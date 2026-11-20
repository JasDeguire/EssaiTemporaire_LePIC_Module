#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_2D_gui.py  --  viewer for the Output/ files of LePIC+ (and the DATA/ files of LePIC2D)

Usage
-----
    python3 GUI/plot_2D_gui.py [path/to/Output]          (default: ./Output, else ./DATA)

Adapted from the LePIC3D_DD 3.4.1 viewer.  The run type is detected from the folder content:
Output/Output_2D/*.mco (LePIC+) or DATA/DATA_2D/*.mco (3-D code) -> plane maps; otherwise the
LePIC2D 2.5D layout (Fortran-unformatted phi.dat, n<k>.dat ...).

Requires: Python >= 3.8, numpy, matplotlib (>= 3.5), tkinter.
Optional: a LaTeX installation (latex + dvipng) for the "Full LaTeX" mode.

LePIC+ (formats from modules/io/mod_output_2d.f90 and modules/simulation/mod_simulation.f90)
--------------------------------------------------------------------------------------------
* Output_2D/it<step>_<name>_<plane>.mco, plane = xy | xz | yz: one set of time-averaged maps per
  save (every nsav steps).  FORMATTED: line 1 = number of nodes "nx+1 ny+1", then the rows from
  the top (largest vertical coordinate) down.  Choose the save with the "Snapshot" box in the top
  bar ("latest" follows the run while auto-update is on).  Horizontal / vertical axes:
  xy -> (x, y), xz -> (x, z), yz -> (y, z).  The planes cross at the middle of the box (the yz plane
  is at the grid when "Holes" = yes).
* Per ptype k (names from Output/particles.out, [e] = 1): n<k> (density), T<k> (temperature, eV),
  u<k>x, u<k>y, u<k>z (mean velocity, m/s), j<k> (|n<v>| with a "vector" block = in-plane flow angle
  -> flux map with streamlines), sour<k>, sink<k> (reaction source / sink rates).
  Fields: phi, Ex, Ey, Ez (time-averaged when the 3D averages are on, ngrid > 0), dr (= dh / lambda_De).
* Files without the it<step>_ prefix: bcnd_*, B_* (|B| in G with field-direction block), Bx/By/Bz_*
  (written once at start-up), and the initial state n<k>_*, phi1_*, Ex/Ey/Ez_* -> snapshot "initial".
* 1D files in Output/: Iw<k>.dat (A) and Pw<k>.dat (W): column 1 = time, then one column per wall label
  igrid = 1..ngrid (igrid = 0 first when xmax < 0); Ptot.dat, Mtot.dat (power and momentum balance),
  skin_depth.dat (RF antenna: <ne>, c/wpe, delta, drive amplitude, P_RF), phi_Te_ne_cntr.dat (probe).
  Times repeated after a restart (the steps between the last backup and the end of a segment are
  redone) are removed when reading: the last record of each time is kept.
* Domain size: first line of ../input_dir/geometry.inp = xmax ymax zmax [cm] (|value| is the length).
"""
import argparse
import os
import re
import shutil
import struct
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from types import SimpleNamespace
from pathlib import Path

import logging
import numpy as np
try:                                     # readers are also used headless by web_gui.py
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
except ImportError:
    tk = ttk = filedialog = messagebox = None

import matplotlib
from matplotlib import cbook, colors as mcolors, mathtext
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
try:
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
except ImportError:
    FigureCanvasTkAgg = NavigationToolbar2Tk = None

logging.getLogger('fontTools').setLevel(logging.ERROR)   # harmless font-timestamp notices on pdf/eps export

# ----------------------------------------------------------------------------
#  Low-level readers
# ----------------------------------------------------------------------------

def _parse_records(buf, endian, msize):
    """Split a Fortran sequential unformatted file into records.
    Returns [(offset, length), ...] or None if the markers are inconsistent."""
    fmt = endian + ('I' if msize == 4 else 'Q')
    recs, pos, n = [], 0, len(buf)
    while pos < n:
        if pos + msize > n:
            return None
        (length,) = struct.unpack_from(fmt, buf, pos)
        start = pos + msize
        end = start + length
        if end + msize > n:
            return None
        (length2,) = struct.unpack_from(fmt, buf, end)
        if length2 != length:
            return None
        recs.append((start, length))
        pos = end + msize
    return recs


def read_field(path):
    """Read a 2D field written by write_data.f90.  Returns (nx, ny, Z) with
    Z[iy, ix] of shape (ny+1, nx+1).  Raises IOError if the file is truncated
    (e.g. while the code is still writing it)."""
    buf = Path(path).read_bytes()
    if len(buf) < 16:
        raise IOError('file empty or incomplete')
    recs = None
    for endian, msize in (('<', 4), ('>', 4), ('<', 8), ('>', 8)):
        recs = _parse_records(buf, endian, msize)
        if recs and len(recs) >= 2:
            break
        recs = None
    if recs is None:
        raise IOError('not a complete Fortran unformatted file (still being written?)')
    s0, l0 = recs[0]
    if l0 == 8:
        n1, n2 = struct.unpack_from(endian + 'ii', buf, s0)
    elif l0 == 16:
        n1, n2 = struct.unpack_from(endian + 'qq', buf, s0)
    else:
        raise IOError('unexpected header record length %d' % l0)
    nelem = (n1 + 1) * (n2 + 1)
    s1, l1 = recs[1]
    if l1 == 8 * nelem:
        dt = np.dtype(endian + 'f8')
    elif l1 == 4 * nelem:
        dt = np.dtype(endian + 'f4')
    else:
        raise IOError('data record has %d bytes, expected %d for a %dx%d grid'
                      % (l1, 8 * nelem, n1 + 1, n2 + 1))
    Z = np.frombuffer(buf, dtype=dt, count=nelem, offset=s1)
    Z = Z.astype(np.float64).reshape((n2 + 1, n1 + 1))  # Fortran (ix fastest)
    return n1, n2, Z


_FIX_EXP = re.compile(r'(?<=\d)([+-]\d{3})(?=\s|$)')   # Fortran "1.0-100" -> "1.0E-100"


def read_table(path):
    """Read an ASCII numeric table with '#' comments (Fortran es16.8 output).
    Tolerates a partially written last line and 3-digit exponents without 'E'."""
    try:
        a = np.loadtxt(path, comments='#', ndmin=2)
        if a.size:
            return a
    except Exception:
        pass
    rows = []
    with open(path, 'r', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in '#!%':
                continue
            try:
                rows.append([float(t) for t in s.split()])
            except ValueError:
                try:
                    rows.append([float(t) for t in _FIX_EXP.sub(r'E\1', s).split()])
                except ValueError:
                    continue
    if not rows:
        raise IOError('no numeric data yet')
    ncol = Counter(len(r) for r in rows).most_common(1)[0][0]
    rows = [r for r in rows if len(r) == ncol]
    return np.array(rows, dtype=float)


def dedupe_time(a):
    """Time series rewritten after a restart: keep the LAST record of each time, sorted by time."""
    if a.ndim != 2 or len(a) < 2 or np.all(np.diff(a[:, 0]) > 0):
        return a
    _, idx = np.unique(a[::-1, 0], return_index=True)
    return a[::-1][idx]


def parse_particles_out(path):
    """Species names in ptype order.  ptype 1 is always the electron."""
    names = []
    if path is not None and Path(path).exists():
        sec = False
        for line in Path(path).read_text(errors='replace').splitlines():
            s = line.strip()
            if not s:
                sec = False
                continue
            if s.endswith(':'):
                sec = s[:-1].strip().lower() in ('ions', 'neutrals')
                continue
            if sec:
                names.append(s.split()[0])
    if not names or names[0].strip('[]').lower() not in ('e', 'e-'):
        names = ['[e]'] + names
    return names


def parse_geometry(path):
    """(Lx, Ly) in cm from geometry.inp:  xmin xmax ymax  (generate_boundary.f90)."""
    for line in Path(path).read_text(errors='replace').splitlines():
        s = line.split('!')[0].strip().replace(',', ' ')
        if not s:
            continue
        xl, xr, yr = (float(t.replace('d', 'e').replace('D', 'e')) for t in s.split()[:3])
        return xr - xl, 2.0 * abs(yr)
    raise ValueError('empty geometry file')


# ---- 3-D code: formatted maps, geometry, flags ------------------------------

def _row_floats(s):
    """One text row -> floats (Fortran 3-digit exponents without 'E', '****' -> nan)."""
    try:
        return np.array(s.split(), dtype=float)
    except ValueError:
        out = []
        for t in _FIX_EXP.sub(r'E\1', s).split():
            try:
                out.append(float(t))
            except ValueError:
                out.append(np.nan)
        return np.array(out, dtype=float)


def read_mco(path):
    """Read a formatted 2D map of the 3-D code (DATA/DATA_2D/*.mco).
    Returns (nx, ny, Z, theta) with Z[iv, ih] of shape (ny+1, nx+1), iv increasing with the
    vertical coordinate (the file lists the rows from the top).  theta is the second block of
    the file (after a line 'vector': in-plane direction angle, j<k>_*.mco and B_*.mco) or None.
    Line 1 is 'nx ny' (cells; LePIC+ writes the node counts nx+1 ny+1) except in Tx/Ty/Tz<k>_*.mco, which have no header: nx, ny are then
    deduced from the data.  Raises IOError for an inconsistent (e.g. half-written) file."""
    lines = Path(path).read_text(errors='replace').splitlines()
    if not lines:
        raise IOError('file empty or incomplete')
    hdr = None
    first = lines[0].split()
    if len(first) == 2:
        try:
            hdr = (int(first[0]), int(first[1]))
        except ValueError:
            hdr = None
    blocks = [[]]
    for ln in lines[1 if hdr else 0:]:
        s = ln.strip()
        if not s:
            continue
        if s.lower().startswith('vector'):
            blocks.append([])
            continue
        blocks[-1].append(_row_floats(s))
    maps = []
    for rows in blocks:
        if not rows:
            continue
        nc = Counter(len(r) for r in rows).most_common(1)[0][0]
        if any(len(r) != nc for r in rows):
            raise IOError('rows of different length (file still being written?)')
        maps.append(np.vstack(rows)[::-1, :])          # first row of the file = top
    if not maps:
        raise IOError('no numeric data yet')
    Z = maps[0]
    ny, nx = Z.shape[0] - 1, Z.shape[1] - 1
    if hdr is not None and hdr not in ((nx, ny), (nx + 1, ny + 1)):   # cells (3-D code) or nodes (LePIC+)
        raise IOError('data are %d x %d nodes but the header says %d x %d cells (incomplete file?)'
                      % (nx + 1, ny + 1, hdr[0], hdr[1]))
    theta = maps[1] if len(maps) > 1 and maps[1].shape == Z.shape else None
    return nx, ny, Z, theta


def read_map(path):
    """(nx, ny, Z, theta) from a 2.5D Fortran-unformatted .dat or a 3-D formatted .mco file."""
    if Path(path).suffix.lower() == '.mco':
        return read_mco(path)
    nx, ny, Z = read_field(path)
    return nx, ny, Z, None


def parse_geometry3d(path):
    """(Lx, Ly, Lz [cm], flag_nmn) from the first line of the 3-D geometry.inp: 'xmax ymax zmax'
    (origin at 0).  A negative value only flags a code option (xmax < 0: Neumann condition on
    the left-hand side, which makes Iw/Pw start at igrid = 0)."""
    for line in Path(path).read_text(errors='replace').splitlines():
        s = line.split('!')[0].strip().replace(',', ' ')
        if not s:
            continue
        v = [float(t.replace('d', 'e').replace('D', 'e')) for t in s.split()[:3]]
        if len(v) < 3:
            raise ValueError('expected "xmax ymax zmax" on the first line')
        return abs(v[0]), abs(v[1]), abs(v[2]), v[0] < 0
    raise ValueError('empty geometry file')


def parse_conditions_flags(path):
    """Flags needed for labelling, from conditions.inp (located through the comment keywords).
    holes: 'Holes (y/n)' = yes -> the yz plane is saved at the grid position instead of xmax/2."""
    out = {'holes': False}
    try:
        for line in Path(path).read_text(errors='replace').splitlines():
            if '!' not in line:
                continue
            body, com = line.split('!', 1)
            toks = body.split()
            if 'holes' in com.lower() and toks:
                out['holes'] = toks[0].lower().startswith('y')
    except OSError:
        pass
    return out


MAPS_DIRS = ('Output_2D', 'DATA_2D')          # LePIC+, 3-D code
_SNAP = re.compile(r'^it(\d+)_(.+)$')
INITIAL = 'initial'                            # unprefixed n<k>/phi1/E maps written at start-up


def maps_dir(datadir):
    """Output/Output_2D (LePIC+) or DATA/DATA_2D (3-D code), or None."""
    for nm in MAPS_DIRS:
        d2 = Path(datadir) / nm
        try:
            if d2.is_dir() and any(p.suffix.lower() == '.mco' for p in d2.iterdir()):
                return d2
        except OSError:
            pass
    return None


def is_3d_dir(datadir):
    return maps_dir(datadir) is not None


def _is_initial(fn):
    return bool(re.match(r'^(?:n\d+|phi1|E[xyz])_(?:xy|xz|yz)\.mco$', fn))


def snapshots(d2):
    """Saved steps (ints, ascending) of Output_2D/it<step>_*.mco, plus INITIAL first if present."""
    steps, init = set(), False
    for nm in os.listdir(d2):
        m = _SNAP.match(nm)
        if m:
            steps.add(int(m.group(1)))
        elif _is_initial(nm):
            init = True
    return ([INITIAL] if init else []) + sorted(steps)


def snapshot_files(d2, step):
    """{name without the it<step>_ prefix: Path} for one snapshot, plus the unprefixed static maps
    (bcnd, B ...).  step = None -> latest; INITIAL -> start-up maps (phi1 is shown as phi)."""
    allf = [p for p in Path(d2).iterdir() if p.is_file() and p.suffix.lower() == '.mco']
    if not any(_SNAP.match(p.name) for p in allf):
        return {p.name: p for p in allf}                 # 3-D code layout: no snapshots
    snaps = snapshots(d2)
    if step is None or (step not in snaps):
        step = snaps[-1] if snaps else None
    out, snap = {}, {}
    for p in allf:
        m = _SNAP.match(p.name)
        if m:
            if step != INITIAL and int(m.group(1)) == step:
                snap[m.group(2)] = p
        elif _is_initial(p.name):
            if step == INITIAL:
                snap[re.sub(r'^phi1_', 'phi_', p.name)] = p
        else:
            out[p.name] = p
    out.update(snap)                                     # the snapshot wins over a stray unprefixed file
    return out


# ----------------------------------------------------------------------------
#  LaTeX / label helpers
# ----------------------------------------------------------------------------

_MT = mathtext.MathTextParser('path')


def check_text(s):
    """None if OK, else a short error message (mathtext validation)."""
    if not s or not cbook.is_math_text(s):
        return None
    try:
        _MT.parse(s, dpi=72, prop=FontProperties(size=12))
        return None
    except Exception as ex:
        return str(ex).strip().splitlines()[0][:140]


def species_math(raw):
    """'[Xe+]' -> r'\\mathrm{Xe}^{+}', '[H2+]' -> r'\\mathrm{H}_{2}^{+}', '[e]' -> 'e'."""
    s = raw.strip()
    if s.startswith('[') and s.endswith(']'):
        s = s[1:-1]
    if s.lower() in ('e', 'e-'):
        return 'e'
    m = re.match(r'^(.*?)([+\-]*)$', s)
    base, ch = m.group(1), m.group(2)
    out = ''
    for tok in re.findall(r'[A-Za-z]+|\d+|.', base):
        if tok.isdigit():
            out += ('_{%s}' % tok) if out else r'\mathrm{%s}' % tok
        else:
            out += r'\mathrm{%s}' % tok
    if ch:
        out += '^{%s}' % ch
    if check_text('$%s$' % out) is not None:
        out = r'\mathrm{%s}' % (re.sub(r'[^A-Za-z0-9]', '', s) or 'p')
    return out


UNIT_FROM_CM = {'m': 1e-2, 'cm': 1.0, 'mm': 10.0}


def axis_coords(n, L_cm, unit):
    """Node coordinates (n+1 points) and the unit label actually used."""
    if unit in UNIT_FROM_CM and L_cm:
        return np.linspace(0.0, L_cm * UNIT_FROM_CM[unit], n + 1), unit
    return np.arange(1, n + 2, dtype=float), 'index'


def gaussian_smooth(Z, sigma, periodic_y=False):
    """Separable Gaussian smoothing of Z[iy, ix] (sigma in cells, 0 = off).
    x: edge-replicated.  y: edge-replicated, or periodic (the last row, iy = ny,
    is the periodic image of row 0 and is kept equal to it)."""
    if not sigma or sigma <= 0:
        return Z
    r = int(np.ceil(3.0 * sigma))
    k = np.arange(-r, r + 1)
    w = np.exp(-0.5 * (k / sigma) ** 2)
    w /= w.sum()

    def conv(A, axis, mode):
        n = A.shape[axis]
        pad = [(0, 0), (0, 0)]
        pad[axis] = (r, r)
        P = np.pad(A, pad, mode=mode)
        out = np.zeros(A.shape, dtype=float)
        for i, wi in enumerate(w):
            sl = [slice(None), slice(None)]
            sl[axis] = slice(i, i + n)
            out += wi * P[tuple(sl)]
        return out

    out = conv(np.asarray(Z, dtype=float), 1, 'edge')
    if periodic_y and out.shape[0] > 2:
        sm = conv(out[:-1, :], 0, 'wrap')
        return np.vstack([sm, sm[:1, :]])
    return conv(out, 0, 'edge')


# ----------------------------------------------------------------------------
#  Dataset registry
# ----------------------------------------------------------------------------

@dataclass
class Entry:
    key: str
    group: str
    label: str
    kind: str                  # field | iw | probe_phi | probe_te | probe_ne
    files: list
    sym: str = ''              # math symbol of the plotted quantity (no $)
    unit: str = ''             # text unit (may contain $..$)
    desc: str = ''
    short: str = ''            # plain-text name for the hover box
    sp: str = ''               # species, math fragment
    extra: dict = field(default_factory=dict)


FIELD_DEFS = [
    ('phi.dat', r'\phi', 'V', 'Electric potential', 'phi'),
    ('Ex.dat', 'E_x', r'V m$^{-1}$', 'Electric field', 'Ex'),
    ('Ey.dat', 'E_y', r'V m$^{-1}$', 'Electric field', 'Ey'),
    ('dr.dat', r'\Delta h/\lambda_{\mathrm{D}e}', '', 'Grid spacing / Debye length', 'dh/lD'),
    ('bcnd.dat', r'\mathrm{bcnd}', '', 'Boundary-condition map', 'bcnd'),
]

SPECIES_DEFS = [   # (file names tried, symbol template, unit, description, hover tag)
    (['n{k}.dat'], 'n_{@}', r'm$^{-3}$', 'Density', 'n'),
    (['T{k}.dat'], 'T_{@}', 'eV', 'Temperature', 'T'),
    (['j{k}x.dat', 'jx{k}.dat'], r'\Gamma_{x,@}', r'm$^{-2}$ s$^{-1}$', 'Flux', 'jx'),
    (['j{k}y.dat', 'jy{k}.dat'], r'\Gamma_{y,@}', r'm$^{-2}$ s$^{-1}$', 'Flux', 'jy'),
    (['j{k}z.dat', 'jz{k}.dat'], r'\Gamma_{z,@}', r'm$^{-2}$ s$^{-1}$', 'Flux', 'jz'),
    (['sour{k}.dat'], r'S^{\mathrm{src}}_{@}', r'm$^{-3}$ s$^{-1}$', 'Source term', 'sour'),
    (['sink{k}.dat'], r'S^{\mathrm{sink}}_{@}', r'm$^{-3}$ s$^{-1}$', 'Sink term', 'sink'),
]

NE_FRACS = [0.5, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]   # columns 4..12 of the probe file
PROBE_NAMES = ('phi_Te_ne_cntr.dat', 'phi_ne_Te_cntr.dat')


# ---- 3-D code: planes, 2D file families, 1D columns -------------------------
PLANES = ('xy', 'xz', 'yz')
TS_KINDS = ('iw', 'probe_phi', 'probe_te', 'probe_ne', 'col')     # entries drawn as curves vs time (or output #)
PLANE_AX = {'xy': ('x', 'y'), 'xz': ('x', 'z'), 'yz': ('y', 'z')}          # (horizontal, vertical)
PLANE_COMPS = {'xy': ('x', 'y', 'z'), 'xz': ('x', 'z', 'y'), 'yz': ('y', 'z', 'x')}  # in-plane h, v, normal
PLANE_NORMAL = {'xy': 'z', 'xz': 'y', 'yz': 'x'}

FIELD_DEFS_3D = [    # file stem, symbol, unit, description, hover tag
    ('phi', r'\phi', 'V', 'Electric potential', 'phi'),
    ('Ex', 'E_x', r'V m$^{-1}$', 'Electric field', 'Ex'),
    ('Ey', 'E_y', r'V m$^{-1}$', 'Electric field', 'Ey'),
    ('Ez', 'E_z', r'V m$^{-1}$', 'Electric field', 'Ez'),
    ('dr', r'\Delta h/\lambda_{\mathrm{D}e}', '', 'Grid spacing / Debye length', 'dh/lD'),
    ('bcnd', r'\mathrm{bcnd}', '', 'Boundary-condition map', 'bcnd'),
]
MAG_DEFS_3D = [      # file stem, symbol, description, hover tag, has direction block
    ('B', 'B', 'Magnetic field magnitude', 'B', True),
    ('Bx', 'B_x', 'Magnetic field', 'Bx', False),
    ('By', 'B_y', 'Magnetic field', 'By', False),
    ('Bz', 'B_z', 'Magnetic field', 'Bz', False),
]
SPECIES_DEFS_3D = [  # file stem template, symbol template, unit, description, hover tag
    ('n{k}', 'n_{@}', r'm$^{-3}$', 'Density', 'n'),
    ('T{k}', 'T_{@}', 'eV', 'Temperature', 'T'),
    ('Tx{k}', 'T_{x,@}', 'eV', 'Temperature x', 'Tx'),
    ('Ty{k}', 'T_{y,@}', 'eV', 'Temperature y', 'Ty'),
    ('Tz{k}', 'T_{z,@}', 'eV', 'Temperature z', 'Tz'),
    ('u{k}x', r'\langle v_x\rangle_{@}', r'm s$^{-1}$', 'Mean velocity x', 'ux'),
    ('u{k}y', r'\langle v_y\rangle_{@}', r'm s$^{-1}$', 'Mean velocity y', 'uy'),
    ('u{k}z', r'\langle v_z\rangle_{@}', r'm s$^{-1}$', 'Mean velocity z', 'uz'),
    ('j{k}x', r'\Gamma_{x,@}', r'm$^{-2}$ s$^{-1}$', 'Flux', 'jx'),
    ('j{k}y', r'\Gamma_{y,@}', r'm$^{-2}$ s$^{-1}$', 'Flux', 'jy'),
    ('j{k}z', r'\Gamma_{z,@}', r'm$^{-2}$ s$^{-1}$', 'Flux', 'jz'),
    ('sour{k}', r'S^{\mathrm{src}}_{@}', r'm$^{-3}$ s$^{-1}$', 'Source term', 'sour'),
    ('sink{k}', r'S^{\mathrm{sink}}_{@}', r'm$^{-3}$ s$^{-1}$', 'Sink term', 'sink'),
]

# 1D files with a time stamp in column 1: (file, group, [(0-based column, name, symbol, unit, description)])
PCT = r'$\%$'
COL_FILES_3D = [
    ('Ptot.dat', 'Totals  (Ptot.dat)', [
        (1, 'Pwall', r'P_{\mathrm{wall}}', 'W', 'Power lost at the walls'),
        (2, 'Pabs', r'P_{\mathrm{abs}}', 'W', 'Absorbed power (heating / RF antenna)'),
        (3, 'Pcoll', r'P_{\mathrm{coll}}', 'W', 'Power lost in collisions'),
        (4, 'Pinj', r'P_{\mathrm{inj}}', 'W', 'Power of the injected particles'),
        (5, 'I_mw', r'I_{\mathrm{mw}}', 'A', 'Electron + negative-ion wall current'),
        (6, 'I_pw', r'I_{\mathrm{pw}}', 'A', 'Positive-ion wall current'),
        (7, 'Ek_ew', r'E_{\mathrm{k,ew}}', 'eV', 'Mean energy of the electrons lost at the walls'),
        (8, 'Ek_iw', r'E_{\mathrm{k,iw}}', 'eV', 'Mean energy of the ions lost at the walls')]),
    ('Mtot.dat', 'Momentum balance  (Mtot.dat)', [
        (1 + 3 * g + c, '%s_%s' % (nm, 'xyz'[c]), r'F_{\mathrm{%s},%s}' % (nm[1:], 'xyz'[c]), un,
         '%s, %s component' % (ds, 'xyz'[c]))
        for g, (nm, un, ds) in enumerate((('Fwall', 'N', 'Force on the walls'),
                                          ('Fabs', 'N', 'Force from the heating / RF field'),
                                          ('Fcoll', 'N', 'Force from collisions'),
                                          ('Finj', r'kg m s$^{-1}$', 'Momentum of the injected particles'),
                                          ('Mtot', r'kg m s$^{-1}$', 'Total momentum of the particles')))
        for c in range(3)]),
    ('skin_depth.dat', 'RF antenna  (skin_depth.dat)', [
        (2, '<ne>', r'\langle n_e\rangle', r'm$^{-3}$', 'Mean electron density over the skin layer / window disk'),
        (3, 'c/wpe', r'c/\omega_{pe}', 'cm', 'Collisionless skin depth from <ne>'),
        (4, 'delta', r'\delta', 'cm', 'Skin depth used for the field'),
        (5, 'E0 / I', r'E_0,\ I', 'V/m or A', 'Drive amplitude (Coil mode: coil current in A)'),
        (6, 'P_RF', r'P_{\mathrm{RF}}', 'W', 'Power absorbed from the RF field over the last period')]),
]
RESID_COLS = [(0, 'res', r'\mathrm{res}', '', 'Poisson solver: residual'),
              (1, 'niter', r'n_{\mathrm{iter}}', '', 'Poisson solver: niter'),
              (2, 'ksor', r'k_{\mathrm{SOR}}', '', 'Poisson solver: ksor')]


def _plane_files(files, stem, used):
    """{plane: Path} for Output_2D/<stem>_<plane>.mco; the names found are added to 'used'."""
    out = {}
    for pl in PLANES:
        fn = '%s_%s.mco' % (stem, pl)
        if fn in files:
            out[pl] = files[fn]
            used.add(fn)
    return out


def scan_mco(d2, names, kmax, step=None):
    """Entries for the formatted 2D maps of Output_2D (one snapshot) or DATA_2D."""
    files = snapshot_files(d2, step)
    used, entries = set(), []

    def avail(pf):
        return '  (%s)' % ', '.join(pl for pl in PLANES if pl in pf)

    def add(key, grp, label, pf, sym, unit, desc, short, sp='', kind='field', **extra):
        first = pf[next(pl for pl in PLANES if pl in pf)]
        fl = first if isinstance(first, list) else [first]
        entries.append(Entry(key, grp, label + avail(pf), kind, fl, sym, unit, desc, short, sp,
                             dict(planes=pf, **extra)))

    for stem, sym, unit, desc, short in FIELD_DEFS_3D:
        pf = _plane_files(files, stem, used)
        if pf:
            add('field:' + stem, 'Fields', '%s  -  %s' % (stem, desc), pf, sym, unit, desc, short)
    for stem, sym, desc, short, vec in MAG_DEFS_3D:
        pf = _plane_files(files, stem, used)
        if pf:
            add('mag:' + stem, 'Magnetic field',
                '%s  -  %s%s' % (stem, desc[:1].lower() + desc[1:], ' + field lines' if vec else ''),
                pf, sym, 'G', desc, short, vec=vec)
    for k in range(1, kmax + 1):
        sp = species_math(names[k - 1])
        grp = 'ptype %d: %s' % (k, names[k - 1])
        tag = names[k - 1].strip('[]')
        for tmpl, sym, unit, desc, short in SPECIES_DEFS_3D:
            stem = tmpl.format(k=k)
            pf = _plane_files(files, stem, used)
            if pf:
                add('sp:%d:%s' % (k, short), grp, '%s  -  %s' % (stem, desc[:1].lower() + desc[1:]), pf,
                    sym.replace('@', sp), unit, desc, '%s[%s]' % (short, tag), sp)
        pf = _plane_files(files, 'j%d' % k, used)  # LePIC+: |n<v>| + in-plane angle block
        if pf:
            add('sp:%d:jvec' % k, grp, 'j%d  -  |flux| with streamlines' % k, pf, r'\Gamma_{%s}' % sp,
                r'm$^{-2}$ s$^{-1}$', 'Flux norm', 'j[%s]' % tag, sp, vec=True)
        fl = {}                                   # |flux|: in-plane h, in-plane v (, normal) component files
        for pl in PLANES:
            c = PLANE_COMPS[pl]
            got = [files.get('j%d%s_%s.mco' % (k, ci, pl)) for ci in c]
            if got[0] is not None and got[1] is not None:
                fl[pl] = got[:2] + ([got[2]] if got[2] is not None else [])
        if fl:
            for pl in PLANES:                     # j<k>_<plane>.mco is then rebuilt from the components
                used.add('j%d_%s.mco' % (k, pl))
            add('sp:%d:jnorm' % k, grp, '|flux| with streamlines', fl, r'\Gamma_{%s}' % sp,
                r'm$^{-2}$ s$^{-1}$', 'Flux norm', 'j[%s]' % tag, sp, kind='flux')
    rest = {}                                     # everything else, grouped by stem
    for fn in files:
        if fn in used:
            continue
        m = re.match(r'^(.*)_(xy|xz|yz)\.mco$', fn)
        rest.setdefault(m.group(1) if m else fn[:-4], {})[m.group(2) if m else 'xy'] = files[fn]
    for stem, pf in sorted(rest.items()):
        add('other:' + stem, 'Other 2D files', '%s  -  not a standard file' % stem, pf,
            r'\mathrm{%s}' % stem.replace('_', r'\_'), '', stem, stem)
    return entries


def scan_cols(files):
    """Time-stamped 1D files of the 3-D code: one entry per column."""
    entries = []
    for fn, grp, cols in COL_FILES_3D:
        if fn not in files:
            continue
        for c, nm, sym, unit, desc in cols:
            entries.append(Entry('col:%s:%d' % (fn, c), grp, '%s  -  %s' % (nm, desc), 'col', [files[fn]],
                                 sym, unit, desc, nm, '', dict(col=c, time=True)))
    if 'residual.dat' in files:
        for c, nm, sym, unit, desc in RESID_COLS:
            entries.append(Entry('col:residual.dat:%d' % c, 'Run diagnostics  (time.dat, residual.dat)',
                                 '%s  -  %s  (vs output #)' % (nm, desc), 'col', [files['residual.dat']],
                                 sym, unit, desc, nm, '', dict(col=c, time=False)))
    return entries


def scan(datadir, step=None):
    """Return (entries, species_names); step = snapshot of Output_2D (None = latest)."""
    datadir = Path(datadir)
    files = {p.name: p for p in datadir.iterdir() if p.is_file()}
    names = parse_particles_out(files.get('particles.out'))
    kmax = len(names)
    for nm in files:
        m = (re.match(r'^(?:n|T|sour|sink|Iw|Pw)(\d+)\.dat$', nm) or re.match(r'^j(\d+)[xyz]\.dat$', nm)
             or re.match(r'^j[xyz](\d+)\.dat$', nm))
        if m:
            kmax = max(kmax, int(m.group(1)))
    while len(names) < kmax:
        names.append('[p%d]' % (len(names) + 1))
    d2 = maps_dir(datadir)
    is3d = d2 is not None
    if is3d:
        for nm in os.listdir(d2):
            nm = _SNAP.sub(r'\2', nm)
            m = re.match(r'^(?:n|T[xyz]?|sour|sink|j|u)(\d+)[xyz]?_(?:xy|xz|yz)\.mco$', nm)
            if m:
                kmax = max(kmax, int(m.group(1)))
        while len(names) < kmax:
            names.append('[p%d]' % (len(names) + 1))

    entries = []
    for fn, sym, unit, desc, short in FIELD_DEFS:
        if fn in files:
            entries.append(Entry('field:' + fn[:-4], 'Fields', '%s  -  %s' % (fn, desc), 'field',
                                 [files[fn]], sym, unit, desc, short))
    for k in range(1, kmax + 1):
        sp = species_math(names[k - 1])
        grp = 'ptype %d: %s' % (k, names[k - 1])
        found = {}
        for cand, sym, unit, desc, short in SPECIES_DEFS:
            for c in cand:
                fn = c.format(k=k)
                if fn in files:
                    entries.append(Entry('sp:%d:%s' % (k, short), grp, '%s  -  %s' % (fn, desc.lower()),
                                         'field', [files[fn]], sym.replace('@', sp), unit, desc,
                                         '%s[%s]' % (short, names[k - 1].strip('[]')), sp))
                    found[short] = files[fn]
                    break
        if 'jx' in found and 'jy' in found:
            fl = [found['jx'], found['jy']] + ([found['jz']] if 'jz' in found else [])
            entries.append(Entry('sp:%d:jnorm' % k, grp, '|flux| with streamlines', 'flux', fl,
                                 r'\Gamma_{%s}' % sp, r'm$^{-2}$ s$^{-1}$', 'Flux norm',
                                 'j[%s]' % names[k - 1].strip('[]'), sp))
    if is3d:
        entries += scan_mco(d2, names, kmax, step)
    for pref, grp, sym, unit, desc in (('Iw', 'Wall currents  (Iw, sum of selected columns)', r'I_{w}', 'A', 'Wall current'),
                                       ('Pw', 'Wall power  (Pw, sum of selected columns)', r'P_{w}',
                                        'W' if is3d else r'eV s$^{-1}$', 'Wall power')):
        for k in range(1, kmax + 1):
            fn = '%s%d.dat' % (pref, k)
            if fn in files:
                entries.append(Entry('%s:%d' % (pref.lower(), k), grp, '%s  -  %s' % (fn, names[k - 1]),
                                     'iw', [files[fn]], sym, unit, desc,
                                     '%s[%s]' % (pref, names[k - 1].strip('[]')), species_math(names[k - 1]),
                                     {'pref': pref}))
    loc = '(xmax/2, ymax/2, zmax/2)' if is3d else '(xmax/2, ymax/2)'
    for pn in PROBE_NAMES:
        if pn in files:
            g = 'Probe  (%s)' % pn
            entries.append(Entry('probe:phi', g, 'phi(t) at ' + loc, 'probe_phi', [files[pn]],
                                 r'\phi', 'V', 'Potential', 'phi'))
            entries.append(Entry('probe:te', g, 'Te(t) at ' + loc, 'probe_te', [files[pn]],
                                 'T_e', 'eV', 'Electron temperature', 'Te'))
            entries.append(Entry('probe:ne', g, 'ne(t), slider on x position', 'probe_ne', [files[pn]],
                                 'n_e', r'm$^{-3}$', 'Electron density', 'ne'))
            break
    if is3d:
        entries += scan_cols(files)
    return entries, names


# ----------------------------------------------------------------------------
#  Application
# ----------------------------------------------------------------------------

TEXT_FIELDS = ('title', 'xlabel', 'ylabel', 'cblabel', 'legend')


def view_labels(plane='xy'):
    """{combobox label: internal mode} for a plane; 'profx' = along the horizontal axis,
    'profy' = along the vertical one, both through the middle node of the other axis
    (for plane 'xy' these are the original 2.5D labels)."""
    a, b = PLANE_AX[plane]
    return {'2D map': 'map',
            'Profile along %s  (at %s = %smax/2)' % (a, b, b): 'profx',
            'Profile along %s  (at %s = %smax/2)' % (b, a, a): 'profy'}


BANDED = 'inferno (10 bands)'
N_BANDS = 10
CMAPS = ['viridis', 'inferno', BANDED, 'plasma', 'magma', 'cividis', 'turbo', 'jet', 'RdBu_r', 'coolwarm', 'gray']
LEGEND_LOCS = ['best', 'upper right', 'upper left', 'lower left', 'lower right',
               'center left', 'center right', 'lower center', 'upper center', 'center']
TEX_MATHTEXT = 'Built-in mathtext'
TEX_FULL = 'Full LaTeX (usetex)'


def usetex_available():
    return bool(shutil.which('latex') and shutil.which('dvipng'))


class App:
    def __init__(self, root, datadir=None):
        self.root = root
        root.title('LePIC+  -  Output viewer')
        root.geometry('1500x900')
        self._ready = False
        self._loading = False
        self._timers = {}
        self._cache = {}
        self.ov = {}                       # user text overrides {(key, mode, field): text}
        self._ok = {}                      # last texts that rendered fine
        self.entries = {}
        self.names = []
        self.cur = None
        self.mode = 'map'
        self.is1d = False
        self.ax = self.cb = self.mesh = self.line = None
        self.hv = None
        self._bg = None
        self._auto_lims = None
        self._sig = None
        self._dirsig = None
        self._note_next = ''
        self.wall_vars = []
        self._wall_built = None
        self.datadir = None
        self.is3d = False                  # Output/Output_2D/*.mco present (LePIC+ / 3-D code)
        self.mapsdir = None
        self.yz_at_grid = False            # 3-D: yz plane saved at the grid position, not at xmax/2
        self.nmn = False                   # 3-D: xmax < 0 in geometry.inp -> Iw/Pw start at igrid = 0
        self.unit_lbl = 'mm'
        self.h_ann = self.h_pt = None
        self._build_ui()
        self._ready = True
        if datadir:
            self.set_datadir(datadir)
        self._poll_id = self.root.after(1000, self.poll)

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        st = ttk.Style()
        st.configure('Bad.TEntry', fieldbackground='#ffd0d0')
        top = ttk.Frame(self.root, padding=(6, 4))
        top.pack(side='top', fill='x')
        ttk.Button(top, text='Open Output folder...', command=self.choose_dir).pack(side='left')
        self.lbl_dir = ttk.Label(top, text='(no folder)')
        self.lbl_dir.pack(side='left', padx=8)
        self.var_snap = tk.StringVar(value='latest')         # LePIC+ it<step>_ snapshots
        self.snap_frame = ttk.Frame(top)
        ttk.Label(self.snap_frame, text='Snapshot:').pack(side='left', padx=(6, 2))
        self.snap_box = ttk.Combobox(self.snap_frame, textvariable=self.var_snap, state='readonly', width=16)
        self.snap_box.pack(side='left')
        self.snap_box.bind('<<ComboboxSelected>>', lambda e: (self._cache.clear(), self.refresh_tree()))
        self.var_plane = tk.StringVar(value='xy')            # 3-D runs only (shown by _apply_mode)
        self.plane_frame = ttk.Frame(top)
        ttk.Label(self.plane_frame, text='Plane:').pack(side='left', padx=(6, 2))
        for pl in PLANES:
            ttk.Radiobutton(self.plane_frame, text=pl, value=pl, variable=self.var_plane).pack(side='left', padx=2)
        ttk.Button(top, text='Reset view', command=lambda: self.render(keep_view=False)).pack(side='right', padx=2)
        ttk.Button(top, text='Reload now', command=self.reload_now).pack(side='right', padx=2)
        self.var_auto = tk.BooleanVar(value=True)
        self.var_int = tk.StringVar(value='60')
        ttk.Label(top, text='s').pack(side='right')
        ttk.Spinbox(top, from_=1, to=3600, increment=5, width=5, textvariable=self.var_int).pack(side='right')
        self.var_int.trace_add('write', lambda *a: self._later('repoll', 600, self._repoll))
        ttk.Checkbutton(top, text='Auto-update every', variable=self.var_auto).pack(side='right', padx=(10, 2))

        pw = ttk.PanedWindow(self.root, orient='horizontal')
        pw.pack(fill='both', expand=True)

        # left: tree
        left = ttk.Frame(pw)
        self.tree = ttk.Treeview(left, show='tree', selectmode='browse')
        self.tree.column('#0', width=330, stretch=True)
        sb = ttk.Scrollbar(left, orient='vertical', command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')
        self.tree.bind('<<TreeviewSelect>>', self.on_select)
        pw.add(left, weight=0)

        # centre: figure
        centre = ttk.Frame(pw)
        centre.rowconfigure(1, weight=1)
        centre.columnconfigure(0, weight=1)
        self.fig = Figure(figsize=(7, 5), dpi=100, layout='constrained')
        self.canvas = FigureCanvasTkAgg(self.fig, master=centre)
        tbf = ttk.Frame(centre)
        tbf.grid(row=0, column=0, sticky='ew')
        self.toolbar = NavigationToolbar2Tk(self.canvas, tbf, pack_toolbar=False)
        self.toolbar.pack(side='left', fill='x')
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky='nsew')
        self.slider_frame = ttk.Frame(centre, padding=(8, 2))
        self.slider_frame.grid(row=2, column=0, sticky='ew')
        self.var_col = tk.IntVar(value=4)
        self.lbl_slider = ttk.Label(self.slider_frame, text='', width=58)
        self.lbl_slider.pack(side='left')
        self.scale = tk.Scale(self.slider_frame, from_=4, to=12, orient='horizontal', resolution=1,
                              showvalue=True, variable=self.var_col, command=self.on_slider)
        self.scale.pack(side='left', fill='x', expand=True)
        self.slider_frame.grid_remove()
        self.lbl_status = ttk.Label(centre, text='', anchor='w')
        self.lbl_status.grid(row=3, column=0, sticky='ew', padx=6, pady=2)
        pw.add(centre, weight=1)

        # right: controls
        right = ttk.Frame(pw, width=330)
        nb = ttk.Notebook(right)
        nb.pack(fill='both', expand=True)
        self._build_tab_plot(nb)
        self._build_tab_text(nb)
        self._build_tab_export(nb)
        pw.add(right, weight=0)

        # keep the figure fitted to its widget (see _sync_fig_size)
        for ev in ('<Map>', '<Configure>'):
            self.canvas.get_tk_widget().bind(ev, lambda e: self._later('fit', 60, self._fit_now), add='+')
        # canvas events
        self.canvas.mpl_connect('draw_event', self._on_draw)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas.mpl_connect('figure_leave_event', lambda ev: self.hide_hover())
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.menu = tk.Menu(self.root, tearoff=0)
        for fmt in ('png', 'pdf', 'eps'):
            self.menu.add_command(label='Save figure as %s...' % fmt.upper(),
                                  command=lambda f=fmt: self.save_as(f))
        self.menu.add_separator()
        self.menu.add_command(label='Reset view', command=lambda: self.render(keep_view=False))

    def _row(self, parent, r, text, widget_factory, **kw):
        lab = ttk.Label(parent, text=text)
        lab.grid(row=r, column=0, sticky='w', padx=4, pady=3)
        w = widget_factory(parent)
        w.grid(row=r, column=1, sticky='ew', padx=4, pady=3, **kw)
        w._lbl = lab                      # so that the label can be renamed (2.5D / 3-D)
        return w

    def _build_tab_plot(self, nb):
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Plot')
        f.columnconfigure(1, weight=1)
        self.var_view = tk.StringVar(value='2D map')
        self.w_view = self._row(f, 0, 'View', lambda p: ttk.Combobox(
            p, textvariable=self.var_view, values=list(view_labels()), state='readonly'))
        self.var_cmap = tk.StringVar(value='viridis')
        self.w_cmap = self._row(f, 1, 'Colormap', lambda p: ttk.Combobox(
            p, textvariable=self.var_cmap, values=CMAPS, state='readonly'))
        self.var_log = tk.BooleanVar(value=False)
        self.w_log = self._row(f, 2, 'Colour scale', lambda p: ttk.Checkbutton(p, text='logarithmic', variable=self.var_log))
        self.var_vmin = tk.StringVar()
        self.var_vmax = tk.StringVar()
        self.w_vmin = self._row(f, 3, 'Colour min (blank=auto)', lambda p: ttk.Entry(p, textvariable=self.var_vmin))
        self.w_vmax = self._row(f, 4, 'Colour max (blank=auto)', lambda p: ttk.Entry(p, textvariable=self.var_vmax))
        self.var_aspect = tk.BooleanVar(value=False)
        self.w_aspect = self._row(f, 5, 'Aspect', lambda p: ttk.Checkbutton(p, text='equal (true geometry)', variable=self.var_aspect))
        ttk.Label(f, text='Smoothing').grid(row=6, column=0, sticky='w', padx=4, pady=3)
        smf = ttk.Frame(f)
        smf.grid(row=6, column=1, sticky='ew', padx=4, pady=3)
        self.var_smon = tk.BooleanVar(value=False)
        self.w_smon = ttk.Checkbutton(smf, text='on', variable=self.var_smon)
        self.w_smon.pack(side='left')
        ttk.Label(smf, text='   sigma [cells]').pack(side='left')
        self.var_smooth = tk.StringVar(value='2')
        self.w_smooth = ttk.Entry(smf, textvariable=self.var_smooth, width=6)
        self.w_smooth.pack(side='left', padx=6)
        self.var_sper = tk.BooleanVar(value=False)
        self.w_sper = self._row(f, 7, 'Smoothing in y', lambda p: ttk.Checkbutton(p, text='periodic (wrap-around)', variable=self.var_sper))
        ttk.Separator(f).grid(row=8, column=0, columnspan=2, sticky='ew', pady=6)
        self.var_unit = tk.StringVar(value='mm')
        self._row(f, 9, 'Length unit', lambda p: ttk.Combobox(
            p, textvariable=self.var_unit, values=['mm', 'cm', 'm', 'grid index'], state='readonly'))
        self.var_Lx = tk.StringVar()
        self.var_Ly = tk.StringVar()
        self.var_Lz = tk.StringVar()
        self.w_Lx = self._row(f, 10, 'Lx = xmax-xmin  [cm]', lambda p: ttk.Entry(p, textvariable=self.var_Lx))
        self.w_Ly = self._row(f, 11, 'Ly = 2|ymax|  [cm]', lambda p: ttk.Entry(p, textvariable=self.var_Ly))
        self.w_Lz = self._row(f, 12, 'Lz = zmax  [cm]', lambda p: ttk.Entry(p, textvariable=self.var_Lz))
        self.w_Lz._lbl.grid_remove()
        self.w_Lz.grid_remove()
        ttk.Button(f, text='Load geometry.inp...', command=self.choose_geometry).grid(
            row=13, column=0, columnspan=2, sticky='ew', padx=4, pady=4)
        # wall-current options (Iw / Pw entries)
        self.iunit_frame = ttk.Frame(f)
        self.iunit_frame.grid(row=14, column=0, columnspan=2, sticky='ew', padx=4)
        ttk.Label(self.iunit_frame, text='Current unit (Iw)').pack(side='left')
        self.var_iunit = tk.StringVar(value='A')
        self.w_iunit = ttk.Combobox(self.iunit_frame, textvariable=self.var_iunit, values=['A', 'mA'],
                                    state='readonly', width=6)
        self.w_iunit.pack(side='left', padx=8)
        self.iunit_frame.grid_remove()
        self.wall_box = ttk.LabelFrame(f, text='Columns summed (Iw / Pw)', padding=4)
        self.wall_box.grid(row=15, column=0, columnspan=2, sticky='ew', padx=4, pady=6)
        self.wall_box.grid_remove()
        # streamline options (flux entries)
        self.stream_box = ttk.LabelFrame(f, text='Flux: streamlines', padding=4)
        self.stream_box.grid(row=15, column=0, columnspan=2, sticky='ew', padx=4, pady=6)
        self.stream_box.columnconfigure(1, weight=1)
        self.var_stream = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.stream_box, text='show streamlines', variable=self.var_stream).grid(
            row=0, column=0, columnspan=2, sticky='w')
        self.var_dens = tk.StringVar(value='1.0')
        ttk.Label(self.stream_box, text='Density').grid(row=1, column=0, sticky='w', pady=2)
        ttk.Entry(self.stream_box, textvariable=self.var_dens, width=8).grid(row=1, column=1, sticky='w', padx=6)
        self.var_scol = tk.StringVar(value='white')
        ttk.Label(self.stream_box, text='Colour').grid(row=2, column=0, sticky='w', pady=2)
        ttk.Combobox(self.stream_box, textvariable=self.var_scol, state='readonly', width=8,
                     values=['white', 'black', 'red', 'orange', 'cyan', 'magenta']).grid(
            row=2, column=1, sticky='w', padx=6)
        self.var_incz = tk.BooleanVar(value=False)
        self.w_incz = ttk.Checkbutton(self.stream_box, text='include jz in the norm', variable=self.var_incz)
        self.w_incz.grid(row=3, column=0, columnspan=2, sticky='w')
        self.stream_box.grid_remove()
        # x-range of the time-series entries (Iw/Pw, probes, Ptot, ...)
        self.tlim_box = ttk.LabelFrame(f, text='Time range (blank = auto)', padding=4)
        self.tlim_box.grid(row=16, column=0, columnspan=2, sticky='ew', padx=4, pady=6)
        self.tlim_box.columnconfigure(1, weight=1)
        self.var_tmin = tk.StringVar()
        self.var_tmax = tk.StringVar()
        self.lbl_tmin = ttk.Label(self.tlim_box, text='min')
        self.lbl_tmin.grid(row=0, column=0, sticky='w', pady=2)
        self.w_tmin = ttk.Entry(self.tlim_box, textvariable=self.var_tmin)
        self.w_tmin.grid(row=0, column=1, sticky='ew', padx=6)
        self.lbl_tmax = ttk.Label(self.tlim_box, text='max')
        self.lbl_tmax.grid(row=1, column=0, sticky='w', pady=2)
        self.w_tmax = ttk.Entry(self.tlim_box, textvariable=self.var_tmax)
        self.w_tmax.grid(row=1, column=1, sticky='ew', padx=6)
        ttk.Button(self.tlim_box, text='Clear', width=6,
                   command=lambda: (self.var_tmin.set(''), self.var_tmax.set(''))).grid(row=0, column=2, rowspan=2, padx=2)
        self.tlim_box.grid_remove()
        for v in (self.var_cmap, self.var_log, self.var_vmin, self.var_vmax, self.var_aspect,
                  self.var_smon, self.var_smooth, self.var_sper, self.var_iunit, self.var_stream, self.var_dens,
                  self.var_scol, self.var_incz):
            v.trace_add('write', lambda *a: self.request_render(keep=True))
        for v in (self.var_view, self.var_unit, self.var_Lx, self.var_Ly, self.var_Lz, self.var_tmin, self.var_tmax):
            v.trace_add('write', lambda *a: self.request_render(keep=False))
        self.var_plane.trace_add('write', lambda *a: self._on_plane())

    def _build_tab_text(self, nb):
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Labels / Fonts')
        f.columnconfigure(1, weight=1)
        self.tv = {k: tk.StringVar() for k in TEXT_FIELDS}
        self.ew = {}
        self.var_show_title = tk.BooleanVar(value=True)
        self.var_show_legend = tk.BooleanVar(value=True)
        ttk.Label(f, text='LaTeX syntax: $...$ for math, e.g.  $n_e$ (m$^{-3}$)',
                  foreground='#555').grid(row=0, column=0, columnspan=2, sticky='w', padx=4)
        rows = [('title', 'Title'), ('xlabel', 'x label'), ('ylabel', 'y label'),
                ('cblabel', 'Colourbar label'), ('legend', 'Legend text')]
        for i, (k, lab) in enumerate(rows, start=1):
            ttk.Label(f, text=lab).grid(row=2 * i - 1, column=0, columnspan=2, sticky='w', padx=4, pady=(6, 0))
            e = ttk.Entry(f, textvariable=self.tv[k])
            e.grid(row=2 * i, column=0, columnspan=2, sticky='ew', padx=4)
            self.ew[k] = e
            self.tv[k].trace_add('write', lambda *a, kk=k: self.on_text_var(kk))
        r = 2 * len(rows) + 1
        self.w_show_title = ttk.Checkbutton(f, text='Show title', variable=self.var_show_title,
                                            command=self.request_restyle)
        self.w_show_title.grid(row=r, column=0, sticky='w', padx=4, pady=(8, 0))
        self.w_show_legend = ttk.Checkbutton(f, text='Show legend', variable=self.var_show_legend,
                                             command=self.request_restyle)
        self.w_show_legend.grid(row=r, column=1, sticky='w', padx=4, pady=(8, 0))
        self.var_legloc = tk.StringVar(value='best')
        self.w_legloc = self._row(f, r + 1, 'Legend position', lambda p: ttk.Combobox(
            p, textvariable=self.var_legloc, values=LEGEND_LOCS, state='readonly'))
        self.var_legloc.trace_add('write', lambda *a: self.request_restyle())
        ttk.Button(f, text='Reset texts to defaults', command=self.reset_texts).grid(
            row=r + 2, column=0, columnspan=2, sticky='ew', padx=4, pady=6)
        ttk.Separator(f).grid(row=r + 3, column=0, columnspan=2, sticky='ew', pady=4)
        self.var_fs = tk.StringVar(value='12')
        self.var_ts = tk.StringVar(value='14')
        self._row(f, r + 4, 'Font size (labels, ticks, legend)', lambda p: ttk.Spinbox(
            p, from_=6, to=40, increment=1, width=6, textvariable=self.var_fs, command=self.request_restyle))
        self._row(f, r + 5, 'Title font size', lambda p: ttk.Spinbox(
            p, from_=6, to=48, increment=1, width=6, textvariable=self.var_ts, command=self.request_restyle))
        for v in (self.var_fs, self.var_ts):
            v.trace_add('write', lambda *a: self.request_restyle())
        self.var_tex = tk.StringVar(value=TEX_MATHTEXT)
        vals = [TEX_MATHTEXT] + ([TEX_FULL] if usetex_available() else [])
        self._row(f, r + 6, 'LaTeX engine', lambda p: ttk.Combobox(
            p, textvariable=self.var_tex, values=vals, state='readonly'))
        self.var_tex.trace_add('write', lambda *a: self.request_render(keep=True))
        if not usetex_available():
            ttk.Label(f, text='(Full LaTeX needs "latex" and "dvipng" in PATH)',
                      foreground='#777').grid(row=r + 7, column=0, columnspan=2, sticky='w', padx=4)

    def _build_tab_export(self, nb):
        f = ttk.Frame(nb, padding=6)
        nb.add(f, text='Export')
        f.columnconfigure(1, weight=1)
        self.var_fixed = tk.BooleanVar(value=True)
        ttk.Checkbutton(f, text='Fixed export size (else as displayed)', variable=self.var_fixed).grid(
            row=0, column=0, columnspan=2, sticky='w', padx=4, pady=4)
        self.var_W = tk.StringVar(value='6.0')
        self.var_H = tk.StringVar(value='4.5')
        self.var_dpi = tk.StringVar(value='300')
        self._row(f, 1, 'Width [inch]', lambda p: ttk.Entry(p, textvariable=self.var_W))
        self._row(f, 2, 'Height [inch]', lambda p: ttk.Entry(p, textvariable=self.var_H))
        self._row(f, 3, 'DPI (png, rasterised maps)', lambda p: ttk.Entry(p, textvariable=self.var_dpi))
        for i, fmt in enumerate(('png', 'pdf', 'eps')):
            ttk.Button(f, text='Save as %s...' % fmt.upper(), command=lambda x=fmt: self.save_as(x)).grid(
                row=4 + i, column=0, columnspan=2, sticky='ew', padx=4, pady=3)
        ttk.Label(f, text='Right-click on the figure for the same menu.\n'
                          'Colour maps are rasterised, everything else\n(axes, text) stays vector in pdf/eps.',
                  foreground='#555').grid(row=8, column=0, columnspan=2, sticky='w', padx=4, pady=8)

    # --------------------------------------------------------------- utils
    def status(self, msg, level='info'):
        col = {'info': '#222', 'warn': '#a60', 'error': '#b00'}[level]
        self.lbl_status.configure(text=msg, foreground=col)

    def _later(self, name, ms, fn):
        if name in self._timers:
            try:
                self.root.after_cancel(self._timers[name])
            except Exception:
                pass
        self._timers[name] = self.root.after(ms, fn)

    def _fit_now(self):
        if self._sync_fig_size():
            self.canvas.draw_idle()

    def request_render(self, keep=True):
        if self._ready and not self._loading:
            self._later('render', 200, lambda: self.render(keep_view=keep))

    def request_restyle(self):
        if self._ready and not self._loading and self.ax is not None:
            self._later('restyle', 700 if self.var_tex.get() == TEX_FULL else 120, self.restyle)

    @staticmethod
    def _float(s):
        s = s.strip().replace(',', '.')
        return float(s) if s else None

    def Lx(self):
        try:
            return self._float(self.var_Lx.get())
        except ValueError:
            return None

    def Ly(self):
        try:
            return self._float(self.var_Ly.get())
        except ValueError:
            return None

    def Lz(self):
        try:
            return self._float(self.var_Lz.get())
        except ValueError:
            return None

    def plane(self):
        """Current plane ('xy' for 2.5D runs)."""
        return self.var_plane.get() if self.is3d else 'xy'

    def ax_names(self):
        return PLANE_AX[self.plane()]

    def plane_L(self):
        """Domain lengths [cm] along the (horizontal, vertical) axes of the current plane."""
        L = {'x': self.Lx(), 'y': self.Ly(), 'z': self.Lz()}
        a, b = self.ax_names()
        return L[a], L[b]

    def view_labels(self):
        return view_labels(self.plane())

    def files_of(self, e):
        """Files of entry e for the current plane ([] if the plane is not available)."""
        pf = e.extra.get('planes')
        if pf is None:
            return e.files
        f = pf.get(self.plane())
        return [] if f is None else (f if isinstance(f, list) else [f])

    @property
    def tkey(self):
        """Key of the text overrides: entry (+ plane in 3-D, since the axis names depend on it)."""
        e = self.cur
        return e.key + ('|' + self.plane() if e.extra.get('planes') else '')

    def _on_plane(self):
        self._refresh_view_labels()
        self.request_render(keep=False)

    def _refresh_view_labels(self):
        old = list(self.w_view['values'])
        idx = old.index(self.var_view.get()) if self.var_view.get() in old else 0
        new = list(self.view_labels())
        self.w_view.configure(values=new)
        self.var_view.set(new[idx])

    def _apply_mode(self):
        """Switch the widgets between the 2.5D and the 3-D layouts."""
        three = self.is3d
        if self.mapsdir is not None and any(_SNAP.match(p) for p in os.listdir(self.mapsdir)):
            self.snap_frame.pack(side='left', after=self.lbl_dir)
        else:
            self.snap_frame.pack_forget()
        if three:
            self.plane_frame.pack(side='left', padx=8)
            self.w_Lz._lbl.grid()
            self.w_Lz.grid()
            self.w_Lx._lbl.configure(text='Lx = xmax  [cm]')
            self.w_Ly._lbl.configure(text='Ly = ymax  [cm]')
            self.w_sper._lbl.configure(text='Smoothing, vertical axis')
            self.w_incz.configure(text='include out-of-plane component in the norm')
        else:
            self.plane_frame.pack_forget()
            self.w_Lz._lbl.grid_remove()
            self.w_Lz.grid_remove()
            self.w_Lx._lbl.configure(text='Lx = xmax-xmin  [cm]')
            self.w_Ly._lbl.configure(text='Ly = 2|ymax|  [cm]')
            self.w_sper._lbl.configure(text='Smoothing in y')
            self.w_incz.configure(text='include jz in the norm')
        self._loading = True
        self.var_incz.set(three)           # 3-D: |j| = full norm, as in the code's own j<k>_<plane>.mco files
        self._loading = False
        self._refresh_view_labels()

    def fs(self):
        try:
            return max(4.0, float(self.var_fs.get()))
        except ValueError:
            return 12.0

    def ts(self):
        try:
            return max(4.0, float(self.var_ts.get()))
        except ValueError:
            return 14.0

    # ------------------------------------------------------------ data dir
    def choose_dir(self):
        d = filedialog.askdirectory(title='Select the Output folder', initialdir=str(self.datadir or '.'))
        if d:
            self.set_datadir(d)

    def set_datadir(self, d):
        d = Path(d)
        if not d.is_dir():
            messagebox.showerror('Error', 'Not a directory: %s' % d)
            return
        if d.name in MAPS_DIRS and not is_3d_dir(d):     # the maps folder itself was selected
            d = d.parent
        self.datadir = d.resolve()
        self.lbl_dir.configure(text=str(self.datadir))
        self._cache.clear()
        self.cur = None
        self.mapsdir = maps_dir(self.datadir)
        self.is3d = self.mapsdir is not None
        self.var_snap.set('latest')
        self.yz_at_grid = False
        self.nmn = False
        self.wall_vars, self._wall_built = [], None
        self._apply_mode()
        for cand in (self.datadir / 'geometry.inp', self.datadir.parent / 'input_dir' / 'geometry.inp'):
            if cand.exists():
                self._load_geometry(cand)
                if self.is3d:
                    self.yz_at_grid = parse_conditions_flags(cand.parent / 'conditions.inp')['holes']
                break
        else:
            self.var_Lx.set('')
            self.var_Ly.set('')
            self.var_Lz.set('')
            self.var_unit.set('grid index')
            self.status('geometry.inp not found: axes in grid index (set Lx, Ly%s in the Plot tab).'
                        % (', Lz' if self.is3d else ''), 'warn')
        self.refresh_tree(select_first=True)

    def _load_geometry(self, path):
        lz = None
        try:
            if self.is3d:
                lx, ly, lz, self.nmn = parse_geometry3d(path)
            else:
                lx, ly = parse_geometry(path)
        except Exception as ex:
            messagebox.showerror('geometry.inp', 'Cannot parse %s:\n%s' % (path, ex))
            return
        self._loading = True
        self.var_Lx.set('%.10g' % lx)
        self.var_Ly.set('%.10g' % ly)
        self.var_Lz.set('%.10g' % lz if lz is not None else '')
        if self.var_unit.get() == 'grid index':
            self.var_unit.set('mm')
        self._loading = False
        self.status('Domain from %s: Lx = %g cm, Ly = %g cm' % (path, lx, ly)
                    + (', Lz = %g cm' % lz if lz is not None else ''))

    def choose_geometry(self):
        p = filedialog.askopenfilename(title='geometry.inp', initialdir=str(self.datadir.parent if self.datadir else '.'))
        if p:
            self._load_geometry(p)
            self.request_render(keep=False)

    def _dir_signature(self):
        try:
            names = tuple(sorted(os.listdir(self.datadir)))
            if self.is3d:
                names += tuple(sorted(os.listdir(self.mapsdir)))
            po = self.datadir / 'particles.out'
            return names, (po.stat().st_mtime_ns if po.exists() else None)
        except OSError:
            return None

    def _snap_step(self):
        """Selected snapshot (None = latest) after refreshing the list of saved steps."""
        snaps = snapshots(self.mapsdir) if self.mapsdir is not None else []
        vals = ['latest'] + [s if s == INITIAL else 'it %d' % s for s in snaps]
        self.snap_box.configure(values=vals)
        cur = self.var_snap.get()
        if cur not in vals:
            self.var_snap.set('latest')
            cur = 'latest'
        if cur == 'latest':
            return None
        return INITIAL if cur == INITIAL else int(cur.split()[1])

    def refresh_tree(self, select_first=False):
        prev = self.cur.key if self.cur else None
        self._dirsig = self._dir_signature()
        try:
            entries, names = scan(self.datadir, self._snap_step())
        except Exception as ex:
            self.status('Scan error: %s' % ex, 'error')
            return
        self.names = names
        self.entries = {e.key: e for e in entries}
        self._loading = True
        self.tree.delete(*self.tree.get_children())
        groups = {}
        for e in entries:
            if e.group not in groups:
                groups[e.group] = self.tree.insert('', 'end', iid='grp:' + e.group, text=e.group, open=True)
            self.tree.insert(groups[e.group], 'end', iid=e.key, text=e.label)
        self._loading = False
        if not entries:
            self.status('No LePIC+ / LePIC2D data files found in %s' % self.datadir, 'warn')
            return
        key = prev if prev in self.entries else (entries[0].key if select_first or prev else None)
        if key:
            if self.cur is not None and key == prev:
                self.cur = self.entries[key]
                self.tree.selection_set(key)    # re-selects silently below
            else:
                self.tree.selection_set(key)
            self.tree.see(key)

    def on_select(self, _=None):
        if self._loading:
            return
        sel = self.tree.selection()
        if not sel or sel[0] not in self.entries:
            return
        new = self.entries[sel[0]]
        if self.cur is not None and new.key == self.cur.key:
            self.cur = new
            return
        self.cur = new
        self._sig = None
        self.render(keep_view=False)

    # ----------------------------------------------------------- loading
    def _cached(self, fn, path):
        st = path.stat()
        sig = (st.st_mtime_ns, st.st_size)
        hit = self._cache.get(path)
        if hit and hit[0] == sig:
            return hit[1]
        val = fn(path)
        self._cache[path] = (sig, val)
        return val

    def load(self, e):
        fl = self.files_of(e)
        if e.kind == 'field':
            nx, ny, Z, th = self._cached(read_map, fl[0])
            return dict(nx=nx, ny=ny, Z=Z, theta=th if e.extra.get('vec') else None)
        if e.kind == 'flux':
            nx, ny, Zx, _ = self._cached(read_map, fl[0])
            _, _, Zy, _ = self._cached(read_map, fl[1])
            if Zx.shape != Zy.shape:
                raise IOError('in-plane flux components have different grids')
            Zz = None
            if self.var_incz.get() and len(fl) > 2:
                _, _, Zz, _ = self._cached(read_map, fl[2])
                if Zz.shape != Zx.shape:
                    raise IOError('out-of-plane flux component has a different grid')
            return dict(nx=nx, ny=ny, Zx=Zx, Zy=Zy, Zz=Zz)
        tab = self._cached(lambda pth: dedupe_time(read_table(pth)), fl[0])
        if e.kind == 'col':
            c = e.extra['col']
            if tab.shape[1] <= c:
                raise IOError('%s has %d columns, column %d requested' % (fl[0].name, tab.shape[1], c + 1))
            if e.extra['time']:
                return dict(t=tab[:, 0] * 1e6, y=tab[:, c], xn='t', xu='\u00b5s')
            return dict(t=np.arange(1, tab.shape[0] + 1, dtype=float), y=tab[:, c], xn='output #', xu='')
        if e.kind == 'iw':
            if tab.shape[1] < 2:
                raise IOError('%s needs at least 2 columns' % e.files[0].name)
            n = tab.shape[1] - 1                       # igrid = 0 .. ngrid
            self._ensure_wall_vars(n)
            sel = np.array([v.get() for v in self.wall_vars], dtype=bool)
            y = tab[:, 1:][:, sel].sum(axis=1) if sel.any() else np.zeros(tab.shape[0])
            if not sel.any():
                self._note_next = 'No column selected: curve set to zero'
            return dict(t=tab[:, 0] * 1e6, y=y, ncols=tab.shape[1], nsel=int(sel.sum()))
        ncols = tab.shape[1]
        t = tab[:, 0] * 1e6
        if e.kind == 'probe_phi':
            if ncols < 2:
                raise IOError('probe file has no column 2')
            return dict(t=t, y=tab[:, 1], ncols=ncols)
        if e.kind == 'probe_te':
            if ncols < 3:
                raise IOError('probe file has no column 3')
            return dict(t=t, y=tab[:, 2], ncols=ncols)
        if ncols < 4:
            raise IOError('probe file has no density column (needs >= 4 columns)')
        col = min(max(self.var_col.get(), 4), ncols)
        return dict(t=t, y=tab[:, col - 1], ncols=ncols, col=col)

    # ------------------------------------------------- Iw / Pw column choice
    def _ensure_wall_vars(self, n):
        if len(self.wall_vars) == n:
            return
        old = [v.get() for v in self.wall_vars]
        # default: wall labels (igrid >= 1) ON, igrid = 0 (3rd dimension) OFF
        self.wall_vars = [tk.BooleanVar(value=(old[i] if i < len(old) else (i >= 1 or self.is3d)))
                          for i in range(n)]
        for v in self.wall_vars:
            v.trace_add('write', lambda *a: self.request_render(keep=True))

    def _build_wall_boxes(self):
        n = len(self.wall_vars)
        if self._wall_built == n:
            return
        for w in self.wall_box.winfo_children():
            w.destroy()
        first = (0 if self.nmn else 1) if self.is3d else 0   # 3-D: Iw/Pw columns are igrid = kn..ngrid, kn = 1 unless xmax < 0
        for i, v in enumerate(self.wall_vars):
            if self.is3d:
                txt = ('col 2: igrid 0 (xmax < 0 in geometry.inp)' if first + i == 0
                       else 'col %d: wall label %d' % (i + 2, first + i))
            else:
                txt = ('col 2: igrid 0, 3rd dimension (unsimulated)' if i == 0
                       else 'col %d: wall label %d' % (i + 2, i))
            ttk.Checkbutton(self.wall_box, text=txt, variable=v).grid(row=i, column=0, sticky='w')
        self._wall_built = n

    # ------------------------------------------------------------ styling
    def apply_rc(self):
        rc = matplotlib.rcParams
        use_tex = self.var_tex.get() == TEX_FULL and usetex_available()
        rc['text.usetex'] = use_tex
        rc['font.family'] = 'serif'
        if use_tex:
            rc['font.serif'] = ['Computer Modern Roman']
            rc['text.latex.preamble'] = r'\usepackage{amsmath}'
        else:
            rc['font.serif'] = ['cmr10', 'DejaVu Serif']
            rc['mathtext.fontset'] = 'cm'
        rc['axes.formatter.use_mathtext'] = True
        rc['axes.formatter.limits'] = (-3, 4)
        rc['axes.unicode_minus'] = False
        rc['axes.linewidth'] = 1.0
        for ax in ('x', 'y'):
            rc[ax + 'tick.direction'] = 'in'
            rc[ax + 'tick.minor.visible'] = True
            rc[ax + 'tick.major.size'] = 5.0
            rc[ax + 'tick.minor.size'] = 2.5
            rc[ax + 'tick.major.width'] = 0.9
            rc[ax + 'tick.minor.width'] = 0.7
        rc['xtick.top'] = True
        rc['ytick.right'] = True
        rc['lines.linewidth'] = 1.6
        rc['pdf.fonttype'] = 42
        rc['ps.fonttype'] = 42
        rc['savefig.dpi'] = 300

    def default_texts(self, e, p):
        u = self.unit_lbl
        a, b = self.ax_names()
        xl = r'$%s$ (%s)' % (a, u) if u != 'index' else r'$%s$ (grid index)' % a
        yl = r'$%s$ (%s)' % (b, u) if u != 'index' else r'$%s$ (grid index)' % b
        q = '$%s$' % e.sym + (' (%s)' % e.unit if e.unit else '')
        pt = ''                                   # plane description appended to the titles (3-D maps)
        if self.is3d and e.extra.get('planes'):
            nrm = PLANE_NORMAL[self.plane()]
            pos = r'x_{\mathrm{grid}}' if (nrm == 'x' and self.yz_at_grid) else r'%s_{\max}/2' % nrm
            pt = r' (%s plane, $%s=%s$)' % (self.plane(), nrm, pos)
        d = dict(title='', xlabel='', ylabel='', cblabel='', legend='')
        if e.kind == 'field':
            d['title'] = ('%s $%s$' % (e.desc, e.sym) if e.sym else e.desc) + pt
            if self.mode == 'map':
                d.update(xlabel=xl, ylabel=yl, cblabel=q)
            else:
                d['xlabel'] = xl if self.mode == 'profx' else yl
                d['ylabel'] = q
                if self.is1d:
                    d['legend'] = '$%s$' % e.sym
                elif self.mode == 'profx':
                    d['legend'] = r'$%s(%s,\,%s_{\max}/2)$' % (e.sym, a, b)
                else:
                    d['legend'] = r'$%s(%s_{\max}/2,\,%s)$' % (e.sym, a, b)
        elif e.kind == 'flux':
            if self.var_incz.get() and len(self.files_of(e)) > 2:
                comp = 'xyz'
            else:
                comp = ''.join(PLANE_COMPS[self.plane()][:2])      # in-plane components (xy for 2.5D)
            sym = r'|\Gamma_{%s,%s}|' % (comp, e.sp)
            qq = '$%s$ (%s)' % (sym, e.unit)
            if self.mode == 'map':
                d.update(title='Flux norm%s, $%s$' % (' and streamlines' if self.var_stream.get() else '', e.sp) + pt,
                         xlabel=xl, ylabel=yl, cblabel=qq)
            else:
                d.update(title='Flux norm, $%s$' % e.sp, xlabel=xl, ylabel=qq, legend='$%s$' % sym)
        elif e.kind == 'iw':
            c = 'P' if e.extra.get('pref') == 'Pw' else 'I'
            unit = 'mA' if (c == 'I' and self.var_iunit.get() == 'mA') else e.unit
            d.update(title='%s, $%s$' % (e.desc, e.sp), xlabel=r'$t$ ($\mu$s)',
                     ylabel=r'$%s_{w}$ (%s)' % (c, unit),
                     legend=r'$%s_{w}^{\mathrm{tot}}$ ($%s$)' % (c, e.sp))
        elif e.kind == 'col':
            d.update(title=e.desc, xlabel=r'$t$ ($\mu$s)' if e.extra['time'] else 'output number',
                     ylabel='$%s$' % e.sym + (' (%s)' % e.unit if e.unit else ''), legend='$%s$' % e.sym)
        elif e.kind == 'probe_phi':
            pos = r'(x_{\max}/2,\,y_{\max}/2' + (r',\,z_{\max}/2)' if self.is3d else ')')
            d.update(title=r'Potential at $%s$' % pos, xlabel=r'$t$ ($\mu$s)',
                     ylabel=r'$\phi$ (V)', legend=r'$\phi(x_{\max}/2)$')
        elif e.kind == 'probe_te':
            pos = r'(x_{\max}/2,\,y_{\max}/2' + (r',\,z_{\max}/2)' if self.is3d else ')')
            d.update(title=r'Electron temperature at $%s$' % pos, xlabel=r'$t$ ($\mu$s)',
                     ylabel=r'$T_e$ (eV)', legend=r'$T_e(x_{\max}/2)$')
        elif e.kind == 'probe_ne':
            col = p['col']
            fr = self.probe_frac(col, p['ncols'])
            yz0 = r'y=y_{\max}/2' + (r',\,z=z_{\max}/2' if self.is3d else '')
            d.update(title=r'Electron density at $%s$' % yz0, xlabel=r'$t$ ($\mu$s)',
                     ylabel=r'$n_e$ (m$^{-3}$)',
                     legend=(r'$n_e(x=%g\,x_{\max})$' % fr) if fr is not None else r'$n_e$ (column %d)' % col)
        return d

    @staticmethod
    def probe_frac(col, ncols):
        i = col - 4
        return NE_FRACS[i] if ncols == 12 and 0 <= i < len(NE_FRACS) else None

    def effective_texts(self):
        d = dict(self.defaults)
        for f in TEXT_FIELDS:
            k = (self.tkey, self.mode, f)
            if k in self.ov:
                d[f] = self.ov[k]
        return d

    def load_text_vars(self):
        self._loading = True
        try:
            for f, t in self.effective_texts().items():
                if self.tv[f].get() != t:
                    self.tv[f].set(t)
        finally:
            self._loading = False

    def on_text_var(self, f):
        if self._loading or self.cur is None:
            return
        self.ov[(self.tkey, self.mode, f)] = self.tv[f].get()
        self.request_restyle()

    def reset_texts(self):
        if self.cur is None:
            return
        for f in TEXT_FIELDS:
            self.ov.pop((self.tkey, self.mode, f), None)
        self.load_text_vars()
        self.restyle()

    def update_controls(self):
        e = self.cur
        is_field = bool(e and e.kind == 'field')
        mapm = bool(e and e.kind in ('field', 'flux')) and self.mode == 'map'

        def st(w, on, ro=False):
            w.configure(state=('readonly' if ro else 'normal') if on else 'disabled')
        st(self.w_view, is_field and not self.is1d, ro=True)
        st(self.w_cmap, mapm, ro=True)
        for w in (self.w_log, self.w_aspect, self.w_vmin, self.w_vmax, self.w_smon):
            st(w, mapm)
        for w in (self.w_smooth, self.w_sper):
            st(w, mapm and self.var_smon.get())
        if e and e.kind in TS_KINDS:
            xu = 'output number' if (e.kind == 'col' and not e.extra['time']) else '\u00b5s'
            self.lbl_tmin.configure(text='min  (%s)' % xu)
            self.lbl_tmax.configure(text='max  (%s)' % xu)
            self.tlim_box.grid()
        else:
            self.tlim_box.grid_remove()
        if e and e.kind == 'flux':
            self.stream_box.grid()
            st(self.w_incz, len(self.files_of(e)) > 2)
        elif e and e.extra.get('vec'):                # 3-D |B| map with field lines from the 'vector' block
            self.stream_box.grid()
            st(self.w_incz, False)
        else:
            self.stream_box.grid_remove()
        if e and e.kind == 'iw' and e.extra.get('pref') == 'Iw':
            self.iunit_frame.grid()
        else:
            self.iunit_frame.grid_remove()
        st(self.ew['cblabel'], mapm)
        for w in (self.ew['legend'], self.w_show_legend):
            st(w, not mapm)
        st(self.w_legloc, not mapm, ro=True)
        if e and e.kind == 'iw':
            self._build_wall_boxes()
            self.wall_box.grid()
        else:
            self.wall_box.grid_remove()
        if e and e.kind == 'probe_ne':
            self.slider_frame.grid()
        else:
            self.slider_frame.grid_remove()

    # ------------------------------------------------------------ building
    def show_message(self, msg):
        self.fig.clear()
        self.ax = self.cb = self.hv = None
        self._bg = None
        self.fig.text(0.5, 0.5, msg, ha='center', va='center', color='#b00', wrap=True)
        self.canvas.draw()

    def _zoomed(self):
        if self.ax is None or self._auto_lims is None:
            return False
        for cur, ref in ((self.ax.get_xlim(), self._auto_lims[0]), (self.ax.get_ylim(), self._auto_lims[1])):
            span = abs(ref[1] - ref[0]) or 1.0
            if abs(cur[0] - ref[0]) > 1e-9 * span or abs(cur[1] - ref[1]) > 1e-9 * span:
                return True
        return False

    def render(self, keep_view=False, auto=False):
        """Full rebuild of the figure from the file(s).  Returns True on success."""
        e = self.cur
        if e is None:
            return False
        fl = self.files_of(e)
        if not fl:                                   # 3-D: quantity not saved in this plane
            avail = ', '.join(pl for pl in PLANES if pl in e.extra.get('planes', {}))
            msg = 'No %s-plane file for %s\n(available planes: %s)' % (self.plane(), e.short, avail or 'none')
            self.show_message(msg)
            self.status(msg.replace('\n', ' '), 'warn')
            self._sig = ()
            return False
        try:
            p = self.load(e)
        except Exception as ex:
            if auto:
                self.status('Waiting for a complete file (%s) ...' % ex, 'warn')
                return False
            msg = '%s:\n%s' % (fl[0].name, ex)
            self.show_message(msg)
            self.status(msg.replace('\n', ' '), 'error')
            return False
        lims = (self.ax.get_xlim(), self.ax.get_ylim()) if (keep_view and self._zoomed()) else None
        self.apply_rc()
        self.fig.clear()
        self._bg = None
        self.ax = self.cb = self.mesh = self.line = None
        self.hv = None
        self.h_ann = self.h_pt = None
        try:
            self.mode = self._build(e, p)
        except Exception as ex:
            self.show_message('Cannot plot %s:\n%s' % (fl[0].name, ex))
            self.status('Cannot plot: %s' % ex, 'error')
            return False
        self.defaults = self.default_texts(e, p)
        self.load_text_vars()
        self.update_controls()
        if e.kind == 'probe_ne':
            self.scale.configure(to=p['ncols'])
            self._slider_label(p)
        self.apply_style()
        self.ax.autoscale_view()
        self._apply_trange(e)
        self._auto_lims = (self.ax.get_xlim(), self.ax.get_ylim())
        if lims:
            self.ax.set_xlim(*lims[0])
            self.ax.set_ylim(*lims[1])
        ok = self.safe_draw()
        self.toolbar.update()
        if fl and ok:
            try:
                self._sig = self.sig_of(e)
            except OSError:
                self._sig = None
            note, self._note_next = self._note_next, ''
            self.status('%s   -   updated %s' % (fl[0].name, time.strftime('%H:%M:%S'))
                        + ('   |   ' + note if note else ''), 'warn' if note else 'info')
        return ok

    def _build(self, e, p):
        ax = self.fig.add_subplot(111)
        self.ax = ax
        unit = self.var_unit.get()
        if e.kind == 'field':
            nx, ny, Z = p['nx'], p['ny'], p['Z']
            self.is1d = (ny == 1)
            a, b = self.ax_names()
            La, Lb = self.plane_L()
            x, ux = axis_coords(nx, La, unit)
            self.unit_lbl = ux
            view = 'profx' if self.is1d else self.view_labels().get(self.var_view.get(), 'map')
            if self.is1d:
                self.line, = ax.plot(x, Z[0, :], color='#1f3a93')
                self.hv = dict(kind='line', x=x, y=Z[0, :], xn=a, xu=ux, yn=e.short)
                return 'profx'
            y, uy = axis_coords(ny, Lb, unit)
            if unit != 'grid index' and ux == 'index':
                self.status('Lx/Ly unknown -> axes in grid index', 'warn')
            if view == 'map':
                sig = self._sigma()
                Z = gaussian_smooth(Z, sig, self.var_sper.get())
                norm = self._norm_for(Z)
                self.mesh = ax.pcolormesh(x, y, Z, shading='nearest', cmap=self._cmap(),
                                          norm=norm, rasterized=True)
                self.cb = self._make_colorbar(ax, Z)
                if p.get('theta') is not None and self.var_stream.get():   # direction lines (3-D |B|)
                    th = np.nan_to_num(p['theta'])
                    self._streamlines(ax, x, y, np.cos(th), np.sin(th), nx, ny, ux == 'index')
                ax.set_xlim(x[0], x[-1])
                ax.set_ylim(y[0], y[-1])
                if self.var_aspect.get():
                    ax.set_aspect('equal')
                self.hv = dict(kind='map', x=x, y=y, Z=Z, xu=ux, yu=uy,
                               yn=e.short + (' (smoothed)' if sig > 0 else ''))
                ax.format_coord = lambda xd, yd: self._coord_text(xd, yd)
                self._make_hover(ax, '+')
                return 'map'
            # cuts through the middle node (same node as the code's probes: n/2 + 1 in Fortran)
            if view == 'profx':
                prof = Z[ny // 2, :]
                xx, xu, xn, yn = x, ux, a, '%s(%s,%smax/2)' % (e.short, a, b)
            else:
                prof = Z[:, nx // 2]
                xx, xu, xn, yn = y, uy, b, '%s(%smax/2,%s)' % (e.short, a, b)
            self.line, = ax.plot(xx, prof, color='#1f3a93')
            self.hv = dict(kind='line', x=xx, y=prof, xn=xn, xu=xu, yn=yn)
            self._make_hover(ax, 'o')
            return view
        if e.kind == 'flux':
            return self._build_flux(e, p, ax, unit)
        # time series
        t, y = p['t'], p['y']
        self.is1d = False
        mode, yu = 'line', ''
        if e.kind == 'iw' and e.extra.get('pref') == 'Iw':
            yu = 'A'
            if self.var_iunit.get() == 'mA':
                y, yu, mode = y * 1e3, 'mA', 'line_mA'
        self.line, = ax.plot(t, y, color='#1f3a93')
        self.hv = dict(kind='line', x=t, y=y, xn=p.get('xn', 't'), xu=p.get('xu', '\u00b5s'), yn=e.short, yu=yu)
        self._make_hover(ax, 'o')
        return mode

    def _apply_trange(self, e):
        """Time-series entries: apply the min/max typed in the 'Time range' box (a blank end stays
        automatic) and rescale y to the points inside the window, so that e.g. the t ~ 0 transient
        of Iw/Ptot does not flatten the rest of the curve."""
        if e.kind not in TS_KINDS or not self.hv or self.hv.get('kind') != 'line':
            return
        lo, hi = self._opt_float(self.var_tmin), self._opt_float(self.var_tmax)
        if lo is None and hi is None:
            return
        x0, x1 = self.ax.get_xlim()
        lo = x0 if lo is None else lo
        hi = x1 if hi is None else hi
        if not lo < hi:
            self._note_next = 'time range ignored (min must be < max)'
            return
        self.ax.set_xlim(lo, hi)
        x, y = np.asarray(self.hv['x'], float), np.asarray(self.hv['y'], float)
        m = (x >= lo) & (x <= hi) & np.isfinite(y)
        if not m.any():
            self._note_next = 'no data inside the requested time range'
            return
        ymin, ymax = float(y[m].min()), float(y[m].max())
        pad = 0.05 * (ymax - ymin) or 0.05 * abs(ymax) or 1.0
        self.ax.set_ylim(ymin - pad, ymax + pad)

    def _sigma(self):
        if not self.var_smon.get():
            return 0.0
        v = self._opt_float(self.var_smooth)
        return max(0.0, v) if v else 0.0

    def _cmap(self):
        """Selected colormap; 'inferno (10 bands)' = inferno sampled at 10 evenly spaced colours."""
        name = self.var_cmap.get()
        if name == BANDED:
            return mcolors.ListedColormap(matplotlib.colormaps['inferno'](np.linspace(0.0, 1.0, N_BANDS)),
                                          name='inferno10')
        return name

    def _make_colorbar(self, ax, Z):
        """Colourbar with end arrows when the user-set limits clip the data; for the banded
        colormap (linear scale) the ticks sit exactly on the band edges."""
        lo, hi = self._opt_float(self.var_vmin), self._opt_float(self.var_vmax)
        zf = Z[np.isfinite(Z)]
        ext_lo = lo is not None and zf.size > 0 and zf.min() < lo
        ext_hi = hi is not None and zf.size > 0 and zf.max() > hi
        extend = {(False, False): 'neither', (True, False): 'min', (False, True): 'max', (True, True): 'both'}[(bool(ext_lo), bool(ext_hi))]
        cb = self.fig.colorbar(self.mesh, ax=ax, pad=0.02, extend=extend)
        if self.var_cmap.get() == BANDED and not isinstance(self.mesh.norm, mcolors.LogNorm):
            a, b = self.mesh.norm.vmin, self.mesh.norm.vmax
            if a is not None and b is not None and b > a:
                cb.set_ticks(np.linspace(a, b, N_BANDS + 1))
                cb.minorticks_off()
        return cb

    def _norm_for(self, Z):
        vmin = self._opt_float(self.var_vmin)
        vmax = self._opt_float(self.var_vmax)
        if vmin is not None and vmax is not None and vmin >= vmax:
            self._note_next = 'Colour min >= colour max: limits ignored'
            vmin = vmax = None
        norm = mcolors.Normalize(vmin, vmax)
        if self.var_log.get():
            pos = Z[np.isfinite(Z) & (Z > 0)]
            if pos.size:
                lo = vmin if (vmin and vmin > 0) else float(pos.min())
                hi = vmax if (vmax and vmax > 0) else float(pos.max())
                if hi > lo:
                    norm = mcolors.LogNorm(lo, hi)
            else:
                self._note_next = 'Log scale ignored: no positive values.'
        return norm

    def _build_flux(self, e, p, ax, unit):
        nx, ny = p['nx'], p['ny']
        self.is1d = (ny == 1)
        La, Lb = self.plane_L()
        x, ux = axis_coords(nx, La, unit)
        self.unit_lbl = ux
        Zx, Zy, Zz = p['Zx'], p['Zy'], p['Zz']
        sig = 0.0 if self.is1d else self._sigma()
        if sig > 0:
            per = self.var_sper.get()
            Zx, Zy = gaussian_smooth(Zx, sig, per), gaussian_smooth(Zy, sig, per)
            if Zz is not None:
                Zz = gaussian_smooth(Zz, sig, per)
        Zn = np.sqrt(Zx ** 2 + Zy ** 2 + (Zz ** 2 if Zz is not None else 0.0))
        nm = '|%s|' % e.short + (' (smoothed)' if sig > 0 else '')
        if self.is1d:
            self.line, = ax.plot(x, Zn[0, :], color='#1f3a93')
            self.hv = dict(kind='line', x=x, y=Zn[0, :], xn=self.ax_names()[0], xu=ux, yn=nm, yu='')
            self._make_hover(ax, 'o')
            return 'profx'
        y, uy = axis_coords(ny, Lb, unit)
        self.mesh = ax.pcolormesh(x, y, Zn, shading='nearest', cmap=self._cmap(),
                                  norm=self._norm_for(Zn), rasterized=True)
        self.cb = self._make_colorbar(ax, Zn)
        if self.var_stream.get():
            self._streamlines(ax, x, y, Zx, Zy, nx, ny, ux == 'index')
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(y[0], y[-1])
        if self.var_aspect.get():
            ax.set_aspect('equal')
        cn = ['j' + c for c in PLANE_COMPS[self.plane()]]              # j<h>, j<v>, j<normal>
        comp = [(cn[0], Zx), (cn[1], Zy)] + ([(cn[2], Zz)] if Zz is not None else [])
        self.hv = dict(kind='map', x=x, y=y, Z=Zn, xu=ux, yu=uy, yn=nm, comp=comp)
        ax.format_coord = lambda xd, yd: self._coord_text(xd, yd)
        self._make_hover(ax, '+')
        return 'map'

    def _streamlines(self, ax, x, y, U, V, nx, ny, index_axes):
        """Streamlines of the in-plane flux direction (drawn in data coordinates)."""
        U, V = np.nan_to_num(U), np.nan_to_num(V)
        La, Lb = self.plane_L()
        if index_axes and La and Lb:                     # axes in cell index: rescale to keep the true direction
            U, V = U / (La / nx), V / (Lb / ny)
        st = max(1, max(nx, ny) // 300)
        xs, ys, U, V = x[::st], y[::st], U[::st, ::st], V[::st, ::st]
        if not np.any(np.hypot(U, V) > 0):
            self._note_next = 'Streamlines skipped: in-plane flux is zero everywhere'
            return
        dens = min(max(self._opt_float(self.var_dens) or 1.0, 0.2), 6.0)
        kw = dict(density=dens, color=self.var_scol.get(), linewidth=0.9, arrowsize=1.1, zorder=3)
        try:
            ax.streamplot(xs, ys, U, V, broken_streamlines=False, **kw)
        except TypeError:                                 # matplotlib < 3.6
            ax.streamplot(xs, ys, U, V, **kw)

    @staticmethod
    def _opt_float(var):
        try:
            s = var.get().strip().replace(',', '.')
            return float(s) if s else None
        except ValueError:
            return None

    def _make_hover(self, ax, marker):
        self.h_pt, = ax.plot([], [], marker, mfc='none', mec='#d00', ms=8, mew=1.4, animated=True,
                             visible=False, zorder=50, label='_hover')
        self.h_ann = ax.annotate('', xy=(0, 0), xytext=(12, 12), textcoords='offset points',
                                 fontsize=9, family='DejaVu Sans', zorder=100, annotation_clip=False,
                                 bbox=dict(boxstyle='round,pad=0.35', fc='#ffffe8', ec='0.25', lw=0.8),
                                 animated=True, visible=False)
        self.h_ann.set_usetex(False)

    def apply_style(self):
        ax = self.ax
        fs, ts = self.fs(), self.ts()
        tx = self.effective_texts()
        bad = {}
        for f in TEXT_FIELDS:
            err = check_text(tx[f]) if not matplotlib.rcParams['text.usetex'] else None
            if err:
                bad[f] = err
                tx[f] = self._ok.get((self.tkey, self.mode, f), '')
        for f in TEXT_FIELDS:
            self.ew[f].configure(style='Bad.TEntry' if f in bad else 'TEntry')
        if bad:
            f0 = next(iter(bad))
            self.status('LaTeX error in "%s": %s' % (f0, bad[f0]), 'error')
        ax.set_title(tx['title'] if self.var_show_title.get() else '', fontsize=ts)
        ax.set_xlabel(tx['xlabel'], fontsize=fs)
        ax.set_ylabel(tx['ylabel'], fontsize=fs)
        ax.tick_params(which='both', labelsize=fs)
        for a in (ax.xaxis, ax.yaxis):
            a.get_offset_text().set_fontsize(fs)
        if self.cb is not None:
            self.cb.set_label(tx['cblabel'], fontsize=fs)
            self.cb.ax.tick_params(which='both', labelsize=fs)
            self.cb.ax.yaxis.get_offset_text().set_fontsize(fs)
        old = ax.get_legend()
        if old is not None:
            old.remove()
        if self.line is not None and self.var_show_legend.get() and tx['legend'].strip():
            ax.legend([self.line], [tx['legend']], loc=self.var_legloc.get(), fontsize=fs,
                      frameon=True, framealpha=1.0, edgecolor='0.2', fancybox=False)
        self._cur_texts = tx

    def restyle(self):
        """Cheap update: texts, fonts, legend only (no data rebuild)."""
        if self.ax is None or self.cur is None:
            return
        self.hide_hover(draw=False)
        self.apply_style()
        self.safe_draw()

    def _sync_fig_size(self):
        """Re-fit the figure to the Tk widget.  matplotlib changes the figure dpi after the
        window is mapped (display scaling); the widget size then does not change, no <Configure>
        event follows, and the figure would stay larger than the widget (bottom/right clipped)."""
        w = self.canvas.get_tk_widget()
        pw, ph = w.winfo_width(), w.winfo_height()
        if pw < 20 or ph < 20:
            return False
        fb = self.fig.bbox
        if abs(fb.width - pw) < 1.5 and abs(fb.height - ph) < 1.5:
            return False
        self.canvas.resize(SimpleNamespace(width=pw, height=ph))
        return True

    def safe_draw(self):
        """Synchronous draw; on a LaTeX failure restore the last good texts."""
        try:
            self._sync_fig_size()
            self.canvas.draw()
        except Exception as ex:
            msg = str(ex).strip().splitlines()
            msg = msg[0] if msg else repr(ex)
            if matplotlib.rcParams['text.usetex']:
                self._note_next = 'Full LaTeX failed (%s) -> switched to built-in mathtext' % msg[:90]
                self.status(self._note_next, 'error')
                self._loading = True
                self.var_tex.set(TEX_MATHTEXT)
                self._loading = False
                self.root.after(10, lambda: self.render(keep_view=True))
                return False
            self._note_next = 'Rendering error: %s (texts reset)' % msg[:120]
            self.status(self._note_next, 'error')
            for f in TEXT_FIELDS:
                self.ov.pop((self.tkey, self.mode, f), None)
            self.root.after(10, lambda: self.render(keep_view=True))
            return False
        for f, t in self._cur_texts.items():
            self._ok[(self.tkey, self.mode, f)] = t
        return True

    # --------------------------------------------------------------- slider
    def _slider_label(self, p):
        col, ncols = p['col'], p['ncols']
        fr = self.probe_frac(col, ncols)
        if fr is None:
            s = 'n_e : column %d of %d' % (col, ncols)
        else:
            s = 'n_e : column %d   x = %g x_max' % (col, fr)
            lx = self.Lx()
            if lx and self.var_unit.get() in UNIT_FROM_CM:
                s += '  (%.4g %s)' % (fr * lx * UNIT_FROM_CM[self.var_unit.get()], self.var_unit.get())
        self.lbl_slider.configure(text=s)

    def on_slider(self, _v=None):
        if self._ready and self.cur is not None and self.cur.kind == 'probe_ne':
            self._later('slider', 40, lambda: self.render(keep_view=False))

    # ---------------------------------------------------------------- hover
    def _on_draw(self, _ev):
        try:
            self._bg = self.canvas.copy_from_bbox(self.fig.bbox)
        except Exception:
            self._bg = None

    def _coord_text(self, xd, yd):
        h = self.hv
        if not h or h['kind'] != 'map':
            return '%.5g, %.5g' % (xd, yd)
        ix = int(np.abs(h['x'] - xd).argmin())
        iy = int(np.abs(h['y'] - yd).argmin())
        return 'x=%.5g y=%.5g  %s=%.5g' % (h['x'][ix], h['y'][iy], h['yn'], h['Z'][iy, ix])

    def hide_hover(self, draw=True):
        if self.h_ann is None or not self.h_ann.get_visible():
            return
        self.h_ann.set_visible(False)
        self.h_pt.set_visible(False)
        if draw and self._bg is not None and self.ax is not None:
            self.canvas.restore_region(self._bg)
            self.canvas.blit(self.fig.bbox)

    def on_motion(self, ev):
        h = self.hv
        if h is None or self.ax is None or self.h_ann is None:
            return
        if ev.inaxes is not self.ax or ev.xdata is None or ev.button is not None:
            self.hide_hover()
            return
        xd, yd = ev.xdata, ev.ydata
        if h['kind'] == 'map':
            ix = int(np.abs(h['x'] - xd).argmin())
            iy = int(np.abs(h['y'] - yd).argmin())
            px, py, val = h['x'][ix], h['y'][iy], h['Z'][iy, ix]
            txt = 'x = %.5g %s\ny = %.5g %s\n%s = %.5g' % (px, h['xu'], py, h['yu'], h['yn'], val)
            for cn, A in h.get('comp', ()):
                txt += '\n%s = %.5g' % (cn, A[iy, ix])
        else:
            i = int(np.abs(h['x'] - xd).argmin())
            px, py = h['x'][i], h['y'][i]
            txt = ('%s = %.5g %s\n%s = %.5g %s' % (h['xn'], px, h['xu'], h['yn'], py, h.get('yu', ''))).replace(' \n', '\n').rstrip()
        self.h_pt.set_data([px], [py])
        self.h_pt.set_visible(True)
        w, hh = self.fig.bbox.width, self.fig.bbox.height
        right = ev.x > 0.62 * w
        top = ev.y > 0.78 * hh
        self.h_ann.xy = (px, py)
        self.h_ann.set_text(txt)
        self.h_ann.set_position((-12 if right else 12, -12 if top else 12))
        self.h_ann.set_ha('right' if right else 'left')
        self.h_ann.set_va('top' if top else 'bottom')
        self.h_ann.set_visible(True)
        if self._bg is None:
            self.canvas.draw()
        self.canvas.restore_region(self._bg)
        self.ax.draw_artist(self.h_pt)
        self.ax.draw_artist(self.h_ann)
        self.canvas.blit(self.fig.bbox)

    def on_press(self, ev):
        if ev.button == 3 and ev.guiEvent is not None:
            try:
                self.menu.tk_popup(ev.guiEvent.x_root, ev.guiEvent.y_root)
            finally:
                self.menu.grab_release()

    # ----------------------------------------------------------------- save
    def save_as(self, fmt):
        if self.ax is None:
            messagebox.showinfo('Save', 'Nothing to save.')
            return
        stem = re.sub(r'[^A-Za-z0-9_.-]+', '_', self.tkey if self.cur else 'figure')
        path = filedialog.asksaveasfilename(
            title='Save figure', defaultextension='.' + fmt, initialfile='%s.%s' % (stem, fmt),
            filetypes=[('%s file' % fmt.upper(), '*.' + fmt)])
        if path:
            self.save_to(path, fmt)

    def save_to(self, path, fmt):
        try:
            W = self._float(self.var_W.get())
            H = self._float(self.var_H.get())
            dpi = self._float(self.var_dpi.get()) or 300
        except ValueError:
            messagebox.showerror('Save', 'Invalid export size / dpi.')
            return False
        self.hide_hover(draw=False)
        old = self.fig.get_size_inches().copy()
        ok = True
        try:
            if self.var_fixed.get() and W and H:
                self.fig.set_size_inches(W, H, forward=False)
            self.fig.savefig(path, format=fmt, dpi=dpi)
            self.status('Saved %s' % path)
        except Exception as ex:
            ok = False
            messagebox.showerror('Save failed', str(ex))
        finally:
            self.fig.set_size_inches(old, forward=False)
            self.safe_draw()
        return ok

    # ------------------------------------------------------- auto-update
    def sig_of(self, e):
        return tuple((p.stat().st_mtime_ns, p.stat().st_size) for p in self.files_of(e))

    def reload_now(self):
        self._cache.clear()
        if self.datadir:
            self.refresh_tree()
        self.render(keep_view=True)

    def poll(self):
        try:
            if self.var_auto.get() and self.datadir is not None:
                if self._dir_signature() != self._dirsig:
                    self.refresh_tree()
                if self.cur is not None:
                    try:
                        sig = self.sig_of(self.cur)
                    except OSError:
                        sig = None
                    if sig is not None and sig != self._sig:
                        self.render(keep_view=True, auto=True)
        except Exception as ex:               # never let the timer die
            self.status('Auto-update error: %s' % ex, 'error')
        finally:
            try:
                ms = int(max(0.2, float(self.var_int.get())) * 1000)
            except ValueError:
                ms = 1000
            self._poll_id = self.root.after(ms, self.poll)

    def _repoll(self):
        """Re-arm the poll timer after the interval was edited."""
        try:
            self.root.after_cancel(self._poll_id)
        except Exception:
            pass
        try:
            ms = int(max(0.2, float(self.var_int.get())) * 1000)
        except ValueError:
            return
        self._poll_id = self.root.after(ms, self.poll)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Viewer for LePIC+ Output/ files (and LePIC2D DATA/ files)')
    ap.add_argument('datadir', nargs='?', default=None,
                    help='Output directory (default: ./Output, else ./DATA): the folder that contains Output_2D/')
    args = ap.parse_args(argv)
    root = tk.Tk()
    d = args.datadir or next((c for c in ('Output', 'DATA') if Path(c).is_dir()), None)
    app = App(root, d)
    if d is None:
        root.after(200, app.choose_dir)
    root.mainloop()


if __name__ == '__main__':
    main()
