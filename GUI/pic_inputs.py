#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pic_inputs.py  --  model layer for the 3D PIC/MCC code input files
===================================================================

Reads, edits and writes  input_dir/conditions.inp, geometry.inp, boundary.inp
(and coil.inp for the Coil RF antenna) exactly the way LePIC+ consumes them:
read_input_cfg / read_coil_turns (modules/config/mod_readConditions.f90),
read_geometry_file / read_boundary_potentials (modules/geometry/mod_boundary.f90)
and generate_boundary (modules/geometry/mod_generateBoundary.f90).
Adapted from the LePIC3D_DD 3.4.1 editor.  It has no GUI dependency; the Tk
interface lives in  pic_input_gui.py.

Design rules (all derived from the Fortran source):

* The Fortran reader is *positional*.  Each READ statement takes the values it
  needs from one line and silently discards the rest of that line (this is why
  "! comments" work).  A line with FEWER values than required makes Fortran
  continue on the next line.  Therefore this module never lets a line shrink:
  every record is always written with all the values the reader needs.
* Lines are modified in place.  A line the user did not touch is written back
  byte-for-byte, so saving without changes produces no diff.
* Sign conventions of the code (e.g. Pabs<0 = off, xr<0 = dielectric) are kept
  in the raw text; the GUI exposes them as check-boxes.
"""
from __future__ import annotations

import difflib
import math
import os
import re
import shutil
import time

# --------------------------------------------------------------------------
#  Tokens
# --------------------------------------------------------------------------
_FLOAT = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$')
_INT = re.compile(r'^[+-]?\d+$')
_TOKEN = re.compile(r'[^\s,]+')


def is_float(t):
    return bool(_FLOAT.match(t or ''))


def is_int(t):
    return bool(_INT.match(t or ''))


def to_f(t):
    return float(t.replace('d', 'e').replace('D', 'e'))


def fmt_float(x):
    """Shortest text that Fortran list-directed input reads back exactly."""
    x = float(x)
    if x == int(x) and abs(x) < 1e15:
        return '%d.' % int(x)
    return repr(x)


def split_comment(raw):
    i = raw.find('!')
    if i < 0:
        return raw, ''
    return raw[:i], raw[i:]


def strip_sign(t):
    return t[1:] if t[:1] in '+-' else t


def is_neg(t):
    return (t or '').startswith('-')


def with_sign(mag, negative):
    mag = strip_sign(mag)
    return ('-' + mag) if negative else mag


class Issue:
    """level: fatal (cannot be saved) > error (code would stop / misbehave)
    > warning > info."""
    __slots__ = ('level', 'msg', 'key', 'tab')
    ORDER = {'fatal': 0, 'error': 1, 'warning': 2, 'info': 3}

    def __init__(self, level, msg, key=None, tab=None):
        self.level, self.msg, self.key, self.tab = level, msg, key, tab

    def __repr__(self):
        return '[%s] %s%s' % (self.level, self.msg, ' (%s)' % self.key if self.key else '')


# --------------------------------------------------------------------------
#  One Fortran record (= one line of the file)
# --------------------------------------------------------------------------
class Rec:
    __slots__ = ('raw', 'toks', 'orig', 'comment', 'new', 'name')

    def __init__(self, raw=None, toks=None, comment='', name=''):
        if raw is not None:
            data, comment = split_comment(raw)
            toks = _TOKEN.findall(data)
        self.raw = raw
        self.toks = list(toks)
        self.orig = list(toks)
        self.comment = comment
        self.new = raw is None
        self.name = name

    def get(self, i):
        return self.toks[i] if i < len(self.toks) else ''

    def set(self, i, v):
        while len(self.toks) <= i:
            self.toks.append('')
        self.toks[i] = v

    def render(self):
        if self.raw is not None and self.toks == self.orig:
            return self.raw
        body = ' '.join(t if t != '' else '?' for t in self.toks)
        return body + (('  ' + self.comment) if self.comment else '')


class InputFile:
    def __init__(self, path):
        self.path = path
        with open(path, 'rb') as fh:
            text = fh.read().decode('latin-1')
        self.eol = '\r\n' if '\r\n' in text else '\n'
        self.final_newline = text.endswith('\n')
        lines = text.replace('\r\n', '\n').split('\n')
        if self.final_newline:
            lines.pop()
        self.original = text.replace('\r\n', '\n')
        self.items = []
        self.notes = []          # parse / legacy-format notes (list of Issue)
        self._parse(lines)

    def render(self):
        out = [it.render() if isinstance(it, Rec) else it for it in self.items]
        return '\n'.join(out) + ('\n' if self.final_newline else '')

    def modified(self):
        return self.render() != self.original

    def _parse(self, lines):
        raise NotImplementedError


def _first_tok(raw):
    d, _ = split_comment(raw)
    t = _TOKEN.findall(d)
    return t[0] if t else None


# --------------------------------------------------------------------------
#  conditions.inp
# --------------------------------------------------------------------------
# (name, number of values, default tokens, comment written for NEW lines)
COND_SPEC = [
    ('rname',  2, ['hydrogen2.dat', '100000'],
     '! name of file storing reactions (in input_dir/), number of time steps'),
    ('Te',     1, ['10.'], '! Te (eV)'),
    ('ng',     1, ['6'], '! # of levels in MG solver'),
    ('epsom',  2, ['1.d-7', '1.6d0'], '! Min. precision in residual, w for SOR [1.6]'),
    ('n',      3, ['128', '64', '64'], '! nx, ny, nz: number of grid CELLS along x, y, z'),
    ('pow',    6, ['0.', '10.', '0.', '100.', '0.', '100.'],
     '! Heating/injection region x<, x>, y<, y>, z<, z> [cm] (yr<0: annular heating)'),
    ('xload',  1, ['10.'], '! Particle loading region [cm]'),
    ('nB',     1, ['0'], '! # of magnetic field entries (<0 positive ions not magnetized)'),
    ('nBgrid', 3, ['1', '1', '1'], '! Grid for the magnetic field map'),
    # -- nB lines  'Bfile' follow --
    ('pabs',   2, ['-1.', '-1.d8'],
     '! Pabs(W) (<0 to turn off), Heating frequency (<0 == off) [, skin depth delta(cm), <0 self-consistent]'),
    ('inj',    2, ['-1.', '2'],
     '! Inject current I(A) (<0 == off), option [1=fixed, 2=equal to ion losses, 3=e-beam, 4=emissive cathode, <0 cosine]'),
    ('tmax',   1, ['1.d-3'], '! t_max(s)'),
    ('kt',     1, ['0.3'], '! time coeff. k ; dt= k*min(dx,dy,dz)/vte'),
    ('n0',     1, ['1.d17'], '! n0 (m-3)'),
    ('ngas',   1, ['1.d20'], '! ng (m-3)'),
    ('npcell', 1, ['20'], '! np_cell: number of particles per cell (<0 for no initial plasma)'),
    ('ngrid',  1, ['4'], '! Number of wall labels (= lines read in boundary.inp; <0 no 3D time averages)'),
    ('nsav',   1, ['1000'], '! frequency for saving data (<0 does not reset averaging)'),
    ('bak',    4, ['1000', '-100', '50.d-6', '100.e-6'],
     '! backup frequency (<0 no averaging), # of profile snapshots (<0 off), time window [t1,t2]'),
    ('keps',   1, ['1.d0'], '! k_eps0 (eps0 -> k_eps0*eps0; <0 iterative P, I_inj or Vgrd)'),
    ('pdf',    7, ['no', '1', '0.', '10.', '0.', '10.', '0'],
     '! Save pdf (y/n), ptype, yl, yr, zl, zr (cm), # of regions n_rg, xl & xr for each region'),
    ('restart', 4, ['no', '1', '1', '1.'],
     '! Restart (yes/no/external), # of OMP threads, # of MPI ranks (of the BACKUP), particle duplication/removal coeff.'),
    ('holes',  6, ['no', '0.8', '0.8', '1', '1', '2'], '! Holes (y/n), Lhy, Lhz (cm), ny*nz holes, grid label'),
    ('neg',    4, ['0.', '1.', '1', '1.d-6'],
     '! jH-(mA/cm2) negative-ion flux, TH(eV), grid label for neg. ion statistics, dielectric capacitance Ca/S (F/m2)'),
    ('sec',    7, ['-0.18', '2', 'No', '0.', '0.', '0.', '0.'],
     '! sec. emission coeff (<0 off), grid label, RF mode (No/Electrode/Wrapped/Planar/Coil), phi0(V), '
     'f0(Hz), phi1(V), f1(Hz) [, E0_RF (V/m; Coil: current I in A) [, R_win(cm)]]'),
    ('coulomb', 1, ['no'], '! Coulomb collisions (no/yes)'),
    ('push',   1, ['momentum'], '! Push scheme (momentum/energy) - EC-PIC (Powis & Kaganovich 2024) if energy'),
]
# Values a line may carry beyond the ones it must have (mod_readConditions reads the
# record as a string and retries with fewer values): record -> number of optional values
OPTIONAL = {'pabs': 1, 'sec': 2}
# RF modes of LePIC+ (case-insensitive; legacy single letters still accepted)
RF_MODES = {'no': 'No', 'n': 'No', 'none': 'No', 'off': 'No',
            'electrode': 'Electrode', 'y': 'Electrode', 'd': 'Electrode',
            'wrapped': 'Wrapped', 'a': 'Wrapped', 'planar': 'Planar', 'p': 'Planar',
            'coil': 'Coil', 'c': 'Coil'}
ANTENNA_MODES = ('Wrapped', 'Planar', 'Coil')


def rf_mode(tok):
    """Canonical RF mode name, or None if the code would stop on it."""
    return RF_MODES.get((tok or '').lower())
HUMAN = {
    'rname': 'reactions file', 'Te': 'electron temperature', 'ng': 'multigrid levels',
    'epsom': 'solver precision / SOR', 'n': 'grid cells', 'pow': 'heating / injection region',
    'xload': 'plasma loading region', 'nB': 'number of magnetic-field entries', 'nBgrid': 'magnetic-field map grid',
    'pabs': 'heating power / frequency', 'inj': 'particle injection', 'tmax': 'simulated time',
    'kt': 'time-step coefficient', 'n0': 'plasma density', 'ngas': 'gas density', 'npcell': 'particles per cell',
    'ngrid': 'number of wall labels', 'nsav': 'saving frequency', 'bak': 'backup frequency / snapshots',
    'keps': 'permittivity scaling', 'pdf': 'velocity-distribution output', 'restart': 'restart',
    'holes': 'grid apertures', 'neg': 'negative-ion emission / dielectric capacitance',
    'sec': 'secondary emission / RF', 'coulomb': 'Coulomb collisions', 'push': 'push scheme',
}
BFILE_DEFAULT = ['s', 'gauss', '1.', '100.', 'By', '1.d10', '0.', '0.', '0.']
BFILE_COMMENT = '! ext. file (y/n/s), name, scaling, B(G), direction, L/x0/y0/z0 (cm)'
BKEYS = ['file', 'name', 'scale', 'B0', 'dir', 'dL', 'x0', 'y0', 'z0']

# key -> (record name, token index)   (GUI key / issue key)
FIELDS = {
    'rname': ('rname', 0), 'nsteps': ('rname', 1), 'Te': ('Te', 0), 'ng': ('ng', 0),
    'eps': ('epsom', 0), 'omega': ('epsom', 1),
    'nx': ('n', 0), 'ny': ('n', 1), 'nz': ('n', 2),
    'xl_pow': ('pow', 0), 'xr_pow': ('pow', 1), 'yl_pow': ('pow', 2),
    'yr_pow': ('pow', 3), 'zl_pow': ('pow', 4), 'zr_pow': ('pow', 5),
    'x_load': ('xload', 0), 'nB': ('nB', 0),
    'nBx': ('nBgrid', 0), 'nBy': ('nBgrid', 1), 'nBz': ('nBgrid', 2),
    'Pabs': ('pabs', 0), 'nu_h': ('pabs', 1), 'delta_heat': ('pabs', 2),
    'I_inj': ('inj', 0), 'opt_inj': ('inj', 1),
    'tmax': ('tmax', 0), 'kt': ('kt', 0), 'n0': ('n0', 0), 'ngas': ('ngas', 0),
    'np_cell': ('npcell', 0), 'ngrid': ('ngrid', 0), 'nsav': ('nsav', 0),
    'nbak': ('bak', 0), 'cnt_plt': ('bak', 1), 't1': ('bak', 2), 't2': ('bak', 3),
    'k_eps0': ('keps', 0),
    'pdf': ('pdf', 0), 'ptype_pdf': ('pdf', 1), 'yl_rg': ('pdf', 2), 'yr_rg': ('pdf', 3),
    'zl_rg': ('pdf', 4), 'zr_rg': ('pdf', 5), 'n_rg': ('pdf', 6),
    'restart': ('restart', 0), 'omp_max': ('restart', 1), 'mpi_max': ('restart', 2),
    'np_dup': ('restart', 3),
    'holes': ('holes', 0), 'Lhy': ('holes', 1), 'Lhz': ('holes', 2), 'nhy': ('holes', 3),
    'nhz': ('holes', 4), 'ind_g': ('holes', 5),
    'jne': ('neg', 0), 'THm': ('neg', 1), 'num_grd': ('neg', 2), 'Ca': ('neg', 3),
    'gam_sec': ('sec', 0), 'igrid_sec': ('sec', 1), 'rf': ('sec', 2), 'phi0': ('sec', 3),
    'f0': ('sec', 4), 'phi1': ('sec', 5), 'f1': ('sec', 6), 'E0': ('sec', 7), 'R_win': ('sec', 8),
    'coulomb': ('coulomb', 0), 'push': ('push', 0),
}
NM_RG = 5          # nm_rg in particle_info.h
MAX_B = 5          # dimension of B_file(5) in particle_info.h
MAX_SEG = 100      # do ig=1,100 in generate_boundary
ME = 9.10938188e-31
QE = 1.60217646e-19
EPS0 = 8.854187817e-12


class Conditions(InputFile):
    """conditions.inp, parsed in the same order as read_input()."""

    def _parse(self, lines):
        self.rec = {}
        self.bfiles = []
        pos = 0
        n_lines = len(lines)

        def skip_blank():
            nonlocal pos
            while pos < n_lines and _first_tok(lines[pos]) is None:
                self.items.append(lines[pos])
                pos += 1
            return pos < n_lines

        def is_end(p):
            t = _first_tok(lines[p])
            return t is not None and t.upper() == 'END'

        def take(name, need, defaults, comment, human=None):
            nonlocal pos
            human = human or HUMAN.get(name, name)
            if not skip_blank() or is_end(pos):
                rec = Rec(toks=defaults[:need], comment=comment, name=name)
                self.items.append(rec)
                self.notes.append(Issue(
                    'warning', 'Line for "%s" is missing in this file (older format?). '
                    'Added with default values - please review it.' % human, tab='Summary'))
                return rec
            rec = Rec(raw=lines[pos], name=name)
            ln = pos + 1
            pos += 1
            self.items.append(rec)
            if len(rec.toks) < need:
                have = len(rec.toks)
                pad = list(defaults[have:need])
                pad += [''] * (need - have - len(pad))
                rec.toks += pad
                self.notes.append(Issue(
                    'warning', 'Line %d ("%s") has %d value(s) but the code reads %d on one line '
                    '(older file format?). Missing values were filled with defaults - please review.'
                    % (ln, human, have, need), tab='Summary'))
            return rec

        spec = {s[0]: s for s in COND_SPEC}

        def get(name, need=None):
            s = spec[name]
            r = take(name, need if need is not None else s[1], s[2], s[3])
            self.rec[name] = r
            return r

        for name in ('rname', 'Te', 'ng', 'epsom', 'n', 'pow', 'xload', 'nB', 'nBgrid'):
            get(name)
        # --- magnetic field entries
        nb = abs(int(self.rec['nB'].get(0))) if is_int(self.rec['nB'].get(0)) else 0
        for k in range(min(nb, MAX_B)):
            r = take('Bfile%d' % k, 9, BFILE_DEFAULT, BFILE_COMMENT, 'magnetic field entry %d' % (k + 1))
            self.bfiles.append(r)
        for name in ('pabs', 'inj', 'tmax', 'kt', 'n0', 'ngas', 'npcell', 'ngrid', 'nsav',
                     'bak', 'keps'):
            get(name)
        # --- pdf line has a variable length (7 + 2*n_rg)
        s = spec['pdf']
        pr = take('pdf', 7, s[2], s[3], 'pdf')
        self.rec['pdf'] = pr
        nrg = int(pr.get(6)) if is_int(pr.get(6)) else 0
        nrg = max(0, min(nrg, 50))
        if len(pr.toks) < 7 + 2 * nrg:
            have = len(pr.toks)
            pr.toks += ['0.'] * (7 + 2 * nrg - have)
            self.notes.append(Issue('warning', 'The diagnostic-region line has fewer values than '
                                    '2 x (number of regions) - missing values set to 0.', tab='Summary'))
        for name in ('restart', 'holes', 'neg', 'sec', 'coulomb', 'push'):
            get(name)
        # --- END marker
        if skip_blank() and is_end(pos):
            rec = Rec(raw=lines[pos], name='END')
            self.items.append(rec)
            pos += 1
        else:
            # look further for END (Fortran keeps reading records until it finds it)
            found = None
            for p in range(pos, n_lines):
                if is_end(p):
                    found = p
                    break
            if found is None:
                rec = Rec(toks=['END'], comment='', name='END')
                self.items.append(rec)
                self.notes.append(Issue('error', 'No "END" line found - the code would fail with an '
                                        'end-of-file error. One was added.', tab='Summary'))
            else:
                for p in range(pos, found):
                    self.items.append(lines[p])
                rec = Rec(raw=lines[found], name='END')
                self.items.append(rec)
                pos = found + 1
        for p in range(pos, n_lines):
            self.items.append(lines[p])
        self.rec['END'] = rec
        # --- LePIC3D_DD files: after the RF line they hold "Save 3D (yes/no) opt" and the
        #     MC/EC flag (0/1), which LePIC+ would read as the Coulomb and push-scheme lines.
        co, pu = self.rec['coulomb'], self.rec['push']
        if not co.new and not pu.new and len(co.toks) >= 2 and is_int(co.get(1)) and pu.get(0) in ('0', '1'):
            ec = pu.get(0)
            co.toks = ['no']
            co.comment = COND_SPEC[[s[0] for s in COND_SPEC].index('coulomb')][3]
            pu.toks = ['energy' if ec == '1' else 'momentum']
            pu.comment = COND_SPEC[[s[0] for s in COND_SPEC].index('push')][3]
            self.notes.append(Issue(
                'warning', 'LePIC3D_DD format found: the "Save 3D data" and MC/EC lines were converted to the '
                'LePIC+ "Coulomb collisions" (set to no) and "push scheme" (%s) lines - please review.'
                % pu.toks[0], key='push', tab='Summary'))
        # tokens beyond what the reader consumes are ignored by Fortran -> report
        extras = []
        for name, r in self.rec.items():
            need = None
            if name in spec:
                need = spec[name][1] + OPTIONAL.get(name, 0)
                if name == 'pdf':
                    nrg = int(r.get(6)) if is_int(r.get(6)) else 0
                    need = 7 + 2 * nrg
            if need is not None and len(r.toks) > need and not r.new:
                extras.append((name, r.toks[need:]))
        for k, r in enumerate(self.bfiles):
            if len(r.toks) > 9 and not r.new:
                extras.append(('B entry %d' % (k + 1), r.toks[9:]))
        self.extras = extras

    # ------------------------------------------------------------ B entries
    def set_nb(self, count):
        """Resize the list of magnetic-field entries (keeps the sign of nB)."""
        count = max(0, min(MAX_B, count))
        while len(self.bfiles) > count:
            r = self.bfiles.pop()
            self.items.remove(r)
        while len(self.bfiles) < count:
            r = Rec(toks=list(BFILE_DEFAULT), comment=BFILE_COMMENT, name='Bfile')
            anchor = self.bfiles[-1] if self.bfiles else self.rec['nBgrid']
            self.items.insert(self.items.index(anchor) + 1, r)
            self.bfiles.append(r)
        old = self.rec['nB'].get(0)
        self.rec['nB'].set(0, ('-' if is_neg(old) else '') + str(count))

    # ------------------------------------------------------------ regions
    def regions(self):
        r = self.rec['pdf']
        n = int(r.get(6)) if is_int(r.get(6)) else 0
        out = []
        for i in range(max(0, n)):
            out.append((r.get(7 + 2 * i), r.get(8 + 2 * i)))
        return out

    def set_regions(self, regs):
        r = self.rec['pdf']
        n_old = int(r.get(6)) if is_int(r.get(6)) else 0
        extras = r.toks[7 + 2 * max(0, n_old):]
        flat = []
        for a, b in regs:
            flat += [a, b]
        r.toks = r.toks[:6] + [str(len(regs))] + flat + extras

    # ---------------------------------------------------------------- tokens
    def tok(self, key):
        rn, i = FIELDS[key]
        return self.rec[rn].get(i)

    def set_tok(self, key, v):
        rn, i = FIELDS[key]
        r = self.rec[rn]
        if rn in OPTIONAL:
            need = [s[1] for s in COND_SPEC if s[0] == rn][0]
            if i >= need and v.strip('+-') == '':
                v = ''                             # cleared optional value: drop it from the line
            r.set(i, v)
            while len(r.toks) > need and r.toks[-1] == '':
                r.toks.pop()
            for j in range(need, len(r.toks)):     # a later optional value keeps its position
                if r.toks[j] == '':
                    r.toks[j] = '0.'
            return
        r.set(i, v)


# --------------------------------------------------------------------------
#  geometry.inp
# --------------------------------------------------------------------------
class Geometry(InputFile):
    def _parse(self, lines):
        self.box = None
        self.segs = []
        pos = 0
        n = len(lines)
        # first non-blank line = box size
        while pos < n and _first_tok(lines[pos]) is None:
            self.items.append(lines[pos])
            pos += 1
        if pos >= n:
            self.box = Rec(toks=['10.', '10.', '10.'], comment='! Box size: xmax, ymax, zmax [cm]', name='box')
            self.items.append(self.box)
            self.notes.append(Issue('error', 'geometry.inp is empty.', tab='Geometry'))
        else:
            self.box = Rec(raw=lines[pos], name='box')
            self.items.append(self.box)
            pos += 1
            if len(self.box.toks) < 3:
                self.notes.append(Issue('error', 'First line of geometry.inp must contain xmax, ymax, zmax.',
                                        tab='Geometry'))
                self.box.toks += ['10.'] * (3 - len(self.box.toks))
        while pos < n:
            if _first_tok(lines[pos]) is None:
                self.items.append(lines[pos])
                pos += 1
                continue
            r = Rec(raw=lines[pos], name='seg')
            if len(r.toks) < 7:
                self.notes.append(Issue(
                    'error', 'geometry.inp line %d has %d value(s); 7 are needed (x<,y<,z<,x>,y>,z>,label).'
                    % (pos + 1, len(r.toks)), tab='Geometry'))
                r.toks += ['0.'] * (7 - len(r.toks))
            self.items.append(r)
            self.segs.append(r)
            pos += 1

    def add_segment(self, after=None, toks=None):
        r = Rec(toks=toks or ['0.', '0.', '0.', '1.', '1.', '1.', '1'],
                comment='! (x<,y<,z<), (x>,y>,z>), label', name='seg')
        if after is None or after not in self.segs:
            anchor = self.segs[-1] if self.segs else self.box
            idx = self.items.index(anchor) + 1
            self.segs.append(r)
        else:
            idx = self.items.index(after) + 1
            self.segs.insert(self.segs.index(after) + 1, r)
        self.items.insert(idx, r)
        return r

    def remove_segment(self, r):
        self.segs.remove(r)
        self.items.remove(r)

    def move_segment(self, r, delta):
        i = self.segs.index(r)
        j = i + delta
        if j < 0 or j >= len(self.segs):
            return
        other = self.segs[j]
        a, b = self.items.index(r), self.items.index(other)
        self.items[a], self.items[b] = other, r
        self.segs[i], self.segs[j] = other, r


# --------------------------------------------------------------------------
#  boundary.inp
# --------------------------------------------------------------------------
class Boundary(InputFile):
    limit = None     # the code reads exactly `ngrid` values; extra lines are not written

    def render(self):
        out = []
        for it in self.items:
            if isinstance(it, Rec):
                if self.limit is not None and it.new and it in self.vals and self.vals.index(it) >= self.limit:
                    continue      # placeholder lines added by the editor beyond `ngrid`
                out.append(it.render())
            else:
                out.append(it)
        return '\n'.join(out) + ('\n' if self.final_newline else '')

    def _parse(self, lines):
        self.vals = []
        for k, raw in enumerate(lines):
            if _first_tok(raw) is None:
                self.items.append(raw)
                continue
            r = Rec(raw=raw, name='V')
            self.items.append(r)
            self.vals.append(r)

    def set_count(self, n):
        """Keep exactly n potential lines (the code reads exactly ngrid)."""
        n = max(0, n)
        while len(self.vals) < n:
            r = Rec(toks=['0.'], comment='! Label #%d: V (volts)' % (len(self.vals) + 1), name='V')
            if self.vals:
                idx = self.items.index(self.vals[-1]) + 1
            else:
                idx = len(self.items)
            self.items.insert(idx, r)
            self.vals.append(r)


# --------------------------------------------------------------------------
#  coil.inp  (RF mode Coil):  N, then N lines  a(cm) x(cm) w
# --------------------------------------------------------------------------
COIL_TURN_COMMENT = '! a(cm) x(cm) relative current'


class CoilFile(InputFile):
    def __init__(self, path):
        if os.path.isfile(path):
            super().__init__(path)
            self.exists = True
        else:                                  # created on save once the Coil mode is chosen
            self.path = path
            self.eol = '\n'
            self.final_newline = True
            self.original = ''
            self.items = []
            self.notes = []
            self.exists = False
            self._parse(['4             ! number of turns',
                         '1.0 -0.5 1.   ' + COIL_TURN_COMMENT, '2.0 -0.5 1.', '3.0 -0.5 1.', '4.0 -0.5 1.'])
            for r in [self.count] + self.turns:
                r.new = True
                r.raw = None

    def _parse(self, lines):
        self.count = None
        self.turns = []
        pos, n = 0, len(lines)
        while pos < n and _first_tok(lines[pos]) is None:
            self.items.append(lines[pos])
            pos += 1
        if pos >= n:
            self.count = Rec(toks=['0'], comment='! number of turns', name='N')
            self.items.append(self.count)
            return
        self.count = Rec(raw=lines[pos], name='N')
        self.items.append(self.count)
        pos += 1
        nt = int(self.count.get(0)) if is_int(self.count.get(0)) else 0
        while pos < n:
            if _first_tok(lines[pos]) is None or len(self.turns) >= max(nt, 0):
                self.items.append(lines[pos])
            else:
                r = Rec(raw=lines[pos], name='turn')
                if len(r.toks) < 3:
                    r.toks += ['1.'] * (3 - len(r.toks))
                    self.notes.append(Issue('error', 'coil.inp line %d needs 3 values: a(cm) x(cm) w.' % (pos + 1),
                                            tab='Emission & RF'))
                self.items.append(r)
                self.turns.append(r)
            pos += 1

    def set_count(self, n):
        n = max(1, n)
        while len(self.turns) > n:
            self.items.remove(self.turns.pop())
        while len(self.turns) < n:
            last = self.turns[-1].toks[:3] if self.turns else ['1.', '-0.5', '1.']
            r = Rec(toks=list(last), comment='', name='turn')
            anchor = self.turns[-1] if self.turns else self.count
            self.items.insert(self.items.index(anchor) + 1, r)
            self.turns.append(r)
        self.count.set(0, str(n))


# --------------------------------------------------------------------------
#  A whole case (the three files, plus coil.inp in Coil mode)
# --------------------------------------------------------------------------
class Case:
    FILES = ('conditions.inp', 'geometry.inp', 'boundary.inp')

    def __init__(self, folder):
        self.folder = os.path.abspath(folder)
        missing = [f for f in self.FILES if not os.path.isfile(os.path.join(self.folder, f))]
        if missing:
            raise FileNotFoundError('Missing in %s: %s' % (self.folder, ', '.join(missing)))
        self.cond = Conditions(os.path.join(self.folder, 'conditions.inp'))
        self.geom = Geometry(os.path.join(self.folder, 'geometry.inp'))
        self.bnd = Boundary(os.path.join(self.folder, 'boundary.inp'))
        self.coil = CoilFile(os.path.join(self.folder, 'coil.inp'))
        self.sync_boundary()
        self.mpi = find_mpi_counts(self.folder)

    def uses_coil(self):
        return rf_mode(self.cond.tok('rf')) == 'Coil'

    @property
    def files(self):
        """Files that are saved: coil.inp only while the Coil mode is selected (or if it exists)."""
        out = [self.cond, self.geom, self.bnd]
        if self.uses_coil() or self.coil.exists:
            out.append(self.coil)
        return out

    def sync_boundary(self):
        """boundary.inp must contain exactly |ngrid| potential values (ngrid < 0: no 3D averages)."""
        t = self.cond.tok('ngrid')
        if is_int(t):
            n = abs(int(t))
            before = len(self.bnd.vals)
            self.bnd.set_count(n)
            self.bnd.limit = n
            if before < n and not getattr(self.bnd, '_noted', False):
                self.bnd._noted = True
                self.bnd.notes.append(Issue(
                    'warning', 'boundary.inp has %d potential value(s) but %d wall labels are declared: the code '
                    'sets the missing potentials to 0 V. Added as 0 V lines - please review.' % (before, n),
                    tab='Potentials'))

    def modified(self):
        self.sync_boundary()
        return any(f.render() != f.original for f in self.files)

    def render_all(self):
        self.sync_boundary()
        return {os.path.basename(f.path): f.render() for f in self.files}

    def diff(self):
        out = []
        new = self.render_all()
        for f in self.files:
            name = os.path.basename(f.path)
            if new[name] != f.original:
                out += list(difflib.unified_diff(
                    f.original.splitlines(), new[name].splitlines(),
                    'a/' + name, 'b/' + name, lineterm='', n=1))
        return out

    def save(self, folder=None):
        """Write the files (backup of the previous ones in <folder>/gui_backups/)."""
        folder = os.path.abspath(folder or self.folder)
        os.makedirs(folder, exist_ok=True)
        new = self.render_all()
        stamp = time.strftime('%Y%m%d_%H%M%S')
        written = []
        for f in self.files:
            name = os.path.basename(f.path)
            dst = os.path.join(folder, name)
            same = os.path.isfile(dst) and open(dst, 'rb').read().decode('latin-1').replace(
                '\r\n', '\n') == new[name]
            if same:
                continue
            if os.path.isfile(dst):
                bdir = os.path.join(folder, 'gui_backups')
                os.makedirs(bdir, exist_ok=True)
                shutil.copy2(dst, os.path.join(bdir, '%s.%s' % (name, stamp)))
            text = new[name]
            if f.eol != '\n':
                text = text.replace('\n', f.eol)
            with open(dst, 'wb') as fh:
                fh.write(text.encode('latin-1', errors='replace'))
            written.append(dst)
        return written


    def input_dir(self):
        """The input_dir/ that LePIC+ reads (./input_dir/*.inp), if this case is a sub-folder of it."""
        parent = os.path.dirname(self.folder)
        return parent if os.path.basename(parent) == 'input_dir' else None

    def activate(self):
        """Copy the saved case files (and its reactions file, if the case has its own) into input_dir/,
        where LePIC+ reads them. Overwritten files go to input_dir/gui_backups/. Returns the copied names."""
        dst = self.input_dir()
        if dst is None:
            raise ValueError('This case is not a sub-folder of input_dir/.')
        names = [os.path.basename(f.path) for f in self.files]
        rn = self.cond.tok('rname')
        if rn and os.path.isfile(os.path.join(self.folder, rn)):
            names.append(rn)
        stamp = time.strftime('%Y%m%d_%H%M%S')
        bdir = os.path.join(dst, 'gui_backups')
        for nm in names:
            src, out = os.path.join(self.folder, nm), os.path.join(dst, nm)
            if not os.path.isfile(src):
                continue
            if os.path.isfile(out):
                os.makedirs(bdir, exist_ok=True)
                shutil.copy2(out, os.path.join(bdir, '%s.%s' % (nm, stamp)))
            shutil.copy2(src, out)
        return names


# --------------------------------------------------------------------------
#  MPI process count, read from the run scripts (the input files do not hold it)
# --------------------------------------------------------------------------
RUN_SCRIPT = 'runModular.sh'      # the only script read for the MPI process count
_LAUNCH = re.compile(r'(?:^|[\s;|&(])(impirun|mpirun|mpiexec(?:\.hydra)?|mpirun_rsh|srun)\b(.*)$')
_NPOPT = re.compile(r'(?:^|\s)(?:-np|--np|-n|--ntasks)(?:=|\s+)(\S+)')


def _slurm_header(text):
    """#SBATCH options -> (ntasks, nodes, tasks_per_node); each int or None."""
    nt = nodes = tpn = None
    for line in text.splitlines():
        m = re.match(r'\s*#SBATCH\s+(.*)$', line)
        if not m:
            continue
        body = m.group(1).split('#')[0].strip()
        for pat, name in ((r'(?:^|\s)(?:-n|--ntasks)(?:=|\s+)(\d+)(?:\s|$)', 'nt'),
                          (r'(?:^|\s)(?:-N|--nodes)(?:=|\s+)(\d+)(?:\s|$)', 'nodes'),
                          (r'--ntasks-per-node(?:=|\s+)(\d+)', 'tpn')):
            mm = re.search(pat, body)
            if mm:
                v = int(mm.group(1))
                if name == 'nt':
                    nt = v
                elif name == 'nodes':
                    nodes = v
                else:
                    tpn = v
    return nt, nodes, tpn


def find_mpi_counts(folder):
    """Read the MPI process count from RUN_SCRIPT in the code folder (the parent of input_dir/).
    Returns a list with at most one dict {script, n, how, note}; n is None when the count is a shell variable."""
    d = os.path.abspath(folder)
    root = None
    for _ in range(3):
        if os.path.basename(d) == 'input_dir':
            root = os.path.dirname(d)
            break
        d = os.path.dirname(d)
    if root is None:
        root = os.path.dirname(os.path.abspath(folder))
    out = []
    names = [RUN_SCRIPT] if os.path.isfile(os.path.join(root, RUN_SCRIPT)) else []
    for nm in names:
        p = os.path.join(root, nm)
        try:
            if not os.path.isfile(p) or os.path.getsize(p) > 200000:
                continue
            raw = open(p, 'rb').read()
        except OSError:
            continue
        if b'\0' in raw[:4096]:
            continue                                   # binary (e.g. the executable itself)
        text = raw.decode('latin-1')
        if not re.search(r'\b(mpirun|impirun|mpiexec|srun|mpirun_rsh)\b|#SBATCH', text):
            continue
        found = None
        for line in text.splitlines():
            if line.lstrip().startswith('#'):
                continue
            m = _LAUNCH.search(line)
            if not m:
                continue
            o = _NPOPT.search(m.group(2))
            if o:
                v = o.group(1)
                if v.isdigit():
                    found = (int(v), '%s %s' % (m.group(1), o.group(0).strip()), None)
                else:
                    found = (None, m.group(1), 'the count is a variable (%s)' % v)
                break
        nt, nodes, tpn = _slurm_header(text)
        slurm = None
        note = None
        if nt:
            slurm = (nt, '#SBATCH -n %d' % nt)
            if nodes and tpn and nodes * tpn != nt:
                note = '-N %d x %d per node = %d differs from -n %d' % (nodes, tpn, nodes * tpn, nt)
        elif nodes and tpn:
            slurm = (nodes * tpn, '#SBATCH -N %d x %d tasks per node' % (nodes, tpn))
        if found and found[0] is not None:
            out.append(dict(script=nm, n=found[0], how=found[1], note=None))
        elif slurm:
            out.append(dict(script=nm, n=slurm[0], how=slurm[1], note=note))
        elif found:
            out.append(dict(script=nm, n=None, how=found[1], note=found[2]))
    return out


def mpi_summary(found):
    if not found:
        return '%s not found next to input_dir/ (nz divisibility is not checked)' % RUN_SCRIPT
    return ',  '.join('%s: %s' % (f['script'], f['n'] if f['n'] else '? (count is a variable, not checked)')
                      for f in found)


# --------------------------------------------------------------------------
#  Decoding  (raw text -> numbers with the semantics of read_input())
# --------------------------------------------------------------------------
class P(dict):
    """Decoded parameters; attribute access, missing/invalid values are None."""
    __getattr__ = dict.get


def _num(rec_get, key, kind, issues, label, tab):
    t = rec_get
    if t == '':
        issues.append(Issue('fatal', '%s: value is missing.' % label, key, tab))
        return None
    if kind == 'int':
        if is_int(t):
            return int(t)
        if is_float(t) and to_f(t) == int(to_f(t)):
            issues.append(Issue('warning', '%s: "%s" is written as a real number where the code reads an integer. '
                                'Some Fortran compilers accept this, others stop; "%d" is safer.'
                                % (label, t, int(to_f(t))), key, tab))
            return int(to_f(t))
        issues.append(Issue('fatal', '%s: "%s" is not an integer (the code reads an integer here).'
                            % (label, t), key, tab))
        return None
    if is_float(t):
        return to_f(t)
    issues.append(Issue('fatal', '%s: "%s" is not a number.' % (label, t), key, tab))
    return None


def _yn(t):
    c = (t or ' ')[0]
    return c in 'yY'


LABELS = {
    'Te': 'Electron temperature', 'ng': 'Multigrid levels', 'eps': 'Solver precision',
    'omega': 'SOR parameter', 'nx': 'Cells in x', 'ny': 'Cells in y', 'nz': 'Cells in z',
    'xl_pow': 'Heating region x_min', 'xr_pow': 'Heating region x_max', 'yl_pow': 'Heating region y_min',
    'yr_pow': 'Heating region y_max', 'zl_pow': 'Heating region z_min', 'zr_pow': 'Heating region z_max',
    'x_load': 'Loading region', 'nB': 'Number of magnetic-field entries',
    'nBx': 'B-map grid x', 'nBy': 'B-map grid y', 'nBz': 'B-map grid z',
    'Pabs': 'Absorbed power', 'nu_h': 'Heating frequency', 'I_inj': 'Injected current',
    'opt_inj': 'Injection option', 'tmax': 'Simulation time', 'kt': 'Time-step coefficient',
    'n0': 'Plasma density', 'ngas': 'Gas density', 'np_cell': 'Particles per cell',
    'ngrid': 'Number of wall labels', 'nsav': 'Saving frequency', 'nbak': 'Backup frequency',
    'cnt_plt': 'Number of snapshots', 't1': 'Snapshot window start', 't2': 'Snapshot window end',
    'k_eps0': 'Permittivity scaling', 'ptype_pdf': 'Species for distribution function',
    'yl_rg': 'Region y_min', 'yr_rg': 'Region y_max', 'zl_rg': 'Region z_min', 'zr_rg': 'Region z_max',
    'n_rg': 'Number of regions', 'omp_max': 'OpenMP threads (backup)', 'mpi_max': 'MPI ranks (backup)',
    'np_dup': 'Duplication/removal coefficient', 'Lhy': 'Hole size Ly', 'Lhz': 'Hole size Lz',
    'nhy': 'Holes along y', 'nhz': 'Holes along z', 'ind_g': 'Grid label',
    'jne': 'Negative-ion flux', 'THm': 'Emission temperature', 'num_grd': 'Grid for ion statistics',
    'Ca': 'Dielectric capacitance', 'gam_sec': 'Secondary emission coefficient',
    'igrid_sec': 'Emitting electrode label', 'phi0': 'RF amplitude 1', 'f0': 'RF frequency 1',
    'phi1': 'RF amplitude 2', 'f1': 'RF frequency 2', 'E0': 'RF antenna field', 'opt3d': '3D data option',
    'nsteps': 'Number of time steps', 'delta_heat': 'Heating skin depth', 'R_win': 'Window radius',
}
INT_KEYS = {'ng', 'nx', 'ny', 'nz', 'nB', 'nBx', 'nBy', 'nBz', 'opt_inj', 'np_cell', 'ngrid', 'nsav',
            'nbak', 'cnt_plt', 'ptype_pdf', 'n_rg', 'omp_max', 'mpi_max', 'nhy', 'nhz', 'ind_g',
            'num_grd', 'igrid_sec', 'nsteps'}
FLOAT_KEYS = {'Te', 'eps', 'omega', 'xl_pow', 'xr_pow', 'yl_pow', 'yr_pow', 'zl_pow', 'zr_pow', 'x_load',
              'Pabs', 'nu_h', 'I_inj', 'tmax', 'kt', 'n0', 'ngas', 't1', 't2', 'k_eps0', 'yl_rg', 'yr_rg',
              'zl_rg', 'zr_rg', 'np_dup', 'Lhy', 'Lhz', 'jne', 'THm', 'Ca', 'gam_sec', 'phi0', 'f0',
              'phi1', 'f1', 'E0', 'delta_heat', 'R_win'}
OPTIONAL_KEYS = {'delta_heat': '0.', 'E0': '0.', 'R_win': '0.'}


def decode(case):
    """Return (P, issues).  P holds typed values; sign conventions are expanded."""
    c = case.cond
    iss = []
    p = P()
    for key in INT_KEYS:
        p[key] = _num(c.tok(key), key, 'int', iss, LABELS.get(key, key), 'Summary')
    for key in FLOAT_KEYS:
        t = c.tok(key)
        if t == '' and key in OPTIONAL_KEYS:
            t = OPTIONAL_KEYS[key]                 # absent optional value: the reader sets it to 0
        p[key] = _num(t, key, 'float', iss, LABELS.get(key, key), 'Summary')
    p['rname'] = c.tok('rname')
    if not p['rname']:
        iss.append(Issue('fatal', 'Reactions file name is missing.', 'rname'))
    p['pdf'] = _yn(c.tok('pdf'))
    rs = c.tok('restart')[:1].lower()
    p['restart'] = 1 if rs == 'y' else (2 if rs == 'e' else 0)
    p['holes'] = _yn(c.tok('holes'))
    p['rf_mode'] = rf_mode(c.tok('rf'))
    p['antenna'] = p['rf_mode'] in ANTENNA_MODES
    p['coulomb'] = _yn(c.tok('coulomb'))
    p['push'] = c.tok('push')
    p['flag_avg3D'] = not is_neg(c.tok('ngrid'))
    if p['ngrid'] is not None:
        p['ngrid'] = abs(p['ngrid'])               # the sign only switches the 3D time averages off
    # sign conventions
    p['flag_ahp'] = p['yr_pow'] is not None and p['yr_pow'] < 0
    p['yr_pow_abs'] = abs(p['yr_pow']) if p['yr_pow'] is not None else None
    p['B_ions_unmagnetized'] = is_neg(c.tok('nB'))
    p['flag_convP'] = p['k_eps0'] is not None and p['k_eps0'] < 0
    # B entries
    be = []
    for k, r in enumerate(c.bfiles):
        d = {'k': k}
        d['file'] = r.get(0)[:1].lower()
        d['name'] = r.get(1)
        for j, kk in enumerate(BKEYS[2:], start=2):
            if kk == 'dir':
                d['dir'] = r.get(4)
            else:
                d[kk] = _num(r.get(j), 'B%d.%s' % (k, kk), 'float', iss, 'Magnetic field entry %d: %s' % (k + 1, kk),
                             'Magnetic field')
        be.append(d)
    p['B'] = be
    return p, iss


def decode_geometry(case):
    g = case.geom
    iss = []
    box = []
    for i, nm in enumerate(('xmax', 'ymax', 'zmax')):
        box.append(_num(g.box.get(i), 'g.' + nm, 'float', iss, 'Box size ' + nm, 'Geometry'))
    segs = []
    for k, r in enumerate(g.segs):
        vals = []
        for i in range(6):
            vals.append(_num(r.get(i), 'g.seg%d.%d' % (k, i), 'float', iss,
                             'Segment %d: coordinate %d' % (k + 1, i + 1), 'Geometry'))
        ind = _num(r.get(6), 'g.seg%d.6' % k, 'int', iss, 'Segment %d: label' % (k + 1), 'Geometry')
        segs.append(vals + [ind])
    return box, segs, iss


def decode_boundary(case):
    iss = []
    v = []
    for k, r in enumerate(case.bnd.vals):
        v.append(_num(r.get(0), 'b.V%d' % (k + 1), 'float', iss, 'Potential of label %d' % (k + 1), 'Potentials'))
    return v, iss


# --------------------------------------------------------------------------
#  Derived numbers
# --------------------------------------------------------------------------
def ns_heat_of(kt):
    if kt >= 0.1:
        return 4
    if kt >= 0.05:
        return 20
    return 40


def is_pow2(n):
    return n > 0 and (n & (n - 1)) == 0


def derived(p, box):
    """Numbers the code prints at start-up (lde, dt, wpe*dt ...). Returns dict or None."""
    try:
        nx, ny, nz = p['nx'], p['ny'], p['nz']
        Te, n0, kt, keps = p['Te'], p['n0'], p['kt'], abs(p['k_eps0'])
        L = [abs(b) for b in box]
        if None in (nx, ny, nz, Te, n0, kt, keps) or None in L or min(nx, ny, nz) < 1:
            return None
        if Te <= 0 or n0 <= 0 or kt <= 0 or keps == 0:
            return None
        h = [L[0] / nx * 1e-2, L[1] / ny * 1e-2, L[2] / nz * 1e-2]          # m
        eps = keps * EPS0
        lde = math.sqrt(eps * Te / (n0 * QE))
        wp = math.sqrt(n0 * QE ** 2 / (eps * ME))
        vt = math.sqrt(2 * QE * Te / ME)
        dt = kt * min(h) / vt
        d = dict(h=h, lde=lde, wp=wp, vt=vt, dt=dt, wpdt=wp * dt, ns_heat=ns_heat_of(kt))
        if p['nu_h'] is not None and p['nu_h'] >= 0:
            d['nudt'] = p['nu_h'] * d['ns_heat'] * dt
        if p['nsteps']:
            d['t_end'] = p['nsteps'] * dt               # LePIC+ runs nsteps steps; tmax is not used
        if p['tmax']:
            d['nsteps_tmax'] = p['tmax'] / dt
        # skin depth from the TRUE eps0 (mod_state.f90 undoes k_eps0), RF update period
        d['skin'] = 2.99792458e8 / math.sqrt(n0 * QE ** 2 / (EPS0 * ME))
        if p['f0'] and p['f0'] > 0:
            d['ns_RF'] = 1.0 / (p['f0'] * dt)
        d['ncells'] = nx * ny * nz
        return d
    except Exception:
        return None


# --------------------------------------------------------------------------
#  Validation (the checks the Fortran code performs, plus sanity checks)
# --------------------------------------------------------------------------
def validate(case):
    p, issues = decode(case)
    box, segs, gi = decode_geometry(case)
    vals, bi = decode_boundary(case)
    issues += gi + bi
    issues += [i for i in case.cond.notes + case.geom.notes + case.bnd.notes]
    if case.uses_coil():
        issues += case.coil.notes
    A = issues.append
    c = case.cond

    def tab(t):
        return t

    # -------- reactions file
    rn = p['rname'] or ''
    if rn and (len(rn) > 20 or re.search(r'[\s,/\'"*]', rn)):
        A(Issue('error', 'Reactions file name must be <= 20 characters without spaces, commas, quotes or "/".',
                'rname', 'Simulation'))
    elif rn:
        d = os.path.join(case.folder, rn)
        d2 = os.path.join(os.path.dirname(case.folder), rn)
        if not (os.path.isfile(d) or os.path.isfile(d2)):
            A(Issue('warning', 'Reactions file "%s" was not found next to the input files or in the parent folder '
                    '(the code opens input_dir/%s).' % (rn, rn), 'rname', 'Simulation'))

    def need(key, cond_, msg, level='error', t='Simulation'):
        if p.get(key) is not None and not cond_(p[key]):
            A(Issue(level, msg, key, t))

    need('Te', lambda v: v > 0, 'Electron temperature must be positive.')
    need('n0', lambda v: v > 0, 'Plasma density must be positive.')
    need('ngas', lambda v: v >= 0, 'Gas density cannot be negative.')
    need('nsteps', lambda v: v > 0, 'The number of time steps must be > 0 (code: "conditions.inp line 1 must be: '
         '<reactions file> <number of time steps (> 0)>").')
    need('kt', lambda v: v > 0, 'The time-step coefficient must be positive.')
    need('kt', lambda v: v <= 1.0, 'Time-step coefficient > 1: an electron crosses more than one cell per '
         'step; results are likely inaccurate.', 'warning')
    need('eps', lambda v: v > 0, 'Solver precision must be positive.')
    need('omega', lambda v: 1.0 <= v < 2.0, 'The SOR parameter should be in [1, 2) (typically 1.4-1.8); '
         'otherwise the Poisson solver may diverge.', 'warning')
    need('np_cell', lambda v: v != 0, 'Particles per cell cannot be 0 (division by zero in the particle weight).')
    need('k_eps0', lambda v: v != 0, 'k_eps0 cannot be 0 ("Please correct k_eps0" - the code stops).')
    need('nsav', lambda v: v != 0, 'Saving frequency cannot be 0.')
    need('nbak', lambda v: v != 0, 'Backup frequency cannot be 0.', t='Output')
    need('cnt_plt', lambda v: v != 0, 'Number of snapshots = 0 causes a division by zero in the code. '
         'Use a negative number to turn snapshots off.', t='Output')

    # -------- grid / multigrid
    n3 = [p['nx'], p['ny'], p['nz']]
    ngl = p['ng']
    if None not in n3 and ngl is not None:
        for key, v in zip(('nx', 'ny', 'nz'), n3):
            if v < 2:
                A(Issue('error', 'Number of cells must be >= 2.', key, 'Simulation'))
        if all(v >= 1 for v in n3):
            # LePIC+ always uses the file value (no power-of-two override) and does not check it itself
            if ngl < 1 or ngl > 12:
                A(Issue('error', 'Number of multigrid levels must be between 1 and 12 (the solver has 12 '
                        'level arrays).', 'ng', 'Simulation'))
            else:
                bad = [nm for nm, v in zip(('nx', 'ny', 'nz'), n3) if v % (2 ** (ngl - 1)) != 0]
                for nm in bad:
                    A(Issue('error', 'With %d multigrid levels every grid size must be divisible by %d: the '
                            'coarse levels are n/2^(level-1). LePIC+ does not check this at start-up.'
                            % (ngl, 2 ** (ngl - 1)), nm, 'Simulation'))
        counts = [f for f in getattr(case, 'mpi', []) if f['n'] and f['n'] > 1]
        bad = [f for f in counts if n3[2] % f['n'] != 0]
        if bad:
            lst = ', '.join('%s (%d processes: %s)' % (f['script'], f['n'], f['how']) for f in bad)
            if len(bad) == len(counts):
                A(Issue('error', 'nz = %d is not divisible by the number of MPI processes set in %s. Each process '
                        'takes exactly nz/N planes of the Poisson domain (code: "nz not divisible by nproc").'
                        % (n3[2], lst), 'nz', 'Simulation'))
            else:
                A(Issue('warning', 'nz = %d is not divisible by the number of MPI processes set in %s. Each process '
                        'takes exactly nz/N planes of the Poisson domain (no remainder handling). Fine if you do '
                        'not use that script; the other scripts use a count that divides nz.' % (n3[2], lst),
                        'nz', 'Simulation'))
    for f in getattr(case, 'mpi', []):
        if f.get('note'):
            A(Issue('info', 'Run script %s: %s.' % (f['script'], f['note']), None, 'Summary'))

    # -------- geometry
    ok_box = None not in box
    if ok_box:
        for nm, b in zip(('xmax', 'ymax', 'zmax'), box):
            if b == 0:
                A(Issue('error', 'Box size cannot be 0.', 'g.' + nm, 'Geometry'))
    ngrid = p['ngrid']
    if ngrid is not None and ngrid < 1:
        A(Issue('error', 'Number of wall labels must be >= 1.', 'ngrid', 'Potentials'))
    if len(segs) > MAX_SEG:
        A(Issue('error', 'The code reads at most %d wall segments; the rest is ignored.' % MAX_SEG, None, 'Geometry'))
    labels_used = set()
    n_dielx = 0
    for k, s in enumerate(segs[:MAX_SEG]):
        if None in s:
            continue
        xl, yl, zl, xr, yr, zr, ind = s
        labels_used.add(abs(ind))
        if ngrid is not None and abs(ind) > ngrid:
            A(Issue('error', 'Segment %d uses label %d but only %d wall labels are defined '
                    '("Insufficient number of wall labels found in file boundary.inp").' % (k + 1, abs(ind), ngrid),
                    'g.seg%d.6' % k, 'Geometry'))
        for i, (v, nm) in enumerate(((xr, 'x>'), (yr, 'y>'), (zr, 'z>'))):
            t = case.geom.segs[k].get(3 + i)
            if is_neg(t) and v == 0:
                A(Issue('warning', 'Segment %d: "-0" cannot carry the dielectric/periodic flag (negative zero is '
                        'not < 0 in Fortran). Use a non-zero %s.' % (k + 1, nm), 'g.seg%d.%d' % (k, 3 + i), 'Geometry'))
        if is_neg(case.geom.segs[k].get(1)) and yl == 0:
            A(Issue('warning', 'Segment %d: "-0" cannot toggle cylinder mode.' % (k + 1), 'g.seg%d.1' % k, 'Geometry'))
        if xr < 0:
            n_dielx += 1
        if ok_box:
            if abs(xr) < xl or abs(yr) < abs(yl) or abs(zr) < zl:
                A(Issue('warning', 'Segment %d: an upper coordinate is smaller than the lower one - the segment '
                        'is empty.' % (k + 1), 'g.seg%d.3' % k, 'Geometry'))
    if n_dielx > 2:
        A(Issue('info', 'More than two x-plane dielectric segments: the code stores only two (type 3 and 4); '
                'extra ones are treated as type 4.', None, 'Geometry'))
    # secondary emission / e-beam need a Z-plane electrode
    gam, ig_sec, opt = p['gam_sec'], p['igrid_sec'], p['opt_inj']
    if None not in (gam, ig_sec, opt):
        if ngrid is not None and not (1 <= ig_sec <= ngrid):
            A(Issue('error', 'Emitting electrode label must be between 1 and the number of wall labels.',
                    'igrid_sec', 'Emission & RF'))
        if gam > 0 or abs(opt) == 4:
            for k, s in enumerate(segs[:MAX_SEG]):
                if None in s or abs(s[6]) != ig_sec:
                    continue
                if abs(s[5]) != abs(s[2]):
                    A(Issue('error', 'Secondary emission / emissive cathode: segment %d (label %d) must be a flat '
                            'plate with z< = z> (code: "secondary particle emission model only along (Oz)").'
                            % (k + 1, ig_sec), 'g.seg%d.5' % k, 'Geometry'))
            if ig_sec not in labels_used:
                A(Issue('warning', 'No geometry segment carries label %d (emitting electrode).' % ig_sec,
                        'igrid_sec', 'Emission & RF'))
    if p['I_inj'] is not None and gam is not None and p['I_inj'] > 0 and gam > 0 and not p['antenna']:
        A(Issue('error', 'Injected current > 0 together with secondary emission > 0 is not allowed '
                '(code: "gam_sec>0 is incompatible with I_inj>0").', 'I_inj', 'Heating & injection'))
    if opt is not None and abs(opt) not in (1, 2, 3, 4):
        A(Issue('error', 'Injection option must be +/-1, 2, 3 or 4.', 'opt_inj', 'Heating & injection'))

    # -------- heating / loading
    if ok_box:
        xm = abs(box[0])
        if p['x_load'] is not None and p['x_load'] > xm:
            A(Issue('info', 'Loading region is larger than the box: it is cut at xmax.', 'x_load', 'Heating & injection'))
        if None not in (p['xl_pow'], p['xr_pow']) and p['xr_pow'] <= p['xl_pow']:
            A(Issue('warning', 'Heating/injection region has x_max <= x_min.', 'xr_pow', 'Heating & injection'))
    d = derived(p, box) if ok_box else None
    if d:
        if d['wpdt'] > 2:
            A(Issue('error', 'wpe*dt = %.2f > 2: the explicit leap-frog scheme is unstable. Reduce the '
                    'time-step coefficient or the density.' % d['wpdt'], 'kt', 'Simulation'))
        elif d['wpdt'] > 0.5:
            A(Issue('warning', 'wpe*dt = %.2f: usually kept below ~0.2-0.5 for accuracy.' % d['wpdt'],
                    'kt', 'Simulation'))
        if 'nudt' in d and d['nudt'] > 1:
            A(Issue('error', 'nu_h * (heating interval) * dt = %.2f > 1 (code: "Warning, nu*dt>1" and stop). '
                    'Reduce the heating frequency.' % d['nudt'], 'nu_h', 'Heating & injection'))
    if p['nu_h'] is not None and p['Pabs'] is not None and not p['antenna']:
        if p['Pabs'] > 0 and p['nu_h'] < 0:
            A(Issue('info', 'Power is set but the heating frequency is negative: heating is OFF.', 'nu_h',
                    'Heating & injection'))
    if (p['delta_heat'] or 0) < 0 and (p['nu_h'] or -1) >= 0 and not p['antenna'] and (p['f0'] or 0) <= 0:
        A(Issue('error', 'Self-consistent skin depth (delta < 0) needs f0 > 0 on the RF line (update period 1/f0) '
                '- the code stops.', 'delta_heat', 'Heating & injection'))

    # -------- magnetic field
    nb = abs(p['nB']) if p['nB'] is not None else 0
    if nb > MAX_B:
        A(Issue('error', 'At most %d magnetic-field entries are supported.' % MAX_B, 'nB', 'Magnetic field'))
    for e in p['B']:
        k = e['k']
        if e['file'] not in ('y', 'n', 's'):
            A(Issue('error', 'Entry %d: first item must be y (map files), n (analytical) or s (skip).' % (k + 1),
                    'B%d.file' % k, 'Magnetic field'))
        if e['file'] == 'n':
            ch = (e['name'] or ' ')[0].lower()
            if ch not in 'gpet':
                A(Issue('warning', 'Entry %d: model name "%s" is not recognised (must start with g, P, e or t); '
                        'the entry has no effect.' % (k + 1, e['name']), 'B%d.name' % k, 'Magnetic field'))
            if ch in 'gp' and (e['dir'] or '').lower() not in ('bx', 'by', 'bz'):
                A(Issue('error', 'Entry %d: direction must be Bx, By or Bz.' % (k + 1), 'B%d.dir' % k,
                        'Magnetic field'))
            if ch == 'g' and e.get('dL') is not None and e['dL'] == 0:
                A(Issue('error', 'Entry %d: Gaussian width cannot be 0.' % (k + 1), 'B%d.dL' % k, 'Magnetic field'))
        if e['file'] == 'y':
            root = os.path.dirname(case.folder)
            for ax in 'xyz':
                fn = os.path.join(os.path.dirname(root), 'mag_field_files', e['name'] + 'B%s.dat' % ax)
                fn2 = os.path.join(root, 'mag_field_files', e['name'] + 'B%s.dat' % ax)
                if not (os.path.isfile(fn) or os.path.isfile(fn2)):
                    A(Issue('warning', 'Entry %d: map file mag_field_files/%sB%s.dat not found from here.'
                            % (k + 1, e['name'], ax), 'B%d.name' % k, 'Magnetic field'))
                    break
    if nb and None not in (p['nBx'], p['nBy'], p['nBz']):
        if min(p['nBx'], p['nBy'], p['nBz']) <= 1 and max(p['nBx'], p['nBy'], p['nBz']) > 1:
            A(Issue('info', 'If any B-map grid size is <= 1 the code resets all three to 1 (uniform field).',
                    'nBx', 'Magnetic field'))
    if nb and all(e['file'] == 's' for e in p['B']):
        A(Issue('info', 'All magnetic-field entries are "skip". (The code still allocates the B arrays because '
                'of a logic quirk in read_input: B_file!=s .or. B_file!=S is always true.)', 'nB', 'Magnetic field'))

    # -------- outputs / restart
    if None not in (p['t1'], p['t2'], p['cnt_plt']) and p['cnt_plt'] > 0 and p['t2'] <= p['t1']:
        A(Issue('warning', 'Snapshot window: end time <= start time.', 't2', 'Output'))
    n_rg = p['n_rg']
    if n_rg is not None and not (0 <= n_rg <= NM_RG):
        A(Issue('error', 'At most %d diagnostic regions are supported (nm_rg).' % NM_RG, 'n_rg', 'Output'))
    if p['restart'] == 2:
        if p['omp_max'] is not None and p['omp_max'] > 100:
            A(Issue('error', 'External restart: OpenMP thread count of the backup must be <= 100.', 'omp_max', 'Output'))
        if None not in (p['mpi_max'],) and p['mpi_max'] < 1:
            A(Issue('error', 'External restart: MPI rank count of the backup must be >= 1.', 'mpi_max', 'Output'))
    npd = p['np_dup']
    if npd is not None:
        if npd <= 0:
            A(Issue('error', 'Duplication/removal coefficient must be non-zero.', 'np_dup', 'Output'))
        elif npd > 1 and abs(npd - round(npd)) > 1e-12 and p['restart'] != 0:
            A(Issue('warning', 'For a coefficient > 1 use a whole number: the particle copies are made with '
                    'INT(coefficient) but the count is multiplied by the exact value (inconsistent otherwise).',
                    'np_dup', 'Output'))
    if p['restart'] == 0 and p['np_dup'] is not None and p['np_dup'] != 1:
        A(Issue('info', 'The duplication coefficient is ignored (set to 1) unless you restart.', 'np_dup', 'Output'))

    # -------- holes / negative ions / dielectric
    if p['holes']:
        for key in ('Lhy', 'Lhz'):
            need(key, lambda v: v > 0, 'Hole size must be positive.', t='Emission & RF')
        for key in ('nhy', 'nhz'):
            need(key, lambda v: v >= 1, 'Number of holes must be >= 1.', t='Emission & RF')
    if p['ind_g'] is not None and ngrid is not None and (p['holes'] or p['jne'] and p['jne'] > 0):
        if not (1 <= p['ind_g'] <= ngrid):
            A(Issue('error', 'Grid label must be between 1 and the number of wall labels.', 'ind_g', 'Emission & RF'))
        elif p['ind_g'] not in labels_used:
            A(Issue('warning', 'No geometry segment carries label %d (the grid to drill holes in).' % p['ind_g'],
                    'ind_g', 'Emission & RF'))
    if p['jne'] is not None and p['jne'] > 0:
        A(Issue('info', 'Negative-ion emission needs a species named [H-] in the reactions file.', 'jne',
                'Emission & RF'))
        if p['THm'] is not None and p['THm'] <= 0:
            A(Issue('error', 'Emission temperature must be positive.', 'THm', 'Emission & RF'))
    if p['THm'] is not None and p['THm'] <= 0 and (abs(opt or 0) in (3, 4) or (gam or 0) > 0):
        A(Issue('error', 'The emission temperature/energy must be > 0 for e-beam / emissive cathode / secondary '
                'emission.', 'THm', 'Emission & RF'))
    if p['Ca'] is not None and p['Ca'] < 0 and p['restart'] == 0:
        A(Issue('warning', 'A negative capacitance is only meaningful with a restart (it tells the code not to '
                'reload the dielectric potentials).', 'Ca', 'Emission & RF'))
    # -------- RF (mod_readConditions.f90 / mod_state.f90)
    mode = p['rf_mode']
    if mode is None:
        A(Issue('error', 'Unknown RF mode "%s" (expected No, Electrode, Wrapped, Planar or Coil) - the code stops.'
                % c.tok('rf'), 'rf', 'Emission & RF'))
    elif mode == 'Electrode':
        A(Issue('error', 'RF mode "Electrode" is not implemented in LePIC+ yet - the code stops.', 'rf',
                'Emission & RF'))
    if p['antenna']:
        A(Issue('info', 'RF antenna (%s): current injection and Maxwellian heating are switched off by the code; '
                'the drive amplitude is rescaled every RF period toward the target power.' % mode,
                'rf', 'Emission & RF'))
        for key, what in (('f0', 'f0 > 0'), ('E0', 'an initial %s > 0' % ('coil current I' if mode == 'Coil'
                                                                          else 'field E0')),
                          ('Pabs', 'a target power Pabs > 0')):
            if p[key] is not None and p[key] <= 0:
                A(Issue('error', 'The RF antenna needs %s - the code stops.' % what, key,
                        'Heating & injection' if key == 'Pabs' else 'Emission & RF'))
        # radius of the <ne> disk (Planar/Coil, R_win > 0) or of the wrapped coil (yr_pow - ymax/2)
        if ok_box and p['yr_pow'] is not None:
            R = p['R_win'] if (mode in ('Planar', 'Coil') and (p['R_win'] or 0) > 0) else \
                abs(p['yr_pow']) - abs(box[1]) / 2.0
            if R <= 0:
                A(Issue('error', 'Antenna radius <= 0: set R_win (Planar/Coil) or y_max of the heating region '
                        '> ymax/2 - the code stops.', 'R_win' if mode in ('Planar', 'Coil') else 'yr_pow',
                        'Emission & RF'))
    elif (p['R_win'] or 0) > 0:
        A(Issue('info', 'R_win is only used by the Planar and Coil antennas; ignored.', 'R_win', 'Emission & RF'))
    if mode == 'Coil':
        A(Issue('info', 'Coil mode: the "initial value" slot is the coil current amplitude I (A, peak); all '
                'species use the Boris push because of B_RF.', 'E0', 'Emission & RF'))
        cf = case.coil
        if not cf.exists:
            A(Issue('warning', 'coil.inp does not exist yet in this folder: the turns shown are defaults and '
                    'will be written on Save.', 'coil.n', 'Emission & RF'))
        nt = cf.count.get(0)
        if not is_int(nt) or int(nt) < 1:
            A(Issue('fatal', 'coil.inp: the number of turns must be an integer >= 1.', 'coil.n', 'Emission & RF'))
        elif len(cf.turns) < int(nt):
            A(Issue('error', 'coil.inp declares %s turns but has %d turn lines - the code stops.'
                    % (nt, len(cf.turns)), 'coil.n', 'Emission & RF'))
        x0 = max(0.0, p['xl_pow']) if p['xl_pow'] is not None else None
        for j, r in enumerate(cf.turns):
            v = [_num(r.get(i), 'coil%d.%d' % (j, i), 'float', issues, 'Coil turn %d: %s' % (j + 1, nm),
                      'Emission & RF') for i, nm in enumerate(('radius a', 'position x', 'weight w'))]
            if v[0] is not None and v[0] <= 0:
                A(Issue('error', 'Coil turn %d: radius must be > 0 - the code stops.' % (j + 1),
                        'coil%d.0' % j, 'Emission & RF'))
            if v[1] is not None and x0 is not None and v[1] >= x0:
                A(Issue('error', 'Coil turn %d is at x = %g cm, not behind the window x < %g cm (x_min of the '
                        'heating region, clamped to >= 0); the field is singular on a wire - the code stops.'
                        % (j + 1, v[1], x0), 'coil%d.1' % j, 'Emission & RF'))
            if v[0] is not None and ok_box and v[0] > max(abs(box[1]), abs(box[2])):
                A(Issue('warning', 'Coil turn %d: radius %g cm is larger than the box cross-section.'
                        % (j + 1, v[0]), 'coil%d.0' % j, 'Emission & RF'))
    if p['flag_convP']:
        A(Issue('info', 'k_eps0 < 0: in LePIC+ the only effect is that a restart reloads the electrode '
                'potentials from Backups/DATA.BAK/Vgrd.bak (the iterative power/current/voltage loop of the '
                'legacy code is not ported).', 'k_eps0', 'Simulation'))

    # -------- push scheme / Coulomb collisions
    ps = p['push'] or ''
    if ps not in ('momentum', 'energy'):
        A(Issue('warning', 'Push scheme "%s" is neither "momentum" nor "energy" (lower case): the code uses the '
                'momentum-conserving push.' % ps, 'push', 'Simulation'))
    if ps == 'energy' and ok_box and box[0] < 0:
        A(Issue('error', 'The energy-conserving push does not support the Neumann boundary at x = 0 (negative '
                'xmax) - the code stops.', 'push', 'Simulation'))
    if (c.tok('coulomb')[:1] or ' ') not in 'yYnN':
        A(Issue('warning', 'Coulomb collisions: anything not starting with y means "no".', 'coulomb', 'Simulation'))

    # extra tokens
    for nm, ex in case.cond.extras:
        A(Issue('info', 'Extra text "%s" after the values of the "%s" line is ignored by this version of the code.'
                % (' '.join(ex[:4]) + (' ...' if len(ex) > 4 else ''), HUMAN.get(nm, nm)), None, 'Summary'))
    # de-duplicate
    seen, out = set(), []
    for i in issues:
        k = (i.level, i.msg, i.key)
        if k not in seen:
            seen.add(k)
            out.append(i)
    out.sort(key=lambda i: Issue.ORDER[i.level])
    return out


# --------------------------------------------------------------------------
#  Geometry rasteriser: a faithful numpy port of the wall-drawing part of
#  generate_boundary()  (used only for the preview)
# --------------------------------------------------------------------------
def _nint(x):
    return int(math.floor(x + 0.5)) if x >= 0 else -int(math.floor(-x + 0.5))


def rasterize(box, segs, n, holes=None):
    """box: signed (xmax,ymax,zmax); segs: list of 7-tuples (signed, as in the file);
    n: (n1,n2,n3) grid cells; holes: dict(Lhy,Lhz,nhy,nhz,ind_g) or None.
    Returns dict(xy, xz, yz arrays of bcnd, dtype dict, L, notes)."""
    import numpy as np
    notes = []
    xs_, ys_, zs_ = box
    flag_nmn = xs_ < 0
    circx_tmp = ys_ < 0
    circz = zs_ < 0
    xmax, ymax, zmax = abs(xs_), abs(ys_), abs(zs_)
    n1, n2, n3 = n
    h1, h2, h3 = xmax / n1, ymax / n2, zmax / n3
    bc = np.zeros((n1 + 3, n2 + 3, n3 + 3), dtype=np.int16)
    dtype = {}
    cnt_yz = 0
    flag_pbc = False
    grid_info = None
    for (xl, yl, zl, xr, yr, zr, ind) in segs[:MAX_SEG]:
        a = abs(ind)
        dtype[a] = 1
        circx = circx_tmp
        if yl < 0:
            circx = not circx_tmp
        yl = abs(yl)
        R = 0.0
        if circx:
            R = min((yr - yl) / 2.0, (zr - zl) / 2.0)
        if circz:
            R = min((xr - xl) / 2.0, (yr - yl) / 2.0)
        if ind == 0:
            flag_pbc = True
        if xr < 0:
            xr = abs(xr)
            dtype[a] = 3 if cnt_yz == 0 else 4
            cnt_yz += 1
        if yr < 0:
            yr = abs(yr)
            dtype[a] = 2
        if zr < 0:
            zr = abs(zr)
            dtype[a] = -dtype[a]
            flag_pbc = True
        ixl = int(xl / h1) + 1
        if ixl <= 1:
            ixl = 0
        ixr = int(xr / h1) + 1
        if ixr >= n1:
            ixr = n1 + 2
        iyl = int(yl / h2) + 1
        if iyl <= 1:
            iyl = 0
        iyr = int(yr / h2) + 1
        if iyr >= n2:
            iyr = n2 + 2
        izl = int(zl / h3) + 1
        if izl <= 1:
            izl = 0
        izr = int(zr / h3) + 1
        if izr >= n3:
            izr = n3 + 2
        if holes and a == holes['ind_g']:
            grid_info = (ixl, ixr)
        if ind == 0:
            ixl = max(ixl, 2)
            ixr = min(ixr, n1)
        if ind >= 0:
            lo, hi = max(ixl, 0), min(ixr, n1 + 2)
            if hi >= lo:
                bc[lo:hi + 1, :, :] = ind
            ixl, ixr = max(ixl, 2), min(ixr, n1)
            iyl, iyr = max(iyl, 2), min(iyr, n2)
            izl, izr = max(izl, 2), min(izr, n3)
        if ixr < ixl or iyr < iyl or izr < izl:
            continue
        xs = (np.arange(ixl, ixr + 1) - 1) * h1
        ys = (np.arange(iyl, iyr + 1) - 1) * h2
        zs = (np.arange(izl, izr + 1) - 1) * h3
        sub = bc[ixl:ixr + 1, iyl:iyr + 1, izl:izr + 1]
        mask = np.ones(sub.shape, dtype=bool)
        if ind > 0 or ind < 0:
            if circz:
                mask &= ((xs[:, None, None] - xmax / 2) ** 2 + (ys[None, :, None] - ymax / 2) ** 2) <= R * R
            if circx:
                mask &= ((ys[None, :, None] - ymax / 2) ** 2 + (zs[None, None, :] - zmax / 2) ** 2) <= R * R
        sub[mask] = -1 if ind >= 0 else a
    if flag_nmn:
        bc[0:2, :, :] = -2
    if flag_pbc:
        negl = [k for k, v in dtype.items() if v < 0]
        if negl:
            m = np.isin(bc[:, :, 1], negl) & (bc[:, :, 2] == -1)
            ixs, iys = np.nonzero(m)
            for ix, iy in zip(ixs, iys):
                bc[ix, iy, 0:2] = 0
                bc[ix, iy, n3 + 1:n3 + 3] = 0
    dtype = {k: abs(v) for k, v in dtype.items()}
    ixg = None
    if holes and holes.get('on') and grid_info:
        igl, igr = grid_info
        Lhy, Lhz = holes['Lhy'], holes['Lhz']
        nhy, nhz = holes['nhy'], holes['nhz']
        iyh = _nint(Lhy / 2.0 / h2) + 1
        izh = _nint(Lhz / 2.0 / h3) + 1
        Lhy, Lhz = min(Lhy, ymax), min(Lhz, zmax)
        disk = (Lhy == Lhz)
        R = Lhy / 2.0
        ixg = _nint((igl + igr) / 2.0)
        igl, igr = max(igl, 0), min(igr, n1 + 2)
        for cz in range(1, nhz + 1):
            zd = -(nhz - 1) * 3.0 * Lhz / 4.0 + (cz - 1) * 3.0 * Lhz / 2.0 + zmax / 2.0
            izl_ = _nint(zd / h3) + 1 - izh
            izr_ = izl_ + 2 * izh
            if (izl_ < 1 or izr_ > n3 + 1) and nhz > 1:
                notes.append('Too many apertures along z: the code would stop.')
                break
            izl_, izr_ = max(izl_, 1), min(izr_, n3 + 1)
            for cy in range(1, nhy + 1):
                yd = -(nhy - 1) * 3.0 * Lhy / 4.0 + (cy - 1) * 3.0 * Lhy / 2.0 + ymax / 2.0
                iyl_ = _nint(yd / h2) + 1 - iyh
                iyr_ = iyl_ + 2 * iyh
                if (iyl_ < 1 or iyr_ > n2 + 1) and nhy > 1:
                    notes.append('Too many apertures along y: the code would stop.')
                    break
                iyl_, iyr_ = max(iyl_, 1), min(iyr_, n2 + 1)
                if izr_ < izl_ or iyr_ < iyl_:
                    continue
                iy = np.arange(iyl_, iyr_ + 1)
                iz = np.arange(izl_, izr_ + 1)
                yy = ((iy - 1) * h2)[:, None]
                zz = ((iz - 1) * h3)[None, :]
                val = -np.ones((len(iy), len(iz)), dtype=np.int16)
                edge = (iy[:, None] == 1) | (iy[:, None] == n2 + 1) | (iz[None, :] == 1) | (iz[None, :] == n3 + 1)
                val[edge] = 0
                m2 = ((yy - yd) ** 2 + (zz - zd) ** 2 <= R * R) if disk else np.ones(val.shape, dtype=bool)
                blk = bc[igl:igr + 1, iyl_:iyr_ + 1, izl_:izr_ + 1]
                blk[:, m2] = val[m2]
    ixm = (n1 // 2 + 1) if not ixg else ixg
    ixm = max(0, min(ixm, n1 + 2))
    return dict(
        xy=bc[1:n1 + 2, 1:n2 + 2, n3 // 2 + 1], xz=bc[1:n1 + 2, n2 // 2 + 1, 1:n3 + 2],
        yz=bc[ixm, 1:n2 + 2, 1:n3 + 2], vol=bc[1:n1 + 2, 1:n2 + 2, 1:n3 + 2], dtype=dtype,
        L=(xmax, ymax, zmax), notes=notes,
        pos=(n3 // 2 * h3, n2 // 2 * h2, (ixm - 1) * h1))



# --------------------------------------------------------------------------
#  3D view: orthographic voxel renderer in pure numpy (no OpenGL / matplotlib)
# --------------------------------------------------------------------------
def view_matrix(az_deg, el_deg):
    """Rows = camera right, up, towards-viewer axes (z of the box points up)."""
    import numpy as np
    el_deg = max(-89.0, min(89.0, el_deg))
    a, e = math.radians(az_deg), math.radians(el_deg)
    d = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
    xc = np.cross([0.0, 0.0, 1.0], d)
    xc /= np.linalg.norm(xc)
    yc = np.cross(d, xc)
    return np.array([xc, yc, d])


class View3D:
    """Projection used both by the renderer and by the GUI overlay (box edges, axes).
    Points are given in index space; `spacing` is the physical size of one cell along each axis."""

    def __init__(self, shape, size, az, el, zoom, spacing=(1.0, 1.0, 1.0)):
        import numpy as np
        self.np = np
        self.R = view_matrix(az, el)
        self.shape = np.array(shape, float)
        self.h = np.asarray(spacing, float)
        self.W, self.H = size
        diag = float(np.linalg.norm(self.shape * self.h))
        self.s = zoom * 0.9 * min(self.W, self.H) / diag          # pixels per physical unit
        self.ctr = (self.shape - 1) / 2.0 * self.h

    def project(self, pts):
        np = self.np
        pc = (np.asarray(pts, float) * self.h - self.ctr) @ self.R.T
        return self.W / 2 + self.s * pc[..., 0], self.H / 2 - self.s * pc[..., 1], pc[..., 2]


def render_voxels(vol, colors, size, az=-50.0, el=28.0, zoom=1.0, clip=None, quality=1.0,
                  bg=(244, 244, 240), max_samples=3.0e6, spacing=(1.0, 1.0, 1.0), plasma=None,
                  plasma_alpha=0.38):
    """vol: int array (nx,ny,nz) of wall labels (>0 walls, -1 plasma, others transparent).
    colors: {label: (r,g,b)} - wall labels not listed are hidden.
    clip: None or (axis, fraction, flip): keep voxels up to `fraction` along the axis (or from it if flip).
    spacing: physical size of a cell along x, y, z (cells need not be cubic).
    plasma: None, or (r,g,b) to draw the plasma region (label -1) as a translucent skin.
    Returns (ppm_bytes, View3D)."""
    import numpy as np
    W, H = size
    h = np.asarray(spacing, float)
    V = View3D(vol.shape, size, az, el, zoom, h)
    img = np.empty((H, W, 3), np.uint8)
    img[:] = bg
    vmax = int(vol.max()) if vol.size else 0
    lut = np.zeros(max(vmax, 0) + 2, bool)
    for k in colors:
        if 0 < k < len(lut):
            lut[k] = True
    vis = (vol > 0) & lut[np.clip(vol, 0, len(lut) - 1)]
    pl = (vol == -1) if plasma is not None else None
    if clip and clip[0] in (0, 1, 2):
        axis, frac, flip = clip
        ids = np.arange(vol.shape[axis]).reshape([-1 if a == axis else 1 for a in range(3)])
        cut = frac * (vol.shape[axis] - 1)
        keep = (ids >= cut) if flip else (ids <= cut)
        vis &= keep
        if pl is not None:
            pl &= keep
    cols = np.zeros((len(lut), 3), float)
    for k, c in colors.items():
        if 0 < k < len(lut):
            cols[k] = c
    L = np.array([0.35, 0.55, 0.75])
    L /= np.linalg.norm(L)

    def exposed_faces(mask):
        out = []
        for axis in range(3):
            for sg in (1, -1):
                nvec = np.zeros(3)
                nvec[axis] = sg
                nc = V.R @ nvec
                if nc[2] <= 1e-6:
                    continue                                  # back face
                nb = np.zeros_like(mask)
                src = [slice(None)] * 3
                dst = [slice(None)] * 3
                if sg > 0:
                    src[axis], dst[axis] = slice(1, None), slice(0, -1)
                else:
                    src[axis], dst[axis] = slice(0, -1), slice(1, None)
                nb[tuple(dst)] = mask[tuple(src)]
                ex = mask & ~nb
                idx = np.argwhere(ex)
                if len(idx):
                    out.append((axis, sg, nc, idx, ex))
        return out

    wf = exposed_faces(vis)
    pf = exposed_faces(pl) if pl is not None else []
    nf = sum(len(f[3]) for f in wf) + sum(len(f[3]) for f in pf)
    head = b'P6\n%d %d\n255\n' % (W, H)
    if not nf:
        return head + img.tobytes(), V
    m = max(2, int(math.ceil(V.s * float(h.max()) * 1.1 * quality)) + 1)
    m = max(2, min(m, int(math.sqrt(max_samples / nf))))
    uu = (np.arange(m) / (m - 1.0) - 0.5)
    U, Vv = np.meshgrid(uu, uu, indexing='ij')
    U, Vv = U.ravel(), Vv.ravel()

    def paint(faces, face_colors):
        """-> (flat pixel index, depth, colour) per sample, sorted far -> near."""
        allx, ally, alld, allc = [], [], [], []
        for axis, sg, nc, idx, ex in faces:
            a1, a2 = [a for a in range(3) if a != axis]
            e1, e2 = np.zeros(3), np.zeros(3)
            e1[a1], e2[a2] = h[a1], h[a2]
            r1, r2 = V.R @ e1, V.R @ e2
            offx = V.s * (r1[0] * U + r2[0] * Vv)
            offy = -V.s * (r1[1] * U + r2[1] * Vv)
            centre = idx.astype(float)
            centre[:, axis] += 0.5 * sg
            sx, sy, dp = V.project(centre)
            shade = 0.5 + 0.5 * max(0.0, float(nc @ L))
            allx.append(sx[:, None] + offx[None, :])
            ally.append(sy[:, None] + offy[None, :])
            alld.append(dp)
            allc.append(np.clip(face_colors(ex) * shade, 0, 255))
        depth = np.concatenate(alld)
        order = np.argsort(depth, kind='stable')
        X = np.concatenate(allx)[order]
        Y = np.concatenate(ally)[order]
        C = np.concatenate(allc)[order].astype(np.uint8)
        D = np.repeat(depth[order], m * m)
        ix = np.rint(X).astype(np.int64).ravel()
        iy = np.rint(Y).astype(np.int64).ravel()
        ok = (ix >= 0) & (ix < W) & (iy >= 0) & (iy < H)
        return (iy * W + ix)[ok], D[ok], np.repeat(C, m * m, axis=0)[ok]

    flat = img.reshape(-1, 3)
    zb = None
    if wf:
        fi, fd, fc = paint(wf, lambda ex: cols[vol[ex]])
        flat[fi] = fc                                      # last write = nearest
        if plasma is not None:
            zb = np.full(W * H, -np.inf)
            zb[fi] = fd
    if pf:
        if zb is None:
            zb = np.full(W * H, -np.inf)
        pc = np.array(plasma, float)
        fi, fd, fc = paint(pf, lambda ex: np.tile(pc, (int(ex.sum()), 1)))
        front = fd > zb[fi]                                # plasma skin only where nothing nearer
        fi, fc = fi[front], fc[front]
        if len(fi):
            layer = np.zeros((W * H, 3), np.float32)
            hit = np.zeros(W * H, bool)
            layer[fi] = fc
            hit[fi] = True
            a = float(plasma_alpha)
            flat[hit] = (flat[hit] * (1 - a) + layer[hit] * a).astype(np.uint8)
    return head + img.tobytes(), V


# --------------------------------------------------------------------------
#  Self test:  python pic_inputs.py  <folder> [<folder> ...]
# --------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    for fol in sys.argv[1:]:
        print('=' * 70)
        print(fol)
        try:
            cs = Case(fol)
        except Exception as e:   # noqa
            print('  cannot load:', e)
            continue
        base = {os.path.basename(f.path): f.render() for f in cs.files}
        print('  round trip (before any edit): ',
              {k: (v == f.original) for (k, v), f in zip(base.items(), cs.files)})
        for i in validate(cs):
            print('  ', i)
