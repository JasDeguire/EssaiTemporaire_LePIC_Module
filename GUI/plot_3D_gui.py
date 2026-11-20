#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_3D_gui.py  --  Viewer for Output/phi_n_3D.dat (LePIC+) and DATA/3D_param.dat (3D PIC code)

Usage
-----
    python3 GUI/plot_3D_gui.py                                  # looks for Output/phi_n_3D.dat
    python3 GUI/plot_3D_gui.py Output/phi_n_3D.dat input_dir/conditions.inp
    python3 GUI/plot_3D_gui.py DATA/3D_param_012.dat

Adapted from the LePIC3D_DD 3.4.1 viewer (plot_3D_param_gui.py).
Requirements: Python 3, numpy, matplotlib (TkAgg backend) and tkinter.

LePIC+ file  (modules/simulation/mod_simulation.f90, Fortran unformatted sequential)
-------------------------------------------------------------------------------------
  record 1 : nx,ny,nz (int32)
  record 2 : map(2,nx+1,ny+1,nz+1) (real8) = time-averaged phi [V] and density [m^-3] of species
             min(2, number of charged species) - the first ion species (avg3D(4:5) of mod_state.f90).
  Written at every save when the 3D averages are on (number of wall labels > 0 in conditions.inp);
  each save overwrites the file.  The cell size is not in the file: it is taken from geometry.inp
  (next to the conditions.inp given below), dx = |xmax|/nx ...
  E-field: computed from phi as below, masked where the density > 0.

3D code file  (Src/main.f90, Fortran unformatted sequential)
-------------------------------------------------------------
Mouse
-----
  hover        : value of the nearest grid node
  right click  : save figure as PNG / PDF / EPS
  3D tab       : left drag = rotate, wheel = zoom, middle drag = pan
"""
import os
import sys
import re
import math
import time
import struct
import atexit
import shutil
import tempfile
import traceback

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import ScalarMappable
from matplotlib.text import Annotation
from matplotlib.lines import Line2D
from matplotlib.transforms import IdentityTransform
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, colorchooser

QE = 1.6e-19          # same value as convert.f90
AMU = 1.66053886e-27  # as in Src/main.f90

# --------------------------------------------------------------------------------------
# Publication style
# --------------------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif", "Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "axes.linewidth": 0.9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 5, "ytick.major.size": 5,
    "xtick.minor.size": 2.5, "ytick.minor.size": 2.5,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "pdf.fonttype": 42, "ps.fonttype": 42,      # embed TrueType (editable text)
    "savefig.dpi": 300,
    "figure.facecolor": "white",
})

UNITS = {"mm": 1.0e3, "cm": 1.0e2, "m": 1.0}
AXN = "xyz"
PLASMA10 = "plasma (10 colors)"      # qualitative (banded) versions: 10 discrete colours
INFERNO10 = "inferno (10 colors)"
DISCRETE = {PLASMA10: "plasma", INFERNO10: "inferno"}
CMAPS = ["viridis", "plasma", PLASMA10, "inferno", INFERNO10, "magma", "cividis", "turbo", "jet",
         "coolwarm", "RdBu_r", "seismic", "gray"]
RES = {"1/8 of grid (default)": 8, "1/2 of grid": 2, "Full grid": 1}


# --------------------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------------------
class Param3D:
    """Lazy (memory-mapped) reader of phi_n_3D.dat (LePIC+) and 3D_param*.dat."""

    def __init__(self, path, geom=None):
        self.path = path
        self._tmp = None
        self.lepic = False
        self.dens_ptype = 1                       # species whose density is stored (option 1)
        size = os.path.getsize(path)
        last = None
        self._try_lepic(path, size, geom)
        for mk in (() if self.lepic else (4, 8)):
            for en in ("<", ">"):
                try:
                    if self._try(path, size, mk, en):
                        break
                except Exception as e:           # noqa
                    last = e
            else:
                continue
            break
        else:
            if not self.lepic:
                raise ValueError("Not a recognised phi_n_3D / 3D_param file (record structure does not "
                                 "match the expected layout)." + (f"\n{last}" if last else ""))
        n = (self.nx + 1, self.ny + 1, self.nz + 1)
        shape = (self.navg,) + n
        if len(self.chunks) == 1:
            self.arr = np.memmap(path, dtype=self.en + "f8", mode="r",
                                 offset=self.chunks[0][0], shape=shape, order="F")
        else:
            self.arr = self._assemble(shape)
        self.n = (self.nx, self.ny, self.nz)
        self.h = (self.dx, self.dy, self.dz)
        self.L = tuple(self.n[a] * self.h[a] for a in range(3))
        # option detection
        if self.navg == 2:
            self.option, self.hasB = 1, False
        elif self.navg == 1 + 7 * self.ntype:
            self.option, self.hasB = 2, False
        elif self.navg == 4 + 7 * self.ntype:
            self.option, self.hasB = 2, True
        else:
            raise ValueError(f"Unexpected number of saved parameters navg3D={self.navg} "
                             f"for ntype={self.ntype}.")
        self._cache = {}

    # ---- low level -----------------------------------------------------------------
    def _try_lepic(self, path, size, geom):
        """LePIC+ phi_n_3D.dat: [nx,ny,nz] then map(2,nx+1,ny+1,nz+1), 4-byte markers (gfortran
        splits records > 2 GB into sub-records with negative markers)."""
        for en in ("<", ">"):
            with open(path, "rb") as f:
                b = f.read(20)
                if len(b) < 20 or struct.unpack(en + "i", b[:4])[0] != 12 or \
                        struct.unpack(en + "i", b[16:20])[0] != 12:
                    continue
                nx, ny, nz = struct.unpack(en + "3i", b[4:16])
                if not all(1 <= v < 100000 for v in (nx, ny, nz)):
                    continue
                expected = 16 * (nx + 1) * (ny + 1) * (nz + 1)
                pos, chunks, tot = 20, [], 0
                while tot < expected:
                    f.seek(pos)
                    bb = f.read(4)
                    if len(bb) < 4:
                        break
                    L = abs(struct.unpack(en + "i", bb)[0])
                    if L == 0 or tot + L > expected:
                        break
                    chunks.append((pos + 4, L))
                    tot += L
                    pos += 8 + L
                if tot != expected or pos > size:
                    continue
            self.nx, self.ny, self.nz, self.navg, self.ntype = nx, ny, nz, 2, 1
            self.mass, self.en, self.chunks, self.lepic = [], en, chunks, True
            self.geom_note = None
            L = None
            if geom and os.path.isfile(geom):
                try:
                    with open(geom, errors="replace") as g:
                        for ln in g:
                            t = ln.split("!")[0].replace(",", " ").split()
                            if t:
                                L = [abs(float(v.replace("d", "e").replace("D", "e"))) * 1e-2 for v in t[:3]]
                                break
                except Exception:          # noqa
                    L = None
            if not L or len(L) < 3:
                L = [nx * 1e-3, ny * 1e-3, nz * 1e-3]
                self.geom_note = "geometry.inp not found: 1 mm cells assumed"
            self.dx, self.dy, self.dz = L[0] / nx, L[1] / ny, L[2] / nz
            po = os.path.join(os.path.dirname(os.path.abspath(path)), "particles.out")
            self.names = parse_particles_out(po)
            ncharged = 1 + sum(1 for nm in self.names[1:] if nm.rstrip("]").endswith(("+", "-")))
            self.dens_ptype = min(2, max(1, ncharged))
            return True
        return False

    def _try(self, path, size, mk, en):
        mi = "i" if mk == 4 else "q"
        with open(path, "rb") as f:
            b = f.read(mk + 20)
            if len(b) < mk + 20:
                return False
            m = struct.unpack(en + mi, b[:mk])[0]
            nx, ny, nz, navg, ntype = struct.unpack(en + "5i", b[mk:mk + 20])
            if not (1 <= ntype <= 64 and 1 <= nx < 100000 and 1 <= ny < 100000
                    and 1 <= nz < 100000 and 1 <= navg < 1000):
                return False
            if m != 20 + 8 * (3 + ntype):
                return False
            f.seek(mk + 20)
            dd = struct.unpack(en + "%dd" % (3 + ntype), f.read(8 * (3 + ntype)))
            if struct.unpack(en + mi, f.read(mk))[0] != m:
                return False
            pos = 2 * mk + m
            expected = 8 * navg * (nx + 1) * (ny + 1) * (nz + 1)
            chunks, tot = [], 0
            while tot < expected:      # follow (sub)records
                f.seek(pos)
                bb = f.read(mk)
                if len(bb) < mk:
                    return False
                m2 = struct.unpack(en + mi, bb)[0]
                L = abs(m2)
                if L == 0 or tot + L > expected:
                    return False
                f.seek(pos + mk + L)
                bb = f.read(mk)
                if len(bb) < mk or struct.unpack(en + mi, bb)[0] != m2:
                    return False
                chunks.append((pos + mk, L))
                tot += L
                pos += 2 * mk + L
            if pos > size:
                return False
        self.nx, self.ny, self.nz, self.navg, self.ntype = nx, ny, nz, navg, ntype
        self.dx, self.dy, self.dz = dd[:3]
        self.mass = list(dd[3:])
        self.en, self.chunks = en, chunks
        return True

    def _assemble(self, shape):
        """Record split in sub-records (>2 GB): build one contiguous temporary file."""
        tmpdir = tempfile.mkdtemp(prefix="p3d_")
        atexit.register(shutil.rmtree, tmpdir, True)
        tmp = os.path.join(tmpdir, "assembled.raw")
        with open(self.path, "rb") as fi, open(tmp, "wb") as fo:
            for off, L in self.chunks:
                fi.seek(off)
                left = L
                while left > 0:
                    blk = fi.read(min(left, 1 << 26))
                    fo.write(blk)
                    left -= len(blk)
        self._tmp = tmp
        return np.memmap(tmp, dtype=self.en + "f8", mode="r", shape=shape, order="F")

    # ---- derived quantities ----------------------------------------------------------
    def plane(self, var, axis, idx, ids, zero_bnd):
        """2D plane of `var` at full-grid index idx along `axis`.
        Returns array [i1,i2] over the two remaining axes (ascending), sampled at ids[a]."""
        key = (var["key"], axis, idx, tuple(len(i) for i in ids), zero_bnd)
        if key in self._cache:
            return self._cache[key]
        rem = [a for a in range(3) if a != axis]
        sel = [None] * 3
        sel[axis] = int(idx)
        sel[rem[0]] = ids[rem[0]][:, None]
        sel[rem[1]] = ids[rem[1]][None, :]
        sel = tuple(sel)
        A = self.arr
        kind, p = var["kind"], var["p"]
        if kind == "phi":
            d = np.array(A[(0,) + sel], dtype=float)
        elif kind == "E":
            d = self.efield(var["c"], axis, idx, ids)
        else:
            base = 7 * p - 6                      # zero-based index of n for ptype p
            n = np.array(A[(base,) + sel], dtype=float)
            if kind == "n":
                d = n
            elif kind == "flux":
                v = np.array(A[(base + 1 + var["c"],) + sel], dtype=float)
                d = n * v
            elif kind == "fnorm":      # |n<v>| = n * sqrt(<vx>^2+<vy>^2+<vz>^2)
                v = np.array(A[(slice(base + 1, base + 4),) + sel], dtype=float)
                d = n * np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
            elif kind == "T":
                v = np.array(A[(slice(base + 1, base + 7),) + sel], dtype=float)
                tp = self.mass[p - 1] / (3.0 * QE) * (
                    v[3] + v[4] + v[5] - v[0] ** 2 - v[1] ** 2 - v[2] ** 2)
                d = np.where(tp > 0.0, tp, 0.0)
            else:
                raise KeyError(kind)
            if zero_bnd:     # convert.f90 only fills interior nodes 2..n (1-based)
                mask = np.zeros(d.shape, bool)
                if idx == 0 or idx == self.n[axis]:
                    mask[:] = True
                for k, a in enumerate(rem):
                    b = (ids[a] == 0) | (ids[a] == self.n[a])
                    mask |= b[:, None] if k == 0 else b[None, :]
                d = np.where(mask, 0.0, d)
        if len(self._cache) > 80:
            self._cache.clear()
        self._cache[key] = d
        return d


def _efield(self, c, axis, idx, ids):
    """E-field component c (0,1,2 = x,y,z) on a plane, exactly as in convert.f90:
         E_x(i,j,k) = -( phi(i+1,j,k) - phi(i-1,j,k) ) / (2 dx)   (same for y, z)
       evaluated on the full-resolution grid at interior nodes only (Fortran 2..n), and only
       where param_3D(2,...) = n(ptype 1) > 0; E = 0 everywhere else."""
    rem = [a for a in range(3) if a != axis]
    shp = (len(ids[rem[0]]), len(ids[rem[1]]))
    # interior-node mask
    inter = np.zeros(shp, bool) if not (1 <= idx <= self.n[axis] - 1) else np.ones(shp, bool)
    for k, a in enumerate(rem):
        ok = (ids[a] >= 1) & (ids[a] <= self.n[a] - 1)
        inter &= ok[:, None] if k == 0 else ok[None, :]
    if not inter.any():
        return np.zeros(shp)

    def sel(shift):
        t = [None] * 3
        t[axis] = int(np.clip(idx + (shift if c == axis else 0), 0, self.n[axis]))
        for k, a in enumerate(rem):
            i = np.clip(ids[a] + (shift if c == a else 0), 0, self.n[a])
            t[a] = i[:, None] if k == 0 else i[None, :]
        return tuple(t)
    A = self.arr
    ph_p = np.array(A[(0,) + sel(+1)], dtype=float)
    ph_m = np.array(A[(0,) + sel(-1)], dtype=float)
    n1 = np.array(A[(1,) + sel(0)], dtype=float)         # param_3D(2,...) = n of ptype 1
    E = -(ph_p - ph_m) / (2.0 * self.h[c])
    return np.where(inter & (n1 > 0.0), E, 0.0)


Param3D.efield = _efield


def _flux_vec(self, p, axis, idx, ids, zero_bnd):
    """Flux vector n<v> of ptype p on a plane: list of the 3 Cartesian components [i1,i2]."""
    rem = [a for a in range(3) if a != axis]
    sel = [None] * 3
    sel[axis] = int(idx)
    sel[rem[0]] = ids[rem[0]][:, None]
    sel[rem[1]] = ids[rem[1]][None, :]
    sel = tuple(sel)
    base = 7 * p - 6
    n = np.array(self.arr[(base,) + sel], dtype=float)
    v = np.array(self.arr[(slice(base + 1, base + 4),) + sel], dtype=float)
    G = [n * v[c] for c in range(3)]
    if zero_bnd:
        mask = np.zeros(n.shape, bool)
        if idx == 0 or idx == self.n[axis]:
            mask[:] = True
        for k, a in enumerate(rem):
            b = (ids[a] == 0) | (ids[a] == self.n[a])
            mask |= b[:, None] if k == 0 else b[None, :]
        G = [np.where(mask, 0.0, g) for g in G]
    return G


Param3D.flux_vec = _flux_vec


def parse_particles_out(path):
    """Species names in ptype order from Output/particles.out ([e] = ptype 1)."""
    names, sec = [], False
    try:
        with open(path, errors="replace") as f:
            for ln in f:
                s = ln.strip()
                if not s:
                    sec = False
                    continue
                if s.endswith(":"):
                    sec = s[:-1].strip().lower() in ("ions", "neutrals")
                    continue
                if sec:
                    names.append(s.split()[0])
    except OSError:
        pass
    return ["[e]"] + names


def parse_ngrid(path):
    """Number of wall labels from LePIC+ conditions.inp, found through its comment (None if absent)."""
    try:
        with open(path, errors="replace") as f:
            for ln in f:
                if "!" in ln and "wall labels" in ln.split("!", 1)[1].lower():
                    return int(ln.split("!")[0].split()[0])
    except Exception:          # noqa
        pass
    return None


def parse_conditions(path):
    """Return dict(save=bool, opt=int) from the 'Save 3D data' line of conditions.inp."""
    try:
        with open(path, errors="replace") as f:
            for ln in f:
                if re.search(r"save\s*3d", ln, re.I):
                    tok = ln.split("!")[0].split()
                    if len(tok) >= 2:
                        return dict(save=tok[0].lower().startswith("y"), opt=int(tok[1]))
                    if tok:
                        return dict(save=tok[0].lower().startswith("y"), opt=None)
    except Exception:
        pass
    return None


# --------------------------------------------------------------------------------------
# GUI helpers
# --------------------------------------------------------------------------------------
class MiniToolbar(NavigationToolbar2Tk):
    toolitems = [t for t in NavigationToolbar2Tk.toolitems
                 if t[0] in ("Home", "Back", "Forward", "Pan", "Zoom")]


class Tab:
    """One figure (canvas, hover tooltip with blitting, right-click save menu)."""

    def __init__(self, app, parent, kind):
        self.app, self.kind = app, kind
        self.frame = ttk.Frame(parent)
        self.fig = Figure(figsize=(6.6, 5.0), dpi=100, facecolor="white")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.toolbar = None
        if kind != "3D":
            self.toolbar = MiniToolbar(self.canvas, self.frame)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.dirty = True
        self.bg = None
        self.ax = None
        self.hover = None          # callable(event) -> (px, py, text) | None
        self.ann = None
        self.mark = None
        self._t_last = 0.0
        self._press = None
        c = self.canvas
        c.mpl_connect("draw_event", self._on_draw)
        c.mpl_connect("motion_notify_event", self._on_move)
        c.mpl_connect("button_press_event", self._on_press)
        c.mpl_connect("button_release_event", self._on_release)
        c.mpl_connect("figure_leave_event", lambda e: self._hide(True))
        self.menu = tk.Menu(self.frame, tearoff=0)
        for fmt in ("png", "pdf", "eps"):
            self.menu.add_command(label=f"Save figure as {fmt.upper()}…",
                                  command=lambda f=fmt: app.save_figure(self, f))
        self.menu.add_separator()
        self.menu.add_command(label="Reset view", command=lambda: app.reset_view(self))

    # ---- tooltip artists ---------------------------------------------------------------
    def new_artists(self):
        self.bg = None
        self.ann = Annotation("", xy=(0, 0), xycoords="figure pixels", xytext=(14, 14),
                              textcoords="offset points", fontsize=9, zorder=100,
                              ha="left", va="bottom", animated=True,
                              bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.25",
                                        alpha=0.96))
        self.mark = Line2D([0], [0], marker="o", ms=7, mfc="none", mec="red", mew=1.3,
                           ls="", transform=IdentityTransform(), animated=True, zorder=101)
        for a in (self.ann, self.mark):
            self.fig.add_artist(a)
            a.set_visible(False)

    def _on_draw(self, ev):
        if self.ann is None or ev.canvas is not self.canvas:
            return          # ignore draw events of temporary save canvases
        self.bg = self.canvas.copy_from_bbox(self.fig.bbox)
        self._draw_anim()

    def _draw_anim(self):
        for a in (self.mark, self.ann):
            if a is not None:
                self.fig.draw_artist(a)

    def _blit(self):
        if self.bg is None:
            return
        self.canvas.restore_region(self.bg)
        self._draw_anim()
        self.canvas.blit(self.fig.bbox)

    def _hide(self, update=False):
        if self.ann is not None and (self.ann.get_visible() or self.mark.get_visible()):
            self.ann.set_visible(False)
            self.mark.set_visible(False)
            if update:
                self._blit()

    def _on_move(self, ev):
        if self.hover is None or self.bg is None:
            return
        now = time.time()
        if now - self._t_last < 0.03:
            return
        self._t_last = now
        if ev.button is not None or ev.x is None:
            self._hide(True)
            return
        try:
            r = self.hover(ev)
        except Exception:
            r = None
        if r is None:
            self._hide(True)
            return
        px, py, text = r
        W, H = self.fig.canvas.get_width_height()
        # flip the box near the right / top edges
        left = px < 0.62 * W
        up = py < 0.78 * H
        self.ann.set_text(text)
        self.ann.xy = (px, py)
        self.ann.set_position((14 if left else -14, 14 if up else -14))
        self.ann.set_ha("left" if left else "right")
        self.ann.set_va("bottom" if up else "top")
        self.mark.set_data([px], [py])
        self.ann.set_visible(True)
        self.mark.set_visible(True)
        self.app.status.set(text.replace("\n", "   "))
        self._blit()

    # ---- right click -------------------------------------------------------------------
    def _on_press(self, ev):
        if ev.button == 3:
            self._press = (ev.x, ev.y)

    def _on_release(self, ev):
        if ev.button == 3 and self._press is not None:
            dx, dy = ev.x - self._press[0], ev.y - self._press[1]
            self._press = None
            if dx * dx + dy * dy < 36 and ev.guiEvent is not None:
                self._hide(True)
                try:
                    self.menu.tk_popup(ev.guiEvent.x_root, ev.guiEvent.y_root)
                finally:
                    self.menu.grab_release()


# --------------------------------------------------------------------------------------
# Main application
# --------------------------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self, data_path=None, cond_path=None):
        super().__init__()
        self.title("LePIC+ 3D viewer")
        self.geometry("1420x880")
        self.minsize(1100, 700)
        self.data = None
        self.vars = []
        self.pos = [0.0, 0.0, 0.0]          # slice coordinates [m]
        self._lock = False
        self._job = None
        self.zoom3d = 1.0
        self.view3d = (25.0, 45.0, 0.0)
        self.status = tk.StringVar(value="Open Output/phi_n_3D.dat (or a 3D_param.dat file).")
        self.ct_store, self._ct_key = {}, None      # contour levels remembered per quantity
        self.ct_unit = tk.StringVar(value="V")
        self._build_controls()
        self._build_tabs()
        ttk.Label(self, textvariable=self.status, anchor="w", relief="sunken").pack(
            side=tk.BOTTOM, fill=tk.X)
        self.left.pack(side=tk.LEFT, fill=tk.Y)
        self.nb.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        # defaults
        if data_path is None:
            for c in ("Output/phi_n_3D.dat", "../Output/phi_n_3D.dat", "phi_n_3D.dat",
                      "DATA/3D_param.dat", "../../DATA/3D_param.dat", "3D_param.dat"):
                if os.path.exists(c):
                    data_path = c
                    break
        if cond_path is None:
            base = os.path.dirname(os.path.abspath(data_path)) if data_path else os.getcwd()
            for c in (os.path.join(base, "..", "input_dir", "conditions.inp"),
                      os.path.join(base, "..", "..", "input_dir", "conditions.inp"),
                      "input_dir/conditions.inp", "conditions.inp"):
                if os.path.exists(c):
                    cond_path = os.path.normpath(c)
                    break
        if data_path:
            self.f_var.set(data_path)
        if cond_path:
            self.c_var.set(cond_path)
        if data_path and os.path.exists(data_path):
            self.after(100, self.load)

    # ---- controls -----------------------------------------------------------------------
    def _build_controls(self):
        outer = ttk.Frame(self, width=330)
        self.left = outer
        cv = tk.Canvas(outer, width=330, highlightthickness=0)
        sb = ttk.Scrollbar(outer, orient="vertical", command=cv.yview)
        cv.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        cv.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        f = ttk.Frame(cv)
        cv.create_window((0, 0), window=f, anchor="nw")
        f.bind("<Configure>", lambda e: cv.configure(scrollregion=cv.bbox("all")))

        def wheel(e):
            cv.yview_scroll(-1 if (e.num == 4 or e.delta > 0) else 1, "units")
        f.bind("<Enter>", lambda e: (cv.bind_all("<MouseWheel>", wheel),
                                     cv.bind_all("<Button-4>", wheel),
                                     cv.bind_all("<Button-5>", wheel)))
        f.bind("<Leave>", lambda e: (cv.unbind_all("<MouseWheel>"),
                                     cv.unbind_all("<Button-4>"), cv.unbind_all("<Button-5>")))

        # files
        lf = ttk.LabelFrame(f, text="Files")
        lf.pack(fill=tk.X, padx=6, pady=4)
        self.f_var, self.c_var = tk.StringVar(), tk.StringVar()
        for r, (lab, var, cmd) in enumerate((("phi_n_3D.dat / 3D_param:", self.f_var, self.browse_data),
                                              ("conditions.inp:", self.c_var, self.browse_cond))):
            ttk.Label(lf, text=lab).grid(row=2 * r, column=0, columnspan=2, sticky="w", padx=4)
            ttk.Entry(lf, textvariable=var, width=34).grid(row=2 * r + 1, column=0, padx=4, pady=1)
            ttk.Button(lf, text="…", width=3, command=cmd).grid(row=2 * r + 1, column=1)
        ttk.Button(lf, text="Load", command=self.load).grid(row=4, column=0, columnspan=2,
                                                            pady=4)
        self.info = tk.Text(f, width=40, height=9, wrap="word", font=("TkFixedFont", 9),
                            relief="flat", bg=self.cget("bg"))
        self.info.pack(fill=tk.X, padx=8)
        self.info.configure(state="disabled")

        # variable
        lf = ttk.LabelFrame(f, text="Quantity")
        lf.pack(fill=tk.X, padx=6, pady=4)
        self.var_cb = ttk.Combobox(lf, state="readonly", width=34)
        self.var_cb.pack(padx=4, pady=3)
        self.var_cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        r = ttk.Frame(lf)
        r.pack(fill=tk.X, padx=4, pady=2)
        ttk.Label(r, text="Species names:").pack(side=tk.LEFT)
        self.sp_var = tk.StringVar()
        e = ttk.Entry(r, textvariable=self.sp_var, width=16)
        e.pack(side=tk.LEFT, padx=4)
        e.bind("<Return>", lambda ev: self.rebuild_vars())
        e.bind("<FocusOut>", lambda ev: self.rebuild_vars())
        ttk.Label(lf, text="e.g.  e, Xe+, Xe   (one per ptype, in order)",
                  foreground="gray40").pack(anchor="w", padx=4)

        # grid
        lf = ttk.LabelFrame(f, text="Grid resolution (per direction)")
        lf.pack(fill=tk.X, padx=6, pady=4)
        self.res_var = tk.StringVar(value=list(RES)[0])
        for k in RES:
            ttk.Radiobutton(lf, text=k, value=k, variable=self.res_var,
                            command=self.refresh).pack(anchor="w", padx=6)
        r = ttk.Frame(lf)
        r.pack(fill=tk.X, padx=6, pady=2)
        self.cap_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(r, text="3D view: at most", variable=self.cap_var,
                        command=self.refresh).pack(side=tk.LEFT)
        self.cap_n = tk.StringVar(value="96")
        sp = ttk.Spinbox(r, from_=24, to=512, increment=8, width=5, textvariable=self.cap_n,
                         command=self.refresh)
        sp.pack(side=tk.LEFT, padx=3)
        sp.bind("<Return>", lambda e: self.refresh())
        ttk.Label(r, text="nodes/axis").pack(side=tk.LEFT)

        # planes
        lf = ttk.LabelFrame(f, text="Plane coordinates")
        lf.pack(fill=tk.X, padx=6, pady=4)
        self.unit_var = tk.StringVar(value="cm")
        r = ttk.Frame(lf)
        r.pack(fill=tk.X, padx=4)
        ttk.Label(r, text="Length unit:").pack(side=tk.LEFT)
        cb = ttk.Combobox(r, textvariable=self.unit_var, values=list(UNITS), width=5,
                          state="readonly")
        cb.pack(side=tk.LEFT, padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self.on_unit())
        self.sl, self.sl_ent, self.sl_var, self.node_lbl = [], [], [], []
        names = ("x  (YZ plane)", "y  (XZ plane)", "z  (XY plane)")
        for a in range(3):
            ttk.Label(lf, text=names[a]).pack(anchor="w", padx=4)
            r = ttk.Frame(lf)
            r.pack(fill=tk.X, padx=4)
            s = ttk.Scale(r, from_=0, to=1, length=190,
                          command=lambda v, a=a: self.on_slider(a, float(v)))
            s.pack(side=tk.LEFT)
            v = tk.StringVar()
            e = ttk.Entry(r, textvariable=v, width=8)
            e.pack(side=tk.LEFT, padx=4)
            e.bind("<Return>", lambda ev, a=a: self.on_entry(a))
            e.bind("<FocusOut>", lambda ev, a=a: self.on_entry(a))
            l = ttk.Label(lf, text="", foreground="gray40")
            l.pack(anchor="w", padx=6)
            self.sl.append(s); self.sl_var.append(v); self.sl_ent.append(e); self.node_lbl.append(l)

        # display
        lf = ttk.LabelFrame(f, text="Display")
        lf.pack(fill=tk.X, padx=6, pady=4)
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=1)
        ttk.Label(r, text="Colormap:").pack(side=tk.LEFT)
        self.cmap_var = tk.StringVar(value="viridis")
        cb = ttk.Combobox(r, textvariable=self.cmap_var, values=CMAPS, width=10,
                          state="readonly")
        cb.pack(side=tk.LEFT, padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=1)
        ttk.Label(r, text="Shading:").pack(side=tk.LEFT)
        self.shade_var = tk.StringVar(value="Nearest node")
        cb = ttk.Combobox(r, textvariable=self.shade_var, width=14, state="readonly",
                          values=["Nearest node", "Smooth (Gouraud)"])
        cb.pack(side=tk.LEFT, padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=1)
        ttk.Label(r, text="Colour scale  min:").pack(side=tk.LEFT)
        self.vmin_var, self.vmax_var = tk.StringVar(), tk.StringVar()
        for v in (self.vmin_var, self.vmax_var):
            e = ttk.Entry(r, textvariable=v, width=9)
            e.pack(side=tk.LEFT, padx=3)
            e.bind("<Return>", lambda ev: self.refresh())
            e.bind("<FocusOut>", lambda ev: self.refresh())
            if v is self.vmin_var:
                ttk.Label(r, text="max:").pack(side=tk.LEFT)
        ttk.Button(r, text="Auto", width=5, command=self.auto_limits).pack(side=tk.LEFT, padx=3)
        self.range_lbl = ttk.Label(lf, text="blank = automatic limits", foreground="gray40")
        self.range_lbl.pack(anchor="w", padx=4)
        self.log_var = tk.BooleanVar(value=False)
        self.eq_var = tk.BooleanVar(value=True)
        self.title_var = tk.BooleanVar(value=True)
        self.zero_var = tk.BooleanVar(value=False)
        for text, v in (("Logarithmic colour scale", self.log_var),
                        ("Equal aspect ratio (true geometry)", self.eq_var),
                        ("Figure title", self.title_var),
                        ("Zero boundary nodes of n, T, flux (as convert.f90)", self.zero_var)):
            ttk.Checkbutton(lf, text=text, variable=v, command=self.refresh).pack(
                anchor="w", padx=4)
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=1)
        ttk.Label(r, text="Title text:").pack(side=tk.LEFT)
        self.ttext_var = tk.StringVar()
        e = ttk.Entry(r, textvariable=self.ttext_var, width=26)
        e.pack(side=tk.LEFT, padx=4)
        e.bind("<Return>", lambda ev: self.refresh())
        e.bind("<FocusOut>", lambda ev: self.refresh())
        ttk.Label(lf, text="(blank = automatic; mathtext such as $T_e$ allowed;\n untick 'Figure title' to remove it)",
                  foreground="gray40", justify="left").pack(anchor="w", padx=4)
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=1)
        ttk.Label(r, text="Font size (pt):").pack(side=tk.LEFT)
        self.fs_var = tk.StringVar(value="11")
        sp = ttk.Spinbox(r, from_=6, to=40, increment=1, width=5, textvariable=self.fs_var,
                         command=self.refresh)
        sp.pack(side=tk.LEFT, padx=4)
        sp.bind("<Return>", lambda ev: self.refresh())
        sp.bind("<FocusOut>", lambda ev: self.refresh())
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=1)
        ttk.Label(r, text="3D mode:").pack(side=tk.LEFT)
        self.mode_var = tk.StringVar(value="Cut-away (octant removed)")
        cb = ttk.Combobox(r, textvariable=self.mode_var, width=24, state="readonly",
                          values=["Cut-away (octant removed)", "Three orthogonal slices"])
        cb.pack(side=tk.LEFT, padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())

        # phi contours (2D tabs)
        lf = ttk.LabelFrame(f, text="Contour lines of the plotted quantity (2D tabs)")
        lf.pack(fill=tk.X, padx=6, pady=4)
        self.ct_master = tk.BooleanVar(value=False)
        ttk.Checkbutton(lf, text="Show contour lines", variable=self.ct_master,
                        command=self.refresh).pack(anchor="w", padx=4)
        self.ct_on, self.ct_lvl, self.ct_lab, self.ct_col, self.ct_btn = [], [], [], [], []
        for i, col in enumerate(("#ffffff", "#ff0000", "#000000", "#ff00ff")):
            r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=(3, 0))
            on = tk.BooleanVar(value=True)
            ttk.Checkbutton(r, text=f"{i+1}", variable=on, command=self.refresh).pack(side=tk.LEFT)
            lv = tk.StringVar()
            e = ttk.Entry(r, textvariable=lv, width=9)
            e.pack(side=tk.LEFT, padx=2)
            ttk.Label(r, textvariable=self.ct_unit, width=7).pack(side=tk.LEFT)
            b = tk.Button(r, width=3, bg=col, activebackground=col, relief="groove",
                          command=lambda i=i: self.pick_color(i))
            b.pack(side=tk.LEFT, padx=5)
            lb = tk.StringVar()
            e2 = ttk.Entry(r, textvariable=lb, width=14)
            e2.pack(side=tk.LEFT)
            for w in (e, e2):
                w.bind("<Return>", lambda ev: self.refresh())
                w.bind("<FocusOut>", lambda ev: self.refresh())
            self.ct_on.append(on); self.ct_lvl.append(lv); self.ct_lab.append(lb)
            self.ct_col.append(col); self.ct_btn.append(b)
        ttk.Label(lf, text="level (units shown) | colour | legend label (blank = auto)",
                  foreground="gray40").pack(anchor="w", padx=4)
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=2)
        ttk.Label(r, text="Legend:").pack(side=tk.LEFT)
        self.leg_var = tk.StringVar(value="best")
        cb = ttk.Combobox(r, textvariable=self.leg_var, width=12, state="readonly",
                          values=["best", "upper right", "upper left", "lower right",
                                  "lower left", "below plot"])
        cb.pack(side=tk.LEFT, padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self.refresh())
        ttk.Label(r, text="width:").pack(side=tk.LEFT)
        self.ctw_var = tk.StringVar(value="1.6")
        sp = ttk.Spinbox(r, from_=0.5, to=6, increment=0.5, width=4, textvariable=self.ctw_var,
                         command=self.refresh)
        sp.pack(side=tk.LEFT, padx=3)
        sp.bind("<Return>", lambda e: self.refresh())

        # flux streamlines
        lf = ttk.LabelFrame(f, text="Flux streamlines (2D tabs, flux quantities)")
        lf.pack(fill=tk.X, padx=6, pady=4)
        self.sl_on = tk.BooleanVar(value=True)
        self.sl_col = "#ffffff"
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=2)
        ttk.Checkbutton(r, text="Show", variable=self.sl_on, command=self.refresh).pack(side=tk.LEFT)
        self.sl_btn = tk.Button(r, width=3, bg=self.sl_col, activebackground=self.sl_col,
                                relief="groove", command=self.pick_sl_color)
        self.sl_btn.pack(side=tk.LEFT, padx=6)
        ttk.Label(r, text="density").pack(side=tk.LEFT)
        self.sl_dens = tk.StringVar(value="1.2")
        sp = ttk.Spinbox(r, from_=0.3, to=5, increment=0.1, width=4, textvariable=self.sl_dens,
                         command=self.refresh)
        sp.pack(side=tk.LEFT, padx=2)
        sp.bind("<Return>", lambda e: self.refresh())
        ttk.Label(r, text="width").pack(side=tk.LEFT)
        self.sl_w = tk.StringVar(value="0.8")
        sp = ttk.Spinbox(r, from_=0.3, to=4, increment=0.1, width=4, textvariable=self.sl_w,
                         command=self.refresh)
        sp.pack(side=tk.LEFT, padx=2)
        sp.bind("<Return>", lambda e: self.refresh())
        ttk.Label(lf, text="Streamlines of the in-plane projection of n<v> of the ptype\n"
                  "(shown for |n<v>| and n<vx>, n<vy>, n<vz>; option 2 only)",
                  foreground="gray40", justify="left").pack(anchor="w", padx=4)

        # saving
        lf = ttk.LabelFrame(f, text="Saving (right click on a figure)")
        lf.pack(fill=tk.X, padx=6, pady=4)
        r = ttk.Frame(lf); r.pack(fill=tk.X, padx=4, pady=2)
        ttk.Label(r, text="Width (in):").pack(side=tk.LEFT)
        self.w_var = tk.StringVar(value="6.5")
        ttk.Entry(r, textvariable=self.w_var, width=5).pack(side=tk.LEFT, padx=3)
        ttk.Label(r, text="PNG dpi:").pack(side=tk.LEFT)
        self.dpi_var = tk.StringVar(value="300")
        ttk.Entry(r, textvariable=self.dpi_var, width=5).pack(side=tk.LEFT, padx=3)
        ttk.Button(lf, text="Reset 3D view", command=lambda: self.reset_view(self.tabs["3D"])
                   ).pack(pady=3)

    def _build_tabs(self):
        self.nb = ttk.Notebook(self)
        self.tabs = {}
        for kind, text in (("XY", "XY plane"), ("XZ", "XZ plane"), ("YZ", "YZ plane"),
                           ("3D", "3D view")):
            t = Tab(self, self.nb, kind)
            self.tabs[kind] = t
            self.nb.add(t.frame, text=text)
        self.nb.bind("<<NotebookTabChanged>>", lambda e: self.draw_current())
        self.tabs["3D"].canvas.mpl_connect("scroll_event", self.on_scroll3d)

    # ---- file handling ----------------------------------------------------------------
    def browse_data(self):
        p = filedialog.askopenfilename(title="phi_n_3D.dat or 3D_param file",
                                       filetypes=[("3D maps", "phi_n_3D*.dat 3D_param*.dat"), ("All", "*")])
        if p:
            self.f_var.set(p)
            self.load()

    def browse_cond(self):
        p = filedialog.askopenfilename(title="conditions.inp",
                                       filetypes=[("conditions.inp", "*.inp"), ("All", "*")])
        if p:
            self.c_var.set(p)
            if self.data:
                self.fill_info()

    def load(self):
        path = self.f_var.get().strip()
        if not os.path.isfile(path):
            messagebox.showerror("File", f"File not found:\n{path}")
            return
        self.config(cursor="watch"); self.update_idletasks()
        try:
            cond = self.c_var.get().strip()
            d = Param3D(path, os.path.join(os.path.dirname(cond), "geometry.inp") if cond else None)
        except Exception as e:
            self.config(cursor="")
            messagebox.showerror("Read error", str(e))
            return
        self.config(cursor="")
        self.data = d
        self.pos = [0.5 * L for L in d.L]
        self.zoom3d = 1.0
        self.view3d = (25.0, 45.0, 0.0)
        self.set_sliders()
        self.rebuild_vars(first=True)
        self.fill_info()
        self.status.set(f"Loaded {os.path.basename(path)}")
        self.refresh()

    def fill_info(self):
        d = self.data
        lines = [f"Grid   : {d.nx} x {d.ny} x {d.nz} cells ({d.nx+1} x {d.ny+1} x {d.nz+1} nodes)",
                 f"dx,dy,dz [mm]: {d.dx*1e3:.4g}, {d.dy*1e3:.4g}, {d.dz*1e3:.4g}",
                 f"Box    : {d.L[0]*100:.4g} x {d.L[1]*100:.4g} x {d.L[2]*100:.4g} cm",
                 f"ptype  : {d.ntype}   parameters/node: {d.navg}"]
        for i, m in enumerate(d.mass):
            lines.append(f"  mass({i+1}) = {m:.4e} kg = {m/AMU:.4g} amu")
        if d.lepic:
            lines.append(f"LePIC+ phi_n_3D: phi and n of ptype {d.dens_ptype} "
                         f"({d.names[d.dens_ptype - 1] if d.dens_ptype <= len(d.names) else '?'})")
            if d.geom_note:
                lines.append(d.geom_note)
            ng = parse_ngrid(self.c_var.get()) if self.c_var.get() else None
            if ng is not None and ng < 0:
                lines.append("conditions.inp: wall labels < 0 -> 3D averages OFF (file may be stale)")
        opt = ("1: phi, n" if d.option == 1 else
               "2: phi, n, <v>, <v^2> per ptype" + (" (+B, not plotted)" if d.hasB else ""))
        lines.append(f"File content : option {opt}")
        c = parse_conditions(self.c_var.get()) if (self.c_var.get() and not d.lepic) else None
        if c:
            txt = f"conditions.inp: Save 3D = {'yes' if c['save'] else 'no'}, opt = {c['opt']}"
            if c["opt"] is not None and c["opt"] != d.option:
                txt += "  (differs from file; file header used)"
            lines.append(txt)
        self.info.configure(state="normal")
        self.info.delete("1.0", "end")
        self.info.insert("end", "\n".join(lines))
        self.info.configure(state="disabled")

    # ---- variables ------------------------------------------------------------------------
    def species(self, p):
        names = [s.strip() for s in self.sp_var.get().split(",")] if self.sp_var.get().strip() else []
        if p - 1 < len(names) and names[p - 1]:
            return r"\mathrm{%s}" % names[p - 1]
        if self.data and self.data.lepic:
            nm = self.data.names[p - 1] if p <= len(self.data.names) else str(p)
            return r"\mathrm{%s}" % (nm.strip("[]").replace("+", "^+").replace("-", "^-") or p)
        if self.data and self.data.option == 1:
            return r"\mathrm{e}"
        return str(p)

    def rebuild_vars(self, first=False):
        if not self.data:
            return
        d = self.data
        keep = self.cur_var()["key"] if (self.vars and not first) else None
        V = [dict(key="phi", kind="phi", p=0, c=0, sym=r"\phi", unit="V", uplain="V", menu="phi   (potential)")]
        for c in range(3):
            V.append(dict(key=f"E{c}", kind="E", p=0, c=c, sym=f"E_{AXN[c]}", unit=r"V m$^{-1}$",
                          uplain="V/m", menu=f"E{AXN[c]}   (electric field, from phi)"))
        plist = [1] if d.option == 1 else list(range(1, d.ntype + 1))
        for p in plist:
            s = self.species(d.dens_ptype if d.lepic else p)
            tag = f"ptype {d.dens_ptype if d.lepic else p}"
            V.append(dict(key=f"n{p}", kind="n", p=p, c=0, sym=f"n_{{{s}}}", unit=r"m$^{-3}$", uplain="m⁻³",
                          menu=f"n   [{tag}]"))
            if d.option == 2:
                V.append(dict(key=f"T{p}", kind="T", p=p, c=0, sym=f"T_{{{s}}}", unit="eV", uplain="eV",
                              menu=f"T   [{tag}]"))
                V.append(dict(key=f"fn{p}", kind="fnorm", p=p, c=0,
                              sym=rf"|n_{{{s}}}\langle\mathbf{{v}}\rangle_{{{s}}}|",
                              unit=r"m$^{-2}$ s$^{-1}$", uplain="m⁻² s⁻¹",
                              menu=f"|n<v>|   [{tag}]  (flux norm)"))
                for c in range(3):
                    V.append(dict(key=f"f{p}{c}", kind="flux", p=p, c=c,
                                  sym=rf"n_{{{s}}}\langle v_{AXN[c]}\rangle_{{{s}}}",
                                  unit=r"m$^{-2}$ s$^{-1}$", uplain="m⁻² s⁻¹", menu=f"n<v{AXN[c]}>   [{tag}]"))
        self.vars = V
        self.var_cb["values"] = [v["menu"] for v in V]
        idx = 0
        if keep:
            for i, v in enumerate(V):
                if v["key"] == keep:
                    idx = i
        self.var_cb.current(idx)
        if not first:
            self.refresh()

    def cur_var(self):
        i = max(0, self.var_cb.current())
        return self.vars[i]

    @staticmethod
    def vlabel(v):
        return rf"${v['sym']}$ ({v['unit']})"

    # ---- planes / sliders ------------------------------------------------------------------
    @property
    def uf(self):
        return UNITS[self.unit_var.get()]

    def set_sliders(self):
        d = self.data
        self._lock = True
        for a in range(3):
            self.sl[a].configure(from_=0.0, to=d.L[a] * self.uf)
            self.sl[a].set(self.pos[a] * self.uf)
            self.sl_var[a].set(f"{self.pos[a]*self.uf:.4g}")
        self._lock = False

    def on_unit(self):
        if self.data:
            self.set_sliders()
            self.refresh()

    def on_slider(self, a, val):
        if self._lock or not self.data:
            return
        self.pos[a] = val / self.uf
        self.sl_var[a].set(f"{val:.4g}")
        self.schedule()

    def on_entry(self, a):
        if not self.data:
            return
        try:
            v = float(self.sl_var[a].get()) / self.uf
        except ValueError:
            return
        v = min(max(v, 0.0), self.data.L[a])
        if abs(v - self.pos[a]) < 1e-15:
            return
        self.pos[a] = v
        self._lock = True
        self.sl[a].set(v * self.uf)
        self._lock = False
        self.schedule()

    def schedule(self):
        if self._job:
            self.after_cancel(self._job)
        self._job = self.after(80, self.refresh)

    @staticmethod
    def make_ids(n, s):
        return np.unique(np.append(np.arange(0, n + 1, s), n)).astype(int)

    def stride(self):
        return RES[self.res_var.get()]

    def stride3d(self):
        s = self.stride()
        if self.cap_var.get():
            try:
                cap = max(8, int(float(self.cap_n.get())))
            except ValueError:
                cap = 96
            s = max(s, int(math.ceil(max(self.data.n) / cap)))
        return s

    def node_index(self, a, ids):
        d = self.data
        return int(ids[a][np.argmin(np.abs(ids[a] * d.h[a] - self.pos[a]))])

    # ---- refresh / drawing ---------------------------------------------------------------------
    def sync_contours(self):
        """Contour levels are in the units of the plotted quantity: keep them per quantity."""
        v = self.cur_var()
        if self._ct_key != v["key"]:
            if self._ct_key is not None:
                self.ct_store[self._ct_key] = [(self.ct_lvl[i].get(), self.ct_lab[i].get())
                                               for i in range(4)]
            for i, (lv, lb) in enumerate(self.ct_store.get(v["key"], [("", "")] * 4)):
                self.ct_lvl[i].set(lv)
                self.ct_lab[i].set(lb)
            self._ct_key = v["key"]
        self.ct_unit.set(v["uplain"])

    def refresh(self):
        self._job = None
        if not self.data:
            return
        self.sync_contours()
        for t in self.tabs.values():
            t.dirty = True
        d = self.data
        ids = [self.make_ids(d.n[a], self.stride()) for a in range(3)]
        for a in range(3):
            i = self.node_index(a, ids)
            self.node_lbl[a].configure(
                text=f"→ node {AXN[a]}index = {i+1}/{d.n[a]+1},  {AXN[a]} = {i*d.h[a]*self.uf:.5g} {self.unit_var.get()}")
        self.draw_current()

    def draw_current(self):
        if not self.data:
            return
        kind = ("XY", "XZ", "YZ", "3D")[self.nb.index(self.nb.select())]
        t = self.tabs[kind]
        if t.dirty:
            self.config(cursor="watch"); self.update_idletasks()
            try:
                if kind == "3D":
                    self.render3d(t)
                else:
                    self.render2d(t, kind)
                t.dirty = False
            except Exception as e:
                traceback.print_exc()
                self.status.set(f"Plot error: {e}")
            finally:
                self.config(cursor="")
            t.canvas.draw()
            if kind == "3D" and self.place_cbar3d(t):
                t.canvas.draw()

    def place_cbar3d(self, t):
        """Put the 3D colour bar right next to the projected bounding box of the plot."""
        ax, cax = t.ax, getattr(t, "cax", None)
        if ax is None or cax is None:
            return False
        xl, yl, zl = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
        C = np.array([[x, y, z] for x in xl for y in yl for z in zl])
        xs, ys, _ = proj3d.proj_transform(C[:, 0], C[:, 1], C[:, 2], ax.get_proj())
        disp = ax.transData.transform(np.column_stack([xs, ys]))
        W, H = t.fig.bbox.width, t.fig.bbox.height
        x1 = disp[:, 0].max() / W
        y0, y1 = disp[:, 1].min() / H, disp[:, 1].max() / H
        h = 0.62 * (y1 - y0)
        cax.set_position([min(x1 + 0.045, 0.94), 0.5 * (y0 + y1) - 0.5 * h, 0.022, h])
        return True

    def reset_view(self, t):
        if t.kind == "3D":
            self.zoom3d = 1.0
            self.view3d = (25.0, 45.0, 0.0)
            t.dirty = True
            self.draw_current()
        else:
            t.toolbar.home()

    # ---- phi contours ----------------------------------------------------------------------
    def pick_color(self, i):
        c = colorchooser.askcolor(color=self.ct_col[i], title=f"Contour {i+1} colour")[1]
        if c:
            self.ct_col[i] = c
            self.ct_btn[i].configure(bg=c, activebackground=c)
            self.refresh()

    @staticmethod
    def _fmt_level(v):
        """mathtext number: 1.8e17 -> 1.8\\times10^{17}"""
        if v != 0 and (abs(v) >= 1e4 or abs(v) < 1e-3):
            m, e = f"{v:.3e}".split("e")
            m = m.rstrip("0").rstrip(".")
            return rf"{m}\times10^{{{int(e)}}}"
        return f"{v:g}"

    def pick_sl_color(self):
        c = colorchooser.askcolor(color=self.sl_col, title="Streamline colour")[1]
        if c:
            self.sl_col = c
            self.sl_btn.configure(bg=c, activebackground=c)
            self.refresh()

    @staticmethod
    def _uniform(c, U, axis):
        """Linear resampling of U along `axis` onto a uniform grid (streamplot needs one)."""
        cu = np.linspace(c[0], c[-1], len(c))
        if np.allclose(c, cu, rtol=1e-9, atol=0.0):
            return c, U
        out = np.empty_like(U)
        if axis == 0:
            for j in range(U.shape[1]):
                out[:, j] = np.interp(cu, c, U[:, j])
        else:
            for i in range(U.shape[0]):
                out[i, :] = np.interp(cu, c, U[i, :])
        return cu, out

    def draw_streamlines(self, ax, var, axis, idx, ids, a1, a2, c1, c2):
        if not self.sl_on.get() or var["kind"] not in ("flux", "fnorm"):
            return ""
        try:
            dens, lw = float(self.sl_dens.get()), float(self.sl_w.get())
        except ValueError:
            dens, lw = 1.2, 0.8
        G = self.data.flux_vec(var["p"], axis, idx, ids, self.zero_var.get())
        U, V = G[a1], G[a2]
        if not (np.any(U) or np.any(V)):
            return "   – no in-plane flux component in this plane: no streamlines"
        x, U = self._uniform(c1, U, 0)
        x, V = self._uniform(c1, V, 0)
        y, U = self._uniform(c2, U, 1)
        y, V = self._uniform(c2, V, 1)
        ax.streamplot(x, y, U.T, V.T, color=self.sl_col, linewidth=lw, density=dens,
                      arrowsize=1.0, arrowstyle="-|>", zorder=3)
        return ""

    def draw_contours(self, fig, ax, D, var, c1, c2):
        if not self.ct_master.get():
            return []
        try:
            lw = float(self.ctw_var.get())
        except ValueError:
            lw = 1.6
        P, handles, skipped = D, [], []
        for i in range(4):
            if not self.ct_on[i].get():
                continue
            txt = self.ct_lvl[i].get().strip()
            if not txt:
                continue
            try:
                lv = float(txt.replace("d", "e").replace("D", "e"))
            except ValueError:
                skipped.append(f"{i+1} (invalid level)")
                continue
            if not (np.nanmin(P) <= lv <= np.nanmax(P)):
                skipped.append(f"{i+1} (level outside the range {np.nanmin(P):.4g}…"
                               f"{np.nanmax(P):.4g} {var['uplain']})")
                continue
            col = self.ct_col[i]
            ax.contour(c1, c2, P.T, levels=[lv], colors=[col], linewidths=lw, linestyles="solid")
            lab = self.ct_lab[i].get().strip() or rf"${var['sym']}$ = ${self._fmt_level(lv)}$ {var['unit']}"
            handles.append(Line2D([], [], color=col, lw=lw, label=lab))
        if handles:
            loc = self.leg_var.get()
            kw = dict(handles=handles, framealpha=0.9, edgecolor="0.3", fancybox=False,
                      fontsize=mpl.rcParams["font.size"] - 1)
            if loc == "below plot":
                fig.subplots_adjust(bottom=0.26)
                ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
                          ncol=min(len(handles), 2), **kw)
            else:
                ax.legend(loc=loc, **kw)
        return skipped

    # ---- fonts / title -------------------------------------------------------------------
    def apply_fonts(self):
        try:
            fs = float(self.fs_var.get())
        except ValueError:
            fs = 11.0
        fs = min(max(fs, 5.0), 48.0)
        mpl.rcParams.update({"font.size": fs, "axes.labelsize": fs + 1, "axes.titlesize": fs + 1,
                             "xtick.labelsize": fs - 0.5, "ytick.labelsize": fs - 0.5})
        return fs

    def title_text(self, auto):
        """None if the title is switched off, else custom text or the automatic one."""
        if not self.title_var.get():
            return None
        return self.ttext_var.get().strip() or auto

    # ---- colour handling ----------------------------------------------------------------------
    def _lim(self, v):
        try:
            return float(v.get().replace("d", "e").replace("D", "e"))
        except ValueError:
            return None

    def auto_limits(self):
        self.vmin_var.set("")
        self.vmax_var.set("")
        self.refresh()

    def make_norm(self, arrays):
        vals = np.concatenate([a[np.isfinite(a)].ravel() for a in arrays]) if arrays else np.zeros(1)
        if vals.size == 0:
            vals = np.array([0.0, 1.0])
        vmin, vmax = self._lim(self.vmin_var), self._lim(self.vmax_var)
        self.range_lbl.configure(text=f"data range in view: {vals.min():.4g} … {vals.max():.4g}"
                                      "   (blank = automatic)")
        for nm, v in (("min", self.vmin_var), ("max", self.vmax_var)):
            if v.get().strip() and self._lim(v) is None:
                self.status.set(f"Colour scale {nm} '{v.get()}' is not a number: ignored.")
        if self.log_var.get():
            pos = vals[vals > 0]
            if pos.size:
                lo = vmin if (vmin is not None and vmin > 0) else pos.min()
                hi = vmax if (vmax is not None and vmax > 0) else pos.max()
                if hi <= lo:
                    hi = lo * 10.0
                return LogNorm(lo, hi)
            self.status.set("No positive values: logarithmic scale not applicable.")
        lo = vmin if vmin is not None else float(vals.min())
        hi = vmax if vmax is not None else float(vals.max())
        if hi <= lo:
            hi = lo + (abs(lo) * 1e-3 if lo != 0 else 1.0)
        return Normalize(lo, hi)

    def cmap(self):
        name = self.cmap_var.get()
        if name in DISCRETE:
            c = mpl.colormaps[DISCRETE[name]].resampled(10)      # 10 discrete colours
        else:
            c = mpl.colormaps[name].copy()
        c.set_bad("0.85")
        return c

    def style_colorbar(self, cb, norm, label):
        cb.set_label(label)
        cb.outline.set_linewidth(0.9)
        if self.cmap_var.get() in DISCRETE and not isinstance(norm, LogNorm):
            cb.set_ticks(norm.inverse(np.linspace(0.0, 1.0, 11)))   # ticks at the colour transitions
            cb.minorticks_off()
        if not isinstance(norm, LogNorm):
            f = ScalarFormatter(useMathText=True)
            f.set_powerlimits((-3, 4))
            cb.formatter = f
            cb.update_ticks()
            cb.ax.yaxis.get_offset_text().set_fontsize(mpl.rcParams["ytick.labelsize"])
        cb.ax.tick_params(direction="in", which="both", right=True)

    # ---- 2D ------------------------------------------------------------------------------------
    def render2d(self, t, kind):
        d, var = self.data, self.cur_var()
        axis = {"XY": 2, "XZ": 1, "YZ": 0}[kind]
        a1, a2 = [a for a in range(3) if a != axis]
        s = self.stride()
        ids = [self.make_ids(d.n[a], s) for a in range(3)]
        idx = self.node_index(axis, ids)
        D = d.plane(var, axis, idx, ids, self.zero_var.get())
        uf, un = self.uf, self.unit_var.get()
        c1, c2 = ids[a1] * d.h[a1] * uf, ids[a2] * d.h[a2] * uf
        norm = self.make_norm([D])
        self.apply_fonts()
        fig = t.fig
        fig.clear()
        fig.subplots_adjust(left=0.13, right=0.80, bottom=0.13, top=0.90)
        ax = fig.add_subplot(111)
        t.ax = ax
        shading = "gouraud" if self.shade_var.get().startswith("Smooth") else "nearest"
        Z = np.ma.masked_invalid(D.T)
        m = ax.pcolormesh(c1, c2, Z, norm=norm, cmap=self.cmap(), shading=shading,
                          rasterized=True)
        ax.set_xlim(c1[0], c1[-1])
        ax.set_ylim(c2[0], c2[-1])
        sl_msg = self.draw_streamlines(ax, var, axis, idx, ids, a1, a2, c1, c2)
        skipped = self.draw_contours(fig, ax, D, var, c1, c2)
        ax.set_xlim(c1[0], c1[-1])
        ax.set_ylim(c2[0], c2[-1])
        ax.set_aspect("equal" if self.eq_var.get() else "auto")
        ax.set_xlabel(f"${AXN[a1]}$ ({un})")
        ax.set_ylabel(f"${AXN[a2]}$ ({un})")
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which="both", direction="in", top=True, right=True)
        cax = make_axes_locatable(ax).append_axes("right", size="4.5%", pad=0.12)
        cb = fig.colorbar(m, cax=cax)
        self.style_colorbar(cb, norm, self.vlabel(var))
        ttl = self.title_text(rf"${var['sym']}$,  {kind} plane,  ${AXN[axis]}$ = "
                              rf"{idx*d.h[axis]*uf:.4g} {un}")
        if ttl is not None:
            ax.set_title(ttl)
        if t.toolbar:
            t.toolbar.update()
        t.new_artists()

        def hover(ev):
            if ev.inaxes is not ax or ev.xdata is None:
                return None
            i = int(np.argmin(np.abs(c1 - ev.xdata)))
            j = int(np.argmin(np.abs(c2 - ev.ydata)))
            val = D[i, j]
            idxs = [0, 0, 0]
            idxs[axis], idxs[a1], idxs[a2] = idx, ids[a1][i], ids[a2][j]
            px, py = ax.transData.transform((c1[i], c2[j]))
            txt = (f"${AXN[a1]}$ = {c1[i]:.5g} {un}\n${AXN[a2]}$ = {c2[j]:.5g} {un}\n"
                   f"${AXN[axis]}$ = {idx*d.h[axis]*uf:.5g} {un}\n"
                   f"${var['sym']}$ = {val:.6g} {var['unit']}\n"
                   f"node (ix,iy,iz) = ({idxs[0]+1}, {idxs[1]+1}, {idxs[2]+1})")
            return px, py, txt
        t.hover = hover
        self.status.set(f"{kind}: {var['menu'].split('  ')[0]}, stride {s} "
                        f"({len(c1)} x {len(c2)} nodes shown)" + sl_msg
                        + (("   – contour(s) skipped: " + "; ".join(skipped)) if skipped else ""))

    # ---- 3D ------------------------------------------------------------------------------------
    def render3d(self, t):
        d, var = self.data, self.cur_var()
        s = self.stride3d()
        ids = [self.make_ids(d.n[a], s) for a in range(3)]
        uf, un = self.uf, self.unit_var.get()
        zero = self.zero_var.get()
        cut = self.mode_var.get().startswith("Cut")
        sidx = [self.node_index(a, ids) for a in range(3)]
        spos = [sidx[a] * d.h[a] * uf for a in range(3)]       # slice coordinates (display unit)
        coords = [ids[a] * d.h[a] * uf for a in range(3)]
        faces = []   # dict(axis, idx, D, mask(cells), nodevis, name)
        for a in range(3):
            a1, a2 = [b for b in range(3) if b != a]
            c1c = 0.5 * (coords[a1][:-1] + coords[a1][1:])
            c2c = 0.5 * (coords[a2][:-1] + coords[a2][1:])
            beyond = (c1c[:, None] > spos[a1]) & (c2c[None, :] > spos[a2])
            full = np.ones((len(c1c), len(c2c)), bool)
            if not cut:
                specs = [(sidx[a], full, "slice")]
            else:
                specs = [(0, full, "lower face"),
                         (d.n[a], ~beyond, "upper face"),
                         (sidx[a], beyond, "cut face")]
            for ix, mask, nm in specs:
                D = d.plane(var, a, ix, ids, zero)
                nv = np.zeros(D.shape, bool)
                nv[:-1, :-1] |= mask; nv[1:, :-1] |= mask
                nv[:-1, 1:] |= mask; nv[1:, 1:] |= mask
                faces.append(dict(a=a, a1=a1, a2=a2, idx=ix, D=D, mask=mask, nv=nv, name=nm))
        norm = self.make_norm([f["D"][f["nv"]] for f in faces])
        cmap = self.cmap()

        # remember view
        old = t.ax
        if old is not None and hasattr(old, "elev"):
            self.view3d = (old.elev, old.azim, getattr(old, "roll", 0.0))
        self.apply_fonts()
        fig = t.fig
        fig.clear()
        fig.subplots_adjust(left=0.06, right=0.90, bottom=0.02, top=0.95)
        ax = fig.add_subplot(111, projection="3d")
        t.ax = ax
        verts, cols, pick_P, pick_V, pick_name = [], [], [], [], []
        for f in faces:
            a, a1, a2 = f["a"], f["a1"], f["a2"]
            c1, c2 = coords[a1], coords[a2]
            P = np.empty((len(c1), len(c2), 3))
            P[..., a] = f["idx"] * d.h[a] * uf
            P[..., a1] = c1[:, None]
            P[..., a2] = c2[None, :]
            D = f["D"]
            q = np.stack([P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]], axis=2)
            cv = 0.25 * (D[:-1, :-1] + D[1:, :-1] + D[1:, 1:] + D[:-1, 1:])
            m = f["mask"]
            verts.append(q[m].reshape(-1, 4, 3))
            rgba = cmap(norm(np.ma.masked_invalid(cv[m])))
            cols.append(rgba)
            nvm = f["nv"]
            pick_P.append(P[nvm])
            pick_V.append(D[nvm])
            pick_name.append(np.full(nvm.sum(), f["name"] + f" ({AXN[a]}-normal)", dtype=object))
        verts = np.concatenate(verts)
        cols = np.concatenate(cols)
        coll = Poly3DCollection(verts, facecolors=cols, edgecolors=cols, linewidths=0.25)
        coll.set_rasterized(True)
        ax.add_collection3d(coll)
        L = [coords[a][-1] for a in range(3)]
        if not cut:   # bounding box for context
            for a in range(3):
                for b in (0, L[(a + 1) % 3]):
                    for c in (0, L[(a + 2) % 3]):
                        p0 = [0, 0, 0]; p1 = [0, 0, 0]
                        p0[(a + 1) % 3] = p1[(a + 1) % 3] = b
                        p0[(a + 2) % 3] = p1[(a + 2) % 3] = c
                        p1[a] = L[a]
                        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]],
                                color="0.25", lw=0.7)
        ax.set_xlim(0, L[0]); ax.set_ylim(0, L[1]); ax.set_zlim(0, L[2])
        self.apply_box(ax)
        ax.view_init(elev=self.view3d[0], azim=self.view3d[1], roll=self.view3d[2])
        ax.mouse_init(rotate_btn=1, pan_btn=2, zoom_btn=[])   # right button -> save menu
        ax.set_xlabel(f"$x$ ({un})", labelpad=8)
        ax.set_ylabel(f"$y$ ({un})", labelpad=8)
        ax.set_zlabel(f"$z$ ({un})", labelpad=6)
        ax.tick_params(which="major", pad=1)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor((1, 1, 1, 1.0))
            axis.pane.set_edgecolor("0.55")
        ax.grid(False)
        cax = fig.add_axes([0.80, 0.30, 0.022, 0.40])    # re-positioned next to the box after drawing
        t.cax = cax
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        self.style_colorbar(cb, norm, self.vlabel(var))
        ttl = self.title_text(rf"${var['sym']}$,  " + ("cut-away" if cut else "slices")
                              + rf" at ($x,y,z$) = ({spos[0]:.3g}, {spos[1]:.3g}, {spos[2]:.3g}) {un}")
        if ttl is not None:
            ax.set_title(ttl, pad=0)
        t.new_artists()

        PP = np.concatenate(pick_P)
        VV = np.concatenate(pick_V)
        NN = np.concatenate(pick_name)

        def hover(ev):
            if ev.inaxes is not ax or ev.x is None:
                return None
            xs, ys, zs = proj3d.proj_transform(PP[:, 0], PP[:, 1], PP[:, 2], ax.get_proj())
            disp = ax.transData.transform(np.column_stack([xs, ys]))
            d2 = (disp[:, 0] - ev.x) ** 2 + (disp[:, 1] - ev.y) ** 2
            cand = np.nonzero(d2 < 14.0 ** 2)[0]
            if cand.size == 0:
                return None
            tol = 0.04 * (zs.max() - zs.min())
            front = cand[zs[cand] <= zs[cand].min() + tol]
            k = front[np.argmin(d2[front])]
            p = PP[k]
            ijk = [int(round(p[a] / (d.h[a] * uf))) + 1 for a in range(3)]
            txt = (f"$x$ = {p[0]:.5g} {un}\n$y$ = {p[1]:.5g} {un}\n$z$ = {p[2]:.5g} {un}\n"
                   f"${var['sym']}$ = {VV[k]:.6g} {var['unit']}\n"
                   f"node (ix,iy,iz) = ({ijk[0]}, {ijk[1]}, {ijk[2]})\n{NN[k]}")
            return disp[k, 0], disp[k, 1], txt
        t.hover = hover
        msg = f"3D: stride {s} ({len(coords[0])} x {len(coords[1])} x {len(coords[2])} nodes)"
        if s != self.stride():
            msg += f" – decimated from stride {self.stride()} for fast rotation (see 'Grid resolution')"
        self.status.set(msg)

    def apply_box(self, ax):
        d = self.data
        asp = tuple(d.L) if self.eq_var.get() else (1.0, 1.0, 1.0)
        m = max(asp)
        ax.set_box_aspect(tuple(v / m for v in asp), zoom=self.zoom3d)

    def on_scroll3d(self, ev):
        t = self.tabs["3D"]
        if t.ax is None or not self.data:
            return
        step = ev.step if ev.step else (1 if ev.button == "up" else -1)
        self.zoom3d = float(np.clip(self.zoom3d * 1.1 ** step, 0.3, 6.0))
        self.apply_box(t.ax)
        t._hide()
        t.canvas.draw_idle()

    # ---- saving --------------------------------------------------------------------------------
    def save_figure(self, t, fmt):
        if not self.data:
            return
        var = self.cur_var()
        name = f"{var['key']}_{t.kind}"
        path = filedialog.asksaveasfilename(
            title=f"Save figure as {fmt.upper()}", defaultextension="." + fmt,
            initialfile=f"{name}.{fmt}", filetypes=[(fmt.upper(), "*." + fmt)])
        if not path:
            return
        try:
            dpi = int(float(self.dpi_var.get()))
        except ValueError:
            dpi = 300
        fig = t.fig
        old = fig.get_size_inches().copy()
        t._hide()
        self.config(cursor="watch"); self.update_idletasks()
        try:
            try:
                w = float(self.w_var.get())
                if w > 0:
                    fig.set_size_inches(w, w * old[1] / old[0], forward=False)
            except ValueError:
                pass
            kw = dict(facecolor="white")
            if t.kind != "3D":      # tight cropping of 3D axes can cut off the z label
                kw.update(bbox_inches="tight", pad_inches=0.06)
            if fmt == "png":
                kw["dpi"] = dpi
            else:
                kw["dpi"] = 300      # resolution of the rasterised field (text/axes stay vector)
            fig.savefig(path, format=fmt, **kw)
            self.status.set(f"Saved {path}")
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save error", str(e))
        finally:
            fig.set_size_inches(old, forward=False)
            self.config(cursor="")
            t.canvas.draw_idle()


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    app = App(args[0] if args else None, args[1] if len(args) > 1 else None)
    app.mainloop()


if __name__ == "__main__":
    main()
