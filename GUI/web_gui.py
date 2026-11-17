#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
web_gui.py  --  browser front end of the LePIC+ GUIs (no tkinter needed)
========================================================================

    python3 GUI/web_gui.py [--port 8765] [--root <LePIC+ folder>]

Starts a small web server on 127.0.0.1 (Python standard library + numpy + matplotlib).
In VS Code Remote-SSH the port is forwarded automatically: open the link VS Code shows
(or the "Ports" panel).  From a plain ssh session:  ssh -L 8765:localhost:8765 <host>,
then open http://localhost:8765 on your own computer.

Three pages:
* Inputs     - the input editor (conditions.inp, geometry.inp, boundary.inp, coil.inp).
               Same parser, checks and save logic as pic_input_gui.py (pic_inputs.py).
* Output 2D  - plane maps and time series of Output/ (readers of plot_2D_gui.py).
* Output 3D  - slices of Output/phi_n_3D.dat.

Only folders inside --root (default: the folder that contains GUI/) can be opened.
The server is single-user: one case is edited at a time.
"""
import argparse
import base64
import io
import json
import math
import os
import re
import struct
import sys
import threading
import traceback
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                      # noqa: E402
import matplotlib                       # noqa: E402
matplotlib.use('Agg')
from matplotlib.figure import Figure    # noqa: E402

import pic_inputs as M                  # noqa: E402
import plot_2D_gui as V                 # noqa: E402
from gui_help import HELP, help_key     # noqa: E402

ROOT = os.path.dirname(HERE)
LOCK = threading.Lock()
STATE = {'case': None}
PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
           '#bcbd22', '#17becf']


def inside_root(path):
    p = os.path.realpath(os.path.join(ROOT, path))
    r = os.path.realpath(ROOT)
    if p != r and not p.startswith(r + os.sep):
        raise ValueError('Only folders inside %s can be opened.' % ROOT)
    return p


def rel(path):
    return os.path.relpath(path, ROOT)


# =====================================================================================
#  Inputs: token access by GUI key
# =====================================================================================
_BK = {nm: i for i, nm in enumerate(M.BKEYS)}


def seg_kind(r):
    t = r.get(6)
    return '0' if t.lstrip('+-') in ('0', '') else ('-' if M.is_neg(t) else '+')


def tok_get(case, key):
    c = case.cond
    if key in M.FIELDS:
        return c.tok(key)
    m = re.match(r'^B(\d+)\.(\w+)$', key)
    if m:
        return c.bfiles[int(m.group(1))].get(_BK[m.group(2)])
    m = re.match(r'^rg(\d+)\.(\d)$', key)
    if m:
        return c.rec['pdf'].get(7 + 2 * int(m.group(1)) + int(m.group(2)))
    m = re.match(r'^b\.V(\d+)$', key)
    if m:
        return case.bnd.vals[int(m.group(1)) - 1].get(0)
    if key == 'coil.n':
        return case.coil.count.get(0)
    m = re.match(r'^coil(\d+)\.(\d)$', key)
    if m:
        return case.coil.turns[int(m.group(1))].get(int(m.group(2)))
    m = re.match(r'^g\.box\.(\d)$', key)
    if m:
        return case.geom.box.get(int(m.group(1)))
    m = re.match(r'^g\.seg(\d+)\.(\d)$', key)
    if m:
        return case.geom.segs[int(m.group(1))].get(int(m.group(2)))
    raise KeyError(key)


def tok_set(case, key, v):
    """Apply one edit. Special keys mirror the Tk editor's compound widgets."""
    c = case.cond
    v = (v or '').strip()
    if key == 'nB.count':
        n = int(v)
        neg = M.is_neg(c.tok('nB'))
        c.set_nb(n)
        c.set_tok('nB', ('-' if neg else '') + str(max(0, min(M.MAX_B, n))))
    elif key == 'nB.neg':
        c.set_tok('nB', M.with_sign(M.strip_sign(c.tok('nB')), v == '1'))
    elif key == 'opt_inj.mode':
        c.set_tok('opt_inj', M.with_sign(v, M.is_neg(c.tok('opt_inj'))))
    elif key == 'opt_inj.cos':
        c.set_tok('opt_inj', M.with_sign(M.strip_sign(c.tok('opt_inj')), v == '1'))
    elif key == 'coil.n':
        case.coil.set_count(int(v))
    elif key == 'ngrid':
        c.set_tok('ngrid', v)
        case.sync_boundary()
    elif key in M.FIELDS:
        c.set_tok(key, v)
    elif re.match(r'^g\.seg\d+\.', key):
        m = re.match(r'^g\.seg(\d+)\.(kind|label|flag\d|\d)$', key)
        r = case.geom.segs[int(m.group(1))]
        what = m.group(2)
        if what == 'kind':
            mag = M.strip_sign(r.get(6)) or '1'
            if mag == '0':
                mag = '1'
            r.set(6, '0' if v == '0' else ('-' + mag if v == '-' else mag))
        elif what == 'label':
            if seg_kind(r) != '0':
                r.set(6, M.with_sign(v, seg_kind(r) == '-'))
        elif what.startswith('flag'):
            i = int(what[4:])
            r.set(i, M.with_sign(r.get(i), v == '1'))
        else:
            i = int(what)
            r.set(i, M.with_sign(v, M.is_neg(r.get(i))))
    else:
        m = re.match(r'^B(\d+)\.(\w+)$', key)
        if m:
            c.bfiles[int(m.group(1))].set(_BK[m.group(2)], v)
            return
        m = re.match(r'^rg(\d+)\.(\d)$', key)
        if m:
            c.rec['pdf'].set(7 + 2 * int(m.group(1)) + int(m.group(2)), v)
            return
        m = re.match(r'^b\.V(\d+)$', key)
        if m:
            case.bnd.vals[int(m.group(1)) - 1].set(0, v)
            return
        m = re.match(r'^coil(\d+)\.(\d)$', key)
        if m:
            case.coil.turns[int(m.group(1))].set(int(m.group(2)), v)
            return
        m = re.match(r'^g\.box\.(\d)$', key)
        if m:
            case.geom.box.set(int(m.group(1)), v)
            return
        raise KeyError(key)


def do_action(case, act, arg):
    c, g = case.cond, case.geom
    if act == 'region_add':
        regs = c.regions()
        if len(regs) < M.NM_RG:
            c.set_regions(regs + [('0.', '10.')])
    elif act == 'region_del':
        regs = c.regions()
        if regs:
            c.set_regions(regs[:-1])
    elif act in ('seg_add', 'seg_dup', 'seg_del', 'seg_up', 'seg_down'):
        k = int(arg) if arg not in (None, '') else None
        cur = g.segs[k] if (k is not None and 0 <= k < len(g.segs)) else None
        if act == 'seg_add':
            g.add_segment(cur)
            return (g.segs.index(cur) + 1) if cur else len(g.segs) - 1
        if act == 'seg_dup' and cur:
            g.add_segment(cur, list(cur.toks[:7]))
            return k + 1
        if act == 'seg_del' and cur:
            g.remove_segment(cur)
            return max(0, min(k, len(g.segs) - 1)) if g.segs else None
        if act in ('seg_up', 'seg_down') and cur:
            g.move_segment(cur, -1 if act == 'seg_up' else 1)
            return g.segs.index(cur)
    return None


# =====================================================================================
#  Inputs: page schema (mirrors pic_input_gui.py)
# =====================================================================================
def F(key, label, unit='', **kw):
    d = dict(key=key, label=label, unit=unit, type='text')
    d.update(kw)
    return d


def signed(key, label, check, unit='', on_when='positive', **kw):
    return F(key, label, unit, type='signed', check=check, on_when=on_when, **kw)


def yesno(key, label, **kw):
    return F(key, label, type='yesno', **kw)


def choice(key, label, options, match='first', **kw):
    return F(key, label, type='choice', options=[list(o) for o in options], match=match, **kw)


def sec(title, fields, desc=None, **kw):
    d = dict(title=title, fields=fields, desc=desc)
    d.update(kw)
    return d


B_TYPES = [('s', 'Not used (skip this entry)'), ('n', 'Analytical profile'), ('y', 'Field map from files')]
B_MODELS = [('gauss', 'Gaussian profile along x'), ('PG', 'PG-current filter field (ELISE fit)'),
            ('EEmagnets', 'Cusp (extraction) magnets'), ('thr', 'Hall-thruster magnet circuit')]
RF_OPTS = [('No', 'No RF'), ('Electrode', 'Electrode - RF voltage (not implemented in LePIC+)'),
           ('Wrapped', 'Wrapped - inductive coil around the driver wall'),
           ('Planar', 'Planar - flat coil under a window (model field)'),
           ('Coil', 'Coil - real coil turns from coil.inp (E and B)')]
INJ_OPTS = [('1', '1 - fixed current, uniform in the box'), ('2', '2 - replace ions lost at the walls'),
            ('3', '3 - electron beam along +z'), ('4', '4 - emissive cathode along z')]


def schema(case):
    c = case.cond
    p, _ = M.decode(case)
    mode = M.rf_mode(c.tok('rf'))
    rs = c.tok('restart')[:1].lower()
    dats = set()
    for d in (case.folder, os.path.dirname(case.folder)):
        try:
            dats |= {x for x in os.listdir(d) if x.lower().endswith('.dat')}
        except OSError:
            pass
    tabs = []
    tabs.append(dict(name='Simulation', sections=[
        sec('Plasma and physics', [
            F('rname', 'Reactions file', '(in input_dir/)', datalist=sorted(dats)),
            F('Te', 'Electron temperature  Te', 'eV'),
            F('n0', 'Plasma density  n0', 'm⁻³'),
            F('ngas', 'Neutral gas density', 'm⁻³'),
            signed('np_cell', 'Macro-particles per cell', 'Start without initial plasma', on_when='negative'),
            signed('k_eps0', 'Permittivity scaling  k_eps0', 'Negative (restart reloads Vgrd.bak)',
                   on_when='negative'),
            yesno('coulomb', 'Coulomb collisions'),
            choice('push', 'Push scheme', [('momentum', 'momentum - standard (momentum-conserving)'),
                                           ('energy', 'energy - EC-PIC (energy-conserving)')], match='exact')]),
        sec('Time', [F('nsteps', 'Number of time steps', '(run length; line 1)'),
                     F('tmax', 'Simulated time  t_max', 's   (read but not used by LePIC+)'),
                     F('kt', 'Time-step coefficient  kt', '(dt = kt * dx_min / vte)')]),
        sec('Grid and Poisson solver', [
            F('nx', 'Cells along x'), F('ny', 'Cells along y'), F('nz', 'Cells along z'),
            F('ng', 'Multigrid levels'), F('eps', 'Solver precision  eps'),
            F('omega', 'SOR relaxation  omega', '(1 ≤ ω < 2)')],
            'The box size is set in the Geometry tab. The number of cells is NOT an exponent.')]))
    inj = c.tok('opt_inj')
    tabs.append(dict(name='Heating & injection', sections=[
        sec('Heating / injection region  (cm)', [
            F('xl_pow', 'x from', 'cm'), F('xr_pow', 'x to', 'cm'), F('yl_pow', 'y from', 'cm'),
            signed('yr_pow', 'y to', 'Annular: heat only OUTSIDE the radius', 'cm', on_when='negative'),
            F('zl_pow', 'z from', 'cm'), F('zr_pow', 'z to', 'cm')],
            'Electron heating uses the x-range; particle injection uses the whole box. '
            'Planar/Coil antennas: the window is at x = "x from".'),
        sec('Initial plasma', [F('x_load', 'Load plasma for x ≤', 'cm')]),
        sec('Electron heating', [
            signed('Pabs', 'Absorbed power  Pabs', 'Power control on (also the RF antenna target)', 'W'),
            signed('nu_h', 'Heating frequency  nu_h', 'Heating on', 's⁻¹'),
            signed('delta_heat', 'Skin depth  delta', 'Self-consistent (c/wpe every 1/f0)',
                   'cm   (0 = uniform heating)', on_when='negative')]),
        sec('Particle injection', [
            signed('I_inj', 'Injected current  I_inj', 'Injection on', 'A'),
            choice('opt_inj.mode', 'Mode', INJ_OPTS, match='exact', ikey='opt_inj', value=M.strip_sign(inj)),
            yesno('opt_inj.cos', 'Cosine profile instead of uniform', ikey='opt_inj', value=M.is_neg(inj))])]))
    bsecs = [sec('Magnetic field', [
        F('nB.count', 'Number of entries', '(0 - %d)' % M.MAX_B, ikey='nB', value=M.strip_sign(c.tok('nB'))),
        yesno('nB.neg', 'Positive ions are not magnetized', ikey='nB', value=M.is_neg(c.tok('nB'))),
        F('nBx', 'Map grid  nx', 'cells'), F('nBy', 'Map grid  ny', 'cells'), F('nBz', 'Map grid  nz', 'cells')],
        'The fields of all entries are added together.')]
    for k, r in enumerate(c.bfiles):
        fm, ch = r.get(0)[:1].lower(), r.get(1)[:1].lower()
        show = {'file'}
        if fm == 'y':
            show |= {'prefix', 'scale'}
        elif fm == 'n':
            show |= {'model'} | {'g': {'B0', 'dir', 'x0', 'dL'}, 'p': {'B0', 'dir'}, 'e': {'B0', 'x0', 'z0', 'dL'},
                                 't': {'B0', 'x0', 'dL'}}.get(ch, {'B0', 'dir', 'x0', 'z0', 'dL'})
        fl = [choice('B%d.file' % k, 'Type', B_TYPES),
              choice('B%d.name' % k, 'Model', B_MODELS, hidden='model' not in show),
              F('B%d.name' % k, 'File name prefix', 'reads <prefix>Bx/By/Bz.dat in mag_field_files/',
                hidden='prefix' not in show, fid='B%d.prefix' % k),
              F('B%d.scale' % k, 'Scale factor', hidden='scale' not in show),
              F('B%d.B0' % k, 'Amplitude  B0', 'G', hidden='B0' not in show),
              choice('B%d.dir' % k, 'Direction', [('Bx', 'Bx'), ('By', 'By'), ('Bz', 'Bz')], match='exact',
                     hidden='dir' not in show),
              F('B%d.x0' % k, 'Position  x0', 'cm', hidden='x0' not in show),
              F('B%d.z0' % k, 'Position  z0', 'cm', hidden='z0' not in show),
              F('B%d.dL' % k, 'Width / length  dL', 'cm', hidden='dL' not in show)]
        bsecs.append(sec('Entry %d' % (k + 1), fl))
    tabs.append(dict(name='Magnetic field', sections=bsecs))
    rfshow = {'Electrode': {'phi0', 'f0', 'phi1', 'f1'}, 'Wrapped': {'f0', 'E0'},
              'Planar': {'f0', 'E0', 'R_win'}, 'Coil': {'f0', 'E0', 'R_win'}}.get(mode, set())
    esecs = [
        sec('Grid apertures (holes)', [yesno('holes', 'Drill apertures in a wall'),
                                       F('Lhy', 'Aperture size  Ly', 'cm'), F('Lhz', 'Aperture size  Lz', 'cm'),
                                       F('nhy', 'Number of apertures along y'), F('nhz', 'Number of apertures along z'),
                                       F('ind_g', 'Label of the wall / grid')]),
        sec('Surface emission of negative ions', [F('jne', 'Negative-ion flux  jH-', 'mA/cm²   (0 = off)'),
                                                  F('num_grd', 'Wall label for impact statistics')]),
        sec('Dielectric walls', [F('Ca', 'Surface capacitance  Ca/S', 'F/m²')]),
        sec('Secondary electron emission / emissive cathode', [
            signed('gam_sec', 'Emission coefficient  gamma', 'Enabled', 'electrons per ion'),
            F('igrid_sec', 'Emitting electrode label'),
            F('THm', 'Emission temperature / energy  TH', 'eV   (also e-beam and negative ions)')]),
        sec('RF', [choice('rf', 'Mode', RF_OPTS, match='exact', value=mode or c.tok('rf')),
                   F('phi0', 'Voltage amplitude 1', 'V', hidden='phi0' not in rfshow),
                   F('f0', 'Frequency  f0', 'Hz', hidden='f0' not in rfshow),
                   F('phi1', 'Voltage amplitude 2', 'V', hidden='phi1' not in rfshow),
                   F('f1', 'Frequency 2', 'Hz', hidden='f1' not in rfshow),
                   F('E0', 'Initial drive amplitude', 'A   (coil current, peak)' if mode == 'Coil' else 'V/m',
                     hidden='E0' not in rfshow),
                   F('R_win', 'Window radius  R_win', 'cm   (0 = y_max - ymax/2)', hidden='R_win' not in rfshow)])]
    if mode == 'Coil':
        fl = [F('coil.n', 'Number of turns')]
        for j in range(len(case.coil.turns)):
            fl += [F('coil%d.0' % j, 'Turn %d: radius a' % (j + 1), 'cm'),
                   F('coil%d.1' % j, 'Turn %d: position x' % (j + 1), 'cm'),
                   F('coil%d.2' % j, 'Turn %d: relative current w' % (j + 1))]
        esecs.append(sec('Coil turns  (coil.inp)', fl, 'Concentric loops coaxial with x. Every turn must be behind '
                         'the window: x < x_min of the heating region (clamped to >= 0).'))
    tabs.append(dict(name='Emission & RF', sections=esecs))
    rg = [F('rg%d.%d' % (i, j), 'Region %d: x %s' % (i + 1, nm), 'cm', ikey='pdf')
          for i in range(len(c.regions())) for j, nm in ((0, 'from'), (1, 'to'))]
    tabs.append(dict(name='Output', sections=[
        sec('Saving data', [
            signed('nsav', 'Save every', 'Keep accumulating averages', 'time steps', on_when='negative'),
            signed('nbak', 'Backup every', 'Do not time-average the saved data', 'time steps', on_when='negative')]),
        sec('Sequence of profile snapshots', [signed('cnt_plt', 'Number of snapshots', 'Enabled'),
                                              F('t1', 'Window start  t1', 's'), F('t2', 'Window end  t2', 's')]),
        sec('Velocity-distribution output and diagnostic regions', [
            yesno('pdf', 'Write velocity distributions'), F('ptype_pdf', 'Species number', '(1 = electrons)'),
            F('yl_rg', 'y from', 'cm'), F('yr_rg', 'y to', 'cm'), F('zl_rg', 'z from', 'cm'), F('zr_rg', 'z to', 'cm')]
            + rg, buttons=[['region_add', '+ Add region'], ['region_del', '− Remove last']]),
        sec('Restart', [
            choice('restart', 'Run type', [('no', 'New run'), ('yes', 'Continue (same number of MPI ranks)'),
                                           ('external', 'Import backup from a different layout')]),
            F('omp_max', 'OpenMP threads of the backup', hidden=rs != 'e'),
            F('mpi_max', 'MPI ranks of the backup', hidden=rs != 'e'),
            F('np_dup', 'Particle duplication / removal', '(1 = unchanged)', hidden=rs not in ('y', 'e'))],
            'These numbers do not set the threads/processes of this run (runModular.sh does, and rewrites '
            'this line).')]))
    # potentials with usage
    box, segs, _ = M.decode_geometry(case)
    use = {}
    for k, s in enumerate(segs):
        if None in s:
            continue
        kind = 'solid block' if s[6] < 0 else 'wall'
        die = ', dielectric' if (s[3] < 0 or s[4] < 0) else ''
        use.setdefault(abs(s[6]), []).append('segment %d (%s%s)' % (k + 1, kind, die))
    pots = []
    for k in range(len(case.bnd.vals)):
        lab = k + 1
        txt = '; '.join(use.get(lab, [])) or 'not used in geometry'
        role = []
        if p['igrid_sec'] == lab and ((p['gam_sec'] or 0) > 0 or abs(p['opt_inj'] or 0) == 4):
            role.append('emitting electrode')
        if p['ind_g'] == lab and (p['holes'] or (p['jne'] or 0) > 0):
            role.append('grid / plasma electrode')
        pots.append(F('b.V%d' % lab, 'Label %d' % lab, 'V', note=txt + ('   [' + ', '.join(role) + ']' if role else '')))
    tabs.append(dict(name='Geometry', geometry=True, sections=[sec('Box size (cm)', [
        signed('g.box.0', 'x size', 'Zero-gradient (Neumann) boundary at x = 0', 'cm', on_when='negative',
               ikey='g.xmax'),
        signed('g.box.1', 'y size', 'Draw wall openings round (cylinders along x)', 'cm', on_when='negative',
               ikey='g.ymax'),
        signed('g.box.2', 'z size', 'Draw round objects along z', 'cm', on_when='negative', ikey='g.zmax')])]))
    tabs.append(dict(name='Potentials', sections=[
        sec('Electrodes and walls (boundary.inp)', [
            signed('ngrid', 'Number of wall labels', 'No 3D time averages (negative)', on_when='negative')],
            'One potential per wall label. Label 0 is reserved for periodic faces and has no potential.'),
        sec('Potentials', pots)]))
    # fill values
    for t in tabs:
        for s in t['sections']:
            for f in s['fields']:
                f.setdefault('ikey', f['key'])
                f.setdefault('fid', f['key'])
                f['help'] = help_key(f['ikey'])
                if 'value' in f:
                    continue
                try:
                    f['value'] = tok_get(case, f['key'])
                except (KeyError, IndexError):
                    f['value'] = ''
    return tabs


def segments(case):
    out = []
    for k, r in enumerate(case.geom.segs):
        t = r.toks + [''] * (7 - len(r.toks))
        kind = seg_kind(r)
        fl = [nm for i, nm in ((3, 'diel. x'), (4, 'diel. y'), (5, 'periodic z'), (1, 'round toggled'))
              if M.is_neg(t[i])]
        out.append(dict(
            label=M.strip_sign(t[6]), kind={'+': 'wall + opening', '-': 'solid block', '0': 'periodic'}[kind],
            kindtok=kind, x='%s – %s' % (t[0], M.strip_sign(t[3])), y='%s – %s' % (M.strip_sign(t[1]), M.strip_sign(t[4])),
            z='%s – %s' % (t[2], M.strip_sign(t[5])), flags=', '.join(fl),
            coords=[M.strip_sign(t[i]) for i in range(6)], flagvals={str(i): M.is_neg(t[i]) for i in (1, 3, 4, 5)}))
    return out


def derived_lines(case):
    p, _ = M.decode(case)
    box, _, _ = M.decode_geometry(case)
    d = M.derived(p, box) if None not in box else None
    if not d:
        return ['(fill in valid numbers to see derived quantities)']
    n3 = (p['nx'], p['ny'], p['nz'])
    L = [abs(b) for b in box]
    lines = [
        'Grid        %d x %d x %d = %.2f M cells      box %.3g x %.3g x %.3g cm      multigrid levels %s' % (
            n3[0], n3[1], n3[2], d['ncells'] / 1e6, L[0], L[1], L[2], p['ng']),
        'Cell size   dx, dy, dz = %.1f, %.1f, %.1f µm' % tuple(h * 1e6 for h in d['h']),
        'Electrons   vte = %.3g m/s      Debye length = %.3g µm      dx/λD = %.2f' % (
            d['vt'], d['lde'] * 1e6, d['h'][0] / d['lde']),
        'Time step   dt = %.3g s      ωpe·dt = %.3g      (ωpe = %.3g rad/s)' % (d['dt'], d['wpdt'], d['wp']),
        'Run length  %s steps × dt = %s      (t_max would need %s steps; not used)      heating every %d steps%s' % (
            p['nsteps'] if p['nsteps'] else '?', ('%.3g s' % d['t_end']) if 't_end' in d else '?',
            ('%.3g' % d['nsteps_tmax']) if 'nsteps_tmax' in d else '?', d['ns_heat'],
            ('      ν_h·dt_heat = %.3g (must be ≤ 1)' % d['nudt']) if 'nudt' in d else '')]
    if p['antenna']:
        lines.append('RF antenna  %s      skin depth c/ωpe(n0) = %.3g cm (true ε0)      %s' % (
            p['rf_mode'], d['skin'] * 100,
            ('update every 1/(f0 dt) = %.0f steps' % d['ns_RF']) if 'ns_RF' in d else 'f0 must be > 0'))
    lines.append('MPI         %s   (nz = %s)' % (M.mpi_summary(case.mpi), n3[2]))
    return lines


def list_cases():
    base = os.path.join(ROOT, 'input_dir')
    out = []
    if all(os.path.isfile(os.path.join(base, f)) for f in M.Case.FILES):
        out.append(dict(name='input_dir  (what LePIC+ runs)', folder=rel(base)))
    try:
        for d in sorted(os.listdir(base)):
            pth = os.path.join(base, d)
            if os.path.isdir(pth) and all(os.path.isfile(os.path.join(pth, f)) for f in M.Case.FILES):
                out.append(dict(name=d, folder=rel(pth)))
    except OSError:
        pass
    return out


def editor_state(case, sel=None):
    issues = M.validate(case)
    return dict(
        folder=rel(case.folder), modified=case.modified(), tabs=schema(case), segments=segments(case),
        issues=[dict(level=i.level, msg=i.msg, key=i.key, tab=i.tab) for i in issues],
        derived=derived_lines(case), sel=sel, can_activate=case.input_dir() is not None,
        files=[os.path.basename(f.path) for f in case.files])


# ---- geometry preview ---------------------------------------------------------------
_GEOM_CACHE = {}


def geometry_raster(case):
    p, _ = M.decode(case)
    box, segs, _ = M.decode_geometry(case)
    if None in box or 0 in box:
        raise ValueError('Fill in a valid box size to see the preview.')
    good = [s for s in segs if None not in s]
    n = [max(4, p['nx'] or 64), max(4, p['ny'] or 64), max(4, p['nz'] or 64)]
    N = n[0] * n[1] * n[2]
    if N > 3e6:
        f = (3e6 / N) ** (1 / 3)
        n = [max(4, int(v * f)) for v in n]
    holes = None
    if p['holes'] and None not in (p['Lhy'], p['Lhz'], p['nhy'], p['nhz'], p['ind_g']):
        holes = dict(on=True, Lhy=p['Lhy'], Lhz=p['Lhz'], nhy=p['nhy'], nhz=p['nhz'], ind_g=p['ind_g'])
    sig = repr((box, good, n, holes))
    if sig not in _GEOM_CACHE:
        _GEOM_CACHE.clear()
        _GEOM_CACHE[sig] = M.rasterize(box, good, n, holes)
    return _GEOM_CACHE[sig]


def bc_color(v, dtype):
    from matplotlib.colors import to_rgb
    if v == -1:
        return to_rgb('#cfe8ff')
    if v == 0:
        return to_rgb('#d7d7d7')
    if v == -2:
        return to_rgb('#e0d0f0')
    c = np.array(to_rgb(PALETTE[(v - 1) % len(PALETTE)]))
    return tuple(c * 0.55 + 0.45) if dtype.get(v, 1) > 1 else tuple(c)


def geometry_png(case, frac):
    """Three cross-sections through the box at the fractions frac = {xy: fz, xz: fy, yz: fx}."""
    r = geometry_raster(case)
    vol, L, dtype = r['vol'], r['L'], r['dtype']
    fig = Figure(figsize=(11, 3.6), dpi=90, layout='constrained')
    planes = (('xy', 2, 'x', 'y', 0, 1), ('xz', 1, 'x', 'z', 0, 2), ('yz', 0, 'y', 'z', 1, 2))
    labels = sorted(int(v) for v in np.unique(vol))
    lut = {v: bc_color(v, dtype) for v in labels}
    for i, (pl, ax_n, hn, vn, ha, va) in enumerate(planes):
        f = frac.get(pl, 0.5)
        idx = int(round(f * (vol.shape[ax_n] - 1)))
        sl = np.take(vol, idx, axis=ax_n)            # [h, v]
        img = np.zeros(sl.shape + (3,))
        for v, col in lut.items():
            img[sl == v] = col
        ax = fig.add_subplot(1, 3, i + 1)
        ax.imshow(img.transpose(1, 0, 2), origin='lower', extent=(0, L[ha], 0, L[va]), interpolation='nearest')
        ax.set_title('%s plane at %s = %.3g cm' % (pl, 'zyx'[i], f * L[ax_n]), fontsize=10)
        ax.set_xlabel('%s (cm)' % hn)
        ax.set_ylabel('%s (cm)' % vn)
        ax.set_aspect('equal')
    names = []
    for v in labels:
        if v > 0:
            names.append((lut[v], 'label %d%s' % (v, '' if dtype.get(v, 1) == 1 else ' (dielectric)')))
    names += [(lut[-1], 'plasma'), (lut.get(0, bc_color(0, {})), 'periodic / unset')]
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c, label=t) for c, t in names], loc='outside lower center',
               ncol=min(8, len(names)), fontsize=9, frameon=False)
    return fig_png(fig), dict(notes=r['notes'])


def fig_png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    return base64.b64encode(buf.getvalue()).decode('ascii')


# =====================================================================================
#  Output 2D
# =====================================================================================
def out_dir(d):
    p = inside_root(d or 'Output')
    if os.path.basename(p) in V.MAPS_DIRS:
        p = os.path.dirname(p)
    return p


def geometry_of(outp):
    for cand in (os.path.join(outp, 'geometry.inp'), os.path.join(os.path.dirname(outp), 'input_dir', 'geometry.inp')):
        if os.path.isfile(cand):
            try:
                lx, ly, lz, nmn = V.parse_geometry3d(cand)
                return dict(L=(lx, ly, lz), nmn=nmn, path=rel(cand))
            except Exception:          # noqa
                pass
    return dict(L=None, nmn=False, path=None)


def snap_arg(s):
    if not s or s == 'latest':
        return None
    return V.INITIAL if s == V.INITIAL else int(s)


def out_scan(d, snap):
    outp = out_dir(d)
    entries, names = V.scan(outp, snap_arg(snap))
    d2 = V.maps_dir(outp)
    snaps = V.snapshots(d2) if d2 is not None else []
    res = []
    for e in entries:
        item = dict(key=e.key, group=e.group, label=e.label, kind=e.kind,
                    planes=sorted((e.extra.get('planes') or {}).keys(), key=V.PLANES.index))
        if e.kind == 'iw':
            try:
                item['ncols'] = V.read_table(e.files[0]).shape[1] - 1
            except Exception:          # noqa
                item['ncols'] = 0
        res.append(item)
    return dict(dir=rel(outp), entries=res, snapshots=[str(s) for s in snaps], names=names,
                geometry=geometry_of(outp), maps=d2 is not None)


def _entry(d, snap, key):
    outp = out_dir(d)
    entries, names = V.scan(outp, snap_arg(snap))
    for e in entries:
        if e.key == key:
            return outp, e, names
    raise KeyError('No dataset %s (in this snapshot?)' % key)


def _tex(sym):
    return '$%s$' % sym if sym else ''


def out_plot(q):
    d, snap, key = q.get('dir'), q.get('snap'), q.get('key')
    plane = q.get('plane') or 'xy'
    mode = q.get('mode') or 'map'
    cmap = q.get('cmap') or 'viridis'
    logc = q.get('log') == '1'
    outp, e, names = _entry(d, snap, key)
    geo = geometry_of(outp)
    fig = Figure(figsize=(float(q.get('w') or 9), float(q.get('h') or 5.2)), dpi=90, layout='constrained')
    ax = fig.add_subplot(111)
    meta = dict(kind='curve')
    tlab = r'$t$ ($\mu$s)'
    if e.kind == 'field' or e.kind == 'flux':
        pf = e.extra.get('planes')
        f = (pf or {}).get(plane) if pf else e.files[0]
        if f is None:
            raise KeyError('%s is not saved in the %s plane' % (e.label, plane))
        f = f[0] if isinstance(f, list) else f
        nx, ny, Z, th = V.read_map(f)
        ha, va = V.PLANE_AX.get(plane, ('x', 'y'))
        L = geo['L']
        if L:
            Lh, Lv = L['xyz'.index(ha)], L['xyz'.index(va)]
            xs, ys, un = np.linspace(0, Lh, nx + 1), np.linspace(0, Lv, ny + 1), 'cm'
        else:
            xs, ys, un = np.arange(nx + 1.0), np.arange(ny + 1.0), 'node'
        title = '%s %s  (%s plane)' % (e.desc, _tex(e.sym), plane)
        unit = e.unit
        if mode == 'map':
            from matplotlib.colors import LogNorm, Normalize
            vmin = float(q['vmin']) if q.get('vmin') else None
            vmax = float(q['vmax']) if q.get('vmax') else None
            Zp = np.where(Z > 0, Z, np.nan) if logc else Z
            norm = LogNorm(vmin=vmin, vmax=vmax) if logc else Normalize(vmin=vmin, vmax=vmax)
            im = ax.pcolormesh(xs, ys, Zp, cmap=cmap, norm=norm, shading='nearest')
            cb = fig.colorbar(im, ax=ax)
            cb.set_label('%s (%s)' % (_tex(e.sym), unit) if unit else _tex(e.sym))
            if th is not None and e.extra.get('vec') and q.get('stream', '1') == '1':
                U, W = np.cos(th), np.sin(th)
                mag = np.where(np.isfinite(Z), np.abs(Z), 0)
                s = max(1, int(math.ceil(max(Z.shape) / 64)))      # streamlines on a <= 64-point grid (speed)
                try:
                    ax.streamplot(xs[::s], ys[::s], (U * (mag > 0))[::s, ::s], (W * (mag > 0))[::s, ::s],
                                  color='white', linewidth=0.7, density=1.0, arrowsize=0.7)
                except Exception:          # noqa
                    pass
            ax.set_aspect('equal' if L else 'auto')
            ax.set_xlim(xs[0], xs[-1])
            ax.set_ylim(ys[0], ys[-1])
            ax.set_xlabel('%s (%s)' % (ha, un))
            ax.set_ylabel('%s (%s)' % (va, un))
            step = max(1, int(math.ceil(max(Z.shape) / 300)))
            meta = dict(kind='map', x=[float(xs[0]), float(xs[-1])], y=[float(ys[0]), float(ys[-1])],
                        z=np.nan_to_num(Z[::step, ::step]).round(6).tolist(), step=step, nx=nx, ny=ny,
                        hn=ha, vn=va, unit=un, zunit=unit)
        else:
            if mode == 'profx':
                yv = Z[ny // 2, :]
                xv, xn = xs, ha
                where = '%s = %.3g %s' % (va, ys[ny // 2], un)
            else:
                yv = Z[:, nx // 2]
                xv, xn = ys, va
                where = '%s = %.3g %s' % (ha, xs[nx // 2], un)
            ax.plot(xv, yv, lw=1.4)
            if logc:
                ax.set_yscale('log')
            ax.set_xlabel('%s (%s)' % (xn, un))
            ax.set_ylabel('%s (%s)' % (_tex(e.sym), unit) if unit else _tex(e.sym))
            title += ', profile at ' + where
            meta = dict(kind='curve', x=np.asarray(xv).tolist(), ys=[np.asarray(yv).tolist()], xn=xn, xu=un)
        ax.set_title(title, fontsize=11)
    else:
        tab = V.dedupe_time(V.read_table(e.files[0]))
        t = tab[:, 0] * 1e6
        series = []
        if e.kind == 'iw':
            first = 0 if geo['nmn'] else 1
            sel = [int(c) for c in (q.get('cols') or '').split(',') if c != '']
            ncol = tab.shape[1] - 1
            sel = [c for c in sel if 1 <= c <= ncol] or list(range(1, ncol + 1))
            tot = np.zeros_like(t)
            for c in sel:
                ax.plot(t, tab[:, c], lw=0.9, alpha=0.75, label='wall label %d' % (c - 1 + first))
                tot += tab[:, c]
                series.append(tab[:, c].tolist())
            if len(sel) > 1:
                ax.plot(t, tot, 'k--', lw=1.4, label='sum')
                series.append(tot.tolist())
            pre = e.extra.get('pref')
            ax.set_ylabel('$%s$ (%s)' % ('I_w' if pre == 'Iw' else 'P_w', 'A' if pre == 'Iw' else
                                          ('W' if V.maps_dir(outp) is not None else 'eV/s')))
            ax.legend(fontsize=8)
            ax.set_title('%s  -  %s' % (e.desc, names[int(key.split(':')[1]) - 1]), fontsize=11)
        elif e.kind == 'col':
            c = e.extra['col']
            if tab.shape[1] <= c:
                raise IOError('%s has %d columns, column %d requested' % (e.files[0].name, tab.shape[1], c + 1))
            x = t if e.extra.get('time', True) else np.arange(1, len(tab) + 1)
            ax.plot(x, tab[:, c], lw=1.2)
            series.append(tab[:, c].tolist())
            ax.set_ylabel('%s (%s)' % (_tex(e.sym), e.unit) if e.unit else _tex(e.sym))
            ax.set_title(e.desc, fontsize=11)
            if not e.extra.get('time', True):
                t, tlab = x, 'output #'
        else:
            if e.kind == 'probe_ne':
                j = max(0, min(int(q.get('probe') or 0), tab.shape[1] - 4))
                y = tab[:, 3 + j]
                frac = V.NE_FRACS[j] if j < len(V.NE_FRACS) else None
                ax.set_title('Electron density at x = %s xmax' % ('%.1f' % frac if frac is not None else '?'),
                             fontsize=11)
            else:
                y = tab[:, 1 if e.kind == 'probe_phi' else 2]
                ax.set_title('%s at the box centre' % e.desc, fontsize=11)
            ax.plot(t, y, lw=1.2)
            series.append(y.tolist())
            ax.set_ylabel('%s (%s)' % (_tex(e.sym), e.unit))
        if logc:
            ax.set_yscale('log')
        ax.set_xlabel(tlab)
        ax.grid(alpha=0.3)
        meta = dict(kind='curve', x=np.asarray(t).tolist(), ys=series, xn='t' if tlab != 'output #' else '#',
                    xu='µs' if tlab != 'output #' else '')
    png = fig_png(fig)
    fig.canvas.draw()
    bb = ax.get_position()
    meta['axes'] = [bb.x0, 1 - bb.y1, bb.width, bb.height]       # fractions of the image, from top-left
    meta['xlim'] = list(ax.get_xlim())
    meta['ylim'] = list(ax.get_ylim())
    meta['ylog'] = ax.get_yscale() == 'log'
    return dict(png=png, meta=meta)


# =====================================================================================
#  Output 3D  (Output/phi_n_3D.dat)
# =====================================================================================
_P3 = {}


def load_phi_n(path):
    st = os.stat(path)
    sig = (path, st.st_mtime_ns, st.st_size)
    if _P3.get('sig') == sig:
        return _P3['data']
    with open(path, 'rb') as f:
        b = f.read(20)
        en = '<' if struct.unpack('<i', b[:4])[0] == 12 else '>'
        if struct.unpack(en + 'i', b[:4])[0] != 12:
            raise ValueError('Not an Output/phi_n_3D.dat file (first record is not nx, ny, nz).')
        nx, ny, nz = struct.unpack(en + '3i', b[4:16])
        expected = 16 * (nx + 1) * (ny + 1) * (nz + 1)
        parts, tot, pos = [], 0, 20
        while tot < expected:                 # sub-records (> 2 GB) have negative markers
            f.seek(pos)
            L = abs(struct.unpack(en + 'i', f.read(4))[0])
            parts.append(np.fromfile(f, dtype=en + 'f8', count=L // 8))
            tot += L
            pos += 8 + L
    arr = np.concatenate(parts).reshape((2, nx + 1, ny + 1, nz + 1), order='F')
    data = dict(n=(nx, ny, nz), phi=arr[0], dens=arr[1])
    _P3.clear()
    _P3.update(sig=sig, data=data)
    return data


def p3_info(q):
    path = inside_root(q.get('file') or 'Output/phi_n_3D.dat')
    dat = load_phi_n(path)
    outp = os.path.dirname(path)
    geo = geometry_of(outp)
    names = V.parse_particles_out(os.path.join(outp, 'particles.out'))
    ncharged = 1 + sum(1 for nm in names[1:] if nm.rstrip(']').endswith(('+', '-')))
    k = min(2, max(1, ncharged))
    return dict(file=rel(path), n=dat['n'], L=geo['L'], geometry=geo['path'],
                species=names[k - 1] if k <= len(names) else 'ptype %d' % k, ptype=k)


def p3_plot(q):
    info = p3_info(q)
    dat = load_phi_n(inside_root(info['file']))
    n = dat['n']
    L = info['L'] or [float(v) for v in n]
    un = 'cm' if info['L'] else 'node'
    axis = 'xyz'.index(q.get('axis') or 'z')
    idx = max(0, min(int(q.get('index') if q.get('index') not in (None, '') else n[axis] // 2), n[axis]))
    qty = q.get('qty') or 'phi'
    h = [L[a] / n[a] * 1e-2 if info['L'] else 1.0 for a in range(3)]
    if qty in ('phi', 'dens'):
        A = np.take(dat[qty], idx, axis=axis)
    else:
        c = 'xyz'.index(qty[1])
        sl = [slice(None)] * 3
        lo, hi = max(0, idx - 1), min(n[axis], idx + 1)
        sl[axis] = slice(lo, hi + 1)
        sub = dat['phi'][tuple(sl)]
        Eall = -np.gradient(sub, h[c], axis=c) if sub.shape[c] > 1 else np.zeros_like(sub)
        A = np.take(Eall, idx - lo, axis=axis)
        A = np.where(np.take(dat['dens'], idx, axis=axis) > 0, A, 0.0)
    rem = [a for a in range(3) if a != axis]
    xs = np.linspace(0, L[rem[0]], n[rem[0]] + 1)
    ys = np.linspace(0, L[rem[1]], n[rem[1]] + 1)
    fig = Figure(figsize=(float(q.get('w') or 9), float(q.get('h') or 5.2)), dpi=90, layout='constrained')
    ax = fig.add_subplot(111)
    logc = q.get('log') == '1' and qty == 'dens'
    from matplotlib.colors import LogNorm
    Z = A.T
    im = ax.pcolormesh(xs, ys, np.where(Z > 0, Z, np.nan) if logc else Z, cmap=q.get('cmap') or 'viridis',
                       shading='nearest', norm=LogNorm() if logc else None)
    sym, unit = {'phi': (r'\phi', 'V'), 'dens': ('n_{%s}' % info['species'].strip('[]'), r'm$^{-3}$'),
                 'Ex': ('E_x', 'V/m'), 'Ey': ('E_y', 'V/m'), 'Ez': ('E_z', 'V/m')}[qty]
    fig.colorbar(im, ax=ax).set_label('$%s$ (%s)' % (sym, unit))
    pos = idx * L[axis] / n[axis]
    ax.set_title('$%s$, %s plane at %s = %.4g %s (node %d/%d), time-averaged' % (
        sym, ''.join('xyz'[a] for a in rem), 'xyz'[axis], pos, un, idx + 1, n[axis] + 1), fontsize=11)
    ax.set_xlabel('%s (%s)' % ('xyz'[rem[0]], un))
    ax.set_ylabel('%s (%s)' % ('xyz'[rem[1]], un))
    ax.set_aspect('equal' if info['L'] else 'auto')
    png = fig_png(fig)
    fig.canvas.draw()
    bb = ax.get_position()
    step = max(1, int(math.ceil(max(Z.shape) / 300)))
    meta = dict(kind='map', axes=[bb.x0, 1 - bb.y1, bb.width, bb.height], xlim=list(ax.get_xlim()),
                ylim=list(ax.get_ylim()), x=[0.0, float(L[rem[0]])], y=[0.0, float(L[rem[1]])],
                z=np.nan_to_num(Z[::step, ::step]).tolist(), nx=n[rem[0]], ny=n[rem[1]], step=step,
                hn='xyz'[rem[0]], vn='xyz'[rem[1]], unit=un, zunit=unit.replace('$', '').replace('^{-3}', '⁻³'))
    return dict(png=png, meta=meta, info=info, index=idx)


# =====================================================================================
#  HTTP
# =====================================================================================
class Handler(BaseHTTPRequestHandler):
    server_version = 'LePICweb/1.0'

    def log_message(self, fmt, *args):          # keep the terminal quiet
        pass

    def _send(self, code, body, ctype='application/json'):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, allow_nan=False, default=str).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(body)

    def _json_safe(self, obj):
        return json.loads(json.dumps(obj, default=str).replace('NaN', 'null').replace('Infinity', 'null'))

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = {k: v[-1] for k, v in urllib.parse.parse_qs(u.query).items()}
        try:
            if u.path in ('/', '/index.html'):
                with open(os.path.join(HERE, 'web', 'index.html'), 'rb') as f:
                    return self._send(200, f.read(), 'text/html; charset=utf-8')
            if u.path == '/app.js':
                with open(os.path.join(HERE, 'web', 'app.js'), 'rb') as f:
                    return self._send(200, f.read(), 'text/javascript; charset=utf-8')
            if u.path == '/api/help':
                return self._send(200, {k: list(v) for k, v in HELP.items()})
            if u.path == '/api/cases':
                return self._send(200, dict(root=ROOT, cases=list_cases()))
            with LOCK:
                case = STATE['case']
                if u.path == '/api/state':
                    if case is None:
                        return self._send(200, dict(folder=None))
                    return self._send(200, editor_state(case))
                if u.path == '/api/diff':
                    return self._send(200, dict(lines=case.diff() if case else []))
                if u.path == '/api/geometry':
                    frac = {k: float(q[k]) for k in ('xy', 'xz', 'yz') if q.get(k)}
                    png, meta = geometry_png(case, frac)
                    return self._send(200, dict(png=png, meta=meta))
            if u.path == '/api/out/scan':
                return self._send(200, self._json_safe(out_scan(q.get('dir'), q.get('snap'))))
            if u.path == '/api/out/plot':
                return self._send(200, self._json_safe(out_plot(q)))
            if u.path == '/api/3d/plot':
                return self._send(200, self._json_safe(p3_plot(q)))
            return self._send(404, dict(error='not found'))
        except Exception as ex:          # noqa
            traceback.print_exc()
            return self._send(400, dict(error=str(ex)))

    def do_POST(self):
        u = urllib.parse.urlparse(self.path)
        try:
            n = int(self.headers.get('Content-Length') or 0)
            body = json.loads(self.rfile.read(n) or b'{}')
            with LOCK:
                case = STATE['case']
                if u.path == '/api/load':
                    STATE['case'] = M.Case(inside_root(body['folder']))
                    return self._send(200, editor_state(STATE['case']))
                if case is None:
                    raise ValueError('No case loaded.')
                if u.path == '/api/set':
                    tok_set(case, body['key'], body.get('value', ''))
                    return self._send(200, editor_state(case, body.get('sel')))
                if u.path == '/api/action':
                    sel = do_action(case, body['action'], body.get('arg'))
                    return self._send(200, editor_state(case, sel))
                if u.path == '/api/save':
                    issues = M.validate(case)
                    fatal = [i.msg for i in issues if i.level == 'fatal']
                    if fatal:
                        raise ValueError('Cannot save, please fix: ' + ' | '.join(fatal[:5]))
                    written = case.save()
                    STATE['case'] = M.Case(case.folder)
                    st = editor_state(STATE['case'])
                    st['message'] = 'Saved %s (previous versions in gui_backups/)' % (
                        ', '.join(os.path.basename(w) for w in written) or 'nothing - files already up to date')
                    return self._send(200, st)
                if u.path == '/api/activate':
                    if case.modified():
                        raise ValueError('Save the case first.')
                    names = case.activate()
                    st = editor_state(case)
                    st['message'] = 'Copied %s to input_dir/ (previous files in input_dir/gui_backups/)' % \
                        ', '.join(names)
                    return self._send(200, st)
                if u.path == '/api/reload':
                    STATE['case'] = M.Case(case.folder)
                    return self._send(200, editor_state(STATE['case']))
            return self._send(404, dict(error='not found'))
        except Exception as ex:          # noqa
            traceback.print_exc()
            return self._send(400, dict(error=str(ex)))


def main():
    global ROOT
    ap = argparse.ArgumentParser(description='Browser GUI for LePIC+ (input editor and Output viewers)')
    ap.add_argument('--port', type=int, default=8765)
    ap.add_argument('--root', default=ROOT, help='LePIC+ folder (contains input_dir/ and Output/)')
    a = ap.parse_args()
    ROOT = os.path.abspath(a.root)
    srv = ThreadingHTTPServer(('127.0.0.1', a.port), Handler)
    print('LePIC+ web GUI for %s' % ROOT)
    print('Open  http://localhost:%d  (VS Code forwards the port; from plain ssh: ssh -L %d:localhost:%d <host>)'
          % (a.port, a.port, a.port))
    print('Ctrl-C to stop.')
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    main()
