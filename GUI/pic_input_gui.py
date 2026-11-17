#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pic_input_gui.py  --  graphical editor for the LePIC+ input files
=================================================================

    python3 GUI/pic_input_gui.py [path/to/input_dir]

Edits  conditions.inp, geometry.inp, boundary.inp  (and coil.inp for the Coil
RF antenna) with plain-language labels, units, contextual help and live checks
that mirror what LePIC+ (mod_readConditions.f90, mod_boundary.f90,
mod_generateBoundary.f90, mod_state.f90) does.  Adapted from the LePIC3D_DD
3.4.1 editor.  Needs only the Python
standard library (tkinter); numpy is optional (geometry preview).

Keep this file next to  pic_inputs.py  (the model layer).
"""
import os
import re
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter import font as tkfont

import pic_inputs as M

try:
    import numpy as np
except Exception:          # preview is optional
    np = None

from gui_help import HELP, help_key  # noqa: E402  (help texts shared with web_gui.py)


LEVEL_COLOR = {'fatal': '#b00020', 'error': '#b00020', 'warning': '#b36b00', 'info': '#1f5f99'}
LEVEL_MARK = {'fatal': '\u2716', 'error': '\u2716', 'warning': '\u26a0', 'info': '\u2139'}
TABS = ['Summary', 'Simulation', 'Heating & injection', 'Magnetic field', 'Emission & RF', 'Output',
        'Geometry', 'Potentials']
PALETTE = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
           '#bcbd22', '#17becf']


def lighten(hexcol, f=0.45):
    r, g, b = [int(hexcol[i:i + 2], 16) for i in (1, 3, 5)]
    return '#%02x%02x%02x' % tuple(int(c * (1 - f) + 255 * f) for c in (r, g, b))


# =========================================================================
class ScrollFrame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0)
        sb = ttk.Scrollbar(self, orient='vertical', command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        self.inner.bind('<Configure>', lambda e: self.canvas.configure(scrollregion=self.canvas.bbox('all')))
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor='nw')
        self.canvas.bind('<Configure>', lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        self.canvas.configure(yscrollcommand=sb.set)
        self.canvas.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')
        for w in (self.canvas, self.inner):
            w.bind('<Enter>', self._bind_wheel)
            w.bind('<Leave>', self._unbind_wheel)

    def _bind_wheel(self, _):
        self.canvas.bind_all('<MouseWheel>', self._wheel)
        self.canvas.bind_all('<Button-4>', lambda e: self.canvas.yview_scroll(-2, 'units'))
        self.canvas.bind_all('<Button-5>', lambda e: self.canvas.yview_scroll(2, 'units'))

    def _unbind_wheel(self, _):
        self.canvas.unbind_all('<MouseWheel>')
        self.canvas.unbind_all('<Button-4>')
        self.canvas.unbind_all('<Button-5>')

    def _wheel(self, e):
        self.canvas.yview_scroll(int(-e.delta / 40) or (-1 if e.delta > 0 else 1), 'units')


# =========================================================================
class App:
    def __init__(self, root, folder=None):
        self.root = root
        self.case = None
        self.widgets = {}        # key -> widget (for "jump to field")
        self.marks = {}          # key -> mark label
        self.pulls = []          # static widgets refresh callbacks
        self.dyn = []            # (key list, pull) of dynamic widgets, rebuilt on demand
        self.dyn_keys = set()
        self._loading = False
        self._job = None
        self._pjob = None
        self.issues = []
        self.cur_seg = None
        self.preview_stale = True
        self.derived_var = tk.StringVar()
        self.status_var = tk.StringVar()
        root.title('LePIC+ input editor')
        root.geometry('1240x860')
        root.minsize(980, 640)
        self._style()
        self._toolbar()
        pw = ttk.PanedWindow(root, orient='vertical')
        pw.pack(fill='both', expand=True, padx=6, pady=(0, 4))
        self.nb = ttk.Notebook(pw)
        pw.add(self.nb, weight=5)
        self._help_panel(pw)
        ttk.Label(root, textvariable=self.status_var, anchor='w', relief='sunken').pack(fill='x', side='bottom')
        self.tabs = {}
        for t in TABS:
            sf = ScrollFrame(self.nb) if t not in ('Summary', 'Geometry') else ttk.Frame(self.nb)
            self.tabs[t] = sf
            self.nb.add(sf, text=' %s ' % t)
        self.nb.bind('<<NotebookTabChanged>>', self._tab_changed)
        def place_sash(tries=[0]):
            h = pw.winfo_height()
            if h < 300 and tries[0] < 40:                 # window not laid out yet: try again shortly
                tries[0] += 1
                self.root.after(100, place_sash)
                return
            if h >= 300:
                pw.sashpos(0, int(h * 0.74))
        self.root.after(100, place_sash)
        self._build_static()
        root.protocol('WM_DELETE_WINDOW', self.quit)
        if folder:
            self.load(folder)

    # ------------------------------------------------------------ styling
    def _style(self):
        st = ttk.Style()
        try:
            st.theme_use('clam')
        except tk.TclError:
            pass
        base = ('TkDefaultFont', 10)
        st.configure('.', font=base)
        st.configure('TLabelframe.Label', font=('TkDefaultFont', 10, 'bold'))
        st.configure('Err.TEntry', fieldbackground='#ffd9d9')
        st.configure('Warn.TEntry', fieldbackground='#fff0c2')
        st.configure('Hint.TLabel', foreground='#555555')
        st.configure('Big.TLabel', font=('TkDefaultFont', 11, 'bold'))
        st.configure('Treeview', rowheight=22)

    def _toolbar(self):
        tb = ttk.Frame(self.root)
        tb.pack(fill='x', padx=6, pady=6)
        ttk.Button(tb, text='Open folder\u2026', command=self.open_folder).pack(side='left')
        ttk.Label(tb, text='  Case:').pack(side='left')
        self.case_var = tk.StringVar()
        self.case_box = ttk.Combobox(tb, textvariable=self.case_var, state='readonly', width=22)
        self.case_box.pack(side='left')
        self.case_box.bind('<<ComboboxSelected>>', self._case_selected)
        ttk.Separator(tb, orient='vertical').pack(side='left', fill='y', padx=8)
        self.btn_save = ttk.Button(tb, text='Save', command=self.save)
        self.btn_save.pack(side='left')
        ttk.Button(tb, text='Save as\u2026', command=self.save_as).pack(side='left', padx=3)
        ttk.Button(tb, text='Review changes\u2026', command=self.show_diff).pack(side='left')
        ttk.Button(tb, text='Reload from disk', command=self.reload).pack(side='left', padx=3)
        ttk.Separator(tb, orient='vertical').pack(side='left', fill='y', padx=8)
        ttk.Button(tb, text='Copy to input_dir/ (run this case)', command=self.activate).pack(side='left')
        self.path_var = tk.StringVar()

    def _help_panel(self, pw):
        f = ttk.LabelFrame(pw, text='What does this do?  (click or hover over a field)')
        pw.add(f, weight=1)
        self.help = tk.Text(f, height=5, wrap='word', borderwidth=0, padx=8, pady=4, state='disabled',
                            background='#fbfbf5', font=('TkDefaultFont', 10))
        sb = ttk.Scrollbar(f, command=self.help.yview)
        self.help.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        self.help.pack(fill='both', expand=True)
        self.help.tag_configure('title', font=('TkDefaultFont', 11, 'bold'))
        self.help.tag_configure('code', foreground='#666666', font='TkFixedFont')
        for lv, col in LEVEL_COLOR.items():
            self.help.tag_configure(lv, foreground=col, font=('TkDefaultFont', 10, 'bold'))

    # ------------------------------------------------------------ helpers
    def show_help(self, key):
        hk = help_key(key)
        h = self.help
        h.configure(state='normal')
        h.delete('1.0', 'end')
        if key is None:
            h.insert('end', 'LePIC+ input editor\n', 'title')
            h.insert('end', 'Click or hover over any field to see what it does in the code, its unit, and any '
                     'problem detected. Fields marked \u2716 / \u26a0 have a problem; the Summary tab lists '
                     'all checks. Nothing is written until you press Save (the previous files are kept in '
                     'gui_backups/).')
            h.configure(state='disabled')
            return
        if hk in HELP:
            title, code, text = HELP[hk]
            h.insert('end', title, 'title')
            h.insert('end', '    code: %s\n' % code, 'code')
        for i in self.issues:
            if i.key == key or (self.cur_seg is not None and key and key.startswith('g.f.') and
                                i.key and i.key.startswith('g.seg%d.' % self.cur_seg) and
                                i.key.split('.')[-1] == key.split('.')[-1]):
                h.insert('end', '%s %s\n' % (LEVEL_MARK[i.level], i.msg), i.level)
        if hk in HELP:
            h.insert('end', HELP[hk][2])
            if hk == 'nz' and self.case is not None:
                h.insert('end', '\nMPI processes in %s: %s' % (M.RUN_SCRIPT, M.mpi_summary(self.case.mpi)))
        h.configure(state='disabled')

    def _reg(self, key, widget, mark, pull, dynamic=False):
        self.widgets[key] = widget
        if mark is not None:
            self.marks[key] = mark
        (self.dyn if dynamic else self.pulls).append((key, pull))
        if dynamic:
            self.dyn_keys.add(key)
        for ev in ('<FocusIn>', '<Enter>'):
            widget.bind(ev, lambda e, k=key: self.show_help(k), add='+')

    def _row(self, f):
        r = getattr(f, '_r', 0)
        f._r = r + 1
        return r

    def _label(self, f, r, text, c=0):
        ttk.Label(f, text=text).grid(row=r, column=c, sticky='w', padx=(8, 8), pady=2)

    def _mark(self, f, r, c):
        m = ttk.Label(f, text='', width=2)
        m.grid(row=r, column=c, sticky='w')
        return m

    def section(self, parent, title, desc=None):
        lf = ttk.LabelFrame(parent, text=title)
        lf.pack(fill='x', padx=10, pady=6)
        lf.columnconfigure(0, minsize=300)
        lf.columnconfigure(4, weight=1)
        lf._r = 0
        if desc:
            ttk.Label(lf, text=desc, style='Hint.TLabel', wraplength=900, justify='left').grid(
                row=0, column=0, columnspan=5, sticky='w', padx=8, pady=(2, 4))
            lf._r = 1
        return lf

    # --- generic field builders -------------------------------------------------
    def _accessors(self, key):
        return (lambda: self.case.cond.tok(key)), (lambda v: self.case.cond.set_tok(key, v))

    def entry(self, f, label, key, unit='', width=12, get=None, put=None, dynamic=False, row=None,
              values=None):
        if get is None:
            get, put = self._accessors(key)
        r = self._row(f) if row is None else row
        self._label(f, r, label)
        var = tk.StringVar()
        if values is not None:
            e = ttk.Combobox(f, textvariable=var, width=width, values=values)
        else:
            e = ttk.Entry(f, textvariable=var, width=width)
        e.grid(row=r, column=1, sticky='w')
        ttk.Label(f, text=unit, style='Hint.TLabel').grid(row=r, column=2, sticky='w', padx=(4, 8))
        mk = self._mark(f, r, 3)

        def pull():
            self._loading = True
            var.set(get())
            self._loading = False

        def push(*_):
            if not self._loading:
                put(var.get().strip())
                self.changed()
        var.trace_add('write', push)
        self._reg(key, e, mk, pull, dynamic)
        return e

    def signed(self, f, check_text, label, key, unit='', width=12, on_when='positive', lock=True,
               get=None, put=None, dynamic=False, stack=False):
        """Check-box + magnitude entry for a number whose SIGN encodes a switch."""
        if get is None:
            get, put = self._accessors(key)
        r = self._row(f)
        self._label(f, r, label)
        vc, ve = tk.BooleanVar(), tk.StringVar()
        if stack:
            e = ttk.Entry(f, textvariable=ve, width=width)
            e.grid(row=r, column=1, sticky='w')
            ttk.Label(f, text=unit, style='Hint.TLabel').grid(row=r, column=2, sticky='w', padx=4)
            cb = ttk.Checkbutton(f, text=check_text, variable=vc)
            cb.grid(row=self._row(f), column=1, columnspan=3, sticky='w')
        else:
            box = ttk.Frame(f)
            box.grid(row=r, column=1, columnspan=2, sticky='w')
            cb = ttk.Checkbutton(box, text=check_text, variable=vc)
            cb.pack(side='left')
            e = ttk.Entry(box, textvariable=ve, width=width)
            e.pack(side='left', padx=(8, 2))
            ttk.Label(box, text=unit, style='Hint.TLabel').pack(side='left', padx=4)
        mk = self._mark(f, r, 3)

        def pull():
            self._loading = True
            t = get()
            neg = M.is_neg(t)
            vc.set((not neg) if on_when == 'positive' else neg)
            ve.set(M.strip_sign(t))
            self._loading = False
            e.configure(state='normal' if (vc.get() or not lock or on_when != 'positive') else 'disabled')

        def push(*_):
            if self._loading:
                return
            on = vc.get()
            neg = (not on) if on_when == 'positive' else on
            put(M.with_sign(ve.get().strip(), neg))
            e.configure(state='normal' if (on or not lock or on_when != 'positive') else 'disabled')
            self.changed()
        vc.trace_add('write', push)
        ve.trace_add('write', push)
        self._reg(key, e, mk, pull, dynamic)
        cb.bind('<FocusIn>', lambda ev, k=key: self.show_help(k), add='+')
        cb.bind('<Enter>', lambda ev, k=key: self.show_help(k), add='+')
        return e

    def yesno(self, f, text, key, get=None, put=None, dynamic=False, after=None):
        if get is None:
            get, put = self._accessors(key)
        r = self._row(f)
        var = tk.BooleanVar()
        cb = ttk.Checkbutton(f, text=text, variable=var)
        cb.grid(row=r, column=0, columnspan=3, sticky='w', padx=8, pady=2)
        mk = self._mark(f, r, 3)

        def pull():
            self._loading = True
            var.set(M._yn(get()))
            self._loading = False
            if after:
                after()

        def push(*_):
            if self._loading:
                return
            if var.get() != M._yn(get()):
                put('yes' if var.get() else 'no')
            if after:
                after()
            self.changed()
        var.trace_add('write', push)
        self._reg(key, cb, mk, pull, dynamic)
        return cb

    def choice(self, f, label, key, options, get=None, put=None, width=30, dynamic=False, after=None, match='first'):
        """options: list of (token to write, label shown)."""
        if get is None:
            get, put = self._accessors(key)
        r = self._row(f)
        self._label(f, r, label)
        var = tk.StringVar()
        cbx = ttk.Combobox(f, textvariable=var, values=[o[1] for o in options], state='readonly', width=width)
        cbx.grid(row=r, column=1, columnspan=2, sticky='w')
        mk = self._mark(f, r, 3)

        def find(tok):
            t = (tok or '').lower()
            for w, lab in options:
                if match == 'first' and t[:1] == w.lower()[:1]:
                    return lab
                if match == 'exact' and t == w.lower():
                    return lab
            return '(unrecognised: %s)' % tok

        def pull():
            self._loading = True
            var.set(find(get()))
            self._loading = False
            if after:
                after()

        def push(*_):
            if self._loading:
                return
            for w, lab in options:
                if lab == var.get():
                    cur = find(get())
                    if cur != lab:
                        put(w)
            if after:
                after()
            self.changed()
        var.trace_add('write', push)
        self._reg(key, cbx, mk, pull, dynamic)
        return cbx

    # ------------------------------------------------------------ change cycle
    def changed(self):
        if self._loading or self.case is None:
            return
        if self._job:
            self.root.after_cancel(self._job)
        self._job = self.root.after(200, self.recheck)
        self.preview_stale = True

    def refresh_all(self):
        for _, p in self.pulls:
            p()
        for _, p in self.dyn:
            p()
        self.recheck()

    def recheck(self):
        self._job = None
        if self.case is None:
            return
        c = self.case
        # structural sync of dynamic widgets
        self._sync_dynamic()
        self.issues = M.validate(c)
        worst = {}
        for i in self.issues:
            if i.key is None:
                continue
            ks = [i.key]
            if self.cur_seg is not None and i.key.startswith('g.seg%d.' % self.cur_seg):
                ks.append('g.f.' + i.key.split('.')[-1])
            for k in ks:
                if k not in worst or M.Issue.ORDER[i.level] < M.Issue.ORDER[worst[k]]:
                    worst[k] = i.level
        for k, mk in self.marks.items():
            try:
                lv = worst.get(k)
                mk.configure(text=LEVEL_MARK[lv] if lv else '', foreground=LEVEL_COLOR[lv] if lv else 'black')
                w = self.widgets.get(k)
                if type(w) is ttk.Entry:
                    w.configure(style='Err.TEntry' if lv in ('fatal', 'error') else
                                ('Warn.TEntry' if lv == 'warning' else 'TEntry'))
            except tk.TclError:
                pass
        self._update_summary()
        self._update_title()
        self._refresh_potential_usage()
        self._refresh_seg_tree()
        try:
            cur = self.nb.tab(self.nb.select(), 'text').strip()
        except tk.TclError:
            cur = ''
        if cur == 'Geometry' and self.preview_stale:
            self.draw_preview()

    @staticmethod
    def _int(v, d=0):
        try:
            return int(v)
        except Exception:
            return d

    def _update_title(self):
        mod = self.case.modified()
        self.root.title('LePIC+ input editor - %s%s' % (self.case.folder, '  *' if mod else ''))
        self.path_var.set(self.case.folder)
        nf = sum(1 for i in self.issues if i.level == 'fatal')
        ne = sum(1 for i in self.issues if i.level == 'error')
        nw = sum(1 for i in self.issues if i.level == 'warning')
        self.status_var.set('  %s    |    %d blocking, %d problem(s), %d warning(s)    |    %s' % (
            'UNSAVED CHANGES' if mod else 'no unsaved changes', nf, ne, nw, self.case.folder))

    # ------------------------------------------------------------ jump
    def jump(self, key, tab=None):
        if tab and tab in self.tabs:
            self.nb.select(self.tabs[tab])
        if key is None:
            return
        m = re.match(r'^g\.seg(\d+)\.(\d+)$', key)
        if m:
            self.nb.select(self.tabs['Geometry'])
            self.select_segment(int(m.group(1)))
            key = 'g.f.' + m.group(2)
        w = self.widgets.get(key)
        if w is not None:
            try:
                w.focus_set()
            except tk.TclError:
                pass
            self.show_help(key)

    def _tab_changed(self, _=None):
        try:
            name = self.nb.tab(self.nb.select(), 'text').strip()
        except tk.TclError:
            return
        if name == 'Geometry' and self.case is not None and self.preview_stale:
            self.root.after(50, self.draw_preview)

    # ------------------------------------------------------------ load / save
    def find_cases(self, folder):
        out = [('(this folder)', folder)]
        try:
            for d in sorted(os.listdir(folder)):
                p = os.path.join(folder, d)
                if os.path.isdir(p) and os.path.isfile(os.path.join(p, 'conditions.inp')):
                    out.append((d, p))
        except OSError:
            pass
        return out

    def load(self, folder):
        try:
            case = M.Case(folder)
        except Exception as e:
            messagebox.showerror('Cannot open', str(e))
            return False
        self.case = case
        self.root_folder = os.path.abspath(folder)
        base = os.path.abspath(folder)
        parent = os.path.dirname(base)
        # case list: siblings if this is a sub-case, children otherwise
        if os.path.basename(parent) == 'input_dir' or os.path.isfile(os.path.join(parent, 'conditions.inp')):
            base = parent
        self.cases = self.find_cases(base)
        self.case_box.configure(values=[c[0] if c[1] != base else '(input_dir itself)' for c in self.cases])
        cur = [i for i, c in enumerate(self.cases) if c[1] == os.path.abspath(folder)]
        self.case_box.current(cur[0] if cur else 0)
        self.path_var.set(case.folder)
        self.cur_seg = 0 if case.geom.segs else None
        self.slice_pos = {'xy': None, 'xz': None, 'yz': None}
        self._seg_cache = None
        self._rebuild_dynamic()
        self.refresh_all()
        self.select_segment(self.cur_seg)
        self.preview_stale = True
        self.show_help(None)
        if self.nb.tab(self.nb.select(), 'text').strip() == 'Geometry':
            self.draw_preview()
        return True

    def _confirm_discard(self):
        if self.case is not None and self.case.modified():
            r = messagebox.askyesnocancel('Unsaved changes', 'Save your changes first?')
            if r is None:
                return False
            if r:
                return self.save()
        return True

    def open_folder(self):
        if not self._confirm_discard():
            return
        d = filedialog.askdirectory(title='Choose a folder containing conditions.inp, geometry.inp, boundary.inp',
                                    initialdir=self.case.folder if self.case else os.getcwd())
        if d:
            self.load(d)

    def _case_selected(self, _):
        i = self.case_box.current()
        if i < 0:
            return
        if not self._confirm_discard():
            return
        self.load(self.cases[i][1])

    def reload(self):
        if self.case is None:
            return
        if self.case.modified() and not messagebox.askyesno('Reload', 'Discard all unsaved changes?'):
            return
        self.load(self.case.folder)

    def _blocking_check(self):
        self.recheck()
        fatal = [i for i in self.issues if i.level == 'fatal']
        if fatal:
            messagebox.showerror('Cannot save', 'Please fix these first:\n\n' +
                                 '\n'.join('- ' + i.msg for i in fatal[:10]))
            return False
        errs = [i for i in self.issues if i.level == 'error']
        if errs:
            return messagebox.askyesno(
                'Problems found',
                'The simulation would stop or misbehave with these settings:\n\n' +
                '\n'.join('- ' + i.msg for i in errs[:8]) + ('\n...' if len(errs) > 8 else '') +
                '\n\nSave anyway?')
        return True

    def save(self, folder=None):
        if self.case is None:
            return False
        if not self._blocking_check():
            return False
        try:
            written = self.case.save(folder)
        except Exception as e:
            messagebox.showerror('Save failed', str(e))
            return False
        target = folder or self.case.folder
        self.load(target)
        self.status_var.set('  Saved %s  (previous versions are in gui_backups/)' % (
            ', '.join(os.path.basename(w) for w in written) or 'nothing to write - files already up to date'))
        return True

    def activate(self):
        if self.case is None:
            return
        dst = self.case.input_dir()
        if dst is None:
            messagebox.showinfo('Copy to input_dir/', 'This case already is input_dir/ (or is not inside it): '
                                'LePIC+ reads ./input_dir/conditions.inp, geometry.inp, boundary.inp directly.')
            return
        if self.case.modified():
            if not messagebox.askyesno('Unsaved changes', 'Save the changes to this case first?'):
                return
            if not self.save():
                return
        if not messagebox.askyesno('Copy to input_dir/',
                                   'LePIC+ only reads the files directly in input_dir/. Copy this case there?\n\n'
                                   'The current input_dir/ files are kept in input_dir/gui_backups/.'):
            return
        try:
            names = self.case.activate()
        except Exception as e:
            messagebox.showerror('Copy failed', str(e))
            return
        self.status_var.set('  Copied %s to %s  (previous files in gui_backups/)' % (', '.join(names), dst))

    def save_as(self):
        if self.case is None:
            return
        d = filedialog.askdirectory(title='Save the three files into this folder (it can be a new folder)',
                                    initialdir=os.path.dirname(self.case.folder), mustexist=False)
        if d:
            self.save(d)

    def show_diff(self):
        if self.case is None:
            return
        lines = self.case.diff()
        w = tk.Toplevel(self.root)
        w.title('Changes that will be written')
        w.geometry('900x500')
        t = tk.Text(w, wrap='none', font='TkFixedFont')
        sb = ttk.Scrollbar(w, command=t.yview)
        t.configure(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        t.pack(fill='both', expand=True)
        t.tag_configure('add', foreground='#0a7a2f', background='#e8f7ec')
        t.tag_configure('del', foreground='#a01020', background='#fbeaec')
        t.tag_configure('hdr', foreground='#555555', font='TkFixedFont')
        if not lines:
            t.insert('end', 'No changes: the files on disk are identical to what the editor would write.')
        for ln in lines:
            tag = 'hdr' if ln.startswith(('---', '+++', '@@')) else ('add' if ln.startswith('+') else
                                                                       ('del' if ln.startswith('-') else None))
            t.insert('end', ln + '\n', tag)
        t.configure(state='disabled')

    def quit(self):
        if self._confirm_discard():
            self.root.destroy()


    # =====================================================================
    #  PAGES
    # =====================================================================
    def _build_static(self):
        self._page_summary()
        self._page_simulation()
        self._page_heating()
        self._page_bfield()
        self._page_emission()
        self._page_output()
        self._page_geometry()
        self._page_potentials()

    def _drop(self, pattern):
        rx = re.compile(pattern)
        self.dyn = [(k, p) for (k, p) in self.dyn if not rx.match(k)]
        for d in (self.widgets, self.marks):
            for k in [k for k in d if rx.match(k)]:
                del d[k]

    def _set_rows(self, frame, rows, visible):
        for r in rows:
            for w in frame.grid_slaves(row=r):
                (w.grid if visible else w.grid_remove)()

    # ------------------------------------------------------------ summary
    def _page_summary(self):
        f = self.tabs['Summary']
        f.columnconfigure(0, weight=1)
        f.rowconfigure(1, weight=1)
        top = ttk.LabelFrame(f, text='Derived quantities (what the code will print at start-up)')
        top.grid(row=0, column=0, sticky='ew', padx=10, pady=6)
        ttk.Label(top, textvariable=self.derived_var, justify='left', font='TkFixedFont').pack(
            anchor='w', padx=8, pady=4)
        box = ttk.LabelFrame(f, text='Checks  -  select a line to read it in full, double-click to go to the field')
        box.grid(row=1, column=0, sticky='nsew', padx=10, pady=6)
        box.columnconfigure(0, weight=1)
        box.rowconfigure(0, weight=1)
        tv = ttk.Treeview(box, columns=('lvl', 'where', 'msg'), show='headings', selectmode='browse')
        for c, t, w, st in (('lvl', 'Level', 80, False), ('where', 'Where', 140, False), ('msg', 'Message', 800, True)):
            tv.heading(c, text=t)
            tv.column(c, width=w, stretch=st)
        sb = ttk.Scrollbar(box, command=tv.yview)
        tv.configure(yscrollcommand=sb.set)
        tv.grid(row=0, column=0, sticky='nsew')
        sb.grid(row=0, column=1, sticky='ns')
        for lv, col in LEVEL_COLOR.items():
            tv.tag_configure(lv, foreground=col)
        tv.bind('<Double-1>', self._issue_jump)
        tv.bind('<<TreeviewSelect>>', self._issue_select)
        self.issue_tv = tv
        self.issue_by_iid = {}
        ttk.Label(f, text='"Blocking" problems prevent saving; "problem" means the Fortran code would stop or '
                  'misbehave (you may still save); warnings and notes are informational.',
                  style='Hint.TLabel', wraplength=1000).grid(row=2, column=0, sticky='w', padx=12, pady=(0, 6))

    def _issue_select(self, _):
        sel = self.issue_tv.selection()
        if not sel:
            return
        i = self.issue_by_iid.get(sel[0])
        if i:
            h = self.help
            h.configure(state='normal')
            h.delete('1.0', 'end')
            h.insert('end', '%s %s' % (LEVEL_MARK[i.level], i.msg), i.level)
            h.configure(state='disabled')

    def _issue_jump(self, _):
        sel = self.issue_tv.selection()
        if sel and sel[0] in self.issue_by_iid:
            i = self.issue_by_iid[sel[0]]
            self.jump(i.key, i.tab)

    def _update_summary(self):
        tv = self.issue_tv
        tv.delete(*tv.get_children())
        self.issue_by_iid = {}
        names = {'fatal': 'BLOCKING', 'error': 'problem', 'warning': 'warning', 'info': 'note'}
        for i in self.issues:
            iid = tv.insert('', 'end', values=(names[i.level], i.tab or '', i.msg), tags=(i.level,))
            self.issue_by_iid[iid] = i
        if not self.issues:
            tv.insert('', 'end', values=('OK', '', 'No problems found.'))
        p, _ = M.decode(self.case)
        box, segs, _ = M.decode_geometry(self.case)
        d = M.derived(p, box) if None not in box else None
        if not d:
            self.derived_var.set('(fill in valid numbers to see derived quantities)')
            return
        n3 = (p['nx'], p['ny'], p['nz'])
        L = [abs(b) for b in box]
        mgs = '%s (from the file)' % p['ng']
        lines = [
            'Grid            %d x %d x %d = %.2f M cells      box %.3g x %.3g x %.3g cm' % (
                n3[0], n3[1], n3[2], d['ncells'] / 1e6, L[0], L[1], L[2]),
            'Cell size       dx, dy, dz = %.1f, %.1f, %.1f um      multigrid levels used: %s' % (
                d['h'][0] * 1e6, d['h'][1] * 1e6, d['h'][2] * 1e6, mgs),
            'Electrons       vte = %.3g m/s      Debye length = %.3g um      dx/lambda_D = %.2f' % (
                d['vt'], d['lde'] * 1e6, d['h'][0] / d['lde']),
            'Time step       dt = %.3g s      wpe*dt = %.3g      (wpe = %.3g rad/s)' % (d['dt'], d['wpdt'], d['wp']),
            'Run length      %s steps x dt = %s      (t_max would need %s steps; not used)      heating every %d steps%s' % (
                p['nsteps'] if p['nsteps'] else '?', ('%.3g s' % d['t_end']) if 't_end' in d else '?',
                ('%.3g' % d['nsteps_tmax']) if 'nsteps_tmax' in d else '?', d['ns_heat'],
                ('      nu_h*dt_heat = %.3g (must be <= 1)' % d['nudt']) if 'nudt' in d else ''),
        ]
        if p['antenna']:
            lines.append('RF antenna      %s      skin depth c/wpe(n0) = %.3g cm (true eps0)      %s' % (
                p['rf_mode'], d['skin'] * 100,
                ('update every 1/(f0 dt) = %.0f steps' % d['ns_RF']) if 'ns_RF' in d else 'f0 must be > 0'))
        lines.append('MPI processes   %s' % M.mpi_summary(self.case.mpi) + ('   (nz = %d)' % n3[2]))
        self.derived_var.set('\n'.join(lines))
        self.short_var.set('dt = %.3g s    wpe*dt = %.3g    Debye length = %.3g um    dx/lambda_D = %.2f' % (
            d['dt'], d['wpdt'], d['lde'] * 1e6, d['h'][0] / d['lde']))

    # ------------------------------------------------------------ simulation
    def _page_simulation(self):
        f = self.tabs['Simulation'].inner
        self.short_var = tk.StringVar()
        ttk.Label(f, textvariable=self.short_var, style='Big.TLabel').pack(anchor='w', padx=14, pady=(8, 0))
        s = self.section(f, 'Plasma and physics')
        self.rname_box = self.entry(s, 'Reactions file', 'rname', '(in input_dir/)', width=24, values=[])
        self.entry(s, 'Electron temperature  Te', 'Te', 'eV')
        self.entry(s, 'Plasma density  n0', 'n0', 'm\u207b\u00b3')
        self.entry(s, 'Neutral gas density', 'ngas', 'm\u207b\u00b3')
        self.signed(s, 'Start without initial plasma', 'Macro-particles per cell', 'np_cell', '',
                    on_when='negative', width=8)
        self.signed(s, 'Iterative mode (auto-adjust power / current / voltage)', 'Permittivity scaling  k_eps0',
                    'k_eps0', '', on_when='negative', width=10)
        self.yesno(s, 'Coulomb collisions', 'coulomb')
        self.choice(s, 'Push scheme', 'push', [('momentum', 'momentum - standard (momentum-conserving)'),
                                              ('energy', 'energy - EC-PIC (energy-conserving)')],
                    width=40, match='exact')
        s = self.section(f, 'Time')
        self.entry(s, 'Number of time steps', 'nsteps', '(run length; line 1)', width=12)
        self.entry(s, 'Simulated time  t_max', 'tmax', 's   (read but not used by LePIC+)')
        self.entry(s, 'Time-step coefficient  kt', 'kt', '(dt = kt * dx_min / vte)')
        s = self.section(f, 'Grid and Poisson solver',
                         'The box size itself is set in the Geometry tab. The number of cells below is NOT an exponent.')
        self.entry(s, 'Cells along x', 'nx', '', width=8)
        self.entry(s, 'Cells along y', 'ny', '', width=8)
        self.entry(s, 'Cells along z', 'nz', '', width=8)
        self.entry(s, 'Multigrid levels', 'ng', '', width=8)
        self.entry(s, 'Solver precision  eps', 'eps', '')
        self.entry(s, 'SOR relaxation  omega', 'omega', '(1 \u2264 \u03c9 < 2)')

    # ------------------------------------------------------------ heating / injection
    def _page_heating(self):
        f = self.tabs['Heating & injection'].inner
        s = self.section(f, 'Heating / injection region  (cm)',
                         'Electron heating uses the x-range; particle injection uses the whole box.')
        self.entry(s, 'x from', 'xl_pow', 'cm')
        self.entry(s, 'x to', 'xr_pow', 'cm')
        self.entry(s, 'y from', 'yl_pow', 'cm')
        self.signed(s, 'Annular: heat only OUTSIDE the radius (round geometry)', 'y to', 'yr_pow', 'cm',
                    on_when='negative', width=10)
        self.entry(s, 'z from', 'zl_pow', 'cm')
        self.entry(s, 'z to', 'zr_pow', 'cm')
        s = self.section(f, 'Initial plasma')
        self.entry(s, 'Load plasma for x \u2264', 'x_load', 'cm')
        s = self.section(f, 'Electron heating')
        self.signed(s, 'Control temperature through absorbed power', 'Absorbed power  Pabs', 'Pabs', 'W')
        self.signed(s, 'Heating on', 'Heating frequency  nu_h', 'nu_h', 's\u207b\u00b9')
        self.signed(s, 'Self-consistent (delta = c/wpe, updated every 1/f0)', 'Skin depth  delta', 'delta_heat',
                    'cm   (0 = uniform heating)', on_when='negative', width=10)
        s = self.section(f, 'Particle injection')
        self.signed(s, 'Injection on', 'Injected current  I_inj', 'I_inj', 'A')
        opts = [('1', '1 - fixed current, uniform in the box'), ('2', '2 - replace ions lost at the walls'),
                ('3', '3 - electron beam along +z'), ('4', '4 - emissive cathode along z')]

        def og():
            return M.strip_sign(self.case.cond.tok('opt_inj'))

        def op(v):
            self.case.cond.set_tok('opt_inj', M.with_sign(v, M.is_neg(self.case.cond.tok('opt_inj'))))
        self.choice(s, 'Mode', 'opt_inj', opts, get=og, put=op, width=40, match='exact')

        def cg():
            return self.case.cond.tok('opt_inj')

        def cp(v):
            self.case.cond.set_tok('opt_inj', v)
        r = self._row(s)
        var = tk.BooleanVar()
        cb = ttk.Checkbutton(s, text='Cosine profile instead of uniform', variable=var)
        cb.grid(row=r, column=1, columnspan=3, sticky='w')

        def pull():
            self._loading = True
            var.set(M.is_neg(cg()))
            self._loading = False

        def push(*_):
            if not self._loading:
                cp(M.with_sign(M.strip_sign(cg()), var.get()))
                self.changed()
        var.trace_add('write', push)
        self.pulls.append(('opt_inj_cos', pull))
        cb.bind('<FocusIn>', lambda e: self.show_help('opt_inj'), add='+')

    # ------------------------------------------------------------ magnetic field
    def _page_bfield(self):
        f = self.tabs['Magnetic field'].inner
        s = self.section(f, 'Magnetic field', 'The fields of all entries are added together.')
        r = self._row(s)
        self._label(s, r, 'Number of entries')
        self.nb_var = tk.StringVar()
        sp = ttk.Spinbox(s, from_=0, to=M.MAX_B, width=4, textvariable=self.nb_var)
        sp.grid(row=r, column=1, sticky='w')
        mk = self._mark(s, r, 3)

        def pull():
            self._loading = True
            t = self.case.cond.tok('nB')
            self.nb_var.set(str(abs(int(t))) if M.is_int(t) else t)
            self._loading = False

        def push(*_):
            if self._loading:
                return
            try:
                n = int(self.nb_var.get())
            except ValueError:
                return
            self.case.cond.set_nb(n)
            self.changed()
        self.nb_var.trace_add('write', push)
        self._reg('nB', sp, mk, pull)
        r = self._row(s)
        var = tk.BooleanVar()
        cb = ttk.Checkbutton(s, text='Positive ions are not magnetized', variable=var)
        cb.grid(row=r, column=0, columnspan=3, sticky='w', padx=8)

        def pull2():
            self._loading = True
            var.set(M.is_neg(self.case.cond.tok('nB')))
            self._loading = False

        def push2(*_):
            if not self._loading:
                t = self.case.cond.tok('nB')
                self.case.cond.set_tok('nB', M.with_sign(M.strip_sign(t), var.get()))
                self.changed()
        var.trace_add('write', push2)
        self.pulls.append(('nB_neg', pull2))
        cb.bind('<FocusIn>', lambda e: self.show_help('nB'), add='+')
        self.entry(s, 'Map grid  nx', 'nBx', 'cells', width=8)
        self.entry(s, 'Map grid  ny', 'nBy', 'cells', width=8)
        self.entry(s, 'Map grid  nz', 'nBz', 'cells', width=8)
        self.b_holder = ttk.Frame(f)
        self.b_holder.pack(fill='x')

    B_TYPES = [('s', 'Not used (skip this entry)'), ('n', 'Analytical profile'), ('y', 'Field map from files')]
    B_MODELS = [('gauss', 'Gaussian profile along x'), ('PG', 'PG-current filter field (ELISE fit)'),
                ('EEmagnets', 'Cusp (extraction) magnets'), ('thr', 'Hall-thruster magnet circuit')]

    def _build_bentries(self):
        for w in self.b_holder.winfo_children():
            w.destroy()
        self._drop(r'^B\d+\.')
        self._brows = {}
        for k, rec in enumerate(self.case.cond.bfiles):
            lf = ttk.LabelFrame(self.b_holder, text='Entry %d' % (k + 1))
            lf.pack(fill='x', padx=10, pady=6)
            lf._r = 0
            rows = {}

            def acc(i, rec=rec):
                return (lambda: rec.get(i)), (lambda v: rec.set(i, v))
            g, p = acc(0)
            self.choice(lf, 'Type', 'B%d.file' % k, self.B_TYPES, get=g, put=p, width=34, dynamic=True,
                        after=lambda k=k: self._b_visibility(k))
            rows['file'] = lf._r - 1
            g, p = acc(1)
            self.choice(lf, 'Model', 'B%d.name' % k, self.B_MODELS, get=g, put=p, width=34, dynamic=True,
                        after=lambda k=k: self._b_visibility(k))
            rows['model'] = lf._r - 1
            self.entry(lf, 'File name prefix', 'B%d.prefix' % k, 'reads <prefix>Bx/By/Bz.dat in mag_field_files/',
                       width=18, get=g, put=p, dynamic=True)
            rows['prefix'] = lf._r - 1
            g, p = acc(2)
            self.entry(lf, 'Scale factor', 'B%d.scale' % k, '', get=g, put=p, dynamic=True)
            rows['scale'] = lf._r - 1
            g, p = acc(3)
            self.entry(lf, 'Amplitude  B0', 'B%d.B0' % k, 'G', get=g, put=p, dynamic=True)
            rows['B0'] = lf._r - 1
            g, p = acc(4)
            self.choice(lf, 'Direction', 'B%d.dir' % k, [('Bx', 'Bx'), ('By', 'By'), ('Bz', 'Bz')], get=g, put=p,
                        width=8, dynamic=True, match='exact')
            rows['dir'] = lf._r - 1
            for nm, i, unit, lab in (('x0', 6, 'cm', 'Position  x0'), ('z0', 8, 'cm', 'Position  z0'),
                                     ('dL', 5, 'cm', 'Width / length  dL')):
                g, p = acc(i)
                self.entry(lf, lab, 'B%d.%s' % (k, nm), unit, get=g, put=p, dynamic=True)
                rows[nm] = lf._r - 1
            self._brows[k] = (lf, rows)

    def _b_visibility(self, k):
        if k not in getattr(self, '_brows', {}) or k >= len(self.case.cond.bfiles):
            return
        lf, rows = self._brows[k]
        rec = self.case.cond.bfiles[k]
        mode = rec.get(0)[:1].lower()
        ch = rec.get(1)[:1].lower()
        show = {'file'}
        if mode == 'y':
            show |= {'prefix', 'scale'}
        elif mode == 'n':
            show |= {'model'}
            by = {'g': {'B0', 'dir', 'x0', 'dL'}, 'p': {'B0', 'dir'}, 'e': {'B0', 'x0', 'z0', 'dL'},
                  't': {'B0', 'x0', 'dL'}}
            show |= by.get(ch, {'B0', 'dir', 'x0', 'z0', 'dL'})
        for nm, r in rows.items():
            self._set_rows(lf, [r], nm in show)

    # ------------------------------------------------------------ emission & RF
    def _page_emission(self):
        f = self.tabs['Emission & RF'].inner
        s = self.section(f, 'Grid apertures (holes)')
        self.yesno(s, 'Drill apertures in a wall', 'holes')
        self.entry(s, 'Aperture size  Ly', 'Lhy', 'cm', width=8)
        self.entry(s, 'Aperture size  Lz', 'Lhz', 'cm', width=8)
        self.entry(s, 'Number of apertures along y', 'nhy', '', width=8)
        self.entry(s, 'Number of apertures along z', 'nhz', '', width=8)
        self.entry(s, 'Label of the wall / grid', 'ind_g', '', width=8)
        s = self.section(f, 'Surface emission of negative ions')
        self.entry(s, 'Negative-ion flux  jH-', 'jne', 'mA/cm\u00b2   (0 = off)')
        self.entry(s, 'Wall label for impact statistics', 'num_grd', '', width=8)
        s = self.section(f, 'Dielectric walls')
        self.entry(s, 'Surface capacitance  Ca/S', 'Ca', 'F/m\u00b2')
        s = self.section(f, 'Secondary electron emission / emissive cathode')
        self.signed(s, 'Enabled', 'Emission coefficient  gamma', 'gam_sec', 'electrons per ion', width=8)
        self.entry(s, 'Emitting electrode label', 'igrid_sec', '', width=8)
        self.entry(s, 'Emission temperature / energy  TH', 'THm', 'eV   (also used by e-beam and negative ions)')
        s = self.section(f, 'RF')
        self.choice(s, 'Mode', 'rf', [('No', 'No RF'),
                                      ('Electrode', 'Electrode - RF voltage (not implemented in LePIC+)'),
                                      ('Wrapped', 'Wrapped - inductive coil around the driver wall'),
                                      ('Planar', 'Planar - flat coil under a window (model field)'),
                                      ('Coil', 'Coil - real coil turns from coil.inp (E and B)')],
                    width=52, match='exact', get=self._rf_get,
                    put=lambda v: self.case.cond.set_tok('rf', v), after=self._rf_visibility)
        self._rf_rows = {}
        self._rf_units = {}
        for key, lab, unit in (('phi0', 'Voltage amplitude 1', 'V'), ('f0', 'Frequency  f0', 'Hz'),
                               ('phi1', 'Voltage amplitude 2', 'V'), ('f1', 'Frequency 2', 'Hz'),
                               ('E0', 'Initial drive amplitude', 'V/m'),
                               ('R_win', 'Window radius  R_win', 'cm   (0 = y_max - ymax/2)')):
            self.entry(s, lab, key, unit)
            self._rf_rows[key] = s._r - 1
        self._rf_frame = s
        self.coil_sec = self.section(f, 'Coil turns  (coil.inp)',
                                     'Concentric loops coaxial with x. Every turn must be behind the window: '
                                     'x < x_min of the heating region (clamped to >= 0).')
        r = self._row(self.coil_sec)
        self._label(self.coil_sec, r, 'Number of turns')
        self.ncoil_var = tk.StringVar()
        sp = ttk.Spinbox(self.coil_sec, from_=1, to=50, width=4, textvariable=self.ncoil_var)
        sp.grid(row=r, column=1, sticky='w')
        mk = self._mark(self.coil_sec, r, 3)

        def pull():
            self._loading = True
            self.ncoil_var.set(self.case.coil.count.get(0))
            self._loading = False

        def push(*_):
            if self._loading:
                return
            try:
                n = int(self.ncoil_var.get())
            except ValueError:
                return
            self.case.coil.set_count(n)
            self.changed()
        self.ncoil_var.trace_add('write', push)
        self._reg('coil.n', sp, mk, pull)
        self.coil_holder = ttk.Frame(self.coil_sec)
        self.coil_holder.grid(row=self._row(self.coil_sec), column=0, columnspan=5, sticky='ew')

    def _rf_get(self):
        m = M.rf_mode(self.case.cond.tok('rf'))
        return m if m else self.case.cond.tok('rf')

    def _rf_visibility(self):
        if self.case is None or not hasattr(self, '_rf_frame'):
            return
        mode = M.rf_mode(self.case.cond.tok('rf'))
        show = {'Electrode': {'phi0', 'f0', 'phi1', 'f1'}, 'Wrapped': {'f0', 'E0'},
                'Planar': {'f0', 'E0', 'R_win'}, 'Coil': {'f0', 'E0', 'R_win'}}.get(mode, set())
        for k, r in self._rf_rows.items():
            self._set_rows(self._rf_frame, [r], k in show)
        unit = 'A   (coil current, peak)' if mode == 'Coil' else 'V/m'
        for w in self._rf_frame.grid_slaves(row=self._rf_rows['E0'], column=2):
            w.configure(text=unit)
        if mode == 'Coil':
            self.coil_sec.pack(fill='x', padx=10, pady=6, after=self._rf_frame)
        else:
            self.coil_sec.pack_forget()

    def _build_coil(self):
        for w in self.coil_holder.winfo_children():
            w.destroy()
        self._drop(r'^coil\d+\.')
        g = ttk.Frame(self.coil_holder)
        g.pack(fill='x')
        g._r = 0
        for j, rec in enumerate(self.case.coil.turns):
            for i, (nm, unit) in enumerate((('radius a', 'cm'), ('position x', 'cm'), ('relative current w', ''))):
                self.entry(g, 'Turn %d: %s' % (j + 1, nm), 'coil%d.%d' % (j, i), unit, width=10,
                           get=lambda rec=rec, i=i: rec.get(i), put=lambda v, rec=rec, i=i: rec.set(i, v),
                           dynamic=True)

    # ------------------------------------------------------------ output / restart
    def _page_output(self):
        f = self.tabs['Output'].inner
        s = self.section(f, 'Saving data')
        self.signed(s, 'Keep accumulating averages (do not reset after saving)', 'Save every', 'nsav', 'time steps',
                    on_when='negative', width=10)
        self.signed(s, 'Do not time-average the saved data', 'Backup every', 'nbak', 'time steps',
                    on_when='negative', width=10)
        s = self.section(f, 'Sequence of profile snapshots')
        self.signed(s, 'Enabled', 'Number of snapshots', 'cnt_plt', '', width=8)
        self.entry(s, 'Window start  t1', 't1', 's')
        self.entry(s, 'Window end  t2', 't2', 's')
        s = self.section(f, 'Velocity-distribution output and diagnostic regions')
        self.yesno(s, 'Write velocity distributions', 'pdf')
        self.entry(s, 'Species number', 'ptype_pdf', '(1 = electrons)', width=8)
        self.entry(s, 'y from', 'yl_rg', 'cm')
        self.entry(s, 'y to', 'yr_rg', 'cm')
        self.entry(s, 'z from', 'zl_rg', 'cm')
        self.entry(s, 'z to', 'zr_rg', 'cm')
        self.rg_holder = ttk.Frame(s)
        r = self._row(s)
        self.rg_holder.grid(row=r, column=0, columnspan=5, sticky='ew')
        bf = ttk.Frame(s)
        bf.grid(row=self._row(s), column=0, columnspan=5, sticky='w', padx=8, pady=4)
        ttk.Button(bf, text='+ Add region', command=lambda: self._regions_add(1)).pack(side='left')
        ttk.Button(bf, text='\u2212 Remove last', command=lambda: self._regions_add(-1)).pack(side='left', padx=6)
        s = self.section(f, 'Restart',
                         'These numbers do not set how many threads or processes this run uses - that is done in the '
                         'run script (OMP_NUM_THREADS, mpirun -np).')
        self.choice(s, 'Run type', 'restart', [('no', 'New run'), ('yes', 'Continue (same parallel layout)'),
                                              ('external', 'Import backup from a different layout')], width=42,
                    after=self._restart_visibility)
        self.entry(s, 'OpenMP threads of the backup', 'omp_max', '', width=8)
        self.entry(s, 'MPI ranks of the backup', 'mpi_max', '', width=8)
        self.entry(s, 'Particle duplication / removal', 'np_dup', '(1 = unchanged)', width=8)
        self._restart_frame = s
        self._restart_rows = {k: s._r - 1 - d for k, d in (('omp_max', 2), ('mpi_max', 1), ('np_dup', 0))}

    def _restart_visibility(self):
        if self.case is None or not hasattr(self, '_restart_frame'):
            return
        m = self.case.cond.tok('restart')[:1].lower()
        self._set_rows(self._restart_frame, [self._restart_rows['omp_max'], self._restart_rows['mpi_max']], m == 'e')
        self._set_rows(self._restart_frame, [self._restart_rows['np_dup']], m in ('y', 'e'))

    def _build_regions(self):
        for w in self.rg_holder.winfo_children():
            w.destroy()
        self._drop(r'^rg\d+\.')
        rg = ttk.Frame(self.rg_holder)
        rg.pack(fill='x')
        rg._r = 0
        pdf = self.case.cond.rec['pdf']
        for i in range(len(self.case.cond.regions())):
            for j, nm in ((0, 'from'), (1, 'to')):
                idx = 7 + 2 * i + j
                self.entry(rg, 'Region %d: x %s' % (i + 1, nm), 'rg%d.%d' % (i, j), 'cm',
                           get=lambda idx=idx: pdf.get(idx), put=lambda v, idx=idx: pdf.set(idx, v), dynamic=True)

    def _regions_add(self, d):
        regs = self.case.cond.regions()
        if d > 0 and len(regs) < M.NM_RG:
            regs.append(('0.', '10.'))
        elif d < 0 and regs:
            regs.pop()
        else:
            return
        self.case.cond.set_regions(regs)
        self.changed()

    # ------------------------------------------------------------ geometry
    def _page_geometry(self):
        f = self.tabs['Geometry']
        pw = ttk.PanedWindow(f, orient='horizontal')
        pw.pack(fill='both', expand=True)
        left = ScrollFrame(pw)
        left.canvas.configure(width=600)
        pw.add(left, weight=3)
        self.geom_pw = pw
        L = left.inner
        s = self.section(L, 'Box size (cm)')
        g = lambda i: (lambda: self.case.geom.box.get(i))
        p = lambda i: (lambda v: self.case.geom.box.set(i, v))
        self.signed(s, 'Zero-gradient (Neumann) boundary at x = 0', 'x size', 'g.xmax', 'cm', on_when='negative',
                    width=9, stack=True, get=g(0), put=p(0))
        self.signed(s, 'Draw wall openings round (cylinders along x)', 'y size', 'g.ymax', 'cm', on_when='negative',
                    width=9, stack=True, get=g(1), put=p(1))
        self.signed(s, 'Draw round objects along z', 'z size', 'g.zmax', 'cm', on_when='negative',
                    width=9, stack=True, get=g(2), put=p(2))
        s = self.section(L, 'Wall segments  (drawn in this order - later ones overwrite earlier ones)')
        cols = ('n', 'label', 'kind', 'x', 'y', 'z', 'flags')
        tv = ttk.Treeview(s, columns=cols, show='headings', height=8, selectmode='browse')
        for c, t, w in (('n', '#', 26), ('label', 'Label', 52), ('kind', 'Kind', 112), ('x', 'x range', 92),
                        ('y', 'y range', 92), ('z', 'z range', 92), ('flags', 'Options', 120)):
            tv.heading(c, text=t)
            tv.column(c, width=w, stretch=(c == 'flags'))
        tv.tag_configure('error', background='#ffd9d9')
        tv.tag_configure('warning', background='#fff0c2')
        tv.grid(row=0, column=0, columnspan=5, sticky='ew', padx=6, pady=4)
        tv.bind('<<TreeviewSelect>>', self._seg_tree_select)
        self.seg_tv = tv
        self._seg_cache = None
        bf = ttk.Frame(s)
        bf.grid(row=1, column=0, columnspan=5, sticky='w', padx=6, pady=2)
        for t, cmd in (('Add', self.seg_add), ('Duplicate', self.seg_dup), ('Delete', self.seg_del),
                       ('\u25b2 Up', lambda: self.seg_move(-1)), ('\u25bc Down', lambda: self.seg_move(1))):
            ttk.Button(bf, text=t, command=cmd, width=9).pack(side='left', padx=2)
        s._r = 2
        form = self.section(L, 'Selected segment')
        self.seg_form = form
        kinds = [('+', 'Wall with an opening (plasma window)'), ('-', 'Solid block / plate'), ('0', 'Periodic faces')]

        def seg():
            return self.case.geom.segs[self.cur_seg] if (self.cur_seg is not None and
                                                         self.cur_seg < len(self.case.geom.segs)) else None

        def kget():
            r = seg()
            t = r.get(6) if r else '1'
            return '0' if (t.lstrip('+-') in ('0', '')) else ('-' if M.is_neg(t) else '+')

        def kput(v):
            r = seg()
            if r:
                mag = M.strip_sign(r.get(6)) or '1'
                if mag == '0':
                    mag = '1'
                r.set(6, '0' if v == '0' else ('-' + mag if v == '-' else mag))
        self.choice(form, 'Kind', 'g.f.6k', kinds, get=kget, put=kput, width=40, match='first',
                    after=lambda: None)

        def lget():
            r = seg()
            return M.strip_sign(r.get(6)) if r else ''

        def lput(v):
            r = seg()
            if r and kget() != '0':
                r.set(6, M.with_sign(v, kget() == '-'))
        self.entry(form, 'Label', 'g.f.6', '(1 ... number of wall labels)', width=8, get=lget, put=lput)
        for i, (lab, unit) in enumerate((('x from', 'cm'), ('y from', 'cm'), ('z from', 'cm'),
                                         ('x to', 'cm'), ('y to', 'cm'), ('z to', 'cm'))):
            def gi(i=i):
                r = seg()
                return M.strip_sign(r.get(i)) if r else ''

            def pi(v, i=i):
                r = seg()
                if r:
                    r.set(i, M.with_sign(v, M.is_neg(r.get(i))))
            self.entry(form, lab, 'g.f.%d' % i, unit, width=10, get=gi, put=pi)
        for i, txt in ((3, 'Dielectric surface in an x = const plane (x-to stored negative)'),
                       (4, 'Dielectric surface in a y = const plane (y-to stored negative)'),
                       (5, 'Periodic in z (z-to stored negative)'),
                       (1, 'Toggle round / rectangular for this segment (y-from stored negative)')):
            def fg(i=i):
                r = seg()
                return M.is_neg(r.get(i)) if r else False

            def fp(v, i=i):
                r = seg()
                if r:
                    r.set(i, M.with_sign(r.get(i), v))
            r_ = self._row(form)
            var = tk.BooleanVar()
            cb = ttk.Checkbutton(form, text=txt, variable=var)
            cb.grid(row=r_, column=0, columnspan=4, sticky='w', padx=8)

            def pull(var=var, fg=fg):
                self._loading = True
                var.set(fg())
                self._loading = False

            def push(*_, var=var, fp=fp):
                if not self._loading and seg() is not None:
                    fp(var.get())
                    self.changed()
            var.trace_add('write', push)
            self.pulls.append(('g.flag%d' % i, pull))
            cb.bind('<FocusIn>', lambda e: self.show_help('g.seg'), add='+')
            cb.bind('<Enter>', lambda e: self.show_help('g.seg'), add='+')
        # preview
        pnb = ttk.Notebook(pw)
        pw.add(pnb, weight=2)
        self.prev_nb = pnb
        rsf = ScrollFrame(pnb)
        rsf.canvas.configure(width=470)
        pnb.add(rsf, text=' Cross-sections ')
        f3 = ScrollFrame(pnb)
        pnb.add(f3, text=' 3D view ')
        self._build_3d(f3.inner)
        pnb.bind('<<NotebookTabChanged>>', lambda e: self._prev_tab_changed())
        right = rsf.inner
        ttk.Label(right, text='Cross-sections of the walls',
                  style='Big.TLabel').pack(anchor='w', padx=8, pady=6)
        self.prev_labels, self.prev_read, self._map = {}, {}, {}
        self._read_font = tkfont.nametofont('TkFixedFont').copy()
        self._read_font.configure(size=8)
        self.nreal = None
        self.slice_pos = {'xy': None, 'xz': None, 'yz': None}      # None = default position
        self.slice_ctl = {}
        self._sl_loading = False
        self.r2 = None
        rb = ttk.Frame(right)
        rb.pack(anchor='w', padx=8)
        ttk.Button(rb, text='Reset planes to default', command=self._slices_reset).pack(side='left')
        ttk.Label(right, text='Move each plane with its slider, the \u25c0 \u25b6 buttons (one cell) or by typing a position '
                  '(default: through the middle of the box). Hover over a picture to read the grid indices ix, iy, iz.',
                  style='Hint.TLabel', wraplength=430,
                  justify='left').pack(anchor='w', padx=8, pady=(2, 0))
        for k, title in (('xy', 'XY plane'), ('xz', 'XZ plane'), ('yz', 'YZ plane')):
            fr = ttk.Frame(right)
            fr.pack(anchor='w', padx=8, pady=3)
            cap = ttk.Label(fr, text=title, style='Hint.TLabel')
            cap.pack(anchor='w')
            ax, nm, _, _ = self.PLANE[k]
            ctl = ttk.Frame(fr)
            ctl.pack(anchor='w', fill='x', pady=(2, 2))
            ttk.Label(ctl, text='%s =' % nm).pack(side='left')
            ttk.Button(ctl, text='\u25c0', width=2, command=lambda k=k: self._slice_step(k, -1)).pack(side='left', padx=(4, 0))
            var = tk.DoubleVar(value=0.0)
            sc = ttk.Scale(ctl, from_=0, to=1, variable=var, orient='horizontal', length=190,
                           command=lambda v, k=k: self._slice_scale(k, v))
            sc.pack(side='left', padx=3)
            ttk.Button(ctl, text='\u25b6', width=2, command=lambda k=k: self._slice_step(k, 1)).pack(side='left')
            ev = tk.StringVar()
            e = ttk.Entry(ctl, textvariable=ev, width=8)
            e.pack(side='left', padx=(6, 2))
            e.bind('<Return>', lambda evt, k=k: self._slice_entry(k))
            e.bind('<FocusOut>', lambda evt, k=k: self._slice_entry(k))
            ttk.Label(ctl, text='cm', style='Hint.TLabel').pack(side='left')
            lab = ttk.Label(fr, relief='solid', borderwidth=1)
            lab.pack(anchor='w')
            rd = tk.Label(fr, text=self.HOVER_HINT, font=self._read_font, anchor='w', justify='left', height=3,
                          width=58, foreground='#333333')
            rd.pack(anchor='w')
            lab.bind('<Motion>', lambda e, k=k: self._hover(k, e))
            lab.bind('<Leave>', lambda e, k=k: self._hover_clear(k))
            self.prev_read[k] = rd
            self.prev_labels[k] = (cap, lab)
            self.slice_ctl[k] = dict(var=var, scale=sc, entry=ev)
        self.legend = ttk.Frame(right)
        self.legend.pack(anchor='w', padx=8, pady=4)
        self.prev_note = ttk.Label(right, text='', style='Hint.TLabel', wraplength=380, justify='left')
        self.prev_note.pack(anchor='w', padx=8)
        self._photos = []

    # --- segment table
    def _seg_row(self, k, rec):
        t = [rec.get(i) for i in range(7)]
        lab = M.strip_sign(t[6])
        kind = 'periodic' if lab in ('0', '') else ('solid block' if M.is_neg(t[6]) else 'wall + opening')
        fl = []
        if M.is_neg(t[3]):
            fl.append('diel. X')
        if M.is_neg(t[4]):
            fl.append('diel. Y')
        if M.is_neg(t[5]):
            fl.append('periodic z')
        if M.is_neg(t[1]):
            fl.append('toggle round')
        s = M.strip_sign
        return (k + 1, lab, kind, '%s \u2192 %s' % (s(t[0]), s(t[3])), '%s \u2192 %s' % (s(t[1]), s(t[4])),
                '%s \u2192 %s' % (s(t[2]), s(t[5])), ', '.join(fl))

    def _refresh_seg_tree(self):
        segs = self.case.geom.segs
        lv = {}
        for i in self.issues:
            m = re.match(r'^g\.seg(\d+)\.', i.key or '')
            if m:
                k = int(m.group(1))
                if k not in lv or M.Issue.ORDER[i.level] < M.Issue.ORDER[lv[k]]:
                    lv[k] = i.level
        rows = [(self._seg_row(k, r), lv.get(k)) for k, r in enumerate(segs)]
        if rows == self._seg_cache:
            return
        self._seg_cache = rows
        tv = self.seg_tv
        tv.unbind('<<TreeviewSelect>>')
        tv.delete(*tv.get_children())
        for vals, level in rows:
            tag = 'error' if level in ('fatal', 'error') else ('warning' if level == 'warning' else '')
            tv.insert('', 'end', values=vals, tags=(tag,) if tag else ())
        if self.cur_seg is not None and self.cur_seg < len(rows):
            tv.selection_set(tv.get_children()[self.cur_seg])
        tv.bind('<<TreeviewSelect>>', self._seg_tree_select)

    def _seg_tree_select(self, _=None):
        sel = self.seg_tv.selection()
        if sel:
            self.select_segment(self.seg_tv.index(sel[0]))

    def select_segment(self, k):
        n = len(self.case.geom.segs)
        self.cur_seg = k if (k is not None and 0 <= k < n) else None
        if self.cur_seg is not None:
            kids = self.seg_tv.get_children()
            if kids and self.seg_tv.selection() != (kids[self.cur_seg],):
                self.seg_tv.selection_set(kids[self.cur_seg])
        for key, p in self.pulls:
            if key.startswith('g.f.') or key.startswith('g.flag'):
                p()
        self.recheck()

    def seg_add(self):
        g = self.case.geom
        cur = g.segs[self.cur_seg] if self.cur_seg is not None and self.cur_seg < len(g.segs) else None
        r = g.add_segment(after=cur)
        self._seg_cache = None
        self.select_segment(g.segs.index(r))
        self.changed()

    def seg_dup(self):
        g = self.case.geom
        if self.cur_seg is None or self.cur_seg >= len(g.segs):
            return
        cur = g.segs[self.cur_seg]
        r = g.add_segment(after=cur, toks=list(cur.toks))
        r.comment = cur.comment
        self._seg_cache = None
        self.select_segment(g.segs.index(r))
        self.changed()

    def seg_del(self):
        g = self.case.geom
        if self.cur_seg is None or self.cur_seg >= len(g.segs):
            return
        g.remove_segment(g.segs[self.cur_seg])
        self._seg_cache = None
        self.select_segment(min(self.cur_seg, len(g.segs) - 1) if g.segs else None)
        self.changed()

    def seg_move(self, d):
        g = self.case.geom
        if self.cur_seg is None or self.cur_seg >= len(g.segs):
            return
        g.move_segment(g.segs[self.cur_seg], d)
        self._seg_cache = None
        self.select_segment(max(0, min(len(g.segs) - 1, self.cur_seg + d)))
        self.changed()

    # --- preview
    def draw_preview(self):
        self.preview_stale = False
        self.stale3 = True
        if self._on_3d_tab():
            self._draw3d()
        for lab in self.legend.winfo_children():
            lab.destroy()
        if np is None:
            self.prev_note.configure(text='Install numpy to enable the preview  (pip install numpy).')
            return
        p, _ = M.decode(self.case)
        box, segs, _ = M.decode_geometry(self.case)
        if None in box or 0 in box:
            self.prev_note.configure(text='Fill in a valid box size to see the preview.')
            return
        good = [s for s in segs if None not in s]
        nr = [p['nx'], p['ny'], p['nz']]
        self.nreal = nr if all(isinstance(v, int) and v >= 1 for v in nr) else None
        n = [p['nx'] or 64, p['ny'] or 64, p['nz'] or 64]
        n = [max(4, v) for v in n]
        N = n[0] * n[1] * n[2]
        if N > 3e6:
            f = (3e6 / N) ** (1 / 3)
            n = [max(4, int(v * f)) for v in n]
        holes = None
        if p['holes'] and None not in (p['Lhy'], p['Lhz'], p['nhy'], p['nhz'], p['ind_g']):
            holes = dict(on=True, Lhy=p['Lhy'], Lhz=p['Lhz'], nhy=p['nhy'], nhz=p['nhz'], ind_g=p['ind_g'])
        try:
            r = M.rasterize(box, good, n, holes)
        except Exception as e:
            self.prev_note.configure(text='Preview failed: %s' % e)
            return
        self.r2, self.n2 = r, n
        self.h2 = [r['L'][i] / n[i] for i in range(3)]
        self._photos, self._used = {}, {}
        for k, ctl in self.slice_ctl.items():
            ctl['scale'].configure(to=r['L'][self.PLANE[k][0]])
        for k in ('xy', 'xz', 'yz'):
            self._render_slice(k)
        note = 'Approximate view on a %dx%dx%d grid, using the same drawing rules as the code.' % tuple(n)
        if self.nreal and list(self.nreal) != list(n):
            note += ('\nThe simulation grid is %dx%dx%d: indices read with the mouse refer to it, while the '\
                     'content shown is that of the nearest preview node.' % tuple(self.nreal))
        if r['notes']:
            note += '\n' + ' '.join(r['notes'])
        self.prev_note.configure(text=note)

    # --- cross-section planes: each can be moved without recomputing the volume
    PLANE = {'xy': (2, 'z', 0, 1), 'xz': (1, 'y', 0, 2), 'yz': (0, 'x', 1, 2)}   # normal axis, name, horizontal, vertical

    def _default_pos(self, k):
        pos = self.r2['pos']
        return {'xy': pos[0], 'xz': pos[1], 'yz': pos[2]}[k]

    def _slice_index(self, k):
        ax = self.PLANE[k][0]
        pos = self.slice_pos[k]
        if pos is None:
            pos = self._default_pos(k)
        pos = min(max(pos, 0.0), self.r2['L'][ax])
        return max(0, min(int(round(pos / self.h2[ax])), self.n2[ax]))

    def _render_slice(self, k, update_scale=True):
        if self.r2 is None or np is None:
            return
        r = self.r2
        vol = r['vol']
        ax, nm, ha, va = self.PLANE[k]
        i = self._slice_index(k)
        arr = vol[:, :, i] if ax == 2 else (vol[:, i, :] if ax == 1 else vol[i, :, :])
        Lh, Lv = r['L'][ha], r['L'][va]
        sc = min(400.0 / Lh, 170.0 / Lv)
        W, H = max(2, int(Lh * sc)), max(2, int(Lv * sc))
        ii = np.clip(np.rint((np.arange(W) + 0.5) / W * (arr.shape[0] - 1)).astype(int), 0, arr.shape[0] - 1)
        jj = np.clip(np.rint((H - 1 - np.arange(H) + 0.5) / H * (arr.shape[1] - 1)).astype(int), 0, arr.shape[1] - 1)
        img = arr[np.ix_(ii, jj)].T
        cmap = {v: self._bc_color(v, r['dtype']) for v in np.unique(img)}
        rows = ['{' + ' '.join(cmap[v] for v in row) + '}' for row in img.tolist()]
        ph = tk.PhotoImage(width=W, height=H)
        ph.put(' '.join(rows))
        self._photos[k] = ph
        self._used[k] = set(int(v) for v in np.unique(img))
        snapped = i * self.h2[ax]
        self._map[k] = dict(ii=ii, jj=jj, W=W, H=H, arr=arr, Lh=Lh, Lv=Lv, fixed=snapped)
        ctl = self.slice_ctl[k]
        self._sl_loading = True
        ctl['entry'].set('%.4g' % snapped)
        if update_scale:
            ctl['var'].set(snapped)
        self._sl_loading = False
        cap, lab = self.prev_labels[k]
        lab.configure(image=ph)
        cap.configure(text='%s plane at %s = %.4g cm   (%s \u2192 right, %s \u2191)   %.3g x %.3g cm' % (
            k.upper(), nm, snapped, 'xyz'[ha], 'xyz'[va], Lh, Lv))
        self._update_legend()

    def _update_legend(self):
        for w in self.legend.winfo_children():
            w.destroy()
        if self.r2 is None:
            return
        dt = self.r2['dtype']
        used = set().union(*self._used.values()) if self._used else set()
        for v, name in ((-1, 'plasma / open'), (0, 'periodic'), (-2, 'zero-gradient (Neumann)')):
            if v in used:
                self._legend_item(self._bc_color(v, dt), name)
        for v in sorted(x for x in used if x > 0):
            t = dt.get(v, 1)
            self._legend_item(self._bc_color(v, dt), 'label %d  (%s)' % (v, 'metal' if t == 1 else 'dielectric, type %d' % t))

    # --- reading grid indices with the mouse
    HOVER_HINT = 'Hover over the picture to read the\nnode indices ix, iy, iz of the simulation grid.'

    def _real_index(self, pos, axis):
        """Node number (1 ... n+1, as in the Fortran code: x = (ix-1)*dx) of the simulation grid at `pos` (cm)."""
        L = self.r2['L'][axis]
        nr = self.nreal[axis] if self.nreal else self.n2[axis]
        return max(1, min(int(round(pos / (L / nr))) + 1, nr + 1)), nr + 1

    def _content_name(self, v):
        v = int(v)
        if v == -1:
            return 'plasma / open'
        if v == 0:
            return 'periodic face'
        if v == -2:
            return 'zero-gradient (Neumann)'
        t = self.r2['dtype'].get(v, 1)
        return 'label %d (%s)' % (v, 'metal' if t == 1 else 'dielectric, type %d' % t)

    def _hover(self, k, ev):
        m = self._map.get(k)
        if not m or self.r2 is None:
            return
        lab = self.prev_labels[k][1]
        px = ev.x - (lab.winfo_width() - m['W']) // 2
        py = ev.y - (lab.winfo_height() - m['H']) // 2
        if not (0 <= px < m['W'] and 0 <= py < m['H']):
            return self._hover_clear(k)
        ax, nm, ha, va = self.PLANE[k]
        pos = [0.0, 0.0, 0.0]
        pos[ha] = (px + 0.5) / m['W'] * m['Lh']
        pos[va] = (m['H'] - py - 0.5) / m['H'] * m['Lv']          # vertical axis points up, row 0 = top
        pos[ax] = m['fixed']
        idx = [self._real_index(pos[a], a) for a in range(3)]
        what = self._content_name(m['arr'][m['ii'][px], m['jj'][py]])
        parts = []
        for a, name in enumerate('xyz'):
            tag = '  (plane)' if a == ax else ''
            parts.append('i%s = %d/%d%s' % (name, idx[a][0], idx[a][1], tag))
        self.prev_read[k].configure(text='   '.join(parts) + '\nx = %.4g   y = %.4g   z = %.4g cm\n%s' % (
            pos[0], pos[1], pos[2], what))

    def _hover_clear(self, k):
        self.prev_read[k].configure(text=self.HOVER_HINT)

    def _slice_scale(self, k, v):
        if self._sl_loading or self.r2 is None:
            return
        self.slice_pos[k] = float(v)
        self._render_slice(k, update_scale=False)

    def _slice_entry(self, k):
        if self._sl_loading or self.r2 is None:
            return
        try:
            v = float(self.slice_ctl[k]['entry'].get().replace('d', 'e').replace('D', 'e'))
        except ValueError:
            self._render_slice(k)                    # restore the previous position
            return
        self.slice_pos[k] = v
        self._render_slice(k)

    def _slice_step(self, k, d):
        if self.r2 is None:
            return
        ax = self.PLANE[k][0]
        i = max(0, min(self._slice_index(k) + d, self.n2[ax]))
        self.slice_pos[k] = i * self.h2[ax]
        self._render_slice(k)

    def _slices_reset(self):
        self.slice_pos = {'xy': None, 'xz': None, 'yz': None}
        for k in ('xy', 'xz', 'yz'):
            self._render_slice(k)

    @staticmethod
    def _bc_color(v, dtype):
        v = int(v)
        if v == -1:
            return '#cfe8ff'
        if v == 0:
            return '#d7d7d7'
        if v == -2:
            return '#e0d0f0'
        c = PALETTE[(v - 1) % len(PALETTE)]
        return lighten(c) if dtype.get(v, 1) > 1 else c

    def _legend_item(self, color, text):
        fr = ttk.Frame(self.legend)
        fr.pack(anchor='w')
        tk.Label(fr, background=color, width=2, relief='solid', borderwidth=1).pack(side='left', padx=2)
        ttk.Label(fr, text=text).pack(side='left')


    # ------------------------------------------------------------ 3D view
    DETAIL_LABELS = ['Low (1/8 of the full grid)', 'Medium (1/2 of the full grid)', 'High (full grid)', 'Custom']
    PLASMA_RGB = (77, 166, 255)
    MAX_NODES3 = 8e7
    CW, CH = 450, 360

    def _on_3d_tab(self):
        try:
            return self.prev_nb.tab(self.prev_nb.select(), 'text').strip() == '3D view'
        except (tk.TclError, AttributeError):
            return False

    def _prev_tab_changed(self):
        if self.case is not None and self._on_3d_tab() and self.stale3:
            self.root.after(30, self._draw3d)

    def _build_3d(self, f):
        self.v3 = dict(az=-50.0, el=28.0, zoom=1.0)
        self.vol3, self.dt3, self.stale3 = None, {}, True
        self.hidden3, self.labels3, self._busy3, self._drag3 = set(), [], False, None
        ttk.Label(f, text='3D view of the walls and electrodes', style='Big.TLabel').pack(anchor='w', padx=8, pady=6)
        self.cv3 = tk.Canvas(f, width=self.CW, height=self.CH, background='#f4f4f0', highlightthickness=1,
                             highlightbackground='#999999', cursor='fleur')
        self.cv3.pack(padx=8, pady=(0, 2))
        self.cv3.bind('<ButtonPress-1>', self._d3_press)
        self.cv3.bind('<B1-Motion>', self._d3_move)
        self.cv3.bind('<ButtonRelease-1>', self._d3_release)
        self.cv3.bind('<MouseWheel>', lambda e: self._d3_zoom(1.1 if e.delta > 0 else 1 / 1.1))
        self.cv3.bind('<Button-4>', lambda e: self._d3_zoom(1.1))
        self.cv3.bind('<Button-5>', lambda e: self._d3_zoom(1 / 1.1))
        ttk.Label(f, text='Drag to rotate, mouse wheel to zoom. Open space is transparent; tick "plasma" to show it.',
                  style='Hint.TLabel', wraplength=430, justify='left').pack(anchor='w', padx=8)
        bf = ttk.Frame(f)
        bf.pack(anchor='w', padx=8, pady=4)
        for t, az, el in (('Isometric', -50, 28), ('Front', -90, 0), ('Side', 0, 0), ('Top', -90, 89)):
            ttk.Button(bf, text=t, width=10, command=lambda az=az, el=el: self._d3_preset(az, el)).pack(
                side='left', padx=2)
        ttk.Button(bf, text='Reset zoom', width=11, command=lambda: self._d3_zoom(None)).pack(side='left', padx=2)
        cf = ttk.LabelFrame(f, text='Options')
        cf.pack(fill='x', padx=8, pady=4)
        r = ttk.Frame(cf)
        r.pack(fill='x', padx=6, pady=3)
        ttk.Label(r, text='Resolution').pack(side='left')
        self.detail_var = tk.StringVar(value=self.DETAIL_LABELS[0])
        cb = ttk.Combobox(r, textvariable=self.detail_var, values=self.DETAIL_LABELS, state='readonly', width=26)
        cb.pack(side='left', padx=6)
        cb.bind('<<ComboboxSelected>>', self._detail_changed)
        self.detail_info = ttk.Label(cf, text='', style='Hint.TLabel', wraplength=420, justify='left')
        self.detail_info.pack(anchor='w', padx=6)
        self.custom_row = ttk.Frame(cf)
        self.custom_vars = [tk.StringVar(value='32'), tk.StringVar(value='32'), tk.StringVar(value='32')]
        ttk.Label(self.custom_row, text='Number of cells').pack(side='left')
        for nm, v in zip(('nx', 'ny', 'nz'), self.custom_vars):
            ttk.Label(self.custom_row, text=nm).pack(side='left', padx=(8, 2))
            ttk.Entry(self.custom_row, textvariable=v, width=6).pack(side='left')
            v.trace_add('write', lambda *a: self._custom_changed())
        self._cjob = None
        r = ttk.Frame(cf)
        r.pack(fill='x', padx=6, pady=3)
        ttk.Label(r, text='Cut-away').pack(side='left')
        self.cut_axis = tk.StringVar(value='Off')
        cb = ttk.Combobox(r, textvariable=self.cut_axis, values=['Off', 'x', 'y', 'z'], state='readonly', width=5)
        cb.pack(side='left', padx=6)
        cb.bind('<<ComboboxSelected>>', lambda e: self._render3())
        self.cut_var = tk.DoubleVar(value=60.0)
        sc = ttk.Scale(r, from_=0, to=100, variable=self.cut_var, orient='horizontal', length=150,
                       command=lambda v: self._render3(0.5))
        sc.pack(side='left', padx=6)
        sc.bind('<ButtonRelease-1>', lambda e: self._render3())
        self.cut_flip = tk.BooleanVar(value=False)
        ttk.Checkbutton(r, text='other side', variable=self.cut_flip, command=self._render3).pack(side='left')
        ttk.Label(cf, text='The cut-away removes everything beyond the chosen position along that axis, '
                  'so that you can look inside closed boxes.', style='Hint.TLabel', wraplength=420,
                  justify='left').pack(anchor='w', padx=6, pady=(0, 4))
        self.leg3 = ttk.LabelFrame(f, text='Show / hide')
        self.leg3.pack(fill='x', padx=8, pady=4)
        pf = ttk.Frame(self.leg3)
        pf.pack(anchor='w', padx=6)
        tk.Label(pf, background='#%02x%02x%02x' % self.PLASMA_RGB, width=2, relief='solid',
                 borderwidth=1).pack(side='left', padx=2)
        self.plasma3 = tk.BooleanVar(value=False)
        ttk.Checkbutton(pf, variable=self.plasma3, command=self._render3,
                        text='plasma  (translucent blue skin)').pack(side='left')
        self.leg3_labels = ttk.Frame(self.leg3)
        self.leg3_labels.pack(anchor='w', fill='x')
        self.note3 = ttk.Label(f, text='', style='Hint.TLabel', wraplength=430, justify='left')
        self.note3.pack(anchor='w', padx=8, pady=(0, 8))

    def _d3_preset(self, az, el):
        self.v3['az'], self.v3['el'] = float(az), float(el)
        self._render3()

    def _d3_zoom(self, f):
        self.v3['zoom'] = 1.0 if f is None else max(0.4, min(5.0, self.v3['zoom'] * f))
        self._render3()

    def _d3_press(self, e):
        self._drag3 = (e.x, e.y, self.v3['az'], self.v3['el'])

    def _d3_move(self, e):
        if not self._drag3:
            return
        x0, y0, a0, e0 = self._drag3
        self.v3['az'] = a0 - (e.x - x0) * 0.5
        self.v3['el'] = max(-89.0, min(89.0, e0 + (e.y - y0) * 0.5))
        if not self._busy3:
            self._busy3 = True
            self.root.after(25, lambda: (self._render3(0.4), setattr(self, '_busy3', False)))

    def _d3_release(self, e):
        self._drag3 = None
        self._render3()

    def _detail_changed(self, _=None):
        if self.detail_var.get() == 'Custom':
            g = getattr(self, 'grid3', None)
            if g:
                for v, n in zip(self.custom_vars, g):
                    v.set(str(n))
            self.custom_row.pack(anchor='w', padx=6, pady=3, after=self.detail_info)
        else:
            self.custom_row.pack_forget()
        self._draw3d()

    def _custom_changed(self):
        if self.detail_var.get() != 'Custom':
            return
        if self._cjob:
            self.root.after_cancel(self._cjob)
        self._cjob = self.root.after(700, self._draw3d)

    def _grid3(self, p, L):
        """Cells (nx, ny, nz) used for the 3D picture, from the resolution choice."""
        mode = self.detail_var.get()
        if mode == 'Custom':
            n = []
            for v in self.custom_vars:
                try:
                    n.append(max(2, min(2000, int(v.get()))))
                except ValueError:
                    return None
            return n
        frac = 1.0 / 8 if mode.startswith('Low') else (0.5 if mode.startswith('Medium') else 1.0)
        full = [p['nx'], p['ny'], p['nz']]
        if all(isinstance(v, int) and v >= 2 for v in full):
            return [v if frac == 1.0 else max(4, int(round(v * frac))) for v in full]
        N = 32 if frac < 0.2 else (64 if frac < 1 else 128)        # grid sizes not valid yet: proportional fallback
        return [max(4, int(round(l / max(L) * N))) for l in L]

    def _draw3d(self):
        """(Re)build the voxel volume from the current geometry, then render."""
        self.stale3 = False
        c = self.cv3
        if np is None:
            c.delete('all')
            c.create_text(self.CW / 2, self.CH / 2, text='Install numpy to enable the 3D view\n(pip install numpy)',
                          justify='center')
            return
        box, segs, _ = M.decode_geometry(self.case)
        p, _ = M.decode(self.case)
        if None in box or 0 in box:
            c.delete('all')
            c.create_text(self.CW / 2, self.CH / 2, text='Fill in a valid box size.')
            return
        good = [s for s in segs if None not in s]
        L = [abs(b) for b in box]
        n = self._grid3(p, L)
        if n is None:
            self.detail_info.configure(text='Enter whole numbers for the number of cells.')
            return
        nodes = (n[0] + 1) * (n[1] + 1) * (n[2] + 1)
        full = [p['nx'], p['ny'], p['nz']]
        fulltxt = ('   (full grid of the simulation: %d x %d x %d)' % tuple(full)) if None not in full else ''
        if nodes > self.MAX_NODES3:
            c.delete('all')
            c.create_text(self.CW / 2, self.CH / 2, justify='center', text='%d x %d x %d cells is too large for the '
                          '3D view.\nChoose a lower resolution.' % tuple(n))
            self.detail_info.configure(text='Too many cells (%.0f million nodes).' % (nodes / 1e6))
            return
        self.grid3 = n
        self.h3 = [L[i] / n[i] for i in range(3)]
        info = '%d x %d x %d cells%s' % (n[0], n[1], n[2], fulltxt)
        if nodes > 4e6:
            info += '\nLarge grid: drawing takes a few seconds; use a lower resolution while rotating.'
        self.detail_info.configure(text=info)
        holes = None
        if p['holes'] and None not in (p['Lhy'], p['Lhz'], p['nhy'], p['nhz'], p['ind_g']):
            holes = dict(on=True, Lhy=p['Lhy'], Lhz=p['Lhz'], nhy=p['nhy'], nhz=p['nhz'], ind_g=p['ind_g'])
        try:
            r = M.rasterize(box, good, n, holes)
        except Exception as e:
            self.note3.configure(text='3D view failed: %s' % e)
            return
        self.vol3, self.dt3, self.L3 = r['vol'].copy(), r['dtype'], L
        labs = sorted(int(v) for v in np.unique(self.vol3) if v > 0)
        if labs != self.labels3:
            self.labels3 = labs
            for w in self.leg3_labels.winfo_children():
                w.destroy()
            self._leg_vars = {}
            for k in labs:
                fr = ttk.Frame(self.leg3_labels)
                fr.pack(anchor='w', padx=6)
                tk.Label(fr, background=self._bc_color(k, self.dt3), width=2, relief='solid',
                         borderwidth=1).pack(side='left', padx=2)
                var = tk.BooleanVar(value=k not in self.hidden3)
                self._leg_vars[k] = var
                t = self.dt3.get(k, 1)
                ttk.Checkbutton(fr, variable=var, command=self._leg_changed,
                                text='label %d  (%s)' % (k, 'metal' if t == 1 else 'dielectric')).pack(side='left')
            if not labs:
                ttk.Label(self.leg3_labels, text='(no walls in this geometry)', style='Hint.TLabel').pack(anchor='w', padx=6)
        self.note3.configure(text='Drawn on %d x %d x %d cells: at coarse resolution thin plates look thicker than they are.' % tuple(n))
        self._render3()

    def _leg_changed(self):
        self.hidden3 = {k for k, v in self._leg_vars.items() if not v.get()}
        self._render3()

    @staticmethod
    def _hex2rgb(h):
        return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))

    def _render3(self, quality=1.0):
        if self.vol3 is None or np is None:
            return
        if self.stale3:
            return self._draw3d()
        colors = {k: self._hex2rgb(self._bc_color(k, self.dt3)) for k in self.labels3 if k not in self.hidden3}
        ax = self.cut_axis.get()
        clip = None
        if ax in ('x', 'y', 'z'):
            clip = ('xyz'.index(ax), self.cut_var.get() / 100.0, self.cut_flip.get())
        ppm, V = M.render_voxels(self.vol3, colors, (self.CW, self.CH), self.v3['az'], self.v3['el'],
                                 self.v3['zoom'], clip, quality, spacing=self.h3,
                                 plasma=self.PLASMA_RGB if self.plasma3.get() else None)
        try:
            self.ph3 = tk.PhotoImage(data=ppm)
        except tk.TclError:
            return
        c = self.cv3
        c.delete('all')
        c.create_image(0, 0, image=self.ph3, anchor='nw')
        sh = self.vol3.shape
        lo, hi = -0.5, [s - 0.5 for s in sh]
        corners = np.array([[x, y, z] for x in (lo, hi[0]) for y in (lo, hi[1]) for z in (lo, hi[2])])
        px, py, _ = V.project(corners)
        for i in range(8):
            for j in range(i + 1, 8):
                if bin(i ^ j).count('1') == 1:
                    c.create_line(px[i], py[i], px[j], py[j], fill='#8a8a8a', dash=(2, 3))
        ox, oy = 38, self.CH - 34
        for k, (nm, col) in enumerate((('x', '#c62828'), ('y', '#2e7d32'), ('z', '#1565c0'))):
            v = V.R[:, k]
            ex, ey = ox + 28 * v[0], oy - 28 * v[1]
            c.create_line(ox, oy, ex, ey, fill=col, width=2, arrow='last')
            c.create_text(ex + (6 if ex >= ox else -6), ey + (0 if ey != oy else 0), text=nm, fill=col,
                          font=('TkDefaultFont', 10, 'bold'))
        Lx, Ly, Lz = self.L3
        c.create_text(self.CW - 6, 8, anchor='ne', fill='#555555', font=('TkDefaultFont', 9),
                      text='box %.3g x %.3g x %.3g cm' % (Lx, Ly, Lz))

    # ------------------------------------------------------------ potentials
    def _page_potentials(self):
        f = self.tabs['Potentials'].inner
        s = self.section(f, 'Electrodes and walls (boundary.inp)',
                         'One potential per wall label. The label numbers are the ones used in the Geometry tab. '
                         'Label 0 is reserved for periodic faces and has no potential.')
        self.signed(s, 'No 3D time averages (negative)', 'Number of wall labels', 'ngrid', '', on_when='negative',
                    width=8)
        self.pot_holder = ttk.Frame(f)
        self.pot_holder.pack(fill='x')
        self.pot_usage = {}

    def _build_potentials(self):
        for w in self.pot_holder.winfo_children():
            w.destroy()
        self._drop(r'^b\.V\d+')
        self.pot_usage = {}
        s = self.section(self.pot_holder, 'Potentials')
        s.columnconfigure(4, weight=1)
        for k, rec in enumerate(self.case.bnd.vals):
            self.entry(s, 'Label %d' % (k + 1), 'b.V%d' % (k + 1), 'V', width=12,
                       get=lambda rec=rec: rec.get(0), put=lambda v, rec=rec: rec.set(0, v), dynamic=True)
            lab = ttk.Label(s, text='', style='Hint.TLabel')
            lab.grid(row=s._r - 1, column=4, sticky='w', padx=8)
            self.pot_usage[k + 1] = lab

    def _refresh_potential_usage(self):
        if not self.pot_usage:
            return
        box, segs, _ = M.decode_geometry(self.case)
        p, _ = M.decode(self.case)
        use = {}
        for k, s in enumerate(segs):
            if None in s:
                continue
            a = abs(s[6])
            kind = 'solid block' if s[6] < 0 else 'wall'
            die = ', dielectric' if (s[3] < 0 or s[4] < 0) else ''
            use.setdefault(a, []).append('segment %d (%s%s)' % (k + 1, kind, die))
        for lab, w in self.pot_usage.items():
            txt = '; '.join(use.get(lab, [])) or 'not used in geometry'
            role = []
            emit = (p['gam_sec'] or 0) > 0 or abs(p['opt_inj'] or 0) == 4
            if p['igrid_sec'] == lab and emit:
                role.append('emitting electrode')
            if p['ind_g'] == lab and (p['holes'] or (p['jne'] or 0) > 0):
                role.append('grid / plasma electrode')
            if role:
                txt += '   [' + ', '.join(role) + ']'
            w.configure(text=txt, foreground='#b36b00' if lab not in use else '#555555')

    # ------------------------------------------------------------ dynamic parts
    def _rebuild_dynamic(self):
        self._build_bentries()
        self._build_regions()
        self._build_potentials()
        self._build_coil()
        self._built = (len(self.case.cond.bfiles), len(self.case.cond.regions()), len(self.case.bnd.vals),
                       len(self.case.coil.turns))
        self._refresh_dat_list()

    def _sync_dynamic(self):
        self.case.sync_boundary()
        cur = (len(self.case.cond.bfiles), len(self.case.cond.regions()), len(self.case.bnd.vals),
               len(self.case.coil.turns))
        if cur != getattr(self, '_built', None):
            old = getattr(self, '_built', (None, None, None, None))
            if cur[0] != old[0]:
                self._build_bentries()
            if cur[1] != old[1]:
                self._build_regions()
            if cur[2] != old[2]:
                self._build_potentials()
            if cur[3] != old[3]:
                self._build_coil()
            self._built = cur
            for _, p in self.dyn:
                p()
        for k in range(len(self.case.cond.bfiles)):
            self._b_visibility(k)
        self._rf_visibility()
        self._restart_visibility()

    def _refresh_dat_list(self):
        vals = set()
        for d in (self.case.folder, os.path.dirname(self.case.folder)):
            try:
                vals |= {x for x in os.listdir(d) if x.lower().endswith('.dat')}
            except OSError:
                pass
        self.rname_box.configure(values=sorted(vals))


def find_default_input_dir():
    """input_dir/ next to / above the working directory or this script (also one level below)."""
    def ok(p):
        return all(os.path.isfile(os.path.join(p, f)) for f in M.Case.FILES)
    starts = [os.getcwd(), os.path.dirname(os.path.abspath(__file__))]
    for st in starts:
        if os.path.basename(st) == 'input_dir' and ok(st):
            return st
        p = st
        for _ in range(5):                       # this folder and up to 4 parents
            c = os.path.join(p, 'input_dir')
            if ok(c):
                return c
            parent = os.path.dirname(p)
            if parent == p:
                break
            p = parent
        try:                                     # e.g. run from a folder that contains the code folder
            for d in sorted(os.listdir(st)):
                c = os.path.join(st, d, 'input_dir')
                if ok(c):
                    return c
        except OSError:
            pass
    return None


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else find_default_input_dir()
    root = tk.Tk()
    app = App(root, folder)
    if folder is None:
        root.after(200, app.open_folder)
    root.mainloop()


if __name__ == '__main__':
    main()
