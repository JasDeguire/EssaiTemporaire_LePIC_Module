// LePIC+ web GUI - front end of web_gui.py (no framework, no external files)
'use strict';

const $ = (s, r = document) => r.querySelector(s);
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const ORDER = { fatal: 0, error: 1, warning: 2, info: 3 };
const MARK = { fatal: '✖', error: '✖', warning: '⚠', info: 'ℹ' };
const LEVEL_NAME = { fatal: 'BLOCKING', error: 'problem', warning: 'warning', info: 'note' };
const CMAPS = ['viridis', 'inferno', 'plasma', 'magma', 'cividis', 'turbo', 'jet', 'RdBu_r', 'coolwarm', 'gray'];

async function api(path, body) {
  const opt = body === undefined ? {} : { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) };
  const r = await fetch(path, opt);
  const j = await r.json();
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}

// ---------------------------------------------------------------- pages
document.querySelectorAll('header .page').forEach(b => b.onclick = () => {
  document.querySelectorAll('header .page').forEach(x => x.classList.toggle('on', x === b));
  document.querySelectorAll('.page-body').forEach(x => x.classList.toggle('on', x.id === 'page-' + b.dataset.page));
  if (b.dataset.page === 'out2d' && !O.scan) oOpen();
  if (b.dataset.page === 'out3d' && !T.loaded) tPlot(true);
});

// =====================================================================
//  INPUT EDITOR
// =====================================================================
let S = null;            // server state
let HELP = {};
let curTab = 'Summary';
let selSeg = 0;
let geoFrac = { xy: 0.5, xz: 0.5, yz: 0.5 };

function status(msg, bad) { const s = $('#status'); s.textContent = msg || ''; s.className = bad ? 'bad' : ''; }

async function initInputs() {
  HELP = await api('/api/help');
  const c = await api('/api/cases');
  const sel = $('#case');
  sel.innerHTML = c.cases.map(x => `<option value="${esc(x.folder)}">${esc(x.name)}</option>`).join('');
  sel.onchange = async () => {
    if (S && S.modified && !confirm('Discard the unsaved changes of ' + S.folder + '?')) { sel.value = S.folder; return; }
    await load(sel.value);
  };
  showHelp(null);
  if (c.cases.length) await load(c.cases[0].folder);
}

async function load(folder) {
  try { setState(await api('/api/load', { folder })); selSeg = 0; status('Loaded ' + folder); }
  catch (e) { status(e.message, true); }
}

function setState(st) {
  S = st;
  if (st.sel !== null && st.sel !== undefined) selSeg = st.sel;
  $('#folder').textContent = st.folder + (st.modified ? '   — UNSAVED CHANGES' : '') + '   [' + st.files.join(', ') + ']';
  $('#btn-activate').disabled = !st.can_activate;
  if (st.message) status(st.message);
  renderTabs();
  renderTab();
}

function worst() {
  const w = {};
  for (const i of S.issues) if (i.key && (!(i.key in w) || ORDER[i.level] < ORDER[w[i.key]])) w[i.key] = i.level;
  return w;
}

function tabNames() { return ['Summary'].concat(S.tabs.map(t => t.name)); }

function renderTabs() {
  const n = {};
  for (const i of S.issues) if (i.level === 'fatal' || i.level === 'error') n[i.tab || 'Summary'] = (n[i.tab || 'Summary'] || 0) + 1;
  $('#tabs').innerHTML = tabNames().map(t => `<button data-t="${esc(t)}" class="${t === curTab ? 'on' : ''}">${esc(t)}${n[t] ? `<span class="cnt">${n[t]}</span>` : ''}</button>`).join('');
  document.querySelectorAll('#tabs button').forEach(b => b.onclick = () => { curTab = b.dataset.t; renderTabs(); renderTab(); });
}

function fieldRow(f, w) {
  if (f.hidden) return '';
  const lv = w[f.ikey];
  const cls = lv === 'fatal' || lv === 'error' ? 'err' : lv === 'warning' ? 'warn' : '';
  const mark = lv ? `<span class="lv-${lv}">${MARK[lv]}</span>` : '';
  const data = `data-key="${esc(f.key)}" data-fid="${esc(f.fid)}" data-help="${esc(f.help || '')}" data-ikey="${esc(f.ikey)}"`;
  let ctl = '';
  if (f.type === 'text') {
    const dl = f.datalist ? `list="dl-${esc(f.fid)}"` : '';
    ctl = `<input type="text" class="${cls}" ${data} ${dl} value="${esc(f.value)}" size="14">` +
      (f.datalist ? `<datalist id="dl-${esc(f.fid)}">${f.datalist.map(d => `<option value="${esc(d)}">`).join('')}</datalist>` : '');
  } else if (f.type === 'signed') {
    const neg = String(f.value).startsWith('-');
    const on = f.on_when === 'positive' ? !neg : neg;
    const mag = String(f.value).replace(/^[+-]/, '');
    const dis = f.on_when === 'positive' && !on ? 'disabled' : '';
    ctl = `<label><input type="checkbox" class="sg-c" ${data} data-on="${f.on_when}" ${on ? 'checked' : ''}> ${esc(f.check)}</label>
           <input type="text" class="sg-v ${cls}" ${data} data-on="${f.on_when}" value="${esc(mag)}" size="11" ${dis}>`;
  } else if (f.type === 'yesno') {
    const v = typeof f.value === 'boolean' ? f.value : /^[yY]/.test(f.value || '');
    ctl = `<input type="checkbox" class="yn" ${data} data-bool="${typeof f.value === 'boolean' ? 1 : 0}" ${v ? 'checked' : ''}>`;
  } else if (f.type === 'choice') {
    const t = String(f.value || '').toLowerCase();
    let hit = null;
    for (const [tok] of f.options) {
      if ((f.match === 'exact' && t === tok.toLowerCase()) || (f.match !== 'exact' && t && t[0] === tok.toLowerCase()[0])) { hit = tok; break; }
    }
    ctl = `<select class="ch" ${data}>` + (hit === null ? `<option selected value="">(unrecognised: ${esc(f.value)})</option>` : '') +
      f.options.map(([tok, lab]) => `<option value="${esc(tok)}" ${tok === hit ? 'selected' : ''}>${esc(lab)}</option>`).join('') + '</select>';
  }
  return `<tr><td class="lab">${esc(f.label)}</td><td>${ctl}</td><td class="unit">${esc(f.unit || '')}</td><td class="mark">${mark}</td>` +
    (f.note ? `<td class="note">${esc(f.note)}</td>` : '') + '</tr>';
}

function sectionHtml(s, w) {
  const rows = s.fields.map(f => fieldRow(f, w)).join('');
  const btn = (s.buttons || []).map(([a, l]) => `<button data-act="${a}">${esc(l)}</button>`).join(' ');
  return `<fieldset><legend>${esc(s.title)}</legend>${s.desc ? `<div class="desc">${esc(s.desc)}</div>` : ''}<table class="form">${rows}</table>${btn ? '<div>' + btn + '</div>' : ''}</fieldset>`;
}

function renderTab() {
  const pane = $('#tabpane');
  const act = document.activeElement && document.activeElement.dataset ? document.activeElement.dataset.fid : null;
  const top = pane.scrollTop;
  const w = worst();
  if (curTab === 'Summary') {
    const rows = S.issues.map((i, k) => `<tr data-k="${k}"><td class="lv-${i.level}">${LEVEL_NAME[i.level]}</td><td>${esc(i.tab || '')}</td><td>${esc(i.msg)}</td></tr>`).join('')
      || '<tr><td colspan="3">No problems found.</td></tr>';
    pane.innerHTML = `<fieldset><legend>Derived quantities (what the code will compute at start-up)</legend><pre class="derived">${esc(S.derived.join('\n'))}</pre></fieldset>
      <fieldset><legend>Checks — click a line to go to the field</legend><table class="issues">${rows}</table>
      <div class="desc">"BLOCKING" prevents saving; "problem" means the Fortran code would stop or misbehave (you may still save).</div></fieldset>`;
    pane.querySelectorAll('table.issues tr[data-k]').forEach(tr => tr.onclick = () => jump(S.issues[+tr.dataset.k]));
  } else {
    const t = S.tabs.find(x => x.name === curTab);
    let html = t.sections.map(s => sectionHtml(s, w)).join('');
    if (t.geometry) html = `<div class="geo"><div>${html}${segHtml(w)}</div><div>${previewHtml()}</div></div>`;
    pane.innerHTML = html;
    if (t.geometry) { bindSegs(); loadPreview(); }
  }
  bindFields(pane);
  pane.scrollTop = top;
  if (act) { const el = pane.querySelector(`[data-fid="${CSS.escape(act)}"]:not(.sg-c)`) || pane.querySelector(`[data-fid="${CSS.escape(act)}"]`); if (el) el.focus(); }
}

function bindFields(root) {
  root.querySelectorAll('[data-key]').forEach(el => {
    const help = () => showHelp(el.dataset.help, el.dataset.ikey);
    el.addEventListener('focus', help);
    el.addEventListener('mouseenter', help);
    if (el.matches('input[type=text]:not(.sg-v)')) el.addEventListener('change', () => send(el.dataset.key, el.value));
    else if (el.matches('.sg-c, .sg-v')) el.addEventListener('change', () => {
      const row = el.closest('td');
      const c = row.querySelector('.sg-c'), v = row.querySelector('.sg-v');
      const neg = c.dataset.on === 'positive' ? !c.checked : c.checked;
      send(el.dataset.key, (neg ? '-' : '') + v.value.trim().replace(/^[+-]/, ''));
    });
    else if (el.matches('.yn')) el.addEventListener('change', () => send(el.dataset.key, el.dataset.bool === '1' ? (el.checked ? '1' : '0') : (el.checked ? 'yes' : 'no')));
    else if (el.matches('.ch')) el.addEventListener('change', () => { if (el.value) send(el.dataset.key, el.value); });
  });
  root.querySelectorAll('[data-act]').forEach(b => b.onclick = () => action(b.dataset.act));
}

async function send(key, value) {
  try { setState(await api('/api/set', { key, value, sel: selSeg })); }
  catch (e) { status(e.message, true); }
}

async function action(act, arg) {
  try { setState(await api('/api/action', { action: act, arg })); }
  catch (e) { status(e.message, true); }
}

function showHelp(hk, ikey) {
  const h = $('#help');
  if (!hk && !ikey) {
    h.innerHTML = '<div class="t">LePIC+ input editor</div>Click or hover over a field to see what it does in the code, its unit and any problem found. ' +
      'Nothing is written until you press Save (previous files are kept in gui_backups/).';
    return;
  }
  let html = '';
  if (HELP[hk]) html += `<div><span class="t">${esc(HELP[hk][0])}</span><span class="code">code: ${esc(HELP[hk][1])}</span></div>`;
  for (const i of (S ? S.issues : [])) if (i.key === ikey) html += `<div class="lv-${i.level}">${MARK[i.level]} ${esc(i.msg)}</div>`;
  if (HELP[hk]) html += `<div>${esc(HELP[hk][2]).replace(/\n/g, '<br>')}</div>`;
  h.innerHTML = html;
}

function jump(i) {
  let tab = i.tab && tabNames().includes(i.tab) ? i.tab : null, key = i.key;
  const m = key && key.match(/^g\.seg(\d+)\.(\d)$/);
  if (m) { tab = 'Geometry'; selSeg = +m[1]; }
  if (!tab && key) for (const t of S.tabs) for (const s of t.sections) for (const f of s.fields) if (f.ikey === key && !tab) tab = t.name;
  if (tab) { curTab = tab; renderTabs(); renderTab(); }
  const el = key && (document.querySelector(`#tabpane [data-ikey="${CSS.escape(m ? key : key)}"]`) || document.querySelector(`#tabpane [data-fid="${CSS.escape(key)}"]`));
  if (el) { el.focus(); el.scrollIntoView({ block: 'center' }); }
  showHelp(el ? el.dataset.help : null, key);
}

// ---- geometry
function segHtml(w) {
  const rows = S.segments.map((s, k) => {
    let lv = null;
    for (const i of S.issues) if (i.key && i.key.startsWith('g.seg' + k + '.') && (!lv || ORDER[i.level] < ORDER[lv])) lv = i.level;
    const cls = (k === selSeg ? 'sel ' : '') + (lv === 'error' || lv === 'fatal' ? 'e' : lv === 'warning' ? 'w' : '');
    return `<tr data-k="${k}" class="${cls}"><td>${k + 1}</td><td>${esc(s.label)}</td><td>${esc(s.kind)}</td><td>${esc(s.x)}</td><td>${esc(s.y)}</td><td>${esc(s.z)}</td><td>${esc(s.flags)}</td></tr>`;
  }).join('');
  let form = '';
  const s = S.segments[selSeg];
  if (s) {
    const k = selSeg, h = 'g.seg';
    const mk = (i) => { const lv = w['g.seg' + k + '.' + i]; return lv ? `<span class="lv-${lv}">${MARK[lv]}</span>` : ''; };
    const cls = (i) => { const lv = w['g.seg' + k + '.' + i]; return lv === 'error' || lv === 'fatal' ? 'err' : lv === 'warning' ? 'warn' : ''; };
    const dat = (key, i) => `data-key="g.seg${k}.${key}" data-fid="g.seg${k}.${key}" data-help="${h}" data-ikey="g.seg${k}.${i}"`;
    const kinds = [['+', 'Wall with an opening (plasma window)'], ['-', 'Solid block / plate'], ['0', 'Periodic faces']];
    form = `<table class="form"><tr><td class="lab">Kind</td><td><select class="ch" ${dat('kind', 6)}>${kinds.map(([t, l]) => `<option value="${t}" ${t === s.kindtok ? 'selected' : ''}>${l}</option>`).join('')}</select></td><td></td><td class="mark">${mk(6)}</td></tr>
      <tr><td class="lab">Label</td><td><input type="text" class="${cls(6)}" size="6" ${dat('label', 6)} value="${esc(s.label)}"></td><td class="unit">(1 … number of wall labels)</td><td></td></tr>` +
      ['x from', 'y from', 'z from', 'x to', 'y to', 'z to'].map((l, i) => `<tr><td class="lab">${l}</td><td><input type="text" class="${cls(i)}" size="10" ${dat(i, i)} value="${esc(s.coords[i])}"></td><td class="unit">cm</td><td class="mark">${mk(i)}</td></tr>`).join('') +
      [[3, 'Dielectric surface in an x = const plane (x-to stored negative)'], [4, 'Dielectric surface in a y = const plane (y-to stored negative)'],
       [5, 'Periodic in z (z-to stored negative)'], [1, 'Toggle round / rectangular for this segment (y-from stored negative)']]
        .map(([i, t]) => `<tr><td colspan="4"><label><input type="checkbox" class="yn" data-bool="1" ${dat('flag' + i, i)} ${s.flagvals[i] ? 'checked' : ''}> ${t}</label></td></tr>`).join('') + '</table>';
  }
  return `<fieldset><legend>Wall segments (drawn in this order — later ones overwrite earlier ones)</legend>
    <table class="segs"><thead><tr><th>#</th><th>Label</th><th>Kind</th><th>x range</th><th>y range</th><th>z range</th><th>Options</th></tr></thead><tbody>${rows}</tbody></table>
    <div>${['seg_add:Add', 'seg_dup:Duplicate', 'seg_del:Delete', 'seg_up:▲ Up', 'seg_down:▼ Down'].map(x => { const [a, l] = x.split(':'); return `<button data-seg="${a}">${l}</button>`; }).join(' ')}</div></fieldset>
    <fieldset><legend>Selected segment ${S.segments.length ? selSeg + 1 : ''}</legend>${form || '<span class="muted">No segment.</span>'}</fieldset>`;
}

function bindSegs() {
  document.querySelectorAll('#tabpane table.segs tbody tr').forEach(tr => tr.onclick = () => { selSeg = +tr.dataset.k; renderTab(); });
  document.querySelectorAll('#tabpane [data-seg]').forEach(b => b.onclick = () => action(b.dataset.seg, selSeg));
}

function previewHtml() {
  const sl = (k, t) => `<label>${t} <input type="range" min="0" max="1" step="0.005" value="${geoFrac[k]}" data-pl="${k}"> <span id="gp-${k}">${(geoFrac[k] * 100).toFixed(0)}%</span></label>`;
  return `<fieldset><legend>Cross-sections of the walls</legend><div class="ctl">${sl('xy', 'xy plane at z')}${sl('xz', 'xz plane at y')}${sl('yz', 'yz plane at x')}
    <button id="gp-reset">Middle</button></div><img id="gp-img" style="max-width:100%"><div class="muted" id="gp-note"></div></fieldset>`;
}

let geoTimer = null;
function loadPreview() {
  document.querySelectorAll('#tabpane [data-pl]').forEach(r => r.oninput = () => {
    geoFrac[r.dataset.pl] = +r.value; $('#gp-' + r.dataset.pl).textContent = (r.value * 100).toFixed(0) + '%';
    clearTimeout(geoTimer); geoTimer = setTimeout(fetchPreview, 250);
  });
  $('#gp-reset').onclick = () => { geoFrac = { xy: 0.5, xz: 0.5, yz: 0.5 }; renderTab(); };
  fetchPreview();
}

let geoSeq = 0;
async function fetchPreview() {
  const seq = ++geoSeq;
  try {
    const r = await api(`/api/geometry?xy=${geoFrac.xy}&xz=${geoFrac.xz}&yz=${geoFrac.yz}`);
    if (seq !== geoSeq) return;
    const img = $('#gp-img'); if (img) img.src = 'data:image/png;base64,' + r.png;
    const n = $('#gp-note'); if (n) n.textContent = (r.meta.notes || []).join(' ');
  } catch (e) { const n = $('#gp-note'); if (n) n.textContent = e.message; }
}

// ---- toolbar
$('#btn-save').onclick = async () => {
  if (!S) return;
  const errs = S.issues.filter(i => i.level === 'error');
  if (errs.length && !confirm('The simulation would stop or misbehave with these settings:\n\n' + errs.slice(0, 8).map(i => '- ' + i.msg).join('\n') + '\n\nSave anyway?')) return;
  try { setState(await api('/api/save', {})); } catch (e) { status(e.message, true); }
};
$('#btn-reload').onclick = async () => {
  if (S && S.modified && !confirm('Discard all unsaved changes?')) return;
  try { setState(await api('/api/reload', {})); status('Reloaded from disk.'); } catch (e) { status(e.message, true); }
};
$('#btn-activate').onclick = async () => {
  if (!S) return;
  if (S.modified) { alert('Save the case first.'); return; }
  if (!confirm('LePIC+ only reads the files directly in input_dir/. Copy this case there?\n\nThe current input_dir/ files are kept in input_dir/gui_backups/.')) return;
  try { setState(await api('/api/activate', {})); } catch (e) { status(e.message, true); }
};
$('#btn-diff').onclick = async () => {
  const r = await api('/api/diff');
  $('#modal-body').innerHTML = '<h3>Changes that will be written</h3>' + (r.lines.length ? '<pre>' + r.lines.map(l =>
    `<span class="${/^(---|\+\+\+|@@)/.test(l) ? 'muted' : l[0] === '+' ? 'add' : l[0] === '-' ? 'del' : ''}">${esc(l)}</span>`).join('\n') + '</pre>'
    : '<p>No changes: the files on disk are identical to what the editor would write.</p>');
  $('#modal').classList.add('on');
};
$('#modal-close').onclick = () => $('#modal').classList.remove('on');
window.addEventListener('beforeunload', e => { if (S && S.modified) { e.preventDefault(); e.returnValue = ''; } });

// =====================================================================
//  hover readout shared by the two viewers
// =====================================================================
function attachReadout(img, out, getMeta) {
  img.addEventListener('mousemove', ev => {
    const m = getMeta(); if (!m) return;
    const r = img.getBoundingClientRect();
    const fx = (ev.clientX - r.left) / r.width, fy = (ev.clientY - r.top) / r.height;
    const [ax, ay, aw, ah] = m.axes;
    if (fx < ax || fx > ax + aw || fy < ay || fy > ay + ah) { out.textContent = ''; return; }
    const x = m.xlim[0] + (fx - ax) / aw * (m.xlim[1] - m.xlim[0]);
    const tt = 1 - (fy - ay) / ah;
    const y = m.ylog ? Math.pow(10, Math.log10(m.ylim[0]) + tt * (Math.log10(m.ylim[1]) - Math.log10(m.ylim[0]))) : m.ylim[0] + tt * (m.ylim[1] - m.ylim[0]);
    if (m.kind === 'map') {
      const ix = Math.round((x - m.x[0]) / (m.x[1] - m.x[0]) * m.nx), iy = Math.round((y - m.y[0]) / (m.y[1] - m.y[0]) * m.ny);
      if (ix < 0 || iy < 0 || ix > m.nx || iy > m.ny) { out.textContent = ''; return; }
      const row = m.z[Math.round(iy / m.step)], v = row ? row[Math.round(ix / m.step)] : undefined;
      out.textContent = `${m.hn} = ${x.toPrecision(4)} ${m.unit}  (node ${ix + 1})    ${m.vn} = ${y.toPrecision(4)} ${m.unit}  (node ${iy + 1})    value = ${v === undefined ? '?' : Number(v).toPrecision(5)} ${m.zunit || ''}`;
    } else {
      const xs = m.x; let k = 0, best = Infinity;
      for (let i = 0; i < xs.length; i++) { const d = Math.abs(xs[i] - x); if (d < best) { best = d; k = i; } }
      out.textContent = `${m.xn} = ${Number(xs[k]).toPrecision(5)} ${m.xu}    ` + m.ys.map((s, j) => (m.ys.length > 1 ? `[${j + 1}] ` : '') + Number(s[k]).toPrecision(5)).join('   ');
    }
  });
  img.addEventListener('mouseleave', () => out.textContent = '');
}

// =====================================================================
//  OUTPUT 2D
// =====================================================================
const O = { scan: null, key: null, meta: null, sig: '' };
$('#o-cmap').innerHTML = CMAPS.map(c => `<option>${c}</option>`).join('');
$('#t-cmap').innerHTML = CMAPS.map(c => `<option>${c}</option>`).join('');

async function oOpen(keepSnap) {
  const snap = keepSnap ? $('#o-snap').value : 'latest';
  try {
    O.scan = await api(`/api/out/scan?dir=${encodeURIComponent($('#o-dir').value)}&snap=${encodeURIComponent(snap)}`);
  } catch (e) { $('#o-status').textContent = e.message; $('#o-tree').innerHTML = ''; return; }
  const sel = $('#o-snap');
  sel.innerHTML = ['latest'].concat(O.scan.snapshots).map(s => `<option value="${s}">${s === 'latest' || s === 'initial' ? s : 'it ' + s}</option>`).join('');
  sel.value = O.scan.snapshots.includes(snap) || snap === 'latest' ? snap : 'latest';
  sel.disabled = !O.scan.snapshots.length;
  $('#o-planes').style.display = O.scan.maps ? '' : 'none';
  const groups = {};
  for (const e of O.scan.entries) (groups[e.group] = groups[e.group] || []).push(e);
  $('#o-tree').innerHTML = Object.entries(groups).map(([g, es]) => `<div class="g">${esc(g)}</div>` +
    es.map(e => `<div class="it ${e.key === O.key ? 'on' : ''}" data-k="${esc(e.key)}">${esc(e.label)}</div>`).join('')).join('')
    || '<span class="muted">No LePIC+ output found here.</span>';
  document.querySelectorAll('#o-tree .it').forEach(d => d.onclick = () => { O.key = d.dataset.k; document.querySelectorAll('#o-tree .it').forEach(x => x.classList.toggle('on', x === d)); oCtl(); oPlot(); });
  const g = O.scan.geometry;
  $('#o-status').textContent = O.scan.dir + (g.path ? '   — box from ' + g.path + ' (must be the case that wrote these files)' : '   — geometry.inp not found: axes in nodes');
  if (!O.key || !O.scan.entries.find(e => e.key === O.key)) O.key = O.scan.entries.length ? O.scan.entries[0].key : null;
  oCtl();
  if (O.key) oPlot();
}

function oEntry() { return O.scan && O.scan.entries.find(e => e.key === O.key); }

function oCtl() {
  const e = oEntry(); if (!e) return;
  const isMap = e.kind === 'field' || e.kind === 'flux';
  $('#o-mode').parentElement.style.display = isMap ? '' : 'none';
  $('#o-cmap').parentElement.style.display = isMap ? '' : 'none';
  $('#o-vmin').parentElement.style.display = isMap ? '' : 'none';
  $('#o-vmax').parentElement.style.display = isMap ? '' : 'none';
  $('#o-stream').parentElement.style.display = isMap && /jvec|mag:B$/.test(e.key) ? '' : 'none';
  $('#o-probe-box').style.display = e.kind === 'probe_ne' ? '' : 'none';
  const cols = $('#o-cols');
  if (e.kind === 'iw') {
    const first = O.scan.geometry.nmn ? 0 : 1;
    if (cols.dataset.k !== e.key) {
      cols.innerHTML = 'Wall labels: ' + Array.from({ length: e.ncols }, (_, i) => `<label><input type="checkbox" value="${i + 1}" checked> ${i + first}</label>`).join(' ');
      cols.dataset.k = e.key;
      cols.querySelectorAll('input').forEach(c => c.onchange = oPlot);
    }
    cols.style.display = '';
  } else cols.style.display = 'none';
  if (isMap && e.planes.length && !e.planes.includes(oPlane())) document.querySelector(`input[name=o-plane][value=${e.planes[0]}]`).checked = true;
}

function oPlane() { return document.querySelector('input[name=o-plane]:checked').value; }

async function oPlot(auto) {
  const e = oEntry(); if (!e) return;
  const box = $('#o-img').parentElement.getBoundingClientRect();
  const w = Math.max(6, Math.min(16, (box.width - 20) / 90)), h = Math.max(4, Math.min(9, (window.innerHeight - 230) / 90));
  const q = new URLSearchParams({
    dir: $('#o-dir').value, snap: $('#o-snap').value, key: e.key, plane: oPlane(), mode: $('#o-mode').value,
    cmap: $('#o-cmap').value, log: $('#o-log').checked ? 1 : 0, vmin: $('#o-vmin').value, vmax: $('#o-vmax').value,
    stream: $('#o-stream').checked ? 1 : 0, probe: $('#o-probe').value, w: w.toFixed(1), h: h.toFixed(1),
    cols: Array.from(document.querySelectorAll('#o-cols input:checked')).map(c => c.value).join(','),
  });
  const seq = O.seq = (O.seq || 0) + 1;            // only the newest request may draw
  try {
    const r = await api('/api/out/plot?' + q);
    if (seq !== O.seq) return;
    $('#o-img').src = 'data:image/png;base64,' + r.png;
    O.meta = r.meta;
    if (!auto) $('#o-read').textContent = '';
  } catch (err) { if (seq === O.seq) $('#o-read').textContent = err.message; }
}

['o-mode', 'o-cmap', 'o-log', 'o-stream', 'o-probe'].forEach(id => $('#' + id).addEventListener('change', () => oPlot()));
['o-vmin', 'o-vmax'].forEach(id => $('#' + id).addEventListener('change', () => oPlot()));
document.querySelectorAll('input[name=o-plane]').forEach(r => r.onchange = () => oPlot());
$('#o-snap').onchange = () => oOpen(true);
$('#o-open').onclick = () => oOpen();
$('#o-reload').onclick = () => oOpen(true);
attachReadout($('#o-img'), $('#o-read'), () => O.meta);
setInterval(() => {
  if (!$('#page-out2d').classList.contains('on') || !$('#o-auto').checked || !O.scan) return;
  const now = Date.now();
  if (now - (O.last || 0) < Math.max(5, +$('#o-int').value) * 1000) return;
  O.last = now; oOpen(true);
}, 1000);

// =====================================================================
//  OUTPUT 3D
// =====================================================================
const T = { loaded: false, meta: null };
async function tPlot(reset) {
  const axis = document.querySelector('input[name=t-axis]:checked').value;
  const q = new URLSearchParams({
    file: $('#t-file').value, qty: $('#t-qty').value, axis, index: reset ? '' : $('#t-idx').value,
    cmap: $('#t-cmap').value, log: $('#t-log').checked ? 1 : 0, w: 10, h: 5.6,
  });
  const seq = T.seq = (T.seq || 0) + 1;
  $('#t-status').textContent = 'Reading ' + $('#t-file').value + ' …';
  try {
    const r = await api('/api/3d/plot?' + q);
    if (seq !== T.seq) return;
    $('#t-img').src = 'data:image/png;base64,' + r.png;
    T.meta = r.meta; T.loaded = true;
    const n = r.info.n['xyz'.indexOf(axis)];
    const sl = $('#t-idx'); sl.max = n; sl.value = r.index;
    $('#t-status').textContent = `${r.info.file}: ${r.info.n.join(' x ')} cells, density of ${r.info.species} (ptype ${r.info.ptype})` +
      (r.info.geometry ? `, box from ${r.info.geometry} (must be the case that wrote the file)` : ', geometry.inp not found: axes in nodes');
  } catch (e) { if (seq === T.seq) $('#t-status').textContent = e.message; }
}
$('#t-open').onclick = () => tPlot(true);
['t-qty', 't-cmap', 't-log'].forEach(id => $('#' + id).addEventListener('change', () => tPlot()));
document.querySelectorAll('input[name=t-axis]').forEach(r => r.onchange = () => tPlot(true));
let tTimer = null;
$('#t-idx').addEventListener('input', () => { clearTimeout(tTimer); tTimer = setTimeout(() => tPlot(), 200); });
attachReadout($('#t-img'), $('#t-read'), () => T.meta);

initInputs().catch(e => status(e.message, true));
