"""`stew serve`: a read-only web page over the catalogue.

It opens the catalogue read-only, answers a few JSON endpoints and serves one
self-contained HTML page (no external scripts, so it works offline). It never
serves dataset files themselves."""

from __future__ import annotations

import json
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

from . import __version__
from .db import Catalog, CatalogError
from .util import Presence, presence

LIST_FIELDS = ("id", "kind", "name", "path", "platform", "format", "organ", "species", "study", "source", "description",
               "owner", "steward", "size_bytes", "n_files", "n_cells", "n_features", "modalities", "tags",
               "added_at", "updated_at")


class Store:
    """One read-only connection per thread."""

    def __init__(self, path: Path):
        self.path = path
        self.local = threading.local()
        Catalog.open_readonly(path)  # fail early if the catalogue can't be read

    @property
    def cat(self) -> Catalog:
        if not hasattr(self.local, "cat"):
            self.local.cat = Catalog.open_readonly(self.path)
        return self.local.cat


def api(store: Store, route: str):
    cat = store.cat
    if route == "/api/datasets":
        return {"datasets": [{k: d.get(k) for k in LIST_FIELDS} for d in cat.list()],
                "relations": cat.all_relations(), "summary": cat.summary(), "catalog": str(store.path),
                "version": __version__}
    if route.startswith("/api/datasets/"):
        key = route.rsplit("/", 1)[1]
        if not key.isdigit():
            return None
        d = cat.by_id(int(key))
        if not d:
            return None
        up, down = cat.relations(d)
        d["relations"] = {"from": [{"kind": k, "id": x["id"], "name": x["name"]} for k, x in up],
                          "to": [{"kind": k, "id": x["id"], "name": x["name"]} for k, x in down]}
        d["links"] = len(cat.links(dataset_id=d["id"]))
        d["exists"] = presence(Path(d["path"])) is not Presence.MISSING
        files = cat.files(d["id"])
        d["checksummed_files"] = len(files)
        return d
    if route == "/api/downloads":
        return {"downloads": cat.downloads(limit=100)}
    return None


class Handler(BaseHTTPRequestHandler):
    store: Store = None
    server_version = f"stew/{__version__}"

    def log_message(self, fmt, *args):
        sys.stderr.write(f"{self.address_string()} {fmt % args}\n")

    def send(self, code: int, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def do_GET(self):
        route = urlsplit(self.path).path
        if route in ("/", "/index.html"):
            return self.send(200, PAGE.encode(), "text/html; charset=utf-8")
        if route.startswith("/api/"):
            try:
                data = api(self.store, route)
            except CatalogError as e:
                return self.send(503, json.dumps({"error": str(e)}).encode(), "application/json")
            if data is None:
                return self.send(404, b'{"error": "not found"}', "application/json")
            return self.send(200, json.dumps(data, default=str).encode(), "application/json")
        self.send(404, b"not found", "text/plain")

    do_HEAD = do_GET


def serve(catalog: Path, host: str, port: int) -> None:
    Handler.store = Store(catalog)
    httpd = ThreadingHTTPServer((host, port), Handler)
    print(f"Serving {catalog} read-only at http://{host}:{port}/  (Ctrl-C to stop)", flush=True)
    if host in ("127.0.0.1", "localhost"):
        print(f"From your own computer: ssh -N -L {port}:localhost:{port} <you>@<this server>, "
              f"then open http://localhost:{port}/", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>stew · dataset catalogue</title>
<style>
:root { --bg:#fbfbfa; --fg:#1d1d1b; --muted:#6b6b66; --line:#e3e2dd; --accent:#2f6f8f; --chip:#eef3f5; --warn:#a4452c; }
@media (prefers-color-scheme: dark) { :root { --bg:#161615; --fg:#e8e6e1; --muted:#9a988f; --line:#2e2d2a;
  --accent:#7fb4cc; --chip:#22303a; --warn:#e08a6e; } }
* { box-sizing:border-box; }
body { margin:0; font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif; background:var(--bg); color:var(--fg); }
header { padding:16px 24px 8px; display:flex; gap:16px; align-items:baseline; flex-wrap:wrap; }
header h1 { font-size:18px; margin:0; font-weight:600; }
header .meta { color:var(--muted); font-size:12px; }
#filters { padding:8px 24px; display:flex; gap:8px; flex-wrap:wrap; align-items:center; border-bottom:1px solid var(--line); }
input, select { font:inherit; padding:5px 8px; border:1px solid var(--line); border-radius:6px; background:var(--bg); color:var(--fg); }
#q { min-width:260px; flex:1 1 260px; }
button { font:inherit; padding:5px 10px; border:1px solid var(--line); border-radius:6px; background:var(--chip); color:var(--fg); cursor:pointer; }
#stats { padding:8px 24px; color:var(--muted); font-size:12px; }
main { display:grid; grid-template-columns:minmax(0,1fr); }
main.open { grid-template-columns:minmax(0,1fr) minmax(320px,38%); }
.tablewrap { overflow:auto; max-height:calc(100vh - 150px); }
table { border-collapse:collapse; width:100%; }
th, td { text-align:left; padding:6px 10px; border-bottom:1px solid var(--line); white-space:nowrap; }
th { position:sticky; top:0; background:var(--bg); font-weight:600; font-size:12px; color:var(--muted); cursor:pointer; user-select:none; }
th.sorted::after { content:" ↓"; } th.sorted.asc::after { content:" ↑"; }
td.num { text-align:right; font-variant-numeric:tabular-nums; }
tr:hover td { background:var(--chip); cursor:pointer; }
tr.sel td { background:var(--chip); }
.chip { display:inline-block; padding:0 6px; margin:0 2px 2px 0; border-radius:10px; background:var(--chip); font-size:12px; }
aside { border-left:1px solid var(--line); padding:16px 20px; overflow:auto; max-height:calc(100vh - 150px); }
aside h2 { margin:0 0 4px; font-size:16px; word-break:break-all; }
aside dl { display:grid; grid-template-columns:110px 1fr; gap:4px 10px; margin:12px 0; }
aside dt { color:var(--muted); } aside dd { margin:0; word-break:break-all; }
code, .mono { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:12px; }
pre { background:var(--chip); padding:8px 10px; border-radius:6px; white-space:pre-wrap; word-break:break-all; }
.muted { color:var(--muted); } .warn { color:var(--warn); }
a { color:var(--accent); }
</style></head>
<body>
<header><h1>Dataset catalogue</h1><span class="meta" id="meta"></span></header>
<div id="filters">
  <input id="q" type="search" placeholder="Search name, description, study, source, path… (all words must match)">
  <select id="f-kind"><option value="dataset">datasets</option><option value="reference">references</option><option value="">everything</option></select>
  <select id="f-organ"><option value="">any organ</option></select>
  <select id="f-species"><option value="">any species</option></select>
  <select id="f-modality"><option value="">any modality</option></select>
  <select id="f-platform"><option value="">any platform</option></select>
  <button id="reset">Clear</button>
  <button id="csv">Download CSV</button>
</div>
<div id="stats"></div>
<main id="main">
  <div class="tablewrap"><table><thead><tr id="head"></tr></thead><tbody id="rows"></tbody></table></div>
  <aside id="detail" hidden></aside>
</main>
<script>
"use strict";
const COLS = [["id","ID",1],["name","Name"],["platform","Platform"],["modalities","Modalities"],["organ","Organ"],
  ["species","Species"],["n_cells","Cells",1],["size_bytes","Size",1],["source","Source"],["owner","Owner"]];
let data = [], rels = [], sortKey = "id", sortAsc = true, selected = null;
const $ = id => document.getElementById(id);
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
function size(n) { if (n == null) return "–"; const u = ["B","KiB","MiB","GiB","TiB","PiB"]; let i = 0, v = n;
  while (v >= 1024 && i < u.length - 1) { v /= 1024; i++; } return i ? v.toFixed(1) + " " + u[i] : n + " B"; }
const num = n => n == null ? "–" : Number(n).toLocaleString();
function cell(d, k) {
  if (k === "modalities") return d.modalities.map(m => `<span class="chip">${esc(m)}</span>`).join("") || "–";
  if (k === "size_bytes") return size(d.size_bytes);
  if (k === "n_cells") return num(d.n_cells);
  return esc(d[k] ?? "–");
}
function haystack(d) { return [d.name,d.description,d.study,d.source,d.path,d.organ,d.species,d.platform,d.format,
  d.owner,d.steward,...d.modalities,...d.tags].join(" ").toLowerCase(); }
function filtered() {
  const words = $("q").value.toLowerCase().split(/\s+/).filter(Boolean);
  const f = k => $("f-" + k).value;
  return data.filter(d => (!f("kind") || d.kind === f("kind")) && (!f("organ") || d.organ === f("organ")) && (!f("species") || d.species === f("species"))
    && (!f("platform") || d.platform === f("platform")) && (!f("modality") || d.modalities.includes(f("modality")))
    && words.every(w => d._hay.includes(w)));
}
function render() {
  const rows = filtered().sort((a, b) => {
    let x = a[sortKey], y = b[sortKey];
    if (Array.isArray(x)) { x = x.join(","); y = y.join(","); }
    if (x == null) return 1; if (y == null) return -1;
    return (x < y ? -1 : x > y ? 1 : 0) * (sortAsc ? 1 : -1);
  });
  $("head").innerHTML = COLS.map(([k, label]) =>
    `<th data-k="${k}" class="${k === sortKey ? "sorted" + (sortAsc ? " asc" : "") : ""}">${label}</th>`).join("");
  $("rows").innerHTML = rows.map(d => `<tr data-id="${d.id}" class="${d.id === selected ? "sel" : ""}">` +
    COLS.map(([k, , n]) => `<td class="${n ? "num" : ""}">${cell(d, k)}</td>`).join("") + "</tr>").join("");
  const total = rows.reduce((s, d) => s + (d.size_bytes || 0), 0);
  $("stats").textContent = `${rows.length} of ${data.length} datasets · ${size(total)}`;
  return rows;
}
function options(id, values) {
  const sel = $(id), cur = sel.value;
  [...new Set(values.filter(Boolean))].sort().forEach(v => { const o = document.createElement("option"); o.value = o.textContent = v; sel.appendChild(o); });
  sel.value = cur;
}
async function showDetail(id) {
  selected = id; render();
  const r = await fetch("api/datasets/" + id); if (!r.ok) return;
  const d = await r.json();
  const rel = (list, word) => list.map(x => `<div>${esc(x.kind)} ${word} <a href="#${x.id}" data-id="${x.id}">#${x.id} ${esc(x.name)}</a></div>`).join("");
  const det = d.details || {};
  const extra = Object.entries(det).filter(([k]) => k !== "obs_columns").map(([k, v]) =>
    `<dt>${esc(k.replace(/_/g, " "))}</dt><dd>${esc(typeof v === "object" ? JSON.stringify(v) : v)}</dd>`).join("");
  $("detail").innerHTML = `<h2>#${d.id} ${esc(d.name)}</h2>
    <div class="muted">${esc(d.description || "")}</div>
    ${d.exists ? "" : '<p class="warn">The path no longer exists.</p>'}
    <pre>stew get ${esc(d.name)} --to .</pre>
    <dl>
      <dt>path</dt><dd class="mono">${esc(d.path)}</dd>
      <dt>source</dt><dd>${esc(d.source || "–")}</dd>
      <dt>platform</dt><dd>${esc(d.platform || "–")}</dd><dt>format</dt><dd>${esc(d.format || "–")}</dd>
      <dt>modalities</dt><dd>${d.modalities.map(m => `<span class="chip">${esc(m)}</span>`).join("") || "–"}</dd>
      <dt>organ</dt><dd>${esc(d.organ || "–")}</dd><dt>species</dt><dd>${esc(d.species || "–")}</dd>
      <dt>cells</dt><dd>${num(d.n_cells)}</dd><dt>features</dt><dd>${num(d.n_features)}</dd>
      <dt>study</dt><dd>${esc(d.study || "–")}</dd>
      <dt>tags</dt><dd>${d.tags.map(t => `<span class="chip">${esc(t)}</span>`).join("") || "–"}</dd>
      <dt>size</dt><dd>${size(d.size_bytes)} in ${num(d.n_files)} file(s)</dd>
      <dt>owner</dt><dd>${esc(d.owner || "–")}</dd><dt>steward</dt><dd>${esc(d.steward || "–")}</dd>
      <dt>added</dt><dd>${esc(d.added_at)} by ${esc(d.added_by)}</dd>
      <dt>checksum</dt><dd class="mono">${esc(d.fingerprint || "not computed")}</dd>
      <dt>links out</dt><dd>${d.links}</dd>
      ${extra}
    </dl>
    ${d.relations.from.length ? "<h3>Comes from</h3>" + rel(d.relations.from, "") : ""}
    ${d.relations.to.length ? "<h3>Used by</h3>" + rel(d.relations.to, "") : ""}`;
  $("detail").hidden = false; $("main").classList.add("open");
  location.hash = id;
}
function csv(rows) {
  const keys = ["id","name","path","platform","format","modalities","organ","species","study","source","n_cells",
    "n_features","size_bytes","owner","steward","tags","description"];
  const q = v => { v = Array.isArray(v) ? v.join(";") : (v ?? ""); v = String(v); return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; };
  const text = [keys.join(","), ...rows.map(d => keys.map(k => q(d[k])).join(","))].join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], {type: "text/csv"})); a.download = "datasets.csv"; a.click();
}
async function load() {
  const r = await fetch("api/datasets");
  const j = await r.json();
  if (!r.ok) { $("stats").textContent = j.error || "could not load the catalogue"; return; }
  data = j.datasets.map(d => ({...d, _hay: haystack(d)})); rels = j.relations;
  $("meta").textContent = `${j.catalog} · stew ${j.version} · read-only`;
  options("f-organ", data.map(d => d.organ)); options("f-species", data.map(d => d.species));
  options("f-platform", data.map(d => d.platform)); options("f-modality", data.flatMap(d => d.modalities));
  render();
  const h = parseInt(location.hash.slice(1)); if (h) showDetail(h);
}
$("head").addEventListener("click", e => { const k = e.target.dataset.k; if (!k) return;
  sortAsc = k === sortKey ? !sortAsc : true; sortKey = k; render(); });
$("rows").addEventListener("click", e => { const tr = e.target.closest("tr"); if (tr) showDetail(+tr.dataset.id); });
$("detail").addEventListener("click", e => { const id = e.target.dataset && e.target.dataset.id; if (id) { e.preventDefault(); showDetail(+id); } });
["q","f-kind","f-organ","f-species","f-modality","f-platform"].forEach(id => $(id).addEventListener("input", render));
$("reset").addEventListener("click", () => { ["q","f-organ","f-species","f-modality","f-platform"].forEach(id => $(id).value = ""); $("f-kind").value = "dataset"; render(); });
$("csv").addEventListener("click", () => csv(render()));
load();
</script></body></html>
"""
