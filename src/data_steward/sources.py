"""Where datasets come from: GEO, SRA/ENA, CELLxGENE, ArrayExpress/BioStudies,
figshare, Zenodo, the Human Cell Atlas, and plain URLs.

`parse` turns what someone typed into a canonical key without touching the
network, so the "already downloaded?" check can run first. `resolve` then asks
the remote side what files there are."""

from __future__ import annotations

import html
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

from . import util, vocab

USER_AGENT = "data-steward (stew; lab dataset catalogue)"
TIMEOUT = 60
UUID = r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
RAW_KINDS = ("geo", "hca")  # sources that have both processed files and raw reads

KINDS_HELP = ("a GEO accession (GSE…/GSM…), an SRA/ENA accession (SRP…, SRR…, PRJNA…, ERP…), "
              "an ArrayExpress/BioStudies accession (E-MTAB-…), cellxgene:<collection id>, "
              "hca:<project id>, figshare:<id>, zenodo:<id>, or an http(s)/ftp URL")


class SourceError(ValueError):
    pass


@dataclass
class Source:
    kind: str                 # geo, sra, cellxgene, biostudies, figshare, zenodo, hca or url
    key: str                  # canonical identity, e.g. geo:GSE123, zenodo:456, url:https://...
    display: str              # what is stored and shown: GSE123, zenodo:456, the URL
    accession: str | None = None
    url: str | None = None
    raw: bool = False         # the raw reads of a GEO series / HCA project rather than its processed files
    related: list[str] = field(default_factory=list)  # keys of collections that contain this one

    @property
    def folder(self) -> str:
        """Folder name under downloads/<kind>/."""
        base = self.accession or util.slug(self.url.rsplit("/", 1)[-1]) or "download"
        return f"{base}-raw" if self.raw else base

    def mention_words(self) -> list[str]:
        """Words that, in a dataset's study, name or path, suggest it is this data."""
        a = self.accession
        if not a:
            return []
        if self.kind in ("zenodo", "figshare"):
            return [f"{self.kind}{sep}{a}" for sep in ("", ".", "-", "_", ":")]
        return [a]


def _make(kind: str, acc: str, raw: bool = False) -> Source:
    if raw and kind not in RAW_KINDS + ("sra",):
        raise SourceError(f"--raw only applies to GEO and HCA (and SRA is always raw reads), not {kind}")
    raw = raw and kind in RAW_KINDS
    display = acc if kind in ("geo", "sra", "biostudies") else f"{kind}:{acc}"
    if raw:
        display += ":raw"
    return Source(kind, f"{kind}:{acc}" + (":raw" if raw else ""), display, accession=acc, raw=raw)


_PATTERNS = (
    ("geo", re.compile(r"^((?:GSE|GSM)\d+)$", re.I), str.upper),
    ("sra", re.compile(r"^((?:[SED]R[PXRS]|PRJ(?:NA|EB|DB))\d+)$", re.I), str.upper),
    ("biostudies", re.compile(r"^(E-[A-Z]{4}-\d+|S-[A-Z]+\d+)$", re.I), str.upper),
    ("zenodo", re.compile(r"^(?:zenodo[:.]?|(?:https?://(?:dx\.)?doi\.org/)?10\.5281/zenodo\.)(\d+)$", re.I), str),
    ("figshare", re.compile(r"^(?:figshare[:.]?|(?:https?://(?:dx\.)?doi\.org/)?10\.6084/m9\.figshare\.)(\d+)(?:\.v\d+)?$",
                            re.I), str),
    ("cellxgene", re.compile(rf"^cellxgene:({UUID})$", re.I), str.lower),
    ("hca", re.compile(rf"^hca:({UUID})$", re.I), str.lower),
)


def normalize_url(url: str) -> str:
    u = urllib.parse.urlsplit(url.strip())
    scheme, host = u.scheme.lower(), (u.hostname or "").lower()
    if host in ("ftp.ncbi.nlm.nih.gov", "ftp.sra.ebi.ac.uk", "ftp.ebi.ac.uk") and scheme == "ftp":
        scheme = "https"  # the same files are served over https
    port = f":{u.port}" if u.port and (scheme, u.port) not in (("http", 80), ("https", 443), ("ftp", 21)) else ""
    path = u.path if u.path in ("", "/") else u.path.rstrip("/")
    return urllib.parse.urlunsplit((scheme, host + port, path, u.query, ""))


def _from_url(url: str, raw: bool) -> Source:
    u = urllib.parse.urlsplit(url)
    host, path = u.hostname or "", u.path
    q = urllib.parse.parse_qs(u.query)
    segs = [s for s in path.split("/") if s]

    def file_in(kind, acc):
        src = Source("url", f"url:{url}", url, url=url)
        src.related.append(f"{kind}:{acc}")
        return src

    if host.endswith("ncbi.nlm.nih.gov"):
        if path.endswith("acc.cgi") and q.get("acc") and re.match(r"^(GSE|GSM)\d+$", q["acc"][0], re.I):
            return _make("geo", q["acc"][0].upper(), raw)
        m = re.match(r"^/geo/(?:series|samples)/[^/]+/((?:GSE|GSM)\d+)(/.*)?$", path, re.I)
        if m:
            rest = (m.group(2) or "").strip("/")
            return _make("geo", m.group(1).upper(), raw) if rest in ("", "suppl") else file_in("geo", m.group(1).upper())
        if len(segs) >= 2 and segs[0] in ("bioproject", "sra") and re.match(r"^(PRJ|[SED]R)", segs[1], re.I):
            return _make("sra", segs[1].upper())
        if segs[:1] == ["sra"] and q.get("term"):
            return _make("sra", q["term"][0].upper())
    if host == "www.ebi.ac.uk":
        if "ena" in segs and segs[-1:] and re.match(r"^(PRJ|[SED]R[PXRS])\w+$", segs[-1], re.I):
            return _make("sra", segs[-1].upper())
        if ("biostudies" in segs or "arrayexpress" in segs) and segs and re.match(r"^(E-[A-Z]{4}-\d+|S-[A-Z]+\d+)$",
                                                                                     segs[-1], re.I):
            return _make("biostudies", segs[-1].upper())
    if host in ("zenodo.org", "www.zenodo.org"):
        m = re.match(r"^/(?:api/)?records?/(\d+)(/.*)?$", path)
        if m:
            return _make("zenodo", m.group(1)) if not (m.group(2) or "").strip("/") else file_in("zenodo", m.group(1))
    if host.endswith("figshare.com") and "articles" in segs:
        nums = [s for s in segs[segs.index("articles") + 1:] if s.isdigit()]
        if nums:
            return _make("figshare", nums[0])
    if host == "cellxgene.cziscience.com" and len(segs) >= 2 and segs[0] == "collections":
        return _make("cellxgene", segs[1].lower())
    if host.endswith("humancellatlas.org") and "projects" in segs:
        i = segs.index("projects")
        if i + 1 < len(segs) and re.match(UUID, segs[i + 1], re.I):
            return _make("hca", segs[i + 1].lower(), raw)
    if host in ("doi.org", "dx.doi.org"):
        return parse(path.lstrip("/"), raw)
    if raw:
        raise SourceError("--raw only applies to GEO and HCA sources")
    return Source("url", f"url:{url}", url, url=url)


def parse(text: str, raw: bool = False) -> Source:
    s = text.strip()
    if s.lower().endswith(":raw"):
        s, raw = s[:-4], True
    for kind, pat, norm in _PATTERNS:
        m = pat.match(s)
        if m:
            return _make(kind, norm(m.group(1)), raw)
    if not re.match(r"^(https?|ftp)://", s, re.I):
        raise SourceError(f"'{text}' is not {KINDS_HELP}")
    return _from_url(normalize_url(s), raw)


# --- remote side ----------------------------------------------------------------------

@dataclass
class RemoteFile:
    name: str
    url: str
    size: int | None = None
    md5: str | None = None
    sha256: str | None = None
    resolve: object = None   # callable giving a fresh download URL (signed URLs), if needed


@dataclass
class Remote:
    title: str | None = None
    species: str | None = None
    organ: str | None = None
    platform: str | None = None
    modalities: set[str] = field(default_factory=set)
    study: str | None = None
    files: list[RemoteFile] = field(default_factory=list)
    parents: list[str] = field(default_factory=list)   # keys of collections containing (or equal to) this
    note: str | None = None


def request(url: str, method: str = "GET", headers: dict | None = None):
    req = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT, **(headers or {})})
    return urllib.request.urlopen(req, timeout=TIMEOUT)


def get_text(url: str) -> str:
    try:
        with request(url) as r:
            return r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        raise SourceError(f"{url} answered {e.code} {e.reason}")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise SourceError(f"cannot reach {url}: {getattr(e, 'reason', e)}")


def get_json(url: str):
    text = get_text(url)
    try:
        return json.loads(text) if text.strip() else None
    except ValueError:
        raise SourceError(f"{url} returned something that isn't JSON")


def head(url: str) -> tuple[int | None, str | None]:
    """Size and file name a URL would download, where the server says."""
    try:
        with request(url, "HEAD") as r:
            size = r.headers.get("Content-Length")
            cd = r.headers.get("Content-Disposition") or ""
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        return None, None
    m = re.search(r"filename\*?=(?:UTF-8'')?\"?([^\";]+)\"?", cd)
    return (int(size) if size and size.isdigit() else None), (urllib.parse.unquote(m.group(1)) if m else None)


def _one(values) -> str | None:
    s = {v for v in values if v}
    return s.pop() if len(s) == 1 else None


def _organ(labels) -> str | None:
    return _one(vocab.infer_organ(vocab.tokens(t)) for t in labels)


_STRATEGY = {"rna-seq": "rna", "scrna-seq": "rna", "snrna-seq": "rna", "atac-seq": "atac", "bisulfite-seq": "methylation",
             "wgs": "dna", "wxs": "dna", "targeted-capture": "dna", "expression profiling": "rna",
             "methylation profiling": "methylation"}


def _strategy_modalities(values) -> set[str]:
    return {m for v in values for k, m in _STRATEGY.items() if v.lower().startswith(k)}


# --- GEO -------------------------------------------------------------------------------------

def _geo_fields(acc: str) -> dict[str, list[str]]:
    text = get_text(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief")
    fields: dict[str, list[str]] = {}
    for line in text.splitlines():
        m = re.match(r"^!(\w+) = (.*)$", line)
        if m:
            fields.setdefault(m.group(1), []).append(m.group(2).strip())
    return fields


def resolve_geo(src: Source) -> Remote:
    acc = src.accession
    fields = _geo_fields(acc)
    pre = "Series" if acc.startswith("GSE") else "Sample"
    if f"{pre}_title" not in fields:
        raise SourceError(f"GEO has no public record {acc} (it may be private or mistyped)")
    r = Remote(title=fields[f"{pre}_title"][0], study=acc)
    if pre == "Series":
        taxa = {vocab.TAXIDS.get(t) for t in fields.get("Series_sample_taxid", [])}
        r.modalities = _strategy_modalities(fields.get("Series_type", []))
    else:
        taxa = {vocab.species(o) for o in fields.get("Sample_organism_ch1", []) + fields.get("Sample_organism_ch2", [])}
        r.modalities = _strategy_modalities(fields.get("Sample_library_strategy", []))
        series = [s.upper() for s in fields.get("Sample_series_id", [])]
        r.parents = [f"geo:{s}" + (":raw" if src.raw else "") for s in series]
        if len(series) == 1:
            r.study = series[0]
    r.species = _one(taxa)
    if not src.raw:
        for url in fields.get(f"{pre}_supplementary_file", []):
            if url.upper() != "NONE":
                url = normalize_url(url)
                r.files.append(RemoteFile(url.rsplit("/", 1)[1], url))
        for f in r.files:
            f.size, _ = head(f.url)
        if not r.files:
            r.note = (f"{acc} has no supplementary (processed) files on GEO. "
                      f"For its raw reads, run `stew download {acc} --raw`.")
        return r
    # Raw reads: GEO points to SRA (or a BioProject), whose runs ENA mirrors.
    sra = []
    for rel in fields.get(f"{pre}_relation", []):
        m = re.search(r"\b((?:SR[PX]|PRJNA)\d+)\b", rel)
        if m and (rel.startswith("SRA") or rel.startswith("BioProject")):
            sra.append(m.group(1))
    if not sra:
        raise SourceError(f"GEO lists no SRA runs for {acc}, so there are no raw reads to fetch")
    acc_sra = next((a for a in sra if a.startswith(("SRP", "SRX"))), sra[0])
    reads = resolve_sra(_make("sra", acc_sra))
    r.files, r.note = reads.files, reads.note
    r.modalities |= reads.modalities
    r.species = r.species or reads.species
    r.parents += [f"sra:{a}" for a in sra]
    return r


# --- SRA via ENA -------------------------------------------------------------------------------

ENA_FIELDS = ("run_accession,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5,"
              "scientific_name,library_strategy,study_title,study_accession")


def resolve_sra(src: Source) -> Remote:
    acc = src.accession
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?" +
           urllib.parse.urlencode({"accession": acc, "result": "read_run", "fields": ENA_FIELDS,
                                   "format": "json", "limit": 0}))
    try:
        rows = get_json(url) or []
    except SourceError:
        rows = []
    if not rows:
        return resolve_sra_ncbi(acc)
    r = Remote(title=rows[0].get("study_title"), study=rows[0].get("study_accession") or acc)
    r.species = _one(vocab.species(x["scientific_name"]) for x in rows if x.get("scientific_name"))
    r.modalities = _strategy_modalities(x.get("library_strategy", "") for x in rows)
    for x in rows:
        kind = "fastq" if x.get("fastq_ftp") else "submitted"
        urls = [u for u in (x.get(f"{kind}_ftp") or "").split(";") if u]
        sizes = (x.get(f"{kind}_bytes") or "").split(";")
        md5s = (x.get(f"{kind}_md5") or "").split(";")
        for i, u in enumerate(urls):
            size = sizes[i] if i < len(sizes) else ""
            r.files.append(RemoteFile(u.rsplit("/", 1)[1], "https://" + u if "://" not in u else u,
                                      int(size) if size.isdigit() else None,
                                      md5s[i] if i < len(md5s) and md5s[i] else None))
    if not r.files:
        r.note = f"ENA lists runs for {acc} but no downloadable files."
    return r


# NCBI's run list, and the full .sra files in the SRA's AWS open-data bucket (plain
# HTTPS with resume, no SRA toolkit needed to download).
NCBI_RUNINFO = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={}"
SRA_ODP = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{0}/{0}"


def resolve_sra_ncbi(acc: str) -> Remote:
    """Runs that ENA doesn't mirror, as .sra files from NCBI's open-data copy."""
    import csv
    import io
    text = get_text(NCBI_RUNINFO.format(urllib.parse.quote(acc)))
    rows = [r for r in csv.DictReader(io.StringIO(text)) if r.get("Run")]
    if not rows:
        raise SourceError(f"neither ENA nor NCBI SRA lists public runs for {acc} (not released, or no reads "
                          "were deposited, e.g. for imaging-based spatial data)")
    r = Remote(title=f"SRA {rows[0].get('SRAStudy') or acc}", study=rows[0].get("SRAStudy") or acc)
    r.species = _one(vocab.species(x["ScientificName"]) for x in rows if x.get("ScientificName"))
    r.modalities = _strategy_modalities(x.get("LibraryStrategy", "") for x in rows)
    for x in rows:
        url = SRA_ODP.format(x["Run"])
        size, _ = head(url)
        r.files.append(RemoteFile(f"{x['Run']}.sra", url, size))
    r.note = (f"ENA doesn't mirror {acc}, so these are NCBI .sra files (reads with qualities, compact). "
              "Convert them with `fasterq-dump --split-3 <file>.sra` from the SRA toolkit when you need FASTQ.")
    return r


# --- CELLxGENE ----------------------------------------------------------------------------------

def resolve_cellxgene(src: Source) -> Remote:
    d = get_json(f"https://api.cellxgene.cziscience.com/curation/v1/collections/{src.accession}")
    if not d or "datasets" not in d:
        raise SourceError(f"CELLxGENE has no public collection {src.accession}")
    from .contents import _platform_from_assay
    r = Remote(title=d.get("name"), study=d.get("doi") or f"cellxgene:{src.accession}")
    organisms, tissues, assays = set(), set(), set()
    for ds in d["datasets"]:
        organisms |= {o["label"] for o in ds.get("organism") or []}
        tissues |= {t["label"] for t in ds.get("tissue") or []}
        assays |= {a["label"] for a in ds.get("assay") or []}
        if ds.get("spatial"):
            r.modalities.add("spatial")
        for a in ds.get("assets") or []:
            if a.get("filetype") == "H5AD":
                name = f"{util.slug(ds.get('title') or 'dataset')[:60]}-{ds['dataset_id'][:8]}.h5ad"
                r.files.append(RemoteFile(name, a["url"], a.get("filesize")))
    r.modalities.add("rna")
    if any("visium" in a.lower() or "slide-seq" in a.lower() for a in assays):
        r.modalities.add("spatial")
    r.species = _one(vocab.species(o) for o in organisms)
    r.organ = _organ(tissues)
    r.platform = _platform_from_assay(sorted(assays))
    if not r.files:
        r.note = f"CELLxGENE collection {src.accession} has no downloadable H5AD files."
    return r


# --- ArrayExpress / BioStudies -----------------------------------------------------------------

def _attrs(node) -> dict[str, str]:
    return {a.get("name", ""): a.get("value", "") for a in node.get("attributes") or []}


def resolve_biostudies(src: Source) -> Remote:
    acc = src.accession
    info = get_json(f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}/info")
    study = get_json(f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}")
    if not info or not study:
        raise SourceError(f"BioStudies/ArrayExpress has no public study {acc}")
    base = info["httpLink"].rstrip("/") + "/Files/"
    r = Remote(title=_attrs(study).get("Title"), study=acc)
    section = study.get("section") or {}
    sattrs = _attrs(section)
    r.species = vocab.species(sattrs["Organism"]) if sattrs.get("Organism") else None
    r.modalities = _strategy_modalities([sattrs.get("Study type", "")])
    if "single" in sattrs.get("Study type", "").lower() and "rna" in sattrs.get("Study type", "").lower():
        r.modalities.add("rna")

    def walk(node):
        for f in node.get("files") or []:
            for x in (f if isinstance(f, list) else [f]):
                p = x.get("path") or x.get("filePath")
                if p and x.get("type") != "directory":
                    r.files.append(RemoteFile(p, base + urllib.parse.quote(p), x.get("size")))
        for s in node.get("subsections") or []:
            for x in (s if isinstance(s, list) else [s]):
                walk(x)

    walk(section)
    if not r.files:  # big studies keep file lists elsewhere; read the folder listing instead
        for name in re.findall(r'href="([^"?/][^"]*)"', get_text(base)):
            if not name.endswith("/"):
                url = base + name
                size, _ = head(url)
                r.files.append(RemoteFile(urllib.parse.unquote(name), url, size))
    if not r.files:
        r.note = f"{acc} has no downloadable files. Its raw reads may be at ENA (try its ERP/PRJEB accession)."
    return r


# --- figshare and Zenodo ------------------------------------------------------------------------

def resolve_figshare(src: Source) -> Remote:
    d = get_json(f"https://api.figshare.com/v2/articles/{src.accession}")
    if not d:
        raise SourceError(f"figshare has no public article {src.accession}")
    r = Remote(title=d.get("title"), study=d.get("doi") or f"figshare:{src.accession}")
    for f in d.get("files") or []:
        if not f.get("is_link_only"):
            r.files.append(RemoteFile(f["name"], f["download_url"], f.get("size"), f.get("computed_md5") or None))
    if not r.files:
        r.note = f"figshare article {src.accession} has no downloadable files."
    return r


def resolve_zenodo(src: Source) -> Remote:
    d = get_json(f"https://zenodo.org/api/records/{src.accession}")
    if not d:
        raise SourceError(f"Zenodo has no public record {src.accession}")
    meta = d.get("metadata", {})
    r = Remote(title=meta.get("title"), study=f"zenodo:{src.accession}")
    for f in d.get("files") or []:
        md5 = f.get("checksum") or ""
        r.files.append(RemoteFile(f["key"], f["links"]["self"], f.get("size"),
                                  md5[4:] if md5.startswith("md5:") else None))
    if not r.files:
        r.note = f"Zenodo record {src.accession} has no downloadable files (access: {meta.get('access_right')})."
    return r


# --- Human Cell Atlas (Azul) ---------------------------------------------------------------------

AZUL = "https://service.azul.data.humancellatlas.org"
RAW_FORMATS = ("fastq", "fastq.gz", "fq.gz", "bam", "cram")


def _azul_url(uuid: str, version: str, catalog: str):
    """HCA files are served from signed URLs that expire, so ask for one just
    before downloading. Azul may ask us to wait while it prepares the file."""
    def resolve():
        url = f"{AZUL}/fetch/repository/files/{uuid}?" + urllib.parse.urlencode({"catalog": catalog,
                                                                                  "version": version})
        for _ in range(30):
            d = get_json(url)
            if d.get("Status") == 302:
                return d["Location"]
            time.sleep(int(d.get("Retry-After", 5)))
            url = d.get("Location", url)
        raise SourceError(f"HCA did not prepare file {uuid} in time")
    return resolve


def resolve_hca(src: Source) -> Remote:
    catalog = (get_json(f"{AZUL}/index/catalogs") or {}).get("default_catalog", "dcp")
    p = get_json(f"{AZUL}/index/projects/{src.accession}?catalog={catalog}")
    if not p or not p.get("projects"):
        raise SourceError(f"HCA has no project {src.accession}")
    proj = p["projects"][0]
    r = Remote(title=proj.get("projectTitle"), study=f"hca:{src.accession}")
    r.species = _one(vocab.species(g) for d in p.get("donorOrganisms") or [] for g in d.get("genusSpecies") or [])
    r.organ = _organ(o for s in p.get("specimens") or [] for o in s.get("organ") or [] if o)
    r.modalities = {"rna"} if any("rna" in str(x).lower() for x in p.get("protocols") or []) else set()
    # Projects can hold tens of thousands of FASTQs; ask only for the file formats wanted.
    base = {"projectId": {"is": [src.accession]}}
    probe = get_json(f"{AZUL}/index/files?" + urllib.parse.urlencode(
        {"catalog": catalog, "filters": json.dumps(base), "size": 1}))
    formats = [t["term"] for t in (probe.get("termFacets", {}).get("fileFormat", {}).get("terms") or [])
               if t.get("term") and (t["term"].lower() in RAW_FORMATS) == src.raw]
    if not formats:
        what = "raw read" if src.raw else "processed (non-FASTQ)"
        raise SourceError(f"HCA project {src.accession} has no {what} files" + ("" if src.raw else "; try --raw"))
    url = f"{AZUL}/index/files?" + urllib.parse.urlencode(
        {"catalog": catalog, "filters": json.dumps({**base, "fileFormat": {"is": formats}}), "size": 1000})
    while url:
        page = get_json(url)
        for hit in page.get("hits") or []:
            f = hit["files"][0]
            is_raw = (f.get("format") or "").lower() in RAW_FORMATS
            if is_raw == src.raw:
                r.files.append(RemoteFile(f["name"], f"{AZUL}/repository/files/{f['uuid']}", f.get("size"),
                                          sha256=f.get("sha256"),
                                          resolve=_azul_url(f["uuid"], f["version"], catalog)))
        url = (page.get("pagination") or {}).get("next")
    if not r.files:
        what = "raw read" if src.raw else "processed (non-FASTQ)"
        r.note = f"HCA project {src.accession} has no {what} files" + ("" if src.raw else "; try --raw")
    return r


# --- plain URLs -------------------------------------------------------------------------------------

def resolve_url(src: Source) -> Remote:
    size, name = head(src.url)
    if not name:
        name = urllib.parse.unquote(urllib.parse.urlsplit(src.url).path.rsplit("/", 1)[-1]) or "download"
    return Remote(files=[RemoteFile(name, src.url, size)])


RESOLVERS = {"geo": resolve_geo, "sra": resolve_sra, "cellxgene": resolve_cellxgene, "biostudies": resolve_biostudies,
             "figshare": resolve_figshare, "zenodo": resolve_zenodo, "hca": resolve_hca, "url": resolve_url}


def resolve(src: Source) -> Remote:
    return RESOLVERS[src.kind](src)


def strip_html(s: str | None) -> str | None:
    return html.unescape(re.sub(r"<[^>]+>", " ", s)).strip() if s else s
