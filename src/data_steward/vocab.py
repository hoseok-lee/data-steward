"""Controlled vocabulary: canonical names and aliases for modalities, organs,
species and platforms, and inference of those from path names and titles."""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Term:
    name: str
    about: str = ""
    aliases: tuple[str, ...] = ()
    # Modalities implied by a platform (empty for other terms).
    implies: tuple[str, ...] = field(default=())


def _t(name, about, *aliases):
    return Term(name, about, tuple(aliases))


def _p(name, aliases, implies):
    return Term(name, "", tuple(aliases), tuple(implies))


MODALITIES = (
    _t("rna", "gene expression (sc/snRNA-seq, spatial transcriptomics)",
       "gex", "expression", "transcriptome", "transcriptomics", "scrna", "snrna", "scrnaseq", "snrnaseq", "rnaseq"),
    _t("atac", "chromatin accessibility", "scatac", "snatac", "atacseq", "accessibility", "chromatin"),
    _t("protein", "protein abundance (CITE-seq ADT, CODEX, IMC, proteomics)",
       "adt", "cite", "citeseq", "proteomics", "codex", "imc"),
    _t("spatial", "spatially resolved measurements", "spatially-resolved"),
    _t("imaging", "fluorescence or other microscopy images", "image", "images", "microscopy", "if"),
    _t("histology", "H&E or other stained whole-slide images", "he", "h&e", "hne", "wsi", "slide"),
    _t("dna", "DNA sequencing (WGS, WES, CNV)", "wgs", "wes", "genome", "exome"),
    _t("methylation", "DNA methylation", "bisulfite", "5mc"),
    _t("vdj", "immune repertoire (TCR, BCR)", "tcr", "bcr", "repertoire"),
)

# Aliases too generic to infer a modality from a path name.
_NO_INFER = {"slide", "image", "images", "genome", "expression", "accessibility", "chromatin"}

ORGANS = (
    _t("lung", "", "lungs", "pulmonary", "luad", "lusc", "nsclc", "airway"),
    _t("liver", "", "hepatic", "hcc", "hepatocellular"),
    _t("breast", "", "mammary", "brca"),
    _t("pancreas", "", "pancreatic", "pdac", "paad"),
    _t("brain", "", "cerebral", "cortex", "cortical", "glioma", "gbm", "glioblastoma", "hippocampus", "cerebellum"),
    _t("kidney", "", "renal", "rcc", "nephron"),
    _t("colon", "", "colorectal", "crc", "coad", "rectum", "rectal"),
    _t("small-intestine", "", "intestine", "intestinal", "ileum", "jejunum", "duodenum", "gut"),
    _t("stomach", "", "gastric"),
    _t("esophagus", "", "oesophagus", "esophageal", "oesophageal"),
    _t("skin", "", "dermal", "epidermis", "melanoma"),
    _t("heart", "", "cardiac", "myocardium"),
    _t("prostate", "", "prad"),
    _t("ovary", "", "ovarian"),
    _t("uterus", "", "endometrium", "endometrial", "uterine"),
    _t("lymph-node", "", "lymph", "ln"),
    _t("spleen", "", "splenic"),
    _t("bone-marrow", "", "marrow", "bm"),
    _t("blood", "", "pbmc", "pbmcs", "peripheral"),
    _t("tonsil", "", "tonsils"),
    _t("thymus", ""),
    _t("bladder", "", "urothelial"),
    _t("thyroid", ""),
    _t("head-and-neck", "", "hnscc", "oral"),
    _t("placenta", "", "placental"),
    _t("embryo", "", "embryonic"),
    _t("eye", "", "retina", "retinal"),
    _t("muscle", "", "skeletal"),
    _t("adipose", "", "fat"),
    _t("testis", "", "testes"),
    _t("bone", ""),
)

SPECIES = (
    _t("human", "", "homo", "sapiens", "homo-sapiens", "hs", "hsa", "hg19", "hg38", "grch37", "grch38",
       "patient", "patients"),
    _t("mouse", "", "mus", "musculus", "mus-musculus", "mm", "mm10", "mm39", "grcm38", "grcm39", "murine", "mice"),
    _t("rat", "", "rattus", "rattus-norvegicus", "rn6", "rn7"),
    _t("zebrafish", "", "danio", "rerio", "danio-rerio", "danrer11"),
    _t("fly", "", "drosophila", "drosophila-melanogaster", "dm6"),
    _t("macaque", "", "macaca", "rhesus", "cynomolgus", "macaca-mulatta"),
    _t("pig", "", "sus", "porcine", "sus-scrofa"),
)

# NCBI taxonomy ids, as GEO reports them.
TAXIDS = {"9606": "human", "10090": "mouse", "10116": "rat", "7955": "zebrafish", "7227": "fly",
          "9544": "macaque", "9823": "pig"}

PLATFORMS = (
    _p("xenium", ["xetg"], ["rna", "spatial"]),  # XETG: Xenium instrument serials in file names
    _p("visium", ["spaceranger"], ["rna", "spatial"]),
    _p("visium-hd", ["visiumhd"], ["rna", "spatial"]),
    _p("cosmx", ["smi"], ["rna", "spatial"]),
    _p("merscope", ["merfish", "vizgen"], ["rna", "spatial"]),
    _p("slide-seq", ["slideseq", "slideseqv2", "curio"], ["rna", "spatial"]),
    _p("stereo-seq", ["stereoseq", "stomics"], ["rna", "spatial"]),
    _p("geomx", ["dsp"], ["rna", "spatial"]),
    _p("codex", ["phenocycler"], ["protein", "spatial", "imaging"]),
    _p("imc", ["hyperion"], ["protein", "spatial", "imaging"]),
    _p("10x-chromium", ["chromium", "cellranger", "10x"], ["rna"]),
    _p("10x-atac", ["cellranger-atac"], ["atac"]),
    _p("10x-multiome", ["multiome"], ["rna", "atac"]),
    _p("smart-seq", ["smartseq", "smartseq2", "smartseq3"], ["rna"]),
)


def _key(s: str) -> str:
    return re.sub(r"[ _]", "-", s.strip().lower())


def _lookup(terms, s: str) -> Term | None:
    k = _key(s)
    for t in terms:
        if t.name == k or k in t.aliases:
            return t
    return None


class VocabError(ValueError):
    pass


def modality(s: str) -> str:
    t = _lookup(MODALITIES, s)
    if t is None:
        known = ", ".join(t.name for t in MODALITIES)
        raise VocabError(f"unknown modality '{s}'. Known: {known} (see `stew vocab`)")
    return t.name


def organ(s: str) -> str:
    t = _lookup(ORGANS, s)
    return t.name if t else _key(s)


def species(s: str) -> str:
    t = _lookup(SPECIES, s)
    return t.name if t else _key(s)


def platform(s: str) -> str:
    t = _lookup(PLATFORMS, s)
    return t.name if t else _key(s)


def platform_modalities(name: str) -> tuple[str, ...]:
    t = _lookup(PLATFORMS, name)
    return t.implies if t else ()


def tokens(s: str) -> list[str]:
    """Lowercase alphanumeric tokens; "Lung5" also yields "lung", and joined
    words are split too ("humanLung" gives "human" and "lung")."""
    out = []
    split = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", s)
    words = re.split(r"[^0-9a-z]+", split.lower())
    if split != s:
        words += [w for w in re.split(r"[^0-9a-z]+", s.lower()) if w not in words]
    for tok in words:
        if not tok:
            continue
        out.append(tok)
        stripped = tok.rstrip("0123456789")
        if stripped and stripped != tok:
            out.append(stripped)
    return out


def _infer(terms, toks) -> str | None:
    """The last token naming a term (the deepest folder, for paths). Tokens
    under three characters are too ambiguous to use."""
    for tok in reversed(toks):
        if len(tok) >= 3:
            t = _lookup(terms, tok)
            if t:
                return t.name
    return None


def infer_organ(toks) -> str | None:
    return _infer(ORGANS, toks)


def infer_species(toks) -> str | None:
    return _infer(SPECIES, toks)


def all_species(toks) -> set[str]:
    return {t.name for tok in toks if len(tok) >= 3 for t in [_lookup(SPECIES, tok)] if t}


def infer_platform(toks) -> str | None:
    if any(a == "visium" and b == "hd" for a, b in zip(toks, toks[1:])):
        return "visium-hd"
    return _infer(PLATFORMS, toks)


def infer_modalities(toks) -> set[str]:
    out = set()
    for tok in toks:
        if len(tok) >= 3 and tok not in _NO_INFER:
            t = _lookup(MODALITIES, tok)
            if t:
                out.add(t.name)
    return out
