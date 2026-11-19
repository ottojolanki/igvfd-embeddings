"""Turn each Set into (a) a readable text document and (b) structured facets.

    python -m igvfsim.describe --variant full            # -> data/descriptions/full.jsonl
    python -m igvfsim.describe --all-variants

Text deliberately leaves out IDs, dates, files, lab, award, collections and publications:
"similar" should mean similar biology/assay, and the link fields are the evaluation proxies.
Sample summaries are kept but scrubbed of batch identifiers (cellular sub-pool codes) unless a
variant asks for them (to measure how much batch text leaks into similarity).
"""
from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass, field

from . import ontology
from .coverage import get_path
from .load import DATA, read

SET_TYPES = ["MeasurementSet", "AnalysisSet", "PredictionSet"]

MAX_LIST = 12  # cap long lists (multiplexed samples, assessed genes) in the text


def uniq(xs):
    return list(dict.fromkeys(x for x in xs if x not in (None, "")))


def norm(s) -> str:
    return re.sub(r"\s+", " ", str(s).strip().lower())


def capped(xs, n=MAX_LIST) -> str:
    xs = uniq(xs)
    return ", ".join(xs[:n]) + (f" and {len(xs) - n} more" if len(xs) > n else "")


_SUBPOOL = re.compile(r"[,;]?\s*\(?cellular sub[- ]pool:?\s*[^,;)]*\)?", re.I)


def scrub_sample_summary(s: str) -> str:
    s = _SUBPOOL.sub("", s)
    return re.sub(r"\s+", " ", s).strip(" ,;")


# ID-like tokens in free text: accessions, PMIDs, and lab batch codes such as
# "donor_SV20_day45_BR_batch2", "pool6:gex_ox_il_2", "igvf_008b".
_IDLIKE = re.compile(
    r"\b(?:ENC[A-Z]{2}\d{3}[A-Z]{3}|IGVF[A-Z]{2}\d{4}[A-Z]{4}|PMID:?\s*\d+|PMCID:?\s*PMC\d+"
    r"|(?=[\w:.-]*\d)(?=[\w:.-]*[_:])[A-Za-z0-9]+(?:[_:.-][A-Za-z0-9]+)+)\b")


def scrub_free_text(s: str) -> str:
    s = _IDLIKE.sub("", s)
    s = re.sub(r"\s+([.;,])", r"\1", s)
    s = re.sub(r"([;,])(\s*[;,.])+", r"\2", s)
    return re.sub(r"\s+", " ", s).strip(" ,;")


@dataclass
class Parsed:
    """Everything describe/facets need, pulled once from the raw object."""
    type: str
    accession: str
    summary: str = ""
    description: str = ""
    assay_titles: list = field(default_factory=list)   # preferred titles (MS/AS); file_set_type (PS)
    assay_terms: list = field(default_factory=list)    # OBI term names (MS/AS); software (PS)
    assay_slims: list = field(default_factory=list)
    assay_term_ids: list = field(default_factory=list)
    sample_terms: list = field(default_factory=list)
    sample_term_ids: list = field(default_factory=list)
    targeted_terms: list = field(default_factory=list)
    classifications: list = field(default_factory=list)
    sample_summaries: list = field(default_factory=list)
    sample_summaries_raw: list = field(default_factory=list)
    slims: list = field(default_factory=list)          # context slims, per-sample logic (see parse)
    description_raw: str = ""
    term_notes: list = field(default_factory=list)     # (name, synonyms, first sentence of definition)
    taxa: list = field(default_factory=list)
    strains: list = field(default_factory=list)
    genes: list = field(default_factory=list)
    readout: list = field(default_factory=list)
    library: list = field(default_factory=list)
    modifications: list = field(default_factory=list)
    treatments: list = field(default_factory=list)
    phenotypes: list = field(default_factory=list)
    workflows: list = field(default_factory=list)
    content_types: list = field(default_factory=list)
    scope: str = ""
    file_set_type: str = ""


def parse(o: dict) -> Parsed:
    t = o["@type"][0]
    p = Parsed(type=t, accession=o["accession"])
    g = lambda path: uniq(get_path(o, path))
    p.summary = (o.get("summary") or "").strip()
    p.description_raw = (o.get("description") or "").strip()
    p.description = scrub_free_text(p.description_raw)
    p.file_set_type = o.get("file_set_type") or ""
    if t == "PredictionSet":
        p.assay_titles = [p.file_set_type] if p.file_set_type else []
        p.assay_terms = g("software.title")
        p.genes = g("assessed_genes.symbol") + g("small_scale_gene_list.symbol")
        p.phenotypes = g("associated_phenotypes.term_name")
        p.scope = o.get("scope") or ""
    else:
        p.assay_titles = g("preferred_assay_titles") or g("assay_titles")
        p.assay_terms = g("assay_term.term_name") or g("assay_titles")
        p.assay_slims = g("assay_term.assay_slims")
        p.assay_term_ids = g("assay_term.@id")
        p.genes = g("targeted_genes.symbol")
    if t == "MeasurementSet":
        p.readout = g("crispr_screen_readout") + g("crispr_screen_biometric.term_name") + \
            g("functional_assay_mechanisms.term_name")
        p.library = g("sequencing_library_types")
    if t == "AnalysisSet":
        p.workflows = g("workflows.name")
        p.readout = g("functional_assay_mechanisms.term_name")
    if t in ("AnalysisSet", "PredictionSet"):
        p.content_types = g("files.content_type")
    p.sample_terms = g("samples.sample_terms.term_name")
    p.sample_term_ids = g("samples.sample_terms.@id")
    p.targeted_terms = g("samples.targeted_sample_term.term_name")
    p.sample_term_ids += g("samples.targeted_sample_term.@id")
    p.classifications = g("samples.classifications")
    p.sample_summaries_raw = g("samples.summary")
    p.sample_summaries = uniq(scrub_sample_summary(s) for s in p.sample_summaries_raw)
    p.taxa = g("donors.taxa") or g("samples.taxa")
    p.strains = g("donors.strain_background")
    p.modifications = g("samples.modifications.modality")
    p.slims = context_slims(o)
    ont = ontology.get()
    for tid in uniq(p.assay_term_ids + p.sample_term_ids)[:4]:
        t = ont.sample.get(tid) or ont.assay.get(tid) or {}
        d = re.split(r"(?<=[a-z0-9)])\.\s", (t.get("definition") or "").strip())[0].rstrip(".")[:220]
        syn = [x for x in (t.get("synonyms") or []) if norm(x) != norm(t.get("term_name", ""))][:4]
        if d or syn:
            p.term_notes.append((t.get("term_name", ""), syn, d))
    p.treatments = g("samples.treatments.treatment_term_name")
    return p


def context_slims(o: dict) -> list[str]:
    """Organ/cell slims describing what the sample *is*, decided per sample.

    Slims on a cell line's term describe its tissue of origin (WTC-11 -> "skin of body",
    HAP-1 -> "bone marrow"), which misleads. So: if the sample was differentiated/targeted,
    use the target's slims; for cell lines use cell slims only; otherwise organ + cell slims.
    System/developmental slims (e.g. "mesoderm") are too broad to help and are left out.
    """
    ont = ontology.get()
    out = []
    for smp in o.get("samples") or []:
        if not isinstance(smp, dict):
            continue
        classes = set(smp.get("classifications") or [])
        tgt = [t.get("@id") for t in _as_list(smp.get("targeted_sample_term")) if isinstance(t, dict)]
        terms = tgt or [t.get("@id") for t in _as_list(smp.get("sample_terms")) if isinstance(t, dict)]
        keys = ("cell_slims",) if ("cell line" in classes and not tgt) else ("organ_slims", "cell_slims")
        for tid in terms:
            t = ont.sample.get(tid) or {}
            for k in keys:
                out += t.get(k) or []
    return uniq(out)


def _as_list(x):
    return [] if x is None else x if isinstance(x, list) else [x]


# ----------------------------------------------------------------------------- text

def _sections(p: Parsed, keep_batch: bool) -> dict[str, str]:
    """Named text sections; variants choose which to include and in what order."""
    s = {}
    if p.summary:
        s["summary"] = p.summary.rstrip(".") + "."
    kind = {"MeasurementSet": "Assay", "AnalysisSet": "Analysis of assay",
            "PredictionSet": "Prediction type"}[p.type]
    if p.assay_titles:
        extra = [a for a in p.assay_terms if norm(a) not in {norm(x) for x in p.assay_titles}]
        label = "software" if p.type == "PredictionSet" else None
        tail = (f" ({label}: {', '.join(extra)})" if label else f" ({', '.join(extra)})") if extra else ""
        s["assay"] = f"{kind}: {', '.join(p.assay_titles)}{tail}."
    if p.assay_slims:
        s["assay_slim"] = f"Assay category: {', '.join(p.assay_slims)}."
    if p.scope:
        s["scope"] = f"Scope: {p.scope}."
    if p.sample_terms:
        kinds = [k for k in p.classifications if norm(k) not in {norm(t) for t in p.sample_terms}]
        s["biosample"] = f"Biosample: {capped(p.sample_terms)}" + (f" ({', '.join(kinds)})" if kinds else "") + "."
    if p.targeted_terms:
        s["targeted"] = f"Differentiated/targeted to: {capped(p.targeted_terms)}."
    if p.slims:
        s["slims"] = f"Tissue/cell context: {capped(p.slims, 10)}."
    if p.taxa:
        s["organism"] = f"Organism: {', '.join(p.taxa)}" + (
            f" (strains: {capped(p.strains, 6)})" if p.strains else "") + "."
    if p.genes:
        s["genes"] = f"{'Assessed' if p.type == 'PredictionSet' else 'Targeted'} genes: {capped(p.genes, 25)}."
    if p.phenotypes:
        s["phenotypes"] = f"Phenotypes: {capped(p.phenotypes)}."
    if p.readout:
        s["readout"] = f"Readout: {', '.join(p.readout)}."
    if p.modifications:
        s["perturbation"] = f"CRISPR modification: {', '.join(p.modifications)}."
    if p.treatments:
        s["treatment"] = f"Treatment: {capped(p.treatments)}."
    if p.library:
        s["library"] = f"Library: {', '.join(p.library)}."
    if p.workflows:
        s["workflow"] = f"Workflow: {', '.join(p.workflows)}."
    if p.content_types:
        s["outputs"] = f"Outputs: {capped(p.content_types, 8)}."
    ss = uniq(p.sample_summaries_raw) if keep_batch else p.sample_summaries
    if ss:
        s["sample_summary"] = "Sample: " + "; ".join(ss[:3]) + (f"; and {len(ss) - 3} more samples" if len(ss) > 3 else "") + "."
    syn = [f"{n} (also: {', '.join(sy)})" for n, sy, _ in p.term_notes if sy]
    if syn:
        s["synonyms"] = "Synonyms: " + "; ".join(syn) + "."
    defs = [f"{n}: {d}" for n, _, d in p.term_notes if d]
    if defs:
        s["definitions"] = "Definitions: " + "; ".join(defs) + "."
    desc = p.description_raw if keep_batch else p.description
    if desc:
        s["description"] = desc
    return s


STRUCTURED = ["assay", "assay_slim", "scope", "biosample", "targeted", "slims", "organism", "genes",
              "phenotypes", "readout", "perturbation", "treatment", "library", "workflow", "outputs"]

# variant name -> (section order, keep batch ids)
VARIANTS: dict[str, tuple[list[str], bool]] = {
    # prototype order: summary, structured fields, then free text
    "full": (["summary"] + STRUCTURED + ["sample_summary", "description"], False),
    "no_description": (["summary"] + STRUCTURED + ["sample_summary"], False),
    "no_summary": (STRUCTURED + ["sample_summary", "description"], False),
    "structured_only": (STRUCTURED, False),
    "structured_first": (STRUCTURED + ["summary", "sample_summary", "description"], False),
    "free_text_only": (["summary", "sample_summary", "description"], False),
    "full_with_batch_ids": (["summary"] + STRUCTURED + ["sample_summary", "description"], True),
    # ontology enrichment: term synonyms / first sentence of term definitions (SampleTerm, AssayTerm)
    "full_syn": (["summary"] + STRUCTURED + ["synonyms", "sample_summary", "description"], False),
    "full_defs": (["summary"] + STRUCTURED + ["synonyms", "definitions", "sample_summary", "description"], False),
}


def text(p: Parsed, variant: str = "full") -> str:
    order, keep_batch = VARIANTS[variant]
    s = _sections(p, keep_batch)
    return " ".join(s[k] for k in order if k in s)


# ----------------------------------------------------------------------------- facets

def facets(p: Parsed) -> dict[str, dict[str, str]]:
    """{facet: {normalized value: display value}}. No links, lab, collections, publications."""
    ont = ontology.get()
    d = lambda xs: {norm(x): str(x) for x in xs if x}
    onto = uniq(a for tid in p.sample_term_ids for a in ont.closure(tid))
    assay_onto = uniq(a for tid in p.assay_term_ids for a in ont.closure(tid))
    return {
        "assay": d(p.assay_titles),
        "assay_term": d(p.assay_terms),
        "assay_slim": d(p.assay_slims),
        "assay_onto": d(assay_onto),
        "biosample": d(p.sample_terms + p.targeted_terms),
        "biosample_slim": d(p.slims + p.classifications),
        "biosample_onto": d(onto),
        "organism": d(p.taxa),
        "genes": {g.upper(): g for g in p.genes},
        "perturbation": d(p.readout + p.modifications + p.treatments + p.library),
        "phenotype": d(p.phenotypes),
        "workflow": d(p.workflows),
    }


# ----------------------------------------------------------------------------- corpus

@dataclass
class Corpus:
    objs: list[dict]
    parsed: list[Parsed]

    @property
    def accessions(self):
        return [p.accession for p in self.parsed]

    @property
    def types(self):
        return [p.type for p in self.parsed]

    def texts(self, variant="full"):
        return [text(p, variant) for p in self.parsed]

    def facets(self):
        return [facets(p) for p in self.parsed]


def load_corpus(types=SET_TYPES) -> Corpus:
    objs = [o for t in types for o in read(t)]
    return Corpus(objs, [parse(o) for o in objs])


def write_jsonl(corpus: Corpus, variant: str, out: str):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        for p, t in zip(corpus.parsed, corpus.texts(variant)):
            f.write(json.dumps({"accession": p.accession, "type": p.type, "text": t}) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", default="full", choices=list(VARIANTS))
    ap.add_argument("--all-variants", action="store_true")
    a = ap.parse_args()
    c = load_corpus()
    for v in (VARIANTS if a.all_variants else [a.variant]):
        out = os.path.join(DATA, "descriptions", f"{v}.jsonl")
        write_jsonl(c, v, out)
        texts = c.texts(v)
        print(f"{v}: {len(texts)} docs, {len(set(texts))} distinct, "
              f"mean {sum(map(len, texts)) / len(texts):.0f} chars -> {out}")


if __name__ == "__main__":
    main()
