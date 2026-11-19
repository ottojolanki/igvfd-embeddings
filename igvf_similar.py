#!/usr/bin/env python3
"""
Similar-experiment recommender + semantic search prototype for the IGVF portal.

Pipeline
  1. LOAD     embedded MeasurementSet JSON (portal API, a JSON dump, or igvfd test inserts)
  2. DESCRIBE each object two ways:
       - text   : a clean, human-readable "document" (no IDs, dates, lab, file lists)
       - facets : sets of normalized values (assay, biosample, organism, genes, ...)
  3. EMBED    the text with an open-source sentence-embedding model
  4. SCORE    hybrid = w_emb * cosine(text embeddings) + w_facet * weighted facet overlap
  5. OUTPUT   top-k similar experiments per accession, each with "why" explanations

Usage
  # from the live portal (run on a machine that can reach api.data.igvf.org)
  python igvf_similar.py build --source portal --out similar.json
  # from a JSON dump (list of embedded objects, or a search response with "@graph")
  python igvf_similar.py build --source dump.json --out similar.json
  # offline smoke test on igvfd's test inserts
  python igvf_similar.py build --source inserts:/path/to/igvfd/src/igvfd/tests/data/inserts --out similar.json
  # free-text semantic search over the same index
  python igvf_similar.py search --index similar.npz "CRISPR screens in K562 measuring proliferation"

Backends
  --backend st     sentence-transformers (default model Qwen/Qwen3-Embedding-0.6B)
  --backend tfidf  scikit-learn TF-IDF; no model download, useful for testing the pipeline
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.request
from collections import defaultdict

import numpy as np

PORTAL = "https://api.data.igvf.org"

# Fields requested from the portal search API. Dotted paths select from embedded objects,
# so one request returns everything needed without pulling the (large) embedded files list.
PORTAL_FIELDS = [
    "@id", "accession", "status", "summary", "description", "file_set_type",
    "preferred_assay_titles", "assay_titles", "assay_slims",
    "assay_term.term_name", "assay_term.term_id",
    "samples.@type", "samples.summary", "samples.classifications",
    "samples.sample_terms.term_name", "samples.sample_terms.term_id",
    "samples.sample_terms.organ_slims", "samples.sample_terms.cell_slims",
    "samples.sample_terms.system_slims",
    "samples.targeted_sample_term.term_name",
    "samples.disease_terms.term_name",
    "samples.treatments.treatment_term_name",
    "samples.modifications.modality",
    "samples.construct_library_sets.file_set_type",
    "donors.taxa", "donors.sex",
    "targeted_genes.symbol",
    "crispr_screen_readout", "crispr_screen_biometric.term_name",
    "functional_assay_mechanisms.term_name",
    "enrichment_designs", "control_types", "sequencing_library_types",
    "related_measurement_sets.measurement_sets",
    "control_file_sets",
    "lab.title",
]

# Facet weights for the structured part of the score. Lab is deliberately absent:
# "similar experiment" should mean similar biology/assay, not "same submitter".
FACET_WEIGHTS = {
    "assay": 0.30,          # preferred assay title(s)
    "assay_term": 0.10,     # OBI assay term
    "assay_slim": 0.05,     # broad assay category (soft match)
    "biosample": 0.25,      # sample term names
    "biosample_slim": 0.05, # organ/cell/system slims, classifications (soft match)
    "organism": 0.05,
    "genes": 0.10,          # targeted genes
    "perturbation": 0.05,   # readouts, mechanisms, modifications, treatments, CLS types
    "disease": 0.05,
}


# ----------------------------------------------------------------------------- helpers

def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def name_of(x, *keys):
    """Return a human-readable name from an embedded object or a raw link string."""
    if isinstance(x, dict):
        for k in keys or ("term_name", "symbol", "title", "name", "accession"):
            if x.get(k):
                return str(x[k])
        return None
    if isinstance(x, str) and not x.startswith("/") and not re.fullmatch(r"[0-9a-f-]{36}", x):
        return x
    return None  # unresolved link: not informative


def names(xs, *keys):
    out = []
    for x in as_list(xs):
        n = name_of(x, *keys)
        if n and n not in out:
            out.append(n)
    return out


def norm(s):
    return re.sub(r"\s+", " ", str(s).strip().lower())


def sample_type(sample):
    t = sample.get("@type")
    t = t[0] if isinstance(t, list) and t else t
    return re.sub(r"(?<!^)(?=[A-Z])", " ", t).lower() if t else None


# ----------------------------------------------------------------------------- describe

def describe(obj):
    """Turn one embedded MeasurementSet into (text, facets)."""
    samples = [s for s in as_list(obj.get("samples")) if isinstance(s, dict)]
    donors = [d for d in as_list(obj.get("donors")) if isinstance(d, dict)]
    for s in samples:  # donors may only be embedded under samples
        donors += [d for d in as_list(s.get("donors")) if isinstance(d, dict)]

    assay_titles = names(obj.get("preferred_assay_titles")) or names(obj.get("assay_titles"))
    assay_term = names(obj.get("assay_term"), "term_name")
    assay_slims = names(obj.get("assay_slims")) or names(
        (obj.get("assay_term") or {}).get("assay_slims") if isinstance(obj.get("assay_term"), dict) else None)

    sample_terms, slims, classes, stypes, targeted, disease, treat, mods, cls = ([] for _ in range(9))
    for s in samples:
        for t in as_list(s.get("sample_terms")):
            n = name_of(t, "term_name")
            if n:
                sample_terms.append(n)
            if isinstance(t, dict):
                for k in ("organ_slims", "cell_slims", "system_slims"):
                    slims += as_list(t.get(k))
        classes += as_list(s.get("classifications"))
        if sample_type(s):
            stypes.append(sample_type(s))
        targeted += names(s.get("targeted_sample_term"), "term_name")
        disease += names(s.get("disease_terms"), "term_name")
        treat += names(s.get("treatments"), "treatment_term_name", "term_name")
        mods += names(s.get("modifications"), "modality", "summary")
        cls += [c.get("file_set_type") for c in as_list(s.get("construct_library_sets"))
                if isinstance(c, dict) and c.get("file_set_type")]

    taxa = sorted({d.get("taxa") for d in donors if d.get("taxa")})
    sexes = sorted({d.get("sex") for d in donors if d.get("sex")})
    genes = names(obj.get("targeted_genes"), "symbol")
    readout = names(obj.get("crispr_screen_readout")) + names(obj.get("crispr_screen_biometric"), "term_name")
    mech = names(obj.get("functional_assay_mechanisms"), "term_name")
    enrich = names(obj.get("enrichment_designs"))
    seqlib = names(obj.get("sequencing_library_types"))

    uniq = lambda xs: list(dict.fromkeys(x for x in xs if x))

    # ---- text document. Order: most discriminative first; models weight early tokens more.
    parts = []
    if obj.get("summary"):
        parts.append(obj["summary"].rstrip(".") + ".")
    if assay_titles:
        extra = [t for t in assay_term if norm(t) not in {norm(a) for a in assay_titles}]
        parts.append(f"Assay: {', '.join(assay_titles)}"
                     + (f" ({', '.join(extra)})" if extra else "") + ".")
    if assay_slims:
        parts.append(f"Assay category: {', '.join(uniq(assay_slims))}.")
    if sample_terms:
        kinds = [k for k in uniq(classes + stypes) if norm(k) not in {norm(t) for t in sample_terms}]
        parts.append(f"Biosample: {', '.join(uniq(sample_terms))}"
                     + (f" ({', '.join(kinds)})" if kinds else "") + ".")
    if slims:
        parts.append(f"Tissue/cell context: {', '.join(uniq(slims))}.")
    if targeted:
        parts.append(f"Differentiated/targeted to: {', '.join(uniq(targeted))}.")
    if taxa:
        parts.append(f"Organism: {', '.join(taxa)}" + (f"; donor sex: {', '.join(sexes)}" if sexes else "") + ".")
    if disease:
        parts.append(f"Disease: {', '.join(uniq(disease))}.")
    if genes:
        parts.append(f"Targeted genes: {', '.join(genes[:25])}.")
    if readout or mech:
        parts.append(f"Readout: {', '.join(uniq(readout + mech))}.")
    if mods or treat or cls:
        parts.append(f"Perturbation: {', '.join(uniq(mods + treat + cls))}.")
    if enrich or seqlib:
        parts.append(f"Library: {', '.join(uniq(enrich + seqlib))}.")
    if obj.get("description"):
        parts.append(obj["description"].strip())
    text = " ".join(parts)

    facets = {
        "assay": {norm(x) for x in assay_titles},
        "assay_term": {norm(x) for x in assay_term},
        "assay_slim": {norm(x) for x in assay_slims},
        "biosample": {norm(x) for x in sample_terms + targeted},
        "biosample_slim": {norm(x) for x in slims + classes},
        "organism": {norm(x) for x in taxa},
        "genes": {x.upper() for x in genes},
        "perturbation": {norm(x) for x in readout + mech + mods + treat + cls},
        "disease": {norm(x) for x in disease},
    }
    return text, facets


# ----------------------------------------------------------------------------- similarity

def facet_matrix(facets_list, weights=FACET_WEIGHTS):
    """Weighted Jaccard over facets, normalized by the weight of facets present in either object.
    Returns (score matrix, per-pair shared-values dict for explanations)."""
    n = len(facets_list)
    num = np.zeros((n, n))
    den = np.zeros((n, n))
    for f, w in weights.items():
        sets = [fs[f] for fs in facets_list]
        # inverted index -> intersection counts, O(sum of postings^2) instead of O(n^2)
        inv = defaultdict(list)
        for i, s in enumerate(sets):
            for v in s:
                inv[v].append(i)
        inter = np.zeros((n, n))
        for idx in inv.values():
            ix = np.array(idx)
            inter[np.ix_(ix, ix)] += 1
        size = np.array([len(s) for s in sets], dtype=float)
        union = size[:, None] + size[None, :] - inter
        present = union > 0
        num += np.where(present, w * inter / np.maximum(union, 1), 0)
        den += np.where(present, w, 0)
    return np.where(den > 0, num / np.maximum(den, 1e-9), 0)


def explain(fa, fb):
    labels = {"assay": "assay", "biosample": "biosample", "genes": "targeted genes",
              "organism": "organism", "perturbation": "perturbation/readout",
              "disease": "disease", "biosample_slim": "tissue/cell type"}
    out = []
    for f, label in labels.items():
        shared = fa[f] & fb[f]
        if shared:
            out.append(f"Same {label}: {', '.join(sorted(shared)[:3])}")
    return out


# ----------------------------------------------------------------------------- embedding

class Embedder:
    def __init__(self, backend="st", model="Qwen/Qwen3-Embedding-0.6B", device=None):
        self.backend, self.model_name = backend, model
        if backend == "st":
            from sentence_transformers import SentenceTransformer
            self.model = SentenceTransformer(model, device=device)
        elif backend == "tfidf":
            from sklearn.feature_extraction.text import TfidfVectorizer
            self.model = TfidfVectorizer(ngram_range=(1, 2), sublinear_tf=True, min_df=1)
        else:
            raise ValueError(backend)

    def _query_kwargs(self):
        # Instruction-tuned models want a query prompt; documents are encoded without one.
        m = self.model_name.lower()
        if "qwen3-embedding" in m:
            return {"prompt": "Instruct: Given a search query, retrieve matching genomics "
                              "experiments from the IGVF data portal\nQuery: "}
        if "e5" in m:
            return {"prompt": "query: "}
        if "bge" in m and "m3" not in m:
            return {"prompt": "Represent this sentence for searching relevant passages: "}
        if "nomic" in m:
            return {"prompt": "search_query: "}
        return {}

    def _doc_kwargs(self):
        m = self.model_name.lower()
        if "e5" in m:
            return {"prompt": "passage: "}
        if "nomic" in m:
            return {"prompt": "search_document: "}
        return {}

    def docs(self, texts):
        if self.backend == "tfidf":
            X = self.model.fit_transform(texts).toarray()
            return X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-9)
        return self.model.encode(texts, batch_size=32, normalize_embeddings=True,
                                 show_progress_bar=len(texts) > 200, **self._doc_kwargs())

    def query(self, text):
        if self.backend == "tfidf":
            raise RuntimeError("query search needs the fitted vectorizer; use backend st")
        return self.model.encode([text], normalize_embeddings=True, **self._query_kwargs())[0]


# ----------------------------------------------------------------------------- loading

def load_portal(item_type="MeasurementSet", status="released"):
    q = [("type", item_type), ("status", status), ("limit", "all"), ("format", "json")]
    q += [("field", f) for f in PORTAL_FIELDS]
    url = f"{PORTAL}/search/?" + "&".join(f"{k}={urllib.request.quote(v)}" for k, v in q)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)["@graph"]


def load_inserts(path):
    """Resolve igvfd test inserts into an approximation of the portal's embedded frame."""
    types = {
        "tissue": "Tissue", "primary_cell": "PrimaryCell", "in_vitro_system": "InVitroSystem",
        "whole_organism": "WholeOrganism", "multiplexed_sample": "MultiplexedSample",
        "technical_sample": "TechnicalSample",
    }
    index = {}
    for fn in os.listdir(path):
        if not fn.endswith(".json"):
            continue
        kind = fn[:-5]
        for o in json.load(open(os.path.join(path, fn))):
            if kind in types:
                o["@type"] = [types[kind]]
            coll = "/" + kind.replace("_", "-") + "s/"
            for k in [o.get("uuid"), o.get("accession"), o.get("name"), o.get("term_id"),
                      *(o.get("aliases") or [])]:
                if k:
                    index[k] = o
                    index[coll + k + "/"] = o
            if kind == "lab" and o.get("name"):
                index["/labs/" + o["name"] + "/"] = o

    def r(x, depth):
        if depth == 0:
            return x
        if isinstance(x, str) and x in index:
            return r(index[x], depth - 1)
        if isinstance(x, list):
            return [r(v, depth) for v in x]
        if isinstance(x, dict):
            return {k: (r(v, depth) if k in RESOLVE else v) for k, v in x.items()}
        return x

    RESOLVE = {"samples", "sample_terms", "donors", "targeted_genes", "assay_term", "lab",
               "targeted_sample_term", "disease_terms", "treatments", "modifications",
               "construct_library_sets", "crispr_screen_biometric", "functional_assay_mechanisms",
               "phenotypic_features"}
    objs = json.load(open(os.path.join(path, "measurement_set.json")))
    return [r(o, 4) for o in objs if o.get("status") in ("released", "in progress", "preview")]


def load(source):
    if source == "portal":
        return load_portal()
    if source.startswith("inserts:"):
        return load_inserts(source.split(":", 1)[1])
    data = json.load(open(source))
    return data["@graph"] if isinstance(data, dict) else data


# ----------------------------------------------------------------------------- commands

def build(args):
    objs = [o for o in load(args.source) if o.get("accession")]
    print(f"loaded {len(objs)} objects", file=sys.stderr)
    described = [describe(o) for o in objs]
    texts = [t for t, _ in described]
    facets = [f for _, f in described]

    emb = Embedder(args.backend, args.model).docs(texts)
    cos = emb @ emb.T
    fac = facet_matrix(facets)
    score = args.w_emb * cos + (1 - args.w_emb) * fac
    np.fill_diagonal(score, -np.inf)

    # Experiments already linked to each other (multiome partners, controls) are "related",
    # not "similar": the portal shows those already, so exclude them here.
    acc_ix = {o["accession"]: i for i, o in enumerate(objs)}
    if not args.keep_related:
        for i, o in enumerate(objs):
            linked = []
            for g in as_list(o.get("related_measurement_sets")):
                linked += as_list(g.get("measurement_sets") if isinstance(g, dict) else g)
            linked += as_list(o.get("control_file_sets"))
            for l in linked:
                a = l.get("accession") if isinstance(l, dict) else str(l).rstrip("/").split("/")[-1]
                if a in acc_ix:
                    score[i, acc_ix[a]] = -np.inf

    out = {}
    for i, o in enumerate(objs):
        top = np.argsort(-score[i])[: args.k]
        out[o["accession"]] = [
            {"accession": objs[j]["accession"], "score": round(float(score[i, j]), 4),
             "embedding_sim": round(float(cos[i, j]), 4), "facet_sim": round(float(fac[i, j]), 4),
             "why": explain(facets[i], facets[j])}
            for j in top if np.isfinite(score[i, j]) and score[i, j] >= args.min_score
        ]
    json.dump({"model": args.model if args.backend == "st" else "tfidf", "w_emb": args.w_emb,
               "similar": out}, open(args.out, "w"), indent=1)
    np.savez_compressed(os.path.splitext(args.out)[0] + ".npz", emb=emb,
                        accessions=np.array([o["accession"] for o in objs]),
                        texts=np.array(texts))
    if args.dump_texts:
        with open(args.dump_texts, "w") as f:
            for o, t in zip(objs, texts):
                f.write(json.dumps({"accession": o["accession"], "text": t}) + "\n")
    print(f"wrote {args.out}", file=sys.stderr)


def search(args):
    z = np.load(args.index)
    q = Embedder("st", args.model).query(args.query)
    sims = z["emb"] @ q
    for j in np.argsort(-sims)[: args.k]:
        print(f"{sims[j]:.3f}  {z['accessions'][j]}  {z['texts'][j][:140]}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--source", default="portal")
    b.add_argument("--out", default="similar.json")
    b.add_argument("--backend", choices=["st", "tfidf"], default="st")
    b.add_argument("--model", default="Qwen/Qwen3-Embedding-0.6B")
    b.add_argument("--w-emb", type=float, default=0.5, help="weight of embedding vs facet score")
    b.add_argument("--k", type=int, default=10)
    b.add_argument("--min-score", type=float, default=0.2)
    b.add_argument("--keep-related", action="store_true")
    b.add_argument("--dump-texts", help="write the generated text documents as JSONL")
    s = sub.add_parser("search")
    s.add_argument("query")
    s.add_argument("--index", default="similar.npz")
    s.add_argument("--model", default="Qwen/Qwen3-Embedding-0.6B")
    s.add_argument("--k", type=int, default=10)
    args = p.parse_args()
    {"build": build, "search": search}[args.cmd](args)


if __name__ == "__main__":
    main()
