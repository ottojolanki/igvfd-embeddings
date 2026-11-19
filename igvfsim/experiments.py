"""Run an experiment sweep from a YAML config; results/<name>/*.json, then EXPERIMENTS.md is rebuilt.

    python -m igvfsim.experiments configs/baseline.yaml
    python -m igvfsim.experiments configs/baseline.yaml --only bge-small   # subset of runs

Config:
    name: baseline
    max_group: 200          # drop proxy groups larger than this
    pool: within            # within | all  (cross-type candidates)
    eyeball: 15             # number of stratified queries to dump rankings for (runs with eyeball: true)
    runs:
      - name: bge-small
        embedder: st        # st | tfidf | bm25 | none | random | oracle | text_ceiling
        model: BAAI/bge-small-en-v1.5
        variant: full       # describe.VARIANTS
        prompt: sym_instruct  # optional doc prompt override
        facets: exact       # exact | onto | onto_only | none
        w_emb: [0.0, 0.5, 1.0]
        eyeball: true
        note: free text
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time

import numpy as np
import yaml

from . import ontology
from .describe import load_corpus
from .embed import SYM_INSTRUCT, embed_texts, tfidf_matrix
from .evaluate import build_positives, evaluate, write_rankings_md
from .sample import stratified_indices
from .score import FacetIndex, Scorer, bm25_sim
from . import report

def oracle_sim(pos, n):
    def f(q, c):
        out = np.zeros((len(q), len(c)), np.float32)
        cix = {int(j): r for r, j in enumerate(c)}
        for r, i in enumerate(q):
            for j in pos[i]:
                if j in cix:
                    out[r, cix[j]] = 1.0
        return out
    return f


def text_ceiling_sim(pos, texts):
    """Best any text-only scorer can do: candidates whose description is identical to some
    positive's description score 1 (they are indistinguishable from it), others 0."""
    tid = {}
    cls = np.array([tid.setdefault(t, len(tid)) for t in texts])

    def f(q, c):
        out = np.zeros((len(q), len(c)), np.float32)
        cc = cls[c]
        for r, i in enumerate(q):
            pc = np.unique(cls[list(pos[i])]) if pos[i] else []
            out[r] = np.isin(cc, pc)
        return out
    return f


def random_sim(seed=0):
    rng = np.random.default_rng(seed)
    return lambda q, c: rng.random((len(q), len(c)), dtype=np.float32)


def run_config(path: str, only: list[str] | None = None, device: str | None = None):
    cfg = yaml.safe_load(open(path))
    name = cfg.get("name", os.path.splitext(os.path.basename(path))[0])
    outdir = os.path.join("results", name)
    os.makedirs(outdir, exist_ok=True)
    corpus = load_corpus()
    facets = corpus.facets()
    pos, pstats = build_positives(corpus, cfg.get("max_group", 200))
    pool = cfg.get("pool", "within")
    onto = ontology.get()
    types = corpus.types
    has_pos = [i for i in range(len(types)) if any(types[j] == types[i] for j in pos[i])]
    sub = [corpus.parsed[i] for i in has_pos]
    eyeball_q = [has_pos[i] for i in stratified_indices(sub, cfg.get("eyeball", 15), seed=1,
                                                         per_type={"MeasurementSet": 8, "AnalysisSet": 4, "PredictionSet": 3})]
    fidx_cache = {}
    for run in cfg["runs"]:
        if only and run["name"] not in only:
            continue
        emb_kind = run.get("embedder", "st")
        variant = run.get("variant", "full")
        fscheme = run.get("facets", "none" if emb_kind in ("random", "oracle", "text_ceiling") else "exact")
        ws = run.get("w_emb", [1.0] if fscheme == "none" else [0.0] if emb_kind == "none" else [0.5])
        ws = ws if isinstance(ws, list) else [ws]
        model, dev, embed_s, text_src = run.get("model"), None, "—", None
        if emb_kind == "st":
            prompt = run.get("prompt")
            prompt = SYM_INSTRUCT if prompt == "sym_instruct" else prompt
            emb, st = embed_texts(corpus.texts(variant), model, prompt=prompt, device=device,
                                  max_len=run.get("max_len", 512))
            text_src, dev = emb, st["device"] or st["original"]["device"]
            embed_s = st["seconds"] if st["encoded"] else st["original"]["seconds"]
            embed_s = f"{embed_s} ({st['unique']} uniq)"
        elif emb_kind == "tfidf":
            t0 = time.time()
            text_src, _ = tfidf_matrix(corpus.texts(variant))
            embed_s, dev, model = round(time.time() - t0, 1), "cpu", "tfidf(1-2gram)"
        elif emb_kind == "bm25":
            t0 = time.time()
            text_src = bm25_sim(corpus.texts(variant))
            embed_s, dev, model = round(time.time() - t0, 1), "cpu", "bm25"
        elif emb_kind == "random":
            text_src, model, variant = random_sim(), "random", None
        elif emb_kind == "oracle":
            text_src, model, variant = oracle_sim(pos, len(types)), "oracle (uses hidden links)", None
        elif emb_kind == "text_ceiling":
            text_src, model = text_ceiling_sim(pos, corpus.texts(variant)), "text ceiling (oracle over identical texts)"
        elif emb_kind == "none":
            text_src, model, variant = None, None, None
        fi = None
        if fscheme != "none":
            fi = fidx_cache.setdefault(fscheme, FacetIndex(facets, fscheme))
        for w in ws:
            run_id = f"{name}/{run['name']}" + (f"@w{w:g}" if len(ws) > 1 else "")
            t0 = time.time()
            scorer = Scorer(text_src, fi, w)
            want_rank = run.get("eyeball") and (w == run.get("eyeball_w", ws[len(ws) // 2]))
            metrics, rankings = evaluate(scorer, corpus, pos, pool=pool, facets=facets,
                                         keep_rankings=eyeball_q if want_rank else None)
            eval_s = time.time() - t0
            row = {"kind": "proxy", "run_id": run_id, "model": model, "device": dev, "variant": variant,
                   "facets": fscheme if fi else None, "w_emb": scorer.w, "pool": pool,
                   "embed_s": embed_s, "eval_s": eval_s, "note": run.get("note", ""), "metrics": metrics,
                   "config": run, "positives": pstats,
                   "date": dt.datetime.now().isoformat(timespec="seconds")}
            fn = os.path.join(outdir, run_id.split("/", 1)[1].replace("/", "_") + ".json")
            json.dump(row, open(fn, "w"), indent=1)
            m = metrics
            print(f"{run_id:45s} macroR@10={m['macro']['R@10']:.3f} MRR={m['macro']['MRR']:.3f} "
                  f"MS={m.get('MeasurementSet', {}).get('R@10', float('nan')):.3f} "
                  f"AS={m.get('AnalysisSet', {}).get('R@10', float('nan')):.3f} "
                  f"PS={m.get('PredictionSet', {}).get('R@10', float('nan')):.3f} "
                  f"lab@10={m['macro']['lab@10']:.2f} ({eval_s:.0f}s)", file=sys.stderr)
            if rankings:
                md = os.path.join(outdir, run_id.split("/", 1)[1].replace("/", "_").replace("@", "_") + "_rankings.md")
                write_rankings_md(md, f"Rankings: {run_id}", corpus, rankings, facets, onto,
                                  variant or "full")
                print(f"   rankings -> {md}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--only", nargs="+")
    ap.add_argument("--device")
    a = ap.parse_args()
    run_config(a.config, a.only, a.device)
    report.build()


if __name__ == "__main__":
    main()
