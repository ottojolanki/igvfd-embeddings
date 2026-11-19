"""Free-text search benchmark: hand-written queries, relevance from structured filters.

    python -m igvfsim.search_eval --model BAAI/bge-small-en-v1.5 --model Qwen/Qwen3-Embedding-0.6B
    python -m igvfsim.search_eval --lexical            # TF-IDF + BM25 baselines
    python -m igvfsim.search_eval --list               # show queries and relevant-set sizes

Documents are the *unique* description texts across all Set types (identical descriptions would
otherwise fill the top-10 with copies). A unique text is relevant if any Set carrying it matches all
filters of the query. Metrics: P@10, nDCG@10, R@100, MRR (macro over queries).
Results go to results/search/<model>.json; EXPERIMENTS.md is rebuilt by igvfsim.report.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os

import numpy as np
import yaml

from .describe import load_corpus, norm
from .embed import embed_texts, prompts_for, slug

FIELDS = {"type", "assay_titles", "assay_terms", "assay_slims", "sample_terms", "targeted_terms", "slims",
          "taxa", "genes", "readout", "modifications", "treatments", "phenotypes", "file_set_type", "workflows"}


def matches(p, rel: dict) -> bool:
    for f, vals in rel.items():
        assert f in FIELDS, f
        have = getattr(p, f)
        have = {norm(x) for x in (have if isinstance(have, list) else [have])}
        if not have & {norm(v) for v in vals}:
            return False
    return True


def load_benchmark(corpus, path="configs/search_queries.yaml", variant="full"):
    qs = yaml.safe_load(open(path))["queries"]
    texts = corpus.texts(variant)
    uniq = list(dict.fromkeys(texts))
    uix = {t: i for i, t in enumerate(uniq)}
    rel = []
    for q in qs:
        r = set()
        for p, t in zip(corpus.parsed, texts):
            if matches(p, q["rel"]):
                r.add(uix[t])
        rel.append(r)
    return qs, uniq, rel


def minmax(S):
    lo, hi = S.min(1, keepdims=True), S.max(1, keepdims=True)
    return (S - lo) / np.maximum(hi - lo, 1e-9)


def bm25_scores(docs, queries, k1=1.2, b=0.75):
    import scipy.sparse as sp
    from sklearn.feature_extraction.text import CountVectorizer
    cv = CountVectorizer(token_pattern=r"(?u)\b[\w+-]+\b")
    tf = cv.fit_transform(docs).astype(np.float32).tocsr()
    dl = np.asarray(tf.sum(1)).ravel()
    df = np.asarray((tf > 0).sum(0)).ravel()
    idf = np.log1p((len(docs) - df + 0.5) / (df + 0.5))
    tf = tf.tocoo()
    w = tf.data * (k1 + 1) / (tf.data + k1 * (1 - b + b * dl[tf.row] / dl.mean())) * idf[tf.col]
    W = sp.csr_matrix((w, (tf.row, tf.col)), shape=tf.shape)
    Q = (cv.transform(queries) > 0).astype(np.float32)
    return (Q @ W.T).toarray()


def metrics(S: np.ndarray, rel: list[set], k=10):
    out = {"P@10": [], "nDCG@10": [], "R@100": [], "MRR": []}
    for s, r in zip(S, rel):
        order = np.argsort(-s)
        hits = np.array([j in r for j in order[:100]], float)
        out["P@10"].append(hits[:k].mean())
        ideal = min(k, len(r))
        out["nDCG@10"].append((hits[:k] / np.log2(np.arange(2, k + 2))).sum()
                              / (1 / np.log2(np.arange(2, ideal + 2))).sum())
        out["R@100"].append(hits.sum() / len(r))
        rr = np.where([j in r for j in order])[0]
        out["MRR"].append(1.0 / (rr[0] + 1) if len(rr) else 0.0)
    per_query = {key: [float(x) for x in v] for key, v in out.items()}
    return {key: float(np.mean(v)) for key, v in out.items()}, per_query


def save(name, model, dev, m, per_query, qs, extra):
    os.makedirs("results/search", exist_ok=True)
    if extra.get("variant", "full") != "full":
        name = f"{name}[{extra['variant']}]"
    row = {"kind": "search", "run_id": f"search/{name}", "model": model, "device": dev, "metrics": m,
           "per_query": {q["q"]: {k: v[i] for k, v in per_query.items()} for i, q in enumerate(qs)},
           "date": dt.datetime.now().isoformat(timespec="seconds"), **extra}
    json.dump(row, open(f"results/search/{slug(name)}.json", "w"), indent=1)
    print(f"{name:40s} P@10={m['P@10']:.3f} nDCG@10={m['nDCG@10']:.3f} R@100={m['R@100']:.3f} MRR={m['MRR']:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", action="append", default=[])
    ap.add_argument("--variant", default="full")
    ap.add_argument("--queries", default="configs/search_queries.yaml")
    ap.add_argument("--device")
    ap.add_argument("--lexical", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--fuse", type=float, nargs="*", default=[], help="also score alpha*dense+(1-alpha)*bm25")
    a = ap.parse_args()
    corpus = load_corpus()
    qs, uniq, rel = load_benchmark(corpus, a.queries, a.variant)
    if a.list:
        for q, r in zip(qs, rel):
            print(f"{len(r):5d}  {q['q']}")
        print(f"{len(uniq)} unique documents")
        return
    for i, r in enumerate(rel):
        assert r, f"query has no relevant documents: {qs[i]['q']}"
    extra = {"variant": a.variant, "n_queries": len(qs), "n_docs": len(uniq)}
    bm = bm25_scores(uniq, [q["q"] for q in qs])
    if a.lexical:
        from sklearn.feature_extraction.text import TfidfVectorizer
        v = TfidfVectorizer(ngram_range=(1, 2), sublinear_tf=True)
        D = v.fit_transform(uniq)
        S = (v.transform([q["q"] for q in qs]) @ D.T).toarray()
        save("tfidf", "tfidf(1-2gram)", "cpu", *metrics(S, rel), qs, extra)
        save("bm25", "bm25", "cpu", *metrics(bm, rel), qs, extra)
    for model in a.model:
        emb, st = embed_texts(uniq, model, device=a.device)
        # queries go through the cache too (python -m igvfsim.embed --config caches them), so
        # models too large for this machine can be scored from embeddings made elsewhere
        Q, qst = embed_texts([q["q"] for q in qs], model, prompt=prompts_for(model)["query"], device=a.device)
        S = Q @ emb.T
        dev = st["device"] or st["original"]["device"]
        save(model, model, dev, *metrics(S, rel), qs, extra)
        for alpha in a.fuse:
            # hybrid lexical+dense: per-query min-max normalized scores, alpha * dense + (1-alpha) * bm25
            save(f"{model}+bm25@{alpha:g}", f"{model} + BM25 (alpha={alpha:g})", dev,
                 *metrics(alpha * minmax(S) + (1 - alpha) * minmax(bm), rel), qs, extra)


if __name__ == "__main__":
    main()
