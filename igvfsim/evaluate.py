"""Proxy-label evaluation.

There are no similarity labels, so positives come from existing portal links:
  publication     Sets citing the same publication
  analysis_input  MeasurementSets that are co-inputs to one AnalysisSet; AnalysisSets sharing an input
  related         MeasurementSet.related_measurement_sets (multiome partners etc.)
  collection      Sets in the same collection
Collections named IGVF_catalog_* are release tags and are ignored. Groups larger than --max-group
(default 200, e.g. "PD single cell multiomics", 1,284 Sets) are dropped: they are batch labels, not similarity signals, and would swamp the metrics.

The links are *hidden* from the scorer: publications, collections and link fields never appear in
the description text or facets, and the "exclude already-linked" filter is off during evaluation.

Metrics (per Set type, macro over queries; a query is any Set with >=1 positive in its pool):
  R@10     |top10 ∩ P| / min(10, |P|)    (capped, so large groups can still reach 1.0)
  MRR      1 / rank of the first positive
  nDCG@10  binary gains
  assay@10 / biosample@10   share of top-10 sharing an assay title / biosample term with the query
  lab@10   share of top-10 from the same lab: lab is NOT a feature, so a high value means
           batch-specific wording is driving the match (leakage diagnostic)
Ties (identical descriptions are common) are broken by a fixed random jitter, which gives the
expected metric under random tie-breaking instead of favoring index order.
"""
from __future__ import annotations

import json
import os
import re
from collections import defaultdict

import numpy as np

from .coverage import get_path
from .describe import Corpus, text
from .score import explain

SOURCES = ["publication", "analysis_input", "related", "collection"]
# Catalog-release tags group Sets by release, not by biology.
IGNORE_COLLECTIONS = re.compile(r"^IGVF_catalog", re.I)
SET_TYPES = ["MeasurementSet", "AnalysisSet", "PredictionSet"]


def build_positives(corpus: Corpus, max_group: int = 200) -> tuple[list[dict[int, set]], dict]:
    """pos[i] = {j: {sources}} over ALL types (filter by pool later). Also returns group stats."""
    acc = {a: i for i, a in enumerate(corpus.accessions)}
    n = len(acc)
    pos = [defaultdict(set) for _ in range(n)]
    stats = {s: {"groups": 0, "dropped_groups": 0, "pairs": 0} for s in SOURCES}

    def add_group(members, source):
        members = sorted(set(members))
        if len(members) < 2:
            return
        if len(members) > max_group:
            stats[source]["dropped_groups"] += 1
            return
        stats[source]["groups"] += 1
        for i in members:
            for j in members:
                if i != j:
                    pos[i][j].add(source)
                    stats[source]["pairs"] += 1

    groups = defaultdict(list)
    for i, o in enumerate(corpus.objs):
        for p in get_path(o, "publications.@id"):
            groups[("publication", p)].append(i)
        for c in get_path(o, "collections"):
            if not IGNORE_COLLECTIONS.match(c):
                groups[("collection", c)].append(i)
        if o["@type"][0] == "AnalysisSet":
            inputs = [acc[a] for a in get_path(o, "input_file_sets.accession") if a in acc]
            groups[("analysis_input", "AS:" + o["accession"])] += [j for j in inputs if corpus.types[j] == "MeasurementSet"]
            for a in get_path(o, "input_file_sets.accession"):  # AS sharing an input (any type)
                groups[("analysis_input", "IN:" + a)].append(i)
        if o["@type"][0] == "MeasurementSet":
            for a in get_path(o, "related_measurement_sets.measurement_sets.accession"):
                if a in acc:
                    groups[("related", o["accession"] + "|" + a)] += [i, acc[a]]
    for (source, _), members in groups.items():
        add_group(members, source)
    return pos, stats


def _agree(fa, fb, f):
    return bool(set(fa.get(f, {})) & set(fb.get(f, {})))


def evaluate(scorer, corpus: Corpus, pos, *, pool: str = "within", k: int = 10, block: int = 512,
             seed: int = 0, n_boot: int = 1000, sources=SOURCES, keep_rankings: list[int] | None = None,
             facets=None):
    """Returns (metrics dict, rankings dict for keep_rankings query indices)."""
    rng = np.random.default_rng(seed)
    types = np.array(corpus.types)
    labs = [o.get("lab", {}).get("title") for o in corpus.objs]
    facets = facets or corpus.facets()
    keep = set(keep_rankings or [])
    out, rankings = {}, {}
    for t in SET_TYPES:
        cand = np.where(types == t)[0] if pool == "within" else np.arange(len(types))
        cpos = {int(c): r for r, c in enumerate(cand)}
        qs = [i for i in np.where(types == t)[0] if any(j in cpos for j in pos[i] if j != i)]
        if not qs:
            continue
        per = defaultdict(list)
        per_src = {s: defaultdict(list) for s in sources}
        for b0 in range(0, len(qs), block):
            qb = np.array(qs[b0:b0 + block])
            S = scorer.block(qb, cand)
            S = S + rng.uniform(0, 1e-6, S.shape).astype(np.float32)  # random tie-breaking
            for r, i in enumerate(qb):
                if i in cpos:
                    S[r, cpos[i]] = -np.inf
            top = np.argpartition(-S, k, axis=1)[:, :k]
            for r, i in enumerate(qb):
                row = S[r]
                tk = top[r][np.argsort(-row[top[r]])]
                tk_idx = cand[tk]
                P = {j: s for j, s in pos[i].items() if j in cpos}
                pr = np.array([cpos[j] for j in P])
                best = row[pr].max()
                first = int((row > best).sum()) + 1
                hits = np.array([j in P for j in tk_idx], float)
                ideal = min(k, len(P))
                per["R@10"].append(hits.sum() / ideal)
                per["MRR"].append(1.0 / first)
                dcg = (hits / np.log2(np.arange(2, k + 2))).sum()
                per["nDCG@10"].append(dcg / (1 / np.log2(np.arange(2, ideal + 2))).sum())
                per["assay@10"].append(np.mean([_agree(facets[i], facets[j], "assay") for j in tk_idx]))
                per["biosample@10"].append(np.mean([_agree(facets[i], facets[j], "biosample") for j in tk_idx]))
                per["lab@10"].append(np.mean([labs[i] == labs[j] for j in tk_idx]))
                for s in sources:
                    Ps = {j for j, ss in P.items() if s in ss}
                    if Ps:
                        prs = np.array([cpos[j] for j in Ps])
                        per_src[s]["R@10"].append(sum(j in Ps for j in tk_idx) / min(k, len(Ps)))
                        per_src[s]["MRR"].append(1.0 / (int((row > row[prs].max()).sum()) + 1))
                if i in keep:
                    rankings[int(i)] = [(int(j), float(row[c]), sorted(P.get(int(j), [])))
                                        for j, c in zip(tk_idx, tk)]
        m = {"n_queries": len(qs), "pool": int(len(cand))}
        for key, vals in per.items():
            m[key] = float(np.mean(vals))
        for key in ("R@10", "MRR"):
            v = np.array(per[key])
            boots = v[rng.integers(0, len(v), (n_boot, len(v)))].mean(1)
            m[key + "_ci"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
        m["by_source"] = {s: {"n": len(d["R@10"]), "R@10": float(np.mean(d["R@10"])),
                              "MRR": float(np.mean(d["MRR"]))} for s, d in per_src.items() if d["R@10"]}
        out[t] = m
    present = [t for t in SET_TYPES if t in out]
    out["macro"] = {key: float(np.mean([out[t][key] for t in present]))
                    for key in ("R@10", "MRR", "nDCG@10", "assay@10", "biosample@10", "lab@10")}
    return out, rankings


def write_rankings_md(path: str, title: str, corpus: Corpus, rankings: dict, facets, onto=None,
                      variant: str = "full"):
    labs = [o.get("lab", {}).get("title") for o in corpus.objs]
    lines = [f"# {title}", "",
             "Top-10 per query. ✓ = proxy positive (source in brackets). Text is the description "
             "the model saw (truncated).", ""]
    for i, rows in rankings.items():
        p = corpus.parsed[i]
        lines += [f"## {p.accession} ({p.type}; {labs[i]})", "",
                  f"> {text(p, variant)[:400]}", "",
                  "| # | accession | score | proxy | why | lab | description |", "|---|---|---:|---|---|---|---|"]
        for r, (j, s, src) in enumerate(rows, 1):
            why = " · ".join(explain(facets[i], facets[j], onto)) or "—"
            mark = "✓ " + ",".join(src) if src else ""
            d = text(corpus.parsed[j], variant)[:140].replace("|", "/")
            lines.append(f"| {r} | {corpus.accessions[j]} | {s:.3f} | {mark} | {why} | {labs[j]} | {d} |")
        lines.append("")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write("\n".join(lines))


def main():
    """Print proxy-positive statistics."""
    import argparse
    from .describe import load_corpus
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-group", type=int, default=200)
    a = ap.parse_args()
    c = load_corpus()
    pos, stats = build_positives(c, a.max_group)
    types = c.types
    print(json.dumps(stats, indent=1))
    for t in SET_TYPES:
        q = [i for i, tt in enumerate(types) if tt == t]
        w = [i for i in q if any(types[j] == t for j in pos[i])]
        sizes = [sum(types[j] == t for j in pos[i]) for i in w]
        print(f"{t}: {len(w)}/{len(q)} have within-type positives; median |P|={np.median(sizes):.0f}, "
              f"p90={np.percentile(sizes, 90):.0f}")
        for s in SOURCES:
            ws = [i for i in q if any(types[j] == t and s in ss for j, ss in pos[i].items())]
            print(f"   {s:15s} {len(ws)} queries")


if __name__ == "__main__":
    main()
