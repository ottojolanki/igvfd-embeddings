"""Rebuild EXPERIMENTS.md from results/**/*.json (+ reports/findings.md, hand-written notes).

    python -m igvfsim.report
"""
from __future__ import annotations

import glob
import json
import os

OUT = "EXPERIMENTS.md"
FINDINGS = "reports/findings.md"
CONFIG_ORDER = ["baselines", "small", "variants", "cross", "large"]

PROXY_INTRO = """## A. Similar-Sets benchmark (proxy links)

Positives = Sets sharing a publication, co-inputs of an AnalysisSet, `related_measurement_sets`, or a
non-release collection (groups >200 dropped); the links are hidden from the scorer
(`igvfsim/evaluate.py`). Within-type pools unless `pool=all`. R@10 is capped recall
(hits / min(10, |P|)); MS/AS/PS = MeasurementSet/AnalysisSet/PredictionSet; macro = mean over types.
assay/bios/lab@10 = share of top-10 with the same assay title / biosample / lab as the query (lab
is not a feature; compare against the oracle's 0.75). Reference rows: `random` (floor), `oracle`
(1.0 by construction), `ceiling-*` (best any text-only scorer can do: identical descriptions are
indistinguishable). `related` = R@10 on multiome-partner positives only (same sample, other assay).
Embed s = wall-clock to encode the unique texts on the device shown (cached reruns cost ~0).
"""
PROXY_HEAD = ("| run | model | dev | variant | facets | w_emb | pool | macro R@10 | macro MRR | MS R@10 [95% CI] "
              "| AS R@10 [95% CI] | PS R@10 | MS related | assay@10 | bios@10 | lab@10 | embed s | eval s | note |\n"
              "|---|---|---|---|---|---:|---|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---|")
SEARCH_INTRO = """## B. Free-text search benchmark

30 hand-written queries (`configs/search_queries.yaml`) phrased as a user would type them; a document
(unique description text, all Set types, 4,297 docs) is relevant if it matches the query's structured
filters. Run with `python -m igvfsim.search_eval`.
"""
SEARCH_HEAD = ("| run | model | dev | variant | P@10 | nDCG@10 [95% CI] | R@100 | MRR |\n|---|---|---|---|---:|---|---:|---:|")


def boot_ci(vals, n=2000, seed=0):
    import numpy as np
    v = np.array(vals, float)
    b = v[np.random.default_rng(seed).integers(0, len(v), (n, len(v)))].mean(1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def f(x, nd=3):
    return "—" if x is None else f"{x:.{nd}f}" if isinstance(x, float) else str(x)


def proxy_row(r):
    m = r["metrics"]
    ms, as_, ps = (m.get(t, {}) for t in ("MeasurementSet", "AnalysisSet", "PredictionSet"))
    ci = lambda d: f"{f(d['R@10'])} [{f(d['R@10_ci'][0])}–{f(d['R@10_ci'][1])}]" if d else "—"
    rel = ms.get("by_source", {}).get("related", {}).get("R@10") if ms else None
    cells = [r["run_id"], r.get("model") or "—", r.get("device") or "—", r.get("variant") or "—",
             r.get("facets") or "—", f(float(r["w_emb"]), 2), r["pool"], f(m["macro"]["R@10"]),
             f(m["macro"]["MRR"]), ci(ms), ci(as_), f(ps.get("R@10")), f(rel),
             f(m["macro"]["assay@10"], 2), f(m["macro"]["biosample@10"], 2), f(m["macro"]["lab@10"], 2),
             r.get("embed_s", "—"), f(float(r["eval_s"]), 0), r.get("note", "")]
    return "| " + " | ".join(str(c).replace("|", "/") for c in cells) + " |"


THROUGHPUT_INTRO = """## C. Embedding throughput

Summed over every cache shard each model wrote (`data/emb_cache/*/shard_*.npz` metadata), so it covers
all variants and both machines. texts/s is wall-clock encode time only (excluding model load).
"""


def throughput_rows():
    import numpy as np
    agg = {}
    for fn in glob.glob("data/emb_cache/*/shard_*.npz"):
        with np.load(fn) as z:
            meta, n = json.loads(str(z["meta"])), len(z["keys"])
        if n < 64:  # query-only shards: too small to time meaningfully
            continue
        host = "mac" if "macOS" in meta.get("host", "") or "Darwin" in meta.get("host", "") else "sandbox"
        k = (meta["model"], host, meta.get("device", "?"), meta.get("dtype", "fp32" if meta.get("device") == "cpu" else "fp16"),
             meta.get("batch_size", "?"))
        a = agg.setdefault(k, [0, 0.0])
        a[0] += n
        a[1] += meta.get("seconds", 0.0)
    rows = ["| model | machine | device | dtype | batch | texts | seconds | texts/s | min per 4.3k texts |",
            "|---|---|---|---|---:|---:|---:|---:|---:|"]
    for (model, host, dev, dt, bs), (n, sec) in sorted(agg.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        rate = n / sec if sec else float("nan")
        rows.append(f"| {model} | {host} | {dev} | {dt} | {bs} | {n} | {sec:.0f} | {rate:.1f} | {4297 / rate / 60:.1f} |")
    return rows


def build():
    rows = [json.load(open(p)) for p in glob.glob("results/*/*.json")]
    proxy = [r for r in rows if r.get("kind") != "search"]
    search = sorted([r for r in rows if r.get("kind") == "search"], key=lambda r: -r["metrics"]["nDCG@10"])
    key = lambda r: (CONFIG_ORDER.index(r["run_id"].split("/")[0]) if r["run_id"].split("/")[0] in CONFIG_ORDER else 99,
                     r.get("date", ""), r["run_id"])
    proxy.sort(key=key)
    parts = ["# Experiments", ""]
    if os.path.exists(FINDINGS):
        parts += [open(FINDINGS).read().strip(), ""]
    parts += [SEARCH_INTRO, SEARCH_HEAD]
    for r in search:
        m = r["metrics"]
        lo, hi = boot_ci([v["nDCG@10"] for v in r["per_query"].values()])
        parts.append(f"| {r['run_id']} | {r['model']} | {r.get('device') or '—'} | {r.get('variant', 'full')} | {f(m['P@10'])} | "
                     f"{f(m['nDCG@10'])} [{f(lo, 2)}–{f(hi, 2)}] | {f(m['R@100'])} | {f(m['MRR'])} |")
    parts += ["", PROXY_INTRO, PROXY_HEAD] + [proxy_row(r) for r in proxy]
    parts += ["", THROUGHPUT_INTRO] + throughput_rows()
    open(OUT, "w").write("\n".join(parts) + "\n")
    print(f"wrote {OUT}: {len(proxy)} proxy rows, {len(search)} search rows")


if __name__ == "__main__":
    build()
