"""Stratified review sample: round-robin over (type, assay title, organism) strata so rare assays
are represented, not just the big batches.

    python -m igvfsim.sample --n 60 --variant full    # -> data/review_sample.jsonl
"""
from __future__ import annotations

import argparse
import json
import os
import random
from collections import defaultdict

from .describe import load_corpus, text
from .load import DATA


def stratified_indices(parsed, n: int, seed: int = 0, per_type: dict | None = None) -> list[int]:
    rng = random.Random(seed)
    by_type = defaultdict(lambda: defaultdict(list))
    for i, p in enumerate(parsed):
        by_type[p.type][(tuple(p.assay_titles), tuple(p.taxa))].append(i)
    types = sorted(by_type)
    per_type = per_type or {t: n // len(types) for t in types}
    out = []
    for t in types:
        strata = list(by_type[t].values())
        for s in strata:
            rng.shuffle(s)
        strata.sort(key=len, reverse=True)
        picked, r = [], 0
        while len(picked) < per_type[t] and any(r < len(s) for s in strata):
            for s in strata:
                if r < len(s) and len(picked) < per_type[t]:
                    picked.append(s[r])
            r += 1
        out += picked
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--variant", default="full")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(DATA, "review_sample.jsonl"))
    a = ap.parse_args()
    c = load_corpus()
    idx = stratified_indices(c.parsed, a.n, a.seed)
    with open(a.out, "w") as f:
        for i in idx:
            p = c.parsed[i]
            f.write(json.dumps({"accession": p.accession, "type": p.type,
                                "assay": p.assay_titles, "lab": c.objs[i]["lab"]["title"],
                                "text": text(p, a.variant)}) + "\n")
    print(f"wrote {len(idx)} -> {a.out}")


if __name__ == "__main__":
    main()
