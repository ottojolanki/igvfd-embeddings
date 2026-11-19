"""Field coverage: % of objects with each (dotted) field populated, per type.

    python -m igvfsim.coverage            # -> reports/coverage.md
"""
from __future__ import annotations

import argparse
import json
import os
from collections import Counter

from .load import FIELDS, read

SET_TYPES = ["MeasurementSet", "AnalysisSet", "PredictionSet"]


def empty(v):
    return v is None or v == "" or v == [] or v == {}


def get_path(obj, path: str) -> list:
    """All non-empty values at a dotted path, flattening lists along the way."""
    cur = [obj]
    for p in path.split("."):
        nxt = []
        for c in cur:
            for x in (c if isinstance(c, list) else [c]):
                if isinstance(x, dict) and p in x:
                    v = x[p]
                    nxt += v if isinstance(v, list) else [v]
        cur = nxt
    return [c for c in cur if not empty(c)]


def leaf_paths(obj, prefix="") -> set[str]:
    out = set()
    for x in (obj if isinstance(obj, list) else [obj]):
        if not isinstance(x, dict):
            if prefix and not empty(x):
                out.add(prefix)
            continue
        for k, v in x.items():
            if k in ("audit",):
                continue
            p = f"{prefix}.{k}" if prefix else k
            if isinstance(v, (dict, list)) and v and any(isinstance(i, dict) for i in (v if isinstance(v, list) else [v])):
                out |= leaf_paths(v, p)
            elif not empty(v):
                out.add(p)
    return out


def coverage(objs: list[dict], requested: list[str]) -> list[tuple[str, float, str]]:
    seen = Counter()
    for o in objs:
        for p in leaf_paths(o):
            seen[p] += 1
    paths = sorted(set(requested) | set(seen))
    rows = []
    for p in paths:
        n = sum(1 for o in objs if get_path(o, p))
        ex = next((get_path(o, p) for o in objs if get_path(o, p)), None)
        rows.append((p, 100.0 * n / max(len(objs), 1), json.dumps(ex[:3])[:80] if ex else ""))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", nargs="+", default=SET_TYPES)
    ap.add_argument("--out", default="reports/coverage.md")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    manifest = json.load(open("data/manifest.json"))
    lines = ["# Field coverage of released Sets", "",
             "Share of released objects where each field (dotted path into the embedded object) "
             "is non-empty. Fields requested from the API but never populated show 0%.", ""]
    for t in a.types:
        objs = read(t)
        lines += [f"## {t} (n={len(objs)}, pulled {manifest.get(t, {}).get('pulled_at')})", "",
                  "| field | % populated | example |", "|---|---:|---|"]
        for p, pct, ex in coverage(objs, FIELDS[t]):
            if p.startswith("@") and p not in ("@id", "@type"):
                continue
            lines.append(f"| `{p}` | {pct:.1f} | {ex.replace('|', '/')} |")
        lines.append("")
    open(a.out, "w").write("\n".join(lines))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
