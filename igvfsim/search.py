"""Free-text semantic search and per-accession "similar Sets" lookup.

    python -m igvfsim.search "CRISPR screens in K562 measuring proliferation"
    python -m igvfsim.search "brain single-nucleus multiome" --type MeasurementSet --k 5
    python -m igvfsim.search --similar IGVFDS7615QNBX --w-emb 0.5 --facets exact

Embeddings come from the cache (data/emb_cache); missing ones are computed on the fly.
--similar excludes Sets already linked to the query (related/controls/shared inputs), because
the portal already shows those; they are "related", not "similar".
"""
from __future__ import annotations

import argparse

import numpy as np

from . import ontology
from .coverage import get_path
from .describe import load_corpus
from .embed import embed_query, embed_texts
from .score import FacetIndex, Scorer, explain


def linked(o: dict) -> set[str]:
    return set(get_path(o, "related_measurement_sets.measurement_sets.accession")
               + get_path(o, "control_file_sets.accession") + get_path(o, "input_for.accession")
               + get_path(o, "input_file_sets.accession"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("query", nargs="?")
    ap.add_argument("--similar", help="accession: list similar Sets instead of text search")
    ap.add_argument("--model", default="Qwen/Qwen3-Embedding-0.6B")
    ap.add_argument("--variant", default="full")
    ap.add_argument("--type", help="restrict results to one Set type")
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--w-emb", type=float, default=0.5)
    ap.add_argument("--facets", default="exact")
    ap.add_argument("--device")
    a = ap.parse_args()
    c = load_corpus()
    texts = c.texts(a.variant)
    emb, _ = embed_texts(texts, a.model, device=a.device)
    types = np.array(c.types)
    if a.similar:
        i = c.accessions.index(a.similar)
        t = a.type or c.types[i]
        cand = np.where(types == t)[0]
        facets = c.facets()
        s = Scorer(emb, FacetIndex(facets, a.facets), a.w_emb).block(np.array([i]), cand)[0]
        skip = linked(c.objs[i]) | {a.similar}
        print(f"{a.similar}: {texts[i][:200]}\n")
        onto, shown = ontology.get(), 0
        for r in np.argsort(-s):
            j = cand[r]
            if c.accessions[j] in skip:
                continue
            why = " · ".join(explain(facets[i], facets[j], onto))
            print(f"{s[r]:.3f}  {c.accessions[j]}  {why}\n        {texts[j][:150]}")
            shown += 1
            if shown == a.k:
                break
        return
    if not a.query:
        ap.error("give a query or --similar ACCESSION")
    q = embed_query(a.query, a.model, a.device)
    sims = emb @ q
    if a.type:
        sims = np.where(types == a.type, sims, -np.inf)
    seen = set()
    shown = 0
    for j in np.argsort(-sims):
        if texts[j] in seen:  # collapse identical descriptions
            continue
        seen.add(texts[j])
        dup = sum(1 for t in texts if t == texts[j])
        print(f"{sims[j]:.3f}  {c.accessions[j]} ({c.types[j]}{f', +{dup - 1} identical' if dup > 1 else ''})"
              f"\n        {texts[j][:180]}")
        shown += 1
        if shown == a.k:
            break


if __name__ == "__main__":
    main()
