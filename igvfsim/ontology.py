"""Ontology lookups: SampleTerm / AssayTerm ancestors and slims, joined by @id.

The portal's search frame does not embed slims/ancestors under samples.sample_terms
(0% coverage), so we pull the term collections separately (load.py) and join locally.
"""
from __future__ import annotations

import math
from collections import Counter
from functools import lru_cache

from .load import read

# Ancestors that carry no information (top of CL/UBERON/EFO); dropped before soft matching.
GENERIC = {
    "continuant", "independent continuant", "material entity", "anatomical entity",
    "anatomical structure", "material anatomical entity", "entity", "thing", "cell",
    "native cell", "experimental factor", "material property", "information entity",
    "cultured cell", "disposition", "quality", "process", "occurrent", "multicellular anatomical structure",
    "anatomical conduit", "anatomical system", "anatomy basic component", "cellular_component",
    "assay", "planned process", "completely executed planned process", "experimental process",
}


class Ontology:
    def __init__(self):
        self.sample = {t["@id"]: t for t in read("SampleTerm")}
        self.assay = {t["@id"]: t for t in read("AssayTerm")}
        anc = Counter()
        for t in self.sample.values():
            anc.update(set(self.closure(t["@id"])))
        n = max(len(self.sample), 1)
        # IDF over the term collection: rare ancestors (e.g. "dopaminergic neuron") count more
        # than common ones (e.g. "eukaryotic cell") in soft matching.
        self.idf = {a: math.log(1 + n / c) for a, c in anc.items()}

    @lru_cache(maxsize=None)
    def closure(self, term_id: str) -> tuple[str, ...]:
        """Term itself + ancestors, as lowercase names, generic roots removed."""
        t = self.sample.get(term_id) or self.assay.get(term_id)
        if not t:
            return ()
        names = [t.get("term_name", "")] + list(t.get("ancestors") or [])
        return tuple(dict.fromkeys(n.lower() for n in names if n and n.lower() not in GENERIC))

    def slims(self, term_id: str) -> list[str]:
        t = self.sample.get(term_id) or self.assay.get(term_id) or {}
        out = []
        for k in ("organ_slims", "cell_slims", "system_slims", "developmental_slims",
                  "assay_slims", "category_slims", "objective_slims"):
            out += t.get(k) or []
        return list(dict.fromkeys(out))


_ONT = None


def get() -> Ontology:
    global _ONT
    if _ONT is None:
        _ONT = Ontology()
    return _ONT
