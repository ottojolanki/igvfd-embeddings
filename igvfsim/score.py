"""Similarity scoring: embedding cosine, weighted facet Jaccard (exact or ontology-aware), hybrid.

score(i, j) = w_emb * cos(e_i, e_j) + (1 - w_emb) * facet_sim(i, j)

facet_sim is the prototype's weighted Jaccard: for each facet present in either object,
Jaccard(values_i, values_j) weighted by the facet weight, normalized by the total weight of
facets present. Values can carry weights (IDF) so that, e.g., sharing the ancestor
"dopaminergic neuron" counts more than sharing "eukaryotic cell".
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np
import scipy.sparse as sp

# Prototype weights (igvf_similar.FACET_WEIGHTS), with "disease" (0% populated) -> "phenotype".
SCHEMES: dict[str, dict[str, float]] = {
    "exact": {"assay": .30, "assay_term": .10, "assay_slim": .05, "biosample": .25,
              "biosample_slim": .05, "organism": .05, "genes": .10, "perturbation": .05,
              "phenotype": .05},
    # Ontology-aware: half of the biosample weight moves to IDF-weighted ancestor-closure overlap
    # (soft match: "putamen" ~ "caudate nucleus" via "striatum"), likewise for assay terms.
    "onto": {"assay": .30, "assay_term": .05, "assay_onto": .05, "assay_slim": .05,
             "biosample": .15, "biosample_onto": .15, "biosample_slim": .0, "organism": .05,
             "genes": .10, "perturbation": .05, "phenotype": .05},
    "onto_only": {"assay": .30, "assay_term": .05, "assay_onto": .05, "assay_slim": .05,
                  "biosample_onto": .30, "organism": .05, "genes": .10, "perturbation": .05,
                  "phenotype": .05},
}
IDF_FACETS = {"biosample_onto", "assay_onto"}


class FacetIndex:
    """Sparse value matrices per facet, for fast blockwise weighted-Jaccard."""

    def __init__(self, facets: list[dict[str, dict[str, str]]], scheme: str | dict = "exact"):
        self.weights = SCHEMES[scheme] if isinstance(scheme, str) else scheme
        self.facets = facets
        n = len(facets)
        self.X, self.wX, self.size = {}, {}, {}
        for f, w in self.weights.items():
            if w <= 0:
                continue
            vocab, rows, cols = {}, [], []
            for i, fs in enumerate(facets):
                for v in fs.get(f, {}):
                    rows.append(i)
                    cols.append(vocab.setdefault(v, len(vocab)))
            X = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(n, max(len(vocab), 1)))
            if f in IDF_FACETS:
                df = np.asarray(X.sum(0)).ravel()
                vw = np.log1p(n / np.maximum(df, 1)).astype(np.float32)
            else:
                vw = np.ones(X.shape[1], np.float32)
            self.X[f] = X
            self.wX[f] = (X @ sp.diags(vw)).tocsr()
            self.size[f] = np.asarray(self.wX[f].sum(1)).ravel()

    def sim(self, q: np.ndarray, c: np.ndarray) -> np.ndarray:
        """Dense (len(q), len(c)) facet similarity between rows q and candidate rows c."""
        num = np.zeros((len(q), len(c)), np.float32)
        den = np.zeros_like(num)
        for f, X in self.X.items():
            w = self.weights[f]
            inter = (self.wX[f][q] @ X[c].T).toarray()
            union = self.size[f][q][:, None] + self.size[f][c][None, :] - inter
            present = union > 0
            num += np.where(present, w * inter / np.maximum(union, 1e-9), 0)
            den += np.where(present, w, 0)
        return np.where(den > 0, num / np.maximum(den, 1e-9), 0).astype(np.float32)


class Scorer:
    """Combines a text similarity source and a facet index.

    text: dense L2-normalized embeddings (n, d), a sparse L2-normalized TF-IDF matrix,
          or a callable (q, c) -> dense block; None for facets-only.
    """

    def __init__(self, text=None, facet_index: FacetIndex | None = None, w_emb: float = 0.5):
        self.text, self.fi, self.w = text, facet_index, w_emb
        if self.text is None:
            self.w = 0.0
        if self.fi is None:
            self.w = 1.0

    def text_sim(self, q, c):
        if callable(self.text):
            return self.text(q, c)
        if sp.issparse(self.text):
            return (self.text[q] @ self.text[c].T).toarray().astype(np.float32)
        return self.text[q] @ self.text[c].T

    def block(self, q: np.ndarray, c: np.ndarray) -> np.ndarray:
        out = np.zeros((len(q), len(c)), np.float32)
        if self.w > 0:
            out += self.w * self.text_sim(q, c)
        if self.w < 1:
            out += (1 - self.w) * self.fi.sim(q, c)
        return out


# ----------------------------------------------------------------------------- BM25

def bm25_sim(texts: list[str], k1: float = 1.2, b: float = 0.75):
    """Returns a callable (q, c) -> BM25 scores using doc q's terms as the query, rescaled per
    query row to [0, 1] (so it is comparable across queries; used as a lexical baseline)."""
    from sklearn.feature_extraction.text import CountVectorizer
    cv = CountVectorizer(ngram_range=(1, 1), token_pattern=r"(?u)\b[\w+-]+\b")
    tf = cv.fit_transform(texts).tocsr().astype(np.float32)
    n = tf.shape[0]
    dl = np.asarray(tf.sum(1)).ravel()
    df = np.asarray((tf > 0).sum(0)).ravel()
    idf = np.log1p((n - df + 0.5) / (df + 0.5)).astype(np.float32)
    tf = tf.tocoo()
    denom = tf.data + k1 * (1 - b + b * dl[tf.row] / dl.mean())
    W = sp.csr_matrix((tf.data * (k1 + 1) / denom * idf[tf.col], (tf.row, tf.col)), shape=tf.shape)
    Q = (W > 0).astype(np.float32).tocsr()

    def f(q, c):
        s = (Q[q] @ W[c].T).toarray().astype(np.float32)
        return s / np.maximum(s.max(1, keepdims=True), 1e-9)
    return f


# ----------------------------------------------------------------------------- explanations

LABELS = [("assay", "Same assay"), ("biosample", "Same biosample"), ("genes", "Same targeted genes"),
          ("perturbation", "Same perturbation/readout"), ("phenotype", "Same phenotype"),
          ("organism", "Same organism"), ("assay_slim", "Same assay category"),
          ("biosample_slim", "Same tissue/cell type")]


def explain(fa: dict, fb: dict, onto=None, max_items: int = 4) -> list[str]:
    """Short human-readable reasons, most specific first.
    e.g. ["Same assay: STARR-seq", "Same biosample: K562"]."""
    out, have_bios = [], False
    for f, label in LABELS:
        shared = sorted(set(fa.get(f, {})) & set(fb.get(f, {})))
        if not shared:
            continue
        if f == "biosample":
            have_bios = True
        if f == "biosample_slim" and have_bios:
            continue
        if f == "assay_slim" and any(r.startswith("Same assay:") for r in out):
            continue
        if f == "organism" and len(out) >= 2:
            continue
        vals = [fa[f][v] for v in shared[:3]]
        out.append(f"{label}: {', '.join(vals)}" + (f" (+{len(shared) - 3})" if len(shared) > 3 else ""))
    if not have_bios and onto is not None:
        # soft biosample match through the ontology: name the most specific shared ancestor
        shared = set(fa.get("biosample_onto", {})) & set(fb.get("biosample_onto", {}))
        if shared:
            best = max(shared, key=lambda a: onto.idf.get(a, 0))
            out.append(f"Related biosample: both {fa['biosample_onto'][best]}")
    return out[:max_items]
