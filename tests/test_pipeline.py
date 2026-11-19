"""Offline unit tests (no network, no model downloads). Run: python -m pytest -q tests"""
import numpy as np
import pytest

from igvfsim import portal
from igvfsim.describe import scrub_free_text, scrub_sample_summary
from igvfsim.score import FacetIndex, Scorer, explain


def test_portal_is_get_only_and_capped():
    p = portal.Portal(rps=1000)
    assert p.min_interval >= 1 / portal.HARD_MAX_RPS
    assert not any(hasattr(p, m) for m in ("post", "put", "patch", "delete"))


@pytest.mark.parametrize("raw,expected", [
    ("source pool6:gex_ox_il_2; export target gex_ox_il_2", "source; export target"),
    ("donor_SV20_day45_BR_batch2; quantify labeled RNA", "quantify labeled RNA"),
    ("data for ENCSR753ZRI were generated", "data for were generated"),
    ("CRISPRa Tet1-dCas9 screen in iPSCs", "CRISPRa Tet1-dCas9 screen in iPSCs"),  # biology kept
])
def test_scrub_free_text(raw, expected):
    assert scrub_free_text(raw) == expected


def test_scrub_sample_summary():
    s = "multiplexed putamen, 8 donors, 8 samples, cellular sub-pool: RP032_G3"
    assert scrub_sample_summary(s) == "multiplexed putamen, 8 donors, 8 samples"


def _f(**kw):
    base = {k: {} for k in ("assay", "assay_term", "assay_slim", "assay_onto", "biosample", "biosample_slim",
                            "biosample_onto", "organism", "genes", "perturbation", "phenotype", "workflow")}
    for k, vals in kw.items():
        base[k] = {v.lower(): v for v in vals}
    return base


def test_facet_sim_identity_and_disjoint():
    fs = [_f(assay=["STARR-seq"], biosample=["K562"]), _f(assay=["STARR-seq"], biosample=["K562"]),
          _f(assay=["MPRA"], biosample=["HepG2"])]
    S = FacetIndex(fs, "exact").sim(np.arange(3), np.arange(3))
    assert S[0, 1] == pytest.approx(1.0)
    assert S[0, 2] == pytest.approx(0.0)


def test_hybrid_weighting():
    emb = np.eye(3, dtype=np.float32)
    fs = [_f(assay=["A"]), _f(assay=["A"]), _f(assay=["B"])]
    s = Scorer(emb, FacetIndex(fs, "exact"), 0.25).block(np.array([0]), np.arange(3))[0]
    assert s[1] == pytest.approx(0.75)  # cos 0, facets 1
    assert s[0] == pytest.approx(1.0)


def test_explain():
    a = _f(assay=["STARR-seq"], biosample=["K562"], organism=["Homo sapiens"])
    assert explain(a, a)[:2] == ["Same assay: STARR-seq", "Same biosample: K562"]
