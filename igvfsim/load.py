"""Pull released Sets (+ ontology terms, publications) from the portal and cache them as JSON.

    python -m igvfsim.load                       # all default types
    python -m igvfsim.load --types PredictionSet --refresh

Raw pages go to data/raw/<Type>/page_NNNN.json (resumable: existing pages are not refetched);
the merged objects go to data/<Type>.jsonl; data/manifest.json records URLs, totals, timestamps.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys

from .portal import API, Portal

DATA = os.environ.get("IGVFSIM_DATA", "data")

_SAMPLE = [
    "samples.@id", "samples.@type", "samples.summary", "samples.classifications",
    "samples.sample_terms.@id", "samples.sample_terms.term_name",
    "samples.targeted_sample_term.@id", "samples.targeted_sample_term.term_name",
    "samples.disease_terms.term_name", "samples.treatments.treatment_term_name",
    "samples.treatments.summary", "samples.modifications.modality",
    "samples.modifications.summary", "samples.construct_library_sets",
    "samples.taxa", "samples.sorted_from", "samples.time_post_change",
    "samples.cellular_sub_pool",
]
_COMMON = [
    "@id", "@type", "accession", "status", "summary", "description", "file_set_type",
    "lab.title", "award.component", "collections", "publications.@id",
    "donors.taxa", "donors.sex", "donors.strain_background",
    "targeted_genes.symbol", "release_timestamp",
] + _SAMPLE

FIELDS: dict[str, list[str]] = {
    "MeasurementSet": _COMMON + [
        "preferred_assay_titles", "assay_titles", "assay_term.@id", "assay_term.term_name",
        "assay_term.assay_slims", "crispr_screen_readout", "primary_crispr_screen_readout",
        "crispr_screen_biometric.term_name", "functional_assay_mechanisms.term_name",
        "sequencing_library_types", "library_construction_platform.term_name",
        "multiome_size", "control_types", "targeted_loci",
        "related_measurement_sets.measurement_sets.accession",
        "control_file_sets.accession", "input_for.accession", "auxiliary_sets.accession",
        "auxiliary_sets.file_set_type",
    ],
    "AnalysisSet": _COMMON + [
        "assay_titles", "preferred_assay_titles", "workflows.name", "files.content_type",
        "input_file_sets.accession", "input_file_sets.@type",
        "input_file_sets.preferred_assay_titles", "input_file_sets.file_set_type",
        "functional_assay_mechanisms.term_name", "input_for.accession",
    ],
    "PredictionSet": _COMMON + [
        "scope", "software.title", "assessed_genes.symbol", "associated_phenotypes.term_name",
        "files.content_type", "input_file_sets.accession", "input_file_sets.@type",
        "small_scale_gene_list.symbol", "large_scale_gene_list", "input_for.accession",
    ],
    "SampleTerm": ["@id", "term_id", "term_name", "ancestors", "organ_slims", "cell_slims",
                   "system_slims", "developmental_slims", "is_a.term_name", "synonyms", "definition"],
    "AssayTerm": ["@id", "term_id", "term_name", "ancestors", "assay_slims", "category_slims",
                  "objective_slims", "preferred_assay_titles", "synonyms", "definition"],
    "Publication": ["@id", "title", "publication_identifiers", "date_published"],
}
SORT = {"SampleTerm": "term_id", "AssayTerm": "term_id", "Publication": ""}
DEFAULT_TYPES = ["MeasurementSet", "AnalysisSet", "PredictionSet", "SampleTerm", "AssayTerm",
                 "Publication"]


def pull(portal: Portal, item_type: str, page_size: int, refresh: bool) -> dict:
    raw_dir = os.path.join(DATA, "raw", item_type)
    os.makedirs(raw_dir, exist_ok=True)
    fields, sort = FIELDS[item_type], SORT.get(item_type, "accession")
    objs, frm, page, total = [], 0, 0, None
    while total is None or frm < total:
        fn = os.path.join(raw_dir, f"page_{page:04d}.json")
        if os.path.exists(fn) and not refresh:
            d = json.load(open(fn))
        else:
            d = portal.search(item_type, fields, frm=frm, limit=page_size, sort=sort)
            d.pop("facets", None)  # facets are large and unused
            json.dump(d, open(fn, "w"))
        total = d["total"]
        objs += d["@graph"]
        print(f"  {item_type}: {len(objs)}/{total}", file=sys.stderr)
        if not d["@graph"]:
            break
        frm += page_size
        page += 1
    ids = {o["@id"] for o in objs}
    if len(ids) != total:
        # The portal's sorted paging is not always a stable total order (seen on AnalysisSet:
        # the same object on two pages, another on none). Fall back to one limit=all request.
        print(f"  {item_type}: paging gave {len(ids)} unique of {total}; refetching with limit=all",
              file=sys.stderr)
        fn = os.path.join(raw_dir, "all.json")
        if os.path.exists(fn) and not refresh:
            d = json.load(open(fn))
        else:
            d = portal.search(item_type, fields, limit="all", sort=sort)
            d.pop("facets", None)
            json.dump(d, open(fn, "w"))
        objs, total = d["@graph"], d["total"]
        ids = {o["@id"] for o in objs}
        if len(ids) != total:
            print(f"  WARNING {item_type}: still {len(ids)} unique of {total}", file=sys.stderr)
    with open(os.path.join(DATA, f"{item_type}.jsonl"), "w") as f:
        for o in {o["@id"]: o for o in objs}.values():
            f.write(json.dumps(o) + "\n")
    return {"total": total, "n_cached": len(ids), "pages": page, "fields": fields,
            "sort": sort, "query": f"{API}/search/?type={item_type}&status=released"}


def read(item_type: str) -> list[dict]:
    fn = os.path.join(DATA, f"{item_type}.jsonl")
    if not os.path.exists(fn):
        raise FileNotFoundError(f"{fn} missing; run `python -m igvfsim.load` first")
    return [json.loads(l) for l in open(fn)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--types", nargs="+", default=DEFAULT_TYPES, choices=list(FIELDS))
    ap.add_argument("--page-size", type=int, default=500)
    ap.add_argument("--rps", type=float, default=4.0)
    ap.add_argument("--refresh", action="store_true", help="refetch pages already on disk")
    a = ap.parse_args()
    portal = Portal(rps=a.rps)
    mpath = os.path.join(DATA, "manifest.json")
    manifest = json.load(open(mpath)) if os.path.exists(mpath) else {}
    for t in a.types:
        before = portal.n_requests
        info = pull(portal, t, a.page_size, a.refresh)
        if portal.n_requests > before or t not in manifest:
            info["pulled_at"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        else:
            info["pulled_at"] = manifest[t].get("pulled_at")
        manifest[t] = info
    json.dump(manifest, open(mpath, "w"), indent=1)
    print(f"done; {portal.n_requests} HTTP requests", file=sys.stderr)


if __name__ == "__main__":
    main()
