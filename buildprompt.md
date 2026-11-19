# Goal
Draft a "similar experiments" recommender and semantic search for the IGVF data portal
(https://data.igvf.org). For each experiment page, I want a ranked list of similar experiments
with a short, readable reason for each match ("Same assay: STARR-seq · Same biosample: K562").
Secondary goal: free-text semantic search over the same objects ("CRISPR screens in K562
measuring proliferation"). Use an open-source embedding model, self-hosted. No paid APIs.

# Data
- Portal API: https://api.data.igvf.org. Objects are JSON (snovault/encoded-style backend).
  - Search: `/search/?type=MeasurementSet&status=released&limit=all&format=json`
  - Add `&field=<dotted.path>` to return only the fields you need from embedded objects.
  - Single object: `/measurement-sets/<ACCESSION>/?format=json` (embedded frame).
- Start with `MeasurementSet` (experiments). Later: AnalysisSet, PredictionSet,
  ConstructLibrarySet.
- Schema reference: https://github.com/IGVF-DACC/igvfd
  - `src/igvfd/schemas/*.json` has the field definitions.
  - `src/igvfd/types/file_set.py` shows which linked objects are embedded and the calculated
    fields (`summary`, `assay_titles`, `preferred_assay_titles`, `assay_slims`,
    `related_measurement_sets`).
  - Ontology terms carry `ancestors`, `organ_slims`, `cell_slims`, `system_slims`, `assay_slims`.
  - `src/igvfd/tests/data/inserts/` has small fixture data for offline tests.
- First, fetch a few real objects and check which fields are actually populated.
  Don't assume the schema matches what's released.

# Starting point
`igvf_similar.py` in this repo is an untested-against-live-data prototype. It has only been run
on the igvfd test inserts with a TF-IDF stand-in, never on live portal data or a real embedding
model. It:
1. Turns each object into a curated text description: assay, biosample, organism, targeted
   genes, perturbation/readout, summary, description. It leaves out IDs, dates, files, lab and
   award, so "similar" means similar biology rather than same submitter.
2. Extracts structured facets: sets of normalized values per field.
3. Scores pairs as w_emb * cosine(embeddings) + (1 - w_emb) * weighted facet Jaccard.
4. Excludes already-linked experiments (multiome partners, controls), since those are
   "related", not "similar".
5. Writes the top-k matches per accession with "why" explanations, plus an .npz index for
   search.

Treat this design as a hypothesis to test, not a spec. Verify its `PORTAL_FIELDS` list against
the live API; it was written from the schema and never run.

# What I want you to do
1. **Data.** Pull all released MeasurementSets, cache them locally as JSON so experiments are
   reproducible, and report field coverage (% of objects with each field populated).
2. **Text descriptions.** Review the generated descriptions on ~20 diverse real objects and
   improve them. Write them to a JSONL file I can read.
3. **Evaluation set.** There are no labels, so build proxy positives from existing links:
   - same publication
   - shared inputs to an AnalysisSet
   - `related_measurement_sets`
   - same series/collection
   For these, hide the link itself when scoring. Report recall@10, MRR, and
   assay/biosample agreement of the top-10. Also produce a small sample of rankings for me to
   eyeball.
4. **Experiments.** Compare:
   - embedding models: start with Qwen/Qwen3-Embedding-0.6B and BAAI/bge-small-en-v1.5; try
     larger ones if the hardware allows
   - embeddings only vs. facets only vs. hybrid, sweeping w_emb
   - text description variants (with/without summary/description, field order)
   - ontology-aware soft matching using term ancestors/slims
   Keep a results table in `EXPERIMENTS.md` with config, metrics, runtime, and a short note.

# How to work
- Use a Python virtualenv. Pin dependencies in requirements.txt.
- Keep the code modular: load / describe / embed / score / evaluate / search. Each step should
  run from the CLI.
- No writing to igvf portal
- At most 10 requests per second to igvf portal. 
- Cache embeddings keyed by model name + text hash so reruns are cheap.
- Check in with me before:
  - changing the core definition of "similar"
  - adding infrastructure (vector DB, GPU services)
  - anything that would write to the portal or the igvfd repo
- Tell me what hardware you're running on and whether a GPU is available before you pick models.
