# igvfd-embeddings

"Similar Sets" recommender + semantic search experiments for the IGVF data portal
(https://data.igvf.org), using self-hosted open-source embedding models.
`igvfsim/` is the modular pipeline; `igvf_similar.py` is the original prototype (kept for reference).

## Setup

```bash
uv venv --python 3.12 .venv            # sandbox (Linux).  On the Mac host use .venv-host instead:
uv pip install --python .venv/bin/python -r requirements.txt
```

## Pipeline (each step is a CLI)

| step | command | output |
|---|---|---|
| load | `python -m igvfsim.load` | `data/raw/<Type>/page_*.json`, `data/<Type>.jsonl`, `data/manifest.json` |
| coverage | `python -m igvfsim.coverage` | `reports/coverage.md` |
| describe | `python -m igvfsim.describe --all-variants` | `data/descriptions/<variant>.jsonl` |
| review sample | `python -m igvfsim.sample --n 60` | `data/review_sample.jsonl` |
| embed | `python -m igvfsim.embed --model BAAI/bge-small-en-v1.5 --variant full` | `data/emb_cache/<model>/shard_*.npz` |
| proxy stats | `python -m igvfsim.evaluate` | stdout |
| experiments | `python -m igvfsim.experiments configs/baselines.yaml` | `EXPERIMENTS.md`, `results/<config>/*.json`, `*_rankings.md` |
| search benchmark | `python -m igvfsim.search_eval --lexical --model BAAI/bge-small-en-v1.5 --fuse 0.5` | `results/search/*.json` |
| report | `python -m igvfsim.report` | `EXPERIMENTS.md` (rebuilt from `results/` + `reports/findings.md`) |
| search | `python -m igvfsim.search "CRISPR screens in K562 measuring proliferation"` | stdout |
| similar | `python -m igvfsim.search --similar IGVFDS7615QNBX` | stdout |

`load` is read-only (GET only), rate-limited (default 4 req/s, hard cap 10), and resumable:
pages already on disk are not refetched (`--refresh` to refetch). Released MS/AS/PS + SampleTerm,
AssayTerm and Publication are about 50 requests in total.

## Hardware: sandbox CPU + host GPU

The sandbox (Linux aarch64, 16 cores, 31 GB, no GPU) runs everything, including models up to
~0.6B params on CPU (Qwen3-Embedding-0.6B: ~6 texts/s, ~12 min for the ~4.3k unique texts of one
variant; bge-small: ~70 texts/s). Larger models run on the Mac host with MPS. Both write to the same
content-addressed cache in `data/emb_cache/` (key = sha1(model | prompt | max_len | text)), so the
evaluation picks up host-computed embeddings automatically.

On the host (macOS, Apple silicon), from the repo directory:

```bash
# 0. one-time setup (separate venv: the sandbox's .venv holds Linux binaries)
uv venv --python 3.12 .venv-host
uv pip install --python .venv-host/bin/python -r requirements.txt

# 1. speed check (~5 min incl. downloads of ~25 GB on first run; writes nothing to the cache)
.venv-host/bin/python -m igvfsim.embed --benchmark 512 --device mps \
    --model BAAI/bge-small-en-v1.5 --model Qwen/Qwen3-Embedding-0.6B \
    --model Qwen/Qwen3-Embedding-4B --model Qwen/Qwen3-Embedding-8B

# 2. full embedding run for every model/variant in configs/large.yaml (+ the 30 search queries)
.venv-host/bin/python -m igvfsim.embed --config configs/large.yaml --device mps 2>&1 | tee results/logs/embed_large_host.log

# 3. evaluate (either machine; reads only the cache)
python -m igvfsim.experiments configs/large.yaml
python -m igvfsim.search_eval --model Qwen/Qwen3-Embedding-4B --model Qwen/Qwen3-Embedding-8B \
    --model BAAI/bge-large-en-v1.5 --fuse 0.5
```

Notes: models load in fp16 on MPS; if a chunk comes out NaN/zero (fp16 overflow) it is redone in fp32
automatically and never cached (`--dtype bf16|fp32` forces a dtype). Runs checkpoint every 256 texts,
so Ctrl-C and rerun resumes. Section C of `EXPERIMENTS.md` (`python -m igvfsim.report`) tabulates
texts/s per model, machine, device and dtype from the cache metadata.

## Evaluation in one paragraph

No labels exist, so positives are existing portal links: same publication, co-inputs to an
AnalysisSet, `related_measurement_sets`, same collection (catalog-release collections and groups >200
are ignored). The links never appear in the text or facets. Metrics are capped R@10, MRR, nDCG@10,
plus assay/biosample agreement of the top-10 and a same-lab rate. Two reference rows matter when
reading results: `oracle` (must be 1.0) and `text_ceiling` (the best any text-only model can do,
because many Sets have identical descriptions). A second benchmark, 30 hand-written free-text
queries whose relevance is defined by structured filters (`configs/search_queries.yaml`), covers the
search use case. See `EXPERIMENTS.md` for results and findings.
