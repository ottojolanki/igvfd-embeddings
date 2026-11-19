#!/bin/bash
# Sequential CPU embedding queue (one model at a time so timings are clean).
set -e
PY=.venv/bin/python
while pgrep -f "igvfsim.embed --model Qwen/Qwen3-Embedding-0.6B --variant full$" >/dev/null; do sleep 10; done
$PY -m igvfsim.embed --model BAAI/bge-small-en-v1.5 --variant no_description --variant no_summary --variant structured_only --variant structured_first --variant free_text_only --variant full_with_batch_ids
$PY -m igvfsim.embed --model BAAI/bge-base-en-v1.5 --variant full
$PY -m igvfsim.embed --model NeuML/pubmedbert-base-embeddings --variant full
$PY -m igvfsim.embed --model nomic-ai/nomic-embed-text-v1.5 --variant full
$PY -m igvfsim.embed --config configs/small_models.yaml
$PY -m igvfsim.embed --config configs/variants.yaml
echo QUEUE_DONE
