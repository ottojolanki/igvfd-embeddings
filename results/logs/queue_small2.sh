#!/bin/bash
set -e
PY=.venv/bin/python
$PY -m igvfsim.embed --config configs/small_models.yaml
$PY -m igvfsim.embed --config configs/variants.yaml
echo QUEUE_DONE
