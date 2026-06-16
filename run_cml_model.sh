#!/usr/bin/env bash
set -euo pipefail

workflow="${1:-figures}"

docker build -t cml-model:paper .

docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper "$workflow"
