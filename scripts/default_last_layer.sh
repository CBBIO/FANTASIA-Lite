#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_FASTA="${1:-$ROOT_DIR/fasta_test/PRUB1_longiso.pep}"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_CMD="$PYTHON_BIN"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_CMD="python3"
else
  PYTHON_CMD="python"
fi

export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

cd "$ROOT_DIR"

exec "$PYTHON_CMD" src/fantasia_pipeline.py \
  --use-gpu-lookup \
  --serial-models \
  --embed-models prot_t5 \
  --sequence-queue-package 100 \
  --embed-batch-size 4 \
  --model-batch-sizes prot_t5=4 \
  --lookup-batch-size 1024 \
  "$INPUT_FASTA"
