#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_FASTA="${1:-$ROOT_DIR/fasta_test/PRUB1_longiso.pep}"
OUTPUT_PATH="${2:-$ROOT_DIR/outputs_layers.npz}"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_CMD="$PYTHON_BIN"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_CMD="python3"
else
  PYTHON_CMD="python"
fi

shift_count=0
if [[ $# -ge 1 ]]; then
  shift_count=$((shift_count + 1))
fi
if [[ $# -ge 2 ]]; then
  shift_count=$((shift_count + 1))
fi
if [[ $shift_count -gt 0 ]]; then
  shift "$shift_count"
fi

export PYTORCH_CUDA_ALLOC_CONF="${PYTORCH_CUDA_ALLOC_CONF:-expandable_segments:True}"

cd "$ROOT_DIR"

layer_args=()
if [[ $# -eq 0 ]]; then
  layer_args+=(--all-layers)
else
  layer_args+=(--layer-indices "$@")
fi

exec "$PYTHON_CMD" src/generate_embeddings.py \
  --fasta "$INPUT_FASTA" \
  --output "$OUTPUT_PATH" \
  --models prot_t5 \
  --queue-batch-size 100 \
  --embed-batch-size 4 \
  --model-batch-sizes prot_t5=4 \
  "${layer_args[@]}"
