#!/bin/bash
# Pipeline TP Bioinformatica - Linux / macOS
# Tambien funciona llamando: python run_pipeline.py
set -euo pipefail

cd "$(dirname "$0")"

if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
fi

export ENTREZ_EMAIL="${ENTREZ_EMAIL:-estudiante@itba.edu.ar}"
export BLAST_MODE="${BLAST_MODE:-remote}"   # remote | local | both
export REQUIRE_EMBOSS="${REQUIRE_EMBOSS:-1}"

echo "[INFO] ENTREZ_EMAIL=$ENTREZ_EMAIL"
echo "[INFO] BLAST_MODE=$BLAST_MODE"
echo "[INFO] REQUIRE_EMBOSS=$REQUIRE_EMBOSS"

python3 run_pipeline.py "$@"
