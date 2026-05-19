#!/bin/bash
# Setup Linux / macOS
set -e
cd "$(dirname "$0")"

python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install biopython

echo ""
echo "Listo. Activar: source .venv/bin/activate"
echo ""
echo "Dependencias del sistema:"
echo "  Linux:   sudo apt install mafft ncbi-blast+"
echo "  macOS:   brew install mafft"
echo "           BLAST+: descargar de NCBI o conda install -c bioconda blast"
echo ""
echo "BLAST local (opcional):"
echo "  python prepare_blast_db.py"
echo "  export BLAST_MODE=both"
echo ""
echo "Ejecutar:"
echo "  export ENTREZ_EMAIL=tu@email.com"
echo "  ./run_pipeline.sh"
