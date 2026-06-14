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
echo "  Linux:   sudo apt install mafft ncbi-blast+ emboss"
echo "  macOS:   brew install mafft emboss"
echo "           BLAST+: descargar de NCBI o conda install -c bioconda blast"
echo "  Conda:   conda install -c bioconda mafft blast emboss"
echo ""
echo "BLAST local (opcional):"
echo "  python prepare_blast_db.py"
echo "  export BLAST_MODE=both"
echo ""
echo "Verificar herramientas instaladas:"
echo "  python check_requirements.py"
echo ""
echo "Ejecutar:"
echo "  export ENTREZ_EMAIL=tu@email.com"
echo "  export REQUIRE_EMBOSS=1"
echo "  ./run_pipeline.sh"
