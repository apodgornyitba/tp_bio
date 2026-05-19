#!/bin/bash
# Crea entorno virtual e instala dependencias
set -e
python3 -m venv .venv
source .venv/bin/activate
pip install biopython
echo "Listo. Activar con: source .venv/bin/activate"
echo "MSA en macOS: brew install mafft"
echo "Luego: export ENTREZ_EMAIL=tu@email.com && ./run_pipeline.sh"
