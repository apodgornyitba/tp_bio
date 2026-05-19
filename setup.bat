@echo off
REM Setup Windows - entorno virtual + BioPython
cd /d "%~dp0"

python -m venv .venv
call .venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install biopython

echo.
echo Listo. Siguiente:
echo   .venv\Scripts\activate.bat
echo   set ENTREZ_EMAIL=tu@email.com
echo.
echo BLAST+ : instalar desde NCBI y agregar al PATH
echo MSA    : conda install -c bioconda mafft
echo.
echo   set BLAST_MODE=both
echo   run_pipeline.bat
