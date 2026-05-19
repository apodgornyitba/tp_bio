@echo off
REM Pipeline TP Bioinformatica - Windows
cd /d "%~dp0"

if exist .venv\Scripts\activate.bat (
    call .venv\Scripts\activate.bat
)

if "%ENTREZ_EMAIL%"=="" set ENTREZ_EMAIL=estudiante@itba.edu.ar
if "%BLAST_MODE%"=="" set BLAST_MODE=remote

echo Usando BLAST_MODE=%BLAST_MODE%
python run_pipeline.py
exit /b %ERRORLEVEL%
