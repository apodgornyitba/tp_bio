#!/bin/bash

# run_pipeline.sh
# Automatización del flujo de tareas para el TP de Bioinformática
# Realiza logueo, control de errores y verificación de formatos.

LOG_FILE="pipeline.log"
exec > >(tee -i $LOG_FILE)
exec 2>&1

# Usar entorno virtual si existe (recomendado en macOS)
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
fi

echo "=================================================="
echo " INICIO DEL PIPELINE DE BIOINFORMÁTICA"
echo " Fecha: $(date)"
echo "=================================================="

# Email requerido por NCBI (cambiar por el del grupo)
export ENTREZ_EMAIL="${ENTREZ_EMAIL:-estudiante@itba.edu.ar}"
echo "[INFO] Usando ENTREZ_EMAIL=$ENTREZ_EMAIL"

# 0. Descarga de datos
echo "[INFO] Paso 0: Obteniendo datos (NM_000207 - INS)"
python3 fetch_data.py
if [ $? -ne 0 ] || [ ! -f "NM_000207.gbk" ]; then
    echo "[ERROR] Falló la descarga del archivo GenBank."
    exit 1
fi
echo "[OK] Archivo GenBank descargado correctamente."

# 1. Ejercicio 1 - Procesamiento de Secuencias (6 marcos de lectura)
echo -e "\n[INFO] Paso 1: Ejecutando Ejercicio 1 (6 marcos de lectura)"
python3 Ex1.py NM_000207.gbk NM_000207_frames.fasta
if [ $? -ne 0 ] || [ ! -f "NM_000207_frames.fasta" ]; then
    echo "[ERROR] Falló la ejecución del Ejercicio 1."
    exit 1
fi
echo "[OK] 6 traducciones guardadas en NM_000207_frames.fasta"

# 2. Ejercicio 2 - BLAST de cada marco (puede tardar 30-60 min)
echo -e "\n[INFO] Paso 2: Ejecutando Ejercicio 2 (BLAST remoto por marco)"
echo "[INFO] Se ejecutará un BLAST por cada una de las 6 secuencias..."
python3 Ex2_a.py NM_000207_frames.fasta blast_results
if [ $? -ne 0 ] || [ ! -f "blast_results.xml" ] || [ ! -f "query_best.fasta" ]; then
    echo "[ERROR] Falló la ejecución del Ejercicio 2 (BLAST)."
    exit 1
fi
echo "[OK] Resultados de BLAST obtenidos. Mejor marco en query_best.fasta"

# 3. Ejercicio 3 - Multiple Sequence Alignment (MSA)
echo -e "\n[INFO] Paso 3: Ejecutando Ejercicio 3 (MSA con MUSCLE)"
python3 Ex3.py blast_results.xml query_best.fasta
if [ $? -ne 0 ] || [ ! -f "msa_output.afa" ]; then
    echo "[ERROR] Falló la ejecución del Ejercicio 3 (MSA)."
    exit 1
fi
echo "[OK] Alineamiento múltiple completado."

echo -e "\n=================================================="
echo " PIPELINE FINALIZADO CORRECTAMENTE"
echo " Archivos generados:"
echo "  - NM_000207.gbk, NM_000207_frames.fasta, frame_annotation.txt"
echo "  - blast_results/ (XML por marco), blast_results.xml, query_best.fasta"
echo "  - msa_input.fasta, msa_output.afa"
echo " Completar interpretacion_blast.md e interpretacion_msa.md con los resultados."
echo "=================================================="
