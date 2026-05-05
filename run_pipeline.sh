#!/bin/bash

# run_pipeline.sh
# Automatización del flujo de tareas para el TP de Bioinformática
# Realiza logueo, control de errores y verificación de formatos.

LOG_FILE="pipeline.log"
exec > >(tee -i $LOG_FILE)
exec 2>&1

echo "=================================================="
echo " INICIO DEL PIPELINE DE BIOINFORMÁTICA"
echo " Fecha: $(date)"
echo "=================================================="

# 0. Descarga de datos
echo "[INFO] Paso 0: Obteniendo datos (NM_000207 - INS)"
python3 fetch_data.py
if [ $? -ne 0 ] || [ ! -f "NM_000207.gbk" ]; then
    echo "[ERROR] Falló la descarga del archivo GenBank."
    exit 1
fi
echo "[OK] Archivo GenBank descargado correctamente."

# 1. Ejercicio 1 - Procesamiento de Secuencias
echo -e "\n[INFO] Paso 1: Ejecutando Ejercicio 1 (Traducción de ORF)"
python3 Ex1.py NM_000207.gbk NM_000207.fasta
if [ $? -ne 0 ] || [ ! -f "NM_000207.fasta" ]; then
    echo "[ERROR] Falló la ejecución del Ejercicio 1."
    exit 1
fi
echo "[OK] Secuencia traducida a formato FASTA."

# 2. Ejercicio 2 - BLAST
echo -e "\n[INFO] Paso 2: Ejecutando Ejercicio 2 (BLAST remoto)"
python3 Ex2_a.py NM_000207.fasta blast_results.xml
if [ $? -ne 0 ] || [ ! -f "blast_results.xml" ]; then
    echo "[ERROR] Falló la ejecución del Ejercicio 2 (BLAST)."
    exit 1
fi
echo "[OK] Resultados de BLAST obtenidos."

# 3. Ejercicio 3 - Multiple Sequence Alignment (MSA)
echo -e "\n[INFO] Paso 3: Ejecutando Ejercicio 3 (MSA con MUSCLE)"
python3 Ex3.py blast_results.xml NM_000207.fasta
if [ $? -ne 0 ] || [ ! -f "msa_output.afa" ]; then
    echo "[ERROR] Falló la ejecución del Ejercicio 3 (MSA)."
    exit 1
fi
echo "[OK] Alineamiento múltiple completado."

echo -e "\n=================================================="
echo " PIPELINE FINALIZADO CORRECTAMENTE"
echo " Los resultados están listos para la interpretación."
echo " Archivos generados: NM_000207.gbk, NM_000207.fasta, blast_results.xml, msa_input.fasta, msa_output.afa"
echo "=================================================="
