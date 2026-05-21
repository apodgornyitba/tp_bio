# TP Cuatrimestral - Parte 1 - Bioinformatica

Este repositorio contiene los scripts requeridos para el TP Cuatrimestral de Bioinformatica.
El trabajo esta centrado en el gen humano INS (Insulina) y su contexto en diabetes mellitus.

## Requisitos previos

El entorno requiere tener instalado Python 3.10+, Biopython y una herramienta de alineamiento multiple (MUSCLE o MAFFT).
Para BLAST local (opcional), tambien se necesita BLAST+.

Ejemplo en Ubuntu/WSL:

```bash
sudo apt-get update
sudo apt-get install -y python3 python3-venv ncbi-blast+ muscle mafft
python3 -m pip install biopython
```

## Archivos y scripts

1. `fetch_data.py`: descarga la secuencia GenBank de NM_000207 desde NCBI.
2. `Ex1.py`: genera traducciones de los 6 marcos de lectura y anota el marco que coincide con la CDS.
3. `Ex2_a.py`: ejecuta BLAST remoto contra Swiss-Prot para cada marco y selecciona la mejor secuencia.
4. `Ex2_local.py`: ejecuta BLAST local (opcional) usando una base Swiss-Prot descargada.
5. `Ex3.py`: toma el mejor resultado BLAST, recupera homologos y corre MSA con MUSCLE o MAFFT.
6. `prepare_blast_db.py`: descarga y prepara la base local para BLAST+.
7. `run_pipeline.py`: orquestador principal multiplataforma.
8. `run_pipeline.sh` y `run_pipeline.bat`: wrappers para Linux/macOS y Windows.

## Modos de BLAST

La variable `BLAST_MODE` controla el modo de ejecucion:

- `remote`: solo BLAST remoto (default).
- `local`: solo BLAST local.
- `both`: ejecuta remoto y local.

## Como ejecutar

### Linux / macOS

```bash
chmod +x run_pipeline.sh
export ENTREZ_EMAIL="tu_email@ejemplo.com"
export BLAST_MODE=remote
./run_pipeline.sh
```

### Windows (CMD o PowerShell)

```cmd
set ENTREZ_EMAIL=tu_email@ejemplo.com
set BLAST_MODE=remote
run_pipeline.bat
```

### BLAST local (opcional)

Setup minimo para modo local:

1. Instalar BLAST+ (`blastp` y `makeblastdb`)
2. Descargar y preparar la base Swiss-Prot
3. Ejecutar el pipeline con `BLAST_MODE=local` o `BLAST_MODE=both`

Linux / macOS:

```bash
sudo apt-get install -y ncbi-blast+
python prepare_blast_db.py
export BLAST_DB="$(pwd)/data/swissprot_db"
export BLAST_MODE=local
./run_pipeline.sh
```

Windows (PowerShell):

```powershell
# Requiere BLAST+ instalado y en PATH
python prepare_blast_db.py
setx BLAST_DB "$PWD\data\swissprot_db"
$env:BLAST_MODE="local"
run_pipeline.bat
```

Para ejecutar remoto y local en la misma corrida:

```bash
export BLAST_MODE=both
./run_pipeline.sh
```

## Salidas del pipeline

Todas las salidas se guardan en la carpeta `results/`:

- `results/NM_000207.gbk`
- `results/NM_000207_frames.fasta`
- `results/frame_annotation.txt`
- `results/blast_results/` y `results/blast_results.xml`
- `results/query_best.fasta`
- `results/blast_results_local/` y `results/blast_results_local.xml` (si corresponde)
- `results/query_best_local.fasta` (si corresponde)
- `results/msa_input.fasta`
- `results/msa_output.afa`
- `results/pipeline.log`

## Notas

- Si hay resultados remotos, Ejercicio 3 prioriza esos XML; si solo hay local, usa el local.
- BLAST remoto puede tardar varios minutos.
- En macOS normalmente se usa MAFFT cuando MUSCLE no esta disponible.
