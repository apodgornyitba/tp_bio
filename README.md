# TP Cuatrimestral - Parte 1 y 2 - Bioinformática

Trabajo sobre el gen **INS (Insulina)** y **diabetes mellitus** (OMIM) - Implementación Completa.

**Compatible con:** Linux, macOS y Windows.

## Requisitos

| Componente | Linux (Debian/Ubuntu) | macOS | Windows |
|------------|----------------------|-------|---------|
| Python 3.10+ | `python3`, `python3-venv` | `python3` | python.org |
| BioPython | `pip install biopython` | igual | igual |
| MSA | `sudo apt install mafft` o `muscle` | `brew install mafft` | `conda install -c bioconda mafft` |
| BLAST+ (opcional local) | `sudo apt install ncbi-blast+` | [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) o `conda install -c bioconda blast` | Instalador NCBI + PATH o `conda install -c bioconda blast` |
| EMBOSS (requerido para entrega final) | `sudo apt install emboss` | `conda install -c bioconda emboss` | `conda install -c bioconda emboss` |
| Internet | NCBI Entrez + BLAST remoto + Descarga de PROSITE | igual | igual |


## Instalación rápida

### Linux / macOS

```bash
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
python3 -m venv .venv
source .venv/bin/activate
pip install biopython

# Herramientas del sistema
# Debian/Ubuntu: sudo apt install mafft ncbi-blast+ emboss
# macOS: brew install mafft && conda install -c bioconda blast emboss

python check_requirements.py
export ENTREZ_EMAIL="tu_email@ejemplo.com"
export REQUIRE_EMBOSS=1
python run_pipeline.py
```

### Windows (PowerShell)

```powershell
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install biopython

# Herramientas del sistema
# conda install -c bioconda mafft blast emboss

python check_requirements.py
$env:ENTREZ_EMAIL="tu_email@ejemplo.com"
$env:REQUIRE_EMBOSS="1"
python run_pipeline.py
```

## Modos de BLAST (`BLAST_MODE`)

| Modo | Descripción |
|------|-------------|
| `remote` | Solo BLAST en NCBI (por defecto) |
| `local` | Solo BLAST+ local contra Swiss-Prot |
| `both` | Remoto **y** local (más puntos en la consigna) |

```bash
# Solo remoto (default)
python run_pipeline.py --blast-mode remote

# Remoto + local
python prepare_blast_db.py   # una vez (~150 MB descarga)
python run_pipeline.py --blast-mode both

# Solo local, usando los XML/local FASTA locales también para el MSA
python run_pipeline.py --blast-mode local --msa-source local
```

**Salidas BLAST:**
- Remoto: `results/blast_results/`, `results/blast_results.xml`, `results/query_best.fasta`
- Local: `results/blast_results_local/`, `results/blast_results_local.xml`, `results/query_best_local.fasta`

El Ejercicio 3 sigue la selección explícita de `--msa-source`. En modo `auto`, usa resultados locales cuando `--blast-mode=local`; en `remote` o `both`, usa los resultados remotos salvo que se indique `--msa-source local`.

## Preparar base BLAST local (opcional)

```bash
python prepare_blast_db.py
export BLAST_DB="$(pwd)/results/data/swissprot_db"   # opcional, es el default del pipeline
```

Descarga Swiss-Prot desde NCBI y ejecuta `makeblastdb`.

## Flujo del pipeline

```
fetch_data.py  →  Ex1.py  →  Ex2_a.py (remoto)  →  Ex3.py  →  Ex4.py (EMBOSS)  →  Ex5.py (Primers)
                              Ex2_local.py (opcional)
```

Orquestador multiplataforma: **`run_pipeline.py`**

Todas las salidas del pipeline se guardan dentro de `results/`.

## Qué hace cada script

| Script | Rol |
|--------|-----|
| `run_pipeline.py` | Orquesta el flujo completo (pasos 0 a 5) y valida outputs. |
| `fetch_data.py` | Descarga el GenBank de `NM_000207` desde NCBI. |
| `Ex1.py` | Traduce los 6 marcos y genera FASTA + anotación de marco CDS. |
| `Ex2_a.py` | Ejecuta BLAST remoto (NCBI) para cada marco. |
| `Ex2_local.py` | Ejecuta BLAST local contra Swiss-Prot (si se habilita modo local). |
| `Ex3.py` | Descarga top 10 hits y construye MSA con query + hits. |
| `Ex4.py` | Ejecuta ORFs y dominios PROSITE con EMBOSS (o fallback si se permite). |
| `Ex5.py` | Diseña primers según `primer_config.json`. |
| `prepare_blast_db.py` | Descarga Swiss-Prot y crea la base BLAST local. |
| `check_requirements.py` | Verifica disponibilidad de Python, MSA, BLAST y EMBOSS. |
| `platform_tools.py` | Utilidades multiplataforma para detectar ejecutables. |
| `blast_common.py` | Funciones compartidas para parseo/resumen de BLAST. |
| `run_pipeline.sh` / `run_pipeline.bat` | Wrappers opcionales para Linux/macOS y Windows. |
| `setup.sh` / `setup.bat` | Setup rápido del entorno virtual + BioPython. |

| Paso | Script | Output principal |
|------|--------|------------------|
| 0 | `fetch_data.py` | `results/NM_000207.gbk` |
| 1 | `Ex1.py` | `results/NM_000207_frames.fasta` (6 marcos) |
| 2a | `Ex2_a.py` | `results/blast_results.xml` (remoto) |
| 2 extra | `Ex2_local.py` | `results/blast_results_local.xml` (local, opcional) |
| 3 | `Ex3.py` | `results/msa_output.afa` (MSA con query CDS recortada) |
| 4 | `Ex4.py` | `results/emboss_results/NM_000207_domains.patmatmotifs` (Dominios PROSITE con EMBOSS nativo si `REQUIRE_EMBOSS=1`) |
| 5 | `Ex5.py` | `results/primer_results/primers.json` (Diseño de primers qPCR sobre variante configurada) |

## Ejecución manual (paso a paso)

```bash
# Paso 0: Descarga GenBank
python fetch_data.py results/NM_000207.gbk
# Paso 1: Traducir 6 marcos
python Ex1.py results/NM_000207.gbk results/NM_000207_frames.fasta
# Paso 2: Ejecutar BLASTp (remoto)
python Ex2_a.py results/NM_000207_frames.fasta results/blast_results
# Paso 3: Alinear hits con la consulta CDS recortada por BLAST
python Ex3.py results/blast_results.xml results/query_best.fasta --msa-input results/msa_input.fasta --msa-output results/msa_output.afa
# Paso 4: Análisis de Dominios con EMBOSS / PROSITE (falla si falta EMBOSS)
python Ex4.py results/NM_000207.gbk results/emboss_results --require-emboss
# Paso 5: Diseño de Primers qPCR parametrizado para la variante configurada
python Ex5.py results/NM_000207.gbk primer_config.json results/primer_results
```

Si existen outputs previos, estos comandos los sobrescriben con la versión actual: query CDS recortada para MSA, EMBOSS nativo obligatorio y primers enfocados en la variante `rs886037863`.

## Variables de entorno

| Variable | Default | Uso |
|----------|---------|-----|
| `ENTREZ_EMAIL` | estudiante@itba.edu.ar | Email obligatorio para NCBI |
| `BLAST_MODE` | `remote` | `remote`, `local` o `both` |
| `MSA_SOURCE` | `auto` | `auto`, `remote` o `local` para elegir qué BLAST alimenta el MSA |
| `BLAST_DB` | `./results/data/swissprot_db` | Prefijo base BLAST local |
| `SKIP_BLAST` | — | `1` para omitir BLAST (debug) |
| `REQUIRE_EMBOSS` | `0` en `python run_pipeline.py` / `1` en wrappers | `1` para exigir `getorf` y `patmatmotifs` nativos |
| `PYTHON` | `python3` / `sys.executable` | Intérprete a usar |

Las mismas opciones también pueden pasarse como argumentos, y tienen prioridad sobre las variables de entorno:

```bash
python run_pipeline.py --blast-mode local --msa-source local --require-emboss
python run_pipeline.py --blast-mode both --msa-source remote --entrez-email "tu_email@ejemplo.com"
python run_pipeline.py --allow-emboss-fallback   # solo desarrollo
```