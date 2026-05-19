# TP Cuatrimestral - Parte 1 - Bioinformática

Trabajo sobre el gen **INS (Insulina)** y **diabetes mellitus** (OMIM).

**Compatible con:** Linux, macOS y Windows.

## Requisitos

| Componente | Linux (Debian/Ubuntu) | macOS | Windows |
|------------|----------------------|-------|---------|
| Python 3.10+ | `python3`, `python3-venv` | `python3` | python.org |
| BioPython | `pip install biopython` | igual | igual |
| MSA | `sudo apt install mafft` o `muscle` | `brew install mafft` | `conda install -c bioconda mafft` |
| BLAST+ (opcional local) | `sudo apt install ncbi-blast+` | [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) | Instalador NCBI + PATH |
| Internet | NCBI Entrez + BLAST remoto | igual | igual |

## Instalación rápida

### Linux / macOS

```bash
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
chmod +x setup.sh run_pipeline.sh
./setup.sh
source .venv/bin/activate
export ENTREZ_EMAIL="tu_email@ejemplo.com"
./run_pipeline.sh
```

### Windows (CMD o PowerShell)

```cmd
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
setup.bat
.venv\Scripts\activate.bat
set ENTREZ_EMAIL=tu_email@ejemplo.com
run_pipeline.bat
```

## Modos de BLAST (`BLAST_MODE`)

| Modo | Descripción |
|------|-------------|
| `remote` | Solo BLAST en NCBI (por defecto) |
| `local` | Solo BLAST+ local contra Swiss-Prot |
| `both` | Remoto **y** local (más puntos en la consigna) |

```bash
# Solo remoto (default)
export BLAST_MODE=remote
./run_pipeline.sh

# Remoto + local
export BLAST_MODE=both
python prepare_blast_db.py   # una vez (~150 MB descarga)
./run_pipeline.sh
```

**Salidas BLAST:**
- Remoto: `blast_results/`, `blast_results.xml`, `query_best.fasta`
- Local: `blast_results_local/`, `blast_results_local.xml`, `query_best_local.fasta`

El Ejercicio 3 usa los resultados **remotos** si existen; si solo corrés local, usa los locales.

## Preparar base BLAST local (opcional)

```bash
python prepare_blast_db.py
export BLAST_DB="$(pwd)/data/swissprot_db"   # opcional, es el default
```

Descarga Swiss-Prot desde NCBI y ejecuta `makeblastdb`.

## Flujo del pipeline

```
fetch_data.py  →  Ex1.py  →  Ex2_a.py (remoto)  →  Ex3.py
                              Ex2_local.py (opcional)
```

Orquestador multiplataforma: **`run_pipeline.py`**

| Paso | Script | Output principal |
|------|--------|------------------|
| 0 | `fetch_data.py` | `NM_000207.gbk` |
| 1 | `Ex1.py` | `NM_000207_frames.fasta` (6 marcos) |
| 2a | `Ex2_a.py` | `blast_results.xml` |
| 2b | `Ex2_local.py` | `blast_results_local.xml` |
| 3 | `Ex3.py` | `msa_output.afa` |

## Ejecución manual (paso a paso)

```bash
python fetch_data.py
python Ex1.py NM_000207.gbk NM_000207_frames.fasta
python Ex2_a.py NM_000207_frames.fasta blast_results
python Ex2_local.py NM_000207_frames.fasta blast_results_local   # opcional
python Ex3.py blast_results.xml query_best.fasta
```

## Variables de entorno

| Variable | Default | Uso |
|----------|---------|-----|
| `ENTREZ_EMAIL` | estudiante@itba.edu.ar | Email obligatorio para NCBI |
| `BLAST_MODE` | `remote` | `remote`, `local` o `both` |
| `BLAST_DB` | `./data/swissprot_db` | Prefijo base BLAST local |
| `SKIP_BLAST` | — | `1` para omitir BLAST (debug) |
| `PYTHON` | `python3` / `sys.executable` | Intérprete a usar |

## Entregables

- Scripts Python + `run_pipeline.py` / `.sh` / `.bat`
- Outputs del pipeline
- `interpretacion_blast.md`, `interpretacion_msa.md`
- `INFORME_PROYECTO.md` (base para informe escrito)
- Presentación oral: investigación OMIM/INS (**sin código**)

## Notas

- BLAST remoto: ~30–60 min (6 consultas).
- BLAST local: más rápido tras tener la base descargada.
- En macOS no hay `muscle` en Homebrew; se usa **MAFFT** automáticamente.
