# TP Cuatrimestral - Parte 1 y 2 - Bioinformática

Trabajo sobre el gen **INS (Insulina)** y **diabetes mellitus** (OMIM) - Implementación Completa.

**Compatible con:** Linux, macOS y Windows.

## Requisitos

| Componente | Linux (Debian/Ubuntu) | macOS | Windows |
|------------|----------------------|-------|---------|
| Python 3.10+ | `python3`, `python3-venv` | `python3` | python.org |
| BioPython | `pip install biopython` | igual | igual |
| MSA | `sudo apt install mafft` o `muscle` | `brew install mafft` | `conda install -c bioconda mafft` |
| BLAST+ (opcional local) | `sudo apt install ncbi-blast+` | [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) | Instalador NCBI + PATH |
| EMBOSS (opcional) | `sudo apt install emboss` | `brew install emboss` | [EMBOSS](https://emboss.sourceforge.net/) |
| Internet | NCBI Entrez + BLAST remoto + Descarga de PROSITE | igual | igual |

> [!NOTE]
> **Robustez del Pipeline:** Si la suite **EMBOSS** no está instalada localmente o falla, el pipeline cuenta con un **fallback en Python puro** de alta precisión para los programas `getorf` y `patmatmotifs` (incluye parser regex de `prosite.dat`). ¡El pipeline se ejecutará correctamente sin dependencias externas complejas!

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
fetch_data.py  →  Ex1.py  →  Ex2_a.py (remoto)  →  Ex3.py  →  Ex4.py (EMBOSS)  →  Ex5.py (Primers)
                              Ex2_local.py (opcional)
```

Orquestador multiplataforma: **`run_pipeline.py`**

| Paso | Script | Output principal |
|------|--------|------------------|
| 0 | `fetch_data.py` | `NM_000207.gbk` |
| 1 | `Ex1.py` | `NM_000207_frames.fasta` (6 marcos) |
| 2a | `Ex2_a.py` | `blast_results.xml` (remoto) |
| 2b | `Ex2_local.py` | `blast_results_local.xml` (local, opcional) |
| 3 | `Ex3.py` | `msa_output.afa` (Alineamiento múltiple) |
| 4 | `Ex4.py` | `emboss_results/NM_000207_domains.patmatmotifs` (Dominios PROSITE) |
| 5 | `Ex5.py` | `primer_results/primers.json` (Diseño de primers qPCR) |

## Ejecución manual (paso a paso)

```bash
# Paso 0: Descarga GenBank
python fetch_data.py
# Paso 1: Traducir 6 marcos
python Ex1.py NM_000207.gbk NM_000207_frames.fasta
# Paso 2: Ejecutar BLASTp (remoto)
python Ex2_a.py NM_000207_frames.fasta blast_results
# Paso 3: Alinear hits con la consulta (MSA)
python Ex3.py blast_results.xml query_best.fasta
# Paso 4: Análisis de Dominios con EMBOSS / PROSITE
python Ex4.py NM_000207.gbk emboss_results
# Paso 5: Diseño de Primers qPCR parametrizado
python Ex5.py NM_000207.gbk primer_config.json primer_results
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
- Archivo de configuración: `primer_config.json`
- Outputs completos del pipeline (incluyendo `emboss_results/` y `primer_results/`)
- Documento local descargado: `prosite.dat`
- `interpretacion_blast.md`, `interpretacion_msa.md`
- `INFORME_PROYECTO.md` (informe completo de reporte integrado, base académica para la entrega)
- Presentación oral: investigación del gen INS y diabetes (**sin código**, 10 min de exposición)

## Notas

- BLAST remoto: ~30–60 min (6 consultas).
- BLAST local: más rápido tras tener la base descargada.
- En macOS no hay `muscle` en Homebrew; se usa **MAFFT** automáticamente.
