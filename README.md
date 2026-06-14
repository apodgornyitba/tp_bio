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

> [!NOTE]
> **Entrega final:** usar `REQUIRE_EMBOSS=1` o `python Ex4.py ... --require-emboss` para demostrar que el Ejercicio 4 corre con **EMBOSS nativo** (`getorf` + `patmatmotifs`). El fallback Python queda solo como modo de desarrollo/portabilidad.
> `patmatmotifs` requiere que PROSITE se prepare con `prosextract`; `Ex4.py` descarga `prosite.dat` y `prosite.doc`, genera el índice local en `.emboss_data/` y evita modificar la instalación de EMBOSS.

## Instalación rápida

### Linux / macOS

```bash
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
chmod +x setup.sh run_pipeline.sh
./setup.sh
source .venv/bin/activate
# Instalar antes herramientas del sistema:
# Linux: sudo apt install mafft ncbi-blast+ emboss
# macOS: brew install mafft && conda install -c bioconda emboss
python check_requirements.py
export ENTREZ_EMAIL="tu_email@ejemplo.com"
export REQUIRE_EMBOSS=1
./run_pipeline.sh
```

### Windows (CMD o PowerShell)

```cmd
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
setup.bat
.venv\Scripts\activate.bat
conda install -c bioconda mafft blast emboss
python check_requirements.py
set ENTREZ_EMAIL=tu_email@ejemplo.com
set REQUIRE_EMBOSS=1
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
./run_pipeline.sh --blast-mode remote

# Remoto + local
python prepare_blast_db.py   # una vez (~150 MB descarga)
./run_pipeline.sh --blast-mode both

# Solo local, usando los XML/local FASTA locales también para el MSA
./run_pipeline.sh --blast-mode local --msa-source local
```

**Salidas BLAST:**
- Remoto: `blast_results/`, `blast_results.xml`, `query_best.fasta`
- Local: `blast_results_local/`, `blast_results_local.xml`, `query_best_local.fasta`

El Ejercicio 3 sigue la selección explícita de `--msa-source`. En modo `auto`, usa resultados locales cuando `--blast-mode=local`; en `remote` o `both`, usa los resultados remotos salvo que se indique `--msa-source local`.

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
| 2 extra | `Ex2_local.py` | `blast_results_local.xml` (local, opcional) |
| 2b | `interpretacion_blast.md` | Interpretación biológica y estadística del BLAST |
| 3 | `Ex3.py` | `msa_output.afa` (MSA con query CDS recortada) |
| 4 | `Ex4.py` | `emboss_results/NM_000207_domains.patmatmotifs` (Dominios PROSITE con EMBOSS nativo si `REQUIRE_EMBOSS=1`) |
| 5 | `Ex5.py` | `primer_results/primers.json` (Diseño de primers qPCR sobre variante configurada) |

## Ejecución manual (paso a paso)

```bash
# Paso 0: Descarga GenBank
python fetch_data.py
# Paso 1: Traducir 6 marcos
python Ex1.py NM_000207.gbk NM_000207_frames.fasta
# Paso 2: Ejecutar BLASTp (remoto)
python Ex2_a.py NM_000207_frames.fasta blast_results
# Paso 3: Alinear hits con la consulta CDS recortada por BLAST
python Ex3.py blast_results.xml query_best.fasta
# Paso 4: Análisis de Dominios con EMBOSS / PROSITE (falla si falta EMBOSS)
python Ex4.py NM_000207.gbk emboss_results --require-emboss
# Paso 5: Diseño de Primers qPCR parametrizado para la variante configurada
python Ex5.py NM_000207.gbk primer_config.json primer_results
```

Si existen outputs previos, estos comandos los sobrescriben con la versión actual: query CDS recortada para MSA, EMBOSS nativo obligatorio y primers enfocados en la variante `rs886037863`.

## Variables de entorno

| Variable | Default | Uso |
|----------|---------|-----|
| `ENTREZ_EMAIL` | estudiante@itba.edu.ar | Email obligatorio para NCBI |
| `BLAST_MODE` | `remote` | `remote`, `local` o `both` |
| `MSA_SOURCE` | `auto` | `auto`, `remote` o `local` para elegir qué BLAST alimenta el MSA |
| `BLAST_DB` | `./data/swissprot_db` | Prefijo base BLAST local |
| `SKIP_BLAST` | — | `1` para omitir BLAST (debug) |
| `REQUIRE_EMBOSS` | `1` en wrappers | `1` para exigir `getorf` y `patmatmotifs` nativos |
| `PYTHON` | `python3` / `sys.executable` | Intérprete a usar |

Las mismas opciones también pueden pasarse como argumentos, y tienen prioridad sobre las variables de entorno:

```bash
python run_pipeline.py --blast-mode local --msa-source local --require-emboss
python run_pipeline.py --blast-mode both --msa-source remote --entrez-email "tu_email@ejemplo.com"
python run_pipeline.py --allow-emboss-fallback   # solo desarrollo
```

## Entregables

- Scripts Python + `run_pipeline.py` / `.sh` / `.bat`
- Archivo de configuración: `primer_config.json`
- Verificación local de dependencias: `check_requirements.py`
- Outputs completos del pipeline (incluyendo `emboss_results/` y `primer_results/`)
- Documentos PROSITE descargados: `prosite.dat`, `prosite.doc`
- `interpretacion_blast.md`, `interpretacion_msa.md`
- `INFORME_PROYECTO.md` (informe completo de reporte integrado, base académica para la entrega)
- Presentación oral: investigación del gen INS y diabetes (**sin código**, 10 min de exposición)

## Notas

- BLAST remoto: ~30–60 min (6 consultas).
- BLAST local: más rápido tras tener la base descargada.
- En macOS no hay `muscle` en Homebrew; se usa **MAFFT** automáticamente.
