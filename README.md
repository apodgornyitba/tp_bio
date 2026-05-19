# TP Cuatrimestral - Parte 1 - Bioinformática

Trabajo práctico sobre el gen humano **INS (Insulina)**, asociado a la **diabetes mellitus** (consulta en OMIM).

## Requisitos previos

- Linux (o Workspace ITBA)
- `python3` + `biopython`
- `muscle` o `mafft` (alineamiento múltiple; en macOS: `brew install mafft`)
- Conexión a internet (NCBI Entrez + BLAST remoto)

```bash
sudo apt-get update && sudo apt-get install -y python3-biopython muscle
# o: pip install biopython
```

Configurar email para NCBI (obligatorio):

```bash
export ENTREZ_EMAIL="tu_email@ejemplo.com"
```

## Flujo del pipeline

```
OMIM/INS → NM_000207.gbk → 6 marcos (FASTA) → BLAST×6 → mejor marco → MSA
```

| Paso | Script | Input | Output |
|------|--------|-------|--------|
| 0 | `fetch_data.py` | — | `NM_000207.gbk` |
| 1 | `Ex1.py` | `.gbk` | `NM_000207_frames.fasta` (6 secuencias) |
| 2 | `Ex2_a.py` | frames FASTA | `blast_results/`, `blast_results.xml`, `query_best.fasta` |
| 3 | `Ex3.py` | XML + query | `msa_input.fasta`, `msa_output.afa` |

## Ejecución

```bash
chmod +x run_pipeline.sh
export ENTREZ_EMAIL="tu_email@ejemplo.com"
./run_pipeline.sh
```

**Nota:** el Ejercicio 2 ejecuta **6 BLAST remotos** (uno por marco). Puede tardar **30–60 minutos**.

Ejecución paso a paso:

```bash
python3 fetch_data.py
python3 Ex1.py NM_000207.gbk NM_000207_frames.fasta
python3 Ex2_a.py NM_000207_frames.fasta blast_results
python3 Ex3.py blast_results.xml query_best.fasta
```

## Entregables

- Scripts: `Ex1.py`, `Ex2_a.py`, `Ex3.py`, `run_pipeline.sh`
- Inputs/outputs generados por el pipeline
- `interpretacion_blast.md` (Ej. 2.b) – completar tras correr BLAST
- `interpretacion_msa.md` (Ej. 3) – completar tras correr MSA
- Presentación oral: investigación OMIM + INS + diabetes (**sin código**)

## Investigación (para la exposición)

1. Buscar **diabetes mellitus** y gen **INS** en [OMIM](https://omim.org/).
2. Explicar función de la insulina (células β del páncreas, glucemia).
3. Por qué se usa **NM_000207** (mRNA maduro, sin intrones).
4. Interpretar resultados BLAST y MSA con los archivos generados.
