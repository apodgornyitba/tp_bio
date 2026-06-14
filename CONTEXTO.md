# Contexto del TP — Bioinformática (Parte 1 + Parte 2)

> Documento de referencia rápida del proyecto completo.  
> **Repo:** https://github.com/apodgornyitba/tp_bio

---

## 1. De qué trata el trabajo

| Campo | Valor |
|-------|-------|
| **Materia** | Introducción a la Bioinformática — TP Cuatrimestral |
| **Enfermedad** | Diabetes mellitus (incl. formas monogénicas MODY10) |
| **Gen** | **INS** (insulina, *Homo sapiens*) |
| **Transcripto principal** | **NM_000207.3** (mRNA maduro, sin intrones) |
| **Proteína** | Preproinsulina → insulina (UniProt **P01308**) |
| **Herramientas** | Python, BioPython, BLAST, MAFFT/MUSCLE, EMBOSS nativo para entrega final, pipeline Bash/Python |

---

## 2. Estado del proyecto

| Componente | Estado |
|------------|--------|
| **Código Ej. 1–5** | ✅ Completo |
| **Pipeline automatizado** | ✅ `run_pipeline.py` (Linux/macOS/Windows) |
| **Outputs Parte 1** | ✅ GenBank, FASTA, BLAST, MSA |
| **Outputs Parte 2** | ✅ EMBOSS nativo + primers con variante configurada |
| **Interpretaciones BLAST/MSA** | ✅ `interpretacion_blast.md`, `interpretacion_msa.md` |
| **Ejercicio 6 (bases de datos)** | ✅ `ejercicio6_bases_datos.md` |
| **Informe base** | ✅ `INFORME_PROYECTO.md` |
| **Presentación oral (Ej. 7)** | ⏳ Pendiente del grupo (slides, 10 min) |

---

## 3. Mapa de ejercicios

```
                    OMIM / investigación
                            │
         ┌──────────────────┼──────────────────┐
         │                  │                  │
    PARTE 1              PARTE 2           NO CÓDIGO
         │                  │                  │
    fetch_data.py       Ex4.py              Ej. 6 → ejercicio6_bases_datos.md
    Ex1.py (6 marcos)   Ex5.py              Ej. 7 → presentación
    Ex2_a.py (BLAST)         │
    Ex2_local.py (opt.)      │
    Ex3.py (MSA)             │
         └──────────────────┘
                    │
              run_pipeline.py
```

| Ej. | Qué hace | Script | Output principal |
|-----|----------|--------|------------------|
| 0 | Baja mRNA INS de NCBI | `fetch_data.py` | `NM_000207.gbk` |
| 1 | 6 marcos de lectura → FASTA | `Ex1.py` | `NM_000207_frames.fasta` |
| 2a | BLAST remoto ×6 marcos | `Ex2_a.py` | `blast_results.xml`, `query_best.fasta` |
| 2b | Interpretación BLAST | — | `interpretacion_blast.md` |
| 2 local | BLAST local (opcional) | `Ex2_local.py` | `blast_results_local.xml` |
| 3 | MSA top 10 + query | `Ex3.py` | `msa_output.afa` |
| 4 | ORFs + dominios PROSITE | `Ex4.py` | `emboss_results/` |
| 5 | Diseño de primers | `Ex5.py` | `primer_results/` |
| 6 | Bases de datos | — | `ejercicio6_bases_datos.md` |
| 7 | Presentación | — | slides del grupo |

---

## 4. Resultados clave (para oral / informe)

### Ejercicio 1
- Solo **`Forward_Frame_3`** coincide con CDS y produce preproinsulina.
- Los otros 5 marcos no tienen sentido biológico.

### Ejercicio 2 (BLAST remoto)
- Marco ganador: **Forward_Frame_3**.
- Mejor hit: insulina humana **P01308**, E-value **~10⁻⁷⁶**, **100% identidad**.
- Top hits: primates (humano, orangután, chimpancé, macacos).

### Ejercicio 3 (MSA)
- 11 secuencias (query + top 10 BLAST).
- Núcleo de insulina **muy conservado**; más divergencia en roedores.

### Ejercicio 4 (EMBOSS / PROSITE)
- 14 ORFs detectados por `getorf` nativo; firma **INSULIN (PS00262)** en `NM_000207.3_4`.
- Para entrega, correr con EMBOSS nativo (`REQUIRE_EMBOSS=1`) para evitar evidencia basada en fallback.

### Ejercicio 5 (Primers)
- 5 primers sobre **NM_000207.3 con la variante c.125T>C aplicada** (18–24 bp, GC 50–60%, Tm ≤ 67°C).
- Config: `primer_config.json` → reporte regenerado en `primer_results/primers_report.txt`.

### Ejercicio 6
- Gen **INS** en NCBI Gene: https://www.ncbi.nlm.nih.gov/gene/3630
- Variante ejemplo: **rs886037863** (V42A) → **MODY10**.
- Detalle completo: `ejercicio6_bases_datos.md`.

---

## 5. Cómo ejecutar

### Setup (una vez)

```bash
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
./setup.sh && source .venv/bin/activate   # Linux/macOS
# setup.bat en Windows
export ENTREZ_EMAIL="tu_email@ejemplo.com"
export REQUIRE_EMBOSS=1
brew install mafft emboss                 # macOS (MSA + EMBOSS)
python check_requirements.py
```

### Pipeline completo

```bash
./run_pipeline.sh                         # Linux/macOS
# run_pipeline.bat en Windows
```

### Variables útiles

| Variable | Default | Uso |
|----------|---------|-----|
| `ENTREZ_EMAIL` | estudiante@itba.edu.ar | Obligatorio NCBI |
| `BLAST_MODE` | `remote` | `remote`, `local`, `both` |
| `BLAST_DB` | `./data/swissprot_db` | Base BLAST local |
| `REQUIRE_EMBOSS` | — | `1` para exigir getorf + patmatmotifs nativos |

### Solo Parte 2 (si Parte 1 ya corrida)

```bash
python Ex4.py NM_000207.gbk emboss_results --require-emboss
python Ex5.py NM_000207.gbk primer_config.json primer_results
```

---

## 6. Archivos del repositorio

### Scripts
```
fetch_data.py    Ex1.py    Ex2_a.py    Ex2_local.py    Ex3.py
Ex4.py           Ex5.py    prepare_blast_db.py
blast_common.py  platform_tools.py
check_requirements.py
run_pipeline.py  run_pipeline.sh  run_pipeline.bat
setup.sh         setup.bat
```

### Datos generados
```
NM_000207.gbk                 NM_000207_frames.fasta
frame_annotation.txt          query_best.fasta
blast_results/                blast_results.xml
msa_input.fasta               msa_output.afa
emboss_results/               prosite.dat
primer_results/               primer_config.json
```

### Documentación
```
README.md                     consigna.md          parte2.md
CONTEXTO.md                   INFORME_PROYECTO.md
interpretacion_blast.md       interpretacion_msa.md
ejercicio6_bases_datos.md
```

---

## 7. Documentos — cuándo usar cada uno

| Archivo | Para qué sirve |
|---------|----------------|
| **CONTEXTO.md** | Vista rápida de todo el TP (este archivo) |
| **README.md** | Instalación y ejecución técnica |
| **INFORME_PROYECTO.md** | Base para redactar informe formal (Claude, Word, etc.) |
| **interpretacion_blast.md** | Entrega Ej. 2.b |
| **interpretacion_msa.md** | Entrega Ej. 3 |
| **ejercicio6_bases_datos.md** | Entrega Ej. 6 |
| **consigna.md** | Parte 1 oficial |
| **parte2.md** | Parte 2 oficial (Ej. 4–7) |

---

## 8. Qué falta hacer (solo el grupo)

1. **Slides presentación (Ej. 7)** — 10 min, biología, no código.
2. **Ensayo oral** — repartir entre integrantes.
3. **Antes de entregar:** correr Ex4 con EMBOSS real (`sudo apt install emboss`, `brew install emboss` o `conda install -c bioconda emboss`).
4. **Opcional:** informe PDF final a partir de `INFORME_PROYECTO.md` + `ejercicio6_bases_datos.md`.

### Estructura sugerida presentación (10 min)

1. Diabetes + OMIM + gen INS  
2. Parte 1: marcos → BLAST → MSA  
3. Parte 2: PROSITE + primers  
4. Ej. 6: NCBI, Ensembl, variante rs886037863  
5. Conclusión  

---

## 9. Links esenciales

| Recurso | URL |
|---------|-----|
| Repo GitHub | https://github.com/apodgornyitba/tp_bio |
| NCBI Gene INS | https://www.ncbi.nlm.nih.gov/gene/3630 |
| OMIM INS | https://omim.org/entry/176730 |
| UniProt P01308 | https://www.uniprot.org/uniprotkb/P01308 |
| Ensembl INS | https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000254647 |
| ClinVar INS | https://www.ncbi.nlm.nih.gov/clinvar/?term=INS[gene] |
| Reactome | https://reactome.org/content/detail/R-HSA-264876 |

---

## 10. Prompt útil para generar informe o presentación

```
Usá CONTEXTO.md, INFORME_PROYECTO.md y ejercicio6_bases_datos.md de este repo.
Redactá [informe académico / guion de presentación 10 min] sobre el TP de
Bioinformática del gen INS y diabetes. No inventes datos. Tono universitario en español.
```

---

*Última actualización: Mayo 2026 — tp_bio ITBA*
