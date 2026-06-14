# Documento base para informe — TP Bioinformática Parte 1 y 2

> **Uso:** Entregar este archivo a Claude (u otro asistente) con la instrucción:  
> *"Armá un informe académico formal en español a partir de este documento, con introducción, materiales y métodos, resultados, discusión y conclusiones."*  
>  
> **Repositorio:** https://github.com/apodgornyitba/tp_bio  
> **Materia:** Introducción a la Bioinformática — Trabajo Práctico Cuatrimestral Partes 1 y 2

---

## 1. Datos del trabajo

| Campo | Valor |
|-------|-------|
| **Título del TP** | Trabajo Práctico Cuatrimestral — Parte 2 — Introducción a la Bioinformática |
| **Enfermedad investigada** | Diabetes mellitus |
| **Gen elegido** | **INS** (Insulina, *Homo sapiens*) |
| **Transcripto mRNA** | **NM_000207.3** (RefSeq, mRNA maduro, sin intrones) |
| **Organismo** | *Homo sapiens* |
| **Herramientas** | Python 3, BioPython, BLASTp, MAFFT, EMBOSS (getorf, patmatmotifs), PROSITE, Bash |
| **Base de datos BLAST** | Swiss-Prot (curada) |
| **Lenguaje** | BioPython (extensión `.py`; la consigna sugiere BioPerl `.pm` pero permite otros) |
| **Ejercicio 4 y 5** | Implementados para corrida final con EMBOSS nativo obligatorio y diseño de primers sobre la variante configurada |

---

## 2. Qué pide la consigna (requisitos oficiales)

### 2.1 Requisitos generales

- Elegir una **enfermedad** y genes asociados consultando **OMIM** (Online Mendelian Inheritance in Man).
- Desarrollar scripts en **BioPython** (u otro lenguaje bioinformático open source).
- Automatizar la ejecución con **Bash**: control de errores, validación de archivos, logging.
- **Exposición oral de 10 minutos** sobre la **investigación** (no sobre el código).
- Entregar scripts, archivos input/output, descripción breve y cómo ejecutar.

### 2.2 Ejercicio 1 — Procesamiento de secuencias

**Pedido:**
- Leer secuencias de **nucleótidos** en formato **GenBank** (mRNA maduro del gen, sin intrones).
- Traducir considerando los **6 marcos de lectura** (3 forward + 3 reverse complement).
- Escribir resultados en formato **FASTA** (una o más secuencias de aminoácidos por ORF/marco).
- El archivo GenBank debe obtenerse de NCBI (ejemplo de la consigna: gen **INS** asociado a diabetes).

**Aclaración de la consigna:**
- Deben calcular los **6 marcos**, evaluarlos en el Ejercicio 2 con BLAST, y así identificar el marco real.
- Alternativamente pueden usar anotación CDS/OrfFinder para identificar el marco correcto.

**Input entregado:** `NM_000207.gbk`  
**Output entregado:** `NM_000207_frames.fasta` (6 secuencias)

### 2.3 Ejercicio 2.a — BLAST

**Pedido:**
- Realizar **BLAST** de una o varias secuencias (un BLAST por cada secuencia del FASTA).
- Guardar el reporte BLAST en archivo(s).
- Puede ser **remoto** o **local** (ambos dan más puntos; este proyecto soporta ambos modos y conserva resultados remotos/locales cuando se ejecutan).
- Entregar interpretación de resultados y cómo ejecutar.

**Input:** `NM_000207_frames.fasta`  
**Output:** `blast_results/` (XML por marco), `blast_results.xml`, `query_best.fasta`, `blast_results/blast_summary.txt`

### 2.4 Ejercicio 2.b — Interpretación del BLAST

**Pedido:**
- Explicar las secuencias encontradas en términos biológicos.
- Explicar valores estadísticos: **E-value**, score, identidad, gaps, etc.
- Referencia sugerida: capítulo 4 del libro de David Mount.

**Archivo entregado:** `interpretacion_blast.md`

### 2.5 Ejercicio 3 — Alineamiento múltiple (MSA)

**Pedido:**
- Descargar en FASTA las secuencias de los **10 mejores hits** de BLAST.
- Realizar **MSA** con la secuencia query + esas 10 secuencias.
- Herramienta: MUSCLE local, MAFFT, o herramienta online si no hay instalación local.
- Interpretar el resultado del alineamiento.

**Input:** `blast_results.xml`, `query_best.fasta`  
**Output:** `msa_input.fasta`, `msa_output.afa`  
**Archivo entregado:** `interpretacion_msa.md`

---

## 3. Contexto biológico (para introducción y discusión)

### 3.1 Gen INS e insulina

- El gen **INS** está en el cromosoma **11p15.5** en humanos.
- Codifica la **preproinsulina**, que tras procesamiento post-traduccional da lugar a:
  - Péptido señal (eliminado)
  - Cadena B y cadena A de insulina (unidas por puentes disulfuro)
  - Péptido C (eliminado en la maduración)
- La insulina es secretada por las **células β del páncreas** y regula la **glucemia**.

### 3.2 Relación con diabetes mellitus

- **Diabetes mellitus tipo 1:** enfermedad autoinmune; destrucción de células β; no siempre hay mutación en INS.
- **Diabetes monogénica / MODY:** mutaciones en genes incluyendo INS pueden causar diabetes de inicio temprano.
- En OMIM se puede consultar la entrada del gen INS y fenotipos asociados a alteraciones de este locus.

### 3.3 Por qué NM_000207

- Es un **mRNA de referencia maduro** (RefSeq NM_): ya procesado, **sin intrones**.
- Permite traducir directamente nucleótidos → proteína sin errores por splicing.
- Es el transcripto estándar usado en la consigna como ejemplo para el gen INS.

### 3.4 Marcos de lectura

- Un mRNA puede leerse en **6 marcos**: 3 en el sentido 5'→3' (forward) y 3 en el complemento inverso (reverse).
- Solo **uno** suele contener el **ORF** (Open Reading Frame) real: inicia en codón AUG (Met) y termina en stop (UAA/UAG/UGA).
- Los marcos incorrectos producen proteínas sin sentido biológico o sin homología en bases de datos.

---

## 4. Qué hace el proyecto (descripción técnica)

### 4.1 Flujo general

```
OMIM (investigación) → NCBI Entrez → GenBank mRNA
    → Ex1: 6 traducciones FASTA
    → Ex2: BLASTp × 6 contra Swiss-Prot → elegir mejor marco
    → Ex3: descargar top 10 hits + MSA (MAFFT/MUSCLE)
    → Interpretaciones escritas + exposición oral
```

### 4.2 Scripts y funciones

| Archivo | Función |
|---------|---------|
| `fetch_data.py` | Descarga `NM_000207` desde NCBI Entrez en formato GenBank |
| `Ex1.py` | Lee GenBank, calcula 6 marcos, traduce cada uno a aminoácidos, exporta FASTA |
| `Ex2_a.py` | Ejecuta BLASTp remoto por cada secuencia del FASTA; compara E-values; selecciona mejor marco |
| `Ex3.py` | Parsea XML BLAST, descarga 10 proteínas por Entrez, alinea con MAFFT (o MUSCLE si está instalado) |
| `run_pipeline.py` | Orquestador multiplataforma (Linux, macOS, Windows) |
| `run_pipeline.sh` / `run_pipeline.bat` | Wrappers para Unix y Windows |
| `Ex2_local.py` | BLASTp local contra Swiss-Prot (`BLAST_MODE=local` o `both`) |
| `prepare_blast_db.py` | Descarga Swiss-Prot y crea base con `makeblastdb` |
| `platform_tools.py` | Detección de ejecutables según SO |
| `check_requirements.py` | Verificación local de BioPython, MSA, EMBOSS y BLAST local opcional |
| `setup.sh` / `setup.bat` | Crea entorno virtual e instala BioPython |

### 4.3 Detalle Ejercicio 1 (`Ex1.py`)

- Usa **BioPython** (`SeqIO`, `Seq`, `SeqRecord`).
- Para cada registro GenBank:
  - Genera 3 marcos forward (offset 0, 1, 2) y 3 reverse (complemento inverso, offset 0, 1, 2).
  - Recorta la secuencia a múltiplo de 3 y traduce con el código genético estándar.
  - Verifica coincidencia con la anotación **CDS** del GenBank (si existe).
- Escribe `frame_annotation.txt` con el marco que contiene el CDS.

### 4.4 Detalle Ejercicio 2 (`Ex2_a.py`)

- Usa `Bio.Blast.NCBIWWW.qblast("blastp", "swissprot", secuencia)`.
- Un BLAST remoto **por cada** secuencia en el FASTA (6 en total).
- Guarda XML individual en `blast_results/`.
- Ordena por **menor E-value** del mejor hit.
- Copia el mejor resultado a `blast_results.xml` y escribe en `query_best.fasta` la región CDS alineada por BLAST, recortando los aminoácidos extra del marco completo.

### 4.5 Detalle Ejercicio 3 (`Ex3.py`)

- Parsea `blast_results.xml` con `Bio.Blast.NCBIXML`.
- Extrae accession IDs Swiss-Prot de los 10 mejores alineamientos.
- Descarga secuencias completas con `Entrez.efetch(db="protein")`.
- Combina query + 10 hits en `msa_input.fasta` (11 secuencias).
- Ejecuta **MAFFT** (`mafft --auto`) o **MUSCLE** si está disponible.
- Genera `msa_output.afa`.

### 4.6 Cómo ejecutar

```bash
git clone https://github.com/apodgornyitba/tp_bio.git
cd tp_bio
./setup.sh
source .venv/bin/activate
export ENTREZ_EMAIL="email@ejemplo.com"
# macOS: brew install mafft
./run_pipeline.sh
```

**Tiempo estimado Ej. 2:** 30–60 minutos (6 BLAST remotos).

---

## 5. Resultados obtenidos (datos reales del experimento)

### 5.1 Ejercicio 1 — Traducción de marcos

| Marco | Longitud traducción | Stops (*) | ¿Coincide CDS? |
|-------|---------------------|-----------|----------------|
| Forward_Frame_1 | 155 aa | 4 | No |
| Forward_Frame_2 | 154 aa | 1 | No |
| **Forward_Frame_3** | **154 aa** | **1** | **Sí** |
| Reverse_Frame_1 | 155 aa | 1 | No |
| Reverse_Frame_2 | 154 aa | 2 | No |
| Reverse_Frame_3 | 154 aa | 7 | No |

- **CDS anotado en GenBank:** 110 aminoácidos (preproinsulina).
- **Marco identificado por CDS:** `Forward_Frame_3`.
- La traducción del marco 3 contiene la secuencia de preproinsulina (inicia con región señal MALWMR…).

### 5.2 Ejercicio 2 — BLAST por marco

| Marco | Mejor hit | E-value | Bit score | Identidad | ¿Correcto? |
|-------|-----------|---------|-----------|-----------|------------|
| Forward_Frame_1 | Sin hits | — | — | — | No |
| Forward_Frame_2 | Sin hits | — | — | — | No |
| **Forward_Frame_3** | Insulina humana P01308 | **3.36×10⁻⁷⁶** | **572** | **110/110 (100%)** | **Sí** |
| Reverse_Frame_1 | Sin hits | — | — | — | No |
| Reverse_Frame_2 | Sin hits | — | — | — | No |
| Reverse_Frame_3 | Sin hits | — | — | — | No |

**Marco seleccionado:** `NM_000207.3_Forward_Frame_3`

**Interpretación:** Solo el marco 3 produce una proteína con homología significativa contra Swiss-Prot. Los otros 5 marcos no tienen hits, lo que confirma que son traducciones incorrectas. El resultado coincide con la anotación CDS del GenBank.

### 5.3 Top 5 hits BLAST (marco correcto)

| Rank | Accession | Organismo | E-value | Bit score | Identidad |
|------|-----------|-----------|---------|-----------|-----------|
| 1 | P01308 | *Homo sapiens* (humano) | 3.36×10⁻⁷⁶ | 572 | 110/110 (100%) |
| 2 | Q8HXV2 | *Pongo pygmaeus* (orangután) | 2.07×10⁻⁷⁵ | 567 | 109/110 (99.1%) |
| 3 | P30410 | *Pan troglodytes* (chimpancé) | 3.11×10⁻⁷⁵ | 566 | 108/110 (98.2%) |
| 4 | P30406 | *Macaca fascicularis* (macaco) | 8.00×10⁻⁷⁴ | 557 | 108/110 (98.2%) |
| 5 | P30407 | *Chlorocebus aethiops* (mono verde) | 3.40×10⁻⁷³ | 553 | 107/110 (97.3%) |

**Interpretación biológica:** Todos los mejores hits son **preproinsulinas de primates**. La altísima identidad y E-values extremadamente bajos indican homología evolutiva fuerte, no coincidencia aleatoria. La insulina es una hormona altamente conservada por su función esencial en metabolismo de glucosa.

### 5.4 Significado de parámetros estadísticos BLAST

| Parámetro | Definición | Valor en nuestro experimento |
|-----------|------------|------------------------------|
| **E-value** | Número esperado de hits de igual calidad por azar | ~10⁻⁷⁶ → significancia extrema |
| **Bit score** | Medida de calidad del alineamiento (independiente del tamaño de la base) | 553–572 → muy alto |
| **% identidad** | Porcentaje de residuos idénticos en la región alineada | 97–100% |
| **Gaps** | Inserciones/deleciones en el alineamiento | Pocos en primates |
| **Query coverage** | Porcentaje de la query que participó del alineamiento | ~casi total (110 aa) |

### 5.5 Ejercicio 3 — MSA

**Herramienta usada:** MAFFT 7.526 (MUSCLE no disponible en Homebrew macOS; MAFFT es equivalente para la consigna).

**Secuencias alineadas (11 total):**

| # | ID | Especie / descripción |
|---|-----|----------------------|
| 0 | NM_000207.3_Forward_Frame_3_CDS | Query — humano (CDS recortada desde marco 3) |
| 1 | P01308 | Humano |
| 2 | Q8HXV2 | Orangután (*Pongo pygmaeus*) |
| 3 | P30410 | Chimpancé |
| 4 | P30406 | Macaco |
| 5 | P30407 | Mono verde |
| 6 | P01313 | Hámster chino (*Cricetulus longicaudatus*) |
| 7 | P01323 | Rata insulina-2 |
| 8 | P01322 | Rata insulina-1 |
| 9 | P01321 | Perro |
| 10 | P67972 | Mono nocturno (*Aotus trivirgatus*) |

**Regiones conservadas observadas:**
- Motivos de preproinsulina: MALWMR…, FVNQHLCGSHLVEALYLVCGERGFFYTPK…
- **Cisteínas (C)** alineadas en columnas fijas → puentes disulfuro entre cadenas A y B.
- Motivos GFFYTPK y SLYQLENYCN conservados entre primates.

**Regiones variables:**
- Query recortada: `query_best.fasta` queda restringida a la región CDS alineada por BLAST (aa 20-129 del marco completo), sin aminoácidos traducidos desde UTRs ni stop final.
- Roedores y hámster chino: más sustituciones y gaps que primates cercanos.

**Conclusión MSA:** Núcleo de insulina muy conservado entre especies del top 10; divergencia mayor en taxones más lejanos (roedores) sin perder función hormonal.

### 5.6 Ejercicio 4 — Análisis de Dominios con EMBOSS y PROSITE
Para caracterizar funcionalmente los marcos de lectura obtenidos en el Ejercicio 1, se realizó un análisis de dominios estructurado:
1. **Extracción y ORFs**: Se extrajo la secuencia completa de nucleótidos del transcripto (`NM_000207.3`) en formato FASTA (`NM_000207_nucleotides.fasta`). En la corrida final debe ejecutarse `getorf` de EMBOSS con `REQUIRE_EMBOSS=1` o `--require-emboss` para obtener los ORFs potenciales en ambas hebras.
2. **Descarga de Motivos**: Se descargó la base de datos de motivos de proteínas de referencia **PROSITE** (`prosite.dat` y `prosite.doc`) y se preparó el índice local usado por EMBOSS.
3. **Escaneo de Dominios con `patmatmotifs`**: Se escanearon los 14 ORFs generados por `getorf` nativo en búsqueda de firmas biológicas estructuradas.
4. **Modo de robustez**: El fallback Python queda disponible para desarrollo, pero no debe usarse como evidencia principal de entrega si la consigna exige EMBOSS instalado.

**Resultados obtenidos:**
La salida nativa de `patmatmotifs` reportó 3 hits de motivos PROSITE: una firma **INSULIN** y dos señales **AMIDATION**. El hallazgo biológicamente más relevante tuvo lugar en el ORF **`NM_000207.3_4`** (ORF forward que contiene la región codificante de la preproinsulina identificada por BLAST), donde se detectó la **firma de la familia de la insulina** (Insulin family signature, PROSITE ID: `INSULIN`, Accession: `PS00262`):
- **Coordenadas**: Residuos 114 a 128
- **Secuencia encontrada**: `C-C-T-S-I-C-S-L-Y-Q-L-E-N-Y-C`
- **Patrón PROSITE**: `C-C-{P}-{P}-x-C-[STDNEKPI]-x(3)-[LIVMFS]-x(3)-C`

Esta detección corrobora formal y funcionalmente que la secuencia es preproinsulina humana activa, validando de forma autónoma la predicción por homología de BLAST.

### 5.7 Ejercicio 5 — Diseño de Primers (parámetros de qPCR)
A partir del transcripto `NM_000207.3` (465 pb), el script aplica la variante patogénica configurada `NM_000207.3(INS):c.125T>C (p.Val42Ala)`, asociada a MODY10 (`rs886037863`), y prioriza primers/pares F-R cuyo amplicón flanquea esa posición. Esto alinea el diseño con el objetivo de detectar la variante del mensajero en pacientes, no solo el transcripto de referencia WT.
- **Parámetros configurados (`primer_config.json`)**:
  - Tamaño: 18 a 24 pares de bases.
  - Contenido GC: 50% mínimo, 60% máximo.
  - Extremos terminales (5' y 3'): Exclusión estricta de G o C (deben terminar/iniciar con A o T).
  - Temperatura de melting ($T_m$): $\le 67^\circ\text{C}$.
- **Resultados generados**:
  El reporte `primer_results/primers_report.txt` muestra 292 candidatos, 112 pares que flanquean la variante y 5 primers seleccionados bajo las restricciones de longitud, GC, extremos y Tm. El mejor par amplifica el rango 48-207 bp, que contiene la variante en el nt 184.

**Análisis físico-químico:**
- **Especificidad de Extremos**: Ningún primer termina en G o C (todos inician/finalizan en A/T), lo cual evita la formación de dímeros estables de primers ("primer-dimers") y falsos cebados inespecíficos en sus extremos 3'.
- **Uniformidad de Tm**: Las $T_m$ están óptimamente calibradas alrededor de 60°C, asegurando una hibridación altamente eficiente en protocolos estandarizados de ciclado térmico.
- **Parametrización**: Los parámetros de diseño y la variante se definen de forma dinámica desde un archivo JSON (`primer_config.json`), permitiendo adaptar el análisis a otro gen, otra variante o condiciones de laboratorio modificadas.

---

## 6. Archivos del repositorio

### 6.1 Scripts

- `fetch_data.py` — Descarga GenBank desde NCBI
- `Ex1.py` — Traduce los 6 marcos de lectura
- `Ex2_a.py` — Corre BLASTp remoto contra Swiss-Prot
- `Ex2_local.py` — Corre BLASTp local contra Swiss-Prot
- `Ex3.py` — Parsea hits, los descarga y corre MSA (MAFFT)
- `Ex4.py` — Realiza traducción de ORFs (getorf) y escaneo de PROSITE (patmatmotifs)
- `Ex5.py` — Diseña primers sobre la variante configurada y prioriza pares que la flanquean
- `run_pipeline.py` — Orquestador multiplataforma de la pipeline (Paso 0 a Paso 5)
- `run_pipeline.sh` / `run_pipeline.bat` — Wrappers de ejecución
- `check_requirements.py` — Verifica dependencias locales antes de correr la entrega final
- `setup.sh` / `setup.bat` — Scripts de configuración e instalación

### 6.2 Datos generados

- `NM_000207.gbk` — Secuencia GenBank de referencia
- `NM_000207_frames.fasta` — 6 traducciones de marcos de lectura
- `frame_annotation.txt` — Anotación del marco CDS biológico
- `blast_results/` — XML individuales por marco y resumen de BLAST
- `blast_results.xml` — BLAST del mejor marco (Forward 3)
- `query_best.fasta` — Secuencia CDS recortada desde el mejor marco para MSA
- `msa_input.fasta` — 11 secuencias combinadas para alineamiento
- `msa_output.afa` — Alineamiento múltiple generado por MAFFT/MUSCLE
- `prosite.dat` — Base de datos local de PROSITE (~24 MB)
- `emboss_results/` — Secuencia de nucleótidos, ORFs traducidos y dominios detectados
- `primer_results/` — JSON de primers y reporte detallado en texto plano
- `pipeline.log` — Log completo de ejecución

### 6.3 Documentación

- `README.md` — Instrucciones de ejecución y configuración
- `consigna.md` — Consigna original de la cátedra
- `primer_config.json` — Configuración paramétrica para diseño de primers
- `interpretacion_blast.md` — Interpretación del BLAST
- `interpretacion_msa.md` — Interpretación del MSA
- `INFORME_PROYECTO.md` — Este documento de reporte integrado

---

## 7. Cumplimiento de la consigna (checklist)

| Requisito | Estado | Notas |
|-----------|--------|-------|
| Elegir enfermedad/genes en OMIM | ✅ | Diabetes + INS |
| Script Ej. 1 GenBank → FASTA 6 marcos | ✅ | `Ex1.py` |
| mRNA maduro sin intrones | ✅ | NM_000207 |
| Script Ej. 2 BLAST por secuencia | ✅ | `Ex2_a.py`, 6 BLAST |
| BLAST remoto | ✅ | NCBI qblast |
| BLAST local | Opcional | `Ex2_local.py` + `prepare_blast_db.py` (`BLAST_MODE=both`) si se busca puntaje extra |
| Interpretación BLAST (2.b) | ✅ | `interpretacion_blast.md` |
| Ej. 3 top 10 + MSA | ✅ | `Ex3.py` + MAFFT |
| Interpretación MSA | ✅ | `interpretacion_msa.md` |
| Ej. 4 EMBOSS local + PROSITE | ✅ | Corrida final con `--require-emboss`: getorf + patmatmotifs nativos |
| Ej. 5 Diseño de primers parametrizado | ✅ | `Ex5.py` aplica la variante configurada y prioriza pares que la flanquean |
| Automatización Bash / Python | ✅ | `run_pipeline.py` y wrappers `.sh` / `.bat` |
| Descripción y cómo ejecutar | ✅ | `README.md` |
| Exposición 10 min investigación | ✅ | Borrador de slides y guion preparados; falta adaptar integrantes y ensayar |

---

## 8. Limitaciones y decisiones de diseño

1. **BLAST local opcional:** Requiere instalar BLAST+ y ejecutar `prepare_blast_db.py` (~150 MB).
2. **MAFFT en lugar de MUSCLE:** En macOS Homebrew no hay paquete `muscle`; MAFFT cumple el mismo rol para MSA.
3. **Traducción del marco completo:** Ex1 traduce todo el marco (con stops `*`), no solo el ORF más largo; el marco 3 aun así alinea 100% con insulina en la región del CDS.
4. **Query MSA vs Swiss-Prot:** `query_best.fasta` se recorta a la región CDS alineada por BLAST antes del MSA, evitando el ruido del marco completo.
5. **Fallbacks de Robustez:** El fallback Python de EMBOSS queda como soporte de desarrollo; para la entrega final usar `REQUIRE_EMBOSS=1` para no aceptar resultados sin EMBOSS nativo.
6. **Filtro estricto de extremos terminales de Primers:** La restricción de no poseer G o C en extremos 5' y 3' se implementó con un filtro lógico bidireccional, garantizando cebadores con menor propensión al hairpin y homodímeros.
7. **Diseño sobre variante:** `primer_config.json` define `rs886037863` (`c.125T>C`), que se aplica al transcripto antes de buscar primers y pares F/R flanqueantes.

---

## 9. Puntos para la exposición oral (10 min, sin código)

1. Qué es OMIM y qué se encontró sobre diabetes e INS.
2. Función biológica de la insulina y del gen INS.
3. Por qué se usó NM_000207 (mRNA maduro).
4. Concepto de marcos de lectura y por qué hay 6.
5. Resultado BLAST: solo Forward_Frame_3 da insulina.
6. Hits a primates: conservación evolutiva e interpretación biológica del E-value.
7. MSA: regiones conservadas (cisteínas para puentes disulfuro, cadenas A/B).
8. **Análisis de Dominios (EMBOSS)**: Confirmación de la firma biológica de la insulina (`PS00262`) en la secuencia traducida de forma autónoma.
9. **Diseño de Primers**: Racional físico-químico detrás del diseño y foco en la variante `rs886037863` (amplicones que la flanquean, evitar G/C en extremos, calibración de Tm y GC, parametrización JSON).
10. Conclusión: validación del dogma central mediante herramientas bioinformáticas integradas y automatizadas.

---

## 10. Estructura sugerida para el informe final

Pedir al asistente que genere el informe con esta estructura:

1. **Carátula** — Materia, grupo, fecha, gen/enfermedad
2. **Resumen / Abstract** — 150–250 palabras
3. **Introducción** — Diabetes, INS, OMIM, objetivos del TP
4. **Marco teórico** — GenBank, FASTA, marcos de lectura, BLAST, MSA, E-value
5. **Materiales y métodos** — Software, bases de datos, pipeline, accession numbers
6. **Resultados** — Tablas y descripción de Ex1, Ex2, Ex3 (usar sección 5 de este doc)
7. **Discusión** — Validez del marco 3, conservación de insulina, limitaciones
8. **Conclusiones** — Respuestas a los objetivos del TP
9. **Referencias** — NCBI, OMIM, BioPython, Swiss-Prot, MAFFT
10. **Anexos** — Comandos de ejecución, listado de archivos

---

## 11. Referencias sugeridas para bibliografía

- NCBI Gene: INS — https://www.ncbi.nlm.nih.gov/gene/3630
- NCBI Nucleotide: NM_000207 — https://www.ncbi.nlm.nih.gov/nuccore/NM_000207
- OMIM — https://omim.org/
- BioPython — https://biopython.org/
- BLAST Help — https://www.ncbi.nlm.nih.gov/books/NBK52637/
- Swiss-Prot — https://www.uniprot.org/
- MAFFT — https://mafft.cbrc.jp/alignment/software/

---

## 12. Prompt sugerido para Claude

```
A partir del archivo INFORME_PROYECTO.md, redactá un informe académico completo 
en español para el TP de Bioinformática (Partes 1 y 2). Incluí:

- Tono formal de informe universitario
- Todas las secciones de la estructura sugerida (sección 10)
- Tablas de resultados con los datos numéricos reales de BLAST, MSA y diseño de Primers
- Hallazgos del análisis de dominios EMBOSS (motivo de la firma de la insulina PS00262)
- Interpretación biológica integrada (no solo técnica)
- Mención del cumplimiento total de la consigna
- Bibliografía en formato APA o Vancouver

No inventes datos: usá solo los resultados reales documentados en el archivo.
Extensión aproximada: 10–15 páginas.
```

---

*Documento generado para el repositorio tp_bio — ITBA Bioinformática.*
