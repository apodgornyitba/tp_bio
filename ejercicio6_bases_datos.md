# Ejercicio 6 — Trabajo con Bases de Datos Biológicas

**Gen:** INS (insulina)  
**Enfermedad investigada:** Diabetes mellitus (incl. MODY10 por mutaciones en INS)  
**Organismo:** *Homo sapiens*  
**Fecha de consulta:** Mayo 2026

---

## a) Entrada en NCBI Gene (Entrez)

### Link Entrez
**https://www.ncbi.nlm.nih.gov/gene/3630**

| Campo | Valor |
|-------|-------|
| **Símbolo** | INS |
| **Nombre completo** | insulin |
| **ID Entrez Gene** | 3630 |
| **HGNC** | HGNC:6081 |
| **Ensembl** | ENSG00000254647 |
| **OMIM** | 176730 |
| **Cromosoma** | 11p15.5 (GRCh38: 2,159,779–2,161,221, hebra complementaria) |
| **Tipo** | Gen codificante de proteína (RefSeq REVIEWED) |

### ¿Qué hace la proteína?

INS codifica la **preproinsulina**, hormona peptídica clave en el metabolismo de **carbohidratos y lípidos**. Tras eliminar el péptido señal, la preproinsulina se procesa en:
- Cadena **B** y cadena **A** (unidas por puentes disulfuro → insulina madura).
- **Péptido C** (eliminado en la maduración).

La insulina se **secreta desde las células β del páncreas** y, al unirse al **receptor de insulina (INSR)**, estimula la captación de glucosa y múltiples vías de señalización metabólica.

### ¿Por qué la elegimos?

- Es el **producto directo del gen central** en nuestro TP de bioinformática.
- Está **asociado clínicamente a diabetes mellitus**: desde diabetes monogénica (**MODY10**) y diabetes neonatal permanente, hasta su rol central en la fisiología de la glucemia.
- Permite integrar todos los ejercicios del pipeline (traducción, BLAST, MSA, EMBOSS, primers) en un caso biológico coherente.
- Es un ejemplo clásico citado en la consigna del TP (gen INS / diabetes).

### Expresión (NCBI Gene)
Expresión **restringida al páncreas** (RPKM ~671.7 en tejido pancreático), coherente con su función endocrina.

---

## b) Genes / proteínas homólogas (HomoloGene, NCBI Orthologs y Ensembl)

### NCBI — HomoloGene y Orthologs

> **Nota:** La base **HomoloGene** fue **retirada por NCBI** (2023). La información de homología actualizada está en **[NCBI Orthologs](https://www.ncbi.nlm.nih.gov/gene/3630/ortholog/)** y en **NCBI Datasets Gene Table**.

- **Link Orthologs:** https://www.ncbi.nlm.nih.gov/gene/3630/ortholog/
- **OrthoDB:** https://www.orthodb.org/?ncbi=3630
- Históricamente, el cluster HomoloGene **30951** agrupaba insulinas de múltiples vertebrados.

**Organismos con ortólogos típicos de INS** (confirmados también por nuestro BLAST del TP):
- *Homo sapiens*, *Pan troglodytes*, *Gorilla gorilla*, *Macaca fascicularis*
- *Rattus norvegicus* (Ins1, Ins2), *Canis lupus*, otros vertebrados

### Ensembl

- **Link:** https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000254647
- **254 ortólogos** reportados en Ensembl.
- **2 parálogos** (incluye relación con locus vecinos como **IGF2** / gen de readthrough **INS-IGF2**).

### ¿Qué tan común es el gen? ¿En qué taxones?

| Aspecto | Conclusión |
|---------|------------|
| **Conservación** | Muy alta entre **vertebrados** (especialmente mamíferos). |
| **Taxones** | Presente en **vertebrados**; no es un gen universal (no se encuentra en bacterias ni plantas). |
| **Función conservada** | Regulación de metabolismo energético / glucemia. |

### Diferencias entre HomoloGene (NCBI) y Ensembl

| Criterio | NCBI (HomoloGene → Orthologs) | Ensembl |
|----------|----------------------------------|---------|
| **Enfoque** | Clusters de homología basados en proteínas anotadas; ahora Orthologs con algoritmo actualizado | Anotación de genoma de referencia + árbol filogenético Ensembl Compara |
| **Cobertura** | Menos especies en la interfaz clásica; Datasets amplía el acceso | **254 ortólogos**, integración con GENCODE/MANE |
| **Parálogos** | Menos destacados en la vista Gene | Reporta explícitamente **2 parálogos** |
| **Actualización** | HomoloGene congelada; usar Orthologs/Datasets | Actualización continua (release 115, 2025) |

**Opinión del grupo:** Ensembl suele ser **más completo para transcriptos y ortología filogenética**; NCBI Gene es **más práctico para integrar RefSeq, ClinVar, OMIM y GEO** en una sola ficha.

---

## c) Transcriptos y splicing alternativo

### NCBI Gene / RefSeq

**Link transcriptos:** https://www.ncbi.nlm.nih.gov/datasets/tables/genes/?table_type=transcripts&ids=3630

| Transcripto | Descripción |
|-------------|-------------|
| **NM_000207.3** | Variante 1 — **MANE Select** / isoforma de referencia (la usamos en el TP) |
| NM_001185097.2 | Variante 2 — difiere en **5' UTR** |
| NM_001185098.2 | Variante 3 — difiere en **5' UTR** |
| NM_001291897.2 | Variante 4 — difiere en **5' UTR** |

**Proteína resultante:** **Todas codifican la misma preproinsulina** (NP_000198.1 → UniProt **P01308**).  
Las diferencias son principalmente en **regiones no traducidas (5' UTR)**, no en isoformas proteicas distintas.

### Ensembl

- Reporta **5 transcriptos** (variantes de splicing).
- Transcripto canónico / MANE: **ENST00000381330** (corresponde a NM_000207).
- Un transcripto marcado como **non-stop decay** (posiblemente no funcional o poco soportado).

### ¿Cuántos splicing alternativos? ¿Se expresan? ¿Funciones alternativas?

| Pregunta | Respuesta |
|----------|-----------|
| **¿Cuántos?** | NCBI: **4 mRNA** RefSeq; Ensembl: **5** variantes de splicing. |
| **¿Se expresan?** | La isoforma principal (**NM_000207**) está muy respaldada (MANE Select, APPRIS P1, TSL 1). Las otras variantes tienen soporte en 5' UTR. |
| **¿Funciones alternativas?** | **No** hay isoformas proteicas alternativas conocidas; el splicing afecta mayormente **regulación transcricional** (UTRs), no dominios proteicos distintos. |

### ¿Por qué difieren NCBI y Ensembl?

- **Anotaciones independientes** (RefSeq vs GENCODE).
- **Criterios de evidencia** distintos para incluir transcriptos de bajo soporte.
- Ensembl incluye **1 variante adicional** y anotaciones de calidad (APPRIS, TSL) más detalladas.

**¿Cuál es más precisa?** Para **uso clínico y nuestro TP**, **NM_000207 (NCBI MANE Select)** es la referencia más aceptada. Ensembl es superior para **comparar splicing a nivel genómico** en estudios de anotación.

---

## d) Interacciones proteína–proteína

### NCBI Gene — Interactions

**Link:** https://www.ncbi.nlm.nih.gov/gene/3630 (sección *Interactions*)

La tabla principal de interacciones génicas de NCBI para INS está **poco poblada** en la ficha estándar. Sin embargo, en **GO (IPI — interacción física)** y en la base **HIV-1 Human Interaction** se documentan interacciones indirectas:

| Interactor | Tipo de interacción |
|------------|---------------------|
| **INSR** (receptor de insulina) | Unión del ligando insulina → señalización downstream |
| **Proteína Vpr (HIV-1)** | Vpr de HIV antagoniza vías de señalización de insulina (efecto patológico) |
| **Auto-interacción** (GO:0042802) | Unión proteína idéntica (dimerización / agregación) |

### UniProt — P01308 (Insulina humana)

**Link:** https://www.uniprot.org/uniprotkb/P01308

| Aspecto | Información UniProt |
|---------|---------------------|
| **Subunidad** | Heterodímero cadena A + cadena B (puentes disulfuro) |
| **Interacción principal** | **Receptor de insulina (INSR)** |
| **Función** | Ligando endocrino; reduce glucemia; aumenta permeabilidad celular a monosacáridos |
| **Localización** | Secretada al espacio extracelular desde células β pancreáticas |

### Comparación NCBI vs UniProt

| Fuente | Fortaleza | Interacciones destacadas |
|--------|-----------|--------------------------|
| **NCBI Gene** | Integración con literatura, GO, ClinVar | INSR (vía GO), interacciones virales (Vpr-HIV) |
| **UniProt** | Curación manual de proteínas | INSR, estructura del complejo, detalles bioquímicos |

**Interacciones únicas por base:**
- **Solo NCBI:** interacciones anotadas en contexto **HIV-1 / Vpr** (efecto sobre vías de insulina).
- **Solo UniProt:** detalle estructural del **complejo insulina–receptor** y anotaciones farmacológicas (Humulin, Novolin).

**Patrón general:** La insulina **no forma redes de interacción masivas** como proteínas de señalización intracelular; su función central es como **ligando del INSR**, iniciando cascadas de fosforilación (PI3K/AKT, etc.).

---

## e) Gene Ontology (GO)

**Fuente:** NCBI Gene / GOA — https://www.ncbi.nlm.nih.gov/gene/3630  
**AmiGO:** https://amigo.geneontology.org/amigo/gene/ENSG00000254647

### Componente celular (CC)

| Término GO | Significado |
|------------|-------------|
| **GO:0005615 — extracellular space** | Insulina madura secretada |
| **GO:0005788 — ER lumen** | Síntesis y plegamiento de preproinsulina |
| **GO:0005796 — Golgi lumen** | Procesamiento y empaquetado en gránulos |
| **GO:0034774 — secretory granule lumen** | Gránulos secretores de células β |
| **GO:0030133 — transport vesicle** | Tráfico hacia la membrana plasmática |

### Proceso biológico (BP)

| Término GO | Significado |
|------------|-------------|
| **GO:0042593 — glucose homeostasis** | Homeostasis de glucosa (función central) |
| **GO:0008286 — insulin receptor signaling pathway** | Vía de señalización del receptor |
| **GO:0006006 — glucose metabolic process** | Metabolismo de glucosa |
| **GO:0045721 — negative regulation of gluconeogenesis** | Supresión de gluconeogénesis hepática |
| **GO:0045821 — positive regulation of glycolytic process** | Estimulación de glucólisis |
| **GO:0046889 — positive regulation of lipid biosynthetic process** | Efectos sobre lípidos |

### Función molecular (MF)

| Término GO | Significado |
|------------|-------------|
| **GO:0005179 — hormone activity** | Actividad hormonal |
| **GO:0005158 — insulin receptor binding** | Unión al receptor INSR |
| **GO:0042802 — identical protein binding** | Dimerización / agregación |

**Resumen:** INS es una **hormona secretada** (extracelular) sintetizada en **ER/Golgi** del páncreas, con función molecular de **ligando del INSR** y procesos biológicos centrados en **homeostasis de glucosa**.

---

## f) Vías metabólicas (pathways)

### Reactome

| Pathway | ID | Rol de INS |
|---------|-----|------------|
| **Regulation of insulin secretion** | [R-HSA-422356](https://reactome.org/content/detail/R-HSA-422356) | INS es el producto final secretado tras estimulación por glucosa |
| **Insulin processing** | [R-HSA-264876](https://reactome.org/content/detail/R-HSA-264876) | Síntesis de preproinsulina, puentes disulfuro, clivaje a insulina + péptido C |
| **Exocytosis of Insulin** | [R-HSA-265166](https://reactome.org/content/detail/R-HSA-265166) | Exocitosis de gránulos de insulina (dependiente de Ca²⁺) |
| **Insulin receptor signalling** | (vías downstream de INSR) | Captación de glucosa, metabolismo |

**Link proteína:** https://reactome.org/content/schema/instance/browser/uniprot:P01308

### KEGG

| Pathway | ID | Relación |
|---------|-----|----------|
| **Type II diabetes mellitus** | hsa04930 | Deficiencia relativa de insulina / resistencia periférica |
| **Maturity onset diabetes of the young** | hsa04950 | Incluye formas monogénicas como **MODY10 (INS)** |
| **Insulin secretion** | hsa04911 | Secreción desde células β |
| **Insulin signaling pathway** | hsa04910 | Señalización downstream del receptor |

**Link KEGG gen:** https://www.genome.jp/dbget-bin/www_bget?hsa:3630

### Interpretación

INS participa en la **cascada de secreción de insulina** (desde glucosa → ATP → cierre de canales KATP → Ca²⁺ → exocitosis) y es el **ligando central** de la **vía de señalización de insulina** que regula metabolismo de glucosa en hígado, músculo y tejido adiposo. En KEGG, mutaciones en INS se vinculan directamente a **MODY** y formas de **diabetes monogénica**.

---

## g) dbSNP / ClinVar — Variante asociada a la patología

### Búsqueda en ClinVar

**Link:** https://www.ncbi.nlm.nih.gov/clinvar/?term=INS[gene]

### Variante elegida para el análisis

| Campo | Valor |
|-------|-------|
| **Nombre HGVS** | NM_000207.3(INS):c.125T>C (p.Val42Ala) |
| **dbSNP (rs)** | **rs886037863** |
| **Proteína** | V42A (valina → alanina en preproinsulina) |
| **Enfermedad** | **Maturity-onset diabetes of the young, type 10 (MODY10)** |
| **OMIM fenotipo** | 613370 |
| **Significado clínico** | **Likely pathogenic / Pathogenic** (según submitter) |
| **Cromosoma** | 11p15.5 |
| **Link ClinVar** | https://www.ncbi.nlm.nih.gov/clinvar/RCV000240176/ |

### ¿Qué variante es?

Es una **mutación missense** en el gen **INS** que altera un aminoácido de la preproinsulina (V42A). Las mutaciones en INS pueden causar **MODY10** por defectos en el procesamiento de proinsulina, estrés del retículo endoplasmico y apoptosis de células β.

### Frecuencia poblacional

- Es una variante **rara** (típica de enfermedades monogénicas).
- En **gnomAD / dbSNP** suele tener frecuencia **< 0.001%** en población general (variante ultra-rara; consultar dbSNP rs886037863 para allele frequency actualizada).
- Algunos casos reportados son **de novo** o **paternos** (ClinVar: submitters clínicos).

### ¿Qué grupo étnico parece más afectado?

MODY10 por mutaciones en INS está descrito en **múltiples etnias** sin predominio claro en las bases públicas; hay reportes en **Europa, Asia y Medio Oriente** (ej. submitters de Italia, China en ClinVar). **No hay un grupo étnico único** — es una enfermedad monogénica rara distribuida globalmente.

### Otras variantes relevantes en INS (referencia)

| Variante | rs | Fenotipo |
|----------|-----|----------|
| c.130G>A (p.Gly44Arg) | rs765512575 | MODY10 (significado clínico discutido) |
| c.188T>C (p.Val63Ala) | — | Reportada en ClinVar |
| c.193C>T (p.Gln65Ter) | — | Variante stop (pérdida de función) |

### Bases complementarias consultadas

| Base | Link | Uso |
|------|------|-----|
| **OMIM — gen INS** | https://omim.org/entry/176730 | Fenotipos MODY10, PNDM, hiperproinsulinemia |
| **OMIM — MODY10** | https://omim.org/entry/613370 | Fenotipo clínico |
| **ClinVar** | https://www.ncbi.nlm.nih.gov/clinvar/?term=INS[gene] | Significado clínico de variantes |
| **GeneCards — INS** | https://www.genecards.org/cgi-bin/carddisp.pl?gene=INS | Resumen integrado |
| **MedlinePlus** | https://medlineplus.gov/genetics/gene/ins | Información para pacientes |
| **GHR** | https://medlineplus.gov/genetics/condition/maturity-onset-diabetes-of-the-young/ | MODY en general |

---

## Conclusión del Ejercicio 6

El gen **INS** es un locus **altamente conservado en vertebrados**, con expresión pancreática específica y anotaciones GO/pathway coherentes con su rol en **homeostasis de glucosa**. Las bases **NCBI** y **Ensembl** concuerdan en la proteína codificada pero difieren levemente en el conteo de transcriptos por criterios de anotación. Las interacciones funcionales más relevantes son con el **receptor de insulina (INSR)**. En **ClinVar**, mutaciones missense como **rs886037863 (V42A)** están asociadas a **MODY10**, vinculando directamente nuestro gen de estudio con la **diabetes mellitus monogénica** investigada en el TP.

---

## Links rápidos (para slides / informe)

| Recurso | URL |
|---------|-----|
| NCBI Gene INS | https://www.ncbi.nlm.nih.gov/gene/3630 |
| Ensembl INS | https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000254647 |
| UniProt P01308 | https://www.uniprot.org/uniprotkb/P01308 |
| OMIM INS | https://omim.org/entry/176730 |
| ClinVar INS | https://www.ncbi.nlm.nih.gov/clinvar/?term=INS[gene] |
| Reactome — Insulin processing | https://reactome.org/content/detail/R-HSA-264876 |
| KEGG hsa:3630 | https://www.genome.jp/dbget-bin/www_bget?hsa:3630 |
