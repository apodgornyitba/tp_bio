# Guion presentación Ejercicio 7 — INS y diabetes

**Archivo de slides:** `outputs/manual-20260614-tp7/presentations/tp-bio-ejercicio7/output/tp-bio-ins-diabetes-ejercicio7.pptx`  
**Duración objetivo:** 9:20–9:40 min, dejando margen para cambio de orador o preguntas.  
**Regla de la consigna:** explicar la investigación y resultados, no hacer una recorrida del código.

## Reparto sugerido

| Bloque | Slides | Tiempo | Foco |
|--------|--------|--------|------|
| Apertura biológica | 1–2 | 2:00 | Diabetes, OMIM, gen INS, transcripto NM_000207.3 |
| Evidencia computacional | 3–6 | 4:40 | Flujo, BLAST, MSA, EMBOSS/PROSITE |
| Variante y cierre clínico | 7–8 | 2:40 | Primers, rs886037863, bases de datos, conclusión |

## Slide por slide

### 1. Tesis del trabajo — 1:00

Mensaje: el TP no es una suma de scripts; es una investigación guiada por una pregunta biológica. Se eligió INS porque conecta diabetes, insulina, transcriptos, homología, dominios y variantes patogénicas.

Decir:
- "Partimos de diabetes mellitus y del gen INS."
- "El objetivo fue sostener esa elección desde la secuencia hasta una variante detectable por primers."
- "La presentación resume la evidencia, no la implementación."

### 2. Elección biológica — 1:00

Mensaje: INS/NM_000207.3 cumple la aclaración de la consigna porque es mRNA maduro, sin intrones, y se usa el mismo gen en ejercicios 1 y 5.

Decir:
- "NM_000207.3 es RefSeq de mRNA maduro."
- "La variante c.125T>C, rs886037863, permite conectar el diseño de primers con MODY10."
- "Esto mantiene coherencia entre OMIM, GenBank, dbSNP/ClinVar y el pipeline."

### 3. Flujo general — 1:10

Mensaje: cada etapa produce la evidencia necesaria para la siguiente.

Decir:
- "Ejercicio 1 traduce seis marcos."
- "Ejercicio 2 usa BLAST para distinguir el marco real."
- "Ejercicio 3 compara conservación."
- "Ejercicio 4 valida dominio con EMBOSS/PROSITE."
- "Ejercicio 5 diseña primers sobre la variante."
- "Ejercicio 6 conecta la proteína con bases de datos biológicas."

### 4. BLAST y marco correcto — 1:15

Mensaje: solo `Forward_Frame_3` tiene homología fuerte con insulina.

Decir:
- "El mejor hit fue P01308, insulina humana."
- "El E-value fue 3.36e-76, prácticamente incompatible con azar."
- "La identidad fue 110/110, 100% en la región alineada."
- "Los otros marcos no generaron hits comparables, por eso se elige este para MSA."

### 5. MSA — 1:10

Mensaje: la conservación no es uniforme; se concentra en regiones funcionales de insulina.

Decir:
- "El alineamiento contiene query más top 10 BLAST."
- "Primates aparecen con identidad muy alta, como se espera por cercanía evolutiva."
- "Las cisteínas conservadas son importantes por puentes disulfuro."
- "Roedores y especies más lejanas divergen más, pero conservan el núcleo funcional."

### 6. EMBOSS / PROSITE — 1:05

Mensaje: hay una validación independiente de familia proteica, no solo homología por BLAST.

Decir:
- "`getorf` generó 14 ORFs."
- "`patmatmotifs` detectó la firma INSULIN de PROSITE."
- "La firma PS00262 apareció en `NM_000207.3_4`, residuos 114–128."
- "Esto confirma funcionalmente que la secuencia contiene una región característica de insulina."

### 7. Primers y variante — 1:20

Mensaje: los primers no se diseñaron genéricamente; se priorizaron pares que flanquean la variante patogénica.

Decir:
- "La variante aplicada fue c.125T>C, p.Val42Ala, rs886037863."
- "Se respetaron 18–24 bp, GC 50–60%, sin G/C terminal y Tm menor o igual a 67 C."
- "El reporte encontró 292 candidatos y 112 pares que flanquean la variante."
- "El mejor amplicón mostrado es 48–207, de 160 bp, conteniendo la posición variante."

### 8. Cierre con bases de datos — 1:30

Mensaje: las bases de datos explican por qué los resultados importan biológicamente.

Decir:
- "NCBI/OMIM dan el marco clínico y genético."
- "Ensembl muestra conservación y transcriptos."
- "UniProt/GO describen actividad hormonal y unión a INSR."
- "Reactome/KEGG ubican la insulina en procesamiento, secreción y señalización."
- "dbSNP/ClinVar conecta rs886037863 con MODY10."
- "Conclusión: cada etapa apoya la misma historia, desde mRNA maduro hasta variante clínicamente relevante."

## Checklist antes de exponer

- Reemplazar o agregar nombres de integrantes si corresponde.
- Ensayar una vez con cronómetro y cortar detalles si pasan de 10 minutos.
- No explicar funciones internas de Python salvo que pregunten.
- Tener a mano `INFORME_PROYECTO.md`, `ejercicio6_bases_datos.md`, `interpretacion_blast.md` y `primer_results/primers_report.txt` para responder preguntas.
- Antes de entrega final, correr el preflight con el PATH que expone EMBOSS/BLAST:

```bash
env PATH=/Users/martin.zahnd/.conda/envs/tp_bio/bin:/opt/homebrew/bin:/usr/bin:/bin .venv/bin/python check_requirements.py
```
