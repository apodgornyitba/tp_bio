# Ejercicio 3 – Interpretación del alineamiento múltiple (MSA)

**Herramienta:** MAFFT (alternativa a MUSCLE; misma función)  
**Input:** Query CDS recortada desde `Forward_Frame_3` + 10 mejores hits de BLAST  
**Output:** `msa_output.afa`

## 1. Secuencias incluidas

| # | ID | Organismo |
|---|-----|-----------|
| 0 | NM_000207.3_Forward_Frame_3_CDS | *Homo sapiens* (query CDS recortada por BLAST) |
| 1 | P01308 | Humano |
| 2 | Q8HXV2 | Orangután (*Pongo pygmaeus*) |
| 3 | P30410 | Chimpancé |
| 4 | P30406 | Macaco |
| 5 | P30407 | Mono verde |
| 6 | P01313 | Hámster chino (*Cricetulus longicaudatus*) |
| 7 | P01323 | Rata (insulina-2) |
| 8 | P01322 | Rata (insulina-1) |
| 9 | P01321 | Perro |
| 10 | P67972 | *Aotus trivirgatus* (mono nocturno) |

## 2. Regiones conservadas

En `msa_output.afa` se observan bloques con **alta identidad** en:

- **Péptido señal y región central de la preproinsulina** (MALWMR…, FVNQHLCGSHLVEALYLVCGERGFFYTPK…).
- **Residuos de cisteína (C)** alineados en columnas fijas → forman **puentes disulfuro** entre cadenas A y B en la insulina madura.
- Motivos **GFFYTPK** y **SLYQLENYCN** muy conservados entre primates.

Estas regiones son esenciales para el plegamiento y la actividad hormonal.

## 3. Regiones variables y gaps

- **Query recortada:** `query_best.fasta` se genera desde la región alineada por BLAST (aa 20-129 del marco completo), eliminando aminoácidos traducidos desde UTRs y el stop final.
- **Rata y hámster chino:** mayor número de sustituciones y pequeños gaps respecto a primates.
- **Mono nocturno y roedores:** muestran más divergencia que los primates cercanos, pero conservan el núcleo funcional.

## 4. Conclusión

El MSA muestra que las insulinas del top 10 de BLAST comparten un **núcleo altamente conservado**, especialmente en las regiones que dan origen a las cadenas A y B. Las diferencias en roedores y especies más lejanas reflejan **divergencia evolutiva** sin perder la función. Esto apoya que el gen **INS** codifica una proteína fundamental para la regulación de la glucemia en vertebrados.
