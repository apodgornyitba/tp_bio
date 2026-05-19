# Ejercicio 3 – Interpretación del alineamiento múltiple (MSA)

**Herramienta:** MAFFT (alternativa a MUSCLE; misma función)  
**Input:** Query (Forward_Frame_3) + 10 mejores hits de BLAST  
**Output:** `msa_output.afa`

## 1. Secuencias incluidas

| # | ID | Organismo |
|---|-----|-----------|
| 0 | NM_000207.3_Forward_Frame_3 | *Homo sapiens* (query) |
| 1 | P01308 | Humano |
| 2 | Q8HXV2 | Gorila |
| 3 | P30410 | Chimpancé |
| 4 | P30406 | Macaco |
| 5 | P30407 | Mono verde |
| 6 | P01313 | Loris |
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

- **Extremo N-terminal de la query:** la traducción del marco completo incluye aminoácidos extra antes de Met (diferencia con Swiss-Prot que empieza en Met1); aparecen gaps en esa zona.
- **Rata y loris:** mayor número de sustituciones y pequeños gaps respecto a primates.
- **Cola C-terminal de la query:** incluye residuos extra por traducir el marco completo hasta el final del mRNA (no solo el CDS de 110 aa).

## 4. Conclusión

El MSA muestra que las insulinas del top 10 de BLAST comparten un **núcleo altamente conservado**, especialmente en las regiones que dan origen a las cadenas A y B. Las diferencias en roedores y especies más lejanas reflejan **divergencia evolutiva** sin perder la función. Esto apoya que el gen **INS** codifica una proteína fundamental para la regulación de la glucemia en vertebrados.
