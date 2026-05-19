# Ejercicio 2.b – Interpretación del BLAST

**Gen:** INS (insulina humana)  
**Transcripto:** NM_000207.3  
**Base de datos:** Swiss-Prot (curada)  
**Programa:** blastp (remoto, NCBI)

## 1. Objetivo

Comparar los **6 marcos de lectura** del Ejercicio 1 y determinar cuál corresponde a la preproinsulina real. El marco correcto debe dar hits de **insulina** con **E-value muy bajo** e **identidad alta**.

## 2. Comparación de marcos de lectura

| Marco | Mejor hit | E-value | % identidad | ¿Correcto? |
|-------|-----------|---------|-------------|------------|
| Forward_Frame_1 | Sin hits | — | — | No |
| Forward_Frame_2 | Sin hits | — | — | No |
| **Forward_Frame_3** | Insulina humana (P01308) | **3.36×10⁻⁷⁶** | **100%** | **Sí** |
| Reverse_Frame_1 | Sin hits | — | — | No |
| Reverse_Frame_2 | Sin hits | — | — | No |
| Reverse_Frame_3 | Sin hits | — | — | No |

**Conclusión:** solo `Forward_Frame_3` produce una secuencia que alinea perfectamente con la insulina anotada en Swiss-Prot. Coincide con el CDS de GenBank (`frame_annotation.txt`). Los otros 5 marcos no tienen homología significativa.

## 3. Top 5 hits del marco correcto (Forward_Frame_3)

1. **sp\|P01308.1 – Insulina humana (Homo sapiens)**
   - E-value: 3.36×10⁻⁷⁶ | Bit score: 572 | Identidad: 110/110 (100%)
   - Es la misma proteína que la anotada para el gen INS; valida la traducción.

2. **sp\|Q8HXV2.1 – Insulina (Gorilla gorilla)**
   - E-value: 2.07×10⁻⁷⁵ | Bit score: 567 | Identidad: 109/110 (99.1%)
   - Homólogo muy cercano; esperable por parentesco evolutivo con humanos.

3. **sp\|P30410.1 – Insulina (Pan troglodytes, chimpancé)**
   - E-value: 3.11×10⁻⁷⁵ | Bit score: 566 | Identidad: 108/110 (98.2%)

4. **sp\|P30406.1 – Insulina (Macaca fascicularis, mono)**
   - E-value: 8.00×10⁻⁷⁴ | Bit score: 557 | Identidad: 108/110 (98.2%)

5. **sp\|P30407.1 – Insulina (Chlorocebus aethiops, mono verde)**
   - E-value: 3.40×10⁻⁷³ | Bit score: 553 | Identidad: 107/110 (97.3%)

Todos los mejores hits son **preproinsulinas de primates**, lo cual confirma que la query es biológicamente correcta.

## 4. Significado de los valores estadísticos

| Parámetro | Interpretación en nuestros resultados |
|-----------|--------------------------------------|
| **E-value** | ~10⁻⁷⁶ → prácticamente imposible que el match sea por azar; homología real. |
| **Bit score** | 553–572 → alineamientos de muy alta calidad. |
| **% identidad** | 97–100% → secuencia muy conservada entre especies. |
| **Longitud** | 107–110 aa alineados → cubre casi toda la preproinsulina. |

## 5. Conclusión biológica

El BLAST confirma que **Forward_Frame_3** es el marco de lectura correcto del mRNA NM_000207. Los hits son insulinas de mamíferos/primates con identidad casi total, coherente con una proteína **funcionalmente crítica** y **evolutivamente conservada**. Esto vincula el gen **INS** con su producto proteico y fundamenta el análisis de alineamiento múltiple del Ejercicio 3.
