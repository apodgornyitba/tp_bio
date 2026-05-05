# TP Cuatrimestral - Parte 1 - Bioinformática

Este repositorio contiene los scripts requeridos para el TP Cuatrimestral de Bioinformática. Se eligió investigar sobre el gen humano **INS (Insulina)**.

## Requisitos previos

El entorno de ejecución requiere tener instalados `python3`, `biopython` y el software `muscle` para los alineamientos múltiples (MSA).
El script de instalación de dependencias corre automáticamente:
```bash
sudo apt-get update && sudo apt-get install -y python3-biopython muscle
```

## Archivos y Scripts

1. **`fetch_data.py`**: Descarga automáticamente la secuencia madura (mRNA) en formato GenBank del gen INS (`NM_000207`) desde NCBI utilizando Entrez.
2. **`Ex1.py`**: *(Ejercicio 1)* Recibe un archivo `.gbk`, busca la secuencia correcta del Open Reading Frame (ORF) según la anotación `CDS` del GenBank (o busca el ORF más largo si no la hay) y la traduce a aminoácidos, exportando el resultado a formato FASTA (`NM_000207.fasta`).
3. **`Ex2_a.py`**: *(Ejercicio 2.a y 2.b)* Lee la secuencia en FASTA generada en el paso anterior y realiza una búsqueda remota de BLAST contra la base de datos `swissprot`. Guarda los resultados en `blast_results.xml` e imprime una interpretación de los mejores alineamientos obtenidos (Explicando el E-value y % de identidad).
4. **`Ex3.py`**: *(Ejercicio 3)* Analiza el archivo de resultados XML de BLAST, extrae los ID de acceso de los 10 mejores resultados (homólogos a la insulina en otras especies) y obtiene sus secuencias completas desde NCBI Entrez. Combina estos 10 secuencias con la secuencia original y realiza un Multiple Sequence Alignment (MSA) local utilizando MUSCLE, exportando los resultados a `msa_output.afa`.
5. **`run_pipeline.sh`**: Es el script de Bash principal que orquesta la ejecución secuencial de todos los pasos, loguea los resultados, realiza chequeo de errores y valida la correcta generación de los archivos.

## Cómo ejecutar

Para correr todos los ejercicios secuencialmente y generar los outputs y logs solicitados:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh
```

El script principal dejará registro de todo el proceso en `pipeline.log` y generará los siguientes archivos:
- `NM_000207.gbk` (Input Ej 1)
- `NM_000207.fasta` (Output Ej 1 / Input Ej 2)
- `blast_results.xml` (Output Ej 2 / Input Ej 3)
- `msa_input.fasta` (Input secundario Ej 3)
- `msa_output.afa` (Output Ej 3 - Alineamiento Múltiple)
