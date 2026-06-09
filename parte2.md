Bioinformatica Trabajo Practico Grupal - Cuatrimestral
pág. 1
Introducción a la Bioinformática
Trabajo Práctico – Parte 2
Ejercicio 4 - EMBOSS. Instalar EMBOSS. Escribir un script que llame a algún programa EMBOSS para
que realice algún análisis sobre la una secuencia de nucleótidos fasta (del Ej. 1). Por ejemplo que
calcule los ORF y obtenga las secuencias de proteínas posibles. Luego bájense los motivos de las bases
de datos PROSITE (archivo prosite.dat) y por medio del llamado a otro programa EMBOSS realizar el
análisis de dominios de las secuencias de aminoácidos obtenidas y escribir los resultados en un archivo
de salida.
− Input : Archivo de secuencias Fasta (ej. Xxxxx.fas con una o más secuencias de aa.
− Output: Archivo de resultados del dominios encontrados en las secuencias de aa.
Ejercicio 5. A partir del transcipto seleccionado para el análisis de su investigación (aquel que presenta
la/las mutaciones que dan lugar a la fenotipia patológica investigada por ustedes), se requiere poder
crear 5 primers que permitan posteriormente realizar análisis cuanti y cualitativos de la presencia de
dicha variante del mensajero en pacientes. Diseñarlos de forma que su tamaño sea de 18 a 24 pares
de bases, posean mínimo un 50% de GC y un máximo de 60%, debemos evitar tener GC en extremos
terminales y su temperatura de melting debe ser igual o inferior a 67 grados celcius.
*Crear un script que lleve a cabo dicho diseño de forma que sea parametizable, reciba como parámetro
de entrada la secuencia del transcipto, y utilice los parámetros de diseño a partir de un archivo de
configuración en formato XML o JSON (o aquel que les resulte mas sencillo de implementar). De no
poder realizarlo a través de Python o emboss o alguna otra distribución o librería de programacion,
pueden utilizar herramientas web y evidenciar los resultados en el informe entregable utilizando
screenshots del proceso realizado
Ejercicio 6. Trabajo con Bases de Datos Biológicas (útil para la presentación de la investigación Ej. 6).
a) A partir del gen o proteína de interés para ustedes dar su link a NCBI-Gene como una entrada de
Entrez, por ej.: http://www.ncbi.nlm.nih.gov/gene/3630 Expliquen brevemente lo que hace la proteína
y por qué la eligieron.
b) ¿Cuántos genes / proteínas homólogas se conocen en otros organismos? Utilicen la información que
está en la base de datos de HomoloGene y en la bases de datos Ensembl . Describan los resultados en
ambas bases de datos, y en qué se diferencian. Mencionen sobre qué tan común creen son estos genes
o proteínas y a qué grupos taxonómicos pertenecen (sólo en las bacterias, en los vertebrados, etc.)
Bioinformatica Trabajo Practico Grupal - Cuatrimestral
pág. 2
c) ¿Cuántos transcriptos y cuántas formas alternativas de splicing son conocidos para este gen /
proteína? ¿Cuáles de estos splicing alternativos se expresan? ¿Tienen funciones alternativas? Buscar
evidencia de esto en las base de datos de NCBI y en los transcriptos de Ensembl ¿Cómo el número de
splicings alternativos diferente entre las dos bases de datos y cuál piensan que es más precisa y por
qué?
d) ¿Con cuántas otras proteínas interactúa el producto génico de su gen? ¿Existe un patrón o relación
entre las interacciones? Mencione las interacciones interesantes o inusuales. Usted encontrará las
interacciones de su gene/proteína tanto en la base de datos NCBI Gene como en la base de datos
UniProt . Compare las dos tablas entre sí. ¿Hay proteínas que interactúan únicas para cada tabla?
e) Expliquen brevemente de qué componente celular forma parte su proteína (pista: se puede estudiar
la información de Gene Ontology - GO), ¿A qué procesos biológicos pertenece (pista idem)? y ¿En qué
función molecular trabaja esta proteína? Los términos ontológicos de genes los pueden encontrar
tanto en NCBI Gene y en la base de datos UniProt como haciendo una búsqueda en AmiGO. 2
f) Discutan brevemente en qué estructura o vías metabólicas específicas (pathways) estaría
participando su gen / proteína? (Reactome, KEGG son algunas bases de datos de pathways).
g) Entrar en la base de datos de variantes genéticas dbSNP e intentar interpretar o encontrar info sobre
alguna variante (reference SNP - rsXXXX) asociada con la patología investigada en su gen de interés.
¿Qué variante es? ¿Hay información sobre la frecuencia que tiene esta variante en la población? ¿Qué
grupo étnico parece ser el más afectado?
NOTA: Para hacer este ejercicio les pueden servir algunas otras bases de datos como:
http://www.genecards.org
http://www.ncbi.nlm.nih.gov/clinvar/ (para obtener información clínica del gen y sus variantes)
https://ghr.nlm.nih.gov
Ejercicio 7. Armar una presentación donde expliquen la enfermedad que investigaron, lo que hicieron
y los resultados que fueron obteniendo en los 5 ejercicios del TP.
Los integrantes de cada grupo tendrán un máximo de 10 minutos para exponer como realizaron el
trabajo práctico y comentar sobre su investigación (y no tanto sobre el código implementado). La
correcta exposición del trabajo realizado por los miembros del grupo también entra en la evaluación