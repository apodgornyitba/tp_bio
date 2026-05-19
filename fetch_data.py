import os
import sys
from Bio import Entrez

Entrez.email = os.environ.get("ENTREZ_EMAIL", "estudiante@itba.edu.ar")
accession = "NM_000207"

try:
    print(f"Fetching {accession} from NCBI (email: {Entrez.email})...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    data = handle.read()
    with open(f"{accession}.gbk", "w") as f:
        f.write(data)
    print(f"Saved {accession}.gbk successfully.")
except Exception as e:
    print(f"Error fetching data: {e}")
    sys.exit(1)
