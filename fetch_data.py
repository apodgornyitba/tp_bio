from Bio import Entrez
import sys

Entrez.email = "test@example.com"
accession = "NM_000207"

try:
    print(f"Fetching {accession} from NCBI...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    data = handle.read()
    with open(f"{accession}.gbk", "w") as f:
        f.write(data)
    print(f"Saved {accession}.gbk successfully.")
except Exception as e:
    print(f"Error fetching data: {e}")
    sys.exit(1)
