import os
import sys
from pathlib import Path
from Bio import Entrez

Entrez.email = os.environ.get("ENTREZ_EMAIL", "estudiante@itba.edu.ar")
accession = "NM_000207"

try:
    output_path = sys.argv[1] if len(sys.argv) > 1 else f"{accession}.gbk"
    output_parent = Path(output_path).parent
    if str(output_parent) not in ("", "."):
        output_parent.mkdir(parents=True, exist_ok=True)
    print(f"Fetching {accession} from NCBI (email: {Entrez.email})...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    data = handle.read()
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(data)
    print(f"Saved {output_path} successfully.")
except Exception as e:
    print(f"Error fetching data: {e}")
    sys.exit(1)
