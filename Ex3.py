import os
import shutil
import sys
import subprocess
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def find_executable(*names, extra_paths=()):
    """Busca un ejecutable en PATH y rutas extra."""
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    for candidate in extra_paths:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def run_muscle(muscle_bin, msa_input, msa_output):
    """Ejecuta MUSCLE probando sintaxis v5 y v3."""
    try:
        subprocess.run(
            [muscle_bin, "-align", msa_input, "-output", msa_output],
            check=True,
            capture_output=True,
            text=True,
        )
        return "muscle"
    except subprocess.CalledProcessError:
        pass

    subprocess.run(
        [muscle_bin, "-in", msa_input, "-out", msa_output],
        check=True,
        capture_output=True,
        text=True,
    )
    return "muscle"


def run_mafft(mafft_bin, msa_input, msa_output):
    """Ejecuta MAFFT como alternativa a MUSCLE."""
    with open(msa_output, "w") as out:
        subprocess.run(
            [mafft_bin, "--auto", msa_input],
            check=True,
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
        )
    return "mafft"


def find_msa_tool():
    """Devuelve (herramienta, ruta_binario) o (None, None)."""
    brew_paths = ("/opt/homebrew/bin", "/usr/local/bin")

    muscle_bin = find_executable(
        "muscle", "muscle5",
        extra_paths=[os.path.join(p, "muscle") for p in brew_paths],
    )
    if muscle_bin:
        return "muscle", muscle_bin

    mafft_bin = find_executable(
        "mafft",
        extra_paths=[os.path.join(p, "mafft") for p in brew_paths],
    )
    if mafft_bin:
        return "mafft", mafft_bin

    return None, None

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 Ex3.py <blast_results.xml> <query.fasta>")
        sys.exit(1)

    blast_xml = sys.argv[1]
    query_fasta = sys.argv[2]

    Entrez.email = os.environ.get("ENTREZ_EMAIL", "estudiante@itba.edu.ar")
    
    # Read query sequence
    try:
        query_record = SeqIO.read(query_fasta, "fasta")
    except Exception as e:
        print(f"Error reading query FASTA: {e}")
        sys.exit(1)

    print("Parsing BLAST XML results...")
    try:
        with open(blast_xml) as f:
            blast_record = NCBIXML.read(f)
    except Exception as e:
        print(f"Error reading BLAST XML: {e}")
        sys.exit(1)

    top10_records = [query_record]
    seen_accessions = set()
    
    print("Fetching top 10 hits from NCBI...")
    count = 0
    for alignment in blast_record.alignments:
        if count >= 10:
            break
        
        # alignment.title is usually something like 'sp|P01308.1| RecName: Full=Insulin...'
        # Let's extract the accession ID
        title_parts = alignment.title.split('|')
        accession = ""
        if len(title_parts) >= 2 and (
            title_parts[0].startswith("sp") or title_parts[0].startswith("tr")
        ):
            accession = title_parts[1].split('.')[0]
        else:
            # If format is different, try to use the first word or alignment.accession
            accession = alignment.accession
            
        if accession in seen_accessions:
            continue
        
        seen_accessions.add(accession)
        print(f"Fetching {accession}...")
        try:
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            top10_records.append(record)
            count += 1
        except Exception as e:
            print(f"Failed to fetch {accession}: {e}")
            # Fallback to hsp.sbjct if Entrez fails
            hsp = alignment.hsps[0]
            seq_str = hsp.sbjct.replace('-', '') # Remove gaps
            record = SeqRecord(Seq(seq_str), id=accession, description=alignment.title)
            top10_records.append(record)
            count += 1
            
    # Save sequences to a fasta file for MSA
    msa_input = "msa_input.fasta"
    SeqIO.write(top10_records, msa_input, "fasta")
    print(f"Saved {len(top10_records)} sequences to {msa_input}")

    msa_output = "msa_output.afa"
    tool, tool_bin = find_msa_tool()
    if not tool:
        print("ERROR: No hay herramienta de MSA instalada (MUSCLE ni MAFFT).")
        print("  macOS:   brew install mafft")
        print("  Ubuntu:  sudo apt-get install muscle")
        print("Luego volver a ejecutar: python3 Ex3.py blast_results.xml query_best.fasta")
        sys.exit(1)

    print(f"Running {tool.upper()} ({tool_bin}) for Multiple Sequence Alignment...")
    try:
        if tool == "muscle":
            run_muscle(tool_bin, msa_input, msa_output)
        else:
            run_mafft(tool_bin, msa_input, msa_output)
        print(f"MSA successfully saved to {msa_output}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {tool.upper()}: {e}")
        if e.stderr:
            print(e.stderr)
        sys.exit(1)
        
    print("\n--- MSA Interpretation Summary (Ejercicio 3) ---")
    print("El alineamiento múltiple se realizó con la secuencia de consulta y los 10 mejores resultados de BLAST.")
    print("A partir del archivo generado (.afa), podemos observar las regiones conservadas a lo largo de la evolución de la insulina en diferentes especies.")
    print("Las columnas sin huecos (gaps) y con los mismos aminoácidos indican regiones funcional y estructuralmente críticas para la proteína.")

if __name__ == "__main__":
    main()
