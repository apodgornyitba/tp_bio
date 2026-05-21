import sys
import subprocess
import shutil
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_muscle(msa_input, msa_output):
    # Support both MUSCLE v3 (-in/-out) and newer syntax (-align/-output).
    if shutil.which("muscle") is None:
        raise RuntimeError("MUSCLE no está instalado o no está en PATH")

    commands = [
        ["muscle", "-align", msa_input, "-output", msa_output],
        ["muscle", "-in", msa_input, "-out", msa_output],
    ]

    last_error = None
    for cmd in commands:
        try:
            subprocess.run(cmd, check=True)
            return
        except subprocess.CalledProcessError as e:
            last_error = e

    raise RuntimeError(f"No se pudo ejecutar MUSCLE con una sintaxis compatible: {last_error}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 Ex3.py <blast_results.xml> <query.fasta>")
        sys.exit(1)

    blast_xml = sys.argv[1]
    query_fasta = sys.argv[2]
    
    Entrez.email = "test@example.com"
    
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
        if len(title_parts) >= 2 and (title_parts[0].startswith('sp') or title_parts[0].startswith('tr')):
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

    # Run Muscle
    msa_output = "msa_output.afa"
    print("Running MUSCLE for Multiple Sequence Alignment...")
    try:
        run_muscle(msa_input, msa_output)
        print(f"MSA successfully saved to {msa_output}")
    except Exception as e:
        print(f"Error running MUSCLE: {e}")
        sys.exit(1)
        
    print("\n--- MSA Interpretation Summary (Ejercicio 3) ---")
    print("El alineamiento múltiple se realizó con la secuencia de consulta y los 10 mejores resultados de BLAST.")
    print("A partir del archivo generado (.afa), podemos observar las regiones conservadas a lo largo de la evolución de la insulina en diferentes especies.")
    print("Las columnas sin huecos (gaps) y con los mismos aminoácidos indican regiones funcional y estructuralmente críticas para la proteína.")

if __name__ == "__main__":
    main()
