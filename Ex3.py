#!/usr/bin/env python3
import os
import subprocess
import sys
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from platform_tools import find_msa_tool, print_install_hints


def run_muscle(muscle_bin, msa_input, msa_output):
    try:
        subprocess.run(
            [muscle_bin, "-align", msa_input, "-output", msa_output],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError:
        subprocess.run(
            [muscle_bin, "-in", msa_input, "-out", msa_output],
            check=True,
            capture_output=True,
            text=True,
        )


def run_mafft(mafft_bin, msa_input, msa_output):
    with open(msa_output, "w", encoding="utf-8") as out:
        subprocess.run(
            [mafft_bin, "--auto", msa_input],
            check=True,
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
        )


def main():
    if len(sys.argv) < 3:
        print("Usage: python Ex3.py <blast_results.xml> <query.fasta>")
        sys.exit(1)

    blast_xml = sys.argv[1]
    query_fasta = sys.argv[2]
    Entrez.email = os.environ.get("ENTREZ_EMAIL", "estudiante@itba.edu.ar")

    try:
        query_record = SeqIO.read(query_fasta, "fasta")
    except Exception as e:
        print(f"Error reading query FASTA: {e}")
        sys.exit(1)

    print("Parsing BLAST XML results...")
    try:
        with open(blast_xml, encoding="utf-8") as f:
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

        title_parts = alignment.title.split("|")
        accession = ""
        if len(title_parts) >= 2 and (
            title_parts[0].startswith("sp") or title_parts[0].startswith("tr")
        ):
            accession = title_parts[1].split(".")[0]
        else:
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
            hsp = alignment.hsps[0]
            seq_str = hsp.sbjct.replace("-", "")
            record = SeqRecord(Seq(seq_str), id=accession, description=alignment.title)
            top10_records.append(record)
            count += 1

    msa_input = "msa_input.fasta"
    SeqIO.write(top10_records, msa_input, "fasta")
    print(f"Saved {len(top10_records)} sequences to {msa_input}")

    msa_output = "msa_output.afa"
    tool, tool_bin = find_msa_tool()
    if not tool:
        print("ERROR: No hay herramienta MSA (MUSCLE ni MAFFT).")
        print_install_hints("msa")
        sys.exit(1)

    print(f"Running {tool.upper()} ({tool_bin})...")
    try:
        if tool == "muscle":
            run_muscle(tool_bin, msa_input, msa_output)
        else:
            run_mafft(tool_bin, msa_input, msa_output)
        print(f"MSA saved to {msa_output}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {tool.upper()}: {e}")
        if e.stderr:
            print(e.stderr)
        sys.exit(1)

    print("\n--- MSA Interpretation Summary (Ejercicio 3) ---")
    print("Alineamiento de query + 10 mejores hits BLAST.")


if __name__ == "__main__":
    main()
