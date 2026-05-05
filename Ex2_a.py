import sys
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 Ex2_a.py <input.fasta> <output.xml>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_xml = sys.argv[2]

    # Read the query sequence
    try:
        record = SeqIO.read(input_fasta, format="fasta")
        query_seq = record.seq
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        sys.exit(1)

    print(f"Running BLAST for {record.id} against swissprot database. This might take a few minutes...")
    try:
        # Run BLAST remotely
        result_handle = NCBIWWW.qblast("blastp", "swissprot", query_seq)
        
        # Save the raw XML result
        with open(output_xml, "w") as out_handle:
            out_handle.write(result_handle.read())
        print(f"BLAST results successfully saved to {output_xml}")

        # Parse the results for Ejercicio 2.b interpretation
        result_handle.seek(0)
        blast_record = NCBIXML.read(result_handle)
        
        print("\n--- BLAST Interpretation Summary (Ejercicio 2.b) ---")
        print("Interpreting the top 5 alignments found:")
        for alignment in blast_record.alignments[:5]:
            for hsp in alignment.hsps:
                print(f"\nSequence: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"E-value: {hsp.expect}")
                print(f"Identities: {hsp.identities}/{hsp.align_length} ({(hsp.identities/hsp.align_length)*100:.2f}%)")
                print("Explanation: The E-value indicates the number of expected hits of similar quality by chance.")
                print("A very low E-value (close to 0) means the match is statistically highly significant.")
                print("High identity % means the sequence is very conserved compared to the query.")
                break # Only show the top HSP for each alignment

    except Exception as e:
        print(f"Error running BLAST: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
