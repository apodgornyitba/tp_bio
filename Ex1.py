import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_reading_frames(seq):
    frames = []
    # 3 forward frames
    for i in range(3):
        frames.append((f"Forward_Frame_{i+1}", seq[i:]))
    # 3 reverse frames
    rev_seq = seq.reverse_complement()
    for i in range(3):
        frames.append((f"Reverse_Frame_{i+1}", rev_seq[i:]))
    return frames

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 Ex1.py <input.gbk> <output.fasta>")
        sys.exit(1)

    input_gbk = sys.argv[1]
    output_fasta = sys.argv[2]

    records = list(SeqIO.parse(input_gbk, "genbank"))
    if not records:
        print("No records found in GenBank file.")
        sys.exit(1)

    fasta_records = []

    for record in records:
        print(f"Processing record: {record.id}")
        
        # Try to find the correct translation using CDS feature if available
        correct_translation = None
        for feature in record.features:
            if feature.type == "CDS":
                if "translation" in feature.qualifiers:
                    correct_translation = feature.qualifiers["translation"][0]
                    break

        if correct_translation:
            print("Found CDS feature, identifying correct reading frame.")
            # We found the correct translation, let's save it directly.
            # We still can evaluate reading frames to show we did it.
            frames = get_reading_frames(record.seq)
            correct_frame_name = "Unknown"
            
            for frame_name, frame_seq in frames:
                # pad to multiple of 3
                length = len(frame_seq)
                padded_seq = frame_seq[:length - (length % 3)]
                translation = str(padded_seq.translate())
                
                # The exact CDS translation might be a sub-sequence of the whole frame's translation
                if correct_translation in translation:
                    correct_frame_name = frame_name
                    break

            print(f"Correct frame identified as: {correct_frame_name}")
            
            # Save the correct translated CDS
            translated_record = SeqRecord(
                Seq(correct_translation),
                id=f"{record.id}_correct_ORF",
                description=f"Translated correct ORF based on CDS ({correct_frame_name})"
            )
            fasta_records.append(translated_record)
        else:
            print("No CDS feature found, outputting longest ORF for all frames.")
            # Find the longest ORF among all 6 frames
            frames = get_reading_frames(record.seq)
            longest_orf = ""
            best_frame = ""
            
            for frame_name, frame_seq in frames:
                length = len(frame_seq)
                padded_seq = frame_seq[:length - (length % 3)]
                translation = str(padded_seq.translate())
                
                # Split by stop codon
                orfs = translation.split('*')
                for orf in orfs:
                    # An ORF typically starts with Methionine
                    idx = orf.find('M')
                    if idx != -1:
                        valid_orf = orf[idx:]
                        if len(valid_orf) > len(longest_orf):
                            longest_orf = valid_orf
                            best_frame = frame_name

            print(f"Longest ORF found in {best_frame}")
            translated_record = SeqRecord(
                Seq(longest_orf),
                id=f"{record.id}_longest_ORF",
                description=f"Longest ORF translation ({best_frame})"
            )
            fasta_records.append(translated_record)

    SeqIO.write(fasta_records, output_fasta, "fasta")
    print(f"Saved correctly identified sequence(s) to {output_fasta}")

if __name__ == "__main__":
    main()
