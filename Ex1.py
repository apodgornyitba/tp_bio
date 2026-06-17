import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def get_reading_frames(seq):
    """Devuelve los 6 marcos de lectura (3 forward + 3 reverse)."""
    frames = []
    for i in range(3):
        frames.append((f"Forward_Frame_{i + 1}", seq[i:]))
    rev_seq = seq.reverse_complement()
    for i in range(3):
        frames.append((f"Reverse_Frame_{i + 1}", rev_seq[i:]))
    return frames


def translate_frame(frame_seq):
    """Traduce un marco recortando a múltiplo de 3."""
    length = len(frame_seq)
    padded = frame_seq[: length - (length % 3)]
    return str(padded.translate())


def find_cds_translation(record):
    """Obtiene la traducción anotada en la feature CDS, si existe."""
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            return feature.qualifiers["translation"][0]
    return None


def identify_frame_with_cds(frames, cds_translation):
    """Identifica en qué marco aparece la traducción del CDS."""
    for frame_name, frame_seq in frames:
        translation = translate_frame(frame_seq)
        if cds_translation in translation:
            return frame_name
    return None


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 Ex1.py <input.gbk> <output.fasta>")
        sys.exit(1)

    input_gbk = sys.argv[1]
    output_fasta = sys.argv[2]
    output_fasta_path = Path(output_fasta)
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(input_gbk, "genbank"))
    if not records:
        print("No records found in GenBank file.")
        sys.exit(1)

    fasta_records = []
    annotation_lines = []

    for record in records:
        print(f"Processing record: {record.id}")
        frames = get_reading_frames(record.seq)
        cds_translation = find_cds_translation(record)
        cds_frame = None

        if cds_translation:
            cds_frame = identify_frame_with_cds(frames, cds_translation)
            print(f"CDS encontrado. Marco de lectura esperado: {cds_frame}")
            annotation_lines.append(f"record={record.id}")
            annotation_lines.append(f"cds_frame={cds_frame}")
            annotation_lines.append(f"cds_length={len(cds_translation)}")

        for frame_name, frame_seq in frames:
            translation = translate_frame(frame_seq)
            seq_id = f"{record.id}_{frame_name}"
            description = f"Traduccion completa del marco {frame_name}"
            if cds_frame == frame_name:
                description += " [coincide con CDS de GenBank]"

            fasta_records.append(
                SeqRecord(Seq(translation), id=seq_id, description=description)
            )
            print(f"  {frame_name}: {len(translation)} aa (stops: {translation.count('*')})")

    SeqIO.write(fasta_records, str(output_fasta_path), "fasta")
    print(f"\nGuardadas {len(fasta_records)} secuencias en {output_fasta}")

    annotation_file = output_fasta_path.with_name("frame_annotation.txt")
    with open(annotation_file, "w") as f:
        f.write("\n".join(annotation_lines) + "\n")
        f.write(f"total_frames={len(fasta_records)}\n")
    print(f"Anotacion de marcos guardada en {annotation_file}")


if __name__ == "__main__":
    main()
