#!/usr/bin/env python3
"""
Ejercicio 2.a - BLAST remoto (NCBI qblast contra Swiss-Prot).
"""
import argparse
import sys
from pathlib import Path

from Bio.Blast import NCBIWWW
from Bio import SeqIO

from blast_common import (
    finalize_best_result,
    parse_blast_xml,
    print_interpretation,
    top_hit_stats,
)


def run_remote_blast(record, output_dir: Path):
    safe_id = record.id.replace("/", "_")
    output_xml = output_dir / f"blast_{safe_id}.xml"

    print(f"\nBLAST REMOTO para {record.id} ({len(record.seq)} aa)...")
    print("Puede tardar varios minutos por secuencia.")

    result_handle = NCBIWWW.qblast("blastp", "swissprot", record.seq)
    xml_data = result_handle.read()
    result_handle.close()

    output_xml.write_text(xml_data, encoding="utf-8")
    blast_record = parse_blast_xml(output_xml)
    stats = top_hit_stats(blast_record)

    if stats:
        print(f"  Mejor hit: E-value={stats['evalue']:.2e}, id={stats['identity_pct']:.1f}%")
    else:
        print("  Sin hits significativos.")
    return output_xml, stats


def main():
    parser = argparse.ArgumentParser(description="BLASTp remoto por marco de lectura")
    parser.add_argument("input_fasta")
    parser.add_argument("output_dir", nargs="?", default="blast_results")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    if not records:
        print("No sequences found in FASTA.")
        sys.exit(1)

    print(f"BLAST REMOTO: {len(records)} secuencia(s) contra swissprot (NCBI).")
    results = []

    for record in records:
        try:
            xml_path, stats = run_remote_blast(record, output_dir)
            results.append({
                "record_id": record.id,
                "record": record,
                "xml_path": str(xml_path),
                "stats": stats,
            })
        except Exception as e:
            print(f"Error en BLAST remoto de {record.id}: {e}")
            sys.exit(1)

    summary_path = output_dir / "blast_summary.txt"
    best = finalize_best_result(results, summary_path)
    print_interpretation(best)
    print(f"\nResumen: {summary_path}")
    print("Salida principal: blast_results.xml, query_best.fasta")


if __name__ == "__main__":
    main()
