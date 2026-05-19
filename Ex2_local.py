#!/usr/bin/env python3
"""
Ejercicio 2.a - BLAST local (BLAST+ contra Swiss-Prot).
Requiere: blastp, makeblastdb y base local (prepare_blast_db.py).
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

from blast_common import finalize_best_result, parse_blast_xml, print_interpretation, top_hit_stats
from platform_tools import default_blast_db, find_blastp, print_install_hints


def run_local_blast(record, blast_db: Path, output_xml: Path, blastp_bin: str):
    print(f"\nBLAST LOCAL para {record.id} ({len(record.seq)} aa)...")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        SeqIO.write(record, tmp, "fasta")
        query_path = tmp.name

    db_str = str(blast_db)
    # Quitar extension si el usuario paso swissprot_db.phr etc.
    for ext in (".phr", ".pin", ".psq", ".pog", ".psd", ".psi"):
        if db_str.endswith(ext):
            db_str = db_str[: -len(ext)]
            break

    cmd = [
        blastp_bin,
        "-query", query_path,
        "-db", db_str,
        "-out", str(output_xml),
        "-outfmt", "5",
        "-evalue", "10",
        "-max_target_seqs", "20",
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    Path(query_path).unlink(missing_ok=True)

    blast_record = parse_blast_xml(output_xml)
    stats = top_hit_stats(blast_record)
    if stats:
        print(f"  Mejor hit: E-value={stats['evalue']:.2e}, id={stats['identity_pct']:.1f}%")
    else:
        print("  Sin hits.")
    return output_xml, stats


def main():
    parser = argparse.ArgumentParser(description="BLASTp local por marco de lectura")
    parser.add_argument("input_fasta")
    parser.add_argument("output_dir", nargs="?", default="blast_results_local")
    parser.add_argument("--db", default=None, help="Prefijo base BLAST (env BLAST_DB)")
    args = parser.parse_args()

    blastp = find_blastp()
    if not blastp:
        print("ERROR: blastp no encontrado.")
        print_install_hints("blast")
        sys.exit(1)

    blast_db = Path(args.db) if args.db else default_blast_db()
    if not Path(f"{blast_db}.phr").exists():
        print(f"ERROR: Base BLAST no encontrada: {blast_db}.phr")
        print("Ejecutar primero: python prepare_blast_db.py")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    if not records:
        print("No hay secuencias en el FASTA.")
        sys.exit(1)

    print(f"BLAST LOCAL: {len(records)} secuencia(s), DB={blast_db}")
    results = []
    for record in records:
        safe_id = record.id.replace("/", "_")
        xml_path = output_dir / f"blast_{safe_id}.xml"
        try:
            xml_path, stats = run_local_blast(record, blast_db, xml_path, blastp)
            results.append({
                "record_id": record.id,
                "record": record,
                "xml_path": str(xml_path),
                "stats": stats,
            })
        except subprocess.CalledProcessError as e:
            print(f"Error BLAST local {record.id}: {e.stderr or e}")
            sys.exit(1)

    summary_path = output_dir / "blast_summary.txt"
    best = finalize_best_result(
        results,
        summary_path,
        xml_copy="blast_results_local.xml",
        fasta_copy="query_best_local.fasta",
    )
    print_interpretation(best)
    print(f"\nResumen: {summary_path}")
    print("Salida principal: blast_results_local.xml")


if __name__ == "__main__":
    main()
