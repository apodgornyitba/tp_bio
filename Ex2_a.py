import sys
from pathlib import Path

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


def top_hit_stats(blast_record):
    """Extrae estadisticas del mejor hit de un resultado BLAST."""
    if not blast_record.alignments:
        return None

    alignment = blast_record.alignments[0]
    if not alignment.hsps:
        return None

    hsp = alignment.hsps[0]
    identity_pct = (hsp.identities / hsp.align_length) * 100 if hsp.align_length else 0
    return {
        "title": alignment.title,
        "evalue": hsp.expect,
        "score": hsp.score,
        "identities": hsp.identities,
        "align_length": hsp.align_length,
        "identity_pct": identity_pct,
    }


def run_blast_for_record(record, output_dir):
    """Ejecuta BLASTp remoto para una secuencia y guarda el XML."""
    safe_id = record.id.replace("/", "_")
    output_xml = output_dir / f"blast_{safe_id}.xml"

    print(f"\nBLAST para {record.id} (longitud {len(record.seq)} aa)...")
    print("Esto puede tardar varios minutos por secuencia.")

    result_handle = NCBIWWW.qblast("blastp", "swissprot", record.seq)
    xml_data = result_handle.read()

    with open(output_xml, "w") as out_handle:
        out_handle.write(xml_data)

    result_handle.close()
    with open(output_xml) as in_handle:
        blast_record = NCBIXML.read(in_handle)

    stats = top_hit_stats(blast_record)
    if stats:
        print(
            f"  Mejor hit: E-value={stats['evalue']:.2e}, "
            f"identidad={stats['identity_pct']:.1f}%"
        )
        print(f"  {stats['title'][:80]}...")
    else:
        print("  Sin hits significativos.")

    return output_xml, stats


def write_summary(summary_path, results):
    """Escribe resumen comparativo de los 6 marcos."""
    lines = [
        "# Resumen BLAST por marco de lectura (Ejercicio 2)",
        "# Menor E-value = mejor candidato al marco correcto",
        "",
    ]

    sorted_results = sorted(
        results,
        key=lambda r: r["stats"]["evalue"] if r["stats"] else float("inf"),
    )

    for rank, item in enumerate(sorted_results, start=1):
        stats = item["stats"]
        lines.append(f"## {rank}. {item['record_id']}")
        if stats:
            lines.append(f"- Archivo: {item['xml_path']}")
            lines.append(f"- Mejor hit: {stats['title']}")
            lines.append(f"- E-value: {stats['evalue']:.2e}")
            lines.append(f"- Bit score: {stats['score']}")
            lines.append(
                f"- Identidad: {stats['identities']}/{stats['align_length']} "
                f"({stats['identity_pct']:.2f}%)"
            )
        else:
            lines.append("- Sin hits")
        lines.append("")

    best = sorted_results[0]
    lines.append(f"MARCO_SELECCIONADO={best['record_id']}")
    lines.append(f"BLAST_XML={best['xml_path']}")

    summary_path.write_text("\n".join(lines), encoding="utf-8")


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 Ex2_a.py <input.fasta> <output_dir>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        records = list(SeqIO.parse(input_fasta, "fasta"))
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        sys.exit(1)

    if not records:
        print("No sequences found in FASTA.")
        sys.exit(1)

    print(f"Se encontraron {len(records)} secuencia(s). BLAST contra swissprot.")
    results = []

    for record in records:
        try:
            xml_path, stats = run_blast_for_record(record, output_dir)
            results.append(
                {
                    "record_id": record.id,
                    "record": record,
                    "xml_path": str(xml_path),
                    "stats": stats,
                }
            )
        except Exception as e:
            print(f"Error en BLAST de {record.id}: {e}")
            sys.exit(1)

    summary_path = output_dir / "blast_summary.txt"
    write_summary(summary_path, results)

    best = min(
        results,
        key=lambda r: r["stats"]["evalue"] if r["stats"] else float("inf"),
    )
    best_xml = Path("blast_results.xml")
    best_fasta = Path("query_best.fasta")

    best_xml.write_text(Path(best["xml_path"]).read_text(encoding="utf-8"), encoding="utf-8")
    SeqIO.write(best["record"], best_fasta, "fasta")

    print("\n--- BLAST Interpretation Summary (Ejercicio 2.b) ---")
    print(f"Marco seleccionado para Ejercicio 3: {best['record_id']}")
    if best["stats"]:
        s = best["stats"]
        print(f"Mejor hit: {s['title']}")
        print(f"E-value: {s['evalue']:.2e} (menor = mas significativo)")
        print(f"Bit score: {s['score']}")
        print(f"Identidad: {s['identity_pct']:.2f}%")
        print(
            "\nInterpretacion: un E-value muy bajo contra proteinas de insulina "
            "indica que este marco de lectura es el correcto."
        )

    print(f"\nResumen guardado en: {summary_path}")
    print(f"BLAST principal para Ej.3: {best_xml}")
    print(f"Query seleccionada: {best_fasta}")


if __name__ == "__main__":
    main()
