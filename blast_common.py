"""Funciones compartidas para BLAST remoto y local."""
from pathlib import Path

from Bio import SeqIO
from Bio.Blast import NCBIXML


def top_hit_stats(blast_record):
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


def parse_blast_xml(xml_path):
    with open(xml_path) as handle:
        return NCBIXML.read(handle)


def write_summary(summary_path, results, title_suffix=""):
    lines = [
        f"# Resumen BLAST por marco de lectura{title_suffix}",
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
    return best


def finalize_best_result(results, summary_path, xml_copy="blast_results.xml",
                         fasta_copy="query_best.fasta"):
    best = write_summary(summary_path, results)
    Path(xml_copy).write_text(Path(best["xml_path"]).read_text(encoding="utf-8"), encoding="utf-8")
    SeqIO.write(best["record"], fasta_copy, "fasta")
    return best


def print_interpretation(best):
    print("\n--- BLAST Interpretation Summary (Ejercicio 2.b) ---")
    print(f"Marco seleccionado para Ejercicio 3: {best['record_id']}")
    if best["stats"]:
        s = best["stats"]
        print(f"Mejor hit: {s['title']}")
        print(f"E-value: {s['evalue']:.2e}")
        print(f"Bit score: {s['score']}")
        print(f"Identidad: {s['identity_pct']:.2f}%")
