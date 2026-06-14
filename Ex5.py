import json
import os
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt


DEFAULT_CONFIG = {
    "min_length": 18,
    "max_length": 24,
    "min_gc": 50.0,
    "max_gc": 60.0,
    "avoid_gc_ends": True,
    "max_tm": 67.0,
    "target_tm": 60.0,
    "num_primers": 5,
    "num_primer_pairs": 5,
    "min_amplicon_size": 70,
    "max_amplicon_size": 180,
    "max_pair_tm_diff": 2.0,
}


def log(msg):
    print(f"[Ex5] {msg}")


def calculate_gc(seq_str):
    """Calcula el porcentaje de GC de una secuencia."""
    g_or_c = sum(1 for c in seq_str if c.upper() in ("G", "C"))
    return (g_or_c / len(seq_str)) * 100.0


def calculate_tm(seq_obj):
    """Calcula la temperatura de melting usando Tm_NN o un fallback de sal."""
    try:
        return mt.Tm_NN(seq_obj)
    except Exception:
        seq_str = str(seq_obj).upper()
        gc_count = seq_str.count("G") + seq_str.count("C")
        return 64.9 + 41.0 * (gc_count - 16.4) / len(seq_str)


def check_terminal_ends(seq_str):
    """Evita tener G o C en los extremos terminales (5' y 3' del primer)."""
    first_char = seq_str[0].upper()
    last_char = seq_str[-1].upper()
    return first_char in ("A", "T") and last_char in ("A", "T")


def load_config(config_path):
    if not config_path.exists():
        log(f"Archivo de configuracion no encontrado en {config_path}. Usando defaults.")
        return DEFAULT_CONFIG.copy()

    with open(config_path, "r", encoding="utf-8") as f:
        loaded = json.load(f)

    config = DEFAULT_CONFIG.copy()
    config.update(loaded)
    return config


def cds_start_for_record(record):
    for feature in record.features:
        if feature.type == "CDS":
            return int(feature.location.start)
    raise ValueError(
        "No se encontro feature CDS en el GenBank; use variant.transcript_position "
        "en la configuracion si el input es FASTA."
    )


def resolve_variant(record, config):
    variant = config.get("variant")
    if not variant:
        return None

    if "transcript_position" in variant:
        transcript_position = int(variant["transcript_position"])
    elif "coding_position" in variant:
        transcript_position = cds_start_for_record(record) + int(variant["coding_position"])
    else:
        raise ValueError(
            "La variante debe definir transcript_position o coding_position en primer_config.json."
        )

    return {
        "name": variant.get("name", "variante configurada"),
        "rsid": variant.get("rsid", ""),
        "disease": variant.get("disease", ""),
        "coding_position": variant.get("coding_position"),
        "transcript_position": transcript_position,
        "ref": variant.get("ref", "").upper(),
        "alt": variant.get("alt", "").upper(),
    }


def apply_variant(seq_obj, variant_info):
    if not variant_info:
        return seq_obj

    seq_chars = list(str(seq_obj).upper())
    index = variant_info["transcript_position"] - 1
    if index < 0 or index >= len(seq_chars):
        raise ValueError(
            f"Posicion de variante fuera del transcripto: {variant_info['transcript_position']}"
        )

    ref_base = variant_info["ref"]
    alt_base = variant_info["alt"]
    if len(ref_base) != 1 or len(alt_base) != 1:
        raise ValueError("Este script espera una variante SNV con ref y alt de una base.")

    observed = seq_chars[index]
    if observed != ref_base:
        raise ValueError(
            f"La base de referencia en posicion transcript {variant_info['transcript_position']} "
            f"es {observed}, pero la configuracion esperaba {ref_base}."
        )

    seq_chars[index] = alt_base
    return Seq("".join(seq_chars))


def find_candidate_primers(seq_obj, config, is_reverse=False):
    """Encuentra candidatos de primers en la secuencia dada."""
    min_len = config["min_length"]
    max_len = config["max_length"]
    min_gc = config["min_gc"]
    max_gc = config["max_gc"]
    max_tm = config["max_tm"]
    avoid_gc_ends = config.get("avoid_gc_ends", True)

    seq_str = str(seq_obj)
    candidates = []

    for length in range(min_len, max_len + 1):
        for i in range(len(seq_str) - length + 1):
            sub_seq_str = seq_str[i : i + length]
            sub_seq_obj = Seq(sub_seq_str)

            gc = calculate_gc(sub_seq_str)
            if not (min_gc <= gc <= max_gc):
                continue

            if avoid_gc_ends and not check_terminal_ends(sub_seq_str):
                continue

            tm = calculate_tm(sub_seq_obj)
            if tm > max_tm:
                continue

            if is_reverse:
                start_pos = len(seq_str) - (i + length) + 1
                end_pos = len(seq_str) - i
            else:
                start_pos = i + 1
                end_pos = i + length

            candidates.append(
                {
                    "sequence": sub_seq_str.upper(),
                    "length": length,
                    "gc_percent": round(gc, 2),
                    "tm": round(tm, 2),
                    "start": start_pos,
                    "end": end_pos,
                    "direction": "Reverse" if is_reverse else "Forward",
                }
            )

    return candidates


def distance_to_variant(candidate, variant_position):
    if candidate["start"] <= variant_position <= candidate["end"]:
        return 0
    return min(abs(candidate["start"] - variant_position), abs(candidate["end"] - variant_position))


def candidate_key(candidate):
    return (
        candidate["direction"],
        candidate["start"],
        candidate["end"],
        candidate["sequence"],
    )


def find_variant_flanking_pairs(forward_candidates, reverse_candidates, variant_info, config):
    if not variant_info:
        return []

    variant_pos = variant_info["transcript_position"]
    min_amplicon = config.get("min_amplicon_size", 70)
    max_amplicon = config.get("max_amplicon_size", 180)
    target_tm = config.get("target_tm", 60.0)
    max_pair_tm_diff = config.get("max_pair_tm_diff", 2.0)

    pairs = []
    for forward in forward_candidates:
        if forward["end"] >= variant_pos:
            continue
        for reverse in reverse_candidates:
            if reverse["start"] <= variant_pos:
                continue

            amplicon_size = reverse["end"] - forward["start"] + 1
            if not (min_amplicon <= amplicon_size <= max_amplicon):
                continue

            tm_diff = abs(forward["tm"] - reverse["tm"])
            if tm_diff > max_pair_tm_diff:
                continue

            score = (
                abs(forward["tm"] - target_tm)
                + abs(reverse["tm"] - target_tm)
                + tm_diff
                + abs(amplicon_size - ((min_amplicon + max_amplicon) / 2)) / 100.0
            )
            pairs.append(
                {
                    "forward": forward,
                    "reverse": reverse,
                    "amplicon_start": forward["start"],
                    "amplicon_end": reverse["end"],
                    "amplicon_size": amplicon_size,
                    "variant_position": variant_pos,
                    "tm_diff": round(tm_diff, 2),
                    "score": round(score, 4),
                }
            )

    pairs.sort(key=lambda item: item["score"])
    return pairs


def select_primers(all_candidates, forward_candidates, reverse_candidates, variant_info, config):
    num_primers = config.get("num_primers", 5)
    target_tm = config.get("target_tm", 60.0)
    pairs = find_variant_flanking_pairs(forward_candidates, reverse_candidates, variant_info, config)

    selected = []
    selected_keys = set()

    for pair in pairs:
        for primer in (pair["forward"], pair["reverse"]):
            key = candidate_key(primer)
            if key in selected_keys:
                continue
            selected.append(primer)
            selected_keys.add(key)
            if len(selected) >= num_primers:
                return selected, pairs

    if variant_info:
        variant_pos = variant_info["transcript_position"]
        fallback_candidates = sorted(
            all_candidates,
            key=lambda item: (
                distance_to_variant(item, variant_pos),
                abs(item["tm"] - target_tm),
            ),
        )
    else:
        fallback_candidates = sorted(all_candidates, key=lambda item: abs(item["tm"] - target_tm))

    for primer in fallback_candidates:
        key = candidate_key(primer)
        if key in selected_keys:
            continue
        selected.append(primer)
        selected_keys.add(key)
        if len(selected) >= num_primers:
            break

    return selected, pairs


def format_primer_row(index, primer):
    range_str = f"{primer['start']}-{primer['end']}"
    return (
        f"{index:<3} | {primer['direction']:<8} | {range_str:<11} | "
        f"{primer['length']:<5} | {primer['gc_percent']:<6} | "
        f"{primer['tm']:<7} | {primer['sequence']:<24}"
    )


def append_primer_table(report, title, primers):
    sequence_header = "Secuencia (5' -> 3')"
    report.append(title)
    report.append("-" * 90)
    report.append(
        f"{'N°':<3} | {'Sentido':<8} | {'Rango (bp)':<11} | {'Largo':<5} | "
        f"{'GC %':<6} | {'Tm (°C)':<7} | {sequence_header:<24}"
    )
    report.append("-" * 90)
    for idx, primer in enumerate(primers, start=1):
        report.append(format_primer_row(idx, primer))
    report.append("-" * 90)
    report.append("")


def main():
    root = Path(__file__).parent.resolve()

    input_file = sys.argv[1] if len(sys.argv) > 1 else "NM_000207.gbk"
    config_file = sys.argv[2] if len(sys.argv) > 2 else "primer_config.json"
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "primer_results"

    input_path = root / input_file
    config_path = root / config_file
    output_path = root / output_dir
    os.makedirs(output_path, exist_ok=True)

    config = load_config(config_path)
    log(f"Parametros de diseno cargados de {config_file}:")
    for key, value in config.items():
        if key != "variant":
            log(f"  {key}: {value}")

    if not input_path.exists():
        log(f"Archivo de transcripto no encontrado en {input_path}.")
        sys.exit(1)

    try:
        fmt = "genbank" if input_file.endswith((".gbk", ".gb", ".gbf")) else "fasta"
        record = SeqIO.read(input_path, fmt)
        log(f"Cargada secuencia '{record.id}' ({len(record.seq)} nucleotidos) desde {input_file}.")
        variant_info = resolve_variant(record, config)
        design_seq = apply_variant(record.seq, variant_info)
    except Exception as exc:
        log(f"Error al preparar la secuencia: {exc}")
        sys.exit(1)

    if variant_info:
        log(
            "Variante aplicada: "
            f"{variant_info['name']} en transcript nt {variant_info['transcript_position']} "
            f"({variant_info['ref']}>{variant_info['alt']})."
        )

    forward_candidates = find_candidate_primers(design_seq, config, is_reverse=False)
    reverse_candidates = find_candidate_primers(design_seq.reverse_complement(), config, is_reverse=True)
    all_candidates = forward_candidates + reverse_candidates
    log(
        f"Se encontraron {len(forward_candidates)} candidatos Forward y "
        f"{len(reverse_candidates)} candidatos Reverse."
    )

    selected_primers, primer_pairs = select_primers(
        all_candidates,
        forward_candidates,
        reverse_candidates,
        variant_info,
        config,
    )
    num_pairs = config.get("num_primer_pairs", 5)
    selected_pairs = primer_pairs[:num_pairs]

    results_json = {
        "config": config,
        "variant": variant_info,
        "summary": {
            "total_forward_candidates": len(forward_candidates),
            "total_reverse_candidates": len(reverse_candidates),
            "total_candidates": len(all_candidates),
            "variant_flanking_pairs": len(primer_pairs),
        },
        "selected_primers": selected_primers,
        "selected_primer_pairs": selected_pairs,
    }

    with open(output_path / "primers.json", "w", encoding="utf-8") as f:
        json.dump(results_json, f, indent=2, ensure_ascii=False)

    report = []
    report.append("=" * 70)
    report.append(" DISEÑO DE PRIMERS PARA TRANSCRIPTO INS (NM_000207.3)")
    report.append("=" * 70)
    report.append(f"Secuencia base: {record.id} ({len(record.seq)} bp)")
    if variant_info:
        report.append(
            "Variante aplicada: "
            f"{variant_info['name']} | transcript nt {variant_info['transcript_position']} | "
            f"{variant_info['ref']}>{variant_info['alt']} | rs={variant_info.get('rsid') or 'N/D'}"
        )
    else:
        report.append("Variante aplicada: ninguna (diseño sobre referencia WT)")
    report.append("Parámetros de Diseño:")
    report.append(f"  - Tamaño: {config['min_length']}-{config['max_length']} bp")
    report.append(f"  - Contenido GC: {config['min_gc']}%-{config['max_gc']}%")
    report.append(f"  - Evitar GC en extremos terminales (5' y 3'): {'Sí' if config.get('avoid_gc_ends', True) else 'No'}")
    report.append(f"  - Temperatura de melting máxima: {config['max_tm']}°C")
    report.append(f"  - Amplicón esperado: {config['min_amplicon_size']}-{config['max_amplicon_size']} bp")
    report.append(f"  - Pares que flanquean la variante: {len(primer_pairs)}")
    report.append(f"  - Candidatos totales encontrados: {len(all_candidates)}")
    report.append("=" * 70)
    report.append("")

    append_primer_table(
        report,
        "TOP PRIMERS SELECCIONADOS (priorizando pares que flanquean la variante):",
        selected_primers,
    )

    if selected_pairs:
        report.append("TOP PARES FORWARD/REVERSE QUE FLANQUEAN LA VARIANTE:")
        report.append("-" * 120)
        report.append(
            f"{'N°':<3} | {'Amplicón':<11} | {'Tm diff':<7} | "
            f"{'Forward (5->3)':<26} | {'Reverse (5->3)':<26}"
        )
        report.append("-" * 120)
        for idx, pair in enumerate(selected_pairs, start=1):
            amplicon = f"{pair['amplicon_start']}-{pair['amplicon_end']} ({pair['amplicon_size']} bp)"
            report.append(
                f"{idx:<3} | {amplicon:<11} | {pair['tm_diff']:<7} | "
                f"{pair['forward']['sequence']:<26} | {pair['reverse']['sequence']:<26}"
            )
        report.append("-" * 120)
        report.append("")
    elif variant_info:
        report.append(
            "ADVERTENCIA: no se encontraron pares F/R que flanqueen la variante con los "
            "parámetros actuales. Relajar min_amplicon_size, max_amplicon_size o max_pair_tm_diff."
        )
        report.append("")

    report_text = "\n".join(report)
    print(report_text)

    with open(output_path / "primers_report.txt", "w", encoding="utf-8") as f:
        f.write(report_text + "\n")

    log(f"Diseño completado. Reporte escrito en {output_path / 'primers_report.txt'}")


if __name__ == "__main__":
    main()
