import json
import sys
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def log(msg):
    print(f"[Ex5] {msg}")

def calculate_gc(seq_str):
    """Calcula el porcentaje de GC de una secuencia."""
    g_or_c = sum(1 for c in seq_str if c in ('G', 'C', 'g', 'c'))
    return (g_or_c / len(seq_str)) * 100.0

def calculate_tm(seq_obj):
    """Calcula la temperatura de melting usando Tm_NN o un fallback de sal."""
    try:
        # standard nearest-neighbor calculation
        return mt.Tm_NN(seq_obj)
    except Exception:
        # Fallback: Salt-adjusted formula
        seq_str = str(seq_obj).upper()
        gc_count = seq_str.count('G') + seq_str.count('C')
        return 64.9 + 41.0 * (gc_count - 16.4) / len(seq_str)

def check_terminal_ends(seq_str):
    """Evita tener G o C en los extremos terminales (5' y 3' del primer)."""
    # Ambos extremos deben ser A o T (no G ni C)
    first_char = seq_str[0].upper()
    last_char = seq_str[-1].upper()
    return first_char in ('A', 'T') and last_char in ('A', 'T')

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
    
    # Recorrer todos los largos posibles
    for length in range(min_len, max_len + 1):
        for i in range(len(seq_str) - length + 1):
            sub_seq_str = seq_str[i : i + length]
            sub_seq_obj = Seq(sub_seq_str)
            
            # 1. Porcentaje GC
            gc = calculate_gc(sub_seq_str)
            if not (min_gc <= gc <= max_gc):
                continue
                
            # 2. Extremos terminales (evitar G/C en extremos 5' y 3')
            if avoid_gc_ends and not check_terminal_ends(sub_seq_str):
                continue
                
            # 3. Temperatura de melting
            tm = calculate_tm(sub_seq_obj)
            if tm > max_tm:
                continue
                
            # Calcular posición en el transcripto original (1-indexed)
            if is_reverse:
                # Si es reversa, las coordenadas de unión en el transcripto original son complementarias:
                # El primer reverse se une a la hebra forward desde el final del primer hasta el inicio
                start_pos = len(seq_str) - (i + length) + 1
                end_pos = len(seq_str) - i
            else:
                start_pos = i + 1
                end_pos = i + length
                
            candidates.append({
                "sequence": sub_seq_str.upper(),
                "length": length,
                "gc_percent": round(gc, 2),
                "tm": round(tm, 2),
                "start": start_pos,
                "end": end_pos,
                "direction": "Reverse" if is_reverse else "Forward"
            })
            
    return candidates

def main():
    ROOT = Path(__file__).parent.resolve()
    
    # 1. Leer argumentos de línea de comando o valores por defecto
    input_file = sys.argv[1] if len(sys.argv) > 1 else "NM_000207.gbk"
    config_file = sys.argv[2] if len(sys.argv) > 2 else "primer_config.json"
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "primer_results"
    
    input_path = ROOT / input_file
    config_path = ROOT / config_file
    output_path = ROOT / output_dir
    os.makedirs(output_path, exist_ok=True)
    
    # 2. Cargar archivo de configuración
    if not config_path.exists():
        log(f"Archivo de configuración no encontrado en {config_path}. Usando valores por defecto.")
        config = {
            "min_length": 18,
            "max_length": 24,
            "min_gc": 50.0,
            "max_gc": 60.0,
            "avoid_gc_ends": True,
            "max_tm": 67.0,
            "num_primers": 5
        }
    else:
        with open(config_path, "r", encoding="utf-8") as f:
            config = json.load(f)
            
    log(f"Parámetros de diseño cargados de {config_file}:")
    for k, v in config.items():
        log(f"  {k}: {v}")
        
    # 3. Leer secuencia de transcripto (GenBank o FASTA)
    if not input_path.exists():
        log(f"Archivo de transcripto no encontrado en {input_path}.")
        sys.exit(1)
        
    try:
        # Detect format
        fmt = "genbank" if input_file.endswith((".gbk", ".gb", ".gbf")) else "fasta"
        record = SeqIO.read(input_path, fmt)
        log(f"Cargada secuencia '{record.id}' ({len(record.seq)} nucleótidos) desde {input_file}.")
    except Exception as e:
        log(f"Error al leer la secuencia: {e}")
        sys.exit(1)
        
    # 4. Encontrar candidatos
    # Hebra Forward
    forward_candidates = find_candidate_primers(record.seq, config, is_reverse=False)
    # Hebra Reverse (usamos el reverse complement de la secuencia original)
    rev_complement_seq = record.seq.reverse_complement()
    reverse_candidates = find_candidate_primers(rev_complement_seq, config, is_reverse=True)
    
    all_candidates = forward_candidates + reverse_candidates
    log(f"Se encontraron {len(forward_candidates)} candidatos Forward y {len(reverse_candidates)} candidatos Reverse.")
    
    # 5. Seleccionar y rankear
    # Criterio de ranking: cercanía a la temperatura de hibridación ideal en PCR (60°C)
    # y que estén bien distribuidos. Ordenamos por |Tm - 60| ascendente.
    all_candidates.sort(key=lambda x: abs(x["tm"] - 60.0))
    
    # Separar en top Forward y top Reverse
    top_forward = [c for c in all_candidates if c["direction"] == "Forward"][:5]
    top_reverse = [c for c in all_candidates if c["direction"] == "Reverse"][:5]
    
    # Seleccionar top 5 general (o 5 pares si es posible, pero mostremos los 5 mejores Forward y Reverse)
    selected_primers = all_candidates[:config.get("num_primers", 5)]
    
    # Guardar resultados en JSON
    results_json = {
        "config": config,
        "summary": {
            "total_forward_candidates": len(forward_candidates),
            "total_reverse_candidates": len(reverse_candidates),
            "total_candidates": len(all_candidates)
        },
        "selected_primers": selected_primers,
        "top_5_forward": top_forward,
        "top_5_reverse": top_reverse
    }
    
    with open(output_path / "primers.json", "w", encoding="utf-8") as f:
        json.dump(results_json, f, indent=2, ensure_ascii=False)
        
    # Generar reporte escrito en texto plano
    report = []
    report.append("=" * 70)
    report.append(" DISEÑO DE PRIMERS PARA TRANSCRIPTO INS (NM_000207.3)")
    report.append("=" * 70)
    report.append(f"Secuencia: {record.id} ({len(record.seq)} bp)")
    report.append(f"Parámetros de Diseño:")
    report.append(f"  - Tamaño: {config['min_length']}-{config['max_length']} bp")
    report.append(f"  - Contenido GC: {config['min_gc']}%-{config['max_gc']}%")
    report.append(f"  - Evitar GC en extremos terminales (5' y 3'): {'Sí' if config.get('avoid_gc_ends', True) else 'No'}")
    report.append(f"  - Temperatura de melting máxima: {config['max_tm']}°C")
    report.append(f"  - Candidatos totales encontrados: {len(all_candidates)}")
    report.append("=" * 70)
    report.append("")
    
    report.append("TOP 5 PRIMERS SELECCIONADOS (Ordenados por Tm más cercana a 60°C):")
    report.append("-" * 90)
    report.append(f"{'N°':<3} | {'Sentido':<8} | {'Rango (bp)':<11} | {'Largo':<5} | {'GC %':<6} | {'Tm (°C)':<7} | {'Secuencia (5\' -> 3\')':<24}")
    report.append("-" * 90)
    for idx, p in enumerate(selected_primers):
        range_str = f"{p['start']}-{p['end']}"
        report.append(f"{idx+1:<3} | {p['direction']:<8} | {range_str:<11} | {p['length']:<5} | {p['gc_percent']:<6} | {p['tm']:<7} | {p['sequence']:<24}")
    report.append("-" * 90)
    report.append("")
    
    report.append("TOP 5 PRIMERS FORWARD CANDIDATOS:")
    report.append("-" * 90)
    report.append(f"{'N°':<3} | {'Sentido':<8} | {'Rango (bp)':<11} | {'Largo':<5} | {'GC %':<6} | {'Tm (°C)':<7} | {'Secuencia (5\' -> 3\')':<24}")
    report.append("-" * 90)
    for idx, p in enumerate(top_forward):
        range_str = f"{p['start']}-{p['end']}"
        report.append(f"{idx+1:<3} | {p['direction']:<8} | {range_str:<11} | {p['length']:<5} | {p['gc_percent']:<6} | {p['tm']:<7} | {p['sequence']:<24}")
    report.append("-" * 90)
    report.append("")
    
    report.append("TOP 5 PRIMERS REVERSE CANDIDATOS:")
    report.append("-" * 90)
    report.append(f"{'N°':<3} | {'Sentido':<8} | {'Rango (bp)':<11} | {'Largo':<5} | {'GC %':<6} | {'Tm (°C)':<7} | {'Secuencia (5\' -> 3\')':<24}")
    report.append("-" * 90)
    for idx, p in enumerate(top_reverse):
        range_str = f"{p['start']}-{p['end']}"
        report.append(f"{idx+1:<3} | {p['direction']:<8} | {range_str:<11} | {p['length']:<5} | {p['gc_percent']:<6} | {p['tm']:<7} | {p['sequence']:<24}")
    report.append("-" * 90)
    
    report_text = "\n".join(report)
    print(report_text)
    
    with open(output_path / "primers_report.txt", "w", encoding="utf-8") as f:
        f.write(report_text + "\n")
        
    log(f"Diseño completado. Reporte escrito en {output_path / 'primers_report.txt'}")

if __name__ == "__main__":
    main()
