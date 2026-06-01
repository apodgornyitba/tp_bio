import os
import sys
import subprocess
import urllib.request
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# URL for PROSITE database
PROSITE_URL = "https://ftp.expasy.org/databases/prosite/prosite.dat"
PROSITE_FILE = "prosite.dat"

def log(msg):
    print(f"[Ex4] {msg}")

def download_prosite(dest_path):
    """Downloads prosite.dat if it does not exist."""
    dest = Path(dest_path)
    if dest.exists() and dest.stat().st_size > 10 * 1024 * 1024:
        log(f"prosite.dat ya existe en {dest} ({dest.stat().st_size / 1024 / 1024:.2f} MB). Saltando descarga.")
        return True
    
    log(f"Descargando prosite.dat de {PROSITE_URL} (esto puede tardar unos segundos)...")
    try:
        def progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(100.0, downloaded * 100.0 / total_size) if total_size > 0 else 0
            sys.stdout.write(f"\r  Progreso: {percent:.1f}% ({downloaded / 1024 / 1024:.1f}/{total_size / 1024 / 1024:.1f} MB)")
            sys.stdout.flush()
        
        urllib.request.urlretrieve(PROSITE_URL, dest, reporthook=progress)
        print() # New line after progress
        log("Descarga completada con éxito.")
        return True
    except Exception as e:
        print()
        log(f"Error al descargar prosite.dat: {e}")
        return False

def extract_nucleotide_fasta(gbk_path, fasta_path):
    """Extrae la secuencia de nucleótidos completa del GenBank y la guarda en FASTA."""
    log(f"Extrayendo secuencia de nucleótidos de {gbk_path} a {fasta_path}...")
    try:
        record = SeqIO.read(gbk_path, "genbank")
        SeqIO.write(record, fasta_path, "fasta")
        log(f"Secuencia extraída con éxito: {len(record.seq)} nucleótidos.")
        return record
    except Exception as e:
        log(f"Error al extraer nucleótidos: {e}")
        sys.exit(1)

def check_emboss_program(name):
    """Verifica si un ejecutable de EMBOSS está disponible en el PATH."""
    try:
        # On Windows, which doesn't exist, we can use shutil.which or just run which on Linux
        import shutil
        return shutil.which(name) is not None
    except Exception:
        return False

def run_getorf_emboss(input_fasta, output_fasta):
    """Ejecuta getorf de EMBOSS para obtener los ORFs."""
    log("Ejecutando getorf de EMBOSS...")
    cmd = ["getorf", "-sequence", str(input_fasta), "-outseq", str(output_fasta), "-find", "0"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log(f"getorf completado con éxito. Resultados en {output_fasta}")
        return True
    except Exception as e:
        log(f"Error al ejecutar getorf de EMBOSS: {e}")
        return False

def run_patmatmotifs_emboss(input_fasta, output_file, prosite_dir):
    """Ejecuta patmatmotifs de EMBOSS para analizar dominios."""
    log("Ejecutando patmatmotifs de EMBOSS...")
    # Set EMBOSS_DATA variable to directory of prosite.dat
    env = os.environ.copy()
    env["EMBOSS_DATA"] = str(prosite_dir)
    
    cmd = ["patmatmotifs", "-sequence", str(input_fasta), "-outfile", str(output_file), "-prune"]
    try:
        subprocess.run(cmd, env=env, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        log(f"patmatmotifs completado con éxito. Resultados en {output_file}")
        return True
    except Exception as e:
        log(f"Error al ejecutar patmatmotifs de EMBOSS: {e}")
        return False

# ================= Fallbacks en Python Puro =================

def get_reading_frames(seq):
    """Devuelve los 6 marcos de lectura (3 forward + 3 reverse complement)."""
    frames = []
    for i in range(3):
        frames.append((f"F{i+1}", i, seq[i:], False)) # Name, offset, seq, is_reverse
    
    rev_seq = seq.reverse_complement()
    for i in range(3):
        frames.append((f"R{i+1}", i, rev_seq[i:], True))
    return frames

def run_getorf_fallback(input_fasta, output_fasta):
    """Implementación fallback en Python para getorf -find 0 (traducción entre STOPs)."""
    log("Iniciando Fallback en Python de getorf (traducción de regiones entre STOPs)...")
    record = SeqIO.read(input_fasta, "fasta")
    seq = record.seq
    
    orf_records = []
    frames = get_reading_frames(seq)
    
    orf_counter = 1
    for frame_name, offset, frame_seq, is_rev in frames:
        # Pad to multiple of 3
        length = len(frame_seq)
        padded = frame_seq[: length - (length % 3)]
        translation = str(padded.translate())
        
        # Split by stop codon '*'
        # getorf includes the last part after the last stop codon, and the first part before the first stop codon.
        parts = translation.split("*")
        
        current_aa_pos = 0
        for idx, part in enumerate(parts):
            # We skip empty ORFs (e.g. consecutive stop codons)
            if not part:
                current_aa_pos += 1
                continue
            
            # Start position in nucleotide sequence (0-indexed)
            nt_start = offset + 3 * current_aa_pos
            nt_end = nt_start + 3 * len(part) + 3 # including stop codon
            
            if is_rev:
                # Reverse coordinates are mapped from the end of the original sequence
                orig_start = len(seq) - nt_end
                orig_end = len(seq) - nt_start
            else:
                orig_start = nt_start
                orig_end = nt_end
            
            # Formatting the ID similar to EMBOSS
            orf_id = f"{record.id}_{orf_counter}"
            description = f"[{orig_start + 1} - {orig_end}] Frame: {frame_name} Translation of ORF between STOPs"
            
            orf_records.append(SeqRecord(Seq(part), id=orf_id, description=description))
            orf_counter += 1
            current_aa_pos += len(part) + 1 # +1 for the stop codon
            
    SeqIO.write(orf_records, output_fasta, "fasta")
    log(f"Fallback getorf finalizado. Guardados {len(orf_records)} ORFs en {output_fasta}")
    return True

def convert_prosite_pattern_to_regex(pattern):
    """Convierte un patrón de PROSITE en una expresión regular Python."""
    # Example: C-C-{P}-x(2)-C-[STDNE]-[NH]-[FYW]-C
    # x -> .
    # x(2) -> .{2}
    # x(2,4) -> .{2,4}
    # [STDNE] -> [STDNE]
    # {P} -> [^P]
    # - -> ''
    # . -> end of pattern (already stripped)
    
    parts = pattern.split("-")
    regex_str = ""
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if part == "x":
            regex_str += "."
        elif part.startswith("x(") and part.endswith(")"):
            inner = part[2:-1]
            regex_str += f".{{{inner}}}"
        elif part.startswith("[") and part.endswith("]"):
            regex_str += part
        elif part.startswith("{") and part.endswith("}"):
            inner = part[1:-1]
            regex_str += f"[^{inner}]"
        else:
            # Literal amino acids
            # Can be multiple like A or just C
            regex_str += part
    
    return re.compile(regex_str)

def parse_prosite_database(prosite_path):
    """Parsea el archivo prosite.dat y extrae los patrones activos."""
    log(f"Parseando base de datos PROSITE ({prosite_path})...")
    motifs = []
    current_motif = {}
    
    with open(prosite_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("ID   "):
                # ID   INSULIN; PATTERN.
                current_motif = {
                    "id": line[5:].split(";")[0].strip(),
                    "type": line[5:].split(";")[1].strip() if ";" in line else ""
                }
            elif line.startswith("AC   "):
                current_motif["ac"] = line[5:].split(";")[0].strip()
            elif line.startswith("DE   "):
                current_motif["de"] = line[5:].strip()
            elif line.startswith("PA   "):
                pa_line = line[5:].strip()
                if "pa" in current_motif:
                    current_motif["pa"] += pa_line
                else:
                    current_motif["pa"] = pa_line
            elif line.startswith("//"):
                # End of entry
                if "pa" in current_motif and "PATTERN" in current_motif.get("type", ""):
                    pattern = current_motif["pa"]
                    if pattern.endswith("."):
                        pattern = pattern[:-1]
                    current_motif["pa_clean"] = pattern
                    try:
                        current_motif["regex"] = convert_prosite_pattern_to_regex(pattern)
                        motifs.append(current_motif)
                    except Exception as e:
                        # Ignore failed patterns (e.g. complex rules)
                        pass
                current_motif = {}
                
    log(f"Parsea completado. Se cargaron {len(motifs)} patrones activos.")
    return motifs

def run_patmatmotifs_fallback(input_fasta, output_file, prosite_path):
    """Implementación fallback en Python para patmatmotifs."""
    log("Iniciando Fallback en Python de patmatmotifs (escaneo de dominios)...")
    motifs = parse_prosite_database(prosite_path)
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    results = []
    results.append("#=======================================")
    results.append("#")
    results.append(f"# Aligned_sequences: {len(records)}")
    results.append("# Database: PROSITE (Python Fallback)")
    results.append("#")
    results.append("#=======================================")
    results.append("")
    
    hits_count = 0
    for record in records:
        seq_str = str(record.seq)
        for motif in motifs:
            regex = motif["regex"]
            # Find all overlapping/non-overlapping matches
            # Using finditer for non-overlapping hits (which matches EMBOSS -prune behavior)
            for match in regex.finditer(seq_str):
                start = match.start() + 1
                end = match.end()
                hit_seq = match.group(0)
                
                # Format hit sequence with hyphens like EMBOSS
                hit_formatted = "-".join(list(hit_seq))
                
                results.append("#=======================================")
                results.append("#")
                results.append(f"# {record.id}: {motif['id']}")
                results.append("#")
                results.append(f"# Sequence: {record.id}")
                results.append(f"# Accession: {motif['ac']}")
                results.append(f"# Description: {motif['de']}")
                results.append(f"# Start: {start}")
                results.append(f"# End: {end}")
                results.append(f"# Hit: {hit_formatted}")
                results.append(f"# Pattern: {motif['pa_clean']}")
                results.append("#")
                results.append("#=======================================")
                results.append("")
                hits_count += 1
                
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(results) + "\n")
        
    log(f"Fallback patmatmotifs finalizado. Se encontraron {hits_count} hits de dominios. Resultados en {output_file}")
    return True

# ================= Main =================

def main():
    if len(sys.argv) < 3:
        input_gbk = "NM_000207.gbk"
        output_dir = "emboss_results"
    else:
        input_gbk = sys.argv[1]
        output_dir = sys.argv[2]
        
    os.makedirs(output_dir, exist_ok=True)
    
    ROOT = Path(__file__).parent.resolve()
    input_gbk_path = ROOT / input_gbk
    output_dir_path = ROOT / output_dir
    
    # 1. Download prosite.dat
    prosite_path = ROOT / PROSITE_FILE
    download_success = download_prosite(prosite_path)
    
    # 2. Extract nucleotide sequence
    nucleotides_fasta = output_dir_path / "NM_000207_nucleotides.fasta"
    extract_nucleotide_fasta(input_gbk_path, nucleotides_fasta)
    
    # 3. Running ORFs step
    orfs_fasta = output_dir_path / "NM_000207_orfs.fasta"
    emboss_getorf_available = check_emboss_program("getorf")
    
    if emboss_getorf_available:
        log("EMBOSS 'getorf' encontrado en el sistema.")
        success = run_getorf_emboss(nucleotides_fasta, orfs_fasta)
        if not success:
            log("Fallo al ejecutar EMBOSS 'getorf'. Usando fallback de Python.")
            run_getorf_fallback(nucleotides_fasta, orfs_fasta)
    else:
        log("EMBOSS 'getorf' NO está instalado en el sistema.")
        log("  [Sugerencia] Para instalarlo: sudo apt-get update && sudo apt-get install -y emboss")
        run_getorf_fallback(nucleotides_fasta, orfs_fasta)
        
    # 4. Running Domains step
    domains_file = output_dir_path / "NM_000207_domains.patmatmotifs"
    emboss_patmatmotifs_available = check_emboss_program("patmatmotifs")
    
    if download_success:
        if emboss_patmatmotifs_available:
            log("EMBOSS 'patmatmotifs' encontrado en el sistema.")
            success = run_patmatmotifs_emboss(orfs_fasta, domains_file, ROOT)
            if not success:
                log("Fallo al ejecutar EMBOSS 'patmatmotifs'. Usando fallback de Python.")
                run_patmatmotifs_fallback(orfs_fasta, domains_file, prosite_path)
        else:
            log("EMBOSS 'patmatmotifs' NO está instalado en el sistema.")
            log("  [Sugerencia] Para instalarlo: sudo apt-get update && sudo apt-get install -y emboss")
            run_patmatmotifs_fallback(orfs_fasta, domains_file, prosite_path)
    else:
        log("No se pudo ejecutar el análisis de dominios porque prosite.dat no está disponible y falló la descarga.")
        sys.exit(1)

    log("Ejercicio 4 finalizado con éxito.")

if __name__ == "__main__":
    main()
