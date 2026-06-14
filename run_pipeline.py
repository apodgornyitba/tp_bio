#!/usr/bin/env python3
"""
Pipeline multiplataforma (Windows, macOS, Linux).
Uso:
  python run_pipeline.py
  python run_pipeline.py --blast-mode both --require-emboss
"""
import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from platform_tools import default_blast_db, find_blastp, python_command, project_root

ROOT = project_root()
LOG_FILE = ROOT / "pipeline.log"
BLAST_MODES = ("remote", "local", "both")
MSA_SOURCES = ("auto", "remote", "local")


def log(msg):
    line = str(msg)
    print(line)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def run_step(name, cmd, required_output=None):
    log(f"\n[INFO] {name}")
    log(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=ROOT)
    if result.returncode != 0:
        log(f"[ERROR] Fallo: {name}")
        sys.exit(result.returncode)
    if required_output and not Path(required_output).exists():
        log(f"[ERROR] No se genero: {required_output}")
        sys.exit(1)
    log(f"[OK] {name}")


def env_flag(name, default=False):
    value = os.environ.get(name)
    if value is None:
        return default
    return value.lower() in ("1", "true", "yes", "on")


def choice(value, valid_values):
    normalized = value.lower()
    if normalized not in valid_values:
        valid = ", ".join(valid_values)
        raise argparse.ArgumentTypeError(f"valor invalido '{value}'. Usar: {valid}")
    return normalized


def parse_args():
    parser = argparse.ArgumentParser(description="Pipeline TP Bioinformatica")
    parser.add_argument(
        "--blast-mode",
        default=os.environ.get("BLAST_MODE", "remote"),
        type=lambda value: choice(value, BLAST_MODES),
        help="remote, local o both. Default: env BLAST_MODE o remote.",
    )
    parser.add_argument(
        "--msa-source",
        default=os.environ.get("MSA_SOURCE", "auto"),
        type=lambda value: choice(value, MSA_SOURCES),
        help="Fuente BLAST para Ejercicio 3: auto, remote o local. Default: env MSA_SOURCE o auto.",
    )
    parser.add_argument(
        "--entrez-email",
        default=os.environ.get("ENTREZ_EMAIL", "estudiante@itba.edu.ar"),
        help="Email para NCBI Entrez. Default: env ENTREZ_EMAIL.",
    )
    parser.add_argument(
        "--blast-db",
        default=os.environ.get("BLAST_DB"),
        help="Prefijo de la base BLAST local. Default: env BLAST_DB o data/swissprot_db.",
    )

    blast_group = parser.add_mutually_exclusive_group()
    blast_group.add_argument(
        "--skip-blast",
        dest="skip_blast",
        action="store_true",
        default=None,
        help="Omitir Ejercicio 2 y reutilizar resultados BLAST existentes.",
    )
    blast_group.add_argument(
        "--run-blast",
        dest="skip_blast",
        action="store_false",
        help="Forzar corrida BLAST aunque SKIP_BLAST este seteado.",
    )

    emboss_group = parser.add_mutually_exclusive_group()
    emboss_group.add_argument(
        "--require-emboss",
        dest="require_emboss",
        action="store_true",
        default=None,
        help="Fallar si getorf/patmatmotifs nativos no estan disponibles.",
    )
    emboss_group.add_argument(
        "--allow-emboss-fallback",
        dest="require_emboss",
        action="store_false",
        help="Permitir fallback Python de EMBOSS para desarrollo.",
    )
    return parser.parse_args()


def selected_msa_artifacts(blast_mode, msa_source):
    source = msa_source
    if source == "auto":
        source = "local" if blast_mode == "local" else "remote"

    if source == "local":
        return source, ROOT / "blast_results_local.xml", ROOT / "query_best_local.fasta"
    return source, ROOT / "blast_results.xml", ROOT / "query_best.fasta"


def main():
    args = parse_args()
    py = python_command()
    email = args.entrez_email
    blast_mode = args.blast_mode
    skip_blast = env_flag("SKIP_BLAST") if args.skip_blast is None else args.skip_blast
    require_emboss = (
        env_flag("REQUIRE_EMBOSS") if args.require_emboss is None else args.require_emboss
    )

    os.environ["ENTREZ_EMAIL"] = email
    os.environ["BLAST_MODE"] = blast_mode
    os.environ["MSA_SOURCE"] = args.msa_source
    os.environ["SKIP_BLAST"] = "1" if skip_blast else "0"
    os.environ["REQUIRE_EMBOSS"] = "1" if require_emboss else "0"
    if args.blast_db:
        os.environ["BLAST_DB"] = args.blast_db

    LOG_FILE.write_text(f"=== Pipeline {datetime.now().isoformat()} ===\n", encoding="utf-8")

    log("=" * 50)
    log(" PIPELINE BIOINFORMATICA (multiplataforma)")
    log(f" Python: {py}")
    log(f" BLAST_MODE: {blast_mode}")
    log(f" MSA_SOURCE: {args.msa_source}")
    log(f" ENTREZ_EMAIL: {email}")
    if args.blast_db:
        log(f" BLAST_DB: {args.blast_db}")
    log(f" SKIP_BLAST: {skip_blast}")
    log(f" REQUIRE_EMBOSS: {require_emboss}")
    log("=" * 50)

    # Paso 0
    run_step("Descarga GenBank", [py, "fetch_data.py"], ROOT / "NM_000207.gbk")

    # Paso 1
    run_step(
        "Ejercicio 1 - 6 marcos",
        [py, "Ex1.py", "NM_000207.gbk", "NM_000207_frames.fasta"],
        ROOT / "NM_000207_frames.fasta",
    )

    # Paso 2
    if not skip_blast:
        if blast_mode in ("remote", "both"):
            log("[INFO] BLAST remoto (6 consultas, ~30-60 min)")
            run_step(
                "Ejercicio 2 - BLAST remoto",
                [py, "Ex2_a.py", "NM_000207_frames.fasta", "blast_results"],
                ROOT / "blast_results.xml",
            )

        if blast_mode in ("local", "both"):
            db = default_blast_db()
            if not Path(f"{db}.phr").exists():
                log("[INFO] Base BLAST local no encontrada. Ejecutando prepare_blast_db.py ...")
                run_step("Preparar base Swiss-Prot", [py, "prepare_blast_db.py"])
            if not find_blastp():
                log("[ERROR] blastp no instalado para modo local.")
                log("       Ver README.md (Windows/Linux/macOS).")
                sys.exit(1)
            run_step(
                "Ejercicio 2 - BLAST local",
                [py, "Ex2_local.py", "NM_000207_frames.fasta", "blast_results_local"],
                ROOT / "blast_results_local.xml",
            )
    else:
        log("[INFO] SKIP_BLAST=1: omitiendo BLAST")

    msa_source, blast_xml, query_fasta = selected_msa_artifacts(blast_mode, args.msa_source)
    log(f"[INFO] Ejercicio 3 usara resultados BLAST '{msa_source}': {blast_xml.name}, {query_fasta.name}")

    missing_inputs = [path.name for path in (blast_xml, query_fasta) if not path.exists()]
    if missing_inputs:
        log(f"[ERROR] Faltan resultados BLAST para Ejercicio 3: {', '.join(missing_inputs)}")
        log("       Ejecutar BLAST correspondiente o elegir --msa-source correcto.")
        sys.exit(1)

    # Paso 3
    run_step(
        "Ejercicio 3 - MSA",
        [py, "Ex3.py", str(blast_xml.name), str(query_fasta.name)],
        ROOT / "msa_output.afa",
    )

    # Paso 4
    ex4_cmd = [py, "Ex4.py", "NM_000207.gbk", "emboss_results"]
    if require_emboss:
        ex4_cmd.append("--require-emboss")
    run_step(
        "Ejercicio 4 - EMBOSS y Dominios",
        ex4_cmd,
        ROOT / "emboss_results" / "NM_000207_domains.patmatmotifs",
    )

    # Paso 5
    run_step(
        "Ejercicio 5 - Diseño de Primers",
        [py, "Ex5.py", "NM_000207.gbk", "primer_config.json", "primer_results"],
        ROOT / "primer_results" / "primers.json",
    )

    log("\n" + "=" * 50)
    log(" PIPELINE FINALIZADO")
    log(" Archivos generados:")
    log("  - NM_000207.gbk, NM_000207_frames.fasta")
    log("  - blast_results* (XML y mejores consultas)")
    log("  - msa_output.afa (alineamiento múltiple)")
    log("  - emboss_results/ (NM_000207_nucleotides.fasta, NM_000207_orfs.fasta, NM_000207_domains.patmatmotifs)")
    log("  - primer_results/ (primers.json, primers_report.txt)")
    if blast_mode == "both":
        log("  Remoto: blast_results/ | Local: blast_results_local/")
    log("=" * 50)


if __name__ == "__main__":
    main()
