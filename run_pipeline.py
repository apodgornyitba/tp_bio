#!/usr/bin/env python3
"""
Pipeline multiplataforma (Windows, macOS, Linux).
Uso:
  python run_pipeline.py
  BLAST_MODE=both python run_pipeline.py
"""
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from platform_tools import default_blast_db, find_blastp, python_command, project_root

ROOT = project_root()
RESULTS = ROOT / "results"
LOG_FILE = RESULTS / "pipeline.log"


def log(msg):
    line = str(msg)
    print(line)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def run_step(name, cmd, required_output=None):
    log(f"\n[INFO] {name}")
    log(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=RESULTS)
    if result.returncode != 0:
        log(f"[ERROR] Fallo: {name}")
        sys.exit(result.returncode)
    if required_output and not Path(required_output).exists():
        log(f"[ERROR] No se genero: {required_output}")
        sys.exit(1)
    log(f"[OK] {name}")


def main():
    py = python_command()
    email = os.environ.get("ENTREZ_EMAIL", "estudiante@itba.edu.ar")
    blast_mode = os.environ.get("BLAST_MODE", "remote").lower()  # remote | local | both
    skip_blast = os.environ.get("SKIP_BLAST", "").lower() in ("1", "true", "yes")

    os.environ["ENTREZ_EMAIL"] = email
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOG_FILE.write_text(f"=== Pipeline {datetime.now().isoformat()} ===\n", encoding="utf-8")

    log("=" * 50)
    log(" PIPELINE BIOINFORMATICA (multiplataforma)")
    log(f" Python: {py}")
    log(f" BLAST_MODE: {blast_mode}")
    log(f" ENTREZ_EMAIL: {email}")
    log("=" * 50)

    # Paso 0
    run_step("Descarga GenBank", [py, str(ROOT / "fetch_data.py")], RESULTS / "NM_000207.gbk")

    # Paso 1
    run_step(
        "Ejercicio 1 - 6 marcos",
        [py, str(ROOT / "Ex1.py"), "NM_000207.gbk", "NM_000207_frames.fasta"],
        RESULTS / "NM_000207_frames.fasta",
    )

    # Paso 2
    if not skip_blast:
        if blast_mode in ("remote", "both"):
            log("[INFO] BLAST remoto (6 consultas, ~30-60 min)")
            run_step(
                "Ejercicio 2 - BLAST remoto",
                [py, str(ROOT / "Ex2_a.py"), "NM_000207_frames.fasta", "blast_results"],
                RESULTS / "blast_results.xml",
            )

        if blast_mode in ("local", "both"):
            db = default_blast_db()
            if not Path(f"{db}.phr").exists():
                log("[INFO] Base BLAST local no encontrada. Ejecutando prepare_blast_db.py ...")
                run_step("Preparar base Swiss-Prot", [py, str(ROOT / "prepare_blast_db.py")])
            if not find_blastp():
                log("[ERROR] blastp no instalado para modo local.")
                log("       Ver README.md (Windows/Linux/macOS).")
                sys.exit(1)
            run_step(
                "Ejercicio 2 - BLAST local",
                [py, str(ROOT / "Ex2_local.py"), "NM_000207_frames.fasta", "blast_results_local"],
                RESULTS / "blast_results_local.xml",
            )
    else:
        log("[INFO] SKIP_BLAST=1: omitiendo BLAST")

    blast_xml = RESULTS / "blast_results.xml"
    query_fasta = RESULTS / "query_best.fasta"
    if blast_mode == "local" and not blast_xml.exists():
        blast_xml = RESULTS / "blast_results_local.xml"
        query_fasta = RESULTS / "query_best_local.fasta"

    if not blast_xml.exists():
        log("[ERROR] No hay resultados BLAST para Ejercicio 3.")
        sys.exit(1)

    # Paso 3
    run_step(
        "Ejercicio 3 - MSA",
        [py, str(ROOT / "Ex3.py"), str(blast_xml.name), str(query_fasta.name)],
        RESULTS / "msa_output.afa",
    )

    log("\n" + "=" * 50)
    log(" PIPELINE FINALIZADO")
    log(" Carpeta de resultados: results/")
    log(" Archivos: NM_000207.gbk, NM_000207_frames.fasta,")
    log("  blast_results*, query_best*.fasta, msa_output.afa")
    if blast_mode == "both":
        log("  Remoto: blast_results/ | Local: blast_results_local/")
    log("=" * 50)


if __name__ == "__main__":
    main()
