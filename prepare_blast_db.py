#!/usr/bin/env python3
"""
Descarga Swiss-Prot (FASTA) y crea base local para BLAST+.
Uso: python prepare_blast_db.py
"""
import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from platform_tools import default_blast_db, find_makeblastdb, print_install_hints, project_root

SWISSPROT_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz"


def download_swissprot(dest_dir: Path):
    gz_path = dest_dir / "swissprot.gz"
    fasta_path = dest_dir / "swissprot"

    dest_dir.mkdir(parents=True, exist_ok=True)
    if fasta_path.exists():
        print(f"[OK] Ya existe {fasta_path}")
        return fasta_path

    print(f"Descargando {SWISSPROT_URL} ...")
    urllib.request.urlretrieve(SWISSPROT_URL, gz_path)
    print("Descomprimiendo...")
    with gzip.open(gz_path, "rb") as src, open(fasta_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    gz_path.unlink(missing_ok=True)
    print(f"[OK] FASTA guardado en {fasta_path}")
    return fasta_path


def build_database(fasta_path: Path, db_prefix: Path):
    makeblastdb = find_makeblastdb()
    if not makeblastdb:
        print("ERROR: makeblastdb no encontrado.")
        print_install_hints("blast")
        sys.exit(1)

    if (db_prefix.parent / f"{db_prefix.name}.phr").exists() or db_prefix.with_suffix(".phr").exists():
        print(f"[OK] Base BLAST ya existe en {db_prefix}")
        return

    print(f"Creando base BLAST con {makeblastdb} ...")
    subprocess.run(
        [
            makeblastdb,
            "-in", str(fasta_path),
            "-dbtype", "prot",
            "-out", str(db_prefix),
            "-title", "swissprot_local",
        ],
        check=True,
    )
    print(f"[OK] Base creada: {db_prefix}")


def main():
    data_dir = project_root() / "data"
    fasta_path = download_swissprot(data_dir)
    db_prefix = default_blast_db()
    build_database(fasta_path, db_prefix)
    print()
    print("Configurar para el pipeline:")
    print(f"  export BLAST_DB={db_prefix}")
    print("  export BLAST_MODE=both   # remoto + local")
    print("  ./run_pipeline.sh")


if __name__ == "__main__":
    main()
