"""
Utilidades multiplataforma para encontrar ejecutables (Windows, macOS, Linux).
"""
import os
import platform
import shutil
import sys
from pathlib import Path


def project_root():
    return Path(__file__).resolve().parent


def extra_bin_dirs():
    """Rutas adicionales segun sistema operativo."""
    dirs = []
    if platform.system() == "Windows":
        program_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        dirs.extend([
            Path(program_files) / "NCBI" / "blast-2.16.0+" / "bin",
            Path(program_files) / "NCBI" / "blast-2.15.0+" / "bin",
            Path(os.environ.get("LOCALAPPDATA", "")) / "NCBI" / "blast" / "bin",
            Path.home() / "blast" / "bin",
        ])
    elif platform.system() == "Darwin":
        dirs.extend([
            Path("/opt/homebrew/bin"),
            Path("/usr/local/bin"),
        ])
    else:
        dirs.extend([
            Path("/usr/bin"),
            Path("/usr/local/bin"),
        ])
    return [d for d in dirs if d and str(d)]


def find_executable(*names):
    """Busca ejecutable en PATH y rutas tipicas del SO."""
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    win_names = []
    if platform.system() == "Windows":
        win_names = [f"{n}.exe" for n in names]
    for name in list(names) + win_names:
        for directory in extra_bin_dirs():
            candidate = directory / name
            if candidate.is_file():
                return str(candidate)
    return None


def find_blastp():
    return find_executable("blastp")


def find_makeblastdb():
    return find_executable("makeblastdb", "formatdb")


def find_msa_tool():
    """Devuelve (nombre, ruta) para muscle o mafft."""
    muscle = find_executable("muscle", "muscle5")
    if muscle:
        return "muscle", muscle
    mafft = find_executable("mafft")
    if mafft:
        return "mafft", mafft
    return None, None


def default_blast_db():
    """Ruta por defecto de la base Swiss-Prot formateada."""
    env = os.environ.get("BLAST_DB")
    if env:
        return Path(env)
    return project_root() / "data" / "swissprot_db"


def python_command():
    """Comando python multiplataforma."""
    return os.environ.get("PYTHON", sys.executable)


def print_install_hints(tool="general"):
    system = platform.system()
    print(f"Sistema detectado: {system} ({platform.machine()})")
    if tool in ("blast", "general"):
        print("BLAST+:")
        if system == "Windows":
            print("  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
            print("  Agregar carpeta bin al PATH")
        elif system == "Darwin":
            print("  Descargar BLAST+ desde NCBI o: conda install -c bioconda blast")
        else:
            print("  sudo apt install ncbi-blast+   # Debian/Ubuntu")
            print("  conda install -c bioconda blast")
    if tool in ("msa", "general"):
        print("MSA (MAFFT o MUSCLE):")
        if system == "Windows":
            print("  conda install -c bioconda mafft")
        elif system == "Darwin":
            print("  brew install mafft")
        else:
            print("  sudo apt install mafft muscle")
