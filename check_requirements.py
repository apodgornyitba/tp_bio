#!/usr/bin/env python3
"""Check local dependencies needed to run the TP pipeline."""

from __future__ import annotations

import importlib.util
import sys

from platform_tools import (
    find_blastp,
    find_emboss_tools,
    find_makeblastdb,
    find_msa_tool,
    print_install_hints,
)


def status(label, value, required=True):
    marker = "OK" if value else ("MISSING" if required else "OPTIONAL")
    detail = value or "no encontrado"
    print(f"[{marker}] {label}: {detail}")
    return bool(value) or not required


def main():
    ok = True

    ok &= status("Python", sys.executable)
    ok &= status("BioPython", "instalado" if importlib.util.find_spec("Bio") else None)

    msa_tool, msa_path = find_msa_tool()
    ok &= status("MSA (MAFFT o MUSCLE)", f"{msa_tool}: {msa_path}" if msa_tool else None)

    emboss_tools = find_emboss_tools()
    ok &= status("EMBOSS getorf", emboss_tools["getorf"])
    ok &= status("EMBOSS prosextract", emboss_tools["prosextract"])
    ok &= status("EMBOSS patmatmotifs", emboss_tools["patmatmotifs"])

    status("BLAST+ blastp (necesario para BLAST local bonus)", find_blastp(), required=False)
    status("BLAST+ makeblastdb (necesario para preparar DB local)", find_makeblastdb(), required=False)

    if ok:
        print("\nListo para correr la entrega final con REQUIRE_EMBOSS=1.")
        return

    print("\nFaltan dependencias requeridas. Instalar segun sistema:")
    print_install_hints("general")
    sys.exit(1)


if __name__ == "__main__":
    main()
