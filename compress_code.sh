#!/usr/bin/env bash
# Create a code-only archive for running the TP Bio pipeline.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
ARCHIVE="${1:-tp_bio_code.tar.gz}"

if [[ "$ARCHIVE" != /* ]]; then
    ARCHIVE="$ROOT/$ARCHIVE"
fi

INCLUDE_FILES=(
    "primer_config.json"
    "setup.sh"
    "setup.bat"
    "run_pipeline.sh"
    "run_pipeline.bat"
    "run_pipeline.py"
    "check_requirements.py"
    "platform_tools.py"
    "fetch_data.py"
    "prepare_blast_db.py"
    "blast_common.py"
    "Ex1.py"
    "Ex2_a.py"
    "Ex2_local.py"
    "Ex3.py"
    "Ex4.py"
    "Ex5.py"
)

FILES_TO_ARCHIVE=()
for file in "${INCLUDE_FILES[@]}"; do
    if [[ "$file" == "$SCRIPT_NAME" || "$file" == "./$SCRIPT_NAME" ]]; then
        continue
    fi
    if [[ ! -f "$ROOT/$file" ]]; then
        echo "[ERROR] Required file missing: $file" >&2
        exit 1
    fi
    FILES_TO_ARCHIVE+=("$file")
done

mkdir -p "$(dirname "$ARCHIVE")"
tar -czf "$ARCHIVE" -C "$ROOT" "${FILES_TO_ARCHIVE[@]}"

echo "[OK] Code archive created: $ARCHIVE"
echo "[INFO] Files included: ${#FILES_TO_ARCHIVE[@]}"
