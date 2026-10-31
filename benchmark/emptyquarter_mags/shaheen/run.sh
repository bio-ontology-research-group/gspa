#!/usr/bin/env bash
# Wrapper around driver.py for printing (default) or submitting (--submit)
# the Shaheen staging and gapsmith batch jobs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/common.sh"

PY="${PYTHON_BIN:-python3}"
if [ ! -x "$PY" ]; then
    PY="$(command -v python3)"
fi

exec "$PY" "$SCRIPT_DIR/driver.py" "$@"
