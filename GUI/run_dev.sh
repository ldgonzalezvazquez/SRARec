#!/usr/bin/env bash
# Launch SRARec Desktop in development (uses the project venv).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
VENV="$REPO/GUI/.venv"
PYTHON="$VENV/bin/python"

if [[ ! -x "$PYTHON" ]]; then
  echo "Virtualenv not found. Create it with:" >&2
  echo "  python3.11 -m venv $VENV && $VENV/bin/pip install -r $REPO/GUI/requirements.txt" >&2
  exit 1
fi

PLUGIN="$("$PYTHON" -c 'import PySide6, os; print(os.path.join(os.path.dirname(PySide6.__file__), "Qt", "plugins"))')"
export QT_PLUGIN_PATH="$PLUGIN"

exec "$PYTHON" "$REPO/GUI/app.py" "$@"
