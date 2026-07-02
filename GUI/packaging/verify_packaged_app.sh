#!/usr/bin/env bash
# Quick smoke test for the packaged SRARec Desktop.app
set -euo pipefail

APP="${1:-$HOME/Desktop/SRARec-desktop-build/dist/SRARec Desktop.app}"
BIN="$APP/Contents/MacOS/SRARec Desktop"

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: Executable not found: $BIN"
  exit 1
fi

echo "Bundle: $APP"
du -sh "$APP"

for rel in \
  "Contents/Frameworks/SRARec1.3.py" \
  "Contents/Frameworks/generate_graphs.py" \
  "Contents/Frameworks/GUI/resources/reference_indexes/sars_cov_2/NC_045512.2.fna" \
  "Contents/Frameworks/GUI/resources/annotations/SC2.gb"
do
  if [[ -f "$APP/$rel" ]]; then
    echo "OK  $rel"
  else
    echo "MISSING $rel"
    exit 1
  fi
done

echo "Launching app (4s) ..."
"$BIN" 2>&1 &
PID=$!
sleep 4
if kill -0 "$PID" 2>/dev/null; then
  echo "OK  App started (pid $PID)"
  kill "$PID" 2>/dev/null || true
  echo "Packaged app verification passed."
else
  echo "ERROR: App exited immediately"
  exit 1
fi
