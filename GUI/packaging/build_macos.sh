#!/usr/bin/env bash
# Build SRARec Desktop .app into ~/Desktop/SRARec-desktop-build
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
VENV="$REPO/GUI/.venv"
BUILD_ROOT="${BUILD_ROOT:-$HOME/Desktop/SRARec-desktop-build}"

if [[ ! -x "$VENV/bin/python" ]]; then
  echo "Creating venv at $VENV ..."
  python3.11 -m venv "$VENV"
  "$VENV/bin/pip" install -r "$REPO/GUI/requirements.txt"
fi

"$VENV/bin/pip" install -q pyinstaller

export BUILD_ROOT
mkdir -p "$BUILD_ROOT"

echo "Generating app icon..."
bash "$REPO/GUI/packaging/make_icon.sh"

echo "Building into $BUILD_ROOT ..."
cd "$REPO"
"$VENV/bin/pyinstaller" \
  --noconfirm \
  --clean \
  --distpath "$BUILD_ROOT/dist" \
  --workpath "$BUILD_ROOT/build" \
  GUI/packaging/SRARecDesktop.spec

echo ""
echo "Done. Application bundle:"
echo "  $BUILD_ROOT/dist/SRARec Desktop.app"
ls -la "$BUILD_ROOT/dist/" || true
