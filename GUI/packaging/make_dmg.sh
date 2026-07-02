#!/usr/bin/env bash
# Create a compressed .dmg from the packaged SRARec Desktop.app
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_ROOT="${BUILD_ROOT:-$HOME/Desktop/SRARec-desktop-build}"
APP_NAME="SRARec Desktop"
APP_PATH="$BUILD_ROOT/dist/${APP_NAME}.app"

if [[ ! -d "$APP_PATH" ]]; then
  echo "ERROR: App bundle not found: $APP_PATH"
  echo "Run GUI/packaging/build_macos.sh first."
  exit 1
fi

VERSION="$(python3 -c "import sys; sys.path.insert(0, '$REPO'); from GUI.config import APP_VERSION; print(APP_VERSION)")"
DMG_NAME="SRARec-Desktop-${VERSION}.dmg"
DMG_PATH="$BUILD_ROOT/dist/$DMG_NAME"
STAGING="$(mktemp -d)"

cleanup() {
  rm -rf "$STAGING"
}
trap cleanup EXIT

echo "Staging DMG contents..."
cp -R "$APP_PATH" "$STAGING/"
ln -s /Applications "$STAGING/Applications"

rm -f "$DMG_PATH"
echo "Creating compressed DMG: $DMG_PATH"
hdiutil create \
  -volname "$APP_NAME" \
  -srcfolder "$STAGING" \
  -ov \
  -format UDZO \
  "$DMG_PATH"

echo ""
echo "Done."
ls -lh "$DMG_PATH"
echo ""
echo "Upload to GitHub Releases:"
echo "  $DMG_PATH"
