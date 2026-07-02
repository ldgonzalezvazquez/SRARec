#!/usr/bin/env bash
# Convert GUI/icon/icon.jpeg to icon.icns for macOS app bundle.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
ICON_SRC="$REPO/GUI/icon/icon.jpeg"
ICON_ICNS="$REPO/GUI/icon/icon.icns"
ICONSET="$REPO/GUI/icon/icon.iconset"

if [[ ! -f "$ICON_SRC" ]]; then
  echo "Error: icon source not found: $ICON_SRC" >&2
  exit 1
fi

rm -rf "$ICONSET"
mkdir -p "$ICONSET"

sips -s format png -z 16 16   "$ICON_SRC" --out "$ICONSET/icon_16x16.png"      >/dev/null
sips -s format png -z 32 32   "$ICON_SRC" --out "$ICONSET/icon_16x16@2x.png"   >/dev/null
sips -s format png -z 32 32   "$ICON_SRC" --out "$ICONSET/icon_32x32.png"      >/dev/null
sips -s format png -z 64 64   "$ICON_SRC" --out "$ICONSET/icon_32x32@2x.png"   >/dev/null
sips -s format png -z 128 128 "$ICON_SRC" --out "$ICONSET/icon_128x128.png"    >/dev/null
sips -s format png -z 256 256 "$ICON_SRC" --out "$ICONSET/icon_128x128@2x.png"  >/dev/null
sips -s format png -z 256 256 "$ICON_SRC" --out "$ICONSET/icon_256x256.png"     >/dev/null
sips -s format png -z 512 512 "$ICON_SRC" --out "$ICONSET/icon_256x256@2x.png"  >/dev/null
sips -s format png -z 512 512   "$ICON_SRC" --out "$ICONSET/icon_512x512.png"     >/dev/null
sips -s format png -z 1024 1024 "$ICON_SRC" --out "$ICONSET/icon_512x512@2x.png" >/dev/null

xattr -cr "$ICONSET"
iconutil -c icns "$ICONSET" -o "$ICON_ICNS"
rm -rf "$ICONSET"

echo "Created $ICON_ICNS"
