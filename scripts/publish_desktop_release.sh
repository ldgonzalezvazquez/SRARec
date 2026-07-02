#!/usr/bin/env bash
# Upload SRARec Desktop zip to GitHub Releases (direct-download link for README).
#
# Prerequisite:
#   brew install gh git-lfs
#   gh auth login
#
# Build the zip on a Mac with the GUI sources, then run from repo root:
#   scripts/publish_desktop_release.sh /path/to/SRARec-Desktop-1.0.0-macos.zip
set -euo pipefail

TAG="${TAG:-desktop-1.0.0}"
ASSET_NAME="${ASSET_NAME:-SRARec-Desktop-1.0.0-macos.zip}"
ZIP_PATH="${1:-}"
GITHUB_REPO="${GITHUB_REPO:-ldgonzalezvazquez/SRARec}"

if [[ -z "$ZIP_PATH" || ! -f "$ZIP_PATH" ]]; then
  echo "Usage: $0 /path/to/${ASSET_NAME}"
  echo ""
  echo "Build the zip with SRARec Desktop packaging (make_release.sh), then upload it here."
  exit 1
fi

if ! command -v gh >/dev/null 2>&1; then
  echo "ERROR: GitHub CLI (gh) is required. Install with: brew install gh"
  exit 1
fi

if ! gh auth status >/dev/null 2>&1; then
  echo "ERROR: Run 'gh auth login' first."
  exit 1
fi

NOTES="SRARec Desktop 1.0.0 for macOS (Apple Silicon and Intel).

Unzip, then double-click Install SRARec Desktop.command.
If macOS blocks the installer: System Settings → Privacy & Security → Open Anyway."

if gh release view "$TAG" -R "$GITHUB_REPO" >/dev/null 2>&1; then
  echo "Updating release $TAG ..."
  gh release upload "$TAG" "$ZIP_PATH#${ASSET_NAME}" --clobber -R "$GITHUB_REPO"
else
  echo "Creating release $TAG ..."
  gh release create "$TAG" "$ZIP_PATH#${ASSET_NAME}" \
    --title "SRARec Desktop 1.0.0 (macOS)" \
    --notes "$NOTES" \
    --latest \
    -R "$GITHUB_REPO"
fi

DOWNLOAD_URL="https://github.com/${GITHUB_REPO}/releases/download/${TAG}/${ASSET_NAME}"
echo ""
echo "Direct download URL (use in README; for pre-releases use tag URL, not /latest/):"
echo "  $DOWNLOAD_URL"
