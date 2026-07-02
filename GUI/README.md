# SRARec Desktop (macOS)

PySide6 graphical interface for read-level recombination analysis from **preprocessed BAM** files.

## Download (end users)

Pre-built macOS installers are published on [GitHub Releases](https://github.com/ldgonzalezvazquez/SRARec/releases):

**[SRARec-Desktop-1.0.0.dmg](https://github.com/ldgonzalezvazquez/SRARec/releases/download/desktop-v1.0.0/SRARec-Desktop-1.0.0.dmg)**

See the main [README](../README.md#srarec-desktop-macos-gui) for installation and usage.

## Requirements

- macOS (Apple Silicon or Intel)
- Python 3.11+ recommended

## Development setup

From the repository root:

```bash
python3.11 -m venv GUI/.venv
source GUI/.venv/bin/activate
pip install -r GUI/requirements.txt
python GUI/app.py
```

Or use the helper script (sets Qt plugin path explicitly):

```bash
GUI/run_dev.sh
```

## BAM input requirements

The BAM file must be ready for SRARec analysis:

1. **Coordinate-sorted** (`@HD SO:coordinate` in the header)
2. **Indexed** (`.bai` alongside the BAM, or allow the app to create it)
3. **Deduplicated** (recommended; duplicates may inflate signals)
4. For paired-end data: **mate information** should be correct (`fixmate` already applied)

Mapping, trimming, and SRA download must be done **before** using this GUI (e.g. with the full `SRARec1.3.py` CLI).

## References

This GUI includes bundled references for:

- SARS-CoV-2 (`NC_045512.2`)
- HIV-1 (`NC_001802.1`)

### Custom reference

Index a custom reference with the CLI from the repository root:

```bash
python3 SRARec1.3.py --index_only --custom_ref /path/to/reference.fasta
```

Analysis with a custom reference must be run with the **full CLI**, not this desktop app.

## Output

Results are written to the selected output folder as `results_ravi.txt` (one line per sample, **1-based** coordinates). Interactive HTML profiles and export options are available from the GUI after a successful run, or from an existing results file via **Save results...** / **Regenerate graph**.

## macOS packaging

Build output goes to **`~/Desktop/SRARec-desktop-build/`** (outside the repository).

```bash
# From repository root
chmod +x GUI/packaging/build_macos.sh
./GUI/packaging/build_macos.sh
```

Result:

```text
~/Desktop/SRARec-desktop-build/dist/SRARec Desktop.app
```

Verify the bundle:

```bash
./GUI/packaging/verify_packaged_app.sh
```

Open the app:

```bash
open ~/Desktop/SRARec-desktop-build/dist/SRARec\ Desktop.app
```

Create a distributable DMG for GitHub Releases:

```bash
./GUI/packaging/make_dmg.sh
# → ~/Desktop/SRARec-desktop-build/dist/SRARec-Desktop-1.0.0.dmg
```

Notes:

- Build uses **Python 3.11 + venv + pip** (not Conda) and **PyInstaller**.
- Bundle size is ~800 MB (Qt + scientific stack + reference files).
- **Code signing and notarization** are still required before distributing to other Macs.
