# Example input files

Place the following files in this directory to try SRARec locally without downloading from SRA:

| File | Purpose |
|---|---|
| `example_1.fastq` | Forward reads (paired-end FASTQ mode) |
| `example_2.fastq` | Reverse reads (paired-end FASTQ mode) |
| `example.bam` | Pre-mapped BAM (BAM mode) |

Download these files directly from the GitHub repository:

**[github.com/ldgonzalezvazquez/SRARec/tree/main/examples](https://github.com/ldgonzalezvazquez/SRARec/tree/main/examples)**

Open each file in that folder and use **Download** (or **Download raw file**). Place the files in this directory before running the examples below.

## BAM naming convention

In `--bam` mode, SRARec expects the BAM to be named `{folder_name}.bam` inside the sample folder. Two options:

```text
# Option A — folder "example", file example/example.bam
example/
└── example.bam

# Option B — folder "examples", file examples/examples.bam
examples/
└── examples.bam
```

If your file is `examples/example.bam`, rename it to `examples/examples.bam` or move it into an `example/` subfolder as shown above.

## Example commands

**Paired-end FASTQ** (slow mode, SARS-CoV-2):

```bash
python3 SRARec1.3.py --fastq examples/example_1.fastq examples/example_2.fastq example_run \
  --slow --platform ILLUMINA
```

**Precomputed BAM** (folder `example` with `example/example.bam`):

```bash
python3 SRARec1.3.py example --bam
```

**SRARec Desktop (macOS GUI)** — download [`example.bam`](https://github.com/ldgonzalezvazquez/SRARec/tree/main/examples) from the link above, then:

1. Install [SRARec Desktop](https://github.com/ldgonzalezvazquez/SRARec/releases/download/desktop-1.0.0/SRARec-Desktop-1.0.0-macos.zip) (see the main [README](../README.md#srarec-desktop-macos)).
2. Launch **SRARec Desktop**.
3. Click **Add BAMs...** and select `example.bam`.
4. Choose platform and reference (SARS-CoV-2 or HIV-1), set an output folder, and click **Run analysis**.
5. When finished, use **Regenerate graph** or **Save graph** to export the recombination profile.

The BAM must be coordinate-sorted; a `.bai` index is optional (the app can create it).

**Plots** (after a run produces `results_slow.txt` or `results_ravi.txt`):

```bash
python3 generate_graphs.py results_slow.txt --virus sars-cov-2 --exclusive-positions
python3 generate_graphs.py results_slow.txt --no-annotation --exclusive-positions
```

## SRA accessions for testing

You can also run directly from NCBI-SRA without local files. Example accessions that work with the default pipeline:

| Accession | Notes |
|---|---|
| `ERR13573889` | Example SRA run |
| `ERR13020098` | Example SRA run |
| `ERR11859871` | Example SRA run |
| `ERR10430152` | Example SRA run |
| `ERR10304387` | Example SRA run |

```bash
python3 SRARec1.3.py ERR13573889 --slow
```

Add `--virus hiv` for HIV-1 analyses. Platform is auto-detected from SRA metadata.
