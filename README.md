# SRARec

[![Python](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#srarec-desktop-macos)
[![Version](https://img.shields.io/badge/version-1.3-green.svg)](#citation)

**SRARec** is an open-source Python tool for detecting read-supported recombination-compatible signals from high-throughput sequencing data.

It accepts **NCBI-SRA accessions**, local **FASTQ** files, or preprocessed **BAM** files, and supports **Illumina**, **Oxford Nanopore (ONT)** and **PacBio** data (single-end and paired-end).

Built-in references: **SARS-CoV-2** (`NC_045512.2`) and **HIV-1** (`NC_001802.1`). Custom references are also supported.

> **macOS — analyse BAM files without the terminal**  
> [Download SRARec Desktop](https://github.com/ldgonzalezvazquez/SRARec/releases/download/desktop-1.0.0/SRARec-Desktop-1.0.0-macos.zip) (~300 MB). See [SRARec Desktop](#srarec-desktop-macos).

> **CoViE study**  
> Large-scale analyses from the CoViE project are kept on the `Exclusive_CoViE` branch. This repository documents the tool only.

---

## Table of contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [Example data](#example-data)
- [How to run](#how-to-run)
- [Interpreting results](#interpreting-results)
- [Plotting recombination profiles](#plotting-recombination-profiles)
- [SRARec Desktop (macOS)](#srarec-desktop-macos)
- [Optional HIV workflow (SPAdes + SHIVER)](#optional-hiv-workflow-spades--shiver)
- [Command-line reference](#command-line-reference)
- [Limitations](#limitations)
- [Citation](#citation)
- [References](#references)
- [License](#license)
- [Contact](#contact)

---

## Quick start

```bash
# 1. Clone and enter the repository
git clone https://github.com/ldgonzalezvazquez/SRARec.git
cd SRARec

# 2. Create the environment (recommended)
conda env create -f environment.yml
conda activate srarec

# 3. Run an analysis (SRA accession or local files in examples/ — see below)
python3 SRARec1.3.py ERR13573889 --slow

# 4. Generate an interactive recombination profile
python3 generate_graphs.py results_slow.txt --virus sars-cov-2 --exclusive-positions
# Without gene annotation (no reference .gb required):
python3 generate_graphs.py results_slow.txt --no-annotation --exclusive-positions
```

For HIV-1, add `--virus hiv`.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ldgonzalezvazquez/SRARec.git
cd SRARec
```

### 2. Create the Conda environment

**Recommended** — installs Python and dependencies:

```bash
conda env create -f environment.yml
conda activate srarec
```

**Minimal** — Python only:

```bash
conda create -n srarec python=3.10
conda activate srarec
pip install pysam tqdm biopython pandas numpy matplotlib plotly
```

### 3. External tools *(only if you did not use `environment.yml`)*

Skip this step if you created the environment with `environment.yml` in step 2 — those dependencies are already included.

If you used the **minimal** installation (Python only), install the external programs SRARec calls for download, preprocessing, mapping and BAM handling:

| Tool | Role |
|---|---|
| [SRA Toolkit](https://github.com/ncbi/sra-tools) | Download and convert SRA runs |
| [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) | Retrieve SRA metadata |
| [Cutadapt](https://cutadapt.readthedocs.io/) | Adapter trimming |
| [BWA-MEM](https://github.com/lh3/bwa) | Map Illumina reads |
| [minimap2](https://github.com/lh3/minimap2) | Map ONT / PacBio reads |
| [SAMtools](https://www.htslib.org/) | BAM processing and pileup |
| [BCFtools](https://samtools.github.io/bcftools/) | Optional VCF output for coinfection screening |
| [Picard](https://broadinstitute.github.io/picard/) | Reference dictionary creation |
| [BBMap](https://bbmap.org/) | Optional paired-end read repair (`repair.sh`) |
| [SPAdes](https://github.com/ablab/spades) | Optional HIV assembly (`--rnaviral` mode) |
| [SHIVER](https://github.com/ChrisHIV/shiver) | Optional sample-specific HIV reference |

Install missing tools with Bioconda:

```bash
conda install -c bioconda -c conda-forge \
  sra-tools entrez-direct cutadapt bwa minimap2 \
  samtools bcftools picard bbmap spades shiver
```

### 4. Verify

```bash
python3 SRARec1.3.py --help
python3 generate_graphs.py --help
```

---

## Example data

The `examples/` directory holds **local test files**. Download them from GitHub:

**[github.com/ldgonzalezvazquez/SRARec/tree/main/examples](https://github.com/ldgonzalezvazquez/SRARec/tree/main/examples)**

| File | Use |
|---|---|
| `examples/example_1.fastq` | Forward reads |
| `examples/example_2.fastq` | Reverse reads |
| `examples/example.bam` | Pre-mapped BAM (also for [SRARec Desktop](#srarec-desktop-macos)) |

See [`examples/README.md`](examples/README.md) for commands and the BAM naming convention (`{folder_name}.bam` inside the sample folder).

**Paired-end FASTQ:**

```bash
python3 SRARec1.3.py --fastq examples/example_1.fastq examples/example_2.fastq example_run \
  --slow --platform ILLUMINA
```

**BAM** (if the BAM is `examples/example.bam`, use folder `example` with `example/example.bam`, or rename to `examples/examples.bam` — details in `examples/README.md`):

```bash
python3 SRARec1.3.py example --bam
```

**SRA accessions** — you can also test without local files. Example runs:

`ERR13573889` · `ERR13020098` · `ERR11859871` · `ERR10430152` · `ERR10304387`

```bash
python3 SRARec1.3.py ERR13573889 --slow
```

Platform is auto-detected from SRA metadata for SRA accessions.

### Repository layout

```text
SRARec/
├── SRARec1.3.py              # Main analysis script
├── generate_graphs.py        # Recombination profile plots
├── environment.yml           # Conda environment
├── examples/                 # Local test data (see examples/README.md)
│   ├── README.md
│   ├── example_1.fastq       # download from GitHub examples/
│   ├── example_2.fastq       # download from GitHub examples/
│   └── example.bam           # download from GitHub examples/
├── reference_indexes/
│   ├── sars_cov_2/           # SARS-CoV-2 FASTA, indexes and NC_045512.2.gb
│   ├── hiv_1/                # HIV-1 FASTA, indexes, NC_001802.1.gb and K03455.1_MOD.gb
│   └── <custom_name>/        # Custom references indexed with --index_only
└── complementary_data/
    └── adapters.fa           # Adapter sequences
```

---

## How to run

Every run requires four choices:

1. **Reference genome** — SARS-CoV-2 (default), HIV-1 (`--virus hiv`), or a custom reference.
2. **Input type** — SRA accession, local FASTQ (`--fastq`), or existing BAM (`--bam`).
3. **Pipeline mode** — slow (`--slow`, recommended) or fast (default).
4. **Search thresholds** — minimum read support and gamete frequency filters used during the recombination search.

### Reference genomes

| Reference | Command flag | Files |
|---|---|---|
| SARS-CoV-2 (default) | *(none)* | `reference_indexes/sars_cov_2/NC_045512.2.fna` |
| HIV-1 | `--virus hiv` | `reference_indexes/hiv_1/NC_001802.1.fna` |
| Custom | `--index_only --custom_ref --virus <name>` then `--virus <name>` | `reference_indexes/<name>/` |

Custom reference — index first, then analyse. The `--virus` name must match the subfolder created under `reference_indexes/`:

```bash
python3 SRARec1.3.py --index_only --custom_ref path/to/reference.fasta --virus my_virus
python3 SRARec1.3.py SRRXXXXXXXX --slow --virus my_virus
```

This creates `reference_indexes/my_virus/reference.fna` and the mapping indexes alongside it.

### Input types

#### SRA accession

```bash
python3 SRARec1.3.py SRRXXXXXXXX --slow
python3 SRARec1.3.py ERRXXXXXXXX --slow --virus hiv
```

Accepted prefixes: `SRR`, `ERR`, `DRR`. The sequencing platform is detected automatically from SRA RunInfo — `--platform` is not required.

#### Local FASTQ (`--fastq`)

Single-end:

```bash
python3 SRARec1.3.py --fastq reads.fastq output_folder \
  --slow --platform ILLUMINA
```

Paired-end:

```bash
python3 SRARec1.3.py --fastq reads_R1.fastq reads_R2.fastq output_folder \
  --slow --platform ILLUMINA --virus hiv
```

Platform values: `ILLUMINA`, `OXFORD_NANOPORE`, `ONT`, `PACBIO`, `PACBIO_SMRT`.

#### Precomputed BAM (`--bam`)

The BAM must be named `{folder_name}.bam` inside the sample folder:

```bash
python3 SRARec1.3.py sample_folder --bam
```

`--platform` is optional in BAM mode: it does not affect the analysis, only the `Platform:` label written to the results file (default: `ILLUMINA`).

Example layout:

```text
sample_folder/
└── {folder_name}.bam
```

### Pipeline modes

| Mode | Flag | When to use |
|---|---|---|
| **Slow** | `--slow` | Standard analyses; trims reads before mapping (recommended) |
| **Fast** | *(default)* | Faster but less conservative — avoid when possible |

### Search thresholds

These parameters control the recombination **search** (which polymorphic sites and gametes are considered). They are passed as comma-separated lists; SRARec tests every combination.

| Argument | Meaning | Default |
|---|---|---:|
| `--MIN_READS_LIST` | Minimum reads supporting a variant at a polymorphic site | `10` |
| `--READS_GAMETOS_MINORITARIOS_LIST` | Minimum reads supporting each minor gamete | `5` |
| `--THRESHOLD_GAMETO_FREQ_LIST` | Minimum frequency for each gamete | `0.05` |

Override defaults when needed, for example:

```bash
python3 SRARec1.3.py SRRXXXXXXXX --slow \
  --MIN_READS_LIST 5 \
  --READS_GAMETOS_MINORITARIOS_LIST 1 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.00
```

To ignore indels, add `--noindel` — recommended for samples with high diversity, where indels can add noise to the search.

### Common examples

```bash
# SARS-CoV-2 from SRA (default search thresholds)
python3 SRARec1.3.py ERR13573889 --slow

# HIV-1 from SRA
python3 SRARec1.3.py ERR13020098 --slow --virus hiv

# Local paired-end FASTQ (files in examples/)
python3 SRARec1.3.py --fastq examples/example_1.fastq examples/example_2.fastq example_run \
  --slow --platform ILLUMINA

# Local ONT reads (--platform required for mapping)
python3 SRARec1.3.py --fastq reads.fastq my_sample --slow --platform ONT

# BAM-based analysis (--platform optional; only written to the Platform field)
python3 SRARec1.3.py example --bam
```

Other example SRA accessions: `ERR11859871`, `ERR10430152`, `ERR10304387`.

---

## Interpreting results

SRARec writes `results_<mode>.txt` in the working directory (`results_slow.txt`, `results_fast.txt`, or `results_bam.txt`). Each line is one sample.

### Result fields

```text
Sra:<sample_id>
Coinf_detected:<True/False>
Min_read_num:<value>
Min_gameto_num:<value>
Min_gameto_freq:<value>
Platform:<platform>
Rec:<recombinant-compatible_pairs>
No_Rec:<non_recombinant-compatible_pairs>
BegFin_Reads:<read_start_end_coordinates>
```

| Field | Meaning |
|---|---|
| `Sra` | Sample identifier |
| `Coinf_detected` | Preliminary screening flag — see below |
| `Min_read_num` / `Min_gameto_num` / `Min_gameto_freq` | Search thresholds used for this line |
| `Platform` | Sequencing platform |
| `Rec` | Pairs classified as recombination-compatible (see example below) |
| `No_Rec` | Pairs that do not meet recombination-compatible criteria |
| `BegFin_Reads` | Read start/end coordinates per read (used for coverage plots) |

> **Coordinates** — From v1.3 onward, positions in `Rec`, `No_Rec` and `BegFin_Reads` are **1-based** (first nucleotide = 1).

### Reading `Rec`, `No_Rec` and `BegFin_Reads` — example

Consider this result line (fields shortened for clarity):

```text
Rec:22877,23020,2,G.T,T.C,mireads:10,minfreq:0.05;24100,24250,1,C.A,G.T,mireads:10,minfreq:0.05;25000,25100,3,A.G,C.T,mireads:10,minfreq:0.05
No_Rec:22877,22880,2
BegFin_Reads:1,299;300,598;599,897
```

**`Rec`** — semicolon-separated recombination-compatible pairs. Each entry has four fields:

```text
position_A, position_B, recovery_count, gamete_detail
```

- `22877,23020,2,G.T,T.C,mireads:10,minfreq:0.05` — positions 22877 and 23020; `recovery_count` is `2` because this pair is supported within two nested compatible ranges in the hierarchical search; gamete pattern `G.T` / `T.C`.
- `24100,24250,1,C.A,G.T,mireads:10,minfreq:0.05` — a second compatible pair with `recovery_count` `1` (supported in one hierarchical context).
- `25000,25100,3,A.G,C.T,mireads:10,minfreq:0.05` — a third compatible pair with `recovery_count` `3` (supported within three nested compatible ranges).

**`No_Rec`** — pairs tested but not compatible. Format: `position_A, position_B, detail`.

- `22877,22880,2` — positions 22877 and 22880; detail `2` indicates partial gamete support (not all four gametes passed the search thresholds).

**`BegFin_Reads`** — start and end of each read used in the analysis, as `start,end` pairs separated by semicolons:

- `1,299` — one read spans positions 1–299
- `300,598` — another read spans 300–598
- `599,897` — another read spans 599–897

These coordinates feed the coverage track in `generate_graphs.py`.

### `Coinf_detected` — when is it `True`?

`Coinf_detected:True` is set when, at **any polymorphic position**, at least **two alleles** each reach a frequency **> 10%** in the read pileup. It is a **preliminary screening flag** used to trigger optional VCF generation.

`Coinf_detected:False` is written when fewer than two polymorphic positions pass the minimum-read filter, so no search is performed.

> **Important** — `Coinf_detected:True` is **not** evidence of coinfection. Do **not** report coinfections in publications based on this flag alone. Always inspect the optional VCF manually, review allele frequencies and read support, and corroborate with independent biological or epidemiological evidence.

### Optional VCF output

When `Coinf_detected:True`, SRARec may write a screening VCF to:

```text
vcf_coinfecciones/<sample_id>.vcf
```

Treat this file as a starting point for manual review, not as a confirmed coinfection call.

### What the results mean

SRARec reports **read-supported recombination-compatible signals**: patterns in the data that are inconsistent with a single non-recombining haplotype, given the read length and coverage available.

These signals are **not** direct proof of recombination. Always interpret them together with sequencing depth, read length, mapping quality, reference choice and biological context.

We recommend analysing results in **context across multiple samples** and generating **genome-wide profiles** with `generate_graphs.py` (using `--exclusive-positions` and, when appropriate, `--aggregate-by-patient`). Aggregating signals across samples helps reveal general patterns and reduces the impact of sequencing artefacts that may appear in individual runs.

To assess a **single candidate event**, inspect the search output carefully: supporting read counts, minor gamete frequency, whether any involved polymorphism is an indel, and `recovery_count` in `Rec` (always ≥ 1 for compatible pairs).

---

## Plotting recombination profiles

Use `generate_graphs.py` to turn a results file into an interactive HTML profile, a PNG figure, and per-position tables.

Bundled GenBank annotations live under `reference_indexes/`:

| `--virus` | Annotation file |
|---|---|
| `sars-cov-2` (default) | `reference_indexes/sars_cov_2/NC_045512.2.gb` |
| `hiv` | `reference_indexes/hiv_1/NC_001802.1.gb` |
| `hxb2` | `reference_indexes/hiv_1/K03455.1_MOD.gb` |

The default HIV-1 mapping reference is **NC_001802.1** (`--virus hiv`). Use **`--virus hxb2`** when positions in your results file are on **HXB2** (K03455) coordinates — for example, after analysing a BAM aligned to HXB2, or after adapting SHIVER output coordinates using SHIVER's correspondence files (see [SHIVER workflow](#optional-hiv-workflow-spades--shiver)).

You can also plot **without any annotation** using **`--no-annotation`**: the recombination and coverage profiles are still drawn, but no gene track or gene legend is shown. This is useful when coordinates do not match a bundled reference (e.g. sample-specific SHIVER references) or when you only need the ρ profile.

### Examples

Always include `--exclusive-positions` to prevent true and false signals from overlapping at the same genomic position.

```bash
# SARS-CoV-2 profile (default annotation)
python3 generate_graphs.py results_slow.txt --virus sars-cov-2 --exclusive-positions

# HIV-1 profile (NC_001802.1 — default built-in reference)
python3 generate_graphs.py results_slow.txt --virus hiv --exclusive-positions

# HIV-1 profile when results are on HXB2 (K03455) coordinates
python3 generate_graphs.py results_slow.txt --virus hxb2 --exclusive-positions

# Profile without gene annotation (no gene track or legend)
python3 generate_graphs.py results_slow.txt --no-annotation --exclusive-positions

# Custom output path
python3 generate_graphs.py results_bam.txt -o output/my_profile.html \
  --virus sars-cov-2 --exclusive-positions

# Explicit annotation file
python3 generate_graphs.py results_slow.txt --exclusive-positions \
  --annotation reference_indexes/hiv_1/NC_001802.1.gb
```

### Output files

| File | Content |
|---|---|
| `*_recomb_profile.html` | Interactive Plotly figure |
| `*_recomb_profile.png` | Static image |
| `*_recombination_values.txt` | Recombination ratio ρ per position |
| `*_coverage_values.txt` | Read coverage from `BegFin_Reads` |

### Main arguments

| Argument | Description |
|---|---|
| `results_txt` | SRARec results file |
| `-o`, `--output` | Output HTML path |
| `--virus` | Built-in annotation selector (`sars-cov-2`, `hiv`, `hxb2`); ignored with `--no-annotation` |
| `--annotation` | Custom GenBank or GFF3 file |
| `--no-annotation` | Skip gene track and gene legend (recombination/coverage only) |
| `--alpha` | Scaling factor for ρ (default: `1.0`) |
| `--exclusive-positions` | **Recommended.** Avoid true/false signal overlap at the same position |
| `--aggregate-by-patient` | Aggregate events by sample before summing |
| `--no-coverage` | Skip the coverage curve |
| `--no-begfin-concentration` | Skip read start/end concentration track |

---

## SRARec Desktop (macOS)

**SRARec Desktop** is a standalone macOS app for BAM-based recombination analysis. It wraps `SRARec1.3.py` and `generate_graphs.py` in a graphical interface.

### Download

**[SRARec Desktop 1.0.0 (macOS zip)](https://github.com/ldgonzalezvazquez/SRARec/releases/download/desktop-1.0.0/SRARec-Desktop-1.0.0-macos.zip)** (~300 MB)

### Install

1. Download and unzip the archive.
2. Double-click `Install SRARec Desktop.command`.
3. If macOS blocks the installer: **System Settings → Privacy & Security → Open Anyway**.
4. The app is installed in **Applications**.

### Quick start

1. Launch **SRARec Desktop**.
2. Click **Add BAMs...** and select one or more coordinate-sorted, indexed BAM files.
3. Choose platform, reference (SARS-CoV-2 or HIV-1) and output folder.
4. Click **Run analysis** — results are written to `results_bam.txt`.
5. Use **Regenerate graph** or **Save graph** to export HTML/PNG profiles.

Mapping, trimming and SRA download must be done beforehand (e.g. with `SRARec1.3.py`).

| Requirement | Detail |
|---|---|
| Platform | macOS (Apple Silicon or Intel) |
| Input | Coordinate-sorted BAM; `.bai` index optional (the app can create it) |
| Disk space | ~800 MB installed |

For installation issues or help incorporating new references: **luisdaniel.gonzalez@uvigo.es** or **luisde365@gmail.com**.

---

## Optional HIV workflow (SPAdes + SHIVER)

For HIV samples where a sample-specific reference is preferred, SRARec can run:

```text
trimming → SPAdes (--rnaviral) → SHIVER → BAM analysis
```

Activate with `--use_shiver_assembly`. Adapter sequences are in `complementary_data/adapters.fa`.

```bash
python3 SRARec1.3.py --fastq --use_shiver_assembly \
  reads_R1.fastq reads_R2.fastq hiv_sample \
  --virus hiv --platform ILLUMINA \
  --spades spades.py \
  --shiver_dir "$CONDA_PREFIX" \
  --shiver_ref_alignment MyRefAlign.fasta \
  --shiver_adapters complementary_data/adapters.fa \
  --shiver_primers complementary_data/adapters.fa
```

**Coordinates and plotting.** SHIVER builds a **sample-specific reference** (`{sample}_ref.fasta`) and SRARec reports recombination positions in `results_*.txt` relative to that reference.

By default, [SHIVER](https://github.com/ChrisHIV/shiver) also writes a **correspondence between the sample-specific reference and HXB2** for coordinate adaptation. The main files are:

| SHIVER output | Role |
|---|---|
| `{sample}_BaseFreqs_WithHXB2.csv` | Base frequencies with an HXB2 coordinate column (pairwise alignment of the mapping reference to HXB2) |
| `{sample}_coords.csv` | Dictionary mapping positions on the sample reference to the global reference alignment |
| `{sample}_BaseFreqs_ForGlobalAln.csv` | Base frequencies with global-alignment coordinates |

These files are produced when `GiveHXB2coords=true` in SHIVER's `config.sh` (the default for HIV). See the [SHIVER repository](https://github.com/ChrisHIV/shiver) and the [SHIVER manual](https://github.com/ChrisHIV/shiver/blob/master/docs/ShiverManual.pdf) for full details on coordinate translation and downstream use.

**SRARec does not automatically remap** `Rec`/`No_Rec`/`BegFin_Reads` in `results_*.txt` to HXB2. When coinfection is detected, SRARec only moves `{sample}_BaseFreqs_WithHXB2.csv` to `csv_with_HXB2/` (along with `{sample}_ref.fasta` to `ref_fasta/`). To plot with HXB2 gene annotations, adapt recombination coordinates using SHIVER's correspondence first, then pass **`--virus hxb2`** to `generate_graphs.py`.

To plot SHIVER output without coordinate adaptation:

- **`--no-annotation`** — profile without genes (recommended if positions are still on the sample-specific reference).

See [Plotting recombination profiles](#plotting-recombination-profiles).

---

## Command-line reference

```text
python3 SRARec1.3.py [options] [--slow] <SRR>
python3 SRARec1.3.py --fastq [--slow] <fastq1> [fastq2] <output_folder> [options]
python3 SRARec1.3.py [--slow] <sample_folder> --bam [options]
```

| Argument | Description |
|---|---|
| `<SRR>` | NCBI-SRA accession |
| `--slow` | Slow pipeline (recommended) |
| `--fastq` | Local FASTQ input |
| `--bam` | Precomputed BAM input |
| `--virus TEXT` | Reference: default SARS-CoV-2; `hiv` for HIV-1; any other value uses `reference_indexes/<name>/` |
| `--custom_ref PATH` | Custom FASTA (with `--index_only` and `--virus`) |
| `--index_only` | Build custom reference index in `reference_indexes/<virus>/` and exit |
| `--platform TEXT` | Platform label for results. Auto-detected from SRA; required for `--fastq` mapping; optional metadata in `--bam` mode (default: `ILLUMINA`) |
| `--noindel` | Ignore indels (recommended for high-diversity samples) |
| `--verbose` | Verbose output |
| `--MIN_READS_LIST` | Minimum reads per variant (default: `10`) |
| `--READS_GAMETOS_MINORITARIOS_LIST` | Minimum reads per gamete (default: `5`) |
| `--THRESHOLD_GAMETO_FREQ_LIST` | Minimum gamete frequency (default: `0.05`) |
| `--use_shiver_assembly` | SPAdes + SHIVER HIV workflow |
| `--spades`, `--spades_threads`, `--spades_mem_gb` | SPAdes options |
| `--shiver_dir`, `--shiver_threads`, `--shiver_ref_alignment`, `--shiver_adapters`, `--shiver_primers` | SHIVER options |

---

## Limitations

- SRARec detects read-supported recombination-compatible signals, not confirmed recombination events.
- Results depend on read length, depth, mapping quality, reference genome, amplicon design, PCR artefacts and search thresholds.
- Prefer multi-sample profiles over isolated per-sample patterns to reduce sequencing noise.
- For individual events, review read support, minor gamete frequency, indel involvement and recovery count before drawing conclusions.
- `Coinf_detected:True` is a screening flag only — not evidence of coinfection.

---

## Citation

If you use SRARec, please cite:

```text
González-Vázquez LD, Iglesias-Rivas P, Arenas M, Martin DP.
SRARec: A program for detecting recombination in sequencing reads and its application to uncover recombination patterns in SARS-CoV-2 and HIV-1
[Journal / preprint / DOI to be added]
```

---

## References

Methodological and software references used by SRARec:

1. Hudson RR, Kaplan NL. Statistical properties of the number of recombination events in the history of a sample of DNA sequences. *Genetics*. 1985;111(1):147–164. doi [https://doi.org/](https://doi.org/10.1093/genetics/111.1.147)

2. Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. *Bioinformatics*. 2009;25(14):1754–1760. doi:[10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324)

3. Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*. 2013;1303.3997. doi:[10.48550/arXiv.1303.3997](
https://doi.org/10.48550/arXiv.1303.3997)

4. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. 2018;34(18):3094–3100. doi:[10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

5. Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience*. 2021;10(2):giab008. doi:[10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)

6. Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*. 2011;17(1):10–12. doi:[10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)

7. Bankevich A, Nurk S, Antipov D, et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. *J Comput Biol*. 2012;19(5):455–477. doi:[10.1089/cmb.2012.0021](https://doi.org/10.1089/cmb.2012.0021)

8. Wymant C, Blanquart F, Golubchik T, et al. Easy and accurate reconstruction of whole HIV genomes from short-read sequence data with shiver. *Virus Evol*. 2018;4(1):vey007. doi:[10.1093/ve/vey007](https://doi.org/10.1093/ve/vey007)

9. Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The sequence read archive. *Nucleic Acids Res*. 2011;39(Database issue):D19–D21. doi:[10.1093/nar/gkq1019](https://doi.org/10.1093/nar/gkq1019)

10. Katz K, Shutov O, Lapoint R, et al. The Sequence Read Archive: a decade more of explosive growth. *Nucleic Acids Res*. 2022;50(D1):D386–D390. doi:[10.1093/nar/gkab1053](https://doi.org/10.1093/nar/gkab1053)

11. Bushnell B, Rood J, Singer E. BBMerge – Accurate paired shotgun read merging via overlap. *PLoS ONE*. 2017;12(10):e0185056. doi:[10.1371/journal.pone.0185056](https://doi.org/10.1371/journal.pone.0185056)

12. Kans J. Entrez Direct: E-utilities on the Unix command line. In: *Entrez Programming Utilities Help*. Bethesda (MD): National Center for Biotechnology Information (US); 2013.

External resources:

- [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [Cutadapt](https://cutadapt.readthedocs.io/)
- [BWA](https://github.com/lh3/bwa)
- [minimap2](https://github.com/lh3/minimap2)
- [SAMtools / HTSlib](https://www.htslib.org/)
- [BCFtools](https://samtools.github.io/bcftools/)
- [Picard](https://broadinstitute.github.io/picard/)
- [SPAdes](https://github.com/ablab/spades)
- [SHIVER](https://github.com/ChrisHIV/shiver)
- [BBMap](https://bbmap.org/)

---

## License

SRARec is distributed under the [MIT License](LICENSE).

---

## Contact

For questions, bug reports or feature requests, open an issue in the repository.

**Maintainers**

- Luis Daniel González-Vázquez — luisdaniel.gonzalez@uvigo.es
- Paula Iglesias-Rivas – paula.iglesias.rivas@uvigo.gal
