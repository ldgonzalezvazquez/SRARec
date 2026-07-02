# SRARec

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#srarec-desktop-macos-gui)
[![Version](https://img.shields.io/badge/version-1.3-green.svg)](#citation)

**SRARec** is an open-source Python tool for detecting read-supported recombination-compatible signals from high-throughput sequencing data. It implements a read-level version of the **four-gamete test**, allowing recombination signal to be evaluated from physically linked variants observed within sequencing reads.

SRARec can analyse data from:

- **NCBI Sequence Read Archive (SRA)** accessions;
- local **FASTQ** files;
- precomputed **BAM** files;
- optional HIV-focused sample-specific assembly workflows using **rnaviralSPAdes** and **SHIVER**.

It supports **Illumina**, **Oxford Nanopore Technologies (ONT)** and **PacBio** data, including both **single-end** and **paired-end** reads.

> **Project note**  
> For analyses associated with the **CoViE study**, use the `Exclusive_CoViE` branch.

> **macOS users — BAM analysis without the terminal**  
> **[Download SRARec Desktop](https://github.com/ldgonzalezvazquez/SRARec/releases/latest/download/SRARec-Desktop-1.0.0-macos.zip)** (~300 MB). Unzip, then double-click `Install SRARec Desktop.command`. See [SRARec Desktop](#srarec-desktop-macos-gui).

---

## Table of contents

- [Overview](#overview)
- [How SRARec works](#how-srarec-works)
- [How to install](#how-to-install)
  - [1. Clone the repository](#1-clone-the-repository)
  - [2. Create the Conda environment](#2-create-the-conda-environment)
  - [3. Install external dependencies](#3-install-external-dependencies)
  - [4. Check the installation](#4-check-the-installation)
- [How to run SRARec](#how-to-run-srarec)
  - [Step 1. Choose the reference genome](#step-1-choose-the-reference-genome)
  - [Step 2. Choose the input type](#step-2-choose-the-input-type)
  - [Step 3. Choose the analysis mode](#step-3-choose-the-analysis-mode)
  - [Step 4. Choose detection thresholds](#step-4-choose-detection-thresholds)
- [Usage examples](#usage-examples)
- [Command-line reference](#command-line-reference)
- [Optional HIV assembly-based workflow](#optional-hiv-assembly-based-workflow)
- [Outputs](#outputs)
- [Plotting recombination profiles](#plotting-recombination-profiles)
- [SRARec Desktop (macOS GUI)](#srarec-desktop-macos-gui)
- [Method summary](#method-summary)
- [Validation](#validation)
- [Large-scale analyses](#large-scale-analyses)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)
- [References and external tools](#references-and-external-tools)
- [Contact](#contact)

---

## Overview

SRARec was developed to detect recombination-compatible patterns directly from sequencing reads. Instead of relying only on consensus genomes, the program examines how variants are linked within individual reads and tests whether pairs of polymorphic positions show all four possible allele combinations.

This makes SRARec useful for viral sequencing datasets where recombination may be present within mixed infections, intra-host diversity, heterogeneous samples or public sequencing datasets.

SRARec reports **read-supported recombination-compatible signals**. These signals should be interpreted together with sequencing depth, read length, mapping quality, genetic diversity, amplicon design and biological context.

---

## How SRARec works

```text
Input data
  ├── NCBI-SRA accession
  ├── Local FASTQ files
  └── Local BAM file
        │
        ▼
Read preprocessing
  ├── SRA download and FASTQ conversion, if required
  ├── paired-end read repair, if required
  └── adapter trimming
        │
        ▼
Mapping to reference genome
  ├── BWA-MEM for Illumina
  └── minimap2 for ONT / PacBio
        │
        ▼
BAM processing and filtering
  ├── sorting and indexing
  ├── mate fixing for paired-end data
  ├── duplicate removal
  ├── mapping-quality filtering
  └── aligned-read-fraction filtering
        │
        ▼
Candidate polymorphic-site detection
  └── positions with at least two variants and enough read support
        │
        ▼
Read-level variant extraction
  └── linked substitutions and/or indels across polymorphic sites
        │
        ▼
Four-gamete test
  ├── recombination-compatible pairs
  └── non-recombinant-compatible pairs
        │
        ▼
Outputs
  ├── results_<mode>.txt
  ├── recombinant-compatible events
  ├── non-recombinant pair log
  ├── read start/end information
  ├── optional VCFs for potential coinfections
  └── genome-wide recombination profiles
```

---

## How to install

### 1. Clone the repository

```bash
git clone <SRARec-repository-url>
cd SRARec
```

### 2. Create the Conda environment

Create and activate a dedicated environment:

```bash
conda create -n srarec python=3.10
conda activate srarec
```


The environment can alternatively be created with the 'environment.yml' file:

```bash
conda env create -f environment.yml
conda activate srarec
```

### 3. Install external dependencies

SRARec uses external command-line tools for SRA download, read preprocessing, mapping, BAM processing, VCF generation and optional assembly-based workflows.

A basic Conda/Bioconda installation can be started with:

```bash
conda install -c bioconda -c conda-forge \
  sra-tools \
  entrez-direct \
  cutadapt \
  bwa \
  minimap2 \
  samtools \
  bcftools \
  picard \
  spades
```

The main dependencies are:

| Tool | Purpose | Official link |
|---|---|---|
| **SRA Toolkit** | Download and convert SRA runs with `prefetch`, `fasterq-dump` and `fastq-dump` | [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools) |
| **Entrez Direct** | Retrieve SRA metadata with `esearch` and `efetch` | [NCBI Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) |
| **Cutadapt** | Adapter trimming | [Cutadapt documentation](https://cutadapt.readthedocs.io/) |
| **BBMap / BBTools** | Optional repair of orphan paired-end reads with `repair.sh` | [BBMap / BBTools](https://bbmap.org/) |
| **BWA-MEM** | Mapping Illumina reads | [BWA GitHub repository](https://github.com/lh3/bwa) |
| **minimap2** | Mapping ONT and PacBio reads | [minimap2 GitHub repository](https://github.com/lh3/minimap2) |
| **SAMtools** | BAM/SAM processing, sorting, indexing and pileup generation | [SAMtools / HTSlib](https://www.htslib.org/) |
| **BCFtools** | Optional VCF generation for potential coinfection analyses | [BCFtools](https://samtools.github.io/bcftools/) |
| **Picard** | Reference dictionary creation in selected workflows | [Picard](https://broadinstitute.github.io/picard/) |
| **SPAdes** | Optional viral assembly before SHIVER-based analysis | [SPAdes GitHub repository](https://github.com/ablab/spades) |
| **SHIVER** | Optional sample-specific HIV reference reconstruction and read mapping | [SHIVER GitHub repository](https://github.com/ChrisHIV/shiver) |
| **NGSNGS** | Optional read simulation for validation analyses | [NGSNGS GitHub repository](https://github.com/RAHenriksen/NGSNGS) |

Some tools, especially **BBMap/BBTools**, **SHIVER**, **NGSNGS**, may need to be installed manually depending on the system.

### 4. Check the installation

Check that SRARec is available:

```bash
python3 SRARec1.3.py --help
```

The previous release `SRARec1.2.py` is kept under `Releases/` for reference.

Check the main external tools:

```bash
prefetch --help
fasterq-dump --help
esearch -help
efetch -help
cutadapt --version
bwa
minimap2 --version
samtools --version
bcftools --version
picard --version
spades.py --help
```


---

## How to run SRARec

Running SRARec requires four decisions:

1. which **reference genome** to use;
2. which **input type** to analyse;
3. which **analysis mode** to run;
4. which **detection thresholds** to apply.

### Step 1. Choose the reference genome

SRARec can use built-in references or a custom reference.

| Reference type | How to use it |
|---|---|
| Default SARS-CoV-2 reference | Run without `--virus`, or use the SARS-CoV-2 option supported by the current script |
| HIV-1 reference | Use `--virus hiv` |
| Custom reference | First create the custom index with `--index_only --custom_ref`; then run the analysis using `--virus <custom_name>` |

In SRARec v1.2, `--virus` should be understood mainly as a **reference selector**. For custom references, the FASTA path is provided during indexing with `--custom_ref`, not directly as the main analysis argument.
SRARec automatically prepares the reference indexes required for mapping. For built-in references, no manual indexing step is needed.

#### Built-in SARS-CoV-2 reference

If no virus is specified, SRARec uses the default SARS-CoV-2 reference:

```text
NC_045512.2.fna
```

Example:

```bash
python3 SRARec1.2.py SRRXXXXXXXX --slow
```

#### Built-in HIV-1 reference

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --virus hiv
```

#### Custom reference

First, create or overwrite the custom reference index:

```bash
python3 SRARec1.2.py \
  --index_only \
  --custom_ref path/to/custom_reference.fasta
```

Then run SRARec using the custom reference index:

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --virus <custom_name>
```

Any `--virus` value that is not recognized as a built-in reference is interpreted as a custom reference selector and uses the most recently created custom reference index.

### Step 2. Choose the input type

SRARec can start from three main input types.

| Input type | Use when | Basic command pattern |
|---|---|---|
| **SRA accession** | The reads are public in NCBI-SRA | `python3 SRARec1.2.py SRRXXXXXXXX [options]` |
| **Local FASTQ** | The reads are already available locally | `python3 SRARec1.2.py --simulacion reads.fastq output_folder [options]` |
| **Precomputed BAM** | Mapping has already been performed | `python3 SRARec1.2.py sample_folder --bam [options]` |

#### SRA accession

Accepted accession formats include:

```text
SRRXXXXXXXX
ERRXXXXXXXX
```

In SRA mode, SRARec retrieves sequencing metadata from NCBI, downloads the SRA file, converts it to FASTQ and then runs the selected pipeline.

#### Local FASTQ files

Use `--simulacion` for local FASTQ files.

Single-end input:

```bash
python3 SRARec1.2.py \
  --simulacion \
  reads.fastq \
  sample_output_folder \
  --slow \
  --platform ILLUMINA \
  --virus <reference_selector>
```

Paired-end input:

```bash
python3 SRARec1.2.py \
  --simulacion \
  reads_R1.fastq \
  reads_R2.fastq \
  sample_output_folder \
  --slow \
  --platform ILLUMINA \
  --virus <reference_selector>
```

Supported platform values include:

```text
ILLUMINA
OXFORD_NANOPORE
ONT
PACBIO
PACBIO_SMRT
```

#### Precomputed BAM file

Use `--bam` when the BAM file already exists.

```bash
python3 SRARec1.2.py \
  sample_folder \
  --bam \
  --platform ILLUMINA \
  --virus <reference_selector>
```

In the current v1.2 workflow, SRARec expects the BAM file to be inside `sample_folder` and to be named:

```text
sample_folder.bam
```

This mode skips FASTQ preprocessing and mapping, then proceeds to the read-level polymorphism and four-gamete-test analysis.

### Step 3. Choose the analysis mode

| Mode | How to activate it | Recommended use |
|---|---|---|
| **Fast mode** | Default; do not add `--slow` | Faster exploratory analyses or high-coverage datasets |
| **Slow mode** | Add `--slow` | Standard analyses; recommended default for publication-oriented analyses |
| **Local/simulation mode** | Add `--simulacion` | Local FASTQ or simulated reads |
| **BAM/RAVI mode** | Add `--bam` | Precomputed BAM files |
| **SHIVER assembly mode** | Add `--use_shiver_assembly` | HIV analyses requiring sample-specific references |

#### Slow mode

The slow mode trims reads before mapping and is recommended for standard analyses:

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --virus <reference_selector>
```

Main steps:

```text
SRA download
  → FASTQ conversion
  → trimming
  → mapping
  → duplicate removal
  → pileup generation
  → read-level polymorphism extraction
  → four-gamete-test analysis
```

#### Fast mode

The fast mode is used by default when `--slow` is not provided:

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --virus <reference_selector>
```

Main steps:

```text
SRA download
  → FASTQ conversion
  → initial mapping
  → duplicate removal
  → initial polymorphism screening
  → extraction of informative reads
  → trimming
  → remapping
  → duplicate removal
  → pileup generation
  → read-level polymorphism extraction
  → four-gamete-test analysis
```

This mode reduces computational load by focusing the final analysis on reads that cover more than one polymorphic position.

### Step 4. Choose detection thresholds

SRARec v1.2 allows one or several threshold combinations to be tested.

| Argument | Meaning | Default |
|---|---|---:|
| `--MIN_READS_LIST` | Minimum number of reads supporting a variant at a polymorphic position | `5` |
| `--READS_GAMETOS_MINORITARIOS_LIST` | Minimum number of reads supporting each minor gamete | `1` |
| `--THRESHOLD_GAMETO_FREQ_LIST` | Minimum frequency required for each gamete | `0.00` |

Each argument accepts a comma-separated list.

#### Recommended conservative settings

For conservative recombination detection, use stricter thresholds than the exploratory defaults:

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05 \
  --virus <reference_selector>
```

These settings require:

```text
at least 10 reads supporting a polymorphic variant
at least 5 reads supporting each gamete
at least 5% frequency for each gamete
```

#### Testing several threshold combinations

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --MIN_READS_LIST 5,10,20 \
  --READS_GAMETOS_MINORITARIOS_LIST 1,5,10 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.00,0.05,0.10 \
  --virus <reference_selector>
```

SRARec runs the final recombination analysis for every combination of the provided thresholds.

#### Ignoring indels

To restrict the analysis to base substitutions and ignore insertions and deletions, add:

```bash
--noindel
```

Example:

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --noindel \
  --virus <reference_selector>
```

---

## Usage examples

### Standard SARS-CoV-2 analysis from SRA

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05
```

### Exploratory SARS-CoV-2 analysis from SRA

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --MIN_READS_LIST 5 \
  --READS_GAMETOS_MINORITARIOS_LIST 1 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.00
```

### HIV-1 analysis from SRA

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --virus hiv \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05
```

### Custom-reference analysis from SRA

Index the reference first:

```bash
python3 SRARec1.2.py \
  --index_only \
  --custom_ref references/my_virus.fasta
```

Then run the analysis:

```bash
python3 SRARec1.2.py SRRXXXXXXXX \
  --slow \
  --virus my_virus \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05
```

### Local paired-end Illumina analysis

```bash
python3 SRARec1.2.py \
  --simulacion \
  sample_R1.fastq \
  sample_R2.fastq \
  sample_output \
  --slow \
  --platform ILLUMINA \
  --virus <reference_selector>
```

### Local ONT analysis

```bash
python3 SRARec1.2.py \
  --simulacion \
  reads.fastq \
  sample_output \
  --slow \
  --platform ONT \
  --virus <reference_selector>
```

### BAM-based analysis

```bash
python3 SRARec1.2.py \
  sample_output \
  --bam \
  --platform ILLUMINA \
  --virus <reference_selector>
```

---

## Command-line reference

```text
usage: SRARec1.2.py [options] <SRR>
       SRARec1.2.py --simulacion <fastq1> [fastq2] <output_folder> [options]
       SRARec1.2.py <sample_folder> --bam [options]
```

### Main arguments

| Argument | Description |
|---|---|
| `<SRR>` | NCBI-SRA accession, for example `SRRXXXXXXXX`, `ERRXXXXXXXX` or `DRRXXXXXXXX` |
| `--slow` | Use the slow pipeline. If omitted, the fast pipeline is used by default |
| `--simulacion` | Use local FASTQ files instead of downloading from SRA |
| `--bam` | Start from an existing BAM file |
| `--virus TEXT` | Reference selector. Built-in options include SARS-CoV-2/default and HIV aliases such as `hiv`/`vih`; any other value uses the latest custom reference index |
| `--custom_ref PATH` | Path to a custom reference FASTA. Used with `--index_only` |
| `--index_only` | Create or overwrite the custom reference index and exit |
| `--platform TEXT` | Sequencing platform for local/simulation or BAM mode |
| `--noindel` | Ignore insertions and deletions |
| `--verbose` | Activate verbose output |

### Detection-threshold arguments

| Argument | Description | Default |
|---|---|---:|
| `--MIN_READS_LIST` | Comma-separated list of minimum read counts per polymorphic variant | `5` |
| `--READS_GAMETOS_MINORITARIOS_LIST` | Comma-separated list of minimum read counts per gamete | `1` |
| `--THRESHOLD_GAMETO_FREQ_LIST` | Comma-separated list of minimum gamete-frequency thresholds | `0.00` |

### SHIVER/SPAdes arguments

| Argument | Description | Default |
|---|---|---|
| `--use_shiver_assembly` | Run SPAdes in RNA-viral mode followed by SHIVER before BAM/RAVI analysis | disabled |
| `--spades` | Path or command name for the SPAdes executable. Use `spades.py`; SRARec adds the RNA-viral mode internally | `spades.py` |
| `--spades_threads` | Threads for SPAdes | `8` |
| `--spades_mem_gb` | Memory in GB for SPAdes | `12` |
| `--shiver_threads` | Threads for SHIVER | `6` |
| `--shiver_dir` | Relative path to the SHIVER directory, or the SHIVER directory available in the current environment | `"$CONDA_PREFIX"` |
| `--shiver_ref_alignment` | Relative path to the SHIVER reference alignment | `MyRefAlign.fasta` |
| `--shiver_adapters` | Relative path to the adapter FASTA used by SHIVER | `adapters.fa` |
| `--shiver_primers` | Relative path to the primer FASTA used by SHIVER. If no separate primer FASTA is used, `adapters.fa` can also be passed here | `adapters.fa` |
| `--spades_conda_env` | Optional Conda environment for SPAdes | `none` |

---


## Optional HIV assembly-based workflow

SRARec v1.2 includes an optional HIV-focused workflow combining:

```text
trimming
  → SPAdes in RNA-viral mode
  → SHIVER
  → BAM generation
  → BAM/RAVI recombination analysis
```

This can be useful when a sample-specific reference is preferred before mapping reads for recombination analysis.

Activate it with:

```bash
--use_shiver_assembly
```

The SPAdes executable should be provided as:

```bash
--spades spades.py
```

SRARec handles the RNA-viral SPAdes mode internally. 

If the project includes a specific primer FASTA, replace the last value with that primer file. If no separate primer file is used, `adapters.fa` can be passed as both the adapter file and the primer file.

### Example using local paired-end FASTQ files

```bash
python3 SRARec1.2.py \
  --simulacion \
  --use_shiver_assembly \
  reads_R1.fastq \
  reads_R2.fastq \
  hiv_sample \
  --virus hiv \
  --platform ILLUMINA \
  --spades spades.py \
  --shiver_dir "$CONDA_PREFIX" \
  --shiver_ref_alignment MyRefAlign.fasta \
  --shiver_adapters adapters.fa \
  --shiver_primers adapters.fa
```

### Example using one SRA accession

```bash
python3 SRARec1.2.py \
  --use_shiver_assembly \
  SRRXXXXXXXX \
  --virus hiv \
  --platform ILLUMINA \
  --spades spades.py \
  --shiver_dir "$CONDA_PREFIX" \
  --shiver_ref_alignment MyRefAlign.fasta \
  --shiver_adapters adapters.fa \
  --shiver_primers adapters.fa
```

---
## Outputs

The main SRARec output is written to:

```text
results_<mode>.txt
```

where `<mode>` corresponds to the selected workflow, for example:

```text
results_fast.txt
results_slow.txt
results_ravi.txt
```

Each line corresponds to one analysed sample and contains fields such as:

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

For each recombination-compatible signal, SRARec reports:

- genomic position pair (**1-based** coordinates in `SRARec1.3` and current `generate_graphs.py` output);
- alleles involved;
- observed gamete combinations;
- gamete frequencies;
- supporting read counts;
- number of times the event was recovered as true in broader comparisons.

> **Coordinate system**  
> From **SRARec v1.3** onward, positions written to `Rec`, `No_Rec` and `BegFin_Reads` use **1-based** indexing (first nucleotide = 1).  
> Older result files that contain position `0` are still supported by `generate_graphs.py`, which detects 0-based legacy files automatically.

### Optional VCF output for potential coinfections

When potential coinfection signal is detected, SRARec can generate a VCF file using `bcftools`.

VCF files are written to:

```text
vcf_coinfecciones/
```

with one file per sample:

```text
vcf_coinfecciones/<sample_id>.vcf
```

---

## Plotting recombination profiles

The repository includes the auxiliary script:

```text
generate_graphs.py
```

at the repository root. It generates genome-wide recombination profiles, interactive HTML figures, PNG plots, and per-position export files from SRARec result files.

It can read SRARec-style result files containing fields such as:

```text
Sra:
Rec:
No_Rec:
BegFin_Reads:
```

or tabular files containing columns such as:

```text
pA
pB
true
Rec
No_Rec
```

For each event with positions `pA` and `pB`, the script assigns the weight:

```text
w = 1 / (pB - pA)
```

Each event contributes this weight across the genomic interval:

```text
[pA, pB)
```

True and false contributions are computed separately. The recombination ratio at position `i` is calculated as:

```text
rho(i) = alpha × true(i) / [true(i) + false(i)]
```

when events are available at that position. Plot axes use **1-based** genomic positions.

### CLI example

```bash
python3 generate_graphs.py results_slow.txt \
  --annotation SC2.gb \
  --alpha 1.0
```

With explicit output HTML:

```bash
python3 generate_graphs.py results_ravi.txt \
  -o output/my_sample_profile.html \
  --annotation SC2.gb
```

### Main arguments

| Argument | Description |
|---|---|
| `results_txt` | Input SRARec result file (positional) |
| `-o`, `--output` | Output HTML path |
| `--annotation` | GenBank annotation file (e.g. `SC2.gb`, `HIV-1.gb`) |
| `--genome-len INT` | Genome length, if not inferred automatically |
| `--alpha FLOAT` | Scaling factor for ρ (default: `1.0`) |
| `--exclusive-positions` | Prevent true/false signals from overlapping at the same position |
| `--aggregate-by-patient` | Aggregate events by SRA/sample before summing |
| `--no-coverage` | Skip the coverage curve |
| `--no-begfin-concentration` | Skip read start/end concentration |

### Generated files

When run from the command line, the script can write:

| File | Content |
|---|---|
| `*_recomb_profile.html` | Interactive Plotly profile |
| `*_recomb_profile.png` | Static matplotlib figure |
| `*_recombination_values.txt` | Recombination ratio ρ per position |
| `*_coverage_values.txt` | Read coverage from `BegFin_Reads` |

The script also computes read start/end concentration curves used in the HTML and PNG plots.

---

## SRARec Desktop (macOS GUI)

**SRARec Desktop** is a standalone macOS application for read-level recombination analysis from **preprocessed BAM** files. It wraps the `SRARec1.3.py` analysis engine and `generate_graphs.py` plotting tools in a graphical workflow — no terminal required for routine BAM-based analyses.

### What it does

- Analyse **one or many BAM files in series** from a single run
- Built-in references: **SARS-CoV-2** and **HIV-1**
- Adjustable detection thresholds (min reads, gamete support, frequency)
- Embedded **recombination profile** plot, updated after each sample in a batch
- Export results in multiple formats (gamete summary, ρ per position, coverage, BegFin, raw lines)
- Regenerate or save **PNG / HTML** profiles for selected samples
- Open existing `results_*.txt` files without re-running an analysis

### Download (macOS)

**[Download SRARec Desktop 1.0.0 (zip)](https://github.com/ldgonzalezvazquez/SRARec/releases/latest/download/SRARec-Desktop-1.0.0-macos.zip)** (~300 MB)

The link starts the download immediately in your browser. The archive contains the app bundle and a one-click installer script (`Install SRARec Desktop.command`).

| Requirement | Detail |
|---|---|
| Platform | macOS (Apple Silicon or Intel) |
| Input | Coordinate-sorted, indexed BAM (`.bam` + `.bai`) |
| Disk space | ~800 MB installed |

### Installation

1. **Download** the zip file using the link above (your browser should start the download immediately).
2. **Unzip** the archive.
3. **Double-click** `Install SRARec Desktop.command`.
4. If macOS blocks the installer, open **System Settings → Privacy & Security** and click **Open Anyway**. Confirm if prompted. The installer will then run automatically in Terminal.
5. The script removes the macOS download block, installs **SRARec Desktop** in **Applications**, and launches the app.
6. Once installation is complete, you can delete the unzipped folder and the downloaded `.zip` file.

If you have any trouble installing, contact us at **luisdaniel.gonzalez@uvigo.es** or **luisde365@gmail.com**.

### Quick start

1. Launch **SRARec Desktop**.
2. Click **Add BAMs...** and select one or more `.bam` files (hold ⌘ to multi-select).
3. Choose **sequencing platform**, **reference**, and **output folder** (default: `~/SRARec_runs`).
4. Adjust thresholds if needed, then click **Run analysis**.
5. Combined results are written to `{output_folder}/results_ravi.txt` (one line per sample).
6. Use **Regenerate graph**, **Save graph**, or **Save results...** to export figures and tables — including from an existing results file without running a new analysis.

### BAM requirements

| Requirement | Notes |
|---|---|
| Coordinate-sorted | `@HD SO:coordinate` in BAM header |
| Indexed | `.bai` file present (the app can create it) |
| Deduplicated | Recommended |
| Paired-end | Mates should be fixed (`samtools fixmate`) |

Mapping, trimming, and SRA download must be done **before** using SRARec Desktop (e.g. with `SRARec1.3.py`).

---

## Method summary

SRARec detects recombination-compatible signals by analysing physically linked variants in individual reads.

### Mapping strategy

SRARec selects the mapper according to the sequencing platform:

| Platform | Mapper | Preset |
|---|---|---|
| `ILLUMINA` | BWA-MEM | `bwa mem -M` |
| `OXFORD_NANOPORE`, `ONT` | minimap2 | `-ax map-ont` |
| `PACBIO`, `PACBIO_SMRT` | minimap2 | `-ax map-pb` |

After mapping, SRARec applies filtering and processing steps including:

```text
MAPQ > 30
minimum aligned fraction of read = 0.9
removal of unmapped reads
removal of secondary alignments
removal of supplementary alignments
sorting
mate fixing for paired-end reads
duplicate removal
BAM indexing
```

### Candidate polymorphic-site detection

SRARec uses `samtools mpileup` to identify candidate polymorphic positions. Candidate positions are retained when at least two variants are supported by the selected minimum number of reads.

The default base-quality filter is:

```text
minimum base quality = 20
```

SRARec can analyse:

- substitutions;
- deletions;
- insertions;
- combined base + insertion events.

Insertions are represented at their leftmost coordinate and are retained only when all inserted nucleotides satisfy the base-quality filter.

### Four-gamete test implementation

For each pair of polymorphic positions, SRARec evaluates whether all four allele combinations are observed among reads spanning both positions.

For two positions, `A` and `B`, and two variants per position, a recombination-compatible signal requires evidence for:

```text
A1-B1
A1-B2
A2-B1
A2-B2
```

A pair is classified as recombination-compatible only when all required gametes pass the selected read-count and frequency thresholds.

### Recursive comparison strategy

Within each analysis group, SRARec starts by comparing the two most distant positions.

- If the test is negative, the complete combination is considered negative.
- If the test is positive, the comparison is refined using the median position of the group.
- Adjacent positions from different groups within the same window are also compared when they do not share a group.

---

## Validation

SRARec was evaluated using three complementary strategies.

### 1. Simulated SARS-CoV-2 sequences

The SARS-CoV-2 reference genome `NC_045512.2` was mutated to generate pairs of sequences with different levels of diversity:

```text
π = 0.005
π = 0.010
π = 0.050
```

Insertion and deletion rates were also varied:

```text
Insertion rates: 0.000 and 0.003
Deletion rates:  0.000 and 0.002
```

Two homologous recombinant sequences were generated using random breakpoints. Paired-end FASTQ files were then simulated with **NGSNGS** using different coverage values, read sizes and quality settings.

SRARec was run under different detection-filter combinations, evaluating true events, missed events and false positives by comparing detected signals to the simulated breakpoint position with a margin of:

```text
±100 nucleotides
```

In total, this validation included:

```text
93,312 simulations
```

### 2. Recombinant simulations using real SARS-CoV-2 VOC sequences

SRARec was also tested using recombinant genomes generated from real sequences belonging to past SARS-CoV-2 variants of concern.

This validation used:

```text
50 sequences per VOC
9,000 simulations
```

### 3. Comparison with published SARS-CoV-2 recombination studies

SRARec was run on NCBI-SRA accessions from previous SARS-CoV-2 recombination studies, including datasets from Wertheim et al. and Pipek et al.

For this comparison, relaxed filters were used:

```text
minimum reads per polymorphism = 5
minimum reads per gamete       = 1
minimum gamete frequency       = 0%
```

Detected signals were not expected to be identical to those obtained with the original methodologies, but they were evaluated in relation to previously reported recombination sites.

---

## Large-scale analyses

SRARec has been applied to large public sequencing datasets from NCBI-SRA.

Project-scale analyses included SARS-CoV-2 and HIV-1 SRA entries downloaded as of:

```text
2025-09-02
```

Analyses were run in parallel at the:

```text
Galician Supercomputing Center (CESGA)
```

using:

```text
default parameters
slow mode
```

### Analysis records

| Dataset | Status/date | Approximate size |
|---|---:|---:|
| SARS-CoV-2 final analysis | 2025-12-12 | ~600,000 SRAs |
| HIV-1 final analysis | 2025-12-12 | ~5,000 SRAs |


---

## Limitations

SRARec detects **read-supported recombination-compatible signals**, not direct molecular recombination events.

Detected signals may be influenced by:

- local genetic diversity;
- read length;
- sequencing depth;
- mapping quality;
- reference choice;
- indels;
- amplicon design;
- PCR chimeras;
- duplicate removal;
- filtering thresholds;
- local genome complexity;
- read start/end distribution;
- assembly or sample-specific reference reconstruction choices in SHIVER-based workflows.

Therefore, SRARec results should be interpreted as evidence of read-level incompatibility with a simple non-recombining history, ideally together with phylogenetic, epidemiological or biological evidence.

---

## Citation

If you use SRARec, please cite:

```text
González-Vázquez LD, Iglesias-Rivas P,  Arenas , Martin DP. ???
SRARec: read-level detection of recombination signals from SRA and local sequencing data.
Version 1.3.
[Journal / preprint / DOI to be added]
```

BibTeX placeholder:

```bibtex
@software{SRARec,
  title        = {SRARec: read-level detection of recombination signals from SRA and local sequencing data},
  author       = {González-Vázquez, Luis Daniel and Iglesias-Rivas, Paula, and Martin, Arenas, Miguel and Darren P.},
  year         = {2026},
  version      = {1.3},
  license      = {MIT}
}
```

> Update this citation once the manuscript, preprint, DOI or Zenodo record is available.

---

## License

SRARec is distributed under the **MIT License**.

See the [`LICENSE`](LICENSE) file for details.

---

## References and external tools

SRARec uses or interfaces with the following tools and resources:

- [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- [Cutadapt](https://cutadapt.readthedocs.io/)
- [BBMap / BBTools](https://bbmap.org/)
- [BWA](https://github.com/lh3/bwa)
- [minimap2](https://github.com/lh3/minimap2)
- [SAMtools / HTSlib](https://www.htslib.org/)
- [BCFtools](https://samtools.github.io/bcftools/)
- [Picard](https://broadinstitute.github.io/picard/)
- [SPAdes](https://github.com/ablab/spades)
- [SHIVER](https://github.com/ChrisHIV/shiver)
- [NGSNGS](https://github.com/RAHenriksen/NGSNGS)

Suggested references to complete in the manuscript/repository:

```text
Kans, 2013
Leinonen et al., 2011
Martin, 2011
Bushnell et al., 2017
Li and Durbin, 2009
Li, 2011
Li, 2013
Li, 2018
Danecek et al., 2021
Henriksen et al., 2023
Wertheim et al., 2022
Pipek et al., 2024
```

---

## Contact

For questions, bug reports or feature requests, please open an issue in the repository.

Maintainers/contributors listed in the project materials:

```text
Luis Daniel González-Vázquez
Paula Iglesias-Rivas
Miguel Arenas
Darren P. Martin
```
