# SRARec

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#graphical-user-interface-gui)
[![Version](https://img.shields.io/badge/version-1.2-green.svg)](#citation)

**SRARec** is an open-source Python tool for detecting read-supported recombination signals from high-throughput sequencing data using a read-level implementation of the **four-gamete test**.

SRARec can analyse sequencing data from the three major sequencing platforms:

- **Illumina**
- **Oxford Nanopore Technologies (ONT)**
- **Pacific Biosciences (PacBio)**

It supports both:

- **single-end (SE)** reads
- **paired-end (PE)** reads

SRARec can be run from local sequencing files, precomputed BAM files, or directly from public **NCBI Sequence Read Archive (SRA)** accessions.

> **Project note**  
> For analyses associated with the **CoViE study**, use the `Exclusive_CoViE` branch.

---

## Table of contents

- [Overview](#overview)
- [Main features](#main-features)
- [Workflow](#workflow)
- [Installation](#installation)
- [External dependencies](#external-dependencies)
- [Quick start](#quick-start)
- [Input data](#input-data)
- [Reference genomes](#reference-genomes)
- [Basic usage](#basic-usage)
- [Detection parameters](#detection-parameters)
- [Recommended conservative settings](#recommended-conservative-settings)
- [Analysis modes](#analysis-modes)
- [HIV assembly-based mode with SHIVER](#hiv-assembly-based-mode-with-shiver)
- [Mapping strategy](#mapping-strategy)
- [Method summary](#method-summary)
- [Polymorphism and indel handling](#polymorphism-and-indel-handling)
- [Four-gamete test implementation](#four-gamete-test-implementation)
- [Outputs](#outputs)
- [Optional VCF output for potential coinfections](#optional-vcf-output-for-potential-coinfections)
- [Graphical user interface GUI](#graphical-user-interface-gui)
- [Command-line reference](#command-line-reference)
- [Example commands](#example-commands)
- [Plotting recombination profiles](#plotting-recombination-profiles)
- [Validation](#validation)
- [Large-scale analyses](#large-scale-analyses)
- [Best-practice recommendations](#best-practice-recommendations)
- [Limitations](#limitations)
- [Recommended repository structure](#recommended-repository-structure)
- [Citation](#citation)
- [License](#license)
- [References and external tools](#references-and-external-tools)
- [Contact](#contact)

---

## Overview

**SRARec** was developed to detect recombination-compatible patterns directly from sequencing reads. Instead of relying only on consensus sequences, SRARec examines how variants are physically linked within individual reads and tests whether pairs of polymorphic positions show all four possible allele combinations.

This read-level strategy makes SRARec suitable for analysing viral populations where recombination may occur within mixed infections, intra-host diversity, or heterogeneous sequencing datasets.

The tool is distributed as open-source software under the **MIT License**.

---

## Main features

SRARec provides:

- Direct analysis from **NCBI-SRA accession IDs**.
- Local analysis from:
  - one FASTQ file for **single-end** data,
  - two FASTQ files for **paired-end** data,
  - or a precomputed **BAM** file.
- Support for **Illumina**, **ONT**, and **PacBio** sequencing data.
- Multiple practical workflows:
  - **fast mode**, used by default;
  - **slow mode**, recommended for standard analyses;
  - **local/simulation mode**, for local FASTQ files or simulated reads;
  - **BAM/RAVI mode**, for analyses starting from a BAM file;
  - **SHIVER/rnaviralSPAdes mode**, for optional HIV sample-specific assembly workflows.
- Automated SRA download, FASTQ conversion, trimming, mapping, filtering, polymorphism extraction and recombination testing.
- Read-level reconstruction of variant combinations, including substitutions and indels.
- Optional exclusion of indels with `--noindel`.
- Adjustable recombination-detection filters.
- Optional VCF generation for potential coinfection analyses.
- A simple **tkinter-based graphical user interface** for moderate analyses.
- Interactive genome-wide recombination summaries.
- Auxiliary plotting through [`generate_graphs.py`](scripts/generate_graphs.py).

---

## Workflow

```text
Input data
    │
    ├── NCBI-SRA accession ID
    ├── Local FASTQ files
    └── Local BAM file
            │
            ▼
Read preprocessing
    │
    ├── SRA download and FASTQ conversion, if needed
    ├── orphan-pair filtering, if needed
    └── adapter trimming
            │
            ▼
Mapping to reference genome
    │
    ├── BWA-MEM for Illumina
    └── minimap2 for ONT / PacBio
            │
            ▼
BAM processing and filtering
    │
    ├── sorting
    ├── mate fixing for paired-end data
    ├── duplicate removal
    ├── mapping-quality filtering
    ├── aligned-read-fraction filtering
    └── indexing
            │
            ▼
Candidate polymorphic-site detection
    │
    └── positions with ≥2 variants and enough read support
            │
            ▼
Read-level variant extraction
    │
    └── linked substitutions and/or indels across polymorphic sites
            │
            ▼
Four-gamete test
    │
    ├── recombinant-compatible pairs
    └── non-recombinant pairs
            │
            ▼
Outputs
    │
    ├── results_<mode>.txt
    ├── recombinant-compatible events
    ├── non-recombinant pair log
    ├── read start/end information
    ├── optional VCFs for potential coinfections
    └── genome-wide recombination profiles
```

---

## Installation

Clone the repository and move into the project directory:

```bash
git clone <SRARec-repository-url>
cd SRARec
```

Create a dedicated Conda environment:

```bash
conda create -n srarec python=3.10
conda activate srarec
```

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

If an `environment.yml` file is provided, the environment can alternatively be created with:

```bash
conda env create -f environment.yml
conda activate srarec
```

Check that the command-line interface is available:

```bash
python3 SRARec.py --help
```

> Replace `SRARec.py` with the actual script name if the executable file has a different name in your repository.

---

## External dependencies

SRARec v1.2 uses several external tools for SRA download, read preprocessing, mapping, BAM processing, variant extraction, optional assembly-based analysis and VCF generation.

Make sure the following tools are installed and available in your `$PATH`.

| Tool | Purpose | Official link |
|---|---|---|
| **SRA Toolkit** | Downloading and converting SRA runs with `prefetch`, `fasterq-dump` and `fastq-dump` | [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools) |
| **Entrez Direct** | Retrieving SRA run metadata with `esearch` and `efetch` | [NCBI Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) |
| **Cutadapt** | Adapter trimming | [Cutadapt documentation](https://cutadapt.readthedocs.io/) |
| **BBMap / BBTools** | Optional repair of orphan paired-end reads with `repair.sh` | [BBMap / BBTools](https://bbmap.org/) |
| **BWA-MEM** | Mapping Illumina reads | [BWA GitHub repository](https://github.com/lh3/bwa) |
| **minimap2** | Mapping ONT and PacBio reads | [minimap2 GitHub repository](https://github.com/lh3/minimap2) |
| **SAMtools** | BAM/SAM processing, sorting, indexing, duplicate removal and pileup generation | [SAMtools / HTSlib](https://www.htslib.org/) |
| **BCFtools** | Optional VCF generation for potential coinfection analyses | [BCFtools](https://samtools.github.io/bcftools/) |
| **Picard** | Reference dictionary creation for some workflows | [Picard](https://broadinstitute.github.io/picard/) |
| **SPAdes / rnaviralSPAdes** | Optional viral assembly before SHIVER-based analysis | [SPAdes GitHub repository](https://github.com/ablab/spades) |
| **SHIVER** | Optional sample-specific HIV reference reconstruction and read mapping | [SHIVER GitHub repository](https://github.com/ChrisHIV/shiver) |
| **NGSNGS** | Optional read simulation for validation analyses | [NGSNGS GitHub repository](https://github.com/RAHenriksen/NGSNGS) |

A basic installation with Conda/Bioconda can be started with:

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

Some tools, especially **BBMap/BBTools**, **SHIVER**, **NGSNGS** and specific **SPAdes/rnaviralSPAdes** installations, may need to be installed manually depending on the system.

After installation, you can test the main dependencies with:

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
rnaviralspades.py --help
```

---

## Quick start

SRARec can be run from:

- an **NCBI-SRA accession**;
- local **FASTQ** files;
- a precomputed **BAM** file;
- or, in HIV-focused workflows, a **SHIVER/rnaviralSPAdes assembly-based mode**.

By default, if no specific reference is provided, SRARec uses the SARS-CoV-2 reference file located at:

```text
SARS_CoV-2/NC_045512.2.fna
```

A custom reference can be provided with:

```bash
--virus path/to/reference.fasta
```

---

## Input data

SRARec supports three main input types.

### 1. NCBI-SRA accession

Use an accession code from the **NCBI Sequence Read Archive**.

Accepted accession formats include:

```text
SRRXXXXXXXX
ERRXXXXXXXX
DRRXXXXXXXX
```

SRARec automatically connects to NCBI-SRA, retrieves the sequencing metadata, downloads the SRA file and converts it to FASTQ.

### 2. Local FASTQ files

Single-end FASTQ:

```text
sample.fastq
sample.fastq.gz
```

Paired-end FASTQ:

```text
sample_R1.fastq
sample_R2.fastq
sample_R1.fastq.gz
sample_R2.fastq.gz
```

Local FASTQ files are analysed with `--simulacion`.

### 3. BAM file

A precomputed BAM file can be analysed using `--bam`.

In the current v1.2 command-line workflow, SRARec expects the BAM file to be located inside the selected sample folder and to be named:

```text
sample_folder.bam
```

---

## Reference genomes

SRARec can use predefined viral references or a user-provided custom reference.

Currently supported/default references include:

| Virus | Default reference |
|---|---|
| **SARS-CoV-2** | `SARS_CoV-2/NC_045512.2.fna` |
| **SARS-CoV-2 reference accession** | `NC_045512.2` |
| **HIV-1 reference accession** | `NC_001802.1` |
| **Custom virus / genome** | user-provided FASTA through `--virus` |

For custom references, indexing is handled by the corresponding mapping and BAM-processing tools when required.

Example:

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --virus path/to/reference.fasta
```

---

## Basic usage

### Run the default pipeline from an SRA accession

```bash
python3 SRARec.py SRRXXXXXXXX \
  --virus path/to/reference.fasta
```

In SRA mode, SRARec automatically:

1. retrieves the sequencing platform and library layout from NCBI;
2. downloads the SRA file with `prefetch`;
3. converts it to FASTQ with `fasterq-dump` or `fastq-dump`;
4. processes the reads according to the selected mode.

---

### Run the slow pipeline

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --virus path/to/reference.fasta
```

The **slow** pipeline performs trimming before mapping and is recommended for standard analyses.

---

### Run the fast pipeline

```bash
python3 SRARec.py SRRXXXXXXXX \
  --virus path/to/reference.fasta
```

The **fast** pipeline is used by default when `--slow` is not provided.

This mode first maps reads, identifies reads covering polymorphic positions, extracts the informative read subset, trims those reads and remaps them before the final recombination analysis.

---

### Run local FASTQ files

Use `--simulacion` for local FASTQ files.

Paired-end example:

```bash
python3 SRARec.py \
  --simulacion \
  reads_R1.fastq \
  reads_R2.fastq \
  sample_output_folder \
  --slow \
  --virus path/to/reference.fasta \
  --platform ILLUMINA
```

Single-end example:

```bash
python3 SRARec.py \
  --simulacion \
  reads.fastq \
  sample_output_folder \
  --slow \
  --virus path/to/reference.fasta \
  --platform ILLUMINA
```

Supported values for `--platform` include:

```text
ILLUMINA
OXFORD_NANOPORE
ONT
PACBIO
PACBIO_SMRT
```

---

### Run from a BAM file

SRARec can also start from a BAM file using `--bam`.

```bash
python3 SRARec.py \
  sample_folder \
  --bam \
  --virus path/to/reference.fasta \
  --platform ILLUMINA
```

In this mode, SRARec expects the BAM file to be located inside `sample_folder` and to be named:

```text
sample_folder.bam
```

The BAM-based workflow is internally treated as `ravi` mode and proceeds directly to the final read-level polymorphism and four-gamete-test analysis.

---

### Ignore indels

To ignore insertions and deletions during polymorphism scanning, use:

```bash
--noindel
```

Example:

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --noindel \
  --virus path/to/reference.fasta
```

This option restricts the analysis to base substitutions and ignores indel events when extracting polymorphic positions.

---

## Detection parameters

SRARec v1.2 allows testing one or several combinations of detection thresholds.

The three main parameters are:

| Argument | Meaning | Default |
|---|---|---:|
| `--MIN_READS_LIST` | Minimum number of reads supporting a variant at a polymorphic position | `5` |
| `--READS_GAMETOS_MINORITARIOS_LIST` | Minimum number of reads supporting each minor gamete | `1` |
| `--THRESHOLD_GAMETO_FREQ_LIST` | Minimum frequency required for each gamete | `0.00` |

Each argument accepts a comma-separated list.

Example with a single parameter combination:

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05 \
  --virus path/to/reference.fasta
```

Example testing several parameter combinations:

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --MIN_READS_LIST 5,10,20 \
  --READS_GAMETOS_MINORITARIOS_LIST 1,5,10 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.00,0.05,0.10 \
  --virus path/to/reference.fasta
```

SRARec will run the final recombination analysis for every combination of the provided thresholds.

---

## Recommended conservative settings

For conservative recombination detection, we recommend using stricter thresholds than the default exploratory settings:

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05 \
  --virus path/to/reference.fasta
```

These settings require:

```text
at least 10 reads supporting a polymorphic variant
at least 5 reads supporting each gamete
at least 5% frequency for each gamete
```

For exploratory analyses, validation against known events or highly restricted datasets, the filters can be relaxed. However, relaxed filters may increase sensitivity at the cost of specificity.

---

## Analysis modes

SRARec v1.2 includes four practical workflows.

### 1. Fast mode

Default mode.

```bash
python3 SRARec.py SRRXXXXXXXX --virus reference.fasta
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

This mode is designed to reduce computational load by focusing the final analysis on reads that cover more than one polymorphic position.

---

### 2. Slow mode

Recommended mode for standard analyses.

```bash
python3 SRARec.py SRRXXXXXXXX --slow --virus reference.fasta
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

This mode is more complete and is recommended when adapter trimming and mapping quality may strongly affect recombination detection.

---

### 3. Local/simulation mode

Used for local FASTQ files or simulated datasets.

```bash
python3 SRARec.py \
  --simulacion \
  reads_R1.fastq \
  reads_R2.fastq \
  sample_output_folder \
  --slow \
  --virus reference.fasta
```

This mode skips SRA download and starts directly from user-provided FASTQ files.

---

### 4. BAM/RAVI mode

Used when a BAM file is already available.

```bash
python3 SRARec.py \
  sample_folder \
  --bam \
  --virus reference.fasta
```

This mode skips FASTQ preprocessing and mapping, then proceeds to sorting, fixing mate information, indexing, pileup generation and recombination analysis.

---

## HIV assembly-based mode with SHIVER

SRARec v1.2 includes an optional workflow that combines:

```text
trimming
  → rnaviralSPAdes
  → SHIVER
  → BAM generation
  → RAVI/BAM-mode recombination analysis
```

This mode can be activated with:

```bash
--use_shiver_assembly
```

Example using local paired-end FASTQ files:

```bash
python3 SRARec.py \
  --simulacion \
  --use_shiver_assembly \
  reads_R1.fastq \
  reads_R2.fastq \
  sample_output_folder \
  --virus path/to/reference.fasta \
  --platform ILLUMINA \
  --shiver_dir software/shiver \
  --rnaviralspades software/SPAdes-4.2.0-Linux/bin/rnaviralspades.py \
  --shiver_ref_alignment MyRefAligment.fasta \
  --shiver_adapters MyAdapters.fasta \
  --shiver_primers MyPrimers.fasta
```

Example using one SRA accession:

```bash
python3 SRARec.py \
  --use_shiver_assembly \
  SRRXXXXXXXX \
  --virus path/to/reference.fasta \
  --platform ILLUMINA \
  --shiver_dir software/shiver \
  --rnaviralspades software/SPAdes-4.2.0-Linux/bin/rnaviralspades.py \
  --shiver_ref_alignment MyRefAligment.fasta \
  --shiver_adapters MyAdapters.fasta \
  --shiver_primers MyPrimers.fasta
```

Additional SHIVER/SPAdes options:

| Argument | Description | Default |
|---|---|---|
| `--shiver_dir` | Path to the SHIVER directory containing `bin/` | `software/shiver` |
| `--rnaviralspades` | Path to `rnaviralspades.py` | `./software/SPAdes-4.2.0-Linux/bin/rnaviralspades.py` |
| `--spades_threads` | Number of threads for rnaviralSPAdes | `8` |
| `--spades_mem_gb` | Memory in GB for rnaviralSPAdes | `12` |
| `--shiver_threads` | Number of threads for SHIVER mapping | `6` |
| `--shiver_ref_alignment` | Reference alignment used by `shiver_init.sh` | `MyRefAligment.fasta` |
| `--shiver_adapters` | Adapter FASTA used by `shiver_init.sh` | `MyAdapters.fasta` |
| `--shiver_primers` | Primer FASTA used by `shiver_init.sh` | `MyPrimers.fasta` |
| `--spades_conda_env` | Optional Conda environment for SPAdes | `none` |

> This workflow is particularly useful when a sample-specific reference is preferred before mapping reads for recombination analysis.

---

## Mapping strategy

SRARec selects the mapper according to the sequencing platform.

| Platform | Mapper | Preset |
|---|---|---|
| `ILLUMINA` | BWA-MEM | `bwa mem -M` |
| `OXFORD NANOPORE`, `ONT`, `OXFORD_NANOPORE` | minimap2 | `-ax map-ont` |
| `PACBIO`, `PACBIO_SMRT` | minimap2 | `-ax map-pb` |

After mapping, SRARec applies several filtering and processing steps:

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

---

## Method summary

SRARec detects recombination-compatible signals by analysing physically linked variants in individual reads.

The main steps are:

### 1. Read preprocessing

For SRA accessions, SRARec obtains platform and layout metadata from NCBI and downloads the sequencing run with `prefetch`.

For paired-end runs, orphan read pairs can be filtered out when necessary.

### 2. Adapter trimming

Trimming is performed with **Cutadapt**.

By default, SRARec can use an adapter database obtained from **BBMap/BBTools**. If the adapter sequences are known, it is recommended to provide a specific adapter file.

### 3. Read mapping

SRARec maps reads to the selected reference using:

- **BWA-MEM** for Illumina reads.
- **minimap2** for ONT and PacBio reads.

### 4. BAM processing and filtering

Duplicates are removed and mapping files are sorted and indexed.

The default mapping filters include:

```text
MAPQ > 30
base quality > 20
mapped length ≥ 90% of read length
```

SRARec also removes:

- unmapped reads,
- secondary alignments,
- supplementary alignments,
- reads with insufficient mapped length.

### 5. Candidate polymorphic-site detection

SRARec extracts pileup information with **SAMtools mpileup** and processes only positions with:

- at least `2` variants,
- at least the selected number of reads supporting each candidate variant.

### 6. Read-level variant extraction

All reads mapping to each polymorphic position are extracted.

For each read, SRARec identifies the exact combination of variants observed across all polymorphic positions covered by that read.

This includes:

- substitutions,
- deletions,
- insertions counted at their leftmost coordinate,
- combined base + insertion events.

Insertions are retained only when all inserted nucleotides have minimum base quality:

```text
base quality > 20
```

### 7. Window construction

Analysis windows are created using the maximum read size mapped from each polymorphic position.

### 8. Dynamic grouping of positions

Within each window, positions are grouped according to the frequency distribution of their variants.

Positions are dynamically assigned to the same group when their distribution differs by less than a defined threshold from the average of at least one variant frequency in the group.

A single position can belong to more than one group.

### 9. Four-gamete testing

SRARec applies the four-gamete test to pairs of polymorphic positions using the variant combinations observed in individual reads.

### 10. Event filtering

Recombinant-compatible events must satisfy the selected thresholds for:

```text
minimum reads per polymorphic variant
minimum reads per gamete
minimum gamete frequency
```

These thresholds can be modified during execution.

---

## Polymorphism and indel handling

SRARec first uses `samtools mpileup` to identify candidate polymorphic positions.

The default base quality filter is:

```text
minimum base quality = 20
```

Candidate positions are retained when at least two variants are supported by the selected minimum number of reads.

SRARec can analyse:

- substitutions,
- deletions,
- insertions,
- combined base + insertion events.

Insertions are represented at their leftmost coordinate. Insertions are filtered using the minimum quality across all inserted bases.

If `--noindel` is used, insertions and deletions are ignored during polymorphism scanning.

---

## Four-gamete test implementation

For each pair of polymorphic positions, SRARec evaluates whether all four allele combinations are observed among reads spanning both positions.

For two positions, `A` and `B`, and two variants per position, a recombination-compatible signal requires evidence for the following four gametes:

```text
A1-B1
A1-B2
A2-B1
A2-B2
```

Equivalent example:

```text
Read 1: A = X, B = X
Read 2: A = X, B = Y
Read 3: A = Y, B = X
Read 4: A = Y, B = Y
```

A pair is classified as recombination-compatible only when all required gametes pass the selected read-count and frequency thresholds.

This pattern is incompatible with a simple non-recombining history under the assumptions of the classical four-gamete test.

For each positive pair, SRARec stores:

```text
position A
position B
number of broader true ranges containing the event
allele combinations involved
minimum number of reads supporting a gamete
minimum gamete frequency
```

### Recursive comparison strategy

Within each analysis group, SRARec starts by comparing the two most distant positions.

- If the test is negative, the complete combination is considered negative.
- If the test is positive, the comparison is refined using the median position of the group.
- For an even number of positions, SRARec compares:
  - the left position with the position to the right of the median,
  - the right position with the position to the left of the median,
  - and the two positions surrounding the median.

Adjacent positions from different groups within the same window are also compared when they do not share a group.

---

## Outputs

SRARec writes its main output to:

```text
results_<mode>.txt
```

where `<mode>` corresponds to the selected workflow, for example:

```text
results_fast.txt
results_slow.txt
results_ravi.txt
```

Each line contains one analysed sample and follows this general structure:

```text
Sra:<sample_id>
Coinf_detected:<True/False>
Min_read_num:<value>
Min_gameto_num:<value>
Min_gameto_freq:<value>
Platform:<platform>
Rec:<recombinant_pairs>
No_Rec:<non_recombinant_pairs>
BegFin_Reads:<read_start_end_coordinates>
```

The `Rec` field contains position pairs that passed the four-gamete test.

The `No_Rec` field contains tested pairs that did not pass the four-gamete test.

The `BegFin_Reads` field stores read start/end coordinates and can be used to evaluate potential coverage or amplicon-related biases.

For each recombinant-compatible signal, SRARec reports:

- genomic position pair,
- alleles involved,
- observed gamete combinations,
- gamete frequencies,
- number of supporting reads,
- number of times the event was recovered as true in broader comparisons.

---

## Optional VCF output for potential coinfections

When potential coinfection signal is detected, SRARec can generate a VCF file using `bcftools`.

VCF files are written to:

```text
vcf_coinfecciones/
```

with one file per sample:

```text
vcf_coinfecciones/<sample_id>.vcf
```

The VCF generation step uses a command equivalent to:

```bash
bcftools mpileup -d 100000 -f reference.fasta deduplicated.bam | \
  bcftools call -mv -Ov --ploidy 1 -o output.vcf
```

---

## Graphical user interface GUI

SRARec includes a simple **tkinter-based graphical user interface** designed to make the tool accessible to users without command-line experience.

The GUI is optimized for:

- **macOS**,
- **Windows**,
- moderate analyses,
- small viral genomes,
- low-diversity datasets.

For large batches or high-throughput analyses, running SRARec on a high-performance computing system is recommended.

The GUI allows users to:

- choose the input data source:
  - NCBI-SRA accession ID,
  - local FASTQ files,
  - local BAM file;
- select a predefined reference:
  - SARS-CoV-2: `NC_045512.2`,
  - HIV-1: `NC_001802.1`;
- load a custom reference sequence;
- adjust detection parameters:
  - minimum number of reads per gamete,
  - minimum gamete frequency,
  - minimum polymorphism support;
- ignore indels if desired;
- switch between `slow` and `fast` modes;
- export detected events in raw or readable format;
- export interactive genome-wide visual summaries.

After each sample is processed, the GUI generates an interactive recombination-ratio summary across the genome. For predefined viral references, gene annotations are retrieved from GenBank and displayed in the figure.

---

## Command-line reference

```text
usage: SRARec.py [options] <SRR>
       SRARec.py --simulacion <fastq1> [fastq2] <output_folder> [options]
       SRARec.py <sample_folder> --bam [options]
```

### Main arguments

| Argument | Description |
|---|---|
| `<SRR>` | NCBI-SRA accession, for example `SRRXXXXXXXX`, `ERRXXXXXXXX` or `DRRXXXXXXXX` |
| `--slow` | Use the slow pipeline. If omitted, the fast pipeline is used by default |
| `--simulacion` | Use local FASTQ files instead of downloading from SRA |
| `--bam` | Start from an existing BAM file |
| `--virus PATH` | Path to the reference FASTA file |
| `--platform TEXT` | Sequencing platform for local/simulation mode |
| `--noindel` | Ignore insertions and deletions |
| `--verbose` | Activate verbose output |

### Detection-threshold arguments

| Argument | Description | Default |
|---|---|---|
| `--MIN_READS_LIST` | Comma-separated list of minimum read counts per polymorphic variant | `5` |
| `--READS_GAMETOS_MINORITARIOS_LIST` | Comma-separated list of minimum read counts per gamete | `1` |
| `--THRESHOLD_GAMETO_FREQ_LIST` | Comma-separated list of minimum gamete-frequency thresholds | `0.00` |

### SHIVER/SPAdes arguments

| Argument | Description | Default |
|---|---|---|
| `--use_shiver_assembly` | Run rnaviralSPAdes + SHIVER before BAM-mode analysis | disabled |
| `--shiver_dir` | Path to SHIVER directory | `software/shiver` |
| `--rnaviralspades` | Path to `rnaviralspades.py` | `./software/SPAdes-4.2.0-Linux/bin/rnaviralspades.py` |
| `--spades_threads` | Threads for rnaviralSPAdes | `8` |
| `--spades_mem_gb` | Memory in GB for rnaviralSPAdes | `12` |
| `--shiver_threads` | Threads for SHIVER | `6` |
| `--shiver_ref_alignment` | Reference alignment for SHIVER | `MyRefAligment.fasta` |
| `--shiver_adapters` | Adapter FASTA for SHIVER | `MyAdapters.fasta` |
| `--shiver_primers` | Primer FASTA for SHIVER | `MyPrimers.fasta` |
| `--spades_conda_env` | Optional Conda environment for SPAdes | `none` |

---

## Example commands

### Conservative SARS-CoV-2 analysis from SRA

```bash
python3 SRARec.py SRRXXXXXXXX \
  --slow \
  --virus SARS_CoV-2/NC_045512.2.fna \
  --MIN_READS_LIST 10 \
  --READS_GAMETOS_MINORITARIOS_LIST 5 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.05
```

### Exploratory SARS-CoV-2 analysis from SRA

```bash
python3 SRARec.py SRRXXXXXXXX \
  --virus SARS_CoV-2/NC_045512.2.fna \
  --MIN_READS_LIST 5 \
  --READS_GAMETOS_MINORITARIOS_LIST 1 \
  --THRESHOLD_GAMETO_FREQ_LIST 0.00
```

### Local paired-end Illumina analysis

```bash
python3 SRARec.py \
  --simulacion \
  sample_R1.fastq \
  sample_R2.fastq \
  sample_output \
  --slow \
  --platform ILLUMINA \
  --virus reference.fasta
```

### Local ONT analysis

```bash
python3 SRARec.py \
  --simulacion \
  reads.fastq \
  sample_output \
  --slow \
  --platform ONT \
  --virus reference.fasta
```

### BAM-based analysis

```bash
python3 SRARec.py \
  sample_output \
  --bam \
  --platform ILLUMINA \
  --virus reference.fasta
```

### HIV analysis with SHIVER assembly

```bash
python3 SRARec.py \
  --simulacion \
  --use_shiver_assembly \
  reads_R1.fastq \
  reads_R2.fastq \
  hiv_sample \
  --platform ILLUMINA \
  --virus HIV-1/NC_001802.1.fasta \
  --shiver_dir software/shiver \
  --rnaviralspades software/SPAdes-4.2.0-Linux/bin/rnaviralspades.py \
  --shiver_ref_alignment MyRefAligment.fasta \
  --shiver_adapters MyAdapters.fasta \
  --shiver_primers MyPrimers.fasta
```

---

## Plotting recombination profiles

The repository includes the auxiliary script:

```text
scripts/generate_graphs.py
```

This script generates genome-wide recombination profiles and interactive figures from SRARec output.

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

### Main calculation

For each event with positions `pA` and `pB`, the script assigns a weight:

```text
w = 1 / (pB - pA)
```

Each event contributes this weight across the genomic interval:

```text
[pA, pB)
```

True and false contributions are computed separately.

The recombination ratio at position `i` is calculated as:

```text
rho(i) = alpha × true(i) / [true(i) + false(i)]
```

when events are available at that position.

If SRA identifiers are present at the start of each line, the script can automatically aggregate results by patient or sample.

The script can also enforce exclusive positions, so true and false events do not overlap at the same genomic coordinate.

### Auxiliary curves from `BegFin_Reads`

The script can also calculate:

| Curve | Description |
|---|---|
| **Coverage** | Number of reads covering each position, computed using a difference array and cumulative sum |
| **BegFin concentration** | Density of read starts and ends across the genome |

### Main arguments

| Argument | Description |
|---|---|
| `--results PATH` | Input SRARec result file |
| `--annotation PATH` | GenBank or GFF annotation file |
| `--genome-len INT` | Genome length, if not inferred automatically |
| `--alpha FLOAT` | Scaling factor for `rho` |
| `--exclusive-positions` | Prevent true and false signals from overlapping at the same position |
| `--aggregate-by-patient` | Aggregate events by SRA/sample before summing |
| `--no-coverage` | Skip the coverage curve |
| `--no-begfin-concentration` | Skip read start/end concentration |
| `--output-dir PATH` | Output directory for figures and files |

Example:

```bash
python3 scripts/generate_graphs.py \
  --results results_slow.txt \
  --annotation annotations/NC_045512.2.gb \
  --genome-len 29903 \
  --alpha 1 \
  --aggregate-by-patient \
  --exclusive-positions \
  --output-dir figures/sars_cov_2/
```

---

## Validation

SRARec was validated using three complementary strategies.

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

Two homologous recombinant sequences were then generated using random breakpoints.

The resulting four-sequence FASTA files were further mutated by increasing:

```text
π += 0.000 or 0.001
Insertion rate += 0.0000 or 0.0003
Deletion rate  += 0.0000 or 0.0002
```

Paired-end FASTQ files were simulated with **NGSNGS** using different:

```text
Coverage values:   200, 350, 500
Read sizes:        300, 500, 1000
Simulated quality: 100
```

SRARec was run under different detection-filter combinations:

```text
minimum reads per polymorphism: 5, 10, 20
minimum reads per gamete:       5, 10, 20
minority gamete frequency:      5%, 10%, 20%
```

Both `fast` and `slow` modes were evaluated.

In total:

```text
93,312 simulations
```

were performed.

For each simulation, true events, missed events and false positives were evaluated by comparing detected signals to the simulated breakpoint position, allowing a margin of:

```text
±100 nucleotides
```

corresponding to approximately `0.3%` of the simulated genome.

### 2. Recombinant simulations using real SARS-CoV-2 VOC sequences

To better approximate the real diversity and substitution distribution observed in SARS-CoV-2, SRARec was also tested using recombinant genomes generated from real sequences belonging to past variants of concern.

This validation used a whole-genome VOC dataset composed of:

```text
50 sequences per VOC
```

The procedure followed the same overall simulation strategy, but diversity filters were removed.

In total:

```text
9,000 simulations
```

were performed.

### 3. Comparison with published SARS-CoV-2 recombination studies

SRARec was also run on NCBI-SRA accessions from two previous SARS-CoV-2 recombination studies:

- Wertheim et al. 2022
- Pipek et al. 2024

For this comparison, filters were relaxed to maximize detection:

```text
minimum reads per polymorphism = 5
minimum reads per gamete       = 1
minimum gamete frequency       = 0%
```

In the Pipek et al. dataset, 118 SRA-NCBI accessions were initially extracted from the supplementary material. Of these, 23 accessions were discarded because they were not SARS-CoV-2 samples or were private experiments.

SRARec results were not expected to be identical to those obtained with the original methodologies, but detected signals were located near recombination sites identified by those studies.

---

## Large-scale analyses

SRARec has been applied to large public sequencing datasets from NCBI-SRA.

The project-scale analyses included SARS-CoV-2 and HIV-1 SRA entries downloaded as of:

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
| SARS-CoV-2 partial analysis | 2025-08-18 | ~350,000 SRAs |
| SARS-CoV-2 final analysis | 2025-12-12 | ~600,000 SRAs |
| HIV-1 final analysis | 2025-12-12 | ~5,000 SRAs |
| Earlier HIV-1 analysis record | draft note | ~1,700 SRAs |

> Update this table before release if the final number of analysed SRAs changes.

### Result folders

```text
results/
├── Partial_res_18Ago/
├── SC2/
└── HIV-1/
```

---

## Best-practice recommendations

For standard analyses, we recommend:

```text
mode = slow
minimum reads per polymorphic variant = 10
minimum reads per gamete = 5
minimum gamete frequency = 5%
```

For high-coverage datasets, the fast mode can reduce runtime by focusing the final analysis on informative read subsets.

For validation, benchmarking or comparison with previous studies, less restrictive settings can be used, but results should be interpreted cautiously.

Recommended use cases:

| Scenario | Recommended mode |
|---|---|
| Standard viral SRA analysis | `slow` |
| High-coverage datasets | `fast` |
| Low-diversity samples | `slow` |
| Large HPC batch analyses | command line |
| Small or moderate analyses | GUI |
| Method comparison / validation | relaxed filters, interpreted cautiously |
| Precomputed alignment/BAM workflows | `--bam` / RAVI mode |
| HIV sample-specific reference workflow | `--use_shiver_assembly` |

---

## Limitations

SRARec detects **read-supported recombination-compatible signals**, not direct molecular recombination events.

Detected signals may be influenced by:

- local genetic diversity,
- read length,
- sequencing depth,
- mapping quality,
- reference choice,
- indels,
- amplicon design,
- PCR chimeras,
- duplicate removal,
- filtering thresholds,
- local genome complexity,
- read start/end distribution,
- assembly or sample-specific reference reconstruction choices in SHIVER-based workflows.

Therefore, SRARec results should be interpreted as evidence of read-level incompatibility with a simple non-recombining history, ideally together with phylogenetic, epidemiological or biological evidence.

---

## Recommended repository structure

```text
SRARec/
├── README.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── SRARec.py
├── scripts/
│   └── generate_graphs.py
├── SARS_CoV-2/
│   └── NC_045512.2.fna
├── HIV-1/
│   └── NC_001802.1.fasta
├── annotations/
│   ├── NC_045512.2.gb
│   └── NC_001802.1.gb
├── examples/
│   ├── local_fastq/
│   ├── local_bam/
│   └── sra_accession/
├── docs/
│   ├── installation.md
│   ├── usage.md
│   ├── gui.md
│   ├── outputs.md
│   └── validation.md
├── software/
│   ├── shiver/
│   └── SPAdes-4.2.0-Linux/
└── tests/
```

> Update this tree if the final repository structure differs from the release layout.

---

## Citation

If you use SRARec, please cite:

```text González-Vázquez LD,Iglesias-Rivas P, Arenas M, Martin DP.
SRARec: read-level detection of recombination signals from SRA and local sequencing data.
Version 1.2.
[Journal / preprint / DOI to be added]
```

BibTeX placeholder:

```bibtex
@software{SRARec,
  title        = {SRARec: read-level detection of recombination signals from SRA and local sequencing data},
  author       = {Iglesias-Rivas, Paula and González-Vázquez, Luis Daniel and Martin, Darren P. and Arenas, Miguel},
  year         = {2026},
  version      = {1.2},
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
Darren P. Martin
Paula Iglesias-Rivas
Miguel Arenas
```
