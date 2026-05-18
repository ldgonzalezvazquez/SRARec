# SRARec

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#graphical-user-interface-gui)
[![Status](https://img.shields.io/badge/status-active%20development-brightgreen.svg)](#)

**SRARec** is an open-source Python tool for detecting read-supported recombination signals from high-throughput sequencing data using a read-level implementation of the **four-gamete test**.

SRARec can analyse sequencing data from:

- **Illumina**
- **Oxford Nanopore Technologies (ONT)**
- **Pacific Biosciences (PacBio)**

and supports both:

- **single-end (SE)** reads
- **paired-end (PE)** reads

SRARec can be run from local sequencing files or directly from public **NCBI Sequence Read Archive (SRA)** accessions.

> **Project note**  
> For analyses associated with the **CoViE study**, use the [`Exclusive_CoViE`](https://github.com/USER/SRARec/tree/Exclusive_CoViE) branch.

---

## Table of contents

- [Overview](#overview)
- [Main features](#main-features)
- [Workflow](#workflow)
- [Installation](#installation)
- [External dependencies](#external-dependencies)
- [Input data](#input-data)
- [Analysis modes](#analysis-modes)
- [Method summary](#method-summary)
- [Four-gamete test implementation](#four-gamete-test-implementation)
- [Outputs](#outputs)
- [Graphical user interface GUI](#graphical-user-interface-gui)
- [Plotting recombination profiles](#plotting-recombination-profiles)
- [Validation](#validation)
- [Large-scale analyses](#large-scale-analyses)
- [Recommended repository structure](#recommended-repository-structure)
- [Citation](#citation)
- [License](#license)
- [References and external tools](#references-and-external-tools)

---

## Overview

**SRARec** was developed to detect recombination-compatible patterns directly from sequencing reads. Instead of relying only on consensus sequences, SRARec examines how variants are physically linked within individual reads and tests whether pairs of polymorphic positions show all four possible allele combinations.

This makes SRARec suitable for analysing viral populations where recombination may occur within mixed infections, intra-host diversity, or heterogeneous sequencing datasets.

The tool is distributed as open-source software under the **MIT License**.

- **Repository:** [SRARec GitHub repository](https://github.com/USER/SRARec)
- **Documentation:** [SRARec documentation](https://github.com/USER/SRARec/wiki)
- **Current publication version:** `v1.1.1`
- **License:** [MIT](LICENSE)

> Replace `USER/SRARec` with the final repository URL before release.

---

## Main features

SRARec provides:

- Direct analysis from **NCBI-SRA accession IDs**.
- Local analysis from:
  - one FASTQ file for **single-end** data,
  - two FASTQ files for **paired-end** data,
  - or a precomputed **BAM** file.
- Support for **Illumina**, **ONT**, and **PacBio** sequencing data.
- Two analysis modes:
  - **slow mode**, recommended for most analyses;
  - **fast mode**, intended for high-coverage datasets.
- Automated trimming, mapping, filtering, polymorphism extraction and recombination testing.
- Read-level reconstruction of variant combinations, including substitutions and indels.
- Adjustable recombination-detection filters.
- Optional exclusion of indels.
- A simple **tkinter-based graphical user interface** for moderate analyses.
- Interactive genome-wide recombination summaries.
- Exportable raw and readable output tables.
- Auxiliary plotting through `generate_graphs.py`.

---

## Workflow

```text
Input data
    │
    ├── Local FASTQ / BAM files
    └── NCBI-SRA accession ID
            │
            ▼
Read preprocessing
    │
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
    ├── duplicate removal
    ├── mapping-quality filtering
    ├── base-quality filtering
    ├── unmapped-read exclusion
    └── multimapping / low-alignment-length filtering
            │
            ▼
Polymorphic-site detection
    │
    └── positions with ≥2 alleles and enough minor-allele support
            │
            ▼
Read-level variant extraction
    │
    └── linked variants across polymorphic sites
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
    ├── recombination event table
    ├── non-recombinant pair log
    ├── read start/end information
    ├── interactive plots
    └── genome-wide recombination profiles
