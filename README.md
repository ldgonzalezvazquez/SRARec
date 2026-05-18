
````markdown
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
````

---

## Installation

Clone the repository:

```bash
git clone https://github.com/USER/SRARec.git
cd SRARec
```

Create a dedicated environment:

```bash
conda create -n srarec python=3.10
conda activate srarec
```

Install Python requirements:

```bash
pip install -r requirements.txt
```

> If the final repository uses a different installation strategy, update this section with the exact command used by the release version.

---

## External dependencies

SRARec internally integrates several established sequencing-processing tools. Make sure these programs are available in your `$PATH` before running large analyses.

| Tool                | Purpose                                                 | Link                                                       |
| ------------------- | ------------------------------------------------------- | ---------------------------------------------------------- |
| **SRA Toolkit**     | Downloading and processing NCBI-SRA runs                | [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools)      |
| **Cutadapt**        | Adapter trimming                                        | [Cutadapt documentation](https://cutadapt.readthedocs.io/) |
| **BBMap / BBTools** | Adapter database and read-processing utilities          | [BBMap / BBTools](https://sourceforge.net/projects/bbmap/) |
| **BWA-MEM**         | Illumina read mapping                                   | [BWA GitHub](https://github.com/lh3/bwa)                   |
| **minimap2**        | ONT and PacBio read mapping                             | [minimap2 GitHub](https://github.com/lh3/minimap2)         |
| **SAMtools**        | BAM processing, indexing and mpileup extraction         | [SAMtools](https://www.htslib.org/)                        |
| **Picard**          | BAM processing and indexing utilities                   | [Picard](https://broadinstitute.github.io/picard/)         |
| **NGSNGS**          | Optional: simulation of sequencing reads for validation | [NGSNGS GitHub](https://github.com/RAHenriksen/NGSNGS)     |

Example installation with Bioconda:

```bash
conda install -c bioconda -c conda-forge \
  sra-tools \
  cutadapt \
  bwa \
  minimap2 \
  samtools \
  picard
```

> BBMap/BBTools and NGSNGS may require separate installation depending on your system.

---

## Input data

SRARec supports two main input modes.

### 1. Local mode

Use local sequencing files:

```text
Single-end FASTQ:
  sample.fastq.gz

Paired-end FASTQ:
  sample_R1.fastq.gz
  sample_R2.fastq.gz

BAM:
  sample.bam
```

### 2. NCBI-SRA mode

Use an accession code from the **NCBI Sequence Read Archive**:

```text
SRRXXXXXXXX
ERRXXXXXXXX
```

SRARec automatically connects to NCBI-SRA and retrieves the corresponding sequencing data.

---

## Reference genomes

SRARec can use predefined viral references or a user-provided custom reference.

Currently supported default references include:

| Virus                     | Default reference   |
| ------------------------- | ------------------- |
| **SARS-CoV-2**            | `NC_045512.2`       |
| **HIV-1**                 | `NC_001802.1`       |
| **Custom virus / genome** | user-provided FASTA |

For custom references, indexing is handled automatically using the corresponding mapping and BAM-processing tools.

---

## Analysis modes

SRARec implements two analysis modes.

### Slow mode

Recommended for most analyses.

```text
trimming
  → mapping
  → polymorphism calling
  → recombination analysis
```

This mode is more complete and is recommended when analysing heterogeneous datasets, moderate coverage, or samples where trimming may strongly affect mapping quality.

### Fast mode

Recommended only for high-coverage datasets.

```text
mapping
  → polymorphism calling
  → read-subset remapping
  → polymorphism calling
  → recombination analysis
```

This mode reduces computational cost by focusing the analysis after an initial mapping and polymorphism-detection step.

---

## Method summary

SRARec detects recombination-compatible signals by analysing physically linked variants in individual reads.

The main steps are:

1. **Read preprocessing**

   For paired-end runs, orphan read pairs are filtered out when necessary.

2. **Adapter trimming**

   Trimming is performed with **Cutadapt**. By default, SRARec can use an adapter database obtained from **BBMap/BBTools**. If the adapter sequences are known, it is recommended to provide a specific adapter file.

3. **Read mapping**

   SRARec maps reads to the selected reference using:

   * **BWA-MEM** for Illumina reads.
   * **minimap2** for ONT and PacBio reads.

4. **BAM processing and filtering**

   Duplicates are removed and mapping files are sorted and indexed.

   The default mapping filters include:

   * `MAPQ > 30`
   * exclusion of unmapped reads
   * base quality `> 20`
   * removal of reads mapping to less than `90%` of their length
   * removal of reads mapping to more than one genomic location

5. **Polymorphic-site detection**

   SRARec extracts pileup information with **SAMtools mpileup** and processes only positions with:

   * at least `2` alleles
   * at least `N` reads supporting the minor allele

   The default value is:

   ```text
   N = 10
   ```

6. **Read-level variant extraction**

   All reads mapping to each polymorphic position are extracted.

   For each read, SRARec identifies the exact combination of variants observed across all polymorphic positions covered by that read.

   This includes:

   * substitutions
   * deletions
   * insertions, counted at their leftmost coordinate

   Insertions are retained only when all inserted nucleotides have minimum base quality `> 20`.

7. **Window construction**

   Analysis windows are created using the maximum read size mapped from each polymorphic position.

8. **Dynamic grouping of positions**

   Within each window, positions are grouped according to the frequency distribution of their variants.

   Positions are dynamically assigned to the same group when their distribution differs by less than `5%` from the average of at least one variant frequency in the group.

   A single position can belong to more than one group.

9. **Four-gamete testing**

   SRARec applies the four-gamete test to pairs of polymorphic positions using the variant combinations observed in individual reads.

10. **Event filtering**

By default, recombinant-compatible events must satisfy:

```text
minimum reads per gamete      = 5
minimum gamete frequency      = 5%
```

These thresholds can be modified during execution.

---

## Four-gamete test implementation

For a pair of polymorphic positions, `A` and `B`, SRARec evaluates whether all four possible combinations of alleles are observed in reads spanning both positions.

For example:

```text
Read 1: A = X, B = X
Read 2: A = X, B = Y
Read 3: A = Y, B = X
Read 4: A = Y, B = Y
```

If all four gametes are observed with enough read support and frequency, the pair is classified as **recombination-compatible**.

This pattern is incompatible with a simple non-recombining history under the assumptions of the classical four-gamete test.

### Recursive comparison strategy

Within each analysis group, SRARec starts by comparing the two most distant positions.

* If the test is negative, the complete combination is considered negative.
* If the test is positive, the comparison is refined using the median position of the group.
* For an even number of positions, SRARec compares:

  * the left position with the position to the right of the median,
  * the right position with the position to the left of the median,
  * and the two positions surrounding the median.

Adjacent positions from different groups within the same window are also compared when they do not share a group.

---

## Outputs

The main output of SRARec is a table of position pairs satisfying the four-gamete test.

For each recombinant-compatible signal, SRARec reports:

* genomic position pair
* alleles involved
* observed gamete combinations
* gamete frequencies
* number of supporting reads
* number of times the event was recovered as true in broader comparisons

SRARec also logs:

* non-recombinant position pairs
* read start and end coordinates from reads contributing to the analysis
* information useful for assessing potential amplicon or read-coverage biases

Typical output files include:

```text
results/
├── recombinant_pairs.tsv
├── non_recombinant_pairs.tsv
├── begfin_reads.tsv
├── readable_events.tsv
├── raw_events.tsv
└── logs/
```

> File names may differ depending on the current release. Update this section with the exact output names used in the final repository version.

---

## Graphical user interface GUI

SRARec includes a simple **tkinter-based graphical user interface** designed to make the tool accessible to users without command-line experience.

The GUI is optimized for:

* **macOS**
* **Windows**
* moderate analyses
* small viral genomes
* low-diversity datasets

For large batches or high-throughput analyses, running SRARec on a high-performance computing system is recommended.

The GUI allows users to:

* choose the input data source:

  * NCBI-SRA accession ID
  * local FASTQ files
  * local BAM file
* select a predefined reference:

  * SARS-CoV-2: `NC_045512.2`
  * HIV-1: `NC_001802.1`
* load a custom reference sequence
* adjust detection parameters:

  * minimum number of reads per gamete
  * minimum gamete frequency
  * minimum polymorphism support
* ignore indels if desired
* switch between `slow` and `fast` modes
* export detected events in raw or readable format
* export interactive genome-wide visual summaries

After each sample is processed, the GUI generates an interactive recombination-ratio summary across the genome. For predefined viral references, gene annotations are retrieved from GenBank and displayed in the figure.

---

## Plotting recombination profiles

The repository includes the auxiliary script:

```text
generate_graphs.py
```

This script generates genome-wide recombination profiles and interactive figures from SRARec output.

It can read:

* SRARec-style result files containing lines such as:

```text
Sra:
Rec:
No_Rec:
BegFin_Reads:
```

* or tabular files containing columns such as:

```text
pA
pB
true
Rec
No_Rec
```

### Main calculation

For each event with positions `pA` and `pB`:

```text
w = 1 / (pB - pA)
```

Each event contributes this weight across the interval:

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

| Curve                    | Description                                                                                  |
| ------------------------ | -------------------------------------------------------------------------------------------- |
| **Coverage**             | Number of reads covering each position, computed using a difference array and cumulative sum |
| **BegFin concentration** | Density of read starts and ends across the genome                                            |

### Main CLI arguments

```bash
python generate_graphs.py \
  --results PATH \
  --annotation PATH \
  --genome-len INT \
  --alpha FLOAT \
  --exclusive-positions \
  --aggregate-by-patient \
  --output-dir PATH
```

Available arguments:

| Argument                    | Description                                                           |
| --------------------------- | --------------------------------------------------------------------- |
| `--results PATH`            | Input SRARec result file                                              |
| `--annotation PATH`         | GenBank or GFF annotation file                                        |
| `--genome-len INT`          | Genome length, if not inferred automatically                          |
| `--alpha FLOAT`             | Scaling factor for rho                                                |
| `--exclusive-positions`     | Prevent true and false signals from co-occurring at the same position |
| `--aggregate-by-patient`    | Aggregate results by SRA/sample before summing                        |
| `--no-coverage`             | Skip the coverage curve                                               |
| `--no-begfin-concentration` | Skip the read start/end concentration curve                           |
| `--output-dir PATH`         | Output directory for figures and files                                |

Example:

```bash
python generate_graphs.py \
  --results results/srarec_results.txt \
  --annotation annotations/NC_045512.2.gb \
  --genome-len 29903 \
  --alpha 1 \
  --aggregate-by-patient \
  --exclusive-positions \
  --output-dir figures/sars_cov_2/
```

---

## Example command-line templates

> The exact command-line entry point should be updated according to the final repository structure.

### Run from an NCBI-SRA accession

```bash
python <SRARec_entrypoint>.py \
  --sra SRRXXXXXXXX \
  --reference references/NC_045512.2.fasta \
  --mode slow \
  --min-polymorphism-reads 10 \
  --min-gamete-reads 5 \
  --min-gamete-frequency 0.05 \
  --output results/SRRXXXXXXXX/
```

### Run from paired-end FASTQ files

```bash
python <SRARec_entrypoint>.py \
  --fastq-r1 data/sample_R1.fastq.gz \
  --fastq-r2 data/sample_R2.fastq.gz \
  --reference references/NC_045512.2.fasta \
  --mode slow \
  --output results/sample/
```

### Run from a BAM file

```bash
python <SRARec_entrypoint>.py \
  --bam data/sample.bam \
  --reference references/NC_045512.2.fasta \
  --mode slow \
  --output results/sample/
```

### Run in fast mode

```bash
python <SRARec_entrypoint>.py \
  --sra SRRXXXXXXXX \
  --reference references/NC_045512.2.fasta \
  --mode fast \
  --output results/SRRXXXXXXXX_fast/
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
Coverage values: 200, 350, 500
Read sizes:      300, 500, 1000
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

* Wertheim et al. 2022
* Pipek et al. 2024

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

| Dataset                       | Status/date | Approximate size |
| ----------------------------- | ----------: | ---------------: |
| SARS-CoV-2 partial analysis   |  2025-08-18 |    ~350,000 SRAs |
| SARS-CoV-2 final analysis     |  2025-12-12 |    ~600,000 SRAs |
| HIV-1 final analysis          |  2025-12-12 |      ~5,000 SRAs |
| Earlier HIV-1 analysis record |  draft note |      ~1,700 SRAs |

> Update this table before release if the final number of analysed SRAs changes.

### Result folders

```text
results/
├── Partial_res_18Ago/
├── SC2/
└── HIV-1/
```

---

## Recommended repository structure

```text
SRARec/
├── README.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── srarec/
│   ├── __init__.py
│   ├── main.py
│   ├── preprocessing.py
│   ├── mapping.py
│   ├── polymorphisms.py
│   ├── four_gamete_test.py
│   ├── outputs.py
│   └── gui.py
├── scripts/
│   └── generate_graphs.py
├── references/
│   ├── NC_045512.2.fasta
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
└── tests/
```

---

## Best-practice recommendations

For conservative analyses, we recommend using:

```text
mode = slow
minimum reads per polymorphism = 10
minimum reads per gamete = 5
minimum gamete frequency = 5%
```

For validation against known or suspected events, filters can be relaxed, but this may increase sensitivity at the cost of specificity.

Recommended use cases:

| Scenario                       | Recommended mode                        |
| ------------------------------ | --------------------------------------- |
| Standard viral SRA analysis    | `slow`                                  |
| High-coverage datasets         | `fast`                                  |
| Low-diversity samples          | `slow`                                  |
| Large HPC batch analyses       | command line                            |
| Small or moderate analyses     | GUI                                     |
| Method comparison / validation | relaxed filters, interpreted cautiously |

---

## Limitations

SRARec detects **read-supported recombination-compatible signals**, not direct molecular recombination events.

Detected signals may be affected by:

* local diversity
* read length
* sequencing depth
* mapping quality
* indels
* amplicon design
* PCR chimeras
* reference choice
* local genome complexity
* filtering thresholds

Therefore, SRARec results should be interpreted as evidence of read-level incompatibility with a simple non-recombining history, ideally combined with additional biological, epidemiological or phylogenetic evidence.

---

## Citation

If you use SRARec, please cite:

```text
Iglesias-Rivas P, González-Vázquez LD, Arenas M.
SRARec: read-level detection of recombination signals from SRA and local sequencing data.
Version 1.1.1.
[Journal / preprint / DOI to be added]
```

BibTeX placeholder:

```bibtex
@software{SRARec,
  title        = {SRARec: read-level detection of recombination signals from SRA and local sequencing data},
  author       = {Iglesias-Rivas, Paula and González-Vázquez, Luis Daniel and Arenas, Miguel},
  year         = {2026},
  version      = {1.1.1},
  license      = {MIT},
  url          = {https://github.com/USER/SRARec}
}
```

---

## License

SRARec is distributed under the **MIT License**.

See the [`LICENSE`](LICENSE) file for details.

---

## References and external tools

SRARec uses or interfaces with several external tools and resources:

* [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra)
* [SRA Toolkit](https://github.com/ncbi/sra-tools)
* [Cutadapt](https://cutadapt.readthedocs.io/)
* [BBMap / BBTools](https://sourceforge.net/projects/bbmap/)
* [BWA](https://github.com/lh3/bwa)
* [minimap2](https://github.com/lh3/minimap2)
* [SAMtools](https://www.htslib.org/)
* [Picard](https://broadinstitute.github.io/picard/)
* [NGSNGS](https://github.com/RAHenriksen/NGSNGS)

Suggested references to complete in the final manuscript/repository:

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

For questions, bug reports or feature requests, please open an issue in the GitHub repository:

[https://github.com/USER/SRARec/issues](https://github.com/USER/SRARec/issues)

Maintainers:

```text
Paula Iglesias-Rivas
Luis Daniel González-Vázquez
Miguel Arenas
```

```
::contentReference[oaicite:2]{index=2}
```

[1]: https://cutadapt.readthedocs.io/?utm_source=chatgpt.com "Cutadapt — Cutadapt 5.2 documentation"
