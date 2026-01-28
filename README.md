# SRARec – Methods

**For CoViE study analyses move to the [Exclusive_CoViE](https://github.com/ldgonzalezvazquez/SRARec/tree/Exclusive_CoViE) branch.** 
---
SRARec is developed in Python 3 and distributed as open-source software (MIT License). The source code, along with detailed documentation and the graphical interface, is available in the authors' GitHub repository (link; version 1.1.1 in the publication) for free use and adaptation. Internally, SRARec integrates several proven sequence alignment and processing tools.
The software is designed to process samples from the three major sequencer manufacturers: Illumina, Oxford Nanopore Technologies (ONT), and Pacific Biosciences (PacBio), with both single-end (SE) and paired-end (PE) samples. It automatically connects to the NCBI Sequence Read Archive (SRA–NCBI) database (Kans, 2013; Leinonen et al., 2011). When starting the tool, you can select local mode by specifying one (SE) or two (PE) FASTQ files, a BAM file, or NCBI–SRA modo by entering a database accession code (standard format: (S/E)RRXXXXXXXX).
SRARec implements two modes: (i) a recommended ‘slow’ mode (trimming, mapping, polymorphism calling and analysis), and (ii) a ‘fast’ mode for high-coverage datasets (mapping, polymorphism calling, read subset remapping, polymorphism calling and analysis). Before trimming, and only if necessary, orphan pairs are filtered out Paired-End runs (Bushnell et al., 2017). Trimming is performed using the cutadapt software (Martin, 2011) and by default a database of adapters obtained from BBMap (Bushnell et al., 2017) is used, although if the adapters are known, when using the tool it is recommended to indicate a specific file containing them. After this, the mapping process begins, using BWA-MEM (Li et Durbin, 2009; Li, 2013) for Illumina and minimap2 (Li et Durbin, 2010; Li, 2018) for ONT and PacBio. Duplicates are then removed and several filters are applied to the mapping (MAPQ > 30; exclusion of unmapped reads; base quality > 20; elimination of reads that map < 90% of their length or map in more than one place) in addition to ordering and indexing (Li, 2011; Danecek et al, 2021). With the resulting file, an annotation of the mapping is extracted with SAMtools mpileup (Li, 2011; Danecek et al, 2021) and then the alignment is processed only in the positions with ≥2 alleles and ≥N supporting reads (default N=10) for the minor allele.
All reads mapping each of the polymorphic positions are extracted. They are then individually analyzed to identify the exact variants in each read from all polymorphic positions that map to that read—including indels, counted at their leftmost coordinate. Insertions are also filtered so that the minimum base quality of all the nucleotides comprising them is >20. Taking the entire combination of variants, analysis windows are created with the maximum read size mapped from each polymorphic position. Within each window, analysis groups are created based on the frequency distribution of the variants mapped to each position, dynamically grouping into the same group those positions that have a percentage of the distribution that deviates by less than 5% from the average of at least one of the percentages of the distribution in each group. The same position can belong to more than one group. For each group in each window, the comparison starts with the pair of positions furthest apart. If the test is negative, it is concluded that all combinations are negative. If the test is positive, it is refined by analyzing these positions with the median position of the group, and so on recursively. In the case of an even number of positions, the left position is compared with the position to the right of the median and vice versa, in addition to comparing both positions on either side of the median with each other. Between groups in the same window, all adjacent positions that do not share a group are compared.
In position comparisons, which are made by peers, it is verified if the rule of the 4 gametes applied to the linking we see between variants in individual reads is followed. Assuming that we did not apply filters, the test would be positive if we find at least 4 reads (R1, R2, R3, R4) that map the same 2 positions (A and B) and that 4 combinations/gametes of possible variants are found in possible incompatibles with a story without recombination (R1-A.X-B.X, R2-A.X-B.Y, R3-A.Y-B. R4-a.y-b.y). By default, filters are applied to guarantee events that have at least 5 reads per gamete and 5% frequency among all possible gametes - these filters can be changed during execution, but we recommend at least these values to avoid false positives due to PCR chimeras according to (Pipek et al., 2024).
The main output of SRARec is a table of positions pairs that satisfy the test and therefore indicate recombination. Non-recombinant pairs are also logged, and indicating for each recombinant signal the alleles involved, the observed frequencies of each gamete, the number of supporting reads, and the number of times that an event was given as true in broader comparisons. The program also generates a list with the start and end of the reads that generated results, in case the user wants to analyze the ampliconic bias of the analysis.

User Graphical Interface (GUI)
We recommend running SRARec on high-performance machines for large batches. But, in order to make SRARec accessible to investigators without experience in the command line and for moderate analyses (small genomes or low diversity), we develop a simple tkinter-based graphical user interface and optimized for MacOS and Windows. The GUI displays a window where the user can: select the data source (by entering an SRA–NCBI ID, or by loading local FASTQ or BAM  files), choose the viral reference genome from preconfigured options for SARS-CoV-2 (NC_045512.2 by default), HIV-1 (NC_001802.1 by default), or load a custom sequence – indexing is handle automatically with BWA/SAMtools/Picard. Filtering and detection parameters (minimum gamete frequency or number of support reads per gamete) can also be adjusted, as well as indels can be ignored in the analysis, or the analysis can be switched between “fast” and “slow” modes. Upon completion of the analysis of each sample, the interface generates an interactive visual summary consisting of the recombination ratio (see above) across the genome, with an annotation of the target virus genes (if this was the default) taken from their annotation in GenBank. This can be exported in various image formats. Tables of detected events are also generated, which can be exported in raw form (all information - see previous section) or in a more readable and organized format.

Simulations validation
We performed simulations for method validation using simulated sequences, real sequences with simulated recombinations, and inputs from other studies involving similar analyses.
We mutated the SARS-CoV-2 reference (NC_045512.2) generating 2 sequences with different diversities (= 0.005, 0.01 and 0.05) and insertion (0.0 and 0.003) and deletion (0.0 and 0.002) rates. From these, we generated 2 homologous recombinant sequences with random breakpoints. We mutated the 4 resulting sequences again increasing the diversity (by pi += 0.0 and 0.001) and the insertion (0.0 and 0.0003) and deletion (0.0 and 0.0002) rates. Using the resulting 4-sequence FASTA, we simulated FASTQ (PE) files with the NGSNGS tool (Henriksen et al., 2023) at different coverages (200, 350, and 500) and read sizes (300, 500, and 1000) with a simulated sequencing quality of 100 (a parameter defined in NGSNGS). We ran SRARec with different detection filters (see above): minimum number of reads per polymorphism (5, 10, and 20), minimum number of reads per gamete (5, 10, and 20), and minimum minority gamete frequency (5%, 10%, and 20%). Thus, we performed two simulations for each possible filter combination across all these steps for each of the modes (“fast” and “slow”), resulting in 93,312 simulations. For each of these simulations, we identified the number of true events, missed events, and false positives by comparing the distance to the simulated breakpoint. We allowed a margin of 100 nucleotides on either side, which represents 0.3% of the simulated genome size. These simulations demonstrated the high performance of the tool, especially in the absence of indels, which increase the number of false positives by worsening mapping (see Results).
To recreate the real diversity and distribution of substitutions in SARS-CoV-2, we recombined with random breakpoints real sequences between past VOCs previously classified by us in (REF JMedVir accepted but not published – Misture whole genome dataset composed of 50 sequences of each VOC). The procedure was the same as in the previous case, but we removed the diversity filters. In total, we performed 9,000 simulations that, although good (see results), reflected a lower detection rate and only 2 false positives.
Alternatively, we ran SRARec on the NCBI SRA entries from two articles that analyzed recombination in SARS-CoV-2 using dissimilar but comparable methodologies. For these simulations, we relaxed the tool's filters to 5 reads per polymorphism, and 1 read and 0% per gamete to maximize detection.
In Wertheim et al., 2022, recombination was analyzed using three tests based on the four-gamete hypothesis, performing multiple sequence alignments of all reads flanking the S and N genes in a coinfected patient who was resequenced four times. We ran SRARec on all accessions deposited at SRA-NCBI for this study.
In Pipek et al., 2024, frequency imbalances and defining mutations of concordant variants with recombination were analyzed. Of the 118 SRA-NCBI accessions extracted from the supplementary material of this article and claimed to be SARS-CoV-2, 23 (i.e.: SRR12061247, SRR12859897 or SRR12825768) were not samples of that virus (verified in the NCBI accession itself) or were from private experiments, and therefore had to be discarded. SRARec was run on the remaining entries.
The results with SRARec, although not exactly identical, were located near the recombination sites identified by the other two methodologies. We also added complementary detection –given that both studies encounter regional limitations, as they study either specific regions or only polymorphisms that define variants (see results).

Real-world data analysis (SARS-CoV-2 and HIV-1)
We downloaded the identifiers for all NCBI SARS-CoV-2 and HIV-1 SRA entries as of September 2, 2025 [update at the end with our HIV-1 results]. We ran the analysis in parallel at the Galician Supercomputing Center (CESGA) with default parameters and in "slow" mode (see above). We obtained results on 500,000 SARS-CoV-2 SRAs and 1,700 HIV-1 SRAs.

Partial results [until August 8, 2025 - ≈350,000 SRAs]:
<img width="990" height="550" alt="Partial_res_18Ago" src="https://github.com/user-attachments/assets/431c1083-6629-4817-876f-233fb5c02741" />

Final results [December 12, 2025 - ≈600,000 SC2 SRAs & ≈ 5,000 HIV-1 SRAs from Ravi Paper]:

[SC2](https://github.com/ldgonzalezvazquez/SRARec/blob/main/SC2_20251214_224201.html)


[HIV-1](https://github.com/ldgonzalezvazquez/SRARec/blob/main/HIV-1_20251212_163459.html)

# [generate_graphs.py](https://github.com/ldgonzalezvazquez/SRARec/blob/main/generate_graphs.py)

This script generates recombination profiles along the genome and interactive
figures (HTML) plus optional PNGs. It reads SRARec-style results (lines with
Sra:, Rec:, No_Rec:, BegFin_Reads:) or tables with pA/pB/true or Rec/No_Rec
columns.

Main calculation
- Events (pA, pB) are extracted from Rec (true=1) and No_Rec (true=0).
- Invalid events are filtered (pB <= pA) and pA/pB are coerced to numeric.
- Each event contributes a weight w = 1/(pB - pA) over the range [pA, pB).
- True/false contributions are separated and rho per position is computed as:
  rho(i) = alpha * true / (true + false) when events exist.
- If SRA is present at the start of each line, the script aggregates by patient
  automatically; it can also enforce exclusive positions so true and false do
  not co-occur at the same position.

Auxiliary curves from BegFin_Reads
- Coverage: counts how many reads cover each position using a difference array
  and cumulative sum (range coverage).
- Begfin concentration: counts read starts and ends, adding only at start and
  end of each read (start/end density).
### Main CLI arguments

- `--results PATH` — input results file  
- `--annotation PATH` — GenBank or GFF for gene track (optional)  
- `--genome-len INT` — genome length (if not inferred)  
- `--alpha FLOAT` — rho scale factor  
- `--exclusive-positions` — true and false cannot overlap at a position  
- `--aggregate-by-patient` — aggregate per SRA before summing  
- `--no-coverage` — skip coverage curve  
- `--no-begfin-concentration` — skip begfin concentration curve  
- `--output-dir PATH` — output folder for figures/files  
