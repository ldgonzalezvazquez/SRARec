# SRARec – Exlusive CoViE Methods 

[See main for more detailed information](https://github.com/ldgonzalezvazquez/SRARec)

Analyses began from coordinate-sorted, indexed BAM files (one per sample), processed with standard tooling for alignment data handling (sorting/indexing and pileup generation with SAMtools/BCFtools). All computations operated directly on mapped reads; no local re-assembly or consensus calling was used to infer phase.
Candidate polymorphic sites were identified by scanning alignment pileups for multi-allelic evidence. For each such position, the pipeline assembled—at read level—the sets of reads supporting each observed allele (single-nucleotide changes and short indels), applying base-quality screening and excluding reads that appeared to support more than one allele at the same position to avoid ambiguous evidence. Deletions were left-anchored so that extended gaps were not double-counted across adjacent coordinates.
To determine which site pairs could be phased on the same molecules, positions were organized into windows whose span did not exceed the maximum mapped read length observed in a sample. Within each window, positions were grouped when their allele-frequency profiles were similar (frequencies within a narrow tolerance of a running group average for at least one component). This yielded compact sets of co-observed positions suitable for read-level recombination testing.
Recombination was assessed with a read-level four-gamete test: for each co-covered site pair, the four possible haplotypes (AB, Ab, aB, ab) were counted from the intersection of read sets supporting each allele. A pair was classified as four-gamete-positive when all four haplotypes were observed directly on reads (or read pairs) with sufficient support; pairs lacking any haplotype were treated as non-positive. Testing proceeded (i) within groups—starting from the most distant positions and refining toward the median—and (ii) between adjacent groups within the same window, while tracking containment to retain only minimal (non-contained) positive ranges.

Graphic explanation:

<img width="842" height="595" alt="FGT" src="https://github.com/user-attachments/assets/b090ae5b-c646-404d-b54b-ecdf2399d0df" />

<img width="842" height="595" alt="Mapping" src="https://github.com/user-attachments/assets/be81ffa6-4119-473c-a5e0-11743b44e70c" />

<img width="842" height="595" alt="SameDistrWindows" src="https://github.com/user-attachments/assets/26fb8425-b1cd-48ad-adee-d1cf6b3544f6" />

<img width="842" height="595" alt="DiffDistrWindows" src="https://github.com/user-attachments/assets/ff1acff5-f689-44b7-b9d2-d1a8938dac63" />

----
For each sample, the pipeline produced a compact record encoding all tested site pairs. For positive pairs, the record includes the genomic coordinates, whether a pair is contained within a broader positive range, and summary support metrics (minimum per-haplotype read count and minimum per-haplotype frequency among co-covering reads). Non-positive pairs were also retained. These per-sample lines constitute the sole inputs to downstream aggregation and visualization.
Pair records were then annotated to coding features (genes, non-structural proteins, individual proteins) using the reference GenBank. Each pair contributes a uniform weight across the genomic interval spanned by its two sites; specifically, a pair spanning (b−a) bases contributes 1/(b−a) to every covered position. 

Graphic explanation:

<img width="842" height="595" alt="Graphs" src="https://github.com/user-attachments/assets/26861375-d7d5-43e1-8bcc-24a101dc6bbb" />

----
For any sample × region (gene/NSP/protein), A recombination ratio ρ = T/(T+F) is computed genome-wide at nucleotide resolution to produce a positional profile. The self-contained interactive HTML renders heatmaps of log10(ρ) by genes and by proteins/NSPs means, and a genome-wide positional track of ρ.

(See [HTML](https://github.com/ldgonzalezvazquez/SRARec/Recombination__partial_9oct25.html) for graphs and partial results)

**References** (Exclusive for CoViE):
- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.
- Hudson RR, Kaplan NL. Statistical properties of the number of recombination events in the history of a sample of DNA sequences. Genetics. 1985 Aug;111(1):147–64. doi: 10.1093/genetics/111.1.147. PMID: 4029609; PMCID: PMC1202594.

-----
-----
Implementations for HIV-1 study:
Figure adapted from: https://github.com/ChrisHIV/shiver?tab=readme-ov-file

<img width="842" height="595" alt="HIV-1meth" src="https://github.com/user-attachments/assets/69fe8957-f261-448e-8a36-fa2293351d32" />

