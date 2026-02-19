# Saccharomycopsis-schoenii-Genome-Assembly
Scripts used in the data analysis and preparation of the 'Chromosome-Scale Genome Assembly and Characterization of Saccharomycopsis schoenii, a Necrotrophic Predatory Yeast' manuscript


# Genome Analysis Scripts for *Saccharomycopsis schoenii*

This repository contains the custom Python scripts used for the genome annotation, structural analysis, and genomic feature characterization of *Saccharomycopsis schoenii*. 

These scripts accompany the manuscript "Chromosome-Scale Genome Assembly and Characterization of *Saccharomycopsis schoenii*, a Necrotrophic Predatory Yeast
" published in *G3: Genes|Genomes|Genetics*:
> **[Divya] [Kriti]**, et al. (2026). "[Chromosome-Scale Genome Assembly and Characterization of Saccharomycopsis schoenii, a Necrotrophic Predatory Yeast
]". *G3: Genes|Genomes|Genetics*. 

## Dependencies
The scripts are written in standard Python 3. Some scripts require the following third-party libraries:
* `pandas`
* `matplotlib`
* `seaborn`
* `biopython`

You can install these dependencies using pip:
`pip install pandas matplotlib seaborn biopython`

## Repository Contents

The scripts are grouped by the specific biological analysis they perform:

### 1. Genome Annotation & Intron Analysis
* **`annotation_stats.py`**: Parses the structural annotation (`.gtf`) to calculate global statistics, including total gene counts, exon counts, transcript numbers, and the proportion of single-exon vs. multi-exon genes.
* **`analyze_introns.py`**: Analyzes the `.gtf` and `.fa` files to extract intron statistics, validate `GT-AG` consensus splice sites, and identify introns within conserved "usual culprit" yeast gene families (e.g., *RPL*, *ACT1*, *DBP2*, *YRA1*).
* **`analyze_introns_w_len.py`**: An expanded version of the intron analysis that outputs the raw length data (`all_intron_lengths.txt`) and automatically generates a histogram plot (`Figure_Intron_Distribution.png`) to visualize the length distribution skew.

### 2. Retrotransposon (TE) Analysis
* **`retrotransposon_analysis.py`**: Parses the `ltrs.gff3` output to generate a multi-panel summary figure (`LTR_Analysis.png`) detailing LTR abundance per chromosome, length distributions, genomic coordinates, and sequence identity.
* **`retrotransposon_stats.py`**: Merges structural data (`ltrs.gff3`) with lineage classification data (`candidates.fasta.gydb.cls.tsv`) to perform deep profiling. It identifies giant elements (>10 kb), checks for centromere-targeting chromoviruses, analyzes mixed domains, and exports a clean summary to `FINAL_TE_DATASET.csv`.
* **`TE_analysis.py`**: A specialized script to analyze TE protein domains. It verifies the structural order of Integrase (INT) and Reverse Transcriptase (RT) to confidently separate Ty1/Copia-like (`INT...RT`) from Ty3/Gypsy-like (`RT...INT`) superfamilies.

### 3. Specific Genomic Features (Centromeres & MAT Locus)
* **`find_centromeres.py`**: Uses regular expressions to scan the genome assembly (`.fa`) for specific centromere DNA binding motifs (`[AG]TCAC[AG]TG...TGT[AT][TG]G[TG]T`) and calculates their global genomic coordinates.
* **`mat_gene_analysis.py`**: Analyzes tBLASTn output (`mat_search.txt`) to locate mating-type (MAT) genes. It specifically searches for synteny by calculating the genomic distance between putative MAT/HMG elements and the highly conserved *SLA2* flanking gene.

## Usage
Scripts are designed to be run from the command line. Ensure that your input files (e.g., `schoenii_annotation.gtf`, `Schoenii_assembly.fa`, `ltrs.gff3`) are located in the same working directory as the scripts, or update the file path variables at the top of the respective script.

**Example:**
```bash
python analyze_introns_w_len.py
