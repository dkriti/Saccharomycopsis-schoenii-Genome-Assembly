# Genome Analysis Scripts for *Saccharomycopsis schoenii*

This repository contains the custom bash and Python scripts used for the *de novo* genome assembly, annotation, structural analysis, and genomic feature characterization of *Saccharomycopsis schoenii*. 

These scripts accompany the manuscript published in *G3: Genes|Genomes|Genetics*:
> **Divya Kriti**, et al. (2026). "Chromosome-Scale Genome Assembly and Characterization of *Saccharomycopsis schoenii*, a Necrotrophic Predatory Yeast". *G3: Genes|Genomes|Genetics*. 

## Dependencies
The workflows in this repository utilize a combination of standard bioinformatics tools and Python 3. 

**Bash Script Dependencies:**
* PacBio SMRT Link (specifically `pbcromwell` and the Improved Phased Assembler, IPA)
* `QUAST`
* `BUSCO` (requires the `saccharomycetes_odb10` lineage dataset)

**Python Library Dependencies:**
* `pandas`
* `matplotlib`
* `seaborn`
* `biopython`

Python dependencies can be installed using pip:
`pip install pandas matplotlib seaborn biopython`

## Repository Contents

The scripts are grouped sequentially by the specific biological pipeline they perform:

### 1. Genome Assembly & Quality Control
* **`denovo_assembly_with_hifiasm.sh`**: A bash script that executes the *de novo* assembly of PacBio HiFi reads using the PacBio Improved Phased Assembler (IPA) via `pbcromwell`. It is configured to downsample coverage to 100x and performs integrated polishing, phasing, and sequence duplicate purging.
* **`Assembly_QC.sh`**: A bash script to evaluate the final purged assembly. It generates contiguity metrics using QUAST, sanitizes the FASTA headers, and runs BUSCO against the Saccharomycetes lineage to assess genome completeness.

### 2. Genome Annotation & Intron Analysis
* **`annotation_stats.py`** / **`stats_sch_annot.py`**: Parses the structural annotation (`.gtf`) to calculate global statistics, including total gene counts, exon counts, transcript numbers, and the proportion of single-exon vs. multi-exon genes.
* **`analyze_introns.py`**: Analyzes the `.gtf` and `.fa` files to extract intron statistics, validate `GT-AG` consensus splice sites, and identify introns within conserved "usual culprit" yeast gene families (e.g., *RPL*, *ACT1*, *DBP2*, *YRA1*).
* **`analyze_introns_w_len.py`**: An expanded version of the intron analysis that outputs the raw length data (`all_intron_lengths.txt`) and automatically generates a histogram plot (`Figure_Intron_Distribution.png`) to visualize the length distribution skew.

### 3. Retrotransposon (TE) Analysis
* **`retrotransposon_analysis.py`**: Parses the `ltrs.gff3` output to generate a multi-panel summary figure (`LTR_Analysis.png`) detailing LTR abundance per chromosome, length distributions, genomic coordinates, and sequence identity.
* **`retrotransposon_stats.py`**: Merges structural data (`ltrs.gff3`) with lineage classification data (`candidates.fasta.gydb.cls.tsv`) to perform deep profiling. It identifies giant elements (>10 kb), checks for centromere-targeting chromoviruses, analyzes mixed domains, and exports a clean summary to `FINAL_TE_DATASET.csv`.
* **`TE_analysis.py`**: A specialized script to analyze TE protein domains. It verifies the structural order of Integrase (INT) and Reverse Transcriptase (RT) to confidently separate Ty1/Copia-like (`INT...RT`) from Ty3/Gypsy-like (`RT...INT`) superfamilies.

### 4. Specific Genomic Features (Centromeres & MAT Locus)
* **`find_centromeres.py`**: Uses regular expressions to scan the genome assembly (`.fa`) for specific centromere DNA binding motifs (`[AG]TCAC[AG]TG...TGT[AT][TG]G[TG]T`) and calculates their global genomic coordinates.
* **`mat_gene_analysis.py`**: Analyzes tBLASTn output (`mat_search.txt`) to locate mating-type (MAT) genes. It specifically searches for synteny by calculating the genomic distance between putative MAT/HMG elements and the highly conserved *SLA2* flanking gene.

## Usage
Scripts are designed to be run from the command line. Ensure that your input files (e.g., `schoenii-revio.xml`, `schoenii_annotation.gtf`, `Schoenii_assembly.fa`) are located in the expected directories, or update the file path variables at the top of the respective script.

**Examples:**
```bash
# Run genome QC
bash Assembly_QC.sh

# Run intron analysis
python analyze_introns_w_len.py
