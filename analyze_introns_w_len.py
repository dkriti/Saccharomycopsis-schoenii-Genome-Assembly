import re
import statistics
from collections import defaultdict
import os

GTF_FILE = "schoenii_annotation.gtf" 
FASTA_FILE = "Schoenii_assembly.fa"
OUTPUT_DATA_FILE = "all_intron_lengths.txt"
OUTPUT_PLOT_FILE = "Figure_Intron_Distribution.png"

# Updated list including regulatory genes and cytoskeleton factors
CULPRIT_KEYWORDS = [
    "RPL", "RPS",       # Ribosomal Proteins
    "ACT1", "COF1",     # Cytoskeleton
    "TEF", "EFB1",      # Translation Factors
    "VMA",              # Vacuolar ATPase
    "YPT", "GPD",       # GTPases and Metabolism
    "YRA1", "SUS1",     # RNA export / Nuclear pore (Regulatory)
    "DBP2",             # RNA Helicase (NMD regulation)
    "GLC7"              # Phosphatase
]
# =================================================

def parse_fasta(fasta_path):
    """Reads genome fasta into a dict for consensus checking"""
    genome = {}
    current_chrom = None
    seq_parts = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_chrom: genome[current_chrom] = "".join(seq_parts)
                    current_chrom = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if current_chrom: genome[current_chrom] = "".join(seq_parts)
        print(f"Loaded genome: {len(genome)} sequences.")
        return genome
    except FileNotFoundError:
        print("Fasta file not found. Skipping consensus sequence checks")
        return None

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def analyze_gtf():
    print(f"Reading {GTF_FILE}...")

    transcripts = defaultdict(list)
    transcript_to_gene = {}
    gene_id_to_name = {}

    # Regex patterns
    re_tid = re.compile(r'transcript_id "([^"]+)"')
    re_gid = re.compile(r'gene_id "([^"]+)"')
    re_gname = re.compile(r'gene_name "([^"]+)"')

    try:
        with open(GTF_FILE, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue

                feature_type = parts[2]
                attributes = parts[8]

                gid_match = re_gid.search(attributes)
                gname_match = re_gname.search(attributes)

                if gid_match and gname_match:
                    gene_id_to_name[gid_match.group(1)] = gname_match.group(1)

                if feature_type == 'exon':
                    tid_match = re_tid.search(attributes)
                    if tid_match and gid_match:
                        t_id = tid_match.group(1)
                        transcripts[t_id].append({
                            'chrom': parts[0],
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'strand': parts[6]
                        })
                        transcript_to_gene[t_id] = gid_match.group(1)
    except FileNotFoundError:
        print(f"ERROR: GTF file '{GTF_FILE}' not found.")
        return

    # --- ANALYZE ---
    intron_lengths = []
    splice_sites = defaultdict(int)
    culprits_found = []

    genome = parse_fasta(FASTA_FILE)
    print(f"\nAnalyzing {len(transcripts)} transcripts...")

    for t_id, exons in transcripts.items():
        exons.sort(key=lambda x: x['start'])
        if len(exons) < 2: continue 

        g_id = transcript_to_gene.get(t_id)
        g_name = gene_id_to_name.get(g_id, "Unknown")
        is_culprit = any(k in g_name for k in CULPRIT_KEYWORDS)

        local_introns = []

        for i in range(len(exons) - 1):
            intron_start = exons[i]['end'] + 1
            intron_end = exons[i+1]['start'] - 1
            length = intron_end - intron_start + 1

            if length < 1: continue

            intron_lengths.append(length)
            local_introns.append(length)

            # Consensus Check
            if genome and exons[0]['chrom'] in genome:
                chrom_seq = genome[exons[0]['chrom']]
                try:
                    seq = chrom_seq[intron_start-1 : intron_end]
                    if exons[0]['strand'] == '-':
                        seq = get_reverse_complement(seq)
                    motif = f"{seq[:2]}-{seq[-2:]}"
                    splice_sites[motif] += 1
                except (IndexError, KeyError):
                    pass

        if is_culprit and local_introns:
             culprits_found.append(f"{g_name} (ID: {g_id}) - Introns: {local_introns} bp")

    if intron_lengths:
        with open(OUTPUT_DATA_FILE, "w") as f:
            for l in intron_lengths:
                f.write(f"{l}\n")
        print(f"\n Exported {len(intron_lengths)} intron lengths to '{OUTPUT_DATA_FILE}'")
        
        try:
            import matplotlib.pyplot as plt
            
            plt.figure(figsize=(10, 6))
            # Histogram
            plt.hist(intron_lengths, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
            
            # Mean and Median
            mean_val = statistics.mean(intron_lengths)
            median_val = statistics.median(intron_lengths)
            plt.axvline(mean_val, color='blue', linestyle='dashed', linewidth=1.5, label=f'Mean: {mean_val:.1f} bp')
            plt.axvline(median_val, color='red', linestyle='dashed', linewidth=1.5, label=f'Median: {median_val:.0f} bp')
            
            plt.xlabel('Intron Length (bp)')
            plt.ylabel('Frequency')
            plt.title('Figure: Distribution of Intron Lengths')
            plt.legend()
            plt.grid(axis='y', alpha=0.5)
            
            plt.savefig(OUTPUT_PLOT_FILE, dpi=300)
            print(f"Generated histogram: '{OUTPUT_PLOT_FILE}'")
            
        except ImportError:
            print("Matplotlib not installed")

    # Print text report
    print("\n" + "="*40)
    print("Intron analysis report")
    print("="*40)
    print(f"Total introns found: {len(intron_lengths)}")
    if intron_lengths:
        print(f"Average Length: {statistics.mean(intron_lengths):.2f} bp")
        print(f"Median Length: {statistics.median(intron_lengths)} bp")
        print(f"Min/Max: {min(intron_lengths)}/{max(intron_lengths)} bp")

    print("\n Consensus splice sites")
    if splice_sites:
        total = sum(splice_sites.values())
        for motif, count in sorted(splice_sites.items(), key=lambda x: x[1], reverse=True)[:5]:
            print(f"  {motif}: {count} ({count/total*100:.1f}%)")
    else:
        print("  (Skipped or none found)")

    print("\n Conserved intron-containing genes")
    culprits_found.sort()
    for c in culprits_found[:100]:
        print(f"  {c}")
    if len(culprits_found) > 100:
        print(f"  ... and {len(culprits_found)-100} more.")

if __name__ == "__main__":
    analyze_gtf()
