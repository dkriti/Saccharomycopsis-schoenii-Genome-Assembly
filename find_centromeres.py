import re
from Bio import SeqIO

FASTA_FILE = "Schoenii_assembly.fa"

def find_centromeres_with_global_coords():
    print(f"Reading {FASTA_FILE} to calculate Global Offsets...")
    
    chrom_offsets = {}
    current_global_position = 0
    
    genome_records = []
    try:
        for record in SeqIO.parse(FASTA_FILE, "fasta"):
            chrom_offsets[record.id] = current_global_position
            genome_records.append(record)
            current_global_position += len(record)
    except FileNotFoundError:
        print("Fasta file not found")
        return

    print(f"Total genome length (Concatenated): {current_global_position:,} bp")
    print("-" * 80)
    print(f"{'Chrom':<10} | {'Local Start':<12} | {'Local End':<12} | {'GLOBAL START':<15} | {'GLOBAL END':<15}")
    print("-" * 80)

    # Scan for motifs
    pattern = re.compile(r"([AG]TCAC[AG]TG)([ATCGN]{70,120})(TGT[AT][TG]G[TG]T)")

    for record in genome_records:
        # Skip mitochondria if present
        if "Mito" in record.id or "ChrM" in record.id or "mt" in record.id.lower():
            continue

        offset = chrom_offsets[record.id]

        # Search Forward Strand
        for match in pattern.finditer(str(record.seq)):
            local_start = match.start() + 1
            local_end = match.end()
            
            global_start = offset + local_start
            global_end = offset + local_end
            
            print(f"{record.id:<10} | {local_start:<12} | {local_end:<12} | {global_start:<15,} | {global_end:<15,}")

        # Search Reverse Strand
        rev_seq = record.seq.reverse_complement()
        for match in pattern.finditer(str(rev_seq)):
            local_start = len(record) - match.end() + 1
            local_end = len(record) - match.start()
            
            global_start = offset + local_start
            global_end = offset + local_end
            
            print(f"{record.id:<10} | {local_start:<12} | {local_end:<12} | {global_start:<15,} | {global_end:<15,} (Rev)")

if __name__ == "__main__":
    find_centromeres_with_global_coords()
