import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

filename = "ltrs.gff3"
data = []

print(f"Reading {filename}...")

try:
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 11: continue

            try:
                data.append({
                    'Start': int(parts[0]),
                    'End': int(parts[1]),
                    'Length': int(parts[2]),
                    'L_LTR_Len': int(parts[5]),
                    'R_LTR_Len': int(parts[8]),
                    'Similarity': float(parts[9]),
                    'Seq_Nr': int(parts[10])
                })
            except ValueError:
                continue
except FileNotFoundError:
    print(f"Error: {filename} not found.")

df = pd.DataFrame(data)

if not df.empty:
    chrom_map = {0: 'ChrI', 1: 'ChrII', 2: 'ChrIII', 3: 'ChrIV', 4: 'ChrV', 5: 'ChrVI'}
    df['Chromosome'] = df['Seq_Nr'].map(chrom_map)

    fig = plt.figure(figsize=(16, 18))
    plt.rcParams.update({'font.size': 13})
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.8, 1], hspace=0.4)

    ax1 = fig.add_subplot(gs[0, 0])
    order = ['ChrI', 'ChrII', 'ChrIII', 'ChrIV', 'ChrV', 'ChrVI']

    sns.countplot(data=df, x='Chromosome', order=order, color='steelblue', ax=ax1, edgecolor='black')

    ax1.set_title("Abundance of LTR Elements per Chromosome", fontsize=15, fontweight='bold', loc='left')
    ax1.set_xlabel("Chromosome", fontsize=14)
    ax1.set_ylabel("Number of Elements", fontsize=14)

    for p in ax1.patches:
        height = p.get_height()
        ax1.annotate(f'{int(height)}', (p.get_x() + p.get_width() / 2., height),
                     ha='center', va='bottom', fontsize=13, fontweight='bold')

    ax4 = fig.add_subplot(gs[2, 1])
    sns.histplot(data=df, x='Similarity', bins=15, color='rebeccapurple', edgecolor='black', ax=ax4)

    ax4.set_title("LTR Pair Sequence Identity", fontsize=16, fontweight='bold', loc='left')
    ax4.set_xlabel("Sequence Identity (%)", fontsize=14)
    ax4.set_ylabel("Frequency (Number of Elements)", fontsize=14)

    ax2 = fig.add_subplot(gs[1, :])
    y_positions = {chrom: i for i, chrom in enumerate(order)}
    max_pos = df['End'].max() * 1.05

    for _, row in df.iterrows():
        chrom = row['Chromosome']
        if pd.isna(chrom): continue
        y = y_positions[chrom]

        ax2.broken_barh([(row['Start'], row['Length'])], (y - 0.25, 0.5),
                        facecolor='#555555', edgecolor='none')

    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels(order, fontsize=14)
    ax2.set_xlim(0, max_pos)
    ax2.set_xlabel("Genomic Coordinate (bp)", fontsize=14)
    ax2.set_ylabel("Chromosome", fontsize=14)
    ax2.set_title("Genomic Distribution of LTR Elements", fontsize=16, fontweight='bold', loc='left')
    ax2.grid(True, axis='x', linestyle='--', alpha=0.5)

    ax3 = fig.add_subplot(gs[0, 1])
    sns.histplot(data=df, x='Length', bins=20, color='teal', edgecolor='black', ax=ax3)

    ax3.set_title("Full Element Length Distribution", fontsize=16, fontweight='bold', loc='left')
    ax3.set_xlabel("Total Element Length (bp)", fontsize=14)
    ax3.set_ylabel("Frequency (Number of Elements)", fontsize=14)
    ax3.axvline(x=5500, color='red', linestyle='--', linewidth=2, label='Canonical Length (~5.5kb)')
    ax3.legend()

    ax5 = fig.add_subplot(gs[2, 0])
    # Plotting both Left and Right LTR lengths to check for symmetry
    ltr_lengths = pd.DataFrame({
        'Length': pd.concat([df['L_LTR_Len'], df['R_LTR_Len']]),
        'Type': ['Left LTR'] * len(df) + ['Right LTR'] * len(df)
    })

    sns.histplot(data=ltr_lengths, x='Length', hue='Type', multiple='dodge', bins=15, shrink=0.8, ax=ax5)
    ax5.set_title("LTR Region Length Distribution", fontsize=16, fontweight='bold', loc='left')
    ax5.set_xlabel("LTR Length (bp)", fontsize=14)
    ax5.set_ylabel("Frequency (Number of LTRs)", fontsize=14)

    plt.tight_layout()
    plt.savefig('LTR_Analysis.png', dpi=300)

else:
    print("No data found.")
