import pandas as pd
import numpy as np
import sys

def analyze_retrotransposons():
    ltr_data = []
    try:
        with open("ltrs.gff3", "r") as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.split()
                if len(parts) < 11: continue
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    length = int(parts[2])
                    similarity = float(parts[9])
                    seq_nr = int(parts[10])
                    
                    # Map seq_nr to chromosome
                    chr_map = {0: 'ChrI', 1: 'ChrII', 2: 'ChrIII', 3: 'ChrIV', 4: 'ChrV', 5: 'ChrVI'}
                    chrom = chr_map.get(seq_nr, f"Chr{seq_nr+1}")
                    
                    # Create Chr:Start-End
                    key = f"{chrom}:{start}-{end}"
                    
                    ltr_data.append({
                        "Key": key,
                        "Chromosome": chrom,
                        "Start": start,
                        "End": end,
                        "Length_bp": length,
                        "LTR_Identity": similarity
                    })
                except ValueError: continue
    except FileNotFoundError:
        print("ltrs.gff3 not found")
        sys.exit()

    df_ltr = pd.DataFrame(ltr_data)
    print(f" -> Found {len(df_ltr)} putative LTR elements.")

    try:
        df_cls = pd.read_csv("candidates.fasta.gydb.cls.tsv", sep="\t")
        # Extract the coordinate Key from #TE column (e.g., LTR_14::ChrI:584115-589719 -> ChrI:584115-589719)
        df_cls['Key'] = df_cls['#TE'].apply(lambda x: x.split("::")[1] if "::" in x else x)
        print(f" -> Found classification for {len(df_cls)} elements.")
    except FileNotFoundError:
        print("candidates.fasta.gydb.cls.tsv not found.")
        df_cls = pd.DataFrame(columns=['Key', 'Superfamily', 'Clade', 'Domains'])

    df_merged = pd.merge(df_ltr, df_cls[['Key', 'Superfamily', 'Clade', 'Domains']], on='Key', how='left')
    
    df_merged['Superfamily'] = df_merged['Superfamily'].fillna('Unclassified')
    df_merged['Clade'] = df_merged['Clade'].fillna('None')

    median_id = df_merged['LTR_Identity'].median()
    perfect_ltrs = df_merged[df_merged['LTR_Identity'] >= 100.0]
    
    print("\n LTR IDENTITY")
    print(f"Median LTR identity: {median_id:.1f}%")
    print(f"Recent insertions (100% Identity): {len(perfect_ltrs)} elements ({len(perfect_ltrs)/len(df_merged)*100:.1f}%)")
    
    if 'Unclassified' in df_merged['Superfamily'].values:
        print("\n   Median identity by superfamily:")
        print(df_merged.groupby('Superfamily')['LTR_Identity'].median().to_string())

    print("\n Size analysis")
    giants = df_merged[df_merged['Length_bp'] > 10000].sort_values('Length_bp', ascending=False)
    
    if not giants.empty:
        print(f"Found {len(giants)} 'Giant' elements (>10 kb):")
        for _, row in giants.iterrows():
            print(f"  * {row['Key']} | Length: {row['Length_bp']} bp | Superfamily: {row['Superfamily']} | Clade: {row['Clade']}")
            if row['Superfamily'] == 'Unclassified':
                print("    -> Note: Lack of domains + large size strongly supports 'Nested Insertion' hypothesis (coding region disrupted).")
            else:
                print("    -> Note: Has domains but size is double canonical -> Nested Insertion likely retaining some host domains.")
    else:
        print("No elements > 10kb found.")

    # Centromere analysis
    print("\n Centromere/Chromovirus check")
    cen_region = df_merged[
        (df_merged['Chromosome'] == 'ChrVI') & 
        (df_merged['Start'] > 650000) & 
        (df_merged['Start'] < 670000)
    ]
    
    if not cen_region.empty:
        print("Element found near ChrVI Centromere (650kb-670kb):")
        for _, row in cen_region.iterrows():
            print(f"  * {row['Key']} | Lineage: {row['Clade']} | Identity: {row['LTR_Identity']}%")
            if 'cer' in str(row['Clade']).lower() or 'crm' in str(row['Clade']).lower():
                print("This is a Chromovirus (Cer/Crm), known to target centromeres.")
    else:
        print("No elements found in the specific ChrVI centromere window.")

    # Domain analysis
    print("\n Domain analysis")
    mixed = df_merged[df_merged['Clade'] == 'mixture']
    print(f"Total Mixed Lineage Elements: {len(mixed)}")
    if not mixed.empty:
        print("Mixed Domains (Biological Insight: Evolutionary Divergence):")
        print(mixed[['Key', 'Superfamily', 'Domains']].head(3).to_string(index=False))

    df_merged.to_csv("FINAL_TE_DATASET.csv", index=False)
    print(" -> Saved merged dataset to 'FINAL_TE_DATASET.csv'. Use this for your tables.")

if __name__ == "__main__":
    analyze_retrotransposons()
