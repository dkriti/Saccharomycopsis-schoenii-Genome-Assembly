import pandas as pd
import sys

cls_file = "candidates.fasta.gydb.cls.tsv"

try:
    df = pd.read_csv(cls_file, sep='\t')
except FileNotFoundError:
    print(f"Could not find {cls_file}")
    sys.exit()

total_rows = len(df)
unique_elements = df['#TE'].nunique()

print(f"\nTotal entries in file: {total_rows}")
print(f"Total unique TEs classified: {unique_elements}")

print("\n By Superfamily")
superfamily_counts = df['Superfamily'].value_counts()
print(superfamily_counts)

print("\n By Clade (lineage)")
clade_counts = df.groupby(['Superfamily', 'Clade']).size()
print(clade_counts)

# Domain check (Ty1/Copia = INT before RT; Ty3/Gypsy = RT before INT)

def check_order(domains_str):
    if not isinstance(domains_str, str): return "No domains"
    
    parts = domains_str.split('|')
    
    # Positions of substrings
    pos_rt = domains_str.find("RT")
    pos_int = domains_str.find("INT")
    pos_rh = domains_str.find("RNaseH")
    
    if pos_rt == -1 or pos_int == -1:
        return "Incomplete (missing RT or INT)"
    
    if pos_int < pos_rt:
        return "Ty1/Copia-like (INT...RT)"
    elif pos_rt < pos_int:
        return "Ty3/Gypsy-like (RT...INT)"
    else:
        return "Ambiguous"

df['Structure_Check'] = df['Domains'].apply(check_order)

print("\n Domain summary:")
print(df.groupby(['Superfamily', 'Structure_Check']).size())

output_summary = "TE_Analysis_Summary.txt"
with open(output_summary, "w") as f:
    f.write(f"Total unique TEs: {unique_elements}\n\n")
    f.write("Superfamily counts:\n")
    f.write(superfamily_counts.to_string())
    f.write("\n\nClade counts:\n")
    f.write(clade_counts.to_string())

print(f"\nSaved summary to {output_summary}")
