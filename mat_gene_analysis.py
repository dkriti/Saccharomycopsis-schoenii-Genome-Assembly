import pandas as pd

def analyze_blast_results(filename):
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    try:
        df = pd.read_csv(filename, sep='\t', names=columns)
        
        unique_queries = df['qseqid'].unique()
        print(f"Queries found: {unique_queries}\n")
        
        # Separate Hits
        sla2_hits = df[df['qseqid'].str.contains("Sla2", case=False)]
        mat_hits = df[df['qseqid'].str.contains("mat|HMG|alpha", case=False)]
        
        if sla2_hits.empty:
            print("No SLA2 hits found")
            return
        
        if mat_hits.empty:
            print("No MAT/HMG/Alpha hits found")
            return

        # Synteny search
        sla2_scaffolds = sla2_hits['sseqid'].unique()
        
        synteny_found = False
        
        for scaffold in sla2_scaffolds:
            # MAT hits
            mats_on_scaffold = mat_hits[mat_hits['sseqid'] == scaffold]
            
            if not mats_on_scaffold.empty:
                synteny_found = True
                print(f"\n Match found on scaffold: {scaffold}")
                
                # Best SLA2 hit
                best_sla2 = sla2_hits[sla2_hits['sseqid'] == scaffold].sort_values('bitscore', ascending=False).iloc[0]
                sla2_start = best_sla2['sstart']
                sla2_end = best_sla2['send']
                print(f"SLA2 Location: {sla2_start} - {sla2_end} (Bitscore: {best_sla2['bitscore']})")
                
                # List all MAT hits & distance
                for idx, row in mats_on_scaffold.iterrows():
                    mat_start = row['sstart']
                    mat_end = row['send']
                    
                    distance = max(0, min(mat_start, mat_end) - max(sla2_start, sla2_end), min(sla2_start, sla2_end) - max(mat_start, mat_end))
                    
                    print(f"MAT Hit ({row['qseqid']}): {mat_start} - {mat_end}")
                    print(f"-> Distance to SLA2: {distance} bp")
                    print(f"-> E-value: {row['evalue']}")
                    print(f"-> Identity: {row['pident']}%")

        if not synteny_found:
            print("\n No scaffolds found containing both SLA2 and MAT hits.")

    except Exception as e:
        print(f"Error processing file: {e}")

if __name__ == "__main__":
    filename = 'mat_search.txt' #output from tblastn 
    analyze_blast_results(filename)
