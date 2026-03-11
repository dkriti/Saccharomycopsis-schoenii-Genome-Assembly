### Assembly QC (QUAST + BUSCO)

cd /assembly/pacbio/dupsTrue-100x/outputs/
quast final_purged_primary.fasta -o quast_final_report-dupsTrue-100x

input_fasta="final_purged_primary.fasta"
output_fasta="final_purged_primary_cleaned-dupsTrue-100x.fasta"
awk '/^>/ \{gsub("/", "_", $0)\}1' $input_fasta > $output_fasta
busco -i final_purged_primary_cleaned-dupsTrue-100x.fasta -o busco_final_report-dupsTrue-100x -l saccharomycetes_odb10 -m genome}
