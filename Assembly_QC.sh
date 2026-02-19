{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red52\green55\blue54;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c26667\c27843\c27451;}
\margl1440\margr1440\vieww28600\viewh17280\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs28 \cf2 \expnd0\expndtw0\kerning0
### Assembly QC (QUAST + BUSCO)\
\
cd /assembly/pacbio/dupsTrue-100x/outputs/\
quast final_purged_primary.fasta -o quast_final_report-dupsTrue-100x\
\
input_fasta="final_purged_primary.fasta"\
output_fasta="final_purged_primary_cleaned-dupsTrue-100x.fasta"\
awk '/^>/ \{gsub("/", "_", $0)\}1' $input_fasta > $output_fasta\
busco -i final_purged_primary_cleaned-dupsTrue-100x.fasta -o busco_final_report-dupsTrue-100x -l saccharomycetes_odb10 -m genome}