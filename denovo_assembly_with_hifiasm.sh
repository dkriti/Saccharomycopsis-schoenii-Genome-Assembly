{\rtf1\ansi\ansicpg1252\cocoartf2867
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red52\green55\blue54;}
{\*\expandedcolortbl;;\cssrgb\c26667\c27843\c27451;}
\margl1440\margr1440\vieww28600\viewh17280\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs28 \cf0 \expnd0\expndtw0\kerning0
#!/bin/bash\
\
INPUT_XML="/genomes/S.schoenii/pacbio/schoenii-revio.xml"\
OUTPUT_DIR="/assembly/pacbio_reassembly/dupsTrue"\
NPROC=128\
\
export LD_PRELOAD=/lib64/libreadline.so.8\
\
/genomes/S.schoenii/pacbio/smrtlink/install/smrtlink-release_13.1.0.221970/bundles/smrttools/current/smrtcmds/bin/pbcromwell run pb_assembly_hifi \\\
-e $INPUT_XML \\\
--task-option ipa2_genome_size=15000000 \\\
--task-option ipa2_downsampled_coverage=100 \\\
--task-option ipa2_run_polishing=True \\\
--task-option ipa2_run_phasing=True \\\
--task-option ipa2_run_purge_dups=True \\\
--task-option ipa2_ctg_prefix="ctg." \\\
--task-option ipa2_reads_db_prefix="reads" \\\
--task-option ipa2_cleanup_intermediate_files=True \\\
--output $OUTPUT_DIR \\\
--nproc $NPROC}