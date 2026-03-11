#!/bin/bash

INPUT_XML="/genomes/S.schoenii/pacbio/schoenii-revio.xml"
OUTPUT_DIR="/assembly/pacbio_reassembly/dupsTrue"

export LD_PRELOAD=/lib64/libreadline.so.8

/genomes/S.schoenii/pacbio/smrtlink/install/smrtlink-release_13.1.0.221970/bundles/smrttools/current/smrtcmds/bin/pbcromwell run pb_assembly_hifi 
-e $INPUT_XML 
--task-option ipa2_genome_size=15000000 
--task-option ipa2_downsampled_coverage=100 
--task-option ipa2_run_polishing=True 
--task-option ipa2_run_phasing=True 
--task-option ipa2_run_purge_dups=True 
--task-option ipa2_ctg_prefix="ctg." 
--task-option ipa2_reads_db_prefix="reads" 
--task-option ipa2_cleanup_intermediate_files=True 
--output $OUTPUT_DIR 
--nproc $NPROC}
