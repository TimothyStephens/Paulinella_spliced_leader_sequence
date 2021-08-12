#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

# cut -f1 ../features_upstream_of_SLaddition_sites/SL_from_mapped_reads.gt10SLreads.txt | sort | uniq > names2run.txt
BAM="All_combined.local_mapping.coordsorted.filtered.bam"
F="names2run.txt"

#### Start Script
run_cmd "python mapped_reads_with_SL_sequences.py -b $BAM -i $F -o $F.SL_from_mapped_reads.txt --info $F.SL_from_mapped_reads.SL_read_info.txt.gz -s CCGGCTTTTCTG 2> $F.SL_from_mapped_reads.log"


