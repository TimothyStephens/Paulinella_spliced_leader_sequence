#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

BAM="SRR3221671.local_mapping.filtered.coordsorted.bam"
F="names2run.txt"

#### Start Script
run_cmd "python mapped_reads_with_SL_sequences.py -b $BAM -i $F -o $F.filtered.SL_from_mapped_reads.txt --info $F.filtered.SL_from_mapped_reads.SL_read_info.txt.gz -s CCGGCTTTTCTG 2> $F.filtered.SL_from_mapped_reads.log"


