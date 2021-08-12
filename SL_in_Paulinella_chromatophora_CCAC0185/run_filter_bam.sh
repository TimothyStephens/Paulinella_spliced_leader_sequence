#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"

NCPUS=6

B1="SRR3221671.local_mapping.coordsorted.bam"
B2="SRR3221671.local_mapping.namesorted.bam"
B3="SRR3221671.local_mapping.filtered.coordsorted.bam"

#### Start Script
run_cmd "samtools sort -@ $NCPUS -n -o $B2 $B1"
run_cmd "./filter_read_pairs.py -b $B2 | samtools sort -@ $NCPUS -o $B3"
run_cmd "samtools index -@ $NCPUS $B3"

