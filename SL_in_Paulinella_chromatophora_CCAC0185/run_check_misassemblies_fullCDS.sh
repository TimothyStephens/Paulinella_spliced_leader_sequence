#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

#### Start Script
R="regions2check.txt"
B="SRR3221671.local_mapping.coordsorted.bam"
run_cmd "python check_4_misassemblies_between_two_points.py --region_size 20 --min_ok_reads 10 -b $B -r $R -o $R.local.check_misassembly --position_info $R.local.position_info > $R.local.log"



