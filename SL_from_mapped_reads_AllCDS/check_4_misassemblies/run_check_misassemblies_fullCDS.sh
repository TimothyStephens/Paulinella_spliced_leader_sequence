#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu


#### Start Script
B="../All_combined.local_mapping.coordsorted.bam"
parallel -j 30 "python check_4_misassemblies_between_two_points.py --region_size 20 --min_ok_reads 10 -b $B -r {} -o {}.local.check_misassembly --position_info {}.local.position_info > {}.local.log" :::: files2run.txt



