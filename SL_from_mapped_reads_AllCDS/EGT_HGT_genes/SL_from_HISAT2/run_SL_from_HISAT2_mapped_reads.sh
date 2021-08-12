#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

FILES2RUN="scaffolds_with_EGT_and_HGT_genes.txt"
BAM="All_combined.HISAT2_RNAseq_mapping.coordSorted.bam"

NCPUS=1
NPARALLEL=30

run_analysis() {
	NCPUS="${1}"; shift
	BAM="${1}"; shift
	FILE="${1}"; shift
	
	exec 1> "${FILE}.log" 2>&1
	echo "## $FILE"
	
	## CODE 2 RUN
	run_cmd "echo \"$FILE\" | python mapped_reads_with_SL_sequences.py -b $BAM -o $FILE.SL_from_mapped_reads.txt --info $FILE.SL_from_mapped_reads.SL_read_info.txt.gz -s CCGGCTTTTCTG"
}


export -f run_analysis
parallel -j $NPARALLEL run_analysis "$NCPUS" "$BAM" :::: <(cut -f 1 "$FILES2RUN")
log "Parallel finished with exit status: $?"

echo -e "## Done running!"

