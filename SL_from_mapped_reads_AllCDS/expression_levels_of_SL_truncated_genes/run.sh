



## 
## Get the average expression level across all samples from CL+HL conditons and link with the number of SL addition sites found in that sequence
## 

M="Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.RSEM.allLibsCombined.isoforms.results.normalized_counts.txt"

##
## CL only
##
OUT="KR01__SeqID_AvgNormExpr.txt.CL"
SUM_OUT="$M.rowsums.txt.CL"
Rscript sum_matrix_rows.R <(cut -f1-22 $M) $SUM_OUT

SL="../names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.insertion_points.noMissassemblies.bed"
awk 'NR>1 {print $1"\t"$2/21}' $SUM_OUT \
  | ~/scripts/add_value_to_table.py -d 0 -a \
    <(cut -f1 $SL | sort | uniq -c | awk '{print $2"\t"$1}') \
  > $OUT


# Get just the Avg. expression levels of transcripts with Zero or Multipl (>=1) SL additon sites
awk '$3==0 {print $2}' $OUT > $OUT.0_SLsites
awk '$3>0 {print $2}'  $OUT > $OUT.multiple_SLsites




##
## Calculate using TPM
##
M="Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.RSEM.allLibsCombined.isoforms.results.tpm"

##
## CL only
##
OUT="KR01__SeqID_AvgTPM.txt.CL"
SUM_OUT="$M.rowsums.txt.CL"
Rscript sum_matrix_rows.R <(cut -f1-22 $M) $SUM_OUT

SL="../names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.insertion_points.noMissassemblies.bed"
awk 'NR>1 {print $1"\t"$2/21}' $SUM_OUT \
  | ~/scripts/add_value_to_table.py -d 0 -a \
    <(cut -f1 $SL | sort | uniq -c | awk '{print $2"\t"$1}') \
  > $OUT


# Get just the Avg. expression levels of transcripts with Zero or Multipl (>=1) SL additon sites
awk '$3==0 {print $2}' $OUT > $OUT.0_SLsites
awk '$3>0 {print $2}'  $OUT > $OUT.multiple_SLsites




