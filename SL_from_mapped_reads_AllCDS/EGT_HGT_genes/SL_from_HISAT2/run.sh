

## Combine results
cat *.SL_from_mapped_reads.txt | awk 'NR==1 || $1!~"^#"' > scaffolds_with_EGT_and_HGT_genes.txt.SL_from_mapped_reads.txt

## Get SL addition sites with >10 reads supporting
awk -F'\t' 'NR==1 || $4>10' scaffolds_with_EGT_and_HGT_genes.txt.SL_from_mapped_reads.txt > scaffolds_with_EGT_and_HGT_genes.txt.SL_from_mapped_reads.gt10SLreads.txt

## Convert to bed format
awk 'NR>1{printf "%s\t%s\t%s\t%s:%s:%.10f\n", $1, $2-1, $2-1, $3, $4, $4/$3}' scaffolds_with_EGT_and_HGT_genes.txt.SL_from_mapped_reads.gt10SLreads.txt > scaffolds_with_EGT_and_HGT_genes.txt.SL_from_mapped_reads.gt10SLreads.insertion_points.bed



