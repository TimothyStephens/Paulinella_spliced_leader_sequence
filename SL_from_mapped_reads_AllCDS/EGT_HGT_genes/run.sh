


wget http://cyanophora.rutgers.edu/P_micropora/EGT_trees.tar.gz
wget http://cyanophora.rutgers.edu/P_micropora/HGT_trees.tar.gz


## Get list of EGT and HGT gene IDs
ls -1 EGT_trees/*.tre | sed -e 's@.*/@@' -e 's/.tre//' -e 's/QUERY-Paulinella_micropora_V._/Paulinella_micropora_KR01_nuclear___/' > EGT_genes.txt
ls -1 HGT_trees/*.tre | sed -e 's@.*/@@' -e 's/.tre//' -e 's/QUERY-Paulinella_micropora_V._/Paulinella_micropora_KR01_nuclear___/' > HGT_genes.txt


## 
cat /scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3 \
 | awk -F'\t' '$3=="mRNA" {print $1"\t"$9}' \
 | sed -e 's/\([^\t]*\)\tID=\([^;]*\);.*/\1\t\2/' \
 | ./grepf_column.py -c 2 -f <(cat EGT_genes.txt HGT_genes.txt) \
 | cut -f1 | sort | uniq \
 > SL_from_HISAT2/scaffolds_with_EGT_and_HGT_genes.txt


## Get SL addition sites
./grepf_column.py -i ../names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.txt -f EGT_genes.txt | awk -F'\t' '{print $1"\tPos:"$2" Read_cov:"$3" SL_read_cov:"$4}'
./grepf_column.py -i ../names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.txt -f HGT_genes.txt | awk -F'\t' '{print $1"\tPos:"$2" Read_cov:"$3" SL_read_cov:"$4}'


##
./grepf_column.py -f EGT_genes.txt -i /scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3.exon_counts
./grepf_column.py -f HGT_genes.txt -i /scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3.exon_counts



