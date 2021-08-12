

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"



## Analysis
cat names2run.txt.filtered.SL_from_mapped_reads.txt | grep -v '#' | awk -F'\t' '$4>10' > names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.txt


# Total transcripts: 60559


awk -F'\t' '{print $1"\t0\t"$2}' CCAC0185-TGI-SCAFF.transdecoder.cds.fasta.fai > regions2check.txt

cat regions2check.txt.local.check_misassembly \
 | awk -F'\t' '$4=="Failed_min_ok_reads"' \
 | awk -F'\t' '{ split($5,s,","); for (x in s) {print $1"\t"s[x]"\t"s[x]+1} }' \
 | bedtools sort | bedtools merge | awk '{print $0"\tmisassembly\t0\t*"}' \
> misassemblies.bed


## CHECK: misassemblies.bed
cut -f1 misassemblies.bed | uniq | wc -l
50260
# 60559 total genes


F="names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads"
awk 'NR>1{printf "%s\t%s\t%s\t%s:%s:%.10f\n", $1, $2-1, $2-1, $3, $4, $4/$3}' $F.txt > $F.insertion_points.bed

# Remove sequences that are <100 bp of a mispredicted region. (NOTE some features in -a will appear multiple times in output if two equal distant features found in -b; use sort|uniq to solve this)
bedtools closest -d \
 -a <(cat $F.insertion_points.bed | bedtools sort) \
 -b <(cut -f1-3 misassemblies.bed | bedtools sort) | awk '$6==-1 || $8>100' | cut -f1-4 \
 | sort | uniq > $F.insertion_points.noMissassemblies.bed





# Number of SL addition sites with >10 supporting reads
cut -f1 $F.insertion_points.noMissassemblies.bed | wc -l
4006

# Number of transcripts with SL addition sites with >10 supporting reads
cut -f1 $F.insertion_points.noMissassemblies.bed | uniq | wc -l
3921

# Number of SL addition sites which are positioned >100 bp from the 5-prime end.
awk -F'\t' '$2>100' $F.insertion_points.noMissassemblies.bed | wc -l
63

# Number of CDS with SL addition sites which are positioned >100 bp from the 5-prime end.
awk -F'\t' '$2>100' $F.insertion_points.noMissassemblies.bed | cut -f1 | sort | uniq | wc -l
59








