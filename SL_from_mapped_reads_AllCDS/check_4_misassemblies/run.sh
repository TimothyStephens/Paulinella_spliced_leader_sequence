

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"

~/scripts/grepf_column.py -i ../Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai -f ../names2run.txt | awk -F'\t' '{print $1"\t0\t"$2}' > regions2check.bed

split --numeric-suffixes=1 -a 3 -l 250 regions2check.bed regions2check.bed.part-

ls -1 --color=none regions2check.bed.part-* > files2run.txt


##
## Remove Paulinella_micropora_KR01_nuclear___g18286.t1 - it have >100,000x coverage and it staking days to analyze. 
##
B="../All_combined.local_mapping.coordsorted.bam"
R="regions2check.bed.part-007.removed1"
python check_4_misassemblies_between_two_points.py --region_size 20 --min_ok_reads 10 -b $B -r ${R} -o ${R}.local.check_misassembly --position_info ${R}.local.position_info > ${R}.local.log



## Get all the misassembled positions and convert them into bed features and then merge bed features. 

# Join results files together and remove headers (seqid   start   end     classification  bad_positions). 
cat regions2check.bed.part-{001..006}.local.check_misassembly \
    regions2check.bed.part-007.removed1.local.check_misassembly \
    regions2check.bed.part-{008..052}.local.check_misassembly \
 | awk -F'\t' '$1!="seqid"' \
 > regions2check.bed.local.check_misassembly

cat regions2check.bed.local.check_misassembly \
 | awk -F'\t' '$4=="Failed_min_ok_reads"' \
 | awk -F'\t' '{ split($5,s,","); for (x in s) {print $1"\t"s[x]"\t"s[x]+1} }' \
 | bedtools sort | bedtools merge | awk '{print $0"\tmisassembly\t0\t*"}' \
> misassemblies.bed


