

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"

awk -F'\t' 'NR==1 || $4>10' names2run.txt.SL_from_mapped_reads.txt > names2run.txt.SL_from_mapped_reads.gt10SLreads.txt
awk -F'\t' 'NR==1 || $4>10' names2run.txt.filtered.SL_from_mapped_reads.txt > names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.txt





F="names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads"

# Total SL positions
awk -F'\t' 'NR>1' $F.txt | wc -l
16037


# Seqs with SL positions
awk -F'\t' 'NR>1 {print $1}' $F.txt | sort | uniq | wc -l
11618



## CHECK: misassemblies.bed
cut -f1 misassemblies.bed | uniq | wc -l
9513
# 32361 total genes
# Had to ignore 'Paulinella_micropora_KR01_nuclear___g18286.t1' as it had such high read cov. (>1000000x)
# that its was taking days to check it using my script.


# Get SL-truncation/insertion points
# Conver 1-based pos to 0-based position between bases (start and end are the same to represent an insertion point between bases)
#         |A|T G C A T G C A|T|G C A T G C A T G|C|A T G C A T G C A|T|G C A T G C A T G C
#         |||               |||                 |||                 |||
# 0:      0|1               9|10               19|20               29|30
# 0 pos:   0                 9                   19                  29
# 1 pos:   1                 10                  20                  30
awk 'NR>1{printf "%s\t%s\t%s\t%s:%s:%.10f\n", $1, $2-1, $2-1, $3, $4, $4/$3}' $F.txt > $F.insertion_points.bed

# Remove sequences that are <100 bp of a mispredicted region. (NOTE some features in -a will appear multiple times in output if two equal distant features found in -b; use sort|uniq to solve this)
bedtools closest -d \
 -a <(cat $F.insertion_points.bed | bedtools sort) \
 -b <(cut -f1-3 misassemblies.bed | bedtools sort) | awk '$6==-1 || $8>100' | cut -f1-4 \
 | sort | uniq > $F.insertion_points.noMissassemblies.bed

## Manually check, using IGV, which features from Paulinella_micropora_KR01_nuclear___g18286.t1 to remove
# Paulinella_micropora_KR01_nuclear___g18286.t1         104     1209786 2575   - remove
# Paulinella_micropora_KR01_nuclear___g18286.t1         108     1320541 11     - remove
# Paulinella_micropora_KR01_nuclear___g18286.t1         3572    1606545 16     - keep






##
## SL positions not near (100 bp) mispredicted regions
##

# Total SL positions - not near mispredictions
awk -F'\t' '{print $1}' $F.insertion_points.noMissassemblies.bed | wc -l
8423


# Seqs with SL positions - not near mispredictions
awk -F'\t' '{print $1}' $F.insertion_points.noMissassemblies.bed | sort | uniq | wc -l
6579


# Total SL positions that are exactly at the 5-prime end of the CDS
awk '$2==0' $F.insertion_points.noMissassemblies.bed | wc -l
4522


# Average proportion of SL reads at each SL-trunction point
awk -F'\t' '{split($4,a,":"); printf "%.8f\n", a[3]}' $F.insertion_points.noMissassemblies.bed | sort | ~/scripts/stdin_Avg.py
0.354691635882 # 35.47%


# Average position along the CDS of SL-truncation points
cat $F.insertion_points.noMissassemblies.bed | cut -f2 | ~/scripts/stdin_Avg.py
176.880802564


# Average position (as a proportion of CDS length) of SL-truncation points
~/scripts/add_value_to_table.py -i $F.insertion_points.noMissassemblies.bed -a Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai | awk -F'\t' '{printf "%.8f\n", $2/$5}' | ~/scripts/stdin_Avg.py
0.105541344355 # 10.55%


# Average position (as a proportion of CDS length) of SL-truncation points with ABOVE average read support
cat $F.insertion_points.noMissassemblies.bed | awk -F'\t' '{ split($4,a,":"); if(a[3]>0.354691635882) {print} }' | ~/scripts/add_value_to_table.py -a Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai | awk -F'\t' '{printf "%.8f\n", $2/$5}' | ~/scripts/stdin_Avg.py
0.0253284258472 # 2.53%


# Average position (as a proportion of CDS length) of SL-truncation points with BELOW average read support
cat $F.insertion_points.noMissassemblies.bed | awk -F'\t' '{ split($4,a,":"); if(a[3]<0.354691635882) {print} }' | ~/scripts/add_value_to_table.py -a Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai | awk -F'\t' '{printf "%.8f\n", $2/$5}' | ~/scripts/stdin_Avg.py
0.194872569533 # 19.49






##
## Overlap between SL addition sites and exon splice sites
##

# Get sites not at start of CDS
awk -F'\t' '$2!=0' $F.insertion_points.noMissassemblies.bed | bedtools sort > $F.insertion_points.noMissassemblies.noStart.bed

# For each SL-addition site get closes (if one exists [single exon genes]) splice site
bedtools closest -d -a $F.insertion_points.noMissassemblies.noStart.bed -b CDS_splice_sites.bed > $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed
# Output: A_seq A_start A_end A_info B_seq B_start B_end dist_btw_A_and_B

cat $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed | wc -l
3906 # 3901 total (5 duplicates [are equal distance between two splice sites])






# Total SL-adition sites not including sites at the start [pos 0] of CDS
cat $F.insertion_points.noMissassemblies.noStart.bed | wc -l
3901

# Total lines returned by bedtools closest (5 duplicates found - equal distance between splice sites)
cat $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed | wc -l
3906

# NOTE: cant use distance in column 8 - positions that are either the same or ofset by 1 are given a distance of 0

# No. SL addition sites that sit exactly at a splice site
awk -F'\t' '$2==$6' $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed | wc -l
2237
# 2237/3901 = 0.5734427 (57.34%)

# No. SL addition sites that are on single-exon genes (i.e. have no valid splice sites in -b)
awk -F'\t' '$6==-1' $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed | wc -l
96
# 96/3901 = 0.024609 (2.46%)




# Seq LOGO for SL-additions sites

# Remove SL-addition sites that are too close to the ends of the CDS (wont have enough seq for the logo)
~/scripts/add_value_to_table.py \
 -i $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed \
 -a Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai \
 | awk -F'\t' '$2>=5 && $3<=($9-5)' \
 > $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed.good4plotting
cat $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed.good4plotting | wc -l
3532 # 3901 total (includes 5 dups)


# Get SL addition sites that DO NOT start exactly at the same position as a splice site [includes single exon genes].
# Remove 5 duplicates at this stage.
awk -F'\t' '$2!=$6' $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed.good4plotting | cut -f1-4 | uniq > SL_addition_sites_not_on_splice_sites.bed
awk -F'\t' '{print $1"\t"$2-5"\t"$3+5"\t"$4}' SL_addition_sites_not_on_splice_sites.bed > SL_addition_sites_not_on_splice_sites.4logo.bed

awk -F'\t' '$2==$6' $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed.good4plotting | cut -f1-4 > SL_addition_sites_on_splice_site.bed
awk -F'\t' '{print $1"\t"$2-5"\t"$3+5"\t"$4}' SL_addition_sites_on_splice_site.bed > SL_addition_sites_on_splice_site.4logo.bed

CDS=Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna
seqkit subseq --bed SL_addition_sites_not_on_splice_sites.4logo.bed $CDS > SL_addition_sites_not_on_splice_sites.4logo.fa
seqkit subseq --bed SL_addition_sites_on_splice_site.4logo.bed $CDS > SL_addition_sites_on_splice_site.4logo.fa



cat SL_addition_sites_not_on_splice_sites.bed | wc -l
1358
grep -c '>' SL_addition_sites_not_on_splice_sites.4logo.fa
1358

cat SL_addition_sites_on_splice_site.bed | wc -l
2169
grep -c '>' SL_addition_sites_on_splice_site.4logo.fa
2169










##
## Find reads that encode SL sequences and take just those with ref positions that had to be adjusted (i.e. SL-seq was in soft-clipped region and ref position had to be "moved" to first mapped position).
##

zcat names2run.txt.filtered.SL_from_mapped_reads.SL_read_info.txt.gz \
 | awk -F'\t' 'NR>1 && $13>0 {print $4"\t"$11}' \
 | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' \
 | sort -k1,1 -k2,2n \
> SL_addition_sites_from_adjusted_leader_positions.txt


# No. SL addition sites that sit exactly at a splice site
awk '$2==$6' $F.insertion_points.noMissassemblies.noStart.dist_to_splice_sites.bed | sort -k1,1 -k2,2n > SL_addition_sites_at_exon_boundaries.bed
cat SL_addition_sites_at_exon_boundaries.bed | wc -l
2237
# 2237/3901 = 0.5734427 (57.34%)

# No of SL-addition sites with reads that had to be adjusted that are at splice sites which overlap with SL addition sites
bedtools intersect -wa -wb -a <(cut -f1-4 SL_addition_sites_at_exon_boundaries.bed | bedtools sort) -b <(awk -F'\t' '$2>0{print $1"\t"$2"\t"$2"\t"$3}' SL_addition_sites_from_adjusted_leader_positions.txt | bedtools sort) | awk '$2==$6' | wc -l
203
# 203/2237 = 0.0907465 = 9.07% [90.93 of SL addition sites that overlap with splice sites have ALL reads supporting that position (i.e. no reads supporting additon site along intron)]
# 2237-203 = 2034; 2034/3901 = 52.14% of SL addition sites that dont start at position 1










