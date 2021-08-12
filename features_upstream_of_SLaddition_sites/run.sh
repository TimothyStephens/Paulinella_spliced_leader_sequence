


export PATH=$PATH:/home/timothy/programs/bedtools-2.29.2/bin

IP="../../../../../Paulinella_micropora_KR01/GeneAnnotation/Paulinella_micropora_KR01_Nuclear_2021_04/InterProScan/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa.interproscan.gff3"

## Get *just* CDS coords of SL-addition sites.
cut -f1-3 names2run.txt.filtered.SL_from_mapped_reads.gt10SLreads.insertion_points.noMissassemblies.bed | bedtools sort > SL_from_mapped_reads.gt10SLreads.cds_coords.bed






correlate_SL_with_features () {
	PREFIX=$1
	F="$PREFIX.SL_addition_sites_intersection"
	
	## Sort domains and merge overlapping features (join domain ids into list)
	##   - seqid(\t)start(\t)end(\t)domain_id_list[comma sep]
	bedtools sort -i $PREFIX.pep_coords.bed | bedtools merge -c 4 -o collapse > $PREFIX.merged.pep_coords.bed
	awk -F'\t' '{print $1"\t"($2*3)"\t"($3*3)"\t"$4}' $PREFIX.merged.pep_coords.bed > $PREFIX.merged.cds_coords.bed
	
	## Get SL-addition sites which are upstream OR overlap with a domain feature
	bedtools closest -t first -D a -id -io -a SL_from_mapped_reads.gt10SLreads.cds_coords.bed -b $PREFIX.merged.cds_coords.bed | awk -F'\t' '$5!=-1' > "$F.upstream.bed"
	bedtools closest -t first -D a -a SL_from_mapped_reads.gt10SLreads.cds_coords.bed -b $PREFIX.merged.cds_coords.bed | awk -F'\t' '$5!=-1' | awk -F'\t' '$2>=$5 && $2<=$6' > "$F.overlap.bed"
	
	## Get SL-addition sites that dont have detectable mispredictions between them and the most upstream part of the freature.
	##   - seqid(\t)region_start(\t)region_end(\t)SL_seqid(\t)SL_start(\t)SL_end
	awk -F'\t' '{print $1"\t"$5"\t"$2"\t"$1"\t"$2"\t"$3}' "$F.upstream.bed" | bedtools intersect -v -a - -b misassemblies.bed > "$F.upstream.without_misassemblies.bed"
	
	## Get SL-addition sites that sit 60 bp (20 aa) from the start of the overlapping feature.
	awk -F'\t' '$2>($5+60)' "$F.overlap.bed" > "$F.overlap.60bp_from_domain_start.bed"
	awk -F'\t' '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' "$F.overlap.60bp_from_domain_start.bed" | bedtools intersect -v -a - -b misassemblies.bed > "$F.overlap.60bp_from_domain_start.without_misassemblies.bed"
}











## Get just the 'domain' results from InterProScan - filter domains by e-value < 1e-5 (or just return all results for ProSiteProfiles as it only has a bitscore)
##
## seqid(\t)start(\t)end(\t)domain_id

cat $IP \
 | awk -F'\t' '$2=="CDD" || $2=="PANTHER" || $2=="Pfam" || $2=="PRINTS" || $2=="SMART" || $2=="SUPERFAMILY" || $2=="TIGRFAM"' \
 | awk -F'\t' '$6<1e-5 {print $1"\t"$4-1"\t"$5"\t"$9}' \
 | sed -e 's/\([^\t]*\)\t\([^\t]*\)\t\([^\t]*\)\t.*;Name=\([^;]*\).*/\1\t\2\t\3\t\4/' \
 > interproscan_domains.pep_coords.bed

cat $IP \
 | awk -F'\t' '$2=="ProSiteProfiles"' \
 | awk -F'\t' '{print $1"\t"$4-1"\t"$5"\t"$9}' \
 | sed -e 's/\([^\t]*\)\t\([^\t]*\)\t\([^\t]*\)\t.*;Name=\([^;]*\).*/\1\t\2\t\3\t\4/' \
 >> interproscan_domains.pep_coords.bed

correlate_SL_with_features "interproscan_domains"

wc -l interproscan_domains.SL_addition_sites_intersection.*.without_misassemblies.bed
#  358 interproscan_domains.SL_addition_sites_intersection.overlap.60bp_from_domain_start.without_misassemblies.bed
#  137 interproscan_domains.SL_addition_sites_intersection.upstream.without_misassemblies.bed
#  495 total








## Get 'transmembrane' domain regions from InterProScan results. 
##
## seqid(\t)start(\t)end(\t)transmembrane_id
cat $IP \
 | awk -F'\t' '$2=="Phobius" || $2=="TMHMM"' \
 | awk -F'\t' '$9~"TRANSMEMBRANE" || $9~"TMhelix" {print $1"\t"$4-1"\t"$5"\t"$9}' \
 | sed -e 's/\([^\t]*\)\t\([^\t]*\)\t\([^\t]*\)\t.*;Name=\([^;]*\).*/\1\t\2\t\3\t\4/' \
 > interproscan_transmembrane.pep_coords.bed

correlate_SL_with_features "interproscan_transmembrane"

wc -l interproscan_transmembrane.SL_addition_sites_intersection.*.without_misassemblies.bed
#   26 interproscan_transmembrane.SL_addition_sites_intersection.overlap.60bp_from_domain_start.without_misassemblies.bed
#  297 interproscan_transmembrane.SL_addition_sites_intersection.upstream.without_misassemblies.bed
#  323 total








## Get 'coil' regions from InterProScan results.
##
## seqid(\t)start(\t)end(\t)coil_id
cat $IP \
 | awk -F'\t' '$2=="Coils"' \
 | awk -F'\t' '{print $1"\t"$4-1"\t"$5"\t"$9}' \
 | sed -e 's/\([^\t]*\)\t\([^\t]*\)\t\([^\t]*\)\t.*;Name=\([^;]*\).*/\1\t\2\t\3\t\4/' \
 > interproscan_coils.pep_coords.bed

correlate_SL_with_features "interproscan_coils"

wc -l interproscan_coils.SL_addition_sites_intersection.*.without_misassemblies.bed
#    6 interproscan_coils.SL_addition_sites_intersection.overlap.60bp_from_domain_start.without_misassemblies.bed
#  128 interproscan_coils.SL_addition_sites_intersection.upstream.without_misassemblies.bed
#  134 total








cut -f1 SL_from_mapped_reads.gt10SLreads.cds_coords.bed | sort | wc -l
# 8423
cut -f1 SL_from_mapped_reads.gt10SLreads.cds_coords.bed | sort | uniq | wc -l
# 6579



##
## Number of SL addition sites
##
cut -f4-6 interproscan_domains.*.without_misassemblies.bed | sort | uniq | wc -l
# 471
cut -f4-6 interproscan_domains.*.without_misassemblies.bed interproscan_transmembrane.*.without_misassemblies.bed | sort | uniq | wc -l
# 681
cut -f4-6 interproscan_domains.*.without_misassemblies.bed interproscan_transmembrane.*.without_misassemblies.bed interproscan_coils.*.without_misassemblies.bed | sort | uniq | wc -l
# 750

# 471 / 8423 = 5.59%
# 681 / 8423 = 8.09%
# 750 / 8423 = 8.90%


##
## Number of genes with SL addition sites
##
cut -f4 interproscan_domains.*.without_misassemblies.bed | sort | uniq | wc -l
# 406
cut -f4 interproscan_domains.*.without_misassemblies.bed interproscan_transmembrane.*.without_misassemblies.bed | sort | uniq | wc -l
# 572
cut -f4 interproscan_domains.*.without_misassemblies.bed interproscan_transmembrane.*.without_misassemblies.bed interproscan_coils.*.without_misassemblies.bed | sort | uniq | wc -l
# 631

# 406 / 6579 = 6.17%
# 572 / 6579 = 8.69%
# 631 / 6579 = 9.59%






