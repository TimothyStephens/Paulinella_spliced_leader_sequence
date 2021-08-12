
export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"

SL="Paulinella_SL_sequence.fa"



##
##
## KR01
##
##
DB="Paulinella_micropora_KR01_nuclear.assembly.fasta"
blastn -query $SL -db $DB -out $DB.blastn-short_SL.outfmt6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -evalue 1000 -task blastn-short
awk -F'\t' '$8==20' $DB.blastn-short_SL.outfmt6 | less

## Get "full" (from pos 1-20) and "partial" (from pos 20 to > pos 1; i.e., missing start)
awk -F'\t' '$7==1 && $8==20' $DB.blastn-short_SL.outfmt6 > $DB.blastn-short_SL.outfmt6.full_SL
awk -F'\t' '$7!=1 && $8==20' $DB.blastn-short_SL.outfmt6 > $DB.blastn-short_SL.outfmt6.partial_SL_start


## 10 bp upstream + 10 bp downstream
seqkit subseq --up-stream 10 --down-stream 10 --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"$9-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.full_SL) \
 $DB > $DB.blastn-short_SL.outfmt6.full_SL.10up_10down.fa

cat $DB.blastn-short_SL.outfmt6.full_SL.10up_10down.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==40' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.full_SL.10up_10down.cleaned.fa


## 30 bp upstream + 10 bp downstream
seqkit subseq --up-stream 30 --down-stream 10 --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"$9-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.full_SL) \
 $DB > $DB.blastn-short_SL.outfmt6.full_SL.30up_10down.fa

cat $DB.blastn-short_SL.outfmt6.full_SL.30up_10down.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==60' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.full_SL.30up_10down.cleaned.fa


## Partial SL sequences - no upstream or downstream as its not conserved and hinders identification of other SL patterns
seqkit subseq --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"($9-($7-1))-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9+($7-1)"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.partial_SL_start) \
 $DB > $DB.blastn-short_SL.outfmt6.partial_SL_start.fa

cat $DB.blastn-short_SL.outfmt6.partial_SL_start.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==20' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa


##
##
## Chromatophora
##
##
DB="Paulinella_chromatophora_CCAC0185_nuclear.assembly.fasta"
blastn -query $SL -db $DB -out $DB.blastn-short_SL.outfmt6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -evalue 1000 -task blastn-short
awk -F'\t' '$8==20' $DB.blastn-short_SL.outfmt6 | less

## Get "full" (from pos 2-20) and "partial" (from pos 20 to > pos 1; i.e., missing start)
awk -F'\t' '$7<=2 && $8==20' $DB.blastn-short_SL.outfmt6 > $DB.blastn-short_SL.outfmt6.full_SL
awk -F'\t' '$7 >2 && $8==20' $DB.blastn-short_SL.outfmt6 > $DB.blastn-short_SL.outfmt6.partial_SL_start


## 10 bp upstream + 10 bp downstream
seqkit subseq --up-stream 10 --down-stream 10 --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"$9-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.full_SL) \
 $DB > $DB.blastn-short_SL.outfmt6.full_SL.10up_10down.fa

cat $DB.blastn-short_SL.outfmt6.full_SL.10up_10down.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==39' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.full_SL.10up_10down.cleaned.fa


## 30 bp upstream + 10 bp downstream
seqkit subseq --up-stream 30 --down-stream 10 --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"$9-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.full_SL) \
 $DB > $DB.blastn-short_SL.outfmt6.full_SL.30up_10down.fa

cat $DB.blastn-short_SL.outfmt6.full_SL.30up_10down.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==59' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.full_SL.30up_10down.cleaned.fa


## Partial SL sequences - no upstream or downstream as its not conserved and hinders identification of other SL patterns
seqkit subseq --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"($9-($7-1))-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9+($7-1)"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.partial_SL_start) \
 $DB > $DB.blastn-short_SL.outfmt6.partial_SL_start.fa

cat $DB.blastn-short_SL.outfmt6.partial_SL_start.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==20' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa






##
##
## Ovalis
##
##
DB="Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta"
blastn -query $SL -db $DB -out $DB.blastn-short_SL.outfmt6 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -evalue 1000 -task blastn-short

## Get "full" (from pos 1-20) and "partial" (from pos 20 to > pos 1; i.e., missing start)
awk -F'\t' '$7==1 && $8==20' $DB.blastn-short_SL.outfmt6 > $DB.blastn-short_SL.outfmt6.full_SL
awk -F'\t' '$7!=1 && $8==20' $DB.blastn-short_SL.outfmt6 > $DB.blastn-short_SL.outfmt6.partial_SL_start

## 10 bp upstream + 10 bp downstream
#seqkit subseq --up-stream 10 --down-stream 10 --bed \
# <(awk -F'\t' '{ if($10>$9){print $2"\t"$9-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.full_SL) \
# $DB > $DB.blastn-short_SL.outfmt6.full_SL.10bp_up-down.fa


## Partial SL sequences - no upstream or downstream as its not conserved and hinders identification of other SL patterns
seqkit subseq --bed \
 <(awk -F'\t' '{ if($10>$9){print $2"\t"($9-($7-1))-1"\t"$10"\tSL\t.\t+"} else{print $2"\t"$10-1"\t"$9+($7-1)"\tSL\t.\t-"} }' $DB.blastn-short_SL.outfmt6.partial_SL_start) \
 $DB > $DB.blastn-short_SL.outfmt6.partial_SL_start.fa

cat $DB.blastn-short_SL.outfmt6.partial_SL_start.fa \
 | sed -e 's/:.*//' \
 | seqkit fx2tab | awk -F'\t' 'length($2)==20' | seqkit tab2fx \
 > $DB.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa

## Build a tree using the SL seq alignment
/home/timothy/programs/iqtree-1.6.12-Linux/bin/iqtree -s Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa -m JC

## Add SL alignment to iTOL annotation file
cat iTOL_alignment_annotation.txt > $DB.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa.iTOL_alignment_annotation.txt
cat $DB.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa >> $DB.blastn-short_SL.outfmt6.partial_SL_start.cleaned.fa.iTOL_alignment_annotation.txt




