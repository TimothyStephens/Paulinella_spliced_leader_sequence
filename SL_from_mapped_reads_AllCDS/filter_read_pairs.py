#!/usr/bin/env python2
DESCRIPTION = '''
Takes a bam file and filters aligned reads. Read must:
	- be aligned
	- have aligned mate
	- have <=5bp 3-prime soft-clipping
	- have aligned mate with <=5bp 3-prime soft-clipping
	- have <=5bp mismatches or indels 

NOTE:
	- Assumes read aligner outputs edit distance in the optional fields section of each read (tag: NM in bowtie2 output)

pysam documentation: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
'''
import sys
import os
import argparse
import logging
import pysam
import gzip

## Pass arguments.
def main():
	# Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-b', '--bam', metavar='aligned_reads.bam', 
		required=True, type=argparse.FileType('r'), 
		help='Aligned BAM file sorted by read name [not by coords] (required)'
	)
	parser.add_argument('-o', '--out', metavar='filtered_reads.bam', 
		required=False, default=sys.stdout, type=argparse.FileType('w'), 
		help='Out BAM file with filtered reads (default: %(default)s)'
	)
	parser.add_argument('--softclipping',
		required=False, default=5, type=int,
		help='Max number of soft-clipped bases allowed at 3-prime ends of read pairs (default: %(default)s)'
	)
	parser.add_argument('--edit_dist',
		required=False, default=5, type=int,
		help='Max number of mismatch or indel bases to allow in a read (default: %(default)s)'
	)
	parser.add_argument('--debug', action='store_true', 
		required=False, 
		help='Print DEBUG info (default: %(default)s)'
	)
	args = parser.parse_args()
	
	# Set up basic debugger
	logFormat = "[%(levelname)s]: %(message)s"
	logging.basicConfig(format=logFormat, stream=sys.stderr, level=logging.INFO)
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	
	logging.debug('%s', args) ## DEBUG
	
	filter_read_pairs(args.bam, args.out, args.softclipping, args.edit_dist)
	
	
def filter_read_pairs(bam_in_fh, bam_out_fh, softclipping, edit_dist):
	bam_in =  pysam.AlignmentFile(bam_in_fh, "rb", require_index=False)
	bam_out = pysam.AlignmentFile(bam_out_fh, "wb", template=bam_in)
	
	# FOR each read aligned to seqid
	for read_1 in bam_in:
		#logging.debug('Found first: %s', read_1) ## DEBUG
		
		# Skip this read if it doesnt have a properly paired mate. 
		if not read_1.is_proper_pair:
			continue
		
		read_2 = next(bam_in)
		#logging.debug('Found second: %s', read_2) ## DEBUG
		logging.debug('Pair found: %s <-> %s', read_1.query_name, read_2.query_name) ## DEBUG
		
		if read_1.query_name != read_2.query_name:
			logging.error('Pair found that dont have matching names: %s <-> %s', read_1.query_name, read_2.query_name) ## ERROR
			logging.error('Was the input BAM file sorted by coords and not by read names?') ## ERROR
			sys.exit(1)
		
		cigar_1 = read_1.cigartuples
		cigar_2 = read_2.cigartuples
		logging.debug('CIGAR: %s <-> %s', cigar_1, cigar_2) ## DEBUG
		
		## skip this pair if either read has > X number of softclipped bases at 3-prime ends
		if (cigar_1[-1][0] == 4 and cigar_1[-1][1] > softclipping) or (cigar_2[-1][0] == 4 and cigar_2[-1][1] > softclipping):
			logging.debug('  - Pair rejected because of soft clipping.') ## DEBUG
			continue
		
		## Cant use CIGAR because if aligner uses 'M' then we dont know if base is match or mismatch.
		#editdist_1 = sum([y for x, y in cigar_1 if x == 1 or x == 2 or x == 8])
		#editdist_2 = sum([y for x, y in cigar_2 if x == 1 or x == 2 or x == 8])
		editdist_1 = read_1.get_tag('NM')
		editdist_2 = read_2.get_tag('NM')
		
		## skip this pair if either read has edit distance > X (indels+mismatches)
		if editdist_1 > edit_dist or editdist_2 > edit_dist:
			logging.debug('  - Pair rejected because of high edit distance.') ## DEBUG
			continue
		
		## Write filtered reads to output file
		bam_out.write(read_1)
		bam_out.write(read_2)
	bam_in.close()
	bam_out.close()


def __parse_file_check_compression(fh, mode='r'):
	'''
	Open stdin/normal/gzip files - check file exists (if mode='r') and open using appropriate function.
	'''
	# Check file exists if mode='r'
	if not os.path.exists(fh) and mode == 'r':
		raise argparse.ArgumentTypeError("The file %s does not exist!" % fh)
	
	## open with gzip if it has the *.gz extension, else open normally (including stdin)
	try:
		if fh.endswith(".gz"):
			return gzip.open(fh, mode+'b')
		else:
			return open(fh, mode)
	except IOError as e:
		raise argparse.ArgumentTypeError('%s' % e)



if __name__ == '__main__':
	main()
