# Scan a genome database using BLAST
# Created by Bryan White, 2015
# input: list of FASTA files
# output: completed BLAST database
#
'''
USAGE
  blastn [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-entrez_query entrez_query]
    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
    [-subject subject_input_file] [-subject_loc range] [-query input_file]
    [-out output_file] [-evalue evalue] [-word_size int_value]
    [-gapopen open_penalty] [-gapextend extend_penalty]
    [-perc_identity float_value] [-qcov_hsp_perc float_value]
    [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
    [-xdrop_gap_final float_value] [-searchsp int_value]
    [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
    [-min_raw_gapped_score int_value] [-template_type type]
    [-template_length int_value] [-dust DUST_options]
    [-filtering_db filtering_database]
    [-window_masker_taxid window_masker_taxid]
    [-window_masker_db window_masker_db] [-soft_masking soft_masking]
    [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
    [-best_hit_score_edge float_value] [-window_size int_value]
    [-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
    [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
    [-outfmt format] [-show_gis] [-num_descriptions int_value]
    [-num_alignments int_value] [-line_length line_length] [-html]
    [-max_target_seqs num_sequences] [-num_threads int_value] [-remote]
    [-version]

DESCRIPTION
   Nucleotide-Nucleotide BLAST 2.2.31+

Use '-help' to print detailed descriptions of command line arguments
'''

import os
import random
import subprocess
import textwrap

## Local libraries
from libs.blst_libs import tstamp

query_file = "query_targets.txt"
db_file = "target_genomes.txt"

output_path = "~/data_analysis/data/genome_assemblies/"
blast_path = "~/data_analysis/apps/ncbi-blast-2.2.31+/bin/"
seq_type = "blastn"

random_scan = 1
r_scan_size = 1000
r_scan_n = 10
merge_contigs = 1

# Load query file names
print(tstamp()+"\tReading queries...")
queries = []
with open(query_file) as inputfile:
        for line in inputfile:
                queries.append(line.strip().split(','))

# Load database file names
print(tstamp()+"\tReading databases...")
databases = []
with open(db_file) as inputfile:
        for line in inputfile:
                databases.append(line.strip().split(','))

# Start main loop of query sequence files
for query_file in queries:
	query_file = ''.join(query_file[0])
	# Begin analysis of query file
	print(tstamp()+"\tLoading query file: "+query_file)

	# Determine query file path
	split_query = query_file.split('/')
	split_query.pop()
	query_path = '/'.join(split_query)
	print(tstamp()+"\tQuery File path: " + query_path)

	# Load sequences into dictionary
	print(tstamp()+"\tLoading sequences...")
	sequences = {}
	qfile = os.path.expanduser(query_file)
	with open(qfile) as inputfiles:
		current_id = ''
		current_seq = ''
		line_i = 0
		for line in inputfiles:
			line = line.strip('\n')
			if(merge_contigs == 0):
				if('>' in line):
					# Set ID
					current_id = line.strip('>')
					print(current_id)
					# Dump previous sequence, start new
					sequences[current_id] = ''
				else:
					# Append sequence to dictionary
					sequences[current_id] += line
			# Only read first ID, thus merging contigs
			if(merge_contigs == 1):
				if('>' in line and line_i == 0):
					# Set ID
					current_id = line.strip('>')
					print(current_id)
					# Dump previous sequence, start new
					sequences[current_id] = ''
				else:
					# Append sequence to dictionary
					sequences[current_id] += line

			line_i = line_i + 1

	print(tstamp()+"\tLoaded "+str(len(sequences))+" sequences")

	#r_scan_sub_file = os.path.expanduser(output_path+"r_scan_contigs.fas")
	r_scan_sub_file = output_path+"r_scan_contigs.fas"
	if(random_scan == 1 and merge_contigs == 1):
		# Access first sequence in dict (should be merged)
		merged_contigs = next (iter (sequences.values()))
		print(tstamp()+"\tGenerating "+str(r_scan_n)+" slices of size "+str(r_scan_size)+"...")
		r_scan_f = open(os.path.expanduser(r_scan_sub_file), 'w')
		for i in range(0,r_scan_n):
			sub_seq = ''
			# Ensure sub sequence size is sufficient (equal to r_scan_size)
			# Extract a substring from the contig
			while True:
				start = random.randint(0, len(merged_contigs))
				end = start+r_scan_size
				sub_seq = merged_contigs[start:end]
				if(len(sub_seq) == r_scan_size):
					break
			r_scan_f.write(">Sub" + str(i) + "\n" + sub_seq + "\n")
		print(tstamp()+"\tPrinted "+str(r_scan_n)+" random sub-sequences to "+r_scan_sub_file)
		r_scan_f.close()
		# Set query file to sub sequence file
		query_file = r_scan_sub_file

	for database in databases:
		database = ''.join(database[0])
		print(tstamp()+"\tTargeting "+database+" with ")
		# Form blastn command
		scandb_command = blast_path + seq_type + " -query " + query_file + " -db " + ''.join(database) +\
			" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\"" +\
			" -out test.tsv"
		print(tstamp()+"\tRunning: ")
		print(scandb_command)
		#sts = subprocess.Popen(scandb_command, shell=True).wait()
		os.system(scandb_command)


'''
 *** Formatting options
 -outfmt <String>
   alignment view options:
     0 = pairwise,
     1 = query-anchored showing identities,
     2 = query-anchored no identities,
     3 = flat query-anchored, show identities,
     4 = flat query-anchored, no identities,
     5 = XML Blast output,
     6 = tabular,
     7 = tabular with comment lines,
     8 = Text ASN.1,
     9 = Binary ASN.1,
    10 = Comma-separated values,
    11 = BLAST archive format (ASN.1),
    12 = JSON Seqalign output,
    13 = JSON Blast output,
    14 = XML2 Blast output

   Options 6, 7, and 10 can be additionally configured to produce
   a custom format specified by space delimited format specifiers.
   The supported format specifiers are:
            qseqid means Query Seq-id
               qgi means Query GI
              qacc means Query accesion
           qaccver means Query accesion.version
              qlen means Query sequence length
            sseqid means Subject Seq-id
         sallseqid means All subject Seq-id(s), separated by a ';'
               sgi means Subject GI
            sallgi means All subject GIs
              sacc means Subject accession
           saccver means Subject accession.version
           sallacc means All subject accessions
              slen means Subject sequence length
            qstart means Start of alignment in query
              qend means End of alignment in query
            sstart means Start of alignment in subject
              send means End of alignment in subject
              qseq means Aligned part of query sequence
              sseq means Aligned part of subject sequence
            evalue means Expect value
          bitscore means Bit score
             score means Raw score
            length means Alignment length
            pident means Percentage of identical matches
            nident means Number of identical matches
          mismatch means Number of mismatches
          positive means Number of positive-scoring matches
           gapopen means Number of gap openings
              gaps means Total number of gaps
              ppos means Percentage of positive-scoring matches
            frames means Query and subject frames separated by a '/'
            qframe means Query frame
            sframe means Subject frame
              btop means Blast traceback operations (BTOP)
           staxids means unique Subject Taxonomy ID(s), separated by a ';'
                         (in numerical order)
         sscinames means unique Subject Scientific Name(s), separated by a ';'
         scomnames means unique Subject Common Name(s), separated by a ';'
        sblastnames means unique Subject Blast Name(s), separated by a ';'
                         (in alphabetical order)
        sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
                         (in alphabetical order)
            stitle means Subject Title
        salltitles means All Subject Title(s), separated by a '<>'
           sstrand means Subject Strand
             qcovs means Query Coverage Per Subject
           qcovhsp means Query Coverage Per HSP
   When not provided, the default value is:
   'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
   evalue bitscore', which is equivalent to the keyword 'std'
   Default = `0'
'''
