# Scan a genome database using BLAST
# Created by Bryan White, 2015
# input: list of FASTA files
# output: completed BLAST database
#

# Copyright (c) 2015-2016 Bryan White, bpcwhite@gmail.com

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import random
import subprocess
import textwrap

## Local libraries
from blst_libs import tistamp, load_sequences, rand_subseq

query_file = "query_targets.txt"
db_file = "target_genomes.txt"

output_path = "~/data_analysis/data/genome_assemblies/"
blast_path = "~/data_analysis/apps/ncbi-blast-2.2.31+/bin/"
seq_type = "blastn"

random_scan = 1
blast_scan = 0
mem_scan = 1
r_scan_size = 500
r_scan_n = 10
merge_contigs = 1
seq_line_limit = 100

'''
-word_size <Integer, >=4>
   Word size for wordfinder algorithm (length of best perfect match)
 -gapopen <Integer>
   Cost to open a gap
 -gapextend <Integer>
   Cost to extend a gap
 -penalty <Integer, <=0>
   Penalty for a nucleotide mismatch
 -reward <Integer, >=0>
   Reward for a nucleotide match
'''

# Load query file names
print(tistamp(1)+"\tReading queries...")
queries = []
with open(query_file) as inputfile:
        for line in inputfile:
                queries.append(line.strip().split(','))

# Load database file names
print(tistamp(1)+"\tReading databases...")
databases = []
with open(db_file) as inputfile:
        for line in inputfile:
                databases.append(line.strip().split(','))

# Start main loop of query sequence files
for query_file in queries:
	query_file = ''.join(query_file[0])
	# Begin analysis of query file
	print(tistamp(1)+"\tLoading query file: "+query_file)

	# Determine query file path
	split_query = query_file.split('/')
	query_file_name = split_query[-1]
	split_query.pop()
	query_path = '/'.join(split_query)
	print(tistamp(1)+"\tQuery File path: " + query_path)

	# Load sequences into dictionary
	print(tistamp(1)+"\tLoading sequences...")
	query_sequences = {}
	query_sequences = load_sequences(query_file, merge_contigs, seq_line_limit)
	print(tistamp(1)+"\tLoaded "+str(len(query_sequences))+" sequences")

	#r_scan_sub_file = os.path.expanduser(output_path+"r_scan_contigs.fas")
	r_scan_sub_file = query_path+"/r_scan_contigs.fas"
	scan_file = query_file # Actual file that will be blasted against the database
	if(random_scan == 1 and merge_contigs == 1 and blast_scan == 1):
		# Access first sequence in dict (should be merged)
		query_sequence = next (iter (query_sequences.values()))
		print(tistamp(1)+"\tGenerating "+str(r_scan_n)+" slices of size "+str(r_scan_size)+"...")
		r_scan_f = open(os.path.expanduser(r_scan_sub_file), 'w')
		for i in range(0,r_scan_n):
			sub_seq = rand_subseq(query_sequence, r_scan_size)
			r_scan_f.write(">Sub" + str(i) + "\n" + sub_seq + "\n")
		print(tistamp(1)+"\tPrinted "+str(r_scan_n)+" random sub-sequences to "+r_scan_sub_file)
		r_scan_f.close()
		# Set query file to sub sequence file
		scan_file = r_scan_sub_file

		for database in databases:
			database = ''.join(database[0])
			# Determine database name
			split_database = database.split('/')
			database_file_name = split_database[-1]
			split_database.pop()
			database_path = '/'.join(split_database)

			# format output path
			print(tistamp(1)+"\tTargeting "+database+" with ")
			blst_results_file = query_path+"/qry_"+query_file_name+"_db_"+database_file_name+"_results_"+tistamp(2)+".tsv"
			print(tistamp(1)+"\tOutputing BLAST results to "+blst_results_file)
			# Form blastn command
			scandb_command = blast_path + seq_type + " -query " + scan_file + " -db " + ''.join(database) +\
				" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\"" +\
				" -out " + blst_results_file
			print(tistamp(1)+"\tRunning: ")
			print(tistamp(1)+"\t"+scandb_command)
			#sts = subprocess.Popen(scandb_command, shell=True).wait()
			os.system(scandb_command)



	if(random_scan == 1 and mem_scan == 1 and merge_contigs == 1):
		# Access first sequence in dict (should be merged)
		query_sequence = next (iter (query_sequences.values()))
		print(tistamp(1)+"\tGenerating "+str(r_scan_n)+" slices of size "+str(r_scan_size)+"...")

		# Generate random query sequences
		rand_query_seqs = {}
		for i in range(0,r_scan_n):
			(sub_seq,start,end) = rand_subseq(query_sequence, r_scan_size)
			midpoint = start + int(start-end)/2
			rand_query_seqs[i] = sub_seq

		print(tistamp(1)+"\tCreated "+str(r_scan_n)+" random sub-sequences.")

		for database in databases:
			database = ''.join(database[0])
			# Determine database name
			split_database = database.split('/')
			database_file_name = split_database[-1]
			split_database.pop()
			database_path = '/'.join(split_database)
			print(database)
			# Load database sequences
			db_sequences = {}
			db_sequences = load_sequences(database, merge_contigs, seq_line_limit)
			db_sequence = next (iter (db_sequences.values()))

			rand_db_seqs = {}
			for i in range(0,r_scan_n):
				(sub_seq,start,end) = rand_subseq(db_sequence, r_scan_size)
				rand_db_seqs[i] = sub_seq

		''' DIST '''
		'''
		ranmem
		'''
