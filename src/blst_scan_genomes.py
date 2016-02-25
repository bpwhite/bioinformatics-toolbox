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
import distance
import time
import argparse
import pickle

## Local libraries
from blst_libs import tistamp, load_sequences, rand_subseq, sl_window, similar

start_time = time.time()

query_file = "query_targets.txt"
db_file = "target_genomes.txt"

output_path = "~/data_analysis/data/genome_assemblies/"
blast_path = "~/data_analysis/apps/ncbi-blast-2.2.31+/bin/"
seq_type = "blastn"
'''
parser = argparse.ArgumentParser(description='Scan genomes')
parser.add_argument('--sum', dest='accumulate',
                   help='sum the integers (default: find the max)')
'''
#args = parser.parse_args()

random_scan = 1
random_pairwise = 1
blast_scan = 0
mem_scan = 1
r_scan_size = 1000
sl_window_size = 1000
r_scan_n = 5000
merge_contigs = 1
seq_line_limit = 9999999999
sl_increment = int(0.75*sl_window_size)
use_fastcomp = 0
dist_pcnt = 0.30

sliding_window = 1
expand_size = 10
#pickle_cache = 'build'
clear_cache = 0

#pickle_cache = 1

print("random_scan\t" 		+ str(random_scan))
print("blast_scan\t" 		+ str(blast_scan))
print("mem_scan\t" 			+ str(mem_scan))
print("r_scan_size\t" 		+ str(r_scan_size))
print("r_scan_n\t" 			+ str(r_scan_n))
print("merge_contigs\t" 	+ str(merge_contigs))
print("seq_line_limit\t" 	+ str(seq_line_limit))
print("sl_window_size\t" 	+ str(sl_window_size))
print("sl_increment\t" 		+ str(sl_increment))
print("Cache: "				+ str(clear_cache))

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

	pickle_query_file = os.path.expanduser(query_file+".pickle")
	if clear_cache == 1:
		if os.path.exists(pickle_query_file):
			print(tistamp(1)+"\tClearing Cache..." + pickle_query_file)
			os.unlink(pickle_query_file)
	if os.path.exists(pickle_query_file):
		print(tistamp(1)+"\tLoading Cache..." + pickle_query_file)
		query_sequences = pickle.load( open( pickle_query_file, "rb" ) )
	else:
		print(tistamp(1)+"\tBuilding Cache..." + pickle_query_file)
		query_sequences = load_sequences(query_file, merge_contigs, seq_line_limit)
		pickle.dump(query_sequences, open( pickle_query_file, "wb" ) )

	print(tistamp(1)+"\tLoaded "+str(len(query_sequences))+" sequences")
	# Access first sequence in dict (should be merged)
	query_sequence = next (iter (query_sequences.values()))

	print(tistamp(1)+"\tGenerating "+str(r_scan_n)+" slices of size "+str(r_scan_size)+"...")

	# Generate random query sequences
	'''
	query_seqs = {}
	if(random_pairwise == 1):
		for i in range(0,r_scan_n):
			(sub_seq,q_pos) = rand_subseq(query_sequence, r_scan_size)
			query_seqs[q_pos] = sub_seq

	print(tistamp(1)+"\tCreated "+str(r_scan_n)+" random sub-sequences.")
	'''
	#elif(sliding_window == 1):
		#query_seqs = sl_window(query_sequence, sl_window_size, sl_increment)


	for database in databases:
		database = ''.join(database[0])
		database_file_name = os.path.expanduser(database)

		dist_cutoff = int(sl_window_size * dist_pcnt)
		print(tistamp(1)+"\tLoading database: "+ database)
		# Load database sequences
		db_sequences = {}
		pickle_db_file = os.path.expanduser(database+".pickle")
		if clear_cache == 1:
			if os.path.exists(pickle_db_file):
				print(tistamp(1)+"\tClearing Cache..." + pickle_db_file)
				os.unlink(pickle_db_file)
		if os.path.exists(pickle_db_file):
			print(tistamp(1)+"\tLoading Cache..." + pickle_db_file)
			db_sequences = pickle.load( open( pickle_db_file, "rb" ) )
		else:
			print(tistamp(1)+"\tBuilding Cache..." + pickle_db_file)
			db_sequences = load_sequences(database, merge_contigs, seq_line_limit)
			pickle.dump(db_sequences, open( pickle_db_file, "wb" ) )

		db_sequence = next (iter (db_sequences.values()))

		search_i = 1
		found_seqs = {}
		while True:
			print(tistamp(1)+"\tSEARCH ROUND: " + str(search_i))

			dup_hits_i = 0
			hits_i = 0
			dup_hits_list = []

			#r_scan_sub_file = os.path.expanduser(output_path+"r_scan_contigs.fas")
			r_scan_sub_file = query_path+"/r_scan_contigs.fas"
			scan_file = query_file # Actual file that will be blasted against the database

			# Access first sequence in dict (should be merged)
			query_sequence = next (iter (query_sequences.values()))

			if search_i == 1:
				print(tistamp(1)+"\tGenerating "+str(r_scan_n)+" slices of size "+str(r_scan_size)+"...")
				r_scan_f = open(os.path.expanduser(r_scan_sub_file), 'w')
				for i in range(0,r_scan_n):
					(sub_seq,seq_start,seq_end) = rand_subseq(query_sequence, r_scan_size)
					r_scan_f.write(">Sub_g:" + str(i) + "_" + str(seq_start) + "_" + str(seq_end) + "\n" + sub_seq + "\n")
				print(tistamp(1)+"\tPrinted "+str(r_scan_n)+" random sub-sequences to "+r_scan_sub_file)
			else:
				print(tistamp(1)+"\tRefining query sequences from: " + str(len(found_seqs)))
				current_expand_size = expand_size*search_i
				current_window_size = current_expand_size*2+r_scan_size
				current_window_increment = int(0.75*current_window_size)
				dist_cutoff = int(current_window_size * dist_pcnt)
				print(tistamp(1)+"\tExpanding: "+ str(current_expand_size))
				print(tistamp(1)+"\tWindow Size: " + str(current_window_size))
				print(tistamp(1)+"\tWindow Increment: " + str(current_window_increment))

				r_scan_f = open(os.path.expanduser(r_scan_sub_file), 'w')
				for q_pos, q_end in found_seqs.items():
					q_left = int(q_pos) - current_expand_size
					q_right = int(q_end) + current_expand_size
					query_seq = query_sequence[q_left:q_right]
					r_scan_f.write(">Sub_g:" + str(i) + "_" + str(seq_start) + "_" + str(seq_end) + "\n" + query_seq + "\n")
					#ldist = distance.levenshtein(q_seq,d_seq)
					#else:
					#ldist = distance.fast_comp(q_seq,d_seq)
				r_scan_f.close()
				# Reset found seqs
				found_seqs = {}

			# Set query file to sub sequence file
			scan_file = r_scan_sub_file
			print(tistamp(1)+"\tTargeting "+database+" with ")
			blst_results_file = query_path+"/qry_"+query_file_name+"_results_"+tistamp(2)+".tsv"
			blst_results_file = os.path.expanduser(blst_results_file)
			print(tistamp(1)+"\tOutputing BLAST results to "+blst_results_file)

			# Form blastn command
			scandb_command = blast_path + seq_type + " -query " + os.path.expanduser(scan_file) \
				+ " -db " + os.path.expanduser(''.join(database)) +\
				" -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\"" +\
				" -out " + blst_results_file
			print(tistamp(1)+"\tRunning: ")
			print(tistamp(1)+"\t"+scandb_command)
			#sts = subprocess.Popen(scandb_command, shell=True).wait()
			os.system(scandb_command)

			print(tistamp(1)+"\t"+str(float(r_scan_size - r_scan_size*dist_pcnt)))
			blst_results = []
			with open(blst_results_file) as inputfile:
			        for line in inputfile:
			                blst_results.append(line.strip().split('\t'))
			print("qseqid\t\t\t\t\tsseqid\t\tpident\t\tlength\t\tmismatch\t\tgapopen\t\tqstart\t\tqend\t\tsstart\t\tsend\t\tevalue\t\tbitscore")
			# 		0				1			2		3			4			5		6			7		8		9		10		11
			# Store blast hits
			blst_hits = {}
			for result in blst_results:
				q_sid = result[0]

				if q_sid in blst_hits:
					blst_hits[q_sid] += 1
				else:
					blst_hits[q_sid] = 1


			# Find hits mapping to only 1 region
			unique_hits = {}
			for q_sid, count in blst_hits.items():
				if count == 1:
					unique_hits[q_sid] = 1

			for unique_q_sid, count in unique_hits.items():
				for result in blst_results:
					q_sid = result[0]
					split_id = q_sid.split('_')
					q_start = split_id[2]
					q_end = split_id[3]
					q_length = result[3]
					for data in result:
						print(data+"\t\t", end="")
					print("")
					if float(q_length) >= float(r_scan_size - r_scan_size*dist_pcnt):
						if q_sid == unique_q_sid:
							print(tistamp(1)+"\tFound at: " + str(q_start) + " => " + str(q_end))
							found_seqs[q_start] = q_end



			# Iterate number of search rounds
			search_i = search_i + 1
			if len(found_seqs) == 0:
				print(tistamp(1)+"\tDone")
				break

for q_start, q_end in found_seqs.items():
	print(q_start + " => " + q_end)


end_time = time.time()
print("Time: " + str(end_time - start_time))
