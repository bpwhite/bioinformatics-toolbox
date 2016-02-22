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
blast_scan = 0
mem_scan = 1
r_scan_size = 11
sl_window_size = 11
r_scan_n = 500
merge_contigs = 1
seq_line_limit = 9999999999
sl_increment = int(0.75*sl_window_size)
dist_cutoff = 2
random_pairwise = 0
sliding_window = 1
expand_size = 10
#pickle_cache = 'build'
clear_cache = 0

#pickle_cache = 1

print(random_scan)
print(blast_scan)
print(mem_scan)
print(r_scan_size)
print(r_scan_n)
print(merge_contigs)
print(seq_line_limit)
print(sl_window_size)
print(sl_increment)
print(dist_cutoff)
print("Cache: "+ str(clear_cache))
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
	query_seqs = {}
	#if(random_pairwise == 1):
	for i in range(0,r_scan_n):
		(sub_seq,q_pos) = rand_subseq(query_sequence, r_scan_size)
		query_seqs[q_pos] = sub_seq

	print(tistamp(1)+"\tCreated "+str(r_scan_n)+" random sub-sequences.")
	#elif(sliding_window == 1):
		#query_seqs = sl_window(query_sequence, sl_window_size, sl_increment)


	for database in databases:
		database = ''.join(database[0])

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

		db_window_seqs = {}
		db_window_seqs = sl_window(db_sequence, sl_window_size, sl_increment)
		print(tistamp(1)+"\tCreated "+str(len(db_window_seqs))+" window sequences.")

		search_i = 1
		found_seqs = {}
		while True:
			print(tistamp(1)+"\tSEARCH ROUND: " + str(search_i))

			dup_hits_i = 0
			hits_i = 0
			dup_hits_list = []

			if search_i > 1:
				print(tistamp(1)+"\tRefining query sequences from: " + str(len(found_seqs)))
				current_expand_size = expand_size*search_i
				current_window_size = current_expand_size*2
				current_window_increment = int(0.75*current_window_size)
				current_dist_cutoff = int(current_window_size*0.25)
				print(tistamp(1)+"\tExpanding: "+ str(current_expand_size))
				print(tistamp(1)+"\tWindow Size: " + str(current_window_size))
				print(tistamp(1)+"\tWindow Increment: " + str(current_window_increment))

				query_seqs = {}
				for q_pos, d_pos in found_seqs.items():
					q_left = q_pos - current_expand_size
					q_right = q_pos + current_expand_size
					query_seqs[q_left] = query_sequence[q_left:q_right]
					#print(query_sequence[q_left:q_right])
					#exit()
				# Validate window Size
				if current_window_size != current_expand_size:
					print(tistamp(1)+"\tQuery sequence size does not match window size.")
				# Clear previously found sequences
				found_seqs = {}
				# Generate new window seqs based on size expand
				db_window_seqs = sl_window(db_sequence, current_window_size, current_window_increment)
				print(tistamp(1)+"\tCreated "+str(len(db_window_seqs))+" window sequences.")

			print(tistamp(1)+"\tScanning Windows...")

			# Query dict
			min_dist = 99999999
			q_scan_i = 0
			q_scan_max = len(query_seqs)
			q_scan_inc = int(q_scan_max * 0.01)
			#print(tistamp(1))
			for q_pos, q_seq in query_seqs.items():
				found_hit = 0
				'''
				if (q_scan_i % q_scan_inc) == 0:
					pcnt_complete = int(q_scan_i/q_scan_max)*100
					print("=", end="")
				'''
				# Database dict
				for d_pos, d_seq in db_window_seqs.items():
					ldist = ''

					if search_i > 1:
						ldist = distance.levenshtein(q_seq,d_seq)
					else:
						ldist = distance.fast_comp(q_seq,d_seq)

					if ldist < min_dist:
						min_dist = ldist

					if ldist != -1 and ldist < dist_cutoff:
						hits_i = hits_i + 1
						found_hit = found_hit + 1
						found_seqs[q_pos] = d_pos
						# Skip seqs matching more than 1 place in the target genome
						if found_hit > 1:
							dup_hits_i = dup_hits_i + 1
							dup_hits_list.append(q_pos)
							continue
				q_scan_i = q_scan_i + 1
			print(tistamp(1)+"\tMinimum distance encountered: " + str(min_dist))
			for dup_q_pos in dup_hits_list:
				if dup_q_pos in found_seqs:
					del found_seqs[dup_q_pos]

			print("Found " + str(hits_i)+ " hits")
			print("Saved " + str(len(found_seqs)) + " hits")
			print("Duplicate hits: " + str(dup_hits_i))
			if search_i > 1:
				for q_pos, d_pos in found_seqs.items():
						print(str(q_pos) + " \t " \
						+ str(d_pos) + " \t " \
						+ str(query_seqs[q_pos]) + " \t " \
						+ db_window_seqs[d_pos])

			# Iterate number of search rounds
			search_i = search_i + 1
			if len(found_seqs) == 0:
				break




end_time = time.time()
print("Time: " + str(end_time - start_time))
