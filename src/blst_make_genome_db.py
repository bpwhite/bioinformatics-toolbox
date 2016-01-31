# Create a blast database from a genome or assembly
# Created by Bryan White, 2015
# input: list of FASTA files
# output: completed BLAST database
#

import os

db_file = "target_genomes.txt"
blast_path = "~/data_analysis/apps/ncbi-blast-2.2.31+/bin/"
db_type = "nucl"

results = []
with open(db_file) as inputfile:
	for line in inputfile:
		results.append(line.strip().split(','))

for blst_db in results:
	makedb_command = blast_path+"makeblastdb"+" -dbtype "+db_type + " -in " + ''.join(blst_db)
	print(makedb_command)
	os.system(makedb_command)
