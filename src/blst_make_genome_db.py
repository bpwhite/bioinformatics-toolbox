# Create a blast database from a genome or assembly
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
