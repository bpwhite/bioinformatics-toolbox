# Create a taxonomy-based file tree structure
# Created by Bryan White, 2016
# input: list of taxa to download

# Copyright (c) 2016 Bryan White, bpcwhite@gmail.com

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

## Local libraries
from blst_libs import tistamp

taxa_file = os.path.expanduser("~/data_analysis/data/genome_assemblies/genomes_euks_02062016.tsv")

output_path = "~/data_analysis/data/genome_assemblies/"
seq_downloader = "seq_convert_genbank.pl"

# Load query file names
print(tistamp(1)+"\tReading taxa...")
taxa = []
with open(taxa_file) as inputfile:
        for line in inputfile:
                taxa.append(line.strip().split(','))
for taxon in taxa:
	print(taxon[0])
	exit()
	#split_taxon = taxon.split('\t')
	#print(split_taxon[0])
