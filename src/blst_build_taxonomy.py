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
import hashlib

## Local libraries
from blst_libs import tistamp,validate_taxa

taxa_file = os.path.expanduser("~/data_analysis/data/genome_assemblies/test_genomes_02062016.tsv")

output_path = "~/data_analysis/data/genome_assemblies/"
seq_downloader = os.path.expanduser("~/data_analysis/code/bioinformatics-toolbox/src/seq_convert_genbank.pl")

# Parameters
parse_headers = 0
taxa_col = 0

# Load taxa file data
print(tistamp(1)+"\tReading taxa...")
taxa = []
with open(taxa_file) as inputfile:
        for line in inputfile:
                taxa.append(line.strip().split(','))

validated_taxa = 0
taxa_i = 0
# Parse data
for taxon in taxa[1:]:
	#print(taxon[0])

	split_tx = taxon[0].split('\t')

	####################################
	# Print header codes
	if(parse_headers == 1):
		headers_i = 0
		for tax_col in split_tx:
			print(tax_col + " \t\t= split_tx["+str(headers_i)+"]")
			headers_i = headers_i + 1
		exit()
	####################################

	####################################
	# Genome database headers
	gen_organism 		= split_tx[0]
	gen_strain 			= split_tx[1]
	gen_biosample 		= split_tx[2]
	gen_bioproject 		= split_tx[3]
	gen_group 			= split_tx[4]
	gen_subgroup 		= split_tx[5]
	gen_size_mb 		= split_tx[6]
	gen_GC_content 		= split_tx[7]
	gen_assembly 		= split_tx[8]
	gen_replicons 		= split_tx[9]
	gen_WGS 			= split_tx[10]
	gen_scaffolds 		= split_tx[11]
	gen_genes 			= split_tx[12]
	gen_proteins 		= split_tx[13]
	gen_release_d		= split_tx[14]
	gen_modified_d 		= split_tx[15]
	gen_level 			= split_tx[16]
	gen_refseq_FTP 		= split_tx[17]
	gen_genbank_FTP 	= split_tx[18]
	gen_final 			= split_tx[19]
	gen_assemblyname 	= split_tx[24]
	gen_genomic 		= split_tx[25]
	gen_protein 		= split_tx[26]
	gen_feature_table 	= split_tx[27]
	####################################

	####################################
	# Build genome key
	gen_key = gen_organism + "_" + gen_strain + "_" + gen_biosample \
		+ "_" + gen_bioproject + "_" + gen_group + "_" + gen_subgroup \
		+ "_" + gen_size_mb
	key_hash = hashlib.sha224(gen_key.encode('utf-8')).hexdigest()
	short_hash = key_hash[0:7]
	print(tistamp(1) + "\t [" + str(taxa_i) + "] " + gen_organism + "_" + short_hash)
	####################################

	####################################
	# Build taxonomy lookup command
	seq_limit = 1
	validation = validate_taxa(seq_downloader, gen_organism, seq_limit)
	print(tistamp(1)+"\tValidation Report For: " + gen_organism)
	print(tistamp(1)+"\t\tTaxonomic Hierarchy: "+ validation[0])
	print(tistamp(1)+"\t\tTaxa Length (req): " + str(validation[1]))
	print(tistamp(1)+"\t\tSame Hierarchy: " + str(validation[2]))
	print(tistamp(1)+"\t\tNuc length: " + str(validation[3]))
	print(tistamp(1)+"\t\tProt length: " + str(validation[4]))

	taxon_hierarchy 		= validation[0]
	tax_length_validation 	= validation[1]
	nuc_validiation 		= validation[3]
	prot_validation 		= validation[4]

	# Must have both a taxonomic hierarchy string and at least 1 type of sequence (prot or nuc)
	if tax_length_validation == 'PASS' and (nuc_validiation == 'PASS' or prot_validation == 'PASS'):
		print(tistamp(1)+"\t**Validation Successful**\n")
		validated_taxa = validated_taxa + 1
	else:
		print(tistamp(1)+"\t**Validation Failed**\n")

	taxa_i = taxa_i + 1
	#exit()

print(tistamp(1)+"\t Successfully validated: " + validated_taxa)

'''
Genome columns as of:
2016-02-07 17:28:56	Reading taxa...
[0] #Organism/Name
[1] Strain
[2] BioSample
[3] BioProject
[4] Group
[5] SubGroup
[6] Size (Mb)
[7] GC%
[8] Assembly
[9] Replicons
[10] WGS
[11] Scaffolds
[12] Genes
[13] Proteins
[14] Release Date
[15] Modify Date
[16] Level
[17] RefSeq FTP
[18] GenBank FTP
[19] Final
[20] 1
[21] 2
[22] 3
[23] 4
[24] AssemblyName
[25] Genomic
[26] Protein
[27] Feature_Table
'''
