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
import textwrap

## Local libraries
from blst_libs import tistamp,validate_taxa

#taxa_file = os.path.expanduser("~/data_analysis/data/genome_assemblies/test_genomes_02062016.tsv")
taxa_file = os.path.expanduser("~/data_analysis/data/genome_assemblies/genome_lists/genomes_euks_Fungi_Protists_02082016.tsv")
#taxa_file = os.path.expanduser("~/data_analysis/data/genome_assemblies/genome_lists/test_genomes_02062016.tsv")
#taxa_file = os.path.expanduser("~/data_analysis/data/genome_assemblies/Sacch_with_headers.tsv")


output_path = os.path.expanduser("~/data_analysis/data/genome_assemblies/assemblies_02152016")
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
reset_validation = 1
taxa_i = 0
skip_to = 0
end_at = 999
# Parse data
for taxon in taxa[1:]:
	taxa_i = taxa_i + 1

	# Skip ahead
	if taxa_i <= skip_to:
		continue
	if taxa_i >= end_at:
		break

	split_tx = ''
	#try:
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
	# Status key
	status_d = {}
	status_d['Genomic'] 		= ''
	status_d['Protein'] 		= ''
	status_d['Features'] 		= ''
	status_d['taxon_hierarchy'] 	= ''
	status_d['tax_length_validation'] 	= ''
	status_d['nuc_validiation'] 	= ''
	status_d['prot_validation'] 	= ''
	####################################

	####################################
	# Genome database headers
	gen_organism 		= split_tx[0]
	gen_organism		= gen_organism.replace("\\","_")
	status_d['gen_organism'] 	= gen_organism			= gen_organism.replace("/","_")
	status_d['gen_strain'] 		= gen_strain 			= split_tx[1]
	status_d['gen_biosample'] 	= gen_biosample 		= split_tx[2]
	status_d['gen_bioproject'] 	= gen_bioproject		= split_tx[3]
	status_d['gen_group'] 		= gen_group 			= split_tx[4]
	status_d['gen_subgroup'] 	= gen_subgroup 			= split_tx[5]
	status_d['gen_size_mb'] 	= gen_size_mb 			= split_tx[6]
	status_d['gen_GC_content'] 	= gen_GC_content 		= split_tx[7]
	status_d['gen_assembly'] 	= gen_assembly 			= split_tx[8]
	status_d['gen_replicons'] 	= gen_replicons			= split_tx[9]
	status_d['gen_WGS'] 		= gen_WGS 				= split_tx[10]
	status_d['gen_scaffolds'] 	= gen_scaffolds 		= split_tx[11]
	status_d['gen_genes'] 		= gen_genes 			= split_tx[12]
	status_d['gen_proteins']	= gen_proteins 			= split_tx[13]
	status_d['gen_release_d'] 	= gen_release_d			= split_tx[14]
	status_d['gen_modified_d'] 	= gen_modified_d 		= split_tx[15]
	status_d['gen_level'] 		= gen_level 			= split_tx[16]
	status_d['gen_refseq_FTP']	= gen_refseq_FTP		= split_tx[17]
	status_d['gen_genbank_FTP'] = gen_genbank_FTP 		= split_tx[18]
	status_d['gen_final'] 		= gen_final 			= split_tx[19]
	status_d['gen_assemblyname']= gen_assemblyname 		= split_tx[24]
	status_d['gen_genomic_ftp'] = gen_genomic_ftp 		= split_tx[25]
	status_d['gen_protein_ftp'] = gen_protein_ftp		= split_tx[26]
	status_d['gen_feature_table_ftp'] = gen_feature_table_ftp = split_tx[27]
	####################################



	####################################
	# Build genome key
	gen_key = gen_organism + "_" + gen_strain + "_" + gen_biosample \
		+ "_" + gen_bioproject + "_" + gen_group + "_" + gen_subgroup \
		+ "_" + gen_size_mb
	key_hash = hashlib.sha224(gen_key.encode('utf-8')).hexdigest()
	short_hash = key_hash[0:7]
	status_d['key_hash'] = key_hash
	status_d['short_hash'] = short_hash
	####################################

	####################################
	# Print assembly data
	print(tistamp(1) + "\t[" + str(taxa_i) + "] " + gen_organism + "_" \
		+ gen_strain + "_" + gen_group + "_" + gen_subgroup)
	print(tistamp(1)+"\tBegin Validation Report")
	print(tistamp(1) \
		+ "\t\tBioSample: " + gen_biosample \
		+ "\t\tBioProject: " + gen_bioproject \
		+ "\t\tGroup: " + gen_group \
		+ "\t\tSubGroup: " + gen_subgroup)
	print(tistamp(1) \
		+ "\t\tSize: " + gen_size_mb + " MB" \
		+ "\t\tGC Content: " + gen_GC_content \
		+ "\t\tAssembly: " + gen_assembly)
	####################################

	####################################
	# Build taxonomy lookup command
	seq_limit = 1
	validation = validate_taxa(seq_downloader, gen_organism, seq_limit)


	print(tistamp(1)+"\t\tTaxonomic Hierarchy: " + validation[0])
	print(tistamp(1) \
		+ "\t\tTaxa Length (req): " + str(validation[1]) \
		+ "\t\tSame Hierarchy: " + str(validation[2]) \
		+ "\t\tNuc length: " + str(validation[3]) \
		+ "\t\tProt length: " + str(validation[4]))

	taxon_hierarchy 		= validation[0]
	tax_length_validation 	= validation[1]
	nuc_validiation 		= validation[3]
	prot_validation 		= validation[4]

	# STATUS #
	status_d['taxon_hierarchy'] 	= taxon_hierarchy
	status_d['tax_length_validation'] 	= tax_length_validation
	status_d['nuc_validiation'] 	= nuc_validiation
	status_d['prot_validation'] 	= prot_validation
	# STATUS #

	# Must have both a taxonomic hierarchy string and at least 1 type of sequence (prot or nuc)
	validation_successful = 0
	if tax_length_validation == 'PASS' and (nuc_validiation == 'PASS' or prot_validation == 'PASS'):
		print(tistamp(1)+"\t**Validation Successful**")
		validated_taxa = validated_taxa + 1
		validation_successful = 1
	else:
		print(tistamp(1)+"\t**Validation FAILED**")
	####################################

	####################################
	# Begin downloads
	gen_directory = output_path + "/" + gen_assemblyname
	# Check if assembly directory exists
	print(tistamp(1)+"\t\tChecking directory...")
	if not os.path.isdir(gen_directory):
		print(tistamp(1)+"\t\tDirectory not found, creating directory")
		print(tistamp(1)+"\t\tCreating: " + gen_directory)
		os.makedirs(gen_directory)
	####################################

	####################################
	if validation_successful == 1:
		file_list = ("Genomic", "Protein", "Features")
		file_types = {}
		file_types['Genomic'] = "_genomic.fna.gz"
		file_types['Protein'] = "_protein.fna.gz"
		file_types['Features'] = "_feature_table.txt.gz"

		for file_type in file_list:
			print(tistamp(1)+"\tBegin " + file_type + " Download")
			filename = file_types[file_type]
			# STATUS #
			status_d[file_type] = "Failed"
			# STATUS #

			ftp = ''
			if file_type == "Genomic":
				ftp = gen_genomic_ftp
			if file_type == "Protein":
				ftp = gen_protein_ftp
			if file_type == "Features":
				ftp = gen_feature_table_ftp



			# STATUS #
			status_d['assembly_directory'] = gen_directory
			# STATUS #
			####################################
			download_filename = gen_directory + "/" + gen_assemblyname + filename
			download_string = "wget -nv -O " + download_filename + " " + ftp

			# STATUS #
			status_d[file_type+'_filename'] = download_filename
			# STATUS #
			if reset_validation == 1 and os.path.exists(download_filename):
				subprocess.check_output("rm "+download_filename, shell=True, stderr=subprocess.STDOUT).decode("utf-8")

			####################################
			# Check size, download if 0
			if (os.path.exists(download_filename) and os.stat(download_filename).st_size < 1) \
				or not os.path.exists(download_filename):
				#print(tistamp(1)+"\t\t"+download_string)
				try:
					print(tistamp(1)+"\t\tDownloading "+ file_type + "...")
					download_cmd = subprocess.check_output(download_string, shell=True, \
						stderr=subprocess.STDOUT).decode("utf-8")
					download_text = textwrap.wrap(download_cmd, width=100)
					for dl_line in download_text:
						print(tistamp(1)+"\t\t\t" + dl_line)
					# Update status
					# STATUS #
					status_d[file_type] = "Downloaded"
					# STATUS #
				except:
					subprocess.check_output("rm "+download_filename, shell=True, stderr=subprocess.STDOUT).decode("utf-8")
				#	print(tistamp(1)+"\tCould not download "+file_type)
			####################################

			####################################
			# Unzip file
			if (os.path.exists(download_filename) and os.stat(download_filename).st_size > 1):
				gunzip_file = download_filename.replace(".gz","")
				gunzip_string = "gunzip -c " + download_filename + " > " + gunzip_file
				print(tistamp(1)+"\t\tUnzipping...")
				# print(tistamp(1)+"\t\t"+gunzip_string)
				# Unzip file
				gunzip_cmd = subprocess.check_output(gunzip_string, \
					shell=True, stderr=subprocess.STDOUT).decode("utf-8")
				bytes_to_mb = 1048576
				unzipped_size = os.path.getsize(gunzip_file)/bytes_to_mb
				if file_type == "Genomic":
					pcnt_error = abs(float(gen_size_mb) - float(unzipped_size))/float(gen_size_mb)*100
					print(tistamp(1)+"\t\tPercent Error Genome Size: " + str(round(pcnt_error,2)) + "%")
					if pcnt_error > 5:
						print(tistamp(1)+"\tWARNING: Download Size Too Small!")
						print(tistamp(1)+"\t" + str(gen_size_mb) + " < " + str(unzipped_size))
				check_unzip_cmd = subprocess.check_output("head -n 1 "+gunzip_file, \
					shell=True, stderr=subprocess.STDOUT).decode("utf-8").rstrip()
				print(tistamp(1)+"\t"+check_unzip_cmd)

				# STATUS #
				status_d[file_type+'_first_line'] = check_unzip_cmd
				status_d[file_type+'_unzipped'] = "Unzipped"
				status_d[file_type+'_unzipped_filename'] = gunzip_file
				status_d[file_type+'_unzipped_size'] = unzipped_size
				status_d[file_type+'_pcnt_error_genome'] = pcnt_error
				# STATUS #
			print(tistamp(1)+"\t\t"+file_type+" Complete")
			####################################
	####################################
	# Print Status
	print(tistamp(1)+"\tPrinting status...")
	status_filename = gen_directory + "/" + gen_assemblyname \
		+ "_status.txt"
	#print(tistamp(1)+status_filename)
	status_f = open(os.path.expanduser(status_filename), 'w')
	for key, value in status_d.items():
		print(tistamp(1)+"\t\t" + key + "=" + str(value))
		status_f.write(key + "=" + str(value) + "\n")
	####################################

	print(tistamp(1)+"\t***********************")
	print(tistamp(1))
	#exit()

print(tistamp(1)+"\tSuccessfully validated: " + str(validated_taxa) + " out of "+ str(taxa_i))

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
