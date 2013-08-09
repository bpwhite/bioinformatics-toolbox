Bioinformatics Toolbox (BTBox)
A multi-purpose collection of bioinformatics tools written in Perl, R, and C.

Authors: Bryan P. White (bpcwhite@gmail.com)

Table of Contents
I. Toolbox Overview
	1. What is this toolbox?
	2. Installation guide.
		2a. Windows
		2b. Linux
		2c. Mac
II. DNA Barcoding
	1. dnab_condense.pl - Condenses sequence data down to either unique haplotypes or 
	locales.
	2. dnab_species_delim.pl - For use in DNA barcode species delimitations
III. Genome Annotation
	1. gen_concatenate_alignments.pl - Concatenates alignments to create a combined 
	data set.
	2. gen_concatenate_contigs.pl - Concatenates contigs to form a single contig.
	3. gen_make_blast_db.pl - Creates a blast database.
	4. gen_extract_genes.pl - Extracts target genes from a BLAST database 
IV. Population Genetics
	1. pop_rename_fasta_arlequin.pl - Converts a FASTA file to an Arlequin file.
V. Sequence Manipulation
	1. seq_convert_fasta_nexus.pl - Converts a FASTA file toa  NEXUS file.
	2. seq_order_specimens.pl - Orders FASTA sequences based on an input list.
	3. seq_print_sequences.pl - Prints FASTA sequences with some info about each 
	sequence.
	4. seq_reading_frame.pl - Attempts to put a protein coding sequence into a 
	reading frame that does not have any stop codons.
	5. seq_rename_fasta.pl - Converts the ID's of a genbank FASTA file into a 
	readable format.
	6. seq_rename_fasta_from_list.pl - Renames the ID's of a Fasta file based on a 
	list.
	7. seq_codon_translate.pl - Translate from a protein alignment back to a
	nucleotide alignment. Input is a protein alignment and the original,
	unaligned nucleotide sequences.
	8. seq_convert_genbank.pl - Downloads most genbank info, including nucleotide
	sequence, for a list of target taxa from NCBI.
VI. Simulations
	1. sim_replicators.pl - Particle physics simulator of basic biomolecules

I. Toolbox Overview
	1. What is this toolbox?
	
	If BioPerl is the glue of bioinformatics, then Bioinformatics Toolbox (BTBox) 
	is a series of already glued together programs, with lots of additions on 
	top of BioPerl.
	
	2. Installation Guide
		2a. Windows
		If you are using programs from the BTBox on Windows, do the following:
			1. Install the latest distribution of ActiveState or Strawberry Perl
			2. Run the "installer.pl" script inside this folder.
		2b. Linux
			Coming soon
		2c. Mac
			1. Download and install Xcode
			2. Go to Preferences > Downloads and install Command Line Tools
			3. Update Perl by running "curl -L http://xrl.us/installperlosx | bash" in the Terminal
			4. Check that your Perl version is 5.16 or greater by "Perl -v"
			5. Run Cpan "Module::Name" until all dependencies are installed.
II. DNA Barcoding
	1. dnab_condense.pl
		
		This script takes an input of aligned FASTSA sequences and outputs 
		representative sequences of distinct haplotypes, and the number of 
		those distinct haplotypes at each location. Input sequences must be formatted
		in the following input format from the example.
		
		Example: 
		
		Input: SampleID|Name|Location
		#####################################################
		>10-SCCWRP-4787|Simulium|19354
		ACTTTA.....
		>10-SCCWRP-4788|Simulium|19354
		ACTTTA.....
		>10-SCCWRP-4789|Simulium|19354
		ACTTTA.....
		>10-SCCWRP-4790|Simulium|19354
		ACTTTA.....
		>10-SCCWRP-4791|Simulium|19355
		ACTTTA.....
		>10-SCCWRP-4792|Simulium|19351
		ACTTTA.....
		>10-SCCWRP-4793|Simulium|19349
		ACTTTA.....
		
		Output: HaplotypeID|RepSeq|RepName|Location|Abundance
		#####################################################
		>27|10-SCCWRP-4787|Simulium|19354|4
		ACTTTA.....
		>27|10-SCCWRP-4787|Simulium|19355|1
		ACTTTA.....
		>27|10-SCCWRP-4787|Simulium|19351|1
		ACTTTA.....
		>27|10-SCCWRP-4787|Simulium|19349|1
		ACTTTA.....
		
		This script implements out the following algorithm:
			1. Reduce sequences to distinct haplotype. Distinct haplotypes are
			those haplotyes that differ by nucleotide differences rather than 
			by gaps, spaces, or extra/special/different characters. In other 
			words, these are truly distinct haplotypes.
			2. Iterates through those distinct haplotypes and determines how 
			many locations that haplotype is found at.
			3. Determines how many individuals share that haplotype at each 
			location.
			4. Outputs in FASTSA a haplotype identification number (HaplotypeID), 
			a representative sequence from the haplotype (RepSeq), a 
			representative name from the haplotype (RepName), the location, and 
			the abundance of individuals sharing that haplotype at that location.
		
	2. dnab_species_delimitation.pl
	
	This script inputs a FASTA file of aligned sequences and outputs species-level
	putative delimitations (AKA. Molecular Operational Taxonomic Units MOTUs). The
	input requirements are that the sequence ID's must follow a particular format
	outlined below.
V. 
	8.  seq_convert_genbank.pl
	8a. Parameters
		-query = COI_full, 16S_full, ND2_full, etc.
		-voucher-only = 1, only download sequences with a voucher tag
	8a. Common commands
		1. Download all COI sequences for a taxa, limit the number of sequences 
		per batch to 500, and give the output file a name:
		seq_convert_genbank.pl 	-query COI_full 
								-voucher-only 1 
								-batch-cap 500 
								-term Gastropoda
								-outp Gastropoda_COI
