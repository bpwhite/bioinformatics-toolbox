# Extracts genes from a blast database.
#
# Copyright (c) 2011, Bryan White

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

use FindBin;
use lib $FindBin::Bin;
use FastaTools;

use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;


print "Blast database file: ";
# chomp(my $input_file = <>);
my $input_file = "bacterial.fasta";

# Fix the fasta
# If first run
print "Fix fasta (remove blanks - do if first run) [y]es or [n]\n";
# chomp(my $fix_fasta = <>);
my $fix_fasta = "no";

if($fix_fasta eq "yes" || $fix_fasta eq "y") {
	FastaTools::fix_fasta($input_file);
}
print "Query file list: \n";
# chomp(my $query_file_list = <>);
my $query_file_list = "query_files.txt";

open QUERY, "< $query_file_list" or die "Can't open $query_file_list : $!";
my @query_files = <QUERY>;
close QUERY;

my @blast_result_files = (); # Store the blast results in these output files.

my $blast_flag = 0;
# Cycle through the files
# my @output_gene_file = ();
foreach my $query_file (@query_files) {
	my @parsed_output_file = split("#",$query_file);
	next if(scalar(@parsed_output_file) > 1);
	
	@parsed_output_file = split(".txt",$query_file);
	my $output_gene_file = $parsed_output_file[0]."_blast.csv";
	print "BLASTing :".$parsed_output_file[0]."\n";
	if($blast_flag == 1) {
		system("blastn 
			-query ".$query_file
			." -db ".$input_file 
			." -word_size 6
			-outfmt 10
			-min_raw_gapped_score 100
			-num_threads 2 
			-gapopen 0 
			-gapextend 2 
			-penalty -1 
			-reward 1 
			-num_alignments 500
			-out ".$output_gene_file."");
	}
			

		push(@blast_result_files, $output_gene_file);
}

# Cycle through the blast output files and extract the genes.
foreach my $blast_file (@blast_result_files) {
	extract_genes($blast_file,$input_file);
}

## End program
## Begin subs

sub extract_genes {
	# Load the gene output file
	my $gene_file = shift;
	my $input_file = shift;
	my @parsed_output_file = split("_",$gene_file);
	my $output_file = '';
	my $output_name = $parsed_output_file[0];
	if(scalar(@parsed_output_file) == 1) {
		$output_file = $parsed_output_file[0]."_output.fasta";
		$output_name = $parsed_output_file[0];
	} else {
		$output_file = $parsed_output_file[0]."_".$parsed_output_file[1]."_output.fasta";
		$output_name = $parsed_output_file[1];
	}
	print "Extracting :".$parsed_output_file[1]."\n";
	# Parse the gene file
	open F, "< $gene_file" or die "Can't open $gene_file : $!";
	my @gene_output = <F>;
	close F;

	my $gene_counter_i = 1;
	my $gene_length = 0; # This will be the maximum gene length.

	
	unlink $output_file;
	open (MYFILE, '>>'.$output_file);
	my $reversals_file = $output_file."_reversals.txt";
	open (REVERSALS, '>>'.$reversals_file);
	unlink $reversals_file;
	my @reversals = ();
	
	foreach my $line (@gene_output) {
		
		my @output_lines = split(",",$line);

		# Blast result variables
		my $seq_id			= $output_lines[1];
		my $alignment_length= $output_lines[3];
		my $q_start			= $output_lines[6];
		my $q_end			= $output_lines[7];
		my $gene_start		= $output_lines[8];
		my $gene_end		= $output_lines[9];
		my $gene_name	 	= $output_lines[0];
		# Set the gene length to the length of the query sequence.
		if($gene_counter_i == 1){
			$gene_length = $alignment_length;
		}
		
		# If the gene is backwards, reverse the read order.
		if($gene_start > $gene_end) {
			print "***".$seq_id." for ".$gene_name." reversed***.\n";
			print REVERSALS "***".$seq_id." for ".$gene_name." reversed***.\n";
			my $gene_start_orig 	= $gene_start;
			my $gene_end_orig		= $gene_end;
			$gene_end = $gene_start_orig;
			$gene_start = $gene_end_orig;
			# Adjust gene start and end for missing data
			# Find the missign start data
			# print "start before ".$gene_start."\n";
			my $missing_end = $gene_length - $q_end;
			$gene_start = $gene_start - $missing_end;
			# print "start ".$gene_start."\n";
			# Find the missing end data
			
			# print "end before ".$gene_end."\n";
			$gene_end 	= $gene_end+$q_start;
			# print "end ".$gene_end."\n";
			my $final_length = $gene_end - $gene_start;
			print "Length ".$final_length."\n";
		} else {
			# Adjust gene start and end for missing data
			# Find the missign start data
			$gene_start = $gene_start - $q_start;
			# Find the missing end data
			my $missing_end = $gene_length - $q_end;
			$gene_end 	= $gene_end+$missing_end;
		}
		
		# Load the genome file
		my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');

		while((my $seqobj = $seqio->next_seq())) {
			
			# Split the ID from the blast output to get the accession #
			my @blast_id = split(/\|/,$seq_id);
			
			# Same thing with the genomic sequences
			my @genome_id = split(/\|/,$seqobj->display_id);
			
			my @gene_id = ();
			my $output_id = '';
			if(scalar(@genome_id) == 1 && scalar(@blast_id) == 1) {
				
				if($blast_id[0] eq $genome_id[0]) {
					# Final gene length modifications
					# Make sure the output length does not exceed the edges of the contig.
					if($gene_end > $seqobj->length) {
						$gene_end = $seqobj->length;
					}
					if($gene_start <= 0) {
						$gene_start = 1;
					}
					$output_id = $gene_name."|".$genome_id[0];
					print MYFILE ">".$output_id."\n";
					print MYFILE $seqobj->subseq($gene_start,$gene_end)."\n";
					last;
				}
			} else {

				if($blast_id[3] eq $genome_id[3]) {
					# Split the genome name further
					@gene_id = split(/\_/,$genome_id[4]);
				
					# Concatenate the final output.
					if($gene_id[3] eq "mitochondrion,") {
						$gene_id[3] = "";
					} else {
						# append sub-species names when available
						$gene_id[3] = "_".$gene_id[3];
					}
					$output_id = $gene_name."|".$genome_id[3]."|".$gene_id[1]."_".$gene_id[2].$gene_id[3];
					
					# Final gene length modifications
					# Make sure the output length does not exceed the edges of the contig.
					if($gene_end > $seqobj->length) {
						$gene_end = $seqobj->length;
					}
					if($gene_start <= 0) {
						$gene_start = 1;
					}
					
					print MYFILE ">".$output_id."\n";
					print MYFILE $seqobj->subseq($gene_start,$gene_end)."\n";
					last;
				}
			}
		}
		
		$gene_counter_i++;
	}
	print "Output ".($gene_counter_i-1)." ".$parsed_output_file[1]." sequences.\n";
	close(MYFILE);
	close(REVERSALS);
}

# Subs