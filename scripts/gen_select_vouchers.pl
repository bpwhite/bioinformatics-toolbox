#!/usr/bin/perl
# Select vouchered specimens from a sql database
#
# Copyright (c) 2013, Bryan White, bpcwhite@gmail.com

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
use lib "$FindBin::Bin/libs/";
use General::Arguments;
use Sequence::PairwiseAlign;
use Sequence::Kimura_Distance;
use DBI;

use strict;
use warnings;

my $params = General::Arguments->new(	arguments_v => \@ARGV,
										option_defs => {'-database' 	=> '', 				# Database name
														'-location'		=> 'localhost',		# IP address/hostname of database
														'-port'			=> '3306',			# DB port
														'-user'			=> '',				# Username
														'-password'		=> '',				# User password	
														'-table'		=> '',				# Table to input data to
														'-voucher-list'	=> '',				# CSV file to upload
														'-outp'			=> 'output',		# Output prefix
														'-max-length'	=> '3000',			# Maximum sequence length allowed
														'-bsreps'		=> '1',				# Number of bootstrap replicates to output
														'-bsrep-size'	=> '1000',			# Default bs rep size
													}
													);

my $database 		= $params->options->{'-database'};
my $location 		= $params->options->{'-location'};
my $port			= $params->options->{'-port'};
my $user 			= $params->options->{'-user'};
my $password 		= $params->options->{'-password'};
my $table 			= $params->options->{'-table'};
my $voucher_file	= $params->options->{'-voucher-list'};
my $output_prefix	= $params->options->{'-outp'};
my $max_seq_length	= $params->options->{'-max-length'};
my $bs_reps			= $params->options->{'-bsreps'};
my $bs_rep_size		= $params->options->{'-bsrep-size'};

my $distance_cutoff = 0.90;
my $comparison_cutoff = 350;
my $validate_alignment = 0;
my $minimum_seq_length = 500;

print "Opening $voucher_file \n";
open VOUCHER, "< $voucher_file";
my @voucher_list = <VOUCHER>;
close VOUCHER or die "Couldn't close $voucher_file!\n";

my $dsn = "DBI:mysql:database=$database;host=$location;port=$port";
my $dbh = DBI->connect($dsn, $user, $password);

# my @gene_list = ("COI", "16S", "18S", "28S", "CYTB");
my @required_gene_list 	= ("COI", "16S");
my @gene_list 			= ("COI", "16S");

my %reference_hash = (
						'COI' => 'TCGCGACAATGATTATTTTCTACAAATCATAAAGATATTGGAACTTTATATTTTATTTTTGGAGCTTGAGCTGGAATAGTTGGAACATCTTTAAGAATTTTAATTCGAGCTGAATTAGGACATCCTGGAGCATTAATTGGAGATGATCAAATTTATAATGTAATTGTAACTGCACATGCTTTTATTATAATTTTTTTTATGGTTATACCTATTATAATTGGTGGATTTGGAAATTGATTAGTGCCTTTAATATTAGGTGCTCCTGATATAGCATTCCCACGAATAAATAATATAAGATTTTGACTACTACCTCCTGCTCTTTCTTTACTATTAGTAAGTAGAATAGTTGAAAATGGAGCTGGAACAGGATGAACTGTTTATCCACCTTTATCCGCTGGAATTGCTCATGGTGGAGCTTCAGTTGATTTAGCTATTTTTTCTCTACATTTAGCAGGGATTTCTTCAATTTTAGGAGCTGTAAATTTTATTACAACTGTAATTAATATACGATCAACAGGAATTTCATTAGATCGTATACCTTTATTTGTTTGATCAGTAGTTATTACTGCTTTATTATTATTATTATCACTTCCAGTACTAGCAGGAGCTATTACTATATTATTAACAGATCGAAATTTAAATACATCATTTTTTGACCCAGCGGGAGGAGGAGATCCTATTTTATATCAACATTTATTTTGATTTTTTGGTCACCCTGAAGTTTATATTTTAATTTTACCTGGATTTGGAATAATTTCTCATATTATTAGACAAGAATCAGGAAAAAAGGAAACTTTTGGTTCTCTAGGAATAATTTATGCTATATTAGCTATTGGATTATTAGGATTTATTGTATGAGCTCATCATATATTTACCGTTGGAATAGATGTAGATACTCGAGCTTATTTTACCTCAGCTACTATAATTATTGCAGTTCCTACTGGAATTAAAATTTTTAGTTGATTAGCTACTTTACATGGAACTCAACTTTCTTATTCTCCAGCTATTTTATGAGCTTTAGGATTTGTTTTTTTATTTACAGTAGGAGGATTAACAGGAGTTGTTTTAGCTAATTCATCAGTAGATATTATTTTACATGATACTTATTATGTAGTAGCTCATTTTCATTATGTTTTATCTATAGGAGCTGTATTTGCTATTATAGCAGGTTTTATTCACTGATACCCCTTATTTACTGGATTAACGTTAAATAATAAATGATTAAAAAGTCATTTCATTATTATATTTATTGGAGTTAATTTAACATTTTTTCCTCAACATTTTTTAGGATTGGCTGGAATACCTCGACGTTATTCAGATTACCCAGATGCTTACACAACATGAAATATTGTATCAACTATTGGATCAACTATTTCATTATTAGGAATCTTATTCTTTTTTTTTATTATTTGAGAAAGTTTAGTATCACAACGACAAGTAATTTACCCAATTCAACTAAATTCATCAATTGAATGATACCAAAATACTCCGCCAGCTGAACATAGATATTCTGAATTACCACTTTTAACAAATTAA',
						'16S' => 'GCTAAACCTAGCCCCAAACCCACTCCACCTTACTACCAGACAACCTTAGCCAAACCATTTACCCAAATAAAGTATAGGCGATAGAAATTGAAACCTGGCGCAATAGATATAGTACCGCAAGGGAAAGATGAAAAATTATAACCAAGCATAATATAGCAAGGACTAACCCCTATACCTTCTGCATAATGAATTAACTAGAAATAACTTTGCAAGGAGAGCCAAAGCTAAGACCCCCGAAACCAGACGAGCTACCTAAGAACAGCTAAAAGAGCACACCCGTCTATGTAGCAAAATAGTGGGAAGATTTATAGGTAGAGGCGACAAACCTACCGAGCCTGGTGATAGCTGGTTGTCCAAGATAGAATCTTAGTTCAACTTTAAATTTGCCCACAGAACCCTCTAAATCCCCTTGTAAATTTAACTGTTAGTCCAAAGAGGAACAGCTCTTTGGACACTAGGAAAAAACCTTGTAGAGAGAGTAAAAAATTTAACACCCATAGTAGGCCTAAAAGCAGCCACCAATTAAGAAAGCGTTCAAGCTCAACACCCACTACCTAAAAAATCCCAAACATATAACTGAACTCCTCACACCCAATTGGACCAATCTATCACCCTATAGAAGAACTAATGTTAGTATAAGTAACATGAAAACATTCTCCTCCGCATAAGCCTGCGTCAGATTAAAACACTGAACTGACAATTAACAGCCCAATATCTACAATCAACCAACAAGTCATTATTACCCTCACTGTCAACCCAACACAGGCATGCTCATAAGGAAAGGTTAAAAAAAGTAAAAGGAACTCGGCAAATCTTACCCCGCCTGTTTACCAAAAACATCACCTCTAGCATCACCAGTATTAGAGGCACCGCCTGCCCAGTGACACATGTTTAACGGCCGCGGTACCCTAACCGTGCAAAGGTAGCATAATCACTTGTTCCTTAAATAGGGACCTGTATGAATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTACCTAACAAACCCACAGGTCCTAAACTACCAAACCTGCATTAAAAATTTCGGTTGGGGCGACCTCGGAGCAGAACCCAACCTCCGAGCAGTACATGCTAAGACTTCACCAGTCAAAGCGAACTACTATACTCAATTGATCCAATAACTTGACCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCTAGAGTCCATATCAACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGATTAAAGTCCTACGTGATCTGAGTTCAGACCGGAGTAATCCAGGTCGGTTTCTATCTACNTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTACTTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTATTATACCCACACCCACCCAAGAACAGGGTTT',
					);
					
my %aligned_seq_hash = ();
my %query_hash = ();
my %used_voucher_hash = ();

foreach my $query_gene (@gene_list) {
	for (my $i = 0; $i < $bs_reps; $i++) {
		my $current_file = $output_prefix.'_'.$query_gene.'_'.$i.'.fas';
		unlink ($current_file);
	}
}

my $num_vouchers = scalar @voucher_list;
my $line_counter = 0;
if ($bs_reps == 1) {
	$bs_rep_size = $num_vouchers;
}
for (my $bs_reps_i = 0; $bs_reps_i < $bs_reps; $bs_reps_i++) {
	print "Creating rep: ".$bs_reps_i."\n";
	my $bs_sample_created = 0;
	my $voucher_counter = 0;
	while ($bs_sample_created < $bs_rep_size) {
	# foreach my $voucher_line (@voucher_list) {
		my $current_voucher = '';
		if($bs_reps == 1) {
			$current_voucher = $voucher_list[$voucher_counter];
		} elsif ($bs_reps > 1) {
			$current_voucher = $voucher_list[int rand ($num_vouchers)];
			if (exists $used_voucher_hash{$current_voucher}) {
				$current_voucher = $current_voucher."_".($used_voucher_hash{$current_voucher}+1);
			} else {
				$used_voucher_hash{$current_voucher} = 1;
			}
		}
		# next;
		# print "[".$line_counter."]\n";
		$line_counter++;
		# print $line."\n";
		$current_voucher =~ s/\n//g; # Strip newlines
		# print $line."\n";
		my ($query_primary_key, $taxon_query, $taxon_id, $accession,
			$seq_query, $seqs_found, $description, $gene,
			$product, $binom, $tax_hierarchy, $pub_title, $pub_authors,
			$pub_abstract_text, $jrn_name, $jrn_doi, $jrn_so,
			$jrn_volume, $jrn_issue, $jrn_pages, $jrn_pubdate, $nuc_seq,
			$prot_seq, $primers, $codon_start, $collection_date, $voucher_id, 
			$collected_by, $identified_by, $organelle, $location, $lat_lon);
		
		my $incomplete_gene_set = 0;
		foreach my $query_gene (@required_gene_list) {
			# my $query = "SELECT * FROM column_names WHERE seq_query LIKE \"\%".$query_gene."\%\" AND voucher_id = \"".$voucher_line."\"";
			my $query = "SELECT nuc_seq FROM column_names WHERE seq_query LIKE \"\%".$query_gene."\%\" AND voucher_id = \"".$current_voucher."\"";
			
			my $statement = $dbh->prepare($query);
			$statement->execute() or die "$@\n";
			if(($statement->rows == 0)) {
				$incomplete_gene_set = 1;
				last;
			}
			# $statement->bind_col(22, \$nuc_seq);
			$statement->bind_col(1, \$nuc_seq);
			$query_hash{$query} = $nuc_seq;
			
			while($statement->fetch()) {
				if((length($nuc_seq) > $max_seq_length) || (length($nuc_seq) < $minimum_seq_length)) {
					$incomplete_gene_set = 1;
					last;
				}
				last; # only check the first fetched sequence for that voucher ID
			}
			
			if($validate_alignment == 1) {
				my ($refseq_aligned, $query_seq_aln) = alignseqs($reference_hash{$query_gene}, $nuc_seq);
				my ($K2P, $transitions, $transversions, $num_comparisons);
				($K2P, $transitions, $transversions, $num_comparisons) = k2p_no_bs(\$refseq_aligned, \$query_seq_aln, 
																					length($query_seq_aln));
				if(($K2P > $distance_cutoff) || ($num_comparisons < $comparison_cutoff)) {
					print $K2P."\n";
					$incomplete_gene_set = 1;
					last;
				}
				$aligned_seq_hash{$nuc_seq} = $query_seq_aln;
			}
		}
		if($incomplete_gene_set == 1) {
			$voucher_counter++;
			if ($bs_reps == 1) {
				$bs_sample_created++;
			}
			next;
		}
		
		# Has a complete gene set.
		foreach my $query_gene (@gene_list) {
			# my $query = "SELECT * FROM column_names WHERE seq_query LIKE \"\%".$query_gene."\%\" AND voucher_id = \"".$voucher_line."\"";
			my $query = "SELECT nuc_seq FROM column_names WHERE seq_query LIKE \"\%".$query_gene."\%\" AND voucher_id = \"".$current_voucher."\"";
		
			my $statement = $dbh->prepare($query);
			$statement->execute() or die "$@\n";
			
			# BIND TABLE COLUMNS TO VARIABLES
			# $statement->bind_col(1, \$query_primary_key);
			# $statement->bind_col(2, \$taxon_query);
			# $statement->bind_col(3, \$taxon_id);
			# $statement->bind_col(4, \$accession);
			# $statement->bind_col(5, \$seq_query);
			# $statement->bind_col(6, \$seqs_found);	
			# $statement->bind_col(7, \$description);	
			# $statement->bind_col(8, \$gene);
			# $statement->bind_col(9, \$product);
			# $statement->bind_col(10, \$binom);
			# $statement->bind_col(11, \$tax_hierarchy);
			# $statement->bind_col(12, \$pub_title);
			# $statement->bind_col(13, \$pub_authors);	
			# $statement->bind_col(14, \$pub_abstract_text);
			# $statement->bind_col(15, \$jrn_name);
			# $statement->bind_col(16, \$jrn_doi);
			# $statement->bind_col(17, \$jrn_so);
			# $statement->bind_col(18, \$jrn_volume);	
			# $statement->bind_col(19, \$jrn_issue);
			# $statement->bind_col(20, \$jrn_pages);
			# $statement->bind_col(21, \$jrn_pubdate);	
			
			# $statement->bind_col(22, \$nuc_seq);
			$statement->bind_col(1, \$nuc_seq);

			# $statement->bind_col(23, \$prot_seq);
			# $statement->bind_col(24, \$primers);
			# $statement->bind_col(25, \$codon_start);
			# $statement->bind_col(26, \$collection_date);
			# $statement->bind_col(27, \$voucher_id);
			# $statement->bind_col(28, \$collected_by);
			# $statement->bind_col(29, \$identified_by);
			# $statement->bind_col(30, \$organelle);
			# $statement->bind_col(31, \$location);
			# $statement->bind_col(32, \$lat_lon);
		
			# Open file for append
			my $aligned_current_seq = $nuc_seq;
			if($validate_alignment == 1) {
				$aligned_seq_hash{$nuc_seq};
			}
			
			my $current_file = $output_prefix.'_'.$query_gene.'_'.$bs_reps_i.'.fas';
			open (FASTAFILE, '>>'.$current_file)  or die "Could not append to $current_file\n";
			
			if($statement->rows == 0) {
				print FASTAFILE ">$query_gene|$current_voucher\n---\n";
				close FASTAFILE;
				next;
			}
			while($statement->fetch()) {
			   print FASTAFILE ">$query_gene|$current_voucher\n$aligned_current_seq\n";
			   last;
			}
			close FASTAFILE;
		}
		# exit if $line_counter == 1000;
		$voucher_counter++;
		$bs_sample_created++;
	}

}
