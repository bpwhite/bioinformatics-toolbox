#!/usr/bin/env perl
# Resamples an alignment or splits it into even parts
#
# Copyright (c) 2011, Bryan White, bpcwhite@gmail.com

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
use lib "$FindBin::Bin/libs/Sequence"; 
use lib "$FindBin::Bin/libs/";

use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::SeqIO::fasta;
use Statistics::Descriptive;
use Math::Random::MT qw(srand rand);

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

print "Enter input filename:";
chomp(my $fasta_input = <>);

print "Resample [1] or Split [2]?";
chomp(my $splitorsample = <>);
# my $splitorsample = 1;

# Parameters
my $sample_size = 0;
my $replicates = 0;
my $chunk_size = 0;
if ($splitorsample == 1) {
	print "Sample size?: ";
	chomp($sample_size = <>);
	print "Replicates?: ";
	chomp($replicates = <>);
	# $replicates = 10;
} elsif($splitorsample == 2) {
	print "Split into what size chunks?";
	chomp($chunk_size = <>);
}

my $seqio  = Bio::SeqIO->new(-file => $fasta_input, '-format' => 'Fasta');

my @output_split = split(/\./,$fasta_input);
my $output_name = $output_split[0];

my %seq_hash = ();
my %unique_seqs = ();

my $num_seqs = 0;
my $num_unique_seqs = 0;



my $seq;
my @seq_array;
while( $seq = $seqio->next_seq() ) {
	# print $seq->display_id."\n";
    push(@seq_array,$seq);
}
my $num_sequences = scalar(@seq_array);
# print $num_sequences."\n";

# Random resampling
if($splitorsample == 1) {

	for(my $rep_number = 1; $rep_number <= $replicates; $rep_number++) {
		# Generate random numbers
		my @rand_numbers = ();
		for(my $i = 0; $i < $sample_size; $i++) {
			my $random_number = int(rand($num_sequences));
			push(@rand_numbers, $random_number);
		}
		

		
		my $new_file = '';
		$new_file = $output_name."_resampled_n_".$sample_size."_rep_".$rep_number."_aln.fas";
		unlink $new_file;
		open (MYFILE, '>>'.$new_file);
		
		foreach my $random_number (@rand_numbers) {
			print MYFILE ">".$seq_array[$random_number]->display_id."\n";
			print MYFILE $seq_array[$random_number]->subseq(1,$seq_array[$random_number]->length)."\n";
		}
		close MYFILE;
	}
	
	print "Done!";
}

# Split an alignment into chunks
if($splitorsample == 2) {
	
	my $number_files = roundup($num_sequences/$chunk_size);
	print $number_files."\n";
	
	my $seq_counter = 0;
	for(my $file_num = 1; $file_num <= $number_files; $file_num++) {

		print "File number: ".$file_num."\n";
		my $new_file = '';
		$new_file = $output_name."_split_".$file_num."_of_".$number_files."_aln.fas";
		unlink $new_file;
		open (MYFILE, '>>'.$new_file);
		
		my $seqs_per_file = roundup($num_sequences/$number_files);
		for(my $split_seq = 0; $split_seq < $seqs_per_file; $split_seq++) {
			if($seq_counter < $num_sequences) {
				
				print  MYFILE ">".$seq_array[$seq_counter]->display_id."\n";
				print  MYFILE $seq_array[$seq_counter]->subseq(1,$seq_array[$seq_counter]->length)."\n";
			}
			$seq_counter++;
		}
		close MYFILE;
	}
	print "Done!";
}


sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}










