#!/usr/bin/env perl
# Condenses haplotypes based on haplotype and location so that remaining sequences
# either differ by location or haplotype.
# 
# Copyright (c) 2013, Bryan White, bpw.informatics.com@gmail.com

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

# Import sequence libs
use Sequence::Fasta;
use Sequence::Kimura_Distance;
use General::Arguments;
use Sequence::Bootstrap;
use Getopt::Long;

# Local modules


use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

# print "    condense.pl  Copyright (c) 2013, Bryan White, bpcwhite\@gmail.com

	# Condenses a FASTA file to distinct haplotype and/or locale.
	# See Readme for instructions.
	
    # This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    # This is free software, and you are welcome to redistribute it
    # under certain conditions; type `show c' for details.\n\n";
my $critical_value = 1.96; # Default
my $alignment_length = 0;
my $input_file = '';
GetOptions ("in=s" 			=> \$input_file,)
or die("Error in command line arguments\n");

fix_fasta($input_file);

my @ambiguous_characters = ('N', 'R', 'Y', 'K', 
							'M', 'S', 'W', 'B', 
							'D', 'H');
							
my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');
my @alignment = ();
my @orig_alignment = ();
my %orig_seq_hash = ();
my %unique_seqs = ();
my %unpacked_sequences = ();
my %unpacked_filtered_sequences = ();

while((my $seq = $seqio->next_seq())) {
	my $seq_gapped = $seq->seq;
	$alignment_length = $seq->length;
	foreach my $ambiguous_character (@ambiguous_characters) {
		$seq_gapped =~ s/$ambiguous_character/-/g;
	}
	my @delimited_seq = split(/\|/,$seq->id);
	my $sample_id = $delimited_seq[0];
	# my $tax_id = $delimited_seq[1];
	# my $location = $delimited_seq[2];
	my $filtered_seq = $seq->seq;
	$filtered_seq =~ s/-/*/g;
	# $seq->description($location);
	$seq->object_id($filtered_seq);
	# $seq->primary_id($sample_id."|".$tax_id);
	push(@alignment,$seq);
	push(@orig_alignment,$seq);
	$orig_seq_hash{$seq->seq} = $seq;
	$unique_seqs{$seq->seq} = $seq->id;
	
	# Unpack sequences
	my @unpacked_seq = unpack("C*", $seq_gapped);
	$unpacked_sequences{$seq_gapped} = \@unpacked_seq;
	my @unpacked_filtered_seq = unpack("C*", $filtered_seq);
	$unpacked_filtered_sequences{$filtered_seq} = \@unpacked_filtered_seq;
}
my $num_seqs = scalar(@alignment);

my $condensed = clean_file_name($input_file)."_condensed_aln.fas";

unlink $condensed;
open (CONDENSED, '>>'.$condensed);

my %distinct_hash = ();
my $matrix_count = ($num_seqs*$num_seqs-$num_seqs)/2;
my $counter_i = 1;
my %unique_seqs1 = %unique_seqs;
my %unique_seqs2 = %unique_seqs;
print "Uniques: ".(keys %unique_seqs1)."\n";

# Add sequences that differ by some genetic distance to a hash.
# foreach my $seq (@alignment) {
	# my $was_found = 0;
	# my $was_in_distincts = 0;
	# for my $distinct_seqs (sort keys %distinct_hash) {
		# my ($transitions,	$transversions,		$bases_compared,
			# $k2p_distance,	$variance,			$stderror,
			# $mink2p,		$maxk2p,			$p_stderror,
			# $p_min,			$p_max,				$p_dist
			# ) = 0;
		# my $search_type = 3;
		# my $cutoff = 0.02;
		# my $filtered_seq = $seq->seq;
		# $filtered_seq =~ s/-/*/g;
		# if(!defined($distinct_hash{$seq->seq})) {
			# print "a\n";
			# $distinct_hash{$seq->seq} = 'a';
			# last;
		# }
		# ($k2p_distance, $transitions,$transversions,$bases_compared) = k2p_no_bs(	\$unpacked_sequences{$distinct_seqs}, 
																					# \$unpacked_filtered_sequences{$filtered_seq}, 
																					# $alignment_length);
		# print $filtered_seq."\n";
		# print $distinct_seqs."\n";
		# print $unpacked_filtered_sequences{$filtered_seq}->[1]."\n";
		# print $unpacked_sequences{$distinct_seqs}->[1]."\n";
		# print $alignment_length."\n";
		# exit;
		# if($k2p_distance <= 0) {
			# $was_in_distincts = 1;
		# }
	# }
	# if($was_in_distincts == 0) {
		# $distinct_hash{$seq->seq} = 'a';
	# }
# }

print "Distincts: ".(keys %unique_seqs)."\n";

my $test = 0;
my $haplo_loc_counter = 1;
my $haplo_abundance_counter = 0;
my @condensed_sequences = ();
# Loop through distincts and count how many sequences match each haplotype
# 1. Loop through distinct hash
# 2. For each distinct haplotype, find all the locations that haplotype is at
# 3. Count the number of haplotypes for each location
# 4. Print a copy of each haplotype for each sequence, it's haplotype number
# and the abundance of that haplotype at a particular location
for my $distinct_hap ( sort keys %unique_seqs ) {
	my $current_haplotype_id = $orig_seq_hash{$distinct_hap}->primary_id;
	my @haplotype_members = ();
	my @haplotype_locations = ();
	foreach my $orig_seq ( @orig_alignment ) {
		my $filtered_orig_seq = $orig_seq->object_id;
		my $orig_seq_location = $orig_seq->description;
		my ($transitions,	$transversions,		$bases_compared,
			$k2p_distance,	$variance,			$stderror,
			$mink2p,		$maxk2p,			$p_stderror,
			$p_min,			$p_max,				$p_dist
			) = 0;
		my $cutoff = 0.02;
		# print $distinct_hap."\n";
		# print $filtered_orig_seq."\n";
		# print $unpacked_sequences{$distinct_hap}->[0]."\n";
		# print $unpacked_filtered_sequences{$filtered_orig_seq}->[0]."\n";
		# exit;
		# die if !defined($unpacked_filtered_sequences{$filtered_orig_seq}->[0]);
		($k2p_distance, $transitions,$transversions,$bases_compared) = k2p_no_bs(	\$unpacked_sequences{$distinct_hap}, 
																					\$unpacked_filtered_sequences{$filtered_orig_seq}, 
																					$alignment_length);
		# print $k2p_distance."\n";
		# exit;
		if($k2p_distance <= 0) {
			next if $orig_seq->id ~~ @condensed_sequences; # Skip found seqs
			push(@condensed_sequences, $orig_seq->id);
			push(@haplotype_members,$orig_seq);
			push(@haplotype_locations,$orig_seq_location);
		}
	}
	my %unique_haplotype_locations = map {$_,1} @haplotype_locations;
	my @unique_haplotype_locations = keys %unique_haplotype_locations;
	foreach my $location (@unique_haplotype_locations) {
		my @members = grep { $_->description eq $location } @haplotype_members;
		print CONDENSED ">".$haplo_loc_counter."|".$current_haplotype_id."|".$location."|".scalar(@members)."\n";
		print CONDENSED $distinct_hap."\n";
		$haplo_abundance_counter += scalar(@members);
	}
	$haplo_loc_counter++;
}

print "Original Sequences:\t".$num_seqs."\n";
print "Condensed Sequences:\t".$haplo_abundance_counter."\n";
print "Distinct haplo_locs:\t".(keys %distinct_hash)."\n";
close(CONDENSED);

