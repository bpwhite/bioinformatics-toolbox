# Condenses haplotypes based on haplotype and location so that remaining sequences
# either differ by location or haplotype.
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

# Must require these functions to work.
require "$FindBin::Bin/libs/Sequence/Kimura_Distance_C.pl";

# use Bioinformatics::General;
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
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

# my $input_file = '18s_mega_aln.fas';
print "Input file name (.fas or .fasta): ";
chomp(my $input_file = <>);
fix_fasta($input_file);

my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');
my @alignment = ();
my @orig_alignment = ();
my %orig_seq_hash = ();
my %unique_seqs = ();
while((my $seq = $seqio->next_seq())) {
	$alignment_length = $seq->length;
	my @delimited_seq = split(/\|/,$seq->id);
	my $sample_id = $delimited_seq[0];
	my $tax_id = $delimited_seq[1];
	my $location = $delimited_seq[2];
	my $filtered_seq = $seq->seq;
	$filtered_seq =~ s/-/*/g;
	$seq->description($location);
	$seq->object_id($filtered_seq);
	$seq->primary_id($sample_id."|".$tax_id);
	push(@alignment,$seq);
	push(@orig_alignment,$seq);
	$orig_seq_hash{$seq->seq} = $seq;
	$unique_seqs{$seq->seq} = $seq->id;
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
# for my $seq1 ( sort keys %unique_seqs1) {
foreach my $seq (@alignment) {
	my $was_found = 0;
	my $was_in_distincts = 0;
	# print $seq1;
	# for my $seq2 ( sort keys %unique_seqs2) {
		# next if $unique_seqs1{$seq1} eq $unique_seqs2{$seq2};
		# my ($transitions,	$transversions,		$bases_compared,
			# $k2p_distance,	$variance,			$stderror,
			# $mink2p,		$maxk2p,			$p_stderror,
			# $p_min,			$p_max,				$p_dist
			# ) = 0;
		# my $search_type = 3;
		# my $cutoff = 0.02;
		# my $filtered_seq = $seq1;
		# $filtered_seq =~ s/-/*/g;
		# c_kimura_distance(	$filtered_seq,		$seq2,		$critical_value,
							# $cutoff, 			$search_type, 		$alignment_length,
							# $transitions,		$transversions,		$bases_compared,
							# $k2p_distance,		$variance,			$stderror,
							# $mink2p,			$maxk2p);
		# if($k2p_distance <= 0) {
			# if(exists($unique_seqs1{$seq1})) { delete $unique_seqs1{$seq1} };
			# if(exists($unique_seqs2{$seq1})) { delete $unique_seqs2{$seq1} };
			
			# last;
		# }
		# if($k2p_distance <= 0) {
			# print "A\n";
			# $was_found = 1;
		# }
		# last if $matrix_count == $counter_i;
		# $counter_i++;
	# }
	for my $distinct_seqs (sort keys %distinct_hash) {
		my ($transitions,	$transversions,		$bases_compared,
			$k2p_distance,	$variance,			$stderror,
			$mink2p,		$maxk2p,			$p_stderror,
			$p_min,			$p_max,				$p_dist
			) = 0;
		my $search_type = 3;
		my $cutoff = 0.02;
		my $filtered_seq = $seq->seq;
		$filtered_seq =~ s/-/*/g;
		c_kimura_distance(	$filtered_seq,		$distinct_seqs,		$critical_value,
							$cutoff, 			$search_type, 		$alignment_length,
							$transitions,		$transversions,		$bases_compared,
							$k2p_distance,		$variance,			$stderror,
							$mink2p,			$maxk2p);
		if($k2p_distance <= 0) {
			$was_in_distincts = 1;
		}
	}
	if($was_in_distincts == 0) {
		# print "D\n";
		$distinct_hash{$seq->seq} = 'a';
	}
	# print $was_found."\n";
	# print $was_in_distincts."\n";
	# print (keys %distinct_hash)."\n";
	# last if $matrix_count == $counter_i;
}
# %distinct_hash = %unique_seqs1;
print "Distincts: ".(keys %distinct_hash)."\n";
# exit;
my $test = 0;
my $haplo_loc_counter = 1;
my $haplo_abundance_counter = 0;
my @condensed_sequences = ();
for my $distinct_hap ( sort keys %distinct_hash ) {
	my $current_haplotype_id = $orig_seq_hash{$distinct_hap}->primary_id;
	# print "[".$haplo_loc_counter."] ".$current_haplotype_id."\n";
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
		my $search_type = 3;
		my $cutoff = 0.02;
		c_kimura_distance(	$filtered_orig_seq,	$distinct_hap,		$critical_value,
							$cutoff, 			$search_type, 		$alignment_length,
							$transitions,		$transversions,		$bases_compared,
							$k2p_distance,		$variance,			$stderror,
							$mink2p,			$maxk2p);
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
		# print "\t".$haplo_loc_counter."|".$current_haplotype_id."|".$location."|".scalar(@members)."\n";
		print CONDENSED ">".$haplo_loc_counter."|".$current_haplotype_id."|".$location."|".scalar(@members)."\n";
		print CONDENSED $distinct_hap."\n";
		$haplo_abundance_counter += scalar(@members);
	}
	$haplo_loc_counter++;
}

# my %unique_condensed_sequences = map {$_->id,1} @condensed_sequences;
# my @unique_condensed_sequences = keys %unique_condensed_sequences;
# foreach my $unique_condensed (@unique_condensed_sequences) {
	# my @sub_condensed = grep { $_->id eq $unique_condensed } @condensed_sequences;
	# if(scalar(@sub_condensed) > 1) {
		# print scalar(@sub_condensed)." ".$sub_condensed[0]->id."\n";
	# }
# }

print "Original Sequences:\t".$num_seqs."\n";
print "Condensed Sequences:\t".$haplo_abundance_counter."\n";
print "Distinct haplo_locs:\t".(keys %distinct_hash)."\n";
close(CONDENSED);

sub fix_fasta {
	my ( $f ) = shift;
	open F, "< $f" or die "Can't open $f : $!";
	my @fasta = <F>;
	close F;
	
	for(my $i = 0; $i < 10; $i++) {
		foreach my $line (@fasta) {
			if($line =~ /^>/) {
				$line =~ s/ /_/; # replace whitespace with _
			}
		}
	}

	unlink $f;
	open (MYFILE, '>>'.$f);
	foreach my $line(@fasta) {
		print MYFILE $line;
	}
	close(MYFILE);
}


sub clean_file_name {
	# Takes a filename and strips its extension off, and returns the file name
	my ($file_name) = @_;
	my @split_file_name = split(/\./,$file_name);
	my $new_file_name = $split_file_name[0];

	return $new_file_name;
}