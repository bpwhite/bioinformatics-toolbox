#!/usr/bin/env perl
# This script converts DNA sequences to Amino Acid sequences.

# Copyright (c) 2013, 2014 Bryan White, bpcwhite@gmail.com

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

use strict;
use warnings;

# Import sequence libs
use Sequence::Fasta;
use Sequence::Kimura_Distance;
use General::Arguments;
use Sequence::Bootstrap;
use Sequence::Garli;

# BioPerl libs
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::SimpleAlign;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;

# Non-bio modules
use Statistics::Descriptive;
use Benchmark qw(:all);
use POSIX qw(ceil);
use Math::Random::MT::Perl qw(srand rand irand);
use Data::Dumper;
use Storable qw(dclone);
use File::Copy;

# Non-bio modules
use Benchmark qw(:all);
use Getopt::Long;
use String::Random;

my $query_aln_file 	= '';

GetOptions ("aln=s" 			=> \$query_aln_file,)
or die("Error in command line arguments\n");

fix_fasta($query_aln_file);

##################################################################
# Import alignment
print "Importing alignment file ".$query_aln_file."...\n";
my $query_aln = Bio::AlignIO->new(-format => 'fasta',
								-file   => $query_aln_file );
my $query_aln_obj = $query_aln->next_aln;
##################################################################

my @query_sequence_array = ();
foreach my $seq ($query_aln_obj->each_seq) {
	push(@query_sequence_array,$seq);
}

my $outp = '_AA.fas';

open (OUTP, '>'.$outp);
print "Translating proteins...\n\n";
foreach my $query_seq (@query_sequence_array) {
	
	my $seq_id = $query_seq->id;
	my $seq_string = $query_seq->seq;
	
	#print $seq_string."\n\n";

	#my $rev_complement = reverse_complement($seq_string);
	#print $rev_complement."\n\n";

	# my $rf1 = translate_protein($seq_string,1);
	# print $rf1."\n";

	#my $rf4 = translate_protein($rev_complement,1);
	#print $rf4."\n";

	#exit;
	#print $seq_id."\n";
	#print $seq_string."\n";
	
	my ($rf, $num_sc, $best_rf, $min_sc, $max_sc) = best_translation($seq_string);

	print $seq_id."\n";
	print $rf.":".$num_sc." => ".$best_rf."\n";

	#exit;
	# open (QUERY, '>'.$temp_query);
	# foreach my $match_seq (@match_aln_array) {
	# 	print QUERY '>Match|'.$match_seq->id."\n";
	# 	print QUERY $match_seq->seq."\n";
	# }
	# print QUERY ">".$seq_id."\n";
	# print QUERY $seq_string."\n";
	# close (QUERY);
	
	
}
print "Done!\n";

close(OUTP);

sub translate_protein {

	my $seq 	= shift;
	my $rfnum 	= shift;

	if ($rfnum > 3 || $rfnum < 1) {
		die("Invalid reading frame");
	}

	my %codons = (
		TTT => 'F',
		TTC => 'F',
		TTA => 'L',
		TTG => 'L',
		CTT => 'L',
		CTC => 'L',
		CTA => 'L',
		CTG => 'L',
		ATT => 'I',
		ATC => 'I',
		ATA => 'I',
		ATG => 'M',
		GTT => 'V',
		GTC => 'V',
		GTA => 'V',
		GTG => 'V',
		TCT => 'S',
		TCC => 'S',
		TCA => 'S',
		TCG => 'S',
		CCT => 'P',
		CCC => 'P',
		CCA => 'P',
		CCG => 'P',
		ACT => 'T',
		ACC => 'T',
		ACA => 'T',
		ACG => 'T',
		GCT => 'A',
		GCC => 'A',
		GCA => 'A',
		GCG => 'A',
		TAT => 'Y',
		TAC => 'Y',
		TAA => '*',
		TAG => '*',
		CAT => 'H',
		CAC => 'H',
		CAA => 'Q',
		CAG => 'Q',
		AAT => 'N',
		AAC => 'N',
		AAA => 'K',
		AAG => 'K',
		GAT => 'D',
		GAC => 'D',
		GAA => 'E',
		GAG => 'E',
		TGT => 'C',
		TGC => 'C',
		TGA => '*',
		TGG => 'W',
		CGT => 'R',
		CGC => 'R',
		CGA => 'R',
		CGG => 'R',
		AGT => 'S',
		AGC => 'S',
		AGA => 'R',
		AGG => 'R',
		GGT => 'G',
		GGC => 'G',
		GGA => 'G',
		GGG => 'G',
		"---" => "---"
	);

	my $length = length($seq);
	my $amino_acid_seq = '';
	for (my $i = $rfnum-1; $i < $length; $i+=3) {
		
		#print substr $seq, $i, $codon_end;
		my $codon = '';
		if(defined(substr $seq, $i, 3)) {
			$codon = substr $seq, $i, 3;
		}
		#print $codon."\n";
		#exit;
		#print $codons{$codon}."\n";
		my $amino_acid = "\?";
		if(exists($codons{$codon})) {
			$amino_acid = $codons{$codon};
		}
		#print $amino_acid;

		$amino_acid_seq .= $amino_acid;
	}

	return $amino_acid_seq;
}

sub reverse_complement {
	my $seq = shift;

	my %complementary_bases = (
		A => 'T',
		T => 'A',
		G => 'C',
		C => 'G',
		"-" => "-"
	);

	my $reverse = reverse $seq;
	my $length = length($reverse);
	#print $reverse."\n";

	my @rev_split = split(//,$reverse);
	#print $rev_split[0]."\n";

	my $rev_complement = '';
	for (my $i = 0; $i < $length; $i++) {
		my $base = $rev_split[$i];
		$rev_complement .= $complementary_bases{$base};
	}

	return $rev_complement;
}

sub best_translation {
	my $seq = shift;

	my %reading_frames = (
		1 => '',
		2 => '',
		3 => '',
		4 => '',
		5 => '',
		6 => ''
	);

	my $rev_complement = reverse_complement($seq);

	for(my $i = 1; $i <= 3; $i++) {
		#print $i."\n";
		$reading_frames{$i} = translate_protein($seq,$i);
		$reading_frames{$i+3} = translate_protein($rev_complement,$i);
	}

	my $max_stop_codons = 0;
	my $min_stop_codons = 999;
	my $most_stop_codons = '';
	my $least_stop_codons = '';
	while (my ($rf,$aaseq) = each (%reading_frames)) {
		my @stop_codons = $aaseq =~ /\*/g;
		my $num_stop_codons = scalar @stop_codons;
		my @start_codons = $aaseq =~ /M/g;
		my $num_start_codons = scalar @start_codons;
		#print scalar $num_stop_codons."\n";
		if($num_stop_codons > $max_stop_codons) {
			$max_stop_codons = $num_stop_codons;
			$most_stop_codons = $rf;
		}
		if($num_stop_codons < $min_stop_codons) {
			$min_stop_codons = $num_stop_codons;
			$least_stop_codons = $rf;
		}

		my @split_aa = split(//,$aaseq);
		my $open_frame = 0;
		my $num_orfs = 0;
		my $orf_length = 0;
		my @orf_lengths = ();
		foreach my $residue (@split_aa) {
			if($residue eq 'M') {
				$open_frame = 1;
			}
			if($open_frame == 1) {
				$orf_length++;
			}
			if(($open_frame == 1) && ($residue eq '*')) {
				$open_frame = 0;
				#print "Orf Length: ".$orf_length."\n";
				push(@orf_lengths, $orf_length);
				$orf_length = 0;
				$num_orfs++;
			}
		}

		my $rounding = "%.1f";
		my $orf_length_stat = Statistics::Descriptive::Full->new();
		$orf_length_stat->add_data(@orf_lengths);
		#print scalar @orf_lengths."\n";
		my $mean_orf_length = $orf_length_stat->mean();
		# MAX #
		###
		###
		
		#print $mean_orf_length."\n";
		$mean_orf_length = sprintf($rounding, $orf_length_stat->mean());
		#print $rf." => ".$aaseq."\n";
		print $rf." => ".$num_stop_codons." => ".$num_start_codons." => ".$num_orfs." => ".$mean_orf_length."\n";
	}
	print "Least stop codons:\n";
	print $min_stop_codons.": RF: ".$least_stop_codons."\n";
	print "Most stop codons\n";
	print $max_stop_codons.": RF: ".$most_stop_codons."\n";

	return ($least_stop_codons,
			$min_stop_codons,
			$reading_frames{$least_stop_codons},
			$min_stop_codons,
			$max_stop_codons);
}


