#!/usr/bin/env perl
# Batch download taxa 
#
# Copyright (c) 2014, Bryan White, bpcwhite@gmail.com

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

use Sequence::QualityCheck;
use Sequence::Fasta;
use Sequence::Kimura_Distance;
use General::Arguments;
use Sequence::Bootstrap;

# BioPerl libs
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::SimpleAlign;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;
use Benchmark qw(:all);

use Data::Dumper;

use Getopt::Long;

print "
******************************************************************

Copyright (c) 2013-2014 Bryan White, bpcwhite\@gmail.com

GNU General Public License, Version 3, 29 June 2007

This program comes with ABSOLUTELY NO WARRANTY; for details type.
This is free software, and you are welcome to redistribute it
under certain conditions. 

A full copy of the GPL 3.0 license should be accompanied with any 
distribution of this software.

******************************************************************
\n\n
";
##################################################################
# Start benchmark
my $t0 = Benchmark->new;
my $k2p1 = 0;
##################################################################

my $fasta_file		= '';
my $output_tag		= 'output';
my $subsample_size	= 100;
my $max_dist		= 0.5;
my $without_replacement = 0;
my $check_species 	= 0;
my $binom_pos		= 1; # Binom position from right
my $inter_intra		= 0; # inter vs. intraspecific distances: 0 = inter, 1 = intra
my $run_tag			= 'def';
my $print_headers	= 0;
my $alignment_length_cutoff = 350;
my $name_size_min	= 4;
my $binom_first		= 0;

GetOptions ("fasta=s" 			=> \$fasta_file,
			"output=s"			=> \$output_tag,
			"samples=s"			=> \$subsample_size,
			"maxdist=s"			=> \$max_dist,
			"replacement=s" 	=> \$without_replacement,
			"checkspecies=s" 	=> \$check_species,
			"binom-pos=s"		=> \$binom_pos,
			"inter-intra=s"		=> \$inter_intra,
			"run-tag=s"			=> \$run_tag,
			"print-headers=s"	=> \$print_headers,
			"binom-first=s"		=> \$binom_first)
or die("Error in command line arguments\n");

fix_fasta($fasta_file);

# Reimport aligned sequences and put in array
print "Loading... ".$fasta_file."\n";
my $alignment = Bio::AlignIO->new(-format => 'fasta',
								-file   => $fasta_file );
my $alignment_obj = $alignment->next_aln;
my @alignment_1 = ();
my @alignment_2 = ();

foreach my $seq ($alignment_obj->each_seq) {

	my $seq_degapped = $seq->seq;
	$seq_degapped =~ s/-/ /g;
	$seq_degapped =~ s/\s+//g;
	if(length($seq_degapped) < $alignment_length_cutoff) {
		next;
	}
	push(@alignment_1,$seq);
}
my $num_sequences = scalar(@alignment_1);

@alignment_2 = @{ dclone(\@alignment_1) };

fisher_yates_shuffle( \@alignment_1 );
fisher_yates_shuffle( \@alignment_2 );

print "Printing to ".$output_tag.".csv\n";
unlink($output_tag.".csv");
open(DISTFILE, '>>'.$output_tag.".csv");
if ($print_headers == 1) {
	print DISTFILE "run_tag,type,id1,id2,dist\n";
}

# Determine alignment length, including gaps
my $aln_length1 = length($alignment_1[0]->seq);
my $aln_length2 = length($alignment_2[0]->seq);
if($aln_length1 != $aln_length2) {
	die("Problem in copying arrays.\n");
}
my @exclude_names = ('sp.','sp','cf.','cf','sl','s.l.');

my $inter_intra_tag = '';
if($inter_intra == 0) {
	$inter_intra_tag = 'Interspecific';
} elsif( $inter_intra == 1) {
	$inter_intra_tag = 'Intraspecific';
}

my $dists_printed = 0;
my @pairs = ();
while ($dists_printed <= $subsample_size) {

	# Select 2 random sequences
	my $sample1 = rand( int( $num_sequences ));
	my $sample2 = rand( int( $num_sequences ));

	# Don't match a sequence to itself
	next if $sample1 == $sample2;
	
	# Defaults to with replacement
	if($without_replacement == 1) {
		my $pair = $sample1.'-'.$sample2;
		next if $pair ~~ @pairs;
		push(@pairs,$pair);
	}

	my $seq1 = $alignment_1[$sample1];
	my $seq2 = $alignment_2[$sample2];

	my $seq_id1 	= $seq1->id;
	my $seq_id2 	= $seq2->id;

	#next if $seq_id1 =~ m/Baetis/;
	#next if $seq_id2 =~ m/Baetis/;

	if($check_species == 1) {
		my @split_id1 = split(/\|/,$seq_id1);
		my @split_id2 = split(/\|/,$seq_id2);
		my @split_taxonid1 = split(/_/, $split_id1[-$binom_pos]);
		my @split_taxonid2 = split(/_/, $split_id2[-$binom_pos]);

		if((scalar(@split_taxonid1) > 1) && (scalar(@split_taxonid2) > 1)) {
			my $binom1 = $split_taxonid1[-2].' '.$split_taxonid1[-1];
			my $binom2 = $split_taxonid2[-2].' '.$split_taxonid2[-1];

			#next if $split_taxonid1[-1] ~~ @exclude_names;
			#next if $split_taxonid1[-2] ~~ @exclude_names;

			#next if $split_taxonid2[-1] ~~ @exclude_names;
			#next if $split_taxonid2[-2] ~~ @exclude_names;

			next if length($split_taxonid1[-1]) <  $name_size_min;
			next if length($split_taxonid1[-2]) <  $name_size_min;

			next if length($split_taxonid2[-1]) <  $name_size_min;
			next if length($split_taxonid2[-2]) <  $name_size_min;

			# Compute interspecific distances
			if($inter_intra == 0) {
				next if $binom1 ne $binom2;
			} elsif ($inter_intra == 1) {
				next if $binom1 eq $binom2;
			}
			#print $binom1." => ".$binom2."\n";
		} else {
			next;
		}
	}

	my $seq1_seq 	= $seq1->seq;
	my $seq2_seq 	= $seq2->seq;

	my ($k2p_distance, $transitions,$transversions,$bases_compared) = k2p_unpack($seq1_seq,$seq2_seq,$aln_length1);
	
	#if($k2p_distance > 0.2) {
	#	print $seq_id1." => ".$seq_id2." => ".$k2p_distance."\n";
	#}
	
	next if $k2p_distance >= $max_dist;
	#print $k2p_distance."\n";

	print DISTFILE $run_tag,",".$inter_intra_tag.",".$seq_id1.",".$seq_id2.",".$k2p_distance."\n";

	$dists_printed++;
}
close(DISTFILE);

my $t1 = Benchmark->new;
my $time_diff = timediff($t1, $t0);
print "\n";
print timestr($time_diff)."\n";