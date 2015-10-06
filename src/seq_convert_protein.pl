#!/usr/bin/env perl
# This script converts DNA sequences to Amino Acid sequences.

# Copyright (c) 2013-2015 Bryan White, bpcwhite@gmail.com

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
use Sequence::Protein;
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
	
	my ($rf, $num_sc, $best_rf, $min_sc, $max_sc, $orfs) = best_translation($seq_string);

	print $seq_id."\n";
	#print $rf.":".$num_sc." => ".$best_rf."\n";

	foreach my $orf (@$orfs) {
		print $orf."\n";
	}
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


