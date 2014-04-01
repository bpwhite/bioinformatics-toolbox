#!/usr/bin/env perl
# This script cleans fasta ID's for bad characters

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

# Import sequence libs
use Sequence::Fasta;
use Sequence::Kimura_Distance;
use General::Arguments;
use Sequence::Bootstrap;
use Sequence::Garli;

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

# Non-bio modules
use Statistics::Descriptive;
use Benchmark qw(:all);
use POSIX qw(ceil);
use File::Copy;
use Math::Random::MT::Perl qw(srand rand irand);
use Digest::SHA qw(sha1 sha1_base64);
use Term::Spinner;
use Data::Dumper;
use Storable qw(dclone);
use File::Copy;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $alignment 	= '';
my $clean 		= '';
my $output 		= '';
my $site_pos	= 2;

GetOptions ("aln=s" 			=> \$alignment,
			"clean=s"			=> \$clean,
			"out=s"				=> \$output,
			"sitepos=s"			=> \$site_pos,)
or die("Error in command line arguments\n");



##################################################################
# Import alignment
print "Importing alignment file ".$alignment."...\n";
my $alignin = Bio::AlignIO->new(-format => 'fasta',
								-file   => $alignment );
my $original_aln = $alignin->next_aln;
##################################################################

my @starting_sequence_array = ();
foreach my $seq ($original_aln->each_seq) {
	push(@starting_sequence_array,$seq);
}

open (OUT, '>', $output);
foreach my $seq (@starting_sequence_array) {
	my $seq_id = $seq->id;

	$seq_text = $seq->seq;
	$seq_text =~ s/\-+$//;
	$seq->seq($seq_text);
	
	$seq_id =~ s/ /_/g;
	$seq->id($seq_id);
	print OUT ">".$seq->id."\n";
	print OUT $seq->seq."\n";
	
}
close (OUT);

