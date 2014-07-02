#!/usr/bin/env perl
# Outputs taxa ID's for sequences based on delimiter position
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
use lib "$FindBin::Bin/libs";

use Sequence::ExemplarGenes;
use Getopt::Long;

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

my $fasta_input 		= '';
my $out_prefix			= 'out';

my $seq_length_cutoff	= 0;

GetOptions ("aln=s" 		=> \$fasta_input,
			"out=s"			=> \$out_prefix,
			"slength"		=> \$seq_length_cutoff,)
or die("Error in command line arguments\n");

my $seqio  = Bio::SeqIO->new(-file => $fasta_input, '-format' => 'Fasta');
my @seq_array = ();
while( my $seq = $seqio->next_seq() ) {
    push(@seq_array,$seq);
}

open (TAXAFILE, '>'.$out_prefix);
foreach my $seq (@seq_array) {
	my $seq_id = $seq->id;
	my @split_id = split(/\|/,$seq_id);
	next if fast_seq_length($seq->seq) < $seq_length_cutoff;
	
	foreach my $split (@split_id) {
		print TAXAFILE $split.",";
	}
	print TAXAFILE "\n";

}
close TAXAFILE;

print "Done!\n";

sub fast_seq_length {
	my $seq = shift;
	
	$seq =~ s/-/ /g;
	$seq =~ s/\s+//g;
	
	return length($seq);
}