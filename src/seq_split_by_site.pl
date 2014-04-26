#!/usr/bin/env perl
# Resamples an alignment or splits it into even parts
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
my $use_coi_exemplar 	= 0;
my $out_prefix			= 'out';

GetOptions ("aln=s" 				=> \$fasta_input,
			"coiex=s"				=> \$use_coi_exemplar,
			"out=s"			=> \$out_prefix,)
or die("Error in command line arguments\n");

my $seqio  = Bio::SeqIO->new(-file => $fasta_input, '-format' => 'Fasta');
my @seq_array = ();
while( my $seq = $seqio->next_seq() ) {
    push(@seq_array,$seq);
}

foreach my $seq (@seq_array) {
	my $seq_id = $seq->id;
	my @split_id = split(/\|/,$seq_id);
	my $site_code = $split_id[4];
	my $accession = $split_id[1];
	my $read_id	  = $split_id[0];
	my $identity  = $split_id[3];
	my $taxa	  = $split_id[2];
	$seq_id = $read_id."|".$accession."|".$identity."|".$taxa;

	my $current_output = $out_prefix."_".$site_code.".fas";
	unless (-e $current_output) {
		open (MYFILE, '>'.$current_output);
		if($use_coi_exemplar == 1) {
			print MYFILE Sequence::ExemplarGenes::print_COI()."\n";
		}
		close MYFILE;
	} else {
		open (MYFILE, '>>'.$current_output);
		print  MYFILE ">".$seq_id."\n";
		print  MYFILE $seq->seq()."\n";
		close MYFILE;
	}
}


print "Done!";