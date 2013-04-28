#!/usr/bin/perl
# Aligns two sequences using MAFFT
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
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::SeqIO::fasta;

use strict;
use warnings;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw( alignseqs );

sub alignseqs {
	my @seqs = @_;
	my $aln_file = 'alnseqs.fas';
	my $output_file = 'out_alnseqs.fas';
	
	unlink $aln_file;
	unlink $output_file;
	# unlink '1';
	open (ALN, '>'.$aln_file);
	my $seq_counter = 0;
	foreach my $seq (@seqs) {
		# print $seq."\n";
		print ALN ">Seq_$seq_counter\n$seq\n";
		$seq_counter++;
	}
	close ALN;
	
	# system('mafft '.$aln_file.' > out_'.$aln_file);
	# system("mafft --preservecase $aln_file > $output_file > nul 2>");
	system("mafft --quiet --retree 1 --maxiterate 0 --preservecase $aln_file > $output_file");

	my $seqio  = Bio::SeqIO->new(-file => $output_file, '-format' => 'Fasta');
	my @aligned_seqs;
	while( my $seq = $seqio->next_seq() ) {
		push(@aligned_seqs,$seq->seq);
	}
	
	unlink $aln_file;
	unlink $output_file;
	# unlink '1';
	return @aligned_seqs;
}