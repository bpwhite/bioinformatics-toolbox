# Inputs a FASTA file and converts it to a nexus file.

# Copyright (c) 2011, Bryan White

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
use Bio::Seq;
use Bio::SeqIO;

use strict;
use warnings;

# my $alignment = 'baetis_test_aln_filtered_aln.fas';
chomp (my $alignment = <>);

my $seqio  = Bio::SeqIO->new(-file => $alignment , '-format' => 'Fasta');

my @delimited_alignment = split(".fas",$alignment);
my $new_alignment = $delimited_alignment[0]."_nex.nexus";

unlink $new_alignment;
open (NEXUS, '>>'.$new_alignment);

my @seq_array = ();
while((my $seqobj = $seqio->next_seq())) {
	push(@seq_array,$seqobj);
}
print NEXUS "#NEXUS\n";
print NEXUS "Begin taxa;\n";
print NEXUS "\tdimensions ntax=".scalar(@seq_array).";\n";
print NEXUS "\ttaxlabels\n";
my $aln_length = 0;
my $current_seq_i = 0;
foreach my $seq (@seq_array) {
	if ($current_seq_i == 0) {
		$aln_length = $seq->length();
	}
	print NEXUS $seq->id."\n";
	$current_seq_i++;
}
print NEXUS ";\n";
print NEXUS "end;\n";
print NEXUS "begin characters;\n";
print NEXUS "\tdimensions nchar=".$aln_length.";\n";
print NEXUS "\tformat missing=? gap=- matchchar=. datatype=dna;\n";
print NEXUS "\toptions gapmode=missing;\n";
print NEXUS "\tmatrix\n\n";
foreach my $seq (@seq_array) {
	# print NEXUS $seq->id()."\n";
	print NEXUS $seq->id()." ";
	print NEXUS $seq->subseq(1,$aln_length)."\n";
	# print NEXUS "\n";
}
print NEXUS ";\n";
print NEXUS "end;\n";