# Simply prints out sequences from a Fasta file
#
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
use 5.12.3;
use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

chomp(my $input_file = <>);
my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');

my %seq_hash = ();
my $num_seqs = 0;
my $num_unique_seqs = 0;

my $output_file = "sequence_list.txt";
unlink($output_file);
open(SEQS, '>>'.$output_file);
while((my $seqobj = $seqio->next_seq())) {
	
	my $was_repeat = 0;
	if(!defined($seq_hash{$seqobj->display_id})) {
		$seq_hash{$seqobj->display_id} = 0;
	}
	my $num_blanks = 0;
	my $max_length = $seqobj->length;
	for ( my $i = 1; $i <= $max_length; $i++) {
		my $nuc = $seqobj->subseq($i,$i) or die("I died at ".$i);
		if($nuc eq "-") {
			$num_blanks++;
		}
	}
	my $length = $seqobj->length-$num_blanks;
	print "Seen sequence ",$seqobj->display_id,", start of seq ",substr($seqobj->seq,1,10)," length: ",$length,"\n";
	print SEQS $seqobj->display_id."\n";
	$num_seqs++;
}
close(SEQS);
print "Total seqs: ".$num_seqs."\n";