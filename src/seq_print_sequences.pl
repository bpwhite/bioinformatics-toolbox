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
use FindBin;
use lib "$FindBin::Bin/libs/Sequence"; 
use lib "$FindBin::Bin/libs/";
use General::Arguments;

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
	
	# my $was_repeat = 0;
	# if(!defined($seq_hash{$seqobj->display_id})) {
		# $seq_hash{$seqobj->display_id} = 0;
	# }
	# my $num_blanks = 0;
	# my $max_length = $seqobj->length;
	# for ( my $i = 1; $i <= $max_length; $i++) {
		# my $nuc = $seqobj->subseq($i,$i) or die("I died at ".$i);
		# if($nuc eq "-") {
			# $num_blanks++;
		# }
	# }
	# my $length = $seqobj->length-$num_blanks;
	# print "Seen sequence ",$seqobj->display_id,", start of seq ",substr($seqobj->seq,1,10)," length: ",$length,"\n";
	# print SEQS $seqobj->display_id."\n";
	my $new_id = filter_one_id($seqobj->display_id);
	$new_id = convert_id_to_name($new_id);
	print SEQS $new_id."\n";
	$num_seqs++;
}
close(SEQS);
print "Total seqs: ".$num_seqs."\n";

sub filter_one_id {
	my $id = shift;
	my $filtered_id = '';

	my @delimited_id = split(/\|/,$id);

	my $num_delimiters = 0;
	$num_delimiters = scalar(@delimited_id);

	# IF there are 5 bold delimited sections, here
	if($num_delimiters == 2) {
		my $new_id = $delimited_id[0]."|".$delimited_id[1];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	}
	elsif($num_delimiters == 5) {
		my $new_id = $delimited_id[2]."|".$delimited_id[3];
		# my $new_id = $delimited_id[2]."|".$delimited_id[3]."|".$delimited_id[4];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	# Only 4 bold delimited sections
	} elsif(($num_delimiters == 4) && ($delimited_id[3] eq "COI_5P")) {
		my $new_id = $delimited_id[1]."|".$delimited_id[2];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	} elsif(($num_delimiters) == 4 && ($delimited_id[2] ne "COI_5P")) {
		my $new_id = $delimited_id[2]."|".$delimited_id[3];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	} 
	elsif($num_delimiters == 3) {
		my $new_id = $delimited_id[1]."|".$delimited_id[2];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	}
	else {
		$filtered_id = $id;
	}
	return $filtered_id;
}


sub convert_id_to_name {
	# Returns a species id from an already filtered bold record
	# Input a node string id
	# Input: ID#|Name
	# Output: Name
	my ($node_id) = @_;
	my @delimited_id = ();
	@delimited_id = split(/\|/,$node_id);
	if(scalar(@delimited_id) == 2) {
		my $new_id = $delimited_id[1];
		return $new_id;
	} else {
		return $node_id;
	}
}