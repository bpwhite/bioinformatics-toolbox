# Concatenates alignments, i.e. it appends a COI alignment to 16S alignment,
# so and so forth, for each of the genes. The final result is a concatenated
# mutli-gene alignment.

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
use Statistics::Descriptive;
use Math::Random::MT qw(srand rand);
use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

my $output_file = "combined_alignment.fas";
unlink $output_file;

my $id_position = 1;
my $replace_hyphens = 'no';

# Get the master list
my $master_alignment = "bopD_aln.fas";
my $master_seqio  = Bio::SeqIO->new(-file => $master_alignment , '-format' => 'Fasta');
my @master_list = ();
while((my $master_seqobj = $master_seqio->next_seq())) {
	push(@master_list, $master_seqobj->display_id);
	# print $master_seqobj->display_id."\n";
	
}

# Set the gene query files.
my @alignments = <*.fas>;
# Now that the list is built, open the new output file.
open (MYFILE, '>>'.$output_file);

my @sequence_position = ();
my %sequence_position_hash = ();

my $taxa_count = 0;
foreach my $master_order_id (@master_list) {
	print "Concatenating ".$master_order_id."\n";
	my $alignment_i = 1;
	foreach my $alignment (@alignments) {
		my $seqio  = Bio::SeqIO->new(-file => $alignment , '-format' => 'Fasta');
		my @delimited_master = split(/\|/,$master_order_id);

		my $found_sequence = 0; # Check if there was a sequence for that ID or not
		my $found_sequence_length = 0; # Use this to determine how many gaps to put
		my $current_gene = '';
		while((my $seqobj = $seqio->next_seq())) {
			my @delimited_current_id = split(/\|/,$seqobj->display_id);
			$current_gene = $delimited_current_id[0];
			if($replace_hyphens eq 'y') {
				$delimited_current_id[1] =~ s/-/_/;
			}
			if($seqobj->length > 0) {
				$found_sequence_length = $seqobj->length; # Get the length of this particular alignment
			}
			
			if($delimited_master[$id_position] eq $delimited_current_id[$id_position]) {
			# print $delimited_master[$id_position]." => ".$delimited_current_id[$id_position]."\n";
				if($alignment_i == 1) {
					if(defined($delimited_current_id[2])) {
						print "\nFound ".$alignment."\n";
						print MYFILE ">".$delimited_current_id[1]."|".$delimited_current_id[2]."\n";
						print  ">".$delimited_current_id[1]."|".$delimited_current_id[2]."\n";
					} else {
						print "\nFound ".$alignment."\n";
						print MYFILE ">".$delimited_current_id[1]."\n";
						# print  ">".$delimited_current_id[1]."\n";

					}
					print MYFILE $seqobj->subseq(1,$seqobj->length);
					# print  $seqobj->subseq(1,$seqobj->length);
					$found_sequence = 1; # Sequence was found
					last;
				} else {
					print "\nFound ".$alignment."\n";
					print MYFILE $seqobj->subseq(1,$seqobj->length);
					# print  $seqobj->subseq(1,$seqobj->length);
					$found_sequence = 1; # Sequence was found
					last;
				}
			}
			
		}
		
		if ($taxa_count == 0) {
			$sequence_position_hash{$current_gene} = $found_sequence_length;
			push(@sequence_position,$found_sequence_length);
		}
		# Sequence wasn't filled, fill it with gaps.
		if ($found_sequence == 0) {
			if($alignment_i == 1) {
				if(defined($delimited_master[2])) {
					print "Missing :".$alignment."\n";
					print MYFILE ">".$delimited_master[1]."|".$delimited_master[2]."\n";
					for(my $i = 1; $i <= $found_sequence_length; $i++) {
						print MYFILE "-";
					}
				} else {
					print "Missing ".$alignment."\n";
					print MYFILE ">".$delimited_master[1]."\n";
					for(my $i = 1; $i <= $found_sequence_length; $i++) {
						print MYFILE "-";
					}
				}
			}
		}
		$alignment_i++;
		# print "\nFound sequence? ".$found_sequence."\n";
	}
	print MYFILE "\n";
	
	$taxa_count++;

}

my $length_counter = 0;
while ( my ($gene, $current_seq_length) = each(%sequence_position_hash) ) {
	my $previous_length = $length_counter + 1;
	$length_counter = $current_seq_length + $length_counter;
	print "charset\t".$gene." = ".($previous_length)."-".$length_counter.";\n";
}

# foreach my $current_seq_length (@sequence_position) {
	# my $previous_length = $length_counter + 1;
	# $length_counter = $current_seq_length + $length_counter;
	# print " ".($previous_length)." - ".$length_counter."\n";
# }

# close(MYFILE);
# unlink(MYFILE);
