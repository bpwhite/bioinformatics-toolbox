# Copyright (c) 2011, Bryan White, bpcwhite@gmail.com

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
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;
use strict;
use warnings;

my $alignment_file = 'Crepipatella_aligned.fas';
# chomp (my $alignment_file = <>);
# $alignment_file = fix_fasta($alignment_file);

my $alignment = Bio::AlignIO->new(-format => 'fasta',
								-file   => $alignment_file );
my $original_aln = $alignment->next_aln;

my @delimited_alignment = split(".fas",$alignment_file);
my $new_alignment = $delimited_alignment[0]."_renamed.fas";


unlink $new_alignment;
open (RENAMED, '>>'.$new_alignment);

my @species_names = ();
my @overall_sequences = ();
my %sequence_hash = ();
foreach my $seq ($original_aln->each_seq) {
	if(defined($sequence_hash{$seq->seq})) {
		$sequence_hash{$seq->seq} = $sequence_hash{$seq->seq} + 1;
	} else {
		$sequence_hash{$seq->seq} = 1;
		push(@overall_sequences,$seq);
	}
}

foreach my $seq (@overall_sequences) {
	$seq->description($sequence_hash{$seq->seq});
}

foreach my $seq (@overall_sequences) {
	
	my $seq_id = $seq->id;
	$seq_id =~ s/\|/_/g;
	my @split_seq_id = split(/\_/,$seq_id);
	# foreach my $split (@split_seq_id) {
		# print $split."\n";
	# }
	my $new_id = $split_seq_id[5]."_".$split_seq_id[6]."_".$split_seq_id[3];
	print ">".$new_id."\n";
	print $seq->seq()."\n";
	print RENAMED ">".$new_id."\n";
	print RENAMED $seq->seq()."\n";
	# $seq->id($new_id);
	# push(@species_names,$split_seq_id[1]."_".$split_seq_id[2]);
	# exit;
	
}
close RENAMED;


sub fix_fasta {
	my ( $f ) = shift;
	open F, "< $f" or die "Can't open $f : $!";
	my @fasta = <F>;
	close F;
	
	for(my $i = 0; $i < 10; $i++) {
		foreach my $line (@fasta) {
			if($line =~ /^>/) {
				$line =~ s/ /_/; # replace whitespace with _
			}
		}
	}

	unlink $f;
	open (MYFILE, '>>'.$f);
	foreach my $line(@fasta) {
		print MYFILE $line;
	}
	close(MYFILE);
}
