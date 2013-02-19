# Copyright (c) 2012, Bryan White, bpcwhite@gmail.com

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
# fix_fasta($input_file);
# my $input_file = 'COI.fas';
my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');
my @alignment = ();
while((my $seq = $seqio->next_seq())) {
	push(@alignment,$seq);
}

my $new_file = "frame_shifted.fas";
print "Creating output file...\n";
unlink $new_file;
open(FRAME, '>>'.$new_file);

my %stop_codons = (
		'TAA' => 'STOP',
		'TAG' => 'STOP'
		# 'TGA' => 'STOP'
	);

my $num_check_stops = 0;
foreach my $seq (@alignment) {
	my $frame_shift = 0;
	
	check_frame:
	my $seq_string = $seq->seq();
	my $seq_length = length($seq_string);
	for(my $i = 0; $i <= $seq_length; $i++) {
		my $seq_base1 = substr($seq_string,$i,1);
		# print $seq_base1."\n";
		last if $i+1 > $seq_length;
		my $seq_base2 = substr($seq_string,($i+1),1);
		# print $seq_base2."\n";
		last if $i+2 > $seq_length;
		my $seq_base3 = substr($seq_string,($i+2),1);
		# print $seq_base3."\n";
		$i = $i+2;
		my $triplet = $seq_base1.$seq_base2.$seq_base3;
		# print "[".$i."] ".$triplet." +".$frame_shift."\n";
		my $check_stop = '';
		if(defined($stop_codons{$triplet})) {
			$check_stop = $stop_codons{$triplet};
		}
		if($frame_shift > 2) {
			# print $triplet."\n";
			print "Could not find frame for ".$seq->id."\n";
			$num_check_stops++;
			last;
		}
		if($check_stop eq 'STOP') {
			$frame_shift++;
			# print $check_stop."\n";
			$seq_string = "-".$seq_string;
			$seq->seq($seq_string);
			# print "YES";
			$check_stop = '';
			goto check_frame;
		}
	}
	# if($num_check_stops > 1) {
		# goto multiple_stop_codons;
	# }
	# my $gaps_string = '';
	# for(my $i = 0; $i<$frame_shift; $i++){
		# my $gaps_string .= '-';
	# }
	print FRAME ">".$seq->id."\n";
	print FRAME $seq->seq()."\n";
# multiple_stop_codons: # Skip the sequence it has more than 1 stop codon
	# exit;
}

close(FRAME);

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