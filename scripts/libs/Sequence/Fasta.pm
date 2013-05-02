# General Fasta tools
#
# Copyright (c) 2013, Bryan White

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
use strict;
use warnings;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw( fix_fasta clean_file_name );

sub fix_fasta {
	my ( $f ) = shift;
	open F, "< $f" or die "Can't open $f : $!";
	my @fasta = <F>;
	close F;

	foreach my $line (@fasta) {
		if($line =~ /^>/) {
			$line =~ s/ /_/g; # replace whitespace with _
		}
	}
	
	unlink $f;
	open (MYFILE, '>>'.$f);
	foreach my $line(@fasta) {
		print MYFILE $line;
	}
	close(MYFILE);
}

sub clean_file_name {
	# Takes a filename and strips its extension off, and returns the file name
	my ($file_name) = @_;
	my @split_file_name = split(/\./,$file_name);
	my $new_file_name = $split_file_name[0];

	return $new_file_name;
}

sub alignment_coverage {
# Reduces an alignment to only positions that are covered by a certain percentage
	my $alignment 		= shift;
	my $coverage_pcnt 	= shift;
	
	
	my %position_counter = ();
	my $min_coverage = 0;
	my $current_seq = 0;
	foreach my $seq ($alignment->each_seq) {
		my $seq_string = $seq->seq;
		if($current_seq == 0) {
			print length($seq_string)."\n";
			$min_coverage = int length($seq_string)*$coverage_pcnt;
		}
		
		$seq_string = '---AAGGCCTT';
		my @unpacked_seq = unpack("C*", $seq_string);
		my $current_position = 0;
		foreach my $letter (@unpacked_seq) {
			# print $letter."\n";
			if($letter != 45) {
				$position_counter{$current_position}++;
			}
			$current_position++;
		}
		$current_seq++;
		# exit;
	}
	print $min_coverage."\n";
	my @valid_positions = ();
	for my $position ( sort keys %position_counter ) {
		print $position." => ".$position_counter{$position}."\n";
		if ($position_counter{$position} >= $min_coverage) {
			push(@valid_positions, $position);
		}
	}
	print "Reducing alignment to ".($coverage_pcnt*100)."% covered positions only.\n";
	print scalar(@valid_positions)." positions remain.\n";
	exit;
	return $alignment;
}


1;