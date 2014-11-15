#!/usr/bin/env perl
# Combine processed genbank files 
#
# Copyright (c) 2014 Bryan White, bpcwhite@gmail.com

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

use Sequence::QualityCheck;
use Sequence::Fasta;
use Sequence::Kimura_Distance;

use Getopt::Long;

print "
******************************************************************

Copyright (c) 2013-2014 Bryan White, bpcwhite\@gmail.com

GNU General Public License, Version 3, 29 June 2007

This program comes with ABSOLUTELY NO WARRANTY; for details type.
This is free software, and you are welcome to redistribute it
under certain conditions. 

A full copy of the GPL 3.0 license should be accompanied with any 
distribution of this software.

******************************************************************
\n\n
";

my $combine_list	= '';
my $combine_columns	= 'output';
my $output_tag		= '';
GetOptions ("list=s" 	=> \$combine_list,
			"columns=s" 		=> \$combine_columns,
			"out=s"				=> \$output_tag,)
or die("Error in command line arguments\n");

my %overall_data_hash = ();
# Load the files to combine
my @combine_list = ();
open(COMBINE, '<'.$combine_list);
@combine_list = <COMBINE>;
close(COMBINE);
my @genbank_lines = ();
# Load the processed genbank CSV files
foreach my $gb_file (@combine_list) {
	my @file_lines = ();
	$gb_file =~ s/\n//g;
	print $gb_file."\n";
	open (GB_LINES, '<'.$gb_file) or die("could not open $_");
	@file_lines = <GB_LINES>;

	# Loop through Genbank CSV lines and output to a FASTA file
	my %params_hash = ();
	my $line_counter = 0;
	my $column_headers = '';
	foreach my $line (@file_lines) {
		#print $line;
		if($line_counter == 0) {
			$column_headers = $line;
		}
		$line =~ s/\n//g;
		$line =~ s/^\"|\"\,$//g;

		my @split_line = split('","',$line);

		if($line_counter == 0) {
			my $column_num = 0;
			foreach my $split (@split_line) {
				$params_hash{$split} = $column_num;
				#print $split." => ".$column_num."\n";
				$column_num++;
			}
			
		} else {
			# Build FASTA format sequences from each line
			my $column_num = 0;
			my $accession = '';
			my $nuc_seq = '';
			my $binom = '';
			foreach my $split (@split_line) {
				#print $column_num." => ".$split."\n";
				if($column_num == $params_hash{'nuc_seq'}) {
					$nuc_seq = $split;
				}
				if($column_num == $params_hash{'accession'}) {
					$accession = $split;
					$overall_data_hash{$accession}->{'gb_file'} = $gb_file;

				}
				if($column_num == $params_hash{'binom'}) {
					#print $params_hash{'binom'}."\n";
					#print $column_num."\n";
					$binom = $split;
					print $binom."\n";
				}

				$column_num++;
			}
		}

		$line_counter++;
	}










	push(@genbank_lines, @file_lines);
	close(GB_LINES);
}

foreach my $gb_line (@genbank_lines) {
	#print $gb_line."\n";
}



