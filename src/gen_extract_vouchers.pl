#!/usr/bin/perl
# extract voucher sequences for multigene analysis
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

#!/usr/bin/perl
use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/libs/";
use General::Arguments;
use General::Voucher;


my $params = General::Arguments->new(	arguments_v => \@ARGV,
										option_defs => {'-csv-file'			=> '',	# Gene file.
														'-voucher-col' 		=> '',	# CSV column containing voucher ID
														'-dna-col'			=> '',  # CSV column containing DNA sequence
														'-species-col'		=> '',  # species column
														'-outp'				=> 'output.fas',	# Output file
													}
													);

my $csv_file	= $params->options->{'-csv-file'};
my $voucher_col = $params->options->{'-voucher-col'}	- 1;
my $dna_col 	= $params->options->{'-dna-col'} 		- 1;
my $species_col = $params->options->{'-species-col'} 	- 1;
my $output_file = $params->options->{'-outp'};

print "Processing...\n";
print $csv_file."\n";
open (VFILE, '<'.$csv_file);
my @current_file = <VFILE>;
close(VFILE);

open (OUTP, '>'.$output_file);
foreach my $line (@current_file) {
	my @split_line = split(/,/,$line);
	foreach my $split (@split_line) {
		$split =~ s/\"//g;
		$split =~ s/^\s+//;
		$split =~ s/\s+$//;
		$split =~ s/ /_/g;
	}
	my $output_string = ">".$split_line[$voucher_col]."|".$split_line[$species_col]."\n".$split_line[$dna_col]."\n";
	
	print $output_string;
	print OUTP $output_string;
}
close(OUTP);
# my $vinfo = General::Voucher->new( voucher_id => 1, sequence => 2, genus => 3, species => 4);
# print $vinfo->voucher_id."\n";