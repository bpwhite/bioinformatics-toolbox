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
										option_defs => {'-csv-file'			=> '',				# Gene file.
														'-voucher-col' 		=> '',				# CSV column containing voucher ID
														'-dna-col'			=> '',  			# CSV column containing DNA sequence
														'-species-col'		=> '',  			# species column
														'-gene-name'		=> 'gene',			# Gene name
														'-outp'				=> 'output.fas',	# Output file
														'-skip'				=> '1',				# Number of lines to skip on import
													}
													);

my $csv_file	= $params->options->{'-csv-file'};
my $voucher_col = $params->options->{'-voucher-col'}	- 1;
my $dna_col 	= $params->options->{'-dna-col'} 		- 1;
my $species_col = $params->options->{'-species-col'} 	- 1;
my $output_file = $params->options->{'-outp'};
my $gene_name	= $params->options->{'-gene-name'};
my $skip_lines	= $params->options->{'-skip'};

print "Processing...\n";
print $csv_file."\n";
open (VFILE, '<'.$csv_file);
my @current_file = <VFILE>;
close(VFILE);

my %voucher_hash = ();
open (OUTP, '>'.$output_file);
my $line_i = 1;
foreach my $line (@current_file) {
	if ($line_i <= $skip_lines) {
		$line_i++;
		next;
	}
	my @split_line = split(/,/,$line);
	foreach my $split (@split_line) {
		$split =~ s/\"//g;
		$split =~ s/^\s+//;
		$split =~ s/\s+$//;
		$split =~ s/\n//;
		$split =~ s/\'//;
		$split =~ s/ /_/g;
	}
	if(exists($voucher_hash{$split_line[$voucher_col]})) {
		$line_i++;
		next;
	} else {
		$voucher_hash{$split_line[$voucher_col]} = 1;
	}
	my $output_string = ">".$gene_name."|".$split_line[$voucher_col]."|".$split_line[$species_col]."\n".$split_line[$dna_col]."\n";
	
	# print $output_string;
	print OUTP $output_string;
	$line_i++;
}
close(OUTP);
# my $vinfo = General::Voucher->new( voucher_id => 1, sequence => 2, genus => 3, species => 4);
# print $vinfo->voucher_id."\n";