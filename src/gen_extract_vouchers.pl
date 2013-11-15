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
use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/libs/";
use General::Arguments;
use General::Voucher;
use Sequence::ExemplarGenes;

my $params = General::Arguments->new(	arguments_v => \@ARGV,
										option_defs => {'-csv-file'			=> '',				# Gene file.
														'-voucher-col' 		=> '',				# CSV column containing voucher ID
														'-dna-col'			=> '',  			# CSV column containing DNA sequence
														'-species-col'		=> '',  			# species column
														'-gene-name'		=> 'gene',			# Gene name
														'-outp'				=> 'output.fas',	# Output file
														'-skip'				=> '1',				# Number of lines to skip on import
														'-gene-exemplar'	=> '',				# Select an exemplar gene from file
														'-validate'			=> '0',				# Validate sequences against an exemplar gene
														'-seq-type'			=> 'nucl',			# Type of sequence
													}
													);
													
my $csv_file		= $params->options->{'-csv-file'};
my $voucher_col 	= $params->options->{'-voucher-col'}	- 1;
my $dna_col 		= $params->options->{'-dna-col'} 		- 1;
my $species_col 	= $params->options->{'-species-col'} 	- 1;
my $output_file 	= $params->options->{'-outp'};
my $gene_name		= $params->options->{'-gene-name'};
my $skip_lines		= $params->options->{'-skip'};
my $gene_exemplar 	= $params->options->{'-gene-exemplar'};

my $validate		= $params->options->{'-validate'};
my $seq_type	 	= $params->options->{'-seq-type'};

my $gene_exemplar_seq = '';
my $file_separator = "\\";
my $detected_os = 'Win32';

if("$^O\n" =~ "Win32") {
	print "Detected Windows\n";
} else {
	print "Detected Linux\n";
	$file_separator = "/";
	$detected_os = 'linux';
}

if($gene_exemplar eq '18S') {
	$gene_exemplar_seq = Sequence::ExemplarGenes->print_18S();
} elsif($gene_exemplar eq '16S') {
	$gene_exemplar_seq = Sequence::ExemplarGenes->print_16S();
} elsif($gene_exemplar eq '28S') {
	$gene_exemplar_seq = Sequence::ExemplarGenes->print_28S();
} elsif($gene_exemplar eq 'COI') {
	$gene_exemplar_seq = Sequence::ExemplarGenes->print_COI();
}

print "Processing...\n";
print $csv_file."\n";
open (VFILE, '<'.$csv_file);
my @current_file = <VFILE>;
close(VFILE);

my %voucher_hash = ();
open (OUTP, '>'.$output_file);
my $line_i = 1;
if($gene_exemplar_seq ne '') {
	print OUTP $gene_exemplar_seq."\n";
}

my $output_path = $output_file.'_temp';
unless(-d ($output_path)) {
	mkdir $output_path;
}

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
		$split =~ s/\'//g;
		$split =~ s/ /_/g;
	}
	if(exists($voucher_hash{$split_line[$voucher_col]})) {
		$line_i++;
		next;
	} else {
		$voucher_hash{$split_line[$voucher_col]} = 1;
	}
	
	my $output_string = ">".$gene_name."|".$split_line[$voucher_col]."|".$split_line[$species_col]."\n".$split_line[$dna_col]."\n";
	
	my $temp_name = $split_line[$voucher_col]."_".;
	print "Creating output folder ".$output_path.$file_separator.$output_file."\n";
	open (BLDB, '>'.$output_path.$file_separator.$output_file);
	print BLDB $gene_exemplar_seq."\n";
	print BLDB $output_string;
	close (BLDB);
	system("makeblastdb -in ".$output_path.$output_file." -dbtype nucl");

	exit;

	# print $output_string;
	print OUTP $output_string;
	$line_i++;
}
close(OUTP);
# my $vinfo = General::Voucher->new( voucher_id => 1, sequence => 2, genus => 3, species => 4);
# print $vinfo->voucher_id."\n";