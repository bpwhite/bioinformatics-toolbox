#!/usr/bin/perl
# This script takes an input of a FASTA sequence and outputs a 
# species delimitation regime given certain parameters.
# It is designed for use with the cytochrome oxidase I (COI) gene
# and to provide species delimitations useful for DNA barcoding.

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
use FindBin;
use lib "$FindBin::Bin/libs/Sequence"; 
use lib "$FindBin::Bin/libs/";
use General::Arguments;



my $params = General::Arguments->new(	arguments_v => \@ARGV,
										option_defs => {'-check-list' 		=> '', 	# List file name
														'-delim-summary'	=> '',	# Delim summary file
														'-delim-results'	=> '',	# Delim contents
														'-delim-aligment'	=> '', 	# Delim alignment
													}
													);
my $check_list 				= $params->options->{'-check-list'};
my $delim_summary_file 		= $params->options->{'-delim-summary'};
my $delim_results_file 		= $params->options->{'-delim-results'};
my $delim_alignment_file 	= $params->options->{'-delim-aligment'};

open (TAXALIST, '<'.$check_list);
my @check_list_lines = <TAXALIST>;
close(TAXALIST);

open (DELIM, '<'.$delim_summary_file);
my @delim_summary_lines = <DELIM>;
close(DELIM);

# open (DELIM, '<'.$delim_results_file);
# my @delim_results_lines = <DELIM>;
# close(DELIM);

# open (DELIM, '<',.$delim_alignment_file);
# my @delim_lines = <DELIM>;
# close(DELIM);

clean_file_lines(\@check_list_lines);
clean_file_lines(\@delim_summary_lines);
# clean_file_lines(\@delim_results_lines);

my %unique_check_list = ();
foreach my $line (@check_list_lines) {
	$unique_check_list{$line} = 'A';
}

print "Taxa,".$delim_summary_lines[0]."\n";
for my $taxa ( sort keys %unique_check_list ) {
	# print $taxa."\n";
	my $found_delim = 0;
	my $delim_fields = 0;
	# my @current_correspondences = grep { $_ == $unique_otu_correspondences } @otu_morpho_correspondence;
	my @current_summary = grep { $_ =~ m/$taxa/ } @delim_summary_lines;
	foreach my $delim_summary (@delim_summary_lines) {
		my @split_summary = split(",",$delim_summary);
		$delim_fields = scalar(@split_summary);
		if($delim_summary =~ m/$taxa/) {
			if($split_summary[1] =~ m/$taxa/) {
				print "\t$taxa,";
				print $delim_summary."\n";
				# foreach my $summary_field (@split_summary) {
					# print $summary_field." => ";
				# }
				# print "\n";
				# print "\t".$split_summary[1]."\n";
				# print "\t".$delim_summary."\n";
			}
			$found_delim = 1;
		}
	}
	if($found_delim == 0) {
		# print "\tNot found.\n";
		print "\t$taxa,";
		for(0..($delim_fields-1)) {
			print "NA,";
		}
		print "\n";
	}
	# my @current_results = grep { $_ == $taxa } @delim_results_lines;
}

sub clean_file_lines {
	my $array_ref = shift;
	foreach my $line (@$array_ref) {
		$line =~ s/\n//g; # replace newlines
		$line =~ s/ /_/g; # replace whitespace
	}
}
