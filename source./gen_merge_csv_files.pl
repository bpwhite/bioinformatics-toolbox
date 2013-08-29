#!/usr/bin/perl
# Merges all CSV files in a directory.
#
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
use lib "$FindBin::Bin/libs/";
use General::Arguments;

my $params = General::Arguments->new(	arguments_v => \@ARGV,
									option_defs => {'-folder' 		=> '', 				# Folder where CSV files are stored.
													'-outp'			=> 'combined.csv',	# File Output name.
													}
													);

my $folder = $params->options->{'-folder'};


my @csv_files = <$folder/*.csv>;
my @output_csv = ();
my $file_counter = 0;
foreach my $csv_file (@csv_files) {
	print "Adding $csv_file\n";
	open F, "< $csv_file";
	my @current_file = <F>;
	close F or die "Couln't close $csv_file!";
	my $line_counter = 0;
	foreach my $line (@current_file) {
		if(($file_counter != 0) && ($line_counter == 0)) { # skip first line
			$line_counter = 1;
			next;
		}
		push(@output_csv, $line);
		$line_counter++;
	}
	$file_counter++;
}

my $output_file = $params->options->{'-outp'};

unlink ($output_file);
open (OUTPUT, '>>'.$output_file) or die "Could not open $f";
foreach my $line(@output_csv) {
	print OUTPUT $line;
}
close(OUTPUT);

print 'Done combining!'."\n";