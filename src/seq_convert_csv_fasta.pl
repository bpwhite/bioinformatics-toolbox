#!/usr/bin/env perl
# 
# Copyright (c) 2013, 2014 Bryan White, bpcwhite@gmail.com

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
use Getopt::Long;

my $csv_file = '';
my $output	 = 'out.fas';

GetOptions ("csv_file=s" 			=> \$csv_file,
			"output=s"				=> \$output,)
	or die("Error in command line arguments\n");
	
my @csv_lines = ();
open (CSVFILE, '<'.$csv_file);
@csv_lines = <CSVFILE>;
close(CSVFILE);

open (FASTA, '>'.$output);
foreach my $line (@csv_lines) {
	my @split_line = split(/\$/,$line);
	print FASTA ">".$split_line[0]."\n";
	print FASTA $split_line[1];
}
close (FASTA);

