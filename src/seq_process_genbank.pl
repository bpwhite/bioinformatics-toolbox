#!/usr/bin/env perl
# Batch download taxa 
#
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
use FindBin;
use lib "$FindBin::Bin/libs";

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


my $genbank_csv		= '';
my $output_tag		= 'output';

GetOptions ("gb=s" 	=> \$genbank_csv,
			"out=s"		=> \$output_tag,)
or die("Error in command line arguments\n");

my @genbank_lines = ();
open (GB_LINES, '<'.$genbank_csv);
@genbank_lines = <GB_LINES>;
close(GB_LINES);

foreach my $line (@genbank_lines) {
	#print $line;
	$line =~ s/\n|\"//g;
	my @split_line = split(',',$line);
	foreach my $split (@split_line) {
		print $split." => ";
	}
	print "\n";
}



