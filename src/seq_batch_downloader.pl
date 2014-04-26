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

my $taxa_list	 	= '';
my $out_prefix		= 'out';

GetOptions ("list=s" 		=> \$taxa_list,
			"out=s"		=> \$out_prefix,)
or die("Error in command line arguments\n");

my @taxa_list = ();
open (TAXALIST, '<'.$taxa_list);
@taxa_list = <TAXALIST>;
close(TAXALIST);

foreach my $taxa (@taxa_list) {
	$taxa =~ s/\n//g; # replace newlines
	$taxa =~ s/ /_/g; # replace whitespace
	my $output_file = $out_prefix."_".$taxa;
	my $search_string = "seq_convert_genbank.pl -term $taxa -query COI -outp $output_file";
	print `$search_string`;
}