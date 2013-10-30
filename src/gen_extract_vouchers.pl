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
										option_defs => {'-csv-files'		=> '',	# First gene file.
														'-voucher-column' 	=> '',	# CSV column containing voucher ID
														'-dna-column'		=> '',  # CSV column containing DNA sequence
													}
													);


my $csv_files	=  $params->options->{'-csv-files'};

my @csv_files = split(/,/,$csv_files);
print "Processing...\n";
foreach my $file (@csv_files) {
	print $file."\n";
	open (VFILE, '<'.$file);
	my @current_file = <VFILE>;
	close(VFILE);
	foreach my $line (@current_file) {
		print $line;
	}
	# my $rect1 = CollisionBlock->new( x =>  1, y =>  1, w => 10, h => 10 );
	my $vinfo = General::Voucher->new( voucher_id => 1, sequence => 2, genus => 3, species => 4);
	# print $vinfo->voucher_id."\n";
}
exit;


