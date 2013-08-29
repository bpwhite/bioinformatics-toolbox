# Order specimen list by a text file of ID's
#
# Copyright (c) 2012, Bryan White

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

open(MASTER, "<master_specimen_order.txt");
my @master_list = <MASTER>;

open(CLUSTER, "<cluster_contents_baetis_500bp.txt");
my @cluster_list = <CLUSTER>;

my $output = 'ordered_specimen_list.txt';
unlink($output);

open (ORDERED, '>>'.$output);
print ORDERED "specimen_id,otu_id\n";

foreach my $master_id (@master_list) {
	# print $master_id."\n";
	chomp($master_id);
	my $specimen_counter = 0;
	my $was_found = 0;
	foreach my $cluster_specimen (@cluster_list) {
		if(($specimen_counter > 0)) {
			# print "\t".$cluster_specimen."\n";
			my @split_specimen_on_comma = split(/,/,$cluster_specimen);
			my @split_specimen_on_bar = split(/\|/,$split_specimen_on_comma[1]);
			my @split_cluster_on_bar = split(/\|/,$split_specimen_on_comma[0]);
			# print $split_specimen_on_comma[1]." \n ";
			# print "|".$split_specimen_on_bar[0]."| = |".$master_id."|\n";
			# foreach my $split(@split_specimen_on_bar) {
				# print $split."\n";
			# }
			if($split_specimen_on_bar[0] eq $master_id) {
				print ORDERED $split_specimen_on_bar[0].",".$split_cluster_on_bar[0]."|".$split_cluster_on_bar[1]."\n";
				$was_found = 1;
			}
			# foreach my $delimited (@delimited_cluster_specimen) {
				# print "\t\t".$delimited." => ";
			# }
			# print "\n";
			# exit;
		}
		$specimen_counter++;
		# exit;
	}
	if($was_found == 0) {
		print ORDERED $master_id.",NA\n";
	# exit;
	}
	
}




