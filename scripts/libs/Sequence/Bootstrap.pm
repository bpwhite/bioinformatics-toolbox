# Class for Bootstrap analysis of sequences
#
# Copyright (c) 2013, Bryan White

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
package Sequence::Bootstrap;

use strict;
use warnings;
use Math::Random::MT::Perl qw(srand rand);

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw(bootstrap_weights);


sub bootstrap_weights {
	my $seq_length = shift;
	
	# Create starting values
	my @unweighted = ();
	my @bs_start_weights = ();

	# Fill a starting weight array with 1's
	# Fill a starting weight array to be
	# bootstrapped with 0's
	my $default_weight = 1;
	my $max_chars = $seq_length-1;
	my $chars_sampled = 0;
	for (0..($max_chars)) {
		push(@unweighted, $default_weight);
		push(@bs_start_weights, 0);
	}
	my @bs_weights = @bs_start_weights;

	# Increase the to-be-returned weight array by
	# one for each time that position is sampled 
	for (0..($max_chars)) {
		my $rand_position = abs int rand ($seq_length);
		# if($rand_position == 653) { die "Found 653" };
		$bs_weights[$rand_position]++;
	}

	# Convert array to string
	my $weights = join("", @bs_weights);
	my $unweighted = join("",@unweighted);
	# Return string of weights
	return ($weights,$unweighted);
}

1;