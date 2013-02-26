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
package Simulation::Alphabet;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;
use Params::Validate qw(:all);

use Simulation::SimulationObject;
use Simulation::Letter;

extends 'Simulation::SimulationObject';

########################################################################
# Class Attributes
class_has 'mean_alphabet_size' => (
	is => 'rw',
	isa => 'Int',
	builder => '_create_mean_alphabet_size',
);

sub _create_mean_alphabet_size {
	return 25;
}
########################################################################

########################################################################
# Attributes
has 'letter_list' => (
	is => 'rw',
	isa => 'ArrayRef[Letter]',
	builder => '_create_letter_list',
	auto_deref => 1,
);
sub _create_letter_list {
	return [];
}
########################################################################

########################################################################
sub BUILD {
	my $self = shift;
	
	for (my $i = 0; $i < $self->mean_alphabet_size; $i++) {
		my $new_letter = Simulation::Letter->new();
		push($self->letter_list, $new_letter);
	}
	# $self->alphabet_stats;
}
########################################################################

########################################################################
sub alphabet_stats {
# Print some alphabet statistics
	my $self = shift;
	foreach my $letter ($self->letter_list) {
		print $letter->num_centroids." => ".$letter->mass."\n";
	}

}

sub alphabet_length {
	my $self = shift;
	my $alphabet_length = 0;
	foreach my $letter ($self->letter_list) {
		$alphabet_length++;
	}
	return $alphabet_length;
}
########################################################################

__PACKAGE__->meta->make_immutable;
1;