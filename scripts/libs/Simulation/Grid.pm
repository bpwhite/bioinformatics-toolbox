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
package Simulation::Grid;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;
use Params::Validate qw(:all);

use Simulation::SimulationObject;
use Simulation::Point;

extends 'Simulation::SimulationObject';
########################################################################
# Class Attributes

########################################################################
# Attributes
has 'xmax' => (
# Max x value of the grid
	is => 'rw',
	isa => 'Int',
	required => 1,
	);
	
has 'ymax' => (
# Max y value of the grid
	is => 'rw',
	isa => 'Int',
	required => 1,
	);

# sub _create_grid_area {
	# my $self = shift;
	# return $self->xmax*$self->ymax;
# }

has 'point_array' => (
# Holds the Point objects for the Grid
	is => 'rw',
	isa => 'ArrayRef[Point]',
	builder => '_create_point_array',
	);
	
has 'hopper_list' => (
# A hash of points that have something in their hopper
	is => 'rw',
	isa => 'ArrayRef[Point]',
	builder => '_create_hopper_list',
	auto_deref => 1,	
	);
########################################################################

########################################################################
sub BUILD {
# Creates array of Points according to the input dimensions.
    my $self = shift;	
	Simulation::SimulationObject->print_to_logfile("Building grid ".$self->xmax." x ".$self->ymax."...");
	# print "Building grid ".$self->xmax." x ".$self->ymax."...\n";
	for (my $x = 0; $x < $self->xmax; $x++) {
        for (my $y = 0; $y < $self->ymax; $y++) {
			$self->point_array->[$x][$y] = Simulation::Point->new(x => $x, y => $y);
		}
	}
}
########################################################################

########################################################################
sub _create_dimensions {
# Builder for dimension attribute
	return [];
}

sub _create_point_array {
# Builder for point_array
	return [];
}

sub _create_hopper_list {
# Builder for hopper list
	return [];
}



sub grid_area {
	my $self = shift;
	return $self->xmax*$self->ymax;
}

sub add_to_hopper_list {
# Adds points that had their hoppers modified to an array.
	my $self = shift;
	my @additions = @_;
	push($self->hopper_list,@additions);
}

sub get_random_point {
# Returns a random Point object.
	my $self = shift;
	my $rand_x = rand int ($self->xmax);
	my $rand_y = rand int ($self->ymax);
	my $current_object = $self->point_array->[$rand_x][$rand_y];
	return $current_object;
}

sub get_random_coords {
# Returns a random coordinate pair
	my $self = shift;
	my $rand_x = int rand ($self->xmax);
	my $rand_y = int rand ($self->ymax);
	return ($rand_x,$rand_y);
}
sub get_point {
# Returns a single point object given a coordinate pair
	my $self = shift;
	my %params = validate(
		@_, {
			x => 1,
			y => 1,
		}
	);
	if(defined($self->point_array->[$params{'x'}][$params{'y'}])) {
		return $self->point_array->[$params{'x'}][$params{'y'}];
	} else {
		die("Tried to access a point that does not exist");
	}
}

sub transfer_to_hopper {
# Transfers the contents of a Point to the hopper another Point
# Deletes the contents of the original point
	my $self = shift;
	my %params = validate(
		@_, {
			x1 => 1, y1 => 1,
			x2 => 1, y2 => 1,
		}
	);
	my $current_point 	= $self->point_array->[$params{'x1'}][$params{'y1'}];
	my $new_point 		= $self->point_array->[$params{'x2'}][$params{'y2'}];
	$new_point->add_to_hopper($current_point->bucket);
	$self->add_to_hopper_list($new_point);
	$current_point->dump_bucket;
}

sub shuffle_hoppers {
# Dumps the buckets, adds the hoppers into buckets, then dumps the hoppers
	my $self = shift;
	foreach my $current_point ($self->hopper_list) {
		$current_point->add_to_bucket($current_point->hopper);
		$current_point->dump_hopper;
	}
	my @flush = ();
	$self->hopper_list(@flush);
}


__PACKAGE__->meta->make_immutable;
1;