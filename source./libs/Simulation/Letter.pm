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
package Simulation::Letter;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;

use Math::Random qw(:all);
 
use Simulation::SimulationObject;
use Simulation::Point;
use Simulation::Grid;
use Simulation::Centroid;

use Params::Validate qw(:all);

extends 'Simulation::SimulationObject';

########################################################################
# Class Attributes
class_has 'distribution' => (
	is      => 'rw',
	isa     => 'Any',
	builder => '_create_distribution',
);
sub _create_distribution {
	return 'gaussian';
}

class_has 'mean_radius' => (
	is      => 'rw',
	isa     => 'Int',
	builder => '_create_mean_radius',
);
sub _create_mean_radius {
	return 25;
}

class_has 'stdev_rate' => (
	is	=> 'rw',
	isa	=> 'Any',
	builder => '_create_stdev_rate',
);
sub _create_stdev_rate {
	return 0.05;
}

class_has 'mean_num_centroids' => (
	is => 'rw',
	isa => 'Int',
	builder => '_create_mean_num_centroids',
);
sub _create_mean_num_centroids {
	return 25;
}

class_has 'centroid_density' => (
	is => 'rw',
	isa => 'Int',
);
########################################################################

########################################################################
# Attributes

has 'grid' => (
	is => 'rw',
	isa => 'Simulation::Grid',
);

has 'radius' => (
# Stores the letter radius
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_radius',
);
sub _create_radius {
# Generates a random letter radius.
	my $self = shift;
	
	my $radius = 1;
	if(Simulation::Letter->distribution eq 'gaussian') {
		$radius =	int random_normal(	1,
							Simulation::Letter->mean_radius, 
							Simulation::Letter->mean_radius*$self->stdev_rate);
		if ($radius > 0) {
			return $radius;
		} else {
			return 1;
		}
	}
}

has 'num_centroids' => (
# Stores the number of centroids
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_num_centroids',
);
sub _create_num_centroids {
# Generates a random number of centroids based on the mean number of centroids.
	my $self = shift;
	
	my $max_centroids = $self->radius*$self->radius; # Grid area.
	my $num_centroids = 1;
	if(Simulation::Letter->distribution eq 'gaussian') {
		$num_centroids = 	int random_normal(	1, 
							Simulation::Letter->mean_num_centroids, 
							Simulation::Letter->mean_num_centroids*$self->stdev_rate);
		if($num_centroids > $max_centroids) {
			$num_centroids = $max_centroids;
		}
		if($num_centroids < 1) {
			$num_centroids = 1;
		}
		return $num_centroids;
	}
}

has 'centroid_list' => (
	is => 'rw',
	isa => 'ArrayRef[Centroid]',
);

has 'mass' => (
# Letter mass should be the sum of the steric radii of its Centroids
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_mass',
);

sub _create_mass {
# Before calling mass, make sure you've counted centroids, and
# summed their radii.
# Letter mass is defined as the sum of the radii of steric portion of
# centroids.
	my $self = shift;
	my $centroids_ref = $self->find_centroids;
	# print $centroids_ref."\n";
	
	my $mass = $self->calculate_mass(centroid_ref => $self->find_centroids);
	return $mass;
}
	
########################################################################

########################################################################
sub BUILD {
	my $self = shift;
	####################################################################
	# Create a Grid for the letter that will hold Centroids
	$self->print_to_logfile("Building Letter...");
	$self->grid(Simulation::Grid->new(	xmax => $self->radius, 
										ymax => $self->radius));
	$self->print_to_logfile("Centroids...".$self->num_centroids);
	$self->print_to_logfile("\n");
	####################################################################
	
	####################################################################
	# Building the Letter consists of adding Centroids, and assigning those
	# Centroids properties from a selected distribution.
	my $previous_point = '';
	my $previous_radius = '';
	for (my $i = 0; $i <= $self->num_centroids; $i++) {
		if($i == 0) {
			# Place the first Centroid randomly
			my $first_point = $self->grid->get_random_point;
			my $new_centroid = Simulation::Centroid->new(max_radius => $self->radius);
			$previous_radius = $new_centroid->steric_radius;
			$first_point->add_to_bucket($new_centroid);
			$previous_point = $first_point;			
		} else {
			# Place the remaining Centroids 
			for (my $j = 0; $j < 9999999; $j++) {
				my @new_centroid_coords = $self->grid->get_random_coords;
				my $x = $new_centroid_coords[0];
				my $y = $new_centroid_coords[1];
				my $distance = $previous_point->euclidean_distance(	x2 => $x, y2 => $y);
				if($distance < $previous_radius) {
					my $new_centroid = Simulation::Centroid->new(max_radius => $self->radius);
					$previous_radius = $new_centroid->steric_radius;
					my $new_point = $self->grid->get_point(x => $x, y => $y);
					$new_point->add_to_bucket($new_centroid);
					goto placed_centroid; # Place the next Centroid.
				}
			}
			die "Could not place new Centroid.\n";
		}
		placed_centroid:

	}
	####################################################################
}
########################################################################

########################################################################
sub find_centroids {
# Search through all the points on the letter grid and return a ref array
# of centroids
	my $self = shift;
	my @centroids = ();
	# Loop through each point and find Centroids.
	my $found_centroids = 0;
	for (my $x = 0; $x < $self->radius; $x++) {
        for (my $y = 0; $y < $self->radius; $y++) {
			# get the point at x,y position
			my $current_point = $self->grid->get_point( x => $x, y => $y);
			# loop through and check the bucket for Centroids
			foreach my $object ($current_point->bucket) {
				if (ref $object eq 'Simulation::Centroid') {
					push(@centroids, $object);
					$found_centroids++;
				}
			}
			last if $found_centroids == $self->num_centroids;
		}
		last if $found_centroids == $self->num_centroids;
	}
	return \@centroids;
}

sub calculate_mass {
# Loops through all Centroids and calculates the mass of the letter
# Mass of the Letter is the sum of all centroid steric radii
# This is analagous to the sum of protons in an attom (atomic mass)
# Or the molecular mass of a molecule (
	my $self = shift;
	my %params = validate(
		@_, {
			centroid_ref => 1,
		}
	);
	my $centroid_ref = $params{'centroid_ref'};
	my @centroid_list = @$centroid_ref;
	my $cumulative_centroid_mass = 0;
	# loop through each centroid and sum its steric radius.
	foreach my $centroid( @centroid_list) {
		$cumulative_centroid_mass += $centroid->steric_radius;
		if($centroid->steric_radius == 0) { exit };
	}
	return $cumulative_centroid_mass;
}


sub rd_pt_within_bounds {
# Returns a random point within the bounds of another point.
# Typically the bounds should be the steric raidius of another
# Centroid.
	my $self = shift;
	my %params = validate(
		@_, {
			xmax => 1,
			ymax => 1,
			current_pt => 1,
		}
	);
	while(my $pt_in_range == 1) {
		my $new_point = $self->grid->get_random_point;
	}
	print $self->xmax."\n";
	print $self->ymax."\n";
	
	
	return 1;
}

########################################################################
__PACKAGE__->meta->make_immutable;
1;