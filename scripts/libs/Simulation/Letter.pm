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
class_has 'letter_distribution' => (
	is      => 'rw',
	isa     => 'Any',
	builder => '_create_letter_distribution',
);
sub _create_letter_distribution {
	return 'gaussian';
}

class_has 'mean_letter_radius' => (
	is      => 'rw',
	isa     => 'Int',
	builder => '_create_mean_letter_radius',
);
sub _create_mean_letter_radius {
	return 25;
}

class_has 'letter_stdev' => (
	is	=> 'rw',
	isa	=> 'Int',
	builder => '_create_letter_stdev',
);
sub _create_letter_stdev {
	return 3;
}

class_has 'mean_letter_num_centroids' => (
	is => 'rw',
	isa => 'Int',
	builder => '_create_mean_num_centroids',
);
sub _create_mean_num_centroids {
	return 25;
}

class_has 'letter_std_dev_rate' => (
	is	=> 'rw',
	isa	=> 'Any',
	builder => '_create_letter_std_dev_rate',
);
sub _create_letter_std_dev_rate {
	return 0.05;
}
########################################################################

########################################################################
# Attributes

has 'letter_grid' => (
	is => 'rw',
	isa => 'Simulation::Grid',
);

has 'letter_radius' => (
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_letter_radius',
);
sub _create_letter_radius {
	my $self = shift;
	
	my $letter_radius = 1;
	if(Simulation::Letter->letter_distribution eq 'gaussian') {
		$letter_radius =	int random_normal(	1,
							Simulation::Letter->mean_letter_radius, 
							Simulation::Letter->mean_letter_radius*$self->letter_std_dev_rate);
		if ($letter_radius > 0) {
			return $letter_radius;
		} else {
			return 1;
		}
	}
}

has 'letter_num_centroids' => (
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_letter_num_centroids',
);
sub _create_letter_num_centroids {
	my $self = shift;
	
	my $max_centroids = $self->letter_radius*$self->letter_radius; # Grid area.
	my $num_centroids = 1;
	if(Simulation::Letter->letter_distribution eq 'gaussian') {
		$num_centroids = 	int random_normal(	1, 
							Simulation::Letter->mean_letter_num_centroids, 
							Simulation::Letter->mean_letter_num_centroids*$self->letter_std_dev_rate);
		if($num_centroids > $max_centroids) {
			$num_centroids = $max_centroids;
		}
		if($num_centroids < 1) {
			$num_centroids = 1;
		}
		return $num_centroids;
	}
}

has 'letter_centroid_list' => (
	is => 'rw',
	isa => 'ArrayRef[Centroid]',
);

has 'letter_mass' => (
# Letter mass should be the sum of the steric radii of its Centroids
	is => 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_letter_mass',
);

sub _create_letter_mass {
# Before calling letter_mass, make sure you've counted centroids, and
# summed their radii.
# Letter mass is defined as the sum of the radii of steric portion of
# centroids.
	my $self = shift;
	my $centroids_ref = $self->find_centroids;
	# print $centroids_ref."\n";
	
	my $letter_mass = $self->calculate_letter_mass(centroid_ref => $self->find_centroids);
	return $letter_mass;
}
	
########################################################################

########################################################################
sub BUILD {
	my $self = shift;
	$self->print_to_logfile("Building Letter...");
	# The letter grid is a square with a radius input at creation time.
	$self->letter_grid(Simulation::Grid->new(xmax => $self->letter_radius, ymax => $self->letter_radius));
	$self->print_to_logfile("Centroids...".$self->letter_num_centroids);
	$self->print_to_logfile("\n");
	
	# Building the Letter consists of adding Centroids, and assigning those
	# Centroids properties from a selected distribution.
	for (my $i = 0; $i <= $self->letter_num_centroids; $i++) {
		my $random_point = $self->letter_grid->get_random_point;
		$random_point->add_to_bucket(Simulation::Centroid->new(max_centroid_radius => $self->letter_radius));
	}
	
}
########################################################################


sub find_centroids {
# Search through all the points on the letter grid
	my $self = shift;
	my @centroids = ();
	# Loop through each point and find centroids.
	my $found_centroids = 0;
	for (my $x = 0; $x < $self->letter_radius; $x++) {
        for (my $y = 0; $y < $self->letter_radius; $y++) {
			# get the point at x,y position
			my $current_point = $self->letter_grid->get_point( x => $x, y => $y);
			# loop through and check the bucket for Centroids
			foreach my $object ($current_point->bucket) {
				if (ref $object eq 'Simulation::Centroid') {
					push(@centroids, $object);
					$found_centroids++;
				}
			}
			last if $found_centroids == $self->letter_num_centroids;
		}
		last if $found_centroids == $self->letter_num_centroids;
	}
	# print ref \@centroids;
	# print "\n";
	# print ref @centroids;
	# print "\n";
	return \@centroids;
}

sub calculate_letter_mass {
# Loops through all Centroids and calculates the mass of the letter
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
		$cumulative_centroid_mass += $centroid->centroid_steric_radius;
	}
	return $cumulative_centroid_mass;
}


__PACKAGE__->meta->make_immutable;
1;