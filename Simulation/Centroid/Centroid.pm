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
package Centroid;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;

extends 'SimulationObject';
########################################################################
# Class Variables
class_has 'centroid_distribution' => (
	is	=> 'rw',
	isa => 'Any',
);

class_has 'centroid_std_dev_rate' => (
	is => 'rw',
	isa => 'Any',
	builder => '_create_centroid_std_dev_rate',
);
sub _create_centroid_std_dev_rate {
	return 0.05;
}

# Steric
class_has 'mean_steric_radius' => (
	is      => 'rw',
	isa     => 'Int',
	lazy	=> 1,
	builder => '_create_mean_steric_radius',
	);
sub _create_mean_steric_radius {
	return 0.25;
}

# Electrostatic
class_has 'mean_electrostatic_radius' => (
	is		=> 'rw',
	isa		=> 'Int',
	lazy	=> 1,
	builder => '_create_mean_electrostatic_radius',
);
sub _create_mean_electrostatic_radius {
	return 0.25;
}

# Hydrogen
class_has 'mean_hydrogen_radius' => (
	is		=> 'rw',
	isa		=> 'Int',
	lazy	=> 1,
	builder => '_create_mean_hydrogen_radius',
);
sub _create_mean_hydrogen_radius {
	return 0.25;
}

# Van der Waals
class_has 'mean_vanderwaals_radius' => (
	is		=> 'rw',
	isa		=> 'Int',
	lazy	=> 1,
	builder => '_create_mean_vanderwaals_radius',
);
sub _create_mean_vanderwaals_radius {
	return 0.25;
}

# Covalent
class_has 'mean_covalent_radius' => (
	is		=> 'rw',
	isa		=> 'Int',
	lazy	=> 1,
	builder => '_create_mean_covalent_radius',
);
sub _create_mean_covalent_radius {
	return 0.25;
}

########################################################################
# Attributes

has 'max_centroid_radius' => ( 
# Centroid force radii are always a function of the max centroid radius,
# which should be the letter radius for this centroid
	is => 'rw',
	isa => 'Int',
	required => 1,
	);
	
# Steric
has 'centroid_steric_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_centroid_steric_radius',
	);
sub _create_centroid_steric_radius {
	my $self = shift;
	my $steric_radius = 1;
	if(Centroid->centroid_distribution eq 'gaussian') {
		$steric_radius =	int random_normal(	1,
							Centroid->mean_steric_radius*$self->max_centroid_radius, 
							Centroid->mean_steric_radius*$self->centroid_std_dev_rate);
	}
}

# Electrostatic
has 'centroid_electrostatic_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_centroid_electrostatic_radius',
	);
sub _create_centroid_electrostatic_radius {
	my $self = shift;
	my $steric_radius = 1;
	if(Centroid->centroid_distribution eq 'gaussian') {
		$steric_radius =	int random_normal(	1,
							Centroid->mean_electrostatic_radius*$self->max_centroid_radius, 
							Centroid->mean_electrostatic_radius*$self->centroid_std_dev_rate);
	}
}

# Hydrogen
has 'centroid_hydrogen_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_centroid_hydrogen_radius',
	);
sub _create_centroid_hydrogen_radius {
	my $self = shift;
	my $steric_radius = 1;
	if(Centroid->centroid_distribution eq 'gaussian') {
		$steric_radius =	int random_normal(	1,
							Centroid->mean_hydrogen_radius*$self->max_centroid_radius, 
							Centroid->mean_hydrogen_radius*$self->centroid_std_dev_rate);
	}
}

# Van der Waals
has 'centroid_vanderwaals_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_centroid_vanderwaals_radius',
	);
sub _create_centroid_vanderwaals_radius {
	my $self = shift;
	my $steric_radius = 1;
	if(Centroid->centroid_distribution eq 'gaussian') {
		$steric_radius =	int random_normal(	1,
							Centroid->mean_vanderwaals_radius*$self->max_centroid_radius, 
							Centroid->mean_vanderwaals_radius*$self->centroid_std_dev_rate);
	}
}

# Covalent
has 'centroid_covalent_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_centroid_covalent_radius',
	);
sub _create_centroid_covalent_radius {
	my $self = shift;
	my $steric_radius = 1;
	if(Centroid->centroid_distribution eq 'gaussian') {
		$steric_radius =	int random_normal(	1,
							Centroid->mean_covalent_radius*$self->max_centroid_radius, 
							Centroid->mean_covalent_radius*$self->centroid_std_dev_rate);
	}
}


__PACKAGE__->meta->make_immutable;
1;