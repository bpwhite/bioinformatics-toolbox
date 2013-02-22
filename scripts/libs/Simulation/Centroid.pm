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
package Simulation::Centroid;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;

use Math::Random qw(:all);

extends 'Simulation::SimulationObject';
########################################################################
# Class Variables
class_has 'distribution' => (
	is	=> 'rw',
	isa => 'Any',
);

class_has 'std_dev_rate' => (
	is => 'rw',
	isa => 'Any',
	builder => '_create_std_dev_rate',
);
sub _create_std_dev_rate {
	return 0.05;
}

# Steric
class_has 'mean_steric_radius' => (
	is      => 'rw',
	isa     => 'Any',
	lazy	=> 1,
	builder => '_create_mean_steric_radius',
	);
sub _create_mean_steric_radius {
	return 0.25;
}

# Electrostatic
class_has 'mean_electrostatic_radius' => (
	is		=> 'rw',
	isa		=> 'Any',
	lazy	=> 1,
	builder => '_create_mean_electrostatic_radius',
);
sub _create_mean_electrostatic_radius {
	return 0.25;
}

# Hydrogen
class_has 'mean_hydrogen_radius' => (
	is		=> 'rw',
	isa		=> 'Any',
	lazy	=> 1,
	builder => '_create_mean_hydrogen_radius',
);
sub _create_mean_hydrogen_radius {
	return 0.25;
}

# Van der Waals
class_has 'mean_vanderwaals_radius' => (
	is		=> 'rw',
	isa		=> 'Any',
	lazy	=> 1,
	builder => '_create_mean_vanderwaals_radius',
);
sub _create_mean_vanderwaals_radius {
	return 0.25;
}

# Covalent
class_has 'mean_covalent_radius' => (
	is		=> 'rw',
	isa		=> 'Any',
	lazy	=> 1,
	builder => '_create_mean_covalent_radius',
);
sub _create_mean_covalent_radius {
	return 0.25;
}

########################################################################
# Attributes

has 'max_radius' => ( 
# Centroid force radii are always a function of the max centroid radius,
# which should be the letter radius for this centroid
	is => 'rw',
	isa => 'Int',
	required => 1,
	);
	
# Steric
has 'steric_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_steric_radius',
	);
sub _create_steric_radius {
	my $self = shift;
	my $steric_radius = 1;
	if(Simulation::Centroid->distribution eq 'gaussian') {
		$steric_radius =	int random_normal(	1,
							Simulation::Centroid->mean_steric_radius*$self->max_radius, 
							Simulation::Centroid->mean_steric_radius*$self->std_dev_rate);
	}
	if($steric_radius < 1) { return 1 };
	return $steric_radius;
}

# Electrostatic
has 'electrostatic_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_electrostatic_radius',
	);
sub _create_electrostatic_radius {
	my $self = shift;
	my $electrostatic_radius = 1;
	if(Simulation::Centroid->distribution eq 'gaussian') {
		$electrostatic_radius =	int random_normal(	1,
							Simulation::Centroid->mean_electrostatic_radius*$self->max_radius, 
							Simulation::Centroid->mean_electrostatic_radius*$self->std_dev_rate);
	}
	if($electrostatic_radius < 1) { return 1 };
	return $electrostatic_radius;
}

# Hydrogen
has 'hydrogen_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_hydrogen_radius',
	);
sub _create_hydrogen_radius {
	my $self = shift;
	my $hydrogen_radius = 1;
	if(Simulation::Centroid->distribution eq 'gaussian') {
		$hydrogen_radius =	int random_normal(	1,
							Simulation::Centroid->mean_hydrogen_radius*$self->max_radius, 
							Simulation::Centroid->mean_hydrogen_radius*$self->std_dev_rate);
	}
	if($hydrogen_radius < 1) { return 1 };
	return $hydrogen_radius;
}

# Van der Waals
has 'vanderwaals_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_vanderwaals_radius',
	);
sub _create_vanderwaals_radius {
	my $self = shift;
	my $vanderwaals_radius = 1;
	if(Simulation::Centroid->distribution eq 'gaussian') {
		$vanderwaals_radius =	int random_normal(	1,
							Simulation::Centroid->mean_vanderwaals_radius*$self->max_radius, 
							Simulation::Centroid->mean_vanderwaals_radius*$self->std_dev_rate);
	}
	if($vanderwaals_radius < 1) { return 1 };
	return $vanderwaals_radius;
}

# Covalent
has 'covalent_radius' => (
	is	=> 'rw',
	isa => 'Int',
	lazy => 1,
	builder => '_create_covalent_radius',
	);
sub _create_covalent_radius {
	my $self = shift;
	my $covalent_radius = 1;
	if(Simulation::Centroid->distribution eq 'gaussian') {
		$covalent_radius =	int random_normal(	1,
							Simulation::Centroid->mean_covalent_radius*$self->max_radius, 
							Simulation::Centroid->mean_covalent_radius*$self->std_dev_rate);
	}
	if($covalent_radius < 1) { return 1 };
	return $covalent_radius;
}


__PACKAGE__->meta->make_immutable;
1;