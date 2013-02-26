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
package Simulation::SimulationInstance;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;

use Simulation::SimulationObject;
use Simulation::Point;
use Simulation::Grid;
use Simulation::Alphabet;
use Simulation::Letter;
use Params::Validate qw(:all);

extends 'Simulation::SimulationObject';
########################################################################
# Class Attributes


########################################################################
# Attributes
has 'grid_xmax' => (
	is => 'rw',
	isa => 'Int',
	required => 1,
	);

has 'grid_ymax' => (
	is => 'rw',
	isa => 'Int',
	required => 1,
	);
	
has 'primary_grid' => (
	is => 'rw',
	isa => 'Simulation::Grid',
	);

has 'current_alphabet' => (
	is => 'rw',
	isa => 'Simulation::Alphabet',
	);
	
has 'alphabet_states' => (
	is => 'rw',
	isa => 'HashRef[Simulation::Alphabet]',
	);
########################################################################

########################################################################
sub BUILD {
	my $self = shift;
	
	####################################################################
	# Build the primary grid for the simulation.
	$self->primary_grid(Simulation::Grid->new(	xmax => $self->grid_xmax, 
												ymax => $self->grid_ymax)
												);

	$self->print_to_logfile("\n");								
	#####################################################################
	
	#####################################################################
	# Build an alphabet for this simulation, or import an already built one.
	if (!defined($self->current_alphabet)) {
		$self->current_alphabet(Simulation::Alphabet->new());
	}
	#####################################################################
	
	#####################################################################
	# Sample the Alphabet and place Letters on the Simulation Grid.
	# 
	#####################################################################
}
########################################################################

########################################################################
sub step_grid {
	my $self = shift;
	$self->primary_grid->shuffle_hoppers;
}

sub point_mover_test {
	my $self = shift;
	my $test1 = $self->primary_grid->get_point( x => 5, y => 10);

	$test1->add_to_bucket("A","B","C");
	print "1\n";
	foreach my $thing ($test1->bucket) {
		print $thing."\n";
	}

	$self->primary_grid->transfer_to_hopper(x1 	=> 5, 	y1 	=> 10,
											x2 	=> 10, 	y2 	=> 15);
	$self->primary_grid->shuffle_hoppers;
	
	print "2\n";
	foreach my $thing ($test1->bucket) {
		print $thing."\n";
	}

	my $test2 = $self->primary_grid->get_point( x => 10, y => 15);
	print "3\n";
	foreach my $thing ($test2->bucket) {
		print $thing."\n";
	}

	print $test2->point_id."\n";
}

__PACKAGE__->meta->make_immutable;
1;