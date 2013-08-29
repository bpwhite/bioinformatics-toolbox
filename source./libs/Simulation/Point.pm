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

package Simulation::Point;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;
use Params::Validate qw(:all);

extends 'Simulation::SimulationObject';
########################################################################
# Class Variables


########################################################################
# Attributes
has 'x' => (
	is => 'rw',
	isa => 'Int',
	);
has 'y' => (
	is => 'rw',
	isa => 'Int',
	);
has 'bucket' => (
	is => 'rw',
	isa => 'ArrayRef[Any]',
	builder => '_create_bucket',
	required => 0,
	auto_deref => 1,
	);
sub _create_bucket {
	return [];
}

has 'hopper' => (
	is => 'rw',
	isa => 'ArrayRef[Any]',
	builder => '_create_hopper',
	required => 0,
	auto_deref => 1,
	);
sub _create_hopper {
	return [];
}
########################################################################

########################################################################
sub BUILD {
	my $self = shift;
}
########################################################################

########################################################################
sub add_to_bucket {
	my $self = shift;
	my @additions = @_;
	push($self->bucket,@additions);
}

sub dump_bucket {
	my $self = shift;
	my @flush = ();
	$self->bucket(@flush);
}

sub add_to_hopper {
	my $self = shift;
	my @additions = @_;
	push($self->hopper,@additions);
}

sub dump_hopper {
	my $self = shift;
	my @flush = ();
	$self->hopper(@flush);
}

sub euclidean_distance {
	my $self = shift;
	my %params = validate(
		@_, {
			x2 => 1,
			y2 => 1,
		}
	);
	my $distance = sqrt(($self->x - $params{'x2'})^2 + ($self->y - $params{'y2'})^2);
	
	return $distance;
}
__PACKAGE__->meta->make_immutable;
1;