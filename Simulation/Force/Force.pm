package Force;
use strict;
use warnings;

use Moose;
use SimulationObject;

extends 'SimulationObject';

has 'strength' => (
	is => 'rw',
	isa => 'Int',
	builder => '_seed_strength',
	required => 1);

has 'direction' => (
	is => 'rw',
	isa => 'Any',
	builder => '_seed_direction',
	required => 1);

sub _seed_strength {
	return int ( rand 3 );
}

sub _seed_direction {
	return int ( rand 3 );
}

__PACKAGE__->meta->make_immutable;

1;