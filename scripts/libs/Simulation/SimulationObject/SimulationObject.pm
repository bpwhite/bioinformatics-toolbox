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

package SimulationObject;
use strict;
use warnings;

use Moose;
use MooseX::ClassAttribute;
use Params::Validate qw(:all);


class_has 'logger_prefix' => (
	is      => 'rw',
	isa     => 'Any',
	default => 'default_log',
);

sub log_message {
	my $self = shift;
}

sub print_to_logfile {
	my $self = shift;
	my $message = shift;
	if(-e $self->logger_prefix.'.txt') {
		# Append message
		open FILE, '>>', $self->logger_prefix.'.txt' or die $!;
		print FILE $message;
	} else {
		# Open message
		open FILE, '>', $self->logger_prefix.'.txt' or die $!;
		print FILE $message;
	}
	close FILE;
}

sub clear_log_file {
	my $self = shift;
	open FILE, '>>', $self->logger_prefix.'.txt' or die $!;
	truncate (FILE, 0);
}

__PACKAGE__->meta->make_immutable;
1;