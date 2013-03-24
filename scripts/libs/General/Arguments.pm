# Interface for command line arguments
#
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

package General::Arguments;

use Moose;

########################################################################
# Class Variables

########################################################################
# Attributes
has 'argument_v' => (
	is	=> 'ro',
	isa => 'ArrayRef',
	required => 1,
	builder => '_build_argument_v',
);
sub _build_argument_v {
	return [];
}
has 'option' => (
  is        => 'rw',
  isa       => 'HashRef',
  default   => sub { {} },
);


########################################################################
sub BUILD {
	my $self = shift;
	
	$self->option->{'test'} = 'blah';
}
########################################################################

1;