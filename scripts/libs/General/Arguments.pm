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
use Moose::Meta::Attribute::Native::Trait::Array;
use Moose::Meta::Attribute::Native::Trait::Hash;

use Hash::Util qw(
                     hash_seed all_keys
                     lock_keys unlock_keys
                     lock_value unlock_value
                     lock_hash unlock_hash
                     lock_keys_plus hash_locked
                     hidden_keys legal_keys
                   );
				   
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
has 'options' => (
  is        => 'rw',
  isa       => 'HashRef',
  default   => sub { {} },
);


########################################################################
sub BUILD {
	my $self = shift;
	
	# my $args = $self->argument_v;

	# Check for proper number of parameters
	die "Odd number of parameters.\n" if scalar $self->argument_v %2;

	# Loop through arguments
	eval { 
		foreach my $argument ($self->argument_v) {
			if ($argument eq '--help') { print "Seq_convert_genbank help info\n" } ;
		}
	};
	if ($@) {
		my @parameter_error_split = split(/ /,$@);
		my @parameter_error = grep $_ =~ /'/, @parameter_error_split;
		print "*************************\n";
		print "Incorrect parameter used: $parameter_error[0].\n";
		print "*************************\n\n\n";
		# print $@."\n";
		die;
	}
}
########################################################################

1;