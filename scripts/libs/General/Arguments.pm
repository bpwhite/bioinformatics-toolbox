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
has 'arguments_v' => (
# Takes ARGV input as an arrayf
	is	=> 'ro',
	isa => 'ArrayRef[Any]',
	builder => '_build_argument_v',
	required => 1,
	# auto_deref => 1,
);
sub _build_argument_v {
	return [];
}

has 'option_defs' => (
# Stores the default parameters for each option
  is        => 'ro',
  isa       => 'HashRef',
  builder   => '_build_option_defs',
  required 	=> 1,
);
sub _build_option_defs {
	return {};
}

has 'options' => (
# Actual list of options
  is        => 'rw',
  isa       => 'HashRef',
  builder   => '_build_options',
);
sub _build_options {
	return {};
}

########################################################################
sub BUILD {
	my $self = shift;
	# print scalar @{$self->arguments_v}."\n";
	# foreach my $test (@{$self->arguments_v}) {
		# print $test."\n";
	# }
	# exit;
	# Check for proper number of parameters
	die "Odd number of parameters.\n" if scalar @{$self->arguments_v} %2;
	
	
	# Loop through arguments
	eval { 
		for (my $arg_i = 0; $arg_i <= (scalar(@{$self->arguments_v})/2); $arg_i+=2) {
			if ($self->arguments_v->[$arg_i] eq '--help') { print "Seq_convert_genbank help info\n"; exit; };
			my $option 	= $self->arguments_v->[$arg_i];
			my $value 	= $self->arguments_v->[$arg_i+1];
			$self->options->{$option} = $value;
		}
	};
	if ($@) {
		my @parameter_error_split = split(/ /,$@);
		my @parameter_error = grep $_ =~ /'/, @parameter_error_split;
		print "*************************\n";
		print "Fatal parameter error: $parameter_error[0].\n";
		print "*************************\n\n\n";
		# print $@."\n";
		exit;
	}
}
########################################################################

1;