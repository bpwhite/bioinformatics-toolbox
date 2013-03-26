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

	# Loop through arguments and set parameters
	eval {
		die "Missing parameters.\n" if scalar @{$self->arguments_v} == 0;

		my @help = grep { $_ eq '--help'} $self->arguments_v;
		if (scalar(@help) > 0) { 
			print $self->print_help; 
			exit;
		}
		die "Odd number of parameters.\n" if scalar @{$self->arguments_v} %2;
		# Set ARGV parameters
		for (my $arg_i = 0; $arg_i < (scalar(@{$self->arguments_v})); $arg_i+=2) {
			my $option = $self->arguments_v->[$arg_i];			
			if(exists($self->option_defs->{$option})) {
				my $value 					= $self->arguments_v->[$arg_i+1];
				$self->options->{$option} 	= $value;
				print "\t".$option." => ".$value."\n";
			} else {
				die "Parameter does not exist: $option";
			}
		}
		# Set default parameters for params not set by ARGV
		while ( my ($key, $value) = each(%{$self->option_defs}) ) {
			if (!exists ($self->options->{$key})) {
				$self->options->{$key} = $self->option_defs->{$key};
				print "\t$key => $value (default)\n";
			}
		}
	};
	if ($@) {
		my @parameter_error_split = split(/ /,$@);
		my @parameter_error = grep $_ =~ /'/, @parameter_error_split;
		print "*************************\n";
		print "Fatal parameter error:\n\n";
		print "$@";
		print "*************************\n";
		exit;
	}
}
########################################################################

sub print_options {
	my $self = shift;
	print "Using parameters...\n";
	while ( my ($key, $value) = each(%{$self->options}) ) {
        print "\t$key => $value\n";
    }
}

sub print_help {
	# my $self = shift;
	
	# print "seq_convert_genbank.pl help\n";
	# while ( my ($key, $value) = each(%{$self->option_defs}) ) {
        # print "\t$key => $value\n";
    # }
}
1;