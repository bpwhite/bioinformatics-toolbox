# Installation script for the Bioinformatics-toolbox

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

print "Detecting operating system...$^O\n";

if($^O =~ m/MSWin32/) {
	print "Found ".$^O."\n";
	print "Detecting distribution...\n";
	my @perl_version = `perl -v`;
	my $perl_distribution = '';
	foreach my $line (@perl_version) {
		if ($line =~ m/ActiveState/) {
			$perl_distribution = 'ActiveState';
			last;
		}
		if ($line =~ m/Strawberry/) {
			$perl_distribution = 'Strawberry';
			last;
		}
	}
	if ($perl_distribution eq 'ActiveState') {
		print "Installing BioPerl...\n";
		print `ppm install BioPerl`;
		print "Done.\n";
		print "Installing MinGW...\n";
		print `ppm install mingw`;
		print "Done.\n";
		print "Installing other modules...\n";
		print `ppm install Moose`;
		print `ppm install MooseX::ClassAttribute`;
		print `ppm install Params::Validate`;
		print `ppm install FindBin`;
		print `ppm install Statistics::Descriptive`;
		print `ppm install Benchmark`;
		print `ppm install POSIX`;
		print `ppm install File::Copy`;
		print "Done.\n";
	} elsif ($perl_distribution eq 'Strawberry') {
		print "Your version of Perl is not currently supported by this installer.\n";
		exit;
	} else {
		print "Your version of Perl is not currently supported by this installer.\n";
		exit;
	}
} elsif ($^O =~ m/linux/) {
	print "Found ".$^O."\n";
	print "Your operating system is not supported by this installer.\n";
	exit;
	print `apt-get install BioPerl`;
} else {
	print "Found ".$^O."\n";
	print "Your operating system is not supported by this installer.\n";
	exit;
}

