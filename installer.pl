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

if($^O =~ m/Win/) {
	print "Installing BioPerl...\n";
	print `ppm install BioPerl`;
	print "Done.\n";
	print "Installing other modules...\n";
	print `ppm install Moose`;
	print `ppm install MooseX::ClassAttribute`;
	print `ppm install Params::Validate`;
	print `ppm install FindBin`;
	print `ppm install Statistics::Descriptive`;
	print `ppm install Benchmark`;
	print `ppm install POSIX`;
	print `ppm install  File::Copy`;
	print "Done.\n";
	
} elsif ($^O =~ m/linux/) {
	print `apt-get install BioPerl`;
}

