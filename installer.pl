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
use Archive::Extract;
use LWP::Simple;
use FindBin;
# use lib "$FindBin::Bin/tmp";
print "Welcome to the Bionformatics ToolBox (BTBox) Installer\n\n";
print "Currently supported operating systems: MSWin32\n";
print "Detecting operating system...$^O\n";

if($^O =~ m/MSWin32/) {
	print "Found ".$^O."\n";
	
	# Install CPAN modules.
	print "Installing CPAN...\n";
	system("perl -MCPAN -e install CPAN");
	system("perl -MCPAN -e o conf prefer_installer MB; o conf prerequisites_policy follow; o conf commit; ");
	system("perl -MCPAN -e o conf prerequisites_policy follow");
	system("perl -MCPAN -e o conf commit");
	system("perl -MCPAN -e install Module::Build");
	system("perl -MCPAN -e install Test::Harness");
	system("perl -MCPAN -e install Inline");
	system("perl -MCPAN -e install Inline::C");
	system("perl -MCPAN -e install Inline::Files");
	system("perl -MCPAN -e install Statistics::Descriptive");
	print "Done!";
	exit;

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

