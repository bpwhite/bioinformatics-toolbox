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
	
	# Download BioPerl from bioperl.org DIST
	print "Downloading BioPerl...\n";
	my $url = 'http://bioperl.org/DIST/BioPerl-1.6.1.zip';
	my $file = 'BioPerl-1.6.1.zip';
	if(-e $file) {
		print "Already downloaded $file\n";
	} else {
		getstore($url,$file);
	}
	
	# Get the clean file name, i.e. no file extension.
	my @delimited_file_name = split(/\./,$file);
	my $clean_file_name = $delimited_file_name[0].".".$delimited_file_name[1].".".$delimited_file_name[2];
	if(-e "$FindBin::Bin/tmp/$clean_file_name") {
		goto skip_extract;
	}

	# Extract BioPerl zip.
	print "Extracting BioPerl to $FindBin::Bin/tmp/$clean_file_name...\n";
	my $ae = Archive::Extract->new( archive => "$FindBin::Bin/$file" );
	my $ok = $ae->extract(to => "$FindBin::Bin/tmp") or die $ae->error;

	skip_extract:
	
	# Install CPAN modules.
	print "Installing CPAN...\n";
	system("perl -MCPAN -e install CPAN");
	system("perl -MCPAN -e o conf prefer_installer MB; o conf prerequisites_policy follow; o conf commit; ");
	system("perl -MCPAN -e o conf prerequisites_policy follow");
	system("perl -MCPAN -e o conf commit");
	system("perl -MCPAN -e install Module::Build");
	system("perl -MCPAN -e install Test::Harness");
	# Change directory to temp working directory to build BioPerl.
	chdir("$FindBin::Bin/tmp/$clean_file_name");
	# exit;
	print "Building BioPerl...\n";
	system("perl Build.PL");
	# exit;
	print "Testing BioPerl...\n";
	system("perl Build test");
	print "Installing BioPerl...\n";
	system("perl Build install");
	
	print "Done!";
	exit;
	# print "Detecting distribution...\n";
	# my @perl_version = `perl -v`;
	# my $perl_distribution = '';
	# foreach my $line (@perl_version) {
		# if ($line =~ m/ActiveState/) {
			# $perl_distribution = 'ActiveState';
			# last;
		# }
		# if ($line =~ m/Strawberry/) {
			# $perl_distribution = 'Strawberry';
			# last;
		# }
	# }
	# if ($perl_distribution eq 'ActiveState') {
		# print `ppm install PPM-Repositories`;
		# print "Installing BioPerl...\n";
		# print `ppm repo add http://bioperl.org/DIST`;
		# print `ppm repo add uwinnipeg`;
		# print `ppm repo add trouchelle http://trouchelle.com/ppm12/`;
		# print `ppm install BioPerl`;
		# exit;
		# print "Done.\n";
		# print "Installing MinGW...\n";
		# print `ppm install mingw`;
		# print "Done.\n";
		# print "Installing other modules...\n";
		# print `ppm install Moose`;
		# print `ppm install MooseX::ClassAttribute`;
		# print `ppm install Params::Validate`;
		# print `ppm install FindBin`;
		# print `ppm install Statistics::Descriptive`;
		# print `ppm install Benchmark`;
		# print `ppm install POSIX`;
		# print `ppm install File::Copy`;
		# print "Done.\n";
	# } elsif ($perl_distribution eq 'Strawberry') {
		# print "Your version of Perl is not currently supported by this installer.\n";
		# exit;
	# } else {
		# print "Your version of Perl is not currently supported by this installer.\n";
		# exit;
	# }
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

