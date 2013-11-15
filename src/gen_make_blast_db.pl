# Inputs a FASTA file and creates a blast database from it.
# 
# Copyright (c) 2013, Bryan White, bpcwhite@gmail.com

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
use strict;
use warnings;

# Create the output file
print "Output file name? ";
chomp(my $output_file = <>);
unlink $output_file;
open (MYFILE, '>>'.$output_file);

# Glob the all the files
my @fasta_files = ();
print "Primary contig folder? [includes all subfolders for 1 level]: ";
chomp(my $contig_folder1 = <>);
@fasta_files = <$contig_folder1\\*.fasta>;

print "Sub folder?: ";
chomp(my $contig_folder2 = <>);
if($contig_folder2 ne "") {
	@fasta_files = <$contig_folder1\\$contig_folder2\\*.fasta>;
}

# Cycle through the files and append them to the output file
foreach my $file (@fasta_files) {
	open F, "< $file" or die "Can't open $file : $!";
	my @fasta_lines = <F>;
	foreach my $line (@fasta_lines) {
		print MYFILE $line."\n";
	}
	close F;
	print "Done appending :".$file."\n";
}
close(MYFILE);
print "\n".scalar(@fasta_files)." files were concatenated into ".$output_file."\n";
print "Creating blast database...\n";

system("makeblastdb -in ".$output_file." -dbtype nucl");

print "Done!";
