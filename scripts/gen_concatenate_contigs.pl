# Concatenates larges FASTA files (contigs) into one single "Gene".
#
# Copyright (c) 2011, Bryan White

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

use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

# chomp(my $input_file = <>);
# my $input_file = "test.fasta";
my @input_files = <*.fasta>;

print "Overwrite?\n";
chomp(my $overwrite_output = <>);

print "Fix fasta (remove blanks - do if first run) [y]es or [n]\n";
# chomp(my $fix_fasta = <>);
my $fix_fasta = "yes";
	
foreach my $input_file (@input_files) {
	# Fix the fasta
	# If first run

	if($fix_fasta eq "yes" || $fix_fasta eq "y") {
		fix_fasta($input_file);
	}
	# Load the genome file
	my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');
	
	# Create the sequence name
	my @output_id_parsed = split(".fasta",$input_file) or die("Input files not .fasta files.");
	my $output_id = $output_id_parsed[0];
	
	# Create the output file
	my $output_file = "appended_".$output_id_parsed[0].".fasta";

	open (OUTPUT, '>>'.$output_file);
	
	my $seq_i = 1;
	while((my $seqobj = $seqio->next_seq())) {
		if($seq_i == 1) {
			print OUTPUT ">".$output_id."\n";
			print OUTPUT $seqobj->subseq(1,$seqobj->length)."\n";
		} else {
			print OUTPUT $seqobj->subseq(1,$seqobj->length)."\n";
		}
		$seq_i++;
	}
	close(OUTPUT);
	
	$seqio = '';
	
	if($overwrite_output eq "y" || $overwrite_output eq "yes") {
		print "Trying to delete...\n";
		rename($output_file,$input_file);
		if(unlink($output_file) == 1) { print "File deleted successfully."};
	}
}

print scalar(@input_files)." sequences output.";
sub fix_fasta {
	my ( $f ) = shift;
	open F, "< $f" or die "Can't open $f : $!";
	my @fasta = <F>;
	close F;
	
	for(my $i = 0; $i < 10; $i++) {
		foreach my $line (@fasta) {
			if($line =~ /^>/) {
				$line =~ s/ /_/; # replace whitespace with _
			}
		}
	}

	unlink $f;
	open (MYFILE, '>>'.$f);
	foreach my $line(@fasta) {
		print MYFILE $line;
	}
	close(MYFILE);
}