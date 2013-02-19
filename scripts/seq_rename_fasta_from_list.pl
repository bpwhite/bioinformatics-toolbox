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
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

# my $input_file = '18s_mega_aln.fas';
chomp(my $input_file = <>);
fix_fasta($input_file);

my $seqio  = Bio::SeqIO->new(-file => $input_file , '-format' => 'Fasta');
my @alignment = ();
while((my $seq = $seqio->next_seq())) {
	push(@alignment,$seq);
}

chomp(my $rename_list = <>);
open F, "< $rename_list";
my @rename_list = <F>;
close F;

my $renamed = clean_file_name($input_file)."_renamed_aln.fas";

unlink $renamed;
open (RENAMED, '>>'.$renamed);

my $rename_length = scalar(@rename_list);
my $align_i = 0;
foreach my $source_name (@rename_list) {
	$source_name =~ s/ /\_/g; 
	print RENAMED ">".$source_name;
	print RENAMED $alignment[$align_i]->seq."\n";
	$align_i++;
}

close(RENAMED);

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


sub clean_file_name {
	# Takes a filename and strips its extension off, and returns the file name
	my ($file_name) = @_;
	my @split_file_name = split(/\./,$file_name);
	my $new_file_name = $split_file_name[0];

	return $new_file_name;
}