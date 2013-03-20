# Converts a nucleotide alignment to a protein alignment, and then 
# allows that protein alignment to be back translated into a nucleotide
# alignment.
# This will be useful for protein coding genes with gapped alignments.
# E.g. heat shock proteins.
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
use FindBin;
use lib "$FindBin::Bin/libs"; 

use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;

use Data::Dumper;

my $nucleotide_seqs_file = 'hsp70_genomic_CDS.fas';
my $protein_aln_file = 'hsp70_CDS_protein_aln.fas';
my $output_file = 'nucleotide_alignment.fas';

my @nucleotide_seqs = ();
my @protein_translation = ();

# Ungapped nucleotide sequences
my $nucleotide_aln_io = Bio::AlignIO->new(-format => 'fasta',
								-file   => $nucleotide_seqs_file );
my $nucleotide_aln = $nucleotide_aln_io->next_aln;
foreach my $seq ($nucleotide_aln->each_seq) {
	push(@nucleotide_seqs, $seq->seq);
}


# Gapped protein alignment
my $protein_aln_io = Bio::AlignIO->new(-format => 'fasta',
								-file   => $protein_aln_file );
my $protein_aln = $protein_aln_io->next_aln;
my $seq_i = 0;

unlink $output_file;
open (OUTPUT, '>>'.$output_file);
foreach my $seq ($protein_aln->each_seq) {
	print OUTPUT ">".$seq->id."\n";
	my $current_nucleotide_seq = $nucleotide_seqs[$seq_i];
	my @nuc_seq_string = split(//,$current_nucleotide_seq);
	my @prot_seq_string = split(//,$seq->seq);
	my $codon_count = 0;
	for (my $i = 0; $i <= length($seq->seq)-2; $i++) {
		my $current_letter = $prot_seq_string[$i];
		if($current_letter ne '-') {
			# print $codon_count." => ";
			print OUTPUT $nuc_seq_string[$codon_count*3].$nuc_seq_string[$codon_count*3+1].$nuc_seq_string[$codon_count*3+2];
			# print $prot_seq_string[$i]."\n";
			$codon_count++;
		} else {
			print OUTPUT "---";
		}
	}
	print OUTPUT "\n";
	# print $current_nucleotide_seq."\n";
	$seq_i++;
}

close(OUTPUT);