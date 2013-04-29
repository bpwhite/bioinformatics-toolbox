#!/usr/bin/perl
# Merges alignments that have the same order.
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
use lib "$FindBin::Bin/libs/";
use General::Arguments;
use Sequence::Fasta;

# BioPerl libs
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::Tree::Node;
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;

# my $params = General::Arguments->new(	arguments_v => \@ARGV,
										# option_defs => {'-outp' 		=> 'merged', 				# Folder where CSV files are stored.
													# }
													# );

# my $output_file = $params->options->{'-outp'};
# my $output_file .= '.fas';

my $output_file = 'merged.fas';
my @alignment_files = <*.fas>;
my @files_to_merge = ();
my $num_alignments = scalar @alignment_files;
while (my $done_selecting == 0) {
	print "Select an alignment file to merge\n";
	for (my $aln_i = 0; $aln_i < $num_alignments; $aln_i++) {
		print "[".($aln_i+1)."] $alignment_files[$aln_i]\n";
	}
	print "Select [0] to finalize selection\n";
	chomp (my $aln_to_add = <>);
	if($aln_to_add == 0) {
		last;
	} else {
		my @split_line = split(",", $aln_to_add);
		foreach my $split (@split_line) {
			print "Adding $alignment_files[$split-1]\n";
			if ($alignment_files[$split-1] ~~ @files_to_merge) {
				print "Already added alignment.\n";
				next;
			}
			push(@files_to_merge, $alignment_files[$split-1]);
			if($aln_to_add == 0) {
				last;
			}
		}
	}
	print "\n";
}

my @merged_output = ();
my $merged_counter = 0;
my $num_files_to_merge = scalar @files_to_merge;
foreach my $file2merge (@files_to_merge) {
	print "Merging $file2merge.\n";
	if($merged_counter == 0) {
		fix_fasta($file2merge);
		my $alignin = Bio::AlignIO->new(-format => 'fasta',
									-file   => $file2merge );
		my $original_aln = $alignin->next_aln;
		foreach my $seq ($original_aln->each_seq) {
			my $merge_line = '>'.$seq->id."\n".$seq->seq;
			push(@merged_output, $merge_line);
		}
	} else {
		fix_fasta($file2merge);
		my $alignin = Bio::AlignIO->new(-format => 'fasta',
									-file   => $file2merge );
		my $original_aln = $alignin->next_aln;
		my $seq_counter = 0;
		foreach my $seq ($original_aln->each_seq) {
			my $merge_line = $seq->seq;
			$merged_output[$seq_counter] .= $merge_line;
			$seq_counter++;
		}
	}
	$merged_counter++;
}

unlink ($output_file);
open (MERGED, '>>'.$output_file) or die "Could not open $f";
foreach my $merged_line (@merged_output) {
	print MERGED $merged_line."\n";
}
close MERGED;
