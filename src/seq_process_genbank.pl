#!/usr/bin/env perl
# Batch download taxa 
#
# Copyright (c) 2011, Bryan White, bpcwhite@gmail.com

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

use Sequence::QualityCheck;
use Sequence::Fasta;
use Sequence::Kimura_Distance;
use General::Arguments;

use Getopt::Long;

print "
******************************************************************

Copyright (c) 2013-2014 Bryan White, bpcwhite\@gmail.com

GNU General Public License, Version 3, 29 June 2007

This program comes with ABSOLUTELY NO WARRANTY; for details type.
This is free software, and you are welcome to redistribute it
under certain conditions. 

A full copy of the GPL 3.0 license should be accompanied with any 
distribution of this software.

******************************************************************
\n\n
";


my $genbank_csv		= '';
my $output_tag		= 'output';
my $match_aln_file  = '';
my $dist_cutoff 	= '0.70';
my $bases_cutoff 	= '0.70';
GetOptions ("gb=s" 				=> \$genbank_csv,
			"out=s"				=> \$output_tag,
			"match=s"			=> \$match_aln_file,
			"cutoff=s"			=> \$dist_cutoff,
			"bases=s"			=> \$bases_cutoff,)
or die("Error in command line arguments\n");

# Load the alignment to match query sequences against
my $match_aln = Bio::AlignIO->new(-format => 'fasta',
								-file   => $match_aln_file );
my $match_aln_obj = $match_aln->next_aln;
		
my @match_aln_array = ();
foreach my $match_seq ($match_aln_obj->each_seq) {
	push(@match_aln_array,$match_seq);
}

# Load the genbank CSV file
my @genbank_lines = ();
open (GB_LINES, '<'.$genbank_csv);
@genbank_lines = <GB_LINES>;
close(GB_LINES);

# Open the output query file for quality checking
my $fasta_output = $output_tag.".fas";
unlink $fasta_output;
open (QUERY_FAS, '>>'.$fasta_output);

# Loop through Genbank CSV lines and output to a FASTA file
my %params_hash = ();
my $line_counter = 0;
foreach my $line (@genbank_lines) {
	#print $line;
	$line =~ s/\n//g;
	$line =~ s/^\"|\"\,$//g;

	my @split_line = split('","',$line);

	if($line_counter == 0) {
		my $column_num = 0;
		foreach my $split (@split_line) {
			$params_hash{$split} = $column_num;
			#print $split." => ".$column_num."\n";
			$column_num++;
		}
	} else {
		# Build FASTA format sequences from each line
		my $column_num = 0;
		my $accession = '';
		my $nuc_seq = '';
		my $binom = '';
		foreach my $split (@split_line) {
			#print $column_num." => ".$split."\n";
			if($column_num == $params_hash{'nuc_seq'}) {
				$nuc_seq = $split;
			}
			if($column_num == $params_hash{'accession'}) {
				$accession = $split;
			}
			if($column_num == $params_hash{'binom'}) {
				#print $params_hash{'binom'}."\n";
				#print $column_num."\n";
				$binom = $split;
			}

			$column_num++;
		}
		print QUERY_FAS ">".$binom."|".$accession."\n";
		print QUERY_FAS $nuc_seq."\n";		
	}

	$line_counter++;
}

# Add match sequence to the output FASTA as first sequence
foreach my $match_seq (@match_aln_array) {
	print QUERY_FAS ">".$match_seq->id."\n";
	print QUERY_FAS $match_seq->seq."\n";
}

close(QUERY_FAS);

my $aligned = $output_tag."_aln.fas";
unlink($aligned);

# Align all query sequences.
my $mafft_string = "mafft --auto --preservecase --adjustdirection --preservecase $fasta_output > $aligned";
print "Calling $mafft_string\n";
system($mafft_string);

# Reimport aligned sequences and put in array
my $qc_alignment = Bio::AlignIO->new(-format => 'fasta',
								-file   => $aligned );
my $qc_alignment_obj = $qc_alignment->next_aln;
my @qc_aln_array = ();
foreach my $qc_seq ($qc_alignment_obj->each_seq) {
	print ">".$qc_seq->id."\n";
	print $qc_seq->seq."\n";
	push(@qc_aln_array,$qc_seq);
}

# Determine alignment length, including gaps
my $aln_length = length($qc_aln_array[0]->seq);

# Loop through aligned sequences and check distance against the match query
my $match_seq_original = $qc_aln_array[0];
my ($k2p_distance, $transitions,$transversions,$bases_compared) = 0;
foreach my $qc_seq (@qc_aln_array) {
	my $seq1 = $qc_seq->seq;
	my $seq2 = $match_seq_original->seq;
	($k2p_distance, $transitions,$transversions,$bases_compared) = k2p_unpack($seq1,$seq2,$aln_length);
	print $k2p_distance."\n";
}


#check_aln($split,$match_aln_file);

