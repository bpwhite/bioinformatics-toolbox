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
my $dist_cutoff 	= '0.55';
my $bases_cutoff 	= '0.70';
my $otu_cutoff = 0;
my $threads = 1;
GetOptions ("gb=s" 				=> \$genbank_csv,
			"out=s"				=> \$output_tag,
			"match=s"			=> \$match_aln_file,
			"cutoff=s"			=> \$dist_cutoff,
			"bases=s"			=> \$bases_cutoff,
			"otu-cutoff=s"		=> \$otu_cutoff,
			"threads=s"			=> \$threads)
or die("Error in command line arguments\n");

# Load the alignment to match query sequences against
my $match_aln = Bio::AlignIO->new(-format => 'fasta',
								-file   => $match_aln_file );
my $match_aln_obj = $match_aln->next_aln;
		
my @match_aln_array = ();
my @match_seq_ids = ();
foreach my $match_seq ($match_aln_obj->each_seq) {
	push(@match_seq_ids, $match_seq->id);
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

# Add match sequence to the output FASTA as first sequence
foreach my $match_seq (@match_aln_array) {
	print QUERY_FAS ">".$match_seq->id."\n";
	print QUERY_FAS $match_seq->seq."\n";
}
# Loop through Genbank CSV lines and output to a FASTA file
my %params_hash = ();
my $line_counter = 0;
my $column_headers = '';
foreach my $line (@genbank_lines) {
	#print $line;
	if($line_counter == 0) {
		$column_headers = $line;
	}
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
		print QUERY_FAS ">".$accession."\n";
		print QUERY_FAS $nuc_seq."\n";		
	}

	$line_counter++;
}

close(QUERY_FAS);

my $alignment_depth = 0;

my $final_depth = recursive_alignment();

print "Final file = ".$final_depth."\n";

my $final_alignment_file = $output_tag."_".$final_depth."_aln.fas";
my $final_alignment_tag = $output_tag."_".$final_depth."_aln";

# Reload finalized alignment for parsing into output file.
print "Loading... ".$final_alignment_file."\n";
my $final_alignment = Bio::AlignIO->new(-format => 'fasta',
								-file   => $final_alignment_file );
my $final_alignment_obj = $final_alignment->next_aln;
my @final_alignment_array = ();

foreach my $final_seq ($final_alignment_obj->each_seq) {
	my $seq_id = $final_seq->id;
	$seq_id =~ s/^_R_//;
	#print $seq_id."\n";
	$final_seq->id($seq_id);
	push(@final_alignment_array,$final_seq);
}

# Assign OTU tags to qc checked sequences.
my %otu_hash = ();
if($otu_cutoff != 0) {
	my $otu_output = "./dnab_otu_delim.pl -aln1 $final_alignment_file -cutoff $otu_cutoff -skip-intra-dist 1 -skip-nn 1";
	print $otu_output."\n";
	my $otu_results = `$otu_output`;
	my $otu_results_path = "./".$final_alignment_tag."_output/".$final_alignment_tag."_otu_results.csv";
	print $otu_results_path."\n";

	open(OTU_RESULTS, '<'.$otu_results_path);
	my @otu_lines = <OTU_RESULTS>;
	close(OTU_RESULTS);
	foreach my $otu_line (@otu_lines) {
		$otu_line =~ s/\n//g;
		my @split_otu_line = split(/\,/,$otu_line);
		#foreach my $split_otu (@split_otu_line) {
		#	print $split_otu." => ";
		#}
		#print "\n";

		$otu_hash{$split_otu_line[1]} = $split_otu_line[0];
	}
}

# Delete all the recursively generated alignment files
for(my $aln_i = 0; $aln_i <= $final_depth; $aln_i++) {
	unlink( $output_tag."_".$aln_i."_aln.fas")
}

print "Matching aligned sequence to genbank data\n";
my $new_genbank_file = $output_tag."_qchecked.csv";
unlink($new_genbank_file);

open(QCDONE, '>>'.$new_genbank_file);
$column_headers =~ s/\n//g;
print QCDONE $column_headers."otu_id,aligned_fasta\n";

my $printed_passed_seqs = 0;
my $genbank_line_i = 0;
foreach my $genbank_line (@genbank_lines) {
	if($genbank_line_i == 0) {
		$genbank_line_i++;	
		next;
	}
	#print $genbank_line."\n";
	my $found_match = 0;
	$genbank_line =~ s/\n//g;
	foreach my $final_seq (@final_alignment_array) {
		#print $final_seq->id."\n";
		my $final_seq_id = $final_seq->id;
		my $otu_id = $otu_hash{$final_seq_id};
		#print $final_seq_id."\n";
		next if ($final_seq_id ~~ @match_seq_ids);

		if($genbank_line =~ /$final_seq_id/) {
			#print $final_seq_id." => ".$genbank_line."\n";
			#print $final_seq_id."\n";
			$printed_passed_seqs++;
			print QCDONE "\"".$genbank_line."\",\"".$otu_id."\",\"\$".$final_seq->seq."\"\n";
			$found_match = 1;
			last;
		}
	}
	if($found_match == 0) {
		print QCDONE "\"".$genbank_line."\",failedqc\n";
	}
	$genbank_line_i++;
}

close(QCDONE);

print "Printed seqs: $printed_passed_seqs\n";
print "Done!\n";

sub recursive_alignment {
	my $aligned = $output_tag."_".$alignment_depth."_aln.fas";
	#unlink($aligned);

	my $last_depth = $alignment_depth-1;
	my $last_alignment = $output_tag."_".$last_depth."_aln.fas";

	my $mafft_string = '';
	if($alignment_depth == 0) {
		unlink($aligned);
		$mafft_string = "mafft --auto --preservecase --adjustdirection --preservecase --thread $threads --quiet --maxiterate 0 --retree 1 --6merpair $fasta_output > $aligned";
		print "Calling $mafft_string at depth $alignment_depth\n";
		my $mafft_output = `$mafft_string`;
		#system($mafft_string);
	} else {
		unlink($aligned);
		$mafft_string = "mafft --auto --preservecase --adjustdirection --preservecase --thread $threads --quiet --maxiterate 0 --retree 1 --6merpair $last_alignment > $aligned";
		print "Calling $mafft_string at depth $alignment_depth\n";
		my $mafft_output = `$mafft_string`;
		#system($mafft_string);
	}

	# Reimport aligned sequences and put in array
	print "Loading... ".$aligned."\n";
	my $qc_alignment = Bio::AlignIO->new(-format => 'fasta',
									-file   => $aligned );
	my $qc_alignment_obj = $qc_alignment->next_aln;
	my @qc_aln_array = ();
	foreach my $qc_seq ($qc_alignment_obj->each_seq) {
		#print $qc_seq->id."\n";
		push(@qc_aln_array,$qc_seq);
	}
	unlink($aligned);
	print "Printing to ".$aligned."\n";
	open(ALIGNED, '>>'.$aligned);
	# Determine alignment length, including gaps
	my $aln_length = length($qc_aln_array[0]->seq);
	my $num_bases_cutoff = $aln_length*$bases_cutoff;
	# Loop through aligned sequences and check distance against the match query
	my $match_seq_original = $qc_aln_array[0];
	my $max_k2p = 0;
	my $seqs_passed = 0;
	foreach my $qc_seq (@qc_aln_array) {
		my $seq1 = $qc_seq->seq;
		my $seq2 = $match_seq_original->seq;

		my ($k2p_distance, $transitions,$transversions,$bases_compared) = k2p_unpack($seq1,$seq2,$aln_length);
		#print $k2p_distance."\n";

		if($k2p_distance > $max_k2p) {
			$max_k2p = $k2p_distance;
			#print $max_k2p."\n";
		}
		#if (($k2p_distance < $dist_cutoff) && ($bases_compared > $num_bases_cutoff)) {
		if ($k2p_distance < $dist_cutoff) {
			$seqs_passed++;
			#print $k2p_distance. " => ".$qc_seq->id."\n";
			print ALIGNED ">".$qc_seq->id."\n";
			print ALIGNED $qc_seq->seq."\n";
		} elsif ($k2p_distance > $dist_cutoff) {
			#print $k2p_distance." with ".$bases_compared." compared at ".$qc_seq->id."\n";
		}
	}
	close(ALIGNED);
	print "Seqs passed: $seqs_passed\n";

	if($max_k2p > $dist_cutoff) {
		print $max_k2p."\n";
		$alignment_depth++;
		recursive_alignment();
	}
	return $alignment_depth;
}

