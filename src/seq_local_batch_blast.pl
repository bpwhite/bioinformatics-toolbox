#!/usr/bin/env perl
# Assign taxon names to a list of sequences from a local blast database

# Copyright (c) 2013, 2014 Bryan White, bpcwhite@gmail.com

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
use lib "$FindBin::Bin/libs/Sequence"; 
use lib "$FindBin::Bin/libs/";

use strict;
use warnings;
# Import sequence libs
use Sequence::Fasta;
use Sequence::Kimura_Distance;
use General::Arguments;
use Sequence::Bootstrap;
use Sequence::Garli;

# BioPerl libs
use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::SimpleAlign;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;

# Non-bio modules
use Benchmark qw(:all);
use Getopt::Long;
use String::Random;
##################################################################
# Start benchmark
my $t0 = Benchmark->new;
my $k2p1 = 0;
##################################################################

my $alignment 	= '';
my $output 		= '';
my $homology_level = 80;
my $alignment_length_pcnt = 0.60;
my $blastdb_name = '';

GetOptions ("aln=s" 			=> \$alignment,
			"out=s"				=> \$output,
			"db=s"				=> \$blastdb_name,
			"homology=s"		=> \$homology_level,
			"aln_length=s"		=> \$alignment_length_pcnt)
or die("Error in command line arguments\n");


##################################################################
# Import alignment
print "Importing alignment file ".$alignment."...\n";
my $alignin = Bio::AlignIO->new(-format => 'fasta',
								-file   => $alignment );
my $original_aln = $alignin->next_aln;
##################################################################

my @starting_sequence_array = ();
foreach my $seq ($original_aln->each_seq) {
	push(@starting_sequence_array,$seq);
}

# Blast output key
# query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
my $seq_pass_i = 0;
my $seq_fail_i = 0;
open (OUT, '>', $output);
open (OUTLOG, '>', $output."_log.txt");
foreach my $seq (@starting_sequence_array) {
	my $seq_id = $seq->id;
	$seq->id($seq_id);
	my $seq_string = $seq->seq;
	
	# Remove primer sequences
	# LCO F GGTCAACAAATCATAAAGATATTGG
	# HCO R TAAACTTCAGGGTGACCAAAAAAT
	# mLepF1 GCTTTCCCACGAATAAATAATA
	# mLepR1 TAAACTTCTGGATGTCCAAAAAATCA

	$seq_string =~ s/GGTCAACAAATCATAAAGATATTGG//g;
	$seq_string =~ s/TAAACTTCAGGGTGACCAAAAAAT//g;
	$seq_string =~ s/GGTCAACAAATCAATAAAGATATTGG//g;

	my $query_length = fast_seq_length($seq_string);
	
    my $temp_query .= 'temp_query_'.int(rand(99999)).'.fas';
	
	open (QUERY, '>'.$temp_query);
	print QUERY ">".$seq_id."\n";
	print QUERY $seq_string."\n";
	close (QUERY);
	
	my $blast_output = `blastn -query $temp_query -db $blastdb_name -outfmt 6 -max_target_seqs 5 -num_threads 4`;
	unlink $temp_query;
	
	my @blast_lines = split(/\n/,$blast_output);
	my @split_blast = ();
	foreach my $blast_line (@blast_lines) {
		@split_blast = split(/\t/,$blast_line);
		if ($blast_output eq '') {
			$seq_fail_i++;
			next;
		}
		$blast_output =~ s/\n//g;
		# BLAST % identity
		if($split_blast[2] < $homology_level) {
			$seq_fail_i++;
			print OUTLOG $blast_output.",Fail homology level\n";
			# print $blast_output."\n";
			# print "Fail homology level\n";
			next;
		}
		# BLAST alignment length
		my $pcnt_query_match = $split_blast[3]/$query_length;
		if($pcnt_query_match < $alignment_length_pcnt) {
			$seq_fail_i++;
			print OUTLOG $blast_output.",Fail alignment length\n";
			# print $blast_output."\n";
			# print "Fail alignment length at $pcnt_query_match\n";
			next;
		}
		$seq_pass_i++;
		last;
	}
	next if $seq_pass_i == 0;
	next if (!defined($split_blast[0]));
	
	$blast_output =~ s/\n//g;
	print $blast_output."\n";
	
	my @seq_site_code = split(/\|/,$seq_id);
	my $site_code = $seq_site_code[-1];
	
	my @seq_read_id = split(/_/,$seq_id);
	my $read_id = $seq_read_id[0];
	
	my $new_seq_id = $read_id."|".$split_blast[1]."|".$split_blast[2]."|".$site_code;

	print OUT ">".$new_seq_id."\n";
	print OUT $seq_string."\n";	

	
	# print OUT ">".$seq->id."\n";
	# print OUT $seq->seq."\n";
}
close (OUT);

print "Failed: ".$seq_fail_i."\n";
print "Passed: ".$seq_pass_i."\n";

my $t1 = Benchmark->new;
my $time_diff = timediff($t1, $t0);
print "\n";
print timestr($time_diff)."\n";

sub fast_seq_length {
	my $seq = shift;
	
	$seq =~ s/-/ /g;
	$seq =~ s/\s+//g;
	
	return length($seq);
}