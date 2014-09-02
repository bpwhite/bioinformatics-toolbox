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
use Sequence::BLAST;
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
print "
******************************************************************

Copyright (c) 2013,2014 Bryan White, bpcwhite\@gmail.com

GNU General Public License, Version 3, 29 June 2007

This program comes with ABSOLUTELY NO WARRANTY; for details type.
This is free software, and you are welcome to redistribute it
under certain conditions. 

A full copy of the GPL 3.0 license should be accompanied with any 
distribution of this software.

******************************************************************
\n\n
";

my $alignment 	= '';
my $output 		= '';
my $identity_level = 80;
my $aln_length_pcnt = 0.6;
my $blastdb_name = '';
my $qc_only = 0;
my $multi_blast = 0;
my $max_target_seqs = 10;
my $short_seq = 0;

GetOptions ("aln=s" 			=> \$alignment,
			"out=s"				=> \$output,
			"db=s"				=> \$blastdb_name,
			"homology=s"		=> \$identity_level,
			"aln_length=s"		=> \$aln_length_pcnt,
			"qc_only=s"			=> \$qc_only,
			"multi_blast=s"		=> \$multi_blast,
			"max_target_seqs=s"	=> \$max_target_seqs,
			"short_seq=s"		=> \$short_seq)
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
my %results_hash = ();
open (OUT, '>', $output);
open (OUTLOG, '>', $output."_log.txt");
close (OUTLOG);
foreach my $seq (@starting_sequence_array) {
	my $seq_id = $seq->id;
	
	my $seq_string = $seq->seq;
	my $query_length = fast_seq_length($seq_string);
	my $blast_output = '';
	if(exists($results_hash{$seq_string})) {
		$blast_output = $results_hash{$seq_string};
	} else {
		$blast_output = blast_output(seq_string 	=> $seq_string,
									seq_id 			=> $seq_id,
									blast_db		=> $blastdb_name,
									identity_level 	=> $identity_level,
									aln_length_pcnt => $aln_length_pcnt,
									log_file		=> $output."_log.txt",
									multi_blast		=> $multi_blast,
									max_target_seqs => $max_target_seqs,
									);
		$results_hash{$seq_string} = $blast_output;
	}
	print $blast_output."\n";
	$seq_fail_i++ if $blast_output eq 'NA';
	next if $blast_output eq 'NA';
	
	$seq_pass_i++;
	
	my @blast_results = split(/\n/,$blast_output);
	
	foreach my $blast_result (@blast_results) {
		$blast_result =~ s/ //g; # replace whitespace
		my @split_blast = split(/\t/,$blast_result);
		
		my @seq_site_code = split(/\|/,$seq_id);
		my $site_code = $seq_site_code[-1];
		
		my @seq_read_id = split(/_/,$seq_id);
		my $read_id = $seq_read_id[0];
		
		my $pcnt_query_match = sprintf( "%.2f",$split_blast[3]/$query_length);
		
		my $new_seq_id = '';
		if($qc_only == 0) {
			$new_seq_id = $read_id."|".$split_blast[1]."|".$split_blast[2]."|".$split_blast[10]."|".$split_blast[11]."|".$split_blast[12]."|".$seq_string."|".$query_length."|".$pcnt_query_match."|".$site_code;
		} elsif ($qc_only == 1) {
			$new_seq_id = $seq_id."|".$split_blast[2]."|".$query_length."|".$pcnt_query_match;
		}

		print OUT ">".$new_seq_id."\n";
		if($short_seq == 0) {
			print OUT $seq_string."\n";
		} elsif($short_seq == 1) {
			print OUT "a\n";
		}
		
	}
	
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
