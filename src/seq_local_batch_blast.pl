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
use Statistics::Descriptive;
use Benchmark qw(:all);
use POSIX qw(ceil);
use File::Copy;
use Math::Random::MT::Perl qw(srand rand irand);
use Digest::SHA qw(sha1 sha1_base64);
use Term::Spinner;
use Data::Dumper;
use Storable qw(dclone);
use File::Copy;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
##################################################################
# Start benchmark
my $t0 = Benchmark->new;
my $k2p1 = 0;
##################################################################

my $alignment 	= '';
my $output 		= '';

GetOptions ("aln=s" 			=> \$alignment,
			"out=s"				=> \$output,
			"db=s"				=> \$blastdb_name,)
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
my $temp_query = 'temp_query.fas';
open (OUT, '>', $output);
foreach my $seq (@starting_sequence_array) {
	my $seq_id = $seq->id;
	$seq->id($seq_id);
	open (QUERY, '>'.$temp_query);
	print QUERY ">".$seq_id."\n";
	print QUERY $seq->seq."\n";
	close (QUERY);
	
	my $blast_output = `blastn -query $temp_query -db $blastdb_name -outfmt 6 -max_target_seqs 1 -num_threads 2`;
	print $blast_output."\n";
	if ($blast_output eq '') {
		$seq_fail_i++;
		next;
	}
	$seq_pass_i++;
	
	my @seq_site_code = split(/\|/,$seq_id);
	my $site_code = $seq_site_code[-1];
	
	
	my @split_blast = split(/\t/,$blast_output);
	my $new_seq_id = $split_blast[1]."|".$site_code;
	
	print OUT ">".$new_seq_id."\n";
	print OUT $seq->seq."\n";	

	unlink $temp_query;
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