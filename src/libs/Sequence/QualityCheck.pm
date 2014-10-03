#!/usr/bin/env perl
# This script takes an input of a FASTA sequence and outputs a 
# species delimitation regime given certain parameters.
# It is designed for use with the cytochrome oxidase I (COI) gene
# and to provide species delimitations useful for DNA barcoding.

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
use Sequence::Fasta;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw( check_aln	
						);
sub check_aln {
	my %params = ( 	aln => '',
					match => '',
					dist_cutoff => '0.50',
					bases => '0.70',
					@_
	);
	# Import sequence libs
	use Sequence::Fasta;
	use Sequence::Kimura_Distance;
	use General::Arguments;
	use Sequence::Bootstrap;
	use Sequence::Garli;

	# BioPerl libs
	use Bio::TreeIO;
	use Bio::Tree::Tree;
	use Bio::Tree::TreeFunctionsI;
	use Bio::Tree::Node;
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
	use Math::Random::MT::Perl qw(srand rand irand);
	use Data::Dumper;
	use Storable qw(dclone);
	use File::Copy;

	# Non-bio modules
	use Benchmark qw(:all);
	use Getopt::Long;
	use String::Random;

	my $query_aln_file 	= '';
	my $match_aln_file 		= '';
	my $dist_cutoff = 0.50;
	my $bases_cutoff = 0.70;
	my $potential_cutoff = 0.70;
	my $potential_bases = 0.50;



	##################################################################
	# Import alignment
	print "Importing alignment file ".$query_aln_file."...\n";
	my $query_aln = Bio::AlignIO->new(-format => 'fasta',
									-file   => $query_aln_file );
	my $query_aln_obj = $query_aln->next_aln;
	##################################################################
	print "Importing alignment file ".$match_aln_file."...\n";
	my $match_aln = Bio::AlignIO->new(-format => 'fasta',
									-file   => $match_aln_file );
	my $match_aln_obj = $match_aln->next_aln;
	##################################################################

	my @query_sequence_array = ();
	foreach my $seq ($query_aln_obj->each_seq) {
		push(@query_sequence_array,$seq);
	}

	my @match_aln_array = ();
	foreach my $seq ($match_aln_obj->each_seq) {
		push(@match_aln_array,$seq);
	}

	my $outp = 'potential.fas';
	my $outp_pristine = 'pristine.fas';

	open (OUTPRISTINE, '>',$outp_pristine);
	open (OUTP, '>'.$outp);
	foreach my $query_seq (@query_sequence_array) {
		
		my $seq_id = $query_seq->id;
		my $seq_string = $query_seq->seq;
		
		my $temp_string = int(rand(99999));
		my $temp_query = 'temp_query_'.$temp_string.'.fas';
		my $temp_output =  'temp_aln_'.$temp_string.'.fas';
		open (QUERY, '>'.$temp_query);
		foreach my $match_seq (@match_aln_array) {
			print QUERY '>Match|'.$match_seq->id."\n";
			print QUERY $match_seq->seq."\n";
		}
		print QUERY ">".$seq_id."\n";
		print QUERY $seq_string."\n";
		close (QUERY);
		
		my $mafft_string = "mafft --auto --preservecase --adjustdirection --preservecase $temp_query > $temp_output";
		print "Calling $mafft_string\n";
		system($mafft_string);
		
		unlink $temp_query;
		
		my $mafft_aln = Bio::AlignIO->new(-format => 'fasta',
										-file   => $temp_output );
		my $mafft_aln_obj = $mafft_aln->next_aln;
		unlink $temp_output;
		
		my @mafft_aln_array = ();
		foreach my $mafft_seq ($mafft_aln_obj->each_seq) {
			push(@mafft_aln_array,$mafft_seq);
		}
		$mafft_aln_obj ='';
		$mafft_aln = '';
		my $aln_length = length($mafft_aln_array[0]->seq);
		unlink $temp_output;
		# print $mafft_aln_array[0]->seq."\n";
		# print $mafft_aln_array[1]->seq."\n";
		
		my $match_seq_final = $mafft_aln_array[0]->seq;
		$match_seq_final =~ s/-/*/g;
		my $query_seq_final = $mafft_aln_array[1]->seq;
		
		my @unpacked_match = unpack("C*", $match_seq_final);
		my $unpacked_match_ref = \@unpacked_match;
		my @unpacked_query = unpack("C*", $query_seq_final);
		my $unpacked_query_ref = \@unpacked_query;
		
		my $k2p_distance = 0;
		my $transitions = 0;
		my $transversions = 0;
		my $bases_compared = 0;
		($k2p_distance, $transitions,$transversions,$bases_compared) = k2p_no_bs2(\$unpacked_match_ref, \$unpacked_query_ref, $aln_length);
		print $k2p_distance."\n";
		print $transitions."\n";
		print $transversions."\n";
		print $bases_compared."\n";
		$k2p_distance = sprintf("%.2f",$k2p_distance);
		
		if(($k2p_distance < $dist_cutoff) && ($bases_compared/$aln_length > $bases_cutoff)) {
			#print OUTPRISTINE ">".$seq_id."|".$k2p_distance."\n";
			#print OUTPRISTINE $query_seq_final."\n";
			if($match_only == 1) {
				return $query_seq_final;
			}
		} elsif(($k2p_distance < $potential_cutoff) && ($bases_compared/$aln_length > $potential_bases)) {
			#print OUTP ">".$seq_id."|".$k2p_distance."\n";
			#print OUTP $query_seq_final."\n";
		}
		unlink $temp_output;
	}

	close(OUTP);
}