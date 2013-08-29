# Convert fasta to Arlequin
# 
# Copyright (c) 2012, Bryan White, bpcwhite@gmail.com

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
# BEGIN {
	# use Cwd;
	# our $directory = cwd;
	# }
# $directory = ;
use lib 'C:\Users\m0ltop\Dropbox\Code\Bio';
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

# Local modules
require 'Kimura_Distance_C.pl';

# Configuration
my $alignment_file = 'Most_DE_COI_in_GB_fas.fas';
# chomp (my $alignment_file = <>);
# $alignment_file = fix_bold_fasta($alignment_file);

my $alignment = Bio::AlignIO->new(-format => 'fasta',
								-file   => $alignment_file );
my $original_aln = $alignment->next_aln;

my @delimited_alignment = split(".fas",$alignment_file);
my $new_alignment = $delimited_alignment[0]."_renamed.arp";

my @alphabet_list = ('A',	'B',	'C',	'D',	'E',	'F',	'G',	'H',	'I',	'J',	'K',	'L',	'M',	'N',	'O',	'P',	'Q',	'R',	'S',	'T',	'U',	'V',	'W',	'X',	'Y',	'Z');

my @species_names = ();
my @sample_names = ();
my @overall_sequences = ();

foreach my $seq ($original_aln->each_seq) {
	push(@overall_sequences,$seq);
}

# MOVE THIS SO HAPLOTYPES ARE CALCULATED FOR EACH SAMPLE.

foreach my $seq (@overall_sequences) {
	my $seq_id = $seq->id;
	my @split_seq_id = split(/\_/,$seq_id);
	# foreach my $split (@split_seq_id) {
		# print $split."\n";
	# }
	# my $new_id = $split_seq_id[1]."_".$split_seq_id[2]."_".$split_seq_id[3];
	# print $new_id."\n";
	# $seq->id($new_id);
	push(@species_names,$split_seq_id[0]."_".$split_seq_id[1]);
	my $sample_code = substr $split_seq_id[2],0,3;
	push(@sample_names,$split_seq_id[0]."_".$split_seq_id[1]."_".$sample_code);
}

my %unique_species_names = map {$_, 1} @species_names;
my @unique_species_names = keys %unique_species_names;

my %unique_sample_names = map {$_, 1} @sample_names;
my @unique_sample_names = keys %unique_sample_names;

foreach my $unique_species (@unique_species_names) {
	my $current_species_file = $unique_species."_".$new_alignment;
	unlink $current_species_file;
	open (RENAMED, '>>'.$current_species_file);


	my @number_samples = grep {$_ =~ m/$unique_species/ } @unique_sample_names;
	print RENAMED 	
	"[Profile]\n
	Title=\"".$unique_species."\"
	NbSamples=".scalar(@number_samples)."
	GenotypicData=0
	DataType=DNA
	LocusSeparator=NONE
	MissingData='-'
	[Data]
	[[Samples]]\n\n";
	my $num_samples = 0;
	foreach my $unique_sample (@unique_sample_names) {
		# print $unique_sample."\n";
		my %sequence_hash = ();
		my @sample_sequences = ();
		foreach my $seq (@overall_sequences) {
			if($seq->id =~ m/$unique_sample/) {
				my $found_haplotype = 0;
				if(exists($sequence_hash{$seq->seq})) {
					$sequence_hash{$seq->seq} = $sequence_hash{$seq->seq} + 1;
				} else {
					for my $hash_seq ( sort keys %sequence_hash ) {
						my $critical_value = 0.05;
						my $cutoff = 0.02;
						my $search_type = 1;
						my $max_seq_length = 350;
						my $query_seq = $seq->seq;
						my $hash_seq_filtered = $hash_seq;
						$hash_seq_filtered =~ s/-/*/g;
						my ($transitions,	$transversions,		$bases_compared,
						$k2p_distance,	$variance,			$stderror,
						$mink2p,		$maxk2p
						) = 0;
						c_kimura_distance(	$query_seq,			$hash_seq_filtered,	$critical_value,
											$cutoff, 			$search_type, 		$max_seq_length,
											$transitions,		$transversions,		$bases_compared,
											$k2p_distance,		$variance,			$stderror,
											$mink2p,			$maxk2p);
						# print $k2p_distance."\n";
						if($k2p_distance <= 0) {
							$sequence_hash{$hash_seq} = $sequence_hash{$hash_seq} + 1;
							$found_haplotype = 1;
							last;
						}
					}
					if($found_haplotype == 0) {
						$sequence_hash{$seq->seq} = 1;
						push(@sample_sequences,$seq);
					}
				}
			}
		}
		foreach my $seq (@sample_sequences) {
			# if($sequence_hash{$seq->seq} > 1) { print $sequence_hash{$seq->seq}."\n" };
			$seq->description($sequence_hash{$seq->seq});
		}
	
		
		if($unique_sample =~ m/$unique_species/) {
			my $num_individuals = 0;
			$num_samples++;
			my @sample_sequence_lines = ();
			# print $unique_sample."\n";
			foreach my $seq (@sample_sequences) {
				# print "\t".$seq->id."\n";
				if($seq->id =~ m/$unique_sample/) {
					push(@sample_sequence_lines,$seq->id."\t\t".$seq->description."   ".$seq->seq);
					$num_individuals++;
					# print "\t".$num_individuals."\n";
				}
			}
			print RENAMED "SampleName=\"".$unique_sample."\"\n";
			print RENAMED "SampleSize=".$num_individuals."\n";
			# if($num_individuals > 2) { print $unique_sample."\n" };
			print RENAMED "SampleData= {\n";
			foreach my $line (@sample_sequence_lines) {
				print RENAMED $line."\n";
			}
			print RENAMED "}\n";
		}			
	}
	print RENAMED 
		"[[Structure]]
			StructureName=\"Group Structure Assignments\"
			NbGroups=".scalar(@number_samples)."\n";
	my $sample_i = 0;		
	foreach my $sample_name (@unique_sample_names) {
		if($sample_name =~ $unique_species) {
			print RENAMED 
		"\n#$alphabet_list[$sample_i]
		\tGroup={
		\t\"$sample_name\"
		\t}
		\t";
			$sample_i++;
		}
	}
	print "Number of samples ".$unique_species.",".scalar(@number_samples).",".$num_samples."\n";
	close (RENAMED);
}

