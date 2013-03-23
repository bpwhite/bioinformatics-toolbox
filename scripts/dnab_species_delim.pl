#!/usr/bin/perl
# This script takes an input of a FASTA sequence and outputs a 
# species delimitation regime given certain parameters.
# It is designed for use with the cytochrome oxidase I (COI) gene
# and to provide species delimitations useful for DNA barcoding.

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
use lib "$FindBin::Bin/libs/Sequence"; 

# Import sequence libs
use Fasta;
require 'Kimura_Distance_C.pl';

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

# Non-bio modules
use Statistics::Descriptive;
# use Math::Random::MT qw(srand rand);
use Benchmark qw(:all);
use POSIX qw(ceil);
use File::Copy;

use strict;
use warnings;


# bit_Test();
# exit;
# my $alignment_file = "preservation_all_08_07_2012.fas";
# my $alignment_file = "just_brachycentridae_aln.fas";
# my $alignment_file = "small_test_aln.fas";
# my $alignment_file = "sccwrp_everything_07_15_2012_aln.fas";
# my $alignment_file = "baetis_07_16_2012_aln.fas";
# my $alignment_file = "trichops_combined_7_16_2012.fas";
# my $alignment_file = "brachycentridae_lepido_serico_07_17_2012_aln.fas";
# my $alignment_file = "leptohyphidae_07_17_2012_aln.fas";
# my $alignment_file = "sccwrp_library_07_17_2012_aln.fas";
# my $alignment_file = "ephemeroptera_sccwrp_public_07_22_2012_aln.fas";
# my $alignment_file = "ephemerella.fas";
# my $alignment_file = "baetis_07_16_2012_unique.fas";
# my $alignment_file = "ephemeroptera_sccwrp_07_26_2012.fas";
# my $alignment_file = "chironomidae_sccwrp_07_26_2012.fas";
# my $alignment_file = "sccwrp_coleoptera_07_27_2012.fas";
# my $alignment_file = "diptera_all_bold_500bp.fas";
# my $alignment_file = "churchill_trichoptera.fas";
# my $alignment_file = "churchill_plecoptera.fas";
# my $alignment_file = "churchill_ephemeroptera.fas";
# my $alignment_file = "churchill_combined.fas";
# my $alignment_file = "ephemeroptera_reference_combined.fas";
# my $alignment_file = "coleoptera_psephenus.fas";
# my $alignment_file = "psephenus_public.fas";
# my $alignment_file = "mopalia.fas";
# my $alignment_file = "BOLD_Ephemeroptera_100bp.fas";
# my $alignment_file = "stroud_combined.fas";
# my $alignment_file = "druiter_trichops.fas";
# my $alignment_file = "sccwrp_hydromod_09_26_2012.fas";
# my $alignment_file = "coleoptera_500bp.fas";
# my $alignment_file = "polyplacophora_500bp.fas";

my @alignment_files = <*.fas>;
my $alignment_number = 0;
foreach my $alignment_file (@alignment_files) {
	print $alignment_number.". ".$alignment_file."\n";
	$alignment_number++;
}
chomp(my $selected_alignment = <>);
my $alignment_file = $alignment_files[($selected_alignment-1)];
# exit;
##################################################################
# Parameters
my $skip_phylo = 'skip';
# my $skip_phylo = 'no';
my $minimum_phylo_alleles = 10;

my %use_tags_hash =  (  1 => 1,	 # Use tags
						2 => 2); # Don't use tags
# my $use_tags = $use_tags_hash{1};
my $use_tags = $use_tags_hash{1};

my %cutoff_hash = ( 1 => '0.001',
					2 => '0.005',
					3 => '0.010',
					4 => '0.020',
					5 => '0.022',
					6 => '0.030',
					7 => '0.040');
# my $cutoff = $cutoff_hash{1};
# my $cutoff = $cutoff_hash{3};
# my $cutoff = $cutoff_hash{4};
# my $cutoff = $cutoff_hash{5};
# my $cutoff = $cutoff_hash{6};
my $cutoff = 0.02;

my %sequence_length_hash = (1 => 25,
							2 => 50,
							3 => 75,
							4 => 125,
							5 => 250,
							6 => 350,
							7 => 400,
							8 => 450,
							9 => 500,
							10 => 650,
							11 => 658);
# my $minimum_sequence_length = $sequence_length_hash{1};
# my $minimum_sequence_length = $sequence_length_hash{4};
my $minimum_sequence_length = 350;

my $phylogenetic_length_cutoff = 350;

my %statistcal_method_hash = ( 1 => 'strict',
								2 => 'analytical',
								3 => 'bootstrap');
# my $statistical_method = $statistcal_method_hash{1};
my $statistical_method = $statistcal_method_hash{2};
# my $statistical_method = $statistcal_method_hash{3};

my 	%bootstrap_reps_hash = (1 => 25,
							2 => 50,
							3 => 75,
							4 => 100,
							5 => 200,
							6 => 300,
							7 => 400,
							8 => 500,
							9 => 1000);
my $num_bootstrap_reps = $bootstrap_reps_hash{5}; # Minimum should be 25

my $minimum_catch_distance = 0.05;
my $num_examplars = 5;

# 2-tailed crit values. Normal distribution
my $alpha_level = 0.05;
my $critical_value = 1.96; # Default
if($alpha_level == 0.01)		{
	$critical_value = 2.575;
}
elsif($alpha_level == 0.025)	{
	$critical_value = 2.2414;
}
elsif($alpha_level == 0.05) 	{
	$critical_value = 1.96;
}
elsif($alpha_level == 0.075) 	{
	$critical_value = 1.7805;
}
elsif($alpha_level == 0.10) 	{ 
	$critical_value = 1.645;
}
elsif($alpha_level == 0.125) 	{ 
	$critical_value = 1.5341;
}
elsif($alpha_level == 0.15) 	{ 
	$critical_value = 1.4395;
}
elsif($alpha_level == 0.175) 	{ 
	$critical_value = 1.3563;
}
elsif($alpha_level == 0.20) 	{
	$critical_value = 1.28;
}

my %intra_otu_distances_hash = ( 1 => 1,
								 2 => 2);
my $intra_otu_distances = $intra_otu_distances_hash{1};

my $interactive_mode = 0;

fix_bold_fasta($alignment_file);

my @alignment_file_split = split(".fas",$alignment_file);
my $alignment_label = $alignment_file_split[0];
my $output_prefix = $alignment_label.'_'.$cutoff.'_'.$minimum_sequence_length.'_'.$statistical_method.'_'.$num_bootstrap_reps;


my $output_path = $alignment_label."\\";
##################################################################

##################################################################
# Check if output path exists, if not, create the output path.
unless(-d ($output_path)) {
	print "Creating output path: ".$output_path."\n";
	mkdir $output_path;
}
##################################################################

##################################################################
# Import alignment
print "Importing alignment file ".$alignment_file."...\n";
my $alignin = Bio::AlignIO->new(-format => 'fasta',
								-file   => $alignment_file );
my $original_aln = $alignin->next_aln;
##################################################################

##################################################################
# Start benchmark
my $t0 = Benchmark->new;
my $k2p1 = 0;
##################################################################

##################################################################
# Tag sequences
my @original_sequence_array = ();
my $current_tag = '';
if ($use_tags == 1) {
	print "Tagging...\n";
	# For each sequence, if the sequence is a >TAG|NAME, begin tagging
	foreach my $seq ($original_aln->each_seq) {
		if($seq->id =~ m/TAG/) {
			my @delimited_tag = split(/\|/,$seq->id);
			$current_tag = $delimited_tag[1];
			print $current_tag."\n";
			next;
		}
		# Split the ID, count the delimiters, attach the tag to the end
		my @delimited_seq_id = split(/\|/,$seq->id);
		my $new_id = '';
		my $num_delimiters = scalar(@delimited_seq_id);
		if ('COI-5P' ~~ @delimited_seq_id) { $num_delimiters-- } ;
		my $delimiter_i = 1;
		foreach my $split_id (@delimited_seq_id) {
			if(($split_id ne 'COI-5P') && ($delimiter_i != $num_delimiters)) {
				$new_id .= $split_id."|";
			} elsif (($split_id ne 'COI-5P') && ($delimiter_i == $num_delimiters)) {
				$new_id .= $split_id."_";
			}
			$delimiter_i++;
		}
		$new_id .= $current_tag;
		$seq->id($new_id);
		# print $seq->id."\n";
		push(@original_sequence_array,$seq);
	}
} else {
	# No tag, just build the initial sequence array.
	foreach my $seq ($original_aln->each_seq) {
		push(@original_sequence_array,$seq);
	}
}
##################################################################

##################################################################
# 1. Delete identical sequences
# 2. Filter sequences
# 3. Remove sequences less than length cutoff
# 4. Replace ambiguous characters with gaps (-)
my %unique_sequences = (); # Don't flush
my @query_seqs_array = ();
my $max_seq_length = 0;
my $alignment_length = 0;
my @ambiguous_characters = ('N', 'R', 'Y', 'K', 
							'M', 'S', 'W', 'B', 
							'D', 'H');
my %outgroup_seqs = ();
foreach my $seq (@original_sequence_array) {
	my $seq_gapped = $seq->seq();
	my $seq_id = $seq->id();
	foreach my $ambiguous_character (@ambiguous_characters) {
		$seq_gapped =~ s/$ambiguous_character/-/g;
	}
	if(fast_seq_length($seq_gapped) < $minimum_sequence_length) { next; };
	if($seq_id =~ m/Outgroup/) {
		$outgroup_seqs{$seq_gapped} = $seq_id;
		next;
	}
	my $seq_degapped = $seq_gapped;
	$seq_degapped =~ s/-/ /g;
	$seq_degapped =~ s/\s+//g;
	
	# Set max sequence length (Alignment length)
	if(length($seq_degapped) > $max_seq_length) { $max_seq_length = length($seq_degapped) } ;
	$alignment_length = length($seq_gapped);
	my $filtered_seq = $seq_gapped;
	$filtered_seq =~ s/-/*/g;
	$seq->seq($seq_gapped);
	$seq->description($seq_degapped);
	$seq->object_id($filtered_seq);
	$unique_sequences{$seq->seq()} = $seq->id;
	push(@query_seqs_array,$seq);
}
##################################################################

##################################################################
my @unique_sequence_array = ();
for my $seq_string ( sort keys %unique_sequences ) {
	if(fast_seq_length($seq_string) < 1) { next; };
	my $seq_id = $unique_sequences{$seq_string};
	my @query_seq = grep { $_->id eq $seq_id } @original_sequence_array;
	push(@unique_sequence_array,$query_seq[0]);
}
print scalar(@unique_sequence_array)." unique sequences remaining.\n";
##################################################################

##################################################################
my %original_seq_hash 		= ();
my $original_seq_hash_ref 	= \%original_seq_hash;
foreach my $seq (@unique_sequence_array) {
	$original_seq_hash_ref->{$seq->id}->{'gapped_seq'} = $seq->seq();
	$original_seq_hash_ref->{$seq->id}->{'filtered_seq'} = $seq->object_id();
}
@unique_sequence_array = (); # Flush
##################################################################

##################################################################
print "Deleting highly similar sequences...\n";
my %seq_hash1 = ();
%seq_hash1 = %original_seq_hash;
my %seq_hash2 = %original_seq_hash;
my $seq_hash1_ref = \%seq_hash1;
my $seq_hash2_ref = \%seq_hash2;
my $seq_to_delete 			= "";
my $seq_i 					= 1;
my $total_seq_comparisons 	= 0;
my $remaining_unique_seqs 	= keys %original_seq_hash;
my $num_unique_seqs 		= keys %seq_hash1;
%original_seq_hash 			= (); # Flush
for my $seq_id1 ( sort keys %seq_hash1 ) {
	my $max_bases_compared = 0;
	my $seq1 = $seq_hash1_ref->{$seq_id1}->{'filtered_seq'};
	my $seq1_gapped = $seq_hash1_ref->{$seq_id1}->{'gapped_seq'};
	my $seq1_length = fast_seq_length($seq1_gapped);
	for my $seq_id2 ( sort keys %seq_hash2 ) {
		next if $seq_id1 eq $seq_id2;
		my $seq2_gapped = $seq_hash1_ref->{$seq_id2}->{'gapped_seq'};
		my ($transitions,	$transversions,		$bases_compared,
			$k2p_distance,	$variance,			$stderror,
			$mink2p,		$maxk2p,			$p_stderror,
			$p_min,			$p_max,				$p_dist
			) = 0;
		my $search_type = 1;
		c_kimura_distance(	$seq1,				$seq2_gapped,		$critical_value,
							$cutoff, 			$search_type, 		$max_seq_length,
							$transitions,		$transversions,		$bases_compared,
							$k2p_distance,		$variance,			$stderror,
							$mink2p,			$maxk2p);
		next if($bases_compared < $minimum_sequence_length);
		if ($k2p_distance <= $cutoff) {
			$seq_to_delete = $seq_id1;
			delete $seq_hash1_ref->{$seq_to_delete};
			delete $seq_hash2_ref->{$seq_to_delete};
			$seq_to_delete = "";
			last;
		}
	}
	$seq_i++;
}

print "\n";
print "Done deleting.\n";
my $remaining_otu = keys %seq_hash1;
print $remaining_otu."\n";

print "Rebuilding sequence arrays...\n";
my $number_remaining_seqs = 0;
##################################################################

##################################################################
# Rebuilt sequence object array for predicted OTU's
my @otu_seqs_array = ();
for my $seq_id (sort keys %seq_hash1) {
	$number_remaining_seqs++;
	foreach my $original_seq (@query_seqs_array) {
		if ($original_seq->id eq $seq_id) {
			push(@otu_seqs_array,$original_seq);
		}
	}
}
%seq_hash1 = (); # Flush
%seq_hash2 = (); # Flush

print "Number predicted OTU's: ".$number_remaining_seqs."\n";


# Begin rebuilding
print "Creating output file...\n";

# my $matches_file = $output_prefix.'_matches.csv';
# unlink $output_path.$matches_file;
# open(MATCHES, '>>'.$output_path.$matches_file);

my $number_query_seqs 	= scalar(@query_seqs_array);
my $query_seq_i 		= 1;

print "Re-Matching deleted sequences to OTU's...\n";
my %query_results_hash 	= (); # Contains each query seq and its corresponding otu matches
my @new_otu_seqs = ();
my $ten_pcnt_search = $number_query_seqs*.10;
foreach my $query_seq (@query_seqs_array) {
	my $query_seq_was_matched 	= 0;
	my $closest_match 			= '';
	my $lowest_distance 		= 2;
	my $lowest_min 				= 2;
	my $lowest_max 				= 2;
	my $lowest_bases_compared 	= 0;
	my $seq1_filtered 			= $query_seq->object_id();
	my $seq1_gapped				= $query_seq->seq();
	foreach my $otu_seq (@otu_seqs_array) {
		my $seq2_gapped = $otu_seq->seq();
		my ($transitions,	$transversions,		$num_bases_compared,
			$current_dist,	$used_ts_shortcut,	$variance,
			$stderror,		$mink2p,			$maxk2p
			) = 0;
		my $calculation_rounding = "%.5f";
		
		if($seq1_gapped eq $seq2_gapped) {
			$mink2p = 0;
			$maxk2p = 0;
			$stderror = 0;
			$current_dist = 0;
			push(@{$query_results_hash{$query_seq->id}},$otu_seq->id);
			$query_seq_was_matched++;
			next;
		} else {
			if($statistical_method eq 'strict') {
					($transitions,		$transversions,		$num_bases_compared,
					$current_dist,		$variance,			$stderror,
					$mink2p,			$maxk2p
					) = 0;
					my $search_type = 2;
					c_kimura_distance(	$seq1_filtered,		$seq2_gapped,		$critical_value,
										$cutoff,			$search_type, 		$max_seq_length,
										$transitions,		$transversions,		$num_bases_compared,
										$current_dist,		$variance,			$stderror,
										$mink2p,			$maxk2p);
					# Update lowest distance
					if($current_dist < $lowest_distance) {
						$lowest_distance 		= $current_dist;
						$lowest_min 			= $mink2p;
						$lowest_max 			= $maxk2p;
						$closest_match 			= $otu_seq->id;
						$lowest_bases_compared 	= $num_bases_compared;
					}
					next if $num_bases_compared < $minimum_sequence_length;
			} elsif($statistical_method eq 'analytical') {
					($transitions,		$transversions,		$num_bases_compared,
					$current_dist,		$variance,			$stderror,
					$mink2p,			$maxk2p
					) = 0;
					my $search_type = 2;
					c_kimura_distance(	$seq1_filtered,		$seq2_gapped,		$critical_value,
										$cutoff, 			$search_type, 		$max_seq_length,
										$transitions,		$transversions,		$num_bases_compared,
										$current_dist,		$variance,			$stderror,
										$mink2p,			$maxk2p);
					
					if($current_dist < $lowest_distance) {
						$lowest_distance 		= $current_dist;
						$lowest_min 			= $mink2p;
						$lowest_max 			= $maxk2p;
						$closest_match 			= $otu_seq->id;
						$lowest_bases_compared 	= $num_bases_compared;
					}
					next if $num_bases_compared < $minimum_sequence_length;
			}
		}
		# If the method is strict and a query has already matched an OTU
		# Strict means only ever match 1 OTU.
		last if (($statistical_method eq 'strict') && ($query_seq_was_matched > 0));

		# If minimum distance is less than cutoff and statistical method is employed, 
		# add the match to the match array for this query sequence.
		if ($mink2p <= $cutoff) {
			push(@{$query_results_hash{$query_seq->id}},$otu_seq->id);
			$query_seq_was_matched++;
		}
	}
	
	# If strict method or no valid minimum match was found, add the closest match to the array.
	if($query_seq_was_matched  == 0) {
		# if($lowest_min < $cutoff) {
			push(@{$query_results_hash{$query_seq->id}},$closest_match);
		# } else {
			# push(@new_otu_seqs,$query_seq);
			# push(@{$query_results_hash{$query_seq->id}},$query_seq->id);
		# }
	}
	$query_seq_i++;
}
push(@otu_seqs_array,@new_otu_seqs);
# close(MATCHES);

# OTU summary and content file names
my $otu_summary_file = $output_prefix.'_otu_summary.csv';
my $otu_summary_excel = $output_prefix.'_otu_summary_excel.csv';
unlink $output_path.$otu_summary_file;
open(OTU_SUMMARY, '>>'.$output_path.$otu_summary_file);
unlink $output_path.$otu_summary_excel;
open(OTU_EXCEL, '>>'.$output_path.$otu_summary_excel);

print OTU_SUMMARY "OTU Results: ".($cutoff*100)."% cutoff\n";
# print OTU_SUMMARY "User: Bryan White\n";
print OTU_SUMMARY "\tAlignment: ".$alignment_file."\n";
print OTU_SUMMARY "\tSequences: ".$number_query_seqs."\n";
print OTU_SUMMARY "\tStatistical method: ".$statistical_method."\n";
if($statistical_method eq 'bootstrap') {
	print OTU_SUMMARY "\tBootstrap reps: ".$num_bootstrap_reps."\n";
}
print OTU_SUMMARY "\tMinimum seq. length: ".$minimum_sequence_length."\n";
print OTU_SUMMARY "\n";

print "OTU Results: ".($cutoff*100)."% cutoff\n";
# print "User: Bryan White\n";
print "\tAlignment: ".$alignment_file."\n";
print "\tSequences: ".$number_query_seqs."\n";
print "\tStatistical method: ".$statistical_method."\n";
if($statistical_method eq 'bootstrap') {
	print "\tBootstrap reps: ".$num_bootstrap_reps."\n";
}
print "\tMinimum seq. length: ".$minimum_sequence_length."\n";
print "\n";

my $otu_i = 1;
# my $num_otus = scalar(@otu_seqs_array);

my $last_otu_seq = '';
my $string_space = "%-60s %-11s %-30s %-50s %-30s %-33s %-33s %-30s %-12s %-12s %-17s %-100s\n";

printf	$string_space,
		"[OTU #] OTU Ref. ID",
		"[# Morpho]",
		"[A - U - D]",
		"[p - C - CR - E - ER - AE]",
		"Avg. Dist SE: (Min - Max)",
		"Avg. Comparisons SE: (Min - Max)",
		"Avg. Length SE: (Min - Max)",		
		"Avg. SE SE: (Min - Max)",
		"[Sub-groups]",
		"[Link Depth]",
		"[Link Strength %]",
		"Nearest Neighbor Dst (Min - Max)";
my $delimiter = ",";
# my $delimiter = '';
print OTU_SUMMARY
		"[OTU #] OTU Ref. ID".$delimiter.
		"[# Morpho]".$delimiter.
		"[A - U - D]".$delimiter.
		"[p - C - CR - E - ER - AE]".$delimiter.
		"Avg. Dist SE: (Min - Max)".$delimiter.
		"Avg. Comparisons SE: (Min - Max)".$delimiter.
		"Avg. Length SE: (Min - Max)".$delimiter.
		"Avg. SE SE: (Min - Max)".$delimiter.
		"[Sub-groups]".$delimiter.
		"[Link Depth]".$delimiter.
		"[Link Strength %]".$delimiter.
		"Nearest Neighbor Dst (Min - Max)".$delimiter."\n";
print OTU_EXCEL
		"OTU #".$delimiter."OTU Ref. ID".$delimiter.
		"# Morpho".$delimiter.
		"Abundance".$delimiter."Unique Alleles".$delimiter."Distinct Alleles".$delimiter.
		"p-value".$delimiter."ml_clusters".$delimiter."ml_conf".$delimiter."ml_entities".$delimiter."entity_conf".$delimiter."adjusted_entites".$delimiter.
		"Avg. Dist".$delimiter."Dist SE ".$delimiter."Min Dist".$delimiter."Max".$delimiter.
		"Avg. Comparisons".$delimiter." Comparisons SE".$delimiter."Min comparisons".$delimiter."Max comparisons".$delimiter.
		"Avg. Length".$delimiter."Length SE".$delimiter."Min length".$delimiter."Max length".$delimiter.
		"Avg. SE".$delimiter. "SE SE".$delimiter."Min SE".$delimiter."Max SE".$delimiter.
		"Sub-groups".$delimiter.
		"Link Depth".$delimiter.
		"Link Strength Percent".$delimiter.
		"Nearest Neighbor".$delimiter."Dist".$delimiter."Min".$delimiter."Max".$delimiter."\n";

# Output the OTU content (match) results
my $otu_results = $output_prefix.'_otu_results.csv';
unlink $output_path.$otu_results;
open(OTU_RESULTS, '>>'.$output_path.$otu_results);
print OTU_RESULTS 'otu_seq,query_seq'."\n";

my $gmyc_results_filename = $output_prefix.'_gmyc_otu_results.csv';
unlink $output_path.$gmyc_results_filename;
open(OTU_GMYC_RESULTS, '>>'.$output_path.$gmyc_results_filename);
print OTU_GMYC_RESULTS 'otu_id,otu_gmyc_seq,query_seq,pvalue'."\n";

my $otu_exemplars = $output_prefix.'_exemplars.fas';
unlink $output_path.$otu_exemplars;
open(EXEMPLARS, '>>'.$output_path.$otu_exemplars);

my @query_seqs_found = ();
my @found_links_OTUs = ();
my $num_exemplars = 0;

my %morpho_name_hash = ();
my $morpho_name_hash_ref = \%morpho_name_hash;
my @otu_morpho_lumps = ();
my @gmyc_morpho_lumps = ();
foreach my $otu_seq (@otu_seqs_array) { # For each OTU
	next if $otu_seq->id ~~ @found_links_OTUs; # Skip found OTU's
	my @overall_query_matches = (); # Store all of the query matches from
									# all linked clusters in this array.
	push(@found_links_OTUs,$otu_seq->id); # Don't search this OTU again.
	
	# For each otu, cycle through the list of query matches
	# Accumulate the sequences that make up this OTU
	my $query_matches_ref 	= find_query_matches(\%query_results_hash,$otu_seq->id);
	@overall_query_matches 	= @$query_matches_ref;
	# Initial number of queries found.
	my @link_strength = ();
	push(@link_strength,scalar(@overall_query_matches));
	# Cycle back through the list of found query seqs that are in this otu
	# Accumulate the OTU's that are shared by these sequences
	# Must also recursively accumulate OTU's
	my @otu_links_array 				= ();
	my $previous_number_queries_found 	= 0;
	my $current_number_queries_found 	= 0;
	
	# Recursively accumulate query sequences for OTUs by following links.
	my $link_depth = 0;
	
accumulate_linked_OTUs:
	# Reduce the overall query matches down to unique matches
	my %unique_overall_query_matches = map {$_, 1} @overall_query_matches;
	my @unique_overall_query_matches = keys %unique_overall_query_matches;
	
	$previous_number_queries_found = scalar(@unique_overall_query_matches);
	
	while (my ($query_id,$matched_otu_array) = each (%query_results_hash)) { # For each query seq
		foreach my $query_match (@overall_query_matches) {
			if ($query_match eq $query_id) {
				foreach my $matched_otu_link (@$matched_otu_array) {
					push(@otu_links_array,$matched_otu_link);
				}
			}
		}
	}
	
	# Reduce linked OTU list down to a unique OTU list
	# Iterate through that unique OTU list and find all of
	# the query result arrays for each of the linked OTUs.
	# Push them all into the overall query matches array.
	my %unique_otu_links   = map { $_, 1 } @otu_links_array;
	my @unique_otu_links = keys %unique_otu_links;
	
	foreach my $unique_linked_otu (@unique_otu_links) {
		push(@found_links_OTUs,$unique_linked_otu); # Don't search this OTU again.
		my $linked_query_matches_ref	= find_query_matches(\%query_results_hash,$unique_linked_otu);
		my @linked_query_matches 		= @$linked_query_matches_ref;
		push(@overall_query_matches,@linked_query_matches);
	}
	%unique_overall_query_matches = map {$_, 1} @overall_query_matches;
	@unique_overall_query_matches = keys %unique_overall_query_matches;
	$current_number_queries_found = scalar(@unique_overall_query_matches);

	# If there were still new queries found, run back again and find more.
	# If no more queries were found, done and script can continue.
	my $added_queries = $current_number_queries_found - $previous_number_queries_found;
	push(@link_strength,$added_queries);
	$link_depth++;
	if ($previous_number_queries_found !=  $current_number_queries_found) { goto accumulate_linked_OTUs } ;
	
	
	##################################################################
	## End accumulation and flush variables ##########################
	@otu_links_array 				= (); # Flush
	@overall_query_matches 			= (); # Flush
	%unique_overall_query_matches 	= (); # Flush
	%unique_otu_links				= (); # Flush
	$query_matches_ref				= ''; # Flush
	##################################################################
	
	##################################################################
	# Collect the sequences for the current OTU into this array, current_otu_sequences
	# Output the OTU content (match) results
	my %current_otu_sequences 		= ();
	my %current_otu_seqs_and_id		= ();
	my %exemplars_hash 				= ();
	my $current_otu_output = $output_prefix.'_'.$otu_i.'.fas';
	unlink $output_path.$current_otu_output;
	open(OTU_FASTA, '>>'.$output_path.$current_otu_output);
	my @seq_lengths = ();
	foreach my $query_seq_id (@unique_overall_query_matches) {
		my @query_seq = grep { $_->id eq $query_seq_id } @query_seqs_array;
		$current_otu_sequences{$query_seq[0]->seq()} = $query_seq[0]->object_id();
		$current_otu_seqs_and_id{$query_seq_id} = $query_seq[0]->seq();
		$exemplars_hash{$query_seq[0]->id} = $query_seq[0]->seq();
		my $current_seq_length = fast_seq_length($query_seq[0]->seq());
		push(@seq_lengths,$current_seq_length);
		my $filtered_otu_id = filter_one_id($otu_seq->id);
		my $filtered_query_id = filter_one_id($query_seq[0]->id());
		
		print OTU_RESULTS $filtered_otu_id.','.$filtered_query_id."\n";
		print OTU_FASTA '>['.$otu_i.']_'.$query_seq_id."\n";
		print OTU_FASTA $query_seq[0]->seq()."\n";
	}
	close(OTU_FASTA);
	##################################################################
	
	##################################################################
	## Begin finding exemplars
	# Sort sequence lengths and print exemplars.
	my @sorted_seq_lengths = (sort { $b <=> $a } @seq_lengths);
	my @printed_exemplars = ();
	my @exemplar_keys = keys %exemplars_hash;
	fisher_yates_shuffle( \@exemplar_keys );    # permutes @array in place
	for (my $exemplar_i = 0; ($exemplar_i < $num_examplars) && ($exemplar_i < scalar(@seq_lengths)); $exemplar_i++) {
		foreach my $seq_id (@exemplar_keys) {
			next if $seq_id ~~ @printed_exemplars;
			if(fast_seq_length($exemplars_hash{$seq_id}) == $sorted_seq_lengths[$exemplar_i]) {
				print EXEMPLARS '>['.$otu_i.']_'.$seq_id."\n";
				print EXEMPLARS $exemplars_hash{$seq_id}."\n";
				push(@printed_exemplars,$seq_id);
				$num_exemplars++;
				last;
			}
		}
	}
	my $first_exemplar_seq_gapped = $exemplars_hash{$printed_exemplars[0]};
	%exemplars_hash 		= (); # Flush
	@exemplar_keys 			= (); # Flush
	@sorted_seq_lengths 	= (); # Flush
	@printed_exemplars 		= (); # Flush
	##################################################################
	
	##################################################################
	# Use first exemplar to find nearest neighbor.
	my $nn_id = '';
	my $nn_dist = 100;
	my $nn_min = 0;
	my $nn_max = 0;
	my $nn_sequence = '';
	foreach my $otu_seq (@otu_seqs_array) { # For each OTU
		next if $otu_seq->id ~~ @unique_otu_links; # Skip found OTU's
		my $neighbor_seq_filtered = $otu_seq->object_id;
		my ($transitions,	$transversions,	$num_bases_compared,	$current_dist,
			$variance,		$stderror,		$mink2p,				$maxk2p,$p_stderror,
			$p_min,			$p_max,			$p_dist
			) = 0;
			my $search_type = 3;
			c_kimura_distance(	$neighbor_seq_filtered,	$first_exemplar_seq_gapped,	$critical_value,
								$cutoff, 				$search_type, 				$max_seq_length,
								$transitions,			$transversions,				$num_bases_compared,
								$current_dist,			$variance,					$stderror,
								$mink2p,				$maxk2p);
			next if $num_bases_compared < $minimum_sequence_length;
		if($current_dist < $nn_dist) {
			$nn_dist = $current_dist;
			$nn_id = $otu_seq->id;
			$nn_min = $mink2p;
			$nn_max = $maxk2p;
			$nn_sequence = $otu_seq->seq();
		}
	}
	my $nn_rounding = "%.3f";
	$nn_id 		= filter_one_id($nn_id);
	$nn_dist 	= sprintf($nn_rounding,$nn_dist)*100;
	$nn_min 	= sprintf($nn_rounding,$nn_min)*100;
	$nn_max 	= sprintf($nn_rounding,$nn_max)*100;
	##################################################################
	## End nearest neighbor search

	##################################################################
	# Collect sequence lengths and distances.
	my @overall_names 			= ();
	my @num_comparisons			= ();
	my %unique_alleles			= ();
	my %distinct_alleles		= ();
	my @distances 				= ();
	my @std_errors				= ();
	my @deflated_dists			= ();
	my @deflated_stderrors		= ();
	my $num_unique_oqm 			= keys %current_otu_sequences;
	my $matrix_count 			= ($num_unique_oqm*$num_unique_oqm-$num_unique_oqm)/2;
	my $maximum_dist			= 0;
	my $minimum_possible_max	= 0;
	my $unique_oqm_i 			= 1;
	for my $seq1_gapped ( sort keys %current_otu_sequences ) {
		my $seq1_filtered = $current_otu_sequences{$seq1_gapped};
		for my $seq2_gapped ( sort keys %current_otu_sequences ) {
			$unique_alleles{$seq1_gapped} = 'a';
			$unique_alleles{$seq2_gapped} = 'a';

			my ($transitions,	$transversions,	$num_bases_compared,$current_dist,
				$used_ts_shortcut,
				$variance,		$stderror,		$mink2p,			$maxk2p,$p_stderror,
				$p_min,			$p_max,			$p_dist) = 0;
			my $search_type = 3;
				c_kimura_distance(	$seq1_filtered,	$seq2_gapped,	$critical_value,
									$cutoff, $search_type, $max_seq_length,
									$transitions,	$transversions,	$num_bases_compared,
									$current_dist,	$variance,		$stderror,
									$mink2p,		$maxk2p);
			if ($num_bases_compared < $minimum_sequence_length) { $unique_oqm_i++; next };
			if($current_dist > 0) {
				# Store distinct allele sequences greater than 0 dist.
				$distinct_alleles{$seq1_gapped} = 'a';
				$distinct_alleles{$seq2_gapped} = 'a';
			}
			push(@distances, $current_dist);
			push(@std_errors, $stderror);
			push(@num_comparisons, $num_bases_compared);
			last if $unique_oqm_i == $matrix_count;
			$unique_oqm_i++;
		}
		last if $unique_oqm_i == $matrix_count;
	}
	##################################################################
	## End intra-OTU distance calculations
	%current_otu_sequences = (); # Flush
	##################################################################
	
	##################################################################
	# Compute overall sub-morpho names
	foreach my $unique_o_q_match1 (@unique_overall_query_matches) {
		my $current_name 	= filter_one_id($unique_o_q_match1);
		$current_name 		= convert_id_to_name($current_name);
		push(@overall_names,$current_name);
	}
	##################################################################
	
	##################################################################
	##################################################################
	##################################################################
	my $p_value = '';
	my $ml_clusters = '';
	my $ml_clusters_conf = '';
	my $ml_entities = '';
	my $ml_entities_conf = '';
	my $adjusted_entities = 0;
	my $skipped_phylo = 0;
	my $num_specimens_assigned_gmyc = 0;
	my %gmyc_name_id_hash = ();
	##################################################################

	# Nexus output
	# Print distinct alleles for phylogenetic analysis
	# Build ML tree -> MrBayes
	# GMYC analysis
	if($skip_phylo eq 'skip') { 
		$skipped_phylo = 1;
		goto skip_phylo;
	}
	if(keys %distinct_alleles < $minimum_phylo_alleles) {
		$skipped_phylo = 1;
		goto skip_phylo;
	}
	my $number_of_phylo_worthy_taxa = 0;
	for my $seq_gapped (sort keys %distinct_alleles) {
		if(fast_seq_length($seq_gapped) > $phylogenetic_length_cutoff) {
			$number_of_phylo_worthy_taxa++;
		}
	}
	$number_of_phylo_worthy_taxa += 1;
	my $current_otu_distincts 	= $output_prefix.'_distinct_'.$otu_i.'.nex';
	unlink $current_otu_distincts;
	open(OTU_DISTINCTS, '>>'.$current_otu_distincts);
	my $distincts_filename = $current_otu_distincts;
	my $log_filename = $output_prefix.'_log_distinct_'.$otu_i.'.txt';
	# Print nexus header
	print OTU_DISTINCTS "# NEXUS\n";
	print OTU_DISTINCTS "begin data;\n";
	print OTU_DISTINCTS	"dimensions ntax=".$number_of_phylo_worthy_taxa." nchar=".$alignment_length.";\n";
	print OTU_DISTINCTS "format datatype=dna interleave=no gap=-;\n";
	print OTU_DISTINCTS "matrix\n";
	
	my $outgroup_id = $nn_id."_Outgroup";
	$outgroup_id =~ s/\|/_id_/g;
	# $outgroup_id =~ s/-/_/g;
	print OTU_DISTINCTS $outgroup_id."\t".$nn_sequence."\n";
	for my $seq_gapped (sort keys %distinct_alleles) {
		if(fast_seq_length($seq_gapped) > $phylogenetic_length_cutoff) {
			my $distinct_id = $unique_sequences{$seq_gapped};
			$distinct_id = filter_one_id($distinct_id);
			$distinct_id =~ s/\|/_id_/g;
			# $distinct_id =~ s/-/_/g;
			print OTU_DISTINCTS $distinct_id."\t";
			print OTU_DISTINCTS $seq_gapped."\n";
		}
	}
	
	print OTU_DISTINCTS ";\n";
	print OTU_DISTINCTS "end;\n";
	
	# my $number_generations = $number_of_phylo_worthy_taxa*100+10000;
	# my $number_generations = 20000000;
	# my $number_generations = 100000;
	my $number_generations = 10000000;
	my $sample_rate = 1000;

	print OTU_DISTINCTS "begin mrbayes;\n";
	print OTU_DISTINCTS "\tset autoclose=yes nowarn=yes;\n";
	print OTU_DISTINCTS "\tlset nst=6 rates=invgamma;\n";
	print OTU_DISTINCTS "\tprset brlenspr=clock:uniform;\n";
	print OTU_DISTINCTS "\tmcmc nruns=1 ngen=".$number_generations." samplefreq=".$sample_rate." savebrlens=yes;\n";
	print OTU_DISTINCTS "\tlog start filename=".$log_filename." replace;\n";
	print OTU_DISTINCTS "\toutgroup ".$outgroup_id.";\n";
	print OTU_DISTINCTS "\tsumt;\n";
	print OTU_DISTINCTS "\tsump;\n";
	print OTU_DISTINCTS "end;\n";
	close(OTU_DISTINCTS);
	
	my $bayes_results = `mrbayes $current_otu_distincts`;
	my $consensus_file = $current_otu_distincts.'.con.tre';
	my $consensus_tree_path = $output_path.$consensus_file;
	my $gmyc_species_list_file = $current_otu_distincts."_species_list.txt";
	my $tree_pdf = $current_otu_distincts."_tree.pdf";
	my $gmyc_results_txt = $current_otu_distincts."_gmyc_results.txt";
	my $gmyc_annotations = $current_otu_distincts."_gmyc_annotations.csv";
	
	unlink $gmyc_species_list_file;
	unlink $output_path.$current_otu_distincts;
	move($current_otu_distincts,$output_path.$current_otu_distincts);
	move($log_filename,$output_path.$log_filename);
	move($consensus_file,$output_path.$consensus_file);
	
	# finished creating tree, move on to coalescent analysis
	$consensus_tree_path =~ s/\\/\\\\/g; # Replace backslashes?
	my $coalescence_script = 'coalesce.r';
	unlink($coalescence_script);
	open(R_COALESCE, '>>'.$coalescence_script);
	# Print part 1 of the coalescence script
	print R_COALESCE	"sink(\"NUL\")\n".
						"require(methods,quietly=TRUE)\n".
						"require(subplex,quietly=TRUE)\n".
						"require(ape,quietly=TRUE)\n".
						"require(splits,quietly=TRUE)\n".
						"require(geiger,quietly=TRUE)\n".
						"source(\"spec.list.1.1.r\")\n".
						"sink()\n".
						"current_tree <- read.nexus(\"$consensus_tree_path\", tree.names = NULL)\n".
						"current_tree <- multi2di(current_tree, random = FALSE)\n".
						"gmyc_results <- gmyc(current_tree, method=\"single\", interval=c(0, 10))\n".
						"sink(\"".$gmyc_results_txt."\")\n".
						"print(summary(gmyc_results))\n".
						"sink()\n".
						"species_list = spec.list(gmyc_results)\n".
						"write.table(species_list,file=\"".$gmyc_species_list_file."\",row.names = FALSE, sep=\",\")\n".
						"tip_colors = c()\n".
						"for(i in 1:length(".'current_tree$tip.label'.")) {\n".
						"	for(j in 1:length(".'species_list$GMYC_spec'.")) {\n".
						"		if(current_tree\$tip.label[i] == ".'toString(species_list$sample_name[j])'.") {\n".
						"			tip_colors = c(tip_colors,".'species_list$GMYC_spec[j])'."\n".
						"			".'current_tree$tip.label[i]'." = paste(\"[\",".'species_list$GMYC_spec[j]'.",\"]_\",".'current_tree$tip.label[i]'.",sep=\"\")\n".
						"			break\n".
						"		}\n".
						"	}\n".
						"}\n".
						"print_to_pdf = TRUE\n".
						"if(print_to_pdf) {\n".
						"	pdffilename = \"".$tree_pdf."\"\n".
						"	pdf(pdffilename,\n".
						"		width = 10, height = 7.5, \n".
						"		pointsize = 12, bg = \"white\")\n".
						"}\n".
						"if(length(current_tree) > 100) {\n".
						"	plot.phylo(	current_tree, \"p\", font=1, cex=0.05,\n".
						"				show.tip.label=TRUE,\n".
						"				tip.color = tip_colors)\n".
						"} else {\n".
						"	plot.phylo(	current_tree, \"p\", font=1, cex=1,\n".
						"				tip.color = tip_colors)\n".
						"}\n".
						"axisPhylo()\n".
						"dev.off()\n";
						
	close(R_COALESCE);
	
	my $R_results = `Rscript $coalescence_script`;
	# print $R_results."\n";
	move($tree_pdf,$output_path.$tree_pdf);

	## Parse GMYC results and extract entity data
	open GMYC_RESULTS, "< ".$gmyc_results_txt;
	my @gmyc_results_lines = <GMYC_RESULTS>;
	my $current_line = 0;
	foreach my $line (@gmyc_results_lines) {
		if ($line =~ m/result of LR test:/) {
			my @delimited_line = split(/ /,$gmyc_results_lines[$current_line]);
			my @split_on_tabs = split(/\t/,$delimited_line[-1]);
			$p_value = $split_on_tabs[-1];
			$p_value =~ s/\n//;
		}
		if ($line =~ m/number of ML clusters:/) {
			my @delimited_line = split(/ /,$gmyc_results_lines[$current_line]);
			my @split_on_tabs = split(/\t/,$delimited_line[-1]);
			$ml_clusters = $split_on_tabs[-1];
			$ml_clusters =~ s/\n//;
			
			my @delimited_line_conf = split(/ /,$gmyc_results_lines[$current_line+1]);
			my @split_on_tabs_conf = split(/\t/,$delimited_line_conf[-1]);
			$ml_clusters_conf = $split_on_tabs_conf[-1];
			$ml_clusters_conf =~ s/\n//;
		}
		if ($line =~ m/number of ML entities:/) {
			my @delimited_line = split(/ /,$gmyc_results_lines[$current_line]);
			my @split_on_tabs = split(/\t/,$delimited_line[-1]);
			$ml_entities = $split_on_tabs[-1];
			$ml_entities =~ s/\n//;
			
			my @delimited_line_conf = split(/ /,$gmyc_results_lines[$current_line+1]);
			my @split_on_tabs_conf = split(/\t/,$delimited_line_conf[-1]);
			$ml_entities_conf = $split_on_tabs_conf[-1];
			$ml_entities_conf =~ s/\n//;
		}
		$current_line++;
	}
	$adjusted_entities = $ml_entities - 1;
	close(GMYC_RESULTS);
	move($gmyc_results_txt,$output_path.$gmyc_results_txt);
	
	##################################################################
	## Parse species list and rebuild haplotypes
	##################################################################
	open SPECIES_LIST, "< ".$gmyc_species_list_file;
	my @species_list_lines = <SPECIES_LIST>;
	close(SPECIES_LIST);
	my %gmyc_hash = ();
	foreach my $species_line (@species_list_lines) {
		my @split_on_comma = split(/,/,$species_line);
		my @split_on_quotes = split(/\"/,$split_on_comma[1]);
		my $gmyc_id = $split_on_quotes[1];
		
		if($gmyc_id =~ m/Outgroup/) { next };
		if($gmyc_id =~ m/sample.name/) { next };
		$species_line =~ s/\n//;
		$gmyc_id =~ s/_id_/\|/g; # Switch the plus sign back to the bar
		# my @split_remainder = split('',$gmyc_id);
		# foreach my $char (@split_remainder) {
			# if ($char eq '_') { $char = '-' };
			# last if $char eq '|';
		# }
		# $gmyc_id = join('',@split_remainder);
		$gmyc_id = filter_one_id($gmyc_id);
		$gmyc_hash{$gmyc_id} = $split_on_comma[0];
	}
	my @gmyc_haplotype_seqs = ();
	my %haplotype_assignment_hash = ();
	for my $gmyc_id (keys %gmyc_hash) {
		my @current_seq = grep { filter_one_id($_->id) eq $gmyc_id } @query_seqs_array;
		# print $gmyc_id."\n";
		# print "id ".$current_seq[0]->id."\n";
		
		# foreach my $id (@query_seqs_array) {
			# if($gmyc_id eq filter_one_id($id->id)) {
				# print "\t".$gmyc_id." => ".filter_one_id($id->id)."\n";
			# } else {
				# print "\t".filter_one_id($id->id)."\n";
			# }
		# }
		
		my $gmyc_haplotype_seq = $current_seq[0];
		foreach my $current_otu_id ( @unique_overall_query_matches ) {
			my $current_otu_seq = $current_otu_seqs_and_id{$current_otu_id};
			my ($transitions,	$transversions,	$num_bases_compared,$current_dist,
				$used_ts_shortcut,
				$variance,		$stderror,		$mink2p,			$maxk2p,$p_stderror,
				$p_min,			$p_max,			$p_dist) = 0;
			my $search_type = 3;
				c_kimura_distance(	$gmyc_haplotype_seq->object_id(),	$current_otu_seq,	$critical_value,
									$cutoff, $search_type, $max_seq_length,
									$transitions,	$transversions,	$num_bases_compared,
									$current_dist,	$variance,		$stderror,
									$mink2p,		$maxk2p);
			if($current_dist == 0) {
				$haplotype_assignment_hash{$current_otu_id} = $gmyc_hash{$gmyc_id};
			}
		}
	}
	for my $specimen_id (keys %haplotype_assignment_hash) {
		my $current_entity_number = $haplotype_assignment_hash{$specimen_id};
		my $filtered_otu_id =  filter_one_id($otu_seq->id);
		my $filtered_specimen_id = filter_one_id($specimen_id);
		if($p_value =~ m/\*/) {
			my $current_gmyc_id = $filtered_otu_id.'_'.$current_entity_number;
			$morpho_name_hash_ref->{$filtered_specimen_id}->{'otu_id'} = $filtered_otu_id;
			$morpho_name_hash_ref->{$filtered_specimen_id}->{'gmyc_id'} = $current_gmyc_id;
			$gmyc_name_id_hash{$filtered_specimen_id} = $current_gmyc_id;
			print OTU_GMYC_RESULTS $filtered_otu_id.','.$current_gmyc_id.','.$filtered_specimen_id.','.$p_value."\n";
			$num_specimens_assigned_gmyc++;
		} else {
			$morpho_name_hash_ref->{$filtered_specimen_id}->{'otu_id'} = $filtered_otu_id;
			$morpho_name_hash_ref->{$filtered_specimen_id}->{'gmyc_id'} = $filtered_otu_id;
			$gmyc_name_id_hash{$filtered_specimen_id} = $filtered_otu_id;
			print OTU_GMYC_RESULTS $filtered_otu_id.','.$filtered_otu_id.','.$filtered_specimen_id.','.$p_value."\n";
			$num_specimens_assigned_gmyc++;
		}
	}
	move($gmyc_species_list_file,$output_path.$gmyc_species_list_file);
	##################################################################
	# End GMYC species list
	# End main GMYC code
	##################################################################
	#*****************************************************************
	##################################################################
	my @removeable_files = ();
	push(@removeable_files,<*$current_otu_distincts*>);
	foreach my $removeable_file (@removeable_files) {
		unlink($removeable_file);
	}
	### Parse GMYC results and rebuild haplotypes
	skip_phylo: # Skip phylogenetic analysis if less than 3 taxa
	if($skipped_phylo == 1) {
		foreach my $query_seq_id (@unique_overall_query_matches) {
			my @query_seq = grep { $_->id eq $query_seq_id } @query_seqs_array;
			my $filtered_otu_id =  filter_one_id($otu_seq->id);
			my $filtered_specimen_id = filter_one_id($query_seq[0]->id());
			$morpho_name_hash_ref->{$filtered_specimen_id}->{'otu_id'} = $filtered_otu_id;
			$morpho_name_hash_ref->{$filtered_specimen_id}->{'gmyc_id'} = $filtered_otu_id;
			$gmyc_name_id_hash{$filtered_specimen_id} = $filtered_otu_id;
			print OTU_GMYC_RESULTS $filtered_otu_id.','.$filtered_otu_id.','.$filtered_specimen_id.',ns'."\n";
		}
	}
	##################################################################
	# Stat summaries
	my $mean_sequence_length 	= 0; # Sequence length stats
	my $min_sequence_length 	= 0;
	my $max_sequence_length 	= 0;
	my $se_sequence_length		= 0;
	my $mean_number_comparisons	= 0; # Num comparisons stats
	my $min_number_comparisons 	= 0;
	my $max_number_comparisons 	= 0;
	my $se_number_comparisons	= 0;
	my $mean_distance 			= 0; # Normal distance stats
	my $min_distance 			= 0;
	my $max_distance			= 0;
	my $se_distance				= 0;
	my $mean_stderrors			= 0; # Standard errors stats
	my $min_stderrors			= 0;
	my $max_stderrors			= 0;
	my $se_stderrors			= 0;
	my $percent_dominant_allele	= 100; # Allele stats
	my $num_unique_alleles		= 0;
	my $num_distinct_alleles	= 0;
	
	my $seq_rounding 			= "%.1f";
	my $seq_length_stat 		= Statistics::Descriptive::Full->new();
	my $num_comparisons_stat 	= Statistics::Descriptive::Full->new();
	if(scalar(@seq_lengths) > 0) {
		$seq_length_stat->add_data(@seq_lengths);
		$mean_sequence_length 	= sprintf($seq_rounding,$seq_length_stat->mean());
		$min_sequence_length 	= sprintf($seq_rounding,$seq_length_stat->min());
		$max_sequence_length 	= sprintf($seq_rounding,$seq_length_stat->max());
		$se_sequence_length		= sprintf($seq_rounding,($seq_length_stat->standard_deviation()/sqrt(scalar(@seq_lengths))));
	}
	if(scalar(@num_comparisons)) {
		$num_comparisons_stat->add_data(@num_comparisons);
		$mean_number_comparisons 	= sprintf($seq_rounding,$num_comparisons_stat->mean());
		$min_number_comparisons 	= sprintf($seq_rounding,$num_comparisons_stat->min());
		$max_number_comparisons 	= sprintf($seq_rounding,$num_comparisons_stat->max());
		$se_number_comparisons		= sprintf($seq_rounding,($num_comparisons_stat->standard_deviation()/sqrt(scalar(@num_comparisons))));
	}
	my $rounding 					= "%.3f";
	my $distances_stat 				= Statistics::Descriptive::Full->new();
	my $stderrors_stat 				= Statistics::Descriptive::Full->new();
	my $deflated_distances_stat 	= Statistics::Descriptive::Full->new();
	my $deflated_std_errors_stat	= Statistics::Descriptive::Full->new();
	my $allele_stats			 	= Statistics::Descriptive::Full->new();
	if(scalar(@distances) > 0) {
		$distances_stat->add_data(@distances);
		$mean_distance 			= sprintf($rounding,$distances_stat->mean())*100;
		$min_distance 			= sprintf($rounding,$distances_stat->min())*100;
		$max_distance 			= sprintf($rounding,$distances_stat->max())*100;
		$se_distance			= sprintf($rounding,($distances_stat->standard_deviation()/sqrt(scalar(@distances))))*100;
		
		$stderrors_stat->add_data(@std_errors);
		$mean_stderrors 		= sprintf($rounding,$stderrors_stat->mean())*100;
		$min_stderrors 			= sprintf($rounding,$stderrors_stat->min())*100;
		$max_stderrors 			= sprintf($rounding,$stderrors_stat->max())*100;
		$se_stderrors			= sprintf($rounding,($stderrors_stat->standard_deviation()/sqrt(scalar(@std_errors))))*100;
	}
	##################################################################
	
	##################################################################
	# Calculate the percentage of query specimens accoutned for in each link step
	my $abundance = scalar(@unique_overall_query_matches);
	my $link_strength_string = '';
	my $link_strength_i = 1;
	my $last_link = scalar(@link_strength);
	my $link_rounding = "%.0f";
	foreach my $link (@link_strength) {
		last if ($link_strength_i == $last_link);
		my $current_link_strength = sprintf($link_rounding,($link/$abundance*100));
		if($link_strength_i != ($last_link-1)) {
			$link_strength_string .= $current_link_strength."_";
		} else {
			$link_strength_string .= $current_link_strength;
		}
		$link_strength_i++;
	}
	##################################################################
	
	##################################################################
	# Number of morpho names
	my %unique_overall_names = map {$_,1} @overall_names;
	my @unique_overall_names = keys %unique_overall_names;
	$num_unique_alleles = keys %unique_alleles;
	# $num_distinct_alleles = keys %unique_distinct_alleles;
	$num_distinct_alleles = keys %distinct_alleles;
	##################################################################
	
	##################################################################
	# Console output
	printf	$string_space,
			"[".$otu_i."] ".filter_one_id($otu_seq->id),
			"[".scalar(@unique_overall_names)."]",
			"[".$abundance." / ".$num_unique_alleles." / ".$num_distinct_alleles."]",
			# "".$num_specimens_assigned_gmyc,
			"[".$p_value." / ".$ml_clusters." / ".$ml_clusters_conf." / ".$ml_entities." / ".$ml_entities_conf." / ".$adjusted_entities."]",
			"".$mean_distance."% SE: ".$se_distance."% (".$min_distance."% - ".$max_distance."%) ",
			"".$mean_number_comparisons." SE: ".$se_number_comparisons." (".$min_number_comparisons." - ".$max_number_comparisons.")",
			"".$mean_sequence_length." SE: ".$se_sequence_length." (".$min_sequence_length." - ".$max_sequence_length.")",
			"".$mean_stderrors."% SE: ".$se_stderrors."% (".$min_stderrors."% - ".$max_stderrors."%) ",
			"[".scalar(@unique_otu_links)."]",
			"[".$link_depth."]",
			"[".$link_strength_string."]",
			"".$nn_id." -> ".$nn_dist."% (".$nn_min."% - ".$nn_max."%)";
	print OTU_SUMMARY
			"[".$otu_i."] ".filter_one_id($otu_seq->id).$delimiter.
			"[".scalar(@unique_overall_names)."]".$delimiter.
			"[".$abundance." / ".$num_unique_alleles." / ".$num_distinct_alleles."]".$delimiter.
			"[".$p_value." / ".$ml_clusters." / ".$ml_clusters_conf." / ".$ml_entities." / ".$ml_entities_conf." / ".$adjusted_entities."]".$delimiter.
			"".$mean_distance."% SE: ".$se_distance."% (".$min_distance."% - ".$max_distance."%) ".$delimiter.
			"".$mean_number_comparisons." SE: ".$se_number_comparisons." (".$min_number_comparisons." - ".$max_number_comparisons.")".$delimiter.
			"".$mean_sequence_length." SE: ".$se_sequence_length." (".$min_sequence_length." - ".$max_sequence_length.")".$delimiter.
			"".$mean_stderrors."% SE: ".$se_stderrors."% (".$min_stderrors."% - ".$max_stderrors."%) ".$delimiter.
			"[".scalar(@unique_otu_links)."]".$delimiter.
			"[".$link_depth."]".$delimiter.
			"[".$link_strength_string."]".$delimiter.
			"".$nn_id." -> ".$nn_dist."% (".$nn_min."% - ".$nn_max."%)".$delimiter."\n";
	print OTU_EXCEL
			$otu_i.$delimiter.filter_one_id($otu_seq->id).$delimiter.
			scalar(@unique_overall_names).$delimiter.
			$abundance.$delimiter.$num_unique_alleles.$delimiter.$num_distinct_alleles.$delimiter.
			$p_value.$delimiter.$ml_clusters.$delimiter.$ml_clusters_conf.$delimiter.$ml_entities.$delimiter.$ml_entities_conf.$delimiter.$adjusted_entities.$delimiter.
			$mean_distance.$delimiter.$se_distance.$delimiter.$min_distance.$delimiter.$max_distance.$delimiter.
			$mean_number_comparisons.$delimiter.$se_number_comparisons.$delimiter.$min_number_comparisons.$delimiter.$max_number_comparisons.$delimiter.
			$mean_sequence_length.$delimiter.$se_sequence_length.$delimiter.$min_sequence_length.$delimiter.$max_sequence_length.$delimiter.
			$mean_stderrors.$delimiter.$se_stderrors.$delimiter.$min_stderrors.$delimiter.$max_stderrors.$delimiter.
			scalar(@unique_otu_links).$delimiter.
			$link_depth.$delimiter.
			$link_strength_string.$delimiter.
			$nn_id.$delimiter.$nn_dist.$delimiter.$nn_min.$delimiter.$nn_max.$delimiter."\n";

	##################################################################
	foreach my $unique_overall_name (@unique_overall_names) {
		my $unique_name_abundance = 0;
		foreach my $unique_name (@overall_names) {
			if($unique_overall_name eq $unique_name) {
				$unique_name_abundance++;
			}
		}
		print "\t->".$unique_overall_name." [".$unique_name_abundance."]\n";
		print OTU_SUMMARY "         -->".$unique_overall_name." [".$unique_name_abundance."]".$delimiter."\n";
		# print "\tp-value: ".$p_value."\n\t# clusters: ".$ml_clusters."\n\tcluster range: ".$ml_clusters_conf."\n\t# entities: ".$ml_entities."\n\tentities range: ".$ml_entities_conf."\n";
	}
	push(@otu_morpho_lumps,scalar(@unique_overall_names));
	##################################################################
	my %unique_gmyc_otu_ids = reverse %gmyc_name_id_hash;
	my @current_gmyc_names = ();
	my $num_gmyc_otu_ids = 1;
	for my $gmyc_otu_id ( sort keys %unique_gmyc_otu_ids ) {
		for my $specimen_id ( sort keys %gmyc_name_id_hash ) {
			if($gmyc_name_id_hash{$specimen_id} eq $gmyc_otu_id) {
				push(@current_gmyc_names, convert_id_to_name($specimen_id));
			}
		}
		my %unique_gmyc_names = map {$_,1} @current_gmyc_names;
		my @unique_gmyc_names = keys %unique_gmyc_names;
		print "\t --> [".$num_gmyc_otu_ids."] ".$gmyc_otu_id."\n";
		foreach my $name (@unique_gmyc_names) {
			print "\t\t".$name."\n";
		}
		push(@gmyc_morpho_lumps,scalar(@unique_gmyc_names));
		$num_gmyc_otu_ids++;
	}
	##################################################################
	push(@query_seqs_found,@unique_overall_query_matches);
	$otu_i++;
}
close(EXEMPLARS);
close(OTU_GMYC_RESULTS);
print OTU_SUMMARY "\n";
print OTU_SUMMARY "Found ".scalar(@query_seqs_found)."\n";
print OTU_SUMMARY "Not matched :\n";
print "\n";
print "Found ".scalar(@query_seqs_found)."\n";
print "Not matched :\n";
foreach my $orig_query_seq (@original_sequence_array) {
	if ($orig_query_seq->id ~~ @query_seqs_found) {
		next;
	} else {
		print OTU_SUMMARY "\t".$orig_query_seq->id." Sequence Length: ".fast_seq_length($orig_query_seq->seq())."\n";
		print "\t".$orig_query_seq->id." Sequence Length: ".fast_seq_length($orig_query_seq->seq())."\n";
	}
}
print "Exemplars printed: ".$num_exemplars."\n";

print "\n\n";
print "Identification Sucess:\n\n";
## Identification analysis
my %id_analysis_morpho_ids = ();
my %id_analysis_otu_ids = ();
my %id_analysis_gmyc_ids = ();
# Build keys for each type of ID (morpho, otu, gmyc)
for my $specimen_id (sort keys %$morpho_name_hash_ref) {
	$id_analysis_morpho_ids{convert_id_to_name($specimen_id)} = 1;
	$id_analysis_otu_ids{$morpho_name_hash_ref->{$specimen_id}->{'otu_id'}} = 1;
	$id_analysis_gmyc_ids{$morpho_name_hash_ref->{$specimen_id}->{'gmyc_id'}} = 1;
}
my @otu_morpho_correspondence = ();
my @gmyc_morpho_correspondence = ();
for my $morpho_id (sort keys %id_analysis_morpho_ids) {
	my @current_morpho_otu_ids = ();
	my @current_morpho_gmyc_ids = ();
	for my $specimen_id (sort keys %$morpho_name_hash_ref) {
		my $current_morpho_name = convert_id_to_name($specimen_id);
		if ($current_morpho_name eq $morpho_id) {
			push(@current_morpho_otu_ids, $morpho_name_hash_ref->{$specimen_id}->{'otu_id'});
			push(@current_morpho_gmyc_ids, $morpho_name_hash_ref->{$specimen_id}->{'gmyc_id'});
		}
	}
	my %current_unique_otu_ids = map {$_,1} @current_morpho_otu_ids;
	my %current_unique_gmyc_ids = map {$_,1} @current_morpho_gmyc_ids;
	my $number_otu_ids = keys %current_unique_otu_ids;
	my $number_gmyc_ids = keys %current_unique_gmyc_ids;
	push(@otu_morpho_correspondence,$number_otu_ids);
	push(@gmyc_morpho_correspondence,$number_gmyc_ids);
	# print $morpho_id."\n";
	# print "\tOTU: \n";
	# for my $unique_otu (sort keys %current_unique_otu_ids) {
		# print "\t\t".$unique_otu."\n";
	# }
	# print "\tGMYC: \n";
	# for my $unique_gmyc (sort keys %current_unique_gmyc_ids) {
		# print "\t\t".$unique_gmyc."\n";
	# }
}
my %unique_otu_correspondences = map {$_,1} @otu_morpho_correspondence;
my %unique_gmyc_correspondences = map {$_,1} @gmyc_morpho_correspondence;
my $identification_rounding = "%.1f";
print "Morpho Correspondence Success rates: \n";
print "\tOTU: # IDs, % Correspondence\n";
for my $unique_otu_correspondences (sort keys %unique_otu_correspondences) {
	my @current_correspondences = grep { $_ == $unique_otu_correspondences } @otu_morpho_correspondence;
	my $percent_correspondence = scalar(@current_correspondences)/scalar(@otu_morpho_correspondence)*100;
	$percent_correspondence = sprintf($identification_rounding,$percent_correspondence);
	print "\t\t".$unique_otu_correspondences.") ".$percent_correspondence."\n";
}
print "\n";
print "\tGMYC: # IDs, % Correspondence\n";
for my $unique_gmyc_correspondences (sort keys %unique_gmyc_correspondences) {
	my @current_correspondences = grep { $_ == $unique_gmyc_correspondences } @gmyc_morpho_correspondence;
	my $percent_correspondence = scalar(@current_correspondences)/scalar(@gmyc_morpho_correspondence)*100;
	$percent_correspondence = sprintf($identification_rounding,$percent_correspondence);
	print "\t\t".$unique_gmyc_correspondences.") ".$percent_correspondence."\n";
}

my %unique_otu_morpho_lumps = map {$_,1} @otu_morpho_lumps;
my %unique_gmyc_morpho_lumps = map {$_,1} @gmyc_morpho_lumps;
print "\tOTU Lumping Rate: # Morpho, % Lumping\n";
for my $unique_otu_morpho_lumps (sort keys %unique_otu_morpho_lumps) {
	my @current_correspondences = grep { $_ == $unique_otu_morpho_lumps } @otu_morpho_lumps;
	my $percent_correspondence = scalar(@current_correspondences)/scalar(@otu_morpho_lumps)*100;
	$percent_correspondence = sprintf($identification_rounding,$percent_correspondence);
	print "\t\t".$unique_otu_morpho_lumps.") ".$percent_correspondence."\n";
}
print "\tGMYC Lumping Rate: # Morpho, % Lumping\n";
for my $unique_gmyc_morpho_lumps (sort keys %unique_gmyc_morpho_lumps) {
	my @current_correspondences = grep { $_ == $unique_gmyc_morpho_lumps } @gmyc_morpho_lumps;
	my $percent_correspondence = scalar(@current_correspondences)/scalar(@gmyc_morpho_lumps)*100;
	$percent_correspondence = sprintf($identification_rounding,$percent_correspondence);
	print "\t\t".$unique_gmyc_morpho_lumps.") ".$percent_correspondence."\n";
}

my $t1 = Benchmark->new;
my $time = timediff($t1, $t0);
print "\n";
print timestr($time)."\n";

print OTU_SUMMARY "\n";
print OTU_SUMMARY timestr($time).$delimiter."\n";

print OTU_SUMMARY "[OTU #] OTU Ref. ID:			Reference ID\n";
print OTU_SUMMARY "[# Morpho]:					# of morphological identifications\n";
print OTU_SUMMARY "[Abundance:					Alleles]	# of sequences: # of unique alleles\n";
print OTU_SUMMARY "Avg. Dist: Max (Adjusted):	# Average K2P Maximum K2P (Minimum Possible K2P)\n";
print OTU_SUMMARY "Upper Limit:					Upper limit of the widest confidence interval\n";
print OTU_SUMMARY "Avg. Length: (Min - Max):	Average sequence length (Minimum - Maximum)\n";
print OTU_SUMMARY "[Sub-groups]:				# of non-statistical sub-groups\n";
print OTU_SUMMARY "[Link Depth]:				Network depth that sub-groups are connected by\n";
print OTU_SUMMARY "[Link Strength %]:			Distribution of sequences within each sub-group\n";

close(OTU_SUMMARY);

########################################################################################################
## Begin subs																						   #
########################################################################################################

sub find_query_matches {
	# Searches a query result hash and returns and array of subsequent 
	# Query sequences for the current OTU.
	my $query_results_hash_ref = shift;
	my $otu_seq_id = shift;
	my %query_results_hash = %$query_results_hash_ref;
	my @query_matches = ();
	while (my ($query_id,$matched_otu_array) = each (%query_results_hash)) { # For each query seq
		foreach my $matched_otu_id (@$matched_otu_array) { # For each matched OTU
			if($matched_otu_id eq $otu_seq_id) {
				push(@query_matches,$query_id);
			}
		}
	}
	
	return \@query_matches;
}

sub tag_array_ids_otu {
	my $id_array_ref = shift;
	my $current_otu_number = shift;
	
	my $num_ids = scalar(@$id_array_ref);
	for (my $i = 1; $i < $num_ids; $i++) {
		@$id_array_ref[$i] .= $current_otu_number;
	}
	
	return $id_array_ref;
}

sub fast_seq_length {
	my $seq = shift;
	
	$seq =~ s/-/ /g;
	$seq =~ s/\s+//g;
	
	return length($seq);
}


# randomly permutate @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i = @$array;
    while ( --$i ) {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}

#################
## General subs
#################

sub fix_bold_fasta {
	my ( $f ) = shift;
	open F, "< $f";
	my @fasta = <F>;
	close F or die "Couldn't close $f!";
	for(my $i = 0; $i <= 15; $i++) {
		foreach my $line (@fasta) {
			if($line =~ /^>/) {
				$line =~ s/ /_/; # replace whitespace with _
				$line =~ s/-/_/; # replace whitespace with _
				$line =~ s/BOLD:/BOLD_/; # replace whitespace with _
			}
		}
	}
		
	unlink ($f, 0) or die "Could not overwrite $f.";
	open (MYFILE, '>>'.$f) or die "Could not open $f";
	foreach my $line(@fasta) {
		print MYFILE $line;
	}
	close(MYFILE);
}

sub filter_one_id {
	my $id = shift;
	my $filtered_id = '';

	my @delimited_id = split(/\|/,$id);

	my $num_delimiters = 0;
	$num_delimiters = scalar(@delimited_id);

	# IF there are 5 bold delimited sections, here
	if($num_delimiters == 2) {
		my $new_id = $delimited_id[0]."|".$delimited_id[1];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	}
	elsif($num_delimiters == 5) {
		my $new_id = $delimited_id[2]."|".$delimited_id[3];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	# Only 4 bold delimited sections
	} elsif(($num_delimiters == 4) && ($delimited_id[3] eq "COI_5P")) {
		my $new_id = $delimited_id[1]."|".$delimited_id[2];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	} elsif(($num_delimiters) == 4 && ($delimited_id[2] ne "COI_5P")) {
		my $new_id = $delimited_id[2]."|".$delimited_id[3];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	} 
	elsif($num_delimiters == 3) {
		my $new_id = $delimited_id[1]."|".$delimited_id[2];
		$new_id =~ s/ /_/;
		$filtered_id = $new_id;
	}
	else {
		$filtered_id = $id;
	}
	return $filtered_id;
}

sub convert_id_to_name {
	# Returns a species id from an already filtered bold record
	# Input a node string id
	# Input: ID#|Name
	# Output: Name
	my ($node_id) = @_;
	my @delimited_id = ();
	@delimited_id = split(/\|/,$node_id);
	if(scalar(@delimited_id) == 2) {
		my $new_id = $delimited_id[1];
		return $new_id;
	} else {
		return $node_id;
	}
}


