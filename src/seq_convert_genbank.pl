#!/usr/bin/env perl
# Downloads and parses genbank files given an input taxon name.
#
# Copyright (c) 2013-2014 Bryan White, bpcwhite@gmail.com

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
use Time::HiRes qw( usleep );
use General::Arguments;
use Sequence::Fasta;

use Bio::Align::AlignI;
use Bio::AlignIO::fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqIO::fasta;
use strict;
use warnings;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Bio::Tools::EUtilities::Query;

use Bio::SeqIO::genbank;
use Bio::Seq::RichSeq;
use Bio::DB::Taxonomy;
use Bio::DB::Taxonomy::list;

use HTTP::Response;

use Data::Dumper;
use Devel::Size qw(size total_size);
use Test::LeakTrace;
use Devel::Size;
use Text::Fuzzy::PP;
use Benchmark qw(:all);

##################################################################
# Start benchmark
my $t0 = Benchmark->new;
my $k2p1 = 0;
##################################################################


print "
******************************************************************

Copyright (c) 2013-2015 Bryan White, bpcwhite\@gmail.com

GNU General Public License, Version 3, 29 June 2007

This program comes with ABSOLUTELY NO WARRANTY; for details type.
This is free software, and you are welcome to redistribute it
under certain conditions.

A full copy of the GPL 3.0 license should be accompanied with any
distribution of this software.

******************************************************************
\n\n
";

my $params = General::Arguments->new(	arguments_v => \@ARGV,
									option_defs => {'-list' 		=> '', 				# List file name
													'-slim' 		=> 999999999,		# Sequence search limit
													'-sflim'		=> 999999999,		# Sequence fetch limit
													'-batch-cap'	=> 350,				# NCBI sequence fetch limit per batch
													'-tlim' 		=> 999999999,		# Taxa limit
													'-user_email' 	=> 'foo@bar.com', 	# User email
													'-outp' 		=> 'output', 		# Output file prefix
													'-query'		=> '',				# Target gene, etc.
													'-pubmed'		=> 0,				# By default skip pubmed download
													'-count-seqs'	=> 0,				# Toggle to only count sequences
													'-voucher-only' => 0,				# Toggle using voucher only sequences
													'-term'			=> '',				# Use a search term instead of a list
													'-cq-file'		=> '',				# Filename for a custom query list
													'-db'			=> 'nucleotide',	# Change to a different ncbi database
													'-pre-dl'		=> '0',				# Use pre-downloaded genbank files
													'-slow-download'=> '0',				# Slow batch download
													'-exemplar-only'=> '0',				# Only get 1 sequence per taxa name
													'-match'		=> '',				# Match against a query sequence
													'-batch-output'	=> '', # Output without headers for merging large lists of taxa
													'-std-out' => 0, # print to std out instead of output files
													}
													);
# Initiate parameters

my $taxa_file = $params->options->{'-list'};
my $ncbi_database = $params->options->{'-db'};
my $pre_downloaded = $params->options->{'-pre-dl'};
my $slow_download = $params->options->{'-slow-download'};
my $exemplar_only = $params->options->{'-exemplar-only'};
my $match_aln_file = $params->options->{'-match'};
my $batch_output = $params->options->{'-batch-output'};
my $std_outp = $params->options->{'-std-out'};

my @taxa_list = ();
my $seq_length_maximum = 3000;

my $match_seq = '';
if($match_aln_file ne '') {
	if ($std_outp != 1) {
		print "Importing alignment file ".$match_aln_file."...\n";
	}
	my $match_aln = Bio::AlignIO->new(-format => 'fasta',
									-file   => $match_aln_file );
	my $match_aln_obj = $match_aln->next_aln;
	my @match_seqs;
	foreach my $seq ($match_aln_obj->each_seq) {
		push(@match_seqs,$seq->seq);
	}
	$match_seq = $match_seqs[0];
}


if ($params->options->{'-term'}) {
	push(@taxa_list, $params->options->{'-term'});
	# $params->options->{'-outp'} = $params->options->{'-term'};
} else {
	open (TAXALIST, '<'.$taxa_file);
	@taxa_list = <TAXALIST>;
	close(TAXALIST);
}
##############################################################################
my @overall_results = ();
my $taxa_counter = 0;
my $overall_output_file = $params->options->{'-outp'}.'.csv';
my $endl = "\n";
my $suppress_output = 'yes';
foreach my $taxa (@taxa_list) {
	$taxa =~ s/\n//g; # replace newlines
	$taxa =~ s/ /_/g; # replace whitespace

	my ($results) = download_target_taxa($taxa,$taxa_counter,$suppress_output,$params);
	push(@overall_results, @$results);

	$taxa_counter++;
}
# Cleaning up GB files...
foreach my $taxa (@taxa_list) {
	# $taxa =~ s/\n//g; # replace newlines
	my $file = 'seqs_'.$taxa.'.gb';
	unlink $file;
}
##############################################################################

##############################################################################
# Print successful output
unlink $overall_output_file;
if($std_outp == 1) {
	foreach my $output_line (@overall_results) {
		print "\"".	$output_line if $output_line ne $endl;
		print $output_line if $output_line eq $endl;
	}
} else {
	open (OVEROUTPUT, '>>'.$overall_output_file);
	foreach my $output_line (@overall_results) {
		print OVEROUTPUT "\"".	$output_line if $output_line ne $endl;
		print OVEROUTPUT 		$output_line if $output_line eq $endl;
	}
	close(OVEROUTPUT);
}
##############################################################################

sub download_target_taxa {
	##############################################################################
	# Subroutine parameters
	my $target_taxon = shift;
	my $taxa_counter = shift;
	my $suppress_output = shift;
	my $params = shift;
	##############################################################################

	##############################################################################
	my $sequence_limit = $params->options->{'-slim'};
	my $skip_pubmed_search = $params->options->{'-pubmed'};
	my $voucher_only = $params->options->{'-voucher-only'};
	# my $skip_pubmed_search = 1;
	my $number_seqs_found = 'NA';

	my $search_options = search_strings($params->options->{'-query'},
										$params->options->{'-cq-file'});
	my $taxa_failed = 0; # Flag for if the taxa lookup fails
	my $failed_taxa = 'NA';
	my $sequence_failed = 0; # Flag for if the sequence lookup fails
	my $failed_sequence = 'NA';
	# my $search_options = '';
	my $taxon_limit = 1;
	my $user_email = 'blah@blah.com';
	my $dlm = ',';
	my $endl = "\n";
	my $output_file = $target_taxon."_output.csv";
	my $sleep_time = 1000; # microsends, pause between seconds
	my $max_num_tries = 2;
	my $max_pubmed_tries = 2;
	my $sequence_file = '';

	##############################################################################
	# Optimization hashes. Swap memory for speed.
	# Stores the info for a given pubmed article to save download time.
	my %pubmed_hash = ();
	my $pubmed_hash_ref = \%pubmed_hash;
	# Stores the failed sequence searches
	my %failed_sequence_hash = ();
	my $failed_sequence_hash_ref = \%failed_sequence_hash;
	# Stores the failed taxa searches.
	my %failed_search_hash = ();
	my $failed_search_hash_ref = \%failed_search_hash;
	##############################################################################

	##############################################################################
	my @output_header = ('taxon_query',$dlm,
		'taxon_id',$dlm,
		'accession',$dlm,
		'seq_query',$dlm,
		'seqs_found',$dlm,
		'description',$dlm,
		'gene',$dlm,
		'product',$dlm,
		'binom',$dlm,
		'tax_hierarchy',$dlm,
		'pub_title',$dlm,
		'pub_authors',$dlm,
		'pub_abstract_text',$dlm,
		'jrn_name',$dlm,
		'jrn_doi',$dlm,
		'jrn_so',$dlm,
		'jrn_volume',$dlm,
		'jrn_issue',$dlm,
		'jrn_pages',$dlm,
		'jrn_pubdate',$dlm,
		'nuc_seq',$dlm,
		'fasta_nuc',$dlm,
		'aln_seq',$dlm,
		'aln_score',$dlm,
		'prot_seq',$dlm,
		'primers',$dlm,
		'codon_start',$dlm,
		'collection_date',$dlm,
		'voucher_id',$dlm,
		'collected_by',$dlm,
		'identified_by',$dlm,
		'organelle',$dlm,
		'location',$dlm,
		'lat_lon',$dlm,
		'isolate',$dlm,
		'bio_material',$dlm,
		$endl);
		# Output file headers.
	my @output_lines = ();
	my @return_output_lines = ();
	if($batch_output ne 1) {
		push(@output_lines, @output_header);
	}
	push(@return_output_lines, @output_lines) if($taxa_counter == 0);
	##############################################################################
	# Retrieve taxon ID

	print "[".$taxa_counter."] Searching for $target_taxon\n";
	my @taxon_ids = ();
	print "\tSearching taxonomy\n";
	my $taxonomy_eutil_tries = 1;

	taxonomy_eutil:
	my $taxonomy_eutil = Bio::DB::EUtilities->new(-eutil    => 'esearch',
												   -db      => 'taxonomy',
												   -retmax  => $taxon_limit,
												   -rettype => 'gb',
												   -email   => $user_email,
												   -term    => $target_taxon);
	eval { @taxon_ids = $taxonomy_eutil->get_ids(); };
	if ($@ || (scalar(@taxon_ids) == 0)) {
		print "\tProblem in taxonomy_eutil. Retrying...$taxonomy_eutil_tries\n";
		if ($taxonomy_eutil_tries == $max_num_tries) {
			$failed_search_hash_ref->{$target_taxon}->{'taxa_id'} = 'NA';
			goto taxa_failed;
		}
		$taxonomy_eutil_tries++;
		usleep($sleep_time); 	# Sleep so you don't overload NCBI's servers.
		goto taxonomy_eutil;
	}

	my $taxon_id = $taxon_ids[0];
	print "\tFound taxon ID: $taxon_id\n";
	##############################################################################

	##############################################################################
	# Retrieve nucleotide sequences
	if($pre_downloaded == 1) {
		$sequence_file = 'seqs_'.$target_taxon.'.gb';
		# unlink($sequence_file);
		goto use_pre_downloaded;
	} else {
		$sequence_file = 'seqs_'.$target_taxon.'.gb';
		unlink($sequence_file);
	}
	my $sequence_search_tries = 1;
	sequence_search:
	print "\tRetrieving sequences...\n";
	my $sequence_term = 'txid'.$taxon_id.'[Organism:exp]';
	print $sequence_term."\n";
	my $sequence_search_limit = ($params->options->{'-slim'});
	$sequence_search_limit = 999999999 if ($params->options->{'-count-seqs'} == 1);

	my $sequence_search = Bio::DB::EUtilities->new(-eutil    => 'esearch',
												   -db      => $ncbi_database,
												   -retmax  => $sequence_search_limit,
												   -rettype => 'gb',
												   -email   => $user_email,
												   -term    => $sequence_term.' '.$search_options);

	my @sequence_ids = ();
	eval { @sequence_ids = $sequence_search->get_ids(); };
	if ($@ || (scalar(@sequence_ids) == 0)) {
		print "\tProblem in sequence_search. Retrying...\n";
		if ($sequence_search_tries == $max_num_tries) {
			$failed_search_hash_ref->{$target_taxon}->{'taxa_id'} = $taxon_id;
			goto sequence_failed;
		}
		$sequence_search_tries++;
		usleep($sleep_time); 	# Sleep so you don't overload NCBI's servers.
		goto sequence_search;
	}
	my %unique_seq_ids = ();
	foreach my $seq_id (@sequence_ids) {
		# print $seq_id."\n";
		if(exists($unique_seq_ids{$seq_id})) {
			$unique_seq_ids{$seq_id}++;
		} else {
			$unique_seq_ids{$seq_id} = 1;
		}
	}
	@sequence_ids = keys %unique_seq_ids;
	$number_seqs_found = scalar @sequence_ids;
	# If only looking to count sequences, go here.
	if ($params->options->{'-count-seqs'} == 1) {
		$failed_search_hash_ref->{$target_taxon}->{'taxa_id'} = $taxon_id;
		usleep($sleep_time); 	# Sleep so you don't overload NCBI's servers.
		goto just_count_seqs;
	}
	print "\tFound ".$number_seqs_found." sequences to download.\n";
	##############################################################################
	# Batch downloader
	my $fetch_limit = $params->options->{'-sflim'}; # Fetch a limited number of seqs
	my $sequence_download_tries = 1;
	my $ncbi_sequence_cap = $params->options->{'-batch-cap'};
	my $int_number_batches = int @sequence_ids / $ncbi_sequence_cap;
	my $mod_number_batches = @sequence_ids % $ncbi_sequence_cap;
	$int_number_batches++ if($mod_number_batches > 0);
	print "\tSequences will be downloaded in ".$int_number_batches." batches\n";
	my $total_batched = 0;
	for (my $batch_i = 0; $batch_i <= $int_number_batches; $batch_i++) {
		my @sub_batch_list = ();
		for (my $sub_i = 0; $sub_i <= $ncbi_sequence_cap; $sub_i++) {
			last if $total_batched == scalar @sequence_ids;
			last if $total_batched == $fetch_limit;
			push(@sub_batch_list, $sequence_ids[$total_batched]);
			$total_batched++;
		}

		sequence_download:
		my $starting_count = 0;
		if ($batch_i == 0) {
			$starting_count = 0;
		} else {
			$starting_count = $total_batched-$ncbi_sequence_cap;
		}
		print "\tDownloading sequences...[".$batch_i."] ".$starting_count." to ".$total_batched."\n";
		my $sub_batch_file = $batch_i."_".$sequence_file;
		my $sequence_fetch = '';
		if($slow_download != 0) {
			usleep($slow_download);
		}
		$sequence_fetch = Bio::DB::EUtilities->new( -eutil   	=> 'efetch',
													   -db      => 'nucleotide',
													   -rettype => 'gb',
													   -email   => $user_email,
													   -id      => \@sub_batch_list);

		eval {
			$sequence_fetch->get_Response(-file => $sub_batch_file);
		};
		if($@) {
			print "\tProblem in sequence download. Retrying...\n";
			print "\t".$@."\n";
			if ($sequence_download_tries == $max_num_tries) {
				$failed_search_hash_ref->{$target_taxon}->{'taxa_id'} = $taxon_id;
				goto sequence_failed;
			}
			$sequence_download_tries++;
			usleep($sleep_time); 	# Sleep so you don't overload NCBI's servers.
			goto sequence_download;
		}
		last if $total_batched == scalar @sequence_ids;
		last if $total_batched == $fetch_limit;
	}

	print "\tSequences downloaded to genbank format.\n";
	##############################################################################

	##############################################################################
	# use_pre_downloaded:
	# print $sequence_file."\n";
	my @batched_files = <*$sequence_file>;
	my @concatenated_batch_file = ();
	foreach my $file (@batched_files) {
		# print $file."\n";
		# exit;
		next if $file eq $sequence_file;
		# print "B";
		open (CURRENT, "<$file") or die $!;
		my @current_file = <CURRENT>;
		# foreach my $current (@current_file) {
			# print $current."\n";
			# exit;
		# }
		push(@concatenated_batch_file, @current_file);
		close CURRENT or die $!;
	}

	# exit;
	open(CONCAT, ">$sequence_file") or die $!;
	foreach my $line (@concatenated_batch_file) {
		print CONCAT $line;
	}
	close(CONCAT) or die $!;
	foreach my $file (@batched_files) {
		unlink($file) or die "|".$file."|";
	}
	use_pre_downloaded:
	##############################################################################


	##############################################################################
	# Open genbank file into needed locations.
	skip_genbank_download:
	open(GENBANK, "<$sequence_file") or die $!;
	my @genbank 						= <GENBANK>;
	close GENBANK;
	my %binomial_name_hash				= ();
	my %exemplar_hash					= ();
	my %accession_hash					= ();
	my %taxonomy_hierarchy_hash 		= ();
	my $current_accession 				= 'NA';
	my @current_hierarchical_taxonomy 	= ();
	my $current_binomial_name 			= 'NA';
	my $found_organism_line 			= 0;
	my $found_accession 				= 0;
	my $organism_line_i 				= 0;
	foreach my $genbank_line (@genbank) {
		if ($genbank_line =~ m/REFERENCE/ && $found_organism_line == 1) {
			$binomial_name_hash{$current_accession} 				= $current_binomial_name;
			push(@{$taxonomy_hierarchy_hash{$current_accession}}, @current_hierarchical_taxonomy);
			# Reset all variables.
			$current_accession = 'NA';
			@current_hierarchical_taxonomy = ();
			$current_binomial_name = 'NA';
			$found_organism_line = 0;
			$found_accession = 0;
			$organism_line_i = 0;
		}
		if($genbank_line =~ m/ACCESSION/) {
			$found_accession = 1;
			$current_accession = $genbank_line;
			$current_accession =~ s/ACCESSION//g; # Parse off the ORGANISM tag
			$current_accession =~ s/^\s+//; # remove leading spaces
			my @accession_split = split(/ +/,$current_accession);
			$current_accession = $accession_split[0];
			$current_accession =~ s/\s+$//; # remove trailing spaces
			$current_accession =~ s/ /_/g; 	# replace inner whitespace
			if(exists($accession_hash{$current_accession})) {
				print $current_accession;
				$accession_hash{$current_accession}++;
			} else {
				$accession_hash{$current_accession} = 1;
			}
		}
		if($genbank_line =~ m/ORGANISM/ && $found_accession == 1) { $found_organism_line = 1 };
		if($found_organism_line == 1) {
			# Parse the first line, which is the binomial name
			if($organism_line_i == 0) {
				$current_binomial_name = $genbank_line;
				$current_binomial_name =~ s/ORGANISM//g; # Parse off the ORGANISM tag
				$current_binomial_name =~ s/^\s+//; # remove leading spaces
				$current_binomial_name =~ s/\s+$//; # remove trailing spaces
				$current_binomial_name =~ s/ /_/g; 	# replace inner whitespace
			} else {
			# Parse the remaining lines which contain the hierarchical taxonomy
				my $hierarchical_taxa_string = $genbank_line;
				$hierarchical_taxa_string =~ s/^\s+//;	# remove leading spaces
				$hierarchical_taxa_string =~ s/\s+$//;	# remove trailing spaces
				$hierarchical_taxa_string =~ s/ //g; 	# replace inner whitespace
				$hierarchical_taxa_string =~ s/\.//g; 	# delete periods
				my @hierarchical_taxa_list = split(/;/,$hierarchical_taxa_string);
				# print $hierarchical_taxa_string."\n";
				push(@current_hierarchical_taxonomy,@hierarchical_taxa_list);
			}
			$organism_line_i++;
		}
	}
	@genbank = (); # FLUSH.

	# close(GENBANK);
	##############################################################################
	my $seqin = Bio::SeqIO->new(-file   => $sequence_file,
								-format => 'genbank');
	##############################################################################

	##############################################################################
	# Loop through the downloaded genbank files and parse all the data
	my $seq_counter = 0;
	while (my $seq = $seqin->next_seq) {
		# Pull these values as you go along and parse the genbank file.
		# NCBI/Taxonomy Variables
		my $taxon_id 				= 'NA'; # Taxon ID for NCBI
		my $accession_number		= 'NA'; # NCBI accession number
		my $long_name 				= 'NA'; # Long description name for the sequence
		my $gene_name				= 'NA'; # Genbank gene name (e.g. COI)
		my $product_name			= 'NA'; # Genbank product name (e.g. cytochrome oxidase I)
		my $binomial_name			= 'NA';
		my @hierarchical_taxonomy	= 'NA';
		my $taxonomy_print_string 	= 'NA'; # Concatenated taxonomy string
		# Reference Variables
		my $publication_title 		= 'NA'; # Most recent publication title
		my $publication_authors		= 'NA'; # Pub authors
		my $abstract_text 			= 'NA'; # Abstract from pubmed article
		my $journal_name			= 'NA'; # Reference journal name
		my $journal_DOI				= 'NA'; # Online DOI access
		my $journal_SO				= 'NA'; # Journal SO (pub info)
		my $journal_volume			= 'NA'; # Journal volume #
		my $journal_issue			= 'NA'; # Journal issue #
		my $journal_pages			= 'NA'; # Journal page #'s
		my $journal_pubdate			= 'NA'; # Journal pubdate
		# Sequence Variables
		my $nucleotide_seq			= 'NA'; # Nucleotide sequence
		my $fasta_nucleotide		= 'NA'; # FASTA output
		my $fasta_to_print			= 'NA'; # FASTA file
		my $amino_acid_seq			= 'NA'; # Amino acid/protein sequence
		my $codon_start				= 'NA'; # Reading frame position
		my @pcr_primers				= 'NA'; # PCR primer list
		my $pcr_primer_print_string = 'NA'; # Concatenated PCR primer string
		my $protein_id				= 'NA'; # Protein NCBI ID
		# Sample Information
		my $collection_date			= 'NA'; #
		my $voucher_id				= 'NA'; #
		my $collected_by			= 'NA'; #
		my $identified_by			= 'NA'; #
		my $organelle				= 'NA'; #
		my $country					= 'NA'; #
		my $lat_lon					= 'NA'; #
		my $isolate_id				= 'NA'; #
		my $bio_material			= 'NA'; #
		##############################################################################
		# Grab some easy variables.
		$long_name 			= $seq->description if defined $seq->description;
		$accession_number 	= $seq->accession_number if defined $seq->accession_number;
		$nucleotide_seq		= $seq->seq if defined $seq->seq;
		if(length($nucleotide_seq) > $seq_length_maximum) {
			$nucleotide_seq = '---';
		}
		$binomial_name		= $binomial_name_hash{$accession_number};

		if(exists($binomial_name_hash{$accession_number})) {
			# Strip punctuation from binomial name
			$binomial_name =~ s/[^\w\s]//g;
			$binomial_name =~ s/\s/_/g;
			$fasta_nucleotide 	= '>'.$binomial_name.'_'.$accession_number.'$'.$nucleotide_seq;
		} else {
			$fasta_nucleotide 	= '>NA_'.$accession_number.'$'.$nucleotide_seq;
		}
		#$fasta_to_print		= '>'.$binomial_name.'_'.$accession_number."\n".$nucleotide_seq;
		##############################################################################

		##############################################################################
		# Obtain sequence feature information.
		# print "\t\t[".$seq_counter."] ".$long_name."\n";
		for my $feat_object ($seq->get_SeqFeatures) {
			for my $tag ($feat_object->get_all_tags) {
				for my $value ($feat_object->get_tag_values($tag)) {
					next if !defined($value); # Make sure value is defined.
					# print "    value: ", $value, "\n";
					if($value =~ m/taxon:/) { my @taxon_split = split(/:/,$value); $taxon_id = $taxon_split[1] };
					push(@pcr_primers, $value) 		 if($tag =~ m/PCR_primers/);
					$collected_by 			= $value if($tag =~ m/collected_by/);
					$collection_date 		= $value if($tag =~ m/collection_date/);
					$country 				= $value if($tag =~ m/country/);
					$identified_by			= $value if($tag =~ m/identified_by/);
					$lat_lon				= $value if($tag =~ m/lat_lon/);
					$organelle				= $value if($tag =~ m/organelle/);
					$codon_start			= $value if($tag =~ m/codon_start/);
					$gene_name				= $value if($tag =~ m/gene/);
					$product_name			= $value if($tag =~ m/product/);
					$protein_id				= $value if($tag =~ m/protein_id/);
					$amino_acid_seq			= $value if($tag =~ m/translation/);
					$voucher_id				= $value if($tag =~ m/specimen_voucher/);
					$isolate_id				= $value if($tag =~ m/isolate/);
					$bio_material			= $value if($tag =~ m/bio_material/);
				}
			}
		}
		# If targetting specimens that only have vouchers.
		if(($voucher_only == 1) && ($voucher_id eq 'NA')) {
			next;
		}
		# exemplar check
		if(!defined($binomial_name)) {
			next;
		}
		if(($exemplar_only == 1) && exists($exemplar_hash{$binomial_name})) {
			next;
		} else {
			$exemplar_hash{$binomial_name} = 1; # add binomial name to exemplar hash
		}
		if($accession_hash{$accession_number} > 1) {
			next;
		}
		##############################################################################
		# Obtain reference annotations.
		# print "\t\t\tRetrieving annotations...\n";
		my $anno_collection = $seq->annotation;
		for my $key ( $anno_collection->get_all_annotation_keys ) {
			my @annotations = $anno_collection->get_Annotations($key);
			for my $value ( @annotations ) {
				if ($value->tagname eq "reference") {
					# Store only the FIRST MOST RECENT reference.
					if(defined($value->authors) && ($publication_authors eq 'NA')) {
						$publication_authors = $value->authors;
					}
					if(defined($value->title) && ($publication_title eq 'NA')) {
						$publication_title = $value->title;
					}
				}
		   }
		}
		##############################################################################

		##############################################################################
		# Search for the pubmed article, get its ID, download the pubmed file.
		# Use the pubmed summary to get
		# Extract the abstract text from the pubmed file.
		my $cleaned_title = $publication_title;
		$cleaned_title =~ s/://g; # Remove colons
		# Check if the pubmed article has already been downloaded.
		# If it has, skip ahead.
		my $pubmed_was_stored = 0;
		if(exists $pubmed_hash_ref->{$cleaned_title}) {
			$pubmed_was_stored = 1;
			# print "\tUsing cached results for pubmed article.\n";
			goto stored_pubmed;
		}

		# Begin pubmed searching
		my $pubmed_esearch_tries = 0;
		my @pubmed_ids = ();
		goto skip_pubmed if $skip_pubmed_search == 0;
		# print "\tSearching pubmed...\n";
		pubmed_search:
		my $seq_pubmed = Bio::DB::EUtilities->new(   	-eutil    	=> 'esearch',
													   -db      	=> 'pubmed',
													   -retmax  	=> 1,
													   -email   	=> $user_email,
													   -term    	=> $cleaned_title,
													   -verbose		=> -1);

		eval {@pubmed_ids = $seq_pubmed->get_ids; };
		$seq_pubmed = undef; # Flush
		if($@ || scalar @pubmed_ids == 0) {
			$pubmed_esearch_tries++;
			print "\tProblem in pubmed search. Retrying...\n";
			if ($pubmed_esearch_tries == $max_pubmed_tries) {
				$pubmed_hash_ref->{$cleaned_title}->{'abstract'}= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'name'} 	= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'DOI'}		= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'SO'}		= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'volume'}	= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'issue'}	= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'pages'}	= 'NA';
				$pubmed_hash_ref->{$cleaned_title}->{'pubdate'} = 'NA';
				goto skip_pubmed;
			}
			goto pubmed_search;
		}

		my $pubmed_summary_tries = 0;
		pubmed_summary:
		my $ds_pubmed = '';
		my $seq_pubmed_summary = Bio::DB::EUtilities->new(   -eutil    	=> 'esummary',
															-db      	=> 'pubmed',
															-retmax  	=> 1,
															-email   	=> $user_email,
															-id    		=> $pubmed_ids[0]);
		# print "\tCould not find article for: ".$cleaned_title."\n" 	if !defined($pubmed_ids[0]);
		# Did find pubmed articles, get the summary
		eval { $ds_pubmed = $seq_pubmed_summary->next_DocSum; };
		if($@) {
			print "\tProblem in pubmed summary. Retrying...\n";

			goto skip_pubmed if $pubmed_summary_tries == $max_pubmed_tries;
			$pubmed_summary_tries++;
			goto pubmed_summary;
		}
		while (my $item = $ds_pubmed->next_Item('flattened'))  {
			next if !defined($item->get_content);
			$journal_name									= $item->get_content if $item->get_name =~ m/FullJournalName/;
			$pubmed_hash_ref->{$cleaned_title}->{'name'} 	= $item->get_content if $item->get_name =~ m/FullJournalName/;
			$journal_DOI 									= $item->get_content if $item->get_name =~ m/DOI/;
			$pubmed_hash_ref->{$cleaned_title}->{'DOI'}		= $item->get_content if $item->get_name =~ m/DOI/;
			$journal_SO 									= $item->get_content if $item->get_name =~ m/SO/;
			$pubmed_hash_ref->{$cleaned_title}->{'SO'}		= $item->get_content if $item->get_name =~ m/SO/;
			$journal_volume 								= $item->get_content if $item->get_name =~ m/Volume/;
			$pubmed_hash_ref->{$cleaned_title}->{'volume'}	= $item->get_content if $item->get_name =~ m/Volume/;
			$journal_issue 									= $item->get_content if $item->get_name =~ m/Issue/;
			$pubmed_hash_ref->{$cleaned_title}->{'issue'}	= $item->get_content if $item->get_name =~ m/Issue/;
			$journal_pages 									= $item->get_content if $item->get_name =~ m/Pages/;
			$pubmed_hash_ref->{$cleaned_title}->{'pages'}	= $item->get_content if $item->get_name =~ m/Pages/;
			$journal_pubdate 								= $item->get_content if $item->get_name =~ m/PubDate/;
			$pubmed_hash_ref->{$cleaned_title}->{'pubdate'} = $item->get_content if $item->get_name =~ m/PubDate/;
		}

		my $pubmed_fetch_tries = 0;
		pubmed_fetch:
		my $pubmed_fetch = Bio::DB::EUtilities->new( -eutil   => 'efetch',
													   -db      => 'pubmed',
													   -email   => $user_email,
													   -id      => $pubmed_ids[0]);
		my $pubmed_file = 'pubmed_'.$taxon_id.'.txt';
		eval { $pubmed_fetch->get_Response(-file => $pubmed_file); };
		if($@) {
			print "\tProblem in pubmed fetch. Retrying...\n";
			$pubmed_fetch_tries++;
			goto skip_pubmed if $pubmed_fetch_tries == $max_pubmed_tries;
			goto pubmed_fetch;
		}
		open (PUBMED, '<'.$pubmed_file);
		my @pubmed_xml_file = <PUBMED>;
		foreach my $line (@pubmed_xml_file) {
			$abstract_text = ''; # If you got this far, clear the NA from abstract text.
			if($line =~ m/AbstractText/) {
				$abstract_text .= $line;
				# Clear XML tags.
				$abstract_text =~ s/\<AbstractText\>//g;
				$abstract_text =~ s/\<\/AbstractText\>//g;
				$abstract_text =~ s/\<.*\>//g; # Remove all XML tags for sure.
				$abstract_text =~ s/^\s+//; # Remove leading space.
				$pubmed_hash_ref->{$cleaned_title}->{'abstract'} = $abstract_text;
				last;
			}
		}
		close(PUBMED);
		unlink($pubmed_file);
		skip_pubmed: # Skipped here because could not find any pubmed articles
		stored_pubmed: # Already had a stored pubmed article in the pubmed hash.
		if($pubmed_was_stored == 1) {
			$abstract_text		= $pubmed_hash_ref->{$cleaned_title}->{'abstract'};
			$journal_name 		= $pubmed_hash_ref->{$cleaned_title}->{'name'};
			$journal_DOI 		= $pubmed_hash_ref->{$cleaned_title}->{'DOI'};
			$journal_SO 		= $pubmed_hash_ref->{$cleaned_title}->{'SO'};
			$journal_volume 	= $pubmed_hash_ref->{$cleaned_title}->{'volume'};
			$journal_issue 		= $pubmed_hash_ref->{$cleaned_title}->{'issue'};
			$journal_pages 		= $pubmed_hash_ref->{$cleaned_title}->{'pages'};
			$journal_pubdate 	= $pubmed_hash_ref->{$cleaned_title}->{'pubdate'};
		}
		##############################################################################

		#############################################################################
		# Concatenate list data like taxonomy hierarchy and primer list.
		if (exists $taxonomy_hierarchy_hash{$accession_number}) {
			my @current_tax_array = @{$taxonomy_hierarchy_hash{$accession_number}};
			foreach my $taxa_i (0 .. $#current_tax_array) { # Print hierarchical taxonomy list
				last if $current_tax_array[$taxa_i] =~ m/COMMENT/;
				$taxonomy_print_string = '' if ($taxa_i == 0);
				$taxonomy_print_string .= $current_tax_array[$taxa_i].';' if ($taxa_i > 0);
			}
		}
		foreach my $primer_i (0 .. $#pcr_primers) { 		# PCR primer list
			$pcr_primer_print_string = '' if ($primer_i == 0);
			$pcr_primer_print_string .= $pcr_primers[$primer_i].';' if ($primer_i > 0);
		}
		#############################################################################

		#############################################################################
		# Prepare output strings

		$target_taxon =~ s/\"|^\s+|\s+$//g;
		$target_taxon =~ s/ |,/_/g;

		$taxon_id =~ s/\"|^\s+|\s+$//g;
		$taxon_id =~ s/ |,/_/g;

		$accession_number =~ s/\"|^\s+|\s+$//g;
		$accession_number =~ s/ |,/_/g;

		$search_options =~ s/\"|^\s+|\s+$//g;
		$search_options =~ s/ |,/_/g;

		$number_seqs_found =~ s/\"|^\s+|\s+$//g;
		$number_seqs_found =~ s/ |,/_/g;

		$long_name =~ s/\"|^\s+|\s+$//g;
		$long_name =~ s/ |,/_/g;

		$gene_name =~ s/\"|^\s+|\s+$//g;
		$gene_name =~ s/ |,/_/g;

		$product_name =~ s/\"|^\s+|\s+$//g;
		$product_name =~ s/ |,/_/g;

		if(defined($binomial_name)) {
			$binomial_name =~ s/\"|^\s+|\s+$//g;
			$binomial_name =~ s/ |,/_/g;
		}

		$taxonomy_print_string =~ s/\"|^\s+|\s+$//g;
		$taxonomy_print_string =~ s/ |,/_/g;

		$publication_title =~ s/\"|^\s+|\s+$//g;
		$publication_title =~ s/ |,/_/g;

		$publication_authors =~ s/\"|^\s+|\s+$//g;
		# $publication_authors =~ s/ |,/_/g;

		$abstract_text =~ s/\"|^\s+|\s+$//g;
		$abstract_text =~ s/ |,/_/g;

		$journal_name =~ s/\"|^\s+|\s+$//g;
		$journal_name =~ s/ |,/_/g;

		$journal_DOI =~ s/\"|^\s+|\s+$//g;
		$journal_DOI =~ s/ |,/_/g;

		$journal_SO =~ s/\"|^\s+|\s+$//g;
		$journal_SO =~ s/ |,/_/g;

		$journal_volume =~ s/\"|^\s+|\s+$//g;
		$journal_volume =~ s/ |,/_/g;

		$journal_issue =~ s/\"|^\s+|\s+$//g;
		$journal_issue =~ s/ |,/_/g;

		$journal_pages =~ s/\"|^\s+|\s+$//g;
		$journal_pages =~ s/ |,/_/g;

		$journal_pubdate =~ s/\"|^\s+|\s+$//g;
		$journal_pubdate =~ s/ |,/_/g;

		$nucleotide_seq =~ s/\"|^\s+|\s+$//g;
		$nucleotide_seq =~ s/ |,/_/g;

		$fasta_nucleotide =~ s/\"|^\s+|\s+$//g;
		$fasta_nucleotide =~ s/ |,/_/g;

		$amino_acid_seq =~ s/\"|^\s+|\s+$//g;
		$amino_acid_seq =~ s/ |,/_/g;

		$pcr_primer_print_string =~ s/\"|^\s+|\s+$//g;
		$pcr_primer_print_string =~ s/ |,/_/g;

		$codon_start =~ s/\"|^\s+|\s+$//g;
		$codon_start =~ s/ |,/_/g;

		$collection_date =~ s/\"|^\s+|\s+$//g;
		$collection_date =~ s/ |,/_/g;

		#$voucher_id =~ s/\"|^\s+|\s+$//g;
		#$voucher_id =~ s/ |,/_/g;

		$collected_by =~ s/\"|^\s+|\s+$//g;
		$collected_by =~ s/ |,/_/g;

		$identified_by =~ s/\"|^\s+|\s+$//g;
		$identified_by =~ s/ |,/_/g;

		$organelle =~ s/\"|^\s+|\s+$//g;
		$organelle =~ s/ |,/_/g;

		$country =~ s/\"|^\s+|\s+$//g;
		$country =~ s/ |,/_/g;

		$lat_lon =~ s/\"|^\s+|\s+$//g;
		$lat_lon =~ s/ |,/_/g;

		## Check nucleotide sequence against query
		my $distance = 0;
		my ($nw, $align) = '';
		my $aligned_seq = "NA";
		my $match_score = "NA";
		my ($k2p_distance, $transitions,$transversions,$bases_compared) = 0;
		my $aln_seq = '';
		if($match_aln_file ne '') {
			($k2p_distance, $transitions,$transversions,$bases_compared, $aln_seq) = aln2seq($nucleotide_seq, $match_seq);
		}

		my @current_output = (	$target_taxon,$dlm,
								$taxon_id,$dlm,
								$accession_number,$dlm,
								$search_options,$dlm,
								$number_seqs_found,$dlm,
								$long_name,$dlm,
								$gene_name,$dlm,
								$product_name,$dlm,
								$binomial_name,$dlm,
								$taxonomy_print_string,$dlm,
								$publication_title,$dlm,
								$publication_authors,$dlm,
								$abstract_text,$dlm,
								$journal_name,$dlm,
								$journal_DOI,$dlm,
								$journal_SO,$dlm,
								$journal_volume,$dlm,
								$journal_issue,$dlm,
								$journal_pages,$dlm,
								$journal_pubdate,$dlm,
								$nucleotide_seq,$dlm,
								$fasta_nucleotide,$dlm,
								$aln_seq,$dlm,
								$k2p_distance,$dlm,
								$amino_acid_seq,$dlm,
								$distance,$dlm,
								$pcr_primer_print_string,$dlm,
								$codon_start,$dlm,
								$collection_date,$dlm,
								$voucher_id,$dlm,
								$collected_by,$dlm,
								$identified_by,$dlm,
								$organelle,$dlm,
								$country,$dlm,
								$lat_lon,$dlm,
								$isolate_id,$dlm,
								$bio_material,$dlm);
		my @cleaned_output = ();
		foreach my $current_output (@current_output) {
			# Catch missing values. Probably should output an error here.
			$current_output = 'NA' if !defined $current_output;
			$current_output =~ s/\n//g; # Clear newlines.
			push(@cleaned_output, $current_output);
		}
		# This will be for the individual output file.
		if($suppress_output ne 'yes') {
			push(@output_lines,@cleaned_output,$endl);
		}
		# This will be returned at the end of the subroutine.
		push(@return_output_lines,@cleaned_output,$endl);
		# print "cleaned_output ".total_size(\@cleaned_output)."\n";
		# print "output_lines ".total_size(\@output_lines)."\n";
		# print "return_output_lines ".total_size(\@return_output_lines)."\n";
		##############################################################################

		#$fasta_to_print
		$seq_counter++;

		###
		###
	}
	$seqin = '';
	# Print to file if output is not suppressed.
	if($suppress_output ne 'yes') {

		unlink $output_file;
		# open (OUTPUT, '>>'.$output_file);
		open(OUTPUT, '>>'.$output_file) or die "Couldn't open: $!";

		foreach my $output_line (@output_lines) {
			print OUTPUT "\"".	$output_line if $output_line ne $endl;
			print OUTPUT 		$output_line if $output_line eq $endl;
		}
		close(OUTPUT);

	}
	#############################################################################

	taxa_failed:
	sequence_failed:
	just_count_seqs:
	if (exists $failed_search_hash_ref->{$target_taxon}) {
		my $failed_seq_taxa_id = $failed_search_hash_ref->{$target_taxon}->{'taxa_id'};
		push(@return_output_lines, 	$target_taxon, $dlm,
									$failed_seq_taxa_id, $dlm,
									'NA',$dlm,
									$search_options,$dlm,
									$number_seqs_found,$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									'NA',$dlm,
									$endl);
	}
	$taxa_counter++;
	# if(open(GENBANK)) {
		# print "File still open.";
		# die;
	# }
	# unlink $sequence_file or warn "Could not unlink $sequence_file\n";
	unlink($sequence_file) if defined $sequence_file;

	print "Duplicate accession Id's eliminated:\n";
	for my $unique_id (sort keys %accession_hash) {
		if($accession_hash{$unique_id} > 1) {
			print $unique_id." => ".$accession_hash{$unique_id}."\n";
		}
	}
	return \@return_output_lines;
}

my $t1 = Benchmark->new;
my $time_diff = timediff($t1, $t0);
print "\n";
print timestr($time_diff)."\n";

# Available database:
# :pubmed, protein, nucleotide, nuccore, nucgss, nucest, structure, genome,
  		        # :biosystems, books, cancerchromosomes, cdd, gap,
		        # :domains, gene, genomeprj, gensat, geo, gds,
		        # :homologene, journals, mesh, ncbisearch, nlmcatalog,
		        # :omia, omim, pepdome, pmc, popset, probe,
		        # :proteinclusters, pcassay, pccompound, pcsubstance,
		        # :seqannot, snp, sra, taxonomy, toolkit, toolkitall,
		        # :unigene, unists

sub search_strings {
	my $short_name = shift;
	my $query_file = shift;

	my %search_string_hash = ();
	$search_string_hash{'COI'} 	= "AND (COI[All Fields] OR \"cytochrome oxidase I\"[All Fields] OR \"cytochrome oxidase subunit I\"[All Fields] OR COX1[All Fields] OR \"COXI\"[All Fields]) NOT (\"complete genome\"[title] OR \"complete DNA\"[title])";
	$search_string_hash{'16S'} 	= "AND (16S[All Fields] OR \"16S ribosomal RNA\"[All Fields] OR \"16S rRNA\"[All Fields]) NOT (\"complete genome\"[title] OR \"complete DNA\"[title])";
	$search_string_hash{'18S'} 	= "AND (18S[All Fields] OR \"18S ribosomal RNA\"[All Fields] OR \"18S rRNA\"[All Fields]) NOT (\"complete genome\"[title] OR \"complete DNA\"[title])";
	$search_string_hash{'28S'} 	= "AND (28S[All Fields] OR \"28S ribosomal RNA\"[All Fields] OR \"28S rRNA\"[All Fields]) NOT (\"complete genome\"[title] OR \"complete DNA\"[title])";
	$search_string_hash{'CYTB'} = "AND (CYTB[All Fields] OR \"cytochrome b\"[All Fields] OR \"cyt b\"[All Fields]) NOT (\"complete genome\"[title] OR \"complete DNA\"[title])";
	$search_string_hash{'ND2'} 	= "AND (ND2[All Fields] OR \"NADH dehydrogenase subunit 2\"[All Fields] OR \"NADH2\"[All Fields]) NOT (\"complete genome\"[title] OR \"complete DNA\"[title])";

	if($query_file ne '') {
		open(QUERY, '<'.$query_file) or die "Couldn't open: $!";
		my @query_lines = <QUERY>;
		foreach my $query_line (@query_lines) {
			my @split_line = split(m/\=/, $query_line);
			print $split_line[0]." => ".$split_line[1]."\n";
			$search_string_hash{$split_line[0]} = $split_line[1];
		}
	}

	if(exists $search_string_hash{$short_name}) {
		print "Executing search string: \n";
		print $search_string_hash{$short_name}."\n";
		return $search_string_hash{$short_name};
	} else {
		return $short_name;
	}
}

sub current_memory {
	use Data::Dumper;
	use PadWalker qw(peek_my);
	use Devel::Size;

	my $hash_ref = peek_my(0);

	# print Dumper($hash_ref);
	print "Total mem: ".total_size($hash_ref)."\n";
}
