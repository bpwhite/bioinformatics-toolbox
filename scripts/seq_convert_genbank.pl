# Downloads and parses genbank files given an input taxon name.
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
use lib "$FindBin::Bin/libs"; 
use Time::HiRes qw( usleep );

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
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Class::Inspector;

use Data::Dumper;

my %arguments = ();
$arguments{'-list'} = ''; # List file name 
$arguments{'-slim'} = '1'; # Sequence limit
$arguments{'-tlim'} = '1'; # Taxa limit
$arguments{'-user_email'} = 'blah@blah.com';

# Loop through arguments
for (my $i = 0; $i <= scalar(@ARGV)/2; $i+=2) {
	my $argument = $ARGV[$i];
	my $parameter = $ARGV[$i+1];
	$arguments{$argument} = $parameter;
	print $argument." => ".$parameter."\n";
}

# Print parameters
print "Parameters:\n";
while ( my ($key, $value) = each(%arguments) ) {
        print "$key => $value\n";
}

exit;

open (TAXALIST, '<'.$taxa_file);
my @taxa_list = <TAXALIST>;

##############################################################################
my @overall_results = ();
my @overall_taxa_failures = ();
my @overall_seq_failures = ();
my $taxa_counter = 0;
my $overall_output_file = 'overall_output.csv';
my $overall_taxa_failures_file = 'overall_taxa_failures.csv';
my $overall_seq_failures_file = 'overall_seq_failures.csv';
my $endl = "\n";
my $suppress_output = 'yes';
foreach my $taxa (@taxa_list) {
	$taxa =~ s/\n//g; # replace newlines
	my ($results, $taxa_failures, $seq_failures) = download_target_taxa($taxa,$taxa_counter,$suppress_output);
	push(@overall_results, @$results);
	push(@overall_taxa_failures, $taxa_failures,$endl);
	push(@overall_seq_failures, $seq_failures,$endl);
	$taxa_counter++;
}
##############################################################################

##############################################################################
# Print successful output
unlink $overall_output_file;
open (OVEROUTPUT, '>>'.$overall_output_file);
foreach my $output_line (@overall_results) {
	print OVEROUTPUT "\"".	$output_line if $output_line ne $endl;
	print OVEROUTPUT 		$output_line if $output_line eq $endl;
}
close(OVEROUTPUT);
##############################################################################

##############################################################################
# Print things that failed taxonomic lookup
unlink $overall_taxa_failures_file;
open (OVERTAXAFAILS, '>>'.$overall_taxa_failures_file);
foreach my $output_line (@overall_taxa_failures) {
	print OVERTAXAFAILS "\"".	$output_line if $output_line ne $endl;
	print OVERTAXAFAILS 		$output_line if $output_line eq $endl;
}
close(OVERTAXAFAILS);
##############################################################################

##############################################################################
# Print things that lacked sequence data
unlink $overall_seq_failures_file;
open (OVERSEQFAILS, '>>'.$overall_seq_failures_file);
foreach my $output_line (@overall_seq_failures) {
	print OVERSEQFAILS "\"".	$output_line if $output_line ne $endl;
	print OVERSEQFAILS 		$output_line if $output_line eq $endl;
}
close(OVERSEQFAILS);
##############################################################################

sub download_target_taxa {
	##############################################################################
	# Subroutine parameters
	my $target_taxon = shift;
	my $taxa_counter = shift;
	my $suppress_output = shift;
	##############################################################################
	
	##############################################################################
	# my $skip_pubmed_search = 0;
	my $skip_pubmed_search = 1;
	my $number_seqs_found = 'NA';
	# my $target_taxon = 'Crassostrea angulata';
	my $search_options = 'AND COI[gene]';
	my $taxa_failed = 0; # Flag for if the taxa lookup fails
	my $failed_taxa = 'NA';
	my $sequence_failed = 0; # Flag for if the sequence lookup fails
	my $failed_sequence = 'NA';
	# my $search_options = '';
	my $taxon_limit = 1;
	my $sequence_limit = 10000;
	my $user_email = 'bpcwhite@gmail.com';
	my $dlm = ',';
	my $endl = "\n";
	my $output_file = $target_taxon."_output.csv";
	my $sleep_time = 1000; # microsends, pause between seconds
	my $max_num_tries = 15;
	my $max_pubmed_tries = 1;
	##############################################################################
	
	##############################################################################
	# Retrieve taxon ID
	my $taxonomy_eutil_tries = 0;
	taxonomy_eutil:
	print "[".$taxa_counter."] Searching for $target_taxon\n";
	my $taxonomy_eutil = Bio::DB::EUtilities->new(-eutil    => 'esearch',
												   -db      => 'taxonomy',
												   -retmax  => $taxon_limit,
												   -rettype => 'gb',
												   -email   => $user_email,
												   -term    => $target_taxon);
										   
	my @taxon_ids = ();
	eval { @taxon_ids = $taxonomy_eutil->get_ids(); };
	if ($@) {
		print "\tProblem in taxonomy_eutil. Retrying...\n";
		$taxonomy_eutil_tries++;
		$taxa_failed = 1 if $taxonomy_eutil_tries == $max_num_tries;
		goto taxa_failed if $taxonomy_eutil_tries == $max_num_tries;
		goto taxonomy_eutil;
	}
	$taxa_failed = 1 if scalar @taxon_ids == 0;
	goto taxa_failed if scalar @taxon_ids == 0;
	my $taxon_id = $taxon_ids[0];
	print "\tFound taxon ID: $taxon_id\n";
	##############################################################################

	##############################################################################
	# Retrieve nucleotide sequences
	my $sequence_search_tries = 0;
	sequence_search:
	print "\tRetrieving sequences...\n";
	my $sequence_term = 'txid'.$taxon_id.'[Organism:exp]';
	my $sequence_search = Bio::DB::EUtilities->new(-eutil    => 'esearch',
												   -db      => 'nucleotide',
												   -retmax  => $sequence_limit,
												   -rettype => 'gb',
												   -email   => $user_email,
												   -term    => $sequence_term.' '.$search_options);
	my @sequence_ids = ();
	eval { @sequence_ids = $sequence_search->get_ids(); };
	if ($@) {
		print "\tProblem in sequence_search. Retrying...\n";
		$sequence_search_tries++;
		$sequence_failed = 1 if $sequence_search_tries == $max_num_tries;
		goto sequence_failed if $sequence_search_tries == $max_num_tries;
		goto sequence_search;
	}
	print "\tFound ".scalar(@sequence_ids)." sequences to download.\n";
	goto sequence_failed if scalar(@sequence_ids) == 0;
	my $sequence_file = 'seqs_'.$taxon_id.'.gb';
	# if(-e $sequence_file) {
		# goto skip_genbank_download;
	# }
	
	my $sequence_download_tries = 0;
	sequence_download:
	print "\tDownloading sequences...\n";
	my $sequence_fetch = Bio::DB::EUtilities->new( -eutil   => 'efetch',
												   -db      => 'nucleotide',
												   -rettype => 'gb',
												   -email   => $user_email,
												   -id      => \@sequence_ids);

	eval { $sequence_fetch->get_Response(-file => $sequence_file); };
	if($@) {
		print "\tProblem in sequence download. Retrying...\n";
		$sequence_download_tries++;
		$sequence_failed = 1 if $sequence_download_tries == $max_num_tries;
		goto sequence_failed if $sequence_download_tries == $max_num_tries;
		goto sequence_download;
	}
	print "\tSequences downloaded to genbank format.\n";
	##############################################################################

	##############################################################################
	# Open genbank file into needed locations.
	skip_genbank_download:
	open (GENBANK, '<'.$sequence_file);
	my @genbank 						= <GENBANK>;
	my %binomial_name_hash				= ();
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
	##############################################################################
	my $seqin = Bio::SeqIO->new(-file   => $sequence_file,
								-format => 'genbank');
	##############################################################################
	# Output file headers.
	my @output_lines = ();
	my @return_output_lines = ();
	push(@output_lines,
		'taxon_query',$dlm,
		'taxon_id',$dlm,
		'accession',$dlm,
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
		$endl);
		push(@return_output_lines, @output_lines) if($taxa_counter == 0);
	##############################################################################
	# Loop through the downloaded genbank files and parse all the data
	my $seq_counter = 0;
	while (my $seq = $seqin->next_seq) {
		print "\tSeq #: ".$seq_counter."\n";
		usleep($sleep_time); 	# Sleep so you don't overload NCBI's servers.
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
		
		
		##############################################################################
		# Grab some easy variables.
		$long_name 			= $seq->description if defined $seq->description;
		$accession_number 	= $seq->accession_number if defined $seq->accession_number;
		$nucleotide_seq		= $seq->seq if defined $seq->seq;
		$binomial_name		= $binomial_name_hash{$accession_number};
		##############################################################################
		
		##############################################################################
		# Obtain sequence feature information.
		print "\tRetrieving sequence features...\n";
		print "\t".$long_name."\n";
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
				}
			}
		}
		##############################################################################
		# Obtain reference annotations.
		print "\tRetrieving annotations...\n";
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
		my $pubmed_esearch_tries = 0;
		my @pubmed_ids = ();
		goto skip_pubmed if $skip_pubmed_search == 1;
		pubmed_search:
		my $seq_pubmed = Bio::DB::EUtilities->new(   	-eutil    	=> 'esearch',
													   -db      	=> 'pubmed',
													   -retmax  	=> 1,
													   -email   	=> $user_email,
													   -term    	=> $cleaned_title,
													   -verbose		=> -1);

		eval {@pubmed_ids = $seq_pubmed->get_ids; };
		if($@) {
			print "\tProblem in pubmed search. Retrying...\n";
			$pubmed_esearch_tries++;
			goto skip_pubmed if $pubmed_esearch_tries == $max_pubmed_tries;
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
		print "\tCould not find article for: ".$cleaned_title."\n" 	if !defined($pubmed_ids[0]);
		# Did find pubmed articles, get the summary
		eval { $ds_pubmed = $seq_pubmed_summary->next_DocSum; };
		if($@) {
			print "\tProblem in pubmed summary. Retrying...\n";
			$pubmed_summary_tries++;
			goto skip_pubmed if $pubmed_summary_tries == $max_pubmed_tries;
			goto pubmed_summary;
		}
		
		while (my $item = $ds_pubmed->next_Item('flattened'))  {
			next if !defined($item->get_content);
			$journal_name 		= $item->get_content if $item->get_name =~ m/FullJournalName/;
			$journal_DOI 		= $item->get_content if $item->get_name =~ m/DOI/;
			$journal_SO 		= $item->get_content if $item->get_name =~ m/SO/;
			$journal_volume 	= $item->get_content if $item->get_name =~ m/Volume/;
			$journal_issue 		= $item->get_content if $item->get_name =~ m/Issue/;
			$journal_pages 		= $item->get_content if $item->get_name =~ m/Pages/;
			$journal_pubdate 	= $item->get_content if $item->get_name =~ m/PubDate/;
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
		$abstract_text = ''; # If you got this far, clear the NA from abstract text.
		open (PUBMED, '<'.$pubmed_file);
		my @pubmed_xml_file = <PUBMED>;
		foreach my $line (@pubmed_xml_file) {
			if($line =~ m/AbstractText/) {
				$abstract_text .= $line;
				# Clear XML tags.
				$abstract_text =~ s/\<AbstractText\>//g;
				$abstract_text =~ s/\<\/AbstractText\>//g;
				$abstract_text =~ s/\<.*\>//g; # Remove all XML tags for sure.
				$abstract_text =~ s/^\s+//; # Remove leading space.
			}
		}
		close(PUBMED);
		unlink($pubmed_file);
		skip_pubmed: # Skipped here because could not find any pubmed articles
		##############################################################################
		
		#############################################################################
		# Concatenate list data like taxonomy hierarchy and primer list.
		if (exists $taxonomy_hierarchy_hash{$accession_number}) {
			my @current_tax_array = @{$taxonomy_hierarchy_hash{$accession_number}};
			foreach my $taxa_i (0 .. $#current_tax_array) { # Print hierarchical taxonomy list
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
		my @current_output = (	$target_taxon,$dlm,
								$taxon_id,$dlm,
								$accession_number,$dlm,
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
								$amino_acid_seq,$dlm,
								$pcr_primer_print_string,$dlm,
								$codon_start,$dlm,
								$collection_date,$dlm,
								$voucher_id,$dlm,
								$collected_by,$dlm,
								$identified_by,$dlm,
								$organelle,$dlm,
								$country,$dlm,
								$lat_lon,$dlm);
		my @cleaned_output = ();
		foreach my $current_output (@current_output) {
			$current_output =~ s/\n//g; # Clear newlines.
			push(@cleaned_output, $current_output);
		}
		# This will be for the individual output file.
		push(@output_lines,@cleaned_output,$endl);
		# This will be returned at the end of the subroutine.
		push(@return_output_lines,@cleaned_output,$endl);
		##############################################################################
		$seq_counter++;
	}
	
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
	##############################################################################
	
	taxa_failed:
	$failed_taxa = $target_taxon if $taxa_failed == 1;
	sequence_failed:
	$failed_sequence = $target_taxon if $sequence_failed == 1;
	$taxa_counter++;
	
	unlink($sequence_file) if defined $sequence_file; # Delete the genbank file.
	return \@return_output_lines,$failed_taxa,$failed_sequence;
}

# Available database:
# :pubmed, protein, nucleotide, nuccore, nucgss, nucest, structure, genome,
  		        # :biosystems, books, cancerchromosomes, cdd, gap,
		        # :domains, gene, genomeprj, gensat, geo, gds,
		        # :homologene, journals, mesh, ncbisearch, nlmcatalog,
		        # :omia, omim, pepdome, pmc, popset, probe,
		        # :proteinclusters, pcassay, pccompound, pcsubstance,
		        # :seqannot, snp, sra, taxonomy, toolkit, toolkitall,
		        # :unigene, unists
