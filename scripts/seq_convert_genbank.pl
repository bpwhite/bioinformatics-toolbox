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

my $target_taxon = 'Tricorythodes';

##############################################################################
# Retrieve taxon ID
print "Searching for $target_taxon\n";
my $taxonomy_eutil = Bio::DB::EUtilities->new(-eutil    => 'esearch',
											   -db      => 'taxonomy',
											   -retmax  => 10,
											   -rettype => 'gb',
											   -email   => 'mymail@foo.bar',
											   -term    => $target_taxon);

my @taxon_ids = $taxonomy_eutil->get_ids();
my $taxon_id = $taxon_ids[0];
print "Found taxon ID: $taxon_id\n";
##############################################################################

##############################################################################
# Retrieve nucleotide sequences
print "Retrieving sequences...\n";
my $sequence_term = 'txid'.$taxon_id.'[Organism:exp]';
my $sequence_search = Bio::DB::EUtilities->new(-eutil    => 'esearch',
											   -db      => 'nucleotide',
											   -retmax  => 1,
											   -rettype => 'gb',
											   -email   => 'mymail@foo.bar',
											   -term    => $sequence_term);
my @sequence_ids = $sequence_search->get_ids();
foreach my $seq_id (@sequence_ids) {
	print $seq_id."\n";
}
my $sequence_fetch = Bio::DB::EUtilities->new( -eutil   => 'efetch',
											   -db      => 'nucleotide',
											   -rettype => 'gb',
											   -email   => 'mymail@foo.bar',
											   -id      => \@sequence_ids);
my $sequence_file = 'seqs_'.$taxon_id.'.gb';
$sequence_fetch->get_Response(-file => $sequence_file);
##############################################################################

##############################################################################
my $seqin = Bio::SeqIO->new(-file   => $sequence_file,
                            -format => 'genbank');

my $dumper_file = "dumper.txt";

# unlink $dumper_file;
# open (DUMP, '>>'.$dumper_file);

##############################################################################
# Loop through the downloaded genbank files and parse all the data
while (my $seq = $seqin->next_seq) {
	# Pull these values
	my $taxon_id 			= ''; # Taxon ID for NCBI
	my $publication_title 	= ''; # Most recent publication title
	# my $publication_authors	= ''; # 
	# my $long_name 			= '';
	# my $publication_title 	= '';
	# my $publication_authors	= '';
	# my $taxon_id 			= '';
	# my $publication_title 	= '';
	# my $publication_authors	= '';
	
	print $seq->display_name."\n";
	print $seq->description."\n";
	print $seq->accession_number."\n";
	
	##############################################################################
	print "\nSequence Features\n\n";
	for my $feat_object ($seq->get_SeqFeatures) {          
	print "primary tag: ", $feat_object->primary_tag, "\n";          
	for my $tag ($feat_object->get_all_tags) {             
		  print "  tag: ", $tag, "\n";
			
		  for my $value ($feat_object->get_tag_values($tag)) {                
			 print "    value: ", $value, "\n";
			if($value =~ m/taxon:/) { 
				my @taxon_split = split(/:/,$value);
				$taxon_id = $taxon_split[1];
			}
		  }          
	   }       
	}
	##############################################################################
	print "\nAnnotations\n\n";
	my $anno_collection = $seq->annotation;
	for my $key ( $anno_collection->get_all_annotation_keys ) {
		my @annotations = $anno_collection->get_Annotations($key);
		for my $value ( @annotations ) {
			print "tagname : ", $value->tagname, "\n";
			  # $value is an Bio::Annotation, and also has an "as_text" method
			print "  annotation value: ", $value->display_text, "\n";
			if ($value->tagname eq "reference") {
				if(defined($value->authors)) { print "author: ",$value->authors(), "\n"};
				if(defined($value->title)) { print "title: ",$value->title(), "\n"};
				if(defined($value->medline)) { print "medline: ",$value->medline(), "\n"};
				if(defined($value->editors)) { print "editors: ",$value->editors(), "\n"};
				if(defined($value->database)) { print "database: ",$value->database(), "\n"};
				if(defined($value->pubmed)) { print "pubmed: ",$value->pubmed(), "\n"};
				if(defined($value->location)) {	print "location: ",$value->location(), "\n"};
				if(defined($value->doi)) {	print "doi: ",$value->doi(), "\n"};
				if(defined($value->gb_reference)) {	print "gb_reference: ",$value->location(), "\n"};
			}
				# my $hash_ref = $value->hash_tree;
				# for my $key (keys %{$hash_ref}) {
					# if(defined($hash_ref->{$key})) {
						# print "\t".$key.": ".$hash_ref->{$key}."\n";
					# }
				# }
			# }
	   }
	}

	##############################################################################
	my $seq_taxonomy = Bio::DB::EUtilities->new(   -eutil    => 'esummary',
												   -db      => 'taxonomy',
												   -retmax  => 1,
												   -rettype => 'gb',
												   -email   => 'mymail@foo.bar',
												   -id    => $taxon_id);
	my ($name)  = $seq_taxonomy->next_DocSum->get_contents_by_name('ScientificName');
	##############################################################################
	
	##############################################################################
	print $name."\n";
	print "\nSequence\n\n";
	print ">$name\n";
	print $seq->seq."\n";
	##############################################################################

}
##############################################################################


