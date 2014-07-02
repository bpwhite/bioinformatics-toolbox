# BLAST Functions
#
# Copyright (c) 2013, 2014 Bryan White

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
use strict;
use warnings;

use Sequence::Fasta;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw( blast_output	
						);

sub blast_output {
	my %p = (
		seq_string 		=> '',
		seq_id			=> '',
		blast_db 		=> '',
		outfmt			=> '',
		max_target_seqs => '1',
		num_threads 	=> '',
		identity_level	=> '',
		aln_length_pcnt	=> '',
		log_file		=> '',
		@_
	);
	
	defined $p{seq_string} or die 'Required parameter "seq_string" not defined';
	
	my $seq_string 		= $p{'seq_string'};
	my $seq_id			= $p{'seq_id'};
	my $blastdb_name 	= $p{'blast_db'};
	my $identity_level 	= $p{'identity_level'};
	my $aln_length_pcnt = $p{'aln_length_pcnt'};
	my $output_log		= $p{'log_file'};
	my $max_target_seqs = $p{'max_target_seqs'};
	
	my $query_length = fast_seq_length($seq_string);
	
    my $temp_query .= 'temp_query_'.int(rand(99999)).'.fas';
	
	open (QUERY, '>'.$temp_query);
	print QUERY ">".$seq_id."\n";
	print QUERY $seq_string."\n";
	close (QUERY);
	# my $blast_output = '';
	my $blast_output = `blastn -query $temp_query -db $blastdb_name -outfmt 6 -max_target_seqs $max_target_seqs -num_threads 2`;
	# print "blastn -query $temp_query -db $blastdb_name -outfmt 6 -max_target_seqs 1 -num_threads 1";
	# print $blast_output."\n";
	unlink $temp_query;
	
	my $seq_fail_i = 0;
	my $seq_pass_i = 0;
	my @blast_lines = split(/\n/,$blast_output);
	my @split_blast = ();
	my $final_blast_output = '';
	foreach my $blast_line (@blast_lines) {
		@split_blast = split(/\t/,$blast_line);
		if ($blast_output eq '') {
			$seq_fail_i++;
			next;
		}
		$blast_output =~ s/\n//g;
		# BLAST % identity
		if($split_blast[2] < $identity_level) {
			$seq_fail_i++;
			open (OUTLOG, '>>', $output_log);
			print OUTLOG $blast_output.",Fail homology level\n";
			close (OUTLOG);
			# print $blast_output."\n";
			# print "Fail homology level\n";
			next;
		}
		# BLAST alignment length cutoff
		# Useful for excluding chimeric sequences
		if ($aln_length_pcnt > 0) {
			my $pcnt_query_match = $split_blast[3]/$query_length;
			if($pcnt_query_match < $aln_length_pcnt) {
				$seq_fail_i++;
				open (OUTLOG, '>>', $output_log);
				print OUTLOG $blast_output.",Fail alignment length\n";
				close (OUTLOG);
				# print $blast_output."\n";
				# print "Fail alignment length at $pcnt_query_match\n";
				next;
			}
		}
		$final_blast_output = $blast_line;
		$seq_pass_i++;
		last;
	}
	return 'NA' if $seq_pass_i == 0;
	return 'NA' if (!defined($split_blast[0]));
	
	$final_blast_output =~ s/\n//g;
	# print $final_blast_output."\n";
	return $final_blast_output;
}

1;