# Calculates the Kimura 2 parameter distance (Kimura 1980) between two sequences.
#
# Copyright (c) 2013, Bryan White

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
# use Inline::Files::Virtual;
use Inline::Files;
use Inline C;

use strict;
use warnings;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw(c_k2p_bootstrap);

sub k2p_bootstrap {
	my $seq1 		= shift;
	my $seq2 		= shift;
	my $length 		= shift;
	my $weights 	= shift;
	
	my %lookup_table = (
						# Transitions
						'AG'	=> 'S',
						'GA'	=> 'S',
						'CT'	=> 'S',
						'TC'	=> 'S',
						# Transversions
						'AC'	=> 'V',
						'CA'	=> 'V',
						'AT'	=> 'V',
						'TA'	=> 'V',
						'CG'	=> 'V',
						'GC'	=> 'V',
						'GT'	=> 'V',
						'TG'	=> 'V',
						# Identical
						'AA'	=> 'I',
						'GG'	=> 'I',
						'CC'	=> 'I',
						'TT'	=> 'I'
					);
	my %results = (
					'S' => 0,
					'V' => 0,
					'I' => 0
					);
	my $comparisons = '';
	my $i = 0;
	while ($i < $length) {
		# print "[".$i."]\n";
		# print $seq1_split[$i].$seq2_split[$i]."\n";
		# print $lookup_table{$seq1_split[$i].$seq2_split[$i]}."\n";
		# exit;
		# my $pair = $lookup_table{$seq1_split[$i].$seq2_split[$i]};
		my $pair = $lookup_table{substr($$seq1, $i,1).substr($$seq2, $i,1)};
		# my $weight = int @{$weights}[$i];
		# $results{$pair} += $weight if defined $pair;
		$results{$pair} += @{$weights}[$i] if defined $pair;
		$i++;
		# if (!exists($lookup_table{substr($$seq1, $i,1).substr($$seq2, $i,1)})) {
			# $i++;
			# next;
		# }
		# my $weight = int @{$weights}[$i];
		# print $weight."\n";
		# if ($weight == 0) {
			# $i++;
			# next;
		# }
		
		
		# print $weight." => ".$pair."\n" if defined $pair;
		# for (my $j = 0; $j < $weight; $j++) {
			# $comparisons .= $pair;
		# }
		# $i++;
		# print $pair."\n";
	}
	# print $comparisons."\n";
	# my $transitions 		= $comparisons =~ tr/S//;
	# my $transversions 		= $comparisons =~ tr/V//;
	# my $num_comparisons 	= length($comparisons);
	my $transitions 		= $results{'S'};
	my $transversions 		= $results{'V'};
	my $num_comparisons 	= $transitions + $transversions + $results{'I'};
	my $P = $transitions/$num_comparisons;
	my $Q = $transversions/$num_comparisons;
	# print $P." => ".$Q."\n";
	my $w1 = 1-2*$P-$Q;
	my $w2 = 1-2*$Q;
	# print $w1." => ".$w2."\n";
	my $K2P = -log($w1)/2-log($w2)/4;
	# print $transitions." => ".$transversions." => ".$num_comparisons."\n";
	# print $K2P."\n";
	return($K2P, $transitions, $transversions, $num_comparisons);
}

__END__
__C__

void c_k2p_bootstrap(char* str1, char* str2, float critical_value,
						float cutoff, float gene_length,
						SV* sv_transitions, SV* sv_transversions,
						SV* sv_comparisons, SV* sv_K2P, char *bs_weights) {
	/* Returns the following:
		transitions
		transversions
		comparisons
		K2P
		variance
		stderror
		min (confidence intervals)
		max (confidence intervals)
		
		XOR values for each transition/transversion
			Transitions
			A->G = 6
			C->T = 23
			Transversions
			A->C = 2
			A->T = 21
			C->G = 4
			G->T = 19
	*/

	int transitions = 0;
	int transversions = 0;
	int comparisons = 0;
	short int i, base_pair = 0;
	float 	P, Q, w1, w2, K2P = 0;
	for(i = 0; i < strlen(str1); i++) {
		int current_weight = bs_weights[i] - '0'; // Convert array of char to single int 
		if(str1[i] == str2[i]) {
			// Skip identicals.
			comparisons += current_weight;
		} else {
			base_pair = str1[i] ^ str2[i]; // XOR base pair together
			// Transitions
			if(base_pair 		== 6) { // A->G
				transitions += current_weight;
				comparisons += current_weight;
			}
			else if(base_pair 	== 23){ // C->T	
				transitions += current_weight;
				comparisons += current_weight;
			} // Transversions
			else if(base_pair 	== 21){ // A->T
					transversions += current_weight;
					comparisons += current_weight;
			}
			else if(base_pair 	== 2){ 	// A->C
					transversions += current_weight;
					comparisons += current_weight;
			}
			else if(base_pair 	== 4){ 	// C->G
					transversions += current_weight;
					comparisons += current_weight;
			}
			else if(base_pair 	== 19){ // C->T
					transversions += current_weight;
					comparisons += current_weight;
			}
		}
	}
	
	// Compute Kimura Distance
	P = (float)transitions/comparisons;
	Q = (float)transversions/comparisons;
	w1 = 1-2*P-Q;
	w2 = 1-2*Q;
	K2P = -log(w1)/2-log(w2)/4;

	// Modify variables.
	sv_setnv(sv_transitions,transitions);
	sv_setnv(sv_transversions,transversions);
	sv_setnv(sv_comparisons,comparisons);
	sv_setnv(sv_K2P,K2P);
}