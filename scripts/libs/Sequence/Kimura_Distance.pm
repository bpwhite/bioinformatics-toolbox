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
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw(k2p_bootstrap k2p_no_bs);

sub k2p_bootstrap {
	my $seq1 		= shift;
	my $seq2 		= shift;
	my $length 		= shift;
	my $weights		= shift;
	
	$$seq1 =~ s/-/*/g;
	# XOR values for each transition/transversion
	# Transitions
	# A->G = 6
	# C->T = 23
	# Transversions
	# A->C = 2
	# A->T = 21
	# C->G = 4
	# G->T = 19
	my @unpacked_seq1 = unpack("C*", $$seq1);
	my @unpacked_seq2 = unpack("C*", $$seq2);
	my $transitions 		= 0;
	my $transversions 		= 0;
	my $num_comparisons 	= 0;
	my $pair = '';
	for (my $i = 0; $i < $length; $i++) {
		if ($unpacked_seq1[$i] == $unpacked_seq2[$i]) {
			$num_comparisons += @{$weights}[$i];
		} else {
			$pair = $unpacked_seq1[$i]^$unpacked_seq2[$i];
			# Transitions
			if($pair 		== 6) { # A->G
				$transitions += @{$weights}[$i];
				$num_comparisons += @{$weights}[$i];
			}
			elsif($pair 	== 23){ # C->T	
				$transitions += @{$weights}[$i];
				$num_comparisons += @{$weights}[$i];
			} 
			# Transversions
			elsif($pair 	== 21){ # A->T
				$transversions += @{$weights}[$i];
				$num_comparisons += @{$weights}[$i];
			}
			elsif($pair 	== 2){ 	# A->C
				$transversions += @{$weights}[$i];
				$num_comparisons += @{$weights}[$i];
			}
			elsif($pair 	== 4){ 	# C->G
				$transversions += @{$weights}[$i];
				$num_comparisons += @{$weights}[$i];
			}
			elsif($pair 	== 19){ # C->T
				$transversions += @{$weights}[$i];
				$num_comparisons += @{$weights}[$i];
			}
		}
	}
	
	my $P = $transitions/$num_comparisons;
	my $Q = $transversions/$num_comparisons;
	my $w1 = 1-2*$P-$Q;
	my $w2 = 1-2*$Q;
	if(($w1 <= 0) || ($w2 <= 0)) { return (2,$transitions,$transversions,$num_comparisons) };
	my $K2P = -log($w1)/2-log($w2)/4;

	return($K2P, $transitions, $transversions, $num_comparisons);
}


sub k2p_no_bs{
	my $seq1 		= shift;
	my $seq2 		= shift;
	my $length 		= shift;
	
	# XOR values for each transition/transversion
	# Transitions
	# A->G = 6
	# C->T = 23
	# Transversions
	# A->C = 2
	# A->T = 21
	# C->G = 4
	# G->T = 19
	my @unpacked_seq1 = unpack("C*", $$seq1);
	my @unpacked_seq2 = unpack("C*", $$seq2);
	my $transitions 		= 0;
	my $transversions 		= 0;
	my $num_comparisons 	= 0;
	my $pair = '';
	for (my $i = 0; $i < $length; $i++) {
		if ($unpacked_seq1[$i] == $unpacked_seq2[$i]) {
			$num_comparisons++;
		} else {
			$pair = $unpacked_seq1[$i]^$unpacked_seq2[$i];
			# Transitions
			if($pair 		== 6) { # A->G
				$transitions++;
				$num_comparisons++;
			}
			elsif($pair 	== 23){ # C->T	
				$transitions++;
				$num_comparisons++;
			} 
			# Transversions
			elsif($pair 	== 21){ # A->T
				$transversions++;
				$num_comparisons++;
			}
			elsif($pair 	== 2){ 	# A->C
				$transversions++;
				$num_comparisons++;
			}
			elsif($pair 	== 4){ 	# C->G
				$transversions++;
				$num_comparisons++;
			}
			elsif($pair 	== 19){ # C->T
				$transversions++;
				$num_comparisons++;
			}
		}
	}
	
	my $P = $transitions/$num_comparisons;
	my $Q = $transversions/$num_comparisons;
	my $w1 = 1-2*$P-$Q;
	my $w2 = 1-2*$Q;
	my $K2P = -log($w1)/2-log($w2)/4;

	return($K2P, $transitions, $transversions, $num_comparisons);
}

sub calc_k2p {
	my $transitions 	= shift;
	my $transversions 	= shift;
	my $num_comparisons = shift;
	
	my $P = $transitions/$num_comparisons;
	my $Q = $transversions/$num_comparisons;
	my $w1 = 1-2*$P-$Q;
	my $w2 = 1-2*$Q;
	my $K2P = -log($w1)/2-log($w2)/4;

	return($K2P, $transitions, $transversions, $num_comparisons);
}

1;