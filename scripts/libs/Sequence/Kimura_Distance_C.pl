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
use Inline::Files;
use Inline C;

use strict;
use warnings;

use Exporter;

our @ISA= qw( Exporter );
our @EXPORT_OK = qw(c_kimura_distance c_k2p_bootstrap);

__C__

#include <stdio.h>
#include <math.h>
#include <string.h>

float bit_Test () {
	int a = 65;
	int g = 71;
	int t = 84;
	int c = 67;
	char hyphen = '-';
	char asterisk = '*';
	int d;
	d = a ^ a;
	printf("a ^ a: %i %d %f %o\n",d,d,d,d);
	d = a ^ g;
	d <<= d;
	printf("a ^ g: %i %d %f %o\n",d,d,d,d);
	d = c ^ t;
	d <<= d;
	printf("c ^ t: %i %d %f %o\n",d,d,d,d);
	d = a ^ c;
	d <<= d;
	printf("a ^ c: %i %d %f %o\n",d,d,d,d);
	d = a ^ t;
	d <<= d;
	printf("a ^ t: %i %d %f %o\n",d,d,d,d);
	d = c ^ g;
	d <<= d;
	printf("c ^ g: %i %d %f %o\n",d,d,d,d);
	d = g ^ t;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = a ^ hyphen;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = g ^ hyphen;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = t ^ hyphen;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = c ^ hyphen;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = a ^ asterisk;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = g ^ asterisk;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = t ^ asterisk;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
	d = c ^ asterisk;
	d <<= d;
	printf("g ^ t: %i %d %f %o\n",d,d,d,d);
}

void c_kimura_distance(char* str1, char* str2, float critical_value,
						float cutoff, int shortcut_type, float gene_length,
						SV* sv_transitions, SV* sv_transversions,
						SV* sv_comparisons, SV* sv_K2P, SV* sv_variance,
						SV* sv_stderror, SV* sv_mink2p, SV* sv_maxk2p) {
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
	float transitions = 0;
	float transversions = 0;
	float comparisons = 0;
	short int i, base_pair = 0;
	float 	P, Q, w1, w2, c1, c2, c3,
			variance, stderror, K2P, 
			mink2p, maxk2p = 0;
	for(i = 0; i < strlen(str1); i++) {
		if(str1[i] == str2[i]) {
			// Skip identicals.
			comparisons++;
		} else {
			base_pair = str1[i] ^ str2[i];
			// Transitions
			if(base_pair 		== 6) { // A->G
				transitions++;
				comparisons++;
			}
			else if(base_pair 	== 23){ // C->T	
				transitions++;
				comparisons++;
			} // Transversions
			else if(base_pair 	== 21){ // A->T
					transversions++;
					comparisons++;
			}
			else if(base_pair 	== 2){ 	// A->C
					transversions++;
					comparisons++;
			}
			else if(base_pair 	== 4){ 	// C->G
					transversions++;
					comparisons++;
			}
			else if(base_pair 	== 19){ // C->T
					transversions++;
					comparisons++;
			}
		}
		if((i % 50) == 0) {
			// Compute Kimura Distance
			P = transitions/gene_length;
			Q = transversions/gene_length;
			w1 = 1-2*P-Q;
			w2 = 1-2*Q;
			c1 = 1/w1;
			c2 = 1/w2;
			c3 = (c1+c2)/2;
			variance = pow(c1,2)*P+pow(c3,2)*Q-pow((c1*P+c3*Q),2);
			stderror = sqrt(variance)/sqrt(gene_length);
			K2P = -log(w1)/2-log(w2)/4;
			// Compute minimum and maximum distance using normal distribution.
			mink2p = K2P - critical_value*stderror;
			if(mink2p < 0) {
				mink2p = 0;
			}
			// maxk2p = K2P + critical_value*stderror;
			
			if(shortcut_type == 1) { // check strict K2P, Delete
				if(K2P > cutoff)
					goto cutoff_shortcut;
			} else if(shortcut_type ==2) { // Check mink2p, Search
				if(mink2p > cutoff)
					goto cutoff_shortcut;
			}
		}
	}
	
	// Compute Kimura Distance
	P = transitions/comparisons;
	Q = transversions/comparisons;
	w1 = 1-2*P-Q;
	w2 = 1-2*Q;
	c1 = 1/w1;
	c2 = 1/w2;
	c3 = (c1+c2)/2;
	variance = pow(c1,2)*P+pow(c3,2)*Q-pow((c1*P+c3*Q),2);
	stderror = sqrt(variance)/sqrt(comparisons);
	K2P = -log(w1)/2-log(w2)/4;
	// Compute minimum and maximum distance using normal distribution.
	mink2p = K2P - critical_value*stderror;
	if(mink2p < 0) {
		mink2p = 0;
	}
	maxk2p = K2P + critical_value*stderror;

	cutoff_shortcut:

	// Modify variables.
	sv_setnv(sv_transitions,transitions);
	sv_setnv(sv_transversions,transversions);
	sv_setnv(sv_comparisons,comparisons);
	sv_setnv(sv_K2P,K2P);
	sv_setnv(sv_variance,variance);
	sv_setnv(sv_stderror,stderror);
	sv_setnv(sv_mink2p,mink2p);
	sv_setnv(sv_maxk2p,maxk2p);
}


void c_k2p_bootstrap(char* str1, char* str2, float critical_value,
						float cutoff, int shortcut_type, float gene_length,
						SV* sv_transitions, SV* sv_transversions,
						SV* sv_comparisons, SV* sv_K2P, SV* sv_variance,
						SV* sv_stderror, SV* sv_mink2p, SV* sv_maxk2p,
						char* bs_weights) {
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
	float transitions = 0;
	float transversions = 0;
	float comparisons = 0;
	short int i, base_pair = 0;
	float 	P, Q, w1, w2, c1, c2, c3,
			variance, stderror, K2P, 
			mink2p, maxk2p = 0;
	for(i = 0; i < strlen(str1); i++) {
		if(str1[i] == str2[i]) {
			// Skip identicals.
			comparisons += bs_weights[i];
		} else {
			base_pair = str1[i] ^ str2[i];
			// Transitions
			if(base_pair 		== 6) { // A->G
				transitions += bs_weights[i];
				comparisons += bs_weights[i];
			}
			else if(base_pair 	== 23){ // C->T	
				transitions += bs_weights[i];
				comparisons += bs_weights[i];
			} // Transversions
			else if(base_pair 	== 21){ // A->T
					transversions += bs_weights[i];
					comparisons += bs_weights[i];
			}
			else if(base_pair 	== 2){ 	// A->C
					transversions += bs_weights[i];
					comparisons += bs_weights[i];
			}
			else if(base_pair 	== 4){ 	// C->G
					transversions += bs_weights[i];
					comparisons += bs_weights[i];
			}
			else if(base_pair 	== 19){ // C->T
					transversions += bs_weights[i];
					comparisons += bs_weights[i];
			}
		}
		if((i % 50) == 0) {
			// Compute Kimura Distance
			P = transitions/gene_length;
			Q = transversions/gene_length;
			w1 = 1-2*P-Q;
			w2 = 1-2*Q;
			c1 = 1/w1;
			c2 = 1/w2;
			c3 = (c1+c2)/2;
			variance = pow(c1,2)*P+pow(c3,2)*Q-pow((c1*P+c3*Q),2);
			stderror = sqrt(variance)/sqrt(gene_length);
			K2P = -log(w1)/2-log(w2)/4;
			// Compute minimum and maximum distance using normal distribution.
			mink2p = K2P - critical_value*stderror;
			if(mink2p < 0) {
				mink2p = 0;
			}
			// maxk2p = K2P + critical_value*stderror;
			
			if(shortcut_type == 1) { // check strict K2P, Delete
				if(K2P > cutoff)
					goto cutoff_shortcut;
			} else if(shortcut_type ==2) { // Check mink2p, Search
				if(mink2p > cutoff)
					goto cutoff_shortcut;
			}
		}
	}
	
	// Compute Kimura Distance
	P = transitions/comparisons;
	Q = transversions/comparisons;
	w1 = 1-2*P-Q;
	w2 = 1-2*Q;
	c1 = 1/w1;
	c2 = 1/w2;
	c3 = (c1+c2)/2;
	variance = pow(c1,2)*P+pow(c3,2)*Q-pow((c1*P+c3*Q),2);
	stderror = sqrt(variance)/sqrt(comparisons);
	K2P = -log(w1)/2-log(w2)/4;
	// Compute minimum and maximum distance using normal distribution.
	mink2p = K2P - critical_value*stderror;
	if(mink2p < 0) {
		mink2p = 0;
	}
	maxk2p = K2P + critical_value*stderror;

	cutoff_shortcut:

	// Modify variables.
	sv_setnv(sv_transitions,transitions);
	sv_setnv(sv_transversions,transversions);
	sv_setnv(sv_comparisons,comparisons);
	sv_setnv(sv_K2P,K2P);
	sv_setnv(sv_variance,variance);
	sv_setnv(sv_stderror,stderror);
	sv_setnv(sv_mink2p,mink2p);
	sv_setnv(sv_maxk2p,maxk2p);
}