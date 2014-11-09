Bioinformatics Toolbox (BTBox)
A multi-purpose collection of bioinformatics tools written in Perl, R, and C.

Contributors: Bryan P. White

Citations:
White, B.P., Pilgrim E.M., Boykin L.M., Stein E.D., Mazor R.D. 2014. Comparing 
four species delimitation methods applied to a DNA barcode data set of insect 
larvae for use in routine bioassessment. Freshwater Science 33(1), 338-348.

Commonly used commands:

Genbank Downloader: seq_convert_genbank.pl

	Download all COI sequences for Gastropoda from NCBI:
	seq_convert_genbank.pl -query COI_full -batch-cap 500 -term Gastropoda

	Only download sequences with a voucher ID:
	seq_convert_genbank.pl -query COI_full -voucher-only 1 -batch-cap 500 
	-term Gastropoda

	Use a list of taxa:
	seq_convert_genbank.pl -list Gastropoda_list.txt -query COI_full 
	-voucher-only 1 -batch-cap 500

	Include pubmed information such as abstracts:
	seq_convert_genbank.pl -list Gastropoda_list.txt -query ND2_full 
	-batch-cap 500 -pubmed 1 -outp OutputFile

Process genbank files: seq_process_genbank.pl
*Requires MAFFT to be in the path already

	Automatically quality check sequences downloaded using the Genbank Downloader
	seq_process_genbank.pl -gb Gastropoda_COI.csv -out Gastropoda_COI_qc 
	-match COI_Match.fas -otu-cutoff 0.01

	Use a specific match sequence, use 3 threads during MAFFT alignment
	seq_process_genbank.pl -gb Gastropoda_COI.csv -out Gastropoda_COI_qc 
	-match COI_Gastropoda_Match.fas -otu-cutoff 0.01 -threads 3

Clustering/OTU Delimitation: dnab_otu_delim.pl

	Delimit clusters at a 2% genetic distance cutoff (Kimura 2-parameter)
	dnab_otu_delim.pl -aln1 sample_baetis_seqs.fas -cutoff 0.02

	Skip calculating intra-OTU pairwise distances
	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 
	-aln1 sample_baetis_seqs.fas

	Just count OTU's, skip intra dist, randomly splice the alignment
	dnab_otu_delim.pl -aln1 sample_baetis_seqs.fas -shortcut-freq 0.05 
	-ran-splice 1 -skip-intra-dist 1 -pseudo-reps 1

	Various bootstrapping methods
	dnab_otu_delim.pl -shortcut-freq 0.05 -ran-splice 1 -skip-intra-dist 1 
	-pseudo-reps 1 -aln1 sample_baetis_seqs.fas

	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -bootstrap 1 
	-bootstrap-size 500 -pseudo-reps 10 -ran-splice 1 -aln1 sample_baetis_seqs.fas

	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -bootstrap 1 
	-bootstrap-size 500 -pseudo-reps 100 -ran-splice 1 -aln1 sample_baetis_seqs.fas

	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -bootstrap 1 
	-bootstrap-size 500 -pseudo-reps 2000 -specific-splice 1:50 
	-aln1 sample_baetis_seqs.fas

	Pseudo-repping correspondence:
	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -bootstrap 1
	 -bootstrap-size 500 -pseudo-reps 100 -skip-nn 1 -aln1 sample_baetis_seqs.fas

	Pseudo-repping splicing:
	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -bootstrap 1 
	-bootstrap-size 500 -pseudo-reps 1000 -skip-nn 1 -min-aln-length 25 
	-ran-splice 1 -aln1 sample_baetis_seqs.fas

	Specific splice for primer:
	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -bootstrap 1 
	-bootstrap-size 500 -pseudo-reps 1000 -skip-nn 1 -min-aln-length 25 
	-specific-splice 1:135 -aln1 sample_baetis_seqs.fas

	Printing short-read simulation (short splice)
	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -skip-nn 1 
	-min-aln-length 654 -aln1 sample_baetis_seqs.fas -print-spliced-aln 1 
	-spliced-aln-size 400 -print-ref-seq 0
	
	dnab_otu_delim.pl -shortcut-freq 0.05 -skip-intra-dist 1 -skip-nn 1 
	-aln1 sample_baetis_seqs.fas -print-spliced-aln 1 -spliced-aln-size 135 
	-print-ref-seq 1

454 Pipeline
*Requires NCBI Standalone BLAST+ to be in the path already

	Parse raw 454 output
	extract_fasta_454.pl -aln 454_output -out 454_output.fas -clean 1

	Create a standalone BLAST database to check sequences against
	makeblastdb -in metazoan_db.fas -dbtype nucl

	BLAST uknown sequences against a BLAST database
	seq_local_batch_blast.pl -aln 454_output.fas -out 454_output_labeled.fas 
	-db metazoan_db.fas -max_target_seqs 10

	Summarize results
	seq_output_taxa.pl -aln TD90_1_test1.fas -out TD90_1_test1.csv


