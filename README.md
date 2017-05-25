# Bioinformatics scripts and utilities

# scripts/split_fasta_to_chunks.pl - Script to split fasta sequence file into user defined chunks

	Usage: perl split_fasta_to_chunks.pl <file.fasta[FILE]> <prefix[STRING]> <no_of_fasta_seqs[INT]>

	Required:
	Positional arguments:
	<file.fasta[FILE]>      - path/to/file.fasta || file.fasta
	<prefix[STRING]>        - prefix name for chunk
	                          Eg., if given prefix is "chunk", files created will be called chunk-1.txt,chunk-2.txt,..,chunk-n.txt
	<no_of_fasta_seqs[INT]> - indicate how many fasta sequences you want in a chunk
	                          Eg., if you want 1000 sequences per chunk then give 1000

	IMPORTANT:
	You should have installed LEAFF (http://kmer.sourceforge.net/wiki/index.php/LEAFF_User's_Guide) and
	the executable leaff to be available in the PATH
