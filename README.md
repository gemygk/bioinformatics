# Bioinformatics scripts and utilities

# scripts/rename_fasta_header.pl -- Script to rename fasta header based on prefix

	Usage: rename_fasta_header.pl <file.fasta[FILE]> <prefix[STRING]>

	Required:
	Positional arguments:
	<file.fasta[FILE]>   - path/to/file.fasta || file.fasta
	<prefix[STRING]>     - prefix name for fasta header
	                       Eg., if given prefix is "peptide",
	                       then header format will be in format peptide_1,peptide_2,..,peptide_n


	NOTE:
	A file will be generated with name <file.fasta>.new-old.id.txt
	with the new id and old id, for future reference.
	Make sure that you have write access in the current folder.
	
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

