#!/usr/bin/env perl
#
#	split_fasta_to_chunks.pl - Script to split fasta sequence file into user defined chunks

#################################################################################
#################################################################################

# split_fasta_to_chunks.pl is released under the MIT license.

# MIT License

# Copyright (c) 2014 Gemy George Kaithakottil

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

#################################################################################
#################################################################################

# 
# AUTHOR:Gemy George Kaithakottil (Gemy.Kaithakottil@gmail.com)

use strict;
use warnings;
use File::Basename;
my $prog = basename($0);

my $usage = "
	Script to split fasta sequence file into user defined chunks

	Usage: perl $prog <file.fasta[FILE]> <prefix[STRING]> <no_of_fasta_seqs[INT]>

	Required:
	Positional arguments:
	<file.fasta[FILE]>      - path/to/file.fasta || file.fasta
	<prefix[STRING]>        - prefix name for chunk
	                          Eg., if given prefix is \"chunk\", files created will be called chunk-1.txt,chunk-2.txt,..,chunk-n.txt
	<no_of_fasta_seqs[INT]> - indicate how many fasta sequences you want in a chunk
	                          Eg., if you want 1000 sequences per chunk then give 1000

	IMPORTANT:
	You should have installed LEAFF (http://kmer.sourceforge.net/wiki/index.php/LEAFF_User's_Guide) and 
	the executable leaff to be available in the PATH

Contact: Gemy.Kaithakottil\@gmail.com
";

my $file = shift or die $usage;
my $prefix = shift or die $usage;
my $num = shift or die $usage;
my $count=0;
unless(-e $file) {
	print "Fatal error: Cannot open $file:$!";
}

# Check that the environment have the executable variables
my @tools = qw (leaff);
foreach my $tool (@tools) {
	my $tool_path = qx(which $tool);
	chomp ($tool_path);
	die "\n## Fatal error: '$tool' executable is not available in the PATH. Please make sure you source the executable '$tool' before executing the script\n" unless ( $tool_path );
}

# get number of fasta sequences in file
my $lines = `grep -c '^>' $file`;
chomp($lines);
$num = $num-1;
for (my $i=0;$i<=($lines-1);$i++) {
	$count++;
	my $j=$i+$num; # like 9,19,29...
	if ($j <= $lines && ($j != $lines)) {
		print "Starting ",($i+1),"-",($j+1)," End\n";
		`leaff -f $file -6 70 -S $i $j > $prefix-$count.txt`;
		$i=$j;
	}
	else {
		$j = ($lines-1);
		print "Starting ",($i+1),"-",($j+1)," End**Last File**\n";
		`leaff -f $file -6 70 -S $i $j > $prefix-$count.txt`;
		last;
	}
}
print "Process completed!\n";

exit;
