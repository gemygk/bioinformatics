#!/usr/bin/env perl
# 
#	01-May-2015
#	Script to rename fasta header 
# - Mainly to shorten fasta headers that are too long, to run tools like transposonPSI, interproscan, Blast2GO

#################################################################################
#################################################################################

# rename_fasta_header.pl is released under the MIT license.

# MIT License

# Copyright (c) 2015 Gemy George Kaithakottil

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
# AUTHOR: Gemy George Kaithakottil (Gemy.Kaithakottil@gmail.com)

use strict;
use warnings;
use File::Basename;
my $prog = basename($0);

my $usage = "
	$prog -- Script to rename fasta header based on prefix

	Usage: perl $prog <file.fasta[FILE]> <prefix[STRING]>

	Required:
	Positional arguments:
	<file.fasta[FILE]>   - path/to/file.fasta || file.fasta
	<prefix[STRING]>     - prefix name for fasta header
	                       Eg., if given prefix is \"peptide\", 
	                       then header format will be in format peptide_1,peptide_2,..,peptide_n


	NOTE:
	A file will be generated with name <file.fasta>.new-old.id.txt
	with the new id and old id, for future reference.
	Make sure that you have write access in the current folder.


Contact: Gemy.Kaithakottil\@gmail.com
";

my $input_file = shift or die $usage;
my $prefix = shift or die $usage;

my $input_file_basename = basename($input_file);

# To write the new id old id relationship
# First removing the file if it already exists.
if (-e "$input_file_basename.new-old.id.txt") {
	system("rm $input_file_basename.new-old.id.txt");
}
open(my $fh1, '>>', "$input_file_basename.new-old.id.txt") or die "Could not open file $input_file_basename.new-old.id.txt $!";

print $fh1 "#1.new_id\t#2.old_id\n";

my %id_hash=(); # to store the fasta ids
my $count=1;    # fasta counter
open(FILE,"< $input_file") or die "$!";
while (<FILE>) {      # Read lines
	# commenting below two lines so that input and output line counts are consistant.
	# next if /^\#/;    # Remove comments from $_.
	# next unless /\S/; # \S matches non-whitespace.  If not found in $_, skip to next line.
	chomp;
	if (/^\>/) {
	    my ($id) = $_ =~ /^>([^\s\n]+)/; # extract fasta id
	    die "Fatal error: Could not extract fasta id from fasta header from line '$_'" unless (defined $id);
		die "Fatal error: Input fasta file '$input_file' contains a duplicate id: $id\n" if (exists $id_hash{$id});
		$id_hash{$id} = 1;
		my $new_header = "$prefix\_$count";
		print $fh1 "$new_header\t$id\n";
		print ">$new_header\n";
		$count++;
	} else {
		print "$_\n";
	}
}
close $fh1;

print STDERR "#INFO: new id and old id relationship file created '$input_file_basename.new-old.id.txt'\n";

exit;