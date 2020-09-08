#!/usr/bin/env perl
#
# Wednesday, 24 May 2017, 07:41PM
# Script to validate GFF3

#################################################################################
#################################################################################

# validate_gff3.pl is released under the MIT license.

# MIT License

# Copyright (c) 2017 Gemy George Kaithakottil

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
# AUTHOR: Gemy George Kaithakottil (Gemy.Kaithakottil@gmail.com

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
my $prog = basename($0);
my $usage = "
	Script to validate GFF3

	Usage: $prog <annotation.gff3[FILE]>

	Required:
	Positional arguments:
	<annotation.gff3[FILE]>   - path/to/annotation.gff3 || annotation.gff3

	IMPORTANT:
	Only cared about gene,mRNA,exon,CDS,five_prime_UTR,three_prime_UTR (case-sensitive) features,
	mainly for protein coding transcripts.

Contact: Gemy.Kaithakottil\@gmail.com
";
my $trans_annotation = shift or die $usage;

# Process GFF3
open(IN, "< $trans_annotation") or die "ERROR: Could not open the file $trans_annotation\n";
my %gene_info=();
my %mrna_info=();
my %mrna_gene_info=();
my %exon_info=();
my %exon_coords_info=();
my %cds_info=();
my %cds_coords_info=();
my %cds_length_info=();
my %utr5_info=();
my %utr3_info=();
while(defined(my $line = <IN>)){
	next if($line =~ /^\#/);
	next unless ($line =~ /\S/); # \S matches non-whitespace.  If not found in $_, skip to next line.
	chomp $line;
	my @f = split("\t", $line);
	my $seqid = $f[0];
	my $type = $f[2];
	next unless ($type =~ /^(gene|mRNA|exon|CDS|five_prime_UTR|three_prime_UTR)$/);
	my $start = $f[3];
	my $end = $f[4];
	warn "WARN: start coordinate is greater than end coordinate [$start > $end], swapping them" if ($start > $end);
	($start,$end) = ($end,$start) if ($start > $end);
	my ($strand) = $f[6];
	warn "WARN: strand is not defined as + or -, assuming it is +, plus strand" unless ($strand =~ /^(\+|\-)$/);
	($strand) = ($strand =~ /^(\+|\-)$/) ? $strand : "+" ; # if a strand is a dot, then change it is a positive
	my ($phase) = ($f[7] =~ /^(0|1|2)$/) ? $f[7] : 0 ; # make sure phase is an integer 
	my $attrib = $f[8];
	if ($type =~ /^(gene)$/) {
		my ($gene) = $attrib =~ /ID\=([^\;\n]+)/ or die "Fatal error: Cannot parse gene ID field";
		if (exists $gene_info{$gene}) {
			die "Fatal error: You have duplicate gene id in the input file '$trans_annotation'\n";
		} else {
			$gene_info{$gene}{start} = $start;
			$gene_info{$gene}{end} = $end;
		}
	}
	elsif ($type =~ /^(mRNA)$/) {
		my ($id) = $attrib =~ /ID\=([^\;\n]+)/ or die "Fatal error: Cannot parse mRNA ID field";
		my ($parent) = $attrib =~ /Parent\=([^\;\n]+)/ or die "Fatal error: Cannot parse mRNA Parent field";
		if (exists $mrna_info{$id}) {
			die "Fatal error: You have duplicate mRNA id in the input file '$trans_annotation'\n";
		} else {
			$mrna_info{$id}{start} = $start;
			$mrna_info{$id}{end} = $end;
			$mrna_info{$id}{strand} = $strand;
			$mrna_info{$id}{gene} = $parent;
		}
		# work out gene span from mRNA line
		if (exists $mrna_gene_info{$parent}) {
			my $prev_start = $mrna_gene_info{$parent}{start};
			my $prev_end = $mrna_gene_info{$parent}{end};
			$mrna_gene_info{$parent}{start} = $start if ($start < $prev_start);
			$mrna_gene_info{$parent}{end} = $end if ($end > $prev_end);
		} else {
			$mrna_gene_info{$parent}{start} = $start;
			$mrna_gene_info{$parent}{end} = $end;
		}
	} elsif ($type =~ /^(exon)$/) {
		my ($id) = $attrib =~ /Parent\=([^\;\n]+)/ or die "Fatal error: Cannot parse exon Parent field";
		push(@{$exon_info{$id}}, [$start, $end]);
		# work out exon start and end
		if (exists $exon_coords_info{$id}) {
			my $prev_start = $exon_coords_info{$id}{start};
			my $prev_end = $exon_coords_info{$id}{end};
			$exon_coords_info{$id}{start} = $start if ($start < $prev_start);
			$exon_coords_info{$id}{end} = $end if ($end > $prev_end);
		} else {
			$exon_coords_info{$id}{start} = $start;
			$exon_coords_info{$id}{end} = $end;
		}
	} elsif ($type =~ /^(CDS)$/) {
		my ($id) = $attrib =~ /Parent\=([^\;\n]+)/ or die "Fatal error: Cannot parse CDS Parent field";
		push(@{$cds_info{$id}}, [$start, $end]);
		# work out cds start and end
		if (exists $cds_coords_info{$id}) {
			my $prev_start = $cds_coords_info{$id}{start};
			my $prev_end = $cds_coords_info{$id}{end};
			$cds_coords_info{$id}{start} = $start if ($start < $prev_start);
			$cds_coords_info{$id}{end} = $end if ($end > $prev_end);
		} else {
			$cds_coords_info{$id}{start} = $start;
			$cds_coords_info{$id}{end} = $end;
		}
		my $mod_end = $end - $phase; # deduct the phase from the end
		push(@{$cds_length_info{$id}}, [$start, $mod_end]); # now add to the hash cds_length_info when calculating CDS length
	} elsif ($type =~ /^(five_prime_UTR|three_prime_UTR)$/) {
		my ($id) = $attrib =~ /Parent\=([^\;\n]+)/ or die "Fatal error: Cannot parse UTR Parent field";
		push(@{$utr5_info{$id}}, [$start, $end]) if ($type =~ /^(five_prime_UTR)$/) ;
		push(@{$utr3_info{$id}}, [$start, $end]) if ($type =~ /^(three_prime_UTR)$/) ;
	}
}
close(IN);

# sort the utr exons
while (my $trans = each %utr5_info) {
	@{$utr5_info{$trans}} = sort {$a->[0] <=> $b->[0]} @{$utr5_info{$trans}};
}
# sort the utr exons
while (my $trans = each %utr3_info) {
	@{$utr3_info{$trans}} = sort {$a->[0] <=> $b->[0]} @{$utr3_info{$trans}};
}
# sort the exons
while (my $trans = each %exon_info) {
	@{$exon_info{$trans}} = sort {$a->[0] <=> $b->[0]} @{$exon_info{$trans}};
}
# sort the CDSs
while (my $trans = each %cds_info) {
	@{$cds_info{$trans}} = sort {$a->[0] <=> $b->[0]} @{$cds_info{$trans}};
}
# sort the CDSs
while (my $trans = each %cds_length_info) {
	@{$cds_length_info{$trans}} = sort {$a->[0] <=> $b->[0]} @{$cds_length_info{$trans}};
}

foreach my $gene (sort keys %mrna_gene_info) {
	if(exists $gene_info{$gene}) {
		my $actual_gene_start = $gene_info{$gene}{start};
		my $actual_gene_end = $gene_info{$gene}{end};
		my $computed_gene_start = $mrna_gene_info{$gene}{start};
		my $computed_gene_end = $mrna_gene_info{$gene}{end};
		warn "Fatal error: Gene start is not consistant to all mRNA spans for gene '$gene'" unless ($actual_gene_start == $computed_gene_start);
		warn "Fatal error: Gene end is not consistant to all mRNA spans for gene '$gene'" unless ($actual_gene_end == $computed_gene_end);
	} else {
		warn "Fatal error: Gene id '$gene' not found in the input '$trans_annotation'";
	}
}

foreach my $trans (sort keys %mrna_info) {
	my $strand = $mrna_info{$trans}{strand};
	# these are genomic coords
	my $actual_cds_start = undef;
	my $actual_cds_end = undef;
	if (exists $cds_coords_info{$trans}) {
		$actual_cds_start = $cds_coords_info{$trans}{start};
		$actual_cds_end = $cds_coords_info{$trans}{end};
	} else {
		$actual_cds_start = 0;
		$actual_cds_end = 0;
	}

	# GET EXON INFO
	my $exon_count=0;
	my $exon_len=0;
	my $exon_start=0;
	my $exon_end=0;
	my $exon_start_boo=0;
	my $exon_end_boo=0;
	if (exists $exon_info{$trans}) {
		my @exons = @{$exon_info{$trans}};
		foreach my $coords (@exons) {
			my $start = $coords->[0];
			my $end = $coords->[1];
			my $span = (($end - $start) + 1);
			$exon_count++;
			$exon_len += $span;
			if ($exon_start_boo) {
				$exon_start = $start if ($start < $exon_start);
			} else {
				$exon_start = $start;
				$exon_start_boo=1;
			}
			if ($exon_end_boo) {
				$exon_end = $end if ($end > $exon_end);
			} else {
				$exon_end = $end;
				$exon_end_boo=1;
			}
		}
	} else { # if transcript not present then give zeros
		$exon_count = 0;
		$exon_len = 0;
		$exon_start = 0;
		$exon_end = 0;
	}
	$mrna_info{$trans}{exon_count} = $exon_count;
	$mrna_info{$trans}{exon_len} = $exon_len;
	$mrna_info{$trans}{exon_start} = $exon_start;
	$mrna_info{$trans}{exon_end} = $exon_end;

	# GET CDS INFO
	my $cds_len_standard=0;
	my $cds_count=0;
	my $cds_complete_match=0;
	my $cds_partial_match=0;
	if (exists $cds_info{$trans}) {
		my @cdss = @{$cds_info{$trans}};
		my $cds_counts = scalar @cdss;
		foreach my $coords (@cdss) {
			my $start = $coords->[0];
			my $end = $coords->[1];
			my $span = (($end - $start) + 1);
			$cds_len_standard += $span;
			$cds_count++;
			# check how many are complete and partial
			if (exists $exon_info{$trans}) {
				my @exons = @{$exon_info{$trans}};
				foreach my $e_coords (@exons) {
					my $e_start = $e_coords->[0];
					my $e_end = $e_coords->[1];
					if ($e_start == $start && $end == $e_end) { # these are complete matches
						$cds_complete_match++;
					}
					elsif ( ($e_start < $start && $end == $e_end) || ($e_start == $start && $end < $e_end) || ($e_start < $start && $end < $e_end)) { # these are patrial matches
						$cds_partial_match++;
					} 
				}
			}
		}
		my $total_cds_count = $cds_complete_match + $cds_partial_match; # total matches should match the cds count
		warn "Fatal error: Not all CDS exons for transcript '$trans' are covered within exon." unless ($total_cds_count == $cds_counts);
	} else { # if transcript not present then give zeros
		$cds_count = 0;
		$cds_len_standard = 0;
		$cds_complete_match=0;
		$cds_partial_match=0;
	}
	$mrna_info{$trans}{cds_len_standard} = $cds_len_standard;
	$mrna_info{$trans}{cds_count} = $cds_count;
	$mrna_info{$trans}{cds_complete_match} = $cds_complete_match;
	$mrna_info{$trans}{cds_partial_match} = $cds_partial_match;

	# GET CDS SPAN INFO AFTER EXCLUDING THE PHASE INFORMATION
	my $cds_len=0;
	if (exists $cds_length_info{$trans}) {
		my @cdss = @{$cds_length_info{$trans}};
		foreach my $coords (@cdss) {
			my $start = $coords->[0];
			my $end = $coords->[1];
			my $span = (($end - $start) + 1);
			$cds_len += $span;
		}
	} else { # if transcript not present then give zeros
		$cds_len = 0;
	}
	$mrna_info{$trans}{cds_len} = $cds_len;

	# GET UTR5 INFO
	my $utr5_count=0;
	my $utr5_len=0;
	my $utr5_complete_match=0;
	my $utr5_partial_match=0;	
	if (exists $utr5_info{$trans}) {
		my @utr5s = @{$utr5_info{$trans}};
		my $utr5s_count= scalar @utr5s;
		foreach my $coords (@utr5s) {
			my $start = $coords->[0];
			my $end = $coords->[1];
			my $span = (($end - $start) + 1);
			$utr5_count++;
			$utr5_len += $span;
			# check how many are complete and partial
			if (exists $exon_info{$trans}) {
				my @exons = @{$exon_info{$trans}};
				foreach my $e_coords (@exons) {
					my $e_start = $e_coords->[0];
					my $e_end = $e_coords->[1];
					if ($e_start == $start && $end == $e_end) { # these are complete matches
						$utr5_complete_match++;
					}
					elsif ( ($strand eq "+") && ($e_start == $start && $end < $e_end) ) { # these are patrial matches
						$utr5_partial_match++;
						warn "Fatal error: utr5 start position not at expected start for transcript '$trans'" unless ( $end == ($actual_cds_start-1) );
					} 
					elsif ( ($strand eq "-") && ($e_start < $start && $end == $e_end)) { # these are patrial matches
						$utr5_partial_match++;
						warn "Fatal error: utr5 start position not at expected start for transcript '$trans'" unless ( $start == ($actual_cds_end+1) );
					} 
				}
			}
		}
		my $total_utr5_count = $utr5_complete_match + $utr5_partial_match; # total matches should match the utr5 count
		warn "Fatal error: Not all utr5 exons for transcript '$trans' are covered within exon." unless ($total_utr5_count == $utr5s_count);
	} 
	else { # if transcript not present then give zeros
		$utr5_count = 0;
		$utr5_len = 0;
		$utr5_complete_match=0;
		$utr5_partial_match=0;
	}

	# GET UTR3 INFO
	my $utr3_count=0;
	my $utr3_len=0;
	my $utr3_complete_match=0;
	my $utr3_partial_match=0;
	if (exists $utr3_info{$trans}) {
		my @utr3s = @{$utr3_info{$trans}};
		my $utr3s_count= scalar @utr3s;
		foreach my $coords (@utr3s) {
			my $start = $coords->[0];
			my $end = $coords->[1];
			my $span = (($end - $start) + 1);
			$utr3_count++;
			$utr3_len += $span;
			# check how many are complete and partial
			if (exists $exon_info{$trans}) {
				my @exons = @{$exon_info{$trans}};
				foreach my $e_coords (@exons) {
					my $e_start = $e_coords->[0];
					my $e_end = $e_coords->[1];
					if ($e_start == $start && $end == $e_end) { # these are complete matches
						$utr3_complete_match++;
					}
					elsif ( ($strand eq "+") && ($e_start < $start && $end == $e_end) ) { # these are patrial matches
						$utr3_partial_match++;
						warn "Fatal error: utr3 start position not at expected start for transcript '$trans'" unless ( $start == ($actual_cds_end+1) );
					} 
					elsif ( ($strand eq "-") && ($e_start == $start && $end < $e_end) ) { # these are patrial matches
						$utr3_partial_match++;
						warn "Fatal error: utr3 start position not at expected start for transcript '$trans'" unless ( $end == ($actual_cds_start-1) );
					}
				}
			}
		}
		my $total_utr3_count = $utr3_complete_match + $utr3_partial_match; # total matches should match the utr3 count
		warn "Fatal error: Not all utr3 exons for transcript '$trans' are covered within exon." unless ($total_utr3_count == $utr3s_count);
	} 
	else { # if transcript not present then give zeros
		$utr3_count = 0;
		$utr3_len = 0;
		$utr3_complete_match=0;
		$utr3_partial_match=0;	
	}
	# check to make sure that we have exon features for transcripts with CDS feature.
	my $mrna_trans_exon_len = $mrna_info{$trans}{exon_len};
	my $mrna_trans_cds_len = $mrna_info{$trans}{cds_len};
	my $mrna_trans_cds_len_std = $mrna_info{$trans}{cds_len_standard};
	my $cds_cdna_ratio = undef;
	if ( ($mrna_trans_exon_len == 0) && ( $mrna_trans_cds_len > 0 || $mrna_trans_cds_len_std > 0 ) ) {
		warn "Fatal error: Transcript '$trans' do not have exon feature, but has CDS feature";
		$cds_cdna_ratio = 0;
	} elsif ( ($mrna_trans_exon_len == 0) && ( $mrna_trans_cds_len == 0 || $mrna_trans_cds_len_std == 0 ) ) {
		warn "Fatal error: Transcript '$trans' do not have exon and CDS features";
		$cds_cdna_ratio = 0;
	} else {
		$cds_cdna_ratio = sprintf("%.2f",( ($mrna_trans_cds_len/$mrna_trans_exon_len) * 100) )
	}
	die "Fatal error: Transcript '$trans' cds to cDNA ratio could not be calculated." unless (defined $cds_cdna_ratio);
	warn "Fatal error: Transcript '$trans' CDS length is greater than exon length" if ( $mrna_trans_cds_len_std > $mrna_trans_exon_len );
	my $actual_trans_start = $mrna_info{$trans}{start};
	my $actual_trans_end = $mrna_info{$trans}{end};
	my $computed_trans_start = $mrna_info{$trans}{exon_start};
	my $computed_trans_end = $mrna_info{$trans}{exon_end};
	warn "Fatal error: Transcript start is not consistant to all exon spans for transcript '$trans'" unless ($actual_trans_start == $computed_trans_start);
	warn "Fatal error: Transcript end is not consistant to all mRNA spans for transcript '$trans'" unless ($actual_trans_end == $computed_trans_end);
}

exit;

