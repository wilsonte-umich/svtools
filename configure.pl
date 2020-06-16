#!/usr/bin/perl
use strict;
use warnings;

our $utilityDescription = "applications for exploring paired-end sequencing data for structural variants";
our @prerequisites = qw(
R
bwa
samtools
bedtools
cat
awk
gunzip
gzip
sort
bgzip
tabix);
require "./lib/configure.pl";

