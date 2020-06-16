
--------------------------------------------------------------------------------
OVERVIEW
--------------------------------------------------------------------------------

svtools explores paired-end/mate-pair high throughput sequencing data for
genomic structural variants, including CNVs as well as inversions and
translocations. svtools implements the analysis logic of VAMP as described here:

  Arlt MF, Ozdemir AC, Birkeland SR, Lyons RH Jr, Glover TW, Wilson TE.
  Comparison of constitutional and replication stress-induced genome structural
  variation by SNP array and mate-pair sequencing. Genetics. 2011 187:675-83.
  
  Birkeland SR, Jin N, Ozdemir AC, Lyons RH Jr, Weisman LS, Wilson TE.
  Discovery of mutations in Saccharomyces cerevisiae by pooled linkage analysis
  and whole-genome sequencing. Genetics. 2010 186:1127-37.

with the following major differences:

1)  svtools does not perform mapping, i.e. read alignment; it takes as input
    SAM/BAM format data created by bwa mem (or possibly other alignment tools).

2)  svtools makes use of long read formats, especially split and clipped
    reads, as a supplement to read pair information.
 
3)  svtools uses only disk files, not an Oracle database.

4)  svtools supports multi-sample comparisons for finding SVs/CNVs unique to
    just one sample.
    
5)  svtools offers multi-sample read depth analysis as a supplement
    to read-pair and split-read information.

--------------------------------------------------------------------------------
USAGE PHASES
-------------------------------------------------------------------------------- 

svtools is a suite of scripts that manage data flow through a series of
open-source third-party programs and applications. 

Command line tools run under a Linux or compatible environment are used to:
- collect anomalous read pairs that nominate structural variant junctions
- create read depth, i.e. coverage, maps using variable-width bins
- call structural variants from one or many samples
  
--------------------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------------------

Analysis Phase

System prerequisites for running svtools in the analysis phase are:

Perl:      http://www.perl.org/

R:         http://www.r-project.org/

bwa:       http://bio-bwa.sourceforge.net/

Samtools:  http://samtools.sourceforge.net/

Bedtools:  http://code.google.com/p/bedtools/

in addition to other standard Linux system utilities. The installation 
utility will apprise you if you are missing a prerequisite. Note that svtools 
cannot be run on Windows in the analysis phase.

Visualization Phase

System prerequisites for running svtools in the analysis phase are:

Perl:      http://www.perl.org/

R:         http://www.r-project.org/

in addition to specific perl packages listed in INSTALL. 

--------------------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------------------

Read the INSTALL file to learn how to use the 'configure.pl' installation script.

--------------------------------------------------------------------------------
ANALYSIS COMMAND SUMMARY
--------------------------------------------------------------------------------

svtools encompasses a series of commands called as follows:

    svtools <command> [options]

svtools commands are as follows:

extract     search SAM/BAM data for read pairs that nominate structural variants
index       created a tabix-indexed bed file of extracted anomalous pairs
find        analyze extracted read pairs for recurring structural variants

coverage    use multiple samples to establish a read depth profile of all samples
    
Use:

    svtools <command> -h/--help

for more information on individual commands and their options.

--------------------------------------------------------------------------------
INPUT DATA
--------------------------------------------------------------------------------

Input data can be provided on STDIN, or using option -f/--file.

extract

'svtools extract' takes as input BAM/SAM format read-pair alignments, sorted by
name. All input data must come from one sequencing library so that only a single
fragment size distribution is expected. Thus, each library must be first be
subjected individually to 'svtools extract' (one sample might have many libraries).

index

'svtools index' takes as input file(s) created by 'svtools extract'. Use the Linux
'cat' command to concatenate multiple extract files into a single indexed bed file.

find

'svtools find' takes as input file(s) created by 'svtools extract'. Multiple
different libraries/samples can be combined, which is why 'extract' and 'find' are
separate commands. 'svtools find' also requires a genome fasta file.

coverage


--------------------------------------------------------------------------------
USAGE EXAMPLES
--------------------------------------------------------------------------------

extract

Place 'svtools extract' in the data stream immediately following read alignment
(note that providing a '.gz' extension results in a gzipped output file):

bwa mem <options and arguments> |
svtools extract -s my_sample -l library_1 -o /path/my_sample.library_1.extract.svt.gz -r |
samtools -Sb - > /path/my_sample.library_1.bam

Alternatively, run 'svtools extract' on an existing coordinate-sorted BAM file:

svtools extract -bx -f /path/my_sample.bam \
                -s my_sample -l library_1 -o /path/my_sample.library_1.extract.svt.gz


index

Make and indexed BED file from multiple 'extract' files:

gunzip -c /path/*.extract.svt.gz |
svtools index -o /path/all_samples.indexed.bed.bgz


find

Call SVs in a single sample:

svtools find -g /path/hg19.fa -s my_sample \
             -f /path/my_sample.extract.svt.gz \
             -o /path/my_sample.find.svt.gz

Merge multiple 'extract' files and submit to SV calling (note that when input is
provided on stdin, svtools will not know to uncompress gzipped data, you must do
it using gunzip or zcat):

gunzip -c /path/*.extract.svt.gz |
svtools find -g /path/hg19.fa -s all_samples -o /path/all_samples.find.svt.gz

