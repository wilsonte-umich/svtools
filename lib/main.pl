use strict;
use warnings;
use Cwd(qw(abs_path));
			
#========================================================================
# 'main.pl' is the command-line interpreter and help generator
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
# global program information
#------------------------------------------------------------------------
use vars qw($utility $libDir);
require "$libDir/command.pl";
require "$libDir/file.pl";
require "$libDir/format.pl";
require "$libDir/common.pl";
our $error = "error";
our ($command, @args) = @ARGV;
$| = 1;
#------------------------------------------------------------------------
# commands and options
#------------------------------------------------------------------------
our %options; # collects the actual options specified by user on command line
#------------------------------------------------------------------------
#my @file    = ("f", "<str>", "path to the input file(s) [STDIN]");
my @bam     = ("b", undef,   "input data is in BAM format [SAM]");
my @tmpDir  = ("T", "<str>", "temporary directory (must exist) [/tmp]");
my @maxMem  = ("m", "<int>", "maximum RAM bytes to use when sorting [1000000000]");
my @pdfDir  = ("f", "<str>", "directory where output from svtools pdf can be found");
my @sample  = ("s", "<str>", "sample name for the data");
my @library = ("l", "<str>", "name of sequencing library (a sample may have many libraries)");
my @group   = ("G", "<str>", "group name for set of input samples/libraries");
my @exclude = ("x", "<str>", "path to a BED file containing genomic regions to exclude");
my @skipChrom=("k", "<str>", "comma-delimited list of chromosomes to skip in CNV calling [none]");
my @output  = ("o", "<str>", "output directory where files will be placed (must exist)");
my @outFile = ("o", "<str>", "output file path");
my @input   = ("i", "<str>", "input directory containing files from previous command");
my @plotDir = ("P", "<str>", "target directory for summary plots (must exist) [plots not created]");
my @repeat  = ("r", undef,   "repeat all input data to STDOUT (requires --output/-o)");
my @nStatPairs = ("n", "<int>", "number of read pairs used to determine insert sizes [1e5]");
my @minTLen = ("l", "<str>", "min TLEN as maxTLen+z*(maxTlen-medTlen), del/dup/inv [1.5/-1/0]");
my @minPair = ("p", "<int>", "only report variants with >=p crossing fragments [4]");
my @unique  = ("u", undef,   "only report variants called in exactly one sample");
my @strict  = ("U", undef,   "only report variants with matching pairs in one sample");
my @junction= ("j", undef,   "only report variants with sequenced junctions");
my @minSize = ("z", "<int>", "only report variants >= z bp on the same chromosome [0]");
my @maxSize = ("Z", "<int>", "only report variants <= Z bp on the same chromosome [1000000000]");
my @chrom   = ("c", "<str>", "chromosome to analyze in this execution");
my @chroms  = ("C", "<str>", "comma-delimited list of all chromosomes to be analyzed");
my @delete  = ("d", undef,   "delete the bin/find source files after processing");
my @minMapQ = ("q", "<int>", "minimum MAPQ of one read for a pair to be considered [5]");
my @gapFile = ("g", "<str>", "BED file of gaps in chromosomes");
my @minMapQ10X = ("q", "<int>", "minimum MAPQ of index read of a pair [30]");
my @binWeight = ("W", "<str>", "average number of fragments per bin per cell [10]");
#------------------------------------------------------------------------
our %commands =  ( # 0=allowed, 1=required
#------------------------------------------------------------------------
# prepare PDF of FR pairs, used in Pairs and MLE approaches
#------------------------------------------------------------------------
    pdf => [
        \&svtools_pdf,
        "calculate the probability density function of FR pairs for predicting fragment sizes", {
		'sample'=>          [1, 1, @sample, "Input Options"],
		'exclude-file'=>    [0, 2, @exclude],         
        'output-dir'=>      [1, 3, @output, "Output Options"],
        'plot-dir'=>        [0, 4, @plotDir],
        'min-qual'=>        [0, 5, @minMapQ, "PDF Options"],
        'max-tLen'=>        [0, 6, "Z", "<int>", "maximum TLEN in bp used to determine library distributions [20000]"],
        'tLen-bin-size'=>   [0, 7, "b", "<int>", "TLEN bin size in bp, i.e. resolution of subsequent calls [200]"],
        'n-stat-pairs'=>    [0, 8, @nStatPairs],
   }],
#------------------------------------------------------------------------
# 'Pairs' approach, i.e. find clusters of anomalous fragments that predict CNVs
#------------------------------------------------------------------------
    extract =>[
        \&svtools_extract,
        "search SAM/BAM data for unique read pairs/fragments that nominate structural variants", {
		'sample'=>          [1, 1, @sample, "Input Options"],
        'library'=>         [1, 2, @library],        
        'bam'=>             [0, 3, @bam],
        'coordinate'=>      [0, 4, "C", undef, "input data are sorted by coordinate (svtools will re-sort)"],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],
        'output-dir'=>      [1, 7, @output, "Output Options"],
        'plot-dir'=>        [0, 8, @plotDir],
        'repeat'=>          [0, 9, @repeat],
        'n-stat-pairs'=>    [0, 10, @nStatPairs, "Extract Options"],
		'exclude-file'=>    [0, 11, @exclude],
		'min-clip'=>        [0, 12, "c", "<int>", "minimum clip length taken as evidence of a junction [8]"],
        'min-qual'=>        [0, 13, @minMapQ],
    }],
    mask =>[
        \&svtools_mask,
        "create a filter to subtract chimeric fragments with observed reference fragment end[s]", {
		'sample'=>          [1, 1, @sample, "Input Options"],
        'library'=>         [1, 2, @library],        
        'bam'=>             [0, 3, @bam],        
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],
        'output-dir'=>      [1, 7, @output, "Output Options"],
        'repeat'=>          [0, 9, @repeat],
		'exclude-file'=>    [0, 11, @exclude, "Mask Options"],
		'min-qual'=>        [0, 13, @minMapQ],
    }],
    find =>[
        \&svtools_find,
        "analyze extracted anomalous pairs for recurring structural variants", {
		'group'=>           [1, 1, @group, "Input Options"],   
        'genome-fasta'=>    [1, 2, "g", "<str>", "path to the genome FASTA file, e.g. /path/hg19.fa"],
        'output-dir'=>      [1, 3, @output, "Output Options"],       
        'tmp-dir'=>         [0, 4, @tmpDir],
        'max-mem'=>         [0, 5, @maxMem],        
		'min-tLen'=>        [0, 6, @minTLen, "Find Options"],
        'min-set-pairs'=>   [0, 7, @minPair],
        'unique'=>          [0, 8, @unique],
        'unique-strict'=>   [0, 9, @strict],
        'jxn-only'=>        [0, 10, @junction],
        'min-size'=>        [0, 11, @minSize],
        'max-size'=>        [0, 12, @maxSize],
        'chrom'=>           [1, 13, @chrom],
		'chroms'=>          [1, 14, @chroms],
        'use-mask'=>        [0, 15, "M", undef, "use mask to suppress anomalies that share an observed proper end"],
    }],
    cat =>[
        \&svtools_cat,
        "concatenate and tabix the find files for a list of chromosomes", {
		'group'=>           [1, 1, @group, "Input Options"],
		'input-dir'=>       [1, 2, @input],
		'delete'=>          [0, 3, @delete],
        'output-dir'=>      [1, 4, @output, "Output Options"],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem], 
		'chroms'=>          [1, 14, @chroms, "Find Options"],
    }],
    index =>[
        \&svtools_index,
        "create a tabix-index of extracted anomalous pairs for use in browser", {
		'group'=>           [1, 1, @group, "Input Options"],
        'output-dir'=>      [1, 2, @output, "Output Options"],    
        'max-mem'=>         [0, 3, @maxMem],
        'min-tLen'=>        [0, 4, @minTLen],
    }],
    pull =>[
        \&svtools_pull,
        "pull the sequences around a CNV junction (or any genome region)", {
        'genome-fasta'=>    [1, 1, "g", "<str>", "path to the genome FASTA file, e.g. /path/hg19.fa"],
        'padding'=>         [0, 2, "p", "<int>", "number of bases to add to the ends of the region [0]"],
        'line-length'=>     [0, 3, "l", "<int>", "number of characters per output line [100]"],
        'clipped-only'=>    [0, 4, "c", undef,   "only return clipped reads that might sequence a junction"],
        'min-qual'=>        [0, 5, "q", "<int>", "minimum MAPQ of reads to be pulled [0]"],
    }],
#------------------------------------------------------------------------
# 'Depth' approach, i.e. find genomic spans with altered fragment coverage
#------------------------------------------------------------------------		  
    map =>[
        \&svtools_map,
        "create a coverage map of each chromosome in a sorted input BAM stream/file", {
		'sample'=>          [1, 1, @sample, "Input Options"],
        'library'=>         [1, 2, @library],       
        'bam'=>             [0, 3, @bam],
        'tmp-dir'=>         [0, 4, @tmpDir],
        'output-dir'=>      [1, 5, @output, "Output Options"],
        'repeat'=>          [0, 6, @repeat],
		'n-stat-pairs'=>    [0, 7, @nStatPairs, "Map Options"],
		'exclude-file'=>    [0, 8, @exclude],
		'min-qual'=>        [0, 9, "q", "<int>", "minimum MAPQ of index read for a pair to be mapped [5]"],
    }],
    bin =>[
        \&svtools_bin,
        "define variable-width bins for a chromosome based on depth across all samples", {
		'group'=>           [1, 1, @group, "Input Options"],
        'gap-file'=>        [0, 2, "g", "<str>", "BED file of gaps in chromosomes"],
        'exclude-file'=>    [0, 3, @exclude],        
        'output-dir'=>      [1, 4, @output, "Output Options"],        
        'tmp-dir'=>         [0, 5, @tmpDir],
		'max-mem'=>         [0, 6, @maxMem], 
        'chrom'=>           [1, 7, @chrom, "Call Options"],
		'chroms'=>          [1, 8, @chroms],
        'bin-count'=>       [0, 9, "n", "<int>", "number of fragments per variable-width bin [100]"],
		'ref-samples'=>     [0, 10, "r", "<str>", "comma-delimited list of samples used to define bins [all samples]"],
    }],
    call =>[
        \&svtools_call,
        "analyze a set of variable-width bins for copy number status based on depth", {
		'group'=>           [1, 1, @group, "Input Options"],
		'input-dir'=>       [1, 2, @input],
		'delete'=>          [0, 3, @delete],
        'output-dir'=>      [1, 4, @output, "Output Options"],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem], 
		'reference'=>       [1, 7, "r", "<str>", "copy number reference region as chrom:start-end/copy_number", "Call Options"],
		'max-copy-number'=> [0, 8, "c", "<int>", "maximum copy number allowed in the HMM [4]"],
		'chroms'=>          [1, 9, @chroms],
		'min-bins'=>        [0, 10, "b", "<int>", "minimum adjacent deviant bins to call a CNV [4]"],
    }],
    exclude =>[
        \&svtools_exclude,
        "create an exclusion file of regions with excessive copy number based on non-excluded call file", {
        'output-file'=>     [1, 2, @outFile, "Exclude Options"], 			
		#'exclude-file'=>    [0, 1, "x", "<str>", "merge this existing exclusions BED file into the output"],
        'copy-number'=>     [0, 3, "c", "<int>", "exclude spans with copy number >=c [5]"],
    }],
    depth_cdf =>[
        \&svtools_depth_cdf,
        "create a CDF plot of bin counts for all samples in a call file", {
		'group'=>           [1, 1, @group, "Input Options"],
		'input-dir'=>       [1, 2, @input],
        'tmp-dir'=>         [0, 3, @tmpDir],
        'output-dir'=>      [1, 4, @output, "Output Options"],
    }],
#------------------------------------------------------------------------
# '10X' approach, i.e. single-cell analysis in 10X scCNV platform
#------------------------------------------------------------------------
#   discover10X =>[
#        \&svtools_discover10X,
#        "identify true cells in a 10X scCNV BAM file, with their associated barcodes", {
#        'expected-cells'=>  [1, 1, "C", "<int>", "approximate number of true cells expected in data", "Discovery Options"],
#		  'min-qual'=>        [0, 2, @minMapQ10X],        
#        #'exclude-file'=>    [0, 3, @exclude],
#        'output-dir'=>      [1, 4, @output, "Output Options"],      
#	     'sample'=>          [1, 5, @sample],
#    }],
#    bin10X =>[
#        \&svtools_bin10X,
#        "create a bin coverage map for true cells in a 10X scCNV BAM file on a chromosome", {
#        'input-dir'=>       [1, 1, @input, "Input Options"],          
#		'sample'=>          [1, 2, @sample],
#        'chrom'=>           [1, 3, @chrom],        
#        'weight-per-cell'=> [0, 4, @binWeight, "Bin Options"],
#	     'min-qual'=>        [0, 5, @minMapQ10X],        
#        #'gap-file'=>        [0, 6, @gapFile],
#        #'exclude-file'=>    [0, 7, @exclude],
#        'output-dir'=>      [0, 8, @output, "Output Options"],
#        'tmp-dir'=>         [0, 9, @tmpDir],
#    }],
    bin10X =>[
        \&svtools_bin10X,
        "create a bin coverage map for true cells in a 10X scCNV BAM file on a chromosome", {
        'input-dir'=>       [1, 1, @input, "Input Options"],     
        'chrom'=>           [1, 2, @chrom],        
        'weight-per-cell'=> [0, 11, @binWeight, "Bin Options"],
	    'min-qual'=>        [0, 12, @minMapQ10X],        
        #'gap-file'=>        [0, 6, @gapFile],
        #'exclude-file'=>    [0, 7, @exclude],
        'output-dir'=>      [1, 21, @output, "Output Options"],
		'sample'=>          [1, 22, @sample],        
        'tmp-dir'=>         [0, 23, @tmpDir],
    }],
    merge10X =>[
        \&svtools_merge10X,
        "combine 10X coverage map over all chromosomes and make aggregate plots", {
        'input-dir'=>       [1, 1, @input, "Input Options"],          
		  'sample'=>          [1, 2, @sample],
        'exclude-file'=>    [0, 3, @exclude],
        'output-dir'=>      [0, 4, @output, "Output Options"],
    }],
    cluster10X =>[
        \&svtools_cluster10X,
        "calculate Z scores for bins x cells and perform hierarchical clustering", {
        'input-dir'=>       [1, 1, @input, "Input Options"],
        'cell-ranger-dir'=> [1, 2, "r", "<str>", "directory containing the CellRanger output files"],
		'sample'=>          [1, 3, @sample],
        'genome'=>          [1, 4, "g", "<str>", "name of the genome for this sample (used by server)"],
        'output-dir'=>      [0, 11, @output, "Output Options"],
        'include-Y'=>       [0, 12, "Y", undef, "include the Y chromosome in the output [FALSE if -y not used]"]
    }],
    
    
    createHMM10X =>[
        \&svtools_createHMM10X,
        "create a Hidden Markov Model for an intended bin weight", {
        'weight-per-cell'=> [0, 1, @binWeight, "HMM Options"],       
    }],
    segment10X => [
        \&svtools_segment10X,
        "find genome spans different than the modal copy number (i.e. CNVs), per cell", {
        'input-dir'=>       [1, 1, @input, "Input Options"],          
	     'sample'=>          [1, 2, @sample],
        'chrom'=>           [1, 3, @chrom],    
        'transition-prob'=> [0, 11, "t", "<dbl>", "HMM transition probability [1e-4]", "Segmentation Options"],
        #'preservation'=>    [0, 12, "P", "<str>", "penalities for CN changes from model [0.99,0.95,0.9,0.8]"],
        'weight-per-cell'=> [0, 13, @binWeight],
        #'bad-probe-freq'=>  [0, 23, "b", "<dbl>", "assume this fraction of probes give unpredictable values [1e-3]"],        
        'output-dir'=>      [0, 21, @output, "Output Options"],
        #'tmp-dir'=>         [0, 40, @tmpDir],
        #'max-mem'=>         [0, 43, @maxMem],        
    }],
#------------------------------------------------------------------------
# 'MLE' approach, i.e. find small CNVs by maximum likelihood estimate of TLENs
#------------------------------------------------------------------------
    scan => [
        \&svtools_scan,
        "calculate deviations of bin cumulative TLEN distributions vs. expected values", {
		'sample'=>          [1, 1, @sample, "Input Options"],
        'pdf-dir'=>         [1, 2, @pdfDir],
		'exclude-file'=>    [0, 3, @exclude],            
        'tmp-dir'=>         [0, 4, @tmpDir],
        'max-mem'=>         [0, 5, @maxMem],        
        'output-dir'=>      [1, 10, @output, "Output Options"],
        'min-qual'=>        [0, 20, @minMapQ, "Scan Options"],
        'genome-bin-size'=> [0, 21, "B", "<int>", "genome fixed-width bin size in bp [1000]"],        
        'max-cnv-size'=>    [0, 24, "Z", "<int>", "largest CNV size in bp that will be interrogated at maximal sensitivity [10000]"],
        'min-scan-pairs'=>  [0, 30, "e", "<int>", "require at least this many pairs to perform scan [10]"],
        'max-scan-pairs'=>  [0, 31, "E", "<int>", "randomly downsample to this many pairs during scan [100]"],
		'min-neg-pairs'=>   [0, 32, "N", "<int>", "ignore pairs with TLEN<0 unless >=N such pairs are present in bin [3]"],
        'chroms'=>          [0, 43, @chroms],
        'region'=>          [0, 44, "R", "<str>", "only look at this region (chr:-start-end); print bins to STDOUT"]
    }],
    crunch => [
        \&svtools_crunch,
        "find runs of scanned bins that are deviant, i.e. are likely to be small CNVs", {
		'sample'=>          [1, 1, @sample, "Input Options"],
		'input-dir'=>       [1, 2, @input],
        'pdf-dir'=>         [1, 2, @pdfDir],
        'tmp-dir'=>         [0, 4, @tmpDir],
        'max-mem'=>         [0, 5, @maxMem], 		
        'output-dir'=>      [1, 10, @output, "Output Options"],		
        'min-cnv-delta=>'=> [0, 40, "d", "<dbl>", "minimum deltaCDF for bin to be considered deviant [5]", "Crunch Options"],
        'min-cnv-bins=>'=>  [0, 41, "n", "<int>", "minimum adjacent deviant bins for CNV to be called [5]"],
		#'min-cnv-pairs=>'=> [0, 42, @minPair],
        'cnv-size-step'=>   [0, 52, "v", "<int>", "query CNV sizes separated by this many bp [500]"],
        'min-cnv-size'=>    [0, 53, "z", "<int>", "smallest CNV size in bp that will be considered non-reference [2000]"],
        'max-cnv-size'=>    [0, 54, "Z", "<int>", "largest CNV size in bp that will be interrogated at maximal sensitivity [10000]"],
        'region'=>          [0, 63, "R", "<str>", "only look at this region (chr:-start-end)"]
    }],
    collapse => [
        \&svtools_collapse,
        "find overlapping/unique CNVs in samples previously analyzed with scan+crunch", {
		'group'=>           [1, 1, @group, "Input Options"],
        'tmp-dir'=>         [0, 2, @tmpDir],
        'max-mem'=>         [0, 3, @maxMem],  
        'output-dir'=>      [1, 4, @output, "Output Options"],
        'min-cnv-delta=>'=> [0, 40, "d", "<dbl>", "minimum deltaCDF for cnv to be called [5]", "Collapse Options"],
    }],    
    
#    mle => [
#        \&svtools_mle,
#        "maxiumum likelihood estimation of small del/dup from coordinate-sorted bam files", {
#		'sample'=>          [1, 1, @sample, "Input Options"],
#        'pdf-dir'=>         [1, 2, @pdfDir],
#		'exclude-file'=>    [0, 3, @exclude],            
#        'tmp-dir'=>         [0, 4, @tmpDir],
#        'max-mem'=>         [0, 5, @maxMem],        
#        'output-dir'=>      [1, 10, @output, "Output Options"],
#        'min-qual'=>        [0, 20, @minMapQ, "MLE Options"],
#        'genome-bin-size'=> [0, 21, "B", "<int>", "genome fixed-width bin size in bp [1000]"],        
#        'cnv-size-step'=>   [0, 22, "v", "<int>", "query CNV sizes separated by this many bp [500]"],
#        'min-cnv-size'=>    [0, 23, "z", "<int>", "smallest CNV size in bp that will be considered non-reference [2000]"],
#        'max-cnv-size'=>    [0, 24, "Z", "<int>", "largest CNV size in bp that will be interrogated at maximal sensitivity [10000]"],
#        'min-mle-pairs'=>   [0, 30, "e", "<int>", "require at least this many pairs to perform MLE [10]"],
#        'max-mle-pairs'=>   [0, 31, "E", "<int>", "randomly downsample to this many pairs before MLE [100]"],
#        'min-cnv-bins=>'=>  [0, 41, "n", "<int>", "minimum adjacent bins for CNV to be called [3]"],
#		'min-cnv-pairs=>'=> [0, 42, @minPair],
#        'region'=>          [0, 43, "R", "<str>", "only look at this region (chr:-start-end); print bins to STDOUT"]
#    }],
#------------------------------------------------------------------------
# 'SNP' approach, i.e. find runs of unidirectional LOH with respect to known SNPs
#------------------------------------------------------------------------	
    genotype => [
        \&svtools_genotype,
        "determine read counts for a sample at a list of predefined SNPs", {
		'sample'=>          [1, 1, @sample, "Input Options"],
		'genome-fasta'=>    [1, 2, "g", "<str>", "path to the genome FASTA file, e.g. /path/hg19.fa"],
		'snp-file'=>        [1, 3, "S", "<str>", "path to the SNP BED file, format chr,pos-1,pos,name,x,x,Ref,Alt"],
		'exclude-file'=>    [0, 4, @exclude],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],  		
        'output-dir'=>      [1, 10, @output, "Output Options"],
        'min-qual'=>        [0, 20, @minMapQ, "Genotype Options"]
    }],
    pool => [
        \&svtools_pool,
        "aggregate the genotypes from a series of samples to determine heterzygous SNPs", {
		'group'=>           [1, 1, @group, "Input Options"],
        'input-dir'=>       [1, 2, @input],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],  		
        'output-dir'=>      [1, 10, @output, "Output Options"],
    }],
    filter => [
        \&svtools_filter,
        "filter heterozygous snps from svtools pool output", {
		'group'=>           [1, 1, @group, "Input Options"],
        'input-dir'=>       [1, 2, @input],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],  		
        'output-dir'=>      [1, 10, @output, "Output Options"],
        'min-reads'=>       [0, 20, "r", "<int>", "only consider SNPs with >=r crossing reads [10]"],
        'max-reads'=>       [0, 21, "R", "<int>", "only consider SNPs with <=R crossing reads [100]"],
        'min-frac-alt'=>    [0, 22, "f", "<dbl>", "only consider SNPs with >=f non-ref reads [0.3]"],
        'max-frac-alt'=>    [0, 23, "F", "<dbl>", "only consider SNPs with <=F non-ref reads [0.7]"], 
    }],
    bias => [
        \&svtools_bias,
        "search for runs of SNP allelic bias indicative of CNVs on one allele", {
		'sample'=>          [1, 1, @sample, "Input Options"],	
        'input-dir'=>       [1, 2, @input],
        'snp-file'=>        [1, 3, "S", "<str>", "path to the SNP BED file, format chr,pos-1,pos,name,x,x,Ref,Alt,nRef,nAlt"],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],		
        'output-dir'=>      [1, 10, @output, "Output Options"],
		'chroms'=>          [1, 20, @chroms, "Bias Options"],
		'zero-prob'=>       [0, 30, "z", "<dbl>", "probability of being in state 0, i.e. heterozygous [0.99]"],
		'persistence'=>     [0, 31, "p", "<dbl>", "probability of remaining in state [0.995]"],
		'base-error-rate'=> [0, 32, "e", "<dbl>", "probability that a base call is erroneous [0.001]"],
    }],
    collapse_bias => [
        \&svtools_collapse_bias,
        "find overlapping/unique CNVs in samples previously analyzed with genotype+bias", {
		'group'=>           [1, 1, @group, "Input Options"],
        'tmp-dir'=>         [0, 2, @tmpDir],
        'max-mem'=>         [0, 3, @maxMem],  
        'output-dir'=>      [1, 4, @output, "Output Options"],
        'min-pRatio=>'=>    [0, 40, "R", "<dbl>", "minimum pRatio for cnv to be called [10]", "Collapse Options"],
    }],
	
#------------------------------------------------------------------------
# filter and sort output from pairs, depth and MLE to find candidate CNVs
#------------------------------------------------------------------------	
    view =>[
        \&svtools_view,
        "list structural variant junctions/pairs from find/index", {
		'command'=>         [0, 1,  "C", "<str>", "svtools command that generated the data (find call collapse collapse_bias)", "Input Options"], 	
		'data-type'=>       [0, 2,  "d", "<str>", "input data type (sets cnvs regions)"], 		
		#------------------
        'output-file'=>     [0, 3,  @outFile, "Output Options"],
        'human'=>           [0, 4,  "h", undef,  "generate a human readable report [svtools format]"],
        'count'=>           [0, 5,  "c", undef,  "show the number of matching junctions/pairs"],
		#------------------
        'min-size'=>        [0, 10,  "z", "<int>", "only report variants >= z bp on the same chromosome", "Common Filter Options"],
        'max-size'=>        [0, 11,  "Z", "<int>", "only report variants <= Z bp on the same chromosome"],
        'min-call'=>        [0, 12,  "l", "<int>", "only report variants with >=l called samples"],	
        'max-call'=>        [0, 13,  "L", "<int>", "only report variants with <=L called samples"],
        'min-pos'=>         [0, 14,  "f", "<int>", "only report variants with >=f positive samples"],	
        'max-pos'=>         [0, 15,  "F", "<int>", "only report variants with <=F positive samples"],			
        'sample'=>          [0, 20,  "s", "<str>", "only report variants where sample s was called"],
		'min-bins'=>        [0, 21,  "b", "<int>", "only report cnvs with >= b bins"],
		'min-set-pairs'=>   [0, 22,  "p", "<int>", "only report variants with >=p crossing fragments"],
        'max-set-pairs'=>   [0, 23,  "P", "<int>", "only report variants with <=P crossing fragments"],
		'min-PRatio'=>      [0, 24,  "R", "<dbl>", "only report cnvs where sample PRatio >= R"],		
		#------------------           
        'jxn-types' =>      [0, 30,  "t", "<str>", "junction/CNV type(s) (del,dup,invF,invR,trans)", "find.sets Options"],        
        'jxn-only'=>        [0, 33, @junction],
        'condensed'=>       [0, 35,  "x", undef,   "output one line per same-chrom CNV, e.g. for plot_prep"],
        #------------------
		'cn-types' =>       [0, 50, "T", "<str>", "copy number type(s) (gain,loss,group)", "call.cnvs Options"],
		'max-uniq'=>        [0, 52, "q", "<int>", "only report cnvs with <= q uniqueness score"],
		'min-delta-CN'=>    [0, 53, "r", "<dbl>", "only report cnvs where sample abs(mean(CN delta)) >= r"],
        #------------------
		'bias-types' =>     [0, 60, "U", "<str>", "bias type(s) (refDup,altDup,refDel,altDel)", "collapse_bias.cnvs Options"],
		'min-snps' =>       [0, 61, "w", "<int>", "minimum number of SNPs in index CNV"],
		'min-reads' =>      [0, 62, "W", "<int>", "minimum number of reads crossing SNPs in index CNV"],
		'min-bias' =>       [0, 63, "v", "<dbl>", "minimum abs(bias) in index CNV"],
		'max-other-bias' => [0, 64, "V", "<dbl>", "maximum abs(bias) in other samples in CNV span"],
        #------------------
		'scan-types' =>     [0, 70, "u", "<str>", "copy number type(s) (del,dup,ins)", "collapse.regions Options"],		
		'min-dlt-sample' => [0, 72, "e", "<dbl>", "minimum deltaCDF for index sample"],
		'max-dlt-other' =>  [0, 73, "E", "<dbl>", "maximum deltaCDF for any other sample"],
		'suppress-gapped'=> [0, 74, "y", undef,   "skip regions where the index sample had a gap"],
		'max-other-gap' =>  [0, 75, "Y", "<int>", "maximum number of other samples allowed to have a gap"],
		'max-low-qual' =>   [0, 76, "Q", "<dbl>", "maximum fraction of low-quality fragments in index sample"],
    }],
    sort =>[ # not implemented yet!
        \&svtools_sort,
        "sort structural variants based on column name", {
		'command'=>         [1, 1, "C", "<str>", "svtools command that generated the data (find call collapse)", "Sort Options"], 	
		'data-type'=>       [1, 2, "d", "<str>", "input data type (sets cnvs regions)"], 	
        'sort-col'=>        [1, 3, "c", "<str>", "column name to sort by"],
        'reverse'=>         [0, 4, "r", undef,   "reverse order, e.g. descending"],
        'tmp-dir'=>         [0, 5, @tmpDir],
        'max-mem'=>         [0, 6, @maxMem],
    }],
#------------------------------------------------------------------------
# make plots of candidate regions for easy screening
#------------------------------------------------------------------------	
    plot_prep => [
        \&svtools_plot_prep,
        "prepare an input stream on STDIN for passing to a plot_XXX command via STDOUT", {
		'group'=>           [1, 1, @group, "Prep Options"],																							  
		'input-type'=>      [0, 2, "i", "<str>", "the type of svtools features on STDIN (pairs depth scan snp)"],
        'id-col'=>  		[0, 3, "f", "<int>", "data column that contains the feature id"],		
        'type-col'=>  		[0, 4, "t", "<int>", "data column that contains the feature type"],
        'size-col'=>  		[0, 5, "z", "<int>", "data column that contains the feature size"],
        'frac-col'=>  		[0, 6, "F", "<int>", "data column that contains the CNV fraction"],
        'sample-col'=>      [0, 7, "s", "<int>", "data column that contains the index sample id"],
    }], 
    plot_scan => [
        \&svtools_plot_scan,
        "make plots of p-values for all spans in a prepared 'collapse' stream on STDIN", {			
		'group'=>           [1, 1, @group, "Plot Options"],
        'plot-dir'=>        [1, 2, @plotDir],
        'pdf-dir'=>         [1, 3, @pdfDir],
		'tmp-dir'=>         [0, 4, @tmpDir],
		'center-width'=>    [0, 10, "c", "<int>", "place BED span in central 1/c of plot width [5]"],
        'bam-glob'=>        [1, 20, "b", "<int>", "how to find bam files, e.g. <path>/SAMPLE.*.bam"],
        'max-cnv-size'=>    [0, 24, "Z", "<int>", "largest CNV size in bp, used to set TLEN limits [10000]"],
		'min-qual'=>        [0, 25, @minMapQ],
    }],
    plot_cnv => [
        \&svtools_plot_cnv,
        "make plots of coverage depth and anomalous pairs in a prepared stream on STDIN", {
		'call-bins-file'=>  [1, 1, "b", "<str>", "call.bins file used to make depth plots", "Input  Options"],
        'idx-anom-file'=>   [1, 2, "a", "<str>", "index.anomalies file used to make fragment plots"],
        'gen-snps-dir'=>    [1, 3, "S", "<str>", "directory with genotype.snps filed used to make SNP plots"],
		'group'=>           [1, 10, @group, "Plot Options"],
        'plot-dir'=>        [1, 11, @plotDir],
		'tmp-dir'=>         [0, 20, @tmpDir],
		'center-width'=>    [0, 30, "c", "<int>", "place BED span in central 1/c of plot width [5]"],
        'gap-file'=>        [0, 40, "g", "<str>", "BED file of gaps in chromosomes"],
        'seg-dup-file'=>    [0, 50, "d", "<str>", "BED file of segmental duplications in chromosomes"],
        'max-snp-count'=>   [0, 60, "s", "<int>", "y limit on SNP count plot [5]"],
		'ref-name'=>        [0, 62, "r", "<int>", "short name for reference SNP alleles [Ref]"],
		'alt-name'=>        [0, 64, "R", "<int>", "short name for alternative SNP alleles [Alt]"],
    }],
#------------------------------------------------------------------------
# pull pairs that flank presumed CNV junctions
#------------------------------------------------------------------------	
	# pending
);
our %longOptions; # for converting short options to long
foreach my $command(keys %commands){
    %{$longOptions{$command}} =
        map { $commands{$command}[2]{$_}[2] => $_ } keys %{$commands{$command}[2]}; 
}
#------------------------------------------------------------------------
# break utility commands into managable help chunks 
#------------------------------------------------------------------------
sub reportCommandChunks {
    reportCommandChunk('Preliminary Analysis', qw(pdf));
    reportCommandChunk('Pairs Approach', qw(extract mask find cat index pull));	
    reportCommandChunk('Depth Approach', qw(map bin call exclude depth_cdf));
    reportCommandChunk('TLEN Approach', qw(scan crunch collapse));
    #reportCommandChunk('MLE Approach', qw(mle collapse));
	 reportCommandChunk('SNP Approach', qw(genotype pool filter bias));
    reportCommandChunk('10X scCNV Depth', qw(bin10X merge10X cluster10X createHMM10X segment10X)); #discover10X
    reportCommandChunk('Filter and Sort', qw(view sort));
	reportCommandChunk('Plot Candidates', qw(plot_prep plot_cnv));
	#reportCommandChunk('Read Pair Recovery', qw(xxx));
}
#========================================================================

1;
