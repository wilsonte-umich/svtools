use strict;
use warnings;

#========================================================================
# 'file.pl' has subs for handling IO
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utility $command $slurp $error
            $tmpDir $outputDir $inputDir);
#========================================================================

#========================================================================
# handle multiple file inputs
#------------------------------------------------------------------------
sub processInputFiles {
    my ($inH, $globs, $allowStdin, $stepsAfter, $isGZip, $sub) = @_;
    foreach my $glob(@$globs){
        if ($glob eq '-' or lc($glob) eq 'stdin') {
            $allowStdin or die "$utility $error: $command does not support stdin as input\n";
            openInputStream(undef, $inH, undef, $stepsAfter, $isGZip);
            &$sub('stdin');
            return;
        } else {
            foreach my $file(glob($glob)){
                -d $file and die "$utility $error: $file is a directory\n";
                openInputStream($file, $inH, undef, $stepsAfter, $isGZip);
                &$sub($file);
            }
        }
    }
}
#========================================================================

#========================================================================
# handle I/O streams
#------------------------------------------------------------------------
sub openInputStream  {
    my ($file, $h, $slurpSize, $stepsAfter, $isGZip) = @_;
    if(!$file or ($file eq "-" or uc($file) eq "STDIN")){
        if($stepsAfter){
            open $$h, "-|", $stepsAfter or die "$error: could not open input stream for reading: $!\n"; 
        } else {
            $$h = *STDIN;        
        }	
    } else {
        $stepsAfter or $stepsAfter = "";	
        $stepsAfter and $stepsAfter = " | $stepsAfter";    
        defined $slurpSize or $slurpSize = "50M";
		$isGZip = $isGZip || $file =~ m|\.gz$|;
        my $gunzip = $isGZip ? "| gunzip -c" : "";
		my $stream = "$slurp -s $slurpSize $file $gunzip $stepsAfter";
        open $$h, "-|", $stream or die "$error: could not open $file for reading: $!\n";
    }    
}
sub openBamStream {
    my ($file, $h, $slurpSize, $minQ, $f, $F, $region, $samtools, $stepsAfter) = @_;
    $minQ or $minQ = "5";
    $f = $f ? "-f $f" : "";
    $F = $F ? "-F $F" : "";
    $region or $region = "";
    $samtools or $samtools = "samtools";
    $stepsAfter or $stepsAfter = "";
    $stepsAfter and $stepsAfter = " | $stepsAfter";    
    defined $slurpSize or $slurpSize = "50M";
    
    #$stepsAfter = " | head -n 100000";
    
    my $stream = "$slurp -s $slurpSize $samtools view -q $minQ $f $F $file $region $stepsAfter";
    open $$h, "-|", $stream or die "$error: could not open $file for reading: $!\n"; 
}
sub openBamHeader {
    my ($file, $h, $samtools, $stepsAfter) = @_;
    $samtools or $samtools = "samtools";
    $stepsAfter or $stepsAfter = "";
    $stepsAfter and $stepsAfter = " | $stepsAfter";    
    my $stream = "$samtools view -H $file $stepsAfter";
    open $$h, "-|", $stream or die "$error: could not open $file for reading: $!\n";    
}
sub openOutputStream {
    my ($file, $h, $gzip, $slurpSize, $stepsBefore) = @_;
    if(!$file or uc($file) eq "STDOUT"){
        if($stepsBefore){
            open $$h, "|-", $stepsBefore or die "$error: could not open output stream for writing: $!\n"; 
        } else {
            $$h = *STDOUT;        
        }
    } else {         
        $stepsBefore or $stepsBefore = "";
        $stepsBefore and $stepsBefore = "$stepsBefore | ";    
        defined $slurpSize or $slurpSize = "50M";
        $gzip = $gzip ? "gzip -c |" : "";
        $file = ($gzip and !($file =~ m|\.gz$|)) ? "$file.gz" : $file;
		my $stream = "$stepsBefore $gzip $slurp -s $slurpSize -o $file";
        open $$h, "|-", $stream or die "$error: could not open $file for writing: $!\n";
		return $stream;
    }
}
sub closeHandles {
    foreach my $h(@_){
        $h and close $h;
    }
}
#========================================================================

#========================================================================
# create consistent file names
#------------------------------------------------------------------------
sub getTmpFile {
	my ($type, $name, $suffix) = @_;
	my $rand = int(rand(1e6));
	return "$tmpDir/$utility.$command.$type.$name.$rand.$suffix";
}
sub getOutFile {
	my ($type, $name, $suffix) = @_;
	return "$outputDir/$utility.$command.$type.$name.$suffix";
}
sub getInFile {
	my ($command, $type, $name, $suffix) = @_;
	return "$inputDir/$utility.$command.$type.$name.$suffix";  
}
sub getFileType {
	my ($file) = @_;
	$file =~ m/$utility\.(\w+\.\w+)/ or die "not recognized as $utility file: $file\n";
	return $1;
}
#========================================================================

#========================================================================
# recover a sample name list
#------------------------------------------------------------------------
sub getSamples {
    my ($smpFile) = @_;
    openInputStream($smpFile, \my $inH);
    my $samples = <$inH>;
    chomp $samples;
    my @samples = split("\t", $samples);
    close $inH;
    return @samples;
}
#========================================================================

1;
