use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $tmpDir
            %findCols %cnvCols %clpCols);
my ($command, $dataType, $sortCol, $reverse, $maxMem);
my ($inH, $outH);
my %cols;

# manage options
sub setOptions_sort {
    setOptionValue(\$command,    'command');
    setOptionValue(\$dataType,   'data-type');
    setOptionValue(\$sortCol,    'sort-col');
    setOptionValue(\$reverse,    'reverse');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000); 
    #-----------------------------------------   
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";
    #-----------------------------------------   
    if ($command eq 'find' and $dataType eq 'sets') {
        %cols = %findCols;
    } elsif ($command eq 'call' and $dataType eq 'cnvs'){
        %cols = %cnvCols;
    } elsif ($command eq 'collapse' and $dataType eq 'regions'){
        %cols = %clpCols;
    } else {
        die "$error: unrecognized data-type: $command $dataType\n";
    } 
}

# main execution block
sub svtools_sort {
    setOptions_sort();

    # open IO streams
    my $col = $cols{$sortCol} + 1;
    my $r = $reverse ? "r" : "";
    my $sort = "sort -T $tmpDir -S $maxMem"."b -k$col,$col"."n$r";
    openInputStream(undef, \$inH, undef, $sort);
    openOutputStream(undef, \$outH, undef, undef);

    # display the results
    while (my $line = <$inH>) { print $outH $line }
}

1;
