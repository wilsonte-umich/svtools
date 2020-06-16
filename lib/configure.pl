#!/usr/bin/perl
use strict;
use warnings;
use Cwd(qw(abs_path));

#========================================================================
# 'configure.pl' checks prerequisites and creates the program target
#========================================================================

#========================================================================
use vars qw($utilityDescription @prerequisites);
$| = 1;
#------------------------------------------------------------------------
# get the current version
#------------------------------------------------------------------------
my $script = abs_path($0);
$script =~ m|(.*)/(.+)-(\d+\.\d+\.\d+)/configure.pl$| or die "fatal error: could not recover version information\n";
my ($wilsonDir, $utility, $version) = ($1, $2, $3);
my $scriptDir = "$wilsonDir/$utility-$version";
print "configuring $utility, version $version\n";
#------------------------------------------------------------------------
# check the license
#------------------------------------------------------------------------
my $getAgain = "please reacquire $utility source code from SourceForge";
my $licenseFile = "$scriptDir/LICENSE";
-e $licenseFile or die "fatal error: missing $utility software LICENSE\n$getAgain\n";
#------------------------------------------------------------------------
# parse the path to the program and lib
#------------------------------------------------------------------------
print "checking prerequisites ";
print ".";
$script = "$scriptDir/$utility";
my $libDir = "$scriptDir/lib";
-d $libDir or die "\nfatal error: missing 'lib' directory\n$getAgain\n";
#------------------------------------------------------------------------
# get path to perl
#------------------------------------------------------------------------
print ".";
my $perlPath = getProgramPath("perl");
chomp $perlPath;
#------------------------------------------------------------------------
# check for prequisites
#------------------------------------------------------------------------
foreach my $command(@prerequisites){ 
    print ".";
    getProgramPath($command);
}
sub getProgramPath {
    my ($command) = @_;
    my $path = whereis($command);
    $path or $path = which($command);
    $path or die "\nfatal error: missing prerequisite: $command\n";
    return $path;
}
sub whereis {
    my ($command) = @_;
    my $path = qx|whereis $command 2>/dev/null|;
    $path or return;
    chomp $path;
    $path or return;
    my @list = split(/\s+/, $path);
    return $list[1];
}
sub which {
    my ($command) = @_;
    my $path = qx|which $command 2>/dev/null|;
    $path or return;
    chomp $path;
    return $path;
}
#------------------------------------------------------------------------
# print the program target script
#------------------------------------------------------------------------
print "\ngenerating $utility program target\n";
open my $outH, ">", $script or die "could not open $script for writing: $!\n";
print $outH
'#!'.$perlPath.'
use strict;
use warnings;
our $utilityDescription = '."'$utilityDescription'".';
our $utility = '."'$utility'".';
our $version = '."'$version'".';
our $perlPath = '."'$perlPath'".';
our $scriptDir = '."'$scriptDir'".';
our $libDir = '."'$libDir'".';
our $slurp = '."'$libDir/slurp'".';
our $glurp = '."'$libDir/glurp'".';
$| = 1;
require "$libDir/main.pl";
require "$libDir/common.pl";
main();
';
close $outH;
#------------------------------------------------------------------------
# make the program target script executable
#------------------------------------------------------------------------
qx|chmod ugo+x $script|;  
print "done\n";
print "created $utility program target:\n  $script\n";   
#========================================================================

