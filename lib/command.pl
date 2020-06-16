use strict;
use warnings;

#========================================================================
# 'command.pl' has subs for parsing and executing the command line
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utilityDescription $utility $version $command @args @optionGroups %useOptionGroupDelimiter  
            %commands %options %longOptions $libDir);
my (@options, %optInfo);
my $commandTabLength = 12; 
my $optionTabLength = 25;
my $separatorLength = 87;
our $error = "error";
#========================================================================

#========================================================================
# parse execution command and options
#-----------------------------------------------------------------------
sub main {
    $error = "$utility error";
    checkCommand();  # provide syntax checking to prevent malformed commands
    $error = "$utility $command error";
	@args = setOptions();
    checkRequiredOptions();
    executeCommand();  # request is valid, proceed with execution
}
sub checkCommand { # check for help request or validity of requested command
    my $versionString = "$utility version $version\n";
    my $descriptionString = "$versionString$utilityDescription";
    $command or reportUsage($descriptionString);	
    ($command eq '-v' or $command eq '--version') and print $versionString and exit;     
    ($command eq '-h' or $command eq '--help') and reportUsage("$descriptionString\n", "all");  
    $commands{$command} or reportUsage("\n$error: '$command' is not a valid command\n", "all", 1);
    %optInfo = %{$commands{$command}[2]};
}
sub executeCommand { # load scripts and execute command
    #foreach my $scriptName(keys %commands){
    #    my $script = "$libDir/$scriptName.pl";     
    #    require $script;
    #}
    require "$libDir/$command.pl";
    &{${$commands{$command}}[0]}(@args);  # add remaining @args
}
#-----------------------------------------------------------------------
sub setOptions { # parse and check validity of options string
    while (my $optionList = shift @args){
        ($optionList and $optionList =~ m/^\-./) or return ($optionList, @args); # last item(s) in list should be masterFile(s)
        push @options, $optionList;    
        if($optionList =~ m/^\-\-(.+)/){ # long option formatted request
            my $longOption = $1;
            defined $optInfo{$longOption} or reportUsage("\n$error: '$longOption' is not a recognized option\n", $command, 1); 
            setOption($longOption);
        } elsif ($optionList =~ m/^\-(.+)/){ # short option formatted request
            foreach my $shortOption(split('', $1)){
                my $longOption = $longOptions{$command}{$shortOption};
                defined $longOption or reportUsage("\n$error: '$shortOption' is not a recognized option\n", $command, 1);  
                setOption($longOption);
            }
        } else {
            reportUsage("$error: malformed option list", undef, 1); 
        }
    }
    return undef;
}           
sub setOption { # check and set option request                
    my ($longOption) = @_;
    $longOption eq 'version' and print "$utility version $version\n" and exit;   
    $longOption eq 'help' and reportUsage("$utility version $version\n", $command); 
    defined $optInfo{$longOption} or reportUsage("\n$error: '$longOption' is not a valid option\n", $command, 1);
    my $value; # boolean options set to value 1, otherwise use supplied value
    if(${$optInfo{$longOption}}[3]){
        my @values;
        #while($args[0] and !($args[0] =~ m/^\-/)){ push @values, shift @args }
        if($args[0] and !($args[0] =~ m/^\-/)){ push @values, shift @args }
        $value = join(" ", @values);
        $value eq "" and $value = undef;
        push @options, $value;    
    } else {
        $value = 1;
    }
    defined $value or reportUsage("\n$error: missing value for option '$longOption'\n", $command, 1);
    $value =~ m/^\-/ and reportUsage("\n$error: missing value for option '$longOption'\n", $command, 1);    
    $options{$longOption} = $value;  
}
sub checkRequiredOptions { # make sure required value options have been supplied
    foreach my $longOption (keys %optInfo){
        $optInfo{$longOption}[0] or next; # option is not required
        defined $options{$longOption} or reportUsage("\n$error: option '$longOption' is required\n", $command, 1);
    }
}
#========================================================================

#========================================================================
# provide help on command line
#------------------------------------------------------------------------
sub reportUsage { 
    my ($message, $command, $die) = @_;
    $message and print "$message\n";  
    print "usage:\t$utility <command> [options] <input file(s)> [chr:start-end]\n";
    if($command){
        if($commands{$command}){ # help for options for known command
            reportOptionHelp($command);
        } else { # help on the set of available commands, organized by topic
            print "\navailable commands (use '$utility <command> --help' for command options):\n\n";
            reportCommandChunks();     
        }
    } else {
        print "use '$utility --help' or '$utility <command> --help' for extended help\n", 
    }   
    my $exitStatus = $die ? 1 : 0;
    exit $exitStatus; 
}
sub reportOptionHelp { 
    my ($command) = @_;
    print "\n";
    print getCommandLine($command);
    my @longOptions = sort {$optInfo{$a}[1] <=> $optInfo{$b}[1]} keys %optInfo;
    if(@longOptions){
        my %parsedOptions;
        foreach my $longOption(@longOptions){
            my ($required, $sort, $shortOption, $type, $help, $header) = @{$optInfo{$longOption}};        
            my $option = "-$shortOption,--$longOption";
            $type and $option .= " $type";
            $required and $help = "**REQUIRED** $help";
            $header and print "\n$header\n";
            print "    $option".(" " x ($optionTabLength - length($option)))."$help\n";
        }
    } else {
        print "    none\n";
    }
    print "\n";
}
sub reportCommandChunk {
    my ($header, @commands) = @_;
    $header and print "$header\n";
    foreach my $command (@commands){
        print "    ", getCommandLine($command);
    }
    print "\n";
}
sub getCommandLine {
    my ($command) = @_;
    return $command, " " x ($commandTabLength - length($command)), ${$commands{$command}}[1], "\n";
}
#========================================================================

#========================================================================
# handle wall time
#------------------------------------------------------------------------
sub getTime { # status update times
    my ($sec, $min, $hr, $day, $month, $year) = localtime(time);
    $year = $year + 1900;
    $month++;
    return "$month/$day/$year $hr:$min:$sec";
}
sub printTime {
    my ($message) = @_;
    $message or $message = "";
    $message and $message .= ": ";
    my $time = getTime();
    print STDERR "$message$time";
}
#========================================================================

#========================================================================
# set options from within command subs
#------------------------------------------------------------------------
sub setOptionValue {
    my ($var, $option, $default) = @_;
    $$var = defined $options{$option} ? $options{$option} : $default;
}
#========================================================================

1;