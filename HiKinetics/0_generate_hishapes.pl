#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use POSIX;
use threads;
use threads::shared;
 
my $inputfilename = $ARGV[0];
my $outputrootname = $ARGV[1];
my $kbest = $ARGV[2];
my $tvalue = $ARGV[3];
my $x_thresh = $ARGV[4];


# my $data_file="files_to_be_processed.txt";
# open(DAT, $data_file) || die("Could not open file!");
# my @raw_data=<DAT>; 
# close(DAT);



# print "Starting main program\n";
my @childs;
# my $wrestler;
# foreach $wrestler (@raw_data)
# {
#     my $outputrootname;
#     my $extension;
#     chomp($wrestler);
#     ($outputrootname,$extension)=split(/\./,$wrestler);

    my $num = 1;

    my $pid = fork();
    if($pid){
	    # parent
	    # print "pid is $pid, parent $$\n";
	    push(@childs, $pid);
    }elsif ($pid == 0){
	    # child            -P ../src/librna/vienna/rna_turner2004.par 
	    system "RNAHeliCes -f $inputfilename -k $kbest -t $tvalue -x$x_thresh > $outputrootname.pre";
            unless (-e "$outputrootname.faas") {
	      system "mkdir $outputrootname.faas";
	    }   
	    exit 0;
    }else{
	    die "couldnt fork: $!\n";
    }
    $num++;

# }
 
#print "waitpid\n";
foreach (@childs) {
        my $tmp = waitpid($_, 0);
        #print "done with pid $tmp\n";
}
# print "End of main program\n";
 
