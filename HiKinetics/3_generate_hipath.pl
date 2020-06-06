#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use POSIX;
use threads;
use threads::shared;
 
my $outputrootname = $ARGV[0];


# my $data_file="files_to_be_processed.txt";
# open(DAT, $data_file) || die("Could not open file!");
# my @raw_data=<DAT>; 
# close(DAT);

#print "Starting main program\n";
my @childs;
# my $wrestler;
# foreach $wrestler (@raw_data)
# {
#     my $outputrootname;
#     my $extension;
#     chomp($wrestler);
#     ($outputrootname,$extension)=split(/\./,$wrestler);


    my $filename = $outputrootname.'.res';
    my $lines = -2; #0;
    my $buffer;
    open(FILE, $filename) or die "Can't open `$filename': $!";
    while (sysread FILE, $buffer, 4096) {
	$lines += ($buffer =~ tr/\n//);
    }
    close FILE;
#LIMITED with 20 
#     if($lines>20){
#       $lines=20;
#     }
    #print "The number of lines in $filename is $lines.\n";

    my $num = 1;
    my $i=0;
    for ($i=0;$i<$lines;$i++)
    {


		    my $pid = fork();
		    if($pid){
			    # parent
			# print "pid is $pid, parent $$\n";
			    push(@childs, $pid);
		    }elsif ($pid == 0){
			    # child
			    #rnaZ($file_path, "RNAzOutput" , $f, $num);
			  # my $inFile = "./".$outputrootname.".faas/".$i."_".$j.".faa";
			  # my $outFile = "./".$outputrootname.".faas/".$i."_".$j.".hip"; 
			  # system "HiPath -f $inFile -t 1 > $outFile";
			    my $j=0;
			    for ($j=0;$j<$lines;$j++)
			    {
                                #print "HiPath -f ${outputrootname}.faas/${i}_${j}.faa -t 1 > ${outputrootname}.faas/${i}_${j}.hip\n";
			        unless (-e "${outputrootname}.faas/${i}_${j}.hip") { 
                                    #TODO: current HiPath only works with turner1999
                                    system "HiPath -f ${outputrootname}.faas/${i}_${j}.faa -t 1 > ${outputrootname}.faas/${i}_${j}.hip";
                                }
			    }
			    exit 0;
		    }else{
			    die "couldnt fork: $!\n";
		    }
		    $num++;

    }

#}
 
#print "waitpid\n";
foreach (@childs) {
        my $tmp = waitpid($_, 0);
        #print "done with pid $tmp\n";
}
 
#print "End.hip calculation\n";
 
