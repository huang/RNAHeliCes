#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use POSIX;
use threads;
use threads::shared;
#use Data::Dump qw(dump);
#use vars qw(@result);
 
my $inputfilename = $ARGV[0];
my $outputrootname = $ARGV[1];
my $initial_hishape = $ARGV[2];


# my $data_file="files_to_be_processed.txt";
# open(DAT, $data_file) || die("Could not open file!");
# my @raw_data=<DAT>; 
# close(DAT);


# print "Starting main program\n";
# my $wrestler;
# foreach $wrestler (@raw_data)
# {
#     my $outputrootname;
#     my $extension;
#     chomp($wrestler);
#     ($outputrootname,$extension)=split(/\./,$wrestler);

    my $filename = $outputrootname.'.pre';
    my $line_no = -2; #0;
    my $line;

    my @childs;
    my @lines;
    #my @hishapes;
    my $line_minus2;
    my $line_minus1;

    #use vars qw($results);

    unless( open(FILE, $filename) ) {
        die "Error: couldn't open file '$filename': $!\n";
    } # if
    while ($line =<FILE>) {
        #print $line;
        chomp($line);
        if ($line_no>=0) {
            $lines[$line_no] = $line;
        } elsif ($line_no==-2) {
            $line_minus2 = $line;
        } elsif ($line_no==-1) {
            $line_minus1 = $line;
        }
        $line_no ++;
    } # foreach

    close FILE;
    #print "The number of line in $filename is $line_no.\n";



    #print $lines[0];
    #print $lines[5];
    #print $lines[$line_no-1];

    unless (-e "$outputrootname.pars") {
        system "mkdir $outputrootname.pars";
    }   
    #my @results;
    my $j;
    # Since hishrep and e are not used in following code, we can set fake ones.
    $lines[$line_no] = "......    0.0    $initial_hishape";
    $line_no ++;
    #print $line_no;
    $initial_hishape =~ s/\[/\\[/g;
    $initial_hishape =~ s/\]/\\]/g;
    for ($j=0;$j<$line_no;$j++)
    {
	my $pid = fork();
	if($pid){
		# parent
		push(@childs, $pid);
	}elsif ($pid == 0){
		# child
		my $hishrep;
		my $e;
		my $related_keywords;
		#print "HiPath -f ${outputrootname}.pre2.Faas/${i}_${j}.faa -t 1 > ${outputrootname}.pre2.Faas/${i}_${j}.hipath\n";
                #print "$lines[$j]\n";
		($hishrep,$e,$related_keywords)=split(/\s+/,$lines[$j]);
                #$hishapes[$j] = $hishape;
                my $hishape = $related_keywords;
		$related_keywords =~ s/\[//g;
		$related_keywords =~ s/\]//g;
                $related_keywords =~ s/\,\(//g;
                $related_keywords =~ s/\,\)//g;
		#my @helixindices=split(/\,/,$hishape);
                #print "RNAHeliCes -f $inputfilename -k 1000 -t 4 -R $related_keywords -e -z\n";
                #http://ubuntuforums.org/showthread.php?t=526553
		my $result = `RNAHeliCes -f $inputfilename -k 1000 -t 4 -R $related_keywords -e -z`;
		#print "result: $result";
		my @result_lines = split /\n/, $result;
		foreach my $result_line (@result_lines) {
                    #print $hishape;                
                    my $hishrep2;
		    my $e2;
		    my $hishape2;
                    my $partition2;
                    ($hishrep2,$e2,$hishape2,$partition2)=split(/\s+/,$result_line);
                    #if($hishape eq "[39]") {
                    #    print "$result_line\n";
                    #    print "$hishape vs. $hishape2: $hishrep eq $hishrep2\n";
                    #}
		    #if($hishrep eq $hishrep2) {
                    if ($result_line =~ m/$initial_hishape/) {
                        system ("echo '$result_line' > $outputrootname.pars/_.par"); 
                    } elsif (defined $hishape2 && $hishape eq $hishape2) {
                        #print "$result_line =~ m/$hishape/\n";#rw2.pars/0.par
                        #print "j= $j; $results[$j] = $result_line;\n";
                        #print $result_line;
                        #$results[$j]="'$result_line'";
                        system ("echo '$result_line' > $outputrootname.pars/$j.par"); 
                    }
		}

		exit 0;
	}else{
		die "couldnt fork: $!\n";
	}
    }


    #print "waitpid\n";
    foreach (@childs) {
	    my $tmp = waitpid($_, 0);
	    #print "done with pid $tmp\n";
    }

    #dump(@results);
    #print $results[0];
    #print join( ',', @results);
    #print $main::results;

    system ("> $outputrootname.pre2"); 
    system ("echo '$line_minus2' >> $outputrootname.pre2");  
    system ("echo '$line_minus1' >> $outputrootname.pre2");  
    #my $containsOpenChain = 0;
    for ($j=0;$j<$line_no;$j++)
    {
	#my $filename = "$outputrootname.pars/$j.par";
	#unless( open(FILE, $filename) ) {
	#    die "Error: couldn't open file '$filename': $!\n";
	#}
	#while ($line =<FILE>) {
	    #print $line;
	    #chomp($line);
             #print "$lines[$j] vs $initial_hishape\n";
	    if(!($lines[$j] =~ m/$initial_hishape/)) {
                system ("cat $outputrootname.pars/$j.par >> $outputrootname.pre2"); 
                #$containsOpenChain = 0;
                #print "true";
	    }
	#}
#	close FILE;
    }
    #if (!$containsOpenChain) {
    system ("cat $outputrootname.pars/_.par >> $outputrootname.pre2");
    #}
    #print "con=$containsOpenChain\n";
# }
 

# print "End of main program\n";
