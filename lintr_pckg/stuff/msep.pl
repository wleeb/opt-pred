#!/usr/bin/perl
#
        use v5.10;
        use strict;
        msep();
#
#
#
#
#
    sub msep{
#
        my @lines, my @ldescr, my @lfunc;
        my $ldescr_ref, my $lines_ref, my $lfunc_ref;
        my $n, my $m1, my $m2, my $nnc, my $nfuncs;
        my $i, my $j1, my $j2;
        my $fout;
#
#        . . . input parameters
#
        my $fin = $ARGV[0];
        my $fkey = $ARGV[1];

        ($lines_ref, $n) = msep_read($fin);
        @lines = @$lines_ref;


        print("\n    lines, total=\n          $n \n");

#
#       find lines beginning and ending the description
#
        ($ldescr_ref, $nnc) = msep_postmatch(\@lines,$n,0,"%%%%%%%%%%%%%");
        @ldescr=@$ldescr_ref;
#
        $m1 = @ldescr[1];
        $m2 = @ldescr[$nnc];

        print("\n    first description, line= \n          $m1 \n");
        print("\n    last description, line=  \n          $m2 \n");

#
#       print testcode and description
#
        msep_prin2file("testcode.m",\@lines,1,$m1-1);
        msep_prin2file("description.m",\@lines,$m1,$m2);

#
#       find all lines with function declarations,
#       after the last description line
#
        ($lfunc_ref, $nfuncs) = msep_postmatch(\@lines,$n,$m2,"      function");
        @lfunc=@$lfunc_ref;
        @lfunc[$nfuncs+1]=$n;

#
#       extract each function's name and print to file
#
        for ($i=1; $i <= $nfuncs; $i=$i+1){
#
        $j1 = @lfunc[$i];
#
#        last line of function
#
        $j2 = @lfunc[$i+1];
        $j2=msep_firstprior(\@lines,$j2,"      end");
        $j2=$j2+1;
#
        $fout = msep_findname(\@lines,$j1,$fkey);
        msep_prin2file($fout,\@lines,$j1,$j2);
        }

        return;
    }
#
#
#
#
#
    sub msep_firstprior{
#
#       Returns the line number containing first occurrence
#       of sym prior to specified line number j
#
        my $lines_ref = $_[0];
        my @lines = @$lines_ref;
#
        my $j = $_[1];
        my $sym = $_[2];

###        print "\n j= $j \n";
###        print "\n symbol= $sym \n";

        my $i;
        for ($i = $j; $i > 0; $i = $i-1){
#
        if ($lines[$i] =~ m/$sym/){
#
        last;
        }
        }

        return($i);

    }
#
#
#
#
#
    sub msep_prin2file{
#
#       prints lines in specified range to file named fun
#
        my $fun=$_[0];
        my $lines_ref=$_[1];
        my @lines = @$lines_ref;
        my $j1=$_[2];
        my $j2=$_[3];

###        print "\n function name = \n $fun \n";
#
#        create new file for function
#
        my $fh;
        open ($fh,">$fun") or die;
#
#       print lines to new file
#
        for (my $i=$j1; $i <= $j2; $i=$i+1){
#
        print $fh "$lines[$i]";
        }
        close $fh;

        return
    }
#
#
#
#
#
    sub msep_findname{

        my $lines_ref = $_[0];
        my @lines = @$lines_ref;
        my $j = $_[1];
        my $fkey = $_[2];
#
#       concatenate function line and next two lines
#
        my $str=@lines[$j] . @lines[$j+1] . @lines[$j+2];
###        print "\n $str \n";
#
#       determine name of the m-file
#
        my $fname;
        if ( $str =~ /.*($fkey.*)\(/ ) {
#
        $fname=$1;
        }
        $fname = $fname . ".m";
###        print "\n $fname \n";

        return($fname);
    }
#
#
#
#
#
    sub msep_postmatch{
#
#       Finds all line numbers after $init which contain
#       the specified string $sym
#
        my $lines_ref=$_[0];
        my @lines = @$lines_ref;
#
        my $n = $_[1];
        my $init = $_[2];
        my $sym = $_[3];

###        print "n = $n \n";
###        print "sym = $sym \n";

        my @lmatch;
        my $nmatch = 0;
        my $i;
        for($i = $init; $i <= $n; $i = $i+1) {
#
        if ($lines[$i] =~ m/$sym/){
#
        $nmatch = $nmatch+1;
        @lmatch[$nmatch] = $i;

###        print "\n line $i: \n";
###        print @lines[$i];
        }

        }

        return(\@lmatch, $nmatch);

    }
#
#
#
#
#
    sub msep_read{
#
#       Reads file with specified name; returns 
#       array with each line, and the line count
#
        my $fname=$_[0];
        my $n;
        my $file;
        my @lines;

        open($file,$fname);
#
#        copy lines to new array
#
        while (<$file>){
#
        $n = $.;
        @lines[$n]=$_;
        }
#
        return (\@lines, $n);
    }
