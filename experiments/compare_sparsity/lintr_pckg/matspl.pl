#!/usr/bin/perl
#
        use v5.10;
        main();
#
#
#
#
#
        sub main
{
        my @lfuns;
        my @lines;

#
#        . . . input parameters
#

###        $infile = "file.m";
###        $mag = "unre_";
###        $del = "       function";




        $infile=$ARGV[0];
        $mag=$ARGV[1];
        $del=$ARGV[2];
#

        open($infile,$infile);

        print "\n";
        print "input file: $infile \n";
        print "magic word: $mag \n";
        print "function delimiter: $del \n";

#
#        copy lines to new array
#
        while (<$infile>)
    {
        $num_line = $.;
        $cur_line=$_;
        @lines[$num_line]=$cur_line;
    }
        $nlines=$num_line;

#
#        find last line of test code
#
        $test_end = "%%%%%%%%%%%%";
        for( my $i=1; $i <= $nlines; $i=$i+1 )
    {
        if ($lines[$i] =~ m/$test_end/){
            $last_test = $i-1;
            last;
        }
    }

#
#        print test code to file
#
        open ($out_test,">test_code.m") or die;
        for ( my $i=1; $i <= $last_test; $i=$i+1 ){
            print $out_test "$lines[$i]";
        }

        close $out_test;



#
#        after that, find line of first function occurrence
#
        open($out_code,">code.m") or die;

        for( my $i=$last_test+1; $i <= $nlines; $i=$i+1 )
    {
        if ($lines[$i] =~ m/$del/){
            last;
        }
        print $out_code "$lines[$i]";
    }


###        print "$last_test \n";
###        print "$first_fun \n";


#
#        find line numbers starting each function
#
        $nfun=0;
        for( my $i=$last_test+1; $i <= $nlines; $i=$i+1 )
    {
        $cur_line=$lines[$i];
        if ($cur_line =~ /$del/){
            $nfun=$nfun+1;
            $lfuns[$nfun]=$i;
        }
    }
        $lfuns[$nfun+1]=$nlines + 1;
        close $infile;


        print "\n";
        print "number of lines = $nlines \n";
        print "number of functions = $nfun \n";

#
#        print each function to a separate file and to global file
#

        $ndone=0;
        for ( my $i=1; $i <= $nfun; $i=$i+1 )
    {
        $k=$lfuns[$i];
        $str = $lines[$k] . $lines[$k+1];
#
#        determine the name and make new file
#
        if ( $str =~ /.*($mag\S*)\(/ ) 
    {
        $name=$1;
        $name = $name . ".m";
        open ($out_file,">$name") or die;

###        print "$name \n";

#
#        print the function text to the new file
#
        $nlines=$lfuns[$i+1]-$lfuns[$i];
        for (my $j=$lfuns[$i]; $j <= $lfuns[$i+1]-1; $j=$j+1) 
    {
        print $out_file "$lines[$j]";
        print $out_code "$lines[$j]";
    }
        $ndone=$ndone+1;
        close $out_file;
    }

    }

        print "number of files created = $ndone \n";

        exit;

}
