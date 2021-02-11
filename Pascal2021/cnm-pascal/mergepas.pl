#!/usr/bin/perl
# mergepas.pl -- 
# J C Nash (C) 2021
# Use includes to build self-contained pascal file

use strict;
use Cwd; # ??
use Cwd 'realpath'; # ?? 
use File::Basename; # module that parses filenames into path, name, extension etc.
use File::Copy;
# use File::HomeDir; # Need to install libfile-homedir-perl

my ($ii, $jj, $kk, $ifn, $tmp, $tfname, $oldpn);
my ($inpath, $fextn, $fnroot);
my ($tfnam, $numfiles, $pname);
my (@lyns, @persname, $d8, @ymd, @tfile);
my ($trln, $titl, $ccwd, $cmd, $nfn, $pfn);
my  $inkstr="\{\$I ";

print "mergepas.pl -- (C) J C Nash 2021\n";
print "\n";
## $tmp =<STDIN>;

my $numclargs = @ARGV; # get the number of arguments on the command line
if ($numclargs==0) { #if no command line parameters, then display a message about the program
   print "This program was written by John C. Nash\n";
   print "     18 Spyglass Ridge, Stittsville, Ontario, K2S 1R6\n";
   print "       nashjc\@uottawa.ca\n";
   print "\n";
   print "Usage:\n";
   print "     mergepas.pl infile\n";
   print "infile name =";
   $ifn = <STDIN>;
} else { 
   $ifn = $ARGV[0]; # filename root is first parameter
}

chomp $ifn;
print "working with input filename = $ifn\n";
$ifn=realpath($ifn); # expand to full path 
($fnroot,$inpath,$fextn) = fileparse($ifn);

print "full input filename = $ifn\n";
# print "filename root = $fnroot\n";
#print "filename extension is $fextn\n";
print "path is $inpath\n";

open (INF, $ifn) or die("could not open $ifn for read");
@lyns = <INF>; # get the data
close INF;

# print @lyns;

open (OUTF, ">./out.pas") or die("Outfile failure");
# now loop over lines
  print "Search for $inkstr \n";
  foreach $tfname (@lyns) {
    $tfname=ltrim($tfname); # remove left white space
    $kk=index($tfname, $inkstr);
#    print $kk," ";
    if ( $kk == 0) {
        print "****Include found: ",$tfname;
        $tfname = substr($tfname, 3); 
#        chomp $tfname;
        $tfname=ltrim($tfname);
        $ii=index($tfname, "\}");
        $tfname=substr($tfname, 0, $ii);
        print "Including:$tfname\n";
        open (INF, $tfname);
        @tfile = <INF>;
        close INF;
        foreach $tmp (@tfile){
           print OUTF $tmp;
        }
    } 
    else {
       print OUTF $tfname;
    }



  }
close(OUTF);

die("Done");

#===============================================
# Perl trim function to remove whitespace from the start and end of the string
sub trim {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim {
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim {
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
#===============================================


