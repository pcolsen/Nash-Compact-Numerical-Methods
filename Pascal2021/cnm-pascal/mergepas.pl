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

my ($ii, $jj, $ofn, $tmp, $oldpn);
my ($outpath, $fextn, $fnroot);
my ($tfnam, $numfiles, $pname);
my (@imfiles, @persname, $d8, @ymd, @tfiles);
my ($trln, $titl, $ccwd, $cmd, $nfn, $pfn);

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
print "working with input filename = $ofn\n";
$ifn=realpath($ifn); # expand to full path 
($fnroot,$inpath,$fextn) = fileparse($ifn);

print "full input filename = $ifn\n";
# print "filename root = $fnroot\n";
#print "filename extension is $fextn\n";
print "path is $inpath\n";

# Now get file names of image files -- jpg only (lower case only)!
print "Image files\n";
@imfiles=();
@tfiles=glob("*.jpg");
$numfiles=@tfiles;
print "$numfiles of type .jpg\n";
@imfiles=(@imfiles,@tfiles);
$numfiles=@imfiles;
print "$numfiles total jpg image files\n";	

# sort the files

@tfiles = sort @imfiles;
@imfiles = @tfiles; 
# Now parse out the filenames to build tex file
 
# loop over files

# copy tex preamble into outfile
my $preamblefile = glob('~/bin/FrohnTexPreamble.txt');
copy($preamblefile, $ofn) or die "Tex preamble copy failed: $!";

# ?? should arrange that page numbers start AFTER preamble

# Note: preamble ASSUMED to end with newpage!

## Now use the title 

if (-f "\@TITLE\@.TXT") { # there is a title -- use it
   open (TFILE, "\@TITLE\@.TXT");
   $titl = <TFILE>;
   print "Title: $titl\n";
   close TFILE;
 } else {
   print "No title file";
   $titl = "   ";
 }

# open file for output
open (OFILE, ">>", $ofn);
print "File $ofn is open\n";

##$tmp = <STDIN>;

$ii = 0; 
print OFILE "\n\n\n";
print OFILE "\\chapter*{$titl}\n\n";

foreach $tfnam (@tfiles) {
      print OFILE "\\begin{figure}[H]\n";
      print OFILE "\\flushleft\n";
      print OFILE "File: $tfnam\\\\\n";
      $trln = `rdjpgcom $tfnam`;
      ## 190713 -- try to fix double \\
      $trln =~ s/\\\\/\\/g; # replace \\ with \ in comments
      $trln =~ s/\\000//g; # get rid of \000 in comments
      print OFILE "\\fbox{\\includegraphics[width=0.45\\textwidth,left]{$tfnam}} \n";

      print OFILE "\\end{figure}\n";
      print OFILE "\n\\vspace{2.5mm}\n";

      print OFILE "$trln\n\n";
      print OFILE "\\vspace{3mm}\n";
#      print OFILE "\\centering\n";
# Need to do stuff for dots and spaces, hence grffile etc. ??
      print OFILE "\n\n";
}

#  $tmp = <STDIN>;

# copy file end to outfile
print OFILE "\n\n\n\\end{document}\n";
close OFILE;

print "rename out.tex to ?";

$nfn = <STDIN>;
chomp($nfn);

if (length($nfn) == 0) { $nfn=$fnroot };
move($fnroot, $nfn) or die("Rename failed"); 

## print "Now process to pdf ";
## $tmp = <STDIN>;
$cmd = "pdflatex ./".$nfn;
print "Command:".$cmd."\n";
## $tmp = <STDIN>;

system($cmd);

$ii=index($nfn, ".");
$pfn = substr($nfn, 0, $ii).".pdf";
print "\n\n";
print "PDF file=",$pfn,"\n";
## $tmp = <STDIN>;


$cmd="evince ./".$pfn;
system($cmd);

print "Upload?\n";
 $tmp = <STDIN>;

$cmd = "scp ".$pfn.' mnash@nashinfo.com:/home/mnash/public_html/FrohnAds/';
system($cmd);

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


