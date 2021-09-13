#!/usr/bin/perl -w

use strict;
use warnings;
use Math::Trig;


$|=1;

my@line = 1;

my $filenameList;
my $datfile = $ARGV[0];
my $Nframes = $ARGV[1];#10000

if(scalar(@ARGV) == 3){
$filenameList=$ARGV[2];
}
else
{
#############################################
use File::Basename;
$filenameList = "list"; # Traj file # $ARGV[0];

my $f = "*"; 
my $dir = dirname ($f);
my $name;
my $file;

my ($outputFile) = "list";   

open(OUTPUT_FILE_HANDLE,">",$outputFile) or die "Error opening output file. \n";

foreach $file (glob ("$dir/*")){
    if (-f $file){
        $name = basename ($file);
        #print $name, "\n";
        if ($name ne "list"){
        my $Line = sprintf("%s", $name);
        print  OUTPUT_FILE_HANDLE "$Line\n";
        }
    } 
}

close (OUTPUT_FILE_HANDLE);

############################################

}

 
 open (INPUT_FILES_HANDLE, $filenameList) or die "can't open list of degrees $!"; #open input data file and degrees file
 
 while (my $filenameLine=<INPUT_FILES_HANDLE>) {

    chomp $filenameLine;

    my $file = $filenameLine; # parameter in command line from angle_conditions.txt
    
      my $commandLine = "~yulian/scripts/Elec_from_AllElectrost_per_bead.pl  $file $datfile $Nframes"; #the string that will be inserted in the command line 
    
    #  print ("$commandLine\n");
     
     `$commandLine`;
  
 }
 close INPUT_FILES_HANDLE;
