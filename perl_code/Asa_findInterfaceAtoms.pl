#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;

my ($inputFile_init) = $ARGV[0];
my ($inputFile) = $ARGV[1];
my ($inputFile_gro) = $ARGV[2];


open(IN,$inputFile) or die "Error opening input file. \n";
open(IN1,$inputFile_init) or die "Error opening input file. \n";
open(IN2,$inputFile_gro) or die "Error opening input file. \n";


my ($FilePrefix) = $inputFile_gro =~ m/(.*)\.gro/g;
my ($frame) = $inputFile_gro =~ m/all_frames_mol_center_\w+\d+_\d_(\d+)\.gro/g; # all_frames_mol_center_tUb94_1_0.gro
$frame = $frame*100;
open(OUT, ">",$FilePrefix."_interface".".gro") or die "Error opening pdb file for writing: $!\n";
print  OUT " t=   $frame.000\n";

my $index_line = 0;
my $index_line1 = 0;
my @arg1;my @arg2;my @arg3;
my @arg4;my @arg5;
my @arg6;my @arg7;
my @arg8;my @arg9;
my @arg1b;my @arg2b;my @arg3b;
my @arg4b;my @arg5b;
my @arg6b;my @arg7b;
my @arg8b;my @arg9b;

#my $currentLine;

while (my$currentLine = <IN>)
  {
          #if ($currentLine =~m/^\s+\d+\s+\w+\s+\w+/gi ){
  ($arg1[$index_line],$arg2[$index_line],$arg3[$index_line],$arg4[$index_line],$arg5[$index_line],$arg6[$index_line],$arg7[$index_line],$arg8[$index_line],$arg9[$index_line]) =
+ $currentLine =~m/^ATOM\s+(\d+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\w*\s*/gi;
#+ $currentLine =~m/^ATOM\s+(\d+)\s+(\w+)\s+(\w+)\s*\w*\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\w*\s*/gi;
  
  #print  "$currentLine  $index_line\n";
  $index_line++;
  #print  "$arg2[$index_line-1] $arg3[$index_line-1]\n";
  }

#print  "$index_line\n";

while (my$currentLine1 = <IN1>)
  {
         #if ($currentLine =~m/^\s+\d+\s+\w+\s+\w+/gi ){
  ($arg1b[$index_line1],$arg2b[$index_line1],$arg3b[$index_line1],$arg4b[$index_line1],$arg5b[$index_line1],$arg6b[$index_line1],$arg7b[$index_line1],$arg8b[$index_line1],$arg9b[$index_line1]) =
+ $currentLine1 =~m/^ATOM\s+(\d+)\s+(\w+)\s+(\w+)\s+\w*\s*(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*\w*\s*/gi;
          #print  OUT " out b $arg9b[$index_line1]\n";
  $index_line1++;
  #print  "$arg2b[$index_line1-1] $arg3b[$index_line1-1]\n"; 
  }

my $z = 0;
#my $i = 93;

#  for (my $i=0;$i<(scalar @arg8b);$i++) { 
   #if ( grep {$_+$arg8b[$i]==$_} @arg8){$z++};

    #if ($arg2[$i] eq $arg2b[$i]){$z++};
#  }

my @arg2ch;my @arg3ch;my @arg4ch;
my $chose_counter = 0;
  
  for (my $i=0;$i<(scalar @arg1b);$i++) { 
    for (my $j=0;$j<(scalar @arg1b);$j++) {
      #if ($arg2[$i] eq $arg2b[$j]  and  $arg3[$i] eq $arg3b[$j] and  $arg4[$i] eq $arg4b[$j] and $arg8b[$j] > 0 and $arg8[$i] == 0){
      #if ($arg2[$i] eq $arg2b[$j]  and  $arg4[$i] eq $arg4b[$j] and $arg8b[$j] > 0 and $arg8[$i] == 0){
     # if ($arg2[$i] eq $arg2b[$j]  and  $arg4[$i] eq $arg4b[$j] and $arg8[$i] == 0 and $arg8b[$j] > 0){
      if ($arg2[$i] eq $arg2b[$j]  and  $arg4[$i] eq $arg4b[$j] and $arg8b[$j] > $arg8[$i]){
       
        #print "$arg2[$i] eq $arg2b[$j]  and  $arg3[$i] eq $arg3b[$j] and  $arg4[$i] eq $arg4b[$j] ___ $arg8b[$j]   $arg8[$i] \n";
        
        $arg2ch[$chose_counter] = $arg2[$i];
        $arg3ch[$chose_counter] = $arg3[$i];
        $arg4ch[$chose_counter] = $arg4[$i];
        #print "$arg2ch[$chose_counter] $arg3ch[$chose_counter] $arg4ch[$chose_counter] \n";
        $chose_counter++;
        $z++;
        #print "$arg8b[$j] $arg8[$i]\n";
      }
      else{
#        print  OUT "$currentLine";
      }
    }
  }

   #print "$z\n";
   
   #foreach (@arg2ch) {  print "$_\n";};
my $zz = 0;
while (my$currentLine2 = <IN2>)
  {
     if ($currentLine2 =~m/^\s+\d+(?!SOL)\w{3}\s+(?!H)\w+\s*\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi 
        or
        + $currentLine2 =~m/^\s+\d+N\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi
        or
        + $currentLine2 =~m/^\s+\d+C\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi
        ){
      #print  "$currentLine2";
      my($resNumber,$resName,$atomName, $X,$Y,$Z) = $currentLine2 =~m/^\s+(\d+)(\w+)\s+(\w+)\s+\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
      
      if (length($resName)>3){$resName = substr($resName, 0, -1);}
      
      #print  "$resNumber $resName $atomName \n";
      
      
   for (my $i=0;$i<(scalar @arg2ch);$i++) { 
    
      if ($atomName eq $arg2ch[$i]  and  $resName eq $arg3ch[$i] and $resNumber eq $arg4ch[$i]){
         $zz++;
         #chomp($currentLine2);
         print  OUT "$currentLine2";
         
      }
   
     }  
     }
     elsif($currentLine2 =~m/^\s*\d+SOL\s+\w+\d*\s*\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi){
       print  OUT "$currentLine2";
      }
     
    
  }
   #print "$zz\n";

#  grep {$_ eq $resNumber} @arg4ch and grep {$_ eq $resName} @arg3ch and grep {$_ eq $atomName} @arg2ch

close IN;
close IN1;
close IN2;
close OUT;

