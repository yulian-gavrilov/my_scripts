#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;


my $start = time;

#require "/home_b/yulian/scripts/dihedrals_sub.pl";
require "/home_b/yulian/scripts/dihedrals_atan_sub.pl";


my ($x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3, $x4, $y4, $z4);

#$x1 = 2.656;   $y1 =    2.433;   $z1 =    2.885;
#$x2 = 2.621;   $y2 =    2.510;   $z2 =    2.762;
#$x3 = 2.744;   $y3 =    2.541;   $z3 =    2.673;
#$x4 = 2.736;   $y4 =    2.678;   $z4 =    2.603;


my ($inputFile) = $ARGV[0];
my $index_line = 0;
my (@X,@Y,@Z);
my $residues = 0;
my @AtomName;
my $frame;
my $resindexName;

my @frameNumber;
my $frameNumber_step = 0; 
my $newFrame = -1;

open(IN,$inputFile) or die "Error opening output file. \n";

my ($groFileName) = $inputFile =~ m/(.*)\.gro/g;

open(OUT, ">",$groFileName."_dihedrals_atan".".dat") or die "Error opening pdb file for writing: $!\n";
open(OUT1, ">",$groFileName."_dihedrals_atan_angles".".dat") or die "Error opening pdb file for writing: $!\n";


while (my$currentLine = <IN>)
  {
  $index_line++;
  
    if ($currentLine =~m/^\s+\d+\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi and $currentLine !~m/CA| C | N | O \s+/){
        
    ($resindexName, $AtomName[$residues], $X[$residues],$Y[$residues],$Z[$residues]) = $currentLine =~m/^\s+(\d+\w+)\s+(\w+)\s+\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
     
               #print  "$currentLine";
          #print  "$resindexName $AtomName[$residues]  $X[$residues]  $Y[$residues]  $Z[$residues]\n";
          #print  "$index_line\n";
         $residues++;
    }
    
    #(my $frame) = $currentLine =~m/\.*t=\s*(\d+\.\d+)$/;
    #(my $box) = $currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi;




    if ($currentLine =~m/\.*t=\s*(\d+\.\d+)$/gi){
        
     ($frame) = $currentLine =~m/\.*t=\s*(\d+\.\d+)$/;
     
        $newFrame++;
        #my ($frameNumber) = $currentLine =~m/\.*t=\s*(\d+\.\d+)$/;
        $frameNumber[$newFrame] = $frame;
        
        if ($frame >= $frameNumber_step){
            
          print "Frame:".$frame.". ";
          $frameNumber_step=$frameNumber_step+10000;
          my $duration_interm = time - $start;
          print "Running for: $duration_interm s\n";
          
        }
        
    }
    
    
 
    
    if ($currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi){
        
     #foreach (@X) {  print "$_\n";};
     
     
     
# PUT A FRAME NUMBER AS A FIRST COLI    
###########################################################
     #if ($frame==0 and $newFrame < 1){
     #  #print OUT "$frame   ";
     #  my $frame_print = sprintf ("%.1f",$frame);
     #  print  OUT "$frame_print       ";
     #  
     # }
     #else{
     # #print OUT "\n$frame   ";
     # my $frame_print = sprintf ("%.1f",$frame);
     # print  OUT "\n$frame_print       ";     
     # 
     #}
########################################################### 
      
     if ($frame==0 and $newFrame < 1){
     print OUT "   ";

       
     }
     else{
      print OUT "\n   ";
     
      
     }
     
########################################################### 

      
      
      
      for (my $i=0;$i<(scalar @X)-3;$i++){
      # $i=1;
      
      if ($frame==2 and $newFrame < 3){
       print OUT1 "$resindexName $AtomName[$i]-$AtomName[$i+1]-$AtomName[$i+2]-$AtomName[$i+3]\n";
      }

       my $out_dihedral = calc_dihedral ($X[$i],$Y[$i],$Z[$i],$X[$i+1],$Y[$i+1],$Z[$i+1],$X[$i+2],$Y[$i+2],$Z[$i+2],$X[$i+3],$Y[$i+3],$Z[$i+3]);

      
      } # for
      
       @AtomName = 0;@X= 0;@Y= 0;@Z= 0;
       $residues=0;
      
    } #if
    
    
  }
  
  
    
close IN;
close OUT;
close OUT1;

my $duration = time - $start;
print "Execution time: $duration s\n";






