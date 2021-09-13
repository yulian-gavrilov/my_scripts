#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;

require "/home_b/yulian/scripts/dihedrals_sub.pl";
#require "/home_b/yulian/scripts/dihedrals_atan_sub.pl";


my ($x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3, $x4, $y4, $z4);

$x1 = 2.656;   $y1 =    2.433;   $z1 =    2.885;
$x2 = 2.621;   $y2 =    2.510;   $z2 =    2.762;
$x3 = 2.744;   $y3 =    2.541;   $z3 =    2.673;
$x4 = 2.736;   $y4 =    2.678;   $z4 =    2.603;


my ($inputFile) = $ARGV[0];
my $index_line = 0;
my (@X,@Y,@Z);
my $residues = 0;
my @AtomName;



open(IN,$inputFile) or die "Error opening output file. \n";

my ($groFileName) = $inputFile =~ m/(.*)\.gro/g;

open(OUT, ">",$groFileName."_dihedrals_cos".".dat") or die "Error opening pdb file for writing: $!\n";


while (my$currentLine = <IN>)
  {
  $index_line++;
  
    if ($currentLine =~m/^\s+\d+\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi and
       +$currentLine !~m/CA| C | N | O \s+/){
        
    ($AtomName[$residues], $X[$residues],$Y[$residues],$Z[$residues]) = $currentLine =~m/^\s+\d+\w+\s+(\w+)\s+\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
     
               #print  "$currentLine";
          #print  "$AtomName[$residues]  $X[$residues]  $Y[$residues]  $Z[$residues]\n";
         
         $residues++;
    }
   }
    
    for (my $i=0;$i<(scalar @X)-3;$i++) {
      
       print OUT "$AtomName[$i]-$AtomName[$i+1]-$AtomName[$i+2]-$AtomName[$i+3]     ";
       my $out_dihedral = calc_dihedral ($X[$i],$Y[$i],$Z[$i],$X[$i+1],$Y[$i+1],$Z[$i+1],$X[$i+2],$Y[$i+2],$Z[$i+2],$X[$i+3],$Y[$i+3],$Z[$i+3]);
       
    }
    
    
close IN;
close OUT;



#my $out = calc_dihedral ($x1, $y1, $z1, $x2, $y2, $z2, $x3, $y3, $z3, $x4, $y4, $z4);





