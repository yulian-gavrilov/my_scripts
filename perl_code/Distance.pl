#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;

my ($inputFileMass) = $ARGV[0];
my ($inputFileCoord) = $ARGV[1];
#my $firstGroup = $ARGV[2]; # L8: 15,16     #44,70: 85,86,135,136
#my $secondGroup =$ARGV[3]; # L8: 161,162   #44.70: 231,232,281,282

my $numberOfBeadsInDomain = 146;
my $numberOfDomains = 2;
open(IN,$inputFileMass) or die "Error opening input (1) file. \n";
open(IN1,$inputFileCoord) or die "Error opening input (2) file. \n";

my ($FilePrefix) = $inputFileCoord =~ m/(.*)\.dat/g;


my $index_line = 0;
my $index_line1 = 0;

my @Mass = 0;

my %hash_res = ("ALA","15","VAL","43","SER","31","LEU","57","THR","45","ILE","57","CYS","47","PHE","91","TYR","107","TRP","130","ASN","58","MET","75","GLN","72","PRO","42","LYS","72","ASP","59","ARG","100","GLU","73","HIS","81");

# Read masses
# If Nrterm of Cterm plus 1 (H) or 17 (OH).
# Domains should be atthached to each other through Cterm
# Lys CB bead in the site of attachement should have 1 a.m.e less. Originally it is 72, so the error is small

while (my$currentLine = <IN>)
  {
  
        if ($currentLine =~m/^\s*\d+\s+CA\s+(?!PRO)(?!GLY)\s*/gi ){
            
            $Mass[$index_line] = 56; 
            
            if ($index_line == 0 or $index_line == $numberOfBeadsInDomain){$Mass[$index_line] = $Mass[$index_line]+1}  # if it is Nterm (domains are attached not through Nterm)       
            if ($index_line == ($numberOfDomains*$numberOfBeadsInDomain)-1){$Mass[$index_line] = $Mass[$index_line]+17} # if it is Cterm of the last domain        
           
            $index_line++;
        }
        
        
        if ($currentLine =~m/^\s*\d+\s+CA\s+PRO\s*/gi ){
            
            $Mass[$index_line] = 55; 
            
            if ($index_line == 0 or $index_line == $numberOfBeadsInDomain){$Mass[$index_line] = $Mass[$index_line]+1}         
            if ($index_line == ($numberOfDomains*$numberOfBeadsInDomain)-1){$Mass[$index_line] = $Mass[$index_line]+17}                                            
           
            $index_line++;                                             
        }
        
        
        if ($currentLine =~m/^\s*\d+\s+CA\s+GLY\s*/gi ){
            
            $Mass[$index_line] = 57; 
            
            if ($index_line == 0 or $index_line == 1*$numberOfBeadsInDomain){$Mass[$index_line] = $Mass[$index_line]+1}         
            if ($index_line == ($numberOfDomains*$numberOfBeadsInDomain)-1){$Mass[$index_line] = $Mass[$index_line]+17}            
            $index_line++;                                            
        }

        if ($currentLine =~m/^\s*\d+\s+CB\s+\w+\s*/gi ){
            
            my ($res) = $currentLine =~m/^\s*\d+\s+CB\s+(\w+)\s*/;
            #print  "$res\n";
            $Mass[$index_line] = $hash_res{$res};
            
            if ($index_line == 0 or $index_line == 1*$numberOfBeadsInDomain){$Mass[$index_line] = $Mass[$index_line]+1}         
            if ($index_line == ($numberOfDomains*$numberOfBeadsInDomain)-1){$Mass[$index_line] = $Mass[$index_line]+17}
           
            $index_line++;
        }

     
              
          
  }

#foreach ( @Mass){print $_."\n"; };
#my $MassSize = scalar @Mass;
#print "Size = $MassSize\n";


# Read Coordinates

my (@Xub,@Yub,@Zub);
my $index_line2 = 0;
my @Mass_group1;
my @Mass_group2;
my (@Xub_group1,@Yub_group1,@Zub_group1);
my (@Xub_group2,@Yub_group2,@Zub_group2);
my (@numeratorX_group1,@numeratorY_group1,@numeratorZ_group1);
my ($numeratorX_group1,$numeratorY_group1,$numeratorZ_group1);
my (@numeratorX_group2,@numeratorY_group2,@numeratorZ_group2);
my ($numeratorX_group2,$numeratorY_group2,$numeratorZ_group2);
my ($denominator_group1);
my ($denominator_group2);
my ($Xcm_group1,$Ycm_group1,$Zcm_group1);
my ($Xcm_group2,$Ycm_group2,$Zcm_group2);
my $Distance;


print "Please enter index numbers of beads in the first group (Format: A_B_C): ";
my $group1 = <STDIN>;
chomp($group1);
my @group1 = split (/_/,$group1);
#my @group1 = (15, 16);


print "Please enter index numbers of beads in the second group (Format: A_B_C): ";
my $group2 = <STDIN>;
chomp($group2);
my @group2 = split (/_/,$group2);
#my @group2 = (161, 162);

open(OUT, ">",$FilePrefix."_distance_g1_".$group1."_g2_".$group2.".dat") or die "Error opening pdb file for writing: $!\n";


#foreach ( @group1){print $_."\n"; };
#foreach ( @group2){print $_."\n"; };


while (my $currentLine = <IN1>)
  {
    
  #print  "$currentLine";
  
        if ($currentLine =~m/^\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*/gi ){
          ($Xub[$index_line2],$Yub[$index_line2],$Zub[$index_line2]) = $currentLine =~m/^\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*/;
          
           $index_line2++;
           
        #my $XubSize = scalar @Xub;
        #print "Size = $XubSize\n";
        
            if ($index_line2>($numberOfDomains*$numberOfBeadsInDomain)-1){
                
               #FIND COORDINATES OF THE CENTER OF MASSES OF TWO GROUPS AND CALCULATE DISTANCE BETWEEN THEM!#
               
                for (my $i=0;$i<(scalar @group1);$i++) {
                    $Mass_group1[$i] = $Mass[$group1[$i]-1];
                    $Xub_group1[$i] = $Xub[$group1[$i]-1];
                    $Yub_group1[$i] = $Yub[$group1[$i]-1];
                    $Zub_group1[$i] = $Zub[$group1[$i]-1];
                    }
                
                for (my $i=0;$i<(scalar @group2);$i++) {
                    $Mass_group2[$i] = $Mass[$group2[$i]-1];
                    $Xub_group2[$i] = $Xub[$group2[$i]-1];
                    $Yub_group2[$i] = $Yub[$group2[$i]-1];
                    $Zub_group2[$i] = $Zub[$group2[$i]-1];
                    }

                for (my $i=0;$i<(scalar @group1);$i++) {
                $numeratorX_group1[$i] = $Mass_group1[$i]*$Xub_group1[$i];
                $numeratorY_group1[$i] = $Mass_group1[$i]*$Yub_group1[$i];
                $numeratorZ_group1[$i] = $Mass_group1[$i]*$Zub_group1[$i];
                }
                
                for (my $i=0;$i<(scalar @group2);$i++) {
                $numeratorX_group2[$i] = $Mass_group2[$i]*$Xub_group2[$i];
                $numeratorY_group2[$i] = $Mass_group2[$i]*$Yub_group2[$i];
                $numeratorZ_group2[$i] = $Mass_group2[$i]*$Zub_group2[$i];
                }                
                
             
                foreach (@numeratorX_group1){$numeratorX_group1 += $_;}
                foreach (@numeratorY_group1){$numeratorY_group1 += $_;}
                foreach (@numeratorZ_group1){$numeratorZ_group1 += $_;}
                
                foreach (@numeratorX_group2){$numeratorX_group2 += $_;}
                foreach (@numeratorY_group2){$numeratorY_group2 += $_;}
                foreach (@numeratorZ_group2){$numeratorZ_group2 += $_;}
                
                foreach (@Mass_group1){$denominator_group1 += $_;}
                foreach (@Mass_group2){$denominator_group2 += $_;}


                $Xcm_group1 =  $numeratorX_group1/$denominator_group1;
                $Ycm_group1 =  $numeratorY_group1/$denominator_group1;
                $Zcm_group1 =  $numeratorZ_group1/$denominator_group1;

                $Xcm_group2 =  $numeratorX_group2/$denominator_group2;
                $Ycm_group2 =  $numeratorY_group2/$denominator_group2;
                $Zcm_group2 =  $numeratorZ_group2/$denominator_group2;


                $Distance = sqrt(($Xcm_group2-$Xcm_group1)**2+($Ycm_group2-$Ycm_group1)**2+($Zcm_group2-$Zcm_group1)**2);
                
                print OUT "$Distance\n";
                
                #print "$Xcm_group1 $Ycm_group1 $Zcm_group1\n";
                #print "$Xcm_group2 $Ycm_group2 $Zcm_group2\n";
                #print "Distance = $Distance\n";
                #print "\n";
                
                #print "$numeratorX_group1\n";
                #print "$denominator_group1\n";                
                
                #foreach (@Mass_group1){print $_."\n"; };
                #print "#\n";
                #foreach (@Xub_group2){print $_."\n"; };
                #print "###\n";
                
                
                @Xub=0;@Yub=0;@Zub=0;$index_line2=0;
                $denominator_group1 = 0;
                $denominator_group2 = 0;
                $numeratorX_group1 = 0;$numeratorY_group1 = 0;$numeratorZ_group1 = 0;
                $numeratorX_group2 = 0;$numeratorY_group2 = 0;$numeratorZ_group2 = 0;
                $Distance = 0;
            }

        }
        
    
#foreach ( @Xub){print $_."\n"; };
#my $XubSize = scalar @Xub;
#print "Size = $XubSize\n";
  }


close IN;
close IN1;
close OUT;