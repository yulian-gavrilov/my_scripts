#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use List::Util qw( min max );
use POSIX;
use Data::Dumper;
use Data::Dumper qw(Dumper);
#my $min = min @numbers;
#my $max = max @numbers;
#Script 
# Input file should be gro file that includes required atoms of the protein and
# water atoms for each frame.
# In gro file coordinates are in nm
#
my $start = time;


my ($inputFile) = $ARGV[0];
my $cutoffDist = $ARGV[1]; #0.3; # 3 Ang
my $cutoffDist2= $ARGV[2]; # 0.4 4 Ang. Max distance between two WM.
#my $last_bin = $ARGV[2]; #100?
#my $side = $ARGV[3]; #box side 10 nm?
#my $max_sphere = $ARGV[2]; #

#my $system = $ARGV[2]; # peptide20_WT 
#my $protein_res = 99;
open(IN,$inputFile) or die "Error opening output file. \n";

my ($groFile) = $inputFile =~ m/(.*)\.gro/g;

open(OUT1, ">",$groFile."_cutoff_".$cutoffDist."_nm_angles".".dat") or die "Error opening gro file for writing: $!\n";

#open(OUT_ALL, ">",$groFile."_allDist.dat") or die "Error opening gro file for writing: $!\n";

#open(OUT_HIST, ">",$groFile."_cutoff_".$cutoffDist."nm_maxSphere".$max_sphere."_rdf.dat") or die "Error opening output file for writing: $!\n";
#open(OUT_BINS, ">",$groFile."_BINS.dat") or die "Error opening output file for writing: $!\n";


# test speed#

my $index_line = 0;
my $index_line1 = 0;
my (@Xprot,@Yprot,@Zprot);
my (@waterIndex,@Xw,@Yw,@Zw);
my $waterNumb= 0;
my @distLow = 0;
my @distHigh = 0;
my @allWater = 0;
my $newFrame = -1;
my @frameNumber;
my $residues;
my $residues1;
my $resNumber = 0;
my $frameNumber_step = 0;  
$residues = 0;
$residues1 = 0;
my $test;
my @waterDist;
my @waterDistALL;
my $index_waterDistALL = 0;
my $index_line_waterDist = 0;
my @XYZ;
my(@Xw_cut_unique,@Yw_cut_unique,@Zw_cut_unique);
my $waterNumbFixed;
my $waterNumbFixed_sum = 0;
my ($boxX, $boxY, $boxZ);
#my $maxbin= $last_bin;
#my $box_half = $side; #*0.5;

##my $dr = $box_half/$maxbin;    # Bin width
#my $dr = 0.01;
#my $maxbin = $max_sphere/$dr;#5;
##my $maxbin = $cutoffDist; #/$dr;#5;



my @hist;
#my $hist = [0]*($maxbin+1);
my $bin;
my (@Xw_cut,@Yw_cut,@Zw_cut, @XYZ_cut);
my$xx; my $yy; my$zz;

my $N_tetra =0;
my @Xtetra_central;
my @Xtetra=0;
my $rij;
my @rij;
my (@calcAngleX,@calcAngleY,@calcAngleZ);
my (@r1_X,@r1_Y,@r1_Z);
my $cos_sum=0;
my @Qtet;
my @keys;
my @waterDistALL_round;
#for (my $i=0;$i<$maxbin;$i++){$hist[$i] = 0;}

my @XYZwa;
my $waterAtomNumb = 0;
my (@Xwa,@Ywa,@Zwa);
my $waterOWNumb = 0;
my @atom;
my $counter1 = 0;
my (@CoordForAngles, @CoordForAnglesX,@CoordForAnglesY,@CoordForAnglesZ);
my @AtomTest;
my @waterIndexTest;
my(@Xw_unique,@Yw_unique,@Zw_unique);
my @CoordForAngles_noRep;
my @CoordForAngles_noRep1;
my @indexToRemove;
my $indexCorrection = 0;
my @secondWater;
my $length;
my @cos_four;
my @angles_four;
my @angle_min;

my $start1;
my $start2;
my $start3;
my $start_run;
my $end_run;
my $run_time;

while (my$currentLine = <IN>)
  {
  $index_line1++;
  
  
    if ($currentLine =~m/^\s+\d+(?!SOL)\w{3}\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi or
        + $currentLine =~m/^\s+\d+N\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi or
        + $currentLine =~m/^\s+\d+C\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi ){


     ($resNumber, $Xprot[$residues],$Yprot[$residues],$Zprot[$residues]) = $currentLine =~m/^\s+(\d+)(?!SOL)\w{3}\w*\s+\w+\s+\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
     
      $residues++;
       
      #print "$Xprot[$residues-1]\n";  
       
              
    }
    
     
      if ($currentLine =~m/^\s*\d+SOL\s+\w+\d*\s*\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi){
             
         
          
            #if ($counter1>100) {last;}
            #print "$currentLine";
            #$counter1++;
             
             ($Xwa[$waterAtomNumb],$Ywa[$waterAtomNumb],$Zwa[$waterAtomNumb]) = $currentLine=~m/^\s*\d+SOL\s+\w+\d*\s*\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
             ($waterIndex[$waterAtomNumb],$atom[$waterAtomNumb] ,$XYZwa[$waterAtomNumb]) = $currentLine=~m/^\s*(\d+)SOL\s+(\w{2})\d*\s*\d+\s+(\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+)$/;
             
              

              #print "$atom $waterIndex[$waterAtomNumb]  $XYZwa[$waterAtomNumb]  \n";

             
            $waterAtomNumb++;
            
      }
                    

  if ($currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi){ # end of a frame    

         ($boxX, $boxY, $boxZ) = $currentLine=~m/^\s*\-?(\d+\.\d+)\s+\-?(\d+\.\d+)\s+\-?(\d+\.\d+)\s*/;


             for (my $i=0;$i<(scalar @Xprot);$i++) { #(scalar @Xprot)
              
              $start_run =  time();
               
               for (my $j=0;$j<(scalar @Xwa);$j=$j+3) { #(scalar @Xprot)

              # choose WM based on cutoff for OW atoms
              
              #print "$waterIndex[0]  $atom[0]\n";
              
              #my $test = sqrt(($Xwa[$waterOWNumb]-$Xprot[$i])**2+($Ywa[$waterOWNumb]-$Yprot[$i])**2+($Zwa[$waterOWNumb]-$Zprot[$i])**2);
              #print "$test\n"; 
              
                if (sqrt(($Xwa[$j]-$Xprot[$i])**2+($Ywa[$j]-$Yprot[$i])**2+($Zwa[$j]-$Zprot[$i])**2)<=$cutoffDist){
              
                   #print "$waterIndex[$j]  $atom[$j] $XYZwa[$j]\n";
                   
                   # Make @CoordForAngles unique later (no repeats of WM coordinates) and spleat on X, Y, Z:
                   
                   push (@CoordForAngles, $XYZwa[$j], $XYZwa[$j+1], $XYZwa[$j+2]);
                   
                   
                   # Just for checking:
                   #push (@CoordForAnglesX, $Xwa[$j], $Xwa[$j+1], $Xwa[$j+2]);
                   #push (@CoordForAnglesY, $Ywa[$j], $Ywa[$j+1], $Ywa[$j+2]);
                   #push (@CoordForAnglesZ, $Zwa[$j], $Zwa[$j+1], $Zwa[$j+2]);
                   
                   #push (@AtomTest, $atom[$j], $atom[$j+1], $atom[$j+2]);
                   #push (@waterIndexTest, $waterIndex[$j], $waterIndex[$j+1], $waterIndex[$j+2]);
                         
     

                   
                   #$CoordForAngles[$k];
                   #$k++;
                 }
               
               }
              #$waterOWNumb = $waterOWNumb+3;
        
             } #for (my $i=0;$i<(scalar @Xprot);$i++)
         
      
 # Take out repetitions of water trajectories
 
my @CoordForAngles_noRep = do { my %seen; grep { !$seen{$_}++ } @CoordForAngles };      
#foreach (@CoordForAngles_noRep ){print "$_\n";}; 
 
      # Slow version of my @CoordForAngles_noRep = do { my %seen; grep { !$seen{$_}++ } @CoordForAngles }; 
      #   @CoordForAngles_noRep1 = @CoordForAngles;
      #              
      #     for (my $i=0;$i<(scalar @CoordForAngles);$i++) { #(scalar @CoordForAngles)
      #        
      #        #print "$CoordForAngles[$i]\n";
      #        
      #        for (my $j=$i+1;$j<(scalar @CoordForAngles);$j++) {
      #        
      #        
      #          if ($CoordForAngles[$i] eq $CoordForAngles[$j]){
      #        
      #           #print "$i $j\n";
      #           #print "$CoordForAngles[$i] | $CoordForAngles[$j]\n";
      #   
      #           
      #           $CoordForAngles_noRep1[$j] = "to delete";
      #           
      #   
      #          }
      #        
      #       #print "$CoordForAngles[$i]\n";
      #     
      #        }
      #      
      #     }
      #     
      #
      #     
      #print "TEST1\n";
      #$end_run = time();
      #$run_time = $end_run - $start_run;
      #print "Job took $run_time seconds\n";
      #
      #
      #for (my $i=0;$i<(scalar @CoordForAngles_noRep1);$i++) {
      # 
      # if ($CoordForAngles_noRep1[$i] ne "to delete"){
      #        
      #        push (@CoordForAngles_noRep,$CoordForAngles_noRep1[$i])
      # }
      #       
      #}
      

      
      #foreach (@CoordForAngles_noRep){print "$_\n";};
      #foreach (@CoordForAngles){print "$_\n";};

      
      for (my $i=0;$i<(scalar @CoordForAngles_noRep);$i++){
            
           ($Xw_cut_unique[$i],$Yw_cut_unique[$i],$Zw_cut_unique[$i]) = split  /\s+/, $CoordForAngles_noRep[$i];
       }
      
  
      
      for (my $i=0;$i<(scalar @CoordForAngles_noRep-3);$i=$i+3) { #(scalar @Xprot)
              
               for (my $j=$i+3;$j<(scalar @CoordForAngles_noRep);$j=$j+3) {
               
                     $xx = ($Xw_cut_unique[$i]-$Xw_cut_unique[$j]);
                     $yy = ($Yw_cut_unique[$i]-$Yw_cut_unique[$j]);
                     $zz = ($Zw_cut_unique[$i]-$Zw_cut_unique[$j]);

                     if ($xx <-$boxX/2) {$xx=$xx+$boxX}
                     if ($xx >$boxX/2) {$xx=$xx-$boxX}

                     if ($yy <-$boxY/2) {$yy=$yy+$boxY}
                     if ($yy >$boxY/2) {$yy=$yy-$boxY}
                     
                     if ($zz <-$boxZ/2) {$zz=$zz+$boxZ}
                     if ($zz >$boxZ/2) {$zz=$zz-$boxZ}
                     
                     my $rij =  sqrt($xx**2 + $yy**2 + $zz**2);
                    
                     push (@rij, $rij);
                     push (@secondWater, $j);
                   
                     #print "$i $j\n";
                     
                } #for (my $j=0;$j<(scalar @CoordForAngles_noRep);$j=$j+3)
               
                                
                my @rij1 = @rij;
                
                #foreach (@rij1){print "$_ | ";};
                #print "\n";
                
                my $rij_min = min (@rij1);
                
                #print "$rij_min\n";

                my($index)= grep { $rij[$_] eq $rij_min } 0..$#rij; # find the index of $rij_min in @rij
                
                my $secondW = $secondWater[$index];
               #print "$secondW\n";
               
               
                #print "index of rij_min in water vector $index\n";
                #my $test = scalar (@Xw_cut_unique);
                #my $test1 = scalar (@rij);
                #print "size of water coord vector (O,H,H) $test\n";
                #print "size of water water dist vector (O only) $test1\n";

                     # we chose the water which is close to the protein and checked w-w distances that are smaller than cutoff
                     # now we take these two WM (OW and HW atoms) and measure the angles
                     
                     
               if ($rij_min<=$cutoffDist2){

                     $calcAngleX[0] =  $Xw_cut_unique[$i]; # O1
                     $calcAngleY[0] =  $Yw_cut_unique[$i]; 
                     $calcAngleZ[0] =  $Zw_cut_unique[$i]; 

                     $calcAngleX[1] =  $Xw_cut_unique[$i+1]; #H11
                     $calcAngleY[1] =  $Yw_cut_unique[$i+1]; 
                     $calcAngleZ[1] =  $Zw_cut_unique[$i+1];
                     
                     $calcAngleX[2] =  $Xw_cut_unique[$i+2]; #H12
                     $calcAngleY[2] =  $Yw_cut_unique[$i+2]; 
                     $calcAngleZ[2] =  $Zw_cut_unique[$i+2]; 


                     $calcAngleX[3] =  $Xw_cut_unique[$secondW]; # O2
                     $calcAngleY[3] =  $Yw_cut_unique[$secondW]; 
                     $calcAngleZ[3] =  $Zw_cut_unique[$secondW]; 

                     $calcAngleX[4] =  $Xw_cut_unique[$secondW+1]; # H21
                     $calcAngleY[4] =  $Yw_cut_unique[$secondW+1]; 
                     $calcAngleZ[4] =  $Zw_cut_unique[$secondW+1];
                     
                     $calcAngleX[5] =  $Xw_cut_unique[$secondW+2]; # H22
                     $calcAngleY[5] =  $Yw_cut_unique[$secondW+2]; 
                     $calcAngleZ[5] =  $Zw_cut_unique[$secondW+2]; 



                    ############################################################
                    #### For the H11 - 01 --- O2 ####
                    # 01 - H11 #
                    $r1_X[0] = ($calcAngleX[0]-$calcAngleX[1]);
                    $r1_Y[0] = ($calcAngleY[0]-$calcAngleY[1]);
                    $r1_Z[0] = ($calcAngleZ[0]-$calcAngleZ[1]);
                    
                                       
                   $length = sqrt($r1_X[0]**2+$r1_Y[0]**2+$r1_Z[0]**2);

                    $r1_X[0] = ($calcAngleX[0]-$calcAngleX[1])/$length; # 01 - H11 #
                    $r1_Y[0] = ($calcAngleY[0]-$calcAngleY[1])/$length;
                    $r1_Z[0] = ($calcAngleZ[0]-$calcAngleZ[1])/$length;
                    $length = 0;
                    
                    # O1 --- O2 # 
                    $r1_X[1] = ($calcAngleX[0]-$calcAngleX[3]);
                    $r1_Y[1] = ($calcAngleY[0]-$calcAngleY[3]);
                    $r1_Z[1] = ($calcAngleZ[0]-$calcAngleZ[3]);
                    
                                       
                   $length = sqrt($r1_X[1]**2+$r1_Y[1]**2+$r1_Z[1]**2);

                    $r1_X[1] = ($calcAngleX[0]-$calcAngleX[3])/$length; # O1 --- O2 #
                    $r1_Y[1] = ($calcAngleY[0]-$calcAngleY[3])/$length;
                    $r1_Z[1] = ($calcAngleZ[0]-$calcAngleZ[3])/$length;
                    $length = 0;
                    
                    $cos_four[0] = $r1_X[0]*$r1_X[1]+$r1_Y[0]*$r1_Y[1]+$r1_Z[0]*$r1_Z[1]; # vectors O1-H11 and O1-O2
                    $angles_four[0] = 180.0/pi*acos ($cos_four[0]);
                    
                    
                    #### For the H12 - 01 --- O2 ####
                    # 01 - H2 #
                    $r1_X[2] = ($calcAngleX[0]-$calcAngleX[2]);
                    $r1_Y[2] = ($calcAngleY[0]-$calcAngleY[2]);
                    $r1_Z[2] = ($calcAngleZ[0]-$calcAngleZ[2]);
                    
                                       
                   $length = sqrt($r1_X[2]**2+$r1_Y[2]**2+$r1_Z[2]**2);

                    $r1_X[2] = ($calcAngleX[0]-$calcAngleX[2])/$length; # 01 - H12 #
                    $r1_Y[2] = ($calcAngleY[0]-$calcAngleY[2])/$length;
                    $r1_Z[2] = ($calcAngleZ[0]-$calcAngleZ[2])/$length;
                    $length = 0;
                    # O1 --- O2 #
                    # the same
                                        
                    $cos_four[1] = $r1_X[2]*$r1_X[1]+$r1_Y[2]*$r1_Y[1]+$r1_Z[2]*$r1_Z[1]; # vectors O1-H12 and O1-O2
                    $angles_four[1] = 180.0/pi*acos ($cos_four[1]);



                    #### For the H21 - 02 --- O1 ####
                    # 02 - H21 #
                    $r1_X[3] = ($calcAngleX[3]-$calcAngleX[4]);
                    $r1_Y[3] = ($calcAngleY[3]-$calcAngleY[4]);
                    $r1_Z[3] = ($calcAngleZ[3]-$calcAngleZ[4]);
                    
                                       
                   $length = sqrt($r1_X[3]**2+$r1_Y[3]**2+$r1_Z[3]**2);

                    $r1_X[3] = ($calcAngleX[3]-$calcAngleX[4])/$length; # 02 - H21 #
                    $r1_Y[3] = ($calcAngleY[3]-$calcAngleY[4])/$length;
                    $r1_Z[3] = ($calcAngleZ[3]-$calcAngleZ[4])/$length;
                    $length = 0;

                    # 02 - O1 #
                    $r1_X[4] = ($calcAngleX[3]-$calcAngleX[0]);
                    $r1_Y[4] = ($calcAngleY[3]-$calcAngleY[0]);
                    $r1_Z[4] = ($calcAngleZ[3]-$calcAngleZ[0]);
                    
                                       
                   $length = sqrt($r1_X[4]**2+$r1_Y[4]**2+$r1_Z[4]**2);

                    $r1_X[4] = ($calcAngleX[3]-$calcAngleX[0])/$length; # 02 - O1 #
                    $r1_Y[4] = ($calcAngleY[3]-$calcAngleY[0])/$length;
                    $r1_Z[4] = ($calcAngleZ[3]-$calcAngleZ[0])/$length;
                    $length = 0;

                    $cos_four[2] = $r1_X[3]*$r1_X[4]+$r1_Y[3]*$r1_Y[4]+$r1_Z[3]*$r1_Z[4]; # vectors O2-H21 and O2-O1
                    $angles_four[2] = 180.0/pi*acos ($cos_four[2]);

                    #### For the H22 - 02 --- O1 ####
                    # O2 - H22 #
                    $r1_X[5] = ($calcAngleX[3]-$calcAngleX[5]);
                    $r1_Y[5] = ($calcAngleY[3]-$calcAngleY[5]);
                    $r1_Z[5] = ($calcAngleZ[3]-$calcAngleZ[5]);
                    
                                       
                   $length = sqrt($r1_X[5]**2+$r1_Y[5]**2+$r1_Z[5]**2);

                    $r1_X[5] = ($calcAngleX[3]-$calcAngleX[5])/$length; # 02 - H22 #
                    $r1_Y[5] = ($calcAngleY[3]-$calcAngleY[5])/$length;
                    $r1_Z[5] = ($calcAngleZ[3]-$calcAngleZ[5])/$length;
                    $length = 0;
                    
                    # O2 - O1 #
                    # the same
                    
                    $cos_four[3] = $r1_X[5]*$r1_X[4]+$r1_Y[5]*$r1_Y[4]+$r1_Z[5]*$r1_Z[4]; # vectors O2-H22 and O2-O1
                    $angles_four[3] = 180.0/pi*acos ($cos_four[3]);

                    ############################################################
                     
                     
                  ## checking   
                  #   my $xx = ($calcAngleX[0]-$calcAngleX[3]);
                  #   my $yy = ($calcAngleY[0]-$calcAngleY[3]);
                  #   my $zz = ($calcAngleZ[0]-$calcAngleZ[3]);
                  #   
                  #   my $r01 =  sqrt($xx**2 + $yy**2 + $zz**2);
                  #   
                  #   
                  # print "Water #$i and its closest neighbor\n";
                  # print "$r01\n";
                  # foreach (@calcAngleX) {print "X $_ "};
                  # print "\n";
                  # foreach (@calcAngleY) {print "Y $_ "};
                  # print "\n";
                  # foreach (@calcAngleZ) {print "Z $_ "};
                  # print "\n";
                
                                     
                     #print "$rij_min\n";
                     
                     #foreach (@angles_four) {print "$_ |"};
                     
                     my $angle_min = min(@angles_four);
                     #print "$angle_min | ";
                     #print "\n";
                     
                     push (@angle_min, $angle_min);
                      
                } #if ($rij_min<=$cutoffDist2)
        
         @rij = (); @rij1 = ();
         @calcAngleX = ();@calcAngleY = ();@calcAngleZ = ();
         @secondWater = ();
         @r1_X = ();@r1_Y = ();@r1_Z = ();
         @cos_four = ();

      #print "NEW Water\n";
       

      
      
       } #for (my $i=0;$i<(scalar @CoordForAngles_noRep-3);$i=$i+3)
      
         @CoordForAngles = ();
         @CoordForAngles_noRep = ();
         @CoordForAngles_noRep1 = ();
         @XYZwa =(); @Xprot = ();@Xwa = ();@Xw_cut_unique = ();
         @Yprot = ();@Ywa = ();@Yw_cut_unique = ();
         @Zprot = ();@Zwa = ();@Zw_cut_unique = ();
         $waterAtomNumb = 0;
      
      #my $test = scalar (@CoordForAngles_noRep);
      #print "Number of considered WM = $test\n";
      

      
     } #if ($currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi)
      

     

    if ($currentLine =~m/\.*t=\s*(\d+\.\d+)$/gi){
    
    # print  OUT1 "$distLow[$newFrame]\n";
      $newFrame++;
      @Xprot= 0;@Yprot= 0;@Zprot= 0; # should become empty
      @Xw= 0;@Yw= 0;@Zw= 0;@waterIndex = 0; # should become empty
      @Xw_cut= 0;@Yw_cut= 0;@Zw_cut= 0;@XYZ_cut = 0;
      
    # foreach (@Xprot) {  print "$_\n";};

     $residues = 0;$waterNumb=0;
      ##$residues1 = 0;
      $resNumber = 0;
        my ($frameNumber) = $currentLine =~m/\.*t=\s*(\d+\.\d+)$/;
        $frameNumber[$newFrame] = $frameNumber;
        if ($frameNumber >= $frameNumber_step)
         {print "Frame:".$frameNumber.". ";
          $frameNumber_step=$frameNumber_step+10000;
          my $duration_interm = time - $start;
          print "Running for: $duration_interm s\n";
         }
    
    # test speed#
    $frameNumber = $frameNumber*(-1);
    }
    
  }
 
do('/home_b/yulian/scripts/histogramPersent.pl');
histogram(1, @angle_min);
 
close IN;
close OUT1;



#foreach (@distLow) {  print "$_\n";};

#close OUT_HIST;
#close OUT_CHECK;
#do('/home_b/yulian/scripts/histogram.pl');
#histogram(0.01, @waterDistALL);

#system ('mv DDF.dat "DDF_".$system.".dat"');

my $duration = time - $start;
print "Execution time: $duration s\n";

