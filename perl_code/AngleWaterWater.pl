#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use List::Util qw( min max );
use POSIX;
use Data::Dumper;

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
#my $last_bin = $ARGV[2]; #100?
#my $side = $ARGV[3]; #box side 10 nm?
#my $max_sphere = $ARGV[2]; #

#my $system = $ARGV[2]; # peptide20_WT 
#my $protein_res = 99;
open(IN,$inputFile) or die "Error opening output file. \n";

my ($groFile) = $inputFile =~ m/(.*)\.gro/g;

open(OUT, ">",$groFile."_cutoff_".$cutoffDist."_nm_Qtet".".dat") or die "Error opening gro file for writing: $!\n";
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
my (@calcTetrahX,@calcTetrahY,@calcTetrahZ);
my (@r1_X,@r1_Y,@r1_Z);
my $cos_sum=0;
my @Qtet;
my @keys;
my @waterDistALL_round;
my @angles;
#for (my $i=0;$i<$maxbin;$i++){$hist[$i] = 0;}



while (my$currentLine = <IN>)
  {
  $index_line1++;
  
    if ($currentLine =~m/^\s+\d+(?!SOL)\w{3}\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi or
        + $currentLine =~m/^\s+\d+N\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi or
        + $currentLine =~m/^\s+\d+C\w+\s+(?!H)\w+\s+\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi ){


     ($resNumber, $Xprot[$residues],$Yprot[$residues],$Zprot[$residues]) = $currentLine =~m/^\s+(\d+)(?!SOL)\w{3}\w*\s+\w+\s+\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
     
      $residues++;
       
    }
    
      if ($currentLine =~m/^\s*\d+SOL\s+OW\d*\s*\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi){
         
         
       ($waterIndex[$waterNumb], $Xw[$waterNumb],$Yw[$waterNumb],$Zw[$waterNumb]) = $currentLine=~m/^\s*(\d+)SOL\s+OW\d*\s*\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)$/;
       ($XYZ[$waterNumb]) = $currentLine=~m/^\s*\d+SOL\s+OW\d*\s*\d+\s+(\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+)$/;

               
         for (my $i=0;$i<(scalar @Xprot);$i++) {  
            
          if (sqrt(($Xw[$waterNumb]-$Xprot[$i])**2+($Yw[$waterNumb]-$Yprot[$i])**2+($Zw[$waterNumb]-$Zprot[$i])**2)<=$cutoffDist){
         
         
            $waterDist[$index_line_waterDist] = sqrt(($Xw[$waterNumb]-$Xprot[$i])**2+($Yw[$waterNumb]-$Yprot[$i])**2+($Zw[$waterNumb]-$Zprot[$i])**2);

            
            
            # make new arrays with coord of OW that <=$cutoffDist
            
            
            #$Xw_cut[$index_line] = $Xw[$waterNumb];
            #$Yw_cut[$index_line] = $Yw[$waterNumb];
            #$Zw_cut[$index_line] = $Zw[$waterNumb];
            $XYZ_cut[$index_line] = $XYZ[$waterNumb];

            #print "$Xw_cut[$index_line]  $Xw[$waterNumb]\n";
            $index_line_waterDist++;
            $index_line++;
            
            
          }
          #else{$waterDist[0] = 0;}
         
         } #for (my $i=0;$i<(scalar @Xprot);$i++)
         

            $waterNumb++;
            
            
            #my $test = scalar(@waterDist);
            
 # ! Important.We need to define the min dist between protein and OW of water molucule. 

            if ( scalar(@waterDist) > 0){
               
            my $min = min @waterDist;
            #print OUT_ALL "$min \n";
            
            $waterDistALL[$index_waterDistALL] = $min;
            $index_waterDistALL++;
            
            @waterDist=();
            $index_line_waterDist=0;
            }
            
            
            
      } #if ($currentLine =~m/^\s*\d+SOL\s+OW\d*\s*\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+$/gi)
      
      
    if (not $distLow[$newFrame]) {$distLow[$newFrame] = 0};
 
    if ($currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi){ # end of a frame
         
         ($boxX, $boxY, $boxZ) = $currentLine=~m/^\s*\-?(\d+\.\d+)\s+\-?(\d+\.\d+)\s+\-?(\d+\.\d+)\s*/;
          
          #print "$boxX, $boxY, $boxZ\n";
          
          $waterNumb=0;
         $index_line=0;
      
          #foreach (@Xw_cut) {  print "$_ ";};
          
 # ! Important.We needed to define the min dist between protein and OW of water molucule.
 # Now we chose the unique coordinates from $XYZ_cut in order to count chosen WM only once.
 # (coordinates were repeted for each case whe WM was close to any heavy atom of the protein)
 
           my %seenXYZ;
           $seenXYZ{$_}++ for @XYZ_cut;
           my @XYZ_cut_unique = keys %seenXYZ;
           %seenXYZ = (0, 0);        
                    

           
           for (my $i=0;$i<(scalar @XYZ_cut_unique);$i++){
            
           ($Xw_cut_unique[$i],$Yw_cut_unique[$i],$Zw_cut_unique[$i]) = split  /\s+/, $XYZ_cut_unique[$i];
           }
           
           $waterNumbFixed = scalar @XYZ_cut_unique;
           
           $waterNumbFixed_sum = $waterNumbFixed_sum+$waterNumbFixed;
           #print "Water Cutoff $waterNumbFixed $waterNumbFixed_sum\n";


# Here we can connect @Xw_cut_unique with @Xprot, etc. in order to include heavy atoms of the protein to
# the calculation of tetrahidrality of water


                 for (my $i=0;$i<(scalar @Xw_cut_unique);$i++) {#$i<(scalar @Xw_cut_unique-1)
                 
                  for (my $j=0;$j<(scalar @Xw_cut_unique);$j++) {#$j<(scalar @Xw_cut_unique )

                    #next if $i==$j;
                     # calculate all OW-OW distances
                     # take into account PBC

                     $xx = ($Xw_cut_unique[$i]-$Xw_cut_unique[$j]);
                     $yy = ($Yw_cut_unique[$i]-$Yw_cut_unique[$j]);
                     $zz = ($Zw_cut_unique[$i]-$Zw_cut_unique[$j]);

                     if ($xx <-$boxX/2) {$xx=$xx+$boxX}
                     if ($xx >$boxX/2) {$xx=$xx-$boxX}
                     
                     if ($yy <-$boxY/2) {$yy=$yy+$boxY}
                     if ($yy >$boxY/2) {$yy=$yy-$boxY}
                     
                     if ($zz <-$boxZ/2) {$zz=$zz+$boxZ}
                     if ($zz >$boxZ/2) {$zz=$zz-$boxZ}
                     
                     $rij =  sqrt($xx**2 + $yy**2 + $zz**2);
                     
                     push (@rij, $rij);
                     #my $test = scalar (@rij);
                     #print "$test |"
                     
                     
                   
                     
                  } #for (my $j=$i+1;$j<(scalar @Xw_cut_unique);$j++)
                   
                  #print "$rij_min\n";
                   
                  # my $test = scalar (@Xw_cut_unique);
                  # my $test1 = scalar (@rij);
                  #print "$test $test1\n";
                    
                  # coordinates of the $i water molecule (OW atom)
                   $calcTetrahX[0] =  $Xw_cut_unique[$i];
                   $calcTetrahY[0] =  $Yw_cut_unique[$i];
                   $calcTetrahZ[0] =  $Zw_cut_unique[$i];

                   
                  # Find four closest neighbors: 
                  my @rij1 = @rij;
                  
                  for (my $k=1-1;$k<5;$k++){
                    
                     my $rij_min = min (@rij1);
                     
                      if ($rij_min==0){
                        @rij1 = grep { $_ != $rij_min } @rij1; # remove $rij_min from @rij
                        next;
                      };

                     #print "$rij_min |";
                     #print "$i |";

                     my($index)= grep { $rij[$_] eq $rij_min } 0..$#rij; # find the index of $rij_min in @rij
                     $calcTetrahX[$k] =  $Xw_cut_unique[$index]; # assuming that @rij and @Xw_cut_unique have the same indexes +$i+1
                     $calcTetrahY[$k] =  $Yw_cut_unique[$index]; # assuming that @rij and @Yw_cut_unique have the same indexes +$i+1
                     $calcTetrahZ[$k] =  $Zw_cut_unique[$index]; # assuming that @rij and @Zw_cut_unique have the same indexes +$i+1

                     #print "$index |";
                     
                     
                     @rij1 = grep { $_ != $rij_min } @rij1; # remove $rij_min from @rij
                  }
                  
                   #print "Water #$i and its four closest neighburs\n";
                   #foreach (@calcTetrahX) {print "X $_ "};
                   #print "\n";
                   #foreach (@calcTetrahY) {print "Y $_ "};
                   #print "\n";
                   #foreach (@calcTetrahZ) {print "Z $_ "};
                   #print "\n";


                   # The unit vectors in the directions of the four bonds:
                   for (my $i=0;$i<4;$i++){
                    
                    $r1_X[$i] = ($calcTetrahX[$i+1]-$calcTetrahX[0]);
                    $r1_Y[$i] = ($calcTetrahY[$i+1]-$calcTetrahY[0]);
                    $r1_Z[$i] = ($calcTetrahZ[$i+1]-$calcTetrahZ[0]);
                    
                    #print "$r1_X[$i] $r1_Y[$i] $r1_Z[$i]\n";
                    
                   my $length = sqrt($r1_X[$i]**2+$r1_Y[$i]**2+$r1_Z[$i]**2);

                    $r1_X[$i] = ($calcTetrahX[$i+1]-$calcTetrahX[0])/$length;
                    $r1_Y[$i] = ($calcTetrahY[$i+1]-$calcTetrahY[0])/$length;
                    $r1_Z[$i] = ($calcTetrahZ[$i+1]-$calcTetrahZ[0])/$length;
                    
                    #print "$r1_X[$i] $r1_Y[$i] $r1_Z[$i]\n";

                   }
                   
                  # calculate Sg (Nutt and Smith JACS 2008, 130, 13066-13073):
                  
                 for (my $j=0;$j<3;$j++){
                    
                    for (my $k=$j+1;$k<4;$k++){
                   
                       #my $length_j = sqrt($r1_X[$j]**2+$r1_Y[$j]**2+$r1_Z[$j]**2);
                       #my $length_k = sqrt($r1_X[$k]**2+$r1_Y[$k]**2+$r1_Z[$k]**2);
                       #print "$length_j $length_k\n";
                       
                       #my $cos_jk = ($r1_X[$j]*$r1_X[$k]+$r1_Y[$j]*$r1_Y[$k]+$r1_Z[$j]*$r1_Z[$k])/($length_j*$length_k);
                       my $cos_jk = ($r1_X[$j]*$r1_X[$k]+$r1_Y[$j]*$r1_Y[$k]+$r1_Z[$j]*$r1_Z[$k]);

                          
                          
                          my $angle = acos ($cos_jk);
                          $angle = $angle*180/3.14;
                          print "$angle\n";
                          push (@angles, $angle);
                          
                          $cos_jk = ($cos_jk + 1/3)**2;
                          
                          #print "$cos_jk\n";
                          
                          $cos_sum = $cos_sum + $cos_jk;
                       
                       #print "$cos_sum\n";
                    }
                 }
                   
                   $Qtet[$i] = 1-3/32*$cos_sum;
                   
                   #if ($Qtet[$i]!=1){print " $Qtet[$i]\n";}
                   
                   
                   @calcTetrahX=();@calcTetrahY=();@calcTetrahZ=();
                   @rij = (); @rij1 = ();$cos_sum = 0;
                   

                  
                  
                 }#for (my $i=0;$i<(scalar @Xw_cut_unique-1);$i++)
                 

               #foreach (@waterDistALL){print "$_\n";};
               #my $waterDistALL_size = scalar (@waterDistALL);
               #print "$waterDistALL_size\n";
               
               for (my $k=0;$k<scalar(@waterDistALL);$k++) {
               
                 my $NewLine = sprintf ("%.3f    %.5f    ",$waterDistALL[$k],$Qtet[$k]);
               
                #print OUT "$waterDistALL[$k]  $Qtet[$k]\n";
                 print OUT "$NewLine\n";
               
               }
               
               #for (my $k=0;$k<scalar(@waterDistALL);$k++) {
               #
               #my $NewCol1 = sprintf ("%.2f    ",$waterDistALL[$k]);
               #
               #print OUT "$NewCol1";
               #
               #} 
               #print OUT "\n";
               #
               #for (my $k=0;$k<scalar(@Qtet);$k++) {
               #
               #my $NewCol2 = sprintf ("%.5f    ",$Qtet[$k]);
               #
               #print OUT "$NewCol2";
               #
               #}  
               #print OUT "\n";
               
               
               
      #########################
               # Make output Qtet as a function of peptide-water distance
               
               #1 make a hash that collects values of $Qtet according to the distance;
               # hash keys: from 0.01 to $cutoffDist with step 0.01 
               
               #my @keys   = qw(root rhino root root root root root root root root root root domainte root);
               #my @values = qw(stam rhino jam onetwo domante ftpsi jay testwp contra raul vnod foos raul bruce);

               #@waterDistALL; # keys
               #@Qtet; # values
               
               

               #for (my $k=0;$k<scalar(@waterDistALL);$k++) {
               #$waterDistALL_round[$k] = sprintf "%.2f", $waterDistALL[$k];
               #}
               # see also round function in histogram.pl in ~yulian/scripts/     
                    
                              
               #foreach (@waterDistALL_round) {print "$_\n";}
               
               #my %hash;
               #@hash{@array1} = @array2;
               #@hash{@waterDistALL_round} = @Qtet;
               #print Dumper \%hash;

               #$keys[0] = 0;
               #for (my$i=1;$i <= $cutoffDist/0.01;$i++){$keys[$i] = $keys[$i-1]+0.01;}
               ##foreach (@keys) {print "$_\n";}
               #my %hash;
               #for my $idx (0 .. $#keys) {push @{ $hash{ $keys[$idx] } }, $Qtet[$idx];}
               ##print Dumper \%hash;
               
              # @Qtet - hash values
               
               #2 Each new frame adds a new value $Qtet to the previous one with the same distance
               
               #3 Devide each hash value by the number of added $Qtet in order to get average values
               
      #########################
        
               
               @waterDistALL =();
               $index_waterDistALL =0;
    #last;
    @XYZ_cut_unique=();
    @Xw_cut_unique=();
    @Yw_cut_unique=();
    @Zw_cut_unique=();
    @Qtet=();
    
    } # if ($currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi)
    
  
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
histogram(1, @angles);

 
close IN;
close OUT;
close OUT1;


#foreach (@distLow) {  print "$_\n";};

#close OUT_HIST;
#close OUT_CHECK;
#do('/home_b/yulian/scripts/histogram.pl');
#histogram(0.01, @waterDistALL);

#system ('mv DDF.dat "DDF_".$system.".dat"');

my $duration = time - $start;
print "Execution time: $duration s\n";

