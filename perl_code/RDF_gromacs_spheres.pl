#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use List::Util qw( min max );
use POSIX;

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
my $max_sphere = $ARGV[2]; #

#my $system = $ARGV[2]; # peptide20_WT 
#my $protein_res = 99;
open(IN,$inputFile) or die "Error opening output file. \n";

my ($groFile) = $inputFile =~ m/(.*)\.gro/g;

#open(OUT1, ">",$groFile."_all_distLower_".$cutoffDist."_nm".".dat") or die "Error opening gro file for writing: $!\n";

#open(OUT_ALL, ">",$groFile."_allDist.dat") or die "Error opening gro file for writing: $!\n";

open(OUT_HIST, ">",$groFile."_cutoff_".$cutoffDist."nm_maxSphere".$max_sphere."_rdf.dat") or die "Error opening output file for writing: $!\n";
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
my $z;my $zz;my $zzz;
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
#my $maxbin= $last_bin;
#my $box_half = $side; #*0.5;

#my $dr = $box_half/$maxbin;    # Bin width
my $dr = 0.01;
my $maxbin = $max_sphere/$dr;#5;
#my $maxbin = $cutoffDist; #/$dr;#5;


my @hist;
#my $hist = [0]*($maxbin+1);
my $bin;
my (@Xw_cut,@Yw_cut,@Zw_cut, @XYZ_cut);
my @water_sphere=0;
my $counter=0;
 

for (my $i=0;$i<$maxbin;$i++){$hist[$i] = 0;}



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

      $z = 0;

         
         
         for (my $i=0;$i<(scalar @Xprot);$i++) {  
            
          if (sqrt(($Xw[$waterNumb]-$Xprot[$i])**2+($Yw[$waterNumb]-$Yprot[$i])**2+($Zw[$waterNumb]-$Zprot[$i])**2)<=$cutoffDist){
         
         
            $waterDist[$index_line_waterDist] = sqrt(($Xw[$waterNumb]-$Xprot[$i])**2+($Yw[$waterNumb]-$Yprot[$i])**2+($Zw[$waterNumb]-$Zprot[$i])**2);

            
            
            # make new arrays with coord of OW that <=$cutoffDist
            
            
            $Xw_cut[$index_line] = $Xw[$waterNumb];
            $Yw_cut[$index_line] = $Yw[$waterNumb];
            $Zw_cut[$index_line] = $Zw[$waterNumb];
            $XYZ_cut[$index_line] = $XYZ[$waterNumb];

            #print "$Xw_cut[$index_line]  $Xw[$waterNumb]\n";
            $index_line_waterDist++;
            $index_line++;
            
            
          }
          #else{$waterDist[0] = 0;}
         
         } #for (my $i=0;$i<(scalar @Xprot);$i++)
         

            $waterNumb++;
            
            
            my $test = scalar(@waterDist);

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
         
          
          
          $waterNumb=0;
         $index_line=0;
      
          #foreach (@Xw_cut) {  print "$_ ";};
          
           my %seenXYZ;
           $seenXYZ{$_}++ for @XYZ_cut;
           my @XYZ_cut_unique = keys %seenXYZ;
           %seenXYZ = (0, 0);        
                    

           
           for (my $i=0;$i<(scalar @XYZ_cut_unique);$i++){
            
           ($Xw_cut_unique[$i],$Yw_cut_unique[$i],$Zw_cut_unique[$i]) = split  /\s+/, $XYZ_cut_unique[$i];
           }
           
           $waterNumbFixed = scalar @XYZ_cut_unique;
           
           #$waterNumbFixed_sum = $waterNumbFixed_sum+$waterNumbFixed;
           #print "Water Cutoff $waterNumbFixed $waterNumbFixed_sum\n";

#my @water_sphere=0;
#my $counter=0;
                 for (my $i=0;$i<(scalar @Xw_cut_unique-1);$i++) {#$i<(scalar @Xw_cut_unique)
                 
                 
                 
                  for (my $j=$i+1;$j<(scalar @Xw_cut_unique);$j++) {#$j<(scalar @Xw_cut_unique )
                  #for (my $j=$i+1;$j<(scalar @Xw);$j++) {
                 
                     
                     #print "$i $j  \n";
                     
                     # calculate all OW-OW distances
                     
                     my $rij =  sqrt(($Xw_cut_unique[$i]-$Xw_cut_unique[$j])**2+($Yw_cut_unique[$i]-$Yw_cut_unique[$j])**2+($Zw_cut_unique[$i]-$Zw_cut_unique[$j])**2);
                     
                                         
                     #my $test = $rij;
                     #my $test1 = $rij/$dr;

                     
                     $bin = int(ceil($rij/$dr)); # determine in which bin the distance falls
                     #$bin = int(floor($rij/$dr)); # determine in which bin the distance falls

                     #print "$test $test1 $bin\n";
                     
                     if ($bin <= $maxbin){$hist[$bin] += 1;}

                       if ($bin == $maxbin){
                        $water_sphere[$counter] = $i;
                        $counter++;
                        }
                                         
                     
                  }
                 
                 }
  
  
  

    #last;
           my %seen_water_sphere;
           $seen_water_sphere{$_}++ for @water_sphere;
           my @water_sphere_unique = keys %seen_water_sphere;
           %seen_water_sphere = (0, 0);    
           my $water_amount = scalar (@water_sphere_unique);
           my $sphereV1 =  pi*($max_sphere)**3*4/3;
           my $IdealDensity2 =  $water_amount/$sphereV1;


    @water_sphere=();
    $counter=0;
    
    print "$IdealDensity2\n";
    
    @XYZ_cut_unique=();
    @Xw_cut_unique=();
    @Yw_cut_unique=();
    @Zw_cut_unique=();

    } # if ($currentLine =~m/^\s*\-?\d+\.\d+\s+\-?\d+\.\d+\s+\-?\d+\.\d+\s*/gi)
    
  
    if ($currentLine =~m/\.*t=\s*(\d+\.\d+)$/gi){
    
    # print  OUT1 "$distLow[$newFrame]\n";
      $newFrame++;
      @Xprot= 0;@Yprot= 0;@Zprot= 0; # should become empty
      @Xw= 0;@Yw= 0;@Zw= 0;@waterIndex = 0; # should become empty
      @Xw_cut= 0;@Yw_cut= 0;@Zw_cut= 0;@XYZ_cut = 0;
      
    # foreach (@Xprot) {  print "$_\n";};

      $z = 0;
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
 
 
 # Normalization
my @val_rdf; 
my $nframes = $newFrame+1;
my $naparticles = $waterNumbFixed_sum/$nframes;


#print "$naparticles $nframes\n";


#my $rho = $naparticles/(($max_sphere)**3); 
#my $norm = 2.0*pi*$dr*$rho*$nframes*$naparticles;
 
my $sphereV =  pi*($max_sphere)**3*4/3;

my $IdealDensity =  $naparticles/$sphereV;# incorrect!!! :

# 1. $naparticles based on $cutoffDist (1nm) and $sphereV based on $max_sphere (5 nm)
# 2. Ideal density = number of water paricles in the ideal sphere WITHOUT a protein !!!

#my $IdealDensity = 13.5428124175756;
#my $IdealDensity = 13.105992;

#my $IdealDensity = 28;

#print "sphereV =  $sphereV\n";
print "IdealDensity =  $IdealDensity\n";
#print "IdealDensity =  $IdealDensity\n";


#my $ntotal = $sphereV*$IdealDensity;
 
#rint "$dr $rho $nframes $naparticles \n"; 
 
#print "$nframes $naparticles $rho  $norm \n"; 
 
#print "$hist[10] $naparticles $nframes $sphereV\n";
 
my @rrr; 
my $rrr; 
my $rrr_next;
my $bin_sum = 0;
open(OUT_CHECK, ">", $groFile."_cutoff_".$cutoffDist."nm_maxSphere".$max_sphere."_rdf_check.dat") or die "Error opening gro file for writing: $!\n";



 for (my $i=0;$i<$maxbin;$i++){
   
   #if ($i ==0){
   #   $rrr[$i] = ($i)*$dr;
   #   $rrr = ($i)*$dr;
   #   $rrr_next = ($i+0.5)*$dr;

   #}
   #else{
      #$rrr[$i] = ($i+1-0.5)*$dr;
      #$rrr = ($i+1-0.5)*$dr;
      #$rrr_next = ($i+1+0.5)*$dr;
      
      $rrr[$i] = ($i+0.5)*$dr;
      $rrr = ($i+0.5)*$dr;
      $rrr_next = ($i+1.5)*$dr;  
      

   #}
   



   my $V_bin = pi*$rrr_next**3*4/3 - pi*$rrr**3*4/3;
   $bin_sum = $bin_sum+$V_bin;
   
   
   my $nideal = $V_bin*$IdealDensity;  # $nideal=$V_bin*$naparticles/$sphereV
   
   #$val_rdf[$i] = $hist[$i]/$nideal;
   $val_rdf[$i] = $hist[$i]/($nframes*$nideal);#*$nideal
   $val_rdf[$i] = 2*$val_rdf[$i]/($naparticles-1);

   #print "$V_bin  $rrr $nideal $val_rdf[$i]\n";
   
   my $currentNewLine = sprintf ("%.5f    %.5f    %.5f    %.5f",$V_bin, $rrr, $nideal, $val_rdf[$i]);
   print OUT_CHECK "$currentNewLine\n";
   



   
   #print "$nideal | $val_rdf[$i] \n";
   
   #$val_rdf[$i] = $hist[$i]/($nframes*$naparticles*$nideal);

   #print "$hist[$i] $norm $rrr $dr\n"; 
   #$val_rdf[$i] = $hist[$i]/$norm/($rrr**2+$dr**2/12.0);
   #$val_rdf[$i] = $hist[$i]/($naparticles*$nframes*$sphereV);
   #$val_rdf[$i] = $hist[$i]/($nframes*$IdealDensity*$naparticles);

   
 }
 
 print OUT_CHECK "sphereV_init  =  $sphereV\n";
 print OUT_CHECK "sphereV_rings =  $bin_sum\n";
 
 
 for (my $i=0;$i<(scalar @val_rdf);$i++){
   
   print OUT_HIST  "$rrr[$i] $val_rdf[$i]\n";
   
 }
 
 
 
 #foreach (@val_rdf) {print OUT_HIST "$_\n"};
 
  
close IN;



#foreach (@distLow) {  print "$_\n";};

close OUT_HIST;
close OUT_CHECK;
#do('/home_b/yulian/scripts/histogram.pl');
#histogram(0.01, @waterDistALL);

#system ('mv DDF.dat "DDF_".$system.".dat"');

my $duration = time - $start;
print "Execution time: $duration s\n";

