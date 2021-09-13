#!/usr/bin/perl
use strict;
use warnings;
use Math::Trig;
use List::Util qw(sum);
use List::Util qw( min max );
use POSIX qw(ceil floor);

my $start = time;

my ($inputFile_allelec) = $ARGV[0];
my ($inputFile_dat) = $ARGV[1];
my ($Nframes) = $ARGV[2];


#my ($NumberContacts) = $ARGV[2]; 

if ($inputFile_allelec =~ /.bz2/){
    `bunzip2 $inputFile_allelec`;
    $inputFile_allelec =~ s/.bz2//g;
}

my $last;
my $frame = 0;
my $cont_line = 0;
my $cont_frame = 0;
my $cont_ndx = 0;
my $cont_ndx1 = 0;
my $elec_total_numb;
my $elec_total_numb_max;
my $beads_total_numb;
my @frame_contacts_ener;
my (@contBeads_dat, @contBeads1, @contBeads2);
my @all_frames_contacts_ener;
my @all_frames_contacts_ener_avg;
my @all_frames_contacts_ener_std;
my $all_frames_contacts_ener_avg;
my $all_frames_contacts_ener_std;
my @array1;
my @elec_matrix;
my @elec_matrix_std;
my $contBeads1;
my $contBeads2;
my $frame_contacts_ener;
my@array2avg;
my@array2std;
my@array3std;

    
open(IN_DAT,$inputFile_dat) or die "Error opening input dat file. \n";
open(IN_ELEC,$inputFile_allelec) or die "Error opening input allCont file. \n";

my ($FilePrefix) = $inputFile_allelec =~ m/(.*)\.dat/g;



open(OUT1, ">",$FilePrefix."_ElecperContAvg".".dat") or die "Error opening pdb file for writing: $!\n";
open(OUT2, ">",$FilePrefix."_ElecperContStd".".dat") or die "Error opening pdb file for writing: $!\n";


while (my $currentLine = <IN_DAT>){

    if ($currentLine=~m/^\s+\d+\s+\d+\s+\-?\d\.\d00$/gi){#    5  303  307    39.915 1.000000 Contacts. to add old contacts in the end
        
       my  ($i) = $currentLine=~m/^\s+\d+\s+(\d+)\s+\-?\d\.\d00$/;
       
       #$contacts_epsilon[$cont_ndx]=$arg1;
       $contBeads_dat[$cont_ndx]=$i;
       #print "$i ";
       #print OUT "$contacts_epsilon[$cont_ndx] ";
       $cont_ndx++;
        
    }
    
    if($currentLine=~m/electrostatic residues/gi){
      

      ($elec_total_numb) = $currentLine =~m/\s+(\d+)/;
      $elec_total_numb_max = ($elec_total_numb*($elec_total_numb-1))/2
      #print "$elec_total_numb\n";
      #next;
    }

    if($currentLine=~m/Hookean bonded pairs/gi){
      

      ($beads_total_numb) = $currentLine =~m/\s+(\d+)/;
      $beads_total_numb = $beads_total_numb+1;
      #print "$beads_total_numb\n";
      #next;
    }


}
close IN_DAT;

#print "\n";
#foreach (@contBeads_dat){print "$_ "; }


   for (my $i=0;$i<($beads_total_numb);$i++){
     for (my $j=0;$j<($beads_total_numb);$j++){

     $elec_matrix[$j][$i][0] = 0;
     
     }
   }


while (my $currentLine = <IN_ELEC>){
 
$last = $currentLine if eof;


   #if ($currentLine=~m/\s+(frame)/gi){
   #   my ($ifr) = $currentLine =~m/\s*(frame)/;
   #   $cont_frame++;
   #   #print "$ifr $cont_frame\n";
   #}
   
   if ($currentLine=~m/\s+(frame)/gi and $cont_line>0){
     #($frame) = $currentLine =~m/\s+(frame)/;
     $frame=1;
     $cont_frame++;
     #print "$frame\n";
     #print "$cont_line\n";
   }
   

 
  if ($currentLine=~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/gi){

    #($contBeads1[$cont_ndx1], $contBeads2[$cont_ndx1], $frame_contacts_ener[$cont_ndx1]) = $currentLine =~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/;
    ($contBeads1, $contBeads2, $frame_contacts_ener) = $currentLine =~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/;

       
     #process data within one frame
     
     $elec_matrix[$contBeads2-1][$contBeads1-1][0] += $frame_contacts_ener; # [$contBeads2-1] # in dat file vs # in the vector which starts from zero, not one
     
  
    if(($frame==1 and $cont_line>0) or $last){
     
     if ($last){$cont_frame++;}


     # report in the end of frame:
     
     # print $elec_matrix[$cont_frame-1][$contBeads1][$contBeads2]

     
     #print "line # $cont_line\n";
     #print "frame# $cont_frame\n";
           

     #my $currentNewLine = sprintf ("%6.3f %6.3f",$sum_NNcont_IDR_IDR_frame, $sum_NNcont_IDR_FOLDED_frame);
     #print OUT "$currentNewLine ";  
     #print OUT "\n";
    
     # end of the frame:
     $cont_ndx1=-1;
     @contBeads1=0;
     @contBeads2=0;
     @frame_contacts_ener=0;
     #$sum_NNcont_IDR_IDR_frame=0;
     #$sum_NNcont_IDR_FOLDED_frame=0;
     #print OUT "\n";
     $frame=0;
    } # if the end of the frame

    $cont_ndx1++;
    
  }  #if ($currentLine=~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/gi)
    
    $cont_line++;
   
}# end while (my $currentLine = <IN_ELEC>)



#print "$elec_matrix[41][30][1]\n";




for (my $i=0;$i<($beads_total_numb);$i++){
  for (my $j=0;$j<($beads_total_numb);$j++){
  
    $array2avg[$i][$j] = $elec_matrix[$i][$j][0]/$Nframes;
   
  } #for (my $j=0...
  
} #for (my $i=0...


foreach my $row(@array2avg){
   foreach my $val(@$row){
    
    #print OUT1 "$val ";
    my $currentNewLine1 = sprintf ("%9.6f",$val); # left justified in a field of 3 characters
    print OUT1 "$currentNewLine1 ";
    
   }
   print OUT1 "\n";
}



###FOR STDV##
close IN_ELEC;
open(IN_ELEC,$inputFile_allelec) or die "Error opening input allCont file. \n";
$cont_ndx1=0;
$cont_line=0;


   for (my $i=0;$i<($beads_total_numb);$i++){
     for (my $j=0;$j<($beads_total_numb);$j++){

     $elec_matrix_std[$j][$i] = 0;
     #$array3std[$i][$j] = 0;
     }
   }
   
   
while (my $currentLine = <IN_ELEC>){
 
$last = $currentLine if eof;


   #if ($currentLine=~m/\s+(frame)/gi){
   #   my ($ifr) = $currentLine =~m/\s*(frame)/;
   #   $cont_frame++;
   #   #print "$ifr $cont_frame\n";
   #}
   
   if ($currentLine=~m/\s+(frame)/gi and $cont_line>0){
     #($frame) = $currentLine =~m/\s+(frame)/;
     $frame=1;
     $cont_frame++;
     #print "$frame\n";
     #print "$cont_line\n";
   }
   
 
  if ($currentLine=~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/gi){

    #($contBeads1[$cont_ndx1], $contBeads2[$cont_ndx1], $frame_contacts_ener[$cont_ndx1]) = $currentLine =~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/;
    ($contBeads1, $contBeads2, $frame_contacts_ener) = $currentLine =~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/;

       
     #process data within one frame
     
     $elec_matrix_std[$contBeads2-1][$contBeads1-1] += ($array2avg[$contBeads2-1][$contBeads1-1]-$frame_contacts_ener)**2; # [$contBeads2-1] # in dat file vs # in the vector which starts from zero, not one
     
  
    if(($frame==1 and $cont_line>0) or $last){
     
     if ($last){$cont_frame++;}

     # end of the frame:
     

     
     
     
     $cont_ndx1=-1;
     @contBeads1=0;
     @contBeads2=0;
     @frame_contacts_ener=0;
     #$sum_NNcont_IDR_IDR_frame=0;
     #$sum_NNcont_IDR_FOLDED_frame=0;
     #print OUT "\n";
     $frame=0;
    } # if the end of the frame

    $cont_ndx1++;
    
  }  #if ($currentLine=~m/\s+(\d+)\s+(\d+)\s+(\-?\d+.\d+E?\-?\d*)/gi)
    
    $cont_line++;
   
}# end while (my $currentLine = <IN_ELEC>)



#foreach my $row(@elec_matrix_std){
#   foreach my $val(@$row){
#      print  "$val ";
#   }
#   print "\n";
#}



     #for (my $i=0;$i<($beads_total_numb);$i++){
     #  for (my $j=0;$j<($beads_total_numb);$j++){
     #
     #   $array2std[$i][$j] += ($array2avg[$i][$j]-$elec_matrix_std[$i][$j])**2;
     #
     #  } #for (my $j=0...
     #
     #} #for (my $i=0..

for (my $i=0;$i<($beads_total_numb);$i++){
  for (my $j=0;$j<($beads_total_numb);$j++){
  
    $array3std[$i][$j] = ($elec_matrix_std[$i][$j]/($Nframes-1))**0.5;
   
  } #for (my $j=0...
  
} #for (my $i=0...


foreach my $row(@array3std){
   foreach my $val(@$row){
    
    #print OUT1 "$val ";
    my $currentNewLine1 = sprintf ("%9.6f",$val); # left justified in a field of 3 characters
    print OUT2 "$currentNewLine1 ";

   }
   print OUT2 "\n";
}

######FOR STDV########






if ($inputFile_allelec =~ /.dat/){
    `bzip2 $inputFile_allelec`;
    #$inputFile_dat =~ s///g;
}


close OUT1;
close OUT2;
close IN_ELEC;

#sub mean {
#    return sum(@_)/@_;
#    }
#
#sub stdv {
#    
#    my $average = &mean(@_);
#    my $sqtotal = 0;
#    foreach(@_) {
#                $sqtotal += ($average-$_) ** 2;
#    }
#    my $std = ($sqtotal / (@_-1)) ** 0.5;
#        return $std;
#    }



my $duration = time - $start;
my $duration_min = $duration/60;
print "Execution time: $duration_min min ($duration s)\n";




