#!/usr/bin/perl
use Date::Calc qw(:all);
use Math::Complex;
use Math::Trig;
$deg_radians=1/57.29577951;
#
#  Input from STDIN to identify project
#
print "Input Name of Project (e.g., River Basin)\n";
print "This script will search for a file of the form: Project.Control\n";
print "Input only the Project name.  Control is appended in the script\n";
chomp($project=@ARGV[0]);
#
#  Open project control file
#
open CONTROL, $project.'.Control' or die "Cannot open Control File \n";
chomp($project_title=<CONTROL>);
print "$project_title  \n";
#
#  Skip descriptive text
#
<CONTROL>;
<CONTROL>;
chomp($zz=<CONTROL>);
#
#  Read start and end dates and model time steps for RBM from control file.
#
(undef,$start_date)=split/:/,$zz;
print "Start Date $_ $start_date\n";
chomp($zz=<CONTROL>);
(undef,$end_date)=split/:/,$zz;
print "End date $_ $end_date\n";
chomp($zz=<CONTROL>);
(undef,$no_time_steps_per_day)=split/:/,$zz;
if ($no_time_steps_per_day>1) {
#  print "Number of times steps = $no_time_steps_per_day\n";
  die "Present version allows only a daily time step\n";}
#
#  Identify input directory for topology file
#  and output directory for RBM network file
#
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
(undef,$topo_file)=split/:/,$zz;
$topo_file=~s/^\s+//;
print "Topo File $topo_file\n";
#
#  Open topology file
#
open TOPO, $topo_file or die  "Cannot open topology file $topo_file\n";
#
# Output Network file
#
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
(undef,$network_file)=split/:/,$zz;
$network_file=~s/^\s+//;
#
#  Open RBM network file for output
#
open NETWORK, ">$network_file" or die "Cannot open network file $network_file\n";
#
# Direct Access file for flow
#
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
(undef,$flow_file)=split/:/,$zz;
$flow_file=~s/^\s+//;
#
# Direct Acces file for heat
#
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
(undef,$heat_file)=split/:/,$zz;
$heat_file=~s/^\s+//;
#
# Print header and file information to NETWORK file
#
print NETWORK "$project_title\n";
print NETWORK "$flow_file\n";
print NETWORK "$heat_file\n";
#
# Files with Mohseni parameters
#
<CONTROL>;
<CONTROL>;
chomp($zz=<CONTROL>);
(undef,$mohseni_file,$grid)=split/:/,$zz;
chomp($mohseni_file);
print "MOHSENI FILE $mohseni_file GRID:$grid\n";
#
#  If the name of the $mohseni_file contains the word 'grid', the Mohseni parameters
#  are read from an array constructed with the same format as the direction file.   This
#  is for the purpose of accommodating a gridded set of parameters, as was the case for
#  work performed at Wageningen University.  In that case, the Mohseni parameters were
#  estimated using kriging.
#
$is_gridded='F';
if (chomp($grid) == 'grid') {
  $mohseni_grid='T';
  $mohseni_alpha_file=$mohseni_file.'.alpha';
  $mohseni_alpha_file=~s/^\s+//;
  print "Mohseni alpha file $mohseni_alpha_file\n";
  &read_mohseni($mohseni_alpha_file);
  for $i (0..$ncol-1){for $j (1..$nrow) {$mohseni_alpha[$i][$j]=$mohseni_parm[$i][$j];}}
#
#@mohseni_alpha=@mohseni_parm;
  $mohseni_beta_file=$mohseni_file.'.beta';
  $mohseni_beta_file=~s/^\s+//;
  &read_mohseni($mohseni_beta_file);
  for $i (0..$ncol-1) {for $j (1..$nrow) {$mohseni_beta[$i][$j]=$mohseni_parm[$i][$j];}}
#
#@mohseni_beta=@$mohseni_parm;
  $mohseni_gamma_file=~s/^\s+//;
  $mohseni_gamma_file=$mohseni_file.'.gamma';
  &read_mohseni($mohseni_gamma_file);
  for $i (0..$ncol-1){for $j (1..$nrow) {$mohseni_gamma[$i][$j]=$mohseni_parm[$i][$j];}}
#
#@mohseni_gamma=@$mohseni_parm;
  $mohseni_mu_file=$mohseni_file.'.mu';
  $mohseni_mu_file=~s/^\s+//;
  &read_mohseni($mohseni_mu_file);
  for $i (0..$ncol-1){for $j (1..$nrow) {$mohseni_mu[$i][$j]=$mohseni_parm[$i][$j];}}
  $mohseni_timelag_file=$mohseni_file.'.timelag';
  $mohseni_timelag_file=~s/^\s+//;
  &read_mohseni($mohseni_timelag_file);
  for $i (0..$ncol-1){for $j (1..$nrow) {$mohseni_timelag[$i][$j]=$mohseni_parm[$i][$j];}}
} else {
#
# This file contains only one set of Mohseni parameters that are applied to all 
# headwater segments
#
    open MOHSENI, "$mohseni_file" or die "Cannot open Mohseni file $mohseni_file\n";
    <MOHSENI>;
    chomp ($timelag=<MOHSENI>);  
    <MOHSENI>;
    chomp ($zz=<MOHSENI>);
    ($alpha,$beta,$gamma,$mu)=split/\s+/,$zz;
    print "$alpha $beta $gamm $mu $timelag\n";
}
#
# Check if there are Heat Dump files (.true. if there are)
#
$heat_dump='F';
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
(undef,$heat_dump)=split/:/,$zz;
if ($heat_dump =~/[Tt]/) {
  $heat_dump='T';
  chomp($zz=<CONTROL>);
  chomp($zz=<CONTROL>);
  chomp($zz=<CONTROL>);
 (undef,$heat_dump_file)=split/:/,$zz;
  $heat_dump_file=~s/^\s+//;
}
#
# Read ndelta-number of computational elements in each segment
#
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
chomp($zz=<CONTROL>);
(undef,$ndelta)=split/:/,$zz;
#
#  Added file open to write files for rout program\n";
#
open FULLDATA,">FullData.Cells";
open ROUTINIT,">Rout.Cells.init";
open ROUT,">Rout.Cells";
printf NETWORK "%10d%10d%10d\n",$start_date,$end_date,$no_time_steps_per_day;
while ($zz=<TOPO>) {
  if ($zz=~/STREAMNAME/) {
     $nh++;
     ($a1,$head_node[$nh],$a2,$head_ndx,$trib_to_ndx[$nh])=split/\s+/,$zz;
     $head_seq[$head_ndx]=$nh;
#print "Stream name $stream_name\n";
     $nn=0;
  }
  if ($zz =~/Node/) {
        $nn++;
        $nseg++;
        ($a1,$node_no[$nh][$nn],$a2,$row,$a3,$col,$a4,$lat[$nh][$nn],$a5,$long[$nh][$nn])=split/\s+/,$zz;
        $mask[$nseg]=$row.'A'.$col;
        $segment_head[$nseg]=$nh;
        $junc_node[$nh]=$mask[$nseg];
        $last_seg[$nh]=$nseg;
        $cell_total[$nh]=$nn;
        $no_head=$nh;
      }
}
$nsegments=$nseg;
#
# Number of forcing cells is total number of segments
# Number of flow cells is number of cells less those that are duplicates
#
$flow_cells=$nseg-$no_head;
$force_cells=$nseg;
#
# Write the number of flow cells and forcing cells at the
# beginning of the rout model station files
#
print ROUTINIT "$flow_cells $force_cells\n";
print ROUT "$flow_cells $force_cells\n";
#print "NSEG $nseg $no_head \n";
$nseg=0;
#
# Default test for existence of a source file is FALSE
#
printf NETWORK "%10d%10d%10d%10s\n", $no_head,$flow_cells,$force_cells,$heat_dump;
if ($heat_dump=~/[Tt]/) {
  print NETWORK "$heat_dump_file\n";
}
#
#  Convert the headwaters number taken from the topology file
#  to a sequential number by matching the row/column
# 
for $nh (1..$no_head) {
  for $ns (1..$nsegments) {
    if ($junc_node[$nh] eq $mask[$ns]) {
      $trib_cell[$nh]=$ns;
      $nseg_last=$ns;
    }
  }
$nh_trib=$head_seq[$trib_to_ndx[$nh]];
#print CELLTABLE "Head $nh Trib To $trib_cell[$nh] Trib Cell $nh_trib\n";
  $lat1=$lat[$nh][1];
  $long1=$long[$nh][1];
  $total_dist=0;
  for $n (1..$cell_total[$nh]) {
    $lat2=$lat[$nh][$n];
    $long2=$long[$nh][$n];
    $dist_x=69.712*cos($lat2*$deg_radians)*($long2-$long1);
    $dist_y=69.712*($lat2-$lat1);
    $dist[$n]=sqrt($dist_x**2+$dist_y**2);
    $total_dist=$total_dist+$dist[$n];
    $lat1=$lat2;
    $long1=$long2;
  }
  $river_mile=$total_dist+$dist[2];
#  
  for $n (1..$cell_total[$nh]) {
    $nseg++;
    ($row,$col)=split/A/,$mask[$nseg];
    if ($n == 1) {
      if($is_gridded == 'T') {
        $alpha=$mohseni_alpha[$col-1][$row];
        $beta=$mohseni_beta[$col-1][$row];
        $gamma=$mohseni_gamma[$col-1][$row];
        $mu=$mohseni_mu[$col-1][$row];
        $timelag=$mohseni_timelag[$col-1][$row];
      }
      printf NETWORK "%5d Headwaters%4d TribCell%6d  Headwaters%8d         R.M. =%10.2f\n"
        ,$cell_total[$nh],$nh,$trib_cell[$nh],$nh_trib,$river_mile;
      printf NETWORK "%8.2f %8.2f %8.4f %8.2f %8.4f\n",$alpha,$beta,$gamma,$mu,$timelag;
    }  
#    $row=int($mask[$nseg]/1000);
#    print "$nseg $mask[$nseg]\n";
#    $col=$mask[$nseg]-1000*$row;
    $river_mile=$river_mile-$dist[$n+1];
    if ($n == $cell_total[$nh]) {
      $river_mile = 0.0;
    }  
    printf NETWORK "Node%6d Row%6d Column%6d  Lat%9s Long%11s R.M. =%10.2f%5d\n"
          ,$nseg,$row,$col,$lat[$nh][$n],$long[$nh][$n],$river_mile,$ndelta;
    $file_ext=$lat[$nh][$n].'_'.$long[$nh][$n];
    if ($n < $cell_total[$nh]) {
#      print "route file $rout_file\n";
      print ROUTINIT "1 $nseg $file_ext $col $row -99\n";
      print ROUTINIT "NONE\n";
      $uh_file=$file_ext.'.uh_s';
      print ROUT "1 $nseg $file_ext $col $row -99\n";
      print ROUT "$uh_file\n";
    }
    else {
      print ROUTINIT "0 $nseg $file_ext $col $row -99\n";
      print ROUT "0 $nseg $file_ext $col $row -99\n";
    }
  }
}
sub read_mohseni {
  ($file_name)=@_;
  open MOHSENI, "$file_name" or die "Cannot open Mohseni file$file_name\n";
  print "File Name $file_name\n";
  chomp($zz=<MOHSENI>);
  (undef,$ncol)=split/\s+/,$zz;
  chomp($zz=<MOHSENI>);
  (undef,$nrow)=split/\s+/,$zz;
  <MOHSENI>;
  <MOHSENI>;
  <MOHSENI>;
  <MOHSENI>;
  for $j (1..$nrow) {
    $jrow=$nrow-$j+1;
    chomp($zz=<MOHSENI>);
    @mohseni_row=split/\s+/,$zz;
    for $i (0..$ncol-1) {
      $mohseni_parm[$i][$jrow]=$mohseni_row[$i];
    }
  }
return 
}   
   
