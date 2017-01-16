# Perl script to build the energy file neeeded by RBM10
# Uses output from VIC energy files where names of files
# are found in WEATHER_data in the local directory
#!/usr/bin/perl
use Date::Calc qw(:all);
#
# Set arrays for column/row moves using numbers 1-8
# from the direction file
@col_move=(0,0,1,1,1,0,-1,-1,-1);
@row_move=(0,1,1,0,-1,-1,-1,0,1);
#
# QUAD is an array with the obverse directions, e.g.
# 5 comes from 1, 6 comes from 2 ...
#
@quad=(0,5,6,7,8,1,2,3,4);
$pi=3.14159;
$m3_to_ft3=3.2808*3.2808*3.2808;
#
# Modification to read in argument from command line (WUR_WF_MvV_2011/04/04) 
chomp($xmask="@ARGV[0]");
$xmask=~s/^\s+//;
print "$xmask\n";
open XMASK, $xmask or die "Cannot open Direction File \n";
print "Name of Output Basin Topology File \n";
print "File name is required for defining the topology\n";
print "in the control file\n"; 
#
# Modification to read in argument from command line (WUR_WF_MvV_2011/04/04) 
chomp($node_file="@ARGV[1]");
open NODES, ">$node_file" or die "Cannot open TestFile \n";
for $i (0..366) {
  $zero[$i] = 0;
}
#
# Read the header information in the direction file
#
chomp($zz=<XMASK>);
($header,$no_cols)=split/\s+/,$zz;
#
chomp($zz=<XMASK>);
($header,$no_rows)=split/\s+/,$zz;
#
chomp($zz=<XMASK>);
($header,$xll)=split/\s+/,$zz;
chomp($zz=<XMASK>);
($header,$yll)=split/\s+/,$zz;
chomp($zz=<XMASK>);
($header,$cell_size)=split/\s+/,$zz;
print "cell size $cell_size\n";
chomp($zz=<XMASK>);
($header,$no_data)=split/\s+/,$zz;
#
# Read the flow direction file
#
for $nr (1..$no_rows) {
  $nrow=$no_rows-$nr+1;
  $zz=<XMASK>;
  @xmask=split/\s+/,$zz;
# If the first column is a blank, shift the array one place to the left
# so the columns will align properly (UW_JRY_2011/05/10)
  if (!$xmask[0]) {shift(@xmask);}
  for $nc (0..$no_cols-1) {
    if($xmask[$nc] >= 1 && $xmask[$nc] <=8) {
#
# set the direction of each cell
#
      $flow_dir=$xmask[$nc];
      $flow_dir[$nc+1][$nrow]=$flow_dir;
#
# Increment the cell number and define the row and column of
# the cell using the cell number
# 
     $no_cells++;
      $cell_row[$no_cells]=$nrow;
      $cell_col[$no_cells]=$nc+1;
#
# Define the cell number for the row/column
#
      $cell_no[$nc+1][$nrow]=$no_cells;
#
# Define the lat/long of the cell using the cell number
#
      $lat[$no_cells]=$yll+($nrow-0.5)*$cell_size;
      $lon[$no_cells]=$xll+($nc+0.5)*$cell_size;
      $nr_cell=$nrow;
      $nc_cell=$nc+1;
      $nr_ds_cell=$nr_cell+$row_move[$flow_dir];
      $nc_ds_cell=$nc_cell+$col_move[$flow_dir];
#
# Identify downstream cell as not being a headwater
#
      $head_cell[$nc_ds_cell][$nr_ds_cell]=-1;
    }
#
# Identify downstream cell as outlet 
#
    if ($xmask[$nc]==9) {
      $outlet_col=$nc+1;
      $outlet_row=$nrow;
      $lat[0]=$yll+($nrow-0.5)*$cell_size;
      $lon[0]=$xll+($nc+0.5)*$cell_size;
    }
  }
}
#
for $ncell (1..$no_cells) {
  $nc=$cell_col[$ncell];
  $nr=$cell_row[$ncell];
  $nr_ds_cell=$nr_cell+$row_move[$flow_dir];
  $nc_ds_cell=$nc_cell+$col_move[$flow_dir];
  $cell_no=$cell_no[$nc_ds_cell][$nr_ds_cell];
  $nr_us_cell[$cell_no][$quad[$flow_dir]]=$nr;
  $nc_us_cell[$cell_no][$quad[$flow_dir]]=$nc;
#
# Test for headwater cell
#
  if ($head_cell[$nc][$nr]!=-1) {
    $no_head++;
#
# Set headwaters array 
#
    $no_head[$no_head]=$no_head;
#
# Identify column/row and flow direction of headwaters
#
    $head_col[$no_head]=$nc;
    $head_row[$no_head]=$nr;
    $dir_head[$no_head]=$flow_dir[$nc][$nr];
    $head_cell[$nc][$nr]=1;
    print "$ncell $nc $nr\n";
  }
}
#
# Total number of cells and headwaters remaining to be included
#
$cell_total=$cell_no;
$no_head_rem=$no_head;
print "Number of Headwaters $no_head\n";
#
# Find the rest of the tributaries
#
$max_seg=-99999;
$stream_order=1;
for $nh (1..$no_head) {
  $stream_order[$nh]=$stream_order;
  $nc=$head_col[$nh];
  $nr=$head_row[$nh];
  $flow_dir=$dir_head[$nh];
  while ($flow_dir >= 1 && $flow_dir <= 8) {
    $no_seg++;
    &route($nc,$nr,$flow_dir);
    $cell_no=$cell_no[$nc_ds_cell][$nr_ds_cell];
    $no_head_waters[$cell_no][$quad[$flow_dir]]++;
    $trib_head_waters[$cell_no][$quad[$flow_dir]]   [$no_head_waters[$cell_no][$quad[$flow_dir]]]=$nh;
    $flow_dir=$flow_dir[$nc_ds_cell][$nr_ds_cell];
    $nc_us[$nh][$no_seg]=$nc;
    $nr_us[$nh][$no_seg]=$nr;
    $no_seg[$nh][$cell_no]=$no_seg;
    $nc=$nc_ds_cell;
    $nr=$nr_ds_cell;
  }
#
# First time through find the Main Stem
#
  if ($no_seg>$max_seg) {
    $max_seg=$no_seg;
    $main_stem=$nh;
  }
  $no_stream_seg[$nh]=$no_seg;
#
# Set the number of segments back to zero for
# the next time through
#
  $no_seg=0;
}
print  "Outlet $outlet_col $outlet_row\n"; 
print "Main Stem $main_stem\n";
#
# Main Stem has only one (1) branch
#
$no_branches[$stream_order]=1;
$branch_headwaters[$no_branches[$stream_order]][$stream_order]=$main_stem;
$nc=$head_col[$main_stem];
$nr=$head_row[$main_stem];
$conf_col[$main_stem]=$outlet_col;
$conf_row[$main_stem]=$outlet_row;
$cell_no[$outlet_col][$outlet_row]=0;
$no_seg[$nh][0]=$max_seg;
$flow_dir=$dir_head[$nh];
$nc_ds_cell=$nc+$col_move[$flow_dir];
$nr_ds_cell=$nr+$row_move[$flow_dir];
for $ns (1..$max_seg-1) {
  &route($nc_ds_cell,$nr_ds_cell,$flow_dir);
  $flow_dir=$flow_dir[$nc_ds_cell][$nr_ds_cell];
}
$nc_us_cell[$no_cell+1][$quad[$flow_dir]]=$nc_ds_cell;
$nr_us_cell[$no_cell+1][$quad[$flow_dir]]=$nr_ds_cell;
$no_head_rem--;
#
# Now find the rest
#
while ($no_head_rem > 0 && $stream_order < 10) {
  $stream_order++;
#
# Look at all the branches at stream order =$stream_order
# where $stream_order is not the conventional, rather it 
# increases with decreasing stream size.
#
  for $nb (1..$no_branches[$stream_order-1]) {
    $nh=$branch_headwaters[$nb][$stream_order-1];
#
# Starting column and row numbers
#
    $start_col=$conf_col[$nh];
    $start_row=$conf_row[$nh];
    $cell_no=$cell_no[$start_col][$start_row];
    $no_seg=$no_seg[$nh][$cell_no];
#    
#    
    for ($ns=$no_seg; $ns>1; $ns--) {
      $nc=$nc_us[$nh][$ns];
      $nr=$nr_us[$nh][$ns];
      $nhc=$head_col[$nh];
      $nhr=$head_row[$nh];
      $cell_no=$cell_no[$nc][$nr];
      $head_quad=$quad[$flow_dir[$nc_us[$nh][$ns-1]][$nr_us[$nh][$ns-1]]];
      for $nq (1..8) {
        if($nq != $head_quad) {
          $no_head_waters=$no_head_waters[$cell_no][$nq];
          if ($no_head_waters > 0) {

            $max_seg=-9999;
            for $nhw (1..$no_head_waters) {
              $nhw=$trib_head_waters[$cell_no][$nq][$nhw];
              $no_seg_test=$no_seg[$nhw][$cell_no];
              if ($no_seg_test > $max_seg) {
                $max_seg=$no_seg_test;
                $main_stem=$nhw;
              }
            }
            $no_branches[$stream_order]++;
            $no_br=$no_branches[$stream_order];
            $branch_headwaters[$no_br][$stream_order]=$main_stem;
            $trib_to[$main_stem]=$nh;
            $conf_col[$main_stem]=$nc;
            $conf_row[$main_stem]=$nr;
            $conf_cell[$main_stem]=$cell_no;
            $nch=$head_col[$main_stem];
            $nrh=$head_row[$main_stem];
            $no_head_rem--;
#            print "No Remaining $no_head_rem\n";
#           print NODES "Stream Order $stream_order Main Stem $nh\n";
#           print NODES "Res $ns $nc $nr $nq $no_branches[$stream_order] HW $main_stem $nch $nrh\n";
          }
        }
      }
    }
#  print NODES "$stream_order $nb $nh $head_col[$nh] $head_row[$nh]\n";
  }
}
print "STREAM LEVEL: $stream_order\n";
#
for ($nord=$stream_order;$nord>=1;$nord--) {
  print NODES "STREAM LEVEL $nord\n";
  for $no_br (1..$no_branches[$nord]) {
    $nh=$branch_headwaters[$no_br][$nord];
    $nc=$head_col[$nh];
    $nr=$head_row[$nh];
    $cell_no=$cell_no[$nc][$nr];
    print NODES "STREAMNAME $cell_no Headwaters $nh $trib_to[$nh]\n";
    print NODES "Node $cell_no  Row $nr Col $nc Lat $lat[$cell_no] Lon $lon[$cell_no] R.M.\n";
    $flow_dir=$dir_head[$nh];
    &route($nc,$nr,$flow_dir);
    $nc=$nc_ds_cell;
    $nr=$nr_ds_cell;
    $flow_dir=$flow_dir[$nc][$nr];
    $cell_no=$cell_no[$nc][$nr];
    print NODES "Node $cell_no  Row $nr Col $nc Lat $lat[$cell_no] Lon $lon[$cell_no] R.M.\n";
    while ($cell_no != $conf_cell[$nh]) {
      &route($nc,$nr,$flow_dir);
      $nc=$nc_ds_cell;
      $nr=$nr_ds_cell;
      $flow_dir=$flow_dir[$nc][$nr];
      $cell_no=$cell_no[$nc][$nr];
      print NODES "Node $cell_no  Row $nr Col $nc Lat $lat[$cell_no] Lon $lon[$cell_no] R.M.\n";
    }
#  $nc=$conf_col[$nh];
#  $nr=$conf_row[$nh];
  &route($nc,$nr,$flow_dir);
  $nc=$nc_ds_cell;
  $nr=$nr_ds_cell;
  $cell_no=$cell_no[$nc][$nr];
#  print NODES "Node $cell_no  Row $nr Col $nc Lat $lat[$cell_no] Lon $lon[$cell_no] R.M.\n";
  print NODES "\n";
  }
  print NODES "\n";print NODES "\n";
}
#
#Calculate     
sub route {
    ($nc,$nr,$flow_dir)=@_;
    $nc_ds_cell=$nc+$col_move[$flow_dir];
    $nr_ds_cell=$nr+$row_move[$flow_dir];
}
