#!C:\perl\bin\perl.exe
use Date::Calc qw(:all);
$out_dir='./Plot';
$no_head=0;
$no_nodes=0;
print "Number of periods per day-nwpd\n";
chomp($nwpd=<STDIN>);
print "Plot control file (.Control appended)\n";
chomp($plot_file=<STDIN>);
$plot_file=$plot_file.'.Control';
open PLOT_CONTROL, "$plot_file";
<PLOT_CONTROL>;
chomp($plot_dir<PLOT_CONTROL>);
print "$plot_dir\n";
chomp($date=<PLOT_CONTROL>);
($start_year,$start_mnth,$start_day)=split/\s+/,$date;
print "$start_year $start_mnth $start_day\n";
$start_date=Date_to_Days($start_year,$start_mnth,$start_day);
chomp($date=<PLOT_CONTROL>);
($stop_year,$stop_mnth,$stop_day)=split/\s+/,$date;
$stop_date=Date_to_Days($stop_year,$stop_mnth,$stop_day);
print "$stop_year $stop_mnth $stop_day\n";
$no_days=$nwpd*($stop_date-$start_date+1);
print "Number of days $no_days\n";
print "1\n";
chomp($in_file=<PLOT_CONTROL>);
print "Temperature files - $in_file\n";
print "2\n";
chomp($nseg=<PLOT_CONTROL>);
print " Segments $nseg\n";
$seg_list="Salmon_0.5.Spat";
while ($station=<PLOT_CONTROL>) {
  ($plot_row,$plot_col,$out_file)=split/\s+/,$station;
open LIST, "$seg_list";
$seg_count=0;
$no_match=0;
while ($no_match == 0) {
  $xx=<LIST>;
# print "LIST $xx\n";
  (undef,$head,$seg,$row,$col,$lat,$long)=split/\s+/,$xx;
   print "$row $col $seg $plot_seg\n";
  if ($row == $plot_row && $col == $plot_col) {
    $plot_seg=$seg;
    print " Matched segment $plot_seg $seg\n";
    $no_match=1;
  }
}
close LIST;
  $out_file=$out_file.'.plot';
  print "Out File $out_file\n";
  open TEMP_OUT, ">$out_file";
  open TEMP_FILE,"$in_file";
  for $nd (1..$no_days) {
    for $ns (1..$nseg) {
#    print "$nd $ns Plot Seg = $plot_seg \n";
      ($xx=<TEMP_FILE>);
      (undef,$time,$nnd,$head,$ns,$nss,$temp,$T_head,$dbt,$flow,$depth)=split/\s+/,$xx;
      $year=int($time);
      $print_date=Date_to_Days($year,1,1)+$nnd-$array_date;
#    print "$nnd $head $nss $nseg\n";
      if ($ns == $plot_seg) {
      print "TEMP_OUT $plot_seg $nd\n";
        printf TEMP_OUT "%11.4f %4d %6d %7.2f %7.2f %7.2f %10.1f %5.1f\n",$time,$nnd,$print_date,$temp,$T_head,$dbt,$flow,$depth;
      }
    }
  }
  close TEMP_FILE;
  close TEMP_OUT;
}
