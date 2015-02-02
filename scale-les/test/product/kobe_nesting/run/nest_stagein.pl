#!/usr/bin/perl
#----------------------------------------------------------------
#   Anaysis for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#

$nmax = 432;
@prc_domains = ( 36, 108, 288 );

$init_num = "00010368000.000";

# make a process table
$color = 0;
$key   = 0;
$nprc  = $prc_domains[$color];
$prc_root = -999;

for($i=0; $i < $nmax; $i++){
   if($key == 0){
      $prc_root[$color] = $i;
   }
   $color_list[$i] = $color;
   $key_list[$i]   = $key;
   $kk = sprintf("%06d",$key);
   $dom = $color + 1;

   print "\#PJM --stgin  \'rank=".$i." ../init/domain_0".$dom."/init_d0".$dom."_".$init_num.".pe".$kk.".nc  \%r:./\' \n";
   print "\#PJM --stgin  \'rank=".$i." ../input/domain_0".$dom."/topo_d0".$dom.".pe".$kk.".nc  \%r:./\' \n";
   print "\#PJM --stgin  \'rank=".$i." ../input/domain_0".$dom."/landuse_d0".$dom.".pe".$kk.".nc  \%r:./\' \n";
   if($dom == 1){
      print "\#PJM --stgin  \'rank=".$i." ../init/domain_0".$dom."/boundary_d0".$dom.".pe".$kk.".nc  \%r:./\' \n";
   }

   $key++;
   if($key >= $nprc){
      $color++;
      $key = 0;
      $nprc  = $prc_domains[$color];
   }
}


#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
#EOF
