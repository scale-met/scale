#!/usr/bin/perl
# Auto Cruising Experiments for scale3
# Ver.0.1 (final)
# 2012/02/01 --- Ryuji Yoshida.
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
#
# Variables Setting
$Drun = 1400;           # interval of run (sec)
$Nrun = 6;            # iteration number of runs; Drun x Nrun = all duration time
$Rrun = 700;           # interval of Restart output (sec); set smaller value than Drun!
$grs =  "100x100x50";  # domain size for 1 PE
$resl = "400";         # resolution (m)
$pes =  "4x4";         # pjsub node shape
$restart = "false";   # restart run signal; true or false
$restart_step = 6;   # starting step of the restart run
#
$year = 2000;
$month = 1;
$day = 1;
$hour = 0;
$mint = 0;
$sec = 0;
#
# Parameters
$homedir = "/work/user0171/scale3";
$bindir = "/work/user0171/scale3/bin/K";
#
$inittype = "init_supercell";
$runtype = "Microphysics";
$user = "user0171";
$num_threads = "8";
$lpgparam = "-s 32MB -d 32MB -h 32MB -t 32MB -p 32MB";
$stopper_head = "end_";
#
#------------------------------------------------------------------------------------------
#

print "|/--- start auto cruising... GO! GO! (^-^*)/\n|\n";
system "touch ".$homedir."/end_dummy";

# cruising setting
$initdir = $homedir."/output/init_supercell/grs".$grs."_pes".$pes;
$sounding = $homedir."/data/supercell/input_sounding.txt";

@tmp = split(/x/, $grs);
$x = $tmp[0];
$y = $tmp[1];
$z = $tmp[2];

@tmp = split(/x/, $pes);
$tofu_x = $tmp[0];
$tofu_y = $tmp[1];
$pes_sum = sprintf("%02d",$tofu_x*$tofu_y);
if($pes_sum <= 3840){
   $rscgrp = "small";
}elsif($pes_sum > 3840){
   $rscgrp = "large";
}else{
   print " !!! Resource Group selection ERROR !!!\n";
   print " ---> PEs SUM:".$pes_sum."\n\n";
   exit;
}

$ticket = $Nrun;
$stopper = $stopper_head."signal";




# initial data make
if( $restart eq "true" ){
   $rest_sig = "true";
   print "|\\\n";
   print "||-- RESTART RUN !!\n";
   print "||-- starting step: ".$restart_step."\n";
   print "|/\n";
}else{
   $restart_step = 0;
   $rest_sig = "false";
   $basename_init = "grid_SPCL-".$resl."_g".$grs."_p".$pes;
   $init_shell = "k_SPCL_init_g".$grs."_p".$pes.".sh";
   &mkinit;
   system "rm -f ".$init_shell." ";
   open(WRT, "> ".$init_shell) || die "File Open Error -mkinit-\n";
   print WRT $lines_init;
   close(WRT);

   print "|---> ".$init_shell."\n";
   $info = `pjsub $init_shell`;
   $r = $info =~ /(.+)(Job .+)/;
   @tmp = split(/ /, $2);
   $pid = $tmp[1];
   print " PJM Process ID: ".$pid."\n";
   &breaf_stop($pid);
   sleep 5;   #5 wait for file handling
   system "mv ".$homedir."/".$init_shell."* ".$initdir;
}

for( $cru=0; $cru<$Nrun; $cru++ ){

      # running setup
      $CR = sprintf("%02d",($cru+$restart_step));
      $CR_1 = sprintf("%02d",($cru+$restart_step -1));
      $outdir = $homedir."/output/grs".$grs."_pes".$pes."_cru".$CR;
      $outdir_1 = $homedir."/output/grs".$grs."_pes".$pes."_cru".$CR_1;
      $basename_run = "grid_SPCL-".$resl."_g".$grs."_p".$pes."_cru".$CR;
      $run_shell = "k_SPCL_g".$grs."_pe".$pes."_cru".$CR.".sh";
      $ptime = "g".$grs."_p".$pes."_cru".$CR;   # basic profiler mode for cruising
      $fpcoll = &fpcoll_maker($ptime);
      $timestamp = &time_maker($cru);
      if( $rest_sig eq "true" ){
         $inputdata = &input_seacher(-99999);
      }else{
         $inputdata = &input_seacher($ticket);
      }

      &mkrun;
      system "rm -f ".$run_shell." ";
      open(WRT, "> ".$run_shell) || die "File Open Error -mkrun-\n";
      print WRT $lines_run;
      close(WRT);

      print "|---> ".$run_shell."\n";
      $info = `pjsub $run_shell`;
      $r = $info =~ /(.+)(Job .+)/;
      @tmp = split(/ /, $2);
      $pid = $tmp[1];
      print " PJM Process ID: ".$pid."\n";
      &breaf_stop($pid);
      sleep 5;   #5 wait for file handling
      system "mv ".$homedir."/".$run_shell."* ".$outdir;

      #&spd2bin();
      if( $rest_sig eq "true" ){
         $rest_sig = "false";
      }
      $ticket--;

}

system "rm ".$homedir."/end_dummy";

print "|--- ticket:".$ticket."\n|\n";
print "|\\--- arrive the goal... hope so...o(^-^)o\n";

############################################################################################
#End of Main Routine

#--------10--------20--------30--------40--------50--------60--------70--------80--------90#

# Make INIT script -------------------------------------------------------------------------
sub mkinit(){
$lines_init = "#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list \"rscgrp=$rscgrp\"
#PJM --rsc-list \"node=$pes\"
#PJM --rsc-list \"elapse=01:00:00\"
#PJM --rsc-list \"node-mem=12Gi\"
##PJM --mpi \"shape=$pes\"
##PJM --mpi \"proc=$pes_sum\"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=$num_threads
export OMP_NUM_THREADS=$num_threads
export FPC=\"fpcoll $fpcoll \"
export LPG=\"lpgparm $lpgparam \"
#
export HMDIR=$homedir
export BIN=$bindir
export EXE=$inittype
export OUTDIR=$initdir
#
# Run Command
export RUN=\"mpiexec -n $pes_sum \$LPG \$BIN/\$EXE $inittype.cnf\"
#
mkdir -p \$OUTDIR
cd \$OUTDIR
#
########################################################################
cat << End_of_SYSIN > \${OUTDIR}/\${EXE}.cnf

#####
#
# Scale3 init_supercell configuration
#
#####

&PARAM_PRC
 PRC_NUM_X       = $tofu_x,
 PRC_NUM_Y       = $tofu_y,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_CONST
/

&PARAM_TIME
 TIME_STARTDATE             = 2000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/

&PARAM_GRID
 GRID_OUT_BASENAME = '$basename_init',
 GRID_DXYZ         = $resl.D0,
 GRID_KMAX         = $z,
 GRID_IMAX         = $x,
 GRID_JMAX         = $y,
 GRID_BUFFER_DZ    = 4.0D3,
 GRID_BUFFER_DX    = 0.0D0,
 GRID_BUFFER_DY    = 0.0D0,
 GRID_BUFFFACT     = 1.0D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = 'fent_fct',
 ATMOS_TYPE_PHY_TB = 'smagorinsky',
 ATMOS_TYPE_PHY_MP = 'NDW6',
 ATMOS_TYPE_PHY_RD = 'mstrnX',
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX            = 11,
 ATMOS_RESTART_OUTPUT       = .true.,
 ATMOS_RESTART_OUT_BASENAME = '$inittype',
/

&PARAM_MKEXP_SUPERCELL
 ENV_IN_SOUNDING_file = '$sounding',
 ZC_BBL =  1.5D3,
 XC_BBL = 80.0D3,
 YC_BBL = 80.0D3,
 ZR_BBL =  1.0D3,
 XR_BBL =  8.0D3,
 YR_BBL =  8.0D3,
/

End_of_SYSIN
########################################################################

# run
echo \"job \${RUNNAME} started at \" `date`
\$RUN > \$OUTDIR/STDOUT_INIT 2>&1
echo \"job \${RUNNAME} end     at \" `date`

touch ".$homedir."/".$stopper."

exit
";
 return;
} #----------------------------------------------------------------------------------------
#
# Make RUN script -------------------------------------------------------------------------
sub mkrun(){
$lines_run = "#! /bin/bash -x
#
# for K Computer
#
#PJM --rsc-list \"rscgrp=$rscgrp\"
#PJM --rsc-list \"node=$pes\"
#PJM --rsc-list \"elapse=01:00:00\"
#PJM --rsc-list \"node-mem=12Gi\"
##PJM --mpi \"shape=$pes\"
##PJM --mpi \"proc=$pes_sum\"
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=$num_threads
export OMP_NUM_THREADS=$num_threads
export FPC=\"fpcoll $fpcoll \"
export LPG=\"lpgparm $lpgparam \"
#
export HMDIR=$homedir
export BIN=$bindir
export EXE=$runtype
export OUTDIR=$outdir
#
# Run Command
export RUN=\" \$FPC mpiexec -n $pes_sum \$LPG \$BIN/\$EXE $runtype.cnf\"
#export RUN=\" mpiexec -n $pes_sum \$LPG \$BIN/\$EXE $runtype.cnf\"
#
mkdir -p \$OUTDIR
cd \$OUTDIR
#
########################################################################
cat << End_of_SYSIN > \${OUTDIR}/\${EXE}.cnf

#####
#
# Scale3 supercell configuration
#
#####

&PARAM_PRC
 PRC_NUM_X       = $tofu_x,
 PRC_NUM_Y       = $tofu_y,
 PRC_PERIODIC_X  = .true.,
 PRC_PERIODIC_Y  = .true.,
/

&PARAM_CONST
/

&PARAM_TIME
 TIME_STARTDATE             = $timestamp
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = $Drun.D0,
 TIME_DURATION_UNIT         = 'SEC',
 TIME_DT                    = 0.7D0,
 TIME_DT_UNIT               = 'SEC',
 TIME_DT_ATMOS_DYN          = 0.7D0,
 TIME_DT_ATMOS_DYN_UNIT     = 'SEC',
 TIME_NSTEP_ATMOS_DYN       = 1,
 TIME_DT_ATMOS_PHY_MP       = 0.7D0,
 TIME_DT_ATMOS_PHY_MP_UNIT  = 'SEC',
 TIME_DT_ATMOS_RESTART      = $Rrun.D0,
 TIME_DT_ATMOS_RESTART_UNIT = 'SEC',
/

&PARAM_GRID
 GRID_OUT_BASENAME = '$basename_run',
 GRID_DXYZ         = $resl.D0,
 GRID_KMAX         = $z,
 GRID_IMAX         = $x,
 GRID_JMAX         = $y,
 GRID_BUFFER_DZ    = 4.0D3,
 GRID_BUFFER_DX    = 0.0D0,
 GRID_BUFFER_DY    = 0.0D0,
 GRID_BUFFFACT     = 1.0D0,
/

&PARAM_ATMOS
 ATMOS_TYPE_DYN    = 'fent_fct',
 ATMOS_TYPE_PHY_TB = 'smagorinsky',
 ATMOS_TYPE_PHY_MP = 'NDW6',
 ATMOS_TYPE_PHY_RD = 'mstrnX',
/

&PARAM_ATMOS_VARS
 ATMOS_QTRC_NMAX              = 11,
 ATMOS_RESTART_IN_BASENAME    = '".$inputdata."',
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = 'restart_supercell',
 ATMOS_RESTART_CHECK          = .false.,
 ATMOS_RESTART_CHECK_BASENAME = '$datadir/check_supercell_63072000040.000',
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TEMP_SFC    =   300.D0     
/
# ATMOS_REFSTATE_OUT_BASENAME = 'refstate',

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_OUT_BASENAME = '',
 ATMOS_BOUNDARY_VALUE_VELX   = 0.D0,
 ATMOS_BOUNDARY_TAUZ = 75.D0,
 ATMOS_BOUNDARY_TAUX =  0.D0,
 ATMOS_BOUNDARY_TAUY =  0.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-2,
/

&PARAM_OCEAN
 OCEAN_TYPE = 'FIXEDSST',
/

&PARAM_OCEAN_VARS
 OCEAN_SST = 293.15D0,
/

&PARAM_HISTORY
 HISTORY_OUT_BASENAME      = 'history',
 HISTORY_DEFAULT_TINTERVAL = 7.0D0,
 HISTORY_DEFAULT_TUNIT     = 'SEC',
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = 'REAL4',
/

&HISTITEM item='DENS' /
#&HISTITEM item='MOMX' /
#&HISTITEM item='MOMY' /
#&HISTITEM item='MOMZ' /
#&HISTITEM item='RHOT' /

#&HISTITEM item='PRES' /
&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
#&HISTITEM item='T'    /
&HISTITEM item='PT'   /

&HISTITEM item='QV'   /
&HISTITEM item='QC'   /
&HISTITEM item='QR'   /
&HISTITEM item='QI'   /
&HISTITEM item='QS'   /
&HISTITEM item='QG'   /

End_of_SYSIN
########################################################################

# run
echo \"job \${RUNNAME} started at \" `date`
\$RUN > \$OUTDIR/STDOUT_RUN 2>&1
echo \"job \${RUNNAME} end     at \" `date`

touch ".$homedir."/".$stopper."

exit
";
 return;
} #----------------------------------------------------------------------------------------
#
# Subroutine "fpcoll_maker" -------------------------------------------------------------------
sub fpcoll_maker($){
   my $fptime = $_[0];
   my $fpline = " ";
#   my $pdir = "pa";    # header name of directory for collected profile data

#   if($fptime < 0){
#      $fpline = "-Ibalance,call,cpu,hwm, -l20 -i20 -o Basic_Profile.txt -m 200000 ";   # for Basic Profiler
      $fpline = "-Ihwm,cpu -l0 -o Basic_Profile_$fptime.txt -m 200000 ";            # for Basic Profiler (simple)
#   }elsif($fptime == 1){
#      $fpline = "-C -d ".$pdir."1 -Usection=range,local_event_number=0,29,29,29,30,5,9,6 -m 200000 ";
#   }elsif($fptime == 2){
#      $fpline = "-C -d ".$pdir."2 -Usection=range,local_event_number=30,30,30,8,29,30,31,0 -m 200000 ";
#   }elsif($fptime == 3){
#      $fpline = "-C -d ".$pdir."3 -Usection=range,local_event_number=31,10,11,30,31,0,30,30 -m 200000 ";
#   }elsif($fptime == 4){
#      $fpline = "-C -d ".$pdir."4 -Usection=range,local_event_number=0,12,48,48,2,32,48,48 -m 200000 ";
#   }elsif($fptime == 5){
#      $fpline = "-C -d ".$pdir."5 -Usection=range,local_event_number=7,7,7,32,0,13,13,22 -m 200000 ";
#   }elsif($fptime == 6){
#      $fpline = "-C -d ".$pdir."6 -Usection=range,local_event_number=0,13,13,13,13,35,35,33 -m 200000 ";
#   }elsif($fptime == 7){
#      $fpline = "-C -d ".$pdir."7 -Usection=range,local_event_number=35,35,26,0,32,7,7,31 -m 200000 ";
#   }else{
#      print " !!! fp_maker ERROR !!!\n";
#      print " ---> fptime:".$fptime."\n\n";
#      exit;
#   }

   return $fpline;
} #----------------------------------------------------------------------------------------
#
# Subroutine "breaf_stop" -------------------------------------------------------------------
sub breaf_stop($){
   my $pjmid = $_[0];
   my $ii = 0;
   my $counter = 0;
   my $targ = $homedir."/".$stopper_head."*";
   my $targ2 = $homedir."/".$stopper;
   my $homedir2 = $homedir."/";
   my $sign = "stop";
   my $sign1 = "NG";
   my $sign2 = "NG";

   # Step to Next, when exist "signal file"
   #                    or process id is deleted from pjstat
   while($sign eq "stop"){
      print "\r|-".$counter;
      $counter++;
      sleep 60;  # default=60 sec
      # Signal File Check
      my $temp = `ls -x $targ `;
      $temp =~ s/\n/ /g;
      $temp =~ s/$homedir2/ /g;
      @tmp = split(/ /, $temp);
      my $ni = @tmp;
      for( $ii=0; $ii<=$ni; $ii++ ){
         if($tmp[$ii] =~ /$stopper/){
            $sign1 = "OK";
         }
      }
      # Status Check
      $info = `pjstat`;
      @tmp = split(/\n/, $info);
      my $ln = @tmp;
      $sign2 = "OK";
      for($ii=0; $ii<=$ln; $ii++){
          if($tmp[$ii] =~ /$user/){
              if($tmp[$ii] =~ /$pjmid/){
                  $sign2 = "NG";
              }
          }
      }
      if($sign1 eq "OK" || $sign2 eq "OK"){
         $sign = "go";
         if($sign1 eq "OK"){
            print "-quit by signal-";
         }elsif($sign2 eq "OK"){
            print "-quit by pjstat-";
         }
      }

   }
   print "->next\n";
   print "|\n";
   system "rm -r $targ2";

   return;
} #----------------------------------------------------------------------------------------
#
# Subroutine "inputdata name & directory search" -------------------------------------------------------------------
sub input_seacher($){
   my $tic = $_[0];
   my $ii = 0;
   my $counter = 0;
   my $input = "0";
   my $datadir_before = $datadir;

   if( $tic == $Nrun ){
      $datadir = $initdir;
      $datadir_before = $datadir;
      $inittype = "init_supercell";
   }else{
      $datadir = $outdir;
      $datadir_before = $outdir_1;
      $inittype = "restart_supercell";
   }

   my $targ = $datadir_before."/LOG.pe000000";
   print $targ."\n";
   open(LOG, "<".$targ) || die "Data Open Error! -LOG.pe000000-\n";
      my @line = <LOG>;
   close(LOG);
   my $ln = @line;
   for( $ii=0; $ii<=$ln; $ii++ ){
      if($line[$ii] =~ /$inittype/ && $line[$ii] =~ /.pe000000/ && $line[$ii] =~ / \*\*\* filename:/){
          $r = $line[$ii] =~ /( \*\*\* filename: )(.+)(.pe.+)/;
          $input = $2;
      }
   }
   print "|--- next input data: ".$input."\n";
   my $input2 = $datadir_before."/".$input;

   return $input2;
} #----------------------------------------------------------------------------------------
#
# Subroutine "time stamp maker" -------------------------------------------------------------------
sub time_maker($){
   my $tt = $_[0];
   my $ii = 0;
   my $time_stp = "0";
   my $yyyy = 0;
   my $mm = 0;
   my $dd = 0;
   my $hh = 0;
   my $mi = 0;
   my $se = 0;

   if( $tt == 0 ){
      $yyyy = sprintf("%4d",$year);
      $mm = sprintf("%2d",$month);
      $dd = sprintf("%2d",$day);
      $hh = sprintf("%2d",$hour);
      $mi = sprintf("%2d",$mint);
      $se = sprintf("%2d",$sec);
      $time_stp = $yyyy.",".$mm.",".$dd.",".$hh.",".$mi.",".$se.",";
   }else{
      # integrate SEC
      $sec = $sec + $Drun;
      while( $sec >= 60 ){
          $sec = $sec - 60;
          # integrate MIN
          $mint = $mint + 1;
          if( $mint >= 60 ){
             $mint = 0;
             # integrate HOUR
             $hour = $hour + 1;
             if( $hour >= 24 ){
                $hour = 0;
                # integrate DAY
                $day = $day + 1;
             }
          }
      }
      $yyyy = sprintf("%4d",$year);
      $mm = sprintf("%2d",$month);
      $dd = sprintf("%2d",$day);
      $hh = sprintf("%2d",$hour);
      $mi = sprintf("%2d",$mint);
      $se = sprintf("%2d",$sec);
      $time_stp = $yyyy.",".$mm.",".$dd.",".$hh.",".$mi.",".$se.",";
   }

   print "|--- next start time [ ".$time_stp." ]\n";

   return $time_stp;
} #----------------------------------------------------------------------------------------
#
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
