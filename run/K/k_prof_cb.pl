#!/usr/bin/perl
# Auto Experiments Shell for scale3
# Ver.0.3
# 2012/01/28 --- Ryuji Yoshida.
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
#
# Variables Setting
@expNUM = ( 1,2,3,4 );                        # experiments number names (e.g. 1,2,3,4,5)
@expPES = ( "1x96", "2x96", "8x96", "32x96", "64x96", "128x96" );
@expGRS = ( "70x100x220" );
$homedir = "/work/user0173/scale/scaling/sc1";
$datadir_head = "/work/user0173/scale/scaling/sc1/data";
$bindir = "\${HMDIR}/bin/K";
#
# Parameters
$user = "user0173";
$num_threads = "8";
$lpgparam = "-s 32MB -d 32MB -h 32MB -t 32MB -p 32MB";
$stopper_head = "end_";
#
#------------------------------------------------------------------------------------------
#
print "|/--- start profiling collection... GO! GO! (^-^*)/\n|\n";
$stopper = $stopper_head."signal";

foreach $pes(@expPES){
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
foreach $grs(@expGRS){
   $datadir = $datadir_head."/grs".$grs."_pes".$pes;
   $basename = "grid_40m_g".$grs."_p".$pes;
   @tmp = split(/x/, $grs);
   $x = $tmp[0];
   $y = $tmp[1];
   $z = $tmp[2];

   $init_shell = "k_cbinit_g".$grs."_p".$pes.".sh";
   &mkinit;
   system "rm -f ".$init_shell." ";
   open(WRT, "> ".$init_shell) || die "File Open Error -mkinit-\n";
   print WRT $lines_init;
   close(WRT);

   print "|---> ".$init_shell."\n";
   #system "pjsub ".$init_shell." ";  [del]
   $info = `pjsub $init_shell`;
   $r = $info =~ /(.+)(Job .+)/;
   @tmp = split(/ /, $2);
   $pid = $tmp[1];
   print " PJM Process ID: ".$pid."\n";
   &breaf_stop($pid);
   sleep 5;   #5 wait for file handling
   system "mv ".$homedir."/".$init_shell."* ".$datadir;

   foreach $exp(@expNUM){
      $EX = sprintf("%02d",$exp);
      $outdir = $homedir."/output/grs".$grs."_pes".$pes."_ex".$EX;
      $basename_run = "grid_40m_g".$grs."_p".$pes."_ex".$EX;

      ###$ptime = "-1";                            # basic profiler mode
      for($ptime = 1; $ptime <= 7; $ptime++){  # detail profiler mode
      $fpcoll = &fpcoll_maker($ptime);

      $run_shell = "k_cb_g".$grs."_pe".$pes."_ex".$EX."_pa".$ptime.".sh";
      &mkrun;
      system "rm -f ".$run_shell." ";
      open(WRT, "> ".$run_shell) || die "File Open Error -mkrun-\n";
      print WRT $lines_run;
      close(WRT);

      print "|---> ".$run_shell."\n";
      #system "pjsub ".$run_shell." ";
      $info = `pjsub $run_shell`;
      $r = $info =~ /(.+)(Job .+)/;
      @tmp = split(/ /, $2);
      $pid = $tmp[1];
      print " PJM Process ID: ".$pid."\n";
      &breaf_stop($pid);
      sleep 5;   #5 wait for file handling
      system "mv ".$homedir."/".$run_shell."* ".$outdir;
      }  # loop-end: detail profiler mode

      &fpcoll_analysis();

   }

}
}

print "|\\--- all of sequences may be well completed... hope so...o(^-^)o\n";

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
export EXE=init_coldbubble
export OUTDIR=$datadir
#
# Run Command
#export RUN=\"mpiexec \$LPG -n $pes_sum \$BIN/\$EXE init_coldbubble.cnf\"
export RUN=\"mpiexec -n $pes_sum \$LPG \$BIN/\$EXE init_coldbubble.cnf\"
#
mkdir -p \$OUTDIR
cd \$OUTDIR
#
########################################################################
cat << End_of_SYSIN > \${OUTDIR}/\${EXE}.cnf

#####
#
# Scale3 benchmark configuration
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
 GRID_OUT_BASENAME = '$basename',
 GRID_DXYZ         = 40.D0,
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
 ATMOS_RESTART_OUT_BASENAME = 'init_coldbubble',
/

&PARAM_MKEXP_COLDBUBBLE
 XC_BBL =   1.4D3,
 YC_BBL =   2.0D3,
 ZC_BBL =   4.0D3,
 XR_BBL =   5.0D2,
 YR_BBL =   5.0D2,
 ZR_BBL =   5.0D2,
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
export EXE=scale3
export OUTDIR=$outdir
#
# Run Command
export RUN=\" \$FPC mpiexec -n $pes_sum \$LPG \$BIN/\$EXE scale3.cnf\"
#
mkdir -p \$OUTDIR
cd \$OUTDIR
#
########################################################################
cat << End_of_SYSIN > \${OUTDIR}/\${EXE}.cnf

#####
#
# Scale3 benchmark configuration
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
 TIME_DURATION              = 80.D0,
 TIME_DURATION_UNIT         = 'SEC',
 TIME_DT                    = 4.0D0,
 TIME_DT_UNIT               = 'SEC',
 TIME_DT_ATMOS_DYN          = 0.8D0,
 TIME_DT_ATMOS_DYN_UNIT     = 'SEC',
 TIME_NSTEP_ATMOS_DYN       = 10,
 TIME_DT_ATMOS_RESTART      = 80.D0,
 TIME_DT_ATMOS_RESTART_UNIT = 'SEC',
/

&PARAM_GRID
 GRID_OUT_BASENAME = '$basename_run',
 GRID_DXYZ         = 40.D0,
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
 ATMOS_RESTART_IN_BASENAME    = '$datadir/init_coldbubble_63072000000.000',
 ATMOS_RESTART_OUTPUT         = .false.,
 ATMOS_RESTART_OUT_BASENAME   = 'restart_coldbubble',
 ATMOS_RESTART_CHECK          = .false.,
 ATMOS_RESTART_CHECK_BASENAME = '$datadir/check_coldbubble_63072000040.000',
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TEMP_SFC    =   300.D0     
/
# ATMOS_REFSTATE_OUT_BASENAME = 'refstate',

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TAUZ = 75.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF = 1.D-3,
/

&PARAM_OCEAN
 OCEAN_TYPE = 'FIXEDSST',
/

&PARAM_OCEAN_VARS
 OCEAN_SST = 293.15D0,
/

&PARAM_HISTORY
 HISTORY_OUT_BASENAME      = 'history',
 HISTORY_DEFAULT_TINTERVAL = 80.D0,
 HISTORY_DEFAULT_TUNIT     = 'SEC',
 HISTORY_DEFAULT_AVERAGE   = .false.,
 HISTORY_DATATYPE          = 'REAL4',
/

#&HISTITEM item='DENS' /
#&HISTITEM item='MOMX' /
#&HISTITEM item='MOMY' /
#&HISTITEM item='MOMZ' /
#&HISTITEM item='RHOT' /
#&HISTITEM item='PRES' /
#&HISTITEM item='U'    /
#&HISTITEM item='V'    /
#&HISTITEM item='W'    /
#&HISTITEM item='T'    /
&HISTITEM item='PT'   /

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
   my $pdir = "pa";    # header name of directory for collected profile data

   if($fptime < 0){
      $fpline = "-Ibalance,call,cpu,hwm, -l20 -i20 -o Basic_Profile.txt -m 200000 ";   # for Basic Profiler
   }elsif($fptime == 1){
      $fpline = "-C -d ".$pdir."1 -Usection=range,local_event_number=0,29,29,29,30,5,9,6 -m 200000 ";
   }elsif($fptime == 2){
      $fpline = "-C -d ".$pdir."2 -Usection=range,local_event_number=30,30,30,8,29,30,31,0 -m 200000 ";
   }elsif($fptime == 3){
      $fpline = "-C -d ".$pdir."3 -Usection=range,local_event_number=31,10,11,30,31,0,30,30 -m 200000 ";
   }elsif($fptime == 4){
      $fpline = "-C -d ".$pdir."4 -Usection=range,local_event_number=0,12,48,48,2,32,48,48 -m 200000 ";
   }elsif($fptime == 5){
      $fpline = "-C -d ".$pdir."5 -Usection=range,local_event_number=7,7,7,32,0,13,13,22 -m 200000 ";
   }elsif($fptime == 6){
      $fpline = "-C -d ".$pdir."6 -Usection=range,local_event_number=0,13,13,13,13,35,35,33 -m 200000 ";
   }elsif($fptime == 7){
      $fpline = "-C -d ".$pdir."7 -Usection=range,local_event_number=35,35,26,0,32,7,7,31 -m 200000 ";
   }else{
      print " !!! fp_maker ERROR !!!\n";
      print " ---> fptime:".$fptime."\n\n";
      exit;
   }

   return $fpline;
} #----------------------------------------------------------------------------------------
#
# Subroutine "fpcoll_maker" -------------------------------------------------------------------
sub fpcoll_analysis(){
   my $pdir = "pa";    # header name of directory for collected profile data

   print "|\\--> [START] analysis collected data\n";

   print "||--> command: fprofpx   target: ".$pdir."1   output: output_prof_1.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."1 -o ".$outdir."/output_prof_1.csv -t csv -Ihwm";

   print "||--> command: fprofpx   target: ".$pdir."2   output: output_prof_2.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."2 -o ".$outdir."/output_prof_2.csv -t csv -Ihwm";

   print "||--> command: fprofpx   target: ".$pdir."3   output: output_prof_3.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."3 -o ".$outdir."/output_prof_3.csv -t csv -Ihwm";

   print "||--> command: fprofpx   target: ".$pdir."4   output: output_prof_4.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."4 -o ".$outdir."/output_prof_4.csv -t csv -Ihwm";

   print "||--> command: fprofpx   target: ".$pdir."5   output: output_prof_5.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."5 -o ".$outdir."/output_prof_5.csv -t csv -Ihwm";

   print "||--> command: fprofpx   target: ".$pdir."6   output: output_prof_6.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."6 -o ".$outdir."/output_prof_6.csv -t csv -Ihwm";

   print "||--> command: fprofpx   target: ".$pdir."7   output: output_prof_7.csv\n";
   system "fprofpx -d ".$outdir."/".$pdir."7 -o ".$outdir."/output_prof_7.csv -t csv -Ihwm";

   system "mkdir -p ".$outdir."/fpcoll_detail";
   system "mv ".$outdir."/".$pdir."* ".$outdir."/fpcoll_detail";
   system "mv ".$outdir."/output_prof_* ".$outdir."/fpcoll_detail";

   print "|/<-- [END] analysis collected data\n|\n";

   return;
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
      sleep 30;  # default=30 sec
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
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
