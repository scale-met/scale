#!/usr/bin/perl
#-------------------------------------------------------------------------------------------
#   AICS Climate Team Submission Tool  (Ver.0.2)
#   2011/12/22 --- Ryuji Yoshida.
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
use strict;
#
#+++++++++++++++++++++++++++++++++++  Do not change here  ++++++++++++++++++++++++++++++++++
our $public = "/work/user0173/public/";
our $user = "";
our $jobs = 0;
our $jid = 0;
our $uid = 0;
our $pid = 0;
our $jobtype;
our $interval = 180;
our @PEs;
our @script;
my  $id = 0;
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Main Routine

&setup();

for( $id=0; $id<$jobs; $id++ ){

    &group_start($id);
    #print "group ok\n";

    &user_start($id);
    #print "user ok\n";

    &casting($id);
    print " casting ok\n";
    sleep $interval;

    &job_quit($id);
    #print "job quit ok\n";
}

exit;


############################################################################################
sub setup{
    #
    my $i = 0;
    my $target = "";
    while( $i <= $#ARGV )
    {
        if( $ARGV[$i] eq "-u" ){
            $i++;
            $user = $ARGV[$i];
        }elsif( $ARGV[$i] eq "--user" ){
            $i++;
            $user = $ARGV[$i];
        }else{
            $target = $ARGV[$i];
        }
        $i++;
    }
    #
    if( $user eq "" || $target eq "" ){
        &script_help();
    }
    #
   if( -f $target ){
       print "user: ".$user."  target: ".$target."\n";
   }else{
       print "\n".$target." is not exist!!\n\n";
       exit;
    }
    #
    my $temp = `cat $target`;
    my @line = split(/\n/, $temp);
    $jobs = @line;
    if( $jobs == 0 ){
        print "\n !! Script is Empty !!\n\n";
        &script_help();
    }
    for($i = 0; $i < $jobs; $i++){
        my @tmp = split(/ /, $line[$i]);
        if($tmp[0] =~ /#/){
            $temp = $tmp[0];
            $temp =~ s/#/ /;
            $PEs[$i] = $temp;
            $script[$i] = $tmp[1];
            #print " PEs:".$PEs[$i]."  script: ".$script[$i]."\n";
	     if( $PEs[$i] eq "" || $PEs[$i] < 1 || $script[$i] eq "" ){
                print "\n !! Script Grammer is incorrect for Climsub !!\n\n";
                &script_help();
             }
        }else{
            print "\n !! Script Grammer is incorrect for Climsub !!\n\n";
            &script_help();
        }
    }
}
#-------------------------------------------------------------------------------------------
sub group_start{
    #
    my $idn = $_[0];
    my $i = 0;
    my $ni = 0;
    my $ctrl = 0;
    my $num = 0;
    my $max = 0;
    my $sign = "stop";
    my $header = "jobid";
    my $threshold = 6;
    my @tmp;
    my @tmp2;

    print " [START] PEs:".$PEs[$idn]."  script: ".$script[$idn]."\n";
    while($sign eq "stop"){
        my $work1 = $public.$header;
        my $temp = `ls -x $work1* `;
        $temp =~ s/\n/ /g;
        $temp =~ s/$public/ /g;
        @tmp = split(/ /, $temp);
        $ni = @tmp;
        $max = 0;
        $num = 0;
        for( $i=0; $i<=$ni; $i++ ){
            if($tmp[$i] =~ /$header/){
                if($tmp[$i] =~ /dummy/){
                }else{
                    $num++;
                    @tmp2 = split(/_/, $tmp[$i]);
                    $jid = $tmp2[1];
                    if($jid > $max){ $max = $jid; }
                }
            }
        }
        if($num >= $threshold){
             if($ctrl == 0){
                 print " --> pending (group)\n";
                 $ctrl = 1;
             }
             sleep $interval;
        }else{
              $sign = "go";
        }
    }
    $jid = $max + 1;
    if( $jid > 99999 ){ $jid = 0; }
    $jid = sprintf("%05d",$jid);
}
#-------------------------------------------------------------------------------------------
sub user_start{
    #
    my $idn = $_[0];
    my $i = 0;
    my $ni = 0;
    my $num = 0;
    my $ctrl = 0;
    my $max = 0;
    my $memo = "";
    my $sign = "stop";
    my $header = "jobid";
    my $threshold = 2;
    my @tmp;
    while($sign eq "stop"){
        my $work1 = $public.$header;
        my $temp = `ls -x $work1* `;
        $temp =~ s/\n/ /g;
        $temp =~ s/$public/ /g;
        @tmp = split(/ /, $temp);
        $ni = @tmp;
        $num = 0;
        for( $i=0; $i<=$ni; $i++ ){
            if($tmp[$i] =~ /$user/){ $num++; }
        }
        if($num >= $threshold){
             if($ctrl == 0){
                 print " --> pending (user)\n";
                 $ctrl = 1;
             }
             sleep $interval;
        }else{
              $sign = "go";
        }
    }
    my $date = `date`;
    my $jname = $header."_".$jid."_".$user;
    $memo = $date."\n JID:".$jid."\n USER:".$user."\n PEs:".$PEs[$idn]."\n SCRIPT: ".$script[$idn]."\n";
    system "rm -fv ".$jname." ";
    open(WRT, "> ".$public.$jname) || die "File Create Error JID\n";
    #print WRT $memo;
    close(WRT);
}
#-------------------------------------------------------------------------------------------
sub casting{
    #
    my $idn = $_[0];
    my $ctrl = 0;
    my $date;
    my $time;
    my $time1;
    my $time2;
    my $time_s1 = "Night Time Sec1 21:00:00 JST 2011";
    my $time_e1 = "Night Time Sec1 24:00:00 JST 2011";
    my $time_s2 = "Night Time Sec2 00:00:00 JST 2011";
    my $time_e2 = "Night Time Sec2 06:00:00 JST 2011";  # for safety, set to 6:00
    my $info;
    my $sign = "stop";
    my $threshold = 2304;
    my @tmp;

    if( $PEs[$idn] < $threshold ){
        $sign = "go";
        $jobtype = "S";
        print "--> Small Size JOB\n";
    }else{
        $jobtype = "L";
        print "--> Large Size JOB\n";
    }

    # Night Time Stopper
    while($sign eq "stop"){
        $date = `date`;
        $time = &to_sec($date);
        $time1 = &to_sec($time_s1);
        $time2 = &to_sec($time_e1);
        if($time1 < $time && $time <= $time2){ $sign = "go"; }
        $time1 = &to_sec($time_s2);
        $time2 = &to_sec($time_e2);
        if($time1 <= $time && $time < $time2){ $sign = "go"; }
        if( $sign ne "go" ){
             if($ctrl == 0){
                 print " Sleeping until Night Time!\n";
                 $ctrl = 1;
             }
             sleep $interval;
        }
    }

    # Setting JOB Process ID
    $info = `pjsub $script[$idn]`;
    my $r = $info =~ /(.+)(Job .+)/;
    @tmp = split(/ /, $2);
    $pid = $tmp[1];
    print " PJM Process ID: ".$pid."\n";
}
#-------------------------------------------------------------------------------------------
sub job_quit{
    #
    my $i;
    my $ni;
    my $sig = "OK";
    my $info;
    my $date;
    my $time;
    my $time1;
    my $time2;
    my $time_s1 = "Day Time Sec1 06:30:00 JST 2011";
    my $time_e1 = "Day Time Sec1 21:00:00 JST 2011";
    my $sign = "stop";
    my @tmp;
    my $header = "jobid";
    my $work1;

    while($sign eq "stop"){

        # Force Quit
        my $work1 = $public."*".$user;
        my $temp = `ls -x $work1* `;
        $temp =~ s/\n/ /g;
        $temp =~ s/$public/ /g;
        @tmp = split(/ /, $temp);
        $ni = @tmp;
        for( $i=0; $i<=$ni; $i++ ){
            if( $tmp[$i] =~ /quit/ || $tmp[$i] =~ /QUIT/ || $tmp[$i] =~ /Quit/ ){
                $date = `date`;
                print "\n !! Force Quit !!\n";
                print $date."\n\n";
                system `pjdel $pid`;
                my $jname = $public.$header."_".$jid."_".$user;
                system "rm -f $jname";
                exit;
            }
        }

        # Status Check
        $info = `pjstat`;
        @tmp = split(/\n/, $info);
        my $ln = @tmp;
        $sig = "OK";
        for($i=0; $i<=$ln; $i++){
            if($tmp[$i] =~ /$user/){
                if($tmp[$i] =~ /$pid/){
                    $sig = "NG";
                }
            }
        }
        if( $sig eq "NG" ){
            # Day Time Killer
            if($jobtype eq "L"){
                $date = `date`;
                $time = &to_sec($date);
                $time1 = &to_sec($time_s1);
                $time2 = &to_sec($time_e1);
                if($time1 < $time && $time < $time2){
                    print " !! Kill Large Size JOB at Day Time !!\n\n";
                    system `pjdel $pid`;
                    goto Label0;
                }
            }
            $sign = "stop";
            sleep $interval;
        }else{
            $sign = "go";
        }
    }

    Label0:print " finish: ".$header."_".$jid."\n";
    my $jname = $public.$header."_".$jid."_".$user;
    system "rm -f $jname";

}
#-------------------------------------------------------------------------------------------
sub to_sec{
    # This routine assumes following format: "Wed Dec 21 18:07:57 JST 2011"
    my $input = $_[0];
    my @tmp;
    my @tmp2;
    my $sec = 0;
    @tmp = split(/ /, $input);
    @tmp2 = split(/:/, $tmp[3]);
    $sec = $tmp2[0]*3600 + $tmp2[1]*60 + $tmp2[3];
    return $sec;
}
#-------------------------------------------------------------------------------------------
sub script_help{
        print "climsub.pl (version 0.1 / 2011-12-20)\n";
        print "\n";
        print "Usage: \n";
        print "  ./climsub.pl -u [user name] [script file] \n";
        print "\n";
        print "  -u, --user     casting user name\n";
        print "\n";
        print "\n";
        print "Example: \n";
        print "  ./climsub.pl -u user0173 job-script \n";
        print "\n";
        print "Job-script grammar:\n";
        print "  #[PEs number] [job script] \n";
        print "\n";
        print "Example: \n";
        print "  #1280 driver_1.sh\n";
        print "  #640  driver_2.sh\n";
        print "  #5120 driver_3.sh\n";
        exit;
}
#--------10--------20--------30--------40--------50--------60--------70--------80--------90#
# Script End
