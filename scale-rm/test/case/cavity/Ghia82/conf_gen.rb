#!/bin/env ruby

require "fileutils"

##################################################################

Ma    = 1.5E-1
PRES0 = 1.0E05
Ulid  = 5.0E01
L     = 2.0E04
KDiff = 1.0E03

Rd    = 287.04E0
Cpd   = 1004.64E0
T0    = Ulid**2 / ( Ma**2 * Cpd/(Cpd-Rd) * Rd )

Re    = Ulid * L / KDiff

#################################################################

TIME_DURATION_SEC       = "864000.D0"
HISTORY_TINTERVAL_SEC   = "7200.0D0"
CONF_GEN_RESOL_HASHLIST = \
[ \
  { "TAG"=>"2000m", "DX"=>2000E0, "DY"=> 2000E0, "DZ"=>2000E0, 
    "IMAX"=>10, "JMAX"=>10, "KMAX"=>1, "DT"=>2.4, "DTDYN"=>2.4E0, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"1000m", "DX"=>1000E0, "DY"=> 1000E0, "DZ"=>1000E0, 
    "IMAX"=>20, "JMAX"=>20, "KMAX"=>1, "DT"=>1.2, "DTDYN"=>1.2E0, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"500m", "DX"=>500E0, "DY"=> 500E0, "DZ"=>500E0, 
    "IMAX"=>40, "JMAX"=>40, "KMAX"=>1, "DT"=>0.6, "DTDYN"=>0.6E0, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"250m", "DX"=>250E0, "DY"=> 250E0, "DZ"=>250.0E0,
    "IMAX"=>40, "JMAX"=>40, "KMAX"=>1, "DT"=>0.3, "DTDYN"=>0.3E0, "NPRCX"=> 2, "NPRCY"=>2}, \
  { "TAG"=>"125m", "DX"=>125E0, "DY"=> 125E0, "DZ"=>125.0E0,
    "IMAX"=>40, "JMAX"=>40, "KMAX"=>1, "DT"=>0.15, "DTDYN"=>0.15E0, "NPRCX"=> 4, "NPRCY"=>4}, \
  { "TAG"=>"063m", "DX"=> 62.5E0, "DY"=>  62.5E0, "DZ"=> 62.5E0,
    "IMAX"=>40, "JMAX"=>40, "KMAX"=>1, "DT"=>0.12, "DTDYN"=>0.12E0, "NPRCX"=> 8, "NPRCY"=>8}, \
]
CONF_GEN_CASE_HASH_LIST = \
[ \
  {"TAG"=>"CTRL"}, \
]
CONF_GEN_NUMERIC_HASHLIST = \
[ \
  {"TAG"=>"FVM_CD2"}, {"TAG"=>"FVM_CD4"}, {"TAG"=>"FVM_CD6"},  \
  {"TAG"=>"FVM_UD1"}, {"TAG"=>"FVM_UD3"}, {"TAG"=>"FVM_UD5"},  \
]


#########################################################

def gen_init_conf( conf_name,
                   nprocx, nprocy, imax, jmax, kmax, dx, dy, dz )
  f = File.open(conf_name, "w")
  f.print <<EOS
######################################################
#
# SCALE-RM mkinit configulation for cavity flow
# (following experimental setup in Ghia et al. 1982)
#
# P0              = #{PRES0}   [Pa]
# T0              = #{T0}      [K]
# Ulid            = #{Ulid}    [m/s]
# Mach number     = #{Ma}
# Reynolds number = #{Re}
#
######################################################

&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
 PRC_PERIODIC_X  = .false.,
 PRC_PERIODIC_Y  = .false.,
/

&PARAM_INDEX
 KMAX = #{kmax}, 
 IMAX = #{imax}, IHALO = 3, 
 JMAX = #{jmax}, JHALO = 3,
/

&PARAM_GRID
 DZ =  #{dz}, 
 DX =  #{dx},  
 DY =  #{dy}, 
/

&PARAM_TIME
 TIME_STARTDATE             = 0000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .true.,
 STATISTICS_use_globalcomm = .true.,
/

&PARAM_TRACER
 TRACER_TYPE = 'DRY',
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init",
/

&PARAM_CONST
 CONST_GRAV      =   0.00000000000000     ,
/

&PARAM_MKINIT
 MKINIT_initname = "CAVITYFLOW",
/

&PARAM_MKINIT_CAVITYFLOW
  Ulid         = #{Ulid}, 
  PRES0        = #{PRES0}, 
  REYNOLDS_NUM = #{Re},
  MACH_NUM     = #{Ma}
/

!&PARAM_ATMOS_HYDROSTATIC
! HYDROSTATIC_uselapserate = .false.,
!/

EOS
  f.close
  
end

def gen_run_conf( conf_name,
                  nprocx, nprocy,
                  imax, jmax, kmax, dx, dy, dz,
                  dtsec, dtsec_dyn,
                  flxEvalType, fctFlag, dataDir )

  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-RM run configulation
#
# P0              = #{PRES0}   [Pa]
# T0              = #{T0}      [K]
# Ulid            = #{Ulid}    [m/s]
# Mach number     = #{Ma}
# Reynolds number = #{Re}
#
#####

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
 PRC_PERIODIC_X  = .false.,
 PRC_PERIODIC_Y  = .false.,
/

&PARAM_INDEX
 KMAX = #{kmax}, 
 IMAX = #{imax}, IHALO = 3, 
 JMAX = #{jmax}, JHALO = 3,
/

&PARAM_GRID
 DZ =  #{dz}, 
 DX =  #{dx},  
 DY =  #{dy}, 
/


&PARAM_TIME
 TIME_STARTDATE             = 0000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = #{TIME_DURATION_SEC},
 TIME_DURATION_UNIT         = "SEC",
 TIME_DT                    = #{dtsec},
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = #{dtsec_dyn}, 
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .false.,
 STATISTICS_use_globalcomm = .false.,
/

&PARAM_TRACER
 TRACER_TYPE = 'DRY',
/

&PARAM_CONST
 CONST_GRAV      =   0.00000000000000     ,
/

!&PARAM_ATMOS_HYDROSTATIC
! HYDROSTATIC_uselapserate = .false.,
!/

&PARAM_ATMOS
 ATMOS_DYN_TYPE    = "FVM-HEVE",
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_IN_BASENAME      = "init_00000101-000000.000",
 ATMOS_RESTART_OUTPUT           = .false.,
 ATMOS_VARS_CHECKRANGE          = .true.,
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE       = "INIT",
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_OUT_BASENAME = "boundary"
 ATMOS_BOUNDARY_TYPE       = "INIT",
 ATMOS_BOUNDARY_USE_VELZ   = .true.,
 ATMOS_BOUNDARY_USE_VELX   = .true.,
 ATMOS_BOUNDARY_USE_VELY   = .true.,
 ATMOS_BOUNDARY_TAUZ       = 1.D0,
 ATMOS_BOUNDARY_TAUX       = 1.D0,
 ATMOS_BOUNDARY_TAUY       = 1.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_TINTEG_LARGE_TYPE = "EULER",
 ATMOS_DYN_TINTEG_SHORT_TYPE = "RK3WS2002",
 ATMOS_DYN_TINTEG_TRACER_TYPE = "RK3WS2002",
 ATMOS_DYN_FVM_FLUX_TYPE        = "#{flxEvalType}",             
 ATMOS_DYN_FVM_FLUX_TRACER_TYPE = "#{flxEvalType}", 
 ATMOS_DYN_NUMERICAL_DIFF_COEF  = 0.D0,
 ATMOS_DYN_DIVDMP_COEF          = 5.D-4,
 ATMOS_DYN_FLAG_FCT_TRACER      = #{fctFlag}, 
/

&PARAM_USER
 USER_do = .true., 
 Ulid    = #{Ulid}, 
 PRES0   = #{PRES0}, 
 KDiff   = #{KDiff}
/

&PARAM_HISTORY
 HISTORY_DEFAULT_BASENAME  = "history",
 HISTORY_DEFAULT_TINTERVAL = #{HISTORY_TINTERVAL_SEC},
 HISTORY_DEFAULT_TUNIT     = "SEC",
 HISTORY_DEFAULT_TAVERAGE  = .false.,
 HISTORY_DEFAULT_DATATYPE  = "REAL4",
 HISTORY_OUTPUT_STEP0      = .true.,
/

&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='VOR'    /
!&HISTITEM item='W'    /
&HISTITEM item='PT'   /
&HISTITEM item='PRES'   /
&HISTITEM item='DENS'   /
!&HISTITEM item='MOMX'   /
!&HISTITEM item='MOMY'   /
!&HISTITEM item='MOMZ'   /
!&HISTITEM item='T'   /


&PARAM_MONITOR
 MONITOR_STEP_INTERVAL = 120,
/

&MONITITEM item='QDRY' /
&MONITITEM item='QTOT' /
&MONITITEM item='ENGT' /
&MONITITEM item='ENGP' /
&MONITITEM item='ENGK' /
&MONITITEM item='ENGI' /
EOS
f.close
end

CONF_GEN_RESOL_HASHLIST.each{|resol_hash|
  CONF_GEN_CASE_HASH_LIST.each{|case_hash|
    CONF_GEN_NUMERIC_HASHLIST.each{|numeric_hash|
      ["F", "T"].each{|fct_flag|

        dataDir = "./#{resol_hash["TAG"]}/#{case_hash["TAG"]}/"
        dataDir += fct_flag=="T" ? "#{numeric_hash["TAG"]}_FCT/" : "#{numeric_hash["TAG"]}/"

        puts "Generate init.conf and run.conf (Dir=#{dataDir}) .."
        if !File.exists?(dataDir) then
          puts "Create directory .."
          FileUtils.mkdir_p(dataDir)
        end
      
        init_conf_name = "#{dataDir}init.conf" 
        gen_init_conf(init_conf_name, 
                      resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["JMAX"], resol_hash["KMAX"], 
                      resol_hash["DX"], resol_hash["DY"], resol_hash["DZ"]  )

        run_conf_name = "#{dataDir}run.conf"
        gen_run_conf(run_conf_name, 
                     resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["JMAX"], resol_hash["KMAX"], 
                     resol_hash["DX"], resol_hash["DY"], resol_hash["DZ"], resol_hash["DT"], resol_hash["DTDYN"], 
                     numeric_hash["TAG"].sub("FVM_",""), fct_flag, dataDir )
      }
    }
  }
}
  
