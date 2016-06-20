#!/bin/env ruby

TIME_DT_SEC             = "5.0D0"
TIME_DURATION_SEC       = "864000.D0"
HISTORY_TINTERVAL_SEC   = "7200.0D0"
CONF_GEN_RESOL_HASHLIST = \
[ \
  { "TAG"=>"500m", "DX"=>500E0, "DZ"=>500E0, 
    "KMAX"=>40, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>5.0, "NPRCX"=> 2, "NPRCY"=>1}, \
  { "TAG"=>"250m", "DX"=>250E0, "DZ"=>250.0E0,
    "KMAX"=>80, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>2.5E0, "NPRCX"=> 4, "NPRCY"=>1}, \
#  { "TAG"=>"100m", "DX"=>100E0, "DZ"=>100.0E0,
#    "KMAX"=>64, "IMAX"=>64, "JMAX"=>3, "DTDYN"=>0.125E0, "NPRCX"=> 8, "NPRCY"=>1}, \
#  { "TAG"=>"050m", "DX"=>50E0, "DZ"=>50.0E0,
#    "KMAX"=>128, "IMAX"=>64, "JMAX"=>3, "DTDYN"=>0.0625E0, "NPRCX"=>16, "NPRCY"=>1}, \
]
CONF_GEN_NUMERIC_HASHLIST = \
[ \
  {"TAG"=>"FDM_CD2"}, {"TAG"=>"FDM_CD4"}, {"TAG"=>"FDM_CD6"},  \
  {"TAG"=>"FDM_UD1"}, {"TAG"=>"FDM_UD3"}, {"TAG"=>"FDM_UD5"},  \
]

#########################################################

def gen_init_conf(conf_name, nprocx, nprocy, imax, kmax, dx, dz)
  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-RM mkinit configulation for coldbubble test
# (following experimental setup in Ghia et al. 1982)
#
#####

&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
 PRC_PERIODIC_X  = .false., 
/

&PARAM_INDEX
 KMAX = #{kmax}, 
 IMAX = #{imax}, 
 JMAX = 3,
/

&PARAM_GRID
 DZ =  #{dz}, 
 DX =  #{dx},  
 DY =  #{dx}, 
 BUFFER_DZ =   0.D0,  
 BUFFFACT  =   1.D0,
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

!&PARAM_ATMOS_HYDROSTATIC
! HYDROSTATIC_uselapserate = .false.,
!/

EOS
  f.close
  
end

def gen_run_conf(conf_name, nprocx, nprocy, imax, kmax, dx, dz, dtsec_dyn, flxEvalType, dataDir)

  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-RM run configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
 PRC_PERIODIC_X  = .false., 
/

&PARAM_INDEX
 KMAX = #{kmax}, 
 IMAX = #{imax}, 
 JMAX = 3,
/

&PARAM_GRID
 DZ =  #{dz}, 
 DX =  #{dx},  
 DY =  #{dx}, 
 BUFFER_DZ =   0.D0,  
 BUFFFACT  =   1.D0,
/
&PARAM_TIME
 TIME_STARTDATE             = 0000, 1, 1, 0, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = #{TIME_DURATION_SEC},
 TIME_DURATION_UNIT         = "SEC",
 TIME_DT                    = #{TIME_DT_SEC},
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = #{dtsec_dyn}, 
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .false.,
 STATISTICS_use_globalcomm = .true.,
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
 ATMOS_DYN_TYPE    = "FDM-HEVE",
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_IN_BASENAME      = "init_00000000000.000",
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
 ATMOS_DYN_NUMERICAL_DIFF_ORDER = 0,     ! Use 2nd-order diffusion!
 ATMOS_DYN_NUMERICAL_DIFF_COEF  = 20.D0,
 ATMOS_DYN_DIVDMP_COEF   = 0.001D0,
 ATMOS_DYN_FLXEVAL_TYPE  = "#{flxEvalType}"
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
&HISTITEM item='W'    /
&HISTITEM item='PT'   /
&HISTITEM item='PRES'   /
&HISTITEM item='DENS'   /
&HISTITEM item='MOMX'   /
&HISTITEM item='MOMY'   /
&HISTITEM item='MOMZ'   /
!&HISTITEM item='T'   /


&PARAM_MONITOR
 MONITOR_STEP_INTERVAL = 12,
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
  CONF_GEN_NUMERIC_HASHLIST.each{|numeric_hash|
    dataDir = "./#{resol_hash["TAG"]}/#{numeric_hash["TAG"]}/"

    puts "generate init.conf and run.conf (Dir=#{dataDir})"

    init_conf_name = "#{dataDir}init.conf" 
    gen_init_conf(init_conf_name, \
                 resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["KMAX"], \
                 resol_hash["DX"], resol_hash["DZ"] )
    run_conf_name = "#{dataDir}run.conf"
    gen_run_conf(run_conf_name, \
                 resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["KMAX"], \
                 resol_hash["DX"], resol_hash["DZ"], resol_hash["DTDYN"], numeric_hash["TAG"], dataDir )
  }
}
