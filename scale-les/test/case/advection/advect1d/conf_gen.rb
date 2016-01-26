#!/bin/env ruby

TIME_DT_SEC             = "0.25D0"
TIME_DURATION_SEC       = "600.D0"
HISTORY_TINTERVAL_SEC   = "50.D0"
CONF_GEN_RESOL_HASHLIST = \
[ \
  { "TAG"=>"500m", "DX"=>500E0, "DZ"=>500.0E0, 
    "KMAX"=>4, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>0.25, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"250m", "DX"=>250E0, "DZ"=>250.0E0, 
    "KMAX"=>4, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>0.125, "NPRCX"=> 2, "NPRCY"=>1}, \
  { "TAG"=>"125m", "DX"=>125E0, "DZ"=>125.0E0, 
    "KMAX"=>4, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>0.0625, "NPRCX"=> 4, "NPRCY"=>1}, \
#  { "TAG"=>"500m", "DX"=>500E0, "DZ"=>500.0E0, 
#    "KMAX"=>3, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>0.6E0, "NPRCX"=> 1, "NPRCY"=>1}, \
]
CONF_GEN_CASE_HASH_LIST = \
[ \
  {"TAG"=>"RECT", "SHAPE_NC"=>"RECT"}, \
  {"TAG"=>"COSBELL", "SHAPE_NC"=>"BUBBLE"} \
]
CONF_GEN_NUMERIC_HASHLIST = \
[ \
  {"TAG"=>"FDM_CD2"}, {"TAG"=>"FDM_CD4"}, {"TAG"=>"FDM_CD6"},  \
  {"TAG"=>"FDM_UD1"}, {"TAG"=>"FDM_UD3"}, {"TAG"=>"FDM_UD5"},  \
]

#########################################################

def gen_init_conf(conf_name, nprocx, nprocy, imax, kmax, dx, dz, shape_nc)
  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-LES mkinit configulation for coldbubble test
# (following experimental setup in Straka et al. 1993)
#
#####

&PARAM_IO
 IO_LOG_BASENAME = 'init_LOG',
/

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
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

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init",
/

&PARAM_TRACER
 TRACER_TYPE = 'SN14',
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init",
/

&PARAM_MKINIT
 MKINIT_initname = "ADVECT",
/

&PARAM_BUBBLE
 BBL_CZ = 10.0D3,
 BBL_CX = 12.0D3,
 BBL_CY = 12.0D3,
 BBL_RZ = 1.0D14,
 BBL_RX = 2.0D3,
 BBL_RY = 1.0D14,
/

&PARAM_RECT
 RCT_CZ = 10.0D3,
 RCT_CX = 12.0D3,
 RCT_CY = 12.0D3,
 RCT_RZ = 1.0D14,
 RCT_RX = 2.0D3,
 RCT_RY = 1.0D14,
/

&PARAM_MKINIT_ADVECT
! ENV_U  = -35.D0,
! ENV_V  = -40.D0,
 ENV_U  = 40.D0,
 ENV_V  = 0.D0, 
 SHAPE_NC  =  '#{shape_nc}', 
 MAXMIN_NC =   1.D0,
/

EOS
  f.close
  
end

def gen_run_conf(conf_name, nprocx, nprocy, imax, kmax, dx, dz, dtsec_dyn, flxEvalType, dataDir)

  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-LES run configulation
#
#####

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
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
 TIME_DT                    = #{dtsec_dyn}, !#{TIME_DT_SEC},
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = #{dtsec_dyn}, 
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .false.,
 STATISTICS_use_globalcomm = .true.,
/

&PARAM_TRACER
 TRACER_TYPE = 'SN14',
/

&PARAM_ATMOS_HYDROSTATIC
 HYDROSTATIC_uselapserate = .true.,
/

&PARAM_ATMOS
 ATMOS_DYN_TYPE    = "FDM-HEVE",
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_IN_BASENAME      = "init_00000000000.000",
 ATMOS_RESTART_OUTPUT           = .false.,
 ATMOS_VARS_CHECKRANGE          = .true.,
/

&PARAM_ATMOS_REFSTATE
! ATMOS_REFSTATE_TYPE       = "INIT",
! ATMOS_REFSTATE_TYPE       = "UNIFORM",
  ATMOS_REFSTATE_POTT_UNIFORM = 300.0D0, 
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE       = "CONST",
 ATMOS_BOUNDARY_USE_VELZ   = .true.,
 ATMOS_BOUNDARY_VALUE_VELZ =  0.D0,
 ATMOS_BOUNDARY_TAUZ       = 10.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF_COEF  = 0.D0,
 ATMOS_DYN_DIVDMP_COEF   = 0.D0,
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
&HISTITEM item='NC'   /
!&HISTITEM item='DENS'   /
!&HISTITEM item='MOMX'   /
!&HISTITEM item='MOMY'   /
!&HISTITEM item='MOMZ'   /
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
  CONF_GEN_CASE_HASH_LIST.each{|case_hash|
    CONF_GEN_NUMERIC_HASHLIST.each{|numeric_hash|
      dataDir = "./#{resol_hash["TAG"]}/#{case_hash["TAG"]}/#{numeric_hash["TAG"]}/"

      puts "generate init.conf and run.conf (Dir=#{dataDir})"
      
      init_conf_name = "#{dataDir}init.conf" 
      gen_init_conf(init_conf_name, \
                    resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["KMAX"], \
                    resol_hash["DX"], resol_hash["DZ"], case_hash["SHAPE_NC"] )
      run_conf_name = "#{dataDir}run.conf"
      gen_run_conf(run_conf_name, \
                   resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["KMAX"], \
                   resol_hash["DX"], resol_hash["DZ"], resol_hash["DTDYN"], numeric_hash["TAG"], dataDir )
    }
  }
}
