#!/bin/env ruby

require 'fileutils'

TIME_DURATION_SEC       = "600.D0"
HISTORY_TINTERVAL_SEC   = "50.D0"
CONF_GEN_RESOL_HASHLIST = \
[ \
  { "TAG"=>"500m", "DX"=>500E0, "DZ"=>500.0E0, 
    "KMAX"=>4, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>1.0, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"250m", "DX"=>250E0, "DZ"=>250.0E0, 
    "KMAX"=>4, "IMAX"=>80, "JMAX"=>3, "DTDYN"=>0.5, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"125m", "DX"=>125E0, "DZ"=>125E0, 
    "KMAX"=>4, "IMAX"=>160, "JMAX"=>3, "DTDYN"=>0.25, "NPRCX"=> 1, "NPRCY"=>1}, \
  { "TAG"=>"063m", "DX"=>62.5E0, "DZ"=>62.5E0, 
    "KMAX"=>4, "IMAX"=>80, "JMAX"=>3, "DTDYN"=>0.125, "NPRCX"=> 4, "NPRCY"=>1} \
]
CONF_GEN_CASE_HASH_LIST = \
[ \
  {"TAG"=>"COS", "SHAPE_NC"=>"COS"},       \
  {"TAG"=>"RECT", "SHAPE_NC"=>"RECT"},     \
  {"TAG"=>"COSBELL", "SHAPE_NC"=>"BUBBLE"} \
]
CONF_GEN_NUMERIC_HASHLIST = \
[ \
  {"TAG"=>"FVM_CD2"}, {"TAG"=>"FVM_CD4"}, {"TAG"=>"FVM_CD6"},  \
  {"TAG"=>"FVM_UD1"}, {"TAG"=>"FVM_UD3"}, {"TAG"=>"FVM_UD5"},  \
]

#########################################################

def gen_init_conf(conf_name, nprocx, nprocy, imax, kmax, dx, dz, shape_nc)
  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-RM mkinit configulation for linear advection test (1D)
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
 IMAX = #{imax}, IHALO = 3, 
 JMAX = 3, JHALO = 3,
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

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init",
/

&PARAM_CONST
 CONST_GRAV      =   0.00000000000000     ,
/

&PARAM_MKINIT
 MKINIT_initname = "TRACERBUBBLE",
/

&PARAM_BUBBLE
 BBL_CZ = 10.0D3,
 BBL_CX = 10.0D3,
 BBL_CY = 12.0D3,
 BBL_RZ = 1.0D14,
 BBL_RX = 3.0D3,
 BBL_RY = 1.0D14,
/

&PARAM_RECT
 RCT_CZ = 10.0D3,
 RCT_CX = 10.0D3,
 RCT_CY = 12.0D3,
 RCT_RZ = 1.0D14,
 RCT_RX = 3.0D3,
 RCT_RY = 1.0D14,
/

&PARAM_MKINIT_TRACERBUBBLE
 ENV_U     = 40.D0,
 ENV_V     = 0.D0, 
 SHAPE_NC  = '#{shape_nc}', 
 BBL_NC    = 1.D0,
/

EOS
  f.close
  
end

def gen_run_conf( conf_name, 
      nprocx, nprocy,
      imax, kmax, dx, dz, dtsec_dyn,
      shape_nc, flxEvalType, fctFlag, dataDir)
1
  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-RM run configulation for linear advection test (1D)
#
#####

&PARAM_PRC
 PRC_NUM_X       = #{nprocx},  
 PRC_NUM_Y       = #{nprocy},
/

&PARAM_INDEX
 KMAX = #{kmax}, 
 IMAX = #{imax}, IHALO = 3, 
 JMAX = 3, JHALO = 3,
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
 TIME_DT                    = #{dtsec_dyn}, 
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = #{dtsec_dyn}, 
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .false.,
 STATISTICS_use_globalcomm = .true.,
/

&PARAM_CONST
 CONST_GRAV      =   0.00000000000000     ,
/
&PARAM_ATMOS_HYDROSTATIC
 HYDROSTATIC_uselapserate = .true.,
/

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
 ATMOS_BOUNDARY_TYPE       = "CONST",
 ATMOS_BOUNDARY_USE_VELZ   = .true.,
 ATMOS_BOUNDARY_VALUE_VELZ =  0.D0,
 ATMOS_BOUNDARY_TAUZ       = 10.D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_TINTEG_LARGE_TYPE = "EULER",
 ATMOS_DYN_TINTEG_SHORT_TYPE = "RK4",
 ATMOS_DYN_TINTEG_TRACER_TYPE = "RK3WS2002",
 ATMOS_DYN_FVM_FLUX_TYPE        = "#{flxEvalType}",             
 ATMOS_DYN_FVM_FLUX_TRACER_TYPE = "#{flxEvalType}", 
 ATMOS_DYN_NUMERICAL_DIFF_COEF  = 0.D0,
 ATMOS_DYN_DIVDMP_COEF          = 0.D0,
 ATMOS_DYN_FLAG_FCT_TRACER      = #{fctFlag}, 
/

&PARAM_USER
 USER_do = .true., 
 InitShape = "#{shape_nc}"
/


&PARAM_HISTORY
 HISTORY_DEFAULT_BASENAME  = "history",
 HISTORY_DEFAULT_TINTERVAL = #{HISTORY_TINTERVAL_SEC},
 HISTORY_DEFAULT_TUNIT     = "SEC",
 HISTORY_DEFAULT_TAVERAGE  = .false.,
 HISTORY_DEFAULT_DATATYPE  = "REAL8",
 HISTORY_OUTPUT_STEP0      = .true.,
/

&HISTITEM item='DENS'    /
&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
&HISTITEM item='NC'   /
&HISTITEM item='l2error'   /
&HISTITEM item='linferror'   /


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
      ["F", "T"].each{|fct_flag|

        dataDir = "./#{resol_hash["TAG"]}/#{case_hash["TAG"]}/"
        dataDir += fct_flag=="T" ? "#{numeric_hash["TAG"]}_FCT/" : "#{numeric_hash["TAG"]}/"

        puts "Generate init.conf and run.conf (Dir=#{dataDir}) .."
        if !File.exists?(dataDir) then
          puts "Create directory .."
          FileUtils.mkdir_p(dataDir)
        end

        #-------------------------------------------------------
        init_conf_name = "#{dataDir}init.conf"
        
        init_shape = case_hash["SHAPE_NC"]
        init_shape = "BUBBLE" if init_shape == "COS"
        
        gen_init_conf(init_conf_name, 
                      resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["KMAX"], 
                      resol_hash["DX"], resol_hash["DZ"], init_shape )

        #-------------------------------------------------------
        run_conf_name = "#{dataDir}run.conf"
        gen_run_conf(run_conf_name, 
                     resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["KMAX"], 
                     resol_hash["DX"], resol_hash["DZ"], resol_hash["DTDYN"], case_hash["SHAPE_NC"], 
                     numeric_hash["TAG"].sub("FVM_",""), fct_flag, dataDir )
      }
    }
  }
}
