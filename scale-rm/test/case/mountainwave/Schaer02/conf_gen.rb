#!/bin/env ruby

TIME_DT_SEC             = "3.0D0"
TIME_DURATION_SEC       = "36000.D0"
HISTORY_TINTERVAL_SEC   = "1800.D0"
CONF_GEN_RESOL_HASHLIST = \
[ \
  { "TAG"=>"2000m", "DX"=>2000E0, "DZ"=>300.0E0, 
    "KMAX"=>65, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>0.5E0, "NPRCX"=> 5, "NPRCY"=>1}, \
  { "TAG"=>"1000m", "DX"=>1000E0, "DZ"=>300.0E0,
    "KMAX"=>65, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>0.5E0, "NPRCX"=> 5, "NPRCY"=>1}, \
  { "TAG"=>"0500m", "DX"=>500E0, "DZ"=>300.0E0,
    "KMAX"=>65, "IMAX"=>40, "JMAX"=>3, "DTDYN"=>0.5E0, "NPRCX"=> 10, "NPRCY"=>1}, \
  { "TAG"=>"0250m", "DX"=>250E0, "DZ"=>300.0E0,
    "KMAX"=>65, "IMAX"=>80, "JMAX"=>3, "DTDYN"=>0.3E0, "NPRCX"=> 10, "NPRCY"=>1}, \
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
 BUFFER_DZ = 9000.D0,  
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

&PARAM_TOPO
 TOPO_OUT_BASENAME = "topo",
/

&PARAM_TRACER
 TRACER_TYPE = 'DRY',
/

&PARAM_ATMOS_VARS
 ATMOS_RESTART_OUTPUT         = .true.,
 ATMOS_RESTART_OUT_BASENAME   = "init",
/

&PARAM_MKTOPO
 MKTOPO_name   = "SCHAER",
/

&PARAM_MKTOPO_SCHAER
 SCHAER_CX     = 50.D3,
 SCHAER_RX     =  5.D3,
 SCHAER_LAMBDA =  4.D3,
 SCHAER_HEIGHT = 250D0,
/

&PARAM_MKINIT
 MKINIT_initname = "MOUNTAINWAVE",
/

&PARAM_MKINIT_MOUNTAINWAVE
 ENV_U        = 10.D0,
 ENV_V        =  0.D0,
 SCORER       = 1.D-3,
/
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
 BUFFER_DZ = 9000.D0,  
 BUFFFACT  =    1.D0,
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

&PARAM_TOPO
 TOPO_IN_BASENAME = "topo",
/

&PARAM_TRACER
 TRACER_TYPE = 'DRY',
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
 ATMOS_REFSTATE_TYPE       = "INIT",
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE       = "CONST",
 ATMOS_BOUNDARY_USE_VELZ   = .true.,
 ATMOS_BOUNDARY_VALUE_VELZ =  0.D0,
 ATMOS_BOUNDARY_TAUZ       = 10.D0,
 ATMOS_BOUNDARY_USE_VELX   = .true.,
 ATMOS_BOUNDARY_VALUE_VELX =  10.D0,
 ATMOS_BOUNDARY_TAUX       = 10.D0,
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
 HISTORY_DEFAULT_ZINTERP   = "HGT",
 HISTORY_OUTPUT_STEP0      = .true.,
/

&HISTITEM item='DENS' /
&HISTITEM item='U'    /
&HISTITEM item='V'    /
&HISTITEM item='W'    /
&HISTITEM item='PT'   /
&HISTITEM item='T'    /
&HISTITEM item='PRES' /

&PARAM_MONITOR
 MONITOR_STEP_INTERVAL = 360,
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
