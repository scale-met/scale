#!/usr/bin/env ruby
# coding: utf-8

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Generate init.conf and run.conf
#
# These configuration files are used for test cases of riging warm bubble, following
# experimental setup in Wicker and Skmarock (1998).
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TIME_DT_SEC             = "5.0D0"
TIME_DURATION_SEC       = "1200.D0"
HISTORY_TINTERVAL_SEC   = "100.D0"

#------------------------------------------------------------------------------------

CONF_GEN_RESOL_HASHLIST = \
[ \
  { "TAG"=>"500m", "DX"=>500E0, "DZ"=>500.0E0,
    "KMAX"=>20, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>0.5E0, "NPRCX"=> 2, "NPRCY"=>1}, \
  { "TAG"=>"250m", "DX"=>250E0, "DZ"=>250.0E0,
    "KMAX"=>40, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>0.25E0, "NPRCX"=> 4, "NPRCY"=>1}, \
  { "TAG"=>"125m", "DX"=>125E0, "DZ"=>125.0E0,
    "KMAX"=>80, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>0.125E0, "NPRCX"=> 8, "NPRCY"=>1}, \
  { "TAG"=>"063m", "DX"=>62.5E0, "DZ"=>62.5E0,
    "KMAX"=>160, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>0.0625E0, "NPRCX"=> 16, "NPRCY"=>1}, \
  { "TAG"=>"031m", "DX"=>31.25E0, "DZ"=>31.25E0,
    "KMAX"=>320, "IMAX"=>20, "JMAX"=>3, "DTDYN"=>0.03125E0, "NPRCX"=> 32, "NPRCY"=>1}, \
  #  { "TAG"=>"050m", "DX"=>50E0, "DZ"=>50.0E0,
#    "KMAX"=>128, "IMAX"=>64, "JMAX"=>3, "DTDYN"=>0.0625E0, "NPRCX"=>16, "NPRCY"=>1}, \
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

require 'fileutils'

def gen_init_conf( conf_name,
                   nprocx, nprocy, imax, jmax, kmax, dx, dy, dz )

  f = File.open(conf_name, "w")
  f.print <<EOS
#####
#
# SCALE-RM mkinit configulation for rising warmbubble test
# (following experimental setup in Wicker and Skmarock (1998))
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
 JMAX = #{jmax}, JHALO = 3,
/

&PARAM_GRID
 DZ =  #{dz},
 DX =  #{dx},
 DY =  #{dy},
 BUFFER_DZ =   0.D0,
 BUFFFACT  =   1.D0,
 GRID_OFFSET_X = -10.D3,
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

&PARAM_MKINIT
 MKINIT_initname = "WARMBUBBLE",
/

&PARAM_BUBBLE
 BBL_CZ =  2.0D3,
 BBL_CX =  0.0D3,
 BBL_CY =  0.0D3,
 BBL_RZ =  2.0D3,
 BBL_RX =  2.0D3,
 BBL_RY =  5.0D13,
/

&PARAM_MKINIT_WARMBUBBLE
 SFC_THETA   = 300.D0,
 ENV_RH      =   0.D0,
 ENV_L1_ZTOP = 10.0D3
 BBL_THETA   =   2.D0,
/

EOS
  f.close

end

def gen_run_conf( conf_name,
                  nprocx, nprocy,
                  imax, jmax, kmax, dx, dy, dz, dtsec_dyn,
                  flxEvalType, fctFlag, dataDir )

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
 IMAX = #{imax}, IHALO = 3,
 JMAX = #{jmax}, JHALO = 3,
/

&PARAM_GRID
 DZ =  #{dz},
 DX =  #{dx},
 DY =  #{dy},
 BUFFER_DZ =   0.D0,
 BUFFFACT  =   1.D0,
 GRID_OFFSET_X = -10.D3,
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
! ATMOS_REFSTATE_TYPE       = "INIT",
! ATMOS_REFSTATE_TYPE       = "UNIFORM",
  ATMOS_REFSTATE_POTT_UNIFORM = 300.0D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_TINTEG_LARGE_TYPE    = "EULER",
 ATMOS_DYN_TINTEG_SHORT_TYPE    = "RK3WS2002",
 ATMOS_DYN_TINTEG_TRACER_TYPE   = "RK3WS2002",
 ATMOS_DYN_FVM_FLUX_TYPE        = "#{flxEvalType}",
 ATMOS_DYN_FVM_FLUX_TRACER_TYPE = "#{flxEvalType}",
 ATMOS_DYN_NUMERICAL_DIFF_COEF  = 0.D0,
 ATMOS_DYN_WDAMP_TAU            = 10.D0,
 ATMOS_DYN_WDAMP_HEIGHT         = 15.D3,
 ATMOS_DYN_DIVDMP_COEF          = 0.1D0,
 ATMOS_DYN_FLAG_FCT_TRACER      = #{fctFlag},
/

!&PARAM_USER
! USER_do = .true.,
!/

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

        init_conf_name = "#{dataDir}init.conf"
        gen_init_conf(init_conf_name,
                      resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["JMAX"], resol_hash["KMAX"],
                      resol_hash["DX"], resol_hash["DX"], resol_hash["DZ"]  )

        run_conf_name = "#{dataDir}run.conf"
        gen_run_conf(run_conf_name,
                     resol_hash["NPRCX"], resol_hash["NPRCY"], resol_hash["IMAX"], resol_hash["JMAX"], resol_hash["KMAX"],
                     resol_hash["DX"], resol_hash["DX"], resol_hash["DZ"], resol_hash["DTDYN"],
                     numeric_hash["TAG"].sub("FVM_",""), fct_flag, dataDir )
      }
    }
  }
}

