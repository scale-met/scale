#!/usr/bin/env ruby
# -*- coding: euc-jp -*-
#
#     

##########################
## Configuration
#########################

# Configuration for creating thumbnail
DCMODEL_THUM_ORIGIN = "/home/ykawai/dcmodel-thum.rb"

# End day of time integration
LAST_DAY= 3300000  #[day]
#LAST_DAY=730000  #[day]
LASTDATE_SUFFIX="_10000yr"
YT_SUFFIX="_0-10000yr"

# Path of parent directory where the simulation results are saved.
DATASTORAGE_DIR="./"

# Path of parent directory where the figures of simulation results are saved.
FIGSTORAGE_DIR="./figs/"

#################################################################

require "numru/gphys"
require "/home/ykawai/collectImgLib.rb"
include NumRu
include CollectImgLib


VarNoAnimFig = Struct.new("VarNoAnimFig", :ncFileSuffix, :name, :cutpos, :intrv, :range, :prefix, :suffix, :gpOpts)
VarAnimFig = Struct.new("VarAnimFig", :name, :cutpos, :intrv, :range, :animTimeInfo, :prefix, :suffix, :gpOpts)

xz_PT_t03600 = VarNoAnimFig.new("", "PT", "time=03600,y=0",  "0.2", "299:303", "xz", "t03600","") 
xz_PT_t07200 = VarNoAnimFig.new("", "PT", "time=07200,y=0",  "0.2", "299:303", "xz", "t07200","") 
xz_PT_t18000 = VarNoAnimFig.new("", "PT", "time=18000,y=0",  "0.2", "299:303", "xz", "t18000","") 
xz_PT_t36000 = VarNoAnimFig.new("", "PT", "time=36000,y=0",  "0.2", "299:303", "xz", "t36000","")

xz_U_t03600_gl = VarNoAnimFig.new("", "U", "time=03600,y=0",  "0.2", "8:12", "xz", "t03600_gl","") 
xz_U_t07200_gl = VarNoAnimFig.new("", "U", "time=07200,y=0",  "0.2", "8:12", "xz", "t07200_gl","") 
xz_U_t18000_gl = VarNoAnimFig.new("", "U", "time=18000,y=0",  "0.2", "8:12", "xz", "t18000_gl","") 
xz_U_t36000_gl = VarNoAnimFig.new("", "U", "time=36000,y=0",  "0.2", "8:12", "xz", "t36000_gl","")

xz_U_t03600 = VarNoAnimFig.new("", "U", "time=03600,y=0,x=40e3:60e3,z=0:10e3",  "0.2", "8:12", "xz", "t03600","") 
xz_U_t07200 = VarNoAnimFig.new("", "U", "time=07200,y=0,x=40e3:60e3,z=0:10e3",  "0.2", "8:12", "xz", "t07200","") 
xz_U_t18000 = VarNoAnimFig.new("", "U", "time=18000,y=0,x=40e3:60e3,z=0:10e3",  "0.2", "8:12", "xz", "t18000","") 
xz_U_t36000 = VarNoAnimFig.new("", "U", "time=36000,y=0,x=40e3:60e3,z=0:10e3",  "0.2", "8:12", "xz", "t36000","")

xz_W_t03600_gl = VarNoAnimFig.new("", "W", "time=03600,y=0,x=40e3:60e3,z=0:10e3",  "0.05", "-1:1", "xz", "t03600","") 
xz_W_t07200_gl = VarNoAnimFig.new("", "W", "time=07200,y=0,x=40e3:60e3,z=0:10e3",  "0.05", "-1:1", "xz", "t07200","") 
xz_W_t18000_gl = VarNoAnimFig.new("", "W", "time=18000,y=0,x=40e3:60e3,z=0:10e3",  "0.05", "-1:1", "xz", "t18000","") 
xz_W_t36000_gl = VarNoAnimFig.new("", "W", "time=36000,y=0,x=40e3:60e3,z=0:10e3",  "0.05", "-1:1", "xz", "t36000","")

xz_W_t03600 = VarNoAnimFig.new("", "W", "time=03600,y=0",  "0.05", "-1:1", "xz", "t03600_gl","") 
xz_W_t07200 = VarNoAnimFig.new("", "W", "time=07200,y=0",  "0.05", "-1:1", "xz", "t07200_gl","") 
xz_W_t18000 = VarNoAnimFig.new("", "W", "time=18000,y=0",  "0.05", "-1:1", "xz", "t18000_gl","") 
xz_W_t36000 = VarNoAnimFig.new("", "W", "time=36000,y=0",  "0.05", "-1:1", "xz", "t36000_gl","")

exps = []
resols = ["0250m"]#"1000m","0500m","0250m"] #"125m"]#["125m", "500m", "1000m", "2000m"]
schemes = ["UD5"]#, "UD5"]#"CD2", "UD3", "CD4", "UD5", "CD6"]
resols.each{|resol|
  schemes.each{|scheme|
    exps.push(Exp.new("exp_#{resol}_FDM_#{scheme}",
                      "#{DATASTORAGE_DIR}/#{resol}/FDM_#{scheme}/", ""))
  }
}

p exps

merge_vars = [ "PT", "U", "W" ]


noAnimFigVars = []
exp_noAnimFigs =  [ \
#                    xz_PT_t03600, xz_PT_t07200, xz_PT_t18000, xz_PT_t36000, \
                    xz_U_t03600, xz_U_t07200, xz_U_t18000, xz_U_t36000, \
                    xz_W_t03600, xz_W_t07200, xz_W_t18000, xz_W_t36000, \
                    xz_U_t03600_gl, xz_U_t07200_gl, xz_U_t18000_gl, xz_U_t36000_gl, \
                    xz_W_t03600_gl, xz_W_t07200_gl, xz_W_t18000_gl, xz_W_t36000_gl, \
                  ]

# 
animTimeInfo = AnimTimeInfo.new(0, LAST_DAY, 1500.0, "day")
animFigFlag = false

animFigVars = []
exp_AnimFigs =  [ 
#                 VarAnimFig.new("U", "lon=0", "0.05", "-1.2:1.2", animTimeInfo, "yz", "anim", ""),  
#                 VarAnimFig.new("PTemp", "lon=0,sig=-1:0", "2", "270:310", animTimeInfo, "yz", "anim", ""),
#                 VarAnimFig.new("Salt", "lon=0,sig=-1:0", "0.2", "33:37", animTimeInfo, "yz", "anim", ""),
#                 VarAnimFig.new("MassStreamFunc", "lat=-90:90,sig=-1:0", "10", "-100:100", animTimeInfo, "yz", "anim", "") 
                ]




#################################

extra_exps = [ 
#              Exp.new("exp_comm", "./common/"), 
#              Exp.new("exp_EOSComp_CARediGM", "./EOSComp/CARediGM/"), 
             ]


####################################

# Register figure objects for each experiment
@expsHash = {}

exps.each_with_index{|exp, i|
  @expsHash[exp.name.gsub("exp_", "")] = exp
  noAnimFigVars[i] = exp_noAnimFigs
  animFigVars[i] = exp_AnimFigs
}

extra_exps.each_with_index{|exp, i|
  @expsHash[exp.name.gsub("exp_", "")] = exp
}


####################################

def getExpDirPath(expHashKey)
p expHashKey
  return @expsHash[expHashKey].dirPath
end


#t_KECompari_EOSComp = NoAnimOverplotFig.new("KEAvg_EOSCompari", 
#                     "#{getExpDirPath("EOSQ_IDIFF_CA_GM")}/KEAvg.nc, \
#                     #{getExpDirPath("EOSL_IDIFF_CA_GM")}/KEAvg.nc, \
#                      #{getExpDirPath("EOSJM95_IDIFF_CA_GM")}/KEAvg.nc", 
#                     "KEAvg,KEAvg,KEAvg", "t=0:1e10", "", "0:0.05", "--title 'K.E.(global mean)'") 
#t_KECompari_EOSComp.createFigure(getExpDirPath("EOSComp_CARediGM"))

#t_KECompari_VDiffComp = NoAnimOverplotFig.new("KEAvg_VDiffCompari", 
#                     "#{getExpDirPath("EOSL_HDIFF")}/KEAvg.nc, \
#                      #{getExpDirPath("EOSL_HDIFF_VDIFF100")}/KEAvg.nc", 
#                     "KEAvg,KEAvg", "t=0:1e10", "", "0:0.05", "--title 'K.E.(global mean)'") 
#t_KECompari_VDiffComp.createFigure(getExpDirPath("VDiffComp_noCARediGM"))


####################################################

extra_exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
#  exp.create_thumb(DCMODEL_THUM_ORIGIN, "png")
}
#exit

exps.each_with_index{|exp, i|
  puts "= experiment: #{exp.name}"

  puts "merge splited data.."
  merge_vars.each{|merge_varname|
    puts "var=#{merge_varname}"
    files = /#{exp.dirPath}history.pe(\d*).nc/
    gp_tmp = GPhys::IO.open(files, merge_varname)
    
    ofile = NetCDF::create("#{exp.dirPath}#{merge_varname}.nc")
    GPhys::NetCDF_IO.write(ofile, gp_tmp)
    ofile.close
  }

  noAnimFigVars[i].each{|figvar|
      exp.add_NoAnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.prefix, figvar.suffix, figvar.gpOpts, "", figvar.ncFileSuffix) 
  }

  if animFigFlag then
    animFigVars[i].each{|figvar|
      exp.add_AnimFig(figvar.name, figvar.cutpos, figvar.intrv, figvar.range, figvar.animTimeInfo, figvar.prefix, figvar.suffix, figvar.gpOpts)
      p figvar.animTimeInfo 
    }
  end


  exp.create_allFig
#  exp.create_thumb(DCMODEL_THUM_ORIGIN, "png") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg") if noAnimFigVars[i].size > 0
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "gif") if animFigVars[i].size > 0 and animFigFlag

  figdirpath = "#{FIGSTORAGE_DIR}/#{exp.name}"
  p "Move created figures from '#{exp.dirPath}' to '#{figdirpath}'"
  if !File.exist?(figdirpath) then
    `mkdir -p #{figdirpath}`
  end
  FileUtils.mv(Dir.glob("#{exp.dirPath}/*.{png,jpg,gif}"), figdirpath)


#  config_nml = Dir.glob("#{exp.dirPath}/config*.nml")[0]
#  p "Move c '#{config_nml}' to '#{FIGSTORAGE_DIR}/#{exp.name}'"
#  FileUtils.cp(config_nml, "#{FIGSTORAGE_DIR}/#{exp.name}/config.nml")


  puts "delete splited data.."
  merge_vars.each{|merge_varname|
    `rm "#{exp.dirPath}#{merge_varname}.nc"`
  }

}

