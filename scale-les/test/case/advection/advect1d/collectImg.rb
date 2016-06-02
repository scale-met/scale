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

xz_NC_t100 = VarNoAnimFig.new("", "NC", "time=100,y=0,z=0",  "0.1", "0:1", "x", "t100","") 
xz_NC_t300 = VarNoAnimFig.new("", "NC", "time=300,y=0,z=0",  "0.1", "0:1", "x", "t300","") 
xz_NC_t500 = VarNoAnimFig.new("", "NC", "time=500,y=0,z=0",  "0.1", "0:1", "x", "t500","")

exps = []
resols = ["125m","250m", "500m"]
shapes = ["COSBELL"]
schemes = ["UD1", "CD2", "UD3", "CD4", "UD5", "CD6"]
resols.each{|resol|
  shapes.each{|shape|
    schemes.each{|scheme|
      exps.push(Exp.new("exp_#{resol}_#{shape}_FDM_#{scheme}",
                        "#{DATASTORAGE_DIR}/#{resol}/#{shape}/FDM_#{scheme}/", ""))
    }
  }
}

p exps

merge_vars = [ "NC" ]


noAnimFigVars = []
exp_noAnimFigs =  [ \
                    xz_NC_t100, xz_NC_t300, xz_NC_t500 \
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

extra_exps.push(Exp.new("exp_FDMComp_COSBELL",
                        "#{FIGSTORAGE_DIR}/exp_FDMComp/COSBELL/", ""))


####################################

# Register figure objects for each experiment
@expsHash = {}
@expsWithMergeDatHash = {}

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

["500m","250m","125m"].each{|resol|
  x_NCCompari_FDMUDComp =  NoAnimOverplotFig.new("NC_#{resol}_FDMUDCompari", \
                     "#{getExpDirPath("125m_COSBELL_FDM_UD5")}/NC.nc, \
                      #{getExpDirPath("#{resol}_COSBELL_FDM_UD1")}/NC.nc, \
                      #{getExpDirPath("#{resol}_COSBELL_FDM_UD3")}/NC.nc, \
                      #{getExpDirPath("#{resol}_COSBELL_FDM_UD5")}/NC.nc", \
                     "NC,NC,NC,NC", "time=500,y=0,z=0", "", "0:1", "--title 'Comparison FMDs(UD1,UD3,UD5)'") 
  x_NCCompari_FDMCDComp =  NoAnimOverplotFig.new("NC_#{resol}_FDMCDCompari", \
                     "#{getExpDirPath("125m_COSBELL_FDM_UD5")}/NC.nc, \
                      #{getExpDirPath("#{resol}_COSBELL_FDM_CD2")}/NC.nc, \
                      #{getExpDirPath("#{resol}_COSBELL_FDM_CD4")}/NC.nc, \
                      #{getExpDirPath("#{resol}_COSBELL_FDM_CD6")}/NC.nc", \
                     "NC,NC,NC,NC", "time=500,y=0,z=0", "", "0:1", "--title 'Comparison FMDs(CD2,CD4,CD6)'") 
  x_NCCompari_FDMUDComp.createFigure(getExpDirPath("FDMComp_COSBELL"))
  x_NCCompari_FDMCDComp.createFigure(getExpDirPath("FDMComp_COSBELL"))

}
extra_exps.each{|exp|
  exp.create_thumb(DCMODEL_THUM_ORIGIN, "jpg")
#  exp.create_thumb(DCMODEL_THUM_ORIGIN, "png")
}
exit
####################################################

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

#  puts "delete merged data.."
#  merge_vars.each{|merge_varname|
#    `rm "#{exp.dirPath}#{merge_varname}.nc"`
#  }  
}

##

