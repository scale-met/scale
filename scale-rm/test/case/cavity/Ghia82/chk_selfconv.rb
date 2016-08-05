require "numru/gphys"
include NumRu

############
#VAR_NAME   = "PT"
#VAR_CUTHASH = {}
#REFSOL_DIR = "025m/CTRL/FVM_UD5"
#RESOL_LIST = ["400m", "200m", "100m", "050m"]
############
#VAR_NAME   = "PT"
#VAR_CUTHASH = {}
#REFSOL_DIR = "031m/CTRL/FVM_UD5"
#RESOL_LIST = ["500m", "250m", "125m", "063m"]
############
#VAR_NAME   = "W"
#VAR_CUTHASH = {}
#REFSOL_DIR = "0125m/CTRL/FVM_UD5"
#RESOL_LIST = ["2000m", "1000m", "0500m", "0250m"]
#TIME_LIST = [0.0, 2.0*3600.0, 4.0*3600.0, 6.0*3600.0, 8.0*3600.0, 10.0*3600.0]
#############
VAR_NAME   = "W"
VAR_CUTHASH = {}
REFSOL_DIR = "025km/CTRL/FVM_UD5"
RESOL_LIST = ["400km", "200km", "100km", "050km"]
TIME_LIST = [0.0, 2.0*24.0*3600.0, 4.0*24.0*3600.0, 6.0*24.0*3600.0, 8.0*24.0*3600.0, 10.0*24.0*3600.0, 12.0*24.0*3600.0, 14.0*24.0*3600.0]
#############

CASE_LIST  = ["CTRL"]
FLXSCHEME_LIST = ["UD1", "UD3", "UD5", "CD2", "CD4", "CD6" ]
timeList = [0.0, 30.0]

GPhys::extrapolation = true

###########################################################

def get_gphys(ncFileDir, varname, cutHash=nil)
  files = /#{ncFileDir}\/history.pe(\d*).nc/
  gphys = GPhys::IO.open(files, varname)
#  p cutHash
  gphys = gphys.cut(cutHash) if cutHash != nil
  return gphys
end

def calc_ErrorNorm(gp, gp_ref, ofname)

  nx = gp.axis('x').pos.length
  ny = gp.axis('y').pos.length
  nz = gp.axis('z').pos.length
  nt = gp.axis('time').pos.length
  nt_ = gp_ref.axis('time').pos.length

  grid = Grid.new(gp.axis('time'))
  gp_l2error = GPhys.new( grid, VArray.new(NArray.sfloat(nt), {"units"=>"1"}, "l2error") )
  gp_linferror = GPhys.new( grid, VArray.new(NArray.sfloat(nt), {"units"=>"1"}, "linferror") )
  
  tidx = 0
  gp_l2error[true] = 1.0
  
  GPhys.each_along_dims(
    [gp, gp_ref, gp_l2error, gp_linferror], "time"){
    |gp_, gp_ref_,l2error,linferror|

    #p gp_ref[true,true,true,tidx..tidx].sum
    #    p "tidx=#{tidx}, linf=" ##{gp_linferror[tidx].to_f}"
    #    gp_l2error[tidx] = Math::sqrt( gp_ref_.val[true,true,true,0..0].sum.to_f )
    p "tidx=#{tidx}"
    gp__ = gp_.val[false,0]
    gp_ref__ = gp_ref_.val[false,0]
    p gp_ref__
    l2error[0] = Math.sqrt( (((gp__ - gp_ref__)**2).sum  / (gp_ref__**2).sum).to_f )
    linferror[0] = ( (gp__ - gp_ref__).abs.max / gp_ref__.abs.max ).to_f 
    
    tidx = tidx + 1
    break if tidx == nt
    
#    p "l2error=#{ gp_l2error.val[tidx].to_f}"
#    p gp_linferror[tidx]

  }

  p "output.."
  ofile = NetCDF::create(ofname)
  GPhys::NetCDF_IO.write(ofile, gp_l2error)
  GPhys::NetCDF_IO.write(ofile, gp_linferror)
  ofile.close
end

RESOL_LIST.each{|resol|
  CASE_LIST.each{|expcase|
    FLXSCHEME_LIST.each{|scheme|
      target_dir = "#{resol}/#{expcase}/FVM_#{scheme}"
      p "target dir=#{target_dir}.."
      
      gp = get_gphys("#{target_dir}", VAR_NAME, VAR_CUTHASH)
      gp_ref = get_gphys("#{REFSOL_DIR}", VAR_NAME, VAR_CUTHASH)

      ofile = NetCDF::create("#{target_dir}/refsol_diff.nc")

      p gp_ref.shape
      p gp.shape
      
      x = gp.axis("x")
      y = gp.axis("y")
      z = gp.axis("z")
      t = gp.axis("time")
      nt = t.pos.length

      tidx = 0
      GPhys::IO.each_along_dims_write(
        [gp], ofile, "time"){
        |gp_|

        p "tidx=#{tidx}"
        break if tidx == nt

        gp_ref_  = gp_ref[false,tidx..tidx].interpolate(x.pos, y.pos, z.pos).rename(VAR_NAME+"_ref")
        
        gp_diff = (gp_ - gp_ref_).rename(VAR_NAME+"_diff")
        tidx = tidx + 1
        [gp_ref_, gp_diff]
      }

      ofile.close

      gp_ref = GPhys::IO.open("#{target_dir}/refsol_diff.nc", VAR_NAME+"_ref")
      calc_ErrorNorm(gp, gp_ref, "#{target_dir}/error_norm.nc")
    }
  }
}
