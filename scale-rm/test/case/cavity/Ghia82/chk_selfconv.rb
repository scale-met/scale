require "numru/gphys"
include NumRu

#############################################################

L     = 20e3
Ulid  = 5.0E01

VAR_NAME   = "V"
VAR_CUTHASH = {}

REFSOL_RESOL = "063m"
REFSOL_DIR   = "063m/CTRL/FVM_UD5"
RESOL_LIST = ["1000m", "500m", "250m", "125m"] #, "063m"]

NPRC_hash ={}
NPRC_hash["2000m"] = {"nprcx"=> 1,"nprcy"=>1}
NPRC_hash["1000m"] = {"nprcx"=> 1,"nprcy"=>1}
NPRC_hash["500m"] = {"nprcx"=>2,"nprcy"=>2}
NPRC_hash["250m"] = {"nprcx"=>4,"nprcy"=>4}
NPRC_hash["125m"] = {"nprcx"=>8,"nprcy"=>8}
NPRC_hash["063m"] = {"nprcx"=>8,"nprcy"=>8}

CASE_LIST  = ["CTRL"]
FLXSCHEME_LIST = ["UD1", "UD3", "UD5", "CD2", "CD4", "CD6" ]
#timeList = [ 0.0, 0.5*86400.0, 1.0*86400.0,
#             2.0*86400.0, 3.0*86400.0, 5.0*86400.0, 7.5*86400.0,
#             10.0*86400.0 ]

#GPhys::extrapolation = true

###########################################################

def get_gphys(ncFileDir, varname, nprcx=1, nprcy=1, cutHash=nil)
  files = Array.new(nprcx){Array.new(nprcy){""}}
  for j in 0..nprcy-1
    for i in 0..nprcx-1
      pe_id = format("%06d", nprcx*j + i) 
      files[i][j] = "#{ncFileDir}/history.pe#{pe_id}.nc"
    end
  end
  gphys = GPhys::IO.open(files, varname)
  gphys = gphys.cut(cutHash) if cutHash != nil
  p 'Finish getting gphys..'
  return gphys
end

def calc_ErrorNorm(gp, gp_ref, ofname)

  nt_ = gp_ref.axis('time').pos.length    
  nt = gp.axis('time').pos.length
  nx = gp.axis('x').pos.length
  ny = gp.axis('y').pos.length
  nz = gp.axis('z').pos.length

  grid = Grid.new(gp_ref.axis('time'))
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
    l2error[0] = Math.sqrt( (((gp__ - gp_ref__)**2).sum  / (Ulid**2 *nx*ny)).to_f )
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

nprc_ref = NPRC_hash[REFSOL_RESOL]
p "nprc_ref info:", nprc_ref
gp_ref = get_gphys("#{REFSOL_DIR}", VAR_NAME,
                   nprc_ref["nprcx"], nprc_ref["nprcy"], 
                   VAR_CUTHASH)

RESOL_LIST.each{|resol|
  CASE_LIST.each{|expcase|
    FLXSCHEME_LIST.each{|scheme|
      target_dir = "#{resol}/#{expcase}/FVM_#{scheme}"
      p "target dir=#{target_dir}.."

      nprc = NPRC_hash[resol]
      p "nprc info:", nprc
      gp = get_gphys("#{target_dir}", VAR_NAME,
                     nprc["nprcx"], nprc["nprcy"], 
                     VAR_CUTHASH)

      nt_ = gp_ref.axis('time').pos.length
      gp = gp.cut("time"=>  0..gp_ref.axis('time').pos.val[nt_-1])


      #-----------------------------------
      
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

        gp_ref_  = gp_ref[false,0..0,tidx..tidx].interpolate(x.pos, y.pos).rename(VAR_NAME+"_ref")
        p gp_ref_.shape
        gp_diff = (gp_ - gp_ref_).rename(VAR_NAME+"_diff")#(gp_[false,0..0,0..0] - gp_ref_).rename(VAR_NAME+"_diff")
        tidx = tidx + 1
        [gp_ref_, gp_diff]
      }
      ofile.close

      
      gp_ref_intp = GPhys::IO.open("#{target_dir}/refsol_diff.nc", VAR_NAME+"_ref")
      calc_ErrorNorm( gp.cut("x"=>0..L).cut("y"=>0..L),
                      gp_ref_intp.cut("x"=>0..L).cut("y"=>0..L),
                      "#{target_dir}/error_norm.nc")

    }
  }
}
