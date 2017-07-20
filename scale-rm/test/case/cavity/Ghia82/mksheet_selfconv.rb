require "numru/gphys"
require "csv"

include NumRu

###########
#VAR_NAME   = "PT"
#VAR_CUTHASH = {}
#REFSOL_DIR = "025m/CTRL/FVM_UD5"
#RESOL_LIST = ["400m", "200m", "100m", "050m"]
#TIME_LIST = [0.0, 300.0, 600.0, 900.0]
############
#VAR_NAME   = "PT"
#VAR_CUTHASH = {}
#REFSOL_DIR = "031m/CTRL/FVM_UD5"
#RESOL_LIST = ["500m", "250m", "125m", "063m"]
#TIME_LIST = [0.0, 250.0, 500.0, 750.0, 1000.0]
############
VAR_NAME   = "W"
VAR_CUTHASH = {}
REFSOL_DIR = "0125m/CTRL/FVM_UD5"
RESOL_LIST = ["1000m", "500m", "250m", "125m"]
TIME_LIST = [0.0, 0.5*86400.0, 1.0*86400.0, 3.0*86400.0, 6.0*86400.0, 10.0*86400.0]
############


CASE_LIST  = ["CTRL"]
FLXSCHEME_LIST = ["UD1", "UD3", "UD5", "CD2", "CD4", "CD6" ]


def get_gphys(ncFileName, varname, cutHash=nil)
  gphys = GPhys::IO.open(ncFileName, varname)
#  p cutHash
  gphys = gphys.cut(cutHash) if cutHash != nil
  return gphys
end

l2error_table = NArray.sfloat(CASE_LIST.length, RESOL_LIST.length, FLXSCHEME_LIST.length, TIME_LIST.length)
linferror_table = NArray.sfloat(CASE_LIST.length, RESOL_LIST.length, FLXSCHEME_LIST.length,TIME_LIST.length)

CASE_LIST.each_with_index{|expcase,c|
  RESOL_LIST.each_with_index{|resol,r|
    FLXSCHEME_LIST.each_with_index{|scheme,s|
      target_dir = "#{resol}/#{expcase}/FVM_#{scheme}"
      p "target dir=#{target_dir}.."

      gp_l2error = get_gphys("#{target_dir}/error_norm.nc", "l2error")
      gp_lInferror = get_gphys("#{target_dir}/error_norm.nc", "linferror")

      TIME_LIST.each_with_index{|time,t|
        l2error_table[c,r,s,t] = gp_l2error.cut("time"=>time).val[0].to_f
        linferror_table[c,r,s,t] = gp_lInferror.cut("time"=>time).val[0].to_f
      }
      
    }
  }

  TIME_LIST.each_with_index{|time,t|
    csv_header = ["# resol #{FLXSCHEME_LIST.join(" ")}" ]
    CSV.open("selfconv/l2error_t#{time.to_i}.csv", "wb", :headers =>csv_header, :write_headers =>true) do |csv|
      RESOL_LIST.each_with_index{|resol,r|
        ary = [resol].push( l2error_table[c,r,true,t].to_a ).flatten
        csv << ary.map(&:to_s)
      }
    end
  TIME_LIST.each_with_index{|time,t|
    csv_header = ["# resol #{FLXSCHEME_LIST.join(" ")}" ]
    CSV.open("selfconv/linferror_t#{time.to_i}.csv", "wb", :headers =>csv_header, :write_headers =>true) do |csv|
      RESOL_LIST.each_with_index{|resol,r|
        ary = [resol].push( linferror_table[c,r,true,t].to_a ).flatten
        csv << ary.map(&:to_s)
      }
    end
  }  

  }  
}
