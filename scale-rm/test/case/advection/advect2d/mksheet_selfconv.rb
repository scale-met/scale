#!/usr/bin/env ruby

require "numru/gphys"
require "fileutils"
require "csv"

include NumRu

################################################################################

CASE_LIST  = ["COSBELL","RECT"]
RESOL_LIST = ["500m", "250m", "125m", "063m"]
FLXSCHEME_LIST = [
  "UD1", "UD3", "UD5", "CD2", "CD4", "CD6",
  "UD3_FCT", "UD5_FCT", "CD2_FCT", "CD4_FCT", "CD6_FCT"  
]
TIME_LIST = [50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]

################################################################################

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

      gp_l2error = get_gphys("#{target_dir}/history.pe000000.nc", "l2error")
      gp_lInferror = get_gphys("#{target_dir}/history.pe000000.nc", "linferror")

      TIME_LIST.each_with_index{|time,t|
        l2error_table[c,r,s,t] = gp_l2error.cut("time"=>time).val[0].to_f
        linferror_table[c,r,s,t] = gp_lInferror.cut("time"=>time).val[0].to_f
      }
      
    }
  }
}

CASE_LIST.each_with_index{|expcase,c|
  dirpath = "#{Dir.pwd}/selfconv/#{expcase}"
  FileUtils.mkdir_p(dirpath) unless FileTest.exists?(dirpath)

  p "Create #{dirpath}/l2error_t*.csv*"  
  TIME_LIST.each_with_index{|time,t|
    csv_header = ["# resol #{FLXSCHEME_LIST.join(" ")}" ]
    CSV.open("#{dirpath}/l2error_t#{time.to_i}.csv", "wb", :headers =>csv_header, :write_headers =>true) do |csv|
      RESOL_LIST.each_with_index{|resol,r|
        ary = [resol].push( l2error_table[c,r,true,t].to_a ).flatten
        csv << ary.map(&:to_s)
      }
    end
  }
  p "Create #{dirpath}/linferror_t*.csv*"  
  TIME_LIST.each_with_index{|time,t|
    csv_header = ["# resol #{FLXSCHEME_LIST.join(" ")}" ]
    CSV.open("#{dirpath}/linferror_t#{time.to_i}.csv", "wb", :headers =>csv_header, :write_headers =>true) do |csv|
      RESOL_LIST.each_with_index{|resol,r|
        ary = [resol].push( linferror_table[c,r,true,t].to_a ).flatten
        csv << ary.map(&:to_s)
      }
    end
  }
  
}
