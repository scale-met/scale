require "numru/gphys"
include NumRu


file1 = ARGV.shift
file2 = ARGV.shift


f1 = NetCDF.open(file1)
f2 = NetCDF.open(file2)

vars = GPhys::IO.var_names(f1)

vars.each do |var|
  v1 = GPhys::IO.open(f1, var)
  v2 = GPhys::IO.open(f2, var)

  unless v1.val == v2.val
    p "#{var} is different"
  end
end

f1.close
f2.close
