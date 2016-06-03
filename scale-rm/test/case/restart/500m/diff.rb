require "numru/gphys"
include NumRu


file1 = ARGV.shift
file2 = ARGV.shift


f1 = NetCDF.open(file1)
f2 = NetCDF.open(file2)

vars = GPhys::IO.var_names(f1)

flag = true
vars.each do |var|
  v1 = GPhys::IO.open(f1, var).val
  v2 = GPhys::IO.open(f2, var).val

  unless v1 == v2
    p v1 - v2
    p "#{var} is different"
    flag = false
  end
end

f1.close
f2.close

if flag
  p "OK"
else
  exit(-1)
end
