#!/usr/bin/env ruby
=begin
 Convert tool from JMA MSM model level data (NuSDaS format) to GrADS data

 Usage:
   ruby convert.rb NuSDaS_top_direcory

 Requires:
  Ruby-NuSDaS (http://ruby.gfd-dennou.org/products/ruby-nusdas/)

 Notes
  * The lowest and heighest level are not physical level in NHM.
    * The values at the heighest level are those of the parent model (GSM).
  * Momentums are at half level (i-1/2 etc), so the last grid data for full level data should be remove.
  * The levels for groud layer are not stored in NuSDaS data file, so the levels are hard coded.
  * Data is stored from north to south in the NuSDaS files.
  * Values are stored at the last validetime.
=end


require "numru/nusdas"
include NumRu
include NMath


# constants
PTref = 300.0 # reference potential temperature
R = 287.05
CP = R * 7.0 * 0.5
RovCP = R / CP
Grav = 9.80665
Gamma = 0.005 # lapse rate K/m

#debug = true
debug = false

#skip_data = true
skip_data = false

class NArray
  alias :__to_s__ :to_s
  def to_s
    # convert to the big endian
    self[true,-1..0,false].hton.__to_s__
  end
end


# file names
prefix = ARGV.shift || "."
fname_air = File.join( prefix, "fcst_mdlA.nus")
fname_land = File.join( prefix, "anal_land0.nus")
basename_out = "MANAL"
grads_namelist = "namelist.grads_boundary"


validtime = -1
zrange = 1..-1
if debug
  xrange = 470..-250
  yrange = 250..-300
  trange = 0..0
else
  xrange = 0..-2
  yrange = 0..-2
  trange = true
end
xrangeu  = xrange.first..(xrange.last+1)
yrangev  = yrange.first..(yrange.last+1)


# get metrics and reference pressure
file_air = NuSDaS.open(fname_air)
zs = file_air.var("ZSsrf")[xrange,yrange,validtime,0]
u = file_air.var("RU")
zrp, zrw, vctrans_p, vctrans_w, dvtrans_p, dvtrans_w = u.dim("z").val(:full)
dvtrans_p = dvtrans_p[zrange].newdim!(0,0)
zrp = zrp[zrange].newdim!(0,0)
g2 = 1.0 + zs * dvtrans_p
z = zs + g2 * zrp
sz = z.to_s


nlon, nlat, nlev = g2.shape

lon = file_air.var("FLONsrf")[xrange,yrange,validtime,0]
lonrad = lon * PI / 180.0
lat = file_air.var("FLATsrf")[xrange,yrange,validtime,0]
latrad = lat * PI / 180.0

meta = file_air.instance_variable_get(:@meta)
basepoint = meta[:basepoint]
basepoint2 = meta[:basepoint2]
standard = meta[:standard]
standard2 = meta[:standard2]
distance = meta[:distance]
lat0 = basepoint2[0]*PI/180.0
lon0 = basepoint2[1]*PI/180.0
lat1 = standard[0]*PI/180.0
lat2 = standard2[0]*PI/180.0
c = log(cos(lat1)/cos(lat2)) / log( tan(PI/4.0-lat1/2.0)/tan(PI/4.0-lat2/2.0) )
mapf = ( cos(latrad)/cos(lat0) )**(c-1.0) * ( (1.0+sin(lat1))/(1.0+sin(latrad)) )**c
mapf = mapf.reshape!(nlon,nlat,1).to_type(NArray::SFLOAT)

dlon = lonrad - lon0
mask = dlon > PI
dlon[mask] = dlon[mask] - PI*2.0
mask = dlon < -PI
dlon[mask] = dlon[mask] + PI*2.0
alpha = c * dlon
alpha = - alpha if latrad[0,0] < 0.0 # southern hemisphere
sin_alpha = sin(alpha).reshape!(nlon,nlat,1)
cos_alpha = cos(alpha).reshape!(nlon,nlat,1)

exer = file_air.var("PAIRF")[xrange,yrange,zrange,validtime,0].to_type(NArray::SFLOAT)
P00 = meta[:subc]["ZHYB"][:presrf]
pbar = P00 * exer**(CP/R)
pbar = pbar.to_type(NArray::SFLOAT)


# open variables
#  atmospheric data
gru     = file_air.var("RU")
grv     = file_air.var("RV")
grhog2  = file_air.var("DNSG2")
gptdev  = file_air.var("PT")
gpdevg2 = file_air.var("PRS")
gqv     = file_air.var("QV")
#gqc     = file_air.var("QC")
#gqr     = file_air.var("QR")
#gqi     = file_air.var("QCI")
#gqs     = file_air.var("QS")
#gqg     = file_air.var("QG")
gslp    = file_air.var("PSEAsrf")
gptsfc  = file_air.var("PTGRDsrf")

#  land data
file_land = NuSDaS.open(fname_land)
gtg      = file_land.var("TUGD")
gsst     = file_land.var("SSTsrf")
gkind    = file_land.var("KINDsrf")
lon_land = gtg.dim("lon").val(:full).to_type(NArray::SFLOAT)
lat_land = gtg.dim("lat").val(:full).to_type(NArray::SFLOAT)
nlon_land, nlat_land = lon_land.shape
nl = gtg.dim("plane").val.length


time = grhog2.dim("basetime").val[trange] / (24*60)
ntime = time.length

lsmask = gkind[true,true,0]
lsmask[ (lsmask-1.0).abs.gt(1e-4) ] = 0.0
lsmask = lsmask.to_s

slon = lon.to_s
slat = lat.to_s
slon_land = lon_land.to_s
slat_land = lat_land.to_s

szs = zs.to_s

if debug

class DFile
  BYTE = 4
  class << self
    def open(name, mode)
      file = DFile.new
      yield(file)
    end
    def nlon=(n); @@nlon = n; end
    def nlat=(n); @@nlat = n; end
    def nlev=(n); @@nlev = n; end
    def nlon_land=(n); @@nlon_land = n; end
    def nlat_land=(n); @@nlat_land = n; end
    def nl=(n); @@nl = n; end
  end

  def write(str)
    len = str.length
    case len
    when @@nlon*@@nlat*BYTE
      dat = NArray.to_na(str,NArray::SFLOAT).reshape!(@@nlon,@@nlat).ntoh
      p [dat.min, dat.max]
    when @@nlon*@@nlat*@@nlev*BYTE
      dat = NArray.to_na(str,NArray::SFLOAT).reshape!(@@nlon,@@nlat,@@nlev)
      dat0 = dat[true,true,0].ntoh
      dat1 = dat[true,true,-1].ntoh
      p [dat0.min, dat0.max, dat1.min, dat1.max]
    when @@nlon_land*@@nlat_land*BYTE
      dat = NArray.to_na(str,NArray::SFLOAT).reshape!(@@nlon_land,@@nlat_land).ntoh
      p [dat.min, dat.max]
    when @@nlon_land*@@nlat_land*@@nl*BYTE
      dat = NArray.to_na(str,NArray::SFLOAT).reshape!(@@nlon_land,@@nlat_land,@@nl)
      dat0 = dat[true,true,0].ntoh
      dat1 = dat[true,true,-1].ntoh
      p [dat0.min, dat0.max, dat1.min, dat1.max]
    else
      p [len, @@nlon, @@nlat, @@nlev, @@nlon_land, @@nlat_land, @@nl]
      raise "error"
    end
    return len
  end

  def print
  end

end
DFile.nlon = nlon
DFile.nlat = nlat
DFile.nlev = nlev
DFile.nlon_land = nlon_land
DFile.nlat_land = nlat_land
DFile.nl = nl

else
  DFile = File
end

t0 = DateTime.new(1801,1,1)


unless skip_data

  ntime.times do |n|
    print "#{n+1}/#{ntime}\n"

    t = (t0 + time[n]).strftime("%Y%m%d%H")
    dir = t[0,6]
    system "mkdir #{dir}" unless File.exist?(dir)

    fname_out = "#{dir}/#{basename_out}_#{t}.grd"
    fname_land_out = "#{dir}/#{basename_out}land_#{t}.grd"
#    File.open(fname_out, "w") do |ofile|
#    File.open(fname_land_out, "w") do |ofile_land|
    DFile.open(fname_out, "w") do |ofile|
    DFile.open(fname_land_out, "w") do |ofile_land|

      # read data at n time step
      ru = gru[xrangeu,yrange ,zrange,validtime,n]
      rv = grv[xrange ,yrangev,zrange,validtime,n]
      rhog2 = grhog2[xrange,yrange,zrange,validtime,n]
      ptdev = gptdev[xrange,yrange,zrange,validtime,n]
      pdevg2 = gpdevg2[xrange,yrange,zrange,validtime,n]
      qv = gqv[xrange,yrange,zrange,validtime,n]
#      qc = gqc[xrange,yrange,zrange,validtime,n]
#      qr = gqr[xrange,yrange,zrange,validtime,n]
#      qi = gqi[xrange,yrange,zrange,validtime,n]
#      qs = gqs[xrange,yrange,zrange,validtime,n]
#      qg = gqg[xrange,yrange,zrange,validtime,n]
      slp = gslp[xrange,yrange,validtime,n] * 100.0 # hPa to Pa
      ptsfc = gptsfc[xrange,yrange,validtime,n] + PTref

      tg   = gtg[true,true,true,n]
      sst = gsst[true,true,n]

      # interpolate to full level
      ru_fl = ( ru[0..-2,true,true] + ru[1..-1,true,true] ) * 0.5
      rv_fl = ( rv[true,0..-2,true] + rv[true,1..-1,true] ) * 0.5

      # remove metrics and density
      u = ru_fl * mapf / rhog2
      v = rv_fl * mapf / rhog2
      # rho = rhog2 / g2
      pdev = pdevg2 / g2

      # rotation
      ulon =   u * cos_alpha + v * sin_alpha
      vlat = - u * sin_alpha + v * cos_alpha

      # total pressure
      pres = pbar + pdev

      # temperature
      pott = ptdev + PTref
      tem  = pott * (pres/P00)**RovCP
      tem = tem.to_type(NArray::SFLOAT)

      # surface pressure
      # estimate from the pressure at the second layer
      tem2 = tem[true,true,1] * ( 1.0 + 0.608 * qv[true,true,1] )
      dz2 = z[true,true,1] - zs
      psfc = pres[true,true,1] * ( 1.0 + Gamma*dz2/tem2 )**(Grav/(Gamma*R))
      psfc = psfc.to_type(NArray::SFLOAT)

      # surface temperature
      tsfc = ptsfc * (psfc/P00)**RovCP
      tsfc = tsfc.to_type(NArray::SFLOAT).to_s

      # write coordinates
      ofile.write slon
      ofile.write slat
      ofile.write pres.to_s

      # write 2D variables
      ofile.write slp.to_s
      ofile.write psfc.to_s
      ofile.write ulon[true,true,0].to_s  # U10
      ofile.write vlat[true,true,0].to_s  # V10
      ofile.write tsfc                 # T2
      ofile.write qv[true,true,0].to_s # Q2
      ofile.write szs                  # topo

      # write 3D variables
      ofile.write sz
      ofile.write ulon.to_s
      ofile.write vlat.to_s
      ofile.write tem.to_s
      ofile.write qv.to_s
#      ofile.write qc.to_s
#      ofile.write qr.to_s
#      ofile.write qc.to_s
#      ofile.write qi.to_s
#      ofile.write qs.to_s
#      ofile.write qg.to_s


      # write land variables
      ofile_land.write slon_land
      ofile_land.write slat_land
      ofile_land.write lsmask
      ofile_land.write tg[true,true,0].to_s # SKINT
      ofile_land.write sst.to_s # SST
      ofile_land.write tg.to_s  # STEMP
    end

  end
  end

end

exit if debug


llev = [0.02, 0.115, 0.39]

vars_air = [
    ["long", 1, "longitude [degree]"],
    ["lati", 1, "latitude [degree]"],
    ["plev", nlev, "pressure level [Pa]"],
    ["MSLP", 1, "mean sea level Pressure [Pa]"],
    ["PSFC", 1, "surface Pressure [Pa]"],
    ["U10", 1, "10 m above ground U-Component of Wind [m/s]"],
    ["V10", 1, "10 m above ground V-Component of Wind [m/s]"],
    ["T2", 1, "2 m above ground Temperature [K]"],
    ["Q2", 1, "2 m above ground Specific Humidity [kg/kg]"],
    ["TOPO", 1, "topographic elevation [m]"],
    ["HGT", nlev, "Altitude [m]"],
    ["U", nlev, "U-Component of Wind [m/s]"],
    ["V", nlev, "V-Component of Wind [m/s]"],
    ["T", nlev, "Tmeperature [K]"],
    ["QV",nlev, "Specific humidity [kg/kg]"]
]

vars_land = [
    ["long_sfc", 1, "longitude [degree]"],
    ["lati_sfc", 1, "latitude [degree]"],
    ["lsmask", 1, "land mask [0-1]"],
    ["SKINT", 1, "skin Temperature [K]"],
    ["SST", 1, "sea surface Temperature [K]"],
    ["STEMP", nl, "soil Temperature [K]"]
]

File.open(grads_namelist, "w") do |file|
  file.print <<-EOL
#
# Dimension
#
&nml_grads_grid
 outer_nx     = #{nlon},
 outer_ny     = #{nlat},
 outer_nz     = #{nlev},
 outer_nl     = #{nl},
 outer_nx_sfc = #{nlon_land},
 outer_ny_sfc = #{nlat_land},
/

#
# Variables
#
  EOL

  trec = nlev * 6 + 8
  i = 1
  vars_air.each do |v,n,desc|
    file.print <<-EOL
&grdvar item='#{"%-9s"%(v+"',")} dtype='map', fname='#{basename_out}', startrec=#{"%3d"%i}, totalrec=#{trec} /
    EOL
    i = i + n
  end

  file.print <<-EOL
&grdvar item='llev',    dtype='levels', lnum=3, lvars=#{llev.join(", ")}, /
  EOL

  trec = nl + 5
  i = 1
  vars_land.each do |v,n,desc|
    file.print <<-EOL
&grdvar item='#{"%-9s"%(v+"',")} dtype='map', fname='#{basename_out}land', startrec=#{"%1d"%i}, totalrec=#{trec} /
    EOL
    i = i + n
  end

end


t00 = t0 + time[0]

grads_ctl = "#{basename_out}#{t00.strftime("%Y%m")}.ctl"
grads_land_ctl = "#{basename_out}land#{t00.strftime("%Y%m")}.ctl"

t0s = (t00).strftime("%Hz%d%b%Y").downcase

nx0 = basepoint[0]
ny0 = basepoint[1]
dx = ( lon[nx0,ny0-1]-lon[nx0-2,ny0-1] ) * 0.5
dy = ( lat[nx0-1,ny0-2]-lat[nx0-1,ny0] ) * 0.5

File.open(grads_ctl, "w") do |file|
  file.print <<-EOL
dset ^%y4%m2/#{basename_out}_%y4%m2%d2%h2.grd
undef 9.999E+20
title #{basename_out}
options template big_endian
pdef #{nlon} #{nlat} lcc #{basepoint2[0]} #{basepoint2[1]} #{nx0} #{nlat-ny0+1} #{standard[0]} #{standard2[0]} #{standard[1]} #{"%.2f"%distance[0]} #{"%.2f"%distance[1]}
xdef #{((lon.max.ceil-lon.min.floor)/dx).ceil} linear #{lon.min.floor} #{"%.5f"%dx}
ydef #{((lat.max.ceil-lat.min.floor)/dy).ceil} linear #{lat.min.floor} #{"%.5f"%dy}
zdef #{nlev} levels #{zrp.to_a.join(" ")}
tdef #{ntime} linear #{t0s} 3hr
VARS #{vars_air.length}
  EOL
  vars_air.each do |vn, n, dec|
    file.print <<-EOL
#{vn} #{n==1 ? 0 : n} 99 #{dec}
    EOL
  end
  file.print <<-EOL
ENDVARS
  EOL
end

File.open(grads_land_ctl, "w") do |file|
  file.print <<-EOL
dset ^\%y4\%m2/#{basename_out}land_\%y4\%m2\%d2\%h2.grd
undef 9.999E+20
title #{basename_out}
options template big_endian
pdef #{nlon_land} #{nlat_land} lcc #{basepoint2[0]} #{basepoint2[1]} #{nx0} #{nlat_land-ny0+1} #{standard[0]} #{standard2[0]} #{standard[1]} #{"%.2f"%distance[0]} #{"%.2f"%distance[1]}
xdef #{((lon_land.max.ceil-lon_land.min.floor)/dx).ceil} linear #{lon_land.min.floor} #{"%.5f"%dx}
ydef #{((lat_land.max.ceil-lat_land.min.floor)/dy).ceil} linear #{lat_land.min.floor} #{"%.5f"%dy}
zdef #{nl} levels #{llev.join(" ")}
tdef #{ntime} linear #{t0s} 3hr
VARS #{vars_land.length}
  EOL
  vars_land.each do |vn, n, dec|
    file.print <<-EOL
#{vn} #{n==1 ? 0 : n} 99 #{dec}
    EOL
  end
  file.print <<-EOL
ENDVARS
  EOL
end
