require "numru/dcl"
include NumRu
include NMath

T0 = 273.15
P0 = 1000.0 # [hPa]
Rdry = 287.04
CPdry = 1004.64
Rvap = 461.50
LHV = 2.501e6
PSAT0 = 6.1078 # [hPa]

filename = ARGV.shift

nz = z = pres = temp = tdew = temp_p = nil
File.open(filename, "r") do |file|

  nz = file.gets.to_i
  z = NArray.sfloat(nz)
  pres = NArray.sfloat(nz)
  temp = NArray.sfloat(nz)
  tdew = NArray.sfloat(nz)
  temp_p = NArray.sfloat(nz)

  file.gets

  n = 0
  while line = file.gets
    dat = line.split.map{|str| str.to_f}
    z[n]      = dat[0]
    pres[n]   = dat[1]
    temp[n]   = dat[2]
    tdew[n]   = dat[3]
    temp_p[n] = dat[6]
    n += 1
  end
end

z = z / 1000.0 # m -> km
temp = temp - T0 # K -> oC
tdew = tdew - T0 # K -> oC
temp_p = temp_p - T0 # K -> oC

kmax = pres.gt(100.0).where[-1]


DCL.swiset("ifl", 1) # PNG output
#DCL.gropn(1)
DCL.gropn(2)
DCL.uzfact(0.7)
DCL.sglset("lclip", true)

vxmin = 0.10
vxmax = 0.93
vymin = 0.05
vymax = 0.68
xmin = -30.0
xmax =  40.0
ymin = pres[0]
ymax = pres[kmax]
DCL.sglset("lfull", true)
DCL.grfrm
DCL.grsvpt(vxmin, vxmax, vymin, vymax)
DCL.grswnd(xmin, xmax, ymin, ymax)
DCL.grstrn(2) # lin-log plot
DCL.grstrf

DCL.uscset("cyside", "l")
DCL.usdaxs
DCL.uysttl("l", "pressure [hPa]", 0)

yfact = log(ymax/ymin) / (vymax - vymin)
xfact = (xmax - xmin) / (vxmax - vxmin)
shift = log(pres/ymin) * xfact / yfact

temp   = temp   + shift
tdew   = tdew   + shift
temp_p = temp_p + shift

DCL.sgplzu(temp,   pres, 1,  2) # black
DCL.sgplzu(tdew,   pres, 1, 42) # blue
DCL.sgplzu(temp_p, pres, 1, 22) # red


# auxiliary lines
## temperature
for t10 in ((xmin/10.0).ceil-10)..((xmax/10.0).floor)
  t0 = t10 * 10.0
  t = NArray.sfloat(nz).fill(t0) + shift
  DCL.sgplzu(t, pres, 3, 1) # dashed black line
end
## dry adiabatic lapslate
for t10 in ((xmin/10.0).ceil)..((xmax/10.0).floor+20)
  t0 = t10 * 10.0
  t = (t0 + T0) * (pres/pres[0])**(Rdry/CPdry) - T0 + shift
  DCL.sgplzu(t, pres, 3, 31) # dashed green line
end
## moist adiabatic lapslate
eps = Rdry / Rvap
a = Rdry / CPdry
b = LHV**2 * eps / ( CPdry * Rdry )
c = LHV / CPdry
ni = 10
for t10 in ((xmin/5.0).ceil)..((xmax/5.0).floor+10)
  t0 = t10 * 5.0
  t = NArray.sfloat(nz)
  t[0] = t0

  for k in 1...nz
    tk = t[k-1] + T0
    pr = pres[k-1]
    dlp = log( pres[k] / pres[k-1] ) / ni
    ni.times do
      es = PSAT0 * exp( LHV / Rvap * ( 1.0/T0 - 1.0/tk ) )
      rs = eps * es / ( pr - es )
      dtdlp = (a*tk + c*rs) / (1.0 + b*rs/tk**2)
      tk = tk + dtdlp * dlp
      pr = pr * exp(dlp)
    end
    t[k] = tk - T0
  end
  t = t + shift
  DCL.sgplzu(t, pres, 3, 791) # dashed orange line
end
## dew point
a = 6.1094
b = 17.625
c = 243.04
for t10 in ((xmin/10.0).ceil-5)..((xmax/10.0).floor)
  t0 = t10 * 10.0
  e = a * exp(b*t0/(t0+c)) * pres / pres[0]
  t = c * log(e/a) / (b-log(e/a)) + shift
  DCL.sgplzu(t, pres, 3, 91) # dashed cyan line
end




# right y axis
DCL.grfig
DCL.grsvpt(vxmin, vxmax, vymin, vymax)
DCL.grswnd(xmin, xmax, z[0], z[kmax])
DCL.grstrn(1)
DCL.grstrf

DCL.uzlset("labelyr", true)
DCL.usyaxs("r")
DCL.uysttl("r", "height [km]", 0)
DCL.uxsttl("t", "Skew-T Log-P diagram", 0)


DCL.grcls
