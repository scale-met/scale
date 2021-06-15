**
** Initial set
**
'set display color white'
'set xlopts 1 4 0.15'
'set ylopts 1 4 0.15'
'set mproj off'
'set xlint 20'
'set ylint 20'
'set gxout grfill'
'c'
'set parea 1 9 1 7.5'
'sdfopen merged-h_history_d01.pe000000.nc'
'set grads off'
**
** MSLP
**
'set t 7'
'set clevs 98000 98500 99000 99500 100000 100500 101000'
'd mslp'
'cbarn'
'draw title MSLP (Pa)'
'printim real_mslp.png'
'c'
**
** PREC
**
'set t 7'
'set clevs 0 5 10 15 20 25 30 35 40'
'd prec*3600'
'cbarn'
'draw title PREC (mm/hr)'
'printim real_prec.png'
'c'
**
** WIND
**
'set lev 850'
'set t 7'
'd sqrt(umet*umet+vmet*vmet)'
'cbarn'
'set clevs 5 10 15 20 25 30 35 40 45 50'
'd umet;skip(vmet,5)'
'draw title Wind (m/s) @ 850hPa'
'printim real_wind.png'
'c'

'quit'
