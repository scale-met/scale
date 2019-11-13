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
**
** MSLP
**
'set parea 1 9 1 7.5'
'open MSLP_d01z-2d_lccr.ctl'
'set grads off'
'set t 7'
'set clevs 98000 98500 99000 99500 100000 100500 101000'
'd mslp'
'cbarn'
'draw title MSLP (Pa)'
'printim real_mslp.png'
'c'
'close 1'
**
** PREC
**
'set parea 1 9 1 7.5'
'open PREC_d01z-2d_lccr.ctl'
'set grads off'
'set t 7'
'set clevs 0 5 10 15 20 25 30 35 40'
'd prec*3600'
'cbarn'
'draw title PREC (mm/hr)'
'printim real_prec.png'
'c'
'close 1'
**
** WIND
**
'set parea 1 9 1 7.5'
'open Umet_d01z-3d_lccr.ctl'
'open Vmet_d01z-3d_lccr.ctl'
'set grads off'
***'set lev 500'
'set lev 850'
'set t 7'
'd sqrt(umet*umet+vmet.2*vmet.2)'
'cbarn'
'set clevs 5 10 15 20 25 30 35 40 45 50'
'd umet;skip(vmet.2,5)'
'draw title Wind (m/s) @ 850hPa'
'printim real_wind.png'
'c'
'close 2'
'close 1'
'quit'

