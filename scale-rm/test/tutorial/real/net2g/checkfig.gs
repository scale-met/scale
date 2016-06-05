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
'open MSLP_d01z-2d_lcc.ctl'
'set grads off'
'set t 25'
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
'open PREC_d01z-2d_lcc.ctl'
'set grads off'
'set t 25'
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
'open U_d01z-3d_lcc.ctl'
'open V_d01z-3d_lcc.ctl'
'set grads off'
'set lev 500'
'set t 25'
'd sqrt(u*u+v.2*v.2)'
'cbarn'
'd u;skip(v.2,5)'
'draw title Wind (m/s) @ 500hPa'
'printim real_wind.png'
'c'
'close 2'
'close 1'
'quit'




