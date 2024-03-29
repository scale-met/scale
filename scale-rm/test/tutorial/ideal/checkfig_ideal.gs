**
** Initial set
**
'set display color white'
'set xlab off'
'set ylopts 1 4 0.15'
'set mproj off'
'set ylint 20'
'c'
**
** Open file
**
'sdfopen merged_history.pe000000.nc'
'set x 1'
'set t 5'
'set z 1 70'
**
** QHYD
**
'set parea 1 9 1 7.5'
'set grads off'
'set gxout shaded'
'set clevs 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5'
'd qhyd*1000'
'set gxout contour'
'set cthick 0.5'
'set clevs 0.5'
'd qhyd*1000'
*
'set ccolor 1'
'set cthick 3'
'd skip(v,2);w'
'cbarn'
'draw title (a) QHYD (10^3 kg/kg) & V;W (m/s)'
'printim ideal_qhyd.png'
'c'
**
** W
**
'set parea 1 9 1 7.5'
'set grads off'
'set gxout shaded'
'set clevs -8 -6.4 -4.8 -3.2 -1.6 0 1.6 3.2 4.8 6.4 8'
'd W'
'set gxout contour'
'set cthick 0.5'
'set clevs 0.5'
'd W'
*
'set ccolor 1'
'set cthick 3'
'd skip(v,2);w'
'cbarn'
'draw title (b) W (m/s) & V;W (m/s)'
'printim ideal_W.png'
'c'
'quit'


