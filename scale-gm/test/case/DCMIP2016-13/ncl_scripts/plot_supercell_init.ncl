;***********************************************
; plot_supercell_init.ncl
;***********************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

;************************************************
begin

;************************************************
; Input parameters
;************************************************

  datafilename = "nicam.163.r200.L40.interp_latlon.nonhydro.nc"

;************************************************
; Initialize workspace
;************************************************

  wks = gsn_open_wks("eps","plot_supercell_init")   

;************************************************
; Plot resources [options]
;************************************************

  res                     = True
 ;res@mpCenterLonF        = 180
  res@cnFillOn            = True       ; color              [default=False]
  res@cnLinesOn           = True       ; draw contour lines [default=True]
  res@lbLabelAutoStride   = True
  res@vpKeepAspect        = True
  res@vpWidthF            = 0.38       ; user specified shape
  res@vpHeightF           = 0.19
  res@gsnSpreadColors     = True       ; use all colors
  res@gsnSpreadColorStart = 2          ; default=2
  res@gsnSpreadColorEnd   = 23        ; final color of "gui-default"
;  res@gsnSpreadColorStart = 2          ; default=2
;  res@gsnSpreadColorEnd   = 11        ; final color of "cosam"

  res@gsnDraw             = False      ; default=True
  res@gsnFrame            = False      ; default=True

  res@trYMaxF = 20.0

  res_p                   = res        ; pressure plot resources
  res_pd                  = res        ; pressure plot resources
  res_u                   = res        ; zonal velocity plot resources
  res_t                   = res        ; temperature plot resources
  res_qv                  = res        ; water vapor plot resources
  res_theta               = res        ; potential temperature plot resources
  res_thetad              = res        ; potential temperature diff plot resources

;************************************************
; Panel plot
;************************************************
  plot = new(6,graphic)
  gsn_merge_colormaps(wks,"gui_default", "BlWhRe")
;  gsn_merge_colormaps(wks,"cosam", "BlRe")

  datafile = addfile(datafilename, "r")

  lev = datafile->lev

;  rho = datafile->Rho
  t = datafile->T
;  theta = datafile->Theta
  u = datafile->U

;  qv = datafile->RhoQv / rho * 1000.0
  qv = datafile->Q * 1000.0
  printVarSummary (qv)
  printMinMax (qv, False)
  copy_VarCoords(t, qv)

;  p = rho * 287.0 * t / 100.0
  p = datafile->P / 100.0
  copy_VarCoords(t, p)
  theta = (t * (1.0 + 0.608*(qv/1000)))/(p/1000)^(0.2875)
  copy_VarCoords(t, theta)
  printVarSummary (theta)
  printMinMax (theta, False)

;  ps = (p(0,0,:,:) * lev(1) - p(0,1,:,:) * lev(0)) / (lev(1) - lev(0))
  ps = datafile->PS / 100.0
;  copy_VarCoords(t(0,0,:,:), ps)

  nlat = dimsizes(theta(0,0,:,0))
  ilateq = nlat/2

  thetad = theta
  do i=0,nlat-1
    thetad(0,:,i,:) = theta(0,:,i,:) - theta(0,:,ilateq,:)
  end do

  pd = p
  do i=0,nlat-1
    pd(0,:,i,:) = p(0,:,i,:) - p(0,:,ilateq,:)
  end do

; ---
;  z = lev * 20.0
  z = lev
  p&lev = z
  pd&lev = z
  t&lev = z
  theta&lev = z
  thetad&lev = z
  qv&lev = z
  u&lev = z

; ---

  res_u@cnLevelSelectionMode= "ManualLevels"
  res_u@cnLevelSpacingF   = 2.5
  res_u@cnMinLevelValF    = -15.0
  res_u@cnMaxLevelValF    = 15.0

  res_u@tiYAxisString     = "Altitude (km)"
  res_u@gsnCenterString   = "Zonal velocity"
  res_u@gsnRightString    = "(m/s)"
  plot(0) = gsn_csm_contour(wks,u(0,:,:,0),res_u) 

  res_qv@cnLevelSelectionMode= "ManualLevels"
  res_qv@cnLevelSpacingF   = 2.0
  res_qv@cnMinLevelValF    = 2.0
  res_qv@cnMaxLevelValF    = 14.0

  res_qv@tiYAxisString     = "Altitude (km)"
  res_qv@gsnCenterString   = "Specific Humidity"
  res_qv@gsnRightString    = "(g/kg)"
  plot(1) = gsn_csm_contour(wks,qv(0,:,:,143),res_qv) 

; ---

  res_theta@cnLevelSelectionMode= "ExplicitLevels"
  res_theta@cnLevels = (/300., 310., 320., 330., 340., 350., 360., 380., 400., 420., 440., 460., 480., 500./)

  res_theta@tiYAxisString     = "Altitude (km)"
  res_theta@gsnCenterString   = "Virtual Potential Temperature (~F5~q~F~~B~v~E~  )"
  res_theta@gsnRightString    = "(K)"
  plot(2) = gsn_csm_contour(wks,theta(0,:,:,0),res_theta)

  res_thetad@cnLevelSelectionMode= "ManualLevels"
  res_thetad@cnLevelSpacingF   = 0.2
  res_thetad@cnMinLevelValF    = -2.0
  res_thetad@cnMaxLevelValF    = 2.0

  res_thetad@trXMinF = 0.
  res_thetad@trXMaxF = 90.

  res_thetad@gsnSpreadColorStart = 24
  res_thetad@gsnSpreadColorEnd = 150

  res_thetad@tiYAxisString     = "Altitude (km)"
  res_thetad@gsnCenterString   = "~F5~q~F~~B~v~E~  - ~F5~q~F~~B~v,eq~E~"
  res_thetad@gsnRightString    = "(K)"
  plot(3) = gsn_csm_contour(wks,thetad(0,:,:,0),res_thetad)

; ---

 res_p@cnLevelSelectionMode= "ExplicitLevels"
 res_p@cnLevels = (/1.e1, 2.e1, 5.e1, 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 7.e2, 8.e2, 9.e2, 1.e3/)

  res_p@tiYAxisString     = "Altitude (km)"
  res_p@gsnCenterString   = "Pressure (p)"
  res_p@gsnRightString    = "(hPa)"
  plot(4) = gsn_csm_contour(wks,p(0,:,:,0),res_p) 

  res_pd@cnLevelSelectionMode= "ManualLevels"
  res_pd@cnLevelSpacingF   = 0.1
  res_pd@cnMinLevelValF    = -1.0
  res_pd@cnMaxLevelValF    = 1.0

  res_pd@trXMinF = 0.
  res_pd@trXMaxF = 90.

  res_pd@gsnSpreadColorStart = 24
  res_pd@gsnSpreadColorEnd = 150

  res_pd@tiYAxisString     = "Altitude (km)"
  res_pd@gsnCenterString   = "p - p~B~eq~E~"
  res_pd@gsnRightString    = "(hPa)"
  plot(5) = gsn_csm_contour(wks,pd(0,:,:,0),res_pd) 

;  res_t@cnLevelSelectionMode= "ManualLevels"
;  res_t@cnLevelSpacingF   = 0.1
;  res_t@cnMinLevelValF    = -1.0
;  res_t@cnMaxLevelValF    =  1.0

;  res_t@tiYAxisString     = "Altitude (km)"
;  res_t@gsnCenterString   = "Temperature"
;  res_t@gsnRightString    = "(K)"
;  plot(0) = gsn_csm_contour(wks,t(0,:,:,0),res_t) 

; ---

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelBottom = 0.1

  gsn_panel(wks,plot,(/3,2/),pres)

  frame(wks)


end
