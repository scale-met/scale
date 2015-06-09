!-------------------------------------------------------------------------------
!
!+  Program driver
!
!-------------------------------------------------------------------------------
program prg_driver
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This program is a driver of non-hydrostatic model based on an
  !       icosahedral grid system.
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !      0.03      04-05-31   Change by addtion of mod[onestep].
  !                05-12-01   M.Satoh add history_setup
  !                05-12-19   S.Iga moved output_timeinfo after output_all
  !                06-04-18   T.Mitsui add sfc_restart
  !                06-04-21   H.Tomita:  remove output_timeinfo due to
  !                                      computational efficeincy.
  !                                      Instead, this process is move to
  !                                      mod[mod_output].
  !                06-08-07   W.Yanase add history_vars
  !                06-09-27   S.Iga add history_vars_cfmip
  !                07-03-23   Y.Niwa add NDG_setup, ndg_do, FLAG_NUDGING
  !                07-06-27   Y.Niwa add history_vars_setup
  !                07-07-24   K.Suzuki: implementing SPRINTARS aerosol model
  !                07-08-06   Y.Niwa: add history_obs, history_vars_obs
  !                07-11-07   T.Mitsui: add option to omit output_all
  !                08-03-10   T.Mitsui: add intermediate output of restart file
  !                08-05-24   T.Mitsui: trivial fix
  !                08-09-09   Y.Niwa : modfied for nudging
  !                09-01-23   H.Tomita: a) abolish mod_output, mod_extdata
  !                                     mod_history_vars_cfmip, mod_o3var.
  !                                     b) introduce mod_extdata.
  !                09-04-14   T.Mitsui: arrange initialization of aerosols
  !                09-07-10   H.Tomita: Add the module [mod_embudget].
  !                09-08-05   S.Iga: remove latlon_setup (suggested by T.Mitsui)
  !                09-08-05   T.Mitsui: add conditioning by ADM_myprc_is_run
  !                                     to keep out extra-processes from main routines.
  !                09-08-18   T.Mitsui: change of 09-08-05 is not enough.
  !                10-03-08   C.Kodama: Modify for overwrite_restart option
  !                10-04-30   M.Satoh: move diagvar_setup
  !                11-09-03   H.Yashiro : New I/O
  !                11-11-28   Y.Yamada : merge Terai-san timer
  !                12-06-07   T.Seiki  : Application to Multi-job System
  !                12-10-12   R.Yoshida  : Modify for Dynamical Core test
  !                12-10-22   R.Yoshida  : add papi instructions
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_MULTI_PRC,      &
     ADM_LOG_FID,        &
     ADM_prc_me,         &
     ADM_prc_run_master, &
     ADM_proc_init,      &
     ADM_proc_finish,    &
     ADM_setup
  use mod_random, only: &
     RANDOM_setup
  use mod_cnst, only: &
     CNST_setup
  use mod_calendar, only: &
     calendar_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_grd, only: &
     GRD_setup
  use mod_gmtr, only: &
     GMTR_setup
  use mod_oprt, only: &
     OPRT_setup
  use mod_vmtr, only: &
     VMTR_setup
  use mod_time, only: &
     TIME_setup,     &
     TIME_report,    &
     TIME_advance,   &
     TIME_LSTEP_MAX, &
     TIME_CSTEP,     &
     TIME_CTIME,     &
     TIME_DTL
  use mod_extdata, only: &
     extdata_setup
  use mod_runconf, only: &
     runconf_setup
  use mod_prgvar, only: &
     prgvar_setup,            &
     restart_input_basename,  &
     restart_output_basename, &
     restart_input,           &
     restart_output
  use mod_dynamics, only: &
     dynamics_setup, &
     dynamics_step
  use mod_forcing_driver, only: &
     forcing_setup, &
     forcing_step
  use mod_history, only: &
     history_setup, &
     history_out,   &
     HIST_output_step0
  use mod_history_vars, only: &
     history_vars_setup, &
     history_vars
  use mod_embudget, only: &
     embudget_setup, &
     embudget_monitor

  !##### OpenACC (for data copy) #####
  use mod_adm, only: &
     ADM_prc_tab,  &
     ADM_rgn_vnum, &
     ADM_IopJop
  use mod_comm, only: &
     sendlist,     sendlist_pl,  &
     sendinfo,     sendinfo_pl,  &
     recvlist,     recvlist_pl,  &
     recvinfo,     recvinfo_pl,  &
     recvlist_r2r, sendlist_r2r, &
     recvlist_r2p, sendlist_r2p, &
     recvlist_p2r, sendlist_p2r, &
     recvlist_sgp, sendlist_sgp, &
     copyinfo_r2r, copyinfo_sgp, &
     copyinfo_r2p, copyinfo_p2r, &
     nsmax,        nsmax_pl,     &
     nrmax,        nrmax_pl,     &
     ncmax_r2r,    ncmax_sgp,    &
     ncmax_r2p,    ncmax_p2r
  use mod_grd, only: &
     GRD_x,     &
     GRD_xt,    &
     GRD_zs,    &
     GRD_rdgz,  &
     GRD_rdgzh, &
     GRD_vz
  use mod_gmtr, only: &
     GMTR_P_var, &
     GMTR_T_var, &
     GMTR_A_var
  use mod_oprt, only: &
     cdiv,        &
     cgrad,       &
     clap,        &
     cinterp_TN,  &
     cinterp_HN,  &
     cinterp_TRA, &
     cinterp_PRA
  use mod_vmtr, only: &
     VMTR_GAM2,      &
     VMTR_GAM2H,     &
     VMTR_GSGAM2,    &
     VMTR_GSGAM2H,   &
     VMTR_RGSQRTH,   &
     VMTR_RGAM,      &
     VMTR_RGAMH,     &
     VMTR_RGSGAM2,   &
     VMTR_RGSGAM2H,  &
     VMTR_W2Cfact,   &
     VMTR_C2Wfact,   &
     VMTR_C2WfactGz, &
     VMTR_PHI
  use mod_runconf, only: &
     CVW
  use mod_prgvar, only: &
     PRG_var,  &
     PRG_var1, &
     DIAG_var
  use mod_bsstate, only: &
     rho_bs, &
     pre_bs, &
     tem_bs
  use mod_numfilter, only: &
     Kh_coef,      &
     Kh_coef_lap1, &
     divdamp_coef
  use mod_vi, only : &
     Mc, &
     Ml, &
     Mu
  use mod_history, only: &
     ksumstr,     &
     cnvpre_klev, &
     cnvpre_fac1, &
     cnvpre_fac2
  !##### OpenACC #####
  implicit none

  character(len=14) :: cdate

  integer :: n
  !-----------------------------------------------------------------------------

  call ADM_proc_init(ADM_MULTI_PRC)

  !---< admin module setup >---
  call ADM_setup('nhm_driver.cnf')

  !#############################################################################

  write(ADM_LOG_FID,*) '##### start  setup     #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### start  setup     #####'
  endif

  call DEBUG_rapstart('Total')
  call DEBUG_rapstart('Setup_ALL')

  !---< radom module setup >---
  call RANDOM_setup

  !---< cnst module setup >---
  call CNST_setup

  !---< calendar module setup >---
  call calendar_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< grid module setup >---
  call GRD_setup

  !---< geometrics module setup >---
  call GMTR_setup

  !---< operator module setup >---
  call OPRT_setup

  !---< vertical metrics module setup >---
  call VMTR_setup

  !---< time module setup >---
  call TIME_setup

  !---< external data module setup >---
  call extdata_setup


  !---< nhm_runconf module setup >---
  call runconf_setup

  !---< prognostic variable module setup >---
  call prgvar_setup
  call restart_input( restart_input_basename )


  !---< dynamics module setup >---
  call dynamics_setup

  !---< forcing module setup >---
  call forcing_setup

  !---< energy&mass budget module setup >---
  call embudget_setup

  !---< history module setup >---
  call history_setup

  !---< history variable module setup >---
  call history_vars_setup

  write(ADM_LOG_FID,*) '##### finish setup     #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### finish setup     #####'
  endif


  call DEBUG_rapend('Setup_ALL')

  !#############################################################################
#ifdef _FIPP_
  call fipp_start()
#endif
  call DEBUG_rapstart('Main_ALL')

  write(ADM_LOG_FID,*) '##### start  main loop #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### start  main loop #####'
  endif

  !$acc data &
  !$acc& pcopyin(ADM_prc_tab,ADM_rgn_vnum,ADM_IopJop) &
  !$acc& pcopyin(sendlist,sendlist_pl) &
  !$acc& pcopyin(sendinfo,sendinfo_pl) &
  !$acc& pcopyin(recvlist,recvlist_pl) &
  !$acc& pcopyin(recvinfo,recvinfo_pl) &
  !$acc& pcopyin(recvlist_r2r,sendlist_r2r) &
  !$acc& pcopyin(recvlist_sgp,sendlist_sgp) &
  !$acc& pcopyin(recvlist_r2p,sendlist_r2p) &
  !$acc& pcopyin(recvlist_p2r,sendlist_p2r) &
  !$acc& pcopyin(copyinfo_r2r,copyinfo_sgp,copyinfo_r2p,copyinfo_p2r) &
  !$acc& pcopyin(nsmax,nsmax_pl,nrmax,nrmax_pl) &
  !$acc& pcopyin(ncmax_r2r,ncmax_sgp,ncmax_r2p,ncmax_p2r) &
  !$acc& pcopyin(GRD_rdgz,GRD_rdgzh,GRD_x,GRD_xt,GRD_vz,GRD_zs) &
  !$acc& pcopyin(GMTR_P_var,GMTR_T_var,GMTR_A_var) &
  !$acc& pcopyin(cdiv,cgrad,clap,cinterp_TN,cinterp_HN,cinterp_TRA,cinterp_PRA) &
  !$acc& pcopyin(VMTR_GAM2,VMTR_GAM2H,VMTR_GSGAM2,VMTR_GSGAM2H) &
  !$acc& pcopyin(VMTR_RGSQRTH,VMTR_RGAM,VMTR_RGAMH,VMTR_RGSGAM2,VMTR_RGSGAM2H) &
  !$acc& pcopyin(VMTR_W2Cfact,VMTR_C2Wfact,VMTR_C2WfactGz,VMTR_PHI) &
  !$acc& pcopyin(CVW) &
  !$acc& pcopyin(rho_bs,pre_bs,tem_bs) &
  !$acc& pcopyin(divdamp_coef,Kh_coef,Kh_coef_lap1) &
  !$acc& pcopyin(Mc,Mu,Ml) &
  !$acc& pcopyin(ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2) &
  !$acc& pcopy  (PRG_var,PRG_var1,DIAG_var)

  !--- history output at initial time
  if ( HIST_output_step0 ) then
     TIME_CSTEP = TIME_CSTEP - 1
     TIME_CTIME = TIME_CTIME - TIME_DTL
     call history_vars
     call TIME_advance
     call history_out
  else
     call TIME_report
  endif

  do n = 1, TIME_LSTEP_MAX

     call DEBUG_rapstart('_Atmos')
     call dynamics_step
     call forcing_step
     call DEBUG_rapend  ('_Atmos')

     call DEBUG_rapstart('_History')
     call history_vars
     call TIME_advance

     !--- budget monitor
     call embudget_monitor
     call history_out

     if (n == TIME_LSTEP_MAX) then
        cdate = ""
        call restart_output( restart_output_basename )
     endif
     call DEBUG_rapend  ('_History')

  enddo

  !$acc end data
  write(ADM_LOG_FID,*) '##### finish main loop #####'
  if ( ADM_prc_me == ADM_prc_run_master ) then
     write(*,*) '##### finish main loop #####'
  endif

  call DEBUG_rapend('Main_ALL')
#ifdef _FIPP_
  call fipp_stop()
#endif
  !#############################################################################

  call DEBUG_rapend('Total')
  call DEBUG_rapreport

  !--- finalize all process
  call ADM_proc_finish

  stop
end program prg_driver
!-------------------------------------------------------------------------------
