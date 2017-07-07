!-------------------------------------------------------------------------------
!> module SCALE-GM (a main routine of global model)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for global scale weather/climate
!!          based on the Non-hydrostatic ICosahedral Atmosphere Model (NICAM)
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_gm_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_log, only: &
     LogInit
  use gtool_file, only: &
     FileCloseAll
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-gm.h"
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: scalegm

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_MID), private, parameter :: MODELNAME = "SCALE-GM ver. "//VERSION

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine scalegm( &
       comm_world,       &
       intercomm_parent, &
       intercomm_child,  &
       cnf_fname         )
    use scale_process, only: &
       PRC_IsMaster,    &
       PRC_MPIstart,    &
       PRC_LOCAL_setup, &
       PRC_MPIfinish
    use scale_const, only: &
       CONST_setup
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use mod_adm, only: &
       ADM_setup
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
    use mod_saturation, only: &
       saturation_setup
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
       ADM_rgn_vnum
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
       GMTR_p, &
       GMTR_t, &
       GMTR_a
    use mod_vmtr, only: &
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

    integer,          intent(in) :: comm_world
    integer,          intent(in) :: intercomm_parent
    integer,          intent(in) :: intercomm_child
    character(len=*), intent(in) :: cnf_fname

    integer :: myrank
    logical :: ismaster

    integer :: n
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

    ! setup standard I/O
    call IO_setup( MODELNAME, .true., cnf_fname )

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    call LogInit( IO_FID_CONF,       &
                  IO_FID_LOG, IO_L,  &
                  IO_FID_NML, IO_NML )

    !---< admin module setup >---
    call ADM_setup

    ! setup PROF
    call PROF_setup

    !###########################################################################
    ! profiler start
    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize',0)


    ! setup constants
    call CONST_setup

    ! setup calendar
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup

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

    !---< saturation module setup >---
    call saturation_setup

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

    call PROF_rapend('Initialize', 0)

    !########## main ##########

#ifdef FIPP
    call fipp_start
#endif
#ifdef PAPI
    call PROF_PAPI_rapstart
#endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    !$acc data &
    !$acc& pcopyin(ADM_prc_tab,ADM_rgn_vnum) &
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
    !$acc& pcopyin(GMTR_p,GMTR_t,GMTR_a) &
    !$acc& pcopyin(cdiv,cgrad,clap,cinterp_TN,cinterp_HN,cinterp_TRA,cinterp_PRA) &
    !$acc& pcopyin(VMTR_GAM2,VMTR_GAM2H,VMTR_GSGAM2,VMTR_GSGAM2H) &
    !$acc& pcopyin(VMTR_RGSQRTH,VMTR_RGAM,VMTR_RGAMH,VMTR_RGSGAM2,VMTR_RGSGAM2H) &
    !$acc& pcopyin(VMTR_W2Cfact,VMTR_C2Wfact,VMTR_C2WfactGz,VMTR_PHI) &
    !$acc& pcopyin(CVW) &
    !$acc& pcopyin(rho_bs,pre_bs,tem_bs) &
    !$acc& pcopyin(divdamp_coef,Kh_coef,Kh_coef_lap1) &
    !$acc& pcopyin(Mc,Mu,Ml) &
    !$acc& pcopyin(ksumstr,cnvpre_klev,cnvpre_fac1,cnvpre_fac2) &
    !$acc& pcopy  (PRG_var,DIAG_var)

    !--- history output at initial time
    if ( HIST_output_step0 ) then
       call TIME_advance( reverse=.true. )
       call TIME_report
       call history_vars
       call TIME_advance
       call history_out
    endif

    do n = 1, TIME_LSTEP_MAX

       call TIME_report

       call PROF_rapstart('_Atmos',1)
       call dynamics_step
       call forcing_step
       call PROF_rapend  ('_Atmos',1)

       call PROF_rapstart('_History',1)
       call history_vars
       call embudget_monitor

       call TIME_advance

       call history_out

       if (n == TIME_LSTEP_MAX) then
          call restart_output( restart_output_basename )
       endif
       call PROF_rapend  ('_History',1)

      if( IO_L ) call flush(IO_FID_LOG)

    enddo

    call PROF_rapend('Main_Loop', 0)

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
    if( IO_L ) write(IO_FID_LOG,*)


    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

#ifdef FIPP
    call fipp_stop
#endif
#ifdef PAPI
    call PROF_PAPI_rapstop
#endif

    !########## Finalize ##########

    call PROF_rapstart('File', 2)
    call FileCloseAll
    call PROF_rapend  ('File', 2)

    call PROF_rapend  ('All', 1)

    call PROF_rapreport
#ifdef PAPI
    call PROF_PAPI_rapreport
#endif

    return
  end subroutine scalegm

end module mod_gm_driver
