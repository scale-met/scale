!-------------------------------------------------------------------------------
!> Program mkrawgrid
!!
!! @par Description
!!          Making horizontal grid systems based on the icosahedral grid configuration
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program mkrawgrid
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use dc_log, only: &
     LogInit
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_process, only: &
     PRC_LOCAL_MPIstart, &
     PRC_MPIfinish
  use scale_const, only: &
     CONST_setup
  use mod_adm, only: &
     ADM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_comm, only: &
     COMM_setup
  use mod_grd, only: &
     GRD_output_hgrid
  use mod_mkgrd, only: &
     MKGRD_setup,        &
     MKGRD_standard,     &
     MKGRD_spring,       &
     MKGRD_OUT_BASENAME, &
     MKGRD_OUT_io_mode
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-gm.h"
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  character(len=H_MID), parameter :: MODELNAME = "SCALE-GM ver. "//VERSION

  integer :: myrank
  logical :: ismaster
  !=============================================================================

  !---< MPI start >---
  call PRC_LOCAL_MPIstart( myrank,  & ! [OUT]
                           ismaster ) ! [OUT]

  !########## Initial setup ##########

  ! setup standard I/O
  call IO_setup( MODELNAME, .false. )
  call IO_LOG_setup( myrank, ismaster )
  call LogInit( IO_FID_CONF,       &
                IO_FID_LOG, IO_L,  &
                IO_FID_NML, IO_NML )

  !---< admin module setup >---
  call ADM_setup

  ! setup PROF
  call PROF_setup

  !#############################################################################
  call PROF_setprefx('INIT')
  call PROF_rapstart('Initialize',0)

  !---< cnst module setup >---
  call CONST_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< mkgrid module setup >---
  call MKGRD_setup

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_MKGRD',0)

  call PROF_rapstart('MKGRD_standard',0)
  call MKGRD_standard
  call PROF_rapend  ('MKGRD_standard',0)

  call PROF_rapstart('MKGRD_spring',0)
  call MKGRD_spring
  call PROF_rapend  ('MKGRD_spring',0)

  call GRD_output_hgrid( basename      = MKGRD_OUT_BASENAME, & ! [IN]
                         output_vertex = .false.,            & ! [IN]
                         io_mode       = MKGRD_OUT_io_mode   ) ! [IN]

  call PROF_rapend('Main_MKGRD',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program mkrawgrid
