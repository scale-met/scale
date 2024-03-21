!-------------------------------------------------------------------------------
!> Program cnvtopo
!!
!! @par Description
!!          Convert topography based on the icosahedral grid configuration
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program ll2ico
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
     PRC_MPIstart,         &
     PRC_SINGLECOM_setup,  &
     PRC_ERRHANDLER_setup, &
     PRC_MPIfinish
  use scale_prc_icoA, only: &
     PRC_ICOA_setup
  use scale_const, only: &
     CONST_setup
  use scale_atmos_grid_icoA_index, only: &
     ATMOS_GRID_icoA_INDEX_setup
  use scale_comm_icoA, only: &
     COMM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_gm_topography, only: &
     TOPOGRAPHY_setup, &
     TOPOGRAPHY_write
  use mod_grd, only: &
     GRD_setup
  use mod_gmtr, only: &
     GMTR_setup
  use mod_oprt, only: &
     OPRT_setup
  use mod_gm_convert, only: &
     CONVERT_setup, &
     CONVERT
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

  character(len=H_LONG) :: cnf_fname ! config file

  integer :: comm
  integer :: nprocs
  integer :: myrank
  logical :: ismaster
  !=============================================================================

  ! start MPI
  call PRC_MPIstart( comm ) ! [OUT]

  ! setup MPI communicator
  call PRC_SINGLECOM_setup( comm,    & ! [IN]
                            nprocs,  & ! [OUT]
                            myrank,  & ! [OUT]
                            ismaster ) ! [OUT]

  ! setup errhandler
  call PRC_ERRHANDLER_setup( .false., & ! [IN]
                             ismaster ) ! [IN]

  !########## Initial setup ##########

  ! setup standard I/O
  cnf_fname = IO_ARG_getfname( ismaster )

  call IO_setup( MODELNAME, cnf_fname )
  call IO_LOG_setup( myrank, ismaster )

  ! setup process
  call PRC_ICOA_setup

  ! setup PROF
  call PROF_setup

  !#############################################################################
  call PROF_setprefx('INIT')
  call PROF_rapstart('Initialize',0)

  !---< cnst module setup >---
  call CONST_setup

  ! setup horizontal/vertical grid coordinates (icosahedral,idealized)
  call ATMOS_GRID_icoA_INDEX_setup

  !---< I/O module setup >---
  call FIO_setup

  !---< comm module setup >---
  call COMM_setup

  !---< topography module setup >---
  call TOPOGRAPHY_setup

  !---< grid module setup >---
  call GRD_setup

  !---< geometrics module setup >---
  call GMTR_setup

  !---< operator module setup >---
  call OPRT_setup

  !---< preprocess converter setup >---
  call CONVERT_setup

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_ll2ico',0)

  ! execute preprocess
  call PROF_rapstart('Convert',1)
  call CONVERT
  call PROF_rapend  ('Convert',1)

  !---< topography output >---
  call TOPOGRAPHY_write

  call PROF_rapend('Main_ll2ico',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program ll2ico
