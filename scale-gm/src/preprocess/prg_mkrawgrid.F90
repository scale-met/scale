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
