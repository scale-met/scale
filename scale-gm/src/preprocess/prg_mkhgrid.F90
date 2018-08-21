!-------------------------------------------------------------------------------
!> Program mkhgrid
!!
!! @par Description
!!          Making horizontal grid systems based on the icosahedral grid configuration
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program mkhgrid
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
     RADIUS => CONST_RADIUS, &
     CONST_setup
  use scale_atmos_grid_icoA_index, only: &
     ATMOS_GRID_icoA_INDEX_setup
  use scale_comm_icoA, only: &
     COMM_setup
  use mod_fio, only: &
     FIO_setup
  use mod_grd, only: &
     GRD_input_hgrid,  &
     GRD_output_hgrid, &
     GRD_makelatlon,   &
     GRD_scaling
  use mod_gmtr, only: &
     GMTR_setup
  use mod_mkgrd, only: &
     MKGRD_setup,        &
     MKGRD_prerotate,    &
     MKGRD_stretch,      &
     MKGRD_shrink,       &
     MKGRD_rotate,       &
     MKGRD_gravcenter,   &
     MKGRD_diagnosis,    &
     MKGRD_IN_BASENAME,  &
     MKGRD_IN_io_mode,   &
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

  call GRD_input_hgrid( basename     = MKGRD_IN_BASENAME, & ! [IN]
                        input_vertex = .false.,           & ! [IN]
                        io_mode      = MKGRD_IN_io_mode   ) ! [IN]

  call PROF_rapstart('MKGRD_prerotate',0)
  call MKGRD_prerotate
  call PROF_rapend  ('MKGRD_prerotate',0)

  call PROF_rapstart('MKGRD_stretch',0)
  call MKGRD_stretch
  call PROF_rapend  ('MKGRD_stretch',0)

  call PROF_rapstart('MKGRD_shrink',0)
  call MKGRD_shrink
  call PROF_rapend  ('MKGRD_shrink',0)

  call PROF_rapstart('MKGRD_rotate',0)
  call MKGRD_rotate
  call PROF_rapend  ('MKGRD_rotate',0)

  call PROF_rapstart('MKGRD_gravcenter',0)
  call MKGRD_gravcenter
  call PROF_rapend  ('MKGRD_gravcenter',0)

  call GRD_output_hgrid( basename      = MKGRD_OUT_BASENAME, & ! [IN]
                         output_vertex = .true.,             & ! [IN]
                         io_mode       = MKGRD_OUT_io_mode   ) ! [IN]


  !---< gmtr module setup >---
  call PROF_rapstart('GMTR_setup',0)
  call GRD_makelatlon
  call GRD_scaling( RADIUS )
  call GMTR_setup
  call PROF_rapend  ('GMTR_setup',0)

  call PROF_rapend('Main_MKGRD',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
end program mkhgrid
