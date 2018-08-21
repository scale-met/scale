!-------------------------------------------------------------------------------
!> Program mkllmap
!!
!! @par Description
!!          Making remapping coefficient between lat-lon and icosahedral grid
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program mkllmap
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
     PRC_abort,            &
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
     GRD_setup
  use mod_latlon, only: &
     LATLON_setup, &
     LATLON_ico_setup
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

  character(len=H_LONG) :: output_dir = './'

  namelist /mkllmap_param/ &
     output_dir

  integer :: ierr
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

  !--- < cnst module setup > ---
  call CONST_setup

  ! setup horizontal/vertical grid coordinates (icosahedral,idealized)
  call ATMOS_GRID_icoA_INDEX_setup

  !---< I/O module setup >---
  call FIO_setup

  !--- < comm module setup > ---
  call COMM_setup

  !--- < grid module setup > ---
  call GRD_setup

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_MKLLMAP',0)

  !--- read parameters
  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '+++ Program[mkllmap]/Category[tool]'
  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=MKLLMAP_PARAM,iostat=ierr)
  if ( ierr < 0 ) then
     if( IO_L ) write(IO_FID_LOG,*) '*** MKLLMAP_PARAM is not specified. use default.'
  elseif( ierr > 0 ) then
     write(*,*) 'xxx Not appropriate names in namelist MKLLMAP_PARAM. STOP.'
     call PRC_abort
  endif
  if( IO_NML ) write(IO_FID_NML,nml=MKLLMAP_PARAM)

  call LATLON_ico_setup

  call LATLON_setup( output_dir )

  call PROF_rapend('Main_MKLLMAP',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

end program mkllmap
