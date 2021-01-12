!-------------------------------------------------------------------------------
!> module TOPOGRAPHY
!!
!! @par Description
!!          Topography module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_gm_topography
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TOPOGRAPHY_setup
  public :: TOPOGRAPHY_fillhalo
  public :: TOPOGRAPHY_write


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,  public :: TOPOGRAPHY_exist    = .false. !< topography exists?

  logical,  public :: TOPOGRAPHY_IN_IDEAL = .false. !< make ideal topo on the fly?

  real(RP), public, allocatable :: TOPOGRAPHY_Zsfc   (:,:,:,:) !< absolute ground height [m]
  real(RP), public, allocatable :: TOPOGRAPHY_Zsfc_pl(:,:,:,:) !< absolute ground height [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPOGRAPHY_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: TOPOGRAPHY_IN_BASENAME  = ''                    !< basename of the input  file
  character(len=H_LONG),  private :: TOPOGRAPHY_OUT_BASENAME = ''                    !< basename of the output file
  character(len=H_MID),   private :: TOPOGRAPHY_OUT_TITLE    = 'SCALE-GM TOPOGRAPHY' !< title    of the output file
  character(len=H_SHORT), private :: TOPOGRAPHY_OUT_DTYPE    = 'DEFAULT'             !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TOPOGRAPHY_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_TOPOGRAPHY / &
       TOPOGRAPHY_IN_IDEAL,     &
       TOPOGRAPHY_IN_BASENAME,  &
       TOPOGRAPHY_OUT_BASENAME, &
       TOPOGRAPHY_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("TOPOGRAPHY_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPOGRAPHY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TOPOGRAPHY_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TOPOGRAPHY_setup",*) 'Not appropriate names in namelist PARAM_TOPOGRAPHY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TOPOGRAPHY)

    allocate( TOPOGRAPHY_Zsfc   (ADM_gall   ,ADM_KNONE,ADM_lall   ,1) )
    allocate( TOPOGRAPHY_Zsfc_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1) )
    TOPOGRAPHY_Zsfc   (:,:,:,:) = 0.0_RP
    TOPOGRAPHY_Zsfc_pl(:,:,:,:) = 0.0_RP

    ! read from file
    call TOPOGRAPHY_read

    return
  end subroutine TOPOGRAPHY_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine TOPOGRAPHY_fillhalo
    use scale_comm_icoA, only: &
       COMM_var
    implicit none
    !---------------------------------------------------------------------------

    call COMM_var( TOPOGRAPHY_Zsfc, TOPOGRAPHY_Zsfc_pl, ADM_KNONE, 1 )

    return
  end subroutine TOPOGRAPHY_fillhalo

  !-----------------------------------------------------------------------------
  !> Read topography
  subroutine TOPOGRAPHY_read
    use mod_fio, only: &
       FIO_input
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("TOPOGRAPHY_read",*) 'Input topography file '

    if ( TOPOGRAPHY_IN_BASENAME /= '' ) then

       call FIO_input(TOPOGRAPHY_Zsfc(:,:,:,1),TOPOGRAPHY_IN_BASENAME,'topo','ZSSFC1',1,1,1)

       call TOPOGRAPHY_fillhalo

       TOPOGRAPHY_exist = .true.

    else
       LOG_INFO_CONT(*) 'topography file is not specified.'

       TOPOGRAPHY_exist = .false.
    endif

    return
  end subroutine TOPOGRAPHY_read

  !-----------------------------------------------------------------------------
  !> Write topography
  subroutine TOPOGRAPHY_write
    use scale_prc, only: &
       PRC_abort
    use mod_io_param, only: &
       IO_REAL4, &
       IO_REAL8
    use mod_fio, only: &
       FIO_output
    implicit none

    integer :: dtype
    !---------------------------------------------------------------------------

    if ( TOPOGRAPHY_OUT_BASENAME /= '' ) then

       LOG_NEWLINE
       LOG_INFO("TOPOGRAPHY_write",*) 'Output topography file '

       call TOPOGRAPHY_fillhalo

       ! dtype is used to define the data type of axis variables in file
       if    ( TOPOGRAPHY_OUT_DTYPE == 'REAL8' ) then
          dtype = IO_REAL8
       elseif( TOPOGRAPHY_OUT_DTYPE == 'REAL4' ) then
          dtype = IO_REAL4
       else
          if    ( RP == 8 ) then
             dtype = IO_REAL8
          elseif( RP == 4 ) then
             dtype = IO_REAL4
          else
             LOG_ERROR("FILE_CARTESC_create",*) 'unsupported data type. Check!', trim(TOPOGRAPHY_OUT_DTYPE)
             call PRC_abort
          endif
       endif

       call FIO_output( TOPOGRAPHY_Zsfc(:,:,:,1),                          & ! [IN]
                        TOPOGRAPHY_OUT_BASENAME, TOPOGRAPHY_OUT_TITLE, '', & ! [IN]
                       'topo', 'topography', '',                           & ! [IN]
                       'm', dtype, 'ZSSFC1', 1, 1, 1, 0.0_DP, 0.0_DP       ) ! [IN]
    endif

    return
  end subroutine TOPOGRAPHY_write

end module mod_gm_topography
