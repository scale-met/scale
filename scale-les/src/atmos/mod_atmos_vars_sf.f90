!-------------------------------------------------------------------------------
!> module ATMOSPHERIC Surface Variables
!!
!! @par Description
!!          Container for atmospheric surface variables
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-03-27 (H.Yashiro)  [new]
!! @li      2013-08-31 (T.Yamaura)  [mod]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_vars_sf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  use gtool_file_h, only: &
     File_HSHORT, &
     File_HMID,   &
     File_HLONG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_vars_sf_setup
  public :: ATMOS_vars_sf_fillhalo
  public :: ATMOS_vars_sf_restart_read
  public :: ATMOS_vars_sf_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
# include "scalelib.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: PREC(:,:) ! surface precipitation rate [kg/m2/s]
  real(RP), public, allocatable :: SWD (:,:) ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: LWD (:,:) ! downward long-wave radiation flux (upward positive) [W/m2]

  real(RP), public, allocatable :: SFLX_MOMZ (:,:) ! momentum z [kg/s/m2]
  real(RP), public, allocatable :: SFLX_MOMX (:,:) ! momentum x [kg/s/m2]
  real(RP), public, allocatable :: SFLX_MOMY (:,:) ! momentum y [kg/s/m2]
  real(RP), public, allocatable :: SFLX_SWU  (:,:) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: SFLX_LWU  (:,:) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: SFLX_SH   (:,:) ! sensible heat flux (upward positive) [W/m2]
  real(RP), public, allocatable :: SFLX_LH   (:,:) ! latent heat flux (upward positive) [W/m2]
  real(RP), public, allocatable :: SFLX_QVAtm(:,:) ! moisture flux for atmosphere [kg/m2/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                   private, save :: ATMOS_SF_RESTART_OUTPUT        = .false.
  character(len=IO_FILECHR), private, save :: ATMOS_SF_RESTART_IN_BASENAME   = ''
  character(len=IO_FILECHR), private, save :: ATMOS_SF_RESTART_OUT_BASENAME  = 'restart_out'
  character(len=File_HLONG), private, save :: ATMOS_SF_RESTART_OUT_TITLE     = 'SCALE-LES PROGNOSTIC VARS.'
  character(len=File_HLONG), private, save :: ATMOS_SF_RESTART_OUT_SOURCE    = 'SCALE-LES ver. '//VERSION
  character(len=File_HLONG), private, save :: ATMOS_SF_RESTART_OUT_INSTITUTE = 'AICS/RIKEN'

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup atmospheric surface variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_sf_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_SF_VARS / &
       ATMOS_SF_RESTART_IN_BASENAME, &
       ATMOS_SF_RESTART_OUTPUT,      &
       ATMOS_SF_RESTART_OUT_BASENAME

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Variables]/Categ[ATMOS]'

    allocate( PREC(IA,JA) )
    allocate( SWD (IA,JA) )
    allocate( LWD (IA,JA) )

    allocate( SFLX_MOMZ (IA,JA) )
    allocate( SFLX_MOMX (IA,JA) )
    allocate( SFLX_MOMY (IA,JA) )
    allocate( SFLX_SWU  (IA,JA) )
    allocate( SFLX_LWU  (IA,JA) )
    allocate( SFLX_SH   (IA,JA) )
    allocate( SFLX_LH   (IA,JA) )
    allocate( SFLX_QVAtm(IA,JA) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SF_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_SF_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_SF_VARS)

    PREC(:,:) = 0.0_RP
    SWD (:,:) = 0.0_RP
    LWD (:,:) = 0.0_RP

    SFLX_MOMZ (:,:) = 0.0_RP
    SFLX_MOMX (:,:) = 0.0_RP
    SFLX_MOMY (:,:) = 0.0_RP
    SFLX_SWU  (:,:) = 0.0_RP
    SFLX_LWU  (:,:) = 0.0_RP
    SFLX_SH   (:,:) = 0.0_RP
    SFLX_LH   (:,:) = 0.0_RP
    SFLX_QVAtm(:,:) = 0.0_RP

    return
  end subroutine ATMOS_vars_sf_setup

  subroutine ATMOS_vars_sf_fillhalo
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    ! fill IHALO & JHALO
    call COMM_vars8( PREC(:,:), 1 )
    call COMM_vars8( SWD (:,:), 2 )
    call COMM_vars8( LWD (:,:), 3 )
    call COMM_wait ( PREC(:,:), 1 )
    call COMM_wait ( SWD (:,:), 2 )
    call COMM_wait ( LWD (:,:), 3 )

    return
  end subroutine ATMOS_vars_sf_fillhalo

  subroutine ATMOS_vars_sf_restart_read
    implicit none
    !---------------------------------------------------------------------------
    return
  end subroutine ATMOS_vars_sf_restart_read

  subroutine ATMOS_vars_sf_restart_write
    implicit none
    !---------------------------------------------------------------------------
    return
  end subroutine ATMOS_vars_sf_restart_write

end module mod_ATMOS_vars_sf
