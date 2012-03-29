!-------------------------------------------------------------------------------
!> module Atmospheric Surface Variables
!!
!! @par Description
!!          Container for atmospheric surface variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-27 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_vars_sf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR,  &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_vars_sf_setup
  public :: ATMOS_vars_sf_restart_read
  public :: ATMOS_vars_sf_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, save :: SFLX_MOMZ(IA,JA) ! momentum z [kg/s/m2]
  real(8), public, save :: SFLX_MOMX(IA,JA) ! momentum x [kg/s/m2]
  real(8), public, save :: SFLX_MOMY(IA,JA) ! momentum y [kg/s/m2]
  real(8), public, save :: SFLX_POTT(IA,JA) ! POTT [K]
  real(8), public, save :: SFLX_QV  (IA,JA) ! tracer mixing ratio [kg/kg]
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                   private, save :: ATMOS_SF_RESTART_OUTPUT       = .false.
  character(len=IO_FILECHR), private, save :: ATMOS_SF_RESTART_IN_BASENAME  = 'restart_in'
  character(len=IO_FILECHR), private, save :: ATMOS_SF_RESTART_OUT_BASENAME = 'restart_out'

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
    integer :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Variables]/Categ[ATMOS]'

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

    return
  end subroutine ATMOS_vars_sf_setup

  !-----------------------------------------------------------------------------
  !> Read restart of atmospheric surface variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_sf_restart_read
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=8)          :: lname

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (atmos) ***'

    bname = ATMOS_SF_RESTART_IN_BASENAME
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

!    call FIO_input( restart_atmos(:,:,:), bname, 'DENS', lname, 1, KMAX, 1 )
!    DENS(KS:KE,IS:IE,JS:JE) = restart_atmos(1:KMAX,1:IMAX,1:JMAX)

    ! fill IHALO & JHALO
!    call COMM_vars8( DENS(:,:,:), 1 )
!    call COMM_wait ( RHOT(:,:,:), 5 )

    return
  end subroutine ATMOS_vars_sf_restart_read

  !-----------------------------------------------------------------------------
  !> Write restart of atmospheric surface variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_sf_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_fileio_h, only: &
       FIO_REAL8
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8) :: restart_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    character(len=8)          :: lname

    integer :: iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos) ***'

    write(bname,'(A,A,F15.3)') trim(ATMOS_SF_RESTART_OUT_BASENAME), '_', NOWSEC
    desc  = 'SCALE3 PROGNOSTIC VARS.'
    write(lname,'(A,I4.4)') 'ZDEF', KMAX

    return
  end subroutine ATMOS_vars_sf_restart_write

end module mod_ATMOS_vars_sf
