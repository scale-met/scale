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
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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
  logical,  public, save :: ATMOS_SF_sw_restart !< output restart?

  real(RP), public, allocatable :: PREC(:,:) ! surface precipitation rate [kg/m2/s]
  real(RP), public, allocatable :: SWD (:,:) ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), public, allocatable :: LWD (:,:) ! downward long-wave radiation flux (upward positive) [W/m2]

  integer,  public, parameter :: PV_NUM = 3

  integer,  public, parameter :: I_PREC = 1
  integer,  public, parameter :: I_SWD  = 2
  integer,  public, parameter :: I_LWD  = 3

  character(len=H_SHORT), public, save :: PV_NAME(PV_NUM) !< name  of the land variables
  character(len=H_MID),   public, save :: PV_DESC(PV_NUM) !< desc. of the land variables
  character(len=H_SHORT), public, save :: PV_UNIT(PV_NUM) !< unit  of the land variables

  data PV_NAME / 'PREC', &
                 'SWD',  &
                 'LWD'   /
  data PV_DESC / 'surface precipitation rate',         &
                 'downward short-wave radiation flux', &
                 'downward long-wave radiation flux'   /
  data PV_UNIT / 'kg/m2/s', &
                 'W/m2',    &
                 'W/m2'     /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,               private, save :: ATMOS_SF_RESTART_OUTPUT        = .false.
  character(len=H_LONG), private, save :: ATMOS_SF_RESTART_IN_BASENAME   = ''
  character(len=H_LONG), private, save :: ATMOS_SF_RESTART_OUT_BASENAME  = ''
  character(len=H_MID),  private, save :: ATMOS_SF_RESTART_OUT_TITLE     = 'SCALE-LES SURFACE VARS.'
  character(len=H_MID),  private, save :: ATMOS_SF_RESTART_OUT_DTYPE     = 'DEFAULT'                 !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup atmospheric surface variables
  !-----------------------------------------------------------------------------
  subroutine ATMOS_vars_sf_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_SF_VARS / &
       ATMOS_SF_RESTART_IN_BASENAME, &
       ATMOS_SF_RESTART_OUTPUT,      &
       ATMOS_SF_RESTART_OUT_BASENAME

    integer :: ip, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[ATMOS-SF VARS]/Categ[ATMOS]'

    allocate( PREC(IA,JA) )
    allocate( SWD (IA,JA) )
    allocate( LWD (IA,JA) )

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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [ATMOS-SF] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,A32,3(A))') &
               '***       |',' VARNAME','|', 'DESCRIPTION                     ','[', 'UNIT            ',']'
    do ip = 1, PV_NUM
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,A,A32,3(A))') &
                  '*** NO.',ip,'|',trim(PV_NAME(ip)),'|', PV_DESC(ip),'[', PV_UNIT(ip),']'
    enddo

    ! restart switch
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( ATMOS_SF_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Atmos-SF restart output : YES'
       ATMOS_SF_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Atmos-SF restart output : NO'
       ATMOS_SF_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_vars_sf_setup

  subroutine ATMOS_vars_sf_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
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
    use scale_fileio, only: &
       FILEIO_read
    use scale_const, only: &
       CONST_UNDEF
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (atmos-sf) ***'

    call PROF_rapstart('FILE I NetCDF')

    if ( ATMOS_SF_RESTART_IN_BASENAME /= '' ) then

       call FILEIO_read( PREC(:,:),                                         & ! [OUT]
                         ATMOS_SF_RESTART_IN_BASENAME, 'PREC', 'XY', step=1 ) ! [IN]
       call FILEIO_read( SWD(:,:),                                          & ! [OUT]
                         ATMOS_SF_RESTART_IN_BASENAME, 'SWD',  'XY', step=1 ) ! [IN]
       call FILEIO_read( LWD(:,:),                                          & ! [OUT]
                         ATMOS_SF_RESTART_IN_BASENAME, 'LWD',  'XY', step=1 ) ! [IN]

       call ATMOS_vars_sf_fillhalo

       call ATMOS_vars_sf_total
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for atmos-sf is not specified.'

       PREC(:,:) = CONST_UNDEF
       SWD (:,:) = CONST_UNDEF
       LWD (:,:) = CONST_UNDEF
    endif

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine ATMOS_vars_sf_restart_read

  subroutine ATMOS_vars_sf_restart_write
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    character(len=H_LONG) :: bname

    integer :: n
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

    if ( ATMOS_SF_RESTART_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (atmos-sf) ***'

       bname = ''
       write(bname(1:15), '(F15.3)') NOWSEC
       do n = 1, 15
          if ( bname(n:n) == ' ' ) bname(n:n) = '0'
       end do
       write(bname,'(A,A,A)') trim(ATMOS_SF_RESTART_OUT_BASENAME), '_', trim(bname)

       call FILEIO_write( PREC(:,:),  bname,                        ATMOS_SF_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(1), PV_DESC(1), PV_UNIT(1), 'XY', ATMOS_SF_RESTART_OUT_DTYPE, .true. ) ! [IN]
       call FILEIO_write( SWD (:,:),  bname,                        ATMOS_SF_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(2), PV_DESC(2), PV_UNIT(2), 'XY', ATMOS_SF_RESTART_OUT_DTYPE, .true. ) ! [IN]
       call FILEIO_write( LWD (:,:),  bname,                        ATMOS_SF_RESTART_OUT_TITLE,        & ! [IN]
                          PV_NAME(3), PV_DESC(3), PV_UNIT(3), 'XY', ATMOS_SF_RESTART_OUT_DTYPE, .true. ) ! [IN]

    endif

    call PROF_rapend  ('FILE O NetCDF')

    return
  end subroutine ATMOS_vars_sf_restart_write

  subroutine ATMOS_vars_sf_total
    use scale_stats, only: &
       STAT_checktotal, &
       STAT_total
    implicit none

    real(RP) :: total
    !---------------------------------------------------------------------------

    if ( STAT_checktotal ) then

!       call STAT_total( total, PREC(:,:), PV_NAME(I_PREC) )
!       call STAT_total( total, SWD (:,:), PV_NAME(I_SWD)  )
!       call STAT_total( total, LWD (:,:), PV_NAME(I_LWD)  )

    endif

    return
  end subroutine ATMOS_vars_sf_total

end module mod_ATMOS_vars_sf
