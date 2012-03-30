!-------------------------------------------------------------------------------
!> module OCEAN VARIABLES
!!
!! @par Description
!!          Container for oceanic variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_vars
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
  public :: OCEAN_vars_setup
  public :: OCEAN_vars_restart_read
  public :: OCEAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public, save :: SST(1,IA,JA) ! sea surface prognostics container (with HALO)

  character(len=FIO_HSHORT), public, save :: OP_NAME(1)
  character(len=FIO_HMID),   public, save :: OP_DESC(1)
  character(len=FIO_HSHORT), public, save :: OP_UNIT(1)

  data OP_NAME / 'SST' /
  data OP_DESC / 'sea surface temp.' /
  data OP_UNIT / 'K' /

  character(len=IO_SYSCHR),  public, save :: OCEAN_TYPE    = 'NONE'

  logical,                   public, save :: OCEAN_sw_sf
  logical,                   public, save :: OCEAN_sw_restart

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,                   private, save :: OCEAN_RESTART_OUTPUT       = .false.
  character(len=IO_FILECHR), public,  save :: OCEAN_RESTART_IN_BASENAME  = ''
  character(len=IO_FILECHR), private, save :: OCEAN_RESTART_OUT_BASENAME = 'restart_out'

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup oceanpheric variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN / &
       OCEAN_TYPE

    NAMELIST / PARAM_OCEAN_VARS / &
       OCEAN_RESTART_IN_BASENAME, &
       OCEAN_RESTART_OUTPUT,      &
       OCEAN_RESTART_OUT_BASENAME

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[OCEAN VARS]/Categ[OCEAN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN)

    if( IO_L ) write(IO_FID_LOG,*) '*** [OCEAN] selected components'

    if ( OCEAN_TYPE == 'FIXEDSST' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocn-Atm Interface : Fixed SST'
       OCEAN_sw_sf = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocn-Atm Interface : NONE'
       OCEAN_sw_sf = .false.
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[OCEAN VARS]/Categ[OCEAN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_VARS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_VARS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN_VARS)

    if( IO_L ) write(IO_FID_LOG,*) '*** [OCEAN] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,*) '***       | VARNAME|DESCRIPTION', &
    '                                                     [UNIT            ]'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,i3,A,A8,5(A))') &
    '*** NO.',1,'|',trim(OP_NAME(1)),'|', OP_DESC(1),'[', OP_UNIT(1),']'

    if( IO_L ) write(IO_FID_LOG,*) 'Output...'
    if ( OCEAN_RESTART_OUTPUT ) then
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output : YES'
       OCEAN_sw_restart = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '  Restart output : NO'
       OCEAN_sw_restart = .false.
    endif
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Read restart of oceanpheric variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_restart_read
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8) :: restart_ocean(1,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file (ocean) ***'

    if ( OCEAN_RESTART_IN_BASENAME /= '' ) then
       bname = OCEAN_RESTART_IN_BASENAME

       call FIO_input( restart_ocean(:,:,:), bname, 'SST', 'ZSFC', 1, 1, 1 )

       SST(1,IS:IE,JS:JE) = restart_ocean(1,1:IMAX,1:JMAX)

       ! fill IHALO & JHALO
       call COMM_vars8( SST(:,:,:), 1 )
       call COMM_wait ( SST(:,:,:), 1 )
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
    endif

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Read restart of oceanpheric variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_restart_write
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_fileio_h, only: &
       FIO_REAL8
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8) :: restart_ocean(1,IMAX,JMAX) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ocean) ***'

    write(bname,'(A,A,F15.3)') trim(OCEAN_RESTART_OUT_BASENAME), '_', NOWSEC
    desc = 'SCALE3 OCEANIC VARS.'

    restart_ocean(1,1:IMAX,1:JMAX) = SST(1,IS:IE,JS:JE)

    call FIO_output( restart_ocean(:,:,:), bname, desc, '',     &
                     'SST', OP_DESC(1), '', OP_UNIT(1),         &
                     FIO_REAL8, 'ZSFC', 1, 1, 1, NOWSEC, NOWSEC )

    return
  end subroutine OCEAN_vars_restart_write

end module mod_ocean_vars
