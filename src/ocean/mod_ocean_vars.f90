!-------------------------------------------------------------------------------
!> module OCEAN VARIABLES
!!
!! @par Description
!!          Container for oceanic variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_SYSCHR, &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID,   &
     FIO_REAL8
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
  public :: OCEAN_vars_put
  public :: OCEAN_vars_get
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,                   public,              save :: O_VA      ! Number of Tracers + 5
  character(len=FIO_HSHORT), public, allocatable, save :: O_NAME(:)
  character(len=FIO_HMID),   public, allocatable, save :: O_DESC(:)
  character(len=FIO_HSHORT), public, allocatable, save :: O_UNIT(:)

  character(len=IO_SYSCHR),  public, save :: OCEAN_TYPE    = 'NONE'

  real(8), public, allocatable, save :: ocean_var(:,:,:,:)      !> prognostics container (with HALO)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private, save :: OCEAN_RESTART_IN_BASENAME      = 'restart_in'
  character(len=IO_FILECHR), private, save :: OCEAN_RESTART_OUT_BASENAME     = 'restart_out'

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup oceanpheric variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF8
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA
    implicit none

    NAMELIST / PARAM_OCEAN / &
       OCEAN_TYPE

    real(8) :: OCEAN_SST = 290.D0

    NAMELIST / PARAM_OCEAN_VARS / &
       OCEAN_RESTART_IN_BASENAME,  &
       OCEAN_RESTART_OUT_BASENAME, &
       OCEAN_SST

    integer :: ierr
    integer :: iv
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [OCEAN] selected components'

    if ( OCEAN_TYPE == 'FIXEDSST' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocn-Atm Interface : Fixed SST'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocn-Atm Interface : NONE'
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

    O_VA = 1
    allocate( O_NAME(O_VA) )
    allocate( O_DESC(O_VA) )
    allocate( O_UNIT(O_VA) )

    O_NAME(1) = 'SST'
    O_DESC(1) = 'sea surface temp.'
    O_UNIT(1) = 'K'

    if( IO_L ) write(IO_FID_LOG,*) '*** [OCEAN] prognostic variables'
    if( IO_L ) write(IO_FID_LOG,*) &
    '***                 : VARNAME         , ', &
    'DESCRIPTION                                                     [UNIT            ]'
    do iv = 1, O_VA
       if( IO_L ) write(IO_FID_LOG,*) '*** NO.',iv,": ",O_NAME(iv),", ", O_DESC(iv),"[", O_UNIT(iv),"]"
    enddo

    allocate( ocean_var(IA,JA,1,O_VA) ); ocean_var(:,:,:,:) = CONST_UNDEF8

    ! tentative: put contstant value
    ocean_var(:,:,:,:) = OCEAN_SST

    return
  end subroutine OCEAN_vars_setup

  !-----------------------------------------------------------------------------
  !> Read restart of oceanpheric variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_restart_read
    use mod_grid, only : &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE
    use mod_comm, only: &
       COMM_vars, &
       COMM_stats
    use mod_fileio, only: &
       FIO_input
    implicit none

    real(8), allocatable :: restart_ocean(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname

    integer :: iv
    !---------------------------------------------------------------------------

    allocate( restart_ocean(IMAX,JMAX,1) )

    bname = OCEAN_RESTART_IN_BASENAME

    do iv = 1, O_VA
       call FIO_input( restart_ocean(:,:,:), bname, O_NAME(iv), 'ZSFC', 1, 1, 1 )

       ocean_var(IS:IE,JS:JE,1,iv) = restart_ocean(1:IMAX,1:JMAX,1)
    enddo

    deallocate( restart_ocean )

    ! fill IHALO & JHALO
    call COMM_vars( ocean_var(:,:,:,:) )

    call COMM_stats( ocean_var(:,:,:,:), O_NAME(:) )

    return
  end subroutine OCEAN_vars_restart_read

  !-----------------------------------------------------------------------------
  !> Read restart of oceanpheric variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_restart_write
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_grid, only : &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE
    use mod_comm, only: &
       COMM_vars, &
       COMM_stats
    use mod_fileio, only: &
       FIO_output
    implicit none

    real(8), allocatable :: restart_ocean(:,:,:) !> restart file (no HALO)

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc

    integer :: iv
    !---------------------------------------------------------------------------

    allocate( restart_ocean(IMAX,JMAX,1) )

    call COMM_stats( ocean_var(:,:,:,:), O_NAME(:) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file (ocean) ***'

    write(bname,'(A,A,F15.3)') trim(OCEAN_RESTART_OUT_BASENAME), '_', NOWSEC
    desc = 'SCALE3 OCEANIC VARS.'

    do iv = 1, O_VA
       restart_ocean(1:IMAX,1:JMAX,1) = ocean_var(IS:IE,JS:JE,1,iv)

       call FIO_output( restart_ocean(:,:,:), bname, desc, '',       &
                        O_NAME(iv), O_DESC(iv), '', O_UNIT(iv),      &
                        FIO_REAL8, 'ZSFC', 1, 1, 1, NOWSEC, NOWSEC )
    enddo

    deallocate( restart_ocean )

    return
  end subroutine OCEAN_vars_restart_write

  !-----------------------------------------------------------------------------
  !> Put and Communicate prognostic variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_put( &
       sst )
    use mod_grid, only: &
       IA   => GRID_IA, &
       JA   => GRID_JA
    use mod_comm, only: &
       COMM_vars
    implicit none

    real(8), intent(in) :: sst(IA,JA,1)
    !---------------------------------------------------------------------------

    ocean_var(:,:,:,1) = sst(:,:,:)

    ! fill IHALO & JHALO
    call COMM_vars( ocean_var(:,:,:,:) )

    return
  end subroutine OCEAN_vars_put

  !-----------------------------------------------------------------------------
  !> Get prognostic variables
  !-----------------------------------------------------------------------------
  subroutine OCEAN_vars_get( &
       sst )
    use mod_grid, only: &
       IA => GRID_IA, &
       JA => GRID_JA
    implicit none

    real(8), intent(out) :: sst(IA,JA,1)
    !---------------------------------------------------------------------------

    sst(:,:,:) = ocean_var(:,:,:,1)

    return
  end subroutine OCEAN_vars_get

end module mod_ocean_vars
