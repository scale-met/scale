!-------------------------------------------------------------------------------
!> module OCEAN Surface
!!
!! @par Description
!!          SST control module (tentative)
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-03-27 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_sf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_FIXEDSST_setup
  public :: OCEAN_FIXEDSST

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(8), private, save :: OCEAN_FIXEDSST_RATE = 2.D-5 ! SST change rate

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup FIXEDSST
  !-----------------------------------------------------------------------------
  subroutine OCEAN_FIXEDSST_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_ocean_vars, only: &
       OCEAN_RESTART_IN_BASENAME, &
       SST
    implicit none

    real(8) :: OCEAN_FIXEDSST_STARTSST = 290.D0  ! SST for initial state
    logical :: OCEAN_FIXEDSST_RESET    = .false. ! reset SST?

    NAMELIST / PARAM_OCEAN_FIXEDSST / &
       OCEAN_FIXEDSST_STARTSST, &
       OCEAN_FIXEDSST_RESET,    &
       OCEAN_FIXEDSST_RATE

    integer :: ierr
    integer :: i,j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FIXEDSST]/Categ[OCEAN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_FIXEDSST,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_FIXEDSST. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN_FIXEDSST)

    if ( OCEAN_RESTART_IN_BASENAME == '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** restart file for ocean is not specified.'
       if( IO_L ) write(IO_FID_LOG,*) '*** default initial SST is used.'
       OCEAN_FIXEDSST_RESET = .true.
    endif

    if ( OCEAN_FIXEDSST_RESET ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** initial SST [K]  :', OCEAN_FIXEDSST_STARTSST

       do j = 1, JA
       do i = 1, IA
          SST(1,i,j) = OCEAN_FIXEDSST_STARTSST
       enddo
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** change rate [K/s]:', OCEAN_FIXEDSST_RATE

    return
  end subroutine OCEAN_FIXEDSST_setup

  !-----------------------------------------------------------------------------
  !> SST change
  !-----------------------------------------------------------------------------
  subroutine OCEAN_FIXEDSST
    use mod_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_history, only: &
       HIST_in
    use mod_ocean_vars, only: &
       OP_DESC, &
       OP_UNIT, &
       SST
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean SST step:', &
               SST(1,IS,JS), '->', SST(1,IS,JS) + OCEAN_FIXEDSST_RATE * dt

    do j = 1, JA
    do i = 1, IA
       SST(1,i,j) = SST(1,i,j) + OCEAN_FIXEDSST_RATE * dt
    enddo
    enddo

    call HIST_in( SST(:,:,:), 'SST', OP_DESC(1), OP_UNIT(1), '2D', dt )

    return
  end subroutine OCEAN_FIXEDSST

end module mod_ocean_sf
