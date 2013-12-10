!-------------------------------------------------------------------------------
!> module OCEAN / Physics Fixed-SST
!!
!! @par Description
!!          ocean physics module, fixed SST
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_phy
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_setup
  public :: OCEAN_PHY

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
  real(RP), private, save :: OCEAN_FIXEDSST_RATE = 2.E-5_RP !< SST change rate

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_ocean_vars, only: &
       OCEAN_RESTART_IN_BASENAME, &
       OCEAN_TYPE_PHY, &
       SST
    implicit none

    real(RP) :: OCEAN_FIXEDSST_STARTSST = 290.0_RP !< SST for initial state
    logical  :: OCEAN_FIXEDSST_RESET    = .false.   !< reset SST?

    NAMELIST / PARAM_OCEAN_FIXEDSST / &
       OCEAN_FIXEDSST_STARTSST, &
       OCEAN_FIXEDSST_RESET,    &
       OCEAN_FIXEDSST_RATE

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FIXEDSST]/Categ[OCEAN]'

    if ( OCEAN_TYPE_PHY /= 'FIXEDSST' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx OCEAN_TYPE_PHY is not FIXEDSST. Check!'
       call PRC_MPIstop
    endif

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
          SST(i,j) = OCEAN_FIXEDSST_STARTSST
       enddo
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** change rate [K/s]:', OCEAN_FIXEDSST_RATE

    return
  end subroutine OCEAN_PHY_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for ocean submodel
  subroutine OCEAN_PHY
    use mod_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_ocean_vars, only: &
       SST
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean SST step:', &
               SST(IS,JS), '->', SST(IS,JS) + OCEAN_FIXEDSST_RATE * dt

    do j = 1, JA
    do i = 1, IA
       SST(i,j) = SST(i,j) + OCEAN_FIXEDSST_RATE * dt
    enddo
    enddo

    return
  end subroutine OCEAN_PHY

end module mod_ocean_phy
