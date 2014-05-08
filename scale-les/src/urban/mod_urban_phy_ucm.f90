!-------------------------------------------------------------------------------
!> module URBAN / Physics Urban Canopy Model (UCM)
!!
!! @par Description
!!          UCM physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_urban_phy_ucm
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
  public :: URBAN_PHY_driver_setup
  public :: URBAN_PHY_driver_first
  public :: URBAN_PHY_driver_final

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_PHY_driver_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_urban_vars, only: &
       URBAN_TYPE_PHY
    implicit none

    logical  :: dummy

    NAMELIST / PARAM_URBAN_UCM / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[UCM]/Categ[URBAN]'

    if ( URBAN_TYPE_PHY /= 'UCM' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx URBAN_TYPE_PHY is not UCM. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_UCM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_UCM. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_URBAN_UCM)

    return
  end subroutine URBAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for urban submodel
  subroutine URBAN_PHY_driver_first
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    implicit none

    ! work
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Urban step: UCM'

    return
  end subroutine URBAN_PHY_driver_first

  subroutine URBAN_PHY_driver_final
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine URBAN_PHY_driver_final

end module mod_urban_phy_ucm
