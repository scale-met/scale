!-------------------------------------------------------------------------------
!> module LAND / Physics Matsiro model
!!
!! @par Description
!!          matsiro-type land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_phy_matsiro
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_matsiro_setup
  public :: LAND_PHY_matsiro

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
  subroutine LAND_PHY_matsiro_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    logical :: dummy

    NAMELIST / PARAM_LAND_MATSIRO / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[MATSIRO] / Categ[LAND PHY] / Origin[SCALE-LES]'

    if ( LAND_TYPE /= 'MATSIRO' ) then
       write(*,*) 'xxx LAND_TYPE is not MATSIRO. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_MATSIRO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_MATSIRO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_MATSIRO)

    return
  end subroutine LAND_PHY_matsiro_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_matsiro( &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_WaterLimit,   &
       LAND_ThermalCond,  &
       LAND_HeatCapacity, &
       LAND_WaterDiff,    &
       FLX_heat,          &
       FLX_precip,        &
       FLX_evap,          &
       CDZ,               &
       LAND_TEMP_t,       &
       LAND_WATER_t       )
    use scale_grid_index
    implicit none

    real(RP), intent(in)  :: LAND_TEMP        (LKMAX,IA,JA)
    real(RP), intent(in)  :: LAND_WATER       (LKMAX,IA,JA)
    real(RP), intent(in)  :: LAND_WaterLimit  (IA,JA)
    real(RP), intent(in)  :: LAND_ThermalCond (IA,JA)
    real(RP), intent(in)  :: LAND_HeatCapacity(IA,JA)
    real(RP), intent(in)  :: LAND_WaterDiff   (IA,JA)
    real(RP), intent(in)  :: FLX_heat         (IA,JA)
    real(RP), intent(in)  :: FLX_precip       (IA,JA)
    real(RP), intent(in)  :: FLX_evap         (IA,JA)
    real(RP), intent(in)  :: CDZ              (LKMAX)
    real(RP), intent(out) :: LAND_TEMP_t      (LKMAX,IA,JA)
    real(RP), intent(out) :: LAND_WATER_t     (LKMAX,IA,JA)

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land step: Matsiro'

    LAND_TEMP_t (:,:,:) = 0.0_RP
    LAND_WATER_t(:,:,:) = 0.0_RP

    return
  end subroutine LAND_PHY_matsiro

end module scale_land_phy_matsiro
