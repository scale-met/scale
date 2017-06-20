!-------------------------------------------------------------------------------
!> module OCEAN / Physics Slab model
!!
!! @par Description
!!          ocean physics module, slab model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_phy_slab
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_SLAB_setup
  public :: OCEAN_PHY_SLAB

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
  real(RP), private :: OCEAN_PHY_SLAB_DEPTH = 10.0_RP !< water depth of slab ocean [m]
  real(RP), private :: OCEAN_PHY_SLAB_HeatCapacity    !< heat capacity of slab ocean [J/K/m2]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_SLAB_setup( OCEAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE

    NAMELIST / PARAM_OCEAN_PHY_SLAB / &
       OCEAN_PHY_SLAB_DEPTH

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[OCEAN PHY] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_PHY_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN_PHY_SLAB)

    OCEAN_PHY_SLAB_HeatCapacity = DWATR * CL * OCEAN_PHY_SLAB_DEPTH

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Slab ocean depth [m]         : ', OCEAN_PHY_SLAB_DEPTH
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean heat capacity [J/K/m2] : ', OCEAN_PHY_SLAB_HeatCapacity

    return
  end subroutine OCEAN_PHY_SLAB_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_PHY_SLAB( &
       OCEAN_TEMP_t,    &
       OCEAN_TEMP,      &
       OCEAN_SFLX_WH,   &
       OCEAN_SFLX_prec, &
       OCEAN_SFLX_evap, &
       dt               )
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    real(RP), intent(out) :: OCEAN_TEMP_t   (IA,JA)
    real(RP), intent(in)  :: OCEAN_TEMP     (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_WH  (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_prec(IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_evap(IA,JA)
    real(DP), intent(in)  :: dt

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean physics step: Slab'

    do j = JS, JE
    do i = IS, IE
      if( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
        OCEAN_TEMP_t(i,j) = - OCEAN_SFLX_WH(i,j) / OCEAN_PHY_SLAB_HeatCapacity
      else
        OCEAN_TEMP_t(i,j) = 0.0_RP
      endif
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_SLAB

end module scale_ocean_phy_slab
