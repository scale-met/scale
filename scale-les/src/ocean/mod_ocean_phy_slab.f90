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
module mod_ocean_phy_slab
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
  public :: OCEAN_PHY_driver_setup
  public :: OCEAN_PHY_driver

  public :: OCEAN_PHY_slab_setup
  public :: OCEAN_PHY_slab

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
  real(RP), private :: OCEAN_PHY_SLAB_DEPTH    = 10.0_RP !< water depth of slab ocean [m]
  logical,  private :: OCEAN_PHY_SLAB_fixedSST = .false. !< Prevent SST change?

  real(RP), private :: OCEAN_PHY_SLAB_HeatCapacity       !< heat capacity of slab ocean [J/K/m2]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_driver_setup
    use mod_ocean_admin, only: &
       OCEAN_TYPE, &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN PHY] / Origin[SCALE-LES]'

    if ( OCEAN_sw ) then

       call OCEAN_PHY_slab_setup( OCEAN_TYPE )

       call OCEAN_PHY_driver( .true., .false. )

    endif

    return
  end subroutine OCEAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine OCEAN_PHY_driver( update_flag, history_flag )
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_history, only: &
       HIST_in
    use mod_ocean_vars, only: &
       OCEAN_TEMP,  &
       OCEAN_TEMP_t
    use mod_cpl_admin, only: &
       CPL_sw_AtmOcn
    use mod_cpl_vars, only: &
       CPL_getOcn
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    real(RP) :: FLX_heat  (IA,JA)
    real(RP) :: FLX_precip(IA,JA)
    real(RP) :: FLX_evap  (IA,JA)
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       if ( CPL_sw_AtmOcn ) then
          call CPL_getOcn( FLX_heat  (:,:), & ! [OUT]
                           FLX_precip(:,:), & ! [OUT]
                           FLX_evap  (:,:)  ) ! [OUT]
       endif

       call OCEAN_PHY_slab( OCEAN_TEMP  (:,:), & ! [IN]
                            FLX_heat    (:,:), & ! [IN]
                            FLX_precip  (:,:), & ! [IN]
                            FLX_evap    (:,:), & ! [IN]
                            OCEAN_TEMP_t(:,:)  ) ! [OUT]

       if ( history_flag ) then
          call HIST_in( OCEAN_TEMP_t(:,:), 'OCEAN_TEMP_t', 'SST tendency', 'K', dt )
       endif

    endif

    return
  end subroutine OCEAN_PHY_driver

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_slab_setup( OCEAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE

    NAMELIST / PARAM_OCEAN_SLAB / &
       OCEAN_PHY_SLAB_DEPTH,   &
       OCEAN_PHY_SLAB_fixedSST

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[OCEAN PHY] / Origin[SCALE-LES]'

    if ( OCEAN_TYPE /= 'SLAB' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx OCEAN_TYPE is not SLAB. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_SLAB)

    OCEAN_PHY_SLAB_HeatCapacity = DWATR * CL * OCEAN_PHY_SLAB_DEPTH

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Prevent SST change?          : ', OCEAN_PHY_SLAB_fixedSST
    if( IO_L ) write(IO_FID_LOG,*) '*** Slab ocean depth [m]         : ', OCEAN_PHY_SLAB_DEPTH
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean heat capacity [J/K/m2] : ', OCEAN_PHY_SLAB_HeatCapacity

    return
  end subroutine OCEAN_PHY_slab_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_PHY_slab( &
       OCEAN_TEMP,   &
       FLX_heat,     &
       FLX_precip,   &
       FLX_evap,     &
       OCEAN_TEMP_t  )
    implicit none

    real(RP), intent(in)  :: OCEAN_TEMP  (IA,JA)
    real(RP), intent(in)  :: FLX_heat    (IA,JA)
    real(RP), intent(in)  :: FLX_precip  (IA,JA)
    real(RP), intent(in)  :: FLX_evap    (IA,JA)
    real(RP), intent(out) :: OCEAN_TEMP_t(IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean step: Slab'

    if( OCEAN_PHY_SLAB_fixedSST ) return

    do j = JS, JE
    do i = IS, IE
       OCEAN_TEMP_t(i,j) = - FLX_heat(i,j) / OCEAN_PHY_SLAB_HeatCapacity
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_slab

end module mod_ocean_phy_slab
