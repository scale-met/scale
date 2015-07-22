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
  real(RP), private :: OCEAN_PHY_SLAB_DEPTH    = 10.0_RP !< water depth of slab ocean [m]
  logical,  private :: OCEAN_PHY_SLAB_fixedSST = .false. !< Prevent SST change?
  real(RP), private :: OCEAN_PHY_SLAB_HeatCapacity       !< heat capacity of slab ocean [J/K/m2]

  logical, allocatable, private :: is_OCN(:,:)

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
    use scale_landuse, only: &
       LANDUSE_fact_ocean
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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_PHY_SLAB)

    if( OCEAN_TYPE == 'CONST' ) then
       OCEAN_PHY_SLAB_fixedSST = .true.
    else if( OCEAN_TYPE == 'SLAB' ) then
       OCEAN_PHY_SLAB_fixedSST = .false.
    else
       write(*,*) 'xxx wrong OCEAN_TYPE. Check!'
       call PRC_MPIstop
    end if

    OCEAN_PHY_SLAB_HeatCapacity = DWATR * CL * OCEAN_PHY_SLAB_DEPTH

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Prevent SST change?          : ', OCEAN_PHY_SLAB_fixedSST
    if( IO_L ) write(IO_FID_LOG,*) '*** Slab ocean depth [m]         : ', OCEAN_PHY_SLAB_DEPTH
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean heat capacity [J/K/m2] : ', OCEAN_PHY_SLAB_HeatCapacity

    ! judge to run slab ocean model
    allocate( is_OCN(IA,JA) )

    do j = JS, JE
    do i = IS, IE
      if( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
        is_OCN(i,j) = .true.
      else
        is_OCN(i,j) = .false.
      end if
    end do
    end do

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
    use scale_grid_index
    implicit none

    real(RP), intent(out) :: OCEAN_TEMP_t   (IA,JA)
    real(RP), intent(in)  :: OCEAN_TEMP     (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_WH  (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_prec(IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_evap(IA,JA)
    real(RP), intent(in)  :: dt

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean step: Slab'

    if( OCEAN_PHY_SLAB_fixedSST )then
       do j = JS, JE
       do i = IS, IE
          OCEAN_TEMP_t(i,j) = 0.0_RP
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          if( is_OCN(i,j) ) then
             OCEAN_TEMP_t(i,j) = - OCEAN_SFLX_WH(i,j) / OCEAN_PHY_SLAB_HeatCapacity
          else
             OCEAN_TEMP_t(i,j) = 0.0_RP
          endif
       enddo
       enddo
    end if

    return
  end subroutine OCEAN_PHY_SLAB

end module scale_ocean_phy_slab
