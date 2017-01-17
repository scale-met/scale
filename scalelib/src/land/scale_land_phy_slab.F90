!-------------------------------------------------------------------------------
!> module LAND / Physics Slab model
!!
!! @par Description
!!          slab-type land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_phy_slab
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
  public :: LAND_PHY_SLAB_setup
  public :: LAND_PHY_SLAB

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
  logical, private :: LAND_PHY_SLAB_const = .false. ! constant condition?

  logical, private :: LAND_PHY_UPDATE_BOTTOM_TEMP  = .false. ! Is LAND_TEMP  updated in the lowest level?
  logical, private :: LAND_PHY_UPDATE_BOTTOM_WATER = .false. ! Is LAND_WATER updated in the lowest level?

  real(RP), private :: WATER_DENSCS !< Heat Capacity (rho*CS) for soil moisture [J/K/m3]

  logical, allocatable, private :: is_LND(:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_PHY_SLAB_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use scale_landuse, only: &
       LANDUSE_fact_land
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    NAMELIST / PARAM_LAND_PHY_SLAB / &
       LAND_PHY_UPDATE_BOTTOM_TEMP,   &
       LAND_PHY_UPDATE_BOTTOM_WATER

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[LAND PHY] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PHY_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_PHY_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_PHY_SLAB)

    if( LAND_TYPE == 'CONST' ) then
       LAND_PHY_SLAB_const = .true.
    else if( LAND_TYPE == 'SLAB' ) then
       LAND_PHY_SLAB_const = .false.
    else
       write(*,*) 'xxx wrong LAND_TYPE. Check!'
       call PRC_MPIstop
    end if

    WATER_DENSCS = DWATR * CL

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Update soil temperature of bottom layer? : ', LAND_PHY_UPDATE_BOTTOM_TEMP
    if( IO_L ) write(IO_FID_LOG,*) '*** Update soil moisture    of bottom layer? : ', LAND_PHY_UPDATE_BOTTOM_WATER

    ! judge to run slab land model
    allocate( is_LND(IA,JA) )

    do j = JS, JE
    do i = IS, IE
      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then
        is_LND(i,j) = .true.
      else
        is_LND(i,j) = .false.
      end if
    end do
    end do

    return
  end subroutine LAND_PHY_SLAB_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_SLAB( &
       TEMP_t,       &
       WATER_t,      &
       TEMP,         &
       WATER,        &
       WaterLimit,   &
       ThermalCond,  &
       HeatCapacity, &
       WaterDiff,    &
       SFLX_GH,      &
       SFLX_prec,    &
       SFLX_evap,    &
       CDZ,          &
       dt            )
    use scale_grid_index
    use scale_const, only: &
       DWATR => CONST_DWATR
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    ! arguments
    real(RP), intent(out) :: TEMP_t      (LKMAX,IA,JA)
    real(RP), intent(out) :: WATER_t     (LKMAX,IA,JA)

    real(RP), intent(in)  :: TEMP        (LKMAX,IA,JA)
    real(RP), intent(in)  :: WATER       (LKMAX,IA,JA)
    real(RP), intent(in)  :: WaterLimit  (IA,JA)
    real(RP), intent(in)  :: ThermalCond (IA,JA)
    real(RP), intent(in)  :: HeatCapacity(IA,JA)
    real(RP), intent(in)  :: WaterDiff   (IA,JA)
    real(RP), intent(in)  :: SFLX_GH     (IA,JA)
    real(RP), intent(in)  :: SFLX_prec   (IA,JA)
    real(RP), intent(in)  :: SFLX_evap   (IA,JA)
    real(RP), intent(in)  :: CDZ         (LKMAX)
    real(DP), intent(in)  :: dt

    ! work
    real(RP) :: TEMP1 (LKMAX,IA,JA)
    real(RP) :: WATER1(LKMAX,IA,JA)
!    real(RP) :: RUNOFF(IA,JA)
    real(RP) :: SOIL_DENSCS(IA,JA)

    real(RP) :: U   (LKMAX-1,IA,JA)
    real(RP) :: M   (LKMAX-1,IA,JA)
    real(RP) :: L   (LKMAX-1,IA,JA)
    real(RP) :: Vin (LKMAX-1,IA,JA)
    real(RP) :: Vout(LKMAX-1,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land  physics step: Slab'

    ! Solve diffusion of soil moisture (tridiagonal matrix)
    do j = JS, JE
    do i = IS, IE
       L(LKS,i,j) = 0.0_RP
       U(LKS,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
       M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = LKS+1, LKE-1
       L(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
       U(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
       M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Vin(LKS  ,i,j) = WATER(LKS,i,j) &
                      + ( SFLX_prec(i,j) - SFLX_evap(i,j) ) / ( CDZ(LKS) * DWATR ) * dt ! input from atmosphere

       do k = LKS+1, LKE-2
          Vin(k,i,j) = WATER(k,i,j)
       enddo

       Vin(LKE-1,i,j) = WATER(LKE-1,i,j) &
                      - U(LKE-1,i,j) * WATER(LKE,i,j)
    enddo
    enddo

    call MATRIX_SOLVER_tridiagonal( LKMAX-1,     & ! [IN]
                                    IA, IS, IE,  & ! [IN]
                                    JA, JS, JE,  & ! [IN]
                                    U   (:,:,:), & ! [IN]
                                    M   (:,:,:), & ! [IN]
                                    L   (:,:,:), & ! [IN]
                                    Vin (:,:,:), & ! [IN]
                                    Vout(:,:,:)  ) ! [OUT]

    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE-1
       WATER1(k,i,j) = Vout(k,i,j)
    enddo
    enddo
    enddo

    ! lowest layer treatment
    if ( LAND_PHY_UPDATE_BOTTOM_WATER ) then
       do j = JS, JE
       do i = IS, IE
          WATER1(LKE,i,j) = Vout(LKE-1,i,j)
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          WATER1(LKE,i,j) = WATER(LKE,i,j)
       enddo
       enddo
    endif



    ! runoff of soil moisture (vertical sum)
    do j = JS, JE
    do i = IS, IE
!       RUNOFF(i,j) = 0.0_RP
       do k = LKS, LKE
!          RUNOFF(i,j)   = RUNOFF(i,j) &
!                             + max( WATER1(k,i,j) - WaterLimit(i,j), 0.0_RP ) * CDZ(k) * DWATR

          WATER1(k,i,j) = min( WATER1(k,i,j), WaterLimit(i,j) )
       enddo
    enddo
    enddo



    ! Solve diffusion of soil temperature (tridiagonal matrix)
    do j = JS, JE
    do i = IS, IE
       SOIL_DENSCS(i,j) = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       L(LKS,i,j) = 0.0_RP
       U(LKS,i,j) = -2.0_RP * ThermalCond(i,j) / ( SOIL_DENSCS(i,j) + WATER(LKS,i,j)*WATER_DENSCS ) &
                    / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
       M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = LKS+1, LKE-1
       L(k,i,j) = -2.0_RP * ThermalCond(i,j) / ( SOIL_DENSCS(i,j) + WATER(k,i,j)*WATER_DENSCS ) &
                  / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
       U(k,i,j) = -2.0_RP * ThermalCond(i,j) / ( SOIL_DENSCS(i,j) + WATER(k,i,j)*WATER_DENSCS ) &
                  / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
       M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Vin(LKS  ,i,j) = TEMP(LKS,i,j) &
                      - SFLX_GH(i,j) / ( SOIL_DENSCS(i,j) + WATER(LKS,i,j)*WATER_DENSCS ) / CDZ(LKS) * dt ! input from atmosphere

       do k = LKS+1, LKE-2
          Vin(k,i,j) = TEMP(k,i,j)
       enddo

       Vin(LKE-1,i,j) = TEMP(LKE-1,i,j) &
                      - U(LKE-1,i,j) * TEMP(LKE,i,j)
    enddo
    enddo

    call MATRIX_SOLVER_tridiagonal( LKMAX-1,     & ! [IN]
                                    IA, IS, IE,  & ! [IN]
                                    JA, JS, JE,  & ! [IN]
                                    U   (:,:,:), & ! [IN]
                                    M   (:,:,:), & ! [IN]
                                    L   (:,:,:), & ! [IN]
                                    Vin (:,:,:), & ! [IN]
                                    Vout(:,:,:)  ) ! [OUT]

    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE-1
       TEMP1(k,i,j) = Vout(k,i,j)
    enddo
    enddo
    enddo

    ! lowest layer treatment
    if ( LAND_PHY_UPDATE_BOTTOM_TEMP ) then
       do j = JS, JE
       do i = IS, IE
          TEMP1(LKE,i,j) = Vout(LKE-1,i,j)
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          TEMP1(LKE,i,j) = TEMP(LKE,i,j)
       enddo
       enddo
    endif


    ! calculate tendency
    if( LAND_PHY_SLAB_const ) then
       do j = JS, JE
       do i = IS, IE
       do k = LKS, LKE
          TEMP_t (k,i,j) = 0.0_RP
          WATER_t(k,i,j) = 0.0_RP
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          if( is_LND(i,j) ) then
             do k = LKS, LKE
                TEMP_t (k,i,j) = ( TEMP1 (k,i,j) - TEMP (k,i,j) ) / dt
                WATER_t(k,i,j) = ( WATER1(k,i,j) - WATER(k,i,j) ) / dt
             enddo
          else
             do k = LKS, LKE
                TEMP_t (k,i,j) = 0.0_RP
                WATER_t(k,i,j) = 0.0_RP
             enddo
          endif
       enddo
       enddo
    endif

    return
  end subroutine LAND_PHY_SLAB

end module scale_land_phy_slab
