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
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    NAMELIST / PARAM_LAND_SLAB / &
       LAND_PHY_UPDATE_BOTTOM_TEMP,   &
       LAND_PHY_UPDATE_BOTTOM_WATER

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[BUCKET] / Categ[LAND PHY] / Origin[SCALE-LES]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_SLAB)

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

    return
  end subroutine LAND_PHY_SLAB_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_SLAB( &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_WaterLimit,   &
       LAND_ThermalCond,  &
       LAND_HeatCapacity, &
       LAND_WaterDiff,    &
       LAND_SFLX_GH,      &
       LAND_SFLX_prec,    &
       LAND_SFLX_evap,    &
       CDZ,               &
       dt                 )
    use scale_grid_index
    use scale_const, only: &
       DWATR => CONST_DWATR
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    ! arguments
    real(RP), intent(out) :: LAND_TEMP_t      (LKMAX,IA,JA)
    real(RP), intent(out) :: LAND_WATER_t     (LKMAX,IA,JA)

    real(RP), intent(in)  :: LAND_TEMP        (LKMAX,IA,JA)
    real(RP), intent(in)  :: LAND_WATER       (LKMAX,IA,JA)
    real(RP), intent(in)  :: LAND_WaterLimit  (IA,JA)
    real(RP), intent(in)  :: LAND_ThermalCond (IA,JA)
    real(RP), intent(in)  :: LAND_HeatCapacity(IA,JA)
    real(RP), intent(in)  :: LAND_WaterDiff   (IA,JA)
    real(RP), intent(in)  :: LAND_SFLX_GH     (IA,JA)
    real(RP), intent(in)  :: LAND_SFLX_prec   (IA,JA)
    real(RP), intent(in)  :: LAND_SFLX_evap   (IA,JA)
    real(RP), intent(in)  :: CDZ              (LKMAX)
    real(RP), intent(in)  :: dt

    ! work
    real(RP) :: LAND_TEMP1 (LKMAX,IA,JA)
    real(RP) :: LAND_WATER1(LKMAX,IA,JA)
!    real(RP) :: LAND_RUNOFF(IA,JA)
    real(RP) :: SOIL_DENSCS(IA,JA)

    real(RP) :: U   (LKMAX-1,IA,JA)
    real(RP) :: M   (LKMAX-1,IA,JA)
    real(RP) :: L   (LKMAX-1,IA,JA)
    real(RP) :: Vin (LKMAX-1,IA,JA)
    real(RP) :: Vout(LKMAX-1,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land step: Slab'

    ! Solve diffusion of soil moisture (tridiagonal matrix)
    do j = JS, JE
    do i = IS, IE
       L(LKS,i,j) = 0.0_RP
       U(LKS,i,j) = -2.0_RP * LAND_WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
       M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = LKS+1, LKE-1
       L(k,i,j) = -2.0_RP * LAND_WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
       U(k,i,j) = -2.0_RP * LAND_WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
       M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Vin(LKS  ,i,j) = LAND_WATER(LKS,i,j) &
                      + ( LAND_SFLX_prec(i,j) - LAND_SFLX_evap(i,j) ) / ( CDZ(LKS) * DWATR ) * dt ! input from atmosphere

       do k = LKS+1, LKE-2
          Vin(k,i,j) = LAND_WATER(k,i,j)
       enddo

       Vin(LKE-1,i,j) = LAND_WATER(LKE-1,i,j) &
                      - U(LKE-1,i,j) * LAND_WATER(LKE,i,j)
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
       LAND_WATER1(k,i,j) = Vout(k,i,j)
    enddo
    enddo
    enddo

    ! lowest layer treatment
    if ( LAND_PHY_UPDATE_BOTTOM_WATER ) then
       do j = JS, JE
       do i = IS, IE
          LAND_WATER1(LKE,i,j) = Vout(LKE-1,i,j)
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          LAND_WATER1(LKE,i,j) = LAND_WATER(LKE,i,j)
       enddo
       enddo
    endif



    ! runoff of soil moisture (vertical sum)
    do j = JS, JE
    do i = IS, IE
!       LAND_RUNOFF(i,j) = 0.0_RP
       do k = LKS, LKE
!          LAND_RUNOFF(i,j)   = LAND_RUNOFF(i,j) &
!                             + max( LAND_WATER1(k,i,j)-LAND_WaterLimit(i,j), 0.0_RP ) * CDZ(k) * DWATR

          LAND_WATER1(k,i,j) = min( LAND_WATER1(k,i,j), LAND_WaterLimit(i,j) )
       enddo
    enddo
    enddo



    ! Solve diffusion of soil temperature (tridiagonal matrix)
    do j = JS, JE
    do i = IS, IE
       SOIL_DENSCS(i,j) = ( 1.0_RP - LAND_WaterLimit(i,j) ) * LAND_HeatCapacity(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       L(LKS,i,j) = 0.0_RP
       U(LKS,i,j) = -2.0_RP * LAND_ThermalCond(i,j) / ( SOIL_DENSCS(i,j) + LAND_WATER(LKS,i,j)*WATER_DENSCS ) &
                    / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
       M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = LKS+1, LKE-1
       L(k,i,j) = -2.0_RP * LAND_ThermalCond(i,j) / ( SOIL_DENSCS(i,j) + LAND_WATER(k,i,j)*WATER_DENSCS ) &
                  / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
       U(k,i,j) = -2.0_RP * LAND_ThermalCond(i,j) / ( SOIL_DENSCS(i,j) + LAND_WATER(k,i,j)*WATER_DENSCS ) &
                  / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
       M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Vin(LKS  ,i,j) = LAND_TEMP(LKS,i,j) &
                      - LAND_SFLX_GH(i,j) / ( SOIL_DENSCS(i,j) + LAND_WATER(LKS,i,j)*WATER_DENSCS ) / CDZ(LKS) * dt ! input from atmosphere

       do k = LKS+1, LKE-2
          Vin(k,i,j) = LAND_TEMP(k,i,j)
       enddo

       Vin(LKE-1,i,j) = LAND_TEMP(LKE-1,i,j) &
                      - U(LKE-1,i,j) * LAND_TEMP(LKE,i,j)
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
       LAND_TEMP1(k,i,j) = Vout(k,i,j)
    enddo
    enddo
    enddo

    ! lowest layer treatment
    if ( LAND_PHY_UPDATE_BOTTOM_TEMP ) then
       do j = JS, JE
       do i = IS, IE
          LAND_TEMP1(LKE,i,j) = Vout(LKE-1,i,j)
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          LAND_TEMP1(LKE,i,j) = LAND_TEMP(LKE,i,j)
       enddo
       enddo
    endif


    ! calculate tendency
    if( LAND_PHY_SLAB_const ) then
       do j = JS, JE
       do i = IS, IE
       do k = LKS, LKE
          LAND_TEMP_t (k,i,j) = 0.0_RP
          LAND_WATER_t(k,i,j) = 0.0_RP
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
       do k = LKS, LKE
          LAND_TEMP_t (k,i,j) = ( LAND_TEMP1 (k,i,j) - LAND_TEMP (k,i,j) ) / dt
          LAND_WATER_t(k,i,j) = ( LAND_WATER1(k,i,j) - LAND_WATER(k,i,j) ) / dt
       enddo
       enddo
       enddo
    endif

    return
  end subroutine LAND_PHY_SLAB

end module scale_land_phy_slab
