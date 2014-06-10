!-------------------------------------------------------------------------------
!> module LAND / Physics Bucket
!!
!! @par Description
!!          bucket-type land physics module
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_phy_bucket
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
  public :: LAND_PHY_driver_setup
  public :: LAND_PHY_driver

  public :: LAND_PHY_bucket_setup
  public :: LAND_PHY_bucket

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: SOLVER_tridiagonal_matrix_1D
  private :: SOLVER_tridiagonal_matrix_3D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: LAND_PHY_UPDATE_BOTTOM_TEMP  = .false. ! Is LAND_TEMP  updated in the lowest level?
  logical, private :: LAND_PHY_UPDATE_BOTTOM_WATER = .false. ! Is LAND_WATER updated in the lowest level?

  real(RP), private :: WATER_DENSCP              !< rho*CP for soil moisture [J/K/m3]

  real(RP), private, parameter :: DSOIL = 1.0_RP !< density of soil [kg/m3]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_PHY_driver_setup
    use mod_land_admin, only: &
       LAND_TYPE, &
       LAND_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[LAND PHY] / Origin[SCALE-LES]'

    if ( LAND_sw ) then

       call LAND_PHY_bucket_setup( LAND_TYPE )

       call LAND_PHY_driver( .true., .false. )

    endif

    return
  end subroutine LAND_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine LAND_PHY_driver( update_flag, history_flag )
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_land_grid, only: &
       GRID_LCDZ
    use scale_history, only: &
       HIST_in
    use mod_land_vars, only: &
       I_WaterLimit,   &
       I_ThermalCond,  &
       I_HeatCapacity, &
       I_WaterDiff,    &
       LAND_TEMP,      &
       LAND_WATER,     &
       LAND_TEMP_t,    &
       LAND_WATER_t,   &
       LAND_PROPERTY
    use mod_cpl_admin, only: &
       CPL_sw_AtmLnd
    use mod_cpl_vars, only: &
       CPL_getLnd
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    real(RP) :: FLX_heat  (IA,JA)
    real(RP) :: FLX_precip(IA,JA)
    real(RP) :: FLX_evap  (IA,JA)
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       if ( CPL_sw_AtmLnd ) then
          call CPL_getLnd( FLX_heat  (:,:), & ! [OUT]
                           FLX_precip(:,:), & ! [OUT]
                           FLX_evap  (:,:)  ) ! [OUT]
       endif

       call LAND_PHY_bucket( LAND_TEMP    (:,:,:),              & ! [IN]
                             LAND_WATER   (:,:,:),              & ! [IN]
                             LAND_PROPERTY(:,:,I_WaterLimit),   & ! [IN]
                             LAND_PROPERTY(:,:,I_ThermalCond),  & ! [IN]
                             LAND_PROPERTY(:,:,I_HeatCapacity), & ! [IN]
                             LAND_PROPERTY(:,:,I_WaterDiff),    & ! [IN]
                             FLX_heat     (:,:),                & ! [IN]
                             FLX_precip   (:,:),                & ! [IN]
                             FLX_evap     (:,:),                & ! [IN]
                             GRID_LCDZ    (:),                  & ! [IN]
                             LAND_TEMP_t  (:,:,:),              & ! [OUT]
                             LAND_WATER_t (:,:,:)               ) ! [OUT]

       if ( history_flag ) then
          call HIST_in( LAND_TEMP_t (:,:,:), 'LAND_TEMP_t',  'Soil temperature tendency', 'K',     dt )
          call HIST_in( LAND_WATER_t(:,:,:), 'LAND_WATER_t', 'Soil moisture    tendency', 'm3/m3', dt )
       endif

    endif

    return
  end subroutine LAND_PHY_driver

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_PHY_bucket_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    NAMELIST / PARAM_LAND_BUCKET / &
       LAND_PHY_UPDATE_BOTTOM_TEMP,   &
       LAND_PHY_UPDATE_BOTTOM_WATER

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[BUCKET] / Categ[LAND PHY] / Origin[SCALE-LES]'

    if ( LAND_TYPE /= 'BUCKET' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx LAND_TYPE is not BUCKET. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_BUCKET,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_BUCKET. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_BUCKET)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Update soil temperature of bottom layer? : ', LAND_PHY_UPDATE_BOTTOM_TEMP
    if( IO_L ) write(IO_FID_LOG,*) '*** Update soil moisture    of bottom layer? : ', LAND_PHY_UPDATE_BOTTOM_WATER

    return
  end subroutine LAND_PHY_bucket_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_bucket( &
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
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
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

    real(RP) :: LAND_TEMP1 (LKMAX,IA,JA)
    real(RP) :: LAND_WATER1(LKMAX,IA,JA)
    real(RP) :: LAND_RUNOFF(IA,JA)
    real(RP) :: SOIL_DENSCP(IA,JA)

    real(RP) :: U   (LKMAX-1,IA,JA)
    real(RP) :: M   (LKMAX-1,IA,JA)
    real(RP) :: L   (LKMAX-1,IA,JA)
    real(RP) :: Vin (LKMAX-1,IA,JA)
    real(RP) :: Vout(LKMAX-1,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land step: Bucket'

    ! Solve diffusion of soil moisture (tridiagonal matrix)
    do j = JS, JE
    do i = IS, IE
       L(LKS,i,j) = 0.0_RP
       U(LKS,i,j) = - 2.0_RP * LAND_WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
       M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = LKS+1, LKE-1
       L(k,i,j) = - 2.0_RP * LAND_WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
       U(k,i,j) = - 2.0_RP * LAND_WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
       M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Vin(LKS  ,i,j) = LAND_WATER(LKS,i,j) &
                      + ( FLX_precip(i,j) + FLX_evap(i,j) ) / ( CDZ(LKS) * DWATR ) * dt ! input from atmosphere

       do k = LKS+1, LKE-2
          Vin(k,i,j) = LAND_WATER(k,i,j)
       enddo

       Vin(LKE-1,i,j) = LAND_WATER(LKE-1,i,j) &
                      - U(LKE-1,i,j) * Vin(LKE,i,j)
    enddo
    enddo

    call SOLVER_tridiagonal_matrix_3D( LKMAX-1,     & ! [IN]
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
    endif



    ! runoff of soil moisture (vertical sum)
    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
       LAND_RUNOFF(i,j)   = LAND_RUNOFF(i,j) &
                          + max( LAND_WATER1(k,i,j)-LAND_WaterLimit(i,j), 0.0_RP ) * CDZ(k) * DWATR

       LAND_WATER1(k,i,j) = min( LAND_WATER1(k,i,j), LAND_WaterLimit(i,j) )
    enddo
    enddo
    enddo



    ! Solve diffusion of soil temperature (tridiagonal matrix)
    do j = JS, JE
    do i = IS, IE
       SOIL_DENSCP(i,j) = ( 1.0_RP - LAND_WaterLimit(i,j) ) * DSOIL * LAND_HeatCapacity(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       L(LKS,i,j) = 0.0_RP
       U(LKS,i,j) = - LAND_ThermalCond(i,j) / ( SOIL_DENSCP(i,j) + LAND_WATER(LKS,i,j)*WATER_DENSCP ) &
                    / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt

       M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = LKS+1, LKE-1
       L(k,i,j) = - 2.0_RP * LAND_ThermalCond(i,j) / ( SOIL_DENSCP(i,j) + LAND_WATER(k,i,j)*WATER_DENSCP ) &
                  / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
       U(k,i,j) = - 2.0_RP * LAND_ThermalCond(i,j) / ( SOIL_DENSCP(i,j) + LAND_WATER(k,i,j)*WATER_DENSCP ) &
                  / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
       M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       Vin(LKS  ,i,j) = LAND_TEMP(LKS,i,j) &
                      + FLX_heat(i,j) / ( SOIL_DENSCP(i,j) + LAND_WATER(LKS,i,j)*WATER_DENSCP ) / CDZ(LKS) * dt ! input from atmosphere

       do k = LKS+1, LKE-2
          Vin(k,i,j) = LAND_TEMP(k,i,j)
       enddo

       Vin(LKE-1,i,j) = LAND_TEMP(LKE-1,i,j) &
                      - U(LKE-1,i,j) * Vin(LKE,i,j)
    enddo
    enddo

    call SOLVER_tridiagonal_matrix_3D( LKMAX-1,     & ! [IN]
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
    endif



    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = ( LAND_TEMP1 (k,i,j) - LAND_TEMP (k,i,j) ) / dt
       LAND_WATER_t(k,i,j) = ( LAND_WATER1(k,i,j) - LAND_WATER(k,i,j) ) / dt
    enddo
    enddo
    enddo

    return
  end subroutine LAND_PHY_bucket

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Thomas's algorithm
  subroutine SOLVER_tridiagonal_matrix_1D( &
       KA, &
       ud,   &
       md,   &
       ld,   &
       iv,   &
       ov    )
    implicit none

    integer,  intent(in)  :: KA     ! array size
    real(RP), intent(in)  :: ud(KA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA) ! input  vector
    real(RP), intent(out) :: ov(KA) ! output vector

    real(RP) :: c(KA)
    real(RP) :: d(KA)

    integer :: k
    !---------------------------------------------------------------------------

    ! foward reduction
    c(1) = ud(1) / md(1)
    d(1) = iv(1) / md(1)
    do k = 2, KA
       c(k) =           ud(k)            / ( md(k) - ld(k) * c(k-1) )
       d(k) = ( iv(k) - ld(k) * d(k-1) ) / ( md(k) - ld(k) * c(k-1) )
    enddo

    ! backward substitution
    ov(KA) = d(KA)
    do k = KA-1, 1, -1
       ov(k) = d(k) - c(k) * ov(k+1)
    enddo

    return
  end subroutine SOLVER_tridiagonal_matrix_1D

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Thomas's algorithm
  subroutine SOLVER_tridiagonal_matrix_3D( &
       KA,         &
       IA, IS, IE, &
       JA, JS, JE, &
       ud,         &
       md,         &
       ld,         &
       iv,         &
       ov          )
    implicit none

    integer,  intent(in)  :: KA           ! array size
    integer,  intent(in)  :: IA, IS, IE   ! array size
    integer,  intent(in)  :: JA, JS, JE   ! array size
    real(RP), intent(in)  :: ud(KA,IA,JA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,IA,JA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,IA,JA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,IA,JA) ! input  vector
    real(RP), intent(out) :: ov(KA,IA,JA) ! output vector

    real(RP) :: c(KA,IA,JA)
    real(RP) :: d(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! foward reduction
    do j = JS, JE
    do i = IS, IE
       c(1,i,j) = ud(1,i,j) / md(1,i,j)
       d(1,i,j) = iv(1,i,j) / md(1,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 2,  KA
       c(k,i,j) =               ud(k,i,j)                &
                / ( md(k,i,j) - ld(k,i,j) * c(k-1,i,j) )
       d(k,i,j) = ( iv(k,i,j) - ld(k,i,j) * d(k-1,i,j) ) &
                / ( md(k,i,j) - ld(k,i,j) * c(k-1,i,j) )
    enddo
    enddo
    enddo

    ! backward substitution
    do j = JS, JE
    do i = IS, IE
       ov(KA,i,j) = d(KA,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KA-1, 1, -1
       ov(k,i,j) = d(k,i,j) - c(k,i,j) * ov(k+1,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine SOLVER_tridiagonal_matrix_3D

end module mod_land_phy_bucket
