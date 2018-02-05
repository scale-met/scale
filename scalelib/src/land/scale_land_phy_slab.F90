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

  use scale_external_input, only: &
     EXTIN_file_limit
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
  logical,               private :: LAND_PHY_update_bottom_temp  = .false. ! Is LAND_TEMP  updated in the lowest level?
  logical,               private :: LAND_PHY_update_bottom_water = .false. ! Is LAND_WATER updated in the lowest level?

  logical,               private :: LAND_PHY_SLAB_nudging                            = .false. ! Is nudging for land physics used?
  real(RP),              private :: LAND_PHY_SLAB_nudging_tau                        = 0.0_RP  ! time constant for nudging [sec]
  character(len=H_LONG), private :: LAND_PHY_SLAB_nudging_basename(EXTIN_file_limit) = ''
  logical,               private :: LAND_PHY_SLAB_nudging_enable_periodic_year       = .false.
  logical,               private :: LAND_PHY_SLAB_nudging_enable_periodic_month      = .false.
  logical,               private :: LAND_PHY_SLAB_nudging_enable_periodic_day        = .false.
  integer,               private :: LAND_PHY_SLAB_nudging_step_fixed                 = 0
  real(RP),              private :: LAND_PHY_SLAB_nudging_offset                     = 0.0_RP
  real(RP),              private :: LAND_PHY_SLAB_nudging_defval                   ! = UNDEF
  logical,               private :: LAND_PHY_SLAB_nudging_check_coordinates          = .true.
  integer,               private :: LAND_PHY_SLAB_nudging_step_limit                 = 0

  real(RP),              private :: WATER_DENSCS !< Heat Capacity (rho*CS) for soil moisture [J/K/m3]

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
    use scale_external_input, only: &
       EXTIN_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    NAMELIST / PARAM_LAND_PHY_SLAB / &
       LAND_PHY_SLAB_nudging,                       &
       LAND_PHY_SLAB_nudging_tau,                   &
       LAND_PHY_SLAB_nudging_basename,              &
       LAND_PHY_SLAB_nudging_enable_periodic_year,  &
       LAND_PHY_SLAB_nudging_enable_periodic_month, &
       LAND_PHY_SLAB_nudging_enable_periodic_day,   &
       LAND_PHY_SLAB_nudging_step_fixed,            &
       LAND_PHY_SLAB_nudging_offset,                &
       LAND_PHY_SLAB_nudging_defval,                &
       LAND_PHY_SLAB_nudging_check_coordinates,     &
       LAND_PHY_SLAB_nudging_step_limit,            &
       LAND_PHY_update_bottom_temp,                 &
       LAND_PHY_update_bottom_water

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[LAND PHY] / Origin[SCALElib]'

    LAND_PHY_SLAB_nudging_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PHY_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_PHY_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_PHY_SLAB)

    if ( LAND_PHY_SLAB_nudging ) then
       if ( LAND_PHY_SLAB_nudging_basename(1) == '' ) then
          write(*,*) 'xxx LAND_PHY_SLAB_nudging_basename is necessary !!'
          call PRC_MPIstop
       end if

       call EXTIN_regist( LAND_PHY_SLAB_nudging_basename(:),           & ! [IN]
                          'LAND_TEMP',                                 & ! [IN]
                          'LXY',                                       & ! [IN]
                          LAND_PHY_SLAB_nudging_enable_periodic_year,  & ! [IN]
                          LAND_PHY_SLAB_nudging_enable_periodic_month, & ! [IN]
                          LAND_PHY_SLAB_nudging_enable_periodic_day,   & ! [IN]
                          LAND_PHY_SLAB_nudging_step_fixed,            & ! [IN]
                          LAND_PHY_SLAB_nudging_offset,                & ! [IN]
                          LAND_PHY_SLAB_nudging_defval,                & ! [IN]
                          LAND_PHY_SLAB_nudging_check_coordinates,     & ! [IN]
                          LAND_PHY_SLAB_nudging_step_limit             ) ! [IN]

       call EXTIN_regist( LAND_PHY_SLAB_nudging_basename(:),           & ! [IN]
                          'LAND_WATER',                                & ! [IN]
                          'LXY',                                       & ! [IN]
                          LAND_PHY_SLAB_nudging_enable_periodic_year,  & ! [IN]
                          LAND_PHY_SLAB_nudging_enable_periodic_month, & ! [IN]
                          LAND_PHY_SLAB_nudging_enable_periodic_day,   & ! [IN]
                          LAND_PHY_SLAB_nudging_step_fixed,            & ! [IN]
                          LAND_PHY_SLAB_nudging_offset,                & ! [IN]
                          LAND_PHY_SLAB_nudging_defval,                & ! [IN]
                          LAND_PHY_SLAB_nudging_check_coordinates,     & ! [IN]
                          LAND_PHY_SLAB_nudging_step_limit             ) ! [IN]

       if( IO_L ) write(IO_FID_LOG,*) '*** Use nudging for Land physics: ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Use nudging for Land physics: OFF'
    end if

    WATER_DENSCS = DWATR * CL

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Update soil temperature of bottom layer? : ', LAND_PHY_update_bottom_temp
    if( IO_L ) write(IO_FID_LOG,*) '*** Update soil moisture    of bottom layer? : ', LAND_PHY_update_bottom_water

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
    use scale_const, only: &
       DWATR => CONST_DWATR
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_process, only: &
       PRC_MPIstop
    use scale_external_input, only: &
       EXTIN_update
    use scale_landuse, only: &
       LANDUSE_fact_land
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
    logical :: solve_matrix
    logical :: error

    real(RP) :: TEMP1 (LKMAX,IA,JA)
    real(RP) :: WATER1(LKMAX,IA,JA)

    real(RP) :: LAND_DENSCS(LKMAX,IA,JA)
    real(RP) :: ThermalDiff(LKMAX,IA,JA)

!    real(RP) :: RUNOFF(IA,JA)

    real(RP) :: U(LKMAX,IA,JA)
    real(RP) :: M(LKMAX,IA,JA)
    real(RP) :: L(LKMAX,IA,JA)
    real(RP) :: V(LKMAX,IA,JA)

    real(RP) :: NDG_TEMP (LKMAX,IA,JA)
    real(RP) :: NDG_WATER(LKMAX,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land physics step: Slab'

    if( LAND_PHY_SLAB_nudging ) then

      call EXTIN_update( &
           TEMP1,       & ! (out)
           'LAND_TEMP', & ! (in)
           NOWDAYSEC,   & ! (in)
           error        ) ! (out)
      if ( error ) then
         write(*,*) 'xxx Requested data is not found!'
         call PRC_MPIstop
      end if

      call EXTIN_update( &
           WATER1,       & ! (out)
           'LAND_WATER', & ! (in)
           NOWDAYSEC,    & ! (in)
           error         ) ! (out)
      if ( error ) then
         write(*,*) 'xxx Requested data is not found!'
         call PRC_MPIstop
      end if

      if( LAND_PHY_SLAB_nudging_tau > 0.0_RP ) then
        ! nudging is used
        solve_matrix = .true.

        NDG_TEMP (:,:,:) = ( TEMP1 (:,:,:) - TEMP (:,:,:) ) / LAND_PHY_SLAB_nudging_tau * dt
        NDG_WATER(:,:,:) = ( WATER1(:,:,:) - WATER(:,:,:) ) / LAND_PHY_SLAB_nudging_tau * dt

      else
        ! replace data to reference
        solve_matrix = .false.

      end if

    else
      ! nudging is NOT used
      solve_matrix = .true.

      NDG_TEMP (:,:,:) = 0.0_RP
      NDG_WATER(:,:,:) = 0.0_RP

    end if

    if( solve_matrix ) then

      ! Solve diffusion of soil moisture (tridiagonal matrix)
      do j = JS, JE
      do i = IS, IE
        L(LKS,i,j) = 0.0_RP
        U(LKS,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
        L(LKE,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
        U(LKE,i,j) = 0.0_RP

        M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
        M(LKE,i,j) = 1.0_RP - L(LKE,i,j) - U(LKE,i,j)
      end do
      end do

      do j = JS, JE
      do i = IS, IE
      do k = LKS+1, LKE-1
        L(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
        U(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
        M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
      end do
      end do
      end do

      ! input from atmosphere
      do j = JS, JE
      do i = IS, IE
        V(LKS,i,j) = WATER(LKS,i,j) + NDG_WATER(LKS,i,j) &
                   + ( SFLX_prec(i,j) - SFLX_evap(i,j) ) / ( CDZ(LKS) * DWATR ) * dt
      end do
      end do

      do j = JS, JE
      do i = IS, IE
      do k = LKS+1, LKE
        V(k,i,j) = WATER(k,i,j) + NDG_WATER(k,i,j)
      end do
      end do
      end do

      call MATRIX_SOLVER_tridiagonal( LKMAX,         & ! [IN]
                                      IA, IS, IE,    & ! [IN]
                                      JA, JS, JE,    & ! [IN]
                                      U     (:,:,:), & ! [IN]
                                      M     (:,:,:), & ! [IN]
                                      L     (:,:,:), & ! [IN]
                                      V     (:,:,:), & ! [IN]
                                      WATER1(:,:,:)  ) ! [OUT]

      if ( .not. LAND_PHY_UPDATE_BOTTOM_WATER ) then
        do j = JS, JE
        do i = IS, IE
          WATER1(LKE,i,j) = WATER(LKE,i,j)
        end do
        end do
      endif

      ! runoff of soil moisture (vertical sum)
      do j = JS, JE
      do i = IS, IE
!        RUNOFF(i,j) = 0.0_RP
        do k = LKS, LKE
!          RUNOFF(i,j) = RUNOFF(i,j) + max( WATER1(k,i,j) - WaterLimit(i,j), 0.0_RP ) * CDZ(k) * DWATR
          WATER1(k,i,j) = min( WATER1(k,i,j), WaterLimit(i,j) )
        end do
      end do
      end do

      ! estimate thermal diffusivity
      do j = JS, JE
      do i = IS, IE
      do k = LKS, LKE
        LAND_DENSCS(k,i,j) = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j) + WATER_DENSCS * WATER1(k,i,j)
        ThermalDiff(k,i,j) = ThermalCond(i,j) / LAND_DENSCS(k,i,j)
      end do
      end do
      end do

      ! Solve diffusion of soil temperature (tridiagonal matrix)
      do j = JS, JE
      do i = IS, IE
        L(LKS,i,j) = 0.0_RP
        U(LKS,i,j) = -2.0_RP * ThermalDiff(LKS,i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
        L(LKE,i,j) = -2.0_RP * ThermalDiff(LKE,i,j) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
        U(LKE,i,j) = 0.0_RP

        M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
        M(LKE,i,j) = 1.0_RP - L(LKE,i,j) - U(LKE,i,j)
      end do
      end do

      do j = JS, JE
      do i = IS, IE
      do k = LKS+1, LKE-1
        L(k,i,j) = -2.0_RP * ThermalDiff(k,i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
        U(k,i,j) = -2.0_RP * ThermalDiff(k,i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
        M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
      end do
      end do
      end do

      ! input from atmosphere
      do j = JS, JE
      do i = IS, IE
        V(LKS,i,j) = TEMP(LKS,i,j) + NDG_TEMP(LKS,i,j) &
                   - SFLX_GH(i,j) / ( LAND_DENSCS(LKS,i,j) * CDZ(LKS) ) * dt
      end do
      end do

      do j = JS, JE
      do i = IS, IE
      do k = LKS+1, LKE
        V(k,i,j) = TEMP(k,i,j) + NDG_TEMP(k,i,j)
      end do
      end do
      end do

      call MATRIX_SOLVER_tridiagonal( LKMAX,        & ! [IN]
                                      IA, IS, IE,   & ! [IN]
                                      JA, JS, JE,   & ! [IN]
                                      U    (:,:,:), & ! [IN]
                                      M    (:,:,:), & ! [IN]
                                      L    (:,:,:), & ! [IN]
                                      V    (:,:,:), & ! [IN]
                                      TEMP1(:,:,:)  ) ! [OUT]

      if ( .not. LAND_PHY_UPDATE_BOTTOM_TEMP ) then
        do j = JS, JE
        do i = IS, IE
          TEMP1(LKE,i,j) = TEMP(LKE,i,j)
        end do
        end do
      endif

    end if

    ! calculate tendency
    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then
        TEMP_t (k,i,j) = ( TEMP1 (k,i,j) - TEMP (k,i,j) ) / dt
        WATER_t(k,i,j) = ( WATER1(k,i,j) - WATER(k,i,j) ) / dt
      else
        TEMP_t (k,i,j) = 0.0_RP
        WATER_t(k,i,j) = 0.0_RP
      end if
    end do
    end do
    end do

    return
  end subroutine LAND_PHY_SLAB

end module scale_land_phy_slab
