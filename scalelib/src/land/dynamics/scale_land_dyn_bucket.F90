!-------------------------------------------------------------------------------
!> module land / dynamics / bucket
!!
!! @par Description
!!          slab-type land model
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_dyn_bucket
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_DYN_BUCKET_setup
  public :: LAND_DYN_BUCKET

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
  logical,                private :: LAND_DYN_BUCKET_update_bottom_temp  = .false. ! Is LAND_TEMP  updated in the lowest level?
  logical,                private :: LAND_DYN_BUCKET_update_bottom_water = .false. ! Is LAND_WATER updated in the lowest level?

  logical,                private :: LAND_DYN_BUCKET_nudging                                          = .false. ! Is nudging for land physics used?
  real(DP),               private :: LAND_DYN_BUCKET_nudging_tau                                      = 0.0_DP  ! time constant for nudging [sec]
  character(len=H_SHORT), private :: LAND_DYN_BUCKET_nudging_tau_unit                                 = "SEC"
  character(len=H_LONG),  private :: LAND_DYN_BUCKET_nudging_basename              = ''
  logical,                private :: LAND_DYN_BUCKET_nudging_basename_add_num      = .false.
  integer,                private :: LAND_DYN_BUCKET_nudging_number_of_files       = 1
  logical,                private :: LAND_DYN_BUCKET_nudging_enable_periodic_year  = .false.
  logical,                private :: LAND_DYN_BUCKET_nudging_enable_periodic_month = .false.
  logical,                private :: LAND_DYN_BUCKET_nudging_enable_periodic_day   = .false.
  integer,                private :: LAND_DYN_BUCKET_nudging_step_fixed            = 0
  real(RP),               private :: LAND_DYN_BUCKET_nudging_offset                = 0.0_RP
  real(RP),               private :: LAND_DYN_BUCKET_nudging_defval                ! = UNDEF
  logical,                private :: LAND_DYN_BUCKET_nudging_check_coordinates     = .true.
  integer,                private :: LAND_DYN_BUCKET_nudging_step_limit            = 0

  real(RP),               private :: WATER_DENSCS !< Heat Capacity (rho*CS) for soil moisture [J/K/m3]
  real(DP),               private :: LAND_DYN_BUCKET_nudging_tausec  !< Relaxation time [sec]
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_DYN_BUCKET_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    namelist / PARAM_LAND_DYN_BUCKET / &
       LAND_DYN_BUCKET_nudging,                       &
       LAND_DYN_BUCKET_nudging_tau,                   &
       LAND_DYN_BUCKET_nudging_tau_unit,              &
       LAND_DYN_BUCKET_nudging_basename,              &
       LAND_DYN_BUCKET_nudging_basename_add_num,      &
       LAND_DYN_BUCKET_nudging_number_of_files,       &
       LAND_DYN_BUCKET_nudging_enable_periodic_year,  &
       LAND_DYN_BUCKET_nudging_enable_periodic_month, &
       LAND_DYN_BUCKET_nudging_enable_periodic_day,   &
       LAND_DYN_BUCKET_nudging_step_fixed,            &
       LAND_DYN_BUCKET_nudging_offset,                &
       LAND_DYN_BUCKET_nudging_defval,                &
       LAND_DYN_BUCKET_nudging_check_coordinates,     &
       LAND_DYN_BUCKET_nudging_step_limit,            &
       LAND_DYN_BUCKET_update_bottom_temp,            &
       LAND_DYN_BUCKET_update_bottom_water

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Setup'

    LAND_DYN_BUCKET_nudging_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_DYN_BUCKET,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_DYN_BUCKET_setup",*) 'Not appropriate names in namelist PARAM_LAND_DYN_BUCKET. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND_DYN_BUCKET)

    if ( LAND_DYN_BUCKET_nudging ) then
       call CALENDAR_unit2sec( LAND_DYN_BUCKET_nudging_tausec, LAND_DYN_BUCKET_nudging_tau, LAND_DYN_BUCKET_nudging_tau_unit )

       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Use nudging for LAND physics  : ON'
       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Relaxation time Tau [sec]     : ', LAND_DYN_BUCKET_nudging_tausec

       if ( LAND_DYN_BUCKET_nudging_tausec == 0.0_RP ) then
          LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Tau=0 means that LST is completely replaced by the external data.'
       endif

       if ( LAND_DYN_BUCKET_nudging_basename == '' ) then
          LOG_ERROR("LAND_DYN_BUCKET_setup",*) 'LAND_DYN_BUCKET_nudging_basename is necessary !!'
          call PRC_abort
       end if

       call FILE_EXTERNAL_INPUT_regist( LAND_DYN_BUCKET_nudging_basename,              & ! [IN]
                                        LAND_DYN_BUCKET_nudging_basename_add_num,      & ! [IN]
                                        LAND_DYN_BUCKET_nudging_number_of_files,       & ! [IN]
                                        'LAND_TEMP',                                   & ! [IN]
                                        'LXY',                                         & ! [IN]
                                        LAND_DYN_BUCKET_nudging_enable_periodic_year,  & ! [IN]
                                        LAND_DYN_BUCKET_nudging_enable_periodic_month, & ! [IN]
                                        LAND_DYN_BUCKET_nudging_enable_periodic_day,   & ! [IN]
                                        LAND_DYN_BUCKET_nudging_step_fixed,            & ! [IN]
                                        LAND_DYN_BUCKET_nudging_offset,                & ! [IN]
                                        LAND_DYN_BUCKET_nudging_defval,                & ! [IN]
                                        check_coordinates = LAND_DYN_BUCKET_nudging_check_coordinates, & ! [IN]
                                        step_limit = LAND_DYN_BUCKET_nudging_step_limit                ) ! [IN]

       call FILE_EXTERNAL_INPUT_regist( LAND_DYN_BUCKET_nudging_basename,              & ! [IN]
                                        LAND_DYN_BUCKET_nudging_basename_add_num,      & ! [IN]
                                        LAND_DYN_BUCKET_nudging_number_of_files,       & ! [IN]
                                        'LAND_WATER',                                  & ! [IN]
                                        'LXY',                                         & ! [IN]
                                        LAND_DYN_BUCKET_nudging_enable_periodic_year,  & ! [IN]
                                        LAND_DYN_BUCKET_nudging_enable_periodic_month, & ! [IN]
                                        LAND_DYN_BUCKET_nudging_enable_periodic_day,   & ! [IN]
                                        LAND_DYN_BUCKET_nudging_step_fixed,            & ! [IN]
                                        LAND_DYN_BUCKET_nudging_offset,                & ! [IN]
                                        LAND_DYN_BUCKET_nudging_defval,                & ! [IN]
                                        check_coordinates = LAND_DYN_BUCKET_nudging_check_coordinates, & ! [IN]
                                        step_limit = LAND_DYN_BUCKET_nudging_step_limit                ) ! [IN]

       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Use nudging for Land physics: ON'
    else
       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Use nudging for Land physics: OFF'
    end if

    WATER_DENSCS = DWATR * CL

    LOG_NEWLINE
    LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Update soil temperature of bottom layer? : ', LAND_DYN_BUCKET_update_bottom_temp
    LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Update soil moisture    of bottom layer? : ', LAND_DYN_BUCKET_update_bottom_water

    return
  end subroutine LAND_DYN_BUCKET_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_DYN_BUCKET( &
       LKMAX, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
       TEMP_t, WATER_t,               &
       WaterLimit,                    &
       ThermalCond,                   &
       HeatCapacity,                  &
       WaterDiff,                     &
       SFLX_GH, SFLX_water, SFLX_ice, &
       fact_land, CDZ,                &
       dt, NOWDAYSEC,                 &
       TEMP, WATER,                   &
       RUNOFF                         )
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       EMELT => CONST_EMELT
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none
    integer, intent(in) :: LKMAX, LKS, LKE
    integer, intent(in) :: LIA, LIS, LIE
    integer, intent(in) :: LJA, LJS, LJE

    real(RP), intent(in) :: TEMP_t      (LKMAX,LIA,LJA)
    real(RP), intent(in) :: WATER_t     (LKMAX,LIA,LJA)
    real(RP), intent(in) :: WaterLimit  (LIA,LJA)
    real(RP), intent(in) :: ThermalCond (LIA,LJA)
    real(RP), intent(in) :: HeatCapacity(LIA,LJA)
    real(RP), intent(in) :: WaterDiff   (LIA,LJA)
    real(RP), intent(in) :: SFLX_GH     (LIA,LJA) ! positive for downward
    real(RP), intent(in) :: SFLX_water  (LIA,LJA) ! positive for downward
    real(RP), intent(in) :: SFLX_ice    (LIA,LJA) ! positive for downward
    real(RP), intent(in) :: fact_land   (LIA,LJA)
    real(RP), intent(in) :: CDZ         (LKMAX)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: NOWDAYSEC

    real(RP), intent(inout) :: TEMP (LKMAX,LIA,LJA)
    real(RP), intent(inout) :: WATER(LKMAX,LIA,LJA)

    real(RP), intent(out) :: RUNOFF(LIA,LJA)

    logical :: solve_matrix
    logical :: error

    real(RP) :: TEMP1 (LKMAX,LIA,LJA)
    real(RP) :: WATER1(LKMAX,LIA,LJA)

    real(RP) :: LAND_DENSCS(LKMAX,LIA,LJA)
    real(RP) :: ThermalDiff(LKMAX,LIA,LJA)

    real(RP) :: U(LKMAX,LIA,LJA)
    real(RP) :: M(LKMAX,LIA,LJA)
    real(RP) :: L(LKMAX,LIA,LJA)
    real(RP) :: V(LKMAX,LIA,LJA)

    real(RP) :: NDG_TEMP (LKMAX,LIA,LJA)
    real(RP) :: NDG_WATER(LKMAX,LIA,LJA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'land / dynamics / bucket'

    if( LAND_DYN_BUCKET_nudging ) then

      call FILE_EXTERNAL_INPUT_update( &
                          'LAND_TEMP', & ! (in)
                          NOWDAYSEC,   & ! (in)
                          TEMP1,       & ! (out)
                          error        ) ! (out)
      if ( error ) then
         LOG_ERROR("LAND_DYN_BUCKET",*) 'Requested data is not found!'
         call PRC_abort
      end if

      call FILE_EXTERNAL_INPUT_update( &
                         'LAND_WATER', & ! (in)
                         NOWDAYSEC,    & ! (in)
                         WATER1,       & ! (out)
                         error         ) ! (out)
      if ( error ) then
         LOG_ERROR("LAND_DYN_BUCKET",*) 'Requested data is not found!'
         call PRC_abort
      end if

      if( LAND_DYN_BUCKET_nudging_tau > 0.0_RP ) then
        ! nudging is used
        solve_matrix = .true.

        !$omp parallel do
        do j = LJS,LJE
        do i = LIS,LIE
        do k = LKS,LKE
           NDG_TEMP (k,i,j) = ( TEMP1 (k,i,j) - TEMP (k,i,j) ) / LAND_DYN_BUCKET_nudging_tausec * dt
           NDG_WATER(k,i,j) = ( WATER1(k,i,j) - WATER(k,i,j) ) / LAND_DYN_BUCKET_nudging_tausec * dt
        end do
        end do
        end do

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
      do j = LJS, LJE
      do i = LIS, LIE
        L(LKS,i,j) = 0.0_RP
        U(LKS,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
        L(LKE,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
        U(LKE,i,j) = 0.0_RP

        M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
        M(LKE,i,j) = 1.0_RP - L(LKE,i,j) - U(LKE,i,j)
      end do
      end do

      do j = LJS, LJE
      do i = LIS, LIE
      do k = LKS+1, LKE-1
        L(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
        U(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
        M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
      end do
      end do
      end do

      ! input from atmosphere
      do j = LJS, LJE
      do i = LIS, LIE
        V(LKS,i,j) = WATER(LKS,i,j) + NDG_WATER(LKS,i,j) &
                   + ( SFLX_water(i,j) + SFLX_ice(i,j) ) / ( CDZ(LKS) * DWATR ) * dt
      end do
      end do

      do j = LJS, LJE
      do i = LIS, LIE
      do k = LKS+1, LKE
        V(k,i,j) = WATER(k,i,j) + NDG_WATER(k,i,j)
      end do
      end do
      end do

      call MATRIX_SOLVER_tridiagonal( LKMAX, 1, LKMAX, & ! [IN]
                                      LIA, LIS, LIE, & ! [IN]
                                      LJA, LJS, LJE, & ! [IN]
                                      U     (:,:,:), & ! [IN]
                                      M     (:,:,:), & ! [IN]
                                      L     (:,:,:), & ! [IN]
                                      V     (:,:,:), & ! [IN]
                                      WATER1(:,:,:)  ) ! [OUT]

      if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
        do j = LJS, LJE
        do i = LIS, LIE
          WATER1(LKE,i,j) = WATER(LKE,i,j)
        end do
        end do
      endif

      ! runoff of soil moisture (vertical sum)
      do j = LJS, LJE
      do i = LIS, LIE
        RUNOFF(i,j) = 0.0_RP
        do k = LKS, LKE
          RUNOFF(i,j) = RUNOFF(i,j) + max( WATER1(k,i,j) - WaterLimit(i,j), 0.0_RP ) * CDZ(k) * DWATR
          WATER1(k,i,j) = min( WATER1(k,i,j), WaterLimit(i,j) )
        end do
      end do
      end do

      ! estimate thermal diffusivity
      do j = LJS, LJE
      do i = LIS, LIE
      do k = LKS, LKE
        LAND_DENSCS(k,i,j) = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j) + WATER_DENSCS * WATER1(k,i,j)
        ThermalDiff(k,i,j) = ThermalCond(i,j) / LAND_DENSCS(k,i,j)
      end do
      end do
      end do

      ! Solve diffusion of soil temperature (tridiagonal matrix)
      do j = LJS, LJE
      do i = LIS, LIE
        L(LKS,i,j) = 0.0_RP
        U(LKS,i,j) = -2.0_RP * ThermalDiff(LKS,i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
        L(LKE,i,j) = -2.0_RP * ThermalDiff(LKE,i,j) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
        U(LKE,i,j) = 0.0_RP

        M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
        M(LKE,i,j) = 1.0_RP - L(LKE,i,j) - U(LKE,i,j)
      end do
      end do

      do j = LJS, LJE
      do i = LIS, LIE
      do k = LKS+1, LKE-1
        L(k,i,j) = -2.0_RP * ThermalDiff(k,i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
        U(k,i,j) = -2.0_RP * ThermalDiff(k,i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
        M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
      end do
      end do
      end do

      ! input from atmosphere
      do j = LJS, LJE
      do i = LIS, LIE
         V(LKS,i,j) = TEMP(LKS,i,j) + NDG_TEMP(LKS,i,j) + TEMP_t(LKS,i,j) &
                    + ( SFLX_GH(i,j) - SFLX_ice(i,j) * EMELT ) / ( LAND_DENSCS(LKS,i,j) * CDZ(LKS) ) * dt
      end do
      end do

      do j = LJS, LJE
      do i = LIS, LIE
      do k = LKS+1, LKE
        V(k,i,j) = TEMP(k,i,j) + NDG_TEMP(k,i,j) + TEMP_t(k,i,j)
      end do
      end do
      end do

      call MATRIX_SOLVER_tridiagonal( LKMAX, 1, LKMAX, & ! [IN]
                                      LIA, LIS, LIE, & ! [IN]
                                      LJA, LJS, LJE, & ! [IN]
                                      U    (:,:,:),  & ! [IN]
                                      M    (:,:,:),  & ! [IN]
                                      L    (:,:,:),  & ! [IN]
                                      V    (:,:,:),  & ! [IN]
                                      TEMP1(:,:,:)   ) ! [OUT]

      if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
        do j = LJS, LJE
        do i = LIS, LIE
          TEMP1(LKE,i,j) = TEMP(LKE,i,j)
        end do
        end do
      endif

    end if

    ! calculate tendency
    do j = LJS, LJE
    do i = LIS, LIE
      if( fact_land(i,j) > 0.0_RP ) then
         do k = LKS, LKE
            TEMP (k,i,j) = TEMP1 (k,i,j)
            WATER(k,i,j) = WATER1(k,i,j)
         end do
      end if
    end do
    end do

    return
  end subroutine LAND_DYN_BUCKET

end module scale_land_dyn_bucket
