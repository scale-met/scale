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
  real(RP),               private :: LAND_DYN_BUCKET_T_frz

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
  real(RP),               private :: ICE_DENSCS   !< Heat Capacity (rho*CS) for ice           [J/K/m3]
  real(DP),               private :: LAND_DYN_BUCKET_nudging_tausec  !< Relaxation time [sec]

  logical,                private :: replace = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_DYN_BUCKET_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       TEM00 => CONST_TEM00, &
       DWATR => CONST_DWATR, &
       DICE  => CONST_DICE,  &
       CL    => CONST_CL
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_regist
    use scale_atmos_hydrometeor, only: &
       CV_WATER, &
       CV_ICE
    implicit none

    namelist / PARAM_LAND_DYN_BUCKET / &
       LAND_DYN_BUCKET_T_frz,                         &
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
    LAND_DYN_BUCKET_T_frz          = TEM00

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

       if ( LAND_DYN_BUCKET_nudging_tausec <= 0.0_RP ) then
          LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Tau<=0 means that LST is completely replaced by the external data.'
          replace = .true.
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
                                        step_limit        = LAND_DYN_BUCKET_nudging_step_limit,        & ! [IN]
                                        allow_missing     = (.not. replace)                            ) ! [IN]

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
                                        step_limit        = LAND_DYN_BUCKET_nudging_step_limit,        & ! [IN]
                                        allow_missing     = (.not. replace)                            ) ! [IN]

       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Use nudging for Land physics: ON'
    else
       LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Use nudging for Land physics: OFF'
    end if

    WATER_DENSCS = DWATR * CV_WATER
    ICE_DENSCS   = DICE  * CV_ICE

    LOG_NEWLINE
    LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Update soil temperature of bottom layer? : ', LAND_DYN_BUCKET_update_bottom_temp
    LOG_INFO("LAND_DYN_BUCKET_setup",*) 'Update soil moisture    of bottom layer? : ', LAND_DYN_BUCKET_update_bottom_water

    return
  end subroutine LAND_DYN_BUCKET_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_DYN_BUCKET( &
       LKMAX, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
       TEMP_t, WATER_t, ICE_t, &
       WaterLimit,             &
       ThermalCond,            &
       HeatCapacity,           &
       WaterDiff,              &
       SFLX_GH, SFLX_water,    &
       SFLX_RHOE,              &
       exists_land, CDZ,       &
       dt, NOWDAYSEC,          &
       TEMP, WATER, ICE,       &
       RUNOFF, RUNOFF_ENGI     )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       DWATR => CONST_DWATR, &
       DICE  => CONST_DICE,  &
       EMELT => CONST_EMELT
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_atmos_hydrometeor, only: &
       CV_WATER, &
       CV_ICE,   &
       LHF
    implicit none
    integer, intent(in) :: LKMAX, LKS, LKE
    integer, intent(in) :: LIA, LIS, LIE
    integer, intent(in) :: LJA, LJS, LJE

    real(RP), intent(in) :: TEMP_t      (LKMAX,LIA,LJA)
    real(RP), intent(in) :: WATER_t     (LKMAX,LIA,LJA)
    real(RP), intent(in) :: ICE_t       (LKMAX,LIA,LJA)
    real(RP), intent(in) :: WaterLimit  (LIA,LJA)
    real(RP), intent(in) :: ThermalCond (LIA,LJA)
    real(RP), intent(in) :: HeatCapacity(LIA,LJA)
    real(RP), intent(in) :: WaterDiff   (LIA,LJA)
    real(RP), intent(in) :: SFLX_GH     (LIA,LJA) ! positive for downward
    real(RP), intent(in) :: SFLX_water  (LIA,LJA) ! positive for downward
    real(RP), intent(in) :: SFLX_RHOE   (LIA,LJA) ! positive for downward
    logical,  intent(in) :: exists_land (LIA,LJA)
    real(RP), intent(in) :: CDZ         (LKMAX)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: NOWDAYSEC

    real(RP), intent(inout) :: TEMP (LKMAX,LIA,LJA)
    real(RP), intent(inout) :: WATER(LKMAX,LIA,LJA)
    real(RP), intent(inout) :: ICE  (LKMAX,LIA,LJA)

    real(RP), intent(out) :: RUNOFF     (LIA,LJA)
    real(RP), intent(out) :: RUNOFF_ENGI(LIA,LJA) ! internal energy of the runoff water

    logical :: error

    real(RP) :: TEMP1 (LKMAX,LIA,LJA)
    real(RP) :: WATER1(LKMAX,LIA,LJA)
    real(RP) :: ICE1  (LKMAX,LIA,LJA)

    real(RP) :: kappa      (LKMAX)

    real(RP) :: U(LKMAX,LIA,LJA)
    real(RP) :: M(LKMAX,LIA,LJA)
    real(RP) :: L(LKMAX,LIA,LJA)
    real(RP) :: V(LKMAX,LIA,LJA)

    real(RP) :: NDG_TEMP (LKMAX,LIA,LJA)
    real(RP) :: NDG_WATER(LKMAX,LIA,LJA)

    real(RP) :: MASS_total(LKMAX)
    real(RP) :: MASS_water(LKMAX)
    real(RP) :: MASS_ice(LKMAX)

    real(RP) :: ENGI(LKMAX,LIA,LJA)
    real(RP) :: CS
    real(RP) :: CL

    real(RP) :: flux(LKS-1:LKE)

    real(RP) :: ro, rw, ri
    real(RP) :: sw

    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'land / dynamics / bucket'

    if ( LAND_DYN_BUCKET_nudging ) then

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

       if ( .not. replace ) then
          ! nudging is used

          !$omp parallel do
          do j = LJS,LJE
          do i = LIS,LIE
          do k = LKS,LKE
             if ( TEMP1(k,i,j) == UNDEF ) then
                NDG_TEMP (k,i,j) = 0.0_RP
             else
                NDG_TEMP (k,i,j) = ( TEMP1 (k,i,j) - TEMP (k,i,j) ) / LAND_DYN_BUCKET_nudging_tausec * dt
             end if
          end do
          end do
          end do

          !$omp parallel do
          do j = LJS,LJE
          do i = LIS,LIE
          do k = LKS,LKE
             if ( WATER1(k,i,j) == UNDEF ) then
                NDG_WATER(k,i,j) = 0.0_RP
             else
                NDG_WATER(k,i,j) = ( WATER1(k,i,j) - ( WATER(k,i,j) + ICE(k,i,j) ) ) / LAND_DYN_BUCKET_nudging_tausec * dt
             end if
          end do
          end do
          end do

       end if

       if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             NDG_WATER(LKE,i,j) = 0.0_RP
          end do
          end do
       end if

       if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             NDG_TEMP(LKE,i,j) = 0.0_RP
          end do
          end do
       end if

    else
       ! nudging is NOT used

       !$omp parallel do
       do j = LJS,LJE
       do i = LIS,LIE
       do k = LKS,LKE
          NDG_TEMP (k,i,j) = 0.0_RP
          NDG_WATER(k,i,j) = 0.0_RP
       end do
       end do
       end do

    end if

    if ( .not. replace ) then

       !$omp parallel do &
       !$omp private(MASS_total, MASS_water, MASS_ice, CS)
       do j = LJS, LJE
       do i = LIS, LIE
          if ( exists_land(i,j) ) then
             do k = LKS, LKE
                MASS_total(k) = DWATR * WATER(k,i,j) + DICE * ICE(k,i,j)
             end do
             MASS_total(LKS) = MASS_total(LKS) + dt * SFLX_water(i,j) / CDZ(LKS)

             CS = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j)
             do k = LKS, LKE
                ENGI(k,i,j) = ( CS + WATER_DENSCS * WATER(k,i,j) + ICE_DENSCS * ICE(k,i,j) ) * TEMP(k,i,j) - LHF * DICE * ICE(k,i,j)
             end do
             ENGI(LKS,i,j) = ENGI(LKS,i,j) + dt * ( SFLX_GH(i,j) + SFLX_RHOE(i,j) ) / CDZ(LKS)

             ! phase change
             do k = LKS, LKE
                MASS_ice(k) = min( MASS_total(k), max( 0.0_RP, &
                     ( ENGI(k,i,j) - ( CS + CV_WATER * MASS_total(k) ) * LAND_DYN_BUCKET_T_frz ) &
                     / ( ( CV_ICE - CV_WATER ) * LAND_DYN_BUCKET_T_frz - LHF ) &
                     ) )
                MASS_water(k) = MASS_total(k) - MASS_ice(k)
                V(k,i,j) = MASS_water(k) / DWATR
                ICE1(k,i,j) = MASS_ice(k) / DICE
                TEMP1(k,i,j) = ( ENGI(k,i,j) + LHF * MASS_ice(k) ) &
                     / ( CS + CV_WATER * MASS_water(k) + CV_ICE * MASS_ice(k) )
             end do
          end if
       end do
       end do


       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          if ( exists_land(i,j) ) then
             L(LKS,i,j) = 0.0_RP
             U(LKS,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
             M(LKS,i,j) = 1.0_RP - L(LKS,i,j) - U(LKS,i,j)
          end if
       end do
       end do

       if ( LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             if ( exists_land(i,j) ) then
                L(LKE,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
             end if
          end do
          end do
       else
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             if ( exists_land(i,j) ) then
                L(LKE,i,j) = 0.0_RP
             end if
          end do
          end do
       end if
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          if ( exists_land(i,j) ) then
             U(LKE,i,j) = 0.0_RP
             M(LKE,i,j) = 1.0_RP - L(LKE,i,j) - U(LKE,i,j)
          end if
       end do
       end do

       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
       do k = LKS+1, LKE-1
          if ( exists_land(i,j) ) then
             L(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
             U(k,i,j) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
             M(k,i,j) = 1.0_RP - L(k,i,j) - U(k,i,j)
          end if
       end do
       end do
       end do

       call MATRIX_SOLVER_tridiagonal( LKMAX, 1, LKMAX, &
                                       LIA, LIS, LIE,   &
                                       LJA, LJS, LJE,   &
                                       U(:,:,:), M(:,:,:), L(:,:,:), & ! [IN]
                                       V(:,:,:),                     & ! [IN]
                                       WATER1(:,:,:),                & ! [OUT]
                                       mask = exists_land(:,:)       ) ! [IN]

       ! temperature

       flux(LKS-1) = 0.0_RP
       flux(LKE)   = 0.0_RP

       !$omp parallel do &
       !$omp private(kappa,CS,CL,sw) &
       !$omp firstprivate(flux)
       do j = LJS, LJE
       do i = LIS, LIE
          if ( exists_land(i,j) ) then

             CS = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j)
             do k = LKS, LKE
                kappa(k) = ThermalCond(i,j) + 0.5_RP * WATER1(k,i,j)**(1.0_RP/3.0_RP)
             end do

             do k = LKS, LKE-1
                flux(k) = - 2.0_RP *  DWATR * WaterDiff(i,j) * ( WATER1(k+1,i,j) - WATER1(k,i,j) ) / ( CDZ(k+1) + CDZ(k) )
                sw = 0.5_RP - sign( 0.5_RP, flux(k) )
                flux(k) = flux(k) * CV_WATER * ( TEMP1(k+1,i,j) * sw + TEMP1(k,i,j) * ( 1.0_RP - sw ) )
             end do
             if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
                flux(LKE) = flux(LKE-1)
             end if

             do k = LKS, LKE
                V(k,i,j) = ENGI(k,i,j) + LHF * DICE * ICE1(k,i,j) &
                         - dt * ( flux(k) - flux(k-1) ) / CDZ(k)
             end do

             CL = CS + WATER_DENSCS * WATER1(LKS,i,j) + ICE_DENSCS * ICE1(LKS,i,j)
             L(LKS,i,j) = 0.0_RP
             U(LKS,i,j) = - ( kappa(LKS) + kappa(LKS+1) ) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
             M(LKS,i,j) = CL - L(LKS,i,j) - U(LKS,i,j)

             CL = CS + WATER_DENSCS * WATER1(LKE,i,j) + ICE_DENSCS * ICE1(LKE,i,j)
             if ( LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then

                L(LKE,i,j) = - ( kappa(LKE) + kappa(LKE-1) ) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
             else
                L(LKE,i,j) = 0.0_RP
             end if
             U(LKE,i,j) = 0.0_RP
             M(LKE,i,j) = CL - L(LKE,i,j) - U(LKE,i,j)

             do k = LKS+1, LKE-1
                CL = CS + WATER_DENSCS * WATER1(k,i,j) + ICE_DENSCS * ICE1(k,i,j)
                L(k,i,j) = - ( kappa(k) + kappa(k-1) ) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
                U(k,i,j) = - ( kappa(k) + kappa(k+1) ) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
                M(k,i,j) = CL - L(k,i,j) - U(k,i,j)
             end do

          end if
       end do
       end do

       call MATRIX_SOLVER_tridiagonal( LKMAX, 1, LKMAX, &
                                       LIA, LIS, LIE,   &
                                       LJA, LJS, LJE,   &
                                       U(:,:,:), M(:,:,:), L(:,:,:), & ! [IN]
                                       V(:,:,:),                     & ! [IN]
                                       TEMP1(:,:,:),                 & ! [OUT]
                                       mask = exists_land(:,:)       ) ! [IN]


       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          if ( exists_land(i,j) ) then
             do k = LKS, LKE
                TEMP1(k,i,j) = TEMP1(k,i,j) + NDG_TEMP(k,i,j)
             end do
          end if
       end do
       end do


       !$omp parallel do &
       !$omp private(ro,rw,ri)
       do j = LJS, LJE
       do i = LIS, LIE
          if ( exists_land(i,j) ) then
             do k = LKS, LKE
                if ( TEMP1(k,i,j) >= LAND_DYN_BUCKET_T_frz ) then
                   WATER1(k,i,j) = WATER1(k,i,j) + NDG_WATER(k,i,j)
                else
                   ICE1(k,i,j) = ICE1(k,i,j) + NDG_WATER(k,i,j)
                end if
             end do

             ! runoff of soil moisture (vertical sum)
             RUNOFF     (i,j) = 0.0_RP
             RUNOFF_ENGI(i,j) = 0.0_RP
             do k = LKS, LKE
                ro = max( WATER1(k,i,j) + ICE1(k,i,j) - WaterLimit(i,j), 0.0_RP )
                rw = min( ro, WATER1(k,i,j) )
                ri = ro - rw
                WATER1(k,i,j) = WATER1(k,i,j) - rw
                ICE1(k,i,j) = ICE1(k,i,j) - ri
                rw = rw * DWATR / dt
                ri = ri * DICE / dt
                RUNOFF(i,j) = RUNOFF(i,j) + ( rw + ri ) * CDZ(k)
                RUNOFF_ENGI(i,j) = RUNOFF_ENGI(i,j) &
                     + ( ( rw * CV_WATER + ri * CV_ICE ) * TEMP1(k,i,j) - ri * LHF ) * CDZ(k)
             end do
          end if
       end do
       end do


       if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             if ( exists_land(i,j) ) then
                WATER1(LKE,i,j) = WATER(LKE,i,j)
                ICE1  (LKE,i,j) = ICE  (LKE,i,j)
             end if
          end do
          end do
       endif

       if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             if ( exists_land(i,j) ) then
                TEMP1(LKE,i,j) = TEMP(LKE,i,j)
             end if
          end do
          end do
       endif

    end if


    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       if( exists_land(i,j) ) then
          do k = LKS, LKE
             TEMP (k,i,j) = TEMP1 (k,i,j)
             WATER(k,i,j) = WATER1(k,i,j)
             ICE  (k,i,j) = ICE1  (k,i,j)
          end do
       end if
    end do
    end do

    return
  end subroutine LAND_DYN_BUCKET

end module scale_land_dyn_bucket
