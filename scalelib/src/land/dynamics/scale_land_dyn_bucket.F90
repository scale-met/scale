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
       MATRIX_SOLVER_tridiagonal, &
       MATRIX_SOLVER_TRIDIAGONAL_1D_CR
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

    real(RP) :: kappa      (LKMAX,LSIZE)

   !  real(RP) :: U(LKMAX,LIA,LJA)
   !  real(RP) :: M(LKMAX,LIA,LJA)
   !  real(RP) :: L(LKMAX,LIA,LJA)
   !  real(RP) :: V(LKMAX,LIA,LJA)
    real(RP) :: F1(LKMAX,LSIZE)
    real(RP) :: F2(LKMAX,LSIZE)
    real(RP) :: F3(LKMAX,LSIZE)
    real(RP) :: V(LKMAX,LSIZE)
    real(RP) :: TEMP2(LKMAX,LSIZE)
    real(RP) :: WATER2(LKMAX,LSIZE)
    real(RP) :: ICE2  (LKMAX,LSIZE)

#ifdef _OPENACC
    real(RP) :: work(LKMAX-1,4) ! for CR
#endif    

    real(RP) :: NDG_TEMP (LKMAX,LIA,LJA)
    real(RP) :: NDG_WATER(LKMAX,LIA,LJA)

    real(RP) :: MASS_total(LKMAX)
    real(RP) :: MASS_water(LKMAX)
    real(RP) :: MASS_ice(LKMAX)

    real(RP) :: ENGI(LKMAX,LSIZE)
    real(RP) :: CS
    real(RP) :: CL

    real(RP) :: flux(LKS-1:LKE,LSIZE)

    real(RP) :: ro, rw, ri
    real(RP) :: sw

    integer :: k, i, j
    integer :: ii
#if LSIZE == 1
    integer, parameter :: l = 1
#else
    integer  :: l
#endif    
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
    end if

    !$acc data copy(TEMP, WATER, ICE, TEMP1, WATER1) &
    !$acc      copyin(TEMP_t, WATER_t, ICE_t, WaterLimit, ThermalCond, HeatCapacity,  &
    !$acc             WaterDiff, SFLX_GH, SFLX_water, SFLX_RHOE, exists_land, CDZ   ) &
    !$acc      copyout(RUNOFF, RUNOFF_ENGI                                          ) &
    !$acc      create(NDG_TEMP, NDG_WATER)

    if ( LAND_DYN_BUCKET_nudging ) then
       if ( .not. replace ) then
          ! nudging is used

          !$omp parallel do
          !$acc kernels
          !$acc loop independent
          do j = LJS,LJE
          !$acc loop independent
          do i = LIS,LIE
          !$acc loop independent
          do k = LKS,LKE
             if ( TEMP1(k,i,j) == UNDEF ) then
                NDG_TEMP (k,i,j) = 0.0_RP
             else
                NDG_TEMP (k,i,j) = ( TEMP1 (k,i,j) - TEMP (k,i,j) ) / LAND_DYN_BUCKET_nudging_tausec * dt
             end if
          end do
          end do
          end do
          !$acc end kernels

          !$omp parallel do
          !$acc kernels
          !$acc loop independent
          do j = LJS,LJE
          !$acc loop independent
          do i = LIS,LIE
          !$acc loop independent
          do k = LKS,LKE
             if ( WATER1(k,i,j) == UNDEF ) then
                NDG_WATER(k,i,j) = 0.0_RP
             else
                NDG_WATER(k,i,j) = ( WATER1(k,i,j) - ( WATER(k,i,j) + ICE(k,i,j) ) ) / LAND_DYN_BUCKET_nudging_tausec * dt
             end if
          end do
          end do
          end do
          !$acc end kernels

       end if

       if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
          !$omp parallel do
          !$acc kernels
          !$acc loop independent
          do j = LJS, LJE
          !$acc loop independent
          do i = LIS, LIE
             NDG_WATER(LKE,i,j) = 0.0_RP
          end do
          end do
          !$acc end kernels
       end if

       if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
          !$omp parallel do
          !$acc kernels
          !$acc loop independent
          do j = LJS, LJE
          !$acc loop independent
          do i = LIS, LIE
             NDG_TEMP(LKE,i,j) = 0.0_RP
          end do
          end do
          !$acc end kernels
       end if

    else
       ! nudging is NOT used
      
       !$omp parallel do
       !$acc kernels
       !$acc loop independent
       do j = LJS,LJE
       !$acc loop independent
       do i = LIS,LIE
       !$acc loop independent
       do k = LKS,LKE
          NDG_TEMP (k,i,j) = 0.0_RP
          NDG_WATER(k,i,j) = 0.0_RP
       end do
       end do
       end do
       !$acc end kernels

    end if

    if ( .not. replace ) then

       !$omp parallel do private( i, ii, l, k, &
       !$omp MASS_total, MASS_water, MASS_ice, &
       !$omp F1, F2, F3, V, flux, kappa, CS, CL, sw )
       !$acc kernels
       !$acc loop independent
       do j = LJS, LJE
#if LSIZE == 1
       !$acc loop independent private( &
       !$acc MASS_total, MASS_water, MASS_ice, &
       !$acc TEMP2, WATER2, ICE2,              &
       !$acc F1, F2, F3, V, work, flux, kappa )
       do i = LIS, LIE
#else
       do ii = LIS, LIE, LSIZE
          do l = 1, LSIZE
             i = ii + l - 1
             if ( i > LIE ) exit
#endif
             if ( exists_land(i,j) ) then

                !$acc loop independent
                do k = LKS, LKE
                   MASS_total(k) = DWATR * WATER(k,i,j) + DICE * ICE(k,i,j)
                end do
                MASS_total(LKS) = MASS_total(LKS) + dt * SFLX_water(i,j) / CDZ(LKS)
  
                CS = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j)
                !$acc loop independent
                do k = LKS, LKE
                   ENGI(k,l) = ( CS + WATER_DENSCS * WATER(k,i,j) + ICE_DENSCS * ICE(k,i,j) ) * TEMP(k,i,j) - LHF * DICE * ICE(k,i,j)
                end do
                ENGI(LKS,l) = ENGI(LKS,l) + dt * ( SFLX_GH(i,j) + SFLX_RHOE(i,j) ) / CDZ(LKS)
  
                ! phase change
                !$acc loop independent
                do k = LKS, LKE
                   MASS_ice(k) = min( MASS_total(k), max( 0.0_RP, &
                       ( ENGI(k,l) - ( CS + CV_WATER * MASS_total(k) ) * LAND_DYN_BUCKET_T_frz ) &
                       / ( ( CV_ICE - CV_WATER ) * LAND_DYN_BUCKET_T_frz - LHF ) &
                       ) )
                   MASS_water(k) = MASS_total(k) - MASS_ice(k)
                   V(k,l) = MASS_water(k) / DWATR
                   ICE2(k,l) = MASS_ice(k) / DICE
                   TEMP2(k,l) = ( ENGI(k,l) + LHF * MASS_ice(k) ) &
                       / ( CS + CV_WATER * MASS_water(k) + CV_ICE * MASS_ice(k) )
                end do

                !--

                F3(LKS,l) = 0.0_RP
                F1(LKS,l) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
                F2(LKS,l) = 1.0_RP - F3(LKS,l) - F1(LKS,l)

                if ( LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
                   F3(LKE,l) = -2.0_RP * WaterDiff(i,j) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
                else
                   F3(LKE,l) = 0.0_RP
                end if

                F1(LKE,l) = 0.0_RP
                F2(LKE,l) = 1.0_RP - F3(LKE,l) - F1(LKE,l)

                do k = LKS+1, LKE-1
                   F3(k,l) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
                   F1(k,l) = -2.0_RP * WaterDiff(i,j) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
                   F2(k,l) = 1.0_RP - F3(k,l) - F1(k,i,j)
                end do

#if LSIZE == 1
                call MATRIX_SOLVER_tridiagonal_1D_CR( LKMAX, 1, LKMAX, &
#ifdef _OPENACC
                                                      work(:,:), &
#endif
                                                      F1(:,1), F2(:,1), F3(:,1), V(:,1), & ! [IN]
                                                      WATER2(:,1)                        ) ! [OUT]
#else
          end do

          call MATRIX_SOLVER_tridiagonal( LKMAX, 1, LKMAX, &
                                          F1(:,:), F2(:,:), F3(:,:), V(:,:), & ! [IN]
                                          WATER2(:,:)                        ) ! [IN]


          do l = 1, LSIZE
             i = ii + l - 1
             if ( i > l ) exit
#endif

                ! temperature

                flux(LKS-1,l) = 0.0_RP
                flux(LKE,l)   = 0.0_RP

                CS = ( 1.0_RP - WaterLimit(i,j) ) * HeatCapacity(i,j)
                !$acc loop independent
                do k = LKS, LKE
                   kappa(k,l) = ThermalCond(i,j) + 0.5_RP * WATER2(k,l)**(1.0_RP/3.0_RP)
                end do

                !$acc loop independent
                do k = LKS, LKE-1
                  flux(k,l) = - 2.0_RP *  DWATR * WaterDiff(i,j) * ( WATER2(k+1,l) - WATER2(k,l) ) / ( CDZ(k+1) + CDZ(k) )
                  sw = 0.5_RP - sign( 0.5_RP, flux(k,l) )
                  flux(k,l) = flux(k,l) * CV_WATER * ( TEMP2(k+1,l) * sw + TEMP2(k,l) * ( 1.0_RP - sw ) )
                end do
                if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
                   flux(LKE,l) = flux(LKE-1,l)
                end if

                !$acc loop independent
                do k = LKS, LKE
                   V(k,l) = ENGI(k,l) + LHF * DICE * ICE2(k,l) &
                          - dt * ( flux(k,l) - flux(k-1,l) ) / CDZ(k)
                end do

                CL = CS + WATER_DENSCS * WATER2(LKS,l) + ICE_DENSCS * ICE2(LKS,l)
                F3(LKS,l) = 0.0_RP
                F1(LKS,l) = - ( kappa(LKS,l) + kappa(LKS+1,l) ) / ( CDZ(LKS) * ( CDZ(LKS) + CDZ(LKS+1) ) ) * dt
                F2(LKS,l) = CL - F3(LKS,l) - F1(LKS,l)

                CL = CS + WATER_DENSCS * WATER2(LKE,l) + ICE_DENSCS * ICE2(LKE,l)
                if ( LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
                   F3(LKE,l) = - ( kappa(LKE,l) + kappa(LKE-1,l) ) / ( CDZ(LKE) * ( CDZ(LKE) + CDZ(LKE-1) ) ) * dt
                else
                   F3(LKE,l) = 0.0_RP
                end if
                F1(LKE,l) = 0.0_RP
                F2(LKE,l) = CL - F3(LKE,l) - F1(LKE,l)

                !$acc loop independent
                do k = LKS+1, LKE-1
                   CL = CS + WATER_DENSCS * WATER2(k,l) + ICE_DENSCS * ICE2(k,l)
                   F3(k,l) = - ( kappa(k,l) + kappa(k-1,l) ) / ( CDZ(k) * ( CDZ(k) + CDZ(k-1) ) ) * dt
                   F1(k,l) = - ( kappa(k,l) + kappa(k+1,l) ) / ( CDZ(k) * ( CDZ(k) + CDZ(k+1) ) ) * dt
                   F2(k,l) = CL - F3(k,l) - F1(k,l)
                end do


#if LSIZE == 1
                call MATRIX_SOLVER_tridiagonal_1D_CR( LKMAX, 1, LKMAX, &
#ifdef _OPENACC
                                                      work(:,:), &
#endif
                                                      F1(:,1), F2(:,1), F3(:,1), V(:,1), & ! [IN]
                                                      TEMP2(:,1)                         ) ! [OUT]
#else
         end do

         call MATRIX_SOLVER_tridiagonal( LKMAX, 1, LKMAX, &
                                         F1(:,:), F2(:,:), F3(:,:), V(:,:), & ! [IN]
                                         TEMP2(:,:)                         ) ! [IN]

         do l = 1, LSIZE
           i = ii + l - 1
           if ( i > l ) exit                          
#endif

             !$acc loop independent
             do k = LKS, LKE
                TEMP2(k,l) = TEMP2(k,l) + NDG_TEMP(k,i,j)

                if ( TEMP2(k,l) >= LAND_DYN_BUCKET_T_frz ) then
                   WATER2(k,l) = WATER2(k,l) + NDG_WATER(k,i,j)
                else
                   ICE2(k,l) = ICE2(k,l) + NDG_WATER(k,i,j)
                end if            
             end do

             ! runoff of soil moisture (vertical sum)
             RUNOFF(i,j)      = 0.0_RP
             RUNOFF_ENGI(i,j) = 0.0_RP

             !$acc loop independent
             do k = LKS, LKE
                ro = max( WATER2(k,l) + ICE2(k,l) - WaterLimit(i,j), 0.0_RP )
                rw = min( ro, WATER2(k,l) )
                ri = ro - rw
                WATER2(k,l) = WATER2(k,l) - rw
                ICE2(k,l) = ICE2(k,l) - ri
                rw = rw * DWATR / dt
                ri = ri * DICE / dt
                RUNOFF(i,j) = RUNOFF(i,j) + ( rw + ri ) * CDZ(k)
                RUNOFF_ENGI(i,j) = RUNOFF_ENGI(i,j) &
                   + ( ( rw * CV_WATER + ri * CV_ICE ) * TEMP2(k,l) - ri * LHF ) * CDZ(k)
             end do

             if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER ) then
                WATER2(LKE,l) = WATER(LKE,i,j)
                ICE2  (LKE,l) = ICE  (LKE,i,j)
             end if
             if ( .not. LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP ) then
                TEMP2(LKE,l) = TEMP(LKE,i,j)      
             end if
#if LSIZE == 1
#else
         end do
#endif
             !$acc loop independent
             do k = LKS, LKE
                TEMP (k,i,j) = TEMP2 (k,l)
                WATER(k,i,j) = WATER2(k,l)
                ICE  (k,i,j) = ICE2  (k,l)
             end do

          end if ! end if exists_land(i,j)
      end do ! end for j
      end do ! end for j

     else  
        ! if replace

        !$acc kernels
        !$acc loop independent
        do j = LJS, LJE
        !$acc loop independent
        do i = LIS, LIE
             !$acc loop independent
             do k = LKS, LKE
              TEMP (k,i,j) = TEMP1 (k,i,j)
              WATER(k,i,j) = WATER1(k,i,j)
              ICE  (k,i,j) = 0.0_RP
            end do          
        end do
        end do
        !$acc end kernels
     end if

    !$acc end data

    return
  end subroutine LAND_DYN_BUCKET

end module scale_land_dyn_bucket
