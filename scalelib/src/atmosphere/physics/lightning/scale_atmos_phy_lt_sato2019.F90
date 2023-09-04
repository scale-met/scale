!-------------------------------------------------------------------------------
!> module atmosphere / physics / lightninh / SATO2019
!!
!! @par Description
!!         Component for sato2019 tracer
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2019-03-19 (Y.Sato) [new] Newly created
!! @li      2021-02-22 (N. Yamashita, T. Iwashita) [add] Add preprocessing subroutine
!! @li      2021-06-24 (T. Iwashita) [add] Add OpenMP option for preprocessing subroutine
!!
!<
!-------------------------------------------------------------------------------
#ifndef COLORING
#ifdef _OPENACC
#define COLORING 2
#elif defined(_OPENMP)
#define COLORING 1
#else
#define COLORING 0
#endif
#endif

#include "scalelib.h"
module scale_atmos_phy_lt_sato2019
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_LT_sato2019_setup
  public :: ATMOS_PHY_LT_sato2019_finalize
  public :: ATMOS_PHY_LT_sato2019_adjustment
  public :: ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT
!  public :: ATMOS_PHY_LT_sato2019_mkinit
!  public :: ATMOS_PHY_LT_electric_field
!  public :: ATMOS_PHY_LT_neutralization_MG2001
!  public :: ATMOS_PHY_LT_neutralization_F2013
!  public :: ATMOS_PHY_LT_judge_absE
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
  !--- for Bi-CGSTAB to solve the poisson equation
  integer,  private                 :: ITMAX = 3000
  real(RP), private                 :: epsilon = 0.1_RP ** (RP*2)
  character(len=64),private         :: NUTR_TYPE='F2013'
  character(len=64),private,save    :: NUTR_qhyd='Each_POLARITY2'
  integer,  private                 :: MAX_NUTR = 100
  real(RP), private                 :: Eint = 150.0E+3_RP     ! [V/m]
  real(RP), private                 :: delEint = 10.0E+3_RP   ! [V/m]
  real(RP), private                 :: Estop = 15.0E+3_RP     ! [V/m]
  real(RP), private                 :: qrho_chan = 0.5_RP     ! [nC/m3]
  real(RP), private                 :: qrho_neut = 0.5_RP     ! [nC/m3]
  real(RP), private                 :: fp = 0.3_RP            ! [-]
  real(RP), private                 :: zcg = 0.0_RP           ! lowest grid point of lightning path for MG2001 [m]
  integer,  private                 :: NUTR_ITMAX = 1000
  integer,  private                 :: FLAG_preprocessing = 2 !> 0: Disable
                                                              !> 1: Gauss-Seidel
                                                              !> 2: Synmetric Gauss-Sidel (Default)
                                                              !> 3: Incomplete LU factorization
  integer,  private                 :: KIJMAXG
  character(len=H_LONG)             :: ATMOS_PHY_LT_LUT_FILENAME !--- LUT file name
  character(len=H_LONG)             :: fname_lut_lt="LUT_TK1978_v.txt" !--- LUT file name
  integer, save                     :: fid_lut_lt
  real(RP), private                 :: R_neut = 2000.0_RP  ! [m]
  real(RP), private                 :: rho0 = 1.225_RP       ! [kg/m3]
  real(RP), save                    :: flg_eint_hgt
  logical,  private                 :: LT_DO_Lightning = .true.

  !--- for history output
  real(RP), allocatable, private    :: d_QCRG_TOT(:,:,:)
  real(RP), allocatable, private    :: LT_PATH_TOT(:,:,:,:)
  real(RP), allocatable, private    :: fls_int_p_tot(:,:,:)
  real(RP), allocatable, private    :: B_F2013_TOT(:,:)
  real(RP), allocatable, private    :: G_F2013(:,:)
  real(RP),              private    :: C_F2013

  real(RP), allocatable, private    :: A(:,:,:,:) !--- A : Laplasian (Coefficient matrix)
  !$acc declare create(d_QCRG_TOT,LT_PATH_TOT,fls_int_p_tot,B_F2013_TOT,G_F2013,A)

  !---
  integer,  parameter, private :: nxlut_lt = 200, nylut_lt = 200
  real(RP), private :: dq_chrg( nxlut_lt,nylut_lt )    !--- charge separation [fC]
  real(RP), private :: grid_lut_t( nxlut_lt )
  real(RP), private :: grid_lut_l( nylut_lt )
!  !$acc declare create(dq_chrg,grid_lut_t,grid_lut_l)
  real(RP), private :: tcrglimit
  logical, private                  :: Hgt_dependency_Eint = .false.

  !--- Index for local variable
  integer, private, parameter :: I_lt_x = 1
  integer, private, parameter :: I_lt_y = 2
  integer, private, parameter :: I_lt_z = 3
  integer, private, parameter :: I_lt_abs = 4

  !--- For history output
  integer, private, parameter :: w_nmax = 11
  integer, private, parameter :: I_Ex = 1
  integer, private, parameter :: I_Ey = 2
  integer, private, parameter :: I_Ez = 3
  integer, private, parameter :: I_Eabs = 4
  integer, private, parameter :: I_Epot = 5
  integer, private, parameter :: I_Qneut = 6
  integer, private, parameter :: I_LTpath = 7
  integer, private, parameter :: I_PosFLASH = 8
  integer, private, parameter :: I_NegFLASH = 9
  integer, private, parameter :: I_FlashPoint = 10
  integer, private, parameter :: I_FOD = 11
  integer,  private              :: HIST_id(w_nmax)
  character(len=H_SHORT), private :: w_name(w_nmax)
  character(len=H_MID),   private :: w_longname(w_nmax)
  character(len=H_SHORT), private :: w_unit(w_nmax)
  data w_name / 'Ex', &
                'Ey', &
                'Ez', &
                'Eabs', &
                'Epot', &
                'Qneut', &
                'LTpath', &
                'PosFLASH', &
                'NegFLASH', &
                'FlashPoint', &
                'FOD' /
  data w_longname / &
                'X component of Electrical Field', &
                'Y component of Electrical Field', &
                'Z component of Electrical Field', &
                'Absolute value of Electrical Field', &
                'Electric Potential', &
                'Cumulative Neutralizated charge', &
                'Cumulative Number of flash path', &
                'Cumulative Number of Positive flash', &
                'Cumulative Number of Negative flash', &
                'Cumulative Number of Flash point', &
                'Flash Origin Density' /
  data w_unit / &
                'kV/m', &
                'kV/m', &
                'kV/m', &
                'kV/m', &
                'V', &
                'nC/m3/s', &
                'num/grid/s', &
                'num/grid/s', &
                'num/grid/s', &
                'num/grid/s', &
                'num/grid/s' /
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_LT_sato2019_setup( KA, KS, KE, &
                                          IA, IS, IE, &
                                          JA, JS, JE, &
                                          IMAXG,      &
                                          JMAXG,      &
                                          KMAX,       &
                                          MP_TYPE,    &
                                          CDX, CDY    )
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster, &
       PRC_myrank
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ
    use scale_const, only: &
       T00 => CONST_TEM00, &
       PI  => CONST_PI
    use scale_file_history, only: &
       FILE_HISTORY_reg
    use scale_atmos_grid_cartesC, only: &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_atmos_grid_cartesC_metric, only: &
       MAPF  => ATMOS_GRID_CARTESC_METRIC_MAPF, &
       J13G  => ATMOS_GRID_CARTESC_METRIC_J13G, &
       J23G  => ATMOS_GRID_CARTESC_METRIC_J23G, &
       J33G  => ATMOS_GRID_CARTESC_METRIC_J33G
    use scale_atmos_grid_cartesC_index, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY, &
       I_UY, &
       I_XV, &
       I_UV
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer, intent(in)  :: IA, IS, IE
    integer, intent(in)  :: JA, JS, JE
    integer, intent(in)  :: IMAXG
    integer, intent(in)  :: JMAXG
    integer, intent(in)  :: KMAX

    character(len=*), intent(in) :: MP_TYPE
    real(RP),         intent(in)  :: CDX(IA)
    real(RP),         intent(in)  :: CDY(JA)

    integer :: n, myu, ip
    integer :: ierr
    integer :: i, j, k

    namelist / PARAM_ATMOS_PHY_LT_SATO2019 / &
         NUTR_TYPE, &
         ATMOS_PHY_LT_LUT_FILENAME, &
         ITMAX, &
         EPSILON, &
         Eint, &
         delEint, &
         Estop, &
         qrho_chan, &
         qrho_neut, &
         fp, &
         NUTR_ITMAX, &
         MAX_NUTR, &
         NUTR_qhyd, &
         zcg, &
         R_neut, &
         Hgt_dependency_Eint, &
         LT_DO_Lightning, &
         FLAG_preprocessing

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_LT_SATO2019,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_LT_sato2019_setup",*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_LT_SATO2019. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_LT_SATO2019)

    !--- If zcg is lower than lowest grid point, zcg set lowest grid point
    if( zcg <= CZ(KS) ) then
       zcg = CZ(KS)
    endif

    if( R_neut <= minval( CDX,1 ) .or. R_neut <= minval( CDY,1 ) ) then
        R_neut = 2.d0 * min( minval( CDX,1 ), minval( CDY,1 ) )
    endif

    !$acc enter data create(dq_chrg,grid_lut_t,grid_lut_l)
    if ( PRC_IsMaster ) then
        fname_lut_lt = ATMOS_PHY_LT_LUT_FILENAME
        fid_lut_lt = IO_get_available_fid()
        !--- open parameter of cloud microphysics
        open ( fid_lut_lt, file = fname_lut_lt, form = 'formatted', status = 'old', iostat=ierr )

        if ( ierr == 0 ) then
          LOG_INFO("ATMOS_PHY_LT_sato2019_setup",'(2A)') 'Read LUT of TK78 table from ', trim(fname_lut_lt)
          read( fid_lut_lt,* )
          do n   = 1, nylut_lt
          do myu = 1, nxlut_lt
             read( fid_lut_lt,* ) grid_lut_t( myu ), grid_lut_l( n ), dq_chrg( myu,n )
          enddo
          enddo
          !$acc update device(dq_chrg,grid_lut_t,grid_lut_l)
        !--- LUT file does not exist
        else
           LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'xxx LUT for LT is requied when ATMOS_PHY_LT_TYPE = SATO2019, stop!'
           call PRC_abort
        endif

        if( NUTR_TYPE /= 'MG2001' .and. NUTR_TYPE /= 'F2013' ) then
           LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'xxx NUTR_TYPE should be MG2001 or F2013 stop!'
           call PRC_abort
        endif
    endif

    call COMM_bcast( nxlut_lt, nylut_lt, dq_chrg )
    call COMM_bcast( nxlut_lt, grid_lut_t )
    call COMM_bcast( nylut_lt, grid_lut_l )

!    KIJMAXG = (IEG-ISG+1)*(JEG-JSG+1)*(KE-KS+1)
    KIJMAXG = IMAXG*JMAXG*KMAX

    ! for history output
    allocate( d_QCRG_TOT(KA,IA,JA) )
    allocate( LT_PATH_TOT(KA,IA,JA,3) )
    allocate( fls_int_p_tot(KA,IA,JA) )
    d_QCRG_TOT(:,:,:) = 0.0_RP
    LT_PATH_TOT(:,:,:,:) = 0.0_RP
    fls_int_p_tot(:,:,:) = 0.0_RP
    !$acc update device(d_QCRG_TOT,LT_PATH_TOT,fls_int_p_tot)

    tcrglimit = -60.0_RP+T00

    if( NUTR_TYPE == 'F2013' ) then
      LOG_INFO("ATMOS_PHY_LT_sato2019_setup",'(A,F15.7,A)') 'Radius of neutralization is ', R_neut, "[m]"
    endif
    allocate( B_F2013_TOT(IA,JA) )
    allocate( G_F2013(IA,JA) )
    B_F2013_TOT(:,:) = 0.0_RP
    do i = 1, IA
    do j = 1, JA
      G_F2013(i,j) = CDX(i)*CDY(j)*1.0E-6_RP       ! [m2] -> [km2]
    enddo
    enddo
    !$acc update device(B_F2013_TOT,G_F2013)
    C_F2013 = PI*R_neut*R_neut*1.0E-6_RP  ! [m2] -> [km2]

    flg_eint_hgt = 0.0_RP
    if( Hgt_dependency_Eint .and. NUTR_TYPE == 'MG2001') then
      flg_eint_hgt = 1.0_RP
    endif

    if( NUTR_qhyd /= 'TOTAL' .and. NUTR_qhyd /= 'Each_POLARITY' .and. NUTR_qhyd /= 'Each_POLARITY2' ) then
      LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'xxx NUTR_qhyd should be TOTAL, or Each_POLARITY, or Each_POLARITY2, stop!'
      call PRC_abort
    endif

!    if( ( NUTR_qhyd == 'Each_POLARITY' .or. NUTR_qhyd /= 'Each_POLARITY2' ) .and. MP_TYPE == 'SUZUKI10' ) then
!       LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'NUTR_qhyd = Each_POLARITY, Each_POLARITY2 is not supported for MP_TYPE SUZUKI10'
!       call PRC_abort
!    endif

    do ip = 1, w_nmax
       if( ip /= I_FOD ) then
          call FILE_HISTORY_reg( w_name(ip), w_longname(ip), w_unit(ip), & ! [IN]
                                 HIST_id(ip)                             ) ! [OUT]
       elseif( ip == I_FOD ) then
          call FILE_HISTORY_reg( w_name(ip), w_longname(ip), w_unit(ip), & ! [IN]
                                 HIST_id(ip), dim_type='XY'              ) ! [OUT]
       endif
    end do

    allocate( A(KA,15,IA,JA) )
    !---- input vector A
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! (k,i,j)
       A(k,1,i,j) = &
                  - MAPF(i,j  ,1,I_XY)*MAPF(i,j  ,1,I_XY)*RCDX(i  )*RFDX(i) &
                  - MAPF(i,j  ,1,I_XY)*MAPF(i,j  ,1,I_XY)*RCDX(i-1)*RFDX(i) &
                  + MAPF(i,j  ,1,I_XY)*J13G(k,i  ,j,I_UYZ)*f2h(k  ,i,j,2)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - MAPF(i,j  ,1,I_XY)*J13G(k,i  ,j,I_UYZ)*f2h(k-1,i,j,1)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - MAPF(i,j  ,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k  ,i,j,2)*0.5_RP*RFDX(i)*RFDZ(k) &
                  + MAPF(i,j  ,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i,j,1)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*J13G(k  ,i,j,I_XYW)*RCDZ(k)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*J13G(k-1,i,j,I_XYW)*RCDZ(k)*RFDZ(k) &
                  - MAPF(i,j  ,2,I_XY)*MAPF(i,j  ,2,I_XY)*RCDY(j  )*RFDY(j) &
                  - MAPF(i,j  ,2,I_XY)*MAPF(i,j  ,2,I_XY)*RCDY(j-1)*RFDY(j) &
                  + MAPF(i,j  ,2,I_XY)*J23G(k,i,j  ,I_XVZ)*f2h(k  ,i,j,2)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - MAPF(i,j  ,2,I_XY)*J23G(k,i,j  ,I_XVZ)*f2h(k-1,i,j,1)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - MAPF(i,j  ,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k  ,i,j,2)*0.5_RP*RFDY(j)*RFDZ(k) &
                  + MAPF(i,j  ,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j,1)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*J23G(k  ,i,j,I_XYW)*RCDZ(k)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*J23G(k-1,i,j,I_XYW)*RCDZ(k)*RFDZ(k) &
                  - J33G*J33G*RFDZ(k)*RCDZ(k) &
                  - J33G*J33G*RFDZ(k)*RCDZ(k-1)

       ! (k-1,i,j)
       A(k,2,i,j) = &
                  - MAPF(i,j,1,I_XY)*J13G(k,i  ,j,I_UYZ)*f2h(k-1,i,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*J13G(k-1,i,j,I_XYW)*RFDZ(k)*RCDZ(k-1) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j  ,I_XVZ)*f2h(k-1,i,j,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*J23G(k-1,i,j,I_XYW)*RFDZ(k)*RCDZ(k-1) &
                  + J33G*J33G*RFDZ(k)*RCDZ(k-1)

       ! (k+1,i,j)
       A(k,3,i,j) = &
                    MAPF(i,j,1,I_XY)*J13G(k,i  ,j,I_UYZ)*f2h(k,i,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k,i,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*J13G(k,i,j,I_XYW)*RFDZ(k)*RCDZ(k) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j  ,I_XVZ)*f2h(k,i,j,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k,i,j,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*J23G(k,i,j,I_XYW)*RFDZ(k)*RCDZ(k) &
                  + J33G*J33G*RFDZ(k)*RCDZ(k)

       ! (k,i-1,j)
       A(k,4,i,j) = &
                    MAPF(i,j,1,I_XY)*MAPF(i-1,j,1,I_XY)*RFDX(i)*RCDX(i-1) &
                  - MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k  ,i-1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i-1,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k  ,i-1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i-1,j,1)*0.50_RP*RFDX(i)*RFDZ(k)

       ! (k,i+1,j)
       A(k,5,i,j) = &
                    MAPF(i,j,1,I_XY)*MAPF(i+1,j,1,I_XY)*RFDX(i)*RCDX(i) &
                  + MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k  ,i+1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k-1,i+1,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k  ,i+1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i+1,j,1)*0.50_RP*RFDX(i)*RFDZ(k)

       ! (k,i,j-1)
       A(k,6,i,j) = &
                    MAPF(i,j,2,I_XY)*MAPF(i,j-1,2,I_XY)*RFDY(j)*RCDY(j-1) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k  ,i,j-1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j-1,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k  ,i,j-1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j-1,1)*0.50_RP*RFDY(j)*RFDZ(k)

       ! (k,i,j+1)
       A(k,7,i,j) = &
                    MAPF(i,j,2,I_XY)*MAPF(i,j+1,2,I_XY)*RFDY(j)*RCDY(j) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k  ,i,j+1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k-1,i,j+1,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k  ,i,j+1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j+1,1)*0.50_RP*RFDY(j)*RFDZ(k)

       ! (k-1,i-1,j)
       A(k,8,i,j) = &
                    MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i-1,j,2)*0.5_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i-1,j,2)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k-1,i+1,j)
       A(k,9,i,j) = &
                  - MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k-1,i+1,j,2)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i+1,j,2)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k-1,i,j-1)
       A(k,10,i,j) = &
                    MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j-1,2)*0.5_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j-1,2)*0.5_RP*RFDY(j)*RFDZ(k)

       ! (k-1,i,j+1)
       A(k,11,i,j) = &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k-1,i,j+1,2)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j+1,2)*0.5_RP*RFDY(j)*RFDZ(k)

       ! (k+1,i-1,j)
       A(k,12,i,j) = &
                  - MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k,i-1,j,1)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k,i-1,j,1)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k+1,i+1,j)
       A(k,13,i,j) = &
                    MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k,i+1,j,1)*0.5_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k,i+1,j,1)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k+1,i,j-1)
       A(k,14,i,j) = &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k,i,j-1,1)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k,i,j-1,1)*0.5_RP*RFDY(j)*RFDZ(k)

       ! (k+1,i,j+1)
       A(k,15,i,j) = &
                    MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k,i,j+1,1)*0.5_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k,i,j+1,1)*0.5_RP*RFDY(j)*RFDZ(k)

    enddo
    enddo
    enddo
    !$acc update device(A)

    return
  end subroutine ATMOS_PHY_LT_sato2019_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_LT_sato2019_finalize

    !$acc exit data delete(dq_chrg,grid_lut_t,grid_lut_l)

    deallocate( d_QCRG_TOT )
    deallocate( LT_PATH_TOT )
    deallocate( fls_int_p_tot )

    deallocate( B_F2013_TOT )
    deallocate( G_F2013 )

    return
  end subroutine ATMOS_PHY_LT_sato2019_finalize

  !-----------------------------------------------------------------------------
  !> Update of charge density
  subroutine ATMOS_PHY_LT_sato2019_adjustment( &
       KA, KS, KE,   &
       IA, IS, IE,   &
       JA, JS, JE,   &
       KIJMAX,   &
       IMAX,     &
       JMAX,     &
       QA_LT,    &
       DENS,     &
       RHOT,     &
       QHYD,     &
       Sarea,    &
       dt_LT,    &
       QTRC,     &
       Epot      )
    use scale_const, only: &
       SMALL => CONST_EPS
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
    use scale_comm_cartesC, only: &
       COMM_datatype, &
       COMM_world, &
       COMM_wait, &
       COMM_vars8
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put
    implicit none

    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE
    integer,  intent(in) :: KIJMAX
    integer,  intent(in) :: IMAX
    integer,  intent(in) :: JMAX
!    character(len=H_SHORT)
    integer,  intent(in) :: QA_LT
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QHYD(KA,IA,JA)
    real(RP), intent(in) :: Sarea(KA,IA,JA,QA_LT)       !--- Surface area and that of each catergory [m2]
    real(DP), intent(in) :: dt_LT
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_LT)
    real(RP), intent(inout) :: Epot(KA,IA,JA) !--- Electrical potential at previous time step [V]

    real(RP) :: QHYD_mass(KA,IA,JA)         !--- Mass of total hydrometeor [kg/m3]
    real(RP) :: QCRG(KA,IA,JA)              !--- Total charge density [nC/m3]
    real(RP) :: Efield(KA,IA,JA,I_lt_abs)   !--- Electrical field (1-3 ->x,y,z, 4->abs. )
    real(RP) :: NUM_end(KA,IA,JA,3)         !--- Number of each flash type (1->negative, 2->ground, 3->positive)
    real(RP) :: d_QCRG(KA,IA,JA)            !--- Change of charge by charge neutralization [fC/m3]
    real(RP) :: LT_PATH(KA,IA,JA)           !--- Lightning path (0-> no path, 1-> path)
    real(RP) :: fls_int_p(KA,IA,JA)
    real(RP) :: Total_Sarea1, Total_Sarea2  !--- Sum of surface area and that of each category [m2]
    real(RP) :: r_totalSarea1, r_totalSarea2
    real(RP) :: neg_crg, pos_crg
    real(RP) :: frac
    real(RP) :: dqneut(KA,IA,JA,QA_LT)
    real(RP) :: dqneut_real(KA,IA,JA,QA_LT)
    real(RP) :: dqneut_real_tot(KA,IA,JA)
    logical  :: flg_chrged(QA_LT)
    real(RP) :: Emax, Emax_old
    logical  :: output_step
    integer  :: flg_lt_neut
    integer  :: i, j, k, m, n, countbin, ip
    real(RP) :: sw, zerosw, positive, negative
    integer  :: count_neut

    logical  :: HIST_sw(w_nmax)
    real(RP) :: w3d(KA,IA,JA)

    real(RP) :: diff_qcrg(0:1), lack(0:1), sum_crg(0:1)
    real(RP) :: crg_rate(QA_LT), qcrg_before(QA_LT)
    real(RP) :: B_F2013(IA,JA)
    integer  :: int_sw
    real(RP) :: iprod, buf, tmp_qcrg, tmp_dqcrg
    integer  :: ierror

    !$acc data &
    !$acc copyin(DENS,RHOT,QHYD,Sarea) &
    !$acc copy(QTRC,Epot) &
    !$acc create(QHYD_mass,QCRG,Efield,NUM_end,d_QCRG,LT_PATH,fls_int_p, &
    !$acc        dqneut,dqneut_real,dqneut_real_tot,flg_chrged,HIST_sw,w3d, &
    !$acc        diff_qcrg,lack,sum_crg,crg_rate,qcrg_before,B_F2013)

    !$acc kernels
    NUM_end(:,:,:,:) = 0.0_RP
    B_F2013(:,:) = 0.0_RP
    !$acc end kernels

    !$omp parallel do private(tmp_qcrg)
    !$acc kernels
    !$acc loop collapse(2)
    do j = JS, JE
    do i = IS, IE

       ! calc total charge density
       do k = KS, KE
          tmp_qcrg = 0.0_RP
          !$acc loop seq
          do n = 1, QA_LT
             tmp_qcrg = tmp_qcrg + QTRC(k,i,j,n)
          end do
          QCRG(k,i,j) = tmp_qcrg * DENS(k,i,j) * 1.E-6_RP ![fC/kg] -> [nc/m3]
       enddo

       do k = KS, KE
          QHYD_mass(k,i,j) = QHYD(k,i,j) * DENS(k,i,j) ![kg/kg] -> [kg/m3]
       enddo

    enddo
    enddo
    !$acc end kernels

    iprod = 0.0_RP
    !$omp parallel do reduction(+:iprod) private(i,j,k)
    !$acc kernels
    !$acc loop collapse(3) reduction(+:iprod)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
        iprod = iprod + abs( QCRG(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    call MPI_AllReduce(iprod, buf, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

    if( buf <= SMALL ) then

       !$omp parallel do
       !$acc kernels
       !$acc loop collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Epot(k,i,j) = 0.0_RP
          Efield(k,i,j,I_lt_x) = 0.0_RP
          Efield(k,i,j,I_lt_y) = 0.0_RP
          Efield(k,i,j,I_lt_z) = 0.0_RP
       end do
       end do
       end do
       !$acc end kernels

    else

       !--- Calculate E field
       call ATMOS_PHY_LT_electric_field( KA, KS, KE,                   &   ! [IN]
                                         IA, IS, IE,                   &   ! [IN]
                                         JA, JS, JE,                   &   ! [IN]
                                         QCRG    (:,:,:),              &   ! [IN]
                                         DENS    (:,:,:),              &   ! [IN]
                                         RHOT    (:,:,:),              &   ! [IN]
                                         Epot    (:,:,:),              &   ! [INOUT]
                                         Efield  (:,:,:,I_lt_x:I_lt_z) )   ! [INOUT]

    endif

    !$omp parallel do
    !$acc kernels
    !$acc loop collapse(3)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Efield(k,i,j,I_lt_abs) = sqrt( Efield(k,i,j,I_lt_x)*Efield(k,i,j,I_lt_x) &
                                    + Efield(k,i,j,I_lt_y)*Efield(k,i,j,I_lt_y) &
                                    + Efield(k,i,j,I_lt_z)*Efield(k,i,j,I_lt_z) )

       LT_PATH(k,i,j) = 0.0_RP

    enddo
    enddo
    enddo
    !$acc end kernels

    if ( LT_DO_Lightning ) then

       !$acc kernels
       d_QCRG_TOT(:,:,:) = 0.0_RP
       fls_int_p_tot(:,:,:) = 0.0_RP
       LT_PATH_TOT(:,:,:,:) = 0.0_RP
       B_F2013_TOT(:,:) = 0.0_RP
       !$acc end kernels

       call ATMOS_PHY_LT_judge_absE( KA, KS, KE,             &   ! [IN]
                                     IA, IS, IE,             &   ! [IN]
                                     JA, JS, JE,             &   ! [IN]
                                     DENS(:,:,:),            &   ! [IN]
                                     Efield(:,:,:,I_lt_abs), &   ! [IN]
                                     Emax,                   &   ! [OUT]
                                     flg_lt_neut             )   ! [OUT]

       count_neut = 0
       do while( flg_lt_neut > 0 )

         Emax_old = Emax

         call COMM_vars8( Efield   (:,:,:,I_lt_x),   1 )
         call COMM_vars8( Efield   (:,:,:,I_lt_y),   2 )
         call COMM_vars8( Efield   (:,:,:,I_lt_z),   3 )
         call COMM_vars8( Efield   (:,:,:,I_lt_abs), 4 )
         call COMM_vars8( QHYD_mass(:,:,:),          5 )
         call COMM_vars8( QCRG     (:,:,:),          6 )
         call COMM_vars8( Epot     (:,:,:),          7 )
         call COMM_wait ( Efield   (:,:,:,I_lt_x),   1 )
         call COMM_wait ( Efield   (:,:,:,I_lt_y),   2 )
         call COMM_wait ( Efield   (:,:,:,I_lt_z),   3 )
         call COMM_wait ( Efield   (:,:,:,I_lt_abs), 4 )
         call COMM_wait ( QHYD_mass(:,:,:),          5 )
         call COMM_wait ( QCRG     (:,:,:),          6 )
         call COMM_wait ( Epot     (:,:,:),          7 )

         !--- Calculate lightning path and charge neutralization
         if( NUTR_TYPE == 'MG2001' ) then
           !$acc update host(Efield,Epot,QCRG,QHYD_mass,NUM_end,LT_PATH)
           call ATMOS_PHY_LT_neutralization_MG2001(              &
                                             KA, KS, KE,         & !  [IN]
                                             IA, IS, IE,         & !  [IN]
                                             JA, JS, JE,         & !  [IN]
                                             KIJMAX, IMAX, JMAX, & !  [IN]
                                             Efield   (:,:,:,:), & !  [IN]
                                             Epot     (:,:,:),   & !  [IN]
                                             DENS     (:,:,:),   & !  [IN]
                                             QCRG     (:,:,:),   & !  [IN]
                                             QHYD_mass(:,:,:),   & !  [IN]
                                             NUM_end  (:,:,:,:), & !  [INOUT]
                                             LT_PATH  (:,:,:),   & !  [INOUT]
                                             fls_int_p(:,:,:),   & !  [OUT]
                                             d_QCRG   (:,:,:)    ) !  [OUT]
           !$acc update device(NUM_end,LT_PATH,fls_int_p,d_QCRG)
         elseif( NUTR_TYPE == 'F2013' ) then
           call ATMOS_PHY_LT_neutralization_F2013(               &
                                             KA, KS, KE,         & !  [IN]
                                             IA, IS, IE,         & !  [IN]
                                             JA, JS, JE,         & !  [IN]
                                             KIJMAX,             & !  [IN]
                                             Efield   (:,:,:,:), & !  [IN]
                                             Epot     (:,:,:),   & !  [IN]
                                             DENS     (:,:,:),   & !  [IN]
                                             QCRG     (:,:,:),   & !  [IN]
                                             QHYD_mass(:,:,:),   & !  [IN]
                                             NUM_end  (:,:,:,:), & !  [INOUT]
                                             LT_PATH  (:,:,:),   & !  [INOUT]
                                             fls_int_p(:,:,:),   & !  [OUT]
                                             d_QCRG   (:,:,:),   & !  [OUT]
                                             B_F2013  (:,:)      ) !  [OUT]
         endif

         call COMM_vars8( LT_path(:,:,:),1 )
         call COMM_vars8( d_QCRG(:,:,:),2 )
         call COMM_wait ( LT_path(:,:,:),1 )
         call COMM_wait ( d_QCRG(:,:,:),2 )

         !$acc kernels
         dqneut(:,:,:,:) = 0.0_RP
         dqneut_real(:,:,:,:) = 0.0_RP
         !$acc end kernels

         !-- Calculate neutralization of each hydrometeor or each category
         select case( NUTR_qhyd )
         case ( 'TOTAL' )

            !$omp parallel do &
            !$omp private(Total_Sarea1,Total_Sarea2,r_totalSarea1,r_totalSarea2,zerosw)
            !$acc kernels
            !$acc loop collapse(3)
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               if( abs( d_QCRG(k,i,j) ) > 0.0_RP ) then
                  Total_Sarea1 = 0.0_RP
                  !$acc loop seq
                  do n = 1, QA_LT
                     Total_Sarea1 = Total_Sarea1 + Sarea(k,i,j,n)
                  enddo
                  zerosw = 0.5_RP - sign( 0.5_RP, Total_Sarea1-SMALL )
                  r_totalSarea1 = 1.0_RP / ( Total_Sarea1 + zerosw ) * ( 1.0_RP - zerosw )
                  !$acc loop seq
                  do n = 1, QA_LT
                     dqneut(k,i,j,n) = d_QCRG(k,i,j)*1.0E+6_RP &
                                     * Sarea(k,i,j,n) * r_totalSarea1 / DENS(k,i,j)
                     QTRC(k,i,j,n) = QTRC(k,i,j,n) + dqneut(k,i,j,n)
                     dqneut_real(k,i,j,n) = dqneut(k,i,j,n)
                  enddo
               endif
            enddo
            enddo
            enddo
            !$acc end kernels

         case ( 'Each_POLARITY', 'Each_POLARITY2' )

            !$omp parallel do &
            !$omp private(Total_Sarea1,Total_Sarea2,r_totalSarea1,r_totalSarea2,flg_chrged,pos_crg,neg_crg,frac,lack, &
            !$omp         diff_qcrg,crg_rate,sum_crg,qcrg_before, &
            !$omp         positive,negative,zerosw,sw,int_sw)
            !$acc parallel
            !$acc loop collapse(2) gang
            do j = JS, JE
            do i = IS, IE
            !$acc loop vector private(lack,flg_chrged,diff_qcrg,crg_rate,sum_crg)
            do k = KS, KE
               lack(:) = 0.0_RP
               if( abs( d_QCRG(k,i,j) ) > 0.0_RP ) then

                  !--- flg whether the charged or not (0.0-> not charged, 1.0->charged)
                  !$acc loop seq
                  do n = 1, QA_LT
                     flg_chrged(n) = abs(QTRC(k,i,j,n)) >= SMALL
                  enddo

                  Total_Sarea1 = 0.0_RP
                  Total_Sarea2 = 0.0_RP
                  pos_crg = 0.0_RP
                  neg_crg = 0.0_RP
                  !$acc loop seq
                  do n = 1, QA_LT
                     if ( flg_chrged(n) ) then
                        positive = 0.5_RP + sign( 0.5_RP, QTRC(k,i,j,n) )
                        negative = 1.0_RP - positive
                        !--- total of positive charge
                        pos_crg = pos_crg + QTRC(k,i,j,n) * positive
                        !--- total of negative charge
                        neg_crg = neg_crg + QTRC(k,i,j,n) * negative
                        !--- Sarea of positively charged hydrometeor
                        Total_Sarea1 = Total_Sarea1 + Sarea(k,i,j,n) * positive
                        !--- Sarea of negatively charged hydrometeor
                        Total_Sarea2 = Total_Sarea2 + Sarea(k,i,j,n) * negative
                     end if
                  end do

                  zerosw = 0.5_RP - sign( 0.5_RP, abs( QCRG(k,i,j) ) - SMALL )
                  frac = d_QCRG(k,i,j) / ( QCRG(k,i,j) + zerosw ) * ( 1.0_RP - zerosw )
                  pos_crg = frac * pos_crg
                  neg_crg = frac * neg_crg

                  !--- remove 0 surface area ( no crg for each porality )
                  zerosw = 0.5_RP - sign( 0.5_RP, Total_Sarea1 - SMALL )
                  r_totalSarea1 = 1.0_RP / ( Total_Sarea1 + zerosw ) * ( 1.0_RP - zerosw )
                  zerosw = 0.5_RP - sign( 0.5_RP, Total_Sarea2 - SMALL )
                  r_totalSarea2 = 1.0_RP / ( Total_Sarea2 + zerosw ) * ( 1.0_RP - zerosw )

                  diff_qcrg(:) = 0.0_RP
                  !$acc loop seq
                  do n = 1, QA_LT
                     crg_rate(n) = 0.0_RP
                  end do
                  sum_crg(:) = 0.0_RP
                  !$acc loop seq
                  do n = 1, QA_LT
                     if ( flg_chrged(n) ) then
                        sw = 0.5_RP + sign( 0.5_RP, QTRC(k,i,j,n) )
                        dqneut(k,i,j,n) = &
                                    + pos_crg * Sarea(k,i,j,n) * r_totalSarea1 &
                                    * sw &   ! sw = 1 positive
                                    + neg_crg * Sarea(k,i,j,n) * r_totalSarea2  &
                                    * ( 1.0_RP - sw ) ! sw = 0 negative
                        qcrg_before(n) = QTRC(k,i,j,n)

                        if( sw == 1.0_RP ) then
                          QTRC(k,i,j,n) = max( QTRC(k,i,j,n) + dqneut(k,i,j,n), 0.0_RP )  !--- limiter
                        elseif( sw == 0.0_RP ) then
                          QTRC(k,i,j,n) = min( QTRC(k,i,j,n) + dqneut(k,i,j,n), 0.0_RP )  !--- limiter
                        endif

                        int_sw = int( sw )   ! 0-> negative, 1-> positive
                        diff_qcrg(int_sw) = diff_qcrg(int_sw) &
                                          + ( QTRC(k,i,j,n) - qcrg_before(n) )
                        sum_crg(int_sw) = sum_crg(int_sw) + QTRC(k,i,j,n)
                     end if
                  enddo

                  if( NUTR_qhyd == 'Each_POLARITY2' ) then ! Adjust
                    lack(0) = neg_crg - diff_qcrg(0) ! negative (Should be Positive)
                    lack(1) = pos_crg - diff_qcrg(1) ! positive (Should be Negative)
#ifdef DEBUG
                    if( lack(0) > 0.0_RP .and. abs(lack(0)/diff_qcrg(0)) > 1.0E-10_RP ) then
                        LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(A,4E15.7)') &
                          "Large negative for lack(0) ", lack(0)/diff_qcrg(0), neg_crg, lack(0), diff_qcrg(0)
                    endif
                    if( lack(1) < 0.0_RP .and. abs(lack(1)/diff_qcrg(1)) > 1.0E-10_RP ) then
                       LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(A,4E15.7)') &
                          "Large positive for lack(1) ", lack(1)/diff_qcrg(1), pos_crg, lack(1), diff_qcrg(1)
                    endif
#endif
                    lack(0) = max( lack(0), 0.0_RP ) ! negative (Should be Positive)
                    lack(1) = min( lack(1), 0.0_RP ) ! positive (Should be Negative)
                    !$acc loop seq
                    do n = 1, QA_LT
                       if ( flg_chrged(n) ) then
                          sw = 0.5_RP + sign( 0.5_RP, QTRC(k,i,j,n) )
                          int_sw = int( sw )   ! 0-> negative, 1-> positive
                          if( sum_crg(int_sw) /= 0.0_RP ) then
                             crg_rate(n) = QTRC(k,i,j,n)/sum_crg(int_sw)
                          else
                             crg_rate(n) = 0.0_RP
                          endif
                          QTRC(k,i,j,n) = QTRC(k,i,j,n) + crg_rate(n) * lack(int_sw)
                       endif
                    enddo
                  endif

                  !$acc loop seq
                  do n = 1, QA_LT
                     if ( flg_chrged(n) ) then
                        dqneut_real(k,i,j,n) = QTRC(k,i,j,n) - qcrg_before(n)
                     endif
                  enddo

               endif

            enddo
            enddo
            enddo
            !$acc end parallel

         end select

         !$acc kernels
         !$acc loop collapse(3)
         do j = JS, JE
         do i = IS, IE

            ! calc total charge density
            do k = KS, KE
               tmp_qcrg = 0.0_RP
               tmp_dqcrg = 0.0_RP
               !$acc loop seq
               do n = 1, QA_LT
                  tmp_qcrg = tmp_qcrg + QTRC(k,i,j,n)
                  tmp_dqcrg = tmp_dqcrg + dqneut_real(k,i,j,n)
               end do
               QCRG(k,i,j) = tmp_qcrg * DENS(k,i,j) * 1.E-6_RP ![fC/kg] -> [nc/m3]
               dqneut_real_tot(k,i,j) = tmp_dqcrg
            enddo

         enddo
         enddo
         !$acc end kernels

         iprod = 0.0_RP
         !$omp parallel do reduction(+:iprod) private(i,j,k)
         !$acc kernels
         !$acc loop collapse(3) reduction(+:iprod)
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
             iprod = iprod + abs( QCRG(k,i,j) )
         enddo
         enddo
         enddo
         !$acc end kernels

         call MPI_AllReduce(iprod, buf, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

         if( buf <= SMALL ) then

            !$omp parallel do
            !$acc kernels
            !$acc loop collapse(3)
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               Epot(k,i,j) = 0.0_RP
               Efield(k,i,j,I_lt_x) = 0.0_RP
               Efield(k,i,j,I_lt_y) = 0.0_RP
               Efield(k,i,j,I_lt_z) = 0.0_RP
            end do
            end do
            end do
            !$acc end kernels

         else

            !--- Calculate E field
            call ATMOS_PHY_LT_electric_field( KA, KS, KE,                   & ! [IN]
                                              IA, IS, IE,                   & ! [IN]
                                              JA, JS, JE,                   & ! [IN]
                                              QCRG    (:,:,:),              & ! [IN]
                                              DENS    (:,:,:),              & ! [IN]
                                              RHOT    (:,:,:),              & ! [IN]
                                              Epot    (:,:,:),              & ! [INOUT]
                                              Efield  (:,:,:,I_lt_x:I_lt_z) ) ! [INOUT]

         endif

         !--- Add Total number of path
         !$omp parallel do
         !$acc kernels
         !$acc loop collapse(3)
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            Efield(k,i,j,I_lt_abs) = sqrt( Efield(k,i,j,I_lt_x)*Efield(k,i,j,I_lt_x) &
                                         + Efield(k,i,j,I_lt_y)*Efield(k,i,j,I_lt_y) &
                                         + Efield(k,i,j,I_lt_z)*Efield(k,i,j,I_lt_z) )
         end do
         end do
         end do
         !$acc end kernels

         !--- Add Total number of charge neutralization and flash point
         if ( HIST_id(I_Qneut) > 0 ) then
            !$omp parallel do
            !$acc kernels
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
                d_QCRG_TOT(k,i,j) = d_QCRG_TOT(k,i,j) + dqneut_real_tot(k,i,j)*1.E-6_RP ![fC/m3]->[nC/m3]
            end do
            end do
            end do
            !$acc end kernels
         end if
         if ( HIST_id(I_FlashPoint) > 0 ) then
            !$omp parallel do
            !$acc kernels
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               fls_int_p_tot(k,i,j) = fls_int_p_tot(k,i,j) + fls_int_p(k,i,j)
            end do
            end do
            end do
            !$acc end kernels
         end if

         if ( HIST_id(I_PosFLASH) > 0 ) then
            !$omp parallel do
            !$acc kernels
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               LT_PATH_TOT(k,i,j,1) = LT_PATH_TOT(k,i,j,1) &
                                    + 0.5_RP + sign( 0.5_RP,-dqneut_real_tot(k,i,j)-SMALL )
            end do
            end do
            end do
            !$acc end kernels
         end if
         if ( HIST_id(I_NegFLASH) > 0 ) then
            !$omp parallel do
            !$acc kernels
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               LT_PATH_TOT(k,i,j,2) = LT_PATH_TOT(k,i,j,2) &
                                    + 0.5_RP + sign( 0.5_RP, dqneut_real_tot(k,i,j)-SMALL )
            end do
            end do
            end do
            !$acc end kernels
         end if
         if ( HIST_id(I_LTpath) > 0 ) then
            !$omp parallel do
            !$acc kernels
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               LT_PATH_TOT(k,i,j,3) = LT_PATH_TOT(k,i,j,3) + LT_PATH(k,i,j)
            end do
            end do
            end do
            !$acc end kernels
         end if
         if ( HIST_id(I_FOD) > 0 ) then
            !$acc kernels
            do j = JS, JE
            do i = IS, IE
               B_F2013_TOT(i,j) = B_F2013_TOT(i,j) + G_F2013(i,j)/C_F2013*B_F2013(i,j)
            enddo
            enddo
            !$acc end kernels
         endif

         call ATMOS_PHY_LT_judge_absE( KA, KS, KE,             &   ! [IN]
                                       IA, IS, IE,             &   ! [IN]
                                       JA, JS, JE,             &   ! [IN]
                                       DENS(:,:,:),            &   ! [IN]
                                       Efield(:,:,:,I_lt_abs), &   ! [IN]
                                       Emax,                   &   ! [OUT]
                                       flg_lt_neut             )   ! [OUT]

         count_neut = count_neut + 1
#ifdef DEBUG
         LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(A,F15.7,A,F15.7,1X,I0)')  &
                   'CHECK', &
                    Emax_old*1.E-3_RP, ' [kV/m] -> ', Emax*1.E-3_RP, count_neut
#endif

         if( flg_lt_neut == 1 .and. Emax == Emax_old ) then
            flg_lt_neut = 0
            if( PRC_IsMaster ) then
               LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(A,2(E15.7,A),1X,I0)')  &
                   'Eabs value after neutralization is same as previous value, Finish', &
                    Emax_old*1.E-3_RP, ' [kV/m] -> ', Emax*1.E-3_RP, count_neut
            endif
         elseif( flg_lt_neut == 1 .and. count_neut == MAX_NUTR ) then
            flg_lt_neut = 0
            if( PRC_IsMaster ) then
               LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(2(E15.7,A),1X,I0)')  &
                   Emax_old*1.E-3_RP, ' [kV/m] -> ', Emax*1.E-3_RP, &
                   ' [kV/m], reach maximum neutralization count, Finish', count_neut
            endif
         elseif( flg_lt_neut == 1 .and. Emax > Emax_old ) then
            flg_lt_neut = 0
            if( PRC_IsMaster ) then
               LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(2(E15.7,A),1X,I0)')  &
                   Emax_old*1.E-3_RP, ' [kV/m] -> ', Emax*1.E-3_RP, &
                   ' [kV/m] larger than previous one by neutralization, back to previous step, Finish', &
                   count_neut
            endif
            !---- Back to Charge density as previous step
            !$acc kernels
            !$acc loop collapse(3)
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               !$acc loop seq
               do n = 1, QA_LT
                  QTRC(k,i,j,n) = QTRC(k,i,j,n) - dqneut_real(k,i,j,n)
               enddo
               d_QCRG_TOT(k,i,j) = d_QCRG_TOT(k,i,j) - dqneut_real_tot(k,i,j)*1.E-6_RP
               d_QCRG(k,i,j) = 0.0_RP
            enddo
            enddo
            enddo
            !$acc end kernels
         elseif( flg_lt_neut == 2 ) then
            flg_lt_neut = 0
            if( PRC_IsMaster ) then
               LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(2(E15.7,A),1X,I0)')  &
                   Emax_old*1.E-3_RP, ' [kV/m] -> ', Emax*1.E-3_RP, &
                   '[kV/m] After neutralization, Finish' , count_neut
            endif
         elseif( flg_lt_neut == 0 ) then
            if( PRC_IsMaster ) then
               LOG_INFO("ATMOS_PHY_LT_sato2019_adjustment",'(2(F15.7,A),1X,I0)')  &
                   Emax_old*1.E-3_RP, ' [kV/m] -> ', Emax*1.E-3_RP, &
                   '[kV/m] After neutralization, Finish', count_neut
            endif
         endif

       enddo

    else

       !$omp parallel do private(i,j,k)
       !$acc kernels
       !$acc loop collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          d_QCRG(k,i,j) = 0.0_RP
       end do
       end do
       end do
       !$acc end kernels

    endif

    !--- For history output
    do ip = 1, w_nmax
       call FILE_HISTORY_query( HIST_id(ip), HIST_sw(ip) )
    end do

    if ( HIST_sw(I_Ex  ) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       !$acc loop collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_x  )*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_Ex  ), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Ey  ) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       !$acc loop collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_y  )*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_Ey  ), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Ez  ) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       !$acc loop collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_z  )*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_Ez  ), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Eabs) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       !$acc loop collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_abs)*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_Eabs), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Epot) ) then
       call FILE_HISTORY_put( HIST_id(I_Epot), Epot(:,:,:) )
    end if
    if ( HIST_sw(I_Qneut) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = d_QCRG_TOT(k,i,j)/dt_LT ![nC/m3/s]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_Qneut), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_LTpath)  ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = LT_PATH_TOT(k,i,j,3)/dt_LT ![num/grid/s]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_LTpath), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_PosFLASH) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = LT_PATH_TOT(k,i,j,1)/dt_LT ![num/grid/s]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_PosFLASH), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_NegFLASH) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = LT_PATH_TOT(k,i,j,2)/dt_LT ![num/grid/s]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_NegFLASH), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_FlashPoint) ) then
       !$omp parallel do private(i,j,k)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = fls_int_p_tot(k,i,j)/dt_LT ![num/grid/s]
       enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_FlashPoint), w3d(:,:,:) )
    end if
    if ( HIST_sw(I_FOD) ) then
       !$omp parallel do private(i,j)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          w3d(1,i,j) = B_F2013_TOT(i,j)/dt_LT ![num/grid/s]
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_put( HIST_id(I_FOD), w3d(1,:,:) )
    end if

    !$acc end data

    return
  end subroutine ATMOS_PHY_LT_sato2019_adjustment
  !-----------------------------------------------------------------------------
  !> calculate electric field from charge density of each grid
  !> temporaly Bi-CGSTAB is used in this component
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_LT_electric_field( &
       KA, KS, KE, & ! [IN]
       IA, IS, IE, & ! [IN]
       JA, JS, JE, & ! [IN]
       QCRG,       & ! [IN]
       DENS,       & ! [IN]
       RHOT,       & ! [IN]
       E_pot,      & ! [INOUT]
       Efield      ) ! [INOUT]
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster, &
       PRC_myrank
    use scale_const, only: &
       EPS    => CONST_EPS, &
       EPSvac => CONST_EPSvac, &
       EPSair => CONST_EPSair
    use scale_atmos_grid_cartesC, only: &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_atmos_grid_cartesC_metric, only: &
       MAPF  => ATMOS_GRID_CARTESC_METRIC_MAPF, &
       GSQRT => ATMOS_GRID_CARTESC_METRIC_GSQRT, &
       J13G  => ATMOS_GRID_CARTESC_METRIC_J13G, &
       J23G  => ATMOS_GRID_CARTESC_METRIC_J23G, &
       J33G  => ATMOS_GRID_CARTESC_METRIC_J33G
    use scale_atmos_grid_cartesC_index, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY, &
       I_UY, &
       I_XV, &
       I_UV
    use scale_prc_cartesC, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_comm_cartesC, only: &
       COMM_datatype, &
       COMM_world, &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: QCRG     (KA,IA,JA)      !-- Charge density [nC/m3]
    real(RP), intent(in)  :: DENS     (KA,IA,JA)      !-- Total density [kg/m3]
    real(RP), intent(in)  :: RHOT     (KA,IA,JA)      !-- density weighted potential temperature [K kg/m3]
    real(RP), intent(inout) :: E_pot (KA,IA,JA)       !-- Electric potential [V]
    real(RP), intent(inout) :: Efield(KA,IA,JA,3)     !-- Electric field [V/m]

    real(RP) :: eps_air
    !--- A x E_pott = - QCRG/epsiron
    real(RP) :: B(KA,IA,JA)               !--- B : -QCRG*DENS/epsiron
    real(RP) :: E_pot_N(KA,IA,JA)         !--- electrical potential calculated by Bi-CGSTAB

    integer :: i, j, k, ijk, ierror
    real(RP) :: iprod, buf

    call PROF_rapstart('LT_E_field', 2)


    !$acc data &
    !$acc copy(E_pot,Efield) &
    !$acc copyin(QCRG,DENS,RHOT) &
    !$acc create(B,E_pot_N) &
    !$acc copyin(GSQRT,MAPF,J13G,J23G,RCDX,RCDY,RFDZ)

    !$omp parallel do private(i,j,k)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       E_pot_N(k,i,j) = E_pot(k,i,j)   !-- initial value -> previous step value
    end do
    end do
    end do
    !$acc end kernels

    call COMM_vars8( E_pot_N, 1 )


    !$omp parallel do private(i,j,k,eps_air)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       eps_air = EPSvac * EPSair   !--- temporary, dependency of epsiron on P and T will be implemented
       B(k,i,j) = - QCRG(k,i,j)/eps_air * 1.0E-9_RP ! [nC/m3] -> [C/m3 * m/F] = [V/m2]
    enddo
    enddo
    enddo
    !$acc end kernels

    call COMM_wait ( E_pot_N, 1 )

    !--- calcuclate counter matrix
    call ATMOS_PHY_LT_solve_bicgstab( &
       KA, KS, KE, & ! (in)
       IA, IS, IE, & ! (in)
       JA, JS, JE, & ! (in)
       E_pot,      & ! (out)
       E_pot_N,    & ! (in)
       A, B        ) ! (in)

    call COMM_vars8( E_pot, 1 )
    call COMM_wait ( E_pot, 1, .true. )

    !$omp parallel do private(i,j)
    !$acc parallel vector_length(32)
    !$acc loop collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = 1, KS-1
          E_pot(k,i,j) = 0.0_RP
       enddo
       do k = KE+1, KA
          E_pot(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo
    !$acc end parallel

    !---- Calculate Electrical Field
    !$omp parallel do private(i,j,k)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Efield(k,i,j,I_lt_x) = - MAPF(i,j,1,I_XYZ)/GSQRT(k,i,j,I_XYZ) * &
                            ( &
                            ( GSQRT(k,i+1,j,I_XYZ)*E_pot(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*E_pot(k,i-1,j) ) &
                            * RCDX(i) * 0.5_RP &
                            + ( J13G(k+1,i,j,I_XYZ)*GSQRT(k+1,i,j,I_XYZ)*E_pot(k+1,i,j) &
                              - J13G(k-1,i,j,I_XYZ)*GSQRT(k-1,i,j,I_XYZ)*E_pot(k-1,i,j) &
                              ) &
                            * RFDZ(k) * 0.5_RP  &
                            )
       Efield(k,i,j,I_lt_y) = - MAPF(i,j,2,I_XYZ)/GSQRT(k,i,j,I_XYZ) * &
                            ( &
                            ( GSQRT(k,i,j+1,I_XYZ)*E_pot(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*E_pot(k,i,j-1) ) &
                            * RCDY(j) * 0.5_RP &
                            + ( J23G(k+1,i,j,I_XYZ)*GSQRT(k+1,i,j,I_XYZ)*E_pot(k+1,i,j) &
                              - J23G(k-1,i,j,I_XYZ)*GSQRT(k-1,i,j,I_XYZ)*E_pot(k-1,i,j) &
                              ) &
                            * RFDZ(k) * 0.5_RP  &
                            )
       Efield(k,i,j,I_lt_z) = - 1.0_RP/GSQRT(k,i,j,I_XYZ) * &
                            ( J33G*E_pot(k+1,i ,j  )*GSQRT(k+1,i,j,I_XYZ) &
                            - J33G*E_pot(k-1,i ,j  )*GSQRT(k-1,i,j,I_XYZ) &
                            ) &
                            * RFDZ(k) * 0.5_RP
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    call PROF_rapend('LT_E_field', 2)

    return
  end subroutine ATMOS_PHY_LT_electric_field

  subroutine ATMOS_PHY_LT_solve_bicgstab( &
       KA, KS, KE, & ! [IN]
       IA, IS, IE, & ! [IN]
       JA, JS, JE, & ! [IN]
       PHI_N,      & ! [OUT]
       PHI,        & ! [IN]
       M, B        ) ! [IN]
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster
    use scale_comm_cartesC, only: &
       COMM_datatype, &
       COMM_world, &
       COMM_vars8, &
       COMM_wait
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(out) :: PHI_N(KA,IA,JA)
    real(RP), intent(in)  :: PHI(KA,IA,JA)
    real(RP), intent(in)  :: M(KA,15,IA,JA)
    real(RP), intent(in)  :: B(KA,IA,JA)

    real(RP) :: r0(KA,IA,JA)

    real(RP) :: p(KA,IA,JA)
    real(RP) :: Mp(KA,IA,JA)
    real(RP) :: s(KA,IA,JA)
    real(RP) :: Ms(KA,IA,JA)
    real(RP) :: al, be, w

    real(RP), pointer :: r(:,:,:), rn(:,:,:), swap(:,:,:)
    real(RP), target :: rbuf1(KA,IA,JA)
    real(RP), target :: rbuf2(KA,IA,JA)
    real(RP) :: v1(KA,IA,JA)

    real(RP) :: r0r
    real(RP) :: norm, error, error2

    real(RP) :: iprod1, iprod2
    real(RP) :: buf(2)

    real(RP):: diag(KA,IA,JA)
    real(RP):: z1(KA,IA,JA)
    real(RP):: z2(KA,IA,JA)
    real(RP):: Mz1(KA,IA,JA)
    real(RP):: Mz2(KA,IA,JA)

    integer :: k, i, j
    integer :: iis, iie, jjs, jje
    integer :: iter
    integer :: ierror

    !$acc data &
    !$acc copyout(PHI_N) &
    !$acc copyin(PHI,M,B) &
    !$acc create(r0,p,Mp,s,Ms,rbuf1,rbuf2,v1,diag,z1,z2,Mz1,Mz2)

    r  => rbuf1
    rn => rbuf2

    if( FLAG_preprocessing == 3 ) then
       !$acc update host(M)
       call ILU_decomp(KA, KS, KE, &
                       IA, IS, IE, &
                       JA, JS, JE, &
                       M,  diag)
       !$acc update device(diag)
    endif

    call mul_matrix( KA, KS, KE, & ! (in)
                     IA, IS, IE, & ! (in)
                     JA, JS, JE, & ! (in)
                     v1, M, PHI  ) ! v1 = M x0

    norm = 0.0_RP
    !$omp parallel do reduction(+:norm) private(i, j, k)
    !$acc kernels
    !$acc loop collapse(3) reduction(+:norm)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       norm = norm + B(k,i,j)**2
    enddo
    enddo
    enddo
    !$acc end kernels

    ! r = b - M x0
    !$omp parallel do private(i, j, k)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       r(k,i,j) = B(k,i,j) - v1(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do private(i, j, k)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       r0(k,i,j) = r(k,i,j)
       p(k,i,j) = r(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    r0r  = 0.0_RP
    !$omp parallel do reduction(+:r0r) private(i, j, k)
    !$acc kernels
    !$acc loop collapse(3) reduction(+:r0r)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       r0r = r0r + r0(k,i,j) * r(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do private(i, j, k)
    !$acc kernels
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       PHI_N(k,i,j) = PHI(k,i,j)
    end do
    end do
    end do
    !$acc end kernels

    buf(1) = r0r
    buf(2) = norm
    call MPI_AllReduce(MPI_IN_PLACE, buf(:), 2, COMM_datatype, MPI_SUM, COMM_world, ierror)
    r0r = buf(1)
    norm = buf(2)
    error2 = norm

    do iter = 1, ITMAX

       call COMM_vars8( p, 1 )

       error = 0.0_RP
       !$omp parallel do reduction(+:error) private(i, j, k)
       !$acc kernels
       !$acc loop collapse(3) reduction(+:error)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          error = error + r(k,i,j)**2
       enddo
       enddo
       enddo
       !$acc end kernels

       call MPI_AllReduce(MPI_IN_PLACE, error, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

       call COMM_wait ( p, 1 )

       if ( sqrt(error/norm) < epsilon ) then
         LOG_INFO("ATMOS_PHY_LT_Efield",'(a,1x,i0,1x,2e15.7)') "Bi-CGSTAB converged:", iter, sqrt(error/norm),norm
         exit
       endif
       error2 = error

       if( FLAG_preprocessing == 0 ) then  !-- No preprocessing

          call mul_matrix( KA, KS, KE, & ! (in)
                           IA, IS, IE, & ! (in)
                           JA, JS, JE, & ! (in)
                           Mp, M, p    )

          iprod1 = 0.0_RP
          !$omp parallel do reduction(+:iprod1) private(i, j, k)
          !$acc kernels
          !$acc loop collapse(3) reduction(+:iprod1)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             iprod1 = iprod1 + r0(k,i,j) * Mp(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

       else   !--- Preprocessing

          if( FLAG_preprocessing == 1 ) then !--- Gauss-Seidel preprocessing

             call gs( KA, KS, KE, & ! (in)
                      IA, IS, IE, & ! (in)
                      JA, JS, JE, & ! (in)
                      z1, M, p    )

          elseif( FLAG_preprocessing == 2 ) then  !--- Synmetric Gauss-Seidel preprocessing (Default)

             call sgs( KA, KS, KE, & ! (in)
                       IA, IS, IE, & ! (in)
                       JA, JS, JE, & ! (in)
                       z1, M, p    )

          elseif( FLAG_preprocessing == 3 ) then  !--- Incomplete Cholesky Factorization preprocessing

             !$acc update host(p,diag) ! M is already updated
             call solve_ILU( KA, KS, KE, & ! (in)
                             IA, IS, IE, & ! (in)
                             JA, JS, JE, & ! (in)
                             z1, M, p, diag)
             !$acc update device(z1)

          endif

          call COMM_vars8( z1, 1 )

          call mul_matrix( KA, KS, KE, & ! (in)
                           IA, IS, IE, & ! (in)
                           JA, JS, JE, & ! (in)
                           Mp, M, p    )

          call COMM_wait ( z1, 1 )

          call mul_matrix( KA, KS, KE, & ! (in)
                           IA, IS, IE, & ! (in)
                           JA, JS, JE, & ! (in)
                           Mz1, M, z1    )

          iprod1 = 0.0_RP
          !$omp parallel do reduction(+:iprod1)  private(i, j, k)
          !$acc kernels
          !$acc loop collapse(3) reduction(+:iprod1)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             iprod1 = iprod1 + r0(k,i,j) * Mz1(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

       endif

       call MPI_AllReduce(MPI_IN_PLACE, iprod1, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

       if ( iprod1 == 0.0_RP ) then
         LOG_INFO("ATMOS_PHY_LT_Efield",'(a,1x,e15.7,1x,i10)') 'Iprod1 is zero(Bi-CGSTAB) skip:', iprod1, iter
         exit
       endif
       al = r0r / iprod1 ! (r0,r) / (r0,Mp)

       if( FLAG_preprocessing == 0 ) then  !-- No preprocessing
          !$omp parallel do
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             s(k,i,j) = r(k,i,j) - al*Mp(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       else  ! Preprocessing
          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             s(k,i,j) = r(k,i,j) - al*Mz1(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       call COMM_vars8( s, 1 )
       call COMM_wait ( s, 1 )
       if( FLAG_preprocessing == 0 ) then  !--- No Preprocessing

          call mul_matrix( KA, KS, KE, & ! (in)
                           IA, IS, IE, & ! (in)
                           JA, JS, JE, & ! (in)
                           Ms, M,  s   )

          iprod1 = 0.0_RP
          iprod2 = 0.0_RP
          !$omp parallel do reduction(+:iprod1,iprod2)
          !$acc kernels
          !$acc loop collapse(3) reduction(+:iprod1,iprod2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             iprod1 = iprod1 + Ms(k,i,j) *  s(k,i,j)
             iprod2 = iprod2 + Ms(k,i,j) * Ms(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

       else  !--- Preprocessing

          if( FLAG_preprocessing == 1 ) then   !--- Gauss-Seidel preprocessing

             call gs( KA, KS, KE, & ! (in)
                      IA, IS, IE, & ! (in)
                      JA, JS, JE, & ! (in)
                      z2, M, s    )

          elseif( FLAG_preprocessing == 2 ) then  !--- Synmetric Gauss-Seidel preprocessing (Default)

             call sgs( KA, IS, KE, & ! (in)
                       IA, IS, IE, & ! (in)
                       JA, JS, JE, & ! (in)
                       z2, M, s    )

          elseif( FLAG_preprocessing == 3 ) then  !--- Incomplete Cholesky Factorization preprocessing

             !$acc update host(s,diag) ! M is already updated
             call solve_ILU( KA, KS, KE, & ! (in)
                             IA, IS, IE, & ! (in)
                             JA, JS, JE, & ! (in)
                             z2, M, s, diag)
             !$acc update device(z2)
          endif

          call COMM_vars8( z2, 1 )
          call COMM_wait ( z2, 1 )
          call mul_matrix( KA, KS, KE, & ! (in)
                           IA, IS, IE, & ! (in)
                           JA, JS, JE, & ! (in)
                           Mz2, M, z2 )

          iprod1 = 0.0_RP
          iprod2 = 0.0_RP
          !$omp parallel do reduction(+:iprod1,iprod2) private(i, j, k)
          !$acc kernels
          !$acc loop collapse(3) reduction(+:iprod1,iprod2)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             iprod1 = iprod1 + Mz2(k,i,j) *  s(k,i,j)
             iprod2 = iprod2 + Mz2(k,i,j) * Mz2(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

       endif

       buf(1) = iprod1
       buf(2) = iprod2
       call MPI_AllReduce(MPI_IN_PLACE, buf(:), 2, COMM_datatype, MPI_SUM, COMM_world, ierror)

       if ( buf(2) == 0.0_RP ) then
         LOG_INFO("ATMOS_PHY_LT_Efield",'(a,1x,e15.7,1x,i10)') 'Buf(2) is zero(Bi-CGSTAB) skip:', buf(2), iter
         exit
       endif
       w = buf(1) / buf(2) ! (Ms,s) / (Ms,Ms)

       if( FLAG_preprocessing == 0 ) then !--- No preprocessing

          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             PHI_N(k,i,j) = PHI_N(k,i,j) + al*p(k,i,j) + w*s(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             rn(k,i,j) = s(k,i,j) - w*Ms(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

       else

          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             PHI_N(k,i,j) = PHI_N(k,i,j) + al*z1(k,i,j) + w*z2(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             rn(k,i,j) = s(k,i,j) - w*Mz2(k,i,j)
          enddo
          enddo
          enddo
          !$acc end kernels

       endif

       iprod1 = 0.0_RP
       !$omp parallel do reduction(+:iprod1) private(i, j, k)
       !$acc kernels
       !$acc loop collapse(3) reduction(+:iprod1)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          iprod1 = iprod1 + r0(k,i,j) * rn(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels

       be = al/w / r0r

       call MPI_AllReduce(iprod1, r0r, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

       be = be * r0r ! al/w * (r0,rn)/(r0,r)

       if( FLAG_preprocessing == 0 ) then !--- No preprocessing
          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             p(k,i,j) = rn(k,i,j) + be * ( p(k,i,j) - w*Mp(k,i,j) )
          enddo
          enddo
          enddo
          !$acc end kernels
       else
          !$omp parallel do private(i, j, k)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             p(k,i,j) = rn(k,i,j) + be * ( p(k,i,j) - w*Mz1(k,i,j) )
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       swap => rn
       rn   => r
       r    => swap

       if ( r0r == 0.0_RP ) then
         LOG_INFO("ATMOS_PHY_LT_Efield",'(a,1x,i0,1x,3e15.7)') "Inner product of r0 and r_itr is zero(Bi-CGSTAB) :", &
                                                                iter, r0r, sqrt(error/norm), norm
         exit
       endif

    enddo

    if ( iter >= ITMAX ) then
       if( PRC_IsMaster ) then
         LOG_WARN("ATMOS_PHY_LT_solve_bicgstab",'(a,1x,2e15.7)') 'Bi-CGSTAB not converged:', error, norm
         LOG_WARN_CONT('(a,1x,2e15.7)') 'Bi-CGSTAB not converged:', epsilon, sqrt(error/norm)
         LOG_WARN_CONT('(a,1x,2e15.7)') 'xxx epsilon(set,last)=', epsilon, sqrt(error/norm)
         if( error /= error ) then
          LOG_ERROR("ATMOS_PHY_LT_solve_bicgstab",*) 'xxx error or norm is NaN Stop!'
          call PRC_abort
         endif
       endif
    endif

    !$acc end data
    return
  end subroutine ATMOS_PHY_LT_solve_bicgstab

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine ILU_decomp(KA, KS, KE, &
                        IA, IS, IE, &
                        JA, JS, JE, &
                        M,  diag)
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: M(KA,15,IA,JA)
    real(RP), intent(out) :: diag(KA,IA,JA)

    integer :: k, i, j

    !$omp parallel do private(i, j, k)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       diag(k,i,j)=M(k,1,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i, j, k)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       if(k /= KE) then
          diag(k+1,i,j) = diag(k+1,i,j) + M(k,3,i,j)*M(k+1,2,i,j)/diag(k,i,j)
          if(i /= IE) then
             diag(k+1,i+1,j) = diag(k+1,i+1,j) + M(k,13,i,j)*M(k+1,8,i+1,j)/diag(k,i,j)
          endif
          if(j /= JE) then
             diag(k+1,i,j+1) = diag(k+1,i,j+1) + M(k,15,i,j)*M(k+1,10,i,j+1)/diag(k,i,j)
          endif
       endif
       if(i /= IE) then
          diag(k,i+1,j) = diag(k,i+1,j) +M(k,5,i,j)*M(k,4,i+1,j)/diag(k,i,j)
          if(k /= KS) then
             diag(k-1,i+1,j) = diag(k-1,i+1,j) + M(k,9,i,j)*M(k-1,12,i+1,j)/diag(k,i,j)
          endif
       endif
       if(j /= JE) then
          diag(k,i,j+1)=diag(k,i,j+1) + M(k,7,i,j)*M(k,6,i,j+1)/diag(k,i,j)
          if(k /= KS) then
             diag(k-1,i,j+1) = diag(k-1,i,j+1) + M(k,11,i,j)*M(k-1,14,i,j+1)/diag(k,i,j)
          endif
       endif
    enddo
    enddo
    enddo

    return
  end subroutine ILU_decomp

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine solve_ILU( KA,KS,KE, &
                        IA,IS,IE, &
                        JA,JS,JE, &
                        Z, M, V, diag)
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: M(KA,15,IA,JA)
    real(RP), intent(in)  :: diag(KA,IA,JA)
    real(RP), intent(in)  :: V(KA,IA,JA)
    real(RP), intent(out) :: Z(KA,IA,JA)

    real(RP):: Y(KA,IA,JA)
    integer :: k,i,j

    call gs_ILU(KA, KS, KE, &
                IA, IS, IE, &
                JA, JS, JE, &
                Y,  M,  V,  diag)

    !$omp parallel do
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
      Y(k,i,j)=diag(k,i,j)*Y(k,i,j)
    enddo
    enddo
    enddo

    call back_sub_ILU(KA, KS, KE, &
                      IA, IS, IE, &
                      JA, JS, JE, &
                      Z,  M,  Y, diag)

    return

  end subroutine solve_ILU

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine back_sub_ILU( KA, KS, KE, &
                           IA, IS, IE, &
                           JA, JS, JE, &
                           Z,  M,  V,diag)
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster
   implicit none
   integer,  intent(in)  :: KA, KS, KE
   integer,  intent(in)  :: IA, IS, IE
   integer,  intent(in)  :: JA, JS, JE
   real(RP), intent(in)  :: V(KA,IA,JA)
   real(RP), intent(in)  :: M(KA,15,IA,JA)
   real(RP), intent(in)  :: diag(KA,IA,JA)
   real(RP), intent(out) :: Z(KA,IA,JA)

   integer :: k, i, j

   Z(:,:,:)=0.0_RP

   do j = JE, JS,-1
   do i = IE, IS,-1

         Z(KE,i,j) = ( V(KE,i,j) &
         -( M(KE,2,i,j) * Z(KE-1,i  ,j  ) &
          + M(KE,4,i,j) * Z(KE  ,i-1,j  ) &
          + M(KE,5,i,j) * Z(KE  ,i+1,j  ) &
          + M(KE,6,i,j) * Z(KE  ,i  ,j-1) &
          + M(KE,7,i,j) * Z(KE  ,i  ,j+1) &
          + M(KE,8,i,j) * Z(KE-1,i-1,j  ) &
          + M(KE,9,i,j) * Z(KE-1,i+1,j  ) &
          + M(KE,10,i,j)* Z(KE-1,i  ,j-1) &
          + M(KE,11,i,j)* Z(KE-1,i  ,j+1) ) )/diag(KE,i,j)

         do k = KE-1, KS+1,-1
            Z(k,i,j) = (V(k,i,j) &
            -( M(k,2,i,j) * Z(k-1,i  ,j  ) &
             + M(k,3,i,j) * Z(k+1,i  ,j  ) &
             + M(k,4,i,j) * Z(k  ,i-1,j  ) &
             + M(k,5,i,j) * Z(k  ,i+1,j  ) &
             + M(k,6,i,j) * Z(k  ,i  ,j-1) &
             + M(k,7,i,j) * Z(k  ,i  ,j+1) &
             + M(k,8,i,j) * Z(k-1,i-1,j  ) &
             + M(k,9,i,j) * Z(k-1,i+1,j  ) &
             + M(k,10,i,j)* Z(k-1,i  ,j-1) &
             + M(k,11,i,j)* Z(k-1,i  ,j+1) &
             + M(k,12,i,j)* Z(k+1,i-1,j  ) &
             + M(k,13,i,j)* Z(k+1,i+1,j  ) &
             + M(k,14,i,j)* Z(k+1,i  ,j-1) &
             + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/diag(k,i,j)
         enddo

         Z(KS,i,j) = (V(KS,i,j) &
         -( M(KS,3,i,j) * Z(KS+1,i  ,j  ) &
          + M(KS,4,i,j) * Z(KS  ,i-1,j  ) &
          + M(KS,5,i,j) * Z(KS  ,i+1,j  ) &
          + M(KS,6,i,j) * Z(KS  ,i  ,j-1) &
          + M(KS,7,i,j) * Z(KS  ,i  ,j+1) &
          + M(KS,12,i,j)* Z(KS+1,i-1,j  ) &
          + M(KS,13,i,j)* Z(KS+1,i+1,j  ) &
          + M(KS,14,i,j)* Z(KS+1,i  ,j-1) &
          + M(KS,15,i,j)* Z(KS+1,i  ,j+1) ) ) /diag(KS,i,j)

   enddo
   enddo

   return
  end subroutine back_sub_ILU

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine gs_ILU( KA, KS, KE, &
                     IA, IS, IE, &
                     JA, JS, JE, &
                     Z,  M,  V , diag) !M*Z=V -> Z=M^{-1} V

   implicit none
   integer,  intent(in)  :: KA, KS, KE
   integer,  intent(in)  :: IA, IS, IE
   integer,  intent(in)  :: JA, JS, JE
   real(RP), intent(in)  :: V(KA,IA,JA)
   real(RP), intent(in)  :: M(KA,15,IA,JA)
   real(RP), intent(in)  :: diag(KA,IA,JA)
   real(RP), intent(out) :: Z(KA,IA,JA)

   integer :: k, i, j

   Z(:,:,:)=0.0_RP

   do j = JS, JE
   do i = IS, IE

         Z(KS,i,j) = (V(KS,i,j) &
         -( M(KS,3,i,j) * Z(KS+1,i  ,j  ) &
          + M(KS,4,i,j) * Z(KS  ,i-1,j  ) &
          + M(KS,5,i,j) * Z(KS  ,i+1,j  ) &
          + M(KS,6,i,j) * Z(KS  ,i  ,j-1) &
          + M(KS,7,i,j) * Z(KS  ,i  ,j+1) &
          + M(KS,12,i,j)* Z(KS+1,i-1,j  ) &
          + M(KS,13,i,j)* Z(KS+1,i+1,j  ) &
          + M(KS,14,i,j)* Z(KS+1,i  ,j-1) &
          + M(KS,15,i,j)* Z(KS+1,i  ,j+1) ) ) /diag(KS,i,j)

         do k = KS+1, KE-1
            Z(k,i,j) = (V(k,i,j) &
            -( M(k,2,i,j) * Z(k-1,i  ,j  ) &
             + M(k,3,i,j) * Z(k+1,i  ,j  ) &
             + M(k,4,i,j) * Z(k  ,i-1,j  ) &
             + M(k,5,i,j) * Z(k  ,i+1,j  ) &
             + M(k,6,i,j) * Z(k  ,i  ,j-1) &
             + M(k,7,i,j) * Z(k  ,i  ,j+1) &
             + M(k,8,i,j) * Z(k-1,i-1,j  ) &
             + M(k,9,i,j) * Z(k-1,i+1,j  ) &
             + M(k,10,i,j)* Z(k-1,i  ,j-1) &
             + M(k,11,i,j)* Z(k-1,i  ,j+1) &
             + M(k,12,i,j)* Z(k+1,i-1,j  ) &
             + M(k,13,i,j)* Z(k+1,i+1,j  ) &
             + M(k,14,i,j)* Z(k+1,i  ,j-1) &
             + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/diag(k,i,j)
         enddo

         Z(KE,i,j) = ( V(KE,i,j) &
         -( M(KE,2,i,j) * Z(KE-1,i  ,j  ) &
          + M(KE,4,i,j) * Z(KE  ,i-1,j  ) &
          + M(KE,5,i,j) * Z(KE  ,i+1,j  ) &
          + M(KE,6,i,j) * Z(KE  ,i  ,j-1) &
          + M(KE,7,i,j) * Z(KE  ,i  ,j+1) &
          + M(KE,8,i,j) * Z(KE-1,i-1,j  ) &
          + M(KE,9,i,j) * Z(KE-1,i+1,j  ) &
          + M(KE,10,i,j)* Z(KE-1,i  ,j-1) &
          + M(KE,11,i,j)* Z(KE-1,i  ,j+1) ) )/diag(KE,i,j)

   enddo
   enddo

   return
  end subroutine gs_ILU

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine sgs( KA, KS, KE, &
                  IA, IS, IE, &
                  JA, JS, JE, &
                  Z,  M,  V)
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: M(KA,15,IA,JA)
    real(RP), intent(in)  :: V(KA,IA,JA)
    real(RP), intent(out) :: Z(KA,IA,JA)
    real(RP):: Y(KA,IA,JA)
    integer :: k,i,j

    !$acc data copyin(M,V) copyout(Z) create(Y)

    call gs( KA, KS, KE, &
             IA, IS, IE, &
             JA, JS, JE, &
             Y,  M,  V )

    !$omp parallel do private(i, j, k)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
            Y(k,i,j)=M(k,1,i,j)*Y(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    call back_sub(KA, KS, KE, &
                  IA, IS, IE, &
                  JA, JS, JE, &
                  Z,  M,  Y )

    !$acc end data

    return
  end subroutine sgs

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine back_sub ( KA, KS, KE, &
                        IA, IS, IE, &
                        JA, JS, JE, &
                        Z,  M,  V   )
   implicit none
   integer,  intent(in)  :: KA, KS, KE
   integer,  intent(in)  :: IA, IS, IE
   integer,  intent(in)  :: JA, JS, JE
   real(RP), intent(in)  :: V(KA,IA,JA)
   real(RP), intent(in)  :: M(KA,15,IA,JA)
   real(RP), intent(out) :: Z(KA,IA,JA)
   integer :: k, i, j

   !$acc data copyin(V,M) copyout(Z)

!$acc kernels
z(:,:,:)=1d30
!$acc end kernels

   !$omp parallel
   !$omp do
   !$acc parallel vector_length(32)
   do i = IS, IE
   do k = KS, KE
      Z(k,i,JS-1) = 0.0_RP
      Z(k,i,JE+1) = 0.0_RP
   end do
   end do
   !$acc end parallel
   !$omp do
   !$acc parallel vector_length(32)
   do j = JS, JE
   do k = KS, KE
      Z(k,IS-1,j) = 0.0_RP
      Z(k,IE+1,j) = 0.0_RP
   end do
   end do
   !$acc end parallel
   !$omp end parallel


#if COLORING == 2
   !$omp parallel do collapse(2)
!   !$acc parallel vector_length(32)
   !$acc parallel vector_length(1)
   !$acc loop collapse(2)
   do j = JE, JS, -1
   do i = IE, IS, -1
      if ( mod((IE-i)+(JE-j),2)==0 ) then

      Z(KE,i,j) = ( V(KE,i,j)             )/M(KE,1,i,j)

      do k = KE-1, KS, -1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,3,i,j) * Z(k+1,i  ,j  ) ) )/M(k,1,i,j)
      enddo

      end if

#elif COLORING == 1
   !$omp parallel do
   !$acc parallel vector_length(1)
   do j = JE, JS, -2
   !$acc loop seq
   do i = IE, IS, -1

      Z(KE,i,j) = ( V(KE,i,j) &
      -( M(KE,5,i,j) * Z(KE  ,i+1,j  ) ) )/M(KE,1,i,j)

      do k = KE-1, KS, -1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,3,i,j) * Z(k+1,i  ,j  ) &
          + M(k,5,i,j) * Z(k  ,i+1,j  ) &
          + M(k,13,i,j)* Z(k+1,i+1,j  ) ) )/M(k,1,i,j)
      enddo
#else
   !$acc parallel
   !$acc loop seq
   do j = JE, JS, -1
   !$acc loop seq
   do i = IE, IS, -1

      Z(KE,i,j) = ( V(KE,i,j) &
      -( M(KE,5,i,j) * Z(KE  ,i+1,j  ) &
       + M(KE,7,i,j) * Z(KE  ,i  ,j+1) ) )/M(KE,1,i,j)

      do k = KE-1, KS, -1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,3,i,j) * Z(k+1,i  ,j  ) &
          + M(k,5,i,j) * Z(k  ,i+1,j  ) &
          + M(k,7,i,j) * Z(k  ,i  ,j+1) &
          + M(k,13,i,j)* Z(k+1,i+1,j  ) &
          + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/M(k,1,i,j)
      enddo
#endif
   end do
   end do
   !$acc end parallel

#if COLORING > 0
#if COLORING == 2
   !$omp parallel do collapse(2)
!   !$acc parallel vector_length(32)
   !$acc parallel vector_length(1)
   !$acc loop collapse(2)
   do j = JE, JS, -1
   do i = IE, IS, -1
      if ( mod((IE-i)+(JE-j),2)==1 ) then

      Z(KE,i,j) = ( V(KE,i,j) &
      -( M(KE,4,i,j) * Z(KE  ,i-1,j  ) &
       + M(KE,5,i,j) * Z(KE  ,i+1,j  ) &
       + M(KE,6,i,j) * Z(KE  ,i  ,j-1) &
       + M(KE,7,i,j) * Z(KE  ,i  ,j+1) &
       + M(KE,8,i,j) * Z(KE-1,i-1,j  ) &
       + M(KE,9,i,j) * Z(KE-1,i+1,j  ) &
       + M(KE,10,i,j)* Z(KE-1,i  ,j-1) &
       + M(KE,11,i,j)* Z(KE-1,i  ,j+1) ) )/M(KE,1,i,j)

      do k = KE-1, KS+1,-1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,3,i,j) * Z(k+1,i  ,j  ) &
          + M(k,4,i,j) * Z(k  ,i-1,j  ) &
          + M(k,5,i,j) * Z(k  ,i+1,j  ) &
          + M(k,6,i,j) * Z(k  ,i  ,j-1) &
          + M(k,7,i,j) * Z(k  ,i  ,j+1) &
          + M(k,8,i,j) * Z(k-1,i-1,j  ) &
          + M(k,9,i,j) * Z(k-1,i+1,j  ) &
          + M(k,10,i,j)* Z(k-1,i  ,j-1) &
          + M(k,11,i,j)* Z(k-1,i  ,j+1) &
          + M(k,12,i,j)* Z(k+1,i-1,j  ) &
          + M(k,13,i,j)* Z(k+1,i+1,j  ) &
          + M(k,14,i,j)* Z(k+1,i  ,j-1) &
          + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/M(k,1,i,j)
      enddo

      Z(KS,i,j) = (V(KS,i,j) &
      -( M(KS,3,i,j) * Z(KS+1,i  ,j  ) &
       + M(KS,4,i,j) * Z(KS  ,i-1,j  ) &
       + M(KS,5,i,j) * Z(KS  ,i+1,j  ) &
       + M(KS,6,i,j) * Z(KS  ,i  ,j-1) &
       + M(KS,7,i,j) * Z(KS  ,i  ,j+1) &
       + M(KS,12,i,j)* Z(KS+1,i-1,j  ) &
       + M(KS,13,i,j)* Z(KS+1,i+1,j  ) &
       + M(KS,14,i,j)* Z(KS+1,i  ,j-1) &
       + M(KS,15,i,j)* Z(KS+1,i  ,j+1) ) ) /M(KS,1,i,j)

      end if
#elif COLORING == 1
   !$omp parallel do
   !$acc parallel vector_length(1)
   do j = JE-1, JS, -2
   !$acc loop seq
   do i = IE, IS, -1
      Z(KE,i,j) = ( V(KE,i,j) &
      -( M(KE,5,i,j) * Z(KE  ,i+1,j  ) &
       + M(KE,6,i,j) * Z(KE  ,i  ,j-1) &
       + M(KE,7,i,j) * Z(KE  ,i  ,j+1) &
       + M(KE,9,i,j) * Z(KE-1,i+1,j  ) &
       + M(KE,10,i,j)* Z(KE-1,i  ,j-1) &
       + M(KE,11,i,j)* Z(KE-1,i  ,j+1) ) )/M(KE,1,i,j)

      do k = KE-1, KS+1,-1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,3,i,j) * Z(k+1,i  ,j  ) &
          + M(k,5,i,j) * Z(k  ,i+1,j  ) &
          + M(k,6,i,j) * Z(k  ,i  ,j-1) &
          + M(k,7,i,j) * Z(k  ,i  ,j+1) &
          + M(k,9,i,j) * Z(k-1,i+1,j  ) &
          + M(k,10,i,j)* Z(k-1,i  ,j-1) &
          + M(k,11,i,j)* Z(k-1,i  ,j+1) &
          + M(k,13,i,j)* Z(k+1,i+1,j  ) &
          + M(k,14,i,j)* Z(k+1,i  ,j-1) &
          + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/M(k,1,i,j)
      enddo

      Z(KS,i,j) = (V(KS,i,j) &
      -( M(KS,3,i,j) * Z(KS+1,i  ,j  ) &
       + M(KS,5,i,j) * Z(KS  ,i+1,j  ) &
       + M(KS,6,i,j) * Z(KS  ,i  ,j-1) &
       + M(KS,7,i,j) * Z(KS  ,i  ,j+1) &
       + M(KS,13,i,j)* Z(KS+1,i+1,j  ) &
       + M(KS,14,i,j)* Z(KS+1,i  ,j-1) &
       + M(KS,15,i,j)* Z(KS+1,i  ,j+1) ) ) /M(KS,1,i,j)
#endif
   end do
   end do
   !$acc end parallel
#endif

   !$acc end data

   return
  end subroutine back_sub

  !----- N. Yamashita and T. Iwashita of Hokkaido Univ. created (2021/2/22)-----
  subroutine gs ( KA, KS, KE, &
                  IA, IS, IE, &
                  JA, JS, JE, &
                  Z,  M,  V   )
   implicit none
   integer,  intent(in)  :: KA, KS, KE
   integer,  intent(in)  :: IA, IS, IE
   integer,  intent(in)  :: JA, JS, JE
   real(RP), intent(in)  :: V(KA,IA,JA)
   real(RP), intent(in)  :: M(KA,15,IA,JA)
   real(RP), intent(out) :: Z(KA,IA,JA)

   integer :: k, i, j, n

   !$acc data copyin(V,M) copyout(Z)

   !$omp parallel
   !$omp do
   !$acc parallel vector_length(32)
   do i = IS, IE
   do k = KS, KE
      Z(k,i,JS-1) = 0.0_RP
      Z(k,i,JE+1) = 0.0_RP
   end do
   end do
   !$acc end parallel
   !$omp do
   !$acc parallel vector_length(32)
   do j = JS, JE
   do k = KS, KE
      Z(k,IS-1,j) = 0.0_RP
      Z(k,IE+1,j) = 0.0_RP
   end do
   end do
   !$acc end parallel
   !$omp end parallel


#if COLORING == 2
   !$omp parallel do collapse(2)
!   !$acc parallel vector_length(32)
   !$acc parallel vector_length(1)
   !$acc loop collapse(2)
   do j = JS, JE
   do i = IS, IE
      if ( mod((i-IS)+(j-JS),2)==0 ) then

      Z(KS,i,j) = (V(KS,i,j)  ) /M(KS,1,i,j)

      do k = KS+1, KE
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,2,i,j) * Z(k-1,i  ,j  ) ) )/M(k,1,i,j)
      enddo

      end if

#elif COLORING == 1
   !$omp parallel do
   !$acc parallel vector_length(1)
   do j = JS, JE, 2
   !$acc loop seq
   do i = IS, IE

      Z(KS,i,j) = (V(KS,i,j) &
      -( M(KS,4,i,j) * Z(KS  ,i-1,j  ) ) ) /M(KS,1,i,j)

      do k = KS+1, KE
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,2,i,j) * Z(k-1,i  ,j  ) &
          + M(k,4,i,j) * Z(k  ,i-1,j  ) &
          + M(k,8,i,j) * Z(k-1,i-1,j  ) ) )/M(k,1,i,j)
      enddo
#else
   !$acc parallel
   !$acc loop seq
   do j = JS, JE
   !$acc loop seq
   do i = IS, IE

      Z(KS,i,j) = (V(KS,i,j) &
      -( M(KS,4,i,j) * Z(KS  ,i-1,j  ) &
       + M(KS,6,i,j) * Z(KS  ,i  ,j-1) ) ) /M(KS,1,i,j)

      do k = KS+1, KE
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,2,i,j) * Z(k-1,i  ,j  ) &
          + M(k,4,i,j) * Z(k  ,i-1,j  ) &
          + M(k,6,i,j) * Z(k  ,i  ,j-1) &
          + M(k,8,i,j) * Z(k-1,i-1,j  ) &
          + M(k,10,i,j)* Z(k-1,i  ,j-1) ) )/M(k,1,i,j)
      enddo
#endif
   end do
   end do
   !$acc end parallel

#if COLORING > 0
#if COLORING == 2
   !$omp parallel do collapse(2)
!   !$acc parallel vector_length(32)
   !$acc parallel vector_length(1)
   !$acc loop collapse(2)
   do j = JS, JE
   do i = IS, IE
      if ( mod((i-IS)+(j-JS),2)==1 ) then

      Z(KS,i,j) = (V(KS,i,j) &
      -( M(KS,4,i,j) * Z(KS  ,i-1,j  )  &
       + M(KS,5,i,j) * Z(KS  ,i+1,j  ) &
       + M(KS,6,i,j) * Z(KS  ,i  ,j-1) &
       + M(KS,7,i,j) * Z(KS  ,i  ,j+1) &
       + M(KS,12,i,j)* Z(KS+1,i-1,j  ) &
       + M(KS,13,i,j)* Z(KS+1,i+1,j  ) &
       + M(KS,14,i,j)* Z(KS+1,i  ,j-1) &
       + M(KS,15,i,j)* Z(KS+1,i  ,j+1) ) ) /M(KS,1,i,j)

      do k = KS+1, KE-1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,2,i,j) * Z(k-1,i  ,j  ) &
          + M(k,4,i,j) * Z(k  ,i-1,j  ) &
          + M(k,5,i,j) * Z(k  ,i+1,j  ) &
          + M(k,6,i,j) * Z(k  ,i  ,j-1) &
          + M(k,7,i,j) * Z(k  ,i  ,j+1) &
          + M(k,8,i,j) * Z(k-1,i-1,j  ) &
          + M(k,9,i,j) * Z(k-1,i+1,j  ) &
          + M(k,10,i,j)* Z(k-1,i  ,j-1) &
          + M(k,11,i,j)* Z(k-1,i  ,j+1) &
          + M(k,12,i,j)* Z(k+1,i-1,j  ) &
          + M(k,13,i,j)* Z(k+1,i+1,j  ) &
          + M(k,14,i,j)* Z(k+1,i  ,j-1) &
          + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/M(k,1,i,j)
      enddo

      Z(KE,i,j) = ( V(KE,i,j) &
      -( M(KE,2,i,j) * Z(KE-1,i  ,j  ) &
       + M(KE,4,i,j) * Z(KE  ,i-1,j  ) &
       + M(KE,5,i,j) * Z(KE  ,i+1,j  ) &
       + M(KE,6,i,j) * Z(KE  ,i  ,j-1) &
       + M(KE,7,i,j) * Z(KE  ,i  ,j+1) &
       + M(KE,8,i,j) * Z(KE-1,i-1,j  ) &
       + M(KE,9,i,j) * Z(KE-1,i+1,j  ) &
       + M(KE,10,i,j)* Z(KE-1,i  ,j-1) &
       + M(KE,11,i,j)* Z(KE-1,i  ,j+1) ) )/M(KE,1,i,j)

      end if
#elif COLORING == 1
   !$omp parallel do
   !$acc parallel vector_length(1)
   do j = JS+1, JE, 2
   !$acc loop seq
   do i = IS, IE
      Z(KS,i,j) = (V(KS,i,j) &
      -( M(KS,4,i,j) * Z(KS  ,i-1,j  )  &
       + M(KS,6,i,j) * Z(KS  ,i  ,j-1) &
       + M(KS,7,i,j) * Z(KS  ,i  ,j+1) &
       + M(KS,12,i,j)* Z(KS+1,i-1,j  ) &
       + M(KS,14,i,j)* Z(KS+1,i  ,j-1) &
       + M(KS,15,i,j)* Z(KS+1,i  ,j+1) ) ) /M(KS,1,i,j)

      do k = KS+1, KE-1
         Z(k,i,j) = (V(k,i,j) &
         -( M(k,2,i,j) * Z(k-1,i  ,j  ) &
          + M(k,4,i,j) * Z(k  ,i-1,j  ) &
          + M(k,6,i,j) * Z(k  ,i  ,j-1) &
          + M(k,7,i,j) * Z(k  ,i  ,j+1) &
          + M(k,8,i,j) * Z(k-1,i-1,j  ) &
          + M(k,10,i,j)* Z(k-1,i  ,j-1) &
          + M(k,11,i,j)* Z(k-1,i  ,j+1) &
          + M(k,12,i,j)* Z(k+1,i-1,j  ) &
          + M(k,14,i,j)* Z(k+1,i  ,j-1) &
          + M(k,15,i,j)* Z(k+1,i  ,j+1) ) )/M(k,1,i,j)
      enddo

      Z(KE,i,j) = ( V(KE,i,j) &
      -( M(KE,2,i,j) * Z(KE-1,i  ,j  ) &
       + M(KE,4,i,j) * Z(KE  ,i-1,j  ) &
       + M(KE,6,i,j) * Z(KE  ,i  ,j-1) &
       + M(KE,7,i,j) * Z(KE  ,i  ,j+1) &
       + M(KE,8,i,j) * Z(KE-1,i-1,j  ) &
       + M(KE,10,i,j)* Z(KE-1,i  ,j-1) &
       + M(KE,11,i,j)* Z(KE-1,i  ,j+1) ) )/M(KE,1,i,j)
#endif
   end do
   end do
   !$acc end parallel
#endif

   !$acc end data

   return
  end subroutine gs

  subroutine mul_matrix( KA,KS,KE, &
                         IA,IS,IE, &
                         JA,JS,JE, &
                         V, M, C)
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(out) :: V(KA,IA,JA)
    real(RP), intent(in)  :: M(KA,15,IA,JA)
    real(RP), intent(in)  :: C(KA,IA,JA)

    integer :: k, i, j

    !$omp parallel do private(i,j,k)
    !$acc kernels copyin(M,C) copyout(V)
    do j = JS, JE
    do i = IS, IE
       V(KS,i,j) = M(KS,1,i,j) * C(KS  ,i  ,j  ) &
                 + M(KS,3,i,j) * C(KS+1,i  ,j  ) &
                 + M(KS,4,i,j) * C(KS  ,i-1,j  ) &
                 + M(KS,5,i,j) * C(KS  ,i+1,j  ) &
                 + M(KS,6,i,j) * C(KS  ,i  ,j-1) &
                 + M(KS,7,i,j) * C(KS  ,i  ,j+1) &
                 + M(KS,12,i,j)* C(KS+1,i-1,j  ) &
                 + M(KS,13,i,j)* C(KS+1,i+1,j  ) &
                 + M(KS,14,i,j)* C(KS+1,i  ,j-1) &
                 + M(KS,15,i,j)* C(KS+1,i  ,j+1)
       do k = KS+1, KE-1
          V(k,i,j) = M(k,1,i,j) * C(k  ,i  ,j  ) &
                   + M(k,2,i,j) * C(k-1,i  ,j  ) &
                   + M(k,3,i,j) * C(k+1,i  ,j  ) &
                   + M(k,4,i,j) * C(k  ,i-1,j  ) &
                   + M(k,5,i,j) * C(k  ,i+1,j  ) &
                   + M(k,6,i,j) * C(k  ,i  ,j-1) &
                   + M(k,7,i,j) * C(k  ,i  ,j+1) &
                   + M(k,8,i,j) * C(k-1,i-1,j  ) &
                   + M(k,9,i,j) * C(k-1,i+1,j  ) &
                   + M(k,10,i,j)* C(k-1,i  ,j-1) &
                   + M(k,11,i,j)* C(k-1,i  ,j+1) &
                   + M(k,12,i,j)* C(k+1,i-1,j  ) &
                   + M(k,13,i,j)* C(k+1,i+1,j  ) &
                   + M(k,14,i,j)* C(k+1,i  ,j-1) &
                   + M(k,15,i,j)* C(k+1,i  ,j+1)
       enddo
       V(KE,i,j) = M(KE,1,i,j) * C(KE  ,i  ,j  ) &
                 + M(KE,2,i,j) * C(KE-1,i  ,j  ) &
                 + M(KE,4,i,j) * C(KE  ,i-1,j  ) &
                 + M(KE,5,i,j) * C(KE  ,i+1,j  ) &
                 + M(KE,6,i,j) * C(KE  ,i  ,j-1) &
                 + M(KE,7,i,j) * C(KE  ,i  ,j+1) &
                 + M(KE,8,i,j) * C(KE-1,i-1,j  ) &
                 + M(KE,9,i,j) * C(KE-1,i+1,j  ) &
                 + M(KE,10,i,j)* C(KE-1,i  ,j-1) &
                 + M(KE,11,i,j)* C(KE-1,i  ,j+1)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine mul_matrix

  !-----------------------------------------------------------------------------
  function f2h( k,i,j,p )

    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ
    use scale_atmos_grid_cartesC_metric, only: &
       GSQRT => ATMOS_GRID_CARTESC_METRIC_GSQRT
    use scale_atmos_grid_cartesC_index, only: &
       I_XYZ

    implicit none

    real(RP) :: f2h
    integer, intent(in) :: k, i, j, p

    f2h = (CDZ(k+p-1)*GSQRT(k+p-1,i,j,I_XYZ) &
        / (CDZ(k)*GSQRT(k,i,j,I_XYZ)+CDZ(k+1)*GSQRT(k+1,i,j,I_XYZ)))

  end function
  !-----------------------------------------------------------------------------
  !> calculate charge neutralization MacGorman et al. (2001)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_LT_neutralization_MG2001( &
       KA, KS, KE, & ! [IN]
       IA, IS, IE, & ! [IN]
       JA, JS, JE, & ! [IN]
       KIJMAX,     & ! [IN]
       IMAX,       & ! [IN]
       JMAX,       & ! [IN]
       Efield,     & ! [IN]
       E_pot,      & ! [IN]
       DENS,       & ! [IN]
       QCRG,       & ! [IN]
       QHYD,       & ! [IN]
       NUM_end,    & ! [INOUT]
       LT_path,    & ! [INOUT]
       fls_int_p,  & ! [OUT]
       d_QCRG      ) ! [OUT]
    use scale_const, only: &
       EPS       => CONST_EPS, &
       EPSvac    => CONST_EPSvac, &
       EPSair    => CONST_EPSair, &
       LARGE_NUM => CONST_HUGE
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ, &
       CX   => ATMOS_GRID_CARTESC_CX, &
       CY   => ATMOS_GRID_CARTESC_CY, &
       CDZ  => ATMOS_GRID_CARTESC_CDZ, &
       CDX  => ATMOS_GRID_CARTESC_CDX, &
       CDY  => ATMOS_GRID_CARTESC_CDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_prc_cartesC, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_N, &
       PRC_HAS_S, &
       PRC_W,     &
       PRC_N,     &
       PRC_E,     &
       PRC_S,     &
       PRC_NW,    &
       PRC_NE,    &
       PRC_SW,    &
       PRC_SE,    &
       PRC_next,  &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_prc, only: &
       PRC_abort,      &
       PRC_IsMaster,   &
       PRC_nprocs,     &
       PRC_masterrank, &
       PRC_MPIbarrier, &
       PRC_myrank
    use scale_comm_cartesC, only: &
       COMM_datatype, &
       COMM_world, &
       COMM_vars8, &
       COMM_wait, &
       COMM_bcast
    use scale_random, only: &
       RANDOM_uniform
    implicit none

    integer,  intent(in)    :: KA, KS, KE
    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    integer,  intent(in)    :: KIJMAX, IMAX, JMAX
    real(RP), intent(in)    :: Efield   (KA,IA,JA,4)    !-- Electric field [V/m] (4->|E|)
    real(RP), intent(in)    :: E_pot    (KA,IA,JA)      !-- Electric potential [V]
    real(RP), intent(in)    :: QCRG     (KA,IA,JA)      !-- Charge density [nC/m3]
    real(RP), intent(in)    :: DENS     (KA,IA,JA)      !-- Total density [kg/m3]
    real(RP), intent(in)    :: QHYD     (KA,IA,JA)      !-- Total hydrometeor mixing ratio
    real(RP), intent(inout) :: NUM_end  (KA,IA,JA,3)    !-- Number of Lightning stop grid
    real(RP), intent(inout) :: LT_path  (KA,IA,JA)      !-- Number of path
    real(RP), intent(out)   :: fls_int_p(KA,IA,JA)      !-- Flash initiation point (0->no flash, 1->flash start at the point)
    real(RP), intent(out)   :: d_QCRG   (KA,IA,JA)      !-- Charge density [nC/m3]

    real(RP) :: E_det, E_x, E_y, E_z, E_max
    real(RP) :: dqrho_pls(KA,IA,JA), dqrho_mns(KA,IA,JA)
    real(RP) :: L_path(KA,IA,JA)
    real(RP) :: sumdqrho_pls, sumdqrho_mns
    real(RP) :: Ncrit, dqrho_cor
    real(RP) :: randnum(1)
    real(RP) :: rprod1, rprod2, rbuf, rbuf2, phi_end(2)
    integer  :: Eovid(4,KIJMAX)
    integer  :: Npls, Nmns, Ntot
    integer  :: count1(PRC_nprocs)
    integer  :: countindx(PRC_nprocs+1)
    integer  :: own_prc_total                         !--- number of grid whose |E| > Eint-dEint for each process
    integer  :: iprod(2), ibuf(6), status(MPI_STATUS_SIZE)
    integer  :: init_point, rank_initpoint, grid_initpoint
    integer  :: current_prc    !--- 1->lightning flash path is include, 0-> nopath in the node
    integer  :: comm_flag, comm_target, stop_flag, corr_flag(2)
    integer  :: end_grid(4), wild_grid(6)
    real(RP) :: pm
    integer  :: k, i, j, ipp, iq, direction, ierr
    integer  :: k1, i1, j1, k2, i2, j2, sign_flash
    integer  :: wild_flag, count_path
    integer  :: flg_num_end(KA,IA,JA,3)
    real(RP) :: Eint_hgt(KA,IA,JA)

    call PROF_rapstart('LT_neut_MG2001', 1)

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       fls_int_p(k,i,j) = 0.0_RP
    end do
    end do
    end do

    do count_path = 1, NUTR_ITMAX


       !--- search grid whose Efield is over threshold of flash inititation
       own_prc_total = 0
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Eint_hgt(k,i,j) = min( Eint*DENS(k,i,j)/rho0,Eint )*flg_eint_hgt &
                          + ( Eint-delEint )*( 1.0_RP-flg_eint_hgt )
          if( flg_eint_hgt == 1.0_RP .and. Eint_hgt(k,i,j) < Estop ) then
             Eint_hgt(k,i,j) = LARGE_NUM !--- ignore when Eint < Estop
          endif
          E_det = Efield(k,i,j,I_lt_abs) - Eint_hgt(k,i,j)
          if( E_det > 0.0_RP ) then
             own_prc_total = own_prc_total + 1
             Eovid(1,own_prc_total) = i
             Eovid(2,own_prc_total) = j
             Eovid(3,own_prc_total) = k
             Eovid(4,own_prc_total) = PRC_myrank
          endif
       enddo
       enddo
       enddo

       call MPI_AllGather( own_prc_total, 1, MPI_integer, &
                           count1, 1, MPI_integer, &
                           COMM_world, ierr )
       countindx(1) = 0

       !**** countindx(PROC_nprocs+1) -> total number of grid with |E|>E_threthold
       do ipp = 1, PRC_nprocs
          countindx(ipp+1) = countindx(ipp) + count1(ipp)
       enddo

       !---- select initial point of flash by random select
       !**** rank_initpoint -> rank number including initpoint
       !**** grid_initpoint -> grid number of init point in rank (rank_initpoint)
       if( PRC_IsMaster ) then
          call RANDOM_uniform(randnum)
          init_point = int( randnum(1) * countindx(PRC_nprocs+1) ) + 1
          do ipp = 1, PRC_nprocs
             if ( countindx(ipp+1) >= init_point ) then
                rank_initpoint = ipp-1
                exit
             end if
          end do
          grid_initpoint = init_point - countindx(rank_initpoint+1)
          ibuf(1) = rank_initpoint
          ibuf(2) = grid_initpoint
       endif

       call COMM_bcast( 2, ibuf )

       rank_initpoint = ibuf(1)
       grid_initpoint = ibuf(2)


       !--- Propagate lightning

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          L_path(k,i,j) = 0.0_RP   !--- +-1 -> path already passed, +-2 -> path calculate current step
          flg_num_end(k,i,j,:) = 0
       end do
       end do
       end do
       wild_flag = 0
       end_grid(:) = 0
       comm_target = -999
       do iq = 1, 2       !--- loop for direction (1-> parallel, 2-> anti-parallel)
          if( iq == 1 ) then
             pm = 1.0_RP
          elseif( iq == 2 ) then
             pm = - 1.0_RP
          endif
          stop_flag = 0
          current_prc = rank_initpoint

          !---- determine initiation point
          if( rank_initpoint == PRC_myrank ) then
             i = Eovid(1,grid_initpoint)
             j = Eovid(2,grid_initpoint)
             k = Eovid(3,grid_initpoint)
             L_path(k,i,j) = pm
             if( iq == 1 ) fls_int_p(k,i,j) = fls_int_p(k,i,j) + 1.0_RP
             sign_flash = 2 + int( sign( 1.0_RP, QCRG(k,i,j) ) )
          else
             i = 0
             j = 0
             k = 0
          endif

          loop_path : do ipp = 1, KIJMAXG  !--- loop for path
             comm_flag = 0
             !---- calculate path of lightning
             if( current_prc == PRC_myrank ) then

                !--- determine direction of path
                E_x = abs( Efield(k,i,j,I_lt_x) )/CDX(i)
                E_y = abs( Efield(k,i,j,I_lt_y) )/CDY(j)
                E_z = abs( Efield(k,i,j,I_lt_z) )/CDZ(k)
                E_det = max( E_x,E_y,E_z )

                i1 = i
                j1 = j
                k1 = k
                if( E_det == E_x ) then
                   direction = int( sign( 1.0_RP,Efield(k,i,j,I_lt_x) ) * pm )
                   i1 = i+direction
                elseif( E_det == E_y ) then
                   direction = int( sign( 1.0_RP,Efield(k,i,j,I_lt_y) ) * pm )
                   j1 = j+direction
                elseif( E_det == E_z ) then
                   direction = int( sign( 1.0_RP,Efield(k,i,j,I_lt_z) ) * pm )
                   k1 = k+direction
                endif


                if( Efield(k1,i1,j1,I_lt_abs) <= Estop ) then
                   !--- stop if |E| < Estop
                   phi_end(iq) = QCRG(k,i,j)
                   wild_flag = 1
                   NUM_end(k,i,j,sign_flash) = NUM_end(k,i,j,sign_flash) + 1.0_RP
                   flg_num_end(k,i,j,sign_flash) = 1
                   L_path(k,i,j) = pm
                   stop_flag = 1
                   end_grid(1) = i
                   end_grid(2) = j
                   end_grid(3) = k
                   end_grid(4) = PRC_myrank
                   if( QHYD(k,i,j) < EPS ) then
                      corr_flag(iq) = 0
                   else
                      corr_flag(iq) = 1
                   endif
                elseif( Efield(k1,i1,j1,I_lt_abs) > Estop ) then
                   !--- propagate lightning path
                   L_path(k, i ,j ) = pm
                   L_path(k1,i1,j1) = pm
                endif

                !--- check whether lightning path reach top or bottom layer
                if( k1 == KE ) then
                   !--- reach at top
                   NUM_end(k,i,j,sign_flash) = NUM_end(k,i,j,sign_flash) + 1.0_RP
                   flg_num_end(k,i,j,sign_flash) = 1
                   L_path(k,i,j) = pm
                   stop_flag = 1
                   end_grid(1) = i
                   end_grid(2) = j
                   end_grid(3) = k
                   end_grid(4) = PRC_myrank
                   if( QHYD(k,i,j) < EPS ) then
                      corr_flag(iq) = 0
                   else
                      corr_flag(iq) = 1
                   endif
                elseif( CZ(k1) < zcg ) then
                   !--- reach at groud
                   NUM_end(k,i,j,2) = NUM_end(k,i,j,2) + 1.0_RP
                   flg_num_end(k,i,j,2) = 1
                   L_path(k,i,j) = pm
                   stop_flag = 1
                   end_grid(1) = i
                   end_grid(2) = j
                   end_grid(3) = k
                   end_grid(4) = PRC_myrank
                   if( QHYD(k,i,j) < EPS ) then
                      corr_flag(iq) = 0
                   else
                      corr_flag(iq) = 1
                   endif
                endif

                if( stop_flag /= 1 ) then
                   !--- check whether lightning path reachs boundary in i-direction
                   if( i1 == IS-1 .and. .not. PRC_HAS_W ) then
                      !--- reach west boundary of global domain(stop)
                      NUM_end(k,i,j,sign_flash) = NUM_end(k,i,j,sign_flash) + 1.0_RP
                      flg_num_end(k,i,j,sign_flash) = 1
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      stop_flag = 1
                      end_grid(1) = i
                      end_grid(2) = j
                      end_grid(3) = k
                      end_grid(4) = PRC_myrank
                      if( QHYD(k,i,j) < EPS ) then
                         corr_flag(iq) = 0
                      else
                         corr_flag(iq) = 1
                      endif
                   elseif( i1 == IS-1 .and. PRC_HAS_W ) then
                      !--- reach west boundary of local domain(propagate)
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      comm_flag = 1
                      comm_target = PRC_next(PRC_W)
                      i1 = IE
                   elseif( i1 == IE+1 .and. .not. PRC_HAS_E ) then
                      !--- reach east boundary of global domain(stop)
                      NUM_end(k,i,j,sign_flash) = NUM_end(k,i,j,sign_flash) + 1.0_RP
                      flg_num_end(k,i,j,sign_flash) = 1
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      stop_flag = 1
                      end_grid(1) = i
                      end_grid(2) = j
                      end_grid(3) = k
                      end_grid(4) = PRC_myrank
                      if( QHYD(k,i,j) < EPS ) then
                         corr_flag(iq) = 0
                      else
                         corr_flag(iq) = 1
                      endif
                   elseif( i1 == IE+1 .and. PRC_HAS_E ) then
                      !--- reach east boundary of local domain(propagate)
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      comm_flag = 1
                      comm_target = PRC_next(PRC_E)
                      i1 = IS
                   endif

                   !--- check whether lightning path reachs boundary in i-direction
                   if( j1 == JS-1 .and. .not. PRC_HAS_S ) then
                      !--- reach south boundary of global domain(stop)
                      NUM_end(k,i,j,sign_flash) = NUM_end(k,i,j,sign_flash) + 1.0_RP
                      flg_num_end(k,i,j,sign_flash) = 1
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      stop_flag = 1
                      end_grid(1) = i
                      end_grid(2) = j
                      end_grid(3) = k
                      end_grid(4) = PRC_myrank
                      if( QHYD(k,i,j) < EPS ) then
                         corr_flag(iq) = 0
                      else
                         corr_flag(iq) = 1
                      endif
                   elseif( j1 == JS-1 .and. PRC_HAS_S ) then
                      !--- reach south boundary of local domain(propagate)
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      comm_flag = 1
                      comm_target = PRC_next(PRC_S)
                      j1 = JE
                   elseif( j1 == JE+1 .and. .not. PRC_HAS_N ) then
                      !--- reach north boundary of global domain(stop)
                      NUM_end(k,i,j,sign_flash) = NUM_end(k,i,j,sign_flash) + 1.0_RP
                      flg_num_end(k,i,j,sign_flash) = 1
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      comm_flag = 1
                      stop_flag = 1
                      end_grid(1) = i
                      end_grid(2) = j
                      end_grid(3) = k
                      end_grid(4) = PRC_myrank
                      if( QHYD(k,i,j) < EPS ) then
                         corr_flag(iq) = 0
                      else
                         corr_flag(iq) = 1
                      endif
                   elseif( j1 == JE+1 .and. PRC_HAS_N ) then
                      !--- reach north boundary of local domain(propagate)
                      L_path(k,i,j) = pm
                      L_path(k1,i1,j1) = pm
                      comm_flag = 1
                      comm_target = PRC_next(PRC_N)
                      j1 = JS
                   endif
                endif

             endif

             !---- determine stop or not
             call MPI_bcast(stop_flag, 1, MPI_integer, current_prc, COMM_world, ierr)

             if( stop_flag == 1 ) then
                !---- send flag wildfire
                call MPI_bcast(wild_flag, 1, MPI_integer, current_prc, COMM_world, ierr)

                if ( PRC_myrank == current_prc ) then
                   ibuf(1:4) = end_grid(:)
                   ibuf(5:6) = corr_flag(:)
                end if
                call MPI_bcast(ibuf, 6, MPI_integer, current_prc, COMM_world, ierr)
                end_grid(:) = ibuf(1:4)
                corr_flag(:) = ibuf(5:6)

                stop_flag = 0
                !---- If lightning path reaches end_grid stop
                exit loop_path
             endif

             !---- determine wether process change occurs or not
             call MPI_bcast(comm_flag, 1, MPI_integer, current_prc, COMM_world, ierr)

             !--- process change occurs
             if( comm_flag == 1 ) then
                call MPI_bcast(comm_target, 1, MPI_integer, current_prc, COMM_world, ierr)

                if ( PRC_myrank == comm_target ) then
                   call MPI_recv(ibuf, 6, MPI_integer, current_prc, 1, COMM_world, status, ierr)
                   k1 = ibuf(1)
                   i1 = ibuf(2)
                   j1 = ibuf(3)
                   corr_flag(:) = ibuf(4:5)
                   sign_flash = ibuf(6)
                end if

                if ( PRC_myrank == current_prc ) then
                   ibuf(1) = k1
                   ibuf(2) = i1
                   ibuf(3) = j1
                   ibuf(4:5) = corr_flag(:)
                   ibuf(6) = sign_flash
                   call MPI_send(ibuf, 6, MPI_integer, comm_target, 1, COMM_world, ierr)
                end if

                !--- change current proc
                current_prc = comm_target
             endif

             if( current_prc == PRC_myrank ) then
                k = k1
                j = j1
                i = i1
             endif

          enddo loop_path

!          call PRC_MPIbarrier

          !---- Wildfire method
          if( wild_flag == 1 .and. end_grid(4) == PRC_myrank ) then
             i1 = end_grid(1)
             j1 = end_grid(2)
             k1 = end_grid(3)
             i2 = end_grid(1)+1
             j2 = end_grid(2)
             k2 = end_grid(3)
             if( abs( QCRG(k2,i2,j2) ) > qrho_chan/DENS(k2,i2,j2) .and. &
                  abs( E_pot(k2,i2,j2) ) > abs( phi_end(iq) ) ) then
                L_path( k2,i2,j2 ) = L_path( k1,i1,j1 )
             endif
             i1 = end_grid(1)
             j1 = end_grid(2)
             k1 = end_grid(3)
             i2 = end_grid(1)-1
             j2 = end_grid(2)
             k2 = end_grid(3)
             if( abs( QCRG(k2,i2,j2) ) > qrho_chan/DENS(k2,i2,j2) .and. &
                  abs( E_pot(k2,i2,j2) ) > abs( phi_end(iq) ) ) then
                L_path( k2,i2,j2 ) = L_path( k1,i1,j1 )
             endif
             i1 = end_grid(1)
             j1 = end_grid(2)
             k1 = end_grid(3)
             i2 = end_grid(1)
             j2 = end_grid(2)+1
             k2 = end_grid(3)
             if( abs( QCRG(k2,i2,j2) ) > qrho_chan/DENS(k2,i2,j2) .and. &
                  abs( E_pot(k2,i2,j2) ) > abs( phi_end(iq) ) ) then
                L_path( k2,i2,j2 ) = L_path( k1,i1,j1 )
             endif
             i1 = end_grid(1)
             j1 = end_grid(2)
             k1 = end_grid(3)
             i2 = end_grid(1)
             j2 = end_grid(2)-1
             k2 = end_grid(3)
             if( abs( QCRG(k2,i2,j2) ) > qrho_chan/DENS(k2,i2,j2) .and. &
                  abs( E_pot(k2,i2,j2) ) > abs( phi_end(iq) ) ) then
                L_path( k2,i2,j2 ) = L_path( k1,i1,j1 )
             endif
             i1 = end_grid(1)
             j1 = end_grid(2)
             k1 = end_grid(3)
             i2 = end_grid(1)
             j2 = end_grid(2)
             k2 = end_grid(3)+1
             if( abs( QCRG(k2,i2,j2) ) > qrho_chan/DENS(k2,i2,j2) .and. &
                  abs( E_pot(k2,i2,j2) ) > abs( phi_end(iq) ) ) then
                L_path( k2,i2,j2 ) = L_path( k1,i1,j1 )
             endif
             i1 = end_grid(1)
             j1 = end_grid(2)
             k1 = end_grid(3)
             i2 = end_grid(1)
             j2 = end_grid(2)
             k2 = end_grid(3)-1
             if( abs( QCRG(k2,i2,j2) ) > qrho_chan/DENS(k2,i2,j2) .and. &
                  abs( E_pot(k2,i2,j2) ) > abs( phi_end(iq) ) ) then
                L_path( k2,i2,j2 ) = L_path( k1,i1,j1 )
             endif
          endif

!          call PRC_MPIbarrier

       enddo  !--- loop for direction

!       call COMM_vars8( L_path,1 )
!       call COMM_wait ( L_path,1 )

       !---- Neutralization
       Npls = 0
       Nmns = 0
       sumdqrho_pls = 0.0_RP
       sumdqrho_mns = 0.0_RP
       !$omp parallel do reduction(+:Npls,Nmns,sumdqrho_pls,sumdqrho_mns)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

          !---- whether the grid is on lightning path or not
          if( L_path(k,i,j) /= 0.0_RP .and. abs( QCRG(k,i,j) ) > qrho_chan ) then
             if ( QCRG(k,i,j) >= 0.0_RP ) then
                Npls = Npls + 1
                dqrho_pls(k,i,j) = fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                sumdqrho_pls = sumdqrho_pls &
                             + fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                dqrho_mns(k,i,j) = 0.0_RP
             else
                Nmns = Nmns + 1
                dqrho_mns(k,i,j) = fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                sumdqrho_mns = sumdqrho_mns &
                             + fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                dqrho_pls(k,i,j) = 0.0_RP
             endif
          else
             dqrho_pls(k,i,j) = 0.0_RP
             dqrho_mns(k,i,j) = 0.0_RP
          endif

       enddo
       enddo
       enddo
       iprod(1) =  Npls
       iprod(2) =  Nmns
       call MPI_AllReduce(iprod, ibuf, 2, MPI_integer, MPI_SUM, COMM_world, ierr)
       Npls = ibuf(1)
       Nmns = ibuf(2)

       Ntot = Npls + Nmns

       if ( Ntot > 0 ) exit

    end do

    if ( count_path > NUTR_ITMAX ) then
       LOG_INFO("ATMOS_PHY_LT_Efield",*) "Reach limit iteration for searching flash path", count_path, Npls, Nmns, current_prc
    endif

    if( corr_flag(1) == 1 .and. corr_flag(2) == 1 ) then
      dqrho_cor = sumdqrho_pls - sumdqrho_mns
      rprod2 = dqrho_cor
      call MPI_AllReduce(rprod2, rbuf, 1, COMM_datatype, MPI_SUM, COMM_world, ierr)
      dqrho_cor = rbuf

      dqrho_cor = dqrho_cor / Ntot    !---- Eq.(2) of MacGorman et al. (2001)
    else
      dqrho_cor = 0.0_RP              !---- No correlation if flash to ground and out of cloud
    endif

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       d_QCRG(k,i,j) = 0.0_RP
    end do
    end do
    end do
    rprod1 = 0.0_RP
    rprod2 = 0.0_RP
    !--- pseud-2D experiment for x-direction
    if( PRC_NUM_X == 1 .and. IMAX == 2 ) then

       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
         if( L_path(k,IS,j) /= 0.0_RP ) then
           L_path(k,IE,j) =   L_path(k,IS,j)   !--- copy IS value to IE (which is regarded as same grid in psued-2D)
           if( dqrho_pls(k,IS,j) /= 0.0_RP ) then
             d_QCRG(k,IS,j) = - dqrho_pls(k,IS,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,IE,j) = d_QCRG(k,IS,j)   !--- copy IS value to IE (which is regarded as same grid in psued-2D)
           elseif( dqrho_mns(k,IS,j) /= 0.0_RP ) then
             d_QCRG(k,IS,j) =   dqrho_mns(k,IS,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,IE,j) = d_QCRG(k,IS,j)   !--- copy IS value to IE (which is regarded as same grid in psued-2D)
           endif
         elseif( L_path(k,IE,j) /= 0.0_RP ) then
           L_path(k,IS,j) =   L_path(k,IE,j)    !--- copy IE value to IS (which is regarded as same grid in psued-2D)
           if( dqrho_pls(k,IE,j) /= 0.0_RP ) then
             d_QCRG(k,IE,j) = - dqrho_pls(k,IE,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,IS,j) = d_QCRG(k,IE,j)   !--- copy IS value to IS (which is regarded as same grid in psued-2D)
           elseif( dqrho_mns(k,IE,j) /= 0.0_RP ) then
             d_QCRG(k,IE,j) =   dqrho_mns(k,IE,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,IS,j) = d_QCRG(k,IE,j)   !--- copy IS value to IS (which is regarded as same grid in psued-2D)
           endif
         endif

         do iq = 1, 3
           if( flg_num_end(k,IS,j,iq) == 1 ) then
               NUM_end(k,IE,j,iq) = NUM_end(k,IE,j,iq) + 1.0_RP
           endif
           if( flg_num_end(k,IE,j,iq) == 1 ) then
               NUM_end(k,IS,j,iq) = NUM_end(k,IS,j,iq) + 1.0_RP
           endif
         enddo

      enddo
      enddo

    !--- pseud-2D experiment for y-direction
    elseif( PRC_NUM_Y == 1 .and. JMAX == 2 ) then

       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
         if( L_path(k,i,JS) /= 0.0_RP ) then
           L_path(k,i,JE) =   L_path(k,i,JS)   !--- copy JS value to JE (which is regarded as same grid in psued-2D)
           if( dqrho_pls(k,i,JS) /= 0.0_RP ) then
             d_QCRG(k,i,JS) = - dqrho_pls(k,i,JS) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,i,JE) = d_QCRG(k,i,JS)  !--- copy JS value to JE (which is regarded as same grid in psued-2D)
           elseif( dqrho_mns(k,i,JS) /= 0.0_RP ) then
             d_QCRG(k,i,JS) =   dqrho_mns(k,i,JS) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,i,JE) = d_QCRG(k,i,JS)  !--- copy JS value to JE (which is regarded as same grid in psued-2D)
           endif
         elseif( L_path(k,i,JE) /= 0.0_RP ) then
           L_path(k,i,JS) =   L_path(k,i,JE)   !--- copy JE value to JS (which is regarded as same grid in psued-2D)
           if( dqrho_pls(k,i,JE) /= 0.0_RP ) then
             d_QCRG(k,i,JE) = - dqrho_pls(k,i,JE) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,i,JS) = d_QCRG(k,i,JE)  !--- copy JS value to JE (which is regarded as same grid in psued-2D)
           elseif( dqrho_mns(k,i,JE) /= 0.0_RP ) then
             d_QCRG(k,i,JE) =   dqrho_mns(k,i,JE) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
             d_QCRG(k,i,JS) = d_QCRG(k,i,JE)  !--- copy JS value to JE (which is regarded as same grid in psued-2D)
           endif
         endif

         do iq = 1, 3
           if( flg_num_end(k,i,JS,iq) == 1 ) then
               NUM_end(k,i,JE,iq) = NUM_end(k,i,JE,iq) + 1.0_RP
           endif
           if( flg_num_end(k,i,JE,iq) == 1 ) then
               NUM_end(k,i,JS,iq) = NUM_end(k,i,JS,iq) + 1.0_RP
           endif
         enddo

      enddo
      enddo

    else

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if( L_path(k,i,j) /= 0.0_RP .and. dqrho_pls(k,i,j) /= 0.0_RP ) then
             d_QCRG(k,i,j) = - dqrho_pls(k,i,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
          elseif( L_path(k,i,j) /= 0.0_RP .and. dqrho_mns(k,i,j) /= 0.0_RP ) then
             d_QCRG(k,i,j) =   dqrho_mns(k,i,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
          endif
       enddo
       enddo
       enddo

    endif

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       LT_path(k,i,j) = LT_path(k,i,j) + L_path(k,i,j)
    enddo
    enddo
    enddo

    call PROF_rapend('LT_neut_MG2001', 1)

    return

  end subroutine ATMOS_PHY_LT_neutralization_MG2001
  !-----------------------------------------------------------------------------
  !> calculate charge neutralization Fierro et al. (2013)
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_LT_neutralization_F2013( &
       KA, KS, KE, & ! [IN]
       IA, IS, IE, & ! [IN]
       JA, JS, JE, & ! [IN]
       KIJMAX,     & ! [IN]
       Efield,     & ! [IN]
       E_pot,      & ! [IN]
       DENS,       & ! [IN]
       QCRG,       & ! [IN]
       QHYD,       & ! [IN]
       NUM_end,    & ! [INOUT]
       LT_path,    & ! [INOUT]
       fls_int_p,  & ! [OUT]
       d_QCRG,     & ! [OUT]
       B_OUT       ) ! [OUT]
    use scale_const, only: &
       EPS    => CONST_EPS, &
       EPSvac => CONST_EPSvac, &
       EPSair => CONST_EPSair
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ, &
       CX   => ATMOS_GRID_CARTESC_CX, &
       CY   => ATMOS_GRID_CARTESC_CY, &
       CDZ  => ATMOS_GRID_CARTESC_CDZ, &
       CDX  => ATMOS_GRID_CARTESC_CDX, &
       CDY  => ATMOS_GRID_CARTESC_CDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_prc_cartesC, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_N, &
       PRC_HAS_S, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S,    &
       PRC_NW,   &
       PRC_NE,   &
       PRC_SW,   &
       PRC_SE,   &
       PRC_next, &
       PRC_NUM_X,&
       PRC_NUM_Y
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster, &
       PRC_nprocs, &
       PRC_masterrank, &
       PRC_MPIbarrier, &
       PRC_myrank
    use scale_comm_cartesC, only: &
       COMM_datatype, &
       COMM_world
!    use scale_random, only: &
!       RANDOM_reset, &
!       RANDOM_uniform
    implicit none

    integer,  intent(in)    :: KA, KS, KE
    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    integer,  intent(in)    :: KIJMAX
    real(RP), intent(in)    :: Efield   (KA,IA,JA,4)    !-- Electric field [V/m] (4->|E|)
    real(RP), intent(in)    :: E_pot    (KA,IA,JA)      !-- Electric potential [V]
    real(RP), intent(in)    :: QCRG     (KA,IA,JA)      !-- Charge density [nC/m3]
    real(RP), intent(in)    :: DENS     (KA,IA,JA)      !-- Total density [kg/m3]
    real(RP), intent(in)    :: QHYD     (KA,IA,JA)      !-- Total hydrometeor mixing ratio
    real(RP), intent(inout) :: NUM_end  (KA,IA,JA,3)    !-- Number of Lightning stop grid
    real(RP), intent(inout) :: LT_path   (KA,IA,JA)     !-- Number of path
    real(RP), intent(out)   :: fls_int_p(KA,IA,JA)      !-- Flash initiation point (0->no flash, 1->flash start at the point)
    real(RP), intent(out)   :: d_QCRG   (KA,IA,JA)      !-- Charge density [nC/m3]
    real(RP), intent(out)   :: B_OUT    (IA,JA)         !-- B value for output

    real(RP), parameter :: q_thre = 0.1_RP ! threshold of discharge zone (Fierro et al. 2013) [nC/m3]

    real(RP) :: B(IA,JA)   !--- B in Fierro et al. 2013
    integer  :: C(IA,JA)   !--- flg for cylinder (in cylinder -> C=1)
    real(RP) :: Q_d, Fpls, Fmns
    real(RP) :: Spls, Smns
    real(RP) :: Spls_g, Smns_g
    real(RP) :: distance
    real(RP) :: Edif(KS:KE), Edif_max, abs_qcrg_max
    real(RP) :: sw

    !integer, allocatable :: proc_num(:), proc_numg(:)
    real(RP),allocatable :: E_exce_x(:), E_exce_x_g(:)  !--- x point of column in which |E|>Eint is included [m] (_g means global attribute)
    real(RP),allocatable :: E_exce_y(:), E_exce_y_g(:)  !--- y point of column in which |E|>Eint is included [m] (_g means global attribute)
    real(RP) :: exce_grid(2,KIJMAX)
    integer  :: count1(PRC_nprocs)
    integer  :: countindx(PRC_nprocs+1)
    integer  :: num_own                               !--- number of column whose |E| > Eint-dEint for each process
    integer  :: num_total                             !--- total number of column whose |E| > Eint-dEint
    real(RP) :: rbuf1, rbuf2
    integer  :: k, i, j, ipp, iq, ierr

    call PROF_rapstart('LT_neut_F2013', 1)

    !$acc data copyin(Efield,E_pot,QCRG,DENS,QHYD,CX,CY) &
    !$acc      copy(NUM_end,LT_path) &
    !$acc      copyout(fls_int_p,d_QCRG,B_OUT) &
    !$acc      create(B,C,Edif,exce_grid)

    !$omp parallel do
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, IE
    do iq = 1, 3
       NUM_end(k,i,j,iq) = 0.0_RP
    end do
    end do
    end do
    end do
    !$acc end kernels

    !--- search grid whose Efield is over threshold of flash inititation
    num_own = 0
    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JS, JE
    do i = IS, IE
       Edif_max = -1.0_RP
       !$acc loop reduction(max:Edif_max)
       do k = KS, KE
          Edif(k) = Efield(k,i,j,I_lt_abs) - ( Eint-delEint )
          Edif_max = max( Edif_max, Edif(k) )
       enddo
       if ( Edif_max > 0.0_RP ) then
          !$acc atomic
          num_own = num_own + 1
          !$acc end atomic
          exce_grid(1,num_own) = CX(i)
          exce_grid(2,num_own) = CY(j)
          do k = KS, KE
             fls_int_p(k,i,j) = 0.5_RP + sign( 0.5_RP, Edif(k)-EPS )
          enddo
       else
          do k = KS, KE
             fls_int_p(k,i,j) = 0.0_RP
          enddo
       endif
    enddo
    enddo
    !$acc end kernels

    !############## calculation by CPU ################
    !$acc update host(exce_grid)
    !**** proc_num(0~) -> process number of each grid with |E|> E_threthold (local)
!    allocate(proc_num(own_prc_total))
    allocate(E_exce_x(num_own))
    allocate(E_exce_y(num_own))

!    proc_num(:) = PRC_myrank
    E_exce_x(1:num_own) = exce_grid(1,1:num_own)
    E_exce_y(1:num_own) = exce_grid(2,1:num_own)

    call MPI_AllGather( num_own, 1, MPI_integer, &
                        count1(:), 1, MPI_integer, &
                        COMM_world, ierr )

    countindx(1) = 0
    do ipp = 1, PRC_nprocs
      countindx(ipp+1) = countindx(ipp) + count1(ipp)
    enddo

    !---- Create global version of proc_num(proc_numg)
    !**** countindx(PROC_nprocs) -> total number of grid with |E|>E_threthold
    !**** proc_numg(0~) -> process number of each grid with |E|> E_threthold (global)
    num_total = countindx(PRC_nprocs+1)
    allocate(E_exce_x_g(num_total))
    allocate(E_exce_y_g(num_total))

    call MPI_AllGatherv( E_exce_x,   num_own, COMM_datatype, &
                         E_exce_x_g, count1, countindx, COMM_datatype, &
                         COMM_world, ierr )

    call MPI_AllGatherv( E_exce_y,   num_own, COMM_datatype, &
                         E_exce_y_g, count1, countindx, COMM_datatype, &
                         COMM_world, ierr )
    !$acc data copyin(E_exce_x_g,E_exce_y_g)
    !############## calculation by CPU ################

    !$acc kernels
    C(:,:) = 0
    do ipp = 1, num_total
       do j = JS, JE
       do i = IS, IE
          distance = sqrt( ( CX(i)-E_exce_x_g(ipp) )**2 &
                         + ( CY(j)-E_exce_y_g(ipp) )**2 )
          if( distance <= R_neut ) then
             C(i,j) = 1
          endif
       enddo
       enddo
    enddo
    !$acc end kernels

    !$acc end data

    Spls = 0.0_RP
    Smns = 0.0_RP
    !$omp parallel do reduction(+:Spls,Smns) &
    !$omp private(abs_qcrg_max,sw)
    !$acc kernels
    !$acc loop collapse(2) reduction(+:Spls,Smns)
    do j = JS, JE
    do i = IS, IE

       B(i,j) = 0.0_RP
       if ( C(i,j) == 1 ) then
          abs_qcrg_max = 0.0_RP
          !$acc loop reduction(max:abs_qcrg_max)
          do k = KS, KE
             abs_qcrg_max = max( abs_qcrg_max, abs( QCRG(k,i,j) ) )
          end do
          if ( abs_qcrg_max >= q_thre ) then
             B(i,j) = 1.0_RP
             do k = KS, KE
                sw = 0.5_RP + sign( 0.5_RP, QCRG(k,i,j) )
                Spls = Spls + QCRG(k,i,j) * sw
                Smns = Smns - QCRG(k,i,j) * ( 1.0_RP - sw )
             enddo
          end if
       end if

    enddo
    enddo
    !$acc end kernels

    rbuf1 = Spls
    call MPI_AllReduce( rbuf1, rbuf2, 1, COMM_datatype, &
                        MPI_SUM, COMM_world, ierr       )
    Spls_g = rbuf2

    rbuf1 = Smns
    call MPI_AllReduce( rbuf1, rbuf2, 1, COMM_datatype, &
                        MPI_SUM, COMM_world, ierr       )
    Smns_g = rbuf2

    if( max( Spls_g,Smns_g )*0.3_RP < min( Spls_g,Smns_g ) ) then
      Q_d = 0.3_RP * max( Spls_g,Smns_g )
    else
      Q_d = min( Spls_g,Smns_g )
    endif

    if( Spls_g /= 0.0_RP ) then
      Fpls = Q_d / Spls_g
    else
      Fpls = 0.0_RP
    endif

    if( Smns_g /= 0.0_RP ) then
      Fmns = Q_d / Smns_g
    else
      Fmns = 0.0_RP
    endif

    !---- select initial point of flash by random select
    !$omp parallel do
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if( QCRG(k,i,j) > q_thre ) then
          d_QCRG(k,i,j) = - Fpls * (   QCRG(k,i,j) - q_thre ) * B(i,j)
          LT_path(k,i,j) = LT_path(k,i,j) + B(i,j)
       elseif( QCRG(k,i,j) < -q_thre ) then
          d_QCRG(k,i,j) = + Fmns * ( - QCRG(k,i,j) - q_thre ) * B(i,j)
          LT_path(k,i,j) = LT_path(k,i,j) - B(i,j)
       else
          d_QCRG(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       B_OUT(i,j) = B(i,j)
    enddo
    enddo
    !$acc end kernels

    deallocate(E_exce_x)
    deallocate(E_exce_y)
    deallocate(E_exce_x_g)
    deallocate(E_exce_y_g)

    !$acc end data

    call PROF_rapend('LT_neut_F2013', 1)

    return
  end subroutine ATMOS_PHY_LT_neutralization_F2013
  !-----------------------------------------------------------------------------
  !> judge max of |E| exceeds E_init
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_LT_judge_absE( &
       KA, KS, KE, & ! [IN]
       IA, IS, IE, & ! [IN]
       JA, JS, JE, & ! [IN]
       DENS,       & ! [IN]
       Efield,     & ! [IN]
       Emax,       & ! [OUT]
       flg_lt_neut ) ! [OUT]
    use scale_comm_cartesC, only: &
       COMM_datatype, &
       COMM_world
    use scale_const, only: &
       LARGE_NUM => CONST_HUGE
    implicit none

    integer,  intent(in)    :: KA, KS, KE
    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    real(RP), intent(in)    :: DENS     (KA,IA,JA)    !-- Total density [kg/m3]
    real(RP), intent(in)    :: Efield   (KA,IA,JA)    !-- Absolute value of Electric field [V/m]
    real(RP), intent(out)   :: Emax                   !-- Maximum value of Electric field [V/m]
    integer,  intent(out)   :: flg_lt_neut            !-- flg 1-> neutralization, 0->no neutralization

    real(RP) :: Eint_hgt(KA,IA,JA)
    real(RP) :: E_det, rbuf, rprod1
    integer  :: own_prc_total, iprod1, buf
    integer  :: k, i, j, ierr

    !$acc data &
    !$acc copyin(DENS,Efield) &
    !$acc create(Eint_hgt)

    !--- search grid whose Efield is over threshold of flash inititation
    own_prc_total = 0
    if( flg_eint_hgt == 1.0_RP ) then
       !$omp parallel do &
       !$omp reduction(+:own_prc_total) &
       !$omp private(E_det)
       !$acc kernels
       !$acc loop collapse(3) reduction(+:own_prc_total)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          Eint_hgt(k,i,j) = min( Eint*DENS(k,i,j)/rho0,Eint )
          if( Eint_hgt(k,i,j) < Estop ) then
             Eint_hgt(k,i,j) = LARGE_NUM !--- ignore when Eint < Estop
          endif
          E_det = Efield(k,i,j) - Eint_hgt(k,i,j)
          if( E_det > 0.0_RP ) then
             own_prc_total = own_prc_total + 1
          endif
       enddo
       enddo
       enddo
       !$acc end kernels
    else
       !$omp parallel do &
       !$omp reduction(+:own_prc_total) &
       !$omp private(E_det)
       !$acc kernels
       !$acc loop collapse(3) reduction(+:own_prc_total)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          E_det = Efield(k,i,j) - ( Eint-delEint )
          if( E_det > 0.0_RP ) then
             own_prc_total = own_prc_total + 1
          endif
       enddo
       enddo
       enddo
       !$acc end kernels
    endif

    !--- Add number of grids, whose |E| are over threshold of flash initiaion, for all process
    iprod1 = own_prc_total
    call MPI_AllReduce(iprod1, buf, 1, MPI_integer, MPI_SUM, COMM_world, ierr)
    rprod1 = 0.0_RP
    !$acc kernels
    !$acc loop collapse(3) reduction(max:rprod1)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rprod1 = max(rprod1,Efield(k,i,j))
    enddo
    enddo
    enddo
    !$acc end kernels
    call MPI_AllReduce(rprod1, rbuf, 1, COMM_datatype, MPI_MAX, COMM_world, ierr)
    !--- exit when no grid point with |E| over threshold of flash initiation exist
    Emax = rbuf
    if( buf == 0 ) then
      flg_lt_neut = 0
    else
      flg_lt_neut = 1
      if( Emax < Eint .and. Emax >= Eint-delEint .and. flg_eint_hgt == 0.0_RP ) then
        flg_lt_neut = 2
      endif
    endif

    !$acc end data

  end subroutine ATMOS_PHY_LT_judge_absE
  !-----------------------------------------------------------------------------
  !> Select cwc-temp point on LUT
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT( &
       KA, KS, KE, & ! [IN]
       IA, IS, IE, & ! [IN]
       JA, JS, JE, & ! [IN]
       NLIQ,      & ! [IN]
       TEMP,       & ! [IN]
       DENS,       & ! [IN]
       QLIQ,       & ! [IN]
       dqcrg,      & ! [OUT]
       beta_crg    ) ! [OUT]
    use scale_const, only: &
       T00 => CONST_TEM00
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    integer,  intent(in)  :: NLIQ
    real(RP), intent(in)  :: TEMP(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QLIQ(KA,IA,JA,NLIQ)
    real(RP), intent(out) :: dqcrg(KA,IA,JA)
    real(RP), intent(out) :: beta_crg(KA,IA,JA)

    integer  :: i, j, k, pp, qq, iq
    integer  :: grid1, grid2
    real(RP) :: cwc
    real(RP) :: diffx, diffy, tmp

    !$acc data copyin(TEMP,DENS,QLIQ) copyout(dqcrg,beta_crg)

    !$omp parallel do &
    !$omp private(cwc,diffx,diffy,tmp,grid1,grid2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if( TEMP(k,i,j) <= T00 .and. TEMP(k,i,j) >= tcrglimit ) then
          cwc = 0.0_RP
          !$acc loop seq
          do iq = 1, NLIQ
             cwc = cwc + QLIQ(k,i,j,iq) * DENS(k,i,j) * 1.0E+3_RP ![g/m3]
          enddo
          diffx = abs( grid_lut_t(1)-TEMP(k,i,j) )
          grid1 = 1
          !$acc loop seq
          do pp = 2, nxlut_lt
             tmp = abs( grid_lut_t(pp)-TEMP(k,i,j) )
             if ( tmp < diffx ) then
                diffx = tmp
                grid1 = pp
             end if
          enddo
          diffy = abs( grid_lut_l(1)-cwc )
          grid2 = 1
          !$acc loop seq
          do qq = 2, nylut_lt
             tmp = abs( grid_lut_l(qq)-cwc )
             if ( tmp < diffy ) then
                diffy = tmp
                grid2 = qq
             end if
          enddo
          dqcrg(k,i,j) = dq_chrg( grid1, grid2 ) &
                       *( 0.5_RP + sign( 0.5_RP,cwc-1.0E-2_RP ) ) !--- no charge separation when cwc < 0.01 [g/m3]
       else
          dqcrg(k,i,j) = 0.0_RP
       endif
       if( TEMP(k,i,j) >= -30.0_RP+T00 ) then
          beta_crg(k,i,j) = 1.0_RP
       elseif( TEMP(k,i,j) < -30.0_RP+T00 .and. TEMP(k,i,j) >= -43.0_RP+T00 ) then
          beta_crg(k,i,j) = 1.0_RP - ( ( TEMP(k,i,j)-T00+30.0_RP )/13.0_RP )**2
!      elseif( TEMP(k,i,j) < -43.0_RP+T00 ) then
       else
          beta_crg(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    return
 end subroutine ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_lt_sato2019
