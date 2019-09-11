!-------------------------------------------------------------------------------
!> module atmosphere / physics / lightninh / SATO2019
!!
!! @par Description
!!         Component for sato2019 tracer
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
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
  public :: ATMOS_PHY_LT_sato2019_tendency
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
  character(len=64),private,save    :: NUTR_qhyd='TOTAL'
  integer,  private                 :: MAX_NUTR = 100
  real(RP), private                 :: Eint = 150.0E+3_RP     ! [V/m]
  real(RP), private                 :: delEint = 10.0E+3_RP   ! [V/m]
  real(RP), private                 :: Estop = 15.0E+3_RP     ! [V/m]
  real(RP), private                 :: qrho_chan = 0.5_RP     ! [nC/m3]
  real(RP), private                 :: qrho_neut = 0.5_RP     ! [nC/m3]
  real(RP), private                 :: fp = 0.3_RP            ! [-]
  real(RP), private                 :: zcg = 500.0_RP         ! [m]
  integer,  private                 :: NUTR_ITMAX = 1000
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

  !---
  integer,  parameter, private :: nxlut_lt = 200, nylut_lt = 200
  real(RP), private :: dq_chrg( nxlut_lt,nylut_lt )    !--- charge separation [fC]
  real(RP), private :: grid_lut_t( nxlut_lt )
  real(RP), private :: grid_lut_l( nylut_lt )
  real(RP), private :: tcrglimit
  logical, private                  :: Hgt_dependency_Eint = .false.

  !--- Index for local variable
  integer, private, parameter :: I_lt_x = 1
  integer, private, parameter :: I_lt_y = 2
  integer, private, parameter :: I_lt_z = 3
  integer, private, parameter :: I_lt_abs = 4
  integer, private :: I_crg_LIQ_s, I_crg_LIQ_e
  integer, private :: I_crg_ICE_s, I_crg_ICE_e

  !--- For history output
  integer, private, parameter :: w_nmax = 13
  integer, private, parameter :: I_QCRG_LIQ = 1
  integer, private, parameter :: I_QCRG_ICE = 2
  integer, private, parameter :: I_QCRG_TOT = 3
  integer, private, parameter :: I_Ex = 4
  integer, private, parameter :: I_Ey = 5
  integer, private, parameter :: I_Ez = 6
  integer, private, parameter :: I_Eabs = 7
  integer, private, parameter :: I_Epot = 8
  integer, private, parameter :: I_Qneut = 9
  integer, private, parameter :: I_LTpath = 10
  integer, private, parameter :: I_PosFLASH = 11
  integer, private, parameter :: I_NegFLASH = 12
  integer, private, parameter :: I_FlashPoint = 13
  integer,  private              :: HIST_id(w_nmax)
  character(len=H_SHORT), private :: w_name(w_nmax)
  character(len=H_MID),   private :: w_longname(w_nmax)
  character(len=H_SHORT), private :: w_unit(w_nmax)
  data w_name / 'CRGD_LIQ', &
                'CRGD_ICE', &
                'CRGD_TOT',&
                'Ex', &
                'Ey', &
                'Ez', &
                'Eabs', &
                'Epot', &
                'Qneut', &
                'LTpath', &
                'PosFLASH', &
                'NegFLASH', &
                'FlashPoint' /
  data w_longname / &
                'Charge density of liquid water', &
                'Charge density of ice water', &
                'Charge density of QHYD', &
                'X component of Electrical Field', &
                'Y component of Electrical Field', &
                'Z component of Electrical Field', &
                'Absolute value of Electrical Field', &
                'Electric Potential', &
                'Cumulative Neutralizated charge', &
                'Cumulative Number of flash', &
                'Cumulative Number of Positive flash', &
                'Cumulative Number of Negative flash', &
                'Cumulative Number of Flash point' /
  data w_unit / 'nC/m3', &
                'nC/m3', &
                'nC/m3', &
                'kV/m', &
                'kV/m', &
                'kV/m', &
                'kV/m', &
                'V', &
                'nC/m3', &
                'num', &
                'num', &
                'num', &
                'num' /
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
                                          nliq, nice, &
                                          CDX, CDY    )
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster, &
       PRC_myrank
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_const, only: &
       T00 => CONST_TEM00
    use scale_file_history, only: &
       FILE_HISTORY_reg
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer, intent(in)  :: IA, IS, IE
    integer, intent(in)  :: JA, JS, JE
    integer, intent(in)  :: IMAXG
    integer, intent(in)  :: JMAXG
    integer, intent(in)  :: KMAX

    character(len=*), intent(in) :: MP_TYPE
    integer,          intent(in)  :: nice
    integer,          intent(in)  :: nliq
    real(RP),         intent(in)  :: CDX(IA)
    real(RP),         intent(in)  :: CDY(JA)

    integer :: n, myu, ip
    integer :: ierr

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
         LT_DO_Lightning

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

    if( R_neut <= minval( CDX,1 ) .or. R_neut <= minval( CDY,1 ) ) then
        R_neut = 2.d0 * min( minval( CDX,1 ), minval( CDY,1 ) )
    endif

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

    call COMM_bcast( dq_chrg, nxlut_lt,nylut_lt )
    call COMM_bcast( grid_lut_t,nxlut_lt )
    call COMM_bcast( grid_lut_l,nylut_lt )

!    KIJMAXG = (IEG-ISG+1)*(JEG-JSG+1)*(KE-KS+1)
    KIJMAXG = IMAXG*JMAXG*KMAX

    ! for history output
    allocate( d_QCRG_TOT(KA,IA,JA) )
    allocate( LT_PATH_TOT(KA,IA,JA,3) )
    allocate( fls_int_p_tot(KA,IA,JA) )
    d_QCRG_TOT(:,:,:) = 0.0_RP
    LT_PATH_TOT(:,:,:,:) = 0.0_RP
    fls_int_p_tot(:,:,:) = 0.0_RP

    tcrglimit = -60.0_RP+T00

    if( NUTR_TYPE == 'F2013' ) then
      LOG_INFO("ATMOS_PHY_LT_sato2019_setup",'(A,F15.7,A)') 'Radius of neutralization is ', R_neut, "[m]"
    endif

    flg_eint_hgt = 0.0_RP
    if( Hgt_dependency_Eint .and. NUTR_TYPE == 'MG2001') then
      flg_eint_hgt = 1.0_RP
    endif

    if( NUTR_qhyd /= 'TOTAL' .and. NUTR_qhyd /= 'Each_POLARITY' ) then
      LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'xxx NUTR_qhyd should be TOTAL or Each_POLARITY, stop!'
      call PRC_abort
    endif

    if( NUTR_qhyd == 'Each_POLARITY' .and. MP_TYPE == 'SUZUKI10' ) then
       LOG_ERROR("ATMOS_PHY_LT_sato2019_setup",*) 'NUTR_qhyd = Each_POLARITY is not supported for MP_TYPE SUZUKI10'
       call PRC_abort
    endif

    I_crg_LIQ_s = 1
    I_crg_LIQ_e = nliq
    I_crg_ICE_s = I_crg_LIQ_e + 1
    I_crg_ICE_e = I_crg_ICE_s + nice - 1

    do ip = 1, w_nmax
       call FILE_HISTORY_reg( w_name(ip), w_longname(ip), w_unit(ip), & ! [IN]
                              HIST_id(ip)                             ) ! [OUT]
    end do

    return
  end subroutine ATMOS_PHY_LT_sato2019_setup

  !-----------------------------------------------------------------------------
  !> Calculate tendency of charge density
  subroutine ATMOS_PHY_LT_sato2019_tendency( &
       KA, KS, KE,   &
       IA, IS, IE,   &
       JA, JS, JE,   &
       KIJMAX,       &
       IMAX,         &
       JMAX,         &
       QA_MP,        &
       QA_LT,        &
       DENS,         &
       RHOT,         &
       QTRC,         &
       QTRC_crg,     &
       dt_MP,        &
       dt_LT,        &
       Sarea,        &
       RHOQ_T_MP,    &
       RHOQ_t_LT_mp, &
       Epot,         &
       RHOQ_t_LT     )
    use scale_const, only: &
       SMALL => CONST_EPS
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
    use scale_comm_cartesC, only: &
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
    integer,  intent(in) :: QA_MP
    integer,  intent(in) :: QA_LT
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: Sarea(KA,IA,JA,QA_LT)       !--- Surface area and that of each catergory [m2]
    real(DP), intent(in) :: dt_MP
    real(DP), intent(in) :: dt_LT
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(in) :: QTRC_crg(KA,IA,JA,QA_LT)
    real(RP), intent(in) :: RHOQ_t_MP(KA,IA,JA,QA_MP)
    real(RP), intent(in) :: RHOQ_t_LT_mp(KA,IA,JA,QA_LT)
    real(RP), intent(inout) :: Epot(KA,IA,JA) !--- Electrical potential at previous time step [V]
    real(RP), intent(inout) :: RHOQ_t_LT(KA,IA,JA,QA_LT)

    real(RP) :: QHYD_mass(KA,IA,JA)         !--- Mass of total hydrometeor [kg/m3]
    real(RP) :: QTRC0(KA,IA,JA,QA_MP)       !--- Tracer after lightning component
    real(RP) :: QTRC0_crg(KA,IA,JA,QA_LT)   !--- Tracer after lightning component
    real(RP) :: QTRC1_crg(KA,IA,JA,QA_LT)   !--- Tracer after lightning component
    real(RP) :: QCRG(KA,IA,JA)              !--- Total charge density [nC/m3]
    real(RP) :: Epot_new(KA,IA,JA)          !--- Electrical potential [V]
    real(RP) :: Efield(KA,IA,JA,I_lt_abs)   !--- Electrical field (1-3 ->x,y,z, 4->abs. )
    real(RP) :: NUM_end(KA,IA,JA,3)         !--- Number of each flash type (1->negative, 2->ground, 3->positive)
    real(RP) :: d_QCRG(KA,IA,JA)            !--- Change of charge by charge neutralization [fC/m3]
    real(RP) :: LT_PATH(KA,IA,JA)           !--- Lightning path (0-> no path, 1-> path)
    real(RP) :: fls_int_p(KA,IA,JA)
    real(RP) :: Total_Sarea(2)              !--- Sum of surface area and that of each catergory [m2]
    real(RP) :: neg_crg, pos_crg, d_pos_crg, d_neg_crg
    real(RP) :: frac, flg_chrged(5), r_totalSarea(2)
    real(RP) :: Emax, Emax_old
    logical  :: output_step
    integer  :: flg_lt_neut
    integer  :: i, j, k, m, n, countbin, ip
    real(RP) :: zerosw, positive, negative

    logical  :: HIST_sw(w_nmax)
    real(RP) :: w3d(KA,IA,JA)

    NUM_end(:,:,:,:) = 0.0_RP

    ! Add tendency of charge density by microphysics
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC0(k,i,j,1:QA_MP) = QTRC(k,i,j,1:QA_MP) + RHOQ_t_MP(k,i,j,1:QA_MP)/DENS(k,i,j) * dt_MP
       QTRC0_crg(k,i,j,1:QA_LT) = QTRC_crg(k,i,j,1:QA_LT) + RHOQ_t_LT_mp(k,i,j,1:QA_LT)/DENS(k,i,j) * dt_MP
       QTRC1_crg(k,i,j,1:QA_LT) = QTRC0_crg(k,i,j,1:QA_LT)
    enddo
    enddo
    enddo

    ! calc total charge density
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QCRG(k,i,j) = 0.0_RP
       do n = 1, QA_LT
          QCRG(k,i,j) = QCRG(k,i,j) + QTRC0_crg(k,i,j,n)
       end do
       QCRG(k,i,j) = QCRG(k,i,j) * DENS(k,i,j) * 1.E-6_RP ![fC/kg] -> [nc/m3]
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QHYD_mass(k,i,j) = 0.0_RP
       do n = 2, QA_MP  ! ignore vapor (I_QV=1)
         QHYD_mass(k,i,j) = QHYD_mass(k,i,j) + QTRC0(k,i,j,n)*DENS(k,i,j) ![kg/kg] -> [kg/m3]
       enddo
    enddo
    enddo
    enddo

    !--- Calculate E field
    Efield(:,:,:,:) = 0.0_RP
    call ATMOS_PHY_LT_electric_field( KA, KS, KE,                   &   ! [IN]
                                      IA, IS, IE,                   &   ! [IN]
                                      JA, JS, JE,                   &   ! [IN]
                                      QCRG    (:,:,:),              &   ! [IN]
                                      DENS    (:,:,:),              &   ! [IN]
                                      RHOT    (:,:,:),              &   ! [IN]
                                      Epot    (:,:,:),              &   ! [IN]
                                      Epot_new(:,:,:),              &   ! [OUT]
                                      Efield  (:,:,:,I_lt_x:I_lt_z) )   ! [OUT]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Efield(k,i,j,I_lt_abs) = sqrt( Efield(k,i,j,I_lt_x)*Efield(k,i,j,I_lt_x) &
                                    + Efield(k,i,j,I_lt_y)*Efield(k,i,j,I_lt_y) &
                                    + Efield(k,i,j,I_lt_z)*Efield(k,i,j,I_lt_z) )
       Epot(k,i,j) = Epot_new(k,i,j)
    enddo
    enddo
    enddo


    d_QCRG(:,:,:) = 0.0_RP
    LT_PATH(:,:,:) = 0.0_RP
    fls_int_p(:,:,:) = 0.0_RP
    if ( LT_DO_Lightning ) then

       call ATMOS_PHY_LT_judge_absE( KA, KS, KE,             &   ! [IN]
                                     IA, IS, IE,             &   ! [IN]
                                     JA, JS, JE,             &   ! [IN]
                                     DENS(:,:,:),            &   ! [IN]
                                     Efield(:,:,:,I_lt_abs), &   ! [IN]
                                     Emax,                   &   ! [OUT]
                                     flg_lt_neut             )   ! [OUT]

       if( flg_lt_neut > 0 .and. PRC_IsMaster ) then
          LOG_INFO("ATMOS_PHY_LT_sato2019_calc_tendency",'(F15.7,A)')  &
              Emax*1.E-3_RP, '[kV/m] Charge Neutralization is calculated'
       endif

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
                                             fls_int_p(:,:,:),   & !  [INOUT]
                                             d_QCRG   (:,:,:)    ) !  [OUT]
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
                                             fls_int_p(:,:,:),   & !  [INOUT]
                                             d_QCRG   (:,:,:)    ) !  [OUT]
         endif

         call COMM_vars8( LT_path(:,:,:),1 )
         call COMM_vars8( d_QCRG(:,:,:),2 )
         call COMM_wait ( LT_path(:,:,:),1 )
         call COMM_wait ( d_QCRG(:,:,:),2 )

         !-- Calculate neutralization of each hydrometeor or each category
         select case( NUTR_qhyd )
         case ( 'TOTAL' )

            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               if( abs( d_QCRG(k,i,j) ) > 0.0_RP ) then
                  Total_Sarea(1) = 0.0_RP
                  do n = 1, QA_LT
                     Total_Sarea(1) = Total_Sarea(1) + Sarea(k,i,j,n)
                  enddo
                  zerosw = 0.5_RP - sign( 0.5_RP, Total_Sarea(1)-SMALL )
                  r_totalSarea(1) = 1.0_RP / ( Total_Sarea(1) + zerosw ) * ( 1.0_RP - zerosw )
                  do n = 1, QA_LT
                     QTRC1_crg(k,i,j,n) = QTRC0_crg(k,i,j,n)  &
                                        + ( d_QCRG(k,i,j)*1.0E+6_RP )  &
                                        * Sarea(k,i,j,n) * r_totalSarea(1) / DENS(k,i,j)
                  enddo
               endif
            enddo
            enddo
            enddo

         case ( 'Each_POLARITY' )

             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                if( abs( d_QCRG(k,i,j) ) > 0.0_RP ) then
                  Total_Sarea(:) = 0.0_RP
                  flg_chrged(:) = 0.0_RP

                  !--- flg whether the charged or not (0.0-> not charged, 1.0->charged)
                  do n = 1, QA_LT
                     flg_chrged(n) = 0.5_RP + sign( 0.5_RP, abs(QTRC0_crg(k,i,j,n))-SMALL )
                  enddo

                  pos_crg = 0.0_RP
                  neg_crg = 0.0_RP
                  do n = 1, QA_LT
                     positive = 0.5_RP + sign( 0.5_RP, QTRC0_crg(k,i,j,n) )
                     negative = 1.0_RP - positive
                     !--- total of positive charge
                     pos_crg = pos_crg + QTRC0_crg(k,i,j,n) * positive * flg_chrged(n)
                     !--- total of negative charge
                     neg_crg = neg_crg - QTRC0_crg(k,i,j,n) * negative * flg_chrged(n)
                     !--- Sarea of positively charged hydrometeor
                     Total_Sarea(1) = Total_Sarea(1) + Sarea(k,i,j,n) * positive * flg_chrged(n)
                     !--- Sarea of negatively charged hydrometeor
                     Total_Sarea(2) = Total_Sarea(2) + Sarea(k,i,j,n) * negative * flg_chrged(n)
                  end do

                  zerosw = 0.5_RP - sign( 0.5_RP, abs( QCRG(k,i,j) ) - SMALL )
                  frac = d_QCRG(k,i,j) / ( QCRG(k,i,j) + zerosw ) * ( 1.0_RP - zerosw )
                  frac = frac * DENS(k,i,j) * 1.E-6_RP ! [fC/kg] -> [nc/m3]
                  d_pos_crg = frac * pos_crg
                  d_neg_crg = frac * neg_crg

                  !--- remove 0 surface area ( no crg for each porality )
                  zerosw = 0.5_RP - sign( 0.5_RP, Total_Sarea(1) - SMALL )
                  r_totalSarea(1) = 1.0_RP / ( Total_Sarea(1) + zerosw ) * ( 1.0_RP - zerosw )
                  zerosw = 0.5_RP - sign( 0.5_RP, Total_Sarea(2) - SMALL )
                  r_totalSarea(2) = 1.0_RP / ( Total_Sarea(2) + zerosw ) * ( 1.0_RP - zerosw )

                  do n = 1, QA_LT
                     QTRC1_crg(k,i,j,n) = QTRC0_crg(k,i,j,n)  &
                                        + ( d_pos_crg*1.0E+6_RP ) / DENS(k,i,j)  &
                                        * Sarea(k,i,j,n) * r_totalSarea(1) &
                                        * ( 0.5_RP + sign( 0.5_RP, QTRC0_crg(k,i,j,n) ) ) &
                                        * flg_chrged(n) &
                                        + ( d_neg_crg*1.0E+6_RP ) / DENS(k,i,j) &
                                        * Sarea(k,i,j,n) * r_totalSarea(2)  &
                                        * ( 0.5_RP - sign( 0.5_RP, QTRC0_crg(k,i,j,n) ) ) &
                                        * flg_chrged(n)
                  enddo
               endif

            enddo
            enddo
            enddo

         end select


         !--- Calculate E field
         Efield(:,:,:,:) = 0.0_RP
         call ATMOS_PHY_LT_electric_field( KA, KS, KE,                   &   ! [IN]
                                           IA, IS, IE,                   &   ! [IN]
                                           JA, JS, JE,                   &   ! [IN]
                                           QCRG    (:,:,:),              &   ! [IN]
                                           DENS    (:,:,:),              &   ! [IN]
                                           RHOT    (:,:,:),              &   ! [IN]
                                           Epot    (:,:,:),              &   ! [IN]
                                           Epot_new(:,:,:),              &   ! [OUT]
                                           Efield  (:,:,:,I_lt_x:I_lt_z) )   ! [OUT]

         !--- Add Total number of charge neutralization and flash point
         if ( HIST_id(I_Qneut) > 0 ) then
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               d_QCRG_TOT(k,i,j) = d_QCRG_TOT(k,i,j) + d_QCRG(k,i,j)
            end do
            end do
            end do
         end if
         if ( HIST_id(I_FlashPoint) > 0 ) then
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               fls_int_p_tot(k,i,j) = fls_int_p_tot(k,i,j) + fls_int_p(k,i,j)
            end do
            end do
            end do
         end if

         !--- Add Total number of path
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            Efield(k,i,j,I_lt_abs) = sqrt( Efield(k,i,j,I_lt_x)*Efield(k,i,j,I_lt_x) &
                                         + Efield(k,i,j,I_lt_y)*Efield(k,i,j,I_lt_y) &
                                         + Efield(k,i,j,I_lt_z)*Efield(k,i,j,I_lt_z) )
            Epot(k,i,j) = Epot_new(k,i,j)
         end do
         end do
         end do

         if ( HIST_id(I_PosFLASH) > 0 ) then
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               LT_PATH_TOT(k,i,j,1) = LT_PATH_TOT(k,i,j,1) &
                                    + 0.5_RP + sign( 0.5_RP,-d_QCRG(k,i,j)-SMALL )
            end do
            end do
            end do
         end if
         if ( HIST_id(I_NegFLASH) > 0 ) then
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               LT_PATH_TOT(k,i,j,2) = LT_PATH_TOT(k,i,j,2) &
                                    + 0.5_RP + sign( 0.5_RP, d_QCRG(k,i,j)-SMALL )
            end do
            end do
            end do
         end if
         if ( HIST_id(I_LTpath) > 0 ) then
            do j = JS, JE
            do i = IS, IE
            do k = KS, KE
               LT_PATH_TOT(k,i,j,3) = LT_PATH_TOT(k,i,j,3) + LT_PATH(k,i,j)
            end do
            end do
            end do
         end if

         call ATMOS_PHY_LT_judge_absE( KA, KS, KE,             &   ! [IN]
                                       IA, IS, IE,             &   ! [IN]
                                       JA, JS, JE,             &   ! [IN]
                                       DENS(:,:,:),            &   ! [IN]
                                       Efield(:,:,:,I_lt_abs), &   ! [IN]
                                       Emax,                   &   ! [OUT]
                                       flg_lt_neut             )   ! [OUT]

         if( flg_lt_neut == 1 .and. Emax == Emax_old ) then
           flg_lt_neut = 0
         elseif( flg_lt_neut == 2 ) then
           flg_lt_neut = 0
         endif

       enddo
    endif

    !--- For history output
    do ip = 1, w_nmax
       call FILE_HISTORY_query( HIST_id(ip), HIST_sw(ip) )
    end do

    if ( HIST_sw(I_QCRG_LIQ) ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = 0.0_RP
          do n = I_crg_LIQ_s, I_crg_LIQ_e
             w3d(k,i,j) = w3d(k,i,j) + QTRC1_crg(k,i,j,n)
          enddo
          w3d(k,i,j) = w3d(k,i,j) * DENS(k,i,j) * 1.E-6_RP ![fC/kg] -> [nc/m3]
       end do
       end do
       end do
       call FILE_HISTORY_put( HIST_id(I_QCRG_LIQ), w3d(:,:,:) )
    end if
    if ( HIST_sw(I_QCRG_ICE) ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = 0.0_RP
          do n = I_crg_ICE_s, I_crg_ICE_e
             w3d(k,i,j) = w3d(k,i,j) + QTRC1_crg(k,i,j,n)
          enddo
          w3d(k,i,j) = w3d(k,i,j) * DENS(k,i,j) * 1.E-6_RP ![fC/kg] -> [nc/m3]
       end do
       end do
       end do
       call FILE_HISTORY_put( HIST_id(I_QCRG_ICE), w3d(:,:,:) )
    end if
    if ( HIST_sw(I_QCRG_TOT) ) then
       call FILE_HISTORY_put( HIST_id(I_QCRG_TOT), QCRG(:,:,:) )
    end if
    if ( HIST_sw(I_Ex  ) ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_x  )*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       call FILE_HISTORY_put( HIST_id(I_Ex  ), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Ey  ) ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_y  )*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       call FILE_HISTORY_put( HIST_id(I_Ey  ), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Ez  ) ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_z  )*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       call FILE_HISTORY_put( HIST_id(I_Ez  ), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Eabs) ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          w3d(k,i,j) = Efield(k,i,j,I_lt_abs)*1.0E-3_RP ![kV/m]
       enddo
       enddo
       enddo
       call FILE_HISTORY_put( HIST_id(I_Eabs), w3d(:,:,:) )
    endif
    if ( HIST_sw(I_Epot) ) then
       call FILE_HISTORY_put( HIST_id(I_Epot), Epot(:,:,:) )
    end if
    if ( HIST_sw(I_Qneut) ) then
       call FILE_HISTORY_put( HIST_id(I_Qneut), d_QCRG_TOT(:,:,:) )
       d_QCRG_tot(:,:,:) = 0.0_RP
    endif
    if ( HIST_sw(I_LTpath)  ) then
       call FILE_HISTORY_put( HIST_id(I_LTpath), LT_PATH_TOT(:,:,:,3) )
       LT_PATH_TOT(:,:,:,3) = 0.0_RP
    endif
    if ( HIST_sw(I_PosFLASH) ) then
       call FILE_HISTORY_put( HIST_id(I_PosFLASH), LT_PATH_TOT(:,:,:,1) )
       LT_PATH_TOT(:,:,:,1) = 0.0_RP
    endif
    if ( HIST_sw(I_NegFLASH) ) then
       call FILE_HISTORY_put( HIST_id(I_NegFLASH), LT_PATH_TOT(:,:,:,2) )
       LT_PATH_TOT(:,:,:,2) = 0.0_RP
    endif
    if ( HIST_sw(I_FlashPoint) ) then
       call FILE_HISTORY_put( HIST_id(I_FlashPoint), fls_int_p_tot(:,:,:) )
       fls_int_p_tot(:,:,:) = 0.0_RP
    end if


    !--- Calculation of tendency
    do n = 1, QA_LT
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOQ_t_LT(k,i,j,n) = ( QTRC1_crg(k,i,j,n)-QTRC0_crg(k,i,j,n) )*DENS(k,i,j)/dt_LT
       enddo
       enddo
       enddo
    enddo

    return
  end subroutine ATMOS_PHY_LT_sato2019_tendency
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
       E_pot_old,  & ! [IN]
       E_pot,      & ! [OUT]
       Efield      ) ! [OUT]
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

    integer,  intent(in)    :: KA, KS, KE
    integer,  intent(in)    :: IA, IS, IE
    integer,  intent(in)    :: JA, JS, JE
    real(RP), intent(in)  :: QCRG     (KA,IA,JA)      !-- Charge density [nC/m3]
    real(RP), intent(in)  :: DENS     (KA,IA,JA)      !-- Total density [kg/m3]
    real(RP), intent(in)  :: RHOT     (KA,IA,JA)      !-- density weighted potential temperature [K kg/m3]
    real(RP), intent(in)  :: E_pot_old(KA,IA,JA)      !-- Electric potential in previous step[V]
    real(RP), intent(out) :: E_pot    (KA,IA,JA)      !-- Electric potential [V]
    real(RP), intent(out) :: Efield   (KA,IA,JA,3)    !-- Electric field [V/m]

    real(RP) :: eps_air(KA,IA,JA)
    !--- A x E_pott = - QCRG/epsiron
    real(RP) :: A(15,KA,IA,JA)            !--- A : Laplasian
    real(RP) :: B(KA,IA,JA)               !--- B : -QCRG*DENS/epsiron
    real(RP) :: E_pot_N(KA,IA,JA)         !--- electrical potential calculated by Bi-CGSTAB

    integer :: i, j, k, ijk, ierror
    real(RP) :: iprod, buf

    call PROF_rapstart('LT_E_field', 1)

    iprod = 0.0_RP
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
        iprod = iprod + abs( QCRG(k,i,j) )
    enddo
    enddo
    enddo

    call MPI_AllReduce(iprod, buf, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

!    if( buf == 0.0_RP ) then
    if( buf <= EPS ) then
      E_pot(:,:,:) = 0.0_RP
      Efield(:,:,:,1:3) = 0.0_RP
      return
    endif

    eps_air(:,:,:) = EPSvac * EPSair
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       eps_air(k,i,j) = EPSvac * EPSair   !--- temporary, dependency of epsiron on P and T will be implemented
       B(k,i,j) = - QCRG(k,i,j)/eps_air(k,i,j) * 1.0E-9_RP ! [nC/m3] -> [C/m3 * m/F] = [V/m2]
    enddo
    enddo
    enddo

    !---- fill halo
    call COMM_vars8( eps_air,1 )
    call COMM_vars8( B,      2 )
    call COMM_wait ( eps_air,1 )
    call COMM_wait ( B,      2 )

    !---- input vector A
    do i = IS, IE
    do j = JS, JE
     do k = KS, KE
       ! (k,i,j)
       A(1,k,i,j) = &
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
       A(2,k,i,j) = &
                  - MAPF(i,j,1,I_XY)*J13G(k,i  ,j,I_UYZ)*f2h(k-1,i,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*J13G(k-1,i,j,I_XYW)*RFDZ(k)*RCDZ(k-1) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j  ,I_XVZ)*f2h(k-1,i,j,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*J23G(k-1,i,j,I_XYW)*RFDZ(k)*RCDZ(k-1) &
                  + J33G*J33G*RFDZ(k)*RCDZ(k-1)

       ! (k+1,i,j)
       A(3,k,i,j) = &
                    MAPF(i,j,1,I_XY)*J13G(k,i  ,j,I_UYZ)*f2h(k,i,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k,i,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*J13G(k,i,j,I_XYW)*RFDZ(k)*RCDZ(k) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j  ,I_XVZ)*f2h(k,i,j,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k,i,j,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*J23G(k,i,j,I_XYW)*RFDZ(k)*RCDZ(k) &
                  + J33G*J33G*RFDZ(k)*RCDZ(k)

       ! (k,i-1,j)
       A(4,k,i,j) = &
                    MAPF(i,j,1,I_XY)*MAPF(i-1,j,1,I_XY)*RFDX(i)*RCDX(i-1) &
                  - MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k  ,i-1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i-1,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k  ,i-1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i-1,j,1)*0.50_RP*RFDX(i)*RFDZ(k)

       ! (k,i+1,j)
       A(5,k,i,j) = &
                    MAPF(i,j,1,I_XY)*MAPF(i+1,j,1,I_XY)*RFDX(i)*RCDX(i) &
                  + MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k  ,i+1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k-1,i+1,j,1)*0.50_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k  ,i+1,j,2)*0.50_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i+1,j,1)*0.50_RP*RFDX(i)*RFDZ(k)

       ! (k,i,j-1)
       A(6,k,i,j) = &
                    MAPF(i,j,2,I_XY)*MAPF(i,j-1,2,I_XY)*RFDY(j)*RCDY(j-1) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k  ,i,j-1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j-1,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k  ,i,j-1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j-1,1)*0.50_RP*RFDY(j)*RFDZ(k)

       ! (k,i,j+1)
       A(7,k,i,j) = &
                    MAPF(i,j,2,I_XY)*MAPF(i,j+1,2,I_XY)*RFDY(j)*RCDY(j) &
                  + MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k  ,i,j+1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k-1,i,j+1,1)*0.50_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k  ,i,j+1,2)*0.50_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j+1,1)*0.50_RP*RFDY(j)*RFDZ(k)

       ! (k-1,i-1,j)
       A(8,k,i,j) = &
                    MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k-1,i-1,j,2)*0.5_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i-1,j,2)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k-1,i+1,j)
       A(9,k,i,j) = &
                  - MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k-1,i+1,j,2)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k-1,i+1,j,2)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k-1,i,j-1)
       A(10,k,i,j) = &
                    MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k-1,i,j-1,2)*0.5_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j-1,2)*0.5_RP*RFDY(j)*RFDZ(k)

       ! (k-1,i,j+1)
       A(11,k,i,j) = &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k-1,i,j+1,2)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k-1,i,j+1,2)*0.5_RP*RFDY(j)*RFDZ(k)

       ! (k+1,i-1,j)
       A(12,k,i,j) = &
                  - MAPF(i,j,1,I_XY)*J13G(k,i-1,j,I_UYZ)*f2h(k,i-1,j,1)*0.5_RP*RFDX(i)*RFDZ(k) &
                  - J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k,i-1,j,1)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k+1,i+1,j)
       A(13,k,i,j) = &
                    MAPF(i,j,1,I_XY)*J13G(k,i,j,I_UYZ)*f2h(k,i+1,j,1)*0.5_RP*RFDX(i)*RFDZ(k) &
                  + J13G(k,i,j,I_XYZ)*MAPF(i,j,1,I_XY)*f2h(k,i+1,j,1)*0.5_RP*RFDX(i)*RFDZ(k)

       ! (k+1,i,j-1)
       A(14,k,i,j) = &
                  - MAPF(i,j,2,I_XY)*J23G(k,i,j-1,I_XVZ)*f2h(k,i,j-1,1)*0.5_RP*RFDY(j)*RFDZ(k) &
                  - J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k,i,j-1,1)*0.5_RP*RFDY(j)*RFDZ(k)

       ! (k+1,i,j+1)
       A(15,k,i,j) = &
                    MAPF(i,j,2,I_XY)*J23G(k,i,j,I_XVZ)*f2h(k,i,j+1,1)*0.5_RP*RFDY(j)*RFDZ(k) &
                  + J23G(k,i,j,I_XYZ)*MAPF(i,j,2,I_XY)*f2h(k,i,j+1,1)*0.5_RP*RFDY(j)*RFDZ(k)

     enddo
    enddo
    enddo

    !---- fill halo
    call COMM_vars8( A(1,:,:,:), 1  )
    call COMM_vars8( A(2,:,:,:), 2  )
    call COMM_vars8( A(3,:,:,:), 3  )
    call COMM_vars8( A(4,:,:,:), 4  )
    call COMM_vars8( A(5,:,:,:), 5  )
    call COMM_vars8( A(6,:,:,:), 6  )
    call COMM_vars8( A(7,:,:,:), 7  )
    call COMM_vars8( A(8,:,:,:), 8  )
    call COMM_vars8( A(9,:,:,:), 9  )
    call COMM_vars8( A(10,:,:,:),10 )
    call COMM_vars8( A(11,:,:,:),11 )
    call COMM_vars8( A(12,:,:,:),12 )
    call COMM_vars8( A(13,:,:,:),13 )
    call COMM_vars8( A(14,:,:,:),14 )
    call COMM_vars8( A(15,:,:,:),15 )
    call COMM_wait ( A(1,:,:,:), 1  )
    call COMM_wait ( A(2,:,:,:), 2  )
    call COMM_wait ( A(3,:,:,:), 3  )
    call COMM_wait ( A(4,:,:,:), 4  )
    call COMM_wait ( A(5,:,:,:), 5  )
    call COMM_wait ( A(6,:,:,:), 6  )
    call COMM_wait ( A(7,:,:,:), 7  )
    call COMM_wait ( A(8,:,:,:), 8  )
    call COMM_wait ( A(9,:,:,:), 9  )
    call COMM_wait ( A(10,:,:,:),10 )
    call COMM_wait ( A(11,:,:,:),11 )
    call COMM_wait ( A(12,:,:,:),12 )
    call COMM_wait ( A(13,:,:,:),13 )
    call COMM_wait ( A(14,:,:,:),14 )
    call COMM_wait ( A(15,:,:,:),15 )

    E_pot_N(:,:,:) = E_pot_old(:,:,:)   !-- initial value -> previous step value
    call COMM_vars8( E_pot_N, 1 )
    call COMM_wait ( E_pot_N, 1 )

    !--- calcuclate counter matrix
    call solve_bicgstab( &
       KA, KS, KE, & ! (in)
       IA, IS, IE, & ! (in)
       JA, JS, JE, & ! (in)
       E_pot,      & ! (out)
       E_pot_n,    & ! (in)
       A, B        ) ! (in)

    call COMM_vars8( E_pot, 1 )
    call COMM_wait ( E_pot, 1 )

    !--- boundary condition
    if( .not. PRC_HAS_W ) then
      do k = 1, KA
      do j = 1, JA
        E_pot(k,1:IS-1,j) = E_pot(k,IS,j)
      enddo
      enddo
    endif
    if( .not. PRC_HAS_E ) then
      do k = 1, KA
      do j = 1, JA
        E_pot(k,IE+1:IA,j) = E_pot(k,IE,j)
      enddo
      enddo
    endif
    if( .not. PRC_HAS_S ) then
      do k = 1, KA
      do i = 1, IA
        E_pot(k,i,1:JS-1) = E_pot(k,i,JS)
      enddo
      enddo
    endif
    if( .not. PRC_HAS_N ) then
      do k = 1, KA
      do i = 1, IA
        E_pot(k,i,JE+1:JA) = E_pot(k,i,JE)
      enddo
      enddo
    endif
    do i = 1, IA
    do j = 1, JA
      E_pot(1:KS-1,i,j) = 0.0_RP
      E_pot(KE+1:KA,i,j) = 0.0_RP
    enddo
    enddo

    !---- Calculate Electrical Field
    Efield(:,:,:,:) = 0.d0
    do i = IS, IE
    do j = JS, JE
     do k = KS, KE
       Efield(k,i,j,I_lt_x) = - MAPF(i,j,1,I_XYZ)/GSQRT(k,i,j,I_XYZ) * &
                            ( &
                            ( GSQRT(k,i+1,j,I_XYZ)*E_pot(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*E_pot(k,i-1,j) ) &
                            * RCDX(i) * 0.5_RP &
                            + ( J13G(k+1,i,j,I_UYW)*GSQRT(k+1,i,j,I_UYW)*E_pot(k+1,i,j) &
                              - J13G(k-1,i,j,I_UYW)*GSQRT(k-1,i,j,I_UYW)*E_pot(k-1,i,j) &
                              ) &
                            * RFDZ(k) * 0.5_RP  &
                            )
       Efield(k,i,j,I_lt_y) = - MAPF(i,j,2,I_XYZ)/GSQRT(k,i,j,I_XYZ) * &
                            ( &
                            ( GSQRT(k,i,j+1,I_XYZ)*E_pot(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*E_pot(k,i,j-1) ) &
                            * RCDY(j) * 0.5_RP &
                            + ( J23G(k+1,i,j,I_UYW)*GSQRT(k+1,i,j,I_UYW)*E_pot(k+1,i,j) &
                              - J23G(k-1,i,j,I_UYW)*GSQRT(k-1,i,j,I_UYW)*E_pot(k-1,i,j) &
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

    call PROF_rapend('LT_E_field', 1)

    return
  end subroutine ATMOS_PHY_LT_electric_field

  subroutine solve_bicgstab( &
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
    real(RP), intent(in)  :: M(15,KA,IA,JA)
    real(RP), intent(in)  :: B(KA,IA,JA)

    real(RP) :: r0(KA,IA,JA)

    real(RP) :: p(KA,IA,JA)
    real(RP) :: Mp(KA,IA,JA)
    real(RP) :: s(KA,IA,JA)
    real(RP) :: Ms(KA,IA,JA)
    real(RP) :: al, be, w

    real(RP), pointer :: r(:,:,:)
    real(RP), pointer :: rn(:,:,:)
    real(RP), pointer :: swap(:,:,:)
    real(RP), target :: v0(KA,IA,JA)
    real(RP), target :: v1(KA,IA,JA)
    real(RP) :: r0r
    real(RP) :: norm, error, error2

    real(RP) :: iprod(2)
    real(RP) :: buf(2)

    integer :: k, i, j
    integer :: iis, iie, jjs, jje
    integer :: iter
    integer :: ierror

    r  => v0
    rn => v1

    call mul_matrix( KA, KS, KE, & ! (in)
                     IA, IS, IE, & ! (in)
                     JA, JS, JE, & ! (in)
                     v1, M, PHI  ) ! v1 = M x0

    norm = 0.0_RP
    r0r  = 0.0_RP

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       norm = norm + B(k,i,j)**2
    enddo
    enddo
    enddo

    ! r = b - M x0
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       r(k,i,j) = B(k,i,j) - v1(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       r0(k,i,j) = r(k,i,j)
       p(k,i,j) = r(k,i,j)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       r0r = r0r + r0(k,i,j) * r(k,i,j)
    enddo
    enddo
    enddo

    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       PHI_N(k,i,j) = PHI(k,i,j)
    end do
    end do
    end do


    iprod(1) = r0r
    iprod(2) = norm
    call MPI_AllReduce(iprod, buf, 2, COMM_datatype, MPI_SUM, COMM_world, ierror)
    r0r = buf(1)
    norm = buf(2)

    error2 = norm

    do iter = 1, ITMAX

       error = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          error = error + r(k,i,j)**2
       enddo
       enddo
       enddo
       call MPI_AllReduce(error, buf, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)
       error = buf(1)

       if ( sqrt(error/norm) < epsilon ) then
         LOG_INFO("ATMOS_PHY_LT_Efield",'(a,1x,i0,1x,2e15.7)') "Bi-CGSTAB converged:", iter, sqrt(error/norm),norm
         exit
       endif
       error2 = error

       call COMM_vars8( p, 1 )
       call COMM_wait ( p, 1 )
       call mul_matrix( KA, KS, KE, & ! (in)
                        IA, IS, IE, & ! (in)
                        JA, JS, JE, & ! (in)
                        Mp, M, p    )

       iprod(1) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          iprod(1) = iprod(1) + r0(k,i,j) * Mp(k,i,j)
       enddo
       enddo
       enddo
       call MPI_AllReduce(iprod, buf, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)
       al = r0r / buf(1) ! (r0,r) / (r0,Mp)

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          s(k,i,j) = r(k,i,j) - al*Mp(k,i,j)
       enddo
       enddo
       enddo

       call COMM_vars8( s, 1 )
       call COMM_wait ( s, 1 )
       call mul_matrix( KA, KS, KE, & ! (in)
                        IA, IS, IE, & ! (in)
                        JA, JS, JE, & ! (in)
                        Ms, M, s )
       iprod(1) = 0.0_RP
       iprod(2) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          iprod(1) = iprod(1) + Ms(k,i,j) *  s(k,i,j)
          iprod(2) = iprod(2) + Ms(k,i,j) * Ms(k,i,j)
       enddo
       enddo
       enddo
       call MPI_AllReduce(iprod, buf, 2, COMM_datatype, MPI_SUM, COMM_world, ierror)
       w = buf(1) / buf(2) ! (Ms,s) / (Ms,Ms)

       iprod(1) = 0.0_RP

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PHI_N(k,i,j) = PHI_N(k,i,j) + al*p(k,i,j) + w*s(k,i,j)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          rn(k,i,j) = s(k,i,j) - w*Ms(k,i,j)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          iprod(1) = iprod(1) + r0(k,i,j) * rn(k,i,j)
       enddo
       enddo
       enddo

       be = al/w / r0r

       call MPI_AllReduce(iprod, r0r, 1, COMM_datatype, MPI_SUM, COMM_world, ierror)

       be = be * r0r ! al/w * (r0,rn)/(r0,r)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          p(k,i,j) = rn(k,i,j) + be * ( p(k,i,j) - w*Mp(k,i,j) )
       enddo
       enddo
       enddo

       swap => rn
       rn => r
       r => swap
    enddo

    if ( iter >= ITMAX ) then
       if( PRC_IsMaster ) then
         write(*,*) 'xxx [atmos_phy_lt] Bi-CGSTAB'
         write(*,*) 'xxx not converged', error, norm
         write(*,*) 'xxx epsilon(set,last)=', epsilon, sqrt(error/norm)
         if( error /= error ) then
          write(*,*) 'xxx error or norm is NaN Stop!'
          call PRC_abort
         endif
       endif
    endif

    return
  end subroutine solve_bicgstab

  subroutine mul_matrix( KA,KS,KE, &
                         IA,IS,IE, &
                         JA,JS,JE, &
                         V, M, C)
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(out) :: V(KA,IA,JA)
    real(RP), intent(in)  :: M(15,KA,IA,JA)
    real(RP), intent(in)  :: C(KA,IA,JA)

    integer :: k, i, j

    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE-1
       V(k,i,j) = M(1,k,i,j) * C(k  ,i  ,j  ) &
                + M(2,k,i,j) * C(k-1,i  ,j  ) &
                + M(3,k,i,j) * C(k+1,i  ,j  ) &
                + M(4,k,i,j) * C(k  ,i-1,j  ) &
                + M(5,k,i,j) * C(k  ,i+1,j  ) &
                + M(6,k,i,j) * C(k  ,i  ,j-1) &
                + M(7,k,i,j) * C(k  ,i  ,j+1) &
                + M(8,k,i,j) * C(k-1,i-1,j  ) &
                + M(9,k,i,j) * C(k-1,i+1,j  ) &
                + M(10,k,i,j)* C(k-1,i  ,j-1) &
                + M(11,k,i,j)* C(k-1,i  ,j+1) &
                + M(12,k,i,j)* C(k+1,i-1,j  ) &
                + M(13,k,i,j)* C(k+1,i+1,j  ) &
                + M(14,k,i,j)* C(k+1,i  ,j-1) &
                + M(15,k,i,j)* C(k+1,i  ,j+1)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       V(KS,i,j) = M(1,KS,i,j) * C(KS  ,i  ,j  ) &
                 + M(3,KS,i,j) * C(KS+1,i  ,j  ) &
                 + M(4,KS,i,j) * C(KS  ,i-1,j  ) &
                 + M(5,KS,i,j) * C(KS  ,i+1,j  ) &
                 + M(6,KS,i,j) * C(KS  ,i  ,j-1) &
                 + M(7,KS,i,j) * C(KS  ,i  ,j+1) &
                 + M(12,KS,i,j)* C(KS+1,i-1,j  ) &
                 + M(13,KS,i,j)* C(KS+1,i+1,j  ) &
                 + M(14,KS,i,j)* C(KS+1,i  ,j-1) &
                 + M(15,KS,i,j)* C(KS+1,i  ,j+1)
       V(KE,i,j) = M(1,KE,i,j) * C(KE  ,i  ,j  ) &
                 + M(2,KE,i,j) * C(KE-1,i  ,j  ) &
                 + M(4,KE,i,j) * C(KE  ,i-1,j  ) &
                 + M(5,KE,i,j) * C(KE  ,i+1,j  ) &
                 + M(6,KE,i,j) * C(KE  ,i  ,j-1) &
                 + M(7,KE,i,j) * C(KE  ,i  ,j+1) &
                 + M(8,KE,i,j) * C(KE-1,i-1,j  ) &
                 + M(9,KE,i,j) * C(KE-1,i+1,j  ) &
                 + M(10,KE,i,j)* C(KE-1,i  ,j-1) &
                 + M(11,KE,i,j)* C(KE-1,i  ,j+1)
    enddo
    enddo

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
    real(RP),allocatable :: randnum(:), randnum_3d(:,:,:)
    real(RP) :: rprod1, rprod2, rbuf, rbuf2, phi_end(2)
    integer, allocatable :: proc_num(:), proc_numg(:)
    integer  :: Eovid(4,KIJMAX)
    integer  :: Npls, Nmns, Ntot
    integer  :: count1(PRC_nprocs)
    integer  :: countindx(PRC_nprocs+1)
    integer  :: own_prc_total                         !--- number of grid whose |E| > Eint-dEint for each process
    integer  :: iprod1, iprod2, iprod3(4), buf
    integer  :: init_point(1), rank_initpoint, grid_initpoint
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

    fls_int_p(:,:,:) = 0.0_RP

    do count_path = 1, NUTR_ITMAX


       Eovid(:,:) = 0
       !--- search grid whose Efield is over threshold of flash inititation
       own_prc_total = 0
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
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

       !**** proc_num(0~) -> process number of each grid with |E|> E_threthold (local)
       allocate(proc_num(own_prc_total))
       proc_num(:) = PRC_myrank
       call MPI_AllGather( own_prc_total, 1, MPI_integer, &
                           count1, 1, MPI_integer, &
                           COMM_world, ierr )
       countindx(1) = 0

       do ipp = 1, PRC_nprocs
          countindx(ipp+1) = countindx(ipp) + count1(ipp)
       enddo

       !---- Create global version of proc_num(proc_numg)
       !**** countindx(PROC_nprocs) -> total number of grid with |E|>E_threthold
       !**** proc_numg(0~) -> process number of each grid with |E|> E_threthold (global)
       allocate(randnum(countindx(PRC_nprocs+1)))
       allocate(proc_numg(countindx(PRC_nprocs+1)))

       call MPI_AllGatherv( proc_num, own_prc_total, MPI_integer, &
                            proc_numg, count1, countindx, MPI_integer, &
                            COMM_world, ierr )

       !---- select initial point of flash by random select
       !**** rank_initpoint -> rank number including initpoint
       !**** grid_initpoint -> grid number of init point in rank (rank_initpoint)
       if( PRC_IsMaster ) then
          call RANDOM_uniform(randnum)
          init_point(1) = minloc( randnum,1 )
          randnum(:) = 0.0_RP
          randnum(init_point(1)) = 1.0_RP
          rank_initpoint = proc_numg(init_point(1))
          grid_initpoint = 0
          do i = init_point(1), 0, -1
             grid_initpoint = grid_initpoint + 1
             if( i == 0 ) then
                grid_initpoint = grid_initpoint - 1
             elseif( i /= 0 .and. proc_numg(i) /= rank_initpoint ) then
                grid_initpoint = grid_initpoint - 1
                exit
             endif
          enddo
       endif

       call COMM_bcast( rank_initpoint )
       call COMM_bcast( grid_initpoint )

       deallocate(randnum)
       deallocate(proc_num)
       deallocate(proc_numg)

       L_path(:,:,:) = 0.0_RP   !--- +-1 -> path already passed, +-2 -> path calculate current step

       !--- Propagate lightning
       flg_num_end(:,:,:,:) = 0
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
             fls_int_p(k,i,j) = fls_int_p(k,i,j) + 1.0_RP
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
             iprod1 = stop_flag
             call MPI_AllReduce(iprod1, iprod2, 1, MPI_integer, MPI_MAX, COMM_world, ierr)
             stop_flag = iprod2

             if( stop_flag == 1 ) then
                !---- send flag wildfire
                iprod1 = wild_flag
                call MPI_AllReduce(iprod1, iprod2, 1, MPI_integer, MPI_MAX, COMM_world, ierr)
                wild_flag = iprod2

                call MPI_bcast(end_grid, 4, MPI_integer, current_prc, COMM_world, ierr)
                call MPI_bcast(corr_flag, 2, MPI_integer, current_prc, COMM_world, ierr)

                stop_flag = 0
                !---- If lightning path reaches end_grid stop
                exit loop_path
             endif

             !---- determine wether process change occurs or not
             iprod1 = comm_flag
             call MPI_AllReduce(iprod1, iprod2, 1, MPI_integer, MPI_MAX, COMM_world, ierr)
             comm_flag = iprod2

             !--- process change occurs
             if( comm_flag == 1 ) then
                call MPI_bcast(comm_target, 1, MPI_integer, current_prc, COMM_world, ierr)
                !--- this part should be changed by 1 to 1 communication
                call MPI_bcast(k1, 1, MPI_integer, current_prc, COMM_world, ierr)
                call MPI_bcast(i1, 1, MPI_integer, current_prc, COMM_world, ierr)
                call MPI_bcast(j1, 1, MPI_integer, current_prc, COMM_world, ierr)
                call MPI_bcast(corr_flag, 2, MPI_integer, current_prc, COMM_world, ierr)
                call MPI_bcast(sign_flash, 1, MPI_integer, current_prc, COMM_world, ierr)
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
       dqrho_pls(:,:,:) = 0.0_RP
       dqrho_mns(:,:,:) = 0.0_RP
       do k = KS, KE
       do i = IS, IE
       do j = JS, JE

          !---- whether the grid is on lightning path or not
          if( L_path(k,i,j) /= 0.0_RP ) then
             if( abs( QCRG(k,i,j) ) > qrho_chan ) then
                if ( QCRG(k,i,j) >= 0.0_RP ) then
                   Npls = Npls + 1
                   dqrho_pls(k,i,j) = fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                   sumdqrho_pls = sumdqrho_pls &
                                + fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                else
                   Nmns = Nmns + 1
                   dqrho_mns(k,i,j) = fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                   sumdqrho_mns = sumdqrho_mns &
                                + fp * ( abs( QCRG(k,i,j) )-qrho_neut )
                endif
             endif
          endif

       enddo
       enddo
       enddo

       iprod1 =  Npls
       call MPI_AllReduce(iprod1, buf, 1, MPI_integer, MPI_SUM, COMM_world, ierr)
       Npls = buf

       iprod1 =  Nmns
       call MPI_AllReduce(iprod1, buf, 1, MPI_integer, MPI_SUM, COMM_world, ierr)
       Nmns = buf

       Ntot = Npls + Nmns

       if ( Ntot > 0 ) exit

    end do

    if ( count_path > NUTR_ITMAX ) then
       LOG_INFO("ATMOS_PHY_LT_Efield",*) "Reach limit iteration for searching flash path", count_path, Npls, Nmns, current_prc
       d_QCRG(:,:,:) = 0.0_RP
       L_path(:,:,:) = 0.0_RP
       return
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

    d_QCRG(:,:,:) = 0.0_RP
    rprod1 = 0.0_RP
    rprod2 = 0.0_RP
    !--- pseud-2D experiment for x-direction
    if( PRC_NUM_X == 1 .and. IMAX == 2 ) then

      do k = KS, KE
      do j = JS, JE
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

      do k = KS, KE
      do i = IS, IE
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

         do pm = 1, 3
           if( flg_num_end(k,i,JS,pm) == 1 ) then
               NUM_end(k,i,JE,pm) = NUM_end(k,i,JE,pm) + 1.0_RP
           endif
           if( flg_num_end(k,i,JE,pm) == 1 ) then
               NUM_end(k,i,JS,pm) = NUM_end(k,i,JS,pm) + 1.0_RP
           endif
         enddo

      enddo
      enddo

    else

      do k = KS, KE
      do i = IS, IE
      do j = JS, JE
         if( L_path(k,i,j) /= 0.0_RP .and. dqrho_pls(k,i,j) /= 0.0_RP ) then
           d_QCRG(k,i,j) = - dqrho_pls(k,i,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
         elseif( L_path(k,i,j) /= 0.0_RP .and. dqrho_mns(k,i,j) /= 0.0_RP ) then
           d_QCRG(k,i,j) =   dqrho_mns(k,i,j) - dqrho_cor  !--- Eq.(3) MacGorman et al. (2001)
         endif
      enddo
      enddo
      enddo

    endif

    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
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
       d_QCRG      ) ! [OUT]
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
       COMM_world, &
       COMM_vars8, &
       COMM_wait, &
       COMM_bcast
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

    real(RP), parameter :: q_thre = 0.1_RP ! threshold of discharge zone (Fierro et al. 2013) [nC/m3]

    integer  :: B(IA,JA)   !--- B in Fierro et al. 2013
    integer  :: C(IA,JA)   !--- flg for cylinder (in cylinder -> C=1)
    real(RP) :: Q_d, Fpls, Fmns
    real(RP) :: Spls, Smns
    real(RP) :: Spls_g, Smns_g
    real(RP) :: distance
    real(RP) :: Edif(KS:KE), abs_qcrg_cyl(KA), abs_qcrg_max

    integer, allocatable :: proc_num(:), proc_numg(:)
    real(RP),allocatable :: E_exce_x(:), E_exce_x_g(:)  !--- x point of column in which |E|>Eint is included [m] (_g means global attribute)
    real(RP),allocatable :: E_exce_y(:), E_exce_y_g(:)  !--- y point of column in which |E|>Eint is included [m] (_g means global attribute)
    real(RP) :: exce_grid(2,KIJMAX)
    integer  :: count1(PRC_nprocs)
    integer  :: countindx(PRC_nprocs+1)
    integer  :: num_own                               !--- number of column whose |E| > Eint-dEint for each process
    integer  :: num_total                             !--- total number of column whose |E| > Eint-dEint
    integer  :: k, i, j, ipp, iq, ierr

    call PROF_rapstart('LT_neut_F2013', 1)

    NUM_end(:,:,:,:) = 0.0_RP
    fls_int_p(:,:,:) = 0.0_RP
    B(:,:) = 0
    C(:,:) = 0
    !--- search grid whose Efield is over threshold of flash inititation
    num_own = 0
    do j = JS, JE
    do i = IS, IE
       Edif(:) = 0.0_RP
       do k = KS, KE
          Edif(k) = Efield(k,i,j,I_lt_abs) - ( Eint-delEint )
       enddo
       if( maxval(Edif,1) > 0.0_RP ) then
         num_own = num_own + 1
         exce_grid(1,num_own) = CX(i)
         exce_grid(2,num_own) = CY(j)
         do k = KS, KE
           if( Edif(k) > 0.0_RP ) then
             fls_int_p(k,i,j) = 1.0_RP
           endif
         enddo
       endif
    enddo
    enddo

    !**** proc_num(0~) -> process number of each grid with |E|> E_threthold (local)
!    allocate(proc_num(own_prc_total))
    allocate(E_exce_x(num_own))
    allocate(E_exce_y(num_own))

!    proc_num(:) = PRC_myrank
    E_exce_x(1:num_own) = exce_grid(1,1:num_own)
    E_exce_y(1:num_own) = exce_grid(2,1:num_own)

    call MPI_AllGather( num_own, 1, MPI_integer, &
                        count1, 1, MPI_integer, &
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


    C(:,:) = 0.0_RP
    do ipp = 1, num_total

      do i = IS, IE
      do j = JS, JE
           distance = sqrt( ( CX(i)-E_exce_x_g(ipp) )**2 &
                          + ( CY(j)-E_exce_y_g(ipp) )**2 )
           C(i,j) = C(i,j) + int(0.5_RP + sign( 0.5_RP,R_neut-distance )) !--- when distance > R_neut C=0, and distance < R_neut C=1
      enddo
      enddo

    enddo

    do i = IS, IE
    do j = JS, JE
      C(i,j) = min( C(i,j),1 )
    enddo
    enddo

    do i = IS, IE
    do j = JS, JE
         abs_qcrg_cyl(:) = 0.0_RP
         abs_qcrg_cyl(KS:KE) = abs( QCRG(KS:KE,i,j) )
         abs_qcrg_max = maxval( abs_qcrg_cyl,1 )
         B(i,j) = int( 0.5_RP + sign( 0.5_RP,abs_qcrg_max-q_thre )*real( C(i,j) ) )
    enddo
    enddo

    Spls = 0.0_RP
    Smns = 0.0_RP
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       if( QCRG(k,i,j) > 0.0_RP .and. B(i,j) == 1 ) then
         Spls = Spls + abs( QCRG(k,i,j) )
       elseif( QCRG(k,i,j) < 0.0_RP .and. B(i,j) == 1) then
         Smns = Smns + abs( QCRG(k,i,j) )
       endif
    enddo
    enddo
    enddo

    call MPI_Allreduce( Spls, Spls_g, 1, COMM_datatype, &
                        MPI_SUM, COMM_world, ierr       )
    call MPI_Allreduce( Smns, Smns_g, 1, COMM_datatype, &
                        MPI_SUM, COMM_world, ierr       )

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
    d_QCRG(:,:,:) = 0.0_RP
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
       if( QCRG(k,i,j) > q_thre ) then
         d_QCRG(k,i,j) = - Fpls * ( abs( QCRG(k,i,j) ) - q_thre ) * real( B(i,j) )
         LT_path(k,i,j) = LT_path(k,i,j) + 1.0_RP * real( B(i,j) )
       elseif( QCRG(k,i,j) < -q_thre ) then
         d_QCRG(k,i,j) = + Fmns * ( abs( QCRG(k,i,j) ) - q_thre ) * real( B(i,j) )
         LT_path(k,i,j) = LT_path(k,i,j) - 1.0_RP * real( B(i,j) )
       endif
    enddo
    enddo
    enddo

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

    !--- search grid whose Efield is over threshold of flash inititation
    own_prc_total = 0
    if( flg_eint_hgt == 1.0_RP ) then
      do k = KS, KE
      do j = JS, JE
      do i = IS, IE
         Eint_hgt(k,i,j) = min( Eint*DENS(k,i,j)/rho0,Eint )*flg_eint_hgt &
                         + ( Eint-delEint )*( 1.0_RP-flg_eint_hgt )
         if( flg_eint_hgt == 1.0_RP .and. Eint_hgt(k,i,j) < Estop ) then
            Eint_hgt(k,i,j) = LARGE_NUM !--- ignore when Eint < Estop
         endif
         E_det = Efield(k,i,j) - Eint_hgt(k,i,j)
         if( E_det > 0.0_RP ) then
           own_prc_total = own_prc_total + 1
         endif
      enddo
      enddo
      enddo
    else
      do k = KS, KE
      do j = JS, JE
      do i = IS, IE
         E_det = Efield(k,i,j) - ( Eint-delEint )
         if( E_det > 0.0_RP ) then
           own_prc_total = own_prc_total + 1
         endif
      enddo
      enddo
      enddo
    endif

    !--- Add number of grids, whose |E| are over threshold of flash initiaion, for all process
    iprod1 = own_prc_total
    call MPI_AllReduce(iprod1, buf, 1, MPI_integer, MPI_SUM, COMM_world, ierr)
    rprod1 = maxval(Efield(:,:,:))
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
    integer  :: grid(2)
    real(RP) :: cwc
    real(RP) :: diffx(nxlut_lt), diffy(nylut_lt)

    dqcrg(:,:,:) = 0.0_RP
    beta_crg(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      if( TEMP(k,i,j) <= T00 .and. TEMP(k,i,j) >= tcrglimit ) then
        cwc = 0.0_RP
        do iq = 1, NLIQ
          cwc = cwc + QLIQ(k,i,j,iq) * DENS(k,i,j) * 1.0E+3_RP ![g/m3]
        enddo
        do pp = 1, nxlut_lt
           diffx( pp ) = abs( grid_lut_t(pp)-TEMP(k,i,j) )
        enddo
        grid(1) = minloc( diffx,1 )
        do qq = 1, nylut_lt
           diffy( qq ) = abs( grid_lut_l(qq)-cwc )
        enddo
        grid(2) = minloc( diffy,1 )
        dqcrg( k,i,j ) = dq_chrg( grid(1), grid(2) ) &
                       *( 0.5_RP + sign( 0.5_RP,cwc-1.0E-2_RP ) ) !--- no charge separation when cwc < 0.01 [g/m3]
      endif
      if( TEMP(k,i,j) >= -30.0_RP+T00 ) then
        beta_crg( k,i,j ) = 1.0_RP
      elseif( TEMP(k,i,j) < -30.0_RP+T00 .and. TEMP(k,i,j) >= -43.0_RP+T00 ) then
        beta_crg( k,i,j ) = 1.0_RP - ( ( TEMP(k,i,j)-T00+30.0_RP )/13.0_RP )**2
      elseif( TEMP(k,i,j) < -43.0_RP+T00 ) then
        beta_crg( k,i,j ) = 0.0_RP
      endif
    enddo
    enddo
    enddo

    return

 end subroutine ATMOS_PHY_LT_sato2019_select_dQCRG_from_LUT
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_lt_sato2019
