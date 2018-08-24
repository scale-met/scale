!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User module for the GABLS2 experiment
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  integer,  private, parameter :: K_ini = 10
  real(RP), private, parameter :: z_ini (K_ini) = &
       (/    0.0_RP,  200.0_RP,  850.0_RP,  900.0_RP,  1000.0_RP, &
          2000.0_RP, 3500.0_RP, 4000.0_RP, 5000.0_RP, 40000.0_RP  /)

  real(RP), private, parameter :: pt_ini(K_ini) = &
       (/  288.0_RP,  286.0_RP,  286.0_RP,  288.0_RP,   292.0_RP, &
           300.0_RP,  310.0_RP,  312.0_RP,  316.0_RP,   800.0_RP  /)

  real(RP), private, parameter :: qv_ini(K_ini) = &
       (/ 0.0025_RP, 0.0025_RP, 0.0025_RP,  0.0025_RP,   0.0005_RP, &
          0.0003_RP, 0.0002_RP, 0.00015_RP, 0.0000_RP,   0.0000_RP  /)

  real(RP), private :: Ps   = 972E+2_RP
  real(RP), private :: Ug   =  3.0_RP
  real(RP), private :: Vg   = -9.0_RP
  real(RP), private :: Z0M  = 0.03_RP
  real(RP), private :: Z0H  = 0.003_RP
  real(RP), private :: ALB_SW = 0.20_RP
  real(RP), private :: ALB_LW = 0.01_RP
  real(RP), private :: LSD  = -0.005_RP ! Large-scale synoptic divergence
  logical,  private :: dry     = .false.

  real(RP), private, allocatable :: fluxf(:) ! large scale sinking (full level)
  real(RP), private, allocatable :: fluxh(:) ! large scale sinking (half level)
  real(RP), private, allocatable :: divf(:)  ! large scale divergence (full level)
  real(RP), private, allocatable :: divh(:)  ! large scale divergence (half level)

  integer,                private, parameter :: QA = 3
  character(len=H_SHORT), private            :: QNAME(QA)
  character(len=H_MID),   private            :: QDESC(QA)
  character(len=H_SHORT), private            :: QUNIT(QA)

  data QNAME / 'QV', &
               'QC', &
               'QI'  /

  data QDESC / 'Ratio of Water Vapor mass to total mass (Specific humidity)',   &
               'Ratio of Cloud Water mass to total mass', &
               'Ratio of Cloud Ice mass to total mass'   /

  data QUNIT / 'kg/kg', &
               'kg/kg', &
               'kg/kg'  /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer Setup
  subroutine USER_tracer_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_USER_qhyd2qtrc
    implicit none

    namelist / PARAM_USER / &
         Ps, &
         Ug, &
         Vg, &
         LSD, &
         dry

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_tracer_setup",*) 'Setup'
    LOG_INFO("USER_tracer_setup",*) '+++ GABLS2 SCM experiment'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_tracer_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_tracer_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    if ( .not. dry ) then
       call ATMOS_HYDROMETEOR_regist( 1, 1,                & ! (in)
                                      QNAME, QDESC, QUNIT, & ! (in)
                                      QS_MP                ) ! (out)
    end if

    QA_MP = 3
    QE_MP = QS_MP + 2

    ATMOS_PHY_MP_USER_qhyd2qtrc => USER_qhyd2qtrc

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ, &
       FZ => ATMOS_GRID_CARTESC_FZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ
    implicit none

    integer  :: k
    !---------------------------------------------------------------------------

    allocate( fluxf(KA), fluxh(KA) )
    allocate( divf(KA), divh(KA) )

    fluxh(KS-1) = 0.0_RP
    do k = KS, KE
       if ( CZ(k) > 1000.0_RP ) then
          fluxf(k) = LSD
       else
          fluxf(k) = LSD * CZ(k) * 1e-3_RP
       end if
       if ( FZ(k) > 1000.0_RP ) then
          fluxh(k) = LSD
       else
          fluxh(k) = LSD * FZ(k) * 1e-3_RP
       end if
    end do

    do k = KS, KE
       divf(k) = ( fluxh(k) - fluxh(k-1) ) * RCDZ(k)
    end do
    do k = KS, KE-1
       divh(k) = ( fluxf(k+1) - fluxf(k) ) * RFDZ(k)
    end do

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_atmos_hydrostatic, only: &
       buildrho => ATMOS_HYDROSTATIC_buildrho
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_SFC_TEMP,   &
       ATMOS_PHY_SF_SFC_albedo, &
       ATMOS_PHY_SF_SFC_Z0M,    &
       ATMOS_PHY_SF_SFC_Z0H,    &
       ATMOS_PHY_SF_SFC_Z0E
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ, &
       FZ => ATMOS_GRID_CARTESC_FZ
    implicit none

    real(RP) :: RHO (KA)
    real(RP) :: TEMP(KA)
    real(RP) :: PRES(KA)
    real(RP) :: PT  (KA)
    real(RP) :: QV  (KA)
    real(RP) :: QC  (KA)
    real(RP) :: Ts
    real(RP) :: QVs

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    call interporate( PT(:), pt_ini )
    if ( dry ) then
       QV(:) = 0.0_RP
       QVs = 0.0_RP
    else
       call interporate( QV(:), qv_ini )
       QVs = qv_ini(1)
    end if
    QC(:) = 0.0_RP

    call buildrho( KA, KS, KE, &
                   PT(:), QV(:), QC(:),         & ! (in)
                   Ps, pt_ini(1), QVs, 0.0_RP,  & ! (in)
                   CZ(:), FZ(:),                & ! (in)
                   RHO(:), TEMP(:), PRES(:), Ts ) ! (out)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = RHO(k) * Ug
       MOMY(k,i,j) = RHO(k) * Vg
       RHOT(k,i,j) = RHO(k) * PT(k)
    end do
    end do
    end do

    QTRC(:,:,:,:) = 0.0_RP
    if ( .not. dry ) then
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          QTRC(k,i,j,I_QV) = QV(k)
       end do
       end do
       end do
    end if

    call set_tg( NOWDAYSEC, Ts )

    ATMOS_PHY_SF_SFC_TEMP(:,:) = Ts

    ATMOS_PHY_SF_SFC_albedo(:,:,:,I_R_IR)  = ALB_LW
    ATMOS_PHY_SF_SFC_albedo(:,:,:,I_R_NIR) = ALB_SW
    ATMOS_PHY_SF_SFC_albedo(:,:,:,I_R_VIS) = ALB_SW

    ATMOS_PHY_SF_SFC_Z0M (:,:) = Z0M
    ATMOS_PHY_SF_SFC_Z0H (:,:) = Z0H
    ATMOS_PHY_SF_SFC_Z0E (:,:) = Z0H


    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine USER_calc_tendency
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC, &
       dt     => TIME_DTSEC
    use scale_atmos_grid_cartesC, only: &
       FDZ  => ATMOS_GRID_CARTESC_FDZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ
    use scale_atmos_dyn, only: &
       CORIOLIS
    use mod_atmos_vars, only: &
       DENS,    &
       MOMZ,    &
       U,       &
       V,       &
       POTT,    &
       QTRC,    &
       MOMZ_tp, &
       RHOU_tp, &
       RHOV_tp, &
       RHOT_tp, &
       RHOQ_tp
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_SFC_TEMP
    implicit none

    real(RP) :: val(KA)

    real(RP) :: Ts

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    val(KS-1) = 0.0_RP

    do j = JS, JE
    do i = IS, IE

       ! large scale divergence

       !! full level
       val(KS) = MOMZ(KS,i,j) / DENS(KS,i,j) * 0.5_RP
       do k = KS+1, KE
          val(k) = ( MOMZ(k,i,j) + MOMZ(k-1,i,j) ) / DENS(k,i,j) * 0.5_RP
       end do
       do k = KS, KE-1
          MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) &
               - ( fluxf(k+1)*val(k+1) - fluxf(k)*val(k) ) * RFDZ(k) &
               + MOMZ(k,i,j) * divh(k) / ( DENS(k,i,j)*FDZ(k+1) + DENS(k+1,i,j)*FDZ(k) ) * ( FDZ(k+1) + FDZ(k) )
       end do

       !! half level
       do k = KS, KE-1
          val(k) = ( U(k,i,j)*FDZ(k+1) + U(k+1,i,j)*FDZ(k) ) &
                 / ( FDZ(k+1) + FDZ(k) )
       end do
       val(KE) = U(KE,i,j)
       do k = KS, KE
          RHOU_tp(k,i,j) = RHOU_tp(k,i,j) &
               - ( fluxh(k)*val(k) - fluxh(k-1)*val(k-1) ) * RCDZ(k) &
               + U(k,i,j) * divf(k)
       end do

       !! half level
       do k = KS, KE-1
          val(k) = ( V(k,i,j)*FDZ(k+1) + V(k+1,i,j)*FDZ(k) ) &
                 / ( FDZ(k+1) + FDZ(k) )
       end do
       val(KE) = V(KE,i,j)
       do k = KS, KE
          RHOV_tp(k,i,j) = RHOV_tp(k,i,j) &
               - ( fluxh(k)*val(k) - fluxh(k-1)*val(k-1) ) * RCDZ(k) &
               + V(k,i,j) * divf(k)
       end do

       !! half level
       do k = KS, KE-1
          val(k) = ( POTT(k,i,j)*FDZ(k+1) + POTT(k+1,i,j)*FDZ(k) ) &
                 / ( FDZ(k+1) + FDZ(k) )
       end do
       val(KE) = POTT(KE,i,j)
       do k = KS, KE
          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) &
               - ( fluxh(k)*val(k) - fluxh(k-1)*val(k-1) ) * RCDZ(k) &
               + POTT(k,i,j) * divf(k)
       end do

       do iq = 1, QA
          !! half level
          do k = KS, KE-1
             val(k) = ( QTRC(k,i,j,iq)*FDZ(k+1) + QTRC(k+1,i,j,iq)*FDZ(k) ) &
                    / ( FDZ(k+1) + FDZ(k) )
          end do
          val(KE) = QTRC(KE,i,j,iq)
          do k = KS, KE
             RHOQ_tp(k,i,j,iq) = RHOQ_tp(k,i,j,iq) &
                  - ( fluxh(k)*val(k) - fluxh(k-1)*val(k-1) ) * RCDZ(k) &
                  + QTRC(k,i,j,iq) * divf(k)
          end do
       end do

       ! geostrophic forcing
       do k = KS, KE
          RHOU_tp(k,i,j) = RHOU_tp(k,i,j) - CORIOLIS(i,j) * Vg * DENS(k,i,j)
          RHOV_tp(k,i,j) = RHOV_tp(k,i,j) + CORIOLIS(i,j) * Ug * DENS(k,i,j)
       end do


    end do
    end do

    call set_tg( NOWDAYSEC, Ts )

    ATMOS_PHY_SF_SFC_TEMP(:,:) = Ts

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_update

  subroutine USER_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QV, QHYD, &
       QTRC,     &
       QNUM      )
    use scale_atmos_hydrometeor, only: &
         N_HYD, &
         I_QV, &
         I_HC, &
         I_HI
    use mod_atmos_phy_mp_vars, only: &
         QA_MP
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: QV   (KA,IA,JA)
    real(RP), intent(in) :: QHYD(KA,IA,JA,N_HYD)

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP)

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)

    QTRC(:,:,:,1) = QV(:,:,:)
    QTRC(:,:,:,2) = QHYD(:,:,:,I_HC)
    QTRC(:,:,:,3) = QHYD(:,:,:,I_HI)

    return
  end subroutine USER_qhyd2qtrc

  !-----------------------------------------------------------------------------
  subroutine interporate( d_out, d_in )
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ
    implicit none

    real(RP), intent(out) :: d_out(KA)
    real(RP), intent(in)  :: d_in(K_ini)

    integer :: k, kk
    !---------------------------------------------------------------------------

    do k = KS, KE

       if ( CZ(k) >= z_ini(K_ini) ) then
          d_out(k) = d_in(K_ini)
       elseif( CZ(k) < z_ini(1) ) then
          d_out(k) = d_in(1)
       else
          do kk = 1, K_ini-1
             if (       z_ini(kk)   <= CZ(k) &
                  .AND. z_ini(kk+1) >  CZ(k) ) then
                d_out(k) = d_in(kk) + ( d_in(kk+1)-d_in(kk) ) * ( CZ(k)-z_ini(kk) ) / ( z_ini(kk+1)-z_ini(kk) )
                exit
             end if
          end do
       end if

    end do

    d_out(   1:KS-1) = d_out(KS)
    d_out(KE+1:KA  ) = d_out(KE)

    return
  end subroutine interporate

  !-----------------------------------------------------------------------------
  subroutine set_tg( &
       current_sec, &
       Ts           )
    implicit none

    real(RP), intent(in)  :: current_sec
    real(RP), intent(out) :: Ts

    real(RP) :: time
    real(RP) :: Tg
    !---------------------------------------------------------------------------

    time = current_sec / 3600.0_RP + 16.0_RP

    if    ( time <= 17.4_RP ) then
       Tg = -10.0_RP - 25.0_RP * cos(time*0.22_RP + 0.2_RP)
    elseif( time <= 30.0_RP ) then
       Tg = -0.54_RP * time + 15.2_RP
    elseif( time <= 41.9_RP ) then
       Tg = -7.0_RP - 25.0_RP * cos(time*0.21_RP + 1.8_RP)
    elseif( time <= 53.3_RP ) then
       Tg = -0.37_RP * time + 18.0_RP
    elseif( time <= 65.6_RP ) then
       Tg = -4.0_RP - 25.0_RP * cos(time * 0.22_RP + 2.5_RP)
    else
       Tg = 4.4_RP
    end if

    Ts = Tg + 273.15_RP

    return
  end subroutine set_tg

end module mod_user
