!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cumulus Parameterization
!!
!! @par Description
!!         Cumulus Convection parameterization
!!
!! @author Team SCALE
!!
!<
module scale_atmos_phy_cp
  !------------------------------------------------------------------------------
  !
  !+++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_tracer
  use scale_atmos_phy_mp
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  abstract interface
     subroutine cp( &
          KA, KS, KE,     &
          IA, IS, IE,     &
          JA, JS, JE,     &
          QA_MP, QS_MP, QE_MP, &
          DENS,           &
!          MOMZ,           &
!          MOMX,           &
!          MOMY,           &
          U,           &
          V,           &
          RHOT,           &
          QTRC,           &
          w0avg,          &
          DENS_t_CP,      &
          MOMZ_t_CP,      &
          MOMX_t_CP,      &
          MOMY_t_CP,      &
          RHOT_t_CP,      &
          RHOQ_t_CP,      &
          MFLX_cloudbase, &
          SFLX_convrain,  &
          cloudtop,       &
          cloudbase,      &
          cldfrac_dp,     &
          cldfrac_sh,     &
          kf_nca          )
       use scale_file_history, only: &
          FILE_HISTORY_in
       use scale_precision
       use scale_tracer
       use scale_const, only: &
          GRAV => CONST_GRAV, &
          R    => CONST_Rdry
       use scale_atmos_hydrometeor, only: &
          I_QV, &
          I_QC, &
          I_QR, &
          I_QI, &
          I_QS
       use scale_time , only :&
          KF_DTSEC => TIME_DTSEC_ATMOS_PHY_CP
       use scale_atmos_grid_cartesC_real, only: &
          FZ => ATMOS_GRID_CARTESC_REAL_FZ
       use scale_atmos_thermodyn, only: &
          THERMODYN_temp_pres   => ATMOS_THERMODYN_temp_pres,   &
          THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,        &
          THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E, &
          THERMODYN_qd          => ATMOS_THERMODYN_qd,          &
          THERMODYN_pott        => ATMOS_THERMODYN_pott
       use scale_atmos_saturation ,only :&

       SATURATION_psat_liq => ATMOS_SATURATION_psat_liq
       integer,  intent(in)    :: KA, KS, KE
       integer,  intent(in)    :: IA, IS, IE
       integer,  intent(in)    :: JA, JS, JE
       integer,  intent(in)    :: QA_MP, QS_MP, QE_MP

       real(RP), intent(in)    :: DENS(KA,IA,JA)
!       real(RP), intent(in)    :: MOMX(KA,IA,JA)
!       real(RP), intent(in)    :: MOMY(KA,IA,JA)
!       real(RP), intent(in)    :: MOMZ(KA,IA,JA)
       real(RP), intent(in)    :: U(KA,IA,JA)
       real(RP), intent(in)    :: V(KA,IA,JA)
       real(RP), intent(in)    :: RHOT(KA,IA,JA)
       real(RP), intent(in)    :: QTRC(KA,IA,JA,QS_MP:QE_MP)
       real(RP), intent(in)    :: w0avg(KA,IA,JA)
       real(RP), intent(inout) :: DENS_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMX_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMY_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: RHOT_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: RHOQ_t_CP(KA,IA,JA,QA_MP)
       real(RP), intent(inout) :: MFLX_cloudbase(IA,JA)
       real(RP), intent(inout) :: SFLX_convrain(IA,JA)
       real(RP), intent(inout) :: cloudtop(IA,JA)
       real(RP), intent(inout) :: cloudbase(IA,JA)
       real(RP), intent(inout) :: cldfrac_dp(KA,IA,JA)
       real(RP), intent(inout) :: cldfrac_sh(KA,IA,JA)
       real(RP), intent(inout) :: kf_nca(IA,JA)
     end subroutine cp
  end interface
  procedure(cp), pointer :: ATMOS_PHY_CP => NULL()
  public :: ATMOS_PHY_CP
  public :: ATMOS_PHY_CP_setup
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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cumulus parameterization
  !------------------------------------------------------------------------------
  subroutine ATMOS_PHY_CP_setup( &
      KA, KS, KE,   &
      IA, IS, IE,   &
      JA, JS, JE,   &
      CP_TYPE       )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_phy_cp_kf, only: &
         ATMOS_PHY_CP_kf_setup, &
         ATMOS_PHY_CP_kf
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    character(len=*), intent(in) :: CP_TYPE
    !------------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** => ', trim(CP_TYPE), ' is selected.'

    select case( CP_TYPE )
    case('OFF')
       ! do nothing
    case('KF')
       call ATMOS_PHY_CP_kf_setup( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                   CP_TYPE )
       ATMOS_PHY_CP => ATMOS_PHY_CP_kf

    case default
       write(*,*) 'xxx invalid Cumulus parameterization type(', trim(CP_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_CP_setup

end module scale_atmos_phy_cp
