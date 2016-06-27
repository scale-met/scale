module scale_atmos_phy_cp
  !------------------------------------------------------------------------------
  !
  !+++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  private

  abstract interface
     subroutine cp( &
          ![IN]
          DENS,      &
          MOMZ,      &
          MOMX,      &
          MOMY,      &
          RHOT,      &
          QTRC,      &
          ![OUT]
          DENS_t_CP, &
          MOMZ_t_CP, &
          MOMX_t_CP, &
          MOMY_t_CP, &
          RHOT_t_CP, &
          RHOQ_t_CP, &
          MFLX_cloudbase, &
          SFLX_convrain,  & ! [INOUT]
          cloudtop,       & ! [INOUT]
          cloudbase,      & ! [INOUT]
          cldfrac_dp,     & ! [INOUT]
          cldfrac_sh,     & ! [INOUT]
          kf_nca,         & ! [INOUT]
          kf_w0avg )        ! [INOUT]
       use scale_precision
       use scale_stdio
       use scale_prof
       use scale_grid_index
       use scale_tracer
       real(RP), intent(in) :: DENS(KA,IA,JA)
       real(RP), intent(in) :: MOMX(KA,IA,JA)
       real(RP), intent(in) :: MOMY(KA,IA,JA)
       real(RP), intent(in) :: MOMZ(KA,IA,JA)
       real(RP), intent(in) :: RHOT(KA,IA,JA)
       real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    ! [out]
       real(RP), intent(inout) :: DENS_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMX_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMY_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: RHOT_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: RHOQ_t_CP(KA,IA,JA,QA)
       real(RP), intent(inout) :: MFLX_cloudbase(IA,JA) ! dammy
       real(RP), intent(inout) :: SFLX_convrain(IA,JA)  ! convective rain
       real(RP), intent(inout) :: cloudtop(IA,JA)       ! cloud top height [m]
       real(RP), intent(inout) :: cloudbase(IA,JA)      ! cloud base height [m]
       real(RP), intent(inout) :: cldfrac_dp(KA,IA,JA)  ! cloud fraction (deep convection)
       real(RP), intent(inout) :: cldfrac_sh(KA,IA,JA)  ! cloud fraction (shallow convection)
       real(RP), intent(inout) :: kf_nca(IA,JA)         ! [step]
       real(RP), intent(inout) :: kf_W0avg(KA,IA,JA)    ! [m/s]
     end subroutine cp
  end interface
  
     procedure(cp), pointer :: ATMOS_PHY_CP => NULL()
     public :: ATMOS_PHY_CP
     public :: ATMOS_PHY_CP_setup
     
   contains
     !------------------------------------------------------------------------------
     !> Setup Cumulus parameterization
     !------------------------------------------------------------------------------
     subroutine ATMOS_PHY_CP_setup( CP_TYPE)
       use scale_process, only: &
            PRC_MPIstop
       use scale_atmos_phy_cp_kf, only: &
            ATMOS_PHY_CP_kf_setup, &
            ATMOS_PHY_CP_kf
       implicit none
       character(len=*),intent(in) :: CP_TYPE
       !------------------------------------------------------------------------------

       select case (CP_TYPE)
       case('KF')
          call ATMOS_PHY_CP_kf_setup( CP_TYPE )
          ATMOS_PHY_CP => ATMOS_PHY_CP_kf
       end select

       return
     end subroutine ATMOS_PHY_CP_setup
          
   end module scale_atmos_phy_cp
