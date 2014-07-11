!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Coversion between atmospheric fluid variables
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-11 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_fluidconv

  implicit none
  private
  public :: sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod

contains
  !-----------------------------------------------------------------------------
  subroutine sdm_rho_qtrc2rhod(dens,qtrc,rhod)
    use scale_precision
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_tracer, only: &
         QAD => QA, QQS, QQE
    ! Input variables
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: qtrc(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
    ! Output variables
    real(RP), intent(out) :: rhod(KA,IA,JA)  ! Density of dry air [kg/m3]

    real(RP) :: qdry(KA,IA,JA)
    integer :: i, j, k,iq  ! index
    !---------------------------------------------------------------------

    qdry(:,:,:) = 1.0_RP
    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
      do iq = QQS, QQE
        qdry(k,i,j) = qdry(k,i,j) - qtrc(k,i,j,iq)
      enddo
      rhod(k,i,j) = dens(k,i,j)*qdry(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine sdm_rho_qtrc2rhod
  !-----------------------------------------------------------------------------
  subroutine sdm_rho_mom2uvw(dens,momx,momy,momz,u,v,w)
    use scale_precision
    use scale_grid_index, only: &
         IA,JA,KA
    ! Input variables
    real(RP), intent(in)  :: dens(KA,IA,JA) ! Density [kg/m3]
    real(RP), intent(in)  :: momx(KA,IA,JA) ! Momentum [kg/s/m2]
    real(RP), intent(in)  :: momy(KA,IA,JA) 
    real(RP), intent(in)  :: momz(KA,IA,JA) 
    ! Output variables
    real(RP), intent(out) :: u(KA,IA,JA)  ! wind velocity [m/s]
    real(RP), intent(out) :: v(KA,IA,JA)  ! 
    real(RP), intent(out) :: w(KA,IA,JA)  ! 

    integer :: i, j, k  ! index
    !---------------------------------------------------------------------

    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
       u(k,i,j) = 2.0_RP * momx(k,i,j) / ( dens(k,i,j)+dens(k,i+1,j) )
       v(k,i,j) = 2.0_RP * momy(k,i,j) / ( dens(k,i,j)+dens(k,i,j+1) )
       w(k,i,j) = 2.0_RP * momz(k,i,j) / ( dens(k,i,j)+dens(k+1,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine sdm_rho_mom2uvw
  !-----------------------------------------------------------------------------
  subroutine sdm_rho_rhot2pt(dens,rhot,pt)
    use scale_precision
    use scale_grid_index, only: &
         IA,JA,KA
    ! Input variables
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: rhot(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    ! Output variables
    real(RP), intent(out) :: pt(KA,IA,JA)         !  Potential temperature [K]
    
    integer :: i, j, k  ! index
    !---------------------------------------------------------------------

    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
       pt(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine sdm_rho_rhot2pt
  !-----------------------------------------------------------------------------
  subroutine sdm_rhot_qtrc2p_t(rhot,qtrc,dens,p,t)
    use scale_precision
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_const, only: &
         P00 => CONST_PRE00, &
         Rdry => CONST_Rdry, &
         Rvap => CONST_Rvap, &
         CPdry => CONST_CPdry
    use scale_atmos_thermodyn, only: &
         CPw => AQ_CP
    use scale_tracer, only: &
         QAD => QA, QQS, QQE, I_QV
    ! Input variables
    real(RP), intent(in) :: rhot(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(in) :: qtrc(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    ! Output variables
    real(RP), intent(out) :: p(KA,IA,JA)          !  Pressure [Pa]
    real(RP), intent(out) :: t(KA,IA,JA)          !  Temperature [K]

    real(RP) :: qdry(KA,IA,JA)
    real(RP) :: rtot(KA,IA,JA)
    real(RP) :: cptot(KA,IA,JA)
    real(RP) :: cpovcv(KA,IA,JA)
    integer :: i, j, k, iq  ! index
    !---------------------------------------------------------------------

    qdry(:,:,:) = 1.0_RP
    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
      do iq = QQS, QQE
        qdry(k,i,j) = qdry(k,i,j) - qtrc(k,i,j,iq)
      enddo
      rtot (k,i,j) = Rdry * qdry(k,i,j) + Rvap * qtrc(k,i,j,I_QV)
      cptot(k,i,j) = CPdry * qdry(k,i,j)
      do iq = QQS, QQE
        cptot(k,i,j) = cptot(k,i,j) + qtrc(k,i,j,iq) * CPw(iq)
      enddo
      cpovcv(k,i,j) = cptot(k,i,j) / ( cptot(k,i,j) - rtot(k,i,j) )

      p(k,i,j) = P00 * ( rhot(k,i,j) * rtot(k,i,j) / P00 )**cpovcv(k,i,j)
      t(k,i,j) = (rhot(k,i,j)/dens(k,i,j)) * (p(k,i,j)/P00)**(rtot(k,i,j)/cptot(k,i,j))
    enddo
    enddo
    enddo

    return
  end subroutine sdm_rhot_qtrc2p_t
end module m_sdm_fluidconv
