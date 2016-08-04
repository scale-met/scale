!-------------------------------------------------------------------------------
!> Module thermodyanics
!!
!! @par Description
!!          This module calculates thermodyanical variables
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_thrmdyn
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID,      &
     kdim => ADM_kall, &
     kmin => ADM_kmin, &
     kmax => ADM_kmax
  use scale_const, only: &
     Rdry  => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     Rvap  => CONST_Rvap,  &
     PRE00 => CONST_PRE00, &
     TEM00 => CONST_TEM00, &
     PSAT0 => CONST_PSAT0, &
     EPSV  => CONST_EPSvap
  use mod_runconf, only: &
     nqmax => TRC_VMAX, &
     NQW_STR,           &
     NQW_END,           &
     I_QV,              &
     I_QC,              &
     I_QR,              &
     I_QI,              &
     I_QS,              &
     I_QG,              &
     LHV,               &
     LHF,               &
     CVW,               &
     CPW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: THRMDYN_qd
  public :: THRMDYN_cv
  public :: THRMDYN_cp
  public :: THRMDYN_rho
  public :: THRMDYN_pre
  public :: THRMDYN_ein
  public :: THRMDYN_tem
  public :: THRMDYN_th
  public :: THRMDYN_eth
  public :: THRMDYN_ent
  public :: THRMDYN_rhoein
  public :: THRMDYN_tempre

  interface THRMDYN_qd
     module procedure THRMDYN_qd_ijk
     module procedure THRMDYN_qd_ijkl
  end interface THRMDYN_qd

  interface THRMDYN_cv
     module procedure THRMDYN_cv_ijk
  end interface THRMDYN_cv

  interface THRMDYN_cp
     module procedure THRMDYN_cp_ijk
  end interface THRMDYN_cp

  interface THRMDYN_rho
     module procedure THRMDYN_rho_ijk
     module procedure THRMDYN_rho_ijkl
  end interface THRMDYN_rho

  interface THRMDYN_pre
     module procedure THRMDYN_pre_ijk
  end interface THRMDYN_pre

  interface THRMDYN_ein
     module procedure THRMDYN_ein_ijk
  end interface THRMDYN_ein

  interface THRMDYN_tem
     module procedure THRMDYN_tem_ijk
  end interface THRMDYN_tem

  interface THRMDYN_th
     module procedure THRMDYN_th_ijk
     module procedure THRMDYN_th_ijkl
  end interface THRMDYN_th

  interface THRMDYN_eth
     module procedure THRMDYN_eth_ijk
     module procedure THRMDYN_eth_ijkl
  end interface THRMDYN_eth

  interface THRMDYN_ent
     module procedure THRMDYN_ent_ijk
  end interface THRMDYN_ent

  interface THRMDYN_rhoein
     module procedure THRMDYN_rhoein_ijkl
  end interface THRMDYN_rhoein

  interface THRMDYN_tempre
     module procedure THRMDYN_tempre_ijkl
  end interface THRMDYN_tempre

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijk( &
       ijdim, &
       kdim,  &
       q,     &
       qd     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(qd) pcopyin(q) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k) = 1.0_RP

       !$acc loop seq
       do nq = NQW_STR,NQW_END
          qd(ij,k) = qd(ij,k) - q(ij,k,nq)
       enddo
       !$acc end loop

    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_qd_ijk

  !-----------------------------------------------------------------------------
  !> calculate dry air
  subroutine THRMDYN_qd_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       q,     &
       qd     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(RP), intent(in)  :: q (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: qd(ijdim,kdim,ldim)       ! dry air mass concentration [kg/kg]

    integer :: ij, k, l,nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(qd) pcopyin(q) async(0)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       qd(ij,k,l) = 1.0_RP

       !$acc loop seq
       do nq = NQW_STR,NQW_END
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       !$acc end loop

    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_qd_ijkl

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cv_ijk( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cv     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: cv(ijdim,kdim)       ! specific heat [J/kg/K]

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(cv) pcopyin(qd,q,CVW) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       !$acc end loop

    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_cv_ijk

  !-----------------------------------------------------------------------------
  !> calculate specific heat
  subroutine THRMDYN_cp_ijk( &
       ijdim, &
       kdim,  &
       qd,    &
       q,     &
       cp     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: qd(ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: cp(ijdim,kdim)       ! specific heat [J/kg/K]

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(cp) pcopyin(qd,q,cpW) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       cp(ij,k) = qd(ij,k) * CPdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cp(ij,k) = cp(ij,k) + q(ij,k,nq) * CPW(nq)
       enddo
       !$acc end loop

    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_cp_ijk

  !-----------------------------------------------------------------------------
  !> calculate density
  subroutine THRMDYN_rho_ijk( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       rho    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(RP), intent(in)  :: pre(ijdim,kdim)       ! pressure    [Pa]
    real(RP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: rho(ijdim,kdim)       ! density     [kg/m3]

    integer :: ij, k
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(rho) pcopyin(pre,tem,qd,q) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = pre(ij,k) / tem(ij,k) / ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_rho_ijk

  !-----------------------------------------------------------------------------
  !> calculate density
  subroutine THRMDYN_rho_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       rho    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(RP), intent(in)  :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(RP), intent(in)  :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]
    real(RP), intent(in)  :: qd (ijdim,kdim,ldim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(rho) pcopyin(pre,tem,qd,q) async(0)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k,l) = pre(ij,k,l) / tem(ij,k,l) / ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap )
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_rho_ijkl

  !-----------------------------------------------------------------------------
  !> calculate pressure
  subroutine THRMDYN_pre_ijk( &
       ijdim, &
       kdim,  &
       rho,   &
       tem,   &
       qd,    &
       q,     &
       pre    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: rho(ijdim,kdim)       ! density     [kg/m3]
    real(RP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(RP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: pre(ijdim,kdim)       ! pressure    [Pa]

    integer :: ij, k
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(pre) pcopyin(rho,tem,qd,q) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       pre(ij,k) = rho(ij,k) * tem(ij,k) * ( qd(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_pre_ijk

  !-----------------------------------------------------------------------------
  !> calculate internal energy
  subroutine THRMDYN_ein_ijk( &
       ijdim, &
       kdim,  &
       tem,   &
       qd,    &
       q,     &
       ein    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: tem(ijdim,kdim)       ! temperature [K]
    real(RP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: ein(ijdim,kdim)       ! internal energy [J]

    real(RP) :: cv(ijdim,kdim)

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(ein,cv) pcopyin(tem,qd,q,CVW) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       !$acc end loop

       ein(ij,k) = tem(ij,k) * cv(ij,k)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_ein_ijk

  !-----------------------------------------------------------------------------
  !> calculate temperature
  subroutine THRMDYN_tem_ijk( &
       ijdim, &
       kdim,  &
       ein,   &
       qd,    &
       q,     &
       tem    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: ein(ijdim,kdim)       ! internal energy [J]
    real(RP), intent(in)  :: qd (ijdim,kdim)       ! dry air mass concentration [kg/kg]
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: tem(ijdim,kdim)       ! temperature [K]

    real(RP) :: cv(ijdim,kdim)

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(tem,cv) pcopyin(ein,qd,q,CVW) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k) = qd(ij,k) * CVdry

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k) = cv(ij,k) + q(ij,k,nq) * CVW(nq)
       enddo
       !$acc end loop

       tem(ij,k) = ein(ij,k) / cv(ij,k)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_tem_ijk

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijk( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       th     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: tem(ijdim,kdim) ! temperature [K]
    real(RP), intent(in)  :: pre(ijdim,kdim) ! pressure    [Pa]
    real(RP), intent(out) :: th (ijdim,kdim) ! potential temperature [K]

    real(RP) :: pre0_kappa, kappa

    integer :: ij, k
    !---------------------------------------------------------------------------

    kappa = Rdry /Cpdry
    pre0_kappa = PRE00**kappa

    !$acc kernels pcopy(th) pcopyin(tem,pre) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k) = tem(ij,k) + pre(ij,k)**kappa * pre0_kappa
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_th_ijk

  !-----------------------------------------------------------------------------
  !> calculate potential temperature
  subroutine THRMDYN_th_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       th     )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(RP), intent(in)  :: tem(ijdim,kdim,ldim) ! temperature [K]
    real(RP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure    [Pa]
    real(RP), intent(out) :: th (ijdim,kdim,ldim) ! potential temperature [K]

    real(RP) :: pre0_kappa, kappa

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    kappa = Rdry /Cpdry
    pre0_kappa = PRE00**kappa

    !$acc kernels pcopy(th) pcopyin(tem,pre) async(0)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       th(ij,k,l) = tem(ij,k,l) + pre(ij,k,l)**kappa * pre0_kappa
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_th_ijkl

  !-----------------------------------------------------------------------------
  !> calculate enthalpy
  subroutine THRMDYN_eth_ijk( &
       ijdim, &
       kdim,  &
       ein,   &
       pre,   &
       rho,   &
       eth    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: ein(ijdim,kdim) ! internal energy [J]
    real(RP), intent(in)  :: pre(ijdim,kdim) ! pressure    [Pa]
    real(RP), intent(in)  :: rho(ijdim,kdim) ! density     [kg/m3]
    real(RP), intent(out) :: eth(ijdim,kdim) ! enthalpy

    integer :: ij, k
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(eth) pcopyin(ein,pre,rho) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       eth(ij,k) = ein(ij,k) + pre(ij,k) / rho(ij,k)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_eth_ijk

  !-----------------------------------------------------------------------------
  !> calculate enthalpy
  subroutine THRMDYN_eth_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       pre,   &
       rho,   &
       eth    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(RP), intent(in)  :: ein(ijdim,kdim,ldim) ! internal energy [J]
    real(RP), intent(in)  :: pre(ijdim,kdim,ldim) ! pressure    [Pa]
    real(RP), intent(in)  :: rho(ijdim,kdim,ldim) ! density     [kg/m3]
    real(RP), intent(out) :: eth(ijdim,kdim,ldim) ! enthalpy

    integer :: ij, k, l
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(eth) pcopyin(ein,pre,rho) async(0)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       eth(ij,k,l) = ein(ij,k,l) + pre(ij,k,l) / rho(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_eth_ijkl

  !-----------------------------------------------------------------------------
  !> calculate entropy
  subroutine THRMDYN_ent_ijk( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       q,     &
       ent    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(in)  :: pre(ijdim,kdim)
    real(RP), intent(in)  :: qd (ijdim,kdim)
    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax)
    real(RP), intent(out) :: ent(ijdim,kdim)

    real(RP) :: Pdry
    real(RP) :: Pvap
    real(RP) :: LH(nqmax)

    real(RP), parameter :: EPS = 1.E-10_RP

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    do nq = NQW_STR, NQW_END
       if ( nq == I_QV ) then
          LH(nq) =  LHV / TEM00
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          LH(nq) = -LHF / TEM00
       else
          LH(nq) = 0.0_RP
       endif
    enddo

    !$acc kernels pcopy(ent) pcopyin(tem,pre,qd,q) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim
       Pdry = max( pre(ij,k) * EPSV*qd(ij,k) / ( EPSV*qd(ij,k) + q(ij,k,I_QV) ), EPS )
       Pvap = max( pre(ij,k) * q(ij,k,I_QV)  / ( EPSV*qd(ij,k) + q(ij,k,I_QV) ), EPS )

       ent(ij,k) = qd(ij,k)      * CPdry * log( tem(ij,k)/TEM00 ) &
                 - qd(ij,k)      * Rdry  * log( Pdry     /PRE00 ) &
                 - q (ij,k,I_QV) * Rvap  * log( Pvap     /PSAT0 )
    enddo
    enddo
    !$acc end kernels

    !$acc kernels pcopy(ent) pcopyin(tem,q,CVW,LH) async(0)
    do k  = 1, kdim
    do ij = 1, ijdim

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          ent(ij,k) = ent(ij,k) + q(ij,k,nq) * CVW(nq) * log( tem(ij,k)/TEM00 ) &
                                + q(ij,k,nq) * LH (nq) / TEM00
       enddo
       !$acc end loop

    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_ent_ijk

  !-----------------------------------------------------------------------------
  !> calculate density & internal energy
  subroutine THRMDYN_rhoein_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       tem,   &
       pre,   &
       q,     &
       rho,   &
       ein    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(RP), intent(in)  :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(RP), intent(in)  :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]
    real(RP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(RP), intent(out) :: ein(ijdim,kdim,ldim)       ! internal energy [J]

    real(RP) :: cv(ijdim,kdim,ldim)
    real(RP) :: qd(ijdim,kdim,ldim)

    integer :: ij, k, l, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(rho,ein,cv,qd) pcopyin(pre,tem,q,CVW) async(0)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_RP
       qd(ij,k,l) = 1.0_RP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       !$acc end loop

       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       rho(ij,k,l) = pre(ij,k,l) / tem(ij,k,l) / ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap )
       ein(ij,k,l) = tem(ij,k,l) * cv(ij,k,l)
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_rhoein_ijkl

  !-----------------------------------------------------------------------------
  !> calculate temperature & pressure
  subroutine THRMDYN_tempre_ijkl( &
       ijdim, &
       kdim,  &
       ldim,  &
       ein,   &
       rho,   &
       q,     &
       tem,   &
       pre    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: ldim
    real(RP), intent(in)  :: ein(ijdim,kdim,ldim)       ! internal energy [J]
    real(RP), intent(in)  :: rho(ijdim,kdim,ldim)       ! density     [kg/m3]
    real(RP), intent(in)  :: q  (ijdim,kdim,ldim,nqmax) ! tracer  mass concentration [kg/kg]
    real(RP), intent(out) :: tem(ijdim,kdim,ldim)       ! temperature [K]
    real(RP), intent(out) :: pre(ijdim,kdim,ldim)       ! pressure    [Pa]

    real(RP) :: cv(ijdim,kdim,ldim)
    real(RP) :: qd(ijdim,kdim,ldim)

    integer :: ij, k, l, nq
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(tem,pre,cv,qd) pcopyin(ein,rho,q,CVW) async(0)
    do l  = 1, ldim
    do k  = 1, kdim
    do ij = 1, ijdim
       cv(ij,k,l) = 0.0_RP
       qd(ij,k,l) = 1.0_RP

       !$acc loop seq
       do nq = NQW_STR, NQW_END
          cv(ij,k,l) = cv(ij,k,l) + q(ij,k,l,nq) * CVW(nq)
          qd(ij,k,l) = qd(ij,k,l) - q(ij,k,l,nq)
       enddo
       !$acc end loop

       cv(ij,k,l) = cv(ij,k,l) + qd(ij,k,l) * CVdry

       tem(ij,k,l) = ein(ij,k,l) / cv(ij,k,l)
       pre(ij,k,l) = rho(ij,k,l) * tem(ij,k,l) * ( qd(ij,k,l)*Rdry + q(ij,k,l,I_QV)*Rvap )
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine THRMDYN_tempre_ijkl

end module mod_thrmdyn
