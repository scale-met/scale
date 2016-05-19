!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Thermodynamics
!!
!! @par Description
!!          Thermodynamics module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)     [new] Import from NICAM
!! @li      2012-02-10 (H.Yashiro)   [mod] Reconstruction
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-12-22 (S.Nishizawa) [mod] Use thermodyn macro set
!!
!<
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
#include "inc_openmp.h"
module scale_atmos_thermodyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_const, only: &
     EPS   => CONST_EPS,   &
     Rdry  => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CPvap => CONST_CPvap, &
     CVdry => CONST_CVdry, &
     CL    => CONST_CL,    &
     CI    => CONST_CI,    &
     Rvap  => CONST_Rvap,  &
     LHV0  => CONST_LHV0,  &
     LHS0  => CONST_LHS0,  &
     LHF0  => CONST_LHF0,  &
     PSAT0 => CONST_PSAT0, &
     PRE00 => CONST_PRE00, &
     TEM00 => CONST_TEM00
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_THERMODYN_setup

  public :: ATMOS_THERMODYN_qd
  public :: ATMOS_THERMODYN_cv
  public :: ATMOS_THERMODYN_cp
  public :: ATMOS_THERMODYN_pott
  public :: ATMOS_THERMODYN_rhoe
  public :: ATMOS_THERMODYN_rhot
  public :: ATMOS_THERMODYN_temp_pres
  public :: ATMOS_THERMODYN_temp_pres_E
  public :: ATMOS_THERMODYN_templhv
  public :: ATMOS_THERMODYN_templhs
  public :: ATMOS_THERMODYN_templhf
  public :: ATMOS_THERMODYN_entr

  interface ATMOS_THERMODYN_qd
     module procedure ATMOS_THERMODYN_qd_0D
     module procedure ATMOS_THERMODYN_qd_3D
  end interface ATMOS_THERMODYN_qd

  interface ATMOS_THERMODYN_cp
     module procedure ATMOS_THERMODYN_cp_0D
     module procedure ATMOS_THERMODYN_cp_3D
  end interface ATMOS_THERMODYN_cp

  interface ATMOS_THERMODYN_cv
     module procedure ATMOS_THERMODYN_cv_0D
     module procedure ATMOS_THERMODYN_cv_3D
  end interface ATMOS_THERMODYN_cv

  interface ATMOS_THERMODYN_rhoe
     module procedure ATMOS_THERMODYN_rhoe_0D
     module procedure ATMOS_THERMODYN_rhoe_3D
  end interface ATMOS_THERMODYN_rhoe

  interface ATMOS_THERMODYN_rhot
     module procedure ATMOS_THERMODYN_rhot_0D
     module procedure ATMOS_THERMODYN_rhot_3D
  end interface ATMOS_THERMODYN_rhot

  interface ATMOS_THERMODYN_pott
     module procedure ATMOS_THERMODYN_pott_0D
     module procedure ATMOS_THERMODYN_pott_3D
  end interface ATMOS_THERMODYN_pott

  interface ATMOS_THERMODYN_temp_pres
     module procedure ATMOS_THERMODYN_temp_pres_0D
     module procedure ATMOS_THERMODYN_temp_pres_3D
  end interface ATMOS_THERMODYN_temp_pres

  interface ATMOS_THERMODYN_temp_pres_E
     module procedure ATMOS_THERMODYN_temp_pres_E_0D
     module procedure ATMOS_THERMODYN_temp_pres_E_3D
  end interface ATMOS_THERMODYN_temp_pres_E

  interface ATMOS_THERMODYN_templhv
     module procedure ATMOS_THERMODYN_templhv_0D
     module procedure ATMOS_THERMODYN_templhv_2D
     module procedure ATMOS_THERMODYN_templhv_3D
  end interface ATMOS_THERMODYN_templhv

  interface ATMOS_THERMODYN_templhs
     module procedure ATMOS_THERMODYN_templhs_0D
     module procedure ATMOS_THERMODYN_templhs_2D
     module procedure ATMOS_THERMODYN_templhs_3D
  end interface ATMOS_THERMODYN_templhs

  interface ATMOS_THERMODYN_templhf
     module procedure ATMOS_THERMODYN_templhf_0D
     module procedure ATMOS_THERMODYN_templhf_2D
     module procedure ATMOS_THERMODYN_templhf_3D
  end interface ATMOS_THERMODYN_templhf

  interface ATMOS_THERMODYN_entr
     module procedure ATMOS_THERMODYN_entr_0D
     module procedure ATMOS_THERMODYN_entr_2D
     module procedure ATMOS_THERMODYN_entr_3D
  end interface ATMOS_THERMODYN_entr

  public :: ATMOS_THERMODYN_tempre
  public :: ATMOS_THERMODYN_tempre2

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public              :: THERMODYN_EMASK = 0.0_RP !< =0: SIMPLE, 1: EXACT

  real(RP), public, allocatable :: AQ_CP(:) !< CP for each hydrometeors [J/kg/K]
  real(RP), public, allocatable :: AQ_CV(:) !< CV for each hydrometeors [J/kg/K]

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
  !> Setup
  subroutine ATMOS_THERMODYN_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CPvap          => CONST_CPvap,          &
       CVvap          => CONST_CVvap,          &
       CL             => CONST_CL,             &
       CI             => CONST_CI,             &
       THERMODYN_TYPE => CONST_THERMODYN_TYPE
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[THERMODYN] / Categ[ATMOS SHARE] / Origin[SCALElib]'

    allocate( AQ_CP(QQA) )
    allocate( AQ_CV(QQA) )

    if ( THERMODYN_TYPE == 'EXACT' ) then
       THERMODYN_EMASK = 1.0_RP

       AQ_CP(I_QV) = CPvap
       AQ_CV(I_QV) = CVvap

       if ( QWS /= 0 ) then
          do n = QWS, QWE
             AQ_CP(n) = CL
             AQ_CV(n) = CL
          enddo
       endif

       if ( QIS /= 0 ) then
          do n = QIS, QIE
             AQ_CP(n) = CI
             AQ_CV(n) = CI
          enddo
       endif
    elseif( THERMODYN_TYPE == 'SIMPLE' ) then
       THERMODYN_EMASK = 0.0_RP

       AQ_CP(I_QV) = CPdry
       AQ_CV(I_QV) = CVdry

       if ( QWS /= 0 ) then
          do n = QWS, QWE
             AQ_CP(n) = CVdry
             AQ_CV(n) = CVdry
          enddo
       endif

       if ( QIS /= 0 ) then
          do n = QIS, QIE
             AQ_CP(n) = CVdry
             AQ_CV(n) = CVdry
          enddo
       endif
    endif

    return
  end subroutine ATMOS_THERMODYN_setup

  !-----------------------------------------------------------------------------
  !> calc dry air mass (0D)
  subroutine ATMOS_THERMODYN_qd_0D( &
       qdry, &
       q     )
    implicit none

    real(RP), intent(out) :: qdry  !< dry mass concentration [kg/kg]
    real(RP), intent(in)  :: q(QA) !< mass concentration     [kg/kg]

    integer :: iqw
    !-----------------------------------------------------------------------------

    qdry = 1.0_RP
#ifndef DRY
    do iqw = QQS, QQE
       qdry = qdry - q(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_qd_0D

  !-----------------------------------------------------------------------------
  !> calc dry air mass (3D)
  subroutine ATMOS_THERMODYN_qd_3D( &
       qdry, &
       q     )
    implicit none

    real(RP), intent(out) :: qdry(KA,IA,JA)    !< dry mass concentration [kg/kg]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration     [kg/kg]

    integer :: k, i, j, iqw
    !-----------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

!       qdry(k,i,j) = 1.0_RP
!       do iqw = QQS, QQE
!          qdry(k,i,j) = qdry(k,i,j) - q(k,i,j,iqw)
!       enddo
       CALC_QDRY( qdry(k,i,j), q, k, i, j, iqw )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_qd_3D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CP,0D)
  subroutine ATMOS_THERMODYN_cp_0D( &
       CPtot, &
       q,     &
       qdry   )
    implicit none

    real(RP), intent(out) :: CPtot !< total specific heat    [J/kg/K]
    real(RP), intent(in)  :: q(QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: qdry  !< dry mass concentration [kg/kg]

    integer :: iqw
    !---------------------------------------------------------------------------

    CPtot = qdry * CPdry
#ifndef DRY
    do iqw = QQS, QQE
       CPtot = CPtot + q(iqw) * AQ_CP(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_cp_0D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CP,3D)
  subroutine ATMOS_THERMODYN_cp_3D( &
       CPtot, &
       q,     &
       qdry   )
    implicit none

    real(RP), intent(out) :: CPtot(KA,IA,JA)    !< total specific heat    [J/kg/K]
    real(RP), intent(in)  :: q    (KA,IA,JA,QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: qdry (KA,IA,JA)    !< dry mass concentration [kg/kg]

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

!       CPtot(k,i,j) = qdry(k,i,j) * CPdry
!       do iqw = QQS, QQE
!          CPtot(k,i,j) = CPtot(k,i,j) + q(k,i,j,iqw) * AQ_CP(iqw)
!       enddo

       CALC_CP(cptot(k,i,j), qdry(k,i,j), q, k, i, j, iqw, CPdry, AQ_CP)

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_cp_3D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CV,0D)
  subroutine ATMOS_THERMODYN_cv_0D( &
       CVtot, &
       q,     &
       qdry   )
    implicit none

    real(RP), intent(out) :: CVtot !< total specific heat    [J/kg/K]
    real(RP), intent(in)  :: q(QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: qdry  !< dry mass concentration [kg/kg]

    integer :: iqw
    !---------------------------------------------------------------------------

    CVtot = qdry * CVdry
#ifndef DRY
    do iqw = QQS, QQE
       CVtot = CVtot + q(iqw) * AQ_CV(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_cv_0D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CV,3D)
  subroutine ATMOS_THERMODYN_cv_3D( &
       CVtot, &
       q,     &
       qdry   )
    implicit none

    real(RP), intent(out) :: CVtot(KA,IA,JA)    !< total specific heat    [J/kg/K]
    real(RP), intent(in)  :: q    (KA,IA,JA,QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: qdry (KA,IA,JA)    !< dry mass concentration [kg/kg]

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

!       CVtot(k,i,j) = qdry(k,i,j) * CVdry
!       do iqw = QQS, QQE
!          CVtot(k,i,j) = CVtot(k,i,j) + q(k,i,j,iqw) * AQ_CV(iqw)
!       enddo

       CALC_CV(cvtot(k,i,j), qdry(k,i,j), q, k, i, j, iqw, CVdry, AQ_CV)

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_cv_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_pott_0D( &
       pott, &
       temp, &
       pres, &
       q     )
    implicit none

    real(RP), intent(out) :: pott  ! potential temperature [K]
    real(RP), intent(in)  :: temp  ! temperature           [K]
    real(RP), intent(in)  :: pres  ! pressure              [Pa]
    real(RP), intent(in)  :: q(QA) ! mass concentration    [kg/kg]

    real(RP) :: qdry
    real(RP) :: Rtot, CVtot, RovCP

    integer  :: iqw
    !---------------------------------------------------------------------------

#ifdef DRY
    CVtot = CVdry
    Rtot  = Rdry
#else
    qdry  = 1.0_RP
    CVtot = 0.0_RP
    do iqw = QQS, QQE
       qdry  = qdry  - q(iqw)
       CVtot = CVtot + q(iqw) * AQ_CV(iqw)
    enddo
    CVtot = CVdry * qdry + CVtot
    Rtot  = Rdry  * qdry + Rvap * q(I_QV)
#endif

    RovCP = Rtot / ( CVtot + Rtot )

    pott  = temp * ( PRE00 / pres )**RovCP

    return
  end subroutine ATMOS_THERMODYN_pott_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_pott_3D( &
       pott, &
       temp, &
       pres, &
       q     )
    implicit none

    real(RP), intent(out) :: pott(KA,IA,JA)    ! potential temperature [K]
    real(RP), intent(in)  :: temp(KA,IA,JA)    ! temperature           [K]
    real(RP), intent(in)  :: pres(KA,IA,JA)    ! pressure              [Pa]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) ! mass concentration    [kg/kg]

    real(RP) :: qdry
    real(RP) :: Rtot, CVtot, RovCP

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
#ifdef DRY
       CVtot = CVdry
       Rtot  = Rdry
#else
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - q(k,i,j,iqw)
          CVtot = CVtot + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * q(k,i,j,I_QV)
#endif

       RovCP = Rtot / ( CVtot + Rtot )

       pott(k,i,j) = temp(k,i,j) * ( PRE00 / pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_pott_3D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> rho * ein (0D)
  subroutine ATMOS_THERMODYN_rhoe_0D( &
       rhoe, &
       rhot, &
       q     )
    implicit none

    real(RP), intent(out) :: rhoe  !< density * internal energy       [J/m3]
    real(RP), intent(in)  :: rhot  !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: q(QA) !< mass concentration              [kg/kg]

    real(RP) :: qdry
    real(RP) :: pres
    real(RP) :: Rtot, CVtot, CPovCV

    integer :: iqw
    !---------------------------------------------------------------------------

#ifdef DRY
    CVtot = CVdry
    Rtot  = Rdry
#else
    qdry  = 1.0_RP
    CVtot = 0.0_RP
    do iqw = QQS, QQE
       qdry  = qdry  - q(iqw)
       CVtot = CVtot + q(iqw) * AQ_CV(iqw)
    enddo
    CVtot = CVdry * qdry + CVtot
    Rtot  = Rdry  * qdry + Rvap * q(I_QV)
#endif

    CPovCV = ( CVtot + Rtot ) / CVtot

    pres = PRE00 * ( rhot * Rtot / PRE00 )**CPovCV

    rhoe = pres / Rtot * CVtot

    return
  end subroutine ATMOS_THERMODYN_rhoe_0D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> rho * ein (3D)
  subroutine ATMOS_THERMODYN_rhoe_3D( &
       rhoe, &
       rhot, &
       q     )
    implicit none

    real(RP), intent(out) :: rhoe(KA,IA,JA)    !< density * internal energy       [J/m3]
    real(RP), intent(in)  :: rhot(KA,IA,JA)    !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration              [kg/kg]

    real(RP) :: qdry
    real(RP) :: pres
    real(RP) :: Rtot, CVtot, CPovCV

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw,qdry,pres,Rtot,CVtot,CPovCV) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
#ifdef DRY
       CVtot = CVdry
       Rtot  = Rdry
#else
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - q(k,i,j,iqw)
          CVtot = CVtot + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * q(k,i,j,I_QV)
#endif


       CPovCV = ( CVtot + Rtot ) / CVtot

       pres = PRE00 * ( rhot(k,i,j) * Rtot / PRE00 )**CPovCV

       rhoe(k,i,j) = pres / Rtot * CVtot
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_rhoe_3D

  !-----------------------------------------------------------------------------
  !> calc rho * ein -> rho * pott (0D)
  subroutine ATMOS_THERMODYN_rhot_0D( &
       rhot, &
       rhoe, &
       q     )
    implicit none

    real(RP), intent(out) :: rhot  !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: rhoe  !< density * internal energy       [J/m3]
    real(RP), intent(in)  :: q(QA) !< mass concentration              [kg/kg]

    real(RP) :: qdry
    real(RP) :: pres
    real(RP) :: Rtot, CVtot, RovCP

    integer :: iqw
    !---------------------------------------------------------------------------

#ifdef DRY
    CVtot = CVdry
    Rtot  = Rdry
#else
    qdry  = 1.0_RP
    CVtot = 0.0_RP
    do iqw = QQS, QQE
       qdry  = qdry  - q(iqw)
       CVtot = CVtot + q(iqw) * AQ_CV(iqw)
    enddo
    CVtot = CVdry * qdry + CVtot
    Rtot  = Rdry  * qdry + Rvap * q(I_QV)
#endif

    RovCP = Rtot / ( CVtot + Rtot )

    pres = rhoe * Rtot / CVtot

    rhot = rhoe / CVtot * ( PRE00 / pres )**RovCP

    return
  end subroutine ATMOS_THERMODYN_rhot_0D

  !-----------------------------------------------------------------------------
  !> calc rho * ein -> rho * pott (3D)
  subroutine ATMOS_THERMODYN_rhot_3D( &
       rhot, &
       rhoe, &
       q     )
    implicit none

    real(RP), intent(out) :: rhot(KA,IA,JA)    !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: rhoe(KA,IA,JA)    !< density * internal energy       [J/m3]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration              [kg/kg]

    real(RP) :: qdry
    real(RP) :: pres
    real(RP) :: Rtot, CVtot, RovCP

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw,qdry,pres,Rtot,CVtot,RovCP) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

#ifdef DRY
       CVtot = CVdry
       Rtot  = Rdry
#else
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - q(k,i,j,iqw)
          CVtot = CVtot + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * q(k,i,j,I_QV)
#endif

       RovCP = Rtot / ( CVtot + Rtot )

       pres = rhoe(k,i,j) * Rtot / CVtot

       rhot(k,i,j) = rhoe(k,i,j) / CVtot * ( PRE00 / pres )**RovCP
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_rhot_3D

  !-----------------------------------------------------------------------------
  !> calc rho, q, rho * pott -> temp & pres (0D)
  subroutine ATMOS_THERMODYN_temp_pres_0D( &
       temp, &
       pres, &
       dens, &
       rhot, &
       q     )
    implicit none

    real(RP), intent(out) :: temp  !< temperature                     [K]
    real(RP), intent(out) :: pres  !< pressure                        [Pa]
    real(RP), intent(in)  :: dens  !< density                         [kg/m3]
    real(RP), intent(in)  :: rhot  !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: q(QA) !< mass concentration              [kg/kg]

    real(RP) :: qdry
    real(RP) :: Rtot, CVtot, CPovCV

    integer :: iqw
    !---------------------------------------------------------------------------

#ifdef DRY
    CVtot = CVdry
    Rtot  = Rdry
#else
    qdry  = 1.0_RP
    CVtot = 0.0_RP
    do iqw = QQS, QQE
       qdry  = qdry  - q(iqw)
       CVtot = CVtot + q(iqw) * AQ_CV(iqw)
    enddo
    CVtot = CVdry * qdry + CVtot
    Rtot  = Rdry  * qdry + Rvap * q(I_QV)
#endif

    CPovCV = ( CVtot + Rtot ) / CVtot

    pres = PRE00 * ( rhot * Rtot / PRE00 )**CPovCV
    temp = pres / ( dens * Rtot )

    return
  end subroutine ATMOS_THERMODYN_temp_pres_0D

  !-----------------------------------------------------------------------------
  !> calc rho, q, rho * pott -> temp & pres (3D)
  subroutine ATMOS_THERMODYN_temp_pres_3D( &
       temp, &
       pres, &
       dens, &
       rhot, &
       q     )
    implicit none

    real(RP), intent(out) :: temp(KA,IA,JA)    !< temperature                     [K]
    real(RP), intent(out) :: pres(KA,IA,JA)    !< pressure                        [Pa]
    real(RP), intent(in)  :: dens(KA,IA,JA)    !< density                         [kg/m3]
    real(RP), intent(in)  :: rhot(KA,IA,JA)    !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration              [kg/kg]

    real(RP) :: qdry
    real(RP) :: Rtot, CVtot, CPovCV

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw,qdry,Rtot,CVtot,CPovCV) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
#ifdef DRY
       CVtot = CVdry
       Rtot  = Rdry
#else
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - q(k,i,j,iqw)
          CVtot = CVtot + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * q(k,i,j,I_QV)
#endif

       CPovCV = ( CVtot + Rtot ) / CVtot

       pres(k,i,j) = PRE00 * ( rhot(k,i,j) * Rtot / PRE00 )**CPovCV
       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rtot )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_temp_pres_3D

  !-----------------------------------------------------------------------------
  !> calc rho, q, rho * pott -> temp & pres (0D)
  subroutine ATMOS_THERMODYN_temp_pres_E_0D( &
       temp, &
       pres, &
       dens, &
       rhoe, &
       q     )
    implicit none

    real(RP), intent(out) :: temp  !< temperature               [K]
    real(RP), intent(out) :: pres  !< pressure                  [Pa]
    real(RP), intent(in)  :: dens  !< density                   [kg/m3]
    real(RP), intent(in)  :: rhoe  !< density * internal energy [J/m3]
    real(RP), intent(in)  :: q(QA) !< mass concentration        [kg/kg]

    real(RP) :: qdry
    real(RP) :: Rtot, CVtot

    integer :: iqw
    !---------------------------------------------------------------------------

#ifdef DRY
    CVtot = CVdry
    Rtot  = Rdry
#else
    qdry  = 1.0_RP
    CVtot = 0.0_RP
    do iqw = QQS, QQE
       qdry  = qdry  - q(iqw)
       CVtot = CVtot + q(iqw) * AQ_CV(iqw)
    enddo
    CVtot = CVdry * qdry + CVtot
    Rtot  = Rdry  * qdry + Rvap * q(I_QV)
#endif

    temp = rhoe / ( dens * CVtot )
    pres = dens * Rtot * temp

    return
  end subroutine ATMOS_THERMODYN_temp_pres_E_0D

  !-----------------------------------------------------------------------------
  !> calc rho, q, rho * pott -> temp & pres (3D)
  subroutine ATMOS_THERMODYN_temp_pres_E_3D( &
       temp, &
       pres, &
       dens, &
       rhoe, &
       q     )
    implicit none

    real(RP), intent(out) :: temp(KA,IA,JA)    !< temperature               [K]
    real(RP), intent(out) :: pres(KA,IA,JA)    !< pressure                  [Pa]
    real(RP), intent(in)  :: dens(KA,IA,JA)    !< density                   [kg/m3]
    real(RP), intent(in)  :: rhoe(KA,IA,JA)    !< density * internal energy [J/m3]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration        [kg/kg]

    real(RP) :: qdry
    real(RP) :: Rtot, CVtot

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw,qdry,Rtot,CVtot) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
#ifdef DRY
       CVtot = CVdry
       Rtot  = Rdry
#else
       qdry  = 1.0_RP
       CVtot = 0.0_RP
       do iqw = QQS, QQE
          qdry  = qdry  - q(k,i,j,iqw)
          CVtot = CVtot + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo
       CVtot = CVdry * qdry + CVtot
       Rtot  = Rdry  * qdry + Rvap * q(k,i,j,I_QV)
#endif

       temp(k,i,j) = rhoe(k,i,j) / ( dens(k,i,j) * CVtot )
       pres(k,i,j) = dens(k,i,j) * Rtot * temp(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_temp_pres_E_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhv_0D( &
       lhv, &
       temp )
    implicit none

    real(RP), intent(out) :: lhv  ! potential temperature [K]
    real(RP), intent(in)  :: temp ! temperature           [K]
    !---------------------------------------------------------------------------

    lhv = LHV0 + ( CPvap - CL ) * ( temp - TEM00 ) * THERMODYN_EMASK

    return
  end subroutine ATMOS_THERMODYN_templhv_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhv_2D( &
       lhv, &
       temp )
    implicit none

    real(RP), intent(out) :: lhv (IA,JA) ! potential temperature [K]
    real(RP), intent(in)  :: temp(IA,JA) ! temperature           [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       lhv(i,j) = LHV0 + ( CPvap - CL ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhv_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhv_3D( &
       lhv, &
       temp )
    implicit none

    real(RP), intent(out) :: lhv (KA,IA,JA) ! potential temperature [K]
    real(RP), intent(in)  :: temp(KA,IA,JA) ! temperature           [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhv(k,i,j) = LHV0 + ( CPvap - CL ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhv_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhs_0D( &
       lhs, &
       temp )
    implicit none

    real(RP), intent(out) :: lhs   ! potential temperature [K]
    real(RP), intent(in)  :: temp  ! temperature           [K]
    !---------------------------------------------------------------------------

    lhs = LHS0 + ( CPvap - CI ) * ( temp - TEM00 ) * THERMODYN_EMASK

    return
  end subroutine ATMOS_THERMODYN_templhs_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhs_2D( &
       lhs, &
       temp )
    implicit none

    real(RP), intent(out) :: lhs (IA,JA) ! potential temperature [K]
    real(RP), intent(in)  :: temp(IA,JA) ! temperature           [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       lhs(i,j) = LHS0 + ( CPvap - CI ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhs_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhs_3D( &
       lhs, &
       temp )
    implicit none

    real(RP), intent(out) :: lhs (KA,IA,JA) ! potential temperature [K]
    real(RP), intent(in)  :: temp(KA,IA,JA) ! temperature           [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhs(k,i,j) = LHS0 + ( CPvap - CI ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhs_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhf_0D( &
       lhf, &
       temp )
    implicit none

    real(RP), intent(out) :: lhf  ! potential temperature [K]
    real(RP), intent(in)  :: temp ! temperature           [K]
    !---------------------------------------------------------------------------

    lhf = LHF0 + ( CL - CI ) * ( temp - TEM00 ) * THERMODYN_EMASK

    return
  end subroutine ATMOS_THERMODYN_templhf_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhf_2D( &
       lhf, &
       temp )
    implicit none

    real(RP), intent(out) :: lhf (IA,JA) ! potential temperature [K]
    real(RP), intent(in)  :: temp(IA,JA) ! temperature           [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
       lhf(i,j) = LHF0 + ( CL - CI ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhf_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhf_3D( &
       lhf, &
       temp )
    implicit none

    real(RP), intent(out) :: lhf (KA,IA,JA) ! potential temperature [K]
    real(RP), intent(in)  :: temp(KA,IA,JA) ! temperature           [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhf(k,i,j) = LHF0 + ( CL - CI ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhf_3D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (0D)
  subroutine ATMOS_THERMODYN_entr_0D( &
       entr, &
       temp, &
       pres, &
       q     )
    implicit none

    real(RP), intent(out) :: entr  !< entropy            [J/K]
    real(RP), intent(in)  :: temp  !< temperature        [K]
    real(RP), intent(in)  :: pres  !< pressure           [Pa]
    real(RP), intent(in)  :: q(QA) !< mass concentration [kg/kg]

    real(RP) :: lhv, lhf
    real(RP) :: qdry, Rtot
    real(RP) :: logT_T0, pres_dry, pres_vap
    real(RP) :: sq

    integer :: iqw
    !---------------------------------------------------------------------------

    lhv = LHV0 + ( CPvap - CL ) * ( temp - TEM00 ) * THERMODYN_EMASK
    lhf = LHF0 + ( CL    - CI ) * ( temp - TEM00 ) * THERMODYN_EMASK

    qdry = 1.0_RP
    do iqw = QQS, QQE
       qdry = qdry - q(iqw)
    enddo
    Rtot = Rdry * qdry + Rvap * q(I_QV)

    logT_T0 = log( temp / TEM00 )

    ! dry air + vapor
    pres_dry = max( pres * qdry    * Rdry / Rtot, EPS )
    pres_vap = max( pres * q(I_QV) * Rvap / Rtot, EPS )

    entr = qdry    * CPdry * logT_T0                 &
         - qdry    * Rdry  * log( pres_dry / PRE00 ) &
         + q(I_QV) * CPvap * logT_T0                 &
         - q(I_QV) * Rvap  * log( pres_vap / PSAT0 ) &
         + q(I_QV) * lhv / TEM00

    ! hydrometeors
    do iqw = QQS, QQE
       sq = 0.0_RP

       if( iqw == I_QC ) sq = q(I_QC) * CL * logT_T0

       if( iqw == I_QR ) sq = q(I_QR) * CL * logT_T0

       if( iqw == I_QI ) sq = q(I_QI) * CI * logT_T0 &
                            - q(I_QI) * lhf / TEM00

       if( iqw == I_QS ) sq = q(I_QS) * CI * logT_T0 &
                            - q(I_QS) * lhf / TEM00

       if( iqw == I_QG ) sq = q(I_QG) * CI * logT_T0 &
                            - q(I_QG) * lhf / TEM00

       entr = entr + sq
    enddo

    return
  end subroutine ATMOS_THERMODYN_entr_0D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (2D)
  subroutine ATMOS_THERMODYN_entr_2D( &
       entr, &
       temp, &
       pres, &
       q     )
    implicit none

    real(RP), intent(out) :: entr(IA,JA)    !< entropy            [J/K]
    real(RP), intent(in)  :: temp(IA,JA)    !< temperature        [K]
    real(RP), intent(in)  :: pres(IA,JA)    !< pressure           [Pa]
    real(RP), intent(in)  :: q   (IA,JA,QA) !< mass concentration [kg/kg]

    real(RP) :: lhv, lhf
    real(RP) :: qdry, Rtot
    real(RP) :: logT_T0, pres_dry, pres_vap
    real(RP) :: sq

    integer :: i, j, iqw
    !---------------------------------------------------------------------------

    ! dry air + vapor
    do j = JSB, JEB
    do i = ISB, IEB
       lhv = LHV0 + ( CPvap - CL ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK

       qdry = 1.0_RP
       do iqw = QQS, QQE
          qdry = qdry - q(i,j,iqw)
       enddo
       Rtot = Rdry * qdry + Rvap * q(i,j,I_QV)

       logT_T0 = log( temp(i,j) / TEM00 )

       pres_dry = max( pres(i,j) * qdry        * Rdry / Rtot, EPS )
       pres_vap = max( pres(i,j) * q(i,j,I_QV) * Rvap / Rtot, EPS )

       entr(i,j) = qdry        * CPdry * logT_T0 &
                 - qdry        * Rdry  * log( pres_dry / PRE00 ) &
                 + q(i,j,I_QV) * CPvap * logT_T0 &
                 - q(i,j,I_QV) * Rvap  * log( pres_vap / PSAT0 ) &
                 + q(i,j,I_QV) * lhv / TEM00
    enddo
    enddo

    ! hydrometeors
    do j = JSB, JEB
    do i = ISB, IEB
       lhf = LHF0 + ( CL - CI ) * ( temp(i,j) - TEM00 ) * THERMODYN_EMASK

       logT_T0 = log( temp(i,j) / TEM00 )

       do iqw = QQS, QQE
          sq = 0.0_RP

          if( iqw == I_QC ) sq = q(i,j,I_QC) *   CL * logT_T0
          if( iqw == I_QR ) sq = q(i,j,I_QR) *   CL * logT_T0
          if( iqw == I_QI ) sq = q(i,j,I_QI) * ( CI * logT_T0 - lhf / TEM00 )
          if( iqw == I_QS ) sq = q(i,j,I_QS) * ( CI * logT_T0 - lhf / TEM00 )
          if( iqw == I_QG ) sq = q(i,j,I_QG) * ( CI * logT_T0 - lhf / TEM00 )

          entr(i,j) = entr(i,j) + sq
       enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_entr_2D

  !-----------------------------------------------------------------------------
  !> calc temp, pres, q -> entropy (3D)
  subroutine ATMOS_THERMODYN_entr_3D( &
       entr, &
       temp, &
       pres, &
       q     )
    implicit none

    real(RP), intent(out) :: entr(KA,IA,JA)    !< entropy            [J/K]
    real(RP), intent(in)  :: temp(KA,IA,JA)    !< temperature        [K]
    real(RP), intent(in)  :: pres(KA,IA,JA)    !< pressure           [Pa]
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration [kg/kg]

    real(RP) :: lhv, lhf
    real(RP) :: qdry, Rtot
    real(RP) :: logT_T0, pres_dry, pres_vap
    real(RP) :: sq

    integer :: k, i, j, iqw
    !---------------------------------------------------------------------------

    ! dry air + vapor
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhv = LHV0 + ( CPvap - CL ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK

       qdry = 1.0_RP
       do iqw = QQS, QQE
          qdry = qdry - q(k,i,j,iqw)
       enddo
       Rtot = Rdry * qdry + Rvap * q(k,i,j,I_QV)

       logT_T0 = log( temp(k,i,j) / TEM00 )

       pres_dry = max( pres(k,i,j) * qdry          * Rdry / Rtot, EPS )
       pres_vap = max( pres(k,i,j) * q(k,i,j,I_QV) * Rvap / Rtot, EPS )

       entr(k,i,j) = qdry          * CPdry * logT_T0 &
                   - qdry          * Rdry  * log( pres_dry / PRE00 ) &
                   + q(k,i,j,I_QV) * CPvap * logT_T0 &
                   - q(k,i,j,I_QV) * Rvap  * log( pres_vap / PSAT0 ) &
                   + q(k,i,j,I_QV) * lhv / TEM00
    enddo
    enddo
    enddo

    ! hydrometeors
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       lhf = LHF0 + ( CL - CI ) * ( temp(k,i,j) - TEM00 ) * THERMODYN_EMASK

       logT_T0 = log( temp(k,i,j) / TEM00 )

       do iqw = QQS, QQE
          sq = 0.0_RP

          if( iqw == I_QC ) sq = q(k,i,j,I_QC) *   CL * logT_T0
          if( iqw == I_QR ) sq = q(k,i,j,I_QR) *   CL * logT_T0
          if( iqw == I_QI ) sq = q(k,i,j,I_QI) * ( CI * logT_T0 - lhf / TEM00 )
          if( iqw == I_QS ) sq = q(k,i,j,I_QS) * ( CI * logT_T0 - lhf / TEM00 )
          if( iqw == I_QG ) sq = q(k,i,j,I_QG) * ( CI * logT_T0 - lhf / TEM00 )

          entr(k,i,j) = entr(k,i,j) + sq
       enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_entr_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_tempre( &
      temp, pres,         &
      Ein,  dens, qdry, q )
    implicit none

    real(RP), intent(out) :: temp(KA,IA,JA)    ! temperature
    real(RP), intent(out) :: pres(KA,IA,JA)    ! pressure
    real(RP), intent(in)  :: Ein (KA,IA,JA)    ! internal energy
    real(RP), intent(in)  :: dens(KA,IA,JA)    ! density
    real(RP), intent(in)  :: qdry(KA,IA,JA)    ! dry concentration
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) ! water concentration

    real(RP) :: cv, Rmoist

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw,cv,Rmoist) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

       CALC_CV(cv, qdry(k,i,j), q, k, i, j, iqw, CVdry, AQ_CV)
       CALC_R(Rmoist, q(k,i,j,I_QV), qdry(k,i,j), Rdry, Rvap)

       temp(k,i,j) = Ein(k,i,j) / cv

       pres(k,i,j) = dens(k,i,j) * Rmoist * temp(k,i,j)

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_tempre

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_tempre2( &
      temp, pres,         &
      dens, pott, qdry, q )
    implicit none

    real(RP), intent(out) :: temp(KA,IA,JA)    ! temperature
    real(RP), intent(out) :: pres(KA,IA,JA)    ! pressure
    real(RP), intent(in)  :: dens(KA,IA,JA)    ! density
    real(RP), intent(in)  :: pott(KA,IA,JA)    ! potential temperature
    real(RP), intent(in)  :: qdry(KA,IA,JA)    ! dry concentration
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) ! water concentration

    real(RP) :: Rmoist, cp

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,iqw,cp,Rmoist) OMP_SCHEDULE_ collapse(2)
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE

       CALC_CP(cp, qdry(k,i,j), q, k, i, j, iqw, CPdry, AQ_CP)
       CALC_R(Rmoist, q(k,i,j,I_QV), qdry(k,i,j), Rdry, Rvap)
       CALC_PRE(pres(k,i,j), dens(k,i,j), pott(k,i,j), Rmoist, cp, PRE00)

       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rmoist )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_tempre2

end module scale_atmos_thermodyn
