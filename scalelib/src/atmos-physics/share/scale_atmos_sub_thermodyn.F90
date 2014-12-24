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
     Rdry  => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     Rvap  => CONST_Rvap,  &
     CPvap => CONST_CPvap, &
     CVvap => CONST_CVvap, &
     CL    => CONST_CL,    &
     CI    => CONST_CI,    &
     LH0   => CONST_LH0,   &
     TEM00 => CONST_TEM00, &
     PRE00 => CONST_PRE00
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
  public :: ATMOS_THERMODYN_rhoe
  public :: ATMOS_THERMODYN_rhot
  public :: ATMOS_THERMODYN_temp_pres
  public :: ATMOS_THERMODYN_temp_pres_E
  public :: ATMOS_THERMODYN_templhv

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

  public :: ATMOS_THERMODYN_tempre
  public :: ATMOS_THERMODYN_tempre2

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[THERMODYN]/Categ[ATMOS]'

    allocate( AQ_CP(QQA) )
    allocate( AQ_CV(QQA) )

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

    RovCP  = Rtot / ( CVtot + Rtot )

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

#ifdef DRY
       CVtot = CVdry
       Rtot  = Rtot
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

       RovCP  = Rtot / ( CVtot + Rtot )

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

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
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA

       CALC_CP(cp, qdry(k,i,j), q, k, i, j, iqw, CPdry, AQ_CP)
       CALC_R(Rmoist, q(k,i,j,I_QV), qdry(k,i,j), Rdry, Rvap)
       CALC_PRE(pres(k,i,j), dens(k,i,j), pott(k,i,j), Rmoist, cp, PRE00)

       temp(k,i,j) = pres(k,i,j) / ( dens(k,i,j) * Rmoist )

    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_tempre2

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhv_0D( &
      lhv,         &
      temp         )
!    use scale_const, only: &
!       THERMODYN_TYPE => CONST_THERMODYN_TYPE
    implicit none

    real(RP), intent(out) :: lhv  ! potential temperature [K]
    real(RP), intent(in)  :: temp  ! temperature           [K]

    !---------------------------------------------------------------------------

!    if ( THERMODYN_TYPE == 'EXACT' ) then
         lhv = LH0 + ( CPvap - CL ) * ( temp - TEM00 )
!    elseif ( THERMODYN_TYPE == 'SIMPLE' ) then
!         lhv = LH0
!    endif

    return
  end subroutine ATMOS_THERMODYN_templhv_0D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhv_2D( &
      lhv,         &
      temp         )
!    use scale_const, only: &
!       THERMODYN_TYPE => CONST_THERMODYN_TYPE
    implicit none

    real(RP), intent(out) :: lhv(IA,JA)     ! potential temperature [K]
    real(RP), intent(in)  :: temp(IA,JA)    ! temperature           [K]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
!       if ( THERMODYN_TYPE == 'EXACT' ) then
         lhv(i,j) = LH0 + ( CPvap - CL ) * ( temp(i,j) - TEM00 )
!       elseif ( THERMODYN_TYPE == 'SIMPLE' ) then
!         lhv(i,j) = LH0
!       endif
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhv_2D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_templhv_3D( &
      lhv,         &
      temp         )
!    use scale_const, only: &
!       THERMODYN_TYPE => CONST_THERMODYN_TYPE
    implicit none

    real(RP), intent(out) :: lhv(KA,IA,JA)     ! potential temperature [K]
    real(RP), intent(in)  :: temp(KA,IA,JA)    ! temperature           [K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
!       if ( THERMODYN_TYPE == 'EXACT' ) then
         lhv(k,i,j) = LH0 + ( CPvap - CL ) * ( temp(k,i,j) - TEM00 )
!       elseif ( THERMODYN_TYPE == 'SIMPLE' ) then
!         lhv(k,i,j) = LH0
!       endif
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_templhv_3D

end module scale_atmos_thermodyn
