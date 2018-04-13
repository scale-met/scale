!-------------------------------------------------------------------------------
!> module atmosphere / thermodyn
!!
!! @par Description
!!          Thermodynamics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_thermodyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  use scale_const, only: &
     Rdry  => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     PRE00 => CONST_PRE00
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_THERMODYN_setup

  public :: ATMOS_THERMODYN_qdry
  public :: ATMOS_THERMODYN_cv
  public :: ATMOS_THERMODYN_cp
  public :: ATMOS_THERMODYN_r
  public :: ATMOS_THERMODYN_specific_heat
  public :: ATMOS_THERMODYN_rhot2pres
  public :: ATMOS_THERMODYN_rhot2rhoe
  public :: ATMOS_THERMODYN_rhoe2rhot
  public :: ATMOS_THERMODYN_rhot2temp_pres
  public :: ATMOS_THERMODYN_rhoe2temp_pres
  public :: ATMOS_THERMODYN_ein2temp_pres
  public :: ATMOS_THERMODYN_pott2temp_pres
  public :: ATMOS_THERMODYN_temp_pres2pott

  interface ATMOS_THERMODYN_qdry
     module procedure ATMOS_THERMODYN_qdry_0D
     module procedure ATMOS_THERMODYN_qdry_3D
  end interface ATMOS_THERMODYN_qdry

  interface ATMOS_THERMODYN_cv
     module procedure ATMOS_THERMODYN_cv_0D
     module procedure ATMOS_THERMODYN_cv_3D
  end interface ATMOS_THERMODYN_cv

  interface ATMOS_THERMODYN_cp
     module procedure ATMOS_THERMODYN_cp_0D
     module procedure ATMOS_THERMODYN_cp_3D
  end interface ATMOS_THERMODYN_cp

  interface ATMOS_THERMODYN_r
     module procedure ATMOS_THERMODYN_r_0D
     module procedure ATMOS_THERMODYN_r_3D
  end interface ATMOS_THERMODYN_r

  interface ATMOS_THERMODYN_specific_heat
     module procedure ATMOS_THERMODYN_specific_heat_0D
     module procedure ATMOS_THERMODYN_specific_heat_1D
     module procedure ATMOS_THERMODYN_specific_heat_3D
  end interface ATMOS_THERMODYN_specific_heat

  interface ATMOS_THERMODYN_rhot2pres
     module procedure ATMOS_THERMODYN_rhot2pres_0D
     module procedure ATMOS_THERMODYN_rhot2pres_3D
  end interface ATMOS_THERMODYN_rhot2pres

  interface ATMOS_THERMODYN_rhot2rhoe
     module procedure ATMOS_THERMODYN_rhot2rhoe_0D
     module procedure ATMOS_THERMODYN_rhot2rhoe_3D
  end interface ATMOS_THERMODYN_rhot2rhoe

  interface ATMOS_THERMODYN_rhoe2rhot
     module procedure ATMOS_THERMODYN_rhoe2rhot_0D
     module procedure ATMOS_THERMODYN_rhoe2rhot_3D
  end interface ATMOS_THERMODYN_rhoe2rhot

  interface ATMOS_THERMODYN_rhoe2temp_pres
     module procedure ATMOS_THERMODYN_rhoe2temp_pres_0D
     module procedure ATMOS_THERMODYN_rhoe2temp_pres_3D
  end interface ATMOS_THERMODYN_rhoe2temp_pres

  interface ATMOS_THERMODYN_rhot2temp_pres
     module procedure ATMOS_THERMODYN_rhot2temp_pres_0D
     module procedure ATMOS_THERMODYN_rhot2temp_pres_1D
     module procedure ATMOS_THERMODYN_rhot2temp_pres_3D
  end interface ATMOS_THERMODYN_rhot2temp_pres

  interface ATMOS_THERMODYN_ein2temp_pres
     module procedure ATMOS_THERMODYN_ein2temp_pres_0D
     module procedure ATMOS_THERMODYN_ein2temp_pres_3D
  end interface ATMOS_THERMODYN_ein2temp_pres

  interface ATMOS_THERMODYN_pott2temp_pres
     module procedure ATMOS_THERMODYN_pott2temp_pres_0D
     module procedure ATMOS_THERMODYN_pott2temp_pres_3D
  end interface ATMOS_THERMODYN_pott2temp_pres

  interface ATMOS_THERMODYN_temp_pres2pott
     module procedure ATMOS_THERMODYN_temp_pres2pott_0D
     module procedure ATMOS_THERMODYN_temp_pres2pott_1D
     module procedure ATMOS_THERMODYN_temp_pres2pott_3D
  end interface ATMOS_THERMODYN_temp_pres2pott


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
  !> Setup
  subroutine ATMOS_THERMODYN_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_THERMODYN_setup",*) 'Setup'
    LOG_INFO("ATMOS_THERMODYN_setup",*) 'No namelists.'

    return
  end subroutine ATMOS_THERMODYN_setup

  !-----------------------------------------------------------------------------
  !> calc dry air mass (0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_THERMODYN_qdry_0D( &
       QA, &
       q, q_mass, &
       qdry       )
    implicit none
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q(QA)      !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: q_mass(QA) !< mass factor 0 or 1

    real(RP), intent(out) :: qdry       !< dry mass concentration [kg/kg]

    integer :: iqw
    !---------------------------------------------------------------------------

    qdry = 1.0_RP
#ifndef DRY
    do iqw = 1, QA
       qdry = qdry - q(iqw)*q_mass(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_qdry_0D

  !-----------------------------------------------------------------------------
  !> calc dry air mass (3D)
  subroutine ATMOS_THERMODYN_qdry_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, QA, &
       q, q_mass, &
       qdry       )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q     (KA,IA,JA,QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: q_mass(QA)          !< mass factor 0 or 1
    real(RP), intent(out) :: qdry  (KA,IA,JA)    !< dry mass concentration [kg/kg]

    integer :: k, i, j
    !-----------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,qdry,q,q_mass,QA)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_qdry( QA, q(k,i,j,:), q_mass(:), qdry(k,i,j) )
    enddo
    enddo
    enddo
    return
  end subroutine ATMOS_THERMODYN_qdry_3D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CV,0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_THERMODYN_cv_0D( &
       QA, &
       q, CVq, qdry, &
       CVtot         )
    implicit none
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q(QA)   !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: CVq(QA) !< specific heat          [J/kg/K]
    real(RP), intent(in)  :: qdry    !< dry mass concentration [kg/kg]

    real(RP), intent(out) :: CVtot   !< total specific heat    [J/kg/K]

    integer :: iqw
    !---------------------------------------------------------------------------

    CVtot = qdry * CVdry
#ifndef DRY
    do iqw = 1, QA
       CVtot = CVtot + q(iqw) * CVq(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_cv_0D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CV,3D)
  subroutine ATMOS_THERMODYN_cv_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, QA, &
       q, CVq, qdry, &
       CVtot         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q    (KA,IA,JA,QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: CVq  (QA)          !< specific heat          [J/kg/K]
    real(RP), intent(in)  :: qdry (KA,IA,JA)    !< dry mass concentration [kg/kg]

    real(RP), intent(out) :: CVtot(KA,IA,JA)    !< total specific heat    [J/kg/K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,cvtot,qdry,q,CVdry,CVq,QA)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_cv( QA, q(k,i,j,:), CVq(:), qdry(k,i,j), CVtot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_cv_3D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CP,0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_THERMODYN_cp_0D( &
       QA, &
       q, CPq, qdry, &
       CPtot         )
    implicit none
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q(QA)   !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: CPq(QA) !< specific heat          [J/kg/K]
    real(RP), intent(in)  :: qdry    !< dry mass concentration [kg/kg]

    real(RP), intent(out) :: CPtot   !< total specific heat    [J/kg/K]

    integer :: iqw
    !---------------------------------------------------------------------------

    CPtot = qdry * CPdry
#ifndef DRY
    do iqw = 1, QA
       CPtot = CPtot + q(iqw) * CPq(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_cp_0D

  !-----------------------------------------------------------------------------
  !> calc total specific heat (CP,3D)
  subroutine ATMOS_THERMODYN_cp_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, QA, &
       q, CPq, qdry, &
       CPtot         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q    (KA,IA,JA,QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: CPq  (QA)          !< specific heat          [J/kg/K]
    real(RP), intent(in)  :: qdry (KA,IA,JA)    !< dry mass concentration [kg/kg]

    real(RP), intent(out) :: CPtot(KA,IA,JA)    !< total specific heat    [J/kg/K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_cp( QA, q(k,i,j,:), CPq(:), qdry(k,i,j), CPtot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_cp_3D

  !-----------------------------------------------------------------------------
  !> calc total gas constant (R,0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_THERMODYN_r_0D( &
       QA, &
       q, Rq, qdry, &
       Rtot         )
    implicit none
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q(QA)   !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: Rq(QA)  !< gas constant           [J/kg/K]
    real(RP), intent(in)  :: qdry    !< dry mass concentration [kg/kg]

    real(RP), intent(out) :: Rtot    !< total gas constant     [J/kg/K]

    integer :: iqw
    !---------------------------------------------------------------------------

    Rtot = qdry * Rdry
#ifndef DRY
    do iqw = 1, QA
       Rtot = Rtot + q(iqw) * Rq(iqw)
    enddo
#endif

    return
  end subroutine ATMOS_THERMODYN_r_0D

  !-----------------------------------------------------------------------------
  !> calc total gas constant (R,3D)
  subroutine ATMOS_THERMODYN_r_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, QA, &
       q, Rq, qdry, &
       Rtot         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q   (KA,IA,JA,QA) !< mass concentration     [kg/kg]
    real(RP), intent(in)  :: Rq  (QA)          !< gas constant           [J/kg/K]
    real(RP), intent(in)  :: qdry(KA,IA,JA)    !< dry mass concentration [kg/kg]

    real(RP), intent(out) :: Rtot(KA,IA,JA)    !< total gas constant     [J/kg/K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_r( QA, q(k,i,j,:), Rq(:), qdry(k,i,j), Rtot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_r_3D

  !-----------------------------------------------------------------------------
  !> ATMOS_THERMODYN_specific_heat_0D
  !! calc heat capacities: qdry, Rtot, CVtot, CPtot
  !<
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_THERMODYN_specific_heat_0D( &
       QA, &
       q, Mq, Rq, CVq, CPq,     &
       Qdry, Rtot, CVtot, CPtot )
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q  (QA) !< mass concentration      [kg/kg]
    real(RP), intent(in)  :: Mq (QA) !< mass of tracer          0 or 1
    real(RP), intent(in)  :: Rq (QA) !< gas constant of tracer  [J/kg/K]
    real(RP), intent(in)  :: CVq(QA) !< specific heat of tracer [J/kg/K]
    real(RP), intent(in)  :: CPq(QA) !< specific heat of tracer [J/kg/K]

    real(RP), intent(out) :: Qdry  !> dry mass concentration [kg/kg]
    real(RP), intent(out) :: Rtot  !> total gas constant     [J/kg/K]
    real(RP), intent(out) :: CVtot !> total specific heat    [J/kg/K]
    real(RP), intent(out) :: CPtot !> total specific heat    [J/kg/K]

    integer :: iqw
    !---------------------------------------------------------------------------

    qdry  = 1.0_RP
    Rtot  = 0.0_RP
    CVtot = 0.0_RP
    CPtot = 0.0_RP
    do iqw = 1, QA
       qdry  = qdry  - q(iqw) * Mq(iqw)
       Rtot  = Rtot  + q(iqw) * Rq(iqw)
       CVtot = CVtot + q(iqw) * CVq(iqw)
       CPtot = CPtot + q(iqw) * CPq(iqw)
    enddo
    Rtot  = Rtot  + qdry * Rdry
    CVtot = CVtot + qdry * CVdry
    CPtot = CPtot + qdry * CPdry

    return
  end subroutine ATMOS_THERMODYN_specific_heat_0D

  !-----------------------------------------------------------------------------
  !> ATMOS_THERMODYN_specific_heat_1D
  !! calc specific heat: qdry, Rtot, CVtot, CPtot
  !<
!OCL SERIAL
  subroutine ATMOS_THERMODYN_specific_heat_1D( &
       KA, KS, KE, &       
       QA,         &
       q, Mq, Rq, CVq, CPq,     &
       Qdry, Rtot, CVtot, CPtot )
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q(KA,QA)  !< mass concentration      [kg/kg]
    real(RP), intent(in)  :: Mq (QA)   !< mass of tracer          0 or 1
    real(RP), intent(in)  :: Rq (QA)   !< gas constant of tracer  [J/kg/K]
    real(RP), intent(in)  :: CVq(QA)   !< specific heat of tracer [J/kg/K]
    real(RP), intent(in)  :: CPq(QA)   !< specific heat of tracer [J/kg/K]

    real(RP), intent(out) :: Qdry (KA) !> dry mass concentration [kg/kg]
    real(RP), intent(out) :: Rtot (KA) !> total gas constant     [J/kg/K]
    real(RP), intent(out) :: CVtot(KA) !> total specific heat    [J/kg/K]
    real(RP), intent(out) :: CPtot(KA) !> total specific heat    [J/kg/K]

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_THERMODYN_specific_heat( QA, &
                                           q(k,:), Mq(:), Rq(:), CVq(:), CPq(:), & ! [IN]
                                           Qdry(k), Rtot(k), CVtot(k), CPtot(k)  ) ! [OUT]
    end do

    return
  end subroutine ATMOS_THERMODYN_specific_heat_1D

  !-----------------------------------------------------------------------------
  !> ATMOS_THERMODYN_specific_heat_3D
  !! calc specific heat: qdry, Rtot, CVtot, CPtot
  !<
  subroutine ATMOS_THERMODYN_specific_heat_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, QA, &
       q,                       &
       Mq, Rq, CVq, CPq,        &
       Qdry, Rtot, CVtot, CPtot )
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: QA

    real(RP), intent(in)  :: q(KA,IA,JA,QA) !< mass concentration      [kg/kg]
    real(RP), intent(in)  :: Mq (QA)        !< mass of tracer          0 or 1
    real(RP), intent(in)  :: Rq (QA)        !< gas constant of tracer  [J/kg/K]
    real(RP), intent(in)  :: CVq(QA)        !< specific heat of tracer [J/kg/K]
    real(RP), intent(in)  :: CPq(QA)        !< specific heat of tracer [J/kg/K]

    real(RP), intent(out) :: Qdry (KA,IA,JA) !> dry mass concentration [kg/kg]
    real(RP), intent(out) :: Rtot (KA,IA,JA) !> total gas constant     [J/kg/K]
    real(RP), intent(out) :: CVtot(KA,IA,JA) !> total specific heat    [J/kg/K]
    real(RP), intent(out) :: CPtot(KA,IA,JA) !> total specific heat    [J/kg/K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_specific_heat( QA, q(k,i,j,:), Mq(:), Rq(:), CVq(:), CPq(:),        & ! [IN]
                                           Qdry(k,i,j), Rtot(k,i,j), CVtot(k,i,j), CPtot(k,i,j) ) ! [OUT]
    end do
    end do
    end do

    return
  end subroutine ATMOS_THERMODYN_specific_heat_3D

  !-----------------------------------------------------------------------------
  !> calc rhot -> pres (0D)
  subroutine ATMOS_THERMODYN_rhot2pres_0D( &
       rhot, CVtot, CPtot, Rtot, &
       pres                )
    implicit none
    real(RP), intent(in) :: rhot
    real(RP), intent(in) :: CVtot
    real(RP), intent(in) :: CPtot
    real(RP), intent(in) :: Rtot

    real(RP), intent(out) :: pres
    !---------------------------------------------------------------------------

    pres = PRE00 * ( rhot * Rtot / PRE00 )**(CPtot/CVtot)

    return
  end subroutine ATMOS_THERMODYN_rhot2pres_0D

  !-----------------------------------------------------------------------------
  !> calc rhot -> pres (3D)
  subroutine ATMOS_THERMODYN_rhot2pres_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       rhot, CVtot, CPtot, Rtot, &
       pres                )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: rhot (KA,IA,JA)
    real(RP), intent(in) :: CVtot(KA,IA,JA)
    real(RP), intent(in) :: CPtot(KA,IA,JA)
    real(RP), intent(in) :: Rtot (KA,IA,JA)

    real(RP), intent(out) :: pres(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE, &
    !$omp        rhot,CVtot,CPtot,Rtot,pres)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_rhot2pres( rtot(k,i,j), CVtot(k,i,j), CPtot(k,i,j), Rtot(k,i,j), pres(k,i,j) )
    end do
    end do
    end do

    return
  end subroutine ATMOS_THERMODYN_rhot2pres_3D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> rho * ein (0D)
  subroutine ATMOS_THERMODYN_rhot2rhoe_0D( &
       rhot, CVtot, CPtot, Rtot, &
       rhoe               )
    implicit none
    real(RP), intent(in)  :: rhot  !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: CVtot !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: CPtot !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: Rtot  !< gas constant                    [J/kg/K]

    real(RP), intent(out) :: rhoe  !< density * internal energy       [J/m3]

    real(RP) :: pres
    !---------------------------------------------------------------------------

    call ATMOS_THERMODYN_rhot2pres( rhot, CVtot, CPtot, Rtot, pres )

    rhoe = pres / Rtot * CVtot

    return
  end subroutine ATMOS_THERMODYN_rhot2rhoe_0D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> rho * ein (3D)
  subroutine ATMOS_THERMODYN_rhot2rhoe_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       rhot, CVtot, CPtot, Rtot, &
       rhoe              )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: rhot (KA,IA,JA) !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: CVtot(KA,IA,JA) !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: CPtot(KA,IA,JA) !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: Rtot (KA,IA,JA) !< gas constant                    [J/kg/K]

    real(RP), intent(out) :: rhoe(KA,IA,JA)  !< density * internal energy       [J/m3]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE, &
    !$omp        rhot,CVtot,CPtot,Rtot,rhoe)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_rhot2rhoe( rhot(k,i,j), CVtot(k,i,j), CPtot(k,i,j), Rtot(k,i,j), rhoe(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_rhot2rhoe_3D

  !-----------------------------------------------------------------------------
  !> calc rho * ein -> rho * pott (0D)
  subroutine ATMOS_THERMODYN_rhoe2rhot_0D( &
       rhoe, CVtot, CPtot, Rtot, &
       rhot                      )
    implicit none
    real(RP), intent(in)  :: rhoe  !< density * internal energy       [J/m3]
    real(RP), intent(in)  :: CVtot !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: CPtot !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: Rtot  !< gas constant                    [J/kg/K]

    real(RP), intent(out) :: rhot  !< density * potential temperature [kg/m3*K]

    real(RP) :: pres
    !---------------------------------------------------------------------------

    pres = rhoe * Rtot / CVtot

    rhot = rhoe / CVtot * ( PRE00 / pres )**(Rtot/CPtot)

    return
  end subroutine ATMOS_THERMODYN_rhoe2rhot_0D

  !-----------------------------------------------------------------------------
  !> calc rho * ein -> rho * pott (3D)
  subroutine ATMOS_THERMODYN_rhoe2rhot_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       rhoe, CVtot, CPtot, Rtot, &
       rhot              )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: rhoe (KA,IA,JA) !< density * internal energy       [J/m3]
    real(RP), intent(in)  :: CVtot(KA,IA,JA) !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: CPtot(KA,IA,JA) !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: Rtot (KA,IA,JA) !< gas constant                    [J/kg/K]

    real(RP), intent(out) :: rhot(KA,IA,JA)  !< density * potential temperature [kg/m3*K]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,&
    !$omp        rhoe,CVtot,CPtot,Rtot,rhot)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_rhoe2rhot( rhoe(k,i,j), CVtot(k,i,j), CPtot(k,i,j), Rtot(k,i,j), rhot(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_rhoe2rhot_3D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> temp & pres (0D)
  subroutine ATMOS_THERMODYN_rhot2temp_pres_0D( &
       dens, rhot, Rtot, CVtot, CPtot, &
       temp, pres                      )
    implicit none
    real(RP), intent(in)  :: dens  !< density                         [kg/m3]
    real(RP), intent(in)  :: rhot  !< density * potential temperature [kg/m3*K]
    real(RP), intent(in)  :: Rtot  !< gass constant                   [kg/kg]
    real(RP), intent(in)  :: CVtot !< specific heat                   [J/kg/K]
    real(RP), intent(in)  :: CPtot !< specific heat                   [J/kg/K]

    real(RP), intent(out) :: temp  !< temperature                     [K]
    real(RP), intent(out) :: pres  !< pressure                        [Pa]
    !---------------------------------------------------------------------------

    pres = PRE00 * ( rhot * Rtot / PRE00 )**(CPtot/CVtot)
    temp = pres / ( dens * Rtot )

    return
  end subroutine ATMOS_THERMODYN_rhot2temp_pres_0D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> temp & pres (1D)
!OCL SERIAL
  subroutine ATMOS_THERMODYN_rhot2temp_pres_1D( &
       KA, KS, KE, &
       dens, rhot, Rtot, CVtot, CPtot, &
       temp, pres                      )
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: dens (KA) !< density              [kg/m3]
    real(RP), intent(in)  :: rhot (KA) !< density * pot. temp. [kg/m3*K]
    real(RP), intent(in)  :: Rtot (KA) !< mass concentration   [kg/kg]
    real(RP), intent(in)  :: CVtot(KA) !< specific heat        [J/kg/K]
    real(RP), intent(in)  :: CPtot(KA) !< specific heat        [J/kg/K]

    real(RP), intent(out) :: temp(KA)  !< temperature          [K]
    real(RP), intent(out) :: pres(KA)  !< pressure             [Pa]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_THERMODYN_rhot2temp_pres( dens(k), rhot(k), Rtot(k), CVtot(k), CPtot(k), &
                                            temp(k), pres(k)                               )
    enddo

    return
  end subroutine ATMOS_THERMODYN_rhot2temp_pres_1D

  !-----------------------------------------------------------------------------
  !> calc rho * pott -> temp & pres (3D)
  subroutine ATMOS_THERMODYN_rhot2temp_pres_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dens, rhot, Rtot, CVtot, CPtot, &
       temp, pres                      )
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: dens (KA,IA,JA) !< density              [kg/m3]
    real(RP), intent(in)  :: rhot (KA,IA,JA) !< density * pot. temp. [kg/m3*K]
    real(RP), intent(in)  :: Rtot (KA,IA,JA) !< mass concentration   [kg/kg]
    real(RP), intent(in)  :: CVtot(KA,IA,JA) !< specific heat        [J/kg/K]
    real(RP), intent(in)  :: CPtot(KA,IA,JA) !< specific heat        [J/kg/K]

    real(RP), intent(out) :: temp(KA,IA,JA)  !< temperature          [K]
    real(RP), intent(out) :: pres(KA,IA,JA)  !< pressure             [Pa]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k) &
    !$omp shared(KS,KE,IS,IE,JS,JE) &
    !$omp shared(temp,pres,dens,rhot,Rtot,CVtot,CPtot)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_rhot2temp_pres( dens(k,i,j), rhot(k,i,j), Rtot(k,i,j), CVtot(k,i,j), CPtot(k,i,j), &
                                            temp(k,i,j), pres(k,i,j)                                           )
    enddo
    enddo
    enddo
    return
  end subroutine ATMOS_THERMODYN_rhot2temp_pres_3D

  !-----------------------------------------------------------------------------
  !> calc rho * ein -> temp & pres (0D)
!OCL SERIAL
!OCL NOSIMD
  subroutine ATMOS_THERMODYN_rhoe2temp_pres_0D( &
       dens, rhoe, CVtot, Rtot, &
       temp, pres               )
    implicit none
    real(RP), intent(in)  :: dens  !< density                   [kg/m3]
    real(RP), intent(in)  :: rhoe  !< density * internal energy [J/m3]
    real(RP), intent(in)  :: CVtot !< specific heat             [J/kg/K]
    real(RP), intent(in)  :: Rtot  !< gas constant              [J/kg/K]

    real(RP), intent(out) :: temp  !< temperature               [K]
    real(RP), intent(out) :: pres  !< pressure                  [Pa]
    !---------------------------------------------------------------------------

    temp = rhoe / ( dens * CVtot )
    pres = dens * Rtot * temp

    return
  end subroutine ATMOS_THERMODYN_rhoe2temp_pres_0D

  !-----------------------------------------------------------------------------
  !> calc rho * ein -> temp & pres (3D)
  subroutine ATMOS_THERMODYN_rhoe2temp_pres_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dens, rhoe, CVtot, Rtot, &
       temp, pres                    )
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: dens (KA,IA,JA) !< density                   [kg/m3]
    real(RP), intent(in)  :: rhoe (KA,IA,JA) !< density * internal energy [J/m3]
    real(RP), intent(in)  :: CVtot(KA,IA,JA) !< specific heat             [J/kg/K]
    real(RP), intent(in)  :: Rtot (KA,IA,JA) !< gas constant              [J/kg/K]

    real(RP), intent(out) :: temp(KA,IA,JA)  !< temperature               [K]
    real(RP), intent(out) :: pres(KA,IA,JA)  !< pressure                  [Pa]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,&
    !$omp        dens,rhoe,CVtot,Rtot,temp,pres)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_rhoe2temp_pres( dens(k,i,j), rhoe(k,i,j), CVtot(k,i,j), Rtot(k,i,j), temp(k,i,j), pres(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_rhoe2temp_pres_3D

  !-----------------------------------------------------------------------------
  !> calc ein -> temp & pres (0D)
  subroutine ATMOS_THERMODYN_ein2temp_pres_0D( &
       Ein,  dens, CVtot, Rtot, &
       temp, pres               )
    implicit none
    real(RP), intent(in)  :: Ein   ! internal energy
    real(RP), intent(in)  :: dens  ! density
    real(RP), intent(in)  :: CVtot ! specific heat
    real(RP), intent(in)  :: Rtot  ! gas constant

    real(RP), intent(out) :: temp  ! temperature
    real(RP), intent(out) :: pres  ! pressure
    !---------------------------------------------------------------------------

    temp = Ein / CVtot
    pres = dens * Rtot * temp

    return
  end subroutine ATMOS_THERMODYN_ein2temp_pres_0D

  !-----------------------------------------------------------------------------
  !> calc ein -> temp & pres (3D)
  subroutine ATMOS_THERMODYN_ein2temp_pres_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Ein,  dens, CVtot, Rtot, &
       temp, pres                    )
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: Ein  (KA,IA,JA) ! internal energy
    real(RP), intent(in)  :: dens (KA,IA,JA) ! density
    real(RP), intent(in)  :: CVtot(KA,IA,JA) ! specific heat
    real(RP), intent(in)  :: Rtot (KA,IA,JA) ! gas constant

    real(RP), intent(out) :: temp(KA,IA,JA)  ! temperature
    real(RP), intent(out) :: pres(KA,IA,JA)  ! pressure

    integer  :: i, j, k
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_ein2temp_pres( Ein(k,i,j), dens(k,i,j), CVtot(k,i,j), Rtot(k,i,j), temp(k,i,j), pres(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_ein2temp_pres_3D

  !-----------------------------------------------------------------------------
  !> calc pott -> temp & pres (0D)
  subroutine ATMOS_THERMODYN_pott2temp_pres_0D( &
       dens, pott, CVtot, CPtot, Rtot, &
       temp, pres                      )
    implicit none

    real(RP), intent(in)  :: dens  ! density
    real(RP), intent(in)  :: pott  ! potential temperature
    real(RP), intent(in)  :: CVtot ! specific heat
    real(RP), intent(in)  :: CPtot ! specific heat
    real(RP), intent(in)  :: Rtot  ! gas constant

    real(RP), intent(out) :: temp  ! temperature
    real(RP), intent(out) :: pres  ! pressure
    !---------------------------------------------------------------------------

    pres = PRE00 * ( dens * pott * Rtot / PRE00 )**(CPtot/CVtot)
    temp = pres / ( dens * Rtot )

    return
  end subroutine ATMOS_THERMODYN_pott2temp_pres_0D

  !-----------------------------------------------------------------------------
  !> calc pott -> temp & pres (3D)
  subroutine ATMOS_THERMODYN_pott2temp_pres_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dens, pott, CVtot, CPtot, Rtot, &
       temp, pres                    )
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: dens (KA,IA,JA) ! density
    real(RP), intent(in)  :: pott (KA,IA,JA) ! potential temperature
    real(RP), intent(in)  :: CVtot(KA,IA,JA) ! specific heat
    real(RP), intent(in)  :: CPtot(KA,IA,JA) ! specific heat
    real(RP), intent(in)  :: Rtot (KA,IA,JA) ! gas constant

    real(RP), intent(out) :: temp(KA,IA,JA)    ! temperature
    real(RP), intent(out) :: pres(KA,IA,JA)    ! pressure

    integer  :: i, j, k
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_pott2temp_pres( dens(k,i,j), pott(k,i,j), CVtot(k,i,j), CPtot(k,i,j), Rtot(k,i,j), temp(k,i,j), pres(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_pott2temp_pres_3D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_temp_pres2pott_0D( &
       temp, pres, CPtot, Rtot, &
       pott                     )
    implicit none
    real(RP), intent(in)  :: temp  ! temperature           [K]
    real(RP), intent(in)  :: pres  ! pressure              [Pa]
    real(RP), intent(in)  :: CPtot ! specific heat         [J/kg/K]
    real(RP), intent(in)  :: Rtot  ! gas constant          [J/kg/K]

    real(RP), intent(out) :: pott  ! potential temperature [K]
    !---------------------------------------------------------------------------

    pott = temp * ( PRE00 / pres )**(Rtot/CPtot)

    return
  end subroutine ATMOS_THERMODYN_temp_pres2pott_0D

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_THERMODYN_temp_pres2pott_1D( &
       KA, KS, KE, &
       temp, pres, CPtot, Rtot, &
       pott              )
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: temp (KA) ! temperature           [K]
    real(RP), intent(in)  :: pres (KA) ! pressure              [Pa]
    real(RP), intent(in)  :: CPtot(KA) ! specific heat         [J/kg/K]
    real(RP), intent(in)  :: Rtot (KA) ! gas constant          [J/kg/K]

    real(RP), intent(out) :: pott(KA)  ! potential temperature [K]

    integer  :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       call ATMOS_THERMODYN_temp_pres2pott( temp(k), pres(k), CPtot(k), Rtot(k), pott(k) )
    enddo

    return
  end subroutine ATMOS_THERMODYN_temp_pres2pott_1D

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_temp_pres2pott_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       temp, pres, CPtot, Rtot, &
       pott              )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: temp (KA,IA,JA) ! temperature           [K]
    real(RP), intent(in)  :: pres (KA,IA,JA) ! pressure              [Pa]
    real(RP), intent(in)  :: CPtot(KA,IA,JA) ! specific heat         [J/kg/K]
    real(RP), intent(in)  :: Rtot (KA,IA,JA) ! gas constant          [J/kg/K]

    real(RP), intent(out) :: pott(KA,IA,JA)  ! potential temperature [K]

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       call ATMOS_THERMODYN_temp_pres2pott( temp(k,i,j), pres(k,i,j), CPtot(k,i,j), Rtot(k,i,j), pott(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_THERMODYN_temp_pres2pott_3D

end module scale_atmos_thermodyn
