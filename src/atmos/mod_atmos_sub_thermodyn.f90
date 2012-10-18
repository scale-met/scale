!-------------------------------------------------------------------------------
!> module Thermodynamics
!!
!! @par Description
!!          Thermodynamics module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)    [new] Import from NICAM
!! @li      2012-02-10 (H.Yashiro)  [mod] Reconstruction
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_thermodyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  use mod_const, only : &
     Rdry  => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     RovCP => CONST_RovCP, &
     PRE00 => CONST_PRE00, &
     Rvap  => CONST_Rvap,  &
     CPvap => CONST_CPvap, &
     CVvap => CONST_CVvap, &
     CL    => CONST_CL,    &
     CI    => CONST_CI
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
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
  public :: ATMOS_THERMODYN_tempre
  public :: ATMOS_THERMODYN_tempre2

  public :: ATMOS_THERMODYN_qd_kij
  public :: ATMOS_THERMODYN_cv_kij
  public :: ATMOS_THERMODYN_cp_kij
  public :: ATMOS_THERMODYN_tempre_kij
  public :: ATMOS_THERMODYN_tempre2_kij
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public,      save :: AQ_CP(QQA) ! CP for each hydrometeors
  real(RP), public,      save :: AQ_CV(QQA) ! CV for each hydrometeors

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
  subroutine ATMOS_THERMODYN_setup
    implicit none
    !---------------------------------------------------------------------------

    AQ_CP(I_QV) = CPvap
    AQ_CV(I_QV) = CVvap

    if ( I_QC > 0 ) then
       AQ_CP(I_QC) = CL
       AQ_CV(I_QC) = CL
    endif
    if ( I_QR > 0 ) then
       AQ_CP(I_QR) = CL
       AQ_CV(I_QR) = CL
    endif
    if ( I_QI > 0 ) then
       AQ_CP(I_QI) = CI
       AQ_CV(I_QI) = CI
    endif
    if ( I_QS > 0 ) then
       AQ_CP(I_QS) = CI
       AQ_CV(I_QS) = CI
    endif
    if ( I_QG > 0 ) then
       AQ_CP(I_QG) = CI
       AQ_CV(I_QG) = CI
    endif

    return
  end subroutine ATMOS_THERMODYN_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_qd( qdry, q )
    implicit none

    real(RP), intent(out) :: qdry(IJA,KA)    ! dry mass concentration
    real(RP), intent(in)  :: q   (IJA,KA,QA) ! mass concentration

    integer :: ij, k, iqw
    !-----------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k  = 1, KA
    do ij = 1, IJA

       qdry(ij,k) = 1.0_RP

       do iqw = QQS, QQE
          qdry(ij,k) = qdry(ij,k) - q(ij,k,iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_qd
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_qd_kij( qdry, q )
    implicit none

    real(RP), intent(out) :: qdry(KA,IA,JA)    ! dry mass concentration
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) ! mass concentration

    integer :: i,j, k, iqw
    !-----------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do j = 1, JA
    do i = 1, IA
    do k  = 1, KA
       qdry(k,i,j) = 1.0_RP
       do iqw = QQS, QQE
          qdry(k,i,j) = qdry(k,i,j) - q(k,i,j,iqw)
       enddo

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_qd_kij
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_cp( cptot, q, qdry )
    implicit none

    real(RP), intent(out) :: cptot(IJA,KA)    ! total specific heat
    real(RP), intent(in)  :: q    (IJA,KA,QA) ! mass concentration
    real(RP), intent(in)  :: qdry (IJA,KA)    ! dry mass concentration

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k  = 1, KA
    do ij = 1, IJA

       cptot(ij,k) = qdry(ij,k) * CPdry

       do iqw = QQS, QQE
          cptot(ij,k) = cptot(ij,k) + q(ij,k,iqw) * AQ_CP(iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_cp
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_cp_kij( cptot, q, qdry )
    implicit none

    real(RP), intent(out) :: cptot(KA,IA,JA)    ! total specific heat
    real(RP), intent(in)  :: q    (KA,IA,JA,QA) ! mass concentration
    real(RP), intent(in)  :: qdry (KA,IA,JA)    ! dry mass concentration

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA

       cptot(k,i,j) = qdry(k,i,j) * CPdry

       do iqw = QQS, QQE
          cptot(k,i,j) = cptot(k,i,j) + q(k,i,j,iqw) * AQ_CP(iqw)
       enddo

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_cp_kij
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_cv( cvtot, q, qdry )
    implicit none

    real(RP), intent(out) :: cvtot(IJA,KA)    ! total specific heat
    real(RP), intent(in)  :: q    (IJA,KA,QA) ! mass concentration
    real(RP), intent(in)  :: qdry (IJA,KA)    ! dry mass concentration

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k  = 1, KA
    do ij = 1, IJA

       cvtot(ij,k) = qdry(ij,k) * CVdry

       do iqw = QQS, QQE
          cvtot(ij,k) = cvtot(ij,k) + q(ij,k,iqw) * AQ_CV(iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_cv
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_cv_kij( cvtot, q, qdry )
    implicit none

    real(RP), intent(out) :: cvtot(KA,IA,JA)    ! total specific heat
    real(RP), intent(in)  :: q    (KA,IA,JA,QA) ! mass concentration
    real(RP), intent(in)  :: qdry (KA,IA,JA)    ! dry mass concentration

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA

       cvtot(k,i,j) = qdry(k,i,j) * CVdry

       do iqw = QQS, QQE
          cvtot(k,i,j) = cvtot(k,i,j) + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_cv_kij
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_tempre( &
      temp, pres,         &
      Ein,  dens, qdry, q )
    implicit none

    real(RP), intent(out) :: temp(IJA,KA)    ! temperature
    real(RP), intent(out) :: pres(IJA,KA)    ! pressure
    real(RP), intent(in)  :: Ein (IJA,KA)    ! internal energy
    real(RP), intent(in)  :: dens(IJA,KA)    ! density
    real(RP), intent(in)  :: qdry(IJA,KA)    ! dry concentration
    real(RP), intent(in)  :: q   (IJA,KA,QA) ! water concentration 

    real(RP) :: cv, Rmoist

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k  = 1, KA
    do ij = 1, IJA

       cv = qdry(ij,k) * CVdry
       do iqw = QQS, QQE
          cv = cv + q(ij,k,iqw) * AQ_CV(iqw)
       enddo
       Rmoist = qdry(ij,k)*Rdry + q(ij,k,I_QV)*Rvap

       temp(ij,k) = Ein(ij,k) / cv

       pres(ij,k) = dens(ij,k) * Rmoist * temp(ij,k)

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_tempre
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_tempre_kij( &
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

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA

       cv = qdry(k,i,j) * CVdry
       do iqw = QQS, QQE
          cv = cv + q(k,i,j,iqw) * AQ_CV(iqw)
       enddo
       Rmoist = qdry(k,i,j)*Rdry + q(k,i,j,I_QV)*Rvap

       temp(k,i,j) = Ein(k,i,j) / cv

       pres(k,i,j) = dens(k,i,j) * Rmoist * temp(k,i,j)

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_tempre_kij
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_tempre2( &
      temp, pres,         &
      dens, pott, qdry, q )
    implicit none

    real(RP), intent(out) :: temp(IJA,KA)    ! temperature
    real(RP), intent(out) :: pres(IJA,KA)    ! pressure
    real(RP), intent(in)  :: dens(IJA,KA)    ! density
    real(RP), intent(in)  :: pott(IJA,KA)    ! potential temperature
    real(RP), intent(in)  :: qdry(IJA,KA)    ! dry concentration
    real(RP), intent(in)  :: q   (IJA,KA,QA) ! water concentration 

    real(RP) :: RPRE00, WKAPPA, rhoRmoist

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    RPRE00   = 1.0_RP / PRE00
    WKAPPA   = 1.0_RP / ( 1.0_RP - RovCP )

    do k  = 1, KA
    do ij = 1, IJA
       rhoRmoist = dens(ij,k) * ( qdry(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )

       temp(ij,k) = ( pott(ij,k) * ( rhoRmoist * RPRE00 )**RovCP )**WKAPPA
       pres(ij,k) = rhoRmoist * temp(ij,k)
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_tempre2
  !-----------------------------------------------------------------------------
  subroutine ATMOS_THERMODYN_tempre2_kij( &
      temp, pres,         &
      dens, pott, qdry, q )
    implicit none

    real(RP), intent(out) :: temp(KA,IA,JA)    ! temperature
    real(RP), intent(out) :: pres(KA,IA,JA)    ! pressure
    real(RP), intent(in)  :: dens(KA,IA,JA)    ! density
    real(RP), intent(in)  :: pott(KA,IA,JA)    ! potential temperature
    real(RP), intent(in)  :: qdry(KA,IA,JA)    ! dry concentration
    real(RP), intent(in)  :: q   (KA,IA,JA,QA) ! water concentration 

    real(RP) :: RPRE00, WKAPPA, rhoRmoist

    integer :: i, j, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_thermodyn')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_thermodyn')
#endif

    RPRE00   = 1.0_RP / PRE00
    WKAPPA   = 1.0_RP / ( 1.0_RP - RovCP )

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA
       rhoRmoist = dens(k,i,j) * ( qdry(k,i,j)*Rdry + q(k,i,j,I_QV)*Rvap )

       temp(k,i,j) = ( pott(k,i,j) * ( rhoRmoist * RPRE00 )**RovCP )**WKAPPA
       pres(k,i,j) = rhoRmoist * temp(k,i,j)
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_thermodyn')
#endif
    call TIME_rapend  ('SUB_thermodyn')

    return
  end subroutine ATMOS_THERMODYN_tempre2_kij
  !-----------------------------------------------------------------------------
end module mod_atmos_thermodyn

