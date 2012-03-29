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
  public :: ATMOS_THRRMODYN_setup
  public :: ATMOS_THRRMODYN_qd
  public :: ATMOS_THRRMODYN_cv
  public :: ATMOS_THRRMODYN_cp
  public :: ATMOS_THRRMODYN_tempre
  public :: ATMOS_THRRMODYN_tempre2

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(8), public,      save :: AQ_CP(QQA) ! CP for each hydrometeors
  real(8), public,      save :: AQ_CV(QQA) ! CV for each hydrometeors

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
  subroutine ATMOS_THRRMODYN_setup
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
  end subroutine ATMOS_THRRMODYN_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THRRMODYN_qd( qdry, q )
    implicit none

    real(8), intent(out) :: qdry(IJA,KA)    ! dry mass concentration
    real(8), intent(in)  :: q   (IJA,KA,QA) ! mass concentration

    integer :: ij, k, iqw
    !-----------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("ATMOS_THRRMODYN_qd")
#endif

    do k  = 1, KA
    do ij = 1, IJA

       qdry(ij,k) = 1.D0

       do iqw = QQS, QQE
          qdry(ij,k) = qdry(ij,k) - q(ij,k,iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("ATMOS_THRRMODYN_qd")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine ATMOS_THRRMODYN_qd

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THRRMODYN_cp( cptot, q, qdry )
    implicit none

    real(8), intent(out) :: cptot(IJA,KA)    ! total specific heat
    real(8), intent(in)  :: q    (IJA,KA,QA) ! mass concentration
    real(8), intent(in)  :: qdry (IJA,KA)    ! dry mass concentration

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("ATMOS_THRRMODYN_cp")
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
call STOP_COLLECTION("ATMOS_THRRMODYN_cp")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine ATMOS_THRRMODYN_cp

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THRRMODYN_cv( cvtot, q, qdry )
    implicit none

    real(8), intent(out) :: cvtot(IJA,KA)    ! total specific heat
    real(8), intent(in)  :: q    (IJA,KA,QA) ! mass concentration
    real(8), intent(in)  :: qdry (IJA,KA)    ! dry mass concentration

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("ATMOS_THRRMODYN_cv")
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
call STOP_COLLECTION("ATMOS_THRRMODYN_cv")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine ATMOS_THRRMODYN_cv

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THRRMODYN_tempre( &
      temp, pres,         &
      Ein,  dens, qdry, q )
    implicit none

    real(8), intent(out) :: temp(IJA,KA)    ! temperature
    real(8), intent(out) :: pres(IJA,KA)    ! pressure
    real(8), intent(in)  :: Ein (IJA,KA)    ! internal energy
    real(8), intent(in)  :: dens(IJA,KA)    ! density
    real(8), intent(in)  :: qdry(IJA,KA)    ! dry concentration
    real(8), intent(in)  :: q   (IJA,KA,QA) ! water concentration 

    real(8) :: cv, Rmoist

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("ATMOS_THRRMODYN_tempre")
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
call STOP_COLLECTION("ATMOS_THRRMODYN_tempre")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine ATMOS_THRRMODYN_tempre

  !-----------------------------------------------------------------------------
  subroutine ATMOS_THRRMODYN_tempre2( &
      temp, pres,         &
      dens, pott, qdry, q )
    implicit none

    real(8), intent(out) :: temp(IJA,KA)    ! temperature
    real(8), intent(out) :: pres(IJA,KA)    ! pressure
    real(8), intent(in)  :: dens(IJA,KA)    ! density
    real(8), intent(in)  :: pott(IJA,KA)    ! potential temperature
    real(8), intent(in)  :: qdry(IJA,KA)    ! dry concentration
    real(8), intent(in)  :: q   (IJA,KA,QA) ! water concentration 

    real(8) :: RPRE00, WKAPPA, rhoRmoist

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("ATMOS_THRRMODYN_tempre2")
#endif

    RPRE00   = 1.D0 / PRE00
    WKAPPA   = 1.D0 / ( 1.D0 - RovCP )

    do k  = 1, KA
    do ij = 1, IJA
       rhoRmoist = dens(ij,k) * ( qdry(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )

       temp(ij,k) = ( pott(ij,k) * ( rhoRmoist * RPRE00 )**RovCP )**WKAPPA
       pres(ij,k) = rhoRmoist * temp(ij,k)
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("ATMOS_THRRMODYN_tempre2")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine ATMOS_THRRMODYN_tempre2

end module mod_atmos_thermodyn

