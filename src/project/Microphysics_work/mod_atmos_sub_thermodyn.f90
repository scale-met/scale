!-------------------------------------------------------------------------------
!> module Thermodynamics
!!
!! @par Description
!!          Thermodynamics module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)   [new] Import from NICAM
!! @li      2012-02-10 (H.Yashiro) [mod] Reconstruction
!!
!<
!-------------------------------------------------------------------------------
module mod_thrmdyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_const, only : &
     CPdry => CONST_CPdry, &
     CVdry => CONST_CVdry, &
     Rdry  => CONST_Rdry,  &
     Rvap  => CONST_Rvap,  &
     RovCP => CONST_RovCP, &
     PRE00 => CONST_PRE00
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use mod_grid, only: &
     KA => GRID_KA, &
     IA => GRID_IA, &
     JA => GRID_JA, &
     KS => GRID_KS, &
     KE => GRID_KE, &
     IS => GRID_IS, &
     IE => GRID_IE, &
     JS => GRID_JS, &
     JE => GRID_JE, &
     IJA => GRID_IJA, &
     IJS => GRID_IJS, &
     IJE => GRID_IJE
  use mod_atmos_vars, only: &
     VA  => A_VA,  &
     QA  => A_QA,  &
     I_QV,         &
     QWA => A_QWA, &
     QWS => A_QWS, &
     QWE => A_QWE, &
     CPw => A_CPw, &
     CVw => A_CVw
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: thrmdyn_qd
  public :: thrmdyn_cv
  public :: thrmdyn_cp
  public :: thrmdyn_tempre
  public :: thrmdyn_tempre2
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
  subroutine thrmdyn_qd( qdry, q  )
    implicit none

    real(8), intent(out) :: qdry(IJA,KA)    ! dry mass concentration
    real(8), intent(in)  :: q   (IJA,KA,QA) ! mass concentration

    integer :: ij, k, iqw
    !-----------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_qd")
#endif

    do k  = KS,  KE
    do ij = IJS, IJE

       qdry(ij,k) = 1.D0

       do iqw = QWS, QWE
          qdry(ij,k) = qdry(ij,k) - q(ij,k,iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_qd")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_qd

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cp( cptot, q, qdry )
    implicit none

    real(8), intent(out) :: cptot(IJA,KA)    ! total specific heat
    real(8), intent(in)  :: q    (IJA,KA,QA) ! mass concentration
    real(8), intent(in)  :: qdry (IJA,KA)    ! dry mass concentration

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_cp")
#endif

    do k  = KS,  KE
    do ij = IJS, IJE

       cptot(ij,k) = qdry(ij,k) * CPdry

       do iqw = QWS, QWE
          cptot(ij,k) = cptot(ij,k) + q(ij,k,iqw) * CPw(iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_cp")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_cp

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cv( cvtot, q, qdry )
    implicit none

    real(8), intent(out) :: cvtot(IJA,KA)    ! total specific heat
    real(8), intent(in)  :: q    (IJA,KA,QA) ! mass concentration
    real(8), intent(in)  :: qdry (IJA,KA)    ! dry mass concentration

    integer :: ij, k, iqw
    !---------------------------------------------------------------------------

    call TIME_rapstart('thrmdyn')
#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_cv")
#endif

    do k  = KS,  KE
    do ij = IJS, IJE

       cvtot(ij,k) = qdry(ij,k) * CVdry

       do iqw = QWS, QWE
          cvtot(ij,k) = cvtot(ij,k) + q(ij,k,iqw) * CVw(iqw)
       enddo

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_cv")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_cv

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempre( &
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
call START_COLLECTION("thrmdyn_tempre")
#endif

    do k  = KS,  KE
    do ij = IJS, IJE

       cv = qdry(ij,k) * CVdry
       do iqw = QWS, QWE
          cv = cv + q(ij,k,QWS+iqw-1) * CVw(iqw)
       enddo
       Rmoist = qdry(ij,k)*Rdry + q(ij,k,I_QV)*Rvap

       temp(ij,k) = Ein(ij,k) / cv

       pres(ij,k) = dens(ij,k) * Rmoist * temp(ij,k)

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempre")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_tempre

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempre2( &
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
call START_COLLECTION("thrmdyn_tempre2")
#endif

    RPRE00   = 1.D0 / PRE00
    WKAPPA   = 1.D0 / ( 1.D0 - RovCP )

    do k  = KS,  KE
    do ij = IJS, IJE
       rhoRmoist = dens(ij,k) * ( qdry(ij,k)*Rdry + q(ij,k,I_QV)*Rvap )

       temp(ij,k) = ( pott(ij,k) * ( rhoRmoist * RPRE00 )**RovCP )**WKAPPA
       pres(ij,k) = rhoRmoist * temp(ij,k)
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempre2")
#endif
    call TIME_rapend  ('thrmdyn')

    return
  end subroutine thrmdyn_tempre2

end module mod_thrmdyn
!-------------------------------------------------------------------------------
