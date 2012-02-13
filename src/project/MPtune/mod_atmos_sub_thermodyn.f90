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
module mod_atmos_thrmdyn
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
  use mod_grid, only: &
     KA => GRID_KA, &
     IA => GRID_IA, &
     JA => GRID_JA, &
     KS => GRID_KS, &
     KE => GRID_KE, &
     IS => GRID_IS, &
     IE => GRID_IE, &
     JS => GRID_JS, &
     JE => GRID_JE
  use mod_atmos_vars, only: &
     VA  => A_VA,  &
     I_QV,         &
     QWA => A_QWA, &
     QWS => A_QWS, &
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

    real(8), intent(out) :: qdry(KA,IA,JA)    ! dry mass concentration
    real(8), intent(in)  :: q   (KA,IA,JA,VA) ! mass concentration

    integer :: i, j, k, iqw
    !-----------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_qd")
#endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       qdry(k,i,j) = 1.D0

       do iqw = 1, QWA
          qdry(k,i,j) = qdry(k,i,j) - q(k,i,j,QWS+iqw-1)
       enddo

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_qd")
#endif

    return
  end subroutine thrmdyn_qd

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cp( cptot, q, qdry )
    implicit none

    real(8), intent(out) :: cptot(KA,IA,JA)    ! total specific heat
    real(8), intent(in)  :: q    (KA,IA,JA,VA) ! mass concentration
    real(8), intent(in)  :: qdry (KA,IA,JA)    ! dry mass concentration

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_cp")
#endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       cptot(k,i,j) = qdry(k,i,j) * CPdry

       do iqw = 1, QWA
          cptot(k,i,j) = cptot(k,i,j) + q(k,i,j,QWS+iqw-1) * CPw(iqw)
       enddo

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_cp")
#endif

    return
  end subroutine thrmdyn_cp

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_cv( cvtot, q, qdry )
    implicit none

    real(8), intent(out) :: cvtot(KA,IA,JA)    ! total specific heat
    real(8), intent(in)  :: q    (KA,IA,JA,VA) ! mass concentration
    real(8), intent(in)  :: qdry (KA,IA,JA)    ! dry mass concentration

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_cv")
#endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       cvtot(k,i,j) = qdry(k,i,j) * CVdry

       do iqw = 1, QWA
          cvtot(k,i,j) = cvtot(k,i,j) + q(k,i,j,QWS+iqw-1) * CVw(iqw)
       enddo

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_cv")
#endif

    return
  end subroutine thrmdyn_cv

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempre( &
      temp, pres,         &
      Ein,  dens, qdry, q )
    implicit none

    real(8), intent(out) :: temp(KA,IA,JA)    ! temperature
    real(8), intent(out) :: pres(KA,IA,JA)    ! pressure
    real(8), intent(in)  :: Ein (KA,IA,JA)    ! internal energy
    real(8), intent(in)  :: dens(KA,IA,JA)    ! density
    real(8), intent(in)  :: qdry(KA,IA,JA)    ! dry concentration
    real(8), intent(in)  :: q   (KA,IA,JA,VA) ! water concentration 

    real(8) :: cv

    integer :: i, j, k, iqw
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_tempre")
#endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       cv = qdry(k,i,j) * CVdry
       do iqw = 1, QWA
          cv = cv + q(k,i,j,QWS+iqw-1) * CVw(iqw)
       enddo

       temp(k,i,j) = Ein(k,i,j) / cv

       pres(k,i,j) = dens(k,i,j) * temp(k,i,j)                 &
                   * ( qdry(k,i,j)*Rdry + q(k,i,j,I_QV)*Rvap )

    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempre")
#endif

    return
  end subroutine thrmdyn_tempre

  !-----------------------------------------------------------------------------
  subroutine thrmdyn_tempre2( &
      temp, pres,         &
      dens, pott, qdry, q )
    implicit none

    real(8), intent(out) :: temp(KA,IA,JA)    ! temperature
    real(8), intent(out) :: pres(KA,IA,JA)    ! pressure
    real(8), intent(in)  :: dens(KA,IA,JA)    ! density
    real(8), intent(in)  :: pott(KA,IA,JA)    ! potential temperature
    real(8), intent(in)  :: qdry(KA,IA,JA)    ! dry concentration
    real(8), intent(in)  :: q   (KA,IA,JA,VA) ! water concentration 

    real(8) :: RPRE00, WKAPPA, rhoR

    integer :: i, j, k
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("thrmdyn_tempre2")
#endif

    RPRE00   = 1.D0 / PRE00
    WKAPPA   = 1.D0 / ( 1.D0 - RovCP )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rhoR = dens(k,i,j) * ( qdry(k,i,j)*Rdry + q(k,i,j,I_QV)*Rvap )

       temp(k,i,j) = ( pott(k,i,j) * ( rhoR * RPRE00 )**RovCP )**WKAPPA
       pres(k,i,j) = rhoR * temp(k,i,j)
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("thrmdyn_tempre2")
#endif

    return
  end subroutine thrmdyn_tempre2

end module mod_atmos_thrmdyn
!-------------------------------------------------------------------------------
