!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Diagnostic
!!
!! @par Description
!!          Calculate diagnostic variables
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2017-02-24 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_diagnostic
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DIAGNOSTIC_get_vel
  public :: ATMOS_DIAGNOSTIC_get_therm
  public :: ATMOS_DIAGNOSTIC_get_phyd
  public :: ATMOS_DIAGNOSTIC_get_potv
  public :: ATMOS_DIAGNOSTIC_get_teml
  public :: ATMOS_DIAGNOSTIC_get_n2

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
  !> ATMOS_DIAGNOSTIC_get_vel
  !! W, U, V
  !<
  subroutine ATMOS_DIAGNOSTIC_get_vel( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       W,    &
       U,    &
       V     )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    use scale_gridtrans, only: &
       GSQRT => GTRANS_GSQRT, &
       J13G  => GTRANS_J13G,  &
       J23G  => GTRANS_J23G,  &
       I_XYW, &
       I_XYZ
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)

    real(RP), intent(out) :: W   (KA,IA,JA)
    real(RP), intent(out) :: U   (KA,IA,JA)
    real(RP), intent(out) :: V   (KA,IA,JA)

    real(RP) :: momws

    integer :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS+1, KE-1
       W(k,i,j) = 0.5_RP * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = 2, JA
    do i = 2, IA
!       W(KS,i,j) = 0.5_RP * ( MOMZ(KS,i,j) ) / DENS(KS,i,j)

       ! at KS+1/2
       momws = MOMZ(KS,i,j) &
             + ( J13G(KS,i,j,I_XYW) * ( MOMX(KS,i,j) + MOMX(KS,i-1,j) + MOMX(KS+1,i,j) + MOMX(KS+1,i-1,j) ) &
               + J23G(KS,i,j,I_XYW) * ( MOMY(KS,i,j) + MOMY(KS,i,j-1) + MOMY(KS+1,i,j) + MOMY(KS+1,i,j-1) ) ) &
               * 0.25_RP / GSQRT(KS,i,j,I_XYW)
       ! at KS
       ! momws at the surface is assumed to be zero
       W(KS,i,j) = momws * 0.5_RP &
                 - ( J13G(KS,i,j,I_XYZ) * ( MOMX(KS,i,j) + MOMX(KS,i-1,j) ) &
                   + J23G(KS,i,j,I_XYZ) * ( MOMY(KS,i,j) + MOMY(KS,i,j-1) ) ) &
                   * 0.5_RP / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
       W(KE,i,j) = 0.5_RP * ( MOMZ(KE-1,i,j) ) / DENS(KE,i,j)
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 2, IA
    do k = KS, KE
       U(k,i,j) = 0.5_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do k = KS, KE
       U(k,1,j) = MOMX(k,1,j) / DENS(k,1,j)
    enddo
    enddo

   !OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = 2, JA
    do i = 1, IA
    do k = KS, KE
       V(k,i,j) = 0.5_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
    do i = 1, IA
    do k = KS, KE
       V(k,i,1) = MOMY(k,i,1) / DENS(k,i,1)
    enddo
    enddo
 
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = 1, JA
    do i  = 1, IA
       W(   1:KS-1,i,j) = W(KS,i,j)
       U(   1:KS-1,i,j) = U(KS,i,j)
       V(   1:KS-1,i,j) = V(KS,i,j)
       W(KE+1:KA,  i,j) = W(KE,i,j)
       U(KE+1:KA,  i,j) = U(KE,i,j)
       V(KE+1:KA,  i,j) = V(KE,i,j)
    enddo
    enddo

    call COMM_vars8( W(:,:,:), 1 )
    call COMM_vars8( U(:,:,:), 2 )
    call COMM_vars8( V(:,:,:), 3 )
    call COMM_wait ( W(:,:,:), 1, .false. )
    call COMM_wait ( U(:,:,:), 2, .false. )
    call COMM_wait ( V(:,:,:), 3, .false. )

    return
  end subroutine ATMOS_DIAGNOSTIC_get_vel

  !-----------------------------------------------------------------------------
  !> ATMOS_DIAGNOSTIC_get_therm
  !! potential temperature, temperature, pressure
  !<
  subroutine ATMOS_DIAGNOSTIC_get_therm( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       DENS, RHOT,             &
       Rtot, CVtot, CPtot,     &
       POTT, TEMP, PRES, EXNER )
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: DENS (KA,IA,JA)
    real(RP), intent(in)  :: RHOT (KA,IA,JA)
    real(RP), intent(in)  :: Rtot (KA,IA,JA)
    real(RP), intent(in)  :: CVtot(KA,IA,JA)
    real(RP), intent(in)  :: CPtot(KA,IA,JA)

    real(RP), intent(out) :: POTT (KA,IA,JA)
    real(RP), intent(out) :: TEMP (KA,IA,JA)
    real(RP), intent(out) :: PRES (KA,IA,JA)
    real(RP), intent(out) :: EXNER(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call THERMODYN_temp_pres( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              DENS(:,:,:), RHOT(:,:,:),                & ! (in)
                              Rtot(:,:,:), CVtot(:,:,:), CPtot(:,:,:), & ! (in)
                              TEMP(:,:,:), PRES(:,:,:)                 ) ! (out)
                              

!OCL XFILL
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k) &
    !$omp shared(POTT,EXNER,RHOT,DENS,TEMP) &
    !$omp shared(KS,KE,IS,IE,JS,JE)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       POTT (k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       EXNER(k,i,j) = TEMP(k,i,j) / POTT(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_DIAGNOSTIC_get_therm

  !-----------------------------------------------------------------------------
  !> ATMOS_DIAGNOSTIC_get_phyd
  !! hydrostatic pressure 
  !<
  subroutine ATMOS_DIAGNOSTIC_get_phyd( &
       KA, KS, KE, &       
       IA, IS, IE, &
       JA, JS, JE, &
       DENS, PRES, &
       CZ, FZ,     &
       PHYD        )
    use scale_const, only: &
       GRAV => CONST_GRAV
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: PRES(KA,IA,JA)

    real(RP), intent(in)  :: CZ(KA,IA,JA)
    real(RP), intent(in)  :: FZ(KA,IA,JA)

    real(RP), intent(out) :: PHYD(KA,IA,JA)

    real(RP) :: ph(KA)  !> hydrostatic pressure at the half level

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,ph) &
    !$omp shared(PHYD,DENS,PRES,CZ,FZ,GRAV) &
    !$omp shared(KS,KE,IS,IE,JS,JE)
    do j = JS, JE
    do i = IS, IE
       ph(KE) = PRES(KE,i,j) - DENS(KE,i,j) * GRAV * ( FZ(KE,i,j) - CZ(KE,i,j) )
       do k = KE, KS, -1
          ph(k-1) = ph(k) + DENS(k,i,j) * GRAV * ( FZ(k,i,j) - FZ(k-1,i,j) )
          PHYD(k,i,j) = ( ph(k) + ph(k-1) ) * 0.5_RP
       end do
    enddo
    enddo

    return
  end subroutine ATMOS_DIAGNOSTIC_get_phyd

  !-----------------------------------------------------------------------------
  !> ATMOS_DIAGNOSTIC_get_n2
  !! N^2
  !<
  subroutine ATMOS_DIAGNOSTIC_get_n2( &
       KA, KS, KE, &       
       IA, IS, IE, &
       JA, JS, JE, &
       POTT, &
       Rtot, &
       CZ,   &
       N2    )
    use scale_const, only: &
       GRAV => CONST_GRAV
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: POTT(KA,IA,JA)
    real(RP), intent(in)  :: Rtot(KA,IA,JA)

    real(RP), intent(in)  :: CZ(KA,IA,JA)

    real(RP), intent(out) :: N2  (KA,IA,JA)

    real(RP) :: RPT(KA) !> Rtot * PT (= Rdry * virtual potential temperature)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k) &
    !$omp private(RPT) &
    !$omp shared(N2,POTT,Rtot,CZ,GRAV) &
    !$omp shared(KS,KE,IS,IE,JS,JE)
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          RPT(k) = Rtot(k,i,j) * POTT(k,i,j)
       end do

       N2(KS,i,j) = GRAV * ( RPT(KS+1) - RPT(KS) ) / ( ( CZ(KS+1,i,j) - CZ(KS,i,j) ) * RPT(KS) )
       do k = KS+1,KE-1
          N2(k,i,j) = GRAV * ( RPT(k+1) - RPT(k-1) ) / ( ( CZ(k+1,i,j) - CZ(k-1,i,j) ) * RPT(k ) )
       end do
       N2(KE,i,j) = GRAV * ( RPT(KE) - RPT(KE-1) ) / ( ( CZ(KE,i,j) - CZ(KE-1,i,j) ) * RPT(KE) )
    end do
    end do

    return
  end subroutine ATMOS_DIAGNOSTIC_get_n2

  !-----------------------------------------------------------------------------
  !> ATMOS_DIAGNOSTIC_get_potv
  !! virtual potential temperature
  !<
  subroutine ATMOS_DIAGNOSTIC_get_potv( &
       KA, KS, KE, &       
       IA, IS, IE, &
       JA, JS, JE, &
       POTT, &
       Rtot, &
       POTV  )
    use scale_const, only: &
       Rdry => CONST_Rdry
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: POTT(KA,IA,JA)
    real(RP), intent(in)  :: Rtot(KA,IA,JA)

    real(RP), intent(out) :: POTV(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k) &
    !$omp shared(POTV,POTT,Rtot,Rdry) &
    !$omp shared(KS,KE,IS,IE,JS,JE)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       POTV(k,i,j) = POTT(k,i,j) * Rtot(k,i,j) / Rdry
    end do
    end do
    end do

    return
  end subroutine ATMOS_DIAGNOSTIC_get_potv

  !-----------------------------------------------------------------------------
  !> ATMOS_DIAGNOSTIC_get_teml
  !! liqued water temperature
  !<
  subroutine ATMOS_DIAGNOSTIC_get_teml( &
       KA, KS, KE, &       
       IA, IS, IE, &
       JA, JS, JE, &
       TEMP,     &
       LHV, LHS, &
       LIQ, ICE, &
       CPtot,    &
       TEML      )
    use scale_const, only: &
       Rdry => CONST_Rdry
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: TEMP (KA,IA,JA)
    real(RP), intent(in)  :: LHV  (KA,IA,JA)
    real(RP), intent(in)  :: LHS  (KA,IA,JA)
    real(RP), intent(in)  :: LIQ  (KA,IA,JA)
    real(RP), intent(in)  :: ICE  (KA,IA,JA)
    real(RP), intent(in)  :: CPtot(KA,IA,JA)

    real(RP), intent(out) :: TEML(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k) &
    !$omp shared(TEML,TEMP,LHV,LHS,LIQ,ICE,CPtot) &
    !$omp shared(KS,KE,IS,IE,JS,JE)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       TEML(k,i,j) = TEMP(k,i,j) &
                   - ( LHV(k,i,j) * LIQ(k,i,j) + LHS(k,i,j) * ICE(k,i,j) ) / CPtot(k,i,j)
    end do
    end do
    end do

    return
  end subroutine ATMOS_DIAGNOSTIC_get_teml

end module scale_atmos_diagnostic
