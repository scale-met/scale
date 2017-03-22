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
  public :: ATMOS_DIAGNOSTIC_get

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

  subroutine ATMOS_DIAGNOSTIC_get( &
       POTT, &
       TEMP, &
       PRES, &
       PHYD, &
       W,    &
       U,    &
       V,    &
       N2,   &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd, &
       THERMODYN_r         => ATMOS_THERMODYN_r,  &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    implicit none
    real(RP), intent(out) :: POTT(KA,IA,JA)
    real(RP), intent(out) :: TEMP(KA,IA,JA)
    real(RP), intent(out) :: PRES(KA,IA,JA)
    real(RP), intent(out) :: PHYD(KA,IA,JA)
    real(RP), intent(out) :: W   (KA,IA,JA)
    real(RP), intent(out) :: U   (KA,IA,JA)
    real(RP), intent(out) :: V   (KA,IA,JA)
    real(RP), intent(out) :: N2  (KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP) :: ph(KA)  !> hydrostatic pressure at the half level
    real(RP) :: RPT(KA) !> Rtot * PT (= Rdry * virtual potential temperature)
    real(RP) :: q(QA)
    real(RP) :: qdry, Rtot

    integer :: k, i, j
    integer :: iq
    !---------------------------------------------------------------------------

    call THERMODYN_temp_pres( TEMP(:,:,:),   & ! [OUT]
                              PRES(:,:,:),   & ! [OUT]
                              DENS(:,:,:),   & ! [IN]
                              RHOT(:,:,:),   & ! [IN]
                              QTRC(:,:,:,:), & ! [IN]
                              TRACER_CV(:),  & ! [IN]
                              TRACER_R(:),   & ! [IN]
                              TRACER_MASS(:) ) ! [IN]


    !$omp parallel do private(i,j,k,ph) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
       ph(KE) = PRES(KE,i,j) - DENS(KE,i,j) * GRAV * ( FZ(KE,i,j) - CZ(KE,i,j) )
       do k = KE, KS, -1
          ph(k-1) = ph(k) + DENS(k,i,j) * GRAV * ( FZ(k,i,j) - FZ(k-1,i,j) )
          PHYD(k,i,j) = ( ph(k) + ph(k-1) ) * 0.5_RP
       end do
    enddo
    enddo


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
    do j = 1, JA
    do i = 1, IA
       W(KS,i,j) = 0.5_RP * (                 MOMZ(KS,i,j) ) / DENS(KS,i,j)
    enddo
    enddo
!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
       W(KE,i,j) = 0.5_RP * ( MOMZ(KE-1,i,j)               ) / DENS(KE,i,j)
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

    call COMM_vars8( U(:,:,:), 1 )
    call COMM_vars8( V(:,:,:), 2 )
    call COMM_wait ( U(:,:,:), 1, .false. )
    call COMM_wait ( V(:,:,:), 2, .false. )

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do OMP_SCHEDULE_ collapse(2)
    !$omp private(i,j,q,rpt,qdry,rtot)
    do j = 1, JA
    do i = 1, IA

       do k = KS, KE
          do iq = 1, QA
             q(iq) = QTRC(k,i,j,iq)
          end do
          call THERMODYN_qd( qdry, q(:), TRACER_MASS(:) )
          call THERMODYN_r ( Rtot, q(:), TRACER_R(:), qdry )
          RPT(k) = Rtot * POTT(k,i,j)
       end do

       N2(KS,i,j) = GRAV * ( RPT(KS+1) - RPT(KS) ) &
                         / ( ( CZ(KS+1,i,j) - CZ(KS,i,j) ) * RPT(KS) )
       do k = KS+1,KE-1
          N2(k,i,j) = GRAV * ( RPT(k+1) - RPT(k-1) ) &
                           / ( ( CZ(k+1,i,j) - CZ(k-1,i,j) ) * RPT(k) )
       end do
       N2(KE,i,j) = GRAV * ( RPT(KE) - RPT(KE-1) ) &
                         / ( ( CZ(KE,i,j) - CZ(KE-1,i,j) ) * RPT(KE) )
    end do
    end do

    return
  end subroutine ATMOS_DIAGNOSTIC_get


end module scale_atmos_diagnostic
