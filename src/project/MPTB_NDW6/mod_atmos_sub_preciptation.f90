!-------------------------------------------------------------------------------
!> module Precipitation Transport
!!
!! @par Description
!!          Precipitation Transport module
!!          (1D advection)
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-08 (H.Yashiro) [mod] New
!!
!<
!-------------------------------------------------------------------------------
module mod_precipitation
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
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

  public :: precipitation

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
  !> precipitation transport
  !-----------------------------------------------------------------------------
  subroutine precipitation( &
       rhog,   &
       rhogvx, &
       rhogvy, &
       rhogw,  &
       rhoge,  &
       rhogq,  &
       tem,    &
       q,      &
       velw,   &
       precip, &
       dt      )
    use mod_const, only : &
       GRAV => CONST_GRAV
    use mod_grid, only: &
       KA   => GRID_KA,   &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       CDZ  => GRID_CDZ,  &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ, &
       IJA  => GRID_IJA,  &
       IJS  => GRID_IJS,  &
       IJE  => GRID_IJE
    use mod_atmos_vars, only: &
       QA  => A_QA,  &
       CVw => A_CVw, &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG, &
       I_NG
    use mod_thrmdyn, only : &
       thrmdyn_qd, &
       thrmdyn_tempre
    implicit none

    real(8), intent(inout) :: rhog  (IJA,KA)
    real(8), intent(inout) :: rhogvx(IJA,KA)
    real(8), intent(inout) :: rhogvy(IJA,KA)
    real(8), intent(inout) :: rhogw (IJA,KA)
    real(8), intent(inout) :: rhoge (IJA,KA)
    real(8), intent(inout) :: rhogq (IJA,KA,QA)
    real(8), intent(in)    :: tem   (IJA,KA)
    real(8), intent(in)    :: q     (IJA,KA,QA)
    real(8), intent(in)    :: velw  (IJA,KA,QA)
    real(8), intent(out)   :: precip(IJA,2)
    real(8), intent(in)    :: dt

    real(8) :: qflx(IJA,KA,QA)
    real(8) :: eflx(IJA,KA)

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

    call TIME_rapstart('precip_transport_nwater')
#ifdef _FPCOLL_
call START_COLLECTION("precip_transport_nwater")
#endif

    do nq = I_QC, I_NG
    do ij = IJS, IJE
       qflx(ij,KE,nq) = 0.D0
    enddo
    enddo
    do ij = IJS, IJE
       eflx(ij,KE) = 0.D0
    enddo

    ! tracer transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative
    do nq = I_QC, I_NG
    do k  = KS-1, KE-1
    do ij = IJS,  IJE
       qflx(ij,k,nq) = velw(ij,k+1,nq) * rhogq(ij,k+1,nq)
    enddo
    enddo
    enddo

    do nq = I_QC, I_QG

       !--- energy
       do k  = KS-1, KE-1
       do ij = IJS,  IJE
          eflx(ij,k) = velw(ij,k+1,nq) * rhogq(ij,k+1,nq) * CVw(nq) * tem(ij,k)
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          rhoge(ij,k) = rhoge(ij,k) - dt * ( ( qflx(ij,k,nq)-qflx(ij,k-1,nq) ) * GRAV * CDZ(k) &           ! potential energy
                                           + ( eflx(ij,k)   -eflx(ij,k-1)    )                 ) * RCDZ(k) ! internal energy
       enddo
       enddo

       !--- momentum z (half level)
       do k  = KS-1, KE-2
       do ij = IJS,  IJE
          eflx(ij,k) = 0.5D0 * ( velw(ij,k+1,nq) + velw(ij,k,nq) ) * rhogw(ij,k) &
                     * 0.5D0 * ( q(ij,k+1,nq)    + q(ij,k,nq)    )
       enddo
       enddo

       do k  = KS,  KE-1
       do ij = IJS, IJE
          rhogw(ij,k) = rhogw(ij,k) - dt * ( eflx(ij,k+1)-eflx(ij,k) ) * RFDZ(k)
       enddo
       enddo

       !--- momentum x
       do k  = KS-1, KE-1
       do ij = IJS,  IJE
          eflx(ij,k) = velw(ij,k+1,nq) * rhogvx(ij,k+1) * q(ij,k+1,nq)
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          rhogvx(ij,k) = rhogvx(ij,k) - dt * ( eflx(ij,k)-eflx(ij,k-1) ) * RCDZ(k)
       enddo
       enddo

       !--- momentum y
       do k  = KS-1, KE-1
       do ij = IJS,  IJE
          eflx(ij,k) = velw(ij,k+1,nq) * rhogvy(ij,k+1) * q(ij,k+1,nq)
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          rhogvy(ij,k) = rhogvy(ij,k) - dt * ( eflx(ij,k)-eflx(ij,k-1) ) * RCDZ(k)
       enddo
       enddo

       !--- update total density
       do k  = KS,  KE
       do ij = IJS, IJE
          rhog(ij,k) = rhog(ij,k) - dt * ( qflx(ij,k,nq)-qflx(ij,k-1,nq) ) * RCDZ(k)
       enddo
       enddo

    enddo

    !--- update tracer
    do nq = I_QC, I_NG
    do k  = KS,  KE
    do ij = IJS, IJE
       rhogq(ij,k,nq) = rhogq(ij,k,nq) - dt * ( qflx(ij,k,nq)-qflx(ij,k-1,nq) ) * RCDZ(k)
    enddo
    enddo
    enddo

    !--- lowermost flux is saved for land process
    do ij = IJS, IJE
       precip(ij,1) = - ( qflx(ij,KS-1,I_QC) &
                        + qflx(ij,KS-1,I_QR) )
       precip(ij,2) = - ( qflx(ij,KS-1,I_QI) &
                        + qflx(ij,KS-1,I_QS) &
                        + qflx(ij,KS-1,I_QG) )
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("precip_transport_nwater")
#endif
    call TIME_rapend  ('precip_transport_nwater')

    return
  end subroutine precipitation

end module mod_precipitation
