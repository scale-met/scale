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
       flux_rain, &
       flux_snow, &
       velw,      &
       dt         )
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
       CZ   => GRID_CZ,   &
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

    real(8), intent(out)   :: flux_rain(IA,JA)
    real(8), intent(out)   :: flux_snow(IA,JA)
    real(8), intent(in)    :: tem      (KA,IA,JA)
    real(8), intent(in)    :: velw     (KA,IA,JA,QA)
    real(8), intent(in)    :: dt

    real(8) :: qflx(KA,IA,JA,QA)
    real(8) :: eflx(KA,IA,JA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('precipitation')
#ifdef _FPCOLL_
call START_COLLECTION("precipitation")
#endif

    do iq = I_QC, I_NG
    do j  = JS, JE
    do i  = IS, IE
       qflx(KE,i,j,iq) = 0.D0
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       eflx(KE,i,j) = 0.D0
    enddo
    enddo

    ! tracer transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative
    do iq = I_QC, I_NG
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS-1, KE-1
       qflx(k,i,j,iq) = velw(k+1,i,j,iq) * var(k+1,i,j,I_DENS) * var(k+1,i,j,5+iq)
    enddo
    enddo
    enddo
    enddo

    do iq = I_QC, I_QG

       !--- internal energy
       do j = JS, JE
       do i = IS, IE
       do k = KS-1, KE-1
          eflx(k,i,j) = velw(k+1,i,j,iq) * rhogq(k+1,i,j,iq) * CVw(iq) * tem(k+1,i,j)
       enddo
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          rhoge(k,i,j) = rhoge(k,i,j) - dt * ( ( eflx(k,i,j)-eflx(ij,k-1) ) ) * RCDZ(k)
       enddo
       enddo

       !--- potential energy
       do k  = KS-1, KE-1
       do ij = IJS,  IJE
          eflx(k,i,j) = velw(k+1,i,j,iq) * rhogq(k+1,i,j,iq) * GRAV * CZ(k+1)
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          rhoge(k,i,j) = rhoge(k,i,j) - dt * ( ( eflx(k,i,j)   -eflx(ij,k-1)    )                 &
                                           - ( qflx(k,i,j,iq)-qflx(ij,k-1,iq) ) * GRAV * CZ(k) ) * RCDZ(k)
       enddo
       enddo

       !--- momentum z (half level)
       do k  = KS-1, KE-2
       do ij = IJS,  IJE
          eflx(k,i,j) = 0.5D0 * ( velw(k+1,i,j,iq) + velw(k,i,j,iq) ) * var(k,i,j,I_MOMZ) &
                     * 0.5D0 * ( q(k+1,i,j,iq)    + q(k,i,j,iq)    )
       enddo
       enddo

       do k  = KS,  KE-1
       do ij = IJS, IJE
          var(k,i,j,I_MOMZ) = var(k,i,j,I_MOMZ) - dt * ( eflx(k+1,i,j)-eflx(k,i,j) ) * RFDZ(k)
       enddo
       enddo

       !--- momentum x
       do k  = KS-1, KE-1
       do ij = IJS,  IJE
          eflx(k,i,j) = velw(k+1,i,j,iq) * var(k+1,i,j,I_MOMX) * q(k+1,i,j,iq)
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          var(k,i,j,I_MOMX) = var(k,i,j,I_MOMX) - dt * ( eflx(k,i,j)-eflx(ij,k-1) ) * RCDZ(k)
       enddo
       enddo

       !--- momentum y
       do k  = KS-1, KE-1
       do ij = IJS,  IJE
          eflx(k,i,j) = velw(k+1,i,j,iq) * var(k+1,i,j,I_MOMY) * q(k+1,i,j,iq)
       enddo
       enddo

       do k  = KS,  KE
       do ij = IJS, IJE
          var(k,i,j,I_MOMY) = var(k,i,j,I_MOMY) - dt * ( eflx(k,i,j)-eflx(ij,k-1) ) * RCDZ(k)
       enddo
       enddo

       !--- update total density
       do k  = KS,  KE
       do ij = IJS, IJE
          rhog(k,i,j) = rhog(k,i,j) - dt * ( qflx(k,i,j,iq)-qflx(ij,k-1,iq) ) * RCDZ(k)
       enddo
       enddo

    enddo

    !--- update tracer
    do iq = I_QC, I_NG
    do k  = KS,  KE
    do ij = IJS, IJE
       rhogq(k,i,j,iq) = rhogq(k,i,j,iq) - dt * ( qflx(k,i,j,iq)-qflx(ij,k-1,iq) ) * RCDZ(k)
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
call STOP_COLLECTION("precipitation")
#endif
    call TIME_rapend  ('precipitation')

    return
  end subroutine precipitation

end module mod_precipitation
