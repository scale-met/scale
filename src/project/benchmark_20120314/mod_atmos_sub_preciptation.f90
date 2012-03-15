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
     IO_FID_LOG, &
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
       rhoq,      &
       rhoe,      &
       temp,      &
       dt         )
    use mod_const, only : &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       RovCP => CONST_RovCP, &
       Rvap  => CONST_Rvap,  &
       PRE00 => CONST_PRE00
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
       RFDZ => GRID_RFDZ
    use mod_atmos_vars, only: &
       var => atmos_var, &
       QA  => A_QA,  &
       CVw => A_CVw, &
       I_DENS, &
       I_MOMZ, &
       I_MOMX, &
       I_MOMY, &
       I_RHOT, &
       I_QV, &
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

    real(8), intent(out)   :: flux_rain(KA,IA,JA)
    real(8), intent(out)   :: flux_snow(KA,IA,JA)
    real(8), intent(in)    :: velw     (KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(8), intent(in)    :: rhoq     (KA,IA,JA,QA) ! rho * q
    real(8), intent(in)    :: rhoe     (KA,IA,JA)
    real(8), intent(in)    :: temp     (KA,IA,JA)
    real(8), intent(in)    :: dt

    real(8) :: qflx    (KA,QA)
    real(8) :: eflx    (KA)
    real(8) :: rhoe_new(KA)
    real(8) :: qdry(KA), CVmoist(KA), Rmoist(KA), pres(KA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('precipitation')
#ifdef _FPCOLL_
call START_COLLECTION("precipitation")
#endif

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative
    !OCL NORECURRENCE,PARALLEL
    do j  = JS, JE
    do i  = IS, IE

       !OCL XFILL
       do iq = I_QC, I_NG
          qflx(KE,iq) = 0.D0
       enddo
       eflx(KE) = 0.D0

       !--- tracer
       do iq = I_QC, I_NG
       do k  = KS-1, KE-1
          qflx(k,iq) = velw(k+1,i,j,iq) * rhoq(k+1,i,j,iq)
       enddo
       enddo ! tracer loop

       !--- lowermost flux is saved for land process
       do k  = KS-1, KE
          flux_rain(k,i,j) = - ( qflx(k,I_QC) &
                               + qflx(k,I_QR) )
          flux_snow(k,i,j) = - ( qflx(k,I_QI) &
                               + qflx(k,I_QS) &
                               + qflx(k,I_QG) )
       enddo 

       !OCL XFILL
       do k = KS, KE
          rhoe_new(k) = 0.D0
       enddo

       do iq = I_QC, I_QG

          !--- internal energy
          do k  = KS-1, KE-1
             eflx(k) = qflx(k,iq) * temp(k+1,i,j) * CVw(iq)
          enddo
          do k  = KS,   KE
             rhoe_new(k) = rhoe(k,i,j) - dt * ( eflx(k)-eflx(k-1) ) * RCDZ(k)
          enddo

          !--- potential energy
          do k  = KS-1, KE-1
             eflx(k) = qflx(k,iq) * GRAV * CZ(k+1)
          enddo
          do k  = KS,   KE
             rhoe_new(k) = rhoe_new(k) - dt * ( ( eflx(k)   -eflx(k-1)    ) &
                                              - ( qflx(k,iq)-qflx(k-1,iq) ) * GRAV * CZ(k) ) * RCDZ(k)
          enddo

          !--- momentum z (half level)
          do k  = KS-1, KE-2
             eflx(k) = 0.5D0 * ( velw(k+1,i,j,iq)   + velw(k,i,j,iq)   ) * var(k,i,j,I_MOMZ) &
                     * 0.5D0 * ( var (k+1,i,j,5+iq) + var (k,i,j,5+iq) )
          enddo
          do k  = KS,  KE-1
             var(k,i,j,I_MOMZ) = var(k,i,j,I_MOMZ) - dt * ( eflx(k+1)-eflx(k) ) * RFDZ(k)
          enddo

          !--- momentum x
          do k  = KS-1, KE-1
             eflx(k) = velw(k+1,i,j,iq) * var(k+1,i,j,I_MOMX) * var(k+1,i,j,5+iq)
          enddo
          do k  = KS,  KE
             var(k,i,j,I_MOMX) = var(k,i,j,I_MOMX) - dt * ( eflx(k)-eflx(k-1) ) * RCDZ(k)
          enddo

          !--- momentum y
          do k  = KS-1, KE-1
             eflx(k) = velw(k+1,i,j,iq) * var(k+1,i,j,I_MOMY) * var(k+1,i,j,5+iq)
          enddo
          do k  = KS,  KE
             var(k,i,j,I_MOMY) = var(k,i,j,I_MOMY) - dt * ( eflx(k)-eflx(k-1) ) * RCDZ(k)
          enddo

          !--- update total density
          do k  = KS,  KE
             var(k,i,j,I_DENS) = var(k,i,j,I_DENS) - dt * ( qflx(k,iq)-qflx(k-1,iq) ) * RCDZ(k)
          enddo

       enddo ! QC-QG loop

       !--- update tracer
       do iq = I_QC, I_NG
       do k  = KS, KE
          var(k,i,j,5+iq) = var(k,i,j,5+iq) - dt * ( qflx(k,iq)-qflx(k-1,iq) ) * RCDZ(k) / var(k,i,j,I_DENS)
       enddo 
       enddo ! QC-QG,NC-NG loop

       do k  = KS,  KE
          qdry(k) = 1.D0
          do iq = I_QV, I_QG
             qdry(k) = qdry(k) - var(k,i,j,5+iq)
          enddo
       enddo

       do k  = KS,  KE
          CVmoist(k) = qdry(k) * CVdry
          do iq = I_QV, I_QG
             CVmoist(k) = CVmoist(k) + var(k,i,j,5+iq) * CVw(iq)
          enddo
       enddo

       do k  = KS,  KE
          Rmoist(k) = qdry(k)*Rdry + var(k,i,j,5+I_QV)*Rvap
          pres  (k) = rhoe_new(k) * Rmoist(k) / CVmoist(k)
       enddo

       do k  = KS,  KE
          var(k,i,j,I_RHOT) = rhoe_new(k) / CVmoist(k) * ( PRE00 / pres(k) )**RovCP
       enddo

    enddo ! I loop
    enddo ! J loop

#ifdef _FPCOLL_
call STOP_COLLECTION("precipitation")
#endif
    call TIME_rapend  ('precipitation')

    return
  end subroutine precipitation

end module mod_precipitation
