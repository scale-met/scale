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
module mod_atmos_precipitation
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
  public :: ATMOS_PRECIPITATION

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
  subroutine ATMOS_PRECIPITATION( &
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
       CZ   => GRID_CZ,   &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_thermodyn, only: &
       CVw => AQ_CV
    implicit none

    real(RP), intent(out)   :: flux_rain(KA,IA,JA)
    real(RP), intent(out)   :: flux_snow(KA,IA,JA)
    real(RP), intent(in)    :: velw     (KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(RP), intent(in)    :: rhoq     (KA,IA,JA,QA) ! rho * q
    real(RP), intent(in)    :: rhoe     (KA,IA,JA)
    real(RP), intent(in)    :: temp     (KA,IA,JA)
    real(RP), intent(in)    :: dt

    real(RP) :: qflx    (KA,QA)
    real(RP) :: eflx    (KA)
    real(RP) :: rhoe_new(KA)
    real(RP) :: qdry(KA), CVmoist(KA), Rmoist(KA), pres(KA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_precipitation')
#ifdef _FPCOLL_
call START_COLLECTION('SUB_precipitation')
#endif

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative
    !OCL NORECURRENCE,PARALLEL
    do j  = JS, JE
    do i  = IS, IE

       !OCL XFILL
       do iq = I_QC, I_NG
          qflx(KE,iq) = 0.0_RP
       enddo
       eflx(KE) = 0.0_RP

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
          rhoe_new(k) = 0.E0_RP
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
             eflx(k) = 0.5E0_RP * ( velw(k+1,i,j,iq) + velw(k,i,j,iq) ) * MOMZ(k,i,j) &
                     * 0.5E0_RP * ( QTRC(k+1,i,j,iq) + QTRC(k,i,j,iq) )
          enddo
          do k  = KS,  KE-1
             MOMZ(k,i,j) = MOMZ(k,i,j) - dt * ( eflx(k+1)-eflx(k) ) * RFDZ(k)
          enddo

          !--- momentum x
          do k  = KS-1, KE-1
             eflx(k) = velw(k+1,i,j,iq) * MOMX(k+1,i,j) * QTRC(k+1,i,j,iq)
          enddo
          do k  = KS,  KE
             MOMX(k,i,j) = MOMX(k,i,j) - dt * ( eflx(k)-eflx(k-1) ) * RCDZ(k)
          enddo

          !--- momentum y
          do k  = KS-1, KE-1
             eflx(k) = velw(k+1,i,j,iq) * MOMY(k+1,i,j) * QTRC(k+1,i,j,iq)
          enddo
          do k  = KS,  KE
             MOMY(k,i,j) = MOMY(k,i,j) - dt * ( eflx(k)-eflx(k-1) ) * RCDZ(k)
          enddo

          !--- update total density
          do k  = KS,  KE
             DENS(k,i,j) = DENS(k,i,j) - dt * ( qflx(k,iq)-qflx(k-1,iq) ) * RCDZ(k)
          enddo

       enddo ! QC-QG loop

       !--- update tracer
       do iq = I_QC, I_NG
       do k  = KS, KE
          QTRC(k,i,j,iq) = QTRC(k,i,j,iq) - dt * ( qflx(k,iq)-qflx(k-1,iq) ) * RCDZ(k) / DENS(k,i,j)
       enddo 
       enddo ! QC-QG,NC-NG loop

       do k  = KS,  KE
          qdry(k) = 1.E0_RP
          do iq = I_QV, I_QG
             qdry(k) = qdry(k) - QTRC(k,i,j,iq)
          enddo
       enddo

       do k  = KS,  KE
          CVmoist(k) = qdry(k) * CVdry
          do iq = I_QV, I_QG
             CVmoist(k) = CVmoist(k) + QTRC(k,i,j,iq) * CVw(iq)
          enddo
       enddo

       do k  = KS,  KE
          Rmoist(k) = qdry(k)*Rdry + QTRC(k,i,j,I_QV)*Rvap
          pres  (k) = rhoe_new(k) * Rmoist(k) / CVmoist(k)
       enddo

       do k  = KS,  KE
          RHOT(k,i,j) = rhoe_new(k) / CVmoist(k) * ( PRE00 / pres(k) )**RovCP
       enddo

    enddo ! I loop
    enddo ! J loop

#ifdef _FPCOLL_
call STOP_COLLECTION('SUB_precipitation')
#endif
    call TIME_rapend  ('SUB_precipitation')

    return
  end subroutine ATMOS_PRECIPITATION

end module mod_atmos_precipitation
