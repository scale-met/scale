!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Precipitation transport
!!
!! @par Description
!!          Precipitation Transport module
!!          (1D advection)
!!
!! @author Team SCALE
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
    use mod_const, only: &
       GRAV  => CONST_GRAV,  &
       Rdry  => CONST_Rdry,  &
       CVdry => CONST_CVdry, &
       Rvap  => CONST_Rvap,  &
       PRE00 => CONST_PRE00
    use mod_grid, only: &
       CZ   => GRID_CZ,   &
       FDZ  => GRID_FDZ,   &
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
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_history, only: &
       HIST_in
    implicit none

    real(RP), intent(out)   :: flux_rain(KA,IA,JA)
    real(RP), intent(out)   :: flux_snow(KA,IA,JA)
    real(RP), intent(inout) :: velw     (KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(RP), intent(in)    :: rhoq     (KA,IA,JA,QA) ! rho * q
    real(RP), intent(in)    :: rhoe     (KA,IA,JA)
    real(RP), intent(in)    :: temp     (KA,IA,JA)
    real(RP), intent(in)    :: dt

    real(RP) :: qflx    (KA,QA)
    real(RP) :: eflx    (KA)
    real(RP) :: rhoe_new(KA)
    real(RP) :: qdry, CVtot, Rtot, pres, RovCP
    real(RP) :: eflx_out(KA,IA,JA)

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call TIME_rapstart('SUB_precipitation')

    do iq = 1, QA
       call COMM_vars8( velw(:,:,:,iq), iq )
    end do
    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), QA+iq )
    end do
    do iq = 1, QA
       call COMM_wait( velw(:,:,:,iq), iq )
    end do
    do iq = 1, QA
       call COMM_wait( QTRC(:,:,:,iq), QA+iq )
    end do

!OCL XFILL
    flux_rain(:,:,:) = 0.0_RP
    flux_snow(:,:,:) = 0.0_RP

    ! tracer/energy transport by falldown
    ! 1st order upwind, forward euler, velocity is always negative
    do j  = JS, JE
    do i  = IS, IE

       !--- tracer
!OCL XFILL
       do iq = I_QC, QA
          do k  = KS-1, KE-1
             qflx(k,iq) = velw(k+1,i,j,iq) * rhoq(k+1,i,j,iq)
          enddo
          qflx(KE,iq) = 0.0_RP
       enddo

       !--- lowermost flux is saved for land process
       do k  = KS-1, KE
           do iq = QWS, QWE
              flux_rain(k,i,j) = flux_rain(k,i,j) - qflx(k,iq)
           enddo
           if( QIS > 0 ) then
             do iq = QIS, QIE
               flux_snow(k,i,j) = flux_snow(k,i,j) - qflx(k,iq)
             enddo
           endif
       enddo 

!OCL XFILL
       do k = KS, KE
          rhoe_new(k) = rhoe(k,i,j)
       enddo
       eflx(KE) = 0.0_RP

       do iq = I_QC, QQE

          !--- internal energy
          do k  = KS-1, KE-1
             eflx(k) = qflx(k,iq) * temp(k+1,i,j) * CVw(iq)
             eflx_out(k,i,j) = qflx(k,iq) * temp(k+1,i,j) * CVw(iq)
          enddo
          do k  = KS,   KE
             rhoe_new(k) = rhoe_new(k) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- potential energy
          do k  = KS, KE-1
             eflx(k) = qflx(k,iq) * GRAV * FDZ(k)
          enddo
          eflx(KS-1) = qflx(KS-1,iq) * GRAV * CZ(KS)
          do k  = KS,   KE
             rhoe_new(k) = rhoe_new(k) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- momentum z (half level)
          do k  = KS-1, KE-1
             eflx(k) = 0.25_RP * ( velw(k+1,i,j,iq) + velw(k,i,j,iq) ) &
                               * ( QTRC(k+1,i,j,iq) + QTRC(k,i,j,iq) ) &
                               * MOMZ(k,i,j)
          enddo
          do k  = KS,  KE-1
             MOMZ(k,i,j) = MOMZ(k,i,j) - dt * ( eflx(k+1) - eflx(k) ) * RFDZ(k)
          enddo

          !--- momentum x
          do k  = KS-1, KE-1
             eflx(k) = 0.25_RP * ( velw(k+1,i,j,iq) + velw(k+1,i+1,j,iq) ) &
                               * ( QTRC(k+1,i,j,iq) + QTRC(k+1,i+1,j,iq) ) &
                               * MOMX(k+1,i,j)
          enddo
          do k  = KS,  KE
             MOMX(k,i,j) = MOMX(k,i,j) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- momentum y
          do k  = KS-1, KE-1
             eflx(k) = 0.25_RP * ( velw(k+1,i,j,iq) + velw(k+1,i,j+1,iq) ) &
                               * ( QTRC(k+1,i,j,iq) + QTRC(k+1,i,j+1,iq) ) &
                               * MOMY(k+1,i,j)
          enddo
          do k  = KS,  KE
             MOMY(k,i,j) = MOMY(k,i,j) - dt * ( eflx(k) - eflx(k-1) ) * RCDZ(k)
          enddo

          !--- update total density
          do k  = KS,  KE
             DENS(k,i,j) = DENS(k,i,j) - dt * ( qflx(k,iq) - qflx(k-1,iq) ) * RCDZ(k)
          enddo

       enddo ! QC-QG loop

       !--- update tracer
       do iq = I_QC, QA
       do k  = KS, KE
          QTRC(k,i,j,iq) = ( rhoq(k,i,j,iq) - dt * ( qflx(k,iq) - qflx(k-1,iq) ) * RCDZ(k) ) / DENS(k,i,j)
       enddo 
       enddo ! QC-QG,NC-NG loop

       do k  = KS,  KE
          qdry  = 1.0_RP
          CVtot = 0.0_RP
          do iq = QQS, QQE
             qdry  = qdry  - QTRC(k,i,j,iq)
             CVtot = CVtot + QTRC(k,i,j,iq) * CVw(iq)
          enddo
          CVtot = CVtot + qdry * CVdry

          Rtot  = qdry * Rdry + QTRC(k,i,j,I_QV) * Rvap
          RovCP = Rtot / ( CVtot + Rtot )

          pres = rhoe_new(k) * Rtot / CVtot
          RHOT(k,i,j) = rhoe_new(k) / CVtot * ( PRE00 / pres )**RovCP
       enddo

    enddo ! I loop
    enddo ! J loop

    call TIME_rapend  ('SUB_precipitation')

    call HIST_in( flux_rain(:,:,:), 'RFLX', 'precipitation flux', 'kg/m2/s', dt)
    call HIST_in( flux_snow(:,:,:), 'SFLX', 'precipitation flux', 'kg/m2/s', dt)
    call HIST_in( eflx_out(:,:,:), 'EFLX', 'energy flux', 'J/kg/s', dt)

    return
  end subroutine ATMOS_PRECIPITATION

end module mod_atmos_precipitation
