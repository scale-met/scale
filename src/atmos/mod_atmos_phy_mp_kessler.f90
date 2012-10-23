!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Kessler parametarization
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-01-14 (Y.Miyamoto) [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP

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
  logical, private, save  :: MP_doreport_tendency = .false. ! report tendency of each process?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doreport_tendency

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** KESSLER-type parametarization'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_const, only: &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       CVdry  => CONST_CVdry,  &
       CPovCV => CONST_CPovCV, &
       Rvap   => CONST_Rvap,   &
       CPvap  => CONST_CPvap,  &
       CVvap  => CONST_CVvap,  &
       EPSvap => CONST_EPSvap, &
       P00    => CONST_PRE00,  &
       T00    => CONST_TEM00,  &
       PSAT0  => CONST_PSAT0,  &
       LH0    => CONST_LH0
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use mod_grid, only : &
       RCDZ => GRID_RCDZ
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total,    &
       DENS, &
       RHOT, &
       QTRC
    implicit none

    ! tendency
    real(RP) :: dq_prcp (KA) ! tendency q (precipitation)
    real(RP) :: dq_cond1(KA) ! tendency q (condensation)
    real(RP) :: dq_cond2(KA) ! tendency q (condensation)
    real(RP) :: dq_evap (KA) ! tendency q (evaporation)
    real(RP) :: dq_auto (KA) ! tendency q (autoconversion)
    real(RP) :: dq_coll (KA) ! tendency q (collection)
    ! for output
    real(RP) :: output_prcp(KA,IA,JA)
    real(RP) :: output_cond(KA,IA,JA)
    real(RP) :: output_evap(KA,IA,JA)
    real(RP) :: output_auto(KA,IA,JA)
    real(RP) :: output_coll(KA,IA,JA)

    ! disgnostics
    real(RP) :: pott(KA)    ! potential temperature [K]
    real(RP) :: pres(KA)    ! pressure [Pa]
    real(RP) :: temp(KA)    ! temperature [K]
    real(RP) :: qvs(KA)     ! saturated water vapor
    real(RP) :: Rmoist(KA)
    real(RP) :: CPmoist(KA)
    real(RP) :: CVmoist(KA)
    real(RP) :: LEovSE(KA)  ! latent heat energy / sensible heat energy 

    real(RP), parameter :: TVf1 = 36.34E0_RP  ! Durran and Klemp (1983)
    real(RP), parameter :: TVf2 = 0.1364E0_RP ! Durran and Klemp (1983)
    real(RP) :: rho_prof(KA) ! averaged profile of rho
    real(RP) :: rho_fact(KA)
    real(RP) :: vel1, vel2, vent

    integer, parameter :: itelim = 1000
    real(RP), parameter :: tt1 = 17.269E0_RP ! Tetens' formula
    real(RP), parameter :: tt2 = 35.85E0_RP  ! Tetens' formula
    real(RP) :: pt_prev
    real(RP) :: efact

    real(RP) :: diffq(KA)

    integer :: k, i, j, iq, ite
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics'

    do j = 1, JA
    do i = 1, IA
       ! total hydrometeor (before correction)
       do k = 1, KA
          diffq(k) = QTRC(k,i,j,I_QV) &
                   + QTRC(k,i,j,I_QC) &
                   + QTRC(k,i,j,I_QR)
       enddo

       ! remove negative value of hydrometeor (mass, number)
       do iq = I_QV, I_QR
       do k  = 1, KA
          if ( QTRC(k,i,j,iq) < 0.E0_RP ) then
             QTRC(k,i,j,iq) = 0.E0_RP
          endif
       enddo
       enddo

       ! apply correction of hydrometeor to total density
       do k  = 1, KA
          DENS(k,i,j) = DENS(k,i,j)        &
                      * ( 1.E0_RP             &
                        + QTRC(k,i,j,I_QV) &
                        + QTRC(k,i,j,I_QC) &
                        + QTRC(k,i,j,I_QR) &
                        - diffq(k)         ) ! after-before
       enddo
    enddo
    enddo

    ! averaged profile of density [g/cc]
    do k = KS, KE
       rho_prof(k) = 0.E0_RP
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rho_prof(k) = rho_prof(k) + DENS(k,i,j)
    enddo
    enddo
    enddo

    do k = KS, KE
       rho_prof(k) = rho_prof(k) / real(IMAX*JMAX,kind=RP) * 1.E-3_RP
       rho_fact(k) = sqrt( rho_prof(KS)/rho_prof(k) )
    enddo

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          pott(k) = RHOT(k,i,j) / DENS(k,i,j)

          Rmoist (k) = Rdry  * (1.E0_RP-QTRC(k,i,j,I_QV)) + Rvap  * QTRC(k,i,j,I_QV)
          CPmoist(k) = CPdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CPvap * QTRC(k,i,j,I_QV)
          CVmoist(k) = CVdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CVvap * QTRC(k,i,j,I_QV)

          pres(k) = P00 * ( DENS(k,i,j) * pott(k) * Rmoist(k) / P00 )**CPovCV
          temp(k) = pres(k) / ( DENS(k,i,j) * Rmoist(k) )
          qvs (k) = EPSvap * PSAT0 / pres(k) * exp( tt1 * (temp(k)-T00) / (temp(k)-tt2) ) ! Tetens' formula

          LEovSE(k) = ( LH0 * pott(k) ) / ( CPdry * temp(k) ) ! use prior value of temp & pt
       enddo

       ! Saturation adjustment
       do k = KS, KE
          if ( QTRC(k,i,j,I_QV) > qvs(k) ) then ! supersaturatd

             do ite = 1, itelim
                pt_prev = pott(k)

                Rmoist (k) = Rdry  * (1.E0_RP-QTRC(k,i,j,I_QV)) + Rvap  * QTRC(k,i,j,I_QV)
                CPmoist(k) = CPdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CPvap * QTRC(k,i,j,I_QV)
                CVmoist(k) = CVdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CVvap * QTRC(k,i,j,I_QV)

                pres(k) = P00 * ( DENS(k,i,j) * pott(k) * Rmoist(k) / P00 )**CPovCV
                temp(k) = pres(k) / ( DENS(k,i,j) * Rmoist(k) )
                qvs (k) = EPSvap * PSAT0 / pres(k) * exp( tt1 * (temp(k)-T00) / (temp(k)-tt2) )

                efact = LEovSE(k) * CVdry/CVmoist(k) &
                      - pott(k) * Rvap/CVmoist(k) * ( 1.E0_RP - (Rdry/Rmoist(k)) / (CPdry/CPmoist(k)) )

                dq_cond1(k) = ( QTRC(k,i,j,I_QV)-qvs(k) ) &
                            / ( 1.E0_RP + qvs(k) * LH0 / CPdry * tt1 * ( T00-tt2 ) / (temp(k)-tt2) )

                pott(k)          = pott(k)          + dq_cond1(k) * efact
                QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) - dq_cond1(k)
                QTRC(k,i,j,I_QC) = QTRC(k,i,j,I_QC) + dq_cond1(k)

                if ( QTRC(k,i,j,I_QC) < 0.E0_RP ) then ! re-adjust negative value
                   pott(k)          = pott(k)          - QTRC(k,i,j,I_QC) * efact
                   QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + QTRC(k,i,j,I_QC)
                   QTRC(k,i,j,I_QC) = 0.E0_RP
                endif

                if( abs( pott(k)-pt_prev ) <= 1.E-4_RP ) exit ! converged
             enddo

             if ( ite > itelim ) then
                if( IO_L ) write(IO_FID_LOG,*) 'xxx not converged!', &
                           k,i,j,pt_prev,pott(k),QTRC(k,i,j,I_QV),QTRC(k,i,j,I_QC),dq_cond1(k)
             endif

          else
             dq_cond1(k) = 0.E0_RP
          endif
       enddo

       ! Evaporation
       do k = KS, KE
          if ( QTRC(k,i,j,I_QR) > 0.E0_RP ) then
             Rmoist(k) = Rdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + Rvap * QTRC(k,i,j,I_QV)

             pres(k) = P00 * ( DENS(k,i,j) * pott(k) * Rmoist(k) / P00 )**CPovCV
             temp(k) = pres(k) / ( DENS(k,i,j) * Rmoist(k) )
             qvs (k) = EPSvap * PSAT0 / pres(k) * exp( tt1 * (temp(k)-T00) / (temp(k)-tt2) )

             if ( qvs(k) > QTRC(k,i,j,I_QV) ) then ! unsaturatd
                vent = 1.6E0_RP + 30.3922E0_RP * ( DENS(k,i,j) * QTRC(k,i,j,I_QR) )**0.2046E0_RP

                ! --- evaporation rate (Ogura and Takahashi, 1971) ---
                dq_evap(k) = ( qvs(k) - QTRC(k,i,j,I_QV) ) &
                           / qvs(k) * vent * ( DENS(k,i,j) * QTRC(k,i,j,I_QR) )**0.525E0_RP   &
                           / ( DENS(k,i,j) * ( 2.03E4_RP + 9.584E6_RP / ( pres(k) * qvs(k) ) ) )
                dq_evap(k) = min( dq_evap(k), QTRC(k,i,j,I_QR)/dt )

                if ( QTRC(k,i,j,I_QV) + dq_evap(k)*dt > qvs(k) ) then
                   dq_evap(k) = ( qvs(k)-QTRC(k,i,j,I_QV) ) / dt
                end if
             else
                dq_evap(k) = 0.E0_RP
             endif
          else
             dq_evap(k) = 0.E0_RP
          endif
       enddo

       do k = KS, KE
          ! Autoconversion (ARPS)
          if ( QTRC(k,i,j,I_QC) > 1.E-3_RP ) then
             dq_auto(k) = 1.E-3_RP * ( QTRC(k,i,j,I_QC) - 1.E-3_RP )
          else
             dq_auto(k) = 0.E0_RP
          endif

          ! Collection
          if ( QTRC(k,i,j,I_QC) > 0.E0_RP .and. QTRC(k,i,j,I_QR) > 0.E0_RP ) then
             dq_coll(k) = 2.2E0_RP * QTRC(k,i,j,I_QC) * QTRC(k,i,j,I_QR)**0.875E0_RP
          else
             dq_coll(k) = 0.E0_RP
          endif

!          if ( ( dq_auto + dq_coll )*dt > QTRC(k,i,j,I_QC) ) then
!             dq_auto = QTRC(k,i,j,I_QC) * dq_auto / ( dq_auto + dq_coll )
!             dq_coll = QTRC(k,i,j,I_QC) * dq_coll / ( dq_auto + dq_coll )
!          endif
       enddo

       ! Precipitation flux
       do k = KS, KE
          vel1 = TVf1 * ( rho_prof(k+1) * QTRC(k+1,i,j,I_QR) )**TVf2 * rho_fact(k+1)
          vel2 = TVf1 * ( rho_prof(k  ) * QTRC(k  ,i,j,I_QR) )**TVf2 * rho_fact(k  )
          if( k == KE ) vel1 = 0.E0_RP

          dq_prcp(k) = ( DENS(k+1,i,j) * QTRC(k+1,i,j,I_QR) * vel1 &
                       - DENS(k  ,i,j) * QTRC(k  ,i,j,I_QR) * vel2 ) * RCDZ(k)
       enddo

       ! Update POTT, QV, QC, QR
       do k = KS, KE
          Rmoist (k) = Rdry  * (1.E0_RP-QTRC(k,i,j,I_QV)) + Rvap  * QTRC(k,i,j,I_QV)
          CPmoist(k) = CPdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CPvap * QTRC(k,i,j,I_QV)
          CVmoist(k) = CVdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CVvap * QTRC(k,i,j,I_QV)

          efact = LEovSE(k) * CVdry/CVmoist(k) &
                - pott(k) * Rvap/CVmoist(k) * ( 1.E0_RP - (Rdry/Rmoist(k)) / (CPdry/CPmoist(k)) )

          pott(k)          = pott(k)          + (                       -dq_evap(k) )*dt * efact
          QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + (                        dq_evap(k) )*dt
          QTRC(k,i,j,I_QC) = QTRC(k,i,j,I_QC) + ( -dq_auto(k)-dq_coll(k)            )*dt
          QTRC(k,i,j,I_QR) = QTRC(k,i,j,I_QR) + (  dq_auto(k)+dq_coll(k)-dq_evap(k) )*dt + dq_prcp(k)/DENS(k,i,j)*dt
       enddo

       ! negtive fixer
       do k = KS, KE
          QTRC(k,i,j,I_QV) = max( QTRC(k,i,j,I_QV), 0.E0_RP )
 
          if ( QTRC(k,i,j,I_QC) < 0.E0_RP ) then
             pott(k)          = pott(k)          - QTRC(k,i,j,I_QC) * efact
             QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + QTRC(k,i,j,I_QC)
             QTRC(k,i,j,I_QC) = 0.E0_RP
          endif

          if ( QTRC(k,i,j,I_QR) < 0.E0_RP ) then
             pott(k)          = pott(k)          - QTRC(k,i,j,I_QR) * efact
             QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + QTRC(k,i,j,I_QR)
             QTRC(k,i,j,I_QR) = 0.E0_RP
          endif
       enddo

       ! Saturation adjustment (again)
       do k = KS, KE
                Rmoist (k) = Rdry  * (1.E0_RP-QTRC(k,i,j,I_QV)) + Rvap  * QTRC(k,i,j,I_QV)
                CPmoist(k) = CPdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CPvap * QTRC(k,i,j,I_QV)
                CVmoist(k) = CVdry * (1.E0_RP-QTRC(k,i,j,I_QV)) + CVvap * QTRC(k,i,j,I_QV)

          pres(k) = P00 * ( DENS(k,i,j) * pott(k) * Rmoist(k) / P00 )**CPovCV
          temp(k) = pres(k) / ( DENS(k,i,j) * Rmoist(k) )
          qvs (k) = EPSvap * PSAT0 / pres(k) * exp( tt1 * (temp(k)-T00) / (temp(k)-tt2) ) ! Tetens' formula

          if ( QTRC(k,i,j,I_QV) > qvs(k) ) then ! supersaturatd

             do ite = 1, itelim
                pt_prev = pott(k)

                pres(k) = P00 * ( DENS(k,i,j) * pott(k) * Rmoist(k) / P00 )**CPovCV
                temp(k) = pres(k) / ( DENS(k,i,j) * Rmoist(k) )
                qvs (k) = EPSvap * PSAT0 / pres(k) * exp( tt1 * (temp(k)-T00) / (temp(k)-tt2) )

                efact = LEovSE(k) * CVdry/CVmoist(k)                                               &
                      - pott(k) * Rvap/CVmoist(k) * ( 1.E0_RP - (Rdry/Rmoist(k)) / (CPdry/CPmoist(k)) )

                dq_cond2(k) = ( QTRC(k,i,j,I_QV)-qvs(k) ) &
                            / ( 1.E0_RP + qvs(k) * LH0 / CPdry * tt1 * ( T00-tt2 ) / (temp(k)-tt2) )

                pott(k)          = pott(k)          + dq_cond2(k) * efact
                QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) - dq_cond2(k)
                QTRC(k,i,j,I_QC) = QTRC(k,i,j,I_QC) + dq_cond2(k)

                if ( QTRC(k,i,j,I_QC) < 0.E0_RP ) then ! re-adjust negative value
                   pott(k)          = pott(k)          - QTRC(k,i,j,I_QC) * efact
                   QTRC(k,i,j,I_QV) = QTRC(k,i,j,I_QV) + QTRC(k,i,j,I_QC)
                   QTRC(k,i,j,I_QC) = 0.E0_RP
                endif

                if( abs( pott(k)-pt_prev ) <= 1.E-4_RP ) exit ! converged
             enddo

             if ( ite > itelim ) then
                if( IO_L ) write(IO_FID_LOG,*) 'xxx not converged!(after)', &
                           k,i,j,pt_prev,pott(k),QTRC(k,i,j,I_QV),QTRC(k,i,j,I_QC),dq_cond2(k)
             endif

          else
             dq_cond2(k) = 0.E0_RP
          endif
       enddo

       do k = KS, KE
          RHOT(k,i,j) = pott(k) * DENS(k,i,j)
       enddo

       if ( MP_doreport_tendency ) then
          do k = KS, KE
             output_prcp(k,i,j) = dq_prcp (k)
             output_cond(k,i,j) = dq_cond1(k) + dq_cond2(k)
             output_evap(k,i,j) = dq_evap (k)
             output_auto(k,i,j) = dq_auto (k)
             output_coll(k,i,j) = dq_coll (k)
          enddo
       endif
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       ! total hydrometeor (before correction)
       do k = 1, KA
          diffq(k) = QTRC(k,i,j,I_QV) &
                   + QTRC(k,i,j,I_QC) &
                   + QTRC(k,i,j,I_QR)
       enddo

       ! remove negative value of hydrometeor (mass, number)
       do iq = I_QV, I_QR
       do k  = 1, KA
          if ( QTRC(k,i,j,iq) < 0.D0 ) then
             QTRC(k,i,j,iq) = 0.D0
          endif
       enddo
       enddo

       ! apply correction of hydrometeor to total density
       do k  = 1, KA
          DENS(k,i,j) = DENS(k,i,j)        &
                      * ( 1.D0             &
                        + QTRC(k,i,j,I_QV) &
                        + QTRC(k,i,j,I_QC) &
                        + QTRC(k,i,j,I_QR) &
                        - diffq(k)         ) ! after-before
       enddo
    enddo
    enddo

    call ATMOS_vars_fillhalo

    call ATMOS_vars_total

    if ( MP_doreport_tendency ) then
       call HIST_in( output_prcp(:,:,:), 'dqprcp', 'tendency by precipitation' , 'kg/kg', '3D', dt )
       call HIST_in( output_cond(:,:,:), 'dqcond', 'tendency by condensation'  , 'kg/kg', '3D', dt )
       call HIST_in( output_evap(:,:,:), 'dqevap', 'tendency by evaporation'   , 'kg/kg', '3D', dt )
       call HIST_in( output_auto(:,:,:), 'dqauto', 'tendency by autoconversion', 'kg/kg', '3D', dt )
       call HIST_in( output_coll(:,:,:), 'dqcoll', 'tendency by collection'    , 'kg/kg', '3D', dt )
    endif

    return
  end subroutine ATMOS_PHY_MP

end module mod_atmos_phy_mp 
