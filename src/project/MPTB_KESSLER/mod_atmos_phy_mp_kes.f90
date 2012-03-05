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
  character(len=32), save :: WLABEL(11)

  logical, private, save  :: doreport_tendency = .false. ! report tendency of each process?

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
       doreport_tendency

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for KESSLER'

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

    WLABEL( 1) = "VAPOR"
    WLABEL( 2) = "CLOUD"
    WLABEL( 3) = "RAIN"
    WLABEL( 4) = "ICE"
    WLABEL( 5) = "SNOW"
    WLABEL( 6) = "GRAUPEL"
    WLABEL( 7) = "CLOUD_NUM"
    WLABEL( 8) = "RAIN_NUM"
    WLABEL( 9) = "ICE_NUM"
    WLABEL(10) = "SNOW_NUM"
    WLABEL(11) = "GRAUPEL_NUM"

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
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait,  &
       COMM_total
    use mod_grid, only : &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX, &
       KA   => GRID_KA,   &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       CDZ  => GRID_CDZ
    use mod_atmos_vars, only: &
       var => atmos_var, &
       A_NAME,           &
       VA  => A_VA,      &
       QA  => A_QA,      &
       I_DENS,           &
       I_RHOT,           &
       I_QV,             &
       I_QC,             &
       I_QR
    use mod_history, only: &
       HIST_in
    implicit none

    ! work
    real(8) :: dq_prcp(KA)        ! tendency q (precipitation)
    real(8) :: dq_cond1, dq_cond2 ! tendency q (condensation)
    real(8) :: dq_evap            ! tendency q (evaporation)
    real(8) :: dq_auto            ! tendency q (autoconversion)
    real(8) :: dq_coll            ! tendency q (collection)

    real(8) :: rho_prof(KA)       ! averaged profile of rho
    real(8) :: rho_fact(KA)
    real(8) :: pt_prev

    real(8) :: rho  ! density
    real(8) :: pt   ! potential temperature
    real(8) :: qv   ! water vapor
    real(8) :: qc   ! cloud water
    real(8) :: qr   ! cloud rain
    real(8) :: pres ! pressure [Pa]
    real(8) :: temp ! temperature

    real(8), parameter :: TVf1 = 36.34D0  ! Durran and Klemp (1983)
    real(8), parameter :: TVf2 = 0.1346D0 ! Durran and Klemp (1983)
    real(8) :: vel1, vel2, vent

    integer, parameter :: itelim = 1000
    real(8), parameter :: tt1 = 17.269D0 ! Tetens' formula
    real(8), parameter :: tt2 = 35.85D0  ! Tetens' formula
    real(8) :: qvs    ! saturated water vapor
    real(8) :: Rmoist, CPmoist, CVmoist
    real(8) :: LEovSE ! latent heat energy / sensible heat energy 
    real(8) :: efact

    ! for output
    real(8) :: output_prcp(KA,IA,JA)
    real(8) :: output_cond(KA,IA,JA)
    real(8) :: output_evap(KA,IA,JA)
    real(8) :: output_auto(KA,IA,JA)
    real(8) :: output_coll(KA,IA,JA)

    integer :: k, i, j, iv, ite
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics'

    ! averaged profile of rhoity [g/cc]
    do k = KS, KE
       rho_prof(k) = 0.D0

       do j = JS, JE
       do i = IS, IE
          rho_prof(k) = rho_prof(k) + var(k,i,j,I_DENS)
       enddo
       enddo

       rho_prof(k) = rho_prof(k) / real(IMAX*JMAX,kind=8) * 1.D-3

       rho_fact(k) = sqrt( rho_prof(KS)/rho_prof(k) )
    enddo

    ! prcpipitation flux
    do j = JS,   JE
    do i = IS,   IE
    do k = KS,   KE
       if ( qr > 0.D0 ) then
          if ( k >= KS .AND. k < KE ) then
             vel1 = TVf1 * ( rho_prof(k+1) * var(k+1,i,j,5+I_QR) )**TVf2 * rho_fact(k+1)
             vel2 = TVf1 * ( rho_prof(k  ) * var(k  ,i,j,5+I_QR) )**TVf2 * rho_fact(k  )

             dq_prcp(k) = ( var(k+1,i,j,I_DENS) * var(k+1,i,j,5+I_QR) * vel1 &
                          - var(k  ,i,j,I_DENS) * var(k  ,i,j,5+I_QR) * vel2 ) / CDZ(k)
          elseif( k == KE ) then
             vel2 = TVf1 * ( rho_prof(k  ) * var(k  ,i,j,5+I_QR) )**TVf2 * rho_fact(k  )

             dq_prcp(k) = ( 0.D0                                             &
                          - var(k  ,i,j,I_DENS) * var(k  ,i,j,5+I_QR) * vel2 ) / CDZ(k)
          endif
       else
          dq_prcp(k) = 0.D0
       endif
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rho = var(k,i,j,I_DENS)
       pt  = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS)
       qv  = var(k,i,j,5+I_QV)
       qc  = var(k,i,j,5+I_QC)
       qr  = var(k,i,j,5+I_QR)

       Rmoist  = Rdry  +  Rvap * qv
       CPmoist = CPdry + CPvap * qv
       CVmoist = CVdry + CVvap * qv

       pres = P00 * ( rho * pt * Rdry / P00 )**CPovCV
       temp = pres / ( rho * Rdry )
       qvs  = EPSvap * PSAT0 / pres * exp( tt1 * (temp-T00) / (temp-tt2) ) ! Tetens' formula

       LEovSE = ( LH0 * pt ) / ( CPdry * temp ) ! use prior value of temp & pt

       ! Saturation adjustment
       if ( qv > qvs ) then ! supersaturatd

          do ite = 1, itelim
             pt_prev = pt

             pres = P00 * ( rho * pt * Rdry / P00 )**CPovCV
             temp = pres / ( rho * Rdry )
             qvs  = EPSvap * PSAT0 / pres * exp( tt1 * (temp-T00) / (temp-tt2) )

             efact = LEovSE * CVdry / CVmoist                                     &
                   - pt * Rvap/CVmoist * ( 1.D0 - (Rdry/Rmoist) / (CPdry/CPmoist) )

             dq_cond1 = ( qv-qvs ) / ( 1.D0 + qvs * LH0 / CPdry * tt1 * ( T00-tt2 ) / (temp-tt2) )

             pt = pt + dq_cond1 * efact
             qv = qv - dq_cond1
             qc = qc + dq_cond1

             if ( qc < 0.D0 ) then ! re-adjust negative value
                pt = pt - qc * efact
                qv = qv + qc
                qc = 0.D0
             endif

             if( abs( pt-pt_prev ) <= 1.D-4 ) exit ! converged
          enddo

          if ( ite > itelim ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx not converged!', k,i,j,pt_prev,pt,qv,qc,dq_cond1
          endif

       else
          dq_cond1 = 0.D0
       endif

       ! Evaporation
       if ( qr > 0.D0 ) then

          pres = P00 * ( rho * pt * Rdry / P00 )**CPovCV
          temp = pres / ( rho * Rdry )
          qvs  = EPSvap * PSAT0 / pres * exp( tt1 * (temp-T00) / (temp-tt2) )

          if ( qvs > qv ) then ! unsaturatd
             vent = 1.6D0 + 30.3922D0 * ( rho * qr )**0.2046D0

             ! --- evaporation rate (Ogura and Takahashi, 1971) ---
             dq_evap = ( qvs - qv ) / qvs * vent * ( rho * qr )**0.525D0   &
                     / ( rho * ( 2.03D4 + 9.584D6 / ( pres * qvs ) ) )
             dq_evap = min( dq_evap, qr/dt )

             if( qv + dq_evap*dt > qvs ) dq_evap = ( qvs-qv ) / dt
          endif
       else
          dq_evap = 0.D0
       endif

       ! Autoconversion (ARPS)
       if ( qc > 1.D-3 ) then
          dq_auto = 1.D-3 * ( qc - 1.D-3 )
       else
          dq_auto = 0.D0
       endif

       ! Collection
       if ( qc > 0.D0 .and. qr > 0.D0 ) then
          dq_coll = 2.2D0 * qc * qr**0.875D0
       else
          dq_coll = 0.D0
       endif

!       if ( ( dq_auto + dq_coll )*dt > qc ) then
!          dq_auto = qc * dq_auto / ( dq_auto + dq_coll )
!          dq_coll = qc * dq_coll / ( dq_auto + dq_coll )
!       endif

       Rmoist  = Rdry  +  Rvap * qv
       CPmoist = CPdry + CPvap * qv
       CVmoist = CVdry + CVvap * qv

       efact = LEovSE * CVdry / CVmoist                                     &
             - pt * Rvap/CVmoist * ( 1.D0 - (Rdry/Rmoist) / (CPdry/CPmoist) )

       pt = pt + (                 -dq_evap )*dt * efact
       qv = qv + (                  dq_evap )*dt
       qc = qc + ( -dq_auto-dq_coll         )*dt
       qr = qr + (  dq_auto+dq_coll-dq_evap )*dt + dq_prcp(k)/rho*dt
       qv = max( qv, 0.D0 )

       if ( qc < 0.D0 ) then
          pt = pt - qc * efact
          qv = qv + qc
          qc = 0.D0
       endif

       if ( qr < 0.D0 ) then
          pt = pt - qr * efact
          qv = qv + qr
          qr = 0.D0
       endif

       pres = P00 * ( rho * pt * Rdry / P00 )**CPovCV
       temp = pres / ( rho * Rdry )
       qvs  = EPSvap * PSAT0 / pres * exp( tt1 * (temp-T00) / (temp-tt2) )

       ! Saturation adjustment (again)
       if ( qv > qvs ) then ! supersaturatd

          do ite = 1, itelim
             pt_prev = pt

             pres = P00 * ( rho * pt * Rdry / P00 )**CPovCV
             temp = pres / ( rho * Rdry )
             qvs  = EPSvap * PSAT0 / pres * exp( tt1 * (temp-T00) / (temp-tt2) )

             efact = LEovSE * CVdry / CVmoist                                     &
                   - pt * Rvap/CVmoist * ( 1.D0 - (Rdry/Rmoist) / (CPdry/CPmoist) )

             dq_cond2 = ( qv-qvs ) / ( 1.D0 + qvs * LH0 / CPdry * tt1 * ( T00-tt2 ) / (temp-tt2) )

             pt = pt + dq_cond2 * efact
             qv = qv - dq_cond2
             qc = qc + dq_cond2

             if ( qc < 0.D0 ) then ! re-adjust negative value
                pt = pt - qc * efact
                qv = qv + qc
                qc = 0.D0
             endif

             if( abs( pt-pt_prev ) <= 1.D-4 ) exit ! converged
          enddo

          if ( ite > itelim ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx not converged!', k,i,j,pt_prev,pt,qv,qc,dq_cond2
          endif

       else
          dq_cond2 = 0.D0
       endif

       var(k,i,j,I_DENS) = rho
       var(k,i,j,I_RHOT) = rho * pt
       var(k,i,j,5+I_QV) = qv
       var(k,i,j,5+I_QC) = qc
       var(k,i,j,5+I_QR) = qr

       if ( doreport_tendency ) then
          output_prcp(k,i,j) = dq_prcp(k)
          output_cond(k,i,j) = dq_cond1 + dq_cond2 
          output_evap(k,i,j) = dq_evap
          output_auto(k,i,j) = dq_auto
          output_coll(k,i,j) = dq_coll
       endif
    enddo
    enddo
    enddo

    if ( doreport_tendency ) then
       call HIST_in( output_prcp(:,:,:), 'dqprcp', 'tendency by precipitation' , 'kg/kg', '3D', dt )
       call HIST_in( output_cond(:,:,:), 'dqcond', 'tendency by condensation'  , 'kg/kg', '3D', dt )
       call HIST_in( output_evap(:,:,:), 'dqevap', 'tendency by evaporation'   , 'kg/kg', '3D', dt )
       call HIST_in( output_auto(:,:,:), 'dqauto', 'tendency by autoconversion', 'kg/kg', '3D', dt )
       call HIST_in( output_coll(:,:,:), 'dqcoll', 'tendency by collection'    , 'kg/kg', '3D', dt )
    endif

    ! fill IHALO & JHALO
    do iv = 1, VA
       call COMM_vars8( var(:,:,:,iv), iv )
       call COMM_wait ( var(:,:,:,iv), iv )
    enddo

    ! check total mass
    call COMM_total( var(:,:,:,:), A_NAME(:) )

    return
  end subroutine ATMOS_PHY_MP

end module mod_atmos_phy_mp 
