!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-03 (Y.Miyamoto) [new]
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE3
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2012-04-10 (Y.Miyamoto) [mod] introduce coefficients for interpolation
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf
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
  public :: ATMOS_PHY_SF_setup
  public :: ATMOS_PHY_SF

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
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
  real(8), private, parameter :: CM0   = 1.0D-3  ! bulk coef. for U*
  real(8), private, parameter :: visck = 1.5D-5  ! kinematic viscosity 

  ! parameters
  real(8), private, parameter :: Z0M0 =   0.0D0      ! base
  real(8), private, parameter :: Z0MR = 0.018D0      ! rough factor
  real(8), private, parameter :: Z0MS =  0.11D0      ! smooth factor
  real(8), private, parameter :: Z0H0 =   1.4D-5
  real(8), private, parameter :: Z0HR =   0.0D0
  real(8), private, parameter :: Z0HS =   0.4D0
  real(8), private, parameter :: Z0E0 =   1.3D-4
  real(8), private, parameter :: Z0ER =   0.0D0
  real(8), private, parameter :: Z0ES =  0.62D0

  ! limiter
  real(8), private, parameter :: Ustar_min =  1.0D-3 ! minimum limit of U*

  real(8), private, parameter :: Z0M_min =    1.0D-5 ! minimum roughness length of u,v,w
  real(8), private, parameter :: Z0H_min =    1.0D-5 !                             T
  real(8), private, parameter :: Z0E_min =    1.0D-5 !                             q

  real(8), private, save      :: CM_min  =    1.0D-5 ! minimum bulk coef. of u,v,w
  real(8), private, save      :: CH_min  =    1.0D-5 !                       T
  real(8), private, save      :: CE_min  =    1.0D-5 !                       q
  real(8), private, parameter :: CM_max  =    2.5D-3 ! maximum bulk coef. of u,v,w
  real(8), private, parameter :: CH_max  =    1.0D0  !                       T
  real(8), private, parameter :: CE_max  =    1.0D0  !                       q

  real(8), private, save      :: U_minM  =    0.0D0  ! minimum U_abs for u,v,w
  real(8), private, save      :: U_minH  =    0.0D0  !                   T
  real(8), private, save      :: U_minE  =    0.0D0  !                   q
  real(8), private, parameter :: U_maxM  =  100.0D0  ! maximum U_abs for u,v,w
  real(8), private, parameter :: U_maxH  =  100.0D0  !                   T
  real(8), private, parameter :: U_maxE  =  100.0D0  !                   q
  logical, private, save      :: flgcnst = .true.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only : &
       CDZ => GRID_CDZ, &
       CZ  => GRID_CZ,  &
       FZ  => GRID_FZ
    implicit none

    real(8) :: ATMOS_PHY_SF_U_minM ! minimum U_abs for u,v,w
    real(8) :: ATMOS_PHY_SF_U_minH !                   T
    real(8) :: ATMOS_PHY_SF_U_minE !                   q
    real(8) :: ATMOS_PHY_SF_CM_min ! minimum bulk coef. of u,v,w
    real(8) :: ATMOS_PHY_SF_CH_min !                       T
    real(8) :: ATMOS_PHY_SF_CE_min !                       q
    logical :: ATMOS_PHY_SF_flgcnst 

    NAMELIST / PARAM_ATMOS_PHY_SF / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_U_minH, &
       ATMOS_PHY_SF_U_minE, &
       ATMOS_PHY_SF_CM_min, &
       ATMOS_PHY_SF_CH_min, &
       ATMOS_PHY_SF_CE_min, &
       ATMOS_PHY_SF_flgcnst

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    ATMOS_PHY_SF_U_minM  = U_minM
    ATMOS_PHY_SF_U_minH  = U_minH
    ATMOS_PHY_SF_U_minE  = U_minE
    ATMOS_PHY_SF_CM_min  = CM_min
    ATMOS_PHY_SF_CH_min  = CH_min
    ATMOS_PHY_SF_CE_min  = CE_min
    ATMOS_PHY_SF_flgcnst = flgcnst

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACE]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF)

    U_minM = ATMOS_PHY_SF_U_minM
    U_minH = ATMOS_PHY_SF_U_minH
    U_minE = ATMOS_PHY_SF_U_minE
    CM_min = ATMOS_PHY_SF_CM_min
    CH_min = ATMOS_PHY_SF_CH_min
    CE_min = ATMOS_PHY_SF_CE_min
    flgcnst = ATMOS_PHY_SF_flgcnst 

    return
  end subroutine ATMOS_PHY_SF_setup

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN, &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       Rvap   => CONST_Rvap,   &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00,  &
       T00    => CONST_TEM00,  &
       LH0    => CONST_LH0,    &
       EPSvap => CONST_EPSvap, &
       PSAT0  => CONST_PSAT0
    use mod_time, only: &
       dttb => TIME_DTSEC_ATMOS_PHY_TB
    use mod_grid, only : &
       CDZ => GRID_CDZ, &
       CZ  => GRID_CZ
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_ocean_vars, only: &
       SST
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    implicit none

    ! monitor
    real(8) :: SHFLX(1,IA,JA) ! sensible heat flux [W/m2]
    real(8) :: LHFLX(1,IA,JA) ! latent   heat flux [W/m2]

    real(8) :: FB  = 9.4d0  ! Louis factor b (bM)
    real(8) :: FBS = 4.7d0  ! Louis factor b' (bM/eM = dE/eE = 9.4/2.0)
    real(8) :: FDM = 7.4d0  ! Louis factor d of u (dM)
    real(8) :: FDH = 5.3d0  ! Louis factor d of T, q (dH)
    real(8) :: THS = 300.0d0

    ! work
    real(8) :: VELZ(IA,JA)
    real(8) :: VELX(IA,JA)
    real(8) :: VELY(IA,JA)

    real(8) :: Uabsu  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(8) :: Uabsv
    real(8) :: Uabsw
    real(8) :: Ustaru ! friction velocity [m/s]
    real(8) :: Ustarv
    real(8) :: Ustarw

    real(8) :: Z0Mu   ! roughness length [m] (momentum,heat,tracer)
    real(8) :: Z0Mv
    real(8) :: Z0Mw
    real(8) :: Z0H
    real(8) :: Z0E

    real(8) :: CMX
    real(8) :: CMMu   ! surface exchange coefficient (momentum,heat,tracer)
    real(8) :: CMMv
    real(8) :: CMMw
    real(8) :: CMH
    real(8) :: CME

    real(8) :: FM, FH
    real(8) :: FCM, FCH
    real(8) :: FZH, FZE
    real(8) :: RIX, RIB
    real(8) :: qdry, Rtot, pres, temp, qvap
    real(8) :: pres_evap ! partial pressure of water vapor at surface [Pa]
    real(8) :: qv_evap   ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j, iw
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    ! momentum -> velocity
    do j = JS-2, JE+2
    do i = IS-2, IE+2
       VELZ(i,j) = MOMZ(KS+1,i,j) / ( DENS(KS+1,i,j)+DENS(KS  ,i,j) ) &
                 + MOMZ(KS  ,i,j) / ( DENS(KS  ,i,j)+DENS(KS  ,i,j) ) 
    enddo
    enddo

    do j = JS-2, JE+2
    do i = IS-2, IE+1
       VELX(i,j) = 2.D0 * MOMX(KS,i,j) / ( DENS(KS,i+1,j)+DENS(KS,i,j) )
    enddo
    enddo

    do j = JS-2, JE+1
    do i = IS-2, IE+2
       VELY(i,j) = 2.D0 * MOMY(KS,i,j) / ( DENS(KS,i,j+1)+DENS(KS,i,j) )
    enddo
    enddo

    do j = JS-1, JE
    do i = IS-1, IE

       ! Gas constant
       qdry = 1.D0
       do iw = QQS, QQE
          qdry = qdry - QTRC(KS,i,j,iw)
       enddo
       Rtot = Rdry*qdry + Rvap*QTRC(KS,i,j,I_QV)

       !--- Qv at sea surface
       pres      = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV
       temp      = ( RHOT(KS,i,j) / DENS(KS,i,j) ) * ( P00 / pres )**RovCP
       pres_evap = PSAT0 * exp( LH0/Rvap * ( 1.D0/T00 - 1.D0/SST(1,i,j) ) )
!       qv_evap   = EPSvap * pres_evap / ( pres - pres_evap )
       qv_evap   = EPSvap * pres_evap / P00

       !--- absolute velocity
       ! at (x, y, layer)
       Uabsw = sqrt( ( ( VELZ(i,j)                 )          )**2 &
                   + ( ( VELX(i,j) + VELX(i-1,j  ) ) * 0.50D0 )**2 &
                   + ( ( VELY(i,j) + VELY(i  ,j-1) ) * 0.50D0 )**2 )
       ! at (u, y, layer)
       Uabsu = sqrt( ( ( VELZ(i,j  ) + VELZ(i+1,j  ) ) * 0.50D0 )**2 &
                   + ( ( VELX(i,j  )                 )          )**2 &
                   + ( ( VELY(i,j  ) + VELY(i+1,j  ) &
                       + VELY(i,j-1) + VELY(i+1,j-1) ) * 0.25D0 )**2 )
       ! at (x, v, layer)
       Uabsv = sqrt( ( ( VELZ(i  ,j) + VELZ(i  ,j+1) ) * 0.50D0 )**2 &
                   + ( ( VELX(i  ,j) + VELX(i  ,j+1) &
                       + VELX(i-1,j) + VELX(i-1,j+1) ) * 0.25D0 )**2 &
                   + ( ( VELY(i  ,j)                 )          )**2 )

       !--- friction velocity at u, v, and w points
!       Ustaru = max ( sqrt ( CM0 ) * Uabsu , Ustar_min )
!       Ustarv = max ( sqrt ( CM0 ) * Uabsv , Ustar_min )
!       Ustarw = max ( sqrt ( CM0 ) * Uabsw , Ustar_min )
       !--- friction velocity at u, v, and w points are fixed as 0.25D0 for DYCOMS-II RF02
       Ustaru = 0.25D0
       Ustarv = 0.25D0
       Ustarw = 0.25D0

       !--- roughness lengths at u, v, and w points
       Z0Mu = max( Z0M0 + Z0MR/GRAV * Ustaru*Ustaru + Z0MS*visck / Ustaru, Z0M_min )
       Z0Mv = max( Z0M0 + Z0MR/GRAV * Ustarv*Ustarv + Z0MS*visck / Ustarv, Z0M_min )
       Z0Mw = max( Z0M0 + Z0MR/GRAV * Ustarw*Ustarw + Z0MS*visck / Ustarw, Z0M_min )
       Z0H  = max( Z0H0 + Z0HR/GRAV * Ustarw*Ustarw + Z0HS*visck / Ustarw, Z0H_min )
       Z0E  = max( Z0E0 + Z0ER/GRAV * Ustarw*Ustarw + Z0ES*visck / Ustarw, Z0E_min )
       FZH  = log( Z0Mw/Z0H ) / log( CZ(KS)/Z0Mw )
       FZE  = log( Z0Mw/Z0E ) / log( CZ(KS)/Z0Mw )

       CMX  = ( KARMAN / log( CZ(KS)/Z0Mw ) )**2
       RIB  = GRAV/THS * CZ(KS) * ( temp*( P00/pres )**RovCP - SST(1,i,j) ) / Uabsw
       FCM  = FDM * CMX * FB * sqrt( CZ(KS)/Z0Mw )
       FCH  = FDH * CMX * FB * sqrt( CZ(KS)/Z0Mw )
       if ( RIB >= 0.D0 ) then
          FM   = 1.D0/( 1.D0 + FBS * RIB*1.D0/( 1.D0 + FZH ) )**2
          FH   = FM
       else
          RIX  = sqrt( abs( RIB*1.D0/( 1.D0 + FZH ) ) )
          FM   = 1.D0 - FB*RIB*1.D0/( 1.D0 + FZH )/( 1.D0 + FCM*RIX )
          FH   = 1.D0 - FB*RIB*1.D0/( 1.D0 + FZH )/( 1.D0 + FCH*RIX )
       endif

       !--- surface exchange coefficients
       CMX   = ( KARMAN / log( CZ(KS)/Z0Mu ) )**2 * FM
       CMMu  = max( min( CMX , CM_max ), CM_min ) * min( max( Uabsu, U_minM ), U_maxM )

       CMX   = ( KARMAN / log( CZ(KS)/Z0Mv ) )**2 * FM
       CMMv  = max( min( CMX , CM_max ), CM_min ) * min( max( Uabsv, U_minM ), U_maxM )

       CMX   = ( KARMAN / log( CZ(KS)/Z0Mw ) )**2 * FM
       CMMw  = max( min( CMX , CM_max ), CM_min ) * min( max( Uabsw, U_minM ), U_maxM )

       CMX   = ( KARMAN*KARMAN / log( CZ(KS)/Z0Mw ) / log( CZ(KS)/Z0H ) ) &
               * FH * 1.D0/( 1.D0 + FZH/sqrt( FM )*FH )
       CMH   = max( min( CMX , CH_max ), CH_min ) * min( max( Uabsw, U_minH ), U_maxH )

       CMX   = ( KARMAN*KARMAN / log( CZ(KS)/Z0Mw ) / log( CZ(KS)/Z0E ) ) &
               * FH * 1.D0/( 1.D0 + FZE/sqrt( FM )*FH )
       CME   = max( min( CMX , CE_max ), CE_min ) * min( max( Uabsw, U_minE ), U_maxE )

       !--- surface fluxes ( at x, y, 10m ) 
       SFLX_MOMX(i,j) = - 0.5D0 * ( DENS(KS,i+1,j)+DENS(KS,i,j) ) * CMMu * VELX(i,j)
       SFLX_MOMY(i,j) = - 0.5D0 * ( DENS(KS,i,j+1)+DENS(KS,i,j) ) * CMMv * VELY(i,j)
       SFLX_MOMZ(i,j) = - DENS(KS,i,j) * CMMw * VELZ(i,j) 
       if( flgcnst ) then
        SFLX_POTT(i,j) =   16.d0/CPdry
        SFLX_QV  (i,j) =   93.d0/LH0
       else
        SFLX_POTT(i,j) =   DENS(KS,i,j) * CMH * ( SST(1,i,j) - RHOT(KS,i,j) / DENS(KS,i,j) )
        SFLX_QV  (i,j) =   DENS(KS,i,j) * CME * ( qv_evap - QTRC(KS,i,j,I_QV) )
       endif

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       SHFLX(1,i,j) = SFLX_POTT(i,j) * CPdry
       LHFLX(1,i,j) = SFLX_QV  (i,j) * LH0
    enddo
    enddo

    call HIST_in( SHFLX(:,:,:), 'SHFLX', 'sensible heat flux', 'W/m2', '2D', dttb )
    call HIST_in( LHFLX(:,:,:), 'LHFLX', 'latent heat flux',   'W/m2', '2D', dttb )

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sf
