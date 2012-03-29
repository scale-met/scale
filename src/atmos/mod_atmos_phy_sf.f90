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

  real(8), private, save      :: U_minM  =    4.0D0  ! minimum U_abs for u,v,w
  real(8), private, save      :: U_minH  =    4.0D0  !                   T
  real(8), private, save      :: U_minE  =    4.0D0  !                   q
  real(8), private, parameter :: U_maxM  = 1000.0D0  ! maximum U_abs for u,v,w
  real(8), private, parameter :: U_maxH  = 1000.0D0  !                   T
  real(8), private, parameter :: U_maxE  = 1000.0D0  !                   q

  real(8), private, save      :: R10M                ! scaling factor for 10m value (momentum)
  real(8), private, parameter :: R10H =  1.D0        ! scaling factor for 10m value (heat)
  real(8), private, parameter :: R10E =  1.D0        ! scaling factor for 10m value (tracer)
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
       CDZ => GRID_CDZ
    implicit none

    real(8) :: ATMOS_PHY_SF_U_minM ! minimum U_abs for u,v,w
    real(8) :: ATMOS_PHY_SF_U_minH !                   T
    real(8) :: ATMOS_PHY_SF_U_minE !                   q
    real(8) :: ATMOS_PHY_SF_CM_min ! minimum bulk coef. of u,v,w
    real(8) :: ATMOS_PHY_SF_CH_min !                       T
    real(8) :: ATMOS_PHY_SF_CE_min !                       q

    NAMELIST / PARAM_ATMOS_PHY_SF / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_U_minH, &
       ATMOS_PHY_SF_U_minE, &
       ATMOS_PHY_SF_CM_min, &
       ATMOS_PHY_SF_CH_min, &
       ATMOS_PHY_SF_CE_min

    integer :: ierr
    !---------------------------------------------------------------------------

    ATMOS_PHY_SF_U_minM = U_minM
    ATMOS_PHY_SF_U_minH = U_minH
    ATMOS_PHY_SF_U_minE = U_minE
    ATMOS_PHY_SF_CM_min = CM_min
    ATMOS_PHY_SF_CH_min = CH_min
    ATMOS_PHY_SF_CE_min = CE_min

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

    R10M = 10.D0 / CDZ(KS) ! scale with height

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
       Rvap   => CONST_Rvap,   &
       RovCP  => CONST_RovCP, &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00,  &
       T00    => CONST_TEM00,  &
       LH0    => CONST_LH0,    &
       EPSvap => CONST_EPSvap, &
       PSAT0  => CONST_PSAT0
    use mod_grid, only : &
       CDZ => GRID_CDZ
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

    real(8) :: qdry, Rtot, pres, temp
    real(8) :: pres_evap ! partial pressure of water vapor at surface [Pa]
    real(8) :: qv_evap   ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j, iw
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    ! momentum -> velocity
    do j = JS-2, JE+2
    do i = IS-2, IE+2
       VELZ(i,j) = 2.D0 * MOMZ(KS,i,j) / ( DENS(KS+1,i,j)+DENS(KS,i,j) )
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
       !--- absolute velocity
       ! at (x, y, layer)
       Uabsw = ( ( VELZ(i,j)                 ) * 0.5D0 )**2 & ! surface is zero
             + ( ( VELX(i,j) + VELX(i-1,j  ) ) * 0.5D0 )**2 &
             + ( ( VELY(i,j) + VELY(i  ,j-1) ) * 0.5D0 )**2
       ! at (u, y, layer)
       Uabsu = ( ( VELZ(i,j  ) + VELZ(i+1,j  ) ) * 0.25D0 )**2 &
             + ( ( VELX(i,j  )                 )          )**2 &
             + ( ( VELY(i,j  ) + VELY(i+1,j  ) &
                 + VELY(i,j-1) + VELY(i+1,j-1) ) * 0.25D0 )**2
       ! at (x, v, layer)
       Uabsv = ( ( VELZ(i  ,j) + VELZ(i  ,j+1) ) * 0.25D0 )**2 &
             + ( ( VELX(i  ,j) + VELX(i  ,j+1) &
                 + VELX(i-1,j) + VELX(i-1,j+1) ) * 0.25D0 )**2 &
             + ( ( VELY(i  ,j)                 )          )**2

       !--- friction velocity
       Ustaru = max ( sqrt ( CM0 * Uabsu ), Ustar_min )
       Ustarv = max ( sqrt ( CM0 * Uabsv ), Ustar_min )
       Ustarw = max ( sqrt ( CM0 * Uabsw ), Ustar_min )

       !--- roughness lengths
       Z0Mu = max( Z0M0 + Z0MR/GRAV * Ustaru*Ustaru + Z0MS*visck / Ustaru, Z0M_min )
       Z0Mv = max( Z0M0 + Z0MR/GRAV * Ustarv*Ustarv + Z0MS*visck / Ustarv, Z0M_min )
       Z0Mw = max( Z0M0 + Z0MR/GRAV * Ustarw*Ustarw + Z0MS*visck / Ustarw, Z0M_min )
       Z0H  = max( Z0H0 + Z0HR/GRAV * Ustarw*Ustarw + Z0HS*visck / Ustarw, Z0H_min )
       Z0E  = max( Z0E0 + Z0ER/GRAV * Ustarw*Ustarw + Z0ES*visck / Ustarw, Z0E_min )

       !--- surface exchange coefficients
       CMX   = KARMAN*KARMAN / log( CDZ(KS)/Z0Mu )
       CMMu  = max( min( CMX / log( CDZ(KS)/Z0Mu ), CM_max ), CM_min ) * min( max( Uabsu, U_minM ), U_maxM )

       CMX   = KARMAN*KARMAN / log( CDZ(KS)/Z0Mv )
       CMMv  = max( min( CMX / log( CDZ(KS)/Z0Mv ), CM_max ), CM_min ) * min( max( Uabsv, U_minM ), U_maxM )

       CMX   = KARMAN*KARMAN / log( CDZ(KS)/Z0Mw )
       CMMw  = max( min( CMX / log( CDZ(KS)/Z0Mw ), CM_max ), CM_min ) * min( max( Uabsw, U_minM ), U_maxM )
       CMH   = max( min( CMX / log( CDZ(KS)/Z0H  ), CH_max ), CH_min ) * min( max( Uabsw, U_minH ), U_maxH )
       CME   = max( min( CMX / log( CDZ(KS)/Z0E  ), CE_max ), CE_min ) * min( max( Uabsw, U_minE ), U_maxE )

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
       qv_evap   = EPSvap * pres_evap / ( pres - pres_evap )

       !--- surface fluxes ( at x, y, 10m ) 
       SFLX_MOMZ(i,j) = DENS(KS,i,j) * Uabsw * R10M * CMMw * VELZ(i,j) * R10M
       SFLX_MOMX(i,j) = 0.5D0 * ( DENS(KS,i+1,j)+DENS(KS,i,j) ) * Uabsu * R10M * CMMu * VELX(i,j) * R10M
       SFLX_MOMY(i,j) = 0.5D0 * ( DENS(KS,i,j+1)+DENS(KS,i,j) ) * Uabsv * R10M * CMMv * VELY(i,j) * R10M
       SFLX_POTT(i,j) = DENS(KS,i,j) * Uabsw * R10M * CMH  * ( SST(1,i,j) - temp*R10H )
       SFLX_QV  (i,j) = DENS(KS,i,j) * Uabsw * R10M * CME  * ( qv_evap - QTRC(KS,i,j,1)*R10E )

    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sf
