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
!! @li      2011-xx-xx (Y.Miyamoto) [new]
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE3
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
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
  real(8), private, save :: CM0   = 1.0D-3  ! bulk coef. for U*
  real(8), private, save :: visck = 1.5D-5  ! kinematic viscosity 

  ! parameters
  real(8), private, save :: Z0M0 =   0.0D0      ! base
  real(8), private, save :: Z0MR = 0.018D0      ! rough factor
  real(8), private, save :: Z0MS =  0.11D0      ! smooth factor
  real(8), private, save :: Z0H0 =   1.4D-5
  real(8), private, save :: Z0HR =   0.0D0
  real(8), private, save :: Z0HS =   0.4D0
  real(8), private, save :: Z0E0 =   1.3D-4
  real(8), private, save :: Z0ER =   0.0D0
  real(8), private, save :: Z0ES =  0.62D0

  ! limiter
  real(8), private, save :: Ustar_min =  1.0D-3 ! minimum limit of U*

  real(8), private, save :: Z0M_min =    1.0D-5 ! minimum roughness length of u,v,w
  real(8), private, save :: Z0H_min =    1.0D-5 !                             T
  real(8), private, save :: Z0E_min =    1.0D-5 !                             q

  real(8), private, save :: CM_min  =    1.0D-5 ! minimum bulk coef. of u,v,w
  real(8), private, save :: CH_min  =    1.0D-5 !                       T
  real(8), private, save :: CE_min  =    1.0D-5 !                       q
  real(8), private, save :: CM_max  =    2.5D-3 ! maximum bulk coef. of u,v,w
  real(8), private, save :: CH_max  =    1.0D0  !                       T
  real(8), private, save :: CE_max  =    1.0D0  !                       q

  real(8), private, save :: U_minM  =    4.0D0  ! minimum U_abs for u,v,w
  real(8), private, save :: U_minH  =    4.0D0  !                   T
  real(8), private, save :: U_minE  =    4.0D0  !                   q
  real(8), private, save :: U_maxM  = 1000.0D0  ! maximum U_abs for u,v,w
  real(8), private, save :: U_maxH  = 1000.0D0  !                   T
  real(8), private, save :: U_maxE  = 1000.0D0  !                   q
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
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

    return
  end subroutine ATMOS_PHY_SF_setup

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF( dens,   pott,   qtrc,          &
                           pres,   velx,   vely,   velz,  &
                           FLXij_sfc, FLXt_sfc, FLXqv_sfc )
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN, &
       Rvap   => CONST_Rvap,   &
       T00    => CONST_TEM00,  &
       LH0    => CONST_LH0,    &
       EPSvap => CONST_EPSvap, &
       PSAT0  => CONST_PSAT0
    use mod_grid, only : &
       IA  => GRID_IA,  &
       JA  => GRID_JA,  &
       KA  => GRID_KA,  &
       IS  => GRID_IS,  &
       IE  => GRID_IE,  &
       JS  => GRID_JS,  &
       JE  => GRID_JE,  &
       KS  => GRID_KS,  &
       WS  => GRID_WS,  &
       DZ  => GRID_DZ
    use mod_atmos_vars, only: &
       QA => A_QA, &
       I_QV
    use mod_ocean_vars, only: &
       OCEAN_vars_get
    implicit none

    ! prognostic value
    real(8), intent(in)  :: dens(IA,JA,KA)      ! density [kg/m3]
    real(8), intent(in)  :: pott(IA,JA,KA)      ! potential temperature [K]
    real(8), intent(in)  :: qtrc(IA,JA,KA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    real(8), intent(in)  :: pres(IA,JA,KA)      ! pressure [Pa]
    real(8), intent(in)  :: velx(IA,JA,KA)      ! velocity(x) [m/s]
    real(8), intent(in)  :: vely(IA,JA,KA)      ! velocity(y) [m/s]
    real(8), intent(in)  :: velz(IA,JA,KA)      ! velocity(z) [m/s]

    ! surface flux
    real(8), intent(out) :: FLXij_sfc(IA,JA,3)  ! => FLXij(1:IA,1:JA,WS,1:3,3)
    real(8), intent(out) :: FLXt_sfc (IA,JA)    ! => FLXt (1:IA,1:JA,WS)
    real(8), intent(out) :: FLXqv_sfc(IA,JA)    ! => FLXq (1:IA,1:JA,WS,1)

    ! from ocean
    real(8) :: sst(IA,JA) ! sea surface temperature

    real(8) :: R10M   ! scaling factor for 10m value (momentum,heat,tracer)
    real(8) :: R10H
    real(8) :: R10E

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

    real(8) :: U10u   ! absolute velocity at the lowermost atmospheric layer [m/s]
    real(8) :: U10v
    real(8) :: U10w

    real(8) :: pres_vap ! partial pressure of water vapor at surface [Pa]
    real(8) :: qv_sfc   ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    call OCEAN_vars_get( sst )

    R10M = 10.D0 / DZ ! scale with height
    R10H =  1.D0      ! assume homogeneous
    R10E =  1.D0      ! assume homogeneous

    do j = JS-1, JE
    do i = IS-1, IE
       !--- absolute velocity ( at x, y, interface )
       Uabsu = ( velx(i,j,KS)                                                             )**2.D0 &
             + ( ( vely(i,j,KS)+vely(i+1,j,KS)+vely(i,j-1,KS)+vely(i+1,j-1,KS) ) * 0.25D0 )**2.D0 &
             + ( ( velz(i,j,WS+1)+velz(i+1,j,WS+1)                             ) * 0.25D0 )**2.D0

       Uabsv = ( ( velx(i,j,KS)+velx(i,j+1,KS)+velx(i-1,j,KS)+velx(i-1,j+1,KS) ) * 0.25D0 )**2.D0 &
             + ( vely(i,j,KS)                                                             )**2.D0 &
             + ( ( velz(i,j,WS+1)+velz(i+1,j,WS+1)                             ) * 0.25D0 )**2.D0

       Uabsw = ( ( velx(i,j,KS)+velx(i-1,j,KS) )  * 0.5D0 )**2.D0 &
             + ( ( vely(i,j,KS)+vely(i,j-1,KS) )  * 0.5D0 )**2.D0 &
             + ( velz(i,j,WS+1)                   * 0.5D0 )**2.D0

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
       CMX   = KARMAN*KARMAN / log( DZ/Z0Mu )
       CMMu  = max( min( CMX / log( DZ/Z0Mu ), CM_max ), CM_min ) * min( max( Uabsu, U_minM ), U_maxM )

       CMX   = KARMAN*KARMAN / log( DZ/Z0Mv )
       CMMv  = max( min( CMX / log( DZ/Z0Mv ), CM_max ), CM_min ) * min( max( Uabsv, U_minM ), U_maxM )

       CMX   = KARMAN*KARMAN / log( DZ/Z0Mw )
       CMMw  = max( min( CMX / log( DZ/Z0Mw ), CM_max ), CM_min ) * min( max( Uabsw, U_minM ), U_maxM )

       CMH   = max( min( CMX / log( DZ/Z0H  ), CH_max ), CH_min ) * min( max( Uabsw, U_minH ), U_maxH )
       CME   = max( min( CMX / log( DZ/Z0E  ), CE_max ), CE_min ) * min( max( Uabsw, U_minE ), U_maxE )

       !--- Qv at sea surface
       pres_vap = PSAT0 * exp( LH0/Rvap * ( 1.D0/T00 - 1.D0/sst(i,j) ) )
       qv_sfc   = EPSvap * pres_vap / ( pres(i,j,KS) - pres_vap )

       !--- velocity at the 10m height
       U10u = Uabsu * R10M
       U10v = Uabsv * R10M
       U10w = Uabsw * R10M

       !--- surface fluxes ( at x, y, 10m ) 
       FLXij_sfc(i,j,1) = 0.5D0 * ( dens(i+1,j,KS)+dens(i,j,KS) ) * U10u * CMMu * velx(i,j,KS)   * R10M
       FLXij_sfc(i,j,2) = 0.5D0 * ( dens(i,j+1,KS)+dens(i,j,KS) ) * U10v * CMMv * vely(i,j,KS)   * R10M
       FLXij_sfc(i,j,3) = dens(i,j,KS) * U10w * CMMw * velz(i,j,WS+1) * R10M

       FLXt_sfc (i,j)   = dens(i,j,KS) * U10w * CMH  * ( sst(i,j) - pott(i,j,KS)*R10H )

       FLXqv_sfc(i,j)   = dens(i,j,KS) * U10w * CME  * ( qv_sfc - qtrc(i,j,KS,I_QV)*R10E )
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sf
