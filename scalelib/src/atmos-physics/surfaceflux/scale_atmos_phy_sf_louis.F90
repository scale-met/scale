!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-03 (Y.Miyamoto)  [new]
!! @li      2011-12-11 (H.Yashiro)   [mod] integrate to SCALE-LES ver.3
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-04-10 (Y.Miyamoto)  [mod] introduce coefficients for interpolation
!! @li      2012-09-11 (S.Nishizawa) [mod] bugfix based on the scale document
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_sf_louis
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_louis_setup
  public :: ATMOS_PHY_SF_louis

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
  real(RP), private, parameter :: Cm0   = 1.0E-3_RP  ! bulk coef. for U*
  real(RP), private, parameter :: visck = 1.5E-5_RP  ! kinematic viscosity

  ! parameters
  real(RP), private, save :: Z00 = 0.0_RP      ! base
  real(RP), private, save :: Z0R = 0.018_RP    ! rough factor
  real(RP), private, save :: Z0S = 0.11_RP     ! smooth factor
  real(RP), private, save :: Zt0 = 1.4E-5_RP
  real(RP), private, save :: ZtR = 0.0_RP
  real(RP), private, save :: ZtS = 0.4_RP
  real(RP), private, save :: Ze0 = 1.3E-4_RP
  real(RP), private, save :: ZeR = 0.0_RP
  real(RP), private, save :: ZeS = 0.62_RP
  real(RP), private, save :: ThS = 300.0_RP

  ! limiter
  real(RP), private, parameter :: Ustar_min =  1.0E-3_RP ! minimum limit of U*

  real(RP), private, parameter :: Z0_min =    1.0E-5_RP ! minimum roughness length of u,v,w
  real(RP), private, parameter :: Zt_min =    1.0E-5_RP !                             T
  real(RP), private, parameter :: Ze_min =    1.0E-5_RP !                             q

  real(RP), private, save      :: Cm_min  =    1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, save      :: Ch_min  =    1.0E-5_RP !                       T
  real(RP), private, save      :: Ce_min  =    1.0E-5_RP !                       q
  real(RP), private, parameter :: Cm_max  =    2.5E-3_RP ! maximum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_max  =    1.0_RP  !                       T
  real(RP), private, parameter :: Ce_max  =    1.0_RP  !                       q

  real(RP), private, save      :: U_minM  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, save      :: U_minH  =    0.0_RP  !                   T
  real(RP), private, save      :: U_minE  =    0.0_RP  !                   q
  real(RP), private, parameter :: U_maxM  =  100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH  =  100.0_RP  !                   T
  real(RP), private, parameter :: U_maxE  =  100.0_RP  !                   q

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_louis_setup( ATMOS_TYPE_PHY_SF )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_PHY_SF

    real(RP) :: ATMOS_PHY_SF_U_minM ! minimum U_abs for u,v,w
    real(RP) :: ATMOS_PHY_SF_U_minH !                   T
    real(RP) :: ATMOS_PHY_SF_U_minE !                   q
    real(RP) :: ATMOS_PHY_SF_CM_min ! minimum bulk coef. of u,v,w
    real(RP) :: ATMOS_PHY_SF_CH_min !                       T
    real(RP) :: ATMOS_PHY_SF_CE_min !                       q
    real(RP) :: ATMOS_PHY_SF_Z00
    real(RP) :: ATMOS_PHY_SF_Z0R
    real(RP) :: ATMOS_PHY_SF_Z0S
    real(RP) :: ATMOS_PHY_SF_Zt0
    real(RP) :: ATMOS_PHY_SF_ZtR
    real(RP) :: ATMOS_PHY_SF_ZtS
    real(RP) :: ATMOS_PHY_SF_Ze0
    real(RP) :: ATMOS_PHY_SF_ZeR
    real(RP) :: ATMOS_PHY_SF_ZeS
    real(RP) :: ATMOS_PHY_SF_ThS

    NAMELIST / PARAM_ATMOS_PHY_SF_LOUIS / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_U_minH, &
       ATMOS_PHY_SF_U_minE, &
       ATMOS_PHY_SF_CM_min, &
       ATMOS_PHY_SF_CH_min, &
       ATMOS_PHY_SF_CE_min, &
       ATMOS_PHY_SF_Z00, &
       ATMOS_PHY_SF_Z0R, &
       ATMOS_PHY_SF_Z0S, &
       ATMOS_PHY_SF_Zt0, &
       ATMOS_PHY_SF_ZtR, &
       ATMOS_PHY_SF_ZtS, &
       ATMOS_PHY_SF_Ze0, &
       ATMOS_PHY_SF_ZeR, &
       ATMOS_PHY_SF_ZeS, &
       ATMOS_PHY_SF_ThS

    integer :: ierr
    !---------------------------------------------------------------------------

    ATMOS_PHY_SF_U_minM = U_minM
    ATMOS_PHY_SF_U_minH = U_minH
    ATMOS_PHY_SF_U_minE = U_minE
    ATMOS_PHY_SF_CM_min = CM_min
    ATMOS_PHY_SF_CH_min = CH_min
    ATMOS_PHY_SF_CE_min = CE_min
    ATMOS_PHY_SF_Z00    = Z00
    ATMOS_PHY_SF_Z0R    = Z0R
    ATMOS_PHY_SF_Z0S    = Z0S
    ATMOS_PHY_SF_Zt0    = Zt0
    ATMOS_PHY_SF_ZtR    = ZtR
    ATMOS_PHY_SF_ZtS    = ZtS
    ATMOS_PHY_SF_Ze0    = Ze0
    ATMOS_PHY_SF_ZeR    = ZeR
    ATMOS_PHY_SF_ZeS    = ZeS
    ATMOS_PHY_SF_ThS    = ThS

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACEFLUX]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Louis parameterization'

    if ( ATMOS_TYPE_PHY_SF /= 'LOUIS' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_SF is not LOUIS. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_LOUIS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_LOUIS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_LOUIS)

    U_minM = ATMOS_PHY_SF_U_minM
    U_minH = ATMOS_PHY_SF_U_minH
    U_minE = ATMOS_PHY_SF_U_minE
    CM_min = ATMOS_PHY_SF_CM_min
    CH_min = ATMOS_PHY_SF_CH_min
    CE_min = ATMOS_PHY_SF_CE_min
    Z00    = ATMOS_PHY_SF_Z00
    Z0R    = ATMOS_PHY_SF_Z0R
    Z0S    = ATMOS_PHY_SF_Z0S
    Zt0    = ATMOS_PHY_SF_Zt0
    ZtR    = ATMOS_PHY_SF_ZtR
    ZtS    = ATMOS_PHY_SF_ZtS
    Ze0    = ATMOS_PHY_SF_Ze0
    ZeR    = ATMOS_PHY_SF_ZeR
    ZeS    = ATMOS_PHY_SF_ZeS
    ThS    = ATMOS_PHY_SF_ThS

    return
  end subroutine ATMOS_PHY_SF_louis_setup

  !-----------------------------------------------------------------------------
  ! calclate surface flux
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_louis( &
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, SST,             & ! (in)
         CZ, ctime                                            ) ! (in)
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN, &
       Rdry   => CONST_Rdry,   &
       Rvap   => CONST_Rvap,   &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00,  &
       T00    => CONST_TEM00,  &
       LH0    => CONST_LH0,    &
       EPSvap => CONST_EPSvap, &
       PSAT0  => CONST_PSAT0
    implicit none

    real(RP), intent(out) :: SFLX_MOMZ(IA,JA)
    real(RP), intent(out) :: SFLX_MOMX(IA,JA)
    real(RP), intent(out) :: SFLX_MOMY(IA,JA)
    real(RP), intent(out) :: SFLX_POTT(IA,JA)
    real(RP), intent(out) :: SFLX_QV  (IA,JA)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: SST (1,IA,JA)

    real(RP), intent(in)  :: CZ(KA)
    real(DP), intent(in)  :: ctime

    real(RP) :: FB  = 9.4_RP  ! Louis factor b (bM)
    real(RP) :: FBS = 4.7_RP  ! Louis factor b' (bM/eM = dE/eE = 9.4/2.0)
    real(RP) :: FDM = 7.4_RP  ! Louis factor d of u (dM)
    real(RP) :: FDH = 5.3_RP  ! Louis factor d of T, q (dH)
    real(RP) :: FR  = 0.74_RP ! turbulent Prandtl number (Businger et al. 1971)

    ! work
    real(RP) :: THETA(IA,JA)

    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Ustar ! friction velocity [m/s]

    real(RP) :: Z0   ! roughness length [m] (momentum,heat,tracer)
    real(RP) :: Zt
    real(RP) :: Ze

    real(RP) :: Cm   ! bulk coefficient (momentum,heat,tracer)
    real(RP) :: Ch
    real(RP) :: Ce

    real(RP) :: a2
    real(RP) :: Fm, Fh, Psih
    real(RP) :: RiB
    real(RP) :: qdry, Rtot, pres, temp
    real(RP) :: pres_evap ! partial pressure of water vapor at surface [Pa]
    real(RP) :: qv_evap   ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j, iw
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    ! rho*theta -> potential temperature at cell centor
    do j = JS, JE
    do i = IS, IE
       THETA(i,j) = RHOT(KS,i,j) / DENS(KS,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE

       ! at cell center

       !--- absolute velocity
       Uabs = sqrt( &
              ( MOMZ(KS,i,j)                  )**2 &
            + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
            + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 &
            ) / DENS(KS,i,j) * 0.5_RP

       !--- friction velocity at u, v, and w points
       Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

       !--- roughness lengths at u, v, and w points
       Z0 = max( Z00 + Z0R/GRAV * Ustar*Ustar + Z0S*visck / Ustar, Z0_min )
       Zt = max( Zt0 + ZtR/GRAV * Ustar*Ustar + ZtS*visck / Ustar, Zt_min )
       Ze = max( Ze0 + ZeR/GRAV * Ustar*Ustar + ZeS*visck / Ustar, Ze_min )

       call get_RiB( &
            RiB, Fm, Fh, Psih,             & ! (out)
            THETA(i,j), SST(1,i,j), Uabs**2, & ! (in)
            CZ(KS), Z0, Zt,                & ! (in)
            KARMAN, FB, FBS, FDM, FDH,     & ! (in)
            ThS, GRAV                    ) ! (in)

       !--- surface exchange coefficients
       a2 = ( KARMAN / log( CZ(KS)/Z0 ) )**2
       Cm = a2 * Fm
       Ch = a2 * Fh / ( FR * ( log( Z0/Zt ) / Psih + 1.0_RP ) )
       Ce = a2 * Fh / ( FR * ( log( Z0/Ze ) / Psih + 1.0_RP ) )

       ! Gas constant
       qdry = 1.0_RP
       do iw = QQS, QQE
          qdry = qdry - QTRC(KS,i,j,iw)
       enddo
       Rtot = Rdry*qdry + Rvap*QTRC(KS,i,j,I_QV)

       !--- saturation at surface
       pres      = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV
       temp      = ( RHOT(KS,i,j) / DENS(KS,i,j) ) * ( P00 / pres )**RovCP
       pres_evap = PSAT0 * exp( LH0/Rvap * ( 1.0_RP/T00 - 1.0_RP/SST(1,i,j) ) )
!       qv_evap   = EPSvap * pres_evap / ( pres - pres_evap )
       qv_evap   = EPSvap * pres_evap / P00

       ! flux
       SFLX_MOMZ(i,j) = - min(max(Cm,Cm_min),Cm_max) * min(max(Uabs,U_minM),U_maxM) &
            * MOMZ(KS,i,j) * 0.5_RP
       SFLX_POTT(i,j) =   min(max(Ch,Ch_min),Ch_max) * min(max(Uabs,U_minH),U_maxH) &
            * ( SST(1,i,j)*DENS(KS,i,j) - RHOT(KS,i,j) )
       SFLX_QV  (i,j) =   min(max(Ce,Ce_min),Ce_max) * min(max(Uabs,U_minE),U_maxE) &
            * DENS(KS,i,j) * ( qv_evap - QTRC(KS,i,j,I_QV) )



       ! at (u, y, layer)
       Uabs = sqrt( &
              ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i+1,j)                                     ) )**2 &
            + ( 2.0_RP *   MOMX(KS,i,j)                                                        )**2 &
            + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
       Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

       Z0 = max( Z00 + Z0R/GRAV * Ustar*Ustar + Z0S*visck / Ustar, Z0_min )
       Zt = max( Zt0 + ZtR/GRAV * Ustar*Ustar + ZtS*visck / Ustar, Zt_min )

       call get_RiB( &
            RiB, Fm, Fh, Psih,                    & ! (out)
            ( THETA(i,j)+THETA(i+1,j) ) * 0.5_RP, & ! (in)
            ( SST(1,i,j)+SST(1,i+1,j) ) * 0.5_RP,     & ! (in)
            Uabs**2,                              & ! (in)
            CZ(KS), Z0, Zt,                       & ! (in)
            KARMAN, FB, FBS, FDM, FDH,            & ! (in)
            ThS, GRAV                             ) ! (in)

       a2 = ( KARMAN / log( CZ(KS)/Z0 ) )**2
       Cm = a2 * Fm

       SFLX_MOMX(i,j) = - min(max(Cm,Cm_min),Cm_min) * min(max(Uabs,U_minM),U_maxM) &
            * MOMX(KS,i,j)



       ! at (x, v, layer)
       Uabs = sqrt( &
              ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i,j+1)                                     ) )**2 &
            + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(KS,i,j)                                                        )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
       Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

       Z0 = max( Z00 + Z0R/GRAV * Ustar*Ustar + Z0S*visck / Ustar, Z0_min )
       Zt = max( Zt0 + ZtR/GRAV * Ustar*Ustar + ZtS*visck / Ustar, Zt_min )

       call get_RiB( &
            RiB, Fm, Fh, Psih,                    & ! (out)
            ( THETA(i,j)+THETA(i,j+1) ) * 0.5_RP, & ! (in)
            ( SST(1,i,j)+SST(1,i,j+1) ) * 0.5_RP,     & ! (in)
            Uabs**2,                              & ! (in)
            CZ(KS), Z0, Zt,                       & ! (in)
            KARMAN, FB, FBS, FDM, FDH,            & ! (in)
            ThS, GRAV                             ) ! (in)

       a2 = ( KARMAN / log( CZ(KS)/Z0 ) )**2
       Cm = a2 * Fm

       SFLX_MOMY(i,j) = - min(max(Cm,Cm_min),Cm_min) * min(max(Uabs,U_minM),U_maxM) &
            * MOMY(KS,i,j)

    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_louis


  subroutine get_RiB( &
       RiB, Fm, Fh, Psih,    &
       theta, theta_sfc, u2, &
       Z, Z0, Zt,            &
       K, FB, FBS, FDM, FDH, &
       ThS, G                )
    use scale_const, only: &
       EPS    => CONST_EPS
    real(RP), intent(out) :: RiB
    real(RP), intent(out) :: Fm
    real(RP), intent(out) :: Fh
    real(RP), intent(out) :: Psih
    real(RP), intent(in)  :: theta
    real(RP), intent(in)  :: theta_sfc
    real(RP), intent(in)  :: u2
    real(RP), intent(in)  :: Z
    real(RP), intent(in)  :: Z0
    real(RP), intent(in)  :: Zt
    real(RP), intent(in)  :: K
    real(RP), intent(in)  :: FB
    real(RP), intent(in)  :: FBS
    real(RP), intent(in)  :: FDM
    real(RP), intent(in)  :: FDH
    real(RP), intent(in)  :: ThS
    real(RP), intent(in)  :: G

    real(RP) :: tmp

    ! the first guess of RiB0 (= RiBt)
    RiB = G/ThS * z * (  theta -  theta_sfc ) / ( u2+EPS )

    ! Fm, Fh, Psi_h/R
    if ( RiB >= 0 ) then
       Fm = 1.0_RP / ( 1.0_RP + FBS * Rib )**2
       Fh = Fm
    else
       tmp = ( K / log( Z/Z0 ) )**2 * FB * sqrt( Z/Z0 * abs(RiB) )
       Fm = 1.0_RP - FB * RiB / ( 1.0_RP + FDM * tmp )
       Fh = 1.0_RP - FB * RiB / ( 1.0_RP + FDH * tmp )
    endif
    Psih = log( Z/Z0 ) * sqrt( Fm ) / Fh

    ! the final estimate of RiB0
    tmp = log( Z0/Zt )
    RiB = RiB - RiB * tmp / ( tmp + Psih )

    ! Fm, Fh, Psih/R
    if ( RiB >= 0.0_RP ) then
       Fm = 1.0_RP / ( 1.0_RP + FBS * Rib )**2
       Fh = Fm
    else
       tmp = ( K / log( Z/Z0 ) )**2 * FB * sqrt( Z/Z0 * abs(RiB) )
       Fm = 1.0_RP - FB * RiB / ( 1.0_RP + FDM * tmp )
       Fh = 1.0_RP - FB * RiB / ( 1.0_RP + FDH * tmp )
    endif
    Psih = log( Z/Z0 ) * sqrt( Fm ) / Fh

    return
  end subroutine get_RiB


end module scale_atmos_phy_sf_louis
