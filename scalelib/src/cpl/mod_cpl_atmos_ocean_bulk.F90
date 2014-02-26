!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Ocean Surface fluxes
!!
!! @par Description
!!          Surface flux between atmosphere and ocean with Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-26 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_ocean_bulk
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_AtmOcn_bulk_setup
  public :: CPL_AtmOcn_bulk

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: bulkflux

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer,  private, save :: nmax = 100 ! maximum iteration number

  real(RP), private, save :: res_min  =   1.0_RP    ! minimum number of residual
  real(RP), private, save :: dTS      =   1.0E-8_RP ! delta TS
  real(RP), private, save :: U_minM   =   0.0_RP    ! minimum U_abs for u,v,w
  real(RP), private, save :: U_minH   =   0.0_RP    !                   T
  real(RP), private, save :: U_minE   =   0.0_RP    !                   q
  real(RP), private, save :: U_maxM   = 100.0_RP    ! maximum U_abs for u,v,w
  real(RP), private, save :: U_maxH   = 100.0_RP    !                   T
  real(RP), private, save :: U_maxE   = 100.0_RP    !                   q

  real(RP), private, save :: EMIT  = 1.0_RP    ! emissivity in long-wave radiation [no unit]
  real(RP), private, save :: TCW   = 0.6_RP    ! thermal conductivity for water [W/m/K]
  real(RP), private, save :: CM0   = 1.0E-3_RP ! bulk coef. for U*
  real(RP), private, save :: visck = 1.5E-5_RP ! kinematic viscosity

  real(RP), private, save :: Ustar_min = 1.0E-3_RP ! minimum fiction velocity
  real(RP), private, save :: Z0M_min   = 1.0E-5_RP ! minimum roughness length for u,v,w
  real(RP), private, save :: Z0H_min   = 1.0E-5_RP !                              T
  real(RP), private, save :: Z0E_min   = 1.0E-5_RP !                              q

  real(RP), private, save :: Z0MI = 0.0E-0_RP ! base roughness rength for u,v,w
  real(RP), private, save :: Z0MR = 1.8E-2_RP ! rough factor for u,v,w
  real(RP), private, save :: Z0MS = 1.1E-1_RP ! smooth factor for u,v,w
  real(RP), private, save :: Z0HI = 1.4E-5_RP ! base roughness rength for T
  real(RP), private, save :: Z0HR = 0.0E-0_RP ! rough factor for T
  real(RP), private, save :: Z0HS = 4.0E-1_RP ! smooth factor for T
  real(RP), private, save :: Z0EI = 1.3E-4_RP ! base roughness rength for q
  real(RP), private, save :: Z0ER = 0.0E-0_RP ! rough factor for q
  real(RP), private, save :: Z0ES = 6.2E-1_RP ! smooth factor for q

contains
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_bulk_setup( CPL_TYPE_AtmOcn )
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_bulkcoef, only: &
       CPL_bulkcoef_setup
    implicit none

    character(len=H_SHORT), intent(in) :: CPL_TYPE_AtmOcn

    integer  :: CPL_AtmOcn_bulk_nmax
    real(RP) :: CPL_AtmOcn_bulk_res_min
    real(RP) :: CPL_AtmOcn_bulk_dTS
    real(RP) :: CPL_AtmOcn_bulk_U_minM
    real(RP) :: CPL_AtmOcn_bulk_U_minH
    real(RP) :: CPL_AtmOcn_bulk_U_minE
    real(RP) :: CPL_AtmOcn_bulk_U_maxM
    real(RP) :: CPL_AtmOcn_bulk_U_maxH
    real(RP) :: CPL_AtmOcn_bulk_U_maxE
    real(RP) :: CPL_AtmOcn_bulk_EMIT
    real(RP) :: CPL_AtmOcn_bulk_TCW
    real(RP) :: CPL_AtmOcn_bulk_CM0
    real(RP) :: CPL_AtmOcn_bulk_visck
    real(RP) :: CPL_AtmOcn_bulk_Z0M_min
    real(RP) :: CPL_AtmOcn_bulk_Z0H_min
    real(RP) :: CPL_AtmOcn_bulk_Z0E_min
    real(RP) :: CPL_AtmOcn_bulk_Z0MI
    real(RP) :: CPL_AtmOcn_bulk_Z0MR
    real(RP) :: CPL_AtmOcn_bulk_Z0MS
    real(RP) :: CPL_AtmOcn_bulk_Z0HI
    real(RP) :: CPL_AtmOcn_bulk_Z0HR
    real(RP) :: CPL_AtmOcn_bulk_Z0HS
    real(RP) :: CPL_AtmOcn_bulk_Z0EI
    real(RP) :: CPL_AtmOcn_bulk_Z0ER
    real(RP) :: CPL_AtmOcn_bulk_Z0ES

    NAMELIST / PARAM_CPL_ATMOCN_BULK / &
       CPL_AtmOcn_bulk_nmax,    &
       CPL_AtmOcn_bulk_res_min, &
       CPL_AtmOcn_bulk_dTS,     &
       CPL_AtmOcn_bulk_U_minM,  &
       CPL_AtmOcn_bulk_U_minH,  &
       CPL_AtmOcn_bulk_U_minE,  &
       CPL_AtmOcn_bulk_U_maxM,  &
       CPL_AtmOcn_bulk_U_maxH,  &
       CPL_AtmOcn_bulk_U_maxE,  &
       CPL_AtmOcn_bulk_EMIT,    &
       CPL_AtmOcn_bulk_TCW,     &
       CPL_AtmOcn_bulk_CM0,     &
       CPL_AtmOcn_bulk_visck,   &
       CPL_AtmOcn_bulk_Z0M_min, &
       CPL_AtmOcn_bulk_Z0H_min, &
       CPL_AtmOcn_bulk_Z0E_min, &
       CPL_AtmOcn_bulk_Z0MI,    &
       CPL_AtmOcn_bulk_Z0MR,    &
       CPL_AtmOcn_bulk_Z0MS,    &
       CPL_AtmOcn_bulk_Z0HI,    &
       CPL_AtmOcn_bulk_Z0HR,    &
       CPL_AtmOcn_bulk_Z0HS,    &
       CPL_AtmOcn_bulk_Z0EI,    &
       CPL_AtmOcn_bulk_Z0ER,    &
       CPL_AtmOcn_bulk_Z0ES

    integer :: ierr
    !---------------------------------------------------------------------------

    CPL_AtmOcn_bulk_nmax    = nmax
    CPL_AtmOcn_bulk_res_min = res_min
    CPL_AtmOcn_bulk_dTS     = dTS
    CPL_AtmOcn_bulk_U_minM  = U_minM
    CPL_AtmOcn_bulk_U_minH  = U_minH
    CPL_AtmOcn_bulk_U_minE  = U_minE
    CPL_AtmOcn_bulk_U_maxM  = U_maxM
    CPL_AtmOcn_bulk_U_maxH  = U_maxH
    CPL_AtmOcn_bulk_U_maxE  = U_maxE
    CPL_AtmOcn_bulk_EMIT    = EMIT
    CPL_AtmOcn_bulk_TCW     = TCW
    CPL_AtmOcn_bulk_CM0     = CM0
    CPL_AtmOcn_bulk_visck   = visck
    CPL_AtmOcn_bulk_Z0M_min = Z0M_min
    CPL_AtmOcn_bulk_Z0H_min = Z0H_min
    CPL_AtmOcn_bulk_Z0E_min = Z0E_min
    CPL_AtmOcn_bulk_Z0MI    = Z0MI
    CPL_AtmOcn_bulk_Z0MR    = Z0MR
    CPL_AtmOcn_bulk_Z0MS    = Z0MS
    CPL_AtmOcn_bulk_Z0HI    = Z0HI
    CPL_AtmOcn_bulk_Z0HR    = Z0HR
    CPL_AtmOcn_bulk_Z0HS    = Z0HS
    CPL_AtmOcn_bulk_Z0EI    = Z0EI
    CPL_AtmOcn_bulk_Z0ER    = Z0ER
    CPL_AtmOcn_bulk_Z0ES    = Z0ES

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean: bulk flux parameter'

    if ( CPL_TYPE_AtmOcn /= 'U95'  .and. &
         CPL_TYPE_AtmOcn /= 'BH91'       ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx CPL_TYPE_AtmOcn is not U95 or BH91. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_ATMOCN_BULK,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_ATMOCN_BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL_ATMOCN_BULK)

    nmax    = CPL_AtmOcn_bulk_nmax
    res_min = CPL_AtmOcn_bulk_res_min
    dTS     = CPL_AtmOcn_bulk_dTS
    U_minM  = CPL_AtmOcn_bulk_U_minM
    U_minH  = CPL_AtmOcn_bulk_U_minH
    U_minE  = CPL_AtmOcn_bulk_U_minE
    U_maxM  = CPL_AtmOcn_bulk_U_maxM
    U_maxH  = CPL_AtmOcn_bulk_U_maxH
    U_maxE  = CPL_AtmOcn_bulk_U_maxE
    EMIT    = CPL_AtmOcn_bulk_EMIT
    TCW     = CPL_AtmOcn_bulk_TCW
    CM0     = CPL_AtmOcn_bulk_CM0
    visck   = CPL_AtmOcn_bulk_visck
    Z0M_min = CPL_AtmOcn_bulk_Z0M_min
    Z0H_min = CPL_AtmOcn_bulk_Z0H_min
    Z0E_min = CPL_AtmOcn_bulk_Z0E_min
    Z0MI    = CPL_AtmOcn_bulk_Z0MI
    Z0MR    = CPL_AtmOcn_bulk_Z0MR
    Z0MS    = CPL_AtmOcn_bulk_Z0MS
    Z0HI    = CPL_AtmOcn_bulk_Z0HI
    Z0HR    = CPL_AtmOcn_bulk_Z0HR
    Z0HS    = CPL_AtmOcn_bulk_Z0HS
    Z0EI    = CPL_AtmOcn_bulk_Z0EI
    Z0ER    = CPL_AtmOcn_bulk_Z0ER
    Z0ES    = CPL_AtmOcn_bulk_Z0ES

    !--- set up bulk coefficient function
    call CPL_bulkcoef_setup( CPL_TYPE_AtmOcn )

    return
  end subroutine CPL_AtmOcn_bulk_setup

  subroutine CPL_AtmOcn_bulk( &
        SST,                                  & ! (inout)
        SST_UPDATE,                           & ! (in)
        XMFLX, YMFLX, ZMFLX,                  & ! (out)
        SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX,  & ! (out)
        DZ, DENS, MOMX, MOMY, MOMZ,           & ! (in)
        RHOS, PRES, ATMP, QV, SWD, LWD,       & ! (in)
        TW, ALBW, DZW                         ) ! (in)
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    ! parameters
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)

    ! works
    integer :: i, j, n

    real(RP), intent(inout) :: SST(IA,JA) ! sea surface temperature [K]
    logical,  intent(in)    :: SST_UPDATE ! is sea surface temperature updated?

    real(RP), intent(out) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP), intent(out) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP), intent(out) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    real(RP), intent(in) :: DZ  (IA,JA) ! height from the surface to the lowest atmospheric layer [m]
    real(RP), intent(in) :: DENS(IA,JA) ! air density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: MOMX(IA,JA) ! momentum x at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: MOMY(IA,JA) ! momentum y at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: MOMZ(IA,JA) ! momentum z at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: RHOS(IA,JA) ! air density at the sruface [kg/m3]
    real(RP), intent(in) :: PRES(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: ATMP(IA,JA) ! air temperature at the surface [K]
    real(RP), intent(in) :: QV  (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

    real(RP), intent(in) :: TW  (IA,JA) ! water temperature [K]
    real(RP), intent(in) :: ALBW(IA,JA) ! surface albedo in short-wave radiation for water [no unit]
    real(RP), intent(in) :: DZW (IA,JA) ! water depth [m]

    real(RP) :: RES   (IA,JA)
    real(RP) :: DRES  (IA,JA)
    real(RP) :: oldRES(IA,JA) ! RES in previous step
    real(RP) :: redf  (IA,JA) ! reduced factor
    !---------------------------------------------------------------------------

    redf  (:,:) = 1.0_RP
    oldRES(:,:) = 1.0E+5_RP

    do n = 1, nmax
      ! calculate surface flux
      call bulkflux( &
        RES, DRES,                           & ! (out)
        XMFLX, YMFLX, ZMFLX,                 & ! (out)
        SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
        SST, DZ, DENS, MOMX, MOMY, MOMZ,     & ! (in)
        RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
        TW, ALBW, DZW                        ) ! (in)

      if( SST_UPDATE ) then

        do j = JS-1, JE+1
        do i = IS-1, IE+1

          if( redf(i,j) < 0.0_RP ) then
            redf(i,j) = 1.0_RP
          end if

          if( abs(RES(i,j)) > abs(oldRES(i,j)) ) then
            redf(i,j) = max( TFa*redf(i,j), redf_min )
          else
            redf(i,j) = min( TFb*redf(i,j), redf_max )
          end if

          if( DRES(i,j) > 0.0_RP ) then
            redf(i,j) = -1.0_RP
          end if

          ! update surface temperature
          SST(i,j)  = SST(i,j) - redf(i,j) * RES(i,j)/DRES(i,j)

          ! put residual in ocean heat flux
          WHFLX(i,j) = WHFLX(i,j) - RES(i,j)

          ! save residual in this step
          oldRES(i,j) = RES(i,j)

        end do
        end do

        if( maxval(abs(RES(IS-1:IE+1,JS-1:JE+1))) < res_min ) then
          ! iteration converged
          exit
        end if

      else
        ! get surface flux without SST updating
        exit

      end if

    end do

    if( n > nmax ) then
      ! not converged and stop program
      if( IO_L ) write(IO_FID_LOG,*) 'Error: surface tempearture is not converged.'
      call PRC_MPIstop
    end if

    return
  end subroutine CPL_AtmOcn_bulk

  subroutine bulkflux( &
      RES, DRES,                           & ! (out)
      XMFLX, YMFLX, ZMFLX,                 & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
      TS, DZ, DENS, MOMX, MOMY, MOMZ,      & ! (in)
      RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
      TW, ALBW, DZW                        ) ! (in)
    use mod_const, only: &
      GRAV   => CONST_GRAV,  &
      CPdry  => CONST_CPdry, &
      Rvap   => CONST_Rvap,  &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0,   &
      P00    => CONST_PRE00
    use mod_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use mod_cpl_bulkcoef, only: &
      CPL_bulkcoef
    implicit none

    ! argument
    real(RP), intent(out) :: RES   (IA,JA) ! residual in the equation of heat balance
    real(RP), intent(out) :: DRES  (IA,JA) ! d(residual) / d(Ts)
    real(RP), intent(out) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP), intent(out) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP), intent(out) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    real(RP), intent(in) :: TS  (IA,JA) ! skin temperature [K]
    real(RP), intent(in) :: DZ  (IA,JA) ! height from the surface to the lowest atmospheric layer [m]

    real(RP), intent(in) :: DENS(IA,JA) ! air density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: MOMX(IA,JA) ! momentum x at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: MOMY(IA,JA) ! momentum y at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: MOMZ(IA,JA) ! momentum z at the lowest atmospheric layer [kg/m2/s]
    real(RP), intent(in) :: RHOS(IA,JA) ! air density at the sruface [kg/m3]
    real(RP), intent(in) :: PRES(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: ATMP(IA,JA) ! air temperature at the surface [K]
    real(RP), intent(in) :: QV  (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

    real(RP), intent(in) :: TW  (IA,JA) ! water temperature [K]
    real(RP), intent(in) :: ALBW(IA,JA) ! surface albedo in short-wave radiation for water [no unit]
    real(RP), intent(in) :: DZW (IA,JA) ! water depth [m]

    ! work
    real(RP) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP) :: Ustar ! friction velocity [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: dCm, dCh, dCe

    real(RP) :: Z0M, Z0H, Z0E ! oceanic modified roughness length [m]
    real(RP) :: SQV, dSQV ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dLWUFLX, dWHFLX, dSHFLX, dLHFLX

    integer :: i, j
    !---------------------------------------------------------------------------

    ! at (u, y, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( MOMZ(i,j) + MOMZ(i+1,j)                               ) )**2 &
           + ( 2.0_RP *   MOMX(i,j)                                               )**2 &
           + ( 0.5_RP * ( MOMY(i,j-1) + MOMY(i,j) + MOMY(i+1,j-1) + MOMY(i+1,j) ) )**2 &
           ) / ( DENS(i,j) + DENS(i+1,j) )

      !--- friction velocity at u, v, and w points
      Ustar = max ( sqrt ( CM0 ) * Uabs , Ustar_min )

      !--- roughness lengths at u, v, and w points
      Z0M = max( Z0MI + Z0MR/GRAV * Ustar*Ustar + Z0MS*visck / Ustar, Z0M_min )
      Z0H = max( Z0HI + Z0HR/GRAV * Ustar*Ustar + Z0HS*visck / Ustar, Z0H_min )
      Z0E = max( Z0EI + Z0ER/GRAV * Ustar*Ustar + Z0ES*visck / Ustar, Z0E_min )

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                           & ! (out)
          ( ATMP(i,j) + ATMP(i+1,j) ) * 0.5_RP, & ! (in)
          ( TS  (i,j) + TS  (i+1,j) ) * 0.5_RP, & ! (in)
          DZ(i,j), Uabs,                        & ! (in)
          Z0M, Z0H, Z0E                         ) ! (in)

      XMFLX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMX(i,j)
    enddo
    enddo

    ! at (x, v, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( MOMZ(i,j) + MOMZ(i,j+1)                               ) )**2 &
           + ( 0.5_RP * ( MOMX(i-1,j) + MOMX(i,j) + MOMX(i-1,j+1) + MOMX(i,j+1) ) )**2 &
           + ( 2.0_RP *   MOMY(i,j)                                               )**2 &
           ) / ( DENS(i,j) + DENS(i,j+1) )

      !--- friction velocity at u, v, and w points
      Ustar = max ( sqrt ( CM0 ) * Uabs , Ustar_min )

      !--- roughness lengths at u, v, and w points
      Z0M = max( Z0MI + Z0MR/GRAV * Ustar*Ustar + Z0MS*visck / Ustar, Z0M_min )
      Z0H = max( Z0HI + Z0HR/GRAV * Ustar*Ustar + Z0HS*visck / Ustar, Z0H_min )
      Z0E = max( Z0EI + Z0ER/GRAV * Ustar*Ustar + Z0ES*visck / Ustar, Z0E_min )

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                           & ! (out)
          ( ATMP(i,j) + ATMP(i,j+1) ) * 0.5_RP, & ! (in)
          ( TS  (i,j) + TS  (i,j+1) ) * 0.5_RP, & ! (in)
          DZ(i,j), Uabs,                        & ! (in)
          Z0M, Z0H, Z0E                         ) ! (in)

      YMFLX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMY(i,j)
    enddo
    enddo

    ! at cell center
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      Uabs = sqrt( &
             ( MOMZ(i,j)               )**2 &
           + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
           + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
           ) / DENS(i,j) * 0.5_RP

      !--- friction velocity at u, v, and w points
      Ustar = max ( sqrt ( CM0 ) * Uabs , Ustar_min )

      !--- roughness lengths at u, v, and w points
      Z0M = max( Z0MI + Z0MR/GRAV * Ustar*Ustar + Z0MS*visck / Ustar, Z0M_min )
      Z0H = max( Z0HI + Z0HR/GRAV * Ustar*Ustar + Z0HS*visck / Ustar, Z0H_min )
      Z0E = max( Z0EI + Z0ER/GRAV * Ustar*Ustar + Z0ES*visck / Ustar, Z0E_min )

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                  & ! (out)
          ATMP(i,j), TS(i,j),          & ! (in)
          DZ(i,j), Uabs,               & ! (in)
          Z0M, Z0H, Z0E                ) ! (in)

      ZMFLX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMZ(i,j) * 0.5_RP

      ! saturation at the surface
      call qsat( SQV, TS(i,j), PRES(i,j) )

      SHFLX (i,j) = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOS(i,j) * Ch * ( TS(i,j) - ATMP(i,j) )
      LHFLX (i,j) = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOS(i,j) * Ce * ( SQV - QV(i,j) )
      WHFLX (i,j) = -2.0_RP * TCW * ( TS(i,j) - TW(i,j)  ) / DZW(i,j)
      SWUFLX(i,j) = ALBW(i,j) * SWD(i,j)
      LWUFLX(i,j) = EMIT * STB * TS(i,j)**4

      ! calculation for residual
      RES(i,j) = SWD(i,j) - SWUFLX(i,j) + LWD(i,j) - LWUFLX(i,j) - SHFLX(i,j) - LHFLX(i,j) + WHFLX(i,j)

      call CPL_bulkcoef( &
          dCm, dCh, dCe,               & ! (out)
          ATMP(i,j), TS(i,j)+dTS,      & ! (in)
          DZ(i,j), Uabs,               & ! (in)
          Z0M, Z0H, Z0E                ) ! (in)

      call qsat( dSQV, TS(i,j)+dTS, PRES(i,j) )

      dSHFLX  = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOS(i,j) &
              * ( (dCh-Ch)/dTS * ( TS(i,j) - ATMP(i,j) ) + Ch )
      dLHFLX  = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOS(i,j) &
              * ( (dCe-Ce)/dTS * ( SQV - QV(i,j) ) + Ce * (dSQV-SQV)/dTS )
      dWHFLX  = -2.0_RP * TCW / DZW(i,j)
      dLWUFLX = 4.0_RP * EMIT * STB * TS(i,j)**3

      ! calculation for d(residual)/dTS
      DRES(i,j) = - dLWUFLX - dSHFLX - dLHFLX + dWHFLX
    enddo
    enddo

    return
  end subroutine bulkflux

end module mod_cpl_atmos_ocean_bulk
