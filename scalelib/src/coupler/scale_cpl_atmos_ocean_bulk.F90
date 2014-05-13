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
module scale_cpl_atmos_ocean_bulk
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
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
  real(RP), private, save :: U_minM   =   0.0_RP    ! minimum U_abs for u,v,w
  real(RP), private, save :: U_minH   =   0.0_RP    !                   T
  real(RP), private, save :: U_minE   =   0.0_RP    !                   q
  real(RP), private, save :: U_maxM   = 100.0_RP    ! maximum U_abs for u,v,w
  real(RP), private, save :: U_maxH   = 100.0_RP    !                   T
  real(RP), private, save :: U_maxE   = 100.0_RP    !                   q

contains
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_bulk_setup( CPL_TYPE_AtmOcn )
    use scale_process, only: &
       PRC_MPIstop
    use scale_ocean_roughness, only: &
       OCEAN_roughness_setup
    use scale_cpl_bulkcoef, only: &
       CPL_bulkcoef_setup
    implicit none

    character(len=H_SHORT), intent(in) :: CPL_TYPE_AtmOcn

    real(RP) :: CPL_AtmOcn_bulk_U_minM
    real(RP) :: CPL_AtmOcn_bulk_U_minH
    real(RP) :: CPL_AtmOcn_bulk_U_minE
    real(RP) :: CPL_AtmOcn_bulk_U_maxM
    real(RP) :: CPL_AtmOcn_bulk_U_maxH
    real(RP) :: CPL_AtmOcn_bulk_U_maxE

    NAMELIST / PARAM_CPL_ATMOCN_BULK / &
       CPL_AtmOcn_bulk_U_minM,  &
       CPL_AtmOcn_bulk_U_minH,  &
       CPL_AtmOcn_bulk_U_minE,  &
       CPL_AtmOcn_bulk_U_maxM,  &
       CPL_AtmOcn_bulk_U_maxH,  &
       CPL_AtmOcn_bulk_U_maxE

    integer :: ierr
    !---------------------------------------------------------------------------

    CPL_AtmOcn_bulk_U_minM  = U_minM
    CPL_AtmOcn_bulk_U_minH  = U_minH
    CPL_AtmOcn_bulk_U_minE  = U_minE
    CPL_AtmOcn_bulk_U_maxM  = U_maxM
    CPL_AtmOcn_bulk_U_maxH  = U_maxH
    CPL_AtmOcn_bulk_U_maxE  = U_maxE

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean: bulk flux parameter'

    if ( CPL_TYPE_AtmOcn /= 'BULK' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx CPL_TYPE_AtmOcn is not BULK. Check!'
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

    U_minM    = CPL_AtmOcn_bulk_U_minM
    U_minH    = CPL_AtmOcn_bulk_U_minH
    U_minE    = CPL_AtmOcn_bulk_U_minE
    U_maxM    = CPL_AtmOcn_bulk_U_maxM
    U_maxH    = CPL_AtmOcn_bulk_U_maxH
    U_maxE    = CPL_AtmOcn_bulk_U_maxE

    !--- set up bulk coefficient function
    call CPL_bulkcoef_setup

    return
  end subroutine CPL_AtmOcn_bulk_setup

  subroutine CPL_AtmOcn_bulk( &
        SST,                                  & ! (inout)
        XMFLX, YMFLX, ZMFLX,                  & ! (out)
        SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX,  & ! (out)
        SST_UPDATE,                           & ! (in)
        DZ, DENS, MOMX, MOMY, MOMZ,           & ! (in)
        RHOS, PRES, ATMP, QV, SWD, LWD,       & ! (in)
        TW, ALB_SW, ALB_LW,                   & ! (in)
        Z0M, Z0H, Z0E                         ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    ! parameters
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)

    ! works
    integer :: i, j, n

    real(RP), intent(inout) :: SST (IA,JA) ! sea surface temperature [K]

    real(RP), intent(out) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP), intent(out) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP), intent(out) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    logical,  intent(in) :: SST_UPDATE  ! is sea surface temperature updated?

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

    real(RP), intent(in) :: TW    (IA,JA) ! water temperature [K]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]

    !---------------------------------------------------------------------------

    if( SST_UPDATE ) then
      ! update surface temperature
      SST(:,:) = TW(:,:)
    else
      ! get surface flux without SST updating
    end if

    ! calculate surface flux
    call bulkflux( &
      XMFLX, YMFLX, ZMFLX,                 & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
      SST, DZ, DENS, MOMX, MOMY, MOMZ,     & ! (in)
      RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
      ALB_SW, ALB_LW, Z0M, Z0H, Z0E        ) ! (in)

    return
  end subroutine CPL_AtmOcn_bulk

  subroutine bulkflux( &
      XMFLX, YMFLX, ZMFLX,                 & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
      TS, DZ, DENS, MOMX, MOMY, MOMZ,      & ! (in)
      RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
      ALB_SW, ALB_LW, Z0M, Z0H, Z0E        ) ! (in)
    use scale_const, only: &
      GRAV   => CONST_GRAV,  &
      CPdry  => CONST_CPdry, &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_cpl_bulkcoef, only: &
      CPL_bulkcoef
    implicit none

    ! argument
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

    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]

    ! work
    real(RP) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP) :: Ustar ! friction velocity [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]

    real(RP) :: SQV ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    ! at cell center
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      Uabs = sqrt( &
             ( MOMZ(i,j)               )**2 &
           + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
           + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
           ) / DENS(i,j) * 0.5_RP

      call CPL_bulkcoef( &
          Cm, Ch, Ce,                  & ! (out)
          ATMP(i,j), TS(i,j),          & ! (in)
          DZ(i,j), Uabs,               & ! (in)
          Z0M(i,j), Z0H(i,j), Z0E(i,j) ) ! (in)

      ! saturation at the surface
      call qsat( SQV, TS(i,j), PRES(i,j) )

      XMFLX (i,j) = -Cm * min(max(Uabs,U_minM),U_maxM) * MOMX(i,j)
      YMFLX (i,j) = -Cm * min(max(Uabs,U_minM),U_maxM) * MOMY(i,j)
      ZMFLX (i,j) = -Cm * min(max(Uabs,U_minM),U_maxM) * MOMZ(i,j)

      SHFLX (i,j) = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOS(i,j) * Ch * ( TS(i,j) - ATMP(i,j) )
      LHFLX (i,j) = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOS(i,j) * Ce * ( SQV - QV(i,j) )
      SWUFLX(i,j) = ALB_SW(i,j) * SWD(i,j)
      LWUFLX(i,j) = ALB_LW(i,j) * LWD(i,j) + ( 1.0_RP - ALB_LW(i,j) ) * STB * TS(i,j)**4

      ! calculation for residual
      WHFLX(i,j) = - SWD(i,j) + SWUFLX(i,j) - LWD(i,j) + LWUFLX(i,j) + SHFLX(i,j) + LHFLX(i,j)
    enddo
    enddo

    return
  end subroutine bulkflux

end module scale_cpl_atmos_ocean_bulk
