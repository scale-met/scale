!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Ocean Surface fluxes
!!
!! @par Description
!!          Surface flux between atmosphere and ocean with constant method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-26 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module scale_cpl_atmos_ocean_const
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
  public :: CPL_AtmOcn_const_setup
  public :: CPL_AtmOcn_const

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
  ! limiter
  real(RP), private, save   :: U_min =   0.0_RP ! minimum U_abs for u,v,w
  real(RP), private, save   :: U_max = 100.0_RP ! maximum U_abs for u,v,w

  real(RP), private, save   :: Cm_min = 1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, save   :: Cm_max = 2.5E-3_RP ! maximum bulk coef. of u,v,w

  real(RP), private, save   :: Const_CM    =   1.1E-3_RP ! constant bulk coef. of u,v,w
  real(RP), private, save   :: Const_SWU   = 200.0_RP    ! constant upward short-wave radiation flux [W/m2]
  real(RP), private, save   :: Const_LWU   = 200.0_RP    ! constant upward long-wave radiation flux [W/m2]
  real(RP), private, save   :: Const_SH    =  15.0_RP    ! constant sensible heat flux [W/m2]
  real(RP), private, save   :: Const_LH    = 115.0_RP    ! constant latent heat flux [W/m2]
  real(RP), private, save   :: Const_WH    =   0.0_RP    ! constant water heat flux [W/m2]
  real(RP), private, save   :: Const_Ustar =   0.25_RP   ! constant friction velocity [m/s]
  real(RP), private, save   :: Const_FREQ  =  24.0_RP    ! frequency of sensible heat flux [hour]

  integer(4), private, save :: CMTYPE = 0 ! 0->Bulk coef. is constant
                                          ! 1->Friction velocity is constant

  logical, private, save    :: DIURNAL = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_const_setup( CPL_TYPE_AtmOcn )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: CPL_TYPE_AtmOcn

    real(RP) :: CPL_AtmOcn_const_U_min
    real(RP) :: CPL_AtmOcn_const_U_max
    real(RP) :: CPL_AtmOcn_const_CM_min
    real(RP) :: CPL_AtmOcn_const_CM_max
    real(RP) :: CPL_AtmOcn_const_CM
    real(RP) :: CPL_AtmOcn_const_SWU
    real(RP) :: CPL_AtmOcn_const_LWU
    real(RP) :: CPL_AtmOcn_const_SH
    real(RP) :: CPL_AtmOcn_const_LH
    real(RP) :: CPL_AtmOcn_const_WH
    real(RP) :: CPL_AtmOcn_const_Ustar
    real(RP) :: CPL_AtmOcn_const_FREQ
    integer  :: CPL_AtmOcn_const_CMTYPE
    logical  :: CPL_AtmOcn_const_DIURNAL

    NAMELIST / PARAM_CPL_ATMOCN_CONST / &
       CPL_AtmOcn_const_U_min,  &
       CPL_AtmOcn_const_U_max,  &
       CPL_AtmOcn_const_CM_min, &
       CPL_AtmOcn_const_CM_max, &
       CPL_AtmOcn_const_CM,     &
       CPL_AtmOcn_const_SWU,    &
       CPL_AtmOcn_const_LWU,    &
       CPL_AtmOcn_const_SH,     &
       CPL_AtmOcn_const_LH,     &
       CPL_AtmOcn_const_WH,     &
       CPL_AtmOcn_const_Ustar,  &
       CPL_AtmOcn_const_FREQ,   &
       CPL_AtmOcn_const_CMTYPE, &
       CPL_AtmOcn_const_DIURNAL

    integer :: ierr
    !---------------------------------------------------------------------------

    CPL_AtmOcn_const_U_min   = U_min
    CPL_AtmOcn_const_U_max   = U_max
    CPL_AtmOcn_const_CM_min  = CM_min
    CPL_AtmOcn_const_CM_max  = CM_max
    CPL_AtmOcn_const_CM      = Const_CM
    CPL_AtmOcn_const_SWU     = Const_SWU
    CPL_AtmOcn_const_LWU     = Const_LWU
    CPL_AtmOcn_const_SH      = Const_SH
    CPL_AtmOcn_const_LH      = Const_LH
    CPL_AtmOcn_const_WH      = Const_WH
    CPL_AtmOcn_const_Ustar   = Const_Ustar
    CPL_AtmOcn_const_FREQ    = Const_FREQ
    CPL_AtmOcn_const_CMTYPE  = CMTYPE
    CPL_AtmOcn_const_DIURNAL = DIURNAL

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean: constant flux parameter'

    if ( CPL_TYPE_AtmOcn /= 'CONST' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx CPL_TYPE_AtmOcn is not CONST. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_ATMOCN_CONST,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_ATMOCN_CONST. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL_ATMOCN_CONST)

    U_min       = CPL_AtmOcn_const_U_min
    U_max       = CPL_AtmOcn_const_U_max
    CM_min      = CPL_AtmOcn_const_CM_min
    CM_max      = CPL_AtmOcn_const_CM_max
    Const_Cm    = CPL_AtmOcn_const_Cm
    Const_SWU   = CPL_AtmOcn_const_SWU
    Const_LWU   = CPL_AtmOcn_const_LWU
    Const_SH    = CPL_AtmOcn_const_SH
    Const_LH    = CPL_AtmOcn_const_LH
    Const_WH    = CPL_AtmOcn_const_WH
    Const_Ustar = CPL_AtmOcn_const_Ustar
    Const_FREQ  = CPL_AtmOcn_const_FREQ
    CMTYPE      = CPL_AtmOcn_const_CMTYPE
    DIURNAL     = CPL_AtmOcn_const_DIURNAL

    return
  end subroutine CPL_AtmOcn_const_setup

  subroutine CPL_AtmOcn_const( &
        SST,                                  & ! (inout)
        XMFLX, YMFLX, ZMFLX,                  & ! (out)
        SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX,  & ! (out)
        SST_UPDATE,                           & ! (in)
        DZ, DENS, MOMX, MOMY, MOMZ,           & ! (in)
        RHOS, PRES, ATMP, QV, SWD, LWD,       & ! (in)
        TW, ALB_SW, ALB_LW,                   & ! (in)
        Z0M, Z0H, Z0E                         ) ! (in)
    use scale_const, only: &
      PI => CONST_PI
    use scale_time, only: &
      TIME => TIME_NOWDAYSEC
    implicit none

    ! argument
    real(RP), intent(inout) :: SST (IA,JA) ! ocean surface temperature [K]

    real(RP), intent(out) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP), intent(out) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP), intent(out) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    logical,  intent(in) :: SST_UPDATE  ! is ocean surface temperature updated?

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

    ! work
    real(RP) :: Uabs  ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP) :: Cm    ! bulk coefficient [no unit]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
      ! at cell center
      Uabs = sqrt( &
             ( MOMZ(i,j)               )**2 &
           + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
           + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
           ) / DENS(i,j) * 0.5_RP

      if( CMTYPE == 1 ) then
        ! friction velocity is constant
        Cm = min( max( Const_Ustar**2 / Uabs**2, Cm_min ), Cm_max )
      else
        ! Bulk coef. is constant
        Cm = Const_Cm
      endif

      XMFLX(i,j) = -Cm * min(max(Uabs,U_min),U_max) * MOMX(i,j)
      YMFLX(i,j) = -Cm * min(max(Uabs,U_min),U_max) * MOMY(i,j)
      ZMFLX(i,j) = -Cm * min(max(Uabs,U_min),U_max) * MOMZ(i,j)

      if( DIURNAL ) then
        ! include diurnal change
        SHFLX (i,j) = Const_SH  * sin( TIME / ( Const_FREQ*3600.0_RP )*2.0_RP*PI )
        SWUFLX(i,j) = max( Const_SWU * sin( TIME / ( Const_FREQ*3600.0_RP )*2.0_RP*PI ), 0.0_RP )
      else
        SHFLX (i,j) = Const_SH
        SWUFLX(i,j) = Const_SWU
      endif

      LHFLX (i,j) = Const_LH
      LWUFLX(i,j) = Const_LWU
      WHFLX (i,j) = Const_WH

    enddo
    enddo

    return
  end subroutine CPL_AtmOcn_const

end module scale_cpl_atmos_ocean_const
