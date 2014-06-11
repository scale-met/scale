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
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  real(RP), private :: U_minM   =   0.0_RP    ! minimum U_abs for u,v,w
  real(RP), private :: U_minH   =   0.0_RP    !                   T
  real(RP), private :: U_minE   =   0.0_RP    !                   q
  real(RP), private :: U_maxM   = 100.0_RP    ! maximum U_abs for u,v,w
  real(RP), private :: U_maxH   = 100.0_RP    !                   T
  real(RP), private :: U_maxE   = 100.0_RP    !                   q

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_AtmOcn_bulk_setup( CPL_TYPE_AtmOcn )
    use scale_process, only: &
       PRC_MPIstop
    use scale_cpl_bulkcoef, only: &
       CPL_bulkcoef_setup
    implicit none

    character(len=*), intent(in) :: CPL_TYPE_AtmOcn

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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[AtmOcn bulk] / Categ[COUPLER] / Origin[SCALElib]'

    if ( CPL_TYPE_AtmOcn /= 'BULK' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx CPL_TYPE_AtmOcn is not BULK. Check!'
       call PRC_MPIstop
    endif

    CPL_AtmOcn_bulk_U_minM  = U_minM
    CPL_AtmOcn_bulk_U_minH  = U_minH
    CPL_AtmOcn_bulk_U_minE  = U_minE
    CPL_AtmOcn_bulk_U_maxM  = U_maxM
    CPL_AtmOcn_bulk_U_maxH  = U_maxH
    CPL_AtmOcn_bulk_U_maxE  = U_maxE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_ATMOCN_BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_ATMOCN_BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CPL_ATMOCN_BULK)

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

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_bulk( &
        SST,        & ! (inout)
        XMFLX,      & ! (out)
        YMFLX,      & ! (out)
        ZMFLX,      & ! (out)
        SHFLX,      & ! (out)
        LHFLX,      & ! (out)
        WHFLX,      & ! (out)
        U10,        & ! (out)
        V10,        & ! (out)
        T2,         & ! (out)
        Q2,         & ! (out)
        SST_UPDATE, & ! (in)
        RHOA,       & ! (in)
        UA,         & ! (in)
        VA,         & ! (in)
        WA,         & ! (in)
        TMPA,       & ! (in)
        PRSA,       & ! (in)
        QVA,        & ! (in)
        PRSS,       & ! (in)
        SWD,        & ! (in)
        LWD,        & ! (in)
        TW,         & ! (in)
        ALB_SW,     & ! (in)
        ALB_LW,     & ! (in)
        Z0M,        & ! (in)
        Z0H,        & ! (in)
        Z0E         ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
      CPdry  => CONST_CPdry, &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0
    use scale_grid_real, only: &
      Z1 => REAL_Z1
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_cpl_bulkcoef, only: &
      CPL_bulkcoef
    implicit none

    ! arguments
    real(RP), intent(inout) :: SST(IA,JA) ! sea surface temperature [K]

    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: WHFLX(IA,JA) ! water heat flux at the surface [W/m2]
    real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

    logical,  intent(in) :: SST_UPDATE  ! is sea surface temperature updated?

    real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface [W/m2]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface [W/m2]

    real(RP), intent(in) :: TW    (IA,JA) ! water temperature [K]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]

    ! works
    real(RP) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: R10m, R02h, R02e ! lapse rate [0-1]
    real(RP) :: SQV ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j, n
    !---------------------------------------------------------------------------

    if( SST_UPDATE ) then
      ! update surface temperature
      SST(:,:) = TW(:,:)
    else
      ! get surface flux without SST updating
    end if

    ! calculate surface flux
    do j = 1, JA
    do i = 1, IA
      Uabs = sqrt( UA(i,j)**2 + VA(i,j)**2 + WA(i,j)**2 )

      call CPL_bulkcoef( &
          Cm,        & ! (out)
          Ch,        & ! (out)
          Ce,        & ! (out)
          R10m,      & ! (out)
          R02h,      & ! (out)
          R02e,      & ! (out)
          TMPA(i,j), & ! (in)
          SST (i,j), & ! (in)
          PRSA(i,j), & ! (in)
          PRSS(i,j), & ! (in)
          Uabs,      & ! (in)
          Z1  (i,j), & ! (in)
          Z0M (i,j), & ! (in)
          Z0H (i,j), & ! (in)
          Z0E (i,j)  ) ! (in)

      XMFLX(i,j) = -Cm * min(max(Uabs,U_minM),U_maxM) * RHOA(i,j) * UA(i,j)
      YMFLX(i,j) = -Cm * min(max(Uabs,U_minM),U_maxM) * RHOA(i,j) * VA(i,j)
      ZMFLX(i,j) = -Cm * min(max(Uabs,U_minM),U_maxM) * RHOA(i,j) * WA(i,j)

      ! saturation at the surface
      call qsat( SQV, SST(i,j), PRSS(i,j) )

      SHFLX (i,j) = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOA(i,j) * Ch * ( SST(i,j) - TMPA(i,j) )
      LHFLX (i,j) = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOA(i,j) * Ce * ( SQV - QVA(i,j) )

      ! calculation for residual
      WHFLX(i,j) = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) * (-1.0_RP) &
                 - ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * SST(i,j)**4 )&
                 + SHFLX(i,j) + LHFLX(i,j)

      ! diagnositc variables
      U10(i,j) = R10m * UA(i,j)
      V10(i,j) = R10m * VA(i,j)

      T2(i,j) = (          R02h ) * TMPA(i,j) &
              + ( 1.0_RP - R02h ) * SST (i,j)
      Q2(i,j) = (          R02e ) * QVA (i,j) &
              + ( 1.0_RP - R02e ) * SQV

    enddo
    enddo

    return
  end subroutine CPL_AtmOcn_bulk

end module scale_cpl_atmos_ocean_bulk
