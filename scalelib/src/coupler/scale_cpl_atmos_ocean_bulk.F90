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
  logical, allocatable, private :: is_FLX(:,:) ! is atmos-land coupler run?

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_AtmOcn_bulk_setup( CPL_TYPE_AtmOcn )
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    use scale_cpl_bulkflux, only: &
       CPL_bulkflux_setup
    implicit none

    character(len=*), intent(in) :: CPL_TYPE_AtmOcn

    logical :: dummy = .false.

    NAMELIST / PARAM_CPL_ATMOCN_BULK / &
      dummy

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[AtmOcn bulk] / Categ[COUPLER] / Origin[SCALElib]'

    if ( CPL_TYPE_AtmOcn /= 'BULK' ) then
       write(*,*) 'xxx CPL_TYPE_AtmOcn is not BULK. Check!'
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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CPL_ATMOCN_BULK)

    !--- judge to run atmos-ocean coupler
    allocate( is_FLX(IA,JA) )

    do j = JS, JE
    do i = IS, IE
      if( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
        is_FLX(i,j) = .true.
      else
        is_FLX(i,j) = .false.
      end if
    end do
    end do

    !--- set up bulk coefficient function
    call CPL_bulkflux_setup

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
        PBL,        & ! (in)
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
      CPdry => CONST_CPdry, &
      STB   => CONST_STB,   &
      LHV0  => CONST_LHV0
    use scale_grid_real, only: &
      Z1 => REAL_Z1
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_cpl_bulkflux, only: &
      CPL_bulkflux
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
    real(RP), intent(in) :: PBL (IA,JA) ! the top of atmospheric mixing layer [m]
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
    real(RP) :: Ustar ! friction velocity [m]
    real(RP) :: Tstar ! friction temperature [K]
    real(RP) :: Qstar ! friction mixing rate [kg/kg]
    real(RP) :: Uabs ! modified absolute velocity [m/s]
    real(RP) :: SQV ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j, n
    !---------------------------------------------------------------------------

    if( SST_UPDATE ) then
      ! update surface temperature
      SST(:,:) = TW(:,:)
    else
      ! get surface flux without SST updating
    end if

    do j = JS, JE
    do i = IS, IE

      if( is_FLX(i,j) ) then
        ! saturation at the surface
        call qsat( SQV, SST(i,j), PRSS(i,j) )

        call CPL_bulkflux( &
            Ustar,     & ! (out)
            Tstar,     & ! (out)
            Qstar,     & ! (out)
            Uabs,      & ! (out)
            TMPA(i,j), & ! (in)
            SST (i,j), & ! (in)
            PRSA(i,j), & ! (in)
            PRSS(i,j), & ! (in)
            QVA (i,j), & ! (in)
            SQV,       & ! (in)
            UA  (i,j), & ! (in)
            VA  (i,j), & ! (in)
            Z1  (i,j), & ! (in)
            PBL (i,j), & ! (in)
            Z0M (i,j), & ! (in)
            Z0H (i,j), & ! (in)
            Z0E (i,j)  ) ! (in)

        XMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * UA(i,j)
        YMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * VA(i,j)
        ZMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * WA(i,j)

        SHFLX (i,j) = -CPdry * RHOA(i,j) * Ustar * Tstar
        LHFLX (i,j) = -LHV0  * RHOA(i,j) * Ustar * Qstar

        ! calculation for residual
        WHFLX(i,j) = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) * (-1.0_RP) &
                   - ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * SST(i,j)**4 )&
                   + SHFLX(i,j) + LHFLX(i,j)

        ! diagnositc variables
        U10(i,j) = UA (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        V10(i,j) = VA (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        T2 (i,j) = SST(i,j) + ( TMPA(i,j) - SST(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                       / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
        Q2 (i,j) = SQV      + (  QVA(i,j) - SQV      ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                       / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

      else
        ! not calculate surface flux
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        ZMFLX(i,j) = 0.0_RP
        SHFLX(i,j) = 0.0_RP
        LHFLX(i,j) = 0.0_RP
        WHFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP
      end if

    enddo
    enddo

    return
  end subroutine CPL_AtmOcn_bulk

end module scale_cpl_atmos_ocean_bulk
