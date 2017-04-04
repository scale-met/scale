!-------------------------------------------------------------------------------
!> module LAND / Surface fluxes with thick-slab land model
!!
!! @par Description
!!          Surface flux with thick-slab land model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_land_sfc_thick_slab
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
  public :: LAND_SFC_THICK_SLAB_setup
  public :: LAND_SFC_THICK_SLAB

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
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_SFC_THICK_SLAB_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    logical :: dummy

    NAMELIST / PARAM_LAND_SFC_THICK_SLAB / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[THICK-SLAB] / Categ[LAND SFC] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SFC_THICK_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_SFC_THICK_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_SFC_THICK_SLAB)

    return
  end subroutine LAND_SFC_THICK_SLAB_setup

  !-----------------------------------------------------------------------------
  subroutine LAND_SFC_THICK_SLAB( &
        LST_t,      &
        ZMFLX,      &
        XMFLX,      &
        YMFLX,      &
        SHFLX,      &
        LHFLX,      &
        GHFLX,      &
        U10,        &
        V10,        &
        T2,         &
        Q2,         &
        TMPA,       &
        PRSA,       &
        WA,         &
        UA,         &
        VA,         &
        RHOA,       &
        QVA,        &
        Z1,         &
        PBL,        &
        PRSS,       &
        LWD,        &
        SWD,        &
        TG,         &
        LST,        &
        QVEF,       &
        ALB_LW,     &
        ALB_SW,     &
        DZG,        &
        TCS,        &
        Z0M,        &
        Z0H,        &
        Z0E,        &
        dt          )
    use scale_process, only: &
      PRC_myrank,  &
      PRC_MPIstop
    use scale_const, only: &
      Rdry  => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      STB   => CONST_STB
    use scale_landuse, only: &
      LANDUSE_fact_land
    use scale_atmos_hydrometeor, only: &
      HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
      BULKFLUX
    implicit none

    ! arguments
    real(RP), intent(out) :: LST_t(IA,JA) ! tendency of LST
    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

    real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Z1  (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL (IA,JA) ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface [J/m2/s]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface [J/m2/s]

    real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
    real(RP), intent(in) :: LST   (IA,JA) ! land surface temperature [K]
    real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
    real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [J/m/K/s]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
    real(DP), intent(in) :: dt            ! delta time

    ! works
    real(RP) :: LST1(IA,JA)

    real(RP) :: res    ! residual
    real(RP) :: dres   ! d(residual)/dLST
    real(RP) :: oldres ! residual in previous step
    real(RP) :: redf   ! reduced factor

    real(RP) :: Ustar, dUstar ! friction velocity [m]
    real(RP) :: Tstar, dTstar ! friction potential temperature [K]
    real(RP) :: Qstar, dQstar ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Uabs,  dUabs  ! modified absolute velocity [m/s]
    real(RP) :: QVsat, dQVsat ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS, dQVS     ! water vapor mixing ratio at surface [kg/kg]

    real(RP) :: FracU10 ! calculation parameter for U10 [-]
    real(RP) :: FracT2  ! calculation parameter for T2 [-]
    real(RP) :: FracQ2  ! calculation parameter for Q2 [-]

    real(RP) :: LHV(IA,JA)    ! latent heat of vaporization [J/kg]

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land surface step: Thick-Slab'

    call HYDROMETEOR_LHV( LHV(:,:), TMPA(:,:) )

    ! update surface temperature
    do j = JS, JE
    do i = IS, IE
      LST1 (i,j) = TG(i,j)
      LST_t(i,j) = ( LST1(i,j) - LST(i,j) ) / dt
    end do
    end do

    ! calculate surface flux
    do j = JS, JE
    do i = IS, IE

      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then

        call qsat( QVsat,     & ! [OUT]
                   LST1(i,j), & ! [IN]
                   PRSS(i,j)  ) ! [IN]

        QVS  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat

        call BULKFLUX( &
            Ustar,     & ! [OUT]
            Tstar,     & ! [OUT]
            Qstar,     & ! [OUT]
            Uabs,      & ! [OUT]
            FracU10,   & ! [OUT]
            FracT2,    & ! [OUT]
            FracQ2,    & ! [OUT]
            TMPA(i,j), & ! [IN]
            LST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVS,       & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            Z1  (i,j), & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ZMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * WA(i,j)
        XMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * UA(i,j)
        YMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * VA(i,j)
        SHFLX(i,j) = -CPdry    * RHOA(i,j) * Ustar * Tstar
        LHFLX(i,j) = -LHV(i,j) * RHOA(i,j) * Ustar * Qstar
        GHFLX(i,j) = -2.0_RP * TCS(i,j) * ( LST1(i,j) - TG(i,j) ) / DZG(i,j)

        ! calculation for residual
        res = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
            + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST1(i,j)**4 ) &
            - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

        ! put residual in ground heat flux
        GHFLX(i,j) = GHFLX(i,j) - res

        ! diagnostic variables considering unstable/stable state
        !U10(i,j) = FracU10 * UA(i,j)
        !V10(i,j) = FracU10 * VA(i,j)
        !T2 (i,j) = ( 1.0_RP - FracT2 ) * LST1(i,j) + FracT2 * TMPA(i,j)
        !Q2 (i,j) = ( 1.0_RP - FracQ2 ) * QVS       + FracQ2 * QVA (i,j)

        ! diagnostic variables for neutral state
        U10(i,j) = UA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        V10(i,j) = VA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        T2 (i,j) = LST1(i,j) + ( TMPA(i,j) - LST1(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
        Q2 (i,j) = QVS       + (  QVA(i,j) - QVS       ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

      else

        ! not calculate surface flux
        ZMFLX(i,j) = 0.0_RP
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        SHFLX(i,j) = 0.0_RP
        LHFLX(i,j) = 0.0_RP
        GHFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP

      end if

    end do
    end do

    return
  end subroutine LAND_SFC_THICK_SLAB

end module scale_land_sfc_thick_slab
