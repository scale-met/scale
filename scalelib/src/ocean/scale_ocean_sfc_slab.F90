!-------------------------------------------------------------------------------
!> module OCEAN / Surface flux with slab ocean model
!!
!! @par Description
!!          Surface flux with slab ocean model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_sfc_slab
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
  public :: OCEAN_SFC_SLAB_setup
  public :: OCEAN_SFC_SLAB
  public :: OCEAN_SFC_SLAB_SimpleAlbedo

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
  logical, private :: SST_UPDATE

  logical, allocatable, private :: is_OCN(:,:)

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_SFC_SLAB_setup( OCEAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLAB] / Categ[OCEAN SFC] / Origin[SCALElib]'

    select case( OCEAN_TYPE )
    case( 'CONST' )
       SST_UPDATE = .false.
    case( 'SLAB', 'FILE' )
       SST_UPDATE = .true.
    case default
       write(*,*) 'xxx wrong OCEAN_TYPE. Check!'
       call PRC_MPIstop
    end select

    ! judge to run slab ocean model
    allocate( is_OCN(IA,JA) )

    do j = JS, JE
    do i = IS, IE
      if( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
        is_OCN(i,j) = .true.
      else
        is_OCN(i,j) = .false.
      end if
    end do
    end do

    return
  end subroutine OCEAN_SFC_SLAB_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_SFC_SLAB( &
        SST_t,  & ! [OUT]
        ZMFLX,  & ! [OUT]
        XMFLX,  & ! [OUT]
        YMFLX,  & ! [OUT]
        SHFLX,  & ! [OUT]
        LHFLX,  & ! [OUT]
        WHFLX,  & ! [OUT]
        U10,    & ! [OUT]
        V10,    & ! [OUT]
        T2,     & ! [OUT]
        Q2,     & ! [OUT]
        TMPA,   & ! [IN]
        PRSA,   & ! [IN]
        WA,     & ! [IN]
        UA,     & ! [IN]
        VA,     & ! [IN]
        RHOA,   & ! [IN]
        QVA,    & ! [IN]
        Z1,     & ! [IN]
        PBL,    & ! [IN]
        PRSS,   & ! [IN]
        LWD,    & ! [IN]
        SWD,    & ! [IN]
        TW,     & ! [IN]
        SST,    & ! [IN]
        ALB_LW, & ! [IN]
        ALB_SW, & ! [IN]
        Z0M,    & ! [IN]
        Z0H,    & ! [IN]
        Z0E,    & ! [IN]
        dt      ) ! [IN]
    use scale_grid_index
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
      Rdry  => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      STB   => CONST_STB
    use scale_bulkflux, only: &
      BULKFLUX
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    implicit none

    ! arguments
    real(RP), intent(out) :: SST_t(IA,JA) ! tendency of sea surface temperature
    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: WHFLX(IA,JA) ! water heat flux at the surface [W/m2]
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
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface [W/m2]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface [W/m2]

    real(RP), intent(in) :: TW    (IA,JA) ! water temperature [K]
    real(RP), intent(in) :: SST   (IA,JA) ! sea surface temperature [K]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
    real(DP), intent(in) :: dt            ! delta time

    ! works
    real(RP) :: SST1(IA,JA)

    real(RP) :: Ustar, Ustar10, Ustar2 ! friction velocity [m]
    real(RP) :: Tstar, Tstar10, Tstar2 ! friction temperature [K]
    real(RP) :: Qstar, Qstar10, Qstar2 ! friction mixing rate [kg/kg]
    real(RP) :: Uabs,  Uabs10,  Uabs2  ! modified absolute velocity [m/s]

    real(RP) :: QVsat        ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: LHV(IA,JA)   ! latent heat of vaporization [J/kg]
    real(RP) :: CVTH

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean surface step: Slab'

    ! update surface temperature
    do j = JS, JE
    do i = IS, IE
      SST1(i,j) = TW(i,j) ! assumed well-mixed condition
    end do
    end do

    call HYDROMETEOR_LHV( LHV, TMPA )

    do j = JS, JE
    do i = IS, IE

      if( is_OCN(i,j) ) then

        ! saturation at the surface
        call qsat( QVsat,     & ! [OUT]
                   SST1(i,j), & ! [IN]
                   PRSS(i,j)  ) ! [IN]

        call BULKFLUX( &
            Ustar,     & ! [OUT]
            Tstar,     & ! [OUT]
            Qstar,     & ! [OUT]
            Uabs,      & ! [OUT]
            TMPA(i,j), & ! [IN]
            SST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVsat,     & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            Z1  (i,j), & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ! for 10m wind
        call BULKFLUX( &
            Ustar10,   & ! [OUT]
            Tstar10,   & ! [OUT]
            Qstar10,   & ! [OUT]
            Uabs10,    & ! [OUT]
            TMPA(i,j), & ! [IN]
            SST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVsat,     & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            10.0_RP,   & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ! for 2m temperature / mixing ratio
        call BULKFLUX( &
            Ustar2,    & ! [OUT]
            Tstar2,    & ! [OUT]
            Qstar2,    & ! [OUT]
            Uabs2,     & ! [OUT]
            TMPA(i,j), & ! [IN]
            SST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVsat,     & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            2.0_RP,    & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ZMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * WA(i,j)
        XMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * UA(i,j)
        YMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * VA(i,j)

        SHFLX (i,j) = -CPdry    * RHOA(i,j) * Ustar * Tstar
        LHFLX (i,j) = -LHV(i,j) * RHOA(i,j) * Ustar * Qstar

        ! calculation for residual
        WHFLX(i,j) = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) * (-1.0_RP) &
                   - ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * SST1(i,j)**4 ) &
                   + SHFLX(i,j) + LHFLX(i,j)

        ! diagnositc variables
        CVTH = ( PRSS(i,j) / PRSA(i,j) ) ** ( Rdry / Cpdry )

        U10(i,j) = UA(i,j) * Ustar / Ustar10 * Uabs10 / Uabs
        V10(i,j) = VA(i,j) * Ustar / Ustar10 * Uabs10 / Uabs
        T2 (i,j) = ( 1.0_RP - Tstar / Tstar2 ) * SST1(i,j) + Tstar / Tstar2 * TMPA(i,j) * CVTH
        Q2 (i,j) = ( 1.0_RP - Qstar / Qstar2 ) * QVsat     + Qstar / Qstar2 * QVA (i,j)

      else
        ! not calculate surface flux
        ZMFLX(i,j) = 0.0_RP
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
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

    ! calculate tendency
    if( SST_UPDATE ) then
      do j = JS, JE
      do i = IS, IE
        SST_t(i,j) = ( SST1(i,j) - SST(i,j) ) / dt
      enddo
      enddo
    else
      do j = JS, JE
      do i = IS, IE
        SST_t(i,j) = 0.0_RP
      enddo
      enddo
    end if

    return
  end subroutine OCEAN_SFC_SLAB

  !-----------------------------------------------------------------------------
  subroutine OCEAN_SFC_SLAB_SimpleAlbedo( &
       SFC_albedo_t, &
       SFC_albedo,   &
       cosSZA,       &
       dt            )
    use scale_grid_index
    use scale_const, only: &
       I_SW  => CONST_I_SW, &
       I_LW  => CONST_I_LW
    implicit none

    ! arguments
    real(RP), intent(out) :: SFC_albedo_t(IA,JA,2) ! tendency of sea surface albedo [0-1]
    real(RP), intent(in)  :: SFC_albedo  (IA,JA,2) ! sea surface                    [0-1]
    real(RP), intent(in)  :: cosSZA      (IA,JA)   ! cos(solar zenith angle)        [0-1]
    real(DP), intent(in)  :: dt                    ! delta time

    ! works
    real(RP) :: SFC_albedo1(IA,JA,2)
    real(RP) :: am1

    real(RP), parameter :: c_ocean_albedo(3) = (/ -0.747900_RP, &
                                                  -4.677039_RP, &
                                                   1.583171_RP  /)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       am1 = max( min( cosSZA(i,j), 0.961_RP ), 0.0349_RP )

       ! SFC_albedo1(i,j,I_LW) = 0.5_RP do nothing
       SFC_albedo1(i,j,I_SW) = exp( ( c_ocean_albedo(3)*am1 + c_ocean_albedo(2) )*am1 + c_ocean_albedo(1) )
    enddo
    enddo

    ! calculate tendency
    do j = JS, JE
    do i = IS, IE
       SFC_albedo_t(i,j,I_LW) = 0.0_RP
       SFC_albedo_t(i,j,I_SW) = ( SFC_albedo1(i,j,I_SW) - SFC_albedo(i,j,I_SW) ) / dt
    enddo
    enddo

    return
  end subroutine OCEAN_SFC_SLAB_SimpleAlbedo

end module scale_ocean_sfc_slab
