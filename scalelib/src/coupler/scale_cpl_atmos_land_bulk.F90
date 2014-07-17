!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Land Surface fluxes
!!
!! @par Description
!!          Surface flux between atmosphere and land with Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module scale_cpl_atmos_land_bulk
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
  public :: CPL_AtmLnd_bulk_setup
  public :: CPL_AtmLnd_bulk

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
  integer,  private :: nmax = 100 ! maximum iteration number

  real(RP), private :: res_min  = 1.0E+0_RP ! minimum number of residual
  real(RP), private :: dres_lim = 1.0E+2_RP ! limit factor of d(residual)/dTS
  real(RP), private :: dTS      = 1.0E-8_RP ! delta surface temp.

  logical, allocatable, private :: is_FLX(:,:) ! is atmos-land coupler run?

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_AtmLnd_bulk_setup( CPL_TYPE_AtmLnd )
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_fact_land
    use scale_cpl_bulkflux, only: &
       CPL_bulkflux_setup
    implicit none

    character(len=*), intent(in) :: CPL_TYPE_AtmLnd

    integer  :: CPL_AtmLnd_bulk_nmax
    real(RP) :: CPL_AtmLnd_bulk_res_min
    real(RP) :: CPL_AtmLnd_bulk_dres_lim
    real(RP) :: CPL_AtmLnd_bulk_dTS

    NAMELIST / PARAM_CPL_ATMLND_BULK / &
       CPL_AtmLnd_bulk_nmax,     &
       CPL_AtmLnd_bulk_res_min,  &
       CPL_AtmLnd_bulk_dres_lim, &
       CPL_AtmLnd_bulk_dTS

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[AtmLnd bulk] / Categ[COUPLER] / Origin[SCALElib]'

    if ( CPL_TYPE_AtmLnd /= 'BULK' ) then
       write(*,*) 'xxx CPL_TYPE_AtmLnd is not BULK. Check!'
       call PRC_MPIstop
    endif

    CPL_AtmLnd_bulk_nmax     = nmax
    CPL_AtmLnd_bulk_res_min  = res_min
    CPL_AtmLnd_bulk_dres_lim = dres_lim
    CPL_AtmLnd_bulk_dTS      = dTS

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_ATMLND_BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_ATMLND_BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CPL_ATMLND_BULK)

    nmax     = CPL_AtmLnd_bulk_nmax
    res_min  = CPL_AtmLnd_bulk_res_min
    dres_lim = CPL_AtmLnd_bulk_dres_lim
    dTS      = CPL_AtmLnd_bulk_dTS

    !--- judge to run atmos-land coupler
    allocate( is_FLX(IA,JA) )

    do j = 1, JA
    do i = 1, IA
      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then
        is_FLX(i,j) = .true.
      else
        is_FLX(i,j) = .false.
      end if
    end do
    end do

    !--- set up bulk coefficient function
    call CPL_bulkflux_setup

    return
  end subroutine CPL_AtmLnd_bulk_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmLnd_bulk( &
        LST,        & ! (inout)
        XMFLX,      & ! (out)
        YMFLX,      & ! (out)
        ZMFLX,      & ! (out)
        SHFLX,      & ! (out)
        LHFLX,      & ! (out)
        GHFLX,      & ! (out)
        U10,        & ! (out)
        V10,        & ! (out)
        T2,         & ! (out)
        Q2,         & ! (out)
        LST_UPDATE, & ! (in)
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
        TG,         & ! (in)
        QVEF,       & ! (in)
        ALB_SW,     & ! (in)
        ALB_LW,     & ! (in)
        TCS,        & ! (in)
        DZG,        & ! (in)
        Z0M,        & ! (in)
        Z0H,        & ! (in)
        Z0E         ) ! (in)
    use scale_const, only: &
      CPdry  => CONST_CPdry, &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0
    use scale_grid_real, only: &
      Z1 => REAL_Z1
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_cpl_bulkflux, only: &
      CPL_bulkflux
    implicit none

    ! parameters
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)

    ! arguments
    real(RP), intent(inout) :: LST(IA,JA) ! land surface temperature [K]

    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [W/m2]
    real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

    logical,  intent(in) :: LST_UPDATE  ! is land surface temperature updated?

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

    real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
    real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation [0-1]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [W/m/K]
    real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]

    ! works
    logical  :: continue_iteration ! does iteration continue?

    real(RP) :: res    ! residual
    real(RP) :: dres   ! d(residual)/dLST
    real(RP) :: oldres ! residual in previous step
    real(RP) :: redf   ! reduced factor

    real(RP) :: Ustar, dUstar ! friction velocity [m]
    real(RP) :: Tstar, dTstar ! friction temperature [K]
    real(RP) :: Qstar, dQstar ! friction mixing rate [kg/kg]
    real(RP) :: Uabs, dUabs ! modified absolute velocity [m/s]
    real(RP) :: SQV, dSQV ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dSHFLX, dLHFLX, dGHFLX

    integer :: i, j, n
    !---------------------------------------------------------------------------

    continue_iteration = LST_UPDATE

    do j = 1, JA
    do i = 1, IA

      if( is_FLX(i,j) ) then
        ! calculate surface flux
        redf   = 1.0_RP
        oldres = 1.0E+5_RP

        ! modified Newton-Raphson method (Tomita 2009)
        do n = 1, nmax
          ! saturation at the surface
          call qsat( SQV, LST(i,j), PRSS(i,j) )

          call CPL_bulkflux( &
              Ustar,     & ! (out)
              Tstar,     & ! (out)
              Qstar,     & ! (out)
              Uabs,      & ! (out)
              TMPA(i,j), & ! (in)
              LST (i,j), & ! (in)
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

          SHFLX(i,j) = -CPdry * RHOA(i,j) * Ustar * Tstar
          LHFLX(i,j) = -LH0   * RHOA(i,j) * Ustar * Qstar * QVEF(i,j)
          GHFLX(i,j) = -2.0_RP * TCS(i,j) * ( LST(i,j) - TG(i,j)  ) / DZG(i,j)

          ! calculation for residual
          res = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
              + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST(i,j)**4 ) &
              - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

          ! d(saturation) at the surface
          call qsat( dSQV, LST(i,j)+dTS, PRSS(i,j) )

          call CPL_bulkflux( &
              dUstar,          & ! (out)
              dTstar,          & ! (out)
              dQstar,          & ! (out)
              dUabs,           & ! (out)
              TMPA(i,j),       & ! (in)
              LST (i,j) + dTS, & ! (in)
              PRSA(i,j),       & ! (in)
              PRSS(i,j),       & ! (in)
              QVA (i,j),       & ! (in)
              dSQV,            & ! (in)
              UA  (i,j),       & ! (in)
              VA  (i,j),       & ! (in)
              Z1  (i,j),       & ! (in)
              PBL (i,j),       & ! (in)
              Z0M (i,j),       & ! (in)
              Z0H (i,j),       & ! (in)
              Z0E (i,j)        ) ! (in)

          dSHFLX  = -CPdry  * RHOA(i,j) * ( (dUstar-Ustar)/dTS * Tstar + Ustar * (dTstar-Tstar)/dTS )
          dLHFLX  = -LH0    * RHOA(i,j) * ( (dUstar-Ustar)/dTS * Qstar + Ustar * (dQstar-Qstar)/dTS ) * QVEF(i,j)
          dGHFLX  = -2.0_RP * TCS(i,j) / DZG(i,j)

          ! calculation for d(residual)/dLST
          dres = -4.0_RP * ( 1.0_RP - ALB_LW(i,j) ) * STB * LST(i,j)**3 &
               - dSHFLX - dLHFLX + dGHFLX

          if( continue_iteration ) then
            if( ( abs(dres) * dres_lim ) < abs(res) ) then
              ! stop iteration to prevent numerical error
              exit
            end if

            if( redf < 0.0_RP ) then
              redf = 1.0_RP
            end if

            if( abs(res) > abs(oldres) ) then
              redf = max( TFa*redf, redf_min )
            else
              redf = min( TFb*redf, redf_max )
            end if

            if( dres > 0.0_RP ) then
              redf = -1.0_RP
            end if

            ! update surface temperature
            LST(i,j) = LST(i,j) - redf * res / dres

            ! put residual in ground heat flux
            GHFLX(i,j) = GHFLX(i,j) - res

            ! save residual in this step
            oldres = res

            if( abs(res) < res_min ) then
              ! iteration converged
              continue_iteration = .false.
            end if

          else
            ! get surface flux without LST updating
            exit
          end if

        end do

        ! diagnositc variables
        U10(i,j) = UA (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        V10(i,j) = VA (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        T2 (i,j) = LST(i,j) + ( TMPA(i,j) - LST(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                       / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
        Q2 (i,j) = SQV      + (  QVA(i,j) - SQV      ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                       / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

        if( n > nmax ) then
          ! not converged and stop program
          if( IO_L ) write(IO_FID_LOG,*) 'Warning: surface tempearture is not converged.'
          if( IO_L ) write(IO_FID_LOG,*) 'Residual [W/m2]', res
        end if

      else
        ! not calculate surface flux
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        ZMFLX(i,j) = 0.0_RP
        SHFLX(i,j) = 0.0_RP
        LHFLX(i,j) = 0.0_RP
        GHFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP
      end if

    enddo
    enddo

    return
  end subroutine CPL_AtmLnd_bulk

end module scale_cpl_atmos_land_bulk
