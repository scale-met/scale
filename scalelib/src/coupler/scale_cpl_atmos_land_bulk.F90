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
  integer,  private, save :: nmax = 100 ! maximum iteration number

  real(RP), private, save :: res_min  =   1.0_RP    ! minimum number of residual
  real(RP), private, save :: dTS      =   1.0E-8_RP ! delta surface temp.
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
  subroutine CPL_AtmLnd_bulk_setup( CPL_TYPE_AtmLnd )
    use scale_process, only: &
       PRC_MPIstop
    use scale_cpl_bulkcoef, only: &
       CPL_bulkcoef_setup
    implicit none

    character(len=H_SHORT), intent(in) :: CPL_TYPE_AtmLnd

    integer  :: CPL_AtmLnd_bulk_nmax
    real(RP) :: CPL_AtmLnd_bulk_res_min
    real(RP) :: CPL_AtmLnd_bulk_dTS
    real(RP) :: CPL_AtmLnd_bulk_U_minM
    real(RP) :: CPL_AtmLnd_bulk_U_minH
    real(RP) :: CPL_AtmLnd_bulk_U_minE
    real(RP) :: CPL_AtmLnd_bulk_U_maxM
    real(RP) :: CPL_AtmLnd_bulk_U_maxH
    real(RP) :: CPL_AtmLnd_bulk_U_maxE

    NAMELIST / PARAM_CPL_ATMLND_BULK / &
       CPL_AtmLnd_bulk_nmax,    &
       CPL_AtmLnd_bulk_res_min, &
       CPL_AtmLnd_bulk_dTS,     &
       CPL_AtmLnd_bulk_U_minM,  &
       CPL_AtmLnd_bulk_U_minH,  &
       CPL_AtmLnd_bulk_U_minE,  &
       CPL_AtmLnd_bulk_U_maxM,  &
       CPL_AtmLnd_bulk_U_maxH,  &
       CPL_AtmLnd_bulk_U_maxE

    integer :: ierr
    !---------------------------------------------------------------------------

    CPL_AtmLnd_bulk_nmax    = nmax
    CPL_AtmLnd_bulk_res_min = res_min
    CPL_AtmLnd_bulk_dTS     = dTS
    CPL_AtmLnd_bulk_U_minM  = U_minM
    CPL_AtmLnd_bulk_U_minH  = U_minH
    CPL_AtmLnd_bulk_U_minE  = U_minE
    CPL_AtmLnd_bulk_U_maxM  = U_maxM
    CPL_AtmLnd_bulk_U_maxH  = U_maxH
    CPL_AtmLnd_bulk_U_maxE  = U_maxE

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land: bulk flux parameter'

    if ( CPL_TYPE_AtmLnd /= 'BULK' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx CPL_TYPE_AtmLnd is not BULK. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_ATMLND_BULK,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_ATMLND_BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL_ATMLND_BULK)

    nmax    = CPL_AtmLnd_bulk_nmax
    res_min = CPL_AtmLnd_bulk_res_min
    dTS     = CPL_AtmLnd_bulk_dTS
    U_minM  = CPL_AtmLnd_bulk_U_minM
    U_minH  = CPL_AtmLnd_bulk_U_minH
    U_minE  = CPL_AtmLnd_bulk_U_minE
    U_maxM  = CPL_AtmLnd_bulk_U_maxM
    U_maxH  = CPL_AtmLnd_bulk_U_maxH
    U_maxE  = CPL_AtmLnd_bulk_U_maxE

    !--- set up bulk coefficient function
    call CPL_bulkcoef_setup

    return
  end subroutine CPL_AtmLnd_bulk_setup

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
    real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

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
    real(RP) :: RES   (IA,JA)
    real(RP) :: DRES  (IA,JA)
    real(RP) :: oldRES(IA,JA) ! RES in previous step
    real(RP) :: redf  (IA,JA) ! reduced factor

    real(RP) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP) :: Cm, Ch, Ce, dCm, dCh, dCe ! bulk transfer coeff. [no unit]
    real(RP) :: R10m, R02h, R02e, dR10m, dR02h, dR02e ! lapse rate [0-1]
    real(RP) :: SQV, dSQV ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dSHFLX, dLHFLX, dGHFLX

    integer :: i, j, n
    !---------------------------------------------------------------------------

    redf  (:,:) = 1.0_RP
    oldRES(:,:) = 1.0E+5_RP

    do n = 1, nmax

      ! calculate surface flux
      do j = JS, JE
      do i = IS, IE
        Uabs = sqrt( UA(i,j)**2 + VA(i,j)**2 + WA(i,j)**2 )

        call CPL_bulkcoef( &
            Cm,        & ! (out)
            Ch,        & ! (out)
            Ce,        & ! (out)
            R10m,      & ! (out)
            R02h,      & ! (out)
            R02e,      & ! (out)
            TMPA(i,j), & ! (in)
            LST (i,j), & ! (in)
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
        call qsat( SQV, LST(i,j), PRSS(i,j) )

        SHFLX (i,j) = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOA(i,j) * Ch * ( LST(i,j) - TMPA(i,j) )
        LHFLX (i,j) = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOA(i,j) * QVEF(i,j) * Ce * ( SQV - QVA(i,j) )
        GHFLX (i,j) = -2.0_RP * TCS(i,j) * ( LST(i,j) - TG(i,j)  ) / DZG(i,j)

        ! calculation for residual
        RES(i,j) = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
                 + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST(i,j)**4 ) &
                 - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

        call CPL_bulkcoef( &
            dCm,             & ! (out)
            dCh,             & ! (out)
            dCe,             & ! (out)
            dR10m,           & ! (out)
            dR02h,           & ! (out)
            dR02e,           & ! (out)
            TMPA(i,j),       & ! (in)
            LST (i,j) + dTS, & ! (in)
            PRSA(i,j),       & ! (in)
            PRSS(i,j),       & ! (in)
            Uabs,            & ! (in)
            Z1  (i,j),       & ! (in)
            Z0M (i,j),       & ! (in)
            Z0H (i,j),       & ! (in)
            Z0E (i,j)        ) ! (in)

        call qsat( dSQV, LST(i,j)+dTS, PRSS(i,j) )

        dSHFLX  = CPdry * min(max(Uabs,U_minH),U_maxH) * RHOA(i,j) &
                * ( (dCh-Ch)/dTS * ( LST(i,j) - TMPA(i,j) ) + Ch )
        dLHFLX  = LH0   * min(max(Uabs,U_minE),U_maxE) * RHOA(i,j) * QVEF(i,j) &
                * ( (dCe-Ce)/dTS * ( SQV - QVA(i,j) ) + Ce * (dSQV-SQV)/dTS )
        dGHFLX  = -2.0_RP * TCS(i,j) / DZG(i,j)

        ! calculation for d(residual)/dTS
        DRES(i,j) = -4.0_RP * ( 1.0_RP - ALB_LW(i,j) ) * STB * LST(i,j)**3 &
                  - dSHFLX - dLHFLX + dGHFLX

        ! diagnositc variables
        U10(i,j) = R10m * UA(i,j)
        V10(i,j) = R10m * VA(i,j)

        T2(i,j) = (          R02h ) * TMPA(i,j) &
                + ( 1.0_RP - R02h ) * LST (i,j)
        Q2(i,j) = (          R02e ) * QVA (i,j) &
                + ( 1.0_RP - R02e ) * QVEF(i,j) * SQV

      enddo
      enddo

      if( LST_UPDATE ) then

        do j = JS, JE
        do i = IS, IE

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
          LST(i,j)  = LST(i,j) - redf(i,j) * RES(i,j)/DRES(i,j)

          ! put residual in ground heat flux
          GHFLX(i,j) = GHFLX(i,j) - RES(i,j)

          ! save residual in this step
          oldRES(i,j) = RES(i,j)

        end do
        end do

        if( maxval(abs(RES(IS:IE,JS:JE))) < res_min ) then
          ! iteration converged
          exit
        end if

      else
        ! get surface flux without LST updating
        exit

      end if

    end do

    if( n > nmax ) then
      ! not converged and stop program
      if( IO_L ) write(IO_FID_LOG,*) 'Error: surface tempearture is not converged.'
      call PRC_MPIstop
    end if

    return
  end subroutine CPL_AtmLnd_bulk

end module scale_cpl_atmos_land_bulk
