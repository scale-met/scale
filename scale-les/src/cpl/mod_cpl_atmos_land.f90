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
module mod_cpl_atmos_land
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
  public :: CPL_AtmLnd_driver_setup
  public :: CPL_AtmLnd_driver

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
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmLnd_driver_setup
    use mod_cpl_vars, only: &
       CPL_flushAtm,            &
       CPL_flushLnd,            &
       CPL_AtmLnd_flushCPL
    implicit none
    !---------------------------------------------------------------------------

    call CPL_flushAtm
    call CPL_flushLnd
    call CPL_AtmLnd_flushCPL

    return
  end subroutine CPL_AtmLnd_driver_setup

  subroutine CPL_AtmLnd_driver( update_flag )
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    use mod_cpl_vars, only: &
       CPL_AtmLnd_putCPL,     &
       CPL_AtmLnd_getAtm2CPL, &
       CPL_AtmLnd_getLnd2CPL, &
       LST
    use mod_cpl_sfcflux, only: &
       CPL_sfcflux
    implicit none

    ! parameters
    integer,  parameter :: nmax     = 100       ! maximum iteration number
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)
    real(RP), parameter :: res_min  = 1.0_RP    ! minimum number of residual

    ! argument
    logical, intent(in) :: update_flag

    ! works
    integer :: i, j, n

    real(RP) :: RES   (IA,JA)
    real(RP) :: DRES  (IA,JA)
    real(RP) :: oldRES(IA,JA) ! RES in previous step
    real(RP) :: redf  (IA,JA) ! reduced factor

    real(RP) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP) :: GHFLX (IA,JA) ! ground heat flux at the surface [W/m2]

    real(RP) :: DZ  (IA,JA) ! height from the surface to the lowest atmospheric layer [m]

    real(RP) :: DENS(IA,JA) ! air density at the lowest atmospheric layer [kg/m3]
    real(RP) :: MOMX(IA,JA) ! momentum x at the lowest atmospheric layer [kg/m2/s]
    real(RP) :: MOMY(IA,JA) ! momentum y at the lowest atmospheric layer [kg/m2/s]
    real(RP) :: MOMZ(IA,JA) ! momentum z at the lowest atmospheric layer [kg/m2/s]
    real(RP) :: RHOS(IA,JA) ! air density at the sruface [kg/m3]
    real(RP) :: PRES(IA,JA) ! pressure at the surface [Pa]
    real(RP) :: ATMP(IA,JA) ! air temperature at the surface [K]
    real(RP) :: QV  (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP) :: PREC(IA,JA) ! precipitaton flux at the surface [kg/m2/s]
    real(RP) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface (upward positive) [W/m2]
    real(RP) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface (upward positive) [W/m2]

    real(RP) :: TG  (IA,JA) ! soil temperature [K]
    real(RP) :: QVEF(IA,JA) ! efficiency of evaporation [no unit]
    real(RP) :: EMIT(IA,JA) ! emissivity in long-wave radiation [no unit]
    real(RP) :: ALB (IA,JA) ! surface albedo in short-wave radiation [no unit]
    real(RP) :: TCS (IA,JA) ! thermal conductivity for soil [W/m/K]
    real(RP) :: DZG (IA,JA) ! soil depth [m]
    real(RP) :: Z0M (IA,JA) ! roughness length for momemtum [m]
    real(RP) :: Z0H (IA,JA) ! roughness length for heat [m]
    real(RP) :: Z0E (IA,JA) ! roughness length for vapor [m]
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Land'

    call CPL_AtmLnd_getAtm2CPL( &
      DENS, MOMX, MOMY, MOMZ, & ! (out)
      RHOS, PRES, ATMP, QV,   & ! (out)
      PREC, SWD, LWD          ) ! (out)

    call CPL_AtmLnd_getLnd2CPL( &
      TG, QVEF, EMIT, & ! (out)
      ALB, TCS, DZG,  & ! (out)
      Z0M, Z0H, Z0E   ) ! (out)

    DZ(:,:) = CZ(KS,:,:) - FZ(KS-1,:,:)

    redf  (:,:) = 1.0_RP
    oldRES(:,:) = 1.0E+5_RP

    do n = 1, nmax
      ! calculate surface flux
      call CPL_sfcflux( &
        RES, DRES,                                              & ! (out)
        XMFLX, YMFLX, ZMFLX,                                    & ! (out)
        SWUFLX, LWUFLX, SHFLX, LHFLX, GHFLX,                    & ! (out)
        DZ, LST,                                                & ! (in)
        DENS, MOMX, MOMY, MOMZ, RHOS, PRES, ATMP, QV, SWD, LWD, & ! (in)
        TG, QVEF, EMIT, ALB, TCS, DZG, Z0M, Z0H, Z0E            ) ! (in)

      if( update_flag ) then

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
          LST(i,j)  = LST(i,j) - redf(i,j) * RES(i,j)/DRES(i,j)

          ! put residual in ground heat flux
          GHFLX(i,j) = GHFLX(i,j) - RES(i,j)

          ! save residual in this step
          oldRES(i,j) = RES(i,j)

        end do
        end do

        if( maxval(abs(RES(IS-1:IE+1,JS-1:JE+1))) < res_min ) then
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

    call CPL_AtmLnd_putCPL( &
      XMFLX, YMFLX, ZMFLX,  &
      SWUFLX, LWUFLX,       &
      SHFLX, LHFLX, GHFLX,  &
      PREC                  )

    return
  end subroutine CPL_AtmLnd_driver

end module mod_cpl_atmos_land
