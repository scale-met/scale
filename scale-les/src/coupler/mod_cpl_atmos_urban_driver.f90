!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Urban Driver
!!
!! @par Description
!!          Coupler driver: atmosphere-urban
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_urban_driver
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
  public :: CPL_AtmUrb_driver_setup
  public :: CPL_AtmUrb_driver

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

  subroutine CPL_AtmUrb_driver_setup
    use mod_atmos_driver, only: &
       ATMOS_SURFACE_SET
    use mod_urban_phy_ucm, only: &
       URBAN_PHY_driver_final
    use scale_cpl_atmos_urban, only: &
       CPL_AtmUrb_setup
    use mod_cpl_vars, only: &
       CPL_TYPE_AtmUrb
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_SURFACE_SET
    call URBAN_PHY_driver_final

    call CPL_AtmUrb_setup( CPL_TYPE_AtmUrb )
    call CPL_AtmUrb_driver( .false. )

    return
  end subroutine CPL_AtmUrb_driver_setup

  subroutine CPL_AtmUrb_driver( update_flag )
    use scale_const, only: &
       LH0  => CONST_LH0
    use scale_grid_real, only: &
       Z1  => REAL_Z1,  &
       LON => REAL_lon, &
       LAT => REAL_lat
    use scale_cpl_atmos_urban, only: &
       CPL_AtmUrb
    use mod_urban_vars, only: &
       TR_URB,  &
       TG_URB,  &
       TB_URB,  &
       TC_URB,  &
       QC_URB,  &
       UC_URB,  &
       TRL_URB, &
       TGL_URB, &
       TBL_URB
    use mod_cpl_vars, only: &
       UST,               &
       DENS => CPL_DENS,  &
       MOMX => CPL_MOMX,  &
       MOMY => CPL_MOMY,  &
       MOMZ => CPL_MOMZ,  &
       TMPA => CPL_TMPA,  &
       QV   => CPL_QV  ,  &
       PREC => CPL_PREC,  &
       SWD  => CPL_SWD ,  &
       LWD  => CPL_LWD ,  &
       CPL_AtmLnd_XMFLX,  & ! tentative
       CPL_AtmLnd_YMFLX,  & ! tentative
       CPL_AtmLnd_ZMFLX,  & ! tentative
       CPL_AtmUrb_XMFLX,  &
       CPL_AtmUrb_YMFLX,  &
       CPL_AtmUrb_ZMFLX,  &
       CPL_AtmUrb_SWUFLX, &
       CPL_AtmUrb_LWUFLX, &
       CPL_AtmUrb_SHFLX,  &
       CPL_AtmUrb_LHFLX,  &
       CPL_AtmUrb_QVFLX,  &
       Urb_GHFLX,         &
       Urb_PRECFLX,       &
       Urb_QVFLX,         &
       CNT_Atm_Urb,       &
       CNT_Urb
    implicit none

    ! argument
    logical, intent(in) :: update_flag

    ! work
    real(RP) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP) :: GHFLX (IA,JA) ! ground heat flux at the surface [W/m2]

    real(RP) :: tmpX(IA,JA) ! temporary XMFLX [kg/m2/s]
    real(RP) :: tmpY(IA,JA) ! temporary YMFLX [kg/m2/s]

    logical  :: LSOLAR = .false.    ! logical [true=both, false=SSG only]
    real(RP) :: QA       ! mixing ratio at 1st atmospheric level  [kg/kg]
    real(RP) :: UA       ! wind speed at 1st atmospheric level    [m/s]
    real(RP) :: U1       ! u at 1st atmospheric level             [m/s]
    real(RP) :: V1       ! v at 1st atmospheric level             [m/s]

    ! parameters for shadow model
    !XX sinDEC, cosDEC ?
    ! real(RP) :: DECLIN   ! solar declination                    [rad]
    !XX  cosSZA     ?
    ! real(RP) :: COSZ     ! sin(fai)*sin(del)+cos(fai)*cos(del)*cos(omg)
    !XX  hourangle  ?
    ! real(RP) :: OMG      ! solar hour angle                       [rad]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Urban'

    do j = JS-1, JE+1
    do i = IS-1, IE+1

      QA = QV(i,j) / (1.0_RP-QV(i,j)) ! mixing ratio at 1st atmospheric level [kg/kg]
                                      ! QV specific humidity                  [kg/kg]
      UA = sqrt(          &
            ( MOMZ(i,j)               )**2 &
          + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
          + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
          ) / DENS(i,j) * 0.5_RP
                                  ! wind speed at 1st atmospheric level   [m/s]
      U1 = 0.5_RP * ( MOMX(i-1,j) + MOMX(i,j) ) / DENS(i,j)
                                  ! u at 1st atmospheric level            [m/s]
      V1 = 0.5_RP * ( MOMY(i,j-1) + MOMY(i,j) ) / DENS(i,j)
                                  ! v at 1st atmospheric level            [m/s]

      call CPL_AtmUrb( &
        TR_URB (i,j),   & ! (inout)
        TB_URB (i,j),   & ! (inout)
        TG_URB (i,j),   & ! (inout)
        TC_URB (i,j),   & ! (inout)
        QC_URB (i,j),   & ! (inout)
        UC_URB (i,j),   & ! (inout)
        TRL_URB(:,i,j), & ! (inout)
        TBL_URB(:,i,j), & ! (inout)
        TGL_URB(:,i,j), & ! (inout)
        UST    (i,j),   & ! (out)
        SHFLX  (i,j),   & ! (out)
        LHFLX  (i,j),   & ! (out)
        SWUFLX (i,j),   & ! (out)
        LWUFLX (i,j),   & ! (out)
        GHFLX  (i,j),   & ! (out)
        LSOLAR,         & ! (in)
        TMPA   (i,j),   & ! (in)
        QA,             & ! (in)
        UA,             & ! (in)
        U1,             & ! (in)
        V1,             & ! (in)
        Z1     (i,j),   & ! (in)
        SWD    (i,j),   & ! (in)
        LWD    (i,j),   & ! (in)
        PREC   (i,j),   & ! (in)
        DENS   (i,j),   & ! (in)
        LON    (i,j),   & ! (in)
        LAT    (i,j)    ) ! (in)

    end do
    end do

    ! --- tentative: caution !! ---
    XMFLX(:,:) = CPL_AtmLnd_XMFLX (:,:)
    YMFLX(:,:) = CPL_AtmLnd_YMFLX (:,:)
    ZMFLX(:,:) = CPL_AtmLnd_ZMFLX (:,:)
    ! --- tentative: caution !! ---

    ! interpolate momentum fluxes
    do j = JS, JE
    do i = IS, IE
      tmpX(i,j) = ( XMFLX(i,j) + XMFLX(i+1,j  ) ) * 0.5_RP ! at u/y-layer
      tmpY(i,j) = ( YMFLX(i,j) + YMFLX(i,  j+1) ) * 0.5_RP ! at x/v-layer
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
      XMFLX(i,j) = tmpX(i,j)
      YMFLX(i,j) = tmpY(i,j)
      ZMFLX(i,j) = ZMFLX(i,j) * 0.5_RP ! at w-layer
    enddo
    enddo

    ! temporal average flux
    CPL_AtmUrb_XMFLX (:,:) = ( CPL_AtmUrb_XMFLX (:,:) * CNT_Atm_Urb + XMFLX (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_YMFLX (:,:) = ( CPL_AtmUrb_YMFLX (:,:) * CNT_Atm_Urb + YMFLX (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_ZMFLX (:,:) = ( CPL_AtmUrb_ZMFLX (:,:) * CNT_Atm_Urb + ZMFLX (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_SWUFLX(:,:) = ( CPL_AtmUrb_SWUFLX(:,:) * CNT_Atm_Urb + SWUFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_LWUFLX(:,:) = ( CPL_AtmUrb_LWUFLX(:,:) * CNT_Atm_Urb + LWUFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_SHFLX (:,:) = ( CPL_AtmUrb_SHFLX (:,:) * CNT_Atm_Urb + SHFLX (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_LHFLX (:,:) = ( CPL_AtmUrb_LHFLX (:,:) * CNT_Atm_Urb + LHFLX (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_QVFLX (:,:) = ( CPL_AtmUrb_QVFLX (:,:) * CNT_Atm_Urb + LHFLX (:,:)/LH0 ) / ( CNT_Atm_Urb + 1.0_RP )

    Urb_GHFLX  (:,:) = ( Urb_GHFLX  (:,:) * CNT_Urb + GHFLX(:,:)     ) / ( CNT_Urb + 1.0_RP )
    Urb_PRECFLX(:,:) = ( Urb_PRECFLX(:,:) * CNT_Urb + PREC (:,:)     ) / ( CNT_Urb + 1.0_RP )
    Urb_QVFLX  (:,:) = ( Urb_QVFLX  (:,:) * CNT_Urb - LHFLX(:,:)/LH0 ) / ( CNT_Urb + 1.0_RP )

    CNT_Atm_Urb = CNT_Atm_Urb + 1.0_RP
    CNT_Urb     = CNT_Urb     + 1.0_RP

    return
  end subroutine CPL_AtmUrb_driver

end module mod_cpl_atmos_urban_driver
