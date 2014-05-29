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
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_driver_setup
    use mod_cpl_admin, only: &
       CPL_sw_AtmUrb,  &
       CPL_TYPE_AtmUrb
    use scale_cpl_atmos_urban, only: &
       CPL_AtmUrb_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[CPL AtmUrb] / Origin[SCALE-LES]'

    if ( CPL_sw_AtmUrb ) then

       call CPL_AtmUrb_setup( CPL_TYPE_AtmUrb )

       call CPL_AtmUrb_driver( .false. )

    endif

    return
  end subroutine CPL_AtmUrb_driver_setup

  !-----------------------------------------------------------------------------
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
       TS_URB,  &
       TRL_URB, &
       TGL_URB, &
       TBL_URB, &
       SHR_URB, &
       SHB_URB, &
       SHG_URB, &
       LHR_URB, &
       LHB_URB, &
       LHG_URB, &
       GHR_URB, &
       GHB_URB, &
       GHG_URB, &
       RnR_URB, &
       RnB_URB, &
       RnG_URB
    use mod_cpl_vars, only: &
       UST,               &
       RHOA => CPL_RHOA,  &
       UA   => CPL_UA,    &
       VA   => CPL_VA,    &
       WA   => CPL_WA,    &
       TMPA => CPL_TMPA,  &
       QVA  => CPL_QVA,   &
       PREC => CPL_PREC,  &
       SWD  => CPL_SWD,   &
       LWD  => CPL_LWD,   &
       CPL_AtmLnd_XMFLX,  & ! tentative
       CPL_AtmLnd_YMFLX,  & ! tentative
       CPL_AtmLnd_ZMFLX,  & ! tentative
       CPL_AtmLnd_U10,    & ! tentative
       CPL_AtmLnd_V10,    & ! tentative
       CPL_AtmLnd_T2,     & ! tentative
       CPL_AtmLnd_Q2,     & ! tentative
       CPL_AtmUrb_XMFLX,  &
       CPL_AtmUrb_YMFLX,  &
       CPL_AtmUrb_ZMFLX,  &
       CPL_AtmUrb_SHFLX,  &
       CPL_AtmUrb_LHFLX,  &
       CPL_AtmUrb_QVFLX,  &
       CPL_AtmUrb_U10,    &
       CPL_AtmUrb_V10,    &
       CPL_AtmUrb_T2,     &
       CPL_AtmUrb_Q2,     &
       Urb_GHFLX,         &
       Urb_PRECFLX,       &
       Urb_QVFLX,         &
       CNT_Atm_Urb,       &
       CNT_Urb
    implicit none

    ! arguments
    logical, intent(in) :: update_flag

    ! works
    real(RP) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP) :: SHFLX(IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP) :: LHFLX(IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP) :: GHFLX(IA,JA) ! ground heat flux at the surface [W/m2]
    real(RP) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

    logical  :: LSOLAR = .false.    ! logical [true=both, false=SSG only]
    real(RP) :: QMA                 ! mixing ratio at the lowest atmospheric level  [kg/kg]
    real(RP) :: Uabs                ! wind speed at the lowest atmospheric level    [m/s]

    ! parameters for shadow model
    !XX sinDEC, cosDEC ?
    ! real(RP) :: DECLIN   ! solar declination                    [rad]
    !XX  cosSZA     ?
    ! real(RP) :: COSZ     ! sin(fai)*sin(del)+cos(fai)*cos(del)*cos(omg)
    !XX  hourangle  ?
    ! real(RP) :: OMG      ! solar hour angle                       [rad]

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Urban'

    do j = 1, JA
    do i = 1, IA

      QMA = QVA(i,j) / (1.0_RP-QVA(i,j)) ! mixing ratio at 1st atmospheric level [kg/kg]
                                         ! QV specific humidity                  [kg/kg]
      Uabs = sqrt( UA(i,j)**2 + VA(i,j)**2 + WA(i,j)**2 )

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
        TS_URB  (i,j),  & ! (out)
        SHR_URB (i,j),  & ! (out)
        SHB_URB (i,j),  & ! (out)
        SHG_URB (i,j),  & ! (out)
        LHR_URB (i,j),  & ! (out)
        LHB_URB (i,j),  & ! (out)
        LHG_URB (i,j),  & ! (out)
        GHR_URB (i,j),  & ! (out)
        GHB_URB (i,j),  & ! (out)
        GHG_URB (i,j),  & ! (out)
        RnR_URB (i,j),  & ! (out)
        RnB_URB (i,j),  & ! (out)
        RnG_URB (i,j),  & ! (out)
        UST    (i,j),   & ! (out)
        SHFLX  (i,j),   & ! (out)
        LHFLX  (i,j),   & ! (out)
        GHFLX  (i,j),   & ! (out)
        LSOLAR,         & ! (in)
        TMPA   (i,j),   & ! (in)
        QMA,            & ! (in)
        Uabs,           & ! (in)
        UA     (i,j),   & ! (in)
        VA     (i,j),   & ! (in)
        Z1     (i,j),   & ! (in)
        SWD    (i,j),   & ! (in)
        LWD    (i,j),   & ! (in)
        PREC   (i,j),   & ! (in)
        RHOA   (i,j),   & ! (in)
        LON    (i,j),   & ! (in)
        LAT    (i,j)    ) ! (in)

    end do
    end do

    ! --- tentative: caution !! ---
    XMFLX(:,:) = CPL_AtmLnd_XMFLX(:,:)
    YMFLX(:,:) = CPL_AtmLnd_YMFLX(:,:)
    ZMFLX(:,:) = CPL_AtmLnd_ZMFLX(:,:)
    U10  (:,:) = CPL_AtmLnd_U10  (:,:)
    V10  (:,:) = CPL_AtmLnd_V10  (:,:)
    T2   (:,:) = CPL_AtmLnd_T2   (:,:)
    Q2   (:,:) = CPL_AtmLnd_Q2   (:,:)
    ! --- tentative: caution !! ---

    ! temporal average flux
    CPL_AtmUrb_XMFLX(:,:) = ( CPL_AtmUrb_XMFLX(:,:) * CNT_Atm_Urb + XMFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_YMFLX(:,:) = ( CPL_AtmUrb_YMFLX(:,:) * CNT_Atm_Urb + YMFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_ZMFLX(:,:) = ( CPL_AtmUrb_ZMFLX(:,:) * CNT_Atm_Urb + ZMFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_SHFLX(:,:) = ( CPL_AtmUrb_SHFLX(:,:) * CNT_Atm_Urb + SHFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_LHFLX(:,:) = ( CPL_AtmUrb_LHFLX(:,:) * CNT_Atm_Urb + LHFLX(:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_QVFLX(:,:) = ( CPL_AtmUrb_QVFLX(:,:) * CNT_Atm_Urb + LHFLX(:,:)/LH0 ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_U10  (:,:) = ( CPL_AtmUrb_U10  (:,:) * CNT_Atm_Urb + U10  (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_V10  (:,:) = ( CPL_AtmUrb_V10  (:,:) * CNT_Atm_Urb + V10  (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_T2   (:,:) = ( CPL_AtmUrb_T2   (:,:) * CNT_Atm_Urb + T2   (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )
    CPL_AtmUrb_Q2   (:,:) = ( CPL_AtmUrb_Q2   (:,:) * CNT_Atm_Urb + Q2   (:,:)     ) / ( CNT_Atm_Urb + 1.0_RP )

    Urb_GHFLX  (:,:) = ( Urb_GHFLX  (:,:) * CNT_Urb + GHFLX(:,:)     ) / ( CNT_Urb + 1.0_RP )
    Urb_PRECFLX(:,:) = ( Urb_PRECFLX(:,:) * CNT_Urb + PREC (:,:)     ) / ( CNT_Urb + 1.0_RP )
    Urb_QVFLX  (:,:) = ( Urb_QVFLX  (:,:) * CNT_Urb - LHFLX(:,:)/LH0 ) / ( CNT_Urb + 1.0_RP )

    CNT_Atm_Urb = CNT_Atm_Urb + 1.0_RP
    CNT_Urb     = CNT_Urb     + 1.0_RP

    return
  end subroutine CPL_AtmUrb_driver

end module mod_cpl_atmos_urban_driver
