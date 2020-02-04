!-------------------------------------------------------------------------------
!> module land / physics / snow / ky90
!!
!! @par Description
!!       Profile model for snowpack by Kondo and Yamazaki (1990)
!!       Kondo and Yamazaki (1990): Journal of Applied Meteorology, 29, pp375-384.
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_phy_snow_ky90
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !--------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_SNOW_KY90_setup
  public :: LAND_PHY_SNOW_KY90

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public          :: W0                    ! Maximum water content [ratio:0-1]
  real(RP), public          :: RHOSNOW   = 400.0_RP  ! Snow density [kg/m3]
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !
  real(RP)                  :: A0, A1, A2, C1, C2, C3
  real(RP)                  :: ESAT, QSAT, DELTAQSAT, CP
  real(RP)                  :: Gres

  integer, parameter        :: data_length_max=10000

  ! model parameters
  real(RP)                  :: LAMBDAS               ! Snow thermal conductivity [W/m/K]
  real(RP)                  :: CSRHOS                ! Heat capacity             [J/m3/K]
                                                     ! (Specific heat of snow [J/kg/K])*(Snow density [kg/m3])
  !real(RP)                  :: RHOSNOW   = 400.0_RP  ! Snow density [kg/m3]
  real(RP)                  :: CS                    ! Specific heat for unit mass [J/kg/K]
  real(RP)                  :: ALBEDO                ! albedo
  real(RP)                  :: ALBEDOMIN = 0.4_RP
  real(RP)                  :: ALBEDOMAX = 0.85_RP
  real(RP)                  :: ZMIN      = 0.01_RP   ! Minimum snow depth [m]

  ! constant variables
  real(RP), parameter       :: epsilon   = 0.97_RP
  real(RP), parameter       :: sigma     = 5.67e-08_RP
  !real(RP), parameter       :: rhoair    = 1.289_RP      ! [kg/m3]
  real(RP), parameter       :: CH        = 0.002_RP
  real(RP), parameter       :: CE        = 0.0021_RP
  real(RP), parameter       :: LV        = 2.5e6_RP       ! [J/kg] at 0 deg.C
  real(RP), parameter       :: LF        = 3.34e5_RP      ! [J/kg] at 0 deg.C

  logical                   :: ALBEDO_const = .true.
  logical                   :: debug = .false.

  integer                   :: ZN_flag, TS_flag, sflag

  !----------------------------------------------------------------------------------------------!
contains
  !----------------------------------------------------------------------------------------------!
  !> Setup
  subroutine LAND_PHY_SNOW_KY90_setup
    use scale_prc, only: &
       PRC_abort
    implicit none
    real(RP) :: snow_conductivity     = 0.42_RP
    real(RP) :: water_content         = 0.1_RP
    real(RP) :: snow_heat_capacityRHO = 8.4e+05_RP
    real(RP) :: snow_rho              = 400.0_RP
    real(RP) :: snowDepth_initial     = 0.0_RP
    real(RP) :: albedo_value          = 0.686_RP

    namelist / PARAM_LAND_PHY_SNOW_KY90 / &
       ALBEDO_const,          &
       snow_conductivity,     &
       water_content,         &
       snow_heat_capacityRHO, &
       snow_rho,              &
       snowDepth_initial,     &
       albedo_value,          &
       debug

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_PHY_SNOW_KY90_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PHY_SNOW_KY90,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_PHY_SNOW_KY90_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_PHY_SNOW_KY90_setup",*) 'Not appropriate names in namelist PARAM_LAND_PHY_SNOW_KY90. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND_PHY_SNOW_KY90)

    LAMBDAS       = snow_conductivity
    W0            = water_content
    CSRHOS        = snow_heat_capacityRHO ! [J/m3/K]
    RHOSNOW       = snow_rho
    CS            = CSRHOS/RHOSNOW        ! [J/kg/K]
    !LOG_WARN("LAND_PHY_SNOW_KY90_setup",*)   "Specific heat capacity of snow [J/kg/K]: ",CS
    ALBEDO        = albedo_value

    return
  endsubroutine LAND_PHY_SNOW_KY90_setup

  !-----------------------------------------------------------------------------
  !> Main routine for land submodel

  subroutine LAND_PHY_SNOW_KY90( &
       LIA, LIS, LIE, LJA, LJS, LJE, &
       SFLX_water, SFLX_ENGI,  & ! [IN]
       PRSA, TA, QA,           & ! [IN]
       WA, UA, VA,             & ! [IN]
       DENS,                   & ! [IN]
       SFLX_RAD_dn,            & ! [IN]
       exists_land, dt,        & ! [IN]
       TSNOW, SWE,             & ! [INOUT]
       SDepth, SDzero,         & ! [INOUT]
       nosnowsec,              & ! [INOUT]
       Salbedo,                & ! [OUT]
       SFLX_SH,                & ! [OUT]
       SFLX_LH, SFLX_QV,       & ! [OUT]
       SFLX_QV_ENGI,           & ! [OUT]
       SFLX_GH, SNOW_LAND_GH,  & ! [OUT]
       SNOW_LAND_Water,        & ! [OUT]
       SNOW_frac               ) ! [OUT]
    use scale_prc, only: &
       PRC_abort
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_saturation, only:  &
       qsatf => ATMOS_SATURATION_psat_all  ! better to  change name from qsatf to qsat
    use scale_const, only:   &
       EPS   => CONST_EPS,   &
       T0    => CONST_TEM00, &
       I_SW  => CONST_I_SW,  &
       I_LW  => CONST_I_LW,  &
       EPSvap => CONST_EPSvap
    implicit none
    integer, intent(in) :: LIA, LIS, LIE
    integer, intent(in) :: LJA, LJS, LJE

    ! input data
    real(RP), intent(in)      :: SFLX_water(LIA,LJA)
    real(RP), intent(in)      :: SFLX_ENGI (LIA,LJA)
    real(RP), intent(in)      :: PRSA      (LIA,LJA)
    real(RP), intent(in)      :: TA        (LIA,LJA)
    real(RP), intent(in)      :: WA        (LIA,LJA)
    real(RP), intent(in)      :: UA        (LIA,LJA)
    real(RP), intent(in)      :: VA        (LIA,LJA)
    real(RP), intent(in)      :: QA        (LIA,LJA)      ! specific humidity [kg/kg]
    real(RP), intent(in)      :: DENS      (LIA,LJA)
    real(RP), intent(in)      :: SFLX_RAD_dn(LIA,LJA,N_RAD_DIR,N_RAD_RGN)
    real(DP), intent(in)      :: dt                     ! dt of land
    logical,  intent(in)      :: exists_land(LIA,LJA)

    ! prognostic variables
    real(RP), intent(inout)   :: TSNOW          (LIA,LJA)   ! snow temperature        [K]
    real(RP), intent(inout)   :: SWE            (LIA,LJA)   ! equivalent water        [kg/m2]
    real(RP), intent(inout)   :: SDepth         (LIA,LJA)   ! depth of melting point  [m]
    real(RP), intent(inout)   :: SDzero         (LIA,LJA)   ! total snow depth        [m]
    real(RP), intent(inout)   :: nosnowsec      (LIA,LJA)   ! elapsed time of no snow condition [s]

    ! updated variables
    real(RP), intent(out)     :: Salbedo        (LIA,LJA,2) ! snow albedo             [-]
    real(RP), intent(out)     :: SFLX_SH        (LIA,LJA) ! sensible heat flux between atmos and snow [W/m2]
    real(RP), intent(out)     :: SFLX_LH        (LIA,LJA) ! latente  heat flux between atmos and snow [W/m2]
    real(RP), intent(out)     :: SFLX_QV        (LIA,LJA) ! evaporation due to LH          [kg/m2/s]
    real(RP), intent(out)     :: SFLX_QV_ENGI   (LIA,LJA) ! internal energy of evaporation [J/m2/s]
    real(RP), intent(out)     :: SFLX_GH        (LIA,LJA) ! whole snowpack Ground flux   [W/m2]
    real(RP), intent(out)     :: SNOW_LAND_GH   (LIA,LJA) ! heat flux from snow to land  [W/m2]
    real(RP), intent(out)     :: SNOW_LAND_water(LIA,LJA) ! water flux from snow to land [W/m2]
    real(RP), intent(out)     :: SNOW_frac      (LIA,LJA) ! snow fraction, defined by time direction [-]

    real(RP)                  :: QCC       (LIA,LJA)
    real(RP)                  :: QFUSION   (LIA,LJA)
    real(RP)                  :: MELT      (LIA,LJA)
    real(RP)                  :: SWEMELT   (LIA,LJA)

    ! works
    real(RP)                  :: TSNOW1           ! updated snow surface temperature [K]
    real(RP)                  :: ZNSNOW1          ! updated freezing depth           [m]
    real(RP)                  :: SWE1             ! updated snow water equivalence   [kg/m2]
    real(RP)                  :: DEPTH1           ! updated snow depth               [m]

    real(RP), parameter       :: Uabs_min = 0.1_RP
    real(RP)                  :: Uabs
    real(RP)                  :: QAsat
    real(RP)                  :: RH
    real(RP)                  :: qdry, psat
    real(RP)                  :: SFLX_SW_dn, SFLX_LW_dn

    real(RP)                  :: w

    integer :: k, i, j
    !---------------------------------------------------------------------------
    LOG_PROGRESS(*) 'Land / physics / snow / KY90'

    do j = LJS, LJE
    do i = LIS, LIE

    if( ( exists_land(i,j) ) .and. &
        ( SWE(i,j)>0. .or. SFLX_water(i,j)>0. ) )then

       Uabs = sqrt( WA(i,j)**2 + UA(i,j)**2 + VA(i,j)**2 )

       !qdry = 1.0_RP - QA(i,j)
       !call qsatf(  TA(i,j), PRSA(i,j), qdry, & ![IN]
       !             QAsat                     ) ![OUT]
       call qsatf(  TA(i,j), psat )
       QAsat = EPSvap * psat / ( PRSA(i,j) - ( 1.0_RP-EPSvap ) * psat )
       RH  = QA(i,j) / QAsat
       if ( debug ) then
          LOG_INFO("LAND_PHY_SNOW_KY90",*) "RH,  ",RH," DENS (1.289 org) ",DENS(i,j)
       end if

       SFLX_LW_dn = SFLX_RAD_dn(i,j,I_R_diffuse,I_R_IR )
       SFLX_SW_dn = SFLX_RAD_dn(i,j,I_R_direct ,I_R_NIR) &
                  + SFLX_RAD_dn(i,j,I_R_diffuse,I_R_NIR) &
                  + SFLX_RAD_dn(i,j,I_R_direct ,I_R_VIS) &
                  + SFLX_RAD_dn(i,j,I_R_diffuse,I_R_VIS)

       TSNOW1  = TSNOW (i,j)
       SWE1    = SWE   (i,j)
       DEPTH1  = SDepth(i,j)
       ZNSNOW1 = SDzero(i,j)

       call SNOW_ky90_main( TSNOW1,                 & ! [INOUT]
                            SWE1,                   & ! [INOUT]
                            DEPTH1,                 & ! [INOUT]
                            ZNSNOW1,                & ! [INOUT]
                            nosnowsec   (i,j),      & ! [INOUT]
                            Salbedo     (i,j,I_SW), & ! [OUT]
                            Salbedo     (i,j,I_LW), & ! [OUT]
                            SFLX_SH     (i,j),      & ! [OUT]
                            SFLX_LH     (i,j),      & ! [OUT]
                            SFLX_GH     (i,j),      & ! [OUT]
                            SFLX_QV     (i,j),      & ! [OUT]
                            SFLX_QV_ENGI(i,j),      & ! [OUT]
                            QCC         (i,j),      & ! [OUT]
                            QFUSION     (i,j),      & ! [OUT]
                            MELT        (i,j),      & ! [OUT]
                            SWEMELT     (i,j),      & ! [OUT]
                            SNOW_LAND_GH(i,j),      & ! [OUT]
                            SFLX_water  (i,j),      & ! [IN]     ! [kg/m2/s]
                            SFLX_ENGI   (i,j),      & ! [IN]
                            TA          (i,j),      & ! [IN]
                            Uabs,                   & ! [IN]
                            RH,                     & ! [IN]
                            DENS        (i,j),      & ! [IN]
                            SFLX_SW_dn,             & ! [IN]
                            SFLX_LW_dn,             & ! [IN]
                            dt                      )


       SNOW_LAND_GH   (i,j) = SNOW_LAND_GH(i,j) / dt ! [J/m2] -> [J/m2/s]
       SNOW_LAND_Water(i,j) = SWEMELT(i,j) / dt      ! [kg/m2] -> [kg/m2/s]

       if ( SWE1 <= 0. .and. SWE(i,j) <= 0. ) then  ! no accumulated snow during the time step
          SNOW_frac      (i,j) = 0.0_RP
       else
          SNOW_frac      (i,j) = 1.0_RP
       endif

       TSNOW   (i,j) = TSNOW1
       SWE     (i,j) = SWE1
       SDepth  (i,j) = DEPTH1
       SDzero  (i,j) = ZNSNOW1

       if ( debug ) then
          LOG_INFO("LAND_PHY_SNOW_KY90",*) "SNOW_frac, SWE, TSNOW", SNOW_frac(i,j), SWE(i,j), TSNOW(i,j)
       end if

    else

       SNOW_LAND_Water(i,j)   = SFLX_water(i,j)
       SNOW_frac      (i,j)   = 0.0_RP

       TSNOW          (i,j)   = T0     !!!
       SWE            (i,j)   = 0.0_RP !!!
       SDepth         (i,j)   = 0.0_RP !!!
       SDzero         (i,j)   = 0.0_RP !!!
       Salbedo        (i,j,:) = 0.0_RP !!!
       SFLX_SH        (i,j)   = 0.0_RP
       SFLX_LH        (i,j)   = 0.0_RP
       SFLX_GH        (i,j)   = 0.0_RP
       SFLX_QV        (i,j)   = 0.0_RP
       QCC            (i,j)   = 0.0_RP
       QFUSION        (i,j)   = 0.0_RP
       MELT           (i,j)   = 0.0_RP
       SWEMELT        (i,j)   = 0.0_RP
       SNOW_LAND_GH   (i,j)   = 0.0_RP

    endif

    call FILE_HISTORY_in( QCC     (:,:), 'LAND_SNOW_QCC',        'Heat used for changing temperature profile', 'J/m2',  dim_type='XY' )
    call FILE_HISTORY_in( QFUSION (:,:), 'LAND_SNOW_QFUSION',    'Heat used for phase change of snow',         'J/m2',  dim_type='XY' )
    call FILE_HISTORY_in( MELT    (:,:), 'LAND_SNOW_MELT',       'Heat used for snow melt',                    'J/m2',  dim_type='XY' )
    call FILE_HISTORY_in( SWEMELT (:,:), 'LAND_SNOW_SWEMELT',    'Equivalent water of melt snow',              'kg/m2', dim_type='XY' )

    end do
    end do

   return
 end subroutine LAND_PHY_SNOW_KY90


  !-----------------------------------------------------------------------------
  !> snow model main routine
  subroutine SNOW_ky90_main ( &
       TSNOW,                 & ! [INOUT]
       SWE,                   & ! [INOUT]
       DEPTH,                 & ! [INOUT]
       ZNSNOW,                & ! [INOUT]
       nosnowsec,             & ! [INOUT]
       ALBEDO_out,            & ! [OUT]
       Emiss,                 & ! [OUT]
       HFLUX,                 & ! [OUT]
       LATENTFLUX,            & ! [OUT]
       GFLUX,                 & ! [OUT]
       EvapFLX,               & ! [OUT]
       Evap_ENGI,             & ! [OUT]
       QCC,                   & ! [OUT]
       QFUSION,               & ! [OUT]
       MELT,                  & ! [OUT]
       SWEMELT,               & ! [OUT]
       Gflux2land,            & ! [OUT]
       SFLX_SNOW,             & ! [IN]
       SFLX_ENGI,             & ! [IN]
       TA,                    & ! [IN]
       UA,                    & ! [IN]
       RH,                    & ! [IN]
       DENS,                  & ! [IN]
       SW,                    & ! [IN]
       LW,                    & ! [IN]
       time                   ) ! [IN]
    use scale_const, only:   &
       T0    => CONST_TEM00
    use scale_atmos_hydrometeor, only: &
       CV_ICE, &
       LHF

    implicit none
    ! prognostic variables
    real(RP), intent(inout)    :: TSNOW         ! snow temperature            [K]
    real(RP), intent(inout)    :: SWE           ! equivalent water [kg/m2]
    real(RP), intent(inout)    :: DEPTH         ! depth of melting point      [m]
    real(RP), intent(inout)    :: ZNSNOW        ! total snow depth            [m]
    ! update variables
    real(RP), intent(inout)    :: nosnowsec
    real(RP), intent(out)      :: ALBEDO_out
    real(RP), intent(out)      :: Emiss

    ! output variables
    real(RP), intent(out)      :: HFLUX         ! HFLUX = whole snow Sensible heat flux    [W/m2]
    real(RP), intent(out)      :: LATENTFLUX    ! LATENTFLUX = whole snow Latent heat flux [W/m2]
    real(RP), intent(out)      :: GFLUX         ! GFLUX = whole snow Ground flux           [W/m2]
    real(RP), intent(out)      :: EvapFLX       ! Evapolation due to LATENTFLUX            [kg/m2/s]
    real(RP), intent(out)      :: Evap_ENGI     ! internal energy flux of evapolation

    real(RP), intent(out)      :: QCC           ! QCC = heat used for change snow condition to isothermal [J m^-2]
    real(RP), intent(out)      :: QFUSION       ! QFUSION = heat used for change snow condition to melt point [J m^-2]
    real(RP), intent(out)      :: MELT          ! MELT = heat used for snow run off production [J m^-2 s^-1]
    real(RP), intent(out)      :: SWEMELT       ! snow water equivalent  ! [kg/m2]
    real(RP), intent(out)      :: Gflux2land    ! Residual heat, goes to land model ! [J/m2]


    ! input data
    real(RP), intent(in)       :: SFLX_SNOW
    real(RP), intent(in)       :: SFLX_ENGI
    real(RP), intent(in)       :: TA
    real(RP), intent(in)       :: UA
    real(RP), intent(in)       :: RH
    real(RP), intent(in)       :: DENS
    real(RP), intent(in)       :: SW
    real(RP), intent(in)       :: LW
    real(DP), intent(in)       :: time

    ! initial value
    real(RP)                   :: TSNOW0           ! Initial time snow surface temperature [K]
    real(RP)                   :: ZNSNOW0          ! ZNSNOW0 = initial freezing depth      [m]
    real(RP)                   :: SWE0             ! SWE0 = snow depth initial value in snow water equivalen [kg/m2]
    real(RP)                   :: DEPTH0           ! DEPTH0 = initial snow depth           [m]

    ! works
    real(RP)                   :: SNOW             ! per dt
    real(RP)                   :: RFLUX            ! RFLUX = whole snow net long wave radiation flux [W/m2]
    real(RP)                   :: LINFLUX          ! LINFLUX = whole snow downward long wave radiation flux [W/m2]
    real(RP)                   :: LOUTFLUX         ! LOUTFLUX = whole snow upward long wave radiation flux [W/m2]
    real(RP)                   :: SFLUX            ! SFLUX = whole snow net short wave radiation flux [W/m2]
    real(RP)                   :: DELTADEPTH

    real(RP)                   :: beta


!---------------------------------------------- !
!        ALL START HERE                         !
!---------------------------------------------- !

    ! initialize
    ZN_flag = 0
    TS_flag = 0

    QCC     = 0.0_RP
    QFUSION = 0.0_RP
    MELT    = 0.0_RP
    SWEMELT = 0.0_RP

    ! snowfall during timestep
    SNOW    = SFLX_SNOW * time

    ! save previous timestep
    TSNOW0  = TSNOW
    ZNSNOW0 = ZNSNOW
    SWE0    = SWE
    DEPTH0  = DEPTH

    ! update
    SWE0    = SWE0    +  SNOW              ! update according to snowfall
    DEPTH0  = DEPTH0  + (SNOW /RHOSNOW)
    ZNSNOW0 = ZNSNOW0 + (SNOW /RHOSNOW)

    if ( debug ) then
       LOG_INFO("SNOW_ky90_main",*) "UA, SNOW,SFLX_SNOW,time : ", UA, SNOW, SFLX_SNOW, time
       LOG_INFO("SNOW_ky90_main",*) "SWE , TSNOW, and TA :     ", SWE0, TSNOW0, TA
       LOG_INFO("SNOW_ky90_main",*) "DEPTH is:                 ", DEPTH0
       LOG_INFO("SNOW_ky90_main",*) "ZN beginning:             ", ZNSNOW0
    end if

!----- Calculating albedo -------------------------------------------!

    if (SNOW > 0.0_RP) then    ! snowfall
       nosnowsec = 0.0_RP
    else
       nosnowsec = nosnowsec + time
    endif
    if ( .not. ALBEDO_const ) then
       ALBEDO = ALBEDOMIN + (ALBEDOMAX-ALBEDOMIN)*exp(-1.0_RP*nosnowsec/real(4.0_RP*24.0_RP*3600.0_RP))
    endif

    ALBEDO_out = ALBEDO
    Emiss      = 1.0_RP - epsilon
    if ( debug ) then
       LOG_INFO("SNOW_ky90_main",*) "Albedo                    ",ALBEDO
    end if

!----- Energy balance at snow surface -------------------------------!

   call groundflux (TSNOW0, TA, UA, RH, DENS, ALBEDO, SW, LW, &  ! [IN]
                    GFLUX, RFLUX, SFLUX, LINFLUX, LOUTFLUX, HFLUX, LATENTFLUX) ! [OUT]

   ! SWE change due to latent heat flux
   EvapFLX = LATENTFLUX / LV                 ! [kg/m2/s] positive => snow decrease

   Evap_ENGI = EvapFLX * ( CV_ICE * TSNOW0 - LHF ) ! internal energy of evapolated water
   GFLUX = GFLUX + SFLX_ENGI - Evap_ENGI           ! add internal energy of precipitation and evapolation

   SWE0        = SWE0    - EvapFLX * time    ! [kg/m2]
   DELTADEPTH  = (EvapFLX * time) /RHOSNOW
   DEPTH0      = DEPTH0  - DELTADEPTH        ! [m]
   ZNSNOW0     = ZNSNOW0 - DELTADEPTH        ! [m]

!! Check whether GFLUX (energy into snowpack) is enough to melt all snow.
!! If GFLUX is enough, the model melts all snow and then go to next timestep.

   call check_allSnowMelt   (GFLUX, TSNOW0, ZNSNOW0, DEPTH0, sflag, time)
   if(sflag .eq. 1)then
      LOG_INFO("SNOW_ky90_main",*) 'LAND/snow: All snow melt'
      QCC        = 0.5_RP*CSRHOS*ZNSNOW0*(T0-TSNOW0)
      QFUSION    = W0*RHOSNOW*LF*ZNSNOW0
      MELT       = (1.0_RP-W0)*RHOSNOW*LF*DEPTH0  ! [J/m2]
      TSNOW      = T0
      ZNSNOW     = 0.0_RP
      DEPTH      = 0.0_RP
      SWE        = 0.0_RP
      !SWEMELT    = MELT /((1.0_RP-W0)*LF)
      SWEMELT    = RHOSNOW*DEPTH0
      Gflux2land = GFLUX*time - (QCC + QFUSION + MELT)   ! [J/m2]

   else

      call cal_param (ZNSNOW0, TSNOW0, GFLUX, TA, UA, RH, DENS, LW, time)

      ! check whether the model has solution
      call check_applicability (GFLUX, TSNOW0, ZNSNOW0, TA, UA, RH, DENS, LW, Gres, beta, time)
      if ((Gres >= 0.0_RP).and.(beta >= 0.0_RP)) then
         LOG_INFO("SNOW_ky90_main",*) 'LAND/snow model is not appropriate',Gres,beta
         QCC             = 0.5_RP*CSRHOS*ZNSNOW0*(T0-TSNOW0)
         QFUSION         = W0*RHOSNOW*LF*ZNSNOW0
         MELT            = Gres
         TSNOW           = T0
         ZNSNOW          = 0.0_RP
      else

         !if (t==1)then
         !   call cal_R1R2 (ZNSNOW0, TSNOW0, GFLUX, TA, UA, RH, rhoair, LW, time)
         !endif

         call snowdepth ( GFLUX, ZNSNOW0, ZNSNOW, time)

         if ( debug ) then
            LOG_INFO("SNOW_ky90_main",*) "ZN is: ", ZNSNOW
         end if
         IF(ZNSNOW < ZMIN) THEN
            ZN_flag = 1
            if ( debug ) then
               LOG_INFO("SNOW_ky90_main",*) "ZN is replaced to: ", ZNSNOW ," to ", ZMIN
            end if
            ZNSNOW = ZMIN
         ELSE IF(ZNSNOW > DEPTH0) THEN
            ZN_flag = 2
            if ( debug ) then
               LOG_INFO("SNOW_ky90_main",*) "ZN is replaced to: ", ZNSNOW ," to ", DEPTH0
            end if
            ZNSNOW = DEPTH0
         END IF

         ! This equation is to calculate TSN
         call  equation415(LAMBDAS, C2, ZNSNOW, RH, QSAT, TSNOW0, ZNSNOW0, GFLUX, TA, UA, DENS, LW, TSNOW, time)

         if ( debug ) then
            LOG_INFO("SNOW_ky90_main",*) 'TSNOW is:       ', TSNOW
         end if

         if (TSNOW > T0) then
            TS_flag = 1
            TSNOW   = T0
            call recalculateZ(ZNSNOW0, TSNOW0, GFLUX, ZNSNOW, time)
            call check_res(ZNSNOW0, ZNSNOW, TSNOW0, TSNOW, GFLUX, TA, UA, RH, DENS, LW, "1", time)
            IF (ZNSNOW < ZMIN) THEN
               ZN_flag = 4
               if ( debug ) then
                  LOG_INFO("SNOW_ky90_main",*) "ZN is updated/replaced to: ", ZNSNOW ," to ", ZMIN
               end if
               ZNSNOW = ZMIN
            ELSE IF(ZNSNOW > DEPTH0) THEN
               ZN_flag = 5
               if ( debug ) then
                  LOG_INFO("SNOW_ky90_main",*) "ZN is updated/replaced to: ", ZNSNOW ," to ", DEPTH0
               end if
               ZNSNOW = DEPTH0
            ELSE
               ZN_flag = 3
            END IF

            if ( debug ) then
               LOG_INFO("SNOW_ky90_main",*) 'TSNOW  is updated:  ', TSNOW
               LOG_INFO("SNOW_ky90_main",*) 'ZNSNOW is updated:  ', ZNSNOW
            end if
         else
            call check_res( ZNSNOW0, ZNSNOW, TSNOW0, TSNOW, GFLUX, TA, UA, RH, DENS, LW, "0", time)
         endif

         IF ( (ZNSNOW-ZMIN)/ZMIN < 0.00001 ) THEN
            call calculationMO(  &
                 GFLUX, CSRHOS, ZNSNOW0, TSNOW0, ZNSNOW, TSNOW, &
                 MELT, QCC, QFUSION, time)
         ELSE
            call calculationNoMO(  &
                 GFLUX, CSRHOS, ZNSNOW0, TSNOW0, ZNSNOW, TSNOW,  &
                 MELT, QCC, QFUSION, time)
         END IF

      endif ! Gres & beta

      Gflux2land = GFLUX*time - (QCC + QFUSION + MELT)
      if ( debug ) then
         LOG_INFO("SNOW_ky90_main",*) "### ZN_flag = ", ZN_flag, "TS_flag = ", TS_flag
         LOG_INFO("SNOW_ky90_main",*) '### Heat flux from snowpack to land surface: ', Gflux2land
      end if

      DELTADEPTH             = MELT / ((1.0_RP-W0)*LF*RHOSNOW)
      SWEMELT                = RHOSNOW*DELTADEPTH
      SWE                    = SWE0        - SWEMELT
      DEPTH                  = DEPTH0      - DELTADEPTH
      if (DEPTH < ZNSNOW) then
         ! NOTICE: energy budget has not been considered thought this process yet.
         if ( debug ) then
            LOG_INFO("SNOW_ky90_main",*) "replace ZNSNOW <= DEPTH"
         end if
         ZNSNOW               = DEPTH
      endif

      if ( debug ) then
         LOG_INFO("SNOW_ky90_main",*) 'MELT in water equivalent is: ', SWEMELT
         LOG_INFO("SNOW_ky90_main",*) 'SWE0 is:                     ', SWE
         LOG_INFO("SNOW_ky90_main",*) 'DELTADEPTH is:               ', DELTADEPTH
         LOG_INFO("SNOW_ky90_main",*) 'DEPTH0 is:                   ', DEPTH
         LOG_INFO("SNOW_ky90_main",*) 'ZNSNOW0 is:                  ', ZNSNOW
      end if

   endif

   return
 end subroutine SNOW_ky90_main

!==============================================================

subroutine groundflux (TS, TA, UA, RH, rhoair, ALPHA, SW, LW, &
                       GFLUX, RFLUX, SFLUX, LINFLUX, LOUTFLUX, HFLUX, LATENTFLUX)
  implicit none

  real(RP), intent(in)     :: TS, TA, UA, RH, rhoair, ALPHA, SW, LW
  real(RP), intent(out)    :: GFLUX, RFLUX, SFLUX, LINFLUX, LOUTFLUX, HFLUX, LATENTFLUX

  ESAT               = 0.6112_RP * exp( (17.67_RP * (TA-273.15_RP)) / (TA-29.66_RP) )
  QSAT               = 0.622_RP * ESAT / 101.325_RP
  DELTAQSAT          = 0.622_RP * LV * QSAT /(287.04_RP * (TA**2))
  CP                 = 1004.67_RP * (1.0_RP + (0.84_RP * QSAT))

  ! Into snowpack is positive
  SFLUX              = (1.0_RP-ALPHA) * SW
  RFLUX              = epsilon * (LW-(sigma*(TS**4)))
  LINFLUX            = epsilon * LW
  LOUTFLUX           = -1.0_RP * (epsilon * sigma * (TS**4))
  HFLUX              = cp * rhoair * CH * UA * (TS-TA)
  LATENTFLUX         = LV * rhoair * CE * UA * ((1-RH)*QSAT + (DELTAQSAT*(TS-TA)))

  GFLUX              = (SFLUX + RFLUX - HFLUX - LATENTFLUX)

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_groundflux",*) "-------------- groundflux --------------"
     LOG_INFO_CONT(*) "GFLUX is:    ", GFLUX
     LOG_INFO_CONT(*) "SFLUX:       ", SFLUX
     LOG_INFO_CONT(*) "RFLUX:       ", RFLUX, LINFLUX+LOUTFLUX
     LOG_INFO_CONT(*) " (LONG in:   ", LINFLUX,")"
     LOG_INFO_CONT(*) " (LONG out:  ", LOUTFLUX,")"
     LOG_INFO_CONT(*) "HFLUX is:    ", HFLUX
     LOG_INFO_CONT(*) "LATENT FLUX: ", LATENTFLUX
  end if

return
end subroutine groundflux

!==============================================================
subroutine check_allSnowMelt (GFLUX, TS1, ZN1, D, sflag, time)
  use scale_const, only:   &
       T0    => CONST_TEM00

  implicit none
  real(RP), intent(in) :: GFLUX, TS1, ZN1, D
  real(DP), intent(in) ::time
  integer, intent(out) :: sflag
  real(RP)             :: energy_in, energy_use
  real(RP)             :: energy_use_ripe, energy_use_melt

  energy_in  = GFLUX * time

  energy_use_ripe = 0.5_RP*CSRHOS*ZN1*(T0-TS1) + W0*RHOSNOW*LF*ZN1
  energy_use_melt = (1.0_RP-W0)*RHOSNOW*LF*D
  energy_use      = energy_use_ripe + energy_use_melt

  if(energy_in >= energy_use)then
     sflag=1
  else
     sflag=0
  endif

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_check_allSnowMelt",*) "Energy in  =",energy_in
     LOG_INFO("LAND_PHY_SNOW_KY90_check_allSnowMelt",*) "Energy use =",energy_use
  end if

  return
end subroutine

!==============================================================
subroutine cal_param (ZN1, TS1, GFLUX, TA, UA, RH, rhoair, LW, time)
  use scale_const, only:   &
       T0    => CONST_TEM00

  implicit none
  real(RP), intent(in) ::  ZN1, TS1, TA, UA, RH, rhoair, LW
  real(RP), intent(in) ::  GFLUX
  real(DP), intent(in) :: time

  C1      = csrhos*0.5_RP
  C2      = (4.0_RP*epsilon*sigma*(TA**3)) + (cp*rhoair*CH*UA) + (LV*rhoair*CE*UA*DELTAQSAT)
  C3      = W0*RHOSNOW*LF

  A0      = LAMBDAS*(C3*ZN1 + (C1*ZN1*(T0-TS1)- GFLUX*time))
  A1      = C2*(C3*ZN1 + C1*ZN1*(T0-TS1) - GFLUX*time) - C3*LAMBDAS
  A2      = C1*  &
            ( epsilon*(LW-(sigma*(TA**4))) + (C2*(TA-T0)) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) ) &
            - C2*C3

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_cal_param",*) "-------------- snowdepth --------------"
     LOG_INFO_CONT(*) "C1",C1
     LOG_INFO_CONT(*) "C2",C2,(4.0_RP*epsilon*sigma*(TA**3)),(cp*rhoair*CH*UA), (LV*rhoair*CE*UA*DELTAQSAT)
     LOG_INFO_CONT(*) "C3",C3
     LOG_INFO_CONT(*) "A0",A0
     LOG_INFO_CONT(*) "A1",A1
     LOG_INFO_CONT(*) "A2",A2
  end if

end subroutine

!==============================================================
subroutine check_applicability (GFLUX, TS1, ZN1, TA, UA, RH, rhoair, LW, GFLUX_res, beta, time)
  use scale_const, only:   &
       T0    => CONST_TEM00

  implicit none
  real(RP), intent(in)  :: GFLUX, TS1, ZN1
  real(RP), intent(in)  :: TA, UA, RH, rhoair, LW
  real(DP), intent(in)  :: time
  real(RP), intent(out) :: GFLUX_res, beta
  real(RP)              :: energy_in, energy_use_max

  energy_in      = GFLUX * time
  energy_use_max = 0.5_RP*CSRHOS*ZN1*(T0-TS1) + W0*RHOSNOW*LF*ZN1

  GFLUX_res      = energy_in - energy_use_max  ! residual energy after being used to melt all snow

  !
  beta = (LAMBDAS/C2)* &
   ( T0 - ( TA*C2 - LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT + epsilon*(LW-sigma*(TA**4)))/C2 )

  return
end subroutine

!==============================================================!
subroutine snowdepth (GFLUX, ZN1, ZN2, time)

  implicit none

  real(RP), intent(in)  :: GFLUX
  real(RP), intent(in)  :: ZN1
  real(DP), intent(in)  ::  time
  real(RP), intent(out) :: ZN2

!calcultaing snowdepth

   !print *, -1.0_RP*A1
   !print *, -1.0_RP*((A1**2.0_RP - 4.0_RP * A2 * A0)**0.5_RP)
   !print *, 2.0_RP*A2
   !print *, A1**2.0_RP, -1.0_RP * 4.0_RP * A2 * A0

   ZN2  = ((-1.0_RP*A1) - ((A1**2.0_RP - 4.0_RP*A2*A0)**0.5_RP)) / (2.0_RP*A2)

   if ( debug ) then
      LOG_INFO("LAND_PHY_SNOW_KY90_snowdepth",*) "ZN old = ",ZN1, "ZN new = ",ZN2
   end if

return
end subroutine snowdepth
!==============================================================
subroutine equation415(LAMBDAS, C2, ZN2, RH, QSAT, TS1, ZN1, GFLUX, TA, UA, rhoair, LW, TS2, time)
  use scale_const, only:   &
       T0    => CONST_TEM00
  implicit none

  real(RP), intent(in)  :: TA, UA, rhoair, LW, ZN2, GFLUX, TS1,ZN1
  real(RP), intent(in)  :: C2, LAMBDAS, RH, QSAT
  real(DP), intent(in)  :: time
  real(RP), intent(out) :: TS2
  real(RP)              :: TS_check

  TS2      = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
           + ZN2*epsilon*(LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_equation415",*) "-------------- equation415 --------------"

     TS_check = GFLUX*time/(0.5*CSRHOS*ZN2) + T0 - ZN1*(T0-TS1)/ZN2 - W0*RHOSNOW*LF*(ZN1-ZN2)/(0.5*CSRHOS*ZN2)

     LOG_INFO_CONT(*) "compare ",TS2, TS_check, TS2-TS_check
     ! When ZN2 is replaced, TS2 does not equal to TS_check.
  end if

 return
end subroutine equation415
!==============================================================
subroutine recalculateZ(ZN1, TS, GFLUX, ZN2, time)
  use scale_const, only:   &
       T0    => CONST_TEM00
  implicit none

  real(RP), intent(in)  :: ZN1
  real(RP), intent(in)  :: TS, GFLUX
  real(DP), intent(in)  :: time
  real(RP), intent(out) :: ZN2

  ZN2 = ZN1 + ((C1*ZN1*(T0-TS) - GFLUX*time) / C3)

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_recalculateZ",*) "-------------- recalculate Z --------------"
  end if

 return
end subroutine recalculateZ

!==============================================================
subroutine calculationMO(GFLUX, CSRHOS, ZN1, TS1, ZN2, TS2, &
                         MELT, QCC, QFUSION, time)
  use scale_const, only:   &
       T0    => CONST_TEM00
  use scale_prc, only: &
       PRC_abort
  implicit none
  real(RP), intent(in)  :: GFLUX, CSRHOS, ZN1, TS1, TS2, ZN2
  real(DP), intent(in)  :: time
  real(RP), intent(out) :: MELT, QCC, QFUSION

  QCC             = 0.5_RP*CSRHOS*(ZN1*(T0-TS1)-ZN2*(T0-TS2))
  QFUSION         = W0*RHOSNOW*LF*(ZN1-ZN2)
  MELT            = ( GFLUX*time - QCC - QFUSION )

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_calculationMO",*) "--------------------------MELT----------------"
     LOG_INFO_CONT(*) "GFLUX*time is: ", GFLUX*time
     LOG_INFO_CONT(*) "QCC is       : ", QCC
     LOG_INFO_CONT(*) "QFUSION is   : ", QFUSION
     LOG_INFO_CONT(*) "QMELT is     : ", MELT

     LOG_INFO_CONT(*) QCC+QFUSION+MELT
     LOG_INFO_CONT(*) "diff= ", QCC + QFUSION + MELT - (GFLUX*time)
  end if

  if ( ABS(QCC+QFUSION+MELT - (GFLUX*time)) > 10.) then
     LOG_ERROR("LAND_PHY_SNOW_KY90_calculationMO",*) "Calculation is fault. Model would include bugs. Please check! Melt"
     call PRC_abort
  endif


 return
end subroutine calculationMO

!==============================================================
subroutine calculationNoMO(GFLUX, CSRHOS, ZN1, TS1, ZN2, TS2, &
                           MELT, QCC, QFUSION, time)
  use scale_const, only:   &
       T0    => CONST_TEM00
 implicit none
 real(RP), intent(in)  :: GFLUX, CSRHOS, ZN1, TS1, TS2, ZN2
 real(DP), intent(in)  :: time
 real(RP), intent(out) :: MELT, QCC, QFUSION

  QCC             = 0.5_RP*CSRHOS*(ZN1*(T0-TS1)-ZN2*(T0-TS2))
  QFUSION         = W0*RHOSNOW*LF*(ZN1-ZN2)
  MELT            = 0.0_RP

  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_calculationNoMO",*) "--------------------------NOMELT----------------"
     LOG_INFO_CONT(*) "GFLUX*time is: ", GFLUX*time
     LOG_INFO_CONT(*) "QCC is       : ", QCC
     LOG_INFO_CONT(*) "QFUSION is   : ", QFUSION
     LOG_INFO_CONT(*) "QMELT is     : ", MELT

     LOG_INFO_CONT(*) QCC+QFUSION+MELT
     LOG_INFO_CONT(*) "diff= ", QCC +QFUSION - (GFLUX*time)
  end if

  !if ( ABS(QCC+QFUSION - (GFLUX*time)) > 10.) then
  !  LOG_ERROR("LAND_PHY_SNOW_KY90_calculationNoMO",*) "Calculation is fault. Model would include bugs. Please check! No Melt"
  !  call PRC_abort
  !endif

 return
end subroutine calculationNoMO

!==============================================================
subroutine check_res(ZN1, ZN2, TS1, TS2, GFLUX, TA, UA, RH, rhoair, LW, flag, time)
  use scale_const, only:   &
       T0    => CONST_TEM00
  implicit none
  real(RP), intent(in) :: ZN1, ZN2, TS1, TS2, GFLUX
  real(RP), intent(in) :: TA, UA, RH, rhoair, LW
  real(DP), intent(in) :: time
  character(len=*)     :: flag
  real(RP)             :: R1 ! Eq. (2)
  real(RP)             :: R2 ! Eq. (8)
  real(RP)             :: R3 ! Eq. (8)

  R3 = -999.
  R1 = C1 * (ZN1*(T0-TS1) - ZN2*(T0-TS2)) + C3*(ZN1-ZN2)-GFLUX*time
  R2 = ZN2*( epsilon*LW-epsilon*sigma*(TA**4) - (C2*(TS2-TA)) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) ) &
       + LAMBDAS*(T0-TS2)
  R3   = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
           + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))-TS2


  if ( debug ) then
     LOG_INFO("LAND_PHY_SNOW_KY90_check_res",*) "R1 is         : ", R1, "flag = ", flag
     LOG_INFO("LAND_PHY_SNOW_KY90_check_res",*) "R2 is         : ", R2, R3,"flag = ", flag

     if(abs(R1)>10000)then
        LOG_INFO("LAND_PHY_SNOW_KY90_check_res",*) C1 * (ZN1*(T0-TS1) - ZN2*(T0-TS2)), C3*(ZN1-ZN2), -GFLUX*time
     endif
  end if

  return
end subroutine check_res

!==============================================================
subroutine cal_R1R2(ZN1, TS1, GFLUX, TA, UA, RH, rhoair, LW, time)
  use scale_const, only:   &
       T0    => CONST_TEM00
  implicit none
  real(RP), intent(in) :: ZN1, TS1
  real(RP), intent(in) :: GFLUX, TA, UA, RH, rhoair, LW
  real(DP), intent(in) :: time
  real(RP)             :: ZN2, TS2
  real(RP)             :: R1 ! Eq. (2)
  real(RP)             :: R2 ! Eq. (8)
  real(RP)             :: R3 ! Eq. (8)

  real(RP)             :: ts0,zn0
  integer              :: i,j
  real(RP)             :: ts_r1,ts_r2,zn_r1,zn_r2
  character(len=3)     :: ttt = ""

  real(RP)               :: a,b,c,d,e,f,g

  ts0 = -50.0_RP + T0   ! -25 to + 25
  zn0 = -10.0_RP        ! -50 to + 50

  !write(ttt,'(i3.3)') t

  open(70,file='check_R1-R2_zn-base'//ttt//'.dat',status='unknown')
  do j=1,10000
     ZN2 = zn0 + 9.0_RP*0.0001_RP*real((j-1),kind=RP)

     if (ZN2/=0.0)then
        ts_r1 = T0 - (C1*ZN1*(T0-TS1) + C3*(ZN1-ZN2) - GFLUX*time)/(C1*ZN2)

        ts_r2 = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
             + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))
     endif
     if ( debug ) write(70,'(3f15.5)') ZN2, ts_r1, ts_r2
  enddo
  do j=1,100000
     ZN2 = -1.0_RP + 2.0_RP*0.00001_RP*real((j-1),kind=RP)

     if (ZN2/=0.0)then
        ts_r1 = T0 - (C1*ZN1*(T0-TS1) + C3*(ZN1-ZN2) - GFLUX*time)/(C1*ZN2)

        ts_r2 = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
             + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))
     endif
     if ( debug ) write(70,'(3f15.5)') ZN2, ts_r1, ts_r2
  enddo
  do j=2,10000
     ZN2 = 1.0_RP + 9.0_RP*0.0001_RP*real((j-1),kind=RP)

     if (ZN2/=0.0)then
        ts_r1 = T0 - (C1*ZN1*(T0-TS1) + C3*(ZN1-ZN2) - GFLUX*time)/(C1*ZN2)

        ts_r2 = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
             + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))
     endif
     if ( debug ) write(70,'(3f15.5)') ZN2, ts_r1, ts_r2
  enddo
  close(70)

  a=T0+C3/C1
  b=-1.0_RP * (C1*ZN1*(T0-TS1) + C3*ZN1 - GFLUX*time)
  c=C1

  d=(LAMBDAS*T0)
  e=( TA*C2 - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) + epsilon*( LW - (sigma*(TA**4)))) ! *ZN2
  f=(LAMBDAS + (C2*ZN2))


  if ( debug ) then
     open(71,file='check_R1-R2-grad'//ttt//'.dat',status='unknown')
     write(71,'(5f15.5)') b/c, d,e,LAMBDAS,C2
     close(71)

     open(70,file='check_R1-R2_ts-base'//ttt//'.dat',status='unknown')
     do i=1,100000
        TS2 = ts0 + 150.0_RP*0.00001_RP*real((i-1),kind=RP)

        zn_r1 = (C1*ZN1*(T0-TS1) + C3*ZN1 - GFLUX*time)/(C1*(T0-TS2)+C3)
        zn_r2 = -1.0_RP*LAMBDAS*(T0-TS2)/(LW*epsilon-epsilon*sigma*(TA**4) - (C2*(TS2-TA)) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT))

        write(70,'(3f15.5)') TS2, zn_r1, zn_r2
     enddo
     close(70)
  end if

  !LOG_INFO("cal_R1R2",*) "aa",(LINFLUX-epsilon*sigma*(TA**4) + (C2*TA) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT))/C2
  return
end subroutine cal_R1R2

!!!!!!!!!!!!!!!!!!!!!!SUBROUTINE END!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module scale_land_phy_snow_ky90
