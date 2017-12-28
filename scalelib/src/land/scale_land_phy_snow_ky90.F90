!-------------------------------------------------------------------------------
!> module LAND / SNOW model for Land Slab model
!!
!! @par Description
!!       Profile model for snowpack by Kondo and Yamazaki (1990)
!!       Kondo and Yamazaki (1990): Journal of Applied Meteorology, 29, pp375-384.
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_phy_snow_KY90
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !use scale_prof
  !use scale_debug
  use scale_grid_index
  !use scale_land_grid_index
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
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !integer, private, parameter        :: RP=8

  !
  real(RP)                  :: A0, A1, A2, C1, C2, C3
  real(RP)                  :: ESAT, QSAT, DELTAQSAT, CP
  !integer                   :: t
  real(RP)                  :: Gres

  integer, parameter        :: data_length_max=10000

  ! model parameters
  real(RP)                  :: LAMBDAS               ! Snow thermal conductivity [W m^-1 K^-1]
  real(RP)                  :: W0                    ! Maximum water content [%]
  real(RP)                  :: CSRHOS                ! Heat capacity [J/m^3/K]
                                                     ! (Specific heat of snow [J kg^-1 K^-1])*(Snow density[kg m^-3])
  real(RP)                  :: RHOSNOW   = 400.0_RP  ! Snow density [kg m^-3]
  real(RP)                  :: CS                    ! Specific heat for unit mass [J/kg/K]
  real(RP)                  :: ALBEDO                ! albedo
  real(RP)                  :: ALBEDOMIN = 0.4_RP
  real(RP)                  :: ALBEDOMAX = 0.85_RP
  real(RP)                  :: ZMIN      = 0.01_RP   ! Minimum snow depth [m]

  ! constant variables
  real(RP), parameter       :: epsilon   = 0.97_RP
  real(RP), parameter       :: sigma     = 5.67e-08_RP
  real(RP), parameter       :: rhoair    = 1.289_RP      ! [kg/m3]
  real(RP), parameter       :: CH        = 0.002_RP
  real(RP), parameter       :: CE        = 0.0021_RP
  real(RP), parameter       :: LV        = 2.5e6_RP
  real(RP), parameter       :: LF        = 3.34e5_RP
  real(RP), parameter       :: T0        = 273.15_RP

  logical                   :: ALBEDO_const = .true.

  integer                   :: ZN_flag, TS_flag, sflag

  !----------------------------------------------------------------------------------------------!
contains
  !----------------------------------------------------------------------------------------------!
  !> Setup
  subroutine LAND_PHY_SNOW_KY90_setup
    use scale_process, only: &
       PRC_MPIstop
    !use mod_land_vars, only: &  ! tentative
    !   SNOW_TEMP,         &
    !   SNOW_SWE,          &
    !   SNOW_Depth,        &
    !   SNOW_Dzero,        &
    !   SNOW_nosnowsec
    implicit none

    ! initial value
    !real(RP), intent(out)      :: TSNOW0   ! Initial time snow surface temperature [K]
    !real(RP), intent(out)      :: ZNSNOW0 ! ZNSNOW0 = initial freezing depth      [m]
    !real(RP), intent(out)      :: SWE0    ! SWE0 = snow depth initial value in snow water equivalen [kg/m2]
    !real(RP), intent(out)      :: DEPTH0  ! DEPTH0 = initial snow depth           [m]
    !real(RP)                   :: nosnowsec                   ! number of hours from latest snow event

    real(RP)                  :: snow_conductivity     = 0.42_RP
    real(RP)                  :: water_content         = 0.1_RP
    real(RP)                  :: snow_heat_capacityRHO = 8.4e+05_RP
    real(RP)                  :: snow_rho              = 400.0_RP
    real(RP)                  :: snowDepth_initial     = 0.0_RP
    real(RP)                  :: albedo_value          = 0.686_RP

    namelist / PARAM_LAND_PHY_SNOW_KY90 /  &
         ALBEDO_const,      &
         snow_conductivity, &
         water_content,     &
         snow_heat_capacityRHO, &
         snow_rho,          &
         snowDepth_initial, &
         albedo_value

    integer :: ierr
    !---------------------------------------------------------------------------

    !open(10,file='variable_snow2.in',status='OLD')
    !read(10,nml=variables)
    !write(6,nml=variables)
    !close(10)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SNOW_KY90] / Categ[LAND PHY] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PHY_SNOW_KY90,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_PHY_SNOW_KY90. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_PHY_SNOW_KY90)

    LAMBDAS       = snow_conductivity
    W0            = water_content
    CSRHOS        = snow_heat_capacityRHO ! [J/m3/K]
    RHOSNOW       = snow_rho
    CS            = CSRHOS/RHOSNOW        ! [J/kg/K]
    !write(*,*)   "Specific heat capacity of snow [J/kg/K]: ",CS
    ALBEDO        = albedo_value


    ! Set initial value: tentative
    !SNOW_TEMP      = T0
    !SNOW_SWE       = snowDepth_initial*RHOSNOW
    !SNOW_Depth     = snowDepth_initial
    !SNOW_Dzero     = SNOW_Depth
    !SNOW_nosnowsec = 0.0_RP

    return
  endsubroutine LAND_PHY_SNOW_KY90_setup

  !-----------------------------------------------------------------------------
  !> Main routine for land submodel

  subroutine LAND_PHY_SNOW_KY90(     &
             TSNOW,                  & ! [INOUT]
             SWE,                    & ! [INOUT]
             SDepth,                 & ! [INOUT]
             SDzero,                 & ! [INOUT]
             Salbedo,                & ! [INOUT]
             nosnowsec,              & ! [INOUT]
             TSNOW_t,                & ! [OUT]
             SFLX_SH,                & ! [OUT]
             SFLX_LH,                & ! [OUT]
             SFLX_GH,                & ! [OUT]
             SNOW_LAND_GH,           & ! [OUT]
             SNOW_LAND_Water,        & ! [OUT]
             SNOW_frac,              & ! [OUT]
             SFLX_rain,              & ! [IN]
             SFLX_snow,              & ! [IN]
             PRSA,                   & ! [IN]
             TA,                     & ! [IN]
             WA,                     & ! [IN]
             UA,                     & ! [IN]
             VA,                     & ! [IN]
             QA,                     & ! [IN]
             SFLX_SW_dn,             & ! [IN]
             SFLX_LW_dn,             & ! [IN]
             dt                      ) ! [IN]
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_saturation, only:  &
       qsatf => ATMOS_SATURATION_pres2qsat_all  ! better to  change name from qsatf to qsat
    use scale_landuse, only: &
       LANDUSE_fact_land

    implicit none
    ! prognostic variables
    real(RP), intent(inout)   :: TSNOW          (IA,JA) ! snow temperature        [K]
    real(RP), intent(inout)   :: SWE            (IA,JA) ! equivalent water        [kg/m2]
    real(RP), intent(inout)   :: SDepth         (IA,JA) ! depth of melting point  [m]
    real(RP), intent(inout)   :: SDzero         (IA,JA) ! total snow depth        [m]
    real(RP), intent(inout)   :: Salbedo        (IA,JA) ! snow albedo             [m]
    real(RP), intent(inout)   :: nosnowsec      (IA,JA) ! elapsed time of no snow condition [s]

    ! updated variables
    real(RP), intent(out)     :: TSNOW_t        (IA,JA) ! updated tendency of snow temperature [K]
    real(RP), intent(out)     :: SFLX_SH        (IA,JA) ! sensible heat flux between atmos and snow [W/m2]
    real(RP), intent(out)     :: SFLX_LH        (IA,JA) ! latente  heat flux between atmos and snow [W/m2]
    real(RP), intent(out)     :: SFLX_GH        (IA,JA) ! whole snowpack Ground flux [W/m2]
    real(RP), intent(out)     :: SNOW_LAND_GH   (IA,JA) ! heat flux from snow to land [W/m2]
    real(RP), intent(out)     :: SNOW_LAND_Water(IA,JA) ! water flux from snow to land [W/m2]
    real(RP), intent(out)     :: SNOW_frac      (IA,JA) ! snow fraction, defined by time direction [-]

    ! input data
    real(RP), intent(in)      :: SFLX_rain (IA,JA)
    real(RP), intent(in)      :: SFLX_snow (IA,JA)
    real(RP), intent(in)      :: PRSA      (IA,JA)
    real(RP), intent(in)      :: TA        (IA,JA)
    real(RP), intent(in)      :: WA        (IA,JA)
    real(RP), intent(in)      :: UA        (IA,JA)
    real(RP), intent(in)      :: VA        (IA,JA)
    real(RP), intent(in)      :: QA        (IA,JA)      ! specific humidity [kg/kg]
    real(RP), intent(in)      :: SFLX_SW_dn(IA,JA)
    real(RP), intent(in)      :: SFLX_LW_dn(IA,JA)
    real(RP), intent(in)      :: dt                     ! dt of land

    real(RP)                  :: QCC       (IA,JA)
    real(RP)                  :: QFUSION   (IA,JA)
    real(RP)                  :: MELT      (IA,JA)
    real(RP)                  :: SWEMELT   (IA,JA)

    ! works
    real(RP)                  :: TSNOW1           ! updated snow surface temperature [K]
    real(RP)                  :: ZNSNOW1          ! updated freezing depth      [m]
    real(RP)                  :: SWE1             ! updated snow water equivalen [kg/m2]
    real(RP)                  :: DEPTH1           ! updated snow depth           [m]

    real(RP), parameter       :: Uabs_min = 0.1_RP
    real(RP)                  :: Uabs
    real(RP)                  :: QAsat
    real(RP)                  :: RH

    integer :: k, i, j
    !---------------------------------------------------------------------------
    if( IO_L ) write(IO_FID_LOG,*) '*** Snow surface physics step: SNOW KY90'

    do j = JS, JE
    do i = IS, IE

    if( ( LANDUSE_fact_land(i,j) > 0.0_RP    ) .and.    &
        ( SWE(i,j)>0. .or. SFLX_snow(i,j)>0. ) )then

       Uabs = max( sqrt( UA(i,j)**2 + VA(i,j)**2 + WA(i,j)**2 ), Uabs_min )
       !Uabs = sqrt( UA(i,j)**2 + VA(i,j)**2 + WA(i,j)**2 )

       call qsatf( QAsat, TA(i,j), PRSA(i,j) )
       RH  = QA(i,j) / QAsat
       write(*,*) "RH,   ",RH

       TSNOW1  = TSNOW (i,j)
       SWE1    = SWE   (i,j)
       DEPTH1  = SDepth(i,j)
       ZNSNOW1 = SDzero(i,j)

       call SNOW_ky90_main( TSNOW1,            & ! [INOUT]
                            SWE1,              & ! [INOUT]
                            DEPTH1,            & ! [INOUT]
                            ZNSNOW1,           & ! [INOUT]
                            nosnowsec   (i,j), & ! [INOUT]
                            Salbedo     (i,j), & ! [INOUT]
                            SFLX_SH     (i,j), & ! [OUT]
                            SFLX_LH     (i,j), & ! [OUT]
                            SFLX_GH     (i,j), & ! [OUT]
                            QCC         (i,j), & ! [OUT]
                            QFUSION     (i,j), & ! [OUT]
                            MELT        (i,j), & ! [OUT]
                            SWEMELT     (i,j), & ! [OUT]
                            SNOW_LAND_GH(i,j), & ! [OUT]
                            SFLX_snow   (i,j), & ! [IN]     ! [kg/m2/s]
                            TA          (i,j), & ! [IN]
                            Uabs,              & ! [IN]
                            RH,                & ! [IN]
                            SFLX_SW_dn  (i,j), & ! [IN]
                            SFLX_LW_dn  (i,j), & ! [IN]
                            dt                 )


       SNOW_LAND_GH   (i,j) = SNOW_LAND_GH(i,j) / dt               ! [J/m2] -> [J/m2/s]
       SNOW_LAND_Water(i,j) = SFLX_rain(i,j) + SWEMELT(i,j) / dt  ! [kg/m2] -> [kg/m2/s]

       if ( SWE1 <= 0. .and. SWE(i,j) <= 0. ) then  ! no accumulated snow during the time step
          SNOW_frac      (i,j) = 0.0_RP
       else
          SNOW_frac      (i,j) = 1.0_RP
       endif

       TSNOW_t (i,j) = ( TSNOW1 - TSNOW (i,j) ) / dt
       TSNOW   (i,j) = TSNOW1
       SWE     (i,j) = SWE1
       SDepth  (i,j) = DEPTH1
       SDzero  (i,j) = ZNSNOW1

       write(*,*) "SNOW_frac, SWE, TSNOW", SNOW_frac(i,j), SWE(i,j), TSNOW(i,j)

    else

       SNOW_LAND_Water(i,j) = SFLX_rain(i,j)
       SNOW_frac      (i,j) = 0.0_RP

       TSNOW_t (i,j)        = 0.0_RP
       TSNOW          (i,j) = T0     !!!
       SWE            (i,j) = 0.0_RP !!!
       SDepth         (i,j) = 0.0_RP !!!
       SDzero         (i,j) = 0.0_RP !!!
       Salbedo        (i,j) = 0.0_RP !!!
       SFLX_SH        (i,j) = 0.0_RP
       SFLX_LH        (i,j) = 0.0_RP
       SFLX_GH        (i,j) = 0.0_RP
       QCC            (i,j) = 0.0_RP
       QFUSION        (i,j) = 0.0_RP
       MELT           (i,j) = 0.0_RP
       SWEMELT        (i,j) = 0.0_RP
       SNOW_LAND_GH   (i,j) = 0.0_RP

    endif

    ! HIST_in

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
       ALBEDO_out,            & ! [INOUT]
       HFLUX,                 & ! [OUT]
       LATENTFLUX,            & ! [OUT]
       GFLUX,                 & ! [OUT]
       QCC,                   & ! [OUT]
       QFUSION,               & ! [OUT]
       MELT,                  & ! [OUT]
       SWEMELT,               & ! [OUT]
       Gflux2land,            & ! [OUT]
       SFLX_SNOW,             & ! [IN]
       TA,                    & ! [IN]
       UA,                    & ! [IN]
       RH,                    & ! [IN]
       SW,                    & ! [IN]
       LW,                    & ! [IN]
       time                   ) ! [IN]

    implicit none
    ! prognostic variables
    real(RP), intent(inout)    :: TSNOW         ! snow temperature            [K]
    real(RP), intent(inout)    :: SWE           ! equivalent water [kg/m2]
    real(RP), intent(inout)    :: DEPTH         ! depth of melting point      [m]
    real(RP), intent(inout)    :: ZNSNOW        ! total snow depth            [m]
    ! update variables
    real(RP), intent(inout)    :: nosnowsec
    real(RP), intent(inout)    :: ALBEDO_out

    ! output variables
    real(RP), intent(out)      :: HFLUX         ! HFLUX = whole snow Sensible heat flux [W/m2]
    real(RP), intent(out)      :: LATENTFLUX    ! LATENTFLUX = whole snow Latent heat flux [W/m2]
    real(RP), intent(out)      :: GFLUX         ! GFLUX = whole snow Ground flux [W/m2]

    real(RP), intent(out)      :: QCC           ! QCC = heat used for change snow condition to isothermal [J m^-2]
    real(RP), intent(out)      :: QFUSION       ! QFUSION = heat used for change snow condition to melt point [J m^-2]
    real(RP), intent(out)      :: MELT          ! MELT = heat used for snow run off production [J m^-2 s^-1]
    real(RP), intent(out)      :: SWEMELT       ! snow water equivalent  ! [kg/m2]
    real(RP), intent(out)      :: Gflux2land    ! Residual heat, goes to land model ! [J/m2]


    ! input data
    real(RP), intent(in)       :: SFLX_SNOW
    real(RP), intent(in)       :: TA
    real(RP), intent(in)       :: UA
    real(RP), intent(in)       :: RH
    real(RP), intent(in)       :: SW
    real(RP), intent(in)       :: LW
    real(RP), intent(in)       :: time

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

    write(*,*) "UA, SNOW,SFLX_SNOW,time : ", UA, SNOW, SFLX_SNOW, time
    write(*,*) "SWE , TSNOW, and TA :     ", SWE0, TSNOW0, TA
    write(*,*) "DEPTH is:                 ", DEPTH0
    write(*,*) "ZN beginning:             ", ZNSNOW0

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
    write(*,*) "Albedo                    ",ALBEDO

!----- Energy balance at snow surface -------------------------------!

   call groundflux (TSNOW0, TA, UA, RH, ALBEDO, SW, LW, &  ! [IN]
                    GFLUX, RFLUX, SFLUX, LINFLUX, LOUTFLUX, HFLUX, LATENTFLUX) ! [OUT]


!! Check whether GFLUX (energy into snowpack) is enough to melt all snow.
!! If GFLUX is enough, the model melts all snow and then go to next timestep.

   call check_allSnowMelt   (GFLUX, TSNOW0, ZNSNOW0, DEPTH0, sflag, time)
   if(sflag .eq. 1)then
      if( IO_L ) write(IO_FID_LOG,*) '*** LAND/snow: All snow melt'
      QCC        = 0.5_RP*CSRHOS*ZNSNOW0*(T0-TSNOW0)
      QFUSION    = W0*RHOSNOW*LF*ZNSNOW0
      MELT       = (1.0_RP-W0)*RHOSNOW*LF*DEPTH0  ! [J/m2]
      TSNOW      = T0
      ZNSNOW     = 0.0_RP
      DEPTH      = 0.0_RP
      SWE        = 0.0_RP
      SWEMELT    = MELT /((1.0_RP-W0)*LF)
      Gflux2land = GFLUX*time - (QCC + QFUSION + MELT)   ! [J/m2]

   else

      call cal_param (ZNSNOW0, TSNOW0, GFLUX, TA, UA, RH, LW, time)

      ! check whether the model has solution
      call check_applicability (GFLUX, TSNOW0, ZNSNOW0, TA, UA, RH, LW, Gres, beta, time)
      if ((Gres >= 0.0_RP).and.(beta >= 0.0_RP)) then
         if( IO_L ) write(IO_FID_LOG,*) '*** LAND/snow model is not appropriate',Gres,beta
         QCC             = 0.5_RP*CSRHOS*ZNSNOW0*(T0-TSNOW0)
         QFUSION         = W0*RHOSNOW*LF*ZNSNOW0
         MELT            = Gres
         TSNOW           = T0
         ZNSNOW          = 0.0_RP
      else

         !if (t==1)then
         !   call cal_R1R2 (ZNSNOW0, TSNOW0,GFLUX,TA,UA,LW, time)
         !endif

         call snowdepth ( GFLUX, ZNSNOW0, ZNSNOW, time)

         write(*,*) "ZN is: ", ZNSNOW
         IF(ZNSNOW < ZMIN) THEN
            ZN_flag = 1
            write(*,*) "ZN is replaced to: ", ZNSNOW ," to ", ZMIN
            ZNSNOW = ZMIN
         ELSE IF(ZNSNOW > DEPTH0) THEN
            ZN_flag = 2
            write(*,*) "ZN is replaced to: ", ZNSNOW ," to ", DEPTH0
            ZNSNOW = DEPTH0
         END IF

         ! This equation is to calculate TSN
         call  equation415(LAMBDAS, C2, ZNSNOW, RH, QSAT, TSNOW0, ZNSNOW0, GFLUX, TA, UA, LW, TSNOW, time)

         write(*,*) 'TSNOW is:       ', TSNOW

         if (TSNOW > T0) then
            TS_flag = 1
            TSNOW   = T0
            call recalculateZ(ZNSNOW0, TSNOW0, GFLUX, ZNSNOW, time)
            call check_res(ZNSNOW0, ZNSNOW, TSNOW0, TSNOW, GFLUX, TA, UA, RH, LW, "1", time)
            IF (ZNSNOW < ZMIN) THEN
               ZN_flag = 4
               write(*,*) "ZN is updated/replaced to: ", ZNSNOW ," to ", ZMIN
               ZNSNOW = ZMIN
            ELSE IF(ZNSNOW > DEPTH0) THEN
               ZN_flag = 5
               write(*,*) "ZN is updated/replaced to: ", ZNSNOW ," to ", DEPTH0
               ZNSNOW = DEPTH0
            ELSE
               ZN_flag = 3
            END IF

            write(*,*) 'TSNOW  is updated:  ', TSNOW
            write(*,*) 'ZNSNOW is updated:  ', ZNSNOW
         else
            call check_res( ZNSNOW0, ZNSNOW, TSNOW0, TSNOW, GFLUX, TA,UA,RH,LW, "0", time)
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
      if( IO_L ) write(IO_FID_LOG,*) "### ZN_flag = ", ZN_flag, "TS_flag = ", TS_flag
      if( IO_L ) write(IO_FID_LOG,*) '### Heat flux from snowpack to land surface: ', Gflux2land

      DELTADEPTH             = MELT / ((1.0_RP-W0)*LF*RHOSNOW)
      SWEMELT                = DELTADEPTH*RHOSNOW
      SWE                    = SWE0        - SWEMELT
      DEPTH                  = DEPTH0      - DELTADEPTH
      if (DEPTH < ZNSNOW) then
         ! NOTICE: energy budget has not been considered thought this process yet.
         ZNSNOW               = DEPTH
      endif

      write(*,*) 'MELT in water equivalent is: ', SWEMELT
      write(*,*) 'SWE0 is:                     ', SWE
      write(*,*) 'DELTADEPTH is:               ', DELTADEPTH
      write(*,*) 'DEPTH0 is:                   ', DEPTH
      write(*,*) 'ZNSNOW0 is:                  ', ZNSNOW

   endif

   return
 end subroutine SNOW_ky90_main

!==============================================================

subroutine groundflux (TS, TA, UA, RH, ALPHA, SW, LW, &
                       GFLUX, RFLUX, SFLUX, LINFLUX, LOUTFLUX, HFLUX, LATENTFLUX)
  implicit none

  real(RP), intent(in)     :: TS, TA, UA, RH, ALPHA, SW, LW
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
  HFLUX              = -1.0_RP * (cp * rhoair * CH * UA * (TS-TA))
  LATENTFLUX         = -1.0_RP * (LV * rhoair * CE * UA * ((1-RH)*QSAT + (DELTAQSAT*(TS-TA))))

  GFLUX              = (SFLUX + RFLUX + HFLUX + LATENTFLUX)

  write(*,*) "-------------- groundflux --------------"
  write(*,*) "GFLUX is:    ", GFLUX
  write(*,*) "SFLUX:       ", SFLUX
  write(*,*) "RFLUX:       ", RFLUX, LINFLUX+LOUTFLUX
  write(*,*) " (LONG in:   ", LINFLUX,")"
  write(*,*) " (LONG out:  ", LOUTFLUX,")"
  write(*,*) "HFLUX is:    ", HFLUX
  write(*,*) "LATENT FLUX: ", LATENTFLUX

return
end subroutine groundflux

!==============================================================
subroutine check_allSnowMelt (GFLUX, TS1, ZN1, D, sflag, time)

  implicit none
  real(RP), intent(in)     :: GFLUX, TS1, ZN1, D, time
  integer, intent(out)     :: sflag
  real(RP)                 :: energy_in, energy_use
  real(RP)                 :: energy_use_ripe, energy_use_melt

  energy_in  = GFLUX * time

  energy_use_ripe = 0.5_RP*CSRHOS*ZN1*(T0-TS1) + W0*RHOSNOW*LF*ZN1
  energy_use_melt = (1.0_RP-W0)*RHOSNOW*LF*D
  energy_use      = energy_use_ripe + energy_use_melt

  if(energy_in >= energy_use)then
     sflag=1
  else
     sflag=0
  endif

  write(*,*) "Energy in  =",energy_in
  write(*,*) "Energy use =",energy_use

  return
end subroutine

!==============================================================
subroutine cal_param (ZN1, TS1, GFLUX, TA, UA, RH, LW, time)

  implicit none
  real(RP), intent(in) ::  ZN1, TS1, TA, UA, RH, LW
  real(RP), intent(in) ::  GFLUX, time

  C1      = csrhos*0.5_RP
  C2      = (4.0_RP*epsilon*sigma*(TA**3)) + (cp*rhoair*CH*UA) + (LV*rhoair*CE*UA*DELTAQSAT)
  C3      = W0*RHOSNOW*LF

  A0      = LAMBDAS*(C3*ZN1 + (C1*ZN1*(T0-TS1)- GFLUX*time))
  A1      = C2*(C3*ZN1 + C1*ZN1*(T0-TS1) - GFLUX*time) - C3*LAMBDAS
  A2      = C1*  &
            ( epsilon*(LW-(sigma*(TA**4))) + (C2*(TA-T0)) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) ) &
            - C2*C3

  write(*,*) "-------------- snowdepth --------------"
  print*, "C1",C1
  print*, "C2",C2,(4.0_RP*epsilon*sigma*(TA**3)),(cp*rhoair*CH*UA), (LV*rhoair*CE*UA*DELTAQSAT)
  print*, "C3",C3
  print*, "A0",A0
  print*, "A1",A1
  print*, "A2",A2

end subroutine

!==============================================================
subroutine check_applicability (GFLUX, TS1, ZN1, TA, UA, RH, LW, GFLUX_res, beta, time)

  implicit none
  real(RP), intent(in)     :: GFLUX, TS1, ZN1, time
  real(RP), intent(in)     :: TA, UA, RH, LW
  real(RP), intent(out)    :: GFLUX_res, beta
  real(RP)                 :: energy_in, energy_use_max

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

  real(RP), intent(in)      :: GFLUX, time
  real(RP), intent(in)      :: ZN1
  real(RP), intent(out)     :: ZN2

!calcultaing snowdepth

   !print *, -1.0_RP*A1
   !print *, -1.0_RP*((A1**2.0_RP - 4.0_RP * A2 * A0)**0.5_RP)
   !print *, 2.0_RP*A2
   !print *, A1**2.0_RP, -1.0_RP * 4.0_RP * A2 * A0

   ZN2  = ((-1.0_RP*A1) - ((A1**2.0_RP - 4.0_RP*A2*A0)**0.5_RP)) / (2.0_RP*A2)

   write(*,*) "ZN old = ",ZN1, "ZN new = ",ZN2

return
end subroutine snowdepth
!==============================================================
subroutine equation415(LAMBDAS, C2, ZN2, RH, QSAT, TS1, ZN1, GFLUX, TA, UA, LW, TS2, time)

  implicit none

  real(RP), intent(in)     :: TA, UA, LW, ZN2, GFLUX, TS1,ZN1
  real(RP), intent(in)     :: C2, LAMBDAS, RH, QSAT, time
  real(RP), intent(out)    :: TS2
  real(RP)                 :: TS_check

  TS2      = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
           + ZN2*epsilon*(LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))

  write(*,*) "-------------- equation415 --------------"

  TS_check = GFLUX*time/(0.5*CSRHOS*ZN2) + T0 - ZN1*(T0-TS1)/ZN2 - W0*RHOSNOW*LF*(ZN1-ZN2)/(0.5*CSRHOS*ZN2)

  write(*,*) "compare ",TS2, TS_check, TS2-TS_check
  ! When ZN2 is replaced, TS2 does not equal to TS_check.

 return
end subroutine equation415
!==============================================================
subroutine recalculateZ(ZN1, TS, GFLUX, ZN2, time)

  implicit none

  real(RP), intent(in)     :: ZN1
  real(RP), intent(in)     :: TS, GFLUX, time
  real(RP), intent(out)    :: ZN2

  ZN2 = ZN1 + ((C1*ZN1*(T0-TS) - GFLUX*time) / C3)

  write(*,*) "-------------- recalculate Z --------------"

 return
end subroutine recalculateZ

!==============================================================
subroutine calculationMO(GFLUX, CSRHOS, ZN1, TS1, ZN2, TS2, &
                         MELT, QCC, QFUSION, time)
  implicit none

  real(RP), intent(in)      :: GFLUX, CSRHOS, ZN1, TS1, TS2, ZN2, time
  real(RP), intent(out)     :: MELT, QCC, QFUSION

  QCC             = 0.5_RP*CSRHOS*(ZN1*(T0-TS1)-ZN2*(T0-TS2))
  QFUSION         = W0*RHOSNOW*LF*(ZN1-ZN2)
  MELT            = ( GFLUX*time - QCC - QFUSION )

  write(*,*) "--------------------------MELT----------------"
  write(*,*) "GFLUX*time is: ", GFLUX*time
  write(*,*) "QCC is       : ", QCC
  write(*,*) "QFUSION is   : ", QFUSION
  write(*,*) "QMELT is     : ", MELT

  write(*,*) QCC+QFUSION+MELT
  write(*,*) "diff= ", QCC + QFUSION + MELT - (GFLUX*time)
  if ( ABS(QCC+QFUSION+MELT - (GFLUX*time)) > 10.) then
    print *, "Calculation is fault. Model would include bugs. Please check! Melt"
    stop
  endif


 return
end subroutine calculationMO

!==============================================================
subroutine calculationNoMO(GFLUX, CSRHOS, ZN1, TS1, ZN2, TS2, &
                           MELT, QCC, QFUSION, time)
 implicit none

 real(RP), intent(in)        :: GFLUX, CSRHOS, ZN1, TS1, TS2, ZN2, time
 real(RP), intent(out)       :: MELT, QCC, QFUSION

  QCC             = 0.5_RP*CSRHOS*(ZN1*(T0-TS1)-ZN2*(T0-TS2))
  QFUSION         = W0*RHOSNOW*LF*(ZN1-ZN2)
  MELT            = 0.0_RP

  write(*,*) "--------------------------NOMELT----------------"
  write(*,*) "GFLUX*time is: ", GFLUX*time
  write(*,*) "QCC is       : ", QCC
  write(*,*) "QFUSION is   : ", QFUSION
  write(*,*) "QMELT is     : ", MELT

  write(*,*) QCC+QFUSION+MELT
  write(*,*) "diff= ", QCC +QFUSION - (GFLUX*time)
  !if ( ABS(QCC+QFUSION - (GFLUX*time)) > 10.) then
  !  print *, "Calculation is fault. Model would include bugs. Please check! No Melt"
  !  stop
  !endif

 return
end subroutine calculationNoMO

!==============================================================
subroutine check_res(ZN1, ZN2, TS1, TS2, GFLUX, TA, UA, RH, LW, flag, time)

  implicit none
  real(RP) , intent(in)  :: ZN1, ZN2, TS1, TS2, GFLUX, time
  real(RP) , intent(in)  :: TA, UA, RH, LW
  character(len=*)       :: flag
  real(RP)               :: R1 ! Eq. (2)
  real(RP)               :: R2 ! Eq. (8)
  real(RP)               :: R3 ! Eq. (8)

  R3 = -999.
  R1 = C1 * (ZN1*(T0-TS1) - ZN2*(T0-TS2)) + C3*(ZN1-ZN2)-GFLUX*time
  R2 = ZN2*( epsilon*LW-epsilon*sigma*(TA**4) - (C2*(TS2-TA)) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) ) &
       + LAMBDAS*(T0-TS2)
  R3   = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
           + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))-TS2


  write(*,*) "R1 is         : ", R1, "flag = ", flag
  write(*,*) "R2 is         : ", R2, R3,"flag = ", flag


  if(abs(R1)>10000)then
     write(*,*) C1 * (ZN1*(T0-TS1) - ZN2*(T0-TS2)), C3*(ZN1-ZN2), -GFLUX*time
  endif

  return
end subroutine check_res

!==============================================================
subroutine cal_R1R2(ZN1, TS1,GFLUX,TA,UA,RH,LW, time)

  implicit none
  real(RP), intent(in)   :: ZN1, TS1
  real(RP), intent(in)   :: GFLUX,TA,UA,RH,LW, time
  real(RP)               :: ZN2, TS2
  real(RP)               :: R1 ! Eq. (2)
  real(RP)               :: R2 ! Eq. (8)
  real(RP)               :: R3 ! Eq. (8)

  real(RP)               :: ts0,zn0
  integer                :: i,j
  real(RP)               :: ts_r1,ts_r2,zn_r1,zn_r2
  character(len=3)       :: ttt

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
     write(70,'(3f15.5)') ZN2, ts_r1, ts_r2
  enddo
  do j=1,100000
     ZN2 = -1.0_RP + 2.0_RP*0.00001_RP*real((j-1),kind=RP)

     if (ZN2/=0.0)then
        ts_r1 = T0 - (C1*ZN1*(T0-TS1) + C3*(ZN1-ZN2) - GFLUX*time)/(C1*ZN2)

        ts_r2 = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
             + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))
     endif
     write(70,'(3f15.5)') ZN2, ts_r1, ts_r2
  enddo
  do j=2,10000
     ZN2 = 1.0_RP + 9.0_RP*0.0001_RP*real((j-1),kind=RP)

     if (ZN2/=0.0)then
        ts_r1 = T0 - (C1*ZN1*(T0-TS1) + C3*(ZN1-ZN2) - GFLUX*time)/(C1*ZN2)

        ts_r2 = ((LAMBDAS*T0) + (TA*C2*ZN2) - (ZN2*LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) &
             + ZN2*epsilon*( LW - (sigma*(TA**4))))/ (LAMBDAS + (C2*ZN2))
     endif
     write(70,'(3f15.5)') ZN2, ts_r1, ts_r2
  enddo
  close(70)

  a=T0+C3/C1
  b=-1.0_RP * (C1*ZN1*(T0-TS1) + C3*ZN1 - GFLUX*time)
  c=C1

  d=(LAMBDAS*T0)
  e=( TA*C2 - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT) + epsilon*( LW - (sigma*(TA**4)))) ! *ZN2
  f=(LAMBDAS + (C2*ZN2))


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

  !print  *,"aa",(LINFLUX-epsilon*sigma*(TA**4) + (C2*TA) - (LV*RHOAIR*CE*UA*(1.0_RP-RH)*QSAT))/C2
  return
end subroutine cal_R1R2

!!!!!!!!!!!!!!!!!!!!!!SUBROUTINE END!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module scale_land_phy_snow_KY90
