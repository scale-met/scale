!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new] Imported from SCALE-RM ver.2
!! @li      2012-01-31 (Y.Miyamoto) [add] Lamb wave test
!! @li      2012-01-31 (Y.Miyamoto) [add] KH wave test
!! @li      2012-01-31 (Y.Miyamoto) [add] turbulence test
!! @li      2011-02-01 (H.Yashiro)  [add] supercell test, follow the supercell test of WRF
!! @li      2012-02-06 (Y.Miyamoto) [add] advection test
!! @li      2012-02-16 (Y.Miyamoto) [mod] added hydrostatic balance calculation
!! @li      2012-03-27 (H.Yashiro)  [mod] change subroutines into one module
!! @li      2012-04-04 (Y.Miyamoto) [new] SQUALLINE test, for GCSS model comparison (Redelsperger et al. 2000)
!! @li      2012-04-06 (H.Yashiro)  [new] uniform state test
!! @li      2012-04-08 (H.Yashiro)  [mod] merge all init programs
!! @li      2012-06-13 (Y.Sato)     [mod] add hbinw option (***HBINW)
!! @li      2013-02-25 (H.Yashiro)  [mod] ISA profile
!! @li      2014-03-27 (A.Noda)     [mod] add DYCOMS2_RF02_DNS
!! @li      2014-04-28 (R.Yoshida)  [add] real case experiment
!! @li      2014-08-26 (A.Noda)     [mod] add GRAYZONE
!! @li      2015-03-27 (Y.Sato)     [mod] add Box aero
!! @li      2015-04-30 (Y.Sato)     [mod] add WARMBUBBLE-AERO
!!
!<
!-------------------------------------------------------------------------------
module mod_mkinit
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_MPIstop
  use scale_const, only: &
     PI    => CONST_PI,    &
     GRAV  => CONST_GRAV,  &
     Pstd  => CONST_Pstd,  &
     Rdry  => CONST_Rdry,  &
     CPdry => CONST_CPdry, &
     LHV   => CONST_LHV,   &
     P00   => CONST_PRE00, &
     I_SW  => CONST_I_SW,  &
     I_LW  => CONST_I_LW
  use scale_random, only: &
     RANDOM_get
  use scale_comm, only: &
     COMM_vars8, &
     COMM_wait
  use scale_grid, only: &
     GRID_CZ,  &
     GRID_CX,  &
     GRID_CY,  &
     GRID_FZ,  &
     GRID_FX,  &
     GRID_FY,  &
     GRID_CXG
  use scale_grid_real, only: &
     REAL_CZ, &
     REAL_FZ
  use scale_atmos_profile, only: &
     PROFILE_isa => ATMOS_PROFILE_isa
  use scale_atmos_hydrostatic, only: &
     HYDROSTATIC_buildrho        => ATMOS_HYDROSTATIC_buildrho,       &
     HYDROSTATIC_buildrho_atmos  => ATMOS_HYDROSTATIC_buildrho_atmos, &
     HYDROSTATIC_buildrho_bytemp => ATMOS_HYDROSTATIC_buildrho_bytemp
  use scale_atmos_saturation, only: &
     SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all, &
     SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
  use mod_atmos_vars, only: &
     DENS, &
     MOMX, &
     MOMY, &
     MOMZ, &
     RHOT, &
     QTRC
  use mod_atmos_phy_ae_vars, only: &
     CCN => ATMOS_PHY_AE_CCN
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKINIT_setup
  public :: MKINIT

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public            :: MKINIT_TYPE        = -1
  integer, public, parameter :: I_IGNORE           =  0

  integer, public, parameter :: I_PLANESTATE       =  1
  integer, public, parameter :: I_TRACERBUBBLE     =  2
  integer, public, parameter :: I_COLDBUBBLE       =  3

  integer, public, parameter :: I_LAMBWAVE         =  4
  integer, public, parameter :: I_GRAVITYWAVE      =  5
  integer, public, parameter :: I_KHWAVE           =  6
  integer, public, parameter :: I_TURBULENCE       =  7
  integer, public, parameter :: I_MOUNTAINWAVE     =  8

  integer, public, parameter :: I_WARMBUBBLE       =  9
  integer, public, parameter :: I_SUPERCELL        = 10
  integer, public, parameter :: I_SQUALLLINE       = 11
  integer, public, parameter :: I_WK1982           = 12
  integer, public, parameter :: I_DYCOMS2_RF01     = 13
  integer, public, parameter :: I_DYCOMS2_RF02     = 14
  integer, public, parameter :: I_RICO             = 15

  integer, public, parameter :: I_INTERPORATION    = 16

  integer, public, parameter :: I_LANDCOUPLE       = 17
  integer, public, parameter :: I_OCEANCOUPLE      = 18
  integer, public, parameter :: I_URBANCOUPLE      = 19
  integer, public, parameter :: I_TRIPLECOUPLE     = 20
  integer, public, parameter :: I_BUBBLECOUPLE     = 21

  integer, public, parameter :: I_SEABREEZE        = 22
  integer, public, parameter :: I_HEATISLAND       = 23

  integer, public, parameter :: I_DYCOMS2_RF02_DNS = 24

  integer, public, parameter :: I_REAL             = 25

  integer, public, parameter :: I_GRAYZONE         = 26
  integer, public, parameter :: I_BOXAERO          = 27
  integer, public, parameter :: I_WARMBUBBLEAERO   = 28

  integer, public, parameter :: I_CAVITYFLOW       = 29

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: BUBBLE_setup
  private :: SBMAERO_setup
  private :: AEROSOL_setup

  private :: MKINIT_planestate
  private :: MKINIT_tracerbubble
  private :: MKINIT_coldbubble
  private :: MKINIT_lambwave
  private :: MKINIT_gravitywave
  private :: MKINIT_khwave
  private :: MKINIT_turbulence
  private :: MKINIT_cavityflow
  private :: MKINIT_mountainwave

  private :: MKINIT_warmbubble
  private :: MKINIT_supercell
  private :: MKINIT_squallline
  private :: MKINIT_wk1982
  private :: MKINIT_DYCOMS2_RF01
  private :: MKINIT_DYCOMS2_RF02
  private :: MKINIT_RICO

  private :: MKINIT_interporation

  private :: MKINIT_landcouple
  private :: MKINIT_oceancouple
  private :: MKINIT_urbancouple
  private :: MKINIT_seabreeze
  private :: MKINIT_heatisland

  private :: MKINIT_DYCOMS2_RF02_DNS

  private :: MKINIT_real

  private :: MKINIT_grayzone

  private :: MKINIT_boxaero
  private :: MKINIT_warmbubbleaero

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter           :: THETAstd = 300.0_RP ! [K]

  real(RP), private, allocatable         :: pres    (:,:,:) ! pressure [Pa]
  real(RP), private, allocatable         :: temp    (:,:,:) ! temperature [K]
  real(RP), private, allocatable         :: pott    (:,:,:) ! potential temperature [K]
  real(RP), private, allocatable         :: qsat    (:,:,:) ! satulated water vapor [kg/kg]
  real(RP), private, allocatable         :: qv      (:,:,:) ! water vapor [kg/kg]
  real(RP), private, allocatable         :: qc      (:,:,:) ! cloud water [kg/kg]
  real(RP), private, allocatable         :: velx    (:,:,:) ! velocity u [m/s]
  real(RP), private, allocatable         :: vely    (:,:,:) ! velocity v [m/s]

  real(RP), private, allocatable         :: pres_sfc(:,:,:) ! surface pressure [Pa]
  real(RP), private, allocatable         :: temp_sfc(:,:,:) ! surface temperature [K]
  real(RP), private, allocatable         :: pott_sfc(:,:,:) ! surface potential temperature [K]
  real(RP), private, allocatable         :: qsat_sfc(:,:,:) ! surface satulated water vapor [kg/kg]
  real(RP), private, allocatable         :: qv_sfc  (:,:,:) ! surface water vapor [kg/kg]
  real(RP), private, allocatable         :: qc_sfc  (:,:,:) ! surface cloud water [kg/kg]

  real(RP), private, allocatable         :: rndm    (:,:,:) ! random    number (0-1)
  real(RP), private, allocatable, target :: bubble  (:,:,:) ! bubble    factor (0-1)
  real(RP), private, allocatable, target :: rect    (:,:,:) ! rectangle factor (0-1)
  real(RP), private, allocatable         :: gan     (:)     ! gamma     factor (0-1)

  logical,  private                      :: flg_intrp = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKINIT_setup
    implicit none

    character(len=H_SHORT) :: MKINIT_initname = 'NONE'

    NAMELIST / PARAM_MKINIT / &
       MKINIT_initname, &
       flg_intrp

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[make init] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT)

    allocate( pres(KA,IA,JA) )
    allocate( temp(KA,IA,JA) )
    allocate( pott(KA,IA,JA) )
    allocate( qsat(KA,IA,JA) )
    allocate( qv  (KA,IA,JA) )
    allocate( qc  (KA,IA,JA) )
    allocate( velx(KA,IA,JA) )
    allocate( vely(KA,IA,JA) )

    allocate( pres_sfc(1,IA,JA) )
    allocate( temp_sfc(1,IA,JA) )
    allocate( pott_sfc(1,IA,JA) )
    allocate( qsat_sfc(1,IA,JA) )
    allocate( qv_sfc  (1,IA,JA) )
    allocate( qc_sfc  (1,IA,JA) )

    allocate( rndm  (KA,IA,JA) )
    allocate( bubble(KA,IA,JA) )
    allocate( rect  (KA,IA,JA) )

    select case(trim(MKINIT_initname))
    case('NONE')
       MKINIT_TYPE = I_IGNORE
    case('PLANESTATE')
       MKINIT_TYPE = I_PLANESTATE
    case('TRACERBUBBLE')
       MKINIT_TYPE = I_TRACERBUBBLE
    case('COLDBUBBLE')
       MKINIT_TYPE = I_COLDBUBBLE
       call BUBBLE_setup
    case('LAMBWAVE')
       MKINIT_TYPE = I_LAMBWAVE
       call BUBBLE_setup
    case('GRAVITYWAVE')
       MKINIT_TYPE = I_GRAVITYWAVE
       call BUBBLE_setup
    case('KHWAVE')
       MKINIT_TYPE = I_KHWAVE
    case('TURBULENCE')
       MKINIT_TYPE = I_TURBULENCE
    case('MOUNTAINWAVE')
       MKINIT_TYPE = I_MOUNTAINWAVE
       call BUBBLE_setup
    case('WARMBUBBLE')
       MKINIT_TYPE = I_WARMBUBBLE
       call BUBBLE_setup
    case('SUPERCELL')
       MKINIT_TYPE = I_SUPERCELL
       call BUBBLE_setup
    case('SQUALLLINE')
       MKINIT_TYPE = I_SQUALLLINE
    case('WK1982')
       MKINIT_TYPE = I_WK1982
       call BUBBLE_setup
    case('DYCOMS2_RF01')
       MKINIT_TYPE = I_DYCOMS2_RF01
    case('DYCOMS2_RF02')
       MKINIT_TYPE = I_DYCOMS2_RF02
    case('RICO')
       MKINIT_TYPE = I_RICO
    case('INTERPORATION')
       MKINIT_TYPE = I_INTERPORATION
    case('LANDCOUPLE')
       MKINIT_TYPE = I_LANDCOUPLE
    case('OCEANCOUPLE')
       MKINIT_TYPE = I_OCEANCOUPLE
    case('URBANCOUPLE')
       MKINIT_TYPE = I_URBANCOUPLE
    case('TRIPLECOUPLE')
       MKINIT_TYPE = I_TRIPLECOUPLE
    case('BUBBLECOUPLE')
       MKINIT_TYPE = I_BUBBLECOUPLE
       call BUBBLE_setup
    case('SEABREEZE')
       MKINIT_TYPE = I_SEABREEZE
    case('HEATISLAND')
       MKINIT_TYPE = I_HEATISLAND
    case('DYCOMS2_RF02_DNS')
       MKINIT_TYPE = I_DYCOMS2_RF02_DNS
    case('REAL')
       MKINIT_TYPE = I_REAL
    case('GRAYZONE')
       MKINIT_TYPE = I_GRAYZONE
    case('BOXAERO')
       MKINIT_TYPE = I_BOXAERO
    case('WARMBUBBLEAERO')
       MKINIT_TYPE = I_WARMBUBBLEAERO
       call BUBBLE_setup
    case('CAVITYFLOW')
       MKINIT_TYPE = I_CAVITYFLOW
    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(MKINIT_initname)
       call PRC_MPIstop
    endselect

    return
  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKINIT
    use scale_const, only: &
       CONST_UNDEF8
    use scale_landuse, only: &
       LANDUSE_write
    use mod_atmos_driver, only: &
       ATMOS_SURFACE_GET
    use mod_atmos_vars, only: &
       ATMOS_sw_restart => ATMOS_RESTART_OUTPUT, &
       ATMOS_vars_restart_write
    use mod_ocean_driver, only: &
       OCEAN_SURFACE_SET
    use mod_ocean_vars, only: &
       OCEAN_sw_restart => OCEAN_RESTART_OUTPUT, &
       OCEAN_vars_restart_write
    use mod_land_driver, only: &
       LAND_SURFACE_SET
    use mod_land_vars, only: &
       LAND_sw_restart => LAND_RESTART_OUTPUT, &
       LAND_vars_restart_write
    use mod_urban_driver, only: &
       URBAN_SURFACE_SET
    use mod_urban_vars, only: &
       URBAN_sw_restart => URBAN_RESTART_OUTPUT, &
       URBAN_vars_restart_write
    implicit none

    integer :: iq
    !---------------------------------------------------------------------------

    if ( MKINIT_TYPE == I_IGNORE ) then
      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  MAKING INITIAL DATA ++++++'
    else
      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'

      !--- Initialize variables
#ifndef DRY
      do iq = 2, QA
         QTRC(:,:,:,iq) = 0.0_RP
      enddo
#endif

      pres(:,:,:) = CONST_UNDEF8
      temp(:,:,:) = CONST_UNDEF8
      pott(:,:,:) = CONST_UNDEF8
      qsat(:,:,:) = CONST_UNDEF8
      qv  (:,:,:) = CONST_UNDEF8
      qc  (:,:,:) = CONST_UNDEF8
      velx(:,:,:) = CONST_UNDEF8
      vely(:,:,:) = CONST_UNDEF8

      rndm  (:,:,:) = CONST_UNDEF8

      pres_sfc(:,:,:) = CONST_UNDEF8
      temp_sfc(:,:,:) = CONST_UNDEF8
      pott_sfc(:,:,:) = CONST_UNDEF8
      qsat_sfc(:,:,:) = CONST_UNDEF8
      qv_sfc  (:,:,:) = CONST_UNDEF8
      qc_sfc  (:,:,:) = CONST_UNDEF8

      if ( TRACER_TYPE == 'SUZUKI10' ) then
         if( IO_L ) write(IO_FID_LOG,*) '*** Aerosols for SBM are included ***'
         call SBMAERO_setup
      endif

      select case(MKINIT_TYPE)
      case(I_PLANESTATE)
         call MKINIT_planestate
      case(I_TRACERBUBBLE)
         call MKINIT_tracerbubble
      case(I_COLDBUBBLE)
         call MKINIT_coldbubble
      case(I_LAMBWAVE)
         call MKINIT_lambwave
      case(I_GRAVITYWAVE)
         call MKINIT_gravitywave
      case(I_KHWAVE)
         call MKINIT_khwave
      case(I_TURBULENCE)
         call MKINIT_turbulence
      case(I_MOUNTAINWAVE)
         call MKINIT_mountainwave
      case(I_WARMBUBBLE)
         call MKINIT_warmbubble
      case(I_SUPERCELL)
         call MKINIT_supercell
      case(I_SQUALLLINE)
         call MKINIT_squallline
      case(I_WK1982)
         call MKINIT_wk1982
      case(I_DYCOMS2_RF01)
         call MKINIT_DYCOMS2_RF01
      case(I_DYCOMS2_RF02)
         call MKINIT_DYCOMS2_RF02
      case(I_RICO)
         call MKINIT_RICO
      case(I_INTERPORATION)
         call MKINIT_INTERPORATION
      case(I_OCEANCOUPLE)
         call MKINIT_planestate
         call MKINIT_oceancouple
      case(I_LANDCOUPLE)
         call MKINIT_planestate
         call MKINIT_landcouple
      case(I_URBANCOUPLE)
         call MKINIT_planestate
         call MKINIT_landcouple ! tentative
         call MKINIT_urbancouple
      case(I_TRIPLECOUPLE)
         call MKINIT_planestate
         call MKINIT_oceancouple
         call MKINIT_landcouple
         call MKINIT_urbancouple
      case(I_BUBBLECOUPLE)
         call MKINIT_planestate
         call MKINIT_warmbubble
         call MKINIT_oceancouple
         call MKINIT_landcouple
         call MKINIT_urbancouple
      case(I_SEABREEZE)
         call MKINIT_planestate
         call MKINIT_seabreeze
      case(I_HEATISLAND)
         call MKINIT_planestate
         call MKINIT_heatisland
      case(I_DYCOMS2_RF02_DNS)
         call MKINIT_DYCOMS2_RF02_DNS
      case(I_REAL)
         call MKINIT_real
      case(I_GRAYZONE)
         call MKINIT_grayzone
      case(I_BOXAERO)
         call MKINIT_boxaero
      case(I_WARMBUBBLEAERO)
         call MKINIT_warmbubbleaero
      case(I_CAVITYFLOW)
         call MKINIT_cavityflow
      case default
         write(*,*) ' xxx Unsupported TYPE:', MKINIT_TYPE
         call PRC_MPIstop
      endselect

      if( IO_L ) write(IO_FID_LOG,*) '++++++ END   MAKING INITIAL  DATA ++++++'

      ! setup surface condition
      call OCEAN_SURFACE_SET( countup = .false. )
      call LAND_SURFACE_SET ( countup = .false. )
      call URBAN_SURFACE_SET( countup = .false. )
      call ATMOS_SURFACE_GET

      ! output boundary file
      call LANDUSE_write

      ! output restart file
      if( ATMOS_sw_restart ) call ATMOS_vars_restart_write
      if( OCEAN_sw_restart ) call OCEAN_vars_restart_write
      if( LAND_sw_restart  ) call LAND_vars_restart_write
      if( URBAN_sw_restart ) call URBAN_vars_restart_write
    endif

    return
  end subroutine MKINIT

  !-----------------------------------------------------------------------------
  !> Bubble
  subroutine BUBBLE_setup
    use scale_const, only: &
       CONST_UNDEF8
    implicit none

    ! Bubble
    logical  :: BBL_eachnode = .false.  ! Arrange bubble at each node? [kg/kg]
    real(RP) :: BBL_CZ       =  2.E3_RP ! center location [m]: z
    real(RP) :: BBL_CX       =  2.E3_RP ! center location [m]: x
    real(RP) :: BBL_CY       =  2.E3_RP ! center location [m]: y
    real(RP) :: BBL_RZ       =  2.E3_RP ! bubble radius   [m]: z
    real(RP) :: BBL_RX       =  2.E3_RP ! bubble radius   [m]: x
    real(RP) :: BBL_RY       =  2.E3_RP ! bubble radius   [m]: y

    NAMELIST / PARAM_BUBBLE / &
       BBL_eachnode, &
       BBL_CZ,       &
       BBL_CX,       &
       BBL_CY,       &
       BBL_RZ,       &
       BBL_RX,       &
       BBL_RY

    real(RP) :: CZ_offset
    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP) :: dist

    integer  :: ierr
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit bubble] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found namelist. Check!'
       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_BUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_BUBBLE)

    bubble(:,:,:) = CONST_UNDEF8

    if ( BBL_eachnode ) then
       CZ_offset = GRID_CZ(KS)
       CX_offset = GRID_CX(IS)
       CY_offset = GRID_CY(JS)
    else
       CZ_offset = 0.0_RP
       CX_offset = 0.0_RP
       CY_offset = 0.0_RP
    endif

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE

       ! make tracer bubble
       dist = ( (GRID_CZ(k)-CZ_offset-BBL_CZ)/BBL_RZ )**2 &
            + ( (GRID_CX(i)-CX_offset-BBL_CX)/BBL_RX )**2 &
            + ( (GRID_CY(j)-CY_offset-BBL_CY)/BBL_RY )**2

       bubble(k,i,j) = cos( 0.5_RP*PI*sqrt( min(dist,1.0_RP) ) )**2

    enddo
    enddo
    enddo

    return
  end subroutine BUBBLE_setup

  !-----------------------------------------------------------------------------
  !> Bubble
  subroutine RECT_setup
    use scale_const, only: &
       CONST_UNDEF8
    implicit none

    ! Bubble
    logical  :: RCT_eachnode = .false.  ! Arrange rectangle at each node? [kg/kg]
    real(RP) :: RCT_CZ       =  2.E3_RP ! center location [m]: z
    real(RP) :: RCT_CX       =  2.E3_RP ! center location [m]: x
    real(RP) :: RCT_CY       =  2.E3_RP ! center location [m]: y
    real(RP) :: RCT_RZ       =  2.E3_RP ! rectangle z width   [m]: z
    real(RP) :: RCT_RX       =  2.E3_RP ! rectangle x width   [m]: x
    real(RP) :: RCT_RY       =  2.E3_RP ! rectangle y width   [m]: y

    NAMELIST / PARAM_RECT / &
       RCT_eachnode, &
       RCT_CZ,       &
       RCT_CX,       &
       RCT_CY,       &
       RCT_RZ,       &
       RCT_RX,       &
       RCT_RY

    real(RP) :: CZ_offset
    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP) :: dist

    integer  :: ierr
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit rectangle] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RECT,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found namelist. Check!'
       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_RECT. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_RECT)

    rect(:,:,:) = CONST_UNDEF8

    if ( RCT_eachnode ) then
       CZ_offset = GRID_CZ(KS)
       CX_offset = GRID_CX(IS)
       CY_offset = GRID_CY(JS)
    else
       CZ_offset = 0.0_RP
       CX_offset = 0.0_RP
       CY_offset = 0.0_RP
    endif

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE

       ! make tracer rectangle
       dist = 2.0_RP * max( &
            abs(GRID_CZ(k) - CZ_offset - RCT_CZ)/RCT_RZ,   &
            abs(GRID_CX(i) - CX_offset - RCT_CX)/RCT_RX,   &
            abs(GRID_CY(j) - CY_offset - RCT_CY)/RCT_RY    &
            & )
       if ( dist <= 1.0_RP ) then
         rect(k,i,j) = 1.0_RP
       else
         rect(k,i,j) = 0.0_RP
       end if
    enddo
    enddo
    enddo

    return
  end subroutine RECT_setup

  !-----------------------------------------------------------------------------
  !> Setup aerosol condition for Kajino 2013 scheme
  subroutine AEROSOL_setup

    use scale_tracer
    use scale_const, only: &
       PI => CONST_PI, &
       CONST_CVdry, &
       CONST_Rdry, &
       CONST_Rvap, &
       CONST_PRE00
    use scale_atmos_thermodyn, only: &
       AQ_CV
    implicit none

    integer, parameter ::  ia_m0  = 1             !1. number conc        [#/m3]
    integer, parameter ::  ia_m2  = 2             !2. 2nd mom conc       [m2/m3]
    integer, parameter ::  ia_m3  = 3             !3. 3rd mom conc       [m3/m3]
    integer, parameter ::  ia_ms  = 4             !4. mass conc          [ug/m3]
    integer, parameter ::  ia_kp  = 5
    integer, parameter ::  ik_out  = 1

    real(RP) :: m0_init = 0.0_RP    ! initial total num. conc. of modes (Atk,Acm,Cor) [#/m3]
    real(RP) :: dg_init = 80.e-9_RP ! initial number equivalen diameters of modes     [m]
    real(RP) :: sg_init = 1.6_RP    ! initial standard deviation                      [-]

    real(RP) :: d_min_inp(3)
    real(RP) :: d_max_inp(3)
    real(RP) :: k_min_inp(3)
    real(RP) :: k_max_inp(3)
    integer  :: n_kap_inp(3)

    real(RP), parameter :: d_min_def = 1.e-9_RP ! default lower bound of 1st size bin
    real(RP), parameter :: d_max_def = 1.e-5_RP ! upper bound of last size bin
    integer,  parameter :: n_kap_def = 1        ! number of kappa bins
    real(RP), parameter :: k_min_def = 0.e0_RP  ! lower bound of 1st kappa bin
    real(RP), parameter :: k_max_def = 1.e0_RP  ! upper bound of last kappa bin
    real(RP) :: c_kappa = 0.3_RP     ! hygroscopicity of condensable mass              [-]
    real(RP), parameter :: cleannumber = 1.e-3_RP

    NAMELIST / PARAM_AERO / &
       m0_init, &
       dg_init, &
       sg_init, &
       d_min_inp, &
       d_max_inp, &
       k_min_inp, &
       k_max_inp, &
       n_kap_inp

    ! local variables
    real(RP), parameter  :: rhod_ae    = 1.83_RP              ! particle density [g/cm3] sulfate assumed
    real(RP), parameter  :: conv_vl_ms = rhod_ae/1.e-12_RP     ! M3(volume)[m3/m3] to mass[m3/m3]

    integer :: ia0, ik, is0, ic, k, i, j, it, iq
    integer :: ierr, n_trans
    real(RP) :: m0t, dgt, sgt, m2t, m3t, mst
    real(RP),allocatable :: aerosol_procs(:,:,:,:) !(n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP),allocatable :: aerosol_activ(:,:,:,:) !(n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP),allocatable :: emis_procs(:,:,:,:)    !(n_atr,n_siz_max,n_kap_max,n_ctg)
    !--- gas
    real(RP) :: conc_gas(GAS_CTG)    !concentration [ug/m3]
    integer :: n_siz_max, n_kap_max, n_ctg
!   integer, allocatable :: it_procs2trans(:,:,:,:) !procs to trans conversion
!   integer, allocatable :: ia_trans2procs(:) !trans to procs conversion
!   integer, allocatable :: is_trans2procs(:) !trans to procs conversion
!   integer, allocatable :: ik_trans2procs(:) !trans to procs conversion
!   integer, allocatable :: ic_trans2procs(:)
    !--- bin settings (lower bound, center, upper bound)
    real(RP),allocatable :: d_lw(:,:), d_ct(:,:), d_up(:,:)  !diameter [m]
    real(RP),allocatable :: k_lw(:,:), k_ct(:,:), k_up(:,:)  !kappa    [-]
    real(RP) :: dlogd, dk                                    !delta log(D), delta K
    real(RP),allocatable :: d_min(:)              !lower bound of 1st size bin   (n_ctg)
    real(RP),allocatable :: d_max(:)              !upper bound of last size bin  (n_ctg)
    integer, allocatable :: n_kap(:)              !number of kappa bins          (n_ctg)
    real(RP),allocatable :: k_min(:)              !lower bound of 1st kappa bin  (n_ctg)
    real(RP),allocatable :: k_max(:)

    real(RP) :: pott, qdry, pres
    real(RP) :: temp, cpa, cva, qsat_tmp, ssliq, Rmoist
    real(RP) :: pi6

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit aerosol] / Categ[preprocess] / Origin[SCALE-RM]'

    pi6   = pi / 6._RP
    n_ctg = AE_CTG

    allocate( d_min(n_ctg) )
    allocate( d_max(n_ctg) )
    allocate( n_kap(n_ctg) )
    allocate( k_min(n_ctg) )
    allocate( k_max(n_ctg) )

    d_min(1:n_ctg) = d_min_def ! lower bound of 1st size bin
    d_max(1:n_ctg) = d_max_def ! upper bound of last size bin
    n_kap(1:n_ctg) = n_kap_def ! number of kappa bins
    k_min(1:n_ctg) = k_min_def ! lower bound of 1st kappa bin
    k_max(1:n_ctg) = k_max_def
    d_min_inp(1:n_ctg) = d_min(1:n_ctg) !def value used if not defined in namelist
    d_max_inp(1:n_ctg) = d_max(1:n_ctg) !def value used if not defined in namelist
    n_kap_inp(1:n_ctg) = n_kap(1:n_ctg) !def value used if not defined in namelist
    k_min_inp(1:n_ctg) = k_min(1:n_ctg) !def value used if not defined in namelist
    k_max_inp(1:n_ctg) = k_max(1:n_ctg) !def value used if not defined in namelist

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_AERO,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found namelist. Default used!'
!       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_AERO. Check!'
       call PRC_MPIstop
    else
       d_min(1:n_ctg) = d_min_inp(1:n_ctg)  ! lower bound of 1st size bin
       d_max(1:n_ctg) = d_max_inp(1:n_ctg)  ! upper bound of last size bin
       n_kap(1:n_ctg) = n_kap_inp(1:n_ctg)  ! number of kappa bins
       k_min(1:n_ctg) = k_min_inp(1:n_ctg)  ! lower bound of 1st kappa bin
       k_max(1:n_ctg) = k_max_inp(1:n_ctg)  ! upper bound of last kappa bin
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_AERO)

    !--- diagnose parameters (n_trans, n_siz_max, n_kap_max)
    n_trans   = 0
    n_siz_max = 0
    n_kap_max = 0
    do ic = 1, n_ctg
      n_trans   = n_trans + NSIZ(ic) * NKAP(ic) * N_ATR
      n_siz_max = max(n_siz_max, NSIZ(ic))
      n_kap_max = max(n_kap_max, NKAP(ic))
    enddo

    allocate( aerosol_procs (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( aerosol_activ (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( emis_procs    (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    aerosol_procs(:,:,:,:) = 0._RP
    aerosol_activ(:,:,:,:) = 0._RP
    emis_procs   (:,:,:,:) = 0._RP

!   allocate( it_procs2trans(N_ATR,n_siz_max,n_kap_max,n_ctg)  )
!   allocate( ia_trans2procs(n_trans) )
!   allocate( is_trans2procs(n_trans) )
!   allocate( ik_trans2procs(n_trans) )
!   allocate( ic_trans2procs(n_trans) )

    !bin setting
    allocate(d_lw(n_siz_max,n_ctg))
    allocate(d_ct(n_siz_max,n_ctg))
    allocate(d_up(n_siz_max,n_ctg))
    allocate(k_lw(n_kap_max,n_ctg))
    allocate(k_ct(n_kap_max,n_ctg))
    allocate(k_up(n_kap_max,n_ctg))
    d_lw(:,:) = 0._RP
    d_ct(:,:) = 0._RP
    d_up(:,:) = 0._RP
    k_lw(:,:) = 0._RP
    k_ct(:,:) = 0._RP
    k_up(:,:) = 0._RP

    do ic = 1, n_ctg

      dlogd = (log(d_max(ic)) - log(d_min(ic)))/real(NSIZ(ic),kind=RP)

      do is0 = 1, NSIZ(ic)  !size bin
        d_lw(is0,ic) = exp(log(d_min(ic))+dlogd* real(is0-1,kind=RP)        )
        d_ct(is0,ic) = exp(log(d_min(ic))+dlogd*(real(is0  ,kind=RP)-0.5_RP))
        d_up(is0,ic) = exp(log(d_min(ic))+dlogd* real(is0  ,kind=RP)        )
      enddo !is (1:n_siz(ic))
      dk    = (k_max(ic) - k_min(ic))/real(n_kap(ic),kind=RP)
      do ik = 1, n_kap(ic)  !size bin
        k_lw(ik,ic) = k_min(ic) + dk  * real(ik-1,kind=RP)
        k_ct(ik,ic) = k_min(ic) + dk  *(real(ik  ,kind=RP)-0.5_RP)
        k_up(ik,ic) = k_min(ic) + dk  * real(ik  ,kind=RP)
      enddo !ik (1:n_kap(ic))

    enddo !ic (1:n_ctg)
!    ik  = 1       !only one kappa bin

    m0t = m0_init !total M0 [#/m3]
    dgt = dg_init ![m]
    sgt = sg_init ![-]

    if ( m0t <= cleannumber ) then
       m0t = cleannumber
       dgt = 0.1E-6_RP
       sgt = 1.3_RP
       if ( IO_L ) then
          write(IO_FID_LOG,*) '*** WARNING! Initial aerosol number is set as ', cleannumber, '[#/m3]'
       endif
    endif

    m2t = m0t*dgt**(2.d0) *dexp(2.0d0 *(dlog(real(sgt,kind=DP))**2.d0)) !total M2 [m2/m3]
    m3t = m0t*dgt**(3.d0) *dexp(4.5d0 *(dlog(real(sgt,kind=DP))**2.d0)) !total M3 [m3/m3]
    mst = m3t*pi6*conv_vl_ms                              !total Ms [ug/m3]

    do ic = 1, n_ctg
    !aerosol_procs initial condition
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, NSIZ(ic)
       if (dgt >= d_lw(is0,ic) .and. dgt < d_up(is0,ic)) then
         aerosol_procs(ia_m0,is0,ik,ic) = aerosol_procs(ia_m0,is0,ik,ic) + m0t          ![#/m3]
         aerosol_procs(ia_m2,is0,ik,ic) = aerosol_procs(ia_m2,is0,ik,ic) + m2t          ![m2/m3]
         aerosol_procs(ia_m3,is0,ik,ic) = aerosol_procs(ia_m3,is0,ik,ic) + m3t          ![m3/m3]
         aerosol_procs(ia_ms,is0,ik,ic) = aerosol_procs(ia_ms,is0,ik,ic) + mst*1.E-9_RP ![kg/m3]
       elseif (dgt < d_lw(1,ic)) then
         aerosol_procs(ia_m0,1 ,ik,ic) = aerosol_procs(ia_m0,1 ,ik,ic) + m0t          ![#/m3]
         aerosol_procs(ia_m2,1 ,ik,ic) = aerosol_procs(ia_m2,1 ,ik,ic) + m2t          ![m2/m3]
         aerosol_procs(ia_m3,1 ,ik,ic) = aerosol_procs(ia_m3,1 ,ik,ic) + m3t          ![m3/m3]
         aerosol_procs(ia_ms,1 ,ik,ic) = aerosol_procs(ia_ms,1 ,ik,ic) + mst*1.E-9_RP ![kg/m3]
       elseif (dgt >= d_up(NSIZ(ic),ic)) then
         aerosol_procs(ia_m0,NSIZ(ic),ik,ic) = aerosol_procs(ia_m0,NSIZ(ic),ik,ic) + m0t          ![#/m3]
         aerosol_procs(ia_m2,NSIZ(ic),ik,ic) = aerosol_procs(ia_m2,NSIZ(ic),ik,ic) + m2t          ![m2/m3]
         aerosol_procs(ia_m3,NSIZ(ic),ik,ic) = aerosol_procs(ia_m3,NSIZ(ic),ik,ic) + m3t          ![m3/m3]
         aerosol_procs(ia_ms,NSIZ(ic),ik,ic) = aerosol_procs(ia_ms,NSIZ(ic),ik,ic) + mst*1.E-9_RP ![kg/m3]
       endif
    enddo
    enddo
    enddo

    CCN(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          pott = RHOT(k,i,j) / DENS(k,i,j)
          qdry = 1.0_RP
          do iq = QQS, QQE
             qdry = qdry - QTRC(k,i,j,iq)
          enddo
          cva = qdry * CONST_CVdry
          do iq = QQS, QQE
             cva = cva + QTRC(k,i,j,iq) * AQ_CV(iq)
          enddo
          Rmoist = CONST_Rdry * qdry + CONST_Rvap * QTRC(k,i,j,I_QV)
          cpa = cva + Rmoist
          pres = CONST_PRE00 &
                * ( DENS(k,i,j)*Rmoist*pott/CONST_PRE00 )**( cpa/(cpa-Rmoist) )
          temp = pres / ( DENS(k,i,j) * Rmoist )
          !--- calculate super saturation of water
          call SATURATION_pres2qsat_liq( qsat_tmp,temp,pres )
          ssliq = QTRC(k,i,j,I_QV)/qsat_tmp - 1.0_RP
          call aerosol_activation_kajino13 &
                                 (c_kappa, ssliq, temp, ia_m0, ia_m2, ia_m3, &
                                  N_ATR,n_siz_max,n_kap_max,n_ctg,NSIZ,n_kap, &
                                  d_ct,aerosol_procs, aerosol_activ)

         do is0 = 1, NSIZ(ic_mix)
           CCN(k,i,j) = CCN(k,i,j) + aerosol_activ(ia_m0,is0,ik_out,ic_mix)
         enddo

       enddo

    enddo
    enddo

    conc_gas(:) = 0.0_RP
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
     it = 0
     do ic = 1, n_ctg       !category
     do ik = 1, NKAP(ic)   !kappa bin
     do is0 = 1, NSIZ(ic)   !size bin
     do ia0 = 1, N_ATR       !attributes
        it = it + 1
        QTRC(k,i,j,QAES-1+it) = aerosol_procs(ia0,is0,ik,ic)/DENS(k,i,j) !#,m2,m3,kg/m3 -> #,m2,m3,kg/kg
     enddo !ia0 (1:N_ATR )
     enddo !is (1:n_siz(ic)  )
     enddo !ik (1:n_kap(ic)  )
     enddo !ic (1:n_ctg      )
     QTRC(k,i,j,QAEE-GAS_CTG+1:QAEE-GAS_CTG+GAS_CTG) = conc_gas(1:GAS_CTG)/DENS(k,i,j)*1.E-9_RP !mixing ratio [kg/kg]
    enddo
    enddo
    enddo

    return

  end subroutine AEROSOL_setup
  !----------------------------------------------------------------------------------------
  ! subroutine 4. aerosol_activation
  !   Abdul-Razzak et al.,   JGR, 1998 [AR98]
  !   Abdul-Razzak and Ghan, JGR, 2000 [AR00]
  !----------------------------------------------------------------------------------------
  subroutine aerosol_activation_kajino13(c_kappa, super, temp_k, &
                                ia_m0, ia_m2, ia_m3, &
                                n_atr,n_siz_max,n_kap_max,n_ctg,n_siz,n_kap, &
                                d_ct, aerosol_procs, aerosol_activ)
  use scale_const, only: &
      CONST_Mvap, &            ! mean molecular weight for water vapor [g/mol]
      CONST_DWATR, &           ! water density   [kg/m3]
      CONST_R, &               ! universal gas constant             [J/mol/K]
      CONST_TEM00
  implicit none
    !i/o variables
      real(RP),intent(in) :: super      ! supersaturation                                 [-]
      real(RP),intent(in) :: c_kappa    ! hygroscopicity of condensable mass              [-]
      real(RP),intent(in) :: temp_k     ! temperature
      integer, intent(in) :: ia_m0, ia_m2, ia_m3
      integer, intent(in) :: n_atr
      integer, intent(in) :: n_siz_max
      integer, intent(in) :: n_kap_max
      integer, intent(in) :: n_ctg
      real(RP),intent(in) :: d_ct(n_siz_max,n_ctg)
      real(RP) :: aerosol_procs(n_atr,n_siz_max,n_kap_max,n_ctg)
      real(RP) :: aerosol_activ(n_atr,n_siz_max,n_kap_max,n_ctg)
      integer, intent(in) :: n_siz(n_ctg), n_kap(n_ctg)
    !local variables
      real(RP),parameter :: two3 = 2._RP/3._RP
      real(RP),parameter :: rt2  = sqrt(2._RP)
      real(RP),parameter :: twort2  = rt2       ! 2/sqrt(2) = sqrt(2)
      real(RP),parameter :: thrrt2  = 3._RP/rt2 ! 3/sqrt(2)
      real(RP) :: smax_inv                ! inverse
      real(RP) :: am,scrit_am,aa,tc,st,bb,ac
      real(RP) :: m0t,m2t,m3t,dgt,sgt,dm2
      real(RP) :: d_crit                  ! critical diameter
      real(RP) :: tmp1, tmp2, tmp3        !
      real(RP) :: ccn_frc,cca_frc,ccv_frc ! activated number,area,volume
      integer  :: is0, ik, ic

      aerosol_activ(:,:,:,:) = 0._RP

      if (super.le.0._RP) return

      smax_inv = 1._RP / super

    !--- surface tension of water
      tc = temp_k - CONST_TEM00
      if (tc >= 0._RP ) then
        st  = 75.94_RP-0.1365_RP*tc-0.3827e-3_RP*tc**2._RP !Gittens, JCIS, 69, Table 4.
      else !t[deg C].lt.0.
        st  = 75.93_RP                +0.115_RP*tc        &  !Eq.5-12, pp.130, PK97
            + 6.818e-2_RP*tc**2._RP+6.511e-3_RP*tc**3._RP &
            + 2.933e-4_RP*tc**4._RP+6.283e-6_RP*tc**5._RP &
            + 5.285e-8_RP*tc**6._RP
      endif
      st      = st * 1.E-3_RP                    ![J/m2]

    !-- Kelvin effect
    !          [J m-2]  [kg mol-1]       [m3 kg-1] [mol K J-1] [K-1]
      aa  = 2._RP * st * CONST_Mvap * 1.E-3_RP / (CONST_DWATR * CONST_R * temp_k ) ![m] Eq.5 in AR98

      do ic = 1, n_ctg
      do ik = 1, n_kap(ic)
      do is0 = 1, n_siz(ic)
        m0t = aerosol_procs(ia_m0,is0,ik,ic)
        m2t = aerosol_procs(ia_m2,is0,ik,ic)
        m3t = aerosol_procs(ia_m3,is0,ik,ic)
        call diag_ds(m0t,m2t,m3t,dgt,sgt,dm2)
        if (dgt <= 0._RP) dgt = d_ct(is0,ic)
        if (sgt <= 0._RP) sgt = 1.3_RP
        am  = dgt * 0.5_RP  !geometric dry mean radius [m]
        bb  = c_kappa
        if (bb > 0._RP .and. am > 0._RP ) then
          scrit_am = 2._RP/sqrt(bb)*(aa/(3._RP*am))**1.5_RP !AR00 Eq.9
        else
          scrit_am = 0._RP
        endif
        ac     = am * (scrit_am * smax_inv)**two3      !AR00 Eq.12
        d_crit = ac * 2._RP
        tmp1   = log(d_crit) - log(dgt)
        tmp2   = 1._RP/(rt2*log(sgt))
        tmp3   = log(sgt)
        ccn_frc= 0.5_RP*(1._RP-erf(tmp1*tmp2))
        cca_frc= 0.5_RP*(1._RP-erf(tmp1*tmp2-twort2*tmp3))
        ccv_frc= 0.5_RP*(1._RP-erf(tmp1*tmp2-thrrt2*tmp3))
        aerosol_activ(ia_m0,is0,ik,ic) = ccn_frc * aerosol_procs(ia_m0,is0,ik,ic)
        aerosol_activ(ia_m2,is0,ik,ic) = cca_frc * aerosol_procs(ia_m2,is0,ik,ic)
        aerosol_activ(ia_m3,is0,ik,ic) = ccv_frc * aerosol_procs(ia_m3,is0,ik,ic)
      enddo !is(1:n_siz(ic))
      enddo !ik(1:n_kap(ic))
      enddo !ic(1:n_ctg)

      return
  end subroutine aerosol_activation_kajino13
  !----------------------------------------------------------------------------------------
  subroutine diag_ds(m0,m2,m3,  & !i
                     dg,sg,dm2)   !o
    implicit none
    real(RP)            :: m0,m2,m3,dg,sg,m3_bar,m2_bar
    real(RP)            :: m2_new,m2_old,dm2
    real(RP), parameter :: sgmax=2.5_RP
    real(RP), parameter :: rk1=2._RP
    real(RP), parameter :: rk2=3._RP
    real(RP), parameter :: ratio  =rk1/rk2
    real(RP), parameter :: rk1_hat=1._RP/(ratio*(rk2-rk1))
    real(RP), parameter :: rk2_hat=ratio/(rk1-rk2)

    dm2=0._RP

    if (m0 <= 0._RP .or. m2 <= 0._RP .or. m3 <= 0._RP) then
      m0=0._RP
      m2=0._RP
      m3=0._RP
      dg=-1._RP
      sg=-1._RP
      return
    endif

    m2_old = m2
    m3_bar = m3/m0
    m2_bar = m2/m0
    dg     = m2_bar**rk1_hat*m3_bar**rk2_hat

    if (m2_bar/m3_bar**ratio < 1._RP) then !stdev > 1.
      sg     = exp(sqrt(2._RP/(rk1*(rk1-rk2))  &
             * log(m2_bar/m3_bar**ratio) ))
    endif

    if (sg > sgmax) then
      print *,'sg=',sg
      sg = sgmax
      call diag_d2(m0,m3,sg, & !i
                   m2,dg     ) !o
  !   print *,'warning! modified as sg exceeded sgmax (diag_ds)'
    endif

    if (m2_bar/m3_bar**ratio >= 1._RP) then !stdev < 1.
      sg = 1._RP
      call diag_d2(m0,m3,sg, & !i
                   m2,dg     ) !o
  !   print *,'warning! modified as sg < 1. (diag_ds)'
    endif

    m2_new = m2
    dm2    = m2_old - m2_new !m2_pres - m2_diag

    return
  end subroutine diag_ds
  !----------------------------------------------------------------------------------------
  ! subroutine 6. diag_d2
  !----------------------------------------------------------------------------------------
  subroutine diag_d2(m0,m3,sg, & !i
                       m2,dg     ) !o
    implicit none
    real(RP)            :: dg,sg,m0,m2,m3,aaa
    real(RP), parameter :: one3=1._RP/3._RP
    aaa = m0               * exp( 4.5_RP * (log(sg)**2._RP) )
    dg  =(m3/aaa)**one3
    m2  = m0 * dg ** 2._RP * exp( 2.0_RP * (log(sg)**2._RP) )

    return
  end subroutine diag_d2
  !--------------------------------------------------------------
  !> Setup aerosol condition for Spectral Bin Microphysics (SBM) model
  subroutine SBMAERO_setup
    use scale_const, only: &
       PI => CONST_PI
    use scale_tracer_suzuki10, only: &
       nccn, nbin
    implicit none

#ifndef DRY
    real(RP) :: xasta, xaend, dxaer
    real(RP), allocatable :: xabnd( : ), xactr( : )

    real(RP) :: F0_AERO      = 1.E+7_RP
    real(RP) :: R0_AERO      = 1.E-7_RP
    real(RP) :: R_MAX        = 1.E-06_RP
    real(RP) :: R_MIN        = 1.E-08_RP
    real(RP) :: A_ALPHA      = 3.0_RP
    real(RP) :: RHO_AERO     = 2.25E+03_RP

    NAMELIST / PARAM_SBMAERO / &
       F0_AERO,      &
       R0_AERO,      &
       R_MAX,        &
       R_MIN,        &
       A_ALPHA,      &
       RHO_AERO

    integer :: ierr
    integer :: iq, i, j, k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit aerobin] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SBMAERO,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found namelist. default value used'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist SBMAERO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_SBMAERO)

    if( nccn /= 0 ) then
      allocate( gan( nccn ) )
      allocate( xactr(nccn) )
      allocate( xabnd(nccn+1) )

      xasta = log( RHO_AERO*4.0_RP/3.0_RP*pi * ( R_MIN )**3 )
      xaend = log( RHO_AERO*4.0_RP/3.0_RP*pi * ( R_MAX )**3 )
      dxaer = ( xaend-xasta )/nccn
      do iq = 1, nccn+1
        xabnd( iq ) = xasta + dxaer*( iq-1 )
      enddo
      do iq = 1, nccn
        xactr( iq ) = ( xabnd( iq )+xabnd( iq+1 ) )*0.5_RP
      enddo
      do iq = 1, nccn
        gan( iq ) = faero( F0_AERO,R0_AERO,xactr( iq ), A_ALPHA, RHO_AERO )*exp( xactr(iq) )
      enddo
    endif

    !--- Hydrometeor is zero at initial time for Bin method
    do iq = 2,  QQA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
        QTRC(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo

    !-- Aerosol distribution
    if( nccn /= 0 ) then
     do iq = QQA+1, QA
     do j  = JS, JE
     do i  = IS, IE
     do k  = KS, KE
       QTRC(k,i,j,iq) = gan(iq-QQA) !/ DENS(k,i,j)
     enddo
     enddo
     enddo
     enddo

     deallocate( xactr )
     deallocate( xabnd )
    endif

#endif
    return
  end subroutine SBMAERO_setup

  !-----------------------------------------------------------------------------
  function faero( f0,r0,x,alpha,rhoa )
    use scale_const, only: &
       pi => CONST_PI
    implicit none

    real(RP), intent(in) ::  x, f0, r0, alpha, rhoa
    real(RP) :: faero
    real(RP) :: rad
    !---------------------------------------------------------------------------

    rad = ( exp(x) * 3.0_RP / 4.0_RP / pi / rhoa )**(1.0_RP/3.0_RP)

    faero = f0 * (rad/r0)**(-alpha)

    return
  end function faero

  !-----------------------------------------------------------------------------
  !> flux setup
  subroutine flux_setup
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain    => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow    => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_up   => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFLX_LW_dn   => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFLX_SW_up   => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFLX_SW_dn   => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOAFLX_LW_up => ATMOS_PHY_RD_TOAFLX_LW_up, &
       TOAFLX_LW_dn => ATMOS_PHY_RD_TOAFLX_LW_dn, &
       TOAFLX_SW_up => ATMOS_PHY_RD_TOAFLX_SW_up, &
       TOAFLX_SW_dn => ATMOS_PHY_RD_TOAFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_downall
    implicit none
    ! Flux from Atmosphere
    real(RP) :: FLX_rain      = 0.0_RP ! surface rain flux                         [kg/m2/s]
    real(RP) :: FLX_snow      = 0.0_RP ! surface snow flux                         [kg/m2/s]
    real(RP) :: FLX_LW_dn     = 0.0_RP ! surface downwad long-wave  radiation flux [J/m2/s]
    real(RP) :: FLX_SW_dn     = 0.0_RP ! surface downwad short-wave radiation flux [J/m2/s]

    NAMELIST / PARAM_MKINIT_FLUX / &
       FLX_rain,      &
       FLX_snow,      &
       FLX_LW_dn,     &
       FLX_SW_dn

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_FLUX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_FLUX. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_FLUX)

    do j = JS, JE
    do i = IS, IE
       SFLX_rain   (i,j) = FLX_rain
       SFLX_snow   (i,j) = FLX_snow

       SFLX_LW_up  (i,j) = 0.0_RP
       SFLX_LW_dn  (i,j) = FLX_LW_dn
       SFLX_SW_up  (i,j) = 0.0_RP
       SFLX_SW_dn  (i,j) = FLX_SW_dn

       TOAFLX_LW_up(i,j) = 0.0_RP
       TOAFLX_LW_dn(i,j) = 0.0_RP
       TOAFLX_SW_up(i,j) = 0.0_RP
       TOAFLX_SW_dn(i,j) = 0.0_RP

       SFLX_rad_dn (i,j,1,1) = FLX_SW_dn
       SFLX_rad_dn (i,j,1,2) = 0.0_RP
       SFLX_rad_dn (i,j,2,1) = FLX_LW_dn
       SFLX_rad_dn (i,j,2,2) = 0.0_RP
    enddo
    enddo

    return
  end subroutine flux_setup

  !-----------------------------------------------------------------------------
  !> Land setup
  subroutine land_setup
    use mod_land_vars, only: &
       LAND_TEMP,       &
       LAND_WATER,      &
       LAND_SFC_TEMP,   &
       LAND_SFC_albedo
    implicit none
    ! Land state
    real(RP) :: LND_TEMP                ! soil temperature           [K]
    real(RP) :: LND_WATER     = 0.15_RP ! soil moisture              [m3/m3]
    real(RP) :: SFC_TEMP                ! land skin temperature      [K]
    real(RP) :: SFC_albedo_LW = 0.01_RP ! land surface albedo for LW [0-1]
    real(RP) :: SFC_albedo_SW = 0.20_RP ! land surface albedo for SW [0-1]

    integer :: i, j
    integer :: ierr

    NAMELIST /PARAM_MKINIT_LAND/ &
       LND_TEMP,      &
       LND_WATER,     &
       SFC_TEMP,      &
       SFC_albedo_LW, &
       SFC_albedo_SW

    LND_TEMP = THETAstd
    SFC_TEMP = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_LAND. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_LAND)

    do j = JS, JE
    do i = IS, IE
       LAND_TEMP      (:,i,j)    = LND_TEMP
       LAND_WATER     (:,i,j)    = LND_WATER
       LAND_SFC_TEMP  (i,j)      = SFC_TEMP
       LAND_SFC_albedo(i,j,I_LW) = SFC_albedo_LW
       LAND_SFC_albedo(i,j,I_SW) = SFC_albedo_SW
    end do
    end do

    return
  end subroutine land_setup

  !-----------------------------------------------------------------------------
  !> Ocean setup
  subroutine ocean_setup
    use mod_ocean_vars, only: &
       OCEAN_TEMP,       &
       OCEAN_SFC_TEMP,   &
       OCEAN_SFC_albedo, &
       OCEAN_SFC_Z0M,    &
       OCEAN_SFC_Z0H,    &
       OCEAN_SFC_Z0E
    implicit none
    ! Ocean state
    real(RP) :: OCN_TEMP                  ! ocean temperature           [K]
    real(RP) :: SFC_TEMP                  ! ocean skin temperature      [K]
    real(RP) :: SFC_albedo_LW = 0.04_RP   ! ocean surface albedo for LW [0-1]
    real(RP) :: SFC_albedo_SW = 0.05_RP   ! ocean surface albedo for SW [0-1]
    real(RP) :: SFC_Z0M       = 1.0e-4_RP ! ocean surface roughness length (momentum) [m]
    real(RP) :: SFC_Z0H       = 1.0e-4_RP ! ocean surface roughness length (heat) [m]
    real(RP) :: SFC_Z0E       = 1.0e-4_RP ! ocean surface roughness length (vapor) [m]

    integer :: i, j
    integer :: ierr

    NAMELIST /PARAM_MKINIT_OCEAN/ &
       OCN_TEMP,      &
       SFC_TEMP,      &
       SFC_albedo_LW, &
       SFC_albedo_SW, &
       SFC_Z0M,       &
       SFC_Z0H,       &
       SFC_Z0E

    OCN_TEMP = THETAstd
    SFC_TEMP = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_OCEAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_OCEAN)


    do j = JS, JE
    do i = IS, IE
       OCEAN_TEMP      (i,j)      = OCN_TEMP
       OCEAN_SFC_TEMP  (i,j)      = SFC_TEMP
       OCEAN_SFC_albedo(i,j,I_LW) = SFC_albedo_LW
       OCEAN_SFC_albedo(i,j,I_SW) = SFC_albedo_SW
       OCEAN_SFC_Z0M   (i,j)      = SFC_Z0M
       OCEAN_SFC_Z0H   (i,j)      = SFC_Z0H
       OCEAN_SFC_Z0E   (i,j)      = SFC_Z0E
    end do
    end do

    return
  end subroutine ocean_setup

  !-----------------------------------------------------------------------------
  !> Urban setup
  subroutine urban_setup
    use mod_urban_vars, only: &
       URBAN_TR,         &
       URBAN_TB,         &
       URBAN_TG,         &
       URBAN_TC,         &
       URBAN_QC,         &
       URBAN_UC,         &
       URBAN_TRL,        &
       URBAN_TBL,        &
       URBAN_TGL,        &
       URBAN_RAINR,      &
       URBAN_RAINB,      &
       URBAN_RAING,      &
       URBAN_ROFF,       &
       URBAN_SFC_TEMP,   &
       URBAN_SFC_albedo
    implicit none
    ! urban state
    real(RP) :: URB_ROOF_TEMP          ! Surface temperature of roof [K]
    real(RP) :: URB_BLDG_TEMP          ! Surface temperature of building [K
    real(RP) :: URB_GRND_TEMP          ! Surface temperature of ground [K]
    real(RP) :: URB_CNPY_TEMP          ! Diagnostic canopy air temperature
    real(RP) :: URB_CNPY_HMDT = 0.0_RP ! Diagnostic canopy humidity [-]
    real(RP) :: URB_CNPY_WIND = 0.0_RP ! Diagnostic canopy wind [m/s]
    real(RP) :: URB_ROOF_LAYER_TEMP    ! temperature in layer of roof [K]
    real(RP) :: URB_BLDG_LAYER_TEMP    ! temperature in layer of building [
    real(RP) :: URB_GRND_LAYER_TEMP    ! temperature in layer of ground [K]
    real(RP) :: URB_ROOF_RAIN = 0.0_RP ! temperature in layer of roof [K]
    real(RP) :: URB_BLDG_RAIN = 0.0_RP ! temperature in layer of building [
    real(RP) :: URB_GRND_RAIN = 0.0_RP ! temperature in layer of ground [K]
    real(RP) :: URB_RUNOFF    = 0.0_RP ! temperature in layer of ground [K]
    real(RP) :: URB_SFC_TEMP           ! Grid average of surface temperature [K]
    real(RP) :: URB_ALB_LW    = 0.0_RP ! Grid average of surface albedo for LW [0-1]
    real(RP) :: URB_ALB_SW    = 0.0_RP ! Grid average of surface albedo for SW [0-1]

    integer :: i, j
    integer :: ierr

    NAMELIST /PARAM_MKINIT_URBAN/ &
       URB_ROOF_TEMP,       &
       URB_BLDG_TEMP,       &
       URB_GRND_TEMP,       &
       URB_CNPY_TEMP,       &
       URB_CNPY_HMDT,       &
       URB_CNPY_WIND,       &
       URB_ROOF_LAYER_TEMP, &
       URB_BLDG_LAYER_TEMP, &
       URB_GRND_LAYER_TEMP, &
       URB_ROOF_RAIN,       &
       URB_BLDG_RAIN,       &
       URB_GRND_RAIN,       &
       URB_RUNOFF,          &
       URB_SFC_TEMP,        &
       URB_ALB_LW,          &
       URB_ALB_SW

    URB_ROOF_TEMP       = THETAstd
    URB_BLDG_TEMP       = THETAstd
    URB_GRND_TEMP       = THETAstd
    URB_CNPY_TEMP       = THETAstd
    URB_ROOF_LAYER_TEMP = THETAstd
    URB_BLDG_LAYER_TEMP = THETAstd
    URB_GRND_LAYER_TEMP = THETAstd

    URB_SFC_TEMP        = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_URBAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_URBAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_URBAN)


    do j = JS, JE
    do i = IS, IE
       URBAN_TR   (i,j)   = URB_ROOF_TEMP
       URBAN_TB   (i,j)   = URB_BLDG_TEMP
       URBAN_TG   (i,j)   = URB_GRND_TEMP
       URBAN_TC   (i,j)   = URB_CNPY_TEMP
       URBAN_QC   (i,j)   = URB_CNPY_HMDT
       URBAN_UC   (i,j)   = URB_CNPY_WIND
       URBAN_TRL  (:,i,j) = URB_ROOF_LAYER_TEMP
       URBAN_TBL  (:,i,j) = URB_BLDG_LAYER_TEMP
       URBAN_TGL  (:,i,j) = URB_GRND_LAYER_TEMP
       URBAN_RAINR(i,j)   = URB_ROOF_RAIN
       URBAN_RAINB(i,j)   = URB_BLDG_RAIN
       URBAN_RAING(i,j)   = URB_GRND_RAIN
       URBAN_ROFF (i,j)   = URB_RUNOFF
       URBAN_SFC_TEMP  (i,j)      = URB_SFC_TEMP
       URBAN_SFC_albedo(i,j,I_LW) = URB_ALB_LW
       URBAN_SFC_albedo(i,j,I_SW) = URB_ALB_SW
    end do
    end do

    return
  end subroutine urban_setup

  !-----------------------------------------------------------------------------
  !> Read sounding data from file
  subroutine read_sounding( &
       DENS, VELX, VELY, POTT, QV )
    implicit none
    real(RP), intent(out) :: DENS(KA)
    real(RP), intent(out) :: VELX(KA)
    real(RP), intent(out) :: VELY(KA)
    real(RP), intent(out) :: POTT(KA)
    real(RP), intent(out) :: QV  (KA)

    real(RP) :: TEMP(KA)
    real(RP) :: PRES(KA)
    real(RP) :: QC  (KA)

    character(len=H_LONG) :: ENV_IN_SOUNDING_file = ''

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(RP) :: SFC_THETA            ! surface potential temperature [K]
    real(RP) :: SFC_PRES             ! surface pressure [hPa]
    real(RP) :: SFC_QV               ! surface watervapor [g/kg]

    real(RP) :: EXP_z   (EXP_klim+1) ! height      [m]
    real(RP) :: EXP_pott(EXP_klim+1) ! potential temperature [K]
    real(RP) :: EXP_qv  (EXP_klim+1) ! water vapor [g/kg]
    real(RP) :: EXP_u   (EXP_klim+1) ! velocity u  [m/s]
    real(RP) :: EXP_v   (EXP_klim+1) ! velocity v  [m/s]

    real(RP) :: fact1, fact2
    integer :: k, kref
    integer :: fid
    integer :: ierr

    NAMELIST /PARAM_MKINIT_SOUNDING/ &
         ENV_IN_SOUNDING_file

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SOUNDING,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SOUNDING. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_SOUNDING)

    !--- prepare sounding profile
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding file:', trim(ENV_IN_SOUNDING_file)
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

       if ( ierr /= 0 ) then
          if( IO_L ) write(*,*) 'xxx Input file not found!'
       endif

       !--- read sounding file till end
       read(fid,*) SFC_PRES, SFC_THETA, SFC_QV

       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pressure [hPa]',     SFC_PRES
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pot. temp  [K]',     SFC_THETA
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface water vapor [g/kg]', SFC_QV

       do k = 2, EXP_klim
          read(fid,*,iostat=ierr) EXP_z(k), EXP_pott(k), EXP_qv(k), EXP_u(k), EXP_v(k)
          if ( ierr /= 0 ) exit
       enddo

       EXP_kmax = k - 1
    close(fid)

    ! Boundary
    EXP_z   (1)          = 0.0_RP
    EXP_pott(1)          = SFC_THETA
    EXP_qv  (1)          = SFC_QV
    EXP_u   (1)          = EXP_u   (2)
    EXP_v   (1)          = EXP_v   (2)
    EXP_z   (EXP_kmax+1) = 100.E3_RP
    EXP_pott(EXP_kmax+1) = EXP_pott(EXP_kmax)
    EXP_qv  (EXP_kmax+1) = EXP_qv  (EXP_kmax)
    EXP_u   (EXP_kmax+1) = EXP_u   (EXP_kmax)
    EXP_v   (EXP_kmax+1) = EXP_v   (EXP_kmax)

    do k = 1, EXP_kmax+1
       EXP_qv(k) = EXP_qv(k) * 1.E-3_RP ! [g/kg]->[kg/kg]
    enddo

    ! calc in dry condition
    pres_sfc = SFC_PRES * 1.E2_RP ! [hPa]->[Pa]
    pott_sfc = SFC_THETA
    qv_sfc   = SFC_QV * 1.E-3_RP ! [g/kg]->[kg/kg]
    qc_sfc   = 0.0_RP

    !--- linear interpolate to model grid
    do k = KS, KE
       qc(k) = 0.0_RP

       do kref = 2, EXP_kmax+1
          if (       GRID_CZ(k) >  EXP_z(kref-1) &
               .AND. GRID_CZ(k) <= EXP_z(kref  ) ) then

             fact1 = ( EXP_z(kref) - GRID_CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( GRID_CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

             pott(k) = EXP_pott(kref-1) * fact1 &
                     + EXP_pott(kref  ) * fact2
             qv  (k) = EXP_qv  (kref-1) * fact1 &
                     + EXP_qv  (kref  ) * fact2
             velx(k) = EXP_u   (kref-1) * fact1 &
                     + EXP_u   (kref  ) * fact2
             vely(k) = EXP_v   (kref-1) * fact1 &
                     + EXP_v   (kref  ) * fact2
          endif
       enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS(:),         & ! [OUT]
                               temp(:),         & ! [OUT]
                               pres(:),         & ! [OUT]
                               pott(:),         & ! [IN]
                               qv  (:),         & ! [IN]
                               qc  (:),         & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1) ) ! [IN]

    return
  end subroutine read_sounding

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform + random disturbance )
  subroutine MKINIT_planestate
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =   0.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_THETA               ! potential temperature of environment [K]
    real(RP) :: ENV_TLAPS    =   0.0_RP ! Lapse rate of THETA [K/m]
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =   0.0_RP ! relative humidity of environment [%]
    ! Disturbance
    real(RP) :: RANDOM_THETA =   0.0_RP ! amplitude of random disturbance theta
    real(RP) :: RANDOM_U     =   0.0_RP ! amplitude of random disturbance u
    real(RP) :: RANDOM_V     =   0.0_RP ! amplitude of random disturbance v
    real(RP) :: RANDOM_RH    =   0.0_RP ! amplitude of random disturbance RH

    NAMELIST / PARAM_MKINIT_PLANESTATE / &
       SFC_THETA,    &
       SFC_PRES,     &
       SFC_RH,       &
       ENV_THETA,    &
       ENV_TLAPS,    &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       RANDOM_THETA, &
       RANDOM_U,     &
       RANDOM_V,     &
       RANDOM_RH

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit Horiz_UNIFORM] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_PLANESTATE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_PLANESTATE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_PLANESTATE)

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE
       pott_sfc(1,i,j) = SFC_THETA
       pres_sfc(1,i,j) = SFC_PRES
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    if ( ENV_THETA < 0.0_RP ) then ! use isa profile

       call PROFILE_isa( KA, KS, KE,      & ! [IN]
                         IA, IS, IE,      & ! [IN]
                         JA, JS, JE,      & ! [IN]
                         pott_sfc(1,:,:), & ! [IN]
                         pres_sfc(1,:,:), & ! [IN]
                         REAL_CZ (:,:,:), & ! [IN]
                         pott    (:,:,:)  ) ! [OUT]

    else

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * REAL_CZ(k,i,j)
       enddo
       enddo
       enddo

    endif

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc(1,:,:), temp_sfc(1,:,:), pres_sfc(1,:,:) )
    call SATURATION_pres2qsat_all( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(1,i,j)

       do k = KS, KE
          qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat(k,i,j)
       enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       pott_sfc(1,i,j) = pott_sfc(1,i,j) + rndm(KS-1,i,j) * RANDOM_THETA

       do k = KS, KE
          pott(k,i,j) = pott(k,i,j) + rndm(k,i,j) * RANDOM_THETA
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_U ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_V ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
    enddo
    enddo
    enddo

    call flux_setup

    call ocean_setup

    return
  end subroutine MKINIT_planestate

  !-----------------------------------------------------------------------------
  !> Make initial state for tracer bubble experiment
  subroutine MKINIT_tracerbubble
    implicit none

#ifndef DRY
    ! Surface state
    real(RP)               :: SFC_THETA            ! surface potential temperature [K]
    real(RP)               :: SFC_PRES             ! surface pressure [Pa]
    ! Environment state
    real(RP)               :: ENV_THETA            ! potential temperature of environment [K]
    real(RP)               :: ENV_U     =   0.0_RP ! velocity u of environment [m/s]
    real(RP)               :: ENV_V     =   0.0_RP ! velocity v of environment [m/s]
    ! Bubble
    character(len=H_SHORT) :: SHAPE_NC  = 'BUBBLE' ! BUBBLE or RECT
    real(RP)               :: BBL_NC    =   1.0_RP ! extremum of NC in bubble [kg/kg]

    NAMELIST / PARAM_MKINIT_TRACERBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       ENV_U,     &
       ENV_V,     &
       SHAPE_NC,  &
       BBL_NC

    real(RP), pointer :: shapeFac(:,:,:) => null()

    integer :: k, i, j, iq
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit TRACERBUBBLE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TRACERBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_TRACERBUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_TRACERBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       pott(k,1,1) = ENV_THETA
       qv  (k,1,1) = 0.0_RP
       qc  (k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U       * DENS(k,1,1)
       MOMY(k,i,j) = ENV_V       * DENS(k,1,1)
       RHOT(k,i,j) = pott(k,1,1) * DENS(k,1,1)

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
    enddo
    enddo
    enddo

    ! make tracer bubble
    if ( I_NC > 0 ) then

       select case(SHAPE_NC)
       case('BUBBLE')
          call BUBBLE_setup
          shapeFac => bubble
       case('RECT')
          call RECT_setup
          shapeFac => rect
       case default
          write(*,*) 'xxx SHAPE_NC=', trim(SHAPE_NC), ' cannot be used on advect. Check!'
          call PRC_MPIstop
       end select

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_NC) = BBL_NC * shapeFac(k,i,j)
       enddo
       enddo
       enddo
    else
       write(*,*) 'xxx tracer I_NC is not defined. Check!'
       call PRC_MPIstop
    endif

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       write(*,*) 'xxx SBM cannot be used on tracerbubble. Check!'
       call PRC_MPIstop
    endif
#endif

    return
  end subroutine MKINIT_tracerbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !! Default values are following by Straka et al. (1993)
  !!  BBL_TEMP = -15.0_RP   ! in temperature  [K]
  !!  BBL_CZ   =   3.0E3_RP ! center location [m]: z
  !!  BBL_CX   =  19.2E3_RP ! center location [m]: x
  !!  BBL_CY   =   1.0E2_RP ! center location [m]: y
  !!  BBL_RZ   =   2.0E3_RP ! bubble radius   [m]: z
  !!  BBL_RX   =   4.0E3_RP ! bubble radius   [m]: x
  !!  BBL_RY   =   1.0E3_RP ! bubble radius   [m]: y
  subroutine MKINIT_coldbubble
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_THETA               ! potential temperature of environment [K]
    ! Bubble
    real(RP) :: BBL_TEMP     = -15.0_RP ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_COLDBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       BBL_TEMP

    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit COLDBUBBLE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_COLDBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_COLDBUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_COLDBUBBLE)

    RovCP = Rdry / CPdry

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       pott(k,1,1) = ENV_THETA
       qv  (k,1,1) = 0.0_RP
       qc  (k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP

       ! make cold bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1)                                           &
                                   + BBL_TEMP * ( P00/pres(k,1,1) )**RovCP * bubble(k,i,j) )

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
    enddo
    enddo
    enddo

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       write(*,*) 'xxx SBM cannot be used on coldbubble. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine MKINIT_coldbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  subroutine MKINIT_lambwave
    implicit none

    ! Surface state
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_TEMP     = 300.0_RP ! temperature of environment [K]
    ! Bubble
    real(RP) :: BBL_PRES     =  100._RP ! extremum of pressure in bubble [Pa]

    NAMELIST / PARAM_MKINIT_LAMBWAVE / &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       ENV_TEMP,  &
       BBL_PRES

    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit LAMBWAVE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_LAMBWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_LAMBWAVE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_LAMBWAVE)

    RovCP = Rdry / CPdry

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = SFC_PRES/(Rdry*ENV_TEMP) * exp( - GRAV/(Rdry*ENV_TEMP) * GRID_CZ(k) )
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make pressure bubble
       pres(k,i,j) = DENS(k,i,j) * ENV_TEMP * Rdry + BBL_PRES * bubble(k,i,j)

       RHOT(k,i,j) = DENS(k,i,j) * ENV_TEMP * ( P00/pres(k,i,j) )**RovCP

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
    enddo
    enddo
    enddo

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       write(*,*) 'xxx SBM cannot be used on lambwave. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine MKINIT_lambwave

  !-----------------------------------------------------------------------------
  !> Make initial state for gravity wave experiment
  !! Default values are following by Skamarock and Klemp (1994)
  subroutine MKINIT_gravitywave
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_U        =  20.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_BVF      =  0.01_RP ! Brunt Vaisala frequencies of environment [1/s]
    ! Bubble
    real(RP) :: BBL_THETA    =  0.01_RP ! extremum of potential temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_GRAVITYWAVE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       ENV_BVF,   &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit GRAVITYWAVE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_GRAVITYWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_GRAVITYWAVE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_GRAVITYWAVE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       pott(k,1,1) = SFC_THETA * exp( ENV_BVF*ENV_BVF / GRAV * GRID_CZ(k) )
       qv  (k,1,1) = 0.0_RP
       qc  (k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,1,1)
       MOMY(k,i,j) = ENV_V * DENS(k,1,1)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
    enddo
    enddo
    enddo

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       write(*,*) 'xxx SBM cannot be used on gravitywave. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine MKINIT_gravitywave

  !-----------------------------------------------------------------------------
  !> Make initial state for Kelvin-Helmholtz wave experiment
  subroutine MKINIT_khwave
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA                  ! surface potential temperature [K]
    real(RP) :: SFC_PRES                   ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_L1_ZTOP    = 1900.0_RP ! top    height of the layer1 (low  THETA) [m]
    real(RP) :: ENV_L3_ZBOTTOM = 2100.0_RP ! bottom height of the layer3 (high THETA) [m]
    real(RP) :: ENV_L1_THETA   =  300.0_RP ! THETA in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_THETA   =  301.0_RP ! THETA in the layer3 (high THETA) [K]
    real(RP) :: ENV_L1_U       =    0.0_RP ! velocity u in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_U       =   20.0_RP ! velocity u in the layer3 (high THETA) [K]
    ! Disturbance
    real(RP) :: RANDOM_U       =    0.0_RP ! amplitude of random disturbance u

    NAMELIST / PARAM_MKINIT_KHWAVE / &
       SFC_THETA,      &
       SFC_PRES,       &
       ENV_L1_ZTOP,    &
       ENV_L3_ZBOTTOM, &
       ENV_L1_THETA,   &
       ENV_L3_THETA,   &
       ENV_L1_U,       &
       ENV_L3_U,       &
       RANDOM_U

    real(RP) :: fact

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit KHWAVE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_KHWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_KHWAVE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_KHWAVE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       fact = ( GRID_CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )
       fact = max( min( fact, 1.0_RP ), 0.0_RP )

       pott(k,1,1) = ENV_L1_THETA * ( 1.0_RP - fact ) &
                   + ENV_L3_THETA * (          fact )

       qv(k,1,1) = 0.0_RP
       qc(k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc(1,1,1), temp_sfc(1,1,1), pres_sfc(1,1,1) )
    call SATURATION_pres2qsat_all( qsat    (:,1,1), temp    (:,1,1), pres    (:,1,1) )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       RHOT(k,i,j) = DENS(k,1,1) * pott(k,1,1)

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       fact = ( GRID_CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )
       fact = max( min( fact, 1.0_RP ), 0.0_RP )

       MOMX(k,i,j) = ( ENV_L1_U * ( 1.0_RP - fact )                 &
                     + ENV_L3_U * (          fact )                 &
                     + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_U &
                     ) * DENS(k,i,j)
    enddo
    enddo
    enddo

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       write(*,*) 'xxx SBM cannot be used on khwave. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine MKINIT_khwave

  !-----------------------------------------------------------------------------
  !> Make initial state for turbulence experiment
  subroutine MKINIT_turbulence
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =   0.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_THETA               ! potential temperature of environment [K]
    real(RP) :: ENV_TLAPS    = 4.E-3_RP ! Lapse rate of THETA [K/m]
    real(RP) :: ENV_U        =   5.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =   0.0_RP ! relative humidity of environment [%]
    ! Disturbance
    real(RP) :: RANDOM_THETA =   1.0_RP ! amplitude of random disturbance theta
    real(RP) :: RANDOM_U     =   0.0_RP ! amplitude of random disturbance u
    real(RP) :: RANDOM_V     =   0.0_RP ! amplitude of random disturbance v
    real(RP) :: RANDOM_RH    =   0.0_RP ! amplitude of random disturbance RH

    NAMELIST / PARAM_MKINIT_TURBULENCE / &
       SFC_THETA,    &
       SFC_PRES,     &
       SFC_RH,       &
       ENV_THETA,    &
       ENV_TLAPS,    &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       RANDOM_THETA, &
       RANDOM_U,     &
       RANDOM_V,     &
       RANDOM_RH

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit TURBULENCE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TURBULENCE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_TURBULENCE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_TURBULENCE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       pott(k,1,1) = ENV_THETA + ENV_TLAPS * GRID_CZ(k)
       qv  (k,1,1) = 0.0_RP
       qc  (k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc(1,1,1), temp_sfc(1,1,1), pres_sfc(1,1,1) )
    call SATURATION_pres2qsat_all( qsat    (:,1,1), temp    (:,1,1), pres    (:,1,1) )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(1,1,1)
       qc_sfc(1,i,j) = 0.0_RP

       do k = KS, KE
          qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat(k,1,1)
          qc(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA + rndm(KS-1,i,j) * RANDOM_THETA

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * GRID_CZ(k) + rndm(k,i,j) * RANDOM_THETA
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_U ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_V ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
    enddo
    enddo
    enddo

#ifndef DRY
    if ( QA >= 2 ) then
       do iq = 2, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
       enddo
       enddo
       enddo
    endif
#endif

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       write(*,*) 'xxx SBM cannot be used on turbulence. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine MKINIT_turbulence

  !-----------------------------------------------------------------------------
  !> Make initial state for cavity flow
  subroutine MKINIT_cavityflow
    implicit none

    ! Nondimenstional numbers for a cavity flow problem
    real(RP) :: REYNOLDS_NUM = 4.D02
    real(RP) :: MACH_NUM     = 3.D-2
    real(RP) :: TEMP0        = 300.D0

    NAMELIST / PARAM_MKINIT_CAVITYFLOW / &
         REYNOLDS_NUM, &
         MACH_NUM,     &
         TEMP0

    real(RP) :: DENS0
    real(RP) :: TEMP
    real(RP) :: Gam, Cs2

    real(RP), parameter :: PRES0 = 1000.E+2_RP
    real(RP), parameter :: Ulid  = 10.0_RP

    integer :: k, i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit CAVITYFLOW] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_CAVITYFLOW,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_CAVITYFLOW. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_CAVITYFLOW)

    Gam   = CPdry / ( CPdry - Rdry )
    Cs2   = ( Ulid / MACH_NUM )**2
    TEMP  = Cs2   / ( Gam * Rdry )
    DENS0 = PRES0 / ( Rdry * TEMP )

    if( IO_L ) write(IO_FID_LOG,*) "DENS = ", DENS0
    if( IO_L ) write(IO_FID_LOG,*) "PRES = ", PRES0
    if( IO_L ) write(IO_FID_LOG,*) "TEMP = ", RHOT(10,10,4)/DENS0, TEMP
    if( IO_L ) write(IO_FID_LOG,*) "Ulid = ", Ulid
    if( IO_L ) write(IO_FID_LOG,*) "Cs   = ", sqrt(Cs2)

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j)   = DENS0
       MOMZ(k,i,j)   = 0.0_RP
       MOMX(k,i,j)   = 0.0_RP
       MOMY(k,i,j)   = 0.0_RP
       PRES(k,i,j)   = PRES0
       RHOT(k,i,j)   = P00/Rdry * (P00/PRES0)**((Rdry - CPdry)/CPdry)
       QTRC(k,i,j,:) = 0.0_RP
    enddo
    enddo
    enddo

    MOMX(KE+1:KA,:,:) = DENS0 * Ulid

    return
  end subroutine MKINIT_cavityflow

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform )
  subroutine MKINIT_mountainwave
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA      ! surface potential temperature [K]
    real(RP) :: SFC_PRES       ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_U = 0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V = 0.0_RP ! velocity v of environment [m/s]

    real(RP) :: SCORER = 2.E-3_RP ! Scorer parameter (~=N/U) [1/m]
    real(RP) :: BBL_NC =   0.0_RP ! extremum of NC in bubble [kg/kg]

    NAMELIST / PARAM_MKINIT_MOUNTAINWAVE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       SCORER,    &
       BBL_NC

    real(RP) :: Ustar2, N2

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit MOUNTAINWAVE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_MOUNTAINWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_MOUNTAINWAVE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_MOUNTAINWAVE)

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Ustar2 = ENV_U * ENV_U + ENV_V * ENV_V
       N2     = Ustar2 * (SCORER*SCORER)

       pott(k,i,j) = SFC_THETA * exp( N2 / GRAV * REAL_CZ(k,i,j) )
       qv  (k,i,j) = 0.0_RP
       qc  (k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = DENS(k,i,j)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U       * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V       * DENS(k,i,j)
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       QTRC(k,i,j,:) = 0.0_RP
    enddo
    enddo
    enddo

    ! optional : add tracer bubble
    if (  BBL_NC > 0.0_RP ) then
       if (  I_NC > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,I_NC) = BBL_NC * bubble(k,i,j)
          enddo
          enddo
          enddo
       else
          write(*,*) 'xxx tracer I_NC is not defined. Check!'
          call PRC_MPIstop
       endif
    endif

    return
  end subroutine MKINIT_mountainwave

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  subroutine MKINIT_warmbubble
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  80.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =  80.0_RP ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =  1.E3_RP ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  = 14.E3_RP ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS = 4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS = 3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(RP) :: BBL_THETA    =   1.0_RP ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit WARMBUBBLE] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_WARMBUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       if    ( GRID_CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( GRID_CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( GRID_CZ(k)-GRID_CZ(k-1) )
       else                                ! Layer 3
          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( GRID_CZ(k)-GRID_CZ(k-1) )
       endif
       qv(k,1,1) = 0.0_RP
       qc(k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc(1,1,1), temp_sfc(1,1,1), pres_sfc(1,1,1) )
    call SATURATION_pres2qsat_all( qsat    (:,1,1), temp    (:,1,1), pres    (:,1,1) )

    qv_sfc(1,1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1,1)
    do k = KS, KE
       if    ( GRID_CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       elseif( GRID_CZ(k) <= ENV_L2_ZTOP ) then ! Layer 2
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       else                                ! Layer 3
          qv(k,1,1) = 0.0_RP
       endif
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       QTRC(k,i,j,I_QV) = qv(k,1,1)
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_warmbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for supercell experiment
  subroutine MKINIT_supercell
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV(KA)

    ! Bubble
    real(RP) :: BBL_THETA = 3.D0 ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_SUPERCELL / &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit SUPERCELL] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SUPERCELL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SUPERCELL. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_SUPERCELL)

    call read_sounding( RHO, VELX, VELY, POTT, QV ) ! (out)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = RHO(k) * VELX(k)
       MOMY(k,i,j) = RHO(k) * VELY(k)

       ! make warm bubble
       RHOT(k,i,j) = RHO(k) * ( POTT(k) + BBL_THETA * bubble(k,i,j) )

       QTRC(k,i,j,I_QV) = QV(k)
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_supercell

  !-----------------------------------------------------------------------------
  !> Make initial state for squallline experiment
  subroutine MKINIT_squallline
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV(KA)

    real(RP) :: RANDOM_THETA =  0.01_RP
    real(RP) :: OFFSET_velx  = 12.0_RP
    real(RP) :: OFFSET_vely  = -2.0_RP

    NAMELIST / PARAM_MKINIT_SQUALLLINE / &
       RANDOM_THETA,         &
       OFFSET_velx,          &
       OFFSET_vely

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit SQUALLLINE] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SQUALLLINE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SQUALLLINE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_SQUALLLINE)

    call read_sounding( RHO, VELX, VELY, POTT, QV ) ! (out)

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ( VELX(k) - OFFSET_velx ) * RHO(k)
       MOMY(k,i,j) = ( VELY(k) - OFFSET_vely ) * RHO(k)
       RHOT(k,i,j) = RHO(k) * ( POTT(k) + rndm(k,i,j) * RANDOM_THETA )

       QTRC(k,i,j,I_QV) = QV(k)
    enddo
    enddo
    enddo

    call flux_setup

    call ocean_setup

    return
  end subroutine MKINIT_squallline

  !-----------------------------------------------------------------------------
  !> Make initial state by Weisman and Klemp (1982)
  subroutine MKINIT_wk1982
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA = 300.0_RP   ! surface pot. temperature [K]
    real(RP) :: SFC_PRES               ! surface pressure         [Pa]
    ! Parameter in Weisman and Klemp (1982)
    real(RP) :: TR_Z      = 12000.0_RP ! height           of tropopause  [m]
    real(RP) :: TR_THETA  =   343.0_RP ! pot. temperature at tropopause  [K]
    real(RP) :: TR_TEMP   =   213.0_RP ! temperature      at tropopause  [K]
    real(RP) :: SHEAR_Z   =  3000.0_RP ! center height of shear layer    [m]
    real(RP) :: SHEAR_U   =    15.0_RP ! velocity u over the shear layer [m/s]
    ! Bubble
    real(RP) :: BBL_THETA = 3.D0 ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_WK1982 / &
       SFC_THETA, &
       SFC_PRES,  &
       TR_Z,      &
       TR_THETA,  &
       TR_TEMP,   &
       SHEAR_Z,   &
       SHEAR_U,   &
       BBL_THETA

    real(RP) :: rh    (KA,IA,JA)
    real(RP) :: rh_sfc(1 ,IA,JA)

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit WK1982] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_PRES  = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WK1982,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_WK1982. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_WK1982)

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if ( REAL_CZ(k,i,j) <= TR_Z ) then ! below initial cloud top
             pott(k,i,j) = pott_sfc(1,i,j) &
                         + ( TR_THETA - pott_sfc(1,i,j) ) * ( REAL_CZ(k,i,j) / TR_Z )**1.25_RP
          else
             pott(k,i,j) = TR_THETA * exp( GRAV * ( REAL_CZ(k,i,j) - TR_Z ) / CPdry / TR_TEMP )
          endif

          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    ! calc QV from RH
    do j = JS, JE
    do i = IS, IE
       rh_sfc(1,i,j) = 1.0_RP - 0.75_RP * ( REAL_FZ(KS-1,i,j) / TR_Z )**1.25_RP

       do k = KS, KE
          if ( REAL_CZ(k,i,j) <= TR_Z ) then ! below initial cloud top
             rh(k,i,j) = 1.0_RP - 0.75_RP * ( REAL_CZ(k,i,j) / TR_Z )**1.25_RP
          else
             rh(k,i,j) = 0.25_RP
          endif
       enddo
    enddo
    enddo

    call SATURATION_pres2qsat_all( qsat_sfc(1,:,:), temp_sfc(1,:,:), pres_sfc(1,:,:) )
    call SATURATION_pres2qsat_all( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = rh_sfc(1,i,j) * qsat_sfc(1,i,j)
       do k = KS, KE
          qv(k,i,j) = rh(k,i,j) * qsat(k,i,j)
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

do k = KS, KE
   if( IO_L ) write(IO_FID_LOG,*) k, REAL_CZ(k,IS,JS), pres(k,IS,JS), pott(k,IS,JS), rh(k,IS,JS), qv(k,IS,JS)*1000
enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = SHEAR_U * tanh( REAL_CZ(k,i,j) / SHEAR_Z ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = 0.0_RP
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,i,j) * ( pott(k,i,j) + BBL_THETA * bubble(k,i,j) )

       QTRC(k,i,j,I_QV) = qv(k,i,j)
    enddo
    enddo
    enddo

#ifndef DRY
    if ( QA >= 2 ) then
       do iq = 2, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = 0.0_RP
       enddo
       enddo
       enddo
       enddo
    endif
#endif

    call flux_setup

    return
  end subroutine MKINIT_wk1982

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF01
    implicit none

#ifndef DRY
    real(RP) :: PERTURB_AMP = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0       ! 0 -> no perturbation
                                       ! 1 -> petrurbation for pt
                                       ! 2 -> perturbation for u, v, w
    logical  :: USE_LWSET    = .false. ! use liq. water. static energy temp.?

    NAMELIST / PARAM_MKINIT_RF01 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG,     &
       USE_LWSET

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: sint
    real(RP) :: GEOP_sw ! switch for geopotential energy correction

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP ! pi/2

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit DYCOMS2RF01] / Categ[preprocess] / Origin[SCALE-RM]'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF01,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_RF01. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_RF01)

    if ( USE_LWSET ) then
       GEOP_sw = 1.0_RP
    else
       GEOP_sw = 0.0_RP
    endif

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1017.8E2_RP ! [Pa]
       pott_sfc(1,i,j) = 289.0_RP    ! [K]
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          velx(k,i,j) =   7.0_RP
          vely(k,i,j) =  -5.5_RP
          if ( GRID_CZ(k) < 820.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 289.0_RP - GRAV / CPdry * GRID_CZ(k) * GEOP_sw
          elseif( GRID_CZ(k) <= 860.0_RP ) then
             sint = sin( pi2 * ( GRID_CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             potl(k,i,j) = ( 289.0_RP - GRAV / CPdry * GRID_CZ(k) * GEOP_sw                          ) * (0.5_RP-sint) &
                         + ( 297.5_RP+sign(abs(GRID_CZ(k)-840.0_RP)**(1.0_RP/3.0_RP),GRID_CZ(k)-840.0_RP) &
                           - GRAV / CPdry * GRID_CZ(k) * GEOP_sw                                     ) * (0.5_RP+sint)
          else
             potl(k,i,j) = 297.5_RP + ( GRID_CZ(k)-840.0_RP )**(1.0_RP/3.0_RP) &
                         - GRAV / CPdry * GRID_CZ(k) * GEOP_sw
          endif

          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP
       enddo

    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               potl    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    ! calc in moist condition
    do j = JS, JE
    do i = IS, IE
       qv_sfc  (1,i,j) = 9.0E-3_RP   ! [kg/kg]
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if    ( GRID_CZ(k) <   820.0_RP ) then ! below initial cloud top
             qall = 9.0E-3_RP
          elseif( GRID_CZ(k) <=  860.0_RP ) then ! boundary
             sint = sin( pi2 * ( GRID_CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             qall = 9.0E-3_RP * (0.5_RP-sint) &
                         + 1.5E-3_RP * (0.5_RP+sint)
          elseif( GRID_CZ(k) <= 5000.0_RP ) then
             qall = 1.5E-3_RP
          else
             qall = 0.0_RP
          endif

          if    ( GRID_CZ(k) <=  600.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( GRID_CZ(k) < 820.0_RP ) then ! in the cloud
             fact = ( GRID_CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact
          elseif( GRID_CZ(k) <= 860.0_RP ) then ! boundary
             sint = sin( pi2 * ( GRID_CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             fact = ( GRID_CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact * (0.5_RP-sint)
          else
             qc(k,i,j) = 0.0_RP
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho_bytemp( DENS    (:,:,:), & ! [OUT]
                                      pott    (:,:,:), & ! [OUT]
                                      pres    (:,:,:), & ! [OUT]
                                      temp    (:,:,:), & ! [IN]
                                      qv      (:,:,:), & ! [IN]
                                      qc      (:,:,:), & ! [IN]
                                      pott_sfc(:,:,:), & ! [OUT]
                                      pres_sfc(:,:,:), & ! [IN]
                                      temp_sfc(:,:,:), & ! [IN]
                                      qv_sfc  (:,:,:), & ! [IN]
                                      qc_sfc  (:,:,:)  ) ! [IN]

    do j  = JS, JE
    do i  = IS, IE
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMZ(k,i,j) = ( 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
       else
          MOMZ(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMX(k,i,j) = ( velx(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       else
          MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMY(k,i,j) = ( vely(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       else
          MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          RHOT(k,i,j) = ( pott(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * DENS(k,i,j)
       else
          RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
       endif
    enddo
    enddo
    enddo

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv(k,i,j) + qc(k,i,j) !--- Super saturated air at initial
          do iq = QQA+1, QA
            QTRC(k,i,j,iq) = QTRC(k,i,j,iq) / DENS(k,i,j)
          enddo
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv(k,i,j)
          QTRC(k,i,j,I_QC) = qc(k,i,j)
       enddo
       enddo
       enddo

       if ( I_NC > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if ( qc(k,i,j) > 0.0_RP ) then
                QTRC(k,i,j,I_NC) = 120.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
             endif
          enddo
          enddo
          enddo
       endif
    endif

#endif
    if ( AETRACER_TYPE == 'KAJINO13' ) then
      call AEROSOL_setup
    endif
    return
  end subroutine MKINIT_DYCOMS2_RF01

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF02
    implicit none

#ifndef DRY
    real(RP) :: PERTURB_AMP  = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0 ! 0 -> no perturbation
                                 ! 1 -> perturbation for PT
                                 ! 2 -> perturbation for u,v,w

    NAMELIST / PARAM_MKINIT_RF02 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: sint

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit DYCOMS2RF02] / Categ[preprocess] / Origin[SCALE-RM]'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_RF02. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_RF02)

    ! calc in dry condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1017.8E2_RP   ! [Pa]
       pott_sfc(1,i,j) = 288.3_RP      ! [K]
       qv_sfc(1,i,j) = 0.0_RP
       qc_sfc(1,i,j) = 0.0_RP

       do k = KS, KE
          velx(k,i,j) =  3.0_RP + 4.3 * GRID_CZ(k)*1.E-3_RP
          vely(k,i,j) = -9.0_RP + 5.6 * GRID_CZ(k)*1.E-3_RP

          if ( GRID_CZ(k) < 775.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
          else if ( GRID_CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (GRID_CZ(k) - 795.0_RP)/20.0_RP )
             potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP &
                         + ( 295.0_RP+sign(abs(GRID_CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),GRID_CZ(k)-795.0_RP) ) &
                         * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( GRID_CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
          endif

          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               potl    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    ! calc in moist condition
    do j = JS, JE
    do i = IS, IE
       qv_sfc  (1,i,j) = 9.45E-3_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if ( GRID_CZ(k) < 775.0_RP ) then ! below initial cloud top
             qall = 9.45E-3_RP ! [kg/kg]
          else if ( GRID_CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (GRID_CZ(k) - 795.0_RP)/20.0_RP )
             qall = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
                   ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-GRID_CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             qall = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-GRID_CZ(k))/500.0_RP ) ) ! [kg/kg]
          endif

          if( GRID_CZ(k) < 400.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( GRID_CZ(k) < 775.0_RP ) then
             fact = ( GRID_CZ(k)-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.65E-3_RP * fact
          elseif( GRID_CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * ( GRID_CZ(k)-795.0_RP )/20.0_RP )
             fact = ( GRID_CZ(k)-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.65E-3_RP * fact * (1.0_RP-sint) * 0.5_RP
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho_bytemp( DENS    (:,:,:), & ! [OUT]
                                      pott    (:,:,:), & ! [OUT]
                                      pres    (:,:,:), & ! [OUT]
                                      temp    (:,:,:), & ! [IN]
                                      qv      (:,:,:), & ! [IN]
                                      qc      (:,:,:), & ! [IN]
                                      pott_sfc(:,:,:), & ! [OUT]
                                      pres_sfc(:,:,:), & ! [IN]
                                      temp_sfc(:,:,:), & ! [IN]
                                      qv_sfc  (:,:,:), & ! [IN]
                                      qc_sfc  (:,:,:)  ) ! [IN]

    do j  = JS, JE
    do i  = IS, IE
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMZ(k,i,j) = ( 0.0_RP + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
     else
       MOMZ(k,i,j) = 0.0_RP
     endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMX(k,i,j) = ( velx(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     else
       MOMX(k,i,j) = ( velx(k,i,j) ) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMY(k,i,j) = ( vely(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     else
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then
       RHOT(k,i,j) = ( pott(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * DENS(k,i,j)
     else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
     endif
    enddo
    enddo
    enddo

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          !--- Super saturated air at initial
          QTRC(k,i,j,I_QV) = qv(k,i,j) + qc(k,i,j)
          do iq = QQA+1, QA
            QTRC(k,i,j,iq) = QTRC(k,i,j,iq) / DENS(k,i,j)
          enddo
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv(k,i,j)
          QTRC(k,i,j,I_QC) = qc(k,i,j)
       enddo
       enddo
       enddo

       if ( I_NC > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if ( qc(k,i,j) > 0.0_RP ) then
                QTRC(k,i,j,I_NC) = 55.0E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
             endif
          enddo
          enddo
          enddo
       endif
    endif

#endif
    if ( AETRACER_TYPE == 'KAJINO13' ) then
      call AEROSOL_setup
    endif
    return
  end subroutine MKINIT_DYCOMS2_RF02

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF02_DNS
    implicit none

#ifndef DRY
    real(RP) :: ZB  = 750.0_RP ! domain bottom
!   real(RP) :: ZT  = 900.0_RP ! domain top
    real(RP) :: CONST_U = 0.0_RP
    real(RP) :: CONST_V = 0.0_RP
    real(RP) :: PRES_ZB = 93060.0_RP
    real(RP) :: PERTURB_AMP  = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0 ! 0 -> no perturbation
                                 ! 1 -> perturbation for PT
                                 ! 2 -> perturbation for u,v,w

    NAMELIST / PARAM_MKINIT_RF02_DNS / &
       ZB, CONST_U, CONST_V,PRES_ZB,&
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit DYCOMS2RF02_DNS] / Categ[preprocess] / Origin[SCALE-RM]'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02_DNS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_RF02_DNS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_RF02_DNS)

    ! calc in dry condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = PRES_ZB
!      pott_sfc(1,i,j) = 288.3_RP      ! [K]
!      qv_sfc  (1,i,j) = 9.45E-3_RP
!      qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE

          velx(k,i,j) = CONST_U
          vely(k,i,j) = CONST_V

!         if ( ZB+GRID_CZ(k) < 775.0_RP ) then ! below initial cloud top
          if ( ZB+GRID_CZ(k) <= 795.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
             qall = 9.45E-3_RP ! [kg/kg]
! necessary?
!         else if ( GRID_CZ(k) <= 815.0_RP ) then
!            sint = sin( pi2 * (GRID_CZ(k) - 795.0_RP)/20.0_RP )
!            potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP + &
!                  ( 295.0_RP+sign(abs(GRID_CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),GRID_CZ(k)-795.0_RP) ) * (1.0_RP+sint)*0.5_RP
!            qall = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
!                  ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-GRID_CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( zb+GRID_CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
             qall = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-(zb+GRID_CZ(k)))/500.0_RP ) ) ! [kg/kg]
          endif

          if( ZB+GRID_CZ(k) < 400.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( ZB+GRID_CZ(k) <= 795.0_RP ) then
             fact = ( (zb+GRID_CZ(k))-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.8E-3_RP * fact
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall - qc(k,i,j)

!if(i==is.and.j==js)write(*,*)'chkk',k,cz(k)+zb,qc(k,i,j),qv(k,i,j)
       enddo
    enddo
    enddo

!write(*,*)'chk3',ks,ke
    ! extrapolation (temtative)
    pott_sfc(1,:,:) = potl(ks,:,:)-0.5*(potl(ks+1,:,:)-potl(ks,:,:))
    qv_sfc  (1,:,:) = qv  (ks,:,:)-0.5*(qv  (ks+1,:,:)-qv  (ks,:,:))
    qc_sfc  (1,:,:) = qc  (ks,:,:)-0.5*(qc  (ks+1,:,:)-qc  (ks,:,:))

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               potl    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

!write(*,*)'chk4.1'
    RovCP = Rdry / CPdry
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LHV / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

!write(*,*)'chk5'
    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               pott    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]

    do j  = JS, JE
    do i  = IS, IE
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

!write(*,*)'chk7'
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMZ(k,i,j) = ( 0.0_RP + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
     else
       MOMZ(k,i,j) = 0.0_RP
     endif
    enddo
    enddo
    enddo

!write(*,*)'chk8'
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMX(k,i,j) = ( velx(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     else
       MOMX(k,i,j) = ( velx(k,i,j) ) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo
!write(*,*)'chk9'

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMY(k,i,j) = ( vely(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     else
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo
!write(*,*)'chk10'

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then
       RHOT(k,i,j) = ( pott(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * DENS(k,i,j)
     else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
     endif
    enddo
    enddo
    enddo

!write(*,*)'chk11'
    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo

!write(*,*)'chk12'
    if ( TRACER_TYPE == 'SUZUKI10' ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          !--- Super saturated air at initial
          QTRC(k,i,j,I_QV) = qv(k,i,j) + qc(k,i,j)

          !--- for aerosol
          do iq = QQA+1, QA
             QTRC(k,i,j,iq) = gan(iq-QQA) / DENS(k,i,j)
          enddo
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv(k,i,j)
          QTRC(k,i,j,I_QC) = qc(k,i,j)
       enddo
       enddo
       enddo

       if ( I_NC > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if ( qc(k,i,j) > 0.0_RP ) then
                QTRC(k,i,j,I_NC) = 55.0E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
             endif
          enddo
          enddo
          enddo
       endif
    endif

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF02_DNS

  !-----------------------------------------------------------------------------
  !> Make initial state for RICO inter comparison
  subroutine MKINIT_RICO
    implicit none

#ifndef DRY
    real(RP):: PERTURB_AMP_PT = 0.1_RP
    real(RP):: PERTURB_AMP_QV = 2.5E-5_RP

    NAMELIST / PARAM_MKINIT_RICO / &
       PERTURB_AMP_PT, &
       PERTURB_AMP_QV

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall ! QV+QC
    real(RP) :: fact

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit RICO] / Categ[preprocess] / Origin[SCALE-RM]'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RICO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_RICO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_RICO)

    ! calc in moist condition
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1015.4E2_RP ! [Pa]
       pott_sfc(1,i,j) = 297.9_RP
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          !--- potential temperature
          if ( GRID_CZ(k) < 740.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 297.9_RP
          else
             fact = ( GRID_CZ(k)-740.0_RP ) * ( 317.0_RP-297.9_RP ) / ( 4000.0_RP-740.0_RP )
             potl(k,i,j) = 297.9_RP + fact
          endif

          !--- horizontal wind velocity
          if ( GRID_CZ(k) <= 4000.0_RP ) then ! below initial cloud top
             fact = ( GRID_CZ(k)-0.0_RP ) * ( -1.9_RP+9.9_RP ) / ( 4000.0_RP-0.0_RP )
             velx(k,i,j) =  -9.9_RP + fact
             vely(k,i,j) =  -3.8_RP
          else
             velx(k,i,j) =  -1.9_RP
             vely(k,i,j) =  -3.8_RP
          endif

          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP

       enddo

    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,:,:), & ! [OUT]
                               temp    (:,:,:), & ! [OUT]
                               pres    (:,:,:), & ! [OUT]
                               potl    (:,:,:), & ! [IN]
                               qv      (:,:,:), & ! [IN]
                               qc      (:,:,:), & ! [IN]
                               temp_sfc(:,:,:), & ! [OUT]
                               pres_sfc(:,:,:), & ! [IN]
                               pott_sfc(:,:,:), & ! [IN]
                               qv_sfc  (:,:,:), & ! [IN]
                               qc_sfc  (:,:,:)  ) ! [IN]


    do j = JS, JE
    do i = IS, IE
       qv_sfc  (1,i,j) = 16.0E-3_RP   ! [kg/kg]
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          !--- mixing ratio of vapor
          if ( GRID_CZ(k) <= 740.0_RP ) then ! below initial cloud top
             fact = ( GRID_CZ(k)-0.0_RP ) * ( 13.8E-3_RP-16.0E-3_RP ) / ( 740.0_RP-0.0_RP )
             qall = 16.0E-3_RP + fact
          elseif ( GRID_CZ(k) <= 3260.0_RP ) then ! boundary
             fact = ( GRID_CZ(k)-740.0_RP ) * ( 2.4E-3_RP-13.8E-3_RP ) / ( 3260.0_RP-740.0_RP )
             qall = 13.8E-3_RP + fact
          elseif( GRID_CZ(k) <= 4000.0_RP ) then
             fact = ( GRID_CZ(k)-3260.0_RP ) * ( 1.8E-3_RP-2.4E-3_RP ) / ( 4000.0_RP-3260.0_RP )
             qall = 2.4E-3_RP + fact
          else
             qall = 0.0_RP
          endif

          qc(k,i,j) = 0.0_RP
          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho_bytemp( DENS    (:,:,:), & ! [OUT]
                                      pott    (:,:,:), & ! [OUT]
                                      pres    (:,:,:), & ! [OUT]
                                      temp    (:,:,:), & ! [IN]
                                      qv      (:,:,:), & ! [IN]
                                      qc      (:,:,:), & ! [IN]
                                      pott_sfc(:,:,:), & ! [OUT]
                                      pres_sfc(:,:,:), & ! [IN]
                                      temp_sfc(:,:,:), & ! [IN]
                                      qv_sfc  (:,:,:), & ! [IN]
                                      qc_sfc  (:,:,:)  ) ! [IN]


    do j = JS, JE
    do i = IS, IE
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA  ,i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT(k,i,j) = ( pott(k,i,j)+2.0_RP*( rndm(k,i,j)-0.5_RP )*PERTURB_AMP_PT ) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random

    if ( TRACER_TYPE == 'SUZUKI10' ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          !--- Super saturated air at initial
          QTRC(k,i,j,I_QV) = qv(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP_QV &
                           + qc(k,i,j)
          do iq = QQA+1, QA
            QTRC(k,i,j,iq) = QTRC(k,i,j,iq) / DENS(k,i,j)
          enddo
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP_QV
          QTRC(k,i,j,I_QC) = qc(k,i,j)
       enddo
       enddo
       enddo

       if ( I_NC > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if ( qc(k,i,j) > 0.0_RP ) then
                QTRC(k,i,j,I_NC) = 70.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
             endif
          enddo
          enddo
          enddo
       endif
    endif

#endif
    if ( AETRACER_TYPE == 'KAJINO13' ) then
      call AEROSOL_setup
    endif
    return
  end subroutine MKINIT_RICO

  !-----------------------------------------------------------------------------
  subroutine MKINIT_interporation
    use gtool_file, only: &
       FileGetShape, &
       FileRead
    implicit none

    real(RP) :: dz(KA,IA,JA)

    real(RP) :: W(KA,IA,JA)
    real(RP) :: U(KA,IA,JA)
    real(RP) :: V(KA,IA,JA)

    real(RP) :: fact_cz0(KA)
    real(RP) :: fact_cz1(KA)
    real(RP) :: fact_fz0(KA)
    real(RP) :: fact_fz1(KA)
    real(RP) :: fact_cx0(IA)
    real(RP) :: fact_cx1(IA)
    real(RP) :: fact_fx0(IA)
    real(RP) :: fact_fx1(IA)
    real(RP) :: fact_cy0(JA)
    real(RP) :: fact_cy1(JA)
    real(RP) :: fact_fy0(JA)
    real(RP) :: fact_fy1(JA)

    integer :: idx_cz0(KA)
    integer :: idx_cz1(KA)
    integer :: idx_fz0(KA)
    integer :: idx_fz1(KA)
    integer :: idx_cx0(IA)
    integer :: idx_cx1(IA)
    integer :: idx_fx0(IA)
    integer :: idx_fx1(IA)
    integer :: idx_cy0(JA)
    integer :: idx_cy1(JA)
    integer :: idx_fy0(JA)
    integer :: idx_fy1(JA)

    real(RP), allocatable :: DENS_ORG(:,:,:)
    real(RP), allocatable :: MOMZ_ORG(:,:,:)
    real(RP), allocatable :: MOMX_ORG(:,:,:)
    real(RP), allocatable :: MOMY_ORG(:,:,:)
    real(RP), allocatable :: RHOT_ORG(:,:,:)
    real(RP), allocatable :: QTRC_ORG(:,:,:,:)

    real(RP), allocatable :: W_ORG(:,:,:)
    real(RP), allocatable :: U_ORG(:,:,:)
    real(RP), allocatable :: V_ORG(:,:,:)
    real(RP), allocatable :: POTT_ORG(:,:,:)

    real(RP), allocatable :: CZ_ORG(:)
    real(RP), allocatable :: FZ_ORG(:)
    real(RP), allocatable :: CX_ORG(:)
    real(RP), allocatable :: FX_ORG(:)
    real(RP), allocatable :: CY_ORG(:)
    real(RP), allocatable :: FY_ORG(:)

    integer :: dims(3)

    character(len=H_LONG) :: BASENAME_ORG = ''

    NAMELIST / PARAM_MKINIT_INTERPORATION / &
         BASENAME_ORG

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit INTERPORATION] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_INTERPORATION,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_INTERPORATION. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_INTERPORATION)

    call FileGetShape( dims(:),                               &
                       BASENAME_ORG, "DENS", 1, single=.true. )

    allocate( dens_org(dims(1),dims(2),dims(3)) )
    allocate( momz_org(dims(1),dims(2),dims(3)) )
    allocate( momx_org(dims(1),dims(2),dims(3)) )
    allocate( momy_org(dims(1),dims(2),dims(3)) )
    allocate( rhot_org(dims(1),dims(2),dims(3)) )
    allocate( qtrc_org(dims(1),dims(2),dims(3),QA) )

    allocate( w_org(dims(1),dims(2),dims(3)) )
    allocate( u_org(dims(1),dims(2),dims(3)) )
    allocate( v_org(dims(1),dims(2),dims(3)) )
    allocate( pott_org(dims(1),dims(2),dims(3)) )

    allocate( cz_org(dims(1)) )
    allocate( fz_org(dims(1)) )
    allocate( cx_org(dims(2)) )
    allocate( fx_org(dims(2)) )
    allocate( cy_org(dims(3)) )
    allocate( fy_org(dims(3)) )

    call FileRead( dens_org(:,:,:),                          &
                   BASENAME_ORG, "DENS", 1, 1, single=.true. )
    call FileRead( momz_org(:,:,:),                          &
                   BASENAME_ORG, "MOMZ", 1, 1, single=.true. )
    call FileRead( momx_org(:,:,:),                          &
                   BASENAME_ORG, "MOMX", 1, 1, single=.true. )
    call FileRead( momy_org(:,:,:),                          &
                   BASENAME_ORG, "MOMY", 1, 1, single=.true. )
    call FileRead( rhot_org(:,:,:),                          &
                   BASENAME_ORG, "RHOT", 1, 1, single=.true. )
    do iq = 1, QA
       call FileRead( qtrc_org(:,:,:,iq),                            &
                      BASENAME_ORG, AQ_NAME(iq), 1, 1, single=.true. )
    end do

    call FileRead( cz_org(:),                              &
                   BASENAME_ORG, "z" , 1, 1, single=.true. )
    call FileRead( cx_org(:),                              &
                   BASENAME_ORG, "x" , 1, 1, single=.true. )
    call FileRead( cy_org(:),                              &
                   BASENAME_ORG, "y" , 1, 1, single=.true. )
    call FileRead( fx_org(:),                              &
                   BASENAME_ORG, "xh", 1, 1, single=.true. )
    call FileRead( fy_org(:),                              &
                   BASENAME_ORG, "yh", 1, 1, single=.true. )

    do k = KS, KE
       call interporation_fact( fact_cz0(k), fact_cz1(k),     & ! (OUT)
                                idx_cz0( k), idx_cz1 (k),     & ! (OUT)
                                GRID_CZ (k), cz_org, dims(1), & ! (IN)
                                .false.                       ) ! (IN)
       call interporation_fact( fact_fz0(k), fact_fz1(k),     & ! (OUT)
                                idx_fz0 (k), idx_fz1 (k),     & ! (OUT)
                                GRID_FZ (k), fz_org, dims(1), & ! (IN)
                                .false.                       ) ! (IN)
    enddo
    do i = IS, IE
       call interporation_fact( fact_cx0(i), fact_cx1(i),     & ! (OUT)
                                idx_cx0 (i), idx_cx1 (i),     & ! (OUT)
                                GRID_CX (i), cx_org, dims(2), & ! (IN)
                                .true.                        ) ! (IN)
       call interporation_fact( fact_fx0(i), fact_fx1(i),     & ! (OUT)
                                idx_fx0 (i), idx_fx1 (i),     & ! (OUT)
                                GRID_FX (i), fx_org, dims(2), & ! (IN)
                                .true.                        ) ! (IN)
    enddo
    do j = JS, JE
       call interporation_fact( fact_cy0(j), fact_cy1(j),     & ! (OUT)
                                idx_cy0 (j), idx_cy1 (j),     & ! (OUT)
                                GRID_CY (j), cy_org, dims(3), & ! (IN)
                                .true.                        ) ! (IN)
       call interporation_fact( fact_fy0(j), fact_fy1(j),     & ! (OUT)
                                idx_fy0 (j), idx_fy1 (j),     & ! (OUT)
                                GRID_FY (j), fy_org, dims(3), & ! (IN)
                                .true.                        ) ! (IN)
    enddo


    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)-1
       w_org(k,i,j) = 2.0_RP * momz_org(k,i,j) / ( dens_org(k+1,i,j) + dens_org(k,i,j) )
    end do
    end do
    end do
    do j = 1, dims(3)
    do i = 1, dims(2)
       w_org(dims(1),i,j) = 0.0_RP
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)-1
    do k = 1, dims(1)
       u_org(k,i,j) = 2.0_RP * momx_org(k,i,j) / ( dens_org(k,i+1,j) + dens_org(k,i,j) )
    end do
    end do
    end do
    do j = 1, dims(3)
    do k = 1, dims(1)
       u_org(k,dims(2),j) = 2.0_RP * momx_org(k,dims(2),j) / ( dens_org(k,1,j) + dens_org(k,dims(2),j) )
    end do
    end do

    do j = 1, dims(3)-1
    do i = 1, dims(2)
    do k = 1, dims(1)
       v_org(k,i,j) = 2.0_RP * momy_org(k,i,j) / ( dens_org(k,i,j+1) + dens_org(k,i,j) )
    end do
    end do
    end do
    do i = 1, dims(2)
    do k = 1, dims(1)
       v_org(k,i,dims(3)) = 2.0_RP * momy_org(k,i,dims(3)) / ( dens_org(k,i,1) + dens_org(k,i,dims(3)) )
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       pott_org(k,i,j) = rhot_org(k,i,j) / dens_org(k,i,j)
    end do
    end do
    end do

    do j = 1, dims(3)
    do i = 1, dims(2)
    do k = 1, dims(1)
       dens_org(k,i,j) = log( dens_org(k,i,j) )
    end do
    end do
    end do

    do j = JS, JE
    do i = IS, IE
    do k = KS, IE
       DENS(k,i,j) = exp( &
                     fact_cz0(k)*fact_cx0(i)*fact_cy0(j)*dens_org(idx_cz0(k),idx_cx0(i),idx_cy0(j)) &
                   + fact_cz1(k)*fact_cx0(i)*fact_cy0(j)*dens_org(idx_cz1(k),idx_cx0(i),idx_cy0(j)) &
                   + fact_cz0(k)*fact_cx1(i)*fact_cy0(j)*dens_org(idx_cz0(k),idx_cx1(i),idx_cy0(j)) &
                   + fact_cz1(k)*fact_cx1(i)*fact_cy0(j)*dens_org(idx_cz1(k),idx_cx1(i),idx_cy0(j)) &
                   + fact_cz0(k)*fact_cx0(i)*fact_cy1(j)*dens_org(idx_cz0(k),idx_cx0(i),idx_cy1(j)) &
                   + fact_cz1(k)*fact_cx0(i)*fact_cy1(j)*dens_org(idx_cz1(k),idx_cx0(i),idx_cy1(j)) &
                   + fact_cz0(k)*fact_cx1(i)*fact_cy1(j)*dens_org(idx_cz0(k),idx_cx1(i),idx_cy1(j)) &
                   + fact_cz1(k)*fact_cx1(i)*fact_cy1(j)*dens_org(idx_cz1(k),idx_cx1(i),idx_cy1(j)) &
                   )

       W(k,i,j) = fact_fz0(k)*fact_cx0(i)*fact_cy0(j)*w_org(idx_fz0(k),idx_cx0(i),idx_cy0(j)) &
                + fact_fz1(k)*fact_cx0(i)*fact_cy0(j)*w_org(idx_fz1(k),idx_cx0(i),idx_cy0(j)) &
                + fact_fz0(k)*fact_cx1(i)*fact_cy0(j)*w_org(idx_fz0(k),idx_cx1(i),idx_cy0(j)) &
                + fact_fz1(k)*fact_cx1(i)*fact_cy0(j)*w_org(idx_fz1(k),idx_cx1(i),idx_cy0(j)) &
                + fact_fz0(k)*fact_cx0(i)*fact_cy1(j)*w_org(idx_fz0(k),idx_cx0(i),idx_cy1(j)) &
                + fact_fz1(k)*fact_cx0(i)*fact_cy1(j)*w_org(idx_fz1(k),idx_cx0(i),idx_cy1(j)) &
                + fact_fz0(k)*fact_cx1(i)*fact_cy1(j)*w_org(idx_fz0(k),idx_cx1(i),idx_cy1(j)) &
                + fact_fz1(k)*fact_cx1(i)*fact_cy1(j)*w_org(idx_fz1(k),idx_cx1(i),idx_cy1(j))

       U(k,i,j) = fact_cz0(k)*fact_fx0(i)*fact_cy0(j)*u_org(idx_cz0(k),idx_fx0(i),idx_cy0(j)) &
                + fact_cz1(k)*fact_fx0(i)*fact_cy0(j)*u_org(idx_cz1(k),idx_fx0(i),idx_cy0(j)) &
                + fact_cz0(k)*fact_fx1(i)*fact_cy0(j)*u_org(idx_cz0(k),idx_fx1(i),idx_cy0(j)) &
                + fact_cz1(k)*fact_fx1(i)*fact_cy0(j)*u_org(idx_cz1(k),idx_fx1(i),idx_cy0(j)) &
                + fact_cz0(k)*fact_fx0(i)*fact_cy1(j)*u_org(idx_cz0(k),idx_fx0(i),idx_cy1(j)) &
                + fact_cz1(k)*fact_fx0(i)*fact_cy1(j)*u_org(idx_cz1(k),idx_fx0(i),idx_cy1(j)) &
                + fact_cz0(k)*fact_fx1(i)*fact_cy1(j)*u_org(idx_cz0(k),idx_fx1(i),idx_cy1(j)) &
                + fact_cz1(k)*fact_fx1(i)*fact_cy1(j)*u_org(idx_cz1(k),idx_fx1(i),idx_cy1(j))

       V(k,i,j) = fact_cz0(k)*fact_cx0(i)*fact_fy0(j)*v_org(idx_cz0(k),idx_cx0(i),idx_fy0(j)) &
                + fact_cz1(k)*fact_cx0(i)*fact_fy0(j)*v_org(idx_cz1(k),idx_cx0(i),idx_fy0(j)) &
                + fact_cz0(k)*fact_cx1(i)*fact_fy0(j)*v_org(idx_cz0(k),idx_cx1(i),idx_fy0(j)) &
                + fact_cz1(k)*fact_cx1(i)*fact_fy0(j)*v_org(idx_cz1(k),idx_cx1(i),idx_fy0(j)) &
                + fact_cz0(k)*fact_cx0(i)*fact_fy1(j)*v_org(idx_cz0(k),idx_cx0(i),idx_fy1(j)) &
                + fact_cz1(k)*fact_cx0(i)*fact_fy1(j)*v_org(idx_cz1(k),idx_cx0(i),idx_fy1(j)) &
                + fact_cz0(k)*fact_cx1(i)*fact_fy1(j)*v_org(idx_cz0(k),idx_cx1(i),idx_fy1(j)) &
                + fact_cz1(k)*fact_cx1(i)*fact_fy1(j)*v_org(idx_cz1(k),idx_cx1(i),idx_fy1(j))

       POTT(k,i,j) = fact_cz0(k)*fact_cx0(i)*fact_cy0(j)*pott_org(idx_cz0(k),idx_cx0(i),idx_cy0(j)) &
                   + fact_cz1(k)*fact_cx0(i)*fact_cy0(j)*pott_org(idx_cz1(k),idx_cx0(i),idx_cy0(j)) &
                   + fact_cz0(k)*fact_cx1(i)*fact_cy0(j)*pott_org(idx_cz0(k),idx_cx1(i),idx_cy0(j)) &
                   + fact_cz1(k)*fact_cx1(i)*fact_cy0(j)*pott_org(idx_cz1(k),idx_cx1(i),idx_cy0(j)) &
                   + fact_cz0(k)*fact_cx0(i)*fact_cy1(j)*pott_org(idx_cz0(k),idx_cx0(i),idx_cy1(j)) &
                   + fact_cz1(k)*fact_cx0(i)*fact_cy1(j)*pott_org(idx_cz1(k),idx_cx0(i),idx_cy1(j)) &
                   + fact_cz0(k)*fact_cx1(i)*fact_cy1(j)*pott_org(idx_cz0(k),idx_cx1(i),idx_cy1(j)) &
                   + fact_cz1(k)*fact_cx1(i)*fact_cy1(j)*pott_org(idx_cz1(k),idx_cx1(i),idx_cy1(j))

       do iq = 1, QA
          QTRC(k,i,j,iq) = fact_cz0(k)*fact_cx0(i)*fact_cy0(j)*qtrc_org(idx_cz0(k),idx_cx0(i),idx_cy0(j),iq) &
                         + fact_cz1(k)*fact_cx0(i)*fact_cy0(j)*qtrc_org(idx_cz1(k),idx_cx0(i),idx_cy0(j),iq) &
                         + fact_cz0(k)*fact_cx1(i)*fact_cy0(j)*qtrc_org(idx_cz0(k),idx_cx1(i),idx_cy0(j),iq) &
                         + fact_cz1(k)*fact_cx1(i)*fact_cy0(j)*qtrc_org(idx_cz1(k),idx_cx1(i),idx_cy0(j),iq) &
                         + fact_cz0(k)*fact_cx0(i)*fact_cy1(j)*qtrc_org(idx_cz0(k),idx_cx0(i),idx_cy1(j),iq) &
                         + fact_cz1(k)*fact_cx0(i)*fact_cy1(j)*qtrc_org(idx_cz1(k),idx_cx0(i),idx_cy1(j),iq) &
                         + fact_cz0(k)*fact_cx1(i)*fact_cy1(j)*qtrc_org(idx_cz0(k),idx_cx1(i),idx_cy1(j),iq) &
                         + fact_cz1(k)*fact_cx1(i)*fact_cy1(j)*qtrc_org(idx_cz1(k),idx_cx1(i),idx_cy1(j),iq)
          enddo
    enddo
    enddo
    enddo

    deallocate( dens_org )
    deallocate( momz_org )
    deallocate( momx_org )
    deallocate( momy_org )
    deallocate( rhot_org )
    deallocate( qtrc_org )

    deallocate( w_org )
    deallocate( u_org )
    deallocate( v_org )
    deallocate( pott_org )

    deallocate( cz_org )
    deallocate( fz_org )
    deallocate( cx_org )
    deallocate( fx_org )
    deallocate( cy_org )
    deallocate( fy_org )

    if ( I_QV > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qv(k,i,j) = QTRC(k,i,j,I_QV)
       end do
       end do
       end do
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qv(k,i,j) = 0.0_RP
       end do
       end do
       end do
    end if

    if ( I_QC > 0 ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qc(k,i,j) = QTRC(k,i,j,I_QC)
       end do
       end do
       end do
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qc(k,i,j) = 0.0_RP
       end do
       end do
       end do
    end if

    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE
       dz(k,i,j) = REAL_CZ(k,i,j) - REAL_CZ(k-1,i,j) ! distance from cell center to cell center
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho_atmos( DENS(:,:,:), & ! [INOUT]
                                     temp(:,:,:), & ! [OUT]
                                     pres(:,:,:), & ! [OUT]
                                     pott(:,:,:), & ! [IN]
                                     qv  (:,:,:), & ! [IN]
                                     qc  (:,:,:), & ! [IN]
                                     dz  (:,:,:)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.5_RP * W(k,i,j) * ( DENS(k,i,j) + DENS(k+1,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = 0.5_RP * U(k,i,j) * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = 0.5_RP * V(k,i,j) * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_interporation

  !-----------------------------------------------------------------------------
  subroutine interporation_fact( &
       fact0, fact1, &
       idx0, idx1,   &
       x, x_org, nx, &
       loop          )
    implicit none

    real(RP), intent(out) :: fact0
    real(RP), intent(out) :: fact1
    integer,  intent(out) :: idx0
    integer,  intent(out) :: idx1
    real(RP), intent(in)  :: x
    integer,  intent(in)  :: nx
    real(RP), intent(in)  :: x_org(nx)
    logical,  intent(in)  :: loop

    real(RP) :: xwork
    integer :: i

    if ( x < x_org(1) ) then
       if ( loop ) then
          xwork = x_org(1) - ( x_org(2) - x_org(1) )**2 / ( x_org(3) - x_org(2) )
          fact0 = ( x_org(1) - x ) / ( x_org(1) - xwork )
          fact1 = ( x - xwork )    / ( x_org(1) - xwork )
          idx0 = nx
          idx1 = 1
       else
          fact0 = ( x_org(2) - x ) / ( x_org(2) - x_org(1) )
          fact1 = ( x - x_org(1) ) / ( x_org(2) - x_org(1) )
          idx0 = 1
          idx1 = 2
       end if
    else if ( x > x_org(nx) ) then
       if ( loop ) then
          xwork = x_org(nx) + ( x_org(nx) - x_org(nx-1) )**2 / ( x_org(nx-1) - x_org(nx-2) )
          fact0 = ( xwork - x )     / ( xwork - x_org(nx) )
          fact1 = ( x - x_org(nx) ) / ( xwork - x_org(nx) )
          idx0 = nx
          idx1 = 1
       else
          fact0 = ( x_org(nx) - x )   / ( x_org(nx) - x_org(nx-1) )
          fact1 = ( x - x_org(nx-1) ) / ( x_org(nx) - x_org(nx-1) )
          idx0 = nx-1
          idx1 = nx
       end if
    else
       do i = 2, nx
          if ( x <= x_org(i) ) then
             fact0 = ( x_org(i) - x )   / ( x_org(i) - x_org(i-1) )
             fact1 = ( x - x_org(i-1) ) / ( x_org(i) - x_org(i-1) )
             idx0 = i-1
             idx1 = i
             exit
          end if
       end do
    end if

    return
  end subroutine interporation_fact
  !-----------------------------------------------------------------------------
  !> Make initial state ( ocean variables )
  subroutine MKINIT_oceancouple
    implicit none

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit OceanCouple] / Categ[preprocess] / Origin[SCALE-RM]'

    call flux_setup

    call ocean_setup

    return
  end subroutine MKINIT_oceancouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( land variables )
  subroutine MKINIT_landcouple
    implicit none

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit LandCouple] / Categ[preprocess] / Origin[SCALE-RM]'

    call flux_setup

    call land_setup

    return
  end subroutine MKINIT_landcouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( urban variables )
  subroutine MKINIT_urbancouple
    implicit none

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit UrbanCouple] / Categ[preprocess] / Origin[SCALE-RM]'

    call flux_setup

    call urban_setup

    return
  end subroutine MKINIT_urbancouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( sea breeze )
  subroutine MKINIT_seabreeze
    use scale_rm_process, only: &
       PRC_NUM_X
    use scale_landuse, only: &
       LANDUSE_calc_fact, &
       LANDUSE_frac_land
    implicit none

    real(RP) :: dist

    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit SEABREEZE] / Categ[preprocess] / Origin[SCALE-RM]'

    call flux_setup

    call land_setup

    call ocean_setup

    ! quartor size of domain
    dist = ( GRID_CXG(IMAX*PRC_NUM_X) - GRID_CXG(1) ) / 8.0_RP

    ! make landuse conditions
    do j = JS, JE
    do i = IS, IE
       if (       GRID_CX(i) >= dist * 3.0_RP &
            .AND. GRID_CX(i) <= dist * 5.0_RP ) then
          LANDUSE_frac_land(i,j) = 1.0_RP
       else
          LANDUSE_frac_land(i,j) = 0.0_RP
       endif
    enddo
    enddo

    call LANDUSE_calc_fact

    return
  end subroutine MKINIT_seabreeze

  !-----------------------------------------------------------------------------
  !> Make initial state ( heat island )
  subroutine MKINIT_heatisland
    use scale_rm_process, only: &
       PRC_NUM_X
    use scale_landuse, only: &
       LANDUSE_calc_fact, &
       LANDUSE_frac_land,    &
       LANDUSE_frac_urban
    implicit none

    real(RP) :: dist

    integer  :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit HEATISLAND] / Categ[preprocess] / Origin[SCALE-RM]'

    call flux_setup

    call land_setup

    call urban_setup

    ! 1/9 size of domain
    dist = ( GRID_CXG(IMAX*PRC_NUM_X) - GRID_CXG(1) ) / 9.0_RP

    ! make landuse conditions
    do j = JS, JE
    do i = IS, IE
       if (       GRID_CX(i) >= dist * 4.0_RP &
            .AND. GRID_CX(i) <  dist * 5.0_RP ) then
          LANDUSE_frac_land(i,j)  = 1.0_RP
          LANDUSE_frac_urban(i,j) = 1.0_RP
       else
          LANDUSE_frac_land(i,j)  = 1.0_RP
          LANDUSE_frac_urban(i,j) = 0.0_RP
       endif
    enddo
    enddo

    call LANDUSE_calc_fact

    return
  end subroutine MKINIT_heatisland

  !-----------------------------------------------------------------------------
  !> Make initial state for grayzone experiment
  subroutine MKINIT_grayzone
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV(KA)

    ! Bubble
    real(RP) :: PERTURB_AMP = 0.0_RP
    integer  :: RANDOM_LIMIT = 0
    integer  :: RANDOM_FLAG  = 0       ! 0 -> no perturbation
                                       ! 1 -> petrurbation for pt
                                       ! 2 -> perturbation for u, v, w

    NAMELIST / PARAM_MKINIT_GRAYZONE / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit GRAYZONE] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_GRAYZONE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_GRAYZONE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_GRAYZONE)

    call read_sounding( RHO, VELX, VELY, POTT, QV ) ! (out)

!   do j = JS, JE
!   do i = IS, IE
    do j = 1, ja
    do i = 1, ia
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
!      MOMZ(k,i,j) = 0.0_RP
!      MOMX(k,i,j) = RHO(k) * VELX(k)
!      MOMY(k,i,j) = RHO(k) * VELY(k)

!      RHOT(k,i,j) = RHO(k) * POTT(k)

       QTRC(k,i,j,I_QV) = QV(k)
    enddo
    enddo
    enddo

    do j  = JS, JE
    do i  = IS, IE
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMZ(k,i,j) = ( 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
       else
          MOMZ(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMX(k,i,j) = ( velx(k) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       else
          MOMX(k,i,j) = velx(k) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMY(k,i,j) = ( vely(k) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       else
          MOMY(k,i,j) = vely(k) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          RHOT(k,i,j) = ( pott(k) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                      * DENS(k,i,j)
       else
          RHOT(k,i,j) = pott(k) * DENS(k,i,j)
       endif
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_grayzone
  !-----------------------------------------------------------------------------
  !> Make initial state of Box model experiment for zerochemical module
  subroutine MKINIT_boxaero

    implicit none

    real(RP) :: init_dens  = 1.12_RP   ![kg/m3]
    real(RP) :: init_temp  = 298.18_RP ![K]
    real(RP) :: init_pres  = 1.E+5_RP  ![Pa]
    real(RP) :: init_ssliq = 0.01_RP   ![%]

    NAMELIST / PARAM_MKINIT_BOXAERO / &
         init_dens, &
         init_temp, &
         init_pres, &
         init_ssliq

    real(RP) :: qsat
    integer :: i, j, k, ierr

    if ( AETRACER_TYPE /= 'KAJINO13' ) then
       if( IO_L ) write(IO_FID_LOG,*) '+++ For [Box model of aerosol],'
       if( IO_L ) write(IO_FID_LOG,*) '+++ AETRACER_TYPE should be KAJINO13. Stop!'
       call PRC_MPIstop
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit Box model of aerosol] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BOXAERO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_BOXAERO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_BOXAERO)

    QTRC(:,:,:,:) = 0.0_RP
    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
      DENS(k,i,j) = init_dens
      MOMX(k,i,j) = 0.0_RP
      MOMY(k,i,j) = 0.0_RP
      MOMZ(k,i,j) = 0.0_RP
      pott(k,i,j) = init_temp * ( P00/init_pres )**( Rdry/CPdry )
      RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)
      call SATURATION_pres2qsat_all( qsat,init_temp,init_pres )
      QTRC(k,i,j,I_QV) = ( init_ssliq + 1.0_RP )*qsat
    enddo
    enddo
    enddo

    if( AETRACER_TYPE /= 'KAJINO13' ) then
      call AEROSOL_setup
    endif

  end subroutine MKINIT_boxaero
  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  subroutine MKINIT_warmbubbleaero
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  80.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =  80.0_RP ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =  1.E3_RP ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  = 14.E3_RP ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS = 4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS = 3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(RP) :: BBL_THETA    =   1.0_RP ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[mkinit WARMBUBBLEAERO] / Categ[preprocess] / Origin[SCALE-RM]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_WARMBUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       if    ( GRID_CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( GRID_CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( GRID_CZ(k)-GRID_CZ(k-1) )
       else                                ! Layer 3
          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( GRID_CZ(k)-GRID_CZ(k-1) )
       endif
       qv(k,1,1) = 0.0_RP
       qc(k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( qsat_sfc(1,1,1), temp_sfc(1,1,1), pres_sfc(1,1,1) )
    call SATURATION_pres2qsat_all( qsat    (:,1,1), temp    (:,1,1), pres    (:,1,1) )

    qv_sfc(1,1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1,1)
    do k = KS, KE
       if    ( GRID_CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       elseif( GRID_CZ(k) <= ENV_L2_ZTOP ) then ! Layer 2
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       else                                ! Layer 3
          qv(k,1,1) = 0.0_RP
       endif
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( DENS    (:,1,1), & ! [OUT]
                               temp    (:,1,1), & ! [OUT]
                               pres    (:,1,1), & ! [OUT]
                               pott    (:,1,1), & ! [IN]
                               qv      (:,1,1), & ! [IN]
                               qc      (:,1,1), & ! [IN]
                               temp_sfc(1,1,1), & ! [OUT]
                               pres_sfc(1,1,1), & ! [IN]
                               pott_sfc(1,1,1), & ! [IN]
                               qv_sfc  (1,1,1), & ! [IN]
                               qc_sfc  (1,1,1)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       QTRC(k,i,j,I_QV) = qv(k,1,1)
    enddo
    enddo
    enddo

    call flux_setup

    if ( AETRACER_TYPE /= 'KAJINO13' ) then
      call AEROSOL_setup
    endif

    return
  end subroutine MKINIT_warmbubbleaero

  !-----------------------------------------------------------------------------
  !> Make initial state ( real case )
  subroutine MKINIT_real
    use mod_realinput, only: &
         REALINPUT_atmos, &
         REALINPUT_surface
    implicit none

    call REALINPUT_atmos( flg_intrp )

    call REALINPUT_surface

    call flux_setup

    return
  end subroutine MKINIT_real

end module mod_mkinit
!-------------------------------------------------------------------------------
