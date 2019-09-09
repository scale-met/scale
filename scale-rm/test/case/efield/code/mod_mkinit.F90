!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mkinit
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index

  use scale_prc, only: &
     PRC_abort
  use scale_const, only: &
     PI     => CONST_PI,     &
     GRAV   => CONST_GRAV,   &
     Pstd   => CONST_Pstd,   &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     P00    => CONST_PRE00,  &
     EPSvap => CONST_EPSvap
  use scale_random, only: &
     RANDOM_uniform
  use scale_comm_cartesC, only: &
     COMM_vars8, &
     COMM_wait
  use scale_atmos_grid_cartesC, only: &
     CZ  => ATMOS_GRID_CARTESC_CZ,  &
     CX  => ATMOS_GRID_CARTESC_CX,  &
     CY  => ATMOS_GRID_CARTESC_CY,  &
     FZ  => ATMOS_GRID_CARTESC_FZ,  &
     FX  => ATMOS_GRID_CARTESC_FX,  &
     FY  => ATMOS_GRID_CARTESC_FY,  &
     CXG => ATMOS_GRID_CARTESC_CXG, &
     FXG => ATMOS_GRID_CARTESC_FXG, &
     FYG => ATMOS_GRID_CARTESC_FYG
  use scale_atmos_grid_cartesC_real, only: &
     REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
     REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
     AREA    => ATMOS_GRID_CARTESC_REAL_AREA
  use scale_atmos_profile, only: &
     PROFILE_isa => ATMOS_PROFILE_isa
  use scale_atmos_hydrometeor, only: &
     HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
  use scale_atmos_hydrostatic, only: &
     HYDROSTATIC_buildrho        => ATMOS_HYDROSTATIC_buildrho,       &
     HYDROSTATIC_buildrho_atmos  => ATMOS_HYDROSTATIC_buildrho_atmos, &
     HYDROSTATIC_buildrho_bytemp => ATMOS_HYDROSTATIC_buildrho_bytemp
  use scale_atmos_saturation, only: &
     SATURATION_psat_all => ATMOS_SATURATION_psat_all, &
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
  integer, public, parameter :: I_BAROCWAVE        = 30
  integer, public, parameter :: I_BOMEX            = 31
  integer, public, parameter :: I_POINTE           = 32

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
  private :: MKINIT_barocwave

  private :: MKINIT_warmbubble
  private :: MKINIT_supercell
  private :: MKINIT_squallline
  private :: MKINIT_wk1982
  private :: MKINIT_DYCOMS2_RF01
  private :: MKINIT_DYCOMS2_RF02
  private :: MKINIT_RICO
  private :: MKINIT_BOMEX
  private :: MKINIT_POINTE

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
  real(RP), private, allocatable         :: qdry    (:,:,:) ! dry air mass ratio [kg/kg]
  real(RP), private, allocatable         :: qsat    (:,:,:) ! satulated water vapor [kg/kg]
  real(RP), private, allocatable         :: qv      (:,:,:) ! water vapor [kg/kg]
  real(RP), private, allocatable         :: qc      (:,:,:) ! cloud water [kg/kg]
  real(RP), private, allocatable         :: nc      (:,:,:) ! cloud water number density [1/kg]
  real(RP), private, allocatable         :: velx    (:,:,:) ! velocity u [m/s]
  real(RP), private, allocatable         :: vely    (:,:,:) ! velocity v [m/s]

  real(RP), private, allocatable         :: pres_sfc(:,:) ! surface pressure [Pa]
  real(RP), private, allocatable         :: temp_sfc(:,:) ! surface temperature [K]
  real(RP), private, allocatable         :: pott_sfc(:,:) ! surface potential temperature [K]
  real(RP), private, allocatable         :: psat_sfc(:,:) ! surface satulated water pressure [Pa]
  real(RP), private, allocatable         :: qsat_sfc(:,:) ! surface satulated water vapor [kg/kg]
  real(RP), private, allocatable         :: qv_sfc  (:,:) ! surface water vapor [kg/kg]
  real(RP), private, allocatable         :: qc_sfc  (:,:) ! surface cloud water [kg/kg]

  real(RP), private, allocatable         :: rndm    (:,:,:) ! random    number (0-1)
  real(RP), private, allocatable, target :: bubble  (:,:,:) ! bubble    factor (0-1)
  real(RP), private, allocatable, target :: rect    (:,:,:) ! rectangle factor (0-1)
  real(RP), private, allocatable         :: gan     (:)     ! gamma     factor (0-1)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKINIT_setup
    implicit none

    character(len=H_SHORT) :: MKINIT_initname = 'NONE'

    namelist / PARAM_MKINIT / &
       MKINIT_initname

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_setup",*) 'Not appropriate names in namelist PARAM_MKINIT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT)

    allocate( pres(KA,IA,JA) )
    allocate( temp(KA,IA,JA) )
    allocate( pott(KA,IA,JA) )
    allocate( qdry(KA,IA,JA) )
    allocate( qsat(KA,IA,JA) )
    allocate( qv  (KA,IA,JA) )
    allocate( qc  (KA,IA,JA) )
    allocate( nc  (KA,IA,JA) )
    allocate( velx(KA,IA,JA) )
    allocate( vely(KA,IA,JA) )

    allocate( pres_sfc(IA,JA) )
    allocate( temp_sfc(IA,JA) )
    allocate( pott_sfc(IA,JA) )
    allocate( psat_sfc(IA,JA) )
    allocate( qsat_sfc(IA,JA) )
    allocate( qv_sfc  (IA,JA) )
    allocate( qc_sfc  (IA,JA) )

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
    case('BOMEX')
       MKINIT_TYPE = I_BOMEX
    case('POINTE')
       MKINIT_TYPE = I_POINTE
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
    case('BAROCWAVE')
       MKINIT_TYPE = I_BAROCWAVE
    case default
       LOG_ERROR("MKINIT_setup",*) 'Unsupported TYPE:', trim(MKINIT_initname)
       call PRC_abort
    endselect

    return
  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKINIT( output )
    use scale_const, only: &
       CONST_UNDEF8
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    implicit none
    logical, intent(out) :: output

    real(RP) :: QHYD(KA,IA,JA,N_HYD)
    real(RP) :: QNUM(KA,IA,JA,N_HYD)

    logical :: convert_qtrc
    !---------------------------------------------------------------------------

    if ( MKINIT_TYPE == I_IGNORE ) then
      LOG_NEWLINE
      LOG_PROGRESS(*) 'skip  making initial data'
      output = .false.
    else
      LOG_NEWLINE
      LOG_PROGRESS(*) 'start making initial data'

      !--- Initialize variables
      pres(:,:,:) = CONST_UNDEF8
      temp(:,:,:) = CONST_UNDEF8
      pott(:,:,:) = CONST_UNDEF8
      qsat(:,:,:) = CONST_UNDEF8
      velx(:,:,:) = CONST_UNDEF8
      vely(:,:,:) = CONST_UNDEF8

      rndm  (:,:,:) = CONST_UNDEF8

      pres_sfc(:,:) = CONST_UNDEF8
      temp_sfc(:,:) = CONST_UNDEF8
      pott_sfc(:,:) = CONST_UNDEF8
      psat_sfc(:,:) = CONST_UNDEF8
      qsat_sfc(:,:) = CONST_UNDEF8

      qv    (:,:,:) = 0.0_RP
      qc    (:,:,:) = 0.0_RP
      nc    (:,:,:) = 0.0_RP
      qv_sfc(:,:) = 0.0_RP
      qc_sfc(:,:) = 0.0_RP

!OCL XFILL
      QTRC(:,:,:,:) = 0.0_RP
!OCL XFILL
      QHYD(:,:,:,:) = 0.0_RP
!OCL XFILL
      QNUM(:,:,:,:) = 0.0_RP

      call PROF_rapstart('_MkInit_main',3)

      convert_qtrc = .true.

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
      case(I_BOMEX)
         call MKINIT_BOMEX
      case(I_POINTE)
         call MKINIT_POINTE
      case(I_OCEANCOUPLE)
         call MKINIT_planestate
         call MKINIT_oceancouple
      case(I_LANDCOUPLE)
         call MKINIT_planestate
         call MKINIT_landcouple
      case(I_URBANCOUPLE)
         call MKINIT_planestate
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
         convert_qtrc = .false.
      case(I_GRAYZONE)
         call MKINIT_grayzone
      case(I_BOXAERO)
         call MKINIT_boxaero
      case(I_WARMBUBBLEAERO)
         call MKINIT_warmbubbleaero
      case(I_CAVITYFLOW)
         call MKINIT_cavityflow
      case(I_BAROCWAVE)
         call MKINIT_barocwave
      case default
         LOG_ERROR("MKINIT",*) 'Unsupported TYPE:', MKINIT_TYPE
         call PRC_abort
      endselect

      call tke_setup

      call AEROSOL_setup

      call SBMAERO_setup( convert_qtrc ) ! [INOUT]

      if ( QA_MP > 0 .AND. convert_qtrc ) then
!OCL XFILL
         QHYD(:,:,:,I_HC) = qc(:,:,:)
!OCL XFILL
         QNUM(:,:,:,I_HC) = nc(:,:,:)
         call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                             qv(:,:,:), QHYD(:,:,:,:), & ! [IN]
                                             QTRC(:,:,:,QS_MP:QE_MP),  & ! [OUT]
                                             QNUM=QNUM(:,:,:,:)        ) ! [IN]
      end if

      call PROF_rapend  ('_MkInit_main',3)

      LOG_PROGRESS(*) 'end   making initial data'

      output = .true.

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
    real(RP) :: BBL_RZ       =  0.0_RP  ! bubble radius   [m]: z
    real(RP) :: BBL_RX       =  0.0_RP  ! bubble radius   [m]: x
    real(RP) :: BBL_RY       =  0.0_RP  ! bubble radius   [m]: y

    namelist / PARAM_BUBBLE / &
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
    real(RP) :: distx, disty, distz

    real(RP) :: Domain_RX, Domain_RY

    integer  :: ierr
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("BUBBLE_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BUBBLE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("BUBBLE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("BUBBLE_setup",*) 'Not appropriate names in namelist PARAM_BUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_BUBBLE)

    if ( abs(BBL_RZ*BBL_RX*BBL_RY) <= 0.0_RP ) then
       LOG_INFO("BUBBLE_setup",*) 'no bubble'
       bubble(:,:,:) = 0.0_RP
    else

       bubble(:,:,:) = CONST_UNDEF8

       if ( BBL_eachnode ) then
          CZ_offset = CZ(KS)
          CX_offset = CX(IS)
          CY_offset = CY(JS)
          Domain_RX = FX(IE) - FX(IS-1)
          Domain_RY = FY(JE) - FY(JS-1)
       else
          CZ_offset = 0.0_RP
          CX_offset = 0.0_RP
          CY_offset = 0.0_RP
          Domain_RX = FXG(IAG-IHALO) - FXG(IHALO)
          Domain_RY = FYG(JAG-JHALO) - FYG(JHALO)
       endif

       ! make bubble coefficient
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE

          distz = ( (CZ(k)-CZ_offset-BBL_CZ)/BBL_RZ )**2

          distx = min( ( (CX(i)-CX_offset-BBL_CX          )/BBL_RX )**2, &
                       ( (CX(i)-CX_offset-BBL_CX-Domain_RX)/BBL_RX )**2, &
                       ( (CX(i)-CX_offset-BBL_CX+Domain_RX)/BBL_RX )**2  )

          disty = min( ( (CY(j)-CY_offset-BBL_CY          )/BBL_RY )**2, &
                       ( (CY(j)-CY_offset-BBL_CY-Domain_RY)/BBL_RY )**2, &
                       ( (CY(j)-CY_offset-BBL_CY+Domain_RY)/BBL_RY )**2  )

          bubble(k,i,j) = cos( 0.5_RP*PI*sqrt( min(distz+distx+disty,1.0_RP) ) )**2

       enddo
       enddo
       enddo
    endif

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

    namelist / PARAM_RECT / &
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

    LOG_NEWLINE
    LOG_INFO("RECT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RECT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_ERROR("RECT_setup",*) 'Not found namelist. Check!'
       call PRC_abort
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RECT_setup",*) 'Not appropriate names in namelist PARAM_RECT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_RECT)

    rect(:,:,:) = CONST_UNDEF8

    if ( RCT_eachnode ) then
       CZ_offset = CZ(KS)
       CX_offset = CX(IS)
       CY_offset = CY(JS)
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
            abs(CZ(k) - CZ_offset - RCT_CZ)/RCT_RZ,   &
            abs(CX(i) - CX_offset - RCT_CX)/RCT_RX,   &
            abs(CY(j) - CY_offset - RCT_CY)/RCT_RY    &
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
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_mkinit
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    implicit none

    real(RP), parameter :: d_min_def = 1.e-9_RP ! default lower bound of 1st size bin
    real(RP), parameter :: d_max_def = 1.e-5_RP ! upper bound of last size bin
    integer,  parameter :: n_kap_def = 1        ! number of kappa bins
    real(RP), parameter :: k_min_def = 0.e0_RP  ! lower bound of 1st kappa bin
    real(RP), parameter :: k_max_def = 1.e0_RP  ! upper bound of last kappa bin

    real(RP) :: m0_init = 0.0_RP    ! initial total num. conc. of modes (Atk,Acm,Cor) [#/m3]
    real(RP) :: dg_init = 80.e-9_RP ! initial number equivalen diameters of modes     [m]
    real(RP) :: sg_init = 1.6_RP    ! initial standard deviation                      [-]

    real(RP) :: d_min_inp(3) = d_min_def
    real(RP) :: d_max_inp(3) = d_max_def
    real(RP) :: k_min_inp(3) = k_min_def
    real(RP) :: k_max_inp(3) = k_max_def
    integer  :: n_kap_inp(3) = n_kap_def

    namelist / PARAM_AERO / &
       m0_init,   &
       dg_init,   &
       sg_init,   &
       d_min_inp, &
       d_max_inp, &
       k_min_inp, &
       k_max_inp, &
       n_kap_inp

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_AE_TYPE /= 'KAJINO13' ) return

    LOG_NEWLINE
    LOG_INFO("AEROSOL_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_AERO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("AEROSOL_setup",*) 'Not found namelist. Default used!'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("AEROSOL_setup",*) 'Not appropriate names in namelist PARAM_AERO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_AERO)

    qdry(:,:,:) = 1.0_RP - qv(:,:,:) - qc(:,:,:)
    call ATMOS_PHY_AE_kajino13_mkinit( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
                                       QA_AE,                   & ! (in)
                                       DENS(:,:,:),             & ! (in)
                                       TEMP(:,:,:),             & ! (in)
                                       PRES(:,:,:),             & ! (in)
                                       QDRY(:,:,:),             & ! (in)
                                       QV  (:,:,:),             & ! (in)
                                       m0_init,                 & ! (in)
                                       dg_init,                 & ! (in)
                                       sg_init,                 & ! (in)
                                       d_min_inp(:),            & ! (in)
                                       d_max_inp(:),            & ! (in)
                                       k_min_inp(:),            & ! (in)
                                       k_max_inp(:),            & ! (in)
                                       n_kap_inp(:),            & ! (in)
                                       QTRC(:,:,:,QS_AE:QE_AE), & ! (out)
                                       CCN(:,:,:)               ) ! (out)

    return
  end subroutine AEROSOL_setup

  !-----------------------------------------------------------------------------
  !> Setup aerosol condition for Spectral Bin Microphysics (SBM) model
  subroutine SBMAERO_setup( convert_qtrc )
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       QHS,  &
       QHE
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use scale_atmos_phy_mp_suzuki10, only: &
       nccn, nbin
    implicit none

    logical, intent(inout) :: convert_qtrc

    real(RP), allocatable :: xabnd(:), xactr(:)

    integer :: iq, i, j, k
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_MP_TYPE /= 'SUZUKI10' ) return

    !--- Super saturated air at initial
    do j = JSB, JEB
    do i = ISB, IEB
    do k  = KS, KE
       QTRC(k,i,j,I_QV) = qv(k,i,j) + qc(k,i,j)
    end do
    end do
    end do

    !-- Aerosol distribution
    if ( nccn /= 0 ) then
       do iq = 1, nccn
       do j = JSB, JEB
       do i = ISB, IEB
       do k  = KS, KE
          QTRC(k,i,j,QHE+iq) = gan(iq) / DENS(k,i,j) ! [note] gan is never set.
       enddo
       enddo
       enddo
       enddo

       deallocate( xactr )
       deallocate( xabnd )
    endif

    convert_qtrc = .false.

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
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down
    implicit none

    ! Flux from Atmosphere
    real(RP) :: FLX_rain   = 0.0_RP ! surface rain flux              [kg/m2/s]
    real(RP) :: FLX_snow   = 0.0_RP ! surface snow flux              [kg/m2/s]
    real(RP) :: FLX_IR_dn  = 0.0_RP ! surface downwad radiation flux [J/m2/s]
    real(RP) :: FLX_NIR_dn = 0.0_RP ! surface downwad radiation flux [J/m2/s]
    real(RP) :: FLX_VIS_dn = 0.0_RP ! surface downwad radiation flux [J/m2/s]

    namelist / PARAM_MKINIT_FLUX / &
       FLX_rain,   &
       FLX_snow,   &
       FLX_IR_dn,  &
       FLX_NIR_dn, &
       FLX_VIS_dn

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_FLUX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("flux_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("flux_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_FLUX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_FLUX)

    do j = JSB, JEB
    do i = ISB, IEB
       SFLX_rain   (i,j) = FLX_rain
       SFLX_snow   (i,j) = FLX_snow

       SFLX_LW_up  (i,j) = 0.0_RP
       SFLX_LW_dn  (i,j) = FLX_IR_dn
       SFLX_SW_up  (i,j) = 0.0_RP
       SFLX_SW_dn  (i,j) = FLX_NIR_dn + FLX_VIS_dn

       TOAFLX_LW_up(i,j) = 0.0_RP
       TOAFLX_LW_dn(i,j) = 0.0_RP
       TOAFLX_SW_up(i,j) = 0.0_RP
       TOAFLX_SW_dn(i,j) = 0.0_RP

       SFLX_rad_dn (i,j,I_R_direct ,I_R_IR)  = 0.0_RP
       SFLX_rad_dn (i,j,I_R_diffuse,I_R_IR)  = FLX_IR_dn
       SFLX_rad_dn (i,j,I_R_direct ,I_R_NIR) = FLX_NIR_dn
       SFLX_rad_dn (i,j,I_R_diffuse,I_R_NIR) = 0.0_RP
       SFLX_rad_dn (i,j,I_R_direct ,I_R_VIS) = FLX_VIS_dn
       SFLX_rad_dn (i,j,I_R_diffuse,I_R_VIS) = 0.0_RP
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

    real(RP) :: LND_TEMP                ! land soil temperature      [K]
    real(RP) :: LND_WATER     = 0.15_RP ! land soil moisture         [m3/m3]
    real(RP) :: SFC_TEMP                ! land skin temperature      [K]
    real(RP) :: SFC_albedo_LW = 0.01_RP ! land surface albedo for LW (0-1)
    real(RP) :: SFC_albedo_SW = 0.20_RP ! land surface albedo for SW (0-1)

    namelist / PARAM_MKINIT_LAND / &
       LND_TEMP,      &
       LND_WATER,     &
       SFC_TEMP,      &
       SFC_albedo_LW, &
       SFC_albedo_SW

    integer  :: ierr
    !---------------------------------------------------------------------------

    LND_TEMP = THETAstd
    SFC_TEMP = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("land_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("land_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_LAND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_LAND)

    LAND_TEMP      (:,:,:)         = LND_TEMP
    LAND_WATER     (:,:,:)         = LND_WATER

    LAND_SFC_TEMP  (:,:)           = SFC_TEMP
    LAND_SFC_albedo(:,:,:,I_R_IR)  = SFC_albedo_LW
    LAND_SFC_albedo(:,:,:,I_R_NIR) = SFC_albedo_SW
    LAND_SFC_albedo(:,:,:,I_R_VIS) = SFC_albedo_SW

    return
  end subroutine land_setup

  !-----------------------------------------------------------------------------
  !> Ocean setup
  subroutine ocean_setup
    use mod_ocean_vars, only: &
       OCEAN_TEMP,       &
       OCEAN_SALT,       &
       OCEAN_UVEL,       &
       OCEAN_VVEL,       &
       OCEAN_OCN_Z0M,    &
       OCEAN_ICE_TEMP,   &
       OCEAN_ICE_MASS,   &
       OCEAN_SFC_TEMP,   &
       OCEAN_SFC_albedo, &
       OCEAN_SFC_Z0M,    &
       OCEAN_SFC_Z0H,    &
       OCEAN_SFC_Z0E
    implicit none

    real(RP) :: OCN_TEMP                 ! ocean temperature                         [K]
    real(RP) :: OCN_SALT      = 0.0_RP   ! ocean salinity                            [psu]
    real(RP) :: OCN_UVEL      = 0.0_RP   ! ocean u-velocity                          [m/s]
    real(RP) :: OCN_VVEL      = 0.0_RP   ! ocean v-velocity                          [m/s]
    real(RP) :: ICE_TEMP                 ! ocean temperature                         [K]
    real(RP) :: ICE_MASS      = 0.0_RP   ! ocean temperature                         [K]
    real(RP) :: SFC_TEMP                 ! ocean skin temperature                    [K]
    real(RP) :: SFC_albedo_LW = 0.04_RP  ! ocean surface albedo for LW               (0-1)
    real(RP) :: SFC_albedo_SW = 0.05_RP  ! ocean surface albedo for SW               (0-1)
    real(RP) :: SFC_Z0M       = 1.E-4_RP ! ocean surface roughness length (momentum) [m]
    real(RP) :: SFC_Z0H       = 1.E-4_RP ! ocean surface roughness length (heat)     [m]
    real(RP) :: SFC_Z0E       = 1.E-4_RP ! ocean surface roughness length (vapor)    [m]

    namelist / PARAM_MKINIT_OCEAN / &
       OCN_TEMP,      &
       OCN_SALT,      &
       OCN_UVEL,      &
       OCN_VVEL,      &
       ICE_TEMP,      &
       ICE_MASS,      &
       SFC_TEMP,      &
       SFC_albedo_LW, &
       SFC_albedo_SW, &
       SFC_Z0M,       &
       SFC_Z0H,       &
       SFC_Z0E

    integer :: ierr
    !---------------------------------------------------------------------------

    OCN_TEMP = THETAstd
    ICE_TEMP = 271.35_RP ! freezing point of the ocean
    SFC_TEMP = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ocean_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ocean_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_OCEAN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_OCEAN)

    OCEAN_TEMP      (:,:,:) = OCN_TEMP
    OCEAN_SALT      (:,:,:) = OCN_SALT
    OCEAN_UVEL      (:,:,:) = OCN_UVEL
    OCEAN_VVEL      (:,:,:) = OCN_VVEL
    OCEAN_OCN_Z0M   (:,:)   = SFC_Z0M
    OCEAN_ICE_TEMP  (:,:)   = ICE_TEMP
    OCEAN_ICE_MASS  (:,:)   = ICE_MASS

    OCEAN_SFC_TEMP  (:,:)           = SFC_TEMP
    OCEAN_SFC_albedo(:,:,:,I_R_IR)  = SFC_albedo_LW
    OCEAN_SFC_albedo(:,:,:,I_R_NIR) = SFC_albedo_SW
    OCEAN_SFC_albedo(:,:,:,I_R_VIS) = SFC_albedo_SW
    OCEAN_SFC_Z0M   (:,:)           = SFC_Z0M
    OCEAN_SFC_Z0H   (:,:)           = SFC_Z0H
    OCEAN_SFC_Z0E   (:,:)           = SFC_Z0E

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

    real(RP) :: URB_ROOF_TEMP                 ! Surface temperature of roof           [K]
    real(RP) :: URB_BLDG_TEMP                 ! Surface temperature of building       [K]
    real(RP) :: URB_GRND_TEMP                 ! Surface temperature of ground         [K]
    real(RP) :: URB_CNPY_TEMP                 ! Diagnostic canopy air temperature     [K]
    real(RP) :: URB_CNPY_HMDT       = 0.0_RP  ! Diagnostic canopy humidity            [kg/kg]
    real(RP) :: URB_CNPY_WIND       = 0.0_RP  ! Diagnostic canopy wind                [m/s]
    real(RP) :: URB_ROOF_LAYER_TEMP           ! temperature in layer of roof          [K]
    real(RP) :: URB_BLDG_LAYER_TEMP           ! temperature in layer of building      [K]
    real(RP) :: URB_GRND_LAYER_TEMP           ! temperature in layer of ground        [K]
    real(RP) :: URB_ROOF_RAIN       = 0.0_RP  ! temperature in layer of roof          [K]
    real(RP) :: URB_BLDG_RAIN       = 0.0_RP  ! temperature in layer of building      [K]
    real(RP) :: URB_GRND_RAIN       = 0.0_RP  ! temperature in layer of ground        [K]
    real(RP) :: URB_RUNOFF          = 0.0_RP  ! temperature in layer of ground        [K]
    real(RP) :: URB_SFC_TEMP                  ! Grid average of surface temperature   [K]
    real(RP) :: URB_ALB_LW          = 0.10_RP ! Grid average of surface albedo for LW (0-1)
    real(RP) :: URB_ALB_SW          = 0.20_RP ! Grid average of surface albedo for SW (0-1)

    namelist / PARAM_MKINIT_URBAN / &
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

    integer :: ierr
    !---------------------------------------------------------------------------

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
       LOG_INFO("urban_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("urban_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_URBAN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_URBAN)

    URBAN_TRL       (:,:,:)         = URB_ROOF_LAYER_TEMP
    URBAN_TBL       (:,:,:)         = URB_BLDG_LAYER_TEMP
    URBAN_TGL       (:,:,:)         = URB_GRND_LAYER_TEMP

    URBAN_TR        (:,:)           = URB_ROOF_TEMP
    URBAN_TB        (:,:)           = URB_BLDG_TEMP
    URBAN_TG        (:,:)           = URB_GRND_TEMP
    URBAN_TC        (:,:)           = URB_CNPY_TEMP
    URBAN_QC        (:,:)           = URB_CNPY_HMDT
    URBAN_UC        (:,:)           = URB_CNPY_WIND
    URBAN_RAINR     (:,:)           = URB_ROOF_RAIN
    URBAN_RAINB     (:,:)           = URB_BLDG_RAIN
    URBAN_RAING     (:,:)           = URB_GRND_RAIN
    URBAN_ROFF      (:,:)           = URB_RUNOFF
    URBAN_SFC_TEMP  (:,:)           = URB_SFC_TEMP
    URBAN_SFC_albedo(:,:,:,I_R_IR)  = URB_ALB_LW
    URBAN_SFC_albedo(:,:,:,I_R_NIR) = URB_ALB_SW
    URBAN_SFC_albedo(:,:,:,I_R_VIS) = URB_ALB_SW

    return
  end subroutine urban_setup

  !-----------------------------------------------------------------------------
  !> TKE setup
  subroutine tke_setup
    use scale_const, only: &
       EPS => CONST_EPS
    use mod_atmos_phy_tb_vars, only: &
       I_TKE
    use mod_atmos_phy_bl_vars, only: &
       QS_BL => QS, &
       QE_BL => QE
    implicit none

    real(RP) :: TKE_CONST

    namelist / PARAM_MKINIT_TKE / &
       TKE_CONST

    integer :: k, i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    TKE_CONST = EPS

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TKE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("tke_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("tke_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_TKE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_TKE)

    if ( I_TKE > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          QTRC(k,i,j,I_TKE) = TKE_CONST
       enddo
       enddo
       enddo
    end if
    if ( QS_BL > 0 ) then
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          QTRC(k,i,j,QS_BL) = TKE_CONST
          QTRC(k,i,j,QS_BL+1:QE_BL) = 0.0_RP
       enddo
       enddo
       enddo
    end if

    return
  end subroutine tke_setup

  !-----------------------------------------------------------------------------
  !> Read sounding data from file
  subroutine read_sounding( &
       DENS, VELX, VELY, POTT, QV )
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_SOUNDING / &
       ENV_IN_SOUNDING_file

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SOUNDING,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("read_sounding",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("read_sounding",*) 'Not appropriate names in namelist PARAM_MKINIT_SOUNDING. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SOUNDING)

    !--- prepare sounding profile
    LOG_INFO("read_sounding",*) 'Input sounding file:', trim(ENV_IN_SOUNDING_file)
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

       if ( ierr /= 0 ) then
          LOG_ERROR("read_sounding",*) '[mod_mkinit/read_sounding] Input file not found!'
       endif

       !--- read sounding file till end
       read(fid,*) SFC_PRES, SFC_THETA, SFC_QV

       LOG_INFO("read_sounding",*) '+ Surface pressure [hPa]',     SFC_PRES
       LOG_INFO("read_sounding",*) '+ Surface pot. temp  [K]',     SFC_THETA
       LOG_INFO("read_sounding",*) '+ Surface water vapor [g/kg]', SFC_QV

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
    pres_sfc(:,:) = SFC_PRES * 1.E2_RP ! [hPa]->[Pa]
    pott_sfc(:,:) = SFC_THETA
    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       qv_sfc(:,:)   = SFC_QV * 1.E-3_RP ! [g/kg]->[kg/kg]
    end if

    !--- linear interpolate to model grid
    do k = KS, KE
       do kref = 2, EXP_kmax+1
          if (       CZ(k) >  EXP_z(kref-1) &
               .AND. CZ(k) <= EXP_z(kref  ) ) then

             fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

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
    if ( ATMOS_HYDROMETEOR_dry ) qv(:) = 0.0_RP

    qc(:) = 0.0_RP

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:), qv(:), qc(:),                                  & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:), temp(:), pres(:), temp_sfc(1,1)                ) ! [OUT]

    return
  end subroutine read_sounding

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform + random disturbance )
  subroutine MKINIT_planestate
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_PLANESTATE / &
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

    LOG_NEWLINE
    LOG_INFO("MKINIT_planestate",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_PLANESTATE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_planestate",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR_CONT(*) 'Not appropriate names in namelist PARAM_MKINIT_PLANESTATE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_PLANESTATE)

    ! calc in dry condition
    do j = JSB, JEB
    do i = ISB, IEB
       pott_sfc(i,j) = SFC_THETA
       pres_sfc(i,j) = SFC_PRES
    enddo
    enddo

    if ( ENV_THETA < 0.0_RP ) then ! use isa profile

       call PROFILE_isa( KA, KS, KE,      & ! [IN]
                         IA, ISB, IEB,    & ! [IN]
                         JA, JSB, JEB,    & ! [IN]
                         pott_sfc(:,:), & ! [IN]
                         pres_sfc(:,:), & ! [IN]
                         REAL_CZ (:,:,:), & ! [IN]
                         pott    (:,:,:)  ) ! [OUT]

    else

       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * REAL_CZ(k,i,j)
       enddo
       enddo
       enddo

    endif

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       ! calc QV from RH
       call SATURATION_psat_all( IA, ISB, IEB, JA, JSB, JEB, &
                                 temp_sfc(:,:), & ! [IN]
                                 psat_sfc(:,:)  ) ! [OUT]
       qdry(:,:,:) = 1.0_RP - qv(:,:,:) - qc(:,:,:)
       call SATURATION_pres2qsat_all( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                      temp(:,:,:), pres(:,:,:), qdry(:,:,:), & ! [IN]
                                      qsat(:,:,:)                            ) ! [OUT]

       call RANDOM_uniform(rndm) ! make random
       do j = JSB, JEB
       do i = ISB, IEB
          qsat_sfc(i,j) = EPSvap * psat_sfc(i,j) / ( pres_sfc(i,j) - ( 1.0_RP-EPSvap ) * psat_sfc(i,j) )
          qv_sfc(i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(i,j)

          do k = KS, KE
             qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat(k,i,j)
          enddo
       enddo
       enddo
    end if

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
       pott_sfc(i,j) = pott_sfc(i,j) + rndm(KS-1,i,j) * RANDOM_THETA

       do k = KS, KE
          pott(k,i,j) = pott(k,i,j) + rndm(k,i,j) * RANDOM_THETA
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, min(IEB,IA-1)
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_U ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, min(JEB,JA-1)
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_V ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call flux_setup

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

    namelist / PARAM_MKINIT_TRACERBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       ENV_U,     &
       ENV_V,     &
       SHAPE_NC,  &
       BBL_NC

    real(RP), pointer :: shapeFac(:,:,:) => null()

    integer :: k, i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_tracerbubble",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TRACERBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_tracerbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_tracerbubble",*) 'Not appropriate names in namelist PARAM_MKINIT_TRACERBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_TRACERBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       pott(k,1,1) = ENV_THETA
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U       * DENS(k,1,1)
       MOMY(k,i,j) = ENV_V       * DENS(k,1,1)
       RHOT(k,i,j) = pott(k,1,1) * DENS(k,1,1)
    enddo
    enddo
    enddo

    ! make tracer bubble
    select case(SHAPE_NC)
    case('BUBBLE')
       call BUBBLE_setup
       shapeFac => bubble
    case('RECT')
       call RECT_setup
       shapeFac => rect
    case default
       LOG_ERROR("MKINIT_tracerbubble",*) 'SHAPE_NC=', trim(SHAPE_NC), ' cannot be used on advect. Check!'
       call PRC_abort
    end select

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       nc(k,i,j) = BBL_NC * shapeFac(k,i,j)
    enddo
    enddo
    enddo

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

    namelist / PARAM_MKINIT_COLDBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       BBL_TEMP

    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_coldbubble",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_COLDBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_coldbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_coldbubble",*) 'Not appropriate names in namelist PARAM_MKINIT_COLDBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_COLDBUBBLE)

    RovCP = Rdry / CPdry

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       pott(k,1,1) = ENV_THETA
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

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
    enddo
    enddo
    enddo

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

    namelist / PARAM_MKINIT_LAMBWAVE / &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       ENV_TEMP,  &
       BBL_PRES

    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_lambwave",*) 'Setup initial state'

    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_LAMBWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_lambwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_lambwave",*) 'Not appropriate names in namelist PARAM_MKINIT_LAMBWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_LAMBWAVE)

    RovCP = Rdry / CPdry

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = SFC_PRES/(Rdry*ENV_TEMP) * exp( - GRAV/(Rdry*ENV_TEMP) * CZ(k) )
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make pressure bubble
       pres(k,i,j) = DENS(k,i,j) * ENV_TEMP * Rdry + BBL_PRES * bubble(k,i,j)

       RHOT(k,i,j) = DENS(k,i,j) * ENV_TEMP * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

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

    namelist / PARAM_MKINIT_GRAVITYWAVE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       ENV_BVF,   &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_gravitywave",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_GRAVITYWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_gravitywave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_gravitywave",*) 'Not appropriate names in namelist PARAM_MKINIT_GRAVITYWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_GRAVITYWAVE)

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       pott(k,1,1) = SFC_THETA * exp( ENV_BVF*ENV_BVF / GRAV * CZ(k) )
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,1,1)
       MOMY(k,i,j) = ENV_V * DENS(k,1,1)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

    enddo
    enddo
    enddo

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

    namelist / PARAM_MKINIT_KHWAVE / &
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
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_khwave",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_KHWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_khwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_khwave",*) 'Not appropriate names in namelist PARAM_MKINIT_KHWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_KHWAVE)

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )
       fact = max( min( fact, 1.0_RP ), 0.0_RP )

       pott(k,1,1) = ENV_L1_THETA * ( 1.0_RP - fact ) &
                   + ENV_L3_THETA * (          fact )
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       RHOT(k,i,j) = DENS(k,1,1) * pott(k,1,1)
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )
       fact = max( min( fact, 1.0_RP ), 0.0_RP )

       MOMX(k,i,j) = ( ENV_L1_U * ( 1.0_RP - fact )                 &
                     + ENV_L3_U * (          fact )                 &
                     + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_U &
                     ) * DENS(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_khwave

  !-----------------------------------------------------------------------------
  !> Make initial state for turbulence experiment
  subroutine MKINIT_turbulence
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_TURBULENCE / &
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

    LOG_NEWLINE
    LOG_INFO("MKINIT_turbulence",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TURBULENCE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_turbulence",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_turbulence",*) 'Not appropriate names in namelist PARAM_MKINIT_TURBULENCE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_TURBULENCE)

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       pott(k,1,1) = ENV_THETA + ENV_TLAPS * CZ(k)
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       ! calc QV from RH
       call SATURATION_psat_all( temp_sfc(1,1), & ! [IN]
                                 psat_sfc(1,1)  ) ! [OUT]
       qdry(:,1,1) = 1.0_RP - qv(:,1,1) - qc(:,1,1)
       call SATURATION_pres2qsat_all( KA, KS, KE, &
                                      temp(:,1,1), pres(:,1,1), qdry(:,1,1), & ! [IN]
                                      qsat(:,1,1)                            ) ! [OUT]

       call RANDOM_uniform(rndm) ! make random
       do j = JSB, JEB
       do i = ISB, IEB
          qsat_sfc(1,1) = EPSvap * psat_sfc(i,j) / ( pres_sfc(i,j) - ( 1.0_RP-EPSvap ) * psat_sfc(i,j) )
          qv_sfc(i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(1,1)

          do k = KS, KE
             qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat(k,1,1)
          enddo
       enddo
       enddo
    end if

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
       pres_sfc(i,j) = SFC_PRES
       pott_sfc(i,j) = SFC_THETA + rndm(KS-1,i,j) * RANDOM_THETA

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * CZ(k) + rndm(k,i,j) * RANDOM_THETA
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_U ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_V ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_turbulence

  !-----------------------------------------------------------------------------
  !> Make initial state for cavity flow
  subroutine MKINIT_cavityflow
    implicit none

    ! Nondimenstional numbers for a cavity flow problem
    real(RP) :: REYNOLDS_NUM = 1.D03
    real(RP) :: MACH_NUM     = 3.D-2
    real(RP) :: Ulid         = 1.D01
    real(RP) :: PRES0        = 1.D05

    namelist / PARAM_MKINIT_CAVITYFLOW / &
         Ulid        , &
         PRES0       , &
         REYNOLDS_NUM, &
         MACH_NUM

    real(RP) :: DENS0
    real(RP) :: TEMP
    real(RP) :: Gam
    real(RP) :: Cs2

    integer :: k, i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_cavityflow",*) 'Setup initial state'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_CAVITYFLOW,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_cavityflow",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_cavityflow",*) 'Not appropriate names in namelist PARAM_MKINIT_CAVITYFLOW. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_CAVITYFLOW)

    Gam   = CPdry / ( CPdry - Rdry )
    Cs2   = ( Ulid / MACH_NUM )**2
    TEMP  = Cs2   / ( Gam * Rdry )
    DENS0 = PRES0 / ( Rdry * TEMP )

    LOG_INFO("MKINIT_cavityflow",*) "DENS = ", DENS0
    LOG_INFO("MKINIT_cavityflow",*) "PRES = ", PRES0
    LOG_INFO("MKINIT_cavityflow",*) "TEMP = ", RHOT(10,10,4)/DENS0, TEMP
    LOG_INFO("MKINIT_cavityflow",*) "Ulid = ", Ulid
    LOG_INFO("MKINIT_cavityflow",*) "Cs   = ", sqrt(Cs2)

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j) = DENS0
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       PRES(k,i,j) = PRES0
       RHOT(k,i,j) = P00/Rdry * (P00/PRES0)**((Rdry - CPdry)/CPdry)
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

    namelist / PARAM_MKINIT_MOUNTAINWAVE / &
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

    LOG_NEWLINE
    LOG_INFO("MKINIT_mountainwave",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_MOUNTAINWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_mountainwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_mountainwave",*) 'Not appropriate names in namelist PARAM_MKINIT_MOUNTAINWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_MOUNTAINWAVE)

    ! calc in dry condition
    do j = JSB, JEB
    do i = ISB, IEB
       pres_sfc(i,j) = SFC_PRES
       pott_sfc(i,j) = SFC_THETA
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       Ustar2 = ENV_U * ENV_U + ENV_V * ENV_V
       N2     = Ustar2 * (SCORER*SCORER)

       pott(k,i,j) = SFC_THETA * exp( N2 / GRAV * REAL_CZ(k,i,j) )
    enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,i,j)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U       * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V       * DENS(k,i,j)
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    ! optional : add tracer bubble
    if ( BBL_NC > 0.0_RP ) then
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          nc(k,i,j) = BBL_NC * bubble(k,i,j)
       enddo
       enddo
       enddo
    endif

    return
  end subroutine MKINIT_mountainwave

  !-----------------------------------------------------------------------------
  !> Make initial state
  !!
  !! The initial state for a baroclinic wave test described by Ullrich and Jablonowski(2012)
  !! is generated.
  subroutine MKINIT_barocwave
    use scale_const, only: &
      OHM => CONST_OHM,        &
      RPlanet => CONST_RADIUS, &
      GRAV    => CONST_GRAV
    use scale_prc
    use scale_atmos_grid_cartesC, only: &
         y0  => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
         FYG => ATMOS_GRID_CARTESC_FYG
    use scale_atmos_hydrometeor, only: &
         ATMOS_HYDROMETEOR_dry

    implicit none

    ! Parameters for global domain size
    real(RP) :: Ly                      ! The domain size in y-direction [m]

    ! Parameters for inital stratification
    real(RP) :: REF_TEMP    = 288.E0_RP ! The reference temperature [K]
    real(RP) :: REF_PRES    = 1.E5_RP   ! The reference pressure [Pa]
    real(RP) :: LAPSE_RATE  = 5.E-3_RP  ! The lapse rate [K/m]

    ! Parameters associated with coriolis parameter on a beta-plane
    real(RP) :: Phi0Deg     = 45.E0_RP  ! The central latitude [degree_north]

    ! Parameters for background zonal jet
    real(RP) :: U0 = 35.E0_RP          ! The parameter associated with zonal jet maximum amplitude  [m/s]
    real(RP) :: b  = 2.E0_RP           ! The vertical half-width [1]

    ! Parameters for inital perturbation of zonal wind with a Gaussian profile
    !
    real(RP) :: Up  = 1.E0_RP         ! The maximum amplitude of zonal wind perturbation [m/s]
    real(RP) :: Lp  = 600.E3_RP       ! The width of Gaussian profile
    real(RP) :: Xc  = 2000.E3_RP      ! The center point (x) of inital perturbation
    real(RP) :: Yc  = 2500.E3_RP      ! The center point (y) of inital perturbation

    namelist / PARAM_MKINIT_BAROCWAVE / &
       REF_TEMP, REF_PRES, LAPSE_RATE, &
       phi0Deg,                        &
       U0, b,                          &
       Up, Lp, Xc, Yc

    real(RP) :: f0, beta0

    real(RP) :: geopot(KA,IA,JA)
    real(RP) :: eta(KA,IA,JA)
    real(RP) :: temp(KA,IA,JA)

    real(RP) :: y
    real(RP) :: ln_eta
    real(RP) :: del_eta
    real(RP) :: yphase
    real(RP) :: temp_vfunc
    real(RP) :: geopot_hvari

    integer :: ierr
    integer :: k, i, j

    integer :: itr

    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_barocwave",*) 'Setup initial state'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BAROCWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_barocwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_barocwave",*) 'Not appropriate names in namelist PARAM_MKINIT_BAROCWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_BAROCWAVE)

    Ly = FYG(JAG-JHALO) - FYG(JHALO)

    ! Set coriolis parameters
    f0 = 2.0_RP*OHM*sin(phi0Deg*PI/180.0_RP)
    beta0 = (2.0_RP*OHM/RPlanet)*cos(phi0Deg*PI/180.0_RP)

    ! Calculate eta(=p/p_s) level corresponding to z level of each (y,z) grid point
    ! using Newton's iteration method

    eta(:,:,:) = 1.0E-8_RP   ! Set first guess of eta

    do j = JSB, JEB
    do i = ISB, IEB            ! Note that initial fields are zonaly symmetric

       y = CY(j)
       yphase  = 2.0_RP*PI*y/Ly

       ! Calc horizontal variation of geopotential height
       geopot_hvari = 0.5_RP*U0*(                                                                          &
            (f0 - beta0*y0)*(y - 0.5_RP*Ly*(1.0_RP + sin(yphase)/PI))                                      &
            + 0.5_RP*beta0*( y**2 - Ly*y/PI*sin(yphase) - 0.5_RP*(Ly/PI)**2*(cos(yphase) + 1.0_RP)         &
                             - Ly**2/3.0_RP                                                        )       &
            )

       ! Set surface pressure and temperature
       pres_sfc(i,j) = REF_PRES
       pott_sfc(i,j) = REF_TEMP - geopot_hvari/Rdry

       do k = KS, KE
          del_eta = 1.0_RP

          !-- The loop for iteration
          itr = 0
          do while( abs(del_eta) > CONV_EPS )
             ln_eta = log(eta(k,i,j))

             temp_vfunc = eta(k,i,j)**(Rdry*LAPSE_RATE/Grav)
             temp(k,i,j) = &
                  REF_TEMP*temp_vfunc &
                  + geopot_hvari/Rdry*(2.0_RP*(ln_eta/b)**2 - 1.0_RP)*exp(-(ln_eta/b)**2)
             geopot(k,i,j) = &
                  REF_TEMP*GRAV/LAPSE_RATE*(1.0_RP - temp_vfunc)  &
                  + geopot_hvari*ln_eta*exp(-(ln_eta/b)**2)

             del_eta = -  ( - Grav*CZ(k) + geopot(k,i,j) )    & ! <- F
                  &      *( - eta(k,i,j)/(Rdry*temp(k,i,j))   ) ! <- (dF/deta)^-1

             eta(k,i,j) = eta(k,i,j) + del_eta
             itr = itr + 1

             if ( itr > ITRMAX ) then
                LOG_ERROR("MKINIT_barocwave",*) "Fail the convergence of iteration. Check!"
                LOG_ERROR_CONT(*) "* (X,Y,Z)=", CX(i), CY(j), CZ(k)
                LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta(k,i,j), "temp=", temp(k,i,j)
                call PRC_abort
             end if
          enddo !- End of loop for iteration ----------------------------

          PRES(k,i,j) = eta(k,i,j)*REF_PRES
          DENS(k,i,j) = PRES(k,i,j)/(Rdry*temp(k,i,j))
          pott(k,i,j) = temp(k,i,j)*eta(k,i,j)**(-Rdry/CPdry)

       enddo

       ! Make density & pressure profile in dry condition using the profile of
       ! potential temperature calculated above.
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  pott(:,i,j), qv(:,i,j), qc(:,i,j),                      & ! [IN]
                                  pres_sfc(i,j), pott_sfc(i,j), qv_sfc(i,j), qc_sfc(i,j), & ! [IN]
                                  REAL_CZ(:,i,j), REAL_FZ(:,i,j),                         & ! [IN]
                                  DENS(:,i,j), temp(:,i,j), pres(:,i,j), temp_sfc(i,j)    ) ! [OUT]
    enddo
    enddo

    !-----------------------------------------------------------------------------------

    do j = JSB, JEB
    do k = KS, KE

       eta(k,IS,j) = pres(k,IS,j)/REF_PRES
       ln_eta = log(eta(k,IS,j))
       yphase = 2.0_RP*PI*CY(j)/Ly
!!$       PRES(k,IS:IE,j) = eta(k,IS,j)*REF_PRES
!!$       DENS(k,IS:IE,j) = PRES(k,IS,j)/(Rdry*temp(k,IS,j))
       DENS(k,IS:IE,j) = DENS(k,IS,j)
       PRES(k,IS:IE,j) = PRES(k,IS,j)
       MOMX(k,IS-1:IE,j) = DENS(k,IS,j)*(-U0*sin(0.5_RP*yphase)**2*ln_eta*exp(-(ln_eta/b)**2))
       RHOT(k,IS:IE,j) = DENS(k,IS,j)*pott(k,IS,j) !temp(k,IS,j)*eta(k,IS,j)**(-Rdry/CPdry)
    enddo
    enddo
    MOMY(:,:,:) = 0.0_RP
    MOMZ(:,:,:) = 0.0_RP

    !---------------------------------------------------------------------------------------

    ! Add the inital perturbation for zonal velocity
    do j = JSB, JEB
    do i = max(ISB-1,1), IEB
       MOMX(KS:kE,i,j) = MOMX(KS:KE,i,j) &
           +  DENS(KS:KE,i,j)* Up*exp( - ((FX(i) - Xc)**2 + (CY(j) - Yc)**2)/Lp**2 )
    enddo
    enddo

    return
  end subroutine MKINIT_barocwave

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  subroutine MKINIT_warmbubble
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_WARMBUBBLE / &
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

    LOG_NEWLINE
    LOG_INFO("MKINIT_warmbubble",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_warmbubble",*) 'QV is not registered'
       call PRC_abort
    end if

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_warmbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_warmbubble",*) 'Not appropriate names in namelist PARAM_MKINIT_WARMBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( CZ(k)-CZ(k-1) )
       else                                ! Layer 3
          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( CZ(k)-CZ(k-1) )
       endif
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    ! calc QV from RH
    call SATURATION_psat_all( temp_sfc(1,1), & ! [IN]
                              psat_sfc(1,1)  ) ! [OUT]
    qsat_sfc(1,1) = EPSvap * psat_sfc(1,1) / ( pres_sfc(1,1) - ( 1.0_RP-EPSvap ) * psat_sfc(1,1) )
    qv_sfc(1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1)
    qdry(:,1,1) = 1.0_RP - qv(:,1,1) - qc(:,1,1)
    call SATURATION_pres2qsat_all( KA, KS, KE, &
                                   temp(:,1,1), pres(:,1,1), qdry(:,1,1), & ! [IN]
                                   qsat(:,1,1)                            ) ! [OUT]
    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       elseif( CZ(k) <= ENV_L2_ZTOP ) then ! Layer 2
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       else                                ! Layer 3
          qv(k,1,1) = 0.0_RP
       endif
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       qv  (k,i,j) = qv(k,1,1)
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_warmbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for supercell experiment
  subroutine MKINIT_supercell
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)

    ! Bubble
    real(RP) :: BBL_THETA = 3.D0 ! extremum of temperature in bubble [K]

    namelist / PARAM_MKINIT_SUPERCELL / &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_supercell",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_supercell",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SUPERCELL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_supercell",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_supercell",*) 'Not appropriate names in namelist PARAM_MKINIT_SUPERCELL. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SUPERCELL)

    call read_sounding( RHO, VELX, VELY, POTT, QV1D ) ! (out)

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = RHO(k) * VELX(k)
       MOMY(k,i,j) = RHO(k) * VELY(k)

       ! make warm bubble
       RHOT(k,i,j) = RHO(k) * ( POTT(k) + BBL_THETA * bubble(k,i,j) )

       qv  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_supercell

  !-----------------------------------------------------------------------------
  !> Make initial state for squallline experiment
  subroutine MKINIT_squallline
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)

    real(RP) :: RANDOM_THETA =  0.01_RP
    real(RP) :: OFFSET_velx  = 12.0_RP
    real(RP) :: OFFSET_vely  = -2.0_RP

    namelist / PARAM_MKINIT_SQUALLLINE / &
       RANDOM_THETA,         &
       OFFSET_velx,          &
       OFFSET_vely

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_squallline",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_squallline",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SQUALLLINE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_squallline",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_squallline",*) 'Not appropriate names in namelist PARAM_MKINIT_SQUALLLINE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SQUALLLINE)

    call read_sounding( RHO, VELX, VELY, POTT, QV1D ) ! (out)

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ( VELX(k) - OFFSET_velx ) * RHO(k)
       MOMY(k,i,j) = ( VELY(k) - OFFSET_vely ) * RHO(k)
       RHOT(k,i,j) = RHO(k) * ( POTT(k) + rndm(k,i,j) * RANDOM_THETA )
       qv  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_squallline

  !-----------------------------------------------------------------------------
  !> Make initial state by Weisman and Klemp (1982)
  subroutine MKINIT_wk1982
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_WK1982 / &
       SFC_THETA, &
       SFC_PRES,  &
       TR_Z,      &
       TR_THETA,  &
       TR_TEMP,   &
       SHEAR_Z,   &
       SHEAR_U,   &
       BBL_THETA

    real(RP) :: rh    (KA,IA,JA)
    real(RP) :: rh_sfc(   IA,JA)

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_wk1982",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_wk1982",*) 'QV is not registered'
       call PRC_abort
    end if

    SFC_PRES  = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WK1982,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_wk1982",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_wk1982",*) 'Not appropriate names in namelist PARAM_MKINIT_WK1982. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_WK1982)

    ! calc in dry condition
    do j = JSB, JEB
    do i = ISB, IEB
       pres_sfc(i,j) = SFC_PRES
       pott_sfc(i,j) = SFC_THETA

       do k = KS, KE
          if ( REAL_CZ(k,i,j) <= TR_Z ) then ! below initial cloud top
             pott(k,i,j) = pott_sfc(i,j) &
                         + ( TR_THETA - pott_sfc(i,j) ) * ( REAL_CZ(k,i,j) / TR_Z )**1.25_RP
          else
             pott(k,i,j) = TR_THETA * exp( GRAV * ( REAL_CZ(k,i,j) - TR_Z ) / CPdry / TR_TEMP )
          endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    ! calc QV from RH
    do j = JSB, JEB
    do i = ISB, IEB
       rh_sfc(i,j) = 1.0_RP - 0.75_RP * ( REAL_FZ(KS-1,i,j) / TR_Z )**1.25_RP

       do k = KS, KE
          if ( REAL_CZ(k,i,j) <= TR_Z ) then ! below initial cloud top
             rh(k,i,j) = 1.0_RP - 0.75_RP * ( REAL_CZ(k,i,j) / TR_Z )**1.25_RP
          else
             rh(k,i,j) = 0.25_RP
          endif
       enddo
    enddo
    enddo

    call SATURATION_psat_all( IA, ISB, IEB, JA, JSB, JEB, &
                              temp_sfc(:,:), & ! [IN]
                              psat_sfc(:,:)  ) ! [OUT]
    qdry(:,:,:) = 1.0_RP - qv(:,:,:) - qc(:,:,:)
    call SATURATION_pres2qsat_all( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                   temp(:,:,:), pres(:,:,:), qdry(:,:,:), & ! [IN]
                                   qsat(:,:,:)                            ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
       qsat_sfc(i,j) = EPSvap * psat_sfc(i,j) / ( pres_sfc(i,j) - ( 1.0_RP-EPSvap ) * psat_sfc(i,j) )
       qv_sfc(i,j) = rh_sfc(i,j) * qsat_sfc(i,j)
       do k = KS, KE
          qv(k,i,j) = rh(k,i,j) * qsat(k,i,j)
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    do k = KS, KE
       LOG_INFO("MKINIT_wk1982",*) k, REAL_CZ(k,IS,JS), pres(k,IS,JS), pott(k,IS,JS), rh(k,IS,JS), qv(k,IS,JS)*1000
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = SHEAR_U * tanh( REAL_CZ(k,i,j) / SHEAR_Z ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = 0.0_RP
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,i,j) * ( pott(k,i,j) + BBL_THETA * bubble(k,i,j) )
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_wk1982

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF01
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP) :: PERTURB_AMP = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0       ! 0 -> no perturbation
                                       ! 1 -> petrurbation for pt
                                       ! 2 -> perturbation for u, v, w
    logical  :: USE_LWSET    = .false. ! use liq. water. static energy temp.?

    namelist / PARAM_MKINIT_RF01 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG,     &
       USE_LWSET

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: sint
    real(RP) :: GEOP_sw ! switch for geopotential energy correction

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP ! pi/2

    LOG_NEWLINE
    LOG_INFO("MKINIT_DYCOMS2_RF01",*) 'Setup initial state'

    rewind(IO_FID_CONF)

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_DYCOMS2_RF01",*) 'QV is not registered'
       call PRC_abort
    end if

    read(IO_FID_CONF,nml=PARAM_MKINIT_RF01,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_DYCOMS2_RF01",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_DYCOMS2_RF01",*) 'Not appropriate names in namelist PARAM_MKINIT_RF01. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RF01)

    if ( USE_LWSET ) then
       GEOP_sw = 1.0_RP
    else
       GEOP_sw = 0.0_RP
    endif

    ! calc in dry condition
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1017.8E2_RP ! [Pa]
       pott_sfc(i,j) = 289.0_RP    ! [K]

       do k = KS, KE
          velx(k,i,j) =   7.0_RP
          vely(k,i,j) =  -5.5_RP
          if ( CZ(k) < 820.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 289.0_RP - GRAV / CPdry * CZ(k) * GEOP_sw
          elseif( CZ(k) <= 860.0_RP ) then
             sint = sin( pi2 * ( CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             potl(k,i,j) = ( 289.0_RP - GRAV / CPdry * CZ(k) * GEOP_sw                          ) * (0.5_RP-sint) &
                         + ( 297.5_RP+sign(abs(CZ(k)-840.0_RP)**(1.0_RP/3.0_RP),CZ(k)-840.0_RP) &
                           - GRAV / CPdry * CZ(k) * GEOP_sw                                     ) * (0.5_RP+sint)
          else
             potl(k,i,j) = 297.5_RP + ( CZ(k)-840.0_RP )**(1.0_RP/3.0_RP) &
                         - GRAV / CPdry * CZ(k) * GEOP_sw
          endif
       enddo

    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    ! calc in moist condition
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc  (i,j) = 9.0E-3_RP   ! [kg/kg]

       do k = KS, KE
          if    ( CZ(k) <   820.0_RP ) then ! below initial cloud top
             qall = 9.0E-3_RP
          elseif( CZ(k) <=  860.0_RP ) then ! boundary
             sint = sin( pi2 * ( CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             qall = 9.0E-3_RP * (0.5_RP-sint) &
                  + 1.5E-3_RP * (0.5_RP+sint)
          elseif( CZ(k) <= 5000.0_RP ) then
             qall = 1.5E-3_RP
          else
             qall = 0.0_RP
          endif

          if    ( CZ(k) <=  600.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( CZ(k) < 820.0_RP ) then ! in the cloud
             fact = ( CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact
          elseif( CZ(k) <= 860.0_RP ) then ! boundary
             sint = sin( pi2 * ( CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             fact = ( CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact * (0.5_RP-sint)
          else
             qc(k,i,j) = 0.0_RP
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 120.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       end if
    enddo
    enddo
    enddo

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF01

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF02
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP) :: PERTURB_AMP  = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0 ! 0 -> no perturbation
                                 ! 1 -> perturbation for PT
                                 ! 2 -> perturbation for u,v,w

    namelist / PARAM_MKINIT_RF02 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: sint

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2
    LOG_NEWLINE
    LOG_INFO("MKINIT_DYCOMS2_RF02",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_DYCOMS2_RF02",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_DYCOMS2_RF02",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_DYCOMS2_RF02",*) 'Not appropriate names in namelist PARAM_MKINIT_RF02. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RF02)

    ! calc in dry condition
    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1017.8E2_RP   ! [Pa]
       pott_sfc(i,j) = 288.3_RP      ! [K]

       do k = KS, KE
          velx(k,i,j) =  3.0_RP + 4.3 * CZ(k)*1.E-3_RP
          vely(k,i,j) = -9.0_RP + 5.6 * CZ(k)*1.E-3_RP

          if ( CZ(k) < 775.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
          else if ( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
             potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP &
                         + ( 295.0_RP+sign(abs(CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),CZ(k)-795.0_RP) ) &
                         * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
          endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    ! calc in moist condition
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc(i,j) = 9.45E-3_RP

       do k = KS, KE
          if ( CZ(k) < 775.0_RP ) then ! below initial cloud top
             qall = 9.45E-3_RP ! [kg/kg]
          else if ( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
             qall = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
                  ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             qall = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ! [kg/kg]
          endif

          if( CZ(k) < 400.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( CZ(k) < 775.0_RP ) then
             fact = ( CZ(k)-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.65E-3_RP * fact
          elseif( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * ( CZ(k)-795.0_RP )/20.0_RP )
             fact = ( CZ(k)-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.65E-3_RP * fact * (1.0_RP-sint) * 0.5_RP
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 55.0E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF02

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF02_DNS
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_RF02_DNS / &
       ZB, CONST_U, CONST_V,PRES_ZB,&
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2

    LOG_NEWLINE
    LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_DYCOMS2_RF02_DNS",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02_DNS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_DYCOMS2_RF02_DNS",*) 'Not appropriate names in namelist PARAM_MKINIT_RF02_DNS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RF02_DNS)

    ! calc in dry condition
    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = PRES_ZB
!      pott_sfc(i,j) = 288.3_RP      ! [K]
!      qv_sfc  (i,j) = 9.45E-3_RP

       do k = KS, KE

          velx(k,i,j) = CONST_U
          vely(k,i,j) = CONST_V

!         if ( ZB+CZ(k) < 775.0_RP ) then ! below initial cloud top
          if ( ZB+CZ(k) <= 795.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
             qall = 9.45E-3_RP ! [kg/kg]
! necessary?
!         else if ( CZ(k) <= 815.0_RP ) then
!            sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
!            potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP + &
!                  ( 295.0_RP+sign(abs(CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),CZ(k)-795.0_RP) ) * (1.0_RP+sint)*0.5_RP
!            qall = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
!                  ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( zb+CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
             qall = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-(zb+CZ(k)))/500.0_RP ) ) ! [kg/kg]
          endif

          if( ZB+CZ(k) < 400.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( ZB+CZ(k) <= 795.0_RP ) then
             fact = ( (zb+CZ(k))-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.8E-3_RP * fact
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall - qc(k,i,j)

          !if(i==is.and.j==js)LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chkk',k,cz(k)+zb,qc(k,i,j),qv(k,i,j)
       enddo
    enddo
    enddo

    !LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chk3',ks,ke
    ! extrapolation (temtative)
    pott_sfc(:,:) = potl(ks,:,:)-0.5*(potl(ks+1,:,:)-potl(ks,:,:))
    qv_sfc  (:,:) = qv  (ks,:,:)-0.5*(qv  (ks+1,:,:)-qv  (ks,:,:))
    qc_sfc  (:,:) = qc  (ks,:,:)-0.5*(qc  (ks+1,:,:)-qc  (ks,:,:))

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    RovCP = Rdry / CPdry
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    !LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chk8'
    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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
    !LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chk9'

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 55.0E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF02_DNS

  !-----------------------------------------------------------------------------
  !> Make initial state for RICO inter comparison
  subroutine MKINIT_RICO
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP):: PERTURB_AMP_PT = 0.1_RP
    real(RP):: PERTURB_AMP_QV = 2.5E-5_RP

    namelist / PARAM_MKINIT_RICO / &
       PERTURB_AMP_PT, &
       PERTURB_AMP_QV

    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall ! QV+QC
    real(RP) :: fact

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_RICO",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_RICO",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RICO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_RICO",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_RICO",*) 'Not appropriate names in namelist PARAM_MKINIT_RICO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RICO)

    ! calc in moist condition
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1015.4E2_RP ! [Pa]
       pott_sfc(i,j) = 297.9_RP

       do k = KS, KE
          !--- potential temperature
          if ( CZ(k) < 740.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 297.9_RP
          else
             fact = ( CZ(k)-740.0_RP ) * ( 317.0_RP-297.9_RP ) / ( 4000.0_RP-740.0_RP )
             potl(k,i,j) = 297.9_RP + fact
          endif

          !--- horizontal wind velocity
          if ( CZ(k) <= 4000.0_RP ) then ! below initial cloud top
             fact = ( CZ(k)-0.0_RP ) * ( -1.9_RP+9.9_RP ) / ( 4000.0_RP-0.0_RP )
             velx(k,i,j) =  -9.9_RP + fact
             vely(k,i,j) =  -3.8_RP
          else
             velx(k,i,j) =  -1.9_RP
             vely(k,i,j) =  -3.8_RP
          endif
       enddo

    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc  (i,j) = 16.0E-3_RP   ! [kg/kg]

       do k = KS, KE
          !--- mixing ratio of vapor
          if ( CZ(k) <= 740.0_RP ) then ! below initial cloud top
             fact = ( CZ(k)-0.0_RP ) * ( 13.8E-3_RP-16.0E-3_RP ) / ( 740.0_RP-0.0_RP )
             qall = 16.0E-3_RP + fact
          elseif ( CZ(k) <= 3260.0_RP ) then ! boundary
             fact = ( CZ(k)-740.0_RP ) * ( 2.4E-3_RP-13.8E-3_RP ) / ( 3260.0_RP-740.0_RP )
             qall = 13.8E-3_RP + fact
          elseif( CZ(k) <= 4000.0_RP ) then
             fact = ( CZ(k)-3260.0_RP ) * ( 1.8E-3_RP-2.4E-3_RP ) / ( 4000.0_RP-3260.0_RP )
             qall = 2.4E-3_RP + fact
          else
             qall = 0.0_RP
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA  ,i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       RHOT(k,i,j) = ( pott(k,i,j)+2.0_RP*( rndm(k,i,j)-0.5_RP )*PERTURB_AMP_PT ) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       qv(k,i,j) = qv(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP_QV
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 70.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo

#endif
    return
  end subroutine MKINIT_RICO

  !-----------------------------------------------------------------------------
  !> Make initial state for BOMEX inter comparison
  subroutine MKINIT_BOMEX
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP):: PERTURB_AMP_PT = 0.1_RP
    real(RP):: PERTURB_AMP_QV = 2.5E-5_RP

    namelist / PARAM_MKINIT_BOMEX / &
       PERTURB_AMP_PT, &
       PERTURB_AMP_QV

    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall ! QV+QC
    real(RP) :: fact

    real(RP) :: qdry, Rtot, CPtot

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_BOMEX",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_BOMEX",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BOMEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_BOMEX",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_BOMEX",*) 'Not appropriate names in namelist PARAM_MKINIT_BOMEX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_BOMEX)

    ! calc in moist condition
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1015.E2_RP ! [Pa]
       pott_sfc(i,j) = 299.1_RP

       do k = KS, KE
          !--- potential temperature
          if ( CZ(k) < 520.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 298.7_RP
          elseif( CZ(k) < 1480.0_RP ) then
             fact = ( CZ(k)-520.0_RP ) * ( 302.4_RP-298.7_RP ) / ( 1480.0_RP-520.0_RP )
             potl(k,i,j) = 298.7_RP + fact
          elseif( CZ(k) < 2000.0_RP ) then
             fact = ( CZ(k)-1480.0_RP ) * ( 308.2_RP-302.4_RP ) / ( 2000.0_RP-1480.0_RP )
             potl(k,i,j) = 302.4_RP + fact
          else
             fact = ( CZ(k)-2000.0_RP ) * 3.65E-3_RP
             potl(k,i,j) = 308.2_RP + fact
          endif

          !--- horizontal wind velocity
          if ( CZ(k) <= 700.0_RP ) then ! below initial cloud top
             velx(k,i,j) =  -8.75_RP
             vely(k,i,j) =   0.0_RP
          else
             fact = 1.8E-3_RP * ( CZ(k)-700.0_RP )
             velx(k,i,j) =  -8.75_RP + fact
             vely(k,i,j) =  0.0_RP
          endif
       enddo

    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc(i,j) = 22.45E-3_RP   ! [kg/kg]

       do k = KS, KE
          !--- mixing ratio of vapor
          if ( CZ(k) <= 520.0_RP ) then ! below initial cloud top
             fact = ( CZ(k)-0.0_RP ) * ( 16.3E-3_RP-17.0E-3_RP ) / ( 520.0_RP-0.0_RP )
             qall = 17.0E-3_RP + fact
          elseif ( CZ(k) <= 1480.0_RP ) then ! boundary
             fact = ( CZ(k)-520.0_RP ) * ( 10.7E-3_RP-16.3E-3_RP ) / ( 1480.0_RP-520.0_RP )
             qall = 16.3E-3_RP + fact
          elseif( CZ(k) <= 2000.0_RP ) then
             fact = ( CZ(k)-1480.0_RP ) * ( 4.2E-3_RP-10.7E-3_RP ) / ( 2000.0_RP-1480.0_RP )
             qall = 10.7E-3_RP + fact
          else
             fact = ( CZ(k)-2000.0_RP ) * ( -1.2E-6_RP )
             qall = 4.2E-3_RP + fact
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       qdry = 1.0_RP - qv(k,i,j) - qc(k,i,j)
       Rtot = Rdry * qdry + Rvap * qv(k,i,j)
       CPtot = CPdry * qdry + CPvap * qv(k,i,j) + CL * qc(k,i,j)
       pott(k,i,j) = ( temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) ) * ( P00 / pres(k,i,j) )**(Rtot/CPtot)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA  ,i,j) = DENS(KE,i,j)
    enddo
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if( CZ(k) <= 1600.0_RP ) then !--- lowest 40 model layer when dz=40m
         RHOT(k,i,j) = ( pott(k,i,j)+2.0_RP*( rndm(k,i,j)-0.5_RP )*PERTURB_AMP_PT ) * DENS(k,i,j)
       else
         RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
       endif
    enddo
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if( CZ(k) <= 1600.0_RP ) then !--- lowest 40 model layer when dz=40m
          qv(k,i,j) = qv(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP_QV
       endif
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 70.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo

#endif
    return
  end subroutine MKINIT_BOMEX

  !-----------------------------------------------------------------------------
  !> Make initial state for BOMEX inter comparison
  subroutine MKINIT_POINTE
    implicit none

    QTRC(:,:,:,:) = 0.0_RP
    DENS(:,:,:) = 1.0_RP
    MOMX(:,:,:) = 0.0_RP
    MOMY(:,:,:) = 0.0_RP
    MOMZ(:,:,:) = 0.0_RP
    pott(:,:,:) = 300.0_RP 
    RHOT(:,:,:) = DENS(:,:,:) * pott(:,:,:)

    return
  end subroutine MKINIT_POINTE
  !-----------------------------------------------------------------------------
  !> Make initial state ( ocean variables )
  subroutine MKINIT_oceancouple
    implicit none

    LOG_NEWLINE
    LOG_INFO("MKINIT_oceancouple",*) 'Setup initial state'

    call flux_setup

    call ocean_setup

    return
  end subroutine MKINIT_oceancouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( land variables )
  subroutine MKINIT_landcouple
    implicit none

    LOG_NEWLINE
    LOG_INFO("MKINIT_landcouple",*) 'Setup initial state'

    call flux_setup

    call land_setup

    return
  end subroutine MKINIT_landcouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( urban variables )
  subroutine MKINIT_urbancouple
    implicit none

    LOG_NEWLINE
    LOG_INFO("MKINIT_urbancouple",*) 'Setup initial state'

    call flux_setup

    call urban_setup

    return
  end subroutine MKINIT_urbancouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( sea breeze )
  subroutine MKINIT_seabreeze
    use scale_prc_cartesC, only: &
       PRC_NUM_X
    use scale_landuse, only: &
       LANDUSE_frac_land, &
       LANDUSE_calc_fact, &
       LANDUSE_fillhalo,  &
       LANDUSE_write
    use scale_atmos_grid_cartesC, only: &
       DOMAIN_CENTER_X => ATMOS_GRID_CARTESC_DOMAIN_CENTER_X
    implicit none

    real(RP) :: LAND_SIZE

    namelist / PARAM_MKINIT_SEABREEZE / &
       LAND_SIZE

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_seabreeze",*) 'Setup initial state'

    LAND_SIZE = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SEABREEZE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_seabreeze",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_seabreeze",*) 'Not appropriate names in namelist PARAM_MKINIT_SEABREEZE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SEABREEZE)

    call flux_setup

    call land_setup

    call ocean_setup

    ! make landuse conditions
    do j = JSB, JEB
    do i = ISB, IEB
       if ( abs( CX(i) - DOMAIN_CENTER_X ) < LAND_SIZE ) then
          LANDUSE_frac_land(i,j) = 1.0_RP
       else
          LANDUSE_frac_land(i,j) = 0.0_RP
       endif
    enddo
    enddo

    ! calculate landuse factors
    call LANDUSE_fillhalo( FILL_BND=.true. )
    call LANDUSE_calc_fact

    ! output landuse file
    call LANDUSE_write

    return
  end subroutine MKINIT_seabreeze

  !-----------------------------------------------------------------------------
  !> Make initial state ( heat island )
  subroutine MKINIT_heatisland
    use scale_prc_cartesC, only: &
       PRC_NUM_X
    use scale_landuse, only: &
       LANDUSE_frac_land,  &
       LANDUSE_frac_urban, &
       LANDUSE_calc_fact,  &
       LANDUSE_fillhalo,   &
       LANDUSE_write
    implicit none

    real(RP) :: dist

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_heatisland",*) 'Setup initial state'

    call flux_setup

    call land_setup

    call urban_setup

    ! 1/9 size of domain
    dist = ( CXG(IMAX*PRC_NUM_X) - CXG(1) ) / 9.0_RP

    ! make landuse conditions
    do j = JSB, JEB
    do i = ISB, IEB
       if (       CX(i) >= dist * 4.0_RP &
            .AND. CX(i) <  dist * 5.0_RP ) then
          LANDUSE_frac_land(i,j)  = 1.0_RP
          LANDUSE_frac_urban(i,j) = 1.0_RP
       else
          LANDUSE_frac_land(i,j)  = 1.0_RP
          LANDUSE_frac_urban(i,j) = 0.0_RP
       endif
    enddo
    enddo

    ! calculate landuse factors
    call LANDUSE_fillhalo( FILL_BND=.true. )
    call LANDUSE_calc_fact

    ! output landuse file
    call LANDUSE_write

    return
  end subroutine MKINIT_heatisland

  !-----------------------------------------------------------------------------
  !> Make initial state for grayzone experiment
  subroutine MKINIT_grayzone
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)

    real(RP) :: PERTURB_AMP = 0.0_RP
    integer  :: RANDOM_LIMIT = 0
    integer  :: RANDOM_FLAG  = 0       ! 0 -> no perturbation
                                       ! 1 -> petrurbation for pt
                                       ! 2 -> perturbation for u, v, w

    namelist / PARAM_MKINIT_GRAYZONE / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_grayzone",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_grayzone",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_GRAYZONE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_grayzone",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_grayzone",*) 'Not appropriate names in namelist PARAM_MKINIT_GRAYZONE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_GRAYZONE)

    call read_sounding( RHO, VELX, VELY, POTT, QV1D ) ! (out)

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
       qv  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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

    call RANDOM_uniform(rndm) ! make random
    do j = JSB, JEB
    do i = ISB, IEB
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
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       CVdry => CONST_CVdry, &
       CVvap => CONST_CVvap, &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_rhot2temp_pres
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE
    implicit none

    real(RP) :: init_dens  = 1.12_RP   ![kg/m3]
    real(RP) :: init_temp  = 298.18_RP ![K]
    real(RP) :: init_pres  = 1.E+5_RP  ![Pa]
    real(RP) :: init_ssliq = 0.01_RP   ![%]

    namelist / PARAM_MKINIT_BOXAERO / &
       init_dens, &
       init_temp, &
       init_pres, &
       init_ssliq

    real(RP) :: rtot (KA,IA,JA)
    real(RP) :: cvtot(KA,IA,JA)
    real(RP) :: cptot(KA,IA,JA)
    real(RP) :: qdry
    real(RP) :: psat, qsat
    integer  :: i, j, k, ierr
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_AE_TYPE /= 'KAJINO13' ) then
       LOG_INFO("MKINIT_boxaero",*) 'For [Box model of aerosol],'
       LOG_INFO("MKINIT_boxaero",*) 'ATMOS_PHY_AE_TYPE should be KAJINO13. Stop! ', trim(ATMOS_PHY_AE_TYPE)
       call PRC_abort
    endif

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_boxaero",*) 'QV is not registered'
       call PRC_abort
    end if

    LOG_NEWLINE
    LOG_INFO("MKINIT_boxaero",*) 'Setup initial state'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BOXAERO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_boxaero",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_boxaero",*) 'Not appropriate names in namelist PARAM_MKINIT_BOXAERO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_BOXAERO)

    call SATURATION_psat_all( init_temp, psat )
    qsat = EPSvap * psat / ( init_pres - ( 1.0_RP-EPSvap ) * psat )

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS(k,i,j) = init_dens
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       MOMZ(k,i,j) = 0.0_RP
       pott(k,i,j) = init_temp * ( P00/init_pres )**(Rdry/CPdry)
       RHOT(k,i,j) = init_dens * pott(k,i,j)

       qv(k,i,j) = ( init_ssliq + 1.0_RP ) * qsat

       qdry = 1.0 - qv(k,i,j)
       rtot (k,i,j) = Rdry  * qdry + Rvap  * qv(i,i,j)
       cvtot(k,i,j) = CVdry * qdry + CVvap * qv(i,i,j)
       cptot(k,i,j) = CPdry * qdry + CPvap * qv(i,i,j)
    enddo
    enddo
    enddo

    call ATMOS_THERMODYN_rhot2temp_pres( KA, 1, KA, IA, 1, IA, JA, 1, JA, &
                                         dens(:,:,:), RHOT(:,:,:),                & ! (in)
                                         rtot(:,:,:), cvtot(:,:,:), cptot(:,:,:), & ! (in)
                                         temp(:,:,:), pres(:,:,:)                 ) ! (out)

    return
  end subroutine MKINIT_boxaero

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  subroutine MKINIT_warmbubbleaero
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
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

    namelist / PARAM_MKINIT_WARMBUBBLE / &
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

    LOG_NEWLINE
    LOG_INFO("MKINIT_warmbubbleaero",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_warmbubbleaero",*) 'QV is not registerd'
       call PRC_abort
    end if


    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_warmbubbleaero",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_warmbubbleaero",*) 'Not appropriate names in namelist PARAM_MKINIT_WARMBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA

    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( CZ(k)-CZ(k-1) )
       else                                ! Layer 3
          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( CZ(k)-CZ(k-1) )
       endif
    enddo

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    ! calc QV from RH
    call SATURATION_psat_all( temp_sfc(1,1), psat_sfc(1,1) ) ! [IN], [OUT]
    qsat_sfc(1,1) = EPSvap * psat_sfc(1,1) / ( pres_sfc(1,1) - ( 1.0_RP-EPSvap ) * psat_sfc(1,1) )

    qdry(:,1,1) = 1.0_RP - qv(:,1,1) - qc(:,1,1)
    call SATURATION_pres2qsat_all( KA, KS, KE, &
                                   temp(:,1,1), pres(:,1,1), qdry(:,1,1), & ! [IN]
                                   qsat(:,1,1)                            ) ! [OUT]
    qv_sfc(1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1)
    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       elseif( CZ(k) <= ENV_L2_ZTOP ) then ! Layer 2
          qv(k,1,1) = ENV_RH * 1.E-2_RP * qsat(k,1,1)
       else                                ! Layer 3
          qv(k,1,1) = 0.0_RP
       endif
    enddo

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1)    ) ! [OUT]

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       qv  (k,i,j) = qv(k,1,1)
    enddo
    enddo
    enddo

    call flux_setup

    return
  end subroutine MKINIT_warmbubbleaero

  !-----------------------------------------------------------------------------
  !> Make initial state ( real case )
  subroutine MKINIT_real
    use mod_realinput, only: &
         REALINPUT_atmos, &
         REALINPUT_surface
    implicit none

    call PROF_rapstart('__Real_Atmos',2)

    call REALINPUT_atmos

    call PROF_rapend  ('__Real_Atmos',2)
    call PROF_rapstart('__Real_Surface',2)

    call REALINPUT_surface

    call PROF_rapend  ('__Real_Surface',2)

    call flux_setup

    return
  end subroutine MKINIT_real

end module mod_mkinit
