!-------------------------------------------------------------------------------
!> module initial
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new] Imported from SCALE-LES ver.2
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
!!
!<
!-------------------------------------------------------------------------------
module mod_mkinit
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_get_available_fid, &
     IO_SYSCHR,            &
     IO_FILECHR,           &
     IO_FID_LOG,  &
     IO_FID_CONF, &
     IO_L
  use mod_process, only: &
     PRC_MPIstop
  use mod_const, only : &
     PI    => CONST_PI,    &
     Pstd  => CONST_Pstd,  &
     CPdry => CONST_CPdry, &
     RovCP => CONST_RovCP, &
     LH0   => CONST_LH0,   &
     P00   => CONST_PRE00
  use mod_random, only: &
     RANDOM_get
  use mod_grid, only : &
     CZ => GRID_CZ, &
     CX => GRID_CX, &
     CY => GRID_CY
  use mod_atmos_vars, only: &
     DENS, &
     MOMX, &
     MOMY, &
     MOMZ, &
     RHOT, &
     QTRC
  use mod_atmos_hydrostatic, only: &
     hydro_buildrho      => ATMOS_HYDRO_buildrho,     &
     hydro_buildrho_1d   => ATMOS_HYDRO_buildrho_1d,  &
     hydro_buildrho_temp => ATMOS_HYDRO_buildrho_temp
  use mod_atmos_saturation, only: &
     SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKINIT_setup
  public :: MKINIT_planestate
  public :: MKINIT_tracerbubble
  public :: MKINIT_coldbubble
  public :: MKINIT_warmbubble
  public :: MKINIT_khwave
  public :: MKINIT_turbulence
  public :: MKINIT_supercell
  public :: MKINIT_squallline
  public :: MKINIT_DYCOMS2_RF01
  public :: MKINIT_DYCOMS2_RF02
  public :: MKINIT_DYCOMS2_RF01_hbinw
  public :: MKINIT_warmbubble_hbinw

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, save      :: MKINIT_TYPE
  integer, public, parameter :: I_PLANESTATE          =  1
  integer, public, parameter :: I_TRACERBUBBLE        =  2
  integer, public, parameter :: I_COLDBUBBLE          =  3
  integer, public, parameter :: I_WARMBUBBLE          =  4
  integer, public, parameter :: I_KHWAVE              =  5
  integer, public, parameter :: I_TURBULENCE          =  6
  integer, public, parameter :: I_SUPERCELL           =  7
  integer, public, parameter :: I_SQUALLLINE          =  8
  integer, public, parameter :: I_DYCOMS2_RF01        =  9
  integer, public, parameter :: I_DYCOMS2_RF02        = 10
  integer, public, parameter :: I_DYCOMS2_RF01_hbinw  = 11
  integer, public, parameter :: I_WARMBUBBLE_hbinw    = 12
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: BUBBLE_setup
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: THETAstd = 300.0_RP ! [K]

  real(RP), private :: pres(KA,IA,JA) ! pressure [Pa]
  real(RP), private :: temp(KA,IA,JA) ! temperature [K]
  real(RP), private :: pott(KA,IA,JA) ! potential temperature [K]
  real(RP), private :: qsat(KA,IA,JA) ! satulated water vapor [kg/kg]
  real(RP), private :: qv  (KA,IA,JA) ! water vapor [kg/kg]
  real(RP), private :: qc  (KA,IA,JA) ! cloud water [kg/kg]
  real(RP), private :: velx(KA,IA,JA) ! velocity u [m/s]
  real(RP), private :: vely(KA,IA,JA) ! velocity v [m/s]

  real(RP), private :: pres_sfc(1,IA,JA)
  real(RP), private :: temp_sfc(1,IA,JA)
  real(RP), private :: pott_sfc(1,IA,JA)
  real(RP), private :: qsat_sfc(1,IA,JA)
  real(RP), private :: qv_sfc  (1,IA,JA)
  real(RP), private :: qc_sfc  (1,IA,JA)

  real(RP), private :: rndm  (KA,IA,JA) ! random number (0-1)
  real(RP), private :: bubble(KA,IA,JA) ! bubble factor (0-1)
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine MKINIT_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF8
    implicit none

    character(len=IO_SYSCHR) :: MKINIT_initname = 'COLDBUBBLE'

    NAMELIST / PARAM_MKINIT / &
       MKINIT_initname

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[MKINIT]/Categ[MKINIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Check!'
       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       pres(k,i,j) = CONST_UNDEF8
       temp(k,i,j) = CONST_UNDEF8
       pott(k,i,j) = CONST_UNDEF8
       qsat(k,i,j) = CONST_UNDEF8
       qv  (k,i,j) = CONST_UNDEF8
       qc  (k,i,j) = CONST_UNDEF8
       velx(k,i,j) = CONST_UNDEF8
       vely(k,i,j) = CONST_UNDEF8

       rndm  (k,i,j) = CONST_UNDEF8
       bubble(k,i,j) = CONST_UNDEF8
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       pres_sfc(1,i,j) = CONST_UNDEF8
       temp_sfc(1,i,j) = CONST_UNDEF8
       pott_sfc(1,i,j) = CONST_UNDEF8
       qsat_sfc(1,i,j) = CONST_UNDEF8
       qv_sfc  (1,i,j) = CONST_UNDEF8
       qc_sfc  (1,i,j) = CONST_UNDEF8
    enddo
    enddo

    select case(trim(MKINIT_initname))
    case('PLANESTATE')
       MKINIT_TYPE = I_PLANESTATE
    case('TRACERBUBBLE')
       MKINIT_TYPE = I_TRACERBUBBLE
       call BUBBLE_setup
    case('COLDBUBBLE')
       MKINIT_TYPE = I_COLDBUBBLE
       call BUBBLE_setup
    case('WARMBUBBLE')
       MKINIT_TYPE = I_WARMBUBBLE
       call BUBBLE_setup
    case('KHWAVE')
       MKINIT_TYPE = I_KHWAVE
    case('TURBULENCE')
       MKINIT_TYPE = I_TURBULENCE
    case('SUPERCELL')
       MKINIT_TYPE = I_SUPERCELL
       call BUBBLE_setup
    case('SQUALLLINE')
       MKINIT_TYPE = I_SQUALLLINE
    case('DYCOMS2_RF01')
       MKINIT_TYPE = I_DYCOMS2_RF01
    case('DYCOMS2_RF02')
       MKINIT_TYPE = I_DYCOMS2_RF02
    case('DYCOMS2_RF01_hbinw')
       MKINIT_TYPE = I_DYCOMS2_RF01_hbinw
    case('WARMBUBBLE_hbinw')
       MKINIT_TYPE = I_WARMBUBBLE_hbinw
    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(MKINIT_initname)
       call PRC_MPIstop
    endselect

    return
  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine BUBBLE_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
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
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[BUBBLE]/Categ[MKINIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Check!'
       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_BUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_BUBBLE)

    if ( BBL_eachnode ) then
       CZ_offset = CZ(KS)
       CX_offset = CX(IS)
       CY_offset = CY(JS)
    else
       CZ_offset = 0.0_RP
       CX_offset = 0.0_RP
       CY_offset = 0.0_RP
    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! make tracer bubble
       dist = ( (CZ(k)-CZ_offset-BBL_CZ)/BBL_RZ )**2 &
            + ( (CX(i)-CX_offset-BBL_CX)/BBL_RX )**2 &
            + ( (CY(j)-CY_offset-BBL_CY)/BBL_RY )**2

       bubble(k,i,j) = cos( 0.5_RP*PI*sqrt( min(dist,1.0_RP) ) )**2

    enddo
    enddo
    enddo

    return
  end subroutine BUBBLE_setup

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform + random disturbance )
  !-----------------------------------------------------------------------------
  subroutine MKINIT_planestate
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA              ! surface potential temperature [K]
    real(RP) :: SFC_PRES               ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  0.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_THETA              ! potential temperature of environment [K]
    real(RP) :: ENV_W        =  0.0_RP ! velocity w of environment [m/s]
    real(RP) :: ENV_U        =  0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =  0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =  0.0_RP ! relative humidity of environment [%]
    ! Disturbance
    real(RP) :: RANDOM_THETA =  0.0_RP ! amplitude of random disturbance theta
    real(RP) :: RANDOM_W     =  0.0_RP ! amplitude of random disturbance w
    real(RP) :: RANDOM_U     =  0.0_RP ! amplitude of random disturbance u
    real(RP) :: RANDOM_V     =  0.0_RP ! amplitude of random disturbance v
    real(RP) :: RANDOM_RH    =  0.0_RP ! amplitude of random disturbance RH

    NAMELIST / PARAM_MKINIT_PLANESTATE / &
       SFC_THETA,    &
       SFC_PRES,     &
       SFC_RH,       &
       ENV_THETA,    &
       ENV_W,        &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       RANDOM_THETA, &
       RANDOM_W,     &
       RANDOM_U,     &
       RANDOM_V,     &
       RANDOM_RH

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Horiz_UNIFORM]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_PLANESTATE)

    ! calc in dry condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA + rndm(KS-1,i,j) * RANDOM_THETA
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + rndm(k,i,j) * RANDOM_THETA
          qv  (k,i,j) = 0.0_RP
          qc  (k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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
    call SATURATION_pres2qsat_liq( qsat_sfc(1,:,:), temp_sfc(1,:,:), pres_sfc(1,:,:) )
    call SATURATION_pres2qsat_liq( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(1,i,j)

       do k = KS, KE
          qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat(k,i,j)
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = ( ENV_W + ( rndm(k,i,j) - 0.5_RP ) * 2.0_RP * RANDOM_W ) &
                   * 0.5_RP * ( DENS(k,i,j) + DENS(k+1,i,j) )
    enddo
    enddo
    enddo

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
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
    enddo
    enddo
    enddo

    if ( QA >= 2 .and. QA <= 11 ) then
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

    return
  end subroutine MKINIT_planestate
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Make initial state for tracer bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_tracerbubble
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA         ! surface potential temperature [K]
    real(RP) :: SFC_PRES          ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_THETA         ! potential temperature of environment [K]
    real(RP) :: ENV_U     =  0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V     =  0.0_RP ! velocity v of environment [m/s]
    ! Bubble
    real(RP) :: BBL_NC    =  1.0_RP ! extremum of NC in bubble [kg/kg]

    NAMELIST / PARAM_MKINIT_TRACERBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       ENV_U,     &
       ENV_V,     &
       BBL_NC

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TRACERBUBBLE]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_TRACERBUBBLE)

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
    call hydro_buildrho_1d( DENS    (:,1,1), & ! [OUT]
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
  !-----------------------------------------------------------------------------
  subroutine MKINIT_coldbubble
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA         ! surface potential temperature [K]
    real(RP) :: SFC_PRES          ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_THETA         ! potential temperature of environment [K]
    ! Bubble
    real(RP) :: BBL_TEMP = -15.0_RP  ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_COLDBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       BBL_TEMP

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[COLDBUBBLE]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_COLDBUBBLE)

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
    call hydro_buildrho_1d( DENS    (:,1,1), & ! [OUT]
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
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo

       ! make cold bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1)                                           &
                                   + BBL_TEMP * ( P00/pres(k,1,1) )**RovCP * bubble(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_coldbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_warmbubble
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  80.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_RH       =  80.0_RP ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =  1.E3_RP ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  = 12.E3_RP ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS = 4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS = 3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(RP) :: BBL_THETA    =   1.0_RP ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       ENV_RH,       &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    pres_sfc(1,1,1) = SFC_PRES
    pott_sfc(1,1,1) = SFC_THETA
    qv_sfc  (1,1,1) = 0.0_RP
    qc_sfc  (1,1,1) = 0.0_RP

    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( CZ(k)-CZ(k-1) )
       else                                ! Layer 3
          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( CZ(k)-CZ(k-1) )
       endif
       qv(k,1,1) = 0.0_RP
       qc(k,1,1) = 0.0_RP
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho_1d( DENS    (:,1,1), & ! [OUT]
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
    call SATURATION_pres2qsat_liq( qsat_sfc(1,1,1), temp_sfc(1,1,1), pres_sfc(1,1,1) )
    call SATURATION_pres2qsat_liq( qsat    (:,1,1), temp    (:,1,1), pres    (:,1,1) )

    qv_sfc(1,1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1,1)
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
    call hydro_buildrho_1d( DENS    (:,1,1), & ! [OUT]
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
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP

       QTRC(k,i,j,I_QV) = qv(k,1,1)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_warmbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for Kelvin Helmholtz experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_khwave
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH         =  0.0_RP  ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_L1_ZTOP    = 1.95E3_RP ! top    height of the layer1 (low  THETA) [m]
    real(RP) :: ENV_L3_ZBOTTOM = 2.05E3_RP ! bottom height of the layer3 (high THETA) [m]
    real(RP) :: ENV_L1_THETA   = 300.0_RP ! THETA in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_THETA   = 305.0_RP ! THETA in the layer3 (high THETA) [K]
    real(RP) :: ENV_L1_U       =   0.0_RP ! velocity u in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_U       =  20.0_RP ! velocity u in the layer3 (high THETA) [K]
    real(RP) :: ENV_L1_RH      =  50.0_RP ! Relative Humidity in the layer1 (low  THETA) [%]
    real(RP) :: ENV_L3_RH      =   0.0_RP ! Relative Humidity in the layer3 (high THETA) [%]

    NAMELIST / PARAM_MKINIT_KHWAVE / &
       SFC_THETA,      &
       SFC_PRES,       &
       ENV_L1_ZTOP,    &
       ENV_L3_ZBOTTOM, &
       ENV_L1_THETA,   &
       ENV_L3_THETA,   &
       ENV_L1_U,       &
       ENV_L3_U,       &
       ENV_L1_RH,      &
       ENV_L3_RH

    real(RP) :: fact
    
    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[KH wave]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_KHWAVE)

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if ( CZ(k) <= ENV_L1_ZTOP ) then       ! Layer 1
             pott(k,i,j) = ENV_L1_THETA
          elseif( CZ(k) <= ENV_L3_ZBOTTOM ) then ! Layer 3
             pott(k,i,j) = ENV_L3_THETA
          else                                   ! Layer 2
             fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

             pott(k,i,j) = ENV_L1_THETA * ( 1.0_RP - fact ) &
                         + ENV_L3_THETA * (        fact )
          endif
          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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
    call SATURATION_pres2qsat_liq( qsat_sfc(1,:,:), temp_sfc(1,:,:), pres_sfc(1,:,:) )
    call SATURATION_pres2qsat_liq( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = SFC_RH * 1.E-2_RP * qsat_sfc(1,i,j)

       do k = KS, KE
          if ( CZ(k) <= ENV_L1_ZTOP ) then    ! Layer 1
             qv(k,i,j) = ENV_L1_RH * 1.E-2_RP * qsat_sfc(k,i,j)
          elseif( CZ(k) <= ENV_L3_ZBOTTOM ) then ! Layer 3
             qv(k,i,j) = ENV_L3_RH * 1.E-2_RP * qsat_sfc(k,i,j)
          else                                ! Layer 2
             fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

             qv(k,i,j) = ( ENV_L1_RH * ( 1.0_RP - fact ) &
                         + ENV_L3_RH * (        fact ) ) * 1.E-2_RP * qsat_sfc(k,i,j)
          endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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

       MOMZ(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo

       fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

       MOMX(k,i,j) = ( ENV_L1_U * ( 1.0_RP - fact ) &
                     + ENV_L3_U * (        fact ) ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_khwave

  !-----------------------------------------------------------------------------
  !> Make initial state for turbulence experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_turbulence
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA                ! surface potential temperature [K]
    real(RP) :: SFC_PRES                 ! surface pressure [Pa]
    real(RP) :: SFC_RH       = 50.0_RP   ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_THETA                ! potential temperature of environment [K]
    real(RP) :: ENV_TLAPS    =  4.E-3_RP ! Lapse rate of THETA [K/m]
    real(RP) :: ENV_U        =  5.0_RP   ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =  0.0_RP   ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       = 50.0_RP   ! relative humidity of environment [%]
    ! Disturbance
    real(RP) :: RANDOM_THETA =  3.0_RP   ! amplitude of random disturbance theta
    real(RP) :: RANDOM_U     =  5.E-2_RP ! amplitude of random disturbance u
    real(RP) :: RANDOM_V     =  0.0_RP   ! amplitude of random disturbance v
    real(RP) :: RANDOM_RH    =  5.E-1_RP ! amplitude of random disturbance RH

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
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TURBULENCE]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_TURBULENCE)

    ! calc in dry condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA 
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * CZ(k) 
          qv  (k,i,j) = 0.0_RP
          qc  (k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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
    call SATURATION_pres2qsat_liq( qsat_sfc(1,:,:), temp_sfc(1,:,:), pres_sfc(1,:,:) )
    call SATURATION_pres2qsat_liq( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(1,i,j)

       do k = KS, KE
          qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.E-2_RP * qsat(k,i,j)
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = pott(k,i,j) + rndm(k,i,j) * RANDOM_THETA
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

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
    enddo
    enddo
    enddo

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

    return
  end subroutine MKINIT_turbulence

  !-----------------------------------------------------------------------------
  !> Make initial state for supercell experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_supercell
    implicit none

    character(len=IO_FILECHR) :: ENV_IN_SOUNDING_file = ''
    ! Bubble
    real(RP) :: BBL_THETA = 3.D0 ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_SUPERCELL / &
       ENV_IN_SOUNDING_file, &
       BBL_THETA

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(RP) :: SFC_THETA          ! surface potential temperature [K]
    real(RP) :: SFC_PRES           ! surface pressure [Pa]
    real(RP) :: SFC_QV             ! surface watervapor [g/kg]

    real(RP) :: EXP_z   (EXP_klim) ! height      [m]
    real(RP) :: EXP_pott(EXP_klim) ! potential temperature [K]
    real(RP) :: EXP_qv  (EXP_klim) ! water vapor [g/kg]
    real(RP) :: EXP_u   (EXP_klim) ! velocity u  [m/s]
    real(RP) :: EXP_v   (EXP_klim) ! velocity v  [m/s]

    real(RP) :: fact1, fact2

    integer :: ierr, fid
    integer :: k, i, j, iq, kref
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SUPERCELL]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SUPERCELL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SUPERCELL. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_SUPERCELL)

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

    EXP_z   (1) = 0.0_RP
    EXP_pott(1) = SFC_THETA
    EXP_qv  (1) = SFC_QV
    EXP_u   (1) = EXP_u(2)
    EXP_v   (1) = EXP_v(2)
    do k = 1, EXP_klim
       EXP_qv(k) = EXP_qv(k) * 1.E-3_RP ![g/kg] -> [kg/kg]
    enddo

    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = SFC_QV
       qc_sfc  (1,i,j) = 0.0_RP
    enddo
    enddo

    !--- linear interpolate to model grid
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       qc(k,i,j) = 0.0_RP

       do kref = 2, EXP_kmax
          if (       CZ(k) >  EXP_z(kref-1) &
               .AND. CZ(k) <= EXP_z(kref)   ) then

             fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

             pott(k,i,j) = EXP_pott(kref-1) * fact1 &
                         + EXP_pott(kref)   * fact2
             velx(k,i,j) = EXP_u   (kref-1) * fact1 &
                         + EXP_u   (kref)   * fact2
             vely(k,i,j) = EXP_v   (kref-1) * fact1 &
                         + EXP_v   (kref)   * fact2
             qv  (k,i,j) = EXP_qv  (kref-1) * fact1 &
                         + EXP_qv  (kref)   * fact2
          endif
       enddo
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,i,j) * ( pott(k,i,j) + BBL_THETA * bubble(k,i,j) )
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_supercell

  !-----------------------------------------------------------------------------
  !> Make initial state for squallline experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_squallline
    implicit none

    character(len=IO_FILECHR) :: ENV_IN_SOUNDING_file = ''

    NAMELIST / PARAM_MKINIT_SQUALLLINE / &
       ENV_IN_SOUNDING_file

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(RP) :: SFC_THETA          ! surface potential temperature [K]
    real(RP) :: SFC_PRES           ! surface pressure [Pa]
    real(RP) :: SFC_QV             ! surface watervapor [g/kg]

    real(RP) :: EXP_z   (EXP_klim) ! height      [m]
    real(RP) :: EXP_pott(EXP_klim) ! potential temperature [K]
    real(RP) :: EXP_qv  (EXP_klim) ! water vapor [g/kg]
    real(RP) :: EXP_u   (EXP_klim) ! velocity u  [m/s]
    real(RP) :: EXP_v   (EXP_klim) ! velocity v  [m/s]

    real(RP) :: fact1, fact2

    integer :: ierr, fid
    integer :: k, i, j, iq, kref
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SQUALLLINE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SQUALLLINE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_SQUALLLINE)

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

    EXP_z   (1) = 0.0_RP
    EXP_pott(1) = SFC_THETA
    EXP_qv  (1) = SFC_QV
    EXP_u   (1) = EXP_u(2)
    EXP_v   (1) = EXP_v(2)
    do k = 1, EXP_klim
       EXP_qv(k) = EXP_qv(k) * 1.E-3_RP ![g/kg] -> [kg/kg]
    enddo

    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = SFC_QV
       qc_sfc  (1,i,j) = 0.0_RP
    enddo
    enddo

    !--- linear interpolate to model grid
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       qc(k,i,j) = 0.0_RP

       do kref = 2, EXP_kmax
          if (       CZ(k) >  EXP_z(kref-1) &
               .AND. CZ(k) <= EXP_z(kref)   ) then

             fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

             pott(k,i,j) = EXP_pott(kref-1) * fact1 &
                         + EXP_pott(kref)   * fact2
             velx(k,i,j) = EXP_u   (kref-1) * fact1 &
                         + EXP_u   (kref)   * fact2
             vely(k,i,j) = EXP_v   (kref-1) * fact1 &
                         + EXP_v   (kref)   * fact2
             qv  (k,i,j) = EXP_qv  (kref-1) * fact1 &
                         + EXP_qv  (kref)   * fact2
          endif
       enddo
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS    (:,:,:), & ! [OUT]
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
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.0_RP
       enddo

       ! make warm bubble
    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_squallline

  !-----------------------------------------------------------------------------
  !> Make initial state for strato cumulus
  !-----------------------------------------------------------------------------
  subroutine MKINIT_stratocumulus
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_L1_ZTOP    = 840.0_RP  ! top    height of the layer1 (low  THETA) [m]
    real(RP) :: ENV_L3_ZBOTTOM = 840.0_RP  ! bottom height of the layer3 (high THETA) [m]
    real(RP) :: ENV_L1_THETA   = 289.0_RP  ! THETA in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_THETA   = 297.5_RP  ! THETA in the layer3 (high THETA) [K]
    real(RP) :: ENV_L1_QV      =   8.5E-3_RP ! Specific Humidity in the layer1 (low  THETA) [kg/kg]
    real(RP) :: ENV_L3_QV      =   1.0E-3_RP ! Specific Humidity in the layer3 (high THETA) [kg/kg]

    real(RP) :: ENV_CL_ZBOTTOM = 600.0_RP  ! bottom height of the cloud layer [m]
    real(RP) :: ENV_CL_ZTOP    = 840.0_RP  ! top    height of the cloud layer [m]
    real(RP) :: ENV_CL_QC      =   0.5E-3_RP ! cloud water mixing ratio in the cloud layer   [kg/kg]
    real(RP) :: ENV_CL_NC      = 120.0E6_RP  ! cloud number concentration in the cloud layer [1/m3]

    real(RP) :: ENV_U          =   7.0_RP  ! velocity u in the layer1 (low  THETA) [K]
    real(RP) :: ENV_V          =  -5.5_RP  ! velocity u in the layer3 (high THETA) [K]

    real(RP) :: RANDOM_AMP     =   1.0E-2_RP ! ratio of random disturbance [0-1]

    NAMELIST / PARAM_MKINIT_STRATOCUMULUS / &
       SFC_THETA,      &
       SFC_PRES,       &
       ENV_L1_ZTOP,    &
       ENV_L3_ZBOTTOM, &
       ENV_L1_THETA,   &
       ENV_L3_THETA,   &
       ENV_L1_QV,      &
       ENV_L3_QV,      &
       ENV_CL_ZBOTTOM, &
       ENV_CL_ZTOP,    &
       ENV_CL_QC,      &
       ENV_CL_NC,      &
       ENV_U,          &
       ENV_V,          &
       RANDOM_AMP

    real(RP) :: fact, disturb

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[STRATOCUMULUS)]/Categ[MKINIT]'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_STRATOCUMULUS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_STRATOCUMULUS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_STRATOCUMULUS)

    ! calc in moist condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       disturb = ( 1.0_RP + 2.0_RP * ( rndm(KS-1,i,j)-0.5_RP ) * RANDOM_AMP )

       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA * disturb
       qv_sfc  (1,i,j) = ENV_L1_QV

       do k = KS, KE
          disturb = ( 1.0_RP + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * RANDOM_AMP )

          if ( CZ(k) <= ENV_L1_ZTOP ) then       ! Layer 1
             pott(k,i,j) = ENV_L1_THETA * disturb
             qv  (k,i,j) = ENV_L1_QV
          elseif( CZ(k) >= ENV_L3_ZBOTTOM ) then ! Layer 3
             pott(k,i,j) = ENV_L3_THETA * disturb &
                         + ( CZ(k) - ENV_L3_ZBOTTOM )**(1.0_RP/3.0_RP) ! increase with height
             qv  (k,i,j) = ENV_L3_QV
          else                                   ! Layer 2
             fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

             pott(k,i,j) = ( ENV_L1_THETA * ( 1.0_RP - fact ) &
                           + ENV_L3_THETA * (        fact ) ) * disturb
             qv  (k,i,j) = ENV_L1_QV    * ( 1.0_RP - fact ) &
                         + ENV_L3_QV    * (        fact )
          endif
          qc(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U * ( 1.0_RP + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * RANDOM_AMP ) ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V * ( 1.0_RP + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * RANDOM_AMP ) ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       QTRC(k,i,j,I_QV) = qv(k,i,j)

       if ( CZ(k) >= ENV_CL_ZBOTTOM .and. CZ(k) <= ENV_CL_ZTOP ) then
          QTRC(k,i,j,I_QC) = ENV_CL_QC
          QTRC(k,i,j,I_NC) = ENV_CL_NC / DENS(k,i,j)
       endif

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_stratocumulus

  !-----------------------------------------------------------------------------
  !> Make initial state for strato cumulus
  !-----------------------------------------------------------------------------
  subroutine MKINIT_DYCOMS2_RF01
    implicit none

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall(KA,IA,JA) ! QV+QC
    real(RP) :: fact

    real(RP) :: pi2 
    real(RP) :: sint

    integer :: ierr
    integer :: k, i, j, iq
    real(RP) :: PERTURB_AMP = 0.0_RP
    integer :: RANDOM_LIMIT = 5
    integer :: RANDOM_FLAG = 0  !- 0 -> no perturbation
                                !- 1 -> petrurbation for pt
                                !- 2 -> perturbation for u, v, w
    real(RP) :: dummy(KA,IA,JA)

    NAMELIST / PARAM_MKINIT_RF01 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG
    !---------------------------------------------------------------------------

    dummy(:,:,:) = 0.d0
    pi2 = atan(1.0_RP) * 2.0_RP ! pi/2
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DYCOMS2_RF01)]/Categ[MKINIT]'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF01,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_RF01. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_RF01)

    call RANDOM_get(rndm) ! make random

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1017.8E2_RP ! [Pa]
       pott_sfc(1,i,j) = 289.0_RP !+ 2.0_RP * ( rndm(KS-1,i,j)-0.50 ) * 0.1_RP ! [K]
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if ( CZ(k) < 820.0_RP ) then ! below initial cloud top
             velx(k,i,j) =   7.0_RP
             vely(k,i,j) =  -5.5_RP
             potl(k,i,j) = 289.0_RP
          else if ( CZ(k) <= 860.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 840.0_RP)/20.0_RP )
             velx(k,i,j) =   7.0_RP
             vely(k,i,j) =  -5.5_RP
             potl(k,i,j) = 289.0_RP * (1.0_RP-sint)*0.5_RP + &
                           (297.5_RP+sign(abs(CZ(k)-840.0_RP)**(1.0_RP/3.0_RP),CZ(k)-840.0_RP)) * (1.0_RP+sint)*0.5_RP
          else
             velx(k,i,j) =   7.0_RP
             vely(k,i,j) =  -5.5_RP
             potl(k,i,j) = 297.5_RP + ( CZ(k)-840.0_RP )**(1.0_RP/3.0_RP) ! [K]
          endif

          qv(k,i,j) = 0.0_RP
          qc(k,i,j) = 0.0_RP
       enddo

    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), potl    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    ! calc in moist condition
    do j = JS, JE
    do i = IS, IE

       qv_sfc  (1,i,j) = 9.0E-3_RP   ! [kg/kg]
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE

          if ( CZ(k) < 820.0_RP ) then ! below initial cloud top
             qall(k,i,j) = 9.0E-3_RP ! [kg/kg]
          elseif ( CZ(k) <= 860.0_RP ) then ! boundary
             sint = sin( pi2 * (CZ(k) - 840.0_RP)/20.0_RP )
             qall(k,i,j) = 9.0E-3_RP * (1.0_RP-sint)*0.5_RP + 1.5E-3_RP * (1.0_RP+sint)*0.5_RP
          elseif( CZ(k) <= 5000.0_RP ) then
             qall(k,i,j) = 1.5E-3_RP ! [kg/kg]
          else
             qall(k,i,j) = 0.0_RP
          endif

          if ( CZ(k) <=  600.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif ( CZ(k) < 820.0_RP ) then ! in the cloud
             fact = ( CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact
          elseif ( CZ(k) <= 860.0_RP ) then ! boundary
             sint = sin( pi2 * (CZ(k) - 840.0_RP)/20.0_RP )
             fact = ( CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact * (1.0_RP-sint)*0.5_RP ! + 0.0_RP * (1.0_RP+sint)*0.5_RP
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall(k,i,j) - qc(k,i,j)
       enddo

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LH0 / CPdry * qc(k,i,j)
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho_temp( DENS(:,:,:), pott    (:,:,:), pres    (:,:,:), temp    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                           pott_sfc(:,:,:), pres_sfc(:,:,:), temp_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

 
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
       MOMZ(k,i,j) = ( 2.0_RP * ( rndm(k,i,j)-0.50 ) * PERTURB_AMP ) &
                     * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
     else
       MOMZ(k,i,j) = 0.0_RP
     endif
     if ( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
       RHOT(k,i,j) = ( pott(k,i,j)+2.0_RP*( rndm(k,i,j)-0.5_RP )*PERTURB_AMP ) * DENS(k,i,j)
     else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
     endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
     if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
       MOMX(k,i,j) = ( velx(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50 ) * PERTURB_AMP ) &
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
     if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
       MOMY(k,i,j) = ( vely(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50 ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     else
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo


    if( QA < 30 ) then !--- for ndw6
     do j = JS, JE
     do i = IS, IE
     do k = KS, KE

        QTRC(k,i,j,I_QV) = qv(k,i,j)
        QTRC(k,i,j,I_QC) = qc(k,i,j)

        if ( qc(k,i,j) > 0.0_RP ) then
           QTRC(k,i,j,I_NC) = 120.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
        endif

     enddo
     enddo
     enddo
    else if ( QA >= 30 ) then !--- for HUCM
     do j = JS, JE
     do i = IS, IE
     do k = KS, KE

        QTRC(k,i,j,I_QV) = qv(k,i,j)+qc(k,i,j) !--- Super saturated air at initial

     enddo
     enddo
     enddo
    end if

    return
  end subroutine MKINIT_DYCOMS2_RF01

  !-----------------------------------------------------------------------------
  !> Make initial state for strato cumulus
  !-----------------------------------------------------------------------------
  subroutine MKINIT_DYCOMS2_RF02
    implicit none

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall(KA,IA,JA) ! QV+QC
    real(RP) :: qc  (KA,IA,JA) ! QC
    real(RP) :: fact !, disturb
   
    real(RP) :: pi2, sint

    integer :: ierr
    integer :: k, i, j, iq
    real(RP) :: PERTURB_AMP = 0.0_RP
    integer :: RANDOM_LIMIT = 5
    integer :: RANDOM_FLAG = 0  !0 -> no perturbation
                                !1 -> perturbation for PT  
                                !2 -> perturbation for u,v,w

    NAMELIST / PARAM_MKINIT_RF02 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DYCOMS2_RF02)]/Categ[MKINIT]'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_RF02. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_RF02)

    ! calc in dry condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1017.8E2_RP   ! [Pa]
       pott_sfc(1,i,j) = 288.3_RP      ! [K]
       qv_sfc  (1,i,j) = 9.45E-3_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          velx(k,i,j) =  3.0_RP + 4.3 * CZ(k)*1.E-3_RP
          vely(k,i,j) = -9.0_RP + 5.6 * CZ(k)*1.E-3_RP

          if ( CZ(k) < 775.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
             qall(k,i,j) = 9.45E-3_RP ! [kg/kg]
          else if ( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
             potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP + &
                   ( 295.0_RP+sign(abs(CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),CZ(k)-795.0_RP) ) * (1.0_RP+sint)*0.5_RP 
             qall(k,i,j) = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
                   ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
             qall(k,i,j) = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ! [kg/kg]
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
          qv(k,i,j) = qall(k,i,j) - qc(k,i,j)

       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), potl    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

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
     if( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then
       RHOT(k,i,j) = ( pott(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.50_RP ) * PERTURB_AMP ) &
                   * DENS(k,i,j)
     else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
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

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       QTRC(k,i,j,I_QC) = qc(k,i,j)

       if ( qc(k,i,j) > 0.0_RP ) then
          QTRC(k,i,j,I_NC) = 55.0E6_RP / DENS(k,i,j)
       endif

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_DYCOMS2_RF02
  !-----------------------------------------------------------------------------
  !> Make initial state for strato cumulus
  !-----------------------------------------------------------------------------
  subroutine MKINIT_DYCOMS2_RF01_hbinw
    implicit none

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall(KA,IA,JA) ! QV+QC
    real(RP) :: qc  (KA,IA,JA) ! QC
    real(RP) :: sint, pi2

    integer :: ierr
    integer :: k, i, j, iq
    real(RP) :: gan, xasta, xaend, dxaer
    real(RP), allocatable :: xabnd( : ), xactr( : )
    real(RP), parameter :: pi = 3.141592_RP, rhoa = 2.25E+03_RP

    real(RP) :: F0_AERO      =  1.E+7_RP ! 
    real(RP) :: R0_AERO      =  1.E-7_RP !
    real(RP) :: R_MAX        =  1.E-06_RP
    real(RP) :: R_MIN        =  1.E-08_RP
    real(RP) :: A_ALPHA      =  3.0_RP
    integer :: nccn_i        =  20
    integer :: nbin_i        =  33

    real(RP) :: PERTURB_AMP = 0.0_RP
    integer :: RANDOM_LIMIT = 5
    integer :: RANDOM_FLAG = 0  !0 -> no perturbation
                                !1 -> perturbation for PT  
                                !2 -> perturbation for u,v,w
 
   NAMELIST / PARAM_MKINIT_HBINW / &
       F0_AERO,      &
       R0_AERO,      &
       R_MAX,        &
       R_MIN,        &
       A_ALPHA,      &
       nccn_i,       &
       nbin_i

    NAMELIST / PARAM_MKINIT_RF01 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    !---------------------------------------------------------------------------
    pi2 = atan(1.0_RP) * 2.0_RP ! pi/2

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DYCOMS2_RF01_hbinw)]/Categ[MKINIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_HBINW,iostat=ierr)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF01,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_HBINW. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_HBINW)
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_RF01)

    allocate( xabnd( nccn_i+1 ) )
    allocate( xactr( nccn_i ) )

    ! calc in moist condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1017.8E2_RP ! [Pa]
       pott_sfc(1,i,j) = 289.0_RP  ! [K]
       qv_sfc  (1,i,j) = 9.0E-3_RP ! [kg/kg]
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if ( CZ(k) < 820.0_RP ) then ! below initial cloud top
             velx(k,i,j) =   7.0_RP
             vely(k,i,j) =  -5.5_RP
             potl(k,i,j) = 289.0_RP
          else if ( CZ(k) <= 860.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 840.0_RP)/20.0_RP )
             velx(k,i,j) =   7.0_RP
             vely(k,i,j) =  -5.5_RP
             potl(k,i,j) = 289.0_RP * (1.0_RP-sint)*0.5_RP + &
                           (297.5_RP+sign(abs(CZ(k)-840.0_RP)**(1.0_RP/3.0_RP),CZ(k)-840.0_RP)) * (1.0_RP+sint)*0.5_RP
          else
             velx(k,i,j) =   7.0_RP
             vely(k,i,j) =  -5.5_RP
             potl(k,i,j) = 297.5_RP + ( CZ(k)-840.0_RP )**(1.0_RP/3.0_RP) ! [K]
          endif

          if ( CZ(k) < 820.0_RP ) then ! below initial cloud top
             qall(k,i,j) = 9.0E-3_RP ! [kg/kg]
          elseif ( CZ(k) <= 860.0_RP ) then ! boundary
             sint = sin( pi2 * (CZ(k) - 840.0_RP)/20.0_RP )
             qall(k,i,j) = 9.0E-3_RP * (1.0_RP-sint)*0.5_RP + 1.5E-3_RP * (1.0_RP+sint)*0.5_RP
          elseif( CZ(k) <= 5000.0_RP ) then
             qall(k,i,j) = 1.5E-3_RP ! [kg/kg]
          else
             qall(k,i,j) = 0.0_RP
          endif

          qc(k,i,j) = 0.0_RP
          qv(k,i,j) = qall(k,i,j) - qc(k,i,j)
       enddo

    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), potl    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )


    call RANDOM_get(rndm)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMZ(k,i,j) = ( 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                     * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
      else
       MOMZ(k,i,j) = 0.0_RP
      end if
      if( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then
       RHOT(k,i,j) = ( pott(k,i,j)+2.0_RP*( rndm(k,i,j)-0.5_RP )*PERTURB_AMP )* DENS(k,i,j)
      else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
     endif
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
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
      if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMY(k,i,j) = ( vely(k,i,j) + 2.0_RP * ( rndm(k,i,j)-0.5_RP ) * PERTURB_AMP ) &
                     * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) ) 
      else
       MOMY(k,i,j) =  vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo


     do j = JS, JE
     do i = IS, IE
     do k = KS, KE

      QTRC(k,i,j,I_QV) = qv(k,i,j)
      do iq = 2, QA-nccn_i
        QTRC(k,i,j,iq) = 0.0_RP
      enddo
      xasta = log( rhoa*4.0_RP/3.0_RP*pi * ( R_MIN )**3 )
      xaend = log( rhoa*4.0_RP/3.0_RP*pi * ( R_MAX )**3 )
      dxaer = ( xaend-xasta )/nccn_i
      do iq = 1, nccn_i+1
        xabnd( iq ) = xasta + dxaer*( iq-1 )
      end do
      do iq = 1, nccn_i
       xactr( iq ) = ( xabnd( iq )+xabnd( iq+1 ) )*0.5_RP
      end do
      do iq = QA-nccn_i+1, QA
       gan = faero( F0_AERO,R0_AERO,xactr( iq-nbin_i-1 ), A_ALPHA )
       QTRC( k,i,j,iq ) = gan*exp( xactr( iq-nbin_i-1 ) )/DENS(k,i,j)
      enddo
     enddo
     enddo
     enddo
    deallocate(xactr)
    deallocate(xabnd)

    return
  end subroutine MKINIT_DYCOMS2_RF01_hbinw
  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_warmbubble_hbinw
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA              ! surface potential temperature [K]
    real(RP) :: SFC_PRES               ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  80.0_RP  ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_RH       =  80.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: SFC_U        =  0.0_RP   ! of environment [%]
    real(RP) :: ENV_L1_U     =  0.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: ENV_L3_U     =  0.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: ENV_L2_U     =  0.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_RH    =  80.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: ENV_L2_RH    =  80.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: ENV_L3_RH    =  80.0_RP  ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =   1.E3_RP  ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  =  12.E3_RP  ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L3_ZTOP  =  12.E3_RP  ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS =   4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS =   3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(RP) :: BBL_THETA    =  1.0_RP ! extremum of temperature in bubble [K]
    real(RP) :: BBL_CZ       =  2.E3_RP ! center location [m]: z
    real(RP) :: BBL_CX       =  2.E3_RP ! center location [m]: x
    real(RP) :: BBL_CY       =  2.E3_RP ! center location [m]: y
    real(RP) :: BBL_RZ       =  2.E3_RP ! bubble radius   [m]: z
    real(RP) :: BBL_RX       =  2.E3_RP ! bubble radius   [m]: x
    real(RP) :: BBL_RY       =  2.E3_RP ! bubble radius   [m]: y

    real(RP) :: F0_AERO      =  1.E+7_RP  ! number concentration of aerosol whose radii is R0_aero [#/m^3]
    real(RP) :: R0_AERO      =  1.E-7_RP  ! center radius of Jounge distribution [m]
    real(RP) :: R_MAX        =  1.E-06_RP ! max radius of aerosol [m]
    real(RP) :: R_MIN        =  1.E-08_RP ! min radius of aerosol [m]
    real(RP) :: A_ALPHA      =  3.0_RP ! index of power low distribution of aerosol (=3 -> Jounge distribution)
    integer :: nbin_i       =  33   ! number of hydrometeors bin
    integer :: nccn_i       =  20   ! number of aerosol bin

    NAMELIST / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       SFC_RH,       &
       ENV_RH,       &
       SFC_U,        &
       ENV_L1_U ,    &
       ENV_L2_U ,    &
       ENV_L3_U ,    &
       ENV_L1_RH,    &
       ENV_L2_RH,    &
       ENV_L3_RH,    &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L3_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA,    &
       BBL_CZ,       &
       BBL_CX,       &
       BBL_CY,       &
       BBL_RZ,       &
       BBL_RX,       &
       BBL_RY,       &
       F0_AERO,      &
       R0_AERO,      &
       R_MAX,        &
       R_MIN,        &
       A_ALPHA,      &
       nccn_i,       &
       nbin_i

    real(RP) :: dist, ENV_U

    integer :: ierr
    integer :: k, i, j, iq
    real(RP) :: gan, xasta, xaend, dxaer
    real(RP), allocatable :: xabnd( : ), xactr( : )
    real(RP), parameter :: pi = 3.141592_RP, rhoa = 2.25e+03_RP
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE_hbin]/Categ[MKINIT]'

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_WARMBUBBLE)

    allocate( xabnd( nccn_i+1 ) )
    allocate( xactr( nccn_i ) )

    ! calc in dry condition
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.0_RP
       qc_sfc  (1,i,j) = 0.0_RP

       do k = KS, KE
          if( k == KS ) then
             pott(k,i,j) = SFC_THETA
          elseif ( CZ(k) < ENV_L2_ZTOP ) then    ! Layer 1
             pott(k,i,j) = SFC_THETA + ENV_L2_TLAPS * ( CZ(k) )
          else 
             pott(k,i,j) =  309.02899_RP + ENV_L3_TLAPS * ( CZ(k)-ENV_L2_ZTOP )
          endif

         ! make warm bubble
          dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
            + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
            + ( (CY(j)-BBL_CY)/BBL_RY )**2
          if( dist < 1.0_RP ) then
           if( dist < ( 500.0_RP/600.0_RP )**2 ) then
            pott(k,i,j) = pott(k,i,j) + BBL_THETA 
           else if ( dist >= ( 500.0_RP/600.0_RP )**2 .and. dist < 1.0_RP ) then
            pott(k,i,j) = pott(k,i,j) + BBL_THETA-BBL_THETA/(1.d0-(500.0_RP/600.0_RP)**2)*(dist-(500.0_RP/600.0_RP)**2)
           end if 
          endif
          qv  (k,i,j) = 0.0_RP
          qc  (k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    ! calc QV from RH
    call SATURATION_pres2qsat_liq( qsat_sfc(1,:,:), temp_sfc(1,:,:), pres_sfc(1,:,:) )
    call SATURATION_pres2qsat_liq( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = SFC_RH * 1.E-2_RP * qsat_sfc(1,i,j)

       do k = KS, KE
       ! make warm bubble
        dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
             + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
             + ( (CY(j)-BBL_CY)/BBL_RY )**2
         if( dist < 1.0_RP ) then
             ENV_RH = ENV_L3_RH 
             qv(k,i,j) = ENV_RH * 1.E-2_RP * qsat(k,i,j)
         else                                ! Layer 3
          if ( CZ(k) <= ENV_L2_ZTOP ) then    ! Layer 1
             ENV_RH = SFC_RH + CZ(k)*( ENV_L2_RH-SFC_RH )/( ENV_L2_ZTOP )
             qv(k,i,j) = ENV_RH * 1.E-2_RP * qsat(k,i,j)
          else                                ! Layer 3
             ENV_RH = ENV_L3_RH 
             qv(k,i,j) = ENV_RH * 1.E-2_RP * qsat(k,i,j)
          endif
         endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), qc    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:), qc_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       MOMZ(k,i,j) = 0.0_RP
!       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
!---- SUzuki Ex.
       if( CZ(k) <= ENV_L3_ZTOP ) then
        ENV_U = SFC_U + CZ(k)*( ENV_L3_U-ENV_L1_U )/( ENV_L3_ZTOP-ENV_L1_ZTOP )
       else
        ENV_U = ENV_L3_U
       end if
       MOMX(k,i,j) = ENV_U * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )

       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA-nccn_i+1
         QTRC(k,i,j,iq) = 0.0_RP
       enddo
       xasta = log( rhoa*4.0_RP/3.0_RP*pi * ( R_MIN )**3 )
       xaend = log( rhoa*4.0_RP/3.0_RP*pi * ( R_MAX )**3 )
       dxaer = ( xaend-xasta )/nccn_i
       do iq = 1, nccn_i+1
        xabnd( iq ) = xasta + dxaer*( iq-1 )
       end do
       do iq = 1, nccn_i
        xactr( iq ) = ( xabnd( iq )+xabnd( iq+1 ) )*0.5_RP
       end do
       do iq = QA-nccn_i+1, QA
        gan = faero( F0_AERO,R0_AERO,xactr( iq-nbin_i-1 ), A_ALPHA )
        QTRC( k,i,j,iq ) = gan*exp( xactr( iq-nbin_i-1 ) )/DENS(k,i,j)
       enddo
       do iq = 2, QA-nccn_i
         QTRC(k,i,j,iq) = 0.0_RP
       enddo

       ! make warm bubble
!       dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
!            + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
!            + ( (CY(j)-BBL_CY)/BBL_RY )**2

!       if ( dist <= 1.0_RP ) then
!          RHOT(k,i,j) = RHOT(k,i,j) + DENS(k,i,j) * BBL_THETA * cos( 0.5_RP*PI*sqrt(dist) )**2
!             ENV_RH = ENV_L3_RH 
!             qv(k,i,j) = ENV_RH * 1.E-2_RP * qsat(k,i,j)
!       endif

    enddo
    enddo
    enddo
    deallocate(xactr)
    deallocate(xabnd)

    return
  end subroutine MKINIT_warmbubble_hbinw

  function faero( f0,r0,x,alpha )

  real(RP), intent(in) ::  x, f0, r0, alpha
  real(RP) :: faero
  real(RP) :: rad
  real(RP), parameter :: pi = 3.141592_RP, rhoa = 2.25E+03_RP

  rad = ( exp( x )*3.0_RP/4.0_RP/pi/rhoa )**( 1.0_RP/3.0_RP )

  faero = f0*( rad/r0 )**( -alpha )

  return

  end function faero

end module mod_mkinit
