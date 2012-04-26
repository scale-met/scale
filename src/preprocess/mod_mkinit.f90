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
     hydro_buildrho      => ATMOS_HYDRO_buildrho
  use mod_atmos_saturation, only: &
     saturation_qsat_sfc   => ATMOS_SATURATION_qsat_sfc, &
     saturation_qsat_water => ATMOS_SATURATION_qsat_water
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
  public :: MKINIT_squalline
  public :: MKINIT_DYCOMS2_RF01
  public :: MKINIT_DYCOMS2_RF02

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, save      :: MKINIT_TYPE
  integer, public, parameter :: I_PLANESTATE    =  1
  integer, public, parameter :: I_TRACERBUBBLE  =  2
  integer, public, parameter :: I_COLDBUBBLE    =  3
  integer, public, parameter :: I_WARMBUBBLE    =  4
  integer, public, parameter :: I_KHWAVE        =  5
  integer, public, parameter :: I_TURBULENCE    =  6
  integer, public, parameter :: I_SUPERCELL     =  7
  integer, public, parameter :: I_SQUALLINE     =  8
  integer, public, parameter :: I_DYCOMS2_RF01  =  9
  integer, public, parameter :: I_DYCOMS2_RF02  = 10
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(8), private, parameter :: THETAstd = 300.D0 ! [K]

  real(8), private :: pres(KA,IA,JA) ! pressure [Pa]
  real(8), private :: temp(KA,IA,JA) ! temperature [K]
  real(8), private :: pott(KA,IA,JA) ! potential temperature [K]
  real(8), private :: qsat(KA,IA,JA) ! satulated water vapor [kg/kg]
  real(8), private :: qv  (KA,IA,JA) ! water vapor [kg/kg]
  real(8), private :: velx(KA,IA,JA) ! velocity u [m/s]
  real(8), private :: vely(KA,IA,JA) ! velocity v [m/s]

  real(8), private :: pres_sfc(1,IA,JA)
  real(8), private :: temp_sfc(1,IA,JA)
  real(8), private :: pott_sfc(1,IA,JA)
  real(8), private :: qsat_sfc(1,IA,JA)
  real(8), private :: qv_sfc  (1,IA,JA)

  real(8), private :: rndm(KA,IA,JA) ! random number (0-1)
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine MKINIT_setup
    implicit none

    character(len=IO_SYSCHR) :: MKINIT_initname = 'COLDBUBBLE'

    NAMELIST / PARAM_MKINIT / &
       MKINIT_initname

    integer :: ierr
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

    select case(trim(MKINIT_initname))
    case('PLANESTATE')
       MKINIT_TYPE = I_PLANESTATE
    case('TRACERBUBBLE')
       MKINIT_TYPE = I_TRACERBUBBLE
    case('COLDBUBBLE')
       MKINIT_TYPE = I_COLDBUBBLE
    case('WARMBUBBLE')
       MKINIT_TYPE = I_WARMBUBBLE
    case('KHWAVE')
       MKINIT_TYPE = I_KHWAVE
    case('TURBULENCE')
       MKINIT_TYPE = I_TURBULENCE
    case('SUPERCELL')
       MKINIT_TYPE = I_SUPERCELL
    case('SQUALLINE')
       MKINIT_TYPE = I_SQUALLINE
    case('DYCOMS2_RF01')
       MKINIT_TYPE = I_DYCOMS2_RF01
    case('DYCOMS2_RF02')
       MKINIT_TYPE = I_DYCOMS2_RF02
    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(MKINIT_initname)
       call PRC_MPIstop
    endselect

    return
  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform + random disturbance )
  !-----------------------------------------------------------------------------
  subroutine MKINIT_planestate
    implicit none

    ! Surface state
    real(8) :: SFC_THETA            ! surface potential temperature [K]
    real(8) :: SFC_PRES             ! surface pressure [Pa]
    real(8) :: SFC_RH       =  0.D0 ! surface relative humidity [%]
    ! Environment state
    real(8) :: ENV_THETA            ! potential temperature of environment [K]
    real(8) :: ENV_W        =  0.D0 ! velocity w of environment [m/s]
    real(8) :: ENV_U        =  0.D0 ! velocity u of environment [m/s]
    real(8) :: ENV_V        =  0.D0 ! velocity v of environment [m/s]
    real(8) :: ENV_RH       =  0.D0 ! relative humidity of environment [%]
    ! Disturbance
    real(8) :: RANDOM_THETA =  0.D0 ! amplitude of random disturbance theta
    real(8) :: RANDOM_W     =  0.D0 ! amplitude of random disturbance w
    real(8) :: RANDOM_U     =  0.D0 ! amplitude of random disturbance u
    real(8) :: RANDOM_V     =  0.D0 ! amplitude of random disturbance v
    real(8) :: RANDOM_RH    =  0.D0 ! amplitude of random disturbance RH
    real(8) :: RANDOM_QTRC  =  0.D0 ! amplitude of random disturbance QTRC

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
       RANDOM_RH,    &
       RANDOM_QTRC

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
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + rndm(k,i,j) * RANDOM_THETA
          qv  (k,i,j) = 0.D0
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    ! calc QV from RH
    call saturation_qsat_sfc  ( qsat_sfc(:,:,:), temp_sfc(:,:,:), pres_sfc(:,:,:) )
    call saturation_qsat_water( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.D-2 * qsat_sfc(1,i,j)
       qv_sfc(1,i,j) = max( qv_sfc(1,i,j), 0.D0 )

       do k = KS, KE
          qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.D-2 * qsat(k,i,j)
          qv(k,i,j) = max( qv(k,i,j), 0.D0 )
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = ( ENV_W + ( rndm(k,i,j) - 0.50 ) * 2.D0 * RANDOM_W ) &
                   * 0.5D0 * ( DENS(k,i,j) + DENS(k+1,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) - 0.50 ) * 2.D0 * RANDOM_U ) &
                   * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) - 0.50 ) * 2.D0 * RANDOM_V ) &
                   * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
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

    if ( QA >= 2 ) then
       do iq = 2, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = rndm(k,i,j) * RANDOM_QTRC
       enddo
       enddo
       enddo
       enddo
    endif

    return
  end subroutine MKINIT_planestate

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_tracerbubble
    implicit none

    ! Surface state
    real(8) :: SFC_THETA            ! surface potential temperature [K]
    real(8) :: SFC_PRES             ! surface pressure [Pa]
    ! Environment state
    real(8) :: ENV_THETA            ! potential temperature of environment [K]
    real(8) :: ENV_U        =  0.D0 ! velocity u of environment [m/s]
    real(8) :: ENV_V        =  0.D0 ! velocity v of environment [m/s]
    ! Bubble
    real(8) :: BBL_NC       =  1.D0 ! extremum of NC in bubble [kg/kg]
    real(8) :: BBL_CZ       =  2.D3 ! center location [m]: z
    real(8) :: BBL_CX       =  2.D3 ! center location [m]: x
    real(8) :: BBL_CY       =  2.D3 ! center location [m]: y
    real(8) :: BBL_RZ       =  2.D3 ! bubble radius   [m]: z
    real(8) :: BBL_RX       =  2.D3 ! bubble radius   [m]: x
    real(8) :: BBL_RY       =  2.D3 ! bubble radius   [m]: y

    NAMELIST / PARAM_MKINIT_TRACERBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       ENV_U,     &
       ENV_V,     &
       BBL_NC,    &
       BBL_CZ,    &
       BBL_CX,    &
       BBL_CY,    &
       BBL_RZ,    &
       BBL_RX,    &
       BBL_RY

    real(8) :: dist

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
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          pott(k,i,j) = ENV_THETA
          qv  (k,i,j) = 0.D0
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = ENV_U * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
       MOMY(k,i,j) = ENV_V * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo

       ! make tracer bubble
       dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
            + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
            + ( (CY(j)-BBL_CY)/BBL_RY )**2

       if ( dist <= 1.D0 ) then
          QTRC(k,i,j,I_NC) = QTRC(k,i,j,I_NC) + BBL_NC * cos( 0.5D0*PI*sqrt(dist) )**2
       endif

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_tracerbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_coldbubble
    implicit none

    ! Surface state
    real(8) :: SFC_THETA            ! surface potential temperature [K]
    real(8) :: SFC_PRES             ! surface pressure [Pa]
    ! Environment state
    real(8) :: ENV_THETA            ! potential temperature of environment [K]
    ! Bubble
    real(8) :: BBL_THETA    = -5.D0 ! extremum of temperature in bubble [K]
    real(8) :: BBL_CZ       =  2.D3 ! center location [m]: z
    real(8) :: BBL_CX       =  2.D3 ! center location [m]: x
    real(8) :: BBL_CY       =  2.D3 ! center location [m]: y
    real(8) :: BBL_RZ       =  2.D3 ! bubble radius   [m]: z
    real(8) :: BBL_RX       =  2.D3 ! bubble radius   [m]: x
    real(8) :: BBL_RY       =  2.D3 ! bubble radius   [m]: y

    NAMELIST / PARAM_MKINIT_COLDBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       BBL_THETA, &
       BBL_CZ,    &
       BBL_CX,    &
       BBL_CY,    &
       BBL_RZ,    &
       BBL_RX,    &
       BBL_RY

    real(8) :: dist

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
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          pott(k,i,j) = ENV_THETA
          qv  (k,i,j) = 0.D0
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = 0.D0
       MOMY(k,i,j) = 0.D0
       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo

       ! make cold bubble
       dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
            + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
            + ( (CY(j)-BBL_CY)/BBL_RY )**2

       if ( dist <= 1.D0 ) then
          RHOT(k,i,j) = RHOT(k,i,j) + DENS(k,i,j) * BBL_THETA * cos( 0.5D0*PI*sqrt(dist) )**2
       endif

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
    real(8) :: SFC_THETA              ! surface potential temperature [K]
    real(8) :: SFC_PRES               ! surface pressure [Pa]
    real(8) :: SFC_RH       =  80.D0  ! surface relative humidity [%]
    ! Environment state
    real(8) :: ENV_RH       =  80.D0  ! Relative Humidity of environment [%]
    real(8) :: ENV_L1_ZTOP  =   1.D3  ! top height of the layer1 (constant THETA)       [m]
    real(8) :: ENV_L2_ZTOP  =  12.D3  ! top height of the layer2 (small THETA gradient) [m]
    real(8) :: ENV_L2_TLAPS =   4.D-3 ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(8) :: ENV_L3_TLAPS =   3.D-2 ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(8) :: BBL_THETA    =  1.D0 ! extremum of temperature in bubble [K]
    real(8) :: BBL_CZ       =  2.D3 ! center location [m]: z
    real(8) :: BBL_CX       =  2.D3 ! center location [m]: x
    real(8) :: BBL_CY       =  2.D3 ! center location [m]: y
    real(8) :: BBL_RZ       =  2.D3 ! bubble radius   [m]: z
    real(8) :: BBL_RX       =  2.D3 ! bubble radius   [m]: x
    real(8) :: BBL_RY       =  2.D3 ! bubble radius   [m]: y

    NAMELIST / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       ENV_RH,       &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA,    &
       BBL_CZ,       &
       BBL_CX,       &
       BBL_CY,       &
       BBL_RZ,       &
       BBL_RX,       &
       BBL_RY

    real(8) :: dist

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
    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          if ( CZ(k) <= ENV_L1_ZTOP ) then    ! Layer 1
             pott(k,i,j) = SFC_THETA
          elseif( CZ(k) <= ENV_L2_ZTOP ) then ! Layer 2
             pott(k,i,j) = pott(k-1,i,j) + ENV_L2_TLAPS * ( CZ(k)-CZ(k-1) )
          else                                ! Layer 3
             pott(k,i,j) = pott(k-1,i,j) + ENV_L3_TLAPS * ( CZ(k)-CZ(k-1) )
          endif
          qv  (k,i,j) = 0.D0
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    ! calc QV from RH
    call saturation_qsat_sfc  ( qsat_sfc(:,:,:), temp_sfc(:,:,:), pres_sfc(:,:,:) )
    call saturation_qsat_water( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = SFC_RH * 1.D-2 * qsat_sfc(1,i,j)

       do k = KS, KE
          if ( CZ(k) <= ENV_L1_ZTOP ) then    ! Layer 1
             qv(k,i,j) = ENV_RH * 1.D-2 * qsat(k,i,j)
          elseif( CZ(k) <= ENV_L2_ZTOP ) then ! Layer 2
             qv(k,i,j) = ENV_RH * 1.D-2 * qsat(k,i,j)
          else                                ! Layer 3
             qv(k,i,j) = 0.D0
          endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = 0.D0
       MOMY(k,i,j) = 0.D0
       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo

       ! make warm bubble
       dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
            + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
            + ( (CY(j)-BBL_CY)/BBL_RY )**2

       if ( dist <= 1.D0 ) then
          RHOT(k,i,j) = RHOT(k,i,j) + DENS(k,i,j) * BBL_THETA * cos( 0.5D0*PI*sqrt(dist) )**2
       endif

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
    real(8) :: SFC_THETA               ! surface potential temperature [K]
    real(8) :: SFC_PRES                ! surface pressure [Pa]
    real(8) :: SFC_RH         =  0.D0  ! surface relative humidity [%]
    ! Environment state
    real(8) :: ENV_L1_ZTOP    = 1.95D3 ! top    height of the layer1 (low  THETA) [m]
    real(8) :: ENV_L3_ZBOTTOM = 2.05D3 ! bottom height of the layer3 (high THETA) [m]
    real(8) :: ENV_L1_THETA   = 300.D0 ! THETA in the layer1 (low  THETA) [K]
    real(8) :: ENV_L3_THETA   = 305.D0 ! THETA in the layer3 (high THETA) [K]
    real(8) :: ENV_L1_U       =   0.D0 ! velocity u in the layer1 (low  THETA) [K]
    real(8) :: ENV_L3_U       =  20.D0 ! velocity u in the layer3 (high THETA) [K]
    real(8) :: ENV_L1_RH      =  50.D0 ! Relative Humidity in the layer1 (low  THETA) [%]
    real(8) :: ENV_L3_RH      =   0.D0 ! Relative Humidity in the layer3 (high THETA) [%]

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

    real(8) :: fact
    
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
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          if ( CZ(k) <= ENV_L1_ZTOP ) then       ! Layer 1
             pott(k,i,j) = ENV_L1_THETA
          elseif( CZ(k) >= ENV_L3_ZBOTTOM ) then ! Layer 3
             pott(k,i,j) = ENV_L3_THETA
          else                                   ! Layer 2
             fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

             pott(k,i,j) = ENV_L1_THETA * ( 1.D0 - fact ) &
                         + ENV_L3_THETA * (        fact )
          endif
          qv  (k,i,j) = 0.D0
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    ! calc QV from RH
    call saturation_qsat_sfc  ( qsat_sfc(:,:,:), temp_sfc(:,:,:), pres_sfc(:,:,:) )
    call saturation_qsat_water( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = SFC_RH * 1.D-2 * qsat_sfc(1,i,j)

       do k = KS, KE
          if ( CZ(k) <= ENV_L1_ZTOP ) then    ! Layer 1
             qv(k,i,j) = ENV_L1_RH * 1.D-2 * qsat_sfc(k,i,j)
          elseif( CZ(k) <= ENV_L3_ZBOTTOM ) then ! Layer 3
             qv(k,i,j) = ENV_L3_RH * 1.D-2 * qsat_sfc(k,i,j)
          else                                ! Layer 2
             fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

             qv(k,i,j) = ( ENV_L1_RH * ( 1.D0 - fact ) &
                         + ENV_L3_RH * (        fact ) ) * 1.D-2 * qsat_sfc(k,i,j)
          endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       MOMZ(k,i,j) = 0.D0
       MOMY(k,i,j) = 0.D0
       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo

       fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

       MOMX(k,i,j) = ( ENV_L1_U * ( 1.D0 - fact ) &
                     + ENV_L3_U * (        fact ) ) &
                   * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )

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
    real(8) :: SFC_THETA             ! surface potential temperature [K]
    real(8) :: SFC_PRES              ! surface pressure [Pa]
    real(8) :: SFC_RH       = 50.D0  ! surface relative humidity [%]
    ! Environment state
    real(8) :: ENV_THETA             ! potential temperature of environment [K]
    real(8) :: ENV_TLAPS    =  4.D-3 ! Lapse rate of THETA [K/m]
    real(8) :: ENV_U        =  5.D0  ! velocity u of environment [m/s]
    real(8) :: ENV_V        =  0.D0  ! velocity v of environment [m/s]
    real(8) :: ENV_RH       = 50.D0  ! relative humidity of environment [%]
    ! Disturbance
    real(8) :: RANDOM_THETA =  3.D0  ! amplitude of random disturbance theta
    real(8) :: RANDOM_U     =  5.D-2 ! amplitude of random disturbance u
    real(8) :: RANDOM_V     =  0.D0  ! amplitude of random disturbance v
    real(8) :: RANDOM_RH    =  5.D-1 ! amplitude of random disturbance RH

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
       pott_sfc(1,i,j) = SFC_THETA + rndm(KS-1,i,j) * RANDOM_THETA
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * CZ(k) + rndm(k,i,j) * RANDOM_THETA
          qv  (k,i,j) = 0.D0
       enddo
    enddo
    enddo

    ! make density & pressure profile in dry condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    ! calc QV from RH
    call saturation_qsat_sfc  ( qsat_sfc(:,:,:), temp_sfc(:,:,:), pres_sfc(:,:,:) )
    call saturation_qsat_water( qsat    (:,:,:), temp    (:,:,:), pres    (:,:,:) )

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       qv_sfc(1,i,j) = ( SFC_RH + rndm(KS-1,i,j) * RANDOM_RH ) * 1.D-2 * qsat_sfc(1,i,j)

       do k = KS, KE
          qv(k,i,j) = ( ENV_RH + rndm(k,i,j) * RANDOM_RH ) * 1.D-2 * qsat(k,i,j)
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.D0
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) - 0.50 ) * 2.D0 * RANDOM_U ) &
                   * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) - 0.50 ) * 2.D0 * RANDOM_V ) &
                   * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    if ( QA >= 2 ) then
       do iq = 2, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = 0.D0
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
    real(8) :: BBL_THETA    = -5.D0 ! extremum of temperature in bubble [K]
    real(8) :: BBL_CZ       =  2.D3 ! center location [m]: z
    real(8) :: BBL_CX       =  2.D3 ! center location [m]: x
    real(8) :: BBL_CY       =  2.D3 ! center location [m]: y
    real(8) :: BBL_RZ       =  2.D3 ! bubble radius   [m]: z
    real(8) :: BBL_RX       =  2.D3 ! bubble radius   [m]: x
    real(8) :: BBL_RY       =  2.D3 ! bubble radius   [m]: y

    NAMELIST / PARAM_MKINIT_SUPERCELL / &
       ENV_IN_SOUNDING_file, &
       BBL_THETA, &
       BBL_CZ,    &
       BBL_CX,    &
       BBL_CY,    &
       BBL_RZ,    &
       BBL_RX,    &
       BBL_RY

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(8) :: SFC_THETA          ! surface potential temperature [K]
    real(8) :: SFC_PRES           ! surface pressure [Pa]
    real(8) :: SFC_QV             ! surface watervapor [g/kg]

    real(8) :: EXP_z   (EXP_klim) ! height      [m]
    real(8) :: EXP_pott(EXP_klim) ! potential temperature [K]
    real(8) :: EXP_qv  (EXP_klim) ! water vapor [g/kg]
    real(8) :: EXP_u   (EXP_klim) ! velocity u  [m/s]
    real(8) :: EXP_v   (EXP_klim) ! velocity v  [m/s]

    real(8) :: dist
    real(8) :: fact1, fact2

    integer :: ierr, fid
    integer :: k, i, j, iq, kref
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[INIT]'

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

    EXP_z   (1) = 0.D0
    EXP_pott(1) = SFC_THETA
    EXP_qv  (1) = SFC_QV
    EXP_u   (1) = EXP_u(2)
    EXP_v   (1) = EXP_v(2)
    do k = 1, EXP_klim
       EXP_qv(k) = EXP_qv(k) * 1.D-3 ![g/kg] -> [kg/kg]
    enddo

    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = SFC_QV
    enddo
    enddo

    !--- linear interpolate to model grid
    do k    = KS, KE
    do kref = 2, EXP_kmax

       if (       CZ(k) >  EXP_z(kref-1) &
            .AND. CZ(k) <= EXP_z(kref)   ) then

          fact1 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )
          fact2 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )

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

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = velx(k,i,j) * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
       MOMY(k,i,j) = vely(k,i,j) * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo

       ! make warm bubble
       dist = ( (CZ(k)-BBL_CZ)/BBL_RZ )**2 &
            + ( (CX(i)-BBL_CX)/BBL_RX )**2 &
            + ( (CY(j)-BBL_CY)/BBL_RY )**2

       if ( dist <= 1.D0 ) then
          RHOT(k,i,j) = RHOT(k,i,j) + DENS(k,i,j) * BBL_THETA * cos( 0.5D0*PI*sqrt(dist) )**2
       endif

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_supercell

  !-----------------------------------------------------------------------------
  !> Make initial state for squalline experiment
  !-----------------------------------------------------------------------------
  subroutine MKINIT_squalline
    implicit none

    character(len=IO_FILECHR) :: ENV_IN_SOUNDING_file = ''

    NAMELIST / PARAM_MKINIT_SQUALLINE / &
       ENV_IN_SOUNDING_file

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(8) :: SFC_THETA          ! surface potential temperature [K]
    real(8) :: SFC_PRES           ! surface pressure [Pa]
    real(8) :: SFC_QV             ! surface watervapor [g/kg]

    real(8) :: EXP_z   (EXP_klim) ! height      [m]
    real(8) :: EXP_pott(EXP_klim) ! potential temperature [K]
    real(8) :: EXP_qv  (EXP_klim) ! water vapor [g/kg]
    real(8) :: EXP_u   (EXP_klim) ! velocity u  [m/s]
    real(8) :: EXP_v   (EXP_klim) ! velocity v  [m/s]

    real(8) :: fact1, fact2

    integer :: ierr, fid
    integer :: k, i, j, iq, kref
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[INIT]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SQUALLINE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SQUALLINE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKINIT_SQUALLINE)

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

    EXP_z   (1) = 0.D0
    EXP_pott(1) = SFC_THETA
    EXP_qv  (1) = SFC_QV
    EXP_u   (1) = EXP_u(2)
    EXP_v   (1) = EXP_v(2)
    do k = 1, EXP_klim
       EXP_qv(k) = EXP_qv(k) * 1.D-3 ![g/kg] -> [kg/kg]
    enddo

    do j = JS, JE
    do i = IS, IE
       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA
       qv_sfc  (1,i,j) = SFC_QV
    enddo
    enddo

    !--- linear interpolate to model grid
    do k    = KS, KE
    do kref = 2, EXP_kmax

       if (       CZ(k) >  EXP_z(kref-1) &
            .AND. CZ(k) <= EXP_z(kref)   ) then

          fact1 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )
          fact2 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )

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

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = DENS(k,i,j) * velx(k,i,j)
       MOMY(k,i,j) = DENS(k,i,j) * vely(k,i,j)
       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       do iq = 2, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_squalline

  !-----------------------------------------------------------------------------
  !> Make initial state for strato cumulus
  !-----------------------------------------------------------------------------
  subroutine MKINIT_stratocumulus
    implicit none

    ! Surface state
    real(8) :: SFC_THETA               ! surface potential temperature [K]
    real(8) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(8) :: ENV_L1_ZTOP    = 840.0D0  ! top    height of the layer1 (low  THETA) [m]
    real(8) :: ENV_L3_ZBOTTOM = 840.0D0  ! bottom height of the layer3 (high THETA) [m]
    real(8) :: ENV_L1_THETA   = 289.0D0  ! THETA in the layer1 (low  THETA) [K]
    real(8) :: ENV_L3_THETA   = 297.5D0  ! THETA in the layer3 (high THETA) [K]
    real(8) :: ENV_L1_QV      =   8.5D-3 ! Specific Humidity in the layer1 (low  THETA) [kg/kg]
    real(8) :: ENV_L3_QV      =   1.0D-3 ! Specific Humidity in the layer3 (high THETA) [kg/kg]

    real(8) :: ENV_CL_ZBOTTOM = 600.0D0  ! bottom height of the cloud layer [m]
    real(8) :: ENV_CL_ZTOP    = 840.0D0  ! top    height of the cloud layer [m]
    real(8) :: ENV_CL_QC      =   0.5D-3 ! cloud water mixing ratio in the cloud layer   [kg/kg]
    real(8) :: ENV_CL_NC      = 120.0D6  ! cloud number concentration in the cloud layer [1/m3]

    real(8) :: ENV_U          =   7.0D0  ! velocity u in the layer1 (low  THETA) [K]
    real(8) :: ENV_V          =  -5.5D0  ! velocity u in the layer3 (high THETA) [K]

    real(8) :: RANDOM_AMP     =   1.0D-2 ! ratio of random disturbance [0-1]

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

    real(8) :: fact, disturb

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
       disturb = ( 1.D0 + 2.D0 * ( rndm(KS-1,i,j)-0.50 ) * RANDOM_AMP )

       pres_sfc(1,i,j) = SFC_PRES
       pott_sfc(1,i,j) = SFC_THETA * disturb
       qv_sfc  (1,i,j) = ENV_L1_QV

       do k = KS, KE
          disturb = ( 1.D0 + 2.D0 * ( rndm(k,i,j)-0.50 ) * RANDOM_AMP )

          if ( CZ(k) <= ENV_L1_ZTOP ) then       ! Layer 1
             pott(k,i,j) = ENV_L1_THETA * disturb
             qv  (k,i,j) = ENV_L1_QV
          elseif( CZ(k) >= ENV_L3_ZBOTTOM ) then ! Layer 3
             pott(k,i,j) = ENV_L3_THETA * disturb &
                         + ( CZ(k) - ENV_L3_ZBOTTOM )**(1.D0/3.D0) ! increase with height
             qv  (k,i,j) = ENV_L3_QV
          else                                   ! Layer 2
             fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )

             pott(k,i,j) = ( ENV_L1_THETA * ( 1.D0 - fact ) &
                           + ENV_L3_THETA * (        fact ) ) * disturb
             qv  (k,i,j) = ENV_L1_QV    * ( 1.D0 - fact ) &
                         + ENV_L3_QV    * (        fact )
          endif
       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.D0
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U * ( 1.D0 + 2.D0 * ( rndm(k,i,j)-0.50 ) * RANDOM_AMP ) ) &
                   * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V * ( 1.D0 + 2.D0 * ( rndm(k,i,j)-0.50 ) * RANDOM_AMP ) ) &
                   * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.D0
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

    real(8) :: potl(KA,IA,JA) ! liquid potential temperature
    real(8) :: qall(KA,IA,JA) ! QV+QC
    real(8) :: qc  (KA,IA,JA) ! QC
    real(8) :: fact

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DYCOMS2_RF01)]/Categ[MKINIT]'

    ! calc in moist condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE

       pres_sfc(1,i,j) = 1017.8D2 ! [Pa]
       pott_sfc(1,i,j) = 292.5D0 + 2.D0 * ( rndm(KS-1,i,j)-0.50 ) * 0.1D0 ! [K]
       qv_sfc  (1,i,j) = 9.0D-3   ! [kg/kg]

       do k = KS, KE
          if ( CZ(k) <= 840.D0 ) then ! below initial cloud top
!             velx(k,i,j) =   6.7D0
!             vely(k,i,j) =  -4.9D0
             velx(k,i,j) =   7.0D0
             vely(k,i,j) =  -5.5D0
             potl(k,i,j) = 289.0D0 + 2.D0 * ( rndm(k,i,j)-0.50 ) * 0.1D0 ! [K]
          else
             velx(k,i,j) =   7.0D0
             vely(k,i,j) =  -5.5D0
             potl(k,i,j) = 297.5D0 + ( CZ(k)-840.D0 )**(1.D0/3.D0) ! [K]
          endif

          if ( CZ(k) <= 840.D0 ) then ! below initial cloud top
             qall(k,i,j) = 9.0D-3 ! [kg/kg]
          elseif(       CZ(k) >   840.D0 &
                  .AND. CZ(k) <= 5000.D0 ) then
             qall(k,i,j) = 1.5D-3 ! [kg/kg]
          else
             qall(k,i,j) = 0.0D0
          endif

          if (       CZ(k) >  600.D0 &
               .AND. CZ(k) <= 840.D0 ) then ! in the cloud
             fact = ( CZ(k)-600.D0 ) / ( 840.D0-600.D0 )

             qc(k,i,j) = 0.45D-3 * fact
          else
             qc(k,i,j) = 0.D0
          endif
          qv(k,i,j) = qall(k,i,j) - qc(k,i,j)
       enddo

    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), potl    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )


    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.D0
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( velx(k,i,j) + 2.D0 * ( rndm(k,i,j)-0.50 ) * 0.1D0 ) &
                   * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( vely(k,i,j) + 2.D0 * ( rndm(k,i,j)-0.50 ) * 0.1D0 ) &
                   * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.D0
    enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       QTRC(k,i,j,I_QC) = qc(k,i,j)

       if ( qc(k,i,j) > 0.D0 ) then
          QTRC(k,i,j,I_NC) = 120.D6 / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_DYCOMS2_RF01

  !-----------------------------------------------------------------------------
  !> Make initial state for strato cumulus
  !-----------------------------------------------------------------------------
  subroutine MKINIT_DYCOMS2_RF02
    implicit none

    real(8) :: potl(KA,IA,JA) ! liquid potential temperature
    real(8) :: qall(KA,IA,JA) ! QV+QC
    real(8) :: qc  (KA,IA,JA) ! QC
    real(8) :: fact, disturb

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[DYCOMS2_RF01)]/Categ[MKINIT]'

    ! calc in dry condition
    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
       disturb = 2.D0 * ( rndm(KS-1,i,j)-0.50 ) ! [-1 to 1]

       pres_sfc(1,i,j) = 1017.8D2                ! [Pa]
       pott_sfc(1,i,j) = 288.3 + disturb * 0.1D0 ! [K]
       qv_sfc  (1,i,j) = 0.D0

       do k = KS, KE
          disturb = 2.D0 * ( rndm(k,i,j)-0.50 ) ! [-1 to 1]

          velx(k,i,j) =  3.D0 + 4.3 * CZ(k)*1.D-3
          vely(k,i,j) = -9.D0 + 5.6 * CZ(k)*1.D-3

          if ( CZ(k) <= 795.D0 ) then ! below initial cloud top
             potl(k,i,j) = 288.3D0 + disturb * 0.1D0 ! [K]
             qall(k,i,j) = 9.45D-3 ! [kg/kg]
          else
             potl(k,i,j) = 295.D0 + ( CZ(k)-795.D0 )**(1.D0/3.D0) + disturb * 0.1D0 ! [K]
             qall(k,i,j) = 5.D-3 - 3.D-3 * ( 1.D0- exp( (CZ(k)-795.D0)/500.D0 ) ) ! [kg/kg]
          endif

          if (       CZ(k) >= 400.D0 &
               .AND. CZ(k) <= 795.D0 ) then ! in the cloud
             fact = ( CZ(k)-400.D0 ) / ( 795.D0-400.D0 )

             qc(k,i,j) = 0.5D-3 * fact
          else
             qc(k,i,j) = 0.D0
          endif
          qv(k,i,j) = qall(k,i,j) - qc(k,i,j)

       enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), potl    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       pott(k,i,j) = potl(k,i,j) + LH0 / CPdry * qc(k,i,j) * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! make density & pressure profile in moist condition
    call hydro_buildrho( DENS(:,:,:), temp    (:,:,:), pres    (:,:,:), pott    (:,:,:), qv    (:,:,:), &
                                      temp_sfc(:,:,:), pres_sfc(:,:,:), pott_sfc(:,:,:), qv_sfc(:,:,:)  )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ(k,i,j) = 0.D0
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = ( velx(k,i,j) + 2.D0 * ( rndm(k,i,j)-0.50 ) * 0.1D0 ) &
                   * 0.5D0 * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    call RANDOM_get(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMY(k,i,j) = ( vely(k,i,j) + 2.D0 * ( rndm(k,i,j)-0.50 ) * 0.1D0 ) &
                   * 0.5D0 * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,iq) = 0.D0
    enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       QTRC(k,i,j,I_QV) = qv(k,i,j)
       QTRC(k,i,j,I_QC) = qc(k,i,j)

       if ( qc(k,i,j) > 0.D0 ) then
          QTRC(k,i,j,I_NC) = 55.D0 / DENS(k,i,j)
       endif

    enddo
    enddo
    enddo

    return
  end subroutine MKINIT_DYCOMS2_RF02

end module mod_mkinit
