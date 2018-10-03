!-------------------------------------------------------------------------------
!> module MKTOPO
!!
!! @par Description
!!          subroutines for preparing topography data (ideal case)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mktopo
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer

  use scale_prc, only: &
     PRC_abort
  use scale_atmos_grid_cartesC, only: &
     CX => ATMOS_GRID_CARTESC_CX, &
     CY => ATMOS_GRID_CARTESC_CY
  use scale_topography, only: &
     TOPO_Zsfc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKTOPO_setup
  public :: MKTOPO

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public            :: MKTOPO_TYPE = -1
  integer, public, parameter :: I_IGNORE    =  0
  integer, public, parameter :: I_FLAT      =  1
  integer, public, parameter :: I_BELLSHAPE =  2
  integer, public, parameter :: I_SCHAER    =  3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MKTOPO_flat
  private :: MKTOPO_bellshape
  private :: MKTOPO_schaer

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKTOPO_setup
    implicit none

    character(len=H_SHORT) :: MKTOPO_name = 'NONE'

    namelist / PARAM_MKTOPO / &
       MKTOPO_name

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKTOPO_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_setup",*) 'Not appropriate names in namelist PARAM_MKTOPO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO)

    select case(MKTOPO_name)
    case('NONE')
       MKTOPO_TYPE = I_IGNORE
    case('FLAT')
       MKTOPO_TYPE = I_FLAT
    case('BELLSHAPE')
       MKTOPO_TYPE = I_BELLSHAPE
    case('SCHAER')
       MKTOPO_TYPE = I_SCHAER
    case default
       LOG_ERROR("MKTOPO_setup",*) 'Unsupported TYPE:', trim(MKTOPO_name)
       call PRC_abort
    endselect

    return
  end subroutine MKTOPO_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKTOPO
    use scale_topography, only: &
       TOPO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( MKTOPO_TYPE == I_IGNORE ) then
       LOG_NEWLINE
       LOG_INFO("MKTOPO",*) 'SKIP  MAKING TOPOGRAPHY DATA'
    else
       LOG_NEWLINE
       LOG_INFO("MKTOP",*) 'START MAKING TOPOGRAPHY DATA'

       select case(MKTOPO_TYPE)
       case(I_FLAT)
          call MKTOPO_flat

       case(I_BELLSHAPE)
          call MKTOPO_bellshape

       case(I_SCHAER)
          call MKTOPO_schaer

       case default
          LOG_ERROR("MKTOPO",*) 'Unsupported TYPE:', MKTOPO_TYPE
          call PRC_abort
       endselect

       LOG_INFO("MKTOPO",*) 'END   MAKING TOPOGRAPHY DATA'

       ! output topography file
       call TOPO_write
    endif

    return
  end subroutine MKTOPO

  !-----------------------------------------------------------------------------
  !> Make flat mountain
  subroutine MKTOPO_flat
    implicit none

    ! flat mountain parameter
    real(RP) :: FLAT_HEIGHT   =  100.0_RP ! height of mountain [m]

    namelist / PARAM_MKTOPO_FLAT / &
       FLAT_HEIGHT

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKTOPO_flat",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_FLAT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_flat",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_flat",*) 'Not appropriate names in namelist PARAM_MKTOPO_FLAT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_FLAT)

    do j = 1, JA
    do i = 1, IA
       TOPO_Zsfc(i,j) = FLAT_HEIGHT
    enddo
    enddo

    return
  end subroutine MKTOPO_flat

  !-----------------------------------------------------------------------------
  !> Make bell-shaped mountain
  subroutine MKTOPO_bellshape
    implicit none

    ! bell-shaped mountain parameter
    logical  :: BELL_eachnode = .false.   ! Arrange mountain at each node? [kg/kg]
    real(RP) :: BELL_CX       =   2.E3_RP ! center location [m]: x
    real(RP) :: BELL_CY       =   2.E3_RP ! center location [m]: y
    real(RP) :: BELL_RX       =   2.E3_RP ! bubble radius   [m]: x
    real(RP) :: BELL_RY       =   2.E3_RP ! bubble radius   [m]: y
    real(RP) :: BELL_HEIGHT   =  100.0_RP ! height of mountain [m]

    namelist / PARAM_MKTOPO_BELLSHAPE / &
       BELL_eachnode, &
       BELL_CX,       &
       BELL_CY,       &
       BELL_RX,       &
       BELL_RY,       &
       BELL_HEIGHT

    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP) :: dist

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKTOPO_bellshape",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_BELLSHAPE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_bellshape",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_bellshape",*) 'Not appropriate names in namelist PARAM_MKTOPO_BELLSHAPE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_BELLSHAPE)

    if ( BELL_eachnode ) then
       CX_offset = CX(IS)
       CY_offset = CY(JS)
    else
       CX_offset = 0.0_RP
       CY_offset = 0.0_RP
    endif

    ! make bell-shaped mountain
    do j = 1, JA
    do i = 1, IA

       dist = ( (CX(i)-CX_offset-BELL_CX)/BELL_RX )**2 &
            + ( (CY(j)-CY_offset-BELL_CY)/BELL_RY )**2

       TOPO_Zsfc(i,j) = BELL_HEIGHT / ( 1.0_RP + dist )

    enddo
    enddo

    return
  end subroutine MKTOPO_bellshape

  !-----------------------------------------------------------------------------
  !> Make Schaer-type mountain
  !> References: Schaer et al, 2002, MWR, Vol.130, 2459-2480
  !>             Klemp et al, 2003,  MWR, Vol.131, 1229-1239
  subroutine MKTOPO_schaer
    use scale_const, only: &
       PI => CONST_PI
    implicit none

    ! Schaer-type mountain parameter
    real(RP) :: SCHAER_CX       =  25.E3_RP ! center location [m]: x
    real(RP) :: SCHAER_RX       =   5.E3_RP ! bubble radius   [m]: x
    real(RP) :: SCHAER_LAMBDA   =   4.E3_RP ! wavelength of wavelike perturbation [m]: x
    real(RP) :: SCHAER_HEIGHT   =  250.0_RP ! height of mountain [m]
    logical  :: SCHAER_SWAPXY   =  .false.

    namelist / PARAM_MKTOPO_SCHAER / &
       SCHAER_CX,     &
       SCHAER_RX,     &
       SCHAER_LAMBDA, &
       SCHAER_HEIGHT, &
       SCHAER_SWAPXY

    real(RP) :: dist

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKTOPO_schaer",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKTOPO_SCHAER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKTOPO_schaer",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKTOPO_schaer",*) 'Not appropriate names in namelist PARAM_MKTOPO_SCHAER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKTOPO_SCHAER)

    ! make bell-shaped mountain
    if ( .NOT. SCHAER_SWAPXY ) then
       do j = 1, JA
       do i = 1, IA

          dist = exp( -( (CX(i)-SCHAER_CX)/SCHAER_RX )**2 )

          TOPO_Zsfc(i,j) = SCHAER_HEIGHT * dist * ( cos( PI*(CX(i)-SCHAER_CX)/SCHAER_LAMBDA ) )**2

       enddo
       enddo
    else
       do j = 1, JA
       do i = 1, IA

          dist = exp( -( (CY(j)-SCHAER_CX)/SCHAER_RX )**2 )

          TOPO_Zsfc(i,j) = SCHAER_HEIGHT * dist * ( cos( PI*(CY(j)-SCHAER_CX)/SCHAER_LAMBDA ) )**2

       enddo
       enddo
    endif

    return
  end subroutine MKTOPO_schaer

end module mod_mktopo
