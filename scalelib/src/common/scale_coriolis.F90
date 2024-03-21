!-------------------------------------------------------------------------------
!> module Coriolis
!!
!! @par Description
!!          Coriolis module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_coriolis
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CORIOLIS_setup
  public :: CORIOLIS_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: CORIOLIS_f(:,:) !> Coriolis parameter


  character(len=H_SHORT), public :: CORIOLIS_type = 'PLANE' !> type of coriolis force: 'PLANE', 'SPHERE'
                                                            !> 'PLANE' : f = CORIOLIS_f0 + CORIOLIS_beta * ( CY - CORIOLIS_y0 )
                                                            !> 'SPHERE': f = 2 * CONST_OHM * sin( lat )
  real(RP),               public :: CORIOLIS_f0   = 0.0_RP
  real(RP),               public :: CORIOLIS_beta = 0.0_RP
  real(RP),               public :: CORIOLIS_y0             !> default is domain center

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
  !> Setup
  subroutine CORIOLIS_setup( &
       IA, JA, &
       LAT,                &
       CY, DOMAIN_CENTER_Y )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       OHM => CONST_OHM
    implicit none
    integer, intent(in) :: IA, JA

    real(RP), intent(in) :: LAT(IA,JA)
    real(RP), intent(in) :: CY (JA)
    real(RP), intent(in) :: DOMAIN_CENTER_Y !< center position of global domain [m]: y

    namelist / PARAM_CORIOLIS / &
       CORIOLIS_type, &
       CORIOLIS_f0,   &
       CORIOLIS_beta, &
       CORIOLIS_y0

    integer :: i, j

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CORIOLIS_setup",*) 'Setup'

    CORIOLIS_y0 = DOMAIN_CENTER_Y

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CORIOLIS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("CORIOLIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CORIOLIS_setup",*) 'Not appropriate names in namelist PARAM_CORIOLIS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CORIOLIS)


    allocate( CORIOLIS_f(IA,JA) )

    ! coriolis parameter
    select case ( CORIOLIS_type )
    case ( 'PLANE' )
       !omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
          CORIOLIS_f(i,j) = CORIOLIS_f0 + CORIOLIS_beta * ( CY(j) - CORIOLIS_y0 )
       end do
       end do
    case ( 'SPHERE' )
       !omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
          CORIOLIS_f(i,j) = 2.0_RP * OHM * sin( LAT(i,j) )
       end do
       end do
    case default
       LOG_ERROR("CORIOLIS_setup",*) 'Coriolis type is invalid: ', trim(coriolis_type)
       LOG_ERROR_CONT(*) 'The type must be PLANE or SPHERE'
       call PRC_abort
    end select

    !$acc enter data copyin(CORIOLIS_f)

    return
  end subroutine CORIOLIS_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine CORIOLIS_finalize
    implicit none
    !---------------------------------------------------------------------------

    !$acc exit data delete(CORIOLIS_f)
    deallocate( CORIOLIS_f )

    return
  end subroutine CORIOLIS_finalize

end module scale_coriolis
