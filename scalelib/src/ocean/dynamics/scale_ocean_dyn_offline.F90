!-------------------------------------------------------------------------------
!> module ocean / dynamics / offline
!!
!! @par Description
!!          ocean offline model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_dyn_offline
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_DYN_OFFLINE_setup
  public :: OCEAN_DYN_OFFLINE

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
  !> Setup
  subroutine OCEAN_DYN_OFFLINE_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_regist
    implicit none

    character(len=H_LONG)  :: OCEAN_DYN_OFFLINE_basename              = ''
    logical                :: OCEAN_DYN_OFFLINE_basename_add_num      = .false.
    integer                :: OCEAN_DYN_OFFLINE_number_of_files       = 1
    logical                :: OCEAN_DYN_OFFLINE_enable_periodic_year  = .false.
    logical                :: OCEAN_DYN_OFFLINE_enable_periodic_month = .false.
    logical                :: OCEAN_DYN_OFFLINE_enable_periodic_day   = .false.
    integer                :: OCEAN_DYN_OFFLINE_step_fixed            = 0
    real(RP)               :: OCEAN_DYN_OFFLINE_defval               != UNDEF
    logical                :: OCEAN_DYN_OFFLINE_check_coordinates     = .true.
    integer                :: OCEAN_DYN_OFFLINE_step_limit            = 0

    namelist / PARAM_OCEAN_DYN_OFFLINE / &
       OCEAN_DYN_OFFLINE_basename,              &
       OCEAN_DYN_OFFLINE_basename_add_num,      &
       OCEAN_DYN_OFFLINE_number_of_files,       &
       OCEAN_DYN_OFFLINE_enable_periodic_year,  &
       OCEAN_DYN_OFFLINE_enable_periodic_month, &
       OCEAN_DYN_OFFLINE_enable_periodic_day,   &
       OCEAN_DYN_OFFLINE_step_fixed,            &
       OCEAN_DYN_OFFLINE_defval,                &
       OCEAN_DYN_OFFLINE_check_coordinates,     &
       OCEAN_DYN_OFFLINE_step_limit

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_DYN_OFFLINE_setup",*) 'Setup'

    OCEAN_DYN_OFFLINE_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_DYN_OFFLINE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_DYN_OFFLINE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_DYN_OFFLINE_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_DYN_OFFLINE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_DYN_OFFLINE)

    LOG_INFO("OCEAN_DYN_OFFLINE_setup",*) 'Use offline ocean'

    if ( OCEAN_DYN_OFFLINE_basename == '' ) then
       LOG_ERROR("OCEAN_DYN_OFFLINE_setup",*) 'OCEAN_DYN_OFFLINE_basename is necessary !!'
       call PRC_abort
    endif

    call FILE_EXTERNAL_INPUT_regist( OCEAN_DYN_OFFLINE_basename,              & ! [IN]
                                     OCEAN_DYN_OFFLINE_basename_add_num,      & ! [IN]
                                     OCEAN_DYN_OFFLINE_number_of_files,       & ! [IN]
                                     'OCEAN_TEMP',                            & ! [IN]
                                     'OXY',                                   & ! [IN]
                                     OCEAN_DYN_OFFLINE_enable_periodic_year,  & ! [IN]
                                     OCEAN_DYN_OFFLINE_enable_periodic_month, & ! [IN]
                                     OCEAN_DYN_OFFLINE_enable_periodic_day,   & ! [IN]
                                     OCEAN_DYN_OFFLINE_step_fixed,            & ! [IN]
                                     OCEAN_DYN_OFFLINE_defval,                & ! [IN]
                                     check_coordinates = OCEAN_DYN_OFFLINE_check_coordinates, & ! [IN]
                                     step_limit        = OCEAN_DYN_OFFLINE_step_limit         ) ! [IN]

    return
  end subroutine OCEAN_DYN_OFFLINE_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_DYN_OFFLINE( &
       OKMAX, OKS, OKE,  &
       OIA,   OIS, OIE,  &
       OJA,   OJS, OJE,  &
       calc_flag,        &
       dt, NOWDAYSEC,    &
       OCEAN_TEMP        )
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    implicit none

    integer,  intent(in)    :: OKMAX, OKS, OKE
    integer,  intent(in)    :: OIA,   OIS, OIE
    integer,  intent(in)    :: OJA,   OJS, OJE
    logical,  intent(in)    :: calc_flag (OIA,OJA) ! to decide calculate or not
    real(DP), intent(in)    :: dt
    real(DP), intent(in)    :: NOWDAYSEC
    real(RP), intent(inout) :: OCEAN_TEMP(OKMAX,OIA,OJA)

    real(RP) :: OCEAN_TEMP_ref(OKMAX,OIA,OJA)

    logical  :: error
    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'ocean / dynamics / offline'

    call FILE_EXTERNAL_INPUT_update( 'OCEAN_TEMP', NOWDAYSEC, OCEAN_TEMP_ref(:,:,:), error )

    if ( error ) then
       LOG_ERROR("OCEAN_DYN_OFFLINE",*) 'Requested data is not found!'
       call PRC_abort
    endif

    do j = OJS, OJE
    do i = OIS, OIE
       if ( calc_flag(i,j) ) then
          OCEAN_TEMP(OKS,i,j) = OCEAN_TEMP_ref(OKS,i,j)
       endif
    enddo
    enddo

    return
  end subroutine OCEAN_DYN_OFFLINE

end module scale_ocean_dyn_offline
