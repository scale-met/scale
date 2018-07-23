!-------------------------------------------------------------------------------
!> module ocean / dynamics / slab model
!!
!! @par Description
!!          ocean slab model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_dyn_slab
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
  public :: OCEAN_DYN_SLAB_setup
  public :: OCEAN_DYN_SLAB

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: OCEAN_DYN_SLAB_DEPTH = 10.0_RP !< water depth of slab ocean [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: OCEAN_DYN_SLAB_HeatCapacity             !< heat capacity of slab ocean [J/K/m2]

  logical,  private :: OCEAN_DYN_SLAB_nudging        = .false. !< SST Nudging is used?
  real(DP), private :: OCEAN_DYN_SLAB_nudging_tausec           !< Relaxation time [sec]
  logical,  private :: OCEAN_DYN_SLAB_offline_mode   = .false. !< Use offline mode?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_DYN_SLAB_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_file_limit, &
       FILE_EXTERNAL_INPUT_regist
    implicit none

    real(DP)               :: OCEAN_DYN_SLAB_nudging_tau                                      = 0.0_DP  ! Relaxation time
    character(len=H_SHORT) :: OCEAN_DYN_SLAB_nudging_tau_unit                                 = "SEC"
    character(len=H_LONG)  :: OCEAN_DYN_SLAB_nudging_basename(FILE_EXTERNAL_INPUT_file_limit) = ''
    logical                :: OCEAN_DYN_SLAB_nudging_enable_periodic_year                     = .false.
    logical                :: OCEAN_DYN_SLAB_nudging_enable_periodic_month                    = .false.
    logical                :: OCEAN_DYN_SLAB_nudging_enable_periodic_day                      = .false.
    integer                :: OCEAN_DYN_SLAB_nudging_step_fixed                               = 0
    real(RP)               :: OCEAN_DYN_SLAB_nudging_offset                                   = 0.0_RP
    real(RP)               :: OCEAN_DYN_SLAB_nudging_defval                                  != UNDEF
    logical                :: OCEAN_DYN_SLAB_nudging_check_coordinates                        = .true.
    integer                :: OCEAN_DYN_SLAB_nudging_step_limit                               = 0

    namelist / PARAM_OCEAN_DYN_SLAB / &
       OCEAN_DYN_SLAB_DEPTH,                         &
       OCEAN_DYN_SLAB_nudging,                       &
       OCEAN_DYN_SLAB_nudging_tau,                   &
       OCEAN_DYN_SLAB_nudging_tau_unit,              &
       OCEAN_DYN_SLAB_nudging_basename,              &
       OCEAN_DYN_SLAB_nudging_enable_periodic_year,  &
       OCEAN_DYN_SLAB_nudging_enable_periodic_month, &
       OCEAN_DYN_SLAB_nudging_enable_periodic_day,   &
       OCEAN_DYN_SLAB_nudging_step_fixed,            &
       OCEAN_DYN_SLAB_nudging_offset,                &
       OCEAN_DYN_SLAB_nudging_defval,                &
       OCEAN_DYN_SLAB_nudging_check_coordinates,     &
       OCEAN_DYN_SLAB_nudging_step_limit

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Setup'

    OCEAN_DYN_SLAB_nudging_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_DYN_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_DYN_SLAB_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_DYN_SLAB. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_DYN_SLAB)

    OCEAN_DYN_SLAB_HeatCapacity = DWATR * CL * OCEAN_DYN_SLAB_DEPTH

    LOG_NEWLINE
    LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Slab ocean depth [m]         : ', OCEAN_DYN_SLAB_DEPTH
    LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Ocean heat capacity [J/K/m2] : ', OCEAN_DYN_SLAB_HeatCapacity

    if ( OCEAN_DYN_SLAB_nudging ) then
       call CALENDAR_unit2sec( OCEAN_DYN_SLAB_nudging_tausec, OCEAN_DYN_SLAB_nudging_tau, OCEAN_DYN_SLAB_nudging_tau_unit )

       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Use nudging for SST : ON'
       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Relaxation time Tau [sec] : ', OCEAN_DYN_SLAB_nudging_tausec

       if ( OCEAN_DYN_SLAB_nudging_tausec == 0.0_RP ) then
          OCEAN_DYN_SLAB_offline_mode = .true.
          LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Tau=0 means that SST is completely replaced by the external data.'
       endif

       if ( OCEAN_DYN_SLAB_nudging_basename(1) == '' ) then
          LOG_ERROR("OCEAN_DYN_SLAB_setup",*) 'OCEAN_DYN_SLAB_nudging_basename is necessary !!'
          call PRC_abort
       endif
    else
       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Use nudging for SST : OFF'
    endif

    if ( OCEAN_DYN_SLAB_nudging ) then
       call FILE_EXTERNAL_INPUT_regist( OCEAN_DYN_SLAB_nudging_basename(:),           & ! [IN]
                                        'OCEAN_TEMP',                                 & ! [IN]
                                        'OXY',                                        & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_enable_periodic_year,  & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_enable_periodic_month, & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_enable_periodic_day,   & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_step_fixed,            & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_offset,                & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_defval,                & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_check_coordinates,     & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_step_limit             ) ! [IN]
    endif

    return
  end subroutine OCEAN_DYN_SLAB_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_DYN_SLAB( &
       OKMAX, OKS, OKE,  &
       OIA,   OIS, OIE,  &
       OJA,   OJS, OJE,  &
       OCEAN_TEMP_t,     &
       OCEAN_SFLX_G,     &
       OCEAN_SFLX_water, &
       OCEAN_SFLX_ice,   &
       calc_flag,        &
       dt, NOWDAYSEC,    &
       OCEAN_TEMP        )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       EMELT => CONST_EMELT
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    implicit none

    integer,  intent(in)    :: OKMAX, OKS, OKE
    integer,  intent(in)    :: OIA,   OIS, OIE
    integer,  intent(in)    :: OJA,   OJS, OJE
    real(RP), intent(in)    :: OCEAN_TEMP_t    (OKMAX,OIA,OJA) ! tendency of ocean temperature
    real(RP), intent(in)    :: OCEAN_SFLX_G    (OIA,OJA)       ! heat         flux from surface to subsurface (open ocean/sea ice)
    real(RP), intent(in)    :: OCEAN_SFLX_water(OIA,OJA)       ! liquid water flux from surface to subsurface (open ocean/sea ice)
    real(RP), intent(in)    :: OCEAN_SFLX_ice  (OIA,OJA)       ! ice    water flux from surface to subsurface (open ocean/sea ice)
    logical,  intent(in)    :: calc_flag       (OIA,OJA)       ! to decide calculate or not
    real(DP), intent(in)    :: dt
    real(DP), intent(in)    :: NOWDAYSEC
    real(RP), intent(inout) :: OCEAN_TEMP      (OKMAX,OIA,OJA)

    real(RP) :: OCEAN_TEMP_t_ndg(OKMAX,OIA,OJA)
    real(RP) :: OCEAN_TEMP_ref  (OKMAX,OIA,OJA)
    real(RP) :: rtau

    logical  :: error
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'ocean / dynamics / slab'

    if ( OCEAN_DYN_SLAB_nudging ) then

       call FILE_EXTERNAL_INPUT_update( 'OCEAN_TEMP', NOWDAYSEC, OCEAN_TEMP_ref(:,:,:), error )

       if ( error ) then
          LOG_ERROR("OCEAN_DYN_SLAB",*) 'Requested data is not found!'
          call PRC_abort
       endif

       ! if OCEAN_DYN_SLAB_nudging_tau < dt, Nudging acts as quasi-prescribed boundary
       rtau = 1.0_RP / max(OCEAN_DYN_SLAB_nudging_tausec,dt)

       do j = OJS, OJE
       do i = OIS, OIE
       do k = OKS, OKE
          OCEAN_TEMP_t_ndg(k,i,j) = ( OCEAN_TEMP_ref(k,i,j) - OCEAN_TEMP(k,i,j) ) * rtau
       enddo
       enddo
       enddo

    else
       OCEAN_TEMP_t_ndg(:,:,:) = 0.0_RP
    endif

    if ( OCEAN_DYN_SLAB_offline_mode ) then

       do j = OJS, OJE
       do i = OIS, OIE
          if ( calc_flag(i,j) ) then
             OCEAN_TEMP(OKS,i,j) = OCEAN_TEMP_ref(OKS,i,j)
          endif
       enddo
       enddo

    else

       do j = OJS, OJE
       do i = OIS, OIE
          if ( calc_flag(i,j) ) then
             ! heat flux from atm/ice at uppermost ocean layer
             OCEAN_TEMP(OKS,i,j) = OCEAN_TEMP(OKS,i,j) + OCEAN_TEMP_t_ndg(OKS,i,j) * dt &
                                 + ( OCEAN_SFLX_G(i,j) - OCEAN_SFLX_ice(i,j) * EMELT ) / OCEAN_DYN_SLAB_HeatCapacity * dt
             do k = OKS, OKE
                OCEAN_TEMP(k,i,j) = OCEAN_TEMP(k,i,j) + OCEAN_TEMP_t_ndg(k,i,j) * dt
             enddo
          endif
       enddo
       enddo

    endif

    return
  end subroutine OCEAN_DYN_SLAB

end module scale_ocean_dyn_slab
