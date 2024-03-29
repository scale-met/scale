!-------------------------------------------------------------------------------
!> module ocean / dynamics / slab
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
  subroutine OCEAN_DYN_SLAB_setup( DEPTH )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       DWATR => CONST_DWATR
    use scale_atmos_hydrometeor, only: &
       CV_WATER
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_regist
    implicit none
    real(RP), intent(in) :: DEPTH

    real(DP)               :: OCEAN_DYN_SLAB_nudging_tau                   = 0.0_DP  ! Relaxation time
    character(len=H_SHORT) :: OCEAN_DYN_SLAB_nudging_tau_unit              = "SEC"
    character(len=H_LONG)  :: OCEAN_DYN_SLAB_nudging_basename              = ''
    logical                :: OCEAN_DYN_SLAB_nudging_basename_add_num      = .false.
    integer                :: OCEAN_DYN_SLAB_nudging_number_of_files       = 1
    logical                :: OCEAN_DYN_SLAB_nudging_enable_periodic_year  = .false.
    logical                :: OCEAN_DYN_SLAB_nudging_enable_periodic_month = .false.
    logical                :: OCEAN_DYN_SLAB_nudging_enable_periodic_day   = .false.
    integer                :: OCEAN_DYN_SLAB_nudging_step_fixed            = 0
    real(RP)               :: OCEAN_DYN_SLAB_nudging_defval                ! = UNDEF
    logical                :: OCEAN_DYN_SLAB_nudging_check_coordinates     = .true.
    integer                :: OCEAN_DYN_SLAB_nudging_step_limit            = 0

    ! obsolete
    real(RP)               :: OCEAN_DYN_SLAB_DEPTH = -1.0_RP

    namelist / PARAM_OCEAN_DYN_SLAB / &
       OCEAN_DYN_SLAB_nudging,                       &
       OCEAN_DYN_SLAB_nudging_tau,                   &
       OCEAN_DYN_SLAB_nudging_tau_unit,              &
       OCEAN_DYN_SLAB_nudging_basename,              &
       OCEAN_DYN_SLAB_nudging_basename_add_num,      &
       OCEAN_DYN_SLAB_nudging_number_of_files,       &
       OCEAN_DYN_SLAB_nudging_enable_periodic_year,  &
       OCEAN_DYN_SLAB_nudging_enable_periodic_month, &
       OCEAN_DYN_SLAB_nudging_enable_periodic_day,   &
       OCEAN_DYN_SLAB_nudging_step_fixed,            &
       OCEAN_DYN_SLAB_nudging_defval,                &
       OCEAN_DYN_SLAB_nudging_check_coordinates,     &
       OCEAN_DYN_SLAB_nudging_step_limit,            &
       OCEAN_DYN_SLAB_DEPTH

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

    if ( OCEAN_DYN_SLAB_DEPTH >= 0.0 ) then
       LOG_ERROR("OCEAN_DYN_SLAB_setup",*) '"OCEAN_DYN_SLAB_DEPTH" is obsolete. USE "ODZ" of "PARAM_OCEAN_GRID_CARTESC"'
       call PRC_abort
    end if

    OCEAN_DYN_SLAB_HeatCapacity = DWATR * CV_WATER * DEPTH

    LOG_NEWLINE
    LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Slab ocean depth [m]         : ', DEPTH
    LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Ocean heat capacity [J/K/m2] : ', OCEAN_DYN_SLAB_HeatCapacity

    if ( OCEAN_DYN_SLAB_nudging ) then
       call CALENDAR_unit2sec( OCEAN_DYN_SLAB_nudging_tausec, OCEAN_DYN_SLAB_nudging_tau, OCEAN_DYN_SLAB_nudging_tau_unit )

       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Use nudging for SST : ON'
       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Relaxation time Tau [sec] : ', OCEAN_DYN_SLAB_nudging_tausec

       if ( OCEAN_DYN_SLAB_nudging_tausec == 0.0_RP ) then
          OCEAN_DYN_SLAB_offline_mode = .true.
          LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Tau=0 means that SST is completely replaced by the external data.'
       endif

       if ( OCEAN_DYN_SLAB_nudging_basename == '' ) then
          LOG_ERROR("OCEAN_DYN_SLAB_setup",*) 'OCEAN_DYN_SLAB_nudging_basename is necessary !!'
          call PRC_abort
       endif
    else
       LOG_INFO("OCEAN_DYN_SLAB_setup",*) 'Use nudging for SST : OFF'
    endif

    if ( OCEAN_DYN_SLAB_nudging ) then
       call FILE_EXTERNAL_INPUT_regist( OCEAN_DYN_SLAB_nudging_basename,              & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_basename_add_num,      & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_number_of_files,       & ! [IN]
                                        'OCEAN_TEMP',                                 & ! [IN]
                                        'OXY',                                        & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_enable_periodic_year,  & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_enable_periodic_month, & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_enable_periodic_day,   & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_step_fixed,            & ! [IN]
                                        OCEAN_DYN_SLAB_nudging_defval,                & ! [IN]
                                        check_coordinates = OCEAN_DYN_SLAB_nudging_check_coordinates, & ! [IN]
                                        step_limit        = OCEAN_DYN_SLAB_nudging_step_limit,        & ! [IN]
                                        allow_missing     = ( .not. OCEAN_DYN_SLAB_offline_mode )     ) ! [IN]
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
       calc_flag,        &
       dt, NOWDAYSEC,    &
       OCEAN_TEMP,       &
       MASS_SUPL,        &
       ENGI_SUPL         )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_atmos_hydrometeor, only: &
       CV_WATER
    implicit none

    integer,  intent(in)    :: OKMAX, OKS, OKE
    integer,  intent(in)    :: OIA,   OIS, OIE
    integer,  intent(in)    :: OJA,   OJS, OJE
    real(RP), intent(in)    :: OCEAN_TEMP_t    (OKMAX,OIA,OJA) ! tendency of ocean temperature
    real(RP), intent(in)    :: OCEAN_SFLX_G    (OIA,OJA)       ! heat flux from surface to subsurface (open ocean/sea ice)
    real(RP), intent(in)    :: OCEAN_SFLX_water(OIA,OJA)       ! mass flux from surface to subsurface (open ocean/sea ice)
    logical,  intent(in)    :: calc_flag       (OIA,OJA)       ! to decide calculate or not
    real(DP), intent(in)    :: dt
    real(DP), intent(in)    :: NOWDAYSEC
    real(RP), intent(inout) :: OCEAN_TEMP      (OKMAX,OIA,OJA)
    real(RP), intent(out)   :: MASS_SUPL       (OIA,OJA)
    real(RP), intent(out)   :: ENGI_SUPL       (OIA,OJA)

    real(RP) :: OCEAN_TEMP_t_ndg(OKMAX,OIA,OJA)
    real(RP) :: OCEAN_TEMP_ref  (OKMAX,OIA,OJA)
    real(RP) :: rtau
    real(RP) :: dCP

    logical  :: error
    integer  :: k, i, j
    !---------------------------------------------------------------------------
    !$acc data copyin (OCEAN_TEMP_t,OCEAN_SFLX_G,OCEAN_SFLX_water,calc_flag) &
    !$acc      copyout(MASS_SUPL,ENGI_SUPL) &
    !$acc      copy(OCEAN_TEMP) &
    !$acc      create (OCEAN_TEMP_t_ndg,OCEAN_TEMP_ref)


    LOG_PROGRESS(*) 'ocean / dynamics / slab'

    if ( OCEAN_DYN_SLAB_nudging ) then

       call FILE_EXTERNAL_INPUT_update( 'OCEAN_TEMP', NOWDAYSEC, OCEAN_TEMP_ref(:,:,:), error )


       if ( error ) then
          LOG_ERROR("OCEAN_DYN_SLAB",*) 'Requested data is not found!'
          call PRC_abort
       endif

       if ( .not. OCEAN_DYN_SLAB_offline_mode ) then

          ! if OCEAN_DYN_SLAB_nudging_tau < dt, Nudging acts as quasi-prescribed boundary
          rtau = 1.0_RP / max(OCEAN_DYN_SLAB_nudging_tausec,dt)

          !$omp parallel do
          !$acc kernels
          do j = OJS, OJE
          do i = OIS, OIE
          do k = OKS, OKE
             if ( OCEAN_TEMP_ref(k,i,j) == UNDEF ) then
                OCEAN_TEMP_t_ndg(k,i,j) = 0.0_RP
             else
                OCEAN_TEMP_t_ndg(k,i,j) = ( OCEAN_TEMP_ref(k,i,j) - OCEAN_TEMP(k,i,j) ) * rtau
             end if
          enddo
          enddo
          enddo
          !$acc end kernels

       end if

    else
       !$omp parallel do
       !$acc kernels
       do j = OJS, OJE
       do i = OIS, OIE
       do k = OKS, OKE
          OCEAN_TEMP_t_ndg(k,i,j) = 0.0_RP
       end do
       end do
       end do
       !$acc end kernels
    endif

    if ( OCEAN_DYN_SLAB_offline_mode ) then

       !$omp parallel do
       !$acc kernels
       do j = OJS, OJE
       do i = OIS, OIE
          if ( calc_flag(i,j) ) then
             OCEAN_TEMP(OKS,i,j) = OCEAN_TEMP_ref(OKS,i,j)
          endif
          MASS_SUPL(i,j) = 0.0_RP
          ENGI_SUPL(i,j) = 0.0_RP
       enddo
       enddo
       !$acc end kernels

    else

       !$omp parallel do private(dCP)
       !$acc kernels
       do j = OJS, OJE
       do i = OIS, OIE
          if ( calc_flag(i,j) ) then
             ! heat flux from atm/ice at uppermost ocean layer
             dCP = CV_WATER * OCEAN_SFLX_water(i,j) * dt
             OCEAN_TEMP(OKS,i,j) = OCEAN_TEMP(OKS,i,j) &
                                 + ( OCEAN_SFLX_G(i,j) * dt - dCP * OCEAN_TEMP(OKS,i,j) ) &
                                   / ( OCEAN_DYN_SLAB_HeatCapacity + dCP ) &
                                 + OCEAN_TEMP_t_ndg(OKS,i,j) * dt
             do k = OKS+1, OKE
                OCEAN_TEMP(k,i,j) = OCEAN_TEMP(k,i,j) + OCEAN_TEMP_t_ndg(k,i,j) * dt
             enddo
             MASS_SUPL(i,j) = - OCEAN_SFLX_water(i,j)
             ENGI_SUPL(i,j) = CV_WATER * MASS_SUPL(i,j) * OCEAN_TEMP(OKS,i,j)
          else
             MASS_SUPL(i,j) = 0.0_RP
             ENGI_SUPL(i,j) = 0.0_RP
          endif
       enddo
       enddo
       !$acc end kernels

    endif
    !$acc end data

    return
  end subroutine OCEAN_DYN_SLAB

end module scale_ocean_dyn_slab
