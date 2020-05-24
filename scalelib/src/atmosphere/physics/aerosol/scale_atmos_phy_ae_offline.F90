!-------------------------------------------------------------------------------
!> module atmosphere / physics / aerosol / offline
!!
!! @par Description
!!          offline aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_offline
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_aerosol, only: &
      N_AE,  &
      I_A01, &
      I_A02, &
      I_A03, &
      I_A04, &
      I_A05, &
      I_A06, &
      I_A07
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_offline_setup
  public :: ATMOS_PHY_AE_offline_tendency
  public :: ATMOS_PHY_AE_offline_effective_radius
  public :: ATMOS_PHY_AE_offline_qtrc2qaero
  public :: ATMOS_PHY_AE_offline_mkinit

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
  integer,  private, parameter :: num_vars_3d = 8

  character(len=H_SHORT) :: vars_3d(num_vars_3d)

  data vars_3d / 'OUTQLD01', 'OUTQLD02', 'OUTQLD03', 'OUTQLD04', 'OUTQLD05', 'OUTQLD06', 'OUTQLD07', 'UNCCN' /

  real(RP), private :: const_value(num_vars_3d)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_offline_setup
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=H_LONG)  :: ATMOS_PHY_AE_offline_basename              = ''
    logical                :: ATMOS_PHY_AE_offline_basename_add_num      = .false.
    integer                :: ATMOS_PHY_AE_offline_number_of_files       = 1
    character(len=H_SHORT) :: ATMOS_PHY_AE_offline_axistype              = 'XYZ'
    logical                :: ATMOS_PHY_AE_offline_enable_periodic_year  = .false.
    logical                :: ATMOS_PHY_AE_offline_enable_periodic_month = .false.
    logical                :: ATMOS_PHY_AE_offline_enable_periodic_day   = .false.
    integer                :: ATMOS_PHY_AE_offline_step_fixed            = 0
    real(RP)               :: ATMOS_PHY_AE_offline_offset                = 0.0_RP
    real(RP)               :: ATMOS_PHY_AE_offline_defval               ! = UNDEF
    logical                :: ATMOS_PHY_AE_offline_check_coordinates     = .true.
    integer                :: ATMOS_PHY_AE_offline_step_limit            = 0

    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld01
    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld02
    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld03
    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld04
    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld05
    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld06
    real(RP)               :: ATMOS_PHY_AE_offline_const_outqld07
    real(RP)               :: ATMOS_PHY_AE_offline_const_unccn

    namelist / PARAM_ATMOS_PHY_AE_OFFLINE / &
       ATMOS_PHY_AE_offline_basename,              &
       ATMOS_PHY_AE_offline_basename_add_num,      &
       ATMOS_PHY_AE_offline_number_of_files,       &
       ATMOS_PHY_AE_offline_axistype,              &
       ATMOS_PHY_AE_offline_enable_periodic_year,  &
       ATMOS_PHY_AE_offline_enable_periodic_month, &
       ATMOS_PHY_AE_offline_enable_periodic_day,   &
       ATMOS_PHY_AE_offline_step_fixed,            &
       ATMOS_PHY_AE_offline_offset,                &
       ATMOS_PHY_AE_offline_defval,                &
       ATMOS_PHY_AE_offline_check_coordinates,     &
       ATMOS_PHY_AE_offline_step_limit,            &
       ATMOS_PHY_AE_offline_const_outqld01,        &
       ATMOS_PHY_AE_offline_const_outqld02,        &
       ATMOS_PHY_AE_offline_const_outqld03,        &
       ATMOS_PHY_AE_offline_const_outqld04,        &
       ATMOS_PHY_AE_offline_const_outqld05,        &
       ATMOS_PHY_AE_offline_const_outqld06,        &
       ATMOS_PHY_AE_offline_const_outqld07,        &
       ATMOS_PHY_AE_offline_const_unccn

    integer :: n, ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_offline_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_AE_offline_setup",*) 'Kajino(2013) scheme'

    ATMOS_PHY_AE_offline_defval = UNDEF

    ATMOS_PHY_AE_offline_const_outqld01 = UNDEF
    ATMOS_PHY_AE_offline_const_outqld02 = UNDEF
    ATMOS_PHY_AE_offline_const_outqld03 = UNDEF
    ATMOS_PHY_AE_offline_const_outqld04 = UNDEF
    ATMOS_PHY_AE_offline_const_outqld05 = UNDEF
    ATMOS_PHY_AE_offline_const_outqld06 = UNDEF
    ATMOS_PHY_AE_offline_const_outqld07 = UNDEF
    ATMOS_PHY_AE_offline_const_unccn    = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_AE_offline,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_AE_offline_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_AE_offline_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_offline. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_offline)

    const_value(1) = ATMOS_PHY_AE_offline_const_outqld01
    const_value(2) = ATMOS_PHY_AE_offline_const_outqld02
    const_value(3) = ATMOS_PHY_AE_offline_const_outqld03
    const_value(4) = ATMOS_PHY_AE_offline_const_outqld04
    const_value(5) = ATMOS_PHY_AE_offline_const_outqld05
    const_value(6) = ATMOS_PHY_AE_offline_const_outqld06
    const_value(7) = ATMOS_PHY_AE_offline_const_outqld07
    const_value(8) = ATMOS_PHY_AE_offline_const_unccn

    do n = 1, num_vars_3d

       if ( const_value(n) < UNDEF*0.1 ) then ! read from external file

          if ( ATMOS_PHY_AE_offline_basename /= '' ) then
             call FILE_EXTERNAL_INPUT_regist( ATMOS_PHY_AE_offline_basename,              & ! [IN]
                                              ATMOS_PHY_AE_offline_basename_add_num,      & ! [IN]
                                              ATMOS_PHY_AE_offline_number_of_files,       & ! [IN]
                                              vars_3d(n),                                 & ! [IN]
                                              ATMOS_PHY_AE_offline_axistype,              & ! [IN]
                                              ATMOS_PHY_AE_offline_enable_periodic_year,  & ! [IN]
                                              ATMOS_PHY_AE_offline_enable_periodic_month, & ! [IN]
                                              ATMOS_PHY_AE_offline_enable_periodic_day,   & ! [IN]
                                              ATMOS_PHY_AE_offline_step_fixed,            & ! [IN]
                                              ATMOS_PHY_AE_offline_offset,                & ! [IN]
                                              ATMOS_PHY_AE_offline_defval,                & ! [IN]
                                              check_coordinates = ATMOS_PHY_AE_offline_check_coordinates, & ! [IN]
                                              step_limit        = ATMOS_PHY_AE_offline_step_limit         ) ! [IN]
          endif

       else ! set constant value

          LOG_INFO("ATMOS_PHY_AE_offline_setup",*) &
            'Constant value is set for ', trim(vars_3d(n)), ', value = ', const_value(n)

       endif

    enddo

    return
  end subroutine ATMOS_PHY_AE_offline_setup

  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics
  subroutine ATMOS_PHY_AE_offline_tendency( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       time_now,   &
       CCN         )
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(DP), intent(in)  :: time_now
    real(RP), intent(out) :: CCN(KA,IA,JA)

    logical  :: error
    integer  :: n
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / aerosol / offline'

    n = 8 ! UNCCN

    if ( const_value(n) < UNDEF*0.1 ) then ! read from external file

       call FILE_EXTERNAL_INPUT_update( vars_3d(n), time_now, CCN(:,:,:), error )
       if ( error ) then
          LOG_ERROR("ATMOS_PHY_AE_offline_flux",*) 'Requested data is not found! ', trim(vars_3d(n))
          call PRC_abort
       endif

    else ! set constant value

       CCN(:,:,:) = const_value(n)

    endif

    return
  end subroutine ATMOS_PHY_AE_offline_tendency

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AE_offline_effective_radius( &
       KA, IA, JA, &
       RH,         &
       Re          )
    implicit none

    integer,  intent(in)  :: KA, IA, JA
    real(RP), intent(in)  :: RH(KA,IA,JA)      ! relative humidity (0-1)
    real(RP), intent(out) :: Re(KA,IA,JA,N_AE) ! effective radius

    real(RP), parameter :: AE_Re(N_AE) = & ! aerosol radius [m]
         (/ 1.6E-6_RP, &  ! Soil dust
              -1.0_RP, &  ! Carbonacerous (BC/OC=0.3)
              -1.0_RP, &  ! Carbonacerous (BC/OC=0.15)
              -1.0_RP, &  ! Carbonacerous (BC/OC=0.)
             4.E-8_RP, &  ! Black carbon
              -1.0_RP, &  ! Sulfate
              -1.0_RP  /) ! Sea salt

    integer  :: iaero
    !---------------------------------------------------------------------------

    do iaero = 1, N_AE

       if ( AE_Re(iaero) < 0.0_RP ) then ! hygroscopic particle : look-up table is based on the RH

          Re(:,:,:,iaero) = RH(:,:,:)

       else                              ! non-hygroscopic particle : look-up table is the effective radius

          Re(:,:,:,iaero) = AE_Re(iaero) * 100.0_RP ! [m=>cm]

       endif

    enddo

    return
  end subroutine ATMOS_PHY_AE_offline_effective_radius

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AE_offline_qtrc2qaero( &
       KA, IA, JA, &
       time_now,   &
       Qe          )
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    integer,  intent(in)  :: KA, IA, JA
    real(DP), intent(in)  :: time_now
    real(RP), intent(out) :: Qe(KA,IA,JA,N_AE) ! aerosol mixing ratio [kg/kg]

    logical  :: error, error_sum
    integer  :: n, iaero
    !---------------------------------------------------------------------------

    error_sum = .false.

    do n = 1, num_vars_3d-1
       iaero = n

       if ( const_value(n) < UNDEF*0.1 ) then ! read from external file

          call FILE_EXTERNAL_INPUT_update( vars_3d(n), time_now, Qe(:,:,:,iaero), error )
          error_sum = ( error .OR. error_sum )

       else ! set constant value

          Qe(:,:,:,iaero) = const_value(n)

       endif

    enddo

    if ( error_sum ) then
       LOG_ERROR("ATMOS_PHY_AE_offline_flux",*) 'Requested data is not found!'
       call PRC_abort
    endif

    return
  end subroutine ATMOS_PHY_AE_offline_qtrc2qaero

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_AE_offline_mkinit( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       ccn_init,   &
       CCN         )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: ccn_init
    real(RP), intent(out) :: CCN(KA,IA,JA)

    integer  :: n
    !---------------------------------------------------------------------------

    call ATMOS_PHY_AE_offline_setup

    n = 8 ! UNCCN

    if ( const_value(n) < UNDEF*0.1 ) then ! read from external file

       CCN(:,:,:) = ccn_init

    else ! set constant value

       CCN(:,:,:) = const_value(n)

    endif

    return
  end subroutine ATMOS_PHY_AE_offline_mkinit

end module scale_atmos_phy_ae_offline
