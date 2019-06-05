!-------------------------------------------------------------------------------
!> module atmosphere / physics / radiation / offline
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          for offline usage (input radiation flux from the file)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_rd_offline
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_offline_setup
  public :: ATMOS_PHY_RD_offline_flux

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
  integer,  private, parameter :: num_vars_3d    = 4
  integer,  private, parameter :: num_vars_2d    = 4
  integer,  private, parameter :: num_vars_2d_op = 4 ! optional

  real(RP), private :: ATMOS_PHY_RD_offline_diffuse_rate = 0.5_RP
  real(RP), private :: ATMOS_PHY_RD_offline_NIR_rate     = 0.5_RP

  logical, private :: vars_2d_exist(num_vars_2d_op)

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_offline_setup
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_file_limit, &
       FILE_EXTERNAL_INPUT_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=H_SHORT) :: vars_3d   (num_vars_3d)
    character(len=H_SHORT) :: vars_2d   (num_vars_2d)
    character(len=H_SHORT) :: vars_2d_op(num_vars_2d_op)

    data vars_3d    / 'RFLX_LW_up', 'RFLX_LW_dn', 'RFLX_SW_up', 'RFLX_SW_dn' /
    data vars_2d    / 'SFLX_LW_up', 'SFLX_LW_dn', 'SFLX_SW_up', 'SFLX_SW_dn' /
    data vars_2d_op / 'SFLX_NIR_dn_dir', 'SFLX_NIR_dn_dif', 'SFLX_VIS_dn_dir', 'SFLX_VIS_dn_dif' /

    character(len=H_LONG)  :: ATMOS_PHY_RD_offline_basename(FILE_EXTERNAL_INPUT_file_limit) = ''
    character(len=H_SHORT) :: ATMOS_PHY_RD_offline_axistype                                 = 'XYZ'
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_year                     = .false.
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_month                    = .false.
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_day                      = .false.
    integer                :: ATMOS_PHY_RD_offline_step_fixed                               = 0
    real(RP)               :: ATMOS_PHY_RD_offline_offset                                   = 0.0_RP
    real(RP)               :: ATMOS_PHY_RD_offline_defval                                 ! = UNDEF
    logical                :: ATMOS_PHY_RD_offline_check_coordinates                        = .true.
    integer                :: ATMOS_PHY_RD_offline_step_limit                               = 0

    namelist / PARAM_ATMOS_PHY_RD_OFFLINE / &
       ATMOS_PHY_RD_offline_basename,              &
       ATMOS_PHY_RD_offline_axistype,              &
       ATMOS_PHY_RD_offline_enable_periodic_year,  &
       ATMOS_PHY_RD_offline_enable_periodic_month, &
       ATMOS_PHY_RD_offline_enable_periodic_day,   &
       ATMOS_PHY_RD_offline_step_fixed,            &
       ATMOS_PHY_RD_offline_offset,                &
       ATMOS_PHY_RD_offline_defval,                &
       ATMOS_PHY_RD_offline_check_coordinates,     &
       ATMOS_PHY_RD_offline_step_limit,            &
       ATMOS_PHY_RD_offline_diffuse_rate,          &
       ATMOS_PHY_RD_offline_NIR_rate

    integer :: n, ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_offline_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_RD_offline_setup",*) 'Offline radiation process'

    ATMOS_PHY_RD_offline_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_OFFLINE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_offline_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_offline_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD_OFFLINE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD_OFFLINE)

    if ( ATMOS_PHY_RD_offline_basename(1) == '' ) then
       LOG_ERROR("ATMOS_PHY_RD_offline_setup",*) 'ATMOS_PHY_RD_offline_basename is necessary'
       call PRC_abort
    end if

    do n = 1, num_vars_3d
       call FILE_EXTERNAL_INPUT_regist( ATMOS_PHY_RD_offline_basename(:),           & ! [IN]
                                        vars_3d(n),                                 & ! [IN]
                                        ATMOS_PHY_RD_offline_axistype,              & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_year,  & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_month, & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_day,   & ! [IN]
                                        ATMOS_PHY_RD_offline_step_fixed,            & ! [IN]
                                        ATMOS_PHY_RD_offline_offset,                & ! [IN]
                                        ATMOS_PHY_RD_offline_defval,                & ! [IN]
                                        check_coordinates = ATMOS_PHY_RD_offline_check_coordinates, & ! [IN]
                                        step_limit = ATMOS_PHY_RD_offline_step_limit                ) ! [IN]
    end do

    do n = 1, num_vars_2d
       call FILE_EXTERNAL_INPUT_regist( ATMOS_PHY_RD_offline_basename(:),           & ! [IN]
                                        vars_2d(n),                                 & ! [IN]
                                        'XY',                                       & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_year,  & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_month, & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_day,   & ! [IN]
                                        ATMOS_PHY_RD_offline_step_fixed,            & ! [IN]
                                        ATMOS_PHY_RD_offline_offset,                & ! [IN]
                                        ATMOS_PHY_RD_offline_defval,                & ! [IN]
                                        check_coordinates = ATMOS_PHY_RD_offline_check_coordinates, & ! [IN]
                                        step_limit = ATMOS_PHY_RD_offline_step_limit                ) ! [IN]
    end do

    do n = 1, num_vars_2d_op
       call FILE_EXTERNAL_INPUT_regist( ATMOS_PHY_RD_offline_basename(:),           & ! [IN]
                                        vars_2d_op(n),                              & ! [IN]
                                        'XY',                                       & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_year,  & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_month, & ! [IN]
                                        ATMOS_PHY_RD_offline_enable_periodic_day,   & ! [IN]
                                        ATMOS_PHY_RD_offline_step_fixed,            & ! [IN]
                                        ATMOS_PHY_RD_offline_offset,                & ! [IN]
                                        ATMOS_PHY_RD_offline_defval,                & ! [IN]
                                        check_coordinates = ATMOS_PHY_RD_offline_check_coordinates, & ! [IN]
                                        step_limit = ATMOS_PHY_RD_offline_step_limit,               & ! [IN]
                                        exist = vars_2d_exist(n)                                    ) ! [OUT]
       if ( vars_2d_exist(n) ) then
          LOG_INFO("ATMOS_PHY_RD_offline_setup",*) '', trim(vars_2d_op(n)), ' found.'
       else
          LOG_INFO("ATMOS_PHY_RD_offline_setup",*) '', trim(vars_2d_op(n)), ' not found.'
       end if
    end do

    return
  end subroutine ATMOS_PHY_RD_offline_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_offline_flux( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       time_now,   &
       flux_rad,   &
       SFLX_rad_dn )
    use scale_prc, only: &
       PRC_abort
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    use scale_atmos_phy_rd_common, only: &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(DP), intent(in)  :: time_now
    real(RP), intent(out) :: flux_rad   (KA,IA,JA,2,2)
    real(RP), intent(out) :: SFLX_rad_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)

    real(RP) :: buffer(IA,JA)
    logical  :: error, error_sum, error_sflx

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / radiation / offline'

    ! [note] external data input is now support only SCALE-RM

    error_sum = .false.

    ! 3D
    call FILE_EXTERNAL_INPUT_update( 'RFLX_LW_up', time_now, flux_rad(:,:,:,I_LW,I_up), error )
    error_sum = ( error .OR. error_sum )

    call FILE_EXTERNAL_INPUT_update( 'RFLX_LW_dn', time_now, flux_rad(:,:,:,I_LW,I_dn), error )
    error_sum = ( error .OR. error_sum )

    call FILE_EXTERNAL_INPUT_update( 'RFLX_SW_up', time_now, flux_rad(:,:,:,I_SW,I_up), error )
    error_sum = ( error .OR. error_sum )

    call FILE_EXTERNAL_INPUT_update( 'RFLX_SW_dn', time_now, flux_rad(:,:,:,I_SW,I_dn), error )
    error_sum = ( error .OR. error_sum )


    ! 2D
    call FILE_EXTERNAL_INPUT_update( 'SFLX_LW_up', time_now, buffer(:,:), error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_LW,I_up) = buffer(i,j)
       end do
       end do
    end if

    call FILE_EXTERNAL_INPUT_update( 'SFLX_LW_dn', time_now, buffer(:,:), error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_LW,I_dn) = buffer(i,j)
       end do
       end do
    end if

    call FILE_EXTERNAL_INPUT_update( 'SFLX_SW_up', time_now, buffer(:,:), error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_SW,I_up) = buffer(i,j)
       end do
       end do
    end if

    call FILE_EXTERNAL_INPUT_update( 'SFLX_SW_dn', time_now, buffer(:,:), error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_SW,I_dn) = buffer(i,j)
       end do
       end do
    end if

    !$omp parallel do default(none) OMP_SCHEDULE_ &
    !$omp private(i,j) &
    !$omp shared(KS,IS,IE,JS,JE,ATMOS_PHY_RD_offline_diffuse_rate) &
    !$omp shared(SFLX_rad_dn,flux_rad)
    do j = JS, JE
    do i = IS, IE
       SFLX_rad_dn(i,j,I_R_direct ,I_R_IR) = 0.0_RP
       SFLX_rad_dn(i,j,I_R_diffuse,I_R_IR) = flux_rad(KS-1,i,j,I_LW,I_dn)
    end do
    end do

    ! 2D optional

    error_sflx = .false.

    if ( vars_2d_exist(1) ) then
       call FILE_EXTERNAL_INPUT_update( 'SFLX_NIR_dn_dir', time_now, SFLX_rad_dn(:,:,I_R_direct,I_R_VIS), error )
       error_sum = ( error .OR. error_sum )
    else
       error_sflx = .true.
    endif

    if ( vars_2d_exist(2) ) then
       call FILE_EXTERNAL_INPUT_update( 'SFLX_NIR_dn_dif', time_now, SFLX_rad_dn(:,:,I_R_direct,I_R_VIS), error )
       error_sum = ( error .OR. error_sum )
    else
       error_sflx = .true.
    endif

    if ( vars_2d_exist(3) ) then
       call FILE_EXTERNAL_INPUT_update( 'SFLX_VIS_dn_dir', time_now, SFLX_rad_dn(:,:,I_R_direct,I_R_VIS), error )
       error_sum = ( error .OR. error_sum )
    else
       error_sflx = .true.
    endif

    if ( vars_2d_exist(4) ) then
       call FILE_EXTERNAL_INPUT_update( 'SFLX_VIS_dn_dif', time_now, SFLX_rad_dn(:,:,I_R_direct,I_R_VIS), error )
       error_sum = ( error .OR. error_sum )
    else
       error_sflx = .true.
    endif

    if ( error_sflx ) then ! reconstruct from lowermost SW flux
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE,ATMOS_PHY_RD_offline_diffuse_rate,ATMOS_PHY_RD_offline_NIR_rate) &
       !$omp shared(SFLX_rad_dn,flux_rad)
       do j = JS, JE
       do i = IS, IE
          SFLX_rad_dn(i,j,I_R_direct ,I_R_NIR) = ( 1.0_RP-ATMOS_PHY_RD_offline_diffuse_rate ) &
                                               * (        ATMOS_PHY_RD_offline_NIR_rate     ) * flux_rad(KS-1,i,j,I_SW,I_dn)
          SFLX_rad_dn(i,j,I_R_diffuse,I_R_NIR) = (        ATMOS_PHY_RD_offline_diffuse_rate ) &
                                               * (        ATMOS_PHY_RD_offline_NIR_rate     ) * flux_rad(KS-1,i,j,I_SW,I_dn)
          SFLX_rad_dn(i,j,I_R_direct ,I_R_VIS) = ( 1.0_RP-ATMOS_PHY_RD_offline_diffuse_rate ) &
                                               * ( 1.0_RP-ATMOS_PHY_RD_offline_NIR_rate     ) * flux_rad(KS-1,i,j,I_SW,I_dn)
          SFLX_rad_dn(i,j,I_R_diffuse,I_R_VIS) = (        ATMOS_PHY_RD_offline_diffuse_rate ) &
                                               * ( 1.0_RP-ATMOS_PHY_RD_offline_NIR_rate     ) * flux_rad(KS-1,i,j,I_SW,I_dn)
       enddo
       enddo
    endif

    if ( error_sum ) then
       LOG_ERROR("ATMOS_PHY_RD_offline_flux",*) 'Requested data is not found!'
       call PRC_abort
    endif

    return
  end subroutine ATMOS_PHY_RD_offline_flux

end module scale_atmos_phy_rd_offline
