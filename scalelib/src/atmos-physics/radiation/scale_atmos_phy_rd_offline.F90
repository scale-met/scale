!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          for offline usage (input radiation flux from the file)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_rd_offline
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_offline_setup
  public :: ATMOS_PHY_RD_offline

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
  integer, private, parameter :: num_vars_3d    = 4
  integer, private, parameter :: num_vars_2d    = 4
  integer, private, parameter :: num_vars_2d_op = 1 ! optional

  real,    private :: ATMOS_PHY_RD_offline_diffuse_rate = 0.5_RP

  logical, private :: vars_2d_exist(num_vars_2d_op)

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_offline_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_external_input, only: &
       EXTIN_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=*), intent(in) :: RD_TYPE

    character(len=H_SHORT) :: vars_3d   (num_vars_3d)
    character(len=H_SHORT) :: vars_2d   (num_vars_2d)
    character(len=H_SHORT) :: vars_2d_op(num_vars_2d_op)

    data vars_3d    / 'RFLX_LW_up', 'RFLX_LW_dn', 'RFLX_SW_up', 'RFLX_SW_dn' /
    data vars_2d    / 'SFLX_LW_up', 'SFLX_LW_dn', 'SFLX_SW_up', 'SFLX_SW_dn' /
    data vars_2d_op / 'SFLX_SW_dn_dir' /

    character(len=H_LONG)  :: ATMOS_PHY_RD_offline_basename              = ''
    character(len=H_SHORT) :: ATMOS_PHY_RD_offline_axistype              = 'XYZ'
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_year  = .false.
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_month = .false.
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_day   = .false.
    integer                :: ATMOS_PHY_RD_offline_step_fixed            = 0
    real(RP)               :: ATMOS_PHY_RD_offline_offset                = 0.0_RP
    real(RP)               :: ATMOS_PHY_RD_offline_defval                !> = UNDEF
    logical                :: ATMOS_PHY_RD_offline_check_coordinates     = .true.
    integer                :: ATMOS_PHY_RD_offline_step_limit            = 0

    NAMELIST / PARAM_ATMOS_PHY_RD_OFFLINE / &
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
       ATMOS_PHY_RD_offline_diffuse_rate

    integer :: n, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[RADIATION] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Offline radiation process'

    if ( RD_TYPE /= 'OFFLINE' ) then
       write(*,*) 'xxx RD_TYPE is not OFFLINE. Check!'
       call PRC_MPIstop
    endif


    ATMOS_PHY_RD_offline_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_OFFLINE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD_OFFLINE. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_RD_OFFLINE)

    if ( ATMOS_PHY_RD_offline_basename == '' ) then
       write(*,*) 'xxx ATMOS_PHY_RD_offline_basename is necessary'
       call PRC_MPIstop
    end if

    do n = 1, num_vars_3d
       call EXTIN_regist( ATMOS_PHY_RD_offline_basename,              & ! [IN]
                          vars_3d(n),                                 & ! [IN]
                          ATMOS_PHY_RD_offline_axistype,              & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_year,  & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_month, & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_day,   & ! [IN]
                          ATMOS_PHY_RD_offline_step_fixed,            & ! [IN]
                          ATMOS_PHY_RD_offline_offset,                & ! [IN]
                          ATMOS_PHY_RD_offline_defval,                & ! [IN]
                          ATMOS_PHY_RD_offline_check_coordinates,     & ! [IN]
                          ATMOS_PHY_RD_offline_step_limit             ) ! [IN]
    end do

    do n = 1, num_vars_2d
       call EXTIN_regist( ATMOS_PHY_RD_offline_basename,              & ! [IN]
                          vars_2d(n),                                 & ! [IN]
                          'XY',                                       & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_year,  & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_month, & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_day,   & ! [IN]
                          ATMOS_PHY_RD_offline_step_fixed,            & ! [IN]
                          ATMOS_PHY_RD_offline_offset,                & ! [IN]
                          ATMOS_PHY_RD_offline_defval,                & ! [IN]
                          ATMOS_PHY_RD_offline_check_coordinates,     & ! [IN]
                          ATMOS_PHY_RD_offline_step_limit             ) ! [IN]
    end do

    do n = 1, num_vars_2d_op
       call EXTIN_regist( ATMOS_PHY_RD_offline_basename,              & ! [IN]
                          vars_2d_op(n),                              & ! [IN]
                          'XY',                                       & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_year,  & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_month, & ! [IN]
                          ATMOS_PHY_RD_offline_enable_periodic_day,   & ! [IN]
                          ATMOS_PHY_RD_offline_step_fixed,            & ! [IN]
                          ATMOS_PHY_RD_offline_offset,                & ! [IN]
                          ATMOS_PHY_RD_offline_defval,                & ! [IN]
                          ATMOS_PHY_RD_offline_check_coordinates,     & ! [IN]
                          ATMOS_PHY_RD_offline_step_limit,            & ! [IN]
                          exist = vars_2d_exist(n)                    ) ! [OUT]
       if ( vars_2d_exist(n) ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** ', trim(vars_2d_op(n)), ' found.'
       else
          if( IO_L ) write(IO_FID_LOG,*) '*** ', trim(vars_2d_op(n)), ' not found.'
       end if
    end do

    return
  end subroutine ATMOS_PHY_RD_offline_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_offline( &
       DENS, RHOT, QTRC,      &
       CZ, FZ,                &
       fact_ocean,            &
       fact_land,             &
       fact_urban,            &
       temp_sfc, albedo_land, &
       solins, cosSZA,        &
       flux_rad,              &
       flux_rad_top,          &
       SFLX_rad_dn            )
!       Jval                   )
    use scale_grid_index
    use scale_tracer
    use scale_process, only: &
       PRC_MPIstop
    use scale_external_input, only: &
       EXTIN_update
    use scale_time, only: &
       TIME_NOWDAYSEC
    use scale_atmos_phy_rd_common, only: &
       I_SW,     &
       I_LW,     &
       I_dn,     &
       I_up,     &
       I_direct, &
       I_diffuse
    implicit none
    real(RP), intent(in)  :: DENS        (KA,IA,JA)
    real(RP), intent(in)  :: RHOT        (KA,IA,JA)
    real(RP), intent(in)  :: QTRC        (KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ          (  KA,IA,JA)    ! UNUSED
    real(RP), intent(in)  :: FZ          (0:KA,IA,JA)
    real(RP), intent(in)  :: fact_ocean  (IA,JA)
    real(RP), intent(in)  :: fact_land   (IA,JA)
    real(RP), intent(in)  :: fact_urban  (IA,JA)
    real(RP), intent(in)  :: temp_sfc    (IA,JA)
    real(RP), intent(in)  :: albedo_land (IA,JA,2)
    real(RP), intent(in)  :: solins      (IA,JA)
    real(RP), intent(in)  :: cosSZA      (IA,JA)
    real(RP), intent(out) :: flux_rad    (KA,IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_top(IA,JA,2,2,2)
    real(RP), intent(out) :: SFLX_rad_dn (IA,JA,2,2)
!    real(RP), intent(out) :: Jval        (KA,IA,JA,CH_QA_photo)

    real(RP) :: buffer(IA,JA)
    logical  :: error, error_sum

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Radiation(offline)'

    ! [note] external data input is now support only SCALE-RM

    error_sum = .false.

    ! 3D
    call EXTIN_update( flux_rad(:,:,:,I_LW,I_up,2), 'RFLX_LW_up', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad(:,:,:,I_LW,I_dn,2), 'RFLX_LW_dn', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad(:,:,:,I_SW,I_up,2), 'RFLX_SW_up', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad(:,:,:,I_SW,I_dn,2), 'RFLX_SW_dn', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )


    ! 2D
    call EXTIN_update( buffer(:,:), 'SFLX_LW_up', TIME_NOWDAYSEC, error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_LW,I_up,2) = buffer(i,j)
       end do
       end do
    end if

    call EXTIN_update( buffer(:,:), 'SFLX_LW_dn', TIME_NOWDAYSEC, error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_LW,I_dn,2) = buffer(i,j)
       end do
       end do
    end if

    call EXTIN_update( buffer(:,:), 'SFLX_SW_up', TIME_NOWDAYSEC, error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_SW,I_up,2) = buffer(i,j)
       end do
       end do
    end if

    call EXTIN_update( buffer(:,:), 'SFLX_SW_dn', TIME_NOWDAYSEC, error )
    if ( error ) then
       error_sum = .true.
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE) &
       !$omp shared(flux_rad,buffer)
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_SW,I_dn,2) = buffer(i,j)
       end do
       end do
    end if


    ! 2D optional
    if ( vars_2d_exist(1) ) then
       call EXTIN_update( SFLX_rad_dn(:,:,I_SW,I_direct), 'SFLX_SW_dn_dir', TIME_NOWDAYSEC, error )
       if ( error ) then
          error_sum = .true.
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(i,j) &
          !$omp shared(KS,IS,IE,JS,JE) &
          !$omp shared(SFLX_rad_dn,flux_rad)
          do j = JS, JE
          do i = IS, IE
             SFLX_rad_dn(i,j,I_SW,I_diffuse) = flux_rad(KS-1,i,j,I_SW,I_dn,2) - SFLX_rad_dn(i,j,I_SW,I_direct)
          end do
          end do
       end if
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp private(i,j) &
       !$omp shared(KS,IS,IE,JS,JE,ATMOS_PHY_RD_offline_diffuse_rate) &
       !$omp shared(SFLX_rad_dn,flux_rad)
       do j = JS, JE
       do i = IS, IE
          SFLX_rad_dn(i,j,I_SW,I_diffuse) = (          ATMOS_PHY_RD_offline_diffuse_rate ) * flux_rad(KS-1,i,j,I_SW,I_dn,2)
          SFLX_rad_dn(i,j,I_SW,I_direct ) = ( 1.0_RP - ATMOS_PHY_RD_offline_diffuse_rate ) * flux_rad(KS-1,i,j,I_SW,I_dn,2)
       end do
       end do
    end if

    if ( error_sum ) then
       write(*,*) 'xxx Requested data is not found!'
       call PRC_MPIstop
    endif

    ! clearsky and TOA value are not defined
    flux_rad    (:,:,:,:,:,1) = 0.0_RP
    flux_rad_top(:,:,:,:,:)   = 0.0_RP

    !$omp parallel do default(none) OMP_SCHEDULE_ &
    !$omp private(i,j) &
    !$omp shared(KS,IS,IE,JS,JE,ATMOS_PHY_RD_offline_diffuse_rate) &
    !$omp shared(SFLX_rad_dn,flux_rad)
    do j = JS, JE
    do i = IS, IE
       SFLX_rad_dn(i,j,I_LW,I_diffuse) = flux_rad(KS-1,i,j,I_LW,I_dn,2) * ATMOS_PHY_RD_offline_diffuse_rate
       SFLX_rad_dn(i,j,I_LW,I_direct ) = flux_rad(KS-1,i,j,I_LW,I_dn,2) - SFLX_rad_dn(i,j,I_LW,I_diffuse)
    end do
    end do

    return
  end subroutine ATMOS_PHY_RD_offline

end module scale_atmos_phy_rd_offline
