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

    integer, parameter     :: num_vars_2d = 6
    integer, parameter     :: num_vars_3d = 4
    character(len=H_SHORT) :: vars_2d(num_vars_2d)
    character(len=H_SHORT) :: vars_3d(num_vars_3d)
    data vars_2d / 'SFLX_LW_up', 'SFLX_LW_dn_dir', 'SFLX_LW_dn_dif', &
                   'SFLX_SW_up', "SFLX_SW_dn_dir", 'SFLX_SW_dn_dif'  /
    data vars_3d / 'RFLX_LW_up', 'RFLX_LW_dn', 'RFLX_SW_up', 'RFLX_SW_dn' /

    character(len=H_LONG)  :: ATMOS_PHY_RD_offline_basename = ''
    character(len=H_SHORT) :: ATMOS_PHY_RD_offline_axistype = 'XYZ'
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_year = .false.
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_month = .false.
    logical                :: ATMOS_PHY_RD_offline_enable_periodic_day = .false.
    integer                :: ATMOS_PHY_RD_offline_step_fixed = 0
    real(RP)               :: ATMOS_PHY_RD_offline_offset = 0.0_RP
    real(RP)               :: ATMOS_PHY_RD_offline_defval  !> = UNDEF
    logical                :: ATMOS_PHY_RD_offline_check_coordinates = .true.
    integer                :: ATMOS_PHY_RD_offline_step_limit = 0

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
            ATMOS_PHY_RD_offline_step_limit

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

    call EXTIN_update( buffer(:,:), 'SFLX_LW_up', TIME_NOWDAYSEC, error )
    if ( error ) then
       error_sum = .true.
    else
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_LW,I_up,2) = buffer(i,j)
       end do
       end do
    end if

    call EXTIN_update( buffer(:,:), 'SFLX_SW_up', TIME_NOWDAYSEC, error )
    if ( error ) then
       error_sum = .true.
    else
       do j = JS, JE
       do i = IS, IE
          flux_rad(KS-1,i,j,I_SW,I_up,2) = buffer(i,j)
       end do
       end do
    end if

    call EXTIN_update( SFLX_rad_dn (:,:,I_LW,I_direct ), 'SFLX_LW_dn_dir', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( SFLX_rad_dn (:,:,I_LW,I_diffuse), 'SFLX_LW_dn_dif', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( SFLX_rad_dn (:,:,I_SW,I_direct ), 'SFLX_SW_dn_dir', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( SFLX_rad_dn (:,:,I_SW,I_diffuse), 'SFLX_SW_dn_dif', TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad    (:,:,:,I_LW,I_up,2),  'RFLX_LW_up',     TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad    (:,:,:,I_LW,I_dn,2),  'RFLX_LW_dn',     TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad    (:,:,:,I_SW,I_up,2),  'RFLX_SW_up',     TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    call EXTIN_update( flux_rad    (:,:,:,I_SW,I_dn,2),  'RFLX_SW_dn',     TIME_NOWDAYSEC, error )
    error_sum = ( error .OR. error_sum )

    if ( error_sum ) then
       write(*,*) 'xxx Requested data is not found!'
       call PRC_MPIstop
    endif

    do j = JS, JE
    do i = IS, IE
       flux_rad(KS-1,i,j,I_LW,I_dn,2) = SFLX_rad_dn (i,j,I_LW,I_direct ) &
                                      + SFLX_rad_dn (i,j,I_LW,I_diffuse)
       flux_rad(KS-1,i,j,I_SW,I_dn,2) = SFLX_rad_dn (i,j,I_SW,I_direct ) &
                                      + SFLX_rad_dn (i,j,I_SW,I_diffuse)
    enddo
    enddo

    flux_rad    (:,:,:,:,:,1) = 0.0_RP
    flux_rad_top(:,:,:,:,:)   = 0.0_RP

    return
  end subroutine ATMOS_PHY_RD_offline

end module scale_atmos_phy_rd_offline
