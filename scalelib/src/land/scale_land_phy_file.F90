!-------------------------------------------------------------------------------
!> module LAND / Physics File
!!
!! @par Description
!!          land physics module, external file input
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_phy_file
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_FILE_setup
  public :: LAND_PHY_FILE

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
  subroutine LAND_PHY_FILE_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_file_limit, &
       FILE_EXTERNAL_INPUT_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    character(len=H_LONG) :: LAND_PHY_FILE_basename(FILE_EXTERNAL_INPUT_file_limit) = ''
    logical               :: LAND_PHY_FILE_enable_periodic_year                     = .false.
    logical               :: LAND_PHY_FILE_enable_periodic_month                    = .false.
    logical               :: LAND_PHY_FILE_enable_periodic_day                      = .false.
    integer               :: LAND_PHY_FILE_step_fixed                               = 0
    real(RP)              :: LAND_PHY_FILE_offset                                   = 0.0_RP
    real(RP)              :: LAND_PHY_FILE_defval                                 ! = UNDEF
    logical               :: LAND_PHY_FILE_check_coordinates                        = .true.
    integer               :: LAND_PHY_FILE_step_limit                               = 0

    NAMELIST / PARAM_LAND_PHY_FILE / &
       LAND_PHY_FILE_basename,              &
       LAND_PHY_FILE_enable_periodic_year,  &
       LAND_PHY_FILE_enable_periodic_month, &
       LAND_PHY_FILE_enable_periodic_day,   &
       LAND_PHY_FILE_step_fixed,            &
       LAND_PHY_FILE_offset,                &
       LAND_PHY_FILE_defval,                &
       LAND_PHY_FILE_check_coordinates,     &
       LAND_PHY_FILE_step_limit

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FILE] / Categ[LAND PHY] / Origin[SCALElib]'

    if ( LAND_TYPE /= 'FILE' ) then
       write(*,*) 'xxx wrong LAND_TYPE. Check!'
       call PRC_MPIstop
    end if

    LAND_PHY_FILE_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_PHY_FILE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_PHY_FILE. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_PHY_FILE)

    if ( LAND_PHY_FILE_basename(1) == '' ) then
       write(*,*) 'xxx LAND_PHY_FILE_basename is necessary'
       call PRC_MPIstop
    end if

    call FILE_EXTERNAL_INPUT_regist( LAND_PHY_FILE_basename(:),           & ! [IN]
                                     'LAND_TEMP',                         & ! [IN]
                                     'LXY',                               & ! [IN]
                                     LAND_PHY_FILE_enable_periodic_year,  & ! [IN]
                                     LAND_PHY_FILE_enable_periodic_month, & ! [IN]
                                     LAND_PHY_FILE_enable_periodic_day,   & ! [IN]
                                     LAND_PHY_FILE_step_fixed,            & ! [IN]
                                     LAND_PHY_FILE_offset,                & ! [IN]
                                     LAND_PHY_FILE_defval,                & ! [IN]
                                     LAND_PHY_FILE_check_coordinates,     & ! [IN]
                                     LAND_PHY_FILE_step_limit             ) ! [IN]

    call FILE_EXTERNAL_INPUT_regist( LAND_PHY_FILE_basename(:),           & ! [IN]
                                     'LAND_WATER',                        & ! [IN]
                                     'LXY',                               & ! [IN]
                                     LAND_PHY_FILE_enable_periodic_year,  & ! [IN]
                                     LAND_PHY_FILE_enable_periodic_month, & ! [IN]
                                     LAND_PHY_FILE_enable_periodic_day,   & ! [IN]
                                     LAND_PHY_FILE_step_fixed,            & ! [IN]
                                     LAND_PHY_FILE_offset,                & ! [IN]
                                     LAND_PHY_FILE_defval,                & ! [IN]
                                     LAND_PHY_FILE_check_coordinates,     & ! [IN]
                                     LAND_PHY_FILE_step_limit             ) ! [IN]

    return
  end subroutine LAND_PHY_FILE_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_FILE( &
       TEMP_t,       &
       WATER_t,      &
       TEMP,         &
       WATER,        &
       WaterLimit,   &
       ThermalCond,  &
       HeatCapacity, &
       WaterDiff,    &
       SFLX_GH,      &
       SFLX_prec,    &
       SFLX_evap,    &
       CDZ,          &
       dt            )
    use scale_landuse, only: &
       LANDUSE_fact_land
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_process, only: &
       PRC_MPIstop
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    implicit none

    ! arguments
    real(RP), intent(out) :: TEMP_t      (LKMAX,IA,JA)
    real(RP), intent(out) :: WATER_t     (LKMAX,IA,JA)

    real(RP), intent(in)  :: TEMP        (LKMAX,IA,JA)
    real(RP), intent(in)  :: WATER       (LKMAX,IA,JA)
    real(RP), intent(in)  :: WaterLimit  (IA,JA)
    real(RP), intent(in)  :: ThermalCond (IA,JA)
    real(RP), intent(in)  :: HeatCapacity(IA,JA)
    real(RP), intent(in)  :: WaterDiff   (IA,JA)
    real(RP), intent(in)  :: SFLX_GH     (IA,JA)
    real(RP), intent(in)  :: SFLX_prec   (IA,JA)
    real(RP), intent(in)  :: SFLX_evap   (IA,JA)
    real(RP), intent(in)  :: CDZ         (LKMAX)
    real(DP), intent(in)  :: dt

    real(RP) :: LAND_TEMP_new (LKMAX,IA,JA)
    real(RP) :: LAND_WATER_new(LKMAX,IA,JA)

    logical :: error

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land physics step: File'

    call FILE_EXTERNAL_INPUT_update( &
         'LAND_TEMP',   & ! (in)
         NOWDAYSEC,     & ! (in)
         LAND_TEMP_new, & ! (out)
         error          ) ! (out)
    if ( error ) then
       write(*,*) 'xxx Requested data is not found!'
       call PRC_MPIstop
    end if

    call FILE_EXTERNAL_INPUT_update( &
         'LAND_WATER',   & ! (in)
         NOWDAYSEC,      & ! (in)
         LAND_WATER_new, & ! (out)
         error           ) ! (out)
    if ( error ) then
       write(*,*) 'xxx Requested data is not found!'
       call PRC_MPIstop
    end if

    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then
        TEMP_t (k,i,j) = ( LAND_TEMP_new (k,i,j) - TEMP (k,i,j) ) / dt
        WATER_t(k,i,j) = ( LAND_WATER_new(k,i,j) - WATER(k,i,j) ) / dt
      else
        TEMP_t (k,i,j) = 0.0_RP
        WATER_t(k,i,j) = 0.0_RP
      end if
    end do
    end do
    end do

    return
  end subroutine LAND_PHY_FILE

end module scale_land_phy_file
