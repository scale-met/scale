!-------------------------------------------------------------------------------
!> module OCEAN / Physics File
!!
!! @par Description
!!          ocean physics module, external file input
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_phy_file
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_FILE_setup
  public :: OCEAN_PHY_FILE

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
  logical, allocatable, private :: is_OCN(:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_FILE_setup( OCEAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    use scale_external_input, only: &
       EXTIN_regist
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE

    character(len=H_LONG) :: OCEAN_PHY_FILE_basename = ''
    logical               :: OCEAN_PHY_FILE_enable_periodic_year = .false.
    logical               :: OCEAN_PHY_FILE_enable_periodic_month = .false.
    logical               :: OCEAN_PHY_FILE_enable_periodic_day = .false.
    integer               :: OCEAN_PHY_FILE_step_fixed = 0
    real(RP)              :: OCEAN_PHY_FILE_offset = 0.0_RP
    real(RP)              :: OCEAN_PHY_FILE_defval  ! = UNDEF
    logical               :: OCEAN_PHY_FILE_check_coordinates = .true.
    integer               :: OCEAN_PHY_FILE_step_limit = 0

    NAMELIST / PARAM_OCEAN_PHY_FILE / &
            OCEAN_PHY_FILE_basename,              &
            OCEAN_PHY_FILE_enable_periodic_year,  &
            OCEAN_PHY_FILE_enable_periodic_month, &
            OCEAN_PHY_FILE_enable_periodic_day,   &
            OCEAN_PHY_FILE_step_fixed,            &
            OCEAN_PHY_FILE_offset,                &
            OCEAN_PHY_FILE_defval,                &
            OCEAN_PHY_FILE_check_coordinates,     &
            OCEAN_PHY_FILE_step_limit

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FILE] / Categ[OCEAN PHY] / Origin[SCALElib]'

    if( OCEAN_TYPE /= 'FILE' ) then
       write(*,*) 'xxx wrong OCEAN_TYPE. Check!'
       call PRC_MPIstop
    end if

    OCEAN_PHY_FILE_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_FILE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_PHY_FILE. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN_PHY_FILE)

    if ( OCEAN_PHY_FILE_basename == '' ) then
       write(*,*) 'xxx OCEAN_PHY_FILE_basename is necessary'
       call PRC_MPIstop
    end if

    call EXTIN_regist( &
         OCEAN_PHY_FILE_basename,              &
         'OCEAN_TEMP',                         &
         'XY',                                 &
         OCEAN_PHY_FILE_enable_periodic_year,  &
         OCEAN_PHY_FILE_enable_periodic_month, &
         OCEAN_PHY_FILE_enable_periodic_day,   &
         OCEAN_PHY_FILE_step_fixed,            &
         OCEAN_PHY_FILE_offset,                &
         OCEAN_PHY_FILE_defval,                &
         OCEAN_PHY_FILE_check_coordinates,     &
         OCEAN_PHY_FILE_step_limit             )

    ! judge to run slab ocean model
    allocate( is_OCN(IA,JA) )

    do j = JS, JE
    do i = IS, IE
       if( LANDUSE_fact_ocean(i,j) > 0.0_RP ) then
          is_OCN(i,j) = .true.
       else
          is_OCN(i,j) = .false.
       end if
    end do
    end do

    return
  end subroutine OCEAN_PHY_FILE_setup

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_PHY_FILE( &
       OCEAN_TEMP_t,    &
       OCEAN_TEMP,      &
       OCEAN_SFLX_WH,   &
       OCEAN_SFLX_prec, &
       OCEAN_SFLX_evap, &
       dt               )
    use scale_grid_index
    use scale_time, only: &
         NOWDAYSEC => TIME_NOWDAYSEC
    use scale_process, only: &
         PRC_MPIstop
    use scale_external_input, only: &
         EXTIN_update
    implicit none

    real(RP), intent(out) :: OCEAN_TEMP_t   (IA,JA)
    real(RP), intent(in)  :: OCEAN_TEMP     (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_WH  (IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_prec(IA,JA)
    real(RP), intent(in)  :: OCEAN_SFLX_evap(IA,JA)
    real(DP), intent(in)  :: dt

    real(RP) :: OCEAN_TEMP_new(IA,JA)

    logical :: error

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean physics step: File'

    call EXTIN_update( &
         OCEAN_TEMP_new, & ! (out)
         'OCEAN_TEMP',   & ! (in)
         NOWDAYSEC,      & ! (in)
         error           ) ! (out)
    if ( error ) then
       write(*,*) 'xxx Requested data is not found!'
       call PRC_MPIstop
    end if

    do j = JS, JE
    do i = IS, IE
       if( is_OCN(i,j) ) then
          OCEAN_TEMP_t(i,j) = ( OCEAN_TEMP_new(i,j) - OCEAN_TEMP(i,j) ) / dt
       else
          OCEAN_TEMP_t(i,j) = 0.0_RP
       endif
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_FILE

end module scale_ocean_phy_file
