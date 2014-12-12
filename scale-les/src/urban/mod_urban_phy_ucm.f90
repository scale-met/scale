!-------------------------------------------------------------------------------
!> module URBAN / Physics Urban Canopy Model (UCM)
!!
!! @par Description
!!          Urban physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_urban_phy_ucm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_PHY_driver_setup
  public :: URBAN_PHY_driver

  public :: URBAN_PHY_ucm_setup
  public :: URBAN_PHY_ucm

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
  subroutine URBAN_PHY_driver_setup
    use mod_urban_admin, only: &
       URBAN_TYPE, &
       URBAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN PHY] / Origin[SCALE-LES]'

    if ( URBAN_sw ) then

       ! setup library component
       call URBAN_PHY_ucm_setup( URBAN_TYPE )

       ! run once (only for the diagnostic value)
       call PROF_rapstart('URB Physics', 1)
       call URBAN_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('URB Physics', 1)

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine URBAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine URBAN_PHY_driver( update_flag )
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use mod_urban_vars, only: &
       URBAN_TEMP,  &
       URBAN_TEMP_t
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call URBAN_PHY_ucm( URBAN_TEMP  (:,:), & ! [IN]
                           URBAN_TEMP_t(:,:)  ) ! [OUT]

       call HIST_in( URBAN_TEMP_t(:,:), 'URBAN_TEMP_t', 'SST tendency', 'K' )

    endif

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, URBAN_TEMP_t(:,:), 'URBAN_TEMP_t' )
    endif

    return
  end subroutine URBAN_PHY_driver

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_PHY_ucm_setup( URBAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: URBAN_TYPE

    logical :: dummy

    NAMELIST / PARAM_URBAN_UCM / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[UCM] / Categ[URBAN PHYSICS] / Origin[SCALE-LES]'

    if ( URBAN_TYPE /= 'UCM' ) then
       write(*,*) 'xxx URBAN_TYPE is not UCM. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_UCM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_UCM. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_URBAN_UCM)

    return
  end subroutine URBAN_PHY_ucm_setup

  !-----------------------------------------------------------------------------
  !>  Urban canopy model
  subroutine URBAN_PHY_ucm( &
       URBAN_TEMP,   &
       URBAN_TEMP_t  )
    implicit none

    real(RP), intent(in)  :: URBAN_TEMP  (IA,JA)
    real(RP), intent(out) :: URBAN_TEMP_t(IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Urban step: UCM'

    do j = JS, JE
    do i = IS, IE
       URBAN_TEMP_t(i,j) = 0.0_RP
    enddo
    enddo

    return
  end subroutine URBAN_PHY_ucm

end module mod_urban_phy_ucm
