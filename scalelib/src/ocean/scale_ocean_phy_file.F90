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
    implicit none

    character(len=*), intent(in) :: OCEAN_TYPE

    character(len=H_LONG) :: OCEAN_PHY_IN_BASENAME

    NAMELIST / PARAM_OCEAN_PHY_FILE / &
         OCEAN_PHY_IN_BASENAME

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FILE] / Categ[OCEAN PHY] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_FILE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_PHY_FILE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_PHY_FILE)

    if( OCEAN_TYPE .ne. 'FILE' ) then
       write(*,*) 'xxx wrong OCEAN_TYPE. Check!'
       call PRC_MPIstop
    end if

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input file name : ', OCEAN_PHY_IN_BASENAME

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

    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean step: File'

    call EXTIN_update( &
         OCEAN_TEMP_new, & ! (out)
         'OCEAN_TEMP',   & ! (in)
         'XY',           & ! (in)
         NOWDAYSEC,      & ! (in)
         error           ) ! (out)
    if ( error ) then
       write(*,*) 'xxx Failed to read "OCEAN_TEMP" data from file!'
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
