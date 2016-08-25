!-------------------------------------------------------------------------------
!> module land grid index
!!
!! @par Description
!!          Grid Index module for land
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_grid_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_GRID_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: LKMAX = 1 ! # of computational cells: z for land

  integer, public :: LKS       ! start point of inner domain: z for land, local
  integer, public :: LKE       ! end   point of inner domain: z for land, local

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
  subroutine LAND_GRID_INDEX_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_LAND_INDEX / &
       LKMAX

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID_INDEX] / Categ[LAND GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_INDEX. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_INDEX)

    LKS  = 1
    LKE  = LKMAX

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Land grid index information ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** z-axis levels :', LKMAX

    return
  end subroutine LAND_GRID_INDEX_setup

end module scale_land_grid_index
