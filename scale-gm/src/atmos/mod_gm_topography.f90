!-------------------------------------------------------------------------------
!> module TOPOGRAPHY
!!
!! @par Description
!!          Topography module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_gm_topography
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TOPO_setup
  public :: TOPO_fillhalo
  public :: TOPO_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: TOPO_exist = .false. !< topography exists?

  real(RP), public, allocatable :: TOPO_Zsfc(:,:) !< absolute ground height [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPO_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: TOPO_IN_BASENAME  = ''                     !< basename of the input  file
  logical,                private :: TOPO_IN_CHECK_COORDINATES = .false.        !> switch for check of coordinates
  character(len=H_LONG),  private :: TOPO_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),   private :: TOPO_OUT_TITLE    = 'SCALE-RM TOPOGRAPHY'  !< title    of the output file
  character(len=H_SHORT), private :: TOPO_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TOPO_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_TOPO / &
       TOPO_IN_BASENAME,          &
       TOPO_IN_CHECK_COORDINATES, &
       TOPO_OUT_BASENAME,         &
       TOPO_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TOPOGRAPHY] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_TOPO. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_TOPO)

    allocate( TOPO_Zsfc(IA,JA) )
    TOPO_Zsfc(:,:) = 0.0_RP

    ! read from file
    call TOPO_read

    return
  end subroutine TOPO_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine TOPO_fillhalo

    return
  end subroutine TOPO_fillhalo

  !-----------------------------------------------------------------------------
  !> Read topography
  subroutine TOPO_read

    return
  end subroutine TOPO_read

  !-----------------------------------------------------------------------------
  !> Write topography
  subroutine TOPO_write

    return
  end subroutine TOPO_write

end module mod_gm_topography
