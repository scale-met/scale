!-------------------------------------------------------------------------------
!> module land / grid / cartesianC
!!
!! @par Description
!!          Grid module for cartesian coordinate for land
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_grid_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_land_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_GRID_CARTESC_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: LAND_GRID_CARTESC_CZ (:) !< center coordinate [m]: z, local=global
  real(RP), public, allocatable :: LAND_GRID_CARTESC_FZ (:) !< face   coordinate [m]: z, local=global
  real(RP), public, allocatable :: LAND_GRID_CARTESC_CDZ(:) !< z-length of control volume [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_GRID_CARTESC_read
  private :: LAND_GRID_CARTESC_generate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: LDZ(100)

  character(len=H_LONG) :: LAND_GRID_CARTESC_IN_BASENAME  = ''
  logical               :: LAND_GRID_CARTESC_IN_AGGREGATE

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_GRID_CARTESC_setup
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_AGGREGATE
    implicit none

    namelist / PARAM_LAND_GRID_CARTESC / &
       LAND_GRID_CARTESC_IN_BASENAME,  &
       LAND_GRID_CARTESC_IN_AGGREGATE, &
       LDZ

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[CartesC] / Categ[LAND GRID] / Origin[SCALElib]'

    LDZ(:) = 0.0_RP

    LAND_GRID_CARTESC_IN_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_GRID_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_GRID_CARTESC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_GRID_CARTESC_setup",*) 'Not appropriate names in namelist PARAM_LAND_GRID_CARTESC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND_GRID_CARTESC)

    allocate( LAND_GRID_CARTESC_CZ (LKS  :LKE) )
    allocate( LAND_GRID_CARTESC_FZ (LKS-1:LKE) )
    allocate( LAND_GRID_CARTESC_CDZ(LKS  :LKE) )

    LOG_NEWLINE
    LOG_INFO("LAND_GRID_CARTESC_setup",*) 'Land grid information '

    if ( LAND_GRID_CARTESC_IN_BASENAME /= '' ) then
       call LAND_GRID_CARTESC_read
    else
       LOG_INFO("LAND_GRID_CARTESC_setup",*) 'Not found input grid file. Grid position is calculated.'

       call LAND_GRID_CARTESC_generate
    endif

    if ( LKE == LKS ) then
       LOG_NEWLINE
       LOG_INFO("LAND_GRID_CARTESC_setup",*) 'Single layer. LDZ = ', LDZ(1)
    else
       LOG_NEWLINE
       LOG_INFO_CONT('(1x,A)') 'Vertical Coordinate'
       LOG_INFO_CONT('(1x,A)')                  '|   k       z      zh      dz   k |'
       LOG_INFO_CONT('(1x,A)')                  '|         [m]     [m]     [m]     |'
       k = LKS-1
       LOG_INFO_CONT('(1x,A,F8.3,A,I4,A)')      '|            ',LAND_GRID_CARTESC_FZ(k),'        ',k,' | Atmosphere interface'
       do k = LKS, LKE-1
       LOG_INFO_CONT('(1x,A,I4,F8.3,A,F8.3,A)') '|',k,LAND_GRID_CARTESC_CZ(k),'        ',LAND_GRID_CARTESC_CDZ(k),'     | '
       LOG_INFO_CONT('(1x,A,F8.3,A,I4,A)')      '|            ',LAND_GRID_CARTESC_FZ(k),'       |',k,' | '
       enddo
       k = LKE
       LOG_INFO_CONT('(1x,A,I4,F8.3,A,F8.3,A)') '|',k,LAND_GRID_CARTESC_CZ(k),'        ',LAND_GRID_CARTESC_CDZ(k),'     | '
       LOG_INFO_CONT('(1x,A,F8.3,A,I4,A)')      '|            ',LAND_GRID_CARTESC_FZ(k),'        ',k,' | bedrock'
       LOG_INFO_CONT('(1x,A)')                  '|=================================|'
    endif

    return
  end subroutine LAND_GRID_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> Read land grid
  subroutine LAND_GRID_CARTESC_read
    use scale_file, only: &
       FILE_open, &
       FILE_read
    use scale_prc, only: &
       PRC_myrank
    implicit none

    integer :: fid
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_GRID_CARTESC_read",*) 'Input land grid file '

    call FILE_open( LAND_GRID_CARTESC_IN_BASENAME, fid, rankid=PRC_myrank, aggregate=LAND_GRID_CARTESC_IN_AGGREGATE )

    call FILE_read( fid, 'LCZ',  LAND_GRID_CARTESC_CZ (:) )
    call FILE_read( fid, 'LCDZ', LAND_GRID_CARTESC_CDZ(:) )
    call FILE_read( fid, 'LFZ',  LAND_GRID_CARTESC_FZ (:) )

    return
  end subroutine LAND_GRID_CARTESC_read

  !-----------------------------------------------------------------------------
  !> Generate land grid
  subroutine LAND_GRID_CARTESC_generate
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    do k = LKS, LKE
       LAND_GRID_CARTESC_CDZ(k) = LDZ(k)
    enddo

    LAND_GRID_CARTESC_FZ(LKS-1) = 0.0_RP

    do k = LKS, LKE
       LAND_GRID_CARTESC_CZ(k) = LAND_GRID_CARTESC_CDZ(k) / 2.0_RP + LAND_GRID_CARTESC_FZ(k-1)
       LAND_GRID_CARTESC_FZ(k) = LAND_GRID_CARTESC_CDZ(k)          + LAND_GRID_CARTESC_FZ(k-1)
    enddo

    return
  end subroutine LAND_GRID_CARTESC_generate

end module scale_land_grid_cartesC
