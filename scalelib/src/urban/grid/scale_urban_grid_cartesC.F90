!-------------------------------------------------------------------------------
!> module urban / grid / cartesianC
!!
!! @par Description
!!          Grid module for cartesian coordinate for urban
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_urban_grid_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_urban_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_GRID_CARTESC_setup
  public :: URBAN_GRID_CARTESC_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: URBAN_GRID_CARTESC_CZ  (:)  !< center coordinate [m]: z, local=global
  real(RP), public, allocatable :: URBAN_GRID_CARTESC_FZ  (:)  !< face   coordinate [m]: z, local=global
  real(RP), public, allocatable :: URBAN_GRID_CARTESC_CDZ (:)  !< z-length of control volume [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: URBAN_GRID_CARTESC_read
  private :: URBAN_GRID_CARTESC_generate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: UDZ(200)

  character(len=H_LONG) :: URBAN_GRID_CARTESC_IN_BASENAME  = ''
  logical               :: URBAN_GRID_CARTESC_IN_AGGREGATE

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_GRID_CARTESC_setup
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_AGGREGATE
    implicit none

    namelist / PARAM_URBAN_GRID_CARTESC / &
       URBAN_GRID_CARTESC_IN_BASENAME,  &
       URBAN_GRID_CARTESC_IN_AGGREGATE, &
       UDZ

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_GRID_CARTESC_setup",*) 'Setup'

    if ( UKMAX < 1 ) then
       LOG_INFO("URBAN_GRID_CARTESC_setup",*) 'Skip because UKMAX < 1'
       return
    end if

    UDZ(:) = 0.0_RP

    URBAN_GRID_CARTESC_IN_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_GRID_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("URBAN_GRID_CARTESC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("URBAN_GRID_CARTESC_setup",*) 'Not appropriate names in namelist PARAM_URBAN_GRID_CARTESC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_URBAN_GRID_CARTESC)

    allocate( URBAN_GRID_CARTESC_CZ (UKS  :UKE) )
    allocate( URBAN_GRID_CARTESC_FZ (UKS-1:UKE) )
    allocate( URBAN_GRID_CARTESC_CDZ(UKS  :UKE) )
    !$acc enter data create(URBAN_GRID_CARTESC_CZ,URBAN_GRID_CARTESC_FZ,URBAN_GRID_CARTESC_CDZ)

    LOG_NEWLINE
    LOG_INFO("URBAN_GRID_CARTESC_setup",*) 'Urban grid information '

    if ( URBAN_GRID_CARTESC_IN_BASENAME /= '' ) then
       call URBAN_GRID_CARTESC_read
    else
       LOG_INFO("URBAN_GRID_CARTESC_setup",*) 'Not found input grid file. Grid position is calculated.'

       call URBAN_GRID_CARTESC_generate
    endif

    if ( UKE == UKS ) then
       LOG_NEWLINE
       LOG_INFO("URBAN_GRID_CARTESC_setup",*) 'Single layer. LDZ = ', UDZ(1)
    else
       LOG_NEWLINE
       LOG_INFO("URBAN_GRID_CARTESC_setup",'(1x,A)') 'Vertical Coordinate'
       LOG_INFO_CONT('(1x,A)')                  '|   k       z      zh      dz   k |'
       LOG_INFO_CONT('(1x,A)')                  '|         [m]     [m]     [m]     |'
       k = UKS-1
       LOG_INFO_CONT('(1x,A,F8.3,A,I4,A)')      '|            ',URBAN_GRID_CARTESC_FZ(k),'        ',k,' | Atmosphere interface'
       do k = UKS, UKE-1
       LOG_INFO_CONT('(1x,A,I4,F8.3,A,F8.3,A)') '|',k,URBAN_GRID_CARTESC_CZ(k),'        ',URBAN_GRID_CARTESC_CDZ(k),'     | '
       LOG_INFO_CONT('(1x,A,F8.3,A,I4,A)')      '|            ',URBAN_GRID_CARTESC_FZ(k),'       |',k,' | '
       enddo
       k = UKE
       LOG_INFO_CONT('(1x,A,I4,F8.3,A,F8.3,A)') '|',k,URBAN_GRID_CARTESC_CZ(k),'        ',URBAN_GRID_CARTESC_CDZ(k),'     | '
       LOG_INFO_CONT('(1x,A,F8.3,A,I4,A)')      '|            ',URBAN_GRID_CARTESC_FZ(k),'        ',k,' | bedrock'
       LOG_INFO_CONT('(1x,A)')                  '|=================================|'
    endif

    !$acc update device(URBAN_GRID_CARTESC_CZ,URBAN_GRID_CARTESC_FZ,URBAN_GRID_CARTESC_CDZ)

    return
  end subroutine URBAN_GRID_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine URBAN_GRID_CARTESC_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_GRID_CARTESC_finalize",*) 'Finalize'

    !$acc exit data delete(URBAN_GRID_CARTESC_CZ,URBAN_GRID_CARTESC_FZ,URBAN_GRID_CARTESC_CDZ)
    deallocate( URBAN_GRID_CARTESC_CZ  )
    deallocate( URBAN_GRID_CARTESC_FZ  )
    deallocate( URBAN_GRID_CARTESC_CDZ )

    return
  end subroutine URBAN_GRID_CARTESC_finalize

  !-----------------------------------------------------------------------------
  !> Read urban grid
  subroutine URBAN_GRID_CARTESC_read
    use scale_file, only: &
       FILE_open, &
       FILE_read
    use scale_prc, only: &
       PRC_myrank
    implicit none

    integer :: fid
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_GRID_CARTESC_read",*) 'Input urban grid file '

    call FILE_open( URBAN_GRID_CARTESC_IN_BASENAME, fid, rankid=PRC_myrank, aggregate=URBAN_GRID_CARTESC_IN_AGGREGATE )

    call FILE_read( fid, 'UCZ',  URBAN_GRID_CARTESC_CZ(:)  )
    call FILE_read( fid, 'UCDZ', URBAN_GRID_CARTESC_CDZ(:) )
    call FILE_read( fid, 'UFZ',  URBAN_GRID_CARTESC_FZ(:)  )

    return
  end subroutine URBAN_GRID_CARTESC_read

  !-----------------------------------------------------------------------------
  !> Generate urban grid
  subroutine URBAN_GRID_CARTESC_generate
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    do k = UKS, UKE
       URBAN_GRID_CARTESC_CDZ(k) = UDZ(k)
    enddo

    URBAN_GRID_CARTESC_FZ(UKS-1) = 0.0_RP

    do k = UKS, UKE
       URBAN_GRID_CARTESC_CZ(k) = URBAN_GRID_CARTESC_CDZ(k) / 2.0_RP + URBAN_GRID_CARTESC_FZ(k-1)
       URBAN_GRID_CARTESC_FZ(k) = URBAN_GRID_CARTESC_CDZ(k)          + URBAN_GRID_CARTESC_FZ(k-1)
    enddo

    return
  end subroutine URBAN_GRID_CARTESC_generate

end module scale_urban_grid_cartesC
