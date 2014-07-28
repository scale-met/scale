!-------------------------------------------------------------------------------
!> module GRID (nesting system)
!!
!! @par Description
!!          Grid module for nesting system
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-28 (R.Yoshida)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_grid_nest
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: NEST_setup
  public :: NEST_domain_relate
  !public :: NEST_intercomm

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public              :: NEST_TILE_NUM_X      !< parent tile number in x-direction
  integer, public              :: NEST_TILE_NUM_Y      !< parent tile number in y-direction
  integer, public, allocatable :: NEST_TILE_ID  (:)    !< parent tile real id
  integer, public, allocatable :: NEST_TILE_ID_X(:)    !< tile id in x-direction
  integer, public, allocatable :: NEST_TILE_ID_Y(:)    !< tile id in y-direction

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  real(RP), private, allocatable :: latlon_catalog(:,:,:)  !< parent latlon catalog [rad]
  real(RP), private              :: corner_loc(4,2)        !< local corner location [rad]

  integer, private :: PARENT_PRC_NUM_X
  integer, private :: PARENT_PRC_NUM_Y
  integer, private :: PARENT_PRC_nmax
  integer, private :: PARENT_KMAX
  integer, private :: PARENT_IMAX
  integer, private :: PARENT_JMAX

  integer, parameter :: I_LON  = 1
  integer, parameter :: I_LAT  = 2
  integer, parameter :: I_NW   = 1
  integer, parameter :: I_NE   = 2
  integer, parameter :: I_SW   = 3
  integer, parameter :: I_SE   = 4

  logical, private :: OFFLINE = .true.

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine NEST_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_LONG) :: LATLON_CATALOGUE_FNAME = 'latlon_domain_catalogue.txt' !< metadata files for lat-lon domain for all processes

    namelist / PARAM_NEST /    &
       PARENT_PRC_NUM_X,       &
       PARENT_PRC_NUM_Y,       &
       PARENT_KMAX,            &
       PARENT_IMAX,            &
       PARENT_JMAX,            &
       OFFLINE,                &
       LATLON_CATALOGUE_FNAME

       debug

    integer :: parent_id
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NEST]/Categ[GRID]'

    call GRID_allocate

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_GRID,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_NEST. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_NEST)

    PARENT_PRC_nmax = PARENT_PRC_NUM_X * PARENT_PRC_NUM_Y
    allocate( latlon_catalog(PARENT_PRC_nmax,4,2) )

    corner_loc(I_NW,I_LON) = REAL_LON(IS,JE) / D2R
    corner_loc(I_NW,I_LAT) = REAL_LAT(IS,JE) / D2R
    corner_loc(I_NE,I_LON) = REAL_LON(IE,JE) / D2R
    corner_loc(I_NE,I_LAT) = REAL_LAT(IE,JE) / D2R
    corner_loc(I_SW,I_LON) = REAL_LON(IS,JS) / D2R
    corner_loc(I_SW,I_LAT) = REAL_LAT(IS,JS) / D2R
    corner_loc(I_SE,I_LON) = REAL_LON(IE,JS) / D2R
    corner_loc(I_SE,I_LAT) = REAL_LAT(IE,JS) / D2R

    if( OFFLINE ) then
       !--- read latlon catalogue
       fid = IO_get_available_fid()
       open( fid,                                    &
             file   = trim(LATLON_CATALOGUE_FNAME),  &
             form   = 'formatted',                   &
             status = 'old',                         &
             iostat = ierr                           )

       if ( ierr /= 0 ) then
          if( IO_L ) write(*,*) 'xxx cannot open latlon-catalogue file!'
          call PRC_MPIstop
       endif

       do i = 1, PRC_nmax
          read(fid,*,iostat=ierr) parent_id, &
                                  latlon_catalog(i,I_NW,I_LON), latlon_catalog(i,I_NE,I_LON), & ! LON: NW, NE
                                  latlon_catalog(i,I_SW,I_LON), latlon_catalog(i,I_SE,I_LON), & ! LON: SW, SE
                                  latlon_catalog(i,I_NW,I_LAT), latlon_catalog(i,I_NE,I_LAT), & ! LAT: NW, NE
                                  latlon_catalog(i,I_SW,I_LAT), latlon_catalog(i,I_SE,I_LAT)    ! LAT: SW, SE
          if ( i /= parent_id ) then
             if( IO_L ) write(*,*) 'xxx internal error: parent mpi id'
             call PRC_MPIstop
          endif
          if ( ierr /= 0 ) exit
       enddo
       close(fid)
    !else
       !
       ! ONLINE RELATIONSHIP
       !
    endif

    return
  end subroutine NEST_setup

  !-----------------------------------------------------------------------------
  !> Solve relationship between ParentDomain & Daughter Domain
  subroutine NEST_domain_relate
    implicit none

    integer, allocatable :: pd_tile_id_x(:,:)
    integer              :: pd_sw_tile, pd_ne_tile
    integer              :: i, j, ii, jj, k, hit
    !---------------------------------------------------------------------------

    allocate( pd_tile_id(PARENT_PRC_nmax,2) )
    jj = 1
    do j = 1, PARENT_PRC_NUM_Y
    ii = 1
    do i = 1, PARENT_PRC_NUM_X
       pd_tile_id(k,1) = i
       pd_tile_id(k,2) = j
    enddo
    enddo

    !--- SW search
    hit = 0
    do i = 1, PARENT_PRC_nmax
       if ( corner_loc(I_SW,I_LON) >= latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_SW,I_LON) <  latlon_catalog(i,I_NE,I_LON) .and. &
            corner_loc(I_SW,I_LAT) >= latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_SW,I_LAT) <  latlon_catalog(i,I_NE,I_LAT)       ) then

          pd_sw_tile = i
          hit = hit + 1
       endif
    enddo
    if ( hit == 0 ) then
       if( IO_L ) write(*,*) 'xxx domain mismatch between parent and daughter: SW search'
       call PRC_MPIstop
    elseif ( hit > 1 ) then
       if( IO_L ) write(*,*) 'xxx internal error: SW search'
       call PRC_MPIstop
    endif

    !--- NE search
    hit = 0
    do i = 1, PARENT_PRC_nmax 
       if ( corner_loc(I_NE,I_LON) >= latlon_catalog(i,I_SW,I_LON) .and. &
            corner_loc(I_NE,I_LON) <  latlon_catalog(i,I_NE,I_LON) .and. &
            corner_loc(I_NE,I_LAT) >= latlon_catalog(i,I_SW,I_LAT) .and. &
            corner_loc(I_NE,I_LAT) <  latlon_catalog(i,I_NE,I_LAT)       ) then

          pd_ne_tile = i
          hit = hit + 1
       endif
    enddo
    if ( hit == 0 ) then
       if( IO_L ) write(*,*) 'xxx domain mismatch between parent and daughter: NE search'
       call PRC_MPIstop
    elseif ( hit > 1 ) then
       if( IO_L ) write(*,*) 'xxx internal error: NE search'
       call PRC_MPIstop
    endif

    NEST_TILE_NUM_X = pd_tile_id(pd_ne_tile,1) - pd_tile_id(pd_sw_tile,1) + 1
    NEST_TILE_NUM_Y = pd_tile_id(pd_ne_tile,2) - pd_tile_id(pd_sw_tile,2) + 1
    allocate( NEST_TILE_ID  (NEST_TILE_NUM_X*NEST_TILE_NUM_Y) )
    allocate( NEST_TILE_ID_X(NEST_TILE_NUM_X*NEST_TILE_NUM_Y) )
    allocate( NEST_TILE_ID_Y(NEST_TILE_NUM_X*NEST_TILE_NUM_Y) )

    k = 1
    jj = 1
    do j = 1, NEST_TILE_NUM_Y
    ii = 1
    do i = 1, NEST_TILE_NUM_X
       if( ii==1 .and. jj==1 ) then
          NEST_TILE_ID(k) = pd_sw_tile
       elseif( ii > NEST_TILE_NUM_X ) then
          NEST_TILE_ID(k) = pd_sw_tile + NEST_TILE_NUM_X*(jj-1)
          ii = 1
       else
          NEST_TILE_ID(k) = NEST_TILE_ID(k-1) + 1
       endif

       NEST_TILE_ID_X(k) = i
       k = k + 1
    enddo
    NEST_TILE_ID_Y(k) = j
    enddo

    return
  end subroutine NEST_domain_relate

end module scale_grid_nest
