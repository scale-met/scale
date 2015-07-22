!-------------------------------------------------------------------------------
!> module GRID (real space)
!!
!! @par Description
!!          Grid module for orthogonal curvelinear, terrain-following coordinate
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-24 (H.Yashiro)  [new] reconstruction from scale_REAL & scale_topography
!!
!<
!-------------------------------------------------------------------------------
module scale_grid_real
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
  public :: REAL_setup
  public :: REAL_update_Z
  public :: REAL_calc_areavol

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: REAL_LON(:,:)       !< longitude [rad,0-2pi]
  real(RP), public, allocatable :: REAL_LAT(:,:)       !< latitude  [rad,-pi,pi]
  real(RP), public, allocatable :: REAL_CZ (:,:,:)     !< geopotential height [m] (cell center)
  real(RP), public, allocatable :: REAL_FZ (:,:,:)     !< geopotential height [m] (cell face  )

  real(RP), public              :: REAL_BASEPOINT_LON  !< position of base point in real world [rad,0-2pi]
  real(RP), public              :: REAL_BASEPOINT_LAT  !< position of base point in real world [rad,-pi,pi]

  real(RP), public, allocatable :: REAL_LONX (:,:)     !< longitude at staggered point (uy) [rad,0-2pi]
  real(RP), public, allocatable :: REAL_LONY (:,:)     !< longitude at staggered point (xv) [rad,0-2pi]
  real(RP), public, allocatable :: REAL_LONXY(:,:)     !< longitude at staggered point (uv) [rad,0-2pi]
  real(RP), public, allocatable :: REAL_LATX (:,:)     !< latitude  at staggered point (uy) [rad,-pi,pi]
  real(RP), public, allocatable :: REAL_LATY (:,:)     !< latitude  at staggered point (xv) [rad,-pi,pi]
  real(RP), public, allocatable :: REAL_LATXY(:,:)     !< latitude  at staggered point (uv) [rad,-pi,pi]
  real(RP), public, allocatable :: REAL_DLON (:,:)     !< delta longitude
  real(RP), public, allocatable :: REAL_DLAT (:,:)     !< delta latitude

  real(RP), public, allocatable :: REAL_Z1  (:,:)      !< Height of the lowermost grid from surface (cell center) [m]
  real(RP), public              :: REAL_ASPECT_MAX     !< maximum aspect ratio of the grid cell
  real(RP), public              :: REAL_ASPECT_MIN     !< minimum aspect ratio of the grid cell

  real(RP), public, allocatable :: REAL_PHI (:,:,:)    !< geopotential [m2/s2] (cell center)

  real(RP), public, allocatable :: REAL_AREA(:,:)      !< horizontal area [m2]
  real(RP), public, allocatable :: REAL_VOL (:,:,:)    !< control volume  [m3]

  real(RP), public, allocatable :: REAL_DOMAIN_CATALOGUE(:,:,:) !< domain latlon catalogue [rad]

  real(RP), public :: REAL_TOTAREA                     !< total area   (local) [m2]
  real(RP), public :: REAL_TOTVOL                      !< total volume (local) [m3]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: REAL_calc_latlon
  private :: REAL_calc_Z

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine REAL_setup
    use scale_process, only: &
       PRC_nmax,    &
       PRC_MPIstop
    use scale_grid, only: &
       GRID_DOMAIN_CENTER_X, &
       GRID_DOMAIN_CENTER_Y
    use scale_mapproj, only: &
       MPRJ_setup
    use scale_fileio, only: &
       FILEIO_set_coordinates
    implicit none

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[REAL]/Categ[GRID]'

    allocate( REAL_LON  (IA,JA) )
    allocate( REAL_LAT  (IA,JA) )
    allocate( REAL_LONX (IA,JA) )
    allocate( REAL_LONY (IA,JA) )
    allocate( REAL_LONXY(IA,JA) )
    allocate( REAL_LATX (IA,JA) )
    allocate( REAL_LATY (IA,JA) )
    allocate( REAL_LATXY(IA,JA) )
    allocate( REAL_DLON (IA,JA) )
    allocate( REAL_DLAT (IA,JA) )

    allocate( REAL_CZ (  KA,IA,JA) )
    allocate( REAL_FZ (0:KA,IA,JA) )
    allocate( REAL_Z1 (     IA,JA) )
    allocate( REAL_PHI(  KA,IA,JA) )

    allocate( REAL_AREA(   IA,JA) )
    allocate( REAL_VOL (KA,IA,JA) )

    allocate( REAL_DOMAIN_CATALOGUE(PRC_nmax,4,2) )

    ! setup map projection
    call MPRJ_setup( GRID_DOMAIN_CENTER_X, GRID_DOMAIN_CENTER_Y )

    ! calc longitude & latitude
    call REAL_calc_latlon

    ! calc real height
    call REAL_calc_Z

    ! calc control area & volume
    ! call REAL_calc_areavol ! must be called after GTRANS_setup

    ! set latlon and z to fileio module
    call FILEIO_set_coordinates( REAL_LON, REAL_LONX, REAL_LONY, REAL_LONXY, &
                                 REAL_LAT, REAL_LATX, REAL_LATY, REAL_LATXY, &
                                 REAL_CZ,  REAL_FZ    )

    return
  end subroutine REAL_setup

  !-----------------------------------------------------------------------------
  !> Re-setup with updated topography
  subroutine REAL_update_Z
    use scale_process, only: &
       PRC_MPIstop
    use scale_fileio, only: &
       FILEIO_set_coordinates
    implicit none
    !---------------------------------------------------------------------------

    ! calc real height
    call REAL_calc_Z

    ! set latlon and z to fileio module
    call FILEIO_set_coordinates( REAL_LON, REAL_LONX, REAL_LONY, REAL_LONXY, &
                                 REAL_LAT, REAL_LATX, REAL_LATY, REAL_LATXY, &
                                 REAL_CZ,  REAL_FZ    )

    return
  end subroutine REAL_update_Z

  !-----------------------------------------------------------------------------
  !> Calc longitude & latitude
  subroutine REAL_calc_latlon
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_grid, only: &
       GRID_DOMAIN_CENTER_X, &
       GRID_DOMAIN_CENTER_Y, &
       GRID_CX, &
       GRID_CY, &
       GRID_FX, &
       GRID_FY
    use scale_mapproj, only: &
       MPRJ_basepoint_lon, &
       MPRJ_basepoint_lat, &
       MPRJ_xy2lonlat
    use scale_process, only: &
       PRC_master,           &
       PRC_myrank,           &
       PRC_nmax,             &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_gather,  &
       COMM_bcast
    implicit none

    real(RP), allocatable :: mine(:,:)    !< send buffer of lon-lat [deg]
    real(RP), allocatable :: whole(:,:)   !< recieve buffer of lon-lat [deg]

    real(RP) :: CX, CY, FX, FY

    integer, parameter :: I_LON  = 1
    integer, parameter :: I_LAT  = 2
    integer, parameter :: I_NW   = 1
    integer, parameter :: I_NE   = 2
    integer, parameter :: I_SW   = 3
    integer, parameter :: I_SE   = 4

    character(len=H_LONG) :: DOMAIN_CATALOGUE_FNAME = 'latlon_domain_catalogue.txt' !< metadata files for lat-lon domain for all processes
    logical               :: DOMAIN_CATALOGUE_OUTPUT = .false.

    NAMELIST / PARAM_DOMAIN_CATALOGUE / &
       DOMAIN_CATALOGUE_FNAME,  &
       DOMAIN_CATALOGUE_OUTPUT

    integer :: i, j
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_DOMAIN_CATALOGUE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_DOMAIN_CATALOGUE. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_DOMAIN_CATALOGUE)

    REAL_BASEPOINT_LON = MPRJ_basepoint_lon * D2R
    REAL_BASEPOINT_LAT = MPRJ_basepoint_lat * D2R

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '*** Base position in the global domain'
    if( IO_L ) write(IO_FID_LOG,*) '->(',REAL_BASEPOINT_LON/D2R,',',REAL_BASEPOINT_LAT/D2R,')'

    do j = 1, JA
    do i = 1, IA
       ! apply offset
       CX = GRID_CX(i) - GRID_DOMAIN_CENTER_X
       CY = GRID_CY(j) - GRID_DOMAIN_CENTER_Y
       FX = GRID_FX(i) - GRID_DOMAIN_CENTER_X
       FY = GRID_FY(j) - GRID_DOMAIN_CENTER_Y

       call MPRJ_xy2lonlat( GRID_CX(i), GRID_CY(j), REAL_LON  (i,j), REAL_LAT  (i,j) )
       call MPRJ_xy2lonlat( GRID_FX(i), GRID_CY(j), REAL_LONX (i,j), REAL_LATX (i,j) )
       call MPRJ_xy2lonlat( GRID_CX(i), GRID_FY(j), REAL_LONY (i,j), REAL_LATY (i,j) )
       call MPRJ_xy2lonlat( GRID_FX(i), GRID_FY(j), REAL_LONXY(i,j), REAL_LATXY(i,j) )
    enddo
    enddo

    REAL_DLON(:,:) = 0.0_RP
    REAL_DLAT(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       REAL_DLON(i,j) = REAL_LONX(i,j) - REAL_LONX(i-1,j)
       REAL_DLAT(i,j) = REAL_LATY(i,j) - REAL_LATY(i,j-1)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '*** Position on the earth (Local)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f10.5,A,f9.5,A,A,f10.5,A,f9.5,A)') &
                                'NW(',REAL_LON(IS,JE)/D2R,',',REAL_LAT(IS,JE)/D2R,')-', &
                                'NE(',REAL_LON(IE,JE)/D2R,',',REAL_LAT(IE,JE)/D2R,')'
    if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
                                '            |                       |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f10.5,A,f9.5,A,A,f10.5,A,f9.5,A)') &
                                'SW(',REAL_LON(IS,JS)/D2R,',',REAL_LAT(IS,JS)/D2R,')-', &
                                'SE(',REAL_LON(IE,JS)/D2R,',',REAL_LAT(IE,JS)/D2R,')'


    allocate( mine (4, 2         ) )
    allocate( whole(4, 2*PRC_nmax) )

    mine(I_NW,I_LON) = REAL_LONXY(IS-1,JE  )/D2R
    mine(I_NE,I_LON) = REAL_LONXY(IE  ,JE  )/D2R
    mine(I_SW,I_LON) = REAL_LONXY(IS-1,JS-1)/D2R
    mine(I_SE,I_LON) = REAL_LONXY(IE  ,JS-1)/D2R
    mine(I_NW,I_LAT) = REAL_LATXY(IS-1,JE  )/D2R
    mine(I_NE,I_LAT) = REAL_LATXY(IE  ,JE  )/D2R
    mine(I_SW,I_LAT) = REAL_LATXY(IS-1,JS-1)/D2R
    mine(I_SE,I_LAT) = REAL_LATXY(IE  ,JS-1)/D2R

    call COMM_gather( whole, mine, 4, 2 ) ! everytime do for online nesting

    if( PRC_myrank == PRC_master )then
       if( DOMAIN_CATALOGUE_OUTPUT ) then
          fid = IO_get_available_fid()
          open( fid,                                    &
                file   = trim(DOMAIN_CATALOGUE_FNAME), &
                form   = 'formatted',                  &
                status = 'replace',                    &
                iostat = ierr                          )

          if ( ierr /= 0 ) then
             if( IO_L ) write(*,*) 'xxx cannot create latlon-catalogue file!'
             call PRC_MPIstop
          endif

          do i = 1, PRC_nmax ! for offline nesting
             write(fid,'(i8,8f32.24)',iostat=ierr) i, whole(I_NW,I_LON+2*(i-1)), whole(I_NE,I_LON+2*(i-1)), & ! LON: NW, NE
                                                      whole(I_SW,I_LON+2*(i-1)), whole(I_SE,I_LON+2*(i-1)), & ! LON: SW, SE
                                                      whole(I_NW,I_LAT+2*(i-1)), whole(I_NE,I_LAT+2*(i-1)), & ! LAT: NW, NE
                                                      whole(I_SW,I_LAT+2*(i-1)), whole(I_SE,I_LAT+2*(i-1))    ! LAT: SW, SE
             if ( ierr /= 0 ) exit
          enddo
          close(fid)
       endif

       do i = 1, PRC_nmax ! for online nesting
          REAL_DOMAIN_CATALOGUE(i,I_NW,I_LON) = whole(I_NW,I_LON+2*(i-1)) ! LON: NW
          REAL_DOMAIN_CATALOGUE(i,I_NE,I_LON) = whole(I_NE,I_LON+2*(i-1)) ! LON: NE
          REAL_DOMAIN_CATALOGUE(i,I_SW,I_LON) = whole(I_SW,I_LON+2*(i-1)) ! LON: SW
          REAL_DOMAIN_CATALOGUE(i,I_SE,I_LON) = whole(I_SE,I_LON+2*(i-1)) ! LON: SE
          REAL_DOMAIN_CATALOGUE(i,I_NW,I_LAT) = whole(I_NW,I_LAT+2*(i-1)) ! LAT: NW
          REAL_DOMAIN_CATALOGUE(i,I_NE,I_LAT) = whole(I_NE,I_LAT+2*(i-1)) ! LAT: NE
          REAL_DOMAIN_CATALOGUE(i,I_SW,I_LAT) = whole(I_SW,I_LAT+2*(i-1)) ! LAT: SW
          REAL_DOMAIN_CATALOGUE(i,I_SE,I_LAT) = whole(I_SE,I_LAT+2*(i-1)) ! LAT: SE
       enddo
    endif

    call COMM_bcast( REAL_DOMAIN_CATALOGUE, PRC_nmax, 4, 2 )

    return
  end subroutine REAL_calc_latlon

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine REAL_calc_Z
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_grid, only: &
       GRID_CZ,  &
       GRID_FZ,  &
       GRID_CDX, &
       GRID_CDY
    use scale_topography, only: &
       Zsfc => TOPO_Zsfc
    implicit none

    real(RP) :: Htop
    real(RP) :: DFZ

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    Htop = GRID_FZ(KE) - GRID_FZ(KS-1)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       REAL_CZ(k,i,j) = ( Htop - Zsfc(i,j) ) / Htop * GRID_CZ(k) + Zsfc(i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 0, KA
       REAL_FZ(k,i,j) = ( Htop - Zsfc(i,j) ) / Htop * GRID_FZ(k) + Zsfc(i,j)
    enddo
    enddo
    enddo

    REAL_Z1(:,:) = REAL_CZ(KS,:,:) - Zsfc(:,:)

    REAL_PHI(:,:,:) = GRAV * REAL_CZ(:,:,:)

    REAL_ASPECT_MAX = -1.E+30_RP
    REAL_ASPECT_MIN =  1.E+30_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DFZ = REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j)
       REAL_ASPECT_MAX = max( REAL_ASPECT_MAX, GRID_CDX(i) / DFZ, GRID_CDY(j) / DFZ )
       REAL_ASPECT_MIN = min( REAL_ASPECT_MIN, GRID_CDX(i) / DFZ, GRID_CDY(j) / DFZ )
    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '*** Minimum & maximum aspect ratio'
    if( IO_L ) write(IO_FID_LOG,*) '->(',REAL_ASPECT_MIN,',',REAL_ASPECT_MAX,')'

    return
  end subroutine REAL_calc_Z

  !-----------------------------------------------------------------------------
  !> Calc control area/volume
  subroutine REAL_calc_areavol( &
       MAPF )
    use scale_const, only: &
       RADIUS => CONST_RADIUS
    use scale_grid, only: &
       DZ, &
       DX, &
       DY
    implicit none
    real(RP), intent(in) :: MAPF(IA,JA,2)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    REAL_TOTAREA   = 0.0_RP
    REAL_AREA(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       REAL_AREA(i,j) = DX * DY / ( MAPF(i,j,1) * MAPF(i,j,2) )
       REAL_TOTAREA = REAL_TOTAREA + REAL_AREA(i,j)
    enddo
    enddo

    REAL_TOTVOL     = 0.0_RP
    REAL_VOL(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       REAL_VOL(k,i,j) = ( REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j) ) * REAL_AREA(i,j)
       REAL_TOTVOL = REAL_TOTVOL + REAL_VOL(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine REAL_calc_areavol

end module scale_grid_real
