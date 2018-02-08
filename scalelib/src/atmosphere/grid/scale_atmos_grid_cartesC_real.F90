!-------------------------------------------------------------------------------
!> module Atmosphere GRID CartesC Real(real space)
!!
!! @par Description
!!          Grid module for orthogonal curvelinear, terrain-following coordinate
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_atmos_grid_cartesC_real
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_GRID_CARTESC_REAL_setup
  public :: ATMOS_GRID_CARTESC_REAL_update_Z
  public :: ATMOS_GRID_CARTESC_REAL_calc_areavol

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LON(:,:)       !< longitude [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LAT(:,:)       !< latitude  [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_CZ (:,:,:)     !< geopotential height [m] (cell center)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_FZ (:,:,:)     !< geopotential height [m] (cell face  )

  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON  !< position of base point in real world [rad,0-2pi]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT  !< position of base point in real world [rad,-pi,pi]

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LONX (:,:)     !< longitude at staggered point (uy) [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LONY (:,:)     !< longitude at staggered point (xv) [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LONXY(:,:)     !< longitude at staggered point (uv) [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LATX (:,:)     !< latitude  at staggered point (uy) [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LATY (:,:)     !< latitude  at staggered point (xv) [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LATXY(:,:)     !< latitude  at staggered point (uv) [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_DLON (:,:)     !< delta longitude
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_DLAT (:,:)     !< delta latitude

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_Z1  (:,:)      !< Height of the lowermost grid from surface (cell center) [m]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_ASPECT_MAX     !< maximum aspect ratio of the grid cell
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_ASPECT_MIN     !< minimum aspect ratio of the grid cell

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_PHI (:,:,:)    !< geopotential [m2/s2] (cell center)

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREA(:,:)      !< horizontal area [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_VOL (:,:,:)    !< control volume  [m3]

  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTAREA        !< total area   (local) [m2]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTVOL         !< total volume (local) [m3]

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(:,:,:) !< domain latlon catalogue [rad]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_GRID_CARTESC_REAL_calc_latlon
  private :: ATMOS_GRID_CARTESC_REAL_calc_Z

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_GRID_CARTESC_REAL_setup
    use scale_process, only: &
       PRC_nprocs,  &
       PRC_MPIstop
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, &
       ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
       ATMOS_GRID_CARTESC_CZ,              &
       ATMOS_GRID_CARTESC_FZ
    use scale_topography, only: &
       TOPO_exist
    use scale_mapproj, only: &
       MPRJ_setup
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates
    use scale_interp_vert, only: &
       INTERP_VERT_setcoef
    implicit none

    character(len=H_LONG) :: DOMAIN_CATALOGUE_FNAME  = 'latlon_domain_catalogue.txt' !< metadata files for lat-lon domain for all processes
    logical               :: DOMAIN_CATALOGUE_OUTPUT = .false.

    NAMELIST / PARAM_DOMAIN_CATALOGUE / &
       DOMAIN_CATALOGUE_FNAME,  &
       DOMAIN_CATALOGUE_OUTPUT

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID_REAL] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_DOMAIN_CATALOGUE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_DOMAIN_CATALOGUE. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_DOMAIN_CATALOGUE)

    allocate( ATMOS_GRID_CARTESC_REAL_LON  (  IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LAT  (  IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LONX (0:IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LONY (  IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LONXY(0:IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LATX (0:IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LATY (  IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LATXY(0:IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_DLON (  IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_DLAT (  IA,  JA) )

    allocate( ATMOS_GRID_CARTESC_REAL_CZ (  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_FZ (0:KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_Z1 (     IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_PHI(  KA,IA,JA) )

    allocate( ATMOS_GRID_CARTESC_REAL_AREA(   IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_VOL (KA,IA,JA) )

    allocate( ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(PRC_nprocs,4,2) )

    ! setup map projection
    call MPRJ_setup( ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y )

    ! calc longitude & latitude
    call ATMOS_GRID_CARTESC_REAL_calc_latlon( DOMAIN_CATALOGUE_FNAME, DOMAIN_CATALOGUE_OUTPUT )

    ! calc real height
    call ATMOS_GRID_CARTESC_REAL_calc_Z

    ! set latlon and z to fileio module
    call FILE_CARTESC_set_coordinates( ATMOS_GRID_CARTESC_REAL_LON, ATMOS_GRID_CARTESC_REAL_LONX, ATMOS_GRID_CARTESC_REAL_LONY, ATMOS_GRID_CARTESC_REAL_LONXY, & ! [IN]
                                       ATMOS_GRID_CARTESC_REAL_LAT, ATMOS_GRID_CARTESC_REAL_LATX, ATMOS_GRID_CARTESC_REAL_LATY, ATMOS_GRID_CARTESC_REAL_LATXY, & ! [IN]
                                       ATMOS_GRID_CARTESC_REAL_CZ,  ATMOS_GRID_CARTESC_REAL_FZ                           ) ! [IN]

    call INTERP_VERT_setcoef( KA, KS, KE,     & ! [IN]
                              IA, 1,  IA,     & ! [IN]
                              JA, 1,  JA,     & ! [IN]
                              TOPO_exist,     & ! [IN]
                              ATMOS_GRID_CARTESC_CZ(:),     & ! [IN]
                              ATMOS_GRID_CARTESC_FZ(:),     & ! [IN]
                              ATMOS_GRID_CARTESC_REAL_CZ(:,:,:), & ! [IN]
                              ATMOS_GRID_CARTESC_REAL_FZ(:,:,:)  ) ! [IN]

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_setup

  !-----------------------------------------------------------------------------
  !> Re-setup with updated topography
  subroutine ATMOS_GRID_CARTESC_REAL_update_Z
    use scale_process, only: &
       PRC_MPIstop
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates
    implicit none
    !---------------------------------------------------------------------------

    ! calc real height
    call ATMOS_GRID_CARTESC_REAL_calc_Z

    ! set latlon and z to fileio module
    call FILE_CARTESC_set_coordinates( ATMOS_GRID_CARTESC_REAL_LON, ATMOS_GRID_CARTESC_REAL_LONX, ATMOS_GRID_CARTESC_REAL_LONY, ATMOS_GRID_CARTESC_REAL_LONXY, & ! [IN]
                                       ATMOS_GRID_CARTESC_REAL_LAT, ATMOS_GRID_CARTESC_REAL_LATX, ATMOS_GRID_CARTESC_REAL_LATY, ATMOS_GRID_CARTESC_REAL_LATXY, & ! [IN]
                                       ATMOS_GRID_CARTESC_REAL_CZ,  ATMOS_GRID_CARTESC_REAL_FZ                           ) ! [IN]

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_update_Z

  !-----------------------------------------------------------------------------
  !> Calc longitude & latitude
  subroutine ATMOS_GRID_CARTESC_REAL_calc_latlon( &
       catalogue_fname, &
       catalogue_output )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_nprocs,  &
       PRC_IsMaster
    use scale_const, only: &
       PI  => CONST_PI, &
       D2R => CONST_D2R
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CX, &
       ATMOS_GRID_CARTESC_CY, &
       ATMOS_GRID_CARTESC_FX, &
       ATMOS_GRID_CARTESC_FY
    use scale_comm, only: &
       COMM_gather,  &
       COMM_bcast
    use scale_mapproj, only: &
       MPRJ_basepoint_lon, &
       MPRJ_basepoint_lat, &
       MPRJ_xy2lonlat
    implicit none

    character(len=*), intent(in) :: catalogue_fname  !< metadata files for lat-lon domain for all processes
    logical,          intent(in) :: catalogue_output

    integer, parameter :: I_LON  = 1
    integer, parameter :: I_LAT  = 2
    integer, parameter :: I_NW   = 1
    integer, parameter :: I_NE   = 2
    integer, parameter :: I_SW   = 3
    integer, parameter :: I_SE   = 4

    real(RP) :: mine (4,2)            !< send    buffer of lon-lat [deg]
    real(RP) :: whole(4,2*PRC_nprocs) !< recieve buffer of lon-lat [deg]

    integer  :: i, j
    integer  :: fid, ierr
    !---------------------------------------------------------------------------

    ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON = MPRJ_basepoint_lon * D2R
    ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT = MPRJ_basepoint_lat * D2R

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Base position in the global domain (lat,lon)'
    if( IO_L ) write(IO_FID_LOG,*) '*** ->(',ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON/D2R,',',ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT/D2R,')'

    do j = 1, JA
    do i = 1, IA
       call MPRJ_xy2lonlat( ATMOS_GRID_CARTESC_CX(i), ATMOS_GRID_CARTESC_CY(j), ATMOS_GRID_CARTESC_REAL_LON  (i,j), ATMOS_GRID_CARTESC_REAL_LAT  (i,j) )
    enddo
    enddo

    do j = 1, JA
    do i = 0, IA
       call MPRJ_xy2lonlat( ATMOS_GRID_CARTESC_FX(i), ATMOS_GRID_CARTESC_CY(j), ATMOS_GRID_CARTESC_REAL_LONX (i,j), ATMOS_GRID_CARTESC_REAL_LATX (i,j) )
    enddo
    enddo

    do j = 0, JA
    do i = 1, IA
       call MPRJ_xy2lonlat( ATMOS_GRID_CARTESC_CX(i), ATMOS_GRID_CARTESC_FY(j), ATMOS_GRID_CARTESC_REAL_LONY (i,j), ATMOS_GRID_CARTESC_REAL_LATY (i,j) )
    enddo
    enddo

    do j = 0, JA
    do i = 0, IA
       call MPRJ_xy2lonlat( ATMOS_GRID_CARTESC_FX(i), ATMOS_GRID_CARTESC_FY(j), ATMOS_GRID_CARTESC_REAL_LONXY(i,j), ATMOS_GRID_CARTESC_REAL_LATXY(i,j) )
    enddo
    enddo

    ATMOS_GRID_CARTESC_REAL_DLON(:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_DLAT(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       ATMOS_GRID_CARTESC_REAL_DLON(i,j) = ATMOS_GRID_CARTESC_REAL_LONX(i,j) - ATMOS_GRID_CARTESC_REAL_LONX(i-1,j)
       ATMOS_GRID_CARTESC_REAL_DLAT(i,j) = ATMOS_GRID_CARTESC_REAL_LATY(i,j) - ATMOS_GRID_CARTESC_REAL_LATY(i,j-1)
       if( ATMOS_GRID_CARTESC_REAL_DLON(i,j) < 0.0_RP ) ATMOS_GRID_CARTESC_REAL_DLON(i,j) = ATMOS_GRID_CARTESC_REAL_DLON(i,j) + 2.0_RP*PI

       if (      ATMOS_GRID_CARTESC_REAL_DLON(i,j) == 0.0_RP &
            .OR. ATMOS_GRID_CARTESC_REAL_DLAT(i,j) == 0.0_RP ) then
          write(*,*) 'xxx Invalid grid distance in lat-lon! i,j=', i,j
          write(*,*) 'xxx Lon(i-1),Lon(i),dlon=', ATMOS_GRID_CARTESC_REAL_LONX(i-1,j)/D2R,ATMOS_GRID_CARTESC_REAL_LONX(i,j)/D2R,ATMOS_GRID_CARTESC_REAL_DLON(i,j)/D2R
          write(*,*) 'xxx Lat(j-1),Lat(j),dlat=', ATMOS_GRID_CARTESC_REAL_LATY(i,j-1)/D2R,ATMOS_GRID_CARTESC_REAL_LATY(i,j)/D2R,ATMOS_GRID_CARTESC_REAL_DLAT(i,j)/D2R
          call PRC_MPIstop
       endif
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Position on the earth (Local)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.5,A,F9.5,A,A,F10.5,A,F9.5,A)') &
                                '*** NW(',ATMOS_GRID_CARTESC_REAL_LON(IS,JE)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IS,JE)/D2R,')', &
                                 ' - NE(',ATMOS_GRID_CARTESC_REAL_LON(IE,JE)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IE,JE)/D2R,')'

    if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
                                '***              |                          |'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F10.5,A,F9.5,A,A,F10.5,A,F9.5,A)') &
                                '*** SW(',ATMOS_GRID_CARTESC_REAL_LON(IS,JS)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IS,JS)/D2R,')', &
                                 ' - SE(',ATMOS_GRID_CARTESC_REAL_LON(IE,JS)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IE,JS)/D2R,')'

    mine(I_NW,I_LON) = ATMOS_GRID_CARTESC_REAL_LONXY(IS-1,JE  )/D2R
    mine(I_NE,I_LON) = ATMOS_GRID_CARTESC_REAL_LONXY(IE  ,JE  )/D2R
    mine(I_SW,I_LON) = ATMOS_GRID_CARTESC_REAL_LONXY(IS-1,JS-1)/D2R
    mine(I_SE,I_LON) = ATMOS_GRID_CARTESC_REAL_LONXY(IE  ,JS-1)/D2R
    mine(I_NW,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATXY(IS-1,JE  )/D2R
    mine(I_NE,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATXY(IE  ,JE  )/D2R
    mine(I_SW,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATXY(IS-1,JS-1)/D2R
    mine(I_SE,I_LAT) = ATMOS_GRID_CARTESC_REAL_LATXY(IE  ,JS-1)/D2R

    call COMM_gather( whole(:,:), mine(:,:), 4, 2 ) ! everytime do for online nesting

    if ( PRC_IsMaster ) then
       if ( catalogue_output ) then

          fid = IO_get_available_fid()
          open( fid,                            &
                file   = trim(catalogue_fname), &
                form   = 'formatted',           &
                status = 'replace',             &
                iostat = ierr                   )

          if ( ierr /= 0 ) then
             write(*,*) 'xxx [ATMOS_GRID_CARTESC_REAL_calc_latlon] cannot create latlon-catalogue file!'
             call PRC_MPIstop
          endif

          do i = 1, PRC_nprocs ! for offline nesting
             write(fid,'(I8,8F32.24)',iostat=ierr) i, whole(I_NW,I_LON+2*(i-1)), whole(I_NE,I_LON+2*(i-1)), & ! LON: NW, NE
                                                      whole(I_SW,I_LON+2*(i-1)), whole(I_SE,I_LON+2*(i-1)), & ! LON: SW, SE
                                                      whole(I_NW,I_LAT+2*(i-1)), whole(I_NE,I_LAT+2*(i-1)), & ! LAT: NW, NE
                                                      whole(I_SW,I_LAT+2*(i-1)), whole(I_SE,I_LAT+2*(i-1))    ! LAT: SW, SE
             if ( ierr /= 0 ) exit
          enddo

          close(fid)

       endif

       do i = 1, PRC_nprocs ! for online nesting
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_NW,I_LON) = whole(I_NW,I_LON+2*(i-1)) ! LON: NW
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_NE,I_LON) = whole(I_NE,I_LON+2*(i-1)) ! LON: NE
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_SW,I_LON) = whole(I_SW,I_LON+2*(i-1)) ! LON: SW
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_SE,I_LON) = whole(I_SE,I_LON+2*(i-1)) ! LON: SE
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_NW,I_LAT) = whole(I_NW,I_LAT+2*(i-1)) ! LAT: NW
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_NE,I_LAT) = whole(I_NE,I_LAT+2*(i-1)) ! LAT: NE
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_SW,I_LAT) = whole(I_SW,I_LAT+2*(i-1)) ! LAT: SW
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_SE,I_LAT) = whole(I_SE,I_LAT+2*(i-1)) ! LAT: SE
       enddo
    endif

    call COMM_bcast( ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(:,:,:), PRC_nprocs, 4, 2 )

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_calc_latlon

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine ATMOS_GRID_CARTESC_REAL_calc_Z
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ,  &
       ATMOS_GRID_CARTESC_FZ,  &
       ATMOS_GRID_CARTESC_CDX, &
       ATMOS_GRID_CARTESC_CDY
    use scale_topography, only: &
       Zsfc => TOPO_Zsfc
    implicit none

    real(RP) :: Htop
    real(RP) :: DFZ

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    Htop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_GRID_CARTESC_REAL_CZ(k,i,j) = ( Htop - Zsfc(i,j) ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zsfc(i,j)
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = 0, KA
       ATMOS_GRID_CARTESC_REAL_FZ(k,i,j) = ( Htop - Zsfc(i,j) ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zsfc(i,j)
    enddo
    enddo
    enddo

    ATMOS_GRID_CARTESC_REAL_Z1(:,:) = ATMOS_GRID_CARTESC_REAL_CZ(KS,:,:) - Zsfc(:,:)

    ATMOS_GRID_CARTESC_REAL_PHI(:,:,:) = GRAV * ATMOS_GRID_CARTESC_REAL_CZ(:,:,:)

    ATMOS_GRID_CARTESC_REAL_ASPECT_MAX = -1.E+30_RP
    ATMOS_GRID_CARTESC_REAL_ASPECT_MIN =  1.E+30_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DFZ = ATMOS_GRID_CARTESC_REAL_FZ(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(k-1,i,j)
       ATMOS_GRID_CARTESC_REAL_ASPECT_MAX = max( ATMOS_GRID_CARTESC_REAL_ASPECT_MAX, ATMOS_GRID_CARTESC_CDX(i) / DFZ, ATMOS_GRID_CARTESC_CDY(j) / DFZ )
       ATMOS_GRID_CARTESC_REAL_ASPECT_MIN = min( ATMOS_GRID_CARTESC_REAL_ASPECT_MIN, ATMOS_GRID_CARTESC_CDX(i) / DFZ, ATMOS_GRID_CARTESC_CDY(j) / DFZ )
    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Minimum & maximum aspect ratio'
    if( IO_L ) write(IO_FID_LOG,*) '*** ->(',ATMOS_GRID_CARTESC_REAL_ASPECT_MIN,',',ATMOS_GRID_CARTESC_REAL_ASPECT_MAX,')'

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_calc_Z

  !-----------------------------------------------------------------------------
  !> Calc control area/volume
  subroutine ATMOS_GRID_CARTESC_REAL_calc_areavol( &
       MAPF )
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CDX, &
       ATMOS_GRID_CARTESC_CDY
    implicit none

    real(RP), intent(in) :: MAPF(IA,JA,2)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ATMOS_GRID_CARTESC_REAL_AREA(:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTAREA   = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       ATMOS_GRID_CARTESC_REAL_AREA(i,j) = ATMOS_GRID_CARTESC_CDX(i) * ATMOS_GRID_CARTESC_CDY(j) / ( MAPF(i,j,1) * MAPF(i,j,2) )
       ATMOS_GRID_CARTESC_REAL_TOTAREA   = ATMOS_GRID_CARTESC_REAL_TOTAREA + ATMOS_GRID_CARTESC_REAL_AREA(i,j)
    enddo
    enddo

    ATMOS_GRID_CARTESC_REAL_VOL(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTVOL     = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ATMOS_GRID_CARTESC_REAL_VOL(k,i,j) = ( ATMOS_GRID_CARTESC_REAL_FZ(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(k-1,i,j) ) * ATMOS_GRID_CARTESC_REAL_AREA(i,j)
       ATMOS_GRID_CARTESC_REAL_TOTVOL     = ATMOS_GRID_CARTESC_REAL_TOTVOL + ATMOS_GRID_CARTESC_REAL_VOL(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_calc_areavol

end module scale_atmos_grid_cartesC_real
