!-------------------------------------------------------------------------------
!> module FILE I/O (netcdf)
!!
!! @par Description
!!          general file I/O module
!!          frontend interface of netCDF I/O routine
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-20 (H.Yashiro)   [new]
!!
!<
module scale_fileio
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FILEIO_setup
  public :: FILEIO_cleanup
  public :: FILEIO_set_coordinates
  public :: FILEIO_set_axes
  public :: FILEIO_read
  public :: FILEIO_write
  public :: FILEIO_flush

  public :: FILEIO_create
  public :: FILEIO_open
  public :: FILEIO_def_var
  public :: FILEIO_enddef
  public :: FILEIO_read_var
  public :: FILEIO_write_var
  public :: FILEIO_close

  interface FILEIO_read
     module procedure FILEIO_read_1D
     module procedure FILEIO_read_2D
     module procedure FILEIO_read_3D
     module procedure FILEIO_read_4D
  end interface FILEIO_read

  interface FILEIO_read_var
     module procedure FILEIO_read_var_1D
     module procedure FILEIO_read_var_2D
     module procedure FILEIO_read_var_3D
     module procedure FILEIO_read_var_4D
  end interface FILEIO_read_var

  interface FILEIO_write
     module procedure FILEIO_write_1D
     module procedure FILEIO_write_2D
     module procedure FILEIO_write_3D
     module procedure FILEIO_write_3D_t
     module procedure FILEIO_write_4D
  end interface FILEIO_write

  interface FILEIO_write_var
     module procedure FILEIO_write_var_1D
     module procedure FILEIO_write_var_2D
     module procedure FILEIO_write_var_3D
     module procedure FILEIO_write_var_3D_t
     module procedure FILEIO_write_var_4D
  end interface FILEIO_write_var

  !-----------------------------------------------------------------------------
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
  real(RP), private, allocatable :: AXIS_LON (:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LONX(:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LONY(:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LONXY(:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LAT (:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LATX(:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LATY(:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LATXY(:,:) ! [deg]

  integer,  private, parameter :: File_nfile_max = 512   ! number limit of file
                                ! Keep consistency with "FILE_MAX" in gtool_netcdf.c
  logical,  private,      save :: File_axes_written(0:File_nfile_max-1)
                                ! whether axes have been written
  !                             ! fid starts from zero so index should start from zero
  logical,  private,      save :: File_nozcoord(0:File_nfile_max-1)
                                ! whether nozcoord is true or false

  integer,  private,      save :: write_buf_amount = 0  ! sum of write buffer amounts

  ! global star and count for XY and XYT
  integer,  private,      save :: startXY(3),           countXY(3)
  integer,  private,      save :: startWestXY(3),       countWestXY(3)
  integer,  private,      save :: startEastXY(3),       countEastXY(3)
  integer,  private,      save :: startSouthXY(3),      countSouthXY(3)
  integer,  private,      save :: startNorthXY(3),      countNorthXY(3)
  integer,  private,      save :: startSouthWestXY(3),  countSouthWestXY(3)
  integer,  private,      save :: startSouthEastXY(3),  countSouthEastXY(3)
  integer,  private,      save :: startNorthWestXY(3),  countNorthWestXY(3)
  integer,  private,      save :: startNorthEastXY(3),  countNorthEastXY(3)

  ! global star and count for ZX
  integer,  private,      save :: startZX(2),           countZX(2)
  integer,  private,      save :: startWestZX(2),       countWestZX(2)
  integer,  private,      save :: startEastZX(2),       countEastZX(2)

  ! global star and count for ZXY, ZXYT
  integer,  private,      save :: startZXY(4),           countZXY(4)
  integer,  private,      save :: startWestZXY(4),       countWestZXY(4)
  integer,  private,      save :: startEastZXY(4),       countEastZXY(4)
  integer,  private,      save :: startSouthZXY(4),      countSouthZXY(4)
  integer,  private,      save :: startNorthZXY(4),      countNorthZXY(4)
  integer,  private,      save :: startSouthWestZXY(4),  countSouthWestZXY(4)
  integer,  private,      save :: startSouthEastZXY(4),  countSouthEastZXY(4)
  integer,  private,      save :: startNorthWestZXY(4),  countNorthWestZXY(4)
  integer,  private,      save :: startNorthEastZXY(4),  countNorthEastZXY(4)

  ! global star and count for Land
  integer,  private,      save :: startLAND(3),           countLAND(3)
  integer,  private,      save :: startWestLAND(3),       countWestLAND(3)
  integer,  private,      save :: startEastLAND(3),       countEastLAND(3)
  integer,  private,      save :: startSouthLAND(3),      countSouthLAND(3)
  integer,  private,      save :: startNorthLAND(3),      countNorthLAND(3)
  integer,  private,      save :: startSouthWestLAND(3),  countSouthWestLAND(3)
  integer,  private,      save :: startSouthEastLAND(3),  countSouthEastLAND(3)
  integer,  private,      save :: startNorthWestLAND(3),  countNorthWestLAND(3)
  integer,  private,      save :: startNorthEastLAND(3),  countNorthEastLAND(3)

  ! global star and count for Urban
  integer,  private,      save :: startURBAN(3),           countURBAN(3)
  integer,  private,      save :: startWestURBAN(3),       countWestURBAN(3)
  integer,  private,      save :: startEastURBAN(3),       countEastURBAN(3)
  integer,  private,      save :: startSouthURBAN(3),      countSouthURBAN(3)
  integer,  private,      save :: startNorthURBAN(3),      countNorthURBAN(3)
  integer,  private,      save :: startSouthWestURBAN(3),  countSouthWestURBAN(3)
  integer,  private,      save :: startSouthEastURBAN(3),  countSouthEastURBAN(3)
  integer,  private,      save :: startNorthWestURBAN(3),  countNorthWestURBAN(3)
  integer,  private,      save :: startNorthEastURBAN(3),  countNorthEastURBAN(3)

  ! MPI element datatype for restart variables
  integer,  private,      save :: dtype

  ! MPI derived datatype for XY, XYT
  integer,  private,      save :: centerTypeXY
  integer,  private,      save :: westTypeXY
  integer,  private,      save :: eastTypeXY
  integer,  private,      save :: southTypeXY
  integer,  private,      save :: northTypeXY
  integer,  private,      save :: southwestTypeXY
  integer,  private,      save :: southeastTypeXY
  integer,  private,      save :: northwestTypeXY
  integer,  private,      save :: northeastTypeXY

  ! MPI derived datatype for ZX
  integer,  private,      save :: centerTypeZX
  integer,  private,      save :: westTypeZX
  integer,  private,      save :: eastTypeZX

  ! MPI derived datatype for ZXY, ZXYT
  integer,  private,      save :: centerTypeZXY
  integer,  private,      save :: westTypeZXY
  integer,  private,      save :: eastTypeZXY
  integer,  private,      save :: southTypeZXY
  integer,  private,      save :: northTypeZXY
  integer,  private,      save :: southwestTypeZXY
  integer,  private,      save :: southeastTypeZXY
  integer,  private,      save :: northwestTypeZXY
  integer,  private,      save :: northeastTypeZXY

  ! MPI derived datatype for Land
  integer,  private,      save :: centerTypeLAND
  integer,  private,      save :: westTypeLAND
  integer,  private,      save :: eastTypeLAND
  integer,  private,      save :: southTypeLAND
  integer,  private,      save :: northTypeLAND
  integer,  private,      save :: southwestTypeLAND
  integer,  private,      save :: southeastTypeLAND
  integer,  private,      save :: northwestTypeLAND
  integer,  private,      save :: northeastTypeLAND

  ! MPI derived datatype for Urban
  integer,  private,      save :: centerTypeURBAN
  integer,  private,      save :: westTypeURBAN
  integer,  private,      save :: eastTypeURBAN
  integer,  private,      save :: southTypeURBAN
  integer,  private,      save :: northTypeURBAN
  integer,  private,      save :: southwestTypeURBAN
  integer,  private,      save :: southeastTypeURBAN
  integer,  private,      save :: northwestTypeURBAN
  integer,  private,      save :: northeastTypeURBAN

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILEIO_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FIELIO] / Categ[ATMOS-RM IO] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** No namelists.'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** NetCDF header information ***'
    if( IO_L ) write(IO_FID_LOG,*) '*** Data source : ', trim(H_SOURCE)
    if( IO_L ) write(IO_FID_LOG,*) '*** Institute   : ', trim(H_INSTITUTE)

    allocate( AXIS_LON (IMAXB,JMAXB) )
    allocate( AXIS_LONX(IMAXB,JMAXB) )
    allocate( AXIS_LONY(IMAXB,JMAXB) )
    allocate( AXIS_LONXY(IMAXB,JMAXB) )
    allocate( AXIS_LAT (IMAXB,JMAXB) )
    allocate( AXIS_LATX(IMAXB,JMAXB) )
    allocate( AXIS_LATY(IMAXB,JMAXB) )
    allocate( AXIS_LATXY(IMAXB,JMAXB) )

    if ( IO_PNETCDF ) call Construct_Derived_Datatype

    return
  end subroutine FILEIO_setup

  !-----------------------------------------------------------------------------
  !> deallocate buffers
  subroutine FILEIO_cleanup
    implicit none
    !---------------------------------------------------------------------------

    deallocate( AXIS_LON   )
    deallocate( AXIS_LONX  )
    deallocate( AXIS_LONY  )
    deallocate( AXIS_LONXY )
    deallocate( AXIS_LAT   )
    deallocate( AXIS_LATX  )
    deallocate( AXIS_LATY  )
    deallocate( AXIS_LATXY )

    call Free_Derived_Datatype

  end subroutine FILEIO_cleanup

  !-----------------------------------------------------------------------------
  !> set latlon and z
  subroutine FILEIO_set_coordinates( &
       LON,  &
       LONX, &
       LONY, &
       LONXY, &
       LAT,  &
       LATX, &
       LATY, &
       LATXY, &
       CZ,   &
       FZ    )
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    real(RP), intent(in) :: LON (IA,JA)
    real(RP), intent(in) :: LONX(IA,JA)
    real(RP), intent(in) :: LONY(IA,JA)
    real(RP), intent(in) :: LONXY(IA,JA)
    real(RP), intent(in) :: LAT (IA,JA)
    real(RP), intent(in) :: LATX(IA,JA)
    real(RP), intent(in) :: LATY(IA,JA)
    real(RP), intent(in) :: LATXY(IA,JA)
    real(RP), intent(in) :: CZ  (  KA,IA,JA)
    real(RP), intent(in) :: FZ  (0:KA,IA,JA)
    !---------------------------------------------------------------------------

    AXIS_LON (:,:) = LON (ISB:IEB,JSB:JEB) / D2R
    AXIS_LONX(:,:) = LONX(ISB:IEB,JSB:JEB) / D2R
    AXIS_LONY(:,:) = LONY(ISB:IEB,JSB:JEB) / D2R
    AXIS_LONXY(:,:) = LONXY(ISB:IEB,JSB:JEB) / D2R
    AXIS_LAT (:,:) = LAT (ISB:IEB,JSB:JEB) / D2R
    AXIS_LATX(:,:) = LATX(ISB:IEB,JSB:JEB) / D2R
    AXIS_LATY(:,:) = LATY(ISB:IEB,JSB:JEB) / D2R
    AXIS_LATXY(:,:) = LATXY(ISB:IEB,JSB:JEB) / D2R

    return
  end subroutine FILEIO_set_coordinates

  !-----------------------------------------------------------------------------
  !> write axis to the file
  subroutine FILEIO_set_axes( &
       fid,   &
       dtype, &
       xy     )
    use gtool_file, only: &
       FilePutAxis,  &
       FileSetTAttr, &
       FilePutAssociatedCoordinates
    use scale_grid, only: &
       GRID_CZ,    &
       GRID_CX,    &
       GRID_CY,    &
       GRID_FZ,    &
       GRID_FX,    &
       GRID_FY,    &
       GRID_CDZ,   &
       GRID_CDX,   &
       GRID_CDY,   &
       GRID_FDZ,   &
       GRID_FDX,   &
       GRID_FDY,   &
       GRID_CBFZ,  &
       GRID_CBFX,  &
       GRID_CBFY,  &
       GRID_FBFZ,  &
       GRID_FBFX,  &
       GRID_FBFY,  &
       GRID_CXG,   &
       GRID_CYG,   &
       GRID_FXG,   &
       GRID_FYG,   &
       GRID_CBFXG, &
       GRID_CBFYG, &
       GRID_FBFXG, &
       GRID_FBFYG
    use scale_land_grid, only: &
       GRID_LCZ, &
       GRID_LFZ
    use scale_urban_grid, only: &
       GRID_UCZ, &
       GRID_UFZ
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: dtype
    logical, intent(in), optional :: xy

    character(len=2) :: AXIS_name(2)
    logical :: xy_
    !---------------------------------------------------------------------------

    if ( present(xy) ) then
       xy_ = xy
    else
       xy_ = .false.
    end if

    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'z',   'Z',               'm', 'z',   dtype, GRID_CZ(KS:KE) )
    end if
    call FilePutAxis( fid, 'x',   'X',               'm', 'x',   dtype, GRID_CX(ISB:IEB) )
    call FilePutAxis( fid, 'y',   'Y',               'm', 'y',   dtype, GRID_CY(JSB:JEB) )
    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'zh',  'Z (half level)',  'm', 'zh',  dtype, GRID_FZ(KS:KE) )
    end if
    call FilePutAxis( fid, 'xh',  'X (half level)',  'm', 'xh',  dtype, GRID_FX(ISB:IEB) )
    call FilePutAxis( fid, 'yh',  'Y (half level)',  'm', 'yh',  dtype, GRID_FY(JSB:JEB) )

    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'lz',  'LZ',              'm', 'lz',  dtype, GRID_LCZ(LKS:LKE) )
       call FilePutAxis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, GRID_LFZ(LKS:LKE) )
       call FilePutAxis( fid, 'uz',  'UZ',              'm', 'uz',  dtype, GRID_UCZ(UKS:UKE) )
       call FilePutAxis( fid, 'uzh', 'UZ (half level)', 'm', 'uzh', dtype, GRID_UFZ(UKS:UKE) )
    end if

    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'CZ',  'Atmos Grid Center Position Z', 'm', 'CZ',  dtype, GRID_CZ )
    end if
    call FilePutAxis( fid, 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  dtype, GRID_CX )
    call FilePutAxis( fid, 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  dtype, GRID_CY )
    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'FZ',  'Atmos Grid Face Position Z',   'm', 'FZ',  dtype, GRID_FZ )
    end if
    call FilePutAxis( fid, 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  dtype, GRID_FX )
    call FilePutAxis( fid, 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  dtype, GRID_FY )

    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'CDZ', 'Grid Cell length Z', 'm', 'CZ',  dtype, GRID_CDZ )
    end if
    call FilePutAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, GRID_CDX )
    call FilePutAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, GRID_CDY )
    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'FDZ', 'Grid distance Z',    'm', 'FDZ', dtype, GRID_FDZ )
    end if
    call FilePutAxis( fid, 'FDX', 'Grid distance X',    'm', 'FDX', dtype, GRID_FDX )
    call FilePutAxis( fid, 'FDY', 'Grid distance Y',    'm', 'FDY', dtype, GRID_FDY )

    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', dtype, GRID_LCZ )
       call FilePutAxis( fid, 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', dtype, GRID_LFZ )
       call FilePutAxis( fid, 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', dtype, GRID_LCZ )

       call FilePutAxis( fid, 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', dtype, GRID_UCZ )
       call FilePutAxis( fid, 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', dtype, GRID_UFZ )
       call FilePutAxis( fid, 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', dtype, GRID_UCZ )
    end if

    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'CBFZ', 'Boundary factor Center Z', '1', 'CZ', dtype, GRID_CBFZ )
    end if
    call FilePutAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, GRID_CBFX )
    call FilePutAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, GRID_CBFY )
    if ( .NOT. xy_ ) then
       call FilePutAxis( fid, 'FBFZ', 'Boundary factor Face Z',   '1', 'CZ', dtype, GRID_FBFZ )
    end if
    call FilePutAxis( fid, 'FBFX', 'Boundary factor Face X',   '1', 'CX', dtype, GRID_FBFX )
    call FilePutAxis( fid, 'FBFY', 'Boundary factor Face Y',   '1', 'CY', dtype, GRID_FBFY )

    call FilePutAxis( fid, 'CXG', 'Grid Center Position X (global)', 'm', 'CXG', dtype, GRID_CXG )
    call FilePutAxis( fid, 'CYG', 'Grid Center Position Y (global)', 'm', 'CYG', dtype, GRID_CYG )
    call FilePutAxis( fid, 'FXG', 'Grid Face Position X (global)',   'm', 'FXG', dtype, GRID_FXG )
    call FilePutAxis( fid, 'FYG', 'Grid Face Position Y (global)',   'm', 'FYG', dtype, GRID_FYG )

    call FilePutAxis( fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', dtype, GRID_CBFXG )
    call FilePutAxis( fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', dtype, GRID_CBFYG )
    call FilePutAxis( fid, 'FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', dtype, GRID_FBFXG )
    call FilePutAxis( fid, 'FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', dtype, GRID_FBFYG )

    ! associate coordinates
    AXIS_name = (/'x ','y '/)
    call FilePutAssociatedCoordinates( fid, 'lon' , 'longitude'             ,            &
                                       'degrees_east' , AXIS_name, dtype, AXIS_LON (:,:) )
    AXIS_name = (/'xh','y '/)
    call FilePutAssociatedCoordinates( fid, 'lon_uy', 'longitude (half level uy)',            &
                                       'degrees_east' , AXIS_name, dtype, AXIS_LONX(:,:) )
    AXIS_name = (/'x ','yh'/)
    call FilePutAssociatedCoordinates( fid, 'lon_xv', 'longitude (half level xv)',            &
                                       'degrees_east' , AXIS_name, dtype, AXIS_LONY(:,:) )
    AXIS_name = (/'xh','yh'/)
    call FilePutAssociatedCoordinates( fid, 'lon_uv', 'longitude (half level uv)',            &
                                       'degrees_east' , AXIS_name, dtype, AXIS_LONXY(:,:) )
    AXIS_name = (/'x ','y '/)
    call FilePutAssociatedCoordinates( fid, 'lat' , 'latitude'              ,            &
                                       'degrees_north', AXIS_name, dtype, AXIS_LAT (:,:) )
    AXIS_name = (/'xh','y '/)
    call FilePutAssociatedCoordinates( fid, 'lat_uy', 'latitude (half level uy)' ,            &
                                       'degrees_north', AXIS_name, dtype, AXIS_LATX(:,:) )
    AXIS_name = (/'x ','yh'/)
    call FilePutAssociatedCoordinates( fid, 'lat_xv', 'latitude (half level xv)' ,            &
                                       'degrees_north', AXIS_name, dtype, AXIS_LATY(:,:) )
    AXIS_name = (/'xh','yh'/)
    call FilePutAssociatedCoordinates( fid, 'lat_uv', 'latitude (half level uv)' ,            &
                                       'degrees_north', AXIS_name, dtype, AXIS_LATXY(:,:) )

    ! attributes
    if ( .NOT. xy_ ) then
       call FileSetTAttr( fid, 'lz',  'positive', 'down' )
       call FileSetTAttr( fid, 'lzh', 'positive', 'down' )
       call FileSetTAttr( fid, 'uz',  'positive', 'down' )
       call FileSetTAttr( fid, 'uzh', 'positive', 'down' )
       call FileSetTAttr( fid, 'LCZ', 'positive', 'down' )
       call FileSetTAttr( fid, 'LFZ', 'positive', 'down' )
       call FileSetTAttr( fid, 'UCZ', 'positive', 'down' )
       call FileSetTAttr( fid, 'UFZ', 'positive', 'down' )
    end if

    return
  end subroutine FILEIO_set_axes

  subroutine Construct_Derived_Datatype
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use MPI
    implicit none

    integer err, order
    integer sizes(3), subsizes(3), sub_off(3)
    integer XSG, YSG
    integer XS, XE, YS, YE
    !---------------------------------------------------------------------------
    order = MPI_ORDER_FORTRAN

    centerTypeXY        = MPI_DATATYPE_NULL
    westTypeXY          = MPI_DATATYPE_NULL
    eastTypeXY          = MPI_DATATYPE_NULL
    southTypeXY         = MPI_DATATYPE_NULL
    northTypeXY         = MPI_DATATYPE_NULL
    southwestTypeXY     = MPI_DATATYPE_NULL
    southeastTypeXY     = MPI_DATATYPE_NULL
    northwestTypeXY     = MPI_DATATYPE_NULL
    northeastTypeXY     = MPI_DATATYPE_NULL
    centerTypeZX        = MPI_DATATYPE_NULL
    westTypeZX          = MPI_DATATYPE_NULL
    eastTypeZX          = MPI_DATATYPE_NULL

    centerTypeZXY       = MPI_DATATYPE_NULL
    westTypeZXY         = MPI_DATATYPE_NULL
    eastTypeZXY         = MPI_DATATYPE_NULL
    southTypeZXY        = MPI_DATATYPE_NULL
    northTypeZXY        = MPI_DATATYPE_NULL
    southwestTypeZXY    = MPI_DATATYPE_NULL
    southeastTypeZXY    = MPI_DATATYPE_NULL
    northwestTypeZXY    = MPI_DATATYPE_NULL
    northeastTypeZXY    = MPI_DATATYPE_NULL

    centerTypeLAND      = MPI_DATATYPE_NULL
    westTypeLAND        = MPI_DATATYPE_NULL
    eastTypeLAND        = MPI_DATATYPE_NULL
    southTypeLAND       = MPI_DATATYPE_NULL
    northTypeLAND       = MPI_DATATYPE_NULL
    southwestTypeLAND   = MPI_DATATYPE_NULL
    southeastTypeLAND   = MPI_DATATYPE_NULL
    northwestTypeLAND   = MPI_DATATYPE_NULL
    northeastTypeLAND   = MPI_DATATYPE_NULL

    centerTypeURBAN     = MPI_DATATYPE_NULL
    westTypeURBAN       = MPI_DATATYPE_NULL
    eastTypeURBAN       = MPI_DATATYPE_NULL
    southTypeURBAN      = MPI_DATATYPE_NULL
    northTypeURBAN      = MPI_DATATYPE_NULL
    southwestTypeURBAN  = MPI_DATATYPE_NULL
    southeastTypeURBAN  = MPI_DATATYPE_NULL
    northwestTypeURBAN  = MPI_DATATYPE_NULL
    northeastTypeURBAN  = MPI_DATATYPE_NULL

    dtype = MPI_FLOAT
    if ( RP .EQ. 8 ) dtype = MPI_DOUBLE_PRECISION

    if ( .NOT. PRC_PERIODIC_X .AND. .NOT. PRC_PERIODIC_Y) then
       ! for axistype == 'XY'
       startXY(1) = IS_inG - IHALO
       startXY(2) = JS_inG - JHALO
       countXY(1) = IA
       countXY(2) = JA

       ! for axistype == 'ZXY'
       startZXY(1)   = 1
       startZXY(2:3) = startXY(1:2)
       countZXY(1)   = KA
       countZXY(2:3) = countXY(1:2)
       ! construct MPI subarray data type
       sizes(1)    = KA
       sizes(2)    = IA
       sizes(3)    = JA
       subsizes(1) = KMAX
       subsizes(2) = IA
       subsizes(3) = JA
       sub_off(1)  = KS - 1   ! MPI start index starts with 0
       sub_off(2)  = 0
       sub_off(3)  = 0
       call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, centerTypeZXY, err)
       call MPI_Type_commit(centerTypeZXY, err)

       ! for axistype == 'Land'
       startLAND(1)   = 1
       startLAND(2:3) = startXY(1:2)
       countLAND(1)   = LKMAX
       countLAND(2:3) = countXY(1:2)
       ! construct MPI subarray data type
       sizes(1)    = LKMAX
       subsizes(1) = LKMAX
       sub_off(1)  = LKS - 1   ! MPI start index starts with 0
       call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, centerTypeLAND, err)
       call MPI_Type_commit(centerTypeLAND, err)

       ! for axistype == 'URBAN'
       startURBAN(1)   = 1
       startURBAN(2:3) = startXY(1:2)
       countURBAN(1)   = UKMAX
       countURBAN(2:3) = countXY(1:2)
       ! construct MPI subarray data type
       sizes(1)    = UKMAX
       subsizes(1) = UKMAX
       sub_off(1)  = UKS - 1   ! MPI start index starts with 0
       call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, centerTypeURBAN, err)
       call MPI_Type_commit(centerTypeURBAN, err)
    else
       sizes(2)    = IA
       sizes(3)    = JA

       if ( PRC_PERIODIC_X ) then
          XSG = ISGB-IHALO
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) XSG = 1
          XS = 1
          XE = IA
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 )           XS = IS
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) XE = IE
          startXY(1) = XSG
          countXY(1) = IA
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 )           countXY(1) = countXY(1) - IHALO
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) countXY(1) = countXY(1) - IHALO
       end if
       if ( PRC_PERIODIC_Y ) then
          YSG = JSGB-JHALO
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) YSG = 1
          YS = 1
          YE = JA
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 )           YS = JS
          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) YE = JE
          startXY(2) = YSG
          countXY(2) = JA
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 )           countXY(2) = countXY(2) - JHALO
          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) countXY(2) = countXY(2) - JHALO
       end if
       ! construct MPI subarray data type
       subsizes(2) = XE-XS+1
       subsizes(3) = YE-YS+1
       sub_off(2)  = XS - 1   ! MPI start index starts with 0
       sub_off(3)  = YS - 1

       ! for axistype == 'XY'
       call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, centerTypeXY, err)
       call MPI_Type_commit(centerTypeXY, err)

       ! for axistype == 'ZXY'
       sizes(1)    = KA
       subsizes(1) = KMAX
       sub_off(1)  = KS - 1   ! MPI start index starts with 0
       startZXY(1)   = 1
       startZXY(2:3) = startXY(1:2)
       countZXY(1)   = KMAX
       countZXY(2:3) = countXY(1:2)
       call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, centerTypeZXY, err)
       call MPI_Type_commit(centerTypeZXY, err)

       ! for axistype == 'Land'
       sizes(1)    = LKMAX
       subsizes(1) = LKMAX
       sub_off(1)  = LKS - 1   ! MPI start index starts with 0
       startLAND(1)   = 1
       startLAND(2:3) = startXY(1:2)
       countLAND(1)   = LKMAX
       countLAND(2:3) = countXY(1:2)
       call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, centerTypeLAND, err)
       call MPI_Type_commit(centerTypeLAND, err)

       ! for axistype == 'Urban'
       sizes(1)    = UKMAX
       subsizes(1) = UKMAX
       sub_off(1)  = UKS - 1   ! MPI start index starts with 0
       startURBAN(1)   = 1
       startURBAN(2:3) = startXY(1:2)
       countURBAN(1)   = UKMAX
       countURBAN(2:3) = countXY(1:2)
       call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, centerTypeURBAN, err)
       call MPI_Type_commit(centerTypeURBAN, err)

       ! boundary processes need 2nd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! west column processes
             startWestXY(1) = IMAXG - IHALO + 1
             startWestXY(2) = JSGB - JHALO
             countWestXY(1) = IHALO
             countWestXY(2) = JA
             XS = 1
             XE = IHALO
             YS = 1
             YE = JA
             if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south-west corner process
                startWestXY(2) = startWestXY(2) + JHALO
                countWestXY(2) = countWestXY(2) - JHALO
                YS = YS + JHALO
             end if
             if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north-west corner process
                countWestXY(2) = countWestXY(2) - JHALO
                YE = YE - JHALO
             end if
             ! construct MPI subarray data type
             sizes(2)    = IA
             sizes(3)    = JA
             subsizes(2) = XE-XS+1
             subsizes(3) = YE-YS+1
             sub_off(2)  = XS - 1   ! MPI start index starts with 0
             sub_off(3)  = YS - 1
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, westTypeXY, err)
             call MPI_Type_commit(westTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startWestZXY(1)   = 1
             startWestZXY(2:3) = startWestXY(1:2)
             countWestZXY(1)   = KMAX
             countWestZXY(2:3) = countWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, westTypeZXY, err)
             call MPI_Type_commit(westTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startWestLAND(1)   = 1
             startWestLAND(2:3) = startWestXY(1:2)
             countWestLAND(1)   = LKMAX
             countWestLAND(2:3) = countWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, westTypeLAND, err)
             call MPI_Type_commit(westTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startWestURBAN(1)   = 1
             startWestURBAN(2:3) = startWestXY(1:2)
             countWestURBAN(1)   = UKMAX
             countWestURBAN(2:3) = countWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, westTypeURBAN, err)
             call MPI_Type_commit(westTypeURBAN, err)
          end if

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! east column processes
             startEastXY(1) = 1
             startEastXY(2) = JSGB - JHALO
             countEastXY(1) = IHALO
             countEastXY(2) = JA
             XS = IE+1
             XE = IA
             YS = 1
             YE = JA
             if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south-east corner process
                startEastXY(2) = startEastXY(2) + JHALO
                countEastXY(2) = countEastXY(2) - JHALO
                YS = YS + JHALO
             end if
             if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north-east corner process
                countEastXY(2) = countEastXY(2) - JHALO
                YE = YE - JHALO
             end if
             ! construct MPI subarray data type
             sizes(2)    = IA
             sizes(3)    = JA
             subsizes(2) = XE-XS+1
             subsizes(3) = YE-YS+1
             sub_off(2)  = XS - 1   ! MPI start index starts with 0
             sub_off(3)  = YS - 1
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, eastTypeXY, err)
             call MPI_Type_commit(eastTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startEastZXY(1)   = 1
             startEastZXY(2:3) = startEastXY(1:2)
             countEastZXY(1)   = KMAX
             countEastZXY(2:3) = countEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, eastTypeZXY, err)
             call MPI_Type_commit(eastTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startEastLAND(1)   = 1
             startEastLAND(2:3) = startEastXY(1:2)
             countEastLAND(1)   = LKMAX
             countEastLAND(2:3) = countEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, eastTypeLAND, err)
             call MPI_Type_commit(eastTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startEastURBAN(1)   = 1
             startEastURBAN(2:3) = startEastXY(1:2)
             countEastURBAN(1)   = UKMAX
             countEastURBAN(2:3) = countEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, eastTypeURBAN, err)
             call MPI_Type_commit(eastTypeURBAN, err)
          end if
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south row processes
             startSouthXY(1) = ISGB - IHALO
             startSouthXY(2) = JMAXG - JHALO + 1
             countSouthXY(1) = IA
             countSouthXY(2) = JHALO
             XS = 1
             XE = IA
             YS = 1
             YE = JHALO
             ! local buffer = (1:IA, 1:JHALO)
             if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! south-west corner process
                startSouthXY(1) = startSouthXY(1) + IHALO
                countSouthXY(1) = countSouthXY(1) - IHALO
                XS = XS + IHALO
             end if
             if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! south-east corner process
                countSouthXY(1) = countSouthXY(1) - IHALO
                XE = XE - IHALO
             end if
             ! construct MPI subarray data type
             sizes(2)    = IA
             sizes(3)    = JA
             subsizes(2) = XE-XS+1
             subsizes(3) = YE-YS+1
             sub_off(2)  = XS - 1   ! MPI start index starts with 0
             sub_off(3)  = YS - 1
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, southTypeXY, err)
             call MPI_Type_commit(southTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startSouthZXY(1)   = 1
             startSouthZXY(2:3) = startSouthXY(1:2)
             countSouthZXY(1)   = KMAX
             countSouthZXY(2:3) = countSouthXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southTypeZXY, err)
             call MPI_Type_commit(southTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startSouthLAND(1)   = 1
             startSouthLAND(2:3) = startSouthXY(1:2)
             countSouthLAND(1)   = LKMAX
             countSouthLAND(2:3) = countSouthXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southTypeLAND, err)
             call MPI_Type_commit(southTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startSouthURBAN(1)   = 1
             startSouthURBAN(2:3) = startSouthXY(1:2)
             countSouthURBAN(1)   = UKMAX
             countSouthURBAN(2:3) = countSouthXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southTypeURBAN, err)
             call MPI_Type_commit(southTypeURBAN, err)
          end if

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north row processes
             startNorthXY(1) = ISGB - IHALO
             startNorthXY(2) = 1
             countNorthXY(1) = IA
             countNorthXY(2) = JHALO
             XS = 1
             XE = IA
             YS = JE+1
             YE = JA
             if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! north-west corner process
                startNorthXY(1) = startNorthXY(1) + IHALO
                countNorthXY(1) = countNorthXY(1) - IHALO
                XS = XS + IHALO
             end if
             if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! north-east corner process
                countNorthXY(1) = countNorthXY(1) - IHALO
                XE = XE - IHALO
             end if
             ! construct MPI subarray data type
             sizes(2)    = IA
             sizes(3)    = JA
             subsizes(2) = XE-XS+1
             subsizes(3) = YE-YS+1
             sub_off(2)  = XS - 1   ! MPI start index starts with 0
             sub_off(3)  = YS - 1
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, northTypeXY, err)
             call MPI_Type_commit(northTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startNorthZXY(1)   = 1
             startNorthZXY(2:3) = startNorthXY(1:2)
             countNorthZXY(1)   = KMAX
             countNorthZXY(2:3) = countNorthXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northTypeZXY, err)
             call MPI_Type_commit(northTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startNorthLAND(1)   = 1
             startNorthLAND(2:3) = startNorthXY(1:2)
             countNorthLAND(1)   = LKMAX
             countNorthLAND(2:3) = countNorthXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northTypeLAND, err)
             call MPI_Type_commit(northTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startNorthURBAN(1)   = 1
             startNorthURBAN(2:3) = startNorthXY(1:2)
             countNorthURBAN(1)   = UKMAX
             countNorthURBAN(2:3) = countNorthXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northTypeURBAN, err)
             call MPI_Type_commit(northTypeURBAN, err)
          end if
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          sizes(2)    = IA
          sizes(3)    = JA
          subsizes(2) = IHALO
          subsizes(3) = JHALO

          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startSouthWestXY(1) = IMAXG - IHALO + 1
             startSouthWestXY(2) = JMAXG - JHALO + 1
             countSouthWestXY(1) = IHALO
             countSouthWestXY(2) = JHALO
             ! construct MPI subarray data type
             sub_off(2)  = 0   ! MPI start index starts with 0
             sub_off(3)  = 0
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, southwestTypeXY, err)
             call MPI_Type_commit(southwestTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startSouthWestZXY(1)   = 1
             startSouthWestZXY(2:3) = startSouthWestXY(1:2)
             countSouthWestZXY(1)   = KMAX
             countSouthWestZXY(2:3) = countSouthWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southwestTypeZXY, err)
             call MPI_Type_commit(southwestTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startSouthWestLAND(1)   = 1
             startSouthWestLAND(2:3) = startSouthWestXY(1:2)
             countSouthWestLAND(1)   = LKMAX
             countSouthWestLAND(2:3) = countSouthWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southwestTypeLAND, err)
             call MPI_Type_commit(southwestTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startSouthWestURBAN(1)   = 1
             startSouthWestURBAN(2:3) = startSouthWestXY(1:2)
             countSouthWestURBAN(1)   = UKMAX
             countSouthWestURBAN(2:3) = countSouthWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southwestTypeURBAN, err)
             call MPI_Type_commit(southwestTypeURBAN, err)
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             startSouthEastXY(1) = 1
             startSouthEastXY(2) = JMAXG - JHALO + 1
             countSouthEastXY(1) = IHALO
             countSouthEastXY(2) = JHALO
             ! construct MPI subarray data type
             sub_off(2)  = IE   ! MPI start index starts with 0
             sub_off(3)  = 0
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, southeastTypeXY, err)
             call MPI_Type_commit(southeastTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startSouthEastZXY(1)   = 1
             startSouthEastZXY(2:3) = startSouthEastXY(1:2)
             countSouthEastZXY(1)   = KMAX
             countSouthEastZXY(2:3) = countSouthEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southeastTypeZXY, err)
             call MPI_Type_commit(southeastTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startSouthEastLAND(1)   = 1
             startSouthEastLAND(2:3) = startSouthEastXY(1:2)
             countSouthEastLAND(1)   = LKMAX
             countSouthEastLAND(2:3) = countSouthEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southeastTypeLAND, err)
             call MPI_Type_commit(southeastTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startSouthEastURBAN(1)   = 1
             startSouthEastURBAN(2:3) = startSouthEastXY(1:2)
             countSouthEastURBAN(1)   = UKMAX
             countSouthEastURBAN(2:3) = countSouthEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, southeastTypeURBAN, err)
             call MPI_Type_commit(southeastTypeURBAN, err)
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             startNorthWestXY(1) = IMAXG - IHALO + 1
             startNorthWestXY(2) = 1
             countNorthWestXY(1) = IHALO
             countNorthWestXY(2) = JHALO
             ! construct MPI subarray data type
             sub_off(2)  = 0   ! MPI start index starts with 0
             sub_off(3)  = JE
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, northwestTypeXY, err)
             call MPI_Type_commit(northwestTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startNorthWestZXY(1)   = 1
             startNorthWestZXY(2:3) = startNorthWestXY(1:2)
             countNorthWestZXY(1)   = KMAX
             countNorthWestZXY(2:3) = countNorthWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northwestTypeZXY, err)
             call MPI_Type_commit(northwestTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startNorthWestLAND(1)   = 1
             startNorthWestLAND(2:3) = startNorthWestXY(1:2)
             countNorthWestLAND(1)   = LKMAX
             countNorthWestLAND(2:3) = countNorthWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northwestTypeLAND, err)
             call MPI_Type_commit(northwestTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startNorthWestURBAN(1)   = 1
             startNorthWestURBAN(2:3) = startNorthWestXY(1:2)
             countNorthWestURBAN(1)   = UKMAX
             countNorthWestURBAN(2:3) = countNorthWestXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northwestTypeURBAN, err)
             call MPI_Type_commit(northwestTypeURBAN, err)
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             startNorthEastXY(1) = 1
             startNorthEastXY(2) = 1
             countNorthEastXY(1) = IHALO
             countNorthEastXY(2) = JHALO
             ! construct MPI subarray data type
             sub_off(2)  = IE   ! MPI start index starts with 0
             sub_off(3)  = JE
             call MPI_Type_create_subarray(2, sizes(2:3), subsizes(2:3), sub_off(2:3), order, dtype, northeastTypeXY, err)
             call MPI_Type_commit(northeastTypeXY, err)
             ! for axistype == 'ZXY'
             sizes(1)    = KA
             subsizes(1) = KMAX
             sub_off(1)  = KS - 1   ! MPI start index starts with 0
             startNorthEastZXY(1)   = 1
             startNorthEastZXY(2:3) = startNorthEastXY(1:2)
             countNorthEastZXY(1)   = KMAX
             countNorthEastZXY(2:3) = countNorthEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northeastTypeZXY, err)
             call MPI_Type_commit(northeastTypeZXY, err)
             ! for axistype == 'Land'
             sizes(1)    = LKMAX
             subsizes(1) = LKMAX
             sub_off(1)  = LKS - 1   ! MPI start index starts with 0
             startNorthEastLAND(1)   = 1
             startNorthEastLAND(2:3) = startNorthEastXY(1:2)
             countNorthEastLAND(1)   = LKMAX
             countNorthEastLAND(2:3) = countNorthEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northeastTypeLAND, err)
             call MPI_Type_commit(northeastTypeLAND, err)
             ! for axistype == 'Urban'
             sizes(1)    = UKMAX
             subsizes(1) = UKMAX
             sub_off(1)  = UKS - 1   ! MPI start index starts with 0
             startNorthEastURBAN(1)   = 1
             startNorthEastURBAN(2:3) = startNorthEastXY(1:2)
             countNorthEastURBAN(1)   = UKMAX
             countNorthEastURBAN(2:3) = countNorthEastXY(1:2)
             call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, dtype, northeastTypeURBAN, err)
             call MPI_Type_commit(northeastTypeURBAN, err)
          end if
       end if
    end if

    ! for axistype == 'ZX'
    sizes(1)    = KA
    subsizes(1) = KMAX
    sizes(1)    = KA
    subsizes(1) = KMAX
    startZX(1)  = KHALO+1
    countZX(1)  = KHALO
    sub_off(1)  = KHALO   ! MPI start index starts with 0
    if ( .NOT. PRC_PERIODIC_X ) then  ! read data and halos in one call
       startZX(2)  = IS_inG - IHALO
       countZX(2)  = IA
       sizes(2)    = IA
       subsizes(2) = IMAXB
       sub_off(2)  = ISB-1   ! MPI start index starts with 0
       call MPI_Type_create_subarray(2, sizes, subsizes, sub_off, order, dtype, centerTypeZX, err)
       call MPI_Type_commit(centerTypeZX, err)
    else
       XSG = ISGB - IHALO
       if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) XSG = 1
       XS = 1
       XE = IA
       if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 )           XS = IS
       if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) XE = IE
       startZX(2) = XSG
       countZX(2) = IA
       if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 )           countZX(2) = countZX(2) - IHALO
       if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) countZX(2) = countZX(2) - IHALO
       sizes(2)    = IA
       subsizes(2) = XE-XS+1
       sub_off(2)  = XS-1   ! MPI start index starts with 0
       call MPI_Type_create_subarray(2, sizes, subsizes, sub_off, order, dtype, centerTypeZX, err)
       call MPI_Type_commit(centerTypeZX, err)
       if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then
          ! west-most process reads its west-side of halos
          startWestZX(1) = KHALO+1
          countWestZX(1) = KHALO
          startWestZX(2) = IMAXG-IHALO+1
          countWestZX(2) = IHALO
          sizes(2)    = IA
          subsizes(2) = IHALO
          sub_off(2)  = 0   ! MPI start index starts with 0
          call MPI_Type_create_subarray(2, sizes, subsizes, sub_off, order, dtype, westTypeZX, err)
          call MPI_Type_commit(westTypeZX, err)
       elseif ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then
          ! east-most process reads its east-side of halos
          startEastZX(1) = KHALO+1
          countEastZX(1) = KHALO
          startEastZX(2) = 1
          countEastZX(2) = IHALO
          sizes(2)    = IA
          subsizes(2) = IHALO
          sub_off(2)  = IE   ! MPI start index starts with 0
          call MPI_Type_create_subarray(2, sizes, subsizes, sub_off, order, dtype, eastTypeZX, err)
          call MPI_Type_commit(eastTypeZX, err)
       end if
    end if

  end subroutine Construct_Derived_Datatype

  subroutine Free_Derived_Datatype
    use MPI
    implicit none
    integer err

    if ( .NOT. IO_PNETCDF ) return

    if (    centerTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(   centerTypeXY, err)
    if (      westTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     westTypeXY, err)
    if (      eastTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     eastTypeXY, err)
    if (     southTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    southTypeXY, err)
    if (     northTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    northTypeXY, err)
    if ( southwestTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southwestTypeXY, err)
    if ( southeastTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southeastTypeXY, err)
    if ( northwestTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northwestTypeXY, err)
    if ( northeastTypeXY      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northeastTypeXY, err)
    if (    centerTypeZX      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(   centerTypeZX, err)
    if (      westTypeZX      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     westTypeZX, err)
    if (      eastTypeZX      .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     eastTypeZX, err)

    if (    centerTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(   centerTypeZXY, err)
    if (      westTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     westTypeZXY, err)
    if (      eastTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     eastTypeZXY, err)
    if (     southTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    southTypeZXY, err)
    if (     northTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    northTypeZXY, err)
    if ( southwestTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southwestTypeZXY, err)
    if ( southeastTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southeastTypeZXY, err)
    if ( northwestTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northwestTypeZXY, err)
    if ( northeastTypeZXY     .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northeastTypeZXY, err)

    if (    centerTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(   centerTypeLAND, err)
    if (      westTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     westTypeLAND, err)
    if (      eastTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     eastTypeLAND, err)
    if (     southTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    southTypeLAND, err)
    if (     northTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    northTypeLAND, err)
    if ( southwestTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southwestTypeLAND, err)
    if ( southeastTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southeastTypeLAND, err)
    if ( northwestTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northwestTypeLAND, err)
    if ( northeastTypeLAND    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northeastTypeLAND, err)

    if (    centerTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(   centerTypeURBAN, err)
    if (      westTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     westTypeURBAN, err)
    if (      eastTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(     eastTypeURBAN, err)
    if (     southTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    southTypeURBAN, err)
    if (     northTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(    northTypeURBAN, err)
    if ( southwestTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southwestTypeURBAN, err)
    if ( southeastTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(southeastTypeURBAN, err)
    if ( northwestTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northwestTypeURBAN, err)
    if ( northeastTypeURBAN   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(northeastTypeURBAN, err)

  end subroutine Free_Derived_Datatype

  !-----------------------------------------------------------------------------
  !> Read 1D data from file
  subroutine FILEIO_read_1D( &
       var,      &
       basename, &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer               :: dim1_max, dim1_S, dim1_E
    real(RP), allocatable :: var1D(:)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var: ', trim(varname)

    if ( axistype == 'Z' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype == 'X' ) then
       dim1_max = IMAXB
       dim1_S   = ISB
       dim1_E   = IEB
    elseif( axistype == 'Y' ) then
       dim1_max = JMAXB
       dim1_S   = JSB
       dim1_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    allocate( var1D(dim1_max) )

    call FileRead( var1D(:), basename, varname, step, PRC_myrank )
    var(dim1_S:dim1_E) = var1D(1:dim1_max)

    deallocate( var1D )

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_1D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILEIO_read_2D( &
       var,      &
       basename, &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    real(RP), allocatable :: var2D(:,:)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(varname)

    if ( axistype == 'XY' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif( axistype == 'ZX' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    allocate( var2D(dim1_max,dim2_max) )

    call FileRead( var2D(:,:), basename, varname, step, PRC_myrank )
    var(dim1_S:dim1_E,dim2_S:dim2_E) = var2D(1:dim1_max,1:dim2_max)

    deallocate( var2D )

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILEIO_read_3D( &
       var,      &
       basename, &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: axistype   !< axis type (Z/X/Y/T)
    integer,          intent(in)  :: step       !< step number

    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    integer               :: dim3_max, dim3_S, dim3_E
    real(RP), allocatable :: var3D(:,:,:)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(varname)

    if ( axistype == 'ZXY' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else if ( axistype == 'XYT' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim3_max = step
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       dim3_S   = 1
       dim3_E   = step
    else if ( axistype == 'Land' ) then
       dim1_max = LKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else if ( axistype == 'Urban' ) then
       dim1_max = UKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = UKS
       dim1_E   = UKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    allocate( var3D(dim1_max,dim2_max,dim3_max) )

    call FileRead( var3D(:,:,:), basename, varname, step, PRC_myrank )
    var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E) = var3D(1:dim1_max,1:dim2_max,1:dim3_max)

    deallocate( var3D )

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_3D

  !-----------------------------------------------------------------------------
  !> Read 4D data from file
  subroutine FILEIO_read_4D( &
       var,      &
       basename, &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename     !< basename of the file
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: axistype     !< axis type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number

    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    integer               :: dim3_max, dim3_S, dim3_E
    integer               :: dim4_max, dim4_S, dim4_E
    real(RP), allocatable :: var4D(:,:,:,:)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 4D var: ', trim(varname)

    if ( axistype == 'ZXYT' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim4_max = step
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       dim4_S   = 1
       dim4_E   = step
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    allocate( var4D(dim1_max,dim2_max,dim3_max,dim4_max) )

    call FileRead( var4D(:,:,:,:), basename, varname, step, PRC_myrank )
    var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,dim4_S:dim4_E) = var4D(1:dim1_max,1:dim2_max,1:dim3_max,1:dim4_max)

    deallocate( var4D )

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_4D

  !-----------------------------------------------------------------------------
  !> Read 1D data from file
  subroutine FILEIO_read_var_1D( &
       var,      &
       fid,      &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use MPI
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer               :: XS, XE, YS, YE  ! local buffer indices
    integer               :: XSG, YSG        ! global array indices
    integer               :: dtype
    integer               :: start(1)   ! start offset of globale variable
    integer               :: count(1)   ! request length to the globale variable
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var: ', trim(varname)

    dtype = MPI_FLOAT
    if ( RP .EQ. 8 ) dtype = MPI_DOUBLE_PRECISION

    if ( axistype .EQ. 'Z' ) then
       start(1) = 1
       count(1) = KMAX
       call FileRead( var, fid, varname, step, KMAX, dtype, start, count )
    elseif( axistype .EQ. 'X' ) then
       if ( .NOT. PRC_PERIODIC_X ) then  ! read data and halos in one call
          start(1) = IS_inG - IHALO
          count(1) = IA
          call FileRead( var, fid, varname, step, IA, dtype, start, count )
       else
          XSG = ISGB - IHALO
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) XSG = 1
          XS = 1
          XE = IA
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 )           XS = IS
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) XE = IE
          start(1) = XSG
          count(1) = IA
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 )           count(1) = count(1) - IHALO
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) count(1) = count(1) - IHALO
          call FileRead( var(XS:XE), fid, varname, step, count(1), dtype, start, count )
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then
             ! west-most process reads its west-side of halos
             start(1) = IMAXG-IHALO+1
             count(1) = IHALO
             call FileRead( var(1:IHALO), fid, varname, step, IHALO, dtype, start, count )
          elseif ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then
             ! east-most process reads its east-side of halos
             start(1) = 1
             count(1) = IHALO
             call FileRead( var(IE+1:IA), fid, varname, step, IHALO, dtype, start, count )
          end if
       end if
    elseif( axistype .EQ. 'Y' ) then
       if ( .NOT. PRC_PERIODIC_Y ) then  ! read data and halos in one call
          start(1) = JS_inG - JHALO
          count(1) = JA
          call FileRead( var, fid, varname, step, JA, dtype, start, count )
       else
          YSG = JSGB - JHALO
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) YSG = 1
          YS = 1
          YE = JA
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 )           YS = JS
          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) YE = JE
          start(1) = YSG
          count(1) = JA
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 )           count(1) = count(1) - JHALO
          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) count(1) = count(1) - JHALO
          call FileRead( var(YS:YE), fid, varname, step, count(1), dtype, start, count )
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             ! south-most process reads its south-side of halos
             start(1) = JMAXG-JHALO+1
             count(1) = JHALO
             call FileRead( var(1:JHALO), fid, varname, step, JHALO, dtype, start, count )
          elseif ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_Y-1 ) then
             ! north-most process reads its north-side of halos
             start(1) = 1
             count(1) = JHALO
             call FileRead( var(JE+1:JA), fid, varname, step, JHALO, dtype, start, count )
          end if
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_var_1D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILEIO_read_var_2D( &
       var,      &
       fid,      &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop, &
       PRC_LOCAL_COMM_WORLD
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use MPI
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(varname)

    if ( axistype .EQ. 'XY' ) then
       if ( .NOT. PRC_PERIODIC_X .AND. .NOT. PRC_PERIODIC_Y) then
          call FileRead( var, fid, varname, step, IA*JA, dtype, startXY, countXY )
       else
          call FileRead( var, fid, varname, step, 1, centerTypeXY, startXY, countXY )
       end if

       ! boundary processes need 2nd or 3rd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) & ! west column processes read west halos
             call FileRead( var, fid, varname, step, 1, westTypeXY, startWestXY, countWestXY )

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) & ! east column processes read east halos
             call FileRead( var, fid, varname, step, 1, eastTypeXY, startEastXY, countEastXY )
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) & ! south row processes read south halos
             call FileRead( var, fid, varname, step, 1, southTypeXY, startSouthXY, countSouthXY )

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) & ! north row processes read north halos
             call FileRead( var, fid, varname, step, 1, northTypeXY, startNorthXY, countNorthXY )
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southwestTypeXY, startSouthWestXY, countSouthWestXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southeastTypeXY, startSouthEastXY, countSouthEastXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northwestTypeXY, startNorthWestXY, countNorthWestXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northeastTypeXY, startNorthEastXY, countNorthEastXY )
          end if
       end if
    elseif( axistype .EQ. 'ZX' ) then
       if ( .NOT. PRC_PERIODIC_X ) then  ! read data and halos in one call
          call FileRead( var, fid, varname, step, 1, centerTypeZX, startZX, countZX )
       else
          call FileRead( var, fid, varname, step, countZX(2), centerTypeZX, startZX, countZX )
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then
             ! west-most process reads its west-side of halos
             call FileRead( var, fid, varname, step, IHALO, westTypeZX, startWestZX, countWestZX )
          elseif ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then
             ! east-most process reads its east-side of halos
             call FileRead( var, fid, varname, step, IHALO, eastTypeZX, startEastZX, countEastZX )
          end if
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_var_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILEIO_read_var_3D( &
       var,      &
       fid,      &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop, &
       PRC_LOCAL_COMM_WORLD
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable
    integer,          intent(in)  :: fid        !< file ID
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: axistype   !< axis type (Z/X/Y/T)
    integer,          intent(in)  :: step       !< step number

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(varname)

    if ( axistype .EQ. 'ZXY' ) then
       call FileRead( var, fid, varname, step, 1, centerTypeZXY, startZXY, countZXY )

       ! boundary processes need 2nd or 3rd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! west column processes read west halos
             call FileRead( var, fid, varname, step, 1, westTypeZXY, startWestZXY, countWestZXY )
          end if

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! east column processes read east halos
             call FileRead( var, fid, varname, step, 1, eastTypeZXY, startEastZXY, countEastZXY )
          end if
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south row processes read south halos
             call FileRead( var, fid, varname, step, 1, southTypeZXY, startSouthZXY, countSouthZXY )
          end if

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north row processes read north halos
             call FileRead( var, fid, varname, step, 1, northTypeZXY, startNorthZXY, countNorthZXY )
          end if
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southwestTypeZXY, startSouthWestZXY, countSouthWestZXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southeastTypeZXY, startSouthEastZXY, countSouthEastZXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northwestTypeZXY, startNorthWestZXY, countNorthWestZXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northeastTypeZXY, startNorthEastZXY, countNorthEastZXY )
          end if
       end if
    else if ( axistype .EQ. 'XYT' ) then
       startXY(3)          = 1
       startWestXY(3)      = 1
       startEastXY(3)      = 1
       startSouthXY(3)     = 1
       startNorthXY(3)     = 1
       startSouthWestXY(3) = 1
       startSouthEastXY(3) = 1
       startNorthWestXY(3) = 1
       startNorthEastXY(3) = 1

       countXY(3)          = step
       countWestXY(3)      = step
       countEastXY(3)      = step
       countSouthXY(3)     = step
       countNorthXY(3)     = step
       countSouthWestXY(3) = step
       countSouthEastXY(3) = step
       countNorthWestXY(3) = step
       countNorthEastXY(3) = step

       if ( .NOT. PRC_PERIODIC_X .AND. .NOT. PRC_PERIODIC_Y) then
          call FileRead( var, fid, varname, step, step*IA*JA, dtype, startXY, countXY )
       else
          call FileRead( var, fid, varname, step, step, centerTypeXY, startXY, countXY )
       end if

       ! boundary processes need 2nd or 3rd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) & ! west column processes read west halos
             call FileRead( var, fid, varname, step, step, westTypeXY, startWestXY, countWestXY )

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) & ! east column processes read east halos
             call FileRead( var, fid, varname, step, step, eastTypeXY, startEastXY, countEastXY )
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) & ! south row processes read south halos
             call FileRead( var, fid, varname, step, step, southTypeXY, startSouthXY, countSouthXY )

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) & ! north row processes read north halos
             call FileRead( var, fid, varname, step, step, northTypeXY, startNorthXY, countNorthXY )
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, step, southwestTypeXY, startSouthWestXY, countSouthWestXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, step, southeastTypeXY, startSouthEastXY, countSouthEastXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, step, northwestTypeXY, startNorthWestXY, countNorthWestXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, step, northeastTypeXY, startNorthEastXY, countNorthEastXY )
          end if
       end if
    else if ( axistype .EQ. 'Land' ) then
       call FileRead( var, fid, varname, step, 1, centerTypeLAND, startLAND, countLAND )

       ! boundary processes need 2nd or 3rd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! west column processes read west halos
             call FileRead( var, fid, varname, step, 1, westTypeLAND, startWestLAND, countWestLAND )
          end if

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! east column processes read east halos
             call FileRead( var, fid, varname, step, 1, eastTypeLAND, startEastLAND, countEastLAND )
          end if
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south row processes read south halos
             call FileRead( var, fid, varname, step, 1, southTypeLAND, startSouthLAND, countSouthLAND )
          end if

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north row processes read north halos
             call FileRead( var, fid, varname, step, 1, northTypeLAND, startNorthLAND, countNorthLAND )
          end if
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southwestTypeLAND, startSouthWestLAND, countSouthWestLAND )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southeastTypeLAND, startSouthEastLAND, countSouthEastLAND )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northwestTypeLAND, startNorthWestLAND, countNorthWestLAND )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northeastTypeLAND, startNorthEastLAND, countNorthEastLAND )
          end if
       end if
    else if ( axistype .EQ. 'Urban' ) then
       call FileRead( var, fid, varname, step, 1, centerTypeURBAN, startURBAN, countURBAN )

       ! boundary processes need 2nd or 3rd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! west column processes read west halos
             call FileRead( var, fid, varname, step, 1, westTypeURBAN, startWestURBAN, countWestURBAN )
          end if

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! east column processes read east halos
             call FileRead( var, fid, varname, step, 1, eastTypeURBAN, startEastURBAN, countEastURBAN )
          end if
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south row processes read south halos
             call FileRead( var, fid, varname, step, 1, southTypeURBAN, startSouthURBAN, countSouthURBAN )
          end if

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north row processes read north halos
             call FileRead( var, fid, varname, step, 1, northTypeURBAN, startNorthURBAN, countNorthURBAN )
          end if
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southwestTypeURBAN, startSouthWestURBAN, countSouthWestURBAN )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, 1, southeastTypeURBAN, startSouthEastURBAN, countSouthEastURBAN )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northwestTypeURBAN, startNorthWestURBAN, countNorthWestURBAN )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, 1, northeastTypeURBAN, startNorthEastURBAN, countNorthEastURBAN )
          end if
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_var_3D

  !-----------------------------------------------------------------------------
  !> Read 4D data from file
  subroutine FILEIO_read_var_4D( &
       var,      &
       fid,      &
       varname,  &
       axistype, &
       step      )
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop, &
       PRC_LOCAL_COMM_WORLD
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable
    integer,          intent(in)  :: fid          !< file ID
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: axistype     !< axis type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 4D var: ', trim(varname)

    if ( axistype .EQ. 'ZXYT' ) then
       startZXY(4)          = 1
       startWestZXY(4)      = 1
       startEastZXY(4)      = 1
       startSouthZXY(4)     = 1
       startNorthZXY(4)     = 1
       startSouthWestZXY(4) = 1
       startSouthEastZXY(4) = 1
       startNorthWestZXY(4) = 1
       startNorthEastZXY(4) = 1

       countZXY(4)          = step
       countWestZXY(4)      = step
       countEastZXY(4)      = step
       countSouthZXY(4)     = step
       countNorthZXY(4)     = step
       countSouthWestZXY(4) = step
       countSouthEastZXY(4) = step
       countNorthWestZXY(4) = step
       countNorthEastZXY(4) = step

       call FileRead( var, fid, varname, step, step, centerTypeZXY, startZXY, countZXY )

       ! boundary processes need 2nd or 3rd read
       if ( PRC_PERIODIC_X ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 ) then ! west column processes read west halos
             call FileRead( var, fid, varname, step, step, westTypeZXY, startWestZXY, countWestZXY )
          end if

          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 ) then ! east column processes read east halos
             call FileRead( var, fid, varname, step, step, eastTypeZXY, startEastZXY, countEastZXY )
          end if
       end if

       if ( PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then ! south row processes read south halos
             call FileRead( var, fid, varname, step, step, southTypeZXY, startSouthZXY, countSouthZXY )
          end if

          if ( PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then ! north row processes read north halos
             call FileRead( var, fid, varname, step, step, northTypeZXY, startNorthZXY, countNorthZXY )
          end if
       end if

       ! 4 corner processes read corner rectangle
       if ( PRC_PERIODIC_X .OR. PRC_PERIODIC_Y ) then
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! south-west
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, step, southwestTypeZXY, startSouthWestZXY, countSouthWestZXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! south-east
               PRC_2Drank(PRC_myrank,2) .EQ. 0 ) then
             call FileRead( var, fid, varname, step, step, southeastTypeZXY, startSouthEastZXY, countSouthEastZXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. 0 .AND. & ! north-west
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, step, northwestTypeZXY, startNorthWestZXY, countNorthWestZXY )
          end if
          if ( PRC_2Drank(PRC_myrank,1) .EQ. PRC_NUM_X-1 .AND. & ! north-east
               PRC_2Drank(PRC_myrank,2) .EQ. PRC_NUM_Y-1 ) then
             call FileRead( var, fid, varname, step, step, northeastTypeZXY, startNorthEastZXY, countNorthEastZXY )
          end if
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILEIO_read_var_4D

  !-----------------------------------------------------------------------------
  !> Write 1D data to file
  subroutine FILEIO_write_1D( &
       var,      &
       basename, &
       title,    &
       varname,  &
       desc,     &
       unit,     &
       axistype, &
       datatype, &
       date,     &
       subsec,   &
       append    )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FileAddVariable, &
       FileSetGlobalAttribute, &
       FileWrite
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in)  :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer, optional, intent(in) :: date(6)    !< ymdhms of the time
    real(DP),optional, intent(in) :: subsec     !< subsec of the time
    logical, optional, intent(in) :: append   !< switch whether append existing file or not (default=false)

    integer               :: dtype
    character(len=2)      :: dims(1)
    integer               :: dim1_max, dim1_S, dim1_E
    real(RP), allocatable :: var1D(:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    character(len=34) :: tunits
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( datatype == 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
    else
       if ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
          dtype = File_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    call FileCreate( fid,            & ! [OUT]
                     fileexisted,    & ! [OUT]
                     basename,       & ! [IN]
                     title,          & ! [IN]
                     H_SOURCE,       & ! [IN]
                     H_INSTITUTE,    & ! [IN]
                     PRC_masterrank, & ! [IN]
                     PRC_myrank,     & ! [IN]
                     rankidx,        & ! [IN]
                     append = append ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
       if ( present( subsec ) ) then
          call FileSetGlobalAttribute( fid, "time", (/subsec/) )
       else
          call FileSetGlobalAttribute( fid, "time", (/NOWMS/) )
       end if
       if ( present( date ) ) then
          call getCFtunits(tunits, date)
       else
          call getCFtunits(tunits, NOWDATE)
       end if
       call FileSetGlobalAttribute( fid, "time_units", tunits )
    endif

    if ( axistype == 'Z' ) then
       dims = (/'z'/)
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype == 'X' ) then
       dims = (/'x'/)
       dim1_max = IMAXB
       dim1_S   = ISB
       dim1_E   = IEB
    elseif( axistype == 'Y' ) then
       dims = (/'y'/)
       dim1_max = JMAXB
       dim1_S   = JSB
       dim1_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype ) ! [IN]

    allocate( var1D(dim1_max) )

    var1D(1:dim1_max) = var(dim1_S:dim1_E)
    call FileWrite( fid, vid, var1D(:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var1D )

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_1D

  !-----------------------------------------------------------------------------
  !> Write 2D data to file
  subroutine FILEIO_write_2D( &
       var,      &
       basename, &
       title,    &
       varname,  &
       desc,     &
       unit,     &
       axistype, &
       datatype, &
       date,     &
       subsec,   &
       append,   &
       nohalo,   &
       nozcoord  )
    use gtool_file, only: &
       RMISS
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FileAddVariable, &
       FileSetGlobalAttribute, &
       FileWrite
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in)  :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)
    integer, optional, intent(in) :: date(6)    !< ymdhms of the time
    real(DP),optional, intent(in) :: subsec     !< subsec of the time
    logical, optional, intent(in) :: append   !< switch whether append existing file or not (default=false)
    logical, optional, intent(in) :: nohalo   !< switch whether include halo data or not (default=false)
    logical, optional, intent(in) :: nozcoord !< switch whether include zcoordinate or not (default=false)

    real(RP)              :: varhalo( size(var(:,1)), size(var(1,:)) )

    integer               :: dtype
    character(len=2)      :: dims(2)
    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    real(RP), allocatable :: var2D(:,:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    integer :: i, j
    logical :: nohalo_
    character(len=34) :: tunits
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( datatype == 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
    else
       if ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
          dtype = File_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    call FileCreate( fid,            & ! [OUT]
                     fileexisted,    & ! [OUT]
                     basename,       & ! [IN]
                     title,          & ! [IN]
                     H_SOURCE,       & ! [IN]
                     H_INSTITUTE,    & ! [IN]
                     PRC_masterrank, & ! [IN]
                     PRC_myrank,     & ! [IN]
                     rankidx,        & ! [IN]
                     append = append ) ! [IN]

    if ( axistype == 'XY' ) then
       dims = (/'x','y'/)
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif ( axistype == 'UY' ) then
       dims = (/'xh','y '/)
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif ( axistype == 'XV' ) then
       dims = (/'x ','yh'/)
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif ( axistype == 'UV' ) then
       dims = (/'xh','yh'/)
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif( axistype == 'ZX' ) then
       dims = (/'z','x'/)
       dim1_max = KMAX
       dim2_max = IMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype, nozcoord ) ! [IN]
       if ( present( subsec ) ) then
          call FileSetGlobalAttribute( fid, "time", (/subsec/) )
       else
          call FileSetGlobalAttribute( fid, "time", (/NOWMS/) )
       end if
       if ( present( date ) ) then
          call getCFtunits(tunits, date)
       else
          call getCFtunits(tunits, NOWDATE)
       end if
       call FileSetGlobalAttribute( fid, "time_units", tunits )
    endif

    varhalo(:,:) = var(:,:)

    if ( nohalo_ ) then
       ! W halo
       do j = 1, JA
       do i = 1, IS-1
          varhalo(i,j) = RMISS
       end do
       end do
       ! E halo
       do j = 1, JA
       do i = IE+1, IA
          varhalo(i,j) = RMISS
       end do
       end do
       ! S halo
       do j = 1, JS-1
       do i = 1, IA
          varhalo(i,j) = RMISS
       end do
       end do
       ! N halo
       do j = JE+1, JA
       do i = 1, IA
          varhalo(i,j) = RMISS
       end do
       end do
    end if

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype ) ! [IN]

    allocate( var2D(dim1_max,dim2_max) )

    var2D(1:dim1_max,1:dim2_max) = varhalo(dim1_S:dim1_E,dim2_S:dim2_E)
    call FileWrite( fid, vid, var2D(:,:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var2D )

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_2D

  !-----------------------------------------------------------------------------
  !> Write 3D data to file
  subroutine FILEIO_write_3D( &
       var,      &
       basename, &
       title,    &
       varname,  &
       desc,     &
       unit,     &
       axistype, &
       datatype, &
       date,     &
       subsec,   &
       append,   &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FileAddVariable, &
       FileSetGlobalAttribute, &
       FileWrite
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    real(RP),          intent(in)  :: var(:,:,:) !< value of the variable
    character(len=*),  intent(in)  :: basename   !< basename of the file
    character(len=*),  intent(in)  :: title      !< title    of the file
    character(len=*),  intent(in)  :: varname    !< name        of the variable
    character(len=*),  intent(in)  :: desc       !< description of the variable
    character(len=*),  intent(in)  :: unit       !< unit        of the variable
    character(len=*),  intent(in)  :: axistype   !< axis type (Z/X/Y)
    character(len=*),  intent(in)  :: datatype   !< data type (REAL8/REAL4/default)
    integer, optional, intent(in)  :: date(6)    !< ymdhms of the time
    real(DP),optional, intent(in)  :: subsec     !< subsec of the time
    logical, optional, intent(in)  :: append     !< append existing (closed) file?
    logical, optional, intent(in)  :: nohalo     !< include halo data?

    real(RP)         :: varhalo( size(var(:,1,1)), size(var(1,:,1)), size(var(1,1,:)) )

    integer          :: dtype
    character(len=2) :: dims(3)
    integer          :: dim1_max, dim1_S, dim1_E
    integer          :: dim2_max, dim2_S, dim2_E
    integer          :: dim3_max, dim3_S, dim3_E

    real(RP), allocatable :: var3D(:,:,:)

    integer :: rankidx(2)
    logical :: append_sw
    logical :: fileexisted
    integer :: fid, vid
    integer :: i, j, k
    logical :: nohalo_
    character(len=34) :: tunits
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( datatype == 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
    else
       if ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
          dtype = File_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    append_sw = .false.
    if ( present(append) ) then
       append_sw = append
    endif

    call FileCreate( fid,                 & ! [OUT]
                     fileexisted,         & ! [OUT]
                     basename,            & ! [IN]
                     title,               & ! [IN]
                     H_SOURCE,            & ! [IN]
                     H_INSTITUTE,         & ! [IN]
                     PRC_masterrank,      & ! [IN]
                     PRC_myrank,          & ! [IN]
                     rankidx,             & ! [IN]
                     append = append_sw   ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
       if ( present( subsec ) ) then
          call FileSetGlobalAttribute( fid, "time", (/subsec/) )
       else
          call FileSetGlobalAttribute( fid, "time", (/NOWMS/) )
       end if
       if ( present( date ) ) then
          call getCFtunits(tunits, date)
       else
          call getCFtunits(tunits, NOWDATE)
       end if
       call FileSetGlobalAttribute( fid, "time_units", tunits )
    endif

    if ( axistype == 'ZXY' ) then
       dims = (/'z','x','y'/)
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype == 'ZHXY' ) then
       dims = (/'zh','x ','y '/)
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype == 'ZXHY' ) then
       dims = (/'z ','xh','y '/)
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype == 'ZXYH' ) then
       dims = (/'z ','x ','yh'/)
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype == 'Land' ) then
       dims = (/'lz','x ','y '/)
       dim1_max = LKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype == 'Urban' ) then
       dims = (/'uz','x ','y '/)
       dim1_max = UKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = UKS
       dim1_E   = UKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    varhalo(:,:,:) = var(:,:,:)

    if ( nohalo_ ) then
       ! W halo
       do k = 1, dim1_max
       do j = 1, JA
       do i = 1, IS-1
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
       ! E halo
       do k = 1, dim1_max
       do j = 1, JA
       do i = IE+1, IA
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
       ! S halo
       do k = 1, dim1_max
       do j = 1, JS-1
       do i = 1, IA
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
       ! N halo
       do k = 1, dim1_max
       do j = JE+1, JA
       do i = 1, IA
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
    end if

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype ) ! [IN]

    allocate( var3D(dim1_max,dim2_max,dim3_max) )

    var3D(1:dim1_max,1:dim2_max,1:dim3_max) = varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E)
    call FileWrite( fid, vid, var3D(:,:,:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var3D )

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_3D

  !-----------------------------------------------------------------------------
  !> Write 3D data with time dimension to file
  subroutine FILEIO_write_3D_t( &
       var,      &
       basename, &
       title,    &
       varname,  &
       desc,     &
       unit,     &
       axistype, &
       datatype, &
       timeintv, &
       tsince,   &
       append,   &
       timetarg, &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FilePutAxis,     &
       FileAddVariable, &
       FileWrite
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    implicit none

    real(RP),          intent(in)  :: var(:,:,:) !< value of the variable
    character(len=*),  intent(in)  :: basename     !< basename of the file
    character(len=*),  intent(in)  :: title        !< title    of the file
    character(len=*),  intent(in)  :: varname      !< name        of the variable
    character(len=*),  intent(in)  :: desc         !< description of the variable
    character(len=*),  intent(in)  :: unit         !< unit        of the variable
    character(len=*),  intent(in)  :: axistype     !< axis type (X/Y/Time)
    character(len=*),  intent(in)  :: datatype     !< data type (REAL8/REAL4/default)
    real(RP),          intent(in)  :: timeintv     !< time interval [sec]
    integer ,          intent(in)  :: tsince(6)    !< start time
    logical, optional, intent(in)  :: append       !< append existing (closed) file?
    integer, optional, intent(in)  :: timetarg     !< target timestep (optional)
    logical, optional, intent(in)  :: nohalo       !< include halo data?

    real(RP)         :: varhalo( size(var(:,1,1)), size(var(1,:,1)) )

    integer          :: dtype
    character(len=2) :: dims(2)
    integer          :: dim1_max, dim1_S, dim1_E
    integer          :: dim2_max, dim2_S, dim2_E

    real(RP), allocatable :: var2D(:,:)
    real(DP) :: time_interval, nowtime

    character(len=34) :: tunits

    integer :: rankidx(2)
    logical :: append_sw
    logical :: fileexisted
    integer :: fid, vid
    integer :: step
    integer :: i, j, n
    logical :: nohalo_
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    time_interval = timeintv
    step = size(var(ISB,JSB,:))

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( datatype == 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
    else
       if ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
          dtype = File_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    append_sw = .false.
    if ( present(append) ) then
       append_sw = append
    endif

    write(tunits,'(a,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') 'seconds since ', tsince

    call FileCreate( fid,                 & ! [OUT]
                     fileexisted,         & ! [OUT]
                     basename,            & ! [IN]
                     title,               & ! [IN]
                     H_SOURCE,            & ! [IN]
                     H_INSTITUTE,         & ! [IN]
                     PRC_masterrank,      & ! [IN]
                     PRC_myrank,          & ! [IN]
                     rankidx,             & ! [IN]
                     time_units = tunits, & ! [IN]
                     append = append_sw   ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
    endif

    if ( axistype == 'XYT' ) then
       dims = (/'x','y'/)
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype, time_interval ) ! [IN]
    allocate( var2D(dim1_max,dim2_max) )

    if ( present(timetarg) ) then
       varhalo(:,:) = var(:,:,timetarg)

       if ( nohalo_ ) then
          ! W halo
          do j = 1, JA
          do i = 1, IS-1
             varhalo(i,j) = RMISS
          end do
          end do
          ! E halo
          do j = 1, JA
          do i = IE+1, IA
             varhalo(i,j) = RMISS
          end do
          end do
          ! S halo
          do j = 1, JS-1
          do i = 1, IA
             varhalo(i,j) = RMISS
          end do
          end do
          ! N halo
          do j = JE+1, JA
          do i = 1, IA
             varhalo(i,j) = RMISS
          end do
          end do
       end if

       nowtime = (timetarg-1) * time_interval
       var2D(1:dim1_max,1:dim2_max) = varhalo(dim1_S:dim1_E,dim2_S:dim2_E)
       call FileWrite( fid, vid, var2D(:,:), nowtime, nowtime ) ! [IN]
    else
       nowtime = 0.0_DP
       do n = 1, step
          varhalo(:,:) = var(:,:,n)

          if ( nohalo_ ) then
             ! W halo
             do j = 1, JA
             do i = 1, IS-1
                varhalo(i,j) = RMISS
             end do
             end do
             ! E halo
             do j = 1, JA
             do i = IE+1, IA
                varhalo(i,j) = RMISS
             end do
             end do
             ! S halo
             do j = 1, JS-1
             do i = 1, IA
                varhalo(i,j) = RMISS
             end do
             end do
             ! N halo
             do j = JE+1, JA
             do i = 1, IA
                varhalo(i,j) = RMISS
             end do
             end do
          end if

          var2D(1:dim1_max,1:dim2_max) = varhalo(dim1_S:dim1_E,dim2_S:dim2_E)
          call FileWrite( fid, vid, var2D(:,:), nowtime, nowtime ) ! [IN]
          nowtime = nowtime + time_interval
       enddo
    endif

    deallocate( var2D )

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_3D_t

  !-----------------------------------------------------------------------------
  !> Write 4D data to file
  subroutine FILEIO_write_4D( &
       var,      &
       basename, &
       title,    &
       varname,  &
       desc,     &
       unit,     &
       axistype, &
       datatype, &
       timeintv, &
       tsince,   &
       append,   &
       timetarg, &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FilePutAxis,     &
       FileAddVariable, &
       FileWrite
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    implicit none

    real(RP),          intent(in)  :: var(:,:,:,:) !< value of the variable
    character(len=*),  intent(in)  :: basename     !< basename of the file
    character(len=*),  intent(in)  :: title        !< title    of the file
    character(len=*),  intent(in)  :: varname      !< name        of the variable
    character(len=*),  intent(in)  :: desc         !< description of the variable
    character(len=*),  intent(in)  :: unit         !< unit        of the variable
    character(len=*),  intent(in)  :: axistype     !< axis type (Z/X/Y/Time)
    character(len=*),  intent(in)  :: datatype     !< data type (REAL8/REAL4/default)
    real(RP),          intent(in)  :: timeintv     !< time interval [sec]
    integer ,          intent(in)  :: tsince(6)    !< start time
    logical, optional, intent(in)  :: append       !< append existing (closed) file?
    integer, optional, intent(in)  :: timetarg     !< target timestep (optional)
    logical, optional, intent(in)  :: nohalo       !< include halo data?

    real(RP)         :: varhalo( size(var(:,1,1,1)), size(var(1,:,1,1)), size(var(1,1,:,1)) )

    integer          :: dtype
    character(len=2) :: dims(3)
    integer          :: dim1_max, dim1_S, dim1_E
    integer          :: dim2_max, dim2_S, dim2_E
    integer          :: dim3_max, dim3_S, dim3_E

    real(RP), allocatable :: var3D(:,:,:)
    real(DP) :: time_interval, nowtime

    character(len=34) :: tunits

    integer :: rankidx(2)
    logical :: append_sw
    logical :: fileexisted
    integer :: fid, vid
    integer :: step
    integer :: i, j, k, n
    logical :: nohalo_
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    time_interval = timeintv
    step = size(var(KS,ISB,JSB,:))

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( datatype == 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
    else
       if ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
          dtype = File_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    append_sw = .false.
    if ( present(append) ) then
       append_sw = append
    endif

    call getCFtunits(tunits, tsince)

    call FileCreate( fid,                 & ! [OUT]
                     fileexisted,         & ! [OUT]
                     basename,            & ! [IN]
                     title,               & ! [IN]
                     H_SOURCE,            & ! [IN]
                     H_INSTITUTE,         & ! [IN]
                     PRC_masterrank,      & ! [IN]
                     PRC_myrank,          & ! [IN]
                     rankidx,             & ! [IN]
                     time_units = tunits, & ! [IN]
                     append = append_sw   ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
    endif

    if ( axistype == 'ZXYT' ) then
       dims = (/'z','x','y'/)
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype, time_interval ) ! [IN]
    allocate( var3D(dim1_max,dim2_max,dim3_max) )

    if ( present(timetarg) ) then
       varhalo(:,:,:) = var(:,:,:,timetarg)

       if ( nohalo_ ) then
          ! W halo
          do k = 1, dim1_max
          do j = 1, JA
          do i = 1, IS-1
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
          ! E halo
          do k = 1, dim1_max
          do j = 1, JA
          do i = IE+1, IA
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
          ! S halo
          do k = 1, dim1_max
          do j = 1, JS-1
          do i = 1, IA
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
          ! N halo
          do k = 1, dim1_max
          do j = JE+1, JA
          do i = 1, IA
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
       end if

       nowtime = (timetarg-1) * time_interval
       var3D(1:dim1_max,1:dim2_max,1:dim3_max) = varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E)
       call FileWrite( fid, vid, var3D(:,:,:), nowtime, nowtime ) ! [IN]
    else
       nowtime = 0.0_DP
       do n = 1, step
          varhalo(:,:,:) = var(:,:,:,n)

          if ( nohalo_ ) then
             ! W halo
             do k = 1, dim1_max
             do j = 1, JA
             do i = 1, IS-1
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
             ! E halo
             do k = 1, dim1_max
             do j = 1, JA
             do i = IE+1, IA
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
             ! S halo
             do k = 1, dim1_max
             do j = 1, JS-1
             do i = 1, IA
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
             ! N halo
             do k = 1, dim1_max
             do j = JE+1, JA
             do i = 1, IA
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
          end if

          var3D(1:dim1_max,1:dim2_max,1:dim3_max) = varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E)
          call FileWrite( fid, vid, var3D(:,:,:), nowtime, nowtime ) ! [IN]
          nowtime = nowtime + time_interval
       enddo
    endif

    deallocate( var3D )

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_4D

  ! private

  subroutine getCFtunits(tunits, date)
    character(len=34), intent(out) :: tunits
    integer,           intent(in)  :: date(6)

    write(tunits,'(a,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') 'seconds since ', date

    return
  end subroutine getCFtunits

  !-----------------------------------------------------------------------------
  !> open a netCDF file for read
  subroutine FILEIO_open( &
       fid,      &
       basename  )
    use gtool_file_h, only: &
       File_FREAD
    use gtool_file, only: &
       FileOpen
    use scale_process, only: &
       PRC_myrank, &
       PRC_LOCAL_COMM_WORLD
    use MPI, only : MPI_COMM_NULL
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of the file

    integer :: comm
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( IO_PNETCDF ) then  ! user input parameter indicates to do PnetCDF I/O
       comm = PRC_LOCAL_COMM_WORLD
    else
       comm = MPI_COMM_NULL
    end if

    call FileOpen( fid,                & ! [OUT]
                   basename,           & ! [IN]
                   File_FREAD,         & ! [IN]
                   comm = comm,        & ! [IN]
                   myrank = PRC_myrank ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_open

  !-----------------------------------------------------------------------------
  !> Create/open a netCDF file
  subroutine FILEIO_create( &
       fid,      &
       basename, &
       title,    &
       datatype, &
       date,     &
       subsec,   &
       append,   &
       nozcoord  )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FileSetGlobalAttribute
    use scale_process, only: &
       PRC_masterrank, &
       PRC_myrank,     &
       PRC_MPIstop,    &
       PRC_LOCAL_COMM_WORLD
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    use MPI, only : MPI_COMM_NULL
    implicit none

    integer,           intent(out) :: fid      !< file ID
    character(len=*),  intent(in)  :: basename !< basename of the file
    character(len=*),  intent(in)  :: title    !< title    of the file
    character(len=*),  intent(in)  :: datatype !< data type (REAL8/REAL4/default)
    integer, optional, intent(in)  :: date(6)  !< ymdhms of the time
    real(DP),optional, intent(in)  :: subsec   !< subsec of the time
    logical, optional, intent(in)  :: append   !< switch whether append existing file or not (default=false)
    logical, optional, intent(in)  :: nozcoord !< switch whether include zcoordinate or not (default=false)

    integer :: rankidx(2)
    integer :: dtype
    logical :: append_sw
    logical :: fileexisted
    integer :: comm
    character(len=34) :: tunits
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    ! dtype is used to define the data type of axis variables in file
    if ( datatype .EQ. 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype .EQ. 'REAL4' ) then
       dtype = File_REAL4
    else
       if ( RP .EQ. 8 ) then
          dtype = File_REAL8
       elseif( RP .EQ. 4 ) then
          dtype = File_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype)
          call PRC_MPIstop
       endif
    endif

    append_sw = .false.
    if ( present(append) ) then
       append_sw = append
    endif

    ! create a netCDF file if not already existed. Otherwise, open it.
    if ( present(date) ) then
      call getCFtunits(tunits, date)
    else
      tunits = 'seconds'
    endif

    if ( IO_PNETCDF ) then  ! user input parameter indicates to do PnetCDF I/O
       comm = PRC_LOCAL_COMM_WORLD
    else
       comm = MPI_COMM_NULL
    end if

    call FileCreate( fid,                 & ! [OUT]
                     fileexisted,         & ! [OUT]
                     basename,            & ! [IN]
                     title,               & ! [IN]
                     H_SOURCE,            & ! [IN]
                     H_INSTITUTE,         & ! [IN]
                     PRC_masterrank,      & ! [IN]
                     PRC_myrank,          & ! [IN]
                     rankidx,             & ! [IN]
                     time_units = tunits, & ! [IN]
                     append = append_sw,  & ! [IN]
                     comm = comm          ) ! [IN]

    if ( .NOT. fileexisted ) then ! do below only once when file is created
       call FILEIO_def_axes( fid, dtype, xy = nozcoord ) ! [IN]
       File_axes_written(fid) = .FALSE.  ! indicating axes have not been written yet
       if ( present( nozcoord ) ) then
          File_nozcoord(fid) = nozcoord
       else
          File_nozcoord(fid) = .FALSE.
       endif

       if ( present( subsec ) ) then
          call FileSetGlobalAttribute( fid, "time", (/subsec/) )
       else
          call FileSetGlobalAttribute( fid, "time", (/NOWMS/) )
       end if
       if ( present( date ) ) then
          call getCFtunits(tunits, date)
       else
          call getCFtunits(tunits, NOWDATE)
       end if
       call FileSetGlobalAttribute( fid, "time_units", tunits )
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF file define mode
  subroutine FILEIO_enddef( &
       fid)
    use gtool_file, only: &
       FileEndDef, &
       FileAttachBuffer
    implicit none

    integer, intent(in) :: fid  !< file ID

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call FileEndDef( fid ) ! [IN]

    ! If this enddef is called the first time, write axis variables
    if ( .NOT. File_axes_written(fid) ) then
       if ( IO_PNETCDF ) then
          call FILEIO_write_axes_par(fid, File_nozcoord(fid))
          ! Tell PnetCDF library to use a buffer of size write_buf_amount to
          ! aggregate write requests to be post in FILEIO_write_var
          call FileAttachBuffer(fid, write_buf_amount)
       else
          call FILEIO_write_axes(fid, File_nozcoord(fid))
       end if
       File_axes_written(fid) = .TRUE.
    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_enddef

  !-----------------------------------------------------------------------------
  !> Flush all pending requests to a netCDF file (PnetCDF only)
  subroutine FILEIO_flush( &
       fid)
    use gtool_file, only: &
       FileFlush
    implicit none

    integer, intent(in) :: fid  !< file ID

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( IO_PNETCDF ) then
       call FileFlush( fid )        ! flush all pending read/write requests
    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_flush

  !-----------------------------------------------------------------------------
  !> Close a netCDF file
  subroutine FILEIO_close( &
       fid)
    use gtool_file, only: &
       FileClose, &
       FileFlush, &
       FileDetachBuffer
    implicit none

    integer, intent(in) :: fid  !< file ID

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( IO_PNETCDF ) then
       call FileFlush( fid )        ! flush all pending read/write requests
       if ( write_buf_amount .GT. 0 ) then
          call FileDetachBuffer( fid ) ! detach PnetCDF aggregation buffer
          write_buf_amount = 0         ! reset write request amount
       end if
    end if

    call FileClose( fid ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_close

  !-----------------------------------------------------------------------------
  !> define axis variables in the file
  subroutine FILEIO_def_axes( &
       fid,   &
       dtype, &
       xy     )
    use gtool_file, only: &
       FileDefAxis,  &
       FileSetTAttr, &
       FileDefAssociatedCoordinates
    use scale_rm_process, only: &
       PRC_NUM_X, &
       PRC_NUM_Y, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: dtype
    logical, intent(in), optional :: xy

    character(len=2) :: AXIS_name(2)
    logical :: xy_

    integer :: IAG, JAG
    integer :: XG, YG
    !---------------------------------------------------------------------------

    ! array size (global domain)
    IAG = IHALO + IMAX*PRC_NUM_X + IHALO
    JAG = JHALO + JMAX*PRC_NUM_Y + JHALO

    ! global size for x and xh
    if ( PRC_PERIODIC_X ) then
       XG = IMAX*PRC_NUM_X
    else
       XG = IAG
    end if
    ! global size for y and yh
    if ( PRC_PERIODIC_Y ) then
       YG = JMAX*PRC_NUM_Y
    else
       YG = JAG
    end if

    if ( present(xy) ) then
       xy_ = xy
    else
       xy_ = .false.
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'z',   'Z',               'm', 'z',   dtype, KE-KS+1 )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'x',   'X',               'm', 'x',   dtype, IEB-ISB+1 )
       call FileDefAxis( fid, 'y',   'Y',               'm', 'y',   dtype, JEB-JSB+1 )
    else
       call FileDefAxis( fid, 'x',   'X',               'm', 'x',   dtype, XG )
       call FileDefAxis( fid, 'y',   'Y',               'm', 'y',   dtype, YG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'zh',  'Z (half level)',  'm', 'zh',  dtype, KE-KS+1 )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'xh',  'X (half level)',  'm', 'xh',  dtype, IEB-ISB+1 )
       call FileDefAxis( fid, 'yh',  'Y (half level)',  'm', 'yh',  dtype, JEB-JSB+1 )
    else
       call FileDefAxis( fid, 'xh',  'X (half level)',  'm', 'xh',  dtype, XG )
       call FileDefAxis( fid, 'yh',  'Y (half level)',  'm', 'yh',  dtype, YG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'lz',  'LZ',              'm', 'lz',  dtype, LKE-LKS+1 )
       call FileDefAxis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, LKE-LKS+1 )
       call FileDefAxis( fid, 'uz',  'UZ',              'm', 'uz',  dtype, UKE-UKS+1 )
       call FileDefAxis( fid, 'uzh', 'UZ (half level)', 'm', 'uzh', dtype, UKE-UKS+1 )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'CZ',  'Atmos Grid Center Position Z', 'm', 'CZ',  dtype, KA )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  dtype, IA )
       call FileDefAxis( fid, 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  dtype, JA )
    else
       call FileDefAxis( fid, 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  dtype, IAG )
       call FileDefAxis( fid, 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'FZ',  'Atmos Grid Face Position Z',   'm', 'FZ',  dtype, KA+1 )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  dtype, IA+1 )
       call FileDefAxis( fid, 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  dtype, JA+1 )
    else
       call FileDefAxis( fid, 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  dtype, IAG+1 )
       call FileDefAxis( fid, 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  dtype, JAG+1 )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'CDZ', 'Grid Cell length Z', 'm', 'CZ',  dtype, KA )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, IA )
       call FileDefAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, JA )
    else
       call FileDefAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, IAG )
       call FileDefAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'FDZ', 'Grid distance Z',    'm', 'FDZ', dtype, KA-1 )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'FDX', 'Grid distance X',    'm', 'FDX', dtype, IA-1 )
       call FileDefAxis( fid, 'FDY', 'Grid distance Y',    'm', 'FDY', dtype, JA-1 )
    else
       call FileDefAxis( fid, 'FDX', 'Grid distance X',    'm', 'FDX', dtype, IAG-1 )
       call FileDefAxis( fid, 'FDY', 'Grid distance Y',    'm', 'FDY', dtype, JAG-1 )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', dtype, LKE-LKS+1 )
       call FileDefAxis( fid, 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', dtype, LKE-LKS+2 )
       call FileDefAxis( fid, 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', dtype, LKE-LKS+1 )

       call FileDefAxis( fid, 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', dtype, UKE-UKS+1 )
       call FileDefAxis( fid, 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', dtype, UKE-UKS+2 )
       call FileDefAxis( fid, 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', dtype, UKE-UKS+1 )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'CBFZ', 'Boundary factor Center Z', '1', 'CZ', dtype, KA )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, IA )
       call FileDefAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, JA )
    else
       call FileDefAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, IAG )
       call FileDefAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'FBFZ', 'Boundary factor Face Z',   '1', 'CZ', dtype, KA )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'FBFX', 'Boundary factor Face X',   '1', 'CX', dtype, IA )
       call FileDefAxis( fid, 'FBFY', 'Boundary factor Face Y',   '1', 'CY', dtype, JA )
    else
       call FileDefAxis( fid, 'FBFX', 'Boundary factor Face X',   '1', 'CX', dtype, IAG )
       call FileDefAxis( fid, 'FBFY', 'Boundary factor Face Y',   '1', 'CY', dtype, JAG )
    end if

    ! TODO: skip 8 axes below when IO_PNETCDF is true, as all axes are now global
    call FileDefAxis( fid, 'CXG', 'Grid Center Position X (global)', 'm', 'CXG', dtype, IAG )
    call FileDefAxis( fid, 'CYG', 'Grid Center Position Y (global)', 'm', 'CYG', dtype, JAG )
    call FileDefAxis( fid, 'FXG', 'Grid Face Position X (global)',   'm', 'FXG', dtype, IAG+1 )
    call FileDefAxis( fid, 'FYG', 'Grid Face Position Y (global)',   'm', 'FYG', dtype, JAG+1 )

    call FileDefAxis( fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', dtype, IAG )
    call FileDefAxis( fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', dtype, JAG )
    call FileDefAxis( fid, 'FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', dtype, IAG )
    call FileDefAxis( fid, 'FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', dtype, JAG )

    ! associate coordinates
    AXIS_name = (/'x ','y '/)
    call FileDefAssociatedCoordinates( fid, 'lon' , 'longitude',                   &
                                       'degrees_east' , AXIS_name, dtype )
    AXIS_name = (/'xh','y '/)
    call FileDefAssociatedCoordinates( fid, 'lon_uy', 'longitude (half level uy)', &
                                       'degrees_east' , AXIS_name, dtype )
    AXIS_name = (/'x ','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lon_xv', 'longitude (half level xv)', &
                                       'degrees_east' , AXIS_name, dtype )
    AXIS_name = (/'xh','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lon_uv', 'longitude (half level uv)', &
                                       'degrees_east' , AXIS_name, dtype )
    AXIS_name = (/'x ','y '/)
    call FileDefAssociatedCoordinates( fid, 'lat' , 'latitude',                    &
                                       'degrees_north', AXIS_name, dtype )
    AXIS_name = (/'xh','y '/)
    call FileDefAssociatedCoordinates( fid, 'lat_uy', 'latitude (half level uy)',  &
                                       'degrees_north', AXIS_name, dtype )
    AXIS_name = (/'x ','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lat_xv', 'latitude (half level xv)',  &
                                       'degrees_north', AXIS_name, dtype )
    AXIS_name = (/'xh','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lat_uv', 'latitude (half level uv)',  &
                                       'degrees_north', AXIS_name, dtype )

    ! attributes
    if ( .NOT. xy_ ) then
       call FileSetTAttr( fid, 'lz',  'positive', 'down' )
       call FileSetTAttr( fid, 'lzh', 'positive', 'down' )
       call FileSetTAttr( fid, 'uz',  'positive', 'down' )
       call FileSetTAttr( fid, 'uzh', 'positive', 'down' )
       call FileSetTAttr( fid, 'LCZ', 'positive', 'down' )
       call FileSetTAttr( fid, 'LFZ', 'positive', 'down' )
       call FileSetTAttr( fid, 'UCZ', 'positive', 'down' )
       call FileSetTAttr( fid, 'UFZ', 'positive', 'down' )
    end if

    return
  end subroutine FILEIO_def_axes

  !-----------------------------------------------------------------------------
  !> write axis to the file
  subroutine FILEIO_write_axes( &
       fid,   &
       xy     )
    use gtool_file, only: &
       FileWriteAxis,  &
       FileWriteAssociatedCoordinates
    use scale_grid, only: &
       GRID_CZ,    &
       GRID_CX,    &
       GRID_CY,    &
       GRID_FZ,    &
       GRID_FX,    &
       GRID_FY,    &
       GRID_CDZ,   &
       GRID_CDX,   &
       GRID_CDY,   &
       GRID_FDZ,   &
       GRID_FDX,   &
       GRID_FDY,   &
       GRID_CBFZ,  &
       GRID_CBFX,  &
       GRID_CBFY,  &
       GRID_FBFZ,  &
       GRID_FBFX,  &
       GRID_FBFY,  &
       GRID_CXG,   &
       GRID_CYG,   &
       GRID_FXG,   &
       GRID_FYG,   &
       GRID_CBFXG, &
       GRID_CBFYG, &
       GRID_FBFXG, &
       GRID_FBFYG
    use scale_land_grid, only: &
       GRID_LCZ, &
       GRID_LFZ
    use scale_urban_grid, only: &
       GRID_UCZ, &
       GRID_UFZ
    implicit none

    integer, intent(in) :: fid
    logical, intent(in), optional :: xy

    logical :: xy_
    !---------------------------------------------------------------------------

    if ( present(xy) ) then
       xy_ = xy
    else
       xy_ = .false.
    end if

    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'z',   GRID_CZ(KS:KE) )
    end if
    call FileWriteAxis( fid, 'x',   GRID_CX(ISB:IEB) )
    call FileWriteAxis( fid, 'y',   GRID_CY(JSB:JEB) )
    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'zh',  GRID_FZ(KS:KE) )
    end if
    call FileWriteAxis( fid, 'xh',  GRID_FX(ISB:IEB) )
    call FileWriteAxis( fid, 'yh',  GRID_FY(JSB:JEB) )

    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'lz',  GRID_LCZ(LKS:LKE) )
       call FileWriteAxis( fid, 'lzh', GRID_LFZ(LKS:LKE) )
       call FileWriteAxis( fid, 'uz',  GRID_UCZ(UKS:UKE) )
       call FileWriteAxis( fid, 'uzh', GRID_UFZ(UKS:UKE) )
    end if

    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'CZ', GRID_CZ )
    end if
    call FileWriteAxis( fid, 'CX', GRID_CX )
    call FileWriteAxis( fid, 'CY', GRID_CY )
    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'FZ', GRID_FZ )
    end if
    call FileWriteAxis( fid, 'FX', GRID_FX )
    call FileWriteAxis( fid, 'FY', GRID_FY )

    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'CDZ', GRID_CDZ )
    end if
    call FileWriteAxis( fid, 'CDX', GRID_CDX )
    call FileWriteAxis( fid, 'CDY', GRID_CDY )
    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'FDZ', GRID_FDZ )
    end if
    call FileWriteAxis( fid, 'FDX', GRID_FDX )
    call FileWriteAxis( fid, 'FDY', GRID_FDY )

    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'LCZ',  GRID_LCZ )
       call FileWriteAxis( fid, 'LFZ',  GRID_LFZ )
       call FileWriteAxis( fid, 'LCDZ', GRID_LCZ )

       call FileWriteAxis( fid, 'UCZ',  GRID_UCZ )
       call FileWriteAxis( fid, 'UFZ',  GRID_UFZ )
       call FileWriteAxis( fid, 'UCDZ', GRID_UCZ )
    end if

    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'CBFZ', GRID_CBFZ )
    end if
    call FileWriteAxis( fid, 'CBFX', GRID_CBFX )
    call FileWriteAxis( fid, 'CBFY', GRID_CBFY )
    if ( .NOT. xy_ ) then
       call FileWriteAxis( fid, 'FBFZ', GRID_FBFZ )
    end if
    call FileWriteAxis( fid, 'FBFX', GRID_FBFX )
    call FileWriteAxis( fid, 'FBFY', GRID_FBFY )

    call FileWriteAxis( fid, 'CXG', GRID_CXG )
    call FileWriteAxis( fid, 'CYG', GRID_CYG )
    call FileWriteAxis( fid, 'FXG', GRID_FXG )
    call FileWriteAxis( fid, 'FYG', GRID_FYG )

    call FileWriteAxis( fid, 'CBFXG', GRID_CBFXG )
    call FileWriteAxis( fid, 'CBFYG', GRID_CBFYG )
    call FileWriteAxis( fid, 'FBFXG', GRID_FBFXG )
    call FileWriteAxis( fid, 'FBFYG', GRID_FBFYG )


    ! associate coordinates
    call FileWriteAssociatedCoordinates( fid, 'lon' ,   AXIS_LON  (:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lon_uy', AXIS_LONX (:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lon_xv', AXIS_LONY (:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lon_uv', AXIS_LONXY(:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lat' ,   AXIS_LAT  (:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lat_uy', AXIS_LATX (:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lat_xv', AXIS_LATY (:,:) )
    call FileWriteAssociatedCoordinates( fid, 'lat_uv', AXIS_LATXY(:,:) )

    return
  end subroutine FILEIO_write_axes

  !-----------------------------------------------------------------------------
  !> write axes to the file in parallel
  subroutine FILEIO_write_axes_par( &
       fid,   &
       xy     )
    use gtool_file, only: &
       FileWriteAxis,  &
       FileWriteAssociatedCoordinates, &
       FileFlush
    use scale_grid, only: &
       GRID_CZ,    &
       GRID_CX,    &
       GRID_CY,    &
       GRID_FZ,    &
       GRID_FX,    &
       GRID_FY,    &
       GRID_CDZ,   &
       GRID_CDX,   &
       GRID_CDY,   &
       GRID_FDZ,   &
       GRID_FDX,   &
       GRID_FDY,   &
       GRID_CBFZ,  &
       GRID_CBFX,  &
       GRID_CBFY,  &
       GRID_FBFZ,  &
       GRID_FBFX,  &
       GRID_FBFY,  &
       GRID_CXG,   &
       GRID_CYG,   &
       GRID_FXG,   &
       GRID_FYG,   &
       GRID_CDXG,  &
       GRID_CDYG,  &
       GRID_FDXG,  &
       GRID_FDYG,  &
       GRID_CBFXG, &
       GRID_CBFYG, &
       GRID_FBFXG, &
       GRID_FBFYG
    use scale_land_grid, only: &
       GRID_LCZ, &
       GRID_LFZ
    use scale_urban_grid, only: &
       GRID_UCZ, &
       GRID_UFZ
    use scale_process, only: &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    integer, intent(in) :: fid
    logical, intent(in), optional :: xy

    logical :: xy_
    integer :: rankidx(2)
    integer :: start(2)
    !---------------------------------------------------------------------------

    if ( present(xy) ) then
       xy_ = xy
    else
       xy_ = .false.
    end if

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    ! For parallel I/O, not all variables are written by all processes.
    ! 1. Let PRC_myrank 0 writes all z axes
    ! 2. Let processes (rankidx(2) == 0) write x axes  (south-most processes)
    !        rankidx(1) == 0           writes west HALO
    !        rankidx(1) == PRC_NUM_X-1 writes east HALO
    !        others                    writes without HALO
    ! 3. Let processes (rankidx(1) == 0) write y axes  (west-most processes)
    !        rankidx(1) == 0           writes south HALO
    !        rankidx(1) == PRC_NUM_Y-1 writes north HALO
    !        others                    writes without HALO

    if ( .NOT. xy_ .AND. PRC_myrank .EQ. 0 ) then
       start(1) = 1
       call FileWriteAxis( fid, 'z',   GRID_CZ(KS:KE),    start )
       call FileWriteAxis( fid, 'zh',  GRID_FZ(KS:KE),    start )
       call FileWriteAxis( fid, 'lz',  GRID_LCZ(LKS:LKE), start )
       call FileWriteAxis( fid, 'lzh', GRID_LFZ(LKS:LKE), start )
       call FileWriteAxis( fid, 'uz',  GRID_UCZ(UKS:UKE), start )
       call FileWriteAxis( fid, 'uzh', GRID_UFZ(UKS:UKE), start )
    end if

    if ( rankidx(2) .EQ. 0 ) then  ! south most row processes write x/xh
       if ( PRC_PERIODIC_X ) then
          start(1) = ISGB
       else
          start(1) = ISGA
       end if
       call FileWriteAxis( fid, 'x',   GRID_CX(ISB:IEB), start )
       call FileWriteAxis( fid, 'xh',  GRID_FX(ISB:IEB), start )
    end if

    if ( rankidx(1) .EQ. 0 ) then  ! west most column processes write y/yh
       if ( PRC_PERIODIC_Y ) then
          start(1) = JSGB
       else
          start(1) = JSGA
       end if
       call FileWriteAxis( fid, 'y',   GRID_CY(JSB:JEB), start )
       call FileWriteAxis( fid, 'yh',  GRID_FY(JSB:JEB), start )
    end if

    if ( .NOT. xy_ .AND. PRC_myrank .EQ. 0 ) then
       start(1) = 1
       call FileWriteAxis( fid, 'CZ',   GRID_CZ,   start )
       call FileWriteAxis( fid, 'CDZ',  GRID_CDZ,  start )
       call FileWriteAxis( fid, 'LCZ',  GRID_LCZ,  start )
       call FileWriteAxis( fid, 'LFZ',  GRID_LFZ,  start )
       call FileWriteAxis( fid, 'LCDZ', GRID_LCZ,  start )
       call FileWriteAxis( fid, 'UCZ',  GRID_UCZ,  start )
       call FileWriteAxis( fid, 'UFZ',  GRID_UFZ,  start )
       call FileWriteAxis( fid, 'UCDZ', GRID_UCZ,  start )
       call FileWriteAxis( fid, 'CBFZ', GRID_CBFZ, start )
       call FileWriteAxis( fid, 'FBFZ', GRID_FBFZ, start )
       call FileWriteAxis( fid, 'FZ',   GRID_FZ,   start )
       call FileWriteAxis( fid, 'FDZ',  GRID_FDZ,  start )
    end if

    if ( PRC_myrank .EQ. 0 ) then ! rank 0 writes entire global axes
      ! these axes always include halos when written to file regardless of PRC_PERIODIC
       start(1) = 1
       call FileWriteAxis( fid, 'CX',   GRID_CXG,   start )
       call FileWriteAxis( fid, 'CDX',  GRID_CDXG,  start )
       call FileWriteAxis( fid, 'CBFX', GRID_CBFXG, start )
       call FileWriteAxis( fid, 'FBFX', GRID_FBFX,  start )
       call FileWriteAxis( fid, 'FDX',  GRID_FDX,   start )
       call FileWriteAxis( fid, 'FX',   GRID_FXG,   start )
       call FileWriteAxis( fid, 'CY',   GRID_CYG,   start )
       call FileWriteAxis( fid, 'CDY',  GRID_CDYG,  start )
       call FileWriteAxis( fid, 'CBFY', GRID_CBFYG, start )
       call FileWriteAxis( fid, 'FBFY', GRID_FBFYG, start )
       call FileWriteAxis( fid, 'FDY',  GRID_FDYG,  start )
       call FileWriteAxis( fid, 'FY',   GRID_FYG,   start )
    end if

    ! global axes: skip 8 axes below when IO_PNETCDF is true, as all axes are now global
    ! call FileWriteAxis( fid, 'CXG', GRID_CXG )
    ! call FileWriteAxis( fid, 'CYG', GRID_CYG )
    ! call FileWriteAxis( fid, 'FXG', GRID_FXG )
    ! call FileWriteAxis( fid, 'FYG', GRID_FYG )
    !
    ! call FileWriteAxis( fid, 'CBFXG', GRID_CBFXG )
    ! call FileWriteAxis( fid, 'CBFYG', GRID_CBFYG )
    ! call FileWriteAxis( fid, 'FBFXG', GRID_FBFXG )
    ! call FileWriteAxis( fid, 'FBFYG', GRID_FBFYG )

    ! associate coordinates (dimensions: x/xh, y/yh)
    ! They are allocated with sizes (IMAXB,JMAXB)
    if ( PRC_PERIODIC_X ) then
       start(1) = ISGB
    else
       start(1) = ISGA
    end if
    if ( PRC_PERIODIC_Y ) then
       start(2) = JSGB
    else
       start(2) = JSGA
    end if
    call FileWriteAssociatedCoordinates( fid, 'lon' ,   AXIS_LON  (:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lon_uy', AXIS_LONX (:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lon_xv', AXIS_LONY (:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lon_uv', AXIS_LONXY(:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lat' ,   AXIS_LAT  (:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lat_uy', AXIS_LATX (:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lat_xv', AXIS_LATY (:,:), start )
    call FileWriteAssociatedCoordinates( fid, 'lat_uv', AXIS_LATXY(:,:), start )

    call FileFlush( fid )  ! PnetCDF only

    return
  end subroutine FILEIO_write_axes_par

  !-----------------------------------------------------------------------------
  !> Define a variable to file
  subroutine FILEIO_def_var( &
       fid,      &
       vid,      &
       varname,  &
       desc,     &
       unit,     &
       axistype, &
       datatype, &
       timeintv, &
       nsteps    )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileDefineVariable
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    integer,            intent(in)  :: fid      !< file ID
    integer,            intent(out) :: vid      !< variable ID
    character(len=*),   intent(in)  :: varname  !< name        of the variable
    character(len=*),   intent(in)  :: desc     !< description of the variable
    character(len=*),   intent(in)  :: unit     !< unit        of the variable
    character(len=*),   intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*),   intent(in)  :: datatype !< data type (REAL8/REAL4/default)
    real(RP), optional, intent(in)  :: timeintv !< time interval [sec]
    integer,  optional, intent(in)  :: nsteps   !< number of time steps

    integer           :: dtype, ndims, elm_size
    character(len=2)  :: dims(3)
    real(DP)          :: time_interval
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( datatype .EQ. 'REAL8' ) then
       dtype = File_REAL8
       elm_size = 4
    elseif( datatype .EQ. 'REAL4' ) then
       dtype = File_REAL4
       elm_size = 8
    else
       if ( RP .EQ. 8 ) then
          dtype = File_REAL8
          elm_size = 8
       elseif( RP .EQ. 4 ) then
          dtype = File_REAL4
          elm_size = 4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    if ( axistype .EQ. 'Z' ) then        ! 1D variable
       ndims = 1
       dims(1) = 'z'
       write_buf_amount = write_buf_amount + KA * elm_size
    elseif( axistype .EQ. 'X' ) then
       ndims = 1
       dims(1) = 'x'
       write_buf_amount = write_buf_amount + IA * elm_size
    elseif( axistype .EQ. 'Y' ) then
       ndims = 1
       dims(1) = 'y'
       write_buf_amount = write_buf_amount + JA * elm_size
    elseif ( axistype .EQ. 'XY' ) then   ! 2D variable
       ndims = 2
       dims(1) = 'x'
       dims(2) = 'y'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif ( axistype .EQ. 'UY' ) then
       ndims = 2
       dims(1) = 'xh'
       dims(2) = 'y'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif ( axistype .EQ. 'XV' ) then
       ndims = 2
       dims(1) = 'x'
       dims(2) = 'yh'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif ( axistype .EQ. 'UV' ) then
       ndims = 2
       dims(1) = 'xh'
       dims(2) = 'yh'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif( axistype .EQ. 'ZX' ) then
       ndims = 2
       dims(1) = 'z'
       dims(2) = 'x'
       write_buf_amount = write_buf_amount + KA * IA * elm_size
    elseif ( axistype .EQ. 'ZXY' ) then  ! 3D variable
       ndims = 3
       dims = (/'z','x','y'/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype .EQ. 'ZHXY' ) then
       ndims = 3
       dims = (/'zh','x ','y '/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype .EQ. 'ZXHY' ) then
       ndims = 3
       dims = (/'z ','xh','y '/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype .EQ. 'ZXYH' ) then
       ndims = 3
       dims = (/'z ','x ','yh'/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype .EQ. 'Land' ) then
       ndims = 3
       dims = (/'lz','x ','y '/)
       write_buf_amount = write_buf_amount + LKMAX * IA * JA * elm_size
    elseif( axistype .EQ. 'Urban' ) then
       ndims = 3
       dims = (/'uz','x ','y '/)
       write_buf_amount = write_buf_amount + UKMAX * IA * JA * elm_size
    elseif ( axistype .EQ. 'XYT' ) then  ! 3D variable with time dimension
       ndims = 2
       dims(1) = 'x'
       dims(2) = 'y'
       if ( present(nsteps) ) then
          write_buf_amount = write_buf_amount + IA * JA * elm_size * nsteps
       else
          write_buf_amount = write_buf_amount + IA * JA * elm_size
       end if
    elseif ( axistype .EQ. 'ZXYT' ) then ! 4D variable
       ndims = 3
       dims = (/'z','x','y'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + KA * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( present(timeintv) ) then  ! 3D/4D variable with time dimension
      time_interval = timeintv
      call FileDefineVariable( fid, vid, varname, desc, unit, ndims, dims, dtype, & ! [IN]
           tint=time_interval ) ! [IN]
    else
      call FileDefineVariable( fid, vid, varname, desc, unit, ndims, dims, dtype ) ! [IN]
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_def_var

  !-----------------------------------------------------------------------------
  !> Write 1D data to file
  subroutine FILEIO_write_var_1D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       axistype  )
    use gtool_file, only: &
       FileWriteVar
    use scale_process, only: &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use scale_time, only: &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    integer,          intent(in)  :: fid      !< file ID
    integer,          intent(in)  :: vid      !< netCDF variable ID
    real(RP),         intent(in)  :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)

    integer :: dim1_max, dim1_S, dim1_E
    integer :: rankidx(2)
    integer :: start(1)
    logical :: exec = .TRUE.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( axistype .EQ. 'Z' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
       start(1) = 1
       if ( IO_PNETCDF .AND. PRC_myrank .GT. 0 ) &
          exec = .FALSE.  ! only rank 0 writes
    elseif( axistype .EQ. 'X' ) then
       dim1_max = IMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       if ( PRC_PERIODIC_X ) then
          start(1) = ISGB
       else
          start(1) = ISGA
       endif
       if ( IO_PNETCDF .AND. rankidx(2) .GT. 0 ) &
          exec = .FALSE.  ! only south most row processes write
    elseif( axistype .EQ. 'Y' ) then
       dim1_max = JMAXB
       dim1_S   = JSB
       dim1_E   = JEB
       if ( PRC_PERIODIC_Y ) then
          start(1) = JSGB
       else
          start(1) = JSGA
       endif
       if ( IO_PNETCDF .AND. rankidx(1) .GT. 0 ) &
          exec = .FALSE.  ! only west most column processes write
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( exec ) call FileWriteVar( vid, var(dim1_S:dim1_E), NOWSEC, NOWSEC, start ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_1D

  !-----------------------------------------------------------------------------
  !> Write 2D data to file
  subroutine FILEIO_write_var_2D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       axistype, &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file, only: &
       FileWriteVar
    use scale_process, only: &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use scale_time, only: &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    integer,          intent(in)  :: fid      !< file ID
    integer,          intent(in)  :: vid      !< netCDF variable ID
    real(RP),         intent(in)  :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    logical, optional, intent(in) :: nohalo   !< switch whether include halo data or not (default=false)

    real(RP)              :: varhalo( size(var(:,1)), size(var(1,:)) )

    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E

    integer :: i, j
    logical :: nohalo_
    integer :: rankidx(2)
    integer :: start(2)
    logical :: exec = .TRUE.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( IO_PNETCDF ) then
       if ( PRC_PERIODIC_X ) then
          start(1) = ISGB
       else
          start(1) = ISGA
       endif
       if ( PRC_PERIODIC_Y ) then
          start(2) = JSGB
       else
          start(2) = JSGA
       endif
    end if

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    if ( axistype .EQ. 'XY' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif ( axistype .EQ. 'UY' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif ( axistype .EQ. 'XV' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif ( axistype .EQ. 'UV' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    elseif( axistype .EQ. 'ZX' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       start(2) = start(1)
       start(1) = 1
       if ( IO_PNETCDF .AND. rankidx(2) .GT. 0 ) &
          exec = .FALSE.  ! only south most row processes write
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( exec ) then
       varhalo(:,:) = var(:,:)

       if ( nohalo_ ) then
          ! W halo
          do j = 1, JA
          do i = 1, IS-1
             varhalo(i,j) = RMISS
          end do
          end do
          ! E halo
          do j = 1, JA
          do i = IE+1, IA
             varhalo(i,j) = RMISS
          end do
          end do
          ! S halo
          do j = 1, JS-1
          do i = 1, IA
             varhalo(i,j) = RMISS
          end do
          end do
          ! N halo
          do j = JE+1, JA
          do i = 1, IA
             varhalo(i,j) = RMISS
          end do
          end do
       end if

       call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), NOWSEC, NOWSEC, start ) ! [IN]
    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_2D

  !-----------------------------------------------------------------------------
  !> Write 3D data to file
  subroutine FILEIO_write_var_3D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       axistype, &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file, only: &
       FileWriteVar
    use scale_process, only: &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use scale_time, only: &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    integer,           intent(in)  :: fid        !< file ID
    integer,           intent(in)  :: vid        !< netCDF variable ID
    real(RP),          intent(in)  :: var(:,:,:) !< value of the variable
    character(len=*),  intent(in)  :: varname    !< name        of the variable
    character(len=*),  intent(in)  :: axistype   !< axis type (Z/X/Y)
    logical, optional, intent(in)  :: nohalo     !< include halo data?

    real(RP)         :: varhalo( size(var(:,1,1)), size(var(1,:,1)), size(var(1,1,:)) )

    integer          :: dim1_max, dim1_S, dim1_E
    integer          :: dim2_max, dim2_S, dim2_E
    integer          :: dim3_max, dim3_S, dim3_E

    integer :: i, j, k
    logical :: nohalo_
    integer :: start(3)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    if ( IO_PNETCDF ) then
       start(1) = 1
       if ( PRC_PERIODIC_X ) then
          start(2) = ISGB
       else
          start(2) = ISGA
       endif
       if ( PRC_PERIODIC_Y ) then
          start(3) = JSGB
       else
          start(3) = JSGA
       endif
    end if

    if ( axistype .EQ. 'ZXY' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype .EQ. 'ZHXY' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype .EQ. 'ZXHY' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype .EQ. 'ZXYH' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype .EQ. 'Land' ) then
       dim1_max = LKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    elseif( axistype .EQ. 'Urban' ) then
       dim1_max = UKMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = UKS
       dim1_E   = UKE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    varhalo(:,:,:) = var(:,:,:)

    if ( nohalo_ ) then
       ! W halo
       do k = 1, dim1_max
       do j = 1, JA
       do i = 1, IS-1
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
       ! E halo
       do k = 1, dim1_max
       do j = 1, JA
       do i = IE+1, IA
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
       ! S halo
       do k = 1, dim1_max
       do j = 1, JS-1
       do i = 1, IA
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
       ! N halo
       do k = 1, dim1_max
       do j = JE+1, JA
       do i = 1, IA
          varhalo(k,i,j) = RMISS
       end do
       end do
       end do
    end if

    call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                       NOWSEC, NOWSEC, start ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_3D

  !-----------------------------------------------------------------------------
  !> Write 3D data with time dimension to file
  subroutine FILEIO_write_var_3D_t( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       axistype, &
       timeintv, &
       timetarg, &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file, only: &
       FileWriteVar
    use scale_process, only: &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    integer,           intent(in)  :: fid          !< file ID
    integer,           intent(in)  :: vid          !< netCDF variable ID
    real(RP),          intent(in)  :: var(:,:,:)   !< value of the variable
    character(len=*),  intent(in)  :: varname      !< name of the variable
    character(len=*),  intent(in)  :: axistype     !< axis type (X/Y/Time)
    real(RP),          intent(in)  :: timeintv     !< time interval [sec]
    integer, optional, intent(in)  :: timetarg     !< target timestep (optional)
    logical, optional, intent(in)  :: nohalo       !< include halo data?

    real(RP)         :: varhalo( size(var(:,1,1)), size(var(1,:,1)) )

    integer          :: dim1_max, dim1_S, dim1_E
    integer          :: dim2_max, dim2_S, dim2_E

    real(DP) :: time_interval, nowtime

    integer :: step
    integer :: i, j, n
    logical :: nohalo_
    integer :: start(3)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    time_interval = timeintv
    step = size(var(ISB,JSB,:))

    if ( axistype .EQ. 'XYT' ) then
       dim1_max = IMAXB
       dim2_max = JMAXB
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( IO_PNETCDF ) then
       if ( PRC_PERIODIC_X ) then
          start(1) = ISGB
       else
          start(1) = ISGA
       endif
       if ( PRC_PERIODIC_Y ) then
          start(2) = JSGB
       else
          start(2) = JSGA
       endif
       ! start(3) will be calculated in file_write_var_par()
    end if

    if ( present(timetarg) ) then
       varhalo(:,:) = var(:,:,timetarg)

       if ( nohalo_ ) then
          ! W halo
          do j = 1, JA
          do i = 1, IS-1
             varhalo(i,j) = RMISS
          end do
          end do
          ! E halo
          do j = 1, JA
          do i = IE+1, IA
             varhalo(i,j) = RMISS
          end do
          end do
          ! S halo
          do j = 1, JS-1
          do i = 1, IA
             varhalo(i,j) = RMISS
          end do
          end do
          ! N halo
          do j = JE+1, JA
          do i = 1, IA
             varhalo(i,j) = RMISS
          end do
          end do
       end if

       nowtime = (timetarg-1) * time_interval
       call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), nowtime, nowtime, start ) ! [IN]
    else
       nowtime = 0.0_DP
       do n = 1, step
          varhalo(:,:) = var(:,:,n)

          if ( nohalo_ ) then
             ! W halo
             do j = 1, JA
             do i = 1, IS-1
                varhalo(i,j) = RMISS
             end do
             end do
             ! E halo
             do j = 1, JA
             do i = IE+1, IA
                varhalo(i,j) = RMISS
             end do
             end do
             ! S halo
             do j = 1, JS-1
             do i = 1, IA
                varhalo(i,j) = RMISS
             end do
             end do
             ! N halo
             do j = JE+1, JA
             do i = 1, IA
                varhalo(i,j) = RMISS
             end do
             end do
          end if

          call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), nowtime, nowtime, start ) ! [IN]
          nowtime = nowtime + time_interval
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_3D_t

  !-----------------------------------------------------------------------------
  !> Write 4D data to file
  subroutine FILEIO_write_var_4D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       axistype, &
       timeintv, &
       timetarg, &
       nohalo    )
    use gtool_file, only: &
       RMISS
    use gtool_file, only: &
       FileWriteVar
    use scale_process, only: &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    integer,           intent(in)  :: fid          !< file ID
    integer,           intent(in)  :: vid          !< netCDF variable ID
    real(RP),          intent(in)  :: var(:,:,:,:) !< value of the variable
    character(len=*),  intent(in)  :: varname      !< name        of the variable
    character(len=*),  intent(in)  :: axistype     !< axis type (Z/X/Y/Time)
    real(RP),          intent(in)  :: timeintv     !< time interval [sec]
    integer, optional, intent(in)  :: timetarg     !< target timestep (optional)
    logical, optional, intent(in)  :: nohalo       !< include halo data?

    real(RP)         :: varhalo( size(var(:,1,1,1)), size(var(1,:,1,1)), size(var(1,1,:,1)) )

    integer          :: dim1_max, dim1_S, dim1_E
    integer          :: dim2_max, dim2_S, dim2_E
    integer          :: dim3_max, dim3_S, dim3_E

    real(DP) :: time_interval, nowtime

    integer :: step
    integer :: i, j, k, n
    logical :: nohalo_
    integer :: start(4), count(4)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    if ( IO_PNETCDF ) then
       start(1) = 1
       if ( PRC_PERIODIC_X ) then
          start(2) = ISGB
       else
          start(2) = ISGA
       endif
       if ( PRC_PERIODIC_Y ) then
          start(3) = JSGB
       else
          start(3) = JSGA
       endif
       ! start(4) will be calculated in file_write_var_par()
    end if

    time_interval = timeintv
    step = size(var(KS,ISB,JSB,:))

    if ( axistype .EQ. 'ZXYT' ) then
       dim1_max = KMAX
       dim2_max = IMAXB
       dim3_max = JMAXB
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( present(timetarg) ) then
       varhalo(:,:,:) = var(:,:,:,timetarg)

       if ( nohalo_ ) then
          ! W halo
          do k = 1, dim1_max
          do j = 1, JA
          do i = 1, IS-1
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
          ! E halo
          do k = 1, dim1_max
          do j = 1, JA
          do i = IE+1, IA
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
          ! S halo
          do k = 1, dim1_max
          do j = 1, JS-1
          do i = 1, IA
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
          ! N halo
          do k = 1, dim1_max
          do j = JE+1, JA
          do i = 1, IA
             varhalo(k,i,j) = RMISS
          end do
          end do
          end do
       end if

       nowtime = (timetarg-1) * time_interval
       call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                          nowtime, nowtime, start ) ! [IN]
    else
       nowtime = 0.0_DP
       do n = 1, step
          varhalo(:,:,:) = var(:,:,:,n)

          if ( nohalo_ ) then
             ! W halo
             do k = 1, dim1_max
             do j = 1, JA
             do i = 1, IS-1
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
             ! E halo
             do k = 1, dim1_max
             do j = 1, JA
             do i = IE+1, IA
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
             ! S halo
             do k = 1, dim1_max
             do j = 1, JS-1
             do i = 1, IA
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
             ! N halo
             do k = 1, dim1_max
             do j = JE+1, JA
             do i = 1, IA
                varhalo(k,i,j) = RMISS
             end do
             end do
             end do
          end if

          call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                             nowtime, nowtime, start ) ! [IN]
          nowtime = nowtime + time_interval
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_4D

end module scale_fileio
