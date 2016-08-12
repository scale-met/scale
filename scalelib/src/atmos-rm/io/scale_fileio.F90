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

  ! global star and count
  integer,  private,      save :: startXY(3),    countXY(3)
  integer,  private,      save :: startZX(2),    countZX(2)
  integer,  private,      save :: startZXY(4),   countZXY(4)
  integer,  private,      save :: startLAND(3),  countLAND(3)
  integer,  private,      save :: startURBAN(3), countURBAN(3)

  ! MPI element datatype for restart variables
  integer,  private,      save :: etype

  ! MPI derived datatypes
  integer,  private,      save :: centerTypeXY
  integer,  private,      save :: centerTypeZX
  integer,  private,      save :: centerTypeZXY
  integer,  private,      save :: centerTypeLAND
  integer,  private,      save :: centerTypeURBAN

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

  !-----------------------------------------------------------------------------
  !> construct MPI derived datatypes for read buffers
  subroutine Construct_Derived_Datatype
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    use MPI
    implicit none

    integer err, order
    integer sizes(3), subsizes(3), sub_off(3)
    integer XSG, YSG
    integer XS, XE, YS, YE
    !---------------------------------------------------------------------------
    order = MPI_ORDER_FORTRAN

    centerTypeXY    = MPI_DATATYPE_NULL
    centerTypeZX    = MPI_DATATYPE_NULL
    centerTypeZXY   = MPI_DATATYPE_NULL
    centerTypeLAND  = MPI_DATATYPE_NULL
    centerTypeURBAN = MPI_DATATYPE_NULL

    etype = MPI_FLOAT
    if ( RP .EQ. 8 ) etype = MPI_DOUBLE_PRECISION

    ! for axistype == 'XY'
    startXY(1) = IS_inG - IHALO
    startXY(2) = JS_inG - JHALO
    countXY(1) = IA
    countXY(2) = JA

    ! for axistype == 'ZXY'
    startZXY(1)   = 1
    startZXY(2:3) = startXY(1:2)
    countZXY(1)   = KMAX
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
    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeZXY, err)
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
    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeLAND, err)
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
    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeURBAN, err)
    call MPI_Type_commit(centerTypeURBAN, err)

    ! for axistype == 'ZX'
    sizes(1)    = KA
    subsizes(1) = KMAX
    sizes(1)    = KA
    subsizes(1) = KMAX
    startZX(1)  = KHALO+1
    countZX(1)  = KHALO
    sub_off(1)  = KHALO   ! MPI start index starts with 0

    startZX(2)  = IS_inG - IHALO
    countZX(2)  = IA
    sizes(2)    = IA
    subsizes(2) = IMAXB
    sub_off(2)  = ISB-1   ! MPI start index starts with 0
    call MPI_Type_create_subarray(2, sizes, subsizes, sub_off, order, etype, centerTypeZX, err)
    call MPI_Type_commit(centerTypeZX, err)

  end subroutine Construct_Derived_Datatype

  !-----------------------------------------------------------------------------
  !> free MPI derived datatypes
  subroutine Free_Derived_Datatype
    use MPI
    implicit none
    integer err

    if ( .NOT. IO_PNETCDF ) return

    if ( centerTypeXY    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeXY,    err)
    if ( centerTypeZX    .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZX,    err)
    if ( centerTypeZXY   .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZXY,   err)
    if ( centerTypeLAND  .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeLAND,  err)
    if ( centerTypeURBAN .NE. MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeURBAN, err)

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
       PRC_NUM_Y
    use MPI
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer               :: XS, XE, YS, YE  ! local buffer indices
    integer               :: XSG, YSG        ! global array indices
    integer               :: start(1)   ! start offset of globale variable
    integer               :: count(1)   ! request length to the globale variable
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var: ', trim(varname)

    if ( axistype .EQ. 'Z' ) then
       start(1) = 1
       count(1) = KMAX
       call FileRead( var, fid, varname, step, KMAX, etype, start, count )
    elseif( axistype .EQ. 'X' ) then
       ! read data and halos in one call
       start(1) = IS_inG - IHALO
       count(1) = IA
       call FileRead( var, fid, varname, step, IA, etype, start, count )
    elseif( axistype .EQ. 'Y' ) then
       ! read data and halos in one call
       start(1) = JS_inG - JHALO
       count(1) = JA
       call FileRead( var, fid, varname, step, JA, etype, start, count )
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
       PRC_NUM_Y
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

    ! read data and halos in one call
    if ( axistype .EQ. 'XY' ) then
       call FileRead( var, fid, varname, step, IA*JA, etype, startXY, countXY )
    elseif( axistype .EQ. 'ZX' ) then
       call FileRead( var, fid, varname, step, 1, centerTypeZX, startZX, countZX )
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
       PRC_NUM_Y
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
    else if ( axistype .EQ. 'XYT' ) then
       startXY(3) = 1
       countXY(3) = step
       call FileRead( var, fid, varname, step, step*IA*JA, etype, startXY, countXY )
    else if ( axistype .EQ. 'Land' ) then
       call FileRead( var, fid, varname, step, 1, centerTypeLAND, startLAND, countLAND )
    else if ( axistype .EQ. 'Urban' ) then
       call FileRead( var, fid, varname, step, 1, centerTypeURBAN, startURBAN, countURBAN )
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
       PRC_NUM_Y
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
       startZXY(4) = 1
       countZXY(4) = step
       call FileRead( var, fid, varname, step, step, centerTypeZXY, startZXY, countZXY )
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
       PRC_2Drank, &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
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
    character(len=8)  :: logical_str
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
       call FileSetGlobalAttribute( fid, "IHALO",  (/IHALO/) )
       call FileSetGlobalAttribute( fid, "JHALO",  (/JHALO/) )
       logical_str = "false"
       if (PRC_PERIODIC_X .AND. .NOT. IO_PNETCDF) logical_str = "true"
       call FileSetGlobalAttribute( fid, "PRC_PERIODIC_X",  trim(logical_str) )
       logical_str = "false"
       if (PRC_PERIODIC_Y .AND. .NOT. IO_PNETCDF) logical_str = "true"
       call FileSetGlobalAttribute( fid, "PRC_PERIODIC_Y",  trim(logical_str) )
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
       PRC_NUM_Y
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
       call FileDefAxis( fid, 'z',   'Z',               'm', 'z',   dtype, KMAX )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'x',   'X',               'm', 'x',   dtype, IMAXB )
       call FileDefAxis( fid, 'y',   'Y',               'm', 'y',   dtype, JMAXB )
    else
       call FileDefAxis( fid, 'x',   'X',               'm', 'x',   dtype, IAG )
       call FileDefAxis( fid, 'y',   'Y',               'm', 'y',   dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'zh',  'Z (half level)',  'm', 'zh',  dtype, KMAX )
    end if
    if ( .NOT. IO_PNETCDF ) then
       call FileDefAxis( fid, 'xh',  'X (half level)',  'm', 'xh',  dtype, IMAXB )
       call FileDefAxis( fid, 'yh',  'Y (half level)',  'm', 'yh',  dtype, JMAXB )
    else
       call FileDefAxis( fid, 'xh',  'X (half level)',  'm', 'xh',  dtype, IAG )
       call FileDefAxis( fid, 'yh',  'Y (half level)',  'm', 'yh',  dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'lz',  'LZ',              'm', 'lz',  dtype, LKMAX )
       call FileDefAxis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, LKMAX )
       call FileDefAxis( fid, 'uz',  'UZ',              'm', 'uz',  dtype, UKMAX )
       call FileDefAxis( fid, 'uzh', 'UZ (half level)', 'm', 'uzh', dtype, UKMAX )
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
       call FileDefAxis( fid, 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', dtype, LKMAX )
       call FileDefAxis( fid, 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', dtype, LKMAX+1 )
       call FileDefAxis( fid, 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', dtype, LKMAX )

       call FileDefAxis( fid, 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', dtype, UKMAX )
       call FileDefAxis( fid, 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', dtype, UKMAX+1 )
       call FileDefAxis( fid, 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', dtype, UKMAX )
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
  !> write axes to a single shared file in parallel
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
       PRC_2Drank
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
       start(1) = ISGA
       call FileWriteAxis( fid, 'x',   GRID_CX(ISB:IEB), start )
       call FileWriteAxis( fid, 'xh',  GRID_FX(ISB:IEB), start )
    end if

    if ( rankidx(1) .EQ. 0 ) then  ! west most column processes write y/yh
       start(1) = JSGA
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
    start(1) = ISGA
    start(2) = JSGA
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
       PRC_2Drank
    use scale_time, only: &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    integer,          intent(in)  :: fid      !< file ID
    integer,          intent(in)  :: vid      !< netCDF variable ID
    real(RP),         intent(in)  :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)

    integer :: dim1_S, dim1_E
    integer :: rankidx(2)
    integer :: start(1)         ! used only when IO_PNETCDF is .true.
    logical :: exec = .TRUE.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( axistype .EQ. 'Z' ) then
       dim1_S   = KS
       dim1_E   = KE
       start(1) = 1
       if ( IO_PNETCDF .AND. PRC_myrank .GT. 0 ) &
          exec = .FALSE.  ! only rank 0 writes
    elseif( axistype .EQ. 'X' ) then
       dim1_S   = ISB
       dim1_E   = IEB
       start(1) = ISGA
       if ( IO_PNETCDF .AND. rankidx(2) .GT. 0 ) &
          exec = .FALSE.  ! only south most row processes write
    elseif( axistype .EQ. 'Y' ) then
       dim1_S   = JSB
       dim1_E   = JEB
       start(1) = JSGA
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
       PRC_NUM_X, &
       PRC_NUM_Y
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

    integer               :: dim1_S, dim1_E
    integer               :: dim2_S, dim2_E

    integer :: i, j
    logical :: nohalo_
    integer :: rankidx(2)
    integer :: start(2)         ! used only when IO_PNETCDF is .true.
    logical :: exec = .TRUE.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = ISGA
    start(2) = JSGA

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    if ( axistype .EQ. 'XY' .OR. &
         axistype .EQ. 'UY' .OR. &
         axistype .EQ. 'XV' .OR. &
         axistype .EQ. 'UV' ) then
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       if ( IO_PNETCDF ) then
          if ( rankidx(1) .EQ. 0             ) dim1_S = 1
          if ( rankidx(1) .EQ. PRC_NUM_X - 1 ) dim1_E = IA
          if ( rankidx(2) .EQ. 0             ) dim2_S = 1
          if ( rankidx(2) .EQ. PRC_NUM_Y - 1 ) dim2_E = JA
       end if
    elseif( axistype .EQ. 'ZX' ) then
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       start(2) = start(1)
       start(1) = 1
       if ( IO_PNETCDF .AND. rankidx(2) .GT. 0 ) then
          exec = .FALSE.  ! only south most row processes write
          if ( rankidx(1) .EQ. 0             ) dim2_S = 1
          if ( rankidx(1) .EQ. PRC_NUM_X - 1 ) dim2_E = IA
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( exec ) then
       if ( nohalo_ ) then ! fill halo cells with RMISS
          varhalo(:,:) = var(:,:)

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

          call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), NOWSEC, NOWSEC, start ) ! [IN]
       else
          call FileWriteVar( vid, var(dim1_S:dim1_E,dim2_S:dim2_E), NOWSEC, NOWSEC, start ) ! [IN]
       end if

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
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
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

    integer          :: dim1_S, dim1_E, dim1_max
    integer          :: dim2_S, dim2_E
    integer          :: dim3_S, dim3_E

    integer :: i, j, k
    logical :: nohalo_
    integer :: rankidx(2)
    integer :: start(3)         ! used only when IO_PNETCDF is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA

    dim2_S   = ISB
    dim2_E   = IEB
    dim3_S   = JSB
    dim3_E   = JEB
    if ( IO_PNETCDF ) then
       if ( rankidx(1) .EQ. 0             ) dim2_S = 1
       if ( rankidx(1) .EQ. PRC_NUM_X - 1 ) dim2_E = IA
       if ( rankidx(2) .EQ. 0             ) dim3_S = 1
       if ( rankidx(2) .EQ. PRC_NUM_Y - 1 ) dim3_E = JA
    end if

    if ( axistype .EQ. 'ZXY'  .OR. &
         axistype .EQ. 'ZHXY' .OR. &
         axistype .EQ. 'ZXHY' .OR. &
         axistype .EQ. 'ZXYH' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype .EQ. 'Land' ) then
       dim1_max = LKMAX
       dim1_S   = LKS
       dim1_E   = LKE
    elseif( axistype .EQ. 'Urban' ) then
       dim1_max = UKMAX
       dim1_S   = UKS
       dim1_E   = UKE
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( nohalo_ ) then
       varhalo(:,:,:) = var(:,:,:)

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

       call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                          NOWSEC, NOWSEC, start ) ! [IN]
    else
       call FileWriteVar( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                          NOWSEC, NOWSEC, start ) ! [IN]
    end if

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
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
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

    integer          :: dim1_S, dim1_E
    integer          :: dim2_S, dim2_E

    real(DP) :: time_interval, nowtime

    integer :: step
    integer :: i, j, n
    logical :: nohalo_
    integer :: rankidx(2)
    integer :: start(3)         ! used only when IO_PNETCDF is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    time_interval = timeintv
    step = size(var(ISB,JSB,:))

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( axistype .EQ. 'XYT' ) then
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       if ( IO_PNETCDF ) then
          if ( rankidx(1) .EQ. 0             ) dim1_S = 1
          if ( rankidx(1) .EQ. PRC_NUM_X - 1 ) dim1_E = IA
          if ( rankidx(2) .EQ. 0             ) dim2_S = 1
          if ( rankidx(2) .EQ. PRC_NUM_Y - 1 ) dim2_E = JA
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    start(1) = ISGA
    start(2) = JSGA
    ! start(3) will be calculated in file_write_var_par()

    if ( present(timetarg) ) then
       nowtime = (timetarg-1) * time_interval

       if ( nohalo_ ) then
          varhalo(:,:) = var(:,:,timetarg)

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

          call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), nowtime, nowtime, start ) ! [IN]
       else
          call FileWriteVar( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,timetarg), nowtime, nowtime, start ) ! [IN]
       end if
    else
       nowtime = 0.0_DP
       do n = 1, step
          if ( nohalo_ ) then
             varhalo(:,:) = var(:,:,n)

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

             call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), nowtime, nowtime, start ) ! [IN]
          else
             call FileWriteVar( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,n), nowtime, nowtime, start ) ! [IN]
          end if
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
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
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

    integer          :: dim1_S, dim1_E, dim1_max
    integer          :: dim2_S, dim2_E
    integer          :: dim3_S, dim3_E

    real(DP) :: time_interval, nowtime

    integer :: step
    integer :: i, j, k, n
    logical :: nohalo_
    integer :: rankidx(2)
    integer :: start(4)         ! used only when IO_PNETCDF is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if ( present(nohalo) ) nohalo_ = nohalo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA
    ! start(4) will be calculated in file_write_var_par()

    time_interval = timeintv
    step = size(var(KS,ISB,JSB,:))

    if ( axistype .EQ. 'ZXYT' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       if ( IO_PNETCDF ) then
          if ( rankidx(1) .EQ. 0             ) dim2_S = 1
          if ( rankidx(1) .EQ. PRC_NUM_X - 1 ) dim2_E = IA
          if ( rankidx(2) .EQ. 0             ) dim3_S = 1
          if ( rankidx(2) .EQ. PRC_NUM_Y - 1 ) dim3_E = JA
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if ( present(timetarg) ) then
       nowtime = (timetarg-1) * time_interval

       if ( nohalo_ ) then
          varhalo(:,:,:) = var(:,:,:,timetarg)

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

          call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                             nowtime, nowtime, start ) ! [IN]
       else
          call FileWriteVar( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,timetarg), &
                             nowtime, nowtime, start ) ! [IN]
       end if
    else
       nowtime = 0.0_DP
       do n = 1, step
          if ( nohalo_ ) then
             varhalo(:,:,:) = var(:,:,:,n)

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

             call FileWriteVar( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                                nowtime, nowtime, start ) ! [IN]
          else
             call FileWriteVar( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,n), &
                                nowtime, nowtime, start ) ! [IN]
          end if
          nowtime = nowtime + time_interval
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_4D

end module scale_fileio
