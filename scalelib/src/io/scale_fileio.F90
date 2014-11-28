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
  public :: FILEIO_set_coordinates
  public :: FILEIO_set_axes
  public :: FILEIO_read
  public :: FILEIO_write

  interface FILEIO_read
     module procedure FILEIO_read_1D
     module procedure FILEIO_read_2D
     module procedure FILEIO_read_3D
     module procedure FILEIO_read_4D
  end interface FILEIO_read

  interface FILEIO_write
     module procedure FILEIO_write_1D
     module procedure FILEIO_write_2D
     module procedure FILEIO_write_3D
     module procedure FILEIO_write_4D
  end interface FILEIO_write

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILEIO_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FIELIO] / Categ[IO] / Origin[SCALElib]'
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

    ! only for register
    call PROF_rapstart('FILE I NetCDF', 2)
    call PROF_rapend  ('FILE I NetCDF')
    call PROF_rapstart('FILE O NetCDF', 2)
    call PROF_rapend  ('FILE O NetCDF')
    call PROF_rapstart('FILE I ASCII', 2)
    call PROF_rapend  ('FILE I ASCII')
    call PROF_rapstart('FILE O ASCII', 2)
    call PROF_rapend  ('FILE O ASCII')
    call PROF_rapstart('FILE O Interpolation', 2)
    call PROF_rapend  ('FILE O Interpolation')

    return
  end subroutine FILEIO_setup

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
       fid,  &
       dtype )
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
       GRID_LFZ, &
       GRID_LCDZ
    use scale_urban_grid, only: &
       GRID_UCZ, &
       GRID_UFZ, &
       GRID_UCDZ
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: dtype

    character(len=2) :: AXIS_name(2)
    !---------------------------------------------------------------------------

    call FilePutAxis( fid, 'z',   'Z',               'm', 'z',   dtype, GRID_CZ(KS:KE) )
    call FilePutAxis( fid, 'x',   'X',               'm', 'x',   dtype, GRID_CX(ISB:IEB) )
    call FilePutAxis( fid, 'y',   'Y',               'm', 'y',   dtype, GRID_CY(JSB:JEB) )
    call FilePutAxis( fid, 'zh',  'Z (half level)',  'm', 'zh',  dtype, GRID_FZ(KS:KE) )
    call FilePutAxis( fid, 'xh',  'X (half level)',  'm', 'xh',  dtype, GRID_FX(ISB:IEB) )
    call FilePutAxis( fid, 'yh',  'Y (half level)',  'm', 'yh',  dtype, GRID_FY(JSB:JEB) )

    call FilePutAxis( fid, 'lz',  'LZ',              'm', 'lz',  dtype, GRID_LCZ(LKS:LKE) )
    call FilePutAxis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, GRID_LFZ(LKS:LKE) )

    call FilePutAxis( fid, 'uz',  'UZ',              'm', 'uz',  dtype, GRID_UCZ(UKS:UKE) )
    call FilePutAxis( fid, 'uzh', 'UZ (half level)', 'm', 'uzh', dtype, GRID_UFZ(UKS:UKE) )

    call FilePutAxis( fid, 'CZ',  'Atmos Grid Center Position Z', 'm', 'CZ',  dtype, GRID_CZ )
    call FilePutAxis( fid, 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  dtype, GRID_CX )
    call FilePutAxis( fid, 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  dtype, GRID_CY )
    call FilePutAxis( fid, 'FZ',  'Atmos Grid Face Position Z',   'm', 'FZ',  dtype, GRID_FZ )
    call FilePutAxis( fid, 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  dtype, GRID_FX )
    call FilePutAxis( fid, 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  dtype, GRID_FY )

    call FilePutAxis( fid, 'CDZ', 'Grid Cell length Z', 'm', 'CZ',  dtype, GRID_CDZ )
    call FilePutAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, GRID_CDX )
    call FilePutAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, GRID_CDY )
    call FilePutAxis( fid, 'FDZ', 'Grid distance Z',    'm', 'FDZ', dtype, GRID_FDZ )
    call FilePutAxis( fid, 'FDX', 'Grid distance X',    'm', 'FDX', dtype, GRID_FDX )
    call FilePutAxis( fid, 'FDY', 'Grid distance Y',    'm', 'FDY', dtype, GRID_FDY )

    call FilePutAxis( fid, 'LCZ',  'Land Grid Center Position Z',  'm', 'LCZ', dtype, GRID_LCZ )
    call FilePutAxis( fid, 'LFZ',  'Land Grid Face Position Z',    'm', 'LFZ', dtype, GRID_LFZ )
    call FilePutAxis( fid, 'LCDZ', 'Land Grid Cell length Z',      'm', 'LCZ', dtype, GRID_LCZ )

    call FilePutAxis( fid, 'UCZ',  'Urban Grid Center Position Z', 'm', 'UCZ', dtype, GRID_UCZ )
    call FilePutAxis( fid, 'UFZ',  'Urban Grid Face Position Z',   'm', 'UFZ', dtype, GRID_UFZ )
    call FilePutAxis( fid, 'UCDZ', 'Urban Grid Cell length Z',     'm', 'UCZ', dtype, GRID_UCZ )

    call FilePutAxis( fid, 'CBFZ', 'Boundary factor Center Z', '1', 'CZ', dtype, GRID_CBFZ )
    call FilePutAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, GRID_CBFX )
    call FilePutAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, GRID_CBFY )
    call FilePutAxis( fid, 'FBFZ', 'Boundary factor Face Z',   '1', 'CZ', dtype, GRID_FBFZ )
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
    call FileSetTAttr( fid, 'lz',  'positive', 'down' )
    call FileSetTAttr( fid, 'lzh', 'positive', 'down' )
    call FileSetTAttr( fid, 'uz',  'positive', 'down' )
    call FileSetTAttr( fid, 'uzh', 'positive', 'down' )
    call FileSetTAttr( fid, 'LCZ', 'positive', 'down' )
    call FileSetTAttr( fid, 'LFZ', 'positive', 'down' )
    call FileSetTAttr( fid, 'UCZ', 'positive', 'down' )
    call FileSetTAttr( fid, 'UFZ', 'positive', 'down' )

    return
  end subroutine FILEIO_set_axes

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

    call PROF_rapstart('FILE I NetCDF', 2)

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

    call PROF_rapend  ('FILE I NetCDF')

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

    call PROF_rapstart('FILE I NetCDF', 2)

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

    call PROF_rapend  ('FILE I NetCDF')

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
    character(len=*), intent(in)  :: axistype   !< axis type (Z/X/Y)
    integer,          intent(in)  :: step       !< step number

    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    integer               :: dim3_max, dim3_S, dim3_E
    real(RP), allocatable :: var3D(:,:,:)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE I NetCDF', 2)

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

    call PROF_rapend  ('FILE I NetCDF')

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

    call PROF_rapstart('FILE I NetCDF', 2)

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

    call PROF_rapend  ('FILE I NetCDF')

    return
  end subroutine FILEIO_read_4D

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
       datatype,  &
       append    )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FileAddVariable, &
       FileWrite
    use scale_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in)  :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    logical, optional, intent(in) :: append   !< switch whether append existing file or not (default=false)

    integer               :: dtype
    character(len=2)      :: dims(1)
    integer               :: dim1_max, dim1_S, dim1_E
    real(RP), allocatable :: var1D(:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF', 2)

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
                     PRC_master,     & ! [IN]
                     PRC_myrank,     & ! [IN]
                     rankidx,        & ! [IN]
                     append = append ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
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
    call FileWrite( vid, var1D(:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var1D )

    call PROF_rapend  ('FILE O NetCDF')

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
       FileWrite
    use scale_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in)  :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)
    logical, optional, intent(in) :: append   !< switch whether append existing file or not (default=false)
    logical, optional, intent(in) :: nohalo   !< switch whether include halo data or not (default=false)

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
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF', 2)

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
                     PRC_master,     & ! [IN]
                     PRC_myrank,     & ! [IN]
                     rankidx,        & ! [IN]
                     append = append ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
    endif

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
    call FileWrite( vid, var2D(:,:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var2D )

    call PROF_rapend  ('FILE O NetCDF')

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
       FileWrite
    use scale_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    implicit none

    real(RP),          intent(in)  :: var(:,:,:) !< value of the variable
    character(len=*),  intent(in)  :: basename   !< basename of the file
    character(len=*),  intent(in)  :: title      !< title    of the file
    character(len=*),  intent(in)  :: varname    !< name        of the variable
    character(len=*),  intent(in)  :: desc       !< description of the variable
    character(len=*),  intent(in)  :: unit       !< unit        of the variable
    character(len=*),  intent(in)  :: axistype   !< axis type (Z/X/Y)
    character(len=*),  intent(in)  :: datatype   !< data type (REAL8/REAL4/default)
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
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF', 2)

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

    call FileCreate( fid,               & ! [OUT]
                     fileexisted,       & ! [OUT]
                     basename,          & ! [IN]
                     title,             & ! [IN]
                     H_SOURCE,          & ! [IN]
                     H_INSTITUTE,       & ! [IN]
                     PRC_master,        & ! [IN]
                     PRC_myrank,        & ! [IN]
                     rankidx,           & ! [IN]
                     append = append_sw ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FILEIO_set_axes( fid, dtype ) ! [IN]
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
    call FileWrite( vid, var3D(:,:,:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var3D )

    call PROF_rapend  ('FILE O NetCDF')

    return
  end subroutine FILEIO_write_3D

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
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
!    use scale_time, only: &
!       NOWSEC => TIME_NOWDAYSEC
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

    integer :: rankidx(2)
    logical :: append_sw
    logical :: fileexisted
    logical :: addvar
    integer :: fid, vid
    integer :: step
    integer :: i, j, k, n
    logical :: nohalo_
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF', 2)

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

    call FileCreate( fid,               & ! [OUT]
                     fileexisted,       & ! [OUT]
                     basename,          & ! [IN]
                     title,             & ! [IN]
                     H_SOURCE,          & ! [IN]
                     H_INSTITUTE,       & ! [IN]
                     PRC_master,        & ! [IN]
                     PRC_myrank,        & ! [IN]
                     rankidx,           & ! [IN]
                     append = append_sw ) ! [IN]

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
       call FileWrite( vid, var3D(:,:,:), nowtime, nowtime ) ! [IN]
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
          call FileWrite( vid, var3D(:,:,:), nowtime, nowtime ) ! [IN]
          nowtime = nowtime + time_interval
       enddo
    endif

    deallocate( var3D )

    call PROF_rapend  ('FILE O NetCDF')

    return
  end subroutine FILEIO_write_4D

end module scale_fileio
!-------------------------------------------------------------------------------
