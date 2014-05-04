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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FILEIO_setup
  public :: FILEIO_set_coordinates
  public :: FILEIO_read
  public :: FILEIO_write

  interface FILEIO_read
     module procedure FILEIO_read_1D
     module procedure FILEIO_read_2D
     module procedure FILEIO_read_3D
  end interface FILEIO_read

  interface FILEIO_write
     module procedure FILEIO_write_1D
     module procedure FILEIO_write_2D
     module procedure FILEIO_write_3D
  end interface FILEIO_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: set_axes

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: AXIS_LON (:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LONX(:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LAT (:,:) ! [deg]
  real(RP), private, allocatable :: AXIS_LATY(:,:) ! [deg]

  real(RP), private, allocatable :: AXIS_CZ(:,:,:) ! [m]
  real(RP), private, allocatable :: AXIS_FZ(:,:,:) ! [m]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILEIO_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FIELIO]/Categ[IO]'

    allocate( AXIS_LON (IA,JA) )
    allocate( AXIS_LONX(IA,JA) )
    allocate( AXIS_LAT (IA,JA) )
    allocate( AXIS_LATY(IA,JA) )

    allocate( AXIS_CZ  (  KA,IA,JA) )
    allocate( AXIS_FZ  (0:KA,IA,JA) )

    ! only for register
    call PROF_rapstart('FILE I NetCDF')
    call PROF_rapend  ('FILE I NetCDF')
    call PROF_rapstart('FILE O NetCDF')
    call PROF_rapend  ('FILE O NetCDF')
    call PROF_rapstart('FILE I ASCII')
    call PROF_rapend  ('FILE I ASCII')
    call PROF_rapstart('FILE O ASCII')
    call PROF_rapend  ('FILE O ASCII')
    call PROF_rapstart('FILE O Interpolation')
    call PROF_rapend  ('FILE O Interpolation')

    return
  end subroutine FILEIO_setup

  !-----------------------------------------------------------------------------
  !> set latlon and z
  subroutine FILEIO_set_coordinates( &
       LON,  &
       LONX, &
       LAT,  &
       LATY, &
       CZ,   &
       FZ    )
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    real(RP), intent(in) :: LON (IA,JA)
    real(RP), intent(in) :: LONX(IA,JA)
    real(RP), intent(in) :: LAT (IA,JA)
    real(RP), intent(in) :: LATY(IA,JA)
    real(RP), intent(in) :: CZ  (  KA,IA,JA)
    real(RP), intent(in) :: FZ  (0:KA,IA,JA)
    !---------------------------------------------------------------------------

    AXIS_LON (:,:) = LON (:,:) / D2R
    AXIS_LONX(:,:) = LONX(:,:) / D2R
    AXIS_LAT (:,:) = LAT (:,:) / D2R
    AXIS_LATY(:,:) = LATY(:,:) / D2R

    AXIS_CZ(:,:,:) = CZ(:,:,:)
    AXIS_FZ(:,:,:) = FZ(:,:,:)

    return
  end subroutine FILEIO_set_coordinates

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

    call PROF_rapstart('FILE I NetCDF')

    if ( axistype == 'Z' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype == 'X' ) then
       dim1_max = IMAX
       dim1_S   = IS
       dim1_E   = IE
    elseif( axistype == 'Y' ) then
       dim1_max = JMAX
       dim1_S   = JS
       dim1_E   = JE
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

    call PROF_rapstart('FILE I NetCDF')

    if ( axistype == 'XY' ) then
       dim1_max = IMAX
       dim2_max = JMAX
       dim1_S   = IS
       dim1_E   = IE
       dim2_S   = JS
       dim2_E   = JE
    elseif( axistype == 'ZX' ) then
       dim1_max = KMAX
       dim2_max = IMAX
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = IS
       dim2_E   = IE
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

    call PROF_rapstart('FILE I NetCDF')

    if ( axistype == 'ZXY' ) then
       dim1_max = KMAX
       dim2_max = IMAX
       dim3_max = JMAX
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = IS
       dim2_E   = IE
       dim3_S   = JS
       dim3_E   = JE
    else if ( axistype == 'Land' ) then
       dim1_max = LKMAX
       dim2_max = IMAX
       dim3_max = JMAX
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = IS
       dim2_E   = IE
       dim3_S   = JS
       dim3_E   = JE
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

    call PROF_rapstart('FILE O NetCDF')

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
       call set_axes( fid, dtype ) ! [IN]
    endif

    if ( axistype == 'Z' ) then
       dims = (/'z'/)
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype == 'X' ) then
       dims = (/'x'/)
       dim1_max = IMAX
       dim1_S   = IS
       dim1_E   = IE
    elseif( axistype == 'Y' ) then
       dims = (/'y'/)
       dim1_max = JMAX
       dim1_S   = JS
       dim1_E   = JE
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

    real(RP),         intent(in)  :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)
    logical, optional, intent(in) :: append   !< switch whether append existing file or not (default=false)

    integer               :: dtype
    character(len=2)      :: dims(2)
    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    real(RP), allocatable :: var2D(:,:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

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
       call set_axes( fid, dtype ) ! [IN]
    endif

    if ( axistype == 'XY' ) then
       dims = (/'x','y'/)
       dim1_max = IMAX
       dim2_max = JMAX
       dim1_S   = IS
       dim1_E   = IE
       dim2_S   = JS
       dim2_E   = JE
    elseif( axistype == 'ZX' ) then
       dims = (/'z','x'/)
       dim1_max = KMAX
       dim2_max = IMAX
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = IS
       dim2_E   = IE
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype ) ! [IN]

    allocate( var2D(dim1_max,dim2_max) )

    var2D(1:dim1_max,1:dim2_max) = var(dim1_S:dim1_E,dim2_S:dim2_E)
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

    real(RP),         intent(in)  :: var(:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: title      !< title    of the file
    character(len=*), intent(in)  :: varname    !< name        of the variable
    character(len=*), intent(in)  :: desc       !< description of the variable
    character(len=*), intent(in)  :: unit       !< unit        of the variable
    character(len=*), intent(in)  :: axistype   !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype   !< data type (REAL8/REAL4/default)
    logical, optional, intent(in) :: append   !< switch whether append existing file or not (default=false)

    integer               :: dtype
    character(len=2)      :: dims(3)
    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    integer               :: dim3_max, dim3_S, dim3_E
    real(RP), allocatable :: var3D(:,:,:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE O NetCDF')

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
       call set_axes( fid, dtype ) ! [IN]
    endif

    if ( axistype == 'ZXY' ) then
       dims = (/'z','x','y'/)
       dim1_max = KMAX
       dim2_max = IMAX
       dim3_max = JMAX
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = IS
       dim2_E   = IE
       dim3_S   = JS
       dim3_E   = JE
    else if ( axistype == 'Land' ) then
       dims = (/'lz','x ','y '/)
       dim1_max = LKMAX
       dim2_max = IMAX
       dim3_max = JMAX
       dim1_S   = LKS
       dim1_E   = LKE
       dim2_S   = IS
       dim2_E   = IE
       dim3_S   = JS
       dim3_E   = JE
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype ) ! [IN]

    allocate( var3D(dim1_max,dim2_max,dim3_max) )

    var3D(1:dim1_max,1:dim2_max,1:dim3_max) = var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E)
    call FileWrite( vid, var3D(:,:,:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var3D )

    call PROF_rapend  ('FILE O NetCDF')

    return
  end subroutine FILEIO_write_3D

  !-----------------------------------------------------------------------------
  subroutine set_axes( &
       fid,  &
       dtype )
    use gtool_file, only: &
       FilePutAxis,     &
       FilePutAssociatedCoordinates
    use scale_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY, &
       GRID_FZ, &
       GRID_FX, &
       GRID_FY
    use scale_land_grid, only: &
       GRID_LCZ, &
       GRID_LFZ
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: dtype
    !---------------------------------------------------------------------------

    call FilePutAxis( fid, 'z', 'Z', 'm', 'z', dtype, GRID_CZ(KS:KE) )
    call FilePutAxis( fid, 'x', 'X', 'm', 'x', dtype, GRID_CX(IS:IE) )
    call FilePutAxis( fid, 'y', 'Y', 'm', 'y', dtype, GRID_CY(JS:JE) )
    call FilePutAxis( fid, 'zh', 'Z (half level)', 'm', 'zh', dtype, GRID_FZ(KS:KE) )
    call FilePutAxis( fid, 'xh', 'X (half level)', 'm', 'xh', dtype, GRID_FX(IS:IE) )
    call FilePutAxis( fid, 'yh', 'Y (half level)', 'm', 'yh', dtype, GRID_FY(JS:JE) )

    call FilePutAxis( fid, 'lz', 'LZ', 'm', 'lz', dtype, GRID_LCZ(LKS:LKE) )
    call FilePutAxis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, GRID_LFZ(LKS:LKE) )

    call FilePutAssociatedCoordinates( fid, 'lon' , 'longitude'             , 'degrees_east' , &
                                       (/'x ','y '/), dtype, AXIS_LON (IS:IE,JS:JE)            )
    call FilePutAssociatedCoordinates( fid, 'lonh', 'longitude (half level)', 'degrees_east' , &
                                       (/'xh','y '/), dtype, AXIS_LONX(IS:IE,JS:JE)            )
    call FilePutAssociatedCoordinates( fid, 'lat' , 'latitude'              , 'degrees_north', &
                                       (/'x ','y '/), dtype, AXIS_LAT (IS:IE,JS:JE)            )
    call FilePutAssociatedCoordinates( fid, 'lath', 'latitude (half level)' , 'degrees_north', &
                                       (/'x ','yh'/), dtype, AXIS_LATY(IS:IE,JS:JE)            )

    return
  end subroutine set_axes

end module scale_fileio
!-------------------------------------------------------------------------------
