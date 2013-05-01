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
module mod_fileio
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
# include "scale-les.h"
  include "inc_precision.h"
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FILEIO_setup
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
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_SYSCHR), private :: FILEIO_H_SOURCE    = 'SCALE-LES ver. '//VERSION !< for header
  character(len=IO_SYSCHR), private :: FILEIO_H_INSTITUTE = 'AICS/RIKEN'      !< for header

  character(len=1), parameter :: FILEIO_dim_name (3) = (/'z','x','y'/)        !< for axis property
  integer,          parameter :: FILEIO_dim_size (3) = (/KMAX,IMAX,JMAX/)     !< for axis property
  character(len=1), parameter :: FILEIO_dim_desc (3) = (/'Z','X','Y'/)        !< for axis property
  character(len=1), parameter :: FILEIO_dim_unit (3) = (/'m','m','m'/)        !< for axis property
  integer,               save :: FILEIO_dim_dtype(3)                          !< for axis property

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILEIO_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_FILEIO / &
       FILEIO_H_SOURCE, &
       FILEIO_H_INSTITUTE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FIELIO]/Categ[IO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_FILEIO,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_FILEIO. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_FILEIO)

    call TIME_rapstart('FILE I NetCDF')
    call TIME_rapend  ('FILE I NetCDF')
    call TIME_rapstart('FILE O NetCDF')
    call TIME_rapend  ('FILE O NetCDF')
    call TIME_rapstart('FILE I ASCII')
    call TIME_rapend  ('FILE I ASCII')
    call TIME_rapstart('FILE O ASCII')
    call TIME_rapend  ('FILE O ASCII')

    return
  end subroutine FILEIO_setup

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
    use mod_process, only: &
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

    call TIME_rapstart('FILE I NetCDF')

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

    call TIME_rapend  ('FILE I NetCDF')

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
    use mod_process, only: &
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

    call TIME_rapstart('FILE I NetCDF')

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

    call TIME_rapend  ('FILE I NetCDF')

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
    use mod_process, only: &
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

    call TIME_rapstart('FILE I NetCDF')

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
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    allocate( var3D(dim1_max,dim2_max,dim3_max) )

    call FileRead( var3D(:,:,:), basename, varname, step, PRC_myrank )
    var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E) = var3D(1:dim1_max,1:dim2_max,1:dim3_max)

    deallocate( var3D )

    call TIME_rapend  ('FILE I NetCDF')

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
       datatype  )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FilePutAxis,     &
       FileAddVariable, &
       FileWrite
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP),         intent(in)  :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer               :: dtype
    character(len=2)      :: dims(1)
    integer               :: dim1_max, dim1_S, dim1_E
    real(RP), allocatable :: var1D(:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

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
    FILEIO_dim_dtype(:) = dtype

    call FileCreate( fid,                & ! [OUT]
                     fileexisted,        & ! [OUT]
                     basename,           & ! [IN]
                     title,              & ! [IN]
                     FILEIO_H_SOURCE,    & ! [IN]
                     FILEIO_H_INSTITUTE, & ! [IN]
                     FILEIO_dim_name,    & ! [IN]
                     FILEIO_dim_size,    & ! [IN]
                     FILEIO_dim_desc,    & ! [IN]
                     FILEIO_dim_unit,    & ! [IN]
                     FILEIO_dim_dtype,   & ! [IN]
                     PRC_master,         & ! [IN]
                     PRC_myrank,         & ! [IN]
                     rankidx             ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FilePutAxis( fid, 'z', GRID_CZ(KS:KE) )
       call FilePutAxis( fid, 'x', GRID_CX(IS:IE) )
       call FilePutAxis( fid, 'y', GRID_CY(JS:JE) )
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

    call TIME_rapend  ('FILE O NetCDF')

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
       datatype  )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FilePutAxis,     &
       FileAddVariable, &
       FileWrite
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP),         intent(in)  :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer               :: dtype
    character(len=2)      :: dims(2)
    integer               :: dim1_max, dim1_S, dim1_E
    integer               :: dim2_max, dim2_S, dim2_E
    real(RP), allocatable :: var2D(:,:)

    integer :: rankidx(2)
    logical :: fileexisted
    integer :: fid, vid
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O NetCDF')

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
    FILEIO_dim_dtype(:) = dtype

    call FileCreate( fid,                & ! [OUT]
                     fileexisted,        & ! [OUT]
                     basename,           & ! [IN]
                     title,              & ! [IN]
                     FILEIO_H_SOURCE,    & ! [IN]
                     FILEIO_H_INSTITUTE, & ! [IN]
                     FILEIO_dim_name,    & ! [IN]
                     FILEIO_dim_size,    & ! [IN]
                     FILEIO_dim_desc,    & ! [IN]
                     FILEIO_dim_unit,    & ! [IN]
                     FILEIO_dim_dtype,   & ! [IN]
                     PRC_master,         & ! [IN]
                     PRC_myrank,         & ! [IN]
                     rankidx             ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FilePutAxis( fid, 'z', GRID_CZ(KS:KE) )
       call FilePutAxis( fid, 'x', GRID_CX(IS:IE) )
       call FilePutAxis( fid, 'y', GRID_CY(JS:JE) )
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

    call TIME_rapend  ('FILE O NetCDF')

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
       datatype  )
    use gtool_file_h, only: &
       File_REAL8, &
       File_REAL4
    use gtool_file, only: &
       FileCreate,      &
       FilePutAxis,     &
       FileAddVariable, &
       FileWrite
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_2Drank, &
       PRC_MPIstop
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    implicit none

    real(RP),         intent(in)  :: var(:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: title      !< title    of the file
    character(len=*), intent(in)  :: varname    !< name        of the variable
    character(len=*), intent(in)  :: desc       !< description of the variable
    character(len=*), intent(in)  :: unit       !< unit        of the variable
    character(len=*), intent(in)  :: axistype   !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype   !< data type (REAL8/REAL4/default)

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

    call TIME_rapstart('FILE O NetCDF')

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
    FILEIO_dim_dtype(:) = dtype

    call FileCreate( fid,                & ! [OUT]
                     fileexisted,        & ! [OUT]
                     basename,           & ! [IN]
                     title,              & ! [IN]
                     FILEIO_H_SOURCE,    & ! [IN]
                     FILEIO_H_INSTITUTE, & ! [IN]
                     FILEIO_dim_name,    & ! [IN]
                     FILEIO_dim_size,    & ! [IN]
                     FILEIO_dim_desc,    & ! [IN]
                     FILEIO_dim_unit,    & ! [IN]
                     FILEIO_dim_dtype,   & ! [IN]
                     PRC_master,         & ! [IN]
                     PRC_myrank,         & ! [IN]
                     rankidx             ) ! [IN]

    if ( .NOT. fileexisted ) then ! only once
       call FilePutAxis( fid, 'z', GRID_CZ(KS:KE) )
       call FilePutAxis( fid, 'x', GRID_CX(IS:IE) )
       call FilePutAxis( fid, 'y', GRID_CY(JS:JE) )
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
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    call FileAddVariable( vid, fid, varname, desc, unit, dims, dtype ) ! [IN]

    allocate( var3D(dim1_max,dim2_max,dim3_max) )

    var3D(1:dim1_max,1:dim2_max,1:dim3_max) = var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E)
    call FileWrite( vid, var3D(:,:,:), NOWSEC, NOWSEC ) ! [IN]

    deallocate( var3D )

    call TIME_rapend  ('FILE O NetCDF')

    return
  end subroutine FILEIO_write_3D

end module mod_fileio
!-------------------------------------------------------------------------------
