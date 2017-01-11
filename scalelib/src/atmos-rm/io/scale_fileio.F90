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
  public :: FILEIO_check_coordinates
  public :: FILEIO_read
  public :: FILEIO_write
  public :: FILEIO_flush

  public :: FILEIO_create
  public :: FILEIO_open
  public :: FILEIO_def_var
  public :: FILEIO_enddef
  public :: FILEIO_write_var
  public :: FILEIO_close

  interface FILEIO_check_coordinates
     module procedure FILEIO_check_coordinates_name
     module procedure FILEIO_check_coordinates_id
  end interface FILEIO_check_coordinates

  interface FILEIO_read
     module procedure FILEIO_read_1D
     module procedure FILEIO_read_2D
     module procedure FILEIO_read_3D
     module procedure FILEIO_read_4D

     module procedure FILEIO_read_var_1D
     module procedure FILEIO_read_var_2D
     module procedure FILEIO_read_var_3D
     module procedure FILEIO_read_var_4D
  end interface FILEIO_read

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
  private :: closeall
  private :: getCFtunits
  private :: check_1d
  private :: check_2d
  private :: check_3d
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: AXIS_LON  (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONX (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONY (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONXY(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LAT  (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATX (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATY (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATXY(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_HZXY (:,:,:)
  real(RP), private, allocatable :: AXIS_HWXY (:,:,:)

  integer,  private, parameter   :: File_nfile_max = 512                  ! number limit of file
                                                                          ! Keep consistency with "FILE_MAX" in gtool_netcdf.c
  logical,  private              :: File_axes_written(0:File_nfile_max-1) ! whether axes have been written
  !                                                                       ! fid starts from zero so index should start from zero
  logical,  private              :: File_closed(0:File_nfile_max-1)       ! whether file has been closed
  logical,  private              :: File_nozcoord(0:File_nfile_max-1)     ! whether nozcoord is true or false
  integer,  private              :: write_buf_amount = 0                  ! sum of write buffer amounts

  ! global star and count
  integer,  private              :: startXY   (3), countXY   (3)
  integer,  private              :: startZX   (2), countZX   (2)
  integer,  private              :: startZXY  (4), countZXY  (4)
  integer,  private              :: startLAND (3), countLAND (3)
  integer,  private              :: startURBAN(3), countURBAN(3)

  ! MPI element datatype for restart variables
  integer,  private              :: etype

  ! MPI derived datatypes
  integer,  private              :: centerTypeXY
  integer,  private              :: centerTypeZX
  integer,  private              :: centerTypeZXY
  integer,  private              :: centerTypeLAND
  integer,  private              :: centerTypeURBAN

  logical,  private              :: set_coordinates = .false.

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

    allocate( AXIS_LON  (IMAXB,JMAXB) )
    allocate( AXIS_LONX (IMAXB,JMAXB) )
    allocate( AXIS_LONY (IMAXB,JMAXB) )
    allocate( AXIS_LONXY(IMAXB,JMAXB) )
    allocate( AXIS_LAT  (IMAXB,JMAXB) )
    allocate( AXIS_LATX (IMAXB,JMAXB) )
    allocate( AXIS_LATY (IMAXB,JMAXB) )
    allocate( AXIS_LATXY(IMAXB,JMAXB) )

    allocate( AXIS_HZXY (KMAX,IMAXB,JMAXB) )
    allocate( AXIS_HWXY (KMAX,IMAXB,JMAXB) )

    if( IO_AGGREGATE ) call Construct_Derived_Datatype

    File_closed(:) = .true.

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
    deallocate( AXIS_HZXY  )
    deallocate( AXIS_HWXY  )

    call Free_Derived_Datatype

    call closeall

    return
  end subroutine FILEIO_cleanup

  !-----------------------------------------------------------------------------
  !> set latlon and z
  subroutine FILEIO_set_coordinates( &
       LON,   &
       LONX,  &
       LONY,  &
       LONXY, &
       LAT,   &
       LATX,  &
       LATY,  &
       LATXY, &
       CZ,    &
       FZ     )
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    real(RP), intent(in) :: LON  (IA,JA)
    real(RP), intent(in) :: LONX (IA,JA)
    real(RP), intent(in) :: LONY (IA,JA)
    real(RP), intent(in) :: LONXY(IA,JA)
    real(RP), intent(in) :: LAT  (IA,JA)
    real(RP), intent(in) :: LATX (IA,JA)
    real(RP), intent(in) :: LATY (IA,JA)
    real(RP), intent(in) :: LATXY(IA,JA)
    real(RP), intent(in) :: CZ   (  KA,IA,JA)
    real(RP), intent(in) :: FZ   (0:KA,IA,JA)
    !---------------------------------------------------------------------------

    AXIS_LON  (:,:)   = LON  (ISB:IEB,JSB:JEB) / D2R
    AXIS_LONX (:,:)   = LONX (ISB:IEB,JSB:JEB) / D2R
    AXIS_LONY (:,:)   = LONY (ISB:IEB,JSB:JEB) / D2R
    AXIS_LONXY(:,:)   = LONXY(ISB:IEB,JSB:JEB) / D2R
    AXIS_LAT  (:,:)   = LAT  (ISB:IEB,JSB:JEB) / D2R
    AXIS_LATX (:,:)   = LATX (ISB:IEB,JSB:JEB) / D2R
    AXIS_LATY (:,:)   = LATY (ISB:IEB,JSB:JEB) / D2R
    AXIS_LATXY(:,:)   = LATXY(ISB:IEB,JSB:JEB) / D2R

    AXIS_HZXY (:,:,:) = CZ(KS:KE,ISB:IEB,JSB:JEB)
    AXIS_HZXY (:,:,:) = FZ(KS:KE,ISB:IEB,JSB:JEB)

    set_coordinates = .true.

    return
  end subroutine FILEIO_set_coordinates

  !-----------------------------------------------------------------------------
  !> check coordinates in the file
  subroutine FILEIO_check_coordinates_name( &
       basename, &
       atmos,    &
       land,     &
       urban     )
    implicit none

    character(len=*), intent(in) :: basename        !< basename of the file
    logical,          intent(in), optional :: atmos !< check atmospheric coordinates
    logical,          intent(in), optional :: land  !< check land coordinates
    logical,          intent(in), optional :: urban !< check urban coordinates

    integer :: fid
    !---------------------------------------------------------------------------

    call FILEIO_open( fid,     & ! [OUT]
                      basename ) ! [IN]

    call FILEIO_check_coordinates_id( fid,               & ! [IN]
                                      atmos, land, urban ) ! [IN]

    return
  end subroutine FILEIO_check_coordinates_name

  !-----------------------------------------------------------------------------
  !> check coordinates in the file
  subroutine FILEIO_check_coordinates_id( &
       fid,   &
       atmos, &
       land,  &
       urban  )
    use scale_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    use scale_land_grid, only: &
       GRID_LCZ
    use scale_urban_grid, only: &
       GRID_UCZ
    implicit none

    integer, intent(in) :: fid
    logical, intent(in), optional :: atmos !< check atmospheric coordinates
    logical, intent(in), optional :: land  !< check land coordinates
    logical, intent(in), optional :: urban !< check urban coordinates

    logical  :: atmos_ = .false.
    logical  :: land_  = .false.
    logical  :: urban_ = .false.

    real(RP) :: buffer_z  (KA)
    real(RP) :: buffer_x  (IA)
    real(RP) :: buffer_y  (JA)
    real(RP) :: buffer_xy (IA,JA)
    real(RP) :: buffer_zxy(KA,IA,JA)
    real(RP) :: buffer_l  (LKMAX)
    real(RP) :: buffer_u  (UKMAX)
    !---------------------------------------------------------------------------

    if( present(atmos) ) atmos_ = atmos
    if( present(land ) ) land_  = land
    if( present(urban) ) urban_ = urban

    call FILEIO_read_var_1D( buffer_x, fid, 'x',  'X', 1 )
    call FILEIO_read_var_1D( buffer_y, fid, 'y',  'Y', 1 )
    call FILEIO_flush( fid )
    call check_1d( GRID_CX(ISB:IEB), buffer_x(ISB:IEB), 'x' )
    call check_1d( GRID_CY(JSB:JEB), buffer_y(JSB:JEB), 'y' )

    if ( set_coordinates ) then
       call FILEIO_read_var_2D( buffer_xy, fid, 'lon', 'XY', 1 )
       call FILEIO_flush( fid )
       call check_2d( AXIS_LON, buffer_xy(ISB:IEB,JSB:JEB), 'lon' )

       call FILEIO_read_var_2D( buffer_xy, fid, 'lat', 'XY', 1 )
       call FILEIO_flush( fid )
       call check_2d( AXIS_LAT, buffer_xy(ISB:IEB,JSB:JEB), 'lat' )
    end if

    if ( atmos_ ) then
       call FILEIO_read_var_1D( buffer_z,   fid, 'z',      'Z',   1 )
       call FILEIO_read_var_3D( buffer_zxy, fid, 'height', 'ZXY', 1 )
       call FILEIO_flush( fid )
       call check_1d( GRID_CZ(KS:KE), buffer_z(KS:KE), 'z' )
       call check_3d( AXIS_HZXY, buffer_zxy(KS:KE,ISB:IEB,JSB:JEB), 'height' )
    end if

    if ( land_ ) then
       call FILEIO_read_var_1D( buffer_l, fid, 'lz',  'LZ', 1 )
       call FILEIO_flush( fid )
       call check_1d( GRID_LCZ(LKS:LKE), buffer_l(LKS:LKE), 'lz' )
    end if

    if ( urban_ ) then
       call FILEIO_read_var_1D( buffer_u, fid, 'uz',  'UZ', 1 )
       call FILEIO_flush( fid )
       call check_1d( GRID_UCZ(UKS:UKE), buffer_u(UKS:UKE), 'uz' )
    end if

    return
  end subroutine FILEIO_check_coordinates_id

  !-----------------------------------------------------------------------------
  !> construct MPI derived datatypes for read buffers
  subroutine Construct_Derived_Datatype
    use mpi
    use scale_process, only: &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X,  &
       PRC_NUM_Y
    implicit none

    integer :: err, order
    integer :: sizes(3), subsizes(3), sub_off(3)
    !---------------------------------------------------------------------------

    order           = MPI_ORDER_FORTRAN

    centerTypeXY    = MPI_DATATYPE_NULL
    centerTypeZX    = MPI_DATATYPE_NULL
    centerTypeZXY   = MPI_DATATYPE_NULL
    centerTypeLAND  = MPI_DATATYPE_NULL
    centerTypeURBAN = MPI_DATATYPE_NULL

    etype           = MPI_FLOAT

    if( RP == 8 ) etype = MPI_DOUBLE_PRECISION

    ! for axistype == 'XY'
    startXY(1)    = IS_inG - IHALO
    startXY(2)    = JS_inG - JHALO
    countXY(1)    = IA
    countXY(2)    = JA
    ! for axistype == 'ZXY'
    startZXY(1)   = 1
    startZXY(2:3) = startXY(1:2)
    countZXY(1)   = KMAX
    countZXY(2:3) = countXY(1:2)
    ! construct MPI subarray data type
    sizes(1)      = KA
    sizes(2)      = IA
    sizes(3)      = JA
    subsizes(1)   = KMAX
    subsizes(2)   = IA
    subsizes(3)   = JA
    sub_off(1)    = KS - 1 ! MPI start index starts with 0
    sub_off(2)    = 0
    sub_off(3)    = 0

    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeZXY, err)
    call MPI_Type_commit(centerTypeZXY, err)

    ! for axistype == 'Land'
    startLAND(1)   = 1
    startLAND(2:3) = startXY(1:2)
    countLAND(1)   = LKMAX
    countLAND(2:3) = countXY(1:2)
    ! construct MPI subarray data type
    sizes(1)       = LKMAX
    subsizes(1)    = LKMAX
    sub_off(1)     = LKS - 1 ! MPI start index starts with 0

    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeLAND, err)
    call MPI_Type_commit(centerTypeLAND, err)

    ! for axistype == 'URBAN'
    startURBAN(1)   = 1
    startURBAN(2:3) = startXY(1:2)
    countURBAN(1)   = UKMAX
    countURBAN(2:3) = countXY(1:2)
    ! construct MPI subarray data type
    sizes(1)        = UKMAX
    subsizes(1)     = UKMAX
    sub_off(1)      = UKS - 1 ! MPI start index starts with 0

    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeURBAN, err)
    call MPI_Type_commit(centerTypeURBAN, err)

    ! for axistype == 'ZX'
    startZX(1)  = KHALO+1
    startZX(2)  = IS_inG - IHALO
    countZX(1)  = KHALO
    countZX(2)  = IA
    ! construct MPI subarray data type
    sizes(1)    = KA
    sizes(2)    = IA
    subsizes(1) = KMAX
    subsizes(2) = IMAXB
    sub_off(1)  = KHALO   ! MPI start index starts with 0
    sub_off(2)  = ISB - 1 ! MPI start index starts with 0

    call MPI_Type_create_subarray(2, sizes, subsizes, sub_off, order, etype, centerTypeZX, err)
    call MPI_Type_commit(centerTypeZX, err)

    return
  end subroutine Construct_Derived_Datatype

  !-----------------------------------------------------------------------------
  !> free MPI derived datatypes
  subroutine Free_Derived_Datatype
    use mpi
    implicit none

    integer :: err
    !---------------------------------------------------------------------------

    if( .NOT. IO_AGGREGATE ) return

    if( centerTypeXY    /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeXY,    err)
    if( centerTypeZX    /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZX,    err)
    if( centerTypeZXY   /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZXY,   err)
    if( centerTypeLAND  /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeLAND,  err)
    if( centerTypeURBAN /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeURBAN, err)

    return
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
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer :: fid
    !---------------------------------------------------------------------------

    call FILEIO_open( fid,      & ! [OUT]
                      basename  ) ! [IN]

    call FILEIO_read_var_1D( var(:),                      & ! [OUT]
                             fid, varname, axistype, step ) ! [IN]

    call FILEIO_close( fid )

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
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer :: fid
    !---------------------------------------------------------------------------

    call FILEIO_open( fid,      & ! [OUT]
                      basename  ) ! [IN]

    call FILEIO_read_var_2D( var(:,:),                    & ! [OUT]
                             fid, varname, axistype, step ) ! [IN]

    call FILEIO_close( fid )

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
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: axistype   !< axis type (Z/X/Y/T)
    integer,          intent(in)  :: step       !< step number

    integer :: fid
    !---------------------------------------------------------------------------

    call FILEIO_open( fid,      & ! [OUT]
                      basename  ) ! [IN]

    call FILEIO_read_var_3D( var(:,:,:),                  & ! [OUT]
                             fid, varname, axistype, step ) ! [IN]

    call FILEIO_close( fid )

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
       PRC_MPIstop
    implicit none

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename     !< basename of the file
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: axistype     !< axis type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number

    integer :: fid
    !---------------------------------------------------------------------------

    call FILEIO_open( fid,      & ! [OUT]
                      basename  ) ! [IN]

    call FILEIO_read_var_4D( var(:,:,:,:),                & ! [OUT]
                             fid, varname, axistype, step ) ! [IN]

    call FILEIO_close( fid )

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
    use mpi
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer :: dim1_S, dim1_E
    integer :: start(1)   ! start offset of globale variable
    integer :: count(1)   ! request length to the globale variable
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 1D var: ', trim(varname)

    if ( IO_AGGREGATE ) then
       ! read data and halos into the local buffer
       if    ( axistype == 'Z' ) then
          start(1) = 1
          count(1) = KMAX
          call FileRead( var(KS:KE), fid, varname, step, PRC_myrank,        &
                         ntypes=KMAX, dtype=etype, start=start, count=count )
       elseif( axistype == 'LZ' ) then
          start(1) = 1
          count(1) = LKMAX
          call FileRead( var, fid, varname, step, PRC_myrank,                &
                         ntypes=LKMAX, dtype=etype, start=start, count=count )
       elseif( axistype == 'UZ' ) then
          start(1) = 1
          count(1) = UKMAX
          call FileRead( var, fid, varname, step, PRC_myrank,                &
                         ntypes=LKMAX, dtype=etype, start=start, count=count )
       elseif( axistype == 'X' .OR. axistype == 'CX' ) then
          start(1) = IS_inG - IHALO
          count(1) = IA
          call FileRead( var, fid, varname, step, PRC_myrank,             &
                         ntypes=IA, dtype=etype, start=start, count=count )
       elseif( axistype == 'Y' .OR. axistype == 'CY' ) then
          start(1) = JS_inG - JHALO
          count(1) = JA
          call FileRead( var, fid, varname, step, PRC_myrank,             &
                         ntypes=JA, dtype=etype, start=start, count=count )
       else
          write(*,*) 'xxx unsupported axis type. Check!: ', trim(axistype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    else
       if    ( axistype == 'Z' ) then
          dim1_S   = KS
          dim1_E   = KE
       elseif( axistype == 'LZ' ) then
          dim1_S   = 1
          dim1_E   = LKMAX
       elseif( axistype == 'UZ' ) then
          dim1_S   = 1
          dim1_E   = UKMAX
       elseif( axistype == 'X' ) then
          dim1_S   = ISB
          dim1_E   = IEB
       elseif( axistype == 'CX' ) then
          dim1_S   = 1
          dim1_E   = IA
       elseif( axistype == 'Y' ) then
          dim1_S   = JSB
          dim1_E   = JEB
       elseif( axistype == 'CY' ) then
          dim1_S   = 1
          dim1_E   = JA
       else
          write(*,*) 'xxx unsupported axis type. Check!: ', trim(axistype), ' item:',trim(varname)
          call PRC_MPIstop
       endif

       call FileRead( var(dim1_S:dim1_E), fid, varname, step, PRC_myrank )
    end if

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
    use mpi
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 2D var: ', trim(varname)

    if ( IO_AGGREGATE ) then
       ! read data and halos into the local buffer
       if    ( axistype == 'XY' ) then
          call FileRead( var, fid, varname, step, PRC_myrank,                    &
                         ntypes=IA*JA, dtype=etype, start=startXY, count=countXY )
       elseif( axistype == 'ZX' ) then
          ! Because KHALO is not saved in files, we use centerTypeZX, an MPI
          ! derived datatype to describe the layout of local read buffer
          call FileRead( var, fid, varname, step, PRC_myrank,                       &
                         ntypes=1, dtype=centerTypeZX, start=startZX, count=countZX )
       else
          write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    else
       if    ( axistype == 'XY' ) then
          dim1_S   = ISB
          dim1_E   = IEB
          dim2_S   = JSB
          dim2_E   = JEB
       elseif( axistype == 'ZX' ) then
          dim1_S   = KS
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
       else
          write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
          call PRC_MPIstop
       endif

       call FileRead( var(dim1_S:dim1_E,dim2_S:dim2_E), fid, varname, step, PRC_myrank )
    end if

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

    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    integer :: dim3_S, dim3_E
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 3D var: ', trim(varname)

    if ( IO_AGGREGATE ) then
       ! read data and halos into the local buffer
       ! Because KHALO is not saved in files, we use mpi derived datatypes to
       ! describe the layout of local read buffer
       if    ( axistype == 'ZXY' ) then
          call FileRead( var, fid, varname, step, PRC_myrank,                          &
                         ntypes=1, dtype=centerTypeZXY, start=startZXY, count=countZXY )
       elseif( axistype == 'XYT' ) then
          startXY(3) = 1
          countXY(3) = step
          call FileRead( var, fid, varname, step, PRC_myrank,                         &
                         ntypes=step*IA*JA, dtype=etype, start=startXY, count=countXY )
       elseif( axistype == 'Land' ) then
          call FileRead( var, fid, varname, step, PRC_myrank,                             &
                         ntypes=1, dtype=centerTypeLAND, start=startLAND, count=countLAND )
       elseif( axistype == 'Urban' ) then
          call FileRead( var, fid, varname, step, PRC_myrank,                                &
                         ntypes=1, dtype=centerTypeURBAN, start=startURBAN, count=countURBAN )
       else
          write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    else
       if    ( axistype == 'ZXY' ) then
          dim1_S   = KS
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       elseif( axistype == 'XYT' ) then
          dim1_S   = ISB
          dim1_E   = IEB
          dim2_S   = JSB
          dim2_E   = JEB
          dim3_S   = 1
          dim3_E   = step
       elseif( axistype == 'Land' ) then
          dim1_S   = LKS
          dim1_E   = LKE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       elseif( axistype == 'Urban' ) then
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

       call FileRead( var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                      fid, varname, step, PRC_myrank                  )
    end if

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

    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    integer :: dim3_S, dim3_E
    integer :: dim4_S, dim4_E
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,A,A15)') '*** Read 4D var: ', trim(varname)

    if ( IO_AGGREGATE ) then
       ! read data and halos into the local buffer
       if ( axistype == 'ZXYT' ) then
          startZXY(4) = 1
          countZXY(4) = step
          call FileRead( var, fid, varname, step, PRC_myrank,                             &
                         ntypes=step, dtype=centerTypeZXY, start=startZXY, count=countZXY )
       else
          write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    else
       if ( axistype == 'ZXYT' ) then
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

       call FileRead( var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,dim4_S:dim4_E), &
                      fid, varname, step, PRC_myrank                                )
    end if

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
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in) :: var(:)   !< value of the variable
    character(len=*), intent(in) :: basename !< basename of the file
    character(len=*), intent(in) :: title    !< title    of the file
    character(len=*), intent(in) :: varname  !< name        of the variable
    character(len=*), intent(in) :: desc     !< description of the variable
    character(len=*), intent(in) :: unit     !< unit        of the variable
    character(len=*), intent(in) :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in) :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6) !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec  !< subsec of the time
    logical,          intent(in), optional :: append  !< switch whether append existing file or not (default=false)

    integer :: fid, vid
    !---------------------------------------------------------------------------

    call FILEIO_create( fid,                                            & ! [OUT]
                        basename, title, datatype, date, subsec, append )

    call FILEIO_def_var( fid, vid, varname, desc, unit, axistype, datatype )

    call FILEIO_enddef( fid )

    call FILEIO_write_var_1D( fid, vid, var, varname, axistype )

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
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in) :: var(:,:) !< value of the variable
    character(len=*), intent(in) :: basename !< basename of the file
    character(len=*), intent(in) :: title    !< title    of the file
    character(len=*), intent(in) :: varname  !< name        of the variable
    character(len=*), intent(in) :: desc     !< description of the variable
    character(len=*), intent(in) :: unit     !< unit        of the variable
    character(len=*), intent(in) :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in) :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6)  !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec   !< subsec of the time
    logical,          intent(in), optional :: append   !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: nohalo   !< switch whether include halo data or not (default=false)
    logical,          intent(in), optional :: nozcoord !< switch whether include zcoordinate or not (default=false)

    integer :: fid, vid
    !---------------------------------------------------------------------------

    call FILEIO_create( fid,                                                      & ! [OUT]
                        basename, title, datatype, date, subsec, append, nozcoord )

    call FILEIO_def_var( fid, vid, varname, desc, unit, axistype, datatype )

    call FILEIO_enddef( fid )

    call FILEIO_write_var_2D( fid, vid, var, varname, axistype, nohalo )

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
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: basename   !< basename of the file
    character(len=*), intent(in) :: title      !< title    of the file
    character(len=*), intent(in) :: varname    !< name        of the variable
    character(len=*), intent(in) :: desc       !< description of the variable
    character(len=*), intent(in) :: unit       !< unit        of the variable
    character(len=*), intent(in) :: axistype   !< axis type (Z/X/Y)
    character(len=*), intent(in) :: datatype   !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6) !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec  !< subsec of the time
    logical,          intent(in), optional :: append  !< append existing (closed) file?
    logical,          intent(in), optional :: nohalo  !< include halo data?

    integer :: fid, vid
    !---------------------------------------------------------------------------

    call FILEIO_create( fid,                                            & ! [OUT]
                        basename, title, datatype, date, subsec, append )

    call FILEIO_def_var( fid, vid, varname, desc, unit, axistype, datatype )

    call FILEIO_enddef( fid )

    call FILEIO_write_var_3D( fid, vid, var, varname, axistype, nohalo )

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
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: basename   !< basename of the file
    character(len=*), intent(in) :: title      !< title    of the file
    character(len=*), intent(in) :: varname    !< name        of the variable
    character(len=*), intent(in) :: desc       !< description of the variable
    character(len=*), intent(in) :: unit       !< unit        of the variable
    character(len=*), intent(in) :: axistype   !< axis type (X/Y/Time)
    character(len=*), intent(in) :: datatype   !< data type (REAL8/REAL4/default)
    real(RP),         intent(in) :: timeintv   !< time interval [sec]
    integer ,         intent(in) :: tsince(6)  !< start time

    logical,          intent(in), optional :: append   !< append existing (closed) file?
    integer,          intent(in), optional :: timetarg !< target timestep (optional)
    logical,          intent(in), optional :: nohalo   !< include halo data?

    integer  :: fid, vid
    integer  :: nsteps

    intrinsic :: size
    !---------------------------------------------------------------------------

    call FILEIO_create( fid,                                             & ! [OUT]
                        basename, title, datatype, tsince, append=append )

    if ( present(timetarg) ) then
       nsteps = 1
    else
       nsteps = size(var,3)
    end if
    call FILEIO_def_var( fid, vid, varname, desc, unit, axistype, datatype, timeintv, nsteps )

    call FILEIO_enddef( fid )

    call FILEIO_write_var_3D_t( fid, vid, var, varname, axistype, timeintv, timetarg, nohalo )

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
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    implicit none

    real(RP),         intent(in) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in) :: basename     !< basename of the file
    character(len=*), intent(in) :: title        !< title    of the file
    character(len=*), intent(in) :: varname      !< name        of the variable
    character(len=*), intent(in) :: desc         !< description of the variable
    character(len=*), intent(in) :: unit         !< unit        of the variable
    character(len=*), intent(in) :: axistype     !< axis type (Z/X/Y/Time)
    character(len=*), intent(in) :: datatype     !< data type (REAL8/REAL4/default)
    real(RP),         intent(in) :: timeintv     !< time interval [sec]
    integer,          intent(in) :: tsince(6)    !< start time

    logical,          intent(in), optional :: append   !< append existing (closed) file?
    integer,          intent(in), optional :: timetarg !< target timestep (optional)
    logical,          intent(in), optional :: nohalo   !< include halo data?

    integer  :: fid, vid
    integer  :: nsteps

    intrinsic :: size
    !---------------------------------------------------------------------------

    call FILEIO_create( fid,                                             & ! [OUT]
                        basename, title, datatype, tsince, append=append )

    if ( present(timetarg) ) then
       nsteps = 1
    else
       nsteps = size(var,3)
    end if
    call FILEIO_def_var( fid, vid, varname, desc, unit, axistype, datatype, timeintv, nsteps )

    call FILEIO_enddef( fid )

    call FILEIO_write_var_4D( fid, vid, var, varname, axistype, timeintv, timetarg, nohalo )

    return
  end subroutine FILEIO_write_4D

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
    use mpi, only : MPI_COMM_NULL
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of the file

    integer :: comm
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( IO_AGGREGATE ) then  ! user input parameter indicates to do PnetCDF I/O
       comm = PRC_LOCAL_COMM_WORLD
    else
       comm = MPI_COMM_NULL
    end if

    call FileOpen( fid,                & ! [OUT]
                   basename,           & ! [IN]
                   File_FREAD,         & ! [IN]
                   comm = comm,        & ! [IN]
                   myrank = PRC_myrank ) ! [IN]

    File_closed(fid) = .false.

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
    use mpi, only: &
       MPI_COMM_NULL
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
       PRC_2Drank,     &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS,   &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6)  !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec   !< subsec of the time
    logical,          intent(in), optional :: append   !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: nozcoord !< switch whether include zcoordinate or not (default=false)

    integer           :: rankidx(2)
    integer           :: dtype
    logical           :: append_sw
    character(len=34) :: tunits
    integer           :: comm
    logical           :: fileexisted
    character(len=8)  :: logical_str
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    ! dtype is used to define the data type of axis variables in file
    if    ( datatype == 'REAL8' ) then
       dtype = File_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
    else
       if    ( RP == 8 ) then
          dtype = File_REAL8
       elseif( RP == 4 ) then
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
       call getCFtunits( tunits, date )
    else
       tunits = 'seconds'
    end if

    if ( IO_AGGREGATE ) then  ! user input parameter indicates to do PnetCDF I/O
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
       File_axes_written(fid) = .false.  ! indicating axes have not been written yet
       if ( present( nozcoord ) ) then
          File_nozcoord(fid) = nozcoord
       else
          File_nozcoord(fid) = .false.
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
       if(PRC_PERIODIC_X .AND. .NOT. IO_AGGREGATE) logical_str = "true"
       call FileSetGlobalAttribute( fid, "PRC_PERIODIC_X",  trim(logical_str) )
       logical_str = "false"
       if(PRC_PERIODIC_Y .AND. .NOT. IO_AGGREGATE) logical_str = "true"
       call FileSetGlobalAttribute( fid, "PRC_PERIODIC_Y",  trim(logical_str) )

       File_closed(fid) = .false.
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF file define mode
  subroutine FILEIO_enddef( &
       fid )
    use gtool_file, only: &
       FileEndDef, &
       FileFlush,  &
       FileAttachBuffer
    implicit none

    integer, intent(in) :: fid  !< file ID
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call FileEndDef( fid ) ! [IN]

    ! If this enddef is called the first time, write axis variables
    if ( .NOT. File_axes_written(fid) ) then
       call FILEIO_write_axes( fid, File_nozcoord(fid) )
       if ( IO_AGGREGATE ) then
          call FileFlush( fid )
          ! Tell PnetCDF library to use a buffer of size write_buf_amount to
          ! aggregate write requests to be post in FILEIO_write_var
          call FileAttachBuffer( fid, write_buf_amount )
       end if
       File_axes_written(fid) = .true.
    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_enddef

  !-----------------------------------------------------------------------------
  !> Flush all pending requests to a netCDF file (PnetCDF only)
  subroutine FILEIO_flush( &
       fid )
    use gtool_file, only: &
       FileFlush
    implicit none

    integer, intent(in) :: fid  !< file ID
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( IO_AGGREGATE ) then
       call FileFlush( fid ) ! flush all pending read/write requests
    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_flush

  !-----------------------------------------------------------------------------
  !> Close a netCDF file
  subroutine FILEIO_close( &
       fid )
    use gtool_file, only: &
       FileClose, &
       FileFlush, &
       FileDetachBuffer
    implicit none

    integer, intent(in) :: fid  !< file ID

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( .NOT. File_closed(fid) ) then

       if ( IO_AGGREGATE ) then
          call FileFlush( fid )        ! flush all pending read/write requests
          if ( write_buf_amount > 0 ) then
             call FileDetachBuffer( fid ) ! detach PnetCDF aggregation buffer
             write_buf_amount = 0         ! reset write request amount
          end if
       end if

       call FileClose( fid ) ! [IN]

       File_closed(fid) = .true.

    end if

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

    character(len=2) :: AXIS_name(3)
    logical          :: xy_
    !---------------------------------------------------------------------------

    if ( present(xy) ) then
       xy_ = xy
    else
       xy_ = .false.
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'z',   'Z',               'm', 'z',   dtype, KMAX )
    end if
    if ( .NOT. IO_AGGREGATE ) then
       call FileDefAxis( fid, 'x',   'X',               'm', 'x',   dtype, IMAXB )
       call FileDefAxis( fid, 'y',   'Y',               'm', 'y',   dtype, JMAXB )
    else
       call FileDefAxis( fid, 'x',   'X',               'm', 'x',   dtype, IAG )
       call FileDefAxis( fid, 'y',   'Y',               'm', 'y',   dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'zh',  'Z (half level)',  'm', 'zh',  dtype, KMAX )
    end if
    if ( .NOT. IO_AGGREGATE ) then
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
    if ( .NOT. IO_AGGREGATE ) then
       call FileDefAxis( fid, 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  dtype, IA )
       call FileDefAxis( fid, 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  dtype, JA )
    else
       call FileDefAxis( fid, 'CX',  'Atmos Grid Center Position X', 'm', 'CX',  dtype, IAG )
       call FileDefAxis( fid, 'CY',  'Atmos Grid Center Position Y', 'm', 'CY',  dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'FZ',  'Atmos Grid Face Position Z',   'm', 'FZ',  dtype, KA+1 )
    end if
    if ( .NOT. IO_AGGREGATE ) then
       call FileDefAxis( fid, 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  dtype, IA+1 )
       call FileDefAxis( fid, 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  dtype, JA+1 )
    else
       call FileDefAxis( fid, 'FX',  'Atmos Grid Face Position X',   'm', 'FX',  dtype, IAG+1 )
       call FileDefAxis( fid, 'FY',  'Atmos Grid Face Position Y',   'm', 'FY',  dtype, JAG+1 )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'CDZ', 'Grid Cell length Z', 'm', 'CZ',  dtype, KA )
    end if
    if ( .NOT. IO_AGGREGATE ) then
       call FileDefAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, IA )
       call FileDefAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, JA )
    else
       call FileDefAxis( fid, 'CDX', 'Grid Cell length X', 'm', 'CX',  dtype, IAG )
       call FileDefAxis( fid, 'CDY', 'Grid Cell length Y', 'm', 'CY',  dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'FDZ', 'Grid distance Z',    'm', 'FDZ', dtype, KA-1 )
    end if
    if ( .NOT. IO_AGGREGATE ) then
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
    if ( .NOT. IO_AGGREGATE ) then
       call FileDefAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, IA )
       call FileDefAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, JA )
    else
       call FileDefAxis( fid, 'CBFX', 'Boundary factor Center X', '1', 'CX', dtype, IAG )
       call FileDefAxis( fid, 'CBFY', 'Boundary factor Center Y', '1', 'CY', dtype, JAG )
    end if

    if ( .NOT. xy_ ) then
       call FileDefAxis( fid, 'FBFZ', 'Boundary factor Face Z',   '1', 'CZ', dtype, KA )
    end if
    if ( .NOT. IO_AGGREGATE ) then
       call FileDefAxis( fid, 'FBFX', 'Boundary factor Face X',   '1', 'CX', dtype, IA )
       call FileDefAxis( fid, 'FBFY', 'Boundary factor Face Y',   '1', 'CY', dtype, JA )
    else
       call FileDefAxis( fid, 'FBFX', 'Boundary factor Face X',   '1', 'CX', dtype, IAG )
       call FileDefAxis( fid, 'FBFY', 'Boundary factor Face Y',   '1', 'CY', dtype, JAG )
    end if

    ! TODO: skip 8 axes below when IO_AGGREGATE is true, as all axes are now global
    call FileDefAxis( fid, 'CXG', 'Grid Center Position X (global)', 'm', 'CXG', dtype, IAG )
    call FileDefAxis( fid, 'CYG', 'Grid Center Position Y (global)', 'm', 'CYG', dtype, JAG )
    call FileDefAxis( fid, 'FXG', 'Grid Face Position X (global)',   'm', 'FXG', dtype, IAG+1 )
    call FileDefAxis( fid, 'FYG', 'Grid Face Position Y (global)',   'm', 'FYG', dtype, JAG+1 )

    call FileDefAxis( fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', dtype, IAG )
    call FileDefAxis( fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', dtype, JAG )
    call FileDefAxis( fid, 'FBFXG', 'Boundary factor Face X (global)',   '1', 'CXG', dtype, IAG )
    call FileDefAxis( fid, 'FBFYG', 'Boundary factor Face Y (global)',   '1', 'CYG', dtype, JAG )

    ! associate coordinates
    AXIS_name(1:2) = (/'x ','y '/)
    call FileDefAssociatedCoordinates( fid, 'lon' , 'longitude',                   &
                                       'degrees_east' , AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'xh','y '/)
    call FileDefAssociatedCoordinates( fid, 'lon_uy', 'longitude (half level uy)', &
                                       'degrees_east' , AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'x ','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lon_xv', 'longitude (half level xv)', &
                                       'degrees_east' , AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'xh','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lon_uv', 'longitude (half level uv)', &
                                       'degrees_east' , AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'x ','y '/)
    call FileDefAssociatedCoordinates( fid, 'lat' , 'latitude',                    &
                                       'degrees_north', AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'xh','y '/)
    call FileDefAssociatedCoordinates( fid, 'lat_uy', 'latitude (half level uy)',  &
                                       'degrees_north', AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'x ','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lat_xv', 'latitude (half level xv)',  &
                                       'degrees_north', AXIS_name(1:2), dtype      )
    AXIS_name(1:2) = (/'xh','yh'/)
    call FileDefAssociatedCoordinates( fid, 'lat_uv', 'latitude (half level uv)',  &
                                       'degrees_north', AXIS_name(1:2), dtype      )

    if ( .NOT. xy_ ) then
       AXIS_name = (/'z', 'x', 'y'/)
       call FileDefAssociatedCoordinates( fid, 'height',     'height above ground level', &
                                          'm', AXIS_name(1:3), dtype                      )
       AXIS_name = (/'zh', 'x ', 'y '/)
       call FileDefAssociatedCoordinates( fid, 'height_wxy', 'height above ground level (half level wxy)', &
                                          'm', AXIS_name(1:3), dtype                                       )
    end if

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
    use scale_process, only: &
       PRC_myrank
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
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
    implicit none

    integer, intent(in) :: fid
    logical, intent(in), optional :: xy

    logical :: xy_
    logical :: put_x
    logical :: put_y
    logical :: put_z

    integer :: rankidx(2)
    integer :: start(3)
    integer :: XSB, XEB, YSB, YEB
    !---------------------------------------------------------------------------

    if ( present(xy) ) then
       xy_ = xy
    else
       xy_ = .false.
    end if

    if ( IO_AGGREGATE ) then
       rankidx(1) = PRC_2Drank(PRC_myrank,1)
       rankidx(2) = PRC_2Drank(PRC_myrank,2)

       ! construct indices independent from PRC_PERIODIC_X/Y
       XSB = 1 + IHALO
       if( rankidx(1) == 0 ) XSB = 1
       XEB = IMAX + IHALO
       if( rankidx(1) == PRC_NUM_X-1 ) XEB = IA

       YSB = 1 + JHALO
       if( rankidx(2) == 0 ) YSB = 1
       YEB = JMAX + JHALO
       if( rankidx(2) == PRC_NUM_Y-1 ) YEB = JA
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

       put_z = ( .NOT. xy_ ) .AND. ( PRC_myrank == 0 ) ! only master process output the vertical coordinates
       put_x = ( rankidx(2) == 0 ) ! only south most row processes write x coordinates
       put_y = ( rankidx(1) == 0 ) ! only west most column processes write y coordinates
    else
       XSB = ISB
       XEB = IEB
       YSB = JSB
       YEB = JEB

       put_z = ( .NOT. xy_ )
       put_x = .true.
       put_y = .true.

       start(:) = 1
    end if

    if ( put_z ) then
       if( IO_AGGREGATE ) start(1) = 1
       call FileWriteAxis( fid, 'z',   GRID_CZ(KS:KE),    start )
       call FileWriteAxis( fid, 'zh',  GRID_FZ(KS:KE),    start )
       call FileWriteAxis( fid, 'lz',  GRID_LCZ(LKS:LKE), start )
       call FileWriteAxis( fid, 'lzh', GRID_LFZ(LKS:LKE), start )
       call FileWriteAxis( fid, 'uz',  GRID_UCZ(UKS:UKE), start )
       call FileWriteAxis( fid, 'uzh', GRID_UFZ(UKS:UKE), start )
    end if
    if ( put_x ) then
       if( IO_AGGREGATE ) start(1) = ISGA
       call FileWriteAxis( fid, 'x',   GRID_CX(XSB:XEB),  start )
       call FileWriteAxis( fid, 'xh',  GRID_FX(XSB:XEB),  start )
    end if
    if ( put_y ) then
       if( IO_AGGREGATE ) start(1) = JSGA
       call FileWriteAxis( fid, 'y',   GRID_CY(YSB:YEB),  start )
       call FileWriteAxis( fid, 'yh',  GRID_FY(YSB:YEB),  start )
    end if

    if ( put_z ) then
       if( IO_AGGREGATE ) start(1) = 1
       call FileWriteAxis( fid, 'CZ',   GRID_CZ,   start )
       call FileWriteAxis( fid, 'FZ',   GRID_FZ,   start )
       call FileWriteAxis( fid, 'CDZ',  GRID_CDZ,  start )
       call FileWriteAxis( fid, 'FDZ',  GRID_FDZ,  start )

       call FileWriteAxis( fid, 'LCZ',  GRID_LCZ,  start )
       call FileWriteAxis( fid, 'LFZ',  GRID_LFZ,  start )
       call FileWriteAxis( fid, 'LCDZ', GRID_LCZ,  start )

       call FileWriteAxis( fid, 'UCZ',  GRID_UCZ,  start )
       call FileWriteAxis( fid, 'UFZ',  GRID_UFZ,  start )
       call FileWriteAxis( fid, 'UCDZ', GRID_UCZ,  start )

       call FileWriteAxis( fid, 'CBFZ', GRID_CBFZ, start )
       call FileWriteAxis( fid, 'FBFZ', GRID_FBFZ, start )
    end if

    if ( IO_AGGREGATE ) then
       if ( PRC_myrank == 0 ) then
          start(1) = 1
          call FileWriteAxis( fid, 'CX',    GRID_CXG,   start )
          call FileWriteAxis( fid, 'CY',    GRID_CYG,   start )
          call FileWriteAxis( fid, 'FX',    GRID_FXG,   start )
          call FileWriteAxis( fid, 'FY',    GRID_FYG,   start )

          call FileWriteAxis( fid, 'CDX',   GRID_CDXG,  start )
          call FileWriteAxis( fid, 'CDY',   GRID_CDYG,  start )
          call FileWriteAxis( fid, 'FDX',   GRID_FDXG,  start )
          call FileWriteAxis( fid, 'FDY',   GRID_FDYG,  start )

          call FileWriteAxis( fid, 'CBFX',  GRID_CBFXG, start )
          call FileWriteAxis( fid, 'CBFY',  GRID_CBFYG, start )
          call FileWriteAxis( fid, 'FBFX',  GRID_FBFXG, start )
          call FileWriteAxis( fid, 'FBFY',  GRID_FBFYG, start )

          call FileWriteAxis( fid, 'CXG',   GRID_CXG,   start )
          call FileWriteAxis( fid, 'CYG',   GRID_CYG,   start )
          call FileWriteAxis( fid, 'FXG',   GRID_FXG,   start )
          call FileWriteAxis( fid, 'FYG',   GRID_FYG,   start )

          call FileWriteAxis( fid, 'CBFXG', GRID_CBFXG, start )
          call FileWriteAxis( fid, 'CBFYG', GRID_CBFYG, start )
          call FileWriteAxis( fid, 'FBFXG', GRID_FBFXG, start )
          call FileWriteAxis( fid, 'FBFYG', GRID_FBFYG, start )
       end if
    else
       call FileWriteAxis( fid, 'CX',    GRID_CX    )
       call FileWriteAxis( fid, 'CY',    GRID_CY    )
       call FileWriteAxis( fid, 'FX',    GRID_FX    )
       call FileWriteAxis( fid, 'FY',    GRID_FY    )

       call FileWriteAxis( fid, 'CDX',   GRID_CDX   )
       call FileWriteAxis( fid, 'CDY',   GRID_CDY   )
       call FileWriteAxis( fid, 'FDX',   GRID_FDX   )
       call FileWriteAxis( fid, 'FDY',   GRID_FDY   )

       call FileWriteAxis( fid, 'CBFX',  GRID_CBFX  )
       call FileWriteAxis( fid, 'CBFY',  GRID_CBFY  )
       call FileWriteAxis( fid, 'FBFX',  GRID_FBFX  )
       call FileWriteAxis( fid, 'FBFY',  GRID_FBFY  )

       call FileWriteAxis( fid, 'CXG',   GRID_CXG   )
       call FileWriteAxis( fid, 'CYG',   GRID_CYG   )
       call FileWriteAxis( fid, 'FXG',   GRID_FXG   )
       call FileWriteAxis( fid, 'FYG',   GRID_FYG   )

       call FileWriteAxis( fid, 'CBFXG', GRID_CBFXG )
       call FileWriteAxis( fid, 'CBFYG', GRID_CBFYG )
       call FileWriteAxis( fid, 'FBFXG', GRID_FBFXG )
       call FileWriteAxis( fid, 'FBFYG', GRID_FBFYG )
    end if

    ! associate coordinates

    if ( IO_AGGREGATE ) then
       start(1) = ISGA
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

    if ( .NOT. xy_ ) then
       if ( IO_AGGREGATE ) then
          start(1) = 1
          start(2) = ISGA
          start(3) = JSGA
       end if
       call FileWriteAssociatedCoordinates( fid, 'height',     AXIS_HZXY(:,:,:), start )
       call FileWriteAssociatedCoordinates( fid, 'height_wxy', AXIS_HWXY(:,:,:), start )
    end if

    return
  end subroutine FILEIO_write_axes

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

    integer,          intent(in)  :: fid      !< file ID
    integer,          intent(out) :: vid      !< variable ID
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: axistype !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    real(RP),         intent(in), optional :: timeintv !< time interval [sec]
    integer,          intent(in), optional :: nsteps   !< number of time steps

    integer          :: dtype, ndims, elm_size
    character(len=2) :: dims(3)
    real(DP)         :: time_interval
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if    ( datatype == 'REAL8' ) then
       dtype = File_REAL8
       elm_size = 4
    elseif( datatype == 'REAL4' ) then
       dtype = File_REAL4
       elm_size = 8
    else
       if    ( RP == 8 ) then
          dtype = File_REAL8
          elm_size = 8
       elseif( RP == 4 ) then
          dtype = File_REAL4
          elm_size = 4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_MPIstop
       endif
    endif

    if    ( axistype == 'Z' ) then        ! 1D variable
       ndims   = 1
       dims(1) = 'z'
       write_buf_amount = write_buf_amount + KA * elm_size
    elseif( axistype == 'X' ) then
       ndims   = 1
       dims(1) = 'x'
       write_buf_amount = write_buf_amount + IA * elm_size
    elseif( axistype == 'Y' ) then
       ndims   = 1
       dims(1) = 'y'
       write_buf_amount = write_buf_amount + JA * elm_size
    elseif( axistype == 'XY' ) then   ! 2D variable
       ndims   = 2
       dims(1) = 'x'
       dims(2) = 'y'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif( axistype == 'UY' ) then
       ndims   = 2
       dims(1) = 'xh'
       dims(2) = 'y'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif( axistype == 'XV' ) then
       ndims   = 2
       dims(1) = 'x'
       dims(2) = 'yh'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif( axistype == 'UV' ) then
       ndims   = 2
       dims(1) = 'xh'
       dims(2) = 'yh'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
    elseif( axistype == 'ZX' ) then
       ndims   = 2
       dims(1) = 'z'
       dims(2) = 'x'
       write_buf_amount = write_buf_amount + KA * IA * elm_size
    elseif( axistype == 'ZXY' ) then  ! 3D variable
       ndims   = 3
       dims    = (/'z','x','y'/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype == 'ZHXY' ) then
       ndims   = 3
       dims    = (/'zh','x ','y '/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype == 'ZXHY' ) then
       ndims   = 3
       dims    = (/'z ','xh','y '/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype == 'ZXYH' ) then
       ndims   = 3
       dims    = (/'z ','x ','yh'/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
    elseif( axistype == 'Land' ) then
       ndims   = 3
       dims    = (/'lz','x ','y '/)
       write_buf_amount = write_buf_amount + LKMAX * IA * JA * elm_size
    elseif( axistype == 'Urban' ) then
       ndims   = 3
       dims    = (/'uz','x ','y '/)
       write_buf_amount = write_buf_amount + UKMAX * IA * JA * elm_size
    elseif( axistype == 'XYT' ) then  ! 3D variable with time dimension
       ndims   = 2
       dims(1) = 'x'
       dims(2) = 'y'
       if ( present(nsteps) ) then
          write_buf_amount = write_buf_amount + IA * JA * elm_size * nsteps
       else
          write_buf_amount = write_buf_amount + IA * JA * elm_size
       end if
    elseif( axistype == 'ZXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'z','x','y'/)
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
                               tint=time_interval                                 ) ! [IN]
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
       FileWrite
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    implicit none

    integer,          intent(in) :: fid      !< file ID
    integer,          intent(in) :: vid      !< netCDF variable ID
    real(RP),         intent(in) :: var(:)   !< value of the variable
    character(len=*), intent(in) :: varname  !< name of the variable
    character(len=*), intent(in) :: axistype !< axis type (Z/X/Y)

    integer :: dim1_S, dim1_E
    integer :: rankidx(2)
    integer :: start(1)         ! used only when IO_AGGREGATE is .true.
    logical :: exec = .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if    ( axistype == 'Z' ) then
       dim1_S   = KS
       dim1_E   = KE
       start(1) = 1
       if( IO_AGGREGATE .AND. PRC_myrank > 0 ) exec = .false.  ! only rank 0 writes
    elseif( axistype == 'X' ) then
       dim1_S   = ISB
       dim1_E   = IEB
       start(1) = ISGA
       if( IO_AGGREGATE .AND. rankidx(2) > 0 ) exec = .false.  ! only south most row processes write
    elseif( axistype == 'Y' ) then
       dim1_S   = JSB
       dim1_E   = JEB
       start(1) = JSGA
       if( IO_AGGREGATE .AND. rankidx(1) > 0 ) exec = .false.  ! only west most column processes write
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    if( exec ) call FileWrite( fid, vid, var(dim1_S:dim1_E), &
                               NOWSEC, NOWSEC, start         ) ! [IN]

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
       FileWrite
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
    logical,          intent(in), optional :: nohalo   !< switch whether include halo data or not (default=false)

    real(RP)              :: varhalo( size(var(:,1)), size(var(1,:)) )

    integer               :: dim1_S, dim1_E
    integer               :: dim2_S, dim2_E

    integer :: i, j
    logical :: nohalo_
    integer :: rankidx(2)
    integer :: start(2)         ! used only when IO_AGGREGATE is .true.
    logical :: exec = .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = ISGA
    start(2) = JSGA

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    if (      axistype == 'XY' &
         .OR. axistype == 'UY' &
         .OR. axistype == 'XV' &
         .OR. axistype == 'UV' ) then
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       if ( IO_AGGREGATE ) then
          if( rankidx(1) == 0             ) dim1_S = 1
          if( rankidx(1) == PRC_NUM_X - 1 ) dim1_E = IA
          if( rankidx(2) == 0             ) dim2_S = 1
          if( rankidx(2) == PRC_NUM_Y - 1 ) dim2_E = JA
       end if
    elseif( axistype == 'ZX' ) then
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       start(2) = start(1)
       start(1) = 1
       if ( IO_AGGREGATE .AND. rankidx(2) > 0 ) then
          exec = .false.  ! only south most row processes write
          if( rankidx(1) == 0             ) dim2_S = 1
          if( rankidx(1) == PRC_NUM_X - 1 ) dim2_E = IA
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

          call FileWrite( fid, vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), &
                          NOWSEC, NOWSEC, start                           ) ! [IN]
       else
          call FileWrite( fid, vid, var(dim1_S:dim1_E,dim2_S:dim2_E), &
                          NOWSEC, NOWSEC, start                       ) ! [IN]
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
       FileWrite
    use scale_process, only: &
       PRC_myrank,  &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X,  &
       PRC_NUM_Y
    use scale_time, only: &
       NOWSEC  => TIME_NOWDAYSEC
    implicit none

    integer,          intent(in) :: fid        !< file ID
    integer,          intent(in) :: vid        !< netCDF variable ID
    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: varname    !< name        of the variable
    character(len=*), intent(in) :: axistype   !< axis type (Z/X/Y)

    logical,          intent(in), optional :: nohalo !< include halo data?

    real(RP) :: varhalo( size(var(:,1,1)), size(var(1,:,1)), size(var(1,1,:)) )

    integer  :: dim1_S, dim1_E, dim1_max
    integer  :: dim2_S, dim2_E
    integer  :: dim3_S, dim3_E

    integer  :: i, j, k
    logical  :: nohalo_
    integer  :: rankidx(2)
    integer  :: start(3) ! used only when IO_AGGREGATE is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA

    dim2_S   = ISB
    dim2_E   = IEB
    dim3_S   = JSB
    dim3_E   = JEB
    if ( IO_AGGREGATE ) then
       if( rankidx(1) == 0             ) dim2_S = 1
       if( rankidx(1) == PRC_NUM_X - 1 ) dim2_E = IA
       if( rankidx(2) == 0             ) dim3_S = 1
       if( rankidx(2) == PRC_NUM_Y - 1 ) dim3_E = JA
    end if

    if (      axistype == 'ZXY'  &
         .OR. axistype == 'ZHXY' &
         .OR. axistype == 'ZXHY' &
         .OR. axistype == 'ZXYH' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif( axistype == 'Land' ) then
       dim1_max = LKMAX
       dim1_S   = LKS
       dim1_E   = LKE
    elseif( axistype == 'Urban' ) then
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

       call FileWrite( fid, vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                       NOWSEC, NOWSEC, start                                         ) ! [IN]
    else
       call FileWrite( fid, vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                       NOWSEC, NOWSEC, start                                     ) ! [IN]
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
       FileWrite
    use scale_process, only: &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    integer,          intent(in) :: fid        !< file ID
    integer,          intent(in) :: vid        !< netCDF variable ID
    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: varname    !< name of the variable
    character(len=*), intent(in) :: axistype   !< axis type (X/Y/Time)
    real(RP),         intent(in) :: timeintv   !< time interval [sec]

    integer,          intent(in), optional :: timetarg !< target timestep (optional)
    logical,          intent(in), optional :: nohalo   !< include halo data?

    real(RP) :: varhalo( size(var(:,1,1)), size(var(1,:,1)) )

    integer  :: dim1_S, dim1_E
    integer  :: dim2_S, dim2_E

    real(DP) :: time_interval, nowtime

    integer  :: step
    integer  :: i, j, n
    logical  :: nohalo_
    integer  :: rankidx(2)
    integer  :: start(3) ! used only when IO_AGGREGATE is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    time_interval = timeintv
    step = size(var(ISB,JSB,:))

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( axistype == 'XYT' ) then
       dim1_S   = ISB
       dim1_E   = IEB
       dim2_S   = JSB
       dim2_E   = JEB
       if ( IO_AGGREGATE ) then
          if( rankidx(1) == 0             ) dim1_S = 1
          if( rankidx(1) == PRC_NUM_X - 1 ) dim1_E = IA
          if( rankidx(2) == 0             ) dim2_S = 1
          if( rankidx(2) == PRC_NUM_Y - 1 ) dim2_E = JA
       end if
    else
       write(*,*) 'xxx unsupported axis type. Check!', trim(axistype), ' item:',trim(varname)
       call PRC_MPIstop
    endif

    start(1) = ISGA
    start(2) = JSGA
    ! start(3) time dimension will be set in file_write_data()

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

          call FileWrite( fid, vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), &
                          nowtime, nowtime, start                         ) ! [IN]
       else
          call FileWrite( fid, vid, var(dim1_S:dim1_E,dim2_S:dim2_E,timetarg), &
                          nowtime, nowtime, start                              ) ! [IN]
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

             call FileWrite( fid, vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), &
                             nowtime, nowtime, start                         ) ! [IN]
          else
             call FileWrite( fid, vid, var(dim1_S:dim1_E,dim2_S:dim2_E,n), &
                             nowtime, nowtime, start                       ) ! [IN]
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
       FileWrite
    use scale_process, only: &
       PRC_myrank,     &
       PRC_MPIstop
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    integer,          intent(in) :: fid          !< file ID
    integer,          intent(in) :: vid          !< netCDF variable ID
    real(RP),         intent(in) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in) :: varname      !< name        of the variable
    character(len=*), intent(in) :: axistype     !< axis type (Z/X/Y/Time)
    real(RP),         intent(in) :: timeintv     !< time interval [sec]

    integer,          intent(in), optional :: timetarg !< target timestep (optional)
    logical,          intent(in), optional :: nohalo   !< include halo data?

    real(RP) :: varhalo( size(var(:,1,1,1)), size(var(1,:,1,1)), size(var(1,1,:,1)) )

    integer  :: dim1_S, dim1_E, dim1_max
    integer  :: dim2_S, dim2_E
    integer  :: dim3_S, dim3_E

    real(DP) :: time_interval, nowtime

    integer  :: step
    integer  :: i, j, k, n
    logical  :: nohalo_
    integer  :: rankidx(2)
    integer  :: start(4) ! used only when IO_AGGREGATE is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    nohalo_ = .false.
    if( present(nohalo) ) nohalo_ = nohalo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA
    ! start(4) time dimension will be set in file_write_data()

    time_interval = timeintv
    step = size(var,4)

    if ( axistype == 'ZXYT' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
       if ( IO_AGGREGATE ) then
          if( rankidx(1) == 0             ) dim2_S = 1
          if( rankidx(1) == PRC_NUM_X - 1 ) dim2_E = IA
          if( rankidx(2) == 0             ) dim3_S = 1
          if( rankidx(2) == PRC_NUM_Y - 1 ) dim3_E = JA
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

          call FileWrite( fid, vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                          nowtime, nowtime, start                                       ) ! [IN]
       else
          call FileWrite( fid, vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,timetarg), &
                          nowtime, nowtime, start                                            ) ! [IN]
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

             call FileWrite( fid, vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                             nowtime, nowtime, start                                       ) ! [IN]
          else
             call FileWrite( fid, vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,n), &
                             nowtime, nowtime, start                                     ) ! [IN]
          end if
          nowtime = nowtime + time_interval
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILEIO_write_var_4D

  !-----------------------------------------------------------------------------
  subroutine closeall
    implicit none

    integer :: fid
    !---------------------------------------------------------------------------

    do fid = 0, File_nfile_max-1
       call FILEIO_close( fid )
    end do

    return
  end subroutine closeall

  !-----------------------------------------------------------------------------
  subroutine getCFtunits(tunits, date)
    implicit none

    character(len=34), intent(out) :: tunits
    integer,           intent(in)  :: date(6)
    !---------------------------------------------------------------------------

    write(tunits,'(a,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') 'seconds since ', date

    return
  end subroutine getCFtunits

  !-----------------------------------------------------------------------------
  subroutine check_1d( &
       expected, buffer, &
       name              )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),         intent(in) :: expected(:)
    real(RP),         intent(in) :: buffer(:)
    character(len=*), intent(in) :: name

    integer :: nmax, n

    intrinsic :: size
    !---------------------------------------------------------------------------

    nmax = size(expected)
    if ( size(buffer) /= nmax ) then
       write(*,*) 'xxx size of coordinate ('//trim(name)//') is different:', nmax, size(buffer)
       call PRC_MPIstop
    end if

    do n=1, nmax
       if ( expected(n) /= buffer(n) ) then
          write(*,*) 'xxx value of coordinate ('//trim(name)//') at ', n, ' is different:', &
                     expected(n), buffer(n)
          call PRC_MPIstop
       end if
    end do

    return
  end subroutine check_1d

  !-----------------------------------------------------------------------------
  subroutine check_2d( &
       expected, buffer, &
       name              )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),         intent(in) :: expected(:,:)
    real(RP),         intent(in) :: buffer(:,:)
    character(len=*), intent(in) :: name

    integer :: imax, jmax
    integer :: i, j

    intrinsic :: size
    !---------------------------------------------------------------------------

    imax = size(expected,1)
    jmax = size(expected,2)
    if ( size(buffer,1) /= imax ) then
       write(*,*) 'xxx the first size of coordinate ('//trim(name)//') is different:', imax, size(buffer,1)
       call PRC_MPIstop
    end if
    if ( size(buffer,2) /= jmax ) then
       write(*,*) 'xxx the second size of coordinate ('//trim(name)//') is different:', jmax, size(buffer,2)
       call PRC_MPIstop
    end if

    do j=1, jmax
    do i=1, imax
       if ( expected(i,j) /= buffer(i,j) ) then
          write(*,*) 'xxx value of coordinate ('//trim(name)//') at (', i, ',', j, ') is different:', &
                     expected(i,j), buffer(i,j)
          call PRC_MPIstop
       end if
    end do
    end do

    return
  end subroutine check_2d

  !-----------------------------------------------------------------------------
  subroutine check_3d( &
       expected, buffer, &
       name              )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),         intent(in) :: expected(:,:,:)
    real(RP),         intent(in) :: buffer(:,:,:)
    character(len=*), intent(in) :: name

    integer :: imax, jmax, kmax
    integer :: i, j, k

    intrinsic :: size
    !---------------------------------------------------------------------------

    kmax = size(expected,1)
    imax = size(expected,2)
    jmax = size(expected,3)
    if ( size(buffer,1) /= kmax ) then
       write(*,*) 'xxx the first size of coordinate ('//trim(name)//') is different:', kmax, size(buffer,1)
       call PRC_MPIstop
    end if
    if ( size(buffer,2) /= imax ) then
       write(*,*) 'xxx the second size of coordinate ('//trim(name)//') is different:', imax, size(buffer,2)
       call PRC_MPIstop
    end if
    if ( size(buffer,3) /= jmax ) then
       write(*,*) 'xxx the second size of coordinate ('//trim(name)//') is different:', jmax, size(buffer,3)
       call PRC_MPIstop
    end if

    do j=1, jmax
    do i=1, imax
    do k=1, kmax
       if ( expected(k,i,j) /= buffer(k,i,j) ) then
          write(*,*) 'xxx value of coordinate ('//trim(name)//') at ', k, ',', i, ',', j, ' is different:', &
                     expected(k,i,j), buffer(k,i,j)
          call PRC_MPIstop
       end if
    end do
    end do
    end do

    return
  end subroutine check_3d

end module scale_fileio
