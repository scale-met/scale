!-------------------------------------------------------------------------------
!> module file / cartesianC
!!
!! @par Description
!!          file I/O handling for the cartesianC grid
!!
!! @author Team SCALE
!!
!<
module scale_file_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  use scale_ocean_grid_index
  use scale_urban_grid_index
  use scale_file_h, only: &
     FILE_FILE_MAX
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FILE_CARTESC_setup
  public :: FILE_CARTESC_cleanup
  public :: FILE_CARTESC_set_coordinates
  public :: FILE_CARTESC_check_coordinates

  public :: FILE_CARTESC_create
  public :: FILE_CARTESC_open
  public :: FILE_CARTESC_def_var
  public :: FILE_CARTESC_enddef
  public :: FILE_CARTESC_write_var
  public :: FILE_CARTESC_read
  public :: FILE_CARTESC_write
  public :: FILE_CARTESC_flush
  public :: FILE_CARTESC_close

  interface FILE_CARTESC_check_coordinates
     module procedure FILE_CARTESC_check_coordinates_name
     module procedure FILE_CARTESC_check_coordinates_id
  end interface FILE_CARTESC_check_coordinates

  interface FILE_CARTESC_read
     module procedure FILE_CARTESC_read_1D
     module procedure FILE_CARTESC_read_2D
     module procedure FILE_CARTESC_read_3D
     module procedure FILE_CARTESC_read_4D

     module procedure FILE_CARTESC_read_var_1D
     module procedure FILE_CARTESC_read_var_2D
     module procedure FILE_CARTESC_read_var_3D
     module procedure FILE_CARTESC_read_var_4D
  end interface FILE_CARTESC_read

  interface FILE_CARTESC_write
     module procedure FILE_CARTESC_write_1D
     module procedure FILE_CARTESC_write_2D
     module procedure FILE_CARTESC_write_3D
     module procedure FILE_CARTESC_write_3D_t
     module procedure FILE_CARTESC_write_4D
  end interface FILE_CARTESC_write

  interface FILE_CARTESC_write_var
     module procedure FILE_CARTESC_write_var_1D
     module procedure FILE_CARTESC_write_var_2D
     module procedure FILE_CARTESC_write_var_3D
     module procedure FILE_CARTESC_write_var_3D_t
     module procedure FILE_CARTESC_write_var_4D
  end interface FILE_CARTESC_write_var

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type, public :: axisattinfo
    integer          :: size_global (1)
    integer          :: start_global(1)
    integer          :: halo_global (2)
    integer          :: halo_local  (2)
    character(len=5) :: periodic
  end type axisattinfo

  type, public :: mappinginfo
    character(len=H_SHORT) :: mapping_name
    real(DP)               :: false_easting                        (1)
    real(DP)               :: false_northing                       (1)
    real(DP)               :: longitude_of_central_meridian        (1)
    real(DP)               :: longitude_of_projection_origin       (1)
    real(DP)               :: latitude_of_projection_origin        (1)
    real(DP)               :: straight_vertical_longitude_from_pole(1)
    real(DP)               :: standard_parallel                    (2)
  end type mappinginfo

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: closeall
  private :: check_1d
  private :: check_2d
  private :: check_3d

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private              :: FILE_CARTESC_datacheck_criteria

  real(RP), private, allocatable :: AXIS_LON  (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONUY(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONXV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONUV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LAT  (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATUY(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATXV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATUV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_HGT   (:,:,:)
  real(RP), private, allocatable :: AXIS_HGTWXY(:,:,:)

  logical,  private              :: File_axes_written(0:FILE_FILE_MAX-1) ! whether axes have been written
  !                                                                       ! fid starts from zero so index should start from zero
  logical,  private              :: File_closed      (0:FILE_FILE_MAX-1) ! whether file has been closed
  logical,  private              :: File_haszcoord   (0:FILE_FILE_MAX-1) ! z-coordinates exist?
  integer,  private              :: write_buf_amount = 0                  ! sum of write buffer amounts

  ! global star and count
  integer,  private              :: startXY   (3), countXY   (3)
  integer,  private              :: startZX   (2), countZX   (2)
  integer,  private              :: startZXY  (4), countZXY  (4)
  integer,  private              :: startZHXY (4), countZHXY (4)
  integer,  private              :: startOCEAN(3), countOCEAN(3)
  integer,  private              :: startLAND (3), countLAND (3)
  integer,  private              :: startURBAN(3), countURBAN(3)
  ! local start and end
  integer,  private              :: ISB2, IEB2, JSB2, JEB2 !> for FILE_AGGREGATE

  ! MPI element datatype for restart variables
  integer,  private              :: etype

  ! MPI derived datatypes
  integer,  private              :: centerTypeXY
  integer,  private              :: centerTypeZX
  integer,  private              :: centerTypeZXY
  integer,  private              :: centerTypeZHXY
  integer,  private              :: centerTypeOCEAN
  integer,  private              :: centerTypeLAND
  integer,  private              :: centerTypeURBAN

  logical,  private              :: set_coordinates = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILE_CARTESC_setup
    use scale_process, only: &
       PRC_myrank, &
       PRC_abort
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_file, only: &
       FILE_setup
    implicit none

    NAMELIST / PARAM_FILE_CARTESC / &
       FILE_CARTESC_datacheck_criteria

    integer :: IM, JM
    integer :: ierr
    !---------------------------------------------------------------------------

    call FILE_setup( PRC_myrank )


    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[FIELIO] / Categ[ATMOS-RM IO] / Origin[SCALElib]'

    FILE_CARTESC_datacheck_criteria = EPS * 10.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_FILE_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_FILE_CARTESC. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_FILE_CARTESC)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** NetCDF header information ***'
    if( IO_L ) write(IO_FID_LOG,*) '*** Data source : ', trim(H_SOURCE)
    if( IO_L ) write(IO_FID_LOG,*) '*** Institute   : ', trim(H_INSTITUTE)
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Data consistency criteria : ', &
                                   '(file-internal)/internal = ', FILE_CARTESC_datacheck_criteria

    ! construct indices independent from PRC_PERIODIC_X/Y
    ISB2 = IS
    if( PRC_2Drank(PRC_myrank,1) == 0           ) ISB2 = 1
    IEB2 = IE
    if( PRC_2Drank(PRC_myrank,1) == PRC_NUM_X-1 ) IEB2 = IA

    JSB2 = JS
    if( PRC_2Drank(PRC_myrank,2) == 0           ) JSB2 = 1
    JEB2 = JE
    if( PRC_2Drank(PRC_myrank,2) == PRC_NUM_Y-1 ) JEB2 = JA

    IM = JEB2 - JSB2 + 1
    JM = JEB2 - JSB2 + 1
    allocate( AXIS_LON  (IM,JM) )
    allocate( AXIS_LONUY(IM,JM) )
    allocate( AXIS_LONXV(IM,JM) )
    allocate( AXIS_LONUV(IM,JM) )
    allocate( AXIS_LAT  (IM,JM) )
    allocate( AXIS_LATUY(IM,JM) )
    allocate( AXIS_LATXV(IM,JM) )
    allocate( AXIS_LATUV(IM,JM) )

    allocate( AXIS_HGT   (KMAX  ,IM,JM) )
    allocate( AXIS_HGTWXY(KMAX+1,IM,JM) )

    call Construct_Derived_Datatype

    File_closed(:) = .true.

    return
  end subroutine FILE_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> deallocate buffers
  subroutine FILE_CARTESC_cleanup
    implicit none
    !---------------------------------------------------------------------------

    deallocate( AXIS_LON   )
    deallocate( AXIS_LONUY )
    deallocate( AXIS_LONXV )
    deallocate( AXIS_LONUV )
    deallocate( AXIS_LAT   )
    deallocate( AXIS_LATUY )
    deallocate( AXIS_LATXV )
    deallocate( AXIS_LATUV )
    deallocate( AXIS_HGT    )
    deallocate( AXIS_HGTWXY )

    call Free_Derived_Datatype

    call closeall

    return
  end subroutine FILE_CARTESC_cleanup

  !-----------------------------------------------------------------------------
  !> set latlon and z
  subroutine FILE_CARTESC_set_coordinates( &
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

    real(RP), intent(in) :: LON  (  IA,  JA)
    real(RP), intent(in) :: LONX (0:IA,  JA)
    real(RP), intent(in) :: LONY (  IA,0:JA)
    real(RP), intent(in) :: LONXY(0:IA,0:JA)
    real(RP), intent(in) :: LAT  (  IA,  JA)
    real(RP), intent(in) :: LATX (0:IA,  JA)
    real(RP), intent(in) :: LATY (  IA,0:JA)
    real(RP), intent(in) :: LATXY(0:IA,0:JA)
    real(RP), intent(in) :: CZ   (  KA,IA,JA)
    real(RP), intent(in) :: FZ   (0:KA,IA,JA)
    !---------------------------------------------------------------------------

    AXIS_LON  (:,:)   = LON  (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LONUY(:,:)   = LONX (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LONXV(:,:)   = LONY (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LONUV(:,:)   = LONXY(ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LAT  (:,:)   = LAT  (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LATUY(:,:)   = LATX (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LATXV(:,:)   = LATY (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LATUV(:,:)   = LATXY(ISB2:IEB2,JSB2:JEB2) / D2R

    AXIS_HGT   (:,:,:) = CZ(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_HGTWXY(:,:,:) = FZ(KS-1:KE,ISB2:IEB2,JSB2:JEB2)

    set_coordinates = .true.

    return
  end subroutine FILE_CARTESC_set_coordinates

  !-----------------------------------------------------------------------------
  !> check coordinates in the file
  subroutine FILE_CARTESC_check_coordinates_name( &
       basename, &
       atmos,    &
       ocean,    &
       land,     &
       urban,    &
       transpose )
    implicit none

    character(len=*), intent(in) :: basename        !< basename of the file
    logical,          intent(in), optional :: atmos !< check atmospheric coordinates
    logical,          intent(in), optional :: ocean !< check ocean coordinates
    logical,          intent(in), optional :: land  !< check land coordinates
    logical,          intent(in), optional :: urban !< check urban coordinates
    logical,          intent(in), optional :: transpose

    logical :: atmos_
    logical :: ocean_
    logical :: land_
    logical :: urban_
    logical :: transpose_

    integer :: fid
    !---------------------------------------------------------------------------

    atmos_ = .false.
    ocean_ = .false.
    land_  = .false.
    urban_ = .false.
    transpose_ = .false.

    if( present(atmos) ) atmos_ = atmos
    if( present(ocean) ) ocean_ = ocean
    if( present(land ) ) land_  = land
    if( present(urban) ) urban_ = urban
    if( present(transpose) ) transpose_ = transpose

    call FILE_CARTESC_open( fid,     & ! [OUT]
                            basename ) ! [IN]

    call FILE_CARTESC_check_coordinates_id( fid,                           & ! [IN]
                                            atmos_, ocean_, land_, urban_, & ! [IN]
                                            transpose_                     ) ! [IN]

    return
  end subroutine FILE_CARTESC_check_coordinates_name

  !-----------------------------------------------------------------------------
  !> check coordinates in the file
  subroutine FILE_CARTESC_check_coordinates_id( &
       fid,   &
       atmos, &
       ocean, &
       land,  &
       urban, &
       transpose )
    use scale_grid, only: &
       GRID_CZ, &
       GRID_CX, &
       GRID_CY
    use scale_ocean_grid, only: &
       GRID_OCZ
    use scale_land_grid, only: &
       GRID_LCZ
    use scale_urban_grid, only: &
       GRID_UCZ
    implicit none

    integer, intent(in) :: fid
    logical, intent(in), optional :: atmos !< check atmospheric coordinates
    logical, intent(in), optional :: ocean !< check ocean coordinates
    logical, intent(in), optional :: land  !< check land coordinates
    logical, intent(in), optional :: urban !< check urban coordinates
    logical, intent(in), optional :: transpose

    logical :: atmos_
    logical :: ocean_
    logical :: land_
    logical :: urban_
    logical :: transpose_

    real(RP) :: buffer_z  (KA)
    real(RP) :: buffer_x  (IA)
    real(RP) :: buffer_y  (JA)
    real(RP) :: buffer_xy (IA,JA)
    real(RP) :: buffer_zxy(KA,IA,JA)
    real(RP) :: buffer_o  (OKMAX)
    real(RP) :: buffer_l  (LKMAX)
    real(RP) :: buffer_u  (UKMAX)

    integer :: XSB, XEB, YSB, YEB
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Check consistency of axis ***'

    atmos_ = .false.
    ocean_ = .false.
    land_  = .false.
    urban_ = .false.
    transpose_ = .false.

    if( present(atmos) ) atmos_ = atmos
    if( present(ocean) ) ocean_ = ocean
    if( present(land ) ) land_  = land
    if( present(urban) ) urban_ = urban
    if( present(transpose) ) transpose_ = transpose


    XSB = ISB - ISB2 + 1
    XEB = IEB - ISB + XSB
    YSB = JSB - JSB2 + 1
    YEB = JEB - JSB + YSB

    call FILE_CARTESC_read_var_1D( buffer_x, fid, 'x',  'X', 1 )
    call FILE_CARTESC_read_var_1D( buffer_y, fid, 'y',  'Y', 1 )
    call FILE_CARTESC_flush( fid )
    call check_1d( GRID_CX(ISB:IEB), buffer_x(ISB:IEB), 'x' )
    call check_1d( GRID_CY(JSB:JEB), buffer_y(JSB:JEB), 'y' )

    if ( set_coordinates ) then
       call FILE_CARTESC_read_var_2D( buffer_xy, fid, 'lon', 'XY', 1 )
       call FILE_CARTESC_flush( fid )
       call check_2d( AXIS_LON(XSB:XEB,YSB:YEB), buffer_xy(ISB:IEB,JSB:JEB), 'lon' )

       call FILE_CARTESC_read_var_2D( buffer_xy, fid, 'lat', 'XY', 1 )
       call FILE_CARTESC_flush( fid )
       call check_2d( AXIS_LAT(XSB:XEB,YSB:YEB), buffer_xy(ISB:IEB,JSB:JEB), 'lat' )
    endif

    if ( atmos_ ) then
       call FILE_CARTESC_read_var_1D( buffer_z,   fid, 'z',      'Z',   1 )
       if ( .not. transpose_ ) then
          call FILE_CARTESC_read_var_3D( buffer_zxy, fid, 'height', 'ZXY', 1 )
       endif
       call FILE_CARTESC_flush( fid )
       call check_1d( GRID_CZ(KS:KE), buffer_z(KS:KE), 'z' )
       if ( .not. transpose_ ) then
          call check_3d( AXIS_HGT(:,XSB:XEB,YSB:YEB), buffer_zxy(KS:KE,ISB:IEB,JSB:JEB), 'height', transpose_ )
       endif
    endif

    if ( ocean_ ) then
       call FILE_CARTESC_read_var_1D( buffer_o, fid, 'oz', 'OZ', 1 )
       call FILE_CARTESC_flush( fid )
       call check_1d( GRID_OCZ(OKS:OKE), buffer_o(OKS:OKE), 'oz' )
    endif

    if ( land_ ) then
       call FILE_CARTESC_read_var_1D( buffer_l, fid, 'lz',  'LZ', 1 )
       call FILE_CARTESC_flush( fid )
       call check_1d( GRID_LCZ(LKS:LKE), buffer_l(LKS:LKE), 'lz' )
    endif

    if ( urban_ ) then
       call FILE_CARTESC_read_var_1D( buffer_u, fid, 'uz',  'UZ', 1 )
       call FILE_CARTESC_flush( fid )
       call check_1d( GRID_UCZ(UKS:UKE), buffer_u(UKS:UKE), 'uz' )
    endif

    return
  end subroutine FILE_CARTESC_check_coordinates_id

  !-----------------------------------------------------------------------------
  !> construct MPI derived datatypes for read buffers
  subroutine Construct_Derived_Datatype
    use mpi
    use scale_process, only: &
       PRC_abort
    use scale_rm_process, only: &
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
    centerTypeZHXY  = MPI_DATATYPE_NULL
    centerTypeOCEAN = MPI_DATATYPE_NULL
    centerTypeLAND  = MPI_DATATYPE_NULL
    centerTypeURBAN = MPI_DATATYPE_NULL

    etype           = MPI_FLOAT

    if( RP == 8 ) etype = MPI_DOUBLE_PRECISION

    ! for dim_type == 'XY'
    startXY(1)    = IS_inG - IHALO
    startXY(2)    = JS_inG - JHALO
    countXY(1)    = IA
    countXY(2)    = JA
    ! for dim_type == 'ZXY'
    startZXY(1)    = 1
    startZXY(2:3)  = startXY(1:2)
    countZXY(1)    = KMAX
    countZXY(2:3)  = countXY(1:2)
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

    ! for dim_type == 'ZHXY'
    startZHXY(1)   = 1
    startZHXY(2:3) = startXY(1:2)
    countZHXY(1)   = KMAX+1
    countZHXY(2:3) = countXY(1:2)
    ! construct MPI subarray data type
    sizes(1)      = KA
    sizes(2)      = IA
    sizes(3)      = JA
    subsizes(1)   = KMAX+1
    subsizes(2)   = IA
    subsizes(3)   = JA
    sub_off(1)    = KS - 2 ! MPI start index starts with 0
    sub_off(2)    = 0
    sub_off(3)    = 0
    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeZHXY, err)
    call MPI_Type_commit(centerTypeZHXY, err)

    ! for dim_type == 'OXY'
    startOCEAN(1)   = 1
    startOCEAN(2:3) = startXY(1:2)
    countOCEAN(1)   = OKMAX
    countOCEAN(2:3) = countXY(1:2)
    ! construct MPI subarray data type
    sizes(1)       = OKMAX
    subsizes(1)    = OKMAX
    sub_off(1)     = OKS - 1 ! MPI start index starts with 0
    call MPI_Type_create_subarray(3, sizes, subsizes, sub_off, order, etype, centerTypeOCEAN, err)
    call MPI_Type_commit(centerTypeOCEAN, err)

    ! for dim_type == 'LXY'
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

    ! for dim_type == 'UXY'
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

    ! for dim_type == 'ZX'
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

    if( centerTypeXY    /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeXY,    err)
    if( centerTypeZX    /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZX,    err)
    if( centerTypeZXY   /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZXY,   err)
    if( centerTypeZHXY  /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeZHXY,  err)
    if( centerTypeOCEAN /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeOCEAN, err)
    if( centerTypeLAND  /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeLAND,  err)
    if( centerTypeURBAN /= MPI_DATATYPE_NULL ) call MPI_Type_free(centerTypeURBAN, err)

    return
  end subroutine Free_Derived_Datatype

  !-----------------------------------------------------------------------------
  !> interface FILE_CARTESC_read
  !!    Read data from file
  !!    This routine is a wrapper of the lower primitive routines
  !<
  !-----------------------------------------------------------------------------
  subroutine FILE_CARTESC_read_1D( &
       var,      &
       basename, &
       varname,  &
       dim_type, &
       step,     &
       aggregate )
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number
    logical,          intent(in), optional :: aggregate

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( fid,      & ! [OUT]
                            basename, & ! [IN]
                            aggregate ) ! [IN]

    call FILE_CARTESC_read_var_1D( var(:),                      & ! [OUT]
                                   fid, varname, dim_type, step ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_1D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILE_CARTESC_read_2D( &
       var,      &
       basename, &
       varname,  &
       dim_type, &
       step,     &
       aggregate )
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value of the variable
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number
    logical,          intent(in), optional :: aggregate

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( fid,      & ! [OUT]
                            basename, & ! [IN]
                            aggregate ) ! [IN]

    call FILE_CARTESC_read_var_2D( var(:,:),                    & ! [OUT]
                                   fid, varname, dim_type, step ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILE_CARTESC_read_3D( &
       var,      &
       basename, &
       varname,  &
       dim_type, &
       step,     &
       aggregate )
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: dim_type   !< dimension type (Z/X/Y/T)
    integer,          intent(in)  :: step       !< step number
    logical,          intent(in), optional :: aggregate

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( fid,      & ! [OUT]
                            basename, & ! [IN]
                            aggregate )

    call FILE_CARTESC_read_var_3D( var(:,:,:),                  & ! [OUT]
                                   fid, varname, dim_type, step ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_3D

  !-----------------------------------------------------------------------------
  !> Read 4D data from file
  subroutine FILE_CARTESC_read_4D( &
       var,      &
       basename, &
       varname,  &
       dim_type, &
       step,     &
       aggregate )
    implicit none

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in)  :: basename     !< basename of the file
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: dim_type     !< dimension type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number
    logical,          intent(in), optional :: aggregate

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( fid,      & ! [OUT]
                            basename, & ! [IN]
                            aggregate )

    call FILE_CARTESC_read_var_4D( var(:,:,:,:),                & ! [OUT]
                                   fid, varname, dim_type, step ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_4D

  !-----------------------------------------------------------------------------
  !> Read 1D data from file
  subroutine FILE_CARTESC_read_var_1D( &
       var,      &
       fid,      &
       varname,  &
       dim_type, &
       step      )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Read
    use scale_process, only: &
       PRC_abort
    use scale_rm_process, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    use mpi
    implicit none

    real(RP),         intent(out) :: var(:)   !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer :: dim1_S, dim1_E
    integer :: start(1)   ! start offset of globale variable
    integer :: count(1)   ! request length to the global variable
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Read from file (1D), name : ', trim(varname)

    if ( FILE_get_aggregate(fid) ) then
       ! read data and halos into the local buffer
       if    ( dim_type == 'Z' ) then
          start(1) = 1
          count(1) = KMAX
          call FILE_Read( var(KS:KE), fid, varname, step,                    &
                          ntypes=KMAX, dtype=etype, start=start, count=count )
       elseif( dim_type == 'OZ' ) then
          start(1) = 1
          count(1) = OKMAX
          call FILE_Read( var, fid, varname, step,                            &
                          ntypes=OKMAX, dtype=etype, start=start, count=count )
       elseif( dim_type == 'LZ' ) then
          start(1) = 1
          count(1) = LKMAX
          call FILE_Read( var, fid, varname, step,                            &
                          ntypes=LKMAX, dtype=etype, start=start, count=count )
       elseif( dim_type == 'UZ' ) then
          start(1) = 1
          count(1) = UKMAX
          call FILE_Read( var, fid, varname, step,                            &
                          ntypes=UKMAX, dtype=etype, start=start, count=count )
       elseif( dim_type == 'X' .OR. dim_type == 'CX' ) then
          start(1) = IS_inG - IHALO
          count(1) = IA
          call FILE_Read( var, fid, varname, step,                         &
                          ntypes=IA, dtype=etype, start=start, count=count )
       elseif( dim_type == 'Y' .OR. dim_type == 'CY' ) then
          start(1) = JS_inG - JHALO
          count(1) = JA
          call FILE_Read( var, fid, varname, step,                         &
                          ntypes=JA, dtype=etype, start=start, count=count )
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_1D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif
    else
       if    ( dim_type == 'Z' ) then
          dim1_S   = KS
          dim1_E   = KE
       elseif( dim_type == 'OZ' ) then
          dim1_S   = 1
          dim1_E   = OKMAX
       elseif( dim_type == 'LZ' ) then
          dim1_S   = 1
          dim1_E   = LKMAX
       elseif( dim_type == 'UZ' ) then
          dim1_S   = 1
          dim1_E   = UKMAX
       elseif( dim_type == 'X' ) then
          dim1_S   = ISB
          dim1_E   = IEB
       elseif( dim_type == 'CX' ) then
          dim1_S   = 1
          dim1_E   = IA
       elseif( dim_type == 'Y' ) then
          dim1_S   = JSB
          dim1_E   = JEB
       elseif( dim_type == 'CY' ) then
          dim1_S   = 1
          dim1_E   = JA
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_1D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       call FILE_Read( var(dim1_S:dim1_E), fid, varname, step )
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_1D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILE_CARTESC_read_var_2D( &
       var,      &
       fid,      &
       varname,  &
       dim_type, &
       step      )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Read
    use scale_process, only: &
       PRC_abort
    use scale_rm_process, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    use mpi
    implicit none

    real(RP),         intent(out) :: var(:,:) !< value of the variable
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)
    integer,          intent(in)  :: step     !< step number

    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Read from file (2D), name : ', trim(varname)

    if ( FILE_get_AGGREGATE(fid) ) then
       ! read data and halos into the local buffer
       if    ( dim_type == 'XY' ) then
          call FILE_Read( var, fid, varname, step,                                &
                          ntypes=IA*JA, dtype=etype, start=startXY, count=countXY )
       elseif( dim_type == 'ZX' ) then
          ! Because KHALO is not saved in files, we use centerTypeZX, an MPI
          ! derived datatype to describe the layout of local read buffer
          call FILE_Read( var, fid, varname, step,                                   &
                          ntypes=1, dtype=centerTypeZX, start=startZX, count=countZX )
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_2D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif
    else
       if    ( dim_type == 'XY' ) then
          dim1_S   = ISB
          dim1_E   = IEB
          dim2_S   = JSB
          dim2_E   = JEB
       elseif( dim_type == 'ZX' ) then
          dim1_S   = KS
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_2D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       call FILE_Read( var(dim1_S:dim1_E,dim2_S:dim2_E), fid, varname, step )
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILE_CARTESC_read_var_3D( &
       var,      &
       fid,      &
       varname,  &
       dim_type, &
       step      )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Read
    use scale_process, only: &
       PRC_abort
    use scale_rm_process, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable
    integer,          intent(in)  :: fid        !< file ID
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: dim_type   !< dimension type (Z/X/Y/T)
    integer,          intent(in)  :: step       !< step number

    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    integer :: dim3_S, dim3_E
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Read from file (3D), name : ', trim(varname)

    if ( FILE_get_AGGREGATE(fid) ) then
       ! read data and halos into the local buffer
       ! Because KHALO is not saved in files, we use mpi derived datatypes to
       ! describe the layout of local read buffer
       if(      dim_type == 'ZXY'  &
           .or. dim_type == 'ZXHY' &
           .or. dim_type == 'ZXYH' ) then
          call FILE_Read( var, fid, varname, step,                                      &
                          ntypes=1, dtype=centerTypeZXY, start=startZXY, count=countZXY )
       elseif( dim_type == 'ZHXY' ) then
          call FILE_Read( var, fid, varname, step,                                      &
                          ntypes=1, dtype=centerTypeZHXY, start=startZHXY, count=countZHXY )
       elseif( dim_type == 'XYT' ) then
          startXY(3) = 1
          countXY(3) = step
          call FILE_Read( var, fid, varname, step,                                     &
                          ntypes=step*IA*JA, dtype=etype, start=startXY, count=countXY )
       elseif( dim_type == 'OXY' ) then
          call FILE_Read( var, fid, varname, step,                                         &
                          ntypes=1, dtype=centerTypeOCEAN, start=startOCEAN, count=countOCEAN )
       elseif( dim_type == 'LXY' ) then
          call FILE_Read( var, fid, varname, step,                                         &
                          ntypes=1, dtype=centerTypeLAND, start=startLAND, count=countLAND )
       elseif( dim_type == 'UXY' ) then
          call FILE_Read( var, fid, varname, step,                                            &
                          ntypes=1, dtype=centerTypeURBAN, start=startURBAN, count=countURBAN )
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_3D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif
    else
       if(      dim_type == 'ZXY'  &
           .or. dim_type == 'ZXHY' &
           .or. dim_type == 'ZXYH' ) then
          dim1_S   = KS
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       elseif( dim_type == 'ZHXY' ) then
          dim1_S   = KS-1
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       elseif( dim_type == 'XYT' ) then
          dim1_S   = ISB
          dim1_E   = IEB
          dim2_S   = JSB
          dim2_E   = JEB
          dim3_S   = 1
          dim3_E   = step
       elseif( dim_type == 'OXY' ) then
          dim1_S   = OKS
          dim1_E   = OKE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       elseif( dim_type == 'LXY' ) then
          dim1_S   = LKS
          dim1_E   = LKE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       elseif( dim_type == 'UXY' ) then
          dim1_S   = UKS
          dim1_E   = UKE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_3D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       call FILE_Read( var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                       fid, varname, step                              )
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_3D

  !-----------------------------------------------------------------------------
  !> Read 4D data from file
  subroutine FILE_CARTESC_read_var_4D( &
       var,      &
       fid,      &
       varname,  &
       dim_type, &
       step      )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Read
    use scale_process, only: &
       PRC_abort
    use scale_rm_process, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable
    integer,          intent(in)  :: fid          !< file ID
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: dim_type     !< dimension type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number

    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    integer :: dim3_S, dim3_E
    integer :: dim4_S, dim4_E
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_I_NetCDF', 2)

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Read from file (4D), name : ', trim(varname)

    if ( FILE_get_AGGREGATE(fid) ) then
       ! read data and halos into the local buffer
       if (      dim_type == 'ZXYT'  &
            .or. dim_type == 'ZXHYT' &
            .or. dim_type == 'ZXYHT' ) then
          startZXY(4) = 1
          countZXY(4) = step
          call FILE_Read( var, fid, varname, step,                                         &
                          ntypes=step, dtype=centerTypeZXY, start=startZXY, count=countZXY )
       elseif ( dim_type == 'ZHXYT' ) then
          startZXY(4) = 1
          countZXY(4) = step
          call FILE_Read( var, fid, varname, step,                                         &
                          ntypes=step, dtype=centerTypeZHXY, start=startZHXY, count=countZHXY )
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_4D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif
    else
       if (      dim_type == 'ZXYT'  &
            .or. dim_type == 'ZXHYT' &
            .or. dim_type == 'ZXYHT' ) then
          dim1_S   = KS
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
          dim4_S   = 1
          dim4_E   = step
       elseif ( dim_type == 'ZHXYT' ) then
          dim1_S   = KS-1
          dim1_E   = KE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
          dim4_S   = 1
          dim4_E   = step
       elseif ( dim_type == 'OXYT' ) then
          dim1_S   = OKS
          dim1_E   = OKE
          dim2_S   = ISB
          dim2_E   = IEB
          dim3_S   = JSB
          dim3_E   = JEB
          dim4_S   = 1
          dim4_E   = step
       else
          write(*,*) 'xxx [FILE_CARTESC_read_var_4D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       call FILE_Read( var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,dim4_S:dim4_E), &
                       fid, varname, step                                            )
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_4D

  !-----------------------------------------------------------------------------
  !> interface FILE_CARTESC_write
  !!   Write data to file
  !!   This routine is a wrapper of the lowere primitive routines
  !<
  !-----------------------------------------------------------------------------
  subroutine FILE_CARTESC_write_1D( &
       var,      &
       basename, &
       title,    &
       varname,  &
       desc,     &
       unit,     &
       dim_type, &
       datatype, &
       date,     &
       subsec,   &
       append,   &
       aggregate )
    use scale_process, only: &
       PRC_abort
    implicit none

    real(RP),         intent(in) :: var(:)   !< value of the variable
    character(len=*), intent(in) :: basename !< basename of the file
    character(len=*), intent(in) :: title    !< title    of the file
    character(len=*), intent(in) :: varname  !< name        of the variable
    character(len=*), intent(in) :: desc     !< description of the variable
    character(len=*), intent(in) :: unit     !< unit        of the variable
    character(len=*), intent(in) :: dim_type !< dimension type (Z/X/Y)
    character(len=*), intent(in) :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6) !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec  !< subsec of the time
    logical,          intent(in), optional :: append  !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: aggregate

    integer :: fid, vid
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Write to file (1D), name : ', trim(varname)

    call FILE_CARTESC_create( fid,                            & ! [OUT]
                              basename, title, datatype,      & ! [IN]
                              date, subsec, append, aggregate ) ! [IN]

    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid                                           ) ! [OUT]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_1D( fid, vid, var, varname, dim_type )

    return
  end subroutine FILE_CARTESC_write_1D

  !-----------------------------------------------------------------------------
  !> Write 2D data to file
  subroutine FILE_CARTESC_write_2D( &
       var,       &
       basename,  &
       title,     &
       varname,   &
       desc,      &
       unit,      &
       dim_type,  &
       datatype,  &
       date,      &
       subsec,    &
       append,    &
       fill_halo, &
       haszcoord, &
       aggregate  )
    use scale_process, only: &
       PRC_abort
    implicit none

    real(RP),         intent(in) :: var(:,:) !< value of the variable
    character(len=*), intent(in) :: basename !< basename of the file
    character(len=*), intent(in) :: title    !< title    of the file
    character(len=*), intent(in) :: varname  !< name        of the variable
    character(len=*), intent(in) :: desc     !< description of the variable
    character(len=*), intent(in) :: unit     !< unit        of the variable
    character(len=*), intent(in) :: dim_type !< dimension type (Z/X/Y)
    character(len=*), intent(in) :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6)   !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec    !< subsec of the time
    logical,          intent(in), optional :: append    !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: fill_halo !< switch whether include halo data or not    (default=false)
    logical,          intent(in), optional :: haszcoord !< switch whether include zcoordinate or not  (default=true)
    logical,          intent(in), optional :: aggregate

    integer :: fid, vid
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Write to file (2D), name : ', trim(varname)

    call FILE_CARTESC_create( fid,                                       & ! [OUT]
                              basename, title, datatype,                 & ! [IN]
                              date, subsec, append, haszcoord, aggregate ) ! [IN]

    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid                                           ) ! [OUT]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_2D( fid, vid, var, varname, dim_type, fill_halo )

    return
  end subroutine FILE_CARTESC_write_2D

  !-----------------------------------------------------------------------------
  !> Write 3D data to file
  subroutine FILE_CARTESC_write_3D( &
       var,       &
       basename,  &
       title,     &
       varname,   &
       desc,      &
       unit,      &
       dim_type,  &
       datatype,  &
       date,      &
       subsec,    &
       append,    &
       fill_halo, &
       aggregate  )
    use scale_process, only: &
       PRC_masterrank, &
       PRC_abort
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: basename   !< basename of the file
    character(len=*), intent(in) :: title      !< title    of the file
    character(len=*), intent(in) :: varname    !< name        of the variable
    character(len=*), intent(in) :: desc       !< description of the variable
    character(len=*), intent(in) :: unit       !< unit        of the variable
    character(len=*), intent(in) :: dim_type   !< dimension type (Z/X/Y)
    character(len=*), intent(in) :: datatype   !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6)   !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec    !< subsec of the time
    logical,          intent(in), optional :: append    !< append existing (closed) file?
    logical,          intent(in), optional :: fill_halo !< include halo data?
    logical,          intent(in), optional :: aggregate

    integer :: fid, vid
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Write to file (3D), name : ', trim(varname)

    call FILE_CARTESC_create( fid,                            & ! [OUT]
                              basename, title, datatype,      & ! [IN]
                              date, subsec, append, aggregate ) ! [IN]

    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid                                           ) ! [OUT]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_3D( fid, vid, var, varname, dim_type, fill_halo )

    return
  end subroutine FILE_CARTESC_write_3D

  !-----------------------------------------------------------------------------
  !> Write 3D data with time dimension to file
  subroutine FILE_CARTESC_write_3D_t( &
       var,       &
       basename,  &
       title,     &
       varname,   &
       desc,      &
       unit,      &
       dim_type,  &
       datatype,  &
       timeintv,  &
       tsince,    &
       append,    &
       timetarg,  &
       timeofs,   &
       fill_halo, &
       aggregate  )
    use scale_process, only: &
       PRC_masterrank, &
       PRC_abort
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: basename   !< basename of the file
    character(len=*), intent(in) :: title      !< title    of the file
    character(len=*), intent(in) :: varname    !< name        of the variable
    character(len=*), intent(in) :: desc       !< description of the variable
    character(len=*), intent(in) :: unit       !< unit        of the variable
    character(len=*), intent(in) :: dim_type   !< dimension type (X/Y/Time)
    character(len=*), intent(in) :: datatype   !< data type (REAL8/REAL4/default)
    real(DP),         intent(in) :: timeintv   !< time interval [sec]
    integer ,         intent(in) :: tsince(6)  !< start time

    logical,          intent(in), optional :: append    !< append existing (closed) file?
    integer,          intent(in), optional :: timetarg  !< target timestep (optional)
    real(DP),         intent(in), optional :: timeofs   !< offset time     (optional)
    logical,          intent(in), optional :: fill_halo !< include halo data?
    logical,          intent(in), optional :: aggregate

    integer  :: fid, vid
    integer  :: nsteps

    intrinsic :: size
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,'(1x,3A)') '*** Write to file (3D), name : ', trim(varname), 'with time dimension'

    call FILE_CARTESC_create( fid,                               & ! [OUT]
                              basename, title, datatype, tsince, & ! [IN]
                              append=append, aggregate=aggregate ) ! [IN]

    if ( present(timetarg) ) then
       nsteps = 1
    else
       nsteps = size(var,3)
    endif
    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               timeintv, nsteps                              ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_3D_t( fid, vid, var, varname, dim_type, timeintv, &
                                      timetarg, timeofs, fill_halo                )

    return
  end subroutine FILE_CARTESC_write_3D_t

  !-----------------------------------------------------------------------------
  !> Write 4D data to file
  subroutine FILE_CARTESC_write_4D( &
       var,       &
       basename,  &
       title,     &
       varname,   &
       desc,      &
       unit,      &
       dim_type,  &
       datatype,  &
       timeintv,  &
       tsince,    &
       append,    &
       timetarg,  &
       timeofs,   &
       fill_halo, &
       aggregate  )
    use scale_process, only: &
       PRC_abort
    implicit none

    real(RP),         intent(in) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in) :: basename     !< basename of the file
    character(len=*), intent(in) :: title        !< title    of the file
    character(len=*), intent(in) :: varname      !< name        of the variable
    character(len=*), intent(in) :: desc         !< description of the variable
    character(len=*), intent(in) :: unit         !< unit        of the variable
    character(len=*), intent(in) :: dim_type     !< dimension type (Z/X/Y/Time)
    character(len=*), intent(in) :: datatype     !< data type (REAL8/REAL4/default)
    real(DP),         intent(in) :: timeintv     !< time interval [sec]
    integer,          intent(in) :: tsince(6)    !< start time

    logical,          intent(in), optional :: append    !< append existing (closed) file?
    integer,          intent(in), optional :: timetarg  !< target timestep (optional)
    real(DP),         intent(in), optional :: timeofs   !< offset time     (optional)
    logical,          intent(in), optional :: fill_halo !< include halo data?
    logical,          intent(in), optional :: aggregate

    integer  :: fid, vid
    integer  :: nsteps

    intrinsic :: size
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** Write to file (4D), name : ', trim(varname)

    call FILE_CARTESC_create( fid,                               & ! [OUT]
                              basename, title, datatype, tsince, & ! [IN]
                              append=append, aggregate=aggregate ) ! [IN]

    if ( present(timetarg) ) then
       nsteps = 1
    else
       nsteps = size(var,3)
    endif
    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               timeintv, nsteps                              ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_4D( fid, vid, var, varname, dim_type, timeintv, &
                                    timetarg, timeofs, fill_halo                )

    return
  end subroutine FILE_CARTESC_write_4D

  !-----------------------------------------------------------------------------
  !> open a netCDF file for read
  subroutine FILE_CARTESC_open( &
       fid,      &
       basename, &
       aggregate )
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_AGGREGATE, &
       FILE_Open
    use scale_process, only: &
       PRC_myrank, &
       PRC_LOCAL_COMM_WORLD
    use mpi, only : MPI_COMM_NULL
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of the file
    logical,          intent(in), optional :: aggregate

    integer :: comm
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    comm = MPI_COMM_NULL
    ! check to do PnetCDF I/O
    if ( present(aggregate) ) then
       if (aggregate) comm = PRC_LOCAL_COMM_WORLD
    else if ( FILE_AGGREGATE ) then
       comm = PRC_LOCAL_COMM_WORLD
    end if

    call FILE_Open( fid,                & ! [OUT]
                    basename,           & ! [IN]
                    FILE_FREAD,         & ! [IN]
                    mpi_comm = comm,    & ! [IN]
                    rankid = PRC_myrank ) ! [IN]

    File_closed(fid) = .false.

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_open

  !-----------------------------------------------------------------------------
  !> Create/open a netCDF file
  subroutine FILE_CARTESC_create( &
       fid,       &
       basename,  &
       title,     &
       datatype,  &
       date,      &
       subsec,    &
       append,    &
       haszcoord, &
       aggregate  )
    use mpi, only: &
       MPI_COMM_NULL
    use scale_file_h, only: &
       FILE_REAL8, &
       FILE_REAL4
    use scale_file, only: &
       FILE_AGGREGATE,    &
       FILE_Create,       &
       FILE_get_CFtunits, &
       FILE_Set_Attribute
    use scale_process, only: &
       PRC_Ismaster, &
       PRC_myrank,   &
       PRC_abort,  &
       PRC_LOCAL_COMM_WORLD
    use scale_rm_process, only: &
       PRC_2Drank,     &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE, &
       NOWMS   => TIME_NOWMS
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(in), optional :: date(6)   !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec    !< subsec of the time
    logical,          intent(in), optional :: append    !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: haszcoord !< switch whether include zcoordinate or not (default=true)
    logical,          intent(in), optional :: aggregate

    character(len=5)  :: periodic_z, periodic_x, periodic_y
    integer           :: dtype
    logical           :: append_sw
    character(len=34) :: tunits
    integer           :: comm
    logical           :: fileexisted
    logical           :: aggregate_
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    ! dtype is used to define the data type of axis variables in file
    if    ( datatype == 'REAL8' ) then
       dtype = FILE_REAL8
    elseif( datatype == 'REAL4' ) then
       dtype = FILE_REAL4
    else
       if    ( RP == 8 ) then
          dtype = FILE_REAL8
       elseif( RP == 4 ) then
          dtype = FILE_REAL4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype)
          call PRC_abort
       endif
    endif

    append_sw = .false.
    if ( present(append) ) then
       append_sw = append
    endif

    ! create a netCDF file if not already existed. Otherwise, open it.
    if ( present(date) ) then
       call FILE_get_CFtunits( tunits, date )
    else
       tunits = 'seconds'
    endif


    ! check to use PnetCDF I/O
    if ( present(aggregate) ) then
       aggregate_ = aggregate
    else
       aggregate_ = FILE_AGGREGATE
    endif

    if ( aggregate_ ) then
       comm = PRC_LOCAL_COMM_WORLD
    else
       comm = MPI_COMM_NULL
    end if

    call FILE_Create( fid,                     & ! [OUT]
                      fileexisted,             & ! [OUT]
                      basename,                & ! [IN]
                      title,                   & ! [IN]
                      H_SOURCE,                & ! [IN]
                      H_INSTITUTE,             & ! [IN]
                      "cartesC",               & ! [IN]
!                      GRID_CARTESC_NAME,       & ! [IN]
                      rankid = PRC_myrank,     & ! [IN]
                      ismaster = PRC_Ismaster, & ! [IN]
                      time_units = tunits,     & ! [IN]
                      append     = append_sw,  & ! [IN]
                      mpi_comm   = comm        ) ! [IN]



    if ( .NOT. fileexisted ) then ! do below only once when file is created

       if ( present( haszcoord ) ) then
          File_haszcoord(fid) = haszcoord
       else
          File_haszcoord(fid) = .true.
       endif

       periodic_z = "false"
       if ( PRC_PERIODIC_X ) then
          periodic_x = "true"
       else
          periodic_x = "false"
       endif
       if ( PRC_PERIODIC_Y ) then
          periodic_y = "true"
       else
          periodic_y = "false"
       endif

       if ( aggregate_ ) then
          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_rank_x", (/0/) ) ! [IN]
          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_rank_y", (/0/) ) ! [IN]

          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_num_x",  (/1/) ) ! [IN]
          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_num_y",  (/1/) ) ! [IN]
       else
          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_rank_x", (/PRC_2Drank(PRC_myrank,1)/) ) ! [IN]
          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_rank_y", (/PRC_2Drank(PRC_myrank,2)/) ) ! [IN]

          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_num_x",  (/PRC_NUM_X/) ) ! [IN]
          call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_num_y",  (/PRC_NUM_Y/) ) ! [IN]
       endif

       call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_periodic_z", periodic_z ) ! [IN]
       call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_periodic_x", periodic_x ) ! [IN]
       call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_periodic_y", periodic_y ) ! [IN]

       call FILE_Set_Attribute( fid, "global", "scale_cartesC_grid_index_kmax",  (/KMAX/)  ) ! [IN]
       call FILE_Set_Attribute( fid, "global", "scale_cartesC_grid_index_imaxg", (/IMAXG/) ) ! [IN]
       call FILE_Set_Attribute( fid, "global", "scale_cartesC_grid_index_jmaxg", (/JMAXG/) ) ! [IN]

       call FILE_Set_Attribute( fid, "global", "scale_cartesC_grid_index_khalo", (/KHALO/) ) ! [IN]
       call FILE_Set_Attribute( fid, "global", "scale_cartesC_grid_index_ihalo", (/IHALO/) ) ! [IN]
       call FILE_Set_Attribute( fid, "global", "scale_cartesC_grid_index_jhalo", (/JHALO/) ) ! [IN]

       call FILE_CARTESC_def_axes( fid,                & ! [IN]
                             dtype,              & ! [IN]
                             File_haszcoord(fid) ) ! [IN]

       if ( present( date ) ) then
          call FILE_get_CFtunits(tunits, date)
       else
          call FILE_get_CFtunits(tunits, NOWDATE)
       endif
       call FILE_Set_Attribute( fid, "global", "time_units", tunits )

       if ( present( subsec ) ) then
          call FILE_Set_Attribute( fid, "global", "time_start", (/subsec/) )
       else
          call FILE_Set_Attribute( fid, "global", "time_start", (/NOWMS/)  )
       endif

       File_axes_written(fid) = .false.  ! indicating axes have not been written yet
       File_closed      (fid) = .false.
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF file define mode
  subroutine FILE_CARTESC_enddef( &
       fid )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_EndDef,    &
       FILE_Flush,     &
       FILE_Attach_Buffer
    use scale_process, only: &
       PRC_myrank
    implicit none

    integer, intent(in) :: fid  !< file ID

    integer :: start(3)
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call FILE_EndDef( fid ) ! [IN]

    ! If this enddef is called the first time, write axis variables
    if ( .NOT. File_axes_written(fid) ) then

       if ( FILE_get_AGGREGATE(fid) ) then
          start(1) = 1
          start(2) = ISGA
          start(3) = JSGA
       else
          start(1) = 1
          start(2) = 1
          start(3) = 1
       endif

       call FILE_CARTESC_write_axes( fid,                 & ! [IN]
                                     File_haszcoord(fid), & ! [IN]
                                     start(:)             ) ! [IN]

       ! Tell PnetCDF library to use a buffer of size write_buf_amount to aggregate write requests to be post in FILE_CARTESC_write_var
       if ( FILE_get_AGGREGATE(fid) ) then
          call FILE_Flush( fid )
          call FILE_Attach_Buffer( fid, write_buf_amount )
       endif

       File_axes_written(fid) = .true.
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_enddef

  !-----------------------------------------------------------------------------
  !> Flush all pending requests to a netCDF file (PnetCDF only)
  subroutine FILE_CARTESC_flush( &
       fid )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Flush
    implicit none

    integer, intent(in) :: fid  !< file ID
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( FILE_get_AGGREGATE(fid) ) then
       call FILE_Flush( fid ) ! flush all pending read/write requests
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_flush

  !-----------------------------------------------------------------------------
  !> Close a netCDF file
  subroutine FILE_CARTESC_close( &
       fid )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Close,     &
       FILE_Flush,     &
       FILE_Detach_Buffer
    implicit none

    integer, intent(in) :: fid  !< file ID

    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( .NOT. File_closed(fid) ) then

       if ( FILE_get_AGGREGATE(fid) ) then
          call FILE_Flush( fid )        ! flush all pending read/write requests
          if ( write_buf_amount > 0 ) then
             call FILE_Detach_Buffer( fid ) ! detach PnetCDF aggregation buffer
             write_buf_amount = 0         ! reset write request amount
          endif
       endif

       call FILE_Close( fid ) ! [IN]

       File_closed(fid) = .true.

    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_close

  !-----------------------------------------------------------------------------
  !> define axis variables in the file
  subroutine FILE_CARTESC_def_axes( &
       fid,   &
       dtype, &
       hasZ   )
    use scale_file, only: &
       FILE_get_AGGREGATE,            &
       FILE_Def_Axis,                 &
       FILE_Set_Attribute,            &
       FILE_Def_AssociatedCoordinate, &
       FILE_Add_AssociatedVariable
    use scale_const, &
       UNDEF => CONST_UNDEF
    use scale_rm_process, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_HAS_W,      &
       PRC_HAS_E,      &
       PRC_HAS_S,      &
       PRC_HAS_N
    use scale_mapproj, only: &
       MPRJ_get_attributes
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: dtype
    logical, intent(in) :: hasZ

    integer :: iall  ! grid size, x-axis, (whole domain or local tile), including halo
    integer :: jall  ! grid size, y-axis, (whole domain or local tile), including halo
    integer :: isize ! grid size, x-axis, (whole domain or local tile), without halo except domain edge
    integer :: jsize ! grid size, y-axis, (whole domain or local tile), without halo except domain edge

    type(axisattinfo) :: ainfo(4) ! x, xh, y, yh
    type(mappinginfo) :: minfo

    character(len=2) :: axisname(3)
    !---------------------------------------------------------------------------

    if ( PRC_PERIODIC_X ) then
       ainfo(1)%periodic = "true"
       ainfo(2)%periodic = "true"
    else
       ainfo(1)%periodic = "false"
       ainfo(2)%periodic = "false"
    endif

    if ( PRC_PERIODIC_Y ) then
       ainfo(3)%periodic = "true"
       ainfo(4)%periodic = "true"
    else
       ainfo(3)%periodic = "false"
       ainfo(4)%periodic = "false"
    endif

    call MPRJ_get_attributes( minfo%mapping_name,                             & ! [OUT]
                              minfo%false_easting                        (1), & ! [OUT]
                              minfo%false_northing                       (1), & ! [OUT]
                              minfo%longitude_of_central_meridian        (1), & ! [OUT]
                              minfo%longitude_of_projection_origin       (1), & ! [OUT]
                              minfo%latitude_of_projection_origin        (1), & ! [OUT]
                              minfo%straight_vertical_longitude_from_pole(1), & ! [OUT]
                              minfo%standard_parallel                    (:)  ) ! [OUT]

    if ( FILE_get_AGGREGATE(fid) ) then
       ! for x = CXG
       ainfo(1)%size_global (1) = IAG
       ainfo(1)%start_global(1) = 1
       ainfo(1)%halo_global (1) = IHALO ! west side
       ainfo(1)%halo_global (2) = IHALO ! east side
       ainfo(1)%halo_local  (1) = IHALO ! west side
       ainfo(1)%halo_local  (2) = IHALO ! east side
       ! for xh = FXG
       ainfo(2)%size_global (1) = IAG+1
       ainfo(2)%start_global(1) = 1
       ainfo(2)%halo_global (1) = IHALO ! west side
       ainfo(2)%halo_global (2) = IHALO ! east side
       ainfo(2)%halo_local  (1) = IHALO ! west side
       ainfo(2)%halo_local  (2) = IHALO ! east side
       ! for y = CYG
       ainfo(3)%size_global (1) = JAG
       ainfo(3)%start_global(1) = 1
       ainfo(3)%halo_global (1) = JHALO ! south side
       ainfo(3)%halo_global (2) = JHALO ! north side
       ainfo(3)%halo_local  (1) = JHALO ! south side
       ainfo(3)%halo_local  (2) = JHALO ! north side
       ! for yh = FYG
       ainfo(4)%size_global (1) = JAG+1
       ainfo(4)%start_global(1) = 1
       ainfo(4)%halo_global (1) = JHALO ! south side
       ainfo(4)%halo_global (2) = JHALO ! north side
       ainfo(4)%halo_local  (1) = JHALO ! south side
       ainfo(4)%halo_local  (2) = JHALO ! north side
    else
       ! for x
       if ( PRC_PERIODIC_X ) then
          ainfo(1)%size_global (1) = IMAX * PRC_NUM_X
          ainfo(1)%start_global(1) = IS_inG - IHALO
          ainfo(1)%halo_global (1) = 0     ! west side
          ainfo(1)%halo_global (2) = 0     ! east side
          ainfo(1)%halo_local  (1) = 0     ! west side
          ainfo(1)%halo_local  (2) = 0     ! east side
       else
          ainfo(1)%size_global (1) = IAG
          ainfo(1)%start_global(1) = ISGA
          ainfo(1)%halo_global (1) = IHALO ! west side
          ainfo(1)%halo_global (2) = IHALO ! east side
          ainfo(1)%halo_local  (1) = IHALO ! west side
          ainfo(1)%halo_local  (2) = IHALO ! east side
          if( PRC_HAS_W ) ainfo(1)%halo_local(1) = 0
          if( PRC_HAS_E ) ainfo(1)%halo_local(2) = 0
       endif
       ! for xh
       ainfo(2) = ainfo(1)
       ! for y
       if ( PRC_PERIODIC_Y ) then
          ainfo(3)%size_global (1) = JMAX * PRC_NUM_Y
          ainfo(3)%start_global(1) = JS_inG - JHALO
          ainfo(3)%halo_global (1) = 0     ! south side
          ainfo(3)%halo_global (2) = 0     ! north side
          ainfo(3)%halo_local  (1) = 0     ! south side
          ainfo(3)%halo_local  (2) = 0     ! north side
       else
          ainfo(3)%size_global (1) = JAG
          ainfo(3)%start_global(1) = JSGA
          ainfo(3)%halo_global (1) = JHALO ! south side
          ainfo(3)%halo_global (2) = JHALO ! north side
          ainfo(3)%halo_local  (1) = JHALO ! south side
          ainfo(3)%halo_local  (2) = JHALO ! north side
          if( PRC_HAS_S ) ainfo(3)%halo_local(1) = 0
          if( PRC_HAS_N ) ainfo(3)%halo_local(2) = 0
       endif
       ! for yh
       ainfo(4) = ainfo(3)
    endif

    if ( FILE_get_AGGREGATE(fid) ) then
       iall  = IAG
       jall  = JAG
       isize = IAG
       jsize = JAG
    else
       iall  = IA
       jall  = JA
       isize = IMAXB
       jsize = JMAXB
    endif

    if( hasZ ) call FILE_Def_Axis( fid, 'z'  , 'Z'              , 'm', 'z'  , dtype, KMAX    )
    if( hasZ ) call FILE_Def_Axis( fid, 'zh' , 'Z (half level)' , 'm', 'zh' , dtype, KMAX+1  )

    if( hasZ ) call FILE_Def_Axis( fid, 'oz' , 'OZ'             , 'm', 'oz' , dtype, OKMAX   )
    if( hasZ ) call FILE_Def_Axis( fid, 'ozh', 'OZ (half level)', 'm', 'ozh', dtype, OKMAX+1 )

    if( hasZ ) call FILE_Def_Axis( fid, 'lz' , 'LZ'             , 'm', 'lz' , dtype, LKMAX   )
    if( hasZ ) call FILE_Def_Axis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, LKMAX+1 )

    if( hasZ ) call FILE_Def_Axis( fid, 'uz' , 'UZ'             , 'm', 'uz' , dtype, UKMAX   )
    if( hasZ ) call FILE_Def_Axis( fid, 'uzh', 'UZ (half level)', 'm', 'uzh', dtype, UKMAX+1 )

               call FILE_Def_Axis( fid, 'x'  , 'X'              , 'm', 'x'  , dtype, isize   )
               call FILE_Def_Axis( fid, 'xh' , 'X (half level)' , 'm', 'xh' , dtype, isize   )
               call FILE_Def_Axis( fid, 'y'  , 'Y'              , 'm', 'y'  , dtype, jsize   )
               call FILE_Def_Axis( fid, 'yh' , 'Y (half level)' , 'm', 'yh' , dtype, jsize   )

    if( hasZ ) call FILE_Def_Axis( fid, 'CZ'   , 'Atmos Grid Center Position Z',      'm', 'CZ',   dtype, KA      )
    if( hasZ ) call FILE_Def_Axis( fid, 'FZ'   , 'Atmos Grid Face   Position Z',      'm', 'FZ',   dtype, KA+1    )
    if( hasZ ) call FILE_Def_Axis( fid, 'CDZ'  , 'Grid Cell length Z',                'm', 'CZ',   dtype, KA      )
    if( hasZ ) call FILE_Def_Axis( fid, 'FDZ'  , 'Grid distance Z',                   'm', 'FDZ',  dtype, KA-1    )
    if( hasZ ) call FILE_Def_Axis( fid, 'CBFZ' , 'Boundary factor Center Z',          '1', 'CZ',   dtype, KA      )
    if( hasZ ) call FILE_Def_Axis( fid, 'FBFZ' , 'Boundary factor Face Z',            '1', 'FZ',   dtype, KA+1    )

    if( hasZ ) call FILE_Def_Axis( fid, 'OCZ'  , 'Ocean Grid Center Position Z',      'm', 'OCZ',  dtype, OKMAX   )
    if( hasZ ) call FILE_Def_Axis( fid, 'OFZ'  , 'Ocean Grid Face   Position Z',      'm', 'OFZ',  dtype, OKMAX+1 )
    if( hasZ ) call FILE_Def_Axis( fid, 'OCDZ' , 'Ocean Grid Cell length Z',          'm', 'OCZ',  dtype, OKMAX   )

    if( hasZ ) call FILE_Def_Axis( fid, 'LCZ'  , 'Land Grid Center Position Z',       'm', 'LCZ',  dtype, LKMAX   )
    if( hasZ ) call FILE_Def_Axis( fid, 'LFZ'  , 'Land Grid Face   Position Z',       'm', 'LFZ',  dtype, LKMAX+1 )
    if( hasZ ) call FILE_Def_Axis( fid, 'LCDZ' , 'Land Grid Cell length Z',           'm', 'LCZ',  dtype, LKMAX   )

    if( hasZ ) call FILE_Def_Axis( fid, 'UCZ'  , 'Urban Grid Center Position Z',      'm', 'UCZ',  dtype, UKMAX   )
    if( hasZ ) call FILE_Def_Axis( fid, 'UFZ'  , 'Urban Grid Face   Position Z',      'm', 'UFZ',  dtype, UKMAX+1 )
    if( hasZ ) call FILE_Def_Axis( fid, 'UCDZ' , 'Urban Grid Cell length Z',          'm', 'UCZ',  dtype, UKMAX   )

               call FILE_Def_Axis( fid, 'CX'   , 'Atmos Grid Center Position X',      'm', 'CX',   dtype, iall   )
               call FILE_Def_Axis( fid, 'CY'   , 'Atmos Grid Center Position Y',      'm', 'CY',   dtype, jall   )
               call FILE_Def_Axis( fid, 'FX'   , 'Atmos Grid Face   Position X',      'm', 'FX',   dtype, iall+1 )
               call FILE_Def_Axis( fid, 'FY'   , 'Atmos Grid Face   Position Y',      'm', 'FY',   dtype, jall+1 )
               call FILE_Def_Axis( fid, 'CDX'  , 'Grid Cell length X',                'm', 'CX',   dtype, iall   )
               call FILE_Def_Axis( fid, 'CDY'  , 'Grid Cell length Y',                'm', 'CY',   dtype, jall   )
               call FILE_Def_Axis( fid, 'FDX'  , 'Grid distance X',                   'm', 'FDX',  dtype, iall-1 )
               call FILE_Def_Axis( fid, 'FDY'  , 'Grid distance Y',                   'm', 'FDY',  dtype, jall-1 )
               call FILE_Def_Axis( fid, 'CBFX' , 'Boundary factor Center X',          '1', 'CX',   dtype, iall   )
               call FILE_Def_Axis( fid, 'CBFY' , 'Boundary factor Center Y',          '1', 'CY',   dtype, jall   )
               call FILE_Def_Axis( fid, 'FBFX' , 'Boundary factor Face X',            '1', 'FX',   dtype, iall+1 )
               call FILE_Def_Axis( fid, 'FBFY' , 'Boundary factor Face Y',            '1', 'FY',   dtype, jall+1 )

               call FILE_Def_Axis( fid, 'CXG'  , 'Grid Center Position X (global)',   'm', 'CXG',  dtype, IAG   )
               call FILE_Def_Axis( fid, 'CYG'  , 'Grid Center Position Y (global)',   'm', 'CYG',  dtype, JAG   )
               call FILE_Def_Axis( fid, 'FXG'  , 'Grid Face   Position X (global)',   'm', 'FXG',  dtype, IAG+1 )
               call FILE_Def_Axis( fid, 'FYG'  , 'Grid Face   Position Y (global)',   'm', 'FYG',  dtype, JAG+1 )
               call FILE_Def_Axis( fid, 'CDXG' , 'Grid Cell length X (global)',       'm', 'CXG',  dtype, IAG   )
               call FILE_Def_Axis( fid, 'CDYG' , 'Grid Cell length Y (global)',       'm', 'CYG',  dtype, JAG   )
               call FILE_Def_Axis( fid, 'FDXG' , 'Grid distance X (global)',          'm', 'FDXG', dtype, IAG-1 )
               call FILE_Def_Axis( fid, 'FDYG' , 'Grid distance Y (global)',          'm', 'FDYG', dtype, JAG-1 )
               call FILE_Def_Axis( fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG',  dtype, IAG   )
               call FILE_Def_Axis( fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG',  dtype, JAG   )
               call FILE_Def_Axis( fid, 'FBFXG', 'Boundary factor Face   X (global)', '1', 'FXG',  dtype, IAG+1 )
               call FILE_Def_Axis( fid, 'FBFYG', 'Boundary factor Face   Y (global)', '1', 'FYG',  dtype, JAG+1 )

    ! associate coordinates
    axisname(1:2) = (/'x ','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'lon'   , 'longitude',                 'degrees_east' , axisname(1:2), dtype )
    axisname(1:2) = (/'xh','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'lon_uy', 'longitude (half level uy)', 'degrees_east' , axisname(1:2), dtype )
    axisname(1:2) = (/'x ','yh'/)
    call FILE_Def_AssociatedCoordinate( fid, 'lon_xv', 'longitude (half level xv)', 'degrees_east' , axisname(1:2), dtype )
    axisname(1:2) = (/'xh','yh'/)
    call FILE_Def_AssociatedCoordinate( fid, 'lon_uv', 'longitude (half level uv)', 'degrees_east' , axisname(1:2), dtype )
    axisname(1:2) = (/'x ','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'lat'   , 'latitude',                  'degrees_north', axisname(1:2), dtype )
    axisname(1:2) = (/'xh','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'lat_uy', 'latitude (half level uy)',  'degrees_north', axisname(1:2), dtype )
    axisname(1:2) = (/'x ','yh'/)
    call FILE_Def_AssociatedCoordinate( fid, 'lat_xv', 'latitude (half level xv)',  'degrees_north', axisname(1:2), dtype )
    axisname(1:2) = (/'xh','yh'/)
    call FILE_Def_AssociatedCoordinate( fid, 'lat_uv', 'latitude (half level uv)',  'degrees_north', axisname(1:2), dtype )

    if ( hasZ ) then
       axisname = (/'z ', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'height',     'height above ground level', &
                                          'm', axisname(1:3), dtype                       )
       axisname = (/'zh', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'height_wxy', 'height above ground level (half level wxy)', &
                                          'm', axisname(1:3), dtype                                        )
    endif

    ! attributes

    if ( hasZ ) then
       call FILE_Set_Attribute( fid, 'oz' , 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'ozh', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'lz' , 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'lzh', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'uz' , 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'uzh', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'OCZ', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'OFZ', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'LCZ', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'LFZ', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'UCZ', 'positive', 'down' )
       call FILE_Set_Attribute( fid, 'UFZ', 'positive', 'down' )
    endif

    call FILE_Set_Attribute( fid, "x" , "size_global" , ainfo(1)%size_global (:) )
    call FILE_Set_Attribute( fid, "x" , "start_global", ainfo(1)%start_global(:) )
    call FILE_Set_Attribute( fid, "x" , "halo_global" , ainfo(1)%halo_global (:) )
    call FILE_Set_Attribute( fid, "x" , "halo_local"  , ainfo(1)%halo_local  (:) )
    call FILE_Set_Attribute( fid, "x" , "periodic"    , ainfo(1)%periodic        )

    call FILE_Set_Attribute( fid, "xh", "size_global" , ainfo(2)%size_global (:) )
    call FILE_Set_Attribute( fid, "xh", "start_global", ainfo(2)%start_global(:) )
    call FILE_Set_Attribute( fid, "xh", "halo_global" , ainfo(2)%halo_global (:) )
    call FILE_Set_Attribute( fid, "xh", "halo_local"  , ainfo(2)%halo_local  (:) )
    call FILE_Set_Attribute( fid, "xh", "periodic"    , ainfo(2)%periodic        )

    call FILE_Set_Attribute( fid, "y" , "size_global" , ainfo(3)%size_global (:) )
    call FILE_Set_Attribute( fid, "y" , "start_global", ainfo(3)%start_global(:) )
    call FILE_Set_Attribute( fid, "y" , "halo_global" , ainfo(3)%halo_global (:) )
    call FILE_Set_Attribute( fid, "y" , "halo_local"  , ainfo(3)%halo_local  (:) )
    call FILE_Set_Attribute( fid, "y" , "periodic"    , ainfo(3)%periodic        )

    call FILE_Set_Attribute( fid, "yh", "size_global" , ainfo(4)%size_global (:) )
    call FILE_Set_Attribute( fid, "yh", "start_global", ainfo(4)%start_global(:) )
    call FILE_Set_Attribute( fid, "yh", "halo_global" , ainfo(4)%halo_global (:) )
    call FILE_Set_Attribute( fid, "yh", "halo_local"  , ainfo(4)%halo_local  (:) )
    call FILE_Set_Attribute( fid, "yh", "periodic"    , ainfo(4)%periodic        )

    ! map projection info

    if ( minfo%mapping_name /= "" ) then
       call FILE_Set_Attribute( fid, "x" , "standard_name", "projection_x_coordinate" )
       call FILE_Set_Attribute( fid, "xh", "standard_name", "projection_x_coordinate" )
       call FILE_Set_Attribute( fid, "y" , "standard_name", "projection_y_coordinate" )
       call FILE_Set_Attribute( fid, "yh", "standard_name", "projection_y_coordinate" )

       call FILE_Add_AssociatedVariable( fid, minfo%mapping_name )
       call FILE_Set_Attribute( fid, minfo%mapping_name, "grid_mapping_name",  minfo%mapping_name )

       if ( minfo%false_easting(1) /= UNDEF ) then
          call FILE_Set_Attribute( fid,                   & ! [IN]
                                   minfo%mapping_name,    & ! [IN]
                                   "false_easting",       & ! [IN]
                                   minfo%false_easting(:) ) ! [IN]
       endif

       if ( minfo%false_northing(1) /= UNDEF ) then
          call FILE_Set_Attribute( fid,                    & ! [IN]
                                   minfo%mapping_name,     & ! [IN]
                                   "false_northing",       & ! [IN]
                                   minfo%false_northing(:) ) ! [IN]
       endif

       if ( minfo%longitude_of_central_meridian(1) /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                   & ! [IN]
                                   minfo%mapping_name,                    & ! [IN]
                                   "longitude_of_central_meridian",       & ! [IN]
                                   minfo%longitude_of_central_meridian(:) ) ! [IN]
       endif

       if ( minfo%longitude_of_projection_origin(1) /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                    & ! [IN]
                                   minfo%mapping_name,                     & ! [IN]
                                   "longitude_of_projection_origin",       & ! [IN]
                                   minfo%longitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( minfo%latitude_of_projection_origin(1) /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                   & ! [IN]
                                   minfo%mapping_name,                    & ! [IN]
                                   "latitude_of_projection_origin",       & ! [IN]
                                   minfo%latitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( minfo%straight_vertical_longitude_from_pole(1) /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                           & ! [IN]
                                   minfo%mapping_name,                            & ! [IN]
                                   "straight_vertical_longitude_from_pole",       & ! [IN]
                                   minfo%straight_vertical_longitude_from_pole(:) ) ! [IN]
       endif

       if ( minfo%standard_parallel(1) /= UNDEF ) then
          if ( minfo%standard_parallel(2) /= UNDEF ) then
             call FILE_Set_Attribute( fid,                         & ! [IN]
                                      minfo%mapping_name,          & ! [IN]
                                      "standard_parallel",         & ! [IN]
                                      minfo%standard_parallel(1:2) ) ! [IN]
          else
             call FILE_Set_Attribute( fid,                         & ! [IN]
                                      minfo%mapping_name,          & ! [IN]
                                      "standard_parallel",         & ! [IN]
                                      minfo%standard_parallel(1:1) ) ! [IN]
          endif
       endif
    endif

    return
  end subroutine FILE_CARTESC_def_axes

  !-----------------------------------------------------------------------------
  !> write axis to the file
  subroutine FILE_CARTESC_write_axes( &
       fid,       &
       haszcoord, &
       start      )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Write_Axis,                 &
       FILE_Write_AssociatedCoordinate
    use scale_process, only: &
       PRC_myrank, &
       PRC_IsMaster
    use scale_rm_process, only: &
       PRC_2Drank
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
    use scale_ocean_grid, only: &
       GRID_OCZ,  &
       GRID_OFZ,  &
       GRID_OCDZ
    use scale_land_grid, only: &
       GRID_LCZ,  &
       GRID_LFZ,  &
       GRID_LCDZ
    use scale_urban_grid, only: &
       GRID_UCZ,  &
       GRID_UFZ,  &
       GRID_UCDZ
    implicit none

    integer, intent(in)  :: fid
    logical, intent(in)  :: haszcoord
    integer, intent(in)  :: start(3)

    logical :: put_z, put_x, put_y
    integer :: XSB, XEB, YSB, YEB
    !---------------------------------------------------------------------------

    if ( FILE_get_AGGREGATE(fid) ) then
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

       put_z = ( PRC_IsMaster                  ) ! only master process output the vertical coordinates
       put_x = ( PRC_2Drank(PRC_myrank,2) == 0 ) ! only south-most row    processes write x coordinates
       put_y = ( PRC_2Drank(PRC_myrank,1) == 0 ) ! only west-most  column processes write y coordinates
    else
       put_z = .true.
       put_x = .true.
       put_y = .true.
    end if

    if ( haszcoord .and. put_z ) then
       call FILE_Write_Axis( fid, 'z'  , GRID_CZ (KS  :KE)  , start(1:1) )
       call FILE_Write_Axis( fid, 'zh' , GRID_FZ (KS-1:KE)  , start(1:1) )

       call FILE_Write_Axis( fid, 'oz' , GRID_OCZ(OKS  :OKE), start(1:1) )
       call FILE_Write_Axis( fid, 'ozh', GRID_OFZ(OKS-1:OKE), start(1:1) )

       call FILE_Write_Axis( fid, 'lz' , GRID_LCZ(LKS  :LKE), start(1:1) )
       call FILE_Write_Axis( fid, 'lzh', GRID_LFZ(LKS-1:LKE), start(1:1) )

       call FILE_Write_Axis( fid, 'uz' , GRID_UCZ(UKS  :UKE), start(1:1) )
       call FILE_Write_Axis( fid, 'uzh', GRID_UFZ(UKS-1:UKE), start(1:1) )
    end if

    if ( put_x ) then
       if ( FILE_get_aggregate(fid) ) then
          call FILE_Write_Axis( fid, 'x' ,  GRID_CX(ISB2:IEB2),  start(2:2) )
          call FILE_Write_Axis( fid, 'xh',  GRID_FX(ISB2:IEB2),  start(2:2) )
       else
          call FILE_Write_Axis( fid, 'x' ,  GRID_CX(ISB:IEB),  start(2:2) )
          call FILE_Write_Axis( fid, 'xh',  GRID_FX(ISB:IEB),  start(2:2) )
       end if
    end if

    if ( put_y ) then
       if ( FILE_get_aggregate(fid) ) then
          call FILE_Write_Axis( fid, 'y' ,  GRID_CY(JSB2:JEB2),  start(3:3) )
          call FILE_Write_Axis( fid, 'yh',  GRID_FY(JSB2:JEB2),  start(3:3) )
       else
          call FILE_Write_Axis( fid, 'y' ,  GRID_CY(JSB:JEB),  start(3:3) )
          call FILE_Write_Axis( fid, 'yh',  GRID_FY(JSB:JEB),  start(3:3) )
       end if
    end if

    ! global coordinates (always including halo)
    if ( haszcoord .and. put_z ) then
       call FILE_Write_Axis( fid, 'CZ'  , GRID_CZ  (:), start(1:1) )
       call FILE_Write_Axis( fid, 'FZ'  , GRID_FZ  (:), start(1:1) )
       call FILE_Write_Axis( fid, 'CDZ' , GRID_CDZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'FDZ' , GRID_FDZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'CBFZ', GRID_CBFZ(:), start(1:1) )
       call FILE_Write_Axis( fid, 'FBFZ', GRID_FBFZ(:), start(1:1) )

       call FILE_Write_Axis( fid, 'OCZ' , GRID_OCZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'OFZ' , GRID_OFZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'OCDZ', GRID_OCDZ(:), start(1:1) )

       call FILE_Write_Axis( fid, 'LCZ' , GRID_LCZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'LFZ' , GRID_LFZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'LCDZ', GRID_LCDZ(:), start(1:1) )

       call FILE_Write_Axis( fid, 'UCZ' , GRID_UCZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'UFZ' , GRID_UFZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'UCDZ', GRID_UCDZ(:), start(1:1) )
    end if

    if ( FILE_get_AGGREGATE(fid) ) then
       if ( PRC_IsMaster ) then
          call FILE_Write_Axis( fid, 'CX',   GRID_CXG  (:) )
          call FILE_Write_Axis( fid, 'CY',   GRID_CYG  (:) )
          call FILE_Write_Axis( fid, 'FX',   GRID_FXG  (:) )
          call FILE_Write_Axis( fid, 'FY',   GRID_FYG  (:) )
          call FILE_Write_Axis( fid, 'CDX',  GRID_CDXG (:) )
          call FILE_Write_Axis( fid, 'CDY',  GRID_CDYG (:) )
          call FILE_Write_Axis( fid, 'FDX',  GRID_FDXG (:) )
          call FILE_Write_Axis( fid, 'FDY',  GRID_FDYG (:) )
          call FILE_Write_Axis( fid, 'CBFX', GRID_CBFXG(:) )
          call FILE_Write_Axis( fid, 'CBFY', GRID_CBFYG(:) )
          call FILE_Write_Axis( fid, 'FBFX', GRID_FBFXG(:) )
          call FILE_Write_Axis( fid, 'FBFY', GRID_FBFYG(:) )
       endif
    else
       call FILE_Write_Axis( fid, 'CX',   GRID_CX  (:) )
       call FILE_Write_Axis( fid, 'CY',   GRID_CY  (:) )
       call FILE_Write_Axis( fid, 'FX',   GRID_FX  (:) )
       call FILE_Write_Axis( fid, 'FY',   GRID_FY  (:) )
       call FILE_Write_Axis( fid, 'CDX',  GRID_CDX (:) )
       call FILE_Write_Axis( fid, 'CDY',  GRID_CDY (:) )
       call FILE_Write_Axis( fid, 'FDX',  GRID_FDX (:) )
       call FILE_Write_Axis( fid, 'FDY',  GRID_FDY (:) )
       call FILE_Write_Axis( fid, 'CBFX', GRID_CBFX(:) )
       call FILE_Write_Axis( fid, 'CBFY', GRID_CBFY(:) )
       call FILE_Write_Axis( fid, 'FBFX', GRID_FBFX(:) )
       call FILE_Write_Axis( fid, 'FBFY', GRID_FBFY(:) )
    endif

    if ( (.not. FILE_get_AGGREGATE(fid)) .or. PRC_IsMaster ) then
       call FILE_Write_Axis( fid, 'CXG',   GRID_CXG  (:) )
       call FILE_Write_Axis( fid, 'CYG',   GRID_CYG  (:) )
       call FILE_Write_Axis( fid, 'FXG',   GRID_FXG  (:) )
       call FILE_Write_Axis( fid, 'FYG',   GRID_FYG  (:) )
       call FILE_Write_Axis( fid, 'CDXG',  GRID_CDXG (:) )
       call FILE_Write_Axis( fid, 'CDYG',  GRID_CDYG (:) )
       call FILE_Write_Axis( fid, 'FDXG',  GRID_FDXG (:) )
       call FILE_Write_Axis( fid, 'FDYG',  GRID_FDYG (:) )
       call FILE_Write_Axis( fid, 'CBFXG', GRID_CBFXG(:) )
       call FILE_Write_Axis( fid, 'CBFYG', GRID_CBFYG(:) )
       call FILE_Write_Axis( fid, 'FBFXG', GRID_FBFXG(:) )
       call FILE_Write_Axis( fid, 'FBFYG', GRID_FBFYG(:) )
    end if

    ! associate coordinates
    if ( FILE_get_AGGREGATE(fid) ) then
       call FILE_Write_AssociatedCoordinate( fid, 'lon'   , AXIS_LON  (:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lon_uy', AXIS_LONUY(:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lon_xv', AXIS_LONXV(:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lon_uv', AXIS_LONUV(:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat'   , AXIS_LAT  (:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat_uy', AXIS_LATUY(:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat_xv', AXIS_LATXV(:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat_uv', AXIS_LATUV(:,:), start(2:3) )
       if ( haszcoord ) then
          call FILE_Write_AssociatedCoordinate( fid, 'height'    , AXIS_HGT   (:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'height_wxy', AXIS_HGTWXY(:,:,:), start(1:3) )
       end if
    else
       XSB = ISB - ISB2 + 1
       XEB = IEB - ISB + XSB
       YSB = JSB - JSB2 + 1
       YEB = JEB - JSB + YSB

       call FILE_Write_AssociatedCoordinate( fid, 'lon'   , AXIS_LON  (XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lon_uy', AXIS_LONUY(XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lon_xv', AXIS_LONXV(XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lon_uv', AXIS_LONUV(XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat'   , AXIS_LAT  (XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat_uy', AXIS_LATUY(XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat_xv', AXIS_LATXV(XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'lat_uv', AXIS_LATUV(XSB:XEB,YSB:YEB), start(2:3) )
       if ( haszcoord ) then
          call FILE_Write_AssociatedCoordinate( fid, 'height'    , AXIS_HGT   (:,XSB:XEB,YSB:YEB), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'height_wxy', AXIS_HGTWXY(:,XSB:XEB,YSB:YEB), start(1:3) )
       end if
    end if

    return
  end subroutine FILE_CARTESC_write_axes

  !-----------------------------------------------------------------------------
  !> Define a variable to file
  subroutine FILE_CARTESC_def_var( &
       fid,      &
       varname,  &
       desc,     &
       unit,     &
       dim_type, &
       datatype, &
       vid,      &
       timeintv, &
       nsteps    )
    use scale_file_h, only: &
       FILE_REAL8, &
       FILE_REAL4
    use scale_file, only: &
       FILE_Def_Variable, &
       FILE_Set_Attribute
    use scale_process, only: &
       PRC_abort
    use scale_mapproj, only: &
       MPRJ_get_attributes
    implicit none

    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: dim_type !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer,          intent(out) :: vid      !< variable ID

    real(DP),         intent(in), optional :: timeintv !< time interval [sec]
    integer,          intent(in), optional :: nsteps   !< number of time steps

    integer          :: dtype, ndims, elm_size
    character(len=2) :: dims(3)
    real(DP)         :: time_interval
    character(len=H_SHORT) :: mapping
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if    ( datatype == 'REAL8' ) then
       dtype = FILE_REAL8
       elm_size = 4
    elseif( datatype == 'REAL4' ) then
       dtype = FILE_REAL4
       elm_size = 8
    else
       if    ( RP == 8 ) then
          dtype = FILE_REAL8
          elm_size = 8
       elseif( RP == 4 ) then
          dtype = FILE_REAL4
          elm_size = 4
       else
          write(*,*) 'xxx unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_abort
       endif
    endif

    if    ( dim_type == 'Z' ) then        ! 1D variable
       ndims   = 1
       dims(1) = 'z'
       write_buf_amount = write_buf_amount + KA * elm_size
       mapping = ""
    elseif( dim_type == 'X' ) then
       ndims   = 1
       dims(1) = 'x'
       write_buf_amount = write_buf_amount + IA * elm_size
       mapping = ""
    elseif( dim_type == 'Y' ) then
       ndims   = 1
       dims(1) = 'y'
       write_buf_amount = write_buf_amount + JA * elm_size
       mapping = ""
    elseif( dim_type == 'XY' ) then   ! 2D variable
       ndims   = 2
       dims(1) = 'x'
       dims(2) = 'y'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'UY' ) then
       ndims   = 2
       dims(1) = 'xh'
       dims(2) = 'y'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'XV' ) then
       ndims   = 2
       dims(1) = 'x'
       dims(2) = 'yh'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'UV' ) then
       ndims   = 2
       dims(1) = 'xh'
       dims(2) = 'yh'
       write_buf_amount = write_buf_amount + IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZX' ) then
       ndims   = 2
       dims(1) = 'z'
       dims(2) = 'x'
       write_buf_amount = write_buf_amount + KA * IA * elm_size
       mapping = ""
    elseif( dim_type == 'ZXY' ) then  ! 3D variable
       ndims   = 3
       dims    = (/'z','x','y'/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZHXY' ) then
       ndims   = 3
       dims    = (/'zh','x ','y '/)
       write_buf_amount = write_buf_amount + (KA+1) * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZXHY' ) then
       ndims   = 3
       dims    = (/'z ','xh','y '/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZXYH' ) then
       ndims   = 3
       dims    = (/'z ','x ','yh'/)
       write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'OXY' ) then
       ndims   = 3
       dims    = (/'oz','x ','y '/)
       write_buf_amount = write_buf_amount + OKMAX * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'LXY' ) then
       ndims   = 3
       dims    = (/'lz','x ','y '/)
       write_buf_amount = write_buf_amount + LKMAX * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'UXY' ) then
       ndims   = 3
       dims    = (/'uz','x ','y '/)
       write_buf_amount = write_buf_amount + UKMAX * IA * JA * elm_size
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'XYT' ) then  ! 3D variable with time dimension
       ndims   = 2
       dims(1) = 'x'
       dims(2) = 'y'
       if ( present(nsteps) ) then
          write_buf_amount = write_buf_amount + IA * JA * elm_size * nsteps
       else
          write_buf_amount = write_buf_amount + IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'z','x','y'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + KA * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + KA * IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZHXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'zh','x ','y '/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + (KA+1) * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + (KA+1) * IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZXHYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'z ','xh','y '/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + KA * (IA+1) * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + KA * (IA+1) * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'ZXYHT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'z ','x ','yh'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + KA * IA * (JA+1) * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + KA * IA * (JA+1) * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'OXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'oz','x','y'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + OKMAX * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + OKMAX * IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'OHXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'ozh','x','y'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + (OKMAX+1) * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + (OKMAX+1) * IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'LXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'lz','x','y'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + LKMAX * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + LKMAX * IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    elseif( dim_type == 'LHXYT' ) then ! 4D variable
       ndims   = 3
       dims    = (/'lzh','x','y'/)
       if ( present(nsteps) ) then
         write_buf_amount = write_buf_amount + (LKMAX+1) * IA * JA * elm_size * nsteps
       else
         write_buf_amount = write_buf_amount + (LKMAX+1) * IA * JA * elm_size
       endif
       call MPRJ_get_attributes( mapping )
    else
       write(*,*) 'xxx [FILE_CARTESC_def_var] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if ( present(timeintv) ) then  ! 3D/4D variable with time dimension
      time_interval = timeintv
      call FILE_Def_Variable( fid, varname, desc, unit, ndims, dims, dtype, & ! [IN]
                              vid,                                          & ! [OUT]
                              timeintv=time_interval                        ) ! [IN]
    else
      call FILE_Def_Variable( fid, varname, desc, unit, ndims, dims, dtype, & ! [IN]
                              vid                                           ) ! [OUT]
    endif

    if ( mapping /= "" ) call FILE_Set_Attribute( fid, varname, "grid_mapping", mapping )

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_def_var

  !-----------------------------------------------------------------------------
  !> Write 1D data to file
  subroutine FILE_CARTESC_write_var_1D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       dim_type  )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Write
    use scale_process, only: &
       PRC_myrank,  &
       PRC_abort
    use scale_rm_process, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    implicit none

    integer,          intent(in) :: fid      !< file ID
    integer,          intent(in) :: vid      !< netCDF variable ID
    real(RP),         intent(in) :: var(:)   !< value of the variable
    character(len=*), intent(in) :: varname  !< name of the variable
    character(len=*), intent(in) :: dim_type !< axis type (Z/X/Y)

    integer :: dim1_S, dim1_E
    integer :: rankidx(2)
    integer :: start(1)         ! used only when AGGREGATE is .true.
    logical :: exec
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if    ( dim_type == 'Z' ) then
       dim1_S   = KS
       dim1_E   = KE
       start(1) = 1
       if( FILE_get_AGGREGATE(fid) .AND. PRC_myrank > 0 ) then
          exec = .false.  ! only rank 0 writes
       else
          exec = .true.
       end if
    elseif( dim_type == 'X' ) then
       if ( FILE_get_AGGREGATE(fid) ) then
          exec = ( rankidx(2) == 0 ) ! only south most row processes write
          dim1_S   = ISB2
          dim1_E   = IEB2
       else
          exec = .true.
          dim1_S   = ISB
          dim1_E   = IEB
       end if
       start(1) = ISGA
    elseif( dim_type == 'Y' ) then
       if ( FILE_get_AGGREGATE(fid) ) then
          exec = ( rankidx(1) == 0 ) ! only west most column processes write
          dim1_S   = JSB2
          dim1_E   = JEB2
       else
          exec = .true.
          dim1_S   = JSB
          dim1_E   = JEB
       end if
       start(1) = JSGA
    else
       write(*,*) 'xxx [FILE_CARTESC_write_var_1D] unsupported dimenstion type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if( exec ) call FILE_Write( vid, var(dim1_S:dim1_E),    & ! [IN]
                                NOWSEC, NOWSEC, start=start ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_write_var_1D

  !-----------------------------------------------------------------------------
  !> Write 2D data to file
  subroutine FILE_CARTESC_write_var_2D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       dim_type, &
       fill_halo )
    use scale_file_h, only: &
       RMISS => FILE_RMISS
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Write
    use scale_process, only: &
       PRC_myrank,     &
       PRC_abort
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
    character(len=*), intent(in)  :: dim_type !< axis type (Z/X/Y)
    logical,          intent(in), optional :: fill_halo !< switch whether include halo data or not (default=false)

    real(RP)              :: varhalo( size(var(:,1)), size(var(1,:)) )

    integer               :: dim1_S, dim1_E
    integer               :: dim2_S, dim2_E

    integer :: i, j
    logical :: fill_halo_
    integer :: rankidx(2)
    integer :: start(2)         ! used only when AGGREGATE is .true.
    logical :: exec
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = ISGA
    start(2) = JSGA

    fill_halo_ = .false.
    if( present(fill_halo) ) fill_halo_ = fill_halo

    if (      dim_type == 'XY' &
         .OR. dim_type == 'UY' &
         .OR. dim_type == 'XV' &
         .OR. dim_type == 'UV' ) then
       if ( FILE_get_AGGREGATE(fid) ) then
          dim1_S   = ISB2
          dim1_E   = IEB2
          dim2_S   = JSB2
          dim2_E   = JEB2
       else
          dim1_S   = ISB
          dim1_E   = IEB
          dim2_S   = JSB
          dim2_E   = JEB
       endif
       exec = .true.
    elseif( dim_type == 'ZX' ) then
       dim1_S   = KS
       dim1_E   = KE
       start(2) = start(1)
       start(1) = 1
       if ( FILE_get_AGGREGATE(fid) ) then
          exec = ( rankidx(2) == 0 )  ! only south most row processes write
          dim2_S = ISB2
          dim2_E = IEB2
       else
          exec = .true.
          dim2_S = ISB
          dim2_E = IEB
       endif
    else
       write(*,*) 'xxx [FILE_CARTESC_write_var_2D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if ( exec ) then
       if ( fill_halo_ ) then ! fill halo cells with RMISS
          varhalo(:,:) = var(:,:)

          ! W halo
          do j = 1, JA
          do i = 1, IS-1
             varhalo(i,j) = RMISS
          enddo
          enddo
          ! E halo
          do j = 1, JA
          do i = IE+1, IA
             varhalo(i,j) = RMISS
          enddo
          enddo
          ! S halo
          do j = 1, JS-1
          do i = 1, IA
             varhalo(i,j) = RMISS
          enddo
          enddo
          ! N halo
          do j = JE+1, JA
          do i = 1, IA
             varhalo(i,j) = RMISS
          enddo
          enddo

          call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), &
                           NOWSEC, NOWSEC, start                      ) ! [IN]
       else
          call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E), &
                           NOWSEC, NOWSEC, start                  ) ! [IN]
       endif

    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_write_var_2D

  !-----------------------------------------------------------------------------
  !> Write 3D data to file
  subroutine FILE_CARTESC_write_var_3D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       dim_type, &
       fill_halo )
    use scale_file_h, only: &
       RMISS => FILE_RMISS
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Write
    use scale_process, only: &
       PRC_myrank,  &
       PRC_abort
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
    character(len=*), intent(in) :: dim_type   !< dimension type (Z/X/Y)

    logical,          intent(in), optional :: fill_halo !< include halo data?

    real(RP) :: varhalo( size(var(:,1,1)), size(var(1,:,1)), size(var(1,1,:)) )

    integer  :: dim1_S, dim1_E, dim1_max
    integer  :: dim2_S, dim2_E
    integer  :: dim3_S, dim3_E

    integer  :: i, j, k
    logical  :: fill_halo_
    integer  :: rankidx(2)
    integer  :: start(3) ! used only when AGGREGATE is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    fill_halo_ = .false.
    if( present(fill_halo) ) fill_halo_ = fill_halo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA

    if (      dim_type == 'ZXY'  &
         .OR. dim_type == 'ZXHY' &
         .OR. dim_type == 'ZXYH' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif ( dim_type == 'ZHXY' ) then
       dim1_max = KMAX+1
       dim1_S   = KS-1
       dim1_E   = KE
    elseif( dim_type == 'OXY' ) then
       dim1_max = OKMAX
       dim1_S   = OKS
       dim1_E   = OKE
    elseif( dim_type == 'LXY' ) then
       dim1_max = LKMAX
       dim1_S   = LKS
       dim1_E   = LKE
    elseif( dim_type == 'UXY' ) then
       dim1_max = UKMAX
       dim1_S   = UKS
       dim1_E   = UKE
    else
       write(*,*) 'xxx [FILE_CARTESC_write_var_3D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if ( FILE_get_AGGREGATE(fid) ) then
       dim2_S   = ISB2
       dim2_E   = IEB2
       dim3_S   = JSB2
       dim3_E   = JEB2
    else
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    endif

    if ( fill_halo_ ) then
       varhalo(:,:,:) = var(:,:,:)

       ! W halo
       do k = 1, dim1_max
       do j = 1, JA
       do i = 1, IS-1
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo
       ! E halo
       do k = 1, dim1_max
       do j = 1, JA
       do i = IE+1, IA
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo
       ! S halo
       do k = 1, dim1_max
       do j = 1, JS-1
       do i = 1, IA
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo
       ! N halo
       do k = 1, dim1_max
       do j = JE+1, JA
       do i = 1, IA
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo

       call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                        NOWSEC, NOWSEC, start                                    ) ! [IN]
    else
       call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                        NOWSEC, NOWSEC, start                                ) ! [IN]
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_write_var_3D

  !-----------------------------------------------------------------------------
  !> Write 3D data with time dimension to file
  subroutine FILE_CARTESC_write_var_3D_t( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       dim_type, &
       timeintv, &
       timetarg, &
       timeofs,  &
       fill_halo )
    use scale_file_h, only: &
       RMISS => FILE_RMISS
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Write
    use scale_process, only: &
       PRC_myrank,     &
       PRC_abort
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    integer,          intent(in) :: fid        !< file ID
    integer,          intent(in) :: vid        !< netCDF variable ID
    real(RP),         intent(in) :: var(:,:,:) !< value of the variable
    character(len=*), intent(in) :: varname    !< name of the variable
    character(len=*), intent(in) :: dim_type   !< dimension type (X/Y/Time)
    real(DP),         intent(in) :: timeintv   !< time interval [sec]

    integer,          intent(in), optional :: timetarg  !< target timestep (optional)
    real(DP),         intent(in), optional :: timeofs   !< offset time     (optional)
    logical,          intent(in), optional :: fill_halo !< include halo data?

    real(RP) :: varhalo( size(var(:,1,1)), size(var(1,:,1)) )

    integer  :: dim1_S, dim1_E
    integer  :: dim2_S, dim2_E

    real(DP) :: time_interval, nowtime

    integer  :: step
    integer  :: i, j, n
    logical  :: fill_halo_
    real(DP) :: timeofs_
    integer  :: rankidx(2)
    integer  :: start(2) ! used only when AGGREGATE is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    fill_halo_ = .false.
    if( present(fill_halo) ) fill_halo_ = fill_halo

    timeofs_ = 0.0_DP
    if( present(timeofs) ) timeofs_ = timeofs

    time_interval = timeintv
    step = size(var(ISB,JSB,:))

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    if ( dim_type == 'XYT' ) then
       if ( FILE_get_AGGREGATE(fid) ) then
          dim1_S   = ISB2
          dim1_E   = IEB2
          dim2_S   = JSB2
          dim2_E   = JEB2
       else
          dim1_S   = ISB
          dim1_E   = IEB
          dim2_S   = JSB
          dim2_E   = JEB
       end if
    else
       write(*,*) 'xxx [FILE_CARTESC_write_var_3D_t] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    start(1) = ISGA
    start(2) = JSGA
    ! start(3) time dimension will be set in file_write_data()

    if ( present(timetarg) ) then
       nowtime = timeofs_ + (timetarg-1) * time_interval

       if ( fill_halo_ ) then
          varhalo(:,:) = var(:,:,timetarg)

          ! W halo
          do j = 1, JA
          do i = 1, IS-1
             varhalo(i,j) = RMISS
          enddo
          enddo
          ! E halo
          do j = 1, JA
          do i = IE+1, IA
             varhalo(i,j) = RMISS
          enddo
          enddo
          ! S halo
          do j = 1, JS-1
          do i = 1, IA
             varhalo(i,j) = RMISS
          enddo
          enddo
          ! N halo
          do j = JE+1, JA
          do i = 1, IA
             varhalo(i,j) = RMISS
          enddo
          enddo

          call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), &
                           nowtime, nowtime, start                    ) ! [IN]
       else
          call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,timetarg), &
                           nowtime, nowtime, start                         ) ! [IN]
       endif
    else
       nowtime = timeofs_
       do n = 1, step
          if ( fill_halo_ ) then
             varhalo(:,:) = var(:,:,n)

             ! W halo
             do j = 1, JA
             do i = 1, IS-1
                varhalo(i,j) = RMISS
             enddo
             enddo
             ! E halo
             do j = 1, JA
             do i = IE+1, IA
                varhalo(i,j) = RMISS
             enddo
             enddo
             ! S halo
             do j = 1, JS-1
             do i = 1, IA
                varhalo(i,j) = RMISS
             enddo
             enddo
             ! N halo
             do j = JE+1, JA
             do i = 1, IA
                varhalo(i,j) = RMISS
             enddo
             enddo

             call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E), &
                              nowtime, nowtime, start                    ) ! [IN]
          else
             call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,n), &
                              nowtime, nowtime, start                  ) ! [IN]
          endif
          nowtime = nowtime + time_interval
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_write_var_3D_t

  !-----------------------------------------------------------------------------
  !> Write 4D data to file
  subroutine FILE_CARTESC_write_var_4D( &
       fid,      &
       vid,      &
       var,      &
       varname,  &
       dim_type, &
       timeintv, &
       timetarg, &
       timeofs,  &
       fill_halo )
    use scale_file_h, only: &
       RMISS => FILE_RMISS
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_Write, &
       FILE_Flush
    use scale_process, only: &
       PRC_myrank,     &
       PRC_abort
    use scale_rm_process, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none

    integer,          intent(in) :: fid          !< file ID
    integer,          intent(in) :: vid          !< netCDF variable ID
    real(RP),         intent(in) :: var(:,:,:,:) !< value of the variable
    character(len=*), intent(in) :: varname      !< name        of the variable
    character(len=*), intent(in) :: dim_type     !< dimension type (Z/X/Y/Time)
    real(DP),         intent(in) :: timeintv     !< time interval [sec]

    integer,          intent(in), optional :: timetarg  !< target timestep (optional)
    real(DP),         intent(in), optional :: timeofs   !< offset time     (optional)
    logical,          intent(in), optional :: fill_halo !< include halo data?

    real(RP) :: varhalo( size(var(:,1,1,1)), size(var(1,:,1,1)), size(var(1,1,:,1)) )

    integer  :: dim1_S, dim1_E, dim1_max
    integer  :: dim2_S, dim2_E
    integer  :: dim3_S, dim3_E

    real(DP) :: time_interval, nowtime

    integer  :: step
    integer  :: i, j, k, n
    logical  :: fill_halo_
    real(DP) :: timeofs_
    integer  :: rankidx(2)
    integer  :: start(3) ! used only when AGGREGATE is .true.
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    fill_halo_ = .false.
    if( present(fill_halo) ) fill_halo_ = fill_halo

    timeofs_ = 0.0_DP
    if( present(timeofs) ) timeofs_ = timeofs

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA
    ! start(4) time dimension will be set in file_write_data()

    time_interval = timeintv
    step = size(var,4)

    if ( FILE_get_AGGREGATE(fid) ) then
       dim2_S   = ISB2
       dim2_E   = IEB2
       dim3_S   = JSB2
       dim3_E   = JEB2
    else
       dim2_S   = ISB
       dim2_E   = IEB
       dim3_S   = JSB
       dim3_E   = JEB
    endif

    if (      dim_type == 'ZXYT'  &
         .OR. dim_type == 'ZXHYT' &
         .OR. dim_type == 'ZXYHT' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif ( dim_type == 'ZHXYT' ) then
       dim1_max = KMAX+1
       dim1_S   = KS-1
       dim1_E   = KE
    elseif( dim_type == 'OXYT' ) then
       dim1_max = OKMAX
       dim1_S   = OKS
       dim1_E   = OKE
    elseif( dim_type == 'OHXYT' ) then
       dim1_max = OKMAX+1
       dim1_S   = OKS-1
       dim1_E   = OKE
    elseif( dim_type == 'LXYT' ) then
       dim1_max = LKMAX
       dim1_S   = LKS
       dim1_E   = LKE
    elseif( dim_type == 'LHXYT' ) then
       dim1_max = LKMAX+1
       dim1_S   = LKS-1
       dim1_E   = LKE
    elseif( dim_type == 'UXYT' ) then
       dim1_max = UKMAX
       dim1_S   = UKS
       dim1_E   = UKE
    elseif( dim_type == 'UHXYT' ) then
       dim1_max = UKMAX+1
       dim1_S   = UKS-1
       dim1_E   = UKE
    else
       write(*,*) 'xxx [FILE_CARTESC_write_var_4D] unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if ( present(timetarg) ) then
       nowtime = timeofs_ + (timetarg-1) * time_interval

       if ( fill_halo_ ) then
          varhalo(:,:,:) = var(:,:,:,timetarg)

          ! W halo
          do k = 1, dim1_max
          do j = 1, JA
          do i = 1, IS-1
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo
          ! E halo
          do k = 1, dim1_max
          do j = 1, JA
          do i = IE+1, IA
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo
          ! S halo
          do k = 1, dim1_max
          do j = 1, JS-1
          do i = 1, IA
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo
          ! N halo
          do k = 1, dim1_max
          do j = JE+1, JA
          do i = 1, IA
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo

          call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                           nowtime, nowtime, start                                  ) ! [IN]
       else
          call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,timetarg), &
                           nowtime, nowtime, start                                       ) ! [IN]
       endif
    else
       nowtime = timeofs_
       do n = 1, step
          if ( fill_halo_ ) then
             varhalo(:,:,:) = var(:,:,:,n)

             ! W halo
             do k = 1, dim1_max
             do j = 1, JA
             do i = 1, IS-1
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo
             ! E halo
             do k = 1, dim1_max
             do j = 1, JA
             do i = IE+1, IA
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo
             ! S halo
             do k = 1, dim1_max
             do j = 1, JS-1
             do i = 1, IA
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo
             ! N halo
             do k = 1, dim1_max
             do j = JE+1, JA
             do i = 1, IA
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo

             call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                              nowtime, nowtime, start                                  ) ! [IN]
          else
             call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,n), &
                              nowtime, nowtime, start                                ) ! [IN]
          endif
          call FILE_Flush( fid )
          nowtime = nowtime + time_interval
       enddo
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_write_var_4D

  !-----------------------------------------------------------------------------
  subroutine closeall
    implicit none

    integer :: fid
    !---------------------------------------------------------------------------

    do fid = 0, FILE_FILE_MAX-1
       call FILE_CARTESC_close( fid )
    enddo

    return
  end subroutine closeall

  !-----------------------------------------------------------------------------
  subroutine check_1d( &
       expected, buffer, &
       name              )
    use scale_process, only: &
       PRC_abort
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP),         intent(in) :: expected(:)
    real(RP),         intent(in) :: buffer(:)
    character(len=*), intent(in) :: name

    real(RP) :: check
    integer  :: nmax, n

    intrinsic :: size
    !---------------------------------------------------------------------------

    nmax = size(expected)
    if ( size(buffer) /= nmax ) then
       write(*,*) 'xxx size of coordinate ('//trim(name)//') is different:', nmax, size(buffer)
       call PRC_abort
    endif

    do n=1, nmax
       if ( abs(expected(n)) > EPS ) then
          check = abs(buffer(n)-expected(n)) / abs(buffer(n)+expected(n)) * 2.0_RP
       else
          check = abs(buffer(n)-expected(n))
       endif

       if ( check > FILE_CARTESC_datacheck_criteria ) then
          write(*,*) 'xxx value of coordinate ('//trim(name)//') at ', n, ' is different:', &
                     expected(n), buffer(n), check
          call PRC_abort
       endif
    enddo

    return
  end subroutine check_1d

  !-----------------------------------------------------------------------------
  subroutine check_2d( &
       expected, buffer, &
       name              )
    use scale_process, only: &
       PRC_abort
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP),         intent(in) :: expected(:,:)
    real(RP),         intent(in) :: buffer(:,:)
    character(len=*), intent(in) :: name

    real(RP) :: check
    integer  :: imax, jmax
    integer  :: i, j

    intrinsic :: size
    !---------------------------------------------------------------------------

    imax = size(expected,1)
    jmax = size(expected,2)
    if ( size(buffer,1) /= imax ) then
       write(*,*) 'xxx the first size of coordinate ('//trim(name)//') is different:', imax, size(buffer,1)
       call PRC_abort
    endif
    if ( size(buffer,2) /= jmax ) then
       write(*,*) 'xxx the second size of coordinate ('//trim(name)//') is different:', jmax, size(buffer,2)
       call PRC_abort
    endif

    do j=1, jmax
    do i=1, imax
       if ( abs(expected(i,j)) > EPS ) then
          check = abs(buffer(i,j)-expected(i,j)) / abs(buffer(i,j)+expected(i,j)) * 2.0_RP
       else
          check = abs(buffer(i,j)-expected(i,j))
       endif

       if ( check > FILE_CARTESC_datacheck_criteria ) then
          write(*,*) 'xxx value of coordinate ('//trim(name)//') at (', i, ',', j, ') is different:', &
                     expected(i,j), buffer(i,j), check
          call PRC_abort
       endif
    enddo
    enddo

    return
  end subroutine check_2d

  !-----------------------------------------------------------------------------
  subroutine check_3d( &
       expected, buffer, &
       name,             &
       transpose         )
    use scale_process, only: &
       PRC_abort
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP),         intent(in) :: expected(:,:,:)
    real(RP),         intent(in) :: buffer(:,:,:)
    character(len=*), intent(in) :: name
    logical,          intent(in) :: transpose

    real(RP) :: check
    integer  :: imax, jmax, kmax
    integer  :: i, j, k

    intrinsic :: size
    !---------------------------------------------------------------------------

    if ( transpose ) then
       kmax = size(expected,3)
       imax = size(expected,1)
       jmax = size(expected,2)
    else
       kmax = size(expected,1)
       imax = size(expected,2)
       jmax = size(expected,3)
    endif
    if ( size(buffer,1) /= kmax ) then
       write(*,*) 'xxx the first size of coordinate ('//trim(name)//') is different:', kmax, size(buffer,1)
       call PRC_abort
    endif
    if ( size(buffer,2) /= imax ) then
       write(*,*) 'xxx the second size of coordinate ('//trim(name)//') is different:', imax, size(buffer,2)
       call PRC_abort
    endif
    if ( size(buffer,3) /= jmax ) then
       write(*,*) 'xxx the third size of coordinate ('//trim(name)//') is different:', jmax, size(buffer,3)
       call PRC_abort
    endif

    if ( transpose ) then
       ! buffer(i,j,k), expected(k,i,j)
       do k=1, kmax
       do j=1, jmax
       do i=1, imax
          if ( abs(expected(k,i,j)) > EPS ) then
             check = abs(buffer(i,j,k)-expected(k,i,j)) / abs(buffer(i,j,k)+expected(k,i,j)) * 2.0_RP
          else
             check = abs(buffer(i,j,k)-expected(k,i,j))
          endif

          if ( check > FILE_CARTESC_datacheck_criteria ) then
             write(*,*) 'xxx value of coordinate ('//trim(name)//') at ', i, ',', j, ',', k, ' is different:', &
                        expected(k,i,j), buffer(i,j,k), check
             call PRC_abort
          endif
       enddo
       enddo
       enddo
    else
       do j=1, jmax
       do i=1, imax
       do k=1, kmax
          if ( abs(expected(k,i,j)) > EPS ) then
             check = abs(buffer(k,i,j)-expected(k,i,j)) / abs(buffer(k,i,j)+expected(k,i,j)) * 2.0_RP
          else
             check = abs(buffer(k,i,j)-expected(k,i,j))
          endif

          if ( check > FILE_CARTESC_datacheck_criteria ) then
             write(*,*) 'xxx value of coordinate ('//trim(name)//') at ', k, ',', i, ',', j, ' is different:', &
                        expected(k,i,j), buffer(k,i,j), check
             call PRC_abort
          endif
       enddo
       enddo
       enddo
    endif

    return
  end subroutine check_3d

end module scale_file_cartesC
