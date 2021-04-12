!-------------------------------------------------------------------------------
!> module file / cartesianC
!!
!! @par Description
!!          file I/O handling for the cartesianC grid
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_file_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_land_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  use scale_urban_grid_cartesC_index
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
  public :: FILE_CARTESC_finalize
  public :: FILE_CARTESC_set_coordinates_atmos
  public :: FILE_CARTESC_set_areavol_atmos
  public :: FILE_CARTESC_set_coordinates_ocean
  public :: FILE_CARTESC_set_coordinates_land
  public :: FILE_CARTESC_set_coordinates_urban
  public :: FILE_CARTESC_check_coordinates

  public :: FILE_CARTESC_get_size
  public :: FILE_CARTESC_create
  public :: FILE_CARTESC_open
  public :: FILE_CARTESC_put_globalAttributes
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

  interface FILE_CARTESC_get_size
     module procedure FILE_CARTESC_get_size_id
     module procedure FILE_CARTESC_get_size_name
  end interface FILE_CARTESC_get_size

  interface FILE_CARTESC_read
     module procedure FILE_CARTESC_read_1D
     module procedure FILE_CARTESC_read_2D
     module procedure FILE_CARTESC_read_3D
     module procedure FILE_CARTESC_read_4D
     module procedure FILE_CARTESC_read_var_1D
     module procedure FILE_CARTESC_read_var_2D
     module procedure FILE_CARTESC_read_var_3D
     module procedure FILE_CARTESC_read_var_4D
     module procedure FILE_CARTESC_read_auto_2D
     module procedure FILE_CARTESC_read_auto_3D
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
    integer :: size_global (1)
    integer :: start_global(1)
    integer :: halo_global (2)
    integer :: halo_local  (2)
    logical :: periodic
  end type axisattinfo

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
  real(RP), private :: FILE_CARTESC_datacheck_criteria

  type dims
     character(len=H_SHORT) :: name
     integer                :: ndims
     character(len=H_SHORT) :: dims(3)
     integer                :: size
     logical                :: mapping
     character(len=H_SHORT) :: area
     character(len=H_SHORT) :: area_x
     character(len=H_SHORT) :: area_y
     character(len=H_SHORT) :: volume
     character(len=H_SHORT) :: location
     character(len=H_SHORT) :: grid
  end type dims
  integer, parameter :: FILE_CARTESC_ndims = 44
  type(dims) :: FILE_CARTESC_dims(FILE_CARTESC_ndims)

  type(axisattinfo) :: FILE_CARTESC_AXIS_info(4) ! x, xh, y, yh


  real(RP), private, allocatable :: AXIS_HGT   (:,:,:)
  real(RP), private, allocatable :: AXIS_HGTWXY(:,:,:)

  real(RP), private, allocatable :: AXIS_LON  (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONUY(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONXV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LONUV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LAT  (:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATUY(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATXV(:,:)   ! [deg]
  real(RP), private, allocatable :: AXIS_LATUV(:,:)   ! [deg]

  real(RP), private, allocatable :: AXIS_TOPO  (:,:)
  real(RP), private, allocatable :: AXIS_LSMASK(:,:)

  real(RP), private, allocatable :: AXIS_AREA     (:,:)
  real(RP), private, allocatable :: AXIS_AREAZUY_X(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAZXV_Y(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAWUY_X(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAWXV_Y(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAUY   (:,:)
  real(RP), private, allocatable :: AXIS_AREAZXY_X(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAZUV_Y(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAXV   (:,:)
  real(RP), private, allocatable :: AXIS_AREAZUV_X(:,:,:)
  real(RP), private, allocatable :: AXIS_AREAZXY_Y(:,:,:)

  real(RP), private, allocatable :: AXIS_VOL   (:,:,:)
  real(RP), private, allocatable :: AXIS_VOLWXY(:,:,:)
  real(RP), private, allocatable :: AXIS_VOLZUY(:,:,:)
  real(RP), private, allocatable :: AXIS_VOLZXV(:,:,:)

  real(RP), private, allocatable :: AXIS_VOLO(:,:,:)
  real(RP), private, allocatable :: AXIS_VOLL(:,:,:)
  real(RP), private, allocatable :: AXIS_VOLU(:,:,:)

  logical,    private :: File_axes_written(0:FILE_FILE_MAX-1)          ! whether axes have been written
  !                                                                    ! fid starts from zero so index should start from zero
  logical,    private :: File_haszcoord   (0:FILE_FILE_MAX-1)          ! z-coordinates exist?
  integer(8), private :: write_buf_amount (0:FILE_FILE_MAX-1)          ! sum of write buffer amounts

  ! global star and count
  integer,  private, target :: startXY   (3), countXY   (3)
  integer,  private, target :: startZX   (2), countZX   (2)
  integer,  private, target :: startZXY  (4), countZXY  (4)
  integer,  private, target :: startZHXY (4), countZHXY (4)
  integer,  private, target :: startOCEAN(4), countOCEAN(4)
  integer,  private, target :: startLAND (4), countLAND (4)
  integer,  private, target :: startURBAN(4), countURBAN(4)
  ! local start and end
  integer,  private :: ISB2, IEB2, JSB2, JEB2 !> for FILE_AGGREGATE

  ! MPI element datatype for restart variables
  integer,  private :: etype

  ! MPI derived datatypes
  integer,  private :: centerTypeXY
  integer,  private :: centerTypeZX
  integer,  private :: centerTypeZXY
  integer,  private :: centerTypeZHXY
  integer,  private :: centerTypeOCEAN
  integer,  private :: centerTypeLAND
  integer,  private :: centerTypeURBAN

  logical,  private :: set_coordinates = .false.

  logical,  private :: prof = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FILE_CARTESC_setup
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_file, only: &
       FILE_setup
    implicit none

    namelist / PARAM_FILE_CARTESC / &
       FILE_CARTESC_datacheck_criteria

    integer :: IM, JM
    integer :: ierr
    !---------------------------------------------------------------------------

    call FILE_setup( PRC_myrank )


    LOG_NEWLINE
    LOG_INFO("FILE_CARTESC_setup",*) 'Setup'

    FILE_CARTESC_datacheck_criteria = 0.1_RP**(RP)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_FILE_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("FILE_CARTESC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("FILE_CARTESC_setup",*) 'Not appropriate names in namelist PARAM_FILE_CARTESC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_FILE_CARTESC)

    LOG_NEWLINE
    LOG_INFO("FILE_CARTESC_setup",*) 'NetCDF header information '
    LOG_INFO_CONT(*) 'Data source : ', trim(H_SOURCE)
    LOG_INFO_CONT(*) 'Institute   : ', trim(H_INSTITUTE)

    LOG_NEWLINE
    LOG_INFO("FILE_CARTESC_setup",*) 'Data consistency criteria : ', &
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

    IM = IEB2 - ISB2 + 1
    JM = JEB2 - JSB2 + 1
    allocate( AXIS_HGT   (KMAX  ,IM,JM) )
    allocate( AXIS_HGTWXY(KMAX+1,IM,JM) )

    allocate( AXIS_LON  (IM,JM) )
    allocate( AXIS_LONUY(IM,JM) )
    allocate( AXIS_LONXV(IM,JM) )
    allocate( AXIS_LONUV(IM,JM) )
    allocate( AXIS_LAT  (IM,JM) )
    allocate( AXIS_LATUY(IM,JM) )
    allocate( AXIS_LATXV(IM,JM) )
    allocate( AXIS_LATUV(IM,JM) )

    allocate( AXIS_TOPO  (IM,JM) )
    allocate( AXIS_LSMASK(IM,JM) )

    allocate( AXIS_AREA     (       IM,JM) )
    allocate( AXIS_AREAZUY_X(KMAX,  IM,JM) )
    allocate( AXIS_AREAZXV_Y(KMAX,  IM,JM) )
    allocate( AXIS_AREAWUY_X(KMAX+1,IM,JM) )
    allocate( AXIS_AREAWXV_Y(KMAX+1,IM,JM) )
    allocate( AXIS_AREAUY   (       IM,JM) )
    allocate( AXIS_AREAZXY_X(KMAX,  IM,JM) )
    allocate( AXIS_AREAZUV_Y(KMAX,  IM,JM) )
    allocate( AXIS_AREAXV   (       IM,JM) )
    allocate( AXIS_AREAZUV_X(KMAX,  IM,JM) )
    allocate( AXIS_AREAZXY_Y(KMAX,  IM,JM) )

    allocate( AXIS_VOL   (KMAX  ,IM,JM) )
    allocate( AXIS_VOLWXY(KMAX+1,IM,JM) )
    allocate( AXIS_VOLZUY(KMAX  ,IM,JM) )
    allocate( AXIS_VOLZXV(KMAX  ,IM,JM) )

    allocate( AXIS_VOLO(OKMAX,IM,JM) )
    allocate( AXIS_VOLL(LKMAX,IM,JM) )
    allocate( AXIS_VOLU(UKMAX,IM,JM) )

    call Construct_Derived_Datatype

    write_buf_amount(:) = 0

    return
  end subroutine FILE_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> deallocate buffers
  subroutine FILE_CARTESC_finalize
    implicit none
    !---------------------------------------------------------------------------

    deallocate( AXIS_HGT    )
    deallocate( AXIS_HGTWXY )

    deallocate( AXIS_LON   )
    deallocate( AXIS_LONUY )
    deallocate( AXIS_LONXV )
    deallocate( AXIS_LONUV )
    deallocate( AXIS_LAT   )
    deallocate( AXIS_LATUY )
    deallocate( AXIS_LATXV )
    deallocate( AXIS_LATUV )

    deallocate( AXIS_TOPO   )
    deallocate( AXIS_LSMASK )

    deallocate( AXIS_AREA      )
    deallocate( AXIS_AREAZUY_X )
    deallocate( AXIS_AREAZXV_Y )
    deallocate( AXIS_AREAWUY_X )
    deallocate( AXIS_AREAWXV_Y )
    deallocate( AXIS_AREAUY    )
    deallocate( AXIS_AREAZXY_X )
    deallocate( AXIS_AREAZUV_Y )
    deallocate( AXIS_AREAXV    )
    deallocate( AXIS_AREAZUV_X )
    deallocate( AXIS_AREAZXY_Y )

    deallocate( AXIS_VOL    )
    deallocate( AXIS_VOLWXY )
    deallocate( AXIS_VOLZUY )
    deallocate( AXIS_VOLZXV )

    deallocate( AXIS_VOLO )
    deallocate( AXIS_VOLL )
    deallocate( AXIS_VOLU )

    call Free_Derived_Datatype

    call closeall

    set_coordinates = .false.

    return
  end subroutine FILE_CARTESC_finalize

  !-----------------------------------------------------------------------------
  !> Get dimension information from file
  !! This subroutine can be called without setup
  !<
  !-----------------------------------------------------------------------------
  subroutine FILE_CARTESC_get_size_name( &
       basename,                  &
       KMAX, OKMAX, LKMAX, UKMAX, &
       IMAXG, JMAXG,              &
       KHALO, IHALO, JHALO,       &
       aggregate                  )
    use scale_file, only: &
       FILE_open
    character(len=*), intent(in) :: basename

    integer, intent(out) :: KMAX, OKMAX, LKMAX, UKMAX
    integer, intent(out) :: IMAXG, JMAXG
    integer, intent(out) :: KHALO, IHALO, JHALO

    logical, intent(in), optional :: aggregate

    integer :: fid

    call FILE_open( basename,           & ! (in)
                    fid,                & ! (out)
                    aggregate=aggregate ) ! (in)

    call FILE_CARTESC_get_size_id( fid,                       & ! (in)
                                   KMAX, OKMAX, LKMAX, UKMAX, & ! (out)
                                   IMAXG, JMAXG,              & ! (out)
                                   KHALO, IHALO, JHALO        ) ! (out)

    return
  end subroutine FILE_CARTESC_get_size_name
  subroutine FILE_CARTESC_get_size_id( &
       fid,                       &
       KMAX, OKMAX, LKMAX, UKMAX, &
       IMAXG, JMAXG,              &
       KHALO, IHALO, JHALO        )
    use scale_file, only: &
       FILE_get_attribute

    integer, intent(in) :: fid

    integer, intent(out) :: KMAX, OKMAX, LKMAX, UKMAX
    integer, intent(out) :: IMAXG, JMAXG
    integer, intent(out) :: KHALO, IHALO, JHALO

    integer :: buf(1)
    logical :: existed

    call FILE_Get_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_kmax",  buf(:)  )
    KMAX = buf(1)
    call FILE_Get_Attribute( fid, "global", "scale_ocean_grid_cartesC_index_kmax",  buf(:), existed=existed  )
    if ( existed ) then
       OKMAX = buf(1)
    else
       OKMAX = -1
    end if
    call FILE_Get_Attribute( fid, "global", "scale_land_grid_cartesC_index_kmax",  buf(:), existed=existed  )
    if ( existed ) then
       LKMAX = buf(1)
    else
       LKMAX = -1
    end if
    call FILE_Get_Attribute( fid, "global", "scale_urban_grid_cartesC_index_kmax",  buf(:), existed=existed  )
    if ( existed ) then
       UKMAX = buf(1)
    else
       UKMAX = -1
    end if

    call FILE_Get_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_imaxg", buf(:)  )
    IMAXG = buf(1)
    call FILE_Get_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_jmaxg", buf(:)  )
    JMAXG = buf(1)

    call FILE_Get_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_khalo", buf(:)  )
    KHALO = buf(1)
    call FILE_Get_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_ihalo", buf(:)  )
    IHALO = buf(1)
    call FILE_Get_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_jhalo", buf(:)  )
    JHALO = buf(1)

    return
  end subroutine FILE_CARTESC_get_size_id

  !-----------------------------------------------------------------------------
  !> set latlon and z for atmosphere
  subroutine FILE_CARTESC_set_coordinates_atmos( &
       CZ, FZ,                       &
       LON, LONUY, LONXV, LONUV,     &
       LAT, LATUY, LATXV, LATUV,     &
       TOPO, LSMASK                 )
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    real(RP), intent(in) :: CZ(  KA,IA,JA)
    real(RP), intent(in) :: FZ(0:KA,IA,JA)
    real(RP), intent(in) :: LON  (  IA,  JA)
    real(RP), intent(in) :: LONUY(0:IA,  JA)
    real(RP), intent(in) :: LONXV(  IA,0:JA)
    real(RP), intent(in) :: LONUV(0:IA,0:JA)
    real(RP), intent(in) :: LAT  (  IA,  JA)
    real(RP), intent(in) :: LATUY(0:IA,  JA)
    real(RP), intent(in) :: LATXV(  IA,0:JA)
    real(RP), intent(in) :: LATUV(0:IA,0:JA)
    real(RP), intent(in) :: TOPO  (  IA,  JA)
    real(RP), intent(in) :: LSMASK(  IA,  JA)
    !---------------------------------------------------------------------------

    AXIS_HGT   (:,:,:) = CZ(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_HGTWXY(:,:,:) = FZ(KS-1:KE,ISB2:IEB2,JSB2:JEB2)

    AXIS_LON  (:,:) = LON  (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LONUY(:,:) = LONUY(ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LONXV(:,:) = LONXV(ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LONUV(:,:) = LONUV(ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LAT  (:,:) = LAT  (ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LATUY(:,:) = LATUY(ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LATXV(:,:) = LATXV(ISB2:IEB2,JSB2:JEB2) / D2R
    AXIS_LATUV(:,:) = LATUV(ISB2:IEB2,JSB2:JEB2) / D2R

    AXIS_TOPO  (:,:) = TOPO  (ISB2:IEB2,JSB2:JEB2)
    AXIS_LSMASK(:,:) = LSMASK(ISB2:IEB2,JSB2:JEB2)

    set_coordinates = .true.

    return
  end subroutine FILE_CARTESC_set_coordinates_atmos

  !-----------------------------------------------------------------------------
  !> set area and volume
  subroutine FILE_CARTESC_set_areavol_atmos( &
       AREA,   AREAZUY_X, AREAZXV_Y, &
               AREAWUY_X, AREAWXV_Y, &
       AREAUY, AREAZXY_X, AREAZUV_Y, &
       AREAXV, AREAZUV_X, AREAZXY_Y, &
       VOL, VOLWXY, VOLZUY, VOLZXV   )
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none
    real(RP), intent(in) :: AREA     (     IA,JA)
    real(RP), intent(in) :: AREAZUY_X(  KA,IA,JA)
    real(RP), intent(in) :: AREAZXV_Y(  KA,IA,JA)
    real(RP), intent(in) :: AREAWUY_X(0:KA,IA,JA)
    real(RP), intent(in) :: AREAWXV_Y(0:KA,IA,JA)
    real(RP), intent(in) :: AREAUY   (     IA,JA)
    real(RP), intent(in) :: AREAZXY_X(  KA,IA,JA)
    real(RP), intent(in) :: AREAZUV_Y(  KA,IA,JA)
    real(RP), intent(in) :: AREAXV   (     IA,JA)
    real(RP), intent(in) :: AREAZUV_X(  KA,IA,JA)
    real(RP), intent(in) :: AREAZXY_Y(  KA,IA,JA)
    real(RP), intent(in) :: VOL   (  KA,IA,JA)
    real(RP), intent(in) :: VOLWXY(0:KA,IA,JA)
    real(RP), intent(in) :: VOLZUY(  KA,IA,JA)
    real(RP), intent(in) :: VOLZXV(  KA,IA,JA)

    AXIS_AREA     (:,:)   = AREA     (        ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAZUY_X(:,:,:) = AREAZUY_X(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAZXV_Y(:,:,:) = AREAZXV_Y(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAWUY_X(:,:,:) = AREAWUY_X(KS-1:KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAWXV_Y(:,:,:) = AREAWXV_Y(KS-1:KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAUY   (:,:)   = AREAUY   (        ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAZXY_X(:,:,:) = AREAZXY_X(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAZUV_Y(:,:,:) = AREAZUV_Y(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAXV   (:,:)   = AREAXV   (        ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAZUV_X(:,:,:) = AREAZUV_X(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_AREAZXY_Y(:,:,:) = AREAZXY_Y(KS  :KE,ISB2:IEB2,JSB2:JEB2)

    AXIS_VOL   (:,:,:) = VOL   (KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_VOLWXY(:,:,:) = VOLWXY(KS-1:KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_VOLZUY(:,:,:) = VOLZUY(KS  :KE,ISB2:IEB2,JSB2:JEB2)
    AXIS_VOLZXV(:,:,:) = VOLZXV(KS  :KE,ISB2:IEB2,JSB2:JEB2)

    return
  end subroutine FILE_CARTESC_set_areavol_atmos

  !-----------------------------------------------------------------------------
  !> set volume for ocean
  subroutine FILE_CARTESC_set_coordinates_ocean( &
       VOL )
    implicit none

    real(RP), intent(in) :: VOL(OKA,OIA,OJA)
    !---------------------------------------------------------------------------

    AXIS_VOLO(:,:,:) = VOL(OKS:OKE,ISB2:IEB2,JSB2:JEB2)

    return
  end subroutine FILE_CARTESC_set_coordinates_ocean

  !-----------------------------------------------------------------------------
  !> set volume for land
  subroutine FILE_CARTESC_set_coordinates_land( &
       VOL )
    implicit none

    real(RP), intent(in) :: VOL(LKA,LIA,LJA)
    !---------------------------------------------------------------------------

    AXIS_VOLL(:,:,:) = VOL(LKS:LKE,ISB2:IEB2,JSB2:JEB2)

    return
  end subroutine FILE_CARTESC_set_coordinates_land

  !-----------------------------------------------------------------------------
  !> set volume for urban
  subroutine FILE_CARTESC_set_coordinates_urban( &
       VOL )
    implicit none

    real(RP), intent(in) :: VOL(UKA,UIA,UJA)
    !---------------------------------------------------------------------------

    AXIS_VOLU(:,:,:) = VOL(UKS:UKE,ISB2:IEB2,JSB2:JEB2)

    return
  end subroutine FILE_CARTESC_set_coordinates_urban

  !-----------------------------------------------------------------------------
  !> check coordinates in the file
  subroutine FILE_CARTESC_check_coordinates_name( &
       basename,                  &
       atmos, ocean, land, urban, &
       transpose                  )
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

    call FILE_CARTESC_open( basename, fid )


    call FILE_CARTESC_check_coordinates_id( fid,                           & ! [IN]
                                            atmos_, ocean_, land_, urban_, & ! [IN]
                                            transpose_                     ) ! [IN]

    return
  end subroutine FILE_CARTESC_check_coordinates_name

  !-----------------------------------------------------------------------------
  !> check coordinates in the file
  subroutine FILE_CARTESC_check_coordinates_id( &
       fid,                       &
       atmos, ocean, land, urban, &
       transpose                  )
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ, &
       ATMOS_GRID_CARTESC_CX, &
       ATMOS_GRID_CARTESC_CY
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_CZ
    use scale_land_grid_cartesC, only: &
       LAND_GRID_CARTESC_CZ
    use scale_urban_grid_cartesC, only: &
       URBAN_GRID_CARTESC_CZ
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

    LOG_NEWLINE
    LOG_INFO("FILE_CARTESC_check_coordinates_id",*) 'Check consistency of axis '

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

    call FILE_CARTESC_read_var_1D( fid, 'x',  'X', buffer_x(:) )
    call FILE_CARTESC_read_var_1D( fid, 'y',  'Y', buffer_y(:) )
    call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
    call check_1d( ATMOS_GRID_CARTESC_CX(ISB:IEB), buffer_x(ISB:IEB), 'x' )
    call check_1d( ATMOS_GRID_CARTESC_CY(JSB:JEB), buffer_y(JSB:JEB), 'y' )

    if ( set_coordinates ) then
       call FILE_CARTESC_read_var_2D( fid, 'lon', 'XY', buffer_xy(:,:) )
       call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
       call check_2d( AXIS_LON(XSB:XEB,YSB:YEB), buffer_xy(ISB:IEB,JSB:JEB), 'lon' )

       call FILE_CARTESC_read_var_2D( fid, 'lat', 'XY', buffer_xy(:,:) )
       call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
       call check_2d( AXIS_LAT(XSB:XEB,YSB:YEB), buffer_xy(ISB:IEB,JSB:JEB), 'lat' )
    endif

    if ( atmos_ ) then
       call FILE_CARTESC_read_var_1D( fid, 'z', 'Z', buffer_z(:) )
       if ( .not. transpose_ ) then
          call FILE_CARTESC_read_var_3D( fid, 'height', 'ZXY', buffer_zxy(:,:,:) )
       endif
       call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
       call check_1d( ATMOS_GRID_CARTESC_CZ(KS:KE), buffer_z(KS:KE), 'z' )
       if ( .not. transpose_ ) then
          call check_3d( AXIS_HGT(:,XSB:XEB,YSB:YEB), buffer_zxy(KS:KE,ISB:IEB,JSB:JEB), 'height', transpose_ )
       endif
    endif

    if ( ocean_ ) then
       call FILE_CARTESC_read_var_1D( fid, 'oz', 'OZ', buffer_o(:) )
       call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
       call check_1d( OCEAN_GRID_CARTESC_CZ(OKS:OKE), buffer_o(OKS:OKE), 'oz' )
    endif

    if ( land_ ) then
       call FILE_CARTESC_read_var_1D( fid, 'lz',  'LZ', buffer_l(:) )
       call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
       call check_1d( LAND_GRID_CARTESC_CZ(LKS:LKE), buffer_l(LKS:LKE), 'lz' )
    endif

    if ( urban_ ) then
       call FILE_CARTESC_read_var_1D( fid, 'uz',  'UZ', buffer_u(:) )
       call FILE_CARTESC_flush( fid ) ! for non-blocking I/O
       call check_1d( URBAN_GRID_CARTESC_CZ(UKS:UKE), buffer_u(UKS:UKE), 'uz' )
    endif

    return
  end subroutine FILE_CARTESC_check_coordinates_id

  !-----------------------------------------------------------------------------
  !> open a netCDF file for read
  subroutine FILE_CARTESC_open( &
       basename, &
       fid,      &
       aggregate )
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_AGGREGATE, &
       FILE_Open
    use scale_prc, only: &
       PRC_myrank
    implicit none

    character(len=*), intent(in)  :: basename !< basename of the file
    integer,          intent(out) :: fid      !< file ID
    logical,          intent(in), optional :: aggregate
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_O_NetCDF', 2)

    call FILE_Open( basename,            & ! [IN]
                    fid,                 & ! [OUT]
                    aggregate=aggregate, & ! [IN]
                    rankid=PRC_myrank    ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_open

  !-----------------------------------------------------------------------------
  !> Create/open a netCDF file
  subroutine FILE_CARTESC_create( &
       basename, title, datatype, &
       fid,                       &
       date, subsec,              &
       haszcoord,                 &
       append, aggregate, single  )
    use scale_file_h, only: &
       FILE_REAL8, &
       FILE_REAL4
    use scale_file, only: &
       FILE_AGGREGATE,    &
       FILE_Create,       &
       FILE_get_CFtunits
    use scale_calendar, only: &
       CALENDAR_get_name
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_time, only: &
       NOWDATE   => TIME_NOWDATE, &
       NOWSUBSEC => TIME_NOWSUBSEC
    use scale_prc_cartesC, only: &
       PRC_2Drank,     &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y
    implicit none

    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: title    !< title    of the file
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)
    integer,          intent(out) :: fid      !< file ID
    integer,          intent(in), optional :: date(6)   !< ymdhms of the time
    real(DP),         intent(in), optional :: subsec    !< subsec of the time
    logical,          intent(in), optional :: append    !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: haszcoord !< switch whether include zcoordinate or not (default=true)
    logical,          intent(in), optional :: aggregate
    logical,          intent(in), optional :: single

    integer                :: dtype
    logical                :: append_sw
    character(len=34)      :: tunits
    character(len=H_SHORT) :: calendar
    real(DP)               :: subsec_
    integer                :: rank_x, rank_y
    integer                :: num_x, num_y
    logical                :: fileexisted
    logical                :: aggregate_
    logical                :: single_
    integer                :: date_(6)
    !---------------------------------------------------------------------------

    prof = .true.
    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    end if

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
          LOG_ERROR("FILE_CARTESC_create",*) 'unsupported data type. Check!', trim(datatype)
          call PRC_abort
       endif
    endif

    append_sw = .false.
    if ( present(append) ) then
       append_sw = append
    endif

    ! create a netCDF file if not already existed. Otherwise, open it.
    if ( present(date) ) then
       date_(:) = date(:)
    else
       date_(:) = NOWDATE(:)
    end if
    if ( date_(1) > 0 ) then
       call FILE_get_CFtunits( date_(:), tunits )
       call CALENDAR_get_name( calendar )
    else
       tunits = 'seconds'
       calendar = ''
    endif


    ! check to use PnetCDF I/O
    if ( present(aggregate) ) then
       aggregate_ = aggregate
    else
       aggregate_ = FILE_AGGREGATE
    endif

    call FILE_Create( basename,                & ! [IN]
                      title,                   & ! [IN]
                      H_SOURCE,                & ! [IN]
                      H_INSTITUTE,             & ! [IN]
                      fid,                     & ! [OUT]
                      fileexisted,             & ! [OUT]
                      rankid     = PRC_myrank, & ! [IN]
                      single     = single_,    & ! [IN]
                      aggregate  = aggregate_, & ! [IN]
                      time_units = tunits,     & ! [IN]
                      calendar   = calendar,   & ! [IN]
                      append     = append_sw   ) ! [IN]


    if ( PRC_myrank /= 0 .and. single_ ) then
       fid = -1
    else if ( .not. fileexisted ) then ! do below only once when file is created

       File_axes_written(fid) = .false.  ! indicating axes have not been written yet

       if ( present( haszcoord ) ) then
          File_haszcoord(fid) = haszcoord
       else
          File_haszcoord(fid) = .true.
       endif

       if ( aggregate_ ) then
          rank_x = 0
          rank_y = 0
          num_x = 1
          num_y = 1
       else
          rank_x = PRC_2Drank(PRC_myrank,1)
          rank_y = PRC_2Drank(PRC_myrank,2)
          num_x = PRC_NUM_X
          num_y = PRC_NUM_Y
       end if

       if ( present( subsec ) ) then
          subsec_ = subsec
       else
          subsec_= NOWSUBSEC
       end if

       call FILE_CARTESC_put_globalAttributes( fid,                            & ! [IN]
                                               rank_x, rank_y,                 & ! [IN]
                                               num_x, num_y,                   & ! [IN]
                                               PRC_PERIODIC_X, PRC_PERIODIC_Y, & ! [IN]
                                               KMAX, OKMAX, LKMAX, UKMAX,      & ! [IN]
                                               IMAXG, JMAXG,                   & ! [IN]
                                               KHALO, IHALO, JHALO,            & ! [IN]
                                               subsec_, tunits, calendar       ) ! [IN]

       call FILE_CARTESC_def_axes( fid,                & ! [IN]
                                   dtype,              & ! [IN]
                                   File_haszcoord(fid) ) ! [IN]

    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)
    prof = .false.

    return
  end subroutine FILE_CARTESC_create

  !-----------------------------------------------------------------------------
  !> Exit netCDF file define mode
  subroutine FILE_CARTESC_enddef( fid )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened,    &
       FILE_EndDef,    &
       FILE_Flush,     &
       FILE_Attach_Buffer
    implicit none

    integer, intent(in) :: fid  !< file ID

    integer :: start(3)
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

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
          call FILE_Attach_Buffer( fid, write_buf_amount(fid) )
       endif

       File_axes_written(fid) = .true.
    endif

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_enddef

  !-----------------------------------------------------------------------------
  !> Flush all pending requests to a netCDF file (PnetCDF only)
  subroutine FILE_CARTESC_flush( fid )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened, &
       FILE_Flush
    implicit none

    integer, intent(in) :: fid  !< file ID
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( FILE_get_AGGREGATE(fid) ) then
       call FILE_Flush( fid ) ! flush all pending read/write requests
    end if

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_flush

  !-----------------------------------------------------------------------------
  !> Close a netCDF file
  subroutine FILE_CARTESC_close( fid )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened,    &
       FILE_Close,     &
       FILE_Flush,     &
       FILE_Detach_Buffer
    implicit none

    integer, intent(in) :: fid  !< file ID
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    if ( FILE_get_AGGREGATE(fid) ) then
       call FILE_Flush( fid )        ! flush all pending read/write requests
       if ( write_buf_amount(fid) > 0 ) then
          call FILE_Detach_Buffer( fid ) ! detach PnetCDF aggregation buffer
          write_buf_amount(fid) = 0      ! reset write request amount
       endif
    endif

    call FILE_Close( fid ) ! [IN]

    call PROF_rapend  ('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_close

  !-----------------------------------------------------------------------------
  !> interface FILE_CARTESC_read
  !!    Read data from file
  !!    This routine is a wrapper of the lower primitive routines
  !<
  !-----------------------------------------------------------------------------
  subroutine FILE_CARTESC_read_1D( &
       basename, varname, &
       dim_type,          &
       var,               &
       step,              &
       aggregate,         &
       allow_missing      )
    implicit none
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)

    real(RP),         intent(out) :: var(:)   !< value of the variable

    integer,          intent(in), optional :: step     !< step number
    logical,          intent(in), optional :: aggregate
    logical,          intent(in), optional :: allow_missing

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( basename,           & ! [IN]
                            fid,                & ! [OUT]
                            aggregate=aggregate ) ! [IN]

    call FILE_CARTESC_read_var_1D( fid, varname, dim_type,     & ! [IN]
                                   var(:),                     & ! [OUT]
                                   step=step,                  & ! [IN]
                                   allow_missing=allow_missing ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_1D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILE_CARTESC_read_2D( &
       basename, varname, &
       dim_type,          &
       var,               &
       step,              &
       aggregate,         &
       allow_missing      )
    implicit none
    character(len=*), intent(in)  :: basename !< basename of the file
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)

    real(RP),         intent(out) :: var(:,:) !< value of the variable

    integer,          intent(in), optional :: step     !< step number
    logical,          intent(in), optional :: aggregate
    logical,          intent(in), optional :: allow_missing

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( basename,           & ! [IN]
                            fid,                & ! [OUT]
                            aggregate=aggregate ) ! [IN]

    call FILE_CARTESC_read_var_2D( fid, varname, dim_type,     & ! [IN]
                                   var(:,:),                   & ! [OUT]
                                   step=step,                  & ! [IN]
                                   allow_missing=allow_missing ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILE_CARTESC_read_3D( &
       basename, varname, &
       dim_type,          &
       var,               &
       step,              &
       aggregate,         &
       allow_missing      )
    implicit none
    character(len=*), intent(in)  :: basename   !< basename of the file
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: dim_type   !< dimension type (Z/X/Y/T)

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable

    integer,          intent(in), optional :: step       !< step number
    logical,          intent(in), optional :: aggregate
    logical,          intent(in), optional :: allow_missing

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( basename,           & ! [IN]
                            fid,                & ! [OUT]
                            aggregate=aggregate ) ! [IN]

    call FILE_CARTESC_read_var_3D( fid, varname, dim_type,     & ! [IN]
                                   var(:,:,:),                 & ! [OUT]
                                   step=step,                  & ! [IN]
                                   allow_missing=allow_missing ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_3D

  !-----------------------------------------------------------------------------
  !> Read 4D data from file
  subroutine FILE_CARTESC_read_4D( &
       basename, varname, &
       dim_type, step,    &
       var,               &
       aggregate,         &
       allow_missing      )
    implicit none
    character(len=*), intent(in)  :: basename     !< basename of the file
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: dim_type     !< dimension type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable

    logical,          intent(in), optional :: aggregate
    logical,          intent(in), optional :: allow_missing

    integer :: fid
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( basename,           & ! [IN]
                            fid,                & ! [OUT]
                            aggregate=aggregate ) ! [IN]

    call FILE_CARTESC_read_var_4D( fid, varname, dim_type, step, & ! [IN]
                                   var(:,:,:,:),                 & ! [OUT]
                                   allow_missing=allow_missing   ) ! [IN]

    call FILE_CARTESC_close( fid )

    return
  end subroutine FILE_CARTESC_read_4D

  !-----------------------------------------------------------------------------
  !> Read 1D data from file
  subroutine FILE_CARTESC_read_var_1D( &
       fid, varname, &
       dim_type,     &
       var,          &
       step,         &
       allow_missing )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened, &
       FILE_Read
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    use mpi
    implicit none
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)

    real(RP),         intent(out) :: var(:)   !< value of the variable

    integer,          intent(in), optional :: step     !< step number
    logical,          intent(in), optional :: allow_missing

    integer :: vsize
    integer :: dim1_S, dim1_E
    integer :: start(1)   ! start offset of globale variable
    integer :: count(1)   ! request length to the global variable
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_I_NetCDF', 2)

    LOG_INFO("FILE_CARTESC_read_var_1D",'(1x,2A)') 'Read from file (1D), name : ', trim(varname)

    if ( FILE_get_aggregate(fid) ) then
       ! read data and halos into the local buffer
       if    ( dim_type == 'Z' ) then
          vsize = KA
          dim1_S = KS
          dim1_E = KE
          start(1) = 1
       elseif( dim_type == 'OZ' ) then
          vsize = OKA
          dim1_S = OKS
          dim1_E = OKE
          start(1) = 1
       elseif( dim_type == 'LZ' ) then
          vsize = LKA
          dim1_S = LKS
          dim1_E = LKE
          start(1) = 1
       elseif( dim_type == 'UZ' ) then
          vsize = UKA
          dim1_S = UKS
          dim1_E = UKE
          start(1) = 1
       elseif( dim_type == 'X' .OR. dim_type == 'CX' ) then
          vsize = IA
          dim1_S = 1
          dim1_E = IA
          start(1) = IS_inG - IHALO
       elseif( dim_type == 'Y' .OR. dim_type == 'CY' ) then
          vsize = JA
          dim1_S = 1
          dim1_E = JA
          start(1) = JS_inG - JHALO
       else
          LOG_ERROR("FILE_CARTESC_read_var_1D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_1D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       count(1) = dim1_E - dim1_S + 1
       call FILE_Read( fid, varname,                               & ! (in)
            var(dim1_S:dim1_E),                                    & ! (out)
            step=step, allow_missing=allow_missing,                & ! (in)
            ntypes=count(1), dtype=etype, start=start, count=count ) ! (in)

    else
       if    ( dim_type == 'Z' ) then
          vsize = KA
          dim1_S = KS
          dim1_E = KE
       elseif( dim_type == 'OZ' ) then
          vsize = OKA
          dim1_S = OKS
          dim1_E = OKE
       elseif( dim_type == 'LZ' ) then
          vsize = LKA
          dim1_S = LKS
          dim1_E = LKE
       elseif( dim_type == 'UZ' ) then
          vsize = UKA
          dim1_S = UKS
          dim1_E = UKE
       elseif( dim_type == 'X' ) then
          vsize = IA
          dim1_S = ISB
          dim1_E = IEB
       elseif( dim_type == 'CX' ) then
          vsize = IA
          dim1_S = 1
          dim1_E = IA
       elseif( dim_type == 'Y' ) then
          vsize = JA
          dim1_S = JSB
          dim1_E = JEB
       elseif( dim_type == 'CY' ) then
          vsize = JA
          dim1_S = 1
          dim1_E = JA
       else
          LOG_ERROR("FILE_CARTESC_read_var_1D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_1D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       call FILE_Read( fid, varname, var(dim1_S:dim1_E), step=step )
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_1D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILE_CARTESC_read_var_2D( &
       fid, varname, &
       dim_type,     &
       var,          &
       step,         &
       allow_missing )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened, &
       FILE_Read
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable
    character(len=*), intent(in)  :: dim_type !< dimension type (Z/X/Y)

    real(RP),         intent(out) :: var(:,:) !< value of the variable

    integer,          intent(in), optional :: step     !< step number
    logical,          intent(in), optional :: allow_missing

    integer :: vsize
    integer :: ntypes, dtype
    integer, pointer :: start(:), count(:)
    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_I_NetCDF', 2)

    LOG_INFO("FILE_CARTESC_read_var_2D",'(1x,2A)') 'Read from file (2D), name : ', trim(varname)

    if ( FILE_get_AGGREGATE(fid) ) then

       ! read data and halos into the local buffer
       if    ( dim_type == 'XY' ) then
          vsize = IA * JA
          ntypes = IA * JA
          dtype = etype
          start => startXY
          count => countXY
       elseif( dim_type == 'ZX' ) then
          ! Because KHALO is not saved in files, we use centerTypeZX, an MPI
          ! derived datatype to describe the layout of local read buffer
          vsize = KA * IA
          ntypes = 1
          dtype = centerTypeZX
          start => startZX
          count => countZX
       else
          LOG_ERROR("FILE_CARTESC_read_var_2D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_2D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       call FILE_Read( fid, varname,                             & ! (in)
            var(:,:),                                            & ! (out)
            step=step, allow_missing=allow_missing,              & ! (in)
            ntypes=ntypes, dtype=dtype, start=start, count=count ) ! (in)

    else

       if    ( dim_type == 'XY' ) then
          vsize = IA * JA
          dim1_S = ISB
          dim1_E = IEB
          dim2_S = JSB
          dim2_E = JEB
       elseif( dim_type == 'ZX' ) then
          vsize = KA * IA
          dim1_S = KS
          dim1_E = KE
          dim2_S = ISB
          dim2_E = IEB
       else
          LOG_ERROR("FILE_CARTESC_read_var_2D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_2D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       call FILE_Read( fid, varname, var(dim1_S:dim1_E,dim2_S:dim2_E), step=step )
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILE_CARTESC_read_var_3D( &
       fid, varname, &
       dim_type,     &
       var,          &
       step,         &
       allow_missing )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened, &
       FILE_Read
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none
    integer,          intent(in)  :: fid        !< file ID
    character(len=*), intent(in)  :: varname    !< name of the variable
    character(len=*), intent(in)  :: dim_type   !< dimension type (Z/X/Y/T)

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable

    integer,          intent(in), optional :: step       !< step number
    logical,          intent(in), optional :: allow_missing

    integer :: vsize
    integer :: ntypes, dtype
    integer, pointer :: start(:), count(:)
    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    integer :: dim3_S, dim3_E
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_I_NetCDF', 2)

    LOG_INFO("FILE_CARTESC_read_var_3D",'(1x,2A)') 'Read from file (3D), name : ', trim(varname)

    if ( FILE_get_AGGREGATE(fid) ) then

       ! read data and halos into the local buffer
       ! Because KHALO is not saved in files, we use mpi derived datatypes to
       ! describe the layout of local read buffer
       if(      dim_type == 'ZXY'  &
           .or. dim_type == 'ZXHY' &
           .or. dim_type == 'ZXYH' ) then
          vsize = KA * IA * JA
          ntypes = 1
          dtype = centerTypeZXY
          start => startZXY
          count => countZXY
       elseif( dim_type == 'ZHXY' ) then
          vsize = KA * IA * JA
          ntypes = 1
          dtype = centerTypeZHXY
          start => startZHXY
          count => countZHXY
       elseif( dim_type == 'XYT' ) then
          if ( .not. present(step) ) then
             LOG_ERROR("FILE_CARTESC_read_var_3D",*) 'step is necessary for "XYT"'
             call PRC_abort
          end if
          vsize = IA * JA * step
          ntypes = IA * JA * step
          dtype = etype
          startXY(3) = 1
          countXY(3) = step
          start => startXY
          count => countXY
       elseif( dim_type == 'OXY' ) then
          vsize = OKA * OIA * OJA
          ntypes = 1
          dtype = centerTypeOCEAN
          start => startOCEAN
          count => countOCEAN
       elseif( dim_type == 'LXY' ) then
          vsize = LKA * LIA * LJA
          ntypes = 1
          dtype = centerTypeLAND
          start => startLAND
          count => countLAND
       elseif( dim_type == 'UXY' ) then
          vsize = UKA * UIA * UJA
          ntypes = 1
          dtype = centerTypeURBAN
          start => startURBAN
          count => countURBAN
       else
          LOG_ERROR("FILE_CARTESC_read_var_3D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_3D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       call FILE_Read( fid, varname,                             & ! (in)
            var(:,:,:),                                          & ! (out)
            step=step, allow_missing=allow_missing,              & ! (in)
            ntypes=ntypes, dtype=dtype, start=start, count=count ) ! (in)

    else
       if(      dim_type == 'ZXY'  &
           .or. dim_type == 'ZXHY' &
           .or. dim_type == 'ZXYH' ) then
          vsize = KA * IA * JA
          dim1_S = KS
          dim1_E = KE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif( dim_type == 'ZHXY' ) then
          vsize = KA * IA * JA
          dim1_S = KS-1
          dim1_E = KE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif( dim_type == 'XYT' ) then
          if ( .not. present(step) ) then
             LOG_ERROR("FILE_CARTESC_read_var_3D",*) 'step is necessary for "XYT"'
             call PRC_abort
          end if
          vsize = IA * JA * step
          dim1_S = ISB
          dim1_E = IEB
          dim2_S = JSB
          dim2_E = JEB
          dim3_S = 1
          dim3_E = step
       elseif( dim_type == 'OXY' ) then
          vsize = OKA * OIA * OJA
          dim1_S = OKS
          dim1_E = OKE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif( dim_type == 'LXY' ) then
          vsize = LKA * LIA * LJA
          dim1_S = LKS
          dim1_E = LKE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif( dim_type == 'UXY' ) then
          vsize = UKA * UIA * UJA
          dim1_S = UKS
          dim1_E = UKE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       else
          LOG_ERROR("FILE_CARTESC_read_var_3D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_3D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       call FILE_Read( fid, varname, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                       step=step, allow_missing=allow_missing                        )

    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_3D

  !-----------------------------------------------------------------------------
  !> Read 4D data from file
  subroutine FILE_CARTESC_read_var_4D( &
       fid, varname, &
       dim_type,     &
       step,         &
       var,          &
       allow_missing )
    use scale_file, only: &
       FILE_get_AGGREGATE, &
       FILE_opened, &
       FILE_Read
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_NUM_X, &
       PRC_NUM_Y
    implicit none
    integer,          intent(in)  :: fid          !< file ID
    character(len=*), intent(in)  :: varname      !< name of the variable
    character(len=*), intent(in)  :: dim_type     !< dimension type (Z/X/Y/Time)
    integer,          intent(in)  :: step         !< step number

    real(RP),         intent(out) :: var(:,:,:,:) !< value of the variable

    logical,          intent(in), optional :: allow_missing

    integer :: vsize
    integer :: dtype
    integer, pointer :: start(:), count(:)
    integer :: dim1_S, dim1_E
    integer :: dim2_S, dim2_E
    integer :: dim3_S, dim3_E
    integer :: dim4_S, dim4_E
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_I_NetCDF', 2)

    LOG_INFO("FILE_CARTESC_read_var_4D",'(1x,2A)') 'Read from file (4D), name : ', trim(varname)

    if ( FILE_get_AGGREGATE(fid) ) then
       ! read data and halos into the local buffer
       if (      dim_type == 'ZXYT'  &
            .or. dim_type == 'ZXHYT' &
            .or. dim_type == 'ZXYHT' ) then
          vsize = KA * IA * JA * step
          dtype = centerTypeZXY
          start => startZXY
          count => countZXY
       elseif ( dim_type == 'ZHXYT' ) then
          vsize = KA * IA * JA * step
          dtype = centerTypeZHXY
          start => startZHXY
          count => countZHXY
       elseif ( dim_type == 'OXYT' ) then
          vsize = OKA * OIA * OJA * step
          dtype = centerTypeOCEAN
          start => startOCEAN
          count => countOCEAN
       elseif ( dim_type == 'LXYT' ) then
          vsize = LKA * LIA * LJA * step
          dtype = centerTypeLAND
          start => startLAND
          count => countLAND
       elseif ( dim_type == 'LXYT' ) then
          vsize = LKA * LIA * LJA * step
          dtype = centerTypeLAND
          start => startLAND
          count => countLAND
       elseif ( dim_type == 'UXYT' ) then
          vsize = UKA * UIA * UJA * step
          dtype = centerTypeURBAN
          start => startURBAN
          count => countURBAN
       else
          LOG_ERROR("FILE_CARTESC_read_var_4D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_4D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       start(4) = 1
       count(4) = step
       call FILE_Read( fid, varname,                           & ! (in)
            var(:,:,:,:),                                      & ! (out)
            allow_missing=allow_missing,                       & ! (in)
            ntypes=step, dtype=dtype, start=start, count=count ) ! (in)

    else
       if (      dim_type == 'ZXYT'  &
            .or. dim_type == 'ZXHYT' &
            .or. dim_type == 'ZXYHT' ) then
          vsize = KA * IA * JA * step
          dim1_S = KS
          dim1_E = KE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif ( dim_type == 'ZHXYT' ) then
          vsize = KA * IA * JA * step
          dim1_S = KS-1
          dim1_E = KE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif ( dim_type == 'OXYT' ) then
          vsize = OKA * OIA * OJA * step
          dim1_S = OKS
          dim1_E = OKE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif ( dim_type == 'LXYT' ) then
          vsize = LKA * LIA * LJA * step
          dim1_S = LKS
          dim1_E = LKE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       elseif ( dim_type == 'OXYT' ) then
          vsize = LKA * LIA * LJA * step
          dim1_S = LKS
          dim1_E = LKE
          dim2_S = ISB
          dim2_E = IEB
          dim3_S = JSB
          dim3_E = JEB
       else
          LOG_ERROR("FILE_CARTESC_read_var_4D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
          call PRC_abort
       endif

       if ( size(var) .ne. vsize ) then
          LOG_ERROR("FILE_CARTESC_read_var_4D",*) 'size of var is invalid: ', trim(varname), size(var), vsize
          call PRC_abort
       end if
       dim4_S   = 1
       dim4_E   = step
       call FILE_Read( fid, varname,                                     & ! (in)
            var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E,dim4_S:dim4_E) ) ! (out)
    endif

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_var_4D

  !-----------------------------------------------------------------------------
  !> Read 2D data from file
  subroutine FILE_CARTESC_read_auto_2D( &
       fid, varname, &
       var,          &
       step, existed )
    use scale_file, only: &
       FILE_opened, &
       FILE_get_shape, &
       FILE_get_dataInfo, &
       FILE_get_attribute, &
       FILE_read
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name of the variable

    real(RP),         intent(out) :: var(:,:) !< value of the variable

    integer, intent(in), optional :: step     !< step number

    logical, intent(out), optional :: existed

    integer :: dims(2)
    integer :: halos(2)
    integer :: start(2)
    integer :: count(2)
    character(len=H_SHORT) :: dnames(2)

    integer :: nx, ny
    integer :: n

    logical :: existed2

    intrinsic size
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_I_NetCDF', 2)

    LOG_INFO("FILE_CARTESC_read_auto_2D",'(1x,2A)') 'Read from file (2D), name : ', trim(varname)

    call FILE_get_dataInfo( fid, varname, dim_name=dnames(:), existed=existed2 )

    if ( present( existed ) ) then
       existed = existed2
       if ( .not. existed2 ) return
    end if

    if ( .not. existed2 ) then
       LOG_ERROR("FILE_CARTESC_read_auto_2D",*) 'variable not found: ', trim(varname)
       call PRC_abort
    end if

    call FILE_get_shape( fid, varname, dims(:) )
    nx = size(var,1)
    ny = size(var,2)

    if ( nx==dims(1) .and. ny==dims(2) ) then
       start(:) = (/1,1/)
    else
       do n = 1, 2
          call FILE_get_attribute( fid, dnames(n), "halo_local", halos(:), existed=existed2 )
          if ( existed2 ) then
             start(n) = halos(1) + 1
          else
             start(n) = 1
          end if
       end do
    end if
    count(:) = (/nx,ny/)

    call FILE_read( fid, varname, var(:,:), step=step, start=start(:), count=count(:) )
    call FILE_CARTESC_flush( fid )

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_auto_2D

  !-----------------------------------------------------------------------------
  !> Read 3D data from file
  subroutine FILE_CARTESC_read_auto_3D( &
       fid, varname, &
       var,          &
       step, existed )
    use scale_file, only: &
       FILE_opened, &
       FILE_get_shape, &
       FILE_get_dataInfo, &
       FILE_get_attribute, &
       FILE_read
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,          intent(in)  :: fid        !< file ID
    character(len=*), intent(in)  :: varname    !< name of the variable

    real(RP),         intent(out) :: var(:,:,:) !< value of the variable

    integer, intent(in), optional :: step       !< step number

    logical, intent(out), optional :: existed

    integer :: dims(3)
    integer :: halos(3)
    integer :: start(3)
    integer :: count(3)
    character(len=H_SHORT) :: dnames(3)

    logical :: existed2

    real(RP), allocatable :: buf(:,:,:)
    integer :: nx, ny, nz
    integer :: n
    integer :: k, i, j

    intrinsic size
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_I_NetCDF', 2)

    LOG_INFO("FILE_CARTESC_read_auto_3D",'(1x,2A)') 'Read from file (3D), name : ', trim(varname)

    call FILE_get_dataInfo( fid, varname, dim_name=dnames(:), existed=existed2 )

    if ( present(existed) ) then
       existed = existed2
       if ( .not. existed2 ) return
    end if

    if ( .not. existed2 ) then
       LOG_ERROR("FILE_CARTESC_read_auto_3D",*) 'variable not found: ', trim(varname)
       call PRC_abort
    end if

    call FILE_get_shape( fid, varname, dims(:) )
    nz = size(var,1)
    nx = size(var,2)
    ny = size(var,3)

    if      ( ( dnames(1)(1:1)=="z" .or. dnames(1)(2:2)=="z" ) .and. dnames(2)(1:1)=="x" .and. dnames(3)(1:1)=="y" ) then
       if ( nz==dims(1) .and. nx==dims(2) .and. ny==dims(3) ) then
          start(:) = (/1,1,1/)
       else if ( dnames(1)=="zh" .and. nz+1==dims(1) .and. nx==dims(2) .and. ny==dims(3) ) then
          start(:) = (/2,1,1/)
       else
          do n = 1, 3
             call FILE_get_attribute( fid, dnames(n), "halo_local", halos(:), existed=existed2 )
             if ( existed2 ) then
                start(n) = halos(1) + 1
             else if ( dnames(n)=="zh" ) then
                start(n) = 2
             else
                start(n) = 1
             end if
          end do
       end if
       count(:) = (/nz,nx,ny/)
       call FILE_read( fid, varname, var(:,:,:), step=step, start=start(:), count=count(:) )
       call FILE_CARTESC_flush( fid )
    else if ( dnames(1)(1:1)=="x" .and. dnames(2)(1:1)=="y" .and. ( dnames(3)(1:1)=="z" .or. dnames(3)(2:2)=="z" ) ) then
       allocate( buf(nx,ny,nz) )
       if ( nx==dims(1) .and. ny==dims(2) .and. nz==dims(3) ) then
          start(:) = (/1,1,1/)
       else if ( nx==dims(1) .and. ny==dims(2) .and. nz+1==dims(3) .and. dnames(3)=="zh" ) then
          start(:) = (/1,1,2/)
       else
          do n = 1, 3
             call FILE_get_attribute( fid, dnames(n), "halo_local", halos(:), existed=existed2 )
             if ( existed2 ) then
                start(n) = halos(1) + 1
             else if ( dnames(n)=="zh" ) then
                start(n) = 2
             else
                start(n) = 1
             end if
          end do
       end if
       count(:) = (/nx,ny,nz/)
       call FILE_read( fid, varname, buf(:,:,:), step=step, start=start(:), count=count(:) )
       call FILE_CARTESC_flush( fid )

       !$omp parallel do
       do j = 1, ny
       do i = 1, nx
       do k = 1, nz
          var(k,i,j) = buf(i,j,k)
       end do
       end do
       end do
       deallocate(buf)
    else
       LOG_ERROR("FILE_CARTESC_read_auto_3D",*) 'invalid dimension'
       call PRC_abort
    end if

    call PROF_rapend  ('FILE_I_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_read_auto_3D

  !-----------------------------------------------------------------------------
  !> interface FILE_CARTESC_write
  !!   Write data to file
  !!   This routine is a wrapper of the lowere primitive routines
  !<
  !-----------------------------------------------------------------------------
  subroutine FILE_CARTESC_write_1D( &
       var,                 &
       basename, title,     &
       varname, desc, unit, &
       dim_type, datatype,  &
       date, subsec,        &
       append, aggregate,   &
       standard_name,       &
       cell_measures        )
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
    character(len=*), intent(in), optional :: standard_name
    character(len=*), intent(in), optional :: cell_measures

    integer :: fid, vid
    !---------------------------------------------------------------------------

    LOG_INFO("FILE_CARTESC_write_1D",'(1x,2A)') 'Write to file (1D), name : ', trim(varname)

    call FILE_CARTESC_create( basename, title, datatype,                        & ! [IN]
                              fid,                                              & ! [OUT]
                              date=date, subsec=subsec,                         & ! [IN]
                              append=append, aggregate=aggregate, single=.true. ) ! [IN]

    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               standard_name=standard_name,                  & ! [IN]
                               cell_measures=cell_measures                   ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_1D( fid, vid, var, varname, dim_type )

    return
  end subroutine FILE_CARTESC_write_1D

  !-----------------------------------------------------------------------------
  !> Write 2D data to file
  subroutine FILE_CARTESC_write_2D( &
       var,                  &
       basename, title,      &
       varname, desc, unit,  &
       dim_type, datatype,   &
       date, subsec,         &
       fill_halo, haszcoord, &
       append, aggregate,    &
       standard_name,        &
       cell_measures         )
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
    logical,          intent(in), optional :: fill_halo !< switch whether include halo data or not    (default=false)
    logical,          intent(in), optional :: haszcoord !< switch whether include zcoordinate or not  (default=true)
    logical,          intent(in), optional :: append    !< switch whether append existing file or not (default=false)
    logical,          intent(in), optional :: aggregate
    character(len=*), intent(in), optional :: standard_name
    character(len=*), intent(in), optional :: cell_measures

    integer :: fid, vid
    !---------------------------------------------------------------------------

    LOG_INFO("FILE_CARTESC_write_2D",'(1x,2A)') 'Write to file (2D), name : ', trim(varname)

    call FILE_CARTESC_create( basename, title, datatype,         & ! [IN]
                              fid,                               & ! [OUT]
                              date=date, subsec=subsec,          & ! [IN]
                              haszcoord=haszcoord,               & ! [IN]
                              append=append, aggregate=aggregate ) ! [IN]

    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               standard_name=standard_name,                  & ! [IN]
                               cell_measures=cell_measures                   ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_2D( fid, vid, var, varname, dim_type, fill_halo )

    return
  end subroutine FILE_CARTESC_write_2D

  !-----------------------------------------------------------------------------
  !> Write 3D data to file
  subroutine FILE_CARTESC_write_3D( &
       var,                 &
       basename, title,     &
       varname, desc, unit, &
       dim_type, datatype,  &
       date, subsec,        &
       fill_halo,           &
       append, aggregate,   &
       standard_name,       &
       cell_measures        )
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
    logical,          intent(in), optional :: fill_halo !< include halo data?
    logical,          intent(in), optional :: append    !< append existing (closed) file?
    logical,          intent(in), optional :: aggregate
    character(len=*), intent(in), optional :: standard_name
    character(len=*), intent(in), optional :: cell_measures

    integer :: fid, vid
    !---------------------------------------------------------------------------

    LOG_INFO("FILE_CARTESC_write_3D",'(1x,2A)') 'Write to file (3D), name : ', trim(varname)

    call FILE_CARTESC_create( basename, title, datatype,         & ! [IN]
                              fid,                               & ! [OUT]
                              date=date, subsec=subsec,          & ! [IN]
                              append=append, aggregate=aggregate ) ! [IN]

    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               standard_name=standard_name,                  & ! [IN]
                               cell_measures=cell_measures                   ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_3D( fid, vid, var, varname, dim_type, fill_halo )


    return
  end subroutine FILE_CARTESC_write_3D

  !-----------------------------------------------------------------------------
  !> Write 3D data with time dimension to file
  subroutine FILE_CARTESC_write_3D_t( &
       var,                 &
       basename, title,     &
       varname, desc, unit, &
       dim_type, datatype,  &
       timeintv, tsince,    &
       timetarg, timeofs,   &
       fill_halo,           &
       append, aggregate,   &
       standard_name,       &
       cell_measures        )
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

    integer,          intent(in), optional :: timetarg  !< target timestep (optional)
    real(DP),         intent(in), optional :: timeofs   !< offset time     (optional)
    logical,          intent(in), optional :: fill_halo !< include halo data?
    logical,          intent(in), optional :: append    !< append existing (closed) file?
    logical,          intent(in), optional :: aggregate
    character(len=*), intent(in), optional :: standard_name
    character(len=*), intent(in), optional :: cell_measures

    integer  :: fid, vid
    integer  :: nsteps

    intrinsic :: size
    !---------------------------------------------------------------------------

    LOG_INFO("FILE_CARTESC_write_3D_t",'(1x,3A)') 'Write to file (3D), name : ', trim(varname), 'with time dimension'

    call FILE_CARTESC_create( basename, title, datatype,         & ! [IN]
                              fid,                               & ! [OUT]
                              date=tsince,                       & ! [IN]
                              append=append, aggregate=aggregate ) ! [IN]

    if ( present(timetarg) ) then
       nsteps = 1
    else
       nsteps = size(var,3)
    endif
    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               standard_name=standard_name,                  & ! [IN]
                               cell_measures=cell_measures,                  & ! [IN]
                               timeintv=timeintv, nsteps=nsteps              ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_3D_t( fid, vid, var, varname, dim_type, timeintv, &
                                      timetarg, timeofs, fill_halo                )

    return
  end subroutine FILE_CARTESC_write_3D_t

  !-----------------------------------------------------------------------------
  !> Write 4D data to file
  subroutine FILE_CARTESC_write_4D( &
       var,                 &
       basename, title,     &
       varname, desc, unit, &
       dim_type, datatype,  &
       timeintv, tsince,    &
       timetarg, timeofs,   &
       fill_halo,           &
       append, aggregate,   &
       standard_name,       &
       cell_measures        )
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

    integer,          intent(in), optional :: timetarg  !< target timestep (optional)
    real(DP),         intent(in), optional :: timeofs   !< offset time     (optional)
    logical,          intent(in), optional :: fill_halo !< include halo data?
    logical,          intent(in), optional :: append    !< append existing (closed) file?
    logical,          intent(in), optional :: aggregate
    character(len=*), intent(in), optional :: standard_name
    character(len=*), intent(in), optional :: cell_measures

    integer  :: fid, vid
    integer  :: nsteps

    intrinsic :: size
    !---------------------------------------------------------------------------

    LOG_INFO("FILE_CARTESC_write_4D",'(1x,2A)') 'Write to file (4D), name : ', trim(varname)

    call FILE_CARTESC_create( basename, title, datatype,         & ! [IN]
                              fid,                               & ! [OUT]
                              date=tsince,                       & ! [IN]
                              append=append, aggregate=aggregate ) ! [IN]

    if ( present(timetarg) ) then
       nsteps = 1
    else
       nsteps = size(var,3)
    endif
    call FILE_CARTESC_def_var( fid, varname, desc, unit, dim_type, datatype, & ! [IN]
                               vid,                                          & ! [OUT]
                               standard_name=standard_name,                  & ! [IN]
                               cell_measures=cell_measures,                  & ! [IN]
                               timeintv=timeintv, nsteps=nsteps              ) ! [IN]

    call FILE_CARTESC_enddef( fid )

    call FILE_CARTESC_write_var_4D( fid, vid, var, varname, dim_type, timeintv, &
                                    timetarg, timeofs, fill_halo                )

    return
  end subroutine FILE_CARTESC_write_4D

  !-----------------------------------------------------------------------------
  !> put global attributes
  subroutine FILE_CARTESC_put_globalAttributes( &
       fid, &
       prc_rank_x, prc_rank_y,         &
       prc_num_x, prc_num_y,           &
       prc_periodic_x, prc_periodic_y, &
       kmax, okmax, lkmax, ukmax,      &
       imaxg, jmaxg,                   &
       khalo, ihalo, jhalo,            &
       time, tunits, calendar          )
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_NAME
    use scale_file, only: &
       FILE_opened, &
       FILE_Set_Attribute

    integer,          intent(in) :: fid
    integer,          intent(in) :: prc_rank_x, prc_rank_y
    integer,          intent(in) :: prc_num_x, prc_num_y
    logical,          intent(in) :: prc_periodic_x, prc_periodic_y
    integer,          intent(in) :: kmax, okmax, lkmax, ukmax
    integer,          intent(in) :: imaxg, jmaxg
    integer,          intent(in) :: khalo, ihalo, jhalo
    real(DP),         intent(in) :: time
    character(len=*), intent(in) :: tunits
    character(len=*), intent(in) :: calendar
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    if ( .not. prof ) call PROF_rapstart('FILE_O_NetCDF', 2)

    call FILE_Set_Attribute( fid, "global", "Conventions", "CF-1.6" ) ! [IN]

    call FILE_Set_Attribute( fid, "global", "grid_name", ATMOS_GRID_CARTESC_NAME ) ! [IN]

    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_rank_x", (/prc_rank_x/) ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_rank_y", (/prc_rank_y/) ) ! [IN]

    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_num_x",  (/prc_num_x/) ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_num_y",  (/prc_num_y/) ) ! [IN]

    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_periodic_z", .false.        ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_periodic_x", prc_periodic_x ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_cartesC_prc_periodic_y", prc_periodic_y ) ! [IN]

    call FILE_Set_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_imaxg", (/imaxg/) ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_jmaxg", (/jmaxg/) ) ! [IN]

                     call FILE_Set_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_kmax", (/kmax/)  ) ! [IN]
    if ( okmax > 0 ) call FILE_Set_Attribute( fid, "global", "scale_ocean_grid_cartesC_index_kmax", (/okmax/) ) ! [IN]
    if ( lkmax > 0 ) call FILE_Set_Attribute( fid, "global", "scale_land_grid_cartesC_index_kmax",  (/lkmax/) ) ! [IN]
    if ( ukmax > 0 ) call FILE_Set_Attribute( fid, "global", "scale_urban_grid_cartesC_index_kmax", (/ukmax/) ) ! [IN]

    call FILE_Set_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_khalo", (/khalo/) ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_ihalo", (/ihalo/) ) ! [IN]
    call FILE_Set_Attribute( fid, "global", "scale_atmos_grid_cartesC_index_jhalo", (/jhalo/) ) ! [IN]

    if ( calendar /= "" ) call FILE_Set_Attribute( fid, "global", "calendar", calendar )
    call FILE_Set_Attribute( fid, "global", "time_units", tunits )
    call FILE_Set_Attribute( fid, "global", "time_start", (/time/) )

    if ( .not. prof ) call PROF_rapend('FILE_O_NetCDF', 2)

    return
  end subroutine FILE_CARTESC_put_globalAttributes

  !-----------------------------------------------------------------------------
  !> define axis variables in the file
  subroutine FILE_CARTESC_def_axes( &
       fid,   &
       dtype, &
       hasZ   )
    use scale_file, only: &
       FILE_opened,                   &
       FILE_get_AGGREGATE,            &
       FILE_Def_Axis,                 &
       FILE_Set_Attribute,            &
       FILE_Def_AssociatedCoordinate, &
       FILE_Add_AssociatedVariable
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_mapprojection, only: &
       MAPPROJECTION_mappinginfo
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: dtype
    logical, intent(in) :: hasZ

    integer :: iall  ! grid size, x-axis, (whole domain or local tile), including halo
    integer :: jall  ! grid size, y-axis, (whole domain or local tile), including halo
    integer :: isize ! grid size, x-axis, (whole domain or local tile), without halo except domain edge
    integer :: jsize ! grid size, y-axis, (whole domain or local tile), without halo except domain edge

    character(len=2) :: axisname(3)

    logical, save :: set_dim = .false.
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

    if ( .not. set_dim ) then
       call set_dimension_informations
       set_dim = .true.
    end if

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

    if ( hasZ ) then
       call FILE_Def_Axis( fid, 'z'  , 'Z'              , 'm', 'z'  , dtype, KMAX,   bounds=.true. )
       call FILE_Def_Axis( fid, 'zh' , 'Z (half level)' , 'm', 'zh' , dtype, KMAX+1, bounds=.true. )

       if ( OKMAX > 0 ) then
          call FILE_Def_Axis( fid, 'oz' , 'OZ'             , 'm', 'oz' , dtype, OKMAX,   bounds=.true. )
          call FILE_Def_Axis( fid, 'ozh', 'OZ (half level)', 'm', 'ozh', dtype, OKMAX+1, bounds=.true. )
       end if

       if ( LKMAX > 0 ) then
          call FILE_Def_Axis( fid, 'lz' , 'LZ'             , 'm', 'lz' , dtype, LKMAX,   bounds=.true. )
          call FILE_Def_Axis( fid, 'lzh', 'LZ (half level)', 'm', 'lzh', dtype, LKMAX+1, bounds=.true. )
       end if

       if ( UKMAX > 0 ) then
          call FILE_Def_Axis( fid, 'uz' , 'UZ'             , 'm', 'uz' , dtype, UKMAX,   bounds=.true. )
          call FILE_Def_Axis( fid, 'uzh', 'UZ (half level)', 'm', 'uzh', dtype, UKMAX+1, bounds=.true. )
       end if
    end if

    call FILE_Def_Axis( fid, 'x'  , 'X'              , 'm', 'x'  , dtype, isize, bounds=.true. )
    call FILE_Def_Axis( fid, 'xh' , 'X (half level)' , 'm', 'xh' , dtype, isize, bounds=.true. )
    call FILE_Def_Axis( fid, 'y'  , 'Y'              , 'm', 'y'  , dtype, jsize, bounds=.true. )
    call FILE_Def_Axis( fid, 'yh' , 'Y (half level)' , 'm', 'yh' , dtype, jsize, bounds=.true. )

    if ( hasZ ) then
       call FILE_Def_Axis( fid, 'CZ'   , 'Atmos Grid Center Position Z',      'm', 'CZ',   dtype, KA      )
       call FILE_Def_Axis( fid, 'FZ'   , 'Atmos Grid Face   Position Z',      'm', 'FZ',   dtype, KA+1    )
       call FILE_Def_Axis( fid, 'CDZ'  , 'Grid Cell length Z',                'm', 'CZ',   dtype, KA      )
       call FILE_Def_Axis( fid, 'FDZ'  , 'Grid distance Z',                   'm', 'FDZ',  dtype, KA-1    )
       call FILE_Def_Axis( fid, 'CBFZ' , 'Boundary factor Center Z',          '1', 'CZ',   dtype, KA      )
       call FILE_Def_Axis( fid, 'FBFZ' , 'Boundary factor Face Z',            '1', 'FZ',   dtype, KA+1    )

       if ( OKMAX > 0 ) then
          call FILE_Def_Axis( fid, 'OCZ'  , 'Ocean Grid Center Position Z',      'm', 'OCZ',  dtype, OKMAX   )
          call FILE_Def_Axis( fid, 'OFZ'  , 'Ocean Grid Face   Position Z',      'm', 'OFZ',  dtype, OKMAX+1 )
          call FILE_Def_Axis( fid, 'OCDZ' , 'Ocean Grid Cell length Z',          'm', 'OCZ',  dtype, OKMAX   )
       end if

       if ( LKMAX > 0 ) then
          call FILE_Def_Axis( fid, 'LCZ'  , 'Land Grid Center Position Z',       'm', 'LCZ',  dtype, LKMAX   )
          call FILE_Def_Axis( fid, 'LFZ'  , 'Land Grid Face   Position Z',       'm', 'LFZ',  dtype, LKMAX+1 )
          call FILE_Def_Axis( fid, 'LCDZ' , 'Land Grid Cell length Z',           'm', 'LCZ',  dtype, LKMAX   )
       end if

       if ( UKMAX > 0 ) then
          call FILE_Def_Axis( fid, 'UCZ'  , 'Urban Grid Center Position Z',      'm', 'UCZ',  dtype, UKMAX   )
          call FILE_Def_Axis( fid, 'UFZ'  , 'Urban Grid Face   Position Z',      'm', 'UFZ',  dtype, UKMAX+1 )
          call FILE_Def_Axis( fid, 'UCDZ' , 'Urban Grid Cell length Z',          'm', 'UCZ',  dtype, UKMAX   )
       end if
    end if

    call FILE_Def_Axis( fid, 'CX'   , 'Atmos Grid Center Position X',      'm', 'CX',  dtype, iall   )
    call FILE_Def_Axis( fid, 'CY'   , 'Atmos Grid Center Position Y',      'm', 'CY',  dtype, jall   )
    call FILE_Def_Axis( fid, 'FX'   , 'Atmos Grid Face   Position X',      'm', 'FX',  dtype, iall+1 )
    call FILE_Def_Axis( fid, 'FY'   , 'Atmos Grid Face   Position Y',      'm', 'FY',  dtype, jall+1 )
    call FILE_Def_Axis( fid, 'CDX'  , 'Grid Cell length X',                'm', 'CX',  dtype, iall   )
    call FILE_Def_Axis( fid, 'CDY'  , 'Grid Cell length Y',                'm', 'CY',  dtype, jall   )
    call FILE_Def_Axis( fid, 'FDX'  , 'Grid distance X',                   'm', 'FX',  dtype, iall+1 )
    call FILE_Def_Axis( fid, 'FDY'  , 'Grid distance Y',                   'm', 'FY',  dtype, jall+1 )
    call FILE_Def_Axis( fid, 'CBFX' , 'Boundary factor Center X',          '1', 'CX',  dtype, iall   )
    call FILE_Def_Axis( fid, 'CBFY' , 'Boundary factor Center Y',          '1', 'CY',  dtype, jall   )
    call FILE_Def_Axis( fid, 'FBFX' , 'Boundary factor Face X',            '1', 'FX',  dtype, iall+1 )
    call FILE_Def_Axis( fid, 'FBFY' , 'Boundary factor Face Y',            '1', 'FY',  dtype, jall+1 )

    call FILE_Def_Axis( fid, 'CXG'  , 'Grid Center Position X (global)',   'm', 'CXG', dtype, IAG   )
    call FILE_Def_Axis( fid, 'CYG'  , 'Grid Center Position Y (global)',   'm', 'CYG', dtype, JAG   )
    call FILE_Def_Axis( fid, 'FXG'  , 'Grid Face   Position X (global)',   'm', 'FXG', dtype, IAG+1 )
    call FILE_Def_Axis( fid, 'FYG'  , 'Grid Face   Position Y (global)',   'm', 'FYG', dtype, JAG+1 )
    call FILE_Def_Axis( fid, 'CDXG' , 'Grid Cell length X (global)',       'm', 'CXG', dtype, IAG   )
    call FILE_Def_Axis( fid, 'CDYG' , 'Grid Cell length Y (global)',       'm', 'CYG', dtype, JAG   )
    call FILE_Def_Axis( fid, 'FDXG' , 'Grid distance X (global)',          'm', 'FXG', dtype, IAG+1 )
    call FILE_Def_Axis( fid, 'FDYG' , 'Grid distance Y (global)',          'm', 'FYG', dtype, JAG+1 )
    call FILE_Def_Axis( fid, 'CBFXG', 'Boundary factor Center X (global)', '1', 'CXG', dtype, IAG   )
    call FILE_Def_Axis( fid, 'CBFYG', 'Boundary factor Center Y (global)', '1', 'CYG', dtype, JAG   )
    call FILE_Def_Axis( fid, 'FBFXG', 'Boundary factor Face   X (global)', '1', 'FXG', dtype, IAG+1 )
    call FILE_Def_Axis( fid, 'FBFYG', 'Boundary factor Face   Y (global)', '1', 'FYG', dtype, JAG+1 )

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
       axisname(1:2) = (/'x ','y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'topo' ,  'topography',                 'm',            axisname(1:2), dtype )
    end if
    axisname(1:2) = (/'x ','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'lsmask', 'fraction for land-sea mask', '1',            axisname(1:2), dtype )

    axisname(1:2) = (/'x ','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'cell_area',    'area of grid cell',                  'm2', axisname(1:2), dtype )
    axisname(1:2) = (/'xh','y '/)
    call FILE_Def_AssociatedCoordinate( fid, 'cell_area_uy', 'area of grid cell (half level uy)',  'm2', axisname(1:2), dtype )
    axisname(1:2) = (/'x ','yh'/)
    call FILE_Def_AssociatedCoordinate( fid, 'cell_area_xv', 'area of grid cell (half level xv)',  'm2', axisname(1:2), dtype )
    axisname(1:2) = (/'xh','yh'/)

    if ( hasZ ) then
       axisname = (/'z ', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'height',     'height above ground level', &
                                           'm', axisname(1:3), dtype                        )
       axisname = (/'zh', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'height_wxy', 'height above ground level (half level wxy)', &
                                           'm', axisname(1:3), dtype                                         )

       axisname = (/'z ', 'xh', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_zuy_x', 'area of grid cell face (half level zuy, normal x)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'z ', 'x ', 'yh'/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_zxv_y', 'area of grid cell face (half level zxv, normal y)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'zh', 'xh', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_wuy_x', 'area of grid cell face (half level wuy, normal x)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'zh', 'x ', 'yh'/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_wxv_y', 'area of grid cell face (half level wxv, normal y)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'z ', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_zxy_x', 'area of grid cell face (half level zxy, normal x)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'z ', 'xh', 'yh'/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_zuv_y', 'area of grid cell face (half level zuv, normal y)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'z ', 'xh', 'yh'/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_zuv_x', 'area of grid cell face (half level zuv, normal x)', &
                                           'm2', axisname(1:3), dtype                       )
       axisname = (/'z ', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_area_zxy_y', 'area of grid cell face (half level zxy, normal y)', &
                                           'm2', axisname(1:3), dtype                       )

       axisname = (/'z ', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_volume',     'volume of grid cell', &
                                           'm3', axisname(1:3), dtype                       )
       axisname = (/'zh', 'x ', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_volume_wxy', 'volume of grid cell (half level wxy)', &
                                           'm3', axisname(1:3), dtype                                       )
       axisname = (/'z ', 'xh', 'y '/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_volume_zuy', 'volume of grid cell (half level zuy)', &
                                           'm3', axisname(1:3), dtype                                       )
       axisname = (/'z ', 'x ', 'yh'/)
       call FILE_Def_AssociatedCoordinate( fid, 'cell_volume_zxv', 'volume of grid cell (half level zxv)', &
                                           'm3', axisname(1:3), dtype                                       )

       if ( OKMAX > 0 ) then
          axisname = (/'oz', 'x ', 'y '/)
          call FILE_Def_AssociatedCoordinate( fid, 'cell_volume_oxy', 'volume of grid cell', &
                                              'm3', axisname(1:3), dtype                      )
       end if
       if ( LKMAX > 0 ) then
          axisname = (/'lz', 'x ', 'y '/)
          call FILE_Def_AssociatedCoordinate( fid, 'cell_volume_lxy', 'volume of grid cell', &
                                              'm3', axisname(1:3), dtype                      )
       end if
       if ( UKMAX > 0 ) then
          axisname = (/'uz', 'x ', 'y '/)
          call FILE_Def_AssociatedCoordinate( fid, 'cell_volume_uxy', 'volume of grid cell', &
                                              'm3', axisname(1:3), dtype                      )
       end if
    endif

    ! attributes

    if ( hasZ ) then
       if ( OKMAX > 0 ) then
          call FILE_Set_Attribute( fid, 'oz' , 'positive', 'down' )
          call FILE_Set_Attribute( fid, 'ozh', 'positive', 'down' )
       end if
       if ( LKMAX > 0 ) then
          call FILE_Set_Attribute( fid, 'lz' , 'positive', 'down' )
          call FILE_Set_Attribute( fid, 'lzh', 'positive', 'down' )
       end if
       if ( UKMAX > 0 ) then
          call FILE_Set_Attribute( fid, 'uz' , 'positive', 'down' )
          call FILE_Set_Attribute( fid, 'uzh', 'positive', 'down' )
       end if
       if ( OKMAX > 0 ) then
          call FILE_Set_Attribute( fid, 'OCZ', 'positive', 'down' )
          call FILE_Set_Attribute( fid, 'OFZ', 'positive', 'down' )
       end if
       if ( LKMAX > 0 ) then
          call FILE_Set_Attribute( fid, 'LCZ', 'positive', 'down' )
          call FILE_Set_Attribute( fid, 'LFZ', 'positive', 'down' )
       end if
       if ( UKMAX > 0 ) then
          call FILE_Set_Attribute( fid, 'UCZ', 'positive', 'down' )
          call FILE_Set_Attribute( fid, 'UFZ', 'positive', 'down' )
       end if
    endif

    if ( FILE_get_AGGREGATE(fid) ) then
       call FILE_Set_Attribute( fid, "x" , "size_global" , (/ IAG /) )
       call FILE_Set_Attribute( fid, "x" , "start_global", (/ 1   /) )
       call FILE_Set_Attribute( fid, "x" , "halo_global" , (/ IHALO, IHALO /) )
       call FILE_Set_Attribute( fid, "x" , "halo_local"  , (/ IHALO, IHALO /) )

       call FILE_Set_Attribute( fid, "xh", "size_global" , (/ IAG+1 /) )
       call FILE_Set_Attribute( fid, "xh", "start_global", (/ 1     /) )
       call FILE_Set_Attribute( fid, "xh", "halo_global" , (/ IHALO, IHALO /) )
       call FILE_Set_Attribute( fid, "xh", "halo_local"  , (/ IHALO, IHALO /) )

       call FILE_Set_Attribute( fid, "y" , "size_global" , (/ JAG /) )
       call FILE_Set_Attribute( fid, "y" , "start_global", (/ 1   /) )
       call FILE_Set_Attribute( fid, "y" , "halo_global" , (/ JHALO, JHALO /) )
       call FILE_Set_Attribute( fid, "y" , "halo_local"  , (/ JHALO, JHALO /) )

       call FILE_Set_Attribute( fid, "yh", "size_global" , (/ JAG+1 /) )
       call FILE_Set_Attribute( fid, "yh", "start_global", (/ 1     /) )
       call FILE_Set_Attribute( fid, "yh", "halo_global" , (/ JHALO, JHALO /) )
       call FILE_Set_Attribute( fid, "yh", "halo_local"  , (/ JHALO, JHALO /) )
    else
       call FILE_Set_Attribute( fid, "x" , "size_global" , FILE_CARTESC_AXIS_info(1)%size_global (:) )
       call FILE_Set_Attribute( fid, "x" , "start_global", FILE_CARTESC_AXIS_info(1)%start_global(:) )
       call FILE_Set_Attribute( fid, "x" , "halo_global" , FILE_CARTESC_AXIS_info(1)%halo_global (:) )
       call FILE_Set_Attribute( fid, "x" , "halo_local"  , FILE_CARTESC_AXIS_info(1)%halo_local  (:) )

       call FILE_Set_Attribute( fid, "xh", "size_global" , FILE_CARTESC_AXIS_info(2)%size_global (:) )
       call FILE_Set_Attribute( fid, "xh", "start_global", FILE_CARTESC_AXIS_info(2)%start_global(:) )
       call FILE_Set_Attribute( fid, "xh", "halo_global" , FILE_CARTESC_AXIS_info(2)%halo_global (:) )
       call FILE_Set_Attribute( fid, "xh", "halo_local"  , FILE_CARTESC_AXIS_info(2)%halo_local  (:) )

       call FILE_Set_Attribute( fid, "y" , "size_global" , FILE_CARTESC_AXIS_info(3)%size_global (:) )
       call FILE_Set_Attribute( fid, "y" , "start_global", FILE_CARTESC_AXIS_info(3)%start_global(:) )
       call FILE_Set_Attribute( fid, "y" , "halo_global" , FILE_CARTESC_AXIS_info(3)%halo_global (:) )
       call FILE_Set_Attribute( fid, "y" , "halo_local"  , FILE_CARTESC_AXIS_info(3)%halo_local  (:) )

       call FILE_Set_Attribute( fid, "yh", "size_global" , FILE_CARTESC_AXIS_info(4)%size_global (:) )
       call FILE_Set_Attribute( fid, "yh", "start_global", FILE_CARTESC_AXIS_info(4)%start_global(:) )
       call FILE_Set_Attribute( fid, "yh", "halo_global" , FILE_CARTESC_AXIS_info(4)%halo_global (:) )
       call FILE_Set_Attribute( fid, "yh", "halo_local"  , FILE_CARTESC_AXIS_info(4)%halo_local  (:) )
    end if

    call FILE_Set_Attribute( fid, "x" , "periodic"    , FILE_CARTESC_AXIS_info(1)%periodic        )
    call FILE_Set_Attribute( fid, "xh", "periodic"    , FILE_CARTESC_AXIS_info(2)%periodic        )
    call FILE_Set_Attribute( fid, "y" , "periodic"    , FILE_CARTESC_AXIS_info(3)%periodic        )
    call FILE_Set_Attribute( fid, "yh", "periodic"    , FILE_CARTESC_AXIS_info(4)%periodic        )

    ! map projection info

    if ( MAPPROJECTION_mappinginfo%mapping_name /= "" ) then
       call FILE_Set_Attribute( fid, "x" , "standard_name", "projection_x_coordinate" )
       call FILE_Set_Attribute( fid, "xh", "standard_name", "projection_x_coordinate" )
       call FILE_Set_Attribute( fid, "y" , "standard_name", "projection_y_coordinate" )
       call FILE_Set_Attribute( fid, "yh", "standard_name", "projection_y_coordinate" )

       call FILE_Add_AssociatedVariable( fid, MAPPROJECTION_mappinginfo%mapping_name )
       call FILE_Set_Attribute( fid, MAPPROJECTION_mappinginfo%mapping_name, "grid_mapping_name",  MAPPROJECTION_mappinginfo%mapping_name )

       if ( MAPPROJECTION_mappinginfo%false_easting /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                    & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name, & ! [IN]
                                   "false_easting",                        & ! [IN]
                                   MAPPROJECTION_mappinginfo%false_easting ) ! [IN]
       endif

       if ( MAPPROJECTION_mappinginfo%false_northing /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                     & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name,  & ! [IN]
                                   "false_northing",                        & ! [IN]
                                   MAPPROJECTION_mappinginfo%false_northing ) ! [IN]
       endif

       if ( MAPPROJECTION_mappinginfo%longitude_of_central_meridian /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                                    & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name,                 & ! [IN]
                                   "longitude_of_central_meridian",                        & ! [IN]
                                   MAPPROJECTION_mappinginfo%longitude_of_central_meridian ) ! [IN]
       endif

       if ( MAPPROJECTION_mappinginfo%longitude_of_projection_origin /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                                     & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name,                  & ! [IN]
                                   "longitude_of_projection_origin",                        & ! [IN]
                                   MAPPROJECTION_mappinginfo%longitude_of_projection_origin ) ! [IN]
       endif

       if ( MAPPROJECTION_mappinginfo%latitude_of_projection_origin /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                                    & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name,                 & ! [IN]
                                   "latitude_of_projection_origin",                        & ! [IN]
                                   MAPPROJECTION_mappinginfo%latitude_of_projection_origin ) ! [IN]
       endif

       if ( MAPPROJECTION_mappinginfo%straight_vertical_longitude_from_pole /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                                            & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name,                         & ! [IN]
                                   "straight_vertical_longitude_from_pole",                        & ! [IN]
                                   MAPPROJECTION_mappinginfo%straight_vertical_longitude_from_pole ) ! [IN]
       endif

       if ( MAPPROJECTION_mappinginfo%standard_parallel(1) /= UNDEF ) then
          if ( MAPPROJECTION_mappinginfo%standard_parallel(2) /= UNDEF ) then
             call FILE_Set_Attribute( fid,                                           & ! [IN]
                                      MAPPROJECTION_mappinginfo%mapping_name,        & ! [IN]
                                      "standard_parallel",                           & ! [IN]
                                      MAPPROJECTION_mappinginfo%standard_parallel(:) ) ! [IN]
          else
             call FILE_Set_Attribute( fid,                                           & ! [IN]
                                      MAPPROJECTION_mappinginfo%mapping_name,        & ! [IN]
                                      "standard_parallel",                           & ! [IN]
                                      MAPPROJECTION_mappinginfo%standard_parallel(1) ) ! [IN]
          endif
       endif

       if ( MAPPROJECTION_mappinginfo%rotation /= UNDEF ) then
          call FILE_Set_Attribute( fid,                                    & ! [IN]
                                   MAPPROJECTION_mappinginfo%mapping_name, & ! [IN]
                                   "rotation",                             & ! [IN]
                                   MAPPROJECTION_mappinginfo%rotation      ) ! [IN]
       endif

    endif

    ! cell measures

    call FILE_Set_Attribute( fid, "cell_area",    "standard_name", "area" ) ! [IN]
    call FILE_Set_Attribute( fid, "cell_area_uy", "standard_name", "area" ) ! [IN]
    call FILE_Set_Attribute( fid, "cell_area_xv", "standard_name", "area" ) ! [IN]

    if ( hasZ ) then
       call FILE_Set_Attribute( fid, "cell_area_zuy_x", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_zxv_y", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_wuy_x", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_wxv_y", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_zxy_x", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_zuv_y", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_zuv_x", "standard_name", "area" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_area_zxy_y", "standard_name", "area" ) ! [IN]

       call FILE_Set_Attribute( fid, "cell_volume",     "standard_name", "volume" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_volume_wxy", "standard_name", "volume" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_volume_zuy", "standard_name", "volume" ) ! [IN]
       call FILE_Set_Attribute( fid, "cell_volume_zxv", "standard_name", "volume" ) ! [IN]

       if ( OKMAX > 0 ) then
          call FILE_Set_Attribute( fid, "cell_volume_oxy", "standard_name", "volume" ) ! [IN]
       end if
       if ( LKMAX > 0 ) then
          call FILE_Set_Attribute( fid, "cell_volume_lxy", "standard_name", "volume" ) ! [IN]
       end if
       if ( UKMAX > 0 ) then
          call FILE_Set_Attribute( fid, "cell_volume_uxy", "standard_name", "volume" ) ! [IN]
       end if
    end if

    ! SGRID
    call FILE_Add_AssociatedVariable( fid, "grid" )
    call FILE_Set_Attribute( fid, "grid", "cf_role",             "grid_topology" )
    call FILE_Set_Attribute( fid, "grid", "topology_dimension",  (/ 2 /) )
    call FILE_Set_Attribute( fid, "grid", "node_dimensions",     "xh yh" )
    call FILE_Set_Attribute( fid, "grid", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_Set_Attribute( fid, "grid", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_Set_Attribute( fid, "grid", "face_coordinates",    "lon lat" )
    call FILE_Set_Attribute( fid, "grid", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_Set_Attribute( fid, "grid", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_Set_Attribute( fid, "grid", "vertical_dimensions", "z: zh (padding: none)" )

    if ( OKMAX > 0 ) then
       call FILE_Add_AssociatedVariable( fid, "grid_ocean" )
       call FILE_Set_Attribute( fid, "grid_ocean", "cf_role",             "grid_topology" )
       call FILE_Set_Attribute( fid, "grid_ocean", "topology_dimension",  (/ 2 /) )
       call FILE_Set_Attribute( fid, "grid_ocean", "node_dimensions",     "xh yh" )
       call FILE_Set_Attribute( fid, "grid_ocean", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
       call FILE_Set_Attribute( fid, "grid_ocean", "node_coordinates",    "lon_uv lat_uv" )
       call FILE_Set_Attribute( fid, "grid_ocean", "face_coordinates",    "lon lat" )
       call FILE_Set_Attribute( fid, "grid_ocean", "edge1_coordinates",   "lon_uy lat_uy" )
       call FILE_Set_Attribute( fid, "grid_ocean", "edge2_coordinates",   "lon_xv lat_xv" )
       call FILE_Set_Attribute( fid, "grid_ocean", "vertical_dimensions", "oz: ozh (padding: none)" )
    end if

    if ( LKMAX > 0 ) then
       call FILE_Add_AssociatedVariable( fid, "grid_land" )
       call FILE_Set_Attribute( fid, "grid_land", "cf_role",             "grid_topology" )
       call FILE_Set_Attribute( fid, "grid_land", "topology_dimension",  (/ 2 /) )
       call FILE_Set_Attribute( fid, "grid_land", "node_dimensions",     "xh yh" )
       call FILE_Set_Attribute( fid, "grid_land", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
       call FILE_Set_Attribute( fid, "grid_land", "node_coordinates",    "lon_uv lat_uv" )
       call FILE_Set_Attribute( fid, "grid_land", "face_coordinates",    "lon lat" )
       call FILE_Set_Attribute( fid, "grid_land", "edge1_coordinates",   "lon_uy lat_uy" )
       call FILE_Set_Attribute( fid, "grid_land", "edge2_coordinates",   "lon_xv lat_xv" )
       call FILE_Set_Attribute( fid, "grid_land", "vertical_dimensions", "lz: lzh (padding: none)" )
    end if

    if ( UKMAX > 0 ) then
       call FILE_Add_AssociatedVariable( fid, "grid_urban" )
       call FILE_Set_Attribute( fid, "grid_urban", "cf_role",             "grid_topology" )
       call FILE_Set_Attribute( fid, "grid_urban", "topology_dimension",  (/ 2 /) )
       call FILE_Set_Attribute( fid, "grid_urban", "node_dimensions",     "xh yh" )
       call FILE_Set_Attribute( fid, "grid_urban", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
       call FILE_Set_Attribute( fid, "grid_urban", "node_coordinates",    "lon_uv lat_uv" )
       call FILE_Set_Attribute( fid, "grid_urban", "face_coordinates",    "lon lat" )
       call FILE_Set_Attribute( fid, "grid_urban", "edge1_coordinates",   "lon_uy lat_uy" )
       call FILE_Set_Attribute( fid, "grid_urban", "edge2_coordinates",   "lon_xv lat_xv" )
       call FILE_Set_Attribute( fid, "grid_urban", "vertical_dimensions", "uz: uzh (padding: none)" )
    end if

    call FILE_Add_AssociatedVariable( fid, "grid_model" )
    call FILE_Set_Attribute( fid, "grid_model", "cf_role",             "grid_topology" )
    call FILE_Set_Attribute( fid, "grid_model", "topology_dimension",  (/ 2 /) )
    call FILE_Set_Attribute( fid, "grid_model", "node_dimensions",     "FX FY" )
    call FILE_Set_Attribute( fid, "grid_model", "face_dimensions",     "CX: FY (padding: none) CY: FY (padding: none)" )
    call FILE_Set_Attribute( fid, "grid_model", "vertical_dimensions", "CZ: FZ (padding: none)" )

    call FILE_Add_AssociatedVariable( fid, "grid_model_global" )
    call FILE_Set_Attribute( fid, "grid_model_global", "cf_role",             "grid_topology" )
    call FILE_Set_Attribute( fid, "grid_model_global", "topology_dimension",  (/ 2 /) )
    call FILE_Set_Attribute( fid, "grid_model_global", "node_dimensions",     "FXG FYG" )
    call FILE_Set_Attribute( fid, "grid_model_global", "face_dimensions",     "CXG: FYG (padding: none) CYG: FYG (padding: none)" )
    call FILE_Set_Attribute( fid, "grid_model_global", "vertical_dimensions", "CZ: FZ (padding: none)" )

    return
  end subroutine FILE_CARTESC_def_axes

  !-----------------------------------------------------------------------------
  !> write axis to the file
  subroutine FILE_CARTESC_write_axes( &
       fid,       &
       haszcoord, &
       start      )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_file, only: &
       FILE_opened, &
       FILE_get_AGGREGATE, &
       FILE_Write_Axis,                 &
       FILE_Write_AssociatedCoordinate
    use scale_prc, only: &
       PRC_myrank, &
       PRC_IsMaster
    use scale_prc_cartesC, only: &
       PRC_2Drank
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ,    &
       ATMOS_GRID_CARTESC_CX,    &
       ATMOS_GRID_CARTESC_CY,    &
       ATMOS_GRID_CARTESC_FZ,    &
       ATMOS_GRID_CARTESC_FX,    &
       ATMOS_GRID_CARTESC_FY,    &
       ATMOS_GRID_CARTESC_CDZ,   &
       ATMOS_GRID_CARTESC_CDX,   &
       ATMOS_GRID_CARTESC_CDY,   &
       ATMOS_GRID_CARTESC_FDZ,   &
       ATMOS_GRID_CARTESC_FDX,   &
       ATMOS_GRID_CARTESC_FDY,   &
       ATMOS_GRID_CARTESC_CBFZ,  &
       ATMOS_GRID_CARTESC_CBFX,  &
       ATMOS_GRID_CARTESC_CBFY,  &
       ATMOS_GRID_CARTESC_FBFZ,  &
       ATMOS_GRID_CARTESC_FBFX,  &
       ATMOS_GRID_CARTESC_FBFY,  &
       ATMOS_GRID_CARTESC_CXG,   &
       ATMOS_GRID_CARTESC_CYG,   &
       ATMOS_GRID_CARTESC_FXG,   &
       ATMOS_GRID_CARTESC_FYG,   &
       ATMOS_GRID_CARTESC_CDXG,  &
       ATMOS_GRID_CARTESC_CDYG,  &
       ATMOS_GRID_CARTESC_FDXG,  &
       ATMOS_GRID_CARTESC_FDYG,  &
       ATMOS_GRID_CARTESC_CBFXG, &
       ATMOS_GRID_CARTESC_CBFYG, &
       ATMOS_GRID_CARTESC_FBFXG, &
       ATMOS_GRID_CARTESC_FBFYG
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_CZ,  &
       OCEAN_GRID_CARTESC_FZ,  &
       OCEAN_GRID_CARTESC_CDZ
    use scale_land_grid_cartesC, only: &
       LAND_GRID_CARTESC_CZ,  &
       LAND_GRID_CARTESC_FZ,  &
       LAND_GRID_CARTESC_CDZ
    use scale_urban_grid_cartesC, only: &
       URBAN_GRID_CARTESC_CZ,  &
       URBAN_GRID_CARTESC_FZ,  &
       URBAN_GRID_CARTESC_CDZ
    implicit none

    integer, intent(in)  :: fid
    logical, intent(in)  :: haszcoord
    integer, intent(in)  :: start(3)

    logical :: put_z, put_x, put_y
    integer :: XSB, XEB, YSB, YEB

    real(RP) :: z_bnds(2,KA), zh_bnds(2,0:KA)
    real(RP) :: oz_bnds(2,OKA), ozh_bnds(2,0:OKA)
    real(RP) :: lz_bnds(2,LKA), lzh_bnds(2,0:LKA)
    real(RP) :: uz_bnds(2,UKA), uzh_bnds(2,0:UKA)
    real(RP) :: x_bnds(2,IA), xh_bnds(2,0:IA)
    real(RP) :: y_bnds(2,JA), yh_bnds(2,0:JA)

    real(RP) :: FDXG(0:IAG), FDYG(0:JAG)
    real(RP) :: FDX(0:IA), FDY(0:JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

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
       ! atmos
       call FILE_Write_Axis( fid, 'z', ATMOS_GRID_CARTESC_CZ(KS:KE), start(1:1) )
       do k = KS, KE
          z_bnds(1,k) = ATMOS_GRID_CARTESC_FZ(k-1)
          z_bnds(2,k) = ATMOS_GRID_CARTESC_FZ(k  )
       end do
       call FILE_Write_AssociatedCoordinate( fid, 'z_bnds', z_bnds(:,KS:KE), start(1:1) )

       call FILE_Write_Axis( fid, 'zh', ATMOS_GRID_CARTESC_FZ(KS-1:KE)  , start(1:1) )
       do k = KS-1, KE
          zh_bnds(1,k) = ATMOS_GRID_CARTESC_CZ(k  )
          zh_bnds(2,k) = ATMOS_GRID_CARTESC_CZ(k+1)
       end do
       call FILE_Write_AssociatedCoordinate( fid, 'zh_bnds', zh_bnds(:,KS-1:KE), start(1:1) )

       ! ocean
       if ( OKMAX > 0 ) then
          call FILE_Write_Axis( fid, 'oz', OCEAN_GRID_CARTESC_CZ(OKS:OKE), start(1:1) )
          do k = OKS, OKE
             oz_bnds(1,k) = OCEAN_GRID_CARTESC_FZ(k-1)
             oz_bnds(2,k) = OCEAN_GRID_CARTESC_FZ(k  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'oz_bnds', oz_bnds(:,OKS:OKE), start(1:1) )

          call FILE_Write_Axis( fid, 'ozh', OCEAN_GRID_CARTESC_FZ(OKS-1:OKE), start(1:1) )
          ozh_bnds(1,OKS-1) = OCEAN_GRID_CARTESC_FZ(OKS-1)
          do k = OKS-1, OKE-1
             ozh_bnds(2,k  ) = OCEAN_GRID_CARTESC_CZ(k+1)
             ozh_bnds(1,k+1) = OCEAN_GRID_CARTESC_CZ(k+1)
          end do
          ozh_bnds(2,OKE) = OCEAN_GRID_CARTESC_FZ(OKE)
          call FILE_Write_AssociatedCoordinate( fid, 'ozh_bnds', ozh_bnds(:,OKS-1:OKE), start(1:1) )
       end if

       ! land
       if ( LKMAX > 0 ) then
          call FILE_Write_Axis( fid, 'lz', LAND_GRID_CARTESC_CZ(LKS:LKE), start(1:1) )
          do k = LKS, LKE
             lz_bnds(1,k) = LAND_GRID_CARTESC_FZ(k-1)
             lz_bnds(2,k) = LAND_GRID_CARTESC_FZ(k  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'lz_bnds', lz_bnds(:,LKS:LKE), start(1:1) )

          call FILE_Write_Axis( fid, 'lzh', LAND_GRID_CARTESC_FZ(LKS-1:LKE), start(1:1) )
          lzh_bnds(1,LKS-1) = LAND_GRID_CARTESC_FZ(LKS-1)
          do k = LKS-1, LKE-1
             lzh_bnds(2,k  ) = LAND_GRID_CARTESC_CZ(k+1)
             lzh_bnds(1,k+1) = LAND_GRID_CARTESC_CZ(k+1)
          end do
          lzh_bnds(2,LKE) = LAND_GRID_CARTESC_FZ(LKE)
          call FILE_Write_AssociatedCoordinate( fid, 'lzh_bnds', lzh_bnds(:,LKS-1:LKE), start(1:1) )
       end if

       ! urban
       if ( UKMAX > 0 ) then
          call FILE_Write_Axis( fid, 'uz', URBAN_GRID_CARTESC_CZ(UKS:UKE), start(1:1) )
          do k = UKS, UKE
             uz_bnds(1,k) = URBAN_GRID_CARTESC_FZ(k-1)
             uz_bnds(2,k) = URBAN_GRID_CARTESC_FZ(k  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'uz_bnds', uz_bnds(:,UKS:UKE), start(1:1) )

          call FILE_Write_Axis( fid, 'uzh', URBAN_GRID_CARTESC_FZ(UKS-1:UKE), start(1:1) )
          uzh_bnds(1,UKS-1) = URBAN_GRID_CARTESC_FZ(UKS-1)
          do k = UKS-1, UKE-1
             uzh_bnds(2,k  ) = URBAN_GRID_CARTESC_CZ(k+1)
             uzh_bnds(1,k+1) = URBAN_GRID_CARTESC_CZ(k+1)
          end do
          uzh_bnds(2,UKE) = URBAN_GRID_CARTESC_FZ(UKE)
          call FILE_Write_AssociatedCoordinate( fid, 'uzh_bnds', uzh_bnds(:,UKS-1:UKE), start(1:1) )
       end if
    end if

    if ( put_x ) then
       if ( FILE_get_aggregate(fid) ) then
          call FILE_Write_Axis( fid, 'x' ,  ATMOS_GRID_CARTESC_CX(ISB2:IEB2),  start(2:2) )
          do i = ISB2, IEB2
             x_bnds(1,i) = ATMOS_GRID_CARTESC_FX(i-1)
             x_bnds(2,i) = ATMOS_GRID_CARTESC_FX(i  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'x_bnds', x_bnds(:,ISB2:IEB2), (/1,start(2)/) )

          call FILE_Write_Axis( fid, 'xh',  ATMOS_GRID_CARTESC_FX(ISB2:IEB2),  start(2:2) )
          do i = ISB2, IEB2-1
             xh_bnds(1,i) = ATMOS_GRID_CARTESC_CX(i  )
             xh_bnds(2,i) = ATMOS_GRID_CARTESC_CX(i+1)
          end do
          xh_bnds(1,IEB2) = ATMOS_GRID_CARTESC_CX(IEB2)
          if ( IEB2 == IA ) then
             xh_bnds(2,IEB2) = ATMOS_GRID_CARTESC_FX(IEB2)
          else
             xh_bnds(2,IEB2) = ATMOS_GRID_CARTESC_CX(IEB2+1)
          end if
          call FILE_Write_AssociatedCoordinate( fid, 'xh_bnds', xh_bnds(:,ISB2:IEB2), (/1,start(2)/) )
       else
          call FILE_Write_Axis( fid, 'x' ,  ATMOS_GRID_CARTESC_CX(ISB:IEB),  start(2:2) )
          do i = ISB2, IEB
             x_bnds(1,i) = ATMOS_GRID_CARTESC_FX(i-1)
             x_bnds(2,i) = ATMOS_GRID_CARTESC_FX(i  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'x_bnds', x_bnds(:,ISB:IEB), (/1,start(2)/) )

          call FILE_Write_Axis( fid, 'xh',  ATMOS_GRID_CARTESC_FX(ISB:IEB),  start(2:2) )
          do i = ISB, IEB-1
             xh_bnds(1,i) = ATMOS_GRID_CARTESC_CX(i  )
             xh_bnds(2,i) = ATMOS_GRID_CARTESC_CX(i+1)
          end do
          xh_bnds(1,IEB) = ATMOS_GRID_CARTESC_CX(IEB)
          if ( IEB == IA ) then
             xh_bnds(2,IEB) = ATMOS_GRID_CARTESC_FX(IEB)
          else
             xh_bnds(2,IEB) = ATMOS_GRID_CARTESC_CX(IEB+1)
          end if
          call FILE_Write_AssociatedCoordinate( fid, 'xh_bnds', xh_bnds(:,ISB:IEB), (/1,start(2)/) )
       end if
    end if

    if ( put_y ) then
       if ( FILE_get_aggregate(fid) ) then
          call FILE_Write_Axis( fid, 'y' ,  ATMOS_GRID_CARTESC_CY(JSB2:JEB2),  start(3:3) )
          do j = JSB2, JEB2
             y_bnds(1,j) = ATMOS_GRID_CARTESC_FY(j-1)
             y_bnds(2,j) = ATMOS_GRID_CARTESC_FY(j  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'y_bnds', y_bnds(:,JSB2:JEB2), (/1,start(2)/) )

          call FILE_Write_Axis( fid, 'yh',  ATMOS_GRID_CARTESC_FY(JSB2:JEB2),  start(3:3) )
          do j = JSB2, JEB2-1
             yh_bnds(1,j) = ATMOS_GRID_CARTESC_CY(j  )
             yh_bnds(2,j) = ATMOS_GRID_CARTESC_CY(j+1)
          end do
          yh_bnds(1,JEB2) = ATMOS_GRID_CARTESC_CY(JEB2)
          if ( JEB2 == JA ) then
             yh_bnds(2,JEB2) = ATMOS_GRID_CARTESC_FY(JEB2)
          else
             yh_bnds(2,JEB2) = ATMOS_GRID_CARTESC_CY(JEB2+1)
          end if
          call FILE_Write_AssociatedCoordinate( fid, 'yh_bnds', yh_bnds(:,JSB2:JEB2), (/1,start(2)/) )
       else
          call FILE_Write_Axis( fid, 'y' ,  ATMOS_GRID_CARTESC_CY(JSB:JEB),  start(3:3) )
          do j = JSB, JEB
             y_bnds(1,j) = ATMOS_GRID_CARTESC_FY(j-1)
             y_bnds(2,j) = ATMOS_GRID_CARTESC_FY(j  )
          end do
          call FILE_Write_AssociatedCoordinate( fid, 'y_bnds', y_bnds(:,JSB:JEB), (/1,start(2)/) )

          call FILE_Write_Axis( fid, 'yh',  ATMOS_GRID_CARTESC_FY(JSB:JEB),  start(3:3) )
          do j = JSB, JEB-1
             yh_bnds(1,j) = ATMOS_GRID_CARTESC_CY(j  )
             yh_bnds(2,j) = ATMOS_GRID_CARTESC_CY(j+1)
          end do
          yh_bnds(1,JEB) = ATMOS_GRID_CARTESC_CY(JEB2)
          if ( JEB == JA ) then
             yh_bnds(2,JEB) = ATMOS_GRID_CARTESC_FY(JEB)
          else
             yh_bnds(2,JEB) = ATMOS_GRID_CARTESC_CY(JEB+1)
          end if
          call FILE_Write_AssociatedCoordinate( fid, 'yh_bnds', yh_bnds(:,JSB:JEB), (/1,start(2)/) )
       end if
    end if

    FDXG(1:IAG-1) = ATMOS_GRID_CARTESC_FDXG(:)
    FDXG(0  ) = UNDEF
    FDXG(IAG) = UNDEF
    FDYG(1:JAG-1) = ATMOS_GRID_CARTESC_FDYG(:)
    FDYG(0  ) = UNDEF
    FDYG(JAG) = UNDEF

    FDX(1:IA-1) = ATMOS_GRID_CARTESC_FDX(:)
    FDX(0 ) = FDXG(IS_inG-IHALO-1)
    FDX(IA) = FDXG(IE_inG+IHALO  )
    FDY(1:JA-1) = ATMOS_GRID_CARTESC_FDY(:)
    FDY(0 ) = FDYG(JS_inG-JHALO-1)
    FDY(JA) = FDYG(JE_inG+JHALO  )

    ! global coordinates (always including halo)
    if ( haszcoord .and. put_z ) then
       call FILE_Write_Axis( fid, 'CZ'  , ATMOS_GRID_CARTESC_CZ  (:), start(1:1) )
       call FILE_Write_Axis( fid, 'FZ'  , ATMOS_GRID_CARTESC_FZ  (:), start(1:1) )
       call FILE_Write_Axis( fid, 'CDZ' , ATMOS_GRID_CARTESC_CDZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'FDZ' , ATMOS_GRID_CARTESC_FDZ (:), start(1:1) )
       call FILE_Write_Axis( fid, 'CBFZ', ATMOS_GRID_CARTESC_CBFZ(:), start(1:1) )
       call FILE_Write_Axis( fid, 'FBFZ', ATMOS_GRID_CARTESC_FBFZ(:), start(1:1) )

       if ( OKMAX > 0 ) then
          call FILE_Write_Axis( fid, 'OCZ' , OCEAN_GRID_CARTESC_CZ (:), start(1:1) )
          call FILE_Write_Axis( fid, 'OFZ' , OCEAN_GRID_CARTESC_FZ (:), start(1:1) )
          call FILE_Write_Axis( fid, 'OCDZ', OCEAN_GRID_CARTESC_CDZ(:), start(1:1) )
       end if

       if ( LKMAX > 0 ) then
          call FILE_Write_Axis( fid, 'LCZ' , LAND_GRID_CARTESC_CZ (:), start(1:1) )
          call FILE_Write_Axis( fid, 'LFZ' , LAND_GRID_CARTESC_FZ (:), start(1:1) )
          call FILE_Write_Axis( fid, 'LCDZ', LAND_GRID_CARTESC_CDZ(:), start(1:1) )
       end if

       if ( UKMAX > 0 ) then
          call FILE_Write_Axis( fid, 'UCZ' , URBAN_GRID_CARTESC_CZ (:), start(1:1) )
          call FILE_Write_Axis( fid, 'UFZ' , URBAN_GRID_CARTESC_FZ (:), start(1:1) )
          call FILE_Write_Axis( fid, 'UCDZ', URBAN_GRID_CARTESC_CDZ(:), start(1:1) )
       end if
    end if

    if ( FILE_get_AGGREGATE(fid) ) then
       if ( PRC_IsMaster ) then
          call FILE_Write_Axis( fid, 'CX',   ATMOS_GRID_CARTESC_CXG  (:) )
          call FILE_Write_Axis( fid, 'CY',   ATMOS_GRID_CARTESC_CYG  (:) )
          call FILE_Write_Axis( fid, 'FX',   ATMOS_GRID_CARTESC_FXG  (:) )
          call FILE_Write_Axis( fid, 'FY',   ATMOS_GRID_CARTESC_FYG  (:) )
          call FILE_Write_Axis( fid, 'CDX',  ATMOS_GRID_CARTESC_CDXG (:) )
          call FILE_Write_Axis( fid, 'CDY',  ATMOS_GRID_CARTESC_CDYG (:) )

          call FILE_Write_Axis( fid, 'FDX',                     FDXG (:) )
          call FILE_Write_Axis( fid, 'FDY',                     FDYG (:) )
          call FILE_Write_Axis( fid, 'CBFX', ATMOS_GRID_CARTESC_CBFXG(:) )
          call FILE_Write_Axis( fid, 'CBFY', ATMOS_GRID_CARTESC_CBFYG(:) )
          call FILE_Write_Axis( fid, 'FBFX', ATMOS_GRID_CARTESC_FBFXG(:) )
          call FILE_Write_Axis( fid, 'FBFY', ATMOS_GRID_CARTESC_FBFYG(:) )
       endif
    else
       call FILE_Write_Axis( fid, 'CX',   ATMOS_GRID_CARTESC_CX  (:) )
       call FILE_Write_Axis( fid, 'CY',   ATMOS_GRID_CARTESC_CY  (:) )
       call FILE_Write_Axis( fid, 'FX',   ATMOS_GRID_CARTESC_FX  (:) )
       call FILE_Write_Axis( fid, 'FY',   ATMOS_GRID_CARTESC_FY  (:) )
       call FILE_Write_Axis( fid, 'CDX',  ATMOS_GRID_CARTESC_CDX (:) )
       call FILE_Write_Axis( fid, 'CDY',  ATMOS_GRID_CARTESC_CDY (:) )
       call FILE_Write_Axis( fid, 'FDX',                     FDX (:) )
       call FILE_Write_Axis( fid, 'FDY',                     FDY (:) )
       call FILE_Write_Axis( fid, 'CBFX', ATMOS_GRID_CARTESC_CBFX(:) )
       call FILE_Write_Axis( fid, 'CBFY', ATMOS_GRID_CARTESC_CBFY(:) )
       call FILE_Write_Axis( fid, 'FBFX', ATMOS_GRID_CARTESC_FBFX(:) )
       call FILE_Write_Axis( fid, 'FBFY', ATMOS_GRID_CARTESC_FBFY(:) )
    endif

    call FILE_Write_Axis( fid, 'CXG',   ATMOS_GRID_CARTESC_CXG  (:) )
    call FILE_Write_Axis( fid, 'CYG',   ATMOS_GRID_CARTESC_CYG  (:) )
    call FILE_Write_Axis( fid, 'FXG',   ATMOS_GRID_CARTESC_FXG  (:) )
    call FILE_Write_Axis( fid, 'FYG',   ATMOS_GRID_CARTESC_FYG  (:) )
    call FILE_Write_Axis( fid, 'CDXG',  ATMOS_GRID_CARTESC_CDXG (:) )
    call FILE_Write_Axis( fid, 'CDYG',  ATMOS_GRID_CARTESC_CDYG (:) )
    call FILE_Write_Axis( fid, 'FDXG',                     FDXG (:) )
    call FILE_Write_Axis( fid, 'FDYG',                     FDYG (:) )
    call FILE_Write_Axis( fid, 'CBFXG', ATMOS_GRID_CARTESC_CBFXG(:) )
    call FILE_Write_Axis( fid, 'CBFYG', ATMOS_GRID_CARTESC_CBFYG(:) )
    call FILE_Write_Axis( fid, 'FBFXG', ATMOS_GRID_CARTESC_FBFXG(:) )
    call FILE_Write_Axis( fid, 'FBFYG', ATMOS_GRID_CARTESC_FBFYG(:) )

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
          call FILE_Write_AssociatedCoordinate( fid, 'topo',   AXIS_TOPO  (:,:), start(2:3) )
       end if
       call FILE_Write_AssociatedCoordinate( fid, 'lsmask', AXIS_LSMASK(:,:), start(2:3) )

       call FILE_Write_AssociatedCoordinate( fid, 'cell_area',    AXIS_AREA  (:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'cell_area_uy', AXIS_AREAUY(:,:), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'cell_area_xv', AXIS_AREAXV(:,:), start(2:3) )

       if ( haszcoord ) then
          call FILE_Write_AssociatedCoordinate( fid, 'height'    , AXIS_HGT   (:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'height_wxy', AXIS_HGTWXY(:,:,:), start(1:3) )

          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zuy_x', AXIS_AREAZUY_X(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zxv_y', AXIS_AREAZXV_Y(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_wuy_x', AXIS_AREAWUY_X(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_wxv_y', AXIS_AREAWXV_Y(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zxy_x', AXIS_AREAZXY_X(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zuv_y', AXIS_AREAZUV_Y(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zuv_x', AXIS_AREAZUV_X(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zxy_y', AXIS_AREAZXY_Y(:,:,:), start(1:3) )

          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume',     AXIS_VOL   (:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_wxy', AXIS_VOLWXY(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_zuy', AXIS_VOLZUY(:,:,:), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_zxv', AXIS_VOLZXV(:,:,:), start(1:3) )

          if ( OKMAX > 0 ) then
             call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_oxy', AXIS_VOLO(:,:,:), start(1:3) )
          end if
          if ( LKMAX > 0 ) then
             call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_lxy', AXIS_VOLL(:,:,:), start(1:3) )
          end if
          if ( UKMAX > 0 ) then
             call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_uxy', AXIS_VOLU(:,:,:), start(1:3) )
          end if
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
          call FILE_Write_AssociatedCoordinate( fid, 'topo',   AXIS_TOPO  (XSB:XEB,YSB:YEB), start(2:3) )
       end if
       call FILE_Write_AssociatedCoordinate( fid, 'lsmask', AXIS_LSMASK(XSB:XEB,YSB:YEB), start(2:3) )

       call FILE_Write_AssociatedCoordinate( fid, 'cell_area',    AXIS_AREA  (XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'cell_area_uy', AXIS_AREAUY(XSB:XEB,YSB:YEB), start(2:3) )
       call FILE_Write_AssociatedCoordinate( fid, 'cell_area_xv', AXIS_AREAXV(XSB:XEB,YSB:YEB), start(2:3) )

       if ( haszcoord ) then
          call FILE_Write_AssociatedCoordinate( fid, 'height'    , AXIS_HGT   (:,XSB:XEB,YSB:YEB), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'height_wxy', AXIS_HGTWXY(:,XSB:XEB,YSB:YEB), start(1:3) )

          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zuy_x', AXIS_AREAZUY_X(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zxv_y', AXIS_AREAZXV_Y(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_wuy_x', AXIS_AREAWUY_X(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_wxv_y', AXIS_AREAWXV_Y(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zxy_x', AXIS_AREAZXY_X(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zuv_y', AXIS_AREAZUV_Y(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zuv_x', AXIS_AREAZUV_X(:,XSB:XEB,YSB:YEB), start(2:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_area_zxy_y', AXIS_AREAZXY_Y(:,XSB:XEB,YSB:YEB), start(2:3) )

          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume',     AXIS_VOL   (:,XSB:XEB,YSB:YEB), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_wxy', AXIS_VOLWXY(:,XSB:XEB,YSB:YEB), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_zuy', AXIS_VOLZUY(:,XSB:XEB,YSB:YEB), start(1:3) )
          call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_zxv', AXIS_VOLZXV(:,XSB:XEB,YSB:YEB), start(1:3) )

          if ( OKMAX > 0 ) then
             call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_oxy', AXIS_VOLO(:,XSB:XEB,YSB:YEB), start(1:3) )
          end if
          if ( LKMAX > 0 ) then
             call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_lxy', AXIS_VOLL(:,XSB:XEB,YSB:YEB), start(1:3) )
          end if
          if ( UKMAX > 0 ) then
             call FILE_Write_AssociatedCoordinate( fid, 'cell_volume_uxy', AXIS_VOLU(:,XSB:XEB,YSB:YEB), start(1:3) )
          end if
       end if
    end if

    return
  end subroutine FILE_CARTESC_write_axes

  !-----------------------------------------------------------------------------
  !> Define a variable to file
  subroutine FILE_CARTESC_def_var( &
       fid,                 &
       varname, desc, unit, &
       dim_type, datatype,  &
       vid,                 &
       standard_name,       &
       timeintv, nsteps,    &
       cell_measures        )
    use scale_file_h, only: &
       FILE_REAL8, &
       FILE_REAL4
    use scale_file, only: &
       FILE_opened, &
       FILE_Def_Variable, &
       FILE_Set_Attribute
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_mapprojection, only: &
       MAPPROJECTION_mappinginfo
    implicit none

    integer,          intent(in)  :: fid      !< file ID
    character(len=*), intent(in)  :: varname  !< name        of the variable
    character(len=*), intent(in)  :: desc     !< description of the variable
    character(len=*), intent(in)  :: unit     !< unit        of the variable
    character(len=*), intent(in)  :: dim_type !< axis type (Z/X/Y)
    character(len=*), intent(in)  :: datatype !< data type (REAL8/REAL4/default)

    integer, intent(out) :: vid      !< variable ID

    character(len=*), intent(in), optional :: standard_name !< standard name of the variable
    real(DP),         intent(in), optional :: timeintv      !< time interval [sec]
    integer,          intent(in), optional :: nsteps        !< number of time steps
    character(len=*), intent(in), optional :: cell_measures !< "area" or "volume"

    character(len=H_MID)   :: standard_name_
    character(len=H_SHORT) :: cell_measures_

    character(len=H_SHORT) :: dimtype

    integer :: dtype, elm_size, ndims
    integer :: dimid, n
    !---------------------------------------------------------------------------

    if ( .not. FILE_opened(fid) ) return

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
          LOG_ERROR("FILE_CARTESC_def_var",*) 'unsupported data type. Check!', trim(datatype), ' item:',trim(varname)
          call PRC_abort
       endif
    endif

    if ( PRC_TwoD ) then
       select case( dim_type )
       case ( "UY" )
          dimtype = "XY"
       case ( "UV" )
          dimtype = "XV"
       case ( "ZXHY" )
          dimtype = "ZXY"
       case ( "ZXHYH")
          dimtype = "ZXYH"
       case ( "ZHXHY" )
          dimtype = "ZHXY"
       case default
          dimtype = dim_type
       end select
    else
       dimtype = dim_type
    end if

    dimid = -1
    do n = 1, FILE_CARTESC_ndims
       if ( FILE_CARTESC_dims(n)%name == dimtype ) then
          dimid = n
          exit
       end if
    end do
    if ( dimid <= -1 ) then
       LOG_ERROR("FILE_CARTESC_def_var",*) 'dim_type is not supported: ', trim(dimtype)
       call PRC_abort
    end if

    if ( present(nsteps) ) then
       write_buf_amount(fid) = write_buf_amount(fid) + FILE_CARTESC_dims(dimid)%size * elm_size * nsteps
    else
       write_buf_amount(fid) = write_buf_amount(fid) + FILE_CARTESC_dims(dimid)%size * elm_size
    end if

    ndims = FILE_CARTESC_dims(dimid)%ndims

    if ( present(standard_name) ) then
       standard_name_ = standard_name
    else
       standard_name_ = ""
    end if

    if ( present(timeintv) ) then  ! 3D/4D variable with time dimension
      call FILE_Def_Variable( fid, varname, desc, unit, standard_name_,             & ! [IN]
                              ndims, FILE_CARTESC_dims(dimid)%dims(1:ndims), dtype, & ! [IN]
                              vid,                                                  & ! [OUT]
                              time_int=timeintv                                     ) ! [IN]
    else
      call FILE_Def_Variable( fid, varname, desc, unit, standard_name_,             & ! [IN]
                              ndims, FILE_CARTESC_dims(dimid)%dims(1:ndims), dtype, & ! [IN]
                              vid                                                   ) ! [OUT]
    endif

    if ( present(cell_measures) ) then
       cell_measures_ = cell_measures
       select case ( cell_measures )
       case ( "area" )
          if ( FILE_CARTESC_dims(dimid)%area == "" ) then
             LOG_ERROR("FILE_CARTESC_def_var",*) 'area is not supported for ', trim(dimtype), ' as cell_measures'
             call PRC_abort
          end if
       case ( "area_z" )
          if ( FILE_CARTESC_dims(dimid)%area == "" ) then
             LOG_ERROR("FILE_CARTESC_def_var",*) 'area_z is not supported for ', trim(dimtype), ' as cell_measures'
             call PRC_abort
          end if
       case ( "area_x" )
          if ( FILE_CARTESC_dims(dimid)%area_x == "" ) then
             LOG_ERROR("FILE_CARTESC_def_var",*) 'area_x is not supported for ', trim(dimtype), ' as cell_measures'
             call PRC_abort
          end if
       case ( "area_y" )
          if ( FILE_CARTESC_dims(dimid)%area_y == "" ) then
             LOG_ERROR("FILE_CARTESC_def_var",*) 'area_y is not supported for ', trim(dimtype), ' as cell_measures'
             call PRC_abort
          end if
       case ( "volume" )
          if ( FILE_CARTESC_dims(dimid)%volume == "" ) then
             LOG_ERROR("FILE_CARTESC_def_var",*) 'volume is not supported for ', trim(dimtype), ' as cell_measures'
             call PRC_abort
          end if
       case default
          LOG_ERROR("FILE_CARTESC_def_var",*) 'cell_measures must be "area" or "volume"'
          call PRC_abort
       end select
    else if ( ndims == 2 ) then
       cell_measures_ = "area"
    else if ( ndims == 3 ) then
       cell_measures_ = "volume"
    else
       cell_measures_ = ""
    end if

    select case( cell_measures_ )
    case ( "area", "area_z" )
       call FILE_Set_Attribute( fid, varname, "cell_measures", "area: "//trim(FILE_CARTESC_dims(dimid)%area) )
    case ( "area_x" )
       call FILE_Set_Attribute( fid, varname, "cell_measures", "area: "//trim(FILE_CARTESC_dims(dimid)%area_x) )
    case ( "area_y" )
       call FILE_Set_Attribute( fid, varname, "cell_measures", "area: "//trim(FILE_CARTESC_dims(dimid)%area_y) )
    case ( "volume" )
       call FILE_Set_Attribute( fid, varname, "cell_measures", "volume: "//trim(FILE_CARTESC_dims(dimid)%volume) )
    end select

    ! mapping
    if ( FILE_CARTESC_dims(dimid)%mapping .and. MAPPROJECTION_mappinginfo%mapping_name /= "" ) then
       call FILE_Set_Attribute( fid, varname, "grid_mapping", MAPPROJECTION_mappinginfo%mapping_name )
    end if

    ! SGRID
    if ( FILE_CARTESC_dims(dimid)%location /= "" ) then
       call FILE_Set_Attribute( fid, varname, "grid",     FILE_CARTESC_dims(dimid)%grid     )
       call FILE_Set_Attribute( fid, varname, "location", FILE_CARTESC_dims(dimid)%location )
    end if

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
       FILE_opened, &
       FILE_Write
    use scale_prc, only: &
       PRC_myrank,  &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_2Drank
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
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

    if ( .not. FILE_opened(fid) ) return

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
       LOG_ERROR("FILE_CARTESC_write_var_1D",*) 'unsupported dimenstion type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if( exec ) call FILE_Write( vid, var(dim1_S:dim1_E),          & ! [IN]
                                NOWDAYSEC, NOWDAYSEC, start=start ) ! [IN]

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
       FILE_opened, &
       FILE_Write
    use scale_prc, only: &
       PRC_myrank,     &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_2Drank, &
       PRC_NUM_X, &
       PRC_NUM_Y
    use scale_time, only: &
       NOWDAYSEC  => TIME_NOWDAYSEC
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

    if ( .not. FILE_opened(fid) ) return

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
       LOG_ERROR("FILE_CARTESC_write_var_2D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if ( exec ) then
       if ( fill_halo_ ) then ! fill halo cells with RMISS
          do j = JS, JE
          do i = IS, IE
             varhalo(i,j) = var(i,j)
          end do
          end do

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
                           NOWDAYSEC, NOWDAYSEC, start                ) ! [IN]
       else
          call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E), &
                           NOWDAYSEC, NOWDAYSEC, start            ) ! [IN]
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
       FILE_opened, &
       FILE_Write
    use scale_prc, only: &
       PRC_myrank,  &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_2Drank, &
       PRC_NUM_X,  &
       PRC_NUM_Y
    use scale_time, only: &
       NOWDAYSEC  => TIME_NOWDAYSEC
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

    if ( .not. FILE_opened(fid) ) return

    call PROF_rapstart('FILE_O_NetCDF', 2)

    fill_halo_ = .false.
    if( present(fill_halo) ) fill_halo_ = fill_halo

    rankidx(1) = PRC_2Drank(PRC_myrank,1)
    rankidx(2) = PRC_2Drank(PRC_myrank,2)

    start(1) = 1
    start(2) = ISGA
    start(3) = JSGA

    if (      dim_type == 'ZXY'   &
         .OR. dim_type == 'ZXHY'  &
         .OR. dim_type == 'ZXYH'  &
         .OR. dim_type == 'ZXHYH' ) then
       dim1_max = KMAX
       dim1_S   = KS
       dim1_E   = KE
    elseif (  dim_type == 'ZHXY'  &
         .OR. dim_type == 'ZHXHY' &
         .OR. dim_type == 'ZHXYH' ) then
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
       LOG_ERROR("FILE_CARTESC_write_var_3D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
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
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = 1, dim1_max
          varhalo(k,i,j) = var(k,i,j)
       enddo
       enddo
       enddo

       ! W halo
       do j = 1, JA
       do i = 1, IS-1
       do k = 1, dim1_max
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo
       ! E halo
       do j = 1, JA
       do i = IE+1, IA
       do k = 1, dim1_max
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo
       ! S halo
       do j = 1, JS-1
       do i = 1, IA
       do k = 1, dim1_max
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo
       ! N halo
       do j = JE+1, JA
       do i = 1, IA
       do k = 1, dim1_max
          varhalo(k,i,j) = RMISS
       enddo
       enddo
       enddo

       call FILE_Write( vid, varhalo(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                        NOWDAYSEC, NOWDAYSEC, start                              ) ! [IN]
    else
       call FILE_Write( vid, var(dim1_S:dim1_E,dim2_S:dim2_E,dim3_S:dim3_E), &
                        NOWDAYSEC, NOWDAYSEC, start                          ) ! [IN]
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
       FILE_opened, &
       FILE_Write
    use scale_prc, only: &
       PRC_myrank,     &
       PRC_abort
    use scale_prc_cartesC, only: &
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

    if ( .not. FILE_opened(fid) ) return

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
       LOG_ERROR("FILE_CARTESC_write_var_3D_t",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    start(1) = ISGA
    start(2) = JSGA
    ! start(3) time dimension will be set in file_write_data()

    if ( present(timetarg) ) then
       nowtime = timeofs_ + (timetarg-1) * time_interval

       if ( fill_halo_ ) then
          do j = JS, JE
          do i = IS, IE
             varhalo(i,j) = var(i,j,timetarg)
          end do
          end do

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
             do j = JS, JE
             do i = IS, IE
                varhalo(i,j) = var(i,j,n)
             end do
             end do

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
       FILE_opened, &
       FILE_Write, &
       FILE_Flush
    use scale_prc, only: &
       PRC_myrank,     &
       PRC_abort
    use scale_prc_cartesC, only: &
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

    if ( .not. FILE_opened(fid) ) return

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
       LOG_ERROR("FILE_CARTESC_write_var_4D",*) 'unsupported dimension type. Check! dim_type:', trim(dim_type), ', item:',trim(varname)
       call PRC_abort
    endif

    if ( present(timetarg) ) then
       nowtime = timeofs_ + (timetarg-1) * time_interval

       if ( fill_halo_ ) then
          do j = JS, JE
          do i = IS, IE
          do k = 1, dim1_max
             varhalo(k,i,j) = var(k,i,j,timetarg)
          end do
          end do
          end do

          ! W halo
          do j = 1, JA
          do i = 1, IS-1
          do k = 1, dim1_max
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo
          ! E halo
          do j = 1, JA
          do i = IE+1, IA
          do k = 1, dim1_max
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo
          ! S halo
          do j = 1, JS-1
          do i = 1, IA
          do k = 1, dim1_max
             varhalo(k,i,j) = RMISS
          enddo
          enddo
          enddo
          ! N halo
          do j = JE+1, JA
          do i = 1, IA
          do k = 1, dim1_max
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
             do j = JS, JE
             do i = IS, IE
             do k = 1, dim1_max
                varhalo(k,i,j) = var(k,i,j,n)
             end do
             end do
             end do

             ! W halo
             do j = 1, JA
             do i = 1, IS-1
             do k = 1, dim1_max
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo
             ! E halo
             do j = 1, JA
             do i = IE+1, IA
             do k = 1, dim1_max
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo
             ! S halo
             do j = 1, JS-1
             do i = 1, IA
             do k = 1, dim1_max
                varhalo(k,i,j) = RMISS
             enddo
             enddo
             enddo
             ! N halo
             do j = JE+1, JA
             do i = 1, IA
             do k = 1, dim1_max
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

  ! private procedures

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
    use scale_prc, only: &
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
       LOG_ERROR("check_1d",*) 'size of coordinate ('//trim(name)//') is different:', nmax, size(buffer)
       call PRC_abort
    endif

    do n=1, nmax
       if ( abs(expected(n)) > EPS ) then
          check = abs(buffer(n)-expected(n)) / abs(buffer(n)+expected(n)) * 2.0_RP
       else
          check = abs(buffer(n)-expected(n))
       endif

       if ( check > FILE_CARTESC_datacheck_criteria ) then
          LOG_ERROR("check_1d",*) 'value of coordinate ('//trim(name)//') at ', n, ' is different:', &
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
    use scale_prc, only: &
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
       LOG_ERROR("check_2d",*) 'the first size of coordinate ('//trim(name)//') is different:', imax, size(buffer,1)
       call PRC_abort
    endif
    if ( size(buffer,2) /= jmax ) then
       LOG_ERROR("check_2d",*) 'the second size of coordinate ('//trim(name)//') is different:', jmax, size(buffer,2)
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
          LOG_ERROR("check_2d",*) 'value of coordinate ('//trim(name)//') at (', i, ',', j, ') is different:', &
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
    use scale_prc, only: &
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
       LOG_ERROR("check_3d",*) 'the first size of coordinate ('//trim(name)//') is different:', kmax, size(buffer,1)
       call PRC_abort
    endif
    if ( size(buffer,2) /= imax ) then
       LOG_ERROR("check_3d",*) 'the second size of coordinate ('//trim(name)//') is different:', imax, size(buffer,2)
       call PRC_abort
    endif
    if ( size(buffer,3) /= jmax ) then
       LOG_ERROR("check_3d",*) 'the third size of coordinate ('//trim(name)//') is different:', jmax, size(buffer,3)
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
             LOG_ERROR("check_3d",*) 'value of coordinate ('//trim(name)//') at ', i, ',', j, ',', k, ' is different:', &
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
             LOG_ERROR("check_3d",*) 'value of coordinate ('//trim(name)//') at ', k, ',', i, ',', j, ' is different:', &
                        expected(k,i,j), buffer(k,i,j), check
             call PRC_abort
          endif
       enddo
       enddo
       enddo
    endif

    return
  end subroutine check_3d

  subroutine set_dimension_informations
    use scale_prc_cartesC, only: &
       PRC_PERIODIC_X, &
       PRC_PERIODIC_Y, &
       PRC_NUM_X,      &
       PRC_NUM_Y,      &
       PRC_HAS_W,      &
       PRC_HAS_E,      &
       PRC_HAS_S,      &
       PRC_HAS_N
    implicit none
    !---------------------------------------------------------------------------

    ! Dimension Information

    ! 1D variable
    call set_dimension( 'Z', 1, 'z', KA )
    call set_dimension( 'X', 1, 'x', IA )
    call set_dimension( 'Y', 1, 'y', JA )

    ! 2D variable
    call set_dimension( 'XY', 2, (/ 'x' , 'y'  /), IA*JA, .true., area='cell_area',    location='face'  )
    call set_dimension( 'UY', 2, (/ 'xh', 'y ' /), IA*JA, .true., area='cell_area_uy', location='edge1' )
    call set_dimension( 'XV', 2, (/ 'x ', 'yh' /), IA*JA, .true., area='cell_area_xv', location='edge2' )
    call set_dimension( 'UV', 2, (/ 'xh', 'yh' /), IA*JA, .true.,                      location='node'  )

    call set_dimension( 'ZX', 2, (/ 'z', 'x' /), KA*IA )
    call set_dimension( 'ZY', 2, (/ 'z', 'y' /), KA*JA )

    ! 3D variable
    call set_dimension( 'ZXY',   3, (/ 'z',  'x',  'y'  /), KA*IA*JA,     .true.,                &
                        area='cell_area',    area_x='cell_area_zxy_x', area_y='cell_area_zxy_y', &
                        volume='cell_volume',     location='face'                                )
    call set_dimension( 'ZHXY',  3, (/ 'zh', 'x ', 'y ' /), (KA+1)*IA*JA, .true.,                &
                        area='cell_area',    area_x='cell_area_wxy_x', area_y='cell_area_wxy_y', &
                        volume='cell_volume_wxy', location='face'                                )
    call set_dimension( 'ZXHY',  3, (/ 'z ', 'xh', 'y ' /), KA*IA*JA,     .true.,                &
                        area='cell_area_uy', area_x='cell_area_zuy_x',                           &
                        volume='cell_volume_zuy', location='edge1'                               )
    call set_dimension( 'ZXYH',  3, (/ 'z ', 'x ', 'yh' /), KA*IA*JA,     .true.,                &
                        area='cell_area_xv',                           area_y='cell_area_zxv_y', &
                        volume='cell_volume_zxv', location='edge2'                               )
    call set_dimension( 'ZXHYH', 3, (/ 'z ', 'xh', 'yh' /), KA*IA*JA,     .true., &
                                             area_x='cell_area_zuv_x', area_y='cell_area_zuv_y', &
                                                  location='node'                                )
    call set_dimension( 'ZHXHY', 3, (/ 'zh', 'xh', 'y ' /), (KA+1)*IA*JA, .true., &
                                             area_x='cell_area_wuy_x',                           &
                                                  location='edge1'                               )
    call set_dimension( 'ZHXYH', 3, (/ 'zh', 'x ', 'yh' /), (KA+1)*IA*JA, .true., &
                                                                       area_y='cell_area_wxv_y', &
                                                  location='edge2'                               )

    if ( OKMAX > 0 ) then
       call set_dimension( 'OXY',  3, (/ 'oz',  'x ',  'y '  /), OKMAX*IA*JA,     .true., area='cell_area', volume='cell_volume_oxy', location='face', grid='ocean' )
       call set_dimension( 'OHXY', 3, (/ 'ozh', 'x  ', 'y  ' /), (OKMAX+1)*IA*JA, .true., area='cell_area', volume='cell_volume_oxy', location='face', grid='ocean' )
    end if
    if ( LKMAX > 0 ) then
       call set_dimension( 'LXY',  3, (/ 'lz',  'x ',  'y '  /), LKMAX*IA*JA,     .true., area='cell_area', volume='cell_volume_lxy', location='face', grid='land'  )
       call set_dimension( 'LHXY', 3, (/ 'lzh', 'x  ', 'y  ' /), (LKMAX+1)*IA*JA, .true., area='cell_area', volume='cell_volume_lxy', location='face', grid='land'  )
    end if
    if ( UKMAX > 0 ) then
       call set_dimension( 'UXY',  3, (/ 'uz',  'x ',  'y '  /), UKMAX*IA*JA,     .true., area='cell_area', volume='cell_volume_uxy', location='face', grid='urban' )
       call set_dimension( 'UHXY', 3, (/ 'uzh', 'x  ', 'y  ' /), (UKMAX+1)*IA*JA, .true., area='cell_area', volume='cell_volume_uxy', location='face', grid='urban' )
    end if


    ! Axis information

    if ( PRC_PERIODIC_X ) then
       FILE_CARTESC_AXIS_info(1)%periodic = .true.
       FILE_CARTESC_AXIS_info(2)%periodic = .true.
    else
       FILE_CARTESC_AXIS_info(1)%periodic = .false.
       FILE_CARTESC_AXIS_info(2)%periodic = .false.
    endif

    if ( PRC_PERIODIC_Y ) then
       FILE_CARTESC_AXIS_info(3)%periodic = .true.
       FILE_CARTESC_AXIS_info(4)%periodic = .true.
    else
       FILE_CARTESC_AXIS_info(3)%periodic = .false.
       FILE_CARTESC_AXIS_info(4)%periodic = .false.
    endif


    ! for x
    if ( PRC_PERIODIC_X ) then
       FILE_CARTESC_AXIS_info(1)%size_global (1) = IMAX * PRC_NUM_X
       FILE_CARTESC_AXIS_info(1)%start_global(1) = IS_inG - IHALO
       FILE_CARTESC_AXIS_info(1)%halo_global (1) = 0     ! west side
       FILE_CARTESC_AXIS_info(1)%halo_global (2) = 0     ! east side
       FILE_CARTESC_AXIS_info(1)%halo_local  (1) = 0     ! west side
       FILE_CARTESC_AXIS_info(1)%halo_local  (2) = 0     ! east side
    else
       FILE_CARTESC_AXIS_info(1)%size_global (1) = IAG
       FILE_CARTESC_AXIS_info(1)%start_global(1) = ISGA
       FILE_CARTESC_AXIS_info(1)%halo_global (1) = IHALO ! west side
       FILE_CARTESC_AXIS_info(1)%halo_global (2) = IHALO ! east side
       FILE_CARTESC_AXIS_info(1)%halo_local  (1) = IHALO ! west side
       FILE_CARTESC_AXIS_info(1)%halo_local  (2) = IHALO ! east side
       if( PRC_HAS_W ) FILE_CARTESC_AXIS_info(1)%halo_local(1) = 0
       if( PRC_HAS_E ) FILE_CARTESC_AXIS_info(1)%halo_local(2) = 0
    endif
    ! for xh
    FILE_CARTESC_AXIS_info(2) = FILE_CARTESC_AXIS_info(1)

    ! for y
    if ( PRC_PERIODIC_Y ) then
       FILE_CARTESC_AXIS_info(3)%size_global (1) = JMAX * PRC_NUM_Y
       FILE_CARTESC_AXIS_info(3)%start_global(1) = JS_inG - JHALO
       FILE_CARTESC_AXIS_info(3)%halo_global (1) = 0     ! south side
       FILE_CARTESC_AXIS_info(3)%halo_global (2) = 0     ! north side
       FILE_CARTESC_AXIS_info(3)%halo_local  (1) = 0     ! south side
       FILE_CARTESC_AXIS_info(3)%halo_local  (2) = 0     ! north side
    else
       FILE_CARTESC_AXIS_info(3)%size_global (1) = JAG
       FILE_CARTESC_AXIS_info(3)%start_global(1) = JSGA
       FILE_CARTESC_AXIS_info(3)%halo_global (1) = JHALO ! south side
       FILE_CARTESC_AXIS_info(3)%halo_global (2) = JHALO ! north side
       FILE_CARTESC_AXIS_info(3)%halo_local  (1) = JHALO ! south side
       FILE_CARTESC_AXIS_info(3)%halo_local  (2) = JHALO ! north side
       if( PRC_HAS_S ) FILE_CARTESC_AXIS_info(3)%halo_local(1) = 0
       if( PRC_HAS_N ) FILE_CARTESC_AXIS_info(3)%halo_local(2) = 0
    endif
    ! for yh
    FILE_CARTESC_AXIS_info(4) = FILE_CARTESC_AXIS_info(3)


    return
  end subroutine set_dimension_informations

  subroutine set_dimension( name, ndims, dims, size, mapping, area, area_x, area_y, volume, location, grid )
    use scale_prc, only: &
       PRC_abort
    character(len=*), intent(in) :: name
    integer,          intent(in) :: ndims
    character(len=*), intent(in) :: dims(ndims)
    integer,          intent(in) :: size
    logical,          intent(in), optional :: mapping
    character(len=*), intent(in), optional :: area
    character(len=*), intent(in), optional :: area_x
    character(len=*), intent(in), optional :: area_y
    character(len=*), intent(in), optional :: volume
    character(len=*), intent(in), optional :: location
    character(len=*), intent(in), optional :: grid

    integer, save :: dimid = 0

    integer :: n

    do n = 1, 2
       dimid = dimid + 1
       if ( dimid > FILE_CARTESC_ndims ) then
          LOG_ERROR("set_dimension",*) 'number of dimensions exceeds the limit', dimid, FILE_CARTESC_ndims
          call PRC_abort
       end if

       if ( n==1 ) then
          FILE_CARTESC_dims(dimid)%name = name
       else
          FILE_CARTESC_dims(dimid)%name = trim(name)//"T"
       end if
       FILE_CARTESC_dims(dimid)%ndims         = ndims
       FILE_CARTESC_dims(dimid)%dims(1:ndims) = dims(:)
       FILE_CARTESC_dims(dimid)%size          = size

       if ( present(mapping) ) then
          FILE_CARTESC_dims(dimid)%mapping  = mapping
       else
          FILE_CARTESC_dims(dimid)%mapping  = .false.
       end if

       if ( present(area) ) then
          FILE_CARTESC_dims(dimid)%area     = area
       else
          FILE_CARTESC_dims(dimid)%area     = ''
       end if
       if ( present(area_x) ) then
          FILE_CARTESC_dims(dimid)%area     = area_x
       else
          FILE_CARTESC_dims(dimid)%area     = ''
       end if
       if ( present(area_y) ) then
          FILE_CARTESC_dims(dimid)%area     = area_y
       else
          FILE_CARTESC_dims(dimid)%area     = ''
       end if

       if ( present(volume) ) then
          FILE_CARTESC_dims(dimid)%volume   = volume
       else
          FILE_CARTESC_dims(dimid)%volume   = ''
       end if

       if ( present(location) ) then
          FILE_CARTESC_dims(dimid)%location = location
          if ( present(grid) ) then
             FILE_CARTESC_dims(dimid)%grid = 'grid_'//trim(grid)
          else
             FILE_CARTESC_dims(dimid)%grid = 'grid'
          end if
       else
          FILE_CARTESC_dims(dimid)%location = ''
          FILE_CARTESC_dims(dimid)%grid     = ''
       end if

    end do

    return
  end subroutine set_dimension

  !-----------------------------------------------------------------------------
  !> construct MPI derived datatypes for read buffers
  subroutine Construct_Derived_Datatype
    use mpi
    use scale_prc_cartesC, only: &
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

    if( RP == 8 ) etype = MPI_DOUBLE

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

    if ( OKMAX > 0 ) then
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
    end if

    if ( LKMAX > 0 ) then
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
    end if

    if ( UKMAX > 0 ) then
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
    end if

    ! for dim_type == 'ZX'
    startZX(1)  = KHALO+1
    startZX(2)  = IS_inG - IHALO
    countZX(1)  = KMAX
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

end module scale_file_cartesC
