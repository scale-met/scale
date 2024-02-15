!-------------------------------------------------------------------------------
!> module REAL input netCDF
!!
!! @par Description
!!          read data from netCDF file for real atmospheric simulations
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module mod_realinput_netcdf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_file_h
  use scale_tracer
  use scale_cpl_sfc_index
  use scale_hash

  use scale_prc, only: &
     myrank => PRC_myrank,  &
     PRC_abort
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ParentAtmosSetupNetCDF
  public :: ParentAtmosOpenNetCDF
  public :: ParentAtmosFinalizeNetCDF
  public :: ParentAtmosInputNetCDF
  public :: ParentLandSetupNetCDF
  public :: ParentLandOpenNetCDF
  public :: ParentLandFinalizeNetCDF
  public :: ParentLandInputNetCDF
  public :: ParentOceanSetupNetCDF
  public :: ParentOceanOpenNetCDF
  public :: ParentOceanFinalizeNetCDF
  public :: ParentOceanInputNetCDF

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
  type :: vinfo
     character(len=32) :: name = ""
     logical  :: zstg = .false.
     logical  :: xstg = .false.
     logical  :: ystg = .false.
     real(RP) :: fact = 1.0_RP
     real(RP) :: offset = 0.0_RP
  end type vinfo

  type(hash_table) :: vars_atmos
  type(hash_table) :: vars_ocean
  type(hash_table) :: vars_land

  real(RP), allocatable :: work3d(:,:,:)
  real(RP), allocatable :: work2d(:,:)


  logical :: SCALE_tile_atm
  logical :: SCALE_tile_lnd
  logical :: SCALE_tile_ocn
  integer :: SCALE_DOMID_atm = -1
  integer :: SCALE_DOMID_lnd = -1
  integer :: SCALE_DOMID_ocn = -1
  integer :: nfiles_atm = 0
  integer :: nfiles_lnd = 0
  integer :: nfiles_ocn = 0
  integer :: fid_atm = -1
  integer :: fid_lnd = -1
  integer :: fid_ocn = -1
  integer, allocatable :: fids_atm(:)
  integer, allocatable :: fids_lnd(:)
  integer, allocatable :: fids_ocn(:)
  integer, allocatable :: tile_id_atm(:)
  integer, allocatable :: tile_id_lnd(:)
  integer, allocatable :: tile_id_ocn(:)

  integer, parameter :: vars_max = 100


  character(len=32) :: zname, zhname
  character(len=32) :: xname, xhname
  character(len=32) :: yname, yhname
  character(len=32) :: tname

  namelist / NetCDF_DIMS / &
       zname, &
       zhname, &
       xname, &
       xhname, &
       yname, &
       yhname, &
       tname

  character(len=32) :: item
  character(len=32) :: name
  logical :: zstg, xstg, ystg
  real(RP) :: fact, offset

  namelist / NetCDF_ITEM / &
       item, &
       name, &
       zstg, &
       xstg, &
       ystg, &
       fact, &
       offset


  logical :: first_atm
  logical :: first_ocn
  logical :: first_lnd


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Atmos Setup
  subroutine ParentAtmosSetupNetCDF( &
      dims,    &
      timelen, &
      mixing_ratio, &
      update_coord, &
      mapping_info, &
      qtrc_flag,    &
      lon_all,      &
      lat_all,      &
      basename_org, &
      basename_num, &
      same_mptype,  &
      pt_dry,       &
      serial,       &
      do_read       )
    use scale_prc, only: &
       PRC_isMaster
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       GRAV  => CONST_GRAV, &
       D2R   => CONST_D2R
    use scale_file, only: &
       FILE_get_dimLength, &
       FILE_get_attribute, &
       FILE_get_dataInfo
    use scale_mapprojection, only: &
       mappinginfo
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_domain_regist_file, &
       COMM_CARTESC_NEST_parent_info
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none
    integer,           intent(out) :: dims(6)
    integer,           intent(out) :: timelen
    logical,           intent(out) :: mixing_ratio
    logical,           intent(out) :: update_coord
    type(mappinginfo), intent(out) :: mapping_info
    logical,           intent(out) :: qtrc_flag(QA)
    real(RP), allocatable, intent(out) :: lon_all(:,:)
    real(RP), allocatable, intent(out) :: lat_all(:,:)

    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: basename_num
    logical,          intent(in) :: same_mptype

    logical, intent(inout) :: pt_dry
    logical, intent(inout) :: serial
    logical, intent(inout) :: do_read

    character(len=8) :: FILE_TYPE ! "SCALE-RM", "WRFARW", "NAMELIST", "AUTO"
    character(len=FILE_HLONG) :: NM_FILE ! namelist
    integer :: SCALE_PARENT_PRC_NUM_X
    integer :: SCALE_PARENT_PRC_NUM_Y
    character(len=FILE_HLONG) :: SCALE_LATLON_CATALOGUE

    namelist / PARAM_MKINIT_REAL_ATMOS_NetCDF / &
         FILE_TYPE, &
         NM_FILE, &
         MIXING_RATIO, &
         SCALE_PARENT_PRC_NUM_X, &
         SCALE_PARENT_PRC_NUM_Y, &
         SCALE_LATLON_CATALOGUE

    character(len=H_SHORT) :: mapping_name
    real(DP) :: false_easting
    real(DP) :: false_northing
    real(DP) :: longitude_of_central_meridian
    real(DP) :: longitude_of_projection_origin
    real(DP) :: latitude_of_projection_origin
    real(DP) :: straight_vertical_longitude_from_pole
    real(DP) :: standard_parallel(2)
    real(DP) :: rotation

    namelist / NetCDF_MAPPROJECTION / &
         mapping_name, &
         false_easting, &
         false_northing, &
         longitude_of_central_meridian, &
         longitude_of_projection_origin, &
         latitude_of_projection_origin, &
         straight_vertical_longitude_from_pole, &
         standard_parallel, &
         rotation

    character(len=32) :: items(vars_max)
    integer :: nvars
    type(vinfo), pointer :: var_info
    class(*), pointer :: v

    character(len=FILE_HLONG) :: basename
    character(len=FILE_HLONG) :: fname
    character(len=16) :: map

    integer :: nmfid
    integer :: i, n, iq
    integer :: ierr
    logical :: exist, error
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentAtmosSetupNetCDF",*) 'Real Case/Atmos Setup'

    FILE_TYPE = "AUTO"
    NM_FILE = ""
    MIXING_RATIO = .false.
    SCALE_PARENT_PRC_NUM_X = -1
    SCALE_PARENT_PRC_NUM_Y = -1
    SCALE_LATLON_CATALOGUE = ""

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_ATMOS_NetCDF,iostat=ierr)
    if( ierr > 0 ) then
       LOG_ERROR("ParentAtmosSetupNetCDF",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_ATMOS_NetCDF. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_ATMOS_NetCDF)

    basename = trim(basename_org) // trim(basename_num)

    fid_atm = -1
    if ( do_read ) then
       call check_filetype(fid_atm, FILE_TYPE, basename, SCALE_tile_atm, "ParentAtmosOpenNetCDF")
    end if

    call COMM_bcast( FILE_TYPE )

    if ( FILE_TYPE == "SCALE-RM" ) then
       call COMM_bcast( SCALE_tile_atm )
       if ( SCALE_tile_atm ) then
          do_read = .true.
          serial = .false.
       end if
    end if


    if ( do_read ) then

       mapping_name = ""
       false_easting = UNDEF
       false_easting = UNDEF
       longitude_of_central_meridian = UNDEF
       longitude_of_projection_origin = UNDEF
       latitude_of_projection_origin = UNDEF
       straight_vertical_longitude_from_pole = UNDEF
       standard_parallel = (/ UNDEF, UNDEF /)
       rotation = UNDEF

       vars_atmos = hash_table()

       select case( FILE_TYPE )
       case ( "SCALE-RM" )
          zname = "z"
          zhname = "zh"
          xname = "x"
          xhname = "xh"
          yname = "y"
          yhname = "yh"
          tname = "time"

          call vars_atmos%put("lon", vinfo("lon"))
          call vars_atmos%put("lat", vinfo("lat"))

          call vars_atmos%put("height",   vinfo("height"))
          call vars_atmos%put("pressure", vinfo("PRES"))

          call vars_atmos%put("DENS", vinfo("DENS"))
          call vars_atmos%put("W",    vinfo("W"))
          call vars_atmos%put("MOMZ", vinfo("MOMZ", zstg=.true.))
          call vars_atmos%put("Umet", vinfo("Umet"))
          call vars_atmos%put("U",    vinfo("U"))
          call vars_atmos%put("MOMX", vinfo("MOMX", xstg=.true.))
          call vars_atmos%put("Vmet", vinfo("Vmet"))
          call vars_atmos%put("V",    vinfo("V"))
          call vars_atmos%put("MOMY", vinfo("MOMY", ystg=.true.))

          call vars_atmos%put("PT", vinfo("PT"))
          call vars_atmos%put("T", vinfo("T"))
          call vars_atmos%put("RHOT", vinfo("RHOT"))

          call vars_atmos%put("QV", vinfo("QV"))
          call vars_atmos%put("RH", vinfo("RH"))

          if ( same_mptype ) then
             do iq = QS_MP, QE_MP
                call vars_atmos%put(TRACER_NAME(iq), vinfo(TRACER_NAME(iq)))
             end do
          else
             call vars_atmos%put("QC", vinfo("QC"))
             call vars_atmos%put("QR", vinfo("QR"))
             call vars_atmos%put("QI", vinfo("QI"))
             call vars_atmos%put("QS", vinfo("QS"))
             call vars_atmos%put("QG", vinfo("QG"))

             call vars_atmos%put("NC", vinfo("NC"))
             call vars_atmos%put("NR", vinfo("NR"))
             call vars_atmos%put("NI", vinfo("NI"))
             call vars_atmos%put("NS", vinfo("NS"))
             call vars_atmos%put("NG", vinfo("NG"))
          end if

          call vars_atmos%put("topo",     vinfo("topo"))
          call vars_atmos%put("MSLP",     vinfo("MSLP"))
          call vars_atmos%put("SFC_PRES", vinfo("SFC_PRES"))
          call vars_atmos%put("U10met",   vinfo("U10met"))
          call vars_atmos%put("U10",      vinfo("U10"))
          call vars_atmos%put("V10met",   vinfo("V10met"))
          call vars_atmos%put("V10",      vinfo("V10"))
          call vars_atmos%put("T2",       vinfo("T2"))
          call vars_atmos%put("Q2",       vinfo("Q2"))
          call vars_atmos%put("RH2",      vinfo("RH2"))

          mixing_ratio = .false.
          update_coord = .false.
          pt_dry       = .false.

          if ( PRC_isMaster ) then
             call FILE_get_attribute( fid_atm, "QV", "grid_mapping", map, existed=exist )
             if ( exist ) then
                call FILE_get_attribute( fid_atm, map, "grid_mapping_name", mapping_name )

                call FILE_get_attribute( fid_atm, map, "false_easting", false_easting, existed=exist )
                call FILE_get_attribute( fid_atm, map, "false_northing", false_northing, existed=exist )
                call FILE_get_attribute( fid_atm, map, "longitude_of_central_meridian", longitude_of_central_meridian, existed=exist )
                call FILE_get_attribute( fid_atm, map, "longitude_of_projection_origin", longitude_of_projection_origin, existed=exist )
                call FILE_get_attribute( fid_atm, map, "latitude_of_projection_origin", latitude_of_projection_origin, existed=exist )
                call FILE_get_attribute( fid_atm, map, "straight_vertical_longitude_from_pole", straight_vertical_longitude_from_pole, existed=exist )
                call FILE_get_attribute( fid_atm, map, "standard_parallel", standard_parallel(:), existed=exist )
                if ( .not. exist ) &
                call FILE_get_attribute( fid_atm, map, "standard_parallel", standard_parallel(1:1), existed=exist )
                call FILE_get_attribute( fid_atm, map, "rotation", rotation, existed=exist )
             end if
          end if

          call COMM_bcast( mapping_name )

          call COMM_bcast( false_easting )
          call COMM_bcast( false_northing )
          call COMM_bcast( longitude_of_central_meridian )
          call COMM_bcast( longitude_of_projection_origin )
          call COMM_bcast( latitude_of_projection_origin )
          call COMM_bcast( straight_vertical_longitude_from_pole )
          call COMM_bcast( 2, standard_parallel )
          call COMM_bcast( rotation )

       case ( "WRFARW" )
          zname = "bottom_top"
          zhname = "bottom_top_stag"
          xname = "west_east"
          xhname = "west_east_stag"
          yname = "south_north"
          yhname = "south_north_stag"
          tname = "Time"

          call vars_atmos%put("lon", vinfo("XLONG"))
          call vars_atmos%put("lat", vinfo("XLAT"))

          call vars_atmos%put("hbar", vinfo("PHB", zstg=.true., fact=1.0_RP/GRAV)) ! geopotential height
          call vars_atmos%put("hdev", vinfo("PH", zstg=.true., fact=1.0_RP/GRAV))

          call vars_atmos%put("pbar", vinfo("PB"))
          call vars_atmos%put("pdev", vinfo("P"))

          call FILE_get_dataInfo( fid_atm, "U", existed=exist )
          if ( exist ) then
             LOG_INFO("ParentAtmosSetupNetCDF",*) 'WRF-ARW FILE-TYPE: WRF History Output'
             call vars_atmos%put("W", vinfo("W",zstg=.true.))
             call vars_atmos%put("U", vinfo("U",xstg=.true.))
             call vars_atmos%put("V", vinfo("V",ystg=.true.))
             call vars_atmos%put("PT", vinfo("T", offset=300.0_RP))
          else
             LOG_INFO("ParentAtmosSetupNetCDF",*) 'WRF-ARW FILE-TYPE: WRF Restart'
             call vars_atmos%put("W", vinfo("W_1"))
             call vars_atmos%put("U", vinfo("U_1"))
             call vars_atmos%put("V", vinfo("V_1"))
             call vars_atmos%put("PT", vinfo("T_1", offset=300.0_RP))
          endif

          if ( same_mptype ) then
             LOG_ERROR("ParentAtmosSetupNetCDF",*) 'same_mptype must be .false. for WRF file'
             call PRC_abort
          end if
          call vars_atmos%put("QV", vinfo("QVAPOR"))
          call vars_atmos%put("QC", vinfo("QCLOUD"))
          call vars_atmos%put("QR", vinfo("QRAIN"))
          call vars_atmos%put("QI", vinfo("QICE"))
          call vars_atmos%put("QS", vinfo("QSNOW"))
          call vars_atmos%put("QG", vinfo("QGRAUP"))
          call vars_atmos%put("NC", vinfo("NC"))
          call vars_atmos%put("NR", vinfo("NR"))
          call vars_atmos%put("NI", vinfo("NI"))
          call vars_atmos%put("NS", vinfo("NS"))
          call vars_atmos%put("NG", vinfo("NG"))

          mixing_ratio = .true.
          pt_dry       = .true.

          call vars_atmos%put("topo", vinfo("HGT"))
          call vars_atmos%put("U10",    vinfo("U10"))
          call vars_atmos%put("V10",    vinfo("V10"))
          call vars_atmos%put("T2",     vinfo("T2"))
          call vars_atmos%put("Q2",     vinfo("Q2"))
          call vars_atmos%put("RH2",    vinfo("RH2"))
          call vars_atmos%put("SFC_PRES",    vinfo("PSFC"))

          call FILE_get_attribute( fid_atm, "global", "MAP_PROJ", i, existed=exist )
          if ( exist ) then
             if ( i == 1 ) then ! Lambert Conformal
                mapping_name = "lambert_conformal_conic"
                call FILE_get_attribute( fid_atm, "global", "TRUELAT1", standard_parallel(1) )
                call FILE_get_attribute( fid_atm, "global", "TRUELAT2", standard_parallel(2) )
                call FILE_get_attribute( fid_atm, "global", "STAND_LON", longitude_of_central_meridian )
             else if ( i >= 3 ) then ! No rotate
                ! do nothing
             else
                LOG_WARN("ParentAtmodSetupNetCDF",*) "This map projection type is not supported: ", i
                LOG_WARN_CONT(*) "Specify map projection parameter manually"
             end if
          end if

          update_coord = .true.

       case ( "NAMELIST" )

          update_coord = .true.

       case default
          LOG_ERROR("ParentAtmosSetupNetCDF",*) 'FILE_TYPE must be "SCALE-RM", "WRFARW", "NAMELIST", or "AUTO", ', trim(FILE_TYPE)
          call PRC_abort
       end select

       !--- read namelist
       nmfid = -1
       if ( NM_FILE /= "" ) then
          nmfid = IO_get_available_fid()
          call IO_get_fname(fname, NM_FILE)
          open(nmfid, file=fname, form="formatted", status="old", action="read", iostat=ierr)
          if ( ierr /= 0 ) then
             LOG_ERROR("ParentAtmosSetupNetCDF",*) 'namelist file is not found! ', trim(fname)
             call PRC_abort
          end if

          read(nmfid, nml=NetCDF_DIMS, iostat=ierr)
          if( ierr > 0 ) then
             LOG_ERROR("ParentAtmosSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_DIMS in ', trim(fname), '. Check!'
             call PRC_abort
          end if

          rewind(nmfid)
          read(nmfid, nml=NetCDF_MAPPROJECTION, iostat=ierr)
          if( ierr > 0 ) then
             LOG_ERROR("ParentAtmosSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_MAPPROJECTION in ', trim(fname), '. Check!'
             call PRC_abort
          end if
          ! items
          rewind(nmfid)
          nvars = 0
          do n = 1, vars_max
             read(nmfid, nml=NetCDF_ITEM, iostat=ierr)
             if( ierr > 0 ) then
                LOG_ERROR("ParentAtmosSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_ITEM in ', trim(fname), '. Check!'
                call PRC_abort
             else if( ierr < 0 ) then
                exit
             end if
             nvars = nvars + 1
             items(nvars) = item
          end do
          if ( nvars > vars_max ) then
             LOG_ERROR("ParentAtmosSetupNetCDF",*) "The number of item in the namelist file exceeds the limit! ", nvars
             call PRC_abort
          end if
          rewind(nmfid)
          do n = 1, nvars
             ! set default
             if ( vars_atmos%has_key(items(n)) ) then
                item = items(n)
                v => vars_atmos%get(item)
                select type( v )
                type is (vinfo)
                   var_info => v
                end select
                name = var_info%name
                zstg = var_info%zstg
                xstg = var_info%xstg
                ystg = var_info%ystg
                fact = var_info%fact
                offset = var_info%offset
             else
                item = items(n)
                name = items(n)
                zstg = .false.
                xstg = .false.
                ystg = .false.
                fact = 1.0_RP
                offset = 0.0_RP
             end if
             read(nmfid, nml=NetCDF_ITEM, iostat=ierr)
             if ( ierr /= 0 ) exit
             ! set params
             call vars_atmos%put(item, vinfo(name=name, zstg=zstg, xstg=xstg, ystg=ystg, fact=fact, offset=offset))
          end do

       else if ( FILE_TYPE == "NAMELIST" ) then
          LOG_ERROR("ParentAtmosSetupNetCDF",*) 'NM_FILE is necessary'
          call PRC_abort
       end if

       mapping_info%mapping_name = mapping_name
       if ( false_easting /= UNDEF ) mapping_info%false_easting = false_easting
       if ( false_northing /= UNDEF ) mapping_info%false_northing = false_northing
       if ( longitude_of_central_meridian /= UNDEF ) mapping_info%longitude_of_central_meridian = longitude_of_central_meridian
       if ( longitude_of_projection_origin /= UNDEF ) mapping_info%longitude_of_projection_origin = longitude_of_projection_origin
       if ( latitude_of_projection_origin /= UNDEF ) mapping_info%latitude_of_projection_origin = latitude_of_projection_origin
       if ( straight_vertical_longitude_from_pole /= UNDEF ) mapping_info%straight_vertical_longitude_from_pole = straight_vertical_longitude_from_pole
       if ( standard_parallel(1) /= UNDEF ) mapping_info%standard_parallel(1) = standard_parallel(1)
       if ( standard_parallel(2) /= UNDEF ) mapping_info%standard_parallel(2) = standard_parallel(2)
       if ( rotation /= UNDEF ) mapping_info%rotation = rotation

    end if

    if ( SCALE_tile_atm ) then

       call COMM_CARTESC_NEST_domain_regist_file( &
            SCALE_DOMID_atm,        & ! [OUT]
            basename,               & ! [IN]
            SCALE_PARENT_PRC_NUM_X, & ! [IN]
            SCALE_PARENT_PRC_NUM_Y, & ! [IN]
            SCALE_LATLON_CATALOGUE  ) ! [IN]

       call COMM_CARTESC_NEST_parent_info( &
            SCALE_DOMID_atm,     & ! [IN]
            KMAX=dims(1),        & ! [OUT]
            IMAXG=dims(2),       & ! [OUT]
            JMAXG=dims(3),       & ! [OUT]
            num_tile=nfiles_atm  ) ! [OUT]

       dims(4) = dims(1)
       dims(5) = dims(2)
       dims(6) = dims(3)

       allocate( fids_atm(nfiles_atm) )
       allocate( tile_id_atm(nfiles_atm) )

       call COMM_CARTESC_NEST_parent_info( &
            SCALE_DOMID_atm,      & ! [IN]
            tile_id = tile_id_atm ) ! [OUT]

       call ParentAtmosOpenNetCDF( basename_org, basename_num )

    else if ( do_read ) then

       call FILE_get_dimLength( fid_atm, zname,  dims(1) )
       call FILE_get_dimLength( fid_atm, xname,  dims(2) )
       call FILE_get_dimLength( fid_atm, yname,  dims(3) )
       call FILE_get_dimLength( fid_atm, zhname, dims(4) )
       call FILE_get_dimLength( fid_atm, xhname, dims(5) )
       call FILE_get_dimLength( fid_atm, yhname, dims(6) )

    end if

    if ( do_read ) then

       do iq = 1, QA
          if ( iq >= QS_MP .and. iq <= QE_MP ) cycle
          qtrc_flag(iq) = .false.
          if ( vars_atmos%has_key( TRACER_NAME(iq) ) ) then
             select type ( v => vars_atmos%get( TRACER_NAME(iq) ) )
             type is ( vinfo )
                if ( v%name .ne. "" ) then
                   call FILE_get_dataInfo( fid_atm, v%name, existed = qtrc_flag(iq) )
                end if
             end select
          end if
       end do

       call FILE_get_dimLength( fid_atm, tname, timelen, error=error )
       if ( error ) timelen = 1

       allocate( lon_all(dims(2), dims(3)) )
       allocate( lat_all(dims(2), dims(3)) )

       call read2d( dims(2), 1, dims(2), dims(3), 1, dims(3), &
            LON_all(:,:), vars_atmos%get("lon"), &
            1, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm )
       lon_all(:,:) = lon_all(:,:) * D2R
       call read2d( dims(2), 1, dims(2), dims(3), 1, dims(3), &
            LAT_all(:,:), vars_atmos%get("lat"), &
            1, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm )
       lat_all(:,:) = lat_all(:,:) * D2R

    end if


    first_atm = .true.

    return
  end subroutine ParentAtmosSetupNetCDF

  !-----------------------------------------------------------------------------
  !> Atmos Open
  subroutine ParentAtmosOpenNetCDF( &
       basename_org, basename_num )
    use scale_file, only: &
       FILE_open
    implicit none
    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: basename_num

    character(len=FILE_HLONG) :: basename

    integer :: n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentAtmosOpenNetCDF",*) 'Real Case/Atmos Open'


    basename = trim(basename_org) // trim(basename_num)

    if ( SCALE_tile_atm ) then

       do n = 1, nfiles_atm
          call FILE_open( &
               basename,              & ! [IN]
               fids_atm(n),           & ! [OUT]
               aggregate=.false.,     & ! [IN]
               rankid=tile_id_atm(n)  ) ! [IN]
       end do

       fid_atm = fids_atm(1)

    else

       call FILE_open(basename, fid_atm, postfix="")

    end if

    return
  end subroutine ParentAtmosOpenNetCDF

  !-----------------------------------------------------------------------------
  !> Atmos Finalize
  subroutine ParentAtmosFinalizeNetCDF
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentAtmosFinalizeNetCDF",*) 'Real Case/Atmos Finalize'

    deallocate( work2d )
    deallocate( work3d )
    if ( allocated(fids_atm)    ) deallocate( fids_atm )
    if ( allocated(tile_id_atm) ) deallocate( tile_id_atm )

    call vars_atmos%destroy()

    SCALE_DOMID_atm = -1
    nfiles_atm = 0
    fid_atm = -1

    return
  end subroutine ParentAtmosFinalizeNetCDF

  !-----------------------------------------------------------------------------
  subroutine ParentAtmosInputNetCDF( &
       KA_org, KS_org, KE_org, &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       QA, &
       cz_org,              &
       w_org, u_org, v_org, &
       pres_org,            &
       dens_org,            &
       temp_org, pt_org,    &
       qtrc_org,            &
       qv_org, rh_org,      &
       qhyd_org, qnum_org,  &
       nopres, nodens,      &
       uvmet,               &
       temp2pt, rh2qv,      &
       qnum_flag,           &
       same_mptype ,        &
       sfc_diagnoses,       &
       update_coord,        &
       dims,                &
       it                   )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       HYD_NAME, &
       NUM_NAME, &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none
    integer, intent(in) :: KA_org, KS_org, KE_org
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org
    integer, intent(in) :: QA

    real(RP), intent(inout) :: cz_org(KA_org,IA_org,JA_org)

    real(RP), intent(out) :: w_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: u_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: v_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: pres_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: dens_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: temp_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: pt_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: qtrc_org(KA_org,IA_org,JA_org,QA)
    real(RP), intent(out) :: qv_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: rh_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: qhyd_org(KA_org,IA_org,JA_org,N_HYD)
    real(RP), intent(out) :: qnum_org(KA_org,IA_org,JA_org,N_HYD)
    logical,  intent(out) :: nopres
    logical,  intent(out) :: nodens
    logical,  intent(out) :: uvmet
    logical,  intent(out) :: temp2pt
    logical,  intent(out) :: rh2qv
    logical,  intent(out) :: qnum_flag

    logical, intent(in)  :: same_mptype
    logical, intent(in)  :: sfc_diagnoses
    logical, intent(in)  :: update_coord
    integer, intent(in)  :: dims(6)
    integer, intent(in)  :: it


    logical :: exist
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( .not. allocated( work2d ) ) then
       allocate( work2d(IA_org,JA_org) )
       allocate( work3d(KA_org-2,IA_org,JA_org) )
    end if

    ! height
    if ( first_atm .or. update_coord ) then
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            cz_org(:,:,:), vars_atmos%get("height"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( .not. exist ) then
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               cz_org(:,:,:), vars_atmos%get("hbar"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             call read3d( KA_org-2, KS_org, KE_org-2, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
                  work3d(:,:,:), vars_atmos%get("hdev"), &
                  it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
                  exist = exist )
          end if
          if ( .not. exist ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) '"height" or "hbar"+"hdev" is necessary'
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             cz_org(k+2,i,j) = cz_org(k+2,i,j) + work3d(k,i,j)
          end do
          end do
          end do
       else
       end if
    end if

    ! tracers other than water contents
    do iq = 1, QA
       if ( iq >= QS_MP .and. iq <= QE_MP ) cycle
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            qtrc_org(:,:,:,iq), vars_atmos%get(TRACER_NAME(iq)), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( .not. exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             qtrc_org(k+2,i,j,iq) = UNDEF
          end do
          end do
          end do
       end if
    end do

    ! water contents
    if ( same_mptype ) then

       do iq = QS_MP, QE_MP
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               qtrc_org(:,:,:,iq), vars_atmos%get(TRACER_NAME(iq)), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( .not. exist ) then
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 1, KA_org-2
                qtrc_org(k+2,i,j,iq) = UNDEF
             end do
             end do
             end do
          end if
       end do

    else

       ! qv
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            qv_org(:,:,:), vars_atmos%get("QV"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          rh2qv = .false.
       else
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               rh_org(:,:,:), vars_atmos%get("RH"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm )
          rh2qv = .true.
       end if

       ! qhyd
       qnum_flag = .false.
       do iq = 1, N_HYD
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               qhyd_org(:,:,:,iq), vars_atmos%get(HYD_NAME(iq)), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( .not. exist ) then
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 1, KA_org-2
                qhyd_org(k+2,i,j,iq) = 0.0_RP
             end do
             end do
             end do
          end if

          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               qnum_org(:,:,:,iq), vars_atmos%get(NUM_NAME(iq)), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             qnum_flag = .true.
          else
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = 1, KA_org-2
                qnum_org(k+2,i,j,iq) = UNDEF
             end do
             end do
             end do
          end if
       end do

    end if

    ! pressure
    call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         pres_org(:,:,:), vars_atmos%get("pressure"), &
         it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
         exist = exist )
    if ( .not. exist ) then
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            pres_org(:,:,:), vars_atmos%get("pbar"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          call read3d( KA_org-2, KS_org, KE_org-2, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work3d(:,:,:), vars_atmos%get("pdev"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( .not. exist ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) '"pdev" is necessary if "pbar" exists'
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             pres_org(k+2,i,j) = pres_org(k+2,i,j) + work3d(k,i,j)
          end do
          end do
          end do
       end if
    end if
    if ( exist ) then
       nopres = .false.
    else
       nopres = .true.
    end if

    ! density
    call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         dens_org(:,:,:), vars_atmos%get("DENS"), &
         it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
         exist = exist )
    if ( exist ) then
       nodens = .false.
    else
       nodens = .true.
    end if

    ! pt
    call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         pt_org(:,:,:), vars_atmos%get("PT"), &
         it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
         exist = exist )
    if ( exist ) then
       temp2pt = .false.
    else
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            pt_org(:,:,:), vars_atmos%get("RHOT"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          if ( nodens ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) "DENS is necessary to calculate PT from RHOT"
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             pt_org(k+2,i,j) = pt_org(k+2,i,j) / dens_org(k+2,i,j)
          end do
          end do
          end do
          temp2pt = .false.
       else
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               temp_org(:,:,:), vars_atmos%get("T"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( .not. exist ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) '"PT", "RHOT", or "T" is necessary'
             call PRC_abort
          end if
          temp2pt = .true.
       end if
    end if

    ! W
    call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         w_org(:,:,:), vars_atmos%get("W"), &
         it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
         exist = exist )
    if ( .not. exist ) then
       if ( nodens ) then
          LOG_ERROR("ParentAtmosInputNetCDF",*) "DENS is necessary to use MOMZ"
          call PRC_abort
       end if
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            w_org(:,:,:), vars_atmos%get("MOMZ"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             w_org(k+2,i,j) =  w_org(k+2,i,j) / dens_org(k+2,i,j)
          end do
          end do
          end do
       else
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             w_org(k+2,i,j) =  0.0_RP
          end do
          end do
          end do
       endif
    end if

    ! U
    call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         u_org(:,:,:), vars_atmos%get("Umet"), &
         it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
         exist = exist )
    if ( exist ) then
       uvmet = .true.
    else
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            u_org(:,:,:), vars_atmos%get("U"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( .not. exist ) then
          if ( nodens ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) "DENS is necessary to use MOMX"
             call PRC_abort
          end if
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               u_org(:,:,:), vars_atmos%get("MOMX"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( .not. exist ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) '"Ument", "U", or "MOMX" is necessary'
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             u_org(k+2,i,j) = u_org(k+2,i,j) / dens_org(k+2,i,j)
          end do
          end do
          end do
       end if
       uvmet = .false.
    end if

    ! V
    if ( uvmet ) then
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            v_org(:,:,:), vars_atmos%get("Vmet"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( .not. exist ) then
          LOG_ERROR("ParentAtmosInputNetCDF",*) "Vmet is required when Umet exists"
          call PRC_abort
       end if
    else
       call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            v_org(:,:,:), vars_atmos%get("V"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( .not. exist ) then
          if ( nodens ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) "DENS is necessary to use MOMY"
             call PRC_abort
          end if
          call read3d( KA_org, KS_org+2, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               v_org(:,:,:), vars_atmos%get("MOMY"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( .not. exist ) then
             LOG_ERROR("ParentAtmosInputNetCDF",*) '"V" or "MOMY" is required when "U" or "MOMX" exists'
             call PRC_abort
          end if
          !$omp parallel do collapse(2)
          do j = 1, JA_org
          do i = 1, IA_org
          do k = 1, KA_org-2
             v_org(k+2,i,j) = v_org(k+2,i,j) / dens_org(k+2,i,j)
          end do
          end do
          end do
       end if
    end if


    if ( sfc_diagnoses ) then

       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          cz_org(1,i,j) = 0.0_RP
       end do
       end do


       if ( first_atm .or. update_coord ) then
          ! topo
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("topo"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                cz_org(2,i,j) = work2d(i,j)
             end do
             end do
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                cz_org(2,i,j) = UNDEF
             end do
             end do
          end if
       end if

       ! MSLP
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            work2d(:,:), vars_atmos%get("MSLP"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(1,i,j) = work2d(i,j)
          end do
          end do
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(1,i,j) = UNDEF
          end do
          end do
       end if
       
       ! SFC_PRES
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            work2d(:,:), vars_atmos%get("SFC_PRES"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(2,i,j) = work2d(i,j)
          end do
          end do
       else
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             pres_org(2,i,j) = UNDEF
          end do
          end do
       end if

       ! U10, V10
       if ( uvmet ) then
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("U10met"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                U_org(2,i,j) = work2d(i,j)
             end do
             end do
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                U_org(2,i,j) = UNDEF
             end do
             end do
          end if
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("V10met"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                V_org(2,i,j) = work2d(i,j)
             end do
             end do
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                V_org(2,i,j) = UNDEF
             end do
             end do
          end if
       else
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("U10"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                U_org(2,i,j) = work2d(i,j)
             end do
             end do
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                U_org(2,i,j) = UNDEF
             end do
             end do
          end if
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("V10"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                V_org(2,i,j) = work2d(i,j)
             end do
             end do
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                V_org(2,i,j) = UNDEF
             end do
             end do
          end if
       end if

       ! T2
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            work2d(:,:), vars_atmos%get("T2"), &
            it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
            exist = exist )
       if ( exist ) then
          !$omp parallel do
          do j = 1, JA_org
          do i = 1, IA_org
             temp_org(2,i,j) = work2d(i,j)
          end do
          end do
       end if

       ! Q2
       if ( rh2qv ) then
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("RH2"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                rh_org(2,i,j) = work2d(i,j)
             end do
             end do
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                rh_org(2,i,j) = UNDEF
             end do
             end do
          end if
       else
          call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
               work2d(:,:), vars_atmos%get("Q2"), &
               it, nfiles_atm, fid_atm, fids_atm, SCALE_tile_atm, SCALE_DOMID_atm, &
               exist = exist )
          if ( exist ) then
             if ( same_mptype ) then
                !$omp parallel do
                do j = 1, JA_org
                do i = 1, IA_org
                   qtrc_org(2,i,j,QS_MP) = work2d(i,j)
                end do
                end do
             else
                !$omp parallel do
                do j = 1, JA_org
                do i = 1, IA_org
                   qv_org(2,i,j) = work2d(i,j)
                end do
                end do
             end if
          else
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                qv_org(2,i,j) = UNDEF
             end do
             end do
          end if
       end if


    else
       
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          cz_org  (1,i,j) = 0.0_RP
          cz_org  (2,i,j) = 0.0_RP
          pres_org(1,i,j) = UNDEF
          pres_org(2,i,j) = UNDEF
          u_org   (2,i,j) = UNDEF
          v_org   (2,i,j) = UNDEF
          temp_org(2,i,j) = UNDEF
          pt_org  (2,i,j) = UNDEF
          qv_org  (2,i,j) = UNDEF
          rh_org  (2,i,j) = UNDEF
       end do
       end do

    end if

    first_atm = .false.

    return
  end subroutine ParentAtmosInputNetCDF

  !-----------------------------------------------------------------------------
  !> Land Setup
  subroutine ParentLandSetupNetCDF( &
       ldims,        &
       timelen,      &
       lon_all,      &
       lat_all,      &
       basename_org, &
       basename_num, &
       use_file_landwater, &
       serial, &
       do_read )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_file, only: &
       FILE_get_dimLength
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_domain_regist_file, &
       COMM_CARTESC_NEST_parent_info

    implicit none

    integer,           intent(out) :: ldims(3)
    integer,           intent(out) :: timelen
    real(RP), allocatable, intent(out) :: lon_all(:,:)
    real(RP), allocatable, intent(out) :: lat_all(:,:)

    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: basename_num
    logical,          intent(in) :: use_file_landwater

    logical, intent(inout) :: serial
    logical, intent(inout) :: do_read

    character(len=8) :: FILE_TYPE = "AUTO" ! "SCALE-RM", "WRFARW", "NAMELIST", "AUTO"
    character(len=FILE_HLONG) :: NM_FILE
    logical :: SCALE_MULTI_FILE = .true.
    integer :: SCALE_PARENT_PRC_NUM_X
    integer :: SCALE_PARENT_PRC_NUM_Y
    character(len=FILE_HLONG) :: SCALE_LATLON_CATALOGUE

    namelist / PARAM_MKINIT_REAL_LAND_NetCDF / &
         FILE_TYPE, &
         NM_FILE, &
         SCALE_PARENT_PRC_NUM_X, &
         SCALE_PARENT_PRC_NUM_Y, &
         SCALE_LATLON_CATALOGUE

    character(len=FILE_HLONG) :: basename
    character(len=FILE_HLONG) :: fname
    integer :: nmfid

    character(len=32) :: items(vars_max)
    integer :: nvars
    type(vinfo), pointer :: var_info
    class(*), pointer :: v

    logical :: error, exist
    integer :: n, i
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_INFO("ParentLandSetupNetCDF",*) 'Real Case/Land Setup'

    FILE_TYPE = "AUTO"
    NM_FILE = ""
    SCALE_PARENT_PRC_NUM_X = -1
    SCALE_PARENT_PRC_NUM_Y = -1
    SCALE_LATLON_CATALOGUE = ""

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_LAND_NetCDF,iostat=ierr)
    if( ierr > 0 ) then
       LOG_ERROR("ParentLandSetupNetCDF",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_LAND_NetCDF. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_LAND_NetCDF)


    basename = trim(basename_org) // trim(basename_num)

    fid_lnd = -1
    if ( do_read ) then
       call check_filetype(fid_lnd, FILE_TYPE, basename, SCALE_tile_lnd, "ParentLandSetupNetCDF")
    end if

    call COMM_bcast( FILE_TYPE )

    if ( FILE_TYPE == "SCALE-RM" ) then
       call COMM_bcast( SCALE_tile_lnd )
       if ( SCALE_tile_lnd ) then
          do_read = .true.
          serial = .false.
       end if
    end if

    if ( do_read ) then
       vars_land = hash_table()

       select case( FILE_TYPE )
       case ( "SCALE-RM" )
          zname = "lz"
          xname = "x"
          yname = "y"
          tname = "time"

          call vars_land%put("lon", vinfo("lon"))
          call vars_land%put("lat", vinfo("lat"))
          call vars_land%put("lz",  vinfo("lz"))

          call vars_land%put("topo", vinfo("topo"))
          call vars_land%put("lsmask", vinfo("lsmask"))

          call vars_land%put("LAND_TEMP", vinfo("LAND_TEMP"))
          if ( use_file_landwater ) then
             call vars_land%put("LAND_WATER", vinfo("LAND_WATER"))
          end if

          call vars_land%put("LAND_SFC_TEMP",  vinfo("LAND_SFC_TEMP"))

          call vars_land%put("LAND_SFC_ALB_IR_dir",  vinfo("LAND_SFC_ALB_IR_dir"))
          call vars_land%put("LAND_SFC_ALB_IR_dif",  vinfo("LAND_SFC_ALB_IR_dif"))
          call vars_land%put("LAND_SFC_ALB_NIR_dir", vinfo("LAND_SFC_ALB_NIR_dir"))
          call vars_land%put("LAND_SFC_ALB_NIR_dif", vinfo("LAND_SFC_ALB_NIR_dif"))
          call vars_land%put("LAND_SFC_ALB_VIS_dir", vinfo("LAND_SFC_ALB_VIS_dir"))
          call vars_land%put("LAND_SFC_ALB_VIS_dif", vinfo("LAND_SFC_ALB_VIS_dif"))

          call vars_land%put("URBAN_SFC_TEMP",  vinfo("URBAN_SFC_TEMP"))

       case ( "WRFARW" )
          zname = "soil_layers_stag"
          xname = "west_east"
          yname = "south_north"
          tname = "Time"

          call vars_land%put("lon", vinfo("XLONG"))
          call vars_land%put("lat", vinfo("XLAT"))
          call vars_land%put("lz",  vinfo("ZS"))

          call vars_land%put("topo",   vinfo("HGT"))
          call vars_land%put("lsmask", vinfo("LANDMASK"))

          call vars_land%put("LAND_TEMP", vinfo("TSLB"))
          if ( use_file_landwater ) then
             call vars_land%put("LAND_WATER", vinfo("SH2O"))
          end if

          call vars_land%put("LAND_SFC_TEMP", vinfo("TSK"))

          call vars_land%put("LAND_SFC_ALB_VIS_dir", vinfo("ALBEDO"))
          call vars_land%put("LAND_SFC_EMIS_IR_dif", vinfo("EMISS"))

          call vars_land%put("URBAN_SFC_TEMP", vinfo("URBAN_SFC_TEMP"))

!       call vars_land%put("SNOW_WATER", vinfo("SNOW"))
!       call vars_land%put("SNOW_TEMP",  vinfo("TSNAV"))

       case ( "NAMELIST" )
       case default
          LOG_ERROR("ParentLandSetupNetCDF",*) 'FILE_TYPE must be "SCALE-RM", "WRFARW", or "AUTO", ', trim(FILE_TYPE)
          call PRC_abort
       end select


       !--- read namelist
       nmfid = -1
       if ( NM_FILE /= "" ) then
          nmfid = IO_get_available_fid()
          call IO_get_fname(fname, NM_FILE)
          open(nmfid, file=fname, form="formatted", status="old", action="read", iostat=ierr)
          if ( ierr /= 0 ) then
             LOG_ERROR("ParentLandSetupNetCDF",*) 'namelist file is not found! ', trim(fname)
             call PRC_abort
          end if

          rewind(nmfid)
          read(nmfid, nml=NetCDF_DIMS, iostat=ierr)
          if( ierr > 0 ) then
             LOG_ERROR("ParentLandSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_DIMS in ', trim(fname), '. Check!'
             call PRC_abort
          end if

          ! items
          rewind(nmfid)
          nvars = 0
          do n = 1, vars_max
             read(nmfid, nml=NetCDF_ITEM, iostat=ierr)
             if( ierr > 0 ) then
                LOG_ERROR("ParentLandSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_ITEM in ', trim(fname), '. Check!'
                call PRC_abort
             else if( ierr < 0 ) then
                exit
             end if
             nvars = nvars + 1
             items(nvars) = item
          end do
          if ( nvars > vars_max ) then
             LOG_ERROR("ParentLandSetupNetCDF",*) "The number of item in the namelist file exceeds the limit! ", nvars
             call PRC_abort
          end if
          rewind(nmfid)
          do n = 1, nvars
             ! set default
             if ( vars_land%has_key(items(n)) ) then
                item = items(n)
                v => vars_land%get(item)
                select type(v)
                type is (vinfo)
                   var_info => v
                end select
                name = var_info%name
                fact = var_info%fact
                offset = var_info%offset
             else
                item = items(n)
                name = items(n)
                fact = 1.0_RP
                offset = 0.0_RP
             end if
             zstg = .false.
             xstg = .false.
             ystg = .false.
             read(nmfid, nml=NetCDF_ITEM, iostat=ierr)
             if ( ierr /= 0 ) exit
             ! set params
             call vars_land%put(item, vinfo(name=name, zstg=zstg, xstg=xstg, ystg=ystg, fact=fact, offset=offset))
          end do

       else if ( FILE_TYPE == "NAMELIST" ) then
          LOG_ERROR("ParentLANDSetupNetCDF",*) 'NM_FILE is necessary'
          call PRC_abort
       end if
       
    end if

    if ( SCALE_tile_lnd ) then

       call COMM_CARTESC_NEST_domain_regist_file( &
            SCALE_DOMID_lnd,        & ! [OUT]
            basename,               & ! [IN]
            SCALE_PARENT_PRC_NUM_X, & ! [IN]
            SCALE_PARENT_PRC_NUM_Y, & ! [IN]
            SCALE_LATLON_CATALOGUE  ) ! [IN]

       call COMM_CARTESC_NEST_parent_info( &
            SCALE_DOMID_lnd,    & ! [IN]
            LKMAX=ldims(1),     & ! [OUT]
            IMAXG=ldims(2),     & ! [OUT]
            JMAXG=ldims(3),     & ! [OUT]
            num_tile=nfiles_lnd ) ! [OUT]

       allocate( fids_lnd   (nfiles_lnd) )
       allocate( tile_id_lnd(nfiles_lnd) )

       call COMM_CARTESC_NEST_parent_info( &
            SCALE_DOMID_lnd,      & ! [IN]
            tile_id = tile_id_lnd ) ! [OUT]

       call ParentLandOpenNetCDF( basename_org, basename_num )

    else if ( do_read ) then

       call FILE_get_dimLength( fid_lnd, zname, ldims(1) )
       call FILE_get_dimLength( fid_lnd, xname, ldims(2) )
       call FILE_get_dimLength( fid_lnd, yname, ldims(3) )

    end if

    if ( do_read ) then

       call FILE_get_dimLength( fid_lnd, tname, timelen, error=error )
       if ( error ) timelen = 1

       allocate( lon_all(ldims(2), ldims(3)) )
       allocate( lat_all(ldims(2), ldims(3)) )

       call read2d( ldims(2), 1, ldims(2), ldims(3), 1, ldims(3), &
            lon_all(:,:), vars_land%get("lon"), &
            1, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )
       lon_all(:,:) = lon_all(:,:) * D2R
       call read2d( ldims(2), 1, ldims(2), ldims(3), 1, ldims(3), &
            lat_all(:,:), vars_land%get("lat"), &
            1, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )
       lat_all(:,:) = lat_all(:,:) * D2R

    end if

    first_lnd = .true.

    return
  end subroutine ParentLandSetupNetCDF

  !-----------------------------------------------------------------------------
  !> Land Open
  subroutine ParentLandOpenNetCDF( &
       basename_org, basename_num )
    use scale_file, only: &
       FILE_open
    implicit none
    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: basename_num

    character(len=FILE_HLONG)  :: basename

    integer :: n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentLandOpenNetCDF",*) 'Real Case/Land Open'

    basename = trim(basename_org) // trim(basename_num)

    if ( SCALE_tile_lnd ) then

       do n = 1, nfiles_lnd
          call FILE_open( &
               basename,             & ! [IN]
               fids_lnd(n),          & ! [OUT]
               aggregate=.false.,    & ! [IN]
               rankid=tile_id_lnd(n) ) ! [IN]
       end do

       fid_lnd = fids_lnd(1)

    else

       call FILE_open(basename, fid_lnd, postfix="")

    end if

    return
  end subroutine ParentLandOpenNetCDF

  !-----------------------------------------------------------------------------
  !> Land Finalize
  subroutine ParentLandFinalizeNetCDF
    implicit none

    !---------------------------------------------------------------------------

    LOG_INFO("ParentLandFinalizeNetCDF",*) 'Real Case/Land Finalize'

    if ( allocated(fids_lnd)    ) deallocate( fids_lnd )
    if ( allocated(tile_id_lnd) ) deallocate( tile_id_lnd )

    call vars_land%destroy()

    SCALE_DOMID_lnd = -1
    nfiles_lnd = 0
    fid_lnd = -1

    return
  end subroutine ParentLandFinalizeNetCDF

  !-----------------------------------------------------------------------------
  subroutine ParentLandInputNetCDF( &
       KA_org, KS_org, KE_org, &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       tg_org,             &
       strg_org,           &
       lst_org,            &
       ust_org,            &
       albg_org,           &
       topo_org,           &
       lmask_org,          &
       lz_org,             &
       use_file_landwater, &
       ldims,              &
       it                  )
    use scale_const, only: &
       D2R => CONST_D2R, &
       UNDEF => CONST_UNDEF
    use scale_file, only: &
       FILE_open, &
       FILE_read
    implicit none
    integer, intent(in) :: KA_org, KS_org, KE_org
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org

    real(RP), intent(out) :: tg_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: strg_org(KA_org,IA_org,JA_org)
    real(RP), intent(out) :: lst_org(IA_org,JA_org)
    real(RP), intent(out) :: ust_org(IA_org,JA_org)
    real(RP), intent(out) :: albg_org(IA_org,JA_org,N_RAD_DIR,N_RAD_RGN)

    real(RP), intent(inout) :: topo_org(IA_org,JA_org)
    real(RP), intent(inout) :: lmask_org(IA_org,JA_org)
    real(RP), intent(inout) :: lz_org(KA_org)

    logical, intent(in) :: use_file_landwater   ! use land water data from files
    integer, intent(in) :: ldims(3)
    integer, intent(in) :: it

    logical :: exist
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( first_lnd ) then
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            topo_org(:,:), vars_land%get("topo"), &
            it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )

       call read1d( KA_org, lz_org(:), vars_land%get("lz"), it, fid_lnd )

       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            lmask_org(:,:), vars_land%get("lsmask"), &
            it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )
    end if

    call read3d( KA_org, KS_org, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         tg_org(:,:,:), vars_land%get("LAND_TEMP"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )


    ! soil liquid water [m3 m-3] (no wrfout-default)
    if( use_file_landwater ) then
       call read3d( KA_org, KS_org, KE_org, IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            strg_org(:,:,:), vars_land%get("LAND_WATER"), &
            it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )
    endif

    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         lst_org(:,:), vars_land%get("LAND_SFC_TEMP"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )

    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         ust_org(:,:), vars_land%get("URBAN_SFC_TEMP"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          ust_org(i,j) = lst_org(i,j)
       end do
       end do
    end if

    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albg_org(:,:,I_R_direct,I_R_VIS), vars_land%get("LAND_SFC_ALB_VIS_dir"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd )
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albg_org(:,:,I_R_diffuse,I_R_VIS), vars_land%get("LAND_SFC_ALB_VIS_dif"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albg_org(i,j,I_R_diffuse,I_R_VIS) = albg_org(i,j,I_R_direct,I_R_VIS)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albg_org(:,:,I_R_direct,I_R_NIR), vars_land%get("LAND_SFC_ALB_NIR_dir"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albg_org(i,j,I_R_direct,I_R_NIR) = albg_org(i,j,I_R_direct,I_R_VIS)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albg_org(:,:,I_R_diffuse,I_R_NIR), vars_land%get("LAND_SFC_ALB_NIR_dif"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albg_org(i,j,I_R_diffuse,I_R_NIR) = albg_org(i,j,I_R_diffuse,I_R_VIS)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albg_org(:,:,I_R_diffuse,I_R_IR), vars_land%get("LAND_SFC_ALB_IR_dif"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
         exist = exist )
    if ( .not. exist ) then
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            albg_org(:,:,I_R_diffuse,I_R_IR), vars_land%get("LAND_SFC_EMIS_IR_dif"), &
            it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
            exist = exist )
       if ( .not. exist ) then
          LOG_ERROR("ParentLandInputNetCDF",*) '"LAND_SFC_ALB_IR_dif" or "LAND_SFC_EMIS_IR_dif" is necessary'
          call PRC_abort
       end if
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albg_org(i,j,I_R_diffuse,I_R_IR) = 1.0_RP - albg_org(i,j,I_R_diffuse,I_R_IR)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albg_org(:,:,I_R_direct,I_R_IR), vars_land%get("LAND_SFC_ALB_IR_dir"), &
         it, nfiles_lnd, fid_lnd, fids_lnd, SCALE_tile_lnd, SCALE_DOMID_lnd, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albg_org(i,j,I_R_direct,I_R_IR) = albg_org(i,j,I_R_diffuse,I_R_IR)
       end do
       end do
    end if

    first_lnd = .false.

    return
  end subroutine ParentLandInputNetCDF

  !-----------------------------------------------------------------------------
  !> Ocean Setup  
  subroutine ParentOceanSetupNetCDF( &
       odims,   &
       timelen, &
       lon_all, &
       lat_all, &
       basename_org, &
       basename_num, &
       serial, &
       do_read )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_file, only: &
       FILE_open, &
       FILE_get_dimLength
    use scale_comm_cartesC, only: &
       COMM_bcast
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_domain_regist_file, &
       COMM_CARTESC_NEST_parent_info

    implicit none

    integer,           intent(out) :: odims(2)
    integer,           intent(out) :: timelen
    real(RP), allocatable, intent(out) :: lon_all(:,:)
    real(RP), allocatable, intent(out) :: lat_all(:,:)

    character(len=*), intent(in)  :: basename_org
    character(len=*), intent(in)  :: basename_num

    logical, intent(inout) :: serial
    logical, intent(inout) :: do_read

    character(len=8) :: FILE_TYPE = "AUTO" ! "SCALE-RM", "WRFARW", "NAMELIST", "AUTO"
    character(len=FILE_HLONG) :: NM_FILE
    logical :: SCALE_MULTI_FILE = .true.
    integer :: SCALE_PARENT_PRC_NUM_X
    integer :: SCALE_PARENT_PRC_NUM_Y
    character(len=FILE_HLONG) :: SCALE_LATLON_CATALOGUE

    namelist / PARAM_MKINIT_REAL_OCEAN_NetCDF / &
         FILE_TYPE, &
         NM_FILE, &
         SCALE_MULTI_FILE, &
         SCALE_PARENT_PRC_NUM_X, &
         SCALE_PARENT_PRC_NUM_Y, &
         SCALE_LATLON_CATALOGUE

    character(len=FILE_HLONG) :: basename
    character(len=FILE_HLONG) :: fname
    integer :: nmfid

    character(len=32) :: items(vars_max)
    integer :: nvars
    type(vinfo), pointer :: var_info
    class(*), pointer :: v

    integer :: n, i
    integer :: ierr
    logical :: error
    !---------------------------------------------------------------------------

    LOG_INFO("ParentOceanSetupNetCDF",*) 'Real Case/Ocean Setup'

    FILE_TYPE = "AUTO"
    NM_FILE = ""
    SCALE_PARENT_PRC_NUM_X = -1
    SCALE_PARENT_PRC_NUM_Y = -1
    SCALE_LATLON_CATALOGUE = ""

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_REAL_OCEAN_NetCDF,iostat=ierr)
    if( ierr > 0 ) then
       LOG_ERROR("ParentOceanSetupNetCDF",*) 'Not appropriate names in namelist PARAM_MKINIT_REAL_OCEAN_NetCDF. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_REAL_OCEAN_NetCDF)


    basename = trim(basename_org) // trim(basename_num)

    fid_ocn = -1
    if ( do_read ) then
       call check_filetype(fid_ocn, FILE_TYPE, basename, SCALE_tile_ocn, "ParentOceanSetupNetCDF")
    end if

    call COMM_bcast( FILE_TYPE )

    if ( FILE_TYPE == "SCALE-RM" ) then
       call COMM_bcast( SCALE_tile_ocn )
       if ( SCALE_tile_ocn ) then
          do_read = .true.
          serial = .false.
       end if
    end if


    if ( do_read ) then
       vars_ocean = hash_table()

       select case( FILE_TYPE )
       case ( "SCALE-RM" )
          xname = "x"
          yname = "y"
          tname = "time"

          call vars_ocean%put("lon", vinfo("lon"))
          call vars_ocean%put("lat", vinfo("lat"))

          call vars_ocean%put("lsmask", vinfo("lsmask"))

          call vars_ocean%put("OCEAN_TEMP", vinfo("OCEAN_TEMP"))

          call vars_ocean%put("OCEAN_SFC_TEMP", vinfo("OCEAN_SFC_TEMP"))
          call vars_ocean%put("OCEAN_SFC_Z0M",  vinfo("OCEAN_SFC_Z0M"))

          call vars_ocean%put("OCEAN_SFC_ALB_IR_dir",   vinfo("OCEAN_SFC_ALB_IR_dir"))
          call vars_ocean%put("OCEAN_SFC_ALB_IR_dif",   vinfo("OCEAN_SFC_ALB_IR_dif"))
          call vars_ocean%put("OCEAN_SFC_ALB_NIR_dir",  vinfo("OCEAN_SFC_ALB_NIR_dir"))
          call vars_ocean%put("OCEAN_SFC_ALB_NIR_dif",  vinfo("OCEAN_SFC_ALB_NIR_dif"))
          call vars_ocean%put("OCEAN_SFC_ALB_VIS_dir",  vinfo("OCEAN_SFC_ALB_VIS_dir"))
          call vars_ocean%put("OCEAN_SFC_ALB_VIS_dif",  vinfo("OCEAN_SFC_ALB_VIS_dif"))

       case ( "WRFARW" )
          xname = "west_east"
          yname = "south_north"
          tname = "Time"

          call vars_ocean%put("lon", vinfo("XLONG"))
          call vars_ocean%put("lat", vinfo("XLAT"))
          call vars_ocean%put("lz",  vinfo("ZS"))

          call vars_ocean%put("topo",   vinfo("HGT"))
          call vars_ocean%put("lsmask", vinfo("LANDMASK"))

          call vars_ocean%put("OCEAN_TEMP", vinfo("OCEAN_TEMP"))

          call vars_ocean%put("OCEAN_SFC_TEMP", vinfo("SST"))
          call vars_ocean%put("OCEAN_SFC_Z0M",  vinfo("ZNT"))

          call vars_ocean%put("OCEAN_SFC_ALB_VIS_dir", vinfo("ALBEDO"))
          call vars_ocean%put("OCEAN_SFC_EMIS_IR_dif",  vinfo("EMISS"))

       case ( "NAMELIST" )
       case default
          LOG_ERROR("ParentOCEANSetupNetCDF",*) 'FILE_TYPE must be "SCALE-RM", "WRFARW", "NAMELIST", or "AUTO", ', trim(FILE_TYPE)
          call PRC_abort
       end select


       !--- read namelist
       nmfid = -1
       if ( NM_FILE /= "" ) then
          nmfid = IO_get_available_fid()
          call IO_get_fname(fname, NM_FILE)
          open(nmfid, file=fname, form="formatted", status="old", action="read", iostat=ierr)
          if ( ierr /= 0 ) then
             LOG_ERROR("ParentOceanSetupNetCDF",*) 'namelist file is not found! ', trim(fname)
             call PRC_abort
          end if

          rewind(nmfid)
          read(nmfid, nml=NetCDF_DIMS, iostat=ierr)
          if( ierr > 0 ) then
             LOG_ERROR("ParentOceanSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_DIMS in ', trim(fname), '. Check!'
             call PRC_abort
          end if

          ! items
          rewind(nmfid)
          nvars = 0
          do n = 1, vars_max
             read(nmfid, nml=NetCDF_ITEM, iostat=ierr)
             if( ierr > 0 ) then
                LOG_ERROR("ParentLandSetupNetCDF",*) 'Not appropriate names in namelist NetCDF_ITEM in ', trim(fname), '. Check!'
                call PRC_abort
             else if( ierr < 0 ) then
                exit
             end if
             nvars = nvars + 1
             items(nvars) = item
          end do
          if ( nvars > vars_max ) then
             LOG_ERROR("ParentLandSetupNetCDF",*) "The number of item in the namelist file exceeds the limit! ", nvars
             call PRC_abort
          end if
          rewind(nmfid)
          do n = 1, nvars
             ! set default
             if ( vars_ocean%has_key(items(n)) ) then
                item = items(n)
                v => vars_ocean%get(item)
                select type(v)
                type is (vinfo)
                   var_info => v
                end select
                name = var_info%name
                fact = var_info%fact
                offset = var_info%offset
             else
                item = items(n)
                name = items(n)
                fact = 1.0_RP
                offset = 0.0_RP
             end if
             zstg = .false.
             xstg = .false.
             ystg = .false.
             read(nmfid, nml=NetCDF_ITEM, iostat=ierr)
             if ( ierr /= 0 ) exit
             ! set params
             call vars_ocean%put(item, vinfo(name=name, zstg=zstg, xstg=xstg, ystg=ystg, fact=fact, offset=offset))
          end do

       else if ( FILE_TYPE == "NAMELIST" ) then
          LOG_ERROR("ParentLANDSetupNetCDF",*) 'NM_FILE is necessary'
          call PRC_abort
       end if

    end if

    if ( SCALE_tile_ocn ) then

       call COMM_CARTESC_NEST_domain_regist_file( &
            SCALE_DOMID_ocn,        & ! [OUT]
            basename,               & ! [IN]
            SCALE_PARENT_PRC_NUM_X, & ! [IN]
            SCALE_PARENT_PRC_NUM_Y, & ! [IN]
            SCALE_LATLON_CATALOGUE  ) ! [IN]

       call COMM_CARTESC_NEST_parent_info( &
            SCALE_DOMID_ocn,    & ! [IN]
            IMAXG=odims(1),     & ! [OUT]
            JMAXG=odims(2),     & ! [OUT]
            num_tile=nfiles_ocn ) ! [OUT]

       allocate( fids_ocn   (nfiles_ocn) )
       allocate( tile_id_ocn(nfiles_ocn) )

       call COMM_CARTESC_NEST_parent_info( &
            SCALE_DOMID_ocn,      & ! [IN]
            tile_id = tile_id_ocn ) ! [OUT]

       call ParentOceanOpenNetCDF( basename_org, basename_num )

    else if ( do_read ) then

       call FILE_get_dimLength( fid_ocn, xname,  odims(1) )
       call FILE_get_dimLength( fid_ocn, yname,  odims(2) )

    end if

    if ( do_read ) then

       call FILE_get_dimLength( fid_ocn, tname, timelen, error=error )
       if ( error ) timelen = 1

       allocate( LON_all(odims(1),odims(2)) )
       allocate( LAT_all(odims(1),odims(2)) )

       call read2d( odims(1), 1, odims(1), odims(2), 1, odims(2), &
            lon_all(:,:), vars_ocean%get("lon"), &
            1, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn )
       lon_all(:,:) = lon_all(:,:) * D2R
       call read2d( odims(1), 1, odims(1), odims(2), 1, odims(2), &
            lat_all(:,:), vars_ocean%get("lat"), &
            1, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn )
       lat_all(:,:) = lat_all(:,:) * D2R

    end if

    first_ocn = .true.

    return
  end subroutine ParentOceanSetupNetCDF

  !-----------------------------------------------------------------------------
  !> Ocean Open
  subroutine ParentOceanOpenNetCDF( &
       basename_org, basename_num )
    use scale_file, only: &
       FILE_open
    implicit none
    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: basename_num

    character(len=FILE_HLONG)  :: basename

    integer :: n
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ParentOceanOpenNetCDF",*) 'Real Case/Ocean Open'

    basename = trim(basename_org) // trim(basename_num)

    if ( SCALE_tile_ocn ) then

       do n = 1, nfiles_ocn
          call FILE_open( &
               basename,             & ! [IN]
               fids_ocn(n),          & ! [OUT]
               aggregate=.false.,    & ! [IN]
               rankid=tile_id_ocn(n) ) ! [IN]
       end do

       fid_ocn = fids_ocn(1)

    else

       call FILE_open(basename, fid_ocn, postfix="")

    end if

    return
  end subroutine ParentOceanOpenNetCDF

  !-----------------------------------------------------------------------------
  !> Ocean Finalize
  subroutine ParentOceanFinalizeNetCDF
    implicit none
    !---------------------------------------------------------------------------

    LOG_INFO("ParentOceanFinalizeNetCDF",*) 'Real Case/Ocean Finalize'

    if ( allocated(fids_ocn)    ) deallocate( fids_ocn )
    if ( allocated(tile_id_ocn) ) deallocate( tile_id_ocn )

    call vars_ocean%destroy()

    SCALE_DOMID_ocn = -1
    nfiles_ocn = 0
    fid_ocn = -1

    return
  end subroutine ParentOceanFinalizeNetCDF

  !-----------------------------------------------------------------------------
  subroutine ParentOceanInputNetCDF( &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       tw_org,    &
       sst_org,   &
       albw_org,  &
       z0w_org,   &
       omask_org, &
       odims,     &
       it         )
    use scale_const, only: &
       D2R => CONST_D2R, &
       UNDEF => CONST_UNDEF
    use scale_file, only: &
       FILE_open, &
       FILE_read
    implicit none
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org

    real(RP), intent(out) :: tw_org(IA_org,JA_org)
    real(RP), intent(out) :: sst_org(IA_org,JA_org)
    real(RP), intent(out) :: albw_org(IA_org,JA_org,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(out) :: z0w_org(IA_org,JA_org)
    real(RP), intent(inout) :: omask_org(IA_org,JA_org)

    integer, intent(in)  :: odims(2)
    integer, intent(in)  :: it

    logical :: exist
    integer :: i, j
    !---------------------------------------------------------------------------

    if ( first_ocn ) then
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            omask_org(:,:), vars_ocean%get("lsmask"), &
            it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn )
    end if

    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         tw_org(:,:), vars_ocean%get("OCEAN_SFC_TEMP"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn )
    sst_org(:,:) = tw_org(:,:)

    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         z0w_org(:,:), vars_ocean%get("OCEAN_SFC_Z0M"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          z0w_org(:,:) = UNDEF
       end do
       end do
    end if

    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albw_org(:,:,I_R_direct,I_R_VIS), vars_ocean%get("OCEAN_SFC_ALB_VIS_dir"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn )
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albw_org(:,:,I_R_diffuse,I_R_VIS), vars_ocean%get("OCEAN_SFC_ALB_VIS_dif"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albw_org(i,j,I_R_diffuse,I_R_VIS) = albw_org(i,j,I_R_direct,I_R_VIS)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albw_org(:,:,I_R_direct,I_R_NIR), vars_ocean%get("OCEAN_SFC_ALB_NIR_dir"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albw_org(i,j,I_R_direct,I_R_NIR) = albw_org(i,j,I_R_direct,I_R_VIS)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albw_org(:,:,I_R_diffuse,I_R_NIR), vars_ocean%get("OCEAN_SFC_ALB_NIR_dif"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albw_org(i,j,I_R_diffuse,I_R_NIR) = albw_org(i,j,I_R_diffuse,I_R_VIS)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albw_org(:,:,I_R_diffuse,I_R_IR), vars_ocean%get("OCEAN_SFC_ALB_IR_dif"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
         exist = exist )
    if ( .not. exist ) then
       call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
            albw_org(:,:,I_R_diffuse,I_R_IR), vars_ocean%get("OCEAN_SFC_EMIS_IR_dif"), &
            it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
            exist = exist )
       if ( .not. exist ) then
          LOG_ERROR("ParentOceanInputNetCDF",*) '"OCEAN_SFC_ALB_IR_dif" or "OCEAN_SFC_EMIS_IR_dif" is necessary'
          call PRC_abort
       end if
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albw_org(i,j,I_R_diffuse,I_R_IR) = 1.0_RP - albw_org(i,j,I_R_diffuse,I_R_IR)
       end do
       end do
    end if
    call read2d( IA_org, IS_org, IE_org, JA_org, JS_org, JE_org, &
         albw_org(:,:,I_R_direct,I_R_IR), vars_ocean%get("OCEAN_SFC_ALB_IR_dir"), &
         it, nfiles_ocn, fid_ocn, fids_ocn, SCALE_tile_ocn, SCALE_DOMID_ocn, &
         exist = exist )
    if ( .not. exist ) then
       !$omp parallel do
       do j = 1, JA_org
       do i = 1, IA_org
          albw_org(i,j,I_R_direct,I_R_IR) = albw_org(i,j,I_R_diffuse,I_R_IR)
       end do
       end do
    end if

    first_ocn = .false.

    return
  end subroutine ParentOceanInputNetCDF

  ! private

  subroutine check_filetype(fid, FILE_TYPE, basename_org, SCALE_tile, subname)
    use scale_file, only: &
       FILE_open, &
       FILE_get_attribute
    integer,          intent(out) :: fid
    character(len=*), intent(inout) :: FILE_TYPE
    logical,          intent(out) :: SCALE_tile
    character(len=*), intent(in) :: basename_org
    character(len=*), intent(in) :: subname

    character(len=FILE_HLONG) :: fname
    character(len=32) :: att
    logical :: exist
    integer :: i

    fname = basename_org
    inquire(file=fname, exist=exist) 
    if ( .not. exist ) then
       fname = trim(basename_org)//".nc"
       inquire(file=fname, exist=exist) 
    end if
    if ( .not. exist ) then
       fname = trim(basename_org)//".pe000000.nc"
       inquire(file=fname, exist=exist) 
    end if
    if ( .not. exist ) then
       LOG_ERROR(subname,*) "file is not found: ", trim(basename_org)
       call PRC_abort
    end if
    call FILE_open(fname, fid, postfix="")
    if ( FILE_TYPE == "AUTO" ) then
       call FILE_get_attribute( &
            fid, "global", "source", &
            att, &
            existed = exist )
       if ( exist .and. att(:8)=="SCALE-RM" ) then
          FILE_TYPE = "SCALE-RM"
          LOG_INFO(subname,*) 'FILE-TYPE SCALE-RM was detected'
       else
          call FILE_get_attribute( &
               fid, "global", "TITLE", &
               att, &
               existed = exist )
          if ( exist .and. index(att, "WRF") > 0 ) then
             FILE_TYPE = "WRFARW"
             LOG_INFO(subname,*) 'FILE-TYPE WRF was detected'
          else
             FILE_TYPE = "NAMELIST"
          end if
       end if
    end if

    scale_tile = .false.
    if ( FILE_TYPE == "SCALE-RM" ) then
       call FILE_get_attribute( &
            fid, "global", "scale_cartesC_prc_num_x", &
            i, &
            existed = exist )
       if ( exist .and. i > 1 ) then
          SCALE_tile = .true.
          LOG_INFO(subname,*) 'Multi files was detected'
       else
          call FILE_get_attribute( &
               fid, "global", "scale_cartesC_prc_num_y", &
               i, &
               existed = exist )
          if ( exist .and. i > 1 ) then
             SCALE_tile = .true.
             LOG_INFO(subname,*) 'Multi files was detected'
          end if
       end if
    end if

    return
  end subroutine check_filetype

  subroutine read3d( &
       KA_org, KS_org, KE_org, &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       val, &
       var, &
       it, &
       nfiles, fid, fids, &
       scale_tile, scale_domid, &
       exist )
    use scale_file, only: &
       FILE_get_dataInfo, &
       FILE_get_shape, &
       FILE_read
    use scale_comm_cartesc_nest, only: &
       COMM_CARTESC_NEST_domain_shape
    integer, intent(in) :: KA_org, KS_org, KE_org
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org

    real(RP), intent(out), target :: val(KA_org,IA_org,JA_org)

    class(*), pointer, intent(in) :: var
    integer,           intent(in) :: it
    integer,           intent(in) :: nfiles
    integer,           intent(in) :: fid, fids(nfiles)
    logical,           intent(in) :: scale_tile
    integer,           intent(in) :: scale_domid

    logical, intent(out), optional :: exist

    real(RP), allocatable :: buf3d(:,:,:)
    real(RP), pointer     :: work(:,:,:)
    real(RP), allocatable, target :: work_t(:,:,:)

    integer :: dims(3)
    integer :: tilei, tilej
    integer :: kmax
    integer :: cxs, cxe, cys, cye
    integer :: pxs, pxe, pys, pye
    logical :: has_tdim
    logical :: transpose
    logical :: exist_
    integer :: i0, i1, j0, j1
    integer :: kst, ist, jst
    integer :: it_
    integer :: k, i, j, n

    if ( .not. associated(var) ) then
       if ( present(exist) ) then
          exist = .false.
       else
          LOG_ERROR("read3d",*) 'data is not found '
          call PRC_abort
       end if
       return
    end if

    select type(var)
    type is (vinfo)
       if ( var%name == "" ) then
          if ( present(exist) ) then
             exist = .false.
          else
             LOG_ERROR("read3d",*) 'data is not found '
             call PRC_abort
          end if
          return
       end if

       call FILE_get_dataInfo( fid, var%name, has_tdim=has_tdim, existed=exist_ )
       if ( .not. exist_ ) then
          if ( present(exist) ) then
             exist = .false.
             return
          else
             LOG_ERROR("read3d",*) 'data is not found: ', trim(var%name)
             call PRC_abort
          end if
       end if

       if ( has_tdim ) then
          it_ = it
       else
          it_ = 1
       end if

       kmax = KE_org - KS_org + 1

       if ( var%zstg ) then
          kst = 1
       else
          kst = 0
       end if
       if ( var%xstg ) then
          ist = 1
       else
          ist = 0
       end if
       if ( var%ystg ) then
          jst = 1
       else
          jst = 0
       end if

       call FILE_get_shape( fid, var%name, dims(:) )
       transpose = dims(1) .ne. kmax+kst

       if ( SCALE_tile ) then
          if ( var%xstg .or. var%ystg ) then
             allocate( work_t(KA_org,IA_org+ist,JA_org+jst) )
             work => work_t
          else
             work => val
          end if
          do n = 1, nfiles
             call COMM_CARTESC_NEST_domain_shape( &
                  tilei, tilej, &
                  cxs, cxe, cys, cye, &
                  pxs, pxe, pys, pye, &
                  SCALE_DOMID, n, &
                  xstg = var%xstg, &
                  ystg = var%ystg )
             i0 = max(IS_org - cxs, 0)
             i1 = max(cxe - IE_org - ist, 0)
             j0 = max(JS_org - cys, 0)
             j1 = max(cye - JE_org - jst, 0)
             if ( transpose ) then
                allocate( buf3d(pxs+i0:pxe-i1,pys+j0:pye-j1,KS_org:KE_org+kst) )
                call FILE_read( fids(n), var%name, buf3d(:,:,:), &
                     step=it_, start=(/pxs+i0,pys+j0,1/), count=(/pxe-pxs+1-i1-i0,pye-pys+1-j1-j0,kmax+kst/) )
                if ( var%zstg ) then
                   !$omp parallel do
                   do j = j0, pye-pys-j1
                   do i = i0, pxe-pxs-i1
                   do k = KS_org, KE_org
                      work(k,cxs+i-IS_org+1,cys+j-JS_org+1) = ( buf3d(pxs+i,pys+j,k) + buf3d(pxs+i,pys+j,k+1) ) * 0.5_RP * var%fact + var%offset
                   end do
                   end do
                   end do
                else
                   !$omp parallel do
                   do j = j0, pye-pys-j1
                   do i = i0, pxe-pxs-i1
                   do k = KS_org, KE_org
                      work(k,cxs+i-IS_org+1,cys+j-JS_org+1) = buf3d(pxs+i,pys+j,k) * var%fact + var%offset
                   end do
                   end do
                   end do
                end if
                deallocate( buf3d )
                if ( var%xstg .and. cxs==2 .and. IS_org==1 ) then ! tentative
                   !$omp parallel do
                   do j = j0, pye-pys-j1
                   do k = KS_org, KE_org
                      work(k,1,cys+j-JS_org+1) = work(k,2,cys+j-JS_org+1)
                   end do
                   end do
                end if
                if ( var%ystg .and. cys==2 .and. JS_org==1 ) then ! tentative
                   !$omp parallel do
                   do i = i0, pxe-pxs-i1
                   do k = KS_org, KE_org
                      work(k,cxs+i-IS_org+1,1) = work(k,cxs+i-IS_org+1,2)
                   end do
                   end do
                end if
             else
                allocate( buf3d(KS_org:KE_org+kst,pxs+i0:pxe-i1,pys+j0:pye-j1) )
                call FILE_read( fids(n), var%name, buf3d(:,:,:), &
                     step=it_, start=(/1,pxs+i0,pys+j0/), count=(/kmax+kst,pxe-pxs+1-i1-i0,pye-pys+1-j1-j0/) )
                if ( var%zstg ) then
                   !$omp parallel do
                   do j = j0, pye-pys-j1
                   do i = i0, pxe-pxs-i1
                   do k = KS_org, KE_org
                      work(k,cxs+i-IS_org+1,cys+j-JS_org+1) = ( buf3d(k,pxs+i,pys+j) + buf3d(k+1,pxs+i,pys+j) ) * 0.5_RP * var%fact + var%offset
                   end do
                   end do
                   end do
                else
                   !$omp parallel do
                   do j = j0, pye-pys-j1
                   do i = i0, pxe-pxs-i1
                   do k = KS_org, KE_org
                      work(k,cxs+i-IS_org+1,cys+j-JS_org+1) = buf3d(k,pxs+i,pys+j) * var%fact + var%offset
                   end do
                   end do
                   end do
                end if
                deallocate( buf3d )
             end if
          end do
          if ( var%xstg ) then
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = KS_org, KE_org
                val(k,i,j) = ( work(k,i,j) + work(k,i+1,j) ) * 0.5_RP
             end do
             end do
             end do
          else if ( var%ystg ) then
             !$omp parallel do collapse(2)
             do j = 1, JA_org
             do i = 1, IA_org
             do k = KS_org, KE_org
                val(k,i,j) = ( work(k,i,j) + work(k,i,j+1) ) * 0.5_RP
             end do
             end do
             end do
          end if
          if ( var%xstg .or. var%ystg ) then
             deallocate( work_t )
          end if
          nullify( work )
       else
          if ( transpose ) then
             allocate( buf3d(IS_org:IE_org+ist,JS_org:JE_org+jst,KS_org:KE_org+kst) )
             call FILE_read( &
                  fid, var%name, &
                  buf3d(:,:,:), &
                  step=it_, &
                  start=(/IS_org,JS_org,1/), &
                  count=(/IA_org+ist,JA_org+jst,kmax+kst/))
             if ( var%zstg ) then
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = ( buf3d(i+IS_org-1,j+JS_org-1,k) + buf3d(i+IS_org-1,j+JS_org-1,k+1) ) * 0.5_RP * var%fact + var%offset
                end do
                end do
                end do
             else if ( var%xstg ) then
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = ( buf3d(i+IS_org-1,j+JS_org-1,k) + buf3d(i+IS_org,j+JS_org-1,k) ) * 0.5_RP * var%fact + var%offset
                end do
                end do
                end do
             else if ( var%ystg ) then
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = ( buf3d(i+IS_org-1,j+JS_org-1,k) + buf3d(i+IS_org-1,j+JS_org,k) ) * 0.5_RP * var%fact + var%offset
                end do
                end do
                end do
             else
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = buf3d(i+IS_org-1,j+JS_org-1,k) * var%fact + var%offset
                end do
                end do
                end do
             end if
             deallocate( buf3d )
          else
             allocate( buf3d(KS_org:KE_org+kst,IS_org:IE_org+ist,JS_org:JE_org+jst) )
             call FILE_read( &
                  fid, var%name, &
                  buf3d(:,:,:), &
                  step=it_, &
                  start=(/1,IS_org,JS_org/), &
                  count=(/kmax+kst,IA_org+ist,JA_org+jst/) )
             if ( var%zstg ) then
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = ( buf3d(k,i+IS_org-1,j+JS_org-1) + buf3d(k+1,i+IS_org-1,j+JS_org-1) ) * 0.5_RP * var%fact + var%offset
                end do
                end do
                end do
             else if ( var%xstg ) then
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = ( buf3d(k,i+IS_org-1,j+JS_org-1) + buf3d(k,i+IS_org,j+JS_org-1) ) * 0.5_RP * var%fact + var%offset
                end do
                end do
                end do
             else if ( var%ystg ) then
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = ( buf3d(k,i+IS_org-1,j+JS_org-1) + buf3d(k,i+IS_org-1,j+JS_org) ) * 0.5_RP * var%fact + var%offset
                end do
                end do
                end do
             else
                !$omp parallel do collapse(2)
                do j = 1, JA_org
                do i = 1, IA_org
                do k = KS_org, KE_org
                   val(k,i,j) = buf3d(k,i+IS_org-1,j+JS_org-1) * var%fact + var%offset
                end do
                end do
                end do
             end if
             deallocate( buf3d )
          end if
       end if

       if ( present(exist) ) exist = .true.
    end select

    return
  end subroutine read3d

  subroutine read2d( &
       IA_org, IS_org, IE_org, &
       JA_org, JS_org, JE_org, &
       val, &
       var, &
       it, &
       nfiles, fid, fids, &
       scale_tile, scale_domid, &
       exist )
    use scale_file, only: &
       FILE_get_dataInfo, &
       FILE_get_shape, &
       FILE_read
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_domain_shape
    integer, intent(in) :: IA_org, IS_org, IE_org
    integer, intent(in) :: JA_org, JS_org, JE_org

    real(RP), intent(out), target :: val(IA_org,JA_org)

    class(*), pointer, intent(in) :: var
    integer,           intent(in) :: it
    integer,           intent(in) :: nfiles
    integer,           intent(in) :: fid, fids(nfiles)
    logical,           intent(in) :: scale_tile
    integer,           intent(in) :: scale_domid

    logical, intent(out), optional :: exist

    real(RP), allocatable :: buf2d(:,:)
    real(RP), pointer     :: work(:,:)
    real(RP), allocatable, target :: work_t(:,:)

    integer :: tilei, tilej
    integer :: cxs, cxe, cys, cye
    integer :: pxs, pxe, pys, pye

    integer :: dims(2)

    integer :: i0, i1, j0, j1
    integer :: ist, jst

    logical :: has_tdim
    logical :: exist_
    integer :: it_
    integer :: n, i, j

    if ( .not. associated(var) ) then
       if ( present(exist) ) then
          exist = .false.
       else
          LOG_ERROR("read3d",*) 'data is not found '
          call PRC_abort
       end if
       return
    end if

    select type(var)
    type is (vinfo)
       if ( var%name == "" ) then
          if ( present(exist) ) then
             exist = .false.
          else
             LOG_ERROR("read2d",*) 'data is not found '
             call PRC_abort
          end if
          return
       end if

       call FILE_get_dataInfo( fid, var%name, has_tdim=has_tdim, existed=exist_ )
       if ( .not. exist_ ) then
          if ( present(exist) ) then
             exist = .false.
             return
          else
             LOG_ERROR("read2d",*) 'data is not found: ', trim(var%name)
             call PRC_abort
          end if
       end if

       if ( has_tdim ) then
          it_ = it
       else
          it_ = 1
       end if

       if ( var%xstg ) then
          ist = 1
       else
          ist = 0
       end if
       if ( var%ystg ) then
          jst = 1
       else
          jst = 0
       end if

       if ( SCALE_DOMID > 0 ) then
          if ( var%xstg .or. var%ystg ) then
             allocate( work_t(IA_org+ist,JA_org+jst) )
             work => work_t
          else
             work => val
          end if
          do n = 1, nfiles
             call COMM_CARTESC_NEST_domain_shape( &
                  tilei, tilej, &
                  cxs, cxe, cys, cye, &
                  pxs, pxe, pys, pye, &
                  SCALE_DOMID, n, &
                  xstg = var%xstg, &
                  ystg = var%ystg )
             i0 = max(IS_org - cxs, 0)
             i1 = max(cxe - IE_org - ist, 0)
             j0 = max(JS_org - cys, 0)
             j1 = max(cye - JE_org - jst, 0)
             allocate( buf2d(pxs+i0:pxe-i1,pys+j0:pye-j1) )
             call FILE_read( fids(n), var%name, buf2d(:,:), &
                  step=it_, start=(/pxs+i0,pys+j0/), count=(/pxe-pxs+1-i1-i0,pye-pys+1-j1-j0/) )
             !$omp parallel do
             do j = j0, pye-pys-j1
             do i = i0, pxe-pxs-i1
                work(cxs+i-IS_org+1,cys+j-JS_org+1) = buf2d(pxs+i,pys+j) * var%fact + var%offset
             end do
             end do
             deallocate( buf2d )
             if ( var%xstg .and. cxs==2 .and. IS_org==1 ) then ! tentative
                !$omp parallel do
                do j = j0, pye-pys-j1
                   work(1,cys+j-JS_org+1) = work(2,cys+j-JS_org+1)
                end do
             end if
             if ( var%ystg .and. cys==2 .and. JS_org==1 ) then ! tentative
                !$omp parallel do
                do i = i0, pxe-pxs-i1
                   work(cxs+i-IS_org+1,1) = work(cxs+i-IS_org+1,2)
                end do
             end if
          end do
          if ( var%xstg ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                val(i,j) = ( work(i,j) + work(i+1,j) ) * 0.5_RP
             end do
             end do
          else if ( var%ystg ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                val(i,j) = ( work(i,j) + work(i,j+1) ) * 0.5_RP
             end do
             end do
          end if
          if ( var%xstg .or. var%ystg ) then
             deallocate( work_t )
          end if
          nullify( work )
       else
          if ( var%xstg .or. var%ystg ) then
             allocate( work_t(IS_org:IE_org+ist,JS_org:JE_org+jst) )
             work => work_t
          else
             work => val
          end if
          call FILE_read( &
               fid, var%name, &
               work(:,:), &
               step=it_, &
               start=(/IS_org,JS_org/), &
               count=(/IA_org+ist,JA_org+jst/) )
          if ( var%xstg ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                val(i,j) = ( work(i,j) + work(i+1,j) ) * 0.5_RP * var%fact + var%offset
             end do
             end do
          else if ( var%ystg ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                val(i,j) = ( work(i,j) + work(i,j+1) ) * 0.5_RP * var%fact + var%offset
             end do
             end do
          else if ( var%fact .ne. 1.0_RP .or. var%offset .ne. 0.0_RP ) then
             !$omp parallel do
             do j = 1, JA_org
             do i = 1, IA_org
                val(i,j) = val(i,j) * var%fact + var%offset
             end do
             end do
          end if
          if ( var%xstg .or. var%ystg ) then
             deallocate( work_t )
          end if
          nullify( work )
       end if

       if ( present(exist) ) exist = .true.
    end select

    return
  end subroutine read2d

  subroutine read1d( &
       KA_org, &
       val, &
       var, &
       it, &
       fid, &
       exist )
    use scale_file, only: &
       FILE_get_dataInfo, &
       FILE_read
    integer, intent(in) :: KA_org

    real(RP), intent(out) :: val(KA_org)

    class(*), pointer, intent(in) :: var
    integer,           intent(in) :: it
    integer,           intent(in) :: fid

    logical, intent(out), optional :: exist

    logical :: has_tdim
    logical :: exist_

    integer :: it_

    if ( .not. associated(var) ) then
       if ( present(exist) ) then
          exist = .false.
       else
          LOG_ERROR("read3d",*) 'data is not found '
          call PRC_abort
       end if
       return
    end if

    select type(var)
    type is ( vinfo )
       if ( var%name == "" ) then
          if ( present(exist) ) then
             exist = .false.
          else
             LOG_ERROR("read1d",*) 'data is not found '
             call PRC_abort
          end if
          return
       end if

       call FILE_get_dataInfo( fid, var%name, has_tdim=has_tdim, existed=exist_ )
       if ( .not. exist_ ) then
          if ( present(exist) ) then
             exist = .false.
             return
          else
             LOG_ERROR("read1d",*) 'data is not found: ', trim(var%name)
             call PRC_abort
          end if
       end if

       if ( has_tdim ) then
          it_ = it
       else
          it_ = 1
       end if

       call FILE_read( fid, var%name, val(:), step=it_ )
       if ( var%fact .ne. 1.0_RP .or. var%offset .ne. 0.0_RP ) then
          val(:) = val(:) * var%fact + var%offset
       end if

       if ( present(exist) ) exist = .true.
    end select

    return
  end subroutine read1d

end module mod_realinput_netcdf
