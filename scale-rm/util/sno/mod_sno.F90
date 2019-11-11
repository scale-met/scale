!-------------------------------------------------------------------------------
!> Module SNO (RM)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_sno
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNO_proc_alloc
  public :: SNO_file_getinfo
  public :: SNO_calc_domainsize
  public :: SNO_read_bcast_1d
  public :: SNO_read_bcast_2d
  public :: SNO_read_map_1d
  public :: SNO_read_map_2d
  public :: SNO_read_map_3d
  public :: SNO_attributes_write

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNO_proc_alloc( &
       nprocs,       &
       myrank,       &
       ismaster,     &
       nprocs_x_out, &
       nprocs_y_out, &
       pstr,         &
       pend          )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)  :: nprocs       ! number of processes               (execution)
    integer, intent(in)  :: myrank       ! my rank                           (execution)
    logical, intent(in)  :: ismaster     ! master process?                   (execution)
    integer, intent(in)  :: nprocs_x_out ! x length of 2D processor topology (output)
    integer, intent(in)  :: nprocs_y_out ! y length of 2D processor topology (output)
    integer, intent(out) :: pstr         ! start index of peXXXXXX to manage (execution)
    integer, intent(out) :: pend         ! end   index of peXXXXXX to manage (execution)

    integer :: nprocs_out        ! number of peXXXXXX files           (output)
    integer :: pe_per_rank
    !---------------------------------------------------------------------------

    nprocs_out  = nprocs_x_out * nprocs_y_out
    pe_per_rank = (nprocs_out-1) / nprocs + 1
    pstr        = myrank * pe_per_rank
    pend        = myrank * pe_per_rank + pe_per_rank - 1
    pend        = min( pend, nprocs_out-1 )

    if ( nprocs > nprocs_out ) then
       LOG_ERROR("SNO_proc_alloc",*) '# of using PEs should be less than # of files for output. Check'
       call PRC_abort
    endif

    LOG_NEWLINE
    LOG_INFO("SNO_proc_alloc",*) 'Process management'
    LOG_INFO_CONT(*)             '# of PEs (this execution)     : ', nprocs
    LOG_INFO_CONT(*)             'Rank     (this execution)     : ', myrank
    LOG_INFO_CONT(*)             'Master?  (this execution)     : ', ismaster
    LOG_INFO_CONT(*)             '# of PEs (output file)        : ', nprocs_out, '(', nprocs_x_out, 'x', nprocs_y_out, ')'
    LOG_INFO_CONT(*)             'Start PE managed by this rank : ', pstr
    LOG_INFO_CONT(*)             'End   PE managed by this rank : ', pend

    return
  end subroutine SNO_proc_alloc

  !-----------------------------------------------------------------------------
  subroutine SNO_file_getinfo( &
       ismaster,     &
       basename,     &
       vars,         &
       nprocs_x_out, &
       nprocs_y_out, &
       nprocs_x_in,  &
       nprocs_y_in,  &
       ngrids_z,     &
       ngrids_x,     &
       ngrids_y,     &
       nhalos_z,     &
       nhalos_x,     &
       nhalos_y,     &
       hinfo,        &
       naxis,        &
       axisname,     &
       nvars,        &
       varname,      &
       plugin_hgridope, &
       debug         )
    use mpi
    use scale_file_h, only: &
       FILE_FREAD
    use scale_file, only: &
       FILE_open,           &
       FILE_get_commonInfo, &
       FILE_get_attribute
    use scale_file_cartesc, only: &
       FILE_CARTESC_get_size
    use scale_atmos_grid_cartesc, only: &
       ATMOS_GRID_CARTESC_NAME
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_abort
    use scale_const, only: &
       CONST_UNDEF
    use mod_sno_h, only: &
       item_limit, &
       commoninfo
    implicit none

    logical,                intent(in)  :: ismaster                 ! master process?                   (execution)
    character(len=*),       intent(in)  :: basename                 ! basename of file                  (input)
    character(len=*),       intent(in)  :: vars(item_limit)         ! name of variables, requested      (execution)
    integer,                intent(in)  :: nprocs_x_out             ! x length of 2D processor topology (output)
    integer,                intent(in)  :: nprocs_y_out             ! y length of 2D processor topology (output)
    integer,                intent(out) :: nprocs_x_in              ! x length of 2D processor topology (input)
    integer,                intent(out) :: nprocs_y_in              ! y length of 2D processor topology (input)
    integer,                intent(out) :: ngrids_z                 ! size of z-axis grids              (global,sometimes including halo)
    integer,                intent(out) :: ngrids_x                 ! size of x-axis grids              (global,sometimes including halo)
    integer,                intent(out) :: ngrids_y                 ! size of y-axis grids              (global,sometimes including halo)
    integer,                intent(out) :: nhalos_z                 ! size of z-axis halo grids         (global,sometimes have a size)
    integer,                intent(out) :: nhalos_x                 ! size of x-axis halo grids         (global,sometimes have a size)
    integer,                intent(out) :: nhalos_y                 ! size of y-axis halo grids         (global,sometimes have a size)
    type(commoninfo),       intent(out) :: hinfo                    ! common information                (input)
    integer,                intent(out) :: naxis                    ! number of axis variables          (input)
    character(len=H_SHORT), intent(out) :: axisname(item_limit)     ! name   of axis variables          (input)
    integer,                intent(out) :: nvars                    ! number of variables               (input)
    character(len=H_SHORT), intent(out) :: varname (item_limit)     ! name   of variables               (input)
    logical,                intent(in)  :: plugin_hgridope          ! if true, some 2d axis var. is treated as normal var.
    logical,                intent(in)  :: debug

    integer                :: procsize(2)                   ! total process size        (x:y)

    integer                :: nvars_file                    ! number of variables from input file
    character(len=H_SHORT) :: varname_file(item_limit) = '' ! name   of variables from input file
    integer                :: nvars_req                     ! number of variables, requested
    logical                :: exist       (item_limit)      ! check flag for requested variables

    integer :: nprocs_in       ! number of peXXXXXX files           (input)

    integer :: ngrids_z_nohalo ! number of z-axis grids             (global,without halo)
    integer :: ngrids_x_nohalo ! number of x-axis grids             (global,without halo)
    integer :: ngrids_y_nohalo ! number of y-axis grids             (global,without halo)
    integer :: ngrids          ! number of        grids             (global)
    integer :: ngrids_x_in     ! number of x-axis grids per process (input,without halo)
    integer :: ngrids_y_in     ! number of y-axis grids per process (input,without halo)
    integer :: ngrids_in       ! number of        grids per process (input,without halo)
    integer :: ngrids_x_out    ! number of x-axis grids per process (output,without halo)
    integer :: ngrids_y_out    ! number of y-axis grids per process (output,without halo)
    integer :: ngrids_out      ! number of        grids per process (output,without halo)

    integer :: nowrank
    integer :: n, nn
    integer :: fid, ierr
    logical :: existed
    !---------------------------------------------------------------------------

    !---  read info from global file

    LOG_NEWLINE
    LOG_INFO("SNO_file_getinfo",*) 'Get process & horizontal grid info from input file'

    if ( ismaster ) then
       nowrank = 0 ! first file

       call FILE_open( basename,        & ! [IN]
                       fid,             & ! [OUT]
                       rankid = nowrank ) ! [IN]

       call FILE_get_commonInfo( fid,             & ! [IN]
                                 item_limit,      & ! [IN]
                                 hinfo%title,     & ! [OUT]
                                 hinfo%source,    & ! [OUT]
                                 hinfo%institute, & ! [OUT]
                                 nvars_file,      & ! [OUT]
                                 varname_file(:)  ) ! [OUT]


       call FILE_get_attribute( fid, "global", "scale_cartesC_prc_num_x", procsize(1:1) )
       call FILE_get_attribute( fid, "global", "scale_cartesC_prc_num_y", procsize(2:2) )

       call FILE_get_attribute( fid, "global", "scale_cartesC_prc_periodic_z", hinfo%periodic(1) )
       call FILE_get_attribute( fid, "global", "scale_cartesC_prc_periodic_x", hinfo%periodic(2) )
       call FILE_get_attribute( fid, "global", "scale_cartesC_prc_periodic_y", hinfo%periodic(3) )

       call FILE_get_attribute( fid, "global", "grid_name", hinfo%grid_name )

       select case ( hinfo%grid_name )
       case ( ATMOS_GRID_CARTESC_NAME )

          call FILE_CARTESC_get_size( fid,                                                     & ! (in)
                                      hinfo%gridsize(1),                                       & ! (out) KMAX
                                      hinfo%gridsize(4), hinfo%gridsize(5), hinfo%gridsize(6), & ! (out) OKMAX, LKMAX, UKMAX
                                      hinfo%gridsize(2), hinfo%gridsize(3),                    & ! (out) IMAXG, JMAXG
                                      hinfo%halosize(1), hinfo%halosize(2), hinfo%halosize(3)  ) ! (out) KHALO, IHALO, JHALO

          call FILE_get_attribute( fid, "global", "time_units", hinfo%time_units    )
          call FILE_get_attribute( fid, "global", "calendar", hinfo%calendar, existed )
          if ( .NOT. existed ) hinfo%calendar = ""
          call FILE_get_attribute( fid, "global", "time_start", hinfo%time_start(:) )

          call FILE_get_attribute( fid, 'x', 'size_global', hinfo%xatt_size_global(:) )
          call FILE_get_attribute( fid, 'x', 'halo_global', hinfo%xatt_halo_global(:) )

          call FILE_get_attribute( fid, 'y', 'size_global', hinfo%yatt_size_global(:) )
          call FILE_get_attribute( fid, 'y', 'halo_global', hinfo%yatt_halo_global(:) )
       case default
          LOG_ERROR("SNO_file_getinfo",*) 'currently only ', trim(ATMOS_GRID_CARTESC_NAME), ' is supported as grid_name'
          call PRC_abort
       end select

    endif

    call MPI_BCAST( hinfo%title              , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%source             , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%institute          , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%grid_name          , H_SHORT, MPI_CHARACTER     , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%periodic        (1), 3    , MPI_LOGICAL         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%gridsize        (1), 6    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%halosize        (1), 3    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%time_units         , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%calendar           , H_SHORT, MPI_CHARACTER     , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%time_start      (1), 1    , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%xatt_size_global(1), 1    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%xatt_halo_global(1), 2    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%yatt_size_global(1), 1    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%yatt_halo_global(1), 2    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( procsize              (1), 2    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    call MPI_BCAST( nvars_file     , 1                 , MPI_INTEGER,   PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( varname_file(1), H_SHORT*item_limit, MPI_CHARACTER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    if ( debug ) then
       LOG_INFO("SNO_file_getinfo",*) "hinfo%title           ", trim(hinfo%title)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%source          ", trim(hinfo%source)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%institute       ", trim(hinfo%institute)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%grid_name       ", trim(hinfo%grid_name)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%periodic        ", hinfo%periodic        (:)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%gridsize        ", hinfo%gridsize        (:)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%halosize        ", hinfo%halosize        (:)
       LOG_INFO("SNO_file_getinfo",*) "procsize              ", procsize              (:)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%xatt_size_global", hinfo%xatt_size_global(:)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%xatt_halo_global", hinfo%xatt_halo_global(:)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%yatt_size_global", hinfo%yatt_size_global(:)
       LOG_INFO("SNO_file_getinfo",*) "hinfo%yatt_halo_global", hinfo%yatt_halo_global(:)
    endif

    !--- process management

    nprocs_x_in = procsize(1)
    nprocs_y_in = procsize(2)
    nprocs_in   = nprocs_x_in * nprocs_y_in

    LOG_NEWLINE
    LOG_INFO("SNO_file_getinfo",*) 'Process info '
    LOG_INFO_CONT(*)               '# of PEs (input  file)        : ', nprocs_in, '(', nprocs_x_in, 'x', nprocs_y_in, ')'

    !--- grid management

    ngrids_z        = hinfo%gridsize(1)
    ngrids_x        = hinfo%xatt_size_global(1)
    ngrids_y        = hinfo%yatt_size_global(1)
    nhalos_z        = hinfo%halosize(1)         ! assume both side have same number
    nhalos_x        = hinfo%xatt_halo_global(1) ! assume both side have same number
    nhalos_y        = hinfo%yatt_halo_global(1) ! assume both side have same number

    ngrids_z_nohalo = ngrids_z - 2*nhalos_z
    ngrids_x_nohalo = ngrids_x - 2*nhalos_x
    ngrids_y_nohalo = ngrids_y - 2*nhalos_y

    ngrids_x_in     = ngrids_x_nohalo / nprocs_x_in
    ngrids_y_in     = ngrids_y_nohalo / nprocs_y_in

    ngrids_x_out    = ngrids_x_nohalo / nprocs_x_out
    ngrids_y_out    = ngrids_y_nohalo / nprocs_y_out

    ngrids          = ngrids_x_nohalo * ngrids_y_nohalo
    ngrids_in       = ngrids_x_in     * ngrids_y_in
    ngrids_out      = ngrids_x_out    * ngrids_y_out

    LOG_NEWLINE
    LOG_INFO("SNO_file_getinfo",*) 'Grid info (horizontal) '
    LOG_INFO_CONT(*)               '# of total grids           without halo : ', ngrids    , '(', ngrids_x_nohalo, 'x', ngrids_y_nohalo, ')'
    LOG_INFO_CONT(*)               '# of grids per input  file without halo : ', ngrids_in , '(', ngrids_x_in    , 'x', ngrids_y_in    , ')'
    LOG_INFO_CONT(*)               '# of grids per output file without halo : ', ngrids_out, '(', ngrids_x_out   , 'x', ngrids_y_out   , ')'
    LOG_INFO_CONT(*)               '# of halo grids (x-axis,one side)       : ', nhalos_x
    LOG_INFO_CONT(*)               '# of halo grids (y-axis,one side)       : ', nhalos_y
    LOG_NEWLINE
    LOG_INFO("SNO_file_getinfo",*) 'Grid info (vertical) '
    LOG_INFO_CONT(*)               '# of vertical grids        without halo : ', ngrids_z_nohalo
    LOG_INFO_CONT(*)               '# of halo grids (z-axis,one side)       : ', nhalos_z

    if ( mod(ngrids_x_nohalo,nprocs_x_out) /= 0 ) then
       LOG_ERROR("SNO_file_getinfo",*) 'The allowable case is that # of the total x-grids is divisible with # of the x-PEs for the output. Stop'
       call PRC_abort
    endif

    if ( mod(ngrids_y_nohalo,nprocs_y_out) /= 0 ) then
       LOG_ERROR("SNO_file_getinfo",*) 'The allowable case is that # of the total y-grids is divisible with # of the y-PEs for the output. Stop'
       call PRC_abort
    endif

    exist(:) = .false.
    do nn = 1, item_limit
       if( vars(nn) == '' ) exit
    enddo
    nvars_req = nn - 1

    if ( debug ) then
       LOG_NEWLINE
       LOG_INFO("SNO_file_getinfo",*) "nvars_req = ", nvars_req
       do nn = 1, nvars_req
          LOG_INFO("SNO_file_getinfo",*) nn, trim(vars(nn))
       enddo
       LOG_NEWLINE
       LOG_INFO("SNO_file_getinfo",*) "nvars_file = ", nvars_file
    endif

    hinfo%minfo_mapping_name                             = ''
    hinfo%minfo_false_easting                        (:) = CONST_UNDEF
    hinfo%minfo_false_northing                       (:) = CONST_UNDEF
    hinfo%minfo_longitude_of_central_meridian        (:) = CONST_UNDEF
    hinfo%minfo_longitude_of_projection_origin       (:) = CONST_UNDEF
    hinfo%minfo_latitude_of_projection_origin        (:) = CONST_UNDEF
    hinfo%minfo_straight_vertical_longitude_from_pole(:) = CONST_UNDEF
    hinfo%minfo_standard_parallel                    (:) = CONST_UNDEF

    nvars = 0
    naxis = 0
    do n = 1, nvars_file

       if ( debug ) then
          LOG_INFO("SNO_file_getinfo",*) n, trim(varname_file(n))
       endif

       select case(varname_file(n))
       case('z','zh','oz','ozh','lz','lzh','uz','uzh','pressure',                                                            &
            'z_bnds','zh_bnds','oz_bnds','ozh_bnds','lz_bnds','lzh_bnds','uz_bnds','uzh_bnds','pressure_bnds',               &
            'CZ','FZ','CDZ','FDZ','CBFZ','FBFZ','OCZ','OFZ','OCDZ','LCZ','LFZ','LCDZ','UCZ','UFZ','UCDZ',                    &
            'x','xh','y','yh',                                                                                               &
            'x_bnds','xh_bnds','y_bnds','yh_bnds',                                                                           &
            'CX','CY','FX','FY','CDX','CDY','FDX','FDY','CBFX','CBFY','FBFX','FBFY',                                         &
            'CXG','CYG','FXG','FYG','CDXG','CDYG','FDXG','FDYG','CBFXG','CBFYG','FBFXG','FBFYG',                             &
            'lon','lon_uy','lon_xv','lon_uv','lat','lat_uy','lat_xv','lat_uv',                               &
            'cell_area','cell_area_uy','cell_area_xv',                                                                       &
            'cell_area_zuy_x','cell_area_zxv_y','cell_area_wuy_x','cell_area_wxv_y',                                         &
            'cell_area_zxy_x','cell_area_zuv_y','cell_area_zuv_x','cell_area_zxy_y',                                         &
            'cell_area_uyz_x','cell_area_xvz_y','cell_area_uyw_x','cell_area_xvw_y',                                         &
            'cell_area_xyz_x','cell_area_uvz_y','cell_area_uvz_x','cell_area_xyz_y',                                         &
            'cell_volume',                                                                                                   &
            'cell_volume_wxy','cell_volume_zuy','cell_volume_zxv','cell_volume_oxy','cell_volume_lxy','cell_volume_uxy',     &
            'cell_volume_xyw','cell_volume_uyz','cell_volume_xvz','cell_volume_xyo','cell_volume_xyl','cell_volume_xyu',     &
            'height','height_wxy','height_xyw','height_uyz','height_xvz','height_uvz','height_uyw','height_xvw','height_uvw' )
          naxis           = naxis + 1
          axisname(naxis) = varname_file(n)
       case('grid','grid_ocean','grid_land','grid_urban','grid_pressure','grid_z','grid_model','grid_model_global')
          ! do nothing
       case('lambert_conformal_conic')
          if ( ismaster ) then
             call FILE_get_attribute( fid, varname_file(n), "grid_mapping_name",             &
                                      hinfo%minfo_mapping_name                               )
             call FILE_get_attribute( fid, varname_file(n), "false_easting",                 &
                                      hinfo%minfo_false_easting                        (:)   )
             call FILE_get_attribute( fid, varname_file(n), "false_northing",                &
                                      hinfo%minfo_false_northing                       (:)   )
             call FILE_get_attribute( fid, varname_file(n), "longitude_of_central_meridian", &
                                      hinfo%minfo_longitude_of_central_meridian        (:)   )
             call FILE_get_attribute( fid, varname_file(n), "latitude_of_projection_origin", &
                                      hinfo%minfo_latitude_of_projection_origin        (:)   )
             call FILE_get_attribute( fid, varname_file(n), "standard_parallel",             &
                                      hinfo%minfo_standard_parallel                    (:)   )
          endif
       case('polar_stereographic')
          if ( ismaster ) then
             call FILE_get_attribute( fid, varname_file(n), "grid_mapping_name",                     &
                                      hinfo%minfo_mapping_name                             )
             call FILE_get_attribute( fid, varname_file(n), "false_easting",                         &
                                      hinfo%minfo_false_easting                        (:) )
             call FILE_get_attribute( fid, varname_file(n), "false_northing",                        &
                                       hinfo%minfo_false_northing                       (:) )
             call FILE_get_attribute( fid, varname_file(n), "latitude_of_projection_origin",         &
                                      hinfo%minfo_latitude_of_projection_origin        (:) )
             call FILE_get_attribute( fid, varname_file(n), "straight_vertical_longitude_from_pole", &
                                      hinfo%minfo_straight_vertical_longitude_from_pole(:) )
             call FILE_get_attribute( fid, varname_file(n), "standard_parallel",                     &
                                      hinfo%minfo_standard_parallel                  (1:1) )
          endif
       case('mercator')
          if ( ismaster ) then
             call FILE_get_attribute( fid, varname_file(n), "grid_mapping_name",              &
                                      hinfo%minfo_mapping_name                             )
             call FILE_get_attribute( fid, varname_file(n), "false_easting",                  &
                                      hinfo%minfo_false_easting                        (:) )
             call FILE_get_attribute( fid, varname_file(n), "false_northing",                 &
                                      hinfo%minfo_false_northing                       (:) )
             call FILE_get_attribute( fid, varname_file(n), "longitude_of_projection_origin", &
                                      hinfo%minfo_longitude_of_projection_origin       (:) )
          endif
       case('equirectangular')
          if ( ismaster ) then
             call FILE_get_attribute( fid, varname_file(n), "grid_mapping_name",             &
                                      hinfo%minfo_mapping_name                             )
             call FILE_get_attribute( fid, varname_file(n), "false_easting",                 &
                                      hinfo%minfo_false_easting                        (:) )
             call FILE_get_attribute( fid, varname_file(n), "false_northing",                &
                                      hinfo%minfo_false_northing                       (:) )
             call FILE_get_attribute( fid, varname_file(n), "longitude_of_central_meridian", &
                                      hinfo%minfo_longitude_of_central_meridian        (:) )
          endif
       case('topo','lsmask')
          if ( plugin_hgridope ) then ! treat as normal variable
             if ( nvars_req == 0 ) then
                nvars          = nvars + 1
                varname(nvars) = varname_file(n)
             else
                do nn = 1, nvars_req
                   if ( varname_file(n) == vars(nn) ) then
                      if ( exist(nn) ) then
                         LOG_ERROR("SNO_file_getinfo",*) 'variable ', trim(vars(nn)), ' is requested two times. check namelist!'
                         call PRC_abort
                      endif

                      nvars          = nvars + 1
                      varname(nvars) = varname_file(n)
                      exist(nn)      = .true.
                   endif
                enddo
             endif
          else ! treat as axis variable
             naxis           = naxis + 1
             axisname(naxis) = varname_file(n)
          endif
       case default
          if( index( varname_file(n), 'time' ) > 0 ) then
             ! do nothing
          else
             ! do nothing
             if ( nvars_req == 0 ) then
                nvars           = nvars + 1
                varname (nvars) = varname_file(n)
             else
                do nn = 1, nvars_req
                   if ( varname_file(n) == vars(nn) ) then
                      if ( exist(nn) ) then
                         LOG_ERROR("SNO_file_getinfo",*) 'variable ', trim(vars(nn)), &
                                                        ' is requested two times. check namelist!'
                         call PRC_abort
                      endif

                      nvars          = nvars + 1
                      varname(nvars) = varname_file(n)
                      exist(nn)      = .true.
                   endif
                enddo
             endif
          endif
       end select

    enddo

    if ( nvars_req > 0 .AND. nvars /= nvars_req ) then
       LOG_NEWLINE
       LOG_ERROR("SNO_file_getinfo",*) 'some requested variables are missing.', nvars, nvars_req
       LOG_ERROR_CONT(*) 'nvars = ', nvars, ', nvars_req = ', nvars_req
       do nn = 1, nvars_req
          LOG_ERROR_CONT(*) 'check:', exist(nn), ', name: ', trim(vars(nn))
       enddo
       call PRC_abort
    endif

    call MPI_BCAST( hinfo%minfo_mapping_name                            , H_SHORT, MPI_CHARACTER,        PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_false_easting                        (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_false_northing                       (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_longitude_of_central_meridian        (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_longitude_of_projection_origin       (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_latitude_of_projection_origin        (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_straight_vertical_longitude_from_pole(1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_standard_parallel                    (1), 2      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    LOG_NEWLINE
    LOG_INFO("SNO_file_getinfo",*) 'Axis list'
    do n = 1, naxis
       LOG_INFO_CONT('(1x,A,I3.3,A,A20)') 'No.', n, ', name: ', trim(axisname(n))
    enddo

    LOG_NEWLINE
    LOG_INFO("SNO_file_getinfo",*) 'Variable list'
    do n = 1, nvars
       LOG_INFO_CONT('(1x,A,I3.3,A,A20)') 'No.', n, ', name: ', trim(varname(n))
    enddo

    return
  end subroutine SNO_file_getinfo

  !-----------------------------------------------------------------------------
  subroutine SNO_calc_domainsize( &
       nprocs_x,       &
       nprocs_y,       &
       px,             &
       py,             &
       ngrids_x_total, &
       ngrids_y_total, &
       nhalos_x,       &
       nhalos_y,       &
       hinfo,          &
       ngrids_x,       &
       ngrids_y,       &
       ngrids_xh,      &
       ngrids_yh       )
    use mod_sno_h, only: &
       commoninfo
    implicit none

    integer,          intent(in)  :: nprocs_x       ! x length of 2D processor topology
    integer,          intent(in)  :: nprocs_y       ! y length of 2D processor topology
    integer,          intent(in)  :: ngrids_x_total ! total grid size of x-axis         (sometimes including halo)
    integer,          intent(in)  :: ngrids_y_total ! total grid size of y-axis         (sometimes including halo)
    integer,          intent(in)  :: nhalos_x       ! halo  grid size of x-axis         (sometimes have a size)
    integer,          intent(in)  :: nhalos_y       ! halo  grid size of y-axis         (sometimes have a size)
    integer,          intent(in)  :: px             ! x index  in 2D processor topology
    integer,          intent(in)  :: py             ! y index  in 2D processor topology
    type(commoninfo), intent(in)  :: hinfo          ! common information                (input)
    integer,          intent(out) :: ngrids_x       ! size of x-axis grids              (sometimes including halo)
    integer,          intent(out) :: ngrids_y       ! size of y-axis grids              (sometimes including halo)
    integer,          intent(out) :: ngrids_xh      ! size of x-axis grids              (sometimes including halo)
    integer,          intent(out) :: ngrids_yh      ! size of y-axis grids              (sometimes including halo)
    !---------------------------------------------------------------------------

    ngrids_y = ( ngrids_y_total - 2*nhalos_y ) / nprocs_y
    if( py == 1        ) ngrids_y = ngrids_y + nhalos_y
    if( py == nprocs_y ) ngrids_y = ngrids_y + nhalos_y

    ngrids_yh = ngrids_y
    if ( nhalos_y == 0 .AND. py == 1 .AND. (.NOT. hinfo%periodic(3)) ) then
       ngrids_yh = ngrids_yh + 1
    endif

    ngrids_x = ( ngrids_x_total - 2*nhalos_x ) / nprocs_x
    if( px == 1        ) ngrids_x = ngrids_x + nhalos_x
    if( px == nprocs_x ) ngrids_x = ngrids_x + nhalos_x

    ngrids_xh = ngrids_x
    if ( nhalos_x == 0 .AND. px == 1 .AND. (.NOT. hinfo%periodic(2)) ) then
       ngrids_xh = ngrids_xh + 1
    endif

    return
  end subroutine SNO_calc_domainsize

  !-----------------------------------------------------------------------------
  subroutine SNO_read_bcast_1d( &
       ismaster,  &
       basename,  &
       varname,   &
       datatype,  &
       varsize,   &
       var        )
    use mpi
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4, &
       UNDEF8 => CONST_UNDEF8, &
       UNDEF  => CONST_UNDEF,  &
       EPS    => CONST_EPS
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8, &
       FILE_dtypelist
    use scale_file, only: &
       FILE_read
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_abort
    implicit none

    logical,          intent(in)  :: ismaster
    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: datatype
    integer,          intent(in)  :: varsize
    real(RP),         intent(out) :: var(varsize)

    real(SP), allocatable :: tmp_SP(:)
    real(DP), allocatable :: tmp_DP(:)

    integer :: nowrank
    integer :: datatype_mpi
    integer :: ierr
    !---------------------------------------------------------------------------

    nowrank = 0 ! first file

    if    ( RP == 4 ) then
       datatype_mpi = MPI_REAL
    elseif( RP == 8 ) then
       datatype_mpi = MPI_DOUBLE_PRECISION
    endif

    if ( datatype == FILE_REAL4 ) then

       if ( ismaster ) then
          allocate( tmp_SP(varsize) )

          call FILE_read( basename, varname, tmp_SP(:), rankid=nowrank )

          where( abs( tmp_SP(:) - UNDEF4 ) > EPS )
             var(:) = real(tmp_SP(:),kind=RP)
          elsewhere
             var(:) = UNDEF
          endwhere

          deallocate( tmp_SP )
       endif

    elseif( datatype == FILE_REAL8 ) then

       if ( ismaster ) then
          allocate( tmp_DP(varsize) )

          call FILE_read( basename, varname, tmp_DP(:), rankid=nowrank )

          where( abs( tmp_DP(:) - UNDEF8 ) > EPS )
             var(:) = real(tmp_DP(:),kind=RP)
          elsewhere
             var(:) = UNDEF
          endwhere

          deallocate( tmp_DP )
       endif

    else
       LOG_ERROR("SNO_read_bcast_1d",*) '[read_bcast_1d] Unsupported datatype : ', trim(FILE_dtypelist(datatype))
       call PRC_abort
    endif

    call MPI_BCAST( var(:), varsize, datatype_mpi, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    return
  end subroutine SNO_read_bcast_1d

  !-----------------------------------------------------------------------------
  subroutine SNO_read_map_1d( &
       basename, &
       nowrank,  &
       nowstep,  &
       varname,  &
       datatype, &
       gin1,     &
       gout1,    &
       localmap, &
       var       )
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4, &
       UNDEF8 => CONST_UNDEF8, &
       UNDEF  => CONST_UNDEF,  &
       EPS    => CONST_EPS
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8, &
       FILE_dtypelist
    use scale_file, only: &
       FILE_read
    use scale_prc, only: &
       PRC_abort
    use mod_sno_h, only: &
       I_map_p, &
       I_map_i
    implicit none

    character(len=*), intent(in)    :: basename
    integer,          intent(in)    :: nowrank
    integer,          intent(in)    :: nowstep
    character(len=*), intent(in)    :: varname
    integer,          intent(in)    :: datatype
    integer,          intent(in)    :: gin1
    integer,          intent(in)    :: gout1
    integer(2),       intent(in)    :: localmap(gout1,2)
    real(RP),         intent(inout) :: var     (gout1)

    real(SP) :: tmp_SP(gin1)
    real(DP) :: tmp_DP(gin1)

    integer  :: i, ii
    !---------------------------------------------------------------------------

    if ( datatype == FILE_REAL4 ) then

       call FILE_read( basename, varname, tmp_SP(:), step=nowstep, rankid=nowrank )

       do i = 1, gout1
          if ( localmap(i,I_map_p) == nowrank ) then
             ii = int(localmap(i,I_map_i))

             if ( abs( tmp_SP(ii) - UNDEF4 ) > EPS ) then
                var(i) = real(tmp_SP(ii),kind=RP)
             else
                var(i) = UNDEF
             endif
          endif
       enddo

    elseif( datatype == FILE_REAL8 ) then

       call FILE_read( basename, varname, tmp_DP(:), step=nowstep, rankid=nowrank )

       do i = 1, gout1
          if ( localmap(i,I_map_p) == nowrank ) then
             ii = int(localmap(i,I_map_i))

             if ( abs( tmp_DP(ii) - UNDEF8 ) > EPS ) then
                var(i) = real(tmp_DP(ii),kind=RP)
             else
                var(i) = UNDEF
             endif
          endif
       enddo

    else
       LOG_ERROR("SNO_read_map_1d",*) '[read_map_1d] Unsupported datatype : ', trim(FILE_dtypelist(datatype))
       call PRC_abort
    endif

    return
  end subroutine SNO_read_map_1d

  !-----------------------------------------------------------------------------
  subroutine SNO_read_bcast_2d( &
       ismaster,  &
       basename,  &
       varname,   &
       datatype,  &
       varsize1,  &
       varsize2,  &
       var        )
    use mpi
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4, &
       UNDEF8 => CONST_UNDEF8, &
       UNDEF  => CONST_UNDEF,  &
       EPS    => CONST_EPS
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8, &
       FILE_dtypelist
    use scale_file, only: &
       FILE_read
    use scale_prc, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_abort
    implicit none

    logical,          intent(in)  :: ismaster
    character(len=*), intent(in)  :: basename
    character(len=*), intent(in)  :: varname
    integer,          intent(in)  :: datatype
    integer,          intent(in)  :: varsize1
    integer,          intent(in)  :: varsize2
    real(RP),         intent(out) :: var(varsize1,varsize2)

    real(SP), allocatable :: tmp_SP(:,:)
    real(DP), allocatable :: tmp_DP(:,:)

    integer :: nowrank
    integer :: datatype_mpi
    integer :: ierr
    !---------------------------------------------------------------------------

    nowrank = 0 ! first file

    if    ( RP == 4 ) then
       datatype_mpi = MPI_REAL
    elseif( RP == 8 ) then
       datatype_mpi = MPI_DOUBLE_PRECISION
    endif

    if ( datatype == FILE_REAL4 ) then

       if ( ismaster ) then
          allocate( tmp_SP(varsize1,varsize2) )

          call FILE_read( basename, varname, tmp_SP(:,:), rankid=nowrank )

          where( abs( tmp_SP(:,:) - UNDEF4 ) > EPS )
             var(:,:) = real(tmp_SP(:,:),kind=RP)
          elsewhere
             var(:,:) = UNDEF
          endwhere

          deallocate( tmp_SP )

       endif

    elseif( datatype == FILE_REAL8 ) then

       if ( ismaster ) then
          allocate( tmp_DP(varsize1,varsize2) )

          call FILE_read( basename, varname, tmp_DP(:,:), rankid=nowrank )

          where( abs( tmp_DP(:,:) - UNDEF8 ) > EPS )
             var(:,:) = real(tmp_DP(:,:),kind=RP)
          elsewhere
             var(:,:) = UNDEF
          endwhere

          deallocate( tmp_DP )

       endif

    else
       LOG_ERROR("SNO_read_bcast_2d",*) '[read_bcast_2d] Unsupported datatype : ', trim(FILE_dtypelist(datatype))
       call PRC_abort
    endif

    call MPI_BCAST( var(:,:), varsize1*varsize2, datatype_mpi, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    return
  end subroutine SNO_read_bcast_2d

  !-----------------------------------------------------------------------------
  subroutine SNO_read_map_2d( &
       basename, &
       nowrank,  &
       nowstep,  &
       varname,  &
       datatype, &
       gin1,     &
       gin2,     &
       gout1,    &
       gout2,    &
       localmap, &
       var       )
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4, &
       UNDEF8 => CONST_UNDEF8, &
       UNDEF  => CONST_UNDEF,  &
       EPS    => CONST_EPS
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8, &
       FILE_dtypelist
    use scale_file, only: &
       FILE_read
    use scale_prc, only: &
       PRC_abort
    use mod_sno_h, only: &
       I_map_p, &
       I_map_i, &
       I_map_j
    implicit none

    character(len=*), intent(in)    :: basename
    integer,          intent(in)    :: nowrank
    integer,          intent(in)    :: nowstep
    character(len=*), intent(in)    :: varname
    integer,          intent(in)    :: datatype
    integer,          intent(in)    :: gin1
    integer,          intent(in)    :: gin2
    integer,          intent(in)    :: gout1
    integer,          intent(in)    :: gout2
    integer(2),       intent(in)    :: localmap(gout1,gout2,3)
    real(RP),         intent(inout) :: var     (gout1,gout2)

    real(SP) :: tmp_SP(gin1,gin2)
    real(DP) :: tmp_DP(gin1,gin2)

    integer  :: i, j, ii, jj
    !---------------------------------------------------------------------------

    if ( datatype == FILE_REAL4 ) then

       call FILE_read( basename, varname, tmp_SP(:,:), step=nowstep, rankid=nowrank )

       do j = 1, gout2
       do i = 1, gout1
          if ( localmap(i,j,I_map_p) == nowrank ) then
             ii = int(localmap(i,j,I_map_i))
             jj = int(localmap(i,j,I_map_j))

             if ( abs( tmp_SP(ii,jj) - UNDEF4 ) > EPS ) then
                var(i,j) = real(tmp_SP(ii,jj),kind=RP)
             else
                var(i,j) = UNDEF
             endif
          endif
       enddo
       enddo

    elseif( datatype == FILE_REAL8 ) then

       call FILE_read( basename, varname, tmp_DP(:,:), step=nowstep, rankid=nowrank )

       do j = 1, gout2
       do i = 1, gout1
          if ( localmap(i,j,I_map_p) == nowrank ) then
             ii = int(localmap(i,j,I_map_i))
             jj = int(localmap(i,j,I_map_j))

             if ( abs( tmp_DP(ii,jj) - UNDEF8 ) > EPS ) then
                var(i,j) = real(tmp_DP(ii,jj),kind=RP)
             else
                var(i,j) = UNDEF
             endif
          endif
       enddo
       enddo

    else
       LOG_ERROR("SNO_read_map_2d",*) '[read_map_2d] Unsupported datatype : ', trim(FILE_dtypelist(datatype))
       call PRC_abort
    endif

    return
  end subroutine SNO_read_map_2d

  !-----------------------------------------------------------------------------
  subroutine SNO_read_map_3d( &
       basename,  &
       nowrank,   &
       nowstep,   &
       varname,   &
       datatype,  &
       transpose, &
       gin1,      &
       gin2,      &
       gin3,      &
       gout1,     &
       gout2,     &
       gout3,     &
       localmap,  &
       var        )
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4, &
       UNDEF8 => CONST_UNDEF8, &
       UNDEF  => CONST_UNDEF,  &
       EPS    => CONST_EPS
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8, &
       FILE_dtypelist
    use scale_file, only: &
       FILE_read
    use scale_prc, only: &
       PRC_abort
    use mod_sno_h, only: &
       I_map_p, &
       I_map_i, &
       I_map_j
    implicit none

    character(len=*), intent(in)    :: basename
    integer,          intent(in)    :: nowrank
    integer,          intent(in)    :: nowstep
    character(len=*), intent(in)    :: varname
    integer,          intent(in)    :: datatype
    logical,          intent(in)    :: transpose
    integer,          intent(in)    :: gin1
    integer,          intent(in)    :: gin2
    integer,          intent(in)    :: gin3
    integer,          intent(in)    :: gout1
    integer,          intent(in)    :: gout2
    integer,          intent(in)    :: gout3
    integer(2),       intent(in)    :: localmap(gout2,gout3,3)
    real(RP),         intent(inout) :: var     (gout1,gout2,gout3)

    real(SP), allocatable :: tmp_SP(:,:,:)
    real(DP), allocatable :: tmp_DP(:,:,:)

    integer  :: k, i, j, ii, jj
    !---------------------------------------------------------------------------

    if ( datatype == FILE_REAL4 ) then

       if ( transpose ) then
          allocate( tmp_SP(gin2,gin3,gin1) ) ! x,y,z
       else
          allocate( tmp_SP(gin1,gin2,gin3) ) ! z,x,y
       endif

       call FILE_read( basename, varname, tmp_SP(:,:,:), step=nowstep, rankid=nowrank )

       if ( transpose ) then
          do j = 1, gout3
          do i = 1, gout2
             if ( localmap(i,j,I_map_p) == nowrank ) then
                ii = int(localmap(i,j,I_map_i))
                jj = int(localmap(i,j,I_map_j))

                do k = 1, gout1
                   if ( abs( tmp_SP(ii,jj,k) - UNDEF4 ) > EPS ) then
                      var(k,i,j) = real(tmp_SP(ii,jj,k),kind=RP)
                   else
                      var(k,i,j) = UNDEF
                   endif
                enddo
             endif
          enddo
          enddo
       else
          do j = 1, gout3
          do i = 1, gout2
             if ( localmap(i,j,I_map_p) == nowrank ) then
                ii = int(localmap(i,j,I_map_i))
                jj = int(localmap(i,j,I_map_j))

                do k = 1, gout1
                   if ( abs( tmp_SP(k,ii,jj) - UNDEF4 ) > EPS ) then
                      var(k,i,j) = real(tmp_SP(k,ii,jj),kind=RP)
                   else
                      var(k,i,j) = UNDEF
                   endif
                enddo
             endif
          enddo
          enddo
       endif

       deallocate( tmp_SP )

    elseif( datatype == FILE_REAL8 ) then

       if ( transpose ) then
          allocate( tmp_DP(gin2,gin3,gin1) ) ! x,y,z
       else
          allocate( tmp_DP(gin1,gin2,gin3) ) ! z,x,y
       endif

       call FILE_read( basename, varname, tmp_DP(:,:,:), step=nowstep, rankid=nowrank )

       if ( transpose ) then
          do j = 1, gout3
          do i = 1, gout2
             if ( localmap(i,j,I_map_p) == nowrank ) then
                ii = int(localmap(i,j,I_map_i))
                jj = int(localmap(i,j,I_map_j))

                do k = 1, gout1
                   if ( abs( tmp_DP(ii,jj,k) - UNDEF8 ) > EPS ) then
                      var(k,i,j) = real(tmp_DP(ii,jj,k),kind=RP)
                   else
                      var(k,i,j) = UNDEF
                   endif
                enddo
             endif
          enddo
          enddo
       else
          do j = 1, gout3
          do i = 1, gout2
             if ( localmap(i,j,I_map_p) == nowrank ) then
                ii = int(localmap(i,j,I_map_i))
                jj = int(localmap(i,j,I_map_j))

                do k = 1, gout1
                   if ( abs( tmp_DP(k,ii,jj) - UNDEF8 ) > EPS ) then
                      var(k,i,j) = real(tmp_DP(k,ii,jj),kind=RP)
                   else
                      var(k,i,j) = UNDEF
                   endif
                enddo
             endif
          enddo
          enddo
       endif

       deallocate( tmp_DP )

    else
       LOG_ERROR("SNO_read_map_3d",*) '[read_map_3d] Unsupported datatype : ', trim(FILE_dtypelist(datatype))
       call PRC_abort
    endif

    return
  end subroutine SNO_read_map_3d

  !-----------------------------------------------------------------------------
  subroutine SNO_attributes_write( &
       fid,          &
       nowrank,      &
       nprocs_x_out, &
       nprocs_y_out, &
       nhalos_x,     &
       nhalos_y,     &
       hinfo,        &
       debug         )
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_set_attribute,         &
       FILE_add_associatedVariable
    use scale_const, &
       UNDEF => CONST_UNDEF
    use scale_atmos_grid_cartesc, only: &
       ATMOS_GRID_CARTESC_NAME
    use scale_file_cartesC, only: &
       FILE_CARTESC_put_globalAttributes, &
       axisattinfo
    use mod_sno_h, only: &
       commoninfo, &
       iteminfo
    implicit none

    integer,          intent(in)  :: fid
    integer,          intent(in)  :: nowrank                               ! current rank                       (output)
    integer,          intent(in)  :: nprocs_x_out                          ! x length of 2D processor topology  (output)
    integer,          intent(in)  :: nprocs_y_out                          ! y length of 2D processor topology  (output)
    integer,          intent(in)  :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)  :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)  :: hinfo                                 ! common information                 (input)
    logical,          intent(in)  :: debug

    type(axisattinfo) :: ainfo(4)

    integer :: rankidx(2) ! my rank in 2D process topology (x:y)
    integer :: IMAX, JMAX
    !---------------------------------------------------------------------------

    rankidx(1) = mod(nowrank,nprocs_x_out)
    rankidx(2) =     nowrank/nprocs_x_out

    select case( hinfo%grid_name )
    case( ATMOS_GRID_CARTESC_NAME )

       call FILE_CARTESC_put_globalAttributes( fid,                                                     &
                                               rankidx(1), rankidx(2),                                  &
                                               nprocs_x_out, nprocs_y_out,                              &
                                               hinfo%periodic(2), hinfo%periodic(3),                    &
                                               hinfo%gridsize(1),                                       &
                                               hinfo%gridsize(4), hinfo%gridsize(5), hinfo%gridsize(6), &
                                               hinfo%gridsize(2), hinfo%gridsize(3),                    &
                                               hinfo%halosize(1), hinfo%halosize(2), hinfo%halosize(3), &
                                               hinfo%time_start(1), hinfo%time_units, hinfo%calendar    )

    case default
       LOG_WARN("SNO_attributes_write",*) 'invalud grid_name: ', trim(hinfo%grid_name)
       call PRC_abort
    end select

    IMAX = hinfo%gridsize(2) / nprocs_x_out
    JMAX = hinfo%gridsize(3) / nprocs_y_out

    ! for x
    ainfo(1)%size_global (1) = hinfo%xatt_size_global(1)
    ainfo(1)%start_global(1) = hinfo%xatt_halo_global(1) + rankidx(1) * IMAX + 1
    ainfo(1)%halo_global (1) = hinfo%xatt_halo_global(1) ! west side
    ainfo(1)%halo_global (2) = hinfo%xatt_halo_global(2) ! east side
    ainfo(1)%halo_local  (1) = 0                         ! west side
    ainfo(1)%halo_local  (2) = 0                         ! east side
    ainfo(1)%periodic        = hinfo%periodic(2)
    if( rankidx(1) == 0              ) ainfo(1)%start_global(1) = ainfo(1)%start_global(1) - nhalos_x
    if( rankidx(1) == 0              ) ainfo(1)%halo_local  (1) = nhalos_x
    if( rankidx(1) == nprocs_x_out-1 ) ainfo(1)%halo_local  (2) = nhalos_x

    ! for xh
    ainfo(2) = ainfo(1)
    if ( (.NOT. ainfo(2)%periodic) .AND. ainfo(2)%halo_global(1) == 0 ) then
       ainfo(2)%size_global(1) = ainfo(2)%size_global(1) + 1
       ainfo(2)%halo_global(1) = ainfo(2)%halo_global(1) + 1
!        if ( rankidx(1) == 0 ) then
!           ainfo(2)%start_global(1) = ainfo(2)%start_global(1) + 1
!        else
          ainfo(2)%halo_local  (1) = ainfo(2)%halo_local  (1) + 1
!        endif
    endif

    ! for y
    ainfo(3)%size_global (1) = hinfo%yatt_size_global(1)
    ainfo(3)%start_global(1) = hinfo%yatt_halo_global(1) + rankidx(2) * JMAX + 1
    ainfo(3)%halo_global (1) = hinfo%yatt_halo_global(1) ! west side
    ainfo(3)%halo_global (2) = hinfo%yatt_halo_global(2) ! east side
    ainfo(3)%halo_local  (1) = 0                         ! west side
    ainfo(3)%halo_local  (2) = 0                         ! east side
    ainfo(3)%periodic        = hinfo%periodic(3)
    if( rankidx(2) == 0              ) ainfo(3)%start_global(1) = ainfo(3)%start_global(1) - nhalos_y
    if( rankidx(2) == 0              ) ainfo(3)%halo_local  (1) = nhalos_y
    if( rankidx(2) == nprocs_y_out-1 ) ainfo(3)%halo_local  (2) = nhalos_y

    ! for yh
    ainfo(4) = ainfo(3)
    if ( (.NOT. ainfo(4)%periodic) .AND. ainfo(4)%halo_global(1) == 0 ) then
       ainfo(4)%size_global(1) = ainfo(4)%size_global(1) + 1
       ainfo(4)%halo_global(1) = ainfo(4)%halo_global(1) + 1
!        if ( rankidx(2) == 0 ) then
!           ainfo(4)%start_global(1) = ainfo(4)%start_global(1) + 1
!        else
          ainfo(4)%halo_local  (1) = ainfo(4)%halo_local  (1) + 1
!        endif
    endif

    call FILE_set_attribute( fid, "x" , "size_global" , ainfo(1)%size_global (:) )
    call FILE_set_attribute( fid, "x" , "start_global", ainfo(1)%start_global(:) )
    call FILE_set_attribute( fid, "x" , "halo_global" , ainfo(1)%halo_global (:) )
    call FILE_set_attribute( fid, "x" , "halo_local"  , ainfo(1)%halo_local  (:) )
    call FILE_set_attribute( fid, "x" , "periodic"    , ainfo(1)%periodic        )

    call FILE_set_attribute( fid, "xh", "size_global" , ainfo(2)%size_global (:) )
    call FILE_set_attribute( fid, "xh", "start_global", ainfo(2)%start_global(:) )
    call FILE_set_attribute( fid, "xh", "halo_global" , ainfo(2)%halo_global (:) )
    call FILE_set_attribute( fid, "xh", "halo_local"  , ainfo(2)%halo_local  (:) )
    call FILE_set_attribute( fid, "xh", "periodic"    , ainfo(2)%periodic        )

    call FILE_set_attribute( fid, "y" , "size_global" , ainfo(3)%size_global (:) )
    call FILE_set_attribute( fid, "y" , "start_global", ainfo(3)%start_global(:) )
    call FILE_set_attribute( fid, "y" , "halo_global" , ainfo(3)%halo_global (:) )
    call FILE_set_attribute( fid, "y" , "halo_local"  , ainfo(3)%halo_local  (:) )
    call FILE_set_attribute( fid, "y" , "periodic"    , ainfo(3)%periodic        )

    call FILE_set_attribute( fid, "yh", "size_global" , ainfo(4)%size_global (:) )
    call FILE_set_attribute( fid, "yh", "start_global", ainfo(4)%start_global(:) )
    call FILE_set_attribute( fid, "yh", "halo_global" , ainfo(4)%halo_global (:) )
    call FILE_set_attribute( fid, "yh", "halo_local"  , ainfo(4)%halo_local  (:) )
    call FILE_set_attribute( fid, "yh", "periodic"    , ainfo(4)%periodic        )

    if ( hinfo%minfo_mapping_name /= "" ) then
       call FILE_set_attribute( fid, "x" , "standard_name", "projection_x_coordinate" )
       call FILE_set_attribute( fid, "xh", "standard_name", "projection_x_coordinate" )
       call FILE_set_attribute( fid, "y" , "standard_name", "projection_y_coordinate" )
       call FILE_set_attribute( fid, "yh", "standard_name", "projection_y_coordinate" )

       call FILE_add_associatedVariable( fid, hinfo%minfo_mapping_name )
       call FILE_set_attribute( fid, hinfo%minfo_mapping_name, "grid_mapping_name",  hinfo%minfo_mapping_name )

       if ( hinfo%minfo_false_easting(1) /= UNDEF ) then
          call FILE_set_attribute( fid,                         & ! [IN]
                                   hinfo%minfo_mapping_name,    & ! [IN]
                                   "false_easting",             & ! [IN]
                                   hinfo%minfo_false_easting(:) ) ! [IN]
       endif

       if ( hinfo%minfo_false_northing(1) /= UNDEF ) then
          call FILE_set_attribute( fid,                          & ! [IN]
                                   hinfo%minfo_mapping_name,     & ! [IN]
                                   "false_northing",             & ! [IN]
                                   hinfo%minfo_false_northing(:) ) ! [IN]
       endif

       if ( hinfo%minfo_longitude_of_central_meridian(1) /= UNDEF ) then
          call FILE_set_attribute( fid,                                         & ! [IN]
                                   hinfo%minfo_mapping_name,                    & ! [IN]
                                   "longitude_of_central_meridian",             & ! [IN]
                                   hinfo%minfo_longitude_of_central_meridian(:) ) ! [IN]
       endif

       if ( hinfo%minfo_longitude_of_projection_origin(1) /= UNDEF ) then
          call FILE_set_attribute( fid,                                          & ! [IN]
                                   hinfo%minfo_mapping_name,                     & ! [IN]
                                   "longitude_of_projection_origin",             & ! [IN]
                                   hinfo%minfo_longitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( hinfo%minfo_latitude_of_projection_origin(1) /= UNDEF ) then
          call FILE_set_attribute( fid,                                         & ! [IN]
                                   hinfo%minfo_mapping_name,                    & ! [IN]
                                   "latitude_of_projection_origin",             & ! [IN]
                                   hinfo%minfo_latitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( hinfo%minfo_straight_vertical_longitude_from_pole(1) /= UNDEF ) then
          call FILE_set_attribute( fid,                                                 & ! [IN]
                                   hinfo%minfo_mapping_name,                            & ! [IN]
                                   "straight_vertical_longitude_from_pole",             & ! [IN]
                                   hinfo%minfo_straight_vertical_longitude_from_pole(:) ) ! [IN]
       endif

       if ( hinfo%minfo_standard_parallel(1) /= UNDEF ) then
          if ( hinfo%minfo_standard_parallel(2) /= UNDEF ) then
             call FILE_set_attribute( fid,                               & ! [IN]
                                      hinfo%minfo_mapping_name,          & ! [IN]
                                      "standard_parallel",               & ! [IN]
                                      hinfo%minfo_standard_parallel(1:2) ) ! [IN]
          else
             call FILE_set_attribute( fid,                               & ! [IN]
                                      hinfo%minfo_mapping_name,          & ! [IN]
                                      "standard_parallel",               & ! [IN]
                                      hinfo%minfo_standard_parallel(1:1) ) ! [IN]
          endif
       endif
    endif

    ! SGRID
    call FILE_add_associatedVariable( fid, "grid" )
    call FILE_set_attribute( fid, "grid", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid", "node_dimensions",     "xh yh" )
    call FILE_set_attribute( fid, "grid", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_set_attribute( fid, "grid", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_set_attribute( fid, "grid", "face_coordinates",    "lon lat" )
    call FILE_set_attribute( fid, "grid", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_set_attribute( fid, "grid", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_set_attribute( fid, "grid", "vertical_dimensions", "z: zh (padding: none)" )

    call FILE_add_associatedVariable( fid, "grid_ocean" )
    call FILE_set_attribute( fid, "grid_ocean", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_ocean", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_ocean", "node_dimensions",     "xh yh" )
    call FILE_set_attribute( fid, "grid_ocean", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_set_attribute( fid, "grid_ocean", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_set_attribute( fid, "grid_ocean", "face_coordinates",    "lon lat" )
    call FILE_set_attribute( fid, "grid_ocean", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_set_attribute( fid, "grid_ocean", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_set_attribute( fid, "grid_ocean", "vertical_dimensions", "oz: ozh (padding: none)" )

    call FILE_add_associatedVariable( fid, "grid_land" )
    call FILE_set_attribute( fid, "grid_land", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_land", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_land", "node_dimensions",     "xh yh" )
    call FILE_set_attribute( fid, "grid_land", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_set_attribute( fid, "grid_land", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_set_attribute( fid, "grid_land", "face_coordinates",    "lon lat" )
    call FILE_set_attribute( fid, "grid_land", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_set_attribute( fid, "grid_land", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_set_attribute( fid, "grid_land", "vertical_dimensions", "lz: lzh (padding: none)" )

    call FILE_add_associatedVariable( fid, "grid_urban" )
    call FILE_set_attribute( fid, "grid_urban", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_urban", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_urban", "node_dimensions",     "xh yh" )
    call FILE_set_attribute( fid, "grid_urban", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_set_attribute( fid, "grid_urban", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_set_attribute( fid, "grid_urban", "face_coordinates",    "lon lat" )
    call FILE_set_attribute( fid, "grid_urban", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_set_attribute( fid, "grid_urban", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_set_attribute( fid, "grid_urban", "vertical_dimensions", "uz: uzh (padding: none)" )

    call FILE_add_associatedVariable( fid, "grid_pressure" )
    call FILE_set_attribute( fid, "grid_pressure", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_pressure", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_pressure", "node_dimensions",     "xh yh" )
    call FILE_set_attribute( fid, "grid_pressure", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_set_attribute( fid, "grid_pressure", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_set_attribute( fid, "grid_pressure", "face_coordinates",    "lon lat" )
    call FILE_set_attribute( fid, "grid_pressure", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_set_attribute( fid, "grid_pressure", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_set_attribute( fid, "grid_pressure", "vertical_dimensions", "pressure" )

    call FILE_add_associatedVariable( fid, "grid_z" )
    call FILE_set_attribute( fid, "grid_z", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_z", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_z", "node_dimensions",     "xh yh" )
    call FILE_set_attribute( fid, "grid_z", "face_dimensions",     "x: xh (padding: none) y: yh (padding: none)" )
    call FILE_set_attribute( fid, "grid_z", "node_coordinates",    "lon_uv lat_uv" )
    call FILE_set_attribute( fid, "grid_z", "face_coordinates",    "lon lat" )
    call FILE_set_attribute( fid, "grid_z", "edge1_coordinates",   "lon_uy lat_uy" )
    call FILE_set_attribute( fid, "grid_z", "edge2_coordinates",   "lon_xv lat_xv" )
    call FILE_set_attribute( fid, "grid_z", "vertical_dimensions", "height_xyw: height (padding: none)" )

    call FILE_add_associatedVariable( fid, "grid_model" )
    call FILE_set_attribute( fid, "grid_model", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_model", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_model", "node_dimensions",     "FX FY" )
    call FILE_set_attribute( fid, "grid_model", "face_dimensions",     "CX: FY (padding: none) CY: FY (padding: none)" )
    call FILE_set_attribute( fid, "grid_model", "vertical_dimensions", "CZ: FZ (padding: none)" )

    call FILE_add_associatedVariable( fid, "grid_model_global" )
    call FILE_set_attribute( fid, "grid_model_global", "cf_role",             "grid_topology" )
    call FILE_set_attribute( fid, "grid_model_global", "topology_dimension",  (/ 2 /) )
    call FILE_set_attribute( fid, "grid_model_global", "node_dimensions",     "FXG FYG" )
    call FILE_set_attribute( fid, "grid_model_global", "face_dimensions",     "CXG: FYG (padding: none) CYG: FYG (padding: none)" )
    call FILE_set_attribute( fid, "grid_model_global", "vertical_dimensions", "CZ: FZ (padding: none)" )

    return
  end subroutine SNO_attributes_write

end module mod_sno
