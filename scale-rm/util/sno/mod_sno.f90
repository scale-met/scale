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
module mod_sno
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
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
  public :: SNO_calc_localsize
  public :: SNO_read_bcast_1d
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
    use scale_process, only: &
       PRC_MPIstop
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
       write(*,*) 'xxx [SNO_proc_alloc] # of using PEs should be less than # of files for output. Check'
       call PRC_MPIstop
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_proc_alloc] Process management'
    if( IO_L ) write(IO_FID_LOG,*) '*** # of PEs (this execution)     : ', nprocs
    if( IO_L ) write(IO_FID_LOG,*) '*** Rank     (this execution)     : ', myrank
    if( IO_L ) write(IO_FID_LOG,*) '*** Master?  (this execution)     : ', ismaster
    if( IO_L ) write(IO_FID_LOG,*) '*** # of PEs (output file)        : ', nprocs_out, '(', nprocs_x_out, 'x', nprocs_y_out, ')'
    if( IO_L ) write(IO_FID_LOG,*) '*** Start PE managed by this rank : ', pstr
    if( IO_L ) write(IO_FID_LOG,*) '*** End   PE managed by this rank : ', pend

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
       ngrids_x,     &
       ngrids_y,     &
       nhalos_x,     &
       nhalos_y,     &
       hinfo,        &
       naxis,        &
       axisname,     &
       nvars,        &
       varname,      &
       debug         )
    use mpi
    use gtool_file_h, only: &
       File_FREAD
    use gtool_file, only: &
       FileOpen,               &
       FileGetCommoninfo,      &
       FileGetGlobalAttribute, &
       FileGetAttribute
    use scale_process, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_MPIstop
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
    integer,                intent(out) :: ngrids_x                 ! size of x-axis grids              (global,sometimes including halo)
    integer,                intent(out) :: ngrids_y                 ! size of y-axis grids              (global,sometimes including halo)
    integer,                intent(out) :: nhalos_x                 ! size of x-axis halo grids         (global,sometimes have a size)
    integer,                intent(out) :: nhalos_y                 ! size of y-axis halo grids         (global,sometimes have a size)
    type(commoninfo),       intent(out) :: hinfo                    ! common information                (input)
    integer,                intent(out) :: naxis                    ! number of axis variables          (input)
    character(len=H_SHORT), intent(out) :: axisname(item_limit)     ! name   of axis variables          (input)
    integer,                intent(out) :: nvars                    ! number of variables               (input)
    character(len=H_SHORT), intent(out) :: varname (item_limit)     ! name   of variables               (input)
    logical,                intent(in)  :: debug

    integer                :: procsize(2)                   ! total process size        (x:y)
    character(len=6)       :: periodic(3)

    integer                :: nvars_file                    ! number of variables from input file
    character(len=H_SHORT) :: varname_file(item_limit) = '' ! name   of variables from input file
    integer                :: nvars_req                     ! number of variables, requested
    logical                :: exist       (item_limit)      ! check flag for requested variables

    integer :: nprocs_in       ! number of peXXXXXX files           (input)

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
    !---------------------------------------------------------------------------

    !---  read info from global file

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [SNO_file_getinfo] Get process & horizontal grid info from input file'
    if( IO_L ) write(IO_FID_LOG,*) '*** Basename of input file : ', trim(basename)

    if ( ismaster ) then
       nowrank = 0 ! first file

       call FileOpen( fid,             & ! [OUT]
                      basename,        & ! [IN]
                      File_FREAD,      & ! [IN]
                      myrank = nowrank ) ! [IN]

       call FileGetCommoninfo( fid,             & ! [IN]
                               nowrank,         & ! [IN]
                               item_limit,      & ! [IN]
                               hinfo%title,     & ! [OUT]
                               hinfo%source,    & ! [OUT]
                               hinfo%institute, & ! [OUT]
                               nvars_file,      & ! [OUT]
                               varname_file(:)  ) ! [OUT]

       call FileGetGlobalAttribute( fid, "scale_rm_prc_num_x", procsize(1:1) )
       call FileGetGlobalAttribute( fid, "scale_rm_prc_num_y", procsize(2:2) )

       call FileGetGlobalAttribute( fid, "scale_rm_prc_periodic_z",   periodic(1) )
       call FileGetGlobalAttribute( fid, "scale_rm_prc_periodic_x",   periodic(2) )
       call FileGetGlobalAttribute( fid, "scale_rm_prc_periodic_y",   periodic(3) )
       hinfo%periodic(1) = trim(periodic(1))
       hinfo%periodic(2) = trim(periodic(2))
       hinfo%periodic(3) = trim(periodic(3))

       call FileGetGlobalAttribute( fid, "scale_rm_grid_index_kmax",  hinfo%gridsize(1:1) )
       call FileGetGlobalAttribute( fid, "scale_rm_grid_index_imaxg", hinfo%gridsize(2:2) )
       call FileGetGlobalAttribute( fid, "scale_rm_grid_index_jmaxg", hinfo%gridsize(3:3) )

       call FileGetGlobalAttribute( fid, "scale_rm_grid_index_khalo", hinfo%halosize(1:1) )
       call FileGetGlobalAttribute( fid, "scale_rm_grid_index_ihalo", hinfo%halosize(2:2) )
       call FileGetGlobalAttribute( fid, "scale_rm_grid_index_jhalo", hinfo%halosize(3:3) )

       call FileGetGlobalAttribute( fid, "time_units", hinfo%time_units    )
       call FileGetGlobalAttribute( fid, "time_start", hinfo%time_start(:) )

       call FileGetAttribute( fid, 'x', 'size_global', hinfo%xatt_size_global(:) )
       call FileGetAttribute( fid, 'x', 'halo_global', hinfo%xatt_halo_global(:) )

       call FileGetAttribute( fid, 'y', 'size_global', hinfo%yatt_size_global(:) )
       call FileGetAttribute( fid, 'y', 'halo_global', hinfo%yatt_halo_global(:) )
    endif

    call MPI_BCAST( hinfo%title              , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%source             , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%institute          , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%periodic        (1), 5*3  , MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%gridsize        (1), 3    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%halosize        (1), 3    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%time_units         , H_MID, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%time_start      (1), 1    , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%xatt_size_global(1), 1    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%xatt_halo_global(1), 2    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%yatt_size_global(1), 1    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%yatt_halo_global(1), 2    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( procsize              (1), 2    , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    call MPI_BCAST( nvars_file     , 1                 , MPI_INTEGER,   PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( varname_file(1), H_SHORT*item_limit, MPI_CHARACTER, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%title           ", trim(hinfo%title)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%source          ", trim(hinfo%source)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%institute       ", trim(hinfo%institute)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%periodic        ", hinfo%periodic        (:)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%gridsize        ", hinfo%gridsize        (:)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%halosize        ", hinfo%halosize        (:)
       if( IO_L ) write(IO_FID_LOG,*) "procsize              ", procsize              (:)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%xatt_size_global", hinfo%xatt_size_global(:)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%xatt_halo_global", hinfo%xatt_halo_global(:)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%yatt_size_global", hinfo%yatt_size_global(:)
       if( IO_L ) write(IO_FID_LOG,*) "hinfo%yatt_halo_global", hinfo%yatt_halo_global(:)
    endif

    !--- process management

    nprocs_x_in = procsize(1)
    nprocs_y_in = procsize(2)
    nprocs_in   = nprocs_x_in * nprocs_y_in

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** # of PEs (input  file)        : ', nprocs_in, '(', nprocs_x_in, 'x', nprocs_y_in, ')'

    !--- horizontal grid management

    ngrids_x        = hinfo%xatt_size_global(1)
    ngrids_y        = hinfo%yatt_size_global(1)
    nhalos_x        = hinfo%xatt_halo_global(1) ! assume both side have same number
    nhalos_y        = hinfo%yatt_halo_global(1) ! assume both side have same number

    ngrids_x_nohalo = ngrids_x - 2*nhalos_x
    ngrids_y_nohalo = ngrids_y - 2*nhalos_y

    ngrids_x_in     = ngrids_x_nohalo / nprocs_x_in
    ngrids_y_in     = ngrids_y_nohalo / nprocs_y_in

    ngrids_x_out    = ngrids_x_nohalo / nprocs_x_out
    ngrids_y_out    = ngrids_y_nohalo / nprocs_y_out

    ngrids          = ngrids_x_nohalo * ngrids_y_nohalo
    ngrids_in       = ngrids_x_in     * ngrids_y_in
    ngrids_out      = ngrids_x_out    * ngrids_y_out

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Grid management : horizontal'
    if( IO_L ) write(IO_FID_LOG,*) '*** # of total grids           without halo : ', &
                                   ngrids    , '(', ngrids_x_nohalo, 'x', ngrids_y_nohalo, ')'
    if( IO_L ) write(IO_FID_LOG,*) '*** # of grids per input  file without halo : ', &
                                   ngrids_in , '(', ngrids_x_in    , 'x', ngrids_y_in    , ')'
    if( IO_L ) write(IO_FID_LOG,*) '*** # of grids per output file without halo : ', &
                                   ngrids_out, '(', ngrids_x_out   , 'x', ngrids_y_out   , ')'
    if( IO_L ) write(IO_FID_LOG,*) '*** # of halo grids (x-axis,one side)       : ', nhalos_x
    if( IO_L ) write(IO_FID_LOG,*) '*** # of halo grids (y-axis,one side)       : ', nhalos_y

    if ( mod(ngrids_x_nohalo,nprocs_x_out) /= 0 ) then
       write(*,*) 'xxx [SNO_file_getinfo] The allowable case is that # of the total x-grids is divisible with # of the x-PEs for the output. Stop'
       call PRC_MPIstop
    endif

    if ( mod(ngrids_y_nohalo,nprocs_y_out) /= 0 ) then
       write(*,*) 'xxx [SNO_file_getinfo] The allowable case is that # of the total y-grids is divisible with # of the y-PEs for the output. Stop'
       call PRC_MPIstop
    endif

    exist(:) = .false.
    do nn = 1, item_limit
       if( vars(nn) == '' ) exit
    enddo
    nvars_req = nn - 1

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) "*** nvars_req = ", nvars_req
       do nn = 1, nvars_req
          if( IO_L ) write(IO_FID_LOG,*) "*** ", nn, trim(vars(nn))
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) "*** nvars_file = ", nvars_file
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
          if( IO_L ) write(IO_FID_LOG,*) "*** ", n, trim(varname_file(n))
       endif

       select case(varname_file(n))
       case('z','zh','lz','lzh','uz','uzh',                                                      &
            'CZ','FZ','CDZ','FDZ','CBFZ','FBFZ','LCZ','LFZ','LCDZ','UCZ','UFZ','UCDZ',           &
            'x','xh','y','yh',                                                                   &
            'CX','CY','FX','FY','CDX','CDY','FDX','FDY','CBFX','CBFY','FBFX','FBFY',             &
            'CXG','CYG','FXG','FYG','CDXG','CDYG','FDXG','FDYG','CBFXG','CBFYG','FBFXG','FBFYG', &
            'lon','lon_uy','lon_xv','lon_uv','lat','lat_uy','lat_xv','lat_uv','topo','lsmask',   &
            'height','height_wxy','height_xyw','height_uyz','height_xvz',                        &
            'height_uvz','height_uyw','height_xvw','height_uvw'                                  )
          naxis           = naxis + 1
          axisname(naxis) = varname_file(n)
       case('time','time_bnds')
          ! do nothing
       case('lambert_conformal_conic')
          if ( ismaster ) then
             call FileGetAttribute( fid, varname_file(n), "grid_mapping_name",                     hinfo%minfo_mapping_name                             )
             call FileGetAttribute( fid, varname_file(n), "false_easting",                         hinfo%minfo_false_easting                        (:) )
             call FileGetAttribute( fid, varname_file(n), "false_northing",                        hinfo%minfo_false_northing                       (:) )
             call FileGetAttribute( fid, varname_file(n), "longitude_of_central_meridian",         hinfo%minfo_longitude_of_central_meridian        (:) )
             call FileGetAttribute( fid, varname_file(n), "latitude_of_projection_origin",         hinfo%minfo_latitude_of_projection_origin        (:) )
             call FileGetAttribute( fid, varname_file(n), "standard_parallel",                     hinfo%minfo_standard_parallel                    (:) )
          endif
       case('polar_stereographic')
          if ( ismaster ) then
             call FileGetAttribute( fid, varname_file(n), "grid_mapping_name",                     hinfo%minfo_mapping_name                             )
             call FileGetAttribute( fid, varname_file(n), "false_easting",                         hinfo%minfo_false_easting                        (:) )
             call FileGetAttribute( fid, varname_file(n), "false_northing",                        hinfo%minfo_false_northing                       (:) )
             call FileGetAttribute( fid, varname_file(n), "latitude_of_projection_origin",         hinfo%minfo_latitude_of_projection_origin        (:) )
             call FileGetAttribute( fid, varname_file(n), "straight_vertical_longitude_from_pole", hinfo%minfo_straight_vertical_longitude_from_pole(:) )
             call FileGetAttribute( fid, varname_file(n), "standard_parallel",                     hinfo%minfo_standard_parallel                  (1:1) )
          endif
       case('mercator')
          if ( ismaster ) then
             call FileGetAttribute( fid, varname_file(n), "grid_mapping_name",                     hinfo%minfo_mapping_name                             )
             call FileGetAttribute( fid, varname_file(n), "false_easting",                         hinfo%minfo_false_easting                        (:) )
             call FileGetAttribute( fid, varname_file(n), "false_northing",                        hinfo%minfo_false_northing                       (:) )
             call FileGetAttribute( fid, varname_file(n), "longitude_of_projection_origin",        hinfo%minfo_longitude_of_projection_origin       (:) )
          endif
       case('equirectangular')
          if ( ismaster ) then
             call FileGetAttribute( fid, varname_file(n), "grid_mapping_name",                     hinfo%minfo_mapping_name                             )
             call FileGetAttribute( fid, varname_file(n), "false_easting",                         hinfo%minfo_false_easting                        (:) )
             call FileGetAttribute( fid, varname_file(n), "false_northing",                        hinfo%minfo_false_northing                       (:) )
             call FileGetAttribute( fid, varname_file(n), "longitude_of_central_meridian",         hinfo%minfo_longitude_of_central_meridian        (:) )
          endif
       case default
          if ( nvars_req == 0 ) then
             nvars           = nvars + 1
             varname (nvars) = varname_file(n)
          else
             do nn = 1, nvars_req
                if ( varname_file(n) == vars(nn) ) then
                   if ( exist(nn) ) then
                      if( IO_L ) write(IO_FID_LOG,*) 'xxx [SNO_file_getinfo] variable ', trim(vars(nn)), &
                                                     ' is requested two times. check namelist!'
                      call PRC_MPIstop
                   endif

                   nvars          = nvars + 1
                   varname(nvars) = varname_file(n)
                   exist(nn)      = .true.
                endif
             enddo
          endif
       end select

    enddo

    if ( nvars_req > 0 .AND. nvars /= nvars_req ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [SNO_file_getinfo] some requested variables are missing.', nvars, nvars_req
       if( IO_L ) write(IO_FID_LOG,*) 'xxx nvars = ', nvars, ', nvars_req = ', nvars_req
       do nn = 1, nvars_req
          if( IO_L ) write(IO_FID_LOG,*) 'xxx check:', exist(nn), ', name: ', trim(vars(nn))
       enddo
       call PRC_MPIstop
    endif

    call MPI_BCAST( hinfo%minfo_mapping_name                            , H_SHORT, MPI_CHARACTER,        PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_false_easting                        (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_false_northing                       (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_longitude_of_central_meridian        (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_longitude_of_projection_origin       (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_latitude_of_projection_origin        (1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_straight_vertical_longitude_from_pole(1), 1      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
    call MPI_BCAST( hinfo%minfo_standard_parallel                    (1), 2      , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Axis to read'
    do n = 1, naxis
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A32)') '*** No.', n, ', name: ', trim(axisname(n))
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Variable to read'
    do n = 1, nvars
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A32)') '*** No.', n, ', name: ', trim(varname(n))
    enddo

    return
  end subroutine SNO_file_getinfo

  !-----------------------------------------------------------------------------
  subroutine SNO_calc_localsize( &
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
    if ( nhalos_y == 0 .AND. py == 1 .AND. hinfo%periodic(3) == 'false' ) then
       ngrids_yh = ngrids_yh + 1
    endif

    ngrids_x = ( ngrids_x_total - 2*nhalos_x ) / nprocs_x
    if( px == 1        ) ngrids_x = ngrids_x + nhalos_x
    if( px == nprocs_x ) ngrids_x = ngrids_x + nhalos_x

    ngrids_xh = ngrids_x
    if ( nhalos_x == 0 .AND. px == 1 .AND. hinfo%periodic(2) == 'false' ) then
       ngrids_xh = ngrids_xh + 1
    endif

    return
  end subroutine SNO_calc_localsize

  !-----------------------------------------------------------------------------
  subroutine SNO_read_bcast_1d( &
       ismaster,  &
       basename,  &
       varname,   &
       datatype,  &
       varsize,   &
       var        )
    use mpi
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8, &
       File_dtypelist
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_MPIstop
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

    if ( datatype == File_REAL4 ) then

       if ( ismaster ) then
          allocate( tmp_SP(varsize) )

          call FileRead( tmp_SP(:), basename, varname, 1, nowrank )

          var(:) = real(tmp_SP(:),kind=RP)

          deallocate( tmp_SP )
       endif

    elseif( datatype == File_REAL8 ) then

       if ( ismaster ) then
          allocate( tmp_DP(varsize) )

          call FileRead( tmp_DP(:), basename, varname, 1, nowrank )

          var(:) = real(tmp_DP(:),kind=RP)

          deallocate( tmp_DP )
       endif

    else
       write(*,*) 'xxx [read_bcast_1d] Unsupported datatype : ', trim(File_dtypelist(datatype))
       call PRC_MPIstop
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
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8, &
       File_dtypelist
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_MPIstop
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

    if ( datatype == File_REAL4 ) then

       call FileRead( tmp_SP(:), basename, varname, nowstep, nowrank )

       do i = 1, gout1
          if ( localmap(i,I_map_p) == nowrank ) then
             ii = int(localmap(i,I_map_i))

             var(i) = real(tmp_SP(ii),kind=RP)
          endif
       enddo

    elseif( datatype == File_REAL8 ) then

       call FileRead( tmp_DP(:), basename, varname, nowstep, nowrank )

       do i = 1, gout1
          if ( localmap(i,I_map_p) == nowrank ) then
             ii = int(localmap(i,I_map_i))

             var(i) = real(tmp_DP(ii),kind=RP)
          endif
       enddo

    else
       write(*,*) 'xxx [read_map_1d] Unsupported datatype : ', trim(File_dtypelist(datatype))
       call PRC_MPIstop
    endif

    return
  end subroutine SNO_read_map_1d

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
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8, &
       File_dtypelist
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_MPIstop
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

    if ( datatype == File_REAL4 ) then

       call FileRead( tmp_SP(:,:), basename, varname, nowstep, nowrank )

       do j = 1, gout2
       do i = 1, gout1
          if ( localmap(i,j,I_map_p) == nowrank ) then
             ii = int(localmap(i,j,I_map_i))
             jj = int(localmap(i,j,I_map_j))

             var(i,j) = real(tmp_SP(ii,jj),kind=RP)
          endif
       enddo
       enddo

    elseif( datatype == File_REAL8 ) then

       call FileRead( tmp_DP(:,:), basename, varname, nowstep, nowrank )

       do j = 1, gout2
       do i = 1, gout1
          if ( localmap(i,j,I_map_p) == nowrank ) then
             ii = int(localmap(i,j,I_map_i))
             jj = int(localmap(i,j,I_map_j))

             var(i,j) = real(tmp_DP(ii,jj),kind=RP)
          endif
       enddo
       enddo

    else
       write(*,*) 'xxx [read_map_2d] Unsupported datatype : ', trim(File_dtypelist(datatype))
       call PRC_MPIstop
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
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8, &
       File_dtypelist
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_MPIstop
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

    if ( datatype == File_REAL4 ) then

       if ( transpose ) then
          allocate( tmp_SP(gin2,gin3,gin1) ) ! x,y,z
       else
          allocate( tmp_SP(gin1,gin2,gin3) ) ! z,x,y
       endif

       call FileRead( tmp_SP(:,:,:), basename, varname, nowstep, nowrank )

       if ( transpose ) then
          do j = 1, gout3
          do i = 1, gout2
             if ( localmap(i,j,I_map_p) == nowrank ) then
                ii = int(localmap(i,j,I_map_i))
                jj = int(localmap(i,j,I_map_j))

                do k = 1, gout1
                   var(k,i,j) = real(tmp_SP(ii,jj,k),kind=RP)
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
                   var(k,i,j) = real(tmp_SP(k,ii,jj),kind=RP)
                enddo
             endif
          enddo
          enddo
       endif

       deallocate( tmp_SP )

    elseif( datatype == File_REAL8 ) then

       if ( transpose ) then
          allocate( tmp_DP(gin2,gin3,gin1) ) ! x,y,z
       else
          allocate( tmp_DP(gin1,gin2,gin3) ) ! z,x,y
       endif

       call FileRead( tmp_DP(:,:,:), basename, varname, nowstep, nowrank )

       if ( transpose ) then
          do j = 1, gout3
          do i = 1, gout2
             if ( localmap(i,j,I_map_p) == nowrank ) then
                ii = int(localmap(i,j,I_map_i))
                jj = int(localmap(i,j,I_map_j))

                do k = 1, gout1
                   var(k,i,j) = real(tmp_DP(ii,jj,k),kind=RP)
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
                   var(k,i,j) = real(tmp_DP(k,ii,jj),kind=RP)
                enddo
             endif
          enddo
          enddo
       endif

       deallocate( tmp_DP )

    else
       write(*,*) 'xxx [read_map_3d] Unsupported datatype : ', trim(File_dtypelist(datatype))
       call PRC_MPIstop
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
       dinfo,        &
       debug         )
    use gtool_file, only: &
       FileSetGlobalAttribute,   &
       FileSetAttribute,         &
       FileAddAssociatedVariable
    use scale_const, &
       UNDEF => CONST_UNDEF
    use scale_fileio, only: &
       axisattinfo, &
       mappinginfo
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
    type(iteminfo),   intent(in)  :: dinfo                                 ! variable information               (input)
    logical,          intent(in)  :: debug

    type(axisattinfo) :: ainfo(4)
    type(mappinginfo) :: minfo

    integer :: rankidx(2) ! my rank in 2D process topology (x:y)
    integer :: IMAX, JMAX
    !---------------------------------------------------------------------------

    rankidx(1) = mod(nowrank,nprocs_x_out)
    rankidx(2) =     nowrank/nprocs_x_out

    call FileSetGlobalAttribute( fid, "scale_rm_prc_rank_x", rankidx(1:1) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_prc_rank_y", rankidx(2:2) ) ! [IN]

    call FileSetGlobalAttribute( fid, "scale_rm_prc_num_x",  (/nprocs_x_out/) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_prc_num_y",  (/nprocs_y_out/) ) ! [IN]

    call FileSetGlobalAttribute( fid, "scale_rm_prc_periodic_z",   hinfo%periodic(1) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_prc_periodic_x",   hinfo%periodic(2) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_prc_periodic_y",   hinfo%periodic(3) ) ! [IN]

    call FileSetGlobalAttribute( fid, "scale_rm_grid_index_kmax",  hinfo%gridsize(1:1) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_grid_index_imaxg", hinfo%gridsize(2:2) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_grid_index_jmaxg", hinfo%gridsize(3:3) ) ! [IN]

    call FileSetGlobalAttribute( fid, "scale_rm_grid_index_khalo", hinfo%halosize(1:1) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_grid_index_ihalo", hinfo%halosize(2:2) ) ! [IN]
    call FileSetGlobalAttribute( fid, "scale_rm_grid_index_jhalo", hinfo%halosize(3:3) ) ! [IN]

    call FileSetGlobalAttribute( fid, "time_units", hinfo%time_units    ) ! [IN]
    call FileSetGlobalAttribute( fid, "time_start", hinfo%time_start(:) ) ! [IN]

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
    if ( ainfo(2)%periodic == "false" .AND. ainfo(2)%halo_global(1) == 0 ) then
       ainfo(2)%size_global(1) = ainfo(2)%size_global(1) + 1
       ainfo(2)%halo_global(1) = ainfo(2)%halo_global(1) + 1
       if ( rankidx(1) == 0 ) then
          ainfo(2)%start_global(1) = ainfo(2)%start_global(1) + 1
       else
          ainfo(2)%halo_local  (1) = ainfo(2)%halo_local  (1) + 1
       endif
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
    if ( ainfo(4)%periodic == "false" .AND. ainfo(4)%halo_global(1) == 0 ) then
       ainfo(4)%size_global(1) = ainfo(4)%size_global(1) + 1
       ainfo(4)%halo_global(1) = ainfo(4)%halo_global(1) + 1
       if ( rankidx(2) == 0 ) then
          ainfo(4)%start_global(1) = ainfo(4)%start_global(1) + 1
       else
          ainfo(4)%halo_local  (1) = ainfo(4)%halo_local  (1) + 1
       endif
    endif

    minfo%mapping_name                             = hinfo%minfo_mapping_name
    minfo%false_easting                        (:) = hinfo%minfo_false_easting                        (:)
    minfo%false_northing                       (:) = hinfo%minfo_false_northing                       (:)
    minfo%longitude_of_central_meridian        (:) = hinfo%minfo_longitude_of_central_meridian        (:)
    minfo%longitude_of_projection_origin       (:) = hinfo%minfo_longitude_of_projection_origin       (:)
    minfo%latitude_of_projection_origin        (:) = hinfo%minfo_latitude_of_projection_origin        (:)
    minfo%straight_vertical_longitude_from_pole(:) = hinfo%minfo_straight_vertical_longitude_from_pole(:)
    minfo%standard_parallel                    (:) = hinfo%minfo_standard_parallel                    (:)

    call FileSetAttribute( fid, "x" , "size_global" , ainfo(1)%size_global (:) )
    call FileSetAttribute( fid, "x" , "start_global", ainfo(1)%start_global(:) )
    call FileSetAttribute( fid, "x" , "halo_global" , ainfo(1)%halo_global (:) )
    call FileSetAttribute( fid, "x" , "halo_local"  , ainfo(1)%halo_local  (:) )
    call FileSetAttribute( fid, "x" , "periodic"    , ainfo(1)%periodic        )

    call FileSetAttribute( fid, "xh", "size_global" , ainfo(2)%size_global (:) )
    call FileSetAttribute( fid, "xh", "start_global", ainfo(2)%start_global(:) )
    call FileSetAttribute( fid, "xh", "halo_global" , ainfo(2)%halo_global (:) )
    call FileSetAttribute( fid, "xh", "halo_local"  , ainfo(2)%halo_local  (:) )
    call FileSetAttribute( fid, "xh", "periodic"    , ainfo(2)%periodic        )

    call FileSetAttribute( fid, "y" , "size_global" , ainfo(3)%size_global (:) )
    call FileSetAttribute( fid, "y" , "start_global", ainfo(3)%start_global(:) )
    call FileSetAttribute( fid, "y" , "halo_global" , ainfo(3)%halo_global (:) )
    call FileSetAttribute( fid, "y" , "halo_local"  , ainfo(3)%halo_local  (:) )
    call FileSetAttribute( fid, "y" , "periodic"    , ainfo(3)%periodic        )

    call FileSetAttribute( fid, "yh", "size_global" , ainfo(4)%size_global (:) )
    call FileSetAttribute( fid, "yh", "start_global", ainfo(4)%start_global(:) )
    call FileSetAttribute( fid, "yh", "halo_global" , ainfo(4)%halo_global (:) )
    call FileSetAttribute( fid, "yh", "halo_local"  , ainfo(4)%halo_local  (:) )
    call FileSetAttribute( fid, "yh", "periodic"    , ainfo(4)%periodic        )

    if ( minfo%mapping_name /= "" ) then
       call FileSetAttribute( fid, "x" , "standard_name", "projection_x_coordinate" )
       call FileSetAttribute( fid, "xh", "standard_name", "projection_x_coordinate" )
       call FileSetAttribute( fid, "y" , "standard_name", "projection_y_coordinate" )
       call FileSetAttribute( fid, "yh", "standard_name", "projection_y_coordinate" )

       call FileAddAssociatedVariable( fid, minfo%mapping_name )
       call FileSetAttribute( fid, minfo%mapping_name, "grid_mapping_name",  minfo%mapping_name )

       if ( minfo%false_easting(1) /= UNDEF ) then
          call FileSetAttribute( fid,                   & ! [IN]
                                 minfo%mapping_name,    & ! [IN]
                                 "false_easting",       & ! [IN]
                                 minfo%false_easting(:) ) ! [IN]
       endif

       if ( minfo%false_northing(1) /= UNDEF ) then
          call FileSetAttribute( fid,                    & ! [IN]
                                 minfo%mapping_name,     & ! [IN]
                                 "false_northing",       & ! [IN]
                                 minfo%false_northing(:) ) ! [IN]
       endif

       if ( minfo%longitude_of_central_meridian(1) /= UNDEF ) then
          call FileSetAttribute( fid,                                   & ! [IN]
                                 minfo%mapping_name,                    & ! [IN]
                                 "longitude_of_central_meridian",       & ! [IN]
                                 minfo%longitude_of_central_meridian(:) ) ! [IN]
       endif

       if ( minfo%longitude_of_projection_origin(1) /= UNDEF ) then
          call FileSetAttribute( fid,                                    & ! [IN]
                                 minfo%mapping_name,                     & ! [IN]
                                 "longitude_of_projection_origin",       & ! [IN]
                                 minfo%longitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( minfo%latitude_of_projection_origin(1) /= UNDEF ) then
          call FileSetAttribute( fid,                                   & ! [IN]
                                 minfo%mapping_name,                    & ! [IN]
                                 "latitude_of_projection_origin",       & ! [IN]
                                 minfo%latitude_of_projection_origin(:) ) ! [IN]
       endif

       if ( minfo%straight_vertical_longitude_from_pole(1) /= UNDEF ) then
          call FileSetAttribute( fid,                                           & ! [IN]
                                 minfo%mapping_name,                            & ! [IN]
                                 "straight_vertical_longitude_from_pole",       & ! [IN]
                                 minfo%straight_vertical_longitude_from_pole(:) ) ! [IN]
       endif

       if ( minfo%standard_parallel(1) /= UNDEF ) then
          if ( minfo%standard_parallel(2) /= UNDEF ) then
             call FileSetAttribute( fid,                         & ! [IN]
                                    minfo%mapping_name,          & ! [IN]
                                    "standard_parallel",         & ! [IN]
                                    minfo%standard_parallel(1:2) ) ! [IN]
          else
             call FileSetAttribute( fid,                         & ! [IN]
                                    minfo%mapping_name,          & ! [IN]
                                    "standard_parallel",         & ! [IN]
                                    minfo%standard_parallel(1:1) ) ! [IN]
          endif
       endif
    endif

    return
  end subroutine SNO_attributes_write

end module mod_sno
