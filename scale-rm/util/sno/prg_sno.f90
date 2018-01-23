!-------------------------------------------------------------------------------
!> Program SNO (RM)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program sno
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use dc_log, only: &
     LogInit
  use scale_file, only: &
     FILE_Close_All
  use scale_precision
  use scale_stdio
  use scale_prof

  use scale_process, only: &
     PRC_MPIstart,       &
     PRC_MPIstop,        &
     PRC_MPIfinish,      &
     PRC_SINGLECOM_setup
  use scale_const, only: &
     CONST_setup
  use scale_calendar, only: &
     CALENDAR_setup
  use mod_sno_h
  use mod_sno, only: &
     SNO_proc_alloc,   &
     SNO_file_getinfo, &
     SNO_calc_localsize
  use mod_sno_map, only: &
     SNO_map_settable_global, &
     SNO_map_settable_local
  use mod_sno_axis, only: &
     SNO_axis_getinfo, &
     SNO_axis_alloc,   &
     SNO_axis_dealloc, &
     SNO_axis_copy,    &
     SNO_axis_read
  use mod_sno_vars, only: &
     SNO_vars_getinfo, &
     SNO_vars_alloc,   &
     SNO_vars_dealloc, &
     SNO_vars_read,    &
     SNO_vars_write
  use mod_sno_grads, only: &
     SNO_grads_setup,    &
     SNO_grads_netcdfctl
  use mod_snoplugin_timeave, only: &
     SNOPLGIN_timeave_setup,   &
     SNOPLGIN_timeave_alloc,   &
     SNOPLGIN_timeave_dealloc, &
     SNOPLGIN_timeave_store
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  ! namelist parameters
  character(len=H_LONG)   :: basename_in      = ''       ! Basename of the input  file
  character(len=H_LONG)   :: dirpath_out      = ''       ! directory path for the output file
  character(len=H_LONG)   :: basename_out     = ''       ! Basename of the output file
  integer                 :: nprocs_x_out     = 1        ! x length of 2D processor topology (output)
  integer                 :: nprocs_y_out     = 1        ! y length of 2D processor topology (output)
  character(len=H_SHORT)  :: vars(item_limit) = ''       ! name of variables
  logical                 :: output_grads     = .false.  ! output grads fortmat file?
  logical                 :: output_gradsctl  = .false.  ! output grads control file for reading single NetCDF file?
  logical                 :: debug            = .false.

  namelist / PARAM_SNO / &
       basename_in,     &
       dirpath_out,     &
       basename_out,    &
       nprocs_x_out,    &
       nprocs_y_out,    &
       vars,            &
       output_grads,    &
       output_gradsctl, &
       debug

  ! MPI parameters
  integer                 :: comm                        ! communicator                      (execution)
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)

  ! file management [SNO_proc_allloc]
  integer                 :: pstr                        ! start index of peXXXXXX to manage (execution)
  integer                 :: pend                        ! end   index of peXXXXXX to manage (execution)

  ! process & grid information from input file [SNO_file_getinfo]
  integer                 :: nprocs_x_in                 ! x length of 2D processor topology (input)
  integer                 :: nprocs_y_in                 ! y length of 2D processor topology (input)
  integer                 :: ngrids_x                    ! size of x-axis grids              (global,sometimes including halo)
  integer                 :: ngrids_y                    ! size of y-axis grids              (global,sometimes including halo)
  integer                 :: nhalos_x                    ! size of x-axis halo grids         (global,sometimes have a size)
  integer                 :: nhalos_y                    ! size of y-axis halo grids         (global,sometimes have a size)

  type(commoninfo)        :: hinfo
  integer                 :: nvars                       ! number of variables               (input)
  character(len=H_SHORT)  :: varname(item_limit)         ! name   of variables               (input)
  integer                 :: naxis                       ! number of axis variables          (input)
  character(len=H_SHORT)  :: axisname(item_limit)        ! name   of axis variables          (input)

  ! axis information from input file [SNO_axis_getinfo]
  type(axisinfo), allocatable :: ainfo(:)

  ! variable information from input file [SNO_vars_getinfo]
  type(iteminfo), allocatable :: dinfo(:)

  ! mapping table [SNO_calc_localsize,SNO_map_settable_global,SNO_map_settable_local]
  integer                 :: ngrids_x_out                ! size of x-axis grids              (output,sometimes including halo)
  integer                 :: ngrids_y_out                ! size of y-axis grids              (output,sometimes including halo)
  integer                 :: ngrids_xh_out               ! size of x-axis grids, staggard    (output,sometimes including halo)
  integer                 :: ngrids_yh_out               ! size of y-axis grids, staggard    (output,sometimes including halo)
  integer(2), allocatable :: globalmap(:,:,:)            ! mapping table                     (global)
  integer(2), allocatable :: localmap (:,:,:)            ! mapping table                     (output)
  logical,    allocatable :: readflag (:,:)              ! flag to read each input file
  integer                 :: ipos                        ! offset of i-index
  integer                 :: jpos                        ! offset of j-index

  ! Plugins
  logical                 :: do_output
  logical                 :: plugin_timeave

  integer :: px, py, p
  integer :: t, v
  integer :: ierr
  !-----------------------------------------------------------------------------

  ! start MPI
  call PRC_MPIstart( comm ) ! [OUT]

  ! setup MPI communicator
  call PRC_SINGLECOM_setup( comm,    & ! [IN]
                            nprocs,  & ! [OUT]
                            myrank,  & ! [OUT]
                            ismaster ) ! [OUT]

  ! setup standard I/O
  call IO_setup( TOOLNAME )

  ! setup Log
  call IO_LOG_setup( myrank, ismaster )
  call LogInit( IO_FID_CONF,       &
                IO_FID_LOG, IO_L,  &
                IO_FID_NML, IO_NML )

  ! setup profiler
  call PROF_setup

  call PROF_rapstart('Main', 0)
  !########## main ##########

  ! setup constants
  call CONST_setup

  ! setup calendar
  call CALENDAR_setup

  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '++++++ Program[SNO] / Categ[Utility] / Origin[SCALE-RM]'

  !--- read namelist
  rewind(IO_FID_CONF)
  read(IO_FID_CONF,nml=PARAM_SNO,iostat=ierr)
  if ( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
  elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_SNO. Check!'
     call PRC_MPIstop
  endif
  if( IO_NML ) write(IO_FID_NML,nml=PARAM_SNO)

  do_output = .true.

  call SNO_grads_setup( nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                        output_grads,               & ! [IN] from namelist
                        output_gradsctl             ) ! [IN] from namelist

  call SNOPLGIN_timeave_setup( plugin_timeave, & ! [OUT]
                               do_output       ) ! [INOUT]

  ! allocate output files to executing processes
  call SNO_proc_alloc( nprocs, myrank, ismaster,   & ! [IN] from MPI
                       nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                       pstr, pend                  ) ! [OUT]

  !#############################################################################
  ! global setting
  !#############################################################################

  ! get common information from input file
  call SNO_file_getinfo( ismaster,                   & ! [IN] from MPI
                         basename_in,                & ! [IN] from namelist
                         vars(:),                    & ! [IN] from namelist
                         nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                         nprocs_x_in,  nprocs_y_in,  & ! [OUT]
                         ngrids_x,     ngrids_y,     & ! [OUT]
                         nhalos_x,     nhalos_y,     & ! [OUT]
                         hinfo,                      & ! [OUT]
                         naxis,                      & ! [OUT]
                         axisname(:),                & ! [OUT]
                         nvars,                      & ! [OUT]
                         varname(:),                 & ! [OUT]
                         debug                       ) ! [IN]

  ! in->out mapping table (global)
  allocate( globalmap(ngrids_x,ngrids_y,3) )

  call SNO_map_settable_global( nprocs_x_in, nprocs_y_in, & ! [IN] from SNO_file_getinfo
                                ngrids_x,    ngrids_y,    & ! [IN] from SNO_file_getinfo
                                nhalos_x,    nhalos_y,    & ! [IN] from SNO_file_getinfo
                                globalmap(:,:,:),         & ! [OUT]
                                debug                     ) ! [IN]

  ! get information of axis from input file
  allocate( ainfo(naxis) )

  call SNO_axis_getinfo( ismaster,          & ! [IN] from MPI
                         basename_in,       & ! [IN] from namelist
                         naxis,             & ! [IN] from SNO_file_getinfo
                         axisname(1:naxis), & ! [IN] from SNO_file_getinfo
                         ainfo   (:),       & ! [OUT]
                         debug              ) ! [IN]

  ! get information of variables from input file
  allocate( dinfo(nvars) )

  call SNO_vars_getinfo( ismaster,         & ! [IN] from MPI
                         basename_in,      & ! [IN] from namelist
                         nvars,            & ! [IN] from SNO_file_getinfo
                         varname(1:nvars), & ! [IN] from SNO_file_getinfo
                         dinfo  (:),       & ! [OUT]
                         debug             ) ! [IN]

  allocate( readflag(nprocs_x_in,nprocs_y_in) )

  !#############################################################################
  ! process each output file
  !#############################################################################

  jpos = 0
  do py = 1, nprocs_y_out

     ipos = 0
     do px = 1, nprocs_x_out

        p = (py-1) * nprocs_x_out + px - 1

        call SNO_calc_localsize( nprocs_x_out,  nprocs_y_out, & ! [IN] from namelist
                                 px,            py,           & ! [IN]
                                 ngrids_x,      ngrids_y,     & ! [IN] from SNO_file_getinfo
                                 nhalos_x,      nhalos_y,     & ! [IN] from SNO_file_getinfo
                                 hinfo,                       & ! [IN]    from SNO_file_getinfo
                                 ngrids_x_out,  ngrids_y_out, & ! [OUT]
                                 ngrids_xh_out, ngrids_yh_out ) ! [OUT]

        if ( p >= pstr .AND. p <= pend ) then
           if( IO_L ) write(IO_FID_LOG,*)
           if( IO_L ) write(IO_FID_LOG,*) '*** now processing p = ', p

           ! in->out mapping table (for one file)
           allocate( localmap(ngrids_x_out,ngrids_y_out,3) )

           call SNO_map_settable_local( nprocs_x_in,  nprocs_y_in,  & ! [IN] from SNO_file_getinfo
                                        ngrids_x,     ngrids_y,     & ! [IN] from SNO_file_getinfo
                                        ngrids_x_out, ngrids_y_out, & ! [IN] from SNO_map_getsize_local
                                        ipos,         jpos,         & ! [IN] offset
                                        globalmap(:,:,:),           & ! [IN] from SNO_map_settable_global
                                        localmap (:,:,:),           & ! [OUT]
                                        readflag (:,:),             & ! [OUT]
                                        debug                       ) ! [IN]

           ! read axis and rearrange
           call SNO_axis_alloc( nprocs_x_out,  nprocs_y_out,  & ! [IN]    from namelist
                                ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                hinfo,                        & ! [IN]    from SNO_file_getinfo
                                naxis,                        & ! [IN]    from SNO_file_getinfo
                                ainfo(:),                     & ! [INOUT] from SNO_axis_getinfo
                                debug                         ) ! [IN]

           call SNO_axis_copy( nprocs_x_out,  nprocs_y_out,  & ! [IN]    from namelist
                               px,            py,            & ! [IN]
                               hinfo,                        & ! [IN]    from SNO_file_getinfo
                               naxis,                        & ! [IN]    from SNO_file_getinfo
                               ainfo(:),                     & ! [INOUT] from SNO_axis_getinfo
                               debug                         ) ! [IN]

           call SNO_axis_read( basename_in,                  & ! [IN]    from namelist
                               nprocs_x_in,   nprocs_y_in,   & ! [IN]    from SNO_file_getinfo
                               ngrids_x,      ngrids_y,      & ! [IN]    from SNO_file_getinfo
                               nhalos_x,      nhalos_y,      & ! [IN]    from SNO_file_getinfo
                               hinfo,                        & ! [IN]    from SNO_file_getinfo
                               ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                               ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                               naxis,                        & ! [IN]    from SNO_file_getinfo
                               ainfo   (:),                  & ! [INOUT] from SNO_axis_getinfo
                               localmap(:,:,:),              & ! [IN]    from SNO_map_settable_local
                               readflag(:,:),                & ! [IN]    from SNO_map_settable_local
                               debug                         ) ! [IN]

           !####################################################################
           ! process each variable
           !####################################################################

           do v = 1, nvars
              if( IO_L ) write(IO_FID_LOG,*)
              if( IO_L ) write(IO_FID_LOG,*) '*** + now processing varname = ', trim(dinfo(v)%varname)

              ! output array allocation

              call SNO_vars_alloc( ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                   ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                   dinfo(v),                     & ! [INOUT] from SNO_vars_getinfo
                                   debug                         ) ! [IN]

              if( plugin_timeave ) call SNOPLGIN_timeave_alloc( dinfo(v), & ! [INOUT] from SNO_vars_getinfo
                                                                debug     ) ! [IN]

              !#################################################################
              ! process each timestep
              !#################################################################

              do t = 1, dinfo(v)%step_nmax
                 if( IO_L ) write(IO_FID_LOG,*)
                 if( IO_L ) write(IO_FID_LOG,*) '*** + + now processing t = ', t

                 call SNO_vars_read( basename_in,                  & ! [IN]    from namelist
                                     t,                            & ! [IN]
                                     nprocs_x_in,   nprocs_y_in,   & ! [IN]    from SNO_file_getinfo
                                     ngrids_x,      ngrids_y,      & ! [IN]    from SNO_file_getinfo
                                     nhalos_x,      nhalos_y,      & ! [IN]    from SNO_file_getinfo
                                     hinfo,                        & ! [IN]    from SNO_file_getinfo
                                     ngrids_x_out,  ngrids_y_out,  & ! [IN]    from SNO_map_getsize_local
                                     ngrids_xh_out, ngrids_yh_out, & ! [IN]    from SNO_map_getsize_local
                                     dinfo(v),                     & ! [INOUT] from SNO_vars_getinfo
                                     localmap(:,:,:),              & ! [IN]    from SNO_map_settable_local
                                     readflag(:,:),                & ! [IN]    from SNO_map_settable_local
                                     debug                         ) ! [IN]

                 if( plugin_timeave ) call SNOPLGIN_timeave_store( dirpath_out,                & ! [IN] from namelist
                                                                   basename_out,               & ! [IN] from namelist
                                                                   output_grads,               & ! [IN] from namelist
                                                                   p,                          & ! [IN]
                                                                   t,                          & ! [IN]
                                                                   nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                                                                   nhalos_x,     nhalos_y,     & ! [IN] from SNO_file_getinfo
                                                                   hinfo,                      & ! [IN] from SNO_file_getinfo
                                                                   naxis,                      & ! [IN] from SNO_file_getinfo
                                                                   ainfo(:),                   & ! [IN] from SNO_axis_getinfo
                                                                   dinfo(v),                   & ! [IN] from SNO_vars_getinfo
                                                                   debug                       ) ! [IN]
                 if ( do_output ) then
                    call SNO_vars_write( dirpath_out,                & ! [IN] from namelist
                                         basename_out,               & ! [IN] from namelist
                                         output_grads,               & ! [IN] from namelist
                                         p,                          & ! [IN]
                                         t,                          & ! [IN]
                                         nprocs_x_out, nprocs_y_out, & ! [IN] from namelist
                                         nhalos_x,     nhalos_y,     & ! [IN] from SNO_file_getinfo
                                         hinfo,                      & ! [IN] from SNO_file_getinfo
                                         naxis,                      & ! [IN] from SNO_file_getinfo
                                         ainfo(:),                   & ! [IN] from SNO_axis_getinfo
                                         dinfo(v),                   & ! [IN] from SNO_vars_getinfo
                                         debug                       ) ! [IN]
                 endif
              enddo ! t loop

              ! output array deallocation

              call SNO_vars_dealloc( dinfo(v), & ! [INOUT] from SNO_vars_getinfo
                                     debug     ) ! [IN]

              if( plugin_timeave ) call SNOPLGIN_timeave_dealloc( debug ) ! [IN]

           enddo ! item loop

           if ( do_output ) then
              if ( output_gradsctl ) then
                 call SNO_grads_netcdfctl( dirpath_out,  & ! [IN] from namelist
                                           basename_out, & ! [IN] from namelist
                                           hinfo,        & ! [IN] from SNO_file_getinfo
                                           naxis,        & ! [IN] from SNO_file_getinfo
                                           ainfo(:),     & ! [IN] from SNO_axis_getinfo
                                           nvars,        & ! [IN] from SNO_file_getinfo
                                           dinfo(:),     & ! [IN] from SNO_vars_getinfo
                                           debug         ) ! [IN]
              endif
           endif

           call SNO_axis_dealloc( naxis,    & ! [IN]    from SNO_file_getinfo
                                  ainfo(:), & ! [INOUT] from SNO_axis_getinfo
                                  debug     ) ! [IN]

           deallocate( localmap )

        endif ! managed by this execution rank?

        ipos = ipos + ngrids_x_out
     enddo ! px loop

     jpos = jpos + ngrids_y_out
  enddo ! py loop

  !########## Finalize ##########
  call PROF_rapend  ('Main', 0)

  call FILE_Close_All

  call PROF_rapreport

  ! stop MPI
  call PRC_MPIfinish

end program sno
