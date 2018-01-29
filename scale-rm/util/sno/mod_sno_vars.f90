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
module mod_sno_vars
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
  public :: SNO_vars_getinfo
  public :: SNO_vars_alloc
  public :: SNO_vars_dealloc
  public :: SNO_vars_read
  public :: SNO_vars_write
  public :: SNO_vars_write_netcdf

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
  subroutine SNO_vars_getinfo( &
       ismaster, &
       basename, &
       nvars,    &
       varname,  &
       dinfo,    &
       debug     )
    use mpi
    use scale_file_h, only: &
       FILE_dtypelist
    use scale_file, only: &
       FILE_Get_All_Datainfo
    use scale_process, only: &
       PRC_masterrank,       &
       PRC_LOCAL_COMM_WORLD, &
       PRC_MPIstop
    use scale_calendar, only: &
       CALENDAR_daysec2date,   &
       CALENDAR_adjust_daysec, &
       CALENDAR_CFunits2sec,   &
       CALENDAR_date2char
    use mod_sno_h, only: &
       step_limit, &
       dim_limit,  &
       iteminfo
    implicit none

    logical,          intent(in)  :: ismaster                 ! master process?                    (execution)
    character(len=*), intent(in)  :: basename                 ! basename of file                   (input)
    integer,          intent(in)  :: nvars                    ! number of variables
    character(len=*), intent(in)  :: varname(nvars)           ! name of variables                  (input)
    type(iteminfo),   intent(out) :: dinfo  (nvars)           ! variable information               (input)
    logical,          intent(in)  :: debug

    integer           :: start_day    !< absolute day    (time=0)
    real(DP)          :: start_sec    !< absolute second (time=0)
    integer           :: now_day      !< absolute day    (time=t)
    real(DP)          :: now_sec      !< absolute second (time=t)
    integer           :: now_date(6)  !< date            (time=t)
    real(DP)          :: now_ms       !< subsecond       (time=t)
    character(len=27) :: now_chardate !< date            (time=t)

    integer :: nowrank
    integer :: ierr
    integer :: v, d, t
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Read information of variables'

    if ( ismaster ) then
       nowrank = 0 ! first file

       do v = 1, nvars
          call FILE_Get_All_Datainfo( step_limit  = step_limit,             & ! [IN]
                                      dim_limit   = dim_limit,              & ! [IN]
                                      basename    = basename,               & ! [IN]
                                      varname     = varname(v),             & ! [IN]
                                      rankid      = nowrank,                & ! [IN]
                                      step_nmax   = dinfo(v)%step_nmax,     & ! [OUT]
                                      description = dinfo(v)%description,   & ! [OUT]
                                      units       = dinfo(v)%units,         & ! [OUT]
                                      datatype    = dinfo(v)%datatype,      & ! [OUT]
                                      dim_rank    = dinfo(v)%dim_rank,      & ! [OUT]
                                      dim_name    = dinfo(v)%dim_name  (:), & ! [OUT]
                                      dim_size    = dinfo(v)%dim_size  (:), & ! [OUT]
                                      time_start  = dinfo(v)%time_start(:), & ! [OUT]
                                      time_end    = dinfo(v)%time_end  (:), & ! [OUT]
                                      time_units  = dinfo(v)%time_units     ) ! [OUT]

          dinfo(v)%varname = varname(v)

          dinfo(v)%transpose = .false.
          if ( dinfo(v)%dim_rank > 2 ) then
             if (       dinfo(v)%dim_name(1) /= 'z'   &
                  .AND. dinfo(v)%dim_name(1) /= 'zh'  &
                  .AND. dinfo(v)%dim_name(1) /= 'oz'  &
                  .AND. dinfo(v)%dim_name(1) /= 'ozh' &
                  .AND. dinfo(v)%dim_name(1) /= 'lz'  &
                  .AND. dinfo(v)%dim_name(1) /= 'lzh' &
                  .AND. dinfo(v)%dim_name(1) /= 'uz'  &
                  .AND. dinfo(v)%dim_name(1) /= 'uzh' ) then
                dinfo(v)%transpose = .true.
             endif
          endif

          if ( dinfo(v)%step_nmax > 1 ) then
             dinfo(v)%dt = dinfo(v)%time_start(2) - dinfo(v)%time_start(1)
          else
             dinfo(v)%dt = 0.0_DP
          endif

       enddo
    endif

    do v = 1, nvars
       call MPI_BCAST( dinfo(v)%varname      , H_SHORT          , MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%description  , H_MID            , MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%units        , H_SHORT          , MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%datatype     , 1                , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%dim_rank     , 1                , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%dim_name(:)  , H_SHORT*dim_limit, MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%dim_size(:)  , dim_limit        , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%transpose    , 1                , MPI_LOGICAL         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%step_nmax    , 1                , MPI_INTEGER         , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%time_start(:), step_limit       , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%time_end  (:), step_limit       , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%dt           , 1                , MPI_DOUBLE_PRECISION, PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )
       call MPI_BCAST( dinfo(v)%time_units   , H_MID            , MPI_CHARACTER       , PRC_masterrank, PRC_LOCAL_COMM_WORLD, ierr )

       if ( debug ) then
          if( IO_L ) write(IO_FID_LOG,*)
          if( IO_L ) write(IO_FID_LOG,*) '*** Var No.', v
          if( IO_L ) write(IO_FID_LOG,*) '*** varname     : ', trim(dinfo(v)%varname)
          if( IO_L ) write(IO_FID_LOG,*) '*** description : ', trim(dinfo(v)%description)
          if( IO_L ) write(IO_FID_LOG,*) '*** units       : ', trim(dinfo(v)%units)
          if( IO_L ) write(IO_FID_LOG,*) '*** datatype    : ', trim(FILE_dtypelist(dinfo(v)%datatype))
          if( IO_L ) write(IO_FID_LOG,*) '*** dim_rank    : ', dinfo(v)%dim_rank
          do d = 1, dinfo(v)%dim_rank
             if( IO_L ) write(IO_FID_LOG,*) '*** dim No.', d
             if( IO_L ) write(IO_FID_LOG,*) '*** + dim_name  : ', trim(dinfo(v)%dim_name(d))
             if( IO_L ) write(IO_FID_LOG,*) '*** + dim_size  : ', dinfo(v)%dim_size(d)
          enddo
          if( IO_L ) write(IO_FID_LOG,*) '*** transpose   : ', dinfo(v)%transpose

          if( IO_L ) write(IO_FID_LOG,*) '*** time_units  : ', trim(dinfo(v)%time_units)
          if( IO_L ) write(IO_FID_LOG,*) '*** step_nmax   : ', dinfo(v)%step_nmax

          start_day = 0
          start_sec = CALENDAR_CFunits2sec( cftime=0.0_DP, cfunits=dinfo(v)%time_units, offset_year=0 )
          call CALENDAR_adjust_daysec( start_day, start_sec )

          do t = 1, dinfo(v)%step_nmax
             now_day = start_day
             now_sec = start_sec + 0.5_DP * ( dinfo(v)%time_start(t) + dinfo(v)%time_end(t) )

             call CALENDAR_adjust_daysec( now_day, now_sec ) ! [INOUT]

             call CALENDAR_daysec2date  ( now_date, now_ms,                & ! [OUT]
                                          now_day,  now_sec, offset_year=0 ) ! [IN]

             call CALENDAR_date2char    ( now_chardate,    & ! [OUT]
                                          now_date, now_ms ) ! [IN]

             if( IO_L ) write(IO_FID_LOG,*) '*** + ', now_chardate
          enddo
       endif
    enddo

    return
  end subroutine SNO_vars_getinfo

  !-----------------------------------------------------------------------------
  subroutine SNO_vars_alloc( &
       ngrids_x_out,  &
       ngrids_y_out,  &
       ngrids_xh_out, &
       ngrids_yh_out, &
       dinfo,         &
       debug          )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    integer,        intent(in)    :: ngrids_x_out             ! number of x-axis grids per process (output,sometimes including halo)
    integer,        intent(in)    :: ngrids_y_out             ! number of y-axis grids per process (output,sometimes including halo)
    integer,        intent(in)    :: ngrids_xh_out            ! number of x-axis grids per process (output,sometimes including halo,staggered)
    integer,        intent(in)    :: ngrids_yh_out            ! number of y-axis grids per process (output,sometimes including halo,staggered)
    type(iteminfo), intent(inout) :: dinfo                    ! variable information               (input)
    logical,        intent(in)    :: debug

    integer  :: gout1, gout2, gout3
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Allocate variable array'

    if ( dinfo%dim_rank == 1 ) then

       gout1 = dinfo%dim_size(1)

       allocate( dinfo%VAR_1d(gout1) )
       dinfo%VAR_1d(:) = 0.0_RP

    elseif( dinfo%dim_rank == 2 ) then

       if    ( dinfo%dim_name(1) == 'x'  ) then
          gout1 = ngrids_x_out
       elseif( dinfo%dim_name(1) == 'xh' ) then
          gout1 = ngrids_xh_out
       else
          gout1 = dinfo%dim_size(1)
       endif

       if    ( dinfo%dim_name(2) == 'y'  ) then
          gout2 = ngrids_y_out
       elseif( dinfo%dim_name(2) == 'yh' ) then
          gout2 = ngrids_yh_out
       else
          gout2 = dinfo%dim_size(2)
       endif

       allocate( dinfo%VAR_2d(gout1,gout2) )
       dinfo%VAR_2d(:,:) = 0.0_RP

    elseif( dinfo%dim_rank == 3 ) then

       if ( dinfo%transpose ) then
          gout1 = dinfo%dim_size(3)

          if    ( dinfo%dim_name(1) == 'x'  ) then
             gout2 = ngrids_x_out
          elseif( dinfo%dim_name(1) == 'xh' ) then
             gout2 = ngrids_xh_out
          else
             gout2 = dinfo%dim_size(1)
          endif

          if    ( dinfo%dim_name(2) == 'y'  ) then
             gout3 = ngrids_y_out
          elseif( dinfo%dim_name(2) == 'yh' ) then
             gout3 = ngrids_yh_out
          else
             gout3 = dinfo%dim_size(2)
          endif
       else
          gout1 = dinfo%dim_size(1)

          if    ( dinfo%dim_name(2) == 'x'  ) then
             gout2 = ngrids_x_out
          elseif( dinfo%dim_name(2) == 'xh' ) then
             gout2 = ngrids_xh_out
          else
             gout2 = dinfo%dim_size(2)
          endif

          if    ( dinfo%dim_name(3) == 'y'  ) then
             gout3 = ngrids_y_out
          elseif( dinfo%dim_name(3) == 'yh' ) then
             gout3 = ngrids_yh_out
          else
             gout3 = dinfo%dim_size(3)
          endif
       endif

       allocate( dinfo%VAR_3d(gout1,gout2,gout3) )
       dinfo%VAR_3d(:,:,:) = 0.0_RP

    endif

    return
  end subroutine SNO_vars_alloc

  !-----------------------------------------------------------------------------
  subroutine SNO_vars_dealloc( &
       dinfo, &
       debug  )
    use mod_sno_h, only: &
       iteminfo
    implicit none

    type(iteminfo), intent(inout) :: dinfo ! variable information               (input)
    logical,        intent(in)    :: debug
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Deallocate variable array'

    if( allocated(dinfo%VAR_1d) ) deallocate( dinfo%VAR_1d )
    if( allocated(dinfo%VAR_2d) ) deallocate( dinfo%VAR_2d )
    if( allocated(dinfo%VAR_3d) ) deallocate( dinfo%VAR_3d )

    return
  end subroutine SNO_vars_dealloc

  !-----------------------------------------------------------------------------
  subroutine SNO_vars_read( &
       basename,      &
       nowstep,       &
       nprocs_x_in,   &
       nprocs_y_in,   &
       ngrids_x,      &
       ngrids_y,      &
       nhalos_x,      &
       nhalos_y,      &
       hinfo,         &
       ngrids_x_out,  &
       ngrids_y_out,  &
       ngrids_xh_out, &
       ngrids_yh_out, &
       dinfo,         &
       localmap,      &
       readflag,      &
       debug          )
    use mod_sno_h, only: &
       I_map_p,    &
       I_map_i,    &
       I_map_j,    &
       commoninfo, &
       iteminfo
    use mod_sno, only: &
       SNO_calc_localsize, &
       SNO_read_map_1d,    &
       SNO_read_map_2d,    &
       SNO_read_map_3d
    implicit none

    character(len=*), intent(in)    :: basename                              ! basename of file                   (input)
    integer,          intent(in)    :: nowstep                               ! current step                       (input)
    integer,          intent(in)    :: nprocs_x_in                           ! x length of 2D processor topology  (input)
    integer,          intent(in)    :: nprocs_y_in                           ! y length of 2D processor topology  (input)
    integer,          intent(in)    :: ngrids_x                              ! number of x-axis grids             (global domain,sometimes including halo)
    integer,          intent(in)    :: ngrids_y                              ! number of y-axis grids             (global domain,sometimes including halo)
    integer,          intent(in)    :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)    :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    integer,          intent(in)    :: ngrids_x_out                          ! number of x-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_y_out                          ! number of y-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_xh_out                         ! number of x-axis grids per process (output,sometimes including halo)
    integer,          intent(in)    :: ngrids_yh_out                         ! number of y-axis grids per process (output,sometimes including halo)
    type(iteminfo),   intent(inout) :: dinfo                                 ! variable information               (input)
    integer(2),       intent(in)    :: localmap(ngrids_x_out,ngrids_y_out,3) ! mapping table from input to output (local domain)
    logical,          intent(in)    :: readflag(nprocs_x_in,nprocs_y_in)     ! local domain requires the input file?
    logical,          intent(in)    :: debug

    integer(2), allocatable :: localmap_1d(:,:)
    integer(2), allocatable :: localmap_2d(:,:,:)

    integer  :: ngrids_x_in       ! number of x-axis grids per process (input,without halo)
    integer  :: ngrids_y_in       ! number of y-axis grids per process (input,without halo)
    integer  :: ngrids_xh_in      ! number of x-axis grids per process (input,without halo)
    integer  :: ngrids_yh_in      ! number of y-axis grids per process (input,without halo)
    integer  :: staggered_x_in
    integer  :: staggered_y_in
    integer  :: staggered_x_out
    integer  :: staggered_y_out

    integer  :: gin1, gin2, gin3
    integer  :: gout1, gout2, gout3
    integer  :: stgin1, stgin2, stgin3
    integer  :: stgout1, stgout2, stgout3
    logical  :: readflag_1d

    integer  :: p, px, py
    integer  :: i, j
    !---------------------------------------------------------------------------

    do py = 1, nprocs_y_in
       do px = 1, nprocs_x_in
          p = (py-1) * nprocs_x_in + px - 1

          call SNO_calc_localsize( nprocs_x_in,  nprocs_y_in, & ! [IN] from namelist
                                   px,           py,          & ! [IN]
                                   ngrids_x,     ngrids_y,    & ! [IN] from SNO_file_getinfo
                                   nhalos_x,     nhalos_y,    & ! [IN] from SNO_file_getinfo
                                   hinfo,                     & ! [IN] from SNO_file_getinfo
                                   ngrids_x_in,  ngrids_y_in, & ! [OUT]
                                   ngrids_xh_in, ngrids_yh_in ) ! [OUT]

          staggered_x_in  = 0
          staggered_y_in  = 0
          staggered_x_out = 0
          staggered_y_out = 0
          if ( ngrids_x_in  /= ngrids_xh_in  ) staggered_x_in  = 1
          if ( ngrids_y_in  /= ngrids_yh_in  ) staggered_y_in  = 1
          if ( ngrids_x_out /= ngrids_xh_out ) staggered_x_out = 1
          if ( ngrids_y_out /= ngrids_yh_out ) staggered_y_out = 1

          if ( readflag(px,py) ) then

             if ( dinfo%dim_rank == 1 ) then

                readflag_1d = .false.

                if    ( dinfo%dim_name(1) == 'x'   ) then
                   gin1    = ngrids_x_in
                   gout1   = ngrids_x_out
                   stgin1  = 0
                   stgout1 = 0

                   allocate( localmap_1d(gout1,2) )

                   do i = 1, ngrids_x_out
                      localmap_1d(i,I_map_p) = localmap(i,1,I_map_p)
                      localmap_1d(i,I_map_i) = localmap(i,1,I_map_i)
                   enddo

                   if( py == 1 ) readflag_1d = .true.

                elseif( dinfo%dim_name(1) == 'xh'  ) then
                   gin1    = ngrids_xh_in
                   gout1   = ngrids_xh_out
                   stgin1  = staggered_x_in
                   stgout1 = staggered_x_out

                   allocate( localmap_1d(gout1,2) )

                   do i = 1, ngrids_x_out
                      localmap_1d(i+stgout1,I_map_p) = localmap(i,1,I_map_p)
                      localmap_1d(i+stgout1,I_map_i) = localmap(i,1,I_map_i) + stgin1
                   enddo

                   if ( stgin1 > 0 .AND. stgout1 > 0 ) then
                      do i = 1, stgout1
                         localmap_1d(i,I_map_p) = localmap(1,1,I_map_p)
                         localmap_1d(i,I_map_i) = i
                      enddo
                   endif

                   if( py == 1 ) readflag_1d = .true.

                elseif( dinfo%dim_name(1) == 'y'   ) then
                   gin1    = ngrids_y_in
                   gout1   = ngrids_y_out
                   stgin1  = 0
                   stgout1 = 0

                   allocate( localmap_1d(gout1,2) )

                   do j = 1, ngrids_y_out
                      localmap_1d(j,I_map_p) = localmap(1,j,I_map_p)
                      localmap_1d(j,I_map_i) = localmap(1,j,I_map_j)
                   enddo

                   if( px == 1 ) readflag_1d = .true.

                elseif( dinfo%dim_name(1) == 'yh'  ) then
                   gin1    = ngrids_yh_in
                   gout1   = ngrids_yh_out
                   stgin1  = staggered_y_in
                   stgout1 = staggered_y_out

                   allocate( localmap_1d(gout1,2) )

                   do j = 1, ngrids_y_out
                      localmap_1d(j,I_map_p) = localmap(1,j,I_map_p)
                      localmap_1d(j,I_map_i) = localmap(1,j,I_map_j) + stgin1
                   enddo

                   if ( stgin1 > 0 .AND. stgout1 > 0 ) then
                      do i = 1, stgout1
                         localmap_1d(i,I_map_p) = localmap(1,1,I_map_p)
                         localmap_1d(i,I_map_i) = i
                      enddo
                   endif

                   if( px == 1 ) readflag_1d = .true.

                endif

                if ( readflag_1d ) then
                   call SNO_read_map_1d( basename, p, nowstep, & ! [IN]
                                         dinfo%varname,        & ! [IN]
                                         dinfo%datatype,       & ! [IN]
                                         gin1,                 & ! [IN]
                                         gout1,                & ! [IN]
                                         localmap_1d (:,:),    & ! [IN]
                                         dinfo%VAR_1d(:)       ) ! [INOUT]
                endif

                deallocate( localmap_1d )

             elseif( dinfo%dim_rank == 2 ) then

                if    ( dinfo%dim_name(1) == 'x'  ) then
                   gin1    = ngrids_x_in
                   gout1   = ngrids_x_out
                   stgin1  = 0
                   stgout1 = 0
                elseif( dinfo%dim_name(1) == 'xh' ) then
                   gin1    = ngrids_xh_in
                   gout1   = ngrids_xh_out
                   stgin1  = staggered_x_in
                   stgout1 = staggered_x_out
                else
                   gin1    = dinfo%dim_size(1)
                   gout1   = dinfo%dim_size(1)
                   stgin1  = 0
                   stgout1 = 0
                endif

                if    ( dinfo%dim_name(2) == 'y'  ) then
                   gin2    = ngrids_y_in
                   gout2   = ngrids_y_out
                   stgin2  = 0
                   stgout2 = 0
                elseif( dinfo%dim_name(2) == 'yh' ) then
                   gin2    = ngrids_yh_in
                   gout2   = ngrids_yh_out
                   stgin2  = staggered_y_in
                   stgout2 = staggered_y_out
                else
                   gin2    = dinfo%dim_size(2)
                   gout2   = dinfo%dim_size(2)
                   stgin2  = 0
                   stgout2 = 0
                endif

                allocate( localmap_2d(gout1,gout2,3) )

                do j = 1, ngrids_y_out
                do i = 1, ngrids_x_out
                   localmap_2d(i+stgout1,j+stgout2,I_map_p) = localmap(i,j,I_map_p)
                   localmap_2d(i+stgout1,j+stgout2,I_map_i) = localmap(i,j,I_map_i) + stgin1
                   localmap_2d(i+stgout1,j+stgout2,I_map_j) = localmap(i,j,I_map_j) + stgin2
                enddo
                enddo

                if ( stgout1 > 0 ) then
                   if ( stgin1 > 0 ) then
                      do j = 1, ngrids_y_out
                      do i = 1, stgout1
                         localmap_2d(i,j+stgout2,I_map_p) = localmap(1,j,I_map_p)
                         localmap_2d(i,j+stgout2,I_map_i) = i
                         localmap_2d(i,j+stgout2,I_map_j) = localmap(1,j,I_map_j) + stgin2
                      enddo
                      enddo
                   else
                      do j = 1, ngrids_y_out
                      do i = 1, stgout1
                         localmap_2d(i,j+stgout2,I_map_p) = -1
                         localmap_2d(i,j+stgout2,I_map_i) = -1
                         localmap_2d(i,j+stgout2,I_map_j) = -1
                      enddo
                      enddo
                   endif
                endif

                if ( stgout2 > 0 ) then
                   if ( stgin2 > 0 ) then
                      do j = 1, stgout2
                      do i = 1, ngrids_x_out
                         localmap_2d(i+stgout1,j,I_map_p) = localmap(i,1,I_map_p)
                         localmap_2d(i+stgout1,j,I_map_i) = localmap(i,1,I_map_i) + stgin1
                         localmap_2d(i+stgout1,j,I_map_j) = j
                      enddo
                      enddo
                   else
                      do j = 1, stgout2
                      do i = 1, ngrids_x_out
                         localmap_2d(i+stgout1,j,I_map_p) = -1
                         localmap_2d(i+stgout1,j,I_map_i) = -1
                         localmap_2d(i+stgout1,j,I_map_j) = -1
                      enddo
                      enddo
                   endif
                endif

                if ( stgout1 > 0 .AND. stgout2 > 0 ) then
                   if ( stgin1 > 0 .AND. stgin2 > 0 ) then
                      do j = 1, stgout2
                      do i = 1, stgout1
                         localmap_2d(i,j,I_map_p) = localmap(1,1,I_map_p)
                         localmap_2d(i,j,I_map_i) = i
                         localmap_2d(i,j,I_map_j) = j
                      enddo
                      enddo
                   else
                      do j = 1, stgout2
                      do i = 1, stgout1
                         localmap_2d(i,j,I_map_p) = -1
                         localmap_2d(i,j,I_map_i) = -1
                         localmap_2d(i,j,I_map_j) = -1
                      enddo
                      enddo
                   endif
                endif

                call SNO_read_map_2d( basename, p, nowstep, & ! [IN]
                                      dinfo%varname,        & ! [IN]
                                      dinfo%datatype,       & ! [IN]
                                      gin1,  gin2,          & ! [IN]
                                      gout1, gout2,         & ! [IN]
                                      localmap_2d (:,:,:),  & ! [IN]
                                      dinfo%VAR_2d(:,:)     ) ! [INOUT]

                deallocate( localmap_2d )

             elseif( dinfo%dim_rank == 3 ) then

                if ( dinfo%transpose ) then
                   gin1  = dinfo%dim_size(3)
                   gout1 = dinfo%dim_size(3)

                   if    ( dinfo%dim_name(1) == 'x'  ) then
                      gin2    = ngrids_x_in
                      gout2   = ngrids_x_out
                      stgin2  = 0
                      stgout2 = 0
                   elseif( dinfo%dim_name(1) == 'xh' ) then
                      gin2    = ngrids_xh_in
                      gout2   = ngrids_xh_out
                      stgin2  = staggered_x_in
                      stgout2 = staggered_x_out
                   else
                      gin2    = dinfo%dim_size(1)
                      gout2   = dinfo%dim_size(1)
                      stgin2  = 0
                      stgout2 = 0
                   endif

                   if    ( dinfo%dim_name(2) == 'y'  ) then
                      gin3    = ngrids_y_in
                      gout3   = ngrids_y_out
                      stgin3  = 0
                      stgout3 = 0
                   elseif( dinfo%dim_name(2) == 'yh' ) then
                      gin3    = ngrids_yh_in
                      gout3   = ngrids_yh_out
                      stgin3  = staggered_y_in
                      stgout3 = staggered_y_out
                   else
                      gin3    = dinfo%dim_size(2)
                      gout3   = dinfo%dim_size(2)
                      stgin3  = 0
                      stgout3 = 0
                   endif
                else
                   gin1  = dinfo%dim_size(1)
                   gout1 = dinfo%dim_size(1)

                   if    ( dinfo%dim_name(2) == 'x'  ) then
                      gin2    = ngrids_x_in
                      gout2   = ngrids_x_out
                      stgin2  = 0
                      stgout2 = 0
                   elseif( dinfo%dim_name(2) == 'xh' ) then
                      gin2    = ngrids_xh_in
                      gout2   = ngrids_xh_out
                      stgin2  = staggered_x_in
                      stgout2 = staggered_x_out
                   else
                      gin2    = dinfo%dim_size(2)
                      gout2   = dinfo%dim_size(2)
                      stgin2  = 0
                      stgout2 = 0
                   endif

                   if    ( dinfo%dim_name(3) == 'y'  ) then
                      gin3    = ngrids_y_in
                      gout3   = ngrids_y_out
                      stgin3  = 0
                      stgout3 = 0
                   elseif( dinfo%dim_name(3) == 'yh' ) then
                      gin3    = ngrids_yh_in
                      gout3   = ngrids_yh_out
                      stgin3  = staggered_y_in
                      stgout3 = staggered_y_out
                   else
                      gin3    = dinfo%dim_size(3)
                      gout3   = dinfo%dim_size(3)
                      stgin3  = 0
                      stgout3 = 0
                   endif
                endif

                allocate( localmap_2d(gout2,gout3,3) )

                do j = 1, ngrids_y_out
                do i = 1, ngrids_x_out
                   localmap_2d(i+stgout2,j+stgout3,I_map_p) = localmap(i,j,I_map_p)
                   localmap_2d(i+stgout2,j+stgout3,I_map_i) = localmap(i,j,I_map_i) + stgin2
                   localmap_2d(i+stgout2,j+stgout3,I_map_j) = localmap(i,j,I_map_j) + stgin3
                enddo
                enddo

                if ( stgout2 > 0 ) then
                   if ( stgin2 > 0 ) then
                      do j = 1, ngrids_y_out
                      do i = 1, stgout2
                         localmap_2d(i,j+stgout3,I_map_p) = localmap(1,j,I_map_p)
                         localmap_2d(i,j+stgout3,I_map_i) = i
                         localmap_2d(i,j+stgout3,I_map_j) = localmap(1,j,I_map_j) + stgin3
                      enddo
                      enddo
                   else
                      do j = 1, ngrids_y_out
                      do i = 1, stgout2
                         localmap_2d(i,j+stgout3,I_map_p) = -1
                         localmap_2d(i,j+stgout3,I_map_i) = -1
                         localmap_2d(i,j+stgout3,I_map_j) = -1
                      enddo
                      enddo
                   endif
                endif

                if ( stgout3 > 0 ) then
                   if ( stgin3 > 0 ) then
                      do j = 1, stgout3
                      do i = 1, ngrids_x_out
                         localmap_2d(i+stgout2,j,I_map_p) = localmap(i,1,I_map_p)
                         localmap_2d(i+stgout2,j,I_map_i) = localmap(i,1,I_map_i) + stgin2
                         localmap_2d(i+stgout2,j,I_map_j) = j
                      enddo
                      enddo
                   else
                      do j = 1, stgout3
                      do i = 1, ngrids_x_out
                         localmap_2d(i+stgout2,j,I_map_p) = -1
                         localmap_2d(i+stgout2,j,I_map_i) = -1
                         localmap_2d(i+stgout2,j,I_map_j) = -1
                      enddo
                      enddo
                   endif
                endif

                if ( stgout2 > 0 .AND. stgout3 > 0 ) then
                   if ( stgin2 > 0 .AND. stgin3 > 0 ) then
                      do j = 1, stgout3
                      do i = 1, stgout2
                         localmap_2d(i,j,I_map_p) = localmap(1,1,I_map_p)
                         localmap_2d(i,j,I_map_i) = i
                         localmap_2d(i,j,I_map_j) = j
                      enddo
                      enddo
                   else
                      do j = 1, stgout3
                      do i = 1, stgout2
                         localmap_2d(i,j,I_map_p) = -1
                         localmap_2d(i,j,I_map_i) = -1
                         localmap_2d(i,j,I_map_j) = -1
                      enddo
                      enddo
                   endif
                endif

                call SNO_read_map_3d( basename, p, nowstep, & ! [IN]
                                      dinfo%varname,        & ! [IN]
                                      dinfo%datatype,       & ! [IN]
                                      dinfo%transpose,      & ! [IN]
                                      gin1,  gin2,  gin3,   & ! [IN]
                                      gout1, gout2, gout3,  & ! [IN]
                                      localmap_2d (:,:,:),  & ! [IN]
                                      dinfo%VAR_3d(:,:,:)   ) ! [INOUT]

                deallocate( localmap_2d )

             endif ! dim rank?

          endif ! readflag?
       enddo ! px
    enddo ! py

    if ( debug ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Varname : ', trim(dinfo%varname)

       if ( allocated(dinfo%VAR_1d) ) then
          if( IO_L ) write(IO_FID_LOG,*) dinfo%VAR_1d(:)
       endif
       if ( allocated(dinfo%VAR_2d) ) then
          if( IO_L ) write(IO_FID_LOG,*) dinfo%VAR_2d(:,:)
       endif
       if ( allocated(dinfo%VAR_3d) ) then
          if( IO_L ) write(IO_FID_LOG,*) dinfo%VAR_3d(:,:,:)
       endif
    endif

    return
  end subroutine SNO_vars_read


  !-----------------------------------------------------------------------------
  subroutine SNO_vars_write( &
       dirpath,       &
       basename,      &
       output_grads,  &
       nowrank,       &
       nowstep,       &
       nprocs_x_out,  &
       nprocs_y_out,  &
       nhalos_x,      &
       nhalos_y,      &
       hinfo,         &
       naxis,         &
       ainfo,         &
       dinfo,         &
       debug          )
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo,   &
       iteminfo
    use mod_sno_grads, only: &
       SNO_grads_write
    implicit none

    character(len=*), intent(in)    :: dirpath                               ! directory path                     (output)
    character(len=*), intent(in)    :: basename                              ! basename of file                   (output)
    logical,          intent(in)    :: output_grads
    integer,          intent(in)    :: nowrank                               ! current rank                       (output)
    integer,          intent(in)    :: nowstep                               ! current step                       (output)
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (output)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (output)
    integer,          intent(in)    :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)    :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    integer,          intent(in)    :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)    :: ainfo(naxis)                          ! axis information                   (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    logical,          intent(in)    :: debug
    !---------------------------------------------------------------------------

    if ( output_grads ) then
       call SNO_grads_write( dirpath,  & ! [IN]
                             nowstep,  & ! [IN]
                             hinfo,    & ! [IN]
                             naxis,    & ! [IN]
                             ainfo(:), & ! [IN]
                             dinfo,    & ! [IN]
                             debug     ) ! [IN]
    else
       call SNO_vars_write_netcdf( dirpath,                    & ! [IN]
                                   basename,                   & ! [IN]
                                   nowrank,                    & ! [IN]
                                   nowstep,                    & ! [IN]
                                   nprocs_x_out, nprocs_y_out, & ! [IN]
                                   nhalos_x,     nhalos_y,     & ! [IN]
                                   hinfo,                      & ! [IN]
                                   naxis,                      & ! [IN]
                                   ainfo(:),                   & ! [IN]
                                   dinfo,                      & ! [IN]
                                   debug                       ) ! [IN]
    endif

    return
  end subroutine SNO_vars_write

  !-----------------------------------------------------------------------------
  subroutine SNO_vars_write_netcdf( &
       dirpath,       &
       basename,      &
       nowrank,       &
       nowstep,       &
       nprocs_x_out,  &
       nprocs_y_out,  &
       nhalos_x,      &
       nhalos_y,      &
       hinfo,         &
       naxis,         &
       ainfo,         &
       dinfo,         &
       debug          )
    use scale_file_h, only: &
       FILE_REAL4, &
       FILE_REAL8
    use scale_file, only: &
       FILE_Create,         &
       FILE_Def_Variable, &
       FILE_Set_Attribute,   &
       FILE_EndDef,         &
       FILE_Write
    use scale_process, only: &
       PRC_MPIstop
    use mod_sno_h, only: &
       commoninfo, &
       axisinfo,   &
       iteminfo
    use mod_sno, only: &
       SNO_attributes_write
    use mod_sno_axis, only: &
       SNO_axis_define, &
       SNO_axis_write
    implicit none

    character(len=*), intent(in)    :: dirpath                               ! directory path                     (output)
    character(len=*), intent(in)    :: basename                              ! basename of file                   (output)
    integer,          intent(in)    :: nowrank                               ! current rank                       (output)
    integer,          intent(in)    :: nowstep                               ! current step                       (output)
    integer,          intent(in)    :: nprocs_x_out                          ! x length of 2D processor topology  (output)
    integer,          intent(in)    :: nprocs_y_out                          ! y length of 2D processor topology  (output)
    integer,          intent(in)    :: nhalos_x                              ! number of x-axis halo grids        (global domain)
    integer,          intent(in)    :: nhalos_y                              ! number of y-axis halo grids        (global domain)
    type(commoninfo), intent(in)    :: hinfo                                 ! common information                 (input)
    integer,          intent(in)    :: naxis                                 ! number of axis variables           (input)
    type(axisinfo),   intent(in)    :: ainfo(naxis)                          ! axis information                   (input)
    type(iteminfo),   intent(in)    :: dinfo                                 ! variable information               (input)
    logical,          intent(in)    :: debug

    character(len=H_LONG) :: basename_mod
    integer               :: fid
    logical               :: fileexisted
    integer               :: vid
    real(SP), allocatable :: VAR_1d_SP(:), VAR_2d_SP(:,:), VAR_3d_SP(:,:,:)
    real(DP), allocatable :: VAR_1d_DP(:), VAR_2d_DP(:,:), VAR_3d_DP(:,:,:)

    integer  :: gout1, gout2, gout3
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( dirpath == '' ) then
       basename_mod = trim(basename)
    else
       basename_mod = trim(dirpath)//'/'//trim(basename)
    endif

    call FILE_Create( basename_mod,                 & ! [IN]
                      hinfo%title,                  & ! [IN]
                      hinfo%source,                 & ! [IN]
                      hinfo%institute,              & ! [IN]
                      hinfo%grid_name,              & ! [IN]
                      fid,                          & ! [OUT]
                      fileexisted,                  & ! [OUT]
                      rankid     = nowrank,         & ! [IN]
                      time_units = dinfo%time_units ) ! [IN]

    if ( .NOT. fileexisted ) then ! do below only once when file is created

       call SNO_axis_define( fid,      & ! [IN]
                             naxis,    & ! [IN]
                             ainfo(:), & ! [IN]
                             debug     ) ! [IN]

       call SNO_attributes_write( fid,          & ! [IN]
                                  nowrank,      & ! [IN]
                                  nprocs_x_out, & ! [IN]
                                  nprocs_y_out, & ! [IN]
                                  nhalos_x,     & ! [IN]
                                  nhalos_y,     & ! [IN]
                                  hinfo,        & ! [IN]
                                  dinfo,        & ! [IN]
                                  debug         ) ! [IN]

    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** + + + define variable'

    if ( dinfo%dt > 0.0_DP ) then
       call FILE_Def_Variable( fid,                & ! [IN]
                               dinfo%varname,      & ! [IN]
                               dinfo%description,  & ! [IN]
                               dinfo%units,        & ! [IN]
                               dinfo%dim_rank,     & ! [IN]
                               dinfo%dim_name,     & ! [IN]
                               dinfo%datatype,     & ! [IN]
                               vid,                & ! [OUT]
                               time_int = dinfo%dt ) ! [IN]
    else
       call FILE_Def_Variable( fid,               & ! [IN]
                               dinfo%varname,     & ! [IN]
                               dinfo%description, & ! [IN]
                               dinfo%units,       & ! [IN]
                               dinfo%dim_rank,    & ! [IN]
                               dinfo%dim_name,    & ! [IN]
                               dinfo%datatype,    & ! [IN]
                               vid                ) ! [OUT]
    endif

    if ( hinfo%minfo_mapping_name /= "" ) then
       call FILE_Set_Attribute( fid, dinfo%varname, "grid_mapping", hinfo%minfo_mapping_name )
    endif

    call FILE_EndDef( fid )

    if ( .NOT. fileexisted ) then ! do below only once when file is created

       call SNO_axis_write( fid,      & ! [IN]
                            naxis,    & ! [IN]
                            ainfo(:), & ! [IN]
                            debug     ) ! [IN]

    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** + + + write variable'

    if ( dinfo%dim_rank == 1 ) then

       gout1 = size(dinfo%VAR_1d(:),1)

       if ( dinfo%datatype == FILE_REAL4 ) then

          allocate( VAR_1d_SP(gout1) )
          VAR_1d_SP(:) = real(dinfo%VAR_1d(:),kind=SP)

          call FILE_Write( vid,                       & ! [IN]
                           VAR_1d_SP(:),              & ! [IN]
                           dinfo%time_start(nowstep), & ! [IN]
                           dinfo%time_end  (nowstep)  ) ! [IN]

          deallocate( VAR_1d_SP )

       elseif( dinfo%datatype == FILE_REAL8 ) then

          allocate( VAR_1d_DP(gout1) )
          VAR_1d_DP(:) = real(dinfo%VAR_1d(:),kind=DP)

          call FILE_Write( vid,                       & ! [IN]
                           VAR_1d_DP(:),              & ! [IN]
                           dinfo%time_start(nowstep), & ! [IN]
                           dinfo%time_end  (nowstep)  ) ! [IN]

          deallocate( VAR_1d_DP )

       endif

    elseif( dinfo%dim_rank == 2 ) then

       gout1 = size(dinfo%VAR_2d(:,:),1)
       gout2 = size(dinfo%VAR_2d(:,:),2)

       if ( dinfo%datatype == FILE_REAL4 ) then

          allocate( VAR_2d_SP(gout1,gout2) )
          VAR_2d_SP(:,:) = real(dinfo%VAR_2d(:,:),kind=SP)

          call FILE_Write( vid,                       & ! [IN]
                           VAR_2d_SP(:,:),            & ! [IN]
                           dinfo%time_start(nowstep), & ! [IN]
                           dinfo%time_end  (nowstep)  ) ! [IN]

          deallocate( VAR_2d_SP )

       elseif( dinfo%datatype == FILE_REAL8 ) then

          allocate( VAR_2d_DP(gout1,gout2) )
          VAR_2d_DP(:,:) = real(dinfo%VAR_2d(:,:),kind=DP)

          call FILE_Write( vid,                       & ! [IN]
                           VAR_2d_DP(:,:),            & ! [IN]
                           dinfo%time_start(nowstep), & ! [IN]
                           dinfo%time_end  (nowstep)  ) ! [IN]

          deallocate( VAR_2d_DP )

       endif

    elseif( dinfo%dim_rank == 3 ) then

       gout1 = size(dinfo%VAR_3d(:,:,:),1)
       gout2 = size(dinfo%VAR_3d(:,:,:),2)
       gout3 = size(dinfo%VAR_3d(:,:,:),3)

       if ( dinfo%transpose ) then
          if ( dinfo%datatype == FILE_REAL4 ) then

             allocate( VAR_3d_SP(gout2,gout3,gout1) )
             do k = 1, gout1
             do j = 1, gout3
             do i = 1, gout2
                VAR_3d_SP(i,j,k) = real(dinfo%VAR_3d(k,i,j),kind=SP)
             enddo
             enddo
             enddo

             call FILE_Write( vid,                       & ! [IN]
                              VAR_3d_SP(:,:,:),          & ! [IN]
                              dinfo%time_start(nowstep), & ! [IN]
                              dinfo%time_end  (nowstep)  ) ! [IN]

             deallocate( VAR_3d_SP )

          elseif( dinfo%datatype == FILE_REAL8 ) then

             allocate( VAR_3d_DP(gout2,gout3,gout1) )
             do k = 1, gout1
             do j = 1, gout3
             do i = 1, gout2
                VAR_3d_DP(i,j,k) = real(dinfo%VAR_3d(k,i,j),kind=DP)
             enddo
             enddo
             enddo

             call FILE_Write( vid,                       & ! [IN]
                              VAR_3d_DP(:,:,:),          & ! [IN]
                              dinfo%time_start(nowstep), & ! [IN]
                              dinfo%time_end  (nowstep)  ) ! [IN]

             deallocate( VAR_3d_DP )

          endif
       else
          if ( dinfo%datatype == FILE_REAL4 ) then

             allocate( VAR_3d_SP(gout1,gout2,gout3) )
             VAR_3d_SP(:,:,:) = real(dinfo%VAR_3d(:,:,:),kind=SP)

             call FILE_Write( vid,                       & ! [IN]
                              VAR_3d_SP(:,:,:),          & ! [IN]
                              dinfo%time_start(nowstep), & ! [IN]
                              dinfo%time_end  (nowstep)  ) ! [IN]

             deallocate( VAR_3d_SP )

          elseif( dinfo%datatype == FILE_REAL8 ) then

             allocate( VAR_3d_DP(gout1,gout2,gout3) )
             VAR_3d_DP(:,:,:) = real(dinfo%VAR_3d(:,:,:),kind=DP)

             call FILE_Write( vid,                       & ! [IN]
                              VAR_3d_DP(:,:,:),          & ! [IN]
                              dinfo%time_start(nowstep), & ! [IN]
                              dinfo%time_end  (nowstep)  ) ! [IN]

             deallocate( VAR_3d_DP )

          endif
       endif

    endif

    return
  end subroutine SNO_vars_write_netcdf

end module mod_sno_vars
