!-------------------------------------------------------------------------------
!> Program FIO ico2ll
!!
!! @par Description
!!          This program converts from data on icosahedral grid (new I/O format)
!!          to that on latitude-longitude grid.
!!          ( packaged NICAM data format : PaNDa )
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program fio_ico2ll
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use mod_io_param
  use scale_stdio
  use scale_prof
  use scale_process, only: &
     PRC_LOCAL_COMM_WORLD, &
     PRC_nprocs,           &
     PRC_mpi_alive
  use scale_const, only: &
     CONST_UNDEF, &
     CONST_UNDEF4
  use scale_calendar, only: &
     CALENDAR_daysec2date,   &
     CALENDAR_adjust_daysec
  use mod_fio, only: &
     headerinfo, &
     datainfo
  use mod_mnginfo_light, only: &
     MNG_mnginfo_input,   &
     MNG_mnginfo_noinput, &
     MNG_PALL,            &
     MNG_prc_tab
  use mod_netcdf, only: & ! [add] 13-04-18 C.Kodama
     NETCDF_handler,        &
     NETCDF_set_logfid,     &
     NETCDF_open_for_write, &
     NETCDF_write,          &
     NETCDF_close
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  integer, parameter :: max_nvar   = 500
  integer, parameter :: max_nstep  = 1500
  integer, parameter :: max_nlayer = 200

  integer, parameter :: flim = 1
  integer,      save :: fmax

  !--- NAMELIST
  integer                :: glevel              = -1
  integer                :: rlevel              = -1
  character(len=H_SHORT) :: grid_topology       = 'ICOSAHEDRON'
                                                ! 'LCP'
                                                ! 'MLCP'
  logical                :: complete            = .false.
  character(len=H_LONG)  :: mnginfo             = ''
  character(len=H_LONG)  :: layerfile_dir       = ''
  character(len=H_LONG)  :: llmap_base          = ''
  character(len=H_LONG)  :: topo_base           = ''
  character(len=H_LONG)  :: infile(flim)        = ''
  integer                :: step_str            = 1
  integer                :: step_end            = max_nstep
  character(len=H_LONG)  :: outfile_dir         = '.'
  character(len=H_SHORT) :: outfile_prefix      = ''
  integer                :: outfile_rec         = 1
  logical                :: lon_swap            = .false.
  logical                :: use_NearestNeighbor = .false.
  logical                :: devide_template     = .false.
  logical                :: output_grads        = .true.
  logical                :: output_gtool        = .false.
  logical                :: output_netcdf       = .false.   ! [add] 13-04-18
  logical                :: datainfo_nodep_pe   = .true.    ! <- can be .true. if data header do not depend on pe.
  character(len=H_SHORT) :: selectvar(max_nvar) = ''
  integer                :: nlim_llgrid         = 10000000  ! limit number of lat-lon grid in 1 ico region
  logical                :: comm_smallchunk     = .true.    ! apply MPI_Allreduce for each k-layer?
  logical                :: dcmip2016           = .false.   ! CF mode for dcmip2016

  logical                :: help = .false.

  namelist /OPTION/ glevel,              &
                    rlevel,              &
                    grid_topology,       &
                    complete,            &
                    mnginfo,             &
                    layerfile_dir,       &
                    llmap_base,          &
                    topo_base,           &
                    infile,              &
                    step_str,            &
                    step_end,            &
                    outfile_dir,         &
                    outfile_prefix,      &
                    outfile_rec,         &
                    lon_swap,            &
                    use_NearestNeighbor, &
                    devide_template,     &
                    output_grads,        &
                    output_gtool,        &
                    output_netcdf,       &  ! [add] 13-04-18
                    datainfo_nodep_pe,   &  ! [add] 13-04-18
                    selectvar,           &
                    nlim_llgrid,         &
                    comm_smallchunk,     &
                    dcmip2016,           &
                    help

  !-----------------------------------------------------------------------------
  character(len=H_LONG) :: infname   = ""
  character(len=H_LONG) :: outbase   = ""
  character(len=H_LONG) :: layerfile = ""
  integer               :: fmode
  integer               :: gtopology
  logical               :: allvar = .true.

  ! ll grid coordinate
  integer              :: imax, jmax
  real(8), allocatable :: lon(:), lat(:)
  real(8), allocatable :: lon_tmp(:) ! [add] 13-04-18

  ! ico2ll weight mapping
  integer              :: num_llgrid
  integer, allocatable :: nmax_llgrid(:,:)
  integer, allocatable :: lon_idx(:,:,:), lat_idx(:,:,:)
  integer, allocatable :: n1(:,:,:), n2(:,:,:), n3(:,:,:)
  real(8), allocatable :: w1(:,:,:), w2(:,:,:), w3(:,:,:)

  ! ico data information
  integer, allocatable :: ifid(:)
  integer, allocatable :: prc_tab_C(:)
  type(headerinfo) hinfo
  type(datainfo)   dinfo

  ! topography data information
  integer, allocatable :: ifid_topo(:)
  real(8), allocatable :: topo(:,:,:)

  integer                             :: num_of_data
  integer                             :: nvar
  character(len=H_SHORT), allocatable :: var_name(:)
  character(len=H_MID),   allocatable :: var_desc(:)
  character(len=H_SHORT), allocatable :: var_unit(:)
  character(len=H_SHORT)              :: var_name_nc
  character(len=H_MID)                :: var_desc_nc
  character(len=H_SHORT)              :: var_unit_nc
  character(len=H_SHORT), allocatable :: var_layername(:)
  integer,                allocatable :: var_datatype(:)
  integer,                allocatable :: var_nlayer(:)
  integer,                allocatable :: var_nstep(:)
  integer(8),             allocatable :: var_time_str(:)
  integer(8),             allocatable :: var_dt(:)
  logical,                allocatable :: var_xi2z(:)
  real(8),                allocatable :: var_ztop(:)
  real(8),                allocatable :: var_zgrid(:,:)
  ! header
  character(len=16),      allocatable :: var_gthead(:,:)
  ! NetCDF handler
  type(netcdf_handler)                :: nc              ! [add] 13-04-18
  character(len=1024)                 :: nc_time_units   ! [add] 13-04-18

  ! ico data
  integer              :: GALL
  integer              :: PALL_global
  integer              :: LALL_global
  integer              :: LALL_local

  real(4), allocatable :: data4allrgn(:)
  real(8), allocatable :: data8allrgn(:)
  real(4), allocatable :: icodata4   (:,:,:)
  real(4), allocatable :: icodata4_z (:,:)

  ! ll data
  real(4), allocatable :: lldata(:,:,:)
  real(4), allocatable :: temp  (:,:)
  real(4), allocatable :: lldata_total(:,:,:)

  ! for MPI
  integer          :: prc_nall, prc_nlocal
  integer          :: prc_myrank
  character(len=6) :: rankstr
  integer          :: pstr, pend, pp

  character(len=H_LONG) :: fname
  character(len=20)     :: tmpl
  character(len=16)     :: gthead(64)
  integer(8)            :: nowsec
  integer(8)            :: recsize ! [mod] 12-04-19 H.Yashiro
  integer               :: kmax, num_of_step, step, date_str(6)

  integer :: histday, offset_year
  real(8) :: histsec, histms

  logical :: addvar
  logical :: exist_topo
  integer :: ntemp
  real(8) :: wtemp
  integer :: fid, did, ofid, irec, ierr
  integer :: v, t, p, l, k, n, i, j, rgnid
  real(8) :: pi
  !=============================================================================

  pi = 4.D0 * atan( 1.D0 ) ! [add] 13-04-18

  !--- read option and preprocess
  call readoption !! set fmax, infile

  if ( glevel==-1 .or. rlevel==-1 ) then
     write(*,*) "xxx Set glevel, rlevel. STOP"
     stop
  endif
  if ( step_str < 1 .or. step_end < 1 ) then
     write(*,*) "xxx step must be >= 1. STOP"
     stop
  elseif( step_str > step_end ) then
     write(*,*) "xxx step_str must be < step_end. STOP"
     stop
  endif

  if ( grid_topology=="ICOSAHEDRON" ) then
     gtopology = IO_ICOSAHEDRON
  elseif( grid_topology=="LCP" ) then
     gtopology = IO_IGA_LCP
  elseif( grid_topology=="MLCP" ) then
     gtopology = IO_IGA_MLCP
  else
     write(*,*) "Unknown type of Grid toporogy:",grid_topology
     stop
  endif

  if    (output_gtool) then
     output_grads    = .false.
     devide_template = .false.
     outfile_rec     = 1
     output_netcdf   = .false.
     call PROF_rapstart('+FILE O GTOOL')
     call PROF_rapend  ('+FILE O GTOOL')
  elseif(output_netcdf) then
     output_grads    = .false.
     devide_template = .false.
     outfile_rec     = 1
     output_gtool    = .false.
     call PROF_rapstart('+FILE O NETCDF')
     call PROF_rapend  ('+FILE O NETCDF')
  elseif(output_grads) then
     call PROF_rapstart('+FILE O GRADS')
     call PROF_rapend  ('+FILE O GRADS')
  endif


  if ( trim(selectvar(1)) /= '' ) then
     allvar = .false.
  endif

  !#########################################################

  !--- prepare region infomation
  if (complete) then ! all region
    fmode = IO_INTEG_FILE
    call MNG_mnginfo_noinput( rlevel )
  else               ! region specified by mnginfo
    fmode = IO_SPLIT_FILE
    call MNG_mnginfo_input( rlevel, trim(mnginfo) )
  endif

  !--- Parallel Excution, No communication
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, prc_nall,   ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, prc_myrank, ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  PRC_mpi_alive        = .true.
  PRC_LOCAL_COMM_WORLD = MPI_COMM_WORLD
  PRC_nprocs           = prc_nall

  call PROF_rapstart('FIO_ICO2LL_MPI')

  ! borrow IO_FID_LOG to share log file id between other module
  IO_FID_LOG = IO_get_available_fid()
  if (output_netcdf) then
     call NETCDF_set_logfid( IO_FID_LOG )
  endif

  write(rankstr,'(I6.6)') prc_myrank
  open(IO_FID_LOG, file='LOG_ico2ll.pe'//trim(rankstr) )
  write(IO_FID_LOG,*) "+++ Parallel Execution, Use MPI"

  PALL_global = MNG_PALL
  LALL_global = 10 * (4**rlevel)
  LALL_local  = LALL_global / PALL_global

  if ( mod( PALL_global, prc_nall) /= 0 ) then
     write(IO_FID_LOG,*) "*** Invalid processor number, STOP:", PALL_global, prc_nall
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
     stop
  endif

  prc_nlocal = PALL_global / prc_nall
  pstr       = prc_myrank*prc_nlocal + 1
  pend       = prc_myrank*prc_nlocal + prc_nlocal
  write(IO_FID_LOG,*) "*** Number of Total .pexxxxxx files: ", PALL_global
  write(IO_FID_LOG,*) "*** Number of PE to packing precess: ", prc_nall
  write(IO_FID_LOG,*) "*** The rank of this process       : ", prc_myrank
  write(IO_FID_LOG,*) "*** Number of files for this rank  : ", prc_nlocal
  write(IO_FID_LOG,*) "*** file ID to pack                : ", pstr-1, " - ", pend-1

  !--- setup
  call fio_syscheck()

  !#########################################################

  write(IO_FID_LOG,*) '*** llmap read start'
  call PROF_rapstart('READ LLMAP')

  !--- Read lat-lon grid information
  fid = IO_get_available_fid()
  open(fid, file=trim(llmap_base)//'.info',form='unformatted',status='old',iostat=ierr)
     if (ierr/=0) then
        write(*,*) 'Cannot open llmap info file!',trim(llmap_base)//'.info'
        stop
     endif

     read(fid) imax
     allocate(lon(imax))
     read(fid) lon(:)
     read(fid) jmax
     allocate(lat(jmax))
     read(fid) lat(:)
  close(fid)

  !--- Read lat-lon weight map
  allocate( nmax_llgrid(LALL_local,prc_nlocal) )



  allocate( lon_idx(nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( lat_idx(nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( n1     (nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( n2     (nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( n3     (nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( w1     (nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( w2     (nlim_llgrid,LALL_local,prc_nlocal) )
  allocate( w3     (nlim_llgrid,LALL_local,prc_nlocal) )

  allocate( temp(imax,jmax) )

  ! read ll-ico relationship
  do p = pstr, pend
     pp = p - pstr + 1

     do l = 1, LALL_local
        rgnid = MNG_prc_tab(l,p)
        call IO_make_idstr(fname,trim(llmap_base),'rgn',rgnid,isrgn=.true.)
        write(IO_FID_LOG,*) 'p=', p, 'l=', l, 'rgnid=', rgnid

        fid = IO_get_available_fid()
        open(fid,file=trim(fname),form='unformatted',status='old',iostat=ierr)
           if (ierr/=0) then
              write(*,*) 'Cannot open llmap file!',trim(fname)
              stop
           endif

           read(fid) num_llgrid
           if ( num_llgrid > nlim_llgrid ) then
              write(*,*) 'less nlim_llgrid, please enlarge.',num_llgrid,' > ',nlim_llgrid
              stop
           endif
           nmax_llgrid(l,pp) = num_llgrid

           if ( num_llgrid /= 0 ) then
              read(fid) lon_idx( 1:num_llgrid,l,pp )
              read(fid) lat_idx( 1:num_llgrid,l,pp )
              read(fid) n1     ( 1:num_llgrid,l,pp )
              read(fid) n2     ( 1:num_llgrid,l,pp )
              read(fid) n3     ( 1:num_llgrid,l,pp )
              read(fid) w1     ( 1:num_llgrid,l,pp )
              read(fid) w2     ( 1:num_llgrid,l,pp )
              read(fid) w3     ( 1:num_llgrid,l,pp )
           endif
        close(fid)

        !--- sort weight (w1>w2>w3)
        if ( use_NearestNeighbor ) then
           do n = 1, nmax_llgrid(l,pp)
              if ( w3(n,l,pp) > w2(n,l,pp) ) then
                 wtemp      = w3(n,l,pp)
                 w3(n,l,pp) = w2(n,l,pp)
                 w2(n,l,pp) = wtemp
                 ntemp      = n3(n,l,pp)
                 n3(n,l,pp) = n2(n,l,pp)
                 n2(n,l,pp) = ntemp
              endif
              if ( w2(n,l,pp) > w1(n,l,pp) ) then
                 wtemp      = w2(n,l,pp)
                 w2(n,l,pp) = w1(n,l,pp)
                 w1(n,l,pp) = wtemp
                 ntemp      = n2(n,l,pp)
                 n2(n,l,pp) = n1(n,l,pp)
                 n1(n,l,pp) = ntemp
              endif
              if ( w3(n,l,pp) > w2(n,l,pp) ) then
                 wtemp      = w3(n,l,pp)
                 w3(n,l,pp) = w2(n,l,pp)
                 w2(n,l,pp) = wtemp
                 ntemp      = n3(n,l,pp)
                 n3(n,l,pp) = n2(n,l,pp)
                 n2(n,l,pp) = ntemp
              endif
           enddo
        endif
     enddo
  enddo

  call PROF_rapend('READ LLMAP')
  write(IO_FID_LOG,*) '*** llmap read end'

  !#########################################################

  write(IO_FID_LOG,*) '*** icodata read start'
  call PROF_rapstart('OPEN ICODATA')

  ! Read icodata information (all process)
  allocate( ifid(prc_nlocal) )

  allocate( prc_tab_C(LALL_local) )

  do p = pstr, pend
     pp = p - pstr + 1
     write(IO_FID_LOG,*) 'p=', pp

     if (complete) then ! all region
        infname = trim(infile(1))//'.rgnall'
     else
        call fio_mk_fname(infname,trim(infile(1)),'pe',p-1,6)
     endif
     prc_tab_C(1:LALL_local) = MNG_prc_tab(1:LALL_local,p)-1

     if ( pp == 1 ) then
        call fio_put_commoninfo( fmode,         &
                                 IO_BIG_ENDIAN, &
                                 gtopology,     &
                                 glevel,        &
                                 rlevel,        &
                                 LALL_local,    &
                                 prc_tab_C      )
     endif

     call fio_register_file(ifid(pp),trim(infname))
     call fio_fopen(ifid(pp),IO_FREAD)

     if ( datainfo_nodep_pe .AND. pp > 1 ) then ! assume that datainfo do not depend on pe.
        call fio_read_pkginfo          ( ifid(pp) )
        call fio_valid_pkginfo_validrgn( ifid(pp), prc_tab_C )
        call fio_copy_datainfo         ( ifid(pp), ifid(1)   )
     else ! normal way to read pkginfo and datainfo
        call fio_read_allinfo_validrgn ( ifid(pp), prc_tab_C )
     endif

  enddo

  call PROF_rapend('OPEN ICODATA')
  write(IO_FID_LOG,*) '*** icodata read end'

  !#########################################################

  write(IO_FID_LOG,*) '*** header check start'
  call PROF_rapstart('CHECK HEADER')

  !--- check all header
  allocate( hinfo%rgnid(LALL_local) )

  allocate( var_nstep    (max_nvar) )
  allocate( var_name     (max_nvar) )
  allocate( var_desc     (max_nvar) )
  allocate( var_unit     (max_nvar) )
  allocate( var_layername(max_nvar) )
  allocate( var_datatype (max_nvar) )
  allocate( var_nlayer   (max_nvar) )
  allocate( var_time_str (max_nvar) )
  allocate( var_dt       (max_nvar) )
  allocate( var_xi2z     (max_nvar) )
  allocate( var_ztop     (max_nvar) )
  allocate( var_zgrid    (max_nlayer, max_nvar) )
  allocate( var_gthead   (64, max_nvar) )

  pp = 1 ! only for first file

  call fio_get_pkginfo(ifid(pp),hinfo)
  num_of_data = hinfo%num_of_data

  nvar = 0
  do did = 0, num_of_data-1
     call fio_get_datainfo(ifid(pp),did,dinfo)

     if (allvar) then ! output all variables
        addvar = .true.
     else             ! select valiables to output
        addvar = .false.

        do v = 1, max_nvar
           if ( trim(selectvar(v)) == trim(dinfo%varname) ) then
              addvar = .true.
              exit
           elseif( trim(selectvar(v)) == '' ) then
              exit
           endif
        enddo
     endif

     do v = 1, nvar
        if ( trim(var_name(v)) == trim(dinfo%varname) ) then
           var_nstep(v) = var_nstep(v) + 1

           if( var_nstep(v) == 2 ) var_dt(v) = dinfo%time_start - var_time_str(v)

           if( var_nstep(v) == step_str ) var_time_str(v) = dinfo%time_start ! [mod] H.Yashiro 20111003

           addvar = .false.
           exit
        endif
     enddo

     if (addvar) then
        nvar = nvar + 1
        var_nstep    (nvar) = 1
        var_name     (nvar) = dinfo%varname
        var_desc     (nvar) = dinfo%description
        var_unit     (nvar) = dinfo%unit
        var_layername(nvar) = dinfo%layername
        var_datatype (nvar) = dinfo%datatype
        var_nlayer   (nvar) = dinfo%num_of_layer
        var_time_str (nvar) = dinfo%time_start
        var_dt       (nvar) = dinfo%time_end - dinfo%time_start
        var_xi2z     (nvar) = .false.
        var_ztop     (nvar) = CONST_UNDEF
        var_zgrid  (:,nvar) = CONST_UNDEF

        if ( prc_myrank == 0 ) then ! ##### only for master process

        if ( dinfo%layername == 'LAYERNM' ) then ! generate dummy
           do k = 1, dinfo%num_of_layer
              var_zgrid(k,nvar) = real(k,kind=8)
           enddo
        else ! read from file
           layerfile = trim(layerfile_dir)//'/'//trim(dinfo%layername)//'.txt'

           fid = IO_get_available_fid()
           open(fid,file=trim(layerfile),form='formatted',status='old',iostat=ierr)
              if ( ierr /= 0 ) then
                 write(*,*) 'xxx layerfile doesnt exist!', trim(layerfile)
                 stop
              endif

              read(fid,*) kmax
              do k = 1, kmax
                 read(fid,'(F16.4)') var_zgrid(k,nvar)
              enddo

              if ( dinfo%layername(1:5) == 'ZSALL' ) then ! check Xi2Z
                 write(IO_FID_LOG,*) '*** Try to convert Xi -> Z : ', dinfo%varname
                 var_xi2z(nvar) = .true.
                 var_ztop(nvar) = 0.5D0 * ( var_zgrid(kmax-1,nvar) + var_zgrid(kmax,nvar) )

                 if ( kmax == dinfo%num_of_layer+2 ) then ! trim HALO
                    write(IO_FID_LOG,*) '*** trim HALO: ', trim(dinfo%layername)
                    do k = 1, kmax-2
                       var_zgrid(k,nvar) = var_zgrid(k+1,nvar)
                    enddo
                 endif
              endif

           close(fid)
        endif

        endif ! ##### master?

     endif
  enddo !--- did LOOP

  GALL = ( (2**(glevel-rlevel))+2 ) &
       * ( (2**(glevel-rlevel))+2 )

  write(IO_FID_LOG,*) '*** get variable informations'
  write(IO_FID_LOG,*) 'num_of_data    : ', num_of_data

  if ( nvar == 0 ) then
     write(*,*) 'No variables to convert. Finish.'
     stop
  endif

  call PROF_rapend('CHECK HEADER')
  write(IO_FID_LOG,*) '*** header check end'

  !#########################################################

  write(IO_FID_LOG,*) '*** topography read start'
  call PROF_rapstart('READ TOPOGRAPHY')

  call PROF_rapstart('+Communication')
  ! broadcast var_xi2z, var_ztop and var_zgrid from master process
  call MPI_Bcast( var_xi2z(1),    &
                  max_nvar,       &
                  MPI_LOGICAL,    &
                  0,              &
                  MPI_COMM_WORLD, &
                  ierr            )

  call MPI_Bcast( var_ztop(1),    &
                  max_nvar,       &
                  MPI_REAL8,      &
                  0,              &
                  MPI_COMM_WORLD, &
                  ierr            )

  call MPI_Bcast( var_zgrid(1,1),      &
                  max_nlayer*max_nvar, &
                  MPI_REAL8,           &
                  0,                   &
                  MPI_COMM_WORLD,      &
                  ierr                 )
  call PROF_rapend  ('+Communication')

  if ( topo_base == '' ) then

     write(IO_FID_LOG,*) '*** topography file is not specified. no vertical conversion.'
     var_xi2z(:) = .false. ! reset flag

  else

     ! Read icodata (topography, all process)
     allocate( ifid_topo (prc_nlocal) )

     allocate( data4allrgn(GALL*LALL_local) )
     allocate( data8allrgn(GALL*LALL_local) )

     allocate( topo(GALL,LALL_local,prc_nlocal) )

     do p = pstr, pend
        pp = p - pstr + 1
        write(IO_FID_LOG,*) 'p=', pp

        prc_tab_C(1:LALL_local) = MNG_prc_tab(1:LALL_local,p)-1

        if (complete) then ! all region
           infname = trim(topo_base)//'.rgnall'
        else
           call fio_mk_fname(infname,trim(topo_base),'pe',p-1,6)
        endif

        call fio_register_file(ifid_topo(pp),trim(infname))
        call fio_fopen(ifid_topo(pp),IO_FREAD)

        if ( datainfo_nodep_pe .AND. pp > 1 ) then ! assume that datainfo do not depend on pe.
           call fio_read_pkginfo          ( ifid_topo(pp) )
           call fio_valid_pkginfo_validrgn( ifid_topo(pp), prc_tab_C    )
           call fio_copy_datainfo         ( ifid_topo(pp), ifid_topo(1) )
        else ! normal way to read pkginfo and datainfo
           call fio_read_allinfo_validrgn ( ifid_topo(pp), prc_tab_C )
        endif

        call fio_get_pkginfo(ifid_topo(pp),hinfo)
        num_of_data = hinfo%num_of_data

        exist_topo = .false.
        do did = 0, num_of_data-1
           call fio_get_datainfo(ifid_topo(pp),did,dinfo)

           if ( dinfo%varname == 'topo' ) then
              exist_topo = .true.

              !--- read from pe000xx file
              if ( dinfo%datatype == IO_REAL4 ) then

                 call fio_read_data(ifid_topo(pp),did,data4allrgn(:))
                 data8allrgn(:) = real(data4allrgn(:),kind=8)

              elseif( dinfo%datatype == IO_REAL8 ) then

                 call fio_read_data(ifid_topo(pp),did,data8allrgn(:))

              endif
              topo(:,:,pp) = reshape( data8allrgn(:), shape(topo(:,:,pp)) )

           endif
        enddo

        if ( .NOT. exist_topo ) then
           write(IO_FID_LOG,*) '*** topography data topo is not found in ', trim(infname), " ! STOP."
           stop
        endif

        call fio_fclose(ifid_topo(pp))
     enddo

     deallocate( ifid_topo )

     deallocate( data4allrgn )
     deallocate( data8allrgn )

  endif

  call PROF_rapend('READ TOPOGRAPHY')
  write(IO_FID_LOG,*) '*** topography read end'

  write(IO_FID_LOG,*) '########## Variable List ########## '
  write(IO_FID_LOG,*) 'ID |NAME            |STEPS|Layername       |START FROM         |DT [sec]|Xi2Z?'
  do v = 1, nvar
     histday     = 0
     histsec     = real(var_time_str(v),kind=RP)
     offset_year = 0

     call CALENDAR_adjust_daysec( histday, histsec ) ! [INOUT]

     call CALENDAR_daysec2date( date_str(:), & ! [OUT]
                                histms,      & ! [OUT]
                                histday,     & ! [IN]
                                histsec,     & ! [IN]
                                offset_year  ) ! [IN]

     write(tmpl,'(I4.4,"/",I2.2,"/",I2.2,1x,I2.2,":",I2.2,":",I2.2)') date_str(:)
     write(IO_FID_LOG,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A19,A1,I8,A1,L5)') &
              v,'|',var_name(v),'|',var_nstep(v),'|',var_layername(v),'|', tmpl,'|', var_dt(v), '|', var_xi2z(v)
  enddo

  !#########################################################

  write(IO_FID_LOG,*) '*** convert start : PaNDa format to lat-lon data'
  call PROF_rapstart('CONVERT')

  !--- start weighting summation
  do v = 1, nvar

     kmax    = var_nlayer(v)
     recsize = int(imax,kind=8)*int(jmax,kind=8)*int(kmax,kind=8)*4_8 ! [mod] 12-04-19 H.Yashiro

     num_of_step = min(step_end,var_nstep(v)) - step_str + 1

     allocate( data4allrgn(GALL*kmax*LALL_local) )
     allocate( data8allrgn(GALL*kmax*LALL_local) )
     allocate( icodata4   (GALL,kmax,LALL_local) )
     allocate( icodata4_z (GALL,kmax)            )

     allocate( lldata(imax,jmax,kmax) ) ! all node have large pallet

     allocate( lldata_total(imax,jmax,kmax) ) ! reduced

     if ( prc_myrank == 0 ) then ! ##### only for master process

     !--- open output file
     outbase = trim(outfile_dir)//'/'//trim(outfile_prefix)//trim(var_name(v))
     ofid    = IO_get_available_fid()

     if ( .NOT. devide_template ) then

        if (output_grads) then ! GrADS Format

           call PROF_rapstart('+FILE O GRADS')
           write(IO_FID_LOG,*)
           write(IO_FID_LOG,*) 'Output: ', trim(outbase)//'.grd', recsize, imax, jmax, kmax
           write(*         ,*) 'Output: ', trim(outbase)//'.grd'

           open( unit   = ofid,                  &
                 file   = trim(outbase)//'.grd', &
                 form   = 'unformatted',         &
                 access = 'direct',              &
                 recl   = recsize,               &
                 status = 'unknown'              )
           irec = 1

           if ( outfile_rec > 1 ) then
              write(IO_FID_LOG,*) 'Change output record position : start from step ', outfile_rec
              irec = outfile_rec
           endif
           call PROF_rapend  ('+FILE O GRADS')

        elseif(output_gtool) then ! GTOOL3 Format

           call PROF_rapstart('+FILE O GTOOL')
           write(IO_FID_LOG,*)
           write(IO_FID_LOG,*) 'Output: ', trim(outbase)//'.gt3', recsize, imax, jmax, kmax
           write(*         ,*) 'Output: ', trim(outbase)//'.gt3'

           open( unit   = ofid,                  &
                 file   = trim(outbase)//'.gt3', &
                 form   = 'unformatted',         &
                 access = 'sequential',          &
                 status = 'unknown'              )

           ! [mod] H.Yashiro 20111003
           call makegtoolheader( var_gthead(:,v),     &
                                 outfile_dir,         &
                                 var_name(v),         &
                                 var_desc(v),         &
                                 var_unit(v),         &
                                 var_layername(v),    &
                                 imax,                &
                                 jmax,                &
                                 var_nlayer(v),       &
                                 lon,                 &
                                 lat,                 &
                                 var_zgrid(1:kmax,v), &
                                 var_dt(v),           &
                                 lon_swap             )
           call PROF_rapend  ('+FILE O GTOOL')

        elseif(output_netcdf) then ! NetCDF format [add] 13-04-18 C.Kodama

           call PROF_rapstart('+FILE O NETCDF')
           write(IO_FID_LOG,*)
           write(IO_FID_LOG,*) 'Output: ', trim(outbase)//'.nc'
           write(*         ,*) 'Output: ', trim(outbase)//'.nc'

           histday     = 0
           histsec     = real(var_time_str(v),kind=RP)
           offset_year = 0

           call CALENDAR_adjust_daysec( histday, histsec ) ! [INOUT]

           call CALENDAR_daysec2date( date_str(:), & ! [OUT]
                                      histms,      & ! [OUT]
                                      histday,     & ! [IN]
                                      histsec,     & ! [IN]
                                      offset_year  ) ! [IN]

           write( nc_time_units,'(A14,I4.4,5(A,I2.2))') "minutes since ", date_str(1), &
                                                                     "-", date_str(2), &
                                                                     "-", date_str(3), &
                                                                     " ", date_str(4), &
                                                                     ":", date_str(5), &
                                                                     ":", date_str(6)

           write(IO_FID_LOG,*) '  nc_time_units = ', trim(nc_time_units)

           allocate( lon_tmp(imax) )

           if ( lon_swap ) then ! low_swap == .true. is not checked yet.
              lon_tmp(1:imax/2)      = ( lon(imax/2+1:imax)   ) * 180.D0 / pi
              lon_tmp(imax/2+1:imax) = ( lon(1:imax/2) + 2*pi ) * 180.D0 / pi
           else
              lon_tmp(:) = lon(:) * 180.D0 / pi
           endif

           if ( dcmip2016 ) then
              call cf_desc_unit( var_name_nc, & ! [OUT]
                                 var_desc_nc, & ! [OUT]
                                 var_unit_nc, & ! [OUT]
                                 var_name(v)  ) ! [IN]
           else
              var_name_nc = trim(var_name(v))
              var_desc_nc = trim(var_desc(v))
              var_unit_nc = trim(var_unit(v))
           endif

           call netcdf_open_for_write( nc,                                       & ! [OUT]
                                       ncfile      = trim(outbase)//'.nc',       & ! [IN]
                                       count       = (/  imax, jmax, kmax, 1 /), & ! [IN]
                                       title       = 'SCALE-GM data output',     & ! [IN]
                                       imax        = imax,                       & ! [IN]
                                       jmax        = jmax,                       & ! [IN]
                                       kmax        = kmax,                       & ! [IN]
                                       tmax        = num_of_step,                & ! [IN]
                                       lon         = lon_tmp,                    & ! [IN]
                                       lat         = (/ ( lat(j)*180.D0/pi, j=1, jmax ) /), & ! [IN]
                                       lev         = var_zgrid(1:kmax,v),        & ! [IN]
                                       time        = (/ (real(t-1,8)*real(var_dt(v),8)/real(60,8),t=1,num_of_step) /), & ! [IN]
                                       lev_units   ='m',                         & ! [IN]
                                       time_units  = trim(nc_time_units),        & ! [IN]
                                       var_name    = trim(var_name_nc),          & ! [IN]
                                       var_desc    = trim(var_desc_nc),          & ! [IN]
                                       var_units   = trim(var_unit_nc),          & ! [IN]
                                       var_missing = CONST_UNDEF4                 ) ! [IN]

           deallocate(lon_tmp)
           call PROF_rapend  ('+FILE O NETCDF')

        endif

     endif

     endif ! ##### master?

     do t = 1, num_of_step

        nowsec = var_time_str(v) + (t-1)*var_dt(v)
        step   = t-1 + step_str

        lldata(:,:,:) = 0.D0 ! cannot be filled by UNDEF because of reducing process

        if ( prc_myrank == 0 ) then ! ##### only for master process

        !--- open output file (every timestep)
        if (devide_template) then
           tmpl = sec2template(nowsec)
           write(IO_FID_LOG,*)
           write(IO_FID_LOG,*) 'Output: ', trim(outbase)//'.'//trim(tmpl)//'.grd'

           open( unit   = ofid,             &
                 file   = trim(outbase)//'.'//trim(tmpl)//'.grd', &
                 form   = 'unformatted',    &
                 access = 'direct',         &
                 recl   = recsize,          &
                 status = 'unknown'         )
           irec = 1
        endif

        endif ! ##### master?

        do p = pstr, pend
           pp = p - pstr + 1

           data4allrgn(:)  = CONST_UNDEF4
           data8allrgn(:)  = CONST_UNDEF
           icodata4(:,:,:) = CONST_UNDEF4

           call PROF_rapstart('+FILE I FIO')
           !--- seek data ID and get information
           call fio_seek_datainfo(did,ifid(pp),var_name(v),step)
           !--- verify
           if ( did == -1 ) then
              write(*,*) 'xxx data not found! varname:',trim(var_name(v)),", step : ",step
              stop
           endif

           !--- read from pe000xx file
           if ( var_datatype(v) == IO_REAL4 ) then

              call fio_read_data(ifid(pp),did,data4allrgn(:))

           elseif( var_datatype(v) == IO_REAL8 ) then

              call fio_read_data(ifid(pp),did,data8allrgn(:))

              data4allrgn(:) = real(data8allrgn(:),kind=4)
              where( data8allrgn(:) < CONST_UNDEF*0.1 )
                 data4allrgn(:) = CONST_UNDEF4
              endwhere

           endif
           icodata4(:,:,:) = reshape( data4allrgn(:), shape(icodata4) )
           call PROF_rapend('+FILE I FIO')

           do l = 1, LALL_local
              if ( t == 1 ) then
                 if ( mod(l,10) == 0 ) then
                    write(IO_FID_LOG,'(1x,I6.6)')              MNG_prc_tab(l,p)
                    write(IO_FID_LOG,'(A)',advance='no')       '          '
                 else
                    write(IO_FID_LOG,'(1x,I6.6)',advance='no') MNG_prc_tab(l,p)
                 endif
              endif

              !--- Zstar(Xi) -> Z coordinate
              if ( var_xi2z(v) ) then
                 call PROF_rapstart('+Xi2Z')

                 if ( var_ztop(v) < 0.D0 ) then
                    write(IO_FID_LOG,*) '*** Ztop is not specified.'
                    write(IO_FID_LOG,*) '*** It will be determined by the vertical axis info in ZSALL**.txt.'
                    stop
                 endif

                 call VINTRPL_Xi2Z ( GALL,               & ! [IN]
                                     kmax,               & ! [IN]
                                     var_zgrid (:,v),    & ! [IN]
                                     var_ztop  (v),      & ! [IN]
                                     topo      (:,l,pp), & ! [IN]
                                     CONST_UNDEF4,        & ! [IN]
                                     icodata4  (:,:,l),  & ! [IN]
                                     icodata4_z(:,:)     ) ! [OUT]

                 call PROF_rapend('+Xi2Z')
              else
                 icodata4_z(:,:) = icodata4(:,:,l)
              endif

              call PROF_rapstart('+Interpolation')
              !--- ico -> lat-lon
              if ( nmax_llgrid(l,pp) /= 0 ) then
                 if ( use_NearestNeighbor ) then ! nearest neighbor
                    do k = 1, kmax
                    do n = 1, nmax_llgrid(l,pp)
                       if    ( icodata4_z(n1(n,l,pp),k) /= CONST_UNDEF4 ) then

                          lldata(lon_idx(n,l,pp),lat_idx(n,l,pp),k) = icodata4_z(n1(n,l,pp),k)

                       elseif( icodata4_z(n2(n,l,pp),k) /= CONST_UNDEF4 ) then

                          lldata(lon_idx(n,l,pp),lat_idx(n,l,pp),k) = icodata4_z(n2(n,l,pp),k)

                       elseif( icodata4_z(n3(n,l,pp),k) /= CONST_UNDEF4 ) then

                          lldata(lon_idx(n,l,pp),lat_idx(n,l,pp),k) = icodata4_z(n3(n,l,pp),k)

                       else

                          lldata(lon_idx(n,l,pp),lat_idx(n,l,pp),k) = CONST_UNDEF4

                       endif
                    enddo
                    enddo
                 else
                    do k = 1, kmax
                    do n = 1, nmax_llgrid(l,pp)
                       if (      icodata4_z(n1(n,l,pp),k) < CONST_UNDEF4*0.1 &
                            .OR. icodata4_z(n2(n,l,pp),k) < CONST_UNDEF4*0.1 &
                            .OR. icodata4_z(n3(n,l,pp),k) < CONST_UNDEF4*0.1 ) then

                          lldata(lon_idx(n,l,pp),lat_idx(n,l,pp),k) = CONST_UNDEF4
                       else
                          lldata(lon_idx(n,l,pp),lat_idx(n,l,pp),k) = w1(n,l,pp) * icodata4_z(n1(n,l,pp),k) &
                                                                    + w2(n,l,pp) * icodata4_z(n2(n,l,pp),k) &
                                                                    + w3(n,l,pp) * icodata4_z(n3(n,l,pp),k)
                       endif
                    enddo
                    enddo
                 endif
              endif
              call PROF_rapend('+Interpolation')

           enddo ! region LOOP

           if ( t==1 ) write(IO_FID_LOG,*)

        enddo ! PE LOOP

        !--- swap longitude
        if (lon_swap) then
           do k = 1, kmax
              temp(1:imax/2,     :) = lldata(imax/2+1:imax,:,k)
              temp(imax/2+1:imax,:) = lldata(1:imax/2     ,:,k)
              lldata(:,:,k)         = temp(:,:)
           enddo
        endif

        call PROF_rapstart('+Communication')
        !--- Gather Lat-Lon data
        if ( comm_smallchunk ) then
           do k = 1, kmax
              call MPI_Allreduce( lldata      (1,1,k), &
                                  lldata_total(1,1,k), &
                                  imax*jmax,           &
                                  MPI_REAL,            &
                                  MPI_SUM,             &
                                  MPI_COMM_WORLD,      &
                                  ierr                 )
           enddo
        else
           call MPI_Allreduce( lldata      (1,1,1), &
                               lldata_total(1,1,1), &
                               imax*jmax*kmax,      &
                               MPI_REAL,            &
                               MPI_SUM,             &
                               MPI_COMM_WORLD,      &
                               ierr                 )
        endif
        call PROF_rapend  ('+Communication')

        if ( prc_myrank == 0 ) then ! ##### only for master process

        !--- output lat-lon data file
        if (output_grads) then

           call PROF_rapstart('+FILE O GRADS')
           write(ofid,rec=irec) lldata_total(:,:,:)
           irec = irec + 1
           call PROF_rapend  ('+FILE O GRADS')

        elseif(output_gtool) then

           call PROF_rapstart('+FILE O GTOOL')
           if ( nowsec < 2*365*24*60*60 ) then ! short term
              write(var_gthead(25,v),'(I16)') int(nowsec,kind=4)
              write(var_gthead(26,v),'(A16)') 'SEC             '
              write(var_gthead(28,v),'(I16)') int(var_dt(v),kind=4)
           else
              write(var_gthead(25,v),'(I16)') int( nowsec/60,kind=4 )
           endif
           write(var_gthead(27,v),'(A16)') CALENDAR_ss2cc_gtool(nowsec)
           gthead(:) = var_gthead(:,v)

           write(ofid) gthead(:)
           write(ofid) lldata_total(:,:,:)
           call PROF_rapend  ('+FILE O GTOOL')

        elseif(output_netcdf) then ! [add] 13.04.18 C.Kodama

           call PROF_rapstart('+FILE O NETCDF')
           call netcdf_write( nc, lldata_total(:,:,:), t=t )
!           do k=1,kmax
!              call netcdf_write( nc, lldata(:,:,k), k=k, t=t)
!           enddo
           call PROF_rapend  ('+FILE O NETCDF')

        endif

        !--- close output file
        if (devide_template) then
           close(ofid)
        endif

        write(IO_FID_LOG,*) ' +append step:', step

        endif ! ##### master?

     enddo ! step LOOP

     if ( prc_myrank == 0 ) then ! ##### only for master process

     !--- close output file
     if (.not. devide_template) then
        close(ofid)
     endif

     if (output_grads) then

        call PROF_rapstart('+FILE O GRADS')
        call makegradsctl( outfile_dir,         &
                           outfile_prefix,      &
                           var_name(v),         &
                           imax,                &
                           jmax,                &
                           kmax,                &
                           lon,                 &
                           lat,                 &
                           var_zgrid(1:kmax,v), &
                           num_of_step,         &
                           var_time_str(v),     &
                           var_dt(v),           &
                           lon_swap,            &
                           devide_template      )
        call PROF_rapend  ('+FILE O GRADS')

     elseif(output_netcdf) then ! [add] 13.04.18 C.Kodama

        call PROF_rapstart('+FILE O NETCDF')
        call netcdf_close( nc )
        call PROF_rapend  ('+FILE O NETCDF')

     endif

     endif ! ##### master?

     deallocate( data4allrgn )
     deallocate( data8allrgn )
     deallocate( icodata4    )
     deallocate( icodata4_z  )

     deallocate( lldata )
     deallocate( lldata_total )

  enddo ! variable LOOP

  do p = pstr, pend
     pp = p - pstr + 1
     call fio_fclose(ifid(pp))
  enddo ! PE LOOP

  call PROF_rapend('CONVERT')
  write(IO_FID_LOG,*) '*** convert finished! '

  call PROF_rapend('FIO_ICO2LL_MPI')
  call PROF_rapreport

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)

  close(IO_FID_LOG)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> read option
  subroutine readoption
    use mod_tool_option, only: &
      OPT_convert, &
      OPT_fid
    implicit none

    integer :: io
    !---------------------------------------------------------------------------

    ! --- Set option
    OPT_fid = IO_get_available_fid()
    open(OPT_fid,status='SCRATCH')

      call OPT_convert( fmax )

      read(OPT_fid,nml=OPTION,iostat=io)

    close(OPT_fid)

    if (      io /= 0     &
         .OR. fmax == 0   &
         .OR. fmax > flim &
         .OR. help        ) call helpoption

  end subroutine readoption

  !-----------------------------------------------------------------------------
  !> display help for option and abort
  subroutine helpoption
    implicit none
    !---------------------------------------------------------------------------

    write(*,OPTION)

    stop
  end subroutine helpoption

  !-----------------------------------------------------------------------------
  subroutine makegradsctl( &
      outfile_dir,    &
      outfile_prefix, &
      varname,        &
      imax,           &
      jmax,           &
      kmax,           &
      lon,            &
      lat,            &
      alt,            &
      nstep,          &
      time_str,       &
      dt,             &
      lon_swap,       &
      devide_template )
    implicit none

    character(len=*),  intent(in) :: outfile_dir
    character(len=16), intent(in) :: outfile_prefix
    character(len=16), intent(in) :: varname
    integer,           intent(in) :: imax
    integer,           intent(in) :: jmax
    integer,           intent(in) :: kmax
    real(8),           intent(in) :: lon(imax)
    real(8),           intent(in) :: lat(jmax)
    real(8),           intent(in) :: alt(kmax)
    integer,           intent(in) :: nstep
    integer(8),        intent(in) :: time_str
    integer(8),        intent(in) :: dt
    logical,           intent(in) :: lon_swap
    logical,           intent(in) :: devide_template

    real(8) :: pi
    real(8) :: temp(imax)

    character(len=32)  :: outfile
    integer            :: fid
    character(len=20)  :: s1, s2
    integer            :: i, j, k
    !---------------------------------------------------------------------------
    pi = 4.D0 * atan( 1.D0 )

    outfile = trim(outfile_prefix)//trim(var_name(v))

    fid = IO_get_available_fid()
    open( unit   = fid,         &
          file   = trim(outfile_dir)//'/'//trim(outfile)//'.ctl', &
          form   = 'formatted', &
          status = 'replace'    )

       if ( devide_template ) then
          write(fid,'(A)') 'DSET ^'//trim(outfile)//'.%y4-%m2-%d2-%h2h%n2m'//'.grd'
          write(fid,'(A)') 'OPTIONS TEMPLATE '
       else
          write(fid,'(A)') 'DSET ^'//trim(outfile)//'.grd'
       endif

       write(fid,'(A)')      'TITLE SCALE-GM data output'
       write(fid,'(A)')      'OPTIONS BIG_ENDIAN '
       write(fid,'(A,ES12.5)') 'UNDEF ', CONST_UNDEF4

       write(fid,'(A,I5,A)') 'XDEF ', imax, ' LEVELS'
       if (lon_swap) then
          temp(1:imax/2)      = lon(imax/2+1:imax)
          temp(imax/2+1:imax) = lon(1:imax/2) + 2*pi
          write(fid,'(10(1x,F9.4))') (temp(i)*180.D0/pi,i=1,imax)
       else
          write(fid,'(10(1x,F9.4))') (lon(i) *180.D0/pi,i=1,imax)
       endif

       write(fid,'(A,I5,A)')    'YDEF ',jmax, ' LEVELS'
       write(fid,'(10(1x,F9.4))') (lat(j)*180.D0/pi,j=1,jmax)

       if ( kmax == 1 ) then
          write(fid,'(A,I5,A,2I5)') 'ZDEF ', kmax, ' LINEAR', 1, 1
       else
          write(fid,'(A,I5,A)') 'ZDEF ', kmax, ' LEVELS'
          write(fid,'(10(1x,F9.2))') (alt(k),k=1,kmax)
       endif

       s1 = trim( sec2initplate(time_str) )
       s2 = trim( timeincrement(int(dt)) )
       write(fid,'(A,I5,2A,1x,A)') 'TDEF ',nstep, ' LINEAR ', trim(s1), trim(s2)

       write(fid,'(A,I5)') 'VARS ', 1
       if ( kmax == 1 ) then
          write(fid,'(A,2I5,1x,A)') trim(varname), 0,    99, 'NONE'
       else
          write(fid,'(A,2I5,1x,A)') trim(varname), kmax, 99, 'NONE'
       endif
       write(fid,'(A)') 'ENDVARS '
    close(fid)

    write(IO_FID_LOG,'(A,A)') 'Generate ',trim(outfile)//'.ctl'
    write(*         ,'(A,A)') 'Generate ',trim(outfile)//'.ctl'

  end subroutine makegradsctl

  !-----------------------------------------------------------------------------
  subroutine makegtoolheader( &
      gthead,      & !--- OUT
      outfile_dir, &
      varname,     &
      description, &
      unit,        &
      layername,   &
      imax,        &
      jmax,        &
      kmax,        &
      lon,         &
      lat,         &
      alt,         &
      dt,          &
      lon_swap     )
    implicit none

    character(len=16),      intent(out) :: gthead(64)
    character(len=H_LONG),  intent(in)  :: outfile_dir
    character(len=H_SHORT), intent(in)  :: varname
    character(len=H_MID),   intent(in)  :: description
    character(len=H_SHORT), intent(in)  :: unit
    character(len=H_SHORT), intent(in)  :: layername
    integer,                intent(in)  :: imax
    integer,                intent(in)  :: jmax
    integer,                intent(in)  :: kmax
    real(8),                intent(in)  :: lon(imax)
    real(8),                intent(in)  :: lat(jmax)
    real(8),                intent(in)  :: alt(kmax)
    integer(8),             intent(in)  :: dt
    logical,                intent(in)  :: lon_swap

    character(len=16) :: axhead(64)
    character(len=16) :: hitem
    character(len=32) :: htitle
    character(len=16) :: gt_axisx
    character(len=16) :: gt_axisy
    character(len=16) :: kdate

    integer           :: ndttm(8)
    character(len=10) :: ndate, ntime, nzone

    real(8) :: pi
    real(8) :: temp(imax)
    real(8) :: dx, lonp1(imax+1)
    integer :: i
    !---------------------------------------------------------------------------
    pi = 4.D0 * atan( 1.D0 )

    hitem  = trim(varname)
    do i=1,16
       if( hitem(i:i)=='_' ) hitem(i:i)  = '-' ! escape underbar
    enddo
    htitle(1:32) = description(1:32) ! trim to 32char
    do i=1,32
       if( htitle(i:i)=='_' ) htitle(i:i) = '-' ! escape underbar
    enddo

    write(gt_axisx,'(A,I4.4)') 'LON', imax
    write(gt_axisy,'(A,I4.4)') 'LAT', jmax

    call date_and_time(ndate, ntime, nzone, ndttm)
    write(kdate,'(I4.4,I2.2,I2.2,1x,I2.2,I2.2,I2.2,1x)') ndttm(1),ndttm(2),ndttm(3),ndttm(5),ndttm(6),ndttm(7)

    gthead(:) = ' '
    write(gthead( 1),'(I16)'  ) 9010
    write(gthead( 2),'(A16)'  ) 'SCALE-GM'
    write(gthead( 3),'(A16)'  ) hitem
    write(gthead(12),'(I16)'  ) 1
    write(gthead(13),'(I16)'  ) 1
    write(gthead(14),'(A16)'  ) htitle(1:16)
    write(gthead(15),'(A16)'  ) htitle(17:32)
    write(gthead(16),'(A16)'  ) unit

    write(gthead(26),'(A16)'  ) 'MIN             '
    write(gthead(28),'(I16)'  ) int(dt/60,kind=4)
    write(gthead(29),'(A16)'  ) gt_axisx ! from info file
    write(gthead(30),'(I16)'  ) 1
    write(gthead(31),'(I16)'  ) imax
    write(gthead(32),'(A16)'  ) gt_axisy ! from info file
    write(gthead(33),'(I16)'  ) 1
    write(gthead(34),'(I16)'  ) jmax
    write(gthead(35),'(A16)'  ) layername
    write(gthead(36),'(I16)'  ) 1
    write(gthead(37),'(I16)'  ) kmax
    write(gthead(38),'(A16)'  ) '             UR4'
    write(gthead(39),'(ES16.7)') CONST_UNDEF4
    write(gthead(40),'(ES16.7)') CONST_UNDEF4
    write(gthead(41),'(ES16.7)') CONST_UNDEF4
    write(gthead(42),'(ES16.7)') CONST_UNDEF4
    write(gthead(43),'(ES16.7)') CONST_UNDEF4
    write(gthead(44),'(I16)'  ) 1
    write(gthead(46),'(I16)'  ) 0
    write(gthead(47),'(ES16.7)') 0.
    write(gthead(48),'(I16)'  ) 0
    write(gthead(60),'(A16)'  ) kdate
    write(gthead(62),'(A16)'  ) kdate
    write(gthead(61),'(A16)'  ) 'SCALE-GM'
    write(gthead(63),'(A16)'  ) 'SCALE-GM'
    write(gthead(64),'(I16)'  ) imax*jmax*kmax

    !--- Generate axis file
    axhead(:) = ' '
    write(axhead( 1),'(I16)'  ) 9010
    write(axhead( 2),'(A16)'  ) 'AXLOC'
    write(axhead(12),'(I16)'  ) 1
    write(axhead(13),'(I16)'  ) 1
    write(axhead(25),'(I16)'  ) 0
    write(axhead(26),'(A16)'  ) 'SEC'
    write(axhead(27),'(A16)'  ) kdate
    write(axhead(28),'(I16)'  ) 1
    write(axhead(30),'(I16)'  ) 1
    write(axhead(33),'(I16)'  ) 1
    write(axhead(34),'(I16)'  ) 1
    write(axhead(36),'(I16)'  ) 1
    write(axhead(37),'(I16)'  ) 1
    write(axhead(38),'(A16)'  ) '             UR4'
    write(axhead(39),'(ES16.7)') -999.0
    write(axhead(44),'(I16)'  ) 1
    write(axhead(46),'(I16)'  ) 0
    write(axhead(47),'(ES16.7)') 0.
    write(axhead(48),'(I16)'  ) 0
    write(axhead(60),'(A16)'  ) kdate
    write(axhead(62),'(A16)'  ) kdate
    write(axhead(61),'(A16)'  ) 'SCALE-GM'
    write(axhead(63),'(A16)'  ) 'SCALE-GM'

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(gt_axisx)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then

       dx = lon(2)-lon(1)

       if (lon_swap) then
          temp(1:imax/2)      = lon(imax/2+1:imax)
          temp(imax/2+1:imax) = lon(1:imax/2) + 2*pi

          lonp1(1) = temp(1) - dx/2
       else
          lonp1(1) = lon(1)  - dx/2
       endif
       if ( abs(lonp1(1)) < 1.D-10 ) lonp1(1) = 0.D0

       do i = 2, imax+1
          lonp1(i) = lonp1(i-1) + dx
       enddo

       write(axhead( 3),'(A16)'  ) trim(gt_axisx)
       write(axhead(29),'(A16)'  ) trim(gt_axisx)
       write(axhead(31),'(I16)'  ) imax+1
       write(axhead(40),'(ES16.7)')   0.E0
       write(axhead(41),'(ES16.7)') 360.E0
       write(axhead(42),'(ES16.7)')  10.E0
       write(axhead(43),'(ES16.7)')  30.E0
       write(axhead(64),'(I16)'  ) imax+1

       write(fid) axhead
       write(fid) real( lonp1(1:imax+1)/pi*180.D0, kind=4 )

       close(fid)
    endif

    open( unit   = fid,          &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(gt_axisy)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then
       write(axhead( 3),'(A16)'  ) trim(gt_axisy)
       write(axhead(29),'(A16)'  ) trim(gt_axisy)
       write(axhead(31),'(I16)'  ) jmax
       write(axhead(40),'(ES16.7)') -90.E0
       write(axhead(41),'(ES16.7)')  90.E0
       write(axhead(42),'(ES16.7)')  10.E0
       write(axhead(43),'(ES16.7)')  30.E0
       write(axhead(64),'(I16)'  ) jmax

       write(fid) axhead
       write(fid) real( lat(1:jmax)/pi*180.D0, kind=4 )

       close(fid)
    endif

    open( unit   = fid,          &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(layername)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then
       write(axhead( 3),'(A16)'  ) trim(layername)
       write(axhead(29),'(A16)'  ) trim(layername)
       write(axhead(31),'(I16)'  ) kmax
       write(axhead(40),'(ES16.7)')     0.E0
       write(axhead(41),'(ES16.7)') real(maxval(alt),kind=4)
       write(axhead(42),'(ES16.7)')  1000.E0
       write(axhead(43),'(ES16.7)') 10000.E0
       if ( layername(1:5) == "STDPL" ) then
          write(axhead(44),'(I16)'  ) -2
       else
          write(axhead(44),'(I16)'  ) 1
       endif
       write(axhead(64),'(I16)'  ) kmax

       write(fid) axhead
       write(fid) real( alt(1:kmax), kind=4 )

       close(fid)
    endif

    return
  end subroutine makegtoolheader

  !-----------------------------------------------------------------------------
  !> output grads-like template part  like 01JAN0000
  function sec2initplate(datesec) result(template)
    implicit none

    integer(8)        :: datesec
    character(len=20) :: template

    integer :: d(6)

    character(len=3) :: nmonth(12)
    data nmonth / 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC' /
    !---------------------------------------------------------------------------

    histday     = 0
    histsec     = real(datesec,kind=RP)
    offset_year = 0

    call CALENDAR_adjust_daysec( histday, histsec ) ! [INOUT]

    call CALENDAR_daysec2date( d(:),        & ! [OUT]
                               histms,      & ! [OUT]
                               histday,     & ! [IN]
                               histsec,     & ! [IN]
                               offset_year  ) ! [IN]

    write(template,'(I2.2,A1,I2.2,A1,I2.2,A3,I4.4)') &
                              d(4), ':', d(5), 'Z', d(3), nmonth(d(2)), d(1)

  end function sec2initplate

  !-----------------------------------------------------------------------------
  !> output grads-like template part  like 2005-12-01-23h50m
  function sec2template(datesec) result(template)
    implicit none

    integer(8)        :: datesec
    character(len=20) :: template

    integer :: d(6)
    !---------------------------------------------------------------------------

    histday     = 0
    histsec     = real(datesec,kind=RP)
    offset_year = 0

    call CALENDAR_adjust_daysec( histday, histsec ) ! [INOUT]

    call CALENDAR_daysec2date( d(:),        & ! [OUT]
                               histms,      & ! [OUT]
                               histday,     & ! [IN]
                               histsec,     & ! [IN]
                               offset_year  ) ! [IN]

    write(template,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1)') &
                          d(1), '-', d(2), '-', d(3), '-', d(4), 'h', d(5), 'm'

  end function sec2template

  !-----------------------------------------------------------------------------
  function timeincrement(isec) result(template)
    implicit none

    integer           :: isec
    character(len=20) :: template

    character(len=18) :: tmp
    !---------------------------------------------------------------------------

    write(tmp,*) max(isec/60, 1)

    template = trim(tmp)//'mn'

  end function timeincrement

  !-----------------------------------------------------------------------------
  !> calendar, sec. -> character (YYYYMMDD HHMMSS)
  function CALENDAR_ss2cc_gtool(datesec) result(template)
    implicit none

    integer(8)        :: datesec
    character(len=16) :: template

    integer :: d(6), i
    !---------------------------------------------------------------------------

    histday     = 0
    histsec     = real(datesec,kind=RP)
    offset_year = 0

    call CALENDAR_adjust_daysec( histday, histsec ) ! [INOUT]

    call CALENDAR_daysec2date( d(:),        & ! [OUT]
                               histms,      & ! [OUT]
                               histday,     & ! [IN]
                               histsec,     & ! [IN]
                               offset_year  ) ! [IN]

    write (template,'(i4.4,i2.2,i2.2,1x,i2.2,i2.2,i2.2,1x)') (d(i),i=1,6)

  end function CALENDAR_ss2cc_gtool

  !-----------------------------------------------------------------------------
  subroutine VINTRPL_Xi2Z( &
      ijdim, &
      kdim,  &
      Xi,    &
      Ztop,  &
      Zsfc,  &
      UNDEF, &
      var_Z, &
      var_Xi )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    real(8), intent(in)  :: Xi       (kdim)
    real(8), intent(in)  :: Ztop
    real(8), intent(in)  :: Zsfc     (ijdim)
    real(4), intent(in)  :: UNDEF

    real(4), intent(in)  :: var_Z (ijdim,kdim)
    real(4), intent(out) :: var_Xi(ijdim,kdim)

    real(8) :: Z(ijdim,kdim)
    integer :: xi2z_idx (2)
    real(4) :: xi2z_coef(3)

    integer :: ij, k, kk
    !---------------------------------------------------------------------------

    do k  = 1, kdim
    do ij = 1, ijdim

       Z(ij,k) = Zsfc(ij) + ( Ztop - Zsfc(ij) ) / Ztop * Xi(k)

    enddo
    enddo

    do k  = 1, kdim
    do ij = 1, ijdim

       if ( Xi(k) <= Zsfc(ij) ) then

          xi2z_idx (1) = 1    ! dummmy
          xi2z_idx (2) = 1    ! dummmy
          xi2z_coef(1) = 0.0
          xi2z_coef(2) = 0.0
          xi2z_coef(3) = 1.0  ! set UNDEF

       elseif( Xi(k) <= Z(ij,1) ) then

          xi2z_idx (1) = 1    ! dummmy
          xi2z_idx (2) = 1
          xi2z_coef(1) = 0.0
          xi2z_coef(2) = 1.0
          xi2z_coef(3) = 0.0

       elseif( Xi(k) > Z(ij,kdim) ) then

          xi2z_idx (1) = kdim
          xi2z_idx (2) = kdim ! dummmy
          xi2z_coef(1) = 1.0
          xi2z_coef(2) = 0.0
          xi2z_coef(3) = 0.0

       elseif( Xi(k) > Ztop ) then

          xi2z_idx (1) = kdim ! dummmy
          xi2z_idx (2) = kdim ! dummmy
          xi2z_coef(1) = 0.0
          xi2z_coef(2) = 0.0
          xi2z_coef(3) = 1.0  ! set UNDEF

       else

          do kk = 2, kdim
             if( Xi(k) <= Z(ij,kk) ) exit
          enddo

          xi2z_idx (1) = kk-1
          xi2z_idx (2) = kk
          xi2z_coef(1) = ( Z (ij,kk) - Xi(k)       ) &
                       / ( Z (ij,kk) - Z (ij,kk-1) )
          xi2z_coef(2) = ( Xi(k)     - Z (ij,kk-1) ) &
                       / ( Z (ij,kk) - Z (ij,kk-1) )
          xi2z_coef(3) = 0.0

       endif

       if (       var_Z(ij,xi2z_idx(1)) <= UNDEF*0.1 &
            .AND. var_Z(ij,xi2z_idx(2)) <= UNDEF*0.1 ) then

          xi2z_coef(1) = 0.0
          xi2z_coef(2) = 0.0
          xi2z_coef(3) = 1.0

       elseif(    var_Z(ij,xi2z_idx(1)) <= UNDEF*0.1 ) then

          xi2z_coef(1) = 0.0
          xi2z_coef(2) = xi2z_coef(1) + xi2z_coef(2)

       elseif(    var_Z(ij,xi2z_idx(1)) <= UNDEF*0.1 ) then

          xi2z_coef(1) = xi2z_coef(1) + xi2z_coef(2)
          xi2z_coef(2) = 0.0

       endif

       var_Xi(ij,k) = xi2z_coef(1) * var_Z(ij,xi2z_idx(1)) &
                    + xi2z_coef(2) * var_Z(ij,xi2z_idx(2)) &
                    + xi2z_coef(3) * UNDEF
    enddo
    enddo

    return
  end subroutine VINTRPL_Xi2Z

  !-----------------------------------------------------------------------------
  subroutine cf_desc_unit( &
      var_name_nc, &
      var_desc_nc, &
      var_unit_nc, &
      var_name     )
    implicit none

    character(len=H_SHORT), intent(out) :: var_name_nc
    character(len=H_MID),   intent(out) :: var_desc_nc
    character(len=H_SHORT), intent(out) :: var_unit_nc
    character(len=H_SHORT), intent(in)  :: var_name
    !---------------------------------------------------------------------------

    select case( trim(var_name) )
    case( 'U', 'u' )
       var_name_nc = "U"
       var_desc_nc = "Zonal wind"
       var_unit_nc = "m/s"
    case( 'V', 'v' )
       var_name_nc = "V"
       var_desc_nc = "Meridional wind"
       var_unit_nc = "m/s"
    case( 'W', 'w' )
       var_name_nc = "W"
       var_desc_nc = "Vertical velocity"
       var_unit_nc = "m/s"
    case( 'PRS', 'prs' )
       var_name_nc = "P"
       var_desc_nc = "Pressure"
       var_unit_nc = "Pa"
    case( 'T', 't' )
       var_name_nc = "T"
       var_desc_nc = "Temperature"
       var_unit_nc = "K"
    case( 'PS', 'ps' )
       var_name_nc = "PS"
       var_desc_nc = "Surface pressure"
       var_unit_nc = "Pa"
    case( 'U500', 'u500' )
       var_name_nc = "U500"
       var_desc_nc = "Zonal wind at 500 hPa"
       var_unit_nc = "m/s"
    case( 'U850', 'u850' )
       var_name_nc = "U850"
       var_desc_nc = "Zonal wind at 850 hPa"
       var_unit_nc = "m/s"
    case( 'V500', 'v500' )
       var_name_nc = "V500"
       var_desc_nc = "Meridional wind at 500 hPa"
       var_unit_nc = "m/s"
    case( 'V850', 'v850' )
       var_name_nc = "V850"
       var_desc_nc = "Meridional wind at 850 hPa"
       var_unit_nc = "m/s"
    case( 'W500', 'w500' )
       var_name_nc = "W500"
       var_desc_nc = "Vertical velocity at 500 hPa"
       var_unit_nc = "m/s"
    case( 'W850', 'w850' )
       var_name_nc = "W850"
       var_desc_nc = "Vertical velocity at 850 hPa"
       var_unit_nc = "m/s"
    case( 'T500', 't500' )
       var_name_nc = "T500"
       var_desc_nc = "Temperature at 500 hPa"
       var_unit_nc = "K"
    case( 'T850', 't850' )
       var_name_nc = "T850"
       var_desc_nc = "Temperature at 850 hPa"
       var_unit_nc = "K"
    case( 'QV', 'qv' )
       var_name_nc = "Q"
       var_desc_nc = "Specific humidity"
       var_unit_nc = "kg/kg"
    case( 'QC', 'qc' )
       var_name_nc = "Qc"
       var_desc_nc = "Cloud water mixing ratio"
       var_unit_nc = "kg/kg"
    case( 'QR', 'qr' )
       var_name_nc = "Qr"
       var_desc_nc = "Rain water mixing ratio"
       var_unit_nc = "kg/kg"
    case( 'PASV1', 'pasv1' )
       var_name_nc = "Q1"
       var_desc_nc = "Singlet chlorine mixing ratio"
       var_unit_nc = "kg/kg"
    case( 'PASV2', 'pasv2' )
       var_name_nc = "Q2"
       var_desc_nc = "Chlorine gas mixing ratio"
       var_unit_nc = "kg/kg"
    case( 'PRCP', 'prcp' )
       var_name_nc = "PRECL"
       var_desc_nc = "Large-scale precipitation rate"
       var_unit_nc = "m/s"
    case( 'CL_COLUMN', 'cl_column' )
       var_name_nc = "Q1c"
       var_desc_nc = "Singlet chlorine mixing ratio (column)"
       var_unit_nc = "kg/kg"
    case( 'CL2_COLUMN', 'cl2_column' )
       var_name_nc = "Q2c"
       var_desc_nc = "Chlorine gas mixing ratio (column)"
       var_unit_nc = "kg/kg"
    case( 'CLY_COLUMN', 'cly_column' )
       var_name_nc = "Cly"
       var_desc_nc = "Cl and Cl2 the weighted sum (column)"
       var_unit_nc = "kg/kg"
    case( 'FORCING_VX', 'forcing_vx' )
       var_name_nc = "Fvx"
       var_desc_nc = "Forcing term of horizontal velocity: vx"
       var_unit_nc = "m/s-2"
    case( 'FORCING_VY', 'forcing_vy' )
       var_name_nc = "Fvy"
       var_desc_nc = "Forcing term of horizontal velocity: vx"
       var_unit_nc = "m/s-2"
    case( 'FORCING_VZ', 'forcing_vz' )
       var_name_nc = "Fvz"
       var_desc_nc = "Forcing term of horizontal velocity: vx"
       var_unit_nc = "m/s-2"
    case( 'FORCING_E', 'forcing_e' )
       var_name_nc = "Fe"
       var_desc_nc = "Forcing term of moist internal energy"
       var_unit_nc = "J/kg/s"
    case( 'FORCING_QV', 'forcing_qv' )
       var_name_nc = "Fqv"
       var_desc_nc = "Forcing term of specific humidity"
       var_unit_nc = "kg/kg/s"
    case( 'FORCING_QC', 'forcing_qc' )
       var_name_nc = "Fqc"
       var_desc_nc = "Forcing term of cloud water mixing ratio"
       var_unit_nc = "mkg/kg/s"
    case( 'FORCING_QR', 'forcing_qr' )
       var_name_nc = "Fqr"
       var_desc_nc = "Forcing term of cloud water mixing ratio"
       var_unit_nc = "kg/kg/s"
    case( 'FORCING_CL', 'forcing_cl' )
       var_name_nc = "Fcl"
       var_desc_nc = "Forcing term of Singlet chlorine mixing ratio"
       var_unit_nc = "kg/kg/s"
    case( 'FORCING_CL2', 'forcing_cl2' )
       var_name_nc = "Fcl2"
       var_desc_nc = "Forcing term of Chlorine gas mixing ratio"
       var_unit_nc = "kg/kg/s"
    case default
       var_name_nc = trim(var_name)
       var_desc_nc = "NIL"
       var_unit_nc = "NIL"
    end select

  end subroutine cf_desc_unit

end program fio_ico2ll
