!-------------------------------------------------------------------------------
!> Program SCALE-LES ver.3
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!!
!<
!-------------------------------------------------------------------------------
program spd2bin
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only : &
    IO_get_available_fid
  use mod_const, only : &
    CONST_UNDEF4, &
    CONST_UNDEF8
  use mod_time, only : &
    TIME_sec2date
  use mod_fileio_h, only : &
    FIO_HSHORT,      &
    FIO_HMID,        &
    FIO_HLONG,        &
    FIO_REAL4,       &
    FIO_REAL8,       &
    FIO_BIG_ENDIAN,  &
    FIO_CARTESIAN,   &
    FIO_SPLIT_FILE,  &
    FIO_MPIIO_NOUSE, &
    FIO_FREAD,       &
    headerinfo,      &
    datainfo
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  integer, parameter :: max_nvar   = 100
  integer, parameter :: max_nstep  = 2000
  integer, parameter :: max_nlayer = 1000

  integer, parameter :: flim = 1
  integer,      save :: fmax

  !--- NAMELIST
  character(LEN=FIO_HLONG)  :: infile(flim)        = ''
  integer                   :: step_str            = 1
  integer                   :: step_end            = max_nstep
  character(LEN=FIO_HLONG)  :: outfile_dir         = '.'
  character(LEN=FIO_HSHORT) :: outfile_prefix      = ''
  integer                   :: outfile_rec         = 1
  logical                   :: devide_template     = .false.
  logical                   :: output_grads        = .true.
  logical                   :: output_gtool        = .false.
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''

  logical                   :: help = .false.


  integer              :: PRC_nmax   = 1 !< total number of processors
  integer              :: PRC_NUM_X  = 1
  integer              :: PRC_NUM_Y  = 1
  integer, allocatable :: PRC_2Drank(:,:)

  integer              :: GRID_IMAX = 60    ! # of computational cells: x
  integer              :: GRID_JMAX = 60    ! # of computational cells: y
  integer              :: GRID_KMAX = 320   ! # of computational cells: z
  real(8)              :: GRID_DX   = 40.D0 ! center/face length [m]: x
  real(8)              :: GRID_DY   = 40.D0 ! center/face length [m]: y
  real(8)              :: GRID_DZ   = 40.D0 ! layer/interface length [m]: z
  real(8), allocatable :: GRID_CX(:)        ! center coordinate [m]: x
  real(8), allocatable :: GRID_CY(:)        ! center coordinate [m]: y
  real(8), allocatable :: GRID_CZ(:)        ! center coordinate [m]: z

  namelist /OPTION/ PRC_nmax,        &
                    PRC_NUM_X,       &
                    PRC_NUM_Y,       &
                    GRID_IMAX,       &
                    GRID_JMAX,       &
                    GRID_KMAX,       &
                    GRID_DX,         &
                    GRID_DY,         &
                    GRID_DZ,         &
                    infile,          &
                    step_str,        &
                    step_end,        &
                    outfile_dir,     &
                    outfile_prefix,  &
                    outfile_rec,     &
                    devide_template, &
                    output_grads,    &
                    output_gtool,    &
                    selectvar,       &
                    help
  !-----------------------------------------------------------------------------
  character(LEN=FIO_HLONG) :: infname   = ""
  character(LEN=FIO_HLONG) :: outbase   = ""
  logical                  :: allvar = .true.

  ! ico data information
  integer, allocatable :: ifid(:)
  type(headerinfo) hinfo 
  type(datainfo)   dinfo 

  integer                                :: num_of_data
  integer                                :: nvar
  character(LEN=FIO_HSHORT), allocatable :: var_name(:)
  character(LEN=FIO_HMID),   allocatable :: var_desc(:)
  character(LEN=FIO_HSHORT), allocatable :: var_unit(:)
  character(LEN=FIO_HSHORT), allocatable :: var_layername(:)
  integer,                   allocatable :: var_datatype(:)
  integer,                   allocatable :: var_nlayer(:)
  integer,                   allocatable :: var_nstep(:)
  integer(8),                allocatable :: var_time_str(:)
  integer(8),                allocatable :: var_dt(:)
  ! header
  character(LEN=16),         allocatable :: var_gthead(:,:)


  integer                  :: ISG, IEG   ! start/end of inner domain: x, global
  integer                  :: JSG, JEG   ! start/end of inner domain: y, global

  real(4), allocatable     :: data4allrgn(:)
  real(8), allocatable     :: data8allrgn(:)
  real(4), allocatable     :: icodata4(:,:,:)

  real(4), allocatable     :: lldata(:,:,:)

  character(LEN=28)        :: tmpl
  character(LEN=16)        :: gthead(64)
  integer(8)               :: nowsec
  integer                  :: kmax, num_of_step, step, date_str(6)
  real(8)                  :: msec_str

  logical :: addvar
  integer :: fid, did, ofid, ierr, irec
  integer :: v, t, p, k, i, j
  !=============================================================================

  !--- read option and preprocess
  call readoption !! set fmax, infile

  if ( step_str < 1 .or. step_end < 1 ) then
     write(*,*) "xxx step must be >= 1. STOP"
     stop
  elseif( step_str > step_end ) then
     write(*,*) "xxx step_str must be < step_end. STOP"
     stop
  endif

  if (output_gtool) then
     output_grads    = .false.
     devide_template = .false.
     outfile_rec = 1
  endif

  if ( trim(selectvar(1)) /= '' ) then
     allvar = .false.
  endif

  allocate( PRC_2Drank(-1:PRC_nmax-1,2) ); PRC_2Drank(:,:) = -1

  do p = 0, PRC_nmax-1
     PRC_2Drank(p,1) = mod(p,PRC_NUM_X)
     PRC_2Drank(p,2) = (p-PRC_2Drank(p,1)) / PRC_NUM_X
  enddo

  allocate( GRID_CX(GRID_IMAX*PRC_NUM_X) )
  allocate( GRID_CY(GRID_JMAX*PRC_NUM_Y) )
  allocate( GRID_CZ(GRID_KMAX          ) )

  ! horizontal coordinate: uniform interval
  do i = 1, GRID_IMAX*PRC_NUM_X
     GRID_CX(i) = GRID_DX * ( dble(i) - 0.5D0 )
  enddo

  do j = 1, GRID_JMAX*PRC_NUM_Y
     GRID_CY(j) = GRID_DY * ( dble(j) - 0.5D0 )
  enddo

  ! vertical coordinate: uniform interval
  do k = 1, GRID_KMAX
     GRID_CZ(k) = GRID_DZ * ( dble(k) - 0.5D0 )
  enddo

  !--- setup
  call fio_syscheck()

  !#########################################################
  ! Read data information

  allocate( ifid(0:PRC_nmax-1) )
  allocate( var_nstep    (max_nvar) )
  allocate( var_name     (max_nvar) )
  allocate( var_desc     (max_nvar) )
  allocate( var_unit     (max_nvar) )
  allocate( var_layername(max_nvar) )
  allocate( var_datatype (max_nvar) )
  allocate( var_nlayer   (max_nvar) )
  allocate( var_time_str (max_nvar) )
  allocate( var_dt       (max_nvar) )
  allocate( var_gthead   (64, max_nvar) )

  do p = 0, PRC_nmax-1
     call fio_mk_fname(infname,trim(infile(1)),'pe',p,6)

     call fio_put_commoninfo( FIO_MPIIO_NOUSE, &
                              FIO_SPLIT_FILE,  &
                              FIO_BIG_ENDIAN,  &
                              FIO_CARTESIAN,   &
                              int(GRID_DX),    &
                              GRID_IMAX,       &
                              1,               &
                              p                )

     call fio_register_file(ifid(p),trim(infname))
     call fio_fopen(ifid(p),FIO_FREAD)
     call fio_read_allinfo(ifid(p))

     if ( p == 0 ) then ! only once
        allocate( hinfo%rgnid(1) )

        call fio_get_pkginfo(ifid(p),hinfo)

        num_of_data = hinfo%num_of_data
        write(*,*) '*** get variable informations'
        write(*,*) 'num_of_data    : ', num_of_data

        nvar = 0
        do did = 0, num_of_data-1
           call fio_get_datainfo(ifid(p),did,dinfo)

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

                 if( var_nstep(v) == 2 ) var_dt(v) = dinfo%time_start * 1.D-3 - var_time_str(v)

                 if( var_nstep(v) == step_str ) var_time_str(v) = dinfo%time_start * 1.D-3 ! [mod] H.Yashiro 20111003

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
              var_time_str (nvar) = dinfo%time_start * 1.D-3
              var_dt       (nvar) = ( dinfo%time_end - dinfo%time_start ) * 1.D-3
           endif

        enddo !--- did LOOP
     endif !--- PE=000000
  enddo !--- PE LOOP

  if ( nvar == 0 ) then
     write(*,*) 'No variables to convert. Finish.'
     stop
  endif

  write(*,*) '########## Variable List ########## '
  write(*,*) 'ID |NAME            |STEPS|Layername       |START FROM                 |DT [sec]'
  do v = 1, nvar
     call TIME_sec2date( date_str(:), msec_str, real(var_time_str(v),kind=8) )
     write(tmpl,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') &
               date_str(1),'/',date_str(2),'/',date_str(3),' ', &
               date_str(4),':',date_str(5),':',date_str(6),' +', &
               msec_str
     write(*,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A27,A1,I8)') &
              v,'|',var_name(v),'|',var_nstep(v),'|',var_layername(v),'|', tmpl,'|', var_dt(v)
  enddo

  write(*,*) '*** convert start : PSD format to binary data'

  !#########################################################
  !--- start weighting summation
  do v = 1, nvar

     kmax = var_nlayer(v)

     !--- open output file
     outbase = trim(outfile_dir)//'/'//trim(outfile_prefix)//trim(var_name(v))
     ofid = IO_get_available_fid()
 
     if (.not. devide_template) then
        if (output_grads) then

           write(*,*) 'Output: ', trim(outbase)//'.grd'
           open( unit   = ofid,                  &
                 file   = trim(outbase)//'.grd', &
                 form   = 'unformatted',         &
                 access = 'direct',              &
                 recl   = GRID_IMAX*PRC_NUM_X*GRID_JMAX*PRC_NUM_Y*GRID_KMAX*4, &
                 status = 'unknown'              )
           irec = 1

           if ( outfile_rec > 1 ) then
              write(*,*) 'Change output record position : start from step ', outfile_rec
              irec = outfile_rec
           endif

        elseif(output_gtool) then

           write(*,*) 'Output: ', trim(outbase)//'.gt3'
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
                                 GRID_IMAX*PRC_NUM_X, &
                                 GRID_JMAX*PRC_NUM_Y, &
                                 GRID_KMAX,           &
                                 GRID_CX,             &
                                 GRID_CY,             &
                                 GRID_CZ,             &
                                 var_dt(v)            )

        endif
     endif

     num_of_step = min(step_end,var_nstep(v)) - step_str + 1

     do t = 1, num_of_step

        allocate( lldata(GRID_IMAX*PRC_NUM_X,GRID_JMAX*PRC_NUM_Y,GRID_KMAX) )
        lldata(:,:,:) = CONST_UNDEF4

        nowsec = var_time_str(v) + (t-1)*var_dt(v)

        !--- open output file (every timestep)
        if (devide_template) then
           tmpl   = sec2template(nowsec)
           write(*,*)
           write(*,*) 'Output: ', trim(outbase)//'.'//trim(tmpl)//'.grd'

           open( unit   = ofid,             &
                 file   = trim(outbase)//'.'//trim(tmpl)//'.grd', &
                 form   = 'unformatted',    &
                 access = 'direct',         &
                 recl   = GRID_IMAX*GRID_JMAX*GRID_KMAX*4, &
                 status = 'unknown'         )
           irec = 1 
        endif

        step = t-1 + step_str

        do p = 0, PRC_nmax-1

           allocate( data4allrgn(GRID_IMAX*GRID_JMAX*GRID_KMAX) )
           allocate( data8allrgn(GRID_IMAX*GRID_JMAX*GRID_KMAX) )
           allocate( icodata4   (GRID_IMAX,GRID_JMAX,GRID_KMAX) )
           data4allrgn(:)  = CONST_UNDEF4
           data8allrgn(:)  = CONST_UNDEF8
           icodata4(:,:,:) = CONST_UNDEF4

           !--- seek data ID and get information
           call fio_seek_datainfo(did,ifid(p),var_name(v),step)
           call fio_get_datainfo(ifid(p),did,dinfo)

           !--- verify
           if ( did == -1 ) then
              write(*,*) 'xxx data not found! varname:',trim(var_name(v)),", step : ",step
              stop
           endif

           !--- read from pe000xx file
           if ( dinfo%datatype == FIO_REAL4 ) then
              call fio_read_data(ifid(p),did,data4allrgn(:))
           elseif( dinfo%datatype == FIO_REAL8 ) then
              call fio_read_data(ifid(p),did,data8allrgn(:))

              data4allrgn(:) = real(data8allrgn(:),kind=4)
              where( data8allrgn(:) == CONST_UNDEF8 )
                 data4allrgn(:) = CONST_UNDEF4
              endwhere

           endif
           icodata4(:,:,:) = reshape( data4allrgn(:), shape(icodata4) )

           ! horizontal index (global domain)
           ISG = 1         + PRC_2Drank(p,1) * GRID_IMAX
           IEG = GRID_IMAX + PRC_2Drank(p,1) * GRID_IMAX
           JSG = 1         + PRC_2Drank(p,2) * GRID_JMAX
           JEG = GRID_JMAX + PRC_2Drank(p,2) * GRID_JMAX

           lldata(ISG:IEG,JSG:JEG,1:GRID_KMAX) = icodata4(1:GRID_IMAX,1:GRID_JMAX,1:GRID_KMAX)

           deallocate( data4allrgn )
           deallocate( data8allrgn )
           deallocate( icodata4    )
        enddo ! PE LOOP

        !--- output lat-lon data file
        if (output_grads) then
           write(ofid,rec=irec) lldata(:,:,:)
           irec = irec + 1
        elseif(output_gtool) then
           write(var_gthead(25,v),'(I16)') int( nowsec/3600,kind=4 )
           write(var_gthead(27,v),'(A16)') calendar_ss2cc_gtool(nowsec)
           gthead(:) = var_gthead(:,v)

           write(ofid) gthead(:)
           write(ofid) lldata(:,:,:)
        endif

        !--- close output file
        if (devide_template) then
           close(ofid)
        endif

        write(*,*) ' +append step:', step
        deallocate( lldata )
     enddo ! step LOOP

     !--- close output file
     if (.not. devide_template) then
        close(ofid)
     endif

     if (output_grads) then
        call makegradsctl( outfile_dir,         &
                           outfile_prefix,      &
                           var_name(v),         &
                           GRID_IMAX*PRC_NUM_X, &
                           GRID_JMAX*PRC_NUM_Y, &
                           GRID_KMAX,           &
                           GRID_CX,             &
                           GRID_CY,             &
                           GRID_CZ,             &
                           var_nstep(v),        &
                           var_time_str(v),     &
                           var_dt(v),           &
                           devide_template      )
     endif

  enddo ! variable LOOP

  do p = 0, PRC_nmax-1
     call fio_fclose(ifid(p))
  enddo

contains
  !-----------------------------------------------------------------------------
  !> read option
  !-----------------------------------------------------------------------------
  subroutine readoption
    use mod_stdio, only : &
       IO_get_available_fid
    use mod_option, only: &
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
  !-----------------------------------------------------------------------------
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
      devide_template )
    implicit none

    character(LEN=128), intent(in) :: outfile_dir
    character(LEN=16),  intent(in) :: outfile_prefix
    character(LEN=16),  intent(in) :: varname
    integer,            intent(in) :: imax
    integer,            intent(in) :: jmax
    integer,            intent(in) :: kmax
    real(8),            intent(in) :: lon(imax)
    real(8),            intent(in) :: lat(jmax)
    real(8),            intent(in) :: alt(kmax)
    integer,            intent(in) :: nstep
    integer(8),         intent(in) :: time_str
    integer(8),         intent(in) :: dt
    logical,            intent(in) :: devide_template
 
    real(8) :: pi

    character(LEN=32)  :: outfile
    integer            :: fid
    character(LEN=20)  :: s1, s2
    integer            :: i, j, k
    !---------------------------------------------------------------------------
    pi = 4.D0 * atan( 1.D0 )

    outfile = trim(outfile_prefix)//trim(var_name(v))

    fid = IO_get_available_fid()
    open( unit   = fid,         &
          file   = trim(outfile_dir)//'/'//trim(outfile)//'.ctl', &
          form   = 'formatted', &
          status = 'replace'    )

       ! S.Iga051226=>
       if ( devide_template ) then
          ! W. Yanase 081008   use %h2 instead of %f2 for GrADS template
          write(fid,'(A)') 'DSET ^'//trim(outfile)//'.%y4-%m2-%d2-%h2h%n2m'//'.grd'
          write(fid,'(A)') 'OPTIONS TEMPLATE '
       else
          write(fid,'(A)') 'DSET ^'//trim(outfile)//'.grd'
       endif
       ! S.Iga051226<=

       write(fid,'(A)')      'TITLE SCALE3 data output'
       write(fid,'(A)')      'OPTIONS BIG_ENDIAN '
       write(fid,'(A,E12.5)') 'UNDEF ', real( -99.9E+33, kind=4 )

       write(fid,'(A,I5,A)') 'XDEF ', imax, ' LEVELS'
       write(fid,'(10(1x,F9.4))') (lon(i),i=1,imax)

       write(fid,'(A,I5,A)')    'YDEF ',jmax, ' LEVELS'
       write(fid,'(10(1x,F9.4))') (lat(j),j=1,jmax)

       if ( kmax == 1 ) then
          write(fid,'(A,I5,A,2I5)') 'ZDEF ', kmax, ' LINEAR', 1, 1
       else
          write(fid,'(A,I5,A)') 'ZDEF ', kmax, ' LEVELS'
          write(fid,'(10(1x,F9.2))') (alt(k),k=1,kmax)
       endif

       s1 = trim( sec2initplate(time_str) ) ! S.Iga060508
       s2 = trim( timeincrement(int(dt)) )  ! S.Iga060508
       write(fid,'(A,I5,2A,1x,A)') 'TDEF ',nstep, ' LINEAR ', trim(s1), trim(s2)

       write(fid,'(a,i5)') 'VARS ', 1
       if ( kmax == 1 ) then
          write(fid,'(a,2i5,1x,a)') trim(varname), 0, 99, 'NONE'
       else
          write(fid,'(a,2i5,1x,a)') trim(varname), kmax, 99, 'NONE'
       endif
       write(fid,'(a)') 'ENDVARS '
    close(fid)

    write(*,'(A,A)') 'Generate ',trim(outfile)//'.ctl'

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
      dt           )
    implicit none

    character(LEN=16),         intent(out) :: gthead(64)
    character(LEN=FIO_HLONG),  intent( in) :: outfile_dir
    character(LEN=FIO_HSHORT), intent( in) :: varname
    character(LEN=FIO_HMID),   intent( in) :: description
    character(LEN=FIO_HSHORT), intent( in) :: unit
    character(LEN=FIO_HSHORT), intent( in) :: layername
    integer,                   intent( in) :: imax
    integer,                   intent( in) :: jmax
    integer,                   intent( in) :: kmax
    real(8),                   intent( in) :: lon(imax)
    real(8),                   intent( in) :: lat(jmax)
    real(8),                   intent( in) :: alt(kmax)
    integer(8),                intent( in) :: dt
 
    character(LEN=16) :: axhead(64)
    character(LEN=16) :: hitem
    character(LEN=32) :: htitle
    character(LEN=16) :: gt_axisx
    character(LEN=16) :: gt_axisy
    character(LEN=16) :: kdate

    integer           :: ndttm(8)
    character(LEN=10) :: ndate, ntime, nzone

    real(8) :: pi
    real(8) :: dx, lonp1(imax+1)
    integer :: i
    !---------------------------------------------------------------------------
    pi = 4.D0 * atan( 1.D0 )

    hitem  = trim(varname)
    do i=1,16
       if( hitem(i:i)=='_' ) hitem(i:i)  = '-' ! escape underbar
    enddo
    htitle = description ! trim to 32char
    do i=1,32
       if( htitle(i:i)=='_' ) htitle(i:i) = '-' ! escape underbar
    enddo

    write(gt_axisx,'(A,I4.4)') 'LON', imax
    write(gt_axisy,'(A,I4.4)') 'LAT', jmax

    call date_and_time(ndate, ntime, nzone, ndttm)
    write(kdate,'(I4.4,I2.2,I2.2,1x,I2.2,I2.2,I2.2,1x)') ndttm(1),ndttm(2),ndttm(3),ndttm(5),ndttm(6),ndttm(7)

    gthead(:) = ' '
    write(gthead( 1),'(I16)'  ) 9010
    write(gthead( 2),'(A16)'  ) 'SCALE3'
    write(gthead( 3),'(A16)'  ) hitem
    write(gthead(12),'(I16)'  ) 1
    write(gthead(13),'(I16)'  ) 1
    write(gthead(14),'(A16)'  ) htitle(1:16)
    write(gthead(15),'(A16)'  ) htitle(17:32)
    write(gthead(16),'(A16)'  ) unit

    write(gthead(26),'(A16)'  ) 'HOUR            '
    write(gthead(28),'(I16)'  ) int(dt/3600,kind=4)
    write(gthead(29),'(A16)'  ) gt_axisx ! from info file
    write(gthead(30),'(I16)'  ) 1
    write(gthead(31),'(I16)'  ) imax
    write(gthead(32),'(A16)'  ) gt_axisy ! from info file
    write(gthead(33),'(I16)'  ) 1
    write(gthead(34),'(I16)'  ) jmax
    write(gthead(35),'(A16)'  ) layername
    write(gthead(36),'(I16)'  ) 1
    write(gthead(37),'(I16)'  ) kmax
    write(gthead(38),'(A16)'  ) 'UR4'
    write(gthead(39),'(E16.7)') real(-99.9E+33,4)
    write(gthead(40),'(E16.7)') real(-99.9E+33,4)
    write(gthead(41),'(E16.7)') real(-99.9E+33,4)
    write(gthead(42),'(E16.7)') real(-99.9E+33,4)
    write(gthead(43),'(E16.7)') real(-99.9E+33,4)
    write(gthead(44),'(I16)'  ) 1
    write(gthead(46),'(I16)'  ) 0
    write(gthead(47),'(E16.7)') 0. 
    write(gthead(48),'(I16)'  ) 0
    write(gthead(60),'(A16)'  ) kdate
    write(gthead(62),'(A16)'  ) kdate
    write(gthead(61),'(A16)'  ) 'SCALE3'
    write(gthead(63),'(A16)'  ) 'SCALE3'
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
    write(axhead(38),'(A16)'  ) 'UR4'
    write(axhead(39),'(E16.7)') -999.0
    write(axhead(44),'(I16)'  ) 1
    write(axhead(46),'(I16)'  ) 0
    write(axhead(47),'(E16.7)') 0. 
    write(axhead(48),'(I16)'  ) 0
    write(axhead(60),'(A16)'  ) kdate
    write(axhead(62),'(A16)'  ) kdate
    write(axhead(61),'(A16)'  ) 'SCALE3'
    write(axhead(63),'(A16)'  ) 'SCALE3'

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(trim(outfile_dir)//'/GTAXLOC.'//trim(gt_axisx)),&
          form   = 'unformatted', &
          access = 'sequential',  &
          status = 'new',         &
          iostat = ierr           )

    if (ierr == 0) then

       dx = lon(2)-lon(1)

       lonp1(1) = lon(1)  - dx/2
       if ( abs(lonp1(1)) < 1.D-10 ) lonp1(1) = 0.D0

       do i = 2, imax+1
          lonp1(i) = lonp1(i-1) + dx
       enddo

       write(axhead( 3),'(A16)'  ) trim(gt_axisx)
       write(axhead(29),'(A16)'  ) trim(gt_axisx)
       write(axhead(31),'(I16)'  ) imax+1
       write(axhead(40),'(E16.7)') lonp1(1)
       write(axhead(41),'(E16.7)') lonp1(imax+1)
       write(axhead(42),'(E16.7)')  40.E0
       write(axhead(43),'(E16.7)') 800.E0
       write(axhead(64),'(I16)'  ) imax+1

       write(fid) axhead
       write(fid) real( lonp1(1:imax+1), kind=4 )

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
       write(axhead(40),'(E16.7)') lat(1)
       write(axhead(41),'(E16.7)') lat(jmax)
       write(axhead(42),'(E16.7)')  40.E0
       write(axhead(43),'(E16.7)') 800.E0
       write(axhead(64),'(I16)'  ) jmax

       write(fid) axhead
       write(fid) real( lat(1:jmax), kind=4 )

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
       write(axhead(40),'(E16.7)')     0.E0
       write(axhead(41),'(E16.7)') real(maxval(alt),kind=4)
       write(axhead(42),'(E16.7)')  1000.E0
       write(axhead(43),'(E16.7)') 10000.E0
       write(axhead(64),'(I16)'  ) kmax

       write(fid) axhead
       write(fid) real( alt(1:kmax), kind=4 )

       close(fid)
    endif

    return
  end subroutine makegtoolheader

  !S.Iga051226 =>
  !-----------------------------------------------------------------------------
  function sec2initplate(datesec) result(template)
    !-- output grads-like template part  like 01JAN0000
    implicit none

    integer(8)        :: datesec
    ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: template
    character(LEN=20) :: template

    integer :: d(6)
    real(8) :: msec

    character(LEN=3) :: nmonth(12)
    data nmonth / 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC' /
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call TIME_sec2date( d(:), msec, real(datesec,kind=8) )

    write(template,'(I2.2,A1,I2.2,A1,I2.2,A3,I4.4)') &
                              d(4), ':', d(5), 'Z', d(3), nmonth(d(2)), d(1)

  end function sec2initplate

  !-----------------------------------------------------------------------------
  function sec2template(datesec) result(template)
    !-- output grads-like template part  like 2005-12-01-23h50m
    implicit none

    integer(8)        :: datesec
    ! [mod] 10/08/03 T.Mitsui, can be compiled by gfortran
!!$  character(*):: template
    character(LEN=20) :: template

    integer :: d(6)
    real(8) :: msec
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call TIME_sec2date( d(:), msec, real(datesec,kind=8) )

    write(template,'(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1)') &
                          d(1), '-', d(2), '-', d(3), '-', d(4), 'h', d(5), 'm'

  end function sec2template

  !-----------------------------------------------------------------------------
  function timeincrement(isec) result(template)
    implicit none

    integer       :: isec
    character(20) :: template

    character(18):: tmp
    !---------------------------------------------------------------------------

    write(tmp,*) isec/60

    template = trim(tmp)//'mn'

  end function timeincrement
  !S.Iga051226 <=

  !-----------------------------------------------------------------------------
  function calendar_ss2cc_gtool(datesec) result(template)
    !--- calendar, sec. -> character (YYYYMMDD HHMMSS)
    implicit none

    integer(8)        :: datesec
    character(LEN=16) :: template

    integer :: d(6), i
    real(8) :: msec
    !---------------------------------------------------------------------------

    ! [Comment] H.Yashiro 20110903
    ! Prefer not to use calendar_dd2ym subroutine
    ! Epoch time is different between calendar_ss2yh and calendar_dd2ym
    ! New I/O stores timestamp, which is generated via calendar_yh2ss
    call TIME_sec2date( d(:), msec, real(datesec,kind=8) )

    write (template,'(i4.4,i2.2,i2.2,1x,i2.2,i2.2,i2.2,1x)') (d(i),i=1,6)

  end function calendar_ss2cc_gtool

end program spd2bin
