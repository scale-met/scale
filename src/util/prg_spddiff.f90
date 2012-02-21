!-------------------------------------------------------------------------------
!> Program SPDDUMP
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          header/data veiwer for formatted data
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
program prg_spddiff
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
  !++ parameters & variables
  !
  integer, parameter :: max_nvar   = 100
  integer, parameter :: max_nstep  = 2000
  integer, parameter :: max_nlayer = 1000

  integer, parameter :: flim = 2
  integer,      save :: fmax

  !--- NAMELIST
  character(LEN=FIO_HLONG)  :: infile(flim)        = ''
  integer                   :: PRC_nmax            = 1 !< total number of processors
  integer                   :: PRC_NUM_X           = 1
  integer                   :: PRC_NUM_Y           = 1
  integer                   :: GRID_IMAX           = 60    ! # of computational cells: x
  integer                   :: GRID_JMAX           = 60    ! # of computational cells: y
  integer                   :: GRID_KMAX           = 320   ! # of computational cells: z
  real(8)                   :: GRID_DX             = 40.D0 ! center/face length [m]: x
  real(8)                   :: GRID_DY             = 40.D0 ! center/face length [m]: y
  real(8)                   :: GRID_DZ             = 40.D0 ! layer/interface length [m]: z
  integer                   :: step_str            = 1
  integer                   :: step_end            = max_nstep
  real(8)                   :: criterion           = 1.D-6 ! difference criterion
  character(LEN=FIO_HSHORT) :: selectvar(max_nvar) = ''

  logical                   :: help = .false.

  namelist /OPTION/ infile,          &
                    PRC_nmax,        &
                    PRC_NUM_X,       &
                    PRC_NUM_Y,       &
                    GRID_IMAX,       &
                    GRID_JMAX,       &
                    GRID_KMAX,       &
                    GRID_DX,         &
                    GRID_DY,         &
                    GRID_DZ,         &
                    step_str,        &
                    step_end,        &
                    criterion,       &
                    selectvar,       &
                    help
  !-----------------------------------------------------------------------------
  character(LEN=FIO_HLONG) :: infname  = ""
  logical                  :: allvar   = .true.

  ! data information
  integer, allocatable :: ifid1(:)
  integer, allocatable :: ifid2(:)
  type(headerinfo) hinfo 
  type(datainfo)   dinfo 

  integer                                :: num_of_data
  integer                                :: nvar1, nvar2
  character(LEN=FIO_HSHORT), allocatable :: var1_name(:), var2_name(:)
  character(LEN=FIO_HMID),   allocatable :: var1_desc(:), var2_desc(:)
  character(LEN=FIO_HSHORT), allocatable :: var1_unit(:), var2_unit(:)
  character(LEN=FIO_HSHORT), allocatable :: var1_layername(:), var2_layername(:)
  integer,                   allocatable :: var1_datatype(:), var2_datatype(:)
  integer,                   allocatable :: var1_nlayer(:), var2_nlayer(:)
  integer,                   allocatable :: var1_nstep(:), var2_nstep(:)
  integer(8),                allocatable :: var1_time_str(:), var2_time_str(:)
  integer(8),                allocatable :: var1_dt(:), var2_dt(:)

  character(LEN=28)    :: tmpl
  integer              :: date_str(6)

  real(4), allocatable     :: data1_r4(:)
  real(4), allocatable     :: data2_r4(:)
  real(4), allocatable     :: data1_r4_3d(:,:,:)
  real(4), allocatable     :: data2_r4_3d(:,:,:)
  real(8), allocatable     :: data1_r8(:)
  real(8), allocatable     :: data2_r8(:)
  real(8), allocatable     :: data1_r8_3d(:,:,:)
  real(8), allocatable     :: data2_r8_3d(:,:,:)

  logical :: addvar, filematch, headercheck, datacheck
  integer :: did, did1, did2
  integer :: v, p, k, i, j, v1, v2
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

  if ( trim(selectvar(1)) /= '' ) then
     allvar = .false.
  endif

  !--- setup
  call fio_syscheck()

  !#########################################################
  ! Read data information

  ! file 1
  allocate( ifid1(0:PRC_nmax-1)      )
  allocate( var1_nstep    (max_nvar) )
  allocate( var1_name     (max_nvar) )
  allocate( var1_desc     (max_nvar) )
  allocate( var1_unit     (max_nvar) )
  allocate( var1_layername(max_nvar) )
  allocate( var1_datatype (max_nvar) )
  allocate( var1_nlayer   (max_nvar) )
  allocate( var1_time_str (max_nvar) )
  allocate( var1_dt       (max_nvar) )

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

     call fio_register_file(ifid1(p),trim(infname))
     call fio_fopen(ifid1(p),FIO_FREAD)
     call fio_read_allinfo(ifid1(p))

     if ( p == 0 ) then ! only once
        allocate( hinfo%rgnid(1) )

        call fio_get_pkginfo(ifid1(p),hinfo)

        num_of_data = hinfo%num_of_data
        write(*,*)
        write(*,*) '*** get variable informations (file1)'
        write(*,*) 'num_of_data    : ', num_of_data

        nvar1 = 0
        do did = 0, num_of_data-1
           call fio_get_datainfo(ifid1(p),did,dinfo)

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

           do v = 1, nvar1
              if ( trim(var1_name(v)) == trim(dinfo%varname) ) then
                 var1_nstep(v) = var1_nstep(v) + 1
                 if( var1_nstep(v) == 2 ) var1_dt(v) = dinfo%time_start * 1.D-3 - var1_time_str(v)
                 if( var1_nstep(v) == step_str ) var1_time_str(v) = dinfo%time_start * 1.D-3
                 addvar = .false.
                 exit
              endif
           enddo

           if (addvar) then
              nvar1 = nvar1 + 1
              var1_nstep    (nvar1) = 1
              var1_name     (nvar1) = dinfo%varname
              var1_desc     (nvar1) = dinfo%description
              var1_unit     (nvar1) = dinfo%unit
              var1_layername(nvar1) = dinfo%layername
              var1_datatype (nvar1) = dinfo%datatype
              var1_nlayer   (nvar1) = dinfo%num_of_layer
              var1_time_str (nvar1) = dinfo%time_start * 1.D-3
              var1_dt       (nvar1) = ( dinfo%time_end - dinfo%time_start ) * 1.D-3
           endif

        enddo !--- did LOOP

        deallocate( hinfo%rgnid )
     endif !--- PE=000000
  enddo !--- PE LOOP

  if ( nvar1 == 0 ) then
     write(*,*) 'No variables to convert. Finish.'
     stop
  endif

  write(*,*) '########## Variable List (file1)', trim(infile(1)), ' ########## '
  write(*,*) 'ID |NAME            |STEPS|Layername       |START FROM                 |DT [sec]'
  do v = 1, nvar1
     call TIME_sec2date( date_str(:), real(var1_time_str(v),kind=8) )
     write(tmpl,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') &
               date_str(1),'/',date_str(2),'/',date_str(3),' ', &
               date_str(4),':',date_str(5),':',date_str(6)
     write(*,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A27,A1,I8)') &
              v,'|',var1_name(v),'|',var1_nstep(v),'|',var1_layername(v),'|', tmpl,'|', var1_dt(v)
  enddo

  ! file 2
  allocate( ifid2(0:PRC_nmax-1)      )
  allocate( var2_nstep    (max_nvar) )
  allocate( var2_name     (max_nvar) )
  allocate( var2_desc     (max_nvar) )
  allocate( var2_unit     (max_nvar) )
  allocate( var2_layername(max_nvar) )
  allocate( var2_datatype (max_nvar) )
  allocate( var2_nlayer   (max_nvar) )
  allocate( var2_time_str (max_nvar) )
  allocate( var2_dt       (max_nvar) )

  do p = 0, PRC_nmax-1
     call fio_mk_fname(infname,trim(infile(2)),'pe',p,6)

     call fio_put_commoninfo( FIO_MPIIO_NOUSE, &
                              FIO_SPLIT_FILE,  &
                              FIO_BIG_ENDIAN,  &
                              FIO_CARTESIAN,   &
                              int(GRID_DX),    &
                              GRID_IMAX,       &
                              1,               &
                              p                )

     call fio_register_file(ifid2(p),trim(infname))
     call fio_fopen(ifid2(p),FIO_FREAD)
     call fio_read_allinfo(ifid2(p))

     if ( p == 0 ) then ! only once
        allocate( hinfo%rgnid(1) )

        call fio_get_pkginfo(ifid2(p),hinfo)

        num_of_data = hinfo%num_of_data
        write(*,*)
        write(*,*) '*** get variable informations (file2)'
        write(*,*) 'num_of_data    : ', num_of_data

        nvar2 = 0
        do did = 0, num_of_data-1
           call fio_get_datainfo(ifid2(p),did,dinfo)

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

           do v = 1, nvar2
              if ( trim(var2_name(v)) == trim(dinfo%varname) ) then
                 var2_nstep(v) = var2_nstep(v) + 1
                 if( var2_nstep(v) == 2 ) var2_dt(v) = dinfo%time_start * 1.D-3 - var2_time_str(v)
                 if( var2_nstep(v) == step_str ) var2_time_str(v) = dinfo%time_start * 1.D-3
                 addvar = .false.
                 exit
              endif
           enddo

           if (addvar) then
              nvar2 = nvar2 + 1
              var2_nstep    (nvar2) = 1
              var2_name     (nvar2) = dinfo%varname
              var2_desc     (nvar2) = dinfo%description
              var2_unit     (nvar2) = dinfo%unit
              var2_layername(nvar2) = dinfo%layername
              var2_datatype (nvar2) = dinfo%datatype
              var2_nlayer   (nvar2) = dinfo%num_of_layer
              var2_time_str (nvar2) = dinfo%time_start * 1.D-3
              var2_dt       (nvar2) = ( dinfo%time_end - dinfo%time_start ) * 1.D-3
           endif

        enddo !--- did LOOP

        deallocate( hinfo%rgnid )
     endif !--- PE=000000
  enddo !--- PE LOOP

  if ( nvar2 == 0 ) then
     write(*,*) 'No variables to convert. Finish.'
     stop
  endif

  write(*,*) '########## Variable List (file2)', trim(infile(2)), ' ########## '
  write(*,*) 'ID |NAME            |STEPS|Layername       |START FROM                 |DT [sec]'
  do v = 1, nvar2
     call TIME_sec2date( date_str(:), real(var2_time_str(v),kind=8) )
     write(tmpl,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,F6.3)') &
               date_str(1),'/',date_str(2),'/',date_str(3),' ', &
               date_str(4),':',date_str(5),':',date_str(6)
     write(*,'(1x,I3,A1,A16,A1,I5,A1,A16,A1,A27,A1,I8)') &
              v,'|',var2_name(v),'|',var2_nstep(v),'|',var2_layername(v),'|', tmpl,'|', var2_dt(v)
  enddo

  !#########################################################

  allocate( data1_r4   (GRID_IMAX*GRID_JMAX*GRID_KMAX) )
  allocate( data2_r4   (GRID_IMAX*GRID_JMAX*GRID_KMAX) )
  allocate( data1_r4_3d(GRID_IMAX,GRID_JMAX,GRID_KMAX) )
  allocate( data2_r4_3d(GRID_IMAX,GRID_JMAX,GRID_KMAX) )

  allocate( data1_r8   (GRID_IMAX*GRID_JMAX*GRID_KMAX) )
  allocate( data2_r8   (GRID_IMAX*GRID_JMAX*GRID_KMAX) )
  allocate( data1_r8_3d(GRID_IMAX,GRID_JMAX,GRID_KMAX) )
  allocate( data2_r8_3d(GRID_IMAX,GRID_JMAX,GRID_KMAX) )

  write(*,*)

  do v1 = 1, nvar1
     write(*,*) '<file1> item:', trim(var1_name(v1)), ', step=', var1_nstep(v1)

     filematch = .false.
     do v2 = 1, nvar2
        if (       var1_name (v1) == var2_name (v2) &
             .AND. var1_nstep(v1) == var2_nstep(v2) ) then
           filematch = .true.
           exit
        endif
     enddo

     if ( filematch ) then
        write(*,*) '   is compared to <file2> item:', trim(var2_name(v2)), ', step=', var2_nstep(v2)

        write(*,*) 
        write(*,*) 'Check Header...'
        headercheck = .true.
        ! header comparison
        if ( var1_desc(v1) /= var2_desc(v2) ) then
           write(*,*) 'xxx mismatch! [Description]', var1_desc(v1), var2_desc(v2)
           headercheck = .false.
        endif
        if ( var1_unit(v1) /= var2_unit(v2) ) then
           write(*,*) 'xxx mismatch! [Unit]', var1_unit(v1), var2_unit(v2)
           headercheck = .false.
        endif
        if ( var1_layername(v1) /= var2_layername(v2) ) then
           write(*,*) 'xxx mismatch! [Layername]', var1_layername(v1), var2_layername(v2)
           headercheck = .false.
        endif
        if ( var1_datatype(v1) /= var2_datatype(v2) ) then
           write(*,*) 'xxx mismatch! [Datatype]', var1_datatype(v1), var2_datatype(v2)
           headercheck = .false.
        endif
        if ( var1_nlayer(v1) /= var2_nlayer(v2) ) then
           write(*,*) 'xxx mismatch! [Num. of Layer]', var1_nlayer(v1), var2_nlayer(v2)
           headercheck = .false.
        endif
        if ( var1_time_str(v1) /= var2_time_str(v2) ) then
           write(*,*) 'xxx mismatch! [Start time]', var1_time_str(v1), var2_time_str(v2)
           headercheck = .false.
        endif
        if ( var1_dt(v1) /= var2_dt(v2) ) then
           write(*,*) 'xxx mismatch! [Delta t]', var1_dt(v1), var2_dt(v2)
           headercheck = .false.
        endif
        if(headercheck) write(*,*) 'Check Header is OK.'

        write(*,*) 
        write(*,*) 'Check Data... criterion=',criterion
        datacheck = .true.
        ! data comparison
        do p = 0, PRC_nmax-1
           call fio_seek_datainfo(did1,ifid1(p),var1_name(v1),var1_nstep(v1))
           call fio_seek_datainfo(did2,ifid2(p),var2_name(v2),var2_nstep(v2))

           !--- read from pe000xx file
           if ( var1_datatype(v1) == FIO_REAL4 ) then
              call fio_read_data(ifid1(p),did1,data1_r4(:))
              data1_r4_3d(:,:,:) = reshape( data1_r4(:), shape(data1_r4_3d) )
              call fio_read_data(ifid2(p),did2,data2_r4(:))
              data2_r4_3d(:,:,:) = reshape( data2_r4(:), shape(data2_r4_3d) )

              do i = 1, GRID_IMAX
              do j = 1, GRID_JMAX
              do k = 1, GRID_KMAX
                 if ( abs( data2_r4_3d(i,j,k) - data1_r4_3d(i,j,k) ) > criterion ) then
                    write(*,*) "there is the diffence!:", p, i, j, k, data2_r4_3d(i,j,k)-data1_r4_3d(i,j,k)
                    datacheck = .false.
                 endif
              enddo
              enddo
              enddo
           elseif( var1_datatype(v1) == FIO_REAL8 ) then
              call fio_read_data(ifid1(p),did1,data1_r8(:))
              data1_r8_3d(:,:,:) = reshape( data1_r8(:), shape(data1_r8_3d) )
              call fio_read_data(ifid2(p),did2,data2_r8(:))
              data2_r8_3d(:,:,:) = reshape( data2_r8(:), shape(data2_r8_3d) )

              do i = 1, GRID_IMAX
              do j = 1, GRID_JMAX
              do k = 1, GRID_KMAX
                 if ( abs( data2_r8_3d(i,j,k) - data1_r8_3d(i,j,k) ) > criterion ) then
                    write(*,*) "there is the diffence!:", p, i, j, k, data2_r8_3d(i,j,k)-data1_r8_3d(i,j,k)
                    datacheck = .false.
                 endif
              enddo
              enddo
              enddo
           endif
        enddo
        if(datacheck) write(*,*) 'Check Data is OK.'
     else
        write(*,*) '   ...doesnt match any item in <file2>'
     endif
  enddo

  deallocate( data1_r4    )
  deallocate( data2_r4    )
  deallocate( data1_r4_3d )
  deallocate( data2_r4_3d )
  deallocate( data1_r8    )
  deallocate( data2_r8    )
  deallocate( data1_r8_3d )
  deallocate( data2_r8_3d )

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

end program prg_spddiff
!-------------------------------------------------------------------------------
