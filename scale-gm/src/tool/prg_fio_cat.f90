!-------------------------------------------------------------------------------
!
!+  Program FIO cat
!
!-------------------------------------------------------------------------------
program fio_cat
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !      combine pe0000x format data
  !      ( packaged NICAM data format : PaNDa )
  !
  !++ Current Corresponding Author: H. Yashiro
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      1.00      13-09-30  H.Yashiro : [NEW]
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use mod_misc, only: &
     MISC_get_available_fid
  use mod_fio, only : &
     FIO_HSHORT,      &
     FIO_HMID,        &
     FIO_HLONG,       &
     FIO_REAL4,       &
     FIO_REAL8,       &
     FIO_BIG_ENDIAN,  &
     FIO_FREAD,       &
     FIO_FWRITE,      &
     headerinfo, &
     datainfo
  use mod_mnginfo_light, only : &
     MNG_mnginfo_input,   &
     MNG_PALL,            &
     MNG_prc_rnum
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ param & variable
  !
  integer, parameter :: max_nvar = 1000

  integer, parameter :: flim = 100
  integer,      save :: fmax

  !--- NAMELIST
  integer                   :: glevel       = -1
  integer                   :: rlevel       = -1
  character(LEN=FIO_HLONG)  :: mnginfo      = ""
  character(LEN=FIO_HLONG)  :: infile(flim) = ""
  character(LEN=FIO_HLONG)  :: outfile      = ""
  logical                   :: use_mpi      = .true.
  integer                   :: pe_str       =  0
  integer                   :: pe_end       = -1
  logical                   :: help         = .false.

  namelist /OPTION/ glevel,  &
                    rlevel,  &
                    mnginfo, &
                    infile,  &
                    outfile, &
                    use_mpi, &
                    pe_str,  &
                    pe_end,  &
                    help

  !-----------------------------------------------------------------------------
  character(LEN=FIO_HLONG)  :: infname  = ""
  character(LEN=FIO_HLONG)  :: outfname = ""

  type(headerinfo)          :: hinfo
  type(datainfo)            :: dinfo

  character(LEN=FIO_HMID)   :: pkg_desc
  character(LEN=FIO_HLONG)  :: pkg_note
  integer                   :: nmax_data

  integer                   :: nvar
  character(LEN=FIO_HSHORT) :: var_name (max_nvar)
  integer                   :: var_nstep(max_nvar)

  integer              :: GALL
  integer              :: KALL
  integer              :: LALL
  real(4), allocatable :: data4_1D(:)
  real(8), allocatable :: data8_1D(:)

  ! for MPI
  integer              :: pe_all
  integer              :: prc_nall, prc_nlocal
  integer              :: prc_myrank, pstr, pend
  integer              :: fid_log
  character(LEN=6)     :: rankstr

  logical :: addvar
  integer :: p, f, v, vid
  integer :: ifid, idid, ofid, odid, ierr
  !=====================================================================

  !--- read option and preprocess
  call readoption !! set fmax, infile

  !--- prepare region infomation
  call MNG_mnginfo_input( rlevel, trim(mnginfo) )

  GALL = ( (2**(glevel-rlevel))+2 ) &
       * ( (2**(glevel-rlevel))+2 )

  fid_log = MISC_get_available_fid()
  if ( use_mpi ) then
     !--- Parallel Excution, No communication
     call MPI_Init(ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, prc_nall,   ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD, prc_myrank, ierr)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     write(rankstr,'(I6.6)') prc_myrank
     open(fid_log, file='msg.pe'//trim(rankstr) )
     write(fid_log,*) "+++ Parallel Execution, Use MPI"

     if( mod( MNG_PALL, prc_nall) /= 0)then
        write(fid_log,*) "*** Invalid processor number, STOP:", MNG_PALL, prc_nall
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
     endif
  else
     open(fid_log, file='msg.serial' )
     write(fid_log,*) "+++ Serial Execution"
     prc_nall   = 1
     prc_myrank = 0
  endif

  if ( pe_end >= 0 ) then
     pe_all = pe_end - pe_str + 1
     write(fid_log,*) "*** pe range is specified. "
     write(fid_log,*) "*** pe(all,start,end)=",pe_all,pe_str,pe_end
  else
     pe_all = MNG_PALL
  endif

  if ( mod( pe_all, prc_nall) /= 0 ) then
     write(fid_log,*) "*** Invalid processor number, STOP:", pe_all, prc_nall
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
     stop
  endif

  prc_nlocal = pe_all / prc_nall
  pstr       = prc_myrank*prc_nlocal + pe_str + 1
  pend       = prc_myrank*prc_nlocal + pe_str + prc_nlocal
  write(fid_log,*) "*** Number of Total .pexxxxxx files: ", MNG_PALL
  write(fid_log,*) "*** Number of PE to packing precess: ", prc_nall
  write(fid_log,*) "*** The rank of this process       : ", prc_myrank
  write(fid_log,*) "*** Number of files for this rank  : ", prc_nlocal
  write(fid_log,*) "*** file ID to pack                : ", pstr-1, " - ", pend-1

  !--- setup
  call fio_syscheck()

  write(fid_log,*) '*** combine start : PaNDa format to PaNDa format data'

  do p = pstr, pend
     write(fid_log,*) '+pe:', p-1
     LALL = MNG_prc_rnum(p)

     call fio_mk_fname(infname, trim(infile(1)),'pe',p-1,6)
     call fio_mk_fname(outfname,trim(outfile),  'pe',p-1,6)
     write(fid_log,*) '++output : ', trim(outfname)

     call fio_register_file(ifid,trim(infname))
     call fio_fopen(ifid,FIO_FREAD)
     ! put information from 1st input file
     call fio_put_commoninfo_fromfile(ifid,FIO_BIG_ENDIAN)

     call fio_read_allinfo(ifid)
     allocate( hinfo%rgnid(LALL) )
     call fio_get_pkginfo(ifid,hinfo)
     pkg_desc  = hinfo%description
     pkg_note  = hinfo%note
     nmax_data = hinfo%num_of_data
     write(fid_log,*) '++input', 1, ' : ', trim(infname), "(n=", nmax_data, ")"

     call fio_register_file(ofid,trim(outfname))
     call fio_fopen(ofid,FIO_FWRITE)
     call fio_put_write_pkginfo(ofid,pkg_desc,pkg_note)

     nvar = 0
     do idid = 0, nmax_data-1
        ! get datainfo from input file
        call fio_get_datainfo(ifid,idid,dinfo)
        KALL = dinfo%num_of_layer

        addvar = .true.
        do v = 1, nvar
           if ( var_name(v) == dinfo%varname ) then
              var_nstep(v) = var_nstep(v) + 1
              vid = v
              addvar = .false.
              exit
           endif
        enddo

        if (addvar) then
           nvar = nvar + 1
           var_nstep(nvar) = 1
           var_name (nvar) = dinfo%varname
           vid = nvar
        endif

        dinfo%step = var_nstep(vid)

        ! read->write data
        if ( dinfo%datatype == FIO_REAL4 ) then
           allocate( data4_1D(GALL*KALL*LALL) )
           call fio_read_data(ifid,idid,data4_1D)
           call fio_put_write_datainfo_data(odid,ofid,dinfo,data4_1D)
           deallocate( data4_1D )
        elseif( dinfo%datatype == FIO_REAL8 ) then
           allocate( data8_1D(GALL*KALL*LALL) )
           call fio_read_data(ifid,idid,data8_1D)
           call fio_put_write_datainfo_data(odid,ofid,dinfo,data8_1D)
           deallocate( data8_1D )
        endif
     enddo

     call fio_fclose(ifid)

     do f = 2, fmax
        call fio_mk_fname(infname, trim(infile(f)),'pe',p-1,6)

        call fio_register_file(ifid,trim(infname))
        call fio_fopen(ifid,FIO_FREAD)
        call fio_read_allinfo(ifid)
        call fio_get_pkginfo(ifid,hinfo)
        nmax_data = hinfo%num_of_data
        write(fid_log,*) '++input', f, ' : ', trim(infname), "(n=", nmax_data, ")"

        do idid = 0, nmax_data-1
           ! get datainfo from input file
           call fio_get_datainfo(ifid,idid,dinfo)
           KALL = dinfo%num_of_layer

           addvar = .true.
           do v = 1, nvar
              if ( var_name(v) == dinfo%varname ) then
                 var_nstep(v) = var_nstep(v) + 1
                 addvar = .false.
                 exit
              endif
           enddo

           if (addvar) then
              nvar = nvar + 1
              var_nstep(nvar) = 1
              var_name (nvar) = dinfo%varname
           endif

           dinfo%step = var_nstep(v)

           ! read->write data
           if ( dinfo%datatype == FIO_REAL4 ) then
              allocate( data4_1D(GALL*KALL*LALL) )
              call fio_read_data(ifid,idid,data4_1D)
              call fio_put_write_datainfo_data(odid,ofid,dinfo,data4_1D)
              deallocate( data4_1D )
           elseif( dinfo%datatype == FIO_REAL8 ) then
              allocate( data8_1D(GALL*KALL*LALL) )
              call fio_read_data(ifid,idid,data8_1D)
              call fio_put_write_datainfo_data(odid,ofid,dinfo,data8_1D)
              deallocate( data8_1D )
           endif
        enddo

        call fio_fclose(ifid)
     enddo

     call fio_fclose(ofid)

     deallocate( hinfo%rgnid )
  enddo ! PE loop

  if ( use_mpi ) then
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
  endif

  close(fid_log)

contains
  !-----------------------------------------------------------------------------
  !> read option
  !-----------------------------------------------------------------------------
  subroutine readoption
    use mod_misc, only : &
      MISC_get_available_fid
    use mod_tool_option, only: &
      OPT_convert, &
      OPT_fid
    implicit none

    integer :: io
    !---------------------------------------------------------------------------

    ! --- Set option
    OPT_fid = MISC_get_available_fid()
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

end program fio_cat
!-------------------------------------------------------------------------------
