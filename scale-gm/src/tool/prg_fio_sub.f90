!-------------------------------------------------------------------------------
!> Program FIO sub
!!
!! @par Description
!!          substitute pe0000x format data
!!          ( packaged NICAM data format : PaNDa )
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
program fio_sub
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use mod_io_param
  use scale_stdio
  use mod_fio, only: &
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

  integer, parameter :: flim = 2
  integer,      save :: fmax

  !--- NAMELIST
  integer               :: glevel       = -1
  integer               :: rlevel       = -1
  character(len=H_LONG) :: mnginfo      = ""
  character(len=H_LONG) :: infile(flim) = ""
  character(len=H_LONG) :: outfile      = ""
  logical               :: use_mpi      = .true.
  integer               :: pe_str       =  0
  integer               :: pe_end       = -1
  logical               :: help         = .false.

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
  character(len=H_LONG)  :: infname1 = ""
  character(len=H_LONG)  :: infname2 = ""
  character(len=H_LONG)  :: outfname = ""

  type(headerinfo)       :: hinfo
  type(datainfo)         :: dinfo1
  type(datainfo)         :: dinfo2

  character(len=H_MID)   :: pkg_desc
  character(len=H_LONG)  :: pkg_note
  integer                :: nmax_data

  integer                :: nvar
  character(len=H_SHORT) :: var_name (max_nvar)
  integer                :: var_nstep(max_nvar)

  integer                :: GALL
  integer                :: KALL
  integer                :: LALL
  real(SP), allocatable  :: data4_1D1(:)
  real(DP), allocatable  :: data8_1D1(:)
  real(SP), allocatable  :: data4_1D2(:)
  real(DP), allocatable  :: data8_1D2(:)

  ! for MPI
  integer                :: pe_all
  integer                :: prc_nall, prc_nlocal
  integer                :: prc_myrank, pstr, pend
  character(len=6)       :: rankstr

  logical :: addvar
  integer :: p, v, vid
  integer :: ifid1, ifid2, idid, ofid, odid, ierr
  !=====================================================================

  !--- read option and preprocess
  call readoption !! set fmax, infile

  !--- prepare region infomation
  call MNG_mnginfo_input( rlevel, trim(mnginfo) )

  GALL = ( (2**(glevel-rlevel))+2 ) &
       * ( (2**(glevel-rlevel))+2 )

  IO_FID_LOG = IO_get_available_fid()
  if ( use_mpi ) then
     !--- Parallel Excution, No communication
     call MPI_Init(ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, prc_nall,   ierr)
     call MPI_Comm_rank(MPI_COMM_WORLD, prc_myrank, ierr)
     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     write(rankstr,'(I6.6)') prc_myrank
     open(IO_FID_LOG, file='LOG_sub.pe'//trim(rankstr) )
     write(IO_FID_LOG,*) "+++ Parallel Execution, Use MPI"

     if( mod( MNG_PALL, prc_nall) /= 0)then
        write(IO_FID_LOG,*) "*** Invalid processor number, STOP:", MNG_PALL, prc_nall
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
     endif
  else
     open(IO_FID_LOG, file='msg.serial' )
     write(IO_FID_LOG,*) "+++ Serial Execution"
     prc_nall   = 1
     prc_myrank = 0
  endif

  if ( pe_end >= 0 ) then
     pe_all = pe_end - pe_str + 1
     write(IO_FID_LOG,*) "*** pe range is specified. "
     write(IO_FID_LOG,*) "*** pe(all,start,end)=",pe_all,pe_str,pe_end
  else
     pe_all = MNG_PALL
  endif

  if ( mod( pe_all, prc_nall) /= 0 ) then
     write(IO_FID_LOG,*) "*** Invalid processor number, STOP:", pe_all, prc_nall
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
     stop
  endif

  prc_nlocal = pe_all / prc_nall
  pstr       = prc_myrank*prc_nlocal + pe_str + 1
  pend       = prc_myrank*prc_nlocal + pe_str + prc_nlocal
  write(IO_FID_LOG,*) "*** Number of Total .pexxxxxx files: ", MNG_PALL
  write(IO_FID_LOG,*) "*** Number of PE to packing precess: ", prc_nall
  write(IO_FID_LOG,*) "*** The rank of this process       : ", prc_myrank
  write(IO_FID_LOG,*) "*** Number of files for this rank  : ", prc_nlocal
  write(IO_FID_LOG,*) "*** file ID to pack                : ", pstr-1, " - ", pend-1

  !--- setup
  call fio_syscheck()
  write(IO_FID_LOG,*) '*** combine start : PaNDa format to PaNDa format data'

  do p = pstr, pend
     write(IO_FID_LOG,*) '+pe:', p-1
     LALL = MNG_prc_rnum(p)

     call fio_mk_fname(infname1,trim(infile(1)),'pe',p-1,6)
     call fio_mk_fname(infname2,trim(infile(2)),'pe',p-1,6)
     call fio_mk_fname(outfname,trim(outfile),  'pe',p-1,6)
     write(IO_FID_LOG,*) '++output : ', trim(outfname)

     call fio_register_file(ifid1,trim(infname1))
     call fio_fopen(ifid1,IO_FREAD)
     ! put information from 1st input file
     call fio_put_commoninfo_fromfile(ifid1,IO_BIG_ENDIAN)

     call fio_read_allinfo(ifid1)
     allocate( hinfo%rgnid(LALL) )
     call fio_get_pkginfo(ifid1,hinfo)
     pkg_desc  = hinfo%description
     pkg_note  = hinfo%note
     nmax_data = hinfo%num_of_data
     write(IO_FID_LOG,*) '++input', 1, ' : ', trim(infname1), "(n=", hinfo%num_of_data, ")"

     call fio_register_file(ifid2,trim(infname2))
     call fio_fopen(ifid2,IO_FREAD)
     call fio_read_allinfo(ifid2)
     call fio_get_pkginfo(ifid2,hinfo)
     write(IO_FID_LOG,*) '++input', 2, ' : ', trim(infname2), "(n=", hinfo%num_of_data, ")"

     if ( hinfo%num_of_data /= nmax_data ) then
        write(IO_FID_LOG,*) "*** Mismatch number of data, STOP:", hinfo%num_of_data, nmax_data
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
     endif

     call fio_register_file(ofid,trim(outfname))
     call fio_fopen(ofid,IO_FWRITE)
     call fio_put_write_pkginfo(ofid,pkg_desc,pkg_note)

     nvar = 0
     do idid = 0, nmax_data-1
        ! get datainfo from input file
        call fio_get_datainfo(ifid1,idid,dinfo1)
        call fio_get_datainfo(ifid2,idid,dinfo2)
        KALL = dinfo1%num_of_layer

        if ( dinfo2%num_of_layer /= KALL ) then
           write(IO_FID_LOG,*) "*** Mismatch number of layer, STOP:", dinfo2%num_of_layer, KALL
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           call MPI_FINALIZE(ierr)
           stop
        endif

        addvar = .true.
        do v = 1, nvar
           if ( var_name(v) == dinfo1%varname ) then
              var_nstep(v) = var_nstep(v) + 1
              vid = v
              addvar = .false.
              exit
           endif
        enddo

        if (addvar) then
           nvar = nvar + 1
           var_nstep(nvar) = 1
           var_name (nvar) = dinfo1%varname
           vid = nvar
        endif

        dinfo1%step = var_nstep(vid)

        ! read->write data
        if ( dinfo1%datatype == IO_REAL4 ) then
           allocate( data4_1D1(GALL*KALL*LALL) )
           call fio_read_data(ifid1,idid,data4_1D1)
           if ( dinfo2%datatype == IO_REAL4 ) then
              allocate( data4_1D2(GALL*KALL*LALL) )
              call fio_read_data(ifid2,idid,data4_1D2)
              data4_1D1(:) = data4_1D1(:) - data4_1D2(:)
              deallocate( data4_1D2 )
           elseif( dinfo2%datatype == IO_REAL8 ) then
              allocate( data8_1D2(GALL*KALL*LALL) )
              call fio_read_data(ifid2,idid,data8_1D2)
              data4_1D1(:) = data4_1D1(:) - real(data8_1D2(:),kind=4)
              deallocate( data8_1D2 )
           endif
           call fio_put_write_datainfo_data(odid,ofid,dinfo1,data4_1D1)
           deallocate( data4_1D1 )
        elseif( dinfo1%datatype == IO_REAL8 ) then
           allocate( data8_1D1(GALL*KALL*LALL) )
           call fio_read_data(ifid1,idid,data8_1D1)
           if ( dinfo2%datatype == IO_REAL4 ) then
              allocate( data4_1D2(GALL*KALL*LALL) )
              call fio_read_data(ifid2,idid,data4_1D2)
              data8_1D1(:) = data8_1D1(:) - real(data4_1D2(:),kind=8)
              deallocate( data4_1D2 )
           elseif( dinfo2%datatype == IO_REAL8 ) then
              allocate( data8_1D2(GALL*KALL*LALL) )
              call fio_read_data(ifid2,idid,data8_1D2)
              data8_1D1(:) = data8_1D1(:) - data8_1D2(:)
              deallocate( data8_1D2 )
           endif
           call fio_put_write_datainfo_data(odid,ofid,dinfo1,data8_1D1)
           deallocate( data8_1D1 )
        endif

     enddo

     call fio_fclose(ifid1)
     call fio_fclose(ifid2)
     call fio_fclose(ofid)

     deallocate( hinfo%rgnid )
  enddo ! PE loop

  if ( use_mpi ) then
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     call MPI_FINALIZE(ierr)
  endif

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

end program fio_sub
