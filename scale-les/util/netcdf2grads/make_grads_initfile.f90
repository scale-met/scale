program convine
  implicit none

  include 'netcdf.inc'

  integer(4), parameter :: nst=1, nen=1  !--- time for average profile
  integer(4), parameter :: nt=nen-nst+1
  integer(4) :: nx, ny, nz, nxp, nyp
  integer(4) :: nzhalo, nxhalo, nyhalo
  integer(4) :: xproc=2, yproc=3  !--- process number 
  real(8) :: dt
  logical :: ofirst = .true.
  logical :: oread = .true.
  integer(4) :: start(4), start2(3)
  integer(4) :: count(4), count2(3)
  character*64 :: cfile
  character*64, allocatable :: vname(:)
  integer(4) :: it, ix, jy, kz, nrec, n
  integer(4) :: ncid, id01, status, iix, jjy, ndim
  integer(4) :: ierr, vcount, prc_num_x, prc_num_y
  real(8),allocatable :: cz(:), cdz(:)
  real(8),allocatable :: cdx(:), cdy(:), cx(:), cy(:)
  real(8),allocatable :: p_cdx(:), p_cdy(:)
  real(8),allocatable :: p_cx(:), p_cy(:)
  real(4),allocatable :: var2d(:,:), var3d(:,:,:)
!  real(8),allocatable :: p_3d(:,:,:,:), p_2d(:,:,:,:)
  real(8),allocatable :: p_3d(:,:,:), p_2d(:,:)
  character*64 :: item
  character*5 :: HISTORY_DEFAULT_TUNIT
  real(8) :: HISTORY_DEFAULT_TINTERVAL
  character*64 :: ATMOS_RESTART_OUT_BASENAME, rfile
  logical :: ATMOS_RESTART_OUTPUT
 
  namelist  / PARAM_PRC / &
    PRC_NUM_X, PRC_NUM_Y
  namelist  / PARAM_HISTORY / &
    HISTORY_DEFAULT_TINTERVAL, &
    HISTORY_DEFAULT_TUNIT
!  namelist  / HISTITEM / &
!    item
  namelist / PARAM_ATMOS_VARS / &
    ATMOS_RESTART_OUTPUT, &
    ATMOS_RESTART_OUT_BASENAME

  write(*,*) "Imput number of variable"
  read(*,*) vcount
  if( vcount == 0 ) then
   write(*,*) "Please determine the variable"
  elseif( vcount /= 0 ) then
   allocate(vname(vcount))
   do n = 1, vcount
    write(*,*) "Imput variable"
    read(*,*) vname(n)
    vname(n) = trim(vname(n))
   enddo
  endif

  !--- read variable calculated by model
  open(10,file="./init.conf", form="formatted", access="sequential",iostat=ierr)
  if( ierr /= 0 )then
   write(*,*) "Fail to open *.conf file"
   stop
  endif

  rewind(10)
  read(10,nml=PARAM_PRC,iostat=ierr)
  xproc = PRC_NUM_X
  yproc = PRC_NUM_Y

  rewind(10)
  read(10,nml=PARAM_ATMOS_VARS,iostat=ierr)
  rfile=trim(ATMOS_RESTART_OUT_BASENAME)

!  if( vcount == 0 ) then
!   rewind(10)
!   do n = 1, 500
!    read(10,nml=HISTITEM,iostat=ierr)
!    if( ierr == -1 ) then
!     vcount = n-1
!     allocate(vname(n-1))
!     exit
!    endif
!   enddo
!
!   rewind(10)
!   do n = 1, vcount
!    read(10,nml=HISTITEM,iostat=ierr)
!    vname(n) = trim(item)
!    vname(n) = trim(vname(n))
!    write(*,'(1x,i3,1x,a)') n, vname(n)
!   enddo
!   close(10)
!  endif

  do it = nst, nen !--- time
  !--- open NetCDF file and read from NetCDF file
   start(1:4) = (/1,1,1,it/)
   start2(1:3) = (/1,1,it/)
   nrec = -1

  do n = 1, vcount
   nrec = -1
   do jy = 1, yproc
   do ix = 1, xproc
    nrec = nrec + 1
    if( nrec < 10 ) then
     write(cfile,'(a,a,a,i1,a)') "./",trim(rfile),"_00000000000.000.pe00000",nrec,".nc"
    elseif( nrec < 100 ) then
     write(cfile,'(a,a,a,i2,a)') "./",trim(rfile),"_00000000000.000.pe0000",nrec,".nc"
    elseif( nrec < 1000 ) then
     write(cfile,'(a,a,a,i3,a)') "./",trim(rfile),"_00000000000.000.pe000",nrec,".nc"
    elseif( nrec < 10000 ) then
     write(cfile,'(a,a,a,i4,a)') "./",trim(rfile),"_00000000000.000.pe00",nrec,".nc"
    elseif( nrec < 100000 ) then
     write(cfile,'(a,a,a,i5,a)') "./",trim(rfile),"_00000000000.000.pe0",nrec,".nc"
    endif

    status = nf_open(cfile,0,ncid)
    if( status /= nf_noerr ) then
     write(*,*) "Stop at nf open"
     stop
    endif

    if( ofirst ) then
      ofirst = .false.

      status = nf_inq_dimid( ncid, 'x', id01 )
      status = nf_inq_dimlen( ncid,id01, nxp )
      nx = nxp * xproc
      status = nf_inq_dimid( ncid, 'CX', id01 )
      status = nf_inq_dimlen( ncid,id01, nxhalo )
      nxhalo = ( nxhalo - nxp )/2

      status = nf_inq_dimid( ncid, 'y', id01 )
      status = nf_inq_dimlen( ncid,id01, nyp )
      ny = nyp * yproc
      status = nf_inq_dimid( ncid, 'CY', id01 )
      status = nf_inq_dimlen( ncid,id01, nyhalo )
      nyhalo = ( nyhalo - nyp )/2

      status = nf_inq_dimid( ncid, 'z', id01 )
      status = nf_inq_dimlen( ncid,id01, nz )
      status = nf_inq_dimid( ncid, 'CZ', id01 )
      status = nf_inq_dimlen( ncid,id01, nzhalo )
      nzhalo = ( nzhalo - nz )/2

!      count(1:4) = (/nxp,nyp,nz,1/)
!      count2(1:3) = (/nxp,nyp,1/)
!      allocate( p_3d(nxp,nyp,nz,1) )
!      allocate( p_2d(nxp,nyp,1) )
      count(1:3) = (/nz,nxp,nyp/)
      count2(1:2) = (/nxp,nyp/)
      allocate( p_3d(nz,nxp,nyp) )
      allocate( p_2d(nxp,nyp) )
      allocate( var3d(nx,ny,nz) )
      allocate( var2d(nx,ny) )
      allocate( cz(nz+2*nzhalo) )
      allocate( cdz(nz+2*nzhalo) )
      allocate( cdx(nx) )
      allocate( cdy(ny) )
      allocate( cx(nx) )
      allocate( cy(ny) )
      allocate( p_cdx(nxp+2*nxhalo) )
      allocate( p_cdy(nyp+2*nyhalo) )
      allocate( p_cx(nxp+2*nxhalo) )
      allocate( p_cy(nyp+2*nyhalo) )

      !--- z grid
      status = nf_inq_varid( ncid,'CZ',id01 )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf inq_varid cz"
       stop
      end if

      status = nf_get_var_double( ncid,id01,cz )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf get_var_double cz"
       stop
      end if

      status = nf_inq_varid( ncid,'CDZ',id01 )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf inq_varid cdz"
       stop
      end if

      status = nf_get_var_double( ncid,id01,cdz )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf get_var_double cdz"
       stop
      end if
    endif

    !--- x grid
    status = nf_inq_varid( ncid,'CX',id01 )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf inq_varid cx"
     stop
    end if

    status = nf_get_var_double( ncid,id01,p_cx )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf get_var_double p_cx"
     stop
    end if

    status = nf_inq_varid( ncid,'CDX',id01 )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf inq_varid cdx"
    stop
    end if

    status = nf_get_var_double( ncid,id01,p_cdx )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf get_var_double cdx"
     stop
    end if

    !--- y grid
    status = nf_inq_varid( ncid,'CY',id01 )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf inq_varid cy"
    stop
    end if

    status = nf_get_var_double( ncid,id01,p_cy )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf get_var_double cy"
     stop
    end if

    status = nf_inq_varid( ncid,'CDY',id01 )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf inq_varid cdy"
    stop
    end if

    status = nf_get_var_double( ncid,id01,p_cdy )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf get_var_double cdy"
     stop
    end if

    !--- read variable
    status = nf_inq_varid( ncid,trim(vname(n)),id01 )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf inq_varid ", trim(vname(n))
     stop
    end if

    status = nf_inq_varndims( ncid,id01,ndim )
    if( status /= nf_noerr) then
     write(*,*) "stop at nf inq_varid dim", trim(vname(n))
     stop
    end if

    if( ndim == 3 ) then
     status = nf_get_vara_double( ncid,id01,start,count,p_3d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double ", trim(vname(n))
      stop
     end if
    elseif( ndim == 2 ) then
     status = nf_get_vara_double( ncid,id01,start2,count2,p_2d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double ", trim(vname(n))
      stop
     end if
    end if

    status = nf_close(ncid)

    !--- conbine variables in each processor to single
    do iix = (ix-1)*nxp+1, (ix-1)*nxp+nxp
      cdx(iix) = p_cdx(iix-(ix-1)*nxp+nxhalo)
      cx(iix) = p_cx(iix-(ix-1)*nxp+nxhalo)
    enddo
    do jjy = (jy-1)*nyp+1, (jy-1)*nyp+nyp
      cdy(jjy) = p_cdy(jjy-(jy-1)*nyp+nyhalo)
      cy(jjy) = p_cy(jjy-(jy-1)*nyp+nyhalo)
    enddo

    if( ndim == 3 ) then
     do iix = (ix-1)*nxp+1, (ix-1)*nxp+nxp
     do jjy = (jy-1)*nyp+1, (jy-1)*nyp+nyp
     do kz = 1, nz
      var3d(iix,jjy,kz) = real(p_3d(kz,iix-(ix-1)*nxp,jjy-(jy-1)*nyp))
     enddo
     enddo
     enddo
    elseif( ndim == 2 ) then
     do iix = (ix-1)*nxp+1, (ix-1)*nxp+nxp
     do jjy = (jy-1)*nyp+1, (jy-1)*nyp+nyp
      var2d(iix,jjy) = real(p_2d(iix-(ix-1)*nxp,jjy-(jy-1)*nyp))
     enddo
     enddo
    endif

   enddo
   enddo

   if( it == nst ) then
     write(*,*) "create ctl file of ", trim(vname(n))
     open(10,file=trim(vname(n))//".ctl", form="formatted", access="sequential" )
     write(10,'(a,1x,a)') "DSET", "^"//trim(vname(n))//".grd"
     write(10,'(a)') "TITLE SCALE-LES data output"
     write(10,'(a)') "OPTIONS BIG_ENDIAN"
     write(10,'(a,1x,e15.7)') "UNDEF", -0.99900E+35
     write(10,'(a,3x,i7,1x,a)') "XDEF", nx, "LEVELS"
     write(10,'(5(1x,e15.7))') cx(1:nx)*1.d-3
     write(10,'(a,3x,i7,1x,a)') "YDEF", ny, "LEVELS"
     write(10,'(5(1x,e15.7))') cy(1:ny)*1.d-3
     if( ndim == 3 ) then
      write(10,'(a,3x,i7,1x,a)') "ZDEF", nz, "LEVELS"
      write(10,'(5(1x,e15.7))') cz(nzhalo+1:nz+nzhalo)*1.d-3
     elseif( ndim == 2 ) then
      write(10,'(a,3x,i7,1x,a,1x,e15.7)') "ZDEF", 1, "LEVELS", cz(nzhalo+1)*1.d-3
     endif
     write(10,'(a,3x,i5,1x,a,1x,a,3x,a)') "TDEF", nen-nst+1, "LINEAR", "00:00Z01JAN2000", "1mn"
     write(10,'(a,3x,i2)') "VARS", 1
     if( ndim == 3 ) then
      write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), nz, 99, "NONE"
     elseif( ndim == 2 ) then
      write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), 0, 99, "NONE"
     endif
     write(10,'(a)') "ENDVARS"
     close(10)
   endif

   if( ndim == 3 ) then
    open(10,file=trim(vname(n))//".grd", &
         form="unformatted", access="direct", recl=4*nx*ny*nz)
    write(10,rec=it-nst+1) (((var3d(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
    close(10)
   elseif( ndim == 2 ) then
    open(10,file=trim(vname(n))//".grd", &
         form="unformatted", access="direct", recl=4*nx*ny)
    write(10,rec=it-nst+1) ((var2d(ix,jy),ix=1,nx),jy=1,ny)
    close(10)
   endif
  enddo !--- for variable (n)
  enddo !--- for time (it)

  write(*,*) "CALCULATION FINISHED"

end program
