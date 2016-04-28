program convine
  implicit none

  include 'netcdf.inc'

  integer(4), parameter :: nst=1, nen=32  !--- time for average profile
  integer(4), parameter :: nt=nen-nst+1
  integer(4) :: nx, ny, nz, nxp, nyp
  integer(4) :: lz
  integer(4) :: nzhalo, nxhalo, nyhalo
  integer(4) :: xproc=3, yproc=3  !--- process number 
  real(8)    :: dt
  logical    :: ofirst = .true.
  logical    :: oread = .true.
  integer(4) :: start(3), start2(2), start4(4)
  integer(4) :: count(3), count2(2), count_land(3) ,count4(4)
  character*64 :: cfile
  character*64, allocatable :: vname(:)
  integer(4) :: it, ix, jy, kz, nrec, n
  integer(4) :: ncid, id01, status, iix, jjy, ndim
  integer(4) :: is, ie, js, je
  integer(4) :: ierr, vcount, prc_num_x, prc_num_y
  real(8),allocatable :: cz(:), cdz(:)
  real(8),allocatable :: cdx(:), cdy(:), cx(:), cy(:)
  real(8),allocatable :: p_cdx(:), p_cdy(:), p_cdz(:) 
  real(8),allocatable :: p_cx(:), p_cy(:), p_cz(:)
  real(4),allocatable :: var2d(:,:), var3d(:,:,:), var_land_3d(:,:,:)
  real(8),allocatable :: p_4d(:,:,:,:), p_3d(:,:,:), p_2d(:,:), p_land_3d(:,:,:)
  character*64 :: item
  character*5  :: HISTORY_DEFAULT_TUNIT
  real(8)      :: HISTORY_DEFAULT_TINTERVAL
  character*64 :: ATMOS_RESTART_OUT_BASENAME, rfile
  logical      :: ATMOS_RESTART_OUTPUT
 
  namelist  / PARAM_PRC / &
    PRC_NUM_X, PRC_NUM_Y
  namelist  / PARAM_HISTORY / &
    HISTORY_DEFAULT_TINTERVAL, &
    HISTORY_DEFAULT_TUNIT
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
  open(10,file="./init.d01.conf", form="formatted", access="sequential",iostat=ierr)
  !open(10,file="./pp.d01.conf", form="formatted", access="sequential",iostat=ierr)
  !open(10,file="./run.d01.conf", form="formatted", access="sequential",iostat=ierr)
  if( ierr /= 0 )then
   write(*,*) "Fail to open *.conf file"
   stop
  endif

  rewind(10)
  read(10,nml=PARAM_PRC,iostat=ierr)
  xproc = PRC_NUM_X
  yproc = PRC_NUM_Y

  do it = nst, nen !--- time
  !--- open NetCDF file and read from NetCDF file
   start(1:3) = (/1,1,it/)
   start2(1:2) = (/1,1/)
   start4(1:4) = (/1,1,1,it/)
   nrec = -1

  do n = 1, vcount
   nrec = -1
   do jy = 1, yproc
   do ix = 1, xproc
    nrec = nrec + 1
    !write(cfile,'(a,i6.6,a)') "./init_00000000000.000.pe",nrec,".nc"
    !write(cfile,'(a,i6.6,a)') "./landuse_d01.pe",nrec,".nc"
    !write(cfile,'(a,i6.6,a)') "./history_d01.pe",nrec,".nc"
     write(cfile,'(a,i6.6,a)') "./boundary_d01.pe",nrec,".nc"

    print *,trim(cfile)
    status = nf_open(cfile,0,ncid)
    if( status /= nf_noerr ) then
     write(*,*) "Stop at nf open"
     stop
    endif

    !if( ofirst ) then
    !  ofirst = .false.

      status = nf_inq_dimid( ncid, 'x', id01 )
      status = nf_inq_dimlen( ncid,id01, nxp )
      status = nf_inq_dimid( ncid, 'CX', id01 )
      status = nf_inq_dimlen( ncid,id01, nxhalo )
      nxhalo = ( nxhalo - nxp )/2

      status = nf_inq_dimid( ncid, 'y', id01 )
      status = nf_inq_dimlen( ncid,id01, nyp )
      status = nf_inq_dimid( ncid, 'CY', id01 )
      status = nf_inq_dimlen( ncid,id01, nyhalo )
      nyhalo = ( nyhalo - nyp )/2

      status = nf_inq_dimid( ncid, 'z', id01 )
      status = nf_inq_dimlen( ncid,id01, nz )
      status = nf_inq_dimid( ncid, 'CZ', id01 )
      status = nf_inq_dimlen( ncid,id01, nzhalo )
      nzhalo = ( nzhalo - nz )/2

      status = nf_inq_dimid( ncid, 'lz', id01 )
      status = nf_inq_dimlen( ncid,id01, lz )

!      count(1:4) = (/nxp,nyp,nz,1/)
!      count2(1:3) = (/nxp,nyp,1/)
!      allocate( p_3d(nxp,nyp,nz,1) )
!      allocate( p_2d(nxp,nyp,1) )
      count(1:3) = (/nz,nxp,nyp/)
      count2(1:2) = (/nxp,nyp/)
      count4(1:4) = (/nz,nxp,nyp,1/)
      count_land(1:3) = (/lz,nxp,nyp/)
      allocate( p_4d(nz,nxp,nyp,1) )
      allocate( p_3d(nz,nxp,nyp) )
      allocate( p_2d(nxp,nyp) )
      allocate( p_land_3d(lz,nxp,nyp) )
      allocate( p_cdx(nxp+2*nxhalo) )
      allocate( p_cdy(nyp+2*nyhalo) )
      allocate( p_cx(nxp+2*nxhalo) )
      allocate( p_cy(nyp+2*nyhalo) )

    if( ofirst ) then
      ofirst = .false.
      nx = (nxp-2) * xproc + 4
      ny = (nyp-2) * yproc + 4
      allocate( var3d(nx,ny,nz) )
      allocate( var2d(nx,ny) )
      allocate( var_land_3d(nx,ny,lz) )
      allocate( cdx(nx) )
      allocate( cdy(ny) )
      allocate( cx(nx) )
      allocate( cy(ny) )
      allocate( cz(nz) )
      allocate( cdz(nz) )
      allocate( p_cdz(nz+2*nzhalo) )
      allocate( p_cz(nz+2*nzhalo) )

      !--- z grid
      status = nf_inq_varid( ncid,'CZ',id01 )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf inq_varid cz"
       stop
      end if

      status = nf_get_var_double( ncid,id01,p_cz )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf get_var_double cz"
       stop
      end if

      status = nf_inq_varid( ncid,'CDZ',id01 )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf inq_varid cdz"
       stop
      end if

      status = nf_get_var_double( ncid,id01,p_cdz )
      if( status /= nf_noerr) then
       write(*,*) "stop at nf get_var_double cdz"
       stop
      end if

      do kz=1,nz
        cz(kz)=p_cz(kz+nzhalo)
        cdz(kz)=p_cdz(kz+nzhalo)
      enddo
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

    if( (trim(vname(n)) == "LAND_TEMP")    .or.  &
        (trim(vname(n)) == "LAND_WATER") ) then 
     status = nf_get_vara_double( ncid,id01,start,count_land,p_land_3d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double land ", trim(vname(n))
      stop
     end if
    elseif( ndim == 4 ) then
     status = nf_get_vara_double( ncid,id01,start4,count4,p_4d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double ", trim(vname(n))
      stop
     end if
    elseif( ndim == 3 ) then
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
    !do iix = (ix-1)*nxp+1, (ix-1)*nxp+nxp
    !  cdx(iix) = p_cdx(iix-(ix-1)*nxp+nxhalo)
    !  cx(iix) = p_cx(iix-(ix-1)*nxp+nxhalo)
    !enddo
    !do jjy = (jy-1)*nyp+1, (jy-1)*nyp+nyp
    !  cdy(jjy) = p_cdy(jjy-(jy-1)*nyp+nyhalo)
    !  cy(jjy) = p_cy(jjy-(jy-1)*nyp+nyhalo)
    !enddo

    if(ix==1)then
      is=1 ; ie=nxp
    else if(ix==xproc)then
      is=(ix-1)*(nxp-2)+2+1 ; ie=(ix-1)*(nxp-2)+2+nxp
    else
      is=(ix-1)*nxp+2+1 ; ie=(ix-1)*nxp+2+nxp
    endif
    if(jy==1)then
      js=1 ; je=nyp
    else if(jy==yproc)then
      js=(jy-1)*(nyp-2)+2+1 ; je=(jy-1)*(nyp-2)+2+nyp
    else
      js=(jy-1)*nyp+2+1 ; je=(jy-1)*nyp+2+nyp
    endif

    if(ix==1)then
     do iix = is, ie
       cdx(iix) = p_cdx(iix)
       cx(iix)  = p_cx(iix)
     enddo
    else if(ix==xproc)then
     do iix = is, ie
       cdx(iix) = p_cdx(iix-((ix-1)*(nxp-2)+2)+2)
       cx(iix)  = p_cx( iix-((ix-1)*(nxp-2)+2)+2)
     enddo
    else
     do iix = is, ie
        cdx(iix) = p_cdx(iix-((ix-1)*nxp+2)+nxhalo)
        cx(iix)  = p_cx( iix-((ix-1)*nxp+2)+nxhalo)
     enddo
    endif
    if(jy==1)then
     do jjy = js, je
       cdy(jjy) = p_cdy(jjy)
       cy(jjy)  = p_cy(jjy)
     enddo
    else if(jy==yproc)then
     do jjy = js, je
       cdy(jjy) = p_cdy(jjy-((jy-1)*(nyp-2)+2)+2)
       cy(jjy)  = p_cy( jjy-((jy-1)*(nyp-2)+2)+2)
     enddo
    else
     do jjy = js, je
       cdy(jjy) = p_cdy(jjy-((jy-1)*nyp+2)+nyhalo)
       cy(jjy)  = p_cy(jjy-((jy-1)*nyp+2)+nyhalo)
     enddo
    endif

    print *,ix,':',is,ie,jy,':',js,je
  
    if( (trim(vname(n)) == "LAND_TEMP")    .or.  &
        (trim(vname(n)) == "LAND_WATER") ) then 
     do iix = is,ie
     do jjy = js,je
         var_land_3d(iix,jjy,1:lz) = real(p_land_3d(1:lz,iix-is+1,jjy-js+1))
     enddo
     enddo
    elseif( ndim == 4 ) then
     do iix = is,ie
     do jjy = js,je
         var3d(iix,jjy,1:nz) = real(p_4d(1:nz,iix-is+1,jjy-js+1,1))
     enddo
     enddo
    elseif( ndim == 3 ) then
     do iix = is,ie
     do jjy = js,je
         var3d(iix,jjy,1:nz) = real(p_3d(1:nz,iix-is+1,jjy-js+1))
     enddo
     enddo
    elseif( ndim == 2 ) then
     do iix = is,ie
     do jjy = js,je
      var2d(iix,jjy) = real(p_2d(iix-is+1,jjy-js+1))
     enddo
     enddo
    endif

     deallocate( p_4d )
     deallocate( p_3d )
     deallocate( p_land_3d )
     deallocate( p_2d )
     deallocate( p_cdx )
     deallocate( p_cdy )
     deallocate( p_cx )
     deallocate( p_cy )

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
     if( (trim(vname(n)) == "LAND_TEMP")    .or.  &
         (trim(vname(n)) == "LAND_WATER") ) then 
      write(10,'(a,3x,i7,1x,a,1x,a)') "ZDEF", lz, "linear", "1 1"
     elseif( ndim == 3 .or. ndim == 4) then
      write(10,'(a,3x,i7,1x,a)') "ZDEF", nz, "LEVELS"
      write(10,'(5(1x,e15.7))') cz(1:nz)*1.d-3
     elseif( ndim == 2 ) then
      write(10,'(a,3x,i7,1x,a,1x,e15.7)') "ZDEF", 1, "LEVELS", cz(1)*1.d-3
     endif
     write(10,'(a,3x,i5,1x,a,1x,a,3x,a)') "TDEF", nen-nst+1, "LINEAR", "00:00Z01JAN2000", "1mn"
     write(10,'(a,3x,i2)') "VARS", 1
     if( (trim(vname(n)) == "LAND_TEMP")    .or.  &
         (trim(vname(n)) == "LAND_WATER") ) then 
      write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), lz, 99, "NONE"
     elseif( ndim == 3 .or. ndim == 4 ) then
      write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), nz, 99, "NONE"
     elseif( ndim == 2 ) then
      write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), 0, 99, "NONE"
     endif
     write(10,'(a)') "ENDVARS"
     close(10)
   endif

   if( (trim(vname(n)) == "LAND_TEMP")    .or.  &
       (trim(vname(n)) == "LAND_WATER") ) then 
    open(10,file=trim(vname(n))//".grd", &
         form="unformatted", access="direct", recl=4*nx*ny*lz)
    write(10,rec=it-nst+1) (((var_land_3d(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,lz)
    close(10)
   elseif( ndim == 3 .or. ndim == 4) then
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
