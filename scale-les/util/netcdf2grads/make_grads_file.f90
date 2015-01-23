program convine
  !---------------------------------------------------------
  ! [for hurricane]
  !ifort -O3 -xHOST -convert big_endian -assume byterecl -I/ap/netcdf4/4.1.3/include -L/ap/netcdf4/4.1.3/lib -lnetcdff -lnetcdf make_grads_file.f90 -o convine
  !---------------------------------------------------------
  ! [for marco]
  !ifort -O3 -xHOST -convert big_endian -assume byterecl -I/ap/netcdf4/4.1.3-nodap/include -L/ap/netcdf4/4.1.3-nodap/lib -lnetcdff -lnetcdf make_grads_file.f90 -o convine
  !---------------------------------------------------------
  ! [for pearlscale]
  !ifort -O3 -xHOST -convert big_endian -assume byterecl -I/opt/ap/netcdf4/include -L/opt/ap/netcdf4/lib -lnetcdff -lnetcdf make_grads_file.f90 -o convine
  !---------------------------------------------------------
  !
  implicit none

  include 'netcdf.inc'
 
  character*2 :: ARG           ! s:silent , i:interactive

  integer(4) :: nst, nen       !--- time for average profile
  integer(4) :: nt
  integer(4) :: vcount
  character*64, allocatable :: vname(:)

  ! for silent
  integer(4)                :: uz,lz
  integer                   :: timestep
  integer                   :: inest
  character*2               :: cnest
  character*100             :: idir,odir,conffile
  character*5               :: delt
  character*15              :: stime
  namelist /info/  timestep,inest,conffile,idir,odir,vcount
  namelist /vari/  vname
  namelist /grads/ delt,stime

  ! for interactive
  data delt  /"1mn"/
  data stime /"00:00Z01JAN2000"/
  data idir  /"./"/
  data odir  /"./"/
  data conffile /"./run.conf"/

  integer(4) :: xproc=2, yproc=3         !--- process number 
  integer(4) :: nx, ny, nz, nxp, nyp
  integer(4) :: nzhalo, nxhalo, nyhalo

  character*100 :: cfile
  integer(4)    :: nrec

  real(8)    :: dt
  logical    :: ofirst = .true.
  logical    :: oread = .true.
  integer(4) :: start(4), start2(3)
  integer(4) :: count(4), count2(3), count_urban(4), count_land(4), count_height(3)

  integer(4) :: it, ix, jy, kz, n
  integer(4) :: ncid, id01, status, iix, jjy, ndim
  integer(4) :: is,ie,js,je
  integer(4) :: ierr, prc_num_x, prc_num_y

  real(8),allocatable :: cz(:), cdz(:)
  real(8),allocatable :: cdx(:), cdy(:), cx(:), cy(:)
  real(8),allocatable :: p_cdx(:), p_cdy(:)
  real(8),allocatable :: p_cx(:), p_cy(:)
  real(4),allocatable :: var2d(:,:), var3d(:,:,:), var_urban_3d(:,:,:), var_land_3d(:,:,:)
  real(8),allocatable :: p_3d(:,:,:,:), p_2d(:,:,:), p_urban_3d(:,:,:,:), p_land_3d(:,:,:,:), p_height(:,:,:)

  character*64 :: HISTORY_DEFAULT_BASENAME
  character*5  :: HISTORY_DEFAULT_TUNIT
  real(8)      :: HISTORY_DEFAULT_TINTERVAL
  logical      :: HISTORY_DEFAULT_ZINTERP = .true.
  character*64 :: item
  logical      :: HIST_BND = .true.

  namelist  / PARAM_PRC / &
    PRC_NUM_X, PRC_NUM_Y
  namelist  / PARAM_HISTORY /  &
    HISTORY_DEFAULT_BASENAME,  &
    HISTORY_DEFAULT_TINTERVAL, &
    HISTORY_DEFAULT_TUNIT,     &
    HISTORY_DEFAULT_ZINTERP
  namelist  /PARAM_HIST/ &
    HIST_BND
  namelist  / HISTITEM / &
    item

!-----------------------------------------------
!  select option
!-----------------------------------------------

  call GETARG(1,ARG)
   
!-----------------------------------------------
!  read information
!-----------------------------------------------

  select case(ARG)
  case("  ") ! empty
     write(*,*) "please specify the option '-i' or '-s'"
     stop
  case("-i") ! interactive
     write(*,*) "path to configure file for run with the quotation mark"
     read(*,*) conffile
     write(*,*) "path to directory of history files with the quotation mark"
     read(*,*) idir
     write(*,*) "path to directory of output files with the quotation mark"
     read(*,*) odir

     write(*,*) "start time of convert data"
     read(*,*) nst
     write(*,*) "end time of convert data"
     read(*,*) nen
     nt=nen-nst+1
     timestep=nt

     write(*,*) "Imput number of variable"
     write(*,*) "0 -> all variable output from model"
     read(*,*) vcount
     if( vcount /= 0 ) then
        allocate(vname(vcount))
        do n = 1, vcount
           write(*,*) "Imput variable"
           read(*,*) vname(n)
           vname(n) = trim(vname(n))
        enddo
     endif

     write(6,nml=info)
     write(6,nml=vari)
     write(6,nml=grads)

  case("-s") ! silent
     open(11,file='namelist.in',status='old')
     read(11,nml=info)
     write(6,nml=info)
     nst=1
     nen=timestep
     nt=nen-nst+1
     if( vcount /= 0 ) then
        allocate(vname(vcount))
     endif
     read(11,nml=vari)
     write(6,nml=vari)
     do n = 1, vcount
        vname(n) = trim(vname(n)) 
     enddo
     read(11,nml=grads)
     write(6,nml=grads)
     close(11)
  end select


  !--- read variable calculated by model
  open(10,file=trim(conffile), form="formatted", access="sequential",iostat=ierr)
  if( ierr /= 0 )then
   write(*,*) "Fail to open *.conf file"
   stop
  endif

  rewind(10)
  read(10,nml=PARAM_PRC,iostat=ierr)
  write(6,nml=PARAM_PRC)
  xproc = PRC_NUM_X
  yproc = PRC_NUM_Y

  rewind(10)
  read(10,nml=PARAM_HIST,iostat=ierr)
  write(6,nml=PARAM_HIST)

  rewind(10)
  read(10,nml=PARAM_HISTORY,iostat=ierr)
  write(6,nml=PARAM_HISTORY)

  if( vcount == 0 ) then
   rewind(10)
   do n = 1, 500
    read(10,nml=HISTITEM,iostat=ierr)
    if( ierr == -1 ) then
     vcount = n-1
     allocate(vname(n-1))
     exit
    endif
   enddo

   rewind(10)
   do n = 1, vcount
    read(10,nml=HISTITEM,iostat=ierr)
    vname(n) = trim(item)
    vname(n) = trim(vname(n))
    write(*,'(1x,i3,1x,a)') n, vname(n)
   enddo
  endif
  close(10)

  do it = nst, nen !--- time
  !--- open NetCDF file and read from NetCDF file
   start(1:4) = (/1,1,1,it/)
   start2(1:3) = (/1,1,it/)
   nrec = -1
   write(*,*) "TIME= ", it

 do n = 1, vcount
   nrec = -1
   do jy = 1, yproc
   do ix = 1, xproc
    nrec = nrec + 1
    write(cfile,'(a,i6.6,a)') trim(HISTORY_DEFAULT_BASENAME)//".pe",nrec,".nc"
  
    cfile=trim(idir)//'/'//trim(cfile)
    !write(6,*) trim(cfile)

    status = nf_open(cfile,0,ncid)
    if( status /= nf_noerr ) then
     write(*,*) "Stop at nf open"
     stop
    endif
 
    !if( ofirst ) then
    !  ofirst = .false.

      status = nf_inq_dimid ( ncid,  'x', id01 )
      status = nf_inq_dimlen( ncid, id01, nxp )
      status = nf_inq_dimid ( ncid, 'CX', id01 )
      status = nf_inq_dimlen( ncid, id01, nxhalo )
      nxhalo = ( nxhalo - nxp )/2

      status = nf_inq_dimid ( ncid,  'y', id01 )
      status = nf_inq_dimlen( ncid, id01, nyp )
      status = nf_inq_dimid ( ncid, 'CY', id01 )
      status = nf_inq_dimlen( ncid, id01, nyhalo )
      nyhalo = ( nyhalo - nyp )/2

      status = nf_inq_dimid ( ncid,  'z', id01 )
      status = nf_inq_dimlen( ncid, id01, nz )
      status = nf_inq_dimid ( ncid, 'CZ', id01 )
      status = nf_inq_dimlen( ncid, id01, nzhalo )
      nzhalo = ( nzhalo - nz )/2

      !urban
      status = nf_inq_dimid ( ncid, 'uz', id01 )
      status = nf_inq_dimlen( ncid, id01, uz )
      !land
      status = nf_inq_dimid ( ncid, 'lz', id01 )
      status = nf_inq_dimlen( ncid, id01, lz )


      count(1:4)        = (/nxp,nyp,nz,1/)
      count2(1:3)       = (/nxp,nyp,1/)
      count_urban(1:4)  = (/nxp,nyp,uz,1/)
      count_land(1:4)   = (/nxp,nyp,lz,1/)
      count_height(1:3) = (/nxp,nyp,nz/)
      allocate( p_3d(nxp,nyp,nz,1) )
      allocate( p_2d(nxp,nyp,1) )
      allocate( p_urban_3d(nxp,nyp,uz,1) )
      allocate( p_land_3d(nxp,nyp,lz,1) )
      allocate( p_height(nxp,nyp,nz) )
      !
      allocate( p_cdx(nxp+2*nxhalo) )
      allocate( p_cdy(nyp+2*nyhalo) )
      allocate( p_cx(nxp+2*nxhalo) )
      allocate( p_cy(nyp+2*nyhalo) )

    if( ofirst ) then
      ofirst = .false.
      if ( HIST_BND ) then  ! .true.
        nx = (nxp-2) * xproc + 4
        ny = (nyp-2) * yproc + 4
      else
        nx = nxp * xproc
        ny = nyp * yproc
      endif
      allocate( var3d(nx,ny,nz) )
      allocate( var2d(nx,ny) )
      allocate( var_urban_3d(nx,ny,uz) )
      allocate( var_land_3d(nx,ny,lz) )
      allocate( cdx(nx) )
      allocate( cdy(ny) )
      allocate( cx(nx) )
      allocate( cy(ny) )
      allocate( cz(nz+2*nzhalo) )
      allocate( cdz(nz+2*nzhalo) )
   endif

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
   !endif

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

    if( (trim(vname(n)) == "TRL_URB").or.  &
        (trim(vname(n)) == "TBL_URB").or.  &
        (trim(vname(n)) == "TGL_URB") ) then 
     status = nf_get_vara_double( ncid,id01,start,count_urban,p_urban_3d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double land ", trim(vname(n))
      stop
     end if
    elseif( (trim(vname(n)) == "LAND_TEMP").or.  &
            (trim(vname(n)) == "LAND_WATER") ) then 
     status = nf_get_vara_double( ncid,id01,start,count_land,p_land_3d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double urban ", trim(vname(n))
      stop
     end if
    elseif( (trim(vname(n)) == "height") )then
     status = nf_get_vara_double( ncid,id01,start,count_height,p_height )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double height ", trim(vname(n))
      stop
     end if
    elseif( ndim == 4 ) then 
     status = nf_get_vara_double( ncid,id01,start,count,p_3d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double 3 ", trim(vname(n))
      stop
     end if
    elseif( ndim == 3 ) then
     status = nf_get_vara_double( ncid,id01,start2,count2,p_2d )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double 2 ", trim(vname(n))
      stop
     end if
    end if
 
    status = nf_close(ncid)

    !--- conbine variables in each processor to single  

    if ( HIST_BND ) then  ! .true.
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
    else
       is=(ix-1)*nxp+1 ; ie=(ix-1)*nxp+nxp
       js=(jy-1)*nyp+1 ; je=(jy-1)*nyp+nyp
    endif


    if ( HIST_BND ) then  ! .true.
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
        cy(jjy)  = p_cy( jjy-((jy-1)*nyp+2)+nyhalo)
      enddo
     endif
    else
      do iix = is, ie
        cdx(iix) = p_cdx(iix-is+1+nxhalo)
        cx(iix)  = p_cx( iix-is+1+nxhalo)
      enddo
      do jjy = js, je
        cdy(jjy) = p_cdy(jjy-js+1+nyhalo)
        cy(jjy)  = p_cy( jjy-js+1+nyhalo)
      enddo
    endif

    if( (trim(vname(n)) == "TRL_URB").or.  &
        (trim(vname(n)) == "TBL_URB").or.  &
        (trim(vname(n)) == "TGL_URB") ) then
     do iix = is,ie
     do jjy = js,je
      var_urban_3d(iix,jjy,1:uz) = real(p_urban_3d(iix-is+1,jjy-js+1,1:uz,1))
     enddo
     enddo
    elseif( (trim(vname(n)) == "LAND_TEMP").or.  &
            (trim(vname(n)) == "LAND_WATER") ) then
     do iix = is,ie
     do jjy = js,je
      var_land_3d(iix,jjy,1:uz) = real(p_land_3d(iix-is+1,jjy-js+1,1:uz,1))
     enddo
     enddo
    elseif( (trim(vname(n)) == "height") )then
     do iix = is,ie
     do jjy = js,je
      var3d(iix,jjy,1:nz) = real(p_height(iix-is+1,jjy-js+1,1:nz))
     enddo
     enddo
    elseif( ndim == 4 ) then 
     do iix = is,ie
     do jjy = js,je
      var3d(iix,jjy,1:nz) = real(p_3d(iix-is+1,jjy-js+1,1:nz,1))
     enddo
     enddo
    elseif( ndim == 3 ) then
     do iix = is,ie
     do jjy = js,je
      var2d(iix,jjy) = real(p_2d(iix-is+1,jjy-js+1,1)) 
     enddo
     enddo
    endif


      deallocate( p_3d )
      deallocate( p_2d )
      deallocate( p_urban_3d )
      deallocate( p_land_3d )
      deallocate( p_height )
      deallocate( p_cdx )
      deallocate( p_cdy )
      deallocate( p_cx )
      deallocate( p_cy )

   enddo
   enddo

   if( it == nst ) then
     write(cnest,'(i2.2)') inest
     write(*,*) "create ctl file of ", trim(vname(n))

     open(10,file=trim(odir)//'/'//trim(vname(n))//'_d'//trim(cnest)//".ctl", form="formatted", access="sequential" )
     write(10,'(a,1x,a)') "DSET", "^"//trim(vname(n))//'_d'//trim(cnest)//".grd"
     write(10,'(a)') "TITLE SCALE3 data output"
     write(10,'(a)') "OPTIONS BIG_ENDIAN"
     write(10,'(a,1x,e15.7)') "UNDEF", -9.9999001E+30
     write(10,'(a,3x,i7,1x,a)') "XDEF", nx, "LEVELS"
     write(10,'(5(1x,e15.7))') cx(1:nx)*1.d-3
     write(10,'(a,3x,i7,1x,a)') "YDEF", ny, "LEVELS"
     write(10,'(5(1x,e15.7))') cy(1:ny)*1.d-3

     if( (trim(vname(n)) == "TRL_URB").or.  &
         (trim(vname(n)) == "TBL_URB").or.  &
         (trim(vname(n)) == "TGL_URB") ) then
        write(10,'(a,3x,i7,1x,a,1x,a)') "ZDEF", uz, "linear", "1 1"
     elseif( (trim(vname(n)) == "LAND_TEMP").or.  &
             (trim(vname(n)) == "LAND_WATER") ) then
        write(10,'(a,3x,i7,1x,a,1x,a)') "ZDEF", lz, "linear", "1 1"
     elseif( (trim(vname(n)) == "height") )then
        write(10,'(a,3x,i7,1x,a)') "ZDEF", nz, "LEVELS"
        write(10,'(5(1x,e15.7))') cz(nzhalo+1:nz+nzhalo)*1.d-3
     elseif( ndim == 4 ) then
        write(10,'(a,3x,i7,1x,a)') "ZDEF", nz, "LEVELS"
        write(10,'(5(1x,e15.7))') cz(nzhalo+1:nz+nzhalo)*1.d-3
     elseif( ndim == 3 ) then
        write(10,'(a,3x,i7,1x,a,1x,e15.7)') "ZDEF", 1, "LEVELS", cz(nzhalo+1)*1.d-3
     endif 

     write(10,'(a,3x,i5,1x,a,1x,a,3x,a)') "TDEF", nen-nst+1, "LINEAR", trim(stime), trim(delt)

     write(10,'(a,3x,i2)') "VARS", 1
     if( (trim(vname(n)) == "TRL_URB").or.  &
         (trim(vname(n)) == "TBL_URB").or.  &
         (trim(vname(n)) == "TGL_URB") ) then
        write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), uz, 99, "NONE"
     elseif( (trim(vname(n)) == "LAND_TEMP").or.  &
             (trim(vname(n)) == "LAND_WATER") ) then
        write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), lz, 99, "NONE"
     elseif( (trim(vname(n)) == "height") )then
        write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), nz, 99, "NONE"
     elseif( ndim == 4 ) then
        write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), nz, 99, "NONE"
     elseif( ndim == 3 ) then
        write(10,'(a,1x,i7,1x,i2,1x,a)') trim(vname(n)), 0, 99, "NONE"
     endif 
     write(10,'(a)') "ENDVARS"

     close(10)
   endif


   if( (trim(vname(n)) == "TRL_URB").or.  &
       (trim(vname(n)) == "TBL_URB").or.  &
       (trim(vname(n)) == "TGL_URB") ) then
      print *, maxval(var_urban_3d), minval(var_urban_3d)

      open(10,file=trim(odir)//'/'//trim(vname(n))//'_d'//trim(cnest)//".grd", &
           form="unformatted", access="direct", recl=4*nx*ny*uz)
      write(10,rec=it-nst+1) (((var_urban_3d(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,uz)
      close(10)
   elseif( (trim(vname(n)) == "LAND_TEMP").or.  &
           (trim(vname(n)) == "LAND_WATER") ) then
      print *, maxval(var_land_3d), minval(var_land_3d)

      open(10,file=trim(odir)//'/'//trim(vname(n))//'_d'//trim(cnest)//".grd", &
           form="unformatted", access="direct", recl=4*nx*ny*lz)
      write(10,rec=it-nst+1) (((var_land_3d(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,lz)
      close(10)
   elseif( (trim(vname(n)) == "height") ) then
      print *, maxval(var3d), minval(var3d)
      open(10,file=trim(odir)//'/'//trim(vname(n))//'_d'//trim(cnest)//".grd", &
           form="unformatted", access="direct", recl=4*nx*ny*nz)
      write(10,rec=it-nst+1) (((var3d(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
      close(10)   
   elseif( ndim == 4 ) then
      print *, maxval(var3d), minval(var3d)
      open(10,file=trim(odir)//'/'//trim(vname(n))//'_d'//trim(cnest)//".grd", &
           form="unformatted", access="direct", recl=4*nx*ny*nz)
      write(10,rec=it-nst+1) (((var3d(ix,jy,kz),ix=1,nx),jy=1,ny),kz=1,nz)
      close(10)   
   elseif( ndim == 3 ) then
      print *, maxval(var2d), minval(var2d)
      open(10,file=trim(odir)//'/'//trim(vname(n))//'_d'//trim(cnest)//".grd", &
           form="unformatted", access="direct", recl=4*nx*ny)
      write(10,rec=it-nst+1) ((var2d(ix,jy),ix=1,nx),jy=1,ny)
      close(10)   
   endif

  enddo !--- for variable (n)
  enddo !--- for time (it)

  write(*,*) "CALCULATION FINISHED"

end program

