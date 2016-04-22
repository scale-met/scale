program convine
 implicit none

  include 'netcdf.inc'

  integer(4), parameter :: xproc=2, yproc=3
  integer(4), parameter :: nst=30, nen=60
  integer(4), parameter :: nt=nen-nst+1
  integer(4), parameter :: nxp=8, nyp=8
  integer(4), parameter :: nzhalo=2
  integer(4), parameter :: nx=nxp*xproc, ny=nyp*yproc, nz=200
  real(8), parameter :: cp=1015.D0, Lp=2.47d+6, dt=60.d0 !sec
  real(8), parameter :: rdry=287.04d0 , p00=101300.d0 , cpdry=1003.5d0 
  real(8), parameter :: lh0=2.5008D+6, rovcp=rdry/cpdry, tint=1.d0*dt*0.d0
  integer(4), parameter :: nzs=200
  real(8), parameter :: s_zmax=1.2d0
  real(8) :: s_hgt(nzs), hgt(nz), tmp2(nzs), rcount(nzs)
  real(8) :: s_qt(nzs), s_ql(nzs), s_tke(nzs), s_pl(nzs)
  real(8) :: s_sgstke(nzs), s_ww(nzs), s_www(nzs)
  real(8) :: s_uu(nzs), s_lwpt(nzs), s_vv(nzs), s_wqt(nzs), s_wpt(nzs)
  real(8) :: s_totwqt(nzs), s_totwpt(nzs)
  integer(4) :: ntm1(1)

  integer(4) :: start(4) 
  integer(4) :: count(4) 
  integer(4) :: count2(3) 
  data count /nxp,nyp,nz,1/
  data count2 /nxp,nyp,1/
  integer(4) :: ix, jy, kz, kz1, nrec, n
  integer(4) :: ncid, id01, status, iix, jjy
  real(8) :: count1
  character(len=H_MID) :: cfile
  real(8),allocatable :: dens(:,:,:), qt(:,:,:), t(:,:,:), pres(:,:,:)
  real(8),allocatable :: pt(:,:,:), w(:,:,:), u(:,:,:), v(:,:,:)
  real(8),allocatable :: qc(:,:,:), qr(:,:,:), ql(:,:,:), qv(:,:,:)
  real(8),allocatable :: lwptp(:,:,:), lwpt(:,:,:), lwp(:,:)
  real(8),allocatable :: tke(:,:,:), tke2d(:,:), sgstke2d(:,:)
  real(8),allocatable :: sgsw(:,:,:), sgsu(:,:,:), sgsv(:,:,:)
  real(8),allocatable :: sgstke(:,:,:)
  real(8),allocatable :: lhflx(:,:), shflx(:,:)
  real(8),allocatable :: ww(:), www(:), uu(:), vv(:)
  real(8),allocatable :: wp(:,:,:), up(:,:,:), vp(:,:,:), qtp(:,:,:)
  real(8),allocatable :: ptp(:,:,:), wpt(:), wqt(:), wlwpt(:)
  real(8),allocatable :: sgswqt(:,:,:), sgswpt(:,:,:), sgsww(:,:,:)
  real(8),allocatable :: sgswqc(:,:,:), sgswqr(:,:,:), sgswqv(:,:,:)
  real(8),allocatable :: zi(:,:), zb(:,:), dens_ap(:,:,:)

  real(8) :: aveqt(nz), avept(nz), avew(nz), aveql(nz), avelwpt(nz)
  real(8) :: avetke(nz), aveww(nz), avewww(nz), avewpt(nz), avewqt(nz), avewlwpt(nz)
  real(8) :: avelwp(nen), avetke2d(nen), avezi(nen), avezb(nen), aveu(nz), avev(nz)
  real(8) :: profileqt(nz), profileql(nz), outavetke(2,nz)
  real(8) :: ccover(nen), avesgstke(nz), avesgstke2d(nen), aveden_a(nz)
  real(8) :: avesgswpt(nz), avesgswqt(nz), totwpt(nz), totwqt(nz), profilept(nz)
  real(8) :: aveuu(nz), avevv(nz), avesgsww(nz)
  real(8) :: cz(nz+2*nzhalo), cdz(nz+2*nzhalo)

  real(8),allocatable :: p_qt(:,:,:,:), p_t(:,:,:,:), p_qv(:,:,:,:)
  real(8),allocatable :: p_pt(:,:,:,:), p_u(:,:,:,:), p_v(:,:,:,:)
  real(8),allocatable :: p_lwpt(:,:,:,:), p_w(:,:,:,:), p_pres(:,:,:,:)
  real(8),allocatable :: p_sgstke(:,:,:,:), p_dens(:,:,:,:)
  real(8),allocatable :: p_qc(:,:,:,:), p_qr(:,:,:,:)
  real(8),allocatable :: p_sgs_wptflx(:,:,:,:), p_sgs_wqtflx(:,:,:,:)
  real(8),allocatable :: p_sgs_wwflx(:,:,:,:), p_sgs_vvflx(:,:,:,:), p_sgs_uuflx(:,:,:,:)
  real(8),allocatable :: p_lhflx(:,:,:), p_shflx(:,:,:)
  real(8),allocatable :: p_sgs_wqcflx(:,:,:,:), p_sgs_wqrflx(:,:,:,:)
  !--------------  
  allocate( lwp(nx,ny) )
  allocate( zb(nx,ny) )
  allocate( zi(nx,ny) )
  allocate( lhflx(nx,ny) )
  allocate( shflx(nx,ny) )
  allocate( qt(nx,ny,nz) )
  allocate( qv(nx,ny,nz) )
  allocate( pt(nx,ny,nz) )
  allocate( pres(nx,ny,nz) )
  allocate( lwpt(nx,ny,nz) )
  allocate( w(nx,ny,nz) )
  allocate( t(nx,ny,nz) )
  allocate( v(nx,ny,nz) )
  allocate( u(nx,ny,nz) )
  allocate( wp(nx,ny,nz) )
  allocate( vp(nx,ny,nz) )
  allocate( up(nx,ny,nz) )
  allocate( ptp(nx,ny,nz) )
  allocate( lwptp(nx,ny,nz) )
  allocate( dens(nx,ny,nz) )
  allocate( dens_ap(nx,ny,nz) )
  allocate( qc(nx,ny,nz) )
  allocate( qr(nx,ny,nz) )
  allocate( ql(nx,ny,nz) )
  allocate( qtp(nx,ny,nz) )
  allocate( tke(nx,ny,nz) )
  allocate( sgstke(nx,ny,nz) )
  allocate( ww(nz) )
  allocate( vv(nz) )
  allocate( uu(nz) )
  allocate( www(nz) )
  allocate( wpt(nz) )
  allocate( wqt(nz) )
  allocate( wlwpt(nz) )
  allocate( tke2d(nx,ny) )
  allocate( sgstke2d(nx,ny) )
  allocate( sgsw(nx,ny,nz) )
  allocate( sgsv(nx,ny,nz) )
  allocate( sgsu(nx,ny,nz) )
  allocate( sgswpt(nx,ny,nz) )
  allocate( sgswqt(nx,ny,nz) )
  allocate( sgsww(nx,ny,nz) )
  allocate( sgswqc(nx,ny,nz) )
  allocate( sgswqv(nx,ny,nz) )
  allocate( sgswqr(nx,ny,nz) )
  !--- for netCDF data
  allocate( p_qt(nxp,nyp,nz,1) )
  allocate( p_qv(nxp,nyp,nz,1) )
  allocate( p_qr(nxp,nyp,nz,1) )
  allocate( p_qc(nxp,nyp,nz,1) )
  allocate( p_pt(nxp,nyp,nz,1) )
  allocate( p_pres(nxp,nyp,nz,1) )
  allocate( p_t(nxp,nyp,nz,1) )
  allocate( p_u(nxp,nyp,nz,1) )
  allocate( p_v(nxp,nyp,nz,1) )
  allocate( p_w(nxp,nyp,nz,1) )
  allocate( p_lwpt(nxp,nyp,nz,1) )
  allocate( p_dens(nxp,nyp,nz,1) )
  allocate( p_sgstke(nxp,nyp,nz,1) )
  allocate( p_sgs_wptflx(nxp,nyp,nz,1) )
  allocate( p_sgs_wqtflx(nxp,nyp,nz,1) )
  allocate( p_sgs_wqcflx(nxp,nyp,nz,1) )
  allocate( p_sgs_wqrflx(nxp,nyp,nz,1) )
  allocate( p_sgs_wwflx(nxp,nyp,nz,1) )
  allocate( p_sgs_uuflx(nxp,nyp,nz,1) )
  allocate( p_sgs_vvflx(nxp,nyp,nz,1) )
  allocate( p_lhflx(nxp,nyp,1) )
  allocate( p_shflx(nxp,nyp,1) )
  !--------------  

  open(50,file="comp_time_evolv.txt")
  write(50,'(1x,a3,6(1x,a15))') "#  ", "LWP", "TKE(GR)", "TKE(SGS)", "Zi", "Zb", "FC"

  do kz = 1, nzs
   s_hgt(kz) = (kz-1)*s_zmax/real(nzs)
  enddo

  !--- 3D value
  profileqt(:) = 0.d0
  profileql(:) = 0.d0
  profilept(:) = 0.d0
  outavetke(:,:)=0.d0
!  avetke(:)=0.d0
!  avesgstke(:)=0.d0
  aveuu(:) = 0.d0
  avevv(:) = 0.d0
  aveww(:) = 0.d0
  avewww(:) = 0.d0
  avewpt(:) = 0.d0
  avewqt(:) = 0.d0
  avewlwpt(:) = 0.d0
  totwpt(:) = 0.d0
  totwqt(:) = 0.d0
  !--- 2D value
  avelwp(:) = 0.d0
  avezi(:) = 0.d0
  avezb(:) = 0.d0
  avetke2d(:) = 0.d0
  avesgstke2d(:) = 0.d0
  do n = 1, nen ! time loop
   !--- open NetCDF file and read from NetCDF file
    start(1:4) = (/1,1,1,n/)
    nrec = -1
    write(*,*) "TIME= ", real(n)*dt, "[sec]"
    do jy = 1, yproc !--- yproc
    do ix = 1, xproc !--- xproc
     nrec = nrec + 1
     if( nrec < 10 ) then
      write(cfile,'(a,i1,a)') "./history.pe00000",nrec,".nc"
     elseif( nrec < 100 ) then
      write(cfile,'(a,i2,a)') "./history.pe0000",nrec,".nc"
     elseif( nrec < 1000 ) then
      write(cfile,'(a,i3,a)') "./history.pe000",nrec,".nc"
     elseif( nrec < 10000 ) then
      write(cfile,'(a,i4,a)') "./history.pe00",nrec,".nc"
     elseif( nrec < 100000 ) then
      write(cfile,'(a,i5,a)') "./history.pe0",nrec,".nc"
      endif
  
     status = nf_open(cfile,0,ncid)
     if( status /= nf_noerr ) then
      write(*,*) "Stop at nf open"
      stop
     endif

     if( n == 1 ) then 
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
 
     status = nf_inq_varid( ncid,'W',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid w"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_w )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double w"
      stop
     end if
 
     status = nf_inq_varid( ncid,'V',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid v"
     stop
     end if
 
     status = nf_get_vara_double( ncid,id01,start,count,p_v )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double v"
      stop
     end if
 
     status = nf_inq_varid( ncid,'U',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid u"
     stop
     end if
 
     status = nf_get_vara_double( ncid,id01,start,count,p_u )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double u"
      stop
     end if
 
     status = nf_inq_varid( ncid,'QR',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid qr"
     stop
     end if
 
     status = nf_get_vara_double( ncid,id01,start,count,p_qr )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double qr"
      stop
     end if
 
     status = nf_inq_varid( ncid,'QC',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid qc"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_qc )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double qc"
      stop
     end if
 
     status = nf_inq_varid( ncid,'QV',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid qv"
     stop
     end if
 
     status = nf_get_vara_double( ncid,id01,start,count,p_qv )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double qv"
      stop
     end if
 
     status = nf_inq_varid( ncid,'DENS',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid dens"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_dens )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double dens"
      stop
     end if
  
     status = nf_inq_varid( ncid,'T',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid t"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_t )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double dens"
      stop
     end if
  
     status = nf_inq_varid( ncid,'TKE',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid sgstke"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_sgstke )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double sgstke"
      stop
     end if
  
     status = nf_inq_varid( ncid,'PT',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid pt"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_pt )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double pt"
      stop
     end if
  
     status = nf_inq_varid( ncid,'PRES',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid pres"
      stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_pres )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double pres"
      stop
     end if
  
     status = nf_inq_varid( ncid,'SGS_ZFLX_RHOT',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid sgs_wptflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_sgs_wptflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_sgs_wptflx"
      stop
     end if
  
     status = nf_inq_varid( ncid,'SGS_ZFLX_QV',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid sgs_wqtflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_sgs_wqtflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_sgs_wqtflx"
      stop
     end if
  
     status = nf_inq_varid( ncid,'SGS_ZFLX_QC',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid sgs_wqcflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_sgs_wqrflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_sgs_wqcflx"
      stop
     end if
  
     status = nf_inq_varid( ncid,'SGS_ZFLX_QR',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid sgs_wqrflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_sgs_wqrflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_sgs_wqrflx"
      stop
     end if
  
     status = nf_inq_varid( ncid,'SGS_ZFLX_MOMZ',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid sgs_wwflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count,p_sgs_wwflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_sgs_wwflx"
      stop
     end if
  
     status = nf_inq_varid( ncid,'LHFLX',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid lhflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count2,p_lhflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_lhflx"
      stop
     end if
  
     status = nf_inq_varid( ncid,'SHFLX',id01 )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf inq_varid shflx"
     stop
     end if
  
     status = nf_get_vara_double( ncid,id01,start,count2,p_shflx )
     if( status /= nf_noerr) then
      write(*,*) "stop at nf get_var_double p_shflx"
      stop
     end if
  
     status = nf_close(ncid)
 
     do iix = (ix-1)*nxp+1, (ix-1)*nxp+nxp
     do jjy = (jy-1)*nyp+1, (jy-1)*nyp+nyp
      w(iix,jjy,1:nz) = real(p_w(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      u(iix,jjy,1:nz) = real(p_u(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      v(iix,jjy,1:nz) = real(p_v(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      t(iix,jjy,1:nz) = real(p_t(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      pt(iix,jjy,1:nz) = real(p_pt(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      pres(iix,jjy,1:nz) = real(p_pres(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      dens(iix,jjy,1:nz) = real(p_dens(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      qc(iix,jjy,1:nz) = real(p_qc(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      qv(iix,jjy,1:nz) = real(p_qv(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      qr(iix,jjy,1:nz) = real(p_qr(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      sgstke(iix,jjy,1:nz) = real(p_sgstke(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      sgsww(iix,jjy,1:nz) = real(p_sgs_wwflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      sgswpt(iix,jjy,1:nz) = real(p_sgs_wptflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      sgswqv(iix,jjy,1:nz) = real(p_sgs_wqtflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      sgswqc(iix,jjy,1:nz) = real(p_sgs_wqcflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      sgswqr(iix,jjy,1:nz) = real(p_sgs_wqrflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1:nz,1))
      lhflx(iix,jjy) = real(p_lhflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1))
      shflx(iix,jjy) = real(p_shflx(iix-(ix-1)*nxp,jjy-(jy-1)*nyp,1))
     enddo
     enddo
    enddo !--- xproc
    enddo !--- yproc

    if( n == 1 ) then
     do kz = 1, nz
      hgt(kz) = cz(kz+nzhalo)
     enddo
    endif
 
    aveden_a(:) = 0.d0
    do kz = 1, nz
    do jy = 1, ny
    do ix = 1, nx 
     aveden_a(kz) = aveden_a(kz)+dens(ix,jy,kz)
    end do
    end do
    end do
    aveden_a(1:nz) = aveden_a(1:nz)/real(nx*ny)
    do kz = 1, nz
    do jy = 1, ny
    do ix = 1, nx
     dens_ap(ix,jy,kz) = dens(ix,jy,kz)-aveden_a(kz)
    end do
    end do
    end do
    do jy = 1, ny
    do ix = 1, nx
    do kz = 2, nz
     sgsww (ix,jy,kz) = ( sgsww (ix,jy,kz-1)+sgsww (ix,jy,kz  ) )*0.5d0
     sgswpt(ix,jy,kz) = ( sgswpt(ix,jy,kz-1)+sgswpt(ix,jy,kz  ) )*0.5d0
     sgswqt(ix,jy,kz) = ( sgswqc(ix,jy,kz-1)+sgswqr(ix,jy,kz-1)+sgswqv(ix,jy,kz-1) &
                        + sgswqc(ix,jy,kz  )+sgswqr(ix,jy,kz  )+sgswqv(ix,jy,kz  ) )*0.5d0
    end do
     sgswpt(ix,jy,1) = 0.d0
     sgswqt(ix,jy,1) = 0.d0
     sgsww (ix,jy,1) = 0.d0
    end do
    end do
 
    zi(:,:) = 0.d0
    do ix = 1, nx
    do jy = 1, ny
     loopzi : do kz = nz, 1, -1
      if( qv(ix,jy,kz)+qc(ix,jy,kz)+qr(ix,jy,kz) > 8.d-3 ) then
       zi(ix,jy) = cz(kz+nzhalo)
       exit loopzi
      end if
     end do loopzi
    end do 
    end do

    zb(:,:) = 0.d0
    do ix = 1, nx
    do jy = 1, ny
     loopzb : do kz = 1, nz
      if( qc(ix,jy,kz)+qr(ix,jy,kz) > 1.d-5 ) then
       zb(ix,jy) = cz(kz+nzhalo)
       exit loopzb
      end if
     end do loopzb
    end do 
    end do
  

   !  average 
   aveql(:)=0.d0
   avelwpt(:)=0.d0
   aveqt(:)=0.d0
   avept(:)=0.d0
   avew(:) = 0.d0
   avev(:) = 0.d0
   aveu(:) = 0.d0
   avesgswpt(:) = 0.d0
   avesgswqt(:) = 0.d0
   tke2d(:,:)=0.d0
   count1 = 0.d0
   avesgsww(:) = 0.d0
   do kz = 1, nz
    do ix = 1, nx
    do jy = 1, ny
     lwpt(ix,jy,kz) = pt(ix,jy,kz) &
                    - ( lh0/cpdry*( qc(ix,jy,kz)+qr(ix,jy,kz) ) ) &
                    * ( p00/pres(ix,jy,kz) )**rovcp
     aveql(kz) = aveql(kz)+(qc(ix,jy,kz)+qr(ix,jy,kz))*dens(ix,jy,kz)
     aveqt(kz) = aveqt(kz)+(qv(ix,jy,kz)+qc(ix,jy,kz)+qr(ix,jy,kz))*dens(ix,jy,kz)
     avept(kz) = avept(kz)+pt(ix,jy,kz)*dens(ix,jy,kz)
     avew(kz)  = avew(kz)+w(ix,jy,kz)*dens(ix,jy,kz)
     avev(kz)  = avev(kz)+v(ix,jy,kz)*dens(ix,jy,kz)
     aveu(kz)  = aveu(kz)+u(ix,jy,kz)*dens(ix,jy,kz)
     avelwpt(kz) = avelwpt(kz)+lwpt(ix,jy,kz)*dens(ix,jy,kz)
     if( kz /= 1 ) then
      avesgswpt(kz) = avesgswpt(kz)+sgswpt(ix,jy,kz)*cp
      avesgswqt(kz) = avesgswqt(kz)+sgswqt(ix,jy,kz)*Lp
     elseif( kz == 1 ) then
      avesgswpt(kz) = avesgswpt(kz)+sgswpt(ix,jy,kz)*cp+shflx(ix,jy)
      avesgswqt(kz) = avesgswqt(kz)+sgswqt(ix,jy,kz)*Lp+lhflx(ix,jy)
     endif
     avesgsww(kz) = avesgsww(kz)+sgsww(ix,jy,kz)
    end do
    end do
    avew(kz)  = avew(kz)/real(nx*ny)/aveden_a(kz)
    avev(kz)  = avev(kz)/real(nx*ny)/aveden_a(kz)
    aveu(kz)  = aveu(kz)/real(nx*ny)/aveden_a(kz)
    avept(kz)  = avept(kz)/real(nx*ny)/aveden_a(kz)
    aveql(kz) = aveql(kz)/real(nx*ny)/aveden_a(kz)
    aveqt(kz) = aveqt(kz)/real(nx*ny)/aveden_a(kz)
    avelwpt(kz) = avelwpt(kz)/real(nx*ny)/aveden_a(kz)
    avesgswpt(kz) = avesgswpt(kz)/real(nx*ny)/aveden_a(kz)
    avesgswqt(kz) = avesgswqt(kz)/real(nx*ny)/aveden_a(kz)
    avesgsww(kz) = avesgsww(kz)/real(nx*ny)/aveden_a(kz)

    !--- calculate prime elements ( written as *p(nx,ny,nz) )
    do ix = 1, nx
    do jy = 1, ny
     wp(ix,jy,kz) = w(ix,jy,kz)-avew(kz)
     vp(ix,jy,kz) = v(ix,jy,kz)-avev(kz)
     up(ix,jy,kz) = u(ix,jy,kz)-aveu(kz)
     ptp(ix,jy,kz) = pt(ix,jy,kz)-avept(kz)
     qtp(ix,jy,kz) = (qv(ix,jy,kz)+qc(ix,jy,kz)+qr(ix,jy,kz))-aveqt(kz)
     lwptp(ix,jy,kz) = lwpt(ix,jy,kz)-avelwpt(kz)
    end do
    end do
   
    do ix = 1, nx
    do jy = 1, ny
     tke(ix,jy,kz) = 0.5d0 * dens(ix,jy,kz) / aveden_a(kz) &
                   * ( up(ix,jy,kz)**2+vp(ix,jy,kz)**2+wp(ix,jy,kz)**2 )
    end do
    end do
   enddo

   !--- variance and skewness 
   vv(:) = 0.d0
   uu(:) = 0.d0
   ww(:) = 0.d0
   www(:) = 0.d0
   wpt(:) = 0.d0
   wqt(:) = 0.d0
   wlwpt(:) = 0.d0
   avetke(:) = 0.d0 
   avesgstke(:) = 0.d0 
   if( n >= nst .and.  n <= nen ) then
    do kz = 1, nz
     do ix = 1, nx
     do jy = 1, ny
      uu(kz) = uu(kz)+up(ix,jy,kz)*up(ix,jy,kz)*dens(ix,jy,kz)
      vv(kz) = vv(kz)+vp(ix,jy,kz)*vp(ix,jy,kz)*dens(ix,jy,kz)
      ww(kz) = ww(kz)+wp(ix,jy,kz)*wp(ix,jy,kz)*dens(ix,jy,kz)
      www(kz) = www(kz)+wp(ix,jy,kz)*wp(ix,jy,kz)*wp(ix,jy,kz)*dens(ix,jy,kz)
      wpt(kz) = wpt(kz)+wp(ix,jy,kz)*ptp(ix,jy,kz)*dens(ix,jy,kz)*cp*dens(ix,jy,kz)
      wqt(kz) = wqt(kz)+wp(ix,jy,kz)*qtp(ix,jy,kz)*dens(ix,jy,kz)*Lp*dens(ix,jy,kz)
      wlwpt(kz) = wlwpt(kz)+wp(ix,jy,kz)*lwptp(ix,jy,kz)*dens(ix,jy,kz)*cp*dens(ix,jy,kz)
      avesgstke(kz) = avesgstke(kz)+sgstke(ix,jy,kz)*dens(ix,jy,kz)
      avetke(kz) = avetke(kz)+tke(ix,jy,kz)*dens(ix,jy,kz)
     end do
     end do
     uu(kz) = uu(kz)/real(nx*ny)/aveden_a(kz)
     vv(kz) = vv(kz)/real(nx*ny)/aveden_a(kz)
     ww(kz) = ww(kz)/real(nx*ny)/aveden_a(kz)
     www(kz) = www(kz)/real(nx*ny)/aveden_a(kz)
     wpt(kz) = wpt(kz)/real(nx*ny)/aveden_a(kz)
     wqt(kz) = wqt(kz)/real(nx*ny)/aveden_a(kz)
     wlwpt(kz) = wlwpt(kz)/real(nx*ny)/aveden_a(kz)
     avesgstke(kz) = avesgstke(kz)/real(nx*ny)/aveden_a(kz)
     avetke(kz) = avetke(kz)/real(nx*ny)/aveden_a(kz)
 
     aveuu(kz) = aveuu(kz)+uu(kz)
     avevv(kz) = avevv(kz)+vv(kz)
     aveww(kz) = aveww(kz)+ww(kz)
     avewww(kz) = avewww(kz)+www(kz)
     avewpt(kz) = avewpt(kz)+wpt(kz)
     avewqt(kz) = avewqt(kz)+wqt(kz)
     avewlwpt(kz) = avewlwpt(kz)+wlwpt(kz)
     totwpt(kz) = totwpt(kz)+avesgswpt(kz)
     totwqt(kz) = totwqt(kz)+avesgswqt(kz)
     profileql(kz) = profileql(kz)+aveql(kz)*1.d+3  ! [g/kg]
     profileqt(kz) = profileqt(kz)+aveqt(kz)*1.d+3  ! [g/kg]
     profilept(kz) = profilept(kz)+avelwpt(kz)
     outavetke(1,kz) = outavetke(1,kz)+avetke(kz)
     outavetke(2,kz) = outavetke(2,kz)+avesgstke(kz)
    end do
   endif

   !--- 2D value
   lwp(:,:) = 0.d0
   tke2d(:,:) = 0.d0
   sgstke2d(:,:) = 0.d0
   count1 = 0.d0
   do jy = 1, ny
   do ix = 1, nx
   do kz = 1, nz 
    lwp(ix,jy) =lwp(ix,jy)+( qc(ix,jy,kz)+qr(ix,jy,kz) )*dens(ix,jy,kz)*1000.d0*cdz(kz+nzhalo) ! [g/m^2]
    tke2d(ix,jy) = tke2d(ix,jy)+tke(ix,jy,kz)*cdz(kz+nzhalo)
    sgstke2d(ix,jy) = sgstke2d(ix,jy)+ (sgstke(ix,jy,kz))*cdz(kz+nzhalo)
   end do

    avelwp(n) = avelwp(n)+lwp(ix,jy)
    avetke2d(n) = avetke2d(n)+tke2d(ix,jy)
    avesgstke2d(n) = avesgstke2d(n)+sgstke2d(ix,jy)
    avezi(n) = avezi(n)+zi(ix,jy)
    if( zb(ix,jy) /= 0.d0 ) then
     count1 = count1+1.d0
     avezb(n) = avezb(n)+zb(ix,jy)
    end if
   end do
   end do

   ccover(n) = count1/real(nx*ny)
   avelwp(n) = avelwp(n)/real(nx*ny)
   avetke2d(n) = avetke2d(n)/real(nx*ny)
   avesgstke2d(n) = avesgstke2d(n)/real(nx*ny)
   avezi(n) = avezi(n)/real(nx*ny)
   if( count1 /= 0.d0 ) then
    avezb(n) = avezb(n)/count1
   else
    avezb(n) = 0.d0
   endif
   write(50,'(1x,f10.3,7(1x,e15.7))') tint+int(real(n)*dt), avelwp(n), avetke2d(n), avesgstke2d(n), avezi(n), avezb(n), ccover(n)

  enddo !--- time loop
 
  aveuu(1:nz) = aveuu(1:nz)/real(nt)
  avevv(1:nz) = avevv(1:nz)/real(nt)
  aveww(1:nz) = aveww(1:nz)/real(nt)
  avewww(1:nz) = avewww(1:nz)/real(nt)
  avewpt(1:nz) = avewpt(1:nz)/real(nt)
  avewqt(1:nz) = avewqt(1:nz)/real(nt)
  avewlwpt(1:nz) = avewlwpt(1:nz)/real(nt)
  totwpt(1:nz) = totwpt(1:nz)/real(nt)
  totwqt(1:nz) = totwqt(1:nz)/real(nt)

  profileql(1:nz) = profileql(1:nz)/real(nt)
  profileqt(1:nz) = profileqt(1:nz)/real(nt)
  profilept(1:nz) = profilept(1:nz)/real(nt)
  outavetke(1,1:nz) = outavetke(1,1:nz)/real(nt) 
  outavetke(2,1:nz) = outavetke(2,1:nz)/real(nt) 
!  avesgstke(1:nz) = avesgstke(1:nz)/real(nt) 
!  avetke(1:nz) = avetke(1:nz)/real(nt)


  rcount(:) = 0.d0
  s_qt(:) = 0.d0
  s_ql(:) = 0.d0
  s_tke(:) = 0.d0
  s_sgstke(:) = 0.d0
  s_ww(:) = 0.d0
  s_www(:) = 0.d0
  s_uu(:) = 0.d0
  s_lwpt(:) = 0.d0
  s_vv(:) = 0.d0
  s_wqt(:) = 0.d0
  s_wpt(:) = 0.d0
  s_totwqt(:) = 0.d0
  s_totwpt(:) = 0.d0
  rcount(:) = 0.d0
  do n = nst, nen
  do kz = 1,  nz
   do kz1 = 1, nzs
    tmp2(kz1) = abs( hgt(kz)/avezi(n)-s_hgt(kz1) )
   enddo
   ntm1(1) = minloc(tmp2,1)
   rcount(ntm1(1)) = rcount(ntm1(1))+1.d0
   s_qt(ntm1(1)) = s_qt(ntm1(1))+profileqt(kz)
   s_ql(ntm1(1)) = s_ql(ntm1(1))+profileql(kz)
   s_pl(ntm1(1)) = s_pl(ntm1(1))+profilept(kz)
   s_tke(ntm1(1)) = s_tke(ntm1(1))+outavetke(1,kz)
   s_sgstke(ntm1(1)) = s_sgstke(ntm1(1))+outavetke(2,kz)
!   s_tke(ntm1(1)) = s_tke(ntm1(1))+avetke(kz)
!   s_sgstke(ntm1(1)) = s_sgstke(ntm1(1))+avesgstke(kz)
   s_uu(ntm1(1)) = s_uu(ntm1(1))+aveuu(kz)
   s_vv(ntm1(1)) = s_vv(ntm1(1))+avevv(kz)
   s_ww(ntm1(1)) = s_ww(ntm1(1))+aveww(kz)
   s_www(ntm1(1)) = s_www(ntm1(1))+avewww(kz)
   s_wqt(ntm1(1)) = s_wqt(ntm1(1))+avewqt(kz)
   s_wpt(ntm1(1)) = s_wpt(ntm1(1))+avewpt(kz)
   s_totwqt(ntm1(1)) = s_totwqt(ntm1(1))+totwqt(kz)
   s_totwpt(ntm1(1)) = s_totwpt(ntm1(1))+totwpt(kz)
  enddo
  enddo

  do kz = 1, nzs
   write(*,*)  kz, rcount(kz)
  enddo

  do kz1 = 1, nzs
!   if( rcount(kz1) /= 0.d0 ) then
    s_qt(kz1) = s_qt(kz1)/rcount(kz1)
    s_ql(kz1) = s_ql(kz1)/rcount(kz1)
    s_pl(kz1) = s_pl(kz1)/rcount(kz1)
    s_tke(kz1) = s_tke(kz1)/rcount(kz1)
    s_sgstke(kz1) = s_sgstke(kz1)/rcount(kz1)
    s_uu(kz1) = s_uu(kz1)/rcount(kz1)
    s_vv(kz1) = s_vv(kz1)/rcount(kz1)
    s_ww(kz1) = s_ww(kz1)/rcount(kz1)
    s_www(kz1) = s_www(kz1)/rcount(kz1)
    s_lwpt(kz1) = s_lwpt(kz1)/rcount(kz1)
    s_wqt(kz1) = s_wqt(kz1)/rcount(kz1)
    s_wpt(kz1) = s_wpt(kz1)/rcount(kz1)
    s_totwqt(kz1) = s_totwqt(kz1)/rcount(kz1)
    s_totwpt(kz1) = s_totwpt(kz1)/rcount(kz1)
!   endif
  enddo

  !--- write
  open(51,file="r_comp_profile.txt")
  open(52,file="comp_profile.txt")

  write(51,'(1x,a10,13(1x,a15))') "# z", "QT",  "QL", "TKE-gr", "TKE-sg-diag", "WW", "WWW", &
                                  "WPT(GR)", "WQT(GR)", "WPT(SGS)", "WQT(SGS)", "PT", "UU", "VV"
  do kz = 1, nzs
   write(51,'(1x,f10.5,14(1x,e15.7))') s_hgt(kz), s_qt(kz), s_ql(kz), s_tke(kz), s_sgstke(kz), s_ww(kz), &
                                       s_www(kz), s_wpt(kz), s_wqt(kz), s_totwpt(kz), s_totwqt(kz), &
                                       s_pl(kz), s_uu(kz), s_vv(kz)
  end do
  write(52,'(1x,a10,13(1x,a15))') "# z", "QT",  "QL", "TKE-gr", "TKE-sg-diag", "WW", "WWW", &
                                  "WPT(GR)", "WQT(GR)", "WPT(SGS)", "WQT(SGS)", "PT", "UU", "VV"
  do kz = 1, nz 
   write(52,'(1x,f10.5,13(1x,e15.7))') cz(kz+nzhalo), profileqt(kz), profileql(kz), outavetke(1,kz), &
                                       outavetke(2,kz), aveww(kz), avewww(kz), avewlwpt(kz), avewqt(kz), &
                                       totwpt(kz), totwqt(kz), profilept(kz), aveuu(kz), avevv(kz)
  end do

  close(50)
  close(51)
  close(52)
  write(*,*) "CALCULATION FINISHED"

end program
