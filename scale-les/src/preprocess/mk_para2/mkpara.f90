!-----------------------------------------------------------------------------------------
!
!+  mkpara.f90 : program for making parameter file used in cloud microphysics
!
!-----------------------------------------------------------------------------------------
!
!--- History
!
!    2006/3/24  K.Suzuki  created
!
!-----------------------------------------------------------------------------------------
program mkpara


implicit none

  include 'setup.h'

!  integer, parameter :: nspc = 7
!  integer, parameter :: nbin = 30, ncld = 10, ndrz = 20 
  integer, parameter :: ndat = 33, icemax = 3, kdeg = 4, ldeg = 4

  real(8) :: xmss( nspc,ndat ), zcap( nspc,ndat ), vtrm( nspc,ndat )
  real(8) :: blkr( nspc,ndat ), blkd( nspc,ndat ), ykrn( nspc,nspc,ndat,ndat )

  real(8) :: ywll( ndat,ndat ), ywli( ndat,ndat,icemax ), ywls( ndat,ndat )
  real(8) :: ywlg( ndat,ndat ), ywlh( ndat,ndat )

  real(8) :: ywil( ndat,ndat,icemax ), ywii( ndat,ndat,icemax,icemax )
  real(8) :: ywis( ndat,ndat,icemax ), ywig( ndat,ndat,icemax )
  real(8) :: ywih( ndat,ndat,icemax )

  real(8) :: ywsl( ndat,ndat ), ywsi( ndat,ndat,icemax ), ywss( ndat,ndat )
  real(8) :: ywsg( ndat,ndat ), ywsh( ndat,ndat )

  real(8) :: ywgl( ndat,ndat ), ywgi( ndat,ndat,icemax ), ywgs( ndat,ndat )
  real(8) :: ywgg( ndat,ndat ), ywgh( ndat,ndat )

  real(8) :: ywhl( ndat,ndat ), ywhi( ndat,ndat,icemax ), ywhs( ndat,ndat )
  real(8) :: ywhg( ndat,ndat ), ywhh( ndat,ndat )

  real(8) :: radc( nbin ), xctr( nbin ), xbnd( nbin+1 ), dxmic
  real(8) :: cctr( nspc,nbin ), cbnd( nspc,nbin+1 )
  real(8) :: ck( nspc,nspc,nbin,nbin )
  real(8) :: vt( nspc,nbin )
  real(8) :: br( nspc,nbin ), bd( nspc,nbin )

  integer :: i, j
  !===== Y.SATO add 2007/10/19 =======================================
!  integer, parameter :: kphase = 2
     !	( Hydro-dynamic = 0, Long = 1, Golovin = 2 )
  !===================================================================
!-----------------------------------------------------------------------------------------

  !--- file reading
  call rdkdat

  !--- grid setting
  call sdfgrid

  !--- capacity
  call getcp

  !--- collection kernel
  call getck

  !--- terminal velocity
  call getvt

  !--- bulk radius
  call getbr

  !--- output
  call paraout

  stop
  !---------------------------------------------------------------------------------------
contains
  !---------------------------------------------------------------------------------------
  subroutine rdkdat

  integer, parameter :: il = 1, ic = 2, ip = 3, id = 4
  integer, parameter :: is = 5, ig = 6, ih = 7
  integer, parameter :: icemax = 3
  
  real(8) :: xl( ndat ), rlec( ndat ), vrl( ndat )
  real(8) :: blkradl( ndat ), blkdnsl( ndat )

  real(8) :: xi( ndat,icemax ), riec( ndat,icemax ), vri( ndat,icemax )
  real(8) :: blkradi( ndat,icemax ), blkdnsi( ndat,icemax )
  real(8) :: xc( ndat ), xp( ndat ), xd( ndat )
  real(8) :: rcec( ndat ), rpec( ndat ), rdec( ndat )
  
  real(8) :: xs( ndat ), rsec( ndat ), vrs( ndat )
  real(8) :: blkrads( ndat ), blkdnss( ndat )

  real(8) :: xg( ndat ), rgec( ndat ), vrg( ndat )
  real(8) :: blkradg( ndat ), blkdnsg( ndat )

  real(8) :: xh( ndat ), rhec( ndat ), vrh( ndat )
  real(8) :: blkradh( ndat ), blkdnsh( ndat )

  integer :: k
 
  !--- file opening
  open ( 51, file = 'masses.asc', form = 'formatted', status = 'old' )
  open ( 52, file = 'capacity.asc', form = 'formatted', status = 'old' )
  open ( 53, file = 'termvels.asc', form = 'formatted', status = 'old' )
  open ( 54, file = 'bulkradii.asc_s_0_03_0_9', form = 'formatted', status = 'old' )
  open ( 55, file = 'bulkdens.asc_s_0_03_0_9', form = 'formatted', status = 'old' )
  open ( 56, file = 'kernels.asc_s_0_03_0_9', form = 'formatted', status = 'old' )

  !--- mass
  read( 51,900 ) xl, xi, xs, xg, xh
  do k = 1, ndat
    xmss( il,k ) = log( xl( k )*1.d-03 )
    xmss( ic,k ) = log( xi( k,1 )*1.d-03 )
    xmss( ip,k ) = log( xi( k,2 )*1.d-03 )
    xmss( id,k ) = log( xi( k,3 )*1.d-03 )
    xmss( is,k ) = log( xs( k )*1.d-03 )
    xmss( ig,k ) = log( xg( k )*1.d-03 )
    xmss( ih,k ) = log( xh( k )*1.d-03 )
  end do

  !--- capacity
  read( 52,900 ) rlec, riec, rsec, rgec, rhec
  do k = 1, ndat ! cm -> m
    zcap( il,k ) = rlec( k )*1.d-02
    zcap( ic,k ) = riec( k,1 )*1.d-02
    zcap( ip,k ) = riec( k,2 )*1.d-02
    zcap( id,k ) = riec( k,3 )*1.d-02
    zcap( is,k ) = rsec( k )*1.d-02
    zcap( ig,k ) = rgec( k )*1.d-02
    zcap( ih,k ) = rhec( k )*1.d-02
  end do

  !--- terminal velocity
  read( 53,900 ) vrl, vri, vrs, vrg, vrh
  do k = 1, ndat ! cm/s -> m/s
    vtrm( il,k ) = vrl( k )*1.d-02
    vtrm( ic,k ) = vri( k,1 )*1.d-02
    vtrm( ip,k ) = vri( k,2 )*1.d-02
    vtrm( id,k ) = vri( k,3 )*1.d-02
    vtrm( is,k ) = vrs( k )*1.d-02
    vtrm( ig,k ) = vrg( k )*1.d-02
    vtrm( ih,k ) = vrh( k )*1.d-02
  end do

  !--- bulk radii
  read( 54,900 ) blkradl, blkradi, blkrads, blkradg, blkradh
  do k = 1, ndat ! cm -> mm
    blkr( il,k ) = blkradl( k )*10.D0
    blkr( ic,k ) = blkradi( k,1 )*10.D0
    blkr( ip,k ) = blkradi( k,2 )*10.D0
    blkr( id,k ) = blkradi( k,3 )*10.D0
    blkr( is,k ) = blkrads( k )*10.D0
    blkr( ig,k ) = blkradg( k )*10.D0
    blkr( ih,k ) = blkradh( k )*10.D0
  end do

  !--- bulk density
  read( 55,900 ) blkdnsl, blkdnsi, blkdnss, blkdnsg, blkdnsh
  do k = 1, ndat ! g/cm^3 -> kg/m^3
    blkd( il,k ) = blkdnsl( k )*1000.D0
    blkd( ic,k ) = blkdnsi( k,1 )*1000.D0
    blkd( ip,k ) = blkdnsi( k,2 )*1000.D0
    blkd( id,k ) = blkdnsi( k,3 )*1000.D0
    blkd( is,k ) = blkdnss( k )*1000.D0
    blkd( ig,k ) = blkdnsg( k )*1000.D0
    blkd( ih,k ) = blkdnsh( k )*1000.D0
  end do

  !--- collection kernel
  read( 56,900 ) ywll, ywli, ywls, ywlg, ywlh, &
                 ywil, ywii, ywis, ywig, ywih, &
                 ywsl, ywsi, ywss, ywsg, ywsh, &
                 ywgl, ywgi, ywgs, ywgg, ywgh, &
                 ywhl, ywhi, ywhs, ywhg, ywhh

  ! cm**3/s -> m**3/s
  do i = 1, ndat
  do j = 1, ndat
    ykrn( il,il,i,j ) = ywll( i,j )*1.D-06
    ykrn( il,ic,i,j ) = ywli( i,j,1 )*1.D-06
    ykrn( il,ip,i,j ) = ywli( i,j,2 )*1.D-06
    ykrn( il,id,i,j ) = ywli( i,j,3 )*1.D-06
    ykrn( il,is,i,j ) = ywls( i,j )*1.D-06
    ykrn( il,ig,i,j ) = ywlg( i,j )*1.D-06
    ykrn( il,ih,i,j ) = ywlh( i,j )*1.D-06
 
    ykrn( ic,il,i,j ) = ywil( i,j,1 )*1.D-06
    ykrn( ip,il,i,j ) = ywil( i,j,2 )*1.D-06
    ykrn( id,il,i,j ) = ywil( i,j,3 )*1.D-06
    ykrn( ic,ic,i,j ) = ywii( i,j,1,1 )*1.D-06
    ykrn( ic,ip,i,j ) = ywii( i,j,1,2 )*1.D-06
    ykrn( ic,id,i,j ) = ywii( i,j,1,3 )*1.D-06
    ykrn( ip,ic,i,j ) = ywii( i,j,2,1 )*1.D-06
    ykrn( ip,ip,i,j ) = ywii( i,j,2,2 )*1.D-06
    ykrn( ip,id,i,j ) = ywii( i,j,2,3 )*1.D-06
    ykrn( id,ic,i,j ) = ywii( i,j,3,1 )*1.D-06
    ykrn( id,ip,i,j ) = ywii( i,j,3,2 )*1.D-06
    ykrn( id,id,i,j ) = ywii( i,j,3,3 )*1.D-06
    ykrn( ic,is,i,j ) = ywis( i,j,1 )*1.D-06
    ykrn( ip,is,i,j ) = ywis( i,j,2 )*1.D-06
    ykrn( id,is,i,j ) = ywis( i,j,3 )*1.D-06
    ykrn( ic,ig,i,j ) = ywig( i,j,1 )*1.D-06
    ykrn( ip,ig,i,j ) = ywig( i,j,2 )*1.D-06
    ykrn( id,ig,i,j ) = ywig( i,j,3 )*1.D-06
    ykrn( ic,ih,i,j ) = ywih( i,j,1 )*1.D-06
    ykrn( ip,ih,i,j ) = ywih( i,j,2 )*1.D-06
    ykrn( id,ih,i,j ) = ywih( i,j,3 )*1.D-06 

    ykrn( is,il,i,j ) = ywsl( i,j )*1.D-06
    ykrn( is,ic,i,j ) = ywsi( i,j,1 )*1.D-06
    ykrn( is,ip,i,j ) = ywsi( i,j,2 )*1.D-06
    ykrn( is,id,i,j ) = ywsi( i,j,3 )*1.D-06
    ykrn( is,is,i,j ) = ywss( i,j )*1.D-06
    ykrn( is,ig,i,j ) = ywsg( i,j )*1.D-06
    ykrn( is,ih,i,j ) = ywsh( i,j )*1.D-06

    ykrn( ig,il,i,j ) = ywgl( i,j )*1.D-06
    ykrn( ig,ic,i,j ) = ywgi( i,j,1 )*1.D-06
    ykrn( ig,ip,i,j ) = ywgi( i,j,2 )*1.D-06
    ykrn( ig,id,i,j ) = ywgi( i,j,3 )*1.D-06
    ykrn( ig,is,i,j ) = ywgs( i,j )*1.D-06
    ykrn( ig,ig,i,j ) = ywgg( i,j )*1.D-06
    ykrn( ig,ig,i,j ) = ywgh( i,j )*1.D-06

    ykrn( ih,il,i,j ) = ywhl( i,j )*1.D-06
    ykrn( ih,ic,i,j ) = ywhi( i,j,1 )*1.D-06
    ykrn( ih,ip,i,j ) = ywhi( i,j,2 )*1.D-06
    ykrn( ih,id,i,j ) = ywhi( i,j,3 )*1.D-06
    ykrn( ih,is,i,j ) = ywhs( i,j )*1.D-06
    ykrn( ih,ig,i,j ) = ywhg( i,j )*1.D-06
    ykrn( ih,ih,i,j ) = ywhh( i,j )*1.D-06
  end do
  end do 

  close ( 51 )
  close ( 52 )
  close ( 53 )
  close ( 54 )
  close ( 55 )
  close ( 56 )

900 format ( 6e13.5 )

  end subroutine rdkdat
  !---------------------------------------------------------------------------------------
  subroutine sdfgrid

  use mod_const, only : rhow => CNST_RHOW, &
                        pi   => CNST_PI

  real(8) :: xsta, xend
  integer :: n

  xsta = log( rhow * 4.D0*pi/3.D0 * ( 3.D-06 )**3 )
  xend = log( rhow * 4.D0*pi/3.D0 * ( 3.D-03 )**3 )

  dxmic = ( xend-xsta )/nbin
  do n = 1, nbin+1
    xbnd( n ) = xsta + dxmic*( n-1 )
  end do
  do n = 1, nbin
    xctr( n ) = ( xbnd( n )+xbnd( n+1 ) )*0.5D0
    radc( n ) = ( exp( xctr( n ) )*3.D0/4.D0/pi/rhow )**( 1.D0/3.D0 )
  end do

  end subroutine sdfgrid
  !---------------------------------------------------------------------------------------
  subroutine getcp

  integer :: myu, n

  do myu = 1, nspc
    do n = 1, nbin
      cctr( myu,n ) = fcpc( myu,xctr( n ) )
    end do
    do n = 1, nbin+1
      cbnd( myu,n ) = fcpc( myu,xbnd( n ) )
    end do
  end do

  end subroutine getcp
  !---------------------------------------------------------------------------------------
  function fcpc( myu,x )

  use mod_spline, only: getknot, getmatrx, getcoef, fspline

  integer, intent(in) :: myu
  real(8), intent(in) :: x
  real(8) :: fcpc

  real(8) :: qknt( ndat+kdeg ), elm( ndat,ndat ), coef( ndat )

  call getknot                        &
         ( ndat, kdeg, xmss( myu,: ), & !--- in
           qknt                       ) !--- out

  call getmatrx                             &
         ( ndat, kdeg, qknt, xmss( myu,: ), & !--- in
           elm                              ) !--- out 

  call getcoef                             &
         ( ndat, kdeg, elm, zcap( myu,: ), & !--- in
           coef                            ) !--- out

  fcpc = fspline ( ndat, kdeg, coef, qknt, xmss( myu,: ), x )

  end function fcpc
  !---------------------------------------------------------------------------------------
  subroutine getck

  integer :: myu, nyu, i, j

  do myu = 1, nspc
  do nyu = 1, nspc
  write( *,* ) ' myu, nyu :', myu, nyu
  do i = 1, nbin
  do j = 1, nbin
    ck( myu,nyu,i,j ) = fckrn( myu,nyu,xctr( i ),xctr( j ) )
  end do
  end do
  end do
  end do 

  return

  end subroutine getck
  !---------------------------------------------------------------------------------------
  function fckrn( myu,nyu,x,y )

  use mod_spline_2d, only: getknot, getcoef2, fspline2
  use mod_const, only: rhow => CNST_RHOW, &
                       pi => CNST_PI

  integer, intent(in) :: myu, nyu
  real(8), intent(in) :: x, y
  real(8) :: fckrn

  real(8) :: qknt( ndat+kdeg ), rknt( ndat+kdeg )
  real(8) :: coef( ndat,ndat )
  
  real(8) :: xlrg, xsml, vlrg, vsml, rlrg

  if( kphase == 0 ) then
   call getknot                        &
         ( ndat, kdeg, xmss( myu,: ), & !--- in
           qknt                       ) !--- out

   rknt( : ) = qknt( : )

   call getcoef2                          &
         ( ndat, ndat, kdeg, kdeg,       & !--- in
           xmss( myu,: ), xmss( nyu,: ), & !--- in 
           qknt, rknt,                   & !--- in
           ykrn( myu,nyu,:,: ),          & !--- in
           coef                          ) !--- out

   fckrn = fspline2                          &
            ( ndat, ndat, kdeg, kdeg,       & !--- in
              coef, qknt, rknt,             & !--- in
              xmss( myu,: ), xmss( nyu,: ), & !--- in
              x, y                          ) !--- in 
  else if( kphase == 1 ) then
   xlrg = max( x, y )
   xsml = min( x, y )
   
   vlrg = (exp( xlrg ) / rhow )*1.D+06
   vsml = (exp( xsml ) / rhow )*1.D+06

   rlrg = ( exp( xlrg )/( 4.D0*pi*rhow )*3.D0 )**(1.D0/3.D0 )*1.D+06

   if( rlrg <=50.D0 ) then
     fckrn = 9.44D+03*( vlrg*vlrg + vsml*vsml )
   else
     fckrn = 5.78D-03*( vlrg+vsml )
   end if
  else if( kphase == 2 ) then 
   fckrn = 1.5D0*( exp(x) +exp(y) )
  end if 

  return

  end function fckrn
  !---------------------------------------------------------------------------------------
  subroutine getvt

  integer :: myu, n

  do myu = 1, nspc
  do n = 1, nbin
    vt( myu,n ) = max( fvterm( myu,xctr( n ) ), 0.D0 )
  end do
  end do

  end subroutine getvt
  !---------------------------------------------------------------------------------------
  function fvterm( myu,x )

  use mod_spline, only: getknot, getmatrx, getcoef, fspline

  integer, intent(in) :: myu
  real(8), intent(in) :: x
  real(8) :: fvterm

  real(8) :: qknt( ndat+kdeg ), elm( ndat,ndat ), coef( ndat )

  call getknot                        &
         ( ndat, kdeg, xmss( myu,: ), & !--- in
           qknt                       ) !--- out

  call getmatrx                             &
         ( ndat, kdeg, qknt, xmss( myu,: ), & !--- in
           elm                              ) !--- out

  call getcoef                             &
         ( ndat, kdeg, elm, vtrm( myu,: ), & !--- in
           coef                            ) !--- out

  fvterm = fspline ( ndat, kdeg, coef, qknt, xmss( myu,: ), x )

  end function fvterm
  !---------------------------------------------------------------------------------------
  subroutine getbr 

  integer :: myu, n

  do myu = 1, nspc
  do n = 1, nbin
    br( myu,n ) = fbulkrad( myu, xctr( n ) ) 
  end do
  end do

  end subroutine getbr
  !---------------------------------------------------------------------------------------
  function fbulkrad( myu,x )

  use mod_spline, only: getknot, getmatrx, getcoef, fspline

  integer, intent(in) :: myu
  real(8), intent(in) :: x
  real(8) :: fbulkrad

  real(8) :: qknt( ndat+kdeg ), elm( ndat,ndat ), coef( ndat )

  call getknot                        &
         ( ndat, kdeg, xmss( myu,: ), & !--- in
           qknt                       ) !--- out

  call getmatrx                             &
         ( ndat, kdeg, qknt, xmss( myu,: ), & !--- in
           elm                              ) !--- out

  call getcoef                             &
         ( ndat, kdeg, elm, blkr( myu,: ), & !--- in
           coef                            ) !--- out

  fbulkrad = fspline ( ndat, kdeg, coef, qknt, xmss( myu,: ), x )

  end function fbulkrad
  !---------------------------------------------------------------------------------------
  subroutine paraout 

  integer :: myu, nyu, i, j, n

  open ( 21, file = 'micpara.dat', form = 'formatted' )

  write( 21,* ) nspc, nbin

  ! grid parameter
  do n = 1, nbin
    write( 21,* ) n, xctr( n ), radc( n )
  end do
  do n = 1, nbin+1
    write( 21,* ) n, xbnd( n )
  end do
  write( 21,* ) dxmic

  ! capacity
  do myu = 1, nspc
    do n = 1, nbin
      write( 21,* ) myu, n, cctr( myu,n )
    end do
    do n = 1, nbin+1
      write( 21,* ) myu, n, cbnd( myu,n )
    end do
  end do

  ! collection kernel
  do myu = 1, nspc
  do nyu = 1, nspc
  do i = 1, nbin
  do j = 1, nbin
    write( 21,* ) myu, nyu, i, j, ck( myu,nyu,i,j )
  end do
  end do
  end do
  end do

  ! falling velocity
  do myu = 1, nspc
  do n = 1, nbin
    write( 21,* ) myu, n, vt( myu,n )
  end do
  end do

  ! bulk radius
  do myu = 1, nspc
  do n = 1, nbin
    write( 21,* ) myu, n, br( myu,n )
  end do
  end do

  if( kphase == 0 ) then
    write( 21,* ) "use hydrodynamic kernel"
  elseif( kphase == 1 ) then
    write( 21,* ) "use LONG's kernel"
  elseif( kphase == 2 ) then
    write( 21,* ) "use GOLOVIN's kernel"
  endif

  close ( 21 )

  end subroutine paraout
  !---------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
end program mkpara
