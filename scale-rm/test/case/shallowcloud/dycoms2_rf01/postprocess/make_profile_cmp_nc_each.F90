program convine
 implicit none

  include 'netcdf.inc'

  integer(4), parameter :: xproc=8, yproc=4
  integer(4), parameter :: nxp=3, nyp=4
  integer(4), parameter :: nx=nxp*xproc, ny=nyp*yproc, nz=276

  integer(4), parameter :: nst=30, nen=60
  integer(4), parameter :: nt=nen-nst+1

  integer(4), parameter :: nzhalo=2

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
  character(len=256) :: cfile
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
