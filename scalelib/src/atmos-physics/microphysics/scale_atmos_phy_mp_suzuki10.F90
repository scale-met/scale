!-------------------------------------------------------------------------------
!> Module Spectran Bin Microphysical Module in SCALE-LES ver. 3
!!
!! @par Description:
!!      This module contains subroutines for the Spectral Bin Model
!!
!! - Reference
!!  - Suzuki et al., 2006
!!    Correlation Pattern between Effective Radius and Optical Thickness of Water Clouds Simulated by a Spectral Bin Microphysics Cloud Model
!!    SOLA, 2: 116-119 doi:10.2151/sola.2006-030
!!  - Suzuki et al., 2010
!!    A Study of Microphysical Mechanisms for Correlation Patterns between Droplet Radius and Optical Thickness of Warm Clouds with a Spectral Bin
!!    J. Atmos. Sci., 67: 1126-1141
!!  - Sato et al., 2009
!!    Application of a Monte Carlo integration method to collision and coagulation growth processes of hydrometeors in a bin-type model
!!    J. Geophy. Res., 114: D09215, doi:10.1029/2008JD011247
!!
!! @author : Team SCALE
!!
!! @par History: Hbinw
!! @li  ver.0.00   2012-06-14 (Y.Sato) [new] Import from version 4.1 of original code
!! @li  ver.0.01   2012-09-14 (Y.Sato) [mod] add a stochastic method (Sato et al. 2009)
!! @li  ver.0.01   2013-02-12 (Y.Sato) [mod] modified for latest version
!! @li  ver.0.01   2013-12-26 (Y.Sato) [mod] mearge all version of Bin scheme
!<
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
module scale_atmos_phy_mp_suzuki10
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer_suzuki10
  use scale_const, only: &
     pi => CONST_PI, &
     CONST_CPdry, &
     CONST_CVdry, &
     CONST_DWATR, &
     CONST_GRAV, &
     CONST_Rvap, &
     CONST_Rdry, &
     CONST_LH0, &
     CONST_LHS0, &
     CONST_EMELT, &
     CONST_TEM00, &
     CONST_TMELT, &
     CONST_PSAT0, &
     CONST_PRE00, &
     CONST_CPovCV, &
     CONST_RovCP
  use scale_history, only: &
     HIST_in
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_suzuki10_setup
  public :: ATMOS_PHY_MP_suzuki10
  public :: ATMOS_PHY_MP_suzuki10_CloudFraction
  public :: ATMOS_PHY_MP_suzuki10_EffectiveRadius
  public :: ATMOS_PHY_MP_suzuki10_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]

# include "kernels.h"
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !--- Indeces
  integer, parameter :: il = 1
  integer, parameter :: ic = 2
  integer, parameter :: ip = 3
  integer, parameter :: id = 4
  integer, parameter :: iss= 5
  integer, parameter :: ig = 6
  integer, parameter :: ih = 7

  !--- bin information of hydrometeors
  real(RP) :: dxmic                !--- d( log(m) ) of hydrometeor bin
  real(RP), allocatable :: xctr( : )         !--- log( m ) value of bin center
  real(RP), allocatable :: xbnd( : )       !--- log( m ) value of bin boundary
  real(RP), allocatable :: radc( : )         !--- radius of hydrometeor at bin center [m]
  real(RP), allocatable :: cctr( :,: )       !--- capacitance of hydrometeor at bin center
  real(RP), allocatable :: cbnd( :,: )     !--- capacitance of hydrometeor at bin boundary
  real(RP), allocatable :: ck( :,:,:,: )  !-- collection kernel
  real(RP), allocatable :: vt( :,: )         !--- terminal velocity of hydrometeor [m/s]
  real(RP), allocatable :: br( :,: )         !--- bulk density of hydrometeor [kg/m^3]
  integer,  allocatable  :: ifrsl( :,: )
  !--- bin information of aerosol (not supported)
  real(RP), allocatable :: xactr( : )        !--- log( ma ) value of bin center
  real(RP), allocatable :: xabnd( : )      !--- log( ma ) value of bin boundary
  real(RP), allocatable :: rada( : )

  real(RP) :: dxaer                !--- d( log(ma) ) of aerosol bin
  real(RP) :: xasta
  real(RP) :: xaend
  real(RP) :: sfc_precp
  real(RP), allocatable, save :: velw( :,:,:,: )
  integer, private, save :: MP_NSTEP_SEDIMENTATION
  real(RP), private, save :: MP_RNSTEP_SEDIMENTATION
  real(RP), private, save :: MP_DTSEC_SEDIMENTATION
  integer, private, save :: ntmax_sedimentation= 1

  !--- constant for bin
  real(RP), parameter :: cldmin = 1.0E-10_RP, eps = 1.0E-30_RP
  real(RP), parameter :: OneovThird = 1.0_RP/3.0_RP, ThirdovForth = 3.0_RP/4.0_RP
  real(RP), parameter :: TwoovThird = 2.0_RP/3.0_RP
  !--- constant for aerosol
  real(RP) :: rhoa   = 2.25E+03_RP         ! density of aerosol ( NaCl )
  real(RP) :: emaer  = 58.0_RP             ! molecular weight of aerosol ( salt )
  real(RP) :: emwtr  = 18.0_RP             ! molecular weight of water
  real(RP) :: rasta  = 1.E-08_RP           ! minimum radius of aerosol (m)
  real(RP) :: raend  = 1.E-06_RP           ! maximum radius of aerosol (m)
  real(RP) :: r0a    = 1.E-07_RP           ! average radius of aerosol (m)
  logical :: flg_regeneration=.false.      ! flag regeneration of aerosol
  logical :: flg_nucl=.false.              ! flag nucleated cloud move into smallest bin
  logical :: flg_icenucl=.false.           ! flag ice nucleation
  logical :: flg_sf_aero =.false.          ! flag surface flux of aerosol
  integer, private, save :: rndm_flgp = 0  ! flag for sthastic integration for coll.-coag.
  logical, private, save :: doautoconversion = .true.  ! apply collision process ?
  logical, private, save :: doprecipitation  = .true.  ! apply sedimentation of hydrometeor ?
  logical, private, save :: donegative_fixer = .true.  ! apply negative fixer?

  real(RP), allocatable :: marate( : )               ! mass rate of each aerosol bin to total aerosol mass
  integer, private, save       :: K10_1, K10_2        ! scaling factor for 10m value (momentum)
  real(RP), private, save      :: R10M1, R10M2        ! scaling factor for 10m value (momentum)
  real(RP), private, save      :: R10H1, R10H2        ! scaling factor for 10m value (heat)
  real(RP), private, save      :: R10E1, R10E2        ! scaling factor for 10m value (tracer)

  character(11),parameter :: fname_micpara="micpara.dat"
  integer(4) :: fid_micpara

  !--- Use for stochastic method
  integer, allocatable :: blrg( :,: ), bsml( :,: )
  real(RP) :: wgtbin
  integer  :: mspc, mbin
  real(RP), private :: rndm(1,1,1)
  !--- use for model without aerosol
  real(RP), private :: c_ccn = 100.E+6_RP
  real(RP), private :: kappa = 0.462_RP
  !--- use for aerosol coupled model
  real(RP), private :: sigma = 7.5E-02_RP  ! water surface tension [ N/m2 ]
  real(RP), private :: vhfct = 2.0_RP    ! van't hoff factor

  !--- for mkpara
  integer, parameter :: ndat = 33, icemax = 3
  integer, parameter :: kdeg = 4, ldeg = 4, nspc_mk = 7
  real(DP) :: dxmic_mk
  real(DP), allocatable :: radc_mk( : ), xctr_mk( : ), xbnd_mk( : )
  real(DP), allocatable :: cctr_mk( :,: ), cbnd_mk( :,: )
  real(DP), allocatable :: ck_mk( :,:,:,: )
  real(DP), allocatable :: vt_mk( :,: )
  real(DP), allocatable :: br_mk( :,: )
  real(DP) :: xmss( nspc_mk,ndat ), zcap( nspc_mk,ndat ), vtrm( nspc_mk,ndat )
  real(DP) :: blkr( nspc_mk,ndat ), blkd( nspc_mk,ndat ), ykrn( nspc_mk,nspc_mk,ndat,ndat )

  real(DP) :: ywll( ndat,ndat ), ywli( ndat,ndat,icemax ), ywls( ndat,ndat )
  real(DP) :: ywlg( ndat,ndat ), ywlh( ndat,ndat )

  real(DP) :: ywil( ndat,ndat,icemax ), ywii( ndat,ndat,icemax,icemax )
  real(DP) :: ywis( ndat,ndat,icemax ), ywig( ndat,ndat,icemax )
  real(DP) :: ywih( ndat,ndat,icemax )

  real(DP) :: ywsl( ndat,ndat ), ywsi( ndat,ndat,icemax ), ywss( ndat,ndat )
  real(DP) :: ywsg( ndat,ndat ), ywsh( ndat,ndat )

  real(DP) :: ywgl( ndat,ndat ), ywgi( ndat,ndat,icemax ), ywgs( ndat,ndat )
  real(DP) :: ywgg( ndat,ndat ), ywgh( ndat,ndat )

  real(DP) :: ywhl( ndat,ndat ), ywhi( ndat,ndat,icemax ), ywhs( ndat,ndat )
  real(DP) :: ywhg( ndat,ndat ), ywhh( ndat,ndat )

  !----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_suzuki10_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_master,  &
       PRC_myrank
    use scale_const, only: &
       CONST_DWATR, &
       CONST_DICE
    use scale_comm, only: &
       COMM_datatype
    use scale_grid, only: &
       CDZ => GRID_CDZ, &
       CZ  => GRID_CZ,  &
       FZ  => GRID_FZ
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    implicit none

    character(len=H_SHORT), intent(in) :: MP_TYPE

    real(RP) :: ATMOS_PHY_MP_RHOA  !--- density of aerosol
    real(RP) :: ATMOS_PHY_MP_EMAER !--- moleculer weight of aerosol
    real(RP) :: ATMOS_PHY_MP_RAMIN !--- minimum radius of aerosol (um)
    real(RP) :: ATMOS_PHY_MP_RAMAX !--- maximum radius of aerosol (um)
    real(RP) :: ATMOS_PHY_MP_R0A   !--- maximum radius of aerosol (um)
    logical :: ATMOS_PHY_MP_FLAG_REGENE  !--- flag of regeneration
    logical :: ATMOS_PHY_MP_FLAG_NUCLEAT !--- flag of regeneration
    logical :: ATMOS_PHY_MP_FLAG_ICENUCLEAT !--- flag of regeneration
    logical :: ATMOS_PHY_MP_FLAG_SFAERO  !--- flag of surface flux of aeorol
    integer :: ATMOS_PHY_MP_RNDM_FLGP  !--- flag of surface flux of aeorol
    integer :: ATMOS_PHY_MP_RNDM_MSPC
    integer :: ATMOS_PHY_MP_RNDM_MBIN

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       ATMOS_PHY_MP_RHOA,  &
       ATMOS_PHY_MP_EMAER, &
       ATMOS_PHY_MP_RAMIN, &
       ATMOS_PHY_MP_RAMAX, &
       ATMOS_PHY_MP_FLAG_REGENE,  &
       ATMOS_PHY_MP_FLAG_NUCLEAT, &
       ATMOS_PHY_MP_FLAG_ICENUCLEAT, &
       ATMOS_PHY_MP_FLAG_SFAERO,  &
       ATMOS_PHY_MP_R0A,   &
       ATMOS_PHY_MP_RNDM_FLGP, &
       ATMOS_PHY_MP_RNDM_MSPC, &
       ATMOS_PHY_MP_RNDM_MBIN, &
       doautoconversion, &
       doprecipitation, &
       donegative_fixer, &
       c_ccn, kappa, &
       sigma, vhfct

    integer :: nnspc, nnbin
    integer :: nn, mm, mmyu, nnyu
    integer :: myu, nyu, i, j, k, n, ierr
    integer :: COMM_world
    !---------------------------------------------------------------------------

    !--- allocation
    allocate( xctr( nbin ) )
    allocate( xbnd( nbin+1 ) )
    allocate( radc( nbin ) )
    allocate( cctr( nspc_mk,nbin ) )
    allocate( cbnd( nspc_mk,nbin+1 ) )
    allocate( ck( nspc_mk,nspc_mk,nbin,nbin ) )
    allocate( vt( nspc_mk,nbin ) )
    allocate( br( nspc_mk,nbin ) )
    allocate( ifrsl( nspc_mk,nspc_mk ) )
    if( nccn /= 0 ) then
      allocate( xactr( nccn ) )
      allocate( xabnd( nccn+1 ) )
      allocate( rada( nccn ) )
    endif

    mbin = nbin/2
    mspc = nspc_mk*nspc_mk

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for SBM (warm cloud)'

    if ( MP_TYPE .ne. 'SUZUKI10' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_MP_TYPE is not SUZUKI10. Check!'
       call PRC_MPIstop
    end if

    ATMOS_PHY_MP_RHOA = rhoa
    ATMOS_PHY_MP_EMAER = emaer
    ATMOS_PHY_MP_RAMIN = rasta
    ATMOS_PHY_MP_RAMAX = raend
    ATMOS_PHY_MP_R0A = r0a
    ATMOS_PHY_MP_FLAG_REGENE = flg_regeneration
    ATMOS_PHY_MP_FLAG_NUCLEAT = flg_nucl
    ATMOS_PHY_MP_FLAG_ICENUCLEAT = flg_icenucl
    ATMOS_PHY_MP_FLAG_SFAERO = flg_sf_aero
    ATMOS_PHY_MP_RNDM_FLGP = rndm_flgp
    ATMOS_PHY_MP_RNDM_MSPC = mspc
    ATMOS_PHY_MP_RNDM_MBIN = mbin

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP, Check!'
     call PRC_MPIstop
    end if
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    if ( nspc /= 1 .and. nspc /= 7 ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx nspc should be set as 1(warm rain) or 7(mixed phase) check!'
       call PRC_MPIstop
    end if

    rhoa = ATMOS_PHY_MP_RHOA
    emaer = ATMOS_PHY_MP_EMAER
    rasta = ATMOS_PHY_MP_RAMIN
    raend = ATMOS_PHY_MP_RAMAX
    r0a   = ATMOS_PHY_MP_R0A
    flg_regeneration = ATMOS_PHY_MP_FLAG_REGENE
    flg_nucl = ATMOS_PHY_MP_FLAG_NUCLEAT
    flg_icenucl = ATMOS_PHY_MP_FLAG_ICENUCLEAT
    flg_sf_aero = ATMOS_PHY_MP_FLAG_SFAERO
    rndm_flgp = ATMOS_PHY_MP_RNDM_FLGP
    mspc = ATMOS_PHY_MP_RNDM_MSPC
    mbin = ATMOS_PHY_MP_RNDM_MBIN

    !--- read micpara.dat (microphysical parameter) and broad cast
    if( PRC_myrank == PRC_master ) then

      fid_micpara = IO_get_available_fid()
      !--- open parameter of cloud microphysics
      open ( fid_micpara, file = fname_micpara, form = 'formatted', status = 'old', iostat=ierr )

      !--- micpara.dat does not exist
      if( ierr == 0 ) then

        read( fid_micpara,* ) nnspc, nnbin

        if( nnbin /= nbin ) then
           if ( IO_L ) write(IO_FID_LOG,*) 'xxx nbin in inc_tracer and nbin in micpara.dat is different check!'
           call PRC_MPIstop
        end if

        ! grid parameter
        if( IO_L ) write(IO_FID_LOG,*)  '*** Radius of cloud ****'
        do n = 1, nbin
          read( fid_micpara,* ) nn, xctr( n ), radc( n )
          if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)') &
                    "Radius of ", n, "th cloud bin (bin center)= ", radc( n ) , "[m]"
        end do
        do n = 1, nbin+1
          read( fid_micpara,* ) nn, xbnd( n )
        end do
        read( fid_micpara,* ) dxmic
        if( IO_L ) write(IO_FID_LOG,*)  '*** Width of Cloud SDF= ', dxmic

        ! capacity
        do myu = 1, nspc_mk
         do n = 1, nbin
          read( fid_micpara,* ) mmyu, nn, cctr( myu,n )
         end do
         do n = 1, nbin+1
          read( fid_micpara,* ) mmyu, nn, cbnd( myu,n )
         end do
        end do

        ! collection kernel
        do myu = 1, nspc_mk
         do nyu = 1, nspc_mk
          do i = 1, nbin
           do j = 1, nbin
            read( fid_micpara,* ) mmyu, nnyu, mm, nn, ck( myu,nyu,i,j )
           enddo
          enddo
         enddo
        enddo

        ! terminal velocity
        do myu = 1, nspc_mk
         do n = 1, nbin
          read( fid_micpara,* ) mmyu, nn, vt( myu,n )
         enddo
        enddo

        ! bulk density
        do myu = 1, nspc_mk
         do n = 1, nbin
          read( fid_micpara,* ) mmyu, nn, br( myu,n )
         enddo
        enddo

        close ( fid_micpara )

      !--- micpara.dat does not exist
      else

        if ( IO_L ) write(IO_FID_LOG,*) 'micpara.dat is created'
        call mkpara

        fid_micpara = IO_get_available_fid()
        !--- open parameter of cloud microphysics
        open ( fid_micpara, file = fname_micpara, form = 'formatted', status = 'old', iostat=ierr )

        read( fid_micpara,* ) nnspc, nnbin

        if( nnbin /= nbin ) then
           if ( IO_L ) write(IO_FID_LOG,*) 'xxx nbin in inc_tracer and nbin in micpara.dat is different check!'
           call PRC_MPIstop
        end if

        ! grid parameter
        if( IO_L ) write(IO_FID_LOG,*)  '*** Radius of cloud ****'
        do n = 1, nbin
          read( fid_micpara,* ) nn, xctr( n ), radc( n )
          if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)') &
                    "Radius of ", n, "th cloud bin (bin center)= ", radc( n ) , "[m]"
        end do
        do n = 1, nbin+1
          read( fid_micpara,* ) nn, xbnd( n )
        end do
        read( fid_micpara,* ) dxmic
        if( IO_L ) write(IO_FID_LOG,*)  '*** Width of Cloud SDF= ', dxmic

        ! capacity
        do myu = 1, nspc_mk
         do n = 1, nbin
          read( fid_micpara,* ) mmyu, nn, cctr( myu,n )
         end do
         do n = 1, nbin+1
          read( fid_micpara,* ) mmyu, nn, cbnd( myu,n )
         end do
        end do

        ! collection kernel
        do myu = 1, nspc_mk
         do nyu = 1, nspc_mk
          do i = 1, nbin
           do j = 1, nbin
            read( fid_micpara,* ) mmyu, nnyu, mm, nn, ck( myu,nyu,i,j )
           enddo
          enddo
         enddo
        enddo

        ! terminal velocity
        do myu = 1, nspc_mk
         do n = 1, nbin
          read( fid_micpara,* ) mmyu, nn, vt( myu,n )
         enddo
        enddo

        ! bulk density
        do myu = 1, nspc_mk
         do n = 1, nbin
          read( fid_micpara,* ) mmyu, nn, br( myu,n )
         enddo
        enddo

        close ( fid_micpara )
     
      endif

    endif

    n = nspc_mk*nspc_mk*nbin*nbin
    COMM_world = MPI_COMM_WORLD
    call MPI_BCAST( xctr, nbin,             COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( dxmic,1,                COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( xbnd, nbin+1,           COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( cctr, nbin*nspc_mk,     COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( cbnd, (nbin+1)*nspc_mk, COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( ck,   n,                COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( br, nbin*nspc_mk,       COMM_datatype, PRC_master, COMM_world, ierr )
    call MPI_BCAST( vt, nbin*nspc_mk,       COMM_datatype, PRC_master, COMM_world, ierr )

    !--- aerosol ( CCN ) (not supported)
    if( nccn /= 0 ) then

    xasta = log( rhoa*4.0_RP/3.0_RP*pi * ( rasta )**3 )
    xaend = log( rhoa*4.0_RP/3.0_RP*pi * ( raend )**3 )

    dxaer = ( xaend-xasta )/nccn

    do n = 1, nccn+1
     xabnd( n ) = xasta + dxaer*( n-1 )
    enddo
    do n = 1, nccn
     xactr( n ) = ( xabnd( n )+xabnd( n+1 ) )*0.50_RP
     rada( n )  = ( exp( xactr( n ) )*ThirdovForth/pi/rhoa )**( OneovThird )
     if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)') &
          "Radius of ", n, "th aerosol bin (bin center)= ", rada( n ) , "[m]"
    enddo

    if( flg_sf_aero ) then
     if ( CZ(KS) >= 10.0_RP ) then
          R10M1 = 10.0_RP / CZ(KS) * 0.50_RP ! scale with height
          R10M2 = 10.0_RP / CZ(KS) * 0.50_RP ! scale with height
          R10H1 = 1.0_RP * 0.50_RP
          R10H2 = 1.0_RP * 0.0_RP
          R10E1 = 1.0_RP * 0.50_RP
          R10E2 = 1.0_RP * 0.50_RP
          K10_1 = KS
          K10_2 = KS
     else
       k = 1
       do while ( CZ(k) < 10.0_RP )
          k = k + 1
          K10_1 = k
          K10_2 = k + 1
          R10M1 = ( CZ(k+1) - 10.0_RP ) / CDZ(k)
          R10M2 = ( 10.0_RP   - CZ(k) ) / CDZ(k)
          R10H1 = ( CZ(k+1) - 10.0_RP ) / CDZ(k)
          R10H2 = ( 10.0_RP   - CZ(k) ) / CDZ(k)
          R10E1 = ( CZ(k+1) - 10.0_RP ) / CDZ(k)
          R10E2 = ( 10.0_RP   - CZ(k) ) / CDZ(k)
       enddo
     endif
    endif

    endif

    ATMOS_PHY_MP_DENS(I_mp_QC)  = CONST_DWATR
    ATMOS_PHY_MP_DENS(I_mp_QCL) = CONST_DICE
    ATMOS_PHY_MP_DENS(I_mp_QD)  = CONST_DICE
    ATMOS_PHY_MP_DENS(I_mp_QS)  = CONST_DICE
    ATMOS_PHY_MP_DENS(I_mp_QG)  = CONST_DICE
    ATMOS_PHY_MP_DENS(I_mp_QH)  = CONST_DICE

    !--- random number setup for stochastic method
    if( rndm_flgp > 0 ) then
     call random_setup( IA*JA*KA )
    endif

    allocate( velw(KA,IA,JA,QA) )
    velw(:,:,:,:) = 0.0_RP
    do myu = 1, nspc
    do n = 1, nbin
      velw(:,:,:,I_QV+(myu-1)*nbin+n) = -vt( myu,n )
    enddo
    enddo

    MP_NSTEP_SEDIMENTATION  = ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(ntmax_sedimentation,kind=RP)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

    return

  end subroutine ATMOS_PHY_MP_suzuki10_setup
  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_suzuki10( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP, &
       TIME_NOWDAYSEC
    use scale_grid, only: &
       GRID_CZ,  &
       GRID_FZ,  &
       GRID_CDZ, &
       GRID_FDZ
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_thermodyn, only: &
       AQ_CV, &
       AQ_CP
    use scale_atmos_phy_mp_common, only: &
       MP_precipitation => ATMOS_PHY_MP_PRECIPITATION
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer => ATMOS_PHY_MP_negative_fixer
    use scale_atmos_saturation, only: &
       pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
       pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)

    real(RP) :: dz (KA)
    real(RP) :: dzh(KA)
    real(RP) :: dt
    real(RP) :: ct

    real(RP) :: temp_mp  (KA,IA,JA)
    real(RP) :: pres_mp  (KA,IA,JA)
    real(RP) :: gdgc (KA,IA,JA,nspc,nbin) !-- SDF of hydrometeors [kg/m^3/unit ln(r)]
    real(RP), allocatable :: gdga (:,:,:,:) !-- SDF of aerosol (not supported)
    real(RP) :: qv_mp    (KA,IA,JA)      !-- Qv [kg/kg]
    real(RP) :: wfall( KA )

    real(RP) :: ssliq(KA,IA,JA)
    real(RP) :: ssice(KA,IA,JA)
    real(RP) :: sum1(KA,IA,JA)
    real(RP) :: sum2
    integer :: m, n, k, i, j, iq, countbin
    logical, save :: ofirst_sdfa = .true.

    real(RP) :: VELX(IA,JA)
    real(RP) :: VELY(IA,JA)
    real(RP) :: SFLX_AERO(IA,JA,nccn)
    real(RP) :: Uabs, bparam
    real(RP) :: AMR(KA,IA,JA)

    real(RP) :: rhogq(KA,IA,JA,QA)
    real(RP) :: rrhog(KA,IA,JA)
    real(RP) :: q(KA,IA,JA,QA)
    real(RP) :: qd(KA,IA,JA)
    real(RP) :: cva(KA,IA,JA)
    real(RP) :: cpa(KA,IA,JA)
    real(RP) :: rhoge(KA,IA,JA)
    real(RP) :: th(KA,IA,JA)
    real(RP) :: Rmoist

    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    real(RP) :: flux_rain (KA,IA,JA)
    real(RP) :: flux_snow (KA,IA,JA)
    real(RP) :: flux_prec (IA,JA)
    integer  :: step
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("MICROPHYSICS")
#endif

     if( donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
       QTRC(:,:,:,QQE+1:QA) = max( QTRC(:,:,:,QQE+1:QA),0.0_RP )
     endif


    if( nspc == 1 ) then
     if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(SBM Liquid water only)'
    elseif( nspc == 7 ) then
     if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(SBM Mixed phase)'
    endif

    if( flg_sf_aero ) then
     do j = JS-2, JE+2
     do i = IS-2, IE+1
       VELX(i,j) = MOMX(K10_1,i,j) / ( DENS(K10_1,i+1,j)+DENS(K10_1,i,j) ) * R10M1 &
                 + MOMX(K10_2,i,j) / ( DENS(K10_2,i+1,j)+DENS(K10_2,i,j) ) * R10M2
     enddo
     enddo

     do j = JS-2, JE+1
     do i = IS-2, IE+2
       VELY(i,j) = MOMY(K10_1,i,j) / ( DENS(K10_1,i,j+1)+DENS(K10_1,i,j) ) * R10M1 &
                 + MOMY(K10_2,i,j) / ( DENS(K10_2,i,j+1)+DENS(K10_2,i,j) ) * R10M2
     enddo
     enddo
    end if

    call PROF_rapstart('MPX ijkconvert')
    dz (:) = GRID_CDZ(:)
    dzh(1) = GRID_FDZ(1)
    dzh(2:KA) = GRID_FDZ(1:KA-1)

    dt = TIME_DTSEC_ATMOS_PHY_MP
    ct = TIME_NOWDAYSEC

    gdgc(:,:,:,:,:) = 0.0_RP
    pres_mp(:,:,:) = 0.0_RP
    temp_mp(:,:,:) = 0.0_RP
    qv_mp(:,:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE

       do k = KS-1, KE+1
          rrhog(k,i,j) = 1.0_RP / DENS(k,i,j)
          q(k,i,j,1:QA) = QTRC(k,i,j,1:QA)
       enddo
       do k = KS, KE
          th(k,i,j) = RHOT(k,i,j) * rrhog(k,i,j)
       enddo
       do k = KS, KE
          CALC_QDRY( qd(k,i,j), q, k, i, j, iq )
       enddo
       do k = KS, KE
          CALC_CV( cva(k,i,j), qd(k,i,j), q, k, i, j, iq, CONST_CVdry, AQ_CV )
       enddo
       do k = KS, KE
          CALC_R( Rmoist, q(k,i,j,I_QV), qd(k,i,j), CONST_Rdry, CONST_Rvap )
          cpa(k,i,j) = cva(k,i,j) + Rmoist
          CALC_PRE( pres_mp(k,i,j), DENS(k,i,j), th(k,i,j), Rmoist, cpa(k,i,j), CONST_PRE00 )
          temp_mp(k,i,j) = pres_mp(k,i,j) / ( DENS(k,i,j) * Rmoist )
          qv_mp(k,i,j) = QTRC(k,i,j,I_QV)

          countbin = 1
          do m = 1, nspc
          do n = 1, nbin
           gdgc( k,i,j,m,n ) = &
               QTRC(k,i,j,countbin+I_QV)*DENS(k,i,j)/dxmic
           countbin = countbin + 1
          end do
          end do
       enddo

    enddo
    enddo

    if( nccn /= 0 ) then

    allocate( gdga (KA,IA,JA,nccn) )
    gdga(:,:,:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          do n = 1, nccn
           gdga( k,i,j,n ) = &
               QTRC(k,i,j,n+nbin+I_QV)*DENS(k,i,j)/dxaer
          end do
          !--- store initial SDF of aerosol
          if( ofirst_sdfa ) then
           allocate( marate( nccn ) )
           sum2 = 0.0_RP
           do n = 1, nccn
             marate( n ) = gdga(k,i,j,n)/exp( xactr( n ) )
             sum2 = sum2 + gdga(k,i,j,n)/exp( xactr( n ) )
           end do
           if( sum2 /= 0.0_RP ) then
            marate( 1:nccn ) = marate( 1:nccn )/sum2
            ofirst_sdfa = .false.
           end if
          end if
       enddo

    enddo
    enddo

    endif


    call PROF_rapend  ('MPX ijkconvert')

    !--- Calculate Super saturation
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE

      call pres2qsat_liq( ssliq(k,i,j),temp_mp(k,i,j),pres_mp(k,i,j) )
      call pres2qsat_ice( ssice(k,i,j),temp_mp(k,i,j),pres_mp(k,i,j) )
      ssliq(k,i,j) = qv_mp(k,i,j)/ssliq(k,i,j)-1.0_RP
      ssice(k,i,j) = qv_mp(k,i,j)/ssice(k,i,j)-1.0_RP

      sum1(k,i,j) = 0.0_RP
      do m = 1, nspc
      do n = 1, nbin
        sum1(k,i,j) = sum1(k,i,j) + gdgc(k,i,j,m,n)*dxmic
      end do
      end do

    enddo
    enddo
    enddo

    !--- Time evolution of SDF
    if( nccn /= 0 ) then    ! with aerosol tracer

     if( nspc == 1 ) then   ! warm rain only with aerosol tracer

     do k = KS, KE
     do j = JS, JE
     do i = IS, IE

      if( ssliq(k,i,j) > 0.0_RP .or. sum1(k,i,j) > cldmin ) then
       call mp_abinw_evolve                    &
            ( pres_mp(k,i,j), DENS(k,i,j), dt, &  !--- in
              gdgc(k,i,j,1,1:nbin),            &  !--- inout
              gdga(k,i,j,1:nccn),              &  !--- inout  for aerosol tracer
              qv_mp(k,i,j), temp_mp(k,i,j)     )  !--- inout
      end if

     end do
     end do
     end do

     elseif( nspc > 1 ) then  ! mixed phase rain with aerosol tracer

     do k = KS, KE
     do j = JS, JE
     do i = IS, IE

      if( ssliq(k,i,j) > 0.0_RP .or. sum1(k,i,j) > cldmin ) then
       call mp_abinf_evolve                    &
            ( pres_mp(k,i,j), DENS(k,i,j), dt, &  !--- in
              gdgc(k,i,j,1:nspc,1:nbin),       &  !--- inout
              gdga(k,i,j,1:nccn),              &  !--- inout  for aerosol tracer
              qv_mp(k,i,j), temp_mp(k,i,j)     )  !--- inout
      end if

     end do
     end do
     end do

     endif

    elseif( nccn == 0 ) then  ! without aerosol tracer

     if( nspc == 1 ) then

     do k = KS, KE
     do j = JS, JE
     do i = IS, IE

      if( ssliq(k,i,j) > 0.0_RP .or. sum1(k,i,j) > cldmin ) then
       call mp_binw_evolve                     &
            ( pres_mp(k,i,j), DENS(k,i,j), dt, &  !--- in
              gdgc(k,i,j,1,1:nbin),            &  !--- inout
              qv_mp(k,i,j), temp_mp(k,i,j)     )  !--- inout
      end if

     end do
     end do
     end do

     elseif( nspc > 1 ) then

     do k = KS, KE
     do j = JS, JE
     do i = IS, IE

      if( ssliq(k,i,j) > 0.0_RP .or. ssice(k,i,j) > 0.0_RP .or. sum1(k,i,j) > cldmin ) then
       call mp_binf_evolve                     &
            ( pres_mp(k,i,j), DENS(k,i,j), dt, &  !--- in
              gdgc(k,i,j,1:nspc,1:nbin),       &  !--- inout
              qv_mp(k,i,j), temp_mp(k,i,j)     )  !--- inout
      end if

     end do
     end do
     end do

     endif

    endif

    !--- SURFACE FLUX by Monahan et al. (1986)
    if( flg_sf_aero .and. nccn /= 0 ) then
     do j = JS-1, JE
     do i = IS-1, IE
       Uabs = sqrt(  ( ( VELX(i,j) + VELX(i-1,j  ) ) * 0.50_RP )**2 &
                   + ( ( VELY(i,j) + VELY(i  ,j-1) ) * 0.50_RP )**2 )
       do n = 1, nccn
        if( rada( n ) <= 2.0E-5_RP .and. rada( n ) >= 3.0E-7_RP ) then
         bparam = ( 0.38_RP - log( rada( n ) ) )/0.65_RP
         SFLX_AERO(i,j,n) = 1.373_RP * Uabs**( 3.41_RP ) * rada( n )**( -3.0_RP ) &
                          * ( 1.0_RP + 0.057_RP * rada( n )**( 1.05_RP ) ) &
                          * 10.0_RP**( 1.19_RP * exp( -bparam*bparam ) )
         ! convert from [#/m^2/um/s] -> [kg/m^3/unit log (m)]
         SFLX_AERO(i,j,n) = SFLX_AERO(i,j,n) / DENS(KS,i,j) &
                          / GRID_CDZ(KS) * rada( n ) / 3.0_RP * dt * exp( xactr( n ) )
         gdga(KS,i,j,n) = gdga(KS,i,j,n)+SFLX_AERO(i,j,n)/dxaer
        end if
       end do
     end do
     end do
    end if

    call PROF_rapstart('MPX ijkconvert')
    AMR(:,:,:) = 0.0_RP
    do j = JS, JE
     do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_QV) = qv_mp(k,i,j)
          countbin = 1
          do m = 1, nspc
          do n = 1, nbin
            countbin = countbin+1
            QTRC(k,i,j,countbin) = gdgc(k,i,j,m,n)/DENS(k,i,j)*dxmic
            if( QTRC(k,i,j,n+I_QV) <= eps ) then
              QTRC(k,i,j,n+I_QV) = 0.0_RP
            end if
          end do
          end do
       enddo

       do k = KS, KE
          CALC_QDRY( qd(k,i,j), q, k, i, j, iq )
       enddo
       do k  = KS, KE
          CALC_CP( cpa(k,i,j), qd(k,i,j), q, k, i, j, iq, CONST_CPdry, AQ_CP )
          CALC_R( Rmoist, QTRC(k,i,j,I_QV), qd(k,i,j), CONST_Rdry, CONST_Rvap )
          RHOT(k,i,j) = temp_mp(k,i,j) * ( CONST_PRE00 / pres_mp(k,i,j) )**(Rmoist/cpa(k,i,j)) &
               * DENS(k,i,j)
       enddo

     enddo
    enddo

    if( nccn /= 0 ) then
    do j = JS, JE
     do i = IS, IE
       do k = KS, KE
          do n = 1, nccn
            QTRC(k,i,j,n+I_QV+nbin*nspc)=gdga(k,i,j,n)/DENS(k,i,j)*dxaer
            AMR(k,i,j) = AMR(k,i,j) + QTRC(k,i,j,n+I_QV+nbin)
          end do
          do n = 1, QA
           q(k,i,j,n) = QTRC(k,i,j,n)
          enddo
       enddo

     enddo
    enddo

    call HIST_in( AMR(:,:,:),  'aerosol', 'aerosol mass', 'kg/m^3', dt)

    deallocate(gdga)

    endif

    !--- gravitational falling
    if ( doprecipitation ) then
    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       flux_rain(k,i,j) = 0.0_RP
       flux_snow(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do step = 1, MP_NSTEP_SEDIMENTATION

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          th(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       do k = KS, KE
          CALC_QDRY( qd(k,i,j), QTRC, k, i, j, iq )
       enddo
       do k = KS, KE
          CALC_CV( cva(k,i,j), qd(k,i,j), QTRC, k, i, j, iq, CONST_CVdry, AQ_CV )
       enddo
       do k = KS, KE
          CALC_R( Rmoist, QTRC(k,i,j,I_QV), qd(k,i,j), CONST_Rdry, CONST_Rvap )
          cpa(k,i,j) = cva(k,i,j) + Rmoist
          CALC_PRE( pres_mp(k,i,j), DENS(k,i,j), th(k,i,j), Rmoist, cpa(k,i,j), CONST_PRE00 )
          temp_mp(k,i,j) = pres_mp(k,i,j) / ( DENS(k,i,j) * Rmoist )
       enddo
       do k = KS, KE
          rhoge(k,i,j)  = DENS(k,i,j) * temp_mp(k,i,j) * cva(k,i,j)
       enddo
       enddo
       enddo

       call MP_precipitation( &
            wflux_rain, wflux_snow, &
            DENS, MOMZ, MOMX, MOMY, &
            rhoge, QTRC, &
            velw, temp_mp, &
            MP_DTSEC_SEDIMENTATION )

       do j = JS, JE
       do i = IS, IE
          do k = KS-1, KE
             flux_rain(k,i,j) = flux_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
             flux_snow(k,i,j) = flux_snow(k,i,j) + wflux_snow(k,i,j) * MP_RNSTEP_SEDIMENTATION
          enddo
          flux_prec(i,j) = flux_rain(KS-1,i,j) + flux_snow(KS-1,i,j)
       enddo
       enddo

    enddo

    endif

    if( donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
       QTRC(:,:,:,QQE+1:QA) = max( QTRC(:,:,:,QQE+1:QA),0.0_RP )
    endif

    call HIST_in( flux_rain(KS-1,:,:), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:), 'SNOW', 'surface snow rate', 'kg/m2/s', dt)
    call HIST_in( flux_prec(:,:),      'PREC', 'surface precipitaion rate', 'kg/m2/s', dt)

    call PROF_rapend  ('MPX ijkconvert')

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_vars8( MOMZ(:,:,:), 2 )
    call COMM_vars8( MOMX(:,:,:), 3 )
    call COMM_vars8( MOMY(:,:,:), 4 )
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( DENS(:,:,:), 1 )
    call COMM_wait ( MOMZ(:,:,:), 2 )
    call COMM_wait ( MOMX(:,:,:), 3 )
    call COMM_wait ( MOMY(:,:,:), 4 )
    call COMM_wait ( RHOT(:,:,:), 5 )

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("MICROPHYSICS")
#endif

    return
  end subroutine ATMOS_PHY_MP_suzuki10
  !-----------------------------------------------------------------------------
  subroutine mp_abinw_evolve        &
      ( pres, dens,                 & !--- in
        dtime,                      & !--- in
        gc,                         & !--- inout
        ga,                         & !--- inout
        qvap, temp                  ) !--- inout

  real(RP), intent(in) :: pres   !  pressure
  real(RP), intent(in) :: dens   !  density of dry air
  real(RP), intent(in) :: dtime  !  time interval

  real(RP), intent(inout) :: gc( nspc,nbin )
  real(RP), intent(inout) :: ga( nccn )  !--- aerosol SDF (not supported)
  real(RP), intent(inout) :: qvap  !  specific humidity
  real(RP), intent(inout) :: temp  !  temperature
  integer :: n
  !
  !
  !--- nucleat from aerosol
  call nucleata                 &
         ( dens, pres, dtime,   & !--- in
           gc(il,1:nbin), ga, qvap, temp   ) !--- inout

  !--- condensation / evaporation
  call cndevpsbla               &
         ( dtime,               & !--- in
           dens, pres,          & !--- in
           gc, ga, qvap, temp   ) !--- inout

  !--- collision-coagulation
  if( doautoconversion ) then
   call  collmain               &
          ( dtime,              & !--- in
            gc                  ) !--- inout
  endif

  return

  end subroutine mp_abinw_evolve
  !-----------------------------------------------------------------------------
  subroutine mp_binw_evolve        &
      ( pres, dens,                 & !--- in
        dtime,                      & !--- in
        gc,                         & !--- inout
        qvap, temp                  ) !--- inout

  real(RP), intent(in) :: pres   !  pressure
  real(RP), intent(in) :: dens   !  density of dry air
  real(RP), intent(in) :: dtime  !  time interval

  real(RP), intent(inout) :: gc( nspc,nbin )
  real(RP), intent(inout) :: qvap  !  specific humidity
  real(RP), intent(inout) :: temp  !  temperature
  integer :: n
  !
  !
  !--- nucleat
  call nucleat                        &
         ( dens, pres, dtime,         & !--- in
           gc(il,1:nbin), qvap, temp  ) !--- inout

  !--- condensation / evaporation
  call cndevpsbl                &
         ( dtime,               & !--- in
           dens, pres,          & !--- in
           gc, qvap, temp       ) !--- inout

  !--- collision-coagulation
  if( doautoconversion ) then
   call  collmain               &
          ( dtime,              & !--- in
            gc                  ) !--- inout
  endif

  return

  end subroutine mp_binw_evolve
  !-----------------------------------------------------------------------------
  subroutine mp_binf_evolve         &
      ( pres, dens,                 & !--- in
        dtime,                      & !--- in
        gc,                         & !--- inout
        qvap, temp                  ) !--- inout

  use scale_const, only: &
     CONST_TEM00
  real(RP), intent(in) :: pres   !  pressure
  real(RP), intent(in) :: dens   !  density
  real(RP), intent(in) :: dtime  !  time interval

  real(RP), intent(inout) :: gc( nspc,nbin )
  real(RP), intent(inout) :: qvap  !  specific humidity
  real(RP), intent(inout) :: temp  !  temperature
  !
  !
  !--- nucleat
  call nucleat                  &
         ( dens, pres, dtime,   & !--- in
           gc(il,1:nbin), qvap, temp )   !--- inout

  !--- freezing / melting
  if ( temp < CONST_TEM00 ) then
    call freezing               &
           ( dtime, dens,       & !--- in
             gc, temp           ) !--- inout
!    call ice_nucleat            &
!           ( dtime, dens, pres, & !--- in
!             gc, qvap, temp     ) !--- inout
  else
    call melting                &
           ( dtime, dens,       & !--- in
             gc, temp           ) !--- inout
  end if

  !--- condensation / evaporation
  call cndevpsbl                &
         ( dtime,               & !--- in
           dens, pres,          & !--- in
           gc, qvap, temp       ) !--- inout

  if( doautoconversion ) then
  !--- collision-coagulation
   call  collmainf              &
           ( temp, dtime,       & !--- in
             gc                 ) !--- inout
  end if

  return

  end subroutine mp_binf_evolve
  !-----------------------------------------------------------------------------
  subroutine mp_abinf_evolve        &
      ( pres, dens,                 & !--- in
        dtime,                      & !--- in
        gc,                         & !--- inout
        ga,                         & !--- inout
        qvap, temp                  ) !--- inout

  use scale_const, only: &
     CONST_TEM00
  real(RP), intent(in) :: pres   !  pressure
  real(RP), intent(in) :: dens   !  density
  real(RP), intent(in) :: dtime  !  time interval

  real(RP), intent(inout) :: gc( nspc,nbin )
  real(RP), intent(inout) :: ga( nccn )  !--- aerosol SDF (not supported)
  real(RP), intent(inout) :: qvap  !  specific humidity
  real(RP), intent(inout) :: temp  !  temperature
  !
  !
  !--- nucleat
  call nucleata                 &
         ( dens, pres, dtime,   & !--- in
           gc(il,1:nbin), ga, qvap, temp )   !--- inout

  !--- freezing / melting
  if ( temp < CONST_TEM00 ) then
    call freezing               &
           ( dtime, dens,       & !--- in
             gc, temp           ) !--- inout
!    call ice_nucleat            &
!           ( dtime, dens, pres, & !--- in
!             gc, qvap, temp     ) !--- inout
  else
    call melting                &
           ( dtime, dens,       & !--- in
             gc, temp           ) !--- inout
  end if

  !--- condensation / evaporation
  call cndevpsbla               &
         ( dtime,               & !--- in
           dens, pres,          & !--- in
           gc, ga, qvap, temp   ) !--- inout

  if( doautoconversion ) then
  !--- collision-coagulation
   call  collmainf              &
           ( temp, dtime,       & !--- in
             gc                 ) !--- inout
  end if

  return

  end subroutine mp_abinf_evolve
  !-----------------------------------------------------------------------------
  subroutine nucleat        &
      ( dens, pres, dtime,  & !--- in
        gc, qvap, temp      ) !--- inout
  !
  !  liquid nucleation from aerosol particle
  !
  !
  use scale_const, only: &
     cp    => CONST_CPdry, &
     rhow  => CONST_DWATR, &
     qlevp => CONST_LH0, &
     rvap  => CONST_Rvap
  use scale_atmos_saturation, only: &
       pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
       pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
  real(RP), intent(in) :: dens   !  density  [ kg/m3 ]
  real(RP), intent(in) :: pres   !  pressure [ Pa ]
  real(RP), intent(in) :: dtime
  !
  real(RP), intent(inout) :: gc( nbin )  !  SDF ( hydrometeors )
  real(RP), intent(inout) :: qvap  !  specific humidity [ kg/kg ]
  real(RP), intent(inout) :: temp  !  temperature [ K ]
  !
  !--- local
  real(RP) :: ssliq, ssice, delcld
  real(RP) :: sumold, sumnew, acoef, bcoef, xcrit, rcrit
  real(RP) :: ractr, rcld, xcld, part, vdmp, dmp
  integer :: n, nc, ncrit
  integer, allocatable, save :: ncld( : )
  logical, save :: ofirst = .true.
  !
  real(RP) :: n_c, sumnum, gcn( nbin )
  !
  !--- supersaturation
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  if ( ssliq <= 0.0_RP ) return
  !--- use for aerosol coupled model
  !--- mass -> number
  do n = 1, nbin
    gcn( n ) = gc( n )/exp( xctr( n ) )
  end do

  sumnum = 0.0_RP
  do n = 1, nbin
    sumnum = sumnum + gcn( n )*dxmic
  enddo
  n_c = c_ccn * ( ssliq * 1.E+2_RP )**( kappa )
  if( n_c > sumnum ) then
    dmp = ( n_c - sumnum ) * exp( xctr( 1 ) )
    dmp = min( dmp,qvap*dens )
    gc( 1 ) = gc( 1 ) + dmp/dxmic
    qvap = qvap - dmp/dens
    qvap = max( qvap,0.0_RP )
    temp = temp + dmp/dens*qlevp/cp
  end if
  !
  return
  !
  end subroutine nucleat
  !-----------------------------------------------------------------------------
  subroutine nucleata       &
      ( dens, pres, dtime,  & !--- in
        gc, ga, qvap, temp  ) !--- inout
  !
  !  liquid nucleation from aerosol particle
  !
  !
  use scale_const, only: &
     cp    => CONST_CPdry, &
     rhow  => CONST_DWATR, &
     qlevp => CONST_LH0, &
     rvap  => CONST_Rvap
  use scale_atmos_saturation, only: &
       pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
       pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
  real(RP), intent(in) :: dens   !  density  [ kg/m3 ]
  real(RP), intent(in) :: pres   !  pressure [ Pa ]
  real(RP), intent(in) :: dtime
  !
  real(RP), intent(inout) :: gc( nbin )  !  SDF ( hydrometeors )
  real(RP), intent(inout) :: ga( nccn )  !  SDF ( aerosol ) : mass
  real(RP), intent(inout) :: qvap  !  specific humidity [ kg/kg ]
  real(RP), intent(inout) :: temp  !  temperature [ K ]
  !
  !--- local
  real(RP) :: gan( nccn )  !  SDF ( aerosol ) : number
  real(RP) :: ssliq, ssice, delcld
  real(RP) :: sumold, sumnew, acoef, bcoef, xcrit, rcrit
  real(RP) :: ractr, rcld, xcld, part, vdmp, dmp
  integer :: n, nc, ncrit
  integer, allocatable, save :: ncld( : )
  logical, save :: ofirst = .true.
  !
  !--- supersaturation
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  if ( ssliq <= 0.0_RP ) return
  !--- use for aerosol coupled model
  !--- mass -> number
  do n = 1, nccn
    gan( n ) = ga( n )/exp( xactr( n ) )
  end do

  acoef = 2.0_RP*sigma/rvap/rhow/temp
  bcoef = vhfct* rhoa/rhow * emwtr/emaer

  !--- relationship of bin number
  if ( ofirst ) then
    allocate ( ncld( 1:nccn ) )
    do n = 1, nccn
      ractr = ( exp( xactr( n ) )*ThirdovForth/pi/rhoa )**( OneovThird )
      rcld  = sqrt( 3.0_RP*bcoef*ractr*ractr*ractr / acoef )
      xcld  = log( rhow * 4.0_RP*pi*OneovThird*rcld*rcld*rcld )
     if( flg_nucl ) then
      ncld( n ) = 1
     else
      ncld( n ) = int( ( xcld-xctr( 1 ) )/dxmic ) + 1
      ncld( n ) = min( max( ncld( n ),1 ),nbin )
     end if
    end do
    ofirst = .false.
  end if

  !--- nucleation
  do n = nccn, 1, -1

      call pres2qsat_liq( ssliq,temp,pres )
      call pres2qsat_ice( ssice,temp,pres )
      ssliq = qvap/ssliq-1.0_RP
      ssice = qvap/ssice-1.0_RP

    if ( ssliq <= 0.0_RP ) exit
    !--- use for aerosol coupled model
    acoef = 2.0_RP*sigma/rvap/rhow/temp
    rcrit = acoef*OneovThird * ( 4.0_RP/bcoef )**( OneovThird ) / ssliq**( TwoovThird )
    xcrit = log( rhoa * 4.0_RP*pi*OneovThird * rcrit*rcrit*rcrit )
    ncrit = int( ( xcrit-xabnd( 1 ) )/dxaer ) + 1

    if ( n == ncrit ) then
      part = ( xabnd( ncrit+1 )-xcrit )/dxaer
    else if ( n > ncrit ) then
      part = 1.0_RP
    else
      exit
    end if

    nc = ncld( n )
    dmp = part*gan( n )*dxaer*exp( xctr( nc ) )
    dmp = min( dmp,qvap*dens )
    gc( nc ) = gc( nc ) + dmp/dxmic
    gan( n ) = gan( n ) - dmp/dxaer/exp( xctr( nc ) )
    gan( n ) = max( gan( n ), 0.0_RP )
    qvap = qvap - dmp/dens
    qvap = max( qvap,0.0_RP )
    temp = temp + dmp/dens*qlevp/cp
  end do

  !--- number -> mass
  do n = 1, nccn
    ga( n ) = gan( n )*exp( xactr( n ) )
  end do
  !
  return
  !
  end subroutine nucleata
  !-----------------------------------------------------------------------------
  subroutine  cndevpsbl     &
      ( dtime,              & !--- in
        dens, pres,         & !--- in
        gc, qvap, temp      ) !--- inout
  !
  real(RP), intent(in) :: dtime
  real(RP), intent(in) :: dens   !  atmospheric density [ kg/m3 ]
  real(RP), intent(in) :: pres   !  atmospheric pressure [ Pa ]
  real(RP), intent(inout) :: gc( nspc,nbin )  ! Size Distribution Function
  real(RP), intent(inout) :: qvap    !  specific humidity [ kg/kg ]
  real(RP), intent(inout) :: temp    !  temperature [ K ]
  !
  !--- local variables
  integer :: iflg( nspc ), n, m, iliq, iice
  real(RP) :: csum( nspc )
  real(RP) :: regene_gcn
  !
  !
  iflg( : ) = 0
  csum( : ) = 0.0_RP
  regene_gcn = 0.0_RP
  do m = 1, nspc
  do n = 1, nbin
    csum( m ) = csum( m )+gc( m,n )*dxmic
  end do
  end do

  do m = 1, nspc
   if ( csum( m ) > cldmin ) iflg( m ) = 1
  enddo

  iliq = iflg( il )
  iice = 0
  do m = 2, nspc
     iice = iice + iflg( m )
  enddo

  if ( iliq == 1 .and. iice == 0 ) then
      call  liqphase            &
              ( dtime, iliq,    & !--- in
                dens, pres,     & !--- in
                gc(il,1:nbin),  & !--- inout
                qvap, temp,     & !--- inout
                regene_gcn      ) !--- out
  elseif ( iliq == 0 .and. iice >= 1 ) then
      call icephase             &
              ( dtime, iflg,    & !--- in
                dens, pres,     & !--- in
                gc, qvap, temp  ) !--- inout
  elseif ( iliq == 1 .and. iice >= 1 ) then
      call mixphase             &
              ( dtime, iflg,    & !--- in
                dens, pres,     & !--- in
                gc, qvap, temp  ) !--- inout
  end if
  !
  end subroutine cndevpsbl
  !-----------------------------------------------------------------------------
  subroutine  cndevpsbla    &
      ( dtime,              & !--- in
        dens, pres,         & !--- in
        gc, ga, qvap, temp  ) !--- inout
  !
  real(RP), intent(in) :: dtime
  real(RP), intent(in) :: dens   !  atmospheric density [ kg/m3 ]
  real(RP), intent(in) :: pres   !  atmospheric pressure [ Pa ]
  real(RP), intent(inout) :: gc( nspc,nbin )  ! Size Distribution Function
  real(RP), intent(inout) :: ga( nccn )  !  SDF ( aerosol ) : mass
  real(RP), intent(inout) :: qvap    !  specific humidity [ kg/kg ]
  real(RP), intent(inout) :: temp    !  temperature [ K ]
  !
  !--- local variables
  integer :: iflg( nspc ), n, m, iliq, iice
  real(RP) :: csum( nspc )
  real(RP) :: regene_gcn
  !
  !
  iflg( : ) = 0
  csum( : ) = 0.0_RP
  regene_gcn = 0.0_RP
  do m = 1, nspc
  do n = 1, nbin
    csum( m ) = csum( m )+gc( m,n )*dxmic
  end do
  end do

  do m = 1, nspc
    if ( csum( m ) > cldmin ) iflg( m ) = 1
  enddo

  iliq = iflg( il )
  iice = 0
  do m = 2, nspc
    iice = iice + iflg( m )
  enddo 

  if ( iliq == 1 .and. iice == 0 ) then
      call  liqphase            &
              ( dtime, iliq,    & !--- in
                dens, pres,     & !--- in
                gc(il,1:nbin),  & !--- inout
                qvap, temp,     & !--- inout
                regene_gcn      ) !--- out
     !--- regeneration of aerosol
      if( flg_regeneration ) then
       call faero( regene_gcn,  & !--- in
                   ga           ) !--- inout
      end if
  elseif ( iliq == 0 .and. iice >= 1 ) then
      call icephase             &
              ( dtime, iflg,    & !--- in
                dens, pres,     & !--- in
                gc, qvap, temp  ) !--- inout
  elseif ( iliq == 1 .and. iice >= 1 ) then
      call mixphase             &
              ( dtime, iflg,    & !--- in
                dens, pres,     & !--- in
                gc, qvap, temp  ) !--- inout
  end if
  !
  end subroutine cndevpsbla
  !-----------------------------------------------------------------------------
  subroutine liqphase   &
      ( dtime, iflg,    & !--- in
        dens, pres,     & !--- in
        gc, qvap, temp, & !--- inout
        regene_gcn      ) !--- out
  !
  use scale_const, only: &
     cp    => CONST_CPdry, &
     rhow  => CONST_DWATR, &
     qlevp => CONST_LH0, &
     rvap  => CONST_Rvap
  use scale_atmos_saturation, only: &
     pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
     pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
  real(RP), intent(in) :: dtime
  integer, intent(in) :: iflg
  real(RP), intent(in) :: dens, pres
  real(RP), intent(inout) :: gc( nbin ), qvap, temp
  real(RP), intent(out) :: regene_gcn
  !
  !--- local variables
  integer :: n, nloop, ncount
  real(RP) :: gclold, gclnew, gtliq, umax, dtcnd
  real(RP) :: sumliq, cefliq, a, sliqtnd, cndmss, ssliq, ssice
  real(RP) :: gcn( nbin ), gdu( nbin+1 ), gcnold( nbin )
  real(RP), parameter :: cflfct = 0.50_RP
  real(RP) :: old_sum_gcn, new_sum_gcn
  !
  gclold = 0.0_RP
  do n = 1, nbin
    gclold = gclold + gc( n )
  end do
  gclold = gclold * dxmic
  !
  !------- mass -> number
  gcn( 1:nbin ) = gc( 1:nbin ) / exp( xctr( 1:nbin ) )

  !
  !------- CFL condition
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  gtliq = gliq( pres,temp )
  umax = cbnd( il,1 )/exp( xbnd( 1 ) )*gtliq*abs( ssliq )
  dtcnd = cflfct*dxmic/umax
  nloop = int( dtime/dtcnd ) + 1
  dtcnd = dtime / nloop
  !
  regene_gcn = 0.0_RP
  !------- loop
  do ncount = 1, nloop

  !----- matrix for supersaturation tendency
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  gtliq = gliq( pres,temp )
  sumliq = 0.0_RP
  old_sum_gcn = 0.0_RP
  do n = 1, nbin
    sumliq = sumliq + gcn( n )*cctr( il,n )
  end do
  sumliq = sumliq * dxmic

  if( qvap /= 0.0_RP ) then
   cefliq = ( ssliq+1.0_RP )*( 1.0_RP/qvap + qlevp*qlevp/cp/rvap/temp/temp )
   a = - cefliq*sumliq*gtliq/dens
   !
   !----- supersaturation tendency
   if ( abs( a*dtcnd ) >= 0.10_RP ) then
     sliqtnd = ssliq*( exp( a*dtcnd )-1.0_RP )/( a*dtcnd )
   else
     sliqtnd = ssliq
   end if
  elseif( qvap == 0.0_RP ) then
   sliqtnd = ssliq
  endif
  !
  !----- change of SDF
  gdu( 1:nbin+1 ) = cbnd( il,1:nbin+1 )/exp( xbnd( 1:nbin+1 ) )*gtliq*sliqtnd
  gcnold( : ) = gcn( : )
  call  advection              &
          ( dtcnd,             & !--- in
            gdu( 1:nbin+1 ),   & !--- in
            gcn( 1:nbin ),     & !--- inout
            regene_gcn         ) !--- inout
  !
  !----- new mass
  gclnew = 0.0_RP
  new_sum_gcn = 0.0_RP
  do n = 1, nbin
    gclnew = gclnew + gcn( n )*exp( xctr( n ) )
    old_sum_gcn = old_sum_gcn + gcnold( n )*dxmic
    new_sum_gcn = new_sum_gcn + gcn( n )*dxmic
  end do

  gclnew = gclnew*dxmic
  !
  !----- change of humidity and temperature
  cndmss = gclnew - gclold
  qvap = qvap - cndmss/dens
  temp = temp + cndmss/dens*qlevp/cp
  !
  gclold = gclnew
  !
  !----- continue/end
  end do
  !
  !------- number -> mass
  do n = 1 , nbin
   gc( n ) = gcn( n )*exp( xctr( n ) )
   if( gc( n ) < 0.0_RP ) then
     cndmss = -gc( n )
     gc( n ) = 0.0_RP
     qvap = qvap + cndmss/dens
     temp = temp - cndmss/dens*qlevp/cp
   endif
  enddo
  !
  end subroutine liqphase
  !-------------------------------------------------------------------------------
  subroutine icephase   &
      ( dtime, iflg,    & !--- in
        dens, pres,     & !--- in
        gc, qvap, temp  ) !--- inout
  !
  use scale_const, only : rhow   => CONST_DWATR,  &
                        qlevp  => CONST_LH0, &
                        qlsbl  => CONST_LHS0, &
                        rvap   => CONST_Rvap,  &
                        cp     => CONST_CPdry
  use scale_atmos_saturation, only: &
     pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
     pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
  !
  real(RP), intent(in) :: dtime
  integer, intent(in) :: iflg( nspc )
  real(RP), intent(in) :: dens, pres
  real(RP), intent(inout) :: gc( nspc,nbin ), qvap, temp
  !
  !--- local variables
  integer :: myu, n, nloop, ncount
  real(RP) :: gciold, gtice, umax, uval, dtcnd
  real(RP) :: sumice, cefice, d, sicetnd, gcinew, sblmss, ssliq, ssice
  real(RP) :: gcn( nspc,nbin ), gdu( nbin+1 )
  real(RP), parameter :: cflfct = 0.50_RP
  real(RP) :: dumm_regene
  !
  !
  !----- old mass
  gciold = 0.0_RP
  do n = 1, nbin
  do myu = 2, nspc
    gciold = gciold + gc( myu,n )
  end do
  end do
  gciold = gciold*dxmic

  !----- mass -> number
  do myu = 2, nspc
    gcn( myu,1:nbin ) = gc( myu,1:nbin )/exp( xctr( 1:nbin ) )
  end do

  !----- CFL condition
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  gtice = gice( pres,temp )
  umax = 0.0_RP
  do myu = 2, nspc
    uval = cbnd( myu,1 )/exp( xbnd( 1 ) )*gtice*abs( ssice )
    umax = max( umax,uval )
  end do

  dtcnd = cflfct*dxmic/umax
  nloop = int( dtime/dtcnd ) + 1
  dtcnd = dtime/nloop

!  ncount = 0
  !----- loop
!  1000 continue
  do ncount = 1, nloop
  !----- matrix for supersaturation tendency
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  gtice = gice( pres,temp )

  sumice = 0.0_RP
  do n = 1, nbin
    sumice = sumice + gcn( ic,n )*cctr( ic,n ) &
                    + gcn( ip,n )*cctr( ip,n ) &
                    + gcn( id,n )*cctr( id,n ) &
                    + gcn( iss,n )*cctr( iss,n ) &
                    + gcn( ig,n )*cctr( ig,n ) &
                    + gcn( ih,n )*cctr( ih,n )
  end do
  sumice = sumice*dxmic
  if( qvap /= 0.0_RP ) then
   cefice = ( ssice+1.0_RP )*( 1.0_RP/qvap + qlsbl*qlsbl/cp/rvap/temp/temp )
   d = - cefice*sumice*gtice/dens

   !----- supersaturation tendency
   if ( abs( d*dtcnd ) >= 0.10_RP ) then
     sicetnd = ssice*( exp( d*dtcnd )-1.0_RP )/( d*dtcnd )
   else
     sicetnd = ssice
   end if
  elseif( qvap == 0.0_RP ) then
   sicetnd = ssice
  endif
  !
  !----- change of SDF
  do myu = 2, nspc
    if ( iflg( myu ) == 1 ) then
      gdu( 1:nbin+1 ) = cbnd( myu,1:nbin+1 )/exp( xbnd( 1:nbin+1 ) )*gtice*sicetnd
      call  advection                                    &
              ( dtcnd, gdu( 1:nbin+1 ),                  & !--- in
                gcn( myu,1:nbin ),                       & !--- inout
                dumm_regene                              ) !--- out
    end if
  end do
  !
  !----- new mass
  gcinew = 0.0_RP
  do n = 1, nbin
  do myu = 2, nspc
    gcinew = gcinew + gcn( myu,n )*exp( xctr( n ) )
  end do
  end do
  gcinew = gcinew * dxmic
  !
  !----- change of humidity and temperature
  sblmss = gcinew - gciold
  qvap = qvap - sblmss/dens
  temp = temp + sblmss/dens*qlsbl/cp

  gciold = gcinew

  !
  !----- continue / end
  end do
!  ncount = ncount + 1
!  if ( ncount < nloop ) go to 1000
  !
  !------- number -> mass
  do myu = 2, nspc
    gc( myu,1:nbin ) = &
                  gcn( myu,1:nbin )*exp( xctr( 1:nbin ) )
  end do
  !
  end subroutine icephase
  !-------------------------------------------------------------------------------
  subroutine mixphase     &
      ( dtime, iflg,      & !--- in
        dens, pres,       & !--- in
        gc, qvap, temp    )!--- inout
  !
  use scale_const, only : rhow   => CONST_DWATR,  &
                        qlevp  => CONST_LH0, &
                        qlsbl  => CONST_LHS0, &
                        rvap   => CONST_Rvap,  &
                        cp     => CONST_CPdry
  use scale_atmos_saturation, only: &
     pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
     pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
  !
  real(RP), intent(in) :: dtime
  integer, intent(in) :: iflg( nspc )
  real(RP), intent(in) :: dens, pres
  !
  real(RP), intent(inout) :: gc( nspc,nbin )
  real(RP), intent(inout) :: qvap, temp
  !
  !--- local variables
  integer :: n, myu, nloop, ncount
  real(RP) :: gclold, gclnew, gciold, gcinew, cndmss, sblmss
  real(RP) :: gtliq, gtice, umax, uval, dtcnd
  real(RP) :: cef1, cef2, cef3, cef4, a, b, c, d
  real(RP) :: rmdplus, rmdmins, ssplus, ssmins, tplus, tmins
  real(RP) :: sliqtnd, sicetnd, ssliq, ssice, sumliq, sumice
  real(RP) :: gcn( nspc,nbin ), gdu( nbin+1 )
  real(RP), parameter :: cflfct = 0.50_RP
  real(RP) :: dumm_regene
  !
  !
  !----- old mass
  gclold = 0.0_RP
  do n = 1, nbin
    gclold = gclold + gc( il,n )
  end do
  gclold = gclold*dxmic

  gciold = 0.0_RP
  do n = 1, nbin
  do myu = 2, nspc
    gciold = gciold + gc( myu,n )
  end do
  end do
  gciold = gciold * dxmic

  !----- mass -> number
  do myu = 1, nspc
    gcn( myu,1:nbin ) = &
         gc( myu,1:nbin ) / exp( xctr( 1:nbin ) )
  end do

  !----- CFL condition
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  gtliq = gliq( pres,temp )
  gtice = gice( pres,temp )

  umax = cbnd( il,1 )/exp( xbnd( 1 ) )*gtliq*abs( ssliq )
  do myu = 2, nspc
    uval = cbnd( myu,1 )/exp( xbnd( 1 ) )*gtice*abs( ssice )
    umax = max( umax,uval )
  end do

  dtcnd = cflfct*dxmic/umax
  nloop = int( dtime/dtcnd ) + 1
  dtcnd = dtime/nloop

!  ncount = 0
  !----- loop
!  1000 continue
  do ncount = 1, nloop

  !-- matrix for supersaturation tendency
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  gtliq = gliq( pres,temp )
  gtice = gice( pres,temp )

  sumliq = 0.0_RP
  do n = 1, nbin
    sumliq = sumliq + gcn( il,n )*cctr( il,n )
  end do
  sumliq = sumliq*dxmic

  sumice = 0.0_RP
  do n = 1, nbin
    sumice = sumice + gcn( ic,n )*cctr( ic,n )  &
                    + gcn( ip,n )*cctr( ip,n )  &
                    + gcn( id,n )*cctr( id,n )  &
                    + gcn( iss,n )*cctr( iss,n )  &
                    + gcn( ig,n )*cctr( ig,n )  &
                    + gcn( ih,n )*cctr( ih,n )
  end do
  sumice = sumice*dxmic

  if( qvap /= 0.0_RP ) then
   cef1 = ( ssliq+1.0_RP )*( 1.0_RP/qvap + qlevp/rvap/temp/temp*qlevp/cp )
   cef2 = ( ssliq+1.0_RP )*( 1.0_RP/qvap + qlevp/rvap/temp/temp*qlsbl/cp )
   cef3 = ( ssice+1.0_RP )*( 1.0_RP/qvap + qlsbl/rvap/temp/temp*qlevp/cp )
   cef4 = ( ssice+1.0_RP )*( 1.0_RP/qvap + qlsbl/rvap/temp/temp*qlsbl/cp )

   a = - cef1*sumliq*gtliq/dens
   b = - cef2*sumice*gtice/dens
   c = - cef3*sumliq*gtliq/dens
   d = - cef4*sumice*gtice/dens

   !--- eigenvalues
   rmdplus = ( ( a+d ) + sqrt( ( a-d )**2 + 4.0_RP*b*c ) ) * 0.50_RP
   rmdmins = ( ( a+d ) - sqrt( ( a-d )**2 + 4.0_RP*b*c ) ) * 0.50_RP

   !--- supersaturation tendency
   ssplus = ( ( rmdmins-a )*ssliq - b*ssice )/b/( rmdmins-rmdplus )
   ssmins = ( ( a-rmdplus )*ssliq + b*ssice )/b/( rmdmins-rmdplus )

   if ( abs( rmdplus*dtcnd ) >= 0.10_RP ) then
     tplus = ( exp( rmdplus*dtcnd )-1.0_RP )/( rmdplus*dtcnd )
   else
     tplus = 1.0_RP
   end if

   if ( abs( rmdmins*dtcnd ) >= 0.10_RP ) then
     tmins = ( exp( rmdmins*dtcnd )-1.0_RP )/( rmdmins*dtcnd )
   else
     tmins = 1.0_RP
   end if
   sliqtnd = b*tplus*ssplus + b*tmins*ssmins
   sicetnd = ( rmdplus-a )*tplus*ssplus  &
           + ( rmdmins-a )*tmins*ssmins
  elseif( qvap == 0.0_RP ) then
   sliqtnd = ssliq
   sicetnd = ssice
  endif

  !--- change of SDF
  !- liquid
  if ( iflg( il ) == 1 ) then
    gdu( 1:nbin+1 ) = cbnd( il,1:nbin+1 )/exp( xbnd( 1:nbin+1 ) )*gtliq*sliqtnd
    call  advection                                  &
            ( dtcnd, gdu( 1:nbin+1 ),                & !--- in
              gcn( il,1:nbin ), & !--- inout
              dumm_regene                            ) !--- out
  end if

  !- ice
  do myu = 2, nspc
    if ( iflg( myu ) == 1 ) then
      gdu( 1:nbin+1 ) = cbnd( myu,1:nbin+1 )/exp( xbnd( 1:nbin+1 ) )*gtice*sicetnd
      call  advection                                    &
              ( dtcnd, gdu( 1:nbin+1 ),                  & !--- in
                gcn( myu,1:nbin ), & !--- inout
                dumm_regene                              ) !--- out
    end if
  end do

  !--- new mass
  gclnew = 0.0_RP
  do n = 1, nbin
    gclnew = gclnew + gcn( il,n )*exp( xctr( n ) )
  end do
  gclnew = gclnew*dxmic

  gcinew = 0.0_RP
  do n = 1, nbin
  do myu = 2, nspc
    gcinew = gcinew + gcn( myu,n )*exp( xctr( n ) )
  end do
  end do
  gcinew = gcinew*dxmic

  !--- change of humidity and temperature
  cndmss = gclnew - gclold
  sblmss = gcinew - gciold

  qvap = qvap - ( cndmss+sblmss )/dens
  temp = temp + ( cndmss*qlevp+sblmss*qlsbl )/dens/cp

  gclold = gclnew
  gciold = gcinew

  !--- continue/end
  end do
!  ncount = ncount + 1
!  if ( ncount < nloop ) go to 1000

  !----- number -> mass
  do myu = 1, nspc
    gc( myu,1:nbin ) = &
         gcn( myu,1:nbin )*exp( xctr( 1:nbin ) )
  end do

  end subroutine mixphase
  !-----------------------------------------------------------------------------
  subroutine advection  &
       ( dtime,         & !--- in
         gdu,           & !--- in
         gdq, regene    ) !--- inout
  !
  real(RP), intent(in) :: dtime, gdu( 1:nbin+1 )
  real(RP), intent(inout) :: gdq( 1:nbin ), regene
  !
  !--- local variables
  real(RP) :: delx
  real(RP) :: qadv( -1:nbin+2 ), uadv( 0:nbin+2 )
  real(RP) :: flq( 1:nbin+1 )
  integer, parameter :: ldeg = 2
  real(RP), parameter :: epsl = 1.E-10_RP
  real(RP) :: acoef( 0:nbin+1,0:ldeg ), sum
  real(RP) :: crn( 0:nbin+2 )
  real(RP) :: aip( 0:nbin+1 ), aim( 0:nbin+1 ), ai( 0:nbin+1 )
  real(RP) :: cmins, cplus
  integer :: i, n
  !
  !
  delx = dxmic
  do i = 1, nbin
    qadv( i ) = gdq( i )
  end do
  qadv( -1 )     = 0.0_RP
  qadv( 0  )     = 0.0_RP
  qadv( nbin+1 ) = 0.0_RP
  qadv( nbin+2 ) = 0.0_RP

  do i = 1, nbin+1
    uadv( i ) = gdu( i )
  end do
  uadv( 0  ) = 0.0_RP
  uadv( nbin+2 ) = 0.0_RP

  !--- flux
  do i = 0, nbin+1
    acoef( i,0 ) = - ( qadv( i+1 )-26.0_RP*qadv( i )+qadv( i-1 ) ) / 24.0_RP
    acoef( i,1 ) = ( qadv( i+1 )-qadv( i-1 ) ) * 0.50_RP
    acoef( i,2 ) = ( qadv( i+1 )-2.0_RP*qadv( i )+qadv( i-1 ) ) * 0.50_RP
  end do

  crn( 0:nbin+2 ) = uadv( 0:nbin+2 )*dtime/delx

  do i = 0, nbin+1
    cplus = ( crn( i+1 ) + abs( crn( i+1 ) ) ) * 0.50_RP
    sum = 0.0_RP
    do n = 0, ldeg
      sum = sum + acoef( i,n )/( n+1 )/2.0_RP**( n+1 )  &
                *( 1.0_RP-( 1.0_RP-2.0_RP*cplus )**( n+1 ) )
    end do
    aip( i ) = max( sum,0.0_RP )
  end do

  do i = 0, nbin+1
    cmins = - ( crn( i ) - abs( crn( i ) ) ) * 0.50_RP
    sum = 0.0_RP
    do n = 0, ldeg
      sum = sum + acoef( i,n )/( n+1 )/2.0_RP**( n+1 ) * (-1)**n &
                *( 1.0_RP-( 1.0_RP-2.0_RP*cmins )**( n+1 ) )
    end do
    aim( i ) = max( sum,0.0_RP )
  end do

  do i = 0, nbin+1
    sum = 0.0_RP
    do n = 0, ldeg
      sum = sum + acoef( i,n )/( n+1 )/2.0_RP**( n+1 ) * ( (-1)**n+1 )
    end do
    ai( i ) = max( sum,aip( i )+aim( i )+epsl )
  end do

  do i = 1, nbin+1
    flq( i ) = ( aip( i-1 )/ai( i-1 )*qadv( i-1 )  &
                -aim( i   )/ai( i   )*qadv( i   ) )*delx/dtime
  end do

  if( flg_regeneration .and. gdu( 1 ) < 0.0_RP ) then
   regene = regene+( -flq( 1 )*dtime/delx )
  else
   regene = 0.0_RP
  end if

  do i = 1, nbin
    gdq( i ) = gdq( i ) - ( flq( i+1 )-flq( i ) )*dtime/delx
  end do

  end subroutine advection
  !-----------------------------------------------------------------------------
  subroutine ice_nucleat        &
           ( dtime, dens, pres, & !--- in
             gc, qvap, temp     ) !--- inout
  !
  use scale_const, only : rhow  => CONST_DWATR, &
                        tmlt  => CONST_TMELT, &
                        qlmlt => CONST_EMELT, &
                        cp    => CONST_CPdry
  use scale_atmos_saturation, only: &
     pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq,   &
     pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice
  !
  real(RP), intent(in) :: dtime, dens, pres
  real(RP), intent(inout) :: gc( nspc,nbin ), temp, qvap
  !--- local variables
  real(RP) :: ssliq, ssice
  real(RP) :: numin, tdel, qdel
  real(RP), parameter :: n0 = 1.E+3_RP
  real(RP), parameter :: acoef = -0.639_RP, bcoef = 12.96_RP
  real(RP), parameter :: tcolmu = 269.0_RP, tcolml = 265.0_RP
  real(RP), parameter :: tdendu = 257.0_RP, tdendl = 255.0_RP
  real(RP), parameter :: tplatu = 250.6_RP
  !
  !--- supersaturation
  call pres2qsat_liq( ssliq,temp,pres )
  call pres2qsat_ice( ssice,temp,pres )
  ssliq = qvap/ssliq-1.0_RP
  ssice = qvap/ssice-1.0_RP

  if( ssice <= 0.0_RP ) return

  numin = bcoef * exp( acoef + bcoef * ssice )
  numin = numin * exp( xctr( 1 ) )/dxmic
  numin = min( numin,qvap*dens )
  !--- -4 [deg] > T >= -8 [deg] and T < -22.4 [deg] -> column
  if( temp <= tplatu .or. ( temp >= tcolml .and. temp < tcolmu ) ) then
   gc( ic,1 ) = gc( ic,1 ) + numin
  !--- -14 [deg] > T >= -18 [deg] -> dendrite
  elseif( temp <= tdendu .and. temp >= tdendl ) then
   gc( id,1 ) = gc( id,1 ) + numin
  !--- else -> plate
  else
   gc( ip,1 ) = gc( ip,1 ) + numin
  endif

  tdel = numin/dens*qlmlt/cp
  temp = temp + tdel
  qdel = numin/dens
  qvap = qvap - qdel

  return
  !
  end subroutine ice_nucleat
  !-----------------------------------------------------------------------------
  subroutine  freezing       &
      ( dtime, dens,  & !--- in
        gc, temp  )            !--- inout
  !
  use scale_const, only : rhow  => CONST_DWATR,  &
                        tmlt  => CONST_TMELT,  &
                        qlmlt => CONST_EMELT, &
                        cp    => CONST_CPdry
  !
  real(RP), intent(in) :: dtime, dens
  real(RP), intent(inout) :: gc( nspc,nbin ), temp
  !
  !--- local variables
  integer :: nbound, n
  real(RP) :: xbound, tc, rate, dmp, frz, sumfrz, tdel
  real(RP), parameter :: coefa = 1.0E-01_RP, coefb = 0.66_RP
  real(RP), parameter :: rbound = 2.0E-04_RP, tthreth = 235.0_RP
  real(RP), parameter :: ncoefim = 1.0E+7_RP, gamm = 3.3_RP
  !
  !
!  if( temp <= tthreth ) then !--- Bigg (1975)
   xbound = log( rhow * 4.0_RP*pi/3.0_RP * rbound**3 )
   nbound = int( ( xbound-xbnd( 1 ) )/dxmic ) + 1

   tc = temp-tmlt
   rate = coefa*exp( -coefb*tc )

   sumfrz = 0.0_RP
   do n = 1, nbin
     dmp = rate*exp( xctr( n ) )
     frz = gc( il,n )*( 1.0_RP-exp( -dmp*dtime ) )

     gc( il,n ) = gc( il,n ) - frz
     gc( il,n ) = max( gc( il,n ),0.0_RP )

     if ( n >= nbound ) then
       gc( ih,n ) = gc( ih,n ) + frz
     else
       gc( ip,n ) = gc( ip,n ) + frz
     end if

     sumfrz = sumfrz + frz
   end do
!  elseif( temp > tthreth ) then !--- Vali (1975)
!
!   tc = temp-tmlt
!   dmp = ncoefim * ( 0.1_RP * ( tmlt-temp )**gamm )
!   sumfrz = 0.0_RP
!   do n = 1, nbin
!    frz = gc( (il-1)*nbin+n )*dmp*xctr( n )/dens
!    frz = max( frz,gc( (il-1)*nbin+n )*dmp*xctr( n )/dens )
!    gc( (il-1)*nbin+n ) = gc( (il-1)*nbin+n ) - frz
!    if( n >= nbound ) then
!     gc( (ih-1)*nbin+n ) = gc( (ih-1)*nbin+n ) + frz
!    else
!     gc( (ip-1)*nbin+n ) = gc( (ip-1)*nbin+n ) + frz
!    end if
!
!    sumfrz = sumfrz + frz
!   end do
!  endif
  sumfrz = sumfrz*dxmic

  tdel = sumfrz/dens*qlmlt/cp
  temp = temp + tdel
  !
  return
  !
  end subroutine freezing
  !-----------------------------------------------------------------------------
  subroutine  melting  &
      ( dtime, dens,   & !--- in
        gc, temp       ) !--- inout
  !
  use scale_const, only : qlmlt => CONST_EMELT, &
                        cp    => CONST_CPdry
  !
  real(RP), intent(in) :: dtime, dens
  !
  real(RP), intent(inout) :: gc( nspc,nbin ), temp
  !
  !--- local variables
  integer :: n, m
  real(RP) :: summlt, sumice, tdel
  !
  !
  summlt = 0.0_RP
  do n = 1, nbin
   sumice = 0.0_RP
   do m = ic, ih
    sumice = sumice + gc( m,n )
    gc( m,n ) = 0.0_RP
   end do
   gc( il,n ) = gc( il,n ) + sumice
   summlt = summlt + sumice
  end do
  summlt = summlt*dxmic

  tdel = - summlt/dens*qlmlt/cp
  temp = temp + tdel
  !
  return
  !
  end subroutine melting
  !-------------------------------------------------------------------------------
  function  gliq( pres,temp )
  !
  use scale_const, only : rair  => CONST_Rdry,  &
                        rvap  => CONST_Rvap,  &
                        qlevp => CONST_LH0
  !
  real(RP), intent(in) :: pres, temp
  real(RP) :: gliq
  !
  real(RP) :: emu, dens, cefd, cefk, f1, f2
  real(RP), parameter :: fct = 1.4E+3_RP
  !
  emu = fmyu( temp )
  dens = pres/rair/temp
  cefd = emu/dens
  cefk = fct*emu

  f1 = rvap*temp/fesatl( temp )/cefd
  f2 = qlevp/cefk/temp*( qlevp/rvap/temp - 1.0_RP )

  gliq = 4.0_RP*pi/( f1+f2 )

  end function gliq
  !-------------------------------------------------------------------------------
  function gice ( pres, temp )
  !
  use scale_const, only : rair  => CONST_Rdry,  &
                        rvap  => CONST_Rvap,  &
                        qlsbl => CONST_LHS0
  !
  real(RP), intent(in) :: pres, temp
  real(RP) :: gice
  !
  real(RP) :: emu, dens, cefd, cefk, f1, f2
  real(RP), parameter :: fct = 1.4E+3_RP
  !
  emu = fmyu( temp )
  dens = pres/rair/temp
  cefd = emu/dens
  cefk = fct*emu

  f1 = rvap*temp/fesati( temp )/cefd
  f2 = qlsbl/cefk/temp*( qlsbl/rvap/temp - 1.0_RP )

  gice = 4.0_RP*pi/( f1+f2 )

  end function gice
  !-------------------------------------------------------------------------------
  function fmyu( temp )
  !
  !
  use scale_const, only: &
     tmlt => CONST_TMELT
  real(RP), intent(in) :: temp
  real(RP) :: fmyu
  !
  real(RP), parameter :: a = 1.72E-05_RP, b = 3.93E+2_RP, c = 1.2E+02_RP
  !
  fmyu = a*( b/( temp+c ) )*( temp/tmlt )**1.50_RP
  !
  end function fmyu
  !------------------------------------------------------------------------
  function fesati( temp )
  !
  use scale_const, only : temp0 => CONST_TEM00, &
                        rvap  => CONST_RVAP,  &
                        qlsbl => CONST_LHS0,  &
                        esat0 => CONST_PSAT0
  !
  real(RP), intent(in) :: temp
  real(RP) :: fesati
  !
  fesati = esat0*exp( qlsbl/rvap*( 1.0_RP/temp0 - 1.0_RP/temp ) )
  !
  return

  end function fesati
  !-------------------------------------------------------------------------------
  function fesatl( temp )
  !
  use scale_const, only: &
     temp0 => CONST_TEM00, &
     esat0 => CONST_PSAT0, &
     qlevp => CONST_LH0, &
     rvap  => CONST_Rvap
  real(RP), intent(in) :: temp
  real(RP) :: fesatl
  !
  fesatl = esat0*exp( qlevp/rvap*( 1.0_RP/temp0 - 1.0_RP/temp ) )
  !
  return

  end function fesatl
  !-----------------------------------------------------------------------
  subroutine collmain &
       ( dtime,       & !--- in
         gc           ) !--- inout
  !
  real(RP), intent(in) :: dtime
  real(RP), intent(inout) :: gc( nspc,nbin )
  !
  !--- local variables
  integer :: iflg( il ), n
  real(RP) :: csum( il )
  !
  !--- judgement of particle existence
    iflg( : ) = 0
    csum( : ) = 0.0_RP
    do n = 1, nbin
      csum( il ) = csum( il ) + gc( il,n )*dxmic
    end do
    if ( csum( il ) > cldmin ) iflg( il ) = 1
  !
  !--- interaction
   if ( iflg( il ) == 1 ) then
      if ( rndm_flgp == 1 ) then  !--- stochastic method

        call r_collcoag        &
            ( il, il, il,      &
              dtime,  wgtbin,  & !--- in
              gc               ) !--- inout

      else  !--- default method

        call  collcoag    &
            ( il, il, il, &
              dtime,      & !--- in
              gc          ) !--- inout

      end if
   end if
  !
  return
  !
  end subroutine collmain
  !-------------------------------------------------------------------------------
  subroutine collmainf &
       ( temp, dtime,  & !--- in
         gc            ) !--- inout

  use scale_process, only: &
     PRC_MPIstop
  !
  real(RP), intent(in) :: dtime, temp
  real(RP), intent(inout) :: gc( nspc,nbin )
  !
  !--- local variables
  integer :: iflg( nspc ), n, irsl, ilrg, isml
  real(RP) :: csum( nspc )
  !
  !--- judgement of particle existence
    iflg( : ) = 0
    csum( : ) = 0.0_RP
    do n = 1, nbin
      csum( il ) = csum( il ) + gc( il,n )*dxmic
      csum( ic ) = csum( ic ) + gc( ic,n )*dxmic
      csum( ip ) = csum( ip ) + gc( ip,n )*dxmic
      csum( id ) = csum( id ) + gc( id,n )*dxmic
      csum( iss ) = csum( iss ) + gc( iss,n )*dxmic
      csum( ig ) = csum( ig ) + gc( ig,n )*dxmic
      csum( ih ) = csum( ih ) + gc( ih,n )*dxmic
    end do
    if ( csum( il ) > cldmin ) iflg( il ) = 1
    if ( csum( ic ) > cldmin ) iflg( ic ) = 1
    if ( csum( ip ) > cldmin ) iflg( ip ) = 1
    if ( csum( id ) > cldmin ) iflg( id ) = 1
    if ( csum( iss ) > cldmin ) iflg( iss ) = 1
    if ( csum( ig ) > cldmin ) iflg( ig ) = 1
    if ( csum( ih ) > cldmin ) iflg( ih ) = 1
  !
  !--- rule of interaction
    call  getrule    &
            ( temp,  & !--- in
              ifrsl  ) !--- out
  !
  !--- interaction
  do isml = 1, nspc
   if ( iflg( isml ) == 1 ) then
    do ilrg = 1, nspc
     if ( iflg( ilrg ) == 1 ) then
      irsl = ifrsl( isml, ilrg )
      if ( rndm_flgp == 1 ) then  !--- stochastic method
        call r_collcoag         &
            ( isml, ilrg, irsl, & !--- in
              dtime,  wgtbin,   & !--- in
              gc                ) !--- inout
      else  !--- default method
        call  collcoag          &
            ( isml, ilrg, irsl, & !--- in
              dtime,            & !--- in
              gc                ) !--- inout

      end if
     end if
    end do
   end if
  end do
  !
  return
  !
  end subroutine collmainf
  !-----------------------------------------------------------------------
  subroutine  collcoag( isml,ilrg,irsl,dtime,gc )
  !------------------------------------------------------------------------------
  !--- reference paper
  !    Bott et al. (1998) J. Atmos. Sci. vol.55, pp. 2284-
  !------------------------------------------------------------------------------
  !
  integer,  intent(in) :: isml, ilrg, irsl
  real(RP), intent(in) :: dtime
  real(RP), intent(inout) :: gc( nspc,nbin )
  !
  !--- local variables
  integer :: i, j, k, l
  real(RP) :: xi, xj, xnew, dmpi, dmpj, frci, frcj
  real(RP) :: gprime, gprimk, wgt, crn, sum, flux
  integer, parameter :: ldeg = 2
  real(RP) :: acoef( 0:ldeg )
  real(RP), parameter :: dmpmin = 1.E-01_RP
  real(RP) :: suri, surj
  !
  small : do i = 1, nbin-1
    if ( gc( isml,i ) <= cldmin ) cycle small
  large : do j = i+1, nbin
    if ( gc( ilrg,j ) <= cldmin ) cycle large

    xi = exp( xctr( i ) )
    xj = exp( xctr( j ) )
    xnew = log( xi+xj )
    k = int( ( xnew-xctr( 1 ) )/dxmic ) + 1
    k = max( max( k,j ),i )
    if ( k >= nbin ) cycle small

    dmpi = ck( isml,ilrg,i,j )*gc( ilrg,j )/xj*dxmic*dtime
    dmpj = ck( ilrg,isml,i,j )*gc( isml,i )/xi*dxmic*dtime

    if ( dmpi <= dmpmin ) then
      frci = gc( isml,i )*dmpi
    else
      frci = gc( isml,i )*( 1.0_RP-exp( -dmpi ) )
    end if

    if ( dmpj <= dmpmin ) then
      frcj = gc( ilrg,j )*dmpj
    else
      frcj = gc( ilrg,j )*( 1.0_RP-exp( -dmpj ) )
    end if

    gprime = frci+frcj
    if ( gprime <= 0.0_RP ) cycle large

    suri = gc( isml,i )
    surj = gc( ilrg,j )
    gc( isml,i ) = gc( isml,i )-frci
    gc( ilrg,j ) = gc( ilrg,j )-frcj
    gc( isml,i ) = max( gc( isml,i )-frci, 0.0_RP )
    gc( ilrg,j ) = max( gc( ilrg,j )-frcj, 0.0_RP )
    frci = suri - gc( isml,i )
    frcj = surj - gc( ilrg,j )
    gprime = frci+frcj

    gprimk = gc( irsl,k ) + gprime
    wgt = gprime / gprimk
    crn = ( xnew-xctr( k ) )/( xctr( k+1 )-xctr( k ) )

    acoef( 0 ) = -( gc( irsl,k+1 )-26.0_RP*gprimk+gc( irsl,k-1 ) )/24.0_RP
    acoef( 1 ) =  ( gc( irsl,k+1 )-gc( irsl,k-1 ) ) *0.50_RP
    acoef( 2 ) =  ( gc( irsl,k+1 )-2.0_RP*gprimk+gc( irsl,k-1 ) ) *0.50_RP

    sum = 0.0_RP
    do l = 0, ldeg
      sum = sum + acoef( l )/( l+1 )/2.0_RP**( l+1 )   &
                *( 1.0_RP-( 1.0_RP-2.0_RP*crn )**( l+1 ) )
    end do

    flux = wgt*sum
    flux = min( max( flux,0.0_RP ),gprime )

    gc( irsl,k ) = gprimk - flux
    gc( irsl,k+1 ) = gc( irsl,k+1 ) + flux

  end do large
  end do small
  !
  return
  !
  end subroutine collcoag
  !-------------------------------------------------------------------
  subroutine  getrule   &
      ( temp,           & !--- in
        ifrsl           ) !--- out
  !
  real(RP), intent(in) :: temp
  integer, intent(out) :: ifrsl( 7,7 )
  !
  real(RP), parameter :: tcrit = 271.15_RP
  !
  !--- liquid + liquid -> liquid
  ifrsl( il,il ) = il
  !
  !--- liquid + column -> ( graupel, hail ) + column
  ifrsl( il,ic ) = ic
  if ( temp < tcrit ) then
    ifrsl( ic,il ) = ig
  else
    ifrsl( ic,il ) = ih
  end if
  !
  !--- liquid + plate -> ( graupel, hail ) + plate
  ifrsl( il,ip ) = ip
  if ( temp < tcrit ) then
    ifrsl( ip,il ) = ig
  else
    ifrsl( ip,il ) = ih
  end if
  !
  !--- liquid + dendrite -> ( graupel, hail ) + dendrite
  ifrsl( il,id ) = id
  if ( temp < tcrit ) then
    ifrsl( id,il ) = ig
  else
    ifrsl( id,il ) = ih
  end if
  !
  !--- liquid + snowflake -> ( graupel, hail ) + snowflake
  ifrsl( il,iss ) = iss
  if ( temp < tcrit ) then
    ifrsl( iss,il ) = ig
  else
    ifrsl( iss,il ) = ih
  end if
  !
  !--- liquid + graupel -> ( graupel, hail )
  ifrsl( il,ig ) = ig
  if ( temp < tcrit ) then
    ifrsl( ig,il ) = ig
  else
    ifrsl( ig,il ) = ih
  end if
  !
  !--- liquid + hail -> ( graupel, hail )
  ifrsl( il,ih ) = ih
  if ( temp < tcrit ) then
    ifrsl( ih,il ) = ig
  else
    ifrsl( ih,il ) = ih
  end if
  !
  !
  !--- column + column -> snowflake
  ifrsl( ic,ic ) = iss
  !
  !--- column + plate -> snowflake
  ifrsl( ic,ip ) = iss
  ifrsl( ip,ic ) = iss
  !
  !--- column + dendrite -> snowflake
  ifrsl( ic,id ) = iss
  ifrsl( id,ic ) = iss
  !
  !--- column + snowflake -> snowflake
  ifrsl( ic,iss ) = iss
  ifrsl( iss,ic ) = iss
  !
  !--- column + graupel -> column + graupel
  ifrsl( ic,ig ) = ig
  ifrsl( ig,ic ) = ic
  !
  !--- column + hail -> column + ( graupel, hail )
  ifrsl( ih,ic ) = ic
  if ( temp < tcrit ) then
    ifrsl( ic,ih ) = ig
  else
    ifrsl( ic,ih ) = ih
  end if
  !
  !
  !--- plate + plate -> snowflake
  ifrsl( ip,ip ) = iss
  !
  !--- plate + dendrite -> snowflake
  ifrsl( ip,id ) = iss
  ifrsl( id,ip ) = iss
  !
  !--- plate + snowflake -> snowflake
  ifrsl( ip,iss ) = iss
  ifrsl( iss,ip ) = iss
  !
  !--- plate + graupel -> plate + graupel
  ifrsl( ip,ig ) = ig
  ifrsl( ig,ip ) = ip
  !
  !--- plate + hail -> plate + ( graupel, hail )
  ifrsl( ih,ip ) = ip
  if ( temp < tcrit ) then
    ifrsl( ip,ih ) = ig
  else
    ifrsl( ip,ih ) = ih
  end if
  !
  !
  !--- dendrite + dendrite -> snowflake
  ifrsl( id,id ) = iss
  !
  !--- dendrite + snowflake -> snowflake
  ifrsl( id,iss ) = iss
  ifrsl( iss,id ) = iss
  !
  !--- dendrite + graupel -> dendrite + graupel
  ifrsl( id,ig ) = ig
  ifrsl( ig,id ) = id
  !
  !--- dendrite + hail -> dendrite + ( graupel, hail )
  ifrsl( ih,id ) = id
  if ( temp < tcrit ) then
    ifrsl( id,ih ) = ig
  else
    ifrsl( id,ih ) = ih
  end if
  !
  !
  !--- snowflake + snowflake -> snowflake
  ifrsl( iss,iss ) = iss
  !
  !--- snowflake + graupel -> snowflake + graupel
  ifrsl( iss,ig ) = ig
  ifrsl( ig,iss ) = iss
  !
  !--- snowflake + hail -> snowflake + ( graupel, hail )
  ifrsl( ih,iss ) = iss
  if ( temp < tcrit ) then
    ifrsl( iss,ih ) = ig
  else
    ifrsl( iss,ih ) = ih
  end if
  !
  !
  !--- graupel + graupel -> graupel
  ifrsl( ig,ig ) = ig
  !
  !--- graupel + hail -> ( graupel, hail )
  if ( temp < tcrit ) then
    ifrsl( ig,ih ) = ig
    ifrsl( ih,ig ) = ig
  else
    ifrsl( ig,ih ) = ih
    ifrsl( ih,ig ) = ih
  end if
  !
  !--- hail + hail -> hail
  ifrsl( ih,ih ) = ih
  !
  !
  return
  !
  end subroutine getrule
  !-----------------------------------------------------------------------------
  subroutine faero( f0,ga )

  real(RP), intent(in) ::  f0
  real(RP), intent(inout) :: ga( nccn )
  real(RP) :: gaero( nccn ) !, f1, radmax, radmin
  real(RP), parameter :: alpha = 3.0_RP
  integer :: n

  do n = 1, nccn
   gaero( n ) = f0*marate( n )*exp( xactr( n ) )/dxaer
   ga( n ) = ga( n )+gaero( n )
  end do

  return

  end subroutine faero
  !-------------------------------------------------------------------------------
  !
  ! + Y. Sato added for stochastic method
  ! + Reference Sato et al. (2009) JGR, doi:10.1029/2008JD011247
  !
  !-------------------------------------------------------------------------------
  subroutine random_setup( mset ) !--- in

   use scale_random, only: &
       RANDOM_get
   use scale_process, only: &
       PRC_MPIstop

   integer, intent(in) :: mset

   !--- local ----
   integer :: n
   real(RP) :: nbinr, tmp1
   real(RP) :: rans( mbin ), ranl( mbin )
   integer  :: pq
   real(RP), allocatable :: ranstmp( : )
   real(RP), allocatable :: ranltmp( : )
   integer :: p, q
   integer :: k, temp
   integer, allocatable :: orderb( : )
   real(RP) :: abq1
   real(RP) :: a
   real(RP), allocatable :: randnum(:,:,:)
  !-------------------------------------------------------
   pq = nbin*(nbin-1)/2
   allocate( blrg( mset, mbin ) )
   allocate( bsml( mset, mbin ) )
   allocate( ranstmp( pq ) )
   allocate( ranltmp( pq ) )
   allocate( orderb( pq ) )
   allocate( randnum(1,1,pq) )

   a = real( nbin )*real( nbin-1 )*0.50_RP
   if( a < mbin ) then
    write(*,*) "mbin should be smaller than {nbin}_C_{2}"
    call PRC_MPIstop
   end if

   wgtbin = a/real( mbin )
   nbinr = real( nbin )

    do p = 1, pq
      orderb( p ) = p
    end do

    do p = 1, nbin-1
      ranstmp( (p-1)*nbin-(p*(p-1))/2+1 : p*nbin-(p*(p+1))/2 ) = p
     do q = 1, nbin-p
        ranltmp( (p-1)*nbin-(p*(p-1))/2+q ) = p+q
      end do
   end do

    do n = 1, mset
      call RANDOM_get( randnum )
       do p = 1, pq
        abq1 = randnum( 1,1,p )
        k = int( abq1*( pq-p-1 ) ) + p
        temp = orderb( p )
        orderb( p ) = orderb( k )
        orderb( k ) = temp
       end do

       do p = 1, mbin
        if( p <= pq ) then
         rans( p ) = ranstmp( orderb( p ) )
         ranl( p ) = ranltmp( orderb( p ) )
        else
         rans( p ) = ranstmp( orderb( p-pq ) )
         ranl( p ) = ranltmp( orderb( p-pq ) )
        end if
         if( rans( p ) >= ranl( p ) ) then
          tmp1 = rans( p )
          rans( p ) = ranl( p )
          ranl( p ) = tmp1
         end if
       end do
         blrg( n,1:mbin ) = int( ranl( 1:mbin ) )
         bsml( n,1:mbin ) = int( rans( 1:mbin ) )
    end do

    deallocate( ranstmp )
    deallocate( ranltmp )
    deallocate( orderb )
    deallocate( randnum )

  end subroutine random_setup
  !-----------------------------------------------------------------------------
  subroutine  r_collcoag( isml, ilrg, irsl, dtime, swgt, gc )
  !-------------------------------------------------------------------------------
  !--- reference paper
  !    Bott et al. (1998) J. Atmos. Sci. vol.55, pp. 2284-
  !    Bott et al. (2000) J. Atmos. Sci. Vol.57, pp. 284-
  !-------------------------------------------------------------------------------
  !
  use scale_random, only: &
      RANDOM_get

  integer,  intent(in) :: isml, ilrg, irsl
  real(RP), intent(in) :: dtime
  real(RP), intent(in) :: swgt
  real(RP), intent(inout) :: gc( nspc,nbin )
  !
  !--- local variables
  integer :: i, j, k, l
  real(RP) :: xi, xj, xnew, dmpi, dmpj, frci, frcj
  real(RP) :: gprime, gprimk, wgt, crn, sum, flux
  integer, parameter :: ldeg = 2
  real(RP), parameter :: dmpmin = 1.E-01_RP, cmin = 1.E-10_RP
  real(RP) :: acoef( 0:ldeg )
  !
  !--- Y.sato added to use code6
  integer :: nums( mbin ), numl( mbin )
  real(RP), parameter :: gt = 1.0_RP
  integer :: s, det
  real(RP) :: nbinr, mbinr        ! use to weight
!  real(RP) :: beta
  real(RP) :: tmpi, tmpj
 !-----------------------------------------------------
  call RANDOM_get( rndm )
  det = int( rndm(1,1,1)*IA*JA*KA )
  nbinr = real( nbin )
  mbinr = real( mbin )
  nums( 1:mbin ) = bsml( det,1:mbin )
  numl( 1:mbin ) = blrg( det,1:mbin )

   do s = 1, mbin
    i = nums( s )
    j = numl( s )

    if ( gc( isml,i ) <= cmin ) cycle !small
    if ( gc( ilrg,j ) <= cmin ) cycle !large

    xi = exp( xctr( i ) )
    xj = exp( xctr( j ) )
    xnew = log( xi+xj )
    k = int( ( xnew-xctr( 1 ) )/dxmic ) + 1
    k = max( max( k,j ),i )
    if( k>= nbin ) cycle

    dmpi = ck( isml,ilrg,i,j )*gc( ilrg,j )/xj*dxmic*dtime
    dmpj = ck( ilrg,isml,i,j )*gc( isml,i )/xi*dxmic*dtime

    if ( dmpi <= dmpmin ) then
      frci = gc( isml,i )*dmpi
    else
      frci = gc( isml,i )*( 1.0_RP-exp( -dmpi ) )
    end if

    if ( dmpj <= dmpmin ) then
      frcj = gc( ilrg,j )*dmpj
    else
      frcj = gc( ilrg,j )*( 1.0_RP-exp( -dmpj ) )
    end if
    tmpi = gc( isml,i )
    tmpj = gc( ilrg,j )

    gc( isml,i ) = gc( isml,i )-frci*swgt
    gc( ilrg,j ) = gc( ilrg,j )-frcj*swgt

    if( j /= k ) then
     gc( ilrg,j ) = max( gc( ilrg,j ), 0.0_RP )
    end if
     gc( isml,i ) = max( gc( isml,i ), 0.0_RP )

    frci = tmpi - gc( isml,i )
    frcj = tmpj - gc( ilrg,j )

    gprime = frci+frcj

    !-----------------------------------------------
    !--- Exponential Flux Method (Bott, 2000, JAS)
    !-----------------------------------------------
!    if ( gprime <= 0.0_RP ) cycle !large
!    gprimk = gc( (irsl-1)*nbin+k ) + gprime
!
!    beta = log( gc( (irsl-1)*nbin+k+1 )/gprimk+1.E-60_RP )
!    crn = ( xnew-xctr( k ) )/( xctr( k+1 )-xctr( k ) )
!
!    flux = ( gprime/beta )*( exp( beta*0.50_RP ) -exp( beta*( 0.50_RP-crn ) ) )
!    flux = min( gprimk ,gprime )
!
!    gc( (irsl-1)*nbin+k ) = gprimk - flux
!    gc( (irsl-1)*nbin+k+1 ) = gc( (irsl-1)*nbin+k+1 ) + flux

    !-----------------------------------------------
    !--- Flux Method (Bott, 1998, JAS)
    !-----------------------------------------------
    if ( gprime <= 0.0_RP ) cycle !large
    gprimk = gc( irsl,k ) + gprime
    wgt = gprime / gprimk
    crn = ( xnew-xctr( k ) )/( xctr( k+1 )-xctr( k ) )

    acoef( 0 ) = -( gc( irsl,k+1 )-26.0_RP*gprimk+gc( irsl,k-1 ) )/24.0_RP
    acoef( 1 ) =  ( gc( irsl,k+1 )-gc( irsl,k-1 ) ) *0.5_RP
    acoef( 2 ) =  ( gc( irsl,k+1 )-2.0_RP*gprimk+gc( irsl,k-1 ) ) *0.50_RP

    sum = 0.0_RP
    do l = 0, ldeg
      sum = sum + acoef( l )/( l+1 )/2.0_RP**( l+1 )   &
                *( 1.0_RP-( 1.0_RP-2.0_RP*crn )**( l+1 ) )
    end do

    flux = wgt*sum
    flux = min( max( flux,0.0_RP ),gprime )

    gc( irsl,k ) = gprimk - flux
    gc( irsl,k+1 ) = gc( irsl,k+1 ) + flux

   end do

  !
  return
  !
  end subroutine r_collcoag
  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_suzuki10_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QAD)

    real(RP) :: qhydro
    integer  :: k, i, j, iq, ihydro
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do ihydro = 1, MP_QA
        do iq = I_QV+nbin*(ihydro-1)+1, I_QV+nbin*ihydro
          qhydro = qhydro + QTRC(k,i,j,iq)
        enddo
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_suzuki10_CloudFraction
  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_suzuki10_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0  )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)        ! density                   [kg/m3]

    real(RP) :: sum2(KA,IA,JA,MP_QA), sum3(KA,IA,JA,MP_QA)
    integer  :: i, j, k, iq, ihydro
    !---------------------------------------------------------------------------

    sum2(:,:,:,:) = 0.0_RP
    sum3(:,:,:,:) = 0.0_RP

    do ihydro = 1, MP_QA
    do k = KS, KE
    do j = JS, JE
    do i = IS, JE
      do iq = I_QV+nbin*(ihydro-1)+1, I_QV+nbin*ihydro
         sum3(k,i,j,ihydro) = sum3(k,i,j,ihydro) + &
                            ( ( QTRC0(k,i,j,iq) * DENS0(k,i,j) ) & !--- [kg/kg] -> [kg/m3]
                            / exp( xctr( iq-(I_QV+nbin*(ihydro-1)+iq) ) ) &   !--- mass -> number
                            * radc( iq-(I_QV+nbin*(ihydro-1)+iq) )**3.0_RP )
         sum2(k,i,j,ihydro) = sum2(k,i,j,ihydro) + &
                            ( ( QTRC0(k,i,j,iq) * DENS0(k,i,j) ) & !--- [kg/kg] -> [kg/m3]
                            / exp( xctr( iq-(I_QV+nbin*(ihydro-1)+iq) ) ) &   !--- mass -> number
                            * radc( iq-(I_QV+nbin*(ihydro-1)+iq) )**2.0_RP )
      enddo
      sum2(k,i,j,ihydro) = 0.5_RP + sign(0.5_RP,sum2(k,i,j,ihydro)-EPS)
      sum3(k,i,j,ihydro) = 0.5_RP + sign(0.5_RP,sum3(k,i,j,ihydro)-EPS)

      if( sum2(k,i,j,ihydro) /= 0.0_RP ) then
       Re(k,i,j,ihydro) = sum3(k,i,j,ihydro) / sum2(k,i,j,ihydro)
      else
       Re(k,i,j,ihydro) = 0.0_RP
      endif
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_suzuki10_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_suzuki10_Mixingratio( &
       Qe,    &
       QTRC0  )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]

    real(RP) :: sum2
    integer  :: i, j, k, iq, ihydro
    !---------------------------------------------------------------------------

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
      do ihydro = 1, MP_QA
        sum2 = 0.0_RP
        do iq = I_QV+nbin*(ihydro-1)+1, I_QV+nbin*ihydro
          sum2 = sum2 + QTRC0(k,i,j,iq)
        enddo
        Qe(k,i,j,ihydro) = sum2
      enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_suzuki10_Mixingratio
  !-----------------------------------------------------------------------------
  !----- mkpara is module to create micpara.dat, which is parameter file of
  !----- micrphyisical proprties of hydrometeors (collision kernel, radius...).
  !----- Imported from preprocess/mk_para2 at 2013/12/26 (Y.Sato)
  !-----------------------------------------------------------------------------
  subroutine mkpara

  implicit none

  integer :: i, j
  !--------------------------------------------------------------------------------------

  allocate( radc_mk( nbin ) )
  allocate( xctr_mk( nbin ) )
  allocate( xbnd_mk( nbin+1 ) )
  allocate( cctr_mk( nspc_mk,nbin ) )
  allocate( cbnd_mk( nspc_mk,nbin+1 ) )
  allocate( ck_mk( nspc_mk,nspc_mk,nbin,nbin ) )
  allocate( vt_mk( nspc_mk,nbin ) )
  allocate( br_mk( nspc_mk,nbin ) )

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

  deallocate( radc_mk )
  deallocate( xctr_mk )
  deallocate( xbnd_mk )
  deallocate( cctr_mk )
  deallocate( cbnd_mk )
  deallocate( ck_mk )
  deallocate( vt_mk )
  deallocate( br_mk )

  end subroutine mkpara
  !---------------------------------------------------------------------------------------
  subroutine rdkdat

  implicit none
  integer, parameter :: il = 1, ic = 2, ip = 3, id = 4
  integer, parameter :: is = 5, ig = 6, ih = 7
  integer, parameter :: icemax = 3

  integer :: k, kk, i, j
  real(DP) :: xl( ndat ), rlec( ndat ), vrl( ndat )
  real(DP) :: blkradl( ndat ), blkdnsl( ndat )

  real(DP) :: xi( ndat,icemax ), riec( ndat,icemax ), vri( ndat,icemax )
  real(DP) :: blkradi( ndat,icemax ), blkdnsi( ndat,icemax )

  real(DP) :: xs( ndat ), rsec( ndat ), vrs( ndat )
  real(DP) :: blkrads( ndat ), blkdnss( ndat )

  real(DP) :: xg( ndat ), rgec( ndat ), vrg( ndat )
  real(DP) :: blkradg( ndat ), blkdnsg( ndat )

  real(DP) :: xh( ndat ), rhec( ndat ), vrh( ndat )
  real(DP) :: blkradh( ndat ), blkdnsh( ndat )

  data xl(1:ndat) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00_DP /
  data xi(1:ndat,1) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00_DP /
  data xi(1:ndat,2) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00 /
  data xi(1:ndat,3) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00_DP /
  data xs(1:ndat) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00_DP /
  data xg(1:ndat) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00_DP /
  data xh(1:ndat) / &
     0.33510E-10_DP,0.67021E-10_DP,0.13404E-09_DP,0.26808E-09_DP,0.53617E-09_DP,0.10723E-08_DP, &
     0.21447E-08_DP,0.42893E-08_DP,0.85786E-08_DP,0.17157E-07_DP,0.34315E-07_DP,0.68629E-07_DP, &
     0.13726E-06_DP,0.27452E-06_DP,0.54903E-06_DP,0.10981E-05_DP,0.21961E-05_DP,0.43923E-05_DP, &
     0.87845E-05_DP,0.17569E-04_DP,0.35138E-04_DP,0.70276E-04_DP,0.14055E-03_DP,0.28110E-03_DP, &
     0.56221E-03_DP,0.11244E-02_DP,0.22488E-02_DP,0.44977E-02_DP,0.89954E-02_DP,0.17991E-01_DP, &
     0.35981E-01_DP,0.71963E-01_DP,0.14393E+00_DP /
  data rlec(1:ndat) / &
     0.20000E-03_DP,0.25198E-03_DP,0.31748E-03_DP,0.40000E-03_DP,0.50397E-03_DP,0.63496E-03_DP, &
     0.80000E-03_DP,0.10079E-02_DP,0.12699E-02_DP,0.16000E-02_DP,0.20159E-02_DP,0.25398E-02_DP, &
     0.32000E-02_DP,0.40317E-02_DP,0.50797E-02_DP,0.64000E-02_DP,0.80635E-02_DP,0.10159E-01_DP, &
     0.12800E-01_DP,0.16127E-01_DP,0.20319E-01_DP,0.25600E-01_DP,0.32254E-01_DP,0.40637E-01_DP, &
     0.51200E-01_DP,0.64508E-01_DP,0.81275E-01_DP,0.10240E+00_DP,0.12902E+00_DP,0.16255E+00_DP, &
     0.20480E+00_DP,0.25803E+00_DP,0.32510E+00_DP /
  data riec(1:ndat,1) / &
     0.31936E-03_DP,0.40397E-03_DP,0.51099E-03_DP,0.64638E-03_DP,0.81764E-03_DP,0.10343E-02_DP, &
     0.13084E-02_DP,0.16551E-02_DP,0.20937E-02_DP,0.26486E-02_DP,0.33506E-02_DP,0.42387E-02_DP, &
     0.64360E-02_DP,0.81426E-02_DP,0.10302E-01_DP,0.13035E-01_DP,0.16494E-01_DP,0.20872E-01_DP, &
     0.26412E-01_DP,0.33426E-01_DP,0.42304E-01_DP,0.53543E-01_DP,0.67770E-01_DP,0.85783E-01_DP, &
     0.10859E+00_DP,0.13746E+00_DP,0.17403E+00_DP,0.22032E+00_DP,0.27895E+00_DP,0.35319E+00_DP, &
     0.44722E+00_DP,0.56630E+00_DP,0.71712E+00_DP /
  data riec(1:ndat,2) / &
     0.13188E-03_DP,0.16615E-03_DP,0.20953E-03_DP,0.27728E-03_DP,0.36694E-03_DP,0.48559E-03_DP, &
     0.64261E-03_DP,0.85040E-03_DP,0.11254E-02_DP,0.14893E-02_DP,0.19709E-02_DP,0.26082E-02_DP, &
     0.34515E-02_DP,0.45676E-02_DP,0.60446E-02_DP,0.79991E-02_DP,0.10586E-01_DP,0.14009E-01_DP, &
     0.18539E-01_DP,0.24533E-01_DP,0.32466E-01_DP,0.42964E-01_DP,0.56857E-01_DP,0.75242E-01_DP, &
     0.99573E-01_DP,0.13177E+00_DP,0.17438E+00_DP,0.23077E+00_DP,0.30539E+00_DP,0.40414E+00_DP, &
     0.53482E+00_DP,0.70775E+00_DP,0.93661E+00_DP /
  data riec(1:ndat,3) /            0.14998E-03_DP,0.18896E-03_DP,0.23808E-03_DP, &
     0.29996E-03_DP,0.37793E-03_DP,0.47616E-03_DP,0.61048E-03_DP,0.81343E-03_DP,0.10839E-02_DP, &
     0.14442E-02_DP,0.19243E-02_DP,0.25640E-02_DP,0.34164E-02_DP,0.45522E-02_DP,0.60656E-02_DP, &
     0.80820E-02_DP,0.10769E-01_DP,0.14349E-01_DP,0.19119E-01_DP,0.25475E-01_DP,0.44576E-01_DP, &
     0.62633E-01_DP,0.88006E-01_DP,0.12366E+00_DP,0.17375E+00_DP,0.24414E+00_DP,0.34304E+00_DP, &
     0.48201E+00_DP,0.67728E+00_DP,0.95164E+00_DP,0.13372E+01_DP,0.18788E+01_DP,0.26400E+01_DP /
  data rsec(1:ndat) / &
     0.92832E-03_DP,0.11696E-02_DP,0.14736E-02_DP,0.18566E-02_DP,0.23392E-02_DP,0.29472E-02_DP, &
     0.37133E-02_DP,0.46784E-02_DP,0.58944E-02_DP,0.74265E-02_DP,0.93569E-02_DP,0.11789E-01_DP, &
     0.14853E-01_DP,0.18714E-01_DP,0.23578E-01_DP,0.29706E-01_DP,0.37427E-01_DP,0.47156E-01_DP, &
     0.59412E-01_DP,0.74855E-01_DP,0.94311E-01_DP,0.11882E+00_DP,0.14971E+00_DP,0.18862E+00_DP, &
     0.23765E+00_DP,0.29942E+00_DP,0.37724E+00_DP,0.47530E+00_DP,0.59884E+00_DP,0.75449E+00_DP, &
     0.95060E+00_DP,0.11977E+01_DP,0.15090E+01_DP /
  data rgec(1:ndat) /              0.27144E-03_DP,0.34200E-03_DP,0.43089E-03_DP, &
     0.54288E-03_DP,0.68399E-03_DP,0.86177E-03_DP,0.10858E-02_DP,0.13680E-02_DP,0.17235E-02_DP, &
     0.21715E-02_DP,0.27360E-02_DP,0.34471E-02_DP,0.43431E-02_DP,0.54719E-02_DP,0.68942E-02_DP, &
     0.86861E-02_DP,0.10944E-01_DP,0.13788E-01_DP,0.17372E-01_DP,0.21888E-01_DP,0.27577E-01_DP, &
     0.34745E-01_DP,0.43775E-01_DP,0.55154E-01_DP,0.69489E-01_DP,0.87551E-01_DP,0.11031E+00_DP, &
     0.13898E+00_DP,0.17510E+00_DP,0.22061E+00_DP,0.27796E+00_DP,0.35020E+00_DP,0.44123E+00_DP /
  data rhec(1:ndat) / &
     0.20715E-03_DP,0.26099E-03_DP,0.32883E-03_DP,0.41430E-03_DP,0.52198E-03_DP,0.65766E-03_DP, &
     0.82860E-03_DP,0.10440E-02_DP,0.13153E-02_DP,0.16572E-02_DP,0.20879E-02_DP,0.26306E-02_DP, &
     0.33144E-02_DP,0.41759E-02_DP,0.52613E-02_DP,0.66288E-02_DP,0.83517E-02_DP,0.10523E-01_DP, &
     0.13258E-01_DP,0.16703E-01_DP,0.21045E-01_DP,0.26515E-01_DP,0.33407E-01_DP,0.42090E-01_DP, &
     0.53030E-01_DP,0.66814E-01_DP,0.84180E-01_DP,0.10606E+00_DP,0.13363E+00_DP,0.16836E+00_DP, &
     0.21212E+00_DP,0.26725E+00_DP,0.33672E+00_DP /
  data vrl(1:ndat) / &
     0.50000E-01_DP,0.78000E-01_DP,0.12000E+00_DP,0.19000E+00_DP,0.31000E+00_DP,0.49000E+00_DP, &
     0.77000E+00_DP,0.12000E+01_DP,0.19000E+01_DP,0.30000E+01_DP,0.48000E+01_DP,0.74000E+01_DP, &
     0.11000E+02_DP,0.17000E+02_DP,0.26000E+02_DP,0.37000E+02_DP,0.52000E+02_DP,0.71000E+02_DP, &
     0.94000E+02_DP,0.12000E+03_DP,0.16000E+03_DP,0.21000E+03_DP,0.26000E+03_DP,0.33000E+03_DP, &
     0.41000E+03_DP,0.48000E+03_DP,0.57000E+03_DP,0.66000E+03_DP,0.75000E+03_DP,0.82000E+03_DP, &
     0.88000E+03_DP,0.90000E+03_DP,0.90000E+03_DP /
  data vri(1:ndat,1) /             0.30000E-01_DP,0.40000E-01_DP,0.60000E-01_DP, &
     0.80000E-01_DP,0.11000E+00_DP,0.15000E+00_DP,0.17000E+00_DP,0.18000E+00_DP,0.20000E+00_DP, &
     0.25000E+00_DP,0.40000E+00_DP,0.60000E+01_DP,0.10000E+02_DP,0.15000E+02_DP,0.20000E+02_DP, &
     0.25000E+02_DP,0.31000E+02_DP,0.37000E+02_DP,0.41000E+02_DP,0.46000E+02_DP,0.51000E+02_DP, &
     0.55000E+02_DP,0.59000E+02_DP,0.62000E+02_DP,0.64000E+02_DP,0.67000E+02_DP,0.68000E+02_DP, &
     0.69000E+02_DP,0.70000E+02_DP,0.71000E+02_DP,0.71500E+02_DP,0.71750E+02_DP,0.72000E+02_DP /
  data vri(1:ndat,2) / &
     0.30000E-01_DP,0.40000E-01_DP,0.50000E-01_DP,0.70000E-01_DP,0.90000E-01_DP,0.12000E+00_DP, &
     0.50000E+00_DP,0.80000E+00_DP,0.16000E+01_DP,0.18000E+01_DP,0.20000E+01_DP,0.30000E+01_DP, &
     0.40000E+01_DP,0.50000E+01_DP,0.80000E+01_DP,0.13000E+02_DP,0.19000E+02_DP,0.26000E+02_DP, &
     0.32000E+02_DP,0.38000E+02_DP,0.47000E+02_DP,0.55000E+02_DP,0.65000E+02_DP,0.73000E+02_DP, &
     0.77000E+02_DP,0.79000E+02_DP,0.80000E+02_DP,0.81000E+02_DP,0.81000E+02_DP,0.82000E+02_DP, &
     0.82000E+02_DP,0.82000E+02_DP,0.82000E+02_DP /
  data vri(1:ndat,3) /             0.35000E-01_DP,0.45000E-01_DP,0.55000E-01_DP, &
     0.75000E-01_DP,0.95000E-01_DP,0.13000E+00_DP,0.60000E+00_DP,0.90000E+00_DP,0.17000E+01_DP, &
     0.20000E+01_DP,0.25000E+01_DP,0.38000E+01_DP,0.50000E+01_DP,0.70000E+01_DP,0.90000E+01_DP, &
     0.11000E+02_DP,0.14000E+02_DP,0.17000E+02_DP,0.21000E+02_DP,0.25000E+02_DP,0.32000E+02_DP, &
     0.38000E+02_DP,0.44000E+02_DP,0.49000E+02_DP,0.53000E+02_DP,0.55000E+02_DP,0.58000E+02_DP, &
     0.59000E+02_DP,0.61000E+02_DP,0.62000E+02_DP,0.63000E+02_DP,0.64000E+02_DP,0.65000E+02_DP /
  data vrs(1:ndat) / &
     0.20000E-01_DP,0.31000E-01_DP,0.49000E-01_DP,0.77000E-01_DP,0.12000E+00_DP,0.19000E+00_DP, &
     0.30000E+00_DP,0.48000E+00_DP,0.76000E+00_DP,0.12000E+01_DP,0.19000E+01_DP,0.30000E+01_DP, &
     0.48000E+01_DP,0.75000E+01_DP,0.11000E+02_DP,0.16000E+02_DP,0.21000E+02_DP,0.26000E+02_DP, &
     0.34000E+02_DP,0.41000E+02_DP,0.49000E+02_DP,0.57000E+02_DP,0.65000E+02_DP,0.73000E+02_DP, &
     0.81000E+02_DP,0.87000E+02_DP,0.93000E+02_DP,0.99000E+02_DP,0.10750E+03_DP,0.11500E+03_DP, &
     0.12500E+03_DP,0.13500E+03_DP,0.14500E+03_DP /
  data vrg(1:ndat) /               0.39000E-01_DP,0.62000E-01_DP,0.97000E-01_DP, &
     0.15000E+00_DP,0.24000E+00_DP,0.38000E+00_DP,0.61000E+00_DP,0.96000E+00_DP,0.15000E+01_DP, &
     0.24000E+01_DP,0.38000E+01_DP,0.61000E+01_DP,0.96000E+01_DP,0.15000E+02_DP,0.23000E+02_DP, &
     0.31000E+02_DP,0.39000E+02_DP,0.49000E+02_DP,0.59000E+02_DP,0.68000E+02_DP,0.79000E+02_DP, &
     0.88000E+02_DP,0.10000E+03_DP,0.11000E+03_DP,0.13000E+03_DP,0.15000E+03_DP,0.17000E+03_DP, &
     0.20000E+03_DP,0.23000E+03_DP,0.26000E+03_DP,0.30000E+03_DP,0.35000E+03_DP,0.40000E+03_DP  /
  data vrh(1:ndat) / &
     0.53000E-01_DP,0.84000E-01_DP,0.13000E+00_DP,0.21000E+00_DP,0.33000E+00_DP,0.52000E+00_DP, &
     0.82000E+00_DP,0.13000E+01_DP,0.21000E+01_DP,0.33000E+01_DP,0.52000E+01_DP,0.82000E+01_DP, &
     0.13000E+02_DP,0.20000E+02_DP,0.28000E+02_DP,0.36000E+02_DP,0.46000E+02_DP,0.56000E+02_DP, &
     0.67000E+02_DP,0.80000E+02_DP,0.97000E+02_DP,0.12000E+03_DP,0.14000E+03_DP,0.17000E+03_DP, &
     0.20000E+03_DP,0.24000E+03_DP,0.29000E+03_DP,0.35000E+03_DP,0.42000E+03_DP,0.51000E+03_DP, &
     0.61000E+03_DP,0.74000E+03_DP,0.89000E+03_DP /
  data blkradl(1:ndat) / &
     0.20000E-03_DP,0.25198E-03_DP,0.31748E-03_DP,0.40000E-03_DP,0.50397E-03_DP,0.63496E-03_DP, &
     0.80000E-03_DP,0.10079E-02_DP,0.12699E-02_DP,0.16000E-02_DP,0.20159E-02_DP,0.25398E-02_DP, &
     0.32000E-02_DP,0.40317E-02_DP,0.50797E-02_DP,0.64000E-02_DP,0.80635E-02_DP,0.10159E-01_DP, &
     0.12800E-01_DP,0.16127E-01_DP,0.20319E-01_DP,0.25600E-01_DP,0.32254E-01_DP,0.40637E-01_DP, &
     0.51200E-01_DP,0.64508E-01_DP,0.81275E-01_DP,0.10240E+00_DP,0.12902E+00_DP,0.16255E+00_DP, &
     0.20480E+00_DP,0.25803E+00_DP,0.32510E+00_DP /
  data blkradi(1:ndat,1) /         0.57452E-03_DP,0.72384E-03_DP,0.91199E-03_DP, &
     0.11490E-02_DP,0.14477E-02_DP,0.18240E-02_DP,0.22981E-02_DP,0.28954E-02_DP,0.36479E-02_DP, &
     0.45961E-02_DP,0.57908E-02_DP,0.72959E-02_DP,0.11572E-01_DP,0.14770E-01_DP,0.18851E-01_DP, &
     0.24060E-01_DP,0.30709E-01_DP,0.39194E-01_DP,0.50025E-01_DP,0.63848E-01_DP,0.81491E-01_DP, &
     0.10401E+00_DP,0.13275E+00_DP,0.16943E+00_DP,0.21625E+00_DP,0.27601E+00_DP,0.35228E+00_DP, &
     0.44962E+00_DP,0.57387E+00_DP,0.73244E+00_DP,0.93484E+00_DP,0.11932E+01_DP,0.15229E+01_DP /
  data blkradi(1:ndat,2) / &
     0.20715E-03_DP,0.26099E-03_DP,0.32912E-03_DP,0.43555E-03_DP,0.57638E-03_DP,0.76276E-03_DP, &
     0.10094E-02_DP,0.13358E-02_DP,0.17677E-02_DP,0.23394E-02_DP,0.30958E-02_DP,0.40969E-02_DP, &
     0.54216E-02_DP,0.71748E-02_DP,0.94948E-02_DP,0.12565E-01_DP,0.16628E-01_DP,0.22005E-01_DP, &
     0.29120E-01_DP,0.38537E-01_DP,0.50998E-01_DP,0.67488E-01_DP,0.89311E-01_DP,0.11819E+00_DP, &
     0.15641E+00_DP,0.20698E+00_DP,0.27391E+00_DP,0.36249E+00_DP,0.47970E+00_DP,0.63481E+00_DP, &
     0.84009E+00_DP,0.11117E+01_DP,0.14712E+01_DP /
  data blkradi(1:ndat,3) /         0.23559E-03_DP,0.29682E-03_DP,0.37397E-03_DP, &
     0.47118E-03_DP,0.59365E-03_DP,0.74795E-03_DP,0.95894E-03_DP,0.12777E-02_DP,0.17025E-02_DP, &
     0.22685E-02_DP,0.30227E-02_DP,0.40275E-02_DP,0.53665E-02_DP,0.71506E-02_DP,0.95278E-02_DP, &
     0.12695E-01_DP,0.16916E-01_DP,0.22539E-01_DP,0.30032E-01_DP,0.40017E-01_DP,0.70019E-01_DP, &
     0.98384E-01_DP,0.13824E+00_DP,0.19424E+00_DP,0.27293E+00_DP,0.38350E+00_DP,0.53885E+00_DP, &
     0.75714E+00_DP,0.10639E+01_DP,0.14948E+01_DP,0.21004E+01_DP,0.29513E+01_DP,0.41469E+01_DP /
  data blkrads(1:ndat) / &
     0.20715E-03_DP,0.26148E-03_DP,0.33067E-03_DP,0.41710E-03_DP,0.52691E-03_DP,0.66640E-03_DP, &
     0.84289E-03_DP,0.10674E-02_DP,0.13513E-02_DP,0.17129E-02_DP,0.21670E-02_DP,0.27521E-02_DP, &
     0.34989E-02_DP,0.44777E-02_DP,0.57347E-02_DP,0.75389E-02_DP,0.99020E-02_DP,0.13161E-01_DP, &
     0.17372E-01_DP,0.23337E-01_DP,0.31058E-01_DP,0.41194E-01_DP,0.55153E-01_DP,0.74854E-01_DP, &
     0.99806E-01_DP,0.13463E+00_DP,0.18136E+00_DP,0.24282E+00_DP,0.32955E+00_DP,0.44123E+00_DP, &
     0.59884E+00_DP,0.77090E+00_DP,0.99387E+00_DP /
  data blkradg(1:ndat) /           0.27144E-03_DP,0.34200E-03_DP,0.43089E-03_DP, &
     0.54288E-03_DP,0.68399E-03_DP,0.86177E-03_DP,0.10858E-02_DP,0.13680E-02_DP,0.17235E-02_DP, &
     0.21715E-02_DP,0.27360E-02_DP,0.34471E-02_DP,0.43431E-02_DP,0.54719E-02_DP,0.68942E-02_DP, &
     0.86861E-02_DP,0.10944E-01_DP,0.13788E-01_DP,0.17372E-01_DP,0.21888E-01_DP,0.27577E-01_DP, &
     0.34745E-01_DP,0.43775E-01_DP,0.55154E-01_DP,0.69489E-01_DP,0.87551E-01_DP,0.11031E+00_DP, &
     0.13898E+00_DP,0.17510E+00_DP,0.22061E+00_DP,0.27796E+00_DP,0.35020E+00_DP,0.44123E+00_DP /
  data blkradh(1:ndat) / &
     0.20715E-03_DP,0.26099E-03_DP,0.32883E-03_DP,0.41430E-03_DP,0.52198E-03_DP,0.65766E-03_DP, &
     0.82860E-03_DP,0.10440E-02_DP,0.13153E-02_DP,0.16572E-02_DP,0.20879E-02_DP,0.26306E-02_DP, &
     0.33144E-02_DP,0.41759E-02_DP,0.52613E-02_DP,0.66288E-02_DP,0.83517E-02_DP,0.10523E-01_DP, &
     0.13258E-01_DP,0.16703E-01_DP,0.21045E-01_DP,0.26515E-01_DP,0.33407E-01_DP,0.42090E-01_DP, &
     0.53030E-01_DP,0.66814E-01_DP,0.84180E-01_DP,0.10606E+00_DP,0.13363E+00_DP,0.16836E+00_DP, &
     0.21212E+00_DP,0.26725E+00_DP,0.33672E+00_DP /
  data blkdnsl(1:ndat) / &
     0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP, &
     0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP, &
     0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP, &
     0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP, &
     0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP, &
     0.10000E+01_DP,0.10000E+01_DP,0.10000E+01_DP /
  data blkdnsi(1:ndat,1) /         0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.87368E+00_DP,0.87072E+00_DP,0.86777E+00_DP, &
     0.86483E+00_DP,0.86189E+00_DP,0.85897E+00_DP,0.85606E+00_DP,0.85316E+00_DP,0.85026E+00_DP, &
     0.84738E+00_DP,0.84451E+00_DP,0.84164E+00_DP,0.83879E+00_DP,0.83595E+00_DP,0.83311E+00_DP, &
     0.83029E+00_DP,0.82747E+00_DP,0.82467E+00_DP,0.82187E+00_DP,0.81908E+00_DP,0.81631E+00_DP /
  data blkdnsi(1:ndat,2) /  &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP /
  data blkdnsi(1:ndat,3) /         0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP, &
     0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP, &
     0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP, &
     0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.61183E+00_DP,0.51790E+00_DP, &
     0.45557E+00_DP,0.40075E+00_DP,0.35252E+00_DP,0.31010E+00_DP,0.27278E+00_DP,0.23995E+00_DP, &
     0.21108E+00_DP,0.18567E+00_DP,0.16333E+00_DP,0.14367E+00_DP,0.12638E+00_DP,0.11118E+00_DP /
  data blkdnss(1:ndat) / &
     0.90000E+00_DP,0.89500E+00_DP,0.88500E+00_DP,0.88200E+00_DP,0.87500E+00_DP,0.86500E+00_DP, &
     0.85500E+00_DP,0.84200E+00_DP,0.83000E+00_DP,0.81500E+00_DP,0.80500E+00_DP,0.78600E+00_DP, &
     0.76500E+00_DP,0.73000E+00_DP,0.69500E+00_DP,0.61183E+00_DP,0.54000E+00_DP,0.46000E+00_DP, &
     0.40000E+00_DP,0.33000E+00_DP,0.28000E+00_DP,0.24000E+00_DP,0.20000E+00_DP,0.16000E+00_DP, &
     0.13500E+00_DP,0.11000E+00_DP,0.90000E-01_DP,0.75000E-01_DP,0.60000E-01_DP,0.50000E-01_DP, &
     0.40000E-01_DP,0.37500E-01_DP,0.35000E-01_DP /
  data blkdnsg(1:ndat) / &
     0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP, &
     0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP, &
     0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP, &
     0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP, &
     0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP, &
     0.40000E+00_DP,0.40000E+00_DP,0.40000E+00_DP /
  data blkdnsh(1:ndat) / &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP, &
     0.90000E+00_DP,0.90000E+00_DP,0.90000E+00_DP /

  !--- mass
  do k = 1, ndat
    xmss( il,k ) = log( dble( xl( k ) )*1.E-03_DP )
    xmss( ic,k ) = log( dble( xi( k,1 ) )*1.E-03_DP )
    xmss( ip,k ) = log( dble( xi( k,2 ) )*1.E-03_DP )
    xmss( id,k ) = log( dble( xi( k,3 ) )*1.E-03_DP )
    xmss( is,k ) = log( dble( xs( k ) )*1.E-03_DP )
    xmss( ig,k ) = log( dble( xg( k ) )*1.E-03_DP )
    xmss( ih,k ) = log( dble( xh( k ) )*1.E-03_DP )
  end do

  !--- capacity
  do k = 1, ndat ! cm -> m
    zcap( il,k ) = dble(rlec( k ))*1.E-02_DP
    zcap( ic,k ) = dble(riec( k,1 ))*1.E-02_DP
    zcap( ip,k ) = dble(riec( k,2 ))*1.E-02_DP
    zcap( id,k ) = dble(riec( k,3 ))*1.E-02_DP
    zcap( is,k ) = dble(rsec( k ))*1.E-02_DP
    zcap( ig,k ) = dble(rgec( k ))*1.E-02_DP
    zcap( ih,k ) = dble(rhec( k ))*1.E-02_DP
  end do

  !--- terminal velocity
  do k = 1, ndat ! cm/s -> m/s
    vtrm( il,k ) = dble(vrl( k ))*1.E-02_DP
    vtrm( ic,k ) = dble(vri( k,1 ))*1.E-02_DP
    vtrm( ip,k ) = dble(vri( k,2 ))*1.E-02_DP
    vtrm( id,k ) = dble(vri( k,3 ))*1.E-02_DP
    vtrm( is,k ) = dble(vrs( k ))*1.E-02_DP
    vtrm( ig,k ) = dble(vrg( k ))*1.E-02_DP
    vtrm( ih,k ) = dble(vrh( k ))*1.E-02_DP
  end do

  !--- bulk radii
  do k = 1, ndat ! cm -> mm
    blkr( il,k ) = dble(blkradl( k ))*10.0_DP
    blkr( ic,k ) = dble(blkradi( k,1 ))*10.0_DP
    blkr( ip,k ) = dble(blkradi( k,2 ))*10.0_DP
    blkr( id,k ) = dble(blkradi( k,3 ))*10.0_DP
    blkr( is,k ) = dble(blkrads( k ))*10.0_DP
    blkr( ig,k ) = dble(blkradg( k ))*10.0_DP
    blkr( ih,k ) = dble(blkradh( k ))*10.0_DP
  end do

  !--- bulk density
  do k = 1, ndat ! g/cm^3 -> kg/m^3
    blkd( il,k ) = dble(blkdnsl( k ))*1000._DP
    blkd( ic,k ) = dble(blkdnsi( k,1 ))*1000._DP
    blkd( ip,k ) = dble(blkdnsi( k,2 ))*1000._DP
    blkd( id,k ) = dble(blkdnsi( k,3 ))*1000._DP
    blkd( is,k ) = dble(blkdnss( k ))*1000._DP
    blkd( ig,k ) = dble(blkdnsg( k ))*1000._DP
    blkd( ih,k ) = dble(blkdnsh( k ))*1000._DP
  end do

  !--- collection kernel
  ! cm**3/s -> m**3/s
  do i = 1, ndat
  do j = 1, ndat
    do k = 1, 7
    do kk = 1, 7
     ykrn( k,kk,i,j ) = kernels( k,kk,i,j )*1.E-06_DP
    enddo
    enddo

  end do
  end do

  end subroutine rdkdat
  !---------------------------------------------------------------------------------------
  subroutine sdfgrid

  real(DP) :: xsta, xend
  integer :: n
  real(DP):: pi_dp = 3.1415920_DP
  real(DP):: rhow_dp = 1.0E+03_DP

  xsta = log( rhow_dp * 4.0_DP*pi_dp/3.0_DP * ( 3.E-06_DP )**3.0_DP )
  xend = log( rhow_dp * 4.0_DP*pi_dp/3.0_DP * ( 3.E-03_DP )**3.0_DP )

  dxmic_mk = ( xend-xsta )/nbin
  dxmic = dxmic_mk
  do n = 1, nbin+1
    xbnd_mk( n ) = xsta + dxmic_mk*( n-1 )
  end do
  do n = 1, nbin
    xctr_mk( n ) = ( xbnd_mk( n )+xbnd_mk( n+1 ) )*0.50_DP
    radc_mk( n ) = ( exp( xctr_mk( n ) )*3.0_DP/4.0_DP/pi_dp/rhow_dp )**( 1.0_DP/3.0_DP )
  end do

  end subroutine sdfgrid
  !---------------------------------------------------------------------------------------
  subroutine getcp

  integer :: myu, n

  do myu = 1, nspc_mk
    do n = 1, nbin
      cctr_mk( myu,n ) = fcpc( myu,xctr_mk( n ) )
    end do
    do n = 1, nbin+1
      cbnd_mk( myu,n ) = fcpc( myu,xbnd_mk( n ) )
    end do
  end do

  end subroutine getcp
  !---------------------------------------------------------------------------------------
  function fcpc( myu,x )

  integer, intent(in) :: myu
  real(DP), intent(in) :: x
  real(DP) :: fcpc

  real(DP) :: qknt( ndat+kdeg ), elm( ndat,ndat ), coef( ndat )

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

  if ( IO_L ) write(IO_FID_LOG,*) 'Create micpara.dat'
  if( kphase == 0 ) then
   if ( IO_L ) write(IO_FID_LOG,*) 'Hydro-dynamic kernel'
  else if( kphase == 1 ) then
   if ( IO_L ) write(IO_FID_LOG,*) 'Long Kernel'
  else if( kphase == 2 ) then
   if ( IO_L ) write(IO_FID_LOG,*) 'Golovin Kernel'
  endif

  do myu = 1, nspc_mk
  do nyu = 1, nspc_mk
  if ( IO_L ) write(IO_FID_LOG,*) ' myu, nyu :', myu, nyu
  do i = 1, nbin
  do j = 1, nbin
    ck_mk( myu,nyu,i,j ) = fckrn( myu,nyu,xctr_mk( i ),xctr_mk( j ) )
  end do
  end do
  end do
  end do

  return

  end subroutine getck
  !---------------------------------------------------------------------------------------
  function fckrn( myu,nyu,x,y )

  integer, intent(in) :: myu, nyu
  real(DP), intent(in) :: x, y
  real(DP) :: fckrn

  real(DP) :: qknt( ndat+kdeg ), rknt( ndat+kdeg )
  real(DP) :: coef( ndat,ndat )

  real(DP) :: xlrg, xsml, vlrg, vsml, rlrg
  real(DP):: pi_dp = 3.1415920_DP
  real(DP):: rhow_dp = 1.0E+03_DP

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

   vlrg = (exp( xlrg ) / rhow_dp )*1.E+06_DP
   vsml = (exp( xsml ) / rhow_dp )*1.E+06_DP

   rlrg = ( exp( xlrg )/( 4.0_DP*pi_dp*rhow_dp )*3.0_DP )**(1.0_DP/3.0_DP )*1.E+06_DP

   if( rlrg <=50.0_DP ) then
     fckrn = 9.44E+03_DP*( vlrg*vlrg + vsml*vsml )
   else
     fckrn = 5.78E-03_DP*( vlrg+vsml )
   end if
  else if( kphase == 2 ) then
   fckrn = 1.5_DP*( exp(x) +exp(y) )
  end if

  return

  end function fckrn
  !---------------------------------------------------------------------------------------
  subroutine getvt

  integer :: myu, n

  do myu = 1, nspc_mk
  do n = 1, nbin
    vt_mk( myu,n ) = max( fvterm( myu,xctr_mk( n ) ), 0.0_DP )
  end do
  end do

  end subroutine getvt
  !---------------------------------------------------------------------------------------
  function fvterm( myu,x )

  integer, intent(in) :: myu
  real(DP), intent(in) :: x
  real(DP) :: fvterm

  real(DP) :: qknt( ndat+kdeg ), elm( ndat,ndat ), coef( ndat )

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

  do myu = 1, nspc_mk
  do n = 1, nbin
    br_mk( myu,n ) = fbulkrad( myu, xctr_mk( n ) )
  end do
  end do

  end subroutine getbr
  !---------------------------------------------------------------------------------------
  function fbulkrad( myu,x )

  integer, intent(in) :: myu
  real(DP), intent(in) :: x
  real(DP) :: fbulkrad

  real(DP) :: qknt( ndat+kdeg ), elm( ndat,ndat ), coef( ndat )

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

  open ( fid_micpara, file = fname_micpara, form = 'formatted', status='new' )

  write( fid_micpara,* ) nspc_mk, nbin

  ! grid parameter
  do n = 1, nbin
    xctr( n ) = xctr_mk( n )
    radc( n ) = radc_mk( n )
    write( fid_micpara,* ) n, xctr( n ), radc( n )
  end do
  do n = 1, nbin+1
    xbnd( n ) = xbnd_mk( n )
    write( fid_micpara,* ) n, xbnd( n )
  end do
  write( fid_micpara,* ) dxmic_mk

  ! capacity
  do myu = 1, nspc_mk
    do n = 1, nbin
      cctr( myu,n ) = cctr_mk( myu,n )
      write( fid_micpara,* ) myu, n, cctr( myu,n )
    end do
    do n = 1, nbin+1
      cbnd( myu,n ) = cbnd_mk( myu,n )
      write( fid_micpara,* ) myu, n, cbnd( myu,n )
    end do
  end do

  ! collection kernel
  do myu = 1, nspc_mk
  do nyu = 1, nspc_mk
  do i = 1, nbin
  do j = 1, nbin
    ck( myu,nyu,i,j ) = ck_mk( myu,nyu,i,j )
    write( fid_micpara,* ) myu, nyu, i, j, ck( myu,nyu,i,j )
  end do
  end do
  end do
  end do

  ! falling velocity
  do myu = 1, nspc_mk
  do n = 1, nbin
    vt( myu,n ) = vt_mk( myu,n )
    write( fid_micpara,* ) myu, n, vt( myu,n )
  end do
  end do

  ! bulk radius
  do myu = 1, nspc_mk
  do n = 1, nbin
    br( myu,n ) = br_mk( myu,n )
    write( fid_micpara,* ) myu, n, br( myu,n )
  end do
  end do

  close ( fid_micpara )

  end subroutine paraout
  !---------------------------------------------------------------------------------------
  !---- unify from other files
  !---------------------------------------------------------------------------------------
  subroutine TINVSS(n,a,dt,e,nn,iw,inder)

  implicit none

  integer, intent(in)    :: n, nn
  integer, intent(inout) :: inder
  real(DP), intent(inout) :: a(nn,n)
  real(DP), intent(inout) :: dt
  real(DP), intent(in)    :: e
  integer, intent(inout) :: iw( 2*n )
  integer :: i, j, k, kk, ij, nnn
  real(DP) :: work, aa, az, eps
  integer :: ipiv, jpiv
  real(DP) :: piv
  !-----------------------------------------

    inder = 0
    if( n < 1 ) then
      goto 910
    elseif( n == 1 ) then
      goto 930
    elseif( n > 1 ) then
      goto 101
    endif

101 continue
    if( n > nn ) then
     inder = -1
     write(6,690) "n= ", n, "nn= ", nn, &
           "n should be less than or equal to nn in TINVSS"
690 format(a8,i5,5x,a4,i5,a55)
     return
    endif

    eps = 0.0_DP
    dt  = 1.0_DP
    do k = 1, n
     piv = 0.0_DP
     do i = k, n
     do j = k, n
      if( abs(a(i,j)) <= abs(piv) ) goto 110   !
      ipiv = i
      jpiv = j
      piv = a(i,j)
110   continue
     enddo
     enddo
     dt = dt * piv
     if( abs(piv) <= eps )  goto 920
     if( k == 1 ) eps = abs(piv)*e
     if( ipiv == k ) goto 130
     dt = -dt
     do j = 1, n
      work = a(ipiv,j)
      a(ipiv,j) = a(k,j)
      a(k,j) = work
     enddo
130  continue
     if( jpiv == k ) goto 150
     dt = -dt
     do i = 1, n
      work = a(i,jpiv)
      a(i,jpiv) = a(i,k)
      a(i,k) = work
     enddo
150  continue
     iw(2*k-1) = ipiv
     aa=1.0_DP/piv
     iw(2*k) = jpiv
     do j = 1, n
      a(k,j) = a(k,j)*aa
     enddo
     do i = 1, n
      if( i == k ) goto 220
      az = a(i,k)
      if( az == 0.0_DP ) goto 220
      do j = 1, n
       a(i,j) = a(i,j)-a(k,j)*az
      enddo
      a(i,k) = -aa*az
220   continue
     enddo
     a(k,k) = aa
    enddo
    do kk = 2, n
      k=n+1-kk
      ij=iw(2*k)
      if( ij == k ) goto 420
      do j = 1, n
       work=a(ij,j)
       a(ij,j) = a(k,j)
       a(k,j) = work
      enddo
420   continue
      ij = iw(2*k-1)
      if( ij == k ) goto 400
      do i = 1, n
       work=a(i,ij)
       a(i,ij)=a(i,k)
       a(i,k)=work
      enddo
400   continue
    enddo

    return

910 continue
    inder = -1
    write(*,691) "n= ", n, "should be positive in TINVSS"
691 format(a8,i5,5x,a30)
    return


920 continue
    dt = 0.0_DP
    inder = n-k+1
    nnn = k-1
    write(*,692) 'given matrix A to TINVSS is ill conditioned, or sigular withrank =', nnn, &
                  'return with no further calculation'
692 format(a,1x,i4,1x,a)
    return

930 continue
    dt=a(1,1)
    k=1
    if( dt == 0.0_DP ) goto 920
    a(1,1) = 1.0_DP/a(1,1)
    return

  end subroutine  TINVSS
  !---------------------------------------------------------------
  subroutine getknot      &
      ( ndat, kdeg, xdat, & !--- in
        qknt              ) !--- out

  integer, intent(in) :: ndat  !  number of data
  integer, intent(in) :: kdeg  !  degree of Spline + 1
  real(DP), intent(in) :: xdat( ndat ) ! data of independent var.

  real(DP), intent(out) :: qknt( ndat+kdeg ) ! knots for B-Spline

  !--- local
  integer :: i

  do i = 1, kdeg
    qknt( i ) = xdat( 1 )
  end do

  do i = 1, ndat-kdeg
    qknt( i+kdeg ) = ( xdat( i )+xdat( i+kdeg ) )*0.50_DP
  end do

  do i = 1, kdeg
    qknt( ndat+i ) = xdat( ndat )
  end do

  return

  end subroutine getknot
  !---------------------------------------------------------------
  recursive function fbspl ( ndat, inum, kdeg, qknt, xdat, x ) &
  result (bspl)

  real(DP) :: bspl

  integer, intent(in) :: ndat  !  number of data
  integer, intent(in) :: inum  !  index of B-Spline
  integer, intent(in) :: kdeg  !  degree of B-Spline + 1
  real(DP), intent(in) :: qknt( ndat+kdeg )  !  knot of B-Spline
  real(DP), intent(in) :: xdat( ndat ) ! data of independent variable
  real(DP), intent(in) :: x     !  interpolation point

  !--- local
  real(DP) :: bsp1, bsp2

  if ( ( inum == 1 .and. x == xdat( 1 ) ) .or. &
       ( inum == ndat .and. x == xdat( ndat ) ) ) then
    bspl = 1.
    return
  end if

  if ( kdeg == 1 ) then
    if ( x >= qknt( inum ) .and. x < qknt( inum+1 ) ) then
      bspl = 1.0_DP
    else
      bspl = 0.0_DP
    end if
  else
    if ( qknt( inum+kdeg-1 ) /= qknt( inum ) ) then
      bsp1 = ( x-qknt( inum ) ) &
            /( qknt( inum+kdeg-1 )-qknt( inum ) ) &
           * fbspl( ndat, inum, kdeg-1, qknt, xdat, x )
    else
      bsp1 = 0.0_DP
    end if
    if ( qknt( inum+kdeg ) /= qknt( inum+1 ) ) then
      bsp2 = ( qknt( inum+kdeg )-x ) &
            /( qknt( inum+kdeg )-qknt( inum+1 ) ) &
           * fbspl( ndat, inum+1, kdeg-1, qknt, xdat, x )
    else
      bsp2 =  0.0_DP
    end if
    bspl = bsp1 + bsp2
  end if

  end function fbspl
  !---------------------------------------------------------------
  function fpb( ndat, i, kdeg, qknt, xdat, elm, x )

  real :: fpb
  integer :: ndat, i, kdeg
  real(DP) :: qknt( ndat+kdeg ), xdat( ndat ), elm( ndat,ndat )
  real(DP) :: x

  integer :: l
  real(DP) :: sum

  sum = 0.0_DP
  do l = 1, ndat
    sum = sum + elm( l,i )*fbspl( ndat, l, kdeg, qknt, xdat, x )
  end do

  fpb = sum

  end function fpb
  !---------------------------------------------------------------
  subroutine getmatrx           &
      ( ndat, kdeg, qknt, xdat, & !--- in
        elm                     ) !--- out

!  use scale_tinvss, only: TINVSS

  integer, intent(in) :: ndat
  integer, intent(in) :: kdeg
  real(DP), intent(in) :: qknt( ndat+kdeg )
  real(DP), intent(in) :: xdat( ndat )

  real(DP), intent(out) :: elm( ndat,ndat )

  !--- local
  real(DP) :: dt
  integer :: iw( 2*ndat ), i, j, inder
  real(DP), parameter :: eps = 0.

  do i = 1, ndat
  do j = 1, ndat
    elm( i,j ) = fbspl( ndat, j, kdeg, qknt, xdat, xdat( i ) )
  end do
  end do

  call TINVSS( ndat, elm, dt, eps, ndat, iw, inder )

  return

  end subroutine getmatrx
  !---------------------------------------------------------------
  subroutine getcoef           &
      ( ndat, kdeg, elm, ydat, & !--- in
        coef                   ) !--- out

  integer, intent(in) :: ndat  !  number of data
  integer, intent(in) :: kdeg  !  degree of Spline + 1
  real(DP), intent(in) :: elm( ndat,ndat ) ! matrix ( inverse )
  real(DP), intent(in) :: ydat( ndat )  ! data of dependent var.

  real(DP), intent(out) :: coef( ndat ) ! expansion coefficient

  !--- local
  integer :: i, j
  real(DP) :: sum

  do i = 1, ndat
    sum = 0.0_DP
    do j = 1, ndat
      sum = sum + elm( i,j )*ydat( j )
    end do
    coef( i ) = sum
  end do

  return

  end subroutine getcoef
  !---------------------------------------------------------------
  function fspline ( ndat, kdeg, coef, qknt, xdat, x )

  integer, intent(in) :: ndat
  integer, intent(in) :: kdeg
  real(DP), intent(in) :: coef( ndat )
  real(DP), intent(in) :: qknt( ndat+kdeg )
  real(DP), intent(in) :: xdat( ndat )
  real(DP), intent(in) :: x

  real(DP) :: fspline

  !--- local
  real(DP) :: sum
  integer :: i

  sum = 0.0_DP
  do i = 1, ndat
    sum = sum + coef( i )*fbspl( ndat, i, kdeg, qknt, xdat, x )
  end do

  fspline = sum

  return

  end function fspline
  !---------------------------------------------------------------
  subroutine getcoef2                 &
      ( mdat, ndat, kdeg, ldeg,       & !--- in
        xdat, ydat, qknt, rknt, zdat, & !--- in
        coef                          ) !--- out

!  use scale_tinvss, only: TINVSS

  integer, intent(in) :: mdat  !  number of data (x-direction)
  integer, intent(in) :: ndat  !  number of data (y-direction)
  integer, intent(in) :: kdeg  !  degree of Spline + 1 (x)
  integer, intent(in) :: ldeg  !  degree of Spline + 1 (y)
  real(DP), intent(in) :: xdat( mdat )  ! data of independent var. (x)
  real(DP), intent(in) :: ydat( ndat )  ! data of independent var. (y)
  real(DP), intent(in) :: qknt( mdat+kdeg ) ! knots of B-Spline (x)
  real(DP), intent(in) :: rknt( ndat+ldeg ) ! knots of B-Spline (y)
  real(DP), intent(in) :: zdat( mdat,ndat ) ! data of dependent var.

  real(DP), intent(out) :: coef( mdat,ndat ) ! expansion coefficient

  !--- local
  real(DP) :: elmx( mdat,mdat ), elmy( ndat,ndat )
  integer :: iw1( 2*mdat ), iw2( 2*ndat )
  real(DP) :: beta( mdat,ndat ), sum, dt
  real(DP), parameter :: eps = 0.0_DP
  integer :: i, j, k, l, inder

  do i = 1, mdat
  do j = 1, mdat
    elmx( i,j ) = fbspl( mdat, j, kdeg, qknt, xdat, xdat( i ) )
  end do
  end do
  call TINVSS( mdat, elmx, dt, eps, mdat, iw1, inder )

  do l = 1, ndat
  do i = 1, mdat
    sum = 0.0_DP
    do j = 1, mdat
      sum = sum + elmx( i,j )*zdat( j,l )
    end do
    beta( i,l ) = sum
  end do
  end do

  do i = 1, ndat
  do j = 1, ndat
    elmy( i,j ) = fbspl( ndat, j, ldeg, rknt, ydat, ydat( i ) )
  end do
  end do
  call TINVSS( ndat, elmy, dt, eps, ndat, iw2, inder )

  do k = 1, mdat
  do i = 1, ndat
    sum = 0.0_DP
    do j = 1, ndat
      sum = sum + elmy( i,j )*beta( k,j )
    end do
    coef( k,i ) = sum
  end do
  end do

  return

  end subroutine getcoef2
  !---------------------------------------------------------------
  function fspline2                          &
             ( mdat, ndat, kdeg, ldeg,       & !--- in
               coef, qknt, rknt, xdat, ydat, & !--- in
               x, y                          ) !--- in

  integer, intent(in) :: mdat
  integer, intent(in) :: ndat
  integer, intent(in) :: kdeg
  integer, intent(in) :: ldeg
  real(DP), intent(in) :: coef( mdat,ndat )
  real(DP), intent(in) :: qknt( mdat+kdeg )
  real(DP), intent(in) :: rknt( ndat+ldeg )
  real(DP), intent(in) :: xdat( mdat ), ydat( ndat )
  real(DP), intent(in) :: x, y

  real(DP) :: fspline2

  !--- local
  real(DP) :: sum, add
  integer :: i, j

  sum = 0.0_DP
  do i = 1, mdat
  do j = 1, ndat
    add = coef( i,j )*fbspl( mdat, i, kdeg, qknt, xdat, x ) &
                     *fbspl( ndat, j, ldeg, rknt, ydat, y )
    sum = sum + add
  end do
  end do

  fspline2 = sum

  return

  end function fspline2
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_mp_suzuki10
!-------------------------------------------------------------------------------
