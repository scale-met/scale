!-------------------------------------------------------------------------------
!> Module Spectran Bin Microphysics
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
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer_suzuki10
  use scale_const, only: &
     pi  => CONST_PI,  &
     EPS => CONST_EPS, &
     CONST_CPvap,  &
     CONST_CPdry,  &
     CONST_CVdry,  &
     CONST_CL,     &
     CONST_CI,     &
     CONST_DWATR,  &
     CONST_GRAV,   &
     CONST_Rvap,   &
     CONST_Rdry,   &
     CONST_LHV0,   &
     CONST_LHS0,   &
     CONST_EMELT,  &
     CONST_TEM00,  &
     CONST_TMELT,  &
     CONST_PSAT0,  &
     CONST_PRE00,  &
     CONST_EPSvap, &
     CONST_THERMODYN_TYPE
  use scale_const, only: &
     CP    => CONST_CPdry, &
     Rvap  => CONST_Rvap,  &
     ESAT0 => CONST_PSAT0, &
     QLMLT => CONST_EMELT, &
     TMLT  => CONST_TMELT, &
     TEMP0 => CONST_TEM00, &
     RHOW  => CONST_DWATR
  use scale_atmos_saturation, only: &
     SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq, &
     SATURATION_pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice, &
     LovR_liq,     &
     LovR_ice,     &
     CPovR_liq,    &
     CPovR_ice
  use scale_atmos_thermodyn, only: &
     THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres, &
     THERMODYN_pott      => ATMOS_THERMODYN_pott
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
  !--- Indeces for determining species of cloud particle
  integer, parameter :: il = 1               !--- index for liquid  water
  integer, parameter :: ic = 2               !--- index for columner ice
  integer, parameter :: ip = 3               !--- index for plate ice
  integer, parameter :: id = 4               !--- index for dendrite ice
  integer, parameter :: iss= 5               !--- index for snow
  integer, parameter :: ig = 6               !--- index for graupel
  integer, parameter :: ih = 7               !--- index for hail

  !--- bin information of hydrometeors
  real(RP) :: dxmic                          !--- d( log(m) ) of hydrometeor bin
  real(RP), allocatable :: xctr( : )         !--- log( m ) value of bin center [xctr = exp( 4/3*pi*DWATER*radc^3 )]
  real(RP), allocatable :: xbnd( : )         !--- log( m ) value of bin boundary
  real(RP), allocatable :: radc( : )         !--- radius of hydrometeor at bin center [m]
  real(RP), allocatable :: cctr( :,: )       !--- capacitance of hydrometeor at bin center (C of equation A.17 in Suzuki 2004)
  real(RP), allocatable :: cbnd( :,: )       !--- capacitance of hydrometeor at bin boundary (C of equation A.17 in Suzuki 2004)
  real(RP), allocatable :: ck( :,:,:,: )     !-- collection kernel (K of equation A.20 in Suzuki 2004)
  real(RP), allocatable :: vt( :,: )         !--- terminal velocity of hydrometeor [m/s]
  real(RP), allocatable :: br( :,: )         !--- bulk density of hydrometeor [kg/m^3]
  integer,  allocatable  :: ifrsl( :,:,: )     !--- type of species after collision
  !--- bin information of aerosol (not supported)
  real(RP), allocatable :: xactr( : )        !--- log( ma ) value of bin center ( ma is mass of aerosol )
  real(RP), allocatable :: xabnd( : )        !--- log( ma ) value of bin boundary ( ma is mass of aerosol )
  real(RP), allocatable :: rada( : )         !--- radius of aerosol at bin center [m]

  real(RP), allocatable :: expxctr( : )      !--- exp( xctr )
  real(RP), allocatable :: expxbnd( : )      !--- exp( xbnd )
  real(RP), allocatable :: expxactr( : )     !--- exp( xactr )
  real(RP), allocatable :: expxabnd( : )     !--- exp( xabnd )
  real(RP), allocatable :: rexpxctr( : )     !--- 1/exp( xctr )
  real(RP), allocatable :: rexpxbnd( : )     !--- 1/exp( xbnd )
  real(RP), allocatable :: rexpxactr( : )    !--- 1/exp( xactr )
  real(RP), allocatable :: rexpxabnd( : )    !--- 1/exp( xabnd1 )

  real(RP) :: dxaer                          !--- d( log(ma) ) of aerosol bin
  real(RP) :: xasta                          !--- exponential of mass of aerosol for smallest aerosol bin
  real(RP) :: xaend                          !--- exponential of mass of aerosol for largest aerosol bin
  real(RP), allocatable, save :: vterm( :,:,:,: )      ! terminal velocity
  integer, private, save :: MP_NSTEP_SEDIMENTATION    ! number of fractional step for sedimentation
  real(RP), private, save :: MP_RNSTEP_SEDIMENTATION  ! 1/MP_NSTEP_SEDIMENTATION
  real(DP), private, save :: MP_DTSEC_SEDIMENTATION   ! DT for sedimentation
  integer, private, save :: MP_ntmax_sedimentation= 1 ! maxinum fractional step
  real(RP), private :: flg_thermodyn         !--- flg for lhv and lhs (0 -> SIMPLE, 1 -> EXACT )
  real(RP), private :: RTEM00                !--- 1/CONST_TEM00

  !--- constant for bin
  real(RP), parameter :: cldmin = 1.0E-10_RP !--- threshould for cloud is regarded existing
  real(RP), parameter :: OneovThird   = 1.0_RP/3.0_RP
  real(RP), parameter :: ThirdovForth = 3.0_RP/4.0_RP
  real(RP), parameter :: TwoovThird   = 2.0_RP/3.0_RP
  !--- constant for aerosol
  real(RP) :: rhoa   = 2.25E+03_RP           ! density of aerosol ( NaCl )
  real(RP) :: emaer  = 58.0_RP               ! molecular weight of aerosol ( salt )
  real(RP) :: emwtr  = 18.0_RP               ! molecular weight of water
  real(RP) :: rasta  = 1.E-08_RP             ! minimum radius of aerosol (m)
  real(RP) :: raend  = 1.E-06_RP             ! maximum radius of aerosol (m)
  real(RP) :: r0a    = 1.E-07_RP             ! average radius of aerosol (m)
  logical :: flg_regeneration=.false.        ! flag regeneration of aerosol
  logical :: flg_nucl=.false.                ! flag nucleated cloud move into smallest bin
  logical :: flg_icenucl=.false.             ! flag ice nucleation
  logical :: flg_sf_aero =.false.            ! flag surface flux of aerosol
  integer, private, save :: rndm_flgp = 0    ! flag for sthastic integration for coll.-coag.
  logical, private, save :: MP_doautoconversion = .true.  ! apply collision process ?
  logical, private, save :: MP_doprecipitation  = .true.  ! apply sedimentation of hydrometeor ?
  logical, private, save :: MP_donegative_fixer = .true.  ! apply negative fixer?

  real(RP), allocatable :: marate( : )                ! mass rate of each aerosol bin to total aerosol mass
  integer, allocatable, save :: ncld( : )             ! bin number of aerosol in bin of hydrometeor
  integer, private, save       :: K10_1, K10_2        ! scaling factor for 10m value (momentum)
  real(RP), private            :: R10M1, R10M2        ! scaling factor for 10m value (momentum)
  real(RP), private            :: R10H1, R10H2        ! scaling factor for 10m value (heat)
  real(RP), private            :: R10E1, R10E2        ! scaling factor for 10m value (tracer)

  character(11),parameter :: fname_micpara="micpara.dat" !--- file name
  integer(4) :: fid_micpara

  !--- Use for stochastic method
  integer, allocatable :: blrg( :,: ), bsml( :,: )
  real(RP) :: wgtbin
  integer  :: mspc, mbin
  real(RP), private :: rndm(1,1,1)

  !--- use for model without aerosol
  real(RP), private :: c_ccn = 100.E+6_RP    ! N0 of Nc = N0*s^kappa
  real(RP), private :: kappa = 0.462_RP      ! kappa of Nc = N0*s^kappa
  !--- use for aerosol coupled model
  real(RP), private :: sigma = 7.5E-02_RP    ! water surface tension [ N/m2 ] (sigma in eq. (A.11) of Suzuki (2004) )
  real(RP), private :: vhfct = 2.0_RP        ! van't hoff factor (i in eq.(A.11) of Suzuki (2004))

  real(RP), parameter :: tcrit = 271.15_RP
  integer, private, allocatable :: kindx( :,: )

  !--- for creating micpara.dat (mkpara)
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
  !> Setup
  subroutine ATMOS_PHY_MP_suzuki10_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop,    &
       PRC_masterrank, &
       PRC_IsMaster
    use scale_const, only: &
       CONST_DWATR, &
       CONST_DICE
    use scale_comm, only: &
       COMM_world,    &
       COMM_datatype
    use scale_grid, only: &
       CDZ => GRID_CDZ, &
       CZ  => GRID_CZ
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use scale_tracer, only: &
       QAD => QA
    implicit none

    character(len=*), intent(in) :: MP_TYPE

    real(RP) :: RHO_AERO  !--- density of aerosol
    real(RP) :: R0_AERO   !--- center radius of aerosol (um)
    real(RP) :: R_MIN     !--- minimum radius of aerosol (um)
    real(RP) :: R_MAX     !--- maximum radius of aerosol (um)
    real(RP) :: S10_EMAER !--- moleculer weight of aerosol
    logical :: S10_FLAG_REGENE  !--- flag of regeneration
    logical :: S10_FLAG_NUCLEAT !--- flag of regeneration
    logical :: S10_FLAG_ICENUCLEAT !--- flag of regeneration
    logical :: S10_FLAG_SFAERO  !--- flag of surface flux of aeorol
    integer :: S10_RNDM_FLGP  !--- flag of surface flux of aeorol
    integer :: S10_RNDM_MSPC
    integer :: S10_RNDM_MBIN

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_ntmax_sedimentation, &
       MP_doautoconversion, &
       MP_doprecipitation, &
       MP_donegative_fixer

    NAMELIST / PARAM_ATMOS_PHY_MP_SUZUKI10 / &
       RHO_AERO,  &
       R_MIN, &
       R_MAX, &
       R0_AERO,   &
       S10_EMAER, &
       S10_FLAG_REGENE,  &
       S10_FLAG_NUCLEAT, &
       S10_FLAG_ICENUCLEAT, &
       S10_FLAG_SFAERO,  &
       S10_RNDM_FLGP, &
       S10_RNDM_MSPC, &
       S10_RNDM_MBIN, &
       c_ccn, kappa, &
       sigma, vhfct

    real(RP), parameter :: max_term_vel = 10.0_RP !-- terminal velocity for calculate dt of sedimentation
    integer :: nstep_max
    integer :: nnspc, nnbin
    integer :: nn, mm, mmyu, nnyu
    integer :: myu, nyu, i, j, k, n, ierr
    !---------------------------------------------------------------------------

    !--- allocation
    allocate( xctr( nbin ) )
    allocate( xbnd( nbin+1 ) )
    allocate( radc( nbin ) )
    allocate( cctr( nbin,nspc_mk ) )
    allocate( cbnd( nbin+1,nspc_mk ) )
    allocate( ck( nspc_mk,nspc_mk,nbin,nbin ) )
    allocate( vt( nspc_mk,nbin ) )
    allocate( br( nspc_mk,nbin ) )
    allocate( ifrsl( 2,nspc_mk,nspc_mk ) )
    allocate( expxctr( nbin ) )
    allocate( expxbnd( nbin+1 ) )
    allocate( rexpxctr( nbin ) )
    allocate( rexpxbnd( nbin+1 ) )
    if ( nccn /= 0 ) then
      allocate( xactr( nccn ) )
      allocate( xabnd( nccn+1 ) )
      allocate( rada( nccn ) )
      allocate( expxactr( nccn ) )
      allocate( expxabnd( nccn+1 ) )
      allocate( rexpxactr( nccn ) )
      allocate( rexpxabnd( nccn+1 ) )
    endif

    mbin = nbin/2
    mspc = nspc_mk*nspc_mk

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for SBM (warm cloud)'

    if ( MP_TYPE /= 'SUZUKI10' ) then
       write(*,*) 'xxx ATMOS_PHY_MP_TYPE is not SUZUKI10. Check!'
       call PRC_MPIstop
    endif

    RHO_AERO = rhoa
    S10_EMAER = emaer
    R_MIN = rasta
    R_MAX = raend
    R0_AERO = r0a
    S10_FLAG_REGENE = flg_regeneration
    S10_FLAG_NUCLEAT = flg_nucl
    S10_FLAG_ICENUCLEAT = flg_icenucl
    S10_FLAG_SFAERO = flg_sf_aero
    S10_RNDM_FLGP = rndm_flgp
    S10_RNDM_MSPC = mspc
    S10_RNDM_MBIN = mbin

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if ( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP, Check!'
     call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SUZUKI10,iostat=ierr)

    if ( ierr < 0 ) then !--- missing
     if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
     write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SUZUKI10, Check!'
     call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP_SUZUKI10)

    if ( nspc /= 1 .AND. nspc /= 7 ) then
       write(*,*) 'xxx nspc should be set as 1(warm rain) or 7(mixed phase) check!'
       call PRC_MPIstop
    endif

    rhoa = RHO_AERO
    emaer = S10_EMAER
    rasta = R_MIN
    raend = R_MAX
    r0a   = R0_AERO
    flg_regeneration = S10_FLAG_REGENE
    flg_nucl = S10_FLAG_NUCLEAT
    flg_icenucl = S10_FLAG_ICENUCLEAT
    flg_sf_aero = S10_FLAG_SFAERO
    rndm_flgp = S10_RNDM_FLGP
    mspc = S10_RNDM_MSPC
    mbin = S10_RNDM_MBIN

    !--- read micpara.dat (microphysical parameter) and broad cast
    if ( PRC_IsMaster ) then

      fid_micpara = IO_get_available_fid()
      !--- open parameter of cloud microphysics
      open ( fid_micpara, file = fname_micpara, form = 'formatted', status = 'old', iostat=ierr )

      !--- micpara.dat does not exist
      if ( ierr == 0 ) then

        read( fid_micpara,* ) nnspc, nnbin

        if ( nnbin /= nbin ) then
           write(*,*) 'xxx nbin in inc_tracer and nbin in micpara.dat is different check!'
           call PRC_MPIstop
        endif

        ! grid parameter
        if( IO_L ) write(IO_FID_LOG,*)  '*** Radius of cloud ****'
        do n = 1, nbin
          read( fid_micpara,* ) nn, xctr( n ), radc( n )
          if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)') &
                    "Radius of ", n, "th cloud bin (bin center)= ", radc( n ) , "[m]"
        enddo
        do n = 1, nbin+1
          read( fid_micpara,* ) nn, xbnd( n )
        enddo
        read( fid_micpara,* ) dxmic
        if( IO_L ) write(IO_FID_LOG,*)  '*** Width of Cloud SDF= ', dxmic

        ! capacity
        do myu = 1, nspc_mk
         do n = 1, nbin
!          read( fid_micpara,* ) mmyu, nn, cctr( myu,n )
          read( fid_micpara,* ) mmyu, nn, cctr( n,myu )
         enddo
         do n = 1, nbin+1
!          read( fid_micpara,* ) mmyu, nn, cbnd( myu,n )
          read( fid_micpara,* ) mmyu, nn, cbnd( n,myu )
         enddo
        enddo

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

        if ( nnbin /= nbin ) then
           write(*,*) 'xxx nbin in inc_tracer and nbin in micpara.dat is different check!'
           call PRC_MPIstop
        endif

        ! grid parameter
        if( IO_L ) write(IO_FID_LOG,*)  '*** Radius of cloud ****'
        do n = 1, nbin
          read( fid_micpara,* ) nn, xctr( n ), radc( n )
          if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)') &
                    "Radius of ", n, "th cloud bin (bin center)= ", radc( n ) , "[m]"
        enddo
        do n = 1, nbin+1
          read( fid_micpara,* ) nn, xbnd( n )
        enddo
        read( fid_micpara,* ) dxmic
        if( IO_L ) write(IO_FID_LOG,*)  '*** Width of Cloud SDF= ', dxmic

        ! capacity
        do myu = 1, nspc_mk
         do n = 1, nbin
!          read( fid_micpara,* ) mmyu, nn, cctr( myu,n )
          read( fid_micpara,* ) mmyu, nn, cctr( n,myu )
         enddo
         do n = 1, nbin+1
!          read( fid_micpara,* ) mmyu, nn, cbnd( myu,n )
          read( fid_micpara,* ) mmyu, nn, cbnd( n,myu )
         enddo
        enddo

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

    call MPI_BCAST( xctr, nbin,                      COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( dxmic,1,                         COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( xbnd, nbin+1,                    COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( cctr, nbin*nspc_mk,              COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( cbnd, (nbin+1)*nspc_mk,          COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( ck,   nspc_mk*nspc_mk*nbin*nbin, COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( br, nbin*nspc_mk,                COMM_datatype, PRC_masterrank, COMM_world, ierr )
    call MPI_BCAST( vt, nbin*nspc_mk,                COMM_datatype, PRC_masterrank, COMM_world, ierr )

    !--- aerosol ( CCN ) (not supported)
    if ( nccn /= 0 ) then

    allocate ( ncld( 1:nccn ) )
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

    if ( flg_sf_aero ) then
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
    if ( rndm_flgp > 0 ) then
     call random_setup( IA*JA*KA )
    endif

    if ( nccn /= 0 ) then
     do n = 1, nccn
      expxactr( n ) = exp( xactr( n ) )
      rexpxactr( n ) = 1.0_RP / exp( xactr( n ) )
     enddo
     do n = 1, nccn+1
      expxabnd( n ) = exp( xabnd( n ) )
      rexpxabnd( n ) = 1.0_RP / exp( xabnd( n ) )
     enddo
    endif

    allocate( vterm(KA,IA,JA,QAD) )
    vterm(:,:,:,:) = 0.0_RP
    do myu = 1, nspc
    do n = 1, nbin
      vterm(:,:,:,I_QV+(myu-1)*nbin+n) = -vt( myu,n )
    enddo
    enddo
    do n = 1, nbin
      expxctr( n ) = exp( xctr( n ) )
      rexpxctr( n ) = 1.0_RP / exp( xctr( n ) )
    enddo
    do n = 1, nbin+1
      expxbnd( n ) = exp( xbnd( n ) )
      rexpxbnd( n ) = 1.0_RP / exp( xbnd( n ) )
    enddo

    allocate( kindx(nbin,nbin) )
    call getrule( ifrsl,kindx )

    if ( CONST_THERMODYN_TYPE == 'EXACT' ) then
      flg_thermodyn = 1.0_RP
    elseif( CONST_THERMODYN_TYPE == 'SIMPLE' ) then
      flg_thermodyn = 0.0_RP
    endif
    RTEM00 = 1.0_RP / CONST_TEM00

    nstep_max = int ( ( TIME_DTSEC_ATMOS_PHY_MP * max_term_vel ) / minval( CDZ ) )
    MP_ntmax_sedimentation = max( MP_ntmax_sedimentation, nstep_max )

    MP_NSTEP_SEDIMENTATION  = MP_ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(MP_ntmax_sedimentation,kind=RP)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

    if ( IO_L ) write(IO_FID_LOG,*)
    if ( IO_L ) write(IO_FID_LOG,*) '*** Timestep of sedimentation is divided into : ', MP_ntmax_sedimentation, ' step'
    if ( IO_L ) write(IO_FID_LOG,*) '*** DT of sedimentation is : ', MP_DTSEC_SEDIMENTATION, '[s]'

    return

  end subroutine ATMOS_PHY_MP_suzuki10_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  subroutine ATMOS_PHY_MP_suzuki10( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       SFLX_rain, &
       SFLX_snow  )
    use scale_grid_index
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_tracer, only: &
       QAD => QA
    use scale_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer => ATMOS_PHY_MP_negative_fixer, &
       MP_precipitation  => ATMOS_PHY_MP_precipitation
    use scale_comm, only: &
       COMM_barrier
!    use scale_grid, only: &
!       GRID_CDZ
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)
    real(RP), intent(in)    :: CCN (KA,IA,JA)
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)

    real(RP) :: RHOE (KA,IA,JA)
    real(RP) :: POTT (KA,IA,JA)
    real(RP) :: TEMP (KA,IA,JA)
    real(RP) :: PRES (KA,IA,JA)
    real(RP) :: qsat (KA,IA,JA)
    real(RP) :: ssliq(KA,IA,JA)
    real(RP) :: ssice(KA,IA,JA)

    integer  :: ijk_index (KIJMAX,3)
    integer  :: index_cld (KIJMAX)
    integer  :: index_cold(KIJMAX)
    integer  :: index_warm(KIJMAX)
    integer  :: ijkcount, ijkcount_cold, ijkcount_warm
    integer  :: ijk, indirect

    real(RP) :: DENS_ijk(KIJMAX)
    real(RP) :: PRES_ijk(KIJMAX)
    real(RP) :: TEMP_ijk(KIJMAX)
    real(RP) :: Qvap_ijk(KIJMAX)
    real(RP) :: Ghyd_ijk(nbin,nspc,KIJMAX)
    real(RP) :: Gaer_ijk(nccn1    ,KIJMAX)
    real(RP) :: cldsum
    integer  :: countbin

!    logical, save :: ofirst_sdfa = .true.
!    real(RP) :: VELX (IA,JA)
!    real(RP) :: VELY (IA,JA)
!    real(RP) :: SFLX_AERO(IA,JA,nccn)
!    real(RP) :: Uabs, bparam
!    real(RP) :: AMR(KA,IA,JA)

    real(RP) :: FLX_rain  (KA,IA,JA)
    real(RP) :: FLX_snow  (KA,IA,JA)
    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    integer  :: step
    integer  :: k, i, j, m, n, iq
    !---------------------------------------------------------------------------

    if    ( nspc == 1 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cloud microphysics(SBM Liquid water only)'
    elseif( nspc >  1 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cloud microphysics(SBM Mixed phase)'
    endif

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    call PROF_rapstart('MP_ijkconvert', 1)

    ijk = 0
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ijk = ijk + 1
       ijk_index(ijk,1) = i
       ijk_index(ijk,2) = j
       ijk_index(ijk,3) = k
    enddo
    enddo
    enddo

    call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                              PRES(:,:,:),  & ! [OUT]
                              DENS(:,:,:),  & ! [IN]
                              RHOT(:,:,:),  & ! [IN]
                              QTRC(:,:,:,:) ) ! [IN]

    call SATURATION_pres2qsat_liq( qsat(:,:,:), & ! [OUT]
                                   TEMP(:,:,:), & ! [IN]
                                   PRES(:,:,:)  ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ssliq(k,i,j) = QTRC(k,i,j,I_QV) / qsat(k,i,j) - 1.0_RP
    enddo
    enddo
    enddo

    if ( nspc == 1 ) then
       ssice(:,:,:) = 0.0_RP
    else
       call SATURATION_pres2qsat_ice( qsat(:,:,:), & ! [OUT]
                                      TEMP(:,:,:), & ! [IN]
                                      PRES(:,:,:)  ) ! [IN]

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ssice(k,i,j) = QTRC(k,i,j,I_QV) / qsat(k,i,j) - 1.0_RP
       enddo
       enddo
       enddo
    endif

!--- store initial SDF of aerosol
!--- this option is not supported
!    if ( ofirst_sdfa ) then
!      allocate( marate( nccn ) )
!      do j = JS, JE
!      do i = IS, IE
!         do k = KS, KE
!           sum2 = 0.0_RP
!           do n = 1, nccn
!             marate( n ) = gdga(k,i,j,n)*rexpxactr( n )
!             sum2 = sum2 + gdga(k,i,j,n)*rexpxactr( n )
!           enddo
!         enddo
!      enddo
!      enddo
!      if ( sum2 /= 0.0_RP ) then
!        marate( 1:nccn ) = marate( 1:nccn )/sum2
!        ofirst_sdfa = .false.
!      endif
!    endif

    !--- Arrange array for microphysics
    ijkcount = 0
    ijkcount_cold = 0
    ijkcount_warm = 0

    ijk = 0
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ijk = ijk + 1

       cldsum   = 0.0_RP
       countbin = QQS
       do m = 1, nspc
       do n = 1, nbin
         countbin = countbin + 1
         cldsum   = cldsum + QTRC(k,i,j,countbin) * DENS(k,i,j) / dxmic
       enddo
       enddo

       if (      cldsum > cldmin       &
            .OR. ssliq(k,i,j) > 0.0_RP &
            .OR. ssice(k,i,j) > 0.0_RP ) then

          ijkcount = ijkcount + 1

          index_cld(ijkcount) = ijk

          DENS_ijk(ijkcount) = DENS(k,i,j)
          PRES_ijk(ijkcount) = PRES(k,i,j)
          TEMP_ijk(ijkcount) = TEMP(k,i,j)
          Qvap_ijk(ijkcount) = QTRC(k,i,j,I_QV)

          countbin = QQS
          do m = 1, nspc
          do n = 1, nbin
             countbin = countbin + 1
             Ghyd_ijk(n,m,ijkcount) = QTRC(k,i,j,countbin) * DENS(k,i,j) / dxmic
          enddo
          enddo

          do n = 1, nccn
             countbin = countbin + 1
             Gaer_ijk(n,ijkcount)   = QTRC(k,i,j,countbin) * DENS(k,i,j) / dxaer
          enddo

          if ( TEMP(k,i,j) < CONST_TEM00 ) then ! cold
            ijkcount_cold = ijkcount_cold + 1
            index_cold(ijkcount_cold) = ijkcount
          else ! warm
            ijkcount_warm = ijkcount_warm + 1
            index_warm(ijkcount_warm) = ijkcount
          endif
       endif

    enddo
    enddo
    enddo

    call PROF_rapend  ('MP_ijkconvert', 1)

    ! tentative timername registration
    call PROF_rapstart('MP_SBM_Main', 1)
    call PROF_rapend  ('MP_SBM_Main', 1)
    call PROF_rapstart('_SBM_Nucleat',    2)
    call PROF_rapend  ('_SBM_Nucleat',    2)
!    call PROF_rapstart('_SBM_NucleatA',   2)
!    call PROF_rapend  ('_SBM_NucleatA',   2)
    call PROF_rapstart('_SBM_Liqphase',   2)
    call PROF_rapend  ('_SBM_Liqphase',   2)
    call PROF_rapstart('_SBM_Icephase',   2)
    call PROF_rapend  ('_SBM_Icephase',   2)
    call PROF_rapstart('_SBM_Mixphase',   2)
    call PROF_rapend  ('_SBM_Mixphase',   2)
    call PROF_rapstart('_SBM_AdvLiq',     3)
    call PROF_rapend  ('_SBM_AdvLiq',     3)
    call PROF_rapstart('_SBM_AdvIce',     3)
    call PROF_rapend  ('_SBM_AdvIce',     3)
    call PROF_rapstart('_SBM_AdvMix',     3)
    call PROF_rapend  ('_SBM_AdvMix',     3)
!    call PROF_rapstart('_SBM_FAero',      2)
!    call PROF_rapend  ('_SBM_FAero',      2)
    call PROF_rapstart('_SBM_Freezing',   2)
    call PROF_rapend  ('_SBM_Freezing',   2)
    call PROF_rapstart('_SBM_IceNucleat', 2)
    call PROF_rapend  ('_SBM_IceNucleat', 2)
    call PROF_rapstart('_SBM_Melting',    2)
    call PROF_rapend  ('_SBM_Melting',    2)
    call PROF_rapstart('_SBM_CollCoag',   2)
    call PROF_rapend  ('_SBM_CollCoag',   2)
!    call PROF_rapstart('_SBM_CollCoagR',  2)
!    call PROF_rapend  ('_SBM_CollCoagR',  2)

    if ( ijkcount > 0 ) then

    call PROF_rapstart('MP_SBM_Main', 1)

    call MP_Suzuki10( ijkcount,                   & ! [IN]
                      ijkcount_cold,              & ! [IN]
                      ijkcount_warm,              & ! [IN]
                      index_cold(    1:ijkcount), & ! [IN]
                      index_warm(    1:ijkcount), & ! [IN]
                      DENS_ijk  (    1:ijkcount), & ! [IN]
                      PRES_ijk  (    1:ijkcount), & ! [IN]
                      TEMP_ijk  (    1:ijkcount), & ! [INOUT]
                      Qvap_ijk  (    1:ijkcount), & ! [INOUT]
                      Ghyd_ijk  (:,:,1:ijkcount), & ! [INOUT]
                      Gaer_ijk  (:,  1:ijkcount), & ! [INOUT]
                      dt                          ) ! [IN]

    call PROF_rapend  ('MP_SBM_Main', 1)

!    if ( flg_sf_aero ) then
!     do j = JS-2, JE+2
!     do i = IS-2, IE+1
!       VELX(i,j) = MOMX(K10_1,i,j) / ( DENS(K10_1,i+1,j)+DENS(K10_1,i,j) ) * R10M1 &
!                 + MOMX(K10_2,i,j) / ( DENS(K10_2,i+1,j)+DENS(K10_2,i,j) ) * R10M2
!     enddo
!     enddo
!
!     do j = JS-2, JE+1
!     do i = IS-2, IE+2
!       VELY(i,j) = MOMY(K10_1,i,j) / ( DENS(K10_1,i,j+1)+DENS(K10_1,i,j) ) * R10M1 &
!                 + MOMY(K10_2,i,j) / ( DENS(K10_2,i,j+1)+DENS(K10_2,i,j) ) * R10M2
!     enddo
!     enddo
!    endif
!
!    !--- SURFACE FLUX by Monahan et al. (1986)
!    if ( flg_sf_aero .AND. nccn /= 0 ) then
!     do j = JS, JE
!     do i = IS, IE
!          ijk = ( j - JS ) * KMAX * IMAX &
!              + ( i - IS ) * KMAX        &
!              + ( KS - KS )              &
!              + 1
!       Uabs = sqrt(  ( ( VELX(i,j) + VELX(i-1,j  ) ) * 0.50_RP )**2 &
!                   + ( ( VELY(i,j) + VELY(i  ,j-1) ) * 0.50_RP )**2 )
!       do n = 1, nccn
!        if ( rada( n ) <= 2.0E-5_RP .AND. rada( n ) >= 3.0E-7_RP ) then
!         bparam = ( 0.38_RP - log( rada( n ) ) )/0.65_RP
!         SFLX_AERO(i,j,n) = 1.373_RP * Uabs**( 3.41_RP ) * rada( n )**( -3.0_RP ) &
!                          * ( 1.0_RP + 0.057_RP * rada( n )**( 1.05_RP ) ) &
!                          * 10.0_RP**( 1.19_RP * exp( -bparam*bparam ) )
!         ! convert from [#/m^2/um/s] -> [kg/m^3/unit log (m)]
!         SFLX_AERO(i,j,n) = SFLX_AERO(i,j,n) / DENS(KS,i,j) &
!                          / GRID_CDZ(KS) * rada( n ) / 3.0_RP * dt * expxactr( n )
!         Gaer_ijk(n,ijk) = Gaer_ijk(n,ijk) + SFLX_AERO(i,j,n)/dxaer
!        endif
!       enddo
!     enddo
!     enddo
!    endif

    call PROF_rapstart('MP_ijkconvert', 1)

    !---- return original array
    do ijk = 1, ijkcount
       indirect = index_cld(ijk)
       i = ijk_index(indirect,1)
       j = ijk_index(indirect,2)
       k = ijk_index(indirect,3)

       TEMP(k,i,j)      = TEMP_ijk(ijk)
       QTRC(k,i,j,I_QV) = Qvap_ijk(ijk)

       countbin = QQS
       do m = 1, nspc
       do n = 1, nbin
          countbin = countbin + 1
          QTRC(k,i,j,countbin) = Ghyd_ijk(n,m,ijk) / DENS(k,i,j) * dxmic
       enddo
       enddo

       do n = 1, nccn
          countbin = countbin + 1
          QTRC(k,i,j,countbin) = Gaer_ijk(n,ijk)   / DENS(k,i,j) * dxaer
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       countbin = QQS
       do m = 1, nspc
       do n = 1, nbin
          countbin = countbin + 1
          if ( QTRC(k,i,j,countbin) < EPS ) then
             QTRC(k,i,j,countbin) = 0.0_RP
          endif
       enddo
       enddo
    enddo
    enddo
    enddo

    call THERMODYN_pott( POTT(:,:,:),  & ! [OUT]
                         TEMP(:,:,:),  & ! [IN]
                         PRES(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT(k,i,j) = POTT(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo

!    if ( nccn /= 0 ) then
!       AMR(:,:,:) = 0.0_RP
!       do j = JS, JE
!       do i = IS, IE
!       do k = KS, KE
!          do n = 1, nccn
!             AMR(k,i,j) = AMR(k,i,j) + QTRC(k,i,j,QQE-1+n)
!          enddo
!       enddo
!       enddo
!       enddo
!    endif

    call PROF_rapend  ('MP_ijkconvert', 1)

    endif

    call PROF_rapstart('MP_barrier', 1)
    call COMM_barrier
    call PROF_rapend  ('MP_barrier', 1)

    FLX_rain(:,:,:) = 0.0_RP
    FLX_snow(:,:,:) = 0.0_RP

    if ( MP_doprecipitation ) then

       call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                            RHOT(:,:,:),  & ! [IN]
                            QTRC(:,:,:,:) ) ! [IN]

       do step = 1, MP_NSTEP_SEDIMENTATION

          call THERMODYN_temp_pres_E( temp(:,:,:),  & ! [OUT]
                                      pres(:,:,:),  & ! [OUT]
                                      DENS(:,:,:),  & ! [IN]
                                      RHOE(:,:,:),  & ! [IN]
                                      QTRC(:,:,:,:) ) ! [IN]

          call MP_precipitation( wflux_rain(:,:,:),     & ! [OUT]
                                 wflux_snow(:,:,:),     & ! [OUT]
                                 DENS    (:,:,:),       & ! [INOUT]
                                 MOMZ    (:,:,:),       & ! [INOUT]
                                 MOMX    (:,:,:),       & ! [INOUT]
                                 MOMY    (:,:,:),       & ! [INOUT]
                                 RHOE    (:,:,:),       & ! [INOUT]
                                 QTRC    (:,:,:,:),     & ! [INOUT]
                                 vterm   (:,:,:,:),     & ! [IN]
                                 temp    (:,:,:),       & ! [IN]
                                 MP_DTSEC_SEDIMENTATION ) ! [IN]

          do j = JS, JE
          do i = IS, IE
          do k = KS-1, KE
             FLX_rain(k,i,j) = FLX_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
             FLX_snow(k,i,j) = FLX_snow(k,i,j) + wflux_snow(k,i,j) * MP_RNSTEP_SEDIMENTATION
          enddo
          enddo
          enddo

       enddo

       call THERMODYN_rhot( RHOT(:,:,:),  & ! [OUT]
                            RHOE(:,:,:),  & ! [IN]
                            QTRC(:,:,:,:) ) ! [IN]
    endif

    !##### END MP Main #####

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    SFLX_rain(:,:) = FLX_rain(KS-1,:,:)
    SFLX_snow(:,:) = FLX_snow(KS-1,:,:)

    return
  end subroutine ATMOS_PHY_MP_suzuki10

  !-----------------------------------------------------------------------------
  subroutine MP_Suzuki10( &
       ijkmax,      &
       ijkmax_cold, &
       ijkmax_warm, &
       index_cold,  &
       index_warm,  &
       dens,        &
       pres,        &
       temp,        &
       qvap,        &
       ghyd,        &
       gaer,        &
       dt           )
    implicit none

    integer,  intent(in)    :: ijkmax
    integer,  intent(in)    :: ijkmax_cold
    integer,  intent(in)    :: ijkmax_warm
    integer , intent(in)    :: index_cold(ijkmax)
    integer , intent(in)    :: index_warm(ijkmax)
    real(RP), intent(in)    :: dens      (ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres      (ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp      (ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap      (ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: ghyd      (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(RP), intent(inout) :: gaer      (nccn1,    ijkmax) ! Mass size distribution function of aerosol
    real(DP), intent(in)    :: dt                           ! Time step interval
    !---------------------------------------------------------------------------

    if ( nccn /= 0 ) then
       if ( nspc == 1 ) then
          !---< warm rain only with aerosol tracer >---

          ! nucleation from aerosol
          call nucleata( ijkmax,      & ! [IN]
                         dens(:),     & ! [IN]
                         pres(:),     & ! [IN]
                         temp(:),     & ! [INOUT]
                         qvap(:),     & ! [INOUT]
                         ghyd(:,:,:), & ! [INOUT]
                         gaer(:,:),   & ! [INOUT]
                         dt           ) ! [IN]

          ! condensation / evaporation
          call cndevpsbla( ijkmax,      & ! [IN]
                           dens(:),     & ! [IN]
                           pres(:),     & ! [IN]
                           temp(:),     & ! [INOUT]
                           qvap(:),     & ! [INOUT]
                           ghyd(:,:,:), & ! [INOUT]
                           gaer(:,:),   & ! [INOUT]
                           dt           ) ! [IN]

          if ( MP_doautoconversion ) then
             ! collision-coagulation
             call collmain( ijkmax,      & ! [IN]
                            temp(:),     & ! [IN]
                            ghyd(:,:,:), & ! [INOUT]
                            dt           ) ! [IN]
          endif

       elseif( nspc > 1 ) then
          !---< mixed phase rain only with aerosol tracer >---

          ! nucleation from aerosol
          call nucleata( ijkmax,      & ! [IN]
                         dens(:),     & ! [IN]
                         pres(:),     & ! [IN]
                         temp(:),     & ! [INOUT]
                         qvap(:),     & ! [INOUT]
                         ghyd(:,:,:), & ! [INOUT]
                         gaer(:,:),   & ! [INOUT]
                         dt           ) ! [IN]

          ! freezing / melting
          call freezing( ijkmax,        & ! [IN]
                         ijkmax_cold,   & ! [IN]
                         index_cold(:), & ! [IN]
                         dens(:),       & ! [IN]
                         temp(:),       & ! [INOUT]
                         ghyd(:,:,:),   & ! [INOUT]
                         dt             ) ! [IN]

!          call ice_nucleat( ijkmax,        & ! [IN]
!                            ijkmax_cold,   & ! [IN]
!                            index_cold(:), & ! [IN]
!                            dens(:),       & ! [IN]
!                            pres(:),       & ! [IN]
!                            temp(:),       & ! [INOUT]
!                            qvap(:),       & ! [INOUT]
!                            ghyd(:,:,:),   & ! [INOUT]
!                            dt             ) ! [IN]

          call melting( ijkmax,        & ! [IN]
                        ijkmax_warm,   & ! [IN]
                        index_warm(:), & ! [IN]
                        dens(:),       & ! [IN]
                        temp(:),       & ! [INOUT]
                        ghyd(:,:,:),   & ! [INOUT]
                        dt             ) ! [IN]

          ! condensation / evaporation
          call cndevpsbla( ijkmax,      & ! [IN]
                           dens(:),     & ! [IN]
                           pres(:),     & ! [IN]
                           temp(:),     & ! [INOUT]
                           qvap(:),     & ! [INOUT]
                           ghyd(:,:,:), & ! [INOUT]
                           gaer(:,:),   & ! [INOUT]
                           dt           ) ! [IN]

          if ( MP_doautoconversion ) then
             ! collision-coagulation
             call collmainf( ijkmax,      & ! [IN]
                             temp(:),     & ! [IN]
                             ghyd(:,:,:), & ! [INOUT]
                             dt           ) ! [IN]
          endif

       endif

    elseif( nccn == 0 ) then

       if ( nspc == 1 ) then
          !---< warm rain only without aerosol tracer >---

          ! nucleation
          call nucleat( ijkmax,      & ! [IN]
                        dens(:),     & ! [IN]
                        pres(:),     & ! [IN]
                        temp(:),     & ! [INOUT]
                        qvap(:),     & ! [INOUT]
                        ghyd(:,:,:), & ! [INOUT]
                        dt           ) ! [IN]

          ! condensation / evaporation
          call cndevpsbl( ijkmax,      & ! [IN]
                          dens(:),     & ! [IN]
                          pres(:),     & ! [IN]
                          temp(:),     & ! [INOUT]
                          qvap(:),     & ! [INOUT]
                          ghyd(:,:,:), & ! [INOUT]
                          dt           ) ! [IN]

          if ( MP_doautoconversion ) then
             ! collision-coagulation
             call collmain( ijkmax,      & ! [IN]
                            temp(:),     & ! [IN]
                            ghyd(:,:,:), & ! [INOUT]
                            dt           ) ! [IN]
          endif

       elseif( nspc > 1 ) then
          !---< mixed phase rain only without aerosol tracer >---

          ! nucleation
          call nucleat( ijkmax,      & ! [IN]
                        dens(:),     & ! [IN]
                        pres(:),     & ! [IN]
                        temp(:),     & ! [INOUT]
                        qvap(:),     & ! [INOUT]
                        ghyd(:,:,:), & ! [INOUT]
                        dt           ) ! [IN]

          ! freezing / melting
          call freezing( ijkmax,        & ! [IN]
                         ijkmax_cold,   & ! [IN]
                         index_cold(:), & ! [IN]
                         dens(:),       & ! [IN]
                         temp(:),       & ! [INOUT]
                         ghyd(:,:,:),   & ! [INOUT]
                         dt             ) ! [IN]

          call ice_nucleat( ijkmax,        & ! [IN]
                            ijkmax_cold,   & ! [IN]
                            index_cold(:), & ! [IN]
                            dens(:),       & ! [IN]
                            pres(:),       & ! [IN]
                            temp(:),       & ! [INOUT]
                            qvap(:),       & ! [INOUT]
                            ghyd(:,:,:),   & ! [INOUT]
                            dt             ) ! [IN]

          call melting( ijkmax,        & ! [IN]
                        ijkmax_warm,   & ! [IN]
                        index_warm(:), & ! [IN]
                        dens(:),       & ! [IN]
                        temp(:),       & ! [INOUT]
                        ghyd(:,:,:),   & ! [INOUT]
                        dt             ) ! [IN]

          ! condensation / evaporation
          call cndevpsbl( ijkmax,      & ! [IN]
                          dens(:),     & ! [IN]
                          pres(:),     & ! [IN]
                          temp(:),     & ! [INOUT]
                          qvap(:),     & ! [INOUT]
                          ghyd(:,:,:), & ! [INOUT]
                          dt           ) ! [IN]

          if ( MP_doautoconversion ) then
             ! collision-coagulation
             call collmainf( ijkmax,      & ! [IN]
                             temp(:),     & ! [IN]
                             ghyd(:,:,:), & ! [INOUT]
                             dt           ) ! [IN]
          endif

       endif

    endif

    return
  end subroutine MP_Suzuki10

  !-----------------------------------------------------------------------------
  subroutine nucleat( &
       ijkmax, &
       dens,   &
       pres,   &
       temp,   &
       qvap,   &
       gc,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

  real(RP) :: ssliq(ijkmax)
  real(RP) :: qlevp(ijkmax)              ! LH
  real(RP) :: dmp
  integer :: n
  !
  real(RP) :: n_c
  real(RP) :: sumnum(ijkmax)
  real(RP) :: gcn( nbin,ijkmax )        ! number of cloud particles in each bin (=gc/exp(xctr))
  real(RP) :: psat
  integer  :: ijk

    call PROF_rapstart('_SBM_Nucleat', 2)

  do ijk = 1, ijkmax

    !--- lhv
    qlevp(ijk) = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
    !--- supersaturation
    psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
          * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
    ssliq(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
    ssliq(ijk) = qvap(ijk)/ssliq(ijk)-1.0_RP

  enddo

  sumnum(:) = 0.0_RP
  do ijk = 1, ijkmax
!    if ( ssliq <= 0.0_RP ) cycle
    if ( ssliq(ijk) > 0.0_RP ) then

     !--- use for aerosol coupled model
     !--- mass -> number
     do n = 1, nbin
       gcn( n,ijk ) = gc( n,il,ijk )*rexpxctr( n )
     enddo

     do n = 1, nbin
       sumnum(ijk) = sumnum(ijk) + gcn( n,ijk )*dxmic
     enddo
     n_c = c_ccn * ( ssliq(ijk) * 1.E+2_RP )**( kappa )
     if ( n_c > sumnum(ijk) ) then
       dmp = ( n_c - sumnum(ijk) ) * expxctr( 1 )
       dmp = min( dmp,qvap(ijk)*dens(ijk) )
       gc( 1,il,ijk ) = gc( 1,il,ijk ) + dmp/dxmic
       qvap(ijk) = qvap(ijk) - dmp/dens(ijk)
       qvap(ijk) = max( qvap(ijk),0.0_RP )
       temp(ijk) = temp(ijk) + dmp/dens(ijk)*qlevp(ijk)/cp
     endif
    endif

  enddo

    call PROF_rapend  ('_SBM_Nucleat', 2)

  return
  end subroutine nucleat

  !-----------------------------------------------------------------------------
  subroutine nucleata( &
       ijkmax, &
       dens,   &
       pres,   &
       temp,   &
       qvap,   &
       gc,     &
       ga,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(RP), intent(inout) :: ga  (nccn     ,ijkmax) ! Mass size distribution function of aerosol
    real(DP), intent(in)    :: dtime                  ! Time step interval

  real(RP) :: gan( nccn )           ! size distribution function ( aerosol ) : number ( gan = ga/exp( xactr ) )
  real(RP) :: ssliq, ssice, qlevp   ! supersaturatioin of liq. and ice, and LH
  real(RP) :: acoef, bcoef          ! A and B in eq. (A.11) of Suzuki (2004)
  real(RP) :: rcrit                 ! critical radius (rcrit, r_N,crit of (A.11) of Suzuki (2004))
  real(RP) :: xcrit                 ! exp of hydrometeror whose radi is corresponding to rcrit (xcrit)
  real(RP) :: ractr, rcld, xcld, part, dmp
  integer :: n, nc, ncrit
!  integer, allocatable, save :: ncld( : )
!  integer, save :: ncld( 1:nccn )
!  logical, save :: ofirst(1:ijkmax) = .true.
  !
  real(RP) :: psat
  integer  :: ijk

    call PROF_rapstart('_SBM_NucleatA', 2)

  do ijk = 1, ijkmax
    !
    !--- lhv
    qlevp = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
    !--- supersaturation
    psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
          * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
    ssliq = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
    psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
          * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
    ssice = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
    ssliq = qvap(ijk)/ssliq-1.0_RP
    ssice = qvap(ijk)/ssice-1.0_RP

    if ( ssliq <= 0.0_RP ) cycle
    !--- use for aerosol coupled model
    !--- mass -> number
    do n = 1, nccn
      gan( n ) = ga( n,ijk )*rexpxactr( n )
    enddo

    acoef = 2.0_RP*sigma/rvap/rhow/temp(ijk)   ! A in (A.11) of Suzuki (2004)
    bcoef = vhfct* rhoa/rhow * emwtr/emaer     ! B in (A.11) of Suzuki (2004)

    !--- relationship of bin number
    do n = 1, nccn
      ractr = ( expxactr( n )*ThirdovForth/pi/rhoa )**( OneovThird )
      rcld  = sqrt( 3.0_RP*bcoef*ractr*ractr*ractr / acoef )
      xcld  = log( rhow * 4.0_RP*pi*OneovThird*rcld*rcld*rcld )
     if ( flg_nucl ) then
      ncld( n ) = 1
     else
      ncld( n ) = int( ( xcld-xctr( 1 ) )/dxmic ) + 1
      ncld( n ) = min( max( ncld( n ),1 ),nbin )
     endif
    enddo

    !--- nucleation
    do n = nccn, 1, -1
        !--- super saturation
        psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
              * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
        ssliq = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
        psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
              * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
        ssice = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
        ssliq = qvap(ijk)/ssliq-1.0_RP
        ssice = qvap(ijk)/ssice-1.0_RP

      if ( ssliq <= 0.0_RP ) exit
      !--- use for aerosol coupled model
      acoef = 2.0_RP*sigma/rvap/rhow/temp(ijk) ! A in (A.11) of Suzuki (2004)
      rcrit = acoef*OneovThird * ( 4.0_RP/bcoef )**( OneovThird ) / ssliq**( TwoovThird ) ! r_{N,crit} in (A.11) of Suzuki (2004)
      xcrit = log( rhoa * 4.0_RP*pi*OneovThird * rcrit*rcrit*rcrit )
      ncrit = int( ( xcrit-xabnd( 1 ) )/dxaer ) + 1

      if ( n == ncrit ) then
        part = ( xabnd( ncrit+1 )-xcrit )/dxaer
      elseif ( n > ncrit ) then
        part = 1.0_RP
      else
        exit
      endif

      !--- calculate mass change
      nc = ncld( n )
      dmp = part*gan( n )*dxaer*expxctr( nc )
      dmp = min( dmp,qvap(ijk)*dens(ijk) )
      gc( nc,il,ijk ) = gc( nc,il,ijk ) + dmp/dxmic
      gan( n ) = gan( n ) - dmp/dxaer*rexpxctr( nc )
      gan( n ) = max( gan( n ), 0.0_RP )
      qvap(ijk) = qvap(ijk) - dmp/dens(ijk)
      qvap(ijk) = max( qvap(ijk),0.0_RP )
      temp(ijk) = temp(ijk) + dmp/dens(ijk)*qlevp/cp
    enddo

    !--- number -> mass
    do n = 1, nccn
      ga( n,ijk ) = gan( n )*expxactr( n )
    enddo

  enddo

    call PROF_rapend  ('_SBM_NucleatA', 2)

    return
  end subroutine nucleata

  !-----------------------------------------------------------------------------
  subroutine cndevpsbl( &
       ijkmax, &
       dens,   &
       pres,   &
       temp,   &
       qvap,   &
       gc,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

    real(RP) :: regene_gcn(ijkmax)
    !---------------------------------------------------------------------------

    call liqphase( ijkmax,        & ! [IN]
                   dens(:),       & ! [IN]
                   pres(:),       & ! [IN]
                   temp(:),       & ! [INOUT]
                   qvap(:),       & ! [INOUT]
                   gc  (:,:,:),   & ! [INOUT]
                   regene_gcn(:), & ! [OUT]
                   dtime          ) ! [IN]

    call icephase( ijkmax,        & ! [IN]
                   dens(:),       & ! [IN]
                   pres(:),       & ! [IN]
                   temp(:),       & ! [INOUT]
                   qvap(:),       & ! [INOUT]
                   gc  (:,:,:),   & ! [INOUT]
                   dtime          ) ! [IN]

    call mixphase( ijkmax,        & ! [IN]
                   dens(:),       & ! [IN]
                   pres(:),       & ! [IN]
                   temp(:),       & ! [INOUT]
                   qvap(:),       & ! [INOUT]
                   gc  (:,:,:),   & ! [INOUT]
                   dtime          ) ! [IN]

    return
  end subroutine cndevpsbl

  !-----------------------------------------------------------------------------
  subroutine cndevpsbla( &
       ijkmax, &
       dens,   &
       pres,   &
       temp,   &
       qvap,   &
       gc,     &
       ga,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(RP), intent(inout) :: ga  (nccn     ,ijkmax) ! Mass size distribution function of aerosol
    real(DP), intent(in)    :: dtime                  ! Time step interval

    real(RP) :: regene_gcn(ijkmax)
    !---------------------------------------------------------------------------

    call liqphase( ijkmax,        & ! [IN]
                   dens(:),       & ! [IN]
                   pres(:),       & ! [IN]
                   temp(:),       & ! [INOUT]
                   qvap(:),       & ! [INOUT]
                   gc  (:,:,:),   & ! [INOUT]
                   regene_gcn(:), & ! [OUT]
                   dtime          ) ! [IN]

    ! regeneration of aerosol
    if ( flg_regeneration ) then
       call faero( ijkmax,        & ! [IN]
                   regene_gcn(:), & ! [IN]
                   ga(:,:)        ) ! [INOUT]
    endif

    call icephase( ijkmax,        & ! [IN]
                   dens(:),       & ! [IN]
                   pres(:),       & ! [IN]
                   temp(:),       & ! [INOUT]
                   qvap(:),       & ! [INOUT]
                   gc  (:,:,:),   & ! [INOUT]
                   dtime          ) ! [IN]

    call mixphase( ijkmax,        & ! [IN]
                   dens(:),       & ! [IN]
                   pres(:),       & ! [IN]
                   temp(:),       & ! [INOUT]
                   qvap(:),       & ! [INOUT]
                   gc  (:,:,:),   & ! [INOUT]
                   dtime          ) ! [IN]

    return
  end subroutine cndevpsbla

  !-----------------------------------------------------------------------------
  subroutine liqphase( &
       ijkmax,     &
       dens,       &
       pres,       &
       temp,       &
       qvap,       &
       gc,         &
       regene_gcn, &
       dtime       )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(RP), intent(out)   :: regene_gcn(ijkmax)     ! mass of regenerated aerosol
    real(DP), intent(in)    :: dtime                  ! Time step interval

  integer  :: n, myu, ncount
  integer  :: nloop(ijkmax)                    ! number of fractional step for condensation
  real(RP) :: gclold(ijkmax)
  real(RP) :: gclnew(ijkmax)
  real(RP) :: cndmss(ijkmax)
  real(RP) :: dtcnd(ijkmax)                    ! dt for condensation with considering CFL condition
  real(RP) :: gtliq(ijkmax)                    ! G of eq. (A.17) of Suzuki (2004)
  real(RP) :: qlevp(ijkmax)                    ! LH
  real(RP) :: cefliq, a, sliqtnd
  real(RP) :: sumliq(ijkmax), umax(ijkmax)
  real(RP) :: ssliq(ijkmax), ssice(ijkmax) ! super saturation
  real(RP) :: gcn( -1:nbin+2,nspc,ijkmax )               ! size distribution function(hydrometeor): number
!  real(RP) :: gcnold( nbin,ijkmax )
  real(RP), parameter :: cflfct = 0.50_RP              ! CFL limiter
!  real(RP) :: old_sum_gcn, new_sum_gcn
  integer  :: iflg( nspc,ijkmax )                ! flag whether calculation is conduct or not
  real(RP) :: csum( nspc,ijkmax )
  real(RP) :: f1, f2, emu, cefd, cefk, festl, psat
  real(RP), parameter :: afmyu = 1.72E-05_RP, bfmyu = 3.93E+2_RP
  real(RP), parameter :: cfmyu = 1.2E+02_RP, fct = 1.4E+3_RP
  real(RP) :: zerosw, qvtmp
  integer :: ijk, nn, mm, loopflg(ijkmax)

  !--- local for advection
  real(RP) :: uadv ( 0:nbin+2,nspc,ijkmax )
  real(RP) :: flq  ( 1:nbin+1,nspc,ijkmax )
  real(RP) :: acoef( 0:2,0:nbin+1,nspc,ijkmax )
  real(RP) :: crn  ( 0:nbin+2,nspc,ijkmax )
  real(RP) :: aip  ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: aim  ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: ai   ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: cmins, cplus
  integer :: nloopmax

    call PROF_rapstart('_SBM_Liqphase', 2)


  iflg(:,:) = 0
  csum(:,:) = 0.0_RP
  do ijk = 1, ijkmax

     do n = 1, nbin
       csum( il,ijk ) = csum( il,ijk ) + gc( n,il,ijk )*dxmic
     enddo

     if( csum( il,ijk ) > cldmin ) iflg( il,ijk ) = 1

  enddo

  nloop(:) = 0
  gclold(:) = 0.0_RP
  do ijk = 1, ijkmax

     do n = 1, nbin
       gclold(ijk) = gclold(ijk) + gc( n,il,ijk ) * dxmic
     enddo
     !
     !------- mass -> number
     do n = 1, nbin
       gcn( n,il,ijk ) = gc( n,il,ijk ) * rexpxctr( n )
     enddo
     gcn( -1,il,ijk ) = 0.0_RP
     gcn(  0,il,ijk ) = 0.0_RP
     gcn( nbin+1,il,ijk ) = 0.0_RP
     gcn( nbin+2,il,ijk ) = 0.0_RP
     !
     !--- lhv
     qlevp(ijk) = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
     !--- super saturation
     psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
           * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
     ssliq(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
     psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
           * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
     ssice(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
     ssliq(ijk) = qvap(ijk)/ssliq(ijk)-1.0_RP
     ssice(ijk) = qvap(ijk)/ssice(ijk)-1.0_RP

     emu = afmyu*( bfmyu/( temp(ijk)+cfmyu ) )*( temp(ijk)/tmlt )**1.50_RP
     cefd = emu/dens(ijk)
     cefk = fct*emu

     festl = esat0*exp( qlevp(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
     f1 = rvap*temp(ijk)/festl/cefd
     f2 = qlevp(ijk)/cefk/temp(ijk)*( qlevp(ijk)/rvap/temp(ijk) - 1.0_RP )
     gtliq(ijk) = 4.0_RP*pi/( f1+f2 )  !--- G of eq. (A.17) of Suzuki (2004)
     !------- CFL condition
     umax(ijk) = cbnd( 1,il )*rexpxbnd( 1 )*gtliq(ijk)*abs( ssliq(ijk) )
     dtcnd(ijk) = cflfct*dxmic/umax(ijk)
     nloop(ijk) = int( dtime/dtcnd(ijk) ) + 1
     dtcnd(ijk) = dtime / nloop(ijk)

     nloop(ijk) = nloop(ijk) * iflg( il,ijk ) !--- for determing trivial loop
  enddo
  nloopmax = maxval(nloop,1)

  !
  regene_gcn(:) = 0.0_RP
!OCL LOOP_NOFISSION
!OCL LOOP_NOINTERCHANGE
  do ncount = 1, nloopmax

    do ijk = 1, ijkmax
       loopflg(ijk) = min( 1, int(nloop(ijk)/ncount) )   ! 0 or 1
    enddo

!OCL LOOP_NOFUSION
    do ijk = 1, ijkmax
     do nn = 1, loopflg(ijk)

       !--- lhv
       qlevp(ijk) = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
       !----- matrix for supersaturation tendency
       psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
             * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
       ssliq(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
       psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
             * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
       ssice(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
       ssliq(ijk) = qvap(ijk)/ssliq(ijk)-1.0_RP
       ssice(ijk) = qvap(ijk)/ssice(ijk)-1.0_RP

       emu = afmyu*( bfmyu/( temp(ijk)+cfmyu ) )*( temp(ijk)/tmlt )**1.50_RP
       cefd = emu/dens(ijk)
       cefk = fct*emu
       festl = esat0*exp( qlevp(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
       f1 = rvap*temp(ijk)/festl/cefd
       f2 = qlevp(ijk)/cefk/temp(ijk)*( qlevp(ijk)/rvap/temp(ijk) - 1.0_RP )
       gtliq(ijk) = 4.0_RP*pi/( f1+f2 ) !--- G of eq. (A.17) of Suzuki (2004)

       sumliq(ijk) = 0.0_RP
!       old_sum_gcn = 0.0_RP
       do n = 1, nbin
         sumliq(ijk) = sumliq(ijk) + gcn( n,il,ijk )*cctr( n,il )*dxmic
       enddo

     enddo
    enddo

!OCL LOOP_NOFUSION
    do ijk = 1, ijkmax
     do nn = 1, loopflg(ijk)

       !----- supersaturation tendency
       zerosw = 0.5_RP + sign( 0.5_RP,qvap(ijk)-EPS )  !--- zerosw = 1 (qv>0), zerosw=0 (qv=0)
       qvtmp = qvap(ijk) * zerosw + ( qvap(ijk)+EPS ) * ( 1.0_RP-zerosw )
       cefliq = ( ssliq(ijk)+1.0_RP )*( 1.0_RP/qvtmp + qlevp(ijk)*qlevp(ijk)/cp/rvap/temp(ijk)/temp(ijk) )
       a = - cefliq*sumliq(ijk)*gtliq(ijk)/dens(ijk)   ! a of eq. (A.19) of Suzuki (2004)

       sliqtnd = zerosw * &
               ( ssliq(ijk) * ( exp( a*dtcnd(ijk) )-1.0_RP )/( a*dtcnd(ijk) ) &
               * ( 0.5_RP + sign( 0.5_RP,abs(a*dtcnd(ijk)-0.1_RP) ) ) &
               + ssliq(ijk) &
               * ( 0.5_RP - sign( 0.5_RP,abs(a*dtcnd(ijk)-0.1_RP) ) ) &
               ) &
               + ssliq(ijk) * ( 1.0_RP - zerosw )
       !
       !----- change of SDF
       do mm = 1, iflg( il,ijk )
         !--- advection speed
         do n = 1, nbin+1
           uadv( n,il,ijk ) = cbnd( n,il )*rexpxbnd( n )*gtliq(ijk)*sliqtnd  ! U of eq. (A.18) of Suzuki (2004)
         enddo
         uadv( 0     ,il,ijk ) = 0.0_RP
         uadv( nbin+2,il,ijk ) = 0.0_RP

!         do n = 1, nbin
!           gcnold( n ) = gcn( n,il,ijk )
!         enddo
       enddo

     enddo
    enddo

       call PROF_rapstart('_SBM_AdvLiq', 3)

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+2
             crn( n,myu,ijk ) = uadv( n,myu,ijk )*dtcnd(ijk)/dxmic
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             acoef(0,n,myu,ijk) = - ( gcn( n+1,myu,ijk )-26.0_RP*gcn( n,myu,ijk )+gcn( n-1,myu,ijk ) ) / 48.0_RP
             acoef(1,n,myu,ijk) =   ( gcn( n+1,myu,ijk )                         -gcn( n-1,myu,ijk ) ) / 16.0_RP
             acoef(2,n,myu,ijk) =   ( gcn( n+1,myu,ijk )- 2.0_RP*gcn( n,myu,ijk )+gcn( n-1,myu,ijk ) ) / 48.0_RP
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             cplus = 1.0_RP - ( crn(n+1,myu,ijk) + abs(crn(n+1,myu,ijk)) )

             aip(n,myu,ijk) = acoef(0,n,myu,ijk) * ( 1.0_RP-cplus**1 ) &
                            + acoef(1,n,myu,ijk) * ( 1.0_RP-cplus**2 ) &
                            + acoef(2,n,myu,ijk) * ( 1.0_RP-cplus**3 )

             aip(n,myu,ijk) = max( aip(n,myu,ijk), 0.0_RP )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             cmins = 1.0_RP - ( abs(crn(n,myu,ijk)) - crn(n,myu,ijk) )

             aim(n,myu,ijk) = acoef(0,n,myu,ijk) * ( 1.0_RP-cmins**1 ) &
                            - acoef(1,n,myu,ijk) * ( 1.0_RP-cmins**2 ) &
                            + acoef(2,n,myu,ijk) * ( 1.0_RP-cmins**3 )

             aim(n,myu,ijk) = max( aim(n,myu,ijk), 0.0_RP )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             ai(n,myu,ijk) = acoef(0,n,myu,ijk) * 2.0_RP &
                           + acoef(2,n,myu,ijk) * 2.0_RP

             ai(n,myu,ijk) = max( ai(n,myu,ijk), aip(n,myu,ijk)+aim(n,myu,ijk)+cldmin )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 1, nbin+1
             flq(n,myu,ijk) = ( aip(n-1,myu,ijk)/ai(n-1,myu,ijk)*gcn( n-1,myu,ijk ) &
                              - aim(n  ,myu,ijk)/ai(n  ,myu,ijk)*gcn( n  ,myu,ijk ) )*dxmic/dtcnd(ijk)
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          regene_gcn(ijk) = regene_gcn(ijk)+( -flq(1,myu,ijk)*dtcnd(ijk)/dxmic )*min( uadv(1,myu,ijk),0.0_RP )/uadv(1,myu,ijk)
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, il
       do mm  = 1, iflg( myu,ijk )
          do n = 1, nbin
             gcn( n,myu,ijk ) = gcn( n,myu,ijk ) - ( flq(n+1,myu,ijk)-flq(n,myu,ijk) )*dtcnd(ijk)/dxmic
          enddo
       enddo
       enddo
       enddo
       enddo

       call PROF_rapend  ('_SBM_AdvLiq', 3)

!OCL LOOP_NOFUSION
    do ijk = 1, ijkmax
     do nn = 1, loopflg(ijk)

       !----- new mass
       gclnew(ijk) = 0.0_RP
!       new_sum_gcn = 0.0_RP
       do n = 1, nbin
         gclnew(ijk) = gclnew(ijk) + gcn( n,il,ijk )*expxctr( n )
!         old_sum_gcn = old_sum_gcn + gcnold( n )*dxmic
!         new_sum_gcn = new_sum_gcn + gcn( n,ijk )*dxmic
       enddo

       gclnew(ijk) = gclnew(ijk)*dxmic
       !
       !----- change of humidity and temperature
       cndmss(ijk) = gclnew(ijk) - gclold(ijk)
       qvap(ijk) = qvap(ijk) - cndmss(ijk)/dens(ijk)
       temp(ijk) = temp(ijk) + cndmss(ijk)/dens(ijk)*qlevp(ijk)/cp
       !
       gclold(ijk) = gclnew(ijk)
       !
       !----- continue/end
     enddo
    enddo

  enddo   !nloop
  !
!OCL NORECURRENCE(gc)
  do ijk = 1, ijkmax
     !------- number -> mass
     do n = 1 , nbin
      gc( n,il,ijk ) = gcn( n,il,ijk )*expxctr( n )
     enddo
  enddo

!OCL NORECURRENCE(gc)
  do ijk = 1, ijkmax
      !--- lhv
      qlevp(ijk) = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
      do n = 1 , nbin
       if ( gc( n,il,ijk ) < 0.0_RP ) then
        cndmss(ijk) = -gc( n,il,ijk )
        gc( n,il,ijk ) = 0.0_RP
        qvap(ijk) = qvap(ijk) + cndmss(ijk)/dens(ijk)
        temp(ijk) = temp(ijk) - cndmss(ijk)/dens(ijk)*qlevp(ijk)/cp
       endif
      enddo
  enddo

    call PROF_rapend  ('_SBM_Liqphase', 2)

    return
  end subroutine liqphase

  !-----------------------------------------------------------------------------
  subroutine icephase( &
       ijkmax, &
       dens,   &
       pres,   &
       temp,   &
       qvap,   &
       gc,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

  integer :: myu, n, ncount
  integer :: nloop(ijkmax)                         !number of fractional step for condensation
  real(RP) :: gciold(ijkmax), gcinew(ijkmax)
  real(RP) :: dtcnd(ijkmax)                        ! dt for condensation with considering CFL condition
  real(RP) :: sblmss(ijkmax)
  real(RP) :: gtice(ijkmax), umax(ijkmax)
  real(RP) :: qlsbl(ijkmax)
  real(RP) :: cefice, d, uval, sicetnd
  real(RP) :: sumice(ijkmax)
  real(RP) :: ssliq(ijkmax), ssice(ijkmax)
  real(RP) :: gcn( -1:nbin+2,nspc,ijkmax )     ! size distribution function (Hydrometeor): number
  real(RP), parameter :: cflfct = 0.50_RP
  real(RP) :: dumm_regene(ijkmax)
  integer :: iflg( nspc,ijkmax )
  integer :: csum( nspc,ijkmax )
  real(RP) :: f1, f2, emu, cefd, cefk, festi, psat
  real(RP), parameter :: afmyu = 1.72E-05_RP, bfmyu = 3.93E+2_RP
  real(RP), parameter :: cfmyu = 1.2E+02_RP, fct = 1.4E+3_RP
  integer :: ijk, mm, nn, loopflg(ijkmax)
  real(RP) :: zerosw, qvtmp

  !--- local for advection
!  real(RP) :: qadv( -1:nbin+2 ), uadv( 0:nbin+2 )
  real(RP) :: uadv ( 0:nbin+2,nspc,ijkmax )
  real(RP) :: flq  ( 1:nbin+1,nspc,ijkmax )
  real(RP) :: acoef( 0:2,0:nbin+1,nspc,ijkmax )
  real(RP) :: crn  ( 0:nbin+2,nspc,ijkmax )
  real(RP) :: aip  ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: aim  ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: ai   ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: cmins, cplus
  integer :: nloopmax
    !---------------------------------------------------------------------------

    call PROF_rapstart('_SBM_Icephase', 2)

  iflg(:,:) = 0
  csum( :,: ) = 0.0_RP
  do ijk = 1, ijkmax

     do myu = 2, nspc
     do n = 1, nbin
       csum( myu,ijk ) = csum( myu,ijk ) + gc( n,myu,ijk )*dxmic
     enddo
     enddo

     do myu = 2, nspc
        if ( csum( myu,ijk ) > cldmin ) iflg( myu,ijk ) = 1
     enddo

  enddo

  gciold(:) = 0.0_RP
  nloop(:) = 0
  do ijk = 1, ijkmax

     !----- old mass
     do myu = 2, nspc
     do n = 1, nbin
       gciold(ijk) = gciold(ijk) + gc( n,myu,ijk )*dxmic
     enddo
     enddo

     !----- mass -> number
     do myu = 2, nspc
       do n = 1, nbin
         gcn( n,myu,ijk ) = gc( n,myu,ijk ) * rexpxctr( n )
       enddo
       gcn( -1,myu,ijk ) = 0.0_RP
       gcn(  0,myu,ijk ) = 0.0_RP
       gcn( nbin+1,myu,ijk ) = 0.0_RP
       gcn( nbin+2,myu,ijk ) = 0.0_RP
     enddo

     !--- lhv
     qlsbl(ijk) = CONST_LHS0 + ( CONST_CPvap - CONST_CI ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
     !--- supersaturation
     psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
           * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
     ssliq(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
     psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
           * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
     ssice(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
     ssliq(ijk) = qvap(ijk)/ssliq(ijk)-1.0_RP
     ssice(ijk) = qvap(ijk)/ssice(ijk)-1.0_RP

     emu = afmyu*( bfmyu/( temp(ijk)+cfmyu ) )*( temp(ijk)/tmlt )**1.50_RP
     cefd = emu/dens(ijk)
     cefk = fct*emu
     festi = esat0*exp( qlsbl(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
     f1 = rvap*temp(ijk)/festi/cefd
     f2 = qlsbl(ijk)/cefk/temp(ijk)*( qlsbl(ijk)/rvap/temp(ijk) - 1.0_RP )
     gtice(ijk) = 4.0_RP*pi/( f1+f2 )
     !----- CFL condition
     umax(ijk) = 0.0_RP
     do myu = 2, nspc
       uval = cbnd( 1,myu )*rexpxbnd( 1 )*gtice(ijk)*abs( ssice(ijk) )
       umax(ijk) = max( umax(ijk),uval )
     enddo

     dtcnd(ijk) = cflfct*dxmic/umax(ijk)
     nloop(ijk) = int( dtime/dtcnd(ijk) ) + 1
     dtcnd(ijk) = dtime/nloop(ijk)

     nloop(ijk) = nloop(ijk) * maxval( iflg( 2:nspc,ijk ) ) !--- for determing trivial loop

  enddo
  nloopmax = maxval(nloop,1)

!OCL LOOP_NOFISSION
!OCL LOOP_NOINTERCHANGE
  do ncount = 1, nloopmax

    do ijk = 1, ijkmax
      loopflg(ijk) = min( 1, int(nloop(ijk)/ncount) )   ! 0 or 1
    enddo

!OCL LOOP_NOFUSION
    do ijk = 1, ijkmax
      do nn = 1, loopflg(ijk)
       !--- lhv
       qlsbl(ijk) = CONST_LHS0 + ( CONST_CPvap - CONST_CI ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
       !----- matrix for supersaturation tendency
       psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
             * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
       ssliq(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
       psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
             * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
       ssice(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
       ssliq(ijk) = qvap(ijk)/ssliq(ijk)-1.0_RP
       ssice(ijk) = qvap(ijk)/ssice(ijk)-1.0_RP

       emu = afmyu*( bfmyu/( temp(ijk)+cfmyu ) )*( temp(ijk)/tmlt )**1.50_RP
       cefd = emu/dens(ijk)
       cefk = fct*emu
       festi = esat0*exp( qlsbl(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
       f1 = rvap*temp(ijk)/festi/cefd
       f2 = qlsbl(ijk)/cefk/temp(ijk)*( qlsbl(ijk)/rvap/temp(ijk) - 1.0_RP )
       gtice(ijk) = 4.0_RP*pi/( f1+f2 )   ! G of (A.17) of Suzuki (2004)

       sumice(ijk) = 0.0_RP
       do myu = 2, nspc
         do n = 1, nbin
           sumice(ijk) = sumice(ijk) + gcn( n,myu,ijk )*cctr( n,myu )*dxmic
         enddo
       enddo

      enddo
    enddo

!OCL LOOP_NOFUSION
    do ijk = 1, ijkmax
      do nn = 1, loopflg(ijk)

        !----- supersaturation tendency
        zerosw = 0.5_RP + sign( 0.5_RP,qvap(ijk)-EPS )  !--- zerosw = 1 (qv>0), zerosw=0 (qv=0)
        qvtmp = qvap(ijk) * zerosw + ( qvap(ijk)+EPS ) * ( 1.0_RP-zerosw )

        cefice = ( ssice(ijk)+1.0_RP )*( 1.0_RP/qvtmp + qlsbl(ijk)*qlsbl(ijk)/cp/rvap/temp(ijk)/temp(ijk) )
        d = - cefice*sumice(ijk)*gtice(ijk)/dens(ijk)  ! d of (A.19) of Suzuki (2004)

        sicetnd = zerosw * &
                ( ssice(ijk) * ( exp( d*dtcnd(ijk) )-1.0_RP )/( d*dtcnd(ijk) ) &
                * ( 0.5_RP + sign( 0.5_RP,abs(d*dtcnd(ijk)-0.1_RP) ) ) &
                + ssice(ijk) &
                * ( 0.5_RP - sign( 0.5_RP,abs(d*dtcnd(ijk)-0.1_RP) ) ) &
                ) &
                + ssice(ijk) * ( 1.0_RP - zerosw )
       !
       !----- change of SDF
       do myu = 2, nspc
       do mm = 1, iflg( myu,ijk )
         !--- advection speed
         do n = 1, nbin+1
           uadv( n,myu,ijk ) = cbnd( n,myu )*rexpxbnd( n )*gtice(ijk)*sicetnd  ! U of eq. (A.18) of Suzuki (2004)
         enddo
         uadv( 0,     myu,ijk ) = 0.0_RP
         uadv( nbin+2,myu,ijk ) = 0.0_RP
       enddo
       enddo

      enddo
    enddo

       call PROF_rapstart('_SBM_AdvIce', 3)

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+2
             crn( n,myu,ijk ) = uadv( n,myu,ijk )*dtcnd(ijk)/dxmic
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             acoef(0,n,myu,ijk) = - ( gcn( n+1,myu,ijk )-26.0_RP*gcn( n,myu,ijk )+gcn( n-1,myu,ijk ) ) / 48.0_RP
             acoef(1,n,myu,ijk) =   ( gcn( n+1,myu,ijk )                         -gcn( n-1,myu,ijk ) ) / 16.0_RP
             acoef(2,n,myu,ijk) =   ( gcn( n+1,myu,ijk )- 2.0_RP*gcn( n,myu,ijk )+gcn( n-1,myu,ijk ) ) / 48.0_RP
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             cplus = 1.0_RP - ( crn(n+1,myu,ijk) + abs(crn(n+1,myu,ijk)) )

             aip(n,myu,ijk) = acoef(0,n,myu,ijk) * ( 1.0_RP-cplus**1 ) &
                            + acoef(1,n,myu,ijk) * ( 1.0_RP-cplus**2 ) &
                            + acoef(2,n,myu,ijk) * ( 1.0_RP-cplus**3 )

             aip(n,myu,ijk) = max( aip(n,myu,ijk), 0.0_RP )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             cmins = 1.0_RP - ( abs(crn(n,myu,ijk)) - crn(n,myu,ijk) )

             aim(n,myu,ijk) = acoef(0,n,myu,ijk) * ( 1.0_RP-cmins**1 ) &
                            - acoef(1,n,myu,ijk) * ( 1.0_RP-cmins**2 ) &
                            + acoef(2,n,myu,ijk) * ( 1.0_RP-cmins**3 )

             aim(n,myu,ijk) = max( aim(n,myu,ijk), 0.0_RP )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             ai(n,myu,ijk) = acoef(0,n,myu,ijk) * 2.0_RP &
                           + acoef(2,n,myu,ijk) * 2.0_RP

             ai(n,myu,ijk) = max( ai(n,myu,ijk), aip(n,myu,ijk)+aim(n,myu,ijk)+cldmin )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 1, nbin+1
             flq(n,myu,ijk) = ( aip(n-1,myu,ijk)/ai(n-1,myu,ijk)*gcn( n-1,myu,ijk ) &
                              - aim(n  ,myu,ijk)/ai(n  ,myu,ijk)*gcn( n  ,myu,ijk ) )*dxmic/dtcnd(ijk)
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
           dumm_regene(ijk) = dumm_regene(ijk)+( -flq(1,myu,ijk)*dtcnd(ijk)/dxmic )*min( uadv(1,myu,ijk),0.0_RP )/uadv(1,myu,ijk)
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 2, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 1, nbin
             gcn( n,myu,ijk ) = gcn( n,myu,ijk ) - ( flq(n+1,myu,ijk)-flq(n,myu,ijk) )*dtcnd(ijk)/dxmic
          enddo
       enddo
       enddo
       enddo
       enddo

       call PROF_rapend  ('_SBM_AdvIce', 3)

!OCL LOOP_NOFUSION
    do ijk = 1, ijkmax
      do nn = 1, loopflg(ijk)

       !----- new mass
       gcinew(ijk) = 0.0_RP
       do n = 1, nbin
       do myu = 2, nspc
         gcinew(ijk) = gcinew(ijk) + gcn( n,myu,ijk )*expxctr( n )*dxmic
       enddo
       enddo
       !
       !----- change of humidity and temperature
       sblmss(ijk) = gcinew(ijk) - gciold(ijk)
       qvap(ijk) = qvap(ijk) - sblmss(ijk)/dens(ijk)
       temp(ijk) = temp(ijk) + sblmss(ijk)/dens(ijk)*qlsbl(ijk)/cp

       gciold(ijk) = gcinew(ijk)
       !
       !----- continue / end
     enddo
    enddo

  enddo ! nloop
  !
!OCL NORECURRENCE(gc)
  do ijk = 1, ijkmax
    !------- number -> mass
    do myu = 2, nspc
    do n = 1, nbin
      gc( n,myu,ijk ) = gcn( n,myu,ijk )*expxctr( n )
    enddo
    enddo
  enddo

    call PROF_rapend  ('_SBM_Icephase', 2)

    return
  end subroutine icephase

  !-----------------------------------------------------------------------------
  subroutine mixphase( &
       ijkmax, &
       dens,   &
       pres,   &
       temp,   &
       qvap,   &
       gc,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

  integer :: n, myu, mm, ncount
  integer :: nloop(ijkmax)
  real(RP) :: gclold(ijkmax), gclnew(ijkmax)
  real(RP) :: gciold(ijkmax), gcinew(ijkmax)
  real(RP) :: cndmss(ijkmax), sblmss(ijkmax)
  real(RP) :: gtliq(ijkmax), gtice(ijkmax)
  real(RP) :: umax(ijkmax), uval, dtcnd(ijkmax)
  real(RP) :: qlevp(ijkmax), qlsbl(ijkmax)
  real(RP) :: cef1, cef2, cef3, cef4, a, b, c, d
  real(RP) :: rmdplus, rmdmins, ssplus, ssmins, tplus, tmins
  real(RP) :: sliqtnd, sicetnd
  real(RP) :: ssliq(ijkmax), ssice(ijkmax)
  real(RP) :: sumliq(ijkmax), sumice(ijkmax)
  real(RP) :: gcn( -1:nbin+2,nspc,ijkmax )
  real(RP), parameter :: cflfct = 0.50_RP
  real(RP) :: dumm_regene(ijkmax)
  real(RP) :: csum( nspc,ijkmax )
  integer :: iflg( nspc,ijkmax )
  real(RP) :: f1, f2, emu, cefd, cefk, festl, festi, psat
  real(RP), parameter :: afmyu = 1.72E-05_RP, bfmyu = 3.93E+2_RP
  real(RP), parameter :: cfmyu = 1.2E+02_RP, fct = 1.4E+3_RP
  real(RP) :: qvtmp, zerosw
  integer :: ijk, nn, loopflg(ijkmax)

  !--- local for advection
!  real(RP) :: qadv( -1:nbin+2 ), uadv( 0:nbin+2 )
  real(RP) :: uadv ( 0:nbin+2,nspc,ijkmax )
  real(RP) :: flq  ( 1:nbin+1,nspc,ijkmax )
  real(RP) :: acoef( 0:2,0:nbin+1,nspc,ijkmax )
  real(RP) :: crn  ( 0:nbin+2,nspc,ijkmax )
  real(RP) :: aip  ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: aim  ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: ai   ( 0:nbin+1,nspc,ijkmax )
  real(RP) :: cmins, cplus
  integer :: nloopmax

    call PROF_rapstart('_SBM_Mixphase', 2)

  dumm_regene( : ) = 0.0_RP
  iflg( :,: ) = 0
  csum( :,: ) = 0.0_RP
  do ijk = 1, ijkmax

      do myu = 1, nspc
      do n = 1, nbin
        csum( myu,ijk ) = csum( myu,ijk )+gc( n,myu,ijk )*dxmic
      enddo
      enddo

      do myu = 1, nspc
        if ( csum( myu,ijk ) > cldmin ) iflg( myu,ijk ) = 1
      enddo
  enddo

  gclold(:) = 0.0_RP
  gciold(:) = 0.0_RP
  nloop(:) = 0
  do ijk = 1, ijkmax
      !----- old mass
      do n = 1, nbin
        gclold(ijk) = gclold(ijk) + gc( n,il,ijk )*dxmic
      enddo

      do myu = 2, nspc
      do n = 1, nbin
        gciold(ijk) = gciold(ijk) + gc( n,myu,ijk )*dxmic
      enddo
      enddo

      !----- mass -> number
      do myu = 1, nspc
       do n = 1, nbin
        gcn( n,myu,ijk ) = gc( n,myu,ijk ) * rexpxctr( n )
       enddo
       gcn( -1,myu,ijk ) = 0.0_RP
       gcn(  0,myu,ijk ) = 0.0_RP
       gcn( nbin+1,myu,ijk ) = 0.0_RP
       gcn( nbin+2,myu,ijk ) = 0.0_RP
      enddo

      !-- thermodyn
      qlevp(ijk) = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
      qlsbl(ijk) = CONST_LHS0 + ( CONST_CPvap - CONST_CI ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
      !-- supersaturation
      psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
            * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
      ssliq(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
      psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
            * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
      ssice(ijk) = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
      ssliq(ijk) = qvap(ijk)/ssliq(ijk)-1.0_RP
      ssice(ijk) = qvap(ijk)/ssice(ijk)-1.0_RP

      emu = afmyu*( bfmyu/( temp(ijk)+cfmyu ) )*( temp(ijk)/tmlt )**1.50_RP
      cefd = emu/dens(ijk)
      cefk = fct*emu

      festl = esat0*exp( qlevp(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
      f1 = rvap*temp(ijk)/festl/cefd
      f2 = qlevp(ijk)/cefk/temp(ijk)*( qlevp(ijk)/rvap/temp(ijk) - 1.0_RP )
      gtliq(ijk) = 4.0_RP*pi/( f1+f2 )

      festi = esat0*exp( qlsbl(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
      f1 = rvap*temp(ijk)/festi/cefd
      f2 = qlsbl(ijk)/cefk/temp(ijk)*( qlsbl(ijk)/rvap/temp(ijk) - 1.0_RP )
      gtice(ijk) = 4.0_RP*pi/( f1+f2 )

      !----- CFL condition
      umax(ijk) = cbnd( 1,il )*rexpxbnd( 1 )*gtliq(ijk)*abs( ssliq(ijk) )
      do myu = 2, nspc
        uval = cbnd( 1,myu )*rexpxbnd( 1 )*gtice(ijk)*abs( ssice(ijk) )
        umax(ijk) = max( umax(ijk),uval )
      enddo

      dtcnd(ijk) = cflfct*dxmic/umax(ijk)
      nloop(ijk) = int( dtime/dtcnd(ijk) ) + 1
      dtcnd(ijk) = dtime/nloop(ijk)

      nloop(ijk) = nloop(ijk) * iflg( il,ijk ) !--- for determing trivial loop
      nloop(ijk) = nloop(ijk) * maxval( iflg( 2:nspc,ijk ) ) !--- for determing trivial loop
!      nloop(ijk) = nloop(ijk) * maxval( iflg( 1:nspc,ijk ) ) !--- for determing trivial loop
  enddo
  nloopmax = maxval(nloop,1)


!OCL LOOP_NOFISSION
!OCL LOOP_NOINTERCHANGE
  do ncount = 1, nloopmax


     do ijk = 1, ijkmax
      loopflg(ijk) = min( 1, int(nloop(ijk)/ncount) )   ! 0 or 1
     enddo

!OCL LOOP_NOFUSION
      do ijk = 1, ijkmax
       do nn = 1, loopflg(ijk)
         !-- thermodyn
         qlevp(ijk) = CONST_LHV0 + ( CONST_CPvap - CONST_CL ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
         qlsbl(ijk) = CONST_LHS0 + ( CONST_CPvap - CONST_CI ) * ( temp(ijk) - CONST_TEM00 ) * flg_thermodyn
         !-- matrix for supersaturation tendency
         psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
               * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
         ssliq(ijk) = qvap(ijk) / ( CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat ) )-1.0_RP
         psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
               * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
         ssice(ijk) = qvap(ijk) / ( CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat ) )-1.0_RP

         emu = afmyu*( bfmyu/( temp(ijk)+cfmyu ) )*( temp(ijk)/tmlt )**1.50_RP
         cefd = emu/dens(ijk)
         cefk = fct*emu
         festl = esat0*exp( qlevp(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
         f1 = rvap*temp(ijk)/festl/cefd
         f2 = qlevp(ijk)/cefk/temp(ijk)*( qlevp(ijk)/rvap/temp(ijk) - 1.0_RP )
         gtliq(ijk) = 4.0_RP*pi/( f1+f2 )
         festi = esat0*exp( qlsbl(ijk)/rvap*( 1.0_RP/temp0 - 1.0_RP/temp(ijk) ) )
         f1 = rvap*temp(ijk)/festi/cefd
         f2 = qlsbl(ijk)/cefk/temp(ijk)*( qlsbl(ijk)/rvap/temp(ijk) - 1.0_RP )
         gtice(ijk) = 4.0_RP*pi/( f1+f2 )

         sumliq(ijk) = 0.0_RP
         do n = 1, nbin
           sumliq(ijk) = sumliq(ijk) + gcn( n,il,ijk )*cctr( n,il )*dxmic
         enddo

         sumice(ijk) = 0.0_RP
         do myu = 2, nspc
         do n = 1, nbin
           sumice(ijk) = sumice(ijk) + gcn( n,myu,ijk )*cctr( n,myu )*dxmic
         enddo
         enddo
       enddo
      enddo


!OCL LOOP_NOFUSION
      do ijk = 1, ijkmax
       do nn = 1, loopflg(ijk)
         zerosw = 0.5_RP + sign( 0.5_RP,qvap(ijk)-EPS )  !--- zerosw = 1 (qv>0), zerosw=0 (qv=0)
         qvtmp = qvap(ijk) * zerosw + ( qvap(ijk)+EPS ) * ( 1.0_RP-zerosw )
         cef1 = ( ssliq(ijk)+1.0_RP )*( 1.0_RP/qvtmp + qlevp(ijk)/rvap/temp(ijk)/temp(ijk)*qlevp(ijk)/cp )
         cef2 = ( ssliq(ijk)+1.0_RP )*( 1.0_RP/qvtmp + qlevp(ijk)/rvap/temp(ijk)/temp(ijk)*qlsbl(ijk)/cp )
         cef3 = ( ssice(ijk)+1.0_RP )*( 1.0_RP/qvtmp + qlsbl(ijk)/rvap/temp(ijk)/temp(ijk)*qlevp(ijk)/cp )
         cef4 = ( ssice(ijk)+1.0_RP )*( 1.0_RP/qvtmp + qlsbl(ijk)/rvap/temp(ijk)/temp(ijk)*qlsbl(ijk)/cp )

         a = - cef1*sumliq(ijk)*gtliq(ijk)/dens(ijk)  ! a of (A.19) of Suzuki (2004)
         b = - cef2*sumice(ijk)*gtice(ijk)/dens(ijk)  ! b of (A.19) of Suzuki (2004)
         c = - cef3*sumliq(ijk)*gtliq(ijk)/dens(ijk)  ! c of (A.19) of Suzuki (2004)
         d = - cef4*sumice(ijk)*gtice(ijk)/dens(ijk)  ! d of (A.19) of Suzuki (2004)

         !--- eigenvalues
         rmdplus = ( ( a+d ) + sqrt( ( a-d )**2 + 4.0_RP*b*c ) ) * 0.50_RP
         rmdmins = ( ( a+d ) - sqrt( ( a-d )**2 + 4.0_RP*b*c ) ) * 0.50_RP

         !--- supersaturation tendency
         ssplus = ( ( rmdmins-a )*ssliq(ijk) - b*ssice(ijk) )/b/( rmdmins-rmdplus )
         ssmins = ( ( a-rmdplus )*ssliq(ijk) + b*ssice(ijk) )/b/( rmdmins-rmdplus )

         tplus = ( exp( rmdplus*dtcnd(ijk) )-1.0_RP )/( rmdplus*dtcnd(ijk) ) &
               * ( 0.5_RP + sign( 0.5_RP,abs(rmdplus*dtcnd(ijk)-0.1_RP) ) ) &
               + 1.0_RP &
               * ( 0.5_RP - sign( 0.5_RP,abs(rmdplus*dtcnd(ijk)-0.1_RP) ) )
         tmins = ( exp( rmdmins*dtcnd(ijk) )-1.0_RP )/( rmdmins*dtcnd(ijk) ) &
               * ( 0.5_RP + sign( 0.5_RP,abs(rmdmins*dtcnd(ijk)-0.1_RP) ) ) &
               + 1.0_RP  &
               * ( 0.5_RP - sign( 0.5_RP,abs(rmdmins*dtcnd(ijk)-0.1_RP) ) )

         sliqtnd = ssliq(ijk) * ( 1.0_RP - zerosw ) &
                 + zerosw * &
                 ( b*tplus*ssplus + b*tmins*ssmins ) ! sliwtnd in page 116 of Suzuki (2004)
         sicetnd = ssice(ijk) * ( 1.0_RP - zerosw ) &
                 + zerosw * &
                 ( ( rmdplus-a )*tplus*ssplus  &   ! sicetnd in page 116 of Suzuki (2004)
                 + ( rmdmins-a )*tmins*ssmins  )

         !--- change of SDF
         do myu = 1, nspc
         do mm = 1, iflg( myu,ijk )
           !--- advection speed
           do n = 1, nbin+1
             !--- myu = 1 -> ssliq, myu > 1 -> ssice
             uadv( n,myu,ijk ) = cbnd( n,myu )*rexpxbnd( n )*gtliq(ijk)*sliqtnd   & ! U of eq. (A.18) of Suzuki (2004)
                           * ( 0.5_RP - sign( 0.5_RP,real(myu)-1.5_RP) ) &
                           + cbnd( n,myu )*rexpxbnd( n )*gtice(ijk)*sicetnd   & ! U of eq. (A.18) of Suzuki (2004)
                           * ( 0.5_RP + sign( 0.5_RP,real(myu)-1.5_RP) )
           enddo
           uadv( 0,     myu,ijk ) = 0.0_RP
           uadv( nbin+2,myu,ijk ) = 0.0_RP
         enddo
         enddo
        enddo
      enddo

       call PROF_rapstart('_SBM_AdvMix', 3)

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+2
             crn( n,myu,ijk ) = uadv( n,myu,ijk )*dtcnd(ijk)/dxmic
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             acoef(0,n,myu,ijk) = - ( gcn( n+1,myu,ijk )-26.0_RP*gcn( n,myu,ijk )+gcn( n-1,myu,ijk ) ) / 48.0_RP
             acoef(1,n,myu,ijk) =   ( gcn( n+1,myu,ijk )                         -gcn( n-1,myu,ijk ) ) / 16.0_RP
             acoef(2,n,myu,ijk) =   ( gcn( n+1,myu,ijk )- 2.0_RP*gcn( n,myu,ijk )+gcn( n-1,myu,ijk ) ) / 48.0_RP
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             cplus = 1.0_RP - ( crn(n+1,myu,ijk) + abs(crn(n+1,myu,ijk)) )

             aip(n,myu,ijk) = acoef(0,n,myu,ijk) * ( 1.0_RP-cplus**1 ) &
                            + acoef(1,n,myu,ijk) * ( 1.0_RP-cplus**2 ) &
                            + acoef(2,n,myu,ijk) * ( 1.0_RP-cplus**3 )

             aip(n,myu,ijk) = max( aip(n,myu,ijk), 0.0_RP )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             cmins = 1.0_RP - ( abs(crn(n,myu,ijk)) - crn(n,myu,ijk) )

             aim(n,myu,ijk) = acoef(0,n,myu,ijk) * ( 1.0_RP-cmins**1 ) &
                            - acoef(1,n,myu,ijk) * ( 1.0_RP-cmins**2 ) &
                            + acoef(2,n,myu,ijk) * ( 1.0_RP-cmins**3 )

             aim(n,myu,ijk) = max( aim(n,myu,ijk), 0.0_RP )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 0, nbin+1
             ai(n,myu,ijk) = acoef(0,n,myu,ijk) * 2.0_RP &
                           + acoef(2,n,myu,ijk) * 2.0_RP

             ai(n,myu,ijk) = max( ai(n,myu,ijk), aip(n,myu,ijk)+aim(n,myu,ijk)+cldmin )
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 1, nbin+1
             flq(n,myu,ijk) = ( aip(n-1,myu,ijk)/ai(n-1,myu,ijk)*gcn( n-1,myu,ijk ) &
                              - aim(n  ,myu,ijk)/ai(n  ,myu,ijk)*gcn( n  ,myu,ijk ) )*dxmic/dtcnd(ijk)
          enddo
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
           dumm_regene(ijk) = dumm_regene(ijk)+( -flq(1,myu,ijk)*dtcnd(ijk)/dxmic )*min( uadv(1,myu,ijk),0.0_RP )/uadv(1,myu,ijk)
       enddo
       enddo
       enddo
       enddo

       do ijk = 1, ijkmax
       do nn  = 1, loopflg(ijk)
       do myu = 1, nspc
       do mm  = 1, iflg( myu,ijk )
          do n = 1, nbin
             gcn( n,myu,ijk ) = gcn( n,myu,ijk ) - ( flq(n+1,myu,ijk)-flq(n,myu,ijk) )*dtcnd(ijk)/dxmic
          enddo
       enddo
       enddo
       enddo
       enddo

       call PROF_rapend  ('_SBM_AdvMix', 3)

!OCL LOOP_NOFUSION
      do ijk = 1, ijkmax
       do nn = 1, loopflg(ijk)
         !--- new mass
         gclnew(ijk) = 0.0_RP
         do n = 1, nbin
           gclnew(ijk) = gclnew(ijk) + gcn( n,il,ijk )*expxctr( n )*dxmic
         enddo

         gcinew(ijk) = 0.0_RP
         do myu = 2, nspc
         do n = 1, nbin
           gcinew(ijk) = gcinew(ijk) + gcn( n,myu,ijk )*expxctr( n )*dxmic
         enddo
         enddo

         !--- change of humidity and temperature
         cndmss(ijk) = gclnew(ijk) - gclold(ijk)
         sblmss(ijk) = gcinew(ijk) - gciold(ijk)

         qvap(ijk) = qvap(ijk) - ( cndmss(ijk)+sblmss(ijk) )/dens(ijk)
         temp(ijk) = temp(ijk) + ( cndmss(ijk)*qlevp(ijk)+sblmss(ijk)*qlsbl(ijk) )/dens(ijk)/cp

         gclold(ijk) = gclnew(ijk)
         gciold(ijk) = gcinew(ijk)
       enddo
      enddo

  enddo   ! ncount

!OCL NORECURRENCE(gc)
  do ijk = 1, ijkmax
    !----- number -> mass
    do myu = 1, nspc
    do n = 1, nbin
      gc( n,myu,ijk ) = gcn( n,myu,ijk )*expxctr( n )
    enddo
    enddo
  enddo

  if ( .not. flg_regeneration ) then
   dumm_regene(:) = 0.0_RP
  endif

    call PROF_rapend  ('_SBM_Mixphase', 2)

    return
  end subroutine mixphase

  !-----------------------------------------------------------------------------
  subroutine ice_nucleat( &
       ijkmax,     &
       num_cold,   &
       index_cold, &
       dens,       &
       pres,       &
       temp,       &
       qvap,       &
       gc,         &
       dtime       )
    implicit none

    integer,  intent(in)    :: ijkmax
    integer,  intent(in)    :: num_cold
    integer,  intent(in)    :: index_cold(ijkmax)
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(in)    :: pres(ijkmax)           ! Pressure          [Pa]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: qvap(ijkmax)           ! Specific humidity [kg/kg]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

  real(RP) :: ssliq, ssice
  real(RP) :: numin, tdel, qdel
!  real(RP), parameter :: n0 = 1.E+3_RP                       ! N_{IN0} in page 112 of Suzuki (2004)
  real(RP), parameter :: acoef = -0.639_RP, bcoef = 12.96_RP ! A and B in paeg 112 of Suzuki (2004)
  ! threshould to determine the type of freezed hydrometeor (the detail is described in page 113 of Suzuki(2004))
  real(RP), parameter :: tcolmu = 269.0_RP, tcolml = 265.0_RP! -4[degC], -8[degC]
  real(RP), parameter :: tdendu = 257.0_RP, tdendl = 255.0_RP! -14[degC], -18[degC]
  real(RP), parameter :: tplatu = 250.6_RP                   ! -22.4[degC]
  !
  real(RP) :: psat
  integer :: ijk, indirect


    call PROF_rapstart('_SBM_IceNucleat', 2)

  do indirect = 1, num_cold
     ijk = index_cold(indirect)

    !--- supersaturation
    psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_liq   &
          * exp( LovR_liq * ( RTEM00 - 1.0_RP/temp(ijk) ) )
    ssliq = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
    psat  = CONST_PSAT0 * ( temp(ijk) * RTEM00 )**CPovR_ice   &
          * exp( LovR_ice * ( RTEM00 - 1.0_RP/temp(ijk) ) )
    ssice = CONST_EPSvap * psat / ( pres(ijk) - ( 1.0_RP-CONST_EPSvap ) * psat )
    ssliq = qvap(ijk)/ssliq-1.0_RP
    ssice = qvap(ijk)/ssice-1.0_RP

    if( ssice <= 0.0_RP ) cycle

    numin = bcoef * exp( acoef + bcoef * ssice )
    numin = numin * expxctr( 1 )/dxmic
    numin = min( numin,qvap(ijk)*dens(ijk) )
    !--- -4 [deg] > T >= -8 [deg] and T < -22.4 [deg] -> column
    if ( temp(ijk) <= tplatu .OR. ( temp(ijk) >= tcolml .AND. temp(ijk) < tcolmu ) ) then
     gc( 1,ic,ijk ) = gc( 1,ic,ijk ) + numin
    !--- -14 [deg] > T >= -18 [deg] -> dendrite
    elseif( temp(ijk) <= tdendu .AND. temp(ijk) >= tdendl ) then
     gc( 1,id,ijk ) = gc( 1,id,ijk ) + numin
    !--- else -> plate
    else
     gc( 1,ip,ijk ) = gc( 1,ip,ijk ) + numin
    endif

    tdel = numin/dens(ijk)*qlmlt/cp
    temp(ijk) = temp(ijk) + tdel
    qdel = numin/dens(ijk)
    qvap(ijk) = qvap(ijk) - qdel

  enddo


    call PROF_rapend  ('_SBM_IceNucleat', 2)

    return
  end subroutine ice_nucleat

  !-----------------------------------------------------------------------------
  subroutine freezing( &
       ijkmax,     &
       num_cold,   &
       index_cold, &
       dens,       &
       temp,       &
       gc,         &
       dtime       )
    implicit none

    integer,  intent(in)    :: ijkmax
    integer,  intent(in)    :: num_cold
    integer,  intent(in)    :: index_cold(ijkmax)
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

  integer :: nbound, n
  real(RP) :: xbound, tc, rate, dmp, frz, sumfrz, tdel
  real(RP), parameter :: coefa = 1.0E-01_RP   ! a_{fr} of eq.(3.19) of Suzuki (2004)
  real(RP), parameter :: coefb = 0.66_RP      ! b_{fr} of eq.(3.19) of Suzuki (2004)
  real(RP), parameter :: rbound = 2.0E-04_RP  ! 200 um
!  real(RP), parameter :: tthreth = 235.0_RP   ! -38 [degC] threshold for using Bigg's parameterization
!  real(RP), parameter :: ncoefim = 1.0E+7_RP  ! N_{im0} of eq.(3.18) of Suzuki (2004)
!  real(RP), parameter :: gamm = 3.3_RP        ! gamma of eq.(3.18) of Suzuki (2004)
  integer :: ijk, indirect

    call PROF_rapstart('_SBM_Freezing', 2)

  do indirect = 1, num_cold
     ijk = index_cold(indirect)

!      if ( temp <= tthreth ) then !--- Bigg (1975)
      xbound = log( rhow * 4.0_RP*pi/3.0_RP * rbound**3 )
      nbound = int( ( xbound-xbnd( 1 ) )/dxmic ) + 1

      tc = temp(ijk)-tmlt
      rate = coefa*exp( -coefb*tc )

      sumfrz = 0.0_RP
      do n = 1, nbound-1
        dmp = rate*expxctr( n )
        frz = gc( n,il,ijk )*( 1.0_RP-exp( -dmp*dtime ) )

        gc( n,il,ijk ) = gc( n,il,ijk ) - frz
        gc( n,il,ijk ) = max( gc( n,il,ijk ),0.0_RP )

        gc( n,ip,ijk ) = gc( n,ip,ijk ) + frz

        sumfrz = sumfrz + frz
      enddo
      do n = nbound, nbin
        dmp = rate*expxctr( n )
        frz = gc( n,il,ijk )*( 1.0_RP-exp( -dmp*dtime ) )

        gc( n,il,ijk ) = gc( n,il,ijk ) - frz
        gc( n,il,ijk ) = max( gc( n,il,ijk ),0.0_RP )

        gc( n,ih,ijk ) = gc( n,ih,ijk ) + frz

        sumfrz = sumfrz + frz
      enddo
!     elseif( temp > tthreth ) then !--- Vali (1975)
!
!      tc = temp-tmlt
!      dmp = ncoefim * ( 0.1_RP * ( tmlt-temp )**gamm )
!      sumfrz = 0.0_RP
!      do n = 1, nbin
!       frz = gc( (il-1)*nbin+n )*dmp*xctr( n )/dens
!       frz = max( frz,gc( (il-1)*nbin+n )*dmp*xctr( n )/dens )
!       gc( (il-1)*nbin+n ) = gc( (il-1)*nbin+n ) - frz
!       if ( n >= nbound ) then
!        gc( (ih-1)*nbin+n ) = gc( (ih-1)*nbin+n ) + frz
!       else
!        gc( (ip-1)*nbin+n ) = gc( (ip-1)*nbin+n ) + frz
!       endif
!
!       sumfrz = sumfrz + frz
!      enddo
!     endif
     sumfrz = sumfrz*dxmic

     tdel = sumfrz/dens(ijk)*qlmlt/cp
     temp(ijk) = temp(ijk) + tdel
  enddo

    call PROF_rapend  ('_SBM_Freezing', 2)

    return
  end subroutine freezing

  !-----------------------------------------------------------------------------
  subroutine melting( &
       ijkmax,     &
       num_warm,   &
       index_warm, &
       dens,       &
       temp,       &
       gc,         &
       dtime       )
    implicit none

    integer,  intent(in)    :: ijkmax
    integer,  intent(in)    :: num_warm
    integer,  intent(in)    :: index_warm(ijkmax)
    real(RP), intent(in)    :: dens(ijkmax)           ! Density           [kg/m3]
    real(RP), intent(inout) :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

  integer :: n, m
  real(RP) :: summlt, sumice, tdel
  integer :: ijk, indirect

    call PROF_rapstart('_SBM_Melting', 2)

  do indirect = 1, num_warm
     ijk = index_warm(indirect)

      summlt = 0.0_RP
      do n = 1, nbin
       sumice = 0.0_RP
       do m = ic, ih
        sumice = sumice + gc( n,m,ijk )
        gc( n,m,ijk ) = 0.0_RP
       enddo
       gc( n,il,ijk ) = gc( n,il,ijk ) + sumice
       summlt = summlt + sumice  !--- All freezed particle melt instantaneously
      enddo
      summlt = summlt*dxmic

      tdel = - summlt/dens(ijk)*qlmlt/cp
      temp(ijk) = temp(ijk) + tdel
      !
  enddo

    call PROF_rapend  ('_SBM_Melting', 2)

  return
  end subroutine melting

  !-----------------------------------------------------------------------------
  subroutine collmain( &
       ijkmax, &
       temp,   &
       ghyd,   &
       dt      )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: ghyd(nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dt                     ! Time step interval
    !---------------------------------------------------------------------------

    if ( rndm_flgp == 1 ) then ! stochastic method
       call r_collcoag( ijkmax,      & ! [IN]
                        wgtbin,      & ! [IN]
                        temp(:),     & ! [IN]
                        ghyd(:,:,:), & ! [INOUT]
                        dt           ) ! [IN]
    else  ! default
       call collcoag( ijkmax,      & ! [IN]
                      temp(:),     & ! [IN]
                      ghyd(:,:,:), & ! [INOUT]
                      dt           ) ! [IN]
    endif

    return
  end subroutine collmain

  !-----------------------------------------------------------------------------
  subroutine collmainf( &
       ijkmax, &
       temp,   &
       ghyd,   &
       dt      )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: ghyd(nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dt                     ! Time step interval
    !---------------------------------------------------------------------------

    if ( rndm_flgp == 1 ) then ! stochastic method
       call r_collcoag( ijkmax,      & ! [IN]
                        wgtbin,      & ! [IN]
                        temp(:),     & ! [IN]
                        ghyd(:,:,:), & ! [INOUT]
                        dt           ) ! [IN]
    else  ! default
       call collcoag( ijkmax,      & ! [IN]
                      temp(:),     & ! [IN]
                      ghyd(:,:,:), & ! [INOUT]
                      dt           ) ! [IN]
    endif

    return
  end subroutine collmainf

  !-----------------------------------------------------------------------------
  !--- reference paper
  !    Bott et al. (1998) J. Atmos. Sci. vol.55, pp. 2284-
  subroutine collcoag( &
       ijkmax, &
       temp,   &
       gc,     &
       dtime   )
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

    integer :: i, j, k, l
    real(RP) :: xi, xj, xnew, dmpi, dmpj, frci, frcj
    real(RP) :: gprime, gprimk, wgt, crn, sum, flux
    integer, parameter :: ldeg = 2
    real(RP) :: acoef( 0:ldeg )
    real(RP), parameter :: dmpmin = 1.E-01_RP
    real(RP) :: suri, surj

    integer :: myu, n, irsl, ilrg, isml
    integer :: ibnd( ijkmax )
    integer :: iflg( nspc,ijkmax )
    integer :: iexst( nbin,nspc,ijkmax )
    real(RP) :: csum( nspc,ijkmax )
    integer :: ijk, nn, mm, pp, qq
    !---------------------------------------------------------------------------

    call PROF_rapstart('_SBM_CollCoag', 2)

  iflg( :,: ) = 0
  iexst( :,:,: ) = 0
  csum( :,: ) = 0.0_RP
  do ijk = 1, ijkmax
     !--- judgement of particle existence
     do myu = 1, nspc
     do n = 1, nbin
       csum( myu,ijk ) = csum( myu,ijk ) + gc( n,myu,ijk )*dxmic
     enddo
     enddo
     do myu = 1, nspc
      if ( csum( myu,ijk ) > cldmin ) iflg( myu,ijk ) = 1
     enddo

     if ( temp(ijk) < tcrit ) then
        ibnd(ijk) = 1
     else
        ibnd(ijk) = 2
     endif

     do myu = 1, nspc
     do n = 1, nbin
       if ( gc( n,myu,ijk ) > cldmin ) then
         iexst( n,myu,ijk ) = 1
       endif
     enddo
     enddo

  enddo

!OCL PARALLEL
  do ijk = 1, ijkmax
    do isml = 1, nspc
    do nn = 1, iflg( isml,ijk )

    do ilrg = 1, nspc
    do mm = 1, iflg( ilrg,ijk )
      !--- rule of interaction
      irsl = ifrsl( ibnd(ijk),isml,ilrg )

      do i = 1, nbin-1  ! small
      do pp  = 1, iexst( i,isml,ijk )

      do j = i+1, nbin  ! large
      do qq  = 1, iexst( j,ilrg,ijk )

        k = kindx( i,j )
        xi = expxctr( i )
        xj = expxctr( j )
        xnew = log( xi+xj )

        dmpi = ck( isml,ilrg,i,j )*gc( j,ilrg,ijk )/xj*dxmic*dtime   ! dg_{i}/dt*dt in page 119 of Suzuki (2004)
        dmpj = ck( ilrg,isml,i,j )*gc( i,isml,ijk )/xi*dxmic*dtime   ! dg_{j}/dt*dt in page 119 of Suzuki (2004)

        if ( dmpi <= dmpmin ) then
          frci = gc( i,isml,ijk )*dmpi                               ! Dg_{i} in page 119 of Suzuki (2004)
        else
          frci = gc( i,isml,ijk )*( 1.0_RP-exp( -dmpi ) )            ! Dg_{i} in page 119 of Suzuki (2004)
        endif

        if ( dmpj <= dmpmin ) then
          frcj = gc( j,ilrg,ijk )*dmpj                               ! Dg_{j} in page 119 of Suzuki (2004)
        else
          frcj = gc( j,ilrg,ijk )*( 1.0_RP-exp( -dmpj ) )            ! Dg_{j} in page 119 of Suzuki (2004)
        endif

        gprime = frci+frcj
!        if ( gprime <= 0.0_RP ) cycle large
        if ( gprime > 0.0_RP .AND. k < nbin ) then

          suri = gc( i,isml,ijk )
          surj = gc( j,ilrg,ijk )
          gc( i,isml,ijk ) = gc( i,isml,ijk )-frci
          gc( j,ilrg,ijk ) = gc( j,ilrg,ijk )-frcj
          gc( i,isml,ijk ) = max( gc( i,isml,ijk )-frci, 0.0_RP )
          gc( j,ilrg,ijk ) = max( gc( j,ilrg,ijk )-frcj, 0.0_RP )
          frci = suri - gc( i,isml,ijk )
          frcj = surj - gc( j,ilrg,ijk )
          gprime = frci+frcj                                       ! g' in page 119 of Suzuki (2004)

          gprimk = gc( k,irsl,ijk ) + gprime                       ! g'_{k} in page 119 of Suzuki (2004)
          wgt = gprime / gprimk                                    ! w in page 119 of Suzuki (2004)
          crn = ( xnew-xctr( k ) )/( xctr( k+1 )-xctr( k ) )       ! c_{k} in page 119 of Suzuki (2004)

          acoef( 0 ) = -( gc( k+1,irsl,ijk )-26.0_RP*gprimk+gc( k-1,irsl,ijk ) )/24.0_RP ! a_{k,0} in page 119 of Suzuki (2004)
          acoef( 1 ) =  ( gc( k+1,irsl,ijk )-gc( k-1,irsl,ijk ) ) *0.50_RP ! a_{k,1} in page 119 of Suzuki (2004)
          acoef( 2 ) =  ( gc( k+1,irsl,ijk )-2.0_RP*gprimk+gc( k-1,irsl,ijk ) ) *0.50_RP ! a_{k,2} in page 119 of Suzuki (2004)

          sum = 0.0_RP
          do l = 0, ldeg
            sum = sum + acoef( l )/( l+1 )/2.0_RP**( l+1 )   &
                      *( 1.0_RP-( 1.0_RP-2.0_RP*crn )**( l+1 ) )
          enddo

          flux = wgt*sum                                           ! f_{k+1/2} in page 119 of Suzuki (2004)
          flux = min( max( flux,0.0_RP ),gprime )

          gc( k,irsl,ijk ) = gprimk - flux                         ! tilda{g_{k}} in page 119 of Suzuki (2004)
          gc( k+1,irsl,ijk ) = gc( k+1,irsl,ijk ) + flux           ! tilda{g_{k+1}} in page 119 of Suzuki (2004)

        endif

      enddo
      enddo !large
      enddo
      enddo !small
    !
    enddo
    enddo

    enddo
    enddo

  enddo

    call PROF_rapend  ('_SBM_CollCoag', 2)

    return
  end subroutine collcoag

  !-------------------------------------------------------------------
  subroutine getrule   &
      ( ifrsl,indx      ) !--- out
  ! subroutine for creating lookup table (Table A.2 of Suzuki 2004)
  integer, intent(out) :: indx( nbin,nbin )
  integer, intent(out) :: ifrsl( 2,nspc_mk,nspc_mk )
  integer :: i, j, k
  real(RP) :: xnew
  !
  do i = 1, nbin
  do j = 1, nbin
    xnew = log( expxctr( i )+expxctr( j ) )
    k = int( ( xnew-xctr( 1 ) )/dxmic ) + 1
    k = max( max( k,j ),i )
    indx( i,j ) = k
  enddo
  enddo
  !
  !
  !--- liquid + liquid -> liquid
  ifrsl( 1:2,il,il ) = il
  !
  !--- liquid + column -> ( graupel, hail ) + column
  ifrsl( 1:2,il,ic ) = ic
  ifrsl( 1,ic,il ) = ig
  ifrsl( 2,ic,il ) = ih
  !
  !--- liquid + plate -> ( graupel, hail ) + plate
  ifrsl( 1:2,il,ip ) = ip
  ifrsl( 1,ip,il ) = ig
  ifrsl( 2,ip,il ) = ih
  !
  !--- liquid + dendrite -> ( graupel, hail ) + dendrite
  ifrsl( 1:2,il,id ) = id
  ifrsl( 1,id,il ) = ig
  ifrsl( 2,id,il ) = ih
  !
  !--- liquid + snowflake -> ( graupel, hail ) + snowflake
  ifrsl( 1:2,il,iss ) = iss
  ifrsl( 1,iss,il ) = ig
  ifrsl( 2,iss,il ) = ih
  !
  !--- liquid + graupel -> ( graupel, hail )
  ifrsl( 1:2,il,ig ) = ig
  ifrsl( 1,ig,il ) = ig
  ifrsl( 2,ig,il ) = ih
  !
  !--- liquid + hail -> ( graupel, hail )
  ifrsl( 1:2,il,ih ) = ih
  ifrsl( 1,ih,il ) = ig
  ifrsl( 2,ih,il ) = ih
  !
  !
  !--- column + column -> snowflake
  ifrsl( 1:2,ic,ic ) = iss
  !
  !--- column + plate -> snowflake
  ifrsl( 1:2,ic,ip ) = iss
  ifrsl( 1:2,ip,ic ) = iss
  !
  !--- column + dendrite -> snowflake
  ifrsl( 1:2,ic,id ) = iss
  ifrsl( 1:2,id,ic ) = iss
  !
  !--- column + snowflake -> snowflake
  ifrsl( 1:2,ic,iss ) = iss
  ifrsl( 1:2,iss,ic ) = iss
  !
  !--- column + graupel -> column + graupel
  ifrsl( 1:2,ic,ig ) = ig
  ifrsl( 1:2,ig,ic ) = ic
  !
  !--- column + hail -> column + ( graupel, hail )
  ifrsl( 1:2,ih,ic ) = ic
  ifrsl( 1,ic,ih ) = ig
  ifrsl( 2,ic,ih ) = ih
  !
  !
  !--- plate + plate -> snowflake
  ifrsl( 1:2,ip,ip ) = iss
  !
  !--- plate + dendrite -> snowflake
  ifrsl( 1:2,ip,id ) = iss
  ifrsl( 1:2,id,ip ) = iss
  !
  !--- plate + snowflake -> snowflake
  ifrsl( 1:2,ip,iss ) = iss
  ifrsl( 1:2,iss,ip ) = iss
  !
  !--- plate + graupel -> plate + graupel
  ifrsl( 1:2,ip,ig ) = ig
  ifrsl( 1:2,ig,ip ) = ip
  !
  !--- plate + hail -> plate + ( graupel, hail )
  ifrsl( 1:2,ih,ip ) = ip
  ifrsl( 1,ip,ih ) = ig
  ifrsl( 2,ip,ih ) = ih
  !
  !
  !--- dendrite + dendrite -> snowflake
  ifrsl( 1:2,id,id ) = iss
  !
  !--- dendrite + snowflake -> snowflake
  ifrsl( 1:2,id,iss ) = iss
  ifrsl( 1:2,iss,id ) = iss
  !
  !--- dendrite + graupel -> dendrite + graupel
  ifrsl( 1:2,id,ig ) = ig
  ifrsl( 1:2,ig,id ) = id
  !
  !--- dendrite + hail -> dendrite + ( graupel, hail )
  ifrsl( 1:2,ih,id ) = id
  ifrsl( 1,id,ih ) = ig
  ifrsl( 2,id,ih ) = ih
  !
  !
  !--- snowflake + snowflake -> snowflake
  ifrsl( 1:2,iss,iss ) = iss
  !
  !--- snowflake + graupel -> snowflake + graupel
  ifrsl( 1:2,iss,ig ) = ig
  ifrsl( 1:2,ig,iss ) = iss
  !
  !--- snowflake + hail -> snowflake + ( graupel, hail )
  ifrsl( 1:2,ih,iss ) = iss
  ifrsl( 1,iss,ih ) = ig
  ifrsl( 2,iss,ih ) = ih
  !
  !
  !--- graupel + graupel -> graupel
  ifrsl( 1:2,ig,ig ) = ig
  !
  !--- graupel + hail -> ( graupel, hail )
  ifrsl( 1,ig,ih ) = ig
  ifrsl( 1,ih,ig ) = ig
  ifrsl( 2,ig,ih ) = ih
  ifrsl( 2,ih,ig ) = ih
  !
  !--- hail + hail -> hail
  ifrsl( 1:2,ih,ih ) = ih
  !
  !
  return
  !
  end subroutine getrule

  !-----------------------------------------------------------------------------
  subroutine faero( &
       ijkmax, &
       f0,     &
       ga      )
    implicit none

    integer , intent(in)    :: ijkmax
    real(RP), intent(in)    :: f0(ijkmax)
    real(RP), intent(inout) :: ga(nccn,ijkmax)

!    real(RP), parameter :: alpha = 3.0_RP

    integer :: ijk, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('_SBM_FAero', 2)

    do ijk = 1, ijkmax
    do n   = 1, nccn
       ga(n,ijk) = ga(n,ijk) + f0(ijk) * marate(n) * expxactr(n) / dxaer
    enddo
    enddo

    call PROF_rapend  ('_SBM_FAero', 2)

    return
  end subroutine faero

  !-----------------------------------------------------------------------------
  ! + Y. Sato added for stochastic method
  ! + Reference Sato et al. (2009) JGR, doi:10.1029/2008JD011247
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
   if ( a < mbin ) then
    write(*,*) "mbin should be smaller than {nbin}_C_{2}"
    call PRC_MPIstop
   endif

   wgtbin = a/real( mbin )
   nbinr = real( nbin )

    do p = 1, pq
      orderb( p ) = p
    enddo

    do p = 1, nbin-1
      ranstmp( (p-1)*nbin-(p*(p-1))/2+1 : p*nbin-(p*(p+1))/2 ) = p
     do q = 1, nbin-p
        ranltmp( (p-1)*nbin-(p*(p-1))/2+q ) = p+q
      enddo
   enddo

    do n = 1, mset
      call RANDOM_get( randnum )
       do p = 1, pq
        abq1 = randnum( 1,1,p )
        k = int( abq1*( pq-p-1 ) ) + p
        temp = orderb( p )
        orderb( p ) = orderb( k )
        orderb( k ) = temp
       enddo

       do p = 1, mbin
        if ( p <= pq ) then
         rans( p ) = ranstmp( orderb( p ) )
         ranl( p ) = ranltmp( orderb( p ) )
        else
         rans( p ) = ranstmp( orderb( p-pq ) )
         ranl( p ) = ranltmp( orderb( p-pq ) )
        endif
         if ( rans( p ) >= ranl( p ) ) then
          tmp1 = rans( p )
          rans( p ) = ranl( p )
          ranl( p ) = tmp1
         endif
       enddo
         blrg( n,1:mbin ) = int( ranl( 1:mbin ) )
         bsml( n,1:mbin ) = int( rans( 1:mbin ) )
    enddo

    deallocate( ranstmp )
    deallocate( ranltmp )
    deallocate( orderb )
    deallocate( randnum )

  end subroutine random_setup

  !-----------------------------------------------------------------------------
  !--- reference paper
  !    Bott et al. (1998) J. Atmos. Sci. vol.55, pp. 2284-
  !    Bott et al. (2000) J. Atmos. Sci. Vol.57, pp. 284-
  subroutine r_collcoag( &
       ijkmax, &
       swgt,   &
       temp,   &
       gc,     &
       dtime   )
    use scale_random, only: &
       RANDOM_get
    implicit none

    integer,  intent(in)    :: ijkmax
    real(RP), intent(in)    :: swgt
    real(RP), intent(in)    :: temp(ijkmax)           ! Temperature       [K]
    real(RP), intent(inout) :: gc  (nbin,nspc,ijkmax) ! Mass size distribution function of hydrometeor
    real(DP), intent(in)    :: dtime                  ! Time step interval

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
!    real(RP) :: beta
    real(RP) :: tmpi, tmpj

    integer :: ibnd( ijkmax )
    integer :: iflg( nspc,ijkmax )
    integer :: iexst( nbin,nspc,ijkmax )
    real(RP) :: csum( nspc,ijkmax )
    integer :: ijk, nn, mm, pp, qq, myu, n, isml, ilrg, irsl
    !---------------------------------------------------------------------------

    call PROF_rapstart('_SBM_CollCoagR', 2)

  iflg( :,: ) = 0
  iexst( :,:,: ) = 0
  csum( :,: ) = 0.0_RP
  do ijk = 1, ijkmax
     !--- judgement of particle existence
     do n = 1, nbin
       csum( il,ijk ) = csum( il,ijk ) + gc( n,il,ijk )*dxmic
       csum( ic,ijk ) = csum( ic,ijk ) + gc( n,ic,ijk )*dxmic
       csum( ip,ijk ) = csum( ip,ijk ) + gc( n,ip,ijk )*dxmic
       csum( id,ijk ) = csum( id,ijk ) + gc( n,id,ijk )*dxmic
       csum( iss,ijk ) = csum( iss,ijk ) + gc( n,iss,ijk )*dxmic
       csum( ig,ijk ) = csum( ig,ijk ) + gc( n,ig,ijk )*dxmic
       csum( ih,ijk ) = csum( ih,ijk ) + gc( n,ih,ijk )*dxmic
     enddo
     if ( csum( il,ijk ) > cldmin ) iflg( il,ijk ) = 1
     if ( csum( ic,ijk ) > cldmin ) iflg( ic,ijk ) = 1
     if ( csum( ip,ijk ) > cldmin ) iflg( ip,ijk ) = 1
     if ( csum( id,ijk ) > cldmin ) iflg( id,ijk ) = 1
     if ( csum( iss,ijk ) > cldmin ) iflg( iss,ijk ) = 1
     if ( csum( ig,ijk ) > cldmin ) iflg( ig,ijk ) = 1
     if ( csum( ih,ijk ) > cldmin ) iflg( ih,ijk ) = 1

     if ( temp(ijk) < tcrit ) then
        ibnd(ijk) = 1
     else
        ibnd(ijk) = 2
     endif

     do myu = 1, nspc
     do n = 1, nbin
       if ( gc( n,myu,ijk ) > cldmin ) then
         iexst( n,myu,ijk ) = 1
       endif
     enddo
     enddo

  enddo

!OCL PARALLEL
  do ijk = 1, ijkmax
    do isml = 1, nspc
    do nn = 1, iflg( isml,ijk )

    do ilrg = 1, nspc
    do mm = 1, iflg( ilrg,ijk )
      !--- rule of interaction
      irsl = ifrsl( ibnd(ijk),isml,ilrg )

      call RANDOM_get( rndm )
      det = int( rndm(1,1,1)*IA*JA*KA )
      nbinr = real( nbin )
      mbinr = real( mbin )
      nums( 1:mbin ) = bsml( det,1:mbin )
      numl( 1:mbin ) = blrg( det,1:mbin )

      do s = 1, mbin
       i = nums( s )
       j = numl( s )

       do pp  = 1, iexst( i,isml,ijk )
       do qq  = 1, iexst( j,ilrg,ijk )

          k = kindx( i,j )
          xi = expxctr( i )
          xj = expxctr( j )
          xnew = log( xi+xj )

          dmpi = ck( isml,ilrg,i,j )*gc( j,ilrg,ijk )/xj*dxmic*dtime
          dmpj = ck( ilrg,isml,i,j )*gc( i,isml,ijk )/xi*dxmic*dtime

          if ( dmpi <= dmpmin ) then
            frci = gc( i,isml,ijk )*dmpi
          else
            frci = gc( i,isml,ijk )*( 1.0_RP-exp( -dmpi ) )
          endif

          if ( dmpj <= dmpmin ) then
            frcj = gc( j,ilrg,ijk )*dmpj
          else
            frcj = gc( j,ilrg,ijk )*( 1.0_RP-exp( -dmpj ) )
          endif
          tmpi = gc( i,isml,ijk )
          tmpj = gc( j,ilrg,ijk )

          gc( i,isml,ijk ) = gc( i,isml,ijk )-frci*swgt
          gc( j,ilrg,ijk ) = gc( j,ilrg,ijk )-frcj*swgt

          if ( j /= k ) then
           gc( j,ilrg,ijk ) = max( gc( j,ilrg,ijk ), 0.0_RP )
          endif
           gc( i,isml,ijk ) = max( gc( i,isml,ijk ), 0.0_RP )

          frci = tmpi - gc( i,isml,ijk )
          frcj = tmpj - gc( j,ilrg,ijk )

          gprime = frci+frcj

          !-----------------------------------------------
          !--- Exponential Flux Method (Bott, 2000, JAS)
          !-----------------------------------------------
      !    if ( gprime <= 0.0_RP ) cycle !large
      !    if ( gprime > 0.0_RP .AND. k < nbin ) then
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
      !    endif

          !-----------------------------------------------
          !--- Flux Method (Bott, 1998, JAS)
          !-----------------------------------------------
!          if ( gprime <= 0.0_RP ) cycle !large
          if ( gprime > 0.0_RP .AND. k < nbin ) then

            gprimk = gc( k,irsl,ijk ) + gprime
            wgt = gprime / gprimk
            crn = ( xnew-xctr( k ) )/( xctr( k+1 )-xctr( k ) )

            acoef( 0 ) = -( gc( k+1,irsl,ijk )-26.0_RP*gprimk+gc( k-1,irsl,ijk ) )/24.0_RP
            acoef( 1 ) =  ( gc( k+1,irsl,ijk )-gc( k-1,irsl,ijk ) ) *0.5_RP
            acoef( 2 ) =  ( gc( k+1,irsl,ijk )-2.0_RP*gprimk+gc( k-1,irsl,ijk ) ) *0.50_RP

            sum = 0.0_RP
            do l = 0, ldeg
              sum = sum + acoef( l )/( l+1 )/2.0_RP**( l+1 )   &
                        *( 1.0_RP-( 1.0_RP-2.0_RP*crn )**( l+1 ) )
            enddo

            flux = wgt*sum
            flux = min( max( flux,0.0_RP ),gprime )

            gc( k,irsl,ijk ) = gprimk - flux
            gc( k+1,irsl,ijk ) = gc( k+1,irsl,ijk ) + flux
          endif

       enddo
       enddo

      enddo ! bin
    !
    enddo
    enddo

    enddo
    enddo

  enddo

    call PROF_rapend  ('_SBM_CollCoagR', 2)

    return
  end subroutine r_collcoag

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_suzuki10_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_grid_index
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
       DENS0, &
       TEMP0  )
    use scale_grid_index
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius          [cm]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)        ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)        ! temperature               [K]

    real(RP), parameter :: um2cm = 100.0_RP

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
                            * rexpxctr( iq-(I_QV+nbin*(ihydro-1)+iq) ) &   !--- mass -> number
                            * radc( iq-(I_QV+nbin*(ihydro-1)+iq) )**3.0_RP )
         sum2(k,i,j,ihydro) = sum2(k,i,j,ihydro) + &
                            ( ( QTRC0(k,i,j,iq) * DENS0(k,i,j) ) & !--- [kg/kg] -> [kg/m3]
                            * rexpxctr( iq-(I_QV+nbin*(ihydro-1)+iq) ) &   !--- mass -> number
                            * radc( iq-(I_QV+nbin*(ihydro-1)+iq) )**2.0_RP )
      enddo
      sum2(k,i,j,ihydro) = 0.5_RP + sign(0.5_RP,sum2(k,i,j,ihydro)-EPS)
      sum3(k,i,j,ihydro) = 0.5_RP + sign(0.5_RP,sum3(k,i,j,ihydro)-EPS)

      if ( sum2(k,i,j,ihydro) /= 0.0_RP ) then
       Re(k,i,j,ihydro) = sum3(k,i,j,ihydro) / sum2(k,i,j,ihydro) * um2cm
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
    use scale_grid_index
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
  subroutine mkpara

  implicit none

  integer :: i, j
  !-----------------------------------------------------------------------------

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
      cctr( n,myu ) = cctr_mk( myu,n )
      write( fid_micpara,* ) myu, n, cctr( n,myu )
    end do
    do n = 1, nbin+1
      cbnd( n,myu ) = cbnd_mk( myu,n )
      write( fid_micpara,* ) myu, n, cbnd( n,myu )
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

end module scale_atmos_phy_mp_suzuki10
!-------------------------------------------------------------------------------
