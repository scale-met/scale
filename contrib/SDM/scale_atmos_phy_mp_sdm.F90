!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Super Droplet Method (SDM)
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-09-30 (Y.Sato)  [new] Implement from Original version of SDM
!! @li      2014-01-22 (Y.Sato)  [rev] Update for scale-0.0.0
!! @li      2014-05-04 (Y.Sato)  [rev] Update for scale-0.0.1
!! @li      2014-06-06 (S.Shima) [rev] Modify several bug 
!! @li      2014-06-07 (Y.Sato)  [rev] Remove dt=max(dt_sdm) and some other changes
!! @li      2014-06-09 (S.Shima) [rev] Check whether dt is the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt. Fixed a hidden bug: "/=" to "==".
!! @li      2014-06-13 (S.Shima) [rev] Common variables are separated into sdm_common.f90
!! @li      2014-06-14 (S.Shima) [rev] Check the initialization of the random number generator.
!! @li      2014-06-24 (S.Shima) [rev] Separated sdm_allocinit
!! @li      2014-06-25 (S.Shima) [rev] Bugfix and improvement of ATMOS_PHY_MP_sdm_restart_in and sdm_iniset
!! @li      2014-06-25 (S.Shima) [rev] Bugfix of sd position initialization, and many modification
!! @li      2014-06-25 (S.Shima) [rev] Bugfix: dx_sdm, dy_sdm, dxiv_sdn, dyiv_sdm restored
!! @li      2014-06-26 (S.Shima) [rev] sdm_getrklu and sdm_z2rk are separated
!! @li      2014-06-27 (S.Shima) [rev] sd data output functionality added
!! @li      2014-07-04 (S.Shima) [rev] Removed comment outputs for debugging
!! @li      2014-07-09 (S.Shima) [rev] Subroutines related to boundary conditions and MPI communications are totally revised
!! @li      2014-07-11 (S.Shima) [rev] Subroutines related to sdm_getvz are revised. Many bug fixes.
!! @li      2014-07-11 (S.Shima) [rev] Subroutines for conversion between fluid variables are separated into module m_sdm_fluidconv
!! @li      2014-07-11 (S.Shima) [rev] Subroutines to impose boundary conditions are separated into module m_sdm_boundary
!! @li      2014-07-11 (S.Shima) [rev] Motion (advection/sedimentation/precipitation) related subroutines are separated into the module m_sdm_motion
!! @li      2014-07-12 (S.Shima) [rev] Add comments concerning when to diagnose QC and QR
!! @li      2014-07-12 (S.Shima) [rev] BUG of random number initialization removed
!<
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
module scale_atmos_phy_mp_sdm
  !-----------------------------------------------------------------------------
  !
  !++ used modules ! For encapsulation, reduce the use of modules here as far as possible. Modules should be called inside the subroutine here.
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer_sdm

  use scale_process, only: &
     mype => PRC_myrank, &
     PRC_MPIstop
  use scale_time, only: &
     TIME_DOATMOS_restart, &
     TIME_NOWSEC
  use gadg_algorithm, only: &
     gadg_count_sort
  use rng_uniform_mt, only: &
     c_rng_uniform_mt, &
     rng_save_state, &
     rng_load_state, &
     gen_rand_array => rng_generate_array
  use m_sdm_common
  use m_sdm_numset
  use m_sdm_memmgr
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  public :: ATMOS_PHY_MP_sdm_setup
  public :: ATMOS_PHY_MP_sdm
  public :: ATMOS_PHY_MP_sdm_CloudFraction
  public :: ATMOS_PHY_MP_sdm_EffectiveRadius
  public :: ATMOS_PHY_MP_sdm_Mixingratio
  public :: ATMOS_PHY_MP_sdm_restart_out
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  logical,  public, save   :: sd_rest_flg_out = .false. ! restart flg of Super Droplet
  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_nmax,    &
       PRC_next,    &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_W,       &
       PRC_E,       &
       PRC_S,       &
       PRC_N,       &
       mype => PRC_myrank, &
       PRC_master
    use scale_const, only: &
       PI     => CONST_PI,    &
       GRAV   => CONST_GRAV,  &
       dens_w => CONST_DWATR,  &
       dens_i => CONST_DICE
    use scale_specfunc, only: &
       SF_gamma
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_grid_real, only: & ! These should be all REAL_xx. Check later.
       REAL_FZ
    use scale_grid, only: &
!!$       CDZ => GRID_CDZ,   &
         GRID_CDX,   &
         GRID_CDY,   &
         GRID_RCDX,   &
         GRID_RCDY,   &
!!$       CZ  => GRID_CZ,    &
!!$       FZ  => GRID_FZ,    &
!!$       FDX => GRID_FDX,   &
!!$       FDY => GRID_FDY,   &
!!$       FDZ => GRID_FDZ,   & 
       GRID_FX,    &
       GRID_FY,    &
       CBFZ => GRID_CBFZ, &
       CBFX => GRID_CBFX, &
       CBFY => GRID_CBFY, &
       ISG, IEG, JSG, JEG,&
       DX,DY,DZ
    use scale_topography, only: &
       TOPO_Zsfc
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE

    logical :: PRC_PERIODIC_X = .true. !< periodic condition or not (X)?
    logical :: PRC_PERIODIC_Y = .true. !< periodic condition or not (Y)?
    namelist / PARAM_PRC / &
         PRC_PERIODIC_X, &
         PRC_PERIODIC_Y

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       docondensation, &
       doautoconversion, &
       doprecipitation

    NAMELIST / PARAM_ATMOS_PHY_MP_SDM / &
       RANDOM_IN_BASENAME, &
       RANDOM_OUT_BASENAME, &
       SD_IN_BASENAME, &
       SD_OUT_BASENAME, &
       domovement, &
       donegative_fixer, &
       sdm_dtcmph, &
       sdm_rdnc, &
       sdm_sdnmlvol, &
       sdm_inisdnc, &
       sdm_aslset, &
       sdm_aslmw, &
       sdm_zlower, &
       sdm_zupper, &
       sdm_extbuf, &
       sdm_rqc2qr, &
       sdm_aslfmrate, &
       sdm_aslfmdt, &
       sdm_aslfmsdnc, &
       sdm_colbrwn, &
       sdm_colkrnl, &
       sdm_mvexchg, &
       sdm_nadjdt, &
       sdm_nadjvar,&
       sdm_dmpvar,&
       sdm_dmpitva,&
       sdm_dmpnskip,& 
       sdm_dmpitvb,& 
       sdm_dmpitvl,&
       sdm_dmpsdsiz

    real(RP) :: dtevl
    real(RP) :: n0, dry_r
    real(RP) :: delta1, delta2 !, sdn_tmp
    real(RP) :: buffact
    integer :: ierr
    integer :: i, j, ip, k, n, s
    integer :: bndsdmdim, bufsiz
    integer :: tmppe, ilcm, igcd
    character(len=17) :: fmt1="(A, '.', A, I*.*)"
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Super Droplet Method'

    if ( MP_TYPE /= 'SDM' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_MP is not SDM. Check!'
       call PRC_MPIstop
    endif

    buffact = 0.0_RP
    do k = KS, KE
      buffact = max( buffact,CBFZ(k) )
    enddo
    do j = JS, JE
      buffact = max( buffact,CBFY(j) )
    enddo
    do i = IS, IE
      buffact = max( buffact,CBFX(i) )
    enddo

    if( buffact > 0.0_RP ) then
       sthopt = 1
    elseif( buffact == 0.0_RP ) then
       sthopt = 0
    endif

    if(sthopt==1) then
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: stretched coordinate is not yet supported!'
       call PRC_MPIstop
    end if

    if( maxval( TOPO_Zsfc ) > 0.0_RP ) then
       trnopt = 2
    elseif( maxval( TOPO_Zsfc ) == 0.0_RP ) then
       trnopt = 0
    endif

    if(trnopt==2) then
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: terrain following coordinate is not yet supported!'
       call PRC_MPIstop
    end if

!! Are these used??
    allocate( KMIN1(KA) )
    allocate( IMIN1(IA) )
    allocate( JMIN1(JA) )
    allocate( KPLS1(KA) )
    allocate( IPLS1(IA) )
    allocate( JPLS1(JA) )

    do k = 2, KA
      KMIN1(k) = k-1
    enddo
    KMIN1(1) = 1
    do i = 2, IA
      IMIN1(i) = i-1
    enddo
    IMIN1(1) = 1
    do j = 2, JA
      JMIN1(j) = j-1
    enddo
    JMIN1(1) = 1

    do k = 1, KA-1
      KPLS1(k) = k+1
    enddo
    KPLS1(KA) = KA
    do i = 1, IA-1
      IPLS1(i) = i+1
    enddo
    IPLS1(IA) = IA
    do j = 1, JA-1
      JPLS1(j) = j+1
    enddo
    JPLS1(JA) = JA

    if( (PRC_PERIODIC_X /= .true.) .or. (PRC_PERIODIC_Y /= .true.))then
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: Only periodic B.C. is supported!'
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: Set PRC_PERIODIC_X=PRC_PERIODIC_Y=.true.'
       call PRC_MPIstop
    else
       wbc=1
       ebc=1
       sbc=1
       nbc=1
    end if

    nsub  = max( PRC_nmax,1 )
    nisub = max( PRC_NUM_X,1 )
    njsub = max( PRC_NUM_Y,1 )

    !--- set process next to the process
    dstw_sub = PRC_next(PRC_W)
    dste_sub = PRC_next(PRC_E)
    srcw_sub = PRC_next(PRC_W)
    srce_sub = PRC_next(PRC_E)

    dsts_sub = PRC_next(PRC_S)
    dstn_sub = PRC_next(PRC_N)
    srcs_sub = PRC_next(PRC_S)
    srcn_sub = PRC_next(PRC_N)

    tag  = 0

! Why overwritten? sdm_zlower and sdm_zupper are defined in run.conf and should not be altered.
!    sdm_zlower = CZ(KS)
!    sdm_zupper = CZ(KE)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SDM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SDM. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP_SDM)

    if( RANDOM_IN_BASENAME == '' ) then
       write(*,*) 'xxx Set Random number file! stop'
       write(*,*) 'To generate random number set, run scale_init with'
       write(*,*) 'flg_sdm = .true. in PARAM_MKINIT of init.conf'
       call PRC_MPIstop
    else
       write(fmt1(14:14),'(I1)') 6
       write(fmt1(16:16),'(I1)') 6
       write(RANDOM_IN_BASENAME,fmt1) trim(RANDOM_IN_BASENAME),'pe',mype
    endif
    fid_random_i = IO_get_available_fid()

    if( SD_IN_BASENAME == '' ) then
       sd_rest_flg_in = .false.
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found S.D. file, generate from random number.'
    else
       sd_rest_flg_in = .true.
       fid_sd_i = IO_get_available_fid()
       write(fmt1(14:14),'(I1)') 6
       write(fmt1(16:16),'(I1)') 6
       write(SD_IN_BASENAME,fmt1) trim(SD_IN_BASENAME),'pe',mype
    endif

    if( SD_OUT_BASENAME == '' ) then
       sd_rest_flg_out = .false.
    else
       sd_rest_flg_out = .true.
    endif

    if( docondensation )   sdm_calvar(1) = .true.
    if( doautoconversion ) sdm_calvar(2) = .true.
    if( domovement )       sdm_calvar(3) = .true.

! zlower+surface height is the lower boundary of SDs.
!!$     if( sdm_zlower < CZ(KS) ) then
!!$      if( mype == PRC_master )  write(*,*) "sdm_zlower was set to CZ(KS) because zlower < CZ(KS)"
!!$      sdm_zlower = CZ(KS)
!!$     endif
! 

!     if( sdm_zupper > CZ(KE) ) then
!      if( mype == PRC_master )  write(*,*) "sdm_zupper was set to CZ(KE) because zupper > CZ(KE)"
!      sdm_zupper = CZ(KE)
!     endif

    if( sdm_zupper > minval(REAL_FZ(KE,IS:IE,JS:JE)) ) then
       if( mype == PRC_master )  write(*,*) "sdm_zupper was set to minval(REAL_FZ(KE)) because zupper > minval(REAL_FZ(KE))"
       sdm_zupper = minval(REAL_FZ(KE,IS:IE,JS:JE))
    endif

     sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
     sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
     sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets

    ! check whether sdm_dtcmph(1:3) > 0
     if(  ( (sdm_dtcmph(1) <= 0.0_RP) .and. docondensation   )   .or. &
         ( (sdm_dtcmph(2) <= 0.0_RP) .and. doautoconversion )   .or. &
         ( (sdm_dtcmph(3) <= 0.0_RP) .and. domovement       )        ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: sdm_dtcmph(1:3) have to be positive'
       call PRC_MPIstop
     end if

    ! check whether dt (dt of mp) is larger than sdm_dtcmph(i).
    if( (dt < sdm_dtcmph(1)) .or. (dt < sdm_dtcmph(2)) .or. (dt < sdm_dtcmph(3)) ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: For now, sdm_dtcmph should be smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_MPIstop
    end if

    ! aerosol nucleation and sd number adjustment functions are not supported yet
    if ( (abs(sdm_aslset) >= 10) .or. (sdm_nadjvar /= 0)) then
       if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: aerosol nucleation and sd number adjustment functions are not supported yet'
       if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: set sdm_aslset < 10 and sdm_nadjvar =0'
       call PRC_MPIstop
    end if

    ! rigorous momentum exchange function is not supported yet
    if ( sdm_mvexchg /= 0) then
       if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: Momentum exchange not yet supported. set sdm_mvexchg = 0'
       call PRC_MPIstop
    end if

    if( docondensation ) then
       nclstp(1)=10*int(1.E+2_RP*(dt+0.0010_RP))            &
            /int(1.E+3_RP*(sdm_dtcmph(1)+0.00010_RP))
       if(mod(10*int(1.E+2_RP*(dt+0.0010_RP)),int(1.E+3_RP*(sdm_dtcmph(1)+0.00010_RP))) /= 0) then 
          if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(1) = 1
    end if

    if( doautoconversion ) then
       nclstp(2)=10*int(1.E+2_RP*(dt+0.0010_RP))            &
            /int(1.E+3_RP*(sdm_dtcmph(2)+0.00010_RP))
       if(mod(10*int(1.E+2_RP*(dt+0.0010_RP)),int(1.E+3_RP*(sdm_dtcmph(2)+0.00010_RP))) /= 0) then 
          if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(2) = 1
    end if

    if( domovement ) then
       nclstp(3)=10*int(1.E+2_RP*(dt+0.0010_RP))            &
            /int(1.E+3_RP*(sdm_dtcmph(3)+0.00010_RP))
       if(mod(10*int(1.E+2_RP*(dt+0.0010_RP)),int(1.E+3_RP*(sdm_dtcmph(3)+0.00010_RP))) /= 0) then 
          if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(3) = 1
    end if
    nclstp(0)=min(nclstp(1),nclstp(2),nclstp(3))

    ilcm=0
    iterate_2: do
       ilcm=ilcm+1
       nclstp(0)=nclstp(0)*ilcm
       if( mod(nclstp(0),nclstp(1)) == 0 .and.           &
           mod(nclstp(0),nclstp(2)) == 0 .and.           &
           mod(nclstp(0),nclstp(3)) == 0   ) then
           exit iterate_2
       else
          nclstp(0)=nclstp(0)/ilcm
       end if
    end do iterate_2

    ! check whether dt is the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt 
    !! find the smallest nclstp(1:3) that satisfies sdm_dtcmph(i)<= dt 
    igcd=maxval(nclstp(1:3))
    if((sdm_dtcmph(1)<= dt).and.docondensation)   igcd=min(nclstp(1),igcd)
    if((sdm_dtcmph(2)<= dt).and.doautoconversion) igcd=min(nclstp(2),igcd)
    if((sdm_dtcmph(3)<= dt).and.domovement)       igcd=min(nclstp(3),igcd)
    !! find the greatest common divisor of nclstp(1:3) that satisfies sdm_dtcmph(i)<= dt 
    do 
       if( igcd == 1) exit
       if( mod(nclstp(1),igcd) == 0 .and.           &
           mod(nclstp(2),igcd) == 0 .and.           &
           mod(nclstp(3),igcd) == 0   ) then
           exit
       end if
       igcd = igcd - 1
    end do
    !! if igcd>1, it meanst dt is not the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt
    if( igcd > 1)then
       if ( IO_L ) write(IO_FID_LOG,*) 'ERROR: TIME_DTSEC_ATMOS_PHY_MP should be the least comon multiple of sdm_dtcmph(1:3) that are smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_MPIstop
    end if

!!$    dx_sdm(1:IA) = CDX(1:IA)
!!$    dy_sdm(1:JA) = CDY(1:JA)
!!$    dz_sdm(1:KA) = CDZ(1:KA)
!!$    dx_sdm(1:IA) = FDX(1:IA)
!!$    dy_sdm(1:JA) = FDY(1:JA)
!!$    dz_sdm(1:KA) = FDZ(1:KA)
    xmax_sdm = GRID_FX(IE)-GRID_FX(IS-1)
    ymax_sdm = GRID_FY(JE)-GRID_FY(JS-1)

    allocate( zph_crs(KA,IA,JA) )

    !--- set number of super droplet etc...
    call sdm_numset(              &
      sdm_extbuf,                 &
      sdm_aslset, sdm_sdnmlvol,   &
      sdm_inisdnc,                &
      sdm_aslfmsdnc,              &
      sdm_zlower, sdm_zupper,     &
      minzph, sdininum_s2c,       &
      sdfmnum_s2c, sdnum_s2c,     &
      ni_s2c, nj_s2c, nk_s2c, zph_crs )

    call sdm_allocinit

    dx_sdm(1:IA) = GRID_CDX(1:IA)
    dy_sdm(1:JA) = GRID_CDY(1:JA)
    dxiv_sdm(1:IA) = GRID_RCDX(1:IA)
    dyiv_sdm(1:JA) = GRID_RCDY(1:JA)

    return
  end subroutine ATMOS_PHY_MP_sdm_setup
  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP !,&
!       TIME_NOWSEC
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    use scale_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_gridtrans, only: &
       I_XYZ, I_XYW,    &
       GTRANS_GSQRT, &
       GTRANS_J13G,  &
       GTRANS_J23G,  &
       GTRANS_J33G
    use scale_const, only: &
       CPdry => CONST_CPdry, &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       P00   => CONST_PRE00
    use scale_atmos_thermodyn, only: &
       CPw => AQ_CP

    use m_sdm_io, only: &
       sdm_outasci
    use m_sdm_coordtrans, only: &
       sdm_rk2z
    use m_sdm_fluidconv, only: &
       sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod
    use m_sdm_motion, only: &
       sdm_getvz
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)        !! Momentum [kg/s/m2]
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]

    real(RP):: dtcl(1:5) ! 1-3 : Time interval of cloud micro physics at current processed time
                         ! 4   : Time interval to form super droplets as aerosol
                         ! 5   : Time interval to adjust number of super droplets in each grid
!    real(RP) :: rst_crs(KA,IA,JA)  ! Base state density x Jacobian
    real(RP) :: pbr_crs(KA,IA,JA)  ! Base state pressure
    real(RP) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
    real(RP) :: ppf_crs(KA,IA,JA)  ! Pressure perturbation at future
!!$    real(RP) :: uf_crs(KA,IA,JA)   ! u components of velocity at future
!!$    real(RP) :: vf_crs(KA,IA,JA)   ! v components of velocity at future
!!$    real(RP) :: wf_crs(KA,IA,JA)   ! w components of velocity at future
!!$!    real(RP) :: wcf_crs(KA,IA,JA)  ! zeta components of contravariant velocity at future
    real(RP) :: ptpf_crs(KA,IA,JA) ! Potential temperature perturbation at future
    real(RP) :: qvf_crs(KA,IA,JA)  ! Water vapor mixing ratio at future
!!$    real(RP) :: prr_crs(IA,JA,1:2) ! Precipitation and accumulation for rain
    real(RP) :: rtmp4(KA,IA,JA)    ! Temporary array
    real(RP) :: rtmp5(KA,IA,JA)    ! Temporary array
    real(RP) :: rtmp6(KA,IA,JA)    ! Temporary array
    ! Output variables
    real(RP) :: exnr_crs(KA,IA,JA) ! Exner function
    integer :: bufsiz1      ! Buffer size for MPI
    integer :: bufsiz2      ! Buffer size for MPI
    ! Work variables
    logical :: lsdmup       ! flag for updating water hydrometeor by SDM
    real(RP) :: sd_nc  ! averaged number concentration in a grid

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)
    real(RP) :: QDRY(KA,IA,JA)
    real(RP) :: CPTOT(KA,IA,JA)
    real(RP) :: RTOT(KA,IA,JA)
    real(RP) :: CPovCV(KA,IA,JA)
    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: flux_tot (KA,IA,JA)
    real(RP) :: flux_rain(KA,IA,JA)
    real(RP) :: flux_snow(KA,IA,JA)
    real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
    real(RP) :: crs_dtmp3(KA,IA,JA), crs_dtmp4(KA,IA,JA)
    real(RP) :: crs_dtmp5(KA,IA,JA), crs_dtmp6(KA,IA,JA)
    integer  :: n, s, k, i, j, iq         ! index
    logical, save ::  TIME_DO_RANDOM_restart = .false.

    real(RP)::tmp_pres, tmp_temp

    real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
    real(RP) :: rhod_scale(KA,IA,JA) ! dry air density
    real(RP) :: t_scale(KA,IA,JA)    ! Temperature

    !---------------------------------------------------------------------------

    ! QTRC except QV (and QDRY) are diagnosed from super-droplets
    ! To make this doubly sure, reset QTRC to zero
    do iq = QQS, QQE
       if(iq /= I_QV) then
          QTRC(:,:,:,iq) = 0.0_RP
       end if
    end do

    if( sd_first ) then
      sd_first = .false.
      if( IO_L ) write(IO_FID_LOG,*) '*** S.D.: setup'
      if( sd_rest_flg_in ) then
         !---- read restart file
         call ATMOS_PHY_MP_sdm_restart_in
      else
         !---- set initial condition
         call sdm_iniset(DENS, RHOT, QTRC,                   &
                      RANDOM_IN_BASENAME, fid_random_i,   &
                      xmax_sdm, ymax_sdm, sdm_dtcmph,     &
                      sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
                      sdm_inisdnc,sdm_zlower,             &
                      sdm_zupper,sdm_calvar,              &
!                      dx_sdm,dy_sdm,dz_sdm,               &
!                      nqw,jcb,                            &
!                      qwtr_crs,zph_crs,                   &
!                      jcb,                                &
                      zph_crs,                            &
                      sdasl_s2c, sdx_s2c, sdy_s2c,        &
                      sdz_s2c, sdr_s2c,                   &
                      sdrk_s2c, sdvz_s2c,                 &
                      sdrkl_s2c, sdrku_s2c                )
         prr_crs(1:IA,1:JA,1:2)=0.0_RP
      end if
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(Super Droplet Method)'

    ! -----
    ! Initialize
    if( .not. sdm_calvar(1) .and. &
       .not. sdm_calvar(2) .and. &
        .not. sdm_calvar(3)       ) return

    dtcl(1:3) = sdm_dtcmph(1:3)
    dtcl(4) =  sdm_aslfmdt
    dtcl(5) =  sdm_nadjdt

    ! Evaluate the diagnostic fluid variables needed for SDM from SCALE intrinsic variables
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Legacy code. Temporarily used for now but finally removed
!    mf(:,:) = 1.0_RP
    QDRY(:,:,:) = 1.0_RP
    do k = 1, KA
    do i = 1, IA
    do j = 1, JA
!      rst_crs(k,i,j) = DENS(k,i,j) * jcb(k,i,j)
      ptbr_crs(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
      do iq = QQS, QQE
        QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
      enddo
      RTOT (k,i,j) = Rdry * QDRY(k,i,j) + Rvap * QTRC(k,i,j,I_QV)
      CPTOT(k,i,j) = CPdry * QDRY(k,i,j)
      do iq = QQS, QQE
        CPTOT(k,i,j) = CPTOT(k,i,j) + QTRC(k,i,j,iq) * CPw(iq)
      enddo
      CPovCV(k,i,j) = CPTOT(k,i,j) / ( CPTOT(k,i,j) - RTOT(k,i,j) )
      pbr_crs(k,i,j) = P00 * ( RHOT(k,i,j) * RTOT(k,i,j) / P00 )**CPovCV(k,i,j)
      ppf_crs(k,i,j) = 0.0_RP
      ptpf_crs(k,i,j) = 0.0_RP
      qvf_crs(k,i,j) = QTRC(k,i,j,I_QV)
      QTRC(k,i,j,I_QC:I_QR) = 0.0_RP
    enddo
    enddo
    enddo
!     jcb8w(KA,:,:) = GTRANS_GSQRT(KA-1,:,:,I_XYW)
!     jcb(KA,:,:) = GTRANS_GSQRT(KA-1,:,:,I_XYW)
!!$    do k = 1, KA
!!$    do i = 1, IA
!!$    do j = 1, JA
!!$      uf_crs(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i,j)+DENS(k,IMIN1(i),j) )
!!$      vf_crs(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j)+DENS(k,i,JMIN1(j)) )
!!$      wf_crs(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k,i,j)+DENS(KMIN1(k),i,j) )
!!$!      wcf_crs(k,i,j) = 2.0_RP * MOMZ(k,i,j) * jcb8w(k,i,j) / ( DENS(k,i,j)+DENS(KMIN1(k),i,j) )
!!$    enddo
!!$    enddo
!!$    enddo
    ! -----
    ! Perform super-droplets method (SDM)
    !== get the exner function at future ==!
    call getexner(pbr_crs,ppf_crs,exnr_crs,rtmp4)
    !== get the zeta components of contravariant velocity ==!
    !== at future                                         ==!
!     call phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,    &
!                  ni,nj,nk,j31,j32,jcb8w,mf,                      &
!                  uf_crs,vf_crs,wf_crs,wcf_crs,rtmp4,rtmp5,rtmp6)

    ! Get dry air density at the start of SDM.
    call sdm_getrhod(pbr_crs,ptbr_crs,ppf_crs,ptpf_crs,       &
                       qvf_crs,rhod_crs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Refactored code here. The legacy code above are now being refactored to this part.
    !!! SCALE intrinsic variables 
    !! DENS:               Density [kg/m3]
    !! MOMZ,MOMX,MOMY:     Momentum [kg/s/m2]
    !! RHOT:               DENS * POTT [K*kg/m3]
    !! QTRC(KA,IA,JA,QAD): Ratio of mass of tracer to total mass[kg/kg]

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! after removing OMP support, this should be also omitted.
    ! Finally, bufsiz1 = bufsiz and bufsiz2 = bndsdmdim, and these should be determined once in sdm_numset or somewhere and should be stored as common variables
    bufsiz1 = nint( sdininum_s2c*(real(sdm_extbuf)*1.E-2_RP) )
    bufsiz1 = nomp * ( int((bufsiz1-1)/nomp) + 1 ) !! being multiple number to 'nomp'
    bufsiz2 = 8 + sdnumasl_s2c    !! n,x,y,rk,u,v,wc(vz),r,asl


    !--
    ! SD data output
    if( (mod(sdm_dmpvar,10)==1) .and. sdm_dmpitva>0.0_RP .and. &
         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
             int(1.E+3_RP*(sdm_dmpitva+0.00010_RP))) == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) ' *** Output Super Droplet Data in ASCII'
       call sdm_rk2z(sdnum_s2c,sdx_s2c,sdy_s2c,sdrk_s2c,sdz_s2c,sdri_s2c,sdrj_s2c)
!!$       ! for testing sdm_getvz
!!$       ppf_crs(:,:,:) = 0.0_RP; ptpf_crs(:,:,:) = 0.0_RP
!!$       tmp_pres=50000.0_RP; tmp_temp=263.15_RP
!!$       pbr_crs(:,:,:) = tmp_pres; rhod_crs(:,:,:) = tmp_pres/tmp_temp/Rdry; ptbr_crs(:,:,:) = tmp_temp*(P00/tmp_pres)**(Rdry/CPdry)
       call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
       call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz(pres_scale,rhod_scale,t_scale,            &
                           sdnum_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c,sdvz_s2c,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'motion' )
       call sdm_outasci(TIME_NOWSEC,                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdm_dmpnskip)
    else if(sdm_dmpvar>1)then
       if( IO_L ) write(IO_FID_LOG,*) 'ERROR: sdm_dmpvar>1 not supported for now. Set sdm_dmpvar=1 (Output Super Droplet Data in ASCII)'
    end if

    !== run SDM at future ==!
     call sdm_calc(MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,                 & 
                   sdm_calvar,sdm_mvexchg,dtcl(1:3), sdm_aslset,  &
                   exnr_crs,pbr_crs,ptbr_crs,         &
                   ppf_crs,ptpf_crs,qvf_crs,&
                   prr_crs,zph_crs,rhod_crs,                      &
                   lsdmup,ni_s2c,nj_s2c,nk_s2c,                   &
                   sdnum_s2c,sdnumasl_s2c,                        &
                   sdn_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,              &
                   sdu_s2c,sdv_s2c,sdvz_s2c,sdr_s2c,sdasl_s2c,    &
                   sdrkl_s2c,sdrku_s2c,                           &
                   rng_s2c,rand_s2c,sortid_s2c,sortkey_s2c,       &
!                   rand_s2c,sortid_s2c,sortkey_s2c,               &
                   sortfreq_s2c,sorttag_s2c,                      &
                   bufsiz1,bufsiz2,sdm_itmp1,sdm_itmp2,           &
                   sd_itmp1,sd_itmp2,sd_itmp3,sd_dtmp1,           &
                   crs_dtmp1,crs_dtmp2,crs_dtmp3,crs_dtmp4,       &
                   crs_dtmp5,crs_dtmp6,rbuf,sbuf)

     !== convert updated contravariant velocity of ==!
     !== super-droplets to {u,v,w} at future       ==!

!     if( sdm_mvexchg>0 ) then

!        call cnt2phy(idsthopt,idtrnopt,idmpopt,idmfcopt,          &
!                     ni,nj,nk,j31,j32,jcb8w,mf,                   &
!                     uf_crs,vf_crs,wcf_crs,wf_crs,                &
!                     rtmp4,rtmp5,rtmp6)
!     end if

     !== convert super-droplets to {qc,qr} at future ==!

!     call sdm_sd2qcqr(nqw,pbr_crs,ptbr_crs,                         &
!                      ppf_crs,ptpf_crs,qvf_crs,qwtrf_crs,zph_crs,   &
     ! We need to rethink whether it's okay to evaluate QC QR here
     ! For example, when we diagnose pressure or rhod during the dynamical process, 
     ! qc and qr should be 0 because they are not included in the total density
     call sdm_sd2qcqr(pbr_crs,ptbr_crs,                             &
                      ppf_crs,ptpf_crs,qvf_crs,zph_crs,             &
                      lsdmup,sdnum_s2c,sdn_s2c,sdx_s2c,sdy_s2c,     &
                      sdr_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,         &
                      rhod_crs,rhoc_sdm,rhor_sdm,                   &
                      !--- Y.Sato add ---
                      rhoa_sdm, sdnumasl_s2c, sdasl_s2c,            &
                      !--- Y.Sato add ---
                      sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2)


! not supported yet. S.Shima
!!$     ! Aerosol formation process of super-droplets
!!$
!!$     if( dtcl(4)>0.0_RP .and. &
!!$         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
!!$             int(1.E+3_RP*(dtcl(4)+0.00010_RP))) == 0 ) then
!!$
!!$        call sdm_aslform(DENS,RHOT,QTRC,                             &   
!!$                         sdm_calvar,sdm_aslset,                      &
!!$                         sdm_aslfmsdnc,sdm_sdnmlvol,                 &
!!$                         sdm_zupper,sdm_zlower,dtcl,                 &
!!$                         jcb,pbr_crs,ptbr_crs,ppf_crs,               &
!!$                         ptpf_crs,qvf_crs,zph_crs,rhod_crs,          &
!!$                         sdnum_s2c,sdnumasl_s2c,sdn_s2c,sdx_s2c,     &
!!$                         sdy_s2c,sdz_s2c,sdrk_s2c,sdu_s2c,           &
!!$                         sdv_s2c,sdvz_s2c,sdr_s2c,sdasl_s2c,         &
!!$                         sdfmnum_s2c,sdn_fm,sdx_fm,sdy_fm,sdz_fm,    &
!!$                         sdri_fm,sdrj_fm,sdrk_fm,sdvz_fm,sdr_fm,sdasl_fm,            &
!!$                         ni_s2c,nj_s2c,nk_s2c,                       &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,        &
!!$                         sorttag_s2c,rng_s2c,                        &
!!$!                         sorttag_s2c,                                &
!!$                         sdm_itmp1,sd_itmp1,sd_itmp2,sd_itmp3)
!!$
!!$     end if

! not supported yet. S.Shima
!!$     ! Adjust number of super-droplets
!!$
!!$     if( dtcl(5)>0.0_RP .and. &
!!$         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
!!$             int(1.E+3_RP*(dtcl(5)+0.00010_RP))) == 0 ) then
!!$
!!$       !== averaged number concentration in a grid ==!
!!$
!!$       sd_nc = sdininum_s2c/real(ni_s2c*nj_s2c*knum_sdm,kind=RP)
!!$
!!$       call sdm_adjsdnum(sdm_nadjvar,ni_s2c,nj_s2c,nk_s2c,          &
!!$                         sdnum_s2c,sdnumasl_s2c,sd_nc,              &
!!$                         sdn_s2c,sdx_s2c,sdy_s2c,sdr_s2c,           &
!!$                         sdasl_s2c,sdvz_s2c,sdrk_s2c,               &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,       &
!!$                         sorttag_s2c,rng_s2c,rand_s2c,                &
!!$!                         sorttag_s2c,rand_s2c,                      &
!!$                         sdm_itmp1,sdm_itmp2,sd_itmp1)
!!$
!!$      end if

! No feedback to atmosphere for debugging.
!!$    !--- update MOMENTUM, RHOT, and QV
!!$    do k = 1, KA
!!$    do i = 1, IA
!!$    do j = 1, JA
!!$      RHOT(k,i,j) = ( ptbr_crs(k,i,j)+ptpf_crs(k,i,j) ) * DENS(k,i,j)
!!$!      MOMX(k,i,j) = uf_crs(k,i,j) * 0.5_RP * ( DENS(k,i,j)+DENS(k,IPLS1(i),j) )
!!$!      MOMY(k,i,j) = vf_crs(k,i,j) * 0.5_RP * ( DENS(k,i,j)+DENS(k,i,JPLS1(j)) )
!!$!      MOMZ(k,i,j) = wcf_crs(k,i,j) * 0.5_RP * ( DENS(k,i,j)+DENS(KPLS1(k),i,j) ) / jcb8w(k,i,j)
!!$      QTRC(k,i,j,I_QV) = qvf_crs(k,i,j)
!!$      QTRC(k,i,j,I_QC) = rhoc_sdm(k,i,j) / DENS(k,i,j)
!!$      QTRC(k,i,j,I_QR) = rhor_sdm(k,i,j) / DENS(k,i,j)
!!$    enddo
!!$    enddo
!!$    enddo

!! Are these correct? Let's check later.
    call HIST_in( rhoa_sdm(:,:,:), 'RAERO', 'aerosol mass conc.', 'kg/m3', dt)
    call HIST_in( prr_crs(:,:,1), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( prr_crs(:,:,1), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt)

    return
  end subroutine ATMOS_PHY_MP_sdm
  !-----------------------------------------------------------------------------
   subroutine sdm_iniset(DENS, RHOT, QTRC,                   &
                         RANDOM_IN_BASENAME, fid_random_i,   &
                         xmax_sdm, ymax_sdm, dtcmph,         &
                         sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
                         sdm_inisdnc,sdm_zlower,             &
                         sdm_zupper,sdm_calvar,              &
!                         dx_sdm,dy_sdm,dz_sdm,               &
!                         nqw,jcb,                            &
!                         qwtr_crs,zph_crs,                   &
!                         jcb,                                &
                         zph_crs,                            &
                         sdasl_s2c, sdx_s2c, sdy_s2c,        &
                         sdz_s2c, sdr_s2c,                   &
                         sdrk_s2c, sdvz_s2c,                 &
                         sdrkl_s2c, sdrku_s2c                )
  !***********************************************************************
  ! Input variables
      use scale_const, only: &
        P00 => CONST_PRE00, &
        Rdry => CONST_Rdry, &
        Rvap => CONST_Rvap, &
        CPdry => CONST_CPdry
      use scale_process, only: &
        PRC_MPIstop, &
        mype => PRC_myrank
      use scale_atmos_thermodyn, only: &
        CPw => AQ_CP
      use scale_tracer, only: &
        QAD => QA
      use scale_grid, only: &
        GRID_FX,    &
        GRID_FY
      use m_sdm_coordtrans, only: &
        sdm_getrklu, &
        sdm_z2rk
      use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod
      use m_sdm_motion, only: &
        sdm_getvz

      real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
      real(RP), intent(in) :: RHOT(KA,IA,JA) ! DENS * POTT [K*kg/m3]
      real(RP), intent(in) :: QTRC(KA,IA,JA,QAD) ! ratio of mass of tracer to total mass[kg/kg]
      character(len=H_LONG), intent(in) :: RANDOM_IN_BASENAME
      integer, intent(in) :: fid_random_i
      real(RP),intent(in) :: xmax_sdm, ymax_sdm
      real(RP),intent(in) :: dtcmph(3)  ! Time interval of cloud micro physics
      logical, intent(in) :: sdm_calvar(3)! Flag for cond./coll/move calculation
      real(RP),intent(in) :: sdm_rdnc     ! Number concentration of real droplets
      real(RP),intent(in) :: sdm_sdnmlvol ! Normal volume for number concentration of super droplets
      integer, intent(in) :: sdm_aslset   ! Option for aerosol species
      real(RP),intent(in) :: sdm_inisdnc  ! Initial number of super droplets per sdm_sdnmlvol
      real(RP),intent(in) :: sdm_zlower   ! Lower limitaion of initial SDs position
      real(RP),intent(in) :: sdm_zupper   ! Upper limitaion of initial SDs position
!      integer, intent(in) :: nqw       ! Number of water hydrometeor array
!      real(RP),intent(in) :: jcb(KA,IA,JA)      ! Jacobian at scalar points
      real(RP),intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinates
!      real(RP),intent(in) :: dx_sdm ! Grid distance in x direction for SDM
!      real(RP),intent(in) :: dy_sdm ! Grid distance in y direction for SDM
!      real(RP),intent(in) :: dz_sdm ! Grid distance in z direction for SDM
!      real(RP),intent(inout) :: qwtr_crs(KA,IA,JA,nqw) ! Water hydrometeor (past)
      real(RP),intent(inout) :: sdasl_s2c(1:sdnum_s2c,1:sdnumasl_s2c)
      real(RP),intent(inout) :: sdx_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdy_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdz_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdr_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdrk_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdvz_s2c(1:sdnum_s2c)
      real(RP),intent(inout) :: sdrkl_s2c(IA,JA)
      real(RP),intent(inout) :: sdrku_s2c(IA,JA)
      ! Work variables
!      real(RP) :: rbr_crs(KA,IA,JA)   ! density ! not used
      real(RP) :: ptbr_crs(KA,IA,JA)  ! potential temperature
      real(RP) :: ptp_crs(KA,IA,JA)   ! Potential temperature perturbation
      real(RP) :: pbr_crs(KA,IA,JA)   ! pressure
      real(RP) :: pp_crs(KA,IA,JA)    ! Pressure perturbation
      real(RP) :: qv_crs(KA,IA,JA)    ! Water vapor mixing ratio
      real(RP) :: n0                            ! number of real droplets per unit volume and per aerosol radius
      real(RP) :: dry_r                         ! aerosol radius
      real(RP) :: delta1, delta2, sdn_tmp       ! temporary
      logical :: lsdmup                         ! flag for updating water hydrometeor by SDM
      integer :: iexced, sdnum_tmp1, sdnum_tmp2 ! temporary
      integer :: i, j, k, n, iq, np             ! index
      real(RP) :: CPTOT(KA,IA,JA), RTOT(KA,IA,JA)
      real(RP) :: QDRY(KA,IA,JA),  CPovCV(KA,IA,JA)
      real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
      integer :: sd_str, sd_end, sd_valid

      real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
      real(RP) :: rhod_scale(KA,IA,JA) ! dry air density
      real(RP) :: t_scale(KA,IA,JA)    ! Temperature
     !---------------------------------------------------------------------

      ! conversion of SCALE variables to CReSS variables: ptbr ptp pbr pp rhod qv 
      ! This part will be omitted in the future.
      CPTOT(:,:,:) = 0.0_RP
      QDRY(:,:,:) = 1.0_RP
      do k = 1, KA
      do i = 1, IA
      do j = 1, JA
        ptbr_crs(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
        do iq = QQS, QQE
          QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
        enddo
        rhod_crs(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) ! dry air density
        RTOT (k,i,j) = Rdry * QDRY(k,i,j) + Rvap * QTRC(k,i,j,I_QV)
        CPTOT(k,i,j) = CPdry * QDRY(k,i,j)
        do iq = QQS, QQE
         CPTOT(k,i,j) = CPTOT(k,i,j) + QTRC(k,i,j,iq) * CPw(iq)
        enddo
        CPovCV(k,i,j) = CPTOT(k,i,j) / ( CPTOT(k,i,j) - RTOT(k,i,j) )
        pbr_crs(k,i,j) = P00 * ( RHOT(k,i,j) * RTOT(k,i,j) / P00 )**CPovCV(k,i,j)
        pp_crs(k,i,j) = 0.0_RP
        ptp_crs(k,i,j) = 0.0_RP
        qv_crs(k,i,j) = QTRC(k,i,j,I_QV) ! Be careful. The definition of qv is different. (SCALE: qv=rhov/rho, CReSS: qv=rhov/rhod)
      enddo
      enddo
      enddo

      if( .not. sdm_calvar(1) .and. &
          .not. sdm_calvar(2) .and. &
          .not. sdm_calvar(3)        ) return

      !### Get index[k/real] at "sdm_zlower" and "sdm_zupper"  ###!
      call sdm_getrklu(sdm_zlower,sdm_zupper,      &
                       sdrkl_s2c,sdrku_s2c)

!!$      !### Get base state density ###!
!!$!      do k=1,nk
!!$!      do j=0,nj+1
!!$!      do i=0,ni+1
!!$      do k = 1, KA
!!$      do j = 1, JA
!!$      do i = 1, IA
!!$         rhod_crs(k,i,j) = rbr_crs(k,i,j)
!!$      end do
!!$      end do
!!$      end do

      ! This should not be called here
!!$      ! restart files are used for S.D.
!!$      if( sd_rest_flg_in ) then
!!$         call ATMOS_PHY_MP_sdm_restart_in
!!$         return
!!$      endif

      !### Get random generator seed ###!
      !! Random number generator has already been initialized in scale-les/src/preprocess/mod_mkinit.f90
      !! Be careful. If unit (=fid_random_i) is specified, filename is ignored and the object is initialized by the unit.
      call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME))
!      call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME), fid_random_i )

      ! Initialized super-droplets.
      !### Get parameter for 3mode log-nomiral distribution ###!

      if( mod(sdm_aslset,10)==1 .or. mod(sdm_aslset,10)==3 ) then

         !! Derksen(2009)
         if( IO_L ) then
          write(IO_FID_LOG,*)  "    number concentration of (NH4)2SO4    : Derksen"
         endif
         n1_amsul   = n1_amsul_derksn
         n2_amsul   = n2_amsul_derksn
         rb1_amsul  = rb1_amsul_derksn
         rb2_amsul  = rb2_amsul_derksn
         sgm1_amsul = sgm1_amsul_derksn
         sgm2_amsul = sgm2_amsul_derksn
         rmax_amsul = rmax_amsul_derksn
         rmin_amsul = rmin_amsul_derksn

         if( mod(sdm_aslset,10)==1 ) then
            n3_nacl    = 1.E-10_RP     !! temporary initialize
            rb3_nacl   = 1.E-10_RP
            sgm3_nacl  = 1.E-10_RP
            rmax_nacl  = 1.E-10_RP
            rmin_nacl  = 1.E-10_RP
         end if

      else if( mod(sdm_aslset,10)==-1 .or. mod(sdm_aslset,10)==-3 ) then

         !! vanZanten(2010)
         if( IO_L ) then
          write(IO_FID_LOG,*) "    number concentration of (NH4)2SO4    : vanZanten"
         endif

         n1_amsul   = n1_amsul_zanten
         n2_amsul   = n2_amsul_zanten
         rb1_amsul  = rb1_amsul_zanten
         rb2_amsul  = rb2_amsul_zanten
         sgm1_amsul = sgm1_amsul_zanten
         sgm2_amsul = sgm2_amsul_zanten
         rmax_amsul = rmax_amsul_zanten
         rmin_amsul = rmin_amsul_zanten

         if( mod(sdm_aslset,10)==-1 ) then
            n3_nacl    = 1.E-10_RP     !! temporary initialize
            rb3_nacl   = 1.E-10_RP
            sgm3_nacl  = 1.E-10_RP
            rmax_nacl  = 1.E-10_RP
            rmin_nacl  = 1.E-10_RP
         end if

      end if

      if( mod(sdm_aslset,10)==2 .or. abs(mod(sdm_aslset,10))==3 ) then

         !! Derksen(2009)
         if( IO_L ) then
          write(IO_FID_LOG,*) "number concentration of NaCl aerosol : Derksen"
         endif
         if( mod(sdm_aslset,10)==2 ) then
            n1_amsul   = 1.E-10_RP     !! temporary initialize
            n2_amsul   = 1.E-10_RP
            rb1_amsul  = 1.E-10_RP
            rb2_amsul  = 1.E-10_RP
            sgm1_amsul = 1.E-10_RP
            sgm2_amsul = 1.E-10_RP
            rmax_amsul = 1.E-10_RP
            rmin_amsul = 1.E-10_RP
         end if
         n3_nacl   = n3_nacl_derksn
         rb3_nacl  = rb3_nacl_derksn
         sgm3_nacl = sgm3_nacl_derksn
         rmax_nacl = rmax_nacl_derksn
         rmin_nacl = rmin_nacl_derksn

      end if

      !### Get random number ###!
      call gen_rand_array( rng_s2c, sdx_s2c  )
      call gen_rand_array( rng_s2c, sdy_s2c  )
      call gen_rand_array( rng_s2c, sdz_s2c  )
      call gen_rand_array( rng_s2c, sdr_s2c  )
      call gen_rand_array( rng_s2c, sdvz_s2c )
! Why are they halved??
!!$      do n=1,sdnum_s2c
!!$       sdx_s2c(n) = sdx_s2c(n)-0.5_RP
!!$       sdy_s2c(n) = sdy_s2c(n)-0.5_RP
!!$       sdz_s2c(n) = sdz_s2c(n)-0.5_RP
!!$       sdr_s2c(n) = sdr_s2c(n)-0.5_RP
!!$       sdvz_s2c(n) = sdvz_s2c(n)-0.5_RP
!!$      enddo

      iexced = 1     !! check for memory size of int*8

      do k=1,sdnumasl_s2c
         call gen_rand_array( rng_s2c, sd_dtmp1 )
         do n=1,sdnum_s2c
            sdasl_s2c(n,k) = sd_dtmp1(n)
         end do
      end do

      !### Aerosol mass, muliplicity ###!
      do k=1,sdnumasl_s2c
         do n=1,sdnum_s2c
            i = mod(n-1,sdnumasl_s2c) + 1    !! select aerosol index
            if( abs(sdm_aslset)==12 ) then
               i = 2                         !! 1 : (NH4)2SO4, 2: NaCl
            end if
            if( k==i ) then
               !! match aerosol type
               if( abs(mod(sdm_aslset,10))==1 .or.                      &
                     ( k==1 .and. abs(mod(sdm_aslset,10))==3 ) ) then

                  !### (NH4)2SO4 [g] ###!
                  delta1 = log(rmax_amsul) - log(rmin_amsul)
                  dry_r  = exp( log(rmin_amsul)+delta1*sdasl_s2c(n,k) )
                  sdasl_s2c(n,k) = F_THRD * ONE_PI                 &
                                * (dry_r*dry_r*dry_r) * rho_amsul
                  !! n0(log(dry_r)) [m-3]
                  !! 2-mode log-noraml distribution for ammonium sulfate
                  delta1 = log(dry_r) - log(rb1_amsul)
                  delta1 = -(delta1*delta1)                             &
                               /(2.0_RP*log(sgm1_amsul)*log(sgm1_amsul))
                  delta2 = log(dry_r) - log(rb2_amsul)
                  delta2 = -(delta2*delta2)                             &
                               /(2.0_RP*log(sgm2_amsul)*log(sgm2_amsul))
                  n0 = (n1_amsul*exp(delta1))                           &
                                 /(sqrt(2.0_RP*ONE_PI)*log(sgm1_amsul))  &
                     + (n2_amsul*exp(delta2))                           &
                                 /(sqrt(2.0_RP*ONE_PI)*log(sgm2_amsul))
                  !! number per unit volume and per aerosol species
                  delta1 = real(sdm_inisdnc,kind=RP)                    &
                                    /real(sdm_sdnmlvol,kind=RP)
                  delta1 = delta1/real(sdnumasl_s2c,kind=RP)
                  !! continuous uniform distribution
                  delta2 = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))
                  !! muliplicity
                  sdn_tmp = n0/(delta1*delta2)
               else if( abs(mod(sdm_aslset,10))==2 .or.                 &
                     ( k==2 .and. abs(mod(sdm_aslset,10))==3 ) ) then
                  !### NaCl(seasalt,[g]) ###!
                  delta1 = log(rmax_nacl) - log(rmin_nacl)
                  dry_r  = exp( log(rmin_nacl) + delta1*sdasl_s2c(n,k) )
                  sdasl_s2c(n,k) = F_THRD * ONE_PI                 &
                                * (dry_r*dry_r*dry_r) * rho_nacl
                  !! n0(log(dry_r)) [m-3]
                  !! 1-mode log-noraml distribution for seasalt
                  delta1 = log(dry_r) - log(rb3_nacl)
                  delta1 = -(delta1*delta1)                             &
                               /(2.0_RP*log(sgm3_nacl)*log(sgm3_nacl))
                  n0 = (n3_nacl*exp(delta1))                            &
                                 /(sqrt(2.0_RP*ONE_PI)*log(sgm3_nacl))
                  !! number per unit volume and per aerosol species
                  delta1 = real(sdm_inisdnc,kind=RP)                    &
                                    /real(sdm_sdnmlvol,kind=RP)
                  delta1 = delta1/real(sdnumasl_s2c,kind=RP)
                  !! continuous uniform distribution
                  delta2 = 1.0_RP/(log(rmax_nacl)-log(rmin_nacl))
                  !! muliplicity
                  sdn_tmp = n0/(delta1*delta2)
               else
                  !### Other ###!
                  sdasl_s2c(n,k) = 1.0E-18_RP                           &
                                 + 1.0E-14_RP                           &
                                 * (log(1.0_RP/(1.0_RP-sdasl_s2c(n,k))))
                  sdn_tmp = real(sdm_rdnc,kind=RP)                      &
                          * real(sdm_sdnmlvol,kind=RP)                  &
                          / real(sdm_inisdnc,kind=RP)
               end if
               !! check muliplicity
               if( sdn_tmp<(2.0_RP**63.0_RP) ) then
                  sdn_s2c(n) = nint( sdn_tmp, kind=DP )
               else
                  iexced = -1
               end if
            else
               sdasl_s2c(n,k) = 0.0_RP
            end if
         end do
      end do

      if( iexced<0 ) then
         if( IO_L ) write(IO_FID_LOG,*) "sdm_iniset, exceeded"
         call PRC_MPIstop
      end if

      !### position of super-droplets in horizontal ###!
      do n=1,sdnum_s2c
         sdx_s2c(n) = xmax_sdm * sdx_s2c(n)+GRID_FX(IS-1)
         sdy_s2c(n) = ymax_sdm * sdy_s2c(n)+GRID_FY(JS-1)
      end do

      !### position of super-droplets in vertical ###!
      !! valid super-droplets
      do n=1,nint(sdininum_s2c)
         if( sdn_s2c(n)>0 ) then
            sdz_s2c(n) = real(minzph+sdm_zlower,kind=RP)             &
                 + sdz_s2c(n)                                  &
                 * real(sdm_zupper-(minzph+sdm_zlower),kind=RP)
         else
            sdz_s2c(n) = INVALID     !!! check muliplicity
         end if
      end do

      !! invalid super-droplets
      do n=nint(sdininum_s2c)+1,sdnum_s2c
         sdz_s2c(n) = INVALID
      end do
      
!!$      sdnum_tmp1 = int( nint(sdininum_s2c)/nomp )
!!$      sdnum_tmp2 = mod( nint(sdininum_s2c),nomp )
!!$
!!$      do np=1,nomp
!!$
!!$         sd_str = int(sdnum_s2c/nomp)*(np-1) + 1
!!$         sd_end = int(sdnum_s2c/nomp)*np
!!$
!!$         if( np<=sdnum_tmp2 ) then
!!$            sd_valid = sd_str + sdnum_tmp1
!!$         else
!!$            sd_valid = sd_str + (sdnum_tmp1-1)
!!$         end if
!!$
!!$         !! valid super-droplets
!!$
!!$         do n=sd_str,sd_valid
!!$            if( sdn_s2c(n)>0 ) then
!!$               sdz_s2c(n) = real(minzph+sdm_zlower,kind=RP)             &
!!$                          + sdz_s2c(n)                                  &
!!$                          * real(sdm_zupper-(minzph+sdm_zlower),kind=RP)
!!$            else
!!$               sdz_s2c(n) = INVALID     !!! check muliplicity
!!$            end if
!!$         end do
!!$
!!$         !! invalid super-droplets
!!$         do n=sd_valid+1,sd_end
!!$            sdz_s2c(n) = INVALID
!!$         end do
!!$      end do

      do n=1,sdnum_s2c
!ORG     sdr_s2c(n) = 1.0e-5 * ( log(1.0/(1.0-sdr_s2c(n))) )**O_THRD
!ORG     sdr_s2c(n) = 1.0d-8

! temporary for test
!!$         sdr_s2c(n) = 1.0E-15_RP
         sdr_s2c(n) = 3.0E-3_RP*sdr_s2c(n)
!         sdr_s2c(n) = exp((log(3.0E-3_RP)-log(1.0E-7_RP))*sdr_s2c(n)+log(1.0E-7_RP))
         
      end do

      !### index[k/real] of super-droplets               ###!
      !### modify position[z] of invalid super-droplets  ###!
      call sdm_z2rk(sdm_zlower,sdm_zupper,            &
                        sdnum_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c )

      !### terminal velocity                            ###!
      !### ( at only condensation/evaporation process ) ###!
      if( sdm_calvar(1) ) then

         call sdm_condevp(sdm_aslset,                                   &
                          sdm_aslmw,sdm_aslion,sdm_dtevl,               &
                          pbr_crs,ptbr_crs,pp_crs,ptp_crs,              &
                          qv_crs,                                       &
                          sdnum_s2c,sdnumasl_s2c,sdx_s2c,sdy_s2c,       &
                          sdr_s2c,sdasl_s2c,sdrk_s2c)

      end if

      call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
      call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
      call sdm_getvz(pres_scale,rhod_scale,t_scale,                    &
                     sdnum_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c, &
                     sdvz_s2c,sd_itmp1,sd_itmp2,sd_itmp3,'motion')
      ! -----
      ! Convert super-droplets to mixing ratio of water hydrometeor
      ! in each grid

      lsdmup = .true.

!      call sdm_sd2qcqr(nqw,pbr_crs,ptbr_crs,                            &
!                       pp_crs,ptp_crs,qv_crs,qwtr_crs,zph_crs,          &
      ! We need to rethink whether it's okay to evaluate QC QR here
      ! For example, when we diagnose pressure or rhod during the dynamical process, 
      ! qc and qr should be 0 because they are not included in the total density
      call sdm_sd2qcqr(pbr_crs,ptbr_crs,                                &
                       pp_crs,ptp_crs,qv_crs,zph_crs,                   &
                       lsdmup,sdnum_s2c,sdn_s2c,sdx_s2c,sdy_s2c,        &
                       sdr_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
                       rhod_crs,rhoc_sdm,rhor_sdm,                      &
                       !--- Y.Sato add ---
                       rhoa_sdm,  sdnumasl_s2c, sdasl_s2c,              &
                       !--- Y.Sato add ---
                       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2            )

!      do k = 1, KA
!      do i = 1, IA
!      do j = 1, JA
!         QTRC(k,i,j,I_QC) = rhoc_sdm(k,i,j) / DENS(k,i,j)
!         QTRC(k,i,j,I_QR) = rhor_sdm(k,i,j) / DENS(k,i,j)
!      enddo
!      enddo
!      enddo

      ! Output logfile about SDM
      if( mype==0 ) then

         if( IO_L ) then
          write(IO_FID_LOG,*)
          write(IO_FID_LOG,'(a)')"  ### [SDM] : Information of super-droplets ###"
          write(IO_FID_LOG,*)
         endif
!          write(moji,'(i25)')sdnum_s2c
!          write(*,'(a43,a25)')                                           &
!      &     "    allocate size for super-droplets     : ", adjustl(moji)

!          write(moji,'(i25)')nint(sdininum_s2c)
!          write(IO_FID_LOG,'(a43,a25)')                                 &
!      &     "    initial number of super-droplets     : ", adjustl(moji)

      end if

    return

  end subroutine sdm_iniset
  !-----------------------------------------------------------------------------
  subroutine sdm_condevp(sdm_aslset,                     &
                         sdm_aslmw,sdm_aslion,sdm_dtevl, &
                         pbr_crs,ptbr_crs,pp_crs,        &
                         ptp_crs,qv_crs,                 &
                         sd_num,sd_numasl,sd_x,sd_y,     &
                         sd_r,sd_asl, sd_rk              )

      use scale_const, only: &
         cp   => CONST_CPdry, &
         p0   => CONST_PRE00, &
         t0   => CONST_TEM00, &
         es0  => CONST_PSAT0, &
         rd   => CONST_Rdry
      use scale_grid, only: &
         FX => GRID_FX, &
         FY => GRID_FY
      ! Input variables
      integer,  intent(in) :: sdm_aslset
      real(RP), intent(in) :: sdm_aslmw(20)
      real(RP), intent(in) :: sdm_aslion(20)
      real(RP), intent(in) :: sdm_dtevl  ! tims step of {condensation/evaporation} process
      real(RP), intent(in) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
      real(RP), intent(in) :: ptp_crs(KA,IA,JA)  ! Potential temperature perturbation
      real(RP), intent(in) :: pbr_crs(KA,IA,JA)  ! Base state pressure
      real(RP), intent(in) :: pp_crs(KA,IA,JA)   ! Pressure perturbation
      real(RP), intent(in) :: qv_crs(KA,IA,JA)   ! Water vapor mixing ratio
      integer,  intent(in) :: sd_num      ! number of super-droplets
      integer,  intent(in) :: sd_numasl   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
      real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
      real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
      ! Input and output variables
      real(RP), intent(inout) :: sd_r(1:sd_num) ! equivalent radius of super-droplets

      ! Internal shared variables
      real(RP):: sd_aslmw(1:22) ! Molecular mass of chemical material contained as water-soluble aerosol in super droplets (default+20)
      real(RP):: sd_aslion(1:22) ! Degree of ion dissociation of chemical material contained as water-soluble aerosol in super droplets (default+20)
!      real :: dx_sdm             ! dx_sdm
!      real :: dy_sdm             ! dy_sdm

      ! Work variables
      real(RP) :: dmask(1:22)  ! mask for vactorization
      real(RP) :: p_sd      ! pressure of the grid contained the SD
      real(RP) :: t_sd      ! temperature of the grid contained the SD
      real(RP) :: pt_sd     ! potential temperature of the grid contained the SD.
      real(RP) :: qv_sd     ! water-vapor of the grid contained the SD
      real(RP) :: qvs_sd    ! Saturation mixing ratio of the grid contained the SD.
      real(RP) :: es_sd     ! Saturation vapor pressure of the grid contained the SD.
      real(RP) :: ss_sd     ! Degree of super-saturation of the grid contained the SD.
      real(RP) :: Fac_dd    ! Fd in growth EQ.(R.R.Rogers)
      real(RP) :: Fac_kk    ! Fk in growth EQ.(R.R.Rogers)
      real(RP) :: Rc        ! Rc
      real(RP) :: ivt_sd    ! 1.d0 / t_sd
      real(RP) :: dtivfdk   ! dt / ( Fk + Td )
      real(RP) :: crd       ! sd_r before interation
      real(RP) :: rdi       ! sd_r after interation
      real(RP) :: a         ! temporary
      real(RP) :: a3        ! temporary
      real(RP) :: b         ! temporary
      real(RP) :: eq_a      ! temporary for growth EQ.
      real(RP) :: eq_b      ! temporary for growth EQ.
      real(RP) :: eq_c      ! temporary for growth EQ.
      real(RP) :: new_rd    ! temporary for interation
      real(RP) :: rd2       ! temporary for interation
      real(RP) :: rd5       ! temporary for interation
      real(RP) :: rd7       ! temporary for interation
      real(RP) :: dterm1    ! temporary for interation
      real(RP) :: dterm2    ! temporary for interation
      real(RP) :: dterm3    ! temporary for interation
      real(RP) :: dtmp      ! temporary for interation
      real(RP) :: rddvcp    ! rd / cp
      integer :: idx_nasl(1:22)  ! index for vactorization

      integer :: i, j, k, n, s, t, it               ! index
      integer :: ix, jy
      ! Parameters
      integer, parameter :: itr_max = 25   ! iteration number
      real(RP), parameter :: epsva = 0.622_RP   ! Molecular weight ratio of vapor/air
      real(RP) :: tmpd
     !---------------------------------------------------------------------

      rddvcp = real(rd,kind=RP)/real(cp,kind=RP)

      !### aerosol type ###!

      if( abs(mod(sdm_aslset,10))==1 ) then

         !### numasl=1 @ init+rest : (NH4)2SO4 ###!

         sd_aslmw(1)  = mass_amsul
         sd_aslion(1) = ion_amsul

      else if( abs(mod(sdm_aslset,10))==2 ) then

         if( abs(sdm_aslset)==2 ) then

            !### numasl=1 @ init : NaCl ###!

            sd_aslmw(1)  = mass_nacl
            sd_aslion(1) = ion_nacl

         else if( abs(sdm_aslset)==12 ) then

            !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

            sd_aslmw(1) = mass_amsul
            sd_aslmw(2) = mass_nacl
            sd_aslion(1) = ion_amsul
            sd_aslion(2) = ion_nacl

         end if

      else if( abs(mod(sdm_aslset,10))==3 ) then

         !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

         sd_aslmw(1) = mass_amsul
         sd_aslmw(2) = mass_nacl

         sd_aslion(1) = ion_amsul
         sd_aslion(2) = ion_nacl

!         do n=1,20
!            call getrname( id_sdm_aslmw  + (n-1), sd_aslmw(n+2)  )
!            call getrname( id_sdm_aslion + (n-1), sd_aslion(n+2) )
!         end do

      end if

      do n=1,22

         if( n<=sd_numasl ) then
            idx_nasl(n) = n
            dmask(n) = 1.0_RP
         else
            idx_nasl(n) = sd_numasl
            dmask(n) = 0.0_RP
         end if

      end do


     ! Convert Super-Droplets to density of cloud-water.
      do n=1,sd_num

         !### Skip invalid super-droplets ###!

         if( sd_rk(n)<VALID2INVALID ) cycle

         !### Get the location and variables of Super-Droplets ###!
         iloop : do ix = IS, IE
          if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
           i = ix
           exit iloop
          endif
         enddo iloop
         jloop : do jy = JS, JE
          if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
           j = jy
           exit jloop
          endif
         enddo jloop
!         i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!         j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
         k = floor( sd_rk(n) )

         p_sd  = real( pbr_crs(k,i,j) + pp_crs(k,i,j), kind=RP )
         pt_sd = real( ptbr_crs(k,i,j) + ptp_crs(k,i,j), kind=RP )
         qv_sd = real( qv_crs(k,i,j), kind=RP )

         !### Calculate degree of super-saturation ###!

         t_sd   = pt_sd * exp( rddvcp*log(p_sd/real(p0,kind=RP)) )
         ivt_sd = 1.0_RP / t_sd

         a = 1.0_RP / ( t_sd - 35.860_RP )
         b = a * ( t_sd - real(t0,kind=RP) )

         es_sd  = real(es0,kind=RP) * exp(17.269_RP*b)   !! Teten-eq.
         qvs_sd = epsva * es_sd / ( p_sd - es_sd )
         ss_sd  = qv_sd / qvs_sd          !! degree of super-saturation

         !### Set parameters for growth EQ. of the radius ###!

         !  LatGas = LatHet / GasV_C
         !  L_RL_K = LatHet * DNS_RL / Heat_C
         !  RLRv_D = DNS_RL * GasV_C / Diff_C
         !  ASL_RR = ASL_FF * ION_asl / WGTasl

         Fac_dd  = RLRv_D * t_sd / es_sd
         Fac_kk  = ( LatGas*ivt_sd - 1.0_RP ) * L_RL_K * ivt_sd
         dtivFdk = real( sdm_dtevl, kind=RP ) / ( Fac_dd + Fac_kk )

         eq_a  = CurveF * ivt_sd

         eq_b = 0.0_RP

         do t=1,22

            s = idx_nasl(t)

            dtmp = sd_asl(n,s) * (real(sd_aslion(s),kind=RP)            &
                                       / real(sd_aslmw(s),kind=RP))
            eq_b = eq_b + dmask(t) * dtmp

         end do

         eq_b = eq_b * ASL_FF

         Rc    = sqrt(eq_b/eq_a)
         eq_c  = ss_sd - 1.0_RP

         a  = eq_a / eq_c
         b  = eq_b / eq_c
         a3 = a * a * a

         !### advance the particle radius ###!

         crd    = sd_r(n)
         new_rd = sd_r(n)

         if( ( ss_sd > 1.0_RP ) .and. ( a3 < b*(27.0_RP/4.0_RP) ) ) then
            new_rd = 1.0E-3_RP
         end if

         if( new_rd<Rc ) then
            new_rd = Rc
         end if

         !== iteration ( newton-raphson method ) ==!

         rdi = new_rd

         do it=1,itr_max

            ! Considering Kohler-curve (R.R.Rogers)

            rd2 = rdi * rdi
            rd5 = rd2 * rd2 * rdi
            rd7 = rd5 * rd2

            dterm1 = dtivFdk * ( rd5*eq_c + (eq_b-eq_a*rd2)*rd2 )
            dterm2 = rd5 * crd * crd
            dterm3 = rd5 - dtivFdk*( -3.0_RP*eq_b + eq_a*rd2 )

            dtmp = rd2 - ( rd7 - 2.0_RP*dterm1 - dterm2 ) / dterm3

            if( dtmp<=0.0_RP ) dtmp = rdi * rdi * 1.E-4_RP

            rdi = sqrt(dtmp)

         end do

         sd_r(n) = rdi   !! particle radius at future

      end do

    return
  end subroutine sdm_condevp
  !-----------------------------------------------------------------------------
!  subroutine sdm_sd2qcqr(nqw,pbr_crs,ptbr_crs,                    &
!                         pp_crs,ptp_crs,qv_crs,qwtr_crs,zph_crs,  &
  subroutine sdm_sd2qcqr(pbr_crs,ptbr_crs,                        &
                         pp_crs,ptp_crs,qv_crs,zph_crs,           &
                         lsdmup,sd_num,sd_n,sd_x,sd_y,sd_r,sd_rk, &
                         sd_rkl,sd_rku,                           &
                         rhod_crs,rhoc_sdm,rhor_sdm,              &
                         !--- Y.Sato add ---
                         rhoa_sdm, sd_numasl, sd_nasl,            &
                         !--- Y.Sato add ---
                         sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2    )
      ! Input variables
!      integer, intent(in) :: nqw ! Number of water hydrometeor array
      real(RP), intent(in) :: pbr_crs(KA,IA,JA)  ! Base state pressure
      real(RP), intent(in) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
      real(RP), intent(in) :: pp_crs(KA,IA,JA)   ! Pressure perturbation
      real(RP), intent(in) :: ptp_crs(KA,IA,JA)  ! Potential temperature perturbation
      real(RP), intent(in) :: qv_crs(KA,IA,JA)   ! Water vapor mixing ratio
      real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
      logical, intent(in) :: lsdmup    ! flag for updating water hydrometeor by SDM
      integer, intent(in) :: sd_num    ! number of super-droplets
      integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
      real(RP), intent(in) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
      real(RP), intent(in) :: sd_r(1:sd_num)      ! equivalent radius of super-droplets
      real(RP), intent(in) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
      real(RP), intent(in) :: sd_rkl(IA,JA)  ! lower boundary index[k/real] at scalar point in SDM calculation area
      real(RP), intent(in) :: sd_rku(IA,JA)  ! upper boundary index[k/real] at scalar point in SDM calculation area
      !--- Y.Sato add ---
      integer, intent(in) :: sd_numasl    ! number of super-droplets
      real(RP), intent(in) :: sd_nasl(1:sd_num,1:sd_numasl)     ! index[k/real] of super-droplets
      !--- Y.Sato add ---
      ! Input and output variables
      real(RP), intent(inout) :: rhod_crs(KA,IA,JA)   ! dry air density
      real(RP), intent(inout) :: rhoc_sdm(KA,IA,JA)   ! density of cloud water
      real(RP), intent(inout) :: rhor_sdm(KA,IA,JA)   ! density of rain water
      !--- Y.Sato add ---
      real(RP), intent(inout) :: rhoa_sdm(KA,IA,JA)   ! density of rain water
      !--- Y.Sato add ---
      ! Output variables
!      real(RP), intent(out) :: qwtr_crs(KA,IA,JA,1:nqw)        ! water hydrometeor
      real(RP), intent(out) :: crs_dtmp1(KA,IA,JA)    ! temporary buffer of CReSS dimension
      real(RP), intent(out) :: crs_dtmp2(KA,IA,JA)    ! temporary buffer of CReSS dimension
      integer, intent(out) :: sd_itmp1(1:sd_num,1:nomp)    ! temporary array of the size of the number of super-droplets.
      integer, intent(out) :: sd_itmp2(1:sd_num,1:nomp)    ! temporary array of the size of the number of super-droplets.
      ! Work variables
      integer :: n, i, j, k   ! index
      !-------------------------------------------------------------------

      ! Get dry air density

      ! bug hides in this subroutine
      ! not used for now
      return

      call sdm_getrhod(pbr_crs,ptbr_crs,pp_crs,ptp_crs,        &
                       qv_crs,rhod_crs)

      ! Get density of water hydrometeor

      if( lsdmup ) then

         !### the case updating water hydrometeor by SDM ###!
         !! convert super-droplets to density of cloud water
         !! and rain water.

         call sdm_sd2rhocr(sdm_rqc2qr,                          &
                           zph_crs,rhoc_sdm,rhor_sdm,           &
                           sd_num,sd_n,sd_x,sd_y,sd_r,sd_rk,    &
                           sd_rkl,sd_rku,                       &
                           !--- Y.Sato add ---
                           sd_numasl,sd_nasl,rhoa_sdm,          &
                           !--- Y.Sato add ---
                           crs_dtmp1,crs_dtmp2,sd_itmp1,sd_itmp2)

      end if

      ! Convert water hydrometeor density to mixing ratio

!      do k=2,nk-1
!      do j=1,nj-1
!      do i=1,ni-1
!      do k=KS,KE
!      do j=JS-1,JE+1
!      do i=IS-1,IE+1
!         qwtr_crs(k,i,j,I_QC) = real( rhoc_sdm(k,i,j)/rhod_crs(k,i,j) )
!         qwtr_crs(k,i,j,I_QR) = real( rhor_sdm(k,i,j)/rhod_crs(k,i,j) )
!      end do
!      end do
!      end do

    return
  end subroutine sdm_sd2qcqr
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
   subroutine sdm_sd2rhocr(sdm_rqc2qr,                        &
                           zph_crs,rhoc_sdm,rhor_sdm,         &
                           sd_num, sd_n,sd_x,sd_y,sd_r,sd_rk, &
                           sd_rkl,sd_rku,                     &
                           !--- Y.Sato add ---
                           sd_numasl,sd_nasl,rhoa_sdm,        &
                           !--- Y.Sato add ---
                           liqc_sdm,liqr_sdm,ilist_c,ilist_r)

     use scale_grid, only: &
      FX => GRID_FX, &
      FY => GRID_FY
     use scale_const, only: &
      rw => CONST_DWATR

      ! Input variables
      real(RP),intent(in) :: sdm_rqc2qr          ! Threshould between qc and qr [m]
      real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
      integer, intent(in) :: sd_num              ! Number of super-droplets
      integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
      real(RP), intent(in) :: sd_r(1:sd_num)     ! equivalent radius of super-droplets
      real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
      real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
      real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
      !--- Y.Sato add ---
      integer, intent(in) :: sd_numasl           ! Number of aeorosol species
      real(RP), intent(in):: sd_nasl(1:sd_num,1:sd_numasl)  ! multiplicity of super-droplets
      !--- Y.Sato add ---
      ! Output variables
      real(RP), intent(out) :: rhoc_sdm(KA,IA,JA)  ! densitiy of cloud water
      real(RP), intent(out) :: rhor_sdm(KA,IA,JA)  ! densitiy of rain water
      real(RP), intent(out) :: liqc_sdm(KA,IA,JA)  ! liquid cloud water of super-droplets
      real(RP), intent(out) :: liqr_sdm(KA,IA,JA)  ! liquid rain water of super-droplets
      !--- Y.Sato add ---
      real(RP), intent(out) :: rhoa_sdm(KA,IA,JA)  ! densitiy of
      !--- Y.Sato add ---
      integer, intent(out) :: ilist_c(1:int(sd_num/nomp),1:nomp)  ! buffer for list vectorization
      integer, intent(out) :: ilist_r(1:int(sd_num/nomp),1:nomp)  ! buffer for list vectorization
      ! Work variables for OpenMP
      integer :: sd_str             ! index of divided loop by OpenMP
      integer :: sd_end             ! index of divided loop by OpenMP
      integer :: np                 ! index for OpenMP
      ! Work variables
      real(RP) :: dcoef(KA,IA,JA)   ! coef.
      real(RP) :: drate             ! temporary
      integer :: nlist_c(1:nomp)    ! list number
      integer :: nlist_r(1:nomp)    ! list number
      integer :: nlist_a(1:nomp)    ! list number
      real(RP):: liqa_sdm(KA,IA,JA)  ! liquid rain water of super-droplets
      integer :: ilist_a(1:int(sd_num/nomp),1:nomp)  ! buffer for list vectorization
      integer :: tlist_c            ! total list number for cloud
      integer :: tlist_r            ! total list number for rain
      integer :: tlist_a            ! total list number for aerosol
      integer :: ccnt               ! counter
      integer :: rcnt               ! counter
      integer :: acnt               ! counter
      integer :: i, ix              ! index
      integer :: j, jy              ! index
      integer :: k                  ! index
      integer :: kl                 ! index
      integer :: ku                 ! index
      integer :: m, mm              ! index
      integer :: n                  ! index

!-----7--------------------------------------------------------------7--

     ! Initialize
      do k = 1, KA
      do i = 1, IA
      do j = 1, JA
       dcoef(k,i,j) = rw * F_THRD * ONE_PI / real(dx_sdm(i)*dy_sdm(j),kind=RP)
!       idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!       idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )
      enddo
      enddo
      enddo

      tlist_c = 0
      tlist_r = 0
      tlist_a = 0

      do np=1,nomp
         nlist_c(np) = 0
         nlist_r(np) = 0
         nlist_a(np) = 0
      end do

!      do k=1,nk
!      do j=0,nj+1
!      do i=0,ni+1
      do k=2,KA
      do j=1,JA
      do i=1,IA
         liqc_sdm(k,i,j) = 0.0_RP
         liqr_sdm(k,i,j) = 0.0_RP
         liqa_sdm(k,i,j) = 0.0_RP
         rhoc_sdm(k,i,j) = 0.0_RP
         rhor_sdm(k,i,j) = 0.0_RP
         rhoa_sdm(k,i,j) = 0.0_RP
      end do
      end do
      end do
! -----

! Get index list for compressing buffer.
      do np=1,nomp

         sd_str = int(sd_num/nomp)*(np-1) + 1
         sd_end = int(sd_num/nomp)*np

         ccnt = 0
         rcnt = 0
         acnt = 0

         do n=sd_str,sd_end

            if( sd_rk(n)<VALID2INVALID ) cycle

            if( sd_r(n)<real(sdm_rqc2qr,kind=RP) ) then

               !### cloud-water ###!

               ccnt = ccnt + 1
               ilist_c(ccnt,np) = n
               acnt = acnt + 1
               ilist_a(acnt,np) = n

            else

               !### rain-water ###!
               rcnt = rcnt + 1
               ilist_r(rcnt,np) = n
               acnt = acnt + 1
               ilist_a(acnt,np) = n

            end if

         end do

         nlist_c(np) = ccnt
         nlist_r(np) = rcnt
         nlist_a(np) = acnt

         tlist_c = ccnt
         tlist_r = rcnt
         tlist_a = acnt

      end do

     ! Get density of cloud-water and rain-water.

      !### cloud-water ###!

      if( tlist_c>0 ) then

         do np=1,nomp

            if( nlist_c(np)>0 ) then

               do m=1,nlist_c(np)

                  n = ilist_c(m,np)

!                  i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!                  j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
                  do ix = IS, IE
                   if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                    i = ix
                    exit
                   endif
                  enddo
                  do jy = JS, JE
                   if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                    j = jy
                    exit
                   endif
                  enddo
                  k = floor( sd_rk(n) )

                  liqc_sdm(k,i,j) = liqc_sdm(k,i,j)                     &
                                  + sd_r(n) * sd_r(n) * sd_r(n)         &
                                  * real(sd_n(n),kind=RP)
                  liqa_sdm(k,i,j) = liqa_sdm(k,i,j)                     &
                                    + real(sd_n(n),kind=RP)
               end do

            end if

         end do

         !=== adjust cloud-water in verical boundary. ===!

!         do j=2,nj-2
!         do i=2,ni-2
         do j=JS, JE
         do i=IS, IE

            !! at lower boundary

            kl    = floor(sd_rkl(i,j))
            drate = real(kl+1,kind=RP) - sd_rkl(i,j)
            if( drate<0.50_RP ) then
               liqc_sdm(kl,i,j) = 0.0_RP           !! <50% in share
            else
               liqc_sdm(kl,i,j) = liqc_sdm(kl,i,j)/drate
            end if

            !! at upper boundary

            ku    = floor(sd_rku(i,j))
            drate = sd_rku(i,j) - real(ku,kind=RP)

            if( drate<0.50_RP ) then
               liqc_sdm(ku,i,j) = 0.0_RP           !! <50% in share
            else
               liqc_sdm(ku,i,j) = liqc_sdm(ku,i,j)/drate
            end if

         end do
         end do

         !=== convert super-droplets to density of cloud-water. ===!

!         do k=2,nk-1
!         do j=1,nj-1
!         do i=1,ni-1
         do k=KS, KE+1
         do j=JS-1, JE+1
         do i=IS-1, IE+1
            rhoc_sdm(k,i,j) = liqc_sdm(k,i,j) * dcoef(k,i,j)       &
                            / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)
!            rhoa_sdm(k,i,j) = liqa_sdm(k,i,j) * dcoef(k,i,j)       &
!                            / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)

         end do
         end do
         end do

      end if


      !### rain-water ###!

      if( tlist_r>0 ) then

         do np=1,nomp

            if( nlist_r(np)>0 ) then

               do m=1,nlist_r(np)

                  n = ilist_r(m,np)

!                  i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!                  j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
                  do ix = IS, IE
                   if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                    i = ix
                    exit
                   endif
                  enddo
                  do jy = JS, JE
                   if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                    j = jy
                    exit
                   endif
                  enddo
                  k = floor( sd_rk(n) )

                  liqr_sdm(k,i,j) = liqr_sdm(k,i,j)                     &
                                  + sd_r(n) * sd_r(n) * sd_r(n)         &
                                  * real(sd_n(n),kind=RP)
               end do

            end if

         end do

         !=== adjust rain-water in verical boundary. ===!

!         do j=2,nj-2
!         do i=2,ni-2
         do j=JS, JE
         do i=IS, IE

            !! at lower boundary

            kl    = floor(sd_rkl(i,j))
            drate = real(kl+1,kind=RP) - sd_rkl(i,j)

            if( drate<0.5d0 ) then
               liqr_sdm(kl,i,j) = 0.e0           !! <50% in share
            else
               liqr_sdm(kl,i,j) = liqr_sdm(kl,i,j)/drate
            end if

            !! at upper boundary

            ku    = floor(sd_rku(i,j))
            drate = sd_rku(i,j) - real(ku,kind=RP)

            if( drate<0.5d0 ) then
               liqr_sdm(ku,i,j) = 0.e0           !! <50% in share
            else
               liqr_sdm(ku,i,j) = liqr_sdm(ku,i,j)/drate
            end if

         end do
         end do

         !=== convert super-droplets to density of rain-water. ===!

!         do k=2,nk-1
!         do j=1,nj-1
!         do i=1,ni-1
         do k=KS,KE+1
         do j=JS-1,JE+1
         do i=IS-1,IE+1

            rhor_sdm(k,i,j) = liqr_sdm(k,i,j) * dcoef(k,i,j)      &
                    / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)

         end do
         end do
         end do


      end if

      !### in-cloud-aerosol ###!

      if( tlist_a>0 ) then

         do np=1,nomp

            if( nlist_a(np)>0 ) then

               do m=1,nlist_a(np)

                  n = ilist_a(m,np)

!                  i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!                  j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
                  do ix = IS, IE
                   if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                    i = ix
                    exit
                   endif
                  enddo
                  do jy = JS, JE
                   if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                    j = jy
                    exit
                   endif
                  enddo
                  k = floor( sd_rk(n) )

                  do mm = 1, sd_numasl
                    liqa_sdm(k,i,j) = liqa_sdm(k,i,j)                     &
                                    + real(sd_nasl(n,mm),kind=RP)
                  enddo
!                    liqa_sdm(k,i,j) = liqa_sdm(k,i,j)                     &
!                                    + real(sd_n(n),kind=RP)
               end do

            end if

         end do


        !=== adjust in-cloud-aeorosl in verical boundary. ===!

!         do j=2,nj-2
!         do i=2,ni-2
         do j=JS, JE
         do i=IS, IE

            !! at lower boundary

            kl    = floor(sd_rkl(i,j))
            drate = real(kl+1,kind=RP) - sd_rkl(i,j)

            if( drate<0.50_RP ) then
               liqa_sdm(kl,i,j) = 0.0_RP           !! <50% in share
            else
               liqa_sdm(kl,i,j) = liqa_sdm(kl,i,j)/drate
            end if

            !! at upper boundary

            ku    = floor(sd_rku(i,j))
            drate = sd_rku(i,j) - real(ku,kind=RP)

            if( drate<0.50_RP ) then
               liqa_sdm(ku,i,j) = 0.0_RP           !! <50% in share
            else
               liqa_sdm(ku,i,j) = liqa_sdm(ku,i,j)/drate
            end if

         end do
         end do

         !=== convert super-droplets to density of rain-water. ===!

!         do k=2,nk-1
!         do j=1,nj-1
!         do i=1,ni-1
         do k=KS,KE+1
         do j=JS-1,JE+1
         do i=IS-1,IE+1

            rhoa_sdm(k,i,j) = liqa_sdm(k,i,j) * dcoef(k,i,j)      &
                    / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)

         end do
         end do
         end do

      end if

!do k=KS,KE
!do j=JS,JE
!do i=IS,IE
!if( rhoc_sdm(k,i,j) == 0.0_RP .and. k < KE/4 ) then
!write(*,*) "hoho", k, i, j, rhoc_sdm(k,i,j), rhoa_sdm(k,i,j)
!endif
!enddo
!enddo
!enddo

    return
  end subroutine sdm_sd2rhocr
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine getexner(pbr,pp,pi,p)

  use scale_const, only: &
      rd => CONST_Rdry, &
      cp => CONST_CPdry, &
      p0 => CONST_PRE00
  implicit none
  ! Input variables
  real(RP),intent(in) :: pbr(KA,IA,JA)  ! Base state pressure
  real(RP),intent(in) :: pp(KA,IA,JA)   ! Pressure perturbation
  ! Output variables
  real(RP),intent(out) :: p(KA,IA,JA)   ! Pressure
  real(RP),intent(out) :: pi(KA,IA,JA)  ! Exner function
  ! Internal shared variables
  real(RP) :: rddvcp      ! rd / cp
  real(RP) :: p0iv        ! 1.0 / p0
  ! Internal private variables
  integer :: i, j, k        ! Array index in x direction
  !--------------------------------------------------------------------
  ! Set the common used variables.
   rddvcp=rd/cp
   p0iv=1.e0/p0

!   do k=1,nk-1
!   do j=1,nj-1
!   do i=1,ni-1
   do k=KS-1,KE+1
   do j=JS-1,JE+1
   do i=IS-1,IE+1
     p(k,i,j)=pbr(k,i,j)+pp(k,i,j)
     pi(k,i,j)=exp(rddvcp*log(p0iv*p(k,i,j)))
   end do
   end do
   end do

    return
  end subroutine getexner
  !-----------------------------------------------------------------------------
  subroutine sdm_calc(MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,              & 
                      sdm_calvar,sdm_mvexchg,dtcl, sdm_aslset,    &
                      exnr_crs,                       &
                      pbr_crs,ptbr_crs,pp_crs, &
                      ptp_crs,qv_crs,prec_crs,zph_crs,rhod_crs,   &
                      lsdmup,ni_sdm,nj_sdm,nk_sdm,                &
                      sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_ri,sd_rj,sd_rk,      &
                      sd_u,sd_v,sd_vz,sd_r,sd_asl,sd_rkl,sd_rku,  &
                      sd_rng,sd_rand,sort_id,sort_key,sort_freq,  &
!                      sd_rand,sort_id,sort_key,sort_freq,         &
                      sort_tag,                                   &
                      bufsiz1,bufsiz2,sdm_itmp1,sdm_itmp2,        &
                      sd_itmp1,sd_itmp2,sd_itmp3,sd_dtmp1,        &
                      crs_val1p,crs_val1c,crs_val2p,crs_val2c,    &
                      crs_val3p,crs_val3c,rbuf,sbuf)

   use scale_process, only: &
       PRC_MPIstop
   use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
   use scale_tracer, only: &
       QAD => QA
   use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod
   use m_sdm_boundary, only: &
        sdm_jdginvdv, sdm_boundary
    use m_sdm_motion, only: &
        sdm_getvz, sdm_getvel, sdm_move
   real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
   real(RP), intent(inout) :: MOMZ(KA,IA,JA)        !! Momentum [kg/s/m2]
   real(RP), intent(inout) :: MOMX(KA,IA,JA)
   real(RP), intent(inout) :: MOMY(KA,IA,JA)
   real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
   real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
   ! Input variables
   logical,  intent(in) :: sdm_calvar(3)
   integer,  intent(in) :: sdm_aslset
   integer,  intent(in) :: sdm_mvexchg
   real(RP), intent(in) :: dtcl(1:3)    ! Time interval of cloud micro physics
!   real(RP), intent(in) :: jcb(KA,IA,JA)      ! Jacobian at scalar points
!   real(RP), intent(in) :: rst_crs(KA,IA,JA)  ! Base state density x Jacobian
   real(RP), intent(in) :: exnr_crs(KA,IA,JA) ! exner function
   real(RP), intent(in) :: pbr_crs(KA,IA,JA)  ! Base state pressure
   real(RP), intent(in) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
   real(RP), intent(in) :: pp_crs(KA,IA,JA)   ! Pressure perturbation
   real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinates
   integer, intent(in) :: ni_sdm                    ! SDM model dimension in x direction
   integer, intent(in) :: nj_sdm                    ! SDM model dimension in y direction
   integer, intent(in) :: nk_sdm                    ! SDM model dimension in z direction
   integer, intent(in) :: sd_num                    ! number of super-droplets
   integer, intent(in) :: sd_numasl     ! number of kind of chemical material contained as water-soluble aerosol in super droplets
   real(RP), intent(in) :: sd_rkl(IA,JA) ! index[k/real] at lower boundary in SDM
   real(RP), intent(in) :: sd_rku(IA,JA) ! index[k/real] at upper boundary in SDM
   integer, intent(in) :: bufsiz1       ! buffer size for MPI
   integer, intent(in) :: bufsiz2       ! buffer size for MPI
   ! Input and output variables
   real(RP), intent(inout) :: rhod_crs(KA,IA,JA)  ! dry air density
   integer(DP), intent(inout) :: sd_n(1:sd_num)    ! multiplicity of super-droplets
   real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
   real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
   real(RP), intent(inout) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
   real(RP), intent(inout) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
   real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets in vertical
   real(RP), intent(inout) :: sd_u(1:sd_num)  ! x-components velocity of super-droplets
   real(RP), intent(inout) :: sd_v(1:sd_num)  ! y-components velocity of super-droplets
   real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets /z velocity of super-droplets
   real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
   real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
   type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
   real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
   integer, intent(inout) :: sort_id(1:sd_num)  ! id that super-droplets sorted by grids
   integer, intent(inout) :: sort_key(1:sd_num) ! sort key
   integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each grid
   integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)  ! accumulated number of super-droplets in each grid
   real(RP), intent(inout) :: ptp_crs(KA,IA,JA)! Potential temperature perturbation
   real(RP), intent(inout) :: qv_crs(KA,IA,JA) ! Water vapor mixing ratio
   real(RP), intent(inout) :: prec_crs(IA,JA,1:2)! [1]:precipitation and [2]:accumlation
   real(DP), intent(inout) :: rbuf(1:bufsiz1,1:bufsiz2,1:2) ! reciving buffer for MPI
   real(DP), intent(inout) :: sbuf(1:bufsiz1,1:bufsiz2,1:2) ! sending buffer for MPI
   ! Output variables
   logical, intent(out) :: lsdmup  ! flag for updating water hydrometeor by SDM
   integer, intent(out) :: sdm_itmp1(1:ni_sdm*nj_sdm*nk_sdm+2) ! temporary array of SDM dimension
   integer, intent(out) :: sdm_itmp2(1:ni_sdm*nj_sdm*nk_sdm+2) ! temporary array of SDM dimension
   integer, intent(out) :: sd_itmp1(1:sd_num,1:nomp) ! temporary array of the size of the number of super-droplets.
   integer, intent(out) :: sd_itmp2(1:sd_num,1:nomp) ! temporary array of the size of the number of super-droplets.
   integer, intent(out) :: sd_itmp3(1:sd_num,1:nomp) ! temporary array of the size of the number o fsuper-droplets.
   real(RP), intent(out) :: sd_dtmp1(1:sd_num) ! temporary array of the size of the number of super-droplets.
   real(RP), intent(out) :: crs_val1p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val1c(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val2p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val2c(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val3p(KA,IA,JA) ! temporary buffer of CReSS dimension
   real(RP), intent(out) :: crs_val3c(KA,IA,JA) ! temporary buffer of CReSS dimension
   ! Internal shared variables
   integer :: istep_sdm        ! step number of SDM
   integer :: istep_evl        ! step number of {condensation/evaporation} process
   integer :: istep_col        ! step number of {stochastic coalescence} process
   integer :: istep_adv        ! step number of {motion of super-droplets} process
   ! Work variables
   integer :: t, n     ! index
   integer :: k,i,j     ! index
   real(RP) :: u_scale(KA,IA,JA)   ! u components of velocity
   real(RP) :: v_scale(KA,IA,JA)   ! v components of velocity
   real(RP) :: w_scale(KA,IA,JA)   ! w components of velocity
   real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
   real(RP) :: rhod_scale(KA,IA,JA) ! dry air density
   real(RP) :: t_scale(KA,IA,JA)    ! Temperature
   !---------------------------------------------------------------------

      ! Initialize and rename variables

      istep_sdm = nclstp(0)                !! maximum step
      istep_evl = nclstp(1)                !! condensation/evaporation
      istep_col = nclstp(2)                !! stochastic coalescence
      istep_adv = nclstp(3)                !! motion of super-droplets

      lsdmup = .false.

      ! Calculate super-droplets process.
      !   1 : condensation / evaporation
      !   2 : stochastic coalescence
      !   3 : motion of super-droplets (advection, terminal velocity)
      do t=1,istep_sdm
         !### Run SDM  ###!

         !=== 1 : condensation / evaporation ===!

         if( sdm_calvar(1) .and.                                 &
             mod(t,istep_sdm/istep_evl)==0 .and. sdm_dtevl>0.d0 ) then

            lsdmup = .true.

            ! get density of liquid-water(qw) before process-1

            call sdm_sd2rhow(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_x,sd_y,sd_r,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { condensation/evaporation } in SDM
            call sdm_condevp(sdm_aslset,            &
                             sdm_aslmw,sdm_aslion,sdm_dtevl,      &
                             pbr_crs,ptbr_crs,pp_crs,             &
                             ptp_crs, qv_crs,                     &
                             sd_num,sd_numasl,sd_x,sd_y,sd_r,sd_asl,&
                             sd_rk)

            ! get density of liquid-water(qw) after process-1
            call sdm_sd2rhow(zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_x,sd_y,sd_r,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! convert the impact of { condensation/evaporation } to
            ! {qv,ptp} in CReSS
            call sdm_sd2qvptp(ptp_crs,qv_crs,exnr_crs,         &
                              rhod_crs,crs_val1p,crs_val1c)

            ! update dry air density

            call sdm_getrhod(pbr_crs,ptbr_crs,pp_crs,ptp_crs,  &
                             qv_crs,rhod_crs)

         end if
         !=== 2 : stochastic coalescence ===!

         if( sdm_calvar(2) .and.                                 &
             mod(t,istep_sdm/istep_col)==0 .and. sdm_dtcol>0.d0 ) then

            lsdmup = .true.

            ! get the terminal velocity of super-droplets

            call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

            call sdm_getvz(pres_scale,rhod_scale,t_scale,            &
                           sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'coales' )
            ! { coalescence } in SDM

            call sdm_coales(sdm_colkrnl,sdm_colbrwn,sdm_aslset,         &
                            sdm_aslrho,sdm_dtcol,                       &
                            pbr_crs,ptbr_crs,pp_crs,ptp_crs,            &
                            zph_crs,                                    &
                            ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl,      &
                            sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_rk,     &
                            sort_id,sort_key,sort_freq,sort_tag,        &
                            sd_rng,sd_rand,                             &
!                            sd_rand,                                    &
                            sdm_itmp1,sdm_itmp2,                        &
                            sd_itmp1(1:sd_num,1),sd_itmp2(1:sd_num,1),  &
                            sd_dtmp1)

         end if


         !=== 3 : motion of super-droplets ===!

         if( sdm_calvar(3) .and.                                 &
             mod(t,istep_sdm/istep_adv)==0 .and. sdm_dtadv>0.d0 ) then

            lsdmup = .true.

            ! get the momentum of super-droplets before moving.

            if( sdm_mvexchg>0 ) then
              write(*,*) "sdm_mvexchg>1 is not applied"
              call PRC_MPIstop
!               call sdm_sd2momnt(ni,nj,nk,zph_crs,                    &
!                                 crs_val1p,crs_val2p,crs_val3p,       &
!                                 sd_num,sd_n,sd_x,sd_y,sd_rk,         &
!                                 sd_u,sd_v,sd_vz,sd_r,sd_itmp1)
!
!               call sdm_adjmomnt(idwbc,idebc,idsbc,idnbc,idxsub,idysub, &
!                                 ni,nj,nk,crs_val1p,crs_val2p,crs_val3p,&
!                                 sd_rkl,sd_rku)
!
            end if

            ! get the terminal velocity of super-droplets
            call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

            call sdm_getvz(pres_scale,rhod_scale,t_scale,            &
                           sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'motion' )

            ! get the moving velocity of super-droplets

            call sdm_rho_mom2uvw(DENS,MOMX,MOMY,MOMZ,u_scale,v_scale,w_scale)

            call sdm_getvel(u_scale,v_scale,w_scale,                   &
                            sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vz)

            ! get the momentum of super-droplets after moving.

            if( sdm_mvexchg>0 ) then
              write(*,*) "sdm_mvexchg>1 is not applied"
              call PRC_MPIstop

!               call s_sdm_sd2momnt(iddx_sdm,iddy_sdm,                   &
!     &                             ni,nj,nk,zph_crs,                    &
!     &                             crs_val1c,crs_val2c,crs_val3c,       &
!     &                             sd_num,sd_n,sd_x,sd_y,sd_rk,         &
!     &                             sd_u,sd_v,sd_vz,sd_r,sd_itmp1)
!
!               call sdm_adjmomnt(idwbc,idebc,idsbc,idnbc,idxsub,idysub, &
!     &                           ni,nj,nk,crs_val1c,crs_val2c,crs_val3c,&
!     &                           sd_rkl,sd_rku)

            end if

            ! { motion } in SDM
            call sdm_move(sdm_dtadv,                         &
                          sd_num,sd_u,sd_v,sd_vz,sd_x,sd_y,sd_rk)

            ! lateral boundary routine in SDM
            ! judge super-droplets as invalid or valid in horizontal

            call sdm_boundary(wbc,ebc,sbc,nbc,                           &
                             sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_rk,     &
                             sd_u,sd_v,sd_vz,sd_r,sd_asl,               &
                             bufsiz1,bufsiz2,sd_itmp1,rbuf,sbuf)

            ! judge super-droplets as invalid or valid in vartical

            call sdm_jdginvdv(sd_rkl,sd_rku,sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk)

            ! convert the impact of { exchange momentum } to
            ! {u,v,w} in CReSS

            if( sdm_mvexchg>0 ) then
              write(*,*) "sdm_mvexchg>1 is not applied"
              call PRC_MPIstop
!               call sdm_sd2uvw(ni,nj,nk,rst_crs,u_crs,v_crs,wc_crs,     &
!     &                         crs_val1p,crs_val2p,crs_val3p,           &
!     &                         crs_val1c,crs_val2c,crs_val3c)
!
            end if

         end if

      end do

      if( .not. lsdmup ) return

      ! Communicate CReSS variables to the neighbor node in horizontal
      ! after running SDM

!      call sdm_commucrs(ni,nj,nk,ptp_crs)
!      call sdm_commucrs(ni,nj,nk,qv_crs)

!      if( sdm_mvexchg>0 ) then
!         call sdm_commucrs(ni,nj,nk,u_crs)
!         call sdm_commucrs(ni,nj,nk,v_crs)
!         call sdm_commucrs(ni,nj,nk,wc_crs)
!      end if

    ! Convert super-droplets to precipitation

      call sdm_sd2prec(dt,                                &
                       prec_crs,                          &
                       sd_num,sd_n,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,  &
                       sd_itmp1(1:sd_num,1),crs_val1c(1,1:IA,1:JA))

! -----

    return
  end subroutine sdm_calc
!---------------------------------------------------------------------------------------
  subroutine sdm_sd2rhow(zph_crs,rhow_sdm,sd_num,sd_n,            &
                         sd_x,sd_y,sd_r,sd_rk,sd_rkl,sd_rku,      &
                         liqw_sdm,ilist)
     use scale_const, only: &
      rw => CONST_DWATR
     use scale_grid, only: &
       FX => GRID_FX, &
       FY => GRID_FY
     ! Input variables
!      integer, intent(in) :: id_dx_sdm
!      integer, intent(in) :: id_dy_sdm  ! Unique index of dy_sdm in namelist table
      real(RP),intent(in) :: zph_crs(KA,IA,JA)    ! z physical coordinate
      integer, intent(in) :: sd_num               ! number of super-droplets
      integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
      real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
      real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
      real(RP), intent(in) :: sd_rk(1:sd_num)! index[k/real] of super-droplets
      real(RP), intent(in) :: sd_rkl(IA,JA) ! lower boundary index[k/real] at scalar point in SDM calculation area
      real(RP), intent(in) :: sd_rku(IA,JA) ! upper boundary index[k/real] at scalar point in SDM calculation area
      ! Output variables
      real(RP), intent(out) :: rhow_sdm(KA,IA,JA) ! densitiy of liquid water
      real(RP), intent(out) :: liqw_sdm(KA,IA,JA) ! liquid water of super-droplets
      integer, intent(out) :: ilist(1:int(sd_num/nomp),1:nomp)   ! buffer for list vectorization
      ! Work variables for OpenMP
      integer :: sd_str        ! index of divided loop by OpenMP
      integer :: sd_end        ! index of divided loop by OpenMP
      integer :: np            ! index for OpenMP
      ! Work variables
      real(RP) :: dcoef(KA,IA,JA)    ! coef.
      real(RP) :: drate        ! temporary
!      integer(kind=i8) :: idx_sdm   ! integer of 'dx_sdm'
!      integer(kind=i8) :: idy_sdm   ! integer of 'dy_sdm'
      integer :: nlist(1:nomp)      ! list number
      integer :: cnt                ! counter
      integer :: i, j, k, kl, ku, m, n, ix, jy    ! index
     !--------------------------------------------------------------------

      ! Initialize
      do k = 1, KA
      do i = 1, IA
      do j = 1, JA
       dcoef(k,i,j) = rw * F_THRD * ONE_PI / real(dx_sdm(i)*dy_sdm(j),kind=RP)
!      idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!      idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )
      enddo
      enddo
      enddo

      do np=1,nomp
         nlist(np) = 0
      end do

!      do k=1,nk
!      do j=0,nj+1
!      do i=0,ni+1
      do k=2,KA
      do j=1,JA
      do i=1,IA
         liqw_sdm(k,i,j) = 0.0_RP
         rhow_sdm(k,i,j) = 0.0_RP
      end do
      end do
      end do


      ! Get index list for compressing buffer.
      do np=1,nomp

         sd_str = int(sd_num/nomp)*(np-1) + 1
         sd_end = int(sd_num/nomp)*np
         cnt    = 0

         do n=sd_str,sd_end
             if( sd_rk(n)>VALID2INVALID ) then
               cnt = cnt + 1
               ilist(cnt,np) = n
            end if
         end do

         nlist(np) = cnt

      end do

      ! Get density of liquid-water.
      !### count voulme of super-droplets ###!
      do np=1,nomp
         if( nlist(np)>0 ) then
            do m=1,nlist(np)
               n = ilist(m,np)
!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
               do ix = IS, IE
                if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                 i = ix
                 exit
                endif
               enddo
               do jy = IS, JE
                if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                 j = jy
                 exit
                endif
               enddo
               k = floor( sd_rk(n) )

               liqw_sdm(k,i,j) = liqw_sdm(k,i,j)                        &
     &                         + sd_r(n) * sd_r(n) * sd_r(n)            &
     &                         * real(sd_n(n),kind=RP)
            end do
         end if
      end do

      ! Adjust liquid-water in verical boundary.
!      do j=2,nj-2
!      do i=2,ni-2
      do j=JS,JE
      do i=IS,IE
         !! at lower boundary
         kl    = floor(sd_rkl(i,j))
         drate = real(kl+1,kind=RP) - sd_rkl(i,j)
         if( drate<0.50_RP ) then
            liqw_sdm(kl,i,j) = 0.0_RP           !! <50% in share
         else
            liqw_sdm(kl,i,j) = liqw_sdm(kl,i,j)/drate
         end if

         !! at upper boundary
         ku    = floor(sd_rku(i,j))
         drate = sd_rku(i,j) - real(ku,kind=RP)
         if( drate<0.50_RP ) then
            liqw_sdm(ku,i,j) = 0.0_RP           !! <50% in share
         else
            liqw_sdm(ku,i,j) = liqw_sdm(ku,i,j)/drate
         end if
      end do
      end do

      ! Convert super-droplets to density of liquid-water.
!      do k=2,nk-1
!      do j=1,nj-1
!      do i=1,ni-1
      do k=KS,KE+1
      do j=JS-1,JE+1
      do i=IS-1,IE+1
         rhow_sdm(k,i,j) = liqw_sdm(k,i,j) * dcoef(k,i,j)             &
                         / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)
      end do
      end do
      end do


    return
  end subroutine sdm_sd2rhow
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2qvptp(ptp_crs,qv_crs,exnr_crs,       &
                          rhod_crs,rhowp_sdm,rhowc_sdm)
     use scale_const, only: &
         cp => CONST_CPdry
     ! Input variables
     real(RP), intent(in) :: exnr_crs(KA,IA,JA) ! Exner function
     real(RP), intent(in) :: rhod_crs(KA,IA,JA) ! dry air density
     real(RP), intent(in) :: rhowp_sdm(KA,IA,JA) ! density of liquid water at preveous step
     real(RP), intent(in) :: rhowc_sdm(KA,IA,JA) ! density of liquid water at current step
     ! Output variables
     real(RP), intent(inout) :: ptp_crs(KA,IA,JA) ! Potential temperature perturbation
     real(RP), intent(inout) :: qv_crs(KA,IA,JA)  ! Water vapor mixing ratio
     ! Work variables
     real(RP) :: sv ! the source term of water vapor associated with condenstation/evaporation process.
     integer :: i, j, k ! index
  !-------------------------------------------------------------------7--
  ! Input the impact of cloud physics using SDM to CReSS.
!      do k=2,nk-1
!      do j=1,nj-1
!      do i=1,ni-1
      do k=KS,KE+1
      do j=JS-1,JE+1
      do i=IS-1,IE+1
         sv = -(rhowc_sdm(k,i,j)-rhowp_sdm(k,i,j))/rhod_crs(k,i,j)
!if( qv_crs(k,i,j) + real(sv) < 0.0_RP ) then
! write(*,*) k,i,j, qv_crs(k,i,j), sv
!endif
!         qv_crs(k,i,j)  = qv_crs(k,i,j) + real(sv)
!ORG     qv_crs(i,j,k)  = max( qv_crs(i,j,k) + real(sv), 1.e-8 )
         qv_crs(k,i,j)  = max( qv_crs(k,i,j) + real(sv), 1.E-8_RP )
         ptp_crs(k,i,j) = ptp_crs(k,i,j)                                &
                        - real( LatHet*sv ) / ( cp*exnr_crs(k,i,j) )
      end do
      end do
      end do


    return
  end subroutine sdm_sd2qvptp
  !-----------------------------------------------------------------------------
  subroutine sdm_getrhod(pbr_crs,ptbr_crs,pp_crs,ptp_crs,&
                         qv_crs,rhod_crs)
    use scale_const, only: &
      rd => CONST_Rdry,  &
      cp => CONST_CPdry, &
      p0 => CONST_PRE00       ! Reference Pressure [Pa]
    ! Input variables
    real(RP), intent(in)  :: pbr_crs(KA,IA,JA)  ! Base state pressure
    real(RP), intent(in)  :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
    real(RP), intent(in)  :: pp_crs(KA,IA,JA)   ! Pressure perturbation
    real(RP), intent(in)  :: ptp_crs(KA,IA,JA)  ! Potential temperature perturbation
    real(RP), intent(in)  :: qv_crs(KA,IA,JA)   ! Water vapor mixing ratio
    ! Output variables
    real(RP), intent(out) :: rhod_crs(KA,IA,JA)  ! dry air densitiy
    ! Work variables
    real(RP) :: rddvcp    ! rd / cp
    real(RP) :: p_crs     ! pressure
    real(RP) :: pt_crs    ! potential temperature
    real(RP):: t_crs     ! temperature
    real(RP):: tv_crs    ! virtual temperature
    integer :: i, j, k, n  ! index
    real(RP), parameter :: epsav = 1.60770_RP   ! Molecular weight ratio of vapor/air
 !---------------------------------------------------------------------
 ! Set constant variable
    rddvcp = rd / cp

 ! Get dry air density.
!    do k=2,nk-1
!    do n=0,(ni-1)*(nj-1)-1
    do k=KS,KE+1
    do n=1,(IE+1)*(JE+1)-1
         !### loop coupling for vectorization ###!
!         i = mod(n,ni-1) + 1        !! do j=1,nj-1
!         j = (n-(i-1))/(ni-1) + 1   !! do i=1,ni-1
         i = mod(n,IE+1) + 1        !! do j=1,nj-1
         j = (n-(i-1))/(IE+1) + 1   !! do i=1,ni-1
         pt_crs = ptbr_crs(k,i,j) + ptp_crs(k,i,j)
         p_crs  = pbr_crs(k,i,j)  + pp_crs(k,i,j)
         t_crs  = pt_crs * exp(rddvcp*log(p_crs/p0))  !! pt => t
         !! Get virtual temperature
         tv_crs = t_crs                                                 &
                * (1.e0+epsav*qv_crs(k,i,j))/(1.e0+qv_crs(k,i,j))

         rhod_crs(k,i,j) = real( (p_crs/(tv_crs*rd))                    &
                                       /(1.e0+qv_crs(k,i,j)), kind=RP )

!DBG     rhod_crs(i,j,k) = real( p_crs / (tv_crs*rd), kind=r8 )
                                   !! density of moist air (rhod+rhov)
    end do
    end do

    return
  end subroutine sdm_getrhod
  !-----------------------------------------------------------------------------
  subroutine sdm_coales(sdm_colkrnl,sdm_colbrwn,               &
                        sdm_aslset,sdm_aslrho,                 &
                        sdm_dtcol,pbr_crs,ptbr_crs,            &
                        pp_crs,ptp_crs,zph_crs,                &
                        ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl, &
                        sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_rk,&
                        sort_id,sort_key,sort_freq,sort_tag,   &
                        sd_rng,sd_rand,                        &
!                        sd_rand,                               &
                        sort_tag0,fsort_id,icp,sd_perm,c_rate  )
      use scale_grid, only: &
         FX => GRID_FX, &
         FY => GRID_FY
      use scale_process, only: &
         PRC_MPIstop
!      use m_gadg_algorithm, only: &
!         gadg_count_sort
      use scale_const, only:  &
         t0 => CONST_TEM00, &    ! Gas Constant of vapor [J/K/kg]
         rw => CONST_Rvap,  &    ! Gas Constant of vapor [J/K/kg]
         rd => CONST_Rdry,  &    ! Gas Constant of dry air [J/K/kg]
         cp => CONST_CPdry, &
         p0 => CONST_PRE00       ! Reference Pressure [Pa]
      !  Input variables
      integer, intent(in) :: sdm_colkrnl   ! Kernel type for coalescence process
      integer, intent(in) :: sdm_colbrwn   ! Control flag of Brownian Coagulation and Scavenging process
      integer, intent(in) :: sdm_aslset    ! Control flag to set species and way of chemical material as water-soluble aerosol
      real(RP),intent(in) :: sdm_aslrho(20)! User specified density of chemical material contained as water-soluble aerosol in super droplets
      real(RP), intent(in) :: sdm_dtcol   ! tims step of {stochastic coalescence} process
      real(RP), intent(in) :: pbr_crs(KA,IA,JA) ! Base state pressure
      real(RP), intent(in) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
      real(RP), intent(in) :: pp_crs(KA,IA,JA) ! Pressure perturbation
      real(RP), intent(in) :: ptp_crs(KA,IA,JA) ! Potential temperature perturbation
      real(RP), intent(in) :: zph_crs(KA,IA,JA) ! z physical coordinate
      integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
      integer, intent(in) :: sd_num  ! number of super-droplets
      integer, intent(in) :: sd_numasl ! Number of kind of chemical material contained as water-soluble aerosol in super droplets
      real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
      real(RP), intent(in) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng ! random number generator
      real(RP),intent(inout) :: sd_rand(1:sd_num) ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num) ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num) ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2) ! accumulated number of super-droplets in each SD-grid
      integer(DP), intent(inout) :: sd_n(1:sd_num) ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)  ! = sort_tag(n) - 1
      integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
      integer, intent(out) :: icp(1:sd_num) ! index of coalescence pair
      integer, intent(out) :: sd_perm(1:sd_num) ! random permutations
      real(RP), intent(out) :: c_rate(1:sd_num) ! coalescence probability
      ! Internal shared variables
      real(RP) :: sd_aslrho(1:22) ! Density of chemical material contained as water-soluble aerosol in super droplets
      integer(RP) :: sd_ncol ! how many times coalescence occurs
      integer :: freq_max ! get the maximum number of super-droplets in each grid
      integer :: hfreq_max ! hfreq_max / 2
      integer :: ipremium  ! premium coef. for coalescence
      ! Work variables
      real(RP) :: dmask(1:22) ! mask for vactorization
      real(RP) :: sd_asl1(1:sd_numasl)
      real(RP) :: sd_asl2(1:sd_numasl) ! aerosol mass of super-droplets with large/small multiplicity
      real(RP) :: sd_cc1  ! slip correction of super-droplets
      real(RP) :: sd_cc2  ! with large/small multiplicity
      real(RP) :: sd_dia1 ! diameter of super-droplets
      real(RP) :: sd_dia2 ! with large/small multiplicity
      real(RP) :: sd_lmd1 ! mean free path of super-droplets
      real(RP) :: sd_lmd2 ! with large/small multiplicity
      real(RP) :: sd_m1   ! mass of super-droplets
      real(RP) :: sd_m2   ! with large/small multiplicity
      real(RP) :: sd_r1   ! radius of super-droplets
      real(RP) :: sd_r2   ! with large/small multiplicity
      real(RP) :: sd_rk1  ! index[k/real] of super-droplets with large multiplicity
      real(RP) :: sd_rw1  ! radius of water parts in super-droplets
      real(RP) :: sd_rw2  ! with large/small multiplicity
      real(RP) :: sd_tmasl1
      real(RP) :: sd_tmasl2 ! total mass of aerosol part in super droplets with large/small multiplicity
      real(RP) :: sd_tvasl1 ! total volume of aerosol part in super
      real(RP) :: sd_tvasl2 ! droplets with large/small multiplicity
      real(RP) :: sd_v1   ! volume of super-droplets
      real(RP) :: sd_v2   ! with large/small multiplicity
      real(RP) :: sd_c1   ! temporary
      real(RP) :: sd_c2
      real(RP) :: sd_d1
      real(RP) :: sd_d2
      real(RP) :: sd_g1
      real(RP) :: sd_g2
      real(RP) :: lmd_crs ! air mean free path
      real(RP) :: p_crs   ! pressure
      real(RP) :: pt_crs  ! potential temperarure
      real(RP) :: t_crs   ! temperarure
      real(RP) :: vis_crs ! dynamic viscosity
      real(RP) :: sumdia  ! sum of variables of a pair of droplets
      real(RP) :: sumd
      real(RP) :: sumc
      real(RP) :: sumg
      real(RP) :: sumr
      real(RP) :: k12     ! brownian coagulation coefficient
      real(RP) :: dvz     ! difference in terminal velocity of a pair of super-droplets
      real(RP) :: dtmp    ! temporary
      real(RP) :: frac    ! fraction parts
      real(RP) :: ivvol(KA,IA,JA)   ! inverse of a grid volume
      real(RP) :: tdeg    ! temperature in degree
      real(RP) :: rq      ! radius ratio of a pair of super-droplets
      real(RP) :: ek      ! temporary
      real(RP) :: p       ! temporary
      real(RP) :: q       ! temporary

      integer(DP) :: idx_sdm  ! integer of 'dx_sdm'
      integer(DP) :: idy_sdm  ! integer of 'dy_sdm'

      integer(DP) :: sd_nmax  ! maximum multiplicity
      integer(DP) :: sd_n1    ! multiplicity of super-droplets with large multiplicity
      integer(DP) :: sd_n2    ! multiplicity of super-droplets  with small multiplicity

      integer, allocatable :: fsort_tag(:) ! buffer for sorting
      integer, allocatable :: fsort_freq(:) ! buffer for sorting

      integer :: idx_nasl(1:22)  ! index for vactorization


      integer :: gnum          ! grid number
      integer :: iexced        ! temporary

      integer :: in1, in2, in3, in4 ! index
      integer :: irr, iqq           ! index of coef.
      integer :: i, j, k, m, n, s   ! index
      integer :: t, tc, tp          ! index
      integer :: ix, jy
!--------------------------------------------------------------------

      ! Initialize
      gnum = ni_sdm * nj_sdm * knum_sdm

!      idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!      idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )

      freq_max = 1
      iexced   = 1   !! check for memory size of int*8

      !### aerosol type ###!

      do n=1,22
         sd_aslrho(n) = 1.0_RP
      end do

      if( abs(mod(sdm_aslset,10))==1 ) then

         !### numasl=1 @ init+rest : (NH4)2SO4 ###!

         sd_aslrho(1) = real(rho_amsul,kind=RP)

      else if( abs(mod(sdm_aslset,10))==2 ) then

         if( abs(sdm_aslset)==2 ) then

            !### numasl=1 @ init : NaCl ###!

            sd_aslrho(1) = real(rho_nacl,kind=RP)

         else if( abs(sdm_aslset)==12 ) then

            !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

            sd_aslrho(1) = real(rho_amsul,kind=RP)
            sd_aslrho(2) = real(rho_nacl,kind=RP)

         end if
      else if( abs(mod(sdm_aslset,10))==3 ) then

         !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

         sd_aslrho(1) = real(rho_amsul,kind=RP)
         sd_aslrho(2) = real(rho_nacl,kind=RP)

!         do n=1,20
!            call getrname( id_sdm_aslrho + (n-1), sd_aslrho(n+2) )
!         end do
         do n=1,20
           sd_aslrho(n+2) = sdm_aslrho(n)
         end do

      end if

    ! Sorting super-droplets.

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_x,sd_y,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

    ! Initialize

      do n=1,22

         if( n<=sd_numasl ) then
            idx_nasl(n) = n
            dmask(n) = 1.0_RP
         else
            idx_nasl(n) = sd_numasl
            dmask(n) = 0.0_RP
         end if

      end do

      do n=1,gnum+2
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      do n=1,sd_num
         c_rate(n) = 0.0_RP
      end do

     ! Get the maximum number of super-droplets in each grid

      do m=1,gnum
         freq_max = max( freq_max, sort_freq(m) )
      end do

      hfreq_max = int(freq_max/2)

     ! Sorting the grids by the number of super-droplets in each grid

      allocate( fsort_tag(0:freq_max+1) )
      allocate( fsort_freq(0:freq_max)  )

      call gadg_count_sort( sort_freq(1:gnum), 0, freq_max, &
                            fsort_freq, fsort_tag, fsort_id )

      fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

    ! Get random number using random number generator

      call gen_rand_array( sd_rng, sd_rand )

    ! Get random permutation layout of super-droplets in each grid

      call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
     &                 sort_tag0,fsort_tag,fsort_id,sd_rand,sd_perm)

    ! Get random number using random number generator

      call gen_rand_array( sd_rng, sd_rand )

    ! Select a pair of super-droples in random permutation layout

      do n=1,hfreq_max

         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            in1 = 2 * n
            in2 = in1 - 1

            !### select pair of super-droplets in each grid ###!

            in3 = sd_perm( sort_tag0(m) + in1 )
            in4 = sd_perm( sort_tag0(m) + in2 )

            !### set the random index ###!
            tc = sort_tag0(m) + n
            tp = tc + int(sort_freq(m)/2)

            icp(tc) = sort_id( sort_tag0(m) + in3 )
            icp(tp) = sort_id( sort_tag0(m) + in4 )

         end do

      end do

     ! Get effective collision of droplets that exceed micrometer-size
     ! ( over 10um )

      if( sdm_colkrnl==0 ) then

         !### Golovin's kernel (m^3/s?) ###!
         do n=1,hfreq_max

            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

               m = fsort_id(t)

               tc = sort_tag0(m) + n
               tp = tc + int(sort_freq(m)/2)

               sd_r1 = sd_r(icp(tc))
               sd_r2 = sd_r(icp(tp))

               c_rate(tc) = 1500.0_RP * 1.333333_RP * ONE_PI              &
                          * ( sd_r1*sd_r1*sd_r1 + sd_r2*sd_r2*sd_r2 )

            end do

         end do

      else if( sdm_colkrnl==1 ) then

         !### Long's kernel (m^2)  ###!

         do n=1,hfreq_max

            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

               m = fsort_id(t)

               tc = sort_tag0(m) + n
               tp = tc + int(sort_freq(m)/2)

               sd_r1 = max( sd_r(icp(tc)), sd_r(icp(tp)) )  !! large
               sd_r2 = min( sd_r(icp(tc)), sd_r(icp(tp)) )  !! small

               if( sd_r1 <= 5.E-5_RP ) then
                  c_rate(tc) = 4.5E+8_RP * ( sd_r1*sd_r1 )                 &
                             * ( 1.0_RP - 3.E-6_RP/(max(3.01E-6_RP,sd_r1)) )
               else
                  c_rate(tc) = 1.d0
               end if

               sumr = sd_r1 + sd_r2

               c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr)

            end do

         end do

      else if( sdm_colkrnl==2 ) then

         !### Hall's kernel (m^2) ###!

         do n=1,hfreq_max

            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

               m = fsort_id(t)

               tc = sort_tag0(m) + n
               tp = tc + int(sort_freq(m)/2)

               sd_r1 = max( sd_r(icp(tc)), sd_r(icp(tp)) )  !! large
               sd_r2 = min( sd_r(icp(tc)), sd_r(icp(tp)) )  !! small

               rq    = sd_r2 / sd_r1
               sd_r1 = sd_r1 * m2micro    !! [m] => [micro-m]
               sd_r2 = sd_r2 * m2micro    !! [m] => [micro-m]

               !! Get index of the array {r0,rat}.

               if( sd_r1 <= r0col(1) ) then
                  irr = 1
               else if( sd_r1 <= r0col(2) ) then
                  irr = 2
               else if( sd_r1 <= r0col(3) ) then
                  irr = 3
               else if( sd_r1 <= r0col(4) ) then
                  irr = 4
               else if( sd_r1 <= r0col(5) ) then
                  irr = 5
               else if( sd_r1 <= r0col(6) ) then
                  irr = 6
               else if( sd_r1 <= r0col(7) ) then
                  irr = 7
               else if( sd_r1 <= r0col(8) ) then
                  irr = 8
               else if( sd_r1 <= r0col(9) ) then
                  irr = 9
               else if( sd_r1 <= r0col(10) ) then
                  irr = 10
               else if( sd_r1 <= r0col(11) ) then
                  irr = 11
               else if( sd_r1 <= r0col(12) ) then
                  irr = 12
               else if( sd_r1 <= r0col(13) ) then
                  irr = 13
               else if( sd_r1 <= r0col(14) ) then
                  irr = 14
               else if( sd_r1 <= r0col(15) ) then
                  irr = 15
               else
                  irr = 16
               end if

               if( rq <= ratcol(2) ) then
                  iqq = 2
               else if( rq <= ratcol(3) ) then
                  iqq = 3
               else if( rq <= ratcol(4) ) then
                  iqq = 4
               else if( rq <= ratcol(5) ) then
                  iqq = 5
               else if( rq <= ratcol(6) ) then
                  iqq = 6
               else if( rq <= ratcol(7) ) then
                  iqq = 7
               else if( rq <= ratcol(8) ) then
                  iqq = 8
               else if( rq <= ratcol(9) ) then
                  iqq = 9
               else if( rq <= ratcol(10) ) then
                  iqq = 10
               else if( rq <= ratcol(11) ) then
                  iqq = 11
               else if( rq <= ratcol(12) ) then
                  iqq = 12
               else if( rq <= ratcol(13) ) then
                  iqq = 13
               else if( rq <= ratcol(14) ) then
                  iqq = 14
               else if( rq <= ratcol(15) ) then
                  iqq = 15
               else if( rq <= ratcol(16) ) then
                  iqq = 16
               else if( rq <= ratcol(17) ) then
                  iqq = 17
               else if( rq <= ratcol(18) ) then
                  iqq = 18
               else if( rq <= ratcol(19) ) then
                  iqq = 19
               else if( rq <= ratcol(20) ) then
                  iqq = 20
               else
                  iqq = 21
               end if
               !! Get c_rate

               if( irr>=16 ) then

                  q  = (rq-ratcol(iqq-1)) / (ratcol(iqq)-ratcol(iqq-1))
                  ek = (1.0_RP-q)*ecoll(15,iqq-1) + q*ecoll(15,iqq)

                  c_rate(tc) = min( ek, 1.0_RP )

               else if( irr>=2 .and. irr<16 ) then

                  p = (sd_r1-r0col(irr-1))/(r0col(irr)-r0col(irr-1))
                  q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

                  c_rate(tc) = (1.0_RP-p)*(1.0_RP-q)*ecoll(irr-1,iqq-1)     &
                             + p*(1.0_RP-q)*ecoll(irr,iqq-1)              &
                             + q*(1.0_RP-p)*ecoll(irr-1,iqq)              &
                             + p*q*ecoll(irr,iqq)

               else

                  q = (rq-ratcol(iqq-1))/(ratcol(iqq)-ratcol(iqq-1))

                  c_rate(tc) = (1.0_RP-q)*ecoll(1,iqq-1) + q*ecoll(1,iqq)

               end if

               sd_r1 = sd_r1 * micro2m    !! [micro-m] => [m]
               sd_r2 = sd_r2 * micro2m    !! [micro-m] => [m]

               !! Get c_rate

               sumr = sd_r1 + sd_r2

               c_rate(tc) = c_rate(tc) * ONE_PI * (sumr*sumr)

            end do

         end do

      else if( sdm_colkrnl==3 ) then

         !### no coalescence effeciency hydrodynamic kernel (m^2) ###!
         do n=1,hfreq_max

            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

               m = fsort_id(t)

               tc = sort_tag0(m) + n
               tp = tc + sort_freq(m)/2

               sumr = sd_r(icp(tc)) + sd_r(icp(tp))

               c_rate(tc) = ONE_PI * (sumr*sumr)

            end do

         end do

      end if

      ! -----

      if( sdm_colkrnl==0 ) then

         !### Golovin's kernel [-] ###!
         do n=1,hfreq_max
            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1
               m = fsort_id(t)

               tc = sort_tag0(m) + n
               tp = tc + sort_freq(m)/2

               !! Get location of Super-Droplets

!               i = int( floor(sd_x(icp(tc))*1.d5,kind=i8)/idx_sdm ) + 2
!               j = int( floor(sd_y(icp(tc))*1.d5,kind=i8)/idy_sdm ) + 2
               do ix = IS, IE
                if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                 i = ix
                 exit
                endif
               enddo
               do jy = JS, JE
                if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                 j = jy
                 exit
                endif
               enddo
               k = floor( sd_rk(icp(tc)) )

               ivvol(k,i,j) = 1.0_RP / real(dx_sdm(i)*dy_sdm(j),kind=RP)               &
                                     / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)

               c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
            end do
         end do

      else

         !### Long's kernel, Hall's kernel,       ###!
         !### no col_effi hydrodynamic kernel [-] ###!

         do n=1,hfreq_max

            do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

               m = fsort_id(t)

               tc = sort_tag0(m) + n
               tp = tc + sort_freq(m)/2

               dvz = abs( sd_vz(icp(tc)) - sd_vz(icp(tp)) )

               !! Get location of Super-Droplets

!               i = int( floor(sd_x(icp(tc))*1.d5,kind=i8)/idx_sdm ) + 2
!               j = int( floor(sd_y(icp(tc))*1.d5,kind=i8)/idy_sdm ) + 2
               do ix = IS, IE
                if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                 i = ix
                 exit
                endif
               enddo
               do jy = JS, JE
                if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                 j = jy
                 exit
                endif
               enddo
               k = floor( sd_rk(icp(tc)) )
               ivvol(k,i,j) = 1.0_RP / real(dx_sdm(i)*dy_sdm(j),kind=RP)               &
                                     / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)

               c_rate(tc) = c_rate(tc) * real(sdm_dtcol,kind=RP)        &
                          * ivvol(k,i,j) * dvz

            end do

         end do

      end if

     ! Get effective collision for "Brownian Coagulation and Scavenging
     ! (Seinfeld & Pandis,2006)" of droplets less than
     ! micrometer-size ( below 1um )
      if( sdm_colbrwn>0 ) then
        do n=1,hfreq_max

          do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            tc = sort_tag0(m) + n
            tp = tc + sort_freq(m)/2

            !### information of super-droplets ###

            !! radius of water parts in droplets

            sd_rw1 = sd_r(icp(tc))    !! [m]
            sd_rw2 = sd_r(icp(tp))

            !! mass and volume of aerosol parts in droplets

            sd_tmasl1 = 0.0_RP
            sd_tmasl2 = 0.0_RP
            sd_tvasl1 = 0.0_RP
            sd_tvasl2 = 0.0_RP

            do k=1,22

               s = idx_nasl(k)

               sd_tmasl1 = sd_tmasl1 + sd_asl(icp(tc),s) * dmask(k)
               sd_tmasl2 = sd_tmasl2 + sd_asl(icp(tp),s) * dmask(k)

               sd_tvasl1 = sd_tvasl1                                    &
                         + sd_asl(icp(tc),s)/sd_aslrho(s) * dmask(k)
               sd_tvasl2 = sd_tvasl2                                    &
                         + sd_asl(icp(tp),s)/sd_aslrho(s) * dmask(k)

            end do

            sd_tmasl1 = sd_tmasl1  * 1.E-3_RP    !! [g]=>[kg]
            sd_tmasl2 = sd_tmasl2  * 1.E-3_RP

            !! diameter and mass and volume of droplets

            dtmp = ONE_PI * F_THRD

            sd_v1 = sd_tvasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1)
            sd_v2 = sd_tvasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2)

            sd_m1 = sd_tmasl1 + dtmp * (sd_rw1*sd_rw1*sd_rw1) * rw
            sd_m2 = sd_tmasl2 + dtmp * (sd_rw2*sd_rw2*sd_rw2) * rw

            sd_dia1 = (6.0_RP*sd_v1/ONE_PI)**O_THRD
            sd_dia2 = (6.0_RP*sd_v2/ONE_PI)**O_THRD

            !### location of super-droplets ###

!            i = int( floor(sd_x(icp(tc))*1.d5,kind=i8)/idx_sdm ) + 2
!            j = int( floor(sd_y(icp(tc))*1.d5,kind=i8)/idy_sdm ) + 2
            do ix = IS, IE
             if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
              i = ix
              exit
             endif
            enddo
            do jy = JS, JE
             if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
              j = jy
              exit
             endif
            enddo
            k = floor( sd_rk(icp(tc)) )

            ivvol(k,i,j) = 1.0_RP / real(dx_sdm(i)*dy_sdm(j),kind=RP)             &
                                  / real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)
            pt_crs = ptbr_crs(k,i,j) + ptp_crs(k,i,j)
            p_crs  = pbr_crs(k,i,j)  + pp_crs(k,i,j)      !! [Pa]
            t_crs  = pt_crs * exp((rd/cp)*log(p_crs/p0))  !! [K]

            !### dynamic viscosity [Pa*s]  ###!
            !### (Pruppacher & Klett,1997) ###!
            tdeg = t_crs - t0     !! [K] => [degC]
            if( tdeg>=0.0_RP ) then
               vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
            else
               vis_crs = ( 1.7180_RP + 4.9E-3_RP*tdeg                        &
                                     -1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
            end if

            !### air mean free path [m] ###!
            dtmp = dsqrt(8.0_RP*mass_air*1.E-3_RP/(ONE_PI*rrst*t_crs))
            lmd_crs = (2.0_RP*vis_crs)/(p_crs*dtmp)

            !### slip correction of droplets [-]  ###!
            dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia1/lmd_crs)
            sd_cc1 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia1)
            dtmp   = 1.2570_RP + 0.40_RP * exp(-0.550_RP*sd_dia2/lmd_crs)
            sd_cc2 = 1.0_RP + (2.0_RP*lmd_crs*dtmp)/(sd_dia2)

            !### diffusion term [m*m/s] ###!
            dtmp = (boltz*t_crs)/(3.0_RP*ONE_PI*vis_crs)
            sd_d1 = dtmp * (sd_cc1/sd_dia1)
            sd_d2 = dtmp * (sd_cc2/sd_dia2)

            !### velocity term [m/s] ###!
            dtmp = (8.0_RP*boltz*t_crs)/ONE_PI
            sd_c1 = dsqrt(dtmp/sd_m1)
            sd_c2 = dsqrt(dtmp/sd_m2)

            !### mean free path of droplets [m] ###!
            dtmp = 8.0_RP/ONE_PI
            sd_lmd1 = dtmp * (sd_d1/sd_c1)
            sd_lmd2 = dtmp * (sd_d2/sd_c2)

            !### length term [m] ###!
            dtmp = (sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)*(sd_dia1+sd_lmd1)&
                 - (sd_dia1*sd_dia1+sd_lmd1*sd_lmd1)**1.50_RP
            sd_g1 = dtmp/(3.0_RP*sd_dia1*sd_lmd1) - sd_dia1
            dtmp = (sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)*(sd_dia2+sd_lmd2)&
                 - (sd_dia2*sd_dia2+sd_lmd2*sd_lmd2)**1.50_RP
            sd_g2 = dtmp/(3.0_RP*sd_dia2*sd_lmd2) - sd_dia2

            !### Brownian Coagulation Coefficient K12 [m3/s] ###!
            sumdia = sd_dia1 + sd_dia2
            sumd   = sd_d1   + sd_d2
            sumc   = dsqrt( sd_c1*sd_c1 + sd_c2*sd_c2 )
            sumg   = dsqrt( sd_g1*sd_g1 + sd_g2*sd_g2 )

            dtmp = sumdia/(sumdia+2.0_RP*sumg) + (8.0_RP*sumd)/(sumdia*sumc)
            k12 = 2.0_RP*ONE_PI * sumdia*sumd/dtmp

            !### add effective collision [-] ###!
            c_rate(tc) = c_rate(tc)                                     &
                       + k12 * real(sdm_dtcol,kind=RP) * ivvol(k,i,j)
          end do

        end do

      end if

     ! Get total effective collision of droplets

      do n=1,hfreq_max
         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            tc = sort_tag0(m) + n
            tp = tc + sort_freq(m)/2

            sd_nmax  = max( sd_n(icp(tc)), sd_n(icp(tp)) )
                               ! maximum multiplicity
            ipremium = sort_freq(m) - 1 + iand(sort_freq(m),1)
                               ! IAND(sort_freq(i),1) => even:0, odd:1

            c_rate(tc) = c_rate(tc) * real( sd_nmax*ipremium, kind=RP )
         end do
      end do

     ! Stochastic coalescence process.
      do n=1,hfreq_max
         do t=fsort_tag(n*2),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            tc = sort_tag0(m) + n
            tp = tc + sort_freq(m)/2
            !### set coalescence count ###!
            sd_ncol = int( c_rate(tc), kind=RP )
            frac = c_rate(tc) - real( sd_ncol, kind=RP )

            ! judge coalescence by random number and fractional part

            if( sd_rand(tc) < frac ) then
               sd_ncol =  sd_ncol + 1
            end if

            if( sd_ncol<=0 ) cycle  !! no coalesecense

            !### coalescence procudure ###!

            if( sd_n(icp(tc)) > sd_n(icp(tp)) ) then

               sd_n1  = sd_n( icp(tc) )
               sd_r1  = sd_r( icp(tc) )
               sd_rk1 = sd_rk( icp(tc) )
               sd_m1  = sd_r1 * sd_r1 * sd_r1

               sd_n2  = sd_n( icp(tp) )
               sd_r2  = sd_r( icp(tp) )
               sd_m2  = sd_r2 * sd_r2 * sd_r2

               do k=1,22
                  s = idx_nasl(k)
                  sd_asl1(s) = sd_asl( icp(tc),s )
                  sd_asl2(s) = sd_asl( icp(tp),s )
               end do

            else

               sd_n1  = sd_n( icp(tp) )
               sd_r1  = sd_r( icp(tp) )
               sd_rk1 = sd_rk( icp(tp) )
               sd_m1  = sd_r1 * sd_r1 * sd_r1

               sd_n2  = sd_n( icp(tc) )
               sd_r2  = sd_r( icp(tc) )
               sd_m2  = sd_r2 * sd_r2 * sd_r2
               do k=1,22
                  s = idx_nasl(k)
                  sd_asl1(s) = sd_asl( icp(tp),s )
                  sd_asl2(s) = sd_asl( icp(tc),s )
               end do

            end if

            sd_ncol = min( sd_ncol, int(sd_n1/sd_n2,kind=RP) )

            if( sd_n1 > sd_n2*sd_ncol ) then

               sd_n1 = sd_n1 - sd_n2*sd_ncol
               sd_r2 = exp( O_THRD                                      &
                            * log(sd_m1*real(sd_ncol,kind=RP)+sd_m2) )
!ORG           sd_r2 = ( sd_m1*real(sd_ncol,kind=r8) + sd_m2 ) ** O_THRD

               do k=1,22
                  s = idx_nasl(k)
                  dtmp = sd_asl1(s) * real(sd_ncol,kind=RP)
                  sd_asl2(s) = sd_asl2(s) + dmask(k) * dtmp
               end do

            else

               !! coalescent SDs with same multiplicity
               !!  - do not change the order of equations for
               !!  - vectorization

               sd_n1 = int( sd_n2/2, kind=RP )
               sd_n2 = sd_n2 - sd_n1

               sd_r1 = exp( O_THRD                                      &
                            * log(sd_m1*real(sd_ncol,kind=RP)+sd_m2) )
!ORG           sd_r1 = ( sd_m1*real(sd_ncol,kind=r8) + sd_m2 ) ** O_THRD
               sd_r2 = sd_r1
               do k=1,22

                  s = idx_nasl(k)
                  dtmp = sd_asl1(s)*real(sd_ncol,kind=RP) + sd_asl2(s)
                  sd_asl1(s) = sd_asl1(s) + dmask(k)*(dtmp-sd_asl1(s))
                  sd_asl2(s) = sd_asl1(s)

               end do

               !! invalid by collisions between SDs with
               !! same muliplicity

               if( sd_n1==0 ) then
                  sd_rk1 = INVALID
               end if

            end if

            !! check muliplicity

            if( sd_n1>(2.0_RP**63._RP) .or. sd_n2>(2.0_RP**63._RP) ) then
               iexced = -1
               cycle
            end if

            if( sd_n(icp(tc)) > sd_n(icp(tp)) ) then
               sd_n( icp(tc) )  = sd_n1
               sd_r( icp(tc) )  = sd_r1
               sd_rk( icp(tc) ) = sd_rk1

               sd_n( icp(tp) )  = sd_n2
               sd_r( icp(tp) )  = sd_r2

               do k=1,22
                  s = idx_nasl(k)
                  sd_asl( icp(tc),s ) = sd_asl1(s)
                  sd_asl( icp(tp),s ) = sd_asl2(s)
               end do

            else

               sd_n( icp(tp) )  = sd_n1
               sd_r( icp(tp) )  = sd_r1
               sd_rk( icp(tp) ) = sd_rk1

               sd_n( icp(tc) )  = sd_n2
               sd_r( icp(tc) )  = sd_r2

               do k=1,22
                  s = idx_nasl(k)
                  sd_asl( icp(tp),s ) = sd_asl1(s)
                  sd_asl( icp(tc),s ) = sd_asl2(s)
               end do

            end if

         end do

      end do

      if( iexced<0 ) then
        call PRC_MPIstop
      end if

     ! Deallocate
      deallocate( fsort_tag  )
      deallocate( fsort_freq )

    return
  end subroutine sdm_coales
  !-----------------------------------------------------------------------------
!  subroutine sdm_sd2momnt(ni,nj,nk,zph_crs,mu_sdm,mv_sdm,mw_sdm,  &
!                          sd_num,sd_n,sd_x,sd_y,sd_rk,            &
!                          sd_u,sd_v,sd_wc,sd_r,ilist)
!     use scale_const, only: &
!        rw => CONST_DWATR
!     use scale_grid, only: &
!        CX => GRID_CX, &
!        CY => GRID_CY
!     ! Input variables
!     integer, intent(in) :: ni  ! Model dimension in x direction
!     integer, intent(in) :: nj  ! Model dimension in y direction
!     integer, intent(in) :: nk  ! Model dimension in z direction
!     real(RP),intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
!     integer, intent(in) :: sd_num  ! number of super-droplets
!     integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
!     real(RP),intent(in) :: sd_x(1:sd_num)    ! x-coordinate of super-droplets
!     real(RP), intent(in) :: sd_y(1:sd_num)    ! y-coordinate of super-droplets
!     real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
!     real(RP), intent(in) :: sd_u(1:sd_num)    ! x-direction velocity of super-droplets
!     real(RP), intent(in) :: sd_v(1:sd_num)    ! y-direction velocity of super-droplets
!     real(RP), intent(in) :: sd_wc(1:sd_num)   ! zeta components of contravariant velocity of super-droplets
!     real(RP), intent(in) :: sd_r(1:sd_num)    ! equivalent radius of super-droplets
!     ! Output variables
!     real(RP), intent(out) :: mu_sdm(KA,IA,JA) ! a total of momentum of super-droplets in x-direction
!     real(RP), intent(out) :: mv_sdm(KA,IA,JA) ! a total of momentum of super-droplets in y-direction
!     real(RP), intent(out) :: mw_sdm(KA,IA,JA) ! a total of momentum of super-droplets in z-direction
!     integer,  intent(out) :: ilist(1:int(sd_num/nomp),1:nomp) ! buffer for list vectorization
!     ! Internal shared variables
!!     real(RP) :: dx_sdm              ! dx_sdm
!!     real(RP) :: dy_sdm              ! dy_sdm
!     ! Work variables for OpenMP
!     integer :: sd_str           ! index of divided loop by OpenMP
!     integer :: sd_end           ! index of divided loop by OpenMP
!     integer :: np               ! index for OpenMP
!     ! Work variables
!     real(RP) :: dcoef(KA,IA,JA)      ! coef.
!     real(RP) :: dist       ! temporary
!     real(RP) :: dtmp       ! temporary
!!     integer(DP) :: idx_sdm ! integer of 'dx_sdm'
!!     integer(DP) :: idy_sdm ! integer of 'dy_sdm'
!     integer :: nlist(1:nomp)    ! list number
!     integer :: cnt              ! counter
!     integer :: i, j, k, m, n    ! index
!     integer :: ix, jy           ! index
!    !---------------------------------------------------------------------
!
!      ! Initialize
!      do k = KS, KE
!      do i = IS, IE
!      do j = JS, JE
!       dcoef(k,i,j) = rw * F_THRD*ONE_PI / real(dx_sdm(i)*dy_sdm(i),kind=RP)
!!       idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!!       idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )
!      enddo
!      enddo
!      enddo
!
!      do np=1,nomp
!         nlist(np) = 0
!      end do
!
!      do k=1,nk
!      do j=0,nj+1
!      do i=0,ni+1
!         mu_sdm(k,i,j) = 0.0_RP
!         mv_sdm(k,i,j) = 0.0_RP
!         mw_sdm(k,i,j) = 0.0_RP
!      end do
!      end do
!      end do
!
!      ! Get index list for compressing buffer.
!      do np=1,nomp
!
!         sd_str = int(sd_num/nomp)*(np-1) + 1
!         sd_end = int(sd_num/nomp)*np
!         cnt    = 0
!         do n=sd_str,sd_end
!             if( sd_rk(n)>VALID2INVALID ) then
!               cnt = cnt + 1
!               ilist(cnt,np) = n
!            end if
!         end do
!         nlist(np) = cnt
!      end do
!
!      ! Get a total of momentum of super-droplets.
!
!      !### count momentum of super-droplets ###!
!      do np=1,nomp
!
!         if( nlist(np)>0 ) then
!            do m=1,nlist(np)
!               n = ilist(m,np)
!
!               dtmp = real(sd_n(n),kind=DP) * sd_r(n)*sd_r(n)*sd_r(n)
!
!               !=== momentum in x-direction ===!
!
!!               dist = sd_x(n) + 0.50_RP * real(dx_sdm,kind=RP)
!               dist = sd_x(n) + 0.50_RP * real(dx_sdm,kind=RP)
!               do ix = IS, IE
!!                if( sd_x(n) < CX(ix) ) then
!                if( dist < CX(ix) ) then
!                 i = ix
!                 exit
!                endif
!               enddo
!               do jy = JS, JE
!                if( sd_y(n) < CY(jy) ) then
!                 j = jy
!                 exit
!                endif
!               enddo
!!               i = int( floor( dist*1.d5   , kind=i8 )/idx_sdm ) + 2
!!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
!               k = floor( sd_rk(n) )
!
!               mu_sdm(i,j,k) = mu_sdm(i,j,k) + sd_u(n) * dtmp
!
!               !=== momentum in y-direction ===!
!
!               dist = sd_y(n) + 0.50_RP * real(dy_sdm,kind=RP)
!               do ix = IS, IE
!                if( sd_x(n) < CX(ix) ) then
!                 i = ix
!                 exit
!                endif
!               enddo
!               do jy = JS, JE
!!                if( sd_y(n) < CY(jy) ) then
!                if( dist < CY(jy) ) then
!                 j = jy
!                 exit
!                endif
!               enddo
!!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!!               j = int( floor( dist*1.d5   , kind=i8 )/idy_sdm ) + 2
!               k = floor( sd_rk(n) )
!
!               mv_sdm(k,i,j) = mv_sdm(k,i,j) + sd_v(n) * dtmp
!
!               !=== momentum in zeta-direction ===!
!
!               do ix = IS, IE
!                if( sd_x(n) < CX(ix) ) then
!                 i = ix
!                 exit
!                endif
!               enddo
!               do jy = JS, JE
!                if( sd_y(n) < CY(jy) ) then
!                 j = jy
!                 exit
!                endif
!               enddo
!!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
!               k = floor( sd_rk(n) + 0.50_RP )
!               mw_sdm(k,i,j) = mw_sdm(k,i,j) + sd_wc(n) * dtmp
!
!            end do
!
!         end if
!
!      end do
!
!      ! Convert super-droplets to density of liquid-water.
!      do k=2,nk-1
!      do j=1,nj-1
!      do i=1,ni-1
!
!         dtmp = dcoef(k,i,j)/real(zph_crs(k+1,i,j)-zph_crs(k,i,j),kind=RP)
!         mu_sdm(k,i,j) = mu_sdm(k,i,j) * dtmp
!         mv_sdm(k,i,j) = mv_sdm(k,i,j) * dtmp
!         mw_sdm(k,i,j) = mw_sdm(k,i,j) * dtmp
!
!      end do
!      end do
!      end do
!
!    return
!  end subroutine sdm_sd2momnt
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine sdm_sd2prec(dtb_crs,                        &
                         prec,sd_num,sd_n,sd_x,sd_y,     &
                         sd_r,sd_ri,sd_rj,sd_rk,ilist,pr_sdm)

!!$      use scale_grid, only: &
!!$          FX => GRID_FX, &
!!$          FY => GRID_FY
      use m_sdm_coordtrans, only: &
           sdm_x2ri, sdm_y2rj

      ! Input variables
      real(RP), intent(in) :: dtb_crs
      integer, intent(in) :: sd_num  ! number of super-droplets
      integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
      real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
      real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
      ! Input and output variables
      real(RP), intent(out) :: sd_ri(1:sd_num)   ! index-i(real) of super-droplets
      real(RP), intent(out) :: sd_rj(1:sd_num)   ! index-j(real) of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
      real(RP), intent(inout) :: prec(IA,JA,1:2) ! precipitation and accumlation
      ! Output variables
!      real(RP), intent(out) :: pr_sdm(KA,IA,JA) ! temporary buffer of CReSS dimension
      real(RP), intent(out) :: pr_sdm(1:IA,1:JA) ! temporary buffer of CReSS dimension
!      integer, intent(out) :: ilist(1:int(sd_num/nomp),1:nomp) ! buffer for list vectorization
      integer, intent(out) :: ilist(1:sd_num) ! buffer for list vectorization
      ! Work variables for OpenMP
      integer :: sd_str        ! index of divided loop by OpenMP
      integer :: sd_end        ! index of divided loop by OpenMP
      integer :: np            ! index for OpenMP
      ! Work variables
      real(RP) :: dcoef(IA,JA) ! coef.
      real(RP) :: dtmp         ! temporary variables
      real(RP) :: dtbiv        ! 1.e0 / time step

!      integer(kind=i8) :: idx_sdm   ! integer of 'dx_sdm'
!      integer(kind=i8) :: idy_sdm   ! integer of 'dy_sdm'

      integer :: nlist(1:nomp)      ! list number
      integer :: tlist              ! total list number
      integer :: cnt                ! counter

      integer :: i, j, m, n, ix, jy ! index
    !-------------------------------------------------------------------

    ! Initialize
      dtbiv = 1.0d0 / dtb_crs
      dcoef(1:IA,1:JA)=0.0d0
     do i = IS, IE
     do j = JS, JE
!!$      dcoef(i,j) = F_THRD * ONE_PI / real(dx_sdm(i)*dy_sdm(j),kind=RP)
      dcoef(i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j)
!      idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!      idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )
     enddo
     enddo

!!$      tlist = 0
!!$
!!$      do np=1,nomp
!!$         nlist(np) = 0
!!$      end do

!      do j=0,nj+1
!      do i=0,ni+1
      do j=1,JA
      do i=1,IA
         pr_sdm(i,j) = 0.d0
      end do
      end do

      ! Get index list for compressing buffer.
      cnt=0
      do n=1,sd_num
         if( sd_rk(n)<VALID2INVALID .and. &
              sd_rk(n)>PREC2INVALID ) then
            cnt = cnt + 1
            ilist(cnt) = n
         end if
      end do
!!$      do np=1,nomp
!!$
!!$         sd_str = int(sd_num/nomp)*(np-1) + 1
!!$         sd_end = int(sd_num/nomp)*np
!!$         cnt    = 0
!!$
!!$         do n=sd_str,sd_end
!!$            if( sd_rk(n)<VALID2INVALID .and. &
!!$                sd_rk(n)>PREC2INVALID ) then
!!$               cnt = cnt + 1
!!$               ilist(cnt,np) = n
!!$            end if
!!$         end do
!!$
!!$         nlist(np) = cnt
!!$         tlist     = cnt
!!$
!!$      end do

      ! Get precipitation
      if( cnt>0 ) then
         !### get horizontal face index(real) of super-droplets ###!
         call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
         call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

         do m=1,cnt
            n=ilist(m)
            i=floor(sd_ri(n))+1
            j=floor(sd_rj(n))+1

            pr_sdm(i,j) = pr_sdm(i,j)                         &
                         + sd_r(n) * sd_r(n) * sd_r(n)           &
                         * real(sd_n(n),kind=RP)

            sd_rk(n) = INVALID     !! convert to invalid
         end do

!!$      if( tlist>0 ) then
!!$
!!$         !## count voulme of super-droplets ###!
!!$         do np=1,nomp
!!$
!!$            if( nlist(np)>0 ) then
!!$
!!$               do m=1,nlist(np)
!!$
!!$                  n = ilist(m,np)
!!$
!!$!                  i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm ) + 2
!!$!                  j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm ) + 2
!!$                  do ix = IS, IE
!!$                   if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
!!$                    i = ix
!!$                    exit
!!$                   endif
!!$                  enddo
!!$                  do jy = JS, JE
!!$                   if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
!!$                    j = jy
!!$                    exit
!!$                   endif
!!$                  enddo
!!$                  pr_sdm(i,j) = pr_sdm(i,j)                         &
!!$                                + sd_r(n) * sd_r(n) * sd_r(n)           &
!!$                                * real(sd_n(n),kind=RP)
!!$
!!$                  sd_rk(n) = INVALID     !! convert to invalid
!!$
!!$               end do
!!$
!!$            end if
!!$
!!$         end do

         !### convert super-droplets to precipitation ###!

!         do j=0,nj+1
!         do i=0,ni+1
         do j=1,JA
         do i=1,IA

            dtmp = real( pr_sdm(i,j) * dcoef(i,j) )

            !! rain fall rate
            prec(i,j,1) = dtmp * dtbiv

            !! accumulation
            prec(i,j,2) = prec(i,j,2) + dtmp

         end do
         end do

      end if

    return
  end subroutine sdm_sd2prec
  !----------------------------------------------------------------------------
!  subroutine sdm_commucrs(ni,nj,nk,val)
!    ! Input variables
!
!    integer, intent(in) :: ni    ! Model dimension in x direction
!    integer, intent(in) :: nj    ! Model dimension in y direction
!    integer, intent(in) :: nk    ! Model dimension in z direction
!
!    ! Input and output variable
!
!    real(RP), intent(inout) ::  val(0:ni+1,0:nj+1,1:nk)
!                                   ! Optional scalar variable

!-----7--------------------------------------------------------------7--

! Exchange the value horizontally.

      !### Exchange the value horizontally between sub domain ###

      !== x direction ==!

!      call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,val,1,1,sbuf)
!
!      call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)
!
!      call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,val,1,1,rbuf)
!
!      !== y direction ==!
!
!      call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,val,1,1,sbuf)
!
!      call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)
!
!      call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,val,1,1,rbuf)


      !### Exchange the value horizontally between grouped domain ###!

!      !== x direction ==!

!     call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,val,1,1,sbuf)
!
!      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)
!
!      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,val,1,1,rbuf)

      !== y direction ==!

!      call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,val,1,1,sbuf)
!
!      call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)
!
!      call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,val,1,1,rbuf)

      !== x direction again ( for L-corner ) ==!
!
!      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,val,1,1,sbuf)
!
!      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)
!
!      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,val,1,1,rbuf)


!      !### Set the periodic boundary conditions ###!

!      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
!     &            2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,val)


      !### Set the boundary conditions at the four corners ###!

!      call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,ni,nj,nk,val)

!    return
!  end subroutine sdm_commucrs
  !----------------------------------------------------------------------------
  subroutine sdm_aslform(DENS,RHOT,QTRC,                          &   
                         sdm_calvar,sdm_aslset,                   &
                         sdm_aslfmsdnc,sdm_sdnmlvol,              &
                         sdm_zupper,sdm_zlower,dtcl,              &
                         pbr_crs,ptbr_crs,pp_crs,         &
                         ptp_crs,qv_crs,zph_crs,rhod_crs,         &
                         sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_z,    &
                         sd_rk,sd_u,sd_v,sd_vz,sd_r,sd_asl,       &
                         sd_fmnum,sd_fmn,sd_fmx,sd_fmy,sd_fmz,    &
                         sd_fmri,sd_fmrj,sd_fmrk,sd_fmvz,sd_fmr,sd_fmasl,         &
                         ni_sdm,nj_sdm,nk_sdm,sort_id,sort_key,   &
                         sort_freq,sort_tag,sd_rng,               &
!                         sort_freq,sort_tag,                      &
                         sort_tag0,sd_itmp1,sd_itmp2,sd_itmp3)

    use scale_tracer, only: &
         QAD => QA
    use m_sdm_coordtrans, only: sdm_z2rk
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod
    use m_sdm_motion, only: &
         sdm_getvz

      ! Input variables
    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(in) :: QTRC(KA,IA,JA,QAD)    !! Ratio of mass of tracer to total mass[kg/kg]
      logical, intent(in) :: sdm_calvar(3) ! Control flag of calculation using SDM
      integer, intent(in) :: sdm_aslset    ! Option for aerosol species
      real(RP),intent(in) :: sdm_sdnmlvol  ! Normal volume for number concentration of super droplets
      real(RP),intent(in) :: sdm_aslfmsdnc ! Number of super droplets at aerosol formation per sdm_sdnmlvol
      real(RP),intent(in) :: sdm_zlower    ! Lower limitaion of initial droplet's position
      real(RP),intent(in) :: sdm_zupper    ! Upper limitaion of initial droplet's position
      real(RP),intent(in) :: dtcl(1:5)     ! Time interval of cloud micro physics
!      real(RP), intent(in) :: jcb_crs(KA,IA,JA)   ! Jacobian at scalar points
      real(RP), intent(in) :: pbr_crs(KA,IA,JA)   ! Base state pressure
      real(RP), intent(in) :: ptbr_crs(KA,IA,JA)  ! Base state potential temperature
      real(RP), intent(in) :: pp_crs(KA,IA,JA)    ! Pressure perturbation
      real(RP), intent(in) :: ptp_crs(KA,IA,JA)   ! Potential temperature perturbation
      real(RP), intent(in) :: qv_crs(KA,IA,JA)    ! Water vapor mixing ratio
      real(RP), intent(in) :: zph_crs(KA,IA,JA)   ! Z-physical coordinates
      real(RP), intent(in) :: rhod_crs(KA,IA,jA)   ! dry air densiy
      integer, intent(in) :: ni_sdm   ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm   ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm   ! SDM model dimension in z direction
      integer, intent(in) :: sd_num   ! number of super-droplets
      integer, intent(in) :: sd_numasl! number of kind of chemical material contained as water-soluble aerosol in super droplets
      integer, intent(in) :: sd_fmnum ! number of super-droplets at aerosol formation
      ! Input and output variables
      integer(DP), intent(inout) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
      real(RP), intent(inout) :: sd_z(1:sd_num)      ! z-coordinate of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
      real(RP), intent(inout) :: sd_u(1:sd_num)      ! x-components velocity of super-droplets
      real(RP), intent(inout) :: sd_v(1:sd_num)      ! y-components velocity of super-droplets
      real(RP), intent(inout) :: sd_vz(1:sd_num)     ! terminal velocity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)      ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
      integer(DP), intent(inout) :: sd_fmn(1:sd_fmnum) ! multiplicity of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmx(1:sd_fmnum)  ! x-coordinate of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmy(1:sd_fmnum)  ! y-coordinate of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmz(1:sd_fmnum)  ! z-coordinate of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmri(1:sd_fmnum) ! index[i/real] of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmrj(1:sd_fmnum) ! index[j/real] of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmrk(1:sd_fmnum) ! index[k/real] of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmvz(1:sd_fmnum) ! terminal velocity of super-droplet at aerosol formation
      real(RP), intent(inout) :: sd_fmr(1:sd_fmnum)  ! equivalent radius of super-droplets at aerosol formation
      real(RP), intent(inout) :: sd_fmasl(1:sd_fmnum,1:sd_numasl)  ! aerosol mass of super-droplets at aerosol formation
      type(c_rng_uniform_mt), intent(inout) :: sd_rng     ! random number generator
      integer, intent(inout) :: sort_id(1:sd_num)         ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num)        ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1) ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)  ! accumulated number of super-droplets in each SD-grid
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)   ! = sort_tag(n) - 1
      integer, intent(out) :: sd_itmp1(1:sd_num,1:nomp)  ! temporary array of the size of the number of super-droplets.
      integer, intent(out) :: sd_itmp2(1:sd_num,1:nomp)  ! temporary array of the size of the number of super-droplets.
      integer, intent(out) :: sd_itmp3(1:sd_num,1:nomp)  ! temporary array of the size of the number of super-droplets.
      ! Work variables
      real(RP) :: sd_fmnc      ! number concentration of super-droplets at aerosol formation
      real(RP) :: sdm_dtevl    ! time step of {condensation/evaporation} process
      real(RP) :: sdm_aslfmdt  ! time step of aerosol formation
      integer :: max_key       ! max key to sort data
      integer :: n_invd        ! number of invalid droplets
      integer :: innum         ! temporary
      integer :: id_invd       ! index
      integer :: k, n

      real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
      real(RP) :: rhod_scale(KA,IA,JA) ! dry air density
      real(RP) :: t_scale(KA,IA,JA)    ! Temperature
    !---------------------------------------------------------------------

      ! Check active

      if( abs(sdm_aslset)<10 ) return
      ! Sorting valid super-droplets

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_x,sd_y,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      max_key = ni_sdm * nj_sdm * knum_sdm + 1

      n_invd = sort_freq(max_key)    !! number of invalid droplets

      do n=1,max_key+1
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      ! Initialize

      sdm_dtevl   = real( dtcl(1), kind=RP )
      sdm_aslfmdt = real( dtcl(4), kind=RP )

      ! Get random number

      call gen_rand_array( sd_rng, sd_fmx  )
      call gen_rand_array( sd_rng, sd_fmy  )
      call gen_rand_array( sd_rng, sd_fmz  )
      call gen_rand_array( sd_rng, sd_fmvz )

      do k=1,sd_numasl

         call gen_rand_array( sd_rng, sd_fmr )

         do n=1,sd_fmnum
            sd_fmasl(n,k) = sd_fmr(n)   !! sd_fmr : temporary
         end do

      end do

      ! Set aerosol mass and muliplicity of super-droplets

      !### Ammonium Sulfate [(NH4)2SO4] ###!

      sd_fmnc = real(sdm_aslfmsdnc,kind=RP)/real(sdm_sdnmlvol,kind=RP)
                                     !! number concentration of super-
                                     !! -droplets at aerosol formation

      call sdm_aslsulf(sdm_aslfmrate,sdm_aslfmdt,                &
                       sd_fmnum,sd_numasl,sd_fmn,sd_fmasl,sd_fmnc)

      ! Set position of super-droplets
      do n=1,sd_fmnum

         sd_fmx(n) = xmax_sdm * sd_fmx(n)
         sd_fmy(n) = ymax_sdm * sd_fmy(n)
         sd_fmz(n) = real(minzph+sdm_zlower,kind=RP) + sd_fmz(n)        &
                   * real(sdm_zupper-(minzph+sdm_zlower),kind=RP)

         sd_fmr(n) = 1.0E-15_RP

      end do

      !### index[k/real] of super-droplets ###!

      call sdm_z2rk(sdm_zlower,sdm_zupper,       &
                        sd_fmnum,sd_fmx,sd_fmy,sd_fmz,sd_fmri,sd_fmrj,sd_fmrk)

      ! Set radius of super-droplets
      if( sdm_calvar(1) ) then

         call sdm_condevp(sdm_aslset,                              &
                          sdm_aslmw,sdm_aslion,sdm_dtevl,           &
                          pbr_crs,ptbr_crs,pp_crs,ptp_crs,          &
                          qv_crs,                                       &
                          sd_fmnum,sd_numasl,sd_fmx,sd_fmy,sd_fmr,      &
                          sd_fmasl,sd_fmrk)

      end if

     ! Set terminal velocity of super-droplets

      call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
      call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

      call sdm_getvz(pres_scale,rhod_scale,t_scale,                         &
                     sd_fmnum,sd_fmx,sd_fmy,sd_fmri,sd_fmrj,sd_fmrk,sd_fmr,sd_fmvz,   &
                     sd_itmp1(1:sd_fmnum,1),sd_itmp2(1:sd_fmnum,1),   &
                     sd_itmp3(1:sd_fmnum,1),'motion')

      ! Add new droplets as formed aerosol

      !### adjust number of super-droplets at aerosol formation ###!

      if( n_invd<sd_fmnum ) then

        if( IO_L ) then
         write(IO_FID_LOG,'(2a)')"  ### [SDM] : warning for buffer size limit",  &
          " @ stop to form super-droplets as aerosol ###"
        endif

         innum = nomp * int(n_invd/nomp)

      else

         innum = sd_fmnum

      end if

      do n=1,innum

         id_invd = sort_id( sort_tag0(max_key) + n )

         sd_n(id_invd)  = sd_fmn(n)

         sd_x(id_invd)  = sd_fmx(n)
         sd_y(id_invd)  = sd_fmy(n)
         sd_z(id_invd)  = sd_fmz(n)
         sd_rk(id_invd) = sd_fmrk(n)

         sd_u(id_invd)  = 0.0_RP
         sd_v(id_invd)  = 0.0_RP
         sd_vz(id_invd) = sd_fmvz(n)

         sd_r(id_invd)  = sd_fmr(n)

      end do

      do k=1,sd_numasl

         do n=1,innum

            id_invd = sort_id( sort_tag0(max_key) + n )

            sd_asl(id_invd,k) = sd_fmasl(n,k)

         end do

      end do


    return
  end subroutine sdm_aslform
  !----------------------------------------------------------------------------
  subroutine sdm_sort(ni_sdm,nj_sdm,nk_sdm,                       &
                      sd_num,sd_n,sd_x,sd_y,sd_rk,                &
                      sort_id,sort_key,sort_freq,sort_tag,jdgtype)
!      use m_gadg_algorithm, only: &
!          gadg_count_sort
      use scale_grid, only: &
          FX => GRID_FX, &
          FY => GRID_FY
      ! Input variables
      integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
      integer, intent(in) :: sd_num  ! number of super-droplets
      integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
      real(RP), intent(in) :: sd_x(1:sd_num)   ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num)   ! y-coordinate of super-droplets
      real(RP), intent(in) :: sd_rk(1:sd_num)  ! index[k/real] of super-droplets
      character(len=5), intent(in) :: jdgtype   ! flag for sorting
      ! Output variables
      integer, intent(out) :: sort_id(1:sd_num)   ! id that super-droplets sorted by sd-grids
      integer, intent(out) :: sort_key(1:sd_num)  ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each sd-grid
      integer, intent(out) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)     ! accumulated number of super-droplets in each sd-grid
      ! Work variables
!      integer(kind=i8) :: idx_sdm   ! integer of 'dx_sdm'
!      integer(kind=i8) :: idy_sdm   ! integer of 'dy_sdm'
      integer :: adr                ! grid id
      integer :: max_key            ! total grid number
      integer :: i, j, k, n         ! index
      integer :: ix, jy
     !--------------------------------------------------------------------

     ! Initialize

      max_key = ni_sdm * nj_sdm * knum_sdm + 1

!      idx_sdm = 10_i8 * floor( 1.e4*(dx_sdm+1.e-5), kind=i8 )
!      idy_sdm = 10_i8 * floor( 1.e4*(dy_sdm+1.e-5), kind=i8 )

      ! Sorting [step-1] -- initialize IKEY

      if( jdgtype .eq. 'valid' ) then

         do n=1,sd_num

            if( sd_rk(n)>VALID2INVALID ) then

               !! get index of super-droplets in physical domain
               !! of CReSS
!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm )
!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm )
               do ix = IS, IE
                   if( sd_x(n) <= ( FX(ix)-FX(IS-1) ) ) then
                    i = ix
                    exit
                   endif
               enddo
               do jy = JS, JE
                   if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                    j = jy
                    exit
                   endif
               enddo
!               k = floor( sd_rk(n) ) - 2
               k = floor( sd_rk(n) ) - KS + 1

               !=========================================
               !--- k-index of SCALE start from 1 but
               !--- that of CRess start from 0
               !=========================================
!               adr = ni_sdm*nj_sdm*k &
               adr = ni_sdm*nj_sdm*(k-1) &
                   + ni_sdm*(j-JS+1) + (i-IS+1) + 1

            else

               adr = max_key    !! invalid super-droplets

            end if

            sort_key(n) = adr
!if( sort_key(n) > max_key .or. sort_key(n) < 1 ) then
!write(*,*) sort_key(n), k, i, j, sd_rk(n), sdz_s2c(n)
!endif
         end do

      else if( jdgtype .eq. 'multi' ) then

         do n=1,sd_num

            if( sd_rk(n)>VALID2INVALID .and. sd_n(n)>1 ) then

               !! get index of super-droplets in physical domain
               !! of CReSS

!               i = int( floor( sd_x(n)*1.d5, kind=i8 )/idx_sdm )
!               j = int( floor( sd_y(n)*1.d5, kind=i8 )/idy_sdm )
               do ix = IS, IE
                   if( sd_x(n) <= ( FX(ix)-FX(IS-1)) ) then
                    i = ix
                    exit
                   endif
               enddo
               do jy = JS, JE
                   if( sd_y(n) <= ( FY(jy)-FY(JS-1) ) ) then
                    j = jy
                    exit
                   endif
               enddo
!               k = floor( sd_rk(n) ) - 2
               k = floor( sd_rk(n) ) - KS + 1

               !=========================================
               !--- k-index of SCALE start from 1 but
               !--- that of CRess start from 0
               !=========================================
!               adr = ni_sdm*nj_sdm*k +ni_sdm*j+i+1
               adr = ni_sdm*nj_sdm*(k-1) &
                   + ni_sdm*(j-JS+1) + (i-IS+1) + 1

            else

               adr = max_key    !! invalid super-droplets

            end if

            sort_key(n) = adr

         end do

      end if

!do n = 1, max_key
! write(102,*) n, sort_freq(n), sort_id(n), sort_key(n), size(sort_key),max_key-1+1
!enddo

      ! Sorting [step-2] -- counting sort
      call gadg_count_sort( sort_key, 1, max_key,                       &
                            sort_freq, sort_tag, sort_id )

!do n = 1, max_key
! write(103,*) n, sort_freq(n), sort_id(n)
!enddo

      sort_tag(max_key+1) = sort_tag(max_key) + sort_freq(max_key)

    return
  end subroutine sdm_sort
  !----------------------------------------------------------------------------
  subroutine sdm_adjsdnum(sdm_nadjvar,ni_sdm,nj_sdm,nk_sdm,    &
                          sd_num,sd_numasl,sd_nc,                 &
                          sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_rk, &
                          sort_id,sort_key,sort_freq,sort_tag,    &
                          sd_rng,sd_rand,                         &
!                          sd_rand,                                &
                          sdm_itmp1,sdm_itmp2,sd_itmp1)
      use scale_process, only: &
           mype => PRC_myrank
      ! Input variables
      integer, intent(in) :: sdm_nadjvar
      integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
      integer, intent(in) :: sd_num  ! number of super-droplets
      integer, intent(in) :: sd_numasl  ! number of kind of chemical material contained as water-soluble aerosol in super droplets
      real(RP), intent(in) :: sd_nc ! averaged number concentration in a grid
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng  ! random number generator
      real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num)  ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num) ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)   ! accumulated number of super-droplets in each SD-grid
      real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
      integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
      real(RP), intent(inout) :: sd_vz(1:sd_num) ! terminal velocity of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sdm_itmp1(1:ni_sdm*nj_sdm*nk_sdm+2)  ! temporary array of SDM dimension
      integer, intent(out) :: sdm_itmp2(1:ni_sdm*nj_sdm*nk_sdm+2)  ! temporary array of SDM dimension
      integer, intent(out) :: sd_itmp1(1:sd_num,1:nomp)  ! temporary array of the size of the number of super-droplets.
      ! Parameter
      real(RP), parameter :: RATE4REMOVE = 1.40_RP ! upper limit rate to averaged number concentration in a grid for removing
      real(RP), parameter :: RATE4ADD  = 0.70_RP  ! lower limit rate to averaged number concentration in a grid for adding
      ! Work variables
      integer :: sdnum_upr ! upper limit number of super-droplets in a grid for removing
      integer :: sdnum_lwr ! lower limit number of super-droplets in a grid for adding
      !------------------------------------------------------------------7--

      if( sdm_nadjvar==0 ) return

      ! Remove super-droplets from grids with large number of super-droplets

      if( sdm_nadjvar==3 .or. sdm_nadjvar==2 ) then

         sdnum_upr = floor( RATE4REMOVE * sd_nc )
         call sdm_sdremove(ni_sdm,nj_sdm,nk_sdm,                        &
                           sdnum_upr,sd_num,sd_n,sd_x,sd_y,sd_rk,       &
                           sort_id,sort_key,sort_freq,sort_tag,         &
                           sd_rng,sd_rand,                              &
!                           sd_rand,                                     &
                           sdm_itmp1,sdm_itmp2,sd_itmp1)

      end if

      ! Add super-droplets to grids with small number of super-droplets

      if( sdm_nadjvar==3 .or. sdm_nadjvar==1 ) then

         sdnum_lwr = floor( RATE4ADD * sd_nc )

         call sdm_sdadd(ni_sdm,nj_sdm,nk_sdm,                           &
                        sdnum_lwr,sd_num,sd_numasl,                     &
                        sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_rk,         &
                        sort_id,sort_key,sort_freq,sort_tag,            &
                        sd_rng,sd_rand,                                 &
!                        sd_rand,                                        &
                        sdm_itmp1,sdm_itmp2,sd_itmp1)

      end if

     ! Message for adjustment

      if( mype==0 .and. sdm_nadjvar==3 ) then
        if( IO_L ) then
         write(IO_FID_LOG,'(a,a,i4,a,i4,a)')                           &
                 "  ### [SDM] : adjust number of super-droplet ",      &
                       "( min",sdnum_lwr," -- max",sdnum_upr," ) ###"
        endif
      else if( mype==0 .and. sdm_nadjvar==1 ) then

        if( IO_L ) then
         write(IO_FID_LOG,'(a,a,i4,a)')                                &
                 "  ### [SDM] : adjust number of super-droplet ",      &
                           "( min",sdnum_lwr," -- max INFINITY ) ###"
        endif

      else if( mype==0 .and. sdm_nadjvar==2 ) then
        if( IO_L ) then
         write(IO_FID_LOG,'(a,a,i4,a)')                                &
                 "  ### [SDM] : adjust number of super-droplet ",      &
                                  "( min 0 -- max",sdnum_upr," ) ###"
        endif

      end if


    return
  end subroutine sdm_adjsdnum
  !----------------------------------------------------------------------------
  subroutine sdm_sdremove(ni_sdm,nj_sdm,nk_sdm,                   &
                          sdnum_upr,sd_num,sd_n,sd_x,sd_y,sd_rk,  &
                          sort_id,sort_key,sort_freq,sort_tag,    &
                          sd_rng,sd_rand,                         &
!                          sd_rand,                                &
                          sort_tag0,fsort_id,isd_perm)

!      use m_gadg_algorithm, only: &
!          gadg_count_sort
      ! Input variables
      integer, intent(in) :: ni_sdm  ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm  ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm  ! SDM model dimension in z direction
      integer, intent(in) :: sdnum_upr ! upper limit number of super-droplets in a grid  for removing
      integer, intent(in) :: sd_num    ! number of super-droplets
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng    ! random number generator
      real(RP), intent(inout) :: sd_rand(1:sd_num)  ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num)        ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num)       ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)  ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)   ! accumulated number of super-droplets in each SD-grid
      integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)  ! = sort_tag(n) - 1
      integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
      integer, intent(out) :: isd_perm(1:sd_num,1:nomp)   ! random permutations
      ! Internal shared variables
      integer :: freq_max   ! get the maximum number of super-droplets in each grid
      ! Work variables
      real(RP) :: drate   ! temporary
      integer, allocatable :: fsort_tag(:)  ! buffer for sorting
      integer, allocatable :: fsort_freq(:) ! buffer for sorting
      integer :: gnum          ! grid number
      integer :: id_vd         ! index of valid super-droplets
      integer :: n_minus       ! number of reducing super-droplets
      integer :: sdnum_valid   ! number of valid super-droplets in each grid
      integer :: ip, m, n, s, t            ! index
     !---------------------------------------------------------------------

      ! Initialize

      gnum = ni_sdm * nj_sdm * knum_sdm

      freq_max = 1

      ! Sorting valid super-droplets

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_x,sd_y,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      ! Initialize
      do n=1,gnum+2
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      ! Get the maximum number of super-droplets in each grid
      do m=1,gnum
         freq_max = max( freq_max, sort_freq(m) )
      end do


!write(*,*) sdnum_upr, freq_max, sort_freq(gnum+1)
      ! Sorting the grids by the number of super-droplets in each grid

      allocate( fsort_tag(0:freq_max+1) )
      allocate( fsort_freq(0:freq_max)  )

      call gadg_count_sort( sort_freq(1:gnum), 0, freq_max,             &
                            fsort_freq, fsort_tag, fsort_id )

      fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

      ! Get random number using random number generator

      call gen_rand_array( sd_rng, sd_rand )

      ! Get random permutation layout of super-droplets in each grid

      call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
                       sort_tag0,fsort_tag,fsort_id,sd_rand,isd_perm)

      ! Remove super-droplets from grids with large number of super-droplets

      do s=1,freq_max
         do t=fsort_tag(sdnum_upr+1),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            !### get number of removable droplets in each grid ###!
            sdnum_valid = sort_freq(m)   !! number of valid droplets

            n_minus = min( int(sdnum_valid/2),                          &
                           max(sdnum_valid-sdnum_upr,0) )

            if( s>sdnum_valid ) cycle    !! cycle at more than number of valid droplets

            !### get index of valid droplets ###!
            ip = isd_perm( sort_tag0(m)  + s ,1)   !! search forward
!SELECT     ip = isd_perm( sort_tag(m+1) - s )   !! search backward

            id_vd = sort_id( sort_tag0(m) + ip )

            !### Remove droplets ###!

            if( s<=n_minus ) then

               !== convert valid droplets to invalid droplet ==!

               sd_rk(id_vd) = INVALID

            else if( s>n_minus .and. s<=sdnum_valid ) then

               !== adjust multiplicity of valid droplets ==!

               drate = real(sdnum_valid,kind=RP)                        &
                                      /real(sdnum_valid-n_minus,kind=RP)

               sd_n(id_vd) = nint( real(sd_n(id_vd),kind=DP)*drate      &
                                                            , kind=DP )

            end if

         end do

      end do


     ! Deallocate

      deallocate( fsort_tag  )
      deallocate( fsort_freq )

    return
  end subroutine sdm_sdremove
  !----------------------------------------------------------------------------
  subroutine sdm_sdadd(ni_sdm,nj_sdm,nk_sdm,                      &
                       sdnum_lwr,sd_num,sd_numasl,                &
                       sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_rk,    &
                       sort_id,sort_key,sort_freq,sort_tag,       &
                       sd_rng,sd_rand,                            &
!                       sd_rand,                                   &
                       sort_tag0,fsort_id,isd_perm)
!      use m_gadg_algorithm, only: &
!          gadg_count_sort
      ! Input variables
      integer, intent(in) :: ni_sdm ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm ! SDM model dimension in z direction
      integer, intent(in) :: sdnum_lwr  ! lower limit number of super-droplets in a grid for adding
      integer, intent(in) :: sd_num ! number of super-droplets
      integer, intent(in) :: sd_numasl  ! number of kind of chemical material contained as water-soluble aerosol in super droplets
      ! Input and output variables
      type(c_rng_uniform_mt), intent(inout) :: sd_rng   ! random number generator
      real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
      integer, intent(inout) :: sort_id(1:sd_num)       ! super-droplets sorted by SD-grids
      integer, intent(inout) :: sort_key(1:sd_num)      ! sort key
      integer, intent(inout) :: sort_freq(1:ni_sdm*nj_sdm*nk_sdm+1)   ! number of super-droplets in each SD-grid
      integer, intent(inout) :: sort_tag(1:ni_sdm*nj_sdm*nk_sdm+2)    ! accumulated number of super-droplets in each SD-grid
      real(RP), intent(inout) :: sd_x(1:sd_num)    ! x-coordinate of super-droplets
      real(RP), intent(inout) :: sd_y(1:sd_num)    ! y-coordinate of super-droplets
      integer(DP), intent(inout) :: sd_n(1:sd_num) ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_r(1:sd_num)    ! equivalent radius of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
      real(RP), intent(inout) :: sd_vz(1:sd_num)   ! terminal velocity of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
      ! Output variables
      integer, intent(out) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2) ! = sort_tag(n) - 1
      integer, intent(out) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
      integer, intent(out) :: isd_perm(1:sd_num,1:nomp)  ! random permutations
      ! Internal shared variables
      integer :: freq_max ! get the maximum number of super-droplets in each grid
      ! Work variables
      integer, allocatable :: fsort_tag(:)  ! buffer for sorting
      integer, allocatable :: fsort_freq(:) ! buffer for sorting
      integer, allocatable :: sort_freq_valid(:) ! number of valid super-droplets in each grid
      integer :: idx_nasl(1:20) ! index for vactorization
      integer :: gnum          ! grid number
      integer :: id_invd       ! index of invalid super-droplets
      integer :: id_vd         ! index of valid super-droplets
      integer :: iwarn         ! warning flag
      integer :: n_invd        ! number of invalid super-droplets
      integer :: n_plus        ! number of adding super-droplets
      integer :: sdnum_valid   ! number of valid super-droplets in each grid
      integer :: sdnum_split   ! number of splitable valid super-droplets in each grid
      integer :: ip            ! index
      integer :: m, n, q, s, t ! index
     !---------------------------------------------------------------------

      ! Initialize
      gnum = ni_sdm * nj_sdm * knum_sdm
      freq_max = 1

      ! Add super-droplets in the grid with small number of super-droplets

      ! Sorting super-droplets

      !### valid droplets ###!

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_x,sd_y,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      !### valid droplets that multiplicity is over 2 ###!

      allocate( sort_freq_valid(1:ni_sdm*nj_sdm*nk_sdm+1) )

      do n=1,ni_sdm*nj_sdm*nk_sdm+1
         sort_freq_valid(n) = sort_freq(n)
      end do

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_x,sd_y,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'multi')

      ! Initialize
      do n=1,20

         if( n<=sd_numasl ) then
            idx_nasl(n) = n
         else
            idx_nasl(n) = sd_numasl
         end if

      end do

      do n=1,gnum+2
         sort_tag0(n) = sort_tag(n) - 1  !! {1-xx} => {0-xx}
      end do

      ! Get the maximum number of super-droplets in each grid

      do m=1,gnum
         freq_max = max( freq_max, sort_freq(m) )
      end do

      ! Sorting the grids by the number of super-droplets in each grid
      allocate( fsort_tag(0:freq_max+1) )
      allocate( fsort_freq(0:freq_max)  )

      call gadg_count_sort( sort_freq(1:gnum), 0, freq_max,             &
                            fsort_freq, fsort_tag, fsort_id )

      fsort_tag(freq_max+1) = fsort_tag(freq_max) + fsort_freq(freq_max)

      ! Get random number using random number generator

      call gen_rand_array( sd_rng, sd_rand )

      ! Get random permutation layout of super-droplets in each grid

      call sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,            &
                       sort_tag0,fsort_tag,fsort_id,sd_rand,isd_perm(1:sd_num,1))

      ! Search invalid super-droplets
      n_invd = sort_freq(gnum+1)   !! number of invalid droplet

      ! Add super-droplets to grids with small number of super-droplets

      q = 1        !! counter
      iwarn = -1   !! warning flg

      do s=1,int(sdnum_lwr/2)

         do t=fsort_tag(1),fsort_tag(sdnum_lwr)-1

            m = fsort_id(t)

            !### get number of addable droplets in each grid ###!

            sdnum_valid = sort_freq_valid(m)
                                        !! number of valid droplets
            sdnum_split = sort_freq(m)
                                        !! number of valid droplets that
                                        !! multiplicity is over 2
                                        !! ( splitable droplets )

            n_plus = min( sdnum_split, max(sdnum_lwr-sdnum_valid,0) )

            if( s>n_plus ) cycle        !! cycle at more than number
                                        !! of addable droplets

            !### get random permutation ###!
            ! 'sdm_getperm' target grid within two or more droplets
            ! ip = 1    : sort_freq(m)=1
            ! ip = perm : sort_freq(m)>1
            ip = min(sort_freq(m),2) - 1                 !! 0 or 1

            ip = (1-ip) + isd_perm(sort_tag0(m) +s,1)*ip   !! forward
!SELECT     ip = (1-ip) + isd_perm(sort_tag(m+1)+s)*ip   !! backward

            !### get index of valid and invalid droplets ###!

            id_vd   = sort_id( sort_tag0(m) + ip )
            id_invd = sort_id( sort_tag0(gnum+1) + min(q,n_invd) )

            !### Add droplets ###!

            !== split multiplicity of valid droplets ==!
            !== and copy other status to invalid one ==!

            sd_n(id_invd)  = sd_n(id_vd)/2
            sd_n(id_vd)    = sd_n(id_vd) - sd_n(id_invd)

            sd_x(id_invd)  = sd_x(id_vd)
            sd_y(id_invd)  = sd_y(id_vd)
            sd_r(id_invd)  = sd_r(id_vd)
            sd_vz(id_invd) = sd_vz(id_vd)
            sd_rk(id_invd) = sd_rk(id_vd)

            do n=1,20
              sd_asl(id_invd,idx_nasl(n)) = sd_asl(id_vd,idx_nasl(n))
            end do

            !### count and check number of invaid droplets ###!

            q = q + 1

            if( q>n_invd ) then
               iwarn = 1
            end if

         end do

      end do

      if( iwarn==1 ) then
        if( IO_L ) then
         write(IO_FID_LOG,'(2a)')"  ### [SDM] : warning for stopping to add",    &
               " super-droplets to avoid buffer limit exceeded ###"
        else
         write(*,'(2a)')"  ### [SDM] : warning for stopping to add",    &
               " super-droplets to avoid buffer limit exceeded ###"
        endif
      end if

      ! Deallocate

      deallocate( sort_freq_valid )

      deallocate( fsort_tag  )
      deallocate( fsort_freq )

    return
  end subroutine sdm_sdadd
  !----------------------------------------------------------------------------
  subroutine sdm_getperm(freq_max,ni_sdm,nj_sdm,nk_sdm,sd_num,    &
                         sort_tag0,fsort_tag,fsort_id,            &
                         sd_rand,sd_perm)

     ! Input variables
      integer, intent(in) :: freq_max ! maximum number of SD in each sd-grid
      integer, intent(in) :: ni_sdm   ! SDM model dimension in x direction
      integer, intent(in) :: nj_sdm   ! SDM model dimension in y direction
      integer, intent(in) :: nk_sdm   ! SDM model dimension in z direction
      integer, intent(in) :: sd_num   ! Number of super-droplets
      integer, intent(in) :: sort_tag0(1:ni_sdm*nj_sdm*nk_sdm+2)
                           ! = sort_tag(n) - 1
                           ! sort_tag(m) : accumulated number of
                           !               super-droplets in each
                           !               SDM-grid
      integer, intent(in) :: fsort_tag(0:freq_max+1)
                           ! accumulated number with respect to
                           ! the number contaiend super-droplets
      integer, intent(in) :: fsort_id(1:ni_sdm*nj_sdm*nk_sdm+2)
                           ! super-droplets sorted by the number
                           ! contaiend super-droplets
      ! Input and output variables
      real(RP), intent(inout) :: sd_rand(1:sd_num) ! random numbers
      ! Output variables
      integer, intent(out) :: sd_perm(1:sd_num)    ! random permutations
     ! Work variables
      integer :: m, n, t, ss, tt         ! index
      !---------------------------------------------------------------

      ! Initialize for calculating random permutations

      do t=fsort_tag(2),fsort_tag(freq_max+1)-1

         m = fsort_id(t)

         tt = sort_tag0(m) + 1
         sd_perm(tt) = 1

      end do

      ! Calculate random permutations

      do n=2,freq_max

         do t=fsort_tag(n),fsort_tag(freq_max+1)-1

            m = fsort_id(t)

            tt = sort_tag0(m) + n
!ORG        ss = sort_tag0(m) + ( int(sd_rand(tt)*n) + 1 )
            ss = sort_tag0(m) + ( int(sd_rand(tt)*real(n,kind=RP)) + 1 )

            !### swap data ###!

            sd_perm(tt) = sd_perm(ss)
            sd_perm(ss) = n

         end do

      end do

    return
  end subroutine sdm_getperm
  !----------------------------------------------------------------------------
  subroutine sdm_aslsulf(sd_aslfmrate,sd_aslfmdt,             &
                         sd_num,sd_numasl,sd_n,sd_asl,sd_fmnc)

      use scale_process, only: &
         PRC_MPIstop
      ! Input variables
      real(RP), intent(in) :: sd_aslfmrate   ! formation rate of aerosol
      real(RP), intent(in) :: sd_aslfmdt     ! time interval to form aerosol
      integer,  intent(in) :: sd_num         ! number of super-droplets
      integer,  intent(in) :: sd_numasl      ! number of kind of chemical material contained as water-soluble aerosol in super droplets
      real(RP), intent(in) :: sd_fmnc        ! number concentration of super-droplets at aerosol formation
      ! Input and output variables
      integer(DP), intent(inout) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
      real(RP), intent(inout) :: sd_asl(1:sd_num,1:sd_numasl)  ! aerosol mass of super-droplets
      ! Work variables
      real(DP) :: n_amsul
      real(RP) :: rb_amsul
      real(RP) :: rmax_amsul
      real(RP) :: rmin_amsul
      real(RP) :: sgm_amsul   ! parameter for 3mode log-nomiral distribution
      real(RP) :: n0          ! number of real droplets per unit volume and per aerosol radius
      real(RP) :: dry_r       ! aerosol radius
      real(RP) :: delta       ! temporary
      real(RP) :: sdn_tmp     ! temporary
      integer :: iexced       ! temporary
      integer :: n, t         ! index
     !-------------------------------------------------------------------

      ! Set ammonium sulfate aerosol mass and muliplicity of super-droplets

      iexced = 1

      rb_amsul   = rb_amsul_zanten_m
      sgm_amsul  = sgm_amsul_zanten_m
      rmax_amsul = rmax_amsul_zanten_m
      rmin_amsul = rmin_amsul_zanten_m

      n_amsul = real(sd_aslfmrate,RP) * real(sd_aslfmdt,RP)

      do n=1,sd_num

         !### 1 : (NH4)2SO4 in modified vanZanten(2010) ###!
         delta = log(rmax_amsul) - log(rmin_amsul)
         dry_r = exp( log(rmin_amsul) + delta*sd_asl(n,1) )
         sd_asl(n,1) = F_THRD * ONE_PI                             &
                              * (dry_r*dry_r*dry_r) * rho_amsul

         !! n0(log(dry_r)) [m-3]
         !! 1st mode log-noraml distribution for ammonium sulfate
         delta = log(dry_r) - log(rb_amsul)
         delta = -(delta*delta)/(2.0_RP*log(sgm_amsul)*log(sgm_amsul))

         n0 = (n_amsul*exp(delta))/(sqrt(2.0_RP*ONE_PI)*log(sgm_amsul))

         !! continuous uniform distribution

         delta = 1.0_RP/(log(rmax_amsul)-log(rmin_amsul))

         !! muliplicity

         sdn_tmp = n0/(sd_fmnc*delta)

         if( sdn_tmp<(2.0_RP**63_RP) ) then
            sd_n(n) = nint( sdn_tmp, DP )
         else
            iexced = -1
         end if

      end do

      if( iexced<0 ) then
         call PRC_MPIstop
      end if

     ! Set other aerosol mass

      do t=2,sd_numasl
      do n=1,sd_num
         sd_asl(n,t) = 0.0_RP
      end do
      end do

    return
  end subroutine sdm_aslsulf
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sdm_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0  )
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QA) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
!
!    real(RP) :: dens
!    real(RP) :: q(QA)
!    real(RP) :: RLMDr, RLMDs, RLMDg

!    real(RP) :: zerosw
!    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

!    Re(:,:,:,I_mp_QC) =   8.E-6_RP
!    Re(:,:,:,I_mp_QI) =  20.E-6_RP

    ! Effective radius is defined by r3m/r2m=1.5/lambda.
!    do j  = JS, JE
!    do i  = IS, IE
!    do k  = KS, KE
!       ! store to work
!       dens = DENS0(k,i,j)
!       do iq = I_QV, I_QG
!          q(iq) = QTRC0(k,i,j,iq)
!       enddo

       ! slope parameter lambda
!       zerosw = 0.5_RP - sign(0.5_RP, q(I_QR) - 1.E-12_RP )
!       RLMDr = sqrt(sqrt( dens * q(I_QR) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )
!
!       zerosw = 0.5_RP - sign(0.5_RP, q(I_QS) - 1.E-12_RP )
!       RLMDs = sqrt(sqrt( dens * q(I_QS) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )
!
!       zerosw = 0.5_RP - sign(0.5_RP, q(I_QG) - 1.E-12_RP )
!       RLMDg = sqrt(sqrt( dens * q(I_QG) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )
!
!       Re(k,i,j,I_mp_QR) = RLMDr * 1.5_RP
!       Re(k,i,j,I_mp_QS) = RLMDs * 1.5_RP
!       Re(k,i,j,I_mp_QG) = RLMDg * 1.5_RP
!    enddo
!    enddo
!    enddo
!
    return
  end subroutine ATMOS_PHY_MP_sdm_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sdm_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)
!
!    real(RP) :: qhydro
!    integer  :: k, i, j, iq
!    !---------------------------------------------------------------------------
!
!    do j  = JS, JE
!    do i  = IS, IE
!    do k  = KS, KE
!       qhydro = 0.D0
!       do iq = 1, MP_QA
!          qhydro = qhydro + QTRC(k,i,j,I_MP2ALL(iq))
!       enddo
!       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
!    enddo
!    enddo
!    enddo
!
    return
  end subroutine ATMOS_PHY_MP_sdm_CloudFraction
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_Mixingratio( &
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

    integer  :: ihydro
    !---------------------------------------------------------------------------

    do ihydro = 1, MP_QA
       Qe(:,:,:,ihydro) = QTRC0(:,:,:,I_MP2ALL(ihydro))
    enddo

    return
  end subroutine ATMOS_PHY_MP_sdm_Mixingratio
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_in
    implicit none

    integer :: sdnum_dum, sdnumasl_dum, sdfmnum_dum
    integer :: dp_dum, rp_dum, n, m, ierr
    real(DP) :: otime

    !### Get random generator seed ###!
    !! Random number generator has already been initialized in scale-les/src/preprocess/mod_mkinit.f90
    !! Be careful. If unit (=fid_random_i) is specified, filename is ignored and the object is initialized by the unit.
    call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME))
!    call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME), fid_random_i )

    open (fid_sd_i, file = trim(SD_IN_BASENAME), &! action = "read", &
          access = "sequential", status = "old", form = "unformatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_restart_in", "read error"
      call PRC_MPIstop
    endif 
    
    !--- read time and precision
    read(fid_sd_i) otime, rp_dum, dp_dum, sdnum_dum, sdnumasl_dum, sdfmnum_dum
    if( IO_L ) write(IO_FID_LOG,*) '*** Input restart file of Super Droplet  '
    if( IO_L ) write(IO_FID_LOG,*) 'SD. restart now time =  ', real(otime,kind=DP)
    if( rp_dum /= RP .or. dp_dum /= DP .or. &
        sdnum_dum /= sdnum_s2c .or. sdnumasl_dum /= sdnumasl_s2c .or. &
        sdfmnum_dum /= sdfmnum_s2c ) then
       write(*,*) 'xxx RP, DP, sdnum_s2c, sdnumasl_s2c, sdfmnum_s2c in'
       write(*,*) 'is different from those in Param file!  stop'
       write(*,*) 'RP(in restart file) = ', rp_dum
       write(*,*) 'DP(in restart file) = ', dp_dum
       write(*,*) 'sdnum(in restart file) = ', sdnum_dum
       write(*,*) 'RP(in restart file) = ', sdnumasl_dum
       write(*,*) 'RP(in restart file) = ', sdfmnum_dum
       call PRC_MPIstop
    endif

    !--- read S.D.
    read(fid_sd_i) (sdn_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdrk_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdx_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdy_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdz_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdr_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdu_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdv_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) (sdvz_s2c(n),n=1,sdnum_s2c)
    read(fid_sd_i) ((sdasl_s2c(n,m),n=1,sdnum_s2c),m=1,sdnumasl_s2c)
    !--- read formation S.D.
    close(fid_sd_i)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed restart file of Super Droplet  '

    return
  end subroutine ATMOS_PHY_MP_sdm_restart_in
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_out(otime)
    implicit none

    real(DP), intent(in) :: otime
    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp, ftmp2
    character(len=H_LONG) :: basename_sd_out
    character(len=H_LONG) :: basename_random
    character(len=H_LONG) :: basename_time
    integer :: n, m, ierr


    !--- output restart file of Super Droplet
    write(basename_time,'(F15.3)') otime
    do n = 1, 15
       if( basename_time(n:n) == ' ' ) basename_time(n:n) = '0'
    enddo
    fid_sd_o = IO_get_available_fid()

    if( SD_OUT_BASENAME == '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found S.D. restart file name. Default used..'
       write(fmt2(14:14),'(I1)') 6
       write(fmt2(16:16),'(I1)') 6
       write(ftmp,fmt3) 'SD_output', '_', trim(basename_time)
!       write(SD_OUT_BASENAME,fmt2) trim(ftmp), 'pe',mype ! Perhaps a bug?
       write(basename_sd_out,fmt2) trim(ftmp),'pe',mype
    else
       write(fmt2(14:14),'(I1)') 6
       write(fmt2(16:16),'(I1)') 6
       write(ftmp,fmt3) trim(SD_OUT_BASENAME), '_', trim(basename_time)
       write(basename_sd_out,fmt2) trim(ftmp),'pe',mype
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Output restart file of Super Droplet  ', trim(basename_sd_out)

    open (fid_sd_o, file = trim(basename_sd_out), & !action = "write", &
          access = "sequential", status = "replace", form = "unformatted", &
          iostat = ierr)

    if( ierr /= 0 ) then
      write(*,*) "sdm_restart_out", "Write error"
      call PRC_MPIstop
    endif 

    !--- write time and precision
    write(fid_sd_o) otime, RP, DP, sdnum_s2c, sdnumasl_s2c, sdfmnum_s2c
    !--- write S.D.
    write(fid_sd_o) (sdn_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdrk_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdx_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdy_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdz_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdr_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdu_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdv_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) (sdvz_s2c(n),n=1,sdnum_s2c)
    write(fid_sd_o) ((sdasl_s2c(n,m),n=1,sdnum_s2c),m=1,sdnumasl_s2c)

    close(fid_sd_o)
    if( IO_L ) write(IO_FID_LOG,*) '*** Closed restart file of Super Droplet  '

    !--- output restart file of Random number
    if( RANDOM_OUT_BASENAME == '' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found random number output file name. Default used..'
       write(ftmp,fmt3) 'random_number_output', '_', trim(basename_time)
       write(basename_random,fmt2) trim(ftmp),'pe',mype
    else
       write(ftmp,fmt3) trim(RANDOM_OUT_BASENAME), '_', trim(basename_time)
       write(basename_random,fmt2) trim(ftmp),'pe',mype
    endif

    fid_random_o = IO_get_available_fid()
    if( IO_L ) write(IO_FID_LOG,*) '*** Output random number for SDM ***', trim(basename_random)
    call rng_save_state( rng_s2c, trim(basename_random))
!    call rng_save_state( rng_s2c, trim(basename_random), fid_random_o )

    return
  end subroutine ATMOS_PHY_MP_sdm_restart_out
end module scale_atmos_phy_mp_sdm
!-------------------------------------------------------------------------------
