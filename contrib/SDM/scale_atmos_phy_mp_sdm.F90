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
!! @li      2014-07-12 (S.Shima) [rev] sdm_sort repaired
!! @li      2014-07-12 (S.Shima) [rev] sdm_coales repaired
!! @li      2014-07-12 (S.Shima) [rev] sdm_coales tuned for FX (use compile option -Kprefetch_indirect)
!! @li      2014-07-12 (S.Shima) [rev] sdm_sort and sdm_getperm are separated into the module m_sdm_idutil
!! @li      2014-07-12 (S.Shima) [rev] sdm_coales is separated into the module m_sdm_coalescence
!! @li      2014-07-13 (S.Shima) [rev] sdm_condevp repaired
!! @li      2014-07-13 (S.Shima) [rev] sdm_sd2qcqr, sdm_sd2rhocr repaired
!! @li      2014-07-14 (S.Shima) [rev] sdm_condevp_updatefluid created modifying and repairing sdm_sd2qvptp
!! @li      2014-07-14 (S.Shima) [rev] Update HALO after call sdm_condevp_updatefluid
!! @li      2014-07-14 (S.Shima) [rev] sdm_condevp tuned for FX (use compile option -Kocl -Kprefetch_indirect)
!! @li      2014-07-14 (S.Shima) [rev] sdm_sd2rhow, sdm_sd2rhocr, sdm_sd2qcqr are separated into m_sdm_sd2fluid
!! @li      2014-07-14 (S.Shima) [rev] sdm_condevp,sdm_condevp_updatefluid are separated into m_sdm_condensation_water
!! @li      2014-07-14 (S.Shima) [rev] Removed unused subroutines, variables
!! @li      2014-07-18 (Y.Sato)  [mod] Modify the modules to calculate variable forradiation
!! @li      2014-07-18 (Y.Sato)  [add] Add history output for QTRC_sdm and input 0 to QTRC at the end of ATMOS_PHY_MP_sdm
!! @li      2014-07-18 (Y.Sato)  [mod] Modify a way to determine whether the grid stretch or not
!! @li      2014-07-18 (Y.Sato)  [add] Add ATMOS_PHY_MP_DENS, which is used in radiation code
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the way to determine whether the grid stretch or not
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the way to determine whether the lateral boundary is periodic or not
!! @li      2014-07-22 (Y.Sato)  [mod] Modify to write error information to normal output
!! @li      2014-07-22 (Y.Sato)  [mod] Modify bugs to calculate variable for radiation process
!! @li      2014-07-22 (Y.Sato)  [mod] Modify timing for assigning to QTRC_sdm and clear QTRC(2:QA)
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the definition of kl and ku for calculating drate
!! @li      2014-07-22 (Y.Sato)  [mod] Modify the order of loop in a part
!! @li      2014-07-24 (Y.Sato)  [mod] Modify a bug for restart
!! @li      2014-07-25 (Y.Sato)  [rev] Move sdm_getrklu from sdm_iniset to ATMOS_PHY_MP_sdm_setup
!! @li      2014-07-25 (Y.Sato)  [rev] Add COMM_var, and COMM_wait for filling u_scale, v_scale, and w_scale
!! @li      2014-12-12 (Y.Sato)  [mod] Modify for using QTRC_sdm in sdm_sd2qcqr in sdm_iniset 
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
         GRID_CDZ,   &
!!$       CZ  => GRID_CZ,    &
!!$       FZ  => GRID_FZ,    &
!!$       FDX => GRID_FDX,   &
!!$       FDY => GRID_FDY,   &
!!$       FDZ => GRID_FDZ,   & 
       GRID_FX,    &
       GRID_FY,    &
!!       CBFZ => GRID_CBFZ, &
!!       CBFX => GRID_CBFX, &
!!       CBFY => GRID_CBFY, &
       ISG, IEG, JSG, JEG,&
       DX,DY,DZ
    use scale_topography, only: &
       TOPO_Zsfc
    use m_sdm_coordtrans, only: &
       sdm_getrklu
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE

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
       write(*,*) 'xxx ATMOS_TYPE_PHY_MP is not SDM. Check!'
       call PRC_MPIstop
    endif

    buffact = 0.0_RP
    do k = KS, KE
      buffact = max( buffact,GRID_CDZ(k)/DZ )
    enddo
    do j = JS, JE
      buffact = max( buffact,GRID_CDY(j)/DY )
    enddo
    do i = IS, IE
      buffact = max( buffact,GRID_CDX(i)/DX )
    enddo

    if( buffact < 1.0_RP .or. buffact > 1.0_RP ) then
       sthopt = 1
    else
       sthopt = 0
    endif

    if(sthopt==1) then
       write(*,*) 'ERROR: stretched coordinate is not yet supported!', buffact
       call PRC_MPIstop
    end if

    if( maxval( TOPO_Zsfc ) > 0.0_RP ) then
       trnopt = 2
    elseif( maxval( TOPO_Zsfc ) == 0.0_RP ) then
       trnopt = 0
    endif

    if(trnopt==2) then
       write(*,*) 'ERROR: terrain following coordinate is not yet supported!'
       call PRC_MPIstop
    end if

    !--- if the lateral boundary is not periodic, PRC_next has negative value
    if( minval( PRC_next,1 ) < 0 ) then
       write(*,*) 'ERROR: Only periodic B.C. is supported!'
       write(*,*) 'ERROR: Set PRC_PERIODIC_X=PRC_PERIODIC_Y=.true.'
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

    if( sdm_zupper > minval(REAL_FZ(KE-1,IS:IE,JS:JE)) ) then
       if( mype == PRC_master )  write(*,*) "sdm_zupper was set to minval(REAL_FZ(KE-1)) because zupper > minval(REAL_FZ(KE-1))"
       sdm_zupper = minval(REAL_FZ(KE-1,IS:IE,JS:JE))
    endif

     sdm_dtevl = real( sdm_dtcmph(1),kind=RP )  !! condensation/evaporation
     sdm_dtcol = real( sdm_dtcmph(2),kind=RP )  !! stochastic coalescence
     sdm_dtadv = real( sdm_dtcmph(3),kind=RP )  !! motion of super-droplets

    ! check whether sdm_dtcmph(1:3) > 0
     if(  ( (sdm_dtcmph(1) <= 0.0_RP) .and. docondensation   )   .or. &
         ( (sdm_dtcmph(2) <= 0.0_RP) .and. doautoconversion )   .or. &
         ( (sdm_dtcmph(3) <= 0.0_RP) .and. domovement       )        ) then
       write(*,*) 'ERROR: sdm_dtcmph(1:3) have to be positive'
       call PRC_MPIstop
     end if

    ! check whether dt (dt of mp) is larger than sdm_dtcmph(i).
    if( (dt < sdm_dtcmph(1)) .or. (dt < sdm_dtcmph(2)) .or. (dt < sdm_dtcmph(3)) ) then
       write(*,*) 'ERROR: For now, sdm_dtcmph should be smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_MPIstop
    end if

    ! aerosol nucleation and sd number adjustment functions are not supported yet
    if ( (abs(sdm_aslset) >= 10) .or. (sdm_nadjvar /= 0)) then
       write(*,*) 'ERROR: aerosol nucleation and sd number adjustment functions are not supported yet'
       write(*,*) 'ERROR: set sdm_aslset < 10 and sdm_nadjvar =0'
       call PRC_MPIstop
    end if

    ! rigorous momentum exchange function is not supported yet
    if ( sdm_mvexchg /= 0) then
       write(*,*) 'ERROR: Momentum exchange not yet supported. set sdm_mvexchg = 0'
       call PRC_MPIstop
    end if

    if( docondensation ) then
       nclstp(1)=10*int(1.E+2_RP*(dt+0.0010_RP))            &
            /int(1.E+3_RP*(sdm_dtcmph(1)+0.00010_RP))
       if(mod(10*int(1.E+2_RP*(dt+0.0010_RP)),int(1.E+3_RP*(sdm_dtcmph(1)+0.00010_RP))) /= 0) then 
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(1) = 1
    end if

    if( doautoconversion ) then
       nclstp(2)=10*int(1.E+2_RP*(dt+0.0010_RP))            &
            /int(1.E+3_RP*(sdm_dtcmph(2)+0.00010_RP))
       if(mod(10*int(1.E+2_RP*(dt+0.0010_RP)),int(1.E+3_RP*(sdm_dtcmph(2)+0.00010_RP))) /= 0) then 
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_MPIstop
       end if
    else
       nclstp(2) = 1
    end if

    if( domovement ) then
       nclstp(3)=10*int(1.E+2_RP*(dt+0.0010_RP))            &
            /int(1.E+3_RP*(sdm_dtcmph(3)+0.00010_RP))
       if(mod(10*int(1.E+2_RP*(dt+0.0010_RP)),int(1.E+3_RP*(sdm_dtcmph(3)+0.00010_RP))) /= 0) then 
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
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
       write(*,*) 'ERROR: TIME_DTSEC_ATMOS_PHY_MP should be the least comon multiple of sdm_dtcmph(1:3) that are smaller than TIME_DTSEC_ATMOS_PHY_MP'
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

    !### Get index[k/real] at "sdm_zlower" and "sdm_zupper"  ###!
    call sdm_getrklu(sdm_zlower,sdm_zupper,      &
                     sdrkl_s2c,sdrku_s2c)

    dx_sdm(1:IA) = GRID_CDX(1:IA)
    dy_sdm(1:JA) = GRID_CDY(1:JA)
    dxiv_sdm(1:IA) = GRID_RCDX(1:IA)
    dyiv_sdm(1:JA) = GRID_RCDY(1:JA)

    ATMOS_PHY_MP_DENS(I_mp_QC) = dens_w ! hydrometeor density [kg/m3]=[g/L]

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
    use m_sdm_sd2fluid, only: &
         sdm_sd2qcqr
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
!!$    real(RP) :: pbr_crs(KA,IA,JA)  ! Base state pressure
!!$    real(RP) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
!!$    real(RP) :: ppf_crs(KA,IA,JA)  ! Pressure perturbation at future
!!$    real(RP) :: uf_crs(KA,IA,JA)   ! u components of velocity at future
!!$    real(RP) :: vf_crs(KA,IA,JA)   ! v components of velocity at future
!!$    real(RP) :: wf_crs(KA,IA,JA)   ! w components of velocity at future
!!$!    real(RP) :: wcf_crs(KA,IA,JA)  ! zeta components of contravariant velocity at future
!!$    real(RP) :: ptpf_crs(KA,IA,JA) ! Potential temperature perturbation at future
!!$    real(RP) :: qvf_crs(KA,IA,JA)  ! Water vapor mixing ratio at future
!!$    real(RP) :: prr_crs(IA,JA,1:2) ! Precipitation and accumulation for rain
    real(RP) :: rtmp4(KA,IA,JA)    ! Temporary array
    real(RP) :: rtmp5(KA,IA,JA)    ! Temporary array
    real(RP) :: rtmp6(KA,IA,JA)    ! Temporary array
    ! Output variables
!!$    real(RP) :: exnr_crs(KA,IA,JA) ! Exner function
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

    real(RP) :: pres_scale(KA,IA,JA) ! Pressure
    real(RP) :: rhod_scale(KA,IA,JA) ! dry air density
    real(RP) :: t_scale(KA,IA,JA)    ! Temperature
    real(RP) :: rhoc_scale(KA,IA,JA) ! cloud water density
    real(RP) :: rhor_scale(KA,IA,JA) ! rain water density
    real(RP) :: QHYD_sdm(KA,IA,JA)

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
       !! Evaluate diagnostic variables
       !!! z
       call sdm_rk2z(sdnum_s2c,sdx_s2c,sdy_s2c,sdrk_s2c,sdz_s2c,sdri_s2c,sdrj_s2c)
       !!! terminal velocity vz
       call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
       call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz(pres_scale,rhod_scale,t_scale,            &
                           sdnum_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c,sdvz_s2c,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'motion' )
       !! Output
       call sdm_outasci(TIME_NOWSEC,                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdm_dmpnskip)
    else if(sdm_dmpvar>1)then
       write(*,*) 'ERROR: sdm_dmpvar>1 not supported for now. Set sdm_dmpvar=1 (Output Super Droplet Data in ASCII)'
    end if

    !== run SDM at future ==!
     call sdm_calc(MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,                 & 
                   sdm_calvar,sdm_mvexchg,dtcl(1:3), sdm_aslset,  &
!!$                   exnr_crs,pbr_crs,ptbr_crs,         &
!!$                   ppf_crs,ptpf_crs,qvf_crs,&
!!$                   prr_crs,zph_crs,rhod_crs,                      &
                   prr_crs,zph_crs,                      &
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

     !== Diagnose QC and QR from super-droplets ==!
     !! note that when SDM is used QC := rhoc/(rhod+rhov), QR := rhor/(rhod+rhov)
     !! We need to rethink whether it's okay to evaluate QC QR here
     !! For example, when we diagnose pressure or rhod during the dynamical process, 
     !! qc and qr should be 0 because they are not included in the total density
     call sdm_sd2qcqr(DENS,QTRC_sdm(:,:,:,I_QC),QTRC_sdm(:,:,:,I_QR),        &
                      zph_crs,                                               &
                      sdnum_s2c,sdn_s2c,sdx_s2c,sdy_s2c,                     &
                      sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,&
                      rhoc_scale,rhor_scale,                                 &
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
!!$                         sdy_s2c,sdz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdu_s2c,           &
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
!!$                         sdasl_s2c,sdvz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,               &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,       &
!!$                         sorttag_s2c,rng_s2c,rand_s2c,                &
!!$!                         sorttag_s2c,rand_s2c,                      &
!!$                         sdm_itmp1,sdm_itmp2,sd_itmp1)
!!$
!!$      end if

    do k = KS, KE
    do i = IS, IE
    do j = JS, JE
        QHYD_sdm(k,i,j) = QTRC_sdm(k,i,j,I_QC)+QTRC_sdm(k,i,j,I_QR)
    enddo
    enddo
    enddo

!! Are these correct? Let's check later.
!!$    call HIST_in( rhoa_sdm(:,:,:), 'RAERO', 'aerosol mass conc.', 'kg/m3', dt)
    call HIST_in( prr_crs(:,:,1), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( prr_crs(:,:,1), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt)
    call HIST_in( QTRC_sdm(:,:,:,I_QC), 'QC_sd', 'mixing ratio of cloud in SDM', 'kg/kg', dt)
    call HIST_in( QTRC_sdm(:,:,:,I_QR), 'QR_sd', 'mixing ratio of rain in SDM', 'kg/kg', dt)
    call HIST_in( QHYD_sdm(:,:,:), 'QHYD_sd', 'mixing ratio of liquid in SDM', 'kg/kg', dt)

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
        sdm_z2rk
      use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod
      use m_sdm_motion, only: &
        sdm_getvz
      use m_sdm_sd2fluid, only: &
           sdm_sd2qcqr
      use m_sdm_condensation_water, only: &
           sdm_condevp

      real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
      real(RP), intent(in) :: RHOT(KA,IA,JA) ! DENS * POTT [K*kg/m3]
      real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD) ! ratio of mass of tracer to total mass[kg/kg]
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
!!$      real(RP) :: ptbr_crs(KA,IA,JA)  ! potential temperature
!!$      real(RP) :: ptp_crs(KA,IA,JA)   ! Potential temperature perturbation
!!$      real(RP) :: pbr_crs(KA,IA,JA)   ! pressure
!!$      real(RP) :: pp_crs(KA,IA,JA)    ! Pressure perturbation
!!$      real(RP) :: qv_crs(KA,IA,JA)    ! Water vapor mixing ratio
      real(RP) :: n0                            ! number of real droplets per unit volume and per aerosol radius
      real(RP) :: dry_r                         ! aerosol radius
      real(RP) :: delta1, delta2, sdn_tmp       ! temporary
      logical :: lsdmup                         ! flag for updating water hydrometeor by SDM
      integer :: iexced, sdnum_tmp1, sdnum_tmp2 ! temporary
      integer :: i, j, k, n, iq, np             ! index
!!$      real(RP) :: CPTOT(KA,IA,JA), RTOT(KA,IA,JA)
!!$      real(RP) :: QDRY(KA,IA,JA),  CPovCV(KA,IA,JA)
      real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
      integer :: sd_str, sd_end, sd_valid

      real(RP) :: pres_scale(KA,IA,JA)  ! Pressure
      real(RP) :: rhod_scale(KA,IA,JA) ! dry air density
      real(RP) :: t_scale(KA,IA,JA)    ! Temperature
      real(RP) :: rhoc_scale(KA,IA,JA) ! cloud water density
      real(RP) :: rhor_scale(KA,IA,JA) ! rain water density
     !---------------------------------------------------------------------

      ! conversion of SCALE variables to CReSS variables: ptbr ptp pbr pp rhod qv 
      ! This part will be omitted in the future.
!!$      CPTOT(:,:,:) = 0.0_RP
!!$      QDRY(:,:,:) = 1.0_RP
!!$      do k = 1, KA
!!$      do i = 1, IA
!!$      do j = 1, JA
!!$        ptbr_crs(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
!!$        do iq = QQS, QQE
!!$          QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
!!$        enddo
!!$        rhod_crs(k,i,j) = DENS(k,i,j) * QDRY(k,i,j) ! dry air density
!!$        RTOT (k,i,j) = Rdry * QDRY(k,i,j) + Rvap * QTRC(k,i,j,I_QV)
!!$        CPTOT(k,i,j) = CPdry * QDRY(k,i,j)
!!$        do iq = QQS, QQE
!!$         CPTOT(k,i,j) = CPTOT(k,i,j) + QTRC(k,i,j,iq) * CPw(iq)
!!$        enddo
!!$        CPovCV(k,i,j) = CPTOT(k,i,j) / ( CPTOT(k,i,j) - RTOT(k,i,j) )
!!$        pbr_crs(k,i,j) = P00 * ( RHOT(k,i,j) * RTOT(k,i,j) / P00 )**CPovCV(k,i,j)
!!$        pp_crs(k,i,j) = 0.0_RP
!!$        ptp_crs(k,i,j) = 0.0_RP
!!$        qv_crs(k,i,j) = QTRC(k,i,j,I_QV) ! Be careful. The definition of qv is different. (SCALE: qv=rhov/rho, CReSS: qv=rhov/rhod)
!!$      enddo
!!$      enddo
!!$      enddo

      if( .not. sdm_calvar(1) .and. &
          .not. sdm_calvar(2) .and. &
          .not. sdm_calvar(3)        ) return

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
         write(*,*) "sdm_iniset, exceeded"
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

      !### index[k/real] of super-droplets               ###!
      !### modify position[z] of invalid super-droplets  ###!
      call sdm_z2rk(sdm_zlower,sdm_zupper,            &
                        sdnum_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c )

      !### initial equivalent radius                     ###!
      do n=1,sdnum_s2c
         sdr_s2c(n) = 1.0E-15_RP

!ORG     sdr_s2c(n) = 1.0e-5 * ( log(1.0/(1.0-sdr_s2c(n))) )**O_THRD
!ORG     sdr_s2c(n) = 1.0d-8

! temporary for test
!         sdr_s2c(n) = 3.0E-3_RP*sdr_s2c(n)
!         sdr_s2c(n) = exp((log(3.0E-3_RP)-log(1.0E-7_RP))*sdr_s2c(n)+log(1.0E-7_RP))
         
      end do
      !### ( at only condensation/evaporation process ) ###!
      if( sdm_calvar(1) ) then

         call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

         call sdm_condevp(sdm_aslset,                                   &
                          sdm_aslmw,sdm_aslion,sdm_dtevl,               &
                          pres_scale,t_scale,QTRC(:,:,:,I_QV),          &
                          sdnum_s2c,sdnumasl_s2c,sdx_s2c,sdy_s2c,       &
                          sdr_s2c,sdasl_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c)

      end if

      !### diagnose terminal velocity (no need evaluate them here) ###!
!!$      call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
!!$      call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
!!$      call sdm_getvz(pres_scale,rhod_scale,t_scale,                    &
!!$                     sdnum_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c, &
!!$                     sdvz_s2c,sd_itmp1,sd_itmp2,sd_itmp3,'motion')

      !### Diagnose QC and QR from super-droplets ###!
      !! note that when SDM is used QC := rhoc/(rhod+rhov), QR := rhor/(rhod+rhov)
!!$      call sdm_sd2qcqr(DENS,QTRC(:,:,:,I_QC),QTRC(:,:,:,I_QR),          &
      call sdm_sd2qcqr(DENS,QTRC_sdm(:,:,:,I_QC),QTRC_sdm(:,:,:,I_QR),          &
                       zph_crs,                   &
                       sdnum_s2c,sdn_s2c,sdx_s2c,sdy_s2c,        &
                       sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
                       rhoc_scale,rhor_scale,                      &
                       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2            )

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
  subroutine sdm_calc(MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,              & 
                      sdm_calvar,sdm_mvexchg,dtcl, sdm_aslset,    &
!!$                      exnr_crs,                       &
!!$                      pbr_crs,ptbr_crs,pp_crs, &
!!$                      ptp_crs,qv_crs,prec_crs,zph_crs,rhod_crs,   &
                      prec_crs,zph_crs,   &
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
   use scale_comm, only: &
        COMM_vars8, &
        COMM_wait
   use m_sdm_fluidconv, only: &
        sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod
   use m_sdm_sd2fluid, only: &
        sdm_sd2rhow, sdm_sd2prec
   use m_sdm_boundary, only: &
        sdm_jdginvdv, sdm_boundary
    use m_sdm_motion, only: &
        sdm_getvz, sdm_getvel, sdm_move
    use m_sdm_coalescence, only: &
        sdm_coales
    use m_sdm_condensation_water, only: &
        sdm_condevp, sdm_condevp_updatefluid
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
!!$   real(RP), intent(in) :: exnr_crs(KA,IA,JA) ! exner function
!!$   real(RP), intent(in) :: pbr_crs(KA,IA,JA)  ! Base state pressure
!!$   real(RP), intent(in) :: ptbr_crs(KA,IA,JA) ! Base state potential temperature
!!$   real(RP), intent(in) :: pp_crs(KA,IA,JA)   ! Pressure perturbation
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
!!$   real(RP), intent(inout) :: rhod_crs(KA,IA,JA)  ! dry air density
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
!!$   real(RP), intent(inout) :: ptp_crs(KA,IA,JA)! Potential temperature perturbation
!!$   real(RP), intent(inout) :: qv_crs(KA,IA,JA) ! Water vapor mixing ratio
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
            !! cres_val1p is the rhow before condevp
            call sdm_sd2rhow(zph_crs,crs_val1p,                       &
                             sd_num,sd_n,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2p,sd_itmp1)

            ! { condensation/evaporation } in SDM
            !! diagnose necessary fluid variables
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! update the equivalent radius of SDs
            call sdm_condevp(sdm_aslset,            &
                             sdm_aslmw,sdm_aslion,sdm_dtevl,      &
                             pres_scale,t_scale,QTRC(:,:,:,I_QV), &
                             sd_num,sd_numasl,sd_x,sd_y,sd_r,sd_asl,&
                             sd_ri,sd_rj,sd_rk)

            ! get density of liquid-water(qw) after process-1
            !! here cres_val1c is the rhow after condevp
            call sdm_sd2rhow(zph_crs,crs_val1c,                       &
                             sd_num,sd_n,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,        &
                             sd_rkl,sd_rku,crs_val2c,sd_itmp1)

            ! exchange the vapor and heat to fluid variables
            call sdm_condevp_updatefluid(RHOT,QTRC,DENS,crs_val1p,crs_val1c)
            !! update the HALO region of the fluid variables
            call COMM_vars8( RHOT(:,:,:), 1 )
            call COMM_vars8( QTRC(:,:,:,I_QV), 2 )
            call COMM_vars8( DENS(:,:,:), 3 )
            call COMM_wait ( RHOT(:,:,:), 1 )
            call COMM_wait ( QTRC(:,:,:,I_QV), 2 )
            call COMM_wait ( DENS(:,:,:), 3 )

         end if
         !=== 2 : stochastic coalescence ===!
         if( sdm_calvar(2) .and.                                 &
             mod(t,istep_sdm/istep_col)==0 .and. sdm_dtcol>0.d0 ) then

            lsdmup = .true.

            ! get the terminal velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz(pres_scale,rhod_scale,t_scale,            &
                           sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'coales' )

            ! { coalescence } in SDM
            call sdm_coales(sdm_colkrnl,sdm_colbrwn,sdm_aslset,         &
                            sdm_aslrho,sdm_dtcol,                       &
                            pres_scale, t_scale,                        &
                            zph_crs,                                    &
                            ni_sdm,nj_sdm,nk_sdm,sd_num,sd_numasl,      &
                            sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,     &
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
              write(*,*) "sdm_mvexchg>1 cannot be applied"
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
            !! diagnose necessary fluid variables
            call sdm_rho_qtrc2rhod(DENS,QTRC,rhod_scale)
            call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)
            !! evaluate the terminal velocity
            call sdm_getvz(pres_scale,rhod_scale,t_scale,            &
                           sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'motion' )

            ! get the moving velocity of super-droplets
            !! diagnose necessary fluid variables
            call sdm_rho_mom2uvw(DENS,MOMX,MOMY,MOMZ,u_scale,v_scale,w_scale)
            !! update the HALO region of the fluid variables
            call COMM_vars8( u_scale(:,:,:), 1 )
            call COMM_vars8( v_scale(:,:,:), 2 )
            call COMM_vars8( w_scale(:,:,:), 3 )
            call COMM_wait ( u_scale(:,:,:), 1 )
            call COMM_wait ( v_scale(:,:,:), 2 )
            call COMM_wait ( w_scale(:,:,:), 3 )

            !! evaluate the velocity of each SD
            call sdm_getvel(u_scale,v_scale,w_scale,                   &
                            sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vz)

            ! get the momentum of super-droplets after moving.
            if( sdm_mvexchg>0 ) then
              write(*,*) "sdm_mvexchg>1 cannnot be applied"
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
            !! evaluate the motion eq.
            call sdm_move(sdm_dtadv,                         &
                          sd_num,sd_u,sd_v,sd_vz,sd_x,sd_y,sd_rk)

            ! lateral boundary routine in SDM
            !! judge super-droplets as invalid or valid in horizontal
            !! do MPI communication to send/receiv SDs
            call sdm_boundary(wbc,ebc,sbc,nbc,                           &
                             sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_rk,     &
                             sd_u,sd_v,sd_vz,sd_r,sd_asl,               &
                             bufsiz1,bufsiz2,sd_itmp1,rbuf,sbuf)

            ! judge super-droplets as invalid or valid in vartical
            call sdm_jdginvdv(sd_rkl,sd_rku,sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk)

            ! convert the impact of { exchange momentum } to
            ! {u,v,w} in CReSS
            if( sdm_mvexchg>0 ) then
              write(*,*) "sdm_mvexchg>1 cannot be applied"
              call PRC_MPIstop
!               call sdm_sd2uvw(ni,nj,nk,rst_crs,u_crs,v_crs,wc_crs,     &
!     &                         crs_val1p,crs_val2p,crs_val3p,           &
!     &                         crs_val1c,crs_val2c,crs_val3c)
!
            end if

         end if

      end do

    ! Convert super-droplets to precipitation

      call sdm_sd2prec(dt,                                &
                       prec_crs,                          &
                       sd_num,sd_n,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,  &
                       sd_itmp1(1:sd_num,1),crs_val1c(1,1:IA,1:JA))

    return
  end subroutine sdm_calc
!---------------------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  subroutine sdm_aslform(DENS,RHOT,QTRC,                          &   
                         sdm_calvar,sdm_aslset,                   &
                         sdm_aslfmsdnc,sdm_sdnmlvol,              &
                         sdm_zupper,sdm_zlower,dtcl,              &
                         pbr_crs,ptbr_crs,pp_crs,         &
                         ptp_crs,qv_crs,zph_crs,rhod_crs,         &
                         sd_num,sd_numasl,sd_n,sd_x,sd_y,sd_z,    &
                         sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vz,sd_r,sd_asl,       &
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
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort
    use m_sdm_condensation_water, only: &
        sdm_condevp
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
      real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
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

      call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

      ! Check active

      if( abs(sdm_aslset)<10 ) return
      ! Sorting valid super-droplets

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
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

         call sdm_rhot_qtrc2p_t(RHOT,QTRC,DENS,pres_scale,t_scale)

         call sdm_condevp(sdm_aslset,                              &
                          sdm_aslmw,sdm_aslion,sdm_dtevl,           &
                          pres_scale,t_scale,QTRC(:,:,:,I_QV),      &
                          sd_fmnum,sd_numasl,sd_fmx,sd_fmy,sd_fmr,      &
                          sd_fmasl,sd_fmri,sd_fmrj,sd_fmrk)

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
  subroutine sdm_adjsdnum(sdm_nadjvar,ni_sdm,nj_sdm,nk_sdm,    &
                          sd_num,sd_numasl,sd_nc,                 &
                          sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk, &
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
      real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
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
                           sdnum_upr,sd_num,sd_n,sd_x,sd_y,sd_ri,sd_rj,sd_rk,       &
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
                        sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,         &
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
                          sdnum_upr,sd_num,sd_n,sd_x,sd_y,sd_ri,sd_rj,sd_rk,  &
                          sort_id,sort_key,sort_freq,sort_tag,    &
                          sd_rng,sd_rand,                         &
!                          sd_rand,                                &
                          sort_tag0,fsort_id,isd_perm)

!      use m_gadg_algorithm, only: &
!          gadg_count_sort
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort,sdm_getperm
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
      real(RP), intent(inout) :: sd_ri(1:sd_num)    ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num)    ! index[j/real] of super-droplets
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

      call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

      ! Initialize

      gnum = ni_sdm * nj_sdm * knum_sdm

      freq_max = 1

      ! Sorting valid super-droplets

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
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
                       sd_n,sd_x,sd_y,sd_r,sd_asl,sd_vz,sd_ri,sd_rj,sd_rk,    &
                       sort_id,sort_key,sort_freq,sort_tag,       &
                       sd_rng,sd_rand,                            &
!                       sd_rand,                                   &
                       sort_tag0,fsort_id,isd_perm)
!      use m_gadg_algorithm, only: &
!          gadg_count_sort
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_idutil, only: &
         sdm_sort,sdm_getperm
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
      real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
      real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
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

      call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

      ! Initialize
      gnum = ni_sdm * nj_sdm * knum_sdm
      freq_max = 1

      ! Add super-droplets in the grid with small number of super-droplets

      ! Sorting super-droplets

      !### valid droplets ###!

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
                    sort_id,sort_key,sort_freq,sort_tag,'valid')

      !### valid droplets that multiplicity is over 2 ###!

      allocate( sort_freq_valid(1:ni_sdm*nj_sdm*nk_sdm+1) )

      do n=1,ni_sdm*nj_sdm*nk_sdm+1
         sort_freq_valid(n) = sort_freq(n)
      end do

      call sdm_sort(ni_sdm,nj_sdm,nk_sdm,sd_num,sd_n,sd_ri,sd_rj,sd_rk,   &
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
    use m_sdm_coordtrans, only: &
       sdm_x2ri, sdm_y2rj
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QA) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]

    real(RP) :: sum3(KA,IA,JA), sum2(KA,IA,JA)
    real(RP) :: drate             ! temporary
    integer  :: k, i, j
    integer  :: tlist_l            ! total list number for cloud
    integer  :: lcnt               ! counter
    integer  :: kl                 ! index
    integer  :: ku                 ! index
    integer  :: m, n               ! index
    integer  :: ilist_l(1:sdnum_s2c)
    !---------------------------------------------------------------------------
    !-- reset
    Re(:,:,:,:) = 0.0_RP
    sum3(:,:,:) = 0.0_RP
    sum2(:,:,:) = 0.0_RP

    call sdm_x2ri(sdnum_s2c,sdx_s2c,sdri_s2c,sdrk_s2c)
    call sdm_y2rj(sdnum_s2c,sdy_s2c,sdrj_s2c,sdrk_s2c)

    ! Effective radius is defined by r3m/r2m=1.5/lambda.
    ! dcoef is cancelled by the division of sum3/sum2
!    do j  = JS, JE
!    do i  = IS, IE
!    do k  = KS, KE
!       dcoef(k,i,j) = dxiv_sdm(i) * dyiv_sdm(j) / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
!    enddo
!    enddo
!    enddo

    tlist_l = 0

    lcnt = 0

    do n=1,sdnum_s2c

       if( sdrk_s2c(n)<VALID2INVALID ) cycle

       lcnt = lcnt + 1
       ilist_l(lcnt) = n

    end do

    tlist_l = lcnt

    if( tlist_l>0 ) then

       do m=1,tlist_l
          n = ilist_l(m)

          i = floor(sdri_s2c(n))+1
          j = floor(sdrj_s2c(n))+1
          k = floor(sdrk_s2c(n))+1

          sum3(k,i,j) = sum3(k,i,j)            &
               + sdr_s2c(n) * sdr_s2c(n) * sdr_s2c(n)   &
               * real(sdn_s2c(n),kind=RP)
          sum2(k,i,j) = sum2(k,i,j)            &
               + sdr_s2c(n) * sdr_s2c(n)             &
               * real(sdn_s2c(n),kind=RP)
       end do

       !=== correction at the verical boundary ===!

       do j=JS, JE
       do i=IS, IE

          !! at lower boundary

          kl    = floor(sdrkl_s2c(i,j))+1
          drate = real(kl,kind=RP) - sdrkl_s2c(i,j)
          if( drate<0.50_RP ) then
             sum3(kl,i,j) = 0.0_RP           !! <50% in share
             sum2(kl,i,j) = 0.0_RP           !! <50% in share
          else
             sum3(kl,i,j) = sum3(kl,i,j)/drate
             sum2(kl,i,j) = sum2(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sdrku_s2c(i,j))+1
          drate = sdrku_s2c(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sum3(ku,i,j) = 0.0_RP           !! <50% in share
             sum2(ku,i,j) = 0.0_RP           !! <50% in share
          else
             sum3(ku,i,j) = sum3(ku,i,j)/drate
             sum2(ku,i,j) = sum2(ku,i,j)/drate
          end if

       end do
       end do

       !=== convert super-droplets to density of cloud-water. ===!

       do k=KS,KE
       do i=IS,IE
       do j=JS,JE
          if( sum2(k,i,j) == 0.0_RP ) then
           Re(k,i,j,I_mp_QC) = 0.0_RP
          else
           Re(k,i,j,I_mp_QC) = sum3(k,i,j)/sum2(k,i,j)
          endif
       end do
       end do
       end do

    end if

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

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do k  = KS, KE
    do i  = IS, IE
    do j  = JS, JE
       qhydro = 0.0_RP
       do iq = QQS, QQE 
          qhydro = qhydro + QTRC_sdm(k,i,j,iq)
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_sdm_CloudFraction
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_Mixingratio( &
       Qe,    &
       QTRC0  )
    use scale_tracer, only: &
       QAD => QA , &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]

    integer  :: ihydro
    !---------------------------------------------------------------------------

!!    do ihydro = 1, 2
!!       Qe(:,:,:,I_mp_QC) = QTRC_sdm(:,:,:,1)+QTRC_sdm(:,:,:,2)
!!    enddo
    Qe(:,:,:,I_mp_QC) = QTRC_sdm(:,:,:,I_QC)+QTRC_sdm(:,:,:,I_QR)

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
