!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          TWP-ICE forcing
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-xx-xx (A.Noda)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_stdio, only: &
     IO_FID_LOG, &
     IO_L
  use scale_precision
  use scale_prof
  use scale_tracer
  use scale_grid_index
  use scale_stdio, only: &
       IO_get_available_fid
  use scale_grid, only: &
       CX => GRID_CX, &
       CY => GRID_CY, &
       CZ => GRID_CZ
  use scale_time, only: &
       TIME_NOWSTEP,&
       TIME_NOWSEC,&
       TIME_DTSEC
  use mod_cpl_vars, only: &
       SST

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_setup
  public :: USER_step
  
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(DP),allocatable :: momz_ls_t(:)
  real(DP),allocatable :: momz_ls_dz_t(:)
  real(DP),allocatable :: z_in(:)
  real(DP),save,allocatable :: time_atm_in(:)
  real(DP),save,allocatable :: time_sst_in(:)
  real(DP),save,allocatable :: sst_in(:)
  real(RP), private, allocatable :: MOMZ_LS(:,:)
  real(RP), private, allocatable :: MOMZ_LS_DZ(:,:)
  real(RP), private, allocatable :: QV_LS(:,:)
  real(RP), private, allocatable :: U_GEOS(:)
  real(RP), private, allocatable :: V_GEOS(:)
  logical,  private, save        :: MOMZ_LS_FLG(6)
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! for surface flux
  real(RP), private, parameter :: Cm0   = 1.0E-3_RP  ! bulk coef. for U*
  real(RP), private, parameter :: visck = 1.5E-5_RP  ! kinematic viscosity

  ! parameters
  real(RP), private, save :: Z00 = 0.0_RP      ! base
  real(RP), private, save :: Z0R = 0.018_RP    ! rough factor
  real(RP), private, save :: Z0S = 0.11_RP     ! smooth factor
  real(RP), private, save :: Zt0 = 1.4E-5_RP
  real(RP), private, save :: ZtR = 0.0_RP
  real(RP), private, save :: ZtS = 0.4_RP
  real(RP), private, save :: Ze0 = 1.3E-4_RP
  real(RP), private, save :: ZeR = 0.0_RP
  real(RP), private, save :: ZeS = 0.62_RP
  real(RP), private, save :: ThS = 300.0_RP

  ! limiter
  real(RP), private, parameter :: Ustar_min =  1.0E-3_RP ! minimum limit of U*

  real(RP), private, parameter :: Z0_min =    1.0E-5_RP ! minimum roughness length of u,v,w
  real(RP), private, parameter :: Zt_min =    1.0E-5_RP ! T
  real(RP), private, parameter :: Ze_min =    1.0E-5_RP ! q

  real(RP), private, save      :: Cm_min  =   1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, save      :: Ch_min  =   1.0E-5_RP ! T
  real(RP), private, save      :: Ce_min  =   1.0E-5_RP ! q
  real(RP), private, parameter :: Cm_max  =   2.5E-3_RP ! maximum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_max  =   1.0_RP    !                       T
  real(RP), private, parameter :: Ce_max  =   1.0_RP    !                       q

  real(RP), private, save      :: U_minM  =    0.0_RP   ! minimum U_abs for u,v,w
  real(RP), private, save      :: U_minH  =    0.0_RP   !                   T
  real(RP), private, save      :: U_minE  =    0.0_RP   !                   q
  real(RP), private, parameter :: U_maxM  =  100.0_RP   ! maximum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH  =  100.0_RP   !                   T
  real(RP), private, parameter :: U_maxE  =  100.0_RP   !                   q

  integer :: mstep=30
  real(RP), save, allocatable :: time_in(:)
  real(RP), save, allocatable :: lhf_in(:)
  real(RP), save, allocatable :: shf_in(:)
  logical, save :: GIVEN_HEAT_FLUX=.false.

  logical,  private, save :: CNST_RAD=.false. ! add constant radiative cooling

  !---
  real(DP), private, save :: TIME0
  real(RP), private, save :: pi2
  integer,  private, save :: Ktop

  logical,  private, save :: USER_do  = .true.

  integer,  private, save :: USER_LS_FLG = 0 !-- 0->no force, 1->TWPICE
  real(RP), private, save :: corioli

  real(RP), private, save :: CNST_SST=276.2_RP

  character(100), private, save :: inbasedir = './'
  character(100), private, save :: fdata_name_atm = 'large_scale_w_force.txt'
  character(100), private, save :: fdata_name_sst = 'sst_force.txt'
  character(100), private, save :: fdata_name_sf = 'srf_flux_force.txt'
  integer, private, save :: fid_data
  integer, private, save :: fid_data_sf
  logical, private, save :: first_in   =.true.

  integer, private, save :: mstep_atm=15
  integer, private, save :: mstep_sst=30
  integer, private, save :: kend=38

  real(RP), allocatable, private, save :: var(:,:)
  real(RP), allocatable, private, save :: wk(:,:)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine USER_setup
    use scale_stdio, only:  &
       IO_FID_CONF
    use scale_process, only:&
       PRC_MPIstop,         &
       PRC_MPIfinish
    use scale_grid, only:   &
       CZ => GRID_CZ
    implicit none


!    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_PHY_SF

    real(RP) :: USER_SF_U_minM ! minimum U_abs for u,v,w
    real(RP) :: USER_SF_U_minH !                   T
    real(RP) :: USER_SF_U_minE !                   q
    real(RP) :: USER_SF_CM_min ! minimum bulk coef. of u,v,w
    real(RP) :: USER_SF_CH_min !                       T
    real(RP) :: USER_SF_CE_min !                       q
    real(RP) :: USER_SF_Z00
    real(RP) :: USER_SF_Z0R
    real(RP) :: USER_SF_Z0S
    real(RP) :: USER_SF_Zt0
    real(RP) :: USER_SF_ZtR
    real(RP) :: USER_SF_ZtS
    real(RP) :: USER_SF_Ze0
    real(RP) :: USER_SF_ZeR
    real(RP) :: USER_SF_ZeS
    real(RP) :: USER_SF_ThS

    namelist / PARAM_USER / &
       USER_do,             &
       inbasedir,           &
       fdata_name_atm,      &
       fdata_name_sst,      &
       fdata_name_sf,       &
       mstep_atm,           &
       mstep_sst,           &
       CNST_RAD,            &
       CNST_SST,            &
       USER_LS_FLG,         &
       !--- for surface flux
       USER_SF_U_minM,      &
       USER_SF_U_minH,      &
       USER_SF_U_minE,      & 
       USER_SF_CM_min,      & 
       USER_SF_CH_min,      & 
       USER_SF_CE_min,      & 
       USER_SF_Z00,         &
       USER_SF_Z0R,         &
       USER_SF_Z0S,         &
       USER_SF_Zt0,         &
       USER_SF_ZtR,         &
       USER_SF_ZtS,         &
       USER_SF_Ze0,         &
       USER_SF_ZeR,         &
       USER_SF_ZeS,         &
       USER_SF_ThS

    integer :: ierr, t
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

! if(io_l)write(IO_FID_LOG,*)'debug stop'  ! ok
!!call PRC_MPIstop
!call PRC_MPIfinish

    allocate( time_sst_in(mstep_sst) )
    allocate( sst_in(mstep_sst) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    !if(io_l)write(IO_FID_LOG,*)'debug stop1',ierr ! ok
    !call PRC_MPIfinish

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)
    !
    U_minM = USER_SF_U_minM
    U_minH = USER_SF_U_minH
    U_minE = USER_SF_U_minE
    CM_min = USER_SF_CM_min
    CH_min = USER_SF_CH_min
    CE_min = USER_SF_CE_min
    Z00    = USER_SF_Z00
    Z0R    = USER_SF_Z0R
    Z0S    = USER_SF_Z0S
    Zt0    = USER_SF_Zt0
    ZtR    = USER_SF_ZtR
    ZtS    = USER_SF_ZtS
    Ze0    = USER_SF_Ze0
    ZeR    = USER_SF_ZeR
    ZeS    = USER_SF_ZeS
    ThS    = USER_SF_ThS
    !
    allocate( time_in(mstep) )
    allocate( lhf_in(mstep) )
    allocate( shf_in(mstep) )
    !
    fid_data = IO_get_available_fid()
    fdata_name_sst=trim(inbasedir)//'/'//trim(fdata_name_sst)
    open(fid_data, file=trim(fdata_name_sst), status='old',iostat=ierr)
    if(ierr /= 0) then
      write(IO_FID_LOG,*) 'Msg : Sub[mod_user_setup]/Mod[uset_setup]'
      write(IO_FID_LOG,*) 'Cannot open the data file for forcing.'
      write(IO_FID_LOG,*) trim(fdata_name_sst)
      write(IO_FID_LOG,*) 'STOP!!'
      call PRC_MPIstop
    endif
    !
    if(IO_L) write(io_fid_log,*) 'Reading external sst'
    read(fid_data,*)
    do t=1, mstep_sst
      read(fid_data,*) time_sst_in(t), sst_in(t)
      if(IO_L) write(io_fid_log,*) t,time_nowsec,time_sst_in(t),sst_in(t)
    enddo
    close(fid_data)
    !
    if( GIVEN_HEAT_FLUX )then
      fid_data_sf = IO_get_available_fid()
      fdata_name_sf=trim(inbasedir)//'/'//trim(fdata_name_sf)
      open(fid_data_sf, file=trim(fdata_name_sf), status='old',iostat=ierr)
      if(ierr /= 0) then
        write(IO_FID_LOG,*) 'Msg : Sub[SF_GRAYZONE_setup]/Mod[sf_grayzone]'
        write(IO_FID_LOG,*) 'Cannot open the data file for forcing.'
        write(IO_FID_LOG,*) trim(fdata_name_sf)
        write(IO_FID_LOG,*) 'STOP!!'
        call PRC_MPIstop
      endif
      read(fid_data_sf,*)
      do t=1, mstep
        read(fid_data_sf,*) time_in(t), lhf_in(t), shf_in(t)
      enddo
      close(fid_data_sf)
    endif

    do t=1, mstep_sst-1
        if( time_nowsec>=time_sst_in(t) )then
          sst(:,:)=( (time_sst_in(t+1)-time_nowsec)*sst_in(t)+(time_nowsec-time_sst_in(t))*sst_in(t+1) )&
                   /(time_sst_in(t+1)-time_sst_in(t))
          exit
        endif
    enddo

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  !-----------------------------------------------------------------------------
  subroutine USER_step
    use scale_stdio, only: &
     IO_get_available_fid, &
     IO_FID_LOG,  &
     IO_L
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
         DENS,    &
         MOMZ,    &
         MOMX,    &
         MOMY,    &
         RHOT,    &
         QTRC,    &
         DENS_tp, &
         MOMZ_tp, &
         MOMX_tp, &
         MOMY_tp, &
         RHOT_tp, &
         QTRC_tp
    use scale_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use scale_time, only: &
        do_phy_sf => TIME_DOATMOS_PHY_SF, &
        dtsf => TIME_DTSEC_ATMOS_PHY_SF, &
        TIME_NOWSEC
    use scale_const, only: &
        GRAV   => CONST_GRAV,   &
        KARMAN => CONST_KARMAN, &
        Rdry   => CONST_Rdry,   &
        Rvap   => CONST_Rvap,   &
        RovCP  => CONST_RovCP,  &
        CPovCV => CONST_CPovCV, &
        CPdry  => CONST_CPdry,  &
        P00    => CONST_PRE00,  &
        T00    => CONST_TEM00,  &
        CPd    => CONST_CPdry,  &
        LH0    => CONST_LH0,    &
        EPSvap => CONST_EPSvap, &
        PSAT0  => CONST_PSAT0
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_history, only: &
         HIST_in

    implicit none

    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: PRES(KA,IA,JA)
    real(RP) :: TEMP(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA), VELY(KA,IA,JA)
    integer :: k, i, j, iq, ierr, t, kk, iv
    integer :: IIS, IIE, JJS, JJE

!   real(RP) :: z_in(kend)=(/ &
!    2.50,     13.33,     33.33,     60.00,     93.33,&
!  133.33,    180.00,    233.33,    293.33,    360.00,&
!  433.33,    513.33,    600.00,    693.33,    793.33,&
!  900.00,   1013.33,   1133.33,   1260.00,   1393.33,&
! 1533.33,   1680.00,   1833.33,   1993.33,   2160.00,&
! 2333.33,   2513.33,   2700.00,   2893.33,   3093.33,&
! 3300.00,   3513.33,   3733.33,   3960.00,   4193.33,&
! 4433.33,   4680.00,   4933.33/)

!   real(RP) :: time_atm_in(mstep_atm)=(/ &
!     0.0,    3600.0,    7200.0,   10800.0,   14400.0,&
! 18000.0,   21600.0,   25200.0,   28800.0,   32400.0,&
! 36000.0,   39600.0,   43200.0,   46800.0,   50400.0)

    !--- for Surface flux
    real(RP) :: FB  = 9.4_RP  ! Louis factor b (bM)
    real(RP) :: FBS = 4.7_RP  ! Louis factor b' (bM/eM = dE/eE = 9.4/2.0)
    real(RP) :: FDM = 7.4_RP  ! Louis factor d of u (dM)
    real(RP) :: FDH = 5.3_RP  ! Louis factor d of T, q (dH)
    real(RP) :: FR  = 0.74_RP ! turbulent Prandtl number (Businger et al. 1971)

    ! work
    real(RP) :: SFLX_MOMZ(IA,JA)
    real(RP) :: SFLX_MOMX(IA,JA)
    real(RP) :: SFLX_MOMY(IA,JA)
    real(RP) :: SFLX_POTT(IA,JA)
    real(RP) :: SFLX_QV  (IA,JA)
    real(RP) :: SHFLX(IA,JA)
    real(RP) :: LHFLX(IA,JA)

    real(RP) :: THETA(IA,JA)

    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Ustar ! friction velocity [m/s]

    real(RP) :: Z0   ! roughness length [m] (momentum,heat,tracer)
    real(RP) :: Zt
    real(RP) :: Ze

    real(RP) :: Cm   ! bulk coefficient (momentum,heat,tracer)
    real(RP) :: Ch
    real(RP) :: Ce

    real(RP) :: a2
    real(RP) :: Fm, Fh, Psih
    real(RP) :: RiB
    real(RP) :: qdry, Rtot, pres_1d, temp_1d
    real(RP) :: pres_evap ! partial pressure of water vapor at surface [Pa]
    real(RP) :: qv_evap   ! saturation water vapor mixing ratio at surface[kg/kg]
    integer :: iw 

    !---------------------------------------------------------------------------

    if ( .not.USER_do ) then
      return
    else

    if( first_in )then
      first_in = .false.

      allocate( MOMZ_LS(KA,2) )
      allocate( MOMZ_LS_DZ(KA,2) )
      allocate( U_GEOS(KA) )
      allocate( V_GEOS(KA) )
      allocate( QV_LS(KA,2) )

      allocate(wk(mstep_atm,1:kend))
      allocate(var(mstep_atm,1:ka))

      allocate( momz_ls_t(ka) )
      allocate( momz_ls_dz_t(ka) )
      allocate( z_in(kend) )
      allocate( time_atm_in(mstep_atm) )
      !
      ! open 1-dim forcing data
      fid_data = IO_get_available_fid()
      fdata_name_atm=trim(inbasedir)//'/'//trim(fdata_name_atm)
      open(fid_data, file=trim(fdata_name_atm), status='old',iostat=ierr)
      if(ierr /= 0) then
        write(IO_FID_LOG,*) 'Msg : Sub[mod_user_setup]/Mod[user_setup]'
        write(IO_FID_LOG,*) 'Cannot open the data file for forcing.'
        write(IO_FID_LOG,*) trim(fdata_name_atm)
        write(IO_FID_LOG,*) 'STOP!!'
        call PRC_MPIstop
      endif
      !
      do iv=1, 3
        read(fid_data,*)
        if(    iv==1)then
          read(fid_data,'(f9.2,4f11.2)') (z_in(k),k=1,kend)
        elseif(iv==2)then
          read(fid_data,'(f9.2,4f11.2)') (time_atm_in(k),k=1,mstep_atm)
        elseif(iv==3)then
          do t=1, mstep_atm
            read(fid_data,'(f9.4,3f11.4)') (wk(t,k),k=1,kend)
          enddo
        endif
      enddo
      close(fid_data)
      !
      if(IO_L)then
        write(io_fid_log,*) 'w forcing height levels:'
        write(*,'(f9.2,4f11.2)') (z_in(k),k=1,kend)
        write(io_fid_log,*) 'w forcing time levels:'
        write(*,'(f9.2,4f11.2)') (time_atm_in(k),k=1,mstep_atm)
        write(io_fid_log,*) 'w forcing:'
        do t=1, mstep_atm
          write(io_fid_log,'(f9.4,3f11.4)')(wk(t,k),k=1,kend)
        enddo
      endif
      !
      do k=ks, ke
        do kk=2, kend
          if( z_in(kk)>cz(k) )then
            var(:,k)=( (z_in(kk)-cz(k))*wk(:,kk-1)+(cz(k)-z_in(kk-1))*wk(:,kk) )&
                    /(z_in(kk)-z_in(kk-1))
            exit
          endif
        enddo
      enddo

    endif

    if( USER_LS_FLG == 0 ) then  ! no large scale sinking

       MOMZ_LS(:,:) = 0.0_RP
       MOMZ_LS_DZ(:,:) = 0.0_RP
       MOMZ_LS_FLG( : ) = .false.
       QV_LS(:,:) = 0.0_RP
       V_GEOS(:) = 0.0_RP
       U_GEOS(:) = 0.0_RP
       corioli = 0.0_RP

    elseif( USER_LS_FLG == 1 ) then 

      do t=1, mstep_atm-1
        if( time_nowsec>time_atm_in(t) )then
          do k=1, ka
            momz_ls_t(k)=( (time_atm_in(t+1)-time_nowsec)*var(t,k)+(time_nowsec-time_atm_in(t))*var(t+1,k) )&
                   /(time_atm_in(t+1)-time_atm_in(t))
          enddo
          do k=2, ka-1
            momz_ls_dz_t(k)=(momz_ls_t(k+1)-momz_ls_t(k-1))/(cz(k+1)-cz(k-1))
          enddo
          momz_ls_dz_t(1) =momz_ls_t(2)/cz(2)
          momz_ls_dz_t(ka)=momz_ls_t(ka-1)
          exit
        endif
      enddo
      if( time_nowsec>time_atm_in(mstep_atm) )then
        write(*,*) 'Integration time exceeds the maximum forcing data length',time_nowsec,time_atm_in(mstep_atm)
        call PRC_MPIstop
      endif

      do t=1, mstep_sst-1
! write(*,*)'chksstuser0',t,time_nowsec,time_sst_in(t),sst_in(t)
        if( time_nowsec>=time_sst_in(t) )then
          sst(:,:)=( (time_sst_in(t+1)-time_nowsec)*sst_in(t)+(time_nowsec-time_sst_in(t))*sst_in(t+1) )&
                   /(time_sst_in(t+1)-time_sst_in(t))
          exit
        endif
      enddo
! write(*,*)'chksstuser1',maxval(sst(:,:)), minval(sst(:,:))
 
      if( time_nowsec>time_sst_in(mstep_sst) )then
        write(*,*) 'Integration time exceeds the maximum forcing data length',time_nowsec,time_sst_in(mstep_sst)
        call PRC_MPIstop
      endif

       MOMZ_LS(:,1)=MOMZ_LS_T(:)
       MOMZ_LS_DZ(:,1)=MOMZ_LS_DZ_T(:)
       MOMZ_LS_FLG(:) = .true.
       U_GEOS(:) = 0.0
       V_GEOS(:) = -15.0-0.0024*cz(:)
       corioli = 1e-5 ! tentative. need to ask stephan 
       MOMZ_LS(:,2)=0.0_RP
       MOMZ_LS_DZ(:,2)=0.0_RP
       do k=KS, KE
         MOMZ_LS(k,2)=(MOMZ_LS(k-1,1)+MOMZ_LS(k,1))*0.5
         MOMZ_LS_DZ(k,2)=(MOMZ_LS_DZ(k-1,1)+MOMZ_LS_DZ(k,1))*0.5
!        Qv_LS(k,2)=(QV_LS(k-1,1)+QV_LS(k,1))*0.5
       enddo
!      QV_LS(:,1) = QV_LS_T(:)
!      QV_LS(:,2)=0.0_RP
       QV_LS(:,:)=0.0_RP

!do k=1,ka
!write(*,*)'chk1',k,momz_ls(k,1),momz_ls_dz(k,1),v_geos(k),momz_ls(k,2),momz_ls_dz(k,2)
!enddo

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       if ( MOMZ_LS_FLG(I_MOMZ) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             WORK(k,i,j) = MOMZ(k,i,j) * 2.0_RP / ( DENS(k+1,i,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-2
             MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) &
                  - MOMZ_LS(k,2) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RCDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMZ_tp(KE-1,i,j) = MOMZ_tp(KE-1,i,j) &
                  - MOMZ_LS(KE-1,2) * (           - WORK(KE-1,i,j) ) * RCDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_MOMX) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE+1
          do k = KS, KE
             VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMX_tp(k,i,j) = MOMX_tp(k,i,j) &
                  + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                  * ( - CORIOLI * V_GEOS(k) &
                      + CORIOLI * 0.25_RP &
                      * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) &
                    ) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMX_tp(KE,i,j) = MOMX_tp(KE,i,j) &
                  + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                  *  ( - CORIOLI * V_GEOS(KE) &
                       + CORIOLI * 0.25_RP &
                       * ( VELY(KE,i,j)+VELY(KE,i+1,j)+VELY(KE,i,j-1)+VELY(KE,i+1,j-1) ) &
                     ) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_MOMY) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = MOMY(k,i,j) * 2.0_RP / ( DENS(k,i,j+1) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE+1
          do i = IIS-1, IIE
          do k = KS, KE
             VELX(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMY_tp(k,i,j) = MOMY_tp(k,i,j) &
                  + 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )  &
                  * ( + CORIOLI * U_GEOS(k) &
                      - CORIOLI * 0.25_RP &
                      * ( VELX(k,i,j)+VELX(k,i,j+1)+VELX(k,i-1,j)+VELX(k,i-1,j+1) ) &
                    ) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE+1
          do i = IIS-1, IIE
             MOMY_tp(KE,i,j) = MOMY_tp(KE,i,j) &
                  + 0.5_RP * ( DENS(KE,i,j+1)+DENS(KE,i,j) ) &
                  * ( + CORIOLI * U_GEOS(KE) &
                      - CORIOLI * 0.25_RP  &
                      * ( VELX(KE,i,j)+VELX(KE,i,j+1)+VELX(KE,i-1,j)+VELX(KE,i-1,j+1) ) &
                    ) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_RHOT) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             RHOT_tp(k,i,j) = RHOT_tp(k,i,j) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             RHOT_tp(KE,i,j) = RHOT_tp(KE,i,j) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo

          if( CNST_RAD )then
            !--- add constant cooling (-2K/dy)
            call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                                      PRES(:,:,:),  & ! [OUT]
                                      DENS(:,:,:),  & ! [IN]
                                      RHOT(:,:,:),  & ! [IN]
                                      QTRC(:,:,:,:) ) ! [IN]
            !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
            do j = JJS, JJE
            do i = IIS, IIE
            do k = KS, KE-1
               RHOT_tp(k,i,j) = RHOT_tp(k,i,j)-2.3e-5*(1.0e5/PRES(k,i,j))**0.28586
            enddo
            enddo
            enddo
          endif

       end if

       if ( MOMZ_LS_FLG(I_QTRC) ) then

          do iq = 1, QA
             !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE-1
                QTRC_tp(k,i,j,iq) = QTRC_tp(k,i,j,iq) &
                     - MOMZ_LS(k,1) * ( QTRC(k+1,i,j,iq) - QTRC(k,i,j,iq) ) * RFDZ(k)
             enddo
             enddo
             enddo
             !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
                QTRC_tp(KE,i,j,iq) = QTRC_tp(KE,i,j,iq) &
                     - MOMZ_LS(KE,1) * ( QTRC(KE,i,j,iq) - QTRC(KE-1,i,j,iq) ) * RFDZ(KE-1)
             enddo
             enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             QTRC_tp(k,i,j,I_QV) = QTRC_tp(k,i,j,I_QV) + QV_LS(k,1)
          enddo
          enddo
          enddo

       end if

      enddo
      enddo

    else
       write(IO_FID_LOG,*)'Not supported user_ls_flg'
       call PRC_MPIstop
    endif
    endif

    if( do_phy_sf ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface paramaterization of GRAYZONE'
!       do t=1, mstep_sst-1
!write(*,*)'chksstuser0',t,time_nowsec,time_sst_in(t),sst_in(t)
!        if( time_nowsec>=time_sst_in(t) )then
!          sst_loc(:,:)=(
!          (time_sst_in(t+1)-time_nowsec)*sst_in(t)+(time_nowsec-time_sst_in(t))*sst_in(t+1)
!          )&
!                   /(time_sst_in(t+1)-time_sst_in(t))
!          exit
!        endif
!      enddo
!write(*,*)'chksstuser1',maxval(sst_loc(:,:)), minval(sst_loc(:,:))
!
!      if( time_nowsec>time_sst_in(mstep_sst) )then
!        write(*,*) 'Integration time exceeds the maximum forcing data
!        length',time_nowsec,time_sst_in(mstep_sst)
!        call PRC_MPIstop
!      endif

       ! rho*theta -> potential temperature at cell centor
       do j = JS, JE
       do i = IS, IE
          THETA(i,j) = RHOT(KS,i,j) / DENS(KS,i,j)
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE

          ! at cell center

          !--- absolute velocity
          Uabs = sqrt( &
                 ( MOMZ(KS,i,j)                  )**2 &
               + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
               + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 &
               ) / DENS(KS,i,j) * 0.5_RP

          !--- friction velocity at u, v, and w points
!         Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

          !--- roughness lengths at u, v, and w points
!         Z0 = max( Z00 + Z0R/GRAV * Ustar*Ustar + Z0S*visck / Ustar, Z0_min )
!         Zt = max( Zt0 + ZtR/GRAV * Ustar*Ustar + ZtS*visck / Ustar, Zt_min )
!         Ze = max( Ze0 + ZeR/GRAV * Ustar*Ustar + ZeS*visck / Ustar, Ze_min )
          Z0 = 6.6E-4_RP
          Zt = 3.7E-6_RP
          Ze = 3.7E-6_RP

          call get_RiB( &
               RiB, Fm, Fh, Psih,                  & ! (out)
               THETA(i,j), SST(i,j), Uabs**2,      & ! (in)
!              THETA(i,j), SST_loc(i,j), Uabs**2,  & ! (in)
               CZ(KS), Z0, Zt,                     & ! (in)
               KARMAN, FB, FBS, FDM, FDH,          & ! (in)
               ThS, GRAV                           ) ! (in)

          !--- surface exchange coefficients
          a2 = ( KARMAN / log( CZ(KS)/Z0 ) )**2
          Cm = a2 * Fm
          Ch = a2 * Fh / ( FR * ( log( Z0/Zt ) / Psih + 1.0_RP ) )
          Ce = a2 * Fh / ( FR * ( log( Z0/Ze ) / Psih + 1.0_RP ) )

          ! Gas constant
          qdry = 1.0_RP
          do iw = QQS, QQE
             qdry = qdry - QTRC(KS,i,j,iw)
          enddo
          Rtot = Rdry*qdry + Rvap*QTRC(KS,i,j,I_QV)

          !--- saturation at surface
          pres_1d   = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV
          temp_1d   = ( RHOT(KS,i,j) / DENS(KS,i,j) ) * ( P00 / pres_1d )**RovCP
          pres_evap = PSAT0 * exp( LH0/Rvap * ( 1.0_RP/T00 - 1.0_RP/SST(i,j) ) )
!          pres_evap = PSAT0 * exp( LH0/Rvap * ( 1.0_RP/T00 - 1.0_RP/SST_loc(i,j) )
!          )
!          qv_evap   = EPSvap * pres_evap / ( pres_1d - pres_evap )
          qv_evap   = EPSvap * pres_evap / P00

          ! flux
          SFLX_MOMZ(i,j) = - min(max(Cm,Cm_min),Cm_max) * min(max(Uabs,U_minM),U_maxM) &
                         * MOMZ(KS,i,j) * 0.5_RP
          SFLX_POTT(i,j) =   min(max(Ch,Ch_min),Ch_max) * min(max(Uabs,U_minH),U_maxH) &
                         * ( SST(i,j) * DENS(KS,i,j) - RHOT(KS,i,j) )
!           * ( SST_loc(i,j)*DENS(KS,i,j) - RHOT(KS,i,j) )
          SFLX_QV  (i,j) =   min(max(Ce,Ce_min),Ce_max) * min(max(Uabs,U_minE),U_maxE) &
                         * DENS(KS,i,j) * ( qv_evap - QTRC(KS,i,j,I_QV) )

          ! at (u, y, layer)
          Uabs = sqrt( &
               ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i+1,j) ) )**2 &
               + ( 2.0_RP *   MOMX(KS,i,j) )**2 &
               + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
               ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
!         Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

!         Z0 = max( Z00 + Z0R/GRAV * Ustar*Ustar + Z0S*visck / Ustar, Z0_min )
!         Zt = max( Zt0 + ZtR/GRAV * Ustar*Ustar + ZtS*visck / Ustar, Zt_min )
          Z0 = 6.6E-4_RP
          Zt = 3.7E-6_RP

          call get_RiB( &
               RiB, Fm, Fh, Psih,                    & ! (out)
               ( THETA(i,j)+THETA(i+1,j) ) * 0.5_RP, & ! (in)
               ( SST(i,j)+SST(i+1,j) ) * 0.5_RP,     & ! (in)
!              ( SST_loc(i,j)+SST(i+1,j) ) * 0.5_RP, & ! (in)
               Uabs**2,                              & ! (in)
               CZ(KS), Z0, Zt,                       & ! (in)
               KARMAN, FB, FBS, FDM, FDH,            & ! (in)
               ThS, GRAV                             ) ! (in)

          a2 = ( KARMAN / log( CZ(KS)/Z0 ) )**2
          Cm = a2 * Fm

          SFLX_MOMX(i,j) = - min(max(Cm,Cm_min),Cm_min) * min(max(Uabs,U_minM),U_maxM) &
                         * MOMX(KS,i,j)


          ! at (x, v, layer)
          Uabs = sqrt( &
                 ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i,j+1) ) )**2 &
               + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
               + ( 2.0_RP *   MOMY(KS,i,j) )**2 &
               ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
!         Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

!         Z0 = max( Z00 + Z0R/GRAV * Ustar*Ustar + Z0S*visck / Ustar, Z0_min )
!         Zt = max( Zt0 + ZtR/GRAV * Ustar*Ustar + ZtS*visck / Ustar, Zt_min )
          Z0 = 6.6E-4_RP
          Zt = 3.7E-6_RP

          call get_RiB( &
               RiB, Fm, Fh, Psih,                       & ! (out)
               ( THETA(i,j)+THETA(i,j+1) ) * 0.5_RP,    & ! (in)
               ( SST(i,j)+SST(i,j+1) ) * 0.5_RP,        & ! (in)
!              ( SST_loc(i,j)+SST_loc(i,j+1) ) * 0.5_RP,& ! (in)
               Uabs**2,                                 & ! (in)
               CZ(KS), Z0, Zt,                          & ! (in)
               KARMAN, FB, FBS, FDM, FDH,               & ! (in)
               ThS, GRAV                                ) ! (in)

          a2 = ( KARMAN / log( CZ(KS)/Z0 ) )**2
          Cm = a2 * Fm

          SFLX_MOMY(i,j) = - min(max(Cm,Cm_min),Cm_min) * min(max(Uabs,U_minM),U_maxM) &
                         * MOMY(KS,i,j)

      enddo
      enddo

      if( GIVEN_HEAT_FLUX )then
        do t=1, mstep-1
          if( time_nowsec>time_in(t) )then
            SFLX_POTT(:,:)=( (time_in(t+1)-time_nowsec)*shf_in(t)+(time_nowsec-time_in(t))*shf_in(t+1) )&
                          /( time_in(t+1)-time_in(t) ) / CPd
            SFLX_QV  (:,:)=( (time_in(t+1)-time_nowsec)*lhf_in(t)+(time_nowsec-time_in(t))*lhf_in(t+1) )&
                          /( time_in(t+1)-time_in(t ))/LH0
          endif
        enddo
        if( time_nowsec>time_in(mstep) )then
          write(*,*) 'Integration time exceeds the maximum forcing data length',time_nowsec,time_in(mstep)
          call PRC_MPIstop
        endif
      endif

!write(*,*)'chksst',dtsf,maxval(sst(1,:,:)), minval(sst(1,:,:))
!     call HIST_in( SST_loc(:,:), 'SST','sst','K',dtsf)
      call HIST_in( SST(:,:), 'SST2','sst','K',dtsf)
      call HIST_in( SFLX_POTT(:,:)*CPd, 'SHF','shf','W/m2',dtsf)
      call HIST_in( SFLX_QV(:,:)*LH0, 'LHF','lhf','W/m2',dtsf)
!write(*,*)'chksst2'

      !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
      do j = JS, JE
      do i = IS, IE

       SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
       LHFLX(i,j) = SFLX_QV  (i,j) * LH0

       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SFLX_POTT(i,j) &
            + SFLX_QV(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
             + SFLX_QV(i,j) * RCDZ(KS)
       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
             + SFLX_MOMZ(i,j) * RFDZ(KS)
       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
             + SFLX_MOMX(i,j) * RCDZ(KS)
       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
             + SFLX_MOMY(i,j) * RCDZ(KS)
       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
             + SFLX_QV(i,j) * RCDZ(KS)
      enddo
      enddo

      call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2', dtsf )
      call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2', dtsf )

    endif

    return
  end subroutine USER_step
  !---------------------------------------------------------------------------------
  subroutine get_RiB( &
       RiB, Fm, Fh, Psih,    &
       theta, theta_sfc, u2, &
       Z, Z0, Zt,            &
       K, FB, FBS, FDM, FDH, &
       ThS, G                )
    use scale_const, only: &
       EPS    => CONST_EPS
    real(RP), intent(out) :: RiB
    real(RP), intent(out) :: Fm
    real(RP), intent(out) :: Fh
    real(RP), intent(out) :: Psih
    real(RP), intent(in)  :: theta
    real(RP), intent(in)  :: theta_sfc
    real(RP), intent(in)  :: u2
    real(RP), intent(in)  :: Z
    real(RP), intent(in)  :: Z0
    real(RP), intent(in)  :: Zt
    real(RP), intent(in)  :: K
    real(RP), intent(in)  :: FB
    real(RP), intent(in)  :: FBS
    real(RP), intent(in)  :: FDM
    real(RP), intent(in)  :: FDH
    real(RP), intent(in)  :: ThS
    real(RP), intent(in)  :: G

    real(RP) :: tmp

    ! the first guess of RiB0 (= RiBt)
    RiB = G/ThS * z * (  theta -  theta_sfc ) / ( u2+EPS )

    ! Fm, Fh, Psi_h/R
    if ( RiB >= 0 ) then
       Fm = 1.0_RP / ( 1.0_RP + FBS * Rib )**2
       Fh = Fm
    else
       tmp = ( K / log( Z/Z0 ) )**2 * FB * sqrt( Z/Z0 * abs(RiB) )
       Fm = 1.0_RP - FB * RiB / ( 1.0_RP + FDM * tmp )
       Fh = 1.0_RP - FB * RiB / ( 1.0_RP + FDH * tmp )
    endif
    Psih = log( Z/Z0 ) * sqrt( Fm ) / Fh

    ! the final estimate of RiB0
    tmp = log( Z0/Zt )
    RiB = RiB - RiB * tmp / ( tmp + Psih )

    ! Fm, Fh, Psih/R
    if ( RiB >= 0.0_RP ) then
       Fm = 1.0_RP / ( 1.0_RP + FBS * Rib )**2
       Fh = Fm
    else
       tmp = ( K / log( Z/Z0 ) )**2 * FB * sqrt( Z/Z0 * abs(RiB) )
       Fm = 1.0_RP - FB * RiB / ( 1.0_RP + FDM * tmp )
       Fh = 1.0_RP - FB * RiB / ( 1.0_RP + FDH * tmp )
    endif
    Psih = log( Z/Z0 ) * sqrt( Fm ) / Fh

    return

  end subroutine get_RiB
  !---------------------------------------------------------------------------------
end module mod_user
