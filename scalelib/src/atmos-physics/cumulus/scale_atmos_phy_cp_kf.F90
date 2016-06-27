module scale_atmos_phy_cp_kf
  !------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !------------------------------------------------------------------------------
  implicit none
  private
  !------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_kf_setup
  public :: ATMOS_PHY_CP_kf

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP),private,save,allocatable :: lifetime(:,:)!(IA,JA) ! convectime lifetime [s]
  integer ,private,save,allocatable :: I_convflag(:,:)!(IA,JA)
  ! convectiontype 0:deep convection 1:shallow convection 2: noconvection allow
  ! kf time controll
  integer,private,save :: TIME_RES_KF   ! time step for kf
  integer,private,save :: TIME_DSTEP_KF ! time interval
  logical,private,save :: TIME_DOKF     ! exclude kf trigger
  ! tuning parameter
  logical,private,save :: PARAM_ATMOS_PHY_CP_kf_wadapt = .true.
  integer,private,save :: PARAM_ATMOS_PHY_CP_kf_w_time = 16
contains
  !------------------------------------------------------------------------------
  !> Setup
  !------------------------------------------------------------------------------
  subroutine ATMOS_PHY_CP_kf_setup (CP_TYPE) !< [IN]
    use scale_process, only: &
         PRC_MPIstop
    use scale_time , only :&
         dt => TIME_DTSEC_ATMOS_PHY_CP ! now assume equal to TIME_DTSEC
    use scale_atmos_phy_cp_kf_sub,only: &
         kf_lutab, &
         CP_kf_param
    implicit none
    character(len=*), intent(in) :: CP_TYPE !< [IN]
    ! tunning parameters
    ! original parameter set  KF 2004(WRF) and JMA-NHM
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_rate      = 0.03_RP  ! ratio of cloud water and precipitation (Ogura and Cho 1973)
    integer  :: PARAM_ATMOS_PHY_CP_kf_trigger   = 1        ! 1 or 3 trigger type WRF:1 JMA:3
    logical  :: PARAM_ATMOS_PHY_CP_kf_qs        = .true.   ! FLAG_OS:  qs is allowed or not
    logical  :: PARAM_ATMOS_PHY_CP_kf_qi        = .true.   ! FLAG_QI:  qi is allowe or not
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dt        = 5._RP    ! KF convection check time interval [min]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dlcape    = 0.1_RP   ! cape decleace rate
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dlifetime = 1800._RP ! minimum lifetimescale of deep convection
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_slifetime = 2400._RP ! shallow convection liftime 
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_DEPTH_USL = 300._RP  ! depth of updraft source layer  [hPa]!!
    integer  :: PARAM_ATMOS_PHY_CP_kf_cond      = 1        ! condload select 1. Ogura and Cho 1973 2. Kessler(leke JMA ver)
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_thres     = 1.e-3_RP ! kessler type autoconversion rate 
    logical  :: PARAM_ATMOS_PHY_CP_kf_warmrain  = .false.  ! QA<=3 then true
    logical  :: PARAM_ATMOS_PHY_CP_kf_LOG       = .false.  ! KF infomation output to log file(not ERROR messeage)
    ! [INTERNAL WORK]
    integer  :: ierr
    !> [Namelist variable]
    NAMELIST / PARAM_ATMOS_PHY_CP_KF / & !> KF tune parameters
         PARAM_ATMOS_PHY_CP_kf_rate, &
         PARAM_ATMOS_PHY_CP_kf_trigger,  & 
         PARAM_ATMOS_PHY_CP_kf_qs, &
         PARAM_ATMOS_PHY_CP_kf_qi, &
         PARAM_ATMOS_PHY_CP_kf_dt, &
         PARAM_ATMOS_PHY_CP_kf_dlcape,&
         PARAM_ATMOS_PHY_CP_kf_dlifetime, &
         PARAM_ATMOS_PHY_CP_kf_slifetime, &
         PARAM_ATMOS_PHY_CP_kf_DEPTH_USL, &
         PARAM_ATMOS_PHY_CP_kf_cond,&
         PARAM_ATMOS_PHY_CP_kf_thres,&
         PARAM_ATMOS_PHY_CP_kf_warmrain, &
         PARAM_ATMOS_PHY_CP_kf_LOG, &
         !
         PARAM_ATMOS_PHY_CP_kf_wadapt, &
         PARAM_ATMOS_PHY_CP_kf_w_time
    !---------------------------------------------------------------------------
    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CUMULUS] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Kain Fritsch(KF) scheme'

    if ( CP_TYPE /= 'KF' ) then
       write(*,*) 'xxx ATMOS_PHY_CP_TYPE is not KF. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_KF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CP_KF. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_CP_KF)

    call kf_lutab ! set kf look up table
    ! out put variables
    allocate(lifetime(IA,JA))
    lifetime = 0._RP
    allocate(I_convflag(IA,JA))
    I_convflag = 2
    ! set kf convection check step
    TIME_DSTEP_KF = nint(PARAM_ATMOS_PHY_CP_kf_dt*60._RP*1.d3)/nint(dt*1.d3) 
    TIME_DSTEP_KF = max(TIME_DSTEP_KF, 1) ! kf time interval step
    TIME_RES_KF = -1   ! initialize to keep consistent for below step
    TIME_DOKF = .true. ! initialize
    call  CP_kf_param( & ! [IN] ! set paremeters in sub module 
         PARAM_ATMOS_PHY_CP_kf_rate, &
         PARAM_ATMOS_PHY_CP_kf_trigger,  & 
         PARAM_ATMOS_PHY_CP_kf_qs, &
         PARAM_ATMOS_PHY_CP_kf_qi, &
         PARAM_ATMOS_PHY_CP_kf_dt, &
         PARAM_ATMOS_PHY_CP_kf_dlcape,&
         PARAM_ATMOS_PHY_CP_kf_dlifetime, &
         PARAM_ATMOS_PHY_CP_kf_slifetime, &
         PARAM_ATMOS_PHY_CP_kf_DEPTH_USL, &
         PARAM_ATMOS_PHY_CP_kf_cond,&
         PARAM_ATMOS_PHY_CP_kf_thres,&
         PARAM_ATMOS_PHY_CP_kf_warmrain, &
         PARAM_ATMOS_PHY_CP_kf_LOG , &
         TIME_DSTEP_KF )
    ! output parameter lists
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) "*** KF checked interval :",TIME_DSTEP_KF,"step"
    if( IO_L ) write(IO_FID_LOG,*) "*** Ogura and Cho condense material convert rate:",PARAM_ATMOS_PHY_CP_kf_rate
    if( IO_L ) write(IO_FID_LOG,*) "*** trigger fanction : ",PARAM_ATMOS_PHY_CP_kf_trigger
    if( IO_L ) write(IO_FID_LOG,*) "*** exist qs ? : ", PARAM_ATMOS_PHY_CP_kf_qs
    if( IO_L ) write(IO_FID_LOG,*) "*** exist qi ? : ", PARAM_ATMOS_PHY_CP_kf_qi
    if( IO_L ) write(IO_FID_LOG,*) "*** cape decrease rate : ",PARAM_ATMOS_PHY_CP_kf_dlcape
    if( IO_L ) write(IO_FID_LOG,*) "*** deep convection minimum lifetime : ",PARAM_ATMOS_PHY_CP_kf_dlifetime
    if( IO_L ) write(IO_FID_LOG,*) "*** shallow convection  lifetime : ",PARAM_ATMOS_PHY_CP_kf_slifetime
    if( IO_L ) write(IO_FID_LOG,*) "*** Updraft souce layer depth [Pa] : ",PARAM_ATMOS_PHY_CP_kf_DEPTH_USL
    if( IO_L ) write(IO_FID_LOG,*) "*** condload type 1. Ogura and Cho (1973) 2. Kessuler :", PARAM_ATMOS_PHY_CP_kf_cond
    if( IO_L ) write(IO_FID_LOG,*) "*** Kessuler type condload's threshold:", PARAM_ATMOS_PHY_CP_kf_thres
    if( IO_L ) write(IO_FID_LOG,*) "*** warmrain ? :", PARAM_ATMOS_PHY_CP_kf_warmrain
    if( IO_L ) write(IO_FID_LOG,*) "*** Runningmean w use adaptive timestep ? :",PARAM_ATMOS_PHY_CP_kf_wadapt
    if( IO_L ) write(IO_FID_LOG,*) "*** Runningmean w use timestep (if not adaptive) :",PARAM_ATMOS_PHY_CP_kf_w_time
    if( IO_L ) write(IO_FID_LOG,*) "*** kf messege out ? :",PARAM_ATMOS_PHY_CP_kf_LOG
    if( IO_L ) write(IO_FID_LOG,*) "+++ DONE KF SETUP +++"
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine ATMOS_PHY_CP_kf_setup
  !------------------------------------------------------------------------------
  ! > main loop
  !------------------------------------------------------------------------------
  subroutine ATMOS_PHY_CP_kf ( &
       ! [IN]
       DENS, & !> [IN]
       MOMZ, & !> [IN]
       MOMX, & !> [IN]
       MOMY, & !> [IN]
       RHOT, & !> [IN]
       QTRC, & !> [IN]
       ! [INOUT]
       DENS_t_CP,       & ! [INOUT]
       MOMZ_t_CP,       & ! DAMMEY(INOUT)
       MOMX_t_CP,       & ! DAMMEY(INOUT)
       MOMY_t_CP,       & ! DAMMEY(INOUT)
       DT_RHOT,         & ! [INOUT]
       DT_RHOQ,         & ! [INOUT]
       MFLX_cloudbase,  & ! DAMMEY(INOUT)
       SFLX_convrain,   & ! [INOUT]
       cloudtop,        & ! [INOUT]
       cloudbase,       & ! [INOUT]
       cldfrac_dp,      & ! [INOUT]
       cldfrac_sh,      & ! [INOUT]
       nca,             & ! [INOUT]
       w0avg )            ! [INOUT]
    use scale_process, only: &
         PRC_MPIstop
    use scale_history, only: &
         HIST_in
    use scale_time , only :&
         dt => TIME_DTSEC_ATMOS_PHY_CP
    use scale_grid
    use scale_const
    use scale_tracer
    use scale_atmos_phy_cp_kf_sub
    implicit none
    ! [in]
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    ! [out]
    real(RP), intent(inout) :: DENS_t_CP(KA,IA,JA)  !
    real(RP), intent(inout) :: MOMZ_t_CP(KA,IA,JA)  ! dammy not change
    real(RP), intent(inout) :: MOMX_t_CP(KA,IA,JA)  ! dammy not change
    real(RP), intent(inout) :: MOMY_t_CP(KA,IA,JA)  ! dammy not change
    real(RP), intent(inout) :: DT_RHOT(KA,IA,JA)    !
    real(RP), intent(inout) :: DT_RHOQ(KA,IA,JA,QA) !
    !
    real(RP), intent(inout) :: MFLX_cloudbase(IA,JA) ! dammy not change
    real(RP), intent(inout) :: SFLX_convrain(IA,JA)  ! convective rain rate [kg/m2/s]
    real(RP), intent(inout) :: cloudtop(IA,JA)       ! cloud top height [m]
    real(RP), intent(inout) :: cloudbase(IA,JA)      ! cloud base height [m]
    real(RP), intent(inout) :: cldfrac_dp(KA,IA,JA)  ! cloud fraction (deep convection)
    real(RP), intent(inout) :: cldfrac_sh(KA,IA,JA)  ! cloud fraction (shallow convection) 
    real(RP), intent(inout) :: nca(IA,JA)            ! convection active time [s] NOTE: must small advective time scale
    real(RP), intent(inout) :: W0avg(KA,IA,JA)       ! running mean vertical velocity[m/s]
    ! 3d variable
    real(RP) :: cz(KA,IA,JA)          ! centerlevel real height [m]
    real(RP) :: u(KA,IA,JA)           ! x-direction velocity [m/s]
    real(RP) :: v(KA,IA,JA)           ! y-direction velocity [m/s]
    real(RP) :: temp(KA,IA,JA)        ! temperature [K]
    real(RP) :: pres(KA,IA,JA)        ! pressure [Pa]
    real(RP) :: qv(KA,IA,JA)          ! water vaper mixing ratio [kg/kg]
    real(RP) :: q_hyd(KA,IA,JA,MP_QA) ! water mixing ratio [kg/kg]
    real(RP) :: qes(KA,IA,JA)         ! saturate water vaper mixing ratio [kg/kg]
    real(RP) :: rh(KA,IA,JA)          ! relative humidity [%]
    real(RP) :: deltap(KA,IA,JA)      ! pressure interval [hPa]
    real(RP) :: deltaz(KA,IA,JA)      ! height interval (center level) [m]
    real(RP) :: deltax                ! delta x [m]
    real(RP) :: cldfrac_kf(KA,2)      ! 1 shallow , 2 deep
    ! do loop vars
    integer :: ii,jj,kk
    !
    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cumulus Parameterization(KF) '
    TIME_RES_KF = TIME_RES_KF + 1
    if (TIME_RES_KF == TIME_DSTEP_KF) then ! mod(TIME_SETP,TIME_DSTEP_KF) == 0 
       TIME_DOKF = .true.
       TIME_RES_KF = 0
    end if
    if (TIME_DOKF ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** KF Convection Check '
       do jj = JS,JE
          do ii = IS,IE
             nca(ii,jj) = nca(ii,jj) - real(TIME_DSTEP_KF,RP)*dt
!             nca(ii,jj) = nca(ii,jj) - dt
          end do
       end do
       call CP_kf_init ( & !> convert variable
            ! [IN]
            DENS(:,:,:),&
            MOMZ(:,:,:),&
            MOMX(:,:,:),&
            MOMY(:,:,:),&
            RHOT(:,:,:),&
            QTRC(:,:,:,:),&
            ! [OUT]
            cz(:,:,:),&
            u(:,:,:),&
            v(:,:,:),&
            temp(:,:,:),&
            pres(:,:,:),&
            qv(:,:,:),&
            q_hyd(:,:,:,:),&
            qes(:,:,:),&
            rh(:,:,:),&
            deltap(:,:,:),&
            deltaz(:,:,:),&
            deltax, &
            W0avg(:,:,:))
       ! KF main loop
       call PROF_rapstart('CP_kf', 3)
       do jj = JS,JE
          do ii = IS,IE
             if ( nca(ii,jj) < 0.5*dt ) then ! start convection check
                ! convert 3d variable to 1d variable
                ! main loop @ scale_atomos_phy_cp_kf_sub.F90
                call CP_kf_main ( &
                     ! [IN]
                     DENS          (:,ii,jj),   &
                     RHOT          (:,ii,jj),   &
                     QTRC          (:,ii,jj,:), &
                     w0avg         (:,ii,jj),   &
                     u             (:,ii,jj),   &
                     v             (:,ii,jj),   &
                     temp          (:,ii,jj),   &
                     pres          (:,ii,jj),   &
                     qv            (:,ii,jj),   &
                     q_hyd         (:,ii,jj,:), & ! inout 
                     qes           (:,ii,jj),   &
                     rh            (:,ii,jj),   &
                     deltap        (:,ii,jj),   &
                     deltaz        (:,ii,jj),   &
                     CZ            (:,ii,jj),   &
                     deltax,                    &
                     ! [INOUT]
                     nca           (ii,jj)  ,   &
                     ! [OUT]
                     DENS_t_CP     (:,ii,jj),   &
                     DT_RHOT       (:,ii,jj),   &
                     DT_RHOQ       (:,ii,jj,:), &
                     SFLX_convrain (ii,jj),     &
                     cldfrac_kf    (:,:),       &
                     lifetime      (ii,jj),     &
                     cloudtop      (ii,jj),     &
                     cloudbase     (ii,jj),     &
                     I_convflag    (ii,jj))
                cldfrac_sh(KS:KE,ii,jj) = cldfrac_KF(KS:KE,1)
                cldfrac_dp(KS:KE,ii,jj) = cldfrac_KF(KS:KE,2)
             end if
          end do
       end do
       call PROF_rapend('CP_kf', 3)
    else
 !      if( IO_L ) write(IO_FID_LOG,*) '*** Ranning mean vertical velocity '
       call kf_wranmean(w0avg,DENS,MOMZ)
    end if
    TIME_DOKF = .false.
    ! OUTPUT
    call HIST_in( lifetime, 'KF_LIFETIME', ' lifetime of KF scheme', 's')
    call HIST_in( real(I_convflag,RP), 'KF_CONVFLAG', 'CONVECTIONFLAG ', 'NONE')
    return
  end subroutine ATMOS_PHY_CP_kf
  !------------------------------------------------------------------------------
  ! private subroutines
  !------------------------------------------------------------------------------
  subroutine CP_kf_init ( & !> convert variables
       ! [IN]
       DENS   ,&
       MOMZ   ,&
       MOMX   ,&
       MOMY   ,&
       RHOT   ,&
       QTRC   ,&
       ! [OUT]
       z_out  ,&
       u      ,&
       v      ,&
       temp   ,&
       pres   ,&
       qv     ,&
       q_hyd  ,&
       qes    ,&
       rh     ,&
       deltap ,&
       deltaz ,&
       deltax ,&
       w0avg)
    use scale_precision
    use scale_grid_index
    use scale_tracer
    use scale_grid,only: &
         DX => DX, &
         DY => DY
    use scale_const, only: &
         GRAV => CONST_GRAV,&
         R    => CONST_Rdry
    use scale_atmos_thermodyn, only: &
         THERMODYN_temp_pres   => ATMOS_THERMODYN_temp_pres, &
         THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,      &
         THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E, &
         THERMODYN_qd          => ATMOS_THERMODYN_qd 
    use scale_atmos_saturation ,only :&
         ATMOS_SATURATION_psat_liq
    use scale_grid_real, only: &
         CZ   => REAL_CZ,&
         FZ   => REAL_FZ
    implicit none
    ! [IN]
    ! scale prognostic variable
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    ! [OUT]
    ! dynamics variable
    real(RP),intent(out) :: z_out(KA,IA,JA)       ! z level (center level)
    real(RP),intent(out) :: u(KA,IA,JA)           ! x-direction velocity
    real(RP),intent(out) :: v(KA,IA,JA)           ! y-direction velocity
    ! thermodynamics variable
    real(RP),intent(out) :: temp(KA,IA,JA)        ! temperature [K]
    real(RP),intent(out) :: pres(KA,IA,JA)        ! pressure [Pa]
    real(RP),intent(out) :: qv(KA,IA,JA)          ! water vaper mixing ratio [kg/kg]
    real(RP),intent(out) :: q_hyd(KA,IA,JA,MP_QA) ! water mixing ratio [kg/kg] not warervaper
    real(RP),intent(out) :: qes(KA,IA,JA)         ! saturation vaper mixing ratio
    real(RP),intent(out) :: rh(KA,IA,JA)          ! relative humidity
    real(RP),intent(out) :: deltap(KA,IA,JA)      ! pressure difference  calculated by hydrostatic balance)
    real(RP),intent(out) :: deltaz(KA,IA,JA)      ! delta z (center level)
    real(RP),intent(out) :: deltax                ! delta x
    real(RP),intent(out) :: w0avg(KA,IA,JA)       ! running mean vertical velocity
    ! [Internal work]
    integer :: kk,ii,jj,iq
    real(RP) :: es(KA,IA,JA)                 ! saturation vaper pressure
    real(RP) :: qd(KA,IA,JA)                 ! dry mixing ratio
    real(RP) :: dfz(KA,IA,JA)                ! z-length of grid(k+1) to grid(k) [m] ! real
    real(RP) :: dxsq
    z_out(:,:,:) = CZ(:,:,:) ! becouse scale_atmos_phy_cp interface ,not use scale_grid
    ! WRF comment out for W0
    !...TST IS THE NUMBER OF TIME STEPS IN 10 MINUTES...W0AVG IS CLOSE TO A    
    !...RUNNING MEAN VERTICAL VELOCITY...NOTE THAT IF YOU CHANGE TST, IT WIL  
    !...CHANGE THE FREQUENCY OF THE CONVECTIVE INTITIATION CHECK (SEE BELOW) 
    !...NOTE THAT THE ORDERING OF VERTICAL LAYERS MUST BE REVERSED FOR W0AVG
    !...BECAUSE THE ORDERING IS REVERSED IN KFPARA...
    ! calculate w0avg
    ! this is adaptive step from WRF
    call kf_wranmean(w0avg,DENS,MOMZ)

    ! calculate u(x-directin velocity ), v(y-direction velocity)
    do jj = JS  , JE
    do ii = IS, IE
    do kk = KS, KE
       u(kk,ii,jj) = 0.5_RP * ( MOMX(kk,ii,jj) + MOMX(kk,ii-1,jj) ) / DENS(kk,ii,jj)
    end do
    end do
    end do

    do jj = JS  , JE
    do ii = IS, IE
    do kk = KS, KE
       v(kk,ii,jj) = 0.5_RP * ( MOMY(kk,ii,jj) + MOMY(kk,ii,jj-1) ) / DENS(kk,ii,jj)
    end do
    end do
    end do

    call THERMODYN_temp_pres( temp(:,:,:), pres(:,:,:), &               ! [out]
                              DENS(:,:,:), RHOT(:,:,:), QTRC(:,:,:,:) ) ! [in]
    call THERMODYN_qd ( qd,QTRC)
    ! calculate water vaper and relative humidity
    call ATMOS_SATURATION_psat_liq(es(:,:,:),temp(:,:,:))
    !OCL XFILL
    do jj = JS, JE
    do ii = IS, IE
    do kk = KS, KE
       qes(kk,ii,jj) = 0.622_RP * es(kk,ii,jj)/( pres(kk,ii,jj) - (1._RP -0.6222_RP)*es(kk,ii,jj) ) 
       qv(kk,ii,jj) = QTRC(kk,ii,jj,I_QV)/qd(kk,ii,jj)
       qv(kk,ii,jj)  = min(qes(kk,ii,jj), qv(kk,ii,jj)) ! conpare qes and qv
       qv(kk,ii,jj)  = max(0.000001_RP, qv(kk,ii,jj))   ! gess lower limit
       rh(kk,ii,jj) = qv(kk,ii,jj)/qes(kk,ii,jj)
    end do
    end do
    end do
    do iq = 1,MP_QA
    do jj = JS, JE
    do ii = IS, IE
    do kk = KS, KE
       q_hyd(kk,ii,jj,iq) = QTRC(kk,ii,jj,I_MP2ALL(iq))/qd(kk,ii,jj)
    end do
    end do
    end do
    end do 
    ! calculate delta z
    deltaz(:,:,:) = 0._RP ! initialize
    do jj = JS,JE
    do ii = IS,IE
    do kk = KS,KE
       ! deltaz is the interval of between model full levels(scalar point )
       deltaz(kk,ii,jj)  = CZ(kk+1,ii,jj) - CZ(kk,ii,jj)
       dfz(kk,ii,jj) = FZ(kk+1,ii,jj) - FZ(kk,ii,jj)
    end do
    end do
    end do
    ! add
    deltaz(KE,:,:) = 0._RP
    ii = IS; jj =JS
    deltax = (DX + DY)*0.5_RP
    dxsq = deltax*deltax
    ! calculate delta P by hydrostatic balance
    !OCL XFILL
    do jj = JS, JE
    do ii = IS, IE
    do kk = KS, KE
       ! deltap is the pressure interval between half levels(face levels) @ SCALE
       deltap(kk,ii,jj) = DENS(kk,ii,jj)*GRAV*dfz(kk,ii,jj) ! rho*g*dz
    end do
    end do
    end do

    return
  end subroutine CP_kf_init

  subroutine kf_wranmean(w0avg,DENS,MOMZ) ! running mean vertical wind speed
    use scale_precision
    use scale_grid_index
    use scale_time , only :&
         DTmodel =>  TIME_DTSEC ! time interval of model
    use scale_atmos_phy_cp_kf_sub,only: &
         DTcp =>  KF_DT  ! time interval of KF time interval
    implicit none
    ![OUT]
    real(RP), intent(out) :: w0avg(KA,IA,JA)
    ![IN]
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    ! [Internal work]
    ! w0avg factors
    real(RP) :: w0
    real(RP) :: w0avgfactor,w0factor,w0den
    ! do loop index
    integer :: kk,ii,jj

    if (PARAM_ATMOS_PHY_CP_kf_wadapt) then
       w0avgfactor = 2._RP*max(DTmodel,DTcp) - DTmodel 
       w0factor = DTmodel
       w0den =  2._RP*max(DTcp,DTmodel)
    else
       ! note w_time is tuning parameter
       w0avgfactor = real(PARAM_ATMOS_PHY_CP_kf_w_time,RP)
       w0factor = 1._RP
       w0den =  1._RP/w0avgfactor
    end if

    do jj = JS, JE
    do ii = IS, IE
    do kk = KS, KE
       w0 = 0.5_RP * ( MOMZ(kk,ii,jj) + MOMZ(kk-1,ii,jj) ) / DENS(kk,ii,jj)
       ! w0avg = (w0avg[old]*cumulustimestep + w0*modeltimestep)/cumulus timestep
       w0avg(kk,ii,jj) = ( w0avg(kk,ii,jj)*w0avgfactor + w0*w0factor )/w0den ! time mean
       ! w0avg = (w0avg*(2*DTcp-DTmodel) + w0*DTmodel)/(2*DTcp) ! time mean
       ! w0avg is module variable
    end do
    end do
    end do

    return
  end subroutine kf_wranmean
end module scale_atmos_phy_cp_kf
