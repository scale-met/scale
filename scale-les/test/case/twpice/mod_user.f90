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
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  use mod_precision
  use mod_prof
  use mod_tracer
  use mod_grid_index

! use dc_types, only: &
!    DP
    use mod_grid, only: &
       CX => GRID_CX, &
       CY => GRID_CY, &
       CZ => GRID_CZ
    use mod_time, only: &
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
! include "inc_precision.h"
! include "inc_index.h"
! include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
! real(DP), public :: momz_ls_t(ka)
! real(DP), public :: momz_ls_dz_t(ka)
! real(DP), public :: u_geos_t(ka)
! real(DP), public :: v_geos_t(ka)
! real(DP), public :: qv_ls_t(ka)
! real(DP), public :: pott_ls_t(ka)
  real(DP),allocatable :: momz_ls_t(:)
  real(DP),allocatable :: momz_ls_dz_t(:)
  real(DP),allocatable :: u_geos_t(:)
  real(DP),allocatable :: v_geos_t(:)
  real(DP),allocatable :: qv_ls_t(:)
  real(DP),allocatable :: pott_ls_t(:)
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
  real(DP), private, save :: TIME0
  real(RP), private, save :: pi2
  integer,  private, save :: Ktop

  logical,  private, save :: USER_do  = .true.
  real(DP), private, save :: FORCE_DURATION = 1200.D0
  real(RP), private, save :: SHIFT_X = 12.0E0_RP
  real(RP), private, save :: SHIFT_Y = -2.0E0_RP
  real(RP), private, save :: DT_MAX  = -6.7e-3_RP
  real(RP), private, save :: DQ_MAX  = -1.675e-6_RP
  real(RP), private, save :: POOL_TOP  = 2.5e3_RP
  real(RP), private, save :: POOL_CX   = 100.e3_RP
  real(RP), private, save :: POOL_CY0  = 100.e3_RP
  real(RP), private, save :: POOL_RX   = 7.e3_RP
  real(RP), private, save :: POOL_RY   = 6.e3_RP
  real(RP), private, save :: POOL_DIST = 15.e3_RP
  integer,  private, save :: POOL_NUM  = 4

  integer,  private, save :: USER_LS_FLG = 0 !-- 0->no force, 1->TWPICE
  real(RP), private, save :: corioli

  real(RP), private, save :: CNST_SST=302.15_RP

  character(100), private, save :: inbasedir = './'
  character(100), private, save :: fdata_name = 'forcing4run.txt'
  integer, private, save :: mstep = 1  ! max time step in fdata_name
  integer, private, save :: intv   = -999
! integer, private, save :: start_hr=0
  real, private, save :: start_hr=0
  logical, public, save :: CNST_RAD=.false. ! add constant radiative cooling
  integer, private, save :: start_step=-999
  integer, private, save :: fid_data
  logical, private, save :: first_in   =.true.
  logical, private, save :: rd1st=.true.

  integer, private, save :: varmax = 5
  integer, private, save :: nU      =1
  integer, private, save :: nV      =2
  integer, private, save :: nPT_tend=3  ! tendency of pot.temp [K/s]
  integer, private, save :: nQV_tend=4  ! tendency of qv [kg/kg/s]
  integer, private, save :: nW_ls   =5  ! large-scale vertical velocity [m/s]

  integer, private, save :: mtnum
  integer, private, save :: stepv1
  integer, private, save :: stepv2

  real(RP), allocatable, private, save :: var1(:,:)
  real(RP), allocatable, private, save :: var2(:,:)
  real(RP), allocatable, private, save :: wk(:,:)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine USER_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop,&
       PRC_MPIfinish
    use mod_grid, only: &
       CZ => GRID_CZ
    implicit none

    namelist / PARAM_USER / &
       USER_do, &
         inbasedir,    &
         fdata_name,   &
         intv,         &
         mstep,        &
         start_hr,     &
         start_step,   &
         CNST_RAD,&
         CNST_SST,&
       USER_LS_FLG,&
       FORCE_DURATION, &
       DT_MAX, &
       DQ_MAX, &
       SHIFT_X, &
       SHIFT_Y, &
       POOL_CX, &
       POOL_CY0, &
       POOL_TOP, &
       POOL_RX, &
       POOL_RY, &
       POOL_DIST, &
       POOL_NUM

!   integer :: k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

if(io_l)write(IO_FID_LOG,*)'debug stop'  ! ok
!!call PRC_MPIstop
!call PRC_MPIfinish

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

if(io_l)write(IO_FID_LOG,*)'debug stop1',ierr ! ok
!call PRC_MPIfinish
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)
    !
    !--- equivalent time step for the current dt
    if( start_step==-999 ) start_step=int(start_hr/TIME_DTSEC)+1
    !--- assume 3hourly external forcing data
    if( intv==-999 ) intv=int(3.0*3600/TIME_DTSEC)

    SST(:,:)=CNST_SST

!if(io_l)write(IO_FID_LOG,*)'debug stop2',ierr ! ok
!call PRC_MPIfinish
!   if ( USER_do ) then
!      if( IO_L ) write(IO_FID_LOG,*) '*** Enable cold pool forcing'
!      TIME0 = NOWSEC

!      pi2 = atan(1.0) * 2.0_RP

!      Ktop = KE
!      do k = KS, KE
!         if ( CZ(k) > POOL_TOP ) then
!            Ktop = k-1
!            exit
!         end if
!      enddo
!   end if

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  !-----------------------------------------------------------------------------
  subroutine USER_step
    use mod_stdio, only: &
     IO_get_available_fid, &
     IO_FID_LOG,  &
     IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       DENS, &
       RHOT, &
       QTRC
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
         DENS,    &
         MOMZ,    &
         MOMX,    &
         MOMY,    &
         RHOT,    &
         QTRC,    &
         MOMZ_tp, &
         MOMX_tp, &
         MOMY_tp, &
         RHOT_tp, &
         QTRC_tp
    use mod_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use mod_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres

    implicit none

    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: PRES(KA,IA,JA)
    real(RP) :: TEMP(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA), VELY(KA,IA,JA)
    integer :: k, i, j, iq, n, ierr, mt, kk
    integer :: IIS, IIE, JJS, JJE

!   real(RP) :: dt, dq
!   real(RP) :: time
!   real(RP) :: fact, dist
!   real(RP) :: POOL_CY
    !---------------------------------------------------------------------------

    SST(:,:)=CNST_SST

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

      allocate(var1(1:ka,varmax))
      allocate(var2(1:ka,varmax))
      allocate(wk(1:ka,varmax))

      allocate( momz_ls_t(ka) )
      allocate( momz_ls_dz_t(ka) )
      allocate( u_geos_t(ka) )
      allocate( v_geos_t(ka) )
      allocate( qv_ls_t(ka) )
      allocate( pott_ls_t(ka) )

      var1  (:,:) = 0.0d0
      var2  (:,:) = 0.0d0
      !
      ! open 1-dim forcing data
      fid_data = IO_get_available_fid()
      fdata_name=trim(inbasedir)//'/'//trim(fdata_name)
      open(fid_data, file=trim(fdata_name), status='old',iostat=ierr)
      if(ierr /= 0) then
        write(IO_FID_LOG,*) 'Msg : Sub[af_twpice_init]/Mod[af_twpice]'
        write(IO_FID_LOG,*) 'Cannot open the data file for forcing.'
        write(IO_FID_LOG,*) trim(fdata_name)
        write(IO_FID_LOG,*) 'STOP!!'
        call PRC_MPIstop
      endif
      !
      !--- Check if whole input data available
      do mt=1, mstep
        read(fid_data,*) kk
        do k=1, ka
!write(*,*)'chkread',k,kk,mt,mstep,ka
          read(fid_data,*,iostat=ierr) (wk(k,n),n=1,varmax)
!if(io_l)write(io_fid_log,*) 'chkread2',k,(wk(k,n),n=1,varmax)
          if( ierr/=0 )then
            write(*,*) 'Not enough data! ',mt,k,mstep,ka,trim(fdata_name)
            call PRC_MPIstop
          endif
        enddo
      end do
      rewind(fid_data)
      !
      mtnum=int(start_hr*3600.0/TIME_DTSEC/intv) + 1
!write(*,*)'chktime',start_hr,time_dtsec,intv,mtnum
      !
!     if(ADM_prc_me == ADM_prc_run_master) then
        if(io_l) write(IO_FID_LOG,*) 'Forcing starts! ',mtnum
!     end if
      !
      ! skip by the adequate time step
      do mt=1, mtnum-1
        read(fid_data,*) kk
        do k=1, ka
          read(fid_data,*) (wk(k,n),n=1,varmax)
        enddo
      enddo
    endif
    !
    call update_var

    endif
!write(*,*)'chkuserstep'

    if( USER_LS_FLG == 0 ) then  ! no large scale sinking

       MOMZ_LS(:,:) = 0.0_RP
       MOMZ_LS_DZ(:,:) = 0.0_RP
       MOMZ_LS_FLG( : ) = .false.
       QV_LS(:,:) = 0.0_RP
       V_GEOS(:) = 0.0_RP
       U_GEOS(:) = 0.0_RP
       corioli = 0.0_RP

    elseif( USER_LS_FLG == 1 ) then ! DYCOMS

       MOMZ_LS(:,1)=MOMZ_LS_T(:)
       MOMZ_LS_DZ(:,1)=MOMZ_LS_DZ_T(:)
       MOMZ_LS_FLG(:) = .true.
       QV_LS(:,1) = QV_LS_T(:)
       U_GEOS(:) = U_GEOS_T(:)
       V_GEOS(:) = V_GEOS_T(:)
       corioli = 7.292115E-5_RP
       MOMZ_LS(:,2)=0.0_RP
       MOMZ_LS_DZ(:,2)=0.0_RP
       QV_LS(:,2)=0.0_RP
       do k=KS, KE
         MOMZ_LS(k,2)=(MOMZ_LS(k-1,1)+MOMZ_LS(k,1))*0.5
         MOMZ_LS_DZ(k,2)=(MOMZ_LS_DZ(k-1,1)+MOMZ_LS_DZ(k,1))*0.5
         Qv_LS(k,2)=(QV_LS(k-1,1)+QV_LS(k,1))*0.5
       enddo

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

    return
  end subroutine USER_step

  !---------------------------------------------------------------------------------
  subroutine update_var
    !
    use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
    !
    implicit none
    !
    integer :: cstep, k
    real(8) :: f
    !
    ! current time step counted from the initial run
    cstep = int(start_hr*3600/TIME_DTSEC) + TIME_NOWSTEP + 1
    !
!write(*,*)'Update af1',cstep,start_step,intv,rd1st,time_nowstep
    !
    if(rd1st) then
       !
       stepv1 = (mtnum-1)*intv + 1
       stepv2 = stepv1 + intv
       !
!write(*,*)'Update af2',mtnum
       call read_var(var1(:,:))
       mtnum = mtnum + 1 
       call read_var(var2(:,:))
       mtnum = mtnum + 1 
       !
       rd1st = .false.
    else
       if(mod((cstep-start_step), intv) == 0) then
!write(*,*)'Update af3',cstep,start_step,intv,mtnum
          var1(:,:)    = var2(:,:)
          call read_var(var2(:,:))
          mtnum = mtnum + 1 
          !
          stepv1 = stepv2
          stepv2 = stepv2 + intv
       end if
    end if
    !
    f = dble(cstep - stepv1)/dble(intv)
!   var_bs(:,:) = var1(:,:)*(1.0d0 - f) + var2(:,:)*f
    u_geos_t (:) = var1(:,nU)      *(1.0d0 - f) + var2(:,nU)      *f
    v_geos_t (:) = var1(:,nV)      *(1.0d0 - f) + var2(:,nV)      *f
    momz_ls_t(:) = var1(:,nW_ls)   *(1.0d0 - f) + var2(:,nW_ls)   *f
    pott_ls_t(:) = var1(:,nPT_tend)*(1.0d0 - f) + var2(:,nPT_tend)*f 
    qv_ls_t  (:) = var1(:,nQV_tend)*(1.0d0 - f) + var2(:,nQV_tend)*f
    do k=2, ka-1
!     momz_ls_dz_t(k) = (momz_ls_t(k+1)-momz_ls_t(k-1))/(cz(k+1)+cz(k-1))
      momz_ls_dz_t(k) = (momz_ls_t(k+1)-momz_ls_t(k-1))/(cz(k+1)-cz(k-1))
    enddo
    momz_ls_dz_t(1) = momz_ls_t(2)/cz(2)
    momz_ls_dz_t(ka) = momz_ls_t(ka-1)

    !
    return
    !
  end subroutine update_var
  !--------------------------------------------------------------------------------
  subroutine read_var(var)
    !
    implicit none
    !
    real(8), intent(out) :: var(1:ka,varmax)
!   real(4) :: tmp(1:ka,varmax)
    !
    integer :: kk, k, n
    !
    var = 0.0d0
    !
    read(fid_data,*) kk
    do k=1, ka
      read(fid_data,*) (wk(k,n),n=1,varmax)
    enddo
    !
    var(:,:)=dble(wk(:,:))

    !
    return
    !
  end subroutine read_var

end module mod_user
