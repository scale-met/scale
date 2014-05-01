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
  use mod_stdio, only: &
     IO_get_available_fid

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
  real(DP),allocatable :: momz_ls_t(:)
  real(DP),allocatable :: momz_ls_dz_t(:)
  real(DP),allocatable :: z_in(:)
  real(DP),save,allocatable :: time_atm_in(:)
  real(DP),save,allocatable :: time_sst_in(:)
  real(DP),save,allocatable :: sst_in(:)
! real(DP),allocatable :: u_geos_t(:)
! real(DP),allocatable :: v_geos_t(:)
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

  integer,  private, save :: USER_LS_FLG = 0 !-- 0->no force, 1->TWPICE
  real(RP), private, save :: corioli

  real(RP), private, save :: CNST_SST=276.2_RP

  character(100), private, save :: inbasedir = './'
  character(100), private, save :: fdata_name_atm = 'large_scale_w_force.txt'
  character(100), private, save :: fdata_name_sst = 'sst_force.txt'
  logical, public, save :: CNST_RAD=.false. ! add constant radiative cooling
  integer, private, save :: fid_data
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
         inbasedir,     &
         fdata_name_atm,&
         fdata_name_sst,&
         mstep_atm,         &
         mstep_sst,         &
         CNST_RAD,&
         CNST_SST,&
       USER_LS_FLG

    integer :: ierr, t
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

if(io_l)write(IO_FID_LOG,*)'debug stop'  ! ok
!!call PRC_MPIstop
!call PRC_MPIfinish

    allocate( time_sst_in(mstep_sst) )
    allocate( sst_in(mstep_sst) )

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

    do t=1, mstep_sst-1
        if( time_nowsec>=time_sst_in(t) )then
          sst(:,:)=( (time_sst_in(t+1)-time_nowsec)*sst_in(t)+(time_nowsec-time_sst_in(t))*sst_in(t+1) )&
                   /(time_sst_in(t+1)-time_sst_in(t))
          exit
        endif
    enddo
!   SST(:,:)=CNST_SST ! initialization

!dbg
!    call PRC_MPIstop

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

    return
  end subroutine USER_step

  !---------------------------------------------------------------------------------
end module mod_user
