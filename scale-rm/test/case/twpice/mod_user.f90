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
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: varmax   = 5
  integer, private, parameter :: nU       = 1
  integer, private, parameter :: nV       = 2
  integer, private, parameter :: nPT_tend = 3  ! tendency of pot.temp [K/s]
  integer, private, parameter :: nQV_tend = 4  ! tendency of qv [kg/kg/s]
  integer, private, parameter :: nW_ls    = 5  ! large-scale vertical velocity [m/s]

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

  real(RP), private, save :: CNST_SST = 302.15_RP

  character(len=H_LONG), private, save :: inbasedir = './'
  character(len=H_LONG), private, save :: fdata_name = 'forcing4run.txt'
  integer, private, save :: mstep = 1  ! max time step in fdata_name
  integer, private, save :: intv   = -999

  real(RP), private, save :: start_hr = 0.0_RP
  logical, public, save :: CNST_RAD=.false. ! add constant radiative cooling
  integer, private, save :: start_step=-999
  integer, private, save :: fid_data
  logical, private, save :: rd1st=.true.


  integer, private, save :: mtnum
  integer, private, save :: stepv1
  integer, private, save :: stepv2

  real(RP), allocatable, private, save :: var1(:,:)
  real(RP), allocatable, private, save :: var2(:,:)
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

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       TIME_DTSEC
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

    real(RP) :: wk(KA,varmax)
    integer  :: kk

    integer  :: ierr
    integer  :: k, n, mt
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    !--- equivalent time step for the current dt
    if( start_step == -999 ) start_step = int(start_hr/TIME_DTSEC)+1
    !--- assume 3hourly external forcing data
    if( intv == -999 ) intv = int(3.0*3600/TIME_DTSEC)

    allocate( MOMZ_LS(KA,2) )
    allocate( MOMZ_LS_DZ(KA,2) )
    allocate( U_GEOS(KA) )
    allocate( V_GEOS(KA) )
    allocate( QV_LS(KA,2) )

    allocate( var1(KA,varmax) )
    allocate( var2(KA,varmax) )

    allocate( momz_ls_t(KA) )
    allocate( momz_ls_dz_t(KA) )
    allocate( u_geos_t(KA) )
    allocate( v_geos_t(KA) )
    allocate( qv_ls_t(KA) )
    allocate( pott_ls_t(KA) )

    var1(:,:) = 0.0_RP
    var2(:,:) = 0.0_RP

    ! open 1-dim forcing data
    fid_data = IO_get_available_fid()
    fdata_name=trim(inbasedir)//'/'//trim(fdata_name)
    open(fid_data, file=trim(fdata_name), status='old',iostat=ierr)
    if ( ierr /= 0 ) then
       write(*,*) 'Cannot open the data file for forcing. STOP! ', trim(fdata_name)
       call PRC_MPIstop
    endif

    !--- Check if whole input data available
    do mt = 1, mstep
       read(fid_data,*) kk

       do k = 1, KA
          read(fid_data,*,iostat=ierr) (wk(k,n),n=1,varmax)

          if ( ierr /= 0 ) then
             write(*,*) 'Not enough data! ',mt, k, mstep, KA, trim(fdata_name)
             call PRC_MPIstop
          endif
       enddo
    enddo
    rewind(fid_data)

    mtnum = int(start_hr*3600.0/TIME_DTSEC/intv) + 1

    if( IO_L ) write(IO_FID_LOG,*) 'Forcing starts! ',mtnum

    ! skip by the adequate time step
    do mt = 1, mtnum-1
       read(fid_data,*) kk
       do k = 1, KA
          read(fid_data,*) (wk(k,n),n=1,varmax)
       enddo
    enddo

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    use mod_cpl_vars, only: &
       SST
    implicit none
    !---------------------------------------------------------------------------

    SST(:,:) = CNST_SST

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
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
    use mod_cpl_vars, only: &
       SST
    implicit none

    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: PRES(KA,IA,JA)
    real(RP) :: TEMP(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA)
    real(RP) :: VELY(KA,IA,JA)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( .not. USER_do ) then
       return
    endif

    SST(:,:) = CNST_SST


    call update_var


    if ( USER_LS_FLG == 0 ) then  ! no large scale sinking

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
          endif

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
          endif

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
          endif

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

             if ( CNST_RAD ) then
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

          endif

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

          endif

       enddo
       enddo

    else
       write(*,*)'Not supported user_ls_flg'
       call PRC_MPIstop
    endif

    return
  end subroutine USER_step

  !---------------------------------------------------------------------------------
  subroutine update_var
    use scale_time, only: &
       TIME_DTSEC, &
       TIME_NOWSTEP
    use scale_grid, only: &
       CZ => GRID_CZ
    implicit none

    real(RP) :: f

    integer :: cstep, k
    !---------------------------------------------------------------------------

    ! current time step counted from the initial run
    cstep = int(start_hr*3600/TIME_DTSEC) + TIME_NOWSTEP + 1

    if ( rd1st ) then

       stepv1 = (mtnum-1)*intv + 1
       stepv2 = stepv1 + intv

       call read_var(var1(:,:))
       mtnum = mtnum + 1
       call read_var(var2(:,:))
       mtnum = mtnum + 1

       rd1st = .false.
    else
       if ( mod((cstep-start_step),intv) == 0 ) then
          var1(:,:)    = var2(:,:)
          call read_var(var2(:,:))
          mtnum = mtnum + 1

          stepv1 = stepv2
          stepv2 = stepv2 + intv
       endif
    endif

    f = real(cstep-stepv1,kind=RP) / real(intv,kind=RP)

!   var_bs(:,:) = var1(:,:)*(1.0_RP - f) + var2(:,:)*f
    u_geos_t (:) = var1(:,nU)       * (1.0_RP-f) + var2(:,nU)       * f
    v_geos_t (:) = var1(:,nV)       * (1.0_RP-f) + var2(:,nV)       * f
    momz_ls_t(:) = var1(:,nW_ls)    * (1.0_RP-f) + var2(:,nW_ls)    * f
    pott_ls_t(:) = var1(:,nPT_tend) * (1.0_RP-f) + var2(:,nPT_tend) * f
    qv_ls_t  (:) = var1(:,nQV_tend) * (1.0_RP-f) + var2(:,nQV_tend) * f

    do k=2, KA-1
      momz_ls_dz_t(k) = (momz_ls_t(k+1)-momz_ls_t(k-1))/(cz(k+1)-cz(k-1))
    enddo
    momz_ls_dz_t(1)  = momz_ls_t(2) / cz(2)
    momz_ls_dz_t(KA) = momz_ls_t(KA-1)

    return
  end subroutine update_var

  !-----------------------------------------------------------------------------
  subroutine read_var(var)
    implicit none

    real(RP), intent(out) :: var(KA,varmax)

    integer :: kk, k, n
    !---------------------------------------------------------------------------

    var = 0.0_RP

    read(fid_data,*) kk
    do k = 1, KA
      read(fid_data,*) (var(k,n),n=1,varmax)
    enddo

    return
  end subroutine read_var

end module mod_user
