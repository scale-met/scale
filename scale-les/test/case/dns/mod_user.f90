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

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

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

    if ( USER_do ) then
! under construction
      return
    endif

    return
  end subroutine USER_step

  !---------------------------------------------------------------------------------

end module mod_user
