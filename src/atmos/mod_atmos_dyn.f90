!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FENT + FCT
!!
!! @par Description
!!          Dynamical core for Atmospheric process
!!          Full explicit, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)   [new] Imported from SCALE-LES ver.2
!! @li      2011-11-11 (H.Yashiro)   [mod] Merged with Y.Miyamoto's
!! @li      2011-12-11 (H.Yashiro)   [mod] Use reference state
!! @li      2011-12-26 (Y.Miyamoto)  [mod] Add numerical diffusion into mass flux calc
!! @li      2012-01-04 (H.Yashiro)   [mod] Nonblocking communication (Y.Ohno)
!! @li      2012-01-25 (H.Yashiro)   [fix] Bugfix (Y.Miyamoto)
!! @li      2011-01-25 (H.Yashiro)   [mod] sprit as "full" FCT (Y.Miyamoto)
!! @li      2012-02-14 (H.Yashiro)   [mod] Cache tiling
!! @li      2012-03-14 (H.Yashiro)   [mod] Bugfix (Y.Miyamoto)
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-04-09 (H.Yashiro)   [mod] Integrate RDMA communication
!! @li      2012-06-10 (Y.Miyamoto)  [mod] large-scale divergence (from H.Yashiro's)
!! @li      2012-07-13 (H.Yashiro)   [mod] prevent conditional branching in FCT
!! @li      2012-07-27 (Y.Miyamoto)  [mod] divegence damping option
!! @li      2012-08-16 (S.Nishizawa) [mod] use FCT for momentum and temperature
!! @li      2012-09-21 (Y.Sato)      [mod] merge DYCOMS-II experimental set
!! @li      2013-03-26 (Y.Sato)      [mod] modify Large scale forcing and corioli forcing
!! @li      2013-04-04 (Y.Sato)      [mod] modify Large scale forcing
!! @li      2013-06-14 (S.Nishizawa) [mod] enable to change order of numerical diffusion
!! @li      2013-06-18 (S.Nishizawa) [mod] split part of RK to other files
!! @li      2013-06-20 (S.Nishizawa) [mod] split large scale sining to other file
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
#ifdef DEBUG
  use mod_debug, only: &
     CHECK
  use mod_const, only: &
     UNDEF => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  use mod_atmos_vars, only: &
     ZDIR, &
     XDIR, &
     YDIR, &
     I_DENS, &
     I_MOMZ, &
     I_MOMX, &
     I_MOMY, &
     I_RHOT
  use mod_atmos_dyn_common, only: &
     FACT_N, &
     FACT_F

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_setup
  public :: ATMOS_DYN
#ifdef DEBUG
  public :: ATMOS_DYN_init
  public :: ATMOS_DYN_main
#endif

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  ! time integration scheme
  integer,  private, parameter :: RK = 3             ! order of Runge-Kutta scheme
  real(RP), private, save      :: DENS_RK1(KA,IA,JA) ! prognostic variables (+1/3 step)
  real(RP), private, save      :: MOMZ_RK1(KA,IA,JA) !
  real(RP), private, save      :: MOMX_RK1(KA,IA,JA) !
  real(RP), private, save      :: MOMY_RK1(KA,IA,JA) !
  real(RP), private, save      :: RHOT_RK1(KA,IA,JA) !
  real(RP), private, save      :: DENS_RK2(KA,IA,JA) ! prognostic variables (+2/3 step)
  real(RP), private, save      :: MOMZ_RK2(KA,IA,JA) !
  real(RP), private, save      :: MOMX_RK2(KA,IA,JA) !
  real(RP), private, save      :: MOMY_RK2(KA,IA,JA) !
  real(RP), private, save      :: RHOT_RK2(KA,IA,JA) !

  ! numerical filter
  integer,  private, save      :: ATMOS_DYN_numerical_diff_order        = 1
  real(RP), private, save      :: ATMOS_DYN_numerical_diff_coef         = 1.0E-4_RP ! nondimensional numerical diffusion
  real(RP), private, save      :: ATMOS_DYN_numerical_diff_sfc_fact     = 1.0_RP
  logical , private, save      :: ATMOS_DYN_numerical_diff_use_refstate = .true.
  real(RP), private, save      :: ATMOS_DYN_DIFF4                                   ! for numerical filter
  real(RP), private, save      :: ATMOS_DYN_CNZ3(3,KA,2)
  real(RP), private, save      :: ATMOS_DYN_CNX3(3,IA,2)
  real(RP), private, save      :: ATMOS_DYN_CNY3(3,JA,2)
  real(RP), private, save      :: ATMOS_DYN_CNZ4(5,KA,2)
  real(RP), private, save      :: ATMOS_DYN_CNX4(5,IA,2)
  real(RP), private, save      :: ATMOS_DYN_CNY4(5,JA,2)

  ! Coriolis force
  logical,  private, save      :: ATMOS_DYN_enable_coriolis = .false. ! enable coriolis force?
  real(RP), private, save      :: ATMOS_DYN_CORIOLI(IA,JA)            ! coriolis term

  ! divergence damping
  real(RP), private, save      :: ATMOS_DYN_divdmp_coef = 0.0_RP      ! Divergence dumping coef

  ! fct
  logical,  private, save      :: ATMOS_DYN_FLAG_FCT_rho      = .false.
  logical,  private, save      :: ATMOS_DYN_FLAG_FCT_momentum = .false.
  logical,  private, save      :: ATMOS_DYN_FLAG_FCT_T        = .false.

  ! work
  real(RP), private, save      :: mflx_hi(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face high order
  real(RP), private, save      :: mflx_av(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face average
  real(RP), private, save      :: VELZ   (KA,IA,JA)    ! velocity w [m/s]
  real(RP), private, save      :: VELX   (KA,IA,JA)    ! velocity u [m/s]
  real(RP), private, save      :: VELY   (KA,IA,JA)    ! velocity v [m/s]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use mod_grid, only : &
       GRID_CDZ, &
       GRID_CDX, &
       GRID_CDY, &
       GRID_FDZ, &
       GRID_FDX, &
       GRID_FDY
    use mod_geometrics, only : &
       GEOMETRICS_lat
    use mod_atmos_dyn_rk, only: &
       ATMOS_DYN_rk_setup
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff_order,        &
       ATMOS_DYN_numerical_diff_coef,         &
       ATMOS_DYN_numerical_diff_sfc_fact,     &
       ATMOS_DYN_numerical_diff_use_refstate, &
       ATMOS_DYN_enable_coriolis,             &
       ATMOS_DYN_divdmp_coef,                 &
       ATMOS_DYN_FLAG_FCT_rho,                &
       ATMOS_DYN_FLAG_FCT_momentum,           &
       ATMOS_DYN_FLAG_FCT_T

    real(RP) :: DT

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Dynamics]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN)

    !-- Block size must be divisible
    if    ( mod(IMAX,IBLOCK) > 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx number of grid size IMAX must be divisible by IBLOCK! ', IMAX, IBLOCK
       call PRC_MPIstop
    elseif( mod(JMAX,JBLOCK) > 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx number of grid size JMAX must be divisible by JBLOCK! ', JMAX, JBLOCK
       call PRC_MPIstop
    endif

    DT = real(TIME_DTSEC_ATMOS_DYN,kind=RP)

    call ATMOS_DYN_init( ATMOS_DYN_DIFF4,                & ! [OUT]
                         ATMOS_DYN_CNZ3,                 & ! [OUT]
                         ATMOS_DYN_CNX3,                 & ! [OUT]
                         ATMOS_DYN_CNY3,                 & ! [OUT]
                         ATMOS_DYN_CNZ4,                 & ! [OUT]
                         ATMOS_DYN_CNX4,                 & ! [OUT]
                         ATMOS_DYN_CNY4,                 & ! [OUT]
                         ATMOS_DYN_CORIOLI,              & ! [OUT]
                         GRID_CDZ, GRID_CDX, GRID_CDY,   & ! [IN]
                         GRID_FDZ, GRID_FDX, GRID_FDY,   & ! [IN]
                         ATMOS_DYN_numerical_diff_order, & ! [IN]
                         ATMOS_DYN_numerical_diff_coef,  & ! [IN]
                         DT,                             & ! [IN]
                         ATMOS_DYN_enable_coriolis,      & ! [IN]
                         GEOMETRICS_lat                  ) ! [IN]

    call ATMOS_DYN_rk_setup

    return
  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_init( &
       DIFF4,            &
       CNZ3, CNX3, CNY3, &
       CNZ4, CNX4, CNY4, &
       CORIOLI,          &
       CDZ, CDX, CDY,    &
       FDZ, FDX, FDY,    &
       numdiff_order,    &
       numdiff_coef,     &
       dt_dyn,           &
       enable_coriolis,  &
       lat               )
    use mod_const, only: &
       PI  => CONST_PI, &
       OHM => CONST_OHM
    implicit none

    real(RP), intent(out) :: DIFF4
    real(RP), intent(out) :: CNZ3(3,KA,2)
    real(RP), intent(out) :: CNX3(3,IA,2)
    real(RP), intent(out) :: CNY3(3,JA,2)
    real(RP), intent(out) :: CNZ4(5,KA,2)
    real(RP), intent(out) :: CNX4(5,IA,2)
    real(RP), intent(out) :: CNY4(5,JA,2)
    real(RP), intent(out) :: CORIOLI(1,IA,JA)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)
    real(RP), intent(in)  :: FDZ(KA-1)
    real(RP), intent(in)  :: FDX(IA-1)
    real(RP), intent(in)  :: FDY(JA-1)
    integer,  intent(in)  :: numdiff_order
    real(RP), intent(in)  :: numdiff_coef
    real(RP), intent(in)  :: dt_dyn
    logical , intent(in)  :: enable_coriolis
    real(RP), intent(in)  :: lat(1,IA,JA)

    integer :: i, j, k
    !---------------------------------------------------------------------------

#ifdef DEBUG
    CNZ3(:,:,:) = UNDEF
    CNX3(:,:,:) = UNDEF
    CNY3(:,:,:) = UNDEF
    CNZ4(:,:,:) = UNDEF
    CNX4(:,:,:) = UNDEF
    CNY4(:,:,:) = UNDEF
    corioli(:,:,:) = UNDEF
#endif

    ! numerical diffusion
    if ( numdiff_coef > 0.0_RP ) then
       DIFF4 = numdiff_coef / ( 2**(4*numdiff_order) * dt_dyn )

       ! z direction
       do k = KS+1, KE-1
          CNZ3(1,k,1) = 1.0_RP / ( FDZ(k) * CDZ(k) * FDZ(k-1) )
          CNZ3(2,k,1) = 1.0_RP / ( FDZ(k  ) * CDZ(k  ) * FDZ(k-1) ) &
                      + 1.0_RP / ( FDZ(k-1) * CDZ(k  ) * FDZ(k-1) ) &
                      + 1.0_RP / ( FDZ(k-1) * CDZ(k-1) * FDZ(k-1) )
       enddo
       do k = KS+2, KE
          CNZ3(3,k,1) = 1.0_RP / ( FDZ(k-1) * CDZ(k  ) * FDZ(k-1) ) &
                      + 1.0_RP / ( FDZ(k-1) * CDZ(k-1) * FDZ(k-1) ) &
                      + 1.0_RP / ( FDZ(k-1) * CDZ(k-1) * FDZ(k-2) )
       enddo
       CNZ3(1,KS-1,1) = 1.0_RP / ( FDZ(KS  ) * CDZ(KS+1) * FDZ(KS  ) )
       CNZ3(1,KS  ,1) = 1.0_RP / ( FDZ(KS  ) * CDZ(KS  ) * FDZ(KS  ) )
       CNZ3(2,KS  ,1) = 1.0_RP / ( FDZ(KS  ) * CDZ(KS  ) * FDZ(KS  ) ) &
                      + 1.0_RP / ( FDZ(KS  ) * CDZ(KS  ) * FDZ(KS  ) ) &
                      + 1.0_RP / ( FDZ(KS  ) * CDZ(KS+1) * FDZ(KS  ) )
       CNZ3(3,KS  ,1) = 1.0_RP / ( FDZ(KS  ) * CDZ(KS  ) * FDZ(KS  ) ) &
                      + 1.0_RP / ( FDZ(KS  ) * CDZ(KS+1) * FDZ(KS  ) ) &
                      + 1.0_RP / ( FDZ(KS  ) * CDZ(KS+1) * FDZ(KS+1) )
       CNZ3(3,KS+1,1) = 1.0_RP / ( FDZ(KS  ) * CDZ(KS+1) * FDZ(KS  ) ) &
                      + 1.0_RP / ( FDZ(KS  ) * CDZ(KS  ) * FDZ(KS  ) ) &
                      + 1.0_RP / ( FDZ(KS  ) * CDZ(KS  ) * FDZ(KS  ) )
       CNZ3(1,KE  ,1) = 1.0_RP / ( FDZ(KE-1) * CDZ(KE  ) * FDZ(KE-1) )
       CNZ3(2,KE  ,1) = 1.0_RP / ( FDZ(KE-1) * CDZ(KE  ) * FDZ(KE-1) ) &
                      + 1.0_RP / ( FDZ(KE-1) * CDZ(KE  ) * FDZ(KE-1) ) &
                      + 1.0_RP / ( FDZ(KE-1) * CDZ(KE-1) * FDZ(KE-1) )
       CNZ3(1,KE+1,1) = 1.0_RP / ( FDZ(KE-2) * CDZ(KE-1) * FDZ(KE-1) )
       CNZ3(2,KE+1,1) = 1.0_RP / ( FDZ(KE-2) * CDZ(KE-1) * FDZ(KE-1) ) &
                      + 1.0_RP / ( FDZ(KE-1) * CDZ(KE-1) * FDZ(KE-1) ) &
                      + 1.0_RP / ( FDZ(KE-1) * CDZ(KE  ) * FDZ(KE-1) )
       CNZ3(3,KE+1,1) = 1.0_RP / ( FDZ(KE-1) * CDZ(KE+1) * FDZ(KE-1) ) &
                      + 1.0_RP / ( FDZ(KE-1) * CDZ(KE  ) * FDZ(KE-1) ) &
                      + 1.0_RP / ( FDZ(KE-1) * CDZ(KE  ) * FDZ(KE-1) )

       do k = KS, KE-1
          CNZ4(1,k,1) = ( CNZ3(1,k+1,1)               ) / CDZ(k)
          CNZ4(2,k,1) = ( CNZ3(2,k+1,1) + CNZ3(1,k,1) ) / CDZ(k)
          CNZ4(3,k,1) = ( CNZ3(3,k+1,1) + CNZ3(2,k,1) ) / CDZ(k)
          CNZ4(4,k,1) = ( CNZ3(1,k  ,1) + CNZ3(3,k,1) ) / CDZ(k)
          CNZ4(5,k,1) = ( CNZ3(1,k-1,1)               ) / CDZ(k)
       enddo

       do k = KS+1, KE-1
          CNZ3(1,k,2) = 1.0_RP / ( CDZ(k+1) * FDZ(k  ) * CDZ(k  ) )
          CNZ3(2,k,2) = 1.0_RP / ( CDZ(k+1) * FDZ(k  ) * CDZ(k  ) ) &
                      + 1.0_RP / ( CDZ(k  ) * FDZ(k  ) * CDZ(k  ) ) &
                      + 1.0_RP / ( CDZ(k  ) * FDZ(k-1) * CDZ(k  ) )
          CNZ3(3,k,2) = 1.0_RP / ( CDZ(k  ) * FDZ(k  ) * CDZ(k  ) ) &
                      + 1.0_RP / ( CDZ(k  ) * FDZ(k-1) * CDZ(k  ) ) &
                      + 1.0_RP / ( CDZ(k  ) * FDZ(k-1) * CDZ(k-1) )
       enddo
       CNZ3(1,KS-1,2) = 1.0_RP / ( CDZ(KS  ) * FDZ(KS  ) * CDZ(KS  ) )
       CNZ3(1,KS  ,2) = 1.0_RP / ( CDZ(KS+1) * FDZ(KS  ) * CDZ(KS  ) )
       CNZ3(2,KS  ,2) = 1.0_RP / ( CDZ(KS+1) * FDZ(KS  ) * CDZ(KS  ) ) &
                      + 1.0_RP / ( CDZ(KS  ) * FDZ(KS  ) * CDZ(KS  ) ) &
                      + 1.0_RP / ( CDZ(KS  ) * FDZ(KS  ) * CDZ(KS  ) )
       CNZ3(3,KS  ,2) = 1.0_RP / ( CDZ(KS  ) * FDZ(KS  ) * CDZ(KS  ) ) &
                      + 1.0_RP / ( CDZ(KS  ) * FDZ(KS  ) * CDZ(KS  ) ) &
                      + 1.0_RP / ( CDZ(KS  ) * FDZ(KS  ) * CDZ(KS+1) )
       CNZ3(1,KE  ,2) = 1.0_RP / ( CDZ(KE-1) * FDZ(KE-1) * CDZ(KE  ) )
       CNZ3(2,KE  ,2) = 1.0_RP / ( CDZ(KE-1) * FDZ(KE-1) * CDZ(KE  ) ) &
                      + 1.0_RP / ( CDZ(KE  ) * FDZ(KE-1) * CDZ(KE  ) ) &
                      + 1.0_RP / ( CDZ(KE  ) * FDZ(KE-1) * CDZ(KE  ) )
       CNZ3(3,KE  ,2) = 1.0_RP / ( CDZ(KE  ) * FDZ(KE-1) * CDZ(KE  ) ) &
                      + 1.0_RP / ( CDZ(KE  ) * FDZ(KE-1) * CDZ(KE  ) ) &
                      + 1.0_RP / ( CDZ(KE  ) * FDZ(KE-1) * CDZ(KE-1) )
       CNZ3(1,KE+1,2) = 1.0_RP / ( CDZ(KE-2) * FDZ(KE-2) * CDZ(KE-1) )
       CNZ3(2,KE+1,2) = 1.0_RP / ( CDZ(KE-2) * FDZ(KE-2) * CDZ(KE-1) ) &
                      + 1.0_RP / ( CDZ(KE-1) * FDZ(KE-2) * CDZ(KE-1) ) &
                      + 1.0_RP / ( CDZ(KE-1) * FDZ(KE-1) * CDZ(KE-1) )
       CNZ3(3,KE+1,2) = 1.0_RP / ( CDZ(KE-1) * FDZ(KE-2) * CDZ(KE-1) ) &
                      + 1.0_RP / ( CDZ(KE-1) * FDZ(KE-1) * CDZ(KE-1) ) &
                      + 1.0_RP / ( CDZ(KE-1) * FDZ(KE-1) * CDZ(KE  ) )

       do k = KS, KE-1
          CNZ4(1,k,2) = ( CNZ3(1,k+1,2)               ) / FDZ(k)
          CNZ4(2,k,2) = ( CNZ3(2,k+1,2) + CNZ3(1,k,2) ) / FDZ(k)
          CNZ4(3,k,2) = ( CNZ3(3,k+1,2) + CNZ3(2,k,2) ) / FDZ(k)
          CNZ4(4,k,2) = ( CNZ3(1,k  ,2) + CNZ3(3,k,2) ) / FDZ(k)
          CNZ4(5,k,2) = ( CNZ3(1,k-1,2)               ) / FDZ(k)
       enddo
!       CNZ4(1,KE,2) = ( CNZ3(1,KE+1,2)                ) / FDZ(KE-1)
       CNZ4(2,KE,2) = ( CNZ3(2,KE+1,2) + CNZ3(1,KE,2) ) / FDZ(KE-1)
       CNZ4(3,KE,2) = ( CNZ3(3,KE+1,2) + CNZ3(2,KE,2) ) / FDZ(KE-1)
       CNZ4(4,KE,2) = ( CNZ3(1,KE  ,2) + CNZ3(3,KE,2) ) / FDZ(KE-1)

       ! x direction
       CNX3(1,IS-1,1) = 1.0_RP / ( FDX(IS-1) * CDX(IS-1) * FDX(IS-2) )
       do i = IS, IE+1
          CNX3(1,i,1) = 1.0_RP / ( FDX(i  ) * CDX(i  ) * FDX(i-1) )
          CNX3(2,i,1) = 1.0_RP / ( FDX(i  ) * CDX(i  ) * FDX(i-1) ) &
                      + 1.0_RP / ( FDX(i-1) * CDX(i  ) * FDX(i-1) ) &
                      + 1.0_RP / ( FDX(i-1) * CDX(i-1) * FDX(i-1) )
          CNX3(3,i,1) = 1.0_RP / ( FDX(i-1) * CDX(i  ) * FDX(i-1) ) &
                      + 1.0_RP / ( FDX(i-1) * CDX(i-1) * FDX(i-1) ) &
                      + 1.0_RP / ( FDX(i-1) * CDX(i-1) * FDX(i-2) )
       enddo

       do i = IS, IE
          CNX4(1,i,1) = ( CNX3(1,i+1,1)               ) / CDX(i)
          CNX4(2,i,1) = ( CNX3(2,i+1,1) + CNX3(1,i,1) ) / CDX(i)
          CNX4(3,i,1) = ( CNX3(3,i+1,1) + CNX3(2,i,1) ) / CDX(i)
          CNX4(4,i,1) = ( CNX3(1,i  ,1) + CNX3(3,i,1) ) / CDX(i)
          CNX4(5,i,1) = ( CNX3(1,i-1,1)               ) / CDX(i)
       enddo

       do i = IS-1, IE+1
          CNX3(1,i,2) = 1.0_RP / ( CDX(i+1) * FDX(i  ) * CDX(i  ) )
          CNX3(2,i,2) = 1.0_RP / ( CDX(i+1) * FDX(i  ) * CDX(i  ) ) &
                      + 1.0_RP / ( CDX(i  ) * FDX(i  ) * CDX(i  ) ) &
                      + 1.0_RP / ( CDX(i  ) * FDX(i-1) * CDX(i  ) )
          CNX3(3,i,2) = 1.0_RP / ( CDX(i  ) * FDX(i  ) * CDX(i  ) ) &
                      + 1.0_RP / ( CDX(i  ) * FDX(i-1) * CDX(i  ) ) &
                      + 1.0_RP / ( CDX(i  ) * FDX(i-1) * CDX(i-1) )
       enddo

       do i = IS, IE
          CNX4(1,i,2) = ( CNX3(1,i+1,2)               ) / FDX(i)
          CNX4(2,i,2) = ( CNX3(2,i+1,2) + CNX3(1,i,2) ) / FDX(i)
          CNX4(3,i,2) = ( CNX3(3,i+1,2) + CNX3(2,i,2) ) / FDX(i)
          CNX4(4,i,2) = ( CNX3(1,i  ,2) + CNX3(3,i,2) ) / FDX(i)
          CNX4(5,i,2) = ( CNX3(1,i-1,2)               ) / FDX(i)
       enddo

       ! y direction
       CNY3(1,JS-1,1) = 1.0_RP / ( FDY(JS-1) * CDY(JS-1) * FDY(JS-2) )
       do j = JS, JE+1
          CNY3(1,j,1) = 1.0_RP / ( FDY(j  ) * CDY(j  ) * FDY(j-1) )
          CNY3(2,j,1) = 1.0_RP / ( FDY(j  ) * CDY(j  ) * FDY(j-1) ) &
                      + 1.0_RP / ( FDY(j-1) * CDY(j  ) * FDY(j-1) ) &
                      + 1.0_RP / ( FDY(j-1) * CDY(j-1) * FDY(j-1) )
          CNY3(3,j,1) = 1.0_RP / ( FDY(j-1) * CDY(j  ) * FDY(j-1) ) &
                      + 1.0_RP / ( FDY(j-1) * CDY(j-1) * FDY(j-1) ) &
                      + 1.0_RP / ( FDY(j-1) * CDY(j-1) * FDY(j-2) )
       enddo

       do j = JS, JE
          CNY4(1,j,1) = ( CNY3(1,j+1,1)               ) / CDY(j)
          CNY4(2,j,1) = ( CNY3(2,j+1,1) + CNY3(1,j,1) ) / CDY(j)
          CNY4(3,j,1) = ( CNY3(3,j+1,1) + CNY3(2,j,1) ) / CDY(j)
          CNY4(4,j,1) = ( CNY3(1,j  ,1) + CNY3(3,j,1) ) / CDY(j)
          CNY4(5,j,1) = ( CNY3(1,j-1,1)               ) / CDY(j)
       enddo

       do j = JS-1, JE+1
          CNY3(1,j,2) = 1.0_RP / ( CDY(j+1) * FDY(j  ) * CDY(j  ) )
          CNY3(2,j,2) = 1.0_RP / ( CDY(j+1) * FDY(j  ) * CDY(j  ) ) &
                      + 1.0_RP / ( CDY(j  ) * FDY(j  ) * CDY(j  ) ) &
                      + 1.0_RP / ( CDY(j  ) * FDY(j-1) * CDY(j  ) )
          CNY3(3,j,2) = 1.0_RP / ( CDY(j  ) * FDY(j  ) * CDY(j  ) ) &
                      + 1.0_RP / ( CDY(j  ) * FDY(j-1) * CDY(j  ) ) &
                      + 1.0_RP / ( CDY(j  ) * FDY(j-1) * CDY(j-1) )
       enddo

       do j = JS, JE
          CNY4(1,j,2) = ( CNY3(1,j+1,2)               ) / FDY(j)
          CNY4(2,j,2) = ( CNY3(2,j+1,2) + CNY3(1,j,2) ) / FDY(j)
          CNY4(3,j,2) = ( CNY3(3,j+1,2) + CNY3(2,j,2) ) / FDY(j)
          CNY4(4,j,2) = ( CNY3(1,j  ,2) + CNY3(3,j,2) ) / FDY(j)
          CNY4(5,j,2) = ( CNY3(1,j-1,2)               ) / FDY(j)
       enddo
    else
       DIFF4 = 0.0_RP

       CNZ3(:,:,:) = 0.0_RP
       CNX3(:,:,:) = 0.0_RP
       CNY3(:,:,:) = 0.0_RP
       CNZ4(:,:,:) = 0.0_RP
       CNX4(:,:,:) = 0.0_RP
       CNY4(:,:,:) = 0.0_RP
    endif

    ! coriolis parameter
    if ( enable_coriolis ) then
       CORIOLI(1,:,:) = 2.0_RP * OHM * sin( lat(1,:,:) )
    else
       CORIOLI(1,:,:) = 0.0_RP
    endif

    return
  end subroutine ATMOS_DYN_init

  !-----------------------------------------------------------------------------
  !> Dynamical Process (Wrapper)
  subroutine ATMOS_DYN
    use mod_time, only: &
       DTSEC           => TIME_DTSEC,           &
       DTSEC_ATMOS_DYN => TIME_DTSEC_ATMOS_DYN, &
       NSTEP_ATMOS_DYN => TIME_NSTEP_ATMOS_DYN
    use mod_grid, only : &
       CDZ  => GRID_CDZ,  &
       CDX  => GRID_CDX,  &
       CDY  => GRID_CDY,  &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_total,  &
       ATMOS_USE_AVERAGE, &
       DENS,    &
       MOMZ,    &
       MOMX,    &
       MOMY,    &
       RHOT,    &
       QTRC,    &
       DENS_av, &
       MOMZ_av, &
       MOMX_av, &
       MOMY_av, &
       RHOT_av, &
       QTRC_av, &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp
    use mod_atmos_thermodyn, only: &
       AQ_CV
    use mod_atmos_refstate, only: &
       REF_dens => ATMOS_REFSTATE_dens, &
       REF_pott => ATMOS_REFSTATE_pott, &
       REF_qv   => ATMOS_REFSTATE_qv
    use mod_atmos_boundary, only: &
       DAMP_var   => ATMOS_BOUNDARY_var,   &
       DAMP_alpha => ATMOS_BOUNDARY_alpha
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics step'

    call ATMOS_DYN_main( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,                   & ! (inout)
         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (inout)
         DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, & ! (in)
         ATMOS_DYN_CNZ3,                                       & ! (in)
         ATMOS_DYN_CNX3,                                       & ! (in)
         ATMOS_DYN_CNY3,                                       & ! (in)
         ATMOS_DYN_CNZ4,                                       & ! (in)
         ATMOS_DYN_CNX4,                                       & ! (in)
         ATMOS_DYN_CNY4,                                       & ! (in)
         CDZ, CDX, CDY, FDZ, FDX, FDY,                         & ! (in)
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   & ! (in)
         AQ_CV,                                                & ! (in)
         REF_dens, REF_pott, REF_qv,                           & ! (in)
         ATMOS_DYN_DIFF4,                                      & ! (in)
         ATMOS_DYN_numerical_diff_order,                       & ! (in)
         ATMOS_DYN_numerical_diff_sfc_fact,                    & ! (in)
         ATMOS_DYN_numerical_diff_use_refstate,                & ! (in)
         ATMOS_DYN_CORIOLI,                                    & ! (in)
         DAMP_var, DAMP_alpha,                                 & ! (in)
         ATMOS_DYN_divdmp_coef,                                & ! (in)
         ATMOS_DYN_FLAG_FCT_rho,                               & ! (in)
         ATMOS_DYN_FLAG_FCT_momentum,                          & ! (in)
         ATMOS_DYN_FLAG_FCT_T,                                 & ! (in)
         ATMOS_USE_AVERAGE,                                    & ! (in)
         DTSEC, DTSEC_ATMOS_DYN, NSTEP_ATMOS_DYN               ) ! (in)

    call ATMOS_vars_total

    return
  end subroutine ATMOS_DYN

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN_main( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,                   &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, QTRC_tp, &
       CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4,                   &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       AQ_CV,                                                &
       REF_dens, REF_pott, REF_qv,                           &
       DIFF4, ND_ORDER, ND_SFC_FACT, ND_USE_RS,              &
       corioli, DAMP_var, DAMP_alpha,                        &
       divdmp_coef,                                          &
       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,          &
       USE_AVERAGE,                                          &
       DTSEC, DTSEC_ATMOS_DYN, NSTEP_ATMOS_DYN               )
    use dc_types, only: &
       DP
    use mod_const, only : &
       Rdry   => CONST_Rdry,  &
       Rvap   => CONST_Rvap,  &
       CVdry  => CONST_CVdry
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_dyn_common, only: &
       ATMOS_DYN_fct
    use mod_atmos_dyn_rk, only: &
       ATMOS_DYN_rk
    use mod_atmos_boundary, only: &
       I_BND_VELZ,  &
       I_BND_VELX,  &
       I_BND_VELY,  &
       I_BND_POTT,  &
       I_BND_QV
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP), intent(inout) :: DENS_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_av(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_av(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_av(KA,IA,JA,QA)

    real(RP), intent(in)    :: DENS_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_tp(KA,IA,JA)
    real(RP), intent(in)    :: QTRC_tp(KA,IA,JA,QA)

    real(RP), intent(in)    :: CNZ3(3,KA,2)
    real(RP), intent(in)    :: CNX3(3,IA,2)
    real(RP), intent(in)    :: CNY3(3,JA,2)
    real(RP), intent(in)    :: CNZ4(5,KA,2)
    real(RP), intent(in)    :: CNX4(5,IA,2)
    real(RP), intent(in)    :: CNY4(5,JA,2)

    real(RP), intent(in)    :: CDZ (KA)
    real(RP), intent(in)    :: CDX (IA)
    real(RP), intent(in)    :: CDY (JA)
    real(RP), intent(in)    :: FDZ (KA-1)
    real(RP), intent(in)    :: FDX (IA-1)
    real(RP), intent(in)    :: FDY (JA-1)
    real(RP), intent(in)    :: RCDZ(KA)
    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)
    real(RP), intent(in)    :: RFDZ(KA-1)
    real(RP), intent(in)    :: RFDX(IA-1)
    real(RP), intent(in)    :: RFDY(JA-1)

    real(RP), intent(in)    :: AQ_CV(QQA)

    real(RP), intent(in)    :: REF_dens(KA)
    real(RP), intent(in)    :: REF_pott(KA)
    real(RP), intent(in)    :: REF_qv  (KA)
    real(RP), intent(in)    :: DIFF4
    integer,  intent(in)    :: ND_ORDER
    real(RP), intent(in)    :: ND_SFC_FACT
    logical,  intent(in)    :: ND_USE_RS
    real(RP), intent(in)    :: CORIOLI(1,IA,JA)
    real(RP), intent(in)    :: DAMP_var  (KA,IA,JA,5)
    real(RP), intent(in)    :: DAMP_alpha(KA,IA,JA,5)
    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_FCT_RHO
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T

    logical,  intent(in)    :: USE_AVERAGE

    real(DP), intent(in)    :: DTSEC
    real(DP), intent(in)    :: DTSEC_ATMOS_DYN
    integer , intent(in)    :: NSTEP_ATMOS_DYN

    ! diagnostic variables
    real(RP) :: PRES  (KA,IA,JA) ! pressure [Pa]
    real(RP) :: POTT  (KA,IA,JA) ! potential temperature [K]
    real(RP) :: dens_s(KA,IA,JA) ! saved density
    real(RP) :: QDRY(KA,IA,JA) ! not used
    real(RP) :: Rtot  (KA,IA,JA) ! R for dry air + vapor
    real(RP) :: CVtot (KA,IA,JA) ! CV

    ! rayleigh damping, numerical diffusion
    real(RP) :: dens_diff(KA,IA,JA)     ! anomary of density
    real(RP) :: pott_diff(KA,IA,JA)     ! anomary of rho * pott
    real(RP) :: qv_diff  (KA,IA,JA)     ! anomary of vapor
    real(RP), target :: num_diff (KA,IA,JA,5,3)

    real(RP) :: qflx_hi  (KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order

    ! tendency
    real(RP) :: DENS_t(KA,IA,JA)
    real(RP) :: MOMZ_t(KA,IA,JA)
    real(RP) :: MOMX_t(KA,IA,JA)
    real(RP) :: MOMY_t(KA,IA,JA)
    real(RP) :: RHOT_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)

    ! For FCT
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux
    real(RP) :: RHOQ     (KA,IA,JA)    ! rho(previous) * phi(previous)

    ! For numerical diffusion
    real(RP), target  :: num_diff_work(KA,IA,JA,5,3)
    real(RP), pointer :: num_diff_pt0(:,:,:,:,:)
    real(RP), pointer :: num_diff_pt1(:,:,:,:,:)
    real(RP), pointer :: tmp_pt      (:,:,:,:,:)


    real(RP) :: dtrk
    integer  :: nd_order4, no
    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k, iq, rko, step
    !---------------------------------------------------------------------------

#ifdef DEBUG
    DENS_RK1(:,:,:) = UNDEF
    MOMZ_RK1(:,:,:) = UNDEF
    MOMX_RK1(:,:,:) = UNDEF
    MOMY_RK1(:,:,:) = UNDEF
    RHOT_RK1(:,:,:) = UNDEF
    DENS_RK2(:,:,:) = UNDEF
    MOMZ_RK2(:,:,:) = UNDEF
    MOMX_RK2(:,:,:) = UNDEF
    MOMY_RK2(:,:,:) = UNDEF
    RHOT_RK2(:,:,:) = UNDEF

    VELZ    (:,:,:) = UNDEF
    VELX    (:,:,:) = UNDEF
    VELY    (:,:,:) = UNDEF
    POTT    (:,:,:) = UNDEF
    PRES    (:,:,:) = UNDEF

    dens_s   (:,:,:)     = UNDEF
    dens_diff(:,:,:)     = UNDEF
    pott_diff(:,:,:)     = UNDEF
    qv_diff  (:,:,:)     = UNDEF
    num_diff (:,:,:,:,:) = UNDEF

    mflx_hi  (:,:,:,:) = UNDEF
    qflx_hi  (:,:,:,:) = UNDEF
    qflx_lo  (:,:,:,:) = UNDEF
    RHOQ     (:,:,:)   = UNDEF
#endif

!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
       MOMZ_RK1( 1:KS-1,i,j) = 0.0_RP
       MOMZ_RK1(KE:KA  ,i,j) = 0.0_RP
       MOMZ_RK2( 1:KS-1,i,j) = 0.0_RP
       MOMZ_RK2(KE:KA  ,i,j) = 0.0_RP
    enddo
    enddo

    ! bottom boundary
!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
       qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
       qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
       mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
    enddo
    enddo


    if ( USE_AVERAGE ) then

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMZ_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMX_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMY_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

#ifndef DRY
!OCL XFILL
    !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(4)
    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QTRC_av(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo
#endif

    endif

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,ZDIR) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,XDIR) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,YDIR) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       dens_s(k,i,j) = DENS(k,i,j)
    enddo
    enddo
    enddo

#ifndef DRY
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          QDRY(k,i,j) = 1.0_RP
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
       do iq = QQS, QQE
          QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          Rtot(k,i,j) = Rdry*QDRY(k,i,j) + Rvap*QTRC(k,i,j,I_QV)
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          CVtot(k,i,j) = CVdry*QDRY(k,i,j)
          do iq = QQS, QQE
             CVtot(k,i,j) = CVtot(k,i,j) + AQ_CV(iq)*QTRC(k,i,j,iq)
          enddo
       enddo
       enddo
       enddo

    enddo
    enddo
#endif

    if ( DIFF4 == 0.0_RP ) then
!OCL XFILL
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS  , JE
       do i = IS  , IE
       do k = KS  , KE
          num_diff(k,i,j,I_DENS,ZDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMZ,ZDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMX,ZDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMY,ZDIR) = 0.0_RP
          num_diff(k,i,j,I_RHOT,ZDIR) = 0.0_RP
       enddo
       enddo
       enddo
!OCL XFILL
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS  , JE
       do i = IS-1, IE
       do k = KS  , KE
          num_diff(k,i,j,I_DENS,XDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMZ,XDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMY,XDIR) = 0.0_RP
          num_diff(k,i,j,I_RHOT,XDIR) = 0.0_RP
       enddo
       enddo
       enddo
!OCL XFILL
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS  , JE
       do i = IS  , IE+1
       do k = KS  , KE
          num_diff(k,i,j,I_MOMX,XDIR) = 0.0_RP
       enddo
       enddo
       enddo
!OCL XFILL
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE
       do i = IS  , IE
       do k = KS  , KE
          num_diff(k,i,j,I_DENS,YDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMZ,YDIR) = 0.0_RP
          num_diff(k,i,j,I_MOMX,YDIR) = 0.0_RP
          num_diff(k,i,j,I_RHOT,YDIR) = 0.0_RP
       enddo
       enddo
       enddo
!OCL XFILL
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS  , JE+1
       do i = IS  , IE
       do k = KS  , KE
          num_diff(k,i,j,I_MOMY,YDIR) = 0.0_RP
       enddo
       enddo
       enddo
    endif



    do step = 1, NSTEP_ATMOS_DYN

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          DENS_t(k,i,j) = DENS_tp(k,i,j)
       enddo
       enddo
       enddo

       !--- prepare rayleigh damping coefficient
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS, KE-1
             MOMZ_t(k,i,j) = MOMZ_tp(k,i,j) &
                           - DAMP_alpha(k,i,j,I_BND_VELZ) * ( MOMZ(k,i,j) - DAMP_var(k,i,j,I_BND_VELZ) &
                                                          * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) ) ) ! rayleigh damping
          enddo
          MOMZ_t(KE,i,j) = 0.0_RP

          do k = KS, KE
             MOMX_t(k,i,j) = MOMX_tp(k,i,j) &
                           - DAMP_alpha(k,i,j,I_BND_VELX) * ( MOMX(k,i,j) - DAMP_var(k,i,j,I_BND_VELX) &
                                                          * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) ) ! rayleigh damping
             MOMY_t(k,i,j) = MOMY_tp(k,i,j) &
                           - DAMP_alpha(k,i,j,I_BND_VELY) * ( MOMY(k,i,j) - DAMP_var(k,i,j,I_BND_VELY) &
                                                          * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) ) ) ! rayleigh damping
             RHOT_t(k,i,j) = RHOT_tp(k,i,j) &
                           - DAMP_alpha(k,i,j,I_BND_POTT) * ( RHOT(k,i,j) - DAMP_var(k,i,j,I_BND_POTT) &
                                                          * DENS(k,i,j) ) ! rayleigh damping
          enddo
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( DENS_t(:,:,:), 1 )
    call COMM_vars8( MOMZ_t(:,:,:), 2 )
    call COMM_vars8( MOMX_t(:,:,:), 3 )
    call COMM_vars8( MOMY_t(:,:,:), 4 )
    call COMM_vars8( RHOT_t(:,:,:), 5 )
    call COMM_wait ( DENS_t(:,:,:), 1 )
    call COMM_wait ( MOMZ_t(:,:,:), 2 )
    call COMM_wait ( MOMX_t(:,:,:), 3 )
    call COMM_wait ( MOMY_t(:,:,:), 4 )
    call COMM_wait ( RHOT_t(:,:,:), 5 )

    !--- prepare numerical diffusion coefficient
    if ( DIFF4 /= 0.0_RP ) then

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j  = JJS-2, JJE+2
          do i  = IIS-2, IIE+2
          do k = KS, KE
             POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          enddo
          enddo
          enddo

          if ( ND_USE_RS ) then
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j  = JJS, JJE
             do i  = IIS, IIE
             do k = KS, KE
                dens_diff(k,i,j) = DENS(k,i,j)               - REF_dens(k)
                pott_diff(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) - REF_pott(k)
             enddo
             enddo
             enddo
          else
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j  = JJS, JJE
             do i  = IIS, IIE
             do k = KS+1, KE-1
                dens_diff(k,i,j) = ( DENS(k,i,j) * 3.0_RP &
                     + ( DENS(k,i+1,j) + DENS(k,i-1,j) + DENS(k,i,j+1) + DENS(k,i,j-1) ) * 2.0_RP &
                     + ( DENS(k,i+2,j) + DENS(k,i-2,j) + DENS(k,i,j+2) + DENS(k,i,j-2) ) &
                     + ( DENS(k+1,i,j) + DENS(k-1,i,j) ) * 2.0_RP &
                     ) / 19.0_RP
                pott_diff(k,i,j) = ( POTT(k,i,j) * 3.0_RP &
                     + ( POTT(k,i+1,j) + POTT(k,i-1,j) + POTT(k,i,j+1) + POTT(k,i,j-1) ) * 2.0_RP &
                     + ( POTT(k,i+2,j) + POTT(k,i-2,j) + POTT(k,i,j+2) + POTT(k,i,j-2) ) &
                     + ( POTT(k+1,i,j) + POTT(k-1,i,j) ) * 2.0_RP &
                     ) / 19.0_RP
             enddo
             enddo
             enddo
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j  = JJS, JJE
             do i  = IIS, IIE
                dens_diff(KS,i,j) = ( DENS(KS,i,j) * 3.0_RP &
                     + ( DENS(KS,i+1,j) + DENS(KS,i-1,j) + DENS(KS,i,j+1) + DENS(KS,i,j-1) ) * 2.0_RP &
                     + ( DENS(KS,i+2,j) + DENS(KS,i-2,j) + DENS(KS,i,j+2) + DENS(KS,i,j-2) ) &
                     + ( DENS(KS+1,i,j) ) * 2.0_RP &
                     ) / 17.0_RP
                dens_diff(KE,i,j) = ( DENS(KE,i,j) * 3.0_RP &
                     + ( DENS(KE,i+1,j) + DENS(KE,i-1,j) + DENS(KE,i,j+1) + DENS(KE,i,j-1) ) * 2.0_RP &
                     + ( DENS(KE,i+2,j) + DENS(KE,i-2,j) + DENS(KE,i,j+2) + DENS(KE,i,j-2) ) &
                     + ( DENS(KE-1,i,j) ) * 2.0_RP &
                     ) / 17.0_RP
                pott_diff(KS,i,j) = ( POTT(KS,i,j) * 3.0_RP &
                     + ( POTT(KS,i+1,j) + POTT(KS,i-1,j) + POTT(KS,i,j+1) + POTT(KS,i,j-1) ) * 2.0_RP &
                     + ( POTT(KS,i+2,j) + POTT(KS,i-2,j) + POTT(KS,i,j+2) + POTT(KS,i,j-2) ) &
                     + ( POTT(KS+1,i,j) ) * 2.0_RP &
                     ) / 17.0_RP
                pott_diff(KE,i,j) = ( POTT(KE,i,j) * 3.0_RP &
                     + ( POTT(KE,i+1,j) + POTT(KE,i-1,j) + POTT(KE,i,j+1) + POTT(KE,i,j-1) ) * 2.0_RP &
                     + ( POTT(KE,i+2,j) + POTT(KE,i-2,j) + POTT(KE,i,j+2) + POTT(KE,i,j-2) ) &
                     + ( POTT(KE-1,i,j) ) * 2.0_RP &
                     ) / 17.0_RP
             enddo
             enddo
          endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+2
          do i = IIS-2, IIE+2
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k,i,j) )
             call CHECK( __LINE__, DENS(k  ,i,j) )
             call CHECK( __LINE__, DENS(k+1,i,j) )
#endif
             VELZ(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+2
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k,i,j) )
             call CHECK( __LINE__, DENS(k,i  ,j) )
             call CHECK( __LINE__, DENS(k,i+1,j) )
#endif
             VELX(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-2, IIE+2
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i,j) )
             call CHECK( __LINE__, DENS(k,i,j  ) )
             call CHECK( __LINE__, DENS(k,i,j+1) )
#endif
             VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,k+1,1) )
          call CHECK( __LINE__, CNZ3(2,k+1,1) )
          call CHECK( __LINE__, CNZ3(3,k+1,1) )
          call CHECK( __LINE__, CNZ3(1,k,1) )
          call CHECK( __LINE__, dens_diff(k+2,i,j) )
          call CHECK( __LINE__, dens_diff(k+1,i,j) )
          call CHECK( __LINE__, dens_diff(k  ,i,j) )
          call CHECK( __LINE__, dens_diff(k-1,i,j) )
#endif
             num_diff(k,i,j,I_DENS,ZDIR) = &
                    ( CNZ3(1,k+1,1) * dens_diff(k+2,i,j) &
                    - CNZ3(2,k+1,1) * dens_diff(k+1,i,j) &
                    + CNZ3(3,k+1,1) * dens_diff(k  ,i,j) &
                    - CNZ3(1,k  ,1) * dens_diff(k-1,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS+1,1) )
          call CHECK( __LINE__, CNZ3(2,KS+1,1) )
          call CHECK( __LINE__, CNZ3(3,KS+1,1) )
          call CHECK( __LINE__, CNZ3(1,KS,1) )
          call CHECK( __LINE__, dens_diff(KS+2,i,j) )
          call CHECK( __LINE__, dens_diff(KS+1,i,j) )
          call CHECK( __LINE__, dens_diff(KS,i,j) )
#endif
             num_diff(KS  ,i,j,I_DENS,ZDIR) = &
                    ( CNZ3(1,KS+1,1) * dens_diff(KS+2,i,j) &
                    - CNZ3(2,KS+1,1) * dens_diff(KS+1,i,j) &
                    + CNZ3(3,KS+1,1) * dens_diff(KS  ,i,j) &
                    - CNZ3(1,KS  ,1) * dens_diff(KS+1,i,j) )
             num_diff(KS-1,i,j,I_DENS,ZDIR) = - num_diff(KS,i,j,I_DENS,ZDIR)
             num_diff(KS-2,i,j,I_DENS,ZDIR) = - num_diff(KS+1,i,j,I_DENS,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE,1) )
          call CHECK( __LINE__, CNZ3(2,KE,1) )
          call CHECK( __LINE__, CNZ3(3,KE,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, dens_diff(KE,i,j) )
          call CHECK( __LINE__, dens_diff(KE-1,i,j) )
          call CHECK( __LINE__, dens_diff(KE-2,i,j) )
#endif
             num_diff(KE-1,i,j,I_DENS,ZDIR) = &
                    ( CNZ3(1,KE  ,1) * dens_diff(KE-1,i,j) &
                    - CNZ3(2,KE  ,1) * dens_diff(KE  ,i,j) &
                    + CNZ3(3,KE  ,1) * dens_diff(KE-1,i,j) &
                    - CNZ3(1,KE-1,1) * dens_diff(KE-2,i,j) )
             num_diff(KE,i,j,I_DENS,ZDIR) = - num_diff(KE-1,i,j,I_DENS,ZDIR)
             num_diff(KE+1,i,j,I_DENS,ZDIR) = - num_diff(KE-2,i,j,I_DENS,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNX3(1,i+1,1) )
          call CHECK( __LINE__, CNX3(2,i+1,1) )
          call CHECK( __LINE__, CNX3(3,i+1,1) )
          call CHECK( __LINE__, CNX3(1,i,1) )
          call CHECK( __LINE__, dens(k,i+2,j) )
          call CHECK( __LINE__, dens(k,i+1,j) )
          call CHECK( __LINE__, dens(k,i  ,j) )
          call CHECK( __LINE__, dens(k,i-1,j) )
#endif
             num_diff(k,i,j,I_DENS,XDIR) = &
                    ( CNX3(1,i+1,1) * DENS(k,i+2,j) &
                    - CNX3(2,i+1,1) * DENS(k,i+1,j) &
                    + CNX3(3,i+1,1) * DENS(k,i  ,j) &
                    - CNX3(1,i  ,1) * DENS(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNY3(1,j+1,1) )
          call CHECK( __LINE__, CNY3(2,j+1,1) )
          call CHECK( __LINE__, CNY3(3,j+1,1) )
          call CHECK( __LINE__, CNY3(1,j,1) )
          call CHECK( __LINE__, dens(k,i,j+2) )
          call CHECK( __LINE__, dens(k,i,j+1) )
          call CHECK( __LINE__, dens(k,i,j  ) )
          call CHECK( __LINE__, dens(k,i,j-1) )
#endif
             num_diff(k,i,j,I_DENS,YDIR) = &
                    ( CNY3(1,j+1,1) * DENS(k,i,j+2) &
                    - CNY3(2,j+1,1) * DENS(k,i,j+1) &
                    + CNY3(3,j+1,1) * DENS(k,i,j  ) &
                    - CNY3(1,j  ,1) * DENS(k,i,j-1) )
          enddo
          enddo
          enddo

          ! z-momentum
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,k,2) )
          call CHECK( __LINE__, CNZ3(2,k,2) )
          call CHECK( __LINE__, CNZ3(3,k,2) )
          call CHECK( __LINE__, CNZ3(1,k-1,2) )
          call CHECK( __LINE__, VELZ(k+1,i,j) )
          call CHECK( __LINE__, VELZ(k  ,i,j) )
          call CHECK( __LINE__, VELZ(k-1,i,j) )
          call CHECK( __LINE__, VELZ(k-2,i,j) )
#endif
             num_diff(k,i,j,I_MOMZ,ZDIR) = &
                    ( CNZ3(1,k  ,2) * VELZ(k+1,i,j) &
                    - CNZ3(2,k  ,2) * VELZ(k  ,i,j) &
                    + CNZ3(3,k  ,2) * VELZ(k-1,i,j) &
                    - CNZ3(1,k-1,2) * VELZ(k-2,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS,2) )
          call CHECK( __LINE__, CNZ3(2,KS,2) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
#endif
             num_diff(KS  ,i,j,I_MOMZ,ZDIR) = &
                    ( CNZ3(1,KS  ,2) * VELZ(KS+1,i,j) &
                    - CNZ3(2,KS  ,2) * VELZ(KS  ,i,j) &
                    - CNZ3(1,KS-1,2) * VELZ(KS+1,i,j) )
             num_diff(KS-1,i,j,I_MOMZ,ZDIR) = - num_diff(KS  ,i,j,I_MOMZ,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS+1,2) )
          call CHECK( __LINE__, CNZ3(2,KS+1,2) )
          call CHECK( __LINE__, CNZ3(3,KS+1,2) )
          call CHECK( __LINE__, VELZ(KS+2,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
#endif
             num_diff(KS+1,i,j,I_MOMZ,ZDIR) = &
                    ( CNZ3(1,KS+1,2) * VELZ(KS+2,i,j) &
                    - CNZ3(2,KS+1,2) * VELZ(KS+1,i,j) &
                    + CNZ3(3,KS+1,2) * VELZ(KS  ,i,j) )
             num_diff(KS-2,i,j,I_MOMZ,ZDIR) = - num_diff(KS+1,i,j,I_MOMZ,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(2,KE-1,2) )
          call CHECK( __LINE__, CNZ3(3,KE-1,2) )
          call CHECK( __LINE__, CNZ3(1,KE-2,2) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-3,i,j) )
#endif
             num_diff(KE-1,i,j,I_MOMZ,ZDIR) = &
                    ( &
                    - CNZ3(2,KE-1,2) * VELZ(KE-1,i,j) &
                    + CNZ3(3,KE-1,2) * VELZ(KE-2,i,j) &
                    - CNZ3(1,KE-2,2) * VELZ(KE-3,i,j) )
             num_diff(KE+2,i,j,I_MOMZ,ZDIR) = - num_diff(KE-1,i,j,I_MOMZ,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(3,KE,2) )
          call CHECK( __LINE__, CNZ3(1,KE-1,2) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
#endif
             num_diff(KE  ,i,j,I_MOMZ,ZDIR) = &
                    ( CNZ3(1,KE  ,2) * VELZ(KE-1,i,j) &
                    + CNZ3(3,KE  ,2) * VELZ(KE-1,i,j) &
                    - CNZ3(1,KE-1,2) * VELZ(KE-2,i,j) )
             num_diff(KE+1,i,j,I_MOMZ,ZDIR) = - num_diff(KE,i,j,I_MOMZ,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, CNX3(1,i+1,1) )
          call CHECK( __LINE__, CNX3(2,i+1,1) )
          call CHECK( __LINE__, CNX3(3,i+1,1) )
          call CHECK( __LINE__, CNX3(1,i  ,1) )
          call CHECK( __LINE__, VELZ(k,i+2,j) )
          call CHECK( __LINE__, VELZ(k,i+1,j) )
          call CHECK( __LINE__, VELZ(k,i  ,j) )
          call CHECK( __LINE__, VELZ(k,i-1,j) )
#endif
             num_diff(k,i,j,I_MOMZ,XDIR) = &
                    ( CNX3(1,i+1,1) * VELZ(k,i+2,j) &
                    - CNX3(2,i+1,1) * VELZ(k,i+1,j) &
                    + CNX3(3,i+1,1) * VELZ(k,i  ,j) &
                    - CNX3(1,i  ,1) * VELZ(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, CNY3(1,j+1,1) )
          call CHECK( __LINE__, CNY3(2,j+1,1) )
          call CHECK( __LINE__, CNY3(3,j+1,1) )
          call CHECK( __LINE__, CNY3(1,j  ,1) )
          call CHECK( __LINE__, VELZ(k,i,j+2) )
          call CHECK( __LINE__, VELZ(k,i,j+1) )
          call CHECK( __LINE__, VELZ(k,i,j  ) )
          call CHECK( __LINE__, VELZ(k,i,j-1) )
#endif
             num_diff(k,i,j,I_MOMZ,YDIR) = &
                    ( CNY3(1,j+1,1) * VELZ(k,i,j+2) &
                    - CNY3(2,j+1,1) * VELZ(k,i,j+1) &
                    + CNY3(3,j+1,1) * VELZ(k,i,j  ) &
                    - CNY3(1,j  ,1) * VELZ(k,i,j-1) )
          enddo
          enddo
          enddo

          ! x-momentum
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS,  JJE
          do i = IIS,  IIE
          do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,k+1,1) )
          call CHECK( __LINE__, CNZ3(2,k+1,1) )
          call CHECK( __LINE__, CNZ3(3,k+1,1) )
          call CHECK( __LINE__, CNZ3(1,k  ,1) )
          call CHECK( __LINE__, VELX(k+2,i,j) )
          call CHECK( __LINE__, VELX(k+1,i,j) )
          call CHECK( __LINE__, VELX(k  ,i,j) )
          call CHECK( __LINE__, VELX(k-1,i,j) )
#endif
             num_diff(k,i,j,I_MOMX,ZDIR) = &
                    ( CNZ3(1,k+1,1) * VELX(k+2,i,j) &
                    - CNZ3(2,k+1,1) * VELX(k+1,i,j) &
                    + CNZ3(3,k+1,1) * VELX(k  ,i,j) &
                    - CNZ3(1,k  ,1) * VELX(k-1,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS+1,1) )
          call CHECK( __LINE__, CNZ3(2,KS+1,1) )
          call CHECK( __LINE__, CNZ3(3,KS+1,1) )
          call CHECK( __LINE__, CNZ3(1,KS  ,1) )
          call CHECK( __LINE__, VELX(KS+2,i,j) )
          call CHECK( __LINE__, VELX(KS+1,i,j) )
          call CHECK( __LINE__, VELX(KS  ,i,j) )
#endif
             num_diff(KS  ,i,j,I_MOMX,ZDIR) = &
                    ( CNZ3(1,KS+1,1) * VELX(KS+2,i,j) &
                    - CNZ3(2,KS+1,1) * VELX(KS+1,i,j) &
                    + CNZ3(3,KS+1,1) * VELX(KS  ,i,j) &
                    - CNZ3(1,KS  ,1) * VELX(KS+1,i,j) )
             num_diff(KS-1,i,j,I_MOMX,ZDIR) = - num_diff(KS  ,i,j,I_MOMX,ZDIR)
             num_diff(KS-2,i,j,I_MOMX,ZDIR) = - num_diff(KS+1,i,j,I_MOMX,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE  ,1) )
          call CHECK( __LINE__, CNZ3(2,KE  ,1) )
          call CHECK( __LINE__, CNZ3(3,KE  ,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, VELX(KE  ,i,j) )
          call CHECK( __LINE__, VELX(KE-1,i,j) )
          call CHECK( __LINE__, VELX(KE-2,i,j) )
#endif
             num_diff(KE-1,i,j,I_MOMX,ZDIR) = &
                    ( CNZ3(1,KE  ,1) * VELX(KE-1,i,j) &
                    - CNZ3(2,KE  ,1) * VELX(KE  ,i,j) &
                    + CNZ3(3,KE  ,1) * VELX(KE-1,i,j) &
                    - CNZ3(1,KE-1,1) * VELX(KE-2,i,j) )
             num_diff(KE  ,i,j,I_MOMX,ZDIR) = - num_diff(KE-1,i,j,I_MOMX,ZDIR)
             num_diff(KE+1,i,j,I_MOMX,ZDIR) = - num_diff(KE-2,i,j,I_MOMX,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNX3(1,i,2) )
          call CHECK( __LINE__, CNX3(2,i,2) )
          call CHECK( __LINE__, CNX3(3,i,2) )
          call CHECK( __LINE__, CNX3(1,i-1,2) )
          call CHECK( __LINE__, VELX(k,i+1,j) )
          call CHECK( __LINE__, VELX(k,i  ,j) )
          call CHECK( __LINE__, VELX(k,i-1,j) )
          call CHECK( __LINE__, VELX(k,i-2,j) )
#endif
             num_diff(k,i,j,I_MOMX,XDIR) = &
                    ( CNX3(1,i  ,2) * VELX(k,i+1,j) &
                    - CNX3(2,i  ,2) * VELX(k,i  ,j) &
                    + CNX3(3,i  ,2) * VELX(k,i-1,j) &
                    - CNX3(1,i-1,2) * VELX(k,i-2,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNY3(1,j+1,1) )
          call CHECK( __LINE__, CNY3(2,j+1,1) )
          call CHECK( __LINE__, CNY3(3,j+1,1) )
          call CHECK( __LINE__, CNY3(1,j  ,1) )
          call CHECK( __LINE__, VELX(k,i,j+2) )
          call CHECK( __LINE__, VELX(k,i,j+1) )
          call CHECK( __LINE__, VELX(k,i,j) )
          call CHECK( __LINE__, VELX(k,i,j-1) )
#endif
             num_diff(k,i,j,I_MOMX,YDIR) = &
                    ( CNY3(1,j+1,1) * VELX(k,i,j+2) &
                    - CNY3(2,j+1,1) * VELX(k,i,j+1) &
                    + CNY3(3,j+1,1) * VELX(k,i,j  ) &
                    - CNY3(1,j  ,1) * VELX(k,i,j-1) )
          enddo
          enddo
          enddo

          ! y-momentum
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,k+1,1) )
          call CHECK( __LINE__, CNZ3(2,k+1,1) )
          call CHECK( __LINE__, CNZ3(3,k+1,1) )
          call CHECK( __LINE__, CNZ3(1,k  ,1) )
          call CHECK( __LINE__, VELY(k+2,i,j) )
          call CHECK( __LINE__, VELY(k+1,i,j) )
          call CHECK( __LINE__, VELY(k  ,i,j) )
          call CHECK( __LINE__, VELY(k-1,i,j) )
#endif
             num_diff(k,i,j,I_MOMY,ZDIR) = &
                    ( CNZ3(1,k+1,1) * VELY(k+2,i,j) &
                    - CNZ3(2,k+1,1) * VELY(k+1,i,j) &
                    + CNZ3(3,k+1,1) * VELY(k  ,i,j) &
                    - CNZ3(1,k  ,1) * VELY(k-1,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS+1,1) )
          call CHECK( __LINE__, CNZ3(2,KS+1,1) )
          call CHECK( __LINE__, CNZ3(3,KS+1,1) )
          call CHECK( __LINE__, CNZ3(1,KS  ,1) )
          call CHECK( __LINE__, VELY(KS+2,i,j) )
          call CHECK( __LINE__, VELY(KS+1,i,j) )
          call CHECK( __LINE__, VELY(KS,i,j) )
#endif
             num_diff(KS  ,i,j,I_MOMY,ZDIR) = &
                    ( CNZ3(1,KS+1,1) * VELY(KS+2,i,j) &
                    - CNZ3(2,KS+1,1) * VELY(KS+1,i,j) &
                    + CNZ3(3,KS+1,1) * VELY(KS  ,i,j) &
                    - CNZ3(1,KS  ,1) * VELY(KS+1,i,j) )
             num_diff(KS-1,i,j,I_MOMY,ZDIR) = - num_diff(KS  ,i,j,I_MOMY,ZDIR)
             num_diff(KS-2,i,j,I_MOMY,ZDIR) = - num_diff(KS+1,i,j,I_MOMY,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE,1) )
          call CHECK( __LINE__, CNZ3(2,KE,1) )
          call CHECK( __LINE__, CNZ3(3,KE,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, VELY(KE,i,j) )
          call CHECK( __LINE__, VELY(KE-1,i,j) )
          call CHECK( __LINE__, VELY(KE-2,i,j) )
#endif
             num_diff(KE-1,i,j,I_MOMY,ZDIR) = &
                    ( CNZ3(1,KE  ,1) * VELY(KE-1,i,j) &
                    - CNZ3(2,KE  ,1) * VELY(KE  ,i,j) &
                    + CNZ3(3,KE  ,1) * VELY(KE-1,i,j) &
                    - CNZ3(1,KE-1,1) * VELY(KE-2,i,j) )
             num_diff(KE  ,i,j,I_MOMY,ZDIR) = - num_diff(KE-1,i,j,I_MOMY,ZDIR)
             num_diff(KE+1,i,j,I_MOMY,ZDIR) = - num_diff(KE-2,i,j,I_MOMY,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNX3(1,i+1,1) )
          call CHECK( __LINE__, CNX3(2,i+1,1) )
          call CHECK( __LINE__, CNX3(3,i+1,1) )
          call CHECK( __LINE__, CNX3(1,i  ,1) )
          call CHECK( __LINE__, VELY(k,i+2,j) )
          call CHECK( __LINE__, VELY(k,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j) )
          call CHECK( __LINE__, VELY(k,i-1,j) )
#endif
             num_diff(k,i,j,I_MOMY,XDIR) = &
                    ( CNX3(1,i+1,1) * VELY(k,i+2,j) &
                    - CNX3(2,i+1,1) * VELY(k,i+1,j) &
                    + CNX3(3,i+1,1) * VELY(k,i  ,j) &
                    - CNX3(1,i  ,1) * VELY(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNY3(1,j  ,2) )
          call CHECK( __LINE__, CNY3(2,j  ,2) )
          call CHECK( __LINE__, CNY3(3,j  ,2) )
          call CHECK( __LINE__, CNY3(1,j-1,2) )
          call CHECK( __LINE__, VELY(k,i,j+1) )
          call CHECK( __LINE__, VELY(k,i,j  ) )
          call CHECK( __LINE__, VELY(k,i,j-1) )
          call CHECK( __LINE__, VELY(k,i,j-2) )
#endif
             num_diff(k,i,j,I_MOMY,YDIR) = &
                    ( CNY3(1,j  ,2) * VELY(k,i,j+1) &
                    - CNY3(2,j  ,2) * VELY(k,i,j  ) &
                    + CNY3(3,j  ,2) * VELY(k,i,j-1) &
                    - CNY3(1,j-1,2) * VELY(k,i,j-2) )
          enddo
          enddo
          enddo

          ! rho * theta
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,k+1,1) )
          call CHECK( __LINE__, CNZ3(2,k+1,1) )
          call CHECK( __LINE__, CNZ3(3,k+1,1) )
          call CHECK( __LINE__, CNZ3(1,k  ,1) )
          call CHECK( __LINE__, pott_diff(k+2,i,j) )
          call CHECK( __LINE__, pott_diff(k+1,i,j) )
          call CHECK( __LINE__, pott_diff(k  ,i,j) )
          call CHECK( __LINE__, pott_diff(k-1,i,j) )
#endif
             num_diff(k,i,j,I_RHOT,ZDIR) = &
                    ( CNZ3(1,k+1,1) * pott_diff(k+2,i,j)   &
                    - CNZ3(2,k+1,1) * pott_diff(k+1,i,j)   &
                    + CNZ3(3,k+1,1) * pott_diff(k  ,i,j)   &
                    - CNZ3(1,k  ,1) * pott_diff(k-1,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS+1,1) )
          call CHECK( __LINE__, CNZ3(2,KS+1,1) )
          call CHECK( __LINE__, CNZ3(3,KS+1,1) )
          call CHECK( __LINE__, CNZ3(1,KS,1) )
          call CHECK( __LINE__, pott_diff(KS+2,i,j) )
          call CHECK( __LINE__, pott_diff(KS+1,i,j) )
          call CHECK( __LINE__, pott_diff(KS,i,j) )
#endif
             num_diff(KS  ,i,j,I_RHOT,ZDIR) = &
                    ( CNZ3(1,KS+1,1) * pott_diff(KS+2,i,j)   &
                    - CNZ3(2,KS+1,1) * pott_diff(KS+1,i,j)   &
                    + CNZ3(3,KS+1,1) * pott_diff(KS  ,i,j)   &
                    - CNZ3(1,KS  ,1) * pott_diff(KS+1,i,j) )
             num_diff(KS-1,i,j,I_RHOT,ZDIR) = - num_diff(KS  ,i,j,I_RHOT,ZDIR)
             num_diff(KS-2,i,j,I_RHOT,ZDIR) = - num_diff(KS+1,i,j,I_RHOT,ZDIR)
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE,1) )
          call CHECK( __LINE__, CNZ3(2,KE,1) )
          call CHECK( __LINE__, CNZ3(3,KE,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, pott_diff(KE,i,j) )
          call CHECK( __LINE__, pott_diff(KE-1,i,j) )
          call CHECK( __LINE__, pott_diff(KE-2,i,j) )
#endif
             num_diff(KE-1,i,j,I_RHOT,ZDIR) = &
                    ( CNZ3(1,KE  ,1) * pott_diff(KE  ,i,j)   &
                    - CNZ3(2,KE  ,1) * pott_diff(KE  ,i,j)   &
                    + CNZ3(3,KE  ,1) * pott_diff(KE-1,i,j)   &
                    - CNZ3(1,KE-1,1) * pott_diff(KE-2,i,j) )
             num_diff(KE  ,i,j,I_RHOT,ZDIR) = - num_diff(KE-1,i,j,I_RHOT,ZDIR)
             num_diff(KE+1,i,j,I_RHOT,ZDIR) = - num_diff(KE-2,i,j,I_RHOT,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNX3(1,i+1,1) )
          call CHECK( __LINE__, CNX3(2,i+1,1) )
          call CHECK( __LINE__, CNX3(3,i+1,1) )
          call CHECK( __LINE__, CNX3(1,i  ,1) )
          call CHECK( __LINE__, pott(k,i+2,j) )
          call CHECK( __LINE__, pott(k,i+1,j) )
          call CHECK( __LINE__, pott(k,i  ,j) )
          call CHECK( __LINE__, pott(k,i-1,j) )
#endif
             num_diff(k,i,j,I_RHOT,XDIR) = &
                    ( CNX3(1,i+1,1) * POTT(k,i+2,j)   &
                    - CNX3(2,i+1,1) * POTT(k,i+1,j)   &
                    + CNX3(3,i+1,1) * POTT(k,i  ,j)   &
                    - CNX3(1,i  ,1) * POTT(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNY3(1,j+1,1) )
          call CHECK( __LINE__, CNY3(2,j+1,1) )
          call CHECK( __LINE__, CNY3(3,j+1,1) )
          call CHECK( __LINE__, CNY3(1,j  ,1) )
          call CHECK( __LINE__, pott(k,i,j+2) )
          call CHECK( __LINE__, pott(k,i,j+1) )
          call CHECK( __LINE__, pott(k,i,j  ) )
          call CHECK( __LINE__, pott(k,i,j-1) )
#endif
             num_diff(k,i,j,I_RHOT,YDIR) = &
                    ( CNY3(1,j+1,1) * POTT(k,i,j+2)   &
                    - CNY3(2,j+1,1) * POTT(k,i,j+1)   &
                    + CNY3(3,j+1,1) * POTT(k,i,j  )   &
                    - CNY3(1,j  ,1) * POTT(k,i,j-1) )
          enddo
          enddo
          enddo


       enddo
       enddo ! end tile


       num_diff_pt0 => num_diff
       num_diff_pt1 => num_diff_work

       do no = 2, nd_order

          call COMM_vars8( num_diff_pt0(:,:,:,I_DENS,ZDIR),  1 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_DENS,XDIR),  2 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_DENS,YDIR),  3 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMZ,ZDIR),  4 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMZ,XDIR),  5 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMZ,YDIR),  6 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMX,ZDIR),  7 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMX,XDIR),  8 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMX,YDIR),  9 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMY,ZDIR), 10 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMY,XDIR), 11 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_MOMY,YDIR), 12 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_RHOT,ZDIR), 13 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_RHOT,XDIR), 14 )
          call COMM_vars8( num_diff_pt0(:,:,:,I_RHOT,YDIR), 15 )

          call COMM_wait( num_diff_pt0(:,:,:,I_DENS,ZDIR), 1 )
          call COMM_wait( num_diff_pt0(:,:,:,I_DENS,XDIR), 2 )
          call COMM_wait( num_diff_pt0(:,:,:,I_DENS,YDIR), 3 )
          call calc_numdiff4(     &
               num_diff_pt1,      & ! (out)
               num_diff_pt0,      & ! (in)
               CNZ4(:,:,1),       & ! (in)
               CNX4(:,:,1),       & ! (in)
               CNY4(:,:,1),       & ! (in)
               I_DENS,            & ! (in)
               KE                 ) ! (in)

          call COMM_wait( num_diff_pt0(:,:,:,I_MOMZ,ZDIR), 4 )
          call COMM_wait( num_diff_pt0(:,:,:,I_MOMZ,XDIR), 5 )
          call COMM_wait( num_diff_pt0(:,:,:,I_MOMZ,YDIR), 6 )
          call calc_numdiff4(     &
               num_diff_pt1,      & ! (out)
               num_diff_pt0,      & ! (in)
               CNZ4(:,:,2),       & ! (in)
               CNX4(:,:,1),       & ! (in)
               CNY4(:,:,1),       & ! (in)
               I_MOMZ,            & ! (in)
               KE-1               ) ! (in)

          call COMM_wait( num_diff_pt0(:,:,:,I_MOMX,ZDIR), 7 )
          call COMM_wait( num_diff_pt0(:,:,:,I_MOMX,XDIR), 8 )
          call COMM_wait( num_diff_pt0(:,:,:,I_MOMX,YDIR), 9 )
          call calc_numdiff4(     &
               num_diff_pt1,      & ! (out)
               num_diff_pt0,      & ! (in)
               CNZ4(:,:,1),       & ! (in)
               CNX4(:,:,2),       & ! (in)
               CNY4(:,:,1),       & ! (in)
               I_MOMX,            & ! (in)
               KE                 ) ! (in)

          call COMM_wait( num_diff_pt0(:,:,:,I_MOMY,ZDIR), 10 )
          call COMM_wait( num_diff_pt0(:,:,:,I_MOMY,XDIR), 11 )
          call COMM_wait( num_diff_pt0(:,:,:,I_MOMY,YDIR), 12 )
          call calc_numdiff4(     &
               num_diff_pt1,      & ! (out)
               num_diff_pt0,      & ! (in)
               CNZ4(:,:,1),       & ! (in)
               CNX4(:,:,1),       & ! (in)
               CNY4(:,:,2),       & ! (in)
               I_MOMY,            & ! (in)
               KE                 ) ! (in)

          call COMM_wait( num_diff_pt0(:,:,:,I_RHOT,ZDIR), 13 )
          call COMM_wait( num_diff_pt0(:,:,:,I_RHOT,XDIR), 14 )
          call COMM_wait( num_diff_pt0(:,:,:,I_RHOT,YDIR), 15 )
          call calc_numdiff4(     &
               num_diff_pt1,      & ! (out)
               num_diff_pt0,      & ! (in)
               CNZ4(:,:,1),       & ! (in)
               CNX4(:,:,1),       & ! (in)
               CNY4(:,:,1),       & ! (in)
               I_RHOT,            & ! (in)
               KE                 ) ! (in)

          ! swap pointer target
          tmp_pt       => num_diff_pt1
          num_diff_pt1 => num_diff_pt0
          num_diff_pt0 => tmp_pt

       enddo

       nd_order4 = nd_order * 4

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS-1, KE
             num_diff(k,i,j,I_DENS,ZDIR) = num_diff_pt0(k,i,j,I_DENS,ZDIR) &
                  * DIFF4 * CDZ(k)**nd_order4
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_DENS,XDIR) = num_diff_pt0(k,i,j,I_DENS,XDIR) &
                  * DIFF4 * CDX(i)**nd_order4
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_DENS,XDIR) = num_DIFF(KS,i,j,I_DENS,XDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_DENS,XDIR) = num_DIFF(KS+1,i,j,I_DENS,XDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_DENS,YDIR) = num_diff_pt0(k,i,j,I_DENS,YDIR) &
                  * DIFF4 * CDY(j)**nd_order4
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_DENS,YDIR) = num_DIFF(KS,i,j,I_DENS,YDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_DENS,YDIR) = num_DIFF(KS+1,i,j,I_DENS,YDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          ! z-momentum
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1,  KE-1
             num_diff(k,i,j,I_MOMZ,ZDIR) = num_diff_pt0(k,i,j,I_MOMZ,ZDIR) &
                  * DIFF4 * ( 0.5_RP*(CDZ(k+1)+CDZ(k)) )**nd_order4 &
                  * DENS(k,i,j)
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             num_diff(k,i,j,I_MOMZ,XDIR) = num_diff_pt0(k,i,j,I_MOMZ,XDIR) &
                  * DIFF4 * CDX(i)**nd_order4 &
                  * 0.25_RP * ( DENS(k+1,i+1,j) + DENS(k+1,i,j) + DENS(k,i+1,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_MOMZ,XDIR) = num_DIFF(KS,i,j,I_MOMZ,XDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_MOMZ,XDIR) = num_DIFF(KS+1,i,j,I_MOMZ,XDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             num_diff(k,i,j,I_MOMZ,YDIR) = num_diff_pt0(k,i,j,I_MOMZ,YDIR) &
                  * DIFF4 * CDY(j)**nd_order4 &
                  * 0.25_RP * ( DENS(k+1,i,j+1) + DENS(k+1,i,j) + DENS(k,i,j+1) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_MOMZ,YDIR) = num_DIFF(KS,i,j,I_MOMZ,YDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_MOMZ,YDIR) = num_DIFF(KS+1,i,j,I_MOMZ,YDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          ! x-momentum
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS,    KE-1
             num_diff(k,i,j,I_MOMX,ZDIR) = num_diff_pt0(k,i,j,I_MOMX,ZDIR) &
                  * DIFF4 * CDZ(k)**nd_order4 &
                  * 0.25_RP * ( DENS(k+1,i+1,j) + DENS(k+1,i,j) + DENS(k,i+1,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMX,XDIR) = num_diff_pt0(k,i,j,I_MOMX,XDIR) &
                  * DIFF4 * ( 0.5_RP*(CDX(i+1)+CDX(i)) )**nd_order4 &
                  * DENS(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_MOMX,XDIR) = num_DIFF(KS,i,j,I_MOMX,XDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_MOMX,XDIR) = num_DIFF(KS+1,i,j,I_MOMX,XDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMX,YDIR) = num_diff_pt0(k,i,j,I_MOMX,YDIR) &
                  * DIFF4 * CDY(j)**nd_order4 &
                  * 0.25_RP * ( DENS(k,i+1,j+1) + DENS(k,i+1,j) + DENS(k,i,j+1) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_MOMX,YDIR) = num_DIFF(KS,i,j,I_MOMX,YDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_MOMX,YDIR) = num_DIFF(KS+1,i,j,I_MOMX,YDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          ! y-momentum
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS,  KE-1
             num_diff(k,i,j,I_MOMY,ZDIR) = num_diff_pt0(k,i,j,I_MOMY,ZDIR) &
                  * DIFF4 * CDZ(k)**nd_order4 &
                  * 0.25_RP * ( DENS(k+1,i,j+1) + DENS(k+1,i,j) + DENS(k,i,j+1) + DENS(k,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMY,XDIR) = num_diff_pt0(k,i,j,I_MOMY,XDIR) &
                  * DIFF4 * CDX(i)**nd_order4 &
                  * 0.25_RP * ( DENS(k,i+1,j+1) + DENS(k,i,j+1) + DENS(k,i+1,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_MOMY,XDIR) = num_DIFF(KS,i,j,I_MOMY,XDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_MOMY,XDIR) = num_DIFF(KS+1,i,j,I_MOMY,XDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMY,YDIR) = num_diff_pt0(k,i,j,I_MOMY,YDIR) &
                  * DIFF4 * ( 0.5_RP*(CDY(j+1)+CDY(j)) )**nd_order4 &
                  * DENS(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_MOMY,YDIR) = num_DIFF(KS,i,j,I_MOMY,YDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_MOMY,YDIR) = num_DIFF(KS+1,i,j,I_MOMY,YDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          ! rho * theta
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS  , KE-1
             num_diff(k,i,j,I_RHOT,ZDIR) = num_diff_pt0(k,i,j,I_RHOT,ZDIR) &
                  * DIFF4 * CDZ(k)**nd_order4 &
                  * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS-1,i,j,I_RHOT,ZDIR) = num_diff_pt0(KS-1,i,j,I_RHOT,ZDIR) &
                  * DIFF4 * CDZ(KS-1)**nd_order4 &
                  * DENS(KS,i,j)
             num_diff(KE  ,i,j,I_RHOT,ZDIR) = num_diff_pt0(KE  ,i,j,I_RHOT,ZDIR) &
                  * DIFF4 * CDZ(KE  )**nd_order4 &
                  * DENS(KE,i,j)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_RHOT,XDIR) = num_diff_pt0(k,i,j,I_RHOT,XDIR) &
                  * DIFF4 * CDX(i)**nd_order4 &
                  * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_RHOT,XDIR) = num_DIFF(KS,i,j,I_RHOT,XDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_RHOT,XDIR) = num_DIFF(KS+1,i,j,I_RHOT,XDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_RHOT,YDIR) = num_diff_pt0(k,i,j,I_RHOT,YDIR) &
                  * DIFF4 * CDY(j)**nd_order4 &
                  * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             num_diff(KS,i,j,I_RHOT,YDIR) = num_DIFF(KS,i,j,I_RHOT,YDIR) &
                  * ND_SFC_FACT
             num_diff(KS+1,i,j,I_RHOT,YDIR) = num_DIFF(KS+1,i,j,I_RHOT,YDIR) &
                  * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          enddo
          enddo

       enddo
       enddo

       call COMM_vars8( num_diff(:,:,:,I_DENS,ZDIR),  1 )
       call COMM_vars8( num_diff(:,:,:,I_DENS,XDIR),  2 )
       call COMM_vars8( num_diff(:,:,:,I_DENS,YDIR),  3 )
       call COMM_vars8( num_diff(:,:,:,I_MOMZ,ZDIR),  4 )
       call COMM_vars8( num_diff(:,:,:,I_MOMZ,XDIR),  5 )
       call COMM_vars8( num_diff(:,:,:,I_MOMZ,YDIR),  6 )
       call COMM_vars8( num_diff(:,:,:,I_MOMX,ZDIR),  7 )
       call COMM_vars8( num_diff(:,:,:,I_MOMX,XDIR),  8 )
       call COMM_vars8( num_diff(:,:,:,I_MOMX,YDIR),  9 )
       call COMM_vars8( num_diff(:,:,:,I_MOMY,ZDIR), 10 )
       call COMM_vars8( num_diff(:,:,:,I_MOMY,XDIR), 11 )
       call COMM_vars8( num_diff(:,:,:,I_MOMY,YDIR), 12 )
       call COMM_vars8( num_diff(:,:,:,I_RHOT,ZDIR), 13 )
       call COMM_vars8( num_diff(:,:,:,I_RHOT,XDIR), 14 )
       call COMM_vars8( num_diff(:,:,:,I_RHOT,YDIR), 15 )

       call COMM_wait( num_diff(:,:,:,I_DENS,ZDIR),  1 )
       call COMM_wait( num_diff(:,:,:,I_DENS,XDIR),  2 )
       call COMM_wait( num_diff(:,:,:,I_DENS,YDIR),  3 )
       call COMM_wait( num_diff(:,:,:,I_MOMZ,ZDIR),  4 )
       call COMM_wait( num_diff(:,:,:,I_MOMZ,XDIR),  5 )
       call COMM_wait( num_diff(:,:,:,I_MOMZ,YDIR),  6 )
       call COMM_wait( num_diff(:,:,:,I_MOMX,ZDIR),  7 )
       call COMM_wait( num_diff(:,:,:,I_MOMX,XDIR),  8 )
       call COMM_wait( num_diff(:,:,:,I_MOMX,YDIR),  9 )
       call COMM_wait( num_diff(:,:,:,I_MOMY,ZDIR), 10 )
       call COMM_wait( num_diff(:,:,:,I_MOMY,XDIR), 11 )
       call COMM_wait( num_diff(:,:,:,I_MOMY,YDIR), 12 )
       call COMM_wait( num_diff(:,:,:,I_RHOT,ZDIR), 13 )
       call COMM_wait( num_diff(:,:,:,I_RHOT,XDIR), 14 )
       call COMM_wait( num_diff(:,:,:,I_RHOT,YDIR), 15 )

    endif ! DIFF4 /= 0.0_RP

    !##### Start RK #####

    !##### RK1 #####
    rko  = 1
    dtrk = DTSEC_ATMOS_DYN / 3.D0

    call ATMOS_DYN_rk( DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (inout)
                       mflx_hi,                                          & ! (out)
                       DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (inout)
                       DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (in)
                       DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                       Rtot, CVtot, CORIOLI,                             & ! (in)
                       num_diff, divdmp_coef,                            & ! (in)
                       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                       CDZ, FDZ, FDX, FDY,                               & ! (in)
                       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                       dtrk, RK, rko,                                    & ! (in)
                       VELZ, VELX, VELY, PRES, POTT                      ) ! (work)

    call COMM_vars8( DENS_RK1(:,:,:), 1 )
    call COMM_vars8( MOMZ_RK1(:,:,:), 2 )
    call COMM_vars8( MOMX_RK1(:,:,:), 3 )
    call COMM_vars8( MOMY_RK1(:,:,:), 4 )
    call COMM_vars8( RHOT_RK1(:,:,:), 5 )
    call COMM_wait ( DENS_RK1(:,:,:), 1 )
    call COMM_wait ( MOMZ_RK1(:,:,:), 2 )
    call COMM_wait ( MOMX_RK1(:,:,:), 3 )
    call COMM_wait ( MOMY_RK1(:,:,:), 4 )
    call COMM_wait ( RHOT_RK1(:,:,:), 5 )

    !##### RK2 #####
    rko  = 2
    dtrk = DTSEC_ATMOS_DYN / 2.D0

    call ATMOS_DYN_rk( DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (inout)
                       mflx_hi,                                          & ! (out)
                       DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (inout)
                       DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1, & ! (in)
                       DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                       Rtot, CVtot, CORIOLI,                             & ! (in)
                       num_diff, divdmp_coef,                            & ! (in)
                       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                       CDZ, FDZ, FDX, FDY,                               & ! (in)
                       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                       dtrk, RK, rko,                                    & ! (in)
                       VELZ, VELX, VELY, PRES, POTT                      ) ! (work)

    call COMM_vars8( DENS_RK2(:,:,:), 1 )
    call COMM_vars8( MOMZ_RK2(:,:,:), 2 )
    call COMM_vars8( MOMX_RK2(:,:,:), 3 )
    call COMM_vars8( MOMY_RK2(:,:,:), 4 )
    call COMM_vars8( RHOT_RK2(:,:,:), 5 )
    call COMM_wait ( DENS_RK2(:,:,:), 1 )
    call COMM_wait ( MOMZ_RK2(:,:,:), 2 )
    call COMM_wait ( MOMX_RK2(:,:,:), 3 )
    call COMM_wait ( MOMY_RK2(:,:,:), 4 )
    call COMM_wait ( RHOT_RK2(:,:,:), 5 )

    !##### RK3 #####
    rko  = 3
    dtrk = DTSEC_ATMOS_DYN

    call ATMOS_DYN_rk( DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (inout)
                       mflx_hi,                                          & ! (out)
                       DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! (inout)
                       DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2, & ! (in)
                       DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
                       Rtot, CVtot, CORIOLI,                             & ! (in)
                       num_diff, divdmp_coef,                            & ! (in)
                       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,      & ! (in)
                       CDZ, FDZ, FDX, FDY,                               & ! (in)
                       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
                       dtrk, RK, rko,                                    & ! (in)
                       VELZ, VELX, VELY, PRES, POTT                      ) ! (work)

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

    if ( USE_AVERAGE ) then

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          DENS_av(k,i,j) = DENS_av(k,i,j) + DENS(k,i,j)
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          MOMZ_av(k,i,j) = MOMZ_av(k,i,j) + MOMZ(k,i,j)
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          MOMX_av(k,i,j) = MOMX_av(k,i,j) + MOMX(k,i,j)
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          MOMY_av(k,i,j) = MOMY_av(k,i,j) + MOMY(k,i,j)
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          RHOT_av(k,i,j) = RHOT_av(k,i,j) + RHOT(k,i,j)
       enddo
       enddo
       enddo

    endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          mflx_av(k,i,j,ZDIR) = mflx_av(k,i,j,ZDIR) + mflx_hi(k,i,j,ZDIR)
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          mflx_av(k,i,j,XDIR) = mflx_av(k,i,j,XDIR) + mflx_hi(k,i,j,XDIR)
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          mflx_av(k,i,j,YDIR) = mflx_av(k,i,j,YDIR) + mflx_hi(k,i,j,YDIR)
       enddo
       enddo
       enddo

    enddo
    enddo

    enddo ! dynamical steps

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       mflx_hi(k,i,j,ZDIR) = mflx_av(k,i,j,ZDIR) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       mflx_hi(k,i,j,XDIR) = mflx_av(k,i,j,XDIR) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       mflx_hi(k,i,j,YDIR) = mflx_av(k,i,j,YDIR) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    call COMM_vars8( mflx_hi(:,:,:,ZDIR), 1 )
    call COMM_vars8( mflx_hi(:,:,:,XDIR), 2 )
    call COMM_vars8( mflx_hi(:,:,:,YDIR), 3 )
    call COMM_wait ( mflx_hi(:,:,:,ZDIR), 1 )
    call COMM_wait ( mflx_hi(:,:,:,XDIR), 2 )
    call COMM_wait ( mflx_hi(:,:,:,YDIR), 3 )

#ifndef DRY
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          QTRC_t(k,i,j,I_QV) = QTRC_tp(k,i,j,I_QV) &
               - DAMP_alpha(k,i,j,I_BND_QV) &
               * ( QTRC(k,i,j,I_QV) - DAMP_var(k,i,j,I_BND_QV) ) * DENS(k,i,j)
       enddo
       enddo
       enddo

       !$omp parallel do private(iq,i,j,k) OMP_SCHEDULE_ collapse(3)
       do iq = I_QV+1, QA
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          QTRC_t(k,i,j,iq) = QTRC_tp(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo

    enddo
    enddo

    do iq = 1, QA
       call COMM_vars8( QTRC_t(:,:,:,iq), 5+iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC_t(:,:,:,iq), 5+iq )
    enddo

    do iq = 1, QA

       if ( DIFF4 /= 0.0_RP ) then

          do JJS = JS, JE, JBLOCK
          JJE = JJS+JBLOCK-1
          do IIS = IS, IE, IBLOCK
          IIE = IIS+IBLOCK-1

          if ( iq == I_QV ) then
             if ( ND_USE_RS ) then
                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
                do j  = JJS, JJE
                do i  = IIS, IIE
                do k = KS, KE
                   qv_diff(k,i,j) = QTRC(k,i,j,I_QV) - REF_qv  (k)
                enddo
                enddo
                enddo
             else
                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
                do j  = JJS, JJE
                do i  = IIS, IIE
                do k = KS+1, KE-1
                   qv_diff(k,i,j) = ( QTRC(k,i,j,I_QV) * 3.0_RP &
                        + ( QTRC(k,i+1,j,I_QV) + QTRC(k,i-1,j,I_QV) + QTRC(k,i,j+1,I_QV) + QTRC(k,i,j-1,I_QV) ) * 2.0_RP &
                        +   QTRC(k,i+2,j,I_QV) + QTRC(k,i-2,j,I_QV) + QTRC(k,i,j+2,I_QV) + QTRC(k,i,j-2,I_QV) &
                        + ( QTRC(k+1,i,j,I_QV) + QTRC(k-1,i,j,I_QV) ) * 2.0_RP &
                        ) / 19.0_RP
                enddo
                enddo
                enddo
                !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
                do j  = JJS, JJE
                do i  = IIS, IIE
                   qv_diff(KS,i,j) = ( QTRC(KS,i,j,I_QV) * 3.0_RP &
                        + ( QTRC(KS,i+1,j,I_QV) + QTRC(KS,i-1,j,I_QV) + QTRC(KS,i,j+1,I_QV) + QTRC(KS,i,j-1,I_QV) ) * 2.0_RP &
                        +   QTRC(KS,i+2,j,I_QV) + QTRC(KS,i-2,j,I_QV) + QTRC(KS,i,j+2,I_QV) + QTRC(KS,i,j-2,I_QV) &
                        +   QTRC(KS+1,i,j,I_QV) * 2.0_RP &
                   ) / 17.0_RP
                   qv_diff(KE,i,j) = ( QTRC(KE,i,j,I_QV) * 3.0_RP &
                        + ( QTRC(KE,i+1,j,I_QV) + QTRC(KE,i-1,j,I_QV) + QTRC(KE,i,j+1,I_QV) + QTRC(KE,i,j-1,I_QV) ) * 2.0_RP &
                        +   QTRC(KE,i+2,j,I_QV) + QTRC(KE,i-2,j,I_QV) + QTRC(KE,i,j+2,I_QV) + QTRC(KE,i,j-2,I_QV) &
                        +   QTRC(KE-1,i,j,I_QV) * 2.0_RP &
                   ) / 17.0_RP
                enddo
                enddo
             endif
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS,   JJE
             do i = IIS,   IIE
             do k = KS+1, KE-2
                num_diff(k,i,j,1,ZDIR) = &
                       ( CNZ3(1,k+1,1) * qv_diff(k+2,i,j)   &
                       - CNZ3(2,k+1,1) * qv_diff(k+1,i,j)   &
                       + CNZ3(3,k+1,1) * qv_diff(k  ,i,j)   &
                       - CNZ3(1,k  ,1) * qv_diff(k-1,i,j) )
             enddo
             enddo
             enddo

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS,   JJE
             do i = IIS,   IIE
                num_diff(KS  ,i,j,1,ZDIR) = &
                       ( CNZ3(1,KS+1,1) * qv_diff(KS+2,i,j)   &
                       - CNZ3(2,KS+1,1) * qv_diff(KS+1,i,j)   &
                       + CNZ3(3,KS+1,1) * qv_diff(KS  ,i,j)   &
                       - CNZ3(1,KS  ,1) * qv_diff(KS+1,i,j) )
                num_diff(KS-1,i,j,1,ZDIR) = - num_diff(KS  ,i,j,1,ZDIR)
                num_diff(KS-2,i,j,1,ZDIR) = - num_diff(KS+1,i,j,1,ZDIR)
                num_diff(KE-1,i,j,1,ZDIR) = &
                       ( CNZ3(1,KE  ,1) * qv_diff(KE  ,i,j)   &
                       - CNZ3(2,KE  ,1) * qv_diff(KE  ,i,j)   &
                       + CNZ3(3,KE  ,1) * qv_diff(KE-1,i,j)   &
                       - CNZ3(1,KE-1,1) * qv_diff(KE-2,i,j) )
                num_diff(KE  ,i,j,1,ZDIR) = - num_diff(KE-1,i,j,1,ZDIR)
                num_diff(KE+1,i,j,1,ZDIR) = - num_diff(KE-2,i,j,1,ZDIR)
             enddo
             enddo

          else

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS,   JJE
             do i = IIS,   IIE
             do k = KS+1, KE-2
                num_diff(k,i,j,1,ZDIR) = &
                       ( CNZ3(1,k+1,1) * QTRC(k+2,i,j,iq)   &
                       - CNZ3(2,k+1,1) * QTRC(k+1,i,j,iq)   &
                       + CNZ3(3,k+1,1) * QTRC(k  ,i,j,iq)   &
                       - CNZ3(1,k  ,1) * QTRC(k-1,i,j,iq) )
             enddo
             enddo
             enddo

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS,   JJE
             do i = IIS,   IIE
#ifdef DEBUG
                call CHECK( __LINE__, CNZ3(1,KS+1,1) )
                call CHECK( __LINE__, CNZ3(2,KS+1,1) )
                call CHECK( __LINE__, CNZ3(3,KS+1,1) )
                call CHECK( __LINE__, CNZ3(1,KS,1) )
                call CHECK( __LINE__, QTRC(KS+2,i,j,iq) )
                call CHECK( __LINE__, QTRC(KS+1,i,j,iq) )
                call CHECK( __LINE__, QTRC(KS,i,j,iq) )
#endif
                num_diff(KS  ,i,j,1,ZDIR) = &
                       ( CNZ3(1,KS+1,1) * QTRC(KS+2,i,j,iq)   &
                       - CNZ3(2,KS+1,1) * QTRC(KS+1,i,j,iq)   &
                       + CNZ3(3,KS+1,1) * QTRC(KS  ,i,j,iq)   &
                       - CNZ3(1,KS  ,1) * QTRC(KS+1,i,j,iq) )
                num_diff(KS-1,i,j,1,ZDIR) = - num_diff(KS  ,i,j,1,ZDIR)
                num_diff(KS-2,i,j,1,ZDIR) = - num_diff(KS+1,i,j,1,ZDIR)
#ifdef DEBUG
                call CHECK( __LINE__, CNZ3(1,KE,1) )
                call CHECK( __LINE__, CNZ3(2,KE,1) )
                call CHECK( __LINE__, CNZ3(3,KE,1) )
                call CHECK( __LINE__, CNZ3(1,KE-1,1) )
                call CHECK( __LINE__, QTRC(KE,i,j,iq) )
                call CHECK( __LINE__, QTRC(KE-1,i,j,iq) )
                call CHECK( __LINE__, QTRC(KE-2,i,j,iq) )
#endif
                num_diff(KE-1,i,j,1,ZDIR) = &
                       ( CNZ3(1,KE  ,1) * QTRC(KE  ,i,j,iq)   &
                       - CNZ3(2,KE  ,1) * QTRC(KE  ,i,j,iq)   &
                       + CNZ3(3,KE  ,1) * QTRC(KE-1,i,j,iq)   &
                       - CNZ3(1,KE-1,1) * QTRC(KE-2,i,j,iq) )
                num_diff(KE  ,i,j,1,ZDIR) = - num_diff(KE-1,i,j,1,ZDIR)
                num_diff(KE+1,i,j,1,ZDIR) = - num_diff(KE-2,i,j,1,ZDIR)
             enddo
             enddo

          endif ! iq == I_QV

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE
             num_diff(k,i,j,1,XDIR) = &
                  ( CNX3(1,i+1,1) * QTRC(k,i+2,j,iq)   &
                  - CNX3(2,i+1,1) * QTRC(k,i+1,j,iq)   &
                  + CNX3(3,i+1,1) * QTRC(k,i  ,j,iq)   &
                  - CNX3(1,i  ,1) * QTRC(k,i-1,j,iq) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE
             num_diff(k,i,j,1,YDIR) = &
                  ( CNY3(1,j+1,1) * QTRC(k,i,j+2,iq)   &
                  - CNY3(2,j+1,1) * QTRC(k,i,j+1,iq)   &
                  + CNY3(3,j+1,1) * QTRC(k,i,j  ,iq)   &
                  - CNY3(1,j  ,1) * QTRC(k,i,j-1,iq) )
          enddo
          enddo
          enddo

          enddo
          enddo


          num_diff_pt0 => num_diff
          num_diff_pt1 => num_diff_work

          do no = 2, nd_order

             call COMM_vars8( num_diff_pt0(:,:,:,1,ZDIR), 1 )
             call COMM_vars8( num_diff_pt0(:,:,:,1,XDIR), 2 )
             call COMM_vars8( num_diff_pt0(:,:,:,1,YDIR), 3 )

             call COMM_wait( num_diff_pt0(:,:,:,1,ZDIR), 1 )
             call COMM_wait( num_diff_pt0(:,:,:,1,XDIR), 2 )
             call COMM_wait( num_diff_pt0(:,:,:,1,YDIR), 3 )

             call calc_numdiff4(     &
                  num_diff_pt1,      & ! (out)
                  num_diff_pt0,      & ! (in)
                  CNZ4(:,:,1),       & ! (in)
                  CNX4(:,:,1),       & ! (in)
                  CNY4(:,:,1),       & ! (in)
                  1,                 & ! (in)
                  KE                 ) ! (in)

             ! swap pointer target
             tmp_pt => num_diff_pt1
             num_diff_pt1 => num_diff_pt0
             num_diff_pt0 => tmp_pt

          enddo ! no

          do JJS = JS, JE, JBLOCK
          JJE = JJS+JBLOCK-1
          do IIS = IS, IE, IBLOCK
          IIE = IIS+IBLOCK-1

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS  , KE-1
                num_diff(k,i,j,1,ZDIR) = num_diff_pt0(k,i,j,1,ZDIR) &
                     * DIFF4 * CDZ(k)**nd_order4 &
                     * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
             enddo
             enddo
             enddo

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
                num_diff(KS-1,i,j,1,ZDIR) = num_diff_pt0(KS-1,i,j,1,ZDIR) &
                     * DIFF4 * CDZ(KS-1)**nd_order4 &
                     * DENS(KS,i,j)
                num_diff(KE  ,i,j,1,ZDIR) = num_diff_pt0(KE  ,i,j,1,ZDIR) &
                     * DIFF4 * CDZ(KE  )**nd_order4 &
                     * DENS(KE,i,j)
             enddo
             enddo

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                num_diff(k,i,j,1,XDIR) = num_diff_pt0(k,i,j,1,XDIR) &
                     * DIFF4 * CDX(i)**nd_order4 &
                     * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
             enddo
             enddo
             enddo

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                num_diff(k,i,j,1,YDIR) = num_diff_pt0(k,i,j,1,YDIR) &
                     * DIFF4 * CDY(j)**nd_order4 &
                     * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
             enddo
             enddo
             enddo

          enddo
          enddo

          call COMM_vars8( num_diff(:,:,:,1,ZDIR), 1 )
          call COMM_vars8( num_diff(:,:,:,1,XDIR), 2 )
          call COMM_vars8( num_diff(:,:,:,1,YDIR), 3 )

          call COMM_wait( num_diff(:,:,:,1,ZDIR), 1 )
          call COMM_wait( num_diff(:,:,:,1,XDIR), 2 )
          call COMM_wait( num_diff(:,:,:,1,YDIR), 3 )



       endif ! DIFF4 .eq. 0.0


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-2
          qflx_lo(k,i,j,ZDIR) = 0.5_RP * (     mflx_hi(k,i,j,ZDIR)  * ( QTRC(k+1,i,j,iq)+QTRC(k,i,j,iq) ) &
                                         - abs(mflx_hi(k,i,j,ZDIR)) * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) )

          qflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( QTRC(k+1,i,j,iq)+QTRC(k  ,i,j,iq) ) &
                                + FACT_F * ( QTRC(k+2,i,j,iq)+QTRC(k-1,i,j,iq) ) ) &
                              + num_diff(k,i,j,1,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          ! just above the bottom boundary
          qflx_lo(KS  ,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KS  ,i,j,ZDIR)  * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                            - abs(mflx_hi(KS  ,i,j,ZDIR)) * ( QTRC(KS+1,i,j,iq)-QTRC(KS,i,j,iq) ) )

          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * mflx_hi(KS  ,i,j,ZDIR) * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                 + num_diff(KS,i,j,1,ZDIR)
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          ! just below the top boundary
          qflx_lo(KE-1,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KE-1,i,j,ZDIR)  * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                            - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) )
          qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP

          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR) * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                 + num_diff(KE-1,i,j,1,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-2, IIE+1
       do k = KS, KE
          qflx_lo(k,i,j,XDIR) = 0.5_RP * (     mflx_hi(k,i,j,XDIR)  * ( QTRC(k,i+1,j,iq)+QTRC(k,i,j,iq) ) &
                                         - abs(mflx_hi(k,i,j,XDIR)) * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( QTRC(k,i+1,j,iq)+QTRC(k,i  ,j,iq) ) &
                                + FACT_F * ( QTRC(k,i+2,j,iq)+QTRC(k,i-1,j,iq) ) ) &
                              + num_diff(k,i,j,1,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          qflx_lo(k,i,j,YDIR) = 0.5_RP * (     mflx_hi(k,i,j,YDIR)  * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j,iq) ) &
                                         - abs(mflx_hi(k,i,j,YDIR)) * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) &
                              * ( FACT_N * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j  ,iq) ) &
                                + FACT_F * ( QTRC(k,i,j+2,iq)+QTRC(k,i,j-1,iq) ) ) &
                              + num_diff(k,i,j,1,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          RHOQ(k,i,j) = QTRC(k,i,j,iq) * dens_s(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    call ATMOS_DYN_fct(qflx_anti,               & ! (out)
                       RHOQ, qflx_hi, qflx_lo,  & ! (in)
                       RCDZ, RCDX, RCDY, Real(DTSEC,kind=RP)  ) ! (in)
#ifdef DEBUG
       qflx_lo  (KS:,:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          QTRC(k,i,j,iq) = ( RHOQ(k,i,j) &
              + DTSEC * ( - ( ( qflx_hi(k  ,i  ,j  ,ZDIR) + qflx_anti(k  ,i  ,j  ,ZDIR) &
                              - qflx_hi(k-1,i  ,j  ,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( qflx_hi(k  ,i  ,j  ,XDIR) + qflx_anti(k  ,i  ,j  ,XDIR) &
                              - qflx_hi(k  ,i-1,j  ,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k  ,i  ,j  ,YDIR) + qflx_anti(k  ,i  ,j  ,YDIR) &
                              - qflx_hi(k  ,i  ,j-1,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                          + QTRC_t(k,i,j,iq) &
                        ) ) / DENS(k,i,j)
       enddo
       enddo
       enddo

    enddo
    enddo
#ifdef DEBUG
    qflx_hi  (KS:,:,:,:) = UNDEF
    qflx_anti(:,:,:,:) = UNDEF
#endif

    enddo ! scalar quantities loop

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo

#endif

!OCL XFILL
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j  = JS, JE
    do i  = IS, IE
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       MOMZ(   1:KS-1,i,j) = MOMZ(KS,i,j)
       MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
       MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
       RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
       MOMZ(KE+1:KA,  i,j) = MOMZ(KE,i,j)
       MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
       MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
       RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
    enddo
    enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
#ifndef DRY
!OCL XFILL
    !$omp parallel do private(iq,i,j) OMP_SCHEDULE_ collapse(3)
    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
       QTRC(   1:KS-1,i,j,iq) = QTRC(KS,i,j,iq)
       QTRC(KE+1:KA,  i,j,iq) = QTRC(KE,i,j,iq)
    enddo
    enddo
    enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF; iq = IUNDEF
#endif
#endif

    if ( USE_AVERAGE ) then

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS_av(k,i,j) = DENS_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMZ_av(k,i,j) = MOMZ_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMX_av(k,i,j) = MOMX_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMY_av(k,i,j) = MOMY_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT_av(k,i,j) = RHOT_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

#ifndef DRY
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(4)
    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QTRC_av(k,i,j,iq) = QTRC(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo
#endif

    endif

    return
  end subroutine ATMOS_DYN_main

  !-----------------------------------------------------------------------------
  subroutine calc_numdiff4( &
       num_diff_pt1, & ! (out)
       num_diff_pt0, & ! (in)
       CNZ4,         & ! (in)
       CNX4,         & ! (in)
       CNY4,         & ! (in)
       I_val,        & ! (in)
       k1            ) ! (in)
    implicit none

    real(RP), intent(out) :: num_diff_pt1(KA,IA,JA,5,3)
    real(RP), intent(in)  :: num_diff_pt0(KA,IA,JA,5,3)
    real(RP), intent(in)  :: CNZ4(5,KA)
    real(RP), intent(in)  :: CNX4(5,IA)
    real(RP), intent(in)  :: CNY4(5,JA)
    integer,  intent(in)  :: I_val
    integer,  intent(in)  :: k1

    integer :: i, j, k
    integer :: iis, iie, jjs, jje
    !---------------------------------------------------------------------------

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1


    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, CNZ4(1,k) )
       call CHECK( __LINE__, CNZ4(2,k) )
       call CHECK( __LINE__, CNZ4(3,k) )
       call CHECK( __LINE__, CNZ4(4,k) )
       call CHECK( __LINE__, CNZ4(5,k) )
       call CHECK( __LINE__, num_diff_pt0(k+2,i,j,I_val,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k+1,i,j,I_val,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k  ,i,j,I_val,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k-1,i,j,I_val,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k-2,i,j,I_val,ZDIR) )
#endif
       num_diff_pt1(k,i,j,I_val,ZDIR) = &
                     ( CNZ4(1,k) * num_diff_pt0(k+2,i,j,I_val,ZDIR) &
                     - CNZ4(2,k) * num_diff_pt0(k+1,i,j,I_val,ZDIR) &
                     + CNZ4(3,k) * num_diff_pt0(k  ,i,j,I_val,ZDIR) &
                     - CNZ4(4,k) * num_diff_pt0(k-1,i,j,I_val,ZDIR) &
                     + CNZ4(5,k) * num_diff_pt0(k-2,i,j,I_val,ZDIR) )
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
       num_diff_pt1(KS-1,i,j,I_val,ZDIR) = - num_diff_pt1(KS  ,i,j,I_val,ZDIR)
       num_diff_pt1(KS-2,i,j,I_val,ZDIR) = - num_diff_pt1(KS+1,i,j,I_val,ZDIR)
       num_diff_pt1(KE  ,i,j,I_val,ZDIR) = - num_diff_pt1(KE-1,i,j,I_val,ZDIR)
       num_diff_pt1(KE+1,i,j,I_val,ZDIR) = - num_diff_pt1(KE-2,i,j,I_val,ZDIR)
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS, K1
#ifdef DEBUG
       call CHECK( __LINE__, CNX4(1,i) )
       call CHECK( __LINE__, CNX4(2,i) )
       call CHECK( __LINE__, CNX4(3,i) )
       call CHECK( __LINE__, CNX4(4,i) )
       call CHECK( __LINE__, CNX4(5,i) )
       call CHECK( __LINE__, num_diff_pt0(k,i-2,j,I_val,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i+1,j,I_val,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i  ,j,I_val,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i-1,j,I_val,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i-2,j,I_val,XDIR) )
#endif
       num_diff_pt1(k,i,j,I_val,XDIR) = &
                    ( CNX4(1,i) * num_diff_pt0(k,i+2,j,I_val,XDIR) &
                    - CNX4(2,i) * num_diff_pt0(k,i+1,j,I_val,XDIR) &
                    + CNX4(3,i) * num_diff_pt0(k,i  ,j,I_val,XDIR) &
                    - CNX4(4,i) * num_diff_pt0(k,i-1,j,I_val,XDIR) &
                    + CNX4(5,i) * num_diff_pt0(k,i-2,j,I_val,XDIR) )
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS, K1
#ifdef DEBUG
       call CHECK( __LINE__, CNY4(1,j) )
       call CHECK( __LINE__, CNY4(2,j) )
       call CHECK( __LINE__, CNY4(3,j) )
       call CHECK( __LINE__, CNY4(4,j) )
       call CHECK( __LINE__, CNY4(5,j) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j-2,I_val,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j+1,I_val,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j  ,I_val,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j-1,I_val,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j-2,I_val,YDIR) )
#endif
       num_diff_pt1(k,i,j,I_val,YDIR) = &
                    ( CNY4(1,j) * num_diff_pt0(k,i,j+2,I_val,YDIR) &
                    - CNY4(2,j) * num_diff_pt0(k,i,j+1,I_val,YDIR) &
                    + CNY4(3,j) * num_diff_pt0(k,i,j  ,I_val,YDIR) &
                    - CNY4(4,j) * num_diff_pt0(k,i,j-1,I_val,YDIR) &
                    + CNY4(5,j) * num_diff_pt0(k,i,j-2,I_val,YDIR) )
    enddo
    enddo
    enddo


    enddo
    enddo

    return
  end subroutine calc_numdiff4

end module mod_atmos_dyn
