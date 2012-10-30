!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FENT + FCT
!!
!! @par Description
!!          Dynamical core for Atmospheric process
!!          Full explicit, no terrain + tracer FCT limiter
!!
!! @author H.Tomita and SCALE developpers
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
!!
!<
!-------------------------------------------------------------------------------
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
#ifdef DEBUG
  real(RP), private, parameter :: UNDEF = -9.999E30_RP
  integer , private, parameter :: IUNDEF = -99999
#endif
  integer, private, parameter :: I_DENS = 1
  integer, private, parameter :: I_MOMZ = 2
  integer, private, parameter :: I_MOMX = 3
  integer, private, parameter :: I_MOMY = 4
  integer, private, parameter :: I_RHOT = 5

  integer, private, parameter :: ZDIR = 1
  integer, private, parameter :: XDIR = 2
  integer, private, parameter :: YDIR = 3

  ! time settings
  integer, private, parameter :: RK = 3 ! order of Runge-Kutta scheme

  ! advection settings
  real(RP), private, parameter :: FACT_N =   7.0_RP / 6.0_RP !  7/6: fourth, 1: second
  real(RP), private, parameter :: FACT_F = - 1.0_RP / 6.0_RP ! -1/6: fourth, 0: second
  real(RP), parameter :: EPSILON = 1.0E-30_RP

  ! numerical filter settings
  real(RP), private, save      :: ATMOS_DYN_numerical_diff = 1.0E-2_RP ! nondimensional numerical diffusion
  real(RP), private, save      :: DIFF4 ! for 4th order numerical filter

  ! coriolis force
  logical, private, save       :: ATMOS_DYN_enable_coriolis = .false. ! enable coriolis force?
  real(RP), private, save      :: CORIOLI(1,IA,JA)                    ! coriolis term

  ! large scale sinking current (option)
  real(RP), private, save      :: ATMOS_DYN_LSsink_D = 0.0_RP           ! Divergence parameter [1/s]
  real(RP), private, save      :: ATMOS_DYN_divdmp_coef = 0.0_RP        ! Divergence dumping coef
  ! fct
  logical, private, save       :: ATMOS_DYN_FLAG_FCT_rho      = .false.
  logical, private, save       :: ATMOS_DYN_FLAG_FCT_momentum = .false.
  logical, private, save       :: ATMOS_DYN_FLAG_FCT_T        = .false.

  ! work
  real(RP), private, save :: DENS_RK1(KA,IA,JA)   ! prognostic variables (+1/3 step)
  real(RP), private, save :: MOMZ_RK1(KA,IA,JA)   !
  real(RP), private, save :: MOMX_RK1(KA,IA,JA)   !
  real(RP), private, save :: MOMY_RK1(KA,IA,JA)   !
  real(RP), private, save :: RHOT_RK1(KA,IA,JA)   !
  real(RP), private, save :: DENS_RK2(KA,IA,JA)   ! prognostic variables (+2/3 step)
  real(RP), private, save :: MOMZ_RK2(KA,IA,JA)   !
  real(RP), private, save :: MOMX_RK2(KA,IA,JA)   !
  real(RP), private, save :: MOMY_RK2(KA,IA,JA)   !
  real(RP), private, save :: RHOT_RK2(KA,IA,JA)   !

  real(RP), private, save :: mflx_hi(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face high order
  real(RP), private, save :: mflx_av(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face average
  real(RP), private, save :: rjpls  (KA,IA,JA)    ! correction factor for incoming antidiffusive flux
  real(RP), private, save :: rjmns  (KA,IA,JA)    ! correction factor for outgoing antidiffusive flux

  real(RP), private, save :: CNDZ(3,KA)
  real(RP), private, save :: CNMZ(3,KA)
  real(RP), private, save :: CNDX(3,IA)
  real(RP), private, save :: CNMX(3,IA)
  real(RP), private, save :: CNDY(3,JA)
  real(RP), private, save :: CNMY(3,JA)

  !-----------------------------------------------------------------------------
contains

#ifdef DEBUG
  subroutine CHECK( line, v )
    integer,  intent(in) :: line
    real(RP), intent(in) :: v
    if ( .not. v .gt. UNDEF ) then
       write(*,*) "use uninitialized value at line ", line
       call abort
    end if
  end subroutine CHECK
#endif

  !-----------------------------------------------------------------------------
  !> Initialize Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only : &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY
    use mod_geometrics, only : &
       lat => GEOMETRICS_lat
    use mod_atmos_vars, only: &
       ATMOS_TYPE_DYN
#ifdef _USE_RDMA
    use mod_comm, only: &
       COMM_set_rdma_variable
#endif
    implicit none

    NAMELIST / PARAM_ATMOS_DYN /    &
       ATMOS_DYN_numerical_diff,    &
       ATMOS_DYN_enable_coriolis,   &
       ATMOS_DYN_divdmp_coef,       &
       ATMOS_DYN_LSsink_D,          &
       ATMOS_DYN_FLAG_FCT_rho,      &
       ATMOS_DYN_FLAG_FCT_momentum, &
       ATMOS_DYN_FLAG_FCT_T

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Dynamics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** FENT + FCT'

    if ( trim(ATMOS_TYPE_DYN) .ne. 'FENT-FCT' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_DYN is not FENT-FCT. Check!'
       call PRC_MPIstop
    end if

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

    call ATMOS_DYN_init( DIFF4, CORIOLI,                      & ! (out)
                         CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY,  & ! (out)
                         CDZ, CDX, CDY,                       & ! (in)
                         lat,                                 & ! (in)
                         ATMOS_DYN_numerical_diff,            & ! (in)
                         ATMOS_DYN_enable_coriolis            ) ! (in)

#ifdef _USE_RDMA
    ! RDMA setting
    call COMM_set_rdma_variable( DENS_RK1(:,:,:), 5+QA+ 1)
    call COMM_set_rdma_variable( MOMZ_RK1(:,:,:), 5+QA+ 2)
    call COMM_set_rdma_variable( MOMY_RK1(:,:,:), 5+QA+ 3)
    call COMM_set_rdma_variable( MOMX_RK1(:,:,:), 5+QA+ 4)
    call COMM_set_rdma_variable( RHOT_RK1(:,:,:), 5+QA+ 5)

    call COMM_set_rdma_variable( DENS_RK2(:,:,:), 5+QA+ 6)
    call COMM_set_rdma_variable( MOMZ_RK2(:,:,:), 5+QA+ 7)
    call COMM_set_rdma_variable( MOMY_RK2(:,:,:), 5+QA+ 8)
    call COMM_set_rdma_variable( MOMX_RK2(:,:,:), 5+QA+ 9)
    call COMM_set_rdma_variable( RHOT_RK2(:,:,:), 5+QA+10)

    call COMM_set_rdma_variable( mflx_hi(:,:,:,ZDIR), 5+QA+11)
    call COMM_set_rdma_variable( mflx_hi(:,:,:,XDIR), 5+QA+12)
    call COMM_set_rdma_variable( mflx_hi(:,:,:,YDIR), 5+QA+13)

    call COMM_set_rdma_variable( rjpls  (:,:,:),      5+QA+14)
    call COMM_set_rdma_variable( rjmns  (:,:,:),      5+QA+15)
#endif

    return

  end subroutine ATMOS_DYN_setup

  subroutine ATMOS_DYN_init( DIFF4, corioli,                     &
                             CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY, &
                             CDZ, CDX, CDY, lat,                 &
                             numerical_diff,                     &
                             enable_coriolis                     )
    use mod_const, only: &
       PI   => CONST_PI, &
       EOHM => CONST_EOHM
    implicit none
    real(RP), intent(out) :: DIFF4
    real(RP), intent(out) :: corioli(1,IA,JA)
    real(RP), intent(out) :: CNDZ(3,KA)
    real(RP), intent(out) :: CNMZ(3,KA)
    real(RP), intent(out) :: CNDX(3,IA)
    real(RP), intent(out) :: CNMX(3,IA)
    real(RP), intent(out) :: CNDY(3,JA)
    real(RP), intent(out) :: CNMY(3,JA)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)
    real(RP), intent(in)  :: lat(1,IA,JA) 
    real(RP), intent(in)  :: numerical_diff
    logical , intent(in)  :: enable_coriolis

    real(RP) :: d2r
    integer :: i, j, k

#ifdef DEBUG
    CNDZ(:,:) = UNDEF
    CNMZ(:,:) = UNDEF
    CNDX(:,:) = UNDEF
    CNMX(:,:) = UNDEF
    CNDY(:,:) = UNDEF
    CNMY(:,:) = UNDEF
    corioli(:,:,:) = UNDEF
#endif

    ! coriolis parameter
    if ( enable_coriolis ) then
       d2r = PI / 180.0_RP
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = 1, JA
       do i = 1, IA
          corioli(1,i,j) = 2.0_RP * EOHM * sin( lat(1,i,j) * d2r )
       enddo
       enddo
    else
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = 1, JA
       do i = 1, IA
          corioli(1,i,j) = 0.0_RP
       enddo
       enddo
    endif


    ! numerical diffusion
    if ( numerical_diff == 0.0_RP ) then
       DIFF4 = 0.0_RP
       CNDZ(:,:) = 0.0_RP
       CNMZ(:,:) = 0.0_RP
       CNDX(:,:) = 0.0_RP
       CNMX(:,:) = 0.0_RP
       CNDY(:,:) = 0.0_RP
       CNMY(:,:) = 0.0_RP
    else
       DIFF4 = - numerical_diff * (-1.0_RP)**( 4/2+1 )

       ! z djrectjon
       do k = KS-1, KE+1
          CNDZ(1,k) = 1.0_RP / ( (CDZ(k+1)+CDZ(k  )) * 0.5_RP * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP )
       enddo
       CNDZ(1,   1:KS-2) = CNDZ(1,KS-1)
       CNDZ(1,KE+2:KA  ) = CNDZ(1,KE+1)

       do k = KS-1, KE+1
          CNDZ(2,k) = 1.0_RP / ( (CDZ(k+1)+CDZ(k  )) * 0.5_RP * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP )
       enddo
       CNDZ(2,   1:KS-2) = CNDZ(2,KS-1)
       CNDZ(2,KE+2:KA  ) = CNDZ(2,KE+1)

       do k = KS, KE+2
          CNDZ(3,k) = 1.0_RP / ( (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k-1) * (CDZ(k-1)+CDZ(k-2)) * 0.5_RP )
       enddo
       CNDZ(3,1   :KS-1) = CNDZ(3,KS  )
       CNDZ(3,KE+2:KA  ) = CNDZ(3,KE+2)

       do k = KS-2, KE+1
          CNMZ(1,k) = 1.0_RP / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5_RP * CDZ(k  ) )
       enddo
       CNMZ(1,   1:KS-2) = CNMZ(1,KS-2)
       CNMZ(1,KE+2:KA  ) = CNMZ(1,KE+1)

       do k = KS-1, KE+1
          CNMZ(2,k) = 1.0_RP / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5_RP * CDZ(k  ) ) &   
                    + 1.0_RP / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5_RP * CDZ(k  ) ) &  
                    + 1.0_RP / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k  ) ) 
       enddo
       CNMZ(2,   1:KS-2) = CNMZ(2,KS-1)
       CNMZ(2,KE+2:KA  ) = CNMZ(2,KE+1)

       do k = KS-1, KE+1
          CNMZ(3,k) = 1.0_RP / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5_RP * CDZ(k  ) ) &
                    + 1.0_RP / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k  ) ) &
                    + 1.0_RP / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5_RP * CDZ(k-1) )
       enddo
       CNMZ(3,   1:KS-2) = CNMZ(3,KS-1)
       CNMZ(3,KE+2:KA  ) = CNMZ(3,KE+1)

       ! x direction
       do i = IS-1, IE+1
          CNDX(1,i) = 1.0_RP / ( (CDX(i+1)+CDX(i  )) * 0.5_RP * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP )
       enddo
       CNDX(1,   1:IS-2) = CNDX(1,IS-1)
       CNDX(1,IE+2:IA  ) = CNDX(1,IE+1)

       do i = IS-1, IE+1
          CNDX(2,i) = 1.0_RP / ( (CDX(i+1)+CDX(i  )) * 0.5_RP * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5_RP )
       enddo
       CNDX(2,   1:IS-2) = CNDX(2,IS-1)
       CNDX(2,IE+2:IA  ) = CNDX(2,IE+1)

       do i = IS, IE+2
          CNDX(3,i) = 1.0_RP / ( (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i-1) * (CDX(i-1)+CDX(i-2)) * 0.5_RP )
       enddo
       CNDX(3,   1:IS-1) = CNDX(3,IS  )
       CNDX(3,IE+2:IA  ) = CNDX(3,IE+2)

       do i = IS-2, IE+1
          CNMX(1,i) = 1.0_RP / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5_RP * CDX(i  ) )
       enddo
       CNMX(1,IE+2:IA  ) = CNMX(1,IE+1)
       CNMX(1,   1:IS-2) = CNMX(1,IS-2)

       do i = IS-1, IE+1
          CNMX(2,i) = 1.0_RP / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5_RP * CDX(i  ) ) & 
                    + 1.0_RP / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5_RP * CDX(i  ) ) & 
                    + 1.0_RP / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i  ) )
       enddo
       CNMX(2,   1:IS-2) = CNMX(2,IS-1)
       CNMX(2,IE+2:IA  ) = CNMX(2,IE+1)

       do i = IS-1, IE+1
          CNMX(3,i) = 1.0_RP / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5_RP * CDX(i  ) ) & 
                    + 1.0_RP / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i  ) ) & 
                    + 1.0_RP / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5_RP * CDX(i-1) )
       enddo
       CNMX(3,   1:IS-2) = CNMX(3,IS-1)
       CNMX(3,IE+2:IA  ) = CNMX(3,IE+1)

       ! y direction
       do j = JS-1, JE+1
          CNDY(1,j) = 1.0_RP / ( (CDY(j+1)+CDY(j  )) * 0.5_RP * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP )
       enddo
       CNDY(1,   1:JS-2) = CNDY(1,JS-1)
       CNDY(1,JE+2:JA  ) = CNDY(1,JE+1)

       do j = JS-1, JE+1
          CNDY(2,j) = 1.0_RP / ( (CDY(j+1)+CDY(j  )) * 0.5_RP * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5_RP )
       enddo
       CNDY(2,   1:JS-2) = CNDY(2,JS-1)
       CNDY(2,JE+2:JA  ) = CNDY(2,JE+1)

       do j = JS, JE+2
          CNDY(3,j) = 1.0_RP / ( (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5_RP ) &
                    + 1.0_RP / ( (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j-1) * (CDY(j-1)+CDY(j-2)) * 0.5_RP )
       enddo
       CNDY(3,   1:JS-1) = CNDY(3,JS  )
       CNDY(3,JE+2:JA  ) = CNDY(3,JE+2)

       do j = JS-2, JE+1
          CNMY(1,j) = 1.0_RP / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5_RP * CDY(j  ) )
       enddo
       CNMY(1,   1:JS-2) = CNMY(1,JS-2)
       CNMY(1,JE+2:JA  ) = CNMY(1,JE+1)

       do j = JS-1, JE+1
          CNMY(2,j) = 1.0_RP / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5_RP * CDY(j  ) ) &   
                    + 1.0_RP / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5_RP * CDY(j  ) ) &  
                    + 1.0_RP / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j  ) ) 
       enddo
       CNMY(2,   1:JS-2) = CNMY(2,JS-1)
       CNMY(2,JE+2:JA  ) = CNMY(2,JE+1)

       do j = JS-1, JE+1
          CNMY(3,j) = 1.0_RP / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5_RP * CDY(j  ) ) &
                    + 1.0_RP / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j  ) ) &
                    + 1.0_RP / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5_RP * CDY(j-1) )
       enddo
       CNMY(3,   1:JS-2) = CNMY(3,JS-1)
       CNMY(3,JE+2:JA  ) = CNMY(3,JE+1)
    end if

    return
  end subroutine ATMOS_DYN_init

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  !-----------------------------------------------------------------------------
!OCL NORECURRENCE
  subroutine ATMOS_DYN
    use mod_time, only: &
       DTSEC           => TIME_DTSEC,           &
       DTSEC_ATMOS_DYN => TIME_DTSEC_ATMOS_DYN, &
       NSTEP_ATMOS_DYN => TIME_NSTEP_ATMOS_DYN
    use mod_grid, only : &
       CZ   => GRID_CZ,   &
       FZ   => GRID_FZ,   &
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
    use mod_comm, only: &
#ifdef _USE_RDMA
       COMM_rdma_vars8, &
#endif
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
       ATMOS_vars_total,  &
       ATMOS_USE_AVERAGE, &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       DENS_av, &
       MOMZ_av, &
       MOMX_av, &
       MOMY_av, &
       RHOT_av, &
       QTRC_av, &
       qflx_sgs_momz, &
       qflx_sgs_momx, &
       qflx_sgs_momy, &
       qflx_sgs_rhot, &
       qflx_sgs_qtrc
    use mod_atmos_refstate, only: &
       REF_dens => ATMOS_REFSTATE_dens, &
       REF_pott => ATMOS_REFSTATE_pott
    use mod_atmos_boundary, only: &
       DAMP_var   => ATMOS_BOUNDARY_var,   &
       DAMP_alpha => ATMOS_BOUNDARY_alpha, &
       I_BND_VELZ,  &
       I_BND_VELX,  &
       I_BND_VELY,  &
       I_BND_POTT
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    implicit none

    real(RP) :: QDRY(KA,IA,JA)      ! dry air mixing ratio [kg/kg]
    real(RP) :: DDIV(KA,IA,JA)      ! divergence

    !---------------------------------------------------------------------------

#ifdef DEBUG
    QDRY(:,:,:) = UNDEF
    DDIV(:,:,:) = UNDEF
#endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamics step'

    call ATMOS_DYN_main( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,        & ! (inout)
         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (out)
         QDRY, DDIV,                                & ! (out)
         CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY,        & ! (in)
         CZ, FZ, CDZ, CDX, CDY, FDZ, FDX, FDY,      & ! (in)
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,        & ! (in)
         qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & !(in)
         qflx_sgs_rhot, qflx_sgs_qtrc,              & !(in)
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY,           & ! (in)
         SFLX_POTT, SFLX_QV,                        & ! (in)
         REF_dens, REF_pott, DIFF4,                 & ! (in)
         CORIOLI, DAMP_var, DAMP_alpha,             & ! (in)
         ATMOS_DYN_divdmp_coef, ATMOS_DYN_LSsink_D, & ! (in)
         ATMOS_DYN_FLAG_FCT_rho,                    & ! (in)
         ATMOS_DYN_FLAG_FCT_momentum,               & ! (in)
         ATMOS_DYN_FLAG_FCT_T,                      & ! (in)
         ATMOS_USE_AVERAGE,                         & ! (in)
         DTSEC, DTSEC_ATMOS_DYN, NSTEP_ATMOS_DYN    ) ! (in)

    call ATMOS_vars_total

    call HIST_in( QDRY(:,:,:), 'QDRY', 'Dry Air mixng ratio', 'kg/kg', '3D', DTSEC )
    call HIST_in( DDIV(:,:,:), 'div',  'Divergence',          's-1',   '3D', DTSEC )
    return

  end subroutine ATMOS_DYN

  subroutine ATMOS_DYN_main( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC,          & ! (inout)
         DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, & ! (out)
         QDRY, DDIV,                                  & ! (out)
         CNDZ, CNMZ, CNDX, CNMX, CNDY, CNMY,          & ! (in)
         CZ, FZ, CDZ, CDX, CDY, FDZ, FDX, FDY,        & ! (in)
         RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          & ! (in)
         qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & !(in)
         qflx_sgs_rhot, qflx_sgs_qtrc,                & !(in)
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY,             & ! (in)
         SFLX_POTT, SFLX_QV,                          & ! (in)
         REF_dens, REF_pott, DIFF4,                   & ! (in)
         corioli, DAMP_var, DAMP_alpha,               & ! (in)
         divdmp_coef, LSsink_D,                       & ! (in)
         FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, & ! (in)
         USE_AVERAGE,                                 & ! (in)
         DTSEC, DTSEC_ATMOS_DYN, NSTEP_ATMOS_DYN      ) ! (in)
    use mod_const, only : &
       Rdry   => CONST_Rdry,   &
       Rvap   => CONST_Rvap,   &
       P00    => CONST_PRE00
    use mod_comm, only: &
#ifdef _USE_RDMA
       COMM_rdma_vars8, &
#endif
       COMM_vars8, &
       COMM_wait
    use mod_atmos_boundary, only: &
       I_BND_VELZ,  &
       I_BND_VELX,  &
       I_BND_VELY,  &
       I_BND_POTT
    use dc_types, only: &
         DP
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

    real(RP), intent(out)   :: QDRY(KA,IA,JA)
    real(RP), intent(out)   :: DDIV(KA,IA,JA)
    real(RP), intent(in)    :: CNDZ(3,KA)
    real(RP), intent(in)    :: CNMZ(3,KA)
    real(RP), intent(in)    :: CNDX(3,IA)
    real(RP), intent(in)    :: CNMX(3,IA)
    real(RP), intent(in)    :: CNDY(3,JA)
    real(RP), intent(in)    :: CNMY(3,JA)
    real(RP), intent(in)    :: CZ(KA)
    real(RP), intent(in)    :: FZ(0:KA)
    real(RP), intent(in)    :: CDZ(KA)
    real(RP), intent(in)    :: CDX(IA)
    real(RP), intent(in)    :: CDY(JA)
    real(RP), intent(in)    :: FDZ(KA-1)
    real(RP), intent(in)    :: FDX(IA-1)
    real(RP), intent(in)    :: FDY(JA-1)
    real(RP), intent(in)    :: RCDZ(KA)
    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)
    real(RP), intent(in)    :: RFDZ(KA-1)
    real(RP), intent(in)    :: RFDX(IA-1)
    real(RP), intent(in)    :: RFDY(JA-1)

    real(RP), intent(in)    :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(in)    :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(in)    :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(in)    :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(in)    :: qflx_sgs_qtrc(KA,IA,JA,QA,3)

    real(RP), intent(in)    :: SFLX_MOMZ(IA,JA)
    real(RP), intent(in)    :: SFLX_MOMX(IA,JA)
    real(RP), intent(in)    :: SFLX_MOMY(IA,JA)
    real(RP), intent(in)    :: SFLX_POTT(IA,JA)
    real(RP), intent(in)    :: SFLX_QV  (IA,JA)

    real(RP), intent(in)    :: REF_dens(KA)
    real(RP), intent(in)    :: REF_pott(KA)
    real(RP), intent(in)    :: DIFF4
    real(RP), intent(in)    :: CORIOLI(1,IA,JA)
    real(RP), intent(in)    :: DAMP_var  (KA,IA,JA,5)
    real(RP), intent(in)    :: DAMP_alpha(KA,IA,JA,5)
    real(RP), intent(in)    :: divdmp_coef
    real(RP), intent(in)    :: LSsink_D

    logical,  intent(in)    :: FLAG_FCT_RHO
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T

    logical,  intent(in)    :: USE_AVERAGE

    real(DP), intent(in)    :: DTSEC
    real(DP), intent(in)    :: DTSEC_ATMOS_DYN
    integer , intent(in)    :: NSTEP_ATMOS_DYN

    ! diagnostic variables
    real(RP) :: PRES  (KA,IA,JA) ! pressure [Pa]
    real(RP) :: VELZ  (KA,IA,JA) ! velocity w [m/s]
    real(RP) :: VELX  (KA,IA,JA) ! velocity u [m/s]
    real(RP) :: VELY  (KA,IA,JA) ! velocity v [m/s]
    real(RP) :: POTT  (KA,IA,JA) ! potential temperature [K]
    real(RP) :: dens_s(KA,IA,JA) ! saved density
    real(RP) :: Rtot  (KA,IA,JA) ! R for dry air + vapor

    ! rayleigh damping, numerical diffusion
    real(RP) :: dens_diff(KA,IA,JA)     ! anomary of density
    real(RP) :: pott_diff(KA,IA,JA)     ! anomary of rho * pott
    real(RP) :: ray_damp (KA,IA,JA,5)
    real(RP) :: num_diff (KA,IA,JA,5,3)

    real(RP) :: qflx_hi  (KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order

    ! For FCT
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux
    real(RP) :: RHOQ     (KA,IA,JA)    ! rho(previous) * phi(previous)

    integer :: IIS, IIE
    integer :: JJS, JJE

    real(RP) :: dtrk
    integer :: i, j, k, iq, rko, step


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

    PRES    (:,:,:) = UNDEF
    VELZ    (:,:,:) = UNDEF
    VELX    (:,:,:) = UNDEF
    VELY    (:,:,:) = UNDEF
    POTT    (:,:,:) = UNDEF

    dens_s   (:,:,:)     = UNDEF
    dens_diff(:,:,:)     = UNDEF
    pott_diff(:,:,:)     = UNDEF
    ray_damp (:,:,:,:)   = UNDEF
    num_diff (:,:,:,:,:) = UNDEF

    mflx_hi  (:,:,:,:)   = UNDEF
    qflx_hi  (:,:,:,:)   = UNDEF
    qflx_lo  (:,:,:,:)   = UNDEF

    RHOQ   (:,:,:) = UNDEF

    rjpls(:,:,:) = UNDEF
    rjmns(:,:,:) = UNDEF
#endif

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
       MOMZ_RK1( 1:KS-1,i,j) = 0.0_RP
       MOMZ_RK1(KE:KA  ,i,j) = 0.0_RP
       MOMZ_RK2( 1:KS-1,i,j) = 0.0_RP
       MOMZ_RK2(KE:KA  ,i,j) = 0.0_RP
    enddo
    enddo

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
    enddo
    enddo

    if ( USE_AVERAGE ) then

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMZ_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMX_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMY_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT_av(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QTRC_av(k,i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    enddo

    end if

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,ZDIR) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,XDIR) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,YDIR) = 0.0_RP
    enddo
    enddo
    enddo

    do step = 1, NSTEP_ATMOS_DYN

#ifdef _FPCOLL_
call TIME_rapstart   ('DYN-set')
call START_COLLECTION("DYN-set")
#endif


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! Gas constant
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          QDRY(k,i,j) = 1.0_RP
       enddo
       enddo
       enddo
       do iq = QQS, QQE
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iq)
       enddo
       enddo
       enddo
       enddo
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          Rtot(k,i,j) = Rdry*QDRY(k,i,j) + Rvap*QTRC(k,i,j,I_QV)
       enddo
       enddo
       enddo
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          dens_s(k,i,j) = DENS(k,i,j)
       enddo
       enddo
       enddo

       !--- prepare rayleigh damping coefficient
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS, KE-1
             ray_damp(k,i,j,I_MOMZ) = - DAMP_alpha(k,i,j,I_BND_VELZ) &
                                    * ( MOMZ(k,i,j) - DAMP_var(k,i,j,I_BND_VELZ) * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) ) )
          enddo
          do k = KS, KE
             ray_damp(k,i,j,I_MOMX) = - DAMP_alpha(k,i,j,I_BND_VELX) &
                                    * ( MOMX(k,i,j) - DAMP_var(k,i,j,I_BND_VELX) * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) )
          enddo
          do k = KS, KE
             ray_damp(k,i,j,I_MOMY) = - DAMP_alpha(k,i,j,I_BND_VELY) &
                                    * ( MOMY(k,i,j) - DAMP_var(k,i,j,I_BND_VELY) * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) ) )
          enddo
          do k = KS, KE
             ray_damp(k,i,j,I_RHOT) = - DAMP_alpha(k,i,j,I_BND_POTT) &
                                    * ( RHOT(k,i,j) - DAMP_var(k,i,j,I_BND_POTT) * DENS(k,i,j) )
          enddo
       enddo
       enddo

       !--- prepare numerical diffusion coefficient
       if ( DIFF4 /= 0.0_RP ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j  = JJS-2, JJE+2
          do i  = IIS-2, IIE+2
          do k = KS, KE
             dens_diff(k,i,j) = DENS(k,i,j)               - REF_dens(k)
             pott_diff(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) - REF_pott(k)
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
          do k = KS+1, KE-2
             num_diff(k,i,j,I_DENS,ZDIR) = DIFF4 * CDZ(k)**4 &
                                         * ( CNDZ(1,k+1) * dens_diff(k+2,i,j) &
                                           - CNDZ(2,k+1) * dens_diff(k+1,i,j) &
                                           + CNDZ(3,k+1) * dens_diff(k  ,i,j) &
                                           - CNDZ(1,k  ) * dens_diff(k-1,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
             num_diff(KS  ,i,j,I_DENS,ZDIR) = DIFF4 * CDZ(KS  )**4 &
                                         * ( CNDZ(1,KS+1) * dens_diff(KS+2,i,j) &
                                           - CNDZ(2,KS+1) * dens_diff(KS+1,i,j) &
                                           + CNDZ(3,KS+1) * dens_diff(KS  ,i,j) &
                                           - CNDZ(1,KS  ) * dens_diff(KS,i,j) )
             num_diff(KE-1,i,j,I_DENS,ZDIR) = DIFF4 * CDZ(KE-1)**4 &
                                         * ( CNDZ(1,KE  ) * dens_diff(KE  ,i,j) &
                                           - CNDZ(2,KE  ) * dens_diff(KE  ,i,j) &
                                           + CNDZ(3,KE  ) * dens_diff(KE-1,i,j) &
                                           - CNDZ(1,KE-1) * dens_diff(KE-2,i,j) )
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE
             num_diff(k,i,j,I_DENS,XDIR) = DIFF4 * CDX(i)**4 &
                                         * ( CNDX(1,i+1) * dens_diff(k,i+2,j) &
                                           - CNDX(2,i+1) * dens_diff(k,i+1,j) &
                                           + CNDX(3,i+1) * dens_diff(k,i  ,j) &
                                           - CNDX(1,i  ) * dens_diff(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE
             num_diff(k,i,j,I_DENS,YDIR) = DIFF4 * CDY(j)**4 &
                                         * ( CNDY(1,j+1) * dens_diff(k,i,j+2) &
                                           - CNDY(2,j+1) * dens_diff(k,i,j+1) &
                                           + CNDY(3,j+1) * dens_diff(k,i,j  ) &
                                           - CNDY(1,j  ) * dens_diff(k,i,j-1) )
          enddo
          enddo
          enddo

          ! z-momentum
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
          do k = KS+1, KE-1
             num_diff(k,i,j,I_MOMZ,ZDIR) = DIFF4 * ( 0.5_RP*(CDZ(k+1)+CDZ(k)) )**4 &
                                         * ( CNMZ(1,k  ) * MOMZ(k+1,i,j) &
                                           - CNMZ(2,k  ) * MOMZ(k  ,i,j) &
                                           + CNMZ(3,k  ) * MOMZ(k-1,i,j) &
                                           - CNMZ(1,k-1) * MOMZ(k-2,i,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE-1
             num_diff(k,i,j,I_MOMZ,XDIR) = DIFF4 * CDX(i)**4 &
                                         * ( CNDX(1,i+1) * MOMZ(k,i+2,j) &
                                           - CNDX(2,i+1) * MOMZ(k,i+1,j) &
                                           + CNDX(3,i+1) * MOMZ(k,i  ,j) &
                                           - CNDX(1,i  ) * MOMZ(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE-1
             num_diff(k,i,j,I_MOMZ,YDIR) = DIFF4 * CDY(j)**4 &
                                         * ( CNDY(1,j+1) * MOMZ(k,i,j+2) &
                                           - CNDY(2,j+1) * MOMZ(k,i,j+1) &
                                           + CNDY(3,j+1) * MOMZ(k,i,j  ) &
                                           - CNDY(1,j  ) * MOMZ(k,i,j-1) )
          enddo
          enddo
          enddo

          ! x-momentum
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
          do k = KS+1, KE-2
             num_diff(k,i,j,I_MOMX,ZDIR) = DIFF4 * CDZ(k)**4 &
                                         * ( CNDZ(1,k+1) * MOMX(k+2,i,j) &
                                           - CNDZ(2,k+1) * MOMX(k+1,i,j) &
                                           + CNDZ(3,k+1) * MOMX(k  ,i,j) &
                                           - CNDZ(1,k  ) * MOMX(k-1,i,j) )
          enddo
          enddo
          enddo
          do j = JJS,   JJE
          do i = IIS,   IIE
             num_diff(KS  ,i,j,I_MOMX,ZDIR) = DIFF4 * CDZ(KS  )**4 &
                                         * ( CNDZ(1,KS+1) * MOMX(KS+2,i,j) &
                                           - CNDZ(2,KS+1) * MOMX(KS+1,i,j) &
                                           + CNDZ(3,KS+1) * MOMX(KS  ,i,j) &
                                           - CNDZ(1,KS  ) * MOMX(KS  ,i,j) )
             num_diff(KE-1,i,j,I_MOMX,ZDIR) = DIFF4 * CDZ(KE-1)**4 &
                                         * ( CNDZ(1,KE  ) * MOMX(KE  ,i,j) &
                                           - CNDZ(2,KE  ) * MOMX(KE  ,i,j) &
                                           + CNDZ(3,KE  ) * MOMX(KE-1,i,j) &
                                           - CNDZ(1,KE-1) * MOMX(KE-2,i,j) )
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE+1
          do k = KS, KE
             num_diff(k,i,j,I_MOMX,XDIR) = DIFF4 * ( 0.5_RP*(CDX(i+1)+CDX(i)) )**4 &
                                         * ( CNMX(1,i  ) * MOMX(k,i+1,j) &
                                           - CNMX(2,i  ) * MOMX(k,i  ,j) &
                                           + CNMX(3,i  ) * MOMX(k,i-1,j) &
                                           - CNMX(1,i-1) * MOMX(k,i-2,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMX,YDIR) = DIFF4 * CDY(j)**4 &
                                         * ( CNDY(1,j+1) * MOMX(k,i,j+2) &
                                           - CNDY(2,j+1) * MOMX(k,i,j+1) &
                                           + CNDY(3,j+1) * MOMX(k,i,j  ) &
                                           - CNDY(1,j  ) * MOMX(k,i,j-1) )
          enddo
          enddo
          enddo

          ! y-momentum
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
          do k = KS+1, KE-2
             num_diff(k,i,j,I_MOMY,ZDIR) = DIFF4 * CDZ(k)**4 &
                                         * ( CNDZ(1,k+1) * MOMY(k+2,i,j) &
                                           - CNDZ(2,k+1) * MOMY(k+1,i,j) &
                                           + CNDZ(3,k+1) * MOMY(k  ,i,j) &
                                           - CNDZ(1,k  ) * MOMY(k-1,i,j) )
          enddo
          enddo
          enddo
          do j = JJS,   JJE
          do i = IIS,   IIE
             num_diff(KS  ,i,j,I_MOMY,ZDIR) = DIFF4 * CDZ(KE  )**4 &
                                         * ( CNDZ(1,KS+1) * MOMY(KS+2,i,j) &
                                           - CNDZ(2,KS+1) * MOMY(KS+1,i,j) &
                                           + CNDZ(3,KS+1) * MOMY(KS  ,i,j) &
                                           - CNDZ(1,KS  ) * MOMY(KS  ,i,j) )
             num_diff(KE-1,i,j,I_MOMY,ZDIR) = DIFF4 * CDZ(KE-1)**4 &
                                         * ( CNDZ(1,KE  ) * MOMY(KE  ,i,j) &
                                           - CNDZ(2,KE  ) * MOMY(KE  ,i,j) &
                                           + CNDZ(3,KE  ) * MOMY(KE-1,i,j) &
                                           - CNDZ(1,KE-1) * MOMY(KE-2,i,j) )
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMY,XDIR) = DIFF4 * CDX(i)**4 &
                                         * ( CNDX(1,i+1) * MOMY(k,i+2,j) &
                                           - CNDX(2,i+1) * MOMY(k,i+1,j) &
                                           + CNDX(3,i+1) * MOMY(k,i  ,j) &
                                           - CNDX(1,i  ) * MOMY(k,i-1,j) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE+1
          do i = IIS, IIE
          do k = KS, KE
             num_diff(k,i,j,I_MOMY,YDIR) = DIFF4 * ( 0.5_RP*(CDY(j+1)+CDY(j)) )**4 &
                                         * ( CNMY(1,j  ) * MOMY(k,i,j+1) &
                                           - CNMY(2,j  ) * MOMY(k,i,j  ) &
                                           + CNMY(3,j  ) * MOMY(k,i,j-1) &
                                           - CNMY(1,j-1) * MOMY(k,i,j-2) )
          enddo
          enddo
          enddo

          ! rho * theta
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
          do k = KS+1, KE-2
             num_diff(k,i,j,I_RHOT,ZDIR) = DIFF4 * CDZ(k)**4 &
                                         * ( CNDZ(1,k+1) * pott_diff(k+2,i,j)   &
                                           - CNDZ(2,k+1) * pott_diff(k+1,i,j)   &
                                           + CNDZ(3,k+1) * pott_diff(k  ,i,j)   &
                                           - CNDZ(1,k  ) * pott_diff(k-1,i,j) ) &
                                         * 0.5_RP * ( FACT_N * ( DENS(k+1,i,j)+DENS(k  ,i,j) ) &
                                                   + FACT_F * ( DENS(k+2,i,j)+DENS(k-1,i,j) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS,   IIE
             num_diff(KS  ,i,j,I_RHOT,ZDIR) = DIFF4 * CDZ(KS  )**4 &
                                         * ( CNDZ(1,KS+1) * pott_diff(KS+2,i,j)   &
                                           - CNDZ(2,KS+1) * pott_diff(KS+1,i,j)   &
                                           + CNDZ(3,KS+1) * pott_diff(KS  ,i,j)   &
                                           - CNDZ(1,KS  ) * pott_diff(KS  ,i,j) ) &
                                         * 0.5_RP * ( DENS(KS+1,i,j)+DENS(KS  ,i,j) )
             num_diff(KE-1,i,j,I_RHOT,ZDIR) = DIFF4 * CDZ(KE-1)**4 &
                                         * ( CNDZ(1,KE  ) * pott_diff(KE  ,i,j)   &
                                           - CNDZ(2,KE  ) * pott_diff(KE  ,i,j)   &
                                           + CNDZ(3,KE  ) * pott_diff(KE-1,i,j)   &
                                           - CNDZ(1,KE-1) * pott_diff(KE-2,i,j) ) &
                                         * 0.5_RP * ( DENS(KE  ,i,j)+DENS(KE-1,i,j) )
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE
          do i = IIS-1, IIE
          do k = KS, KE
             num_diff(k,i,j,I_RHOT,XDIR) = DIFF4 * CDX(i)**4 &
                                         * ( CNDX(1,i+1) * pott_diff(k,i+2,j)   &
                                           - CNDX(2,i+1) * pott_diff(k,i+1,j)   &
                                           + CNDX(3,i+1) * pott_diff(k,i  ,j)   &
                                           - CNDX(1,i  ) * pott_diff(k,i-1,j) ) &
                                         * 0.5_RP * ( FACT_N * ( DENS(k,i+1,j)+DENS(k,i  ,j) ) &
                                                    + FACT_F * ( DENS(k,i+2,j)+DENS(k,i-1,j) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE
          do k = KS, KE
             num_diff(k,i,j,I_RHOT,YDIR) = DIFF4 * CDY(j)**4 &
                                         * ( CNDY(1,j+1) * pott_diff(k,i,j+2)   &
                                           - CNDY(2,j+1) * pott_diff(k,i,j+1)   &
                                           + CNDY(3,j+1) * pott_diff(k,i,j  )   &
                                           - CNDY(1,j  ) * pott_diff(k,i,j-1) ) &
                                         * 0.5_RP * ( FACT_N * ( DENS(k,i,j+1)+DENS(k,i,j  ) ) &
                                                    + FACT_F * ( DENS(k,i,j+2)+DENS(k,i,j-1) ) )
          enddo
          enddo
          enddo

       else
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS  , JJE
          do i = IIS  , IIE
          do k = KS   , KE
             num_diff(k,i,j,I_DENS,ZDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMZ,ZDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMX,ZDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMY,ZDIR) = 0.0_RP
             num_diff(k,i,j,I_RHOT,ZDIR) = 0.0_RP
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS  , JJE
          do i = IIS-1, IIE
          do k = KS   , KE
             num_diff(k,i,j,I_DENS,XDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMZ,XDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMY,XDIR) = 0.0_RP
             num_diff(k,i,j,I_RHOT,XDIR) = 0.0_RP
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS  , JJE
          do i = IIS  , IIE+1
          do k = KS   , KE
             num_diff(k,i,j,I_MOMX,XDIR) = 0.0_RP
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS  , IIE
          do k = KS   , KE
             num_diff(k,i,j,I_DENS,YDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMZ,YDIR) = 0.0_RP
             num_diff(k,i,j,I_MOMX,YDIR) = 0.0_RP
             num_diff(k,i,j,I_RHOT,YDIR) = 0.0_RP
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS  , JJE+1
          do i = IIS  , IIE
          do k = KS   , KE
             num_diff(k,i,j,I_MOMY,YDIR) = 0.0_RP
          enddo
          enddo
          enddo
       end if

    enddo
    enddo ! end tile

#ifdef _FPCOLL_
call STOP_COLLECTION ("DYN-set")
call TIME_rapend     ('DYN-set')
call TIME_rapstart   ('DYN-rk3')
call START_COLLECTION("DYN-rk3")
#endif

    !##### Start RK #####


    !##### RK1 #####
    rko = 1
    dtrk  = DTSEC_ATMOS_DYN / (RK - rko + 1)
    call calc_rk(DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1,  & ! (out)
                 DDIV,                                              & ! (out)
                 mflx_hi,                                           & ! (inout)
                 DENS,     MOMZ,     MOMX,     MOMY,     RHOT,      & ! (in)
                 DENS,     MOMZ,     MOMX,     MOMY,     RHOT,      & ! (in)
                 qflx_sgs_momz, qflx_sgs_momx,                      & ! (in)
                 qflx_sgs_momy, qflx_sgs_rhot,                      & ! (in)
                 SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT,        & ! (in)
                 Rtot, CORIOLI,                                     & ! (in)
                 num_diff, ray_damp, divdmp_coef, LSsink_d,         & ! (in)
                 FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,       & ! (in)
                 CZ, FZ, CDZ,                                       & ! (in)
                 FDZ, FDX, FDY, RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY, & ! (in)
                 dtrk, rko,                                         & ! (in)
                 VELZ, VELX, VELY, PRES, POTT,                      & ! (work)
                 qflx_hi, qflx_lo, qflx_anti, rjpls, rjmns          ) ! (work)
#ifdef _USE_RDMA
    call COMM_rdma_vars8( 5+QA+1, 5 )
#else
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
#endif

    !##### RK2 #####
    rko = 2
    dtrk  = DTSEC_ATMOS_DYN / (RK - rko + 1)
    call calc_rk(DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2,  & ! (out)
                 DDIV,                                              & ! (out)
                 mflx_hi,                                           & ! (inout)
                 DENS,     MOMZ,     MOMX,     MOMY,     RHOT,      & ! (in)
                 DENS_RK1, MOMZ_RK1, MOMX_RK1, MOMY_RK1, RHOT_RK1,  & ! (in)
                 qflx_sgs_momz, qflx_sgs_momx,                      & ! (in)
                 qflx_sgs_momy, qflx_sgs_rhot,                      & ! (in)
                 SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT,        & ! (in)
                 Rtot, CORIOLI,                                     & ! (in)
                 num_diff, ray_damp, divdmp_coef, LSsink_d,         & ! (in)
                 FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,       & ! (in)
                 CZ, FZ, CDZ,                                       & ! (in)
                 FDZ, FDX, FDY, RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY, & ! (in)
                 dtrk, rko,                                         & ! (in)
                 VELZ, VELX, VELY, PRES, POTT,                      & ! (work)
                 qflx_hi, qflx_lo, qflx_anti, rjpls, rjmns          ) ! (work)
#ifdef _USE_RDMA
    call COMM_rdma_vars8( 5+QA+6, 5 )
#else
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
#endif

    !##### RK3 #####
    rko = 3
    dtrk  = DTSEC_ATMOS_DYN / (RK - rko + 1)
    call calc_rk(DENS,     MOMZ,     MOMX,     MOMY,     RHOT,      & ! (out)
                 DDIV,                                              & ! (out)
                 mflx_hi,                                           & ! (inout)
                 DENS,     MOMZ,     MOMX,     MOMY,     RHOT,      & ! (in)
                 DENS_RK2, MOMZ_RK2, MOMX_RK2, MOMY_RK2, RHOT_RK2,  & ! (in)
                 qflx_sgs_momz, qflx_sgs_momx,                      & ! (in)
                 qflx_sgs_momy, qflx_sgs_rhot,                      & ! (in)
                 SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT,        & ! (in)
                 Rtot, CORIOLI,                                     & ! (in)
                 num_diff, ray_damp, divdmp_coef, LSsink_d,         & ! (in)
                 FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T,       & ! (in)
                 CZ, FZ, CDZ,                                       & ! (in)
                 FDZ, FDX, FDY, RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY, & ! (in)
                 dtrk, rko,                                         & ! (in)
                 VELZ, VELX, VELY, PRES, POTT,                      & ! (work)
                 qflx_hi, qflx_lo, qflx_anti, rjpls, rjmns          ) ! (work)
#ifdef _USE_RDMA
    call COMM_rdma_vars8( 1, 5 )
#else
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
#endif

    if ( USE_AVERAGE ) then

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS_av(k,i,j) = DENS_av(k,i,j) + DENS(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMZ_av(k,i,j) = MOMZ_av(k,i,j) + MOMZ(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMX_av(k,i,j) = MOMX_av(k,i,j) + MOMX(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMY_av(k,i,j) = MOMY_av(k,i,j) + MOMY(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT_av(k,i,j) = RHOT_av(k,i,j) + RHOT(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(4)
    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QTRC_av(k,i,j,iq) = QTRC_av(k,i,j,iq) + QTRC(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    end if

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,ZDIR) = mflx_av(k,i,j,ZDIR) + mflx_hi(k,i,j,ZDIR)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,XDIR) = mflx_av(k,i,j,XDIR) + mflx_hi(k,i,j,XDIR)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_av(k,i,j,YDIR) = mflx_av(k,i,j,YDIR) + mflx_hi(k,i,j,YDIR)
    enddo
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION ("DYN-rk3")
call TIME_rapend     ('DYN-rk3')
#endif

    enddo ! dynamical steps

#ifdef _FPCOLL_
call TIME_rapstart   ('DYN-fct')
call START_COLLECTION("DYN-fct")
#endif

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_hi(k,i,j,ZDIR) = mflx_av(k,i,j,ZDIR) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_hi(k,i,j,XDIR) = mflx_av(k,i,j,XDIR) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       mflx_hi(k,i,j,YDIR) = mflx_av(k,i,j,YDIR) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !##### advection of scalar quantity #####
    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-2
          qflx_lo(k,i,j,ZDIR) = 0.5_RP * (     mflx_hi(k,i,j,ZDIR)  * ( QTRC(k+1,i,j,iq)+QTRC(k,i,j,iq) ) &
                                         - abs(mflx_hi(k,i,j,ZDIR)) * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) ) &
                              + qflx_sgs_qtrc(k,i,j,iq,ZDIR)

          qflx_hi(k,i,j,ZDIR) = 0.5_RP * mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( QTRC(k+1,i,j,iq)+QTRC(k  ,i,j,iq) ) &
                                + FACT_F * ( QTRC(k+2,i,j,iq)+QTRC(k-1,i,j,iq) ) ) &
                              + qflx_sgs_qtrc(k,i,j,iq,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! Surface QV Flux
       if ( iq == I_QV ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             qflx_lo(KS-1,i,j,ZDIR) = SFLX_QV(i,j)
             qflx_hi(KS-1,i,j,ZDIR) = SFLX_QV(i,j)
          enddo
          enddo
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_lo(KS  ,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KS  ,i,j,ZDIR)  * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) & ! just above the bottom boundary
                                            - abs(mflx_hi(KS  ,i,j,ZDIR)) * ( QTRC(KS+1,i,j,iq)-QTRC(KS,i,j,iq) ) ) &
                                 + qflx_sgs_qtrc(KS,i,j,iq,ZDIR)

          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * mflx_hi(KS  ,i,j,ZDIR) * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) &
                                 + qflx_sgs_qtrc(KS,i,j,iq,ZDIR)
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          qflx_lo(KE-1,i,j,ZDIR) = 0.5_RP * (     mflx_hi(KE-1,i,j,ZDIR)  * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) & ! just below the top boundary
                                            - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) ) &
                                 + qflx_sgs_qtrc(KE-1,i,j,iq,ZDIR)

          qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary

          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR) * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) &
                                 + qflx_sgs_qtrc(KE-1,i,j,iq,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = mflx_hi(KE,i,j,ZDIR) &
               * ( QTRC(KE,i,j,iq) - 0.5_RP * ( QTRC(KE,i,j,iq) - QTRC(KE-1,i,j,iq) ) * FZ(KE) / FZ(KE-1) ) ! top boundary ( 0.0 if LSsink_D == 0 )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-2, IIE+1
       do k = KS, KE
          qflx_lo(k,i,j,XDIR) = 0.5_RP * (     mflx_hi(k,i,j,XDIR)  * ( QTRC(k,i+1,j,iq)+QTRC(k,i,j,iq) ) &
                                         - abs(mflx_hi(k,i,j,XDIR)) * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) ) &
                              + qflx_sgs_qtrc(k,i,j,iq,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_hi(k,i,j,XDIR) = 0.5_RP * mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( QTRC(k,i+1,j,iq)+QTRC(k,i  ,j,iq) ) &
                                + FACT_F * ( QTRC(k,i+2,j,iq)+QTRC(k,i-1,j,iq) ) ) &
                              + qflx_sgs_qtrc(k,i,j,iq,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-2, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          qflx_lo(k,i,j,YDIR) = 0.5_RP * (     mflx_hi(k,i,j,YDIR)  * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j,iq) ) &
                                         - abs(mflx_hi(k,i,j,YDIR)) * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) ) &
                              + qflx_sgs_qtrc(k,i,j,iq,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.5_RP * mflx_hi(k,i,j,YDIR) &
                              * ( FACT_N * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j  ,iq) ) &
                                + FACT_F * ( QTRC(k,i,j+2,iq)+QTRC(k,i,j-1,iq) ) ) &
                              + qflx_sgs_qtrc(k,i,j,iq,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update rho*Q with monotone(diffusive) flux divergence
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

    call fct(qflx_anti,               & ! (out)
             RHOQ, qflx_hi, qflx_lo,  & ! (in)
             RCDZ, RCDX, RCDY, DTSEC, & ! (in)
             rjpls, rjmns             ) ! (work)
#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       rjpls(:,:,:) = UNDEF
       rjmns(:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
                          - LSsink_D * QTRC(k,i,j,iq) ) & ! part of large scale sinking
                        ) / DENS(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo
#ifdef DEBUG
    qflx_hi  (:,:,:,:) = UNDEF
    qflx_anti(:,:,:,:) = UNDEF
#endif

    enddo ! scalar quantities loop

#ifdef _USE_RDMA
    call COMM_rdma_vars8( 6, QA )
#else
    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo
#endif

#ifdef _FPCOLL_
call STOP_COLLECTION ("DYN-fct")
call TIME_rapend     ('DYN-fct')
#endif

!OCL XFILL
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
!OCL XFILL
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

    if ( USE_AVERAGE ) then

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS_av(k,i,j) = DENS_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMZ_av(k,i,j) = MOMZ_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMX_av(k,i,j) = MOMX_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMY_av(k,i,j) = MOMY_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT_av(k,i,j) = RHOT_av(k,i,j) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(4)
    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QTRC_av(k,i,j,iq) = QTRC_av(k,i,j,iq) / NSTEP_ATMOS_DYN
    enddo
    enddo
    enddo
    enddo

    else

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS_av(k,i,j) = DENS(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMZ_av(k,i,j) = MOMZ(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMX_av(k,i,j) = MOMX(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       MOMY_av(k,i,j) = MOMY(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       RHOT_av(k,i,j) = RHOT(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) schedule(static,1) collapse(4)
    do iq = 1, QA
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       QTRC_av(k,i,j,iq) = QTRC(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    endif

    return
  end subroutine ATMOS_DYN_main


  subroutine calc_rk(DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
                     DDIV,                                        &
                     mflx_hi,                                     &
                     DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
                     DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
                     qflx_sgs_momz, qflx_sgs_momx,                &
                     qflx_sgs_momy, qflx_sgs_rhot,                &
                     SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT,  &
                     Rtot, CORIOLI,                               &
                     num_diff, ray_damp, divdmp_coef, LSsink_d,   &
                     FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
                     CZ, FZ, CDZ, FDZ, FDX, FDY,                  &
                     RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
                     dtrk, rko,                                   &
                     VELZ, VELX, VELY, PRES, POTT,                &
                     qflx_hi, qflx_lo, qflx_anti, rjpls, rjmns    )
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       Rvap   => CONST_Rvap,   &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00
    use mod_comm, only: &
#ifdef _USE_RDMA
       COMM_rdma_vars8, &
#endif
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out) :: DDIV(KA,IA,JA)      ! divergence

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! rho * vel(x,y,z)

    real(RP), intent(in),target :: DENS0(KA,IA,JA)   ! prognostic variables
    real(RP), intent(in),target :: MOMZ0(KA,IA,JA)   ! at previous dynamical time step
    real(RP), intent(in),target :: MOMX0(KA,IA,JA)   !
    real(RP), intent(in),target :: MOMY0(KA,IA,JA)   !
    real(RP), intent(in),target :: RHOT0(KA,IA,JA)   !

    real(RP), intent(in) :: DENS(KA,IA,JA)   ! prognostic variables
    real(RP), intent(in) :: MOMZ(KA,IA,JA)   ! at previous RK step
    real(RP), intent(in) :: MOMX(KA,IA,JA)   !
    real(RP), intent(in) :: MOMY(KA,IA,JA)   !
    real(RP), intent(in) :: RHOT(KA,IA,JA)   !

    real(RP), intent(in) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_sgs_rhot(KA,IA,JA,3)

    real(RP), intent(in) :: SFLX_MOMZ(IA,JA)
    real(RP), intent(in) :: SFLX_MOMX(IA,JA)
    real(RP), intent(in) :: SFLX_MOMY(IA,JA)
    real(RP), intent(in) :: SFLX_POTT(IA,JA)

    real(RP), intent(in) :: Rtot(KA,IA,JA) ! R for dry air + vapor
    real(RP), intent(in) :: CORIOLI(1,IA,JA)
    real(RP), intent(in) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in) :: ray_damp(KA,IA,JA,5)
    real(RP), intent(in) :: divdmp_coef
    real(RP), intent(in) :: LSsink_D

    logical,  intent(in) :: FLAG_FCT_RHO
    logical,  intent(in) :: FLAG_FCT_MOMENTUM
    logical,  intent(in) :: FLAG_FCT_T

    real(RP), intent(in) :: CZ(KA)
    real(RP), intent(in) :: FZ(0:KA)

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: RFDZ(KA-1)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RFDY(JA-1)
    real(RP), intent(in) :: dtrk
    integer , intent(in) :: rko

    ! diagnostic variables (work space)
    real(RP), intent(out) :: PRES(KA,IA,JA) ! pressure [Pa]
    real(RP), intent(out) :: VELZ(KA,IA,JA) ! velocity w [m/s]
    real(RP), intent(out) :: VELX(KA,IA,JA) ! velocity u [m/s]
    real(RP), intent(out) :: VELY(KA,IA,JA) ! velocity v [m/s]
    real(RP), intent(out) :: POTT(KA,IA,JA) ! potential temperature [K]

    ! flux (work space)
    real(RP), intent(out)   :: qflx_hi  (KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_lo  (KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_anti(KA,IA,JA,3)

    ! factor for FCT (work space)
    real(RP), intent(out) :: rjpls(KA,IA,JA)
    real(RP), intent(out) :: rjmns(KA,IA,JA)

    real(RP), pointer :: org(:,:,:)
    real(RP), target  :: work(KA,IA,JA)

    real(RP) :: rdtrk
    real(RP) :: vel
    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    mflx_hi(KS-1,:,:,ZDIR) = 0.0_RP

    qflx_lo(:,:,:,:) = UNDEF
    qflx_hi(:,:,:,:) = UNDEF

    rjpls(:,:,:) = UNDEF
    rjmns(:,:,:) = UNDEF

    work(:,:,:) = UNDEF
#endif

    rdtrk = 1.0_RP / dtrk

    if ( FLAG_FCT_RHO ) then
       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = DENS0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => DENS0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+2
       do i = IIS-1, IIE+2
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
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+2
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
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-2, JJE+1
       do i = IIS-1, IIE+2
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
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! pressure, pott. temp.
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
#if defined(__INTEL_COMPILER)
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, Rtot(k,i,j) )
#endif
             PRES(k,i,j) = RHOT(k,i,j) * Rtot(k,i,j) / P00
          enddo
          if ( RP == 8 ) then
#ifdef DEBUG
             call vdpowx( KE-KS+1, PRES(KS:KE,i,j), CPovCV, PRES(KS:KE,i,j) )
#else
             call vdpowx( KE, PRES(:,i,j), CPovCV, PRES(:,i,j) )
#endif
          else
#ifdef DEBUG
             call vspowx( KE-KS+1, PRES(KS:KE,i,j), CPovCV, PRES(KS:KE,i,j) )
#else
             call vspowx( KE, PRES(:,i,j), CPovCV, PRES(:,i,j) )
#endif
          end if
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS(k,i,j) )
#endif
             PRES(k,i,j) = PRES(k,i,j) * P00
             POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          enddo
#else
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, Rtot(k,i,j) )
#endif
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**CPovCV
          enddo
          do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, DENS(k,i,j) )
#endif
             POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          enddo
#endif
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! 3D divergence for damping
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k  ,i  ,j  ) )
          call CHECK( __LINE__, MOMZ(k-1,i  ,j  ) )
          call CHECK( __LINE__, MOMX(k  ,i  ,j  ) )
          call CHECK( __LINE__, MOMX(k  ,i-1,j  ) )
          call CHECK( __LINE__, MOMY(k  ,i  ,j  ) )
          call CHECK( __LINE__, MOMY(k  ,i  ,j-1) )
#endif
            DDIV(k,i,j) = ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * RCDZ(k) &
                        + ( MOMX(k,i,j) - MOMX(k  ,i-1,j  ) ) * RCDX(i) &
                        + ( MOMY(k,i,j) - MOMY(k  ,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### continuity equation #####

       ! monotonic flux
       if ( FLAG_FCT_RHO ) then
          ! at (x, y, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, DENS(k  ,i,j) )
             call CHECK( __LINE__, DENS(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  VELZ(k,i,j)  * ( DENS(k+1,i,j)+DENS(k,i,j) ) &
             -abs(VELZ(k,i,j)) * ( DENS(k+1,i,j)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, y, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k,i,j) )
             call CHECK( __LINE__, DENS(k,i  ,j) )
             call CHECK( __LINE__, DENS(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.5_RP * ( &
                  VELX(k,i,j)  * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
             -abs(VELX(k,i,j)) * ( DENS(k,i+1,j)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELY(k,i,j) )
             call CHECK( __LINE__, DENS(k,i,j  ) )
             call CHECK( __LINE__, DENS(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.5_RP * ( &
                  VELY(k,i,j)  * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
             -abs(VELY(k,i,j)) * ( DENS(k,i,j+1)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       endif

       ! at (x, y, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,ZDIR) )
#endif
          mflx_hi(k,i,j,ZDIR) = ( -MOMZ(k+1,i,j) + 8.0_RP*MOMZ(k,i,j) - MOMZ(k-1,i,j) ) / 6.0_RP &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KS,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_DENS,ZDIR) )
#endif
          mflx_hi(KS,i,j,ZDIR) = MOMZ(KS,i,j) &
                              + num_diff(KS,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          mflx_hi(KE,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k,i+1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
          mflx_hi(k,i,j,XDIR) = ( -MOMX(k,i+1,j) + 8.0_RP*MOMX(k,i,j) -MOMX(k,i-1,j) ) / 6.0_RP &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = ( -MOMY(k,i,j+1) + 8.0_RP*MOMY(k,i,j) -MOMY(k,i,j-1) ) / 6.0_RP &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update density
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS0(k,i,j) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, mflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j-1,YDIR) )
#endif
          DENS_RK(k,i,j) = DENS0(k,i,j) &
               + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                            + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                            + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) ) ! divergence
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    enddo
    enddo

    if ( FLAG_FCT_RHO ) then
       call fct(qflx_anti,          & ! (out)
            org, mflx_hi, qflx_lo,  & ! (in)
            RCDZ, RCDX, RCDY, dtrk, & ! (in)
            rjpls, rjmns            ) ! (work)
#ifdef DEBUG
       qflx_lo(:,:,:,:) = UNDEF
       rjpls(:,:,:) = UNDEF
       rjmns(:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update rho
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, DENS_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             DENS_RK(k,i,j) = DENS_RK(k,i,j) &
                  + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                               + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                               + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !--- update mflx_hi
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k,i,j,ZDIR) )
             call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
             call CHECK( __LINE__, qflx_anti(k,i,j,XDIR) )
             call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
             call CHECK( __LINE__, qflx_anti(k,i,j,YDIR) )
#endif
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + qflx_anti(k,i,j,ZDIR)
             mflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) + qflx_anti(k,i,j,XDIR)
             mflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) + qflx_anti(k,i,j,YDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       enddo
       enddo
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
#endif

#ifdef _USE_RDMA
    call COMM_rdma_vars8( 5+QA+11, 3 )
#else
    call COMM_vars8( mflx_hi(:,:,:,ZDIR), 1 )
    call COMM_vars8( mflx_hi(:,:,:,XDIR), 2 )
    call COMM_vars8( mflx_hi(:,:,:,YDIR), 3 )
    call COMM_wait ( mflx_hi(:,:,:,ZDIR), 1 )
    call COMM_wait ( mflx_hi(:,:,:,XDIR), 2 )
    call COMM_wait ( mflx_hi(:,:,:,YDIR), 3 )
#endif

    ! add momentum flux corresponding to large scale sinking
    if ( LSsink_D .ne. 0.0_RP ) then
       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
#endif
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                  - LSsink_D * FZ(k) ! large scale sinking
          enddo
          enddo
          enddo
       enddo
       enddo
    end if

#ifdef DEBUG
    work(:,:,:) = UNDEF
#endif
    if ( FLAG_FCT_MOMENTUM ) then
       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = MOMZ0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => MOMZ0
       end if
    end if


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (z) #####

       if ( FLAG_FCT_MOMENTUM ) then

          ! monotonic flux
          ! at (x, y, layer)
          ! note than z-index is added by -1
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS+1, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k  ,i,j) )
             call CHECK( __LINE__, DENS(k,i,j) )
#endif
             vel = ( ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) * 0.5_RP &
                     - LSsink_D * CZ(k) & ! part of large scale sinking
                   ) / DENS(k,i,j)
             qflx_lo(k-1,i,j,ZDIR) = 0.5_RP * ( &
                       ( vel ) * ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) &
                  - abs( vel ) * ( MOMZ(k,i,j)-MOMZ(k-1,i,j) ) ) &
                                   + qflx_sgs_momz(k,i,j,ZDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             ! k = KS
             qflx_lo(KS-1,i,j,ZDIR) = SFLX_MOMZ(i,j)
             ! k = KE
             qflx_lo(KE-1,i,j,ZDIR) = 0.0_RP
             ! k = KE+1
             qflx_lo(KE,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, y, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k  ,i,j) )
             call CHECK( __LINE__, VELX(k+1,i,j) )
             call CHECK( __LINE__, MOMZ(k,i  ,j) )
             call CHECK( __LINE__, MOMZ(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.25_RP * ( &
                       ( VELX(k+1,i,j)+VELX(k,i,j) ) * ( MOMZ(k,i+1,j)+MOMZ(k,i,j) ) &
                  - abs( VELX(k+1,i,j)+VELX(k,i,j) ) * ( MOMZ(k,i+1,j)-MOMZ(k,i,j) ) ) &
                                 + qflx_sgs_momz(k,i,j,XDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
             qflx_lo(KE,i,j,XDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, v, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELY(k  ,i,j) )
             call CHECK( __LINE__, VELY(k+1,i,j) )
             call CHECK( __LINE__, MOMZ(k,i,j  ) )
             call CHECK( __LINE__, MOMZ(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.25_RP * ( &
                      ( VELY(k+1,i,j)+VELY(k,i,j) ) * ( MOMZ(k,i,j+1)+MOMZ(k,i,j) ) &
                  -abs( VELY(k+1,i,j)+VELY(k,i,j) ) * ( MOMZ(k,i,j+1)-MOMZ(k,i,j) ) ) &
                                 + qflx_sgs_momz(k,i,j,YDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE,i,j,YDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       ! high order flux
       ! at (x, y, layer)
       ! note than z-index is added by -1
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k  ,i,j) )
          call CHECK( __LINE__, VELZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k-2,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,ZDIR) )
#endif
          vel = ( ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) * 0.5_RP &
                  - LSsink_D * CZ(k) & ! part of large scale sinking
                ) / DENS(k,i,j)
          qflx_hi(k-1,i,j,ZDIR) = 0.5_RP * vel  &
                              * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                                + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) ) &
                              + qflx_sgs_momz(k,i,j,ZDIR) &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, num_diff(KS+1,i,j,I_MOMZ,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMZ,ZDIR) )
#endif
          ! k = KS
          qflx_hi(KS-1,i,j,ZDIR) = SFLX_MOMZ(i,j) ! surface flux
          ! k = KS+1
          vel = ( ( MOMZ(KS+1,i,j)+MOMZ(KS,i,j) ) * 0.5_RP &
                  - LSsink_D * CZ(KS+1) & ! part of large scale sinking
                ) / DENS(KS,i,j)
          qflx_hi(KS,i,j,ZDIR) = 0.5_RP * vel &
                              * ( MOMZ(KS+1,i,j)+MOMZ(KS,i,j) ) &
                              + qflx_sgs_momz(KS+1,i,j,ZDIR) &
                              + num_diff(KS+1,i,j,I_MOMZ,ZDIR) * rdtrk
          ! k = KE-1
          vel = ( ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) * 0.5_RP &
                  - LSsink_D * CZ(KE-1) & ! part of large scale sinking
                ) / DENS(KE-1,i,j)
          qflx_hi(KE-2,i,j,ZDIR) = 0.5_RP * vel &
                              * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                              + qflx_sgs_momz(KE-1,i,j,ZDIR) &
                              + num_diff(KE-1,i,j,I_MOMZ,ZDIR) * rdtrk
          ! k = KE
          qflx_hi(KE-1,i,j,ZDIR) = - LSsink_D &
               * ( 0.5_RP * CZ(KE-1) * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) / DENS(KE-1,i,j) &
                 + 2.0_RP * FDZ(KE-1) * MOMZ(KE-1,i,j) / ( DENS(KE,i,j)+DENS(KE-1,i,j) ) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k  ,i,j) )
          call CHECK( __LINE__, VELX(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k,i-1,j) )
          call CHECK( __LINE__, MOMZ(k,i  ,j) )
          call CHECK( __LINE__, MOMZ(k,i+1,j) )
          call CHECK( __LINE__, MOMZ(k,i+2,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = 0.25_RP * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMZ(k,i+1,j)+MOMZ(k,i  ,j) ) &
                                + FACT_F * ( MOMZ(k,i+2,j)+MOMZ(k,i-1,j) ) ) &
                              + qflx_sgs_momz(k,i,j,XDIR) &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
          qflx_hi(KE,i,j,XDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k  ,i,j) )
          call CHECK( __LINE__, VELY(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k,i,j-1) )
          call CHECK( __LINE__, MOMZ(k,i,j  ) )
          call CHECK( __LINE__, MOMZ(k,i,j+1) )
          call CHECK( __LINE__, MOMZ(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = 0.25_RP * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMZ(k,i,j+1)+MOMZ(k,i,j  ) ) &
                                + FACT_F * ( MOMZ(k,i,j+2)+MOMZ(k,i,j-1) ) ) &
                              + qflx_sgs_momz(k,i,j,YDIR) &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
          qflx_hi(KE,i,j,YDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(z)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, PRES(k  ,i,j) )
          call CHECK( __LINE__, PRES(k+1,i,j) )
          call CHECK( __LINE__, DENS(k  ,i,j) )
          call CHECK( __LINE__, DENS(k+1,i,j) )
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, ray_damp(k,i,j,i_MOMZ) )
#endif
          MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
               + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                            + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                          - ( PRES(k+1,i,j)-PRES(k,i,j) ) * RFDZ(k)                      & ! pressure gradient force
                          - ( DENS(k+1,i,j)+DENS(k,i,j) ) * 0.5_RP * GRAV                & ! gravity force
                          + divdmp_coef * dtrk  * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) & ! divergence damping
                          - 2.0_RP * LSsink_D * MOMZ(k,i,j) / ( DENS(k+1,i,j) + DENS(k,i,j) ) & ! part of large scale sinking
                          + ray_damp(k,i,j,I_MOMZ)                                       ) ! additional damping
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       call fct(qflx_anti,              & ! (out)
                org, qflx_hi, qflx_lo,  & ! (in)
                RFDZ, RCDX, RCDY, dtrk, & ! (in)
                rjpls, rjmns            ) ! (work)
#ifdef DEBUG
       qflx_lo(:,:,:,:) = UNDEF
       rjpls(:,:,:) = UNDEF
       rjmns(:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(z)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             MOMZ_RK(k,i,j) = MOMZ_RK(k,i,j) &
                  + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                               + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                               + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if

#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
    work(:,:,:) = UNDEF
#endif
    if ( FLAG_FCT_MOMENTUM ) then

       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = MOMX0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => MOMX0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (x) #####

       ! monotonic flux
       if ( FLAG_FCT_MOMENTUM ) then
          ! at (u, y, interface)
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i  ,j) )
             call CHECK( __LINE__, VELZ(k,i+1,j) )
             call CHECK( __LINE__, MOMX(k  ,i,j) )
             call CHECK( __LINE__, MOMX(k+1,i,j) )
             call CHECK( __LINE__, DENS(k+1,i+1,j) )
             call CHECK( __LINE__, DENS(k+1,i  ,j) )
             call CHECK( __LINE__, DENS(k  ,i+1,j) )
             call CHECK( __LINE__, DENS(k  ,i  ,j) )
#endif
             vel = 0.5_RP * ( VELZ(k,i+1,j) + VELZ(k,i,j) ) &
                 - 4.0_RP * LSsink_D * FZ(k) / ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) ) ! large scale sinking
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  ( vel ) * ( MOMX(k+1,i,j)+MOMX(k,i,j) ) &
              -abs( vel ) * ( MOMX(k+1,i,j)-MOMX(k,i,j) ) ) &
                                 + qflx_sgs_momx(k,i,j,ZDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = SFLX_MOMX(i,j)
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, y, layer)
          ! note that x-index is added by -1
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+2
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k,i-1,j) )
             call CHECK( __LINE__, MOMX(k,i  ,j) )
             call CHECK( __LINE__, DENS(k,i,j) )
#endif
             vel = ( MOMX(k,i,j)+MOMX(k,i-1,j) ) / DENS(k,i,j)
             qflx_lo(k,i-1,j,XDIR) = 0.25_RP * ( &
                  ( vel ) * ( MOMX(k,i,j)+MOMX(k,i-1,j) ) &
              -abs( vel ) * ( MOMX(k,i,j)-MOMX(k,i-1,j) ) ) &
                                   + qflx_sgs_momx(k,i,j,XDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELY(k,i+1,j) )
             call CHECK( __LINE__, VELY(k,i  ,j) )
             call CHECK( __LINE__, MOMX(k,i,j  ) )
             call CHECK( __LINE__, MOMX(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.25_RP * ( &
                  ( VELY(k,i+1,j)+VELY(k,i,j) ) * ( MOMX(k,i,j+1)+MOMX(k,i,j) ) &
              -abs( VELY(k,i+1,j)+VELY(k,i,j) ) * ( MOMX(k,i,j+1)-MOMX(k,i,j) ) ) &
                                 + qflx_sgs_momx(k,i,j,YDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if
       ! high order flux
       ! at (u, y, interface)
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k,i  ,j) )
          call CHECK( __LINE__, VELZ(k,i+1,j) )
          call CHECK( __LINE__, DENS(k+1,i+1,j) )
          call CHECK( __LINE__, DENS(k+1,i  ,j) )
          call CHECK( __LINE__, DENS(k  ,i+1,j) )
          call CHECK( __LINE__, DENS(k  ,i  ,j) )
          call CHECK( __LINE__, MOMX(k-1,i,j) )
          call CHECK( __LINE__, MOMX(k  ,i,j) )
          call CHECK( __LINE__, MOMX(k+1,i,j) )
          call CHECK( __LINE__, MOMX(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,ZDIR) )
#endif
          vel = 0.5_RP * ( VELZ(k,i+1,j) + VELZ(k,i,j) ) &
              - 4.0_RP * LSsink_D * FZ(k) / ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) ) ! large scale sinking
          qflx_hi(k,i,j,ZDIR) = 0.5_RP * vel &
                              * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                                + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) &
                              + qflx_sgs_momx(k,i,j,ZDIR) &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS,i  ,j) )
          call CHECK( __LINE__, VELZ(KS,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i  ,j) )
          call CHECK( __LINE__, DENS(KS  ,i+1,j) )
          call CHECK( __LINE__, DENS(KS  ,i  ,j) )
          call CHECK( __LINE__, MOMX(KS+1,i,j) )
          call CHECK( __LINE__, MOMX(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMX,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i  ,j) )
          call CHECK( __LINE__, VELZ(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i  ,j) )
          call CHECK( __LINE__, DENS(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMX(KE-1,i,j) )
          call CHECK( __LINE__, MOMX(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMX,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = SFLX_MOMX(i,j)
          vel = 0.5_RP * ( VELZ(KS,i+1,j) + VELZ(KS,i,j) ) &
              - 4.0_RP * LSsink_D * FZ(KS) / ( DENS(KS+1,i+1,j)+DENS(KS+1,i,j)+DENS(KS,i+1,j)+DENS(KS,i,j) ) ! large scale sinking
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * vel & ! just above the bottom boundary
                                 * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) &
                                 + qflx_sgs_momx(KS  ,i,j,ZDIR) &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          vel = 0.5_RP * ( VELZ(KE-1,i+1,j) + VELZ(KE-1,i,j) ) &
              - 4.0_RP * LSsink_D * FZ(KE-1) / ( DENS(KE,i+1,j)+DENS(KE,i,j)+DENS(KE-1,i+1,j)+DENS(KE-1,i,j) ) ! large scale sinking
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * vel & ! just below the top boundary
                                 * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) ) &
                                 + qflx_sgs_momx(KE-1,i,j,ZDIR) &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE,i,j,ZDIR) = - 2.0_RP * LSsink_D &
               * ( FZ(KE-1) * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) ) / ( DENS(KE,i+1,j)+DENS(KE,i,j)+DENS(KE-1,i+1,j)+DENS(KE-1,i,j) ) &
                 + CDZ(KE) * MOMX(KE,i,j) / ( DENS(KE,i+1,j)+DENS(KE,i,j) ) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, y, layer)
       ! note that x-index is added by -1
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k,i  ,j) )
          call CHECK( __LINE__, VELX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i-2,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,XDIR) )
#endif
          qflx_hi(k,i-1,j,XDIR) = 0.25_RP * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                              * ( FACT_N * ( MOMX(k,i  ,j)+MOMX(k,i-1,j) ) &
                                + FACT_F * ( MOMX(k,i+1,j)+MOMX(k,i-2,j) ) ) &
                              + qflx_sgs_momx(k,i,j,XDIR) &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i,j-1) )
          call CHECK( __LINE__, MOMX(k,i,j  ) )
          call CHECK( __LINE__, MOMX(k,i,j+1) )
          call CHECK( __LINE__, MOMX(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = 0.25_RP * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMX(k,i,j+1)+MOMX(k,i,j  ) ) &
                                + FACT_F * ( MOMX(k,i,j+2)+MOMX(k,i,j-1) ) ) &
                              + qflx_sgs_momx(k,i,j,YDIR) &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(x)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, PRES(k,i+1,j) )
          call CHECK( __LINE__, PRES(k,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j  ) )
          call CHECK( __LINE__, VELY(k,i+1,j  ) )
          call CHECK( __LINE__, VELY(k,i  ,j-1) )
          call CHECK( __LINE__, VELY(k,i+1,j-1) )
          call CHECK( __LINE__, DDIV(k,i+1,j) )
          call CHECK( __LINE__, DDIV(k,i  ,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
          call CHECK( __LINE__, ray_damp(k,i,j,I_MOMX) )
#endif
          MOMX_RK(k,i,j) = MOMX0(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                                       + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                                     - ( PRES(k,i+1,j)-PRES(k,i,j) ) * RFDX(i)                     & ! pressure gradient force
                                     + 0.125_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i+1,j) )              &
                                     * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) & ! coriolis force
                                     + divdmp_coef * dtrk * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) * FDX(i) & ! divergence damping
                                     - 2.0_RP * LSsink_D * MOMX(k,i,j) / ( DENS(k,i+1,j) + DENS(k,i,j) ) & ! part of large scale sinking
                                     + ray_damp(k,i,j,I_MOMX)                                      ) ! additional damping
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       call fct(qflx_anti,              & ! (out)
                org, qflx_hi, qflx_lo,  & ! (in)
                RCDZ, RFDX, RCDY, dtrk, & ! (in)
                rjpls, rjmns            ) ! (work)
#ifdef DEBUG
       qflx_lo(:,:,:,:) = UNDEF
       rjpls(:,:,:) = UNDEF
       rjmns(:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(x)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMX_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             MOMX_RK(k,i,j) = MOMX_RK(k,i,j) &
                  + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                               + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                               + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
    work(:,:,:) = UNDEF
#endif
    if ( FLAG_FCT_MOMENTUM ) then

       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = MOMY0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => MOMY0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (y) #####
       ! monotonic flux
       if ( FLAG_FCT_MOMENTUM ) then
          ! at (x, v, interface)
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i,j+1) )
             call CHECK( __LINE__, VELZ(k,i,j  ) )
             call CHECK( __LINE__, DENS(k+1,i+1,j) )
             call CHECK( __LINE__, DENS(k+1,i  ,j) )
             call CHECK( __LINE__, DENS(k  ,i+1,j) )
             call CHECK( __LINE__, DENS(k  ,i  ,j) )
             call CHECK( __LINE__, MOMY(k  ,i,j) )
             call CHECK( __LINE__, MOMY(k+1,i,j) )
#endif
             vel = 0.5_RP * ( VELZ(k,i,j+1) + VELZ(k,i,j) ) &
                 - 4.0_RP * LSsink_D * FZ(k) / ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) ) ! large scale sinking
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  ( vel ) * ( MOMY(k+1,i,j)+MOMY(k,i,j) ) &
              -abs( vel ) * ( MOMY(k+1,i,j)-MOMY(k,i,j) ) ) &
                                 + qflx_sgs_momy(k,i,j,ZDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = SFLX_MOMY(i,j)
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k,i,j+1) )
             call CHECK( __LINE__, VELX(k,i,j  ) )
             call CHECK( __LINE__, MOMY(k,i  ,j) )
             call CHECK( __LINE__, MOMY(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.25_RP * ( &
                  ( VELX(k,i,j+1)+VELX(k,i,j) ) * ( MOMY(k,i+1,j)+MOMY(k,i,j) ) &
              -abs( VELX(k,i,j+1)+VELX(k,i,j) ) * ( MOMY(k,i+1,j)-MOMY(k,i,j) ) ) &
                                 + qflx_sgs_momy(k,i,j,XDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, y, layer)
          ! note that y-index is added by -1
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+2
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i,j-1) )
             call CHECK( __LINE__, MOMY(k,i,j  ) )
             call CHECK( __LINE__, DENS(k,i,j) )
#endif
             vel = ( MOMY(k,i,j)+MOMY(k,i,j-1) ) / DENS(k,i,j)
             qflx_lo(k,i,j-1,YDIR) = 0.25_RP * ( &
                  ( vel ) * ( MOMY(k,i,j)+MOMY(k,i,j-1) ) &
              -abs( vel ) * ( MOMY(k,i,j)-MOMY(k,i,j-1) ) ) &
                                   + qflx_sgs_momy(k,i,j,YDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       ! high order flux
       ! at (x, v, interface)
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k,i,j+1) )
          call CHECK( __LINE__, VELZ(k,i,j  ) )
          call CHECK( __LINE__, DENS(k+1,i+1,j) )
          call CHECK( __LINE__, DENS(k+1,i  ,j) )
          call CHECK( __LINE__, DENS(k  ,i+1,j) )
          call CHECK( __LINE__, DENS(k  ,i  ,j) )
          call CHECK( __LINE__, MOMY(k-1,i,j) )
          call CHECK( __LINE__, MOMY(k  ,i,j) )
          call CHECK( __LINE__, MOMY(k+1,i,j) )
          call CHECK( __LINE__, MOMY(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,ZDIR) )
#endif
          vel = 0.5_RP * ( VELZ(k,i,j+1) + VELZ(k,i,j) ) &
              - 4.0_RP * LSsink_D * FZ(k) / ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) ) ! large scale sinking
          qflx_hi(k,i,j,ZDIR) = 0.5_RP * vel &
                              * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                                + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) ) &
                              + qflx_sgs_momy(k,i,j,ZDIR) &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j+1) )
          call CHECK( __LINE__, VELZ(KS  ,i,j  ) )
          call CHECK( __LINE__, DENS(KS+1,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i  ,j) )
          call CHECK( __LINE__, DENS(KS  ,i+1,j) )
          call CHECK( __LINE__, DENS(KS  ,i  ,j) )
          call CHECK( __LINE__, MOMY(KS+1,i,j) )
          call CHECK( __LINE__, MOMY(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMY,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i,j+1) )
          call CHECK( __LINE__, VELZ(KE-1,i,j  ) )
          call CHECK( __LINE__, DENS(KE  ,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i  ,j) )
          call CHECK( __LINE__, DENS(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMY(KE-1,i,j) )
          call CHECK( __LINE__, MOMY(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMY,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = SFLX_MOMY(i,j)
          vel = 0.5_RP * ( VELZ(KS,i,j+1) + VELZ(KS,i,j) ) &
              - 4.0_RP * LSsink_D * FZ(KS) / ( DENS(KS+1,i,j+1)+DENS(KS+1,i,j)+DENS(KS,i,j+1)+DENS(KS,i,j) ) ! large scale sinking
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * vel & ! just above the bottom boundary
                                 * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) )              &
                                 + qflx_sgs_momy(KS  ,i,j,ZDIR) &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          vel = 0.5_RP * ( VELZ(KE-1,i,j+1) + VELZ(KE-1,i,j) ) &
              - 4.0_RP * LSsink_D * FZ(KE-1) / ( DENS(KE,i,j+1)+DENS(KE,i,j)+DENS(KE-1,i,j+1)+DENS(KE-1,i,j) ) ! large scale sinking
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * vel & ! just below the top boundary
                                 * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) )              &
                                 + qflx_sgs_momy(KE-1,i,j,ZDIR) &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = - 2.0_RP * LSsink_D &
               * ( FZ(KE-1) * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) ) / ( DENS(KE,i,j+1)+DENS(KE,i,j)+DENS(KE-1,i,j+1)+DENS(KE-1,i,j) ) &
                 + CDZ(KE) * MOMY(KE,i,j) / ( DENS(KE,i,j+1)+DENS(KE,i,j) ) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k,i,j+1) )
          call CHECK( __LINE__, VELX(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i-1,j) )
          call CHECK( __LINE__, MOMY(k,i  ,j) )
          call CHECK( __LINE__, MOMY(k,i+1,j) )
          call CHECK( __LINE__, MOMY(k,i+2,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = 0.25_RP * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMY(k,i+1,j)+MOMY(k,i  ,j) ) &
                                + FACT_F * ( MOMY(k,i+2,j)+MOMY(k,i-1,j) ) ) &
                              + qflx_sgs_momy(k,i,j,XDIR) &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, y, layer)
       ! note that y-index is added by -1
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k,i,j  ) )
          call CHECK( __LINE__, VELY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j-2) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,YDIR) )
#endif
          qflx_hi(k,i,j-1,YDIR) = 0.25_RP * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                              * ( FACT_N * ( MOMY(k,i,j  )+MOMY(k,i,j-1) ) &
                                + FACT_F * ( MOMY(k,i,j+1)+MOMY(k,i,j-2) ) ) &
                              + qflx_sgs_momy(k,i,j,YDIR) &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(y)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, PRES(k,i,j  ) )
          call CHECK( __LINE__, PRES(k,i,j+1) )
          call CHECK( __LINE__, CORIOLI(1,i,j  ) )
          call CHECK( __LINE__, CORIOLI(1,i,j+1) )
          call CHECK( __LINE__, VELX(k,i  ,j  ) )
          call CHECK( __LINE__, VELX(k,i  ,j+1) )
          call CHECK( __LINE__, VELX(k,i-1,j  ) )
          call CHECK( __LINE__, VELX(k,i-1,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j  ) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
          call CHECK( __LINE__, ray_damp(k,i,j,I_MOMY) )
#endif
          MOMY_RK(k,i,j) = MOMY0(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RFDY(j) ) &
                                     - ( PRES(k,i,j+1)-PRES(k,i,j) ) * RFDY(j)                     & ! pressure gradient force
                                     - 0.125_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i,j+1) )              &
                                     * ( VELX(k,i,j)+VELX(k,i,j+1)+VELX(k,i-1,j)+VELX(k,i-1,j+1) ) & ! coriolis force
                                     + divdmp_coef * dtrk * ( DDIV(k,i,j+1)-DDIV(k,i,j) ) * FDY(j) & ! divergence damping
                                     - 2.0_RP * LSsink_D * MOMY(k,i,j) / ( DENS(k,i,j+1) + DENS(k,i,j) ) & ! part of large scale sinking
                                     + ray_damp(k,i,j,I_MOMY)                                      ) ! additional damping
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       call fct(qflx_anti,              & ! (out)
                org, qflx_hi, qflx_lo,  & ! (in)
                RCDZ, RCDX, RFDY, dtrk, & ! (in)
                rjpls, rjmns            ) ! (work)
#ifdef DEBUG
       qflx_lo(:,:,:,:) = UNDEF
       rjpls(:,:,:) = UNDEF
       rjmns(:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(y)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             MOMY_RK(k,i,j) = MOMY_RK(k,i,j) &
                  + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                               + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                               + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RFDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
    work(:,:,:) = UNDEF
#endif

    if ( FLAG_FCT_T ) then
       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => RHOT0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### Thermodynamic Equation #####

       ! monotonic flux
       if ( FLAG_FCT_T ) then
          ! at (x, y, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
             call CHECK( __LINE__, POTT(k  ,i,j) )
             call CHECK( __LINE__, POTT(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  mflx_hi(k,i,j,ZDIR)  * ( POTT(k+1,i,j)+POTT(k,i,j) ) &
             -abs(mflx_hi(k,i,j,ZDIR)) * ( POTT(k+1,i,j)-POTT(k,i,j) ) ) &
                                 + qflx_sgs_rhot(k,i,j,ZDIR)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = SFLX_POTT(i,j)
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo

#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, y, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
             call CHECK( __LINE__, POTT(k,i  ,j) )
             call CHECK( __LINE__, POTT(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.5_RP * ( &
                  mflx_hi(k,i,j,XDIR)  * ( POTT(k,i+1,j)+POTT(k,i,j) ) &
             -abs(mflx_hi(k,i,j,XDIR)) * ( POTT(k,i+1,j)-POTT(k,i,j) ) )  &
                                 + qflx_sgs_rhot(k,i,j,XDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
             call CHECK( __LINE__, POTT(k,i,j  ) )
             call CHECK( __LINE__, POTT(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.5_RP * ( &
                  mflx_hi(k,i,j,YDIR)  * ( POTT(k,i,j+1)+POTT(k,i,j) ) &
             -abs(mflx_hi(k,i,j,YDIR)) * ( POTT(k,i,j+1)-POTT(k,i,j) ) ) &
                                 + qflx_sgs_rhot(k,i,j,YDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       ! at (x, y, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(k-1,i,j) )
          call CHECK( __LINE__, POTT(k  ,i,j) )
          call CHECK( __LINE__, POTT(k+1,i,j) )
          call CHECK( __LINE__, POTT(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,ZDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = 0.5_RP * mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j)+POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j)+POTT(k-1,i,j) ) ) &
                              + qflx_sgs_rhot(k,i,j,ZDIR) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(KS,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_RHOT,ZDIR) )
          call CHECK( __LINE__, mflx_hi(KE-1,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_RHOT,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = SFLX_POTT(i,j)
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * mflx_hi(KS  ,i,j,ZDIR)  &      ! just above the bottom boundary
                                 * ( POTT(KS+1,i,j)+POTT(KS,i,j) ) &
                                 + qflx_sgs_rhot(KS  ,i,j,ZDIR) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR)  &      ! just below the top boundary
                                 * ( POTT(KE,i,j)+POTT(KE-1,i,j) ) &
                                 + qflx_sgs_rhot(KE-1,i,j,ZDIR) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = mflx_hi(KE,i,j,ZDIR) &
               * ( POTT(KE,i,j) - 0.5_RP * ( POTT(KE,i,j) - POTT(KE-1,i,j) ) * FZ(KE) / FZ(KE-1) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = 0.5_RP * mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) &
                              + qflx_sgs_rhot(k,i,j,XDIR) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = 0.5_RP * mflx_hi(k,i,j,YDIR) &
                              * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) &
                              + qflx_sgs_rhot(k,i,j,YDIR) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update rho*theta
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, RHOT0(k,i,j) )
          call CHECK( __LINE__, ray_damp(k,i,j,i_RHOT) )
#endif
          RHOT_RK(k,i,j) = RHOT0(k,i,j) &
               + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                          - LSsink_D * RHOT(k,i,j) / DENS(k,i,j) & ! part of large scale sinking
                          + ray_damp(k,i,j,I_RHOT)  ) ! additional damping
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    enddo
    enddo

    if ( FLAG_FCT_T ) then
       call fct(qflx_anti,               & ! (out)
                org, qflx_hi, qflx_lo, & ! (in)
                RCDZ, RCDX, RCDY, dtrk,  & ! (in)
                rjpls, rjmns             ) ! (work)
#ifdef DEBUG
       qflx_lo(:,:,:,:) = UNDEF
       rjpls(:,:,:) = UNDEF
       rjmns(:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update rho*theta
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             RHOT_RK(k,i,j) = RHOT_RK(k,i,j) &
                  + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                               + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                               + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
       qflx_hi(:,:,:,:) = UNDEF
#endif

  end subroutine calc_rk

  subroutine fct( qflx_anti,                &
                  phi_in, qflx_hi, qflx_lo, &
                  rdz, rdx, rdy, dtrk,      &
                  rjpls, rjmns              )
    use mod_comm, only: &
#ifdef _USE_RDMA
       COMM_rdma_vars8, &
#endif
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: qflx_anti(KA,IA,JA,3)

    real(RP), intent(in) :: phi_in(KA,IA,JA) ! physical quantity
    real(RP), intent(in) :: qflx_hi(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_lo(KA,IA,JA,3)

    real(RP), intent(in) :: RDZ(:)
    real(RP), intent(in) :: RDX(:)
    real(RP), intent(in) :: RDY(:)
    real(RP), intent(in) :: dtrk

    ! factor for FCT (work space)
    real(RP), intent(inout) :: rjpls(KA,IA,JA)
    real(RP), intent(inout) :: rjmns(KA,IA,JA)


    ! work for FCT
    real(RP) :: phi_lo(KA,IA,JA)
    real(RP) :: pjpls(KA,IA,JA)
    real(RP) :: pjmns(KA,IA,JA)
    real(RP) :: qjpls(KA,IA,JA)
    real(RP) :: qjmns(KA,IA,JA)

    real(RP) :: zerosw, dirsw

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j, ijs

#ifdef DEBUG
    qflx_anti(:,:,:,:) = UNDEF

    pjpls(:,:,:) = UNDEF
    pjmns(:,:,:) = UNDEF
    qjpls(:,:,:) = UNDEF
    qjmns(:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,ZDIR) )
#endif
          qflx_anti(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_lo(k,i,j,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,XDIR) )
#endif
          qflx_anti(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_lo(k,i,j,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,YDIR) )
#endif
          qflx_anti(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) - qflx_lo(k,i,j,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update monotone scheme
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, phi_in(k,i,j) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j-1,YDIR) )
#endif
          phi_lo(k,i,j) = phi_in(k,i,j) &
                          + dtrk * ( - ( ( qflx_lo(k,i,j,ZDIR)-qflx_lo(k-1,i  ,j  ,ZDIR) ) * RDZ(k) &
                                       + ( qflx_lo(k,i,j,XDIR)-qflx_lo(k  ,i-1,j  ,XDIR) ) * RDX(i) &
                                       + ( qflx_lo(k,i,j,YDIR)-qflx_lo(k  ,i  ,j-1,YDIR) ) * RDY(j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net incoming quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
          pjpls(k,i,j) = dtrk * ( ( max(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) - min(0.0_RP,qflx_anti(k,i,j,ZDIR)) ) * RDZ(k) &
                                + ( max(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) - min(0.0_RP,qflx_anti(k,i,j,XDIR)) ) * RDX(i) &
                                + ( max(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) - min(0.0_RP,qflx_anti(k,i,j,YDIR)) ) * RDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net outgoing quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
          pjmns(k,i,j) = dtrk * ( ( max(0.0_RP,qflx_anti(k,i,j,ZDIR)) - min(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) ) * RDZ(k) &
                                + ( max(0.0_RP,qflx_anti(k,i,j,XDIR)) - min(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) ) * RDX(i) &
                                + ( max(0.0_RP,qflx_anti(k,i,j,YDIR)) - min(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) ) * RDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc allowable range or quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, phi_in(k  ,i  ,j  ) )
          call CHECK( __LINE__, phi_in(k-1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(k+1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(k  ,i-1,j  ) )
          call CHECK( __LINE__, phi_in(k  ,i+1,j  ) )
          call CHECK( __LINE__, phi_in(k  ,i  ,j+1) )
          call CHECK( __LINE__, phi_in(k  ,i  ,j-1) )
          call CHECK( __LINE__, phi_lo(k  ,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(k-1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(k+1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(k  ,i-1,j  ) )
          call CHECK( __LINE__, phi_lo(k  ,i+1,j  ) )
          call CHECK( __LINE__, phi_lo(k  ,i  ,j+1) )
          call CHECK( __LINE__, phi_lo(k  ,i  ,j-1) )
#endif
          qjpls(k,i,j) = max( phi_in(k  ,i  ,j  ), &
                              phi_in(k+1,i  ,j  ), &
                              phi_in(k-1,i  ,j  ), &
                              phi_in(k  ,i+1,j  ), &
                              phi_in(k  ,i-1,j  ), &
                              phi_in(k  ,i  ,j+1), &
                              phi_in(k  ,i  ,j-1), &
                              phi_lo(k  ,i  ,j  ), &
                              phi_lo(k+1,i  ,j  ), &
                              phi_lo(k-1,i  ,j  ), &
                              phi_lo(k  ,i+1,j  ), &
                              phi_lo(k  ,i-1,j  ), &
                              phi_lo(k  ,i  ,j+1), &
                              phi_lo(k  ,i  ,j-1) ) &
                       - phi_lo(k,i,j)
          qjmns(k,i,j) = phi_lo(k,i,j) &
                       - min( phi_in(k  ,i  ,j  ), &
                              phi_in(k+1,i  ,j  ), &
                              phi_in(k-1,i  ,j  ), &
                              phi_in(k  ,i-1,j  ), &
                              phi_in(k  ,i+1,j  ), &
                              phi_in(k  ,i  ,j+1), &
                              phi_in(k  ,i  ,j-1), &
                              phi_lo(k  ,i  ,j  ), &
                              phi_lo(k+1,i  ,j  ), &
                              phi_lo(k-1,i  ,j  ), &
                              phi_lo(k  ,i-1,j  ), &
                              phi_lo(k  ,i+1,j  ), &
                              phi_lo(k  ,i  ,j+1), &
                              phi_lo(k  ,i  ,j-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, phi_in(KS  ,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KS+1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KS  ,i-1,j  ) )
          call CHECK( __LINE__, phi_in(KS  ,i+1,j  ) )
          call CHECK( __LINE__, phi_in(KS  ,i  ,j+1) )
          call CHECK( __LINE__, phi_in(KS  ,i  ,j-1) )
          call CHECK( __LINE__, phi_lo(KS  ,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KS+1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KS  ,i-1,j  ) )
          call CHECK( __LINE__, phi_lo(KS  ,i+1,j  ) )
          call CHECK( __LINE__, phi_lo(KS  ,i  ,j+1) )
          call CHECK( __LINE__, phi_lo(KS  ,i  ,j-1) )
          call CHECK( __LINE__, phi_in(KE  ,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KE-1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KE  ,i-1,j  ) )
          call CHECK( __LINE__, phi_in(KE  ,i+1,j  ) )
          call CHECK( __LINE__, phi_in(KE  ,i  ,j+1) )
          call CHECK( __LINE__, phi_in(KE  ,i  ,j-1) )
          call CHECK( __LINE__, phi_lo(KE  ,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KE-1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KE  ,i-1,j  ) )
          call CHECK( __LINE__, phi_lo(KE  ,i+1,j  ) )
          call CHECK( __LINE__, phi_lo(KE  ,i  ,j+1) )
          call CHECK( __LINE__, phi_lo(KE  ,i  ,j-1) )
#endif
          qjpls(KS,i,j) = max( phi_in(KS  ,i  ,j  ), &
                               phi_in(KS+1,i  ,j  ), &
                               phi_in(KS  ,i+1,j  ), &
                               phi_in(KS  ,i-1,j  ), &
                               phi_in(KS  ,i  ,j+1), &
                               phi_in(KS  ,i  ,j-1), &
                               phi_lo(KS  ,i  ,j  ), &
                               phi_lo(KS+1,i  ,j  ), &
                               phi_lo(KS  ,i+1,j  ), &
                               phi_lo(KS  ,i-1,j  ), &
                               phi_lo(KS  ,i  ,j+1), &
                               phi_lo(KS  ,i  ,j-1) ) &
                        - phi_lo(KS,i,j)
          qjmns(KS,i,j) = phi_lo(KS,i,j) &
                        - min( phi_in(KS  ,i  ,j  ), &
                               phi_in(KS+1,i  ,j  ), &
                               phi_in(KS  ,i-1,j  ), &
                               phi_in(KS  ,i+1,j  ), &
                               phi_in(KS  ,i  ,j+1), &
                               phi_in(KS  ,i  ,j-1), &
                               phi_lo(KS  ,i  ,j  ), &
                               phi_lo(KS+1,i  ,j  ), &
                               phi_lo(KS  ,i-1,j  ), &
                               phi_lo(KS  ,i+1,j  ), &
                               phi_lo(KS  ,i  ,j+1), &
                               phi_lo(KS  ,i  ,j-1) )
          qjpls(KE,i,j) = max( phi_in(KE  ,i  ,j  ), &
                               phi_in(KE-1,i  ,j  ), &
                               phi_in(KE  ,i+1,j  ), &
                               phi_in(KE  ,i-1,j  ), &
                               phi_in(KE  ,i  ,j+1), &
                               phi_in(KE  ,i  ,j-1), &
                               phi_lo(KE  ,i  ,j  ), &
                               phi_lo(KE-1,i  ,j  ), &
                               phi_lo(KE  ,i+1,j  ), &
                               phi_lo(KE  ,i-1,j  ), &
                               phi_lo(KE  ,i  ,j+1), &
                               phi_lo(KE  ,i  ,j-1) ) &
                        - phi_lo(KE,i,j)
          qjmns(KE,i,j) = phi_lo(KE,i,j) &
                        - min( phi_in(KE  ,i  ,j  ), &
                               phi_in(KE-1,i  ,j  ), &
                               phi_in(KE  ,i-1,j  ), &
                               phi_in(KE  ,i+1,j  ), &
                               phi_in(KE  ,i  ,j+1), &
                               phi_in(KE  ,i  ,j-1), &
                               phi_lo(KE  ,i  ,j  ), &
                               phi_lo(KE-1,i  ,j  ), &
                               phi_lo(KE  ,i-1,j  ), &
                               phi_lo(KE  ,i+1,j  ), &
                               phi_lo(KE  ,i  ,j+1), &
                               phi_lo(KE  ,i  ,j-1) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- incoming flux limitation factor [0-1]
       !$omp parallel do private(i,j,k,zerosw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, pjpls(k,i,j) )
          call CHECK( __LINE__, qjpls(k,i,j) )
#endif
          ! if pjpls == 0, zerosw = 1 and rjpls = 0
          zerosw = 0.5_RP - sign( 0.5_RP, pjpls(k,i,j)-EPSILON )
          rjpls(k,i,j) = min( 1.0_RP, qjpls(k,i,j) * ( 1.0_RP-zerosw ) / ( pjpls(k,i,j)-zerosw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- outgoing flux limitation factor [0-1]
       !$omp parallel do private(i,j,k,zerosw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, pjmns(k,i,j) )
          call CHECK( __LINE__, qjmns(k,i,j) )
#endif
          ! if pjmns == 0, zerosw = 1 and rjmns = 0
          zerosw = 0.5_RP - sign( 0.5_RP, pjmns(k,i,j)-EPSILON )
          rjmns(k,i,j) = min( 1.0_RP, qjmns(k,i,j) * ( 1.0_RP-zerosw ) / ( pjmns(k,i,j)-zerosw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

#ifdef _USE_RDMA
    call COMM_rdma_vars8( 5+QA+14, 2 )
#else
    call COMM_vars8( rjpls(:,:,:), 1 )
    call COMM_vars8( rjmns(:,:,:), 2 )
    call COMM_wait ( rjpls(:,:,:), 1 )
    call COMM_wait ( rjmns(:,:,:), 2 )
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !--- update high order flux with antidiffusive flux
       !$omp parallel do private(i,j,k,dirsw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS , KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,ZDIR) )
          call CHECK( __LINE__, rjpls(k  ,i,j) )
          call CHECK( __LINE__, rjpls(k+1,i,j) )
          call CHECK( __LINE__, rjmns(k  ,i,j) )
          call CHECK( __LINE__, rjmns(k+1,i,j) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,ZDIR) )
          qflx_anti(k,i,j,ZDIR) = qflx_anti(k,i,j,ZDIR) &
                 * ( min( rjpls(k+1,i,j),rjmns(k  ,i,j) ) * (          dirsw ) &
                   + min( rjpls(k  ,i,j),rjmns(k+1,i,j) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(KE,i,j,ZDIR) )
          call CHECK( __LINE__, rjpls(KE  ,i,j) )
          call CHECK( __LINE__, rjmns(KE  ,i,j) )
#endif
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP ! top    boundary
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       if ( IIS == IS ) then
          ijs = IIS-1
       else
          ijs = IIS
       end if
       !$omp parallel do private(i,j,k,dirsw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = ijs, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,XDIR) )
          call CHECK( __LINE__, rjpls(k,i  ,j) )
          call CHECK( __LINE__, rjpls(k,i+1,j) )
          call CHECK( __LINE__, rjmns(k,i  ,j) )
          call CHECK( __LINE__, rjmns(k,i+1,j) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,XDIR) )
          qflx_anti(k,i,j,XDIR) = qflx_anti(k,i,j,XDIR) &
                 * ( min( rjpls(k,i+1,j),rjmns(k,i  ,j) ) * (          dirsw ) &
                   + min( rjpls(k,i  ,j),rjmns(k,i+1,j) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       if ( JJS == JS ) then
          ijs = JJS-1
       else
          ijs = JJS
       end if
       !$omp parallel do private(i,j,k,dirsw) schedule(static,1) collapse(2)
       do j = ijs, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,YDIR) )
          call CHECK( __LINE__, rjpls(k,i,j+1) )
          call CHECK( __LINE__, rjpls(k,i,j  ) )
          call CHECK( __LINE__, rjmns(k,i,j  ) )
          call CHECK( __LINE__, rjmns(k,i,j+1) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,YDIR) )
          qflx_anti(k,i,j,YDIR) = qflx_anti(k,i,j,YDIR) &
                 * ( min( rjpls(k,i,j+1),rjmns(k,i,j  ) ) * (          dirsw ) &
                   + min( rjpls(k,i,j  ),rjmns(k,i,j+1) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    return
  end subroutine fct
end module mod_atmos_dyn
