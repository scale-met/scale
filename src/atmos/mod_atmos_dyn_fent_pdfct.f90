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
!! @li      2011-11-11 (H.Yashiro)  [new] Imported from SCALE-LES ver.2
!! @li      2011-11-11 (H.Yashiro)  [mod] Merged with Y.Miyamoto's
!! @li      2011-12-11 (H.Yashiro)  [mod] Use reference state
!! @li      2011-12-26 (Y.Miyamoto) [mod] Add numerical diffusion into mass flux calc
!! @li      2012-01-04 (H.Yashiro)  [mod] Nonblocking communication (Y.Ohno)
!! @li      2012-01-25 (H.Yashiro)  [fix] Bugfix (Y.Miyamoto)
!! @li      2012-01-25 (H.Yashiro)  [mod] Positive definite FCT (Y.Miyamoto)
!! @li      2012-02-14 (H.Yashiro)  [mod] Cache tiling
!! @li      2012-03-14 (H.Yashiro)  [mod] Bugfix (Y.Miyamoto)
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
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

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
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
  real(8), private, parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), private, parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second

  ! numerical filter settings
  real(8), private, save      :: ATMOS_DYN_numerical_diff = 1.D-2 ! nondimensional numerical diffusion
  real(8), private, save      :: DIFF4 ! for 4th order numerical filter
  real(8), private, save      :: DIFF2 ! for 2nd order numerical filter

  ! work
  real(8), private, save :: DENS_RK1(KA,IA,JA)   ! prognostic variables (+1/3 step)
  real(8), private, save :: MOMZ_RK1(KA,IA,JA)   !
  real(8), private, save :: MOMX_RK1(KA,IA,JA)   !
  real(8), private, save :: MOMY_RK1(KA,IA,JA)   !
  real(8), private, save :: RHOT_RK1(KA,IA,JA)   !
  real(8), private, save :: DENS_RK2(KA,IA,JA)   ! prognostic variables (+2/3 step)
  real(8), private, save :: MOMZ_RK2(KA,IA,JA)   !
  real(8), private, save :: MOMX_RK2(KA,IA,JA)   !
  real(8), private, save :: MOMY_RK2(KA,IA,JA)   !
  real(8), private, save :: RHOT_RK2(KA,IA,JA)   !

  real(8), private, save :: rjmns   (KA,IA,JA,3) ! minus in (x,y,z)-direction

  real(8), private, save :: CNDZ(3,KA)
  real(8), private, save :: CNMZ(3,KA)
  real(8), private, save :: CNDX(3,IA)
  real(8), private, save :: CNMX(3,IA)
  real(8), private, save :: CNDY(3,JA)
  real(8), private, save :: CNMY(3,JA)

  !-----------------------------------------------------------------------------
contains

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
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff

    integer :: ierr
    integer :: k, i, j
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

    DIFF4 = - ATMOS_DYN_numerical_diff * (-1.D0)**( 4/2+1 )
    DIFF2 = - ATMOS_DYN_numerical_diff * (-1.D0)**( 2/2+1 )

    !OCL XFILL
    do j = 1, JA
    do i = 1, IA
       rjmns(KS-1,i,j,ZDIR) = 0.D0
       rjmns(KS-1,i,j,XDIR) = 0.D0
       rjmns(KS-1,i,j,YDIR) = 0.D0
       rjmns(KE+1,i,j,ZDIR) = 0.D0
       rjmns(KE+1,i,j,XDIR) = 0.D0
       rjmns(KE+1,i,j,YDIR) = 0.D0
    enddo
    enddo

    ! z djrectjon
    do k = KS-1, KE+1
       CNDZ(1,k) = 1.D0 / ( (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 )
    enddo
    CNDZ(1,   1:KS-2) = CNDZ(1,KS-1)
    CNDZ(1,KE+2:KA  ) = CNDZ(1,KE+1)

    do k = KS-1, KE+1
       CNDZ(2,k) = 1.D0 / ( (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 )
    enddo
    CNDZ(2,   1:KS-2) = CNDZ(2,KS-1)
    CNDZ(2,KE+2:KA  ) = CNDZ(2,KE+1)

    do k = KS, KE+2
       CNDZ(3,k) = 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k-1)+CDZ(k-2)) * 0.5D0 )
    enddo
    CNDZ(3,1   :KS-1) = CNDZ(3,KS  )
    CNDZ(3,KE+2:KA  ) = CNDZ(3,KE+2)

    do k = KS-2, KE+1
       CNMZ(1,k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) )
    enddo
    CNMZ(1,   1:KS-2) = CNMZ(1,KS-2)
    CNMZ(1,KE+2:KA  ) = CNMZ(1,KE+1)

    do k = KS-1, KE+1
       CNMZ(2,k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &   
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &  
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) 
    enddo
    CNMZ(2,   1:KS-2) = CNMZ(2,KS-1)
    CNMZ(2,KE+2:KA  ) = CNMZ(2,KE+1)

    do k = KS-1, KE+1
       CNMZ(3,k) = 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) &
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) )
    enddo
    CNMZ(3,   1:KS-2) = CNMZ(3,KS-1)
    CNMZ(3,KE+2:KA  ) = CNMZ(3,KE+1)

    ! x direction
    do i = IS-1, IE+1
       CNDX(1,i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
    enddo
    CNDX(1,   1:IS-2) = CNDX(1,IS-1)
    CNDX(1,IE+2:IA  ) = CNDX(1,IE+1)

    do i = IS-1, IE+1
       CNDX(2,i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
    enddo
    CNDX(2,   1:IS-2) = CNDX(2,IS-1)
    CNDX(2,IE+2:IA  ) = CNDX(2,IE+1)

    do i = IS, IE+2
       CNDX(3,i) = 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i-1)+CDX(i-2)) * 0.5D0 )
    enddo
    CNDX(3,   1:IS-1) = CNDX(3,IS  )
    CNDX(3,IE+2:IA  ) = CNDX(3,IE+2)

    do i = IS-2, IE+1
       CNMX(1,i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) )
    enddo
    CNMX(1,IE+2:IA  ) = CNMX(1,IE+1)
    CNMX(1,   1:IS-2) = CNMX(1,IS-2)

    do i = IS-1, IE+1
       CNMX(2,i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) )
    enddo
    CNMX(2,   1:IS-2) = CNMX(2,IS-1)
    CNMX(2,IE+2:IA  ) = CNMX(2,IE+1)

    do i = IS-1, IE+1
       CNMX(3,i) = 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) )
    enddo
    CNMX(3,   1:IS-2) = CNMX(3,IS-1)
    CNMX(3,IE+2:IA  ) = CNMX(3,IE+1)

    ! y direction
    do j = JS-1, JE+1
       CNDY(1,j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
    enddo
    CNDY(1,   1:JS-2) = CNDY(1,JS-1)
    CNDY(1,JE+2:JA  ) = CNDY(1,JE+1)

    do j = JS-1, JE+1
       CNDY(2,j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
    enddo
    CNDY(2,   1:JS-2) = CNDY(2,JS-1)
    CNDY(2,JE+2:JA  ) = CNDY(2,JE+1)

    do j = JS, JE+2
       CNDY(3,j) = 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j-1)+CDY(j-2)) * 0.5D0 )
    enddo
    CNDY(3,   1:JS-1) = CNDY(3,JS  )
    CNDY(3,JE+2:JA  ) = CNDY(3,JE+2)

    do j = JS-2, JE+1
       CNMY(1,j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) )
    enddo
    CNMY(1,   1:JS-2) = CNMY(1,JS-2)
    CNMY(1,JE+2:JA  ) = CNMY(1,JE+1)

    do j = JS-1, JE+1
       CNMY(2,j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &   
                 + 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &  
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) 
    enddo
    CNMY(2,   1:JS-2) = CNMY(2,JS-1)
    CNMY(2,JE+2:JA  ) = CNMY(2,JE+1)

    do j = JS-1, JE+1
       CNMY(3,j) = 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) &
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) )
    enddo
    CNMY(3,   1:JS-2) = CNMY(3,JS-1)
    CNMY(3,JE+2:JA  ) = CNMY(3,JE+1)

  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       Rvap   => CONST_Rvap,   &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00
    use mod_time, only: &
       TIME_DTSEC_ATMOS_DYN, &
       TIME_NSTEP_ATMOS_DYN
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait,  &
       COMM_total
    use mod_grid, only : &
       CDZ  => GRID_CDZ,  &
       CDX  => GRID_CDX,  &
       CDY  => GRID_CDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use mod_atmos_vars, only: &
       ATMOS_vars_total,   &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
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
    implicit none

    ! diagnostic variables
    real(8) :: PRES  (KA,IA,JA) ! pressure [Pa]
    real(8) :: VELZ  (KA,IA,JA) ! velocity w [m/s]
    real(8) :: VELX  (KA,IA,JA) ! velocity u [m/s]
    real(8) :: VELY  (KA,IA,JA) ! velocity v [m/s]
    real(8) :: POTT  (KA,IA,JA) ! potential temperature [K]
    real(8) :: dens_s(KA,IA,JA) ! saved density
    real(8) :: QDRY  (KA,IA,JA) ! dry air mixing ratio [kg/kg]
    real(8) :: Rtot  (KA,IA,JA) ! R for dry air + vapor

    ! rayleigh damping, numerical diffusion
    real(8) :: dens_diff(KA,IA,JA)     ! anomary of density
    real(8) :: pott_diff(KA,IA,JA)     ! anomary of rho * pott
    real(8) :: ray_damp (KA,IA,JA,5)
    real(8) :: num_diff (KA,IA,JA,5,3)

    ! mass flux
    real(8) :: mflx_hi  (KA,IA,JA,3)   ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8) :: qflx_hi  (KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order

    ! For FCT
    real(8) :: qflx_lo  (KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face low  order
    real(8) :: qflx_anti(KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face antidiffusive
    real(8) :: pjmns, tmp

    integer :: IIS, IIE
    integer :: JJS, JJE

    real(8) :: dtrk, rdtrk
    integer :: i, j, k, iq, iw, rko, step
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("DYNAMICS")
#endif

    !OCL XFILL
    do j = JS, JE+1
    do i = IS, IE+1
       VELZ(KS-1,i,j) = 0.D0
       VELZ(KE  ,i,j) = 0.D0
    enddo
    enddo

    !OCL XFILL
    do j = 1, JA
    do i = 1, IA
       MOMZ_RK1( 1:KS-1,i,j) = 0.D0
       MOMZ_RK1(KE:KA  ,i,j) = 0.D0
       MOMZ_RK2( 1:KS-1,i,j) = 0.D0
       MOMZ_RK2(KE:KA  ,i,j) = 0.D0
    enddo
    enddo

    !OCL XFILL
    do j = JS, JE
    do i = IS, IE
       mflx_hi(KS-1,i,j,ZDIR) = 0.D0 ! bottom boundary
       mflx_hi(KE  ,i,j,ZDIR) = 0.D0 ! top    boundary
    enddo
    enddo

    do step = 1, TIME_NSTEP_ATMOS_DYN
!
!    DENS_RK1(:,:,:) = -9.999D30
!    MOMZ_RK1(:,:,:) = -9.999D30
!    MOMZ_RK1(1:KS-1,:,:) = 0.D0
!    MOMZ_RK1(KE:KA ,:,:) = 0.D0
!    MOMX_RK1(:,:,:) = -9.999D30
!    MOMY_RK1(:,:,:) = -9.999D30
!    RHOT_RK1(:,:,:) = -9.999D30
!    DENS_RK2(:,:,:) = -9.999D30
!    MOMZ_RK2(:,:,:) = -9.999D30
!    MOMZ_RK2(1:KS-1,:,:) = 0.D0
!    MOMZ_RK2(KE:KA ,:,:) = 0.D0
!    MOMX_RK2(:,:,:) = -9.999D30
!    MOMY_RK2(:,:,:) = -9.999D30
!    RHOT_RK2(:,:,:) = -9.999D30
!    PRES    (:,:,:) = -9.999D30
!    VELZ    (:,:,:) = -9.999D30
!    VELZ    (1:KS-1,:,:) = 0.D0
!    VELZ    (KE:KA ,:,:) = 0.D0
!    VELX    (:,:,:) = -9.999D30
!    VELY    (:,:,:) = -9.999D30
!    POTT    (:,:,:) = -9.999D30
!
!    dens_s   (:,:,:)     = -9.999D30
!    dens_diff(:,:,:)     = -9.999D30
!    pott_diff(:,:,:)     = -9.999D30
!    ray_damp (:,:,:,:)   = -9.999D30
!    num_diff (:,:,:,:,:) = -9.999D30
!
!    mflx_hi  (:,:,:,:)   = -9.999D30
!    mflx_hi  (KS-1,:,:,ZDIR) = 0.D0
!    mflx_hi  (KE  ,:,:,ZDIR) = 0.D0
!    qflx_hi  (:,:,:,:)   = -9.999D30
!    qflx_lo  (:,:,:,:)   = -9.999D30
!    qflx_anti(:,:,:,:)   = -9.999D30

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical small step:', step

#ifdef _FPCOLL_
call START_COLLECTION("SET")
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !--- prepare rayleigh damping coefficient
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS, KE-1
             ray_damp(k,i,j,I_MOMZ) = - DAMP_alpha(k,i,j,I_BND_VELZ) &
                                    * ( MOMZ(k,i,j) - DAMP_var(k,i,j,I_BND_VELZ) * 0.5D0 * ( DENS(k+1,i,j)+DENS(k,i,j) ) )
          enddo
          do k = KS, KE
             ray_damp(k,i,j,I_MOMX) = - DAMP_alpha(k,i,j,I_BND_VELX) &
                                    * ( MOMX(k,i,j) - DAMP_var(k,i,j,I_BND_VELX) * 0.5D0 * ( DENS(k,i+1,j)+DENS(k,i,j) ) )
             ray_damp(k,i,j,I_MOMY) = - DAMP_alpha(k,i,j,I_BND_VELY) &
                                    * ( MOMY(k,i,j) - DAMP_var(k,i,j,I_BND_VELY) * 0.5D0 * ( DENS(k,i,j+1)+DENS(k,i,j) ) )
             ray_damp(k,i,j,I_RHOT) = - DAMP_alpha(k,i,j,I_BND_POTT) &
                                    * ( RHOT(k,i,j) - DAMP_var(k,i,j,I_BND_POTT) * DENS(k,i,j) )
          enddo 
       enddo
       enddo

       !--- prepare numerical diffusion coefficient
       do j  = JJS-2, JJE+2
       do i  = IIS-2, IIE+2
          do k = KS, KE
             dens_s(k,i,j)    = DENS(k,i,j)
             dens_diff(k,i,j) = DENS(k,i,j)               - REF_dens(k)
             pott_diff(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) - REF_pott(k)
          enddo
          !OCL XFILL
          do k = 1, KS-1 
             dens_diff(k,i,j) = dens_diff(KS,i,j)
             pott_diff(k,i,j) = pott_diff(KS,i,j)
          enddo
          !OCL XFILL
          do k = KE+1, KA
             dens_diff(k,i,j) = dens_diff(KE,i,j)
             pott_diff(k,i,j) = pott_diff(KE,i,j)
          enddo
       enddo
       enddo

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

       do j = JJS,   JJE
       do i = IIS,   IIE
          num_diff(KS  ,i,j,I_DENS,ZDIR) = DIFF2 * CDZ(KS) &
                                         * 4.0D0 * ( dens_diff(KS+1,i,j)-dens_diff(KS,i,j) ) 
          num_diff(KE-1,i,j,I_DENS,ZDIR) = DIFF2 * CDZ(KE-1) &
                                         * 4.0D0 * ( dens_diff(KE,i,j)-dens_diff(KE-1,i,j) )
       enddo
       enddo

       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          num_diff(k,i,j,I_DENS,XDIR) = DIFF4 * CDX(i)**4 &
                                      * ( CNDX(1,i+1) * dens_diff(k,i+2,j) &
                                        - CNDX(2,i+1) * dens_diff(k,i+1,j) &
                                        + CNDX(3,i+1) * dens_diff(k,i  ,j) &
                                        - CNDX(1,i  ) * dens_diff(k,i-1,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          num_diff(k,i,j,I_DENS,YDIR) = DIFF4 * CDY(j)**4 &
                                      * ( CNDY(1,j+1) * dens_diff(k,i,j+2) &
                                        - CNDY(2,j+1) * dens_diff(k,i,j+1) &
                                        + CNDY(3,j+1) * dens_diff(k,i,j  ) &
                                        - CNDY(1,j  ) * dens_diff(k,i,j-1) )
       enddo
       enddo
       enddo

       ! z-momentum
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          num_diff(k,i,j,I_MOMZ,ZDIR) = DIFF4 * ( 0.5D0*(CDZ(k+1)+CDZ(k)) )**4 &
                                      * ( CNMZ(1,k  ) * MOMZ(k+1,i,j) &
                                        - CNMZ(2,k  ) * MOMZ(k  ,i,j) &
                                        + CNMZ(3,k  ) * MOMZ(k-1,i,j) &
                                        - CNMZ(1,k-1) * MOMZ(k-2,i,j) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS-1, KE
          num_diff(k,i,j,I_MOMZ,XDIR) = DIFF4 * CDX(i)**4 &
                                      * ( CNDX(1,i+1) * MOMZ(k,i+2,j) &
                                        - CNDX(2,i+1) * MOMZ(k,i+1,j) &
                                        + CNDX(3,i+1) * MOMZ(k,i  ,j) &
                                        - CNDX(1,i  ) * MOMZ(k,i-1,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS-1, KE
          num_diff(k,i,j,I_MOMZ,YDIR) = DIFF4 * CDY(j)**4 &
                                      * ( CNDY(1,j+1) * MOMZ(k,i,j+2) &
                                        - CNDY(2,j+1) * MOMZ(k,i,j+1) &
                                        + CNDY(3,j+1) * MOMZ(k,i,j  ) &
                                        - CNDY(1,j  ) * MOMZ(k,i,j-1) )
       enddo
       enddo
       enddo

       ! x-momentum
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          num_diff(k,i,j,I_MOMX,ZDIR) = DIFF4 * CDZ(k)**4 &
                                      * ( CNDZ(1,k+1) * MOMX(k+2,i,j) &
                                        - CNDZ(2,k+1) * MOMX(k+1,i,j) &
                                        + CNDZ(3,k+1) * MOMX(k  ,i,j) &
                                        - CNDZ(1,k  ) * MOMX(k-1,i,j) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          num_diff(k,i,j,I_MOMX,XDIR) = DIFF4 * ( 0.5D0*(CDX(i+1)+CDX(i)) )**4 &
                                      * ( CNMX(1,i  ) * MOMX(k,i+1,j) &
                                        - CNMX(2,i  ) * MOMX(k,i  ,j) &
                                        + CNMX(3,i  ) * MOMX(k,i-1,j) &
                                        - CNMX(1,i-1) * MOMX(k,i-2,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          num_diff(k,i,j,I_MOMX,YDIR) = DIFF4 * CDY(j)**4 &
                                      * ( CNDY(1,j+1) * MOMX(k,i,j+2) &
                                        - CNDY(2,j+1) * MOMX(k,i,j+1) &
                                        + CNDY(3,j+1) * MOMX(k,i,j  ) &
                                        - CNDY(1,j  ) * MOMX(k,i,j-1) )
       enddo
       enddo
       enddo

       ! y-momentum
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS, KE-1
          num_diff(k,i,j,I_MOMY,ZDIR) = DIFF4 * CDZ(k)**4 &
                                      * ( CNDZ(1,k+1) * MOMY(k+2,i,j) &
                                        - CNDZ(2,k+1) * MOMY(k+1,i,j) &
                                        + CNDZ(3,k+1) * MOMY(k  ,i,j) &
                                        - CNDZ(1,k  ) * MOMY(k-1,i,j) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          num_diff(k,i,j,I_MOMY,XDIR) = DIFF4 * CDX(i)**4 &
                                      * ( CNDX(1,i+1) * MOMY(k,i+2,j) &
                                        - CNDX(2,i+1) * MOMY(k,i+1,j) &
                                        + CNDX(3,i+1) * MOMY(k,i  ,j) &
                                        - CNDX(1,i  ) * MOMY(k,i-1,j) )
       enddo
       enddo
       enddo

       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
          num_diff(k,i,j,I_MOMY,YDIR) = DIFF4 * ( 0.5D0*(CDY(j+1)+CDY(j)) )**4 &
                                      * ( CNMY(1,j  ) * MOMY(k,i,j+1) &
                                        - CNMY(2,j  ) * MOMY(k,i,j  ) &
                                        + CNMY(3,j  ) * MOMY(k,i,j-1) &
                                        - CNMY(1,j-1) * MOMY(k,i,j-2) )
       enddo
       enddo
       enddo

       ! rho * theta
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          num_diff(k,i,j,I_RHOT,ZDIR) = DIFF4 * CDZ(k)**4 &
                                      * ( CNDZ(1,k+1) * pott_diff(k+2,i,j)   &
                                        - CNDZ(2,k+1) * pott_diff(k+1,i,j)   &
                                        + CNDZ(3,k+1) * pott_diff(k  ,i,j)   &
                                        - CNDZ(1,k  ) * pott_diff(k-1,i,j) ) &
                                      * 0.5D0 * ( FACT_N * ( DENS(k+1,i,j)+DENS(k  ,i,j) ) &
                                                + FACT_F * ( DENS(k+2,i,j)+DENS(k-1,i,j) ) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE
       do i = IIS,   IIE
          num_diff(KS  ,i,j,I_RHOT,ZDIR) = DIFF2 * CDZ(KS) & 
                                         * 4.0D0 * ( pott_diff(KS+1,i,j)-pott_diff(KS,i,j) ) &
                                         * 0.5D0 * ( DENS(KS+1,i,j)+DENS(KS,i,j) )
          num_diff(KE-1,i,j,I_RHOT,ZDIR) = DIFF2 * CDZ(KE-1) &
                                         * 4.0D0 * ( pott_diff(KE,i,j)-pott_diff(KE-1,i,j) ) &
                                         * 0.5D0 * ( DENS(KE,i,j)+DENS(KE-1,i,j) )
       enddo
       enddo

       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          num_diff(k,i,j,I_RHOT,XDIR) = DIFF4 * CDX(i)**4 &
                                      * ( CNDX(1,i+1) * pott_diff(k,i+2,j)   &
                                        - CNDX(2,i+1) * pott_diff(k,i+1,j)   &
                                        + CNDX(3,i+1) * pott_diff(k,i  ,j)   &
                                        - CNDX(1,i  ) * pott_diff(k,i-1,j) ) &
                                      * 0.5D0 * ( FACT_N * ( DENS(k,i+1,j)+DENS(k,i  ,j) ) &
                                                + FACT_F * ( DENS(k,i+2,j)+DENS(k,i-1,j) ) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          num_diff(k,i,j,I_RHOT,YDIR) = DIFF4 * CDY(j)**4 &
                                      * ( CNDY(1,j+1) * pott_diff(k,i,j+2)   &
                                        - CNDY(2,j+1) * pott_diff(k,i,j+1)   &
                                        + CNDY(3,j+1) * pott_diff(k,i,j  )   &
                                        - CNDY(1,j  ) * pott_diff(k,i,j-1) ) &
                                      * 0.5D0 * ( FACT_N * ( DENS(k,i,j+1)+DENS(k,i,j  ) ) &
                                                + FACT_F * ( DENS(k,i,j+2)+DENS(k,i,j-1) ) )
       enddo
       enddo
       enddo

       ! Gas constant
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
          do k = KS, KE
             QDRY(k,i,j) = 1.D0

             do iw = QQS, QQE
                QDRY(k,i,j) = QDRY(k,i,j) - QTRC(k,i,j,iw)
             enddo
          enddo
          do k = KS, KE
             Rtot(k,i,j) = Rdry!*QDRY(k,i,j) + Rvap*QTRC(k,i,j,I_QV)
          enddo
       enddo
       enddo

    enddo
    enddo ! end tile

#ifdef _FPCOLL_
call STOP_COLLECTION("SET")
call START_COLLECTION("RK3")
#endif

    !##### Start RK #####

    !##### RK1 #####
    rko = 1
    dtrk  = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)
    rdtrk = 1.D0 / dtrk

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS,  KE-1
          VELZ(k,i,j) = 2.D0 * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          VELX(k,i,j) = 2.D0 * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS,   IIE+1
       do k = KS, KE
          VELY(k,i,j) = 2.D0 * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
       enddo
       enddo
       enddo

       ! pressure, pott. temp.
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
          do k = KS, KE
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**CPovCV
          enddo
          do k = KS, KE
             POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) 
          enddo
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          mflx_hi(k,i,j,ZDIR) = MOMZ(k,i,j) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = MOMX(k,i,j) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = MOMY(k,i,j) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update density
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          DENS_RK1(k,i,j) = DENS(k,i,j) &
                          + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                       + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                       + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) ) ! divergence
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                              * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                                + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMZ(k,i+1,j)+MOMZ(k,i  ,j) ) &
                                + FACT_F * ( MOMZ(k,i+2,j)+MOMZ(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMZ(k,i,j+1)+MOMZ(k,i,j  ) ) &
                                + FACT_F * ( MOMZ(k,i,j+2)+MOMZ(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          MOMZ_RK1(k,i,j) = MOMZ(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RFDZ(k) &
                                       + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                     - ( PRES(k+1,i,j)-PRES(k,i,j) ) * RFDZ(k)                         & ! pressure gradient force
                                     - ( DENS(k+1,i,j)+DENS(k,i,j) ) * 0.5D0 * GRAV                    & ! gravity force
                                     + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i+1,j)+VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                                + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                           ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( VELZ(KS  ,i+1,j)+VELZ(KS  ,i,j) ) & ! just above the bottom boundary
                                 * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( VELZ(KE-1,i+1,j)+VELZ(KE-1,i,j) ) & ! just below the top boundary
                                 * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                           ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                              * ( FACT_N * ( MOMX(k,i  ,j)+MOMX(k,i-1,j) ) &
                                + FACT_F * ( MOMX(k,i+1,j)+MOMX(k,i-2,j) ) ) &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMX(k,i,j+1)+MOMX(k,i,j  ) ) &
                                + FACT_F * ( MOMX(k,i,j+2)+MOMX(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMX_RK1(k,i,j) = MOMX(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RFDX(i) &
                                       + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                     - ( PRES(k,i+1,j)-PRES(k,i,j) ) * RFDX(i)                         & ! pressure gradient force
                                     + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i,j+1)+VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                                + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                           ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( VELZ(KS  ,i,j+1)+VELZ(KS  ,i,j) ) & ! just above the bottom boundary
                                 * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) )              &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( VELZ(KE-1,i,j+1)+VELZ(KE-1,i,j) ) & ! just below the top boundary
                                 * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) )              &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                           ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMY(k,i+1,j)+MOMY(k,i  ,j) ) &
                                + FACT_F * ( MOMY(k,i+2,j)+MOMY(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                              * ( FACT_N * ( MOMY(k,i,j  )+MOMY(k,i,j-1) ) &
                                + FACT_F * ( MOMY(k,i,j+1)+MOMY(k,i,j-2) ) ) &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMY_RK1(k,i,j) = MOMY(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k)   &
                                       + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i)   &
                                       + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) & ! flux divergence
                                     - ( PRES(k,i,j+1)-PRES(k,i,j) ) * RFDY(j)                         & ! pressure gradient force
                                     + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j)+POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j)+POTT(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                   ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR)  &      ! just above the bottom boundary
                                 * ( POTT(KS+1,i,j)+POTT(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR)  &      ! just below the top boundary
                                 * ( POTT(KE,i,j)+POTT(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                   ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR) &
                              * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT_RK1(k,i,j) = RHOT(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) & ! divergence
                                     + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping
       enddo
       enddo
       enddo

    enddo
    enddo

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
    rko = 2
    dtrk  = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)
    rdtrk = 1.D0 / dtrk

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS,  KE-1
          VELZ(k,i,j) = 2.D0 * MOMZ_RK1(k,i,j) / ( DENS_RK1(k+1,i,j)+DENS_RK1(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          VELX(k,i,j) = 2.D0 * MOMX_RK1(k,i,j) / ( DENS_RK1(k,i+1,j)+DENS_RK1(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS,   IIE+1
       do k = KS, KE
          VELY(k,i,j) = 2.D0 * MOMY_RK1(k,i,j) / ( DENS_RK1(k,i,j+1)+DENS_RK1(k,i,j) )
       enddo
       enddo
       enddo

       ! pressure, pott. temp.
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
          do k = KS, KE
             PRES(k,i,j) = P00 * ( RHOT_RK1(k,i,j) * Rtot(k,i,j) / P00 )**CPovCV
          enddo
          do k = KS, KE
             POTT(k,i,j) = RHOT_RK1(k,i,j) / DENS_RK1(k,i,j) 
          enddo
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          mflx_hi(k,i,j,ZDIR) = MOMZ_RK1(k,i,j) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = MOMX_RK1(k,i,j) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = MOMY_RK1(k,i,j) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update density
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          DENS_RK2(k,i,j) = DENS(k,i,j) &
                          + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                       + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                       + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) ) ! divergence
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                              * ( FACT_N * ( MOMZ_RK1(k  ,i,j)+MOMZ_RK1(k-1,i,j) ) &
                                + FACT_F * ( MOMZ_RK1(k+1,i,j)+MOMZ_RK1(k-2,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMZ_RK1(k,i+1,j)+MOMZ_RK1(k,i  ,j) ) &
                                + FACT_F * ( MOMZ_RK1(k,i+2,j)+MOMZ_RK1(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMZ_RK1(k,i,j+1)+MOMZ_RK1(k,i,j  ) ) &
                                + FACT_F * ( MOMZ_RK1(k,i,j+2)+MOMZ_RK1(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          MOMZ_RK2(k,i,j) = MOMZ(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RFDZ(k)   &
                                       + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RCDX(i)   &
                                       + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                     - ( PRES(k+1,i,j)-PRES(k,i,j) ) * RFDZ(k)                         & ! pressure gradient force
                                     - ( DENS_RK1(k+1,i,j)+DENS_RK1(k,i,j) ) * 0.5D0 * GRAV            & ! gravity force
                                     + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i+1,j)+VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMX_RK1(k+1,i,j)+MOMX_RK1(k  ,i,j) ) &
                                + FACT_F * ( MOMX_RK1(k+2,i,j)+MOMX_RK1(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                           ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( VELZ(KS  ,i+1,j)+VELZ(KS  ,i,j) ) & ! just above the bottom boundary
                                 * ( MOMX_RK1(KS+1,i,j)+MOMX_RK1(KS,i,j) )      &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( VELZ(KE-1,i+1,j)+VELZ(KE-1,i,j) ) & ! just below the top boundary
                                 * ( MOMX_RK1(KE,i,j)+MOMX_RK1(KE-1,i,j) )      &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                           ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                              * ( FACT_N * ( MOMX_RK1(k,i  ,j)+MOMX_RK1(k,i-1,j) ) &
                                + FACT_F * ( MOMX_RK1(k,i+1,j)+MOMX_RK1(k,i-2,j) ) ) &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMX_RK1(k,i,j+1)+MOMX_RK1(k,i,j  ) ) &
                                + FACT_F * ( MOMX_RK1(k,i,j+2)+MOMX_RK1(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMX_RK2(k,i,j) = MOMX(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RFDX(i) &
                                       + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                     - ( PRES(k,i+1,j)-PRES(k,i,j) ) * RFDX(i)                         & ! pressure gradient force
                                     + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i,j+1)+VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMY_RK1(k+1,i,j)+MOMY_RK1(k  ,i,j) ) &
                                + FACT_F * ( MOMY_RK1(k+2,i,j)+MOMY_RK1(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( VELZ(KS  ,i,j+1)+VELZ(KS  ,i,j) ) & ! just above the bottom boundary
                                 * ( MOMY_RK1(KS+1,i,j)+MOMY_RK1(KS,i,j) )                      &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( VELZ(KE-1,i,j+1)+VELZ(KE-1,i,j) ) & ! just below the top boundary
                                 * ( MOMY_RK1(KE,i,j)+MOMY_RK1(KE-1,i,j) )                      &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMY_RK1(k,i+1,j)+MOMY_RK1(k,i  ,j) )     &
                                + FACT_F * ( MOMY_RK1(k,i+2,j)+MOMY_RK1(k,i-1,j) ) )   &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                              * ( FACT_N * ( MOMY_RK1(k,i,j  )+MOMY_RK1(k,i,j-1) )     &
                                + FACT_F * ( MOMY_RK1(k,i,j+1)+MOMY_RK1(k,i,j-2) ) )   &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMY_RK2(k,i,j) = MOMY(k,i,j) &
                                + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k)   &
                                             + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i)   &
                                             + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) & ! flux divergence
                                           - ( PRES(k,i,j+1)-PRES(k,i,j) ) * RFDY(j)     & ! pressure gradient force
                                           + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j)+POTT(k  ,i,j) )   &
                                + FACT_F * ( POTT(k+2,i,j)+POTT(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                   ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR)  &      ! just above the bottom boundary
                                 * ( POTT(KS+1,i,j)+POTT(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR)  &      ! just below the top boundary
                                 * ( POTT(KE,i,j)+POTT(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                   ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) )   &
                                + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR) &
                              * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) )   &
                                + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT_RK2(k,i,j) = RHOT(k,i,j) &
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) & ! divergence
                                     + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping
       enddo
       enddo
       enddo

    enddo
    enddo

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
    rko = 3
    dtrk  = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)
    rdtrk = 1.D0 / dtrk

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS,  KE-1
          VELZ(k,i,j) = 2.D0 * MOMZ_RK2(k,i,j) / ( DENS_RK2(k+1,i,j)+DENS_RK2(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          VELX(k,i,j) = 2.D0 * MOMX_RK2(k,i,j) / ( DENS_RK2(k,i+1,j)+DENS_RK2(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS,   IIE+1
       do k = KS, KE
          VELY(k,i,j) = 2.D0 * MOMY_RK2(k,i,j) / ( DENS_RK2(k,i,j+1)+DENS_RK2(k,i,j) )
       enddo
       enddo
       enddo

       ! pressure, pott. temp.
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
          do k = KS, KE
             PRES(k,i,j) = P00 * ( RHOT_RK2(k,i,j) * Rtot(k,i,j) / P00 )**CPovCV
          enddo
          do k = KS, KE
             POTT(k,i,j) = RHOT_RK2(k,i,j) / DENS_RK2(k,i,j) 
          enddo
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          mflx_hi(k,i,j,ZDIR) = MOMZ_RK2(k,i,j) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = MOMX_RK2(k,i,j) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = MOMY_RK2(k,i,j) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update density
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          DENS(k,i,j) = DENS(k,i,j) &
                            + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k)   &
                                         + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i)   &
                                         + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) ) ! divergence
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                              * ( FACT_N * ( MOMZ_RK2(k  ,i,j)+MOMZ_RK2(k-1,i,j) )     &
                                + FACT_F * ( MOMZ_RK2(k+1,i,j)+MOMZ_RK2(k-2,i,j) ) )   &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMZ_RK2(k,i+1,j)+MOMZ_RK2(k,i  ,j) )     &
                                + FACT_F * ( MOMZ_RK2(k,i+2,j)+MOMZ_RK2(k,i-1,j) ) )   &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMZ_RK2(k,i,j+1)+MOMZ_RK2(k,i,j  ) )     &
                                + FACT_F * ( MOMZ_RK2(k,i,j+2)+MOMZ_RK2(k,i,j-1) ) )   &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          MOMZ(k,i,j) = MOMZ(k,i,j) &
                      + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RFDZ(k)   &
                                   + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RCDX(i)   &
                                   + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                 - ( PRES(k+1,i,j)-PRES(k,i,j) ) * RFDZ(k)                         & ! pressure gradient force
                                 - ( DENS_RK2(k+1,i,j)+DENS_RK2(k,i,j) ) * 0.5D0 * GRAV            & ! gravity force
                                 + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i+1,j)+VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMX_RK2(k+1,i,j)+MOMX_RK2(k  ,i,j) )     &
                                + FACT_F * ( MOMX_RK2(k+2,i,j)+MOMX_RK2(k-1,i,j) ) )   &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( VELZ(KS  ,i+1,j)+VELZ(KS  ,i,j) ) & ! just above the bottom boundary
                                 * ( MOMX_RK2(KS+1,i,j)+MOMX_RK2(KS,i,j) )                      &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( VELZ(KE-1,i+1,j)+VELZ(KE-1,i,j) ) & ! just below the top boundary
                                 * ( MOMX_RK2(KE,i,j)+MOMX_RK2(KE-1,i,j) )                      &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                              * ( FACT_N * ( MOMX_RK2(k,i  ,j)+MOMX_RK2(k,i-1,j) )     &
                                + FACT_F * ( MOMX_RK2(k,i+1,j)+MOMX_RK2(k,i-2,j) ) )   &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMX_RK2(k,i,j+1)+MOMX_RK2(k,i,j  ) )     &
                                + FACT_F * ( MOMX_RK2(k,i,j+2)+MOMX_RK2(k,i,j-1) ) )   &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMX(k,i,j) = MOMX(k,i,j) &
                            + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RCDZ(k)   &
                                         + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RFDX(i)   &
                                         + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                       - ( PRES(k,i+1,j)-PRES(k,i,j) ) * RFDX(i)     & ! pressure gradient force
                                       + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( VELZ(k,i,j+1)+VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMY_RK2(k+1,i,j)+MOMY_RK2(k  ,i,j) )     &
                                + FACT_F * ( MOMY_RK2(k+2,i,j)+MOMY_RK2(k-1,i,j) ) )   &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( VELZ(KS  ,i,j+1)+VELZ(KS  ,i,j) ) & ! just above the bottom boundary
                                 * ( MOMY_RK2(KS+1,i,j)+MOMY_RK2(KS,i,j) )                      &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( VELZ(KE-1,i,j+1)+VELZ(KE-1,i,j) ) & ! just below the top boundary
                                 * ( MOMY_RK2(KE,i,j)+MOMY_RK2(KE-1,i,j) )                      &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMY_RK2(k,i+1,j)+MOMY_RK2(k,i  ,j) ) &
                                + FACT_F * ( MOMY_RK2(k,i+2,j)+MOMY_RK2(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                              * ( FACT_N * ( MOMY_RK2(k,i,j  )+MOMY_RK2(k,i,j-1) ) &
                                + FACT_F * ( MOMY_RK2(k,i,j+1)+MOMY_RK2(k,i,j-2) ) ) &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMY(k,i,j) = MOMY(k,i,j) &
                      + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k) &
                                   + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i) &
                                   + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) & ! flux divergence
                                 - ( PRES(k,i,j+1)-PRES(k,i,j) ) * RFDY(j)                         & ! pressure gradient force
                                 + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j)+POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j)+POTT(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                   ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR)  &      ! just above the bottom boundary
                                 * ( POTT(KS+1,i,j)+POTT(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR)  &      ! just below the top boundary
                                 * ( POTT(KE,i,j)+POTT(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                   ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR) &
                              * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT(k,i,j) = RHOT(k,i,j) &
                      + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                   + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                   + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) & ! divergence
                                 + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping
       enddo
       enddo
       enddo

    enddo
    enddo

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

#ifdef _FPCOLL_
call STOP_COLLECTION("RK3")
call START_COLLECTION("FCT")
#endif

    !##### advection of scalar quantity #####

    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
          qflx_lo(k,i,j,ZDIR) = 0.5D0 * (     mflx_hi(k,i,j,ZDIR)  * ( QTRC(k+1,i,j,iq)+QTRC(k,i,j,iq) ) &
                                        - abs(mflx_hi(k,i,j,ZDIR)) * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) )

          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                      &
                              * ( FACT_N * ( QTRC(k+1,i,j,iq)+QTRC(k  ,i,j,iq) ) &
                                + FACT_F * ( QTRC(k+2,i,j,iq)+QTRC(k-1,i,j,iq) ) )

          qflx_anti(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_lo(k,i,j,ZDIR)
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_lo(KS-1,i,j,ZDIR) = 0.D0                                                                          ! bottom boundary
          qflx_lo(KS  ,i,j,ZDIR) = 0.5D0 * (     mflx_hi(KS  ,i,j,ZDIR)  * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) ) & ! just above the bottom boundary
                                           - abs(mflx_hi(KS  ,i,j,ZDIR)) * ( QTRC(KS+1,i,j,iq)-QTRC(KS,i,j,iq) ) )
          qflx_lo(KE-1,i,j,ZDIR) = 0.5D0 * (     mflx_hi(KE-1,i,j,ZDIR)  * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) ) & ! just below the top boundary
                                           - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( QTRC(KE,i,j,iq)-QTRC(KE-1,i,j,iq) ) )
          qflx_lo(KE  ,i,j,ZDIR) = 0.D0                                                                          ! top boundary

          qflx_hi(KS-1,i,j,ZDIR) = 0.D0 
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR) * ( QTRC(KS+1,i,j,iq)+QTRC(KS,i,j,iq) )
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR) * ( QTRC(KE,i,j,iq)+QTRC(KE-1,i,j,iq) )
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0 

          qflx_anti(KS-1,i,j,ZDIR) = 0.D0
          qflx_anti(KS  ,i,j,ZDIR) = qflx_hi(KS  ,i,j,ZDIR) - qflx_lo(KS  ,i,j,ZDIR)
          qflx_anti(KE-1,i,j,ZDIR) = qflx_hi(KE-1,i,j,ZDIR) - qflx_lo(KE-1,i,j,ZDIR)
          qflx_anti(KE  ,i,j,ZDIR) = 0.D0
       enddo
       enddo
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_lo(k,i,j,XDIR) = 0.5D0 * (     mflx_hi(k,i,j,XDIR)  * ( QTRC(k,i+1,j,iq)+QTRC(k,i,j,iq) ) &
                                        - abs(mflx_hi(k,i,j,XDIR)) * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) )

          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                      &
                              * ( FACT_N * ( QTRC(k,i+1,j,iq)+QTRC(k,i  ,j,iq) ) &
                                + FACT_F * ( QTRC(k,i+2,j,iq)+QTRC(k,i-1,j,iq) ) )

          qflx_anti(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_lo(k,i,j,XDIR)
       enddo
       enddo
       enddo
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_lo(k,i,j,YDIR) = 0.5D0 * (     mflx_hi(k,i,j,YDIR)  * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j,iq) ) &
                                        - abs(mflx_hi(k,i,j,YDIR)) * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) )

          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                      &
                              * ( FACT_N * ( QTRC(k,i,j+1,iq)+QTRC(k,i,j  ,iq) ) &
                                + FACT_F * ( QTRC(k,i,j+2,iq)+QTRC(k,i,j-1,iq) ) )

          qflx_anti(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) - qflx_lo(k,i,j,YDIR)
       enddo
       enddo
       enddo

       ! --- STEP C: compute the outgoing fluxes in each cell ---
       !OCL SIMD
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE

          pjmns = max( 0.D0, qflx_hi(k,i,j,ZDIR) ) - min( 0.D0, qflx_hi(k-1,i  ,j  ,ZDIR) ) &
                + max( 0.D0, qflx_hi(k,i,j,XDIR) ) - min( 0.D0, qflx_hi(k  ,i-1,j  ,XDIR) ) &
                + max( 0.D0, qflx_hi(k,i,j,YDIR) ) - min( 0.D0, qflx_hi(k  ,i  ,j-1,YDIR) )

          if ( pjmns > 0 ) then
             tmp = QTRC(k,i,j,iq) / pjmns * rdtrk
             rjmns(k,i,j,ZDIR) = tmp * CDZ(k)
             rjmns(k,i,j,XDIR) = tmp * CDX(i)
             rjmns(k,i,j,YDIR) = tmp * CDY(j)
          else
             rjmns(k,i,j,ZDIR) = 0.D0
             rjmns(k,i,j,XDIR) = 0.D0
             rjmns(k,i,j,YDIR) = 0.D0
          endif

       enddo
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( rjmns(:,:,:,ZDIR), ZDIR )
    call COMM_vars8( rjmns(:,:,:,XDIR), XDIR )
    call COMM_vars8( rjmns(:,:,:,YDIR), YDIR )
    call COMM_wait ( rjmns(:,:,:,ZDIR), ZDIR )
    call COMM_wait ( rjmns(:,:,:,XDIR), XDIR )
    call COMM_wait ( rjmns(:,:,:,YDIR), YDIR )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! --- [STEP 7S] limit the antidiffusive flux ---
       !OCL norecurrence
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS-1, KE
          if ( qflx_anti(k,i,j,ZDIR) >= 0 ) then
             if ( rjmns(k  ,i,j,ZDIR) >= 0.D0 .AND. rjmns(k  ,i,j,ZDIR) < 1.D0 ) then
                qflx_hi(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) * rjmns(k  ,i,j,ZDIR)
             endif
          else
             if ( rjmns(k+1,i,j,ZDIR) >= 0.D0 .AND. rjmns(k+1,i,j,ZDIR) < 1.D0 ) then
                qflx_hi(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) * rjmns(k+1,i,j,ZDIR)
             endif
          endif
       enddo
       enddo
       enddo
       !OCL norecurrence
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          if ( qflx_anti(k,i,j,XDIR) >= 0 ) then
             if ( rjmns(k,i  ,j,XDIR) >= 0.D0 .AND. rjmns(k,i  ,j,XDIR) < 1.D0 ) then
                qflx_hi(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) * rjmns(k,i  ,j,XDIR)
             endif
          else
             if ( rjmns(k,i+1,j,XDIR) >= 0.D0 .AND. rjmns(k,i+1,j,XDIR) < 1.D0 ) then
                qflx_hi(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) * rjmns(k,i+1,j,XDIR)
             endif
          endif
       enddo
       enddo
       enddo
       !OCL norecurrence
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          if ( qflx_anti(k,i,j,YDIR) >= 0 ) then
             if ( rjmns(k,i,j  ,YDIR) >= 0.D0 .AND. rjmns(k,i,j  ,YDIR) < 1.D0 ) then
                qflx_hi(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) * rjmns(k,i,j  ,YDIR)
             endif
          else
             if ( rjmns(k,i,j+1,YDIR) >= 0.D0 .AND. rjmns(k,i,j+1,YDIR) < 1.D0 ) then
                qflx_hi(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) * rjmns(k,i,j+1,YDIR)
             endif
          endif
       enddo
       enddo
       enddo

       !--- modify value with antidiffusive fluxes
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          QTRC(k,i,j,iq) = ( QTRC(k,i,j,iq) * dens_s(k,i,j) &
                           + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                        + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                        + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) ) &
                           ) / DENS(k,i,j)
       enddo
       enddo
       enddo

    enddo
    enddo

    enddo ! scalar quantities loop

    do iq = 1, QA
       call COMM_vars8( QTRC(:,:,:,iq), iq )
    enddo
    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), iq )
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("FCT")
#endif

    enddo ! dynamical steps

#ifdef _FPCOLL_
call STOP_COLLECTION("DYNAMICS")
#endif

    call ATMOS_vars_total

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
