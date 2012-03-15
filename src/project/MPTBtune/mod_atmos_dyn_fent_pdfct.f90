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
!! @li      2011-01-04 (H.Yashiro)  [mod] Nonblocking communication (Y.Ohno)
!! @li      2011-01-25 (H.Yashiro)  [mod] Bugfix (Y.Miyamoto)
!! @li      2011-01-25 (H.Yashiro)  [mod] Positive definite FCT (Y.Miyamoto)
!! @li      2011-02-14 (H.Yashiro)  [mod] Expand RK3, put directive (fasttune2)
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
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
  integer, parameter :: KMAX = 336
  integer, parameter :: IMAX =  64
  integer, parameter :: JMAX =  64
  integer, parameter :: KA   = 340
  integer, parameter :: IA   =  68
  integer, parameter :: JA   =  68
  integer, parameter :: KS   =   3
  integer, parameter :: KE   = 338
  integer, parameter :: IS   =   3
  integer, parameter :: IE   =  66
  integer, parameter :: JS   =   3
  integer, parameter :: JE   =  66
  integer, parameter :: VA   = 16
  integer, parameter :: QA   = 11
  integer, parameter :: I_DENS = 1
  integer, parameter :: I_MOMX = 2
  integer, parameter :: I_MOMY = 3
  integer, parameter :: I_MOMZ = 4
  integer, parameter :: I_RHOT = 5

  integer, parameter :: IBLOCK = 8
  integer, parameter :: JBLOCK = 8

  integer, parameter :: I_PRES = 1
  integer, parameter :: I_VELX = 2
  integer, parameter :: I_VELY = 3
  integer, parameter :: I_VELZ = 4
  integer, parameter :: I_POTT = 5

  integer, parameter :: ZDIR   = 1
  integer, parameter :: XDIR   = 2
  integer, parameter :: YDIR   = 3

  ! time settings
  integer, parameter :: RK = 3 ! order of Runge-Kutta scheme

  ! advection settings
  real(8), parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second

  ! numerical filter settings
  real(8), save      :: ATMOS_DYN_numerical_diff = 1.D-2 ! nondimensional numerical diffusion
  real(8), save      :: DIFF4 ! for 4th order numerical filter
  real(8), save      :: DIFF2 ! for 2nd order numerical filter

  real(8), allocatable, save :: CNDZ(:,:)
  real(8), allocatable, save :: CNMZ(:,:)
  real(8), allocatable, save :: CNDX(:,:)
  real(8), allocatable, save :: CNMX(:,:)
  real(8), allocatable, save :: CNDY(:,:)
  real(8), allocatable, save :: CNMY(:,:)

!  integer, private, save :: IBLOCK, JBLOCK

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only : &
       CDX  => GRID_CDX, &
       CDY  => GRID_CDY, &
       CDZ  => GRID_CDZ
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff
!       IBLOCK, &
!       JBLOCK

    integer :: ierr
    integer :: i, j, k
    !---------------------------------------------------------------------------

!    IBLOCK = IMAX
!    JBLOCK = JMAX

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
!    if    ( mod(IMAX,IBLOCK) > 0 ) then
!       if( IO_L ) write(IO_FID_LOG,*) 'xxx number of grid size IMAX must be divisible by IBLOCK! ', IMAX, IBLOCK
!       call PRC_MPIstop
!    elseif( mod(JMAX,JBLOCK) > 0 ) then
!       if( IO_L ) write(IO_FID_LOG,*) 'xxx number of grid size JMAX must be divisible by JBLOCK! ', JMAX, JBLOCK
!       call PRC_MPIstop
!    endif

    DIFF4 = - ATMOS_DYN_numerical_diff * (-1.D0)**dble( 4/2+1 )
    DIFF2 = - ATMOS_DYN_numerical_diff * (-1.D0)**dble( 2/2+1 )

    allocate( CNDX(3,IA) )
    allocate( CNMX(3,IA) )
    allocate( CNDY(3,JA) )
    allocate( CNMY(3,JA) )
    allocate( CNDZ(3,KA) )
    allocate( CNMZ(3,KA) )

    ! z djrectjon
    do k = KS-1, KE+1
       CNDZ(1,k) = 1.D0 / ( (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 )
       CNDZ(2,k) = 1.D0 / ( (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 )
    enddo
    CNDZ(1,1)  = CNDZ(1,KS-1)
    CNDZ(2,1)  = CNDZ(2,KS-1)
    CNDZ(1,KA) = CNDZ(1,KE+1)
    CNDZ(2,KA) = CNDZ(2,KE+1)

    do k = KS, KE+2
       CNDZ(3,k) = 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k-1)+CDZ(k-2)) * 0.5D0 )
    enddo
    CNDZ(3,1)    = CNDZ(3,KS)
    CNDZ(3,KS-1) = CNDZ(3,KS)

    do k = KS-2, KE+1
       CNMZ(1,k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) )
    enddo
    CNMZ(1,KA) = CNMZ(1,KE+1)

    do k = KS-1, KE+1
       CNMZ(2,k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &   
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &  
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) 
       CNMZ(3,k) = 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) &
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) )
    enddo
    CNMZ(2,1)  = CNMZ(2,KS-1)
    CNMZ(3,1)  = CNMZ(3,KS-1)
    CNMZ(2,KA) = CNMZ(2,KE+1)
    CNMZ(3,KA) = CNMZ(3,KE+1)

    ! x direction
    do i = IS-1, IE+1
       CNDX(1,i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
       CNDX(2,i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
    enddo
    CNDX(1,1)  = CNDX(1,IS-1)
    CNDX(2,1)  = CNDX(2,IS-1)
    CNDX(1,IA) = CNDX(1,IE+1)
    CNDX(2,IA) = CNDX(2,IE+1)

    do i = IS, IE+2
       CNDX(3,i) = 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i-1)+CDX(i-2)) * 0.5D0 )
    enddo
    CNDX(3,1)    = CNDX(3,IS)
    CNDX(3,IS-1) = CNDX(3,IS)

    do i = IS-2, IE+1
       CNMX(1,i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) )
    enddo
    CNMX(1,IA) = CNMX(1,IE+1)

    do i = IS-1, IE+1
       CNMX(2,i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) )
       CNMX(3,i) = 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) )
    enddo
    CNMX(2,1)  = CNMX(2,IS-1)
    CNMX(3,1)  = CNMX(3,IS-1)
    CNMX(2,IA) = CNMX(2,IE+1)
    CNMX(3,IA) = CNMX(3,IE+1)

    ! y direction
    do j = JS-1, JE+1
       CNDY(1,j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
       CNDY(2,j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
    enddo
    CNDY(1,1)  = CNDY(1,JS-1)
    CNDY(2,1)  = CNDY(2,JS-1)
    CNDY(1,JA) = CNDY(1,JE+1)
    CNDY(2,JA) = CNDY(2,JE+1)

    do j = JS, JE+2
       CNDY(3,j) = 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j-1)+CDY(j-2)) * 0.5D0 )
    enddo
    CNDY(3,1)    = CNDY(3,JS)
    CNDY(3,JS-1) = CNDY(3,JS)

    do j = JS-2, JE+1
       CNMY(1,j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) )
    enddo
    CNMY(1,JA) = CNMY(1,JE+1)

    do j = JS-1, JE+1
       CNMY(2,j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &   
                 + 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &  
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) 
       CNMY(3,j) = 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) &
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) )
    enddo
    CNMY(2,1)  = CNMY(2,JS-1)
    CNMY(3,1)  = CNMY(3,JS-1)
    CNMY(2,JA) = CNMY(2,JE+1)
    CNMY(3,JA) = CNMY(3,JE+1)

  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd
    use mod_time, only: &
       TIME_DTSEC_ATMOS_DYN, &
       TIME_NSTEP_ATMOS_DYN
    use mod_comm, only: &
       COMM_vars, &
       COMM_vars8, &
       COMM_wait, &
       COMM_total
    use mod_grid, only : &
       CDX  => GRID_CDX,  &
       CDY  => GRID_CDY,  &
       CDZ  => GRID_CDZ,  &
       RDXC => GRID_RCDX, &
       RDYC => GRID_RCDY, &
       RDZC => GRID_RCDZ, &
       RDXF => GRID_RFDX, &
       RDYF => GRID_RFDY, &
       RDZF => GRID_RFDZ
    use mod_atmos_vars, only: &
       var => atmos_var, &
       A_NAME
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

    ! work
    real(8) :: var_RK1  (KA,IA,JA,5)   ! prognostic variables (+1/3 step)
    real(8) :: var_RK2  (KA,IA,JA,5)   ! prognostic variables (+2/3 step)
    real(8) :: diagvar  (KA,IA,JA,5)   ! diagnostic variables (work)

    ! rayleigh damping, numerical diffusion
    real(8) :: dens_s   (KA,IA,JA)     ! saved density
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
    real(8) :: rjmns    (KA,IA,JA,3)   ! minus in (x,y,z)-direction
    real(8) :: pjmns, tmp

    integer :: IIS, IIE
    integer :: JJS, JJE

    real(8) :: dtrk, rdtrk
    integer :: i, j, k, iq, rko, step
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("DYNAMICS")
#endif

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

    do step = 1, TIME_NSTEP_ATMOS_DYN

!    diagvar  (:,:,:,:)   = -9.999D30
!    dens_diff(:,:,:)     = -9.999D30
!    pott_diff(:,:,:)     = -9.999D30
!    ray_damp (:,:,:,:)   = -9.999D30
!    num_diff (:,:,:,:,:) = -9.999D30
!    mflx_hi  (:,:,:,:)   = -9.999D30
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
                                   * ( var(k,i,j,I_MOMZ)          &
                                     - DAMP_var(k,i,j,I_BND_VELZ) &
                                     * 0.5D0 * ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) ) )
       enddo
       do k = KS, KE
          ray_damp(k,i,j,I_MOMX) = - DAMP_alpha(k,i,j,I_BND_VELX) &
                              * ( var(k,i,j,I_MOMX) &
                                - DAMP_var(k,i,j,I_BND_VELX) * 0.5D0 * ( var(k,i+1,j,I_DENS)+var(k,i,j,I_DENS) ) )
          ray_damp(k,i,j,I_MOMY) = - DAMP_alpha(k,i,j,I_BND_VELY) &
                              * ( var(k,i,j,I_MOMY) &
                                - DAMP_var(k,i,j,I_BND_VELY) * 0.5D0 * ( var(k,i,j+1,I_DENS)+var(k,i,j,I_DENS) ) )
          ray_damp(k,i,j,I_RHOT) = - DAMP_alpha(k,i,j,I_BND_POTT) &
                              * ( var(k,i,j,I_RHOT) &
                                - DAMP_var(k,i,j,I_BND_POTT) * var(k,i,j,I_DENS) )
       enddo 
    enddo
    enddo

    !--- prepare numerical diffusion coefficient
    do j  = JJS-2, JJE+2
    do i  = IIS-2, IIE+2
       do k = KS, KE
          dens_s(k,i,j)    = var(k,i,j,I_DENS)
          dens_diff(k,i,j) = var(k,i,j,I_DENS)                     - REF_dens(k)
          pott_diff(k,i,j) = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS) - REF_pott(k)
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
       num_diff(k,i,j,I_DENS,ZDIR) = DIFF4 * CDZ(k)**4                  &
                                   * ( CNDZ(1,k+1) * dens_diff(k+2,i,j) &
                                     - CNDZ(2,k+1) * dens_diff(k+1,i,j) &
                                     + CNDZ(3,k+1) * dens_diff(k  ,i,j) &
                                     - CNDZ(1,k  ) * dens_diff(k-1,i,j) )

       num_diff(k,i,j,I_RHOT,ZDIR) = DIFF4 * CDZ(k)**4                    &
                                   * ( CNDZ(1,k+1) * pott_diff(k+2,i,j)   &
                                     - CNDZ(2,k+1) * pott_diff(k+1,i,j)   &
                                     + CNDZ(3,k+1) * pott_diff(k  ,i,j)   &
                                     - CNDZ(1,k  ) * pott_diff(k-1,i,j) ) &
                                   * 0.5D0 * ( FACT_N * ( var(k+1,i,j,I_DENS)+var(k  ,i,j,I_DENS) ) &
                                             + FACT_F * ( var(k+2,i,j,I_DENS)+var(k-1,i,j,I_DENS) ) )
    enddo
    enddo
    enddo

    do j = JJS,   JJE
    do i = IIS,   IIE
       num_diff(KS  ,i,j,I_DENS,ZDIR) = DIFF2 * CDZ(KS)                                 &
                                      * 4.0D0 * ( dens_diff(KS+1,i,j)-dens_diff(KS,i,j) ) 

       num_diff(KS  ,i,j,I_RHOT,ZDIR) = DIFF2 * CDZ(KS)                                   & 
                                      * 4.0D0 * ( pott_diff(KS+1,i,j)-pott_diff(KS,i,j) ) &
                                      * 0.5D0 * ( var(KS+1,i,j,I_DENS)+var(KS,i,j,I_DENS) )

       num_diff(KE-1,i,j,I_DENS,ZDIR) = DIFF2 * CDZ(KE-1)                               &
                                      * 4.0D0 * ( dens_diff(KE,i,j)-dens_diff(KE-1,i,j) )

       num_diff(KE-1,i,j,I_RHOT,ZDIR) = DIFF2 * CDZ(KE-1)                                 &
                                      * 4.0D0 * ( pott_diff(KE,i,j)-pott_diff(KE-1,i,j) ) &
                                      * 0.5D0 * ( var(KE,i,j,I_DENS)+var(KE-1,i,j,I_DENS) )
    enddo
    enddo

    do j = JJS,   JJE
    do i = IIS-1, IIE
    do k = KS,   KE
       num_diff(k,i,j,I_DENS,XDIR) = DIFF4 * CDX(i)**4                  &
                                   * ( CNDX(1,i+1) * dens_diff(k,i+2,j) &
                                     - CNDX(2,i+1) * dens_diff(k,i+1,j) &
                                     + CNDX(3,i+1) * dens_diff(k,i  ,j) &
                                     - CNDX(1,i  ) * dens_diff(k,i-1,j) )

       num_diff(k,i,j,I_RHOT,XDIR) = DIFF4 * CDX(i)**4                    &
                                   * ( CNDX(1,i+1) * pott_diff(k,i+2,j)   &
                                     - CNDX(2,i+1) * pott_diff(k,i+1,j)   &
                                     + CNDX(3,i+1) * pott_diff(k,i  ,j)   &
                                     - CNDX(1,i  ) * pott_diff(k,i-1,j) ) &
                                   * 0.5D0 * ( FACT_N * ( var(k,i+1,j,I_DENS)+var(k,i  ,j,I_DENS) ) &
                                             + FACT_F * ( var(k,i+2,j,I_DENS)+var(k,i-1,j,I_DENS) ) )
    enddo
    enddo
    enddo

    !OCL SERIAL
    do j = JJS-1, JJE
    !OCL PARALLEL
    do i = IIS,   IIE
    do k = KS,   KE
       num_diff(k,i,j,I_DENS,YDIR) = DIFF4 * CDY(j)**4                  &
                                   * ( CNDY(1,j+1) * dens_diff(k,i,j+2) &
                                     - CNDY(2,j+1) * dens_diff(k,i,j+1) &
                                     + CNDY(3,j+1) * dens_diff(k,i,j  ) &
                                     - CNDY(1,j  ) * dens_diff(k,i,j-1) )

       num_diff(k,i,j,I_RHOT,YDIR) = DIFF4 * CDY(j)**4                    &
                                   * ( CNDY(1,j+1) * pott_diff(k,i,j+2)   &
                                     - CNDY(2,j+1) * pott_diff(k,i,j+1)   &
                                     + CNDY(3,j+1) * pott_diff(k,i,j  )   &
                                     - CNDY(1,j  ) * pott_diff(k,i,j-1) ) &
                                   * 0.5D0 * ( FACT_N * ( var(k,i,j+1,I_DENS)+var(k,i,j  ,I_DENS) ) &
                                             + FACT_F * ( var(k,i,j+2,I_DENS)+var(k,i,j-1,I_DENS) ) )
    enddo
    enddo
    enddo

    ! z-momentum
    !OCL PARALLEL
    do j = JJS,   JJE
    do i = IIS,   IIE
    do k = KS,   KE
       num_diff(k,i,j,I_MOMZ,ZDIR) = DIFF4 * ( 0.5D0*(CDZ(k+1)+CDZ(k)) )**4 &
                                   * ( CNMZ(1,k  ) * var(k+1,i,j,I_MOMZ) &
                                     - CNMZ(2,k  ) * var(k  ,i,j,I_MOMZ) &
                                     + CNMZ(3,k  ) * var(k-1,i,j,I_MOMZ) &
                                     - CNMZ(1,k-1) * var(k-2,i,j,I_MOMZ) )
    enddo
    enddo
    enddo
    do j = JJS,   JJE
    do i = IIS-1, IIE
    do k = KS-1, KE
       num_diff(k,i,j,I_MOMZ,XDIR) = DIFF4 * CDX(i)**4 &
                                   * ( CNDX(1,i+1) * var(k,i+2,j,I_MOMZ) &
                                     - CNDX(2,i+1) * var(k,i+1,j,I_MOMZ) &
                                     + CNDX(3,i+1) * var(k,i  ,j,I_MOMZ) &
                                     - CNDX(1,i  ) * var(k,i-1,j,I_MOMZ) )
    enddo
    enddo
    enddo
    do j = JJS-1, JJE
    do i = IIS,   IIE
    do k = KS-1, KE
       num_diff(k,i,j,I_MOMZ,YDIR) = DIFF4 * CDY(j)**4 &
                                   * ( CNDY(1,j+1) * var(k,i,j+2,I_MOMZ) &
                                     - CNDY(2,j+1) * var(k,i,j+1,I_MOMZ) &
                                     + CNDY(3,j+1) * var(k,i,j  ,I_MOMZ) &
                                     - CNDY(1,j  ) * var(k,i,j-1,I_MOMZ) )
    enddo
    enddo
    enddo

    ! x-momentum
    do j = JJS,   JJE
    do i = IIS,   IIE
    do k = KS,   KE-1
       num_diff(k,i,j,I_MOMX,ZDIR) = DIFF4 * CDZ(k)**4 &
                                   * ( CNDZ(1,k+1) * var(k+2,i,j,I_MOMX) &
                                     - CNDZ(2,k+1) * var(k+1,i,j,I_MOMX) &
                                     + CNDZ(3,k+1) * var(k  ,i,j,I_MOMX) &
                                     - CNDZ(1,k  ) * var(k-1,i,j,I_MOMX) )
    enddo
    enddo
    enddo
    do j = JJS,   JJE
    do i = IIS,   IIE+1
    do k = KS,   KE
       num_diff(k,i,j,I_MOMX,XDIR) = DIFF4 * ( 0.5D0*(CDX(i+1)+CDX(i)) )**4 &
                                   * ( CNMX(1,i  ) * var(k,i+1,j,I_MOMX) &
                                     - CNMX(2,i  ) * var(k,i  ,j,I_MOMX) &
                                     + CNMX(3,i  ) * var(k,i-1,j,I_MOMX) &
                                     - CNMX(1,i-1) * var(k,i-2,j,I_MOMX) )
    enddo
    enddo
    enddo
    !OCL SERIAL
    do j = JJS-1, JJE
    !OCL PARALLEL
    do i = IIS,   IIE
    do k = KS,   KE
       num_diff(k,i,j,I_MOMX,YDIR) = DIFF4 * CDY(j)**4                   &
                                   * ( CNDY(1,j+1) * var(k,i,j+2,I_MOMX) &
                                     - CNDY(2,j+1) * var(k,i,j+1,I_MOMX) &
                                     + CNDY(3,j+1) * var(k,i,j  ,I_MOMX) &
                                     - CNDY(1,j  ) * var(k,i,j-1,I_MOMX) )
    enddo
    enddo
    enddo

    ! y-momentum
    do j = JJS,   JJE
    do i = IIS,   IIE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMY,ZDIR) = DIFF4 * CDZ(k)**4 &
                                   * ( CNDZ(1,k+1) * var(k+2,i,j,I_MOMY) &
                                     - CNDZ(2,k+1) * var(k+1,i,j,I_MOMY) &
                                     + CNDZ(3,k+1) * var(k  ,i,j,I_MOMY) &
                                     - CNDZ(1,k  ) * var(k-1,i,j,I_MOMY) )
    enddo
    enddo
    enddo
    do j = JJS,   JJE
    do i = IIS-1, IIE
    do k = KS,   KE
       num_diff(k,i,j,I_MOMY,XDIR) = DIFF4 * CDX(i)**4 &
                                   * ( CNDX(1,i+1) * var(k,i+2,j,I_MOMY) &
                                     - CNDX(2,i+1) * var(k,i+1,j,I_MOMY) &
                                     + CNDX(3,i+1) * var(k,i  ,j,I_MOMY) &
                                     - CNDX(1,i  ) * var(k,i-1,j,I_MOMY) )
    enddo
    enddo
    enddo
    !OCL SERIAL
    do j = JJS, JJE+1
    !OCL PARALLEL
    do i = IIS, IIE
    do k = KS, KE
       num_diff(k,i,j,I_MOMY,YDIR) = DIFF4 * ( 0.5D0*(CDY(j+1)+CDY(j)) )**4 &
                                   * ( CNMY(1,j  ) * var(k,i,j+1,I_MOMY) &
                                     - CNMY(2,j  ) * var(k,i,j  ,I_MOMY) &
                                     + CNMY(3,j  ) * var(k,i,j-1,I_MOMY) &
                                     - CNMY(1,j-1) * var(k,i,j-2,I_MOMY) )
    enddo
    enddo
    enddo

    enddo
    enddo

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
          diagvar(k,i,j,I_VELZ) = 2.D0 * var(k,i,j,I_MOMZ) &
                                / ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) )
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          diagvar(KS-1,i,j,I_VELZ) = 0.D0
          diagvar(KE  ,i,j,I_VELZ) = 0.D0
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          diagvar(k,i,j,I_VELX) = 2.D0 * var(k,i,j,I_MOMX) &
                                / ( var(k,i+1,j,I_DENS)+var(k,i,j,I_DENS) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS,   IIE+1
       do k = KS, KE
          diagvar(k,i,j,I_VELY) = 2.D0 * var(k,i,j,I_MOMY) &
                                / ( var(k,i,j+1,I_DENS)+var(k,i,j,I_DENS) )
       enddo
       enddo
       enddo

       ! pressure, pott. temp.
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          diagvar(k,i,j,I_PRES) = Pstd * ( var(k,i,j,I_RHOT) * Rdry / Pstd )**CPovCV
       enddo
       do k = KS, KE
          diagvar(k,i,j,I_POTT) = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS) 
       enddo
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          mflx_hi(k,i,j,ZDIR) = var(k,i,j,I_MOMZ) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS,   JJE
       do i = IIS,   IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.D0 ! bottom boundary
          mflx_hi(KE  ,i,j,ZDIR) = 0.D0 ! top    boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = var(k,i,j,I_MOMX) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = var(k,i,j,I_MOMY) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update density
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK1(k,i,j,I_DENS) = var(k,i,j,I_DENS) &
                                + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                             + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                             + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) ) ! divergence
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELZ)+diagvar(k-1,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k  ,i,j,I_MOMZ)+var(k-1,i,j,I_MOMZ) )     &
                                + FACT_F * ( var(k+1,i,j,I_MOMZ)+var(k-2,i,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i+1,j,I_MOMZ)+var(k,i  ,j,I_MOMZ) )     &
                                + FACT_F * ( var(k,i+2,j,I_MOMZ)+var(k,i-1,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j+1,I_MOMZ)+var(k,i,j  ,I_MOMZ) )     &
                                + FACT_F * ( var(k,i,j+2,I_MOMZ)+var(k,i,j-1,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          var_RK1(k,i,j,I_MOMZ) = var(k,i,j,I_MOMZ) &
                                + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RDZF(k)   &
                                             + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RDXC(i)   &
                                             + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                           - ( diagvar(k+1,i,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDZF(k)     & ! pressure gradient force
                                           - ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) ) * 0.5D0 * GRAV        & ! gravity force
                                           + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k+1,i,j,I_MOMX)+var(k  ,i,j,I_MOMX) )     &
                                + FACT_F * ( var(k+2,i,j,I_MOMX)+var(k-1,i,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( diagvar(KS  ,i+1,j,I_VELZ)+diagvar(KS  ,i,j,I_VELZ) ) & ! just above the bottom boundary
                                 * ( var(KS+1,i,j,I_MOMX)+var(KS,i,j,I_MOMX) )                      &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(KE-1,i+1,j,I_VELZ)+diagvar(KE-1,i,j,I_VELZ) ) & ! just below the top boundary
                                 * ( var(KE,i,j,I_MOMX)+var(KE-1,i,j,I_MOMX) )                      &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELX)+diagvar(k,i-1,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i  ,j,I_MOMX)+var(k,i-1,j,I_MOMX) )     &
                                + FACT_F * ( var(k,i+1,j,I_MOMX)+var(k,i-2,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j+1,I_MOMX)+var(k,i,j  ,I_MOMX) )     &
                                + FACT_F * ( var(k,i,j+2,I_MOMX)+var(k,i,j-1,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK1(k,i,j,I_MOMX) = var(k,i,j,I_MOMX) &
                                + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RDZC(k)   &
                                             + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RDXF(i)   &
                                             + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                           - ( diagvar(k,i+1,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDXF(i)     & ! pressure gradient force
                                           + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k+1,i,j,I_MOMY)+var(k  ,i,j,I_MOMY) )     &
                                + FACT_F * ( var(k+2,i,j,I_MOMY)+var(k-1,i,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( diagvar(KS  ,i,j+1,I_VELZ)+diagvar(KS  ,i,j,I_VELZ) ) & ! just above the bottom boundary
                                 * ( var(KS+1,i,j,I_MOMY)+var(KS,i,j,I_MOMY) )                      &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(KE-1,i,j+1,I_VELZ)+diagvar(KE-1,i,j,I_VELZ) ) & ! just below the top boundary
                                 * ( var(KE,i,j,I_MOMY)+var(KE-1,i,j,I_MOMY) )                      &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i+1,j,I_MOMY)+var(k,i  ,j,I_MOMY) )     &
                                + FACT_F * ( var(k,i+2,j,I_MOMY)+var(k,i-1,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       !OCL SERIAL
       do j = JJS, JJE+1
       !OCL PARALLEL
       do i = IIS, IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELY)+diagvar(k,i,j-1,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j  ,I_MOMY)+var(k,i,j-1,I_MOMY) )     &
                                + FACT_F * ( var(k,i,j+1,I_MOMY)+var(k,i,j-2,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK1(k,i,j,I_MOMY) = var(k,i,j,I_MOMY) &
                                + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RDZC(k)   &
                                             + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RDXC(i)   &
                                             + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RDYF(j) ) & ! flux divergence
                                           - ( diagvar(k,i,j+1,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDYF(j)     & ! pressure gradient force
                                           + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                                      &
                              * ( FACT_N * ( diagvar(k+1,i,j,I_POTT)+diagvar(k  ,i,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k+2,i,j,I_POTT)+diagvar(k-1,i,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                  ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR)                      & ! just above the bottom boundary
                                 * ( diagvar(KS+1,i,j,I_POTT)+diagvar(KS,i,j,I_POTT) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR)                      & ! just below the top boundary
                                 * ( diagvar(KE,i,j,I_POTT)+diagvar(KE-1,i,j,I_POTT) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                  ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i+1,j,I_POTT)+diagvar(k,i  ,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i+2,j,I_POTT)+diagvar(k,i-1,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i,j+1,I_POTT)+diagvar(k,i,j  ,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i,j+2,I_POTT)+diagvar(k,i,j-1,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK1(k,i,j,I_RHOT) = var(k,i,j,I_RHOT) &
                                + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                             + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                             + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) & ! divergence
                                           + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping force
       enddo
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( var_RK1(:,:,:,1), 1 )
    call COMM_vars8( var_RK1(:,:,:,2), 2 )
    call COMM_vars8( var_RK1(:,:,:,3), 3 )
    call COMM_vars8( var_RK1(:,:,:,4), 4 )
    call COMM_vars8( var_RK1(:,:,:,5), 5 )
    call COMM_wait ( var_RK1(:,:,:,1), 1 )
    call COMM_wait ( var_RK1(:,:,:,2), 2 )
    call COMM_wait ( var_RK1(:,:,:,3), 3 )
    call COMM_wait ( var_RK1(:,:,:,4), 4 )
    call COMM_wait ( var_RK1(:,:,:,5), 5 )

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
       do k = KS,   KE-1
          diagvar(k,i,j,I_VELZ) = 2.D0 * var_RK1(k,i,j,I_MOMZ) &
                                / ( var_RK1(k+1,i,j,I_DENS)+var_RK1(k,i,j,I_DENS) )
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          diagvar(KS-1,i,j,I_VELZ) = 0.D0
          diagvar(KE  ,i,j,I_VELZ) = 0.D0
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          diagvar(k,i,j,I_VELX) = 2.D0 * var_RK1(k,i,j,I_MOMX) &
                                / ( var_RK1(k,i+1,j,I_DENS)+var_RK1(k,i,j,I_DENS) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS,   IIE+1
       do k = KS, KE
          diagvar(k,i,j,I_VELY) = 2.D0 * var_RK1(k,i,j,I_MOMY) &
                                / ( var_RK1(k,i,j+1,I_DENS)+var_RK1(k,i,j,I_DENS) )
       enddo
       enddo
       enddo

       ! pressure, pott. temp.
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          diagvar(k,i,j,I_PRES) = Pstd * ( var_RK1(k,i,j,I_RHOT) * Rdry / Pstd )**CPovCV
       enddo
       do k = KS, KE
          diagvar(k,i,j,I_POTT) = var_RK1(k,i,j,I_RHOT) / var_RK1(k,i,j,I_DENS) 
       enddo
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          mflx_hi(k,i,j,ZDIR) = var_RK1(k,i,j,I_MOMZ) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS,   JJE
       do i = IIS,   IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.D0 ! bottom boundary
          mflx_hi(KE  ,i,j,ZDIR) = 0.D0 ! top    boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = var_RK1(k,i,j,I_MOMX) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = var_RK1(k,i,j,I_MOMY) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update density
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK2(k,i,j,I_DENS) = var(k,i,j,I_DENS) &
                                + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                             + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                             + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) ) ! divergence
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELZ)+diagvar(k-1,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var_RK1(k  ,i,j,I_MOMZ)+var_RK1(k-1,i,j,I_MOMZ) )     &
                                + FACT_F * ( var_RK1(k+1,i,j,I_MOMZ)+var_RK1(k-2,i,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var_RK1(k,i+1,j,I_MOMZ)+var_RK1(k,i  ,j,I_MOMZ) )     &
                                + FACT_F * ( var_RK1(k,i+2,j,I_MOMZ)+var_RK1(k,i-1,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var_RK1(k,i,j+1,I_MOMZ)+var_RK1(k,i,j  ,I_MOMZ) )     &
                                + FACT_F * ( var_RK1(k,i,j+2,I_MOMZ)+var_RK1(k,i,j-1,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          var_RK2(k,i,j,I_MOMZ) = var(k,i,j,I_MOMZ) &
                                + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RDZF(k)   &
                                             + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RDXC(i)   &
                                             + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                           - ( diagvar(k+1,i,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDZF(k)     & ! pressure gradient force
                                           - ( var_RK1(k+1,i,j,I_DENS)+var_RK1(k,i,j,I_DENS) ) * 0.5D0 * GRAV        & ! gravity force
                                           + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var_RK1(k+1,i,j,I_MOMX)+var_RK1(k  ,i,j,I_MOMX) )     &
                                + FACT_F * ( var_RK1(k+2,i,j,I_MOMX)+var_RK1(k-1,i,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( diagvar(KS  ,i+1,j,I_VELZ)+diagvar(KS  ,i,j,I_VELZ) ) & ! just above the bottom boundary
                                 * ( var_RK1(KS+1,i,j,I_MOMX)+var_RK1(KS,i,j,I_MOMX) )                      &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(KE-1,i+1,j,I_VELZ)+diagvar(KE-1,i,j,I_VELZ) ) & ! just below the top boundary
                                 * ( var_RK1(KE,i,j,I_MOMX)+var_RK1(KE-1,i,j,I_MOMX) )                      &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELX)+diagvar(k,i-1,j,I_VELX) ) &
                              * ( FACT_N * ( var_RK1(k,i  ,j,I_MOMX)+var_RK1(k,i-1,j,I_MOMX) )     &
                                + FACT_F * ( var_RK1(k,i+1,j,I_MOMX)+var_RK1(k,i-2,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var_RK1(k,i,j+1,I_MOMX)+var_RK1(k,i,j  ,I_MOMX) )     &
                                + FACT_F * ( var_RK1(k,i,j+2,I_MOMX)+var_RK1(k,i,j-1,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK2(k,i,j,I_MOMX) = var(k,i,j,I_MOMX) &
                                + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RDZC(k)   &
                                             + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RDXF(i)   &
                                             + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                           - ( diagvar(k,i+1,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDXF(i)     & ! pressure gradient force
                                           + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var_RK1(k+1,i,j,I_MOMY)+var_RK1(k  ,i,j,I_MOMY) )     &
                                + FACT_F * ( var_RK1(k+2,i,j,I_MOMY)+var_RK1(k-1,i,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( diagvar(KS  ,i,j+1,I_VELZ)+diagvar(KS  ,i,j,I_VELZ) ) & ! just above the bottom boundary
                                 * ( var_RK1(KS+1,i,j,I_MOMY)+var_RK1(KS,i,j,I_MOMY) )                      &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(KE-1,i,j+1,I_VELZ)+diagvar(KE-1,i,j,I_VELZ) ) & ! just below the top boundary
                                 * ( var_RK1(KE,i,j,I_MOMY)+var_RK1(KE-1,i,j,I_MOMY) )                      &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var_RK1(k,i+1,j,I_MOMY)+var_RK1(k,i  ,j,I_MOMY) )     &
                                + FACT_F * ( var_RK1(k,i+2,j,I_MOMY)+var_RK1(k,i-1,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       !OCL SERIAL
       do j = JJS, JJE+1
       !OCL PARALLEL
       do i = IIS, IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELY)+diagvar(k,i,j-1,I_VELY) ) &
                              * ( FACT_N * ( var_RK1(k,i,j  ,I_MOMY)+var_RK1(k,i,j-1,I_MOMY) )     &
                                + FACT_F * ( var_RK1(k,i,j+1,I_MOMY)+var_RK1(k,i,j-2,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK2(k,i,j,I_MOMY) = var(k,i,j,I_MOMY) &
                                + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RDZC(k)   &
                                             + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RDXC(i)   &
                                             + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RDYF(j) ) & ! flux divergence
                                           - ( diagvar(k,i,j+1,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDYF(j)     & ! pressure gradient force
                                           + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                                      &
                              * ( FACT_N * ( diagvar(k+1,i,j,I_POTT)+diagvar(k  ,i,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k+2,i,j,I_POTT)+diagvar(k-1,i,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                  ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR)                      & ! just above the bottom boundary
                                 * ( diagvar(KS+1,i,j,I_POTT)+diagvar(KS,i,j,I_POTT) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR)                      & ! just below the top boundary
                                 * ( diagvar(KE,i,j,I_POTT)+diagvar(KE-1,i,j,I_POTT) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                  ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i+1,j,I_POTT)+diagvar(k,i  ,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i+2,j,I_POTT)+diagvar(k,i-1,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i,j+1,I_POTT)+diagvar(k,i,j  ,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i,j+2,I_POTT)+diagvar(k,i,j-1,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var_RK2(k,i,j,I_RHOT) = var(k,i,j,I_RHOT) &
                                + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                             + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                             + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) & ! divergence
                                           + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping force
       enddo
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( var_RK2(:,:,:,1), 1 )
    call COMM_vars8( var_RK2(:,:,:,2), 2 )
    call COMM_vars8( var_RK2(:,:,:,3), 3 )
    call COMM_vars8( var_RK2(:,:,:,4), 4 )
    call COMM_vars8( var_RK2(:,:,:,5), 5 )
    call COMM_wait ( var_RK2(:,:,:,1), 1 )
    call COMM_wait ( var_RK2(:,:,:,2), 2 )
    call COMM_wait ( var_RK2(:,:,:,3), 3 )
    call COMM_wait ( var_RK2(:,:,:,4), 4 )
    call COMM_wait ( var_RK2(:,:,:,5), 5 )

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
       do k = KS,   KE-1
          diagvar(k,i,j,I_VELZ) = 2.D0 * var_RK2(k,i,j,I_MOMZ) &
                                / ( var_RK2(k+1,i,j,I_DENS)+var_RK2(k,i,j,I_DENS) )
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          diagvar(KS-1,i,j,I_VELZ) = 0.D0
          diagvar(KE  ,i,j,I_VELZ) = 0.D0
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          diagvar(k,i,j,I_VELX) = 2.D0 * var_RK2(k,i,j,I_MOMX) &
                                / ( var_RK2(k,i+1,j,I_DENS)+var_RK2(k,i,j,I_DENS) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS,   IIE+1
       do k = KS, KE
          diagvar(k,i,j,I_VELY) = 2.D0 * var_RK2(k,i,j,I_MOMY) &
                                / ( var_RK2(k,i,j+1,I_DENS)+var_RK2(k,i,j,I_DENS) )
       enddo
       enddo
       enddo

       ! pressure, pott. temp.
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE
          diagvar(k,i,j,I_PRES) = Pstd * ( var_RK2(k,i,j,I_RHOT) * Rdry / Pstd )**CPovCV
       enddo
       do k = KS, KE
          diagvar(k,i,j,I_POTT) = var_RK2(k,i,j,I_RHOT) / var_RK2(k,i,j,I_DENS) 
       enddo
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          mflx_hi(k,i,j,ZDIR) = var_RK2(k,i,j,I_MOMZ) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS,   JJE
       do i = IIS,   IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.D0 ! bottom boundary
          mflx_hi(KE  ,i,j,ZDIR) = 0.D0 ! top    boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = var_RK2(k,i,j,I_MOMX) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = var_RK2(k,i,j,I_MOMY) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update density
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var(k,i,j,I_DENS) = var(k,i,j,I_DENS) &
                            + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                         + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                         + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) ) ! divergence
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELZ)+diagvar(k-1,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var_RK2(k  ,i,j,I_MOMZ)+var_RK2(k-1,i,j,I_MOMZ) )     &
                                + FACT_F * ( var_RK2(k+1,i,j,I_MOMZ)+var_RK2(k-2,i,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var_RK2(k,i+1,j,I_MOMZ)+var_RK2(k,i  ,j,I_MOMZ) )     &
                                + FACT_F * ( var_RK2(k,i+2,j,I_MOMZ)+var_RK2(k,i-1,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var_RK2(k,i,j+1,I_MOMZ)+var_RK2(k,i,j  ,I_MOMZ) )     &
                                + FACT_F * ( var_RK2(k,i,j+2,I_MOMZ)+var_RK2(k,i,j-1,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS,   KE-1
          var(k,i,j,I_MOMZ) = var(k,i,j,I_MOMZ) &
                            + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RDZF(k)   &
                                         + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RDXC(i)   &
                                         + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                       - ( diagvar(k+1,i,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDZF(k)     & ! pressure gradient force
                                       - ( var_RK2(k+1,i,j,I_DENS)+var_RK2(k,i,j,I_DENS) ) * 0.5D0 * GRAV        & ! gravity force
                                       + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var_RK2(k+1,i,j,I_MOMX)+var_RK2(k  ,i,j,I_MOMX) )     &
                                + FACT_F * ( var_RK2(k+2,i,j,I_MOMX)+var_RK2(k-1,i,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( diagvar(KS  ,i+1,j,I_VELZ)+diagvar(KS  ,i,j,I_VELZ) ) & ! just above the bottom boundary
                                 * ( var_RK2(KS+1,i,j,I_MOMX)+var_RK2(KS,i,j,I_MOMX) )                      &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(KE-1,i+1,j,I_VELZ)+diagvar(KE-1,i,j,I_VELZ) ) & ! just below the top boundary
                                 * ( var_RK2(KE,i,j,I_MOMX)+var_RK2(KE-1,i,j,I_MOMX) )                      &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS,   JJE
       do i = IIS,   IIE+1
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELX)+diagvar(k,i-1,j,I_VELX) ) &
                              * ( FACT_N * ( var_RK2(k,i  ,j,I_MOMX)+var_RK2(k,i-1,j,I_MOMX) )     &
                                + FACT_F * ( var_RK2(k,i+1,j,I_MOMX)+var_RK2(k,i-2,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var_RK2(k,i,j+1,I_MOMX)+var_RK2(k,i,j  ,I_MOMX) )     &
                                + FACT_F * ( var_RK2(k,i,j+2,I_MOMX)+var_RK2(k,i,j-1,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var(k,i,j,I_MOMX) = var(k,i,j,I_MOMX) &
                            + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RDZC(k)   &
                                         + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RDXF(i)   &
                                         + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                       - ( diagvar(k,i+1,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDXF(i)     & ! pressure gradient force
                                       + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var_RK2(k+1,i,j,I_MOMY)+var_RK2(k  ,i,j,I_MOMY) )     &
                                + FACT_F * ( var_RK2(k+2,i,j,I_MOMY)+var_RK2(k-1,i,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.25D0 * ( diagvar(KS  ,i,j+1,I_VELZ)+diagvar(KS  ,i,j,I_VELZ) ) & ! just above the bottom boundary
                                 * ( var_RK2(KS+1,i,j,I_MOMY)+var_RK2(KS,i,j,I_MOMY) )                      &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(KE-1,i,j+1,I_VELZ)+diagvar(KE-1,i,j,I_VELZ) ) & ! just below the top boundary
                                 * ( var_RK2(KE,i,j,I_MOMY)+var_RK2(KE-1,i,j,I_MOMY) )                      &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var_RK2(k,i+1,j,I_MOMY)+var_RK2(k,i  ,j,I_MOMY) )     &
                                + FACT_F * ( var_RK2(k,i+2,j,I_MOMY)+var_RK2(k,i-1,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       !OCL SERIAL
       do j = JJS, JJE+1
       !OCL PARALLEL
       do i = IIS, IIE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELY)+diagvar(k,i,j-1,I_VELY) ) &
                              * ( FACT_N * ( var_RK2(k,i,j  ,I_MOMY)+var_RK2(k,i,j-1,I_MOMY) )     &
                                + FACT_F * ( var_RK2(k,i,j+1,I_MOMY)+var_RK2(k,i,j-2,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var(k,i,j,I_MOMY) = var(k,i,j,I_MOMY) &
                            + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RDZC(k)   &
                                         + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RDXC(i)   &
                                         + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RDYF(j) ) & ! flux divergence
                                       - ( diagvar(k,i,j+1,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDYF(j)     & ! pressure gradient force
                                       + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                                      &
                              * ( FACT_N * ( diagvar(k+1,i,j,I_POTT)+diagvar(k  ,i,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k+2,i,j,I_POTT)+diagvar(k-1,i,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0                                                  ! bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR)                      & ! just above the bottom boundary
                                 * ( diagvar(KS+1,i,j,I_POTT)+diagvar(KS,i,j,I_POTT) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR)                      & ! just below the top boundary
                                 * ( diagvar(KE,i,j,I_POTT)+diagvar(KE-1,i,j,I_POTT) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(KE  ,i,j,ZDIR) = 0.D0                                                  ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i+1,j,I_POTT)+diagvar(k,i  ,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i+2,j,I_POTT)+diagvar(k,i-1,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i,j+1,I_POTT)+diagvar(k,i,j  ,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i,j+2,I_POTT)+diagvar(k,i,j-1,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo

       !--- update rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          var(k,i,j,I_RHOT) = var(k,i,j,I_RHOT) &
                            + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                         + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                         + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) & ! divergence
                                       + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping force
       enddo
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( var(:,:,:,1), 1 )
    call COMM_vars8( var(:,:,:,2), 2 )
    call COMM_vars8( var(:,:,:,3), 3 )
    call COMM_vars8( var(:,:,:,4), 4 )
    call COMM_vars8( var(:,:,:,5), 5 )
    call COMM_wait ( var(:,:,:,1), 1 )
    call COMM_wait ( var(:,:,:,2), 2 )
    call COMM_wait ( var(:,:,:,3), 3 )
    call COMM_wait ( var(:,:,:,4), 4 )
    call COMM_wait ( var(:,:,:,5), 5 )

#ifdef _FPCOLL_
call STOP_COLLECTION("RK3")
call START_COLLECTION("FCT")
#endif

    !##### advection of scalar quantity #####

    if ( QA > 0 ) then
    do iq = 6, 5+QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
          qflx_lo(k,i,j,ZDIR) = 0.5D0 * (     mflx_hi(k,i,j,ZDIR)  * ( var(k+1,i,j,iq)+var(k,i,j,iq) ) &
                                        - abs(mflx_hi(k,i,j,ZDIR)) * ( var(k+1,i,j,iq)-var(k,i,j,iq) ) )

          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                    &
                              * ( FACT_N * ( var(k+1,i,j,iq)+var(k  ,i,j,iq) ) &
                                + FACT_F * ( var(k+2,i,j,iq)+var(k-1,i,j,iq) ) )

          qflx_anti(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_lo(k,i,j,ZDIR)
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_lo(KS-1,i,j,ZDIR) = 0.D0                                                                          ! bottom boundary
          qflx_lo(KS  ,i,j,ZDIR) = 0.5D0 * (     mflx_hi(KS  ,i,j,ZDIR)  * ( var(KS+1,i,j,iq)+var(KS,i,j,iq) ) & ! just above the bottom boundary
                                           - abs(mflx_hi(KS  ,i,j,ZDIR)) * ( var(KS+1,i,j,iq)-var(KS,i,j,iq) ) )
          qflx_lo(KE-1,i,j,ZDIR) = 0.5D0 * (     mflx_hi(KE-1,i,j,ZDIR)  * ( var(KE,i,j,iq)+var(KE-1,i,j,iq) ) & ! just below the top boundary
                                           - abs(mflx_hi(KE-1,i,j,ZDIR)) * ( var(KE,i,j,iq)-var(KE-1,i,j,iq) ) )
          qflx_lo(KE  ,i,j,ZDIR) = 0.D0                                                                          ! top boundary

          qflx_hi(KS-1,i,j,ZDIR) = 0.D0 
          qflx_hi(KS  ,i,j,ZDIR) = 0.5D0 * mflx_hi(KS  ,i,j,ZDIR) * ( var(KS+1,i,j,iq)+var(KS,i,j,iq) )
          qflx_hi(KE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(KE-1,i,j,ZDIR) * ( var(KE,i,j,iq)+var(KE-1,i,j,iq) )
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
          qflx_lo(k,i,j,XDIR) = 0.5D0 * (     mflx_hi(k,i,j,XDIR)  * ( var(k,i+1,j,iq)+var(k,i,j,iq) ) &
                                        - abs(mflx_hi(k,i,j,XDIR)) * ( var(k,i+1,j,iq)-var(k,i,j,iq) ) )

          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                    &
                              * ( FACT_N * ( var(k,i+1,j,iq)+var(k,i  ,j,iq) ) &
                                + FACT_F * ( var(k,i+2,j,iq)+var(k,i-1,j,iq) ) )

          qflx_anti(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_lo(k,i,j,XDIR)
       enddo
       enddo
       enddo
       !OCL SERIAL
       do j = JJS-1, JJE
       !OCL PARALLEL
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_lo(k,i,j,YDIR) = 0.5D0 * (     mflx_hi(k,i,j,YDIR)  * ( var(k,i,j+1,iq)+var(k,i,j,iq) ) &
                                        - abs(mflx_hi(k,i,j,YDIR)) * ( var(k,i,j+1,iq)-var(k,i,j,iq) ) )

          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                    &
                              * ( FACT_N * ( var(k,i,j+1,iq)+var(k,i,j  ,iq) ) &
                                + FACT_F * ( var(k,i,j+2,iq)+var(k,i,j-1,iq) ) )

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
             tmp = var(k,i,j,iq) / pjmns * rdtrk
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
       do j = JJS-1, JJE
       do i = IIS-1, IIE
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
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS-1, KE
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
       do i = IIS-1, IIE
       do k = KS-1, KE
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
          var(k,i,j,iq) = ( var(k,i,j,iq) * dens_s(k,i,j)                                                &
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)     &
                                       + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)     &
                                       + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) ) &
                          ) / var(k,i,j,I_DENS)
       enddo
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( var(:,:,:,iq), iq )
    call COMM_wait ( var(:,:,:,iq), iq )

    enddo ! scalar quantities loop
    endif

#ifdef _FPCOLL_
call STOP_COLLECTION("FCT")
#endif

    enddo ! dynamical steps

#ifdef _FPCOLL_
call STOP_COLLECTION("DYNAMICS")
#endif

    ! check total mass
    call COMM_total( var(:,:,:,:), A_NAME(:) )

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
