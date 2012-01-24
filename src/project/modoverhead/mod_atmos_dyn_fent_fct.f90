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
  integer, private, parameter :: XDIR   = 1
  integer, private, parameter :: YDIR   = 2
  integer, private, parameter :: ZDIR   = 3

  integer, private, parameter :: I_DENS = 1
  integer, private, parameter :: I_MOMX = 2
  integer, private, parameter :: I_MOMY = 3
  integer, private, parameter :: I_MOMZ = 4
  integer, private, parameter :: I_RHOT = 5

  integer, private, parameter :: I_PRES = 1
  integer, private, parameter :: I_VELX = 2
  integer, private, parameter :: I_VELY = 3
  integer, private, parameter :: I_VELZ = 4
  integer, private, parameter :: I_POTT = 5

  ! time settings
  integer, private, parameter :: RK = 3 ! order of Runge-Kutta scheme

  ! advection settings
  real(8), private, parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), private, parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second
  real(8), private, parameter :: FACT_M =   9.D0 / 8.D0 !  7/6: fourth, 1: second
  real(8), private, parameter :: FACT_G = - 1.D0 / 8.D0 ! -1/6: fourth, 0: second

  ! numerical filter settings
  integer, private, parameter :: DF   = 4     ! order of numerical filter

  real(8), private, save      :: ATMOS_DYN_numerical_diff = 1.D-3 ! nondimensional numerical diffusion
  real(8), private, save      :: DIFF

  real(8), private, allocatable, save :: CNDx1(:), CNDx2(:), CNDx3(:)
  real(8), private, allocatable, save :: CNMx1(:), CNMx2(:), CNMx3(:)
  real(8), private, allocatable, save :: CNDy1(:), CNDy2(:), CNDy3(:)
  real(8), private, allocatable, save :: CNMy1(:), CNMy2(:), CNMy3(:)
  real(8), private, allocatable, save :: CNDz1(:), CNDz2(:), CNDz3(:)
  real(8), private, allocatable, save :: CNMz1(:), CNMz2(:), CNMz3(:)

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
       KA  => GRID_KA, &
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KS  => GRID_KS, &
       KE  => GRID_KE, &
       WS  => GRID_WS, &
       WE  => GRID_WE, &
       IS  => GRID_IS, &
       IE  => GRID_IE, &
       JS  => GRID_JS, &
       JE  => GRID_JE, &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY
    implicit none

    NAMELIST / PARAM_ATMOS_DYN / &
       ATMOS_DYN_numerical_diff

    integer :: ierr
    integer :: i, j, k
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

    DIFF = - ATMOS_DYN_numerical_diff * (-1.D0)**dble( DF/2+1 )

    allocate( CNDx1(IA) )
    allocate( CNDx2(IA) )
    allocate( CNDx3(IA) )
    allocate( CNMx1(IA) )
    allocate( CNMx2(IA) )
    allocate( CNMx3(IA) )

    allocate( CNDy1(JA) )
    allocate( CNDy2(JA) )
    allocate( CNDy3(JA) )
    allocate( CNMy1(JA) )
    allocate( CNMy2(JA) )
    allocate( CNMy3(JA) )

    allocate( CNDz1(KA) )
    allocate( CNDz2(KA) )
    allocate( CNDz3(KA) )
    allocate( CNMz1(KA) )
    allocate( CNMz2(KA) )
    allocate( CNMz3(KA) )

    ! z djrectjon
    do k = KS-1, KE+1
       CNDz1(k) = 1.D0 / ( (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 )
       CNDz2(k) = 1.D0 / ( (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 )
    enddo
    CNDz1(1)  = CNDz1(KS)
    CNDz2(1)  = CNDz2(KS)
    CNDz1(KA) = CNDz1(KE)
    CNDz2(KA) = CNDz2(KE)

    do k = KS, KE+2
       CNDz3(k) = 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) * (CDZ(k-1)+CDZ(k-2)) * 0.5D0 )
    enddo
    CNDz3(1)    = CNDz3(KS)
    CNDz3(KS-1) = CNDz3(KS)

    do k = KS-2, KE+1
       CNMz1(k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) )
    enddo
    CNMz1(KA) = CNMz1(KE)

    do k = KS-1, KE+1
       CNMz2(k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &   
                + 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &  
                + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) 
       CNMz3(k) = 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &
                + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) &
                + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) )
    enddo
    CNMz2(1)  = CNMz2(KS)
    CNMz3(1)  = CNMz3(KS)
    CNMz2(KA) = CNMz2(KE)
    CNMz3(KA) = CNMz3(KE)

    ! x direction
    do i = IS-1, IE+1
       CNDx1(i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
       CNDx2(i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
    enddo
    CNDx1(1)  = CNDx1(IS)
    CNDx2(1)  = CNDx2(IS)
    CNDx1(IA) = CNDx1(IE)
    CNDx2(IA) = CNDx2(IE)

    do i = IS, IE+2
       CNDx3(i) = 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i-1)+CDX(i-2)) * 0.5D0 )
    enddo
    CNDx3(1)    = CNDx3(IS)
    CNDx3(IS-1) = CNDx3(IS)

    do i = IS-2, IE+1
       CNMx1(i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) )
    enddo
    CNMx1(IA) = CNMx1(IE)

    do i = IS-1, IE+1
       CNMx2(i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                + 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) )
       CNMx3(i) = 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) ) & 
                + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) )
    enddo
    CNMx2(1)  = CNMx2(IS)
    CNMx3(1)  = CNMx3(IS)
    CNMx2(IA) = CNMx2(IE)
    CNMx3(IA) = CNMx3(IE)

    ! y direction
    do j = JS, JE
       CNDy1(j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
       CNDy2(j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
       CNDy3(j) = 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j-1)+CDY(j-2)) * 0.5D0 )
       CNMy1(j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) )
       CNMy2(j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &   
                + 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &  
                + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) 
       CNMy3(j) = 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &
                + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) &
                + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) )
    enddo
    do j = 1, JS-1
       CNDy1(j) = CNDy1(JS)
       CNDy2(j) = CNDy2(JS)
       CNDy3(j) = CNDy3(JS)
       CNMy1(j) = CNMy1(JS)
       CNMy2(j) = CNMy2(JS)
       CNMy3(j) = CNMy3(JS)
    enddo
    do j = JE+1, JA
       CNDy1(j) = CNDy1(JE)
       CNDy2(j) = CNDy2(JE)
       CNDy3(j) = CNDy3(JE)
       CNMy1(j) = CNMy1(JE)
       CNMy2(j) = CNMy2(JE)
       CNMy3(j) = CNMy3(JE)
    enddo

  end subroutine ATMOS_DYN_setup

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN( &
     var, var_s, diagvar, &
     CDZ, CDX, CDY, RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY, &
     ray_damp, DAMP_var, DAMP_alpha, &
     num_diff, dens_diff, pott_diff, REF_dens, REF_pott, &
     mflx_hi, mflx_lo, qflx_hi, qflx_lo, qflx_anti, ddiv, rjpls, rjmns, &
     KA, IA, JA, VA, QA, &
     KS, KE, WS, WE, IS, IE, JS, JE )
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
       COMM_wait, &
       COMM_total
    use mod_atmos_vars, only: &
       A_NAME
    use mod_atmos_boundary, only: &
       I_BND_VELX, &
       I_BND_VELY, &
       I_BND_VELZ, &
       I_BND_POTT
    implicit none

    ! array size
    integer, intent(in)    :: KA
    integer, intent(in)    :: IA
    integer, intent(in)    :: JA
    integer, intent(in)    :: VA
    integer, intent(in)    :: QA
    ! variables
    real(8), intent(inout) :: var    (KA,IA,JA,VA)  ! prognostic variables
    real(8), intent(inout) :: var_s  (KA,IA,JA,VA)  ! prognostic variables (previous step)
    real(8), intent(inout) :: diagvar(KA,IA,JA,5)   ! diagnostic variables (work)
    ! grid
    real(8), intent(in)    :: CDZ (KA)
    real(8), intent(in)    :: CDX (IA)
    real(8), intent(in)    :: CDY (JA)
    real(8), intent(in)    :: RCDZ(KA)
    real(8), intent(in)    :: RCDX(IA)
    real(8), intent(in)    :: RCDY(JA)
    real(8), intent(in)    :: RFDZ(KA)
    real(8), intent(in)    :: RFDX(IA)
    real(8), intent(in)    :: RFDY(JA)
    ! dumping
    real(8), intent(inout) :: ray_damp  (KA,IA,JA,5)
    real(8), intent(in)    :: DAMP_var  (KA,IA,JA,5)
    real(8), intent(in)    :: DAMP_alpha(KA,IA,JA,5)
    ! numerical diffusion
    real(8), intent(inout) :: num_diff (KA,IA,JA,5,3)
    real(8), intent(inout) :: dens_diff(KA,IA,JA)     ! anomary of density
    real(8), intent(inout) :: pott_diff(KA,IA,JA)     ! anomary of rho * pott
    real(8), intent(in)    :: REF_dens (KA)
    real(8), intent(in)    :: REF_pott (KA)
    ! mass flux
    real(8), intent(inout) :: mflx_hi  (KA,IA,JA,3)   ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8), intent(inout) :: mflx_lo  (KA,IA,JA,3)   ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8), intent(inout) :: qflx_hi  (KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(8), intent(inout) :: qflx_lo  (KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face low  order
    real(8), intent(inout) :: qflx_anti(KA,IA,JA,3)   ! rho * vel(x,y,z) * phi @ (u,v,w)-face antidiffusive
    real(8), intent(inout) :: ddiv     (KA,IA,JA)     ! density divergence
    ! flux correction
    real(8), intent(inout) :: rjpls    (KA,IA,JA,3)   ! plus  in (x,y,z)-direction
    real(8), intent(inout) :: rjmns    (KA,IA,JA,3)   ! minus in (x,y,z)-direction

    ! array index
    integer, intent(in)    :: KS, KE
    integer, intent(in)    :: WS, WE
    integer, intent(in)    :: IS, IE
    integer, intent(in)    :: JS, JE

    real(8) :: pjpls, pjmns, pjmax, pjmin, var_l
    real(8) :: dtrk, rdtrk
    integer :: i, j, k, iq, iv, rko, step
    !---------------------------------------------------------------------------

    do step = 1, TIME_NSTEP_ATMOS_DYN

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical small step:', step

    call START_COLLECTION("prologue")

    do iv = 1, VA
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       var_s(k,i,j,iv) = var(k,i,j,iv)
    enddo 
    enddo
    enddo
    enddo

    !--- prepare rayleigh damping coefficient
    do j = JS, JE
    do i = IS, IE
       do k = WS+1, WE-1
          ray_damp(k,i,j,I_MOMZ) = - DAMP_alpha(k,i,j,I_BND_VELZ) * ( var(k,i,j,I_MOMZ)  &
                                 - DAMP_var(k,i,j,I_BND_VELZ) * 0.5D0 * ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) ) )
       enddo
       do k = KS, KE
          ray_damp(k,i,j,I_MOMX) = - DAMP_alpha(k,i,j,I_BND_VELX) * ( var(k,i,j,I_MOMX)  &
                                 - DAMP_var(k,i,j,I_BND_VELX) * 0.5D0 * ( var(k,i+1,j,I_DENS)+var(k,i,j,I_DENS) ) )
          ray_damp(k,i,j,I_MOMY) = - DAMP_alpha(k,i,j,I_BND_VELY) * ( var(k,i,j,I_MOMY)  &
                                 - DAMP_var(k,i,j,I_BND_VELY) * 0.5D0 * ( var(k,i,j+1,I_DENS)+var(k,i,j,I_DENS) ) )
          ray_damp(k,i,j,I_RHOT) = - DAMP_alpha(k,i,j,I_BND_POTT) * ( var(k,i,j,I_RHOT)  &
                                 - DAMP_var(k,i,j,I_BND_POTT) * var(k,i,j,I_DENS) )
       enddo 
    enddo
    enddo

    !--- prepare numerical diffusion coefficient
    ! note: the value must divide by dtrk.
    do j  = 1, JA
    do i  = 1, IA
       do k = KS, KE
          dens_diff(k,i,j) = var(k,i,j,I_DENS) - REF_dens(k)
          pott_diff(k,i,j) = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS) - REF_pott(k)
       enddo
       do k = 1, KS-1 
          dens_diff(k,i,j) = dens_diff(KS,i,j)
          pott_diff(k,i,j) = pott_diff(KS,i,j)
       enddo
       do k = KE+1, KA
          dens_diff(k,i,j) = dens_diff(KE,i,j)
          pott_diff(k,i,j) = pott_diff(KE,i,j)
       enddo
    enddo
    enddo

    do j = JS,   JE
    do i = IS,   IE
    do k = WS+1, WE-1
       num_diff(k,i,j,I_DENS,ZDIR) = DIFF * CDZ(k)**DF                 &
                                   * ( CNDz1(k+1) * dens_diff(k+2,i,j) &
                                     - CNDz2(k+1) * dens_diff(k+1,i,j) &
                                     + CNDz3(k+1) * dens_diff(k  ,i,j) &
                                     - CNDz1(k  ) * dens_diff(k-1,i,j) )
       num_diff(k,i,j,I_RHOT,ZDIR) = DIFF * CDZ(k)**DF                 &
                                   * ( CNDz1(k+1) * pott_diff(k+2,i,j) &
                                     - CNDz2(k+1) * pott_diff(k+1,i,j) &
                                     + CNDz3(k+1) * pott_diff(k  ,i,j) &
                                     - CNDz1(k  ) * pott_diff(k-1,i,j) ) * var(k,i,j,I_DENS)
    enddo
    enddo
    enddo

    do j = JS,   JE
    do i = IS-1, IE
    do k = KS,   KE
       num_diff(k,i,j,I_DENS,XDIR) = DIFF * CDX(i)**DF                 &
                                   * ( CNDx1(i+1) * dens_diff(k,i+2,j) &
                                     - CNDx2(i+1) * dens_diff(k,i+1,j) &
                                     + CNDx3(i+1) * dens_diff(k,i  ,j) &
                                     - CNDx1(i  ) * dens_diff(k,i-1,j) )
       num_diff(k,i,j,I_RHOT,XDIR) = DIFF * CDX(i)**DF                 &
                                   * ( CNDx1(i+1) * pott_diff(k,i+2,j) &
                                     - CNDx2(i+1) * pott_diff(k,i+1,j) &
                                     + CNDx3(i+1) * pott_diff(k,i  ,j) &
                                     - CNDx1(i  ) * pott_diff(k,i-1,j) ) * var(k,i,j,I_DENS)

    enddo
    enddo
    enddo

    do j = JS-1, JE
    do i = IS,   IE
    do k = KS,   KE
       num_diff(k,i,j,I_DENS,YDIR) = DIFF * CDY(j)**DF                 &
                                   * ( CNDy1(j+1) * dens_diff(k,i,j+2) &
                                     - CNDy2(j+1) * dens_diff(k,i,j+1) &
                                     + CNDy3(j+1) * dens_diff(k,i,j  ) &
                                     - CNDy1(j  ) * dens_diff(k,i,j-1) )
       num_diff(k,i,j,I_RHOT,YDIR) = DIFF * CDY(j)**DF                 &
                                   * ( CNDy1(j+1) * pott_diff(k,i,j+2) &
                                     - CNDy2(j+1) * pott_diff(k,i,j+1) &
                                     + CNDy3(j+1) * pott_diff(k,i,j  ) &
                                     - CNDy1(j  ) * pott_diff(k,i,j-1) ) * var(k,i,j,I_DENS)
    enddo
    enddo
    enddo

    ! z-momentum
    do j = JS,   JE
    do i = IS,   IE
    do k = KS,   KE
       num_diff(k,i,j,I_MOMZ,ZDIR) = DIFF * ( 0.5D0*(CDZ(k+1)+CDZ(k)) )**DF &
                                   * ( CNMz1(k  ) * var(k+1,i,j,I_MOMZ) &
                                     - CNMz2(k  ) * var(k  ,i,j,I_MOMZ) &
                                     + CNMz3(k  ) * var(k-1,i,j,I_MOMZ) &
                                     - CNMz1(k-1) * var(k-2,i,j,I_MOMZ) )
    enddo
    enddo
    enddo
    do j = JS,   JE
    do i = IS-1, IE
    do k = WS,   WE
       num_diff(k,i,j,I_MOMZ,XDIR) = DIFF * CDX(i)**DF                  &
                                   * ( CNDx1(i+1) * var(k,i+2,j,I_MOMZ) &
                                     - CNDx2(i+1) * var(k,i+1,j,I_MOMZ) &
                                     + CNDx3(i+1) * var(k,i  ,j,I_MOMZ) &
                                     - CNDx1(i  ) * var(k,i-1,j,I_MOMZ) )
    enddo
    enddo
    enddo
    do j = JS-1, JE
    do i = IS,   IE
    do k = WS,   WE
       num_diff(k,i,j,I_MOMZ,YDIR) = DIFF * CDY(j)**DF                  &
                                   * ( CNDy1(j+1) * var(k,i,j+2,I_MOMZ) &
                                     - CNDy2(j+1) * var(k,i,j+1,I_MOMZ) &
                                     + CNDy3(j+1) * var(k,i,j  ,I_MOMZ) &
                                     - CNDy1(j  ) * var(k,i,j-1,I_MOMZ) )
    enddo
    enddo
    enddo

    ! x-momentum
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+1, WE-1
       num_diff(k,i,j,I_MOMX,ZDIR) = DIFF * CDZ(k)**DF                  &
                                   * ( CNDz1(k+1) * var(k+2,i,j,I_MOMX) &
                                     - CNDz2(k+1) * var(k+1,i,j,I_MOMX) &
                                     + CNDz3(k+1) * var(k  ,i,j,I_MOMX) &
                                     - CNDz1(k  ) * var(k-1,i,j,I_MOMX) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE+1
    do k = KS, KE
       num_diff(k,i,j,I_MOMX,XDIR) = DIFF * ( 0.5D0*(CDX(i+1)+CDX(i)) )**DF &
                                   * ( CNMx1(i  ) * var(k,i+1,j,I_MOMX) &
                                     - CNMx2(i  ) * var(k,i  ,j,I_MOMX) &
                                     + CNMx3(i  ) * var(k,i-1,j,I_MOMX) &
                                     - CNMx1(i-1) * var(k,i-2,j,I_MOMX) )
    enddo
    enddo
    enddo
    do j = JS-1, JE
    do i = IS,   IE
    do k = KS,   KE
       num_diff(k,i,j,I_MOMX,YDIR) = DIFF * CDY(j)**DF                  &
                                   * ( CNDy1(j+1) * var(k,i,j+2,I_MOMX) &
                                     - CNDy2(j+1) * var(k,i,j+1,I_MOMX) &
                                     + CNDy3(j+1) * var(k,i,j  ,I_MOMX) &
                                     - CNDy1(j  ) * var(k,i,j-1,I_MOMX) )
    enddo
    enddo
    enddo

    ! y-momentum
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+1, WE-1
       num_diff(k,i,j,I_MOMY,ZDIR) = DIFF * CDZ(k)**DF                  &
                                   * ( CNDz1(k+1) * var(k+2,i,j,I_MOMY) &
                                     - CNDz2(k+1) * var(k+1,i,j,I_MOMY) &
                                     + CNDz3(k+1) * var(k  ,i,j,I_MOMY) &
                                     - CNDz1(k  ) * var(k-1,i,j,I_MOMY) )
    enddo
    enddo
    enddo
    do j = JS,   JE
    do i = IS-1, IE
    do k = KS,   KE
       num_diff(k,i,j,I_MOMY,XDIR) = DIFF * CDX(i)**DF                  &
                                   * ( CNDx1(i+1) * var(k,i+2,j,I_MOMY) &
                                     - CNDx2(i+1) * var(k,i+1,j,I_MOMY) &
                                     + CNDx3(i+1) * var(k,i  ,j,I_MOMY) &
                                     - CNDx1(i  ) * var(k,i-1,j,I_MOMY) )
    enddo
    enddo
    enddo
    do j = JS, JE+1
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_MOMY,YDIR) = DIFF * ( 0.5D0*(CDY(j+1)+CDY(j)) )**DF &
                                   * ( CNMy1(j  ) * var(k,i,j+1,I_MOMY) &
                                     - CNMy2(j  ) * var(k,i,j  ,I_MOMY) &
                                     + CNMy3(j  ) * var(k,i,j-1,I_MOMY) &
                                     - CNMy1(j-1) * var(k,i,j-2,I_MOMY) )
    enddo
    enddo
    enddo

    call STOP_COLLECTION("prologue")
    call START_COLLECTION("RK")

    !##### Start RK #####
    do rko = 1, RK
       dtrk  = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)
       rdtrk = 1.D0 / dtrk

       if ( rko > 1 ) then
          call COMM_wait( var(:,:,:,I_DENS), I_DENS )
          call COMM_wait( var(:,:,:,I_MOMZ), I_MOMZ )
          call COMM_wait( var(:,:,:,I_MOMX), I_MOMX )
          call COMM_wait( var(:,:,:,I_MOMY), I_MOMY )
       endif

       ! momentum -> velocity
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA-1
          diagvar(k,i,j,I_VELZ) = 2.D0 * var(k,i,j,I_MOMZ) / ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) )
       enddo
       enddo
       enddo
       do j = 1, JA
       do i = 1, IA
          diagvar( 1:WS,i,j,I_VELZ) = 0.D0 ! bottom boundary
          diagvar(WE:KA,i,j,I_VELZ) = 0.D0 ! top    boundary
       enddo
       enddo

       do j = 1, JA
       do i = 1, IA-1
       do k = 1, KA
          diagvar(k,i,j,I_VELX) = 2.D0 * var(k,i,j,I_MOMX) / ( var(k,i+1,j,I_DENS)+var(k,i,j,I_DENS) )
       enddo
       enddo
       enddo
       do j = 1, JA
       do k = 1, KA
          diagvar(k,IA,j,I_VELX) = diagvar(k,IA-1,j,I_VELX)
       enddo
       enddo

       do j = 1, JA-1
       do i = 1, IA
       do k = 1, KA
          diagvar(k,i,j,I_VELY) = 2.D0 * var(k,i,j,I_MOMY) / ( var(k,i,j+1,I_DENS)+var(k,i,j,I_DENS) )
       enddo 
       enddo
       enddo
       do i = 1, IA
       do k = 1, KA
          diagvar(k,i,JA,I_VELY) = diagvar(k,i,JA-1,I_VELY)
       enddo
       enddo

       !##### continuity equation #####
       ! at (x, y, interface)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          mflx_hi(k,i,j,ZDIR) = 0.5D0 * diagvar(k,i,j,I_VELZ)                            &
                              * ( FACT_N * ( var(k+1,i,j,I_DENS)+var(k  ,i,j,I_DENS) )   &
                                + FACT_F * ( var(k+2,i,j,I_DENS)+var(k-1,i,j,I_DENS) ) ) &
                              + num_diff(k,i,j,I_DENS,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JS, JE
       do i = IS, IE
          mflx_hi(WS  ,i,j,ZDIR) = 0.D0                                            ! bottom boundary
          mflx_hi(WS+1,i,j,ZDIR) = 0.5D0 * diagvar(WS+1,i,j,I_VELZ)              &
                                 * ( var(WS+2,i,j,I_DENS)+var(WS+1,i,j,I_DENS) ) & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,I_DENS,ZDIR) * rdtrk
          mflx_hi(WE-1,i,j,ZDIR) = 0.5D0 * diagvar(WE-1,i,j,I_VELZ)              &
                                 * ( var(WE  ,i,j,I_DENS)+var(WE-1,i,j,I_DENS) ) & ! just below the top boundary
                                 + num_diff(WE-1,i,j,I_DENS,ZDIR) * rdtrk
          mflx_hi(WE  ,i,j,ZDIR) = 0.D0                                            ! top boundary
       enddo
       enddo
       ! at (u, y, layer)
       do j = JS,   JE
       do i = IS-1, IE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = 0.5D0 * diagvar(k,i,j,I_VELX)                            &
                              * ( FACT_N * ( var(k,i+1,j,I_DENS)+var(k,i  ,j,I_DENS) )   &
                                + FACT_F * ( var(k,i+2,j,I_DENS)+var(k,i-1,j,I_DENS) ) ) &
                              + num_diff(k,i,j,I_DENS,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JS-1, JE
       do i = IS,   IE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = 0.5D0 * diagvar(k,i,j,I_VELY)                            &
                              * ( FACT_N * ( var(k,i,j+1,I_DENS)+var(k,i,j  ,I_DENS) )   &
                                + FACT_F * ( var(k,i,j+2,I_DENS)+var(k,i,j-1,I_DENS) ) ) &
                              + num_diff(k,i,j,I_DENS,YDIR) * rdtrk
       enddo
       enddo
       enddo

       if ( rko == RK .AND. QA > 0 ) then
          call COMM_vars( mflx_hi(:,:,:,ZDIR),VA+ZDIR )
          call COMM_vars( mflx_hi(:,:,:,XDIR),VA+XDIR )
          call COMM_vars( mflx_hi(:,:,:,YDIR),VA+YDIR )
       endif

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JS,   JE
       do i = IS,   IE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELZ)+diagvar(k-1,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k  ,i,j,I_MOMZ)+var(k-1,i,j,I_MOMZ) )     &
                                + FACT_F * ( var(k+1,i,j,I_MOMZ)+var(k-2,i,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JS,   JE
       do i = IS,   IE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0 ! bottom cell center
          qflx_hi(KE+1,i,j,ZDIR) = 0.D0 ! top    cell center
       enddo
       enddo
       ! at (u, y, interface)
       do j = JS,   JE
       do i = IS-1, IE
       do k = WS,   WE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i+1,j,I_MOMZ)+var(k,i  ,j,I_MOMZ) )     &
                                + FACT_F * ( var(k,i+2,j,I_MOMZ)+var(k,i-1,j,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JS-1, JE
       do i = IS,   IE
       do k = WS,   WE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j+1,I_MOMZ)+var(k,i,j  ,I_MOMZ) )     &
                                + FACT_F * ( var(k,i,j+2,I_MOMZ)+var(k,i,j-1,I_MOMZ) ) )   &
                              + num_diff(k,i,j,I_MOMZ,YDIR) * rdtrk
       enddo
       enddo
       enddo

       if ( rko > 1 ) then
          call COMM_wait( var(:,:,:,I_RHOT), I_RHOT )
       endif

       ! pressure
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          diagvar(k,i,j,I_PRES) = Pstd * ( var(k,i,j,I_RHOT) * Rdry / Pstd )**CPovCV
       enddo
       enddo
       enddo

       !--- update momentum(z)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+1, WE-1
          var(k,i,j,I_MOMZ) = var_s(k,i,j,I_MOMZ) &
                            + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RFDZ(k)   &
                                         + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RCDX(i)   &
                                         + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                       - ( diagvar(k+1,i,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RFDZ(k)     & ! pressure gradient force
                                       - ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) ) * 0.5D0 * GRAV        & ! gravity force
                                       + ray_damp(k,i,j,I_MOMZ)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       !--- update density
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,I_DENS) = var_s(k,i,j,I_DENS) &
                            + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k)   &
                                         + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i)   &
                                         + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) ) ! divergence
       enddo
       enddo
       enddo

       call COMM_vars( var(:,:,:,I_DENS),I_DENS )
       call COMM_vars( var(:,:,:,I_MOMZ),I_MOMZ )

       !##### momentum equation (x) #####
       ! at (x, y, layer)
       do j = JS, JE
       do i = IS, IE+1
       do k = KS, KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELX)+diagvar(k,i-1,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i  ,j,I_MOMX)+var(k,i-1,j,I_MOMX) )     &
                                + FACT_F * ( var(k,i+1,j,I_MOMX)+var(k,i-2,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do j = JS-1, JE
       do i = IS,   IE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j+1,I_MOMX)+var(k,i,j  ,I_MOMX) )     &
                                + FACT_F * ( var(k,i,j+2,I_MOMX)+var(k,i,j-1,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,YDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k+1,i,j,I_MOMX)+var(k  ,i,j,I_MOMX) )     &
                                + FACT_F * ( var(k+2,i,j,I_MOMX)+var(k-1,i,j,I_MOMX) ) )   &
                              + num_diff(k,i,j,I_MOMX,ZDIR) * rdtrk
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_hi(WS  ,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(WS+1,i,j,ZDIR) = 0.25D0 * ( diagvar(WS+1,i+1,j,I_VELZ)+diagvar(WS+1,i,j,I_VELZ) ) &
                                 * ( var(WS+2,i,j,I_MOMX)+var(WS+1,i,j,I_MOMX) )                    & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(WE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(WE-1,i+1,j,I_VELZ)+diagvar(WE-1,i,j,I_VELZ) ) &
                                 * ( var(WE  ,i,j,I_MOMX)+var(WE-1,i,j,I_MOMX) )                    & ! just below the top boundary
                                 + num_diff(WE-1,i,j,I_MOMX,ZDIR) * rdtrk
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       !--- update momentum(x)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,I_MOMX) = var_s(k,i,j,I_MOMX) &
                            + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RCDZ(k)   &
                                         + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RFDX(i)   &
                                         + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) & ! flux divergence
                                       - ( diagvar(k,i+1,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RFDX(i)     & ! pressure gradient force
                                       + ray_damp(k,i,j,I_MOMX)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       call COMM_vars( var(:,:,:,I_MOMX),I_MOMX )

       !##### momentum equation (y) #####
       ! at (u, v, layer)
       do j = JS,   JE
       do i = IS-1, IE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i+1,j,I_MOMY)+var(k,i  ,j,I_MOMY) )     &
                                + FACT_F * ( var(k,i+2,j,I_MOMY)+var(k,i-1,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, y, layer)
       do j = JS, JE+1
       do i = IS, IE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELY)+diagvar(k,i,j-1,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j  ,I_MOMY)+var(k,i,j-1,I_MOMY) )     &
                                + FACT_F * ( var(k,i,j+1,I_MOMY)+var(k,i,j-2,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,YDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k+1,i,j,I_MOMY)+var(k  ,i,j,I_MOMY) )     &
                                + FACT_F * ( var(k+2,i,j,I_MOMY)+var(k-1,i,j,I_MOMY) ) )   &
                              + num_diff(k,i,j,I_MOMY,ZDIR) * rdtrk
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_hi(WS  ,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(WS+1,i,j,ZDIR) = 0.25D0 * ( diagvar(WS+1,i,j+1,I_VELZ)+diagvar(WS+1,i,j,I_VELZ) ) &
                                 * ( var(WS+2,i,j,I_MOMY)+var(WS+1,i,j,I_MOMY) )                    & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(WE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(WE-1,i,j+1,I_VELZ)+diagvar(WE-1,i,j,I_VELZ) ) &
                                 * ( var(WE  ,i,j,I_MOMY)+var(WE-1,i,j,I_MOMY) )                    & ! just below the top boundary
                                 + num_diff(WE-1,i,j,I_MOMY,ZDIR) * rdtrk
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

       !--- update momentum(y)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,I_MOMY) = var_s(k,i,j,I_MOMY) &
                            + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k)   &
                                         + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i)   &
                                         + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) & ! flux divergence
                                       - ( diagvar(k,i,j+1,I_PRES)-diagvar(k,i,j,I_PRES) ) * RFDY(j)     & ! pressure gradient force
                                       + ray_damp(k,i,j,I_MOMY)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       call COMM_vars( var(:,:,:,I_MOMY),I_MOMY )

       !##### Thermodynamic Equation #####

       do j = 1, JA !JS, JE
       do i = 1, IA !IS, IE
       do k = 1, KA !KS, KE
          diagvar(k,i,j,I_POTT) = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS) 
       enddo
       enddo
       enddo

       ! at (u, y, layer)
       do j = JS,   JE
       do i = IS-1, IE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i+1,j,I_POTT)+diagvar(k,i  ,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i+2,j,I_POTT)+diagvar(k,i-1,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JS-1, JE
       do i = IS,   IE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                                      &
                              * ( FACT_N * ( diagvar(k,i,j+1,I_POTT)+diagvar(k,i,j  ,I_POTT) )   &
                                + FACT_F * ( diagvar(k,i,j+2,I_POTT)+diagvar(k,i,j-1,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,YDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, y, interface)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                                      &
                              * ( FACT_N * ( diagvar(k+1,i,j,I_POTT)+diagvar(k  ,i,j,I_POTT) )   &
                                + FACT_F * ( diagvar(k+2,i,j,I_POTT)+diagvar(k-1,i,j,I_POTT) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JS, JE
       do i = IS, IE
          qflx_hi(WS  ,i,j,ZDIR) = 0.D0                                                    ! bottom boundary
          qflx_hi(WS+1,i,j,ZDIR) = 0.5D0 * mflx_hi(WS+1,i,j,ZDIR)                        &
                                 * ( diagvar(WS+2,i,j,I_POTT)+diagvar(WS+1,i,j,I_POTT) ) & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(WE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(WE-1,i,j,ZDIR)                        &
                                 * ( diagvar(WE  ,i,j,I_POTT)+diagvar(WE-1,i,j,I_POTT) ) & ! just below the top boundary
                                 + num_diff(WE-1,i,j,I_RHOT,ZDIR) * rdtrk
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0                                                    ! top boundary
       enddo
       enddo

       !--- update rho*theta
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,I_RHOT) = var_s(k,i,j,I_RHOT) &
                            + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k)   &
                                         + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i)   &
                                         + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) & ! divergence
                                       + ray_damp(k,i,j,I_RHOT)                                          ) ! additional damping force
       enddo
       enddo
       enddo

       call COMM_vars( var(:,:,:,I_RHOT),I_RHOT )

    enddo ! RK loop

    call STOP_COLLECTION("RK")
    call START_COLLECTION("FCT")

    !##### advection of scalar quantity #####

    if ( QA > 0 ) then

    ! calc low-order mass flux and high-low difference (at last step of RK)
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+2, WE-2
       mflx_lo(k,i,j,ZDIR) = 0.5D0 * (     diagvar(k,i,j,I_VELZ)  * ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) ) &
                                     - abs(diagvar(k,i,j,I_VELZ)) * ( var(k+1,i,j,I_DENS)-var(k,i,j,I_DENS) ) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       mflx_lo(WS  ,i,j,ZDIR) = 0.D0                                                                         ! bottom boundary
       mflx_lo(WS+1,i,j,ZDIR) = 0.5D0 * diagvar(WS+1,i,j,I_VELZ) * ( var(WS+2,i,j,I_DENS)+var(WS+1,i,j,I_DENS) ) ! just above the bottom boundary
       mflx_lo(WE-1,i,j,ZDIR) = 0.5D0 * diagvar(WE-1,i,j,I_VELZ) * ( var(WE  ,i,j,I_DENS)+var(WE-1,i,j,I_DENS) ) ! just below the top boundary
       mflx_lo(WE  ,i,j,ZDIR) = 0.D0                                                                         ! top boundary 
    enddo
    enddo

    do j = JS,   JE
    do i = IS-1, IE
    do k = KS,   KE
       mflx_lo(k,i,j,XDIR) = 0.5D0 * (     diagvar(k,i,j,I_VELX)  * ( var(k,i+1,j,I_DENS)+var(k,i,j,I_DENS) ) &
                                     - abs(diagvar(k,i,j,I_VELX)) * ( var(k,i+1,j,I_DENS)-var(k,i,j,I_DENS) ) )
    enddo
    enddo
    enddo

    do j = JS-1, JE
    do i = IS,   IE
    do k = KS,   KE
       mflx_lo(k,i,j,YDIR) = 0.5D0 * (     diagvar(k,i,j,I_VELY)  * ( var(k,i,j+1,I_DENS)+var(k,i,j,I_DENS) ) &
                                     - abs(diagvar(k,i,j,I_VELY)) * ( var(k,i,j+1,I_DENS)-var(k,i,j,I_DENS) ) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ddiv(k,i,j) = ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                   + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                   + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j)
    enddo
    enddo
    enddo

    endif

    call COMM_wait( var(:,:,:,I_DENS), I_DENS )
    call COMM_wait( var(:,:,:,I_MOMZ), I_MOMZ )
    call COMM_wait( var(:,:,:,I_MOMX), I_MOMX )
    call COMM_wait( var(:,:,:,I_MOMY), I_MOMY )

    if ( QA > 0 ) then

    call COMM_wait( var(:,:,:,VA+ZDIR), VA+ZDIR )
    call COMM_wait( var(:,:,:,VA+XDIR), VA+XDIR )
    call COMM_wait( var(:,:,:,VA+YDIR), VA+YDIR )

    do iq = 6, 5+QA

       call COMM_wait( var(:,:,:,iq-1), iq-1 )

       do j = JS,   JE
       do i = IS-1, IE
       do k = KS,   KE
          qflx_lo(k,i,j,XDIR) = 0.5D0                                                          &
                              * (     mflx_lo(k,i,j,XDIR)  * ( var(k,i+1,j,iq)+var(k,i,j,iq) ) &
                                - abs(mflx_lo(k,i,j,XDIR)) * ( var(k,i+1,j,iq)-var(k,i,j,iq) ) )

          qflx_hi(k,i,j,XDIR) = 0.5D0 * mflx_hi(k,i,j,XDIR)                    &
                              * ( FACT_N * ( var(k,i+1,j,iq)+var(k,i  ,j,iq) ) &
                                + FACT_F * ( var(k,i+2,j,iq)+var(k,i-1,j,iq) ) )

          qflx_anti(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_lo(k,i,j,XDIR)
       enddo
       enddo
       enddo

       do j = JS-1, JE
       do i = IS,   IE
       do k = KS,   KE
          qflx_lo(k,i,j,YDIR) = 0.5D0                                                          &
                              * (     mflx_lo(k,i,j,YDIR)  * ( var(k,i,j+1,iq)+var(k,i,j,iq) ) &
                                - abs(mflx_lo(k,i,j,YDIR)) * ( var(k,i,j+1,iq)-var(k,i,j,iq) ) )

          qflx_hi(k,i,j,YDIR) = 0.5D0 * mflx_hi(k,i,j,YDIR)                    &
                              * ( FACT_N * ( var(k,i,j+1,iq)+var(k,i,j  ,iq) ) &
                                + FACT_F * ( var(k,i,j+2,iq)+var(k,i,j-1,iq) ) )

          qflx_anti(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) - qflx_lo(k,i,j,YDIR)
       enddo
       enddo
       enddo

       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          qflx_lo(k,i,j,ZDIR) = 0.5D0                                                          &
                              * (     mflx_lo(k,i,j,ZDIR)  * ( var(k+1,i,j,iq)+var(k,i,j,iq) ) &
                                - abs(mflx_lo(k,i,j,ZDIR)) * ( var(k+1,i,j,iq)-var(k,i,j,iq) ) )

          qflx_hi(k,i,j,ZDIR) = 0.5D0 * mflx_hi(k,i,j,ZDIR)                    &
                              * ( FACT_N * ( var(k+1,i,j,iq)+var(k  ,i,j,iq) ) &
                                + FACT_F * ( var(k+2,i,j,iq)+var(k-1,i,j,iq) ) )

          qflx_anti(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_lo(k,i,j,ZDIR)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_lo(WS  ,i,j,ZDIR) = 0.D0                                                                   ! bottom boundary
          qflx_lo(WS+1,i,j,ZDIR) = 0.5D0 * mflx_lo(WS+1,i,j,ZDIR) * ( var(WS+2,i,j,iq)+var(WS+1,i,j,iq) ) ! just above the bottom boundary
          qflx_lo(WE-1,i,j,ZDIR) = 0.5D0 * mflx_lo(WE-1,i,j,ZDIR) * ( var(WE  ,i,j,iq)+var(WE-1,i,j,iq) ) ! just below the top boundary
          qflx_lo(WE  ,i,j,ZDIR) = 0.D0                                                                   ! top boundary

          qflx_hi(WS  ,i,j,ZDIR) = 0.D0 
          qflx_hi(WS+1,i,j,ZDIR) = qflx_lo(WS+1,i,j,ZDIR)
          qflx_hi(WE-1,i,j,ZDIR) = qflx_lo(WE-1,i,j,ZDIR)
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0 

          qflx_anti(WS  ,i,j,ZDIR) = 0.D0
          qflx_anti(WS+1,i,j,ZDIR) = 0.D0
          qflx_anti(WE-1,i,j,ZDIR) = 0.D0
          qflx_anti(WE  ,i,j,ZDIR) = 0.D0
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          !--- update value with flux-divergence from the monotone scheme
          var(k,i,j,iq) = var_s(k,i,j,iq) &
                        + dtrk * ( - ( ( qflx_lo(k,i,j,XDIR)-qflx_lo(k  ,i-1,j,  XDIR) ) * RCDX(i)   &
                                     + ( qflx_lo(k,i,j,YDIR)-qflx_lo(k  ,i,  j-1,YDIR) ) * RCDY(j)   &
                                     + ( qflx_lo(k,i,j,ZDIR)-qflx_lo(k-1,i,  j,  ZDIR) ) * RCDZ(k) ) &
                                   + var_s(k,i,j,iq) * ddiv(k,i,j)                                   &
                                 ) / var(k,i,j,I_DENS)

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( var_s(k  ,i  ,j  ,iq), &
                       var_s(k  ,i+1,j  ,iq), &
                       var_s(k  ,i-1,j  ,iq), &
                       var_s(k  ,i  ,j+1,iq), &
                       var_s(k  ,i  ,j-1,iq), &
                       var_s(k+1,i  ,j  ,iq), &
                       var_s(k-1,i  ,j  ,iq)  )

          pjmin = min( var_s(k  ,i  ,j  ,iq), &
                       var_s(k  ,i+1,j  ,iq), &
                       var_s(k  ,i-1,j  ,iq), &
                       var_s(k  ,i  ,j+1,iq), &
                       var_s(k  ,i  ,j-1,iq), &
                       var_s(k+1,i  ,j  ,iq), &
                       var_s(k-1,i  ,j  ,iq)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(k  ,i-1,j  ,XDIR) ) - min( 0.D0, qflx_anti(k  ,i  ,j  ,XDIR) ) &
                + max( 0.D0, qflx_anti(k  ,i  ,j-1,YDIR) ) - min( 0.D0, qflx_anti(k  ,i  ,j  ,YDIR) ) &
                + max( 0.D0, qflx_anti(k-1,i  ,j  ,ZDIR) ) - min( 0.D0, qflx_anti(k  ,i  ,j  ,ZDIR) )
          pjmns = max( 0.D0, qflx_anti(k  ,i  ,j  ,XDIR) ) - min( 0.D0, qflx_anti(k  ,i-1,j  ,XDIR) ) &
                + max( 0.D0, qflx_anti(k  ,i  ,j  ,YDIR) ) - min( 0.D0, qflx_anti(k  ,i  ,j-1,YDIR) ) &
                + max( 0.D0, qflx_anti(k  ,i  ,j  ,ZDIR) ) - min( 0.D0, qflx_anti(k-1,i  ,j  ,ZDIR) )
          ! --- incoming fluxes ---
          if ( pjpls > 0 ) then
             rjpls(k,i,j,XDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(k,i,j,XDIR)+mflx_lo(k  ,i-1,j  ,XDIR)) * 0.5D0)
             rjpls(k,i,j,YDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(k,i,j,YDIR)+mflx_lo(k  ,i  ,j-1,YDIR)) * 0.5D0)
             rjpls(k,i,j,ZDIR) = (pjmax-var_l) / pjpls * abs((mflx_lo(k,i,j,ZDIR)+mflx_lo(k-1,i  ,j  ,ZDIR)) * 0.5D0)
          else
             rjpls(k,i,j,XDIR:ZDIR) = 0.D0
          endif
          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             rjmns(k,i,j,XDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(k,i,j,XDIR)+mflx_lo(k  ,i-1,j  ,XDIR)) * 0.5D0)
             rjmns(k,i,j,YDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(k,i,j,YDIR)+mflx_lo(k  ,i  ,j-1,YDIR)) * 0.5D0)
             rjmns(k,i,j,ZDIR) = (var_l-pjmin) / pjmns * abs((mflx_lo(k,i,j,ZDIR)+mflx_lo(k-1,i  ,j  ,ZDIR)) * 0.5D0)
          else
             rjmns(k,i,j,XDIR:ZDIR) = 0.D0
          endif
       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux ---
       do j = JS-1, JE 
       do i = IS-1, IE
       do k = KS-1, KE 
          if ( qflx_anti(k,i,j,XDIR) >= 0 ) then
             qflx_anti(k,i,j,XDIR) = qflx_anti(k,i,j,XDIR) * min( rjpls(k,i+1,j,XDIR), rjmns(k,i  ,j,XDIR), 1.D0 )
          else
             qflx_anti(k,i,j,XDIR) = qflx_anti(k,i,j,XDIR) * min( rjpls(k,i  ,j,XDIR), rjmns(k,i+1,j,XDIR), 1.D0 )
          endif
          if ( qflx_anti(k,i,j,YDIR) >= 0 ) then
             qflx_anti(k,i,j,YDIR) = qflx_anti(k,i,j,YDIR) * min( rjpls(k,i,j+1,YDIR), rjmns(k,i,j  ,YDIR), 1.D0 )
          else
             qflx_anti(k,i,j,YDIR) = qflx_anti(k,i,j,YDIR) * min( rjpls(k,i,j  ,YDIR), rjmns(k,i,j+1,YDIR), 1.D0 )
          endif
          if ( qflx_anti(k,i,j,ZDIR) >= 0 ) then
             qflx_anti(k,i,j,ZDIR) = qflx_anti(k,i,j,ZDIR) * min( rjpls(k+1,i,j,ZDIR), rjmns(k  ,i,j,ZDIR), 1.D0 )
          else
             qflx_anti(k,i,j,ZDIR) = qflx_anti(k,i,j,ZDIR) * min( rjpls(k  ,i,j,ZDIR), rjmns(k+1,i,j,ZDIR), 1.D0 )
          endif
       enddo
       enddo
       enddo

       !--- modify value with antidiffusive fluxes
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,iq) = var(k,i,j,iq) &
                        + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR)-qflx_anti(k-1,i,  j,  ZDIR) ) * RCDZ(k)   &
                                     + ( qflx_anti(k,i,j,XDIR)-qflx_anti(k  ,i-1,j,  XDIR) ) * RCDX(i)   &
                                     + ( qflx_anti(k,i,j,YDIR)-qflx_anti(k  ,i,  j-1,YDIR) ) * RCDY(j) ) &
                                 ) / var(k,i,j,I_DENS)
       enddo
       enddo
       enddo

       call COMM_vars( var(:,:,:,iq),iq )

    enddo ! scalar quantities loop

    call COMM_wait( var(:,:,:,iq-1), iq-1 )

    else

    call COMM_wait( var(:,:,:,I_RHOT), I_RHOT )

    endif

    call STOP_COLLECTION("FCT")

    ! check total mass
    call COMM_total( var(:,:,:,1:5), A_NAME(1:5) )

    enddo ! dynamical steps

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
