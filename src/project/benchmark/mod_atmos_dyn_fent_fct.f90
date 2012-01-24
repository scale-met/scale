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
  integer, parameter :: I_PRES = 1
  integer, parameter :: I_VELX = 2
  integer, parameter :: I_VELY = 3
  integer, parameter :: I_VELZ = 4
  integer, parameter :: I_POTT = 5

  integer, parameter :: XDIR   = 1
  integer, parameter :: YDIR   = 2
  integer, parameter :: ZDIR   = 3

  ! time settings
  integer, parameter :: RK = 3 ! order of Runge-Kutta scheme

  ! advection settings
  real(8), parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second
  real(8), parameter :: FACT_M =   9.D0 / 8.D0 !  7/6: fourth, 1: second
  real(8), parameter :: FACT_G = - 1.D0 / 8.D0 ! -1/6: fourth, 0: second

  ! numerical filter settings
  integer, parameter :: DF   = 4     ! order of numerical filter

  real(8), save      :: ATMOS_DYN_numerical_diff = 1.D-3 ! nondimensional numerical diffusion
  real(8), save      :: DIFF

  real(8), allocatable, save :: CNDX(:,:)
  real(8), allocatable, save :: CNMX(:,:)
  real(8), allocatable, save :: CNDY(:,:)
  real(8), allocatable, save :: CNMY(:,:)
  real(8), allocatable, save :: CNDZ(:,:)
  real(8), allocatable, save :: CNMZ(:,:)

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
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KA  => GRID_KA, &
       IS  => GRID_IS, &
       IE  => GRID_IE, &
       JS  => GRID_JS, &
       JE  => GRID_JE, &
       KS  => GRID_KS, &
       KE  => GRID_KE, &
       WS  => GRID_WS, &
       WE  => GRID_WE, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY, &
       CDZ => GRID_CDZ
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
    CNDZ(1,1)  = CNDZ(1,KS)
    CNDZ(2,1)  = CNDZ(2,KS)
    CNDZ(1,KA) = CNDZ(1,KE)
    CNDZ(2,KA) = CNDZ(2,KE)

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
    CNMZ(1,KA) = CNMZ(1,KE)

    do k = KS-1, KE+1
       CNMZ(2,k) = 1.D0 / ( CDZ(k+1) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &   
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &  
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) 
       CNMZ(3,k) = 1.D0 / ( CDZ(k  ) * (CDZ(k+1)+CDZ(k  )) * 0.5D0 * CDZ(k  ) ) &
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k  ) ) &
                 + 1.D0 / ( CDZ(k  ) * (CDZ(k  )+CDZ(k-1)) * 0.5D0 * CDZ(k-1) )
    enddo
    CNMZ(2,1)  = CNMZ(2,KS)
    CNMZ(3,1)  = CNMZ(3,KS)
    CNMZ(2,KA) = CNMZ(2,KE)
    CNMZ(3,KA) = CNMZ(3,KE)

    ! x direction
    do i = IS-1, IE+1
       CNDX(1,i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
       CNDX(2,i) = 1.D0 / ( (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) * (CDX(i  )+CDX(i-1)) * 0.5D0 )
    enddo
    CNDX(1,1)  = CNDX(1,IS)
    CNDX(2,1)  = CNDX(2,IS)
    CNDX(1,IA) = CNDX(1,IE)
    CNDX(2,IA) = CNDX(2,IE)

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
    CNMX(1,IA) = CNMX(1,IE)

    do i = IS-1, IE+1
       CNMX(2,i) = 1.D0 / ( CDX(i+1) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) )
       CNMX(3,i) = 1.D0 / ( CDX(i  ) * (CDX(i+1)+CDX(i  )) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i  ) ) & 
                 + 1.D0 / ( CDX(i  ) * (CDX(i  )+CDX(i-1)) * 0.5D0 * CDX(i-1) )
    enddo
    CNMX(2,1)  = CNMX(2,IS)
    CNMX(3,1)  = CNMX(3,IS)
    CNMX(2,IA) = CNMX(2,IE)
    CNMX(3,IA) = CNMX(3,IE)

    ! y direction
    do j = JS, JE
       CNDY(1,j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
       CNDY(2,j) = 1.D0 / ( (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 )
       CNDY(3,j) = 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j  )+CDY(j-1)) * 0.5D0 ) &
                 + 1.D0 / ( (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) * (CDY(j-1)+CDY(j-2)) * 0.5D0 )
       CNMY(1,j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) )
       CNMY(2,j) = 1.D0 / ( CDY(j+1) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &   
                 + 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &  
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) 
       CNMY(3,j) = 1.D0 / ( CDY(j  ) * (CDY(j+1)+CDY(j  )) * 0.5D0 * CDY(j  ) ) &
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j  ) ) &
                 + 1.D0 / ( CDY(j  ) * (CDY(j  )+CDY(j-1)) * 0.5D0 * CDY(j-1) )
    enddo
    do j = 1, JS-1
       CNDY(1,j) = CNDY(1,JS)
       CNDY(2,j) = CNDY(2,JS)
       CNDY(3,j) = CNDY(3,JS)
       CNMY(1,j) = CNMY(1,JS)
       CNMY(2,j) = CNMY(2,JS)
       CNMY(3,j) = CNMY(3,JS)
    enddo
    do j = JE+1, JA
       CNDY(1,j) = CNDY(1,JE)
       CNDY(2,j) = CNDY(2,JE)
       CNDY(3,j) = CNDY(3,JE)
       CNMY(1,j) = CNMY(1,JE)
       CNMY(2,j) = CNMY(2,JE)
       CNMY(3,j) = CNMY(3,JE)
    enddo

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
       COMM_wait, &
       COMM_total
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       WS   => GRID_WS,   &
       WE   => GRID_WE,   &
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
       A_NAME,      &
       VA  => A_VA, &
       QA  => A_QA
    use mod_atmos_refstate, only: &
       REF_dens => ATMOS_REFSTATE_dens, &
       REF_pott => ATMOS_REFSTATE_pott
    use mod_atmos_boundary, only: &
       DAMP_var   => ATMOS_BOUNDARY_var,   &
       DAMP_alpha => ATMOS_BOUNDARY_alpha
    implicit none

    ! work
    real(8) :: var_s  (KA,IA,JA,VA)            ! prognostic variables (previous step)
    real(8) :: diagvar(KA,IA,JA,I_PRES:I_POTT) ! diagnostic variables (work)

    real(8) :: dens_diff(KA,IA,JA)             ! anomary of density
    real(8) :: pott_diff(KA,IA,JA)             ! anomary of rho * pott
    real(8) :: ray_damp(KA,IA,JA,1:5)
    real(8) :: num_diff(KA,IA,JA,1:5,XDIR:ZDIR)

    ! mass flux
    real(8) :: mflx_hi  (KA,IA,JA,XDIR:ZDIR) ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8) :: mflx_lo  (KA,IA,JA,XDIR:ZDIR) ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8) :: qflx_hi  (KA,IA,JA,XDIR:ZDIR) ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(8) :: qflx_lo  (KA,IA,JA,XDIR:ZDIR) ! rho * vel(x,y,z) * phi @ (u,v,w)-face low  order
    real(8) :: qflx_anti(KA,IA,JA,XDIR:ZDIR) ! rho * vel(x,y,z) * phi @ (u,v,w)-face antidiffusive
    real(8) :: ddiv     (KA,IA,JA)           ! density divergence
    ! flux correction
    real(8) :: rjpls    (KA,IA,JA,XDIR:ZDIR) ! plus  in (x,y,z)-direction
    real(8) :: rjmns    (KA,IA,JA,XDIR:ZDIR) ! minus in (x,y,z)-direction
    real(8) :: pjpls, pjmns, pjmax, pjmin

    real(8) :: dtrk, rdtrk
    integer :: i, j, k, iq, iv, rko, step
    !---------------------------------------------------------------------------

    call fpcoll_start

    do step = 1, TIME_NSTEP_ATMOS_DYN

    if( IO_L ) write(IO_FID_LOG,*) '*** Dynamical small step:', step

call START_COLLECTION("SET")

call START_COLLECTION("PRL1")
    do iv = 1, VA
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       var_s(k,i,j,iv) = var(k,i,j,iv)
    enddo 
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL1")

call START_COLLECTION("PRL2")
    !--- prepare rayleigh damping coefficient
    do j = JS, JE
    do i = IS, IE
       do k = WS+1, WE-1
          ray_damp(k,i,j,4) = - DAMP_alpha(k,i,j,1) * ( var(k,i,j,4) &
                              - DAMP_var(k,i,j,1) * 0.5D0 * ( var(k+1,i,j,1)+var(k,i,j,1) ) )
       enddo
       do k = KS, KE
          ray_damp(k,i,j,2) = - DAMP_alpha(k,i,j,2) * ( var(k,i,j,2) &
                              - DAMP_var(k,i,j,2) * 0.5D0 * ( var(k,i+1,j,1)+var(k,i,j,1) ) )
          ray_damp(k,i,j,3) = - DAMP_alpha(k,i,j,3) * ( var(k,i,j,3) &
                              - DAMP_var(k,i,j,3) * 0.5D0 * ( var(k,i,j+1,1)+var(k,i,j,1) ) )
          ray_damp(k,i,j,5) = - DAMP_alpha(k,i,j,4) * ( var(k,i,j,5) &
                              - DAMP_var(k,i,j,4) * var(k,i,j,1) )
       enddo 
    enddo
    enddo
call STOP_COLLECTION("PRL2")

    !--- prepare numerical diffusion coefficient
    ! note: the value must divide by dtrk.
call START_COLLECTION("PRL3")
    do j  = 1, JA
    do i  = 1, IA
       do k = KS, KE
          dens_diff(k,i,j) = var(k,i,j,1) - REF_dens(k)
          pott_diff(k,i,j) = var(k,i,j,5) / var(k,i,j,1) - REF_pott(k)
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
call STOP_COLLECTION("PRL3")

call START_COLLECTION("PRL4")
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+1, WE-1
       num_diff(k,i,j,1,ZDIR) = DIFF * CDZ(k)**DF                  &
                              * ( CNDZ(1,k+1) * dens_diff(k+2,i,j) &
                                - CNDZ(2,k+1) * dens_diff(k+1,i,j) &
                                + CNDZ(3,k+1) * dens_diff(k  ,i,j) &
                                - CNDZ(1,k  ) * dens_diff(k-1,i,j) )
       num_diff(k,i,j,5,ZDIR) = DIFF * CDZ(k)**DF                  &
                              * ( CNDZ(1,k+1) * pott_diff(k+2,i,j) &
                                - CNDZ(2,k+1) * pott_diff(k+1,i,j) &
                                + CNDZ(3,k+1) * pott_diff(k  ,i,j) &
                                - CNDZ(1,k  ) * pott_diff(k-1,i,j) ) * var(k,i,j,1)
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL4")

call START_COLLECTION("PRL5")
    do j = JS,   JE
    do i = IS-1, IE
    do k = KS,   KE
       num_diff(k,i,j,1,XDIR) = DIFF * CDX(i)**DF                  &
                              * ( CNDX(1,i+1) * dens_diff(k,i+2,j) &
                                - CNDX(2,i+1) * dens_diff(k,i+1,j) &
                                + CNDX(3,i+1) * dens_diff(k,i  ,j) &
                                - CNDX(1,i  ) * dens_diff(k,i-1,j) )
       num_diff(k,i,j,5,XDIR) = DIFF * CDX(i)**DF                  &
                              * ( CNDX(1,i+1) * pott_diff(k,i+2,j) &
                                - CNDX(2,i+1) * pott_diff(k,i+1,j) &
                                + CNDX(3,i+1) * pott_diff(k,i  ,j) &
                                - CNDX(1,i  ) * pott_diff(k,i-1,j) ) * var(k,i,j,1)

    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL5")

call START_COLLECTION("PRL6")
    do j = JS-1, JE
    do i = IS,   IE
    do k = KS,   KE
       num_diff(k,i,j,1,YDIR) = DIFF * CDY(j)**DF                  &
                              * ( CNDY(1,j+1) * dens_diff(k,i,j+2) &
                                - CNDY(2,j+1) * dens_diff(k,i,j+1) &
                                + CNDY(3,j+1) * dens_diff(k,i,j  ) &
                                - CNDY(1,j  ) * dens_diff(k,i,j-1) )
       num_diff(k,i,j,5,YDIR) = DIFF * CDY(j)**DF                  &
                              * ( CNDY(1,j+1) * pott_diff(k,i,j+2) &
                                - CNDY(2,j+1) * pott_diff(k,i,j+1) &
                                + CNDY(3,j+1) * pott_diff(k,i,j  ) &
                                - CNDY(1,j  ) * pott_diff(k,i,j-1) ) * var(k,i,j,1)
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL6")

call START_COLLECTION("PRL7")
    ! z-momentum
    do j = JS,   JE
    do i = IS,   IE
    do k = KS,   KE
       num_diff(k,i,j,4,ZDIR) = DIFF * ( 0.5D0*(CDZ(k+1)+CDZ(k)) )**DF &
                              * ( CNMZ(1,k  ) * var(k+1,i,j,4) &
                                - CNMZ(2,k  ) * var(k  ,i,j,4) &
                                + CNMZ(3,k  ) * var(k-1,i,j,4) &
                                - CNMZ(1,k-1) * var(k-2,i,j,4) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL7")
call START_COLLECTION("PRL8")
    do j = JS,   JE
    do i = IS-1, IE
    do k = WS,   WE
       num_diff(k,i,j,4,XDIR) = DIFF * CDX(i)**DF                   &
                              * ( CNDX(1,i+1) * var(k,i+2,j,4) &
                                - CNDX(2,i+1) * var(k,i+1,j,4) &
                                + CNDX(3,i+1) * var(k,i  ,j,4) &
                                - CNDX(1,i  ) * var(k,i-1,j,4) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL8")
call START_COLLECTION("PRL9")
    do j = JS-1, JE
    do i = IS,   IE
    do k = WS,   WE
       num_diff(k,i,j,4,YDIR) = DIFF * CDY(j)**DF              &
                              * ( CNDY(1,j+1) * var(k,i,j+2,4) &
                                - CNDY(2,j+1) * var(k,i,j+1,4) &
                                + CNDY(3,j+1) * var(k,i,j  ,4) &
                                - CNDY(1,j  ) * var(k,i,j-1,4) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL9")

call START_COLLECTION("PRL10")
    ! x-momentum
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+1, WE-1
       num_diff(k,i,j,2,ZDIR) = DIFF * CDZ(k)**DF              &
                              * ( CNDZ(1,k+1) * var(k+2,i,j,2) &
                                - CNDZ(2,k+1) * var(k+1,i,j,2) &
                                + CNDZ(3,k+1) * var(k  ,i,j,2) &
                                - CNDZ(1,k  ) * var(k-1,i,j,2) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL10")
call START_COLLECTION("PRL11")
    do j = JS, JE
    do i = IS, IE+1
    do k = KS, KE
       num_diff(k,i,j,2,XDIR) = DIFF * ( 0.5D0*(CDX(i+1)+CDX(i)) )**DF &
                              * ( CNMX(1,i  ) * var(k,i+1,j,2) &
                                - CNMX(2,i  ) * var(k,i  ,j,2) &
                                + CNMX(3,i  ) * var(k,i-1,j,2) &
                                - CNMX(1,i-1) * var(k,i-2,j,2) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL11")
call START_COLLECTION("PRL12")
    do j = JS-1, JE
    do i = IS,   IE
    do k = KS,   KE
       num_diff(k,i,j,2,YDIR) = DIFF * CDY(j)**DF                  &
                                   * ( CNDY(1,j+1) * var(k,i,j+2,2) &
                                     - CNDY(2,j+1) * var(k,i,j+1,2) &
                                     + CNDY(3,j+1) * var(k,i,j  ,2) &
                                     - CNDY(1,j  ) * var(k,i,j-1,2) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL12")

call START_COLLECTION("PRL13")
    ! y-momentum
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+1, WE-1
       num_diff(k,i,j,3,ZDIR) = DIFF * CDZ(k)**DF                  &
                                   * ( CNDZ(1,k+1) * var(k+2,i,j,3) &
                                     - CNDZ(2,k+1) * var(k+1,i,j,3) &
                                     + CNDZ(3,k+1) * var(k  ,i,j,3) &
                                     - CNDZ(1,k  ) * var(k-1,i,j,3) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL13")
call START_COLLECTION("PRL14")
    do j = JS,   JE
    do i = IS-1, IE
    do k = KS,   KE
       num_diff(k,i,j,3,XDIR) = DIFF * CDX(i)**DF                  &
                                   * ( CNDX(1,i+1) * var(k,i+2,j,3) &
                                     - CNDX(2,i+1) * var(k,i+1,j,3) &
                                     + CNDX(3,i+1) * var(k,i  ,j,3) &
                                     - CNDX(1,i  ) * var(k,i-1,j,3) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL14")
call START_COLLECTION("PRL15")
    do j = JS, JE+1
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,3,YDIR) = DIFF * ( 0.5D0*(CDY(j+1)+CDY(j)) )**DF &
                                   * ( CNMY(1,j  ) * var(k,i,j+1,3) &
                                     - CNMY(2,j  ) * var(k,i,j  ,3) &
                                     + CNMY(3,j  ) * var(k,i,j-1,3) &
                                     - CNMY(1,j-1) * var(k,i,j-2,3) )
    enddo
    enddo
    enddo
call STOP_COLLECTION("PRL15")

call STOP_COLLECTION("SET")
call START_COLLECTION("RK3")

    !##### Start RK #####
    do rko = 1, RK
       dtrk  = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)
       rdtrk = 1.D0 / dtrk

call START_COLLECTION("PRL16")
       ! momentum -> velocity
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA-1
          diagvar(k,i,j,I_VELZ) = 2.D0 * var(k,i,j,4) / ( var(k+1,i,j,1)+var(k,i,j,1) )
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL16")
call START_COLLECTION("PRL17")
       do j = 1, JA
       do i = 1, IA
          diagvar( 1:WS,i,j,I_VELZ) = 0.D0 ! bottom boundary
          diagvar(WE:KA,i,j,I_VELZ) = 0.D0 ! top    boundary
       enddo
       enddo
call STOP_COLLECTION("PRL17")

call START_COLLECTION("PRL18")
       do j = 1, JA
       do i = 1, IA-1
       do k = 1, KA
          diagvar(k,i,j,I_VELX) = 2.D0 * var(k,i,j,2) / ( var(k,i+1,j,1)+var(k,i,j,1) )
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL18")
call START_COLLECTION("PRL19")
       do j = 1, JA
       do k = 1, KA
          diagvar(k,IA,j,I_VELX) = diagvar(k,IA-1,j,I_VELX)
       enddo
       enddo
call STOP_COLLECTION("PRL19")

call START_COLLECTION("PRL20")
       do j = 1, JA-1
       do i = 1, IA
       do k = 1, KA
          diagvar(k,i,j,I_VELY) = 2.D0 * var(k,i,j,3) / ( var(k,i,j+1,1)+var(k,i,j,1) )
       enddo 
       enddo
       enddo
call STOP_COLLECTION("PRL20")
call START_COLLECTION("PRL21")
       do i = 1, IA
       do k = 1, KA
          diagvar(k,i,JA,I_VELY) = diagvar(k,i,JA-1,I_VELY)
       enddo
       enddo
call STOP_COLLECTION("PRL21")

       !##### continuity equation #####
       ! at (x, y, interface)
call START_COLLECTION("PRL22")
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          mflx_hi(k,i,j,ZDIR) = 0.5D0 * diagvar(k,i,j,I_VELZ)                            &
                              * ( FACT_N * ( var(k+1,i,j,1)+var(k  ,i,j,1) )   &
                                + FACT_F * ( var(k+2,i,j,1)+var(k-1,i,j,1) ) ) &
                              + num_diff(k,i,j,1,ZDIR) * rdtrk
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL22")
call START_COLLECTION("PRL23")
       do j = JS, JE
       do i = IS, IE
          mflx_hi(WS  ,i,j,ZDIR) = 0.D0                                            ! bottom boundary
          mflx_hi(WS+1,i,j,ZDIR) = 0.5D0 * diagvar(WS+1,i,j,I_VELZ)              &
                                 * ( var(WS+2,i,j,1)+var(WS+1,i,j,1) ) & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,1,ZDIR) * rdtrk
          mflx_hi(WE-1,i,j,ZDIR) = 0.5D0 * diagvar(WE-1,i,j,I_VELZ)              &
                                 * ( var(WE  ,i,j,1)+var(WE-1,i,j,1) ) & ! just below the top boundary
                                 + num_diff(WE-1,i,j,1,ZDIR) * rdtrk
          mflx_hi(WE  ,i,j,ZDIR) = 0.D0                                            ! top boundary
       enddo
       enddo
call STOP_COLLECTION("PRL23")
       ! at (u, y, layer)
call START_COLLECTION("PRL24")
       do j = JS,   JE
       do i = IS-1, IE
       do k = KS,   KE
          mflx_hi(k,i,j,XDIR) = 0.5D0 * diagvar(k,i,j,I_VELX)                            &
                              * ( FACT_N * ( var(k,i+1,j,1)+var(k,i  ,j,1) )   &
                                + FACT_F * ( var(k,i+2,j,1)+var(k,i-1,j,1) ) ) &
                              + num_diff(k,i,j,1,XDIR) * rdtrk
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL24")
       ! at (x, v, layer)
call START_COLLECTION("PRL25")
       do j = JS-1, JE
       do i = IS,   IE
       do k = KS,   KE
          mflx_hi(k,i,j,YDIR) = 0.5D0 * diagvar(k,i,j,I_VELY)                            &
                              * ( FACT_N * ( var(k,i,j+1,1)+var(k,i,j  ,1) )   &
                                + FACT_F * ( var(k,i,j+2,1)+var(k,i,j-1,1) ) ) &
                              + num_diff(k,i,j,1,YDIR) * rdtrk
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL25")


       !##### momentum equation (z) #####
       ! at (x, y, layer)
call START_COLLECTION("PRL26")
       do j = JS,   JE
       do i = IS,   IE
       do k = KS,   KE
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELZ)+diagvar(k-1,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k  ,i,j,4)+var(k-1,i,j,4) )     &
                                + FACT_F * ( var(k+1,i,j,4)+var(k-2,i,j,4) ) )   &
                              + num_diff(k,i,j,4,ZDIR) * rdtrk
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL26")
call START_COLLECTION("PRL27")
       do j = JS,   JE
       do i = IS,   IE
          qflx_hi(KS-1,i,j,ZDIR) = 0.D0 ! bottom cell center
          qflx_hi(KE+1,i,j,ZDIR) = 0.D0 ! top    cell center
       enddo
       enddo
call STOP_COLLECTION("PRL27")
       ! at (u, y, interface)
call START_COLLECTION("PRL28")
       do j = JS,   JE
       do i = IS-1, IE
       do k = WS,   WE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i+1,j,4)+var(k,i  ,j,4) )     &
                                + FACT_F * ( var(k,i+2,j,4)+var(k,i-1,j,4) ) )   &
                              + num_diff(k,i,j,4,XDIR) * rdtrk
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL28")
       ! at (x, v, interface)
call START_COLLECTION("PRL29")
       do j = JS-1, JE
       do i = IS,   IE
       do k = WS,   WE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k+1,i,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j+1,4)+var(k,i,j  ,4) )     &
                                + FACT_F * ( var(k,i,j+2,4)+var(k,i,j-1,4) ) )   &
                              + num_diff(k,i,j,4,YDIR) * rdtrk
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL29")

       ! pressure
call START_COLLECTION("PRL30")
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          diagvar(k,i,j,I_PRES) = Pstd * ( var(k,i,j,5) * Rdry / Pstd )**CPovCV
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL30")


call START_COLLECTION("PRL31")
       !--- update momentum(z)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+1, WE-1
          var(k,i,j,4) = var_s(k,i,j,4) &
                            + dtrk * ( - ( ( qflx_hi(k+1,i,j,ZDIR)-qflx_hi(k,i  ,j  ,ZDIR) ) * RDZF(k)   &
                                         + ( qflx_hi(k  ,i,j,XDIR)-qflx_hi(k,i-1,j  ,XDIR) ) * RDXC(i)   &
                                         + ( qflx_hi(k  ,i,j,YDIR)-qflx_hi(k,i  ,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                       - ( diagvar(k+1,i,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDZF(k)     & ! pressure gradient force
                                       - ( var(k+1,i,j,1)+var(k,i,j,1) ) * 0.5D0 * GRAV        & ! gravity force
                                       + ray_damp(k,i,j,4)                                          ) ! additional damping force
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL31")

call START_COLLECTION("PRL32")
       !--- update density
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,1) = var_s(k,i,j,1) &
                            + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                         + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                         + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) ) ! divergence
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL32")

       !##### momentum equation (x) #####
       ! at (x, y, layer)
       do j = JS, JE
       do i = IS, IE+1
       do k = KS, KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELX)+diagvar(k,i-1,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i  ,j,2)+var(k,i-1,j,2) )     &
                                + FACT_F * ( var(k,i+1,j,2)+var(k,i-2,j,2) ) )   &
                              + num_diff(k,i,j,2,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do j = JS-1, JE
       do i = IS,   IE
       do k = KS,   KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELY)+diagvar(k,i,j,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j+1,2)+var(k,i,j  ,2) )     &
                                + FACT_F * ( var(k,i,j+2,2)+var(k,i,j-1,2) ) )   &
                              + num_diff(k,i,j,2,YDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i+1,j,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k+1,i,j,2)+var(k  ,i,j,2) )     &
                                + FACT_F * ( var(k+2,i,j,2)+var(k-1,i,j,2) ) )   &
                              + num_diff(k,i,j,2,ZDIR) * rdtrk
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_hi(WS  ,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(WS+1,i,j,ZDIR) = 0.25D0 * ( diagvar(WS+1,i+1,j,I_VELZ)+diagvar(WS+1,i,j,I_VELZ) ) &
                                 * ( var(WS+2,i,j,2)+var(WS+1,i,j,2) )                    & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,2,ZDIR) * rdtrk
          qflx_hi(WE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(WE-1,i+1,j,I_VELZ)+diagvar(WE-1,i,j,I_VELZ) ) &
                                 * ( var(WE  ,i,j,2)+var(WE-1,i,j,2) )                    & ! just below the top boundary
                                 + num_diff(WE-1,i,j,2,ZDIR) * rdtrk
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

call START_COLLECTION("PRL37")
       !--- update momentum(x)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,2) = var_s(k,i,j,2) &
                            + dtrk * ( - ( ( qflx_hi(k,i  ,j,ZDIR)-qflx_hi(k-1,i,j,  ZDIR) ) * RDZC(k)   &
                                         + ( qflx_hi(k,i+1,j,XDIR)-qflx_hi(k  ,i,j,  XDIR) ) * RDXF(i)   &
                                         + ( qflx_hi(k,i  ,j,YDIR)-qflx_hi(k  ,i,j-1,YDIR) ) * RDYC(j) ) & ! flux divergence
                                       - ( diagvar(k,i+1,j,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDXF(i)     & ! pressure gradient force
                                       + ray_damp(k,i,j,2)                                          ) ! additional damping force
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL37")

       !##### momentum equation (y) #####
       ! at (u, v, layer)
       do j = JS,   JE
       do i = IS-1, IE
       do k = KS,   KE
          qflx_hi(k,i,j,XDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELX)+diagvar(k,i,j,I_VELX) ) &
                              * ( FACT_N * ( var(k,i+1,j,3)+var(k,i  ,j,3) )     &
                                + FACT_F * ( var(k,i+2,j,3)+var(k,i-1,j,3) ) )   &
                              + num_diff(k,i,j,3,XDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, y, layer)
       do j = JS, JE+1
       do i = IS, IE
       do k = KS, KE
          qflx_hi(k,i,j,YDIR) = 0.25D0 * ( diagvar(k,i,j,I_VELY)+diagvar(k,i,j-1,I_VELY) ) &
                              * ( FACT_N * ( var(k,i,j  ,3)+var(k,i,j-1,3) )     &
                                + FACT_F * ( var(k,i,j+1,3)+var(k,i,j-2,3) ) )   &
                              + num_diff(k,i,j,3,YDIR) * rdtrk
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JS,   JE
       do i = IS,   IE
       do k = WS+2, WE-2
          qflx_hi(k,i,j,ZDIR) = 0.25D0 * ( diagvar(k,i,j+1,I_VELZ)+diagvar(k,i,j,I_VELZ) ) &
                              * ( FACT_N * ( var(k+1,i,j,3)+var(k  ,i,j,3) )     &
                                + FACT_F * ( var(k+2,i,j,3)+var(k-1,i,j,3) ) )   &
                              + num_diff(k,i,j,3,ZDIR) * rdtrk
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_hi(WS  ,i,j,ZDIR) = 0.D0                                                               ! bottom boundary
          qflx_hi(WS+1,i,j,ZDIR) = 0.25D0 * ( diagvar(WS+1,i,j+1,I_VELZ)+diagvar(WS+1,i,j,I_VELZ) ) &
                                 * ( var(WS+2,i,j,3)+var(WS+1,i,j,3) )                    & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,3,ZDIR) * rdtrk
          qflx_hi(WE-1,i,j,ZDIR) = 0.25D0 * ( diagvar(WE-1,i,j+1,I_VELZ)+diagvar(WE-1,i,j,I_VELZ) ) &
                                 * ( var(WE  ,i,j,3)+var(WE-1,i,j,3) )                    & ! just below the top boundary
                                 + num_diff(WE-1,i,j,3,ZDIR) * rdtrk
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0                                                               ! top boundary
       enddo
       enddo

call START_COLLECTION("PRL42")
       !--- update momentum(y)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,3) = var_s(k,i,j,3) &
                            + dtrk * ( - ( ( qflx_hi(k,i,j  ,ZDIR)-qflx_hi(k-1,i  ,j,ZDIR) ) * RDZC(k)   &
                                         + ( qflx_hi(k,i,j  ,XDIR)-qflx_hi(k  ,i-1,j,XDIR) ) * RDXC(i)   &
                                         + ( qflx_hi(k,i,j+1,YDIR)-qflx_hi(k  ,i  ,j,YDIR) ) * RDYF(j) ) & ! flux divergence
                                       - ( diagvar(k,i,j+1,I_PRES)-diagvar(k,i,j,I_PRES) ) * RDYF(j)     & ! pressure gradient force
                                       + ray_damp(k,i,j,3)                                          ) ! additional damping force
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL42")

       !##### Thermodynamic Equation #####

       do j = 1, JA !JS, JE
       do i = 1, IA !IS, IE
       do k = 1, KA !KS, KE
          diagvar(k,i,j,I_POTT) = var(k,i,j,5) / var(k,i,j,1) 
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
                              + num_diff(k,i,j,5,XDIR) * rdtrk
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
                              + num_diff(k,i,j,5,YDIR) * rdtrk
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
                              + num_diff(k,i,j,5,ZDIR) * rdtrk
       enddo
       enddo
       enddo
       do j = JS, JE
       do i = IS, IE
          qflx_hi(WS  ,i,j,ZDIR) = 0.D0                                                    ! bottom boundary
          qflx_hi(WS+1,i,j,ZDIR) = 0.5D0 * mflx_hi(WS+1,i,j,ZDIR)                        &
                                 * ( diagvar(WS+2,i,j,I_POTT)+diagvar(WS+1,i,j,I_POTT) ) & ! just above the bottom boundary
                                 + num_diff(WS+1,i,j,5,ZDIR) * rdtrk
          qflx_hi(WE-1,i,j,ZDIR) = 0.5D0 * mflx_hi(WE-1,i,j,ZDIR)                        &
                                 * ( diagvar(WE  ,i,j,I_POTT)+diagvar(WE-1,i,j,I_POTT) ) & ! just below the top boundary
                                 + num_diff(WE-1,i,j,5,ZDIR) * rdtrk
          qflx_hi(WE  ,i,j,ZDIR) = 0.D0                                                    ! top boundary
       enddo
       enddo

call START_COLLECTION("PRL48")
       !--- update rho*theta
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,5) = var_s(k,i,j,5) &
                            + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR)-qflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                         + ( qflx_hi(k,i,j,XDIR)-qflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                         + ( qflx_hi(k,i,j,YDIR)-qflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j) ) & ! divergence
                                       + ray_damp(k,i,j,5)                                          ) ! additional damping force
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL48")

       call COMM_vars( var(:,:,:,1), 1 )
       call COMM_vars( var(:,:,:,2), 2 )
       call COMM_vars( var(:,:,:,3), 3 )
       call COMM_vars( var(:,:,:,4), 4 )
       call COMM_vars( var(:,:,:,5), 5 )
       call COMM_wait( var(:,:,:,1), 1 )
       call COMM_wait( var(:,:,:,2), 2 )
       call COMM_wait( var(:,:,:,3), 3 )
       call COMM_wait( var(:,:,:,4), 4 )
       call COMM_wait( var(:,:,:,5), 5 )

    enddo ! RK loop

    call COMM_vars( mflx_hi(:,:,:,ZDIR), VA+ZDIR )
    call COMM_vars( mflx_hi(:,:,:,XDIR), VA+XDIR )
    call COMM_vars( mflx_hi(:,:,:,YDIR), VA+YDIR )
    call COMM_wait( mflx_hi(:,:,:,ZDIR), VA+ZDIR )
    call COMM_wait( mflx_hi(:,:,:,XDIR), VA+XDIR )
    call COMM_wait( mflx_hi(:,:,:,YDIR), VA+YDIR )

call STOP_COLLECTION("RK3")
call START_COLLECTION("FCT")
    !##### advection of scalar quantity #####

    if ( QA > 0 ) then

    ! calc low-order mass flux and high-low difference (at last step of RK)
    do j = JS,   JE
    do i = IS,   IE
    do k = WS+2, WE-2
       mflx_lo(k,i,j,ZDIR) = 0.5D0 * (     diagvar(k,i,j,I_VELZ)  * ( var(k+1,i,j,1)+var(k,i,j,1) ) &
                                     - abs(diagvar(k,i,j,I_VELZ)) * ( var(k+1,i,j,1)-var(k,i,j,1) ) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       mflx_lo(WS  ,i,j,ZDIR) = 0.D0                                                                         ! bottom boundary
       mflx_lo(WS+1,i,j,ZDIR) = 0.5D0 * diagvar(WS+1,i,j,I_VELZ) * ( var(WS+2,i,j,1)+var(WS+1,i,j,1) ) ! just above the bottom boundary
       mflx_lo(WE-1,i,j,ZDIR) = 0.5D0 * diagvar(WE-1,i,j,I_VELZ) * ( var(WE  ,i,j,1)+var(WE-1,i,j,1) ) ! just below the top boundary
       mflx_lo(WE  ,i,j,ZDIR) = 0.D0                                                                         ! top boundary 
    enddo
    enddo

    do j = JS,   JE
    do i = IS-1, IE
    do k = KS,   KE
       mflx_lo(k,i,j,XDIR) = 0.5D0 * (     diagvar(k,i,j,I_VELX)  * ( var(k,i+1,j,1)+var(k,i,j,1) ) &
                                     - abs(diagvar(k,i,j,I_VELX)) * ( var(k,i+1,j,1)-var(k,i,j,1) ) )
    enddo
    enddo
    enddo

    do j = JS-1, JE
    do i = IS,   IE
    do k = KS,   KE
       mflx_lo(k,i,j,YDIR) = 0.5D0 * (     diagvar(k,i,j,I_VELY)  * ( var(k,i,j+1,1)+var(k,i,j,1) ) &
                                     - abs(diagvar(k,i,j,I_VELY)) * ( var(k,i,j+1,1)-var(k,i,j,1) ) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ddiv(k,i,j) = ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RDZC(k) &
                   + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RDXC(i) &
                   + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RDYC(j)
    enddo
    enddo
    enddo

    endif

    if ( QA > 0 ) then

    do iq = 6, 5+QA

call START_COLLECTION("PRL50")
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
call STOP_COLLECTION("PRL50")

call START_COLLECTION("PRL51")
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
call STOP_COLLECTION("PRL51")

call START_COLLECTION("PRL52")
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
call STOP_COLLECTION("PRL52")

call START_COLLECTION("PRL53")
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
call STOP_COLLECTION("PRL53")

call START_COLLECTION("PRL54")
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          !--- update value with flux-divergence from the monotone scheme
          var(k,i,j,iq) = var_s(k,i,j,iq) &
                        + dtrk * ( - ( ( qflx_lo(k,i,j,XDIR)-qflx_lo(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                     + ( qflx_lo(k,i,j,YDIR)-qflx_lo(k  ,i,  j-1,YDIR) ) * RDYC(j)   &
                                     + ( qflx_lo(k,i,j,ZDIR)-qflx_lo(k-1,i,  j,  ZDIR) ) * RDZC(k) ) &
                                   + var_s(k,i,j,iq) * ddiv(k,i,j)                                   &
                                 ) / var(k,i,j,1)

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
             rjpls(k,i,j,XDIR) = (pjmax-var(k,i,j,iq)) / pjpls * abs((mflx_lo(k,i,j,XDIR)+mflx_lo(k  ,i-1,j  ,XDIR)) * 0.5D0)
             rjpls(k,i,j,YDIR) = (pjmax-var(k,i,j,iq)) / pjpls * abs((mflx_lo(k,i,j,YDIR)+mflx_lo(k  ,i  ,j-1,YDIR)) * 0.5D0)
             rjpls(k,i,j,ZDIR) = (pjmax-var(k,i,j,iq)) / pjpls * abs((mflx_lo(k,i,j,ZDIR)+mflx_lo(k-1,i  ,j  ,ZDIR)) * 0.5D0)
          else
             rjpls(k,i,j,XDIR:ZDIR) = 0.D0
          endif
          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             rjmns(k,i,j,XDIR) = (var(k,i,j,iq)-pjmin) / pjmns * abs((mflx_lo(k,i,j,XDIR)+mflx_lo(k  ,i-1,j  ,XDIR)) * 0.5D0)
             rjmns(k,i,j,YDIR) = (var(k,i,j,iq)-pjmin) / pjmns * abs((mflx_lo(k,i,j,YDIR)+mflx_lo(k  ,i  ,j-1,YDIR)) * 0.5D0)
             rjmns(k,i,j,ZDIR) = (var(k,i,j,iq)-pjmin) / pjmns * abs((mflx_lo(k,i,j,ZDIR)+mflx_lo(k-1,i  ,j  ,ZDIR)) * 0.5D0)
          else
             rjmns(k,i,j,XDIR:ZDIR) = 0.D0
          endif
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL54")

call START_COLLECTION("PRL55")
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
call STOP_COLLECTION("PRL55")

call START_COLLECTION("PRL56")
       !--- modify value with antidiffusive fluxes
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          var(k,i,j,iq) = var(k,i,j,iq) &
                        + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR)-qflx_anti(k-1,i,  j,  ZDIR) ) * RDZC(k)   &
                                     + ( qflx_anti(k,i,j,XDIR)-qflx_anti(k  ,i-1,j,  XDIR) ) * RDXC(i)   &
                                     + ( qflx_anti(k,i,j,YDIR)-qflx_anti(k  ,i,  j-1,YDIR) ) * RDYC(j) ) &
                                 ) / var(k,i,j,1)
       enddo
       enddo
       enddo
call STOP_COLLECTION("PRL56")

       call COMM_vars( var(:,:,:,iq), iq )
       call COMM_wait( var(:,:,:,iq), iq )

    enddo ! scalar quantities loop

    endif

call STOP_COLLECTION("FCT")

    enddo ! dynamical steps

    call fpcoll_stop

    ! check total mass
    call COMM_total( var(:,:,:,1:6), A_NAME(1:6) )

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
