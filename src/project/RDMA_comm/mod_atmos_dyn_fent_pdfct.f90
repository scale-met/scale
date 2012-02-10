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
    use mod_time, only: &
       TIME_NSTEP_ATMOS_DYN
    use mod_comm, only: &
       COMM_vars, &
       COMM_wait
    use mod_atmos_vars, only: &
       var => atmos_var, &
       QA  => A_QA, &
       I_DENS,      &
       I_MOMX,      &
       I_MOMY,      &
       I_MOMZ,      &
       I_RHOT
    implicit none

    integer :: iq, rko, step
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("DYNAMICS")
#endif

    do step = 1, TIME_NSTEP_ATMOS_DYN
    !##### Start RK #####
    do rko = 1, RK
       if ( rko > 1 ) then
          call COMM_wait( var(:,:,:,I_DENS), I_DENS )
          call COMM_wait( var(:,:,:,I_MOMZ), I_MOMZ )
          call COMM_wait( var(:,:,:,I_MOMX), I_MOMX )
          call COMM_wait( var(:,:,:,I_MOMY), I_MOMY )
          call COMM_wait( var(:,:,:,I_RHOT), I_RHOT )
       endif

       call COMM_vars( var(:,:,:,I_DENS), I_DENS )
       call COMM_vars( var(:,:,:,I_MOMZ), I_MOMZ )
       call COMM_vars( var(:,:,:,I_MOMX), I_MOMX )
       call COMM_vars( var(:,:,:,I_MOMY), I_MOMY )
       call COMM_vars( var(:,:,:,I_RHOT), I_RHOT )
    enddo ! RK loop

    call COMM_wait( var(:,:,:,I_DENS), I_DENS )
    call COMM_wait( var(:,:,:,I_MOMZ), I_MOMZ )
    call COMM_wait( var(:,:,:,I_MOMX), I_MOMX )
    call COMM_wait( var(:,:,:,I_MOMY), I_MOMY )

    if ( QA > 0 ) then
       do iq = 6, 5+QA
          call COMM_wait( var(:,:,:,iq-1), iq-1 )
          call COMM_vars( var(:,:,:,iq), iq )
       enddo ! scalar quantities loop
       call COMM_wait( var(:,:,:,iq-1), iq-1 )
    else
       call COMM_wait( var(:,:,:,I_RHOT), I_RHOT )
    endif

    enddo ! dynamical steps

#ifdef _FPCOLL_
call STOP_COLLECTION("DYNAMICS")
#endif

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
