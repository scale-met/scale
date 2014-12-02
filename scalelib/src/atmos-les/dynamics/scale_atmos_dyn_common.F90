!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          common subroutines for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] splited from dynamical core
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_filter_setup
  public :: ATMOS_DYN_numfilter_coef
  public :: ATMOS_DYN_numfilter_coef_q
  public :: ATMOS_DYN_filter_tend
  public :: ATMOS_DYN_fct

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  ! advection settings
  real(RP), public, parameter :: FACT_N =   7.0_RP / 12.0_RP !  7/12: fourth, 1: second
  real(RP), public, parameter :: FACT_F = - 1.0_RP / 12.0_RP ! -1/12: fourth, 0: second

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  real(RP), allocatable :: CNZ3(:,:,:)
  real(RP), allocatable :: CNX3(:,:,:)
  real(RP), allocatable :: CNY3(:,:,:)
  real(RP), allocatable :: CNZ4(:,:,:)
  real(RP), allocatable :: CNX4(:,:,:)
  real(RP), allocatable :: CNY4(:,:,:)

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_filter_setup( &
       CDZ, CDX, CDY, FDZ, FDX, FDY        )
    implicit none

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)
    real(RP), intent(in)  :: FDZ(KA-1)
    real(RP), intent(in)  :: FDX(IA-1)
    real(RP), intent(in)  :: FDY(JA-1)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! allocation
    allocate( CNZ3(3,KA,2) )
    allocate( CNX3(3,IA,2) )
    allocate( CNY3(3,JA,2) )
    allocate( CNZ4(5,KA,2) )
    allocate( CNX4(5,IA,2) )
    allocate( CNY4(5,JA,2) )

#ifdef DEBUG
    CNZ3(:,:,:) = UNDEF
    CNX3(:,:,:) = UNDEF
    CNY3(:,:,:) = UNDEF
    CNZ4(:,:,:) = UNDEF
    CNX4(:,:,:) = UNDEF
    CNY4(:,:,:) = UNDEF
#endif

    ! z direction
    do k = KS+1, KE-1
       CNZ3(1,k,1) = 1.0_RP / ( FDZ(k  ) * CDZ(k  ) * FDZ(k-1) )
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



  end subroutine ATMOS_DYN_filter_setup

  !-----------------------------------------------------------------------------
  !> Calc coefficient of numerical filter
  subroutine ATMOS_DYN_numfilter_coef( &
       num_diff,                               &
       DENS, MOMZ, MOMX, MOMY, RHOT,           &
       CDZ, CDX, CDY, FDZ, FDX, FDY, DT,       &
       REF_dens, REF_pott,                     &
       ND_COEF, ND_ORDER, ND_SFC_FACT, ND_USE_RS )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: num_diff(KA,IA,JA,5,3)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)
    real(RP), intent(in)  :: FDZ(KA-1)
    real(RP), intent(in)  :: FDX(IA-1)
    real(RP), intent(in)  :: FDY(JA-1)

    real(RP), intent(in)  :: DT

    real(RP), intent(in)  :: REF_dens(KA,IA,JA)
    real(RP), intent(in)  :: REF_pott(KA,IA,JA)

    real(RP), intent(in)  :: ND_COEF
    integer,  intent(in)  :: ND_ORDER
    real(RP), intent(in)  :: ND_SFC_FACT
    logical,  intent(in)  :: ND_USE_RS

    ! diagnostic variables
    real(RP) :: VELZ (KA,IA,JA) ! velocity w [m/s]
    real(RP) :: VELX (KA,IA,JA) ! velocity u [m/s]
    real(RP) :: VELY (KA,IA,JA) ! velocity v [m/s]
    real(RP) :: POTT (KA,IA,JA) ! potential temperature [K]

    real(RP) :: dens_diff(KA,IA,JA)     ! anomary of density
    real(RP) :: pott_diff(KA,IA,JA)     ! anomary of rho * pott

    real(RP) :: work(KA,IA,JA,3,2)
    integer  :: iwork

    real(RP) :: DIFF4
    integer  :: nd_order4
    real(RP) :: nd_coef_cdz(KA)
    real(RP) :: nd_coef_cdx(IA)
    real(RP) :: nd_coef_cdy(JA)
    real(RP) :: nd_coef_fdz(KA-1)
    real(RP) :: nd_coef_fdx(IA-1)
    real(RP) :: nd_coef_fdy(JA-1)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! numerical diffusion
    nd_order4 = nd_order * 4
    DIFF4 = ND_COEF / ( 2**(nd_order4) * DT )
    do k = KS-1, KE
       nd_coef_cdz(k) = DIFF4 * CDZ(k)**nd_order4
    end do
    do k = KS+1, KE-1
       nd_coef_fdz(k) = DIFF4 * FDZ(k)**nd_order4
    end do
    do i = IS, IE
       nd_coef_cdx(i) = DIFF4 * CDX(i)**nd_order4
       nd_coef_fdx(i) = DIFF4 * FDX(i)**nd_order4
    end do
    do j = JS, JE
       nd_coef_cdy(j) = DIFF4 * CDY(j)**nd_order4
       nd_coef_fdy(j) = DIFF4 * FDY(j)**nd_order4
    end do


    !###########################################################################
    ! 1st order coefficients
    !###########################################################################

    if ( .not. ND_USE_RS ) then

       call PROF_rapstart("NumFilter Main", 3)

       do j = JS-2, JE+2
       do i = IS-2, IE+2
       do k = KS, KE
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       end do
       end do
       end do

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS+1, KE-1
          dens_diff(k,i,j) = ( ( DENS(k,i,j)                                             ) * 3.0_RP &
                             + ( DENS(k,i+1,j)+DENS(k,i-1,j)+DENS(k,i,j+1)+DENS(k,i,j-1) ) * 2.0_RP &
                             + ( DENS(k,i+2,j)+DENS(k,i-2,j)+DENS(k,i,j+2)+DENS(k,i,j-2) ) &
                             + ( DENS(k+1,i,j)+DENS(k-1,i,j)                             ) * 2.0_RP &
                             ) / 19.0_RP

          pott_diff(k,i,j) = ( ( POTT(k,i,j)                                             ) * 3.0_RP &
                             + ( POTT(k,i+1,j)+POTT(k,i-1,j)+POTT(k,i,j+1)+POTT(k,i,j-1) ) * 2.0_RP &
                             + ( POTT(k,i+2,j)+POTT(k,i-2,j)+POTT(k,i,j+2)+POTT(k,i,j-2) ) &
                             + ( POTT(k+1,i,j)+POTT(k-1,i,j)                             ) * 2.0_RP &
                             ) / 19.0_RP
       enddo
       enddo
       enddo

       do j  = JS, JE
       do i  = IS, IE
          dens_diff(KS,i,j) = ( ( DENS(KS,i,j)                                                ) * 3.0_RP &
                              + ( DENS(KS,i+1,j)+DENS(KS,i-1,j)+DENS(KS,i,j+1)+DENS(KS,i,j-1) ) * 2.0_RP &
                              + ( DENS(KS,i+2,j)+DENS(KS,i-2,j)+DENS(KS,i,j+2)+DENS(KS,i,j-2) ) &
                              + ( DENS(KS+1,i,j)                                              ) * 2.0_RP &
                              ) / 17.0_RP
          dens_diff(KE,i,j) = ( ( DENS(KE,i,j)                                                ) * 3.0_RP &
                              + ( DENS(KE,i+1,j)+DENS(KE,i-1,j)+DENS(KE,i,j+1)+DENS(KE,i,j-1) ) * 2.0_RP &
                              + ( DENS(KE,i+2,j)+DENS(KE,i-2,j)+DENS(KE,i,j+2)+DENS(KE,i,j-2) ) &
                              + ( DENS(KE-1,i,j)                                              ) * 2.0_RP &
                              ) / 17.0_RP

          pott_diff(KS,i,j) = ( ( POTT(KS,i,j)                                                ) * 3.0_RP &
                              + ( POTT(KS,i+1,j)+POTT(KS,i-1,j)+POTT(KS,i,j+1)+POTT(KS,i,j-1) ) * 2.0_RP &
                              + ( POTT(KS,i+2,j)+POTT(KS,i-2,j)+POTT(KS,i,j+2)+POTT(KS,i,j-2) ) &
                              + ( POTT(KS+1,i,j)                                              ) * 2.0_RP &
                              ) / 17.0_RP
          pott_diff(KE,i,j) = ( ( POTT(KE,i,j)                                                ) * 3.0_RP &
                              + ( POTT(KE,i+1,j)+POTT(KE,i-1,j)+POTT(KE,i,j+1)+POTT(KE,i,j-1) ) * 2.0_RP &
                              + ( POTT(KE,i+2,j)+POTT(KE,i-2,j)+POTT(KE,i,j+2)+POTT(KE,i,j-2) ) &
                              + ( POTT(KE-1,i,j)                                              ) * 2.0_RP &
                              ) / 17.0_RP
       end do
       end do

       call PROF_rapend  ("NumFilter Main")

       call PROF_rapstart("NumFilter Comm", 3)

       call COMM_vars8( dens_diff,  1 )
       call COMM_vars8( pott_diff,  2 )

       call COMM_wait ( dens_diff,  1 )
       call COMM_wait ( pott_diff,  2 )

       call PROF_rapend  ("NumFilter Comm")

    end if

    call PROF_rapstart("NumFilter Main", 3)



    !-----< density >-----

    if ( ND_USE_RS ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          dens_diff(k,i,j) = DENS(k,i,j) - REF_dens(k,i,j)
       enddo
       enddo
       enddo
    endif

    call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                       dens_diff, & ! (in)
                       nd_order, & ! (in)
                       0, 0, 0, KE )

    call PROF_rapstart("NumFilter Main", 3)

    !-----< density >-----

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       num_diff(k,i,j,I_DENS,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-2,i,j,I_DENS,ZDIR) = 0.0_RP
       num_diff(KE+1:KA  ,i,j,I_DENS,ZDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_DENS,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_DENS,XDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_DENS,XDIR) = num_diff(KS  ,i,j,I_DENS,XDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_DENS,XDIR) = num_diff(KS+1,i,j,I_DENS,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_DENS,XDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_DENS,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_DENS,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_DENS,YDIR) = num_diff(KS  ,i,j,I_DENS,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_DENS,YDIR) = num_diff(KS+1,i,j,I_DENS,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_DENS,YDIR) = 0.0_RP
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call PROF_rapstart("NumFilter Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_DENS,ZDIR),  1 )
    call COMM_vars8( num_diff(:,:,:,I_DENS,XDIR),  2 )
    call COMM_vars8( num_diff(:,:,:,I_DENS,YDIR),  3 )

    call PROF_rapend  ("NumFilter Comm")


    !-----< z-momentum >-----

    call PROF_rapstart("NumFilter Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS, KE-1
       VELZ(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                       VELZ, & ! (in)
                       nd_order, & ! (in)
                       1, 0, 0, KE-1 )

    call PROF_rapstart("NumFilter Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE-1
       num_diff(k,i,j,I_MOMZ,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_fdz(k) &
                                   * DENS(k,i,j)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS,i,j,I_MOMZ,ZDIR) = 0.0_RP
       num_diff(KE:KA,i,j,I_MOMZ,ZDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMZ,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                   * 0.25_RP * ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMZ,XDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMZ,XDIR) = num_diff(KS  ,i,j,I_MOMZ,XDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMZ,XDIR) = num_diff(KS+1,i,j,I_MOMZ,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE:KA  ,i,j,I_MOMZ,XDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMZ,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                                   * 0.25_RP * ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMZ,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMZ,YDIR) = num_diff(KS  ,i,j,I_MOMZ,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMZ,YDIR) = num_diff(KS+1,i,j,I_MOMZ,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

       num_diff(KE:KA  ,i,j,I_MOMZ,YDIR) = 0.0_RP
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call PROF_rapstart("NumFilter Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_MOMZ,ZDIR),  4 )
    call COMM_vars8( num_diff(:,:,:,I_MOMZ,XDIR),  5 )
    call COMM_vars8( num_diff(:,:,:,I_MOMZ,YDIR),  6 )

    call PROF_rapend  ("NumFilter Comm")


    !-----< x-momentum >-----

    call PROF_rapstart("NumFilter Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = KS, KE
       VELX(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                       VELX, & ! (in)
                       nd_order, & ! (in)
                       0, 1, 0, KE )


    call PROF_rapstart("NumFilter Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMX,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                   * 0.25_RP * ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMX,ZDIR) = 0.0_RP
       num_diff(KE:KA  ,i,j,I_MOMX,ZDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_MOMX,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_fdx(i) &
                                   * DENS(k,i,j)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_MOMX,XDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMX,XDIR) = num_diff(KS  ,i,j,I_MOMX,XDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMX,XDIR) = num_diff(KS+1,i,j,I_MOMX,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

       num_diff(KE+1:KA  ,i,j,I_MOMX,XDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_MOMX,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                                   * 0.25_RP * ( DENS(k,i+1,j+1)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_MOMX,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMX,YDIR) = num_diff(KS  ,i,j,I_MOMX,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMX,YDIR) = num_diff(KS+1,i,j,I_MOMX,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_MOMX,YDIR) = 0.0_RP
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call PROF_rapstart("NumFilter Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_MOMX,ZDIR),  7 )
    call COMM_vars8( num_diff(:,:,:,I_MOMX,XDIR),  8 )
    call COMM_vars8( num_diff(:,:,:,I_MOMX,YDIR),  9 )

    call PROF_rapend  ("NumFilter Comm")


    !-----< y-momentum >-----

    call PROF_rapstart("NumFilter Comm", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = KS, KE
       VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo

    call PROF_rapend  ("NumFilter Comm")

    call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                       VELY, & ! (in)
                       nd_order, & ! (in)
                       0, 0, 1, KE )

    call PROF_rapstart("NumFilter Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMY,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                   * 0.25_RP * ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMY,ZDIR) = 0.0_RP
       num_diff(KE:KA  ,i,j,I_MOMY,ZDIR) = 0.0_RP
    end do
    end do

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_MOMY,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                   * 0.25_RP * ( DENS(k,i+1,j+1)+DENS(k,i,j+1)+DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_MOMY,XDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMY,XDIR) = num_diff(KS  ,i,j,I_MOMY,XDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMY,XDIR) = num_diff(KS+1,i,j,I_MOMY,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_MOMY,XDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_MOMY,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_fdy(j) &
                                   * DENS(k,i,j)
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_MOMY,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMY,YDIR) = num_diff(KS  ,i,j,I_MOMY,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMY,YDIR) = num_diff(KS+1,i,j,I_MOMY,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_MOMY,YDIR) = 0.0_RP
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call PROF_rapstart("NumFilter Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_MOMY,ZDIR), 10 )
    call COMM_vars8( num_diff(:,:,:,I_MOMY,XDIR), 11 )
    call COMM_vars8( num_diff(:,:,:,I_MOMY,YDIR), 12 )

    call PROF_rapend  ("NumFilter Comm")

    !-----< rho * theta >-----

    call PROF_rapstart("NumFilter Main", 3)

    if ( ND_USE_RS ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE
          pott_diff(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) - REF_pott(k,i,j)
       enddo
       enddo
       enddo
    endif

    call PROF_rapend  ("NumFilter Main")

    call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                       pott_diff, & ! (in)
                       nd_order, & ! (in)
                       0, 0, 0, KE )

    call PROF_rapstart("NumFilter Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_RHOT,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                   * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-2,i,j,I_RHOT,ZDIR) = 0.0_RP
       num_diff(KS-1,i,j,I_RHOT,ZDIR) = work(KS-1,i,j,ZDIR,iwork) * nd_coef_cdz(KS-1) &
                                      * DENS(KS,i,j)
       num_diff(KE  ,i,j,I_RHOT,ZDIR) = work(KE  ,i,j,ZDIR,iwork) * nd_coef_cdz(KE  ) &
                                      * DENS(KE,i,j)
       num_diff(KE+1:KA  ,i,j,I_RHOT,ZDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_RHOT,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                   * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_RHOT,XDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_RHOT,XDIR) = num_diff(KS  ,i,j,I_RHOT,XDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_RHOT,XDIR) = num_diff(KS+1,i,j,I_RHOT,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_RHOT,XDIR) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_RHOT,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                                   * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_RHOT,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_RHOT,YDIR) = num_diff(KS  ,i,j,I_RHOT,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_RHOT,YDIR) = num_diff(KS+1,i,j,I_RHOT,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_RHOT,YDIR) = 0.0_RP
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call PROF_rapstart("NumFilter Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_RHOT,ZDIR), 13 )
    call COMM_vars8( num_diff(:,:,:,I_RHOT,XDIR), 14 )
    call COMM_vars8( num_diff(:,:,:,I_RHOT,YDIR), 15 )

    call COMM_wait ( num_diff(:,:,:,I_DENS,ZDIR),  1 )
    call COMM_wait ( num_diff(:,:,:,I_DENS,XDIR),  2 )
    call COMM_wait ( num_diff(:,:,:,I_DENS,YDIR),  3 )
    call COMM_wait ( num_diff(:,:,:,I_MOMZ,ZDIR),  4 )
    call COMM_wait ( num_diff(:,:,:,I_MOMZ,XDIR),  5 )
    call COMM_wait ( num_diff(:,:,:,I_MOMZ,YDIR),  6 )
    call COMM_wait ( num_diff(:,:,:,I_MOMX,ZDIR),  7 )
    call COMM_wait ( num_diff(:,:,:,I_MOMX,XDIR),  8 )
    call COMM_wait ( num_diff(:,:,:,I_MOMX,YDIR),  9 )
    call COMM_wait ( num_diff(:,:,:,I_MOMY,ZDIR), 10 )
    call COMM_wait ( num_diff(:,:,:,I_MOMY,XDIR), 11 )
    call COMM_wait ( num_diff(:,:,:,I_MOMY,YDIR), 12 )
    call COMM_wait ( num_diff(:,:,:,I_RHOT,ZDIR), 13 )
    call COMM_wait ( num_diff(:,:,:,I_RHOT,XDIR), 14 )
    call COMM_wait ( num_diff(:,:,:,I_RHOT,YDIR), 15 )

    call PROF_rapend  ("NumFilter Comm")

    return
  end subroutine ATMOS_DYN_numfilter_coef

  !-----------------------------------------------------------------------------
  !> Calc coefficient of numerical filter
  subroutine ATMOS_DYN_numfilter_coef_q( &
       num_diff_q,                             &
       DENS, QTRC,                             &
       CDZ, CDX, CDY, dt,                      &
       REF_qv, iq,                             &
       ND_COEF, ND_ORDER, ND_SFC_FACT, ND_USE_RS )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: num_diff_q(KA,IA,JA,3)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA)

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)

    real(RP), intent(in)  :: dt

    real(RP), intent(in)  :: REF_qv(KA,IA,JA)
    integer,  intent(in)  :: iq

    real(RP), intent(in)  :: ND_COEF
    integer,  intent(in)  :: ND_ORDER
    real(RP), intent(in)  :: ND_SFC_FACT
    logical,  intent(in)  :: ND_USE_RS

    real(RP) :: qv_diff(KA,IA,JA) ! anomary of water vapor

    real(RP) :: work(KA,IA,JA,3,2)
    integer  :: iwork

    real(RP) :: DIFF4
    integer  :: nd_order4
    real(RP) :: nd_coef_cdz(KA)
    real(RP) :: nd_coef_cdx(IA)
    real(RP) :: nd_coef_cdy(JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !###########################################################################
    ! 1st order coefficients
    !###########################################################################

    nd_order4 = nd_order * 4
    DIFF4 = ND_COEF / ( 2**(nd_order4) * DT )
    do k = KS-1, KE
       nd_coef_cdz(k) = DIFF4 * CDZ(k)**nd_order4
    end do
    do i = IS, IE
       nd_coef_cdx(i) = DIFF4 * CDX(i)**nd_order4
    end do
    do j = JS, JE
       nd_coef_cdy(j) = DIFF4 * CDY(j)**nd_order4
    end do

    if ( iq == I_QV .and. (.not. ND_USE_RS) ) then

       call PROF_rapstart("NumFilter Main", 3)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS+1, KE-1
          qv_diff(k,i,j) = ( ( QTRC(k,i,j)                                             ) * 3.0_RP &
                           + ( QTRC(k,i+1,j)+QTRC(k,i-1,j)+QTRC(k,i,j+1)+QTRC(k,i,j-1) ) * 2.0_RP &
                           + ( QTRC(k,i+2,j)+QTRC(k,i-2,j)+QTRC(k,i,j+2)+QTRC(k,i,j-2) ) &
                           + ( QTRC(k+1,i,j)+QTRC(k-1,i,j)                             ) * 2.0_RP &
                           ) / 19.0_RP
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS-1, JE+2
       do i  = IS-1, IE+2
          qv_diff(KS,i,j) = ( ( QTRC(KS,i,j)                                                ) * 3.0_RP &
                            + ( QTRC(KS,i+1,j)+QTRC(KS,i-1,j)+QTRC(KS,i,j+1)+QTRC(KS,i,j-1) ) * 2.0_RP &
                            + ( QTRC(KS,i+2,j)+QTRC(KS,i-2,j)+QTRC(KS,i,j+2)+QTRC(KS,i,j-2) ) &
                            + ( QTRC(KS+1,i,j)                                              ) * 2.0_RP &
                            ) / 17.0_RP
          qv_diff(KE,i,j) = ( ( QTRC(KE,i,j)                                                ) * 3.0_RP &
                            + ( QTRC(KE,i+1,j)+QTRC(KE,i-1,j)+QTRC(KE,i,j+1)+QTRC(KE,i,j-1) ) * 2.0_RP &
                            + ( QTRC(KE,i+2,j)+QTRC(KE,i-2,j)+QTRC(KE,i,j+2)+QTRC(KE,i,j-2) ) &
                            + ( QTRC(KE-1,i,j)                                              ) * 2.0_RP &
                            ) / 17.0_RP
       end do
       end do

       call PROF_rapend  ("NumFilter Main")

       call PROF_rapstart("NumFilter Comm", 3)

       call COMM_vars8(qv_diff, 1)
       call COMM_wait (qv_diff, 1)

       call PROF_rapend  ("NumFilter Comm")

    end if

    if ( iq == I_QV ) then

       if ( ND_USE_RS ) then

          call PROF_rapstart("NumFilter Main", 3)

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JS-1, JE+2
          do i = IS-1, IE+2
          do k = KS, KE
             qv_diff(k,i,j) = QTRC(k,i,j) - REF_qv(k,i,j)
          enddo
          enddo
          enddo

          call PROF_rapend  ("NumFilter Main")

       endif

       call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                          qv_diff, & ! (in)
                          nd_order, & ! (in)
                          0, 0, 0, KE )

    else ! iq /= I_QV

       call calc_numdiff( work(:,:,:,:,1), iwork, & ! (out)
                          QTRC, & ! (in)
                          nd_order, & ! (in)
                          0, 0, 0, KE )

    endif ! QV or not?


    call PROF_rapstart("NumFilter Main", 3)


    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff_q(k,i,j,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                              * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff_q(KS-1,i,j,ZDIR) = work(KS-1,i,j,ZDIR,iwork) * nd_coef_cdz(KS-1) &
                                 * DENS(KS,i,j)
       num_diff_q(KE  ,i,j,ZDIR) = work(KE  ,i,j,ZDIR,iwork) * nd_coef_cdz(KE  ) &
                                 * DENS(KE,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff_q(k,i,j,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                              * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff_q(KS  ,i,j,XDIR) = num_diff_q(KS  ,i,j,XDIR) * ND_SFC_FACT
       num_diff_q(KS+1,i,j,XDIR) = num_diff_q(KS+1,i,j,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff_q(k,i,j,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                              * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       num_diff_q(KS  ,i,j,YDIR) = num_diff_q(KS  ,i,j,YDIR) * ND_SFC_FACT
       num_diff_q(KS+1,i,j,YDIR) = num_diff_q(KS+1,i,j,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
    enddo
    enddo

    call PROF_rapend  ("NumFilter Main")

    call PROF_rapstart("NumFilter Comm", 3)

    call COMM_vars8( num_diff_q(:,:,:,ZDIR), 1 )
    call COMM_vars8( num_diff_q(:,:,:,XDIR), 2 )
    call COMM_vars8( num_diff_q(:,:,:,YDIR), 3 )

    call COMM_wait ( num_diff_q(:,:,:,ZDIR), 1 )
    call COMM_wait ( num_diff_q(:,:,:,XDIR), 2 )
    call COMM_wait ( num_diff_q(:,:,:,YDIR), 3 )

    call PROF_rapend  ("NumFilter Comm")

    return
  end subroutine ATMOS_DYN_numfilter_coef_q

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_filter_tend( &
       phi_t, &
       phi, &
       rdz, rdx, rdy, &
       KO, IO, JO )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    real(RP), intent(out) :: phi_t(KA,IA,JA)
    real(RP), intent(in ) :: phi  (KA,IA,JA)
    real(RP), intent(in ) :: rdz(:)
    real(RP), intent(in ) :: rdx(:)
    real(RP), intent(in ) :: rdy(:)
    integer , intent(in ) :: KO
    integer , intent(in ) :: IO
    integer , intent(in ) :: JO

    real(RP) :: flux(KA,IA,JA,3)

    integer :: k, i, j

    call calc_diff3( flux, & ! (out)
                     phi, & ! (in)
                     KO, IO, JO )

    call COMM_vars8( flux(:,:,:,XDIR), 1 )
    call COMM_vars8( flux(:,:,:,YDIR), 2 )
    call COMM_wait ( flux(:,:,:,XDIR), 1 )
    call COMM_wait ( flux(:,:,:,YDIR), 2 )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       phi_t(k,i,j) = ( flux(k+KO,i,j,ZDIR) - flux(k-1+KO,i,j,ZDIR) ) * RDZ(k) &
                    + ( flux(k,i+IO,j,XDIR) - flux(k,i-1+IO,j,XDIR) ) * RDX(i) &
                    + ( flux(k,i,j+JO,YDIR) - flux(k,i,j-1+JO,YDIR) ) * RDY(j)
    end do
    end do
    end do

    return
  end subroutine ATMOS_DYN_filter_tend

  !-----------------------------------------------------------------------------
  subroutine calc_numdiff(&
       work, iwork, &
       data, &
       nd_order, &
       KO, IO, JO, KEE )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    real(RP), intent(out) :: work(KA,IA,JA,3,2)
    integer,  intent(out) :: iwork
    real(RP), intent(in)  :: data(KA,IA,JA)
    integer,  intent(in)  :: nd_order
    integer,  intent(in)  :: KO
    integer,  intent(in)  :: IO
    integer,  intent(in)  :: JO
    integer,  intent(in)  :: KEE

    integer :: i_in, i_out, i_tmp

    integer :: no

    call PROF_rapstart("NumFilter Main", 3)

    call calc_diff3( work(:,:,:,:,1), & ! (out)
                     data, & ! (in)
                     KO, IO, JO ) ! (in)

    call PROF_rapend  ("NumFilter Main")

    !###########################################################################
    ! High order coefficients
    !###########################################################################

    i_in  = 1
    i_out = 2

    do no = 2, nd_order

       call PROF_rapstart("NumFilter Comm", 3)

       call COMM_vars8( work(:,:,:,ZDIR,i_in),  1 )
       call COMM_vars8( work(:,:,:,XDIR,i_in),  2 )
       call COMM_vars8( work(:,:,:,YDIR,i_in),  3 )

       call COMM_wait ( work(:,:,:,ZDIR,i_in),  1 )
       call COMM_wait ( work(:,:,:,XDIR,i_in),  2 )
       call COMM_wait ( work(:,:,:,YDIR,i_in),  3 )

       call PROF_rapend  ("NumFilter Comm")

       call PROF_rapstart("NumFilter Main", 3)

       call calc_diff4( work(:,:,:,:,i_out), & ! (out)
                        work(:,:,:,:,i_in), & ! (in)
                        CNZ4(:,:,1),             & ! (in)
                        CNX4(:,:,1),             & ! (in)
                        CNY4(:,:,1),             & ! (in)
                        KE                       ) ! (in)

       call PROF_rapend  ("NumFilter Main")

       ! swap pointer target
       i_tmp = i_in
       i_in  = i_out
       i_out = i_tmp
    enddo

    iwork = i_in

    return
  end subroutine calc_numdiff

  !-----------------------------------------------------------------------------
  subroutine calc_diff3( &
       diff, &
       phi, &
       KO, IO, JO )
    implicit none
    real(RP), intent(out) :: diff(KA,IA,JA,3)
    real(RP), intent(in ) :: phi(KA,IA,JA)
    integer , intent(in ) :: KO
    integer , intent(in ) :: IO
    integer , intent(in ) :: JO

    integer :: k, i, j

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE-2-KO
#ifdef DEBUG
       call CHECK( __LINE__, phi(k+2,i,j) )
       call CHECK( __LINE__, phi(k+1,i,j) )
       call CHECK( __LINE__, phi(k  ,i,j) )
       call CHECK( __LINE__, phi(k-1,i,j) )
#endif
       diff(k+KO,i,j,ZDIR) = ( + CNZ3(1,k+1,1+KO) * phi(k+2,i,j) &
                               - CNZ3(2,k+1,1+KO) * phi(k+1,i,j) &
                               + CNZ3(3,k+1,1+KO) * phi(k  ,i,j) &
                               - CNZ3(1,k  ,1+KO) * phi(k-1,i,j) )
    enddo
    enddo
    enddo

    if ( KO==1 ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
#ifdef DEBUG
          call CHECK( __LINE__, phi(KS+1,i,j) )
          call CHECK( __LINE__, phi(KS  ,i,j) )
#endif
          diff(KS,i,j,ZDIR) = ( + CNZ3(1,KS  ,2) * phi(KS+1,i,j) &
                                - CNZ3(2,KS  ,2) * phi(KS  ,i,j) &
                                - CNZ3(1,KS-1,2) * phi(KS+1,i,j) )
       end do
       end do
    end if
    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
#ifdef DEBUG
       call CHECK( __LINE__, phi(KS+2,i,j) )
       call CHECK( __LINE__, phi(KS+1,i,j) )
       call CHECK( __LINE__, phi(KS,i,j) )
#endif
       diff(KS+KO,i,j,ZDIR) = ( + CNZ3(1,KS+1,1+KO) * phi(KS+2,i,j) &
                                - CNZ3(2,KS+1,1+KO) * phi(KS+1,i,j) &
                                + CNZ3(3,KS+1,1+KO) * phi(KS  ,i,j) &
                                - CNZ3(1,KS  ,1)    * phi(KS+1,i,j) * (1-KO) )
       diff(KS-1,i,j,ZDIR) = - diff(KS  ,i,j,ZDIR)
       diff(KS-2,i,j,ZDIR) = - diff(KS+1,i,j,ZDIR)
    enddo
    enddo

    if ( KO==0 ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
#ifdef DEBUG
       call CHECK( __LINE__, phi(KE,i,j) )
       call CHECK( __LINE__, phi(KE-1,i,j) )
       call CHECK( __LINE__, phi(KE-2,i,j) )
#endif
          diff(KE-1,i,j,ZDIR) = ( + CNZ3(1,KE  ,1+KO) * phi(KE-1,i,j) &
                                  - CNZ3(2,KE  ,1+KO) * phi(KE  ,i,j) &
                                  + CNZ3(3,KE  ,1+KO) * phi(KE-1,i,j) &
                                  - CNZ3(1,KE-1,1+KO) * phi(KE-2,i,j) )
          diff(KE+2,i,j,ZDIR) = 0.0_RP
       end do
       end do
    else
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
#ifdef DEBUG
          call CHECK( __LINE__, phi(KE-1,i,j) )
          call CHECK( __LINE__, phi(KE-2,i,j) )
          call CHECK( __LINE__, phi(KE-3,i,j) )
#endif
          diff(KE-1,i,j,ZDIR) = ( - CNZ3(2,KE-1,2) * phi(KE-1,i,j) &
                                  + CNZ3(3,KE-1,2) * phi(KE-2,i,j) &
                                  - CNZ3(1,KE-2,2) * phi(KE-3,i,j) )
          diff(KE  ,i,j,ZDIR) = ( + CNZ3(1,KE  ,2) * phi(KE-1,i,j) &
                                  + CNZ3(3,KE  ,2) * phi(KE-1,i,j) &
                                  - CNZ3(1,KE-1,2) * phi(KE-2,i,j) )
       end do
       end do
    end if
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       diff(KE  +KO,i,j,ZDIR) = - diff(KE-1+KO,i,j,ZDIR)
       diff(KE+1+KO,i,j,ZDIR) = - diff(KE-2+KO,i,j,ZDIR)
    enddo
    enddo
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS-IO, IE-IO
    do k = KS, KE-KO
#ifdef DEBUG
       call CHECK( __LINE__, phi(k,i+2,j) )
       call CHECK( __LINE__, phi(k,i+1,j) )
       call CHECK( __LINE__, phi(k,i  ,j) )
       call CHECK( __LINE__, phi(k,i-1,j) )
#endif
       diff(k,i+IO,j,XDIR) = ( + CNX3(1,i+1,1) * phi(k,i+2,j) &
                               - CNX3(2,i+1,1) * phi(k,i+1,j) &
                               + CNX3(3,i+1,1) * phi(k,i  ,j) &
                               - CNX3(1,i  ,1) * phi(k,i-1,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       diff(   1:KS-1,i,j,XDIR) = 0.0_RP
       diff(KE+1:KA  ,i,j,XDIR) = 0.0_RP
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-JO, JE-JO
    do i = IS, IE
    do k = KS, KE-KO
#ifdef DEBUG
       call CHECK( __LINE__, phi(k,i,j+2) )
       call CHECK( __LINE__, phi(k,i,j+1) )
       call CHECK( __LINE__, phi(k,i,j  ) )
       call CHECK( __LINE__, phi(k,i,j-1) )
#endif
       diff(k,i,j+JO,YDIR) = ( + CNY3(1,j+1,1) * phi(k,i,j+2) &
                               - CNY3(2,j+1,1) * phi(k,i,j+1) &
                               + CNY3(3,j+1,1) * phi(k,i,j  ) &
                               - CNY3(1,j  ,1) * phi(k,i,j-1) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       diff(   1:KS-1,i,j,YDIR) = 0.0_RP
       diff(KE+1:KA  ,i,j,YDIR) = 0.0_RP
    enddo
    enddo

    return
  end subroutine calc_diff3

  !-----------------------------------------------------------------------------
  subroutine calc_diff4( &
       num_diff_pt1, &
       num_diff_pt0, &
       CNZ4,         &
       CNX4,         &
       CNY4,         &
       k1            )
    implicit none

    real(RP), intent(out) :: num_diff_pt1(KA,IA,JA,3)
    real(RP), intent(in)  :: num_diff_pt0(KA,IA,JA,3)
    real(RP), intent(in)  :: CNZ4(5,KA)
    real(RP), intent(in)  :: CNX4(5,IA)
    real(RP), intent(in)  :: CNY4(5,JA)
    integer,  intent(in)  :: k1

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, CNZ4(1,k) )
       call CHECK( __LINE__, CNZ4(2,k) )
       call CHECK( __LINE__, CNZ4(3,k) )
       call CHECK( __LINE__, CNZ4(4,k) )
       call CHECK( __LINE__, CNZ4(5,k) )
       call CHECK( __LINE__, num_diff_pt0(k+2,i,j,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k+1,i,j,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k  ,i,j,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k-1,i,j,ZDIR) )
       call CHECK( __LINE__, num_diff_pt0(k-2,i,j,ZDIR) )
#endif
       num_diff_pt1(k,i,j,ZDIR) = &
                     ( CNZ4(1,k) * num_diff_pt0(k+2,i,j,ZDIR) &
                     - CNZ4(2,k) * num_diff_pt0(k+1,i,j,ZDIR) &
                     + CNZ4(3,k) * num_diff_pt0(k  ,i,j,ZDIR) &
                     - CNZ4(4,k) * num_diff_pt0(k-1,i,j,ZDIR) &
                     + CNZ4(5,k) * num_diff_pt0(k-2,i,j,ZDIR) )
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff_pt1(KS-1,i,j,ZDIR) = - num_diff_pt1(KS  ,i,j,ZDIR)
       num_diff_pt1(KS-2,i,j,ZDIR) = - num_diff_pt1(KS+1,i,j,ZDIR)
       num_diff_pt1(KE  ,i,j,ZDIR) = - num_diff_pt1(KE-1,i,j,ZDIR)
       num_diff_pt1(KE+1,i,j,ZDIR) = - num_diff_pt1(KE-2,i,j,ZDIR)
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, K1
#ifdef DEBUG
       call CHECK( __LINE__, CNX4(1,i) )
       call CHECK( __LINE__, CNX4(2,i) )
       call CHECK( __LINE__, CNX4(3,i) )
       call CHECK( __LINE__, CNX4(4,i) )
       call CHECK( __LINE__, CNX4(5,i) )
       call CHECK( __LINE__, num_diff_pt0(k,i-2,j,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i+1,j,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i  ,j,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i-1,j,XDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i-2,j,XDIR) )
#endif
       num_diff_pt1(k,i,j,XDIR) = &
                    ( CNX4(1,i) * num_diff_pt0(k,i+2,j,XDIR) &
                    - CNX4(2,i) * num_diff_pt0(k,i+1,j,XDIR) &
                    + CNX4(3,i) * num_diff_pt0(k,i  ,j,XDIR) &
                    - CNX4(4,i) * num_diff_pt0(k,i-1,j,XDIR) &
                    + CNX4(5,i) * num_diff_pt0(k,i-2,j,XDIR) )
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, K1
#ifdef DEBUG
       call CHECK( __LINE__, CNY4(1,j) )
       call CHECK( __LINE__, CNY4(2,j) )
       call CHECK( __LINE__, CNY4(3,j) )
       call CHECK( __LINE__, CNY4(4,j) )
       call CHECK( __LINE__, CNY4(5,j) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j-2,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j+1,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j  ,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j-1,YDIR) )
       call CHECK( __LINE__, num_diff_pt0(k,i,j-2,YDIR) )
#endif
       num_diff_pt1(k,i,j,YDIR) = &
                    ( CNY4(1,j) * num_diff_pt0(k,i,j+2,YDIR) &
                    - CNY4(2,j) * num_diff_pt0(k,i,j+1,YDIR) &
                    + CNY4(3,j) * num_diff_pt0(k,i,j  ,YDIR) &
                    - CNY4(4,j) * num_diff_pt0(k,i,j-1,YDIR) &
                    + CNY4(5,j) * num_diff_pt0(k,i,j-2,YDIR) )
    enddo
    enddo
    enddo

    return
  end subroutine calc_diff4

  !-----------------------------------------------------------------------------
  !> Flux Correction Transport Limiter
  subroutine ATMOS_DYN_fct( &
       qflx_anti,           &
       phi_in, DENS0, DENS, &
       qflx_hi, qflx_lo,    &
       mflx_hi,             &
       rdz, rdx, rdy,       &
       GSQRT, MAPF, dt,     &
       flag_vect )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       IUNDEF => CONST_UNDEF2, &
       EPSILON => CONST_EPS
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_gridtrans, only: &
       I_XYZ
    implicit none

    real(RP), intent(out) :: qflx_anti(KA,IA,JA,3)

    real(RP), intent(in) :: phi_in(KA,IA,JA) ! physical quantity
    real(RP), intent(in) :: DENS0(KA,IA,JA)
    real(RP), intent(in) :: DENS (KA,IA,JA)

    real(RP), intent(in) :: qflx_hi(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_lo(KA,IA,JA,3)
    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3)

    real(RP), intent(in) :: RDZ(:)
    real(RP), intent(in) :: RDX(:)
    real(RP), intent(in) :: RDY(:)

    real(RP), intent(in) :: GSQRT(KA,IA,JA) !< vertical metrics {G}^1/2
    real(RP), intent(in) :: MAPF(IA,JA,2)   !< map factor

    real(RP), intent(in) :: dt

    logical, intent(in) :: flag_vect

    ! work for FCT
    real(RP) :: phi_lo(KA,IA,JA)
    real(RP) :: pjpls(KA,IA,JA)
    real(RP) :: pjmns(KA,IA,JA)
    real(RP) :: qjpls(KA,IA,JA)
    real(RP) :: qjmns(KA,IA,JA)
    real(RP) :: rjpls(KA,IA,JA)
    real(RP) :: rjmns(KA,IA,JA)

    real(RP) :: qmin, qmax
    real(RP) :: zerosw, dirsw

    real(RP) :: rw, ru, rv
    real(RP) :: qa_in, qb_in
    real(RP) :: qa_lo, qb_lo
    real(RP) :: x, y

    integer :: k, i, j, ijs
    integer :: IIS, IIE, JJS, JJE
    !---------------------------------------------------------------------------

#ifdef DEBUG
    qflx_anti(:,:,:,:) = UNDEF

    pjpls(:,:,:) = UNDEF
    pjmns(:,:,:) = UNDEF
    qjpls(:,:,:) = UNDEF
    qjmns(:,:,:) = UNDEF
    rjpls(:,:,:) = UNDEF
    rjmns(:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
          phi_lo(k,i,j) = ( phi_in(k,i,j) * DENS0(k,i,j) &
                          + dt * ( - ( ( qflx_lo(k,i,j,ZDIR)-qflx_lo(k-1,i  ,j  ,ZDIR) ) * RDZ(k) &
                                     + ( qflx_lo(k,i,j,XDIR)-qflx_lo(k  ,i-1,j  ,XDIR) ) * RDX(i) &
                                     + ( qflx_lo(k,i,j,YDIR)-qflx_lo(k  ,i  ,j-1,YDIR) ) * RDY(j) &
                                     ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)                 ) &
                          ) / DENS(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net incoming quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
          pjpls(k,i,j) = dt * ( ( max(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) - min(0.0_RP,qflx_anti(k,i,j,ZDIR)) ) * RDZ(k) &
                              + ( max(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) - min(0.0_RP,qflx_anti(k,i,j,XDIR)) ) * RDX(i) &
                              + ( max(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) - min(0.0_RP,qflx_anti(k,i,j,YDIR)) ) * RDY(j) &
                              ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net outgoing quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
          pjmns(k,i,j) = dt * ( ( max(0.0_RP,qflx_anti(k,i,j,ZDIR)) - min(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) ) * RDZ(k) &
                              + ( max(0.0_RP,qflx_anti(k,i,j,XDIR)) - min(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) ) * RDX(i) &
                              + ( max(0.0_RP,qflx_anti(k,i,j,YDIR)) - min(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) ) * RDY(j) &
                              ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc allowable range or quantity change by antidiffusive flux

       if (flag_vect) then

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
          rw = (mflx_hi(k,i,j,ZDIR)+mflx_hi(k-1,i  ,j  ,ZDIR)) * RDZ(k) ! 2 * rho * w / dz
          ru = (mflx_hi(k,i,j,XDIR)+mflx_hi(k  ,i-1,j  ,XDIR)) * RDX(i) ! 2 * rho * u / dx
          rv = (mflx_hi(k,i,j,YDIR)+mflx_hi(k  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy
          if ( abs(ru) < EPSILON .and. abs(rv) < EPSILON .and. abs(rw) < EPSILON ) then
             qa_in = phi_in(k,i,j)
             qa_lo = phi_lo(k,i,j)
          elseif ( abs(ru) .ge. abs(rv) .and. abs(ru) .ge. abs(rw) ) then
             x = rv / ru
             y = rw / ru
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j+1) * x ) * y
                   qb_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j+1) * x ) * y
                   qb_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j-1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j-1) * x ) * y
                   qb_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j-1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j-1) * x ) * y
                   qb_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j-1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j+1) * x ) * y
                end if
             end if
          else if ( abs(rv) .ge. abs(ru) .and. abs(rv) .ge. abs(rw) ) then
             x = rw / rv
             y = ru / rv
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j+1) * x ) * y
                   qb_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j+1) * x ) * y
                   qb_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j-1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j+1) * x ) * y
                   qb_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j+1) * x ) * y
                   qb_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j-1) * x ) * y
                end if
             end if
          else
             x = ru / rw
             y = rv / rw
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j+1) * x ) * y
                   qb_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j-1) * x ) * y
                   qb_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j-1) * x ) * y
                   qb_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j+1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j+1) * x ) * y
                   qb_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j+1) * x ) * y
                   qb_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j-1) * x ) * y
                   qb_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j-1) * x ) * y
                   qb_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j+1) * x ) * y
                end if
             end if
          end if
          qmax = max( &
               phi_in(k,i,j), qa_in, qb_in, &
               phi_lo(k,i,j), qa_lo, qb_lo  )
          qmin = min( &
               phi_in(k,i,j), qa_in, qb_in, &
               phi_lo(k,i,j), qa_lo, qb_lo  )
          qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
          qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
       end do
       end do
       end do

       do j = JJS, JJE
       do i = IIS, IIE
          k = KS
          rw = mflx_hi(k,i,j,ZDIR) * RDZ(k) ! rho * w / dz
          ru = (mflx_hi(k,i,j,XDIR)+mflx_hi(k  ,i-1,j  ,XDIR)) * RDX(i) ! rho * u / dx
          rv = (mflx_hi(k,i,j,YDIR)+mflx_hi(k  ,i  ,j-1,YDIR)) * RDY(j) ! rho * v / dy
          if ( abs(ru) < EPSILON .and. abs(rv) < EPSILON .and. abs(rw) < EPSILON ) then
             qa_in = phi_in(k,i,j)
             qa_lo = phi_lo(k,i,j)
          else if ( abs(ru) .ge. abs(rv) .and. abs(ru) .ge. abs(rw) ) then
             x = rv / ru
             y = rw / ru
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j-1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j+1) * x ) * y
                end if
             end if
          else if ( abs(rv) .ge. abs(ru) .and. abs(rv) .ge. abs(rw) ) then
             x = rw / rv
             y = ru / rv
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j+1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j-1) * x ) * y
                end if
             end if
          else
             x = ru / rw
             y = rv / rw
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i+1,j-1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k+1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k+1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k+1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k+1,i-1,j-1) * x ) * y
                end if
             end if
          end if
          qmax = max( &
               phi_in(k,i,j), qa_in, &
               phi_lo(k,i,j), qa_lo )
          qmin = min( &
               phi_in(k,i,j), qa_in, &
               phi_lo(k,i,j), qa_lo )
          qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
          qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
       end do
       end do

       do j = JJS, JJE
       do i = IIS, IIE
          k = KE
          rw = mflx_hi(k-1,i,j,ZDIR) * RDZ(k) ! rho * w / dz
          ru = (mflx_hi(k,i,j,XDIR)+mflx_hi(k  ,i-1,j  ,YDIR)) * RDX(i) ! rho * u / dx
          rv = (mflx_hi(k,i,j,YDIR)+mflx_hi(k  ,i  ,j-1,YDIR)) * RDY(j) ! rho * v / dy
          if ( abs(ru) < EPSILON .and. abs(rv) < EPSILON .and. abs(rw) < EPSILON ) then
             qa_in = phi_in(k,i,j)
             qa_lo = phi_lo(k,i,j)
          else if ( abs(ru) .ge. abs(rv) .and. abs(ru) .ge. abs(rw) ) then
             x = rv / ru
             y = rw / ru
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j+1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i-1,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i-1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k  ,i+1,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i+1,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j-1) * x ) * y
                end if
             end if
          else if ( abs(rv) .ge. abs(ru) .and. abs(rv) .ge. abs(rw) ) then
             x = rw / rv
             y = ru / rv
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j-1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j-1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i+1,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j+1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_in(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k  ,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i  ,j+1) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k  ,i-1,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j+1) * x ) * y
                end if
             end if
          else
             x = ru / rw
             y = rv / rw
             if ( x .ge. 0.0_RP ) then
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k-1 ,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1 ,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i-1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i-1,j+1) * x ) * y
                end if
             else ! x < 0
                x = -x
                if ( y .ge. 0.0_RP ) then
                   qa_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j  ) * x ) * y &
                         + ( phi_in(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j-1) * x ) * y
                   qa_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j  ) * x ) * y &
                         + ( phi_lo(k-1,i  ,j-1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j-1) * x ) * y
                else ! y < 0
                   y = -y
                   qa_in = ( phi_in(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_in(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_in(k-1,i+1,j+1) * x ) * y
                   qa_lo = ( phi_lo(k-1,i  ,j  ) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j  ) * x ) * (1.0_RP-y) &
                         + ( phi_lo(k-1,i  ,j+1) * (1.0_RP-x) &
                           + phi_lo(k-1,i+1,j+1) * x ) * y
                end if
             end if
          end if
          qmax = max( &
               phi_in(k,i,j), qa_in, &
               phi_lo(k,i,j), qa_lo )
          qmin = min( &
               phi_in(k,i,j), qa_in, &
               phi_lo(k,i,j), qa_lo )
          qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
          qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
       end do
       end do

       else

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
          qmax = max( phi_in(k  ,i  ,j  ), &
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
                      phi_lo(k  ,i  ,j-1) )
          qmin = min( phi_in(k  ,i  ,j  ), &
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
          qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
          qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
          qmax = max( phi_in(KS  ,i  ,j  ), &
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
                      phi_lo(KS  ,i  ,j-1) )
          qmin = min( phi_in(KS  ,i  ,j  ), &
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
                      phi_lo(KS  ,i  ,j-1) )
          qjmns(KS,i,j) = ( phi_lo(KS,i,j) - qmin ) * DENS(KS,i,j)
          qjpls(KS,i,j) = ( qmax - phi_lo(KS,i,j) ) * DENS(KS,i,j)

          qmax = max( phi_in(KE  ,i  ,j  ), &
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
                      phi_lo(KE  ,i  ,j-1) )
          qmin = min( phi_in(KE  ,i  ,j  ), &
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
          qjpls(KE,i,j) = ( qmax - phi_lo(KE,i,j) ) * DENS(KE,i,j)
          qjmns(KE,i,j) = ( phi_lo(KE,i,j) - qmin ) * DENS(KE,i,j)
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       !--- incoming flux limitation factor [0-1]
       !$omp parallel do private(i,j,k,zerosw) OMP_SCHEDULE_ collapse(2)
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
       !$omp parallel do private(i,j,k,zerosw) OMP_SCHEDULE_ collapse(2)
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

    call COMM_vars8( rjpls(:,:,:), 1 )
    call COMM_vars8( rjmns(:,:,:), 2 )
    call COMM_wait ( rjpls(:,:,:), 1 )
    call COMM_wait ( rjmns(:,:,:), 2 )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !--- update high order flux with antidiffusive flux
       !$omp parallel do private(i,j,k,dirsw) OMP_SCHEDULE_ collapse(2)
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
                 * ( 1.0_RP &
                   - min( rjpls(k+1,i,j),rjmns(k  ,i,j) ) * (          dirsw ) &
                   - min( rjpls(k  ,i,j),rjmns(k+1,i,j) ) * ( 1.0_RP - dirsw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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

       !$omp parallel do private(i,j,k,dirsw) OMP_SCHEDULE_ collapse(2)
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
                 * ( 1.0_RP &
                   - min( rjpls(k,i+1,j),rjmns(k,i  ,j) ) * (          dirsw ) &
                   - min( rjpls(k,i  ,j),rjmns(k,i+1,j) ) * ( 1.0_RP - dirsw ) )
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
       !$omp parallel do private(i,j,k,dirsw) OMP_SCHEDULE_ collapse(2)
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
                 * ( 1.0_RP &
                   - min( rjpls(k,i,j+1),rjmns(k,i,j  ) ) * (          dirsw ) &
                   - min( rjpls(k,i,j  ),rjmns(k,i,j+1) ) * ( 1.0_RP - dirsw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    return
  end subroutine ATMOS_DYN_fct

end module scale_atmos_dyn_common
