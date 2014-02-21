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
  public :: ATMOS_DYN_numfilter_setup
  public :: ATMOS_DYN_numfilter_coef
  public :: ATMOS_DYN_numfilter_coef_q
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
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_numfilter_setup( &
       DIFF4,                              &
       CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4, &
       CDZ, CDX, CDY, FDZ, FDX, FDY,       &
       ND_ORDER, ND_COEF,                  &
       DT                                  )
    implicit none

    real(RP), intent(out) :: DIFF4
    real(RP), intent(out) :: CNZ3(3,KA,2)
    real(RP), intent(out) :: CNX3(3,IA,2)
    real(RP), intent(out) :: CNY3(3,JA,2)
    real(RP), intent(out) :: CNZ4(5,KA,2)
    real(RP), intent(out) :: CNX4(5,IA,2)
    real(RP), intent(out) :: CNY4(5,JA,2)

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)
    real(RP), intent(in)  :: FDZ(KA-1)
    real(RP), intent(in)  :: FDX(IA-1)
    real(RP), intent(in)  :: FDY(JA-1)

    integer,  intent(in)  :: ND_ORDER
    real(RP), intent(in)  :: ND_COEF
    real(RP), intent(in)  :: DT

    integer :: k, i, j
    !---------------------------------------------------------------------------

#ifdef DEBUG
    CNZ3(:,:,:) = UNDEF
    CNX3(:,:,:) = UNDEF
    CNY3(:,:,:) = UNDEF
    CNZ4(:,:,:) = UNDEF
    CNX4(:,:,:) = UNDEF
    CNY4(:,:,:) = UNDEF
#endif

    ! numerical diffusion
    if ( ND_COEF > 0.0_RP ) then
       DIFF4 = ND_COEF / ( 2**(4*ND_ORDER) * DT )

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
    else
       DIFF4 = 0.0_RP

       CNZ3(:,:,:) = 0.0_RP
       CNX3(:,:,:) = 0.0_RP
       CNY3(:,:,:) = 0.0_RP
       CNZ4(:,:,:) = 0.0_RP
       CNX4(:,:,:) = 0.0_RP
       CNY4(:,:,:) = 0.0_RP
    endif

    return
  end subroutine ATMOS_DYN_numfilter_setup

  !-----------------------------------------------------------------------------
  !> Calc coefficient of numerical filter
  subroutine ATMOS_DYN_numfilter_coef( &
       num_diff,                               &
       DENS, MOMZ, MOMX, MOMY, RHOT,           &
       CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4,     &
       CDZ, CDX, CDY, FDZ, FDX, FDY,           &
       REF_dens, REF_pott,                     &
       DIFF4, ND_ORDER, ND_SFC_FACT, ND_USE_RS )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out), target :: num_diff(KA,IA,JA,5,3)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)

    real(RP), intent(in)  :: CNZ3(3,KA,2)
    real(RP), intent(in)  :: CNX3(3,IA,2)
    real(RP), intent(in)  :: CNY3(3,JA,2)
    real(RP), intent(in)  :: CNZ4(5,KA,2)
    real(RP), intent(in)  :: CNX4(5,IA,2)
    real(RP), intent(in)  :: CNY4(5,JA,2)

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)
    real(RP), intent(in)  :: FDZ(KA-1)
    real(RP), intent(in)  :: FDX(IA-1)
    real(RP), intent(in)  :: FDY(JA-1)

    real(RP), intent(in)  :: REF_dens(KA,IA,JA)
    real(RP), intent(in)  :: REF_pott(KA,IA,JA)

    real(RP), intent(in)  :: DIFF4
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

    real(RP), target  :: num_diff_work(KA,IA,JA,5,3)
    real(RP), pointer :: num_diff_pt0 (:,:,:,:,:)
    real(RP), pointer :: num_diff_pt1 (:,:,:,:,:)
    real(RP), pointer :: tmp_pt       (:,:,:,:,:)

    integer  :: nd_order4, no

    integer :: IIS, IIE
    integer :: JJS, JJE
    integer :: k, i, j
    !---------------------------------------------------------------------------

    POTT(:,:,:) = RHOT(:,:,:) / DENS(:,:,:)

    !###########################################################################
    ! 1st order coefficients
    !###########################################################################

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !-----< Density & Potential temperature >-----

       if ( ND_USE_RS ) then
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             dens_diff(k,i,j) = DENS(k,i,j) - REF_dens(k,i,j)
             pott_diff(k,i,j) = POTT(k,i,j) - REF_pott(k,i,j)
          enddo
          enddo
          enddo
       else
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
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

          do j  = JJS, JJE
          do i  = IIS, IIE
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
          enddo
          enddo
       endif

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
          num_diff(k,i,j,I_DENS,ZDIR) = ( + CNZ3(1,k+1,1) * dens_diff(k+2,i,j) &
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
          num_diff(KS  ,i,j,I_DENS,ZDIR) = ( + CNZ3(1,KS+1,1) * dens_diff(KS+2,i,j) &
                                             - CNZ3(2,KS+1,1) * dens_diff(KS+1,i,j) &
                                             + CNZ3(3,KS+1,1) * dens_diff(KS  ,i,j) &
                                             - CNZ3(1,KS  ,1) * dens_diff(KS+1,i,j) )
          num_diff(KS-1,i,j,I_DENS,ZDIR) = - num_diff(KS  ,i,j,I_DENS,ZDIR)
          num_diff(KS-2,i,j,I_DENS,ZDIR) = - num_diff(KS+1,i,j,I_DENS,ZDIR)
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE,1) )
          call CHECK( __LINE__, CNZ3(2,KE,1) )
          call CHECK( __LINE__, CNZ3(3,KE,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, dens_diff(KE,i,j) )
          call CHECK( __LINE__, dens_diff(KE-1,i,j) )
          call CHECK( __LINE__, dens_diff(KE-2,i,j) )
#endif
          num_diff(KE-1,i,j,I_DENS,ZDIR) = ( + CNZ3(1,KE  ,1) * dens_diff(KE-1,i,j) &
                                             - CNZ3(2,KE  ,1) * dens_diff(KE  ,i,j) &
                                             + CNZ3(3,KE  ,1) * dens_diff(KE-1,i,j) &
                                             - CNZ3(1,KE-1,1) * dens_diff(KE-2,i,j) )
          num_diff(KE  ,i,j,I_DENS,ZDIR) = - num_diff(KE-1,i,j,I_DENS,ZDIR)
          num_diff(KE+1,i,j,I_DENS,ZDIR) = - num_diff(KE-2,i,j,I_DENS,ZDIR)
          num_diff(KE+2,i,j,I_DENS,ZDIR) = 0.0_RP
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
          call CHECK( __LINE__, DENS(k,i+2,j) )
          call CHECK( __LINE__, DENS(k,i+1,j) )
          call CHECK( __LINE__, DENS(k,i  ,j) )
          call CHECK( __LINE__, DENS(k,i-1,j) )
#endif
          num_diff(k,i,j,I_DENS,XDIR) = ( + CNX3(1,i+1,1) * DENS(k,i+2,j) &
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
          call CHECK( __LINE__, DENS(k,i,j+2) )
          call CHECK( __LINE__, DENS(k,i,j+1) )
          call CHECK( __LINE__, DENS(k,i,j  ) )
          call CHECK( __LINE__, DENS(k,i,j-1) )
#endif
          num_diff(k,i,j,I_DENS,YDIR) = ( + CNY3(1,j+1,1) * DENS(k,i,j+2) &
                                          - CNY3(2,j+1,1) * DENS(k,i,j+1) &
                                          + CNY3(3,j+1,1) * DENS(k,i,j  ) &
                                          - CNY3(1,j  ,1) * DENS(k,i,j-1) )
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          num_diff(   1:KS-1,i,j,I_DENS,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_DENS,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< z-momentum >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
       do k = KS, KE-1
          VELZ(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo

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
          num_diff(k,i,j,I_MOMZ,ZDIR) = ( + CNZ3(1,k  ,2) * VELZ(k+1,i,j) &
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
          call CHECK( __LINE__, CNZ3(1,KS+1,2) )
          call CHECK( __LINE__, CNZ3(2,KS+1,2) )
          call CHECK( __LINE__, CNZ3(3,KS+1,2) )
          call CHECK( __LINE__, VELZ(KS+2,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
#endif
          num_diff(KS+1,i,j,I_MOMZ,ZDIR) = ( + CNZ3(1,KS+1,2) * VELZ(KS+2,i,j) &
                                             - CNZ3(2,KS+1,2) * VELZ(KS+1,i,j) &
                                             + CNZ3(3,KS+1,2) * VELZ(KS  ,i,j) )
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KS,2) )
          call CHECK( __LINE__, CNZ3(2,KS,2) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
#endif
          num_diff(KS  ,i,j,I_MOMZ,ZDIR) = ( + CNZ3(1,KS  ,2) * VELZ(KS+1,i,j) &
                                             - CNZ3(2,KS  ,2) * VELZ(KS  ,i,j) &
                                             - CNZ3(1,KS-1,2) * VELZ(KS+1,i,j) )
          num_diff(KS-1,i,j,I_MOMZ,ZDIR) = - num_diff(KS  ,i,j,I_MOMZ,ZDIR)
          num_diff(KS-2,i,j,I_MOMZ,ZDIR) = - num_diff(KS+1,i,j,I_MOMZ,ZDIR)
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(2,KE-1,2) )
          call CHECK( __LINE__, CNZ3(3,KE-1,2) )
          call CHECK( __LINE__, CNZ3(1,KE-2,2) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-3,i,j) )
#endif
          num_diff(KE-1,i,j,I_MOMZ,ZDIR) = ( - CNZ3(2,KE-1,2) * VELZ(KE-1,i,j) &
                                             + CNZ3(3,KE-1,2) * VELZ(KE-2,i,j) &
                                             - CNZ3(1,KE-2,2) * VELZ(KE-3,i,j) )
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(3,KE,2) )
          call CHECK( __LINE__, CNZ3(1,KE-1,2) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
#endif
          num_diff(KE  ,i,j,I_MOMZ,ZDIR) = ( + CNZ3(1,KE  ,2) * VELZ(KE-1,i,j) &
                                             + CNZ3(3,KE  ,2) * VELZ(KE-1,i,j) &
                                             - CNZ3(1,KE-1,2) * VELZ(KE-2,i,j) )
          num_diff(KE+1,i,j,I_MOMZ,ZDIR) = - num_diff(KE,i,j,I_MOMZ,ZDIR)
          num_diff(KE+2,i,j,I_MOMZ,ZDIR) = - num_diff(KE-1,i,j,I_MOMZ,ZDIR)
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
          num_diff(k,i,j,I_MOMZ,XDIR) = ( + CNX3(1,i+1,1) * VELZ(k,i+2,j) &
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
          num_diff(k,i,j,I_MOMZ,YDIR) = ( + CNY3(1,j+1,1) * VELZ(k,i,j+2) &
                                          - CNY3(2,j+1,1) * VELZ(k,i,j+1) &
                                          + CNY3(3,j+1,1) * VELZ(k,i,j  ) &
                                          - CNY3(1,j  ,1) * VELZ(k,i,j-1) )
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          num_diff( 1:KS-1,i,j,I_MOMZ,XDIR) = 0.0_RP
          num_diff(KE:KA  ,i,j,I_MOMZ,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< x-momentum >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+1
       do k = KS, KE
          VELX(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo

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
          num_diff(k,i,j,I_MOMX,ZDIR) = ( + CNZ3(1,k+1,1) * VELX(k+2,i,j) &
                                          - CNZ3(2,k+1,1) * VELX(k+1,i,j) &
                                          + CNZ3(3,k+1,1) * VELX(k  ,i,j) &
                                          - CNZ3(1,k  ,1) * VELX(k-1,i,j) )
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
          call CHECK( __LINE__, CNZ3(1,KS  ,1) )
          call CHECK( __LINE__, VELX(KS+2,i,j) )
          call CHECK( __LINE__, VELX(KS+1,i,j) )
          call CHECK( __LINE__, VELX(KS  ,i,j) )
#endif
          num_diff(KS  ,i,j,I_MOMX,ZDIR) = ( + CNZ3(1,KS+1,1) * VELX(KS+2,i,j) &
                                             - CNZ3(2,KS+1,1) * VELX(KS+1,i,j) &
                                             + CNZ3(3,KS+1,1) * VELX(KS  ,i,j) &
                                             - CNZ3(1,KS  ,1) * VELX(KS+1,i,j) )
          num_diff(KS-1,i,j,I_MOMX,ZDIR) = - num_diff(KS  ,i,j,I_MOMX,ZDIR)
          num_diff(KS-2,i,j,I_MOMX,ZDIR) = - num_diff(KS+1,i,j,I_MOMX,ZDIR)
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE  ,1) )
          call CHECK( __LINE__, CNZ3(2,KE  ,1) )
          call CHECK( __LINE__, CNZ3(3,KE  ,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, VELX(KE  ,i,j) )
          call CHECK( __LINE__, VELX(KE-1,i,j) )
          call CHECK( __LINE__, VELX(KE-2,i,j) )
#endif
          num_diff(KE-1,i,j,I_MOMX,ZDIR) = ( + CNZ3(1,KE  ,1) * VELX(KE-1,i,j) &
                                             - CNZ3(2,KE  ,1) * VELX(KE  ,i,j) &
                                             + CNZ3(3,KE  ,1) * VELX(KE-1,i,j) &
                                             - CNZ3(1,KE-1,1) * VELX(KE-2,i,j) )
          num_diff(KE  ,i,j,I_MOMX,ZDIR) = - num_diff(KE-1,i,j,I_MOMX,ZDIR)
          num_diff(KE+1,i,j,I_MOMX,ZDIR) = - num_diff(KE-2,i,j,I_MOMX,ZDIR)
          num_diff(KE+2,i,j,I_MOMX,ZDIR) = 0.0_RP
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
          num_diff(k,i,j,I_MOMX,XDIR) = ( + CNX3(1,i  ,2) * VELX(k,i+1,j) &
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
          num_diff(k,i,j,I_MOMX,YDIR) = ( + CNY3(1,j+1,1) * VELX(k,i,j+2) &
                                          - CNY3(2,j+1,1) * VELX(k,i,j+1) &
                                          + CNY3(3,j+1,1) * VELX(k,i,j  ) &
                                          - CNY3(1,j  ,1) * VELX(k,i,j-1) )
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          num_diff(   1:KS-1,i,j,I_MOMX,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_MOMX,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< y-momentum >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+1
       do i = IIS-2, IIE+2
       do k = KS, KE
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
          call CHECK( __LINE__, CNZ3(1,k  ,1) )
          call CHECK( __LINE__, VELY(k+2,i,j) )
          call CHECK( __LINE__, VELY(k+1,i,j) )
          call CHECK( __LINE__, VELY(k  ,i,j) )
          call CHECK( __LINE__, VELY(k-1,i,j) )
#endif
          num_diff(k,i,j,I_MOMY,ZDIR) = ( + CNZ3(1,k+1,1) * VELY(k+2,i,j) &
                                          - CNZ3(2,k+1,1) * VELY(k+1,i,j) &
                                          + CNZ3(3,k+1,1) * VELY(k  ,i,j) &
                                          - CNZ3(1,k  ,1) * VELY(k-1,i,j) )
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
          call CHECK( __LINE__, CNZ3(1,KS  ,1) )
          call CHECK( __LINE__, VELY(KS+2,i,j) )
          call CHECK( __LINE__, VELY(KS+1,i,j) )
          call CHECK( __LINE__, VELY(KS,i,j) )
#endif
          num_diff(KS  ,i,j,I_MOMY,ZDIR) = ( + CNZ3(1,KS+1,1) * VELY(KS+2,i,j) &
                                             - CNZ3(2,KS+1,1) * VELY(KS+1,i,j) &
                                             + CNZ3(3,KS+1,1) * VELY(KS  ,i,j) &
                                             - CNZ3(1,KS  ,1) * VELY(KS+1,i,j) )
          num_diff(KS-1,i,j,I_MOMY,ZDIR) = - num_diff(KS  ,i,j,I_MOMY,ZDIR)
          num_diff(KS-2,i,j,I_MOMY,ZDIR) = - num_diff(KS+1,i,j,I_MOMY,ZDIR)
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE,1) )
          call CHECK( __LINE__, CNZ3(2,KE,1) )
          call CHECK( __LINE__, CNZ3(3,KE,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, VELY(KE,i,j) )
          call CHECK( __LINE__, VELY(KE-1,i,j) )
          call CHECK( __LINE__, VELY(KE-2,i,j) )
#endif
          num_diff(KE-1,i,j,I_MOMY,ZDIR) = ( + CNZ3(1,KE  ,1) * VELY(KE-1,i,j) &
                                             - CNZ3(2,KE  ,1) * VELY(KE  ,i,j) &
                                             + CNZ3(3,KE  ,1) * VELY(KE-1,i,j) &
                                             - CNZ3(1,KE-1,1) * VELY(KE-2,i,j) )
          num_diff(KE  ,i,j,I_MOMY,ZDIR) = - num_diff(KE-1,i,j,I_MOMY,ZDIR)
          num_diff(KE+1,i,j,I_MOMY,ZDIR) = - num_diff(KE-2,i,j,I_MOMY,ZDIR)
          num_diff(KE+2,i,j,I_MOMY,ZDIR) = 0.0_RP
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
          num_diff(k,i,j,I_MOMY,XDIR) = ( + CNX3(1,i+1,1) * VELY(k,i+2,j) &
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
          num_diff(k,i,j,I_MOMY,YDIR) = ( + CNY3(1,j  ,2) * VELY(k,i,j+1) &
                                          - CNY3(2,j  ,2) * VELY(k,i,j  ) &
                                          + CNY3(3,j  ,2) * VELY(k,i,j-1) &
                                          - CNY3(1,j-1,2) * VELY(k,i,j-2) )
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          num_diff(   1:KS-1,i,j,I_MOMY,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_MOMY,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< rho * theta >-----

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
          num_diff(k,i,j,I_RHOT,ZDIR) = ( + CNZ3(1,k+1,1) * pott_diff(k+2,i,j) &
                                          - CNZ3(2,k+1,1) * pott_diff(k+1,i,j) &
                                          + CNZ3(3,k+1,1) * pott_diff(k  ,i,j) &
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
          num_diff(KS  ,i,j,I_RHOT,ZDIR) = ( + CNZ3(1,KS+1,1) * pott_diff(KS+2,i,j) &
                                             - CNZ3(2,KS+1,1) * pott_diff(KS+1,i,j) &
                                             + CNZ3(3,KS+1,1) * pott_diff(KS  ,i,j) &
                                             - CNZ3(1,KS  ,1) * pott_diff(KS+1,i,j) )
          num_diff(KS-1,i,j,I_RHOT,ZDIR) = - num_diff(KS  ,i,j,I_RHOT,ZDIR)
          num_diff(KS-2,i,j,I_RHOT,ZDIR) = - num_diff(KS+1,i,j,I_RHOT,ZDIR)
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, CNZ3(1,KE,1) )
          call CHECK( __LINE__, CNZ3(2,KE,1) )
          call CHECK( __LINE__, CNZ3(3,KE,1) )
          call CHECK( __LINE__, CNZ3(1,KE-1,1) )
          call CHECK( __LINE__, pott_diff(KE,i,j) )
          call CHECK( __LINE__, pott_diff(KE-1,i,j) )
          call CHECK( __LINE__, pott_diff(KE-2,i,j) )
#endif
          num_diff(KE-1,i,j,I_RHOT,ZDIR) = ( + CNZ3(1,KE  ,1) * pott_diff(KE  ,i,j) &
                                             - CNZ3(2,KE  ,1) * pott_diff(KE  ,i,j) &
                                             + CNZ3(3,KE  ,1) * pott_diff(KE-1,i,j) &
                                             - CNZ3(1,KE-1,1) * pott_diff(KE-2,i,j) )
          num_diff(KE  ,i,j,I_RHOT,ZDIR) = - num_diff(KE-1,i,j,I_RHOT,ZDIR)
          num_diff(KE+1,i,j,I_RHOT,ZDIR) = - num_diff(KE-2,i,j,I_RHOT,ZDIR)
          num_diff(KE+2,i,j,I_RHOT,ZDIR) = 0.0_RP
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
          call CHECK( __LINE__, POTT(k,i+2,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
#endif
          num_diff(k,i,j,I_RHOT,XDIR) = ( + CNX3(1,i+1,1) * POTT(k,i+2,j) &
                                          - CNX3(2,i+1,1) * POTT(k,i+1,j) &
                                          + CNX3(3,i+1,1) * POTT(k,i  ,j) &
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
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
#endif
          num_diff(k,i,j,I_RHOT,YDIR) = ( + CNY3(1,j+1,1) * POTT(k,i,j+2) &
                                          - CNY3(2,j+1,1) * POTT(k,i,j+1) &
                                          + CNY3(3,j+1,1) * POTT(k,i,j  ) &
                                          - CNY3(1,j  ,1) * POTT(k,i,j-1) )
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          num_diff(   1:KS-1,i,j,I_RHOT,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_RHOT,YDIR) = 0.0_RP
       enddo
       enddo

    enddo
    enddo ! end tile

    !###########################################################################
    ! High order coefficients
    !###########################################################################

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

       call COMM_wait ( num_diff_pt0(:,:,:,I_DENS,ZDIR),  1 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_DENS,XDIR),  2 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_DENS,YDIR),  3 )

       call calc_numdiff4( num_diff_pt1(:,:,:,:,:), & ! (out)
                           num_diff_pt0(:,:,:,:,:), & ! (in)
                           CNZ4(:,:,1),             & ! (in)
                           CNX4(:,:,1),             & ! (in)
                           CNY4(:,:,1),             & ! (in)
                           I_DENS,                  & ! (in)
                           KE                       ) ! (in)

       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMZ,ZDIR),  4 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMZ,XDIR),  5 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMZ,YDIR),  6 )

       call calc_numdiff4( num_diff_pt1(:,:,:,:,:), & ! (out)
                           num_diff_pt0(:,:,:,:,:), & ! (in)
                           CNZ4(:,:,2),             & ! (in)
                           CNX4(:,:,1),             & ! (in)
                           CNY4(:,:,1),             & ! (in)
                           I_MOMZ,                  & ! (in)
                           KE-1                     ) ! (in)

       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMX,ZDIR),  7 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMX,XDIR),  8 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMX,YDIR),  9 )

       call calc_numdiff4( num_diff_pt1(:,:,:,:,:), & ! (out)
                           num_diff_pt0(:,:,:,:,:), & ! (in)
                           CNZ4(:,:,1),             & ! (in)
                           CNX4(:,:,2),             & ! (in)
                           CNY4(:,:,1),             & ! (in)
                           I_MOMX,                  & ! (in)
                           KE                       ) ! (in)

       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMY,ZDIR), 10 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMY,XDIR), 11 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_MOMY,YDIR), 12 )

       call calc_numdiff4( num_diff_pt1(:,:,:,:,:), & ! (out)
                           num_diff_pt0(:,:,:,:,:), & ! (in)
                           CNZ4(:,:,1),             & ! (in)
                           CNX4(:,:,1),             & ! (in)
                           CNY4(:,:,2),             & ! (in)
                           I_MOMY,                  & ! (in)
                           KE                       ) ! (in)

       call COMM_wait ( num_diff_pt0(:,:,:,I_RHOT,ZDIR), 13 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_RHOT,XDIR), 14 )
       call COMM_wait ( num_diff_pt0(:,:,:,I_RHOT,YDIR), 15 )

       call calc_numdiff4( num_diff_pt1(:,:,:,:,:), & ! (out)
                           num_diff_pt0(:,:,:,:,:), & ! (in)
                           CNZ4(:,:,1),             & ! (in)
                           CNX4(:,:,1),             & ! (in)
                           CNY4(:,:,1),             & ! (in)
                           I_RHOT,                  & ! (in)
                           KE                       ) ! (in)

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

       !-----< density >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS-1, KE
             num_diff(k,i,j,I_DENS,ZDIR) = num_diff_pt0(k,i,j,I_DENS,ZDIR) * DIFF4 * CDZ(k)**nd_order4
          enddo
          num_diff(   1:KS-2,i,j,I_DENS,ZDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_DENS,ZDIR) = 0.0_RP

          do k = KS+2, KE
             num_diff(k,i,j,I_DENS,XDIR) = num_diff_pt0(k,i,j,I_DENS,XDIR) * DIFF4 * CDX(i)**nd_order4
          enddo
          num_diff(KS  ,i,j,I_DENS,XDIR) = num_DIFF(KS  ,i,j,I_DENS,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_DENS,XDIR) = num_DIFF(KS+1,i,j,I_DENS,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_DENS,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_DENS,XDIR) = 0.0_RP

          do k = KS+2, KE
             num_diff(k,i,j,I_DENS,YDIR) = num_diff_pt0(k,i,j,I_DENS,YDIR) * DIFF4 * CDY(j)**nd_order4
          enddo
          num_diff(KS  ,i,j,I_DENS,YDIR) = num_DIFF(KS  ,i,j,I_DENS,YDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_DENS,YDIR) = num_DIFF(KS+1,i,j,I_DENS,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_DENS,YDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_DENS,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< z-momentum >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS+1, KE-1
             num_diff(k,i,j,I_MOMZ,ZDIR) = num_diff_pt0(k,i,j,I_MOMZ,ZDIR) * DIFF4 * FDZ(k)**nd_order4 &
                                         * DENS(k,i,j)
          enddo

          num_diff( 1:KS,i,j,I_MOMZ,ZDIR) = 0.0_RP
          num_diff(KE:KA,i,j,I_MOMZ,ZDIR) = 0.0_RP

          do k = KS, KE-1
             num_diff(k,i,j,I_MOMZ,XDIR) = num_diff_pt0(k,i,j,I_MOMZ,XDIR) * DIFF4 * CDX(i)**nd_order4 &
                                         * 0.25_RP * ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          num_diff(KS  ,i,j,I_MOMZ,XDIR) = num_DIFF(KS  ,i,j,I_MOMZ,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMZ,XDIR) = num_DIFF(KS+1,i,j,I_MOMZ,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff( 1:KS-1,i,j,I_MOMZ,XDIR) = 0.0_RP
          num_diff(KE:KA  ,i,j,I_MOMZ,XDIR) = 0.0_RP

          do k = KS, KE-1
             num_diff(k,i,j,I_MOMZ,YDIR) = num_diff_pt0(k,i,j,I_MOMZ,YDIR) * DIFF4 * CDY(j)**nd_order4 &
                                         * 0.25_RP * ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          num_diff(KS  ,i,j,I_MOMZ,YDIR) = num_DIFF(KS  ,i,j,I_MOMZ,YDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMZ,YDIR) = num_DIFF(KS+1,i,j,I_MOMZ,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff( 1:KS-1,i,j,I_MOMZ,YDIR) = 0.0_RP
          num_diff(KE:KA  ,i,j,I_MOMZ,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< x-momentum >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS,    KE-1
             num_diff(k,i,j,I_MOMX,ZDIR) = num_diff_pt0(k,i,j,I_MOMX,ZDIR) * DIFF4 * CDZ(k)**nd_order4 &
                                         * 0.25_RP * ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) )
          enddo

          num_diff( 1:KS-1,i,j,I_MOMX,ZDIR) = 0.0_RP
          num_diff(KE:KA  ,i,j,I_MOMX,ZDIR) = 0.0_RP

          do k = KS, KE
             num_diff(k,i,j,I_MOMX,XDIR) = num_diff_pt0(k,i,j,I_MOMX,XDIR) * DIFF4 * FDX(i)**nd_order4 &
                                         * DENS(k,i,j)
          enddo
          num_diff(KS  ,i,j,I_MOMX,XDIR) = num_DIFF(KS  ,i,j,I_MOMX,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMX,XDIR) = num_DIFF(KS+1,i,j,I_MOMX,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_MOMX,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_MOMX,XDIR) = 0.0_RP

          do k = KS, KE
             num_diff(k,i,j,I_MOMX,YDIR) = num_diff_pt0(k,i,j,I_MOMX,YDIR) * DIFF4 * CDY(j)**nd_order4 &
                                         * 0.25_RP * ( DENS(k,i+1,j+1)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          num_diff(KS  ,i,j,I_MOMX,YDIR) = num_DIFF(KS  ,i,j,I_MOMX,YDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMX,YDIR) = num_DIFF(KS+1,i,j,I_MOMX,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_MOMX,YDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_MOMX,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< y-momentum >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS,  KE-1
             num_diff(k,i,j,I_MOMY,ZDIR) = num_diff_pt0(k,i,j,I_MOMY,ZDIR) * DIFF4 * CDZ(k)**nd_order4 &
                                         * 0.25_RP * ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) )
          enddo

          num_diff( 1:KS-1,i,j,I_MOMY,ZDIR) = 0.0_RP
          num_diff(KE:KA  ,i,j,I_MOMY,ZDIR) = 0.0_RP

          do k = KS, KE
             num_diff(k,i,j,I_MOMY,XDIR) = num_diff_pt0(k,i,j,I_MOMY,XDIR) * DIFF4 * CDX(i)**nd_order4 &
                                         * 0.25_RP * ( DENS(k,i+1,j+1)+DENS(k,i,j+1)+DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          num_diff(KS  ,i,j,I_MOMY,XDIR) = num_DIFF(KS  ,i,j,I_MOMY,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMY,XDIR) = num_DIFF(KS+1,i,j,I_MOMY,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_MOMY,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_MOMY,XDIR) = 0.0_RP

          do k = KS, KE
             num_diff(k,i,j,I_MOMY,YDIR) = num_diff_pt0(k,i,j,I_MOMY,YDIR) * DIFF4 * FDY(j)**nd_order4 &
                                         * DENS(k,i,j)
          enddo
          num_diff(KS  ,i,j,I_MOMY,YDIR) = num_DIFF(KS  ,i,j,I_MOMY,YDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMY,YDIR) = num_DIFF(KS+1,i,j,I_MOMY,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_MOMY,YDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_MOMY,YDIR) = 0.0_RP
       enddo
       enddo

       !-----< rho * theta >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS  , KE-1
             num_diff(k,i,j,I_RHOT,ZDIR) = num_diff_pt0(k,i,j,I_RHOT,ZDIR) * DIFF4 * CDZ(k)**nd_order4 &
                                         * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
          enddo
          num_diff(KS-1,i,j,I_RHOT,ZDIR) = num_diff_pt0(KS-1,i,j,I_RHOT,ZDIR) * DIFF4 * CDZ(KS-1)**nd_order4 &
                                         * DENS(KS,i,j)
          num_diff(KE  ,i,j,I_RHOT,ZDIR) = num_diff_pt0(KE  ,i,j,I_RHOT,ZDIR) * DIFF4 * CDZ(KE  )**nd_order4 &
                                         * DENS(KE,i,j)

          num_diff(   1:KS-2,i,j,I_RHOT,ZDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_RHOT,ZDIR) = 0.0_RP

          do k = KS, KE
             num_diff(k,i,j,I_RHOT,XDIR) = num_diff_pt0(k,i,j,I_RHOT,XDIR) * DIFF4 * CDX(i)**nd_order4 &
                                         * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          num_diff(KS  ,i,j,I_RHOT,XDIR) = num_DIFF(KS  ,i,j,I_RHOT,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_RHOT,XDIR) = num_DIFF(KS+1,i,j,I_RHOT,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_RHOT,XDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_RHOT,XDIR) = 0.0_RP

          do k = KS, KE
             num_diff(k,i,j,I_RHOT,YDIR) = num_diff_pt0(k,i,j,I_RHOT,YDIR) * DIFF4 * CDY(j)**nd_order4 &
                                         * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          num_diff(KS  ,i,j,I_RHOT,YDIR) = num_DIFF(KS  ,i,j,I_RHOT,YDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_RHOT,YDIR) = num_DIFF(KS+1,i,j,I_RHOT,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          num_diff(   1:KS-1,i,j,I_RHOT,YDIR) = 0.0_RP
          num_diff(KE+1:KA  ,i,j,I_RHOT,YDIR) = 0.0_RP
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

    return
  end subroutine ATMOS_DYN_numfilter_coef

  !-----------------------------------------------------------------------------
  !> Calc coefficient of numerical filter
  subroutine ATMOS_DYN_numfilter_coef_q( &
       num_diff_q,                             &
       DENS, QTRC,                             &
       CNZ3, CNX3, CNY3, CNZ4, CNX4, CNY4,     &
       CDZ, CDX, CDY,                          &
       REF_qv, iq,                             &
       DIFF4, ND_ORDER, ND_SFC_FACT, ND_USE_RS )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: num_diff_q(KA,IA,JA,3)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA)

    real(RP), intent(in)  :: CNZ3(3,KA,2)
    real(RP), intent(in)  :: CNX3(3,IA,2)
    real(RP), intent(in)  :: CNY3(3,JA,2)
    real(RP), intent(in)  :: CNZ4(5,KA,2)
    real(RP), intent(in)  :: CNX4(5,IA,2)
    real(RP), intent(in)  :: CNY4(5,JA,2)

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)

    real(RP), intent(in)  :: REF_qv(KA,IA,JA)
    integer,  intent(in)  :: iq

    real(RP), intent(in)  :: DIFF4
    integer,  intent(in)  :: ND_ORDER
    real(RP), intent(in)  :: ND_SFC_FACT
    logical,  intent(in)  :: ND_USE_RS

    real(RP) :: qv_diff(KA,IA,JA) ! anomary of water vapor

    real(RP), target  :: num_diff     (KA,IA,JA,5,3)
    real(RP), target  :: num_diff_work(KA,IA,JA,5,3)
    real(RP), pointer :: num_diff_pt0 (:,:,:,:,:)
    real(RP), pointer :: num_diff_pt1 (:,:,:,:,:)
    real(RP), pointer :: tmp_pt       (:,:,:,:,:)

    integer  :: nd_order4, no

    integer :: IIS, IIE
    integer :: JJS, JJE
    integer :: k, i, j
    !---------------------------------------------------------------------------

    !###########################################################################
    ! 1st order coefficients
    !###########################################################################

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       if ( iq == I_QV ) then

          if ( ND_USE_RS ) then
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                qv_diff(k,i,j) = QTRC(k,i,j) - REF_qv(k,i,j)
             enddo
             enddo
             enddo
          else
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
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
             do j  = JJS, JJE
             do i  = IIS, IIE
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
             enddo
             enddo
          endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, CNZ3(1,k+1,1) )
             call CHECK( __LINE__, CNZ3(2,k+1,1) )
             call CHECK( __LINE__, CNZ3(3,k+1,1) )
             call CHECK( __LINE__, CNZ3(1,k,1) )
             call CHECK( __LINE__, qv_diff(k+2,i,j) )
             call CHECK( __LINE__, qv_diff(k+1,i,j) )
             call CHECK( __LINE__, qv_diff(k  ,i,j) )
             call CHECK( __LINE__, qv_diff(k-1,i,j) )
#endif
             num_diff(k,i,j,1,ZDIR) = ( + CNZ3(1,k+1,1) * qv_diff(k+2,i,j) &
                                        - CNZ3(2,k+1,1) * qv_diff(k+1,i,j) &
                                        + CNZ3(3,k+1,1) * qv_diff(k  ,i,j) &
                                        - CNZ3(1,k  ,1) * qv_diff(k-1,i,j) )
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
             call CHECK( __LINE__, qv_diff(KS+2,i,j) )
             call CHECK( __LINE__, qv_diff(KS+1,i,j) )
             call CHECK( __LINE__, qv_diff(KS,i,j) )
#endif
             num_diff(KS  ,i,j,1,ZDIR) = ( + CNZ3(1,KS+1,1) * qv_diff(KS+2,i,j) &
                                           - CNZ3(2,KS+1,1) * qv_diff(KS+1,i,j) &
                                           + CNZ3(3,KS+1,1) * qv_diff(KS  ,i,j) &
                                           - CNZ3(1,KS  ,1) * qv_diff(KS+1,i,j) )
             num_diff(KS-1,i,j,1,ZDIR) = - num_diff(KS  ,i,j,1,ZDIR)
             num_diff(KS-2,i,j,1,ZDIR) = - num_diff(KS+1,i,j,1,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
             call CHECK( __LINE__, CNZ3(1,KE,1) )
             call CHECK( __LINE__, CNZ3(2,KE,1) )
             call CHECK( __LINE__, CNZ3(3,KE,1) )
             call CHECK( __LINE__, CNZ3(1,KE-1,1) )
             call CHECK( __LINE__, qv_diff(KE,i,j) )
             call CHECK( __LINE__, qv_diff(KE-1,i,j) )
             call CHECK( __LINE__, qv_diff(KE-2,i,j) )
#endif
             num_diff(KE-1,i,j,1,ZDIR) = ( + CNZ3(1,KE  ,1) * qv_diff(KE-1,i,j) &
                                           - CNZ3(2,KE  ,1) * qv_diff(KE  ,i,j) &
                                           + CNZ3(3,KE  ,1) * qv_diff(KE-1,i,j) &
                                           - CNZ3(1,KE-1,1) * qv_diff(KE-2,i,j) )
             num_diff(KE  ,i,j,1,ZDIR) = - num_diff(KE-1,i,j,1,ZDIR)
             num_diff(KE+1,i,j,1,ZDIR) = - num_diff(KE-2,i,j,1,ZDIR)
          enddo
          enddo

       else ! iq /= I_QV

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-2
#ifdef DEBUG
             call CHECK( __LINE__, CNZ3(1,k+1,1) )
             call CHECK( __LINE__, CNZ3(2,k+1,1) )
             call CHECK( __LINE__, CNZ3(3,k+1,1) )
             call CHECK( __LINE__, CNZ3(1,k,1) )
             call CHECK( __LINE__, QTRC(k+2,i,j) )
             call CHECK( __LINE__, QTRC(k+1,i,j) )
             call CHECK( __LINE__, QTRC(k  ,i,j) )
             call CHECK( __LINE__, QTRC(k-1,i,j) )
#endif
             num_diff(k,i,j,1,ZDIR) = ( + CNZ3(1,k+1,1) * QTRC(k+2,i,j) &
                                        - CNZ3(2,k+1,1) * QTRC(k+1,i,j) &
                                        + CNZ3(3,k+1,1) * QTRC(k  ,i,j) &
                                        - CNZ3(1,k  ,1) * QTRC(k-1,i,j) )
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
             call CHECK( __LINE__, QTRC(KS+2,i,j) )
             call CHECK( __LINE__, QTRC(KS+1,i,j) )
             call CHECK( __LINE__, QTRC(KS,i,j) )
#endif
             num_diff(KS  ,i,j,1,ZDIR) = ( + CNZ3(1,KS+1,1) * QTRC(KS+2,i,j) &
                                           - CNZ3(2,KS+1,1) * QTRC(KS+1,i,j) &
                                           + CNZ3(3,KS+1,1) * QTRC(KS  ,i,j) &
                                           - CNZ3(1,KS  ,1) * QTRC(KS+1,i,j) )
             num_diff(KS-1,i,j,1,ZDIR) = - num_diff(KS  ,i,j,1,ZDIR)
             num_diff(KS-2,i,j,1,ZDIR) = - num_diff(KS+1,i,j,1,ZDIR)
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
#ifdef DEBUG
             call CHECK( __LINE__, CNZ3(1,KE,1) )
             call CHECK( __LINE__, CNZ3(2,KE,1) )
             call CHECK( __LINE__, CNZ3(3,KE,1) )
             call CHECK( __LINE__, CNZ3(1,KE-1,1) )
             call CHECK( __LINE__, QTRC(KE,i,j) )
             call CHECK( __LINE__, QTRC(KE-1,i,j) )
             call CHECK( __LINE__, QTRC(KE-2,i,j) )
#endif
             num_diff(KE-1,i,j,1,ZDIR) = ( + CNZ3(1,KE  ,1) * QTRC(KE-1,i,j) &
                                           - CNZ3(2,KE  ,1) * QTRC(KE  ,i,j) &
                                           + CNZ3(3,KE  ,1) * QTRC(KE-1,i,j) &
                                           - CNZ3(1,KE-1,1) * QTRC(KE-2,i,j) )
             num_diff(KE  ,i,j,1,ZDIR) = - num_diff(KE-1,i,j,1,ZDIR)
             num_diff(KE+1,i,j,1,ZDIR) = - num_diff(KE-2,i,j,1,ZDIR)
          enddo
          enddo

       endif ! QV or not?

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, CNX3(1,i+1,1) )
          call CHECK( __LINE__, CNX3(2,i+1,1) )
          call CHECK( __LINE__, CNX3(3,i+1,1) )
          call CHECK( __LINE__, CNX3(1,i  ,1) )
          call CHECK( __LINE__, QTRC(k,i+2,j) )
          call CHECK( __LINE__, QTRC(k,i+1,j) )
          call CHECK( __LINE__, QTRC(k,i  ,j) )
          call CHECK( __LINE__, QTRC(k,i-1,j) )
#endif
          num_diff(k,i,j,1,XDIR) = ( + CNX3(1,i+1,1) * QTRC(k,i+2,j) &
                                     - CNX3(2,i+1,1) * QTRC(k,i+1,j) &
                                     + CNX3(3,i+1,1) * QTRC(k,i  ,j) &
                                     - CNX3(1,i  ,1) * QTRC(k,i-1,j) )
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
          call CHECK( __LINE__, QTRC(k,i,j+2) )
          call CHECK( __LINE__, QTRC(k,i,j+1) )
          call CHECK( __LINE__, QTRC(k,i,j  ) )
          call CHECK( __LINE__, QTRC(k,i,j-1) )
#endif
          num_diff(k,i,j,1,YDIR) = ( + CNY3(1,j+1,1) * QTRC(k,i,j+2) &
                                     - CNY3(2,j+1,1) * QTRC(k,i,j+1) &
                                     + CNY3(3,j+1,1) * QTRC(k,i,j  ) &
                                     - CNY3(1,j  ,1) * QTRC(k,i,j-1) )
       enddo
       enddo
       enddo

    enddo
    enddo ! end tile

    !###########################################################################
    ! High order coefficients
    !###########################################################################

    num_diff_pt0 => num_diff
    num_diff_pt1 => num_diff_work

    do no = 2, nd_order

       call COMM_vars8( num_diff_pt0(:,:,:,1,ZDIR), 1 )
       call COMM_vars8( num_diff_pt0(:,:,:,1,XDIR), 2 )
       call COMM_vars8( num_diff_pt0(:,:,:,1,YDIR), 3 )

       call COMM_wait ( num_diff_pt0(:,:,:,1,ZDIR), 1 )
       call COMM_wait ( num_diff_pt0(:,:,:,1,XDIR), 2 )
       call COMM_wait ( num_diff_pt0(:,:,:,1,YDIR), 3 )

       call calc_numdiff4( num_diff_pt1(:,:,:,:,:), & ! (out)
                           num_diff_pt0(:,:,:,:,:), & ! (in)
                           CNZ4(:,:,1),             & ! (in)
                           CNX4(:,:,1),             & ! (in)
                           CNY4(:,:,1),             & ! (in)
                           1,                       & ! (in)
                           KE                       ) ! (in)

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
          do k = KS  , KE-1
             num_diff_q(k,i,j,ZDIR) = num_diff_pt0(k,i,j,1,ZDIR) * DIFF4 * CDZ(k)**nd_order4 &
                                    * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
          enddo
          num_diff_q(KS-1,i,j,ZDIR) = num_diff_pt0(KS-1,i,j,1,ZDIR) * DIFF4 * CDZ(KS-1)**nd_order4 &
                                    * DENS(KS,i,j)
          num_diff_q(KE  ,i,j,ZDIR) = num_diff_pt0(KE  ,i,j,1,ZDIR) * DIFF4 * CDZ(KE  )**nd_order4 &
                                    * DENS(KE,i,j)

          do k = KS, KE
             num_diff_q(k,i,j,XDIR) = num_diff_pt0(k,i,j,1,XDIR) * DIFF4 * CDX(i)**nd_order4 &
                                    * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          num_diff_q(KS  ,i,j,XDIR) = num_DIFF(KS  ,i,j,1,XDIR) * ND_SFC_FACT
          num_diff_q(KS+1,i,j,XDIR) = num_DIFF(KS+1,i,j,1,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

          do k = KS, KE
             num_diff_q(k,i,j,YDIR) = num_diff_pt0(k,i,j,1,YDIR) * DIFF4 * CDY(j)**nd_order4 &
                                    * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          num_diff_q(KS  ,i,j,YDIR) = num_DIFF(KS  ,i,j,1,YDIR) * ND_SFC_FACT
          num_diff_q(KS+1,i,j,YDIR) = num_DIFF(KS+1,i,j,1,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( num_diff_q(:,:,:,ZDIR), 1 )
    call COMM_vars8( num_diff_q(:,:,:,XDIR), 2 )
    call COMM_vars8( num_diff_q(:,:,:,YDIR), 3 )

    call COMM_wait ( num_diff_q(:,:,:,ZDIR), 1 )
    call COMM_wait ( num_diff_q(:,:,:,XDIR), 2 )
    call COMM_wait ( num_diff_q(:,:,:,YDIR), 3 )

    return
  end subroutine ATMOS_DYN_numfilter_coef_q

  !-----------------------------------------------------------------------------
  subroutine calc_numdiff4( &
       num_diff_pt1, &
       num_diff_pt0, &
       CNZ4,         &
       CNX4,         &
       CNY4,         &
       I_val,        &
       k1            )
    implicit none

    real(RP), intent(out) :: num_diff_pt1(KA,IA,JA,5,3)
    real(RP), intent(in)  :: num_diff_pt0(KA,IA,JA,5,3)
    real(RP), intent(in)  :: CNZ4(5,KA)
    real(RP), intent(in)  :: CNX4(5,IA)
    real(RP), intent(in)  :: CNY4(5,JA)
    integer,  intent(in)  :: I_val
    integer,  intent(in)  :: k1

    integer :: i, j, k
    integer :: IIS, IIE, JJS, JJE
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

  !-----------------------------------------------------------------------------
  !> Flux Correction Transport Limiter
  subroutine ATMOS_DYN_fct(      &
       qflx_anti,                &
       phi_in, qflx_hi, qflx_lo, &
       rdz, rdx, rdy, dtrk       )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       IUNDEF => CONST_UNDEF2, &
       EPSILON => CONST_EPS
    use scale_comm, only: &
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

    ! work for FCT
    real(RP) :: phi_lo(KA,IA,JA)
    real(RP) :: pjpls(KA,IA,JA)
    real(RP) :: pjmns(KA,IA,JA)
    real(RP) :: qjpls(KA,IA,JA)
    real(RP) :: qjmns(KA,IA,JA)
    real(RP) :: rjpls(KA,IA,JA)
    real(RP) :: rjmns(KA,IA,JA)

    real(RP) :: zerosw, dirsw

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
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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
                 * ( min( rjpls(k+1,i,j),rjmns(k  ,i,j) ) * (          dirsw ) &
                   + min( rjpls(k  ,i,j),rjmns(k+1,i,j) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
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
  end subroutine ATMOS_DYN_fct

end module scale_atmos_dyn_common
