!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          numerical filter with FVM for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_fvm_numfilter
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
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
  public :: ATMOS_DYN_fvm_numfilter_setup
  public :: ATMOS_DYN_fvm_numfilter_flux
  public :: ATMOS_DYN_fvm_numfilter_flux_q
  public :: ATMOS_DYN_fvm_numfilter_tend
  public :: ATMOS_DYN_fvm_apply_numfilter

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  private :: calc_diff1
  private :: calc_diff3
  private :: calc_diff4
  private :: calc_numdiff

  private :: fvm_add_FluxDiv_xyz
  private :: fvm_add_FluxDiv_xyw
  private :: fvm_add_FluxDiv_uyz
  private :: fvm_add_FluxDiv_xvz

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  real(RP), allocatable :: CNZ1(:,:)
  real(RP), allocatable :: CNX1(:,:)
  real(RP), allocatable :: CNY1(:,:)
  real(RP), allocatable :: CNZ3(:,:,:)
  real(RP), allocatable :: CNX3(:,:,:)
  real(RP), allocatable :: CNY3(:,:,:)
  real(RP), allocatable :: CNZ4(:,:,:)
  real(RP), allocatable :: CNX4(:,:,:)
  real(RP), allocatable :: CNY4(:,:,:)

  integer :: I_COMM_DENS_Z = 1
  integer :: I_COMM_DENS_X = 2
  integer :: I_COMM_DENS_Y = 3
  integer :: I_COMM_MOMZ_Z = 4
  integer :: I_COMM_MOMZ_X = 5
  integer :: I_COMM_MOMZ_Y = 6
  integer :: I_COMM_MOMX_Z = 7
  integer :: I_COMM_MOMX_X = 8
  integer :: I_COMM_MOMX_Y = 9
  integer :: I_COMM_MOMY_Z = 10
  integer :: I_COMM_MOMY_X = 11
  integer :: I_COMM_MOMY_Y = 12
  integer :: I_COMM_RHOT_Z = 13
  integer :: I_COMM_RHOT_X = 14
  integer :: I_COMM_RHOT_Y = 15
  integer :: I_COMM_QTRC_Z = 1
  integer :: I_COMM_QTRC_X = 2
  integer :: I_COMM_QTRC_Y = 3

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_fvm_numfilter_setup( &
       num_diff, num_diff_q,        &
       CDZ, CDX, CDY, FDZ, FDX, FDY )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    implicit none
    real(RP), intent(inout) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(inout) :: num_diff_q(KA,IA,JA,3)
    real(RP), intent(in)    :: CDZ(KA)
    real(RP), intent(in)    :: CDX(IA)
    real(RP), intent(in)    :: CDY(JA)
    real(RP), intent(in)    :: FDZ(KA-1)
    real(RP), intent(in)    :: FDX(IA-1)
    real(RP), intent(in)    :: FDY(JA-1)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( ( ( .not. PRC_TwoD ) .and. IHALO < 2 ) .or. JHALO < 2 .or. KHALO < 2 ) then
       LOG_ERROR("ATMOS_DYN_filter_setup",*) 'number of HALO must be at least 2 for numrical filter'
       call PRC_abort
    end if

    ! allocation
    allocate( CNZ1(KA,2) )
    allocate( CNX1(IA,2) )
    allocate( CNY1(JA,2) )    
    allocate( CNZ3(3,KA,2) )
    allocate( CNX3(3,IA,2) )
    allocate( CNY3(3,JA,2) )
    allocate( CNZ4(5,KA,2) )
    allocate( CNX4(5,IA,2) )
    allocate( CNY4(5,JA,2) )


    call COMM_vars8_init( 'num_diff_DENS_Z', num_diff(:,:,:,I_DENS,ZDIR), I_COMM_DENS_Z )
    call COMM_vars8_init( 'num_diff_DENS_X', num_diff(:,:,:,I_DENS,XDIR), I_COMM_DENS_X )
    call COMM_vars8_init( 'num_diff_DENS_Y', num_diff(:,:,:,I_DENS,YDIR), I_COMM_DENS_Y )
    call COMM_vars8_init( 'num_diff_MOMZ_Z', num_diff(:,:,:,I_MOMZ,ZDIR), I_COMM_MOMZ_Z )
    call COMM_vars8_init( 'num_diff_MOMZ_X', num_diff(:,:,:,I_MOMZ,XDIR), I_COMM_MOMZ_X )
    call COMM_vars8_init( 'num_diff_MOMZ_Y', num_diff(:,:,:,I_MOMZ,YDIR), I_COMM_MOMZ_Y )
    call COMM_vars8_init( 'num_diff_MOMX_Z', num_diff(:,:,:,I_MOMX,ZDIR), I_COMM_MOMX_Z )
    call COMM_vars8_init( 'num_diff_MOMX_X', num_diff(:,:,:,I_MOMX,XDIR), I_COMM_MOMX_X )
    call COMM_vars8_init( 'num_diff_MOMX_Y', num_diff(:,:,:,I_MOMX,YDIR), I_COMM_MOMX_Y )
    call COMM_vars8_init( 'num_diff_MOMY_Z', num_diff(:,:,:,I_MOMY,ZDIR), I_COMM_MOMY_Z )
    call COMM_vars8_init( 'num_diff_MOMY_X', num_diff(:,:,:,I_MOMY,XDIR), I_COMM_MOMY_X )
    call COMM_vars8_init( 'num_diff_MOMY_Y', num_diff(:,:,:,I_MOMY,YDIR), I_COMM_MOMY_Y )
    call COMM_vars8_init( 'num_diff_RHOT_Z', num_diff(:,:,:,I_RHOT,ZDIR), I_COMM_RHOT_Z )
    call COMM_vars8_init( 'num_diff_RHOT_X', num_diff(:,:,:,I_RHOT,XDIR), I_COMM_RHOT_X )
    call COMM_vars8_init( 'num_diff_RHOT_Y', num_diff(:,:,:,I_RHOT,YDIR), I_COMM_RHOT_Y )

    call COMM_vars8_init( 'num_diff_QTRC_Z', num_diff_q(:,:,:,ZDIR), I_COMM_QTRC_Z )
    call COMM_vars8_init( 'num_diff_QTRC_X', num_diff_q(:,:,:,XDIR), I_COMM_QTRC_X )
    call COMM_vars8_init( 'num_diff_QTRC_Y', num_diff_q(:,:,:,YDIR), I_COMM_QTRC_Y )

#ifdef DEBUG
    CNX1(:,:)   = UNDEF
    CNY1(:,:)   = UNDEF
    CNZ1(:,:)   = UNDEF
    CNZ3(:,:,:) = UNDEF
    CNX3(:,:,:) = UNDEF
    CNY3(:,:,:) = UNDEF
    CNZ4(:,:,:) = UNDEF
    CNX4(:,:,:) = UNDEF
    CNY4(:,:,:) = UNDEF
#endif

    !* z direction **********************************************

    do k = KS-1, KE
      CNZ1(k,1) = 1.0_RP / FDZ(k)
    end do

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

    !-

    do k = KS, KE
      CNZ1(k,2) = 1.0_RP / CDZ(k)
    end do

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
!    CNZ4(1,KE,2) = ( CNZ3(1,KE+1,2)                ) / FDZ(KE-1)
    CNZ4(2,KE,2) = ( CNZ3(2,KE+1,2) + CNZ3(1,KE,2) ) / FDZ(KE-1)
    CNZ4(3,KE,2) = ( CNZ3(3,KE+1,2) + CNZ3(2,KE,2) ) / FDZ(KE-1)
    CNZ4(4,KE,2) = ( CNZ3(1,KE  ,2) + CNZ3(3,KE,2) ) / FDZ(KE-1)

    !* x direction *************************************************
    if ( .not. PRC_TwoD ) then

       do i = IS-1, IE
          CNX1(i,1) = 1.0_RP / FDX(i)
       end do

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

       !-

       do i = IS, IE
          CNX1(i,2) = 1.0_RP / CDX(i)
       end do

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
    end if

    !* y direction ************************************************    

    do j = JS-1, JE
      CNY1(j,1) = 1.0_RP / FDY(j)
    end do

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

    !-

    do j = JS, JE
      CNY1(j,2) = 1.0_RP / CDY(j)
    end do

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

    return
  end subroutine ATMOS_DYN_fvm_numfilter_setup

  !-----------------------------------------------------------------------------
  !> Calculate fluxes with numerical filter for prognostic variables of dynamical core
  !
  ! The num_diff (for example, in the x direction and for the variable, RHOT) 
  ! is calculated in the following manner: 
  !   num_diff = - (-1)^(n+1) * (ND_COEF * (Delta x)^2n / (2^2n Delta t))
  !              * DENS * (d^2n-1 THETA / dx^2n-1)
  ! where n is ND_LAPLACIAN_NUM. 
  ! Note that we mulitply a minus sign in the above expression since num_diff fluxes calculated here are 
  ! directly added to advective fluxes. 
  !
  subroutine ATMOS_DYN_fvm_numfilter_flux( &
       num_diff,                                         &
       DENS, MOMZ, MOMX, MOMY, RHOT,                     &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                     &
       TwoD, DT,                                         &
       REF_dens, REF_pott,                               &
       ND_COEF, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS )
    use scale_comm_cartesC, only: &
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

    logical,  intent(in)  :: TwoD
    real(RP), intent(in)  :: DT

    real(RP), intent(in)  :: REF_dens(KA,IA,JA)
    real(RP), intent(in)  :: REF_pott(KA,IA,JA)

    real(RP), intent(in)  :: ND_COEF
    integer,  intent(in)  :: ND_LAPLACIAN_NUM
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

    real(RP) :: diff_coef_tmp
    integer  :: nd_order
    real(RP) :: nd_coef_cdz(KA)
    real(RP) :: nd_coef_cdx(IA)
    real(RP) :: nd_coef_cdy(JA)
    real(RP) :: nd_coef_fdz(KA-1)
    real(RP) :: nd_coef_fdx(IA-1)
    real(RP) :: nd_coef_fdy(JA-1)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! calculate the coefficient numerical diffusion

    nd_order = ND_LAPLACIAN_NUM * 2
    diff_coef_tmp = - (-1)**(mod(ND_LAPLACIAN_NUM+1,2)) &
                    * ND_COEF / ( 2**(nd_order) * DT )
    do k = KS-1, KE
       nd_coef_cdz(k) = diff_coef_tmp * CDZ(k)**nd_order
    end do
    do k = KS+1, KE-1
       nd_coef_fdz(k) = diff_coef_tmp * FDZ(k)**nd_order
    end do
    if ( .not. TwoD ) then
       do i = IS, IE
          nd_coef_cdx(i) = diff_coef_tmp * CDX(i)**nd_order
          nd_coef_fdx(i) = diff_coef_tmp * FDX(i)**nd_order
       end do
    end if
    do j = JS, JE
       nd_coef_cdy(j) = diff_coef_tmp * CDY(j)**nd_order
       nd_coef_fdy(j) = diff_coef_tmp * FDY(j)**nd_order
    end do

    if ( .NOT. ND_USE_RS ) then

       ! In order to relax the diffusion of the reference state, 
       ! we smooth the density and potential temperature fields 
       ! given for calculating fluxes of numerical diffusion. 

       call PROF_rapstart("NumFilter_Main", 3)

       do j = JS-2, JE+2
       do i = max(IS-2,1), min(IE+2,IA)
       do k = KS, KE
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       end do
       end do
       end do

       if ( TwoD ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = JS, JE
          do k = KS+1, KE-1
             dens_diff(k,IS,j) = ( ( DENS(k,IS,j)                  ) * 3.0_RP &
                                 + ( DENS(k,IS,j+1)+DENS(k,IS,j-1) ) * 2.0_RP &
                                 + ( DENS(k,IS,j+2)+DENS(k,IS,j-2) ) &
                                 + ( DENS(k+1,IS,j)+DENS(k-1,IS,j) ) * 2.0_RP &
                                 ) / 13.0_RP

             pott_diff(k,IS,j) = ( ( POTT(k,IS,j)                  ) * 3.0_RP &
                                 + ( POTT(k,IS,j+1)+POTT(k,IS,j-1) ) * 2.0_RP &
                                 + ( POTT(k,IS,j+2)+POTT(k,IS,j-2) ) &
                                 + ( POTT(k+1,IS,j)+POTT(k-1,IS,j) ) * 2.0_RP &
                                 ) / 13.0_RP
          enddo
          enddo

          !$omp parallel do
          do j  = JS, JE
             dens_diff(KS,IS,j) = ( ( DENS(KS,IS,j)                   ) * 3.0_RP &
                                  + ( DENS(KS,IS,j+1)+DENS(KS,IS,j-1) ) * 2.0_RP &
                                  + ( DENS(KS,IS,j+2)+DENS(KS,IS,j-2) ) &
                                  + ( DENS(KS+1,IS,j)                 ) * 2.0_RP &
                                  ) / 11.0_RP
             dens_diff(KE,IS,j) = ( ( DENS(KE,IS,j)                   ) * 3.0_RP &
                                  + ( DENS(KE,IS,j+1)+DENS(KE,IS,j-1) ) * 2.0_RP &
                                  + ( DENS(KE,IS,j+2)+DENS(KE,IS,j-2) ) &
                                  + ( DENS(KE-1,IS,j)                 ) * 2.0_RP &
                                  ) / 11.0_RP

             pott_diff(KS,IS,j) = ( ( POTT(KS,IS,j)                   ) * 3.0_RP &
                                  + ( POTT(KS,IS,j+1)+POTT(KS,IS,j-1) ) * 2.0_RP &
                                  + ( POTT(KS,IS,j+2)+POTT(KS,IS,j-2) ) &
                                  + ( POTT(KS+1,IS,j)                 ) * 2.0_RP &
                                  ) / 11.0_RP
             pott_diff(KE,IS,j) = ( ( POTT(KE,IS,j)                  ) * 3.0_RP &
                                 + ( POTT(KE,IS,j+1)+POTT(KE,IS,j-1) ) * 2.0_RP &
                                 + ( POTT(KE,IS,j+2)+POTT(KE,IS,j-2) ) &
                                 + ( POTT(KE-1,IS,j)                 ) * 2.0_RP &
                                 ) / 11.0_RP
          end do
       else
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

          !$omp parallel do
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
       end if

       call PROF_rapend  ("NumFilter_Main", 3)

       call PROF_rapstart("NumFilter_Comm", 3)

       call COMM_vars8( dens_diff,  1 )
       call COMM_vars8( pott_diff,  2 )

       call COMM_wait ( dens_diff,  1 )
       call COMM_wait ( pott_diff,  2 )

       call PROF_rapend  ("NumFilter_Comm", 3)

    end if


    !-----< density >-----

    if ( ND_USE_RS ) then

       call PROF_rapstart("NumFilter_Main", 3)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = max(IS-1,1), min(IE+2,IA)
       do k = KS, KE
          dens_diff(k,i,j) = DENS(k,i,j) - REF_dens(k,i,j)
       enddo
       enddo
       enddo

       call PROF_rapend("NumFilter_Main", 3)

    endif

    call calc_numdiff( work, iwork,      & ! (out)
                       dens_diff,        & ! (in)
                       TwoD,             & ! (in)
                       ND_LAPLACIAN_NUM, & ! (in)
                       0, 0, 0, KE )

    call PROF_rapstart("NumFilter_Main", 3)

    !-----< density >-----

    !$omp parallel private(i,j,k) 
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       num_diff(k,i,j,I_DENS,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k)
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-2,i,j,I_DENS,ZDIR) = 0.0_RP
       num_diff(KE+1:KA  ,i,j,I_DENS,ZDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait
    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          num_diff(k,i,j,I_DENS,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          num_diff(   1:KS-1,i,j,I_DENS,XDIR) = 0.0_RP
          num_diff(KS  ,i,j,I_DENS,XDIR) = num_diff(KS  ,i,j,I_DENS,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_DENS,XDIR) = num_diff(KS+1,i,j,I_DENS,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          num_diff(KE+1:KA  ,i,j,I_DENS,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait
    end if
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_DENS,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j)
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_DENS,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_DENS,YDIR) = num_diff(KS  ,i,j,I_DENS,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_DENS,YDIR) = num_diff(KS+1,i,j,I_DENS,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_DENS,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do
    !$omp end parallel
    call PROF_rapend  ("NumFilter_Main", 3)

    call PROF_rapstart("NumFilter_Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_DENS,ZDIR),  I_COMM_DENS_Z )
    if ( .not. TwoD ) &
    call COMM_vars8( num_diff(:,:,:,I_DENS,XDIR),  I_COMM_DENS_X )
    call COMM_vars8( num_diff(:,:,:,I_DENS,YDIR),  I_COMM_DENS_Y )

    call PROF_rapend  ("NumFilter_Comm", 3)


    !-----< z-momentum >-----

    call PROF_rapstart("NumFilter_Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+2
    do i = max(IS-2,1), min(IE+2,IA)
    do k = KS, KE-1
       VELZ(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo

    call PROF_rapend  ("NumFilter_Main", 3)

    call calc_numdiff( work, iwork,      & ! (out)
                       VELZ,             & ! (in)
                       TwoD,             & ! (in)
                       ND_LAPLACIAN_NUM, & ! (in)
                       1, 0, 0, KE-1 )

    call PROF_rapstart("NumFilter_Main", 3)

    !$omp parallel private(i,j,k)

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE-1
       num_diff(k,i,j,I_MOMZ,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_fdz(k) &
                                   * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS,i,j,I_MOMZ,ZDIR) = 0.0_RP
       num_diff(KE:KA,i,j,I_MOMZ,ZDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          num_diff(k,i,j,I_MOMZ,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                      * 0.25_RP * ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)    
       do j = JS, JE
       do i = IS, IE
          num_diff( 1:KS-1,i,j,I_MOMZ,XDIR) = 0.0_RP
          num_diff(KS  ,i,j,I_MOMZ,XDIR) = num_diff(KS  ,i,j,I_MOMZ,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMZ,XDIR) = num_diff(KS+1,i,j,I_MOMZ,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          num_diff(KE:KA  ,i,j,I_MOMZ,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait
    end if

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMZ,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                                   * 0.25_RP * ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)    
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMZ,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMZ,YDIR) = num_diff(KS  ,i,j,I_MOMZ,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMZ,YDIR) = num_diff(KS+1,i,j,I_MOMZ,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP

       num_diff(KE:KA  ,i,j,I_MOMZ,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do

    !$omp end parallel
    call PROF_rapend  ("NumFilter_Main", 3)

    call PROF_rapstart("NumFilter_Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_MOMZ,ZDIR), I_COMM_MOMZ_Z )
    if ( .not. TwoD ) &
    call COMM_vars8( num_diff(:,:,:,I_MOMZ,XDIR), I_COMM_MOMZ_X )
    call COMM_vars8( num_diff(:,:,:,I_MOMZ,YDIR), I_COMM_MOMZ_Y )

    call PROF_rapend  ("NumFilter_Comm", 3)


    !-----< x-momentum >-----

    call PROF_rapstart("NumFilter_Main", 3)

    if ( TwoD ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_
       do j = JS-2, JE+2
       do k = KS, KE
          VELX(k,IS,j) = MOMX(k,IS,j) / DENS(k,IS,j)
       enddo
       enddo
    else
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-2, JE+2
       do i = IS-2, IE+1
       do k = KS, KE
          VELX(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
    end if

    call PROF_rapend  ("NumFilter_Main", 3)

    call calc_numdiff( work, iwork,      & ! (out)
                       VELX,             & ! (in)
                       TwoD,             & ! (in)
                       ND_LAPLACIAN_NUM, & ! (in)
                       0, 1, 0, KE )


    call PROF_rapstart("NumFilter_Main", 3)

    !$omp parallel private(i,j,k) 

    if ( TwoD ) then
       !$omp do OMP_SCHEDULE_
       do j = JS, JE
       do k = KS, KE-1
          num_diff(k,IS,j,I_MOMX,ZDIR) = work(k,IS,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                       * 0.5_RP * ( DENS(k+1,IS,j)+DENS(k,IS,j) )
       enddo
       enddo
       !$omp end do nowait
    else
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE-1
          num_diff(k,i,j,I_MOMX,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                      * 0.25_RP * ( DENS(k+1,i+1,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
    end if
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMX,ZDIR) = 0.0_RP
       num_diff(KE:KA  ,i,j,I_MOMX,ZDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          num_diff(k,i,j,I_MOMX,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_fdx(i) &
                                      * DENS(k,i,j)
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          num_diff(   1:KS-1,i,j,I_MOMX,XDIR) = 0.0_RP
          num_diff(KS  ,i,j,I_MOMX,XDIR) = num_diff(KS  ,i,j,I_MOMX,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMX,XDIR) = num_diff(KS+1,i,j,I_MOMX,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          num_diff(KE+1:KA  ,i,j,I_MOMX,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait
    end if

    if ( TwoD ) then
       !$omp do OMP_SCHEDULE_
       do j = JS, JE
       do k = KS, KE
          num_diff(k,IS,j,I_MOMX,YDIR) = work(k,IS,j,YDIR,iwork) * nd_coef_cdy(j) &
                                       * 0.5_RP * ( DENS(k,IS,j+1)+DENS(k,IS,j) )
       enddo
       enddo
       !$omp end do nowait
    else
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          num_diff(k,i,j,I_MOMX,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                                      * 0.25_RP * ( DENS(k,i+1,j+1)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
    end if
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_MOMX,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMX,YDIR) = num_diff(KS  ,i,j,I_MOMX,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMX,YDIR) = num_diff(KS+1,i,j,I_MOMX,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_MOMX,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    !$omp end parallel
    call PROF_rapend  ("NumFilter_Main", 3)

    call PROF_rapstart("NumFilter_Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_MOMX,ZDIR), I_COMM_MOMX_Z )
    if ( .not. TwoD ) &
    call COMM_vars8( num_diff(:,:,:,I_MOMX,XDIR), I_COMM_MOMX_X )
    call COMM_vars8( num_diff(:,:,:,I_MOMX,YDIR), I_COMM_MOMX_Y )

    call PROF_rapend  ("NumFilter_Comm", 3)


    !-----< y-momentum >-----

    call PROF_rapstart("NumFilter_Main", 3)

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+1
    do i = max(IS-2,1), min(IE+2,IA)
    do k = KS, KE
       VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo

    call PROF_rapend  ("NumFilter_Main", 3)

    call calc_numdiff( work, iwork,      & ! (out)
                       VELY,             & ! (in)
                       TwoD,             & ! (in)
                       ND_LAPLACIAN_NUM, & ! (in)
                       0, 0, 1, KE )

    call PROF_rapstart("NumFilter_Main", 3)

    !$omp parallel private(i,j,k) 

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_MOMY,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                   * 0.25_RP * ( DENS(k+1,i,j+1)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)    
    do j = JS, JE
    do i = IS, IE
       num_diff( 1:KS-1,i,j,I_MOMY,ZDIR) = 0.0_RP
       num_diff(KE:KA  ,i,j,I_MOMY,ZDIR) = 0.0_RP
    end do
    end do
    !$omp end do nowait

    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          num_diff(k,i,j,I_MOMY,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                      * 0.25_RP * ( DENS(k,i+1,j+1)+DENS(k,i,j+1)+DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          num_diff(   1:KS-1,i,j,I_MOMY,XDIR) = 0.0_RP
          num_diff(KS  ,i,j,I_MOMY,XDIR) = num_diff(KS  ,i,j,I_MOMY,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_MOMY,XDIR) = num_diff(KS+1,i,j,I_MOMY,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          num_diff(KE+1:KA  ,i,j,I_MOMY,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait
    end if

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_MOMY,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_fdy(j) &
                                   * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_MOMY,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_MOMY,YDIR) = num_diff(KS  ,i,j,I_MOMY,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_MOMY,YDIR) = num_diff(KS+1,i,j,I_MOMY,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_MOMY,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    !$omp end parallel 
    call PROF_rapend  ("NumFilter_Main", 3)

    call PROF_rapstart("NumFilter_Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_MOMY,ZDIR), I_COMM_MOMY_Z )
    if ( .not. TwoD ) &
    call COMM_vars8( num_diff(:,:,:,I_MOMY,XDIR), I_COMM_MOMY_X )
    call COMM_vars8( num_diff(:,:,:,I_MOMY,YDIR), I_COMM_MOMY_Y )

    call PROF_rapend  ("NumFilter_Comm", 3)

    !-----< rho * theta >-----

    call PROF_rapstart("NumFilter_Main", 3)

    if ( ND_USE_RS ) then
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = max(IS-1,1), min(IE+2,IA)
       do k = KS, KE
          pott_diff(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) - REF_pott(k,i,j)
       enddo
       enddo
       enddo
    endif

    call PROF_rapend  ("NumFilter_Main", 3)

    call calc_numdiff( work, iwork,      & ! (out)
                       pott_diff,        & ! (in)
                       TwoD,             & ! (in)
                       ND_LAPLACIAN_NUM, & ! (in)
                       0, 0, 0, KE )

    call PROF_rapstart("NumFilter_Main", 3)

    !$omp parallel private(i,j,k) 

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff(k,i,j,I_RHOT,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                                   * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
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
    !$omp end do nowait

    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          num_diff(k,i,j,I_RHOT,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                      * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          num_diff(   1:KS-1,i,j,I_RHOT,XDIR) = 0.0_RP
          num_diff(KS  ,i,j,I_RHOT,XDIR) = num_diff(KS  ,i,j,I_RHOT,XDIR) * ND_SFC_FACT
          num_diff(KS+1,i,j,I_RHOT,XDIR) = num_diff(KS+1,i,j,I_RHOT,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          num_diff(KE+1:KA  ,i,j,I_RHOT,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait
    end if

    !$omp do OMP_SCHEDULE_ collapse(2)    
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff(k,i,j,I_RHOT,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                                   * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)    
    do j = JS, JE
    do i = IS, IE
       num_diff(   1:KS-1,i,j,I_RHOT,YDIR) = 0.0_RP
       num_diff(KS  ,i,j,I_RHOT,YDIR) = num_diff(KS  ,i,j,I_RHOT,YDIR) * ND_SFC_FACT
       num_diff(KS+1,i,j,I_RHOT,YDIR) = num_diff(KS+1,i,j,I_RHOT,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff(KE+1:KA  ,i,j,I_RHOT,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait
    !$omp end parallel
    call PROF_rapend  ("NumFilter_Main", 3)

    call PROF_rapstart("NumFilter_Comm", 3)

    call COMM_vars8( num_diff(:,:,:,I_RHOT,ZDIR), I_COMM_RHOT_Z )
    if ( .not. TwoD ) &
    call COMM_vars8( num_diff(:,:,:,I_RHOT,XDIR), I_COMM_RHOT_X )
    call COMM_vars8( num_diff(:,:,:,I_RHOT,YDIR), I_COMM_RHOT_Y )

    call COMM_wait ( num_diff(:,:,:,I_DENS,ZDIR), I_COMM_DENS_Z )
    if ( .not. TwoD ) &
    call COMM_wait ( num_diff(:,:,:,I_DENS,XDIR), I_COMM_DENS_X )
    call COMM_wait ( num_diff(:,:,:,I_DENS,YDIR), I_COMM_DENS_Y )
    call COMM_wait ( num_diff(:,:,:,I_MOMZ,ZDIR), I_COMM_MOMZ_Z )
    if ( .not. TwoD ) &
    call COMM_wait ( num_diff(:,:,:,I_MOMZ,XDIR), I_COMM_MOMZ_X )
    call COMM_wait ( num_diff(:,:,:,I_MOMZ,YDIR), I_COMM_MOMZ_Y )
    call COMM_wait ( num_diff(:,:,:,I_MOMX,ZDIR), I_COMM_MOMX_Z )
    if ( .not. TwoD ) &
    call COMM_wait ( num_diff(:,:,:,I_MOMX,XDIR), I_COMM_MOMX_X )
    call COMM_wait ( num_diff(:,:,:,I_MOMX,YDIR), I_COMM_MOMX_Y )
    call COMM_wait ( num_diff(:,:,:,I_MOMY,ZDIR), I_COMM_MOMY_Z )
    if ( .not. TwoD ) &
    call COMM_wait ( num_diff(:,:,:,I_MOMY,XDIR), I_COMM_MOMY_X )
    call COMM_wait ( num_diff(:,:,:,I_MOMY,YDIR), I_COMM_MOMY_Y )
    call COMM_wait ( num_diff(:,:,:,I_RHOT,ZDIR), I_COMM_RHOT_Z )
    if ( .not. TwoD ) &
    call COMM_wait ( num_diff(:,:,:,I_RHOT,XDIR), I_COMM_RHOT_X )
    call COMM_wait ( num_diff(:,:,:,I_RHOT,YDIR), I_COMM_RHOT_Y )

    call PROF_rapend  ("NumFilter_Comm", 3)

    return
  end subroutine ATMOS_DYN_fvm_numfilter_flux

  !-----------------------------------------------------------------------------
  !> Calculate fluxes with numerical filter for tracer variables
  subroutine ATMOS_DYN_fvm_numfilter_flux_q( &
       num_diff_q,                                       &
       DENS, QTRC, is_qv,                                &
       CDZ, CDX, CDY, TwoD, dt,                          &
       REF_qv, iq,                                       &
       ND_COEF, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: num_diff_q(KA,IA,JA,3)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA)
    logical,  intent(in)  :: is_qv

    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: CDX(IA)
    real(RP), intent(in)  :: CDY(JA)

    logical,  intent(in)  :: TwoD
    real(RP), intent(in)  :: dt

    real(RP), intent(in)  :: REF_qv(KA,IA,JA)
    integer,  intent(in)  :: iq

    real(RP), intent(in)  :: ND_COEF
    integer,  intent(in)  :: ND_LAPLACIAN_NUM
    real(RP), intent(in)  :: ND_SFC_FACT
    logical,  intent(in)  :: ND_USE_RS

    real(RP) :: qv_diff(KA,IA,JA) ! anomary of water vapor

    real(RP) :: work(KA,IA,JA,3,2)
    integer  :: iwork

    real(RP) :: diff_coef_tmp
    integer  :: nd_order
    real(RP) :: nd_coef_cdz(KA)
    real(RP) :: nd_coef_cdx(IA)
    real(RP) :: nd_coef_cdy(JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !###########################################################################
    ! 1st order coefficients
    !###########################################################################

    nd_order = ND_LAPLACIAN_NUM * 2
    diff_coef_tmp = - (-1)**(mod(ND_LAPLACIAN_NUM+1,2)) &
                    * ND_COEF / ( 2**(nd_order) * DT )

    do k = KS-1, KE
       nd_coef_cdz(k) = diff_coef_tmp * CDZ(k)**nd_order
    end do
    if ( .not. TwoD ) then
    do i = IS, IE
       nd_coef_cdx(i) = diff_coef_tmp * CDX(i)**nd_order
    end do
    end if
    do j = JS, JE
       nd_coef_cdy(j) = diff_coef_tmp * CDY(j)**nd_order
    end do

    if ( is_qv .AND. (.NOT. ND_USE_RS) ) then

       call PROF_rapstart("NumFilter_Main", 3)

       if ( TwoD ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = JS-1, JE+2
          do k = KS+1, KE-1
             qv_diff(k,IS,j) = ( ( QTRC(k,IS,j)                  ) * 3.0_RP &
                               + ( QTRC(k,IS,j+1)+QTRC(k,IS,j-1) ) * 2.0_RP &
                               + ( QTRC(k,IS,j+2)+QTRC(k,IS,j-2) ) &
                               + ( QTRC(k+1,IS,j)+QTRC(k-1,IS,j) ) * 2.0_RP &
                               ) / 13.0_RP
          enddo
          enddo

          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j  = JS-1, JE+2
             qv_diff(KS,IS,j) = ( ( QTRC(KS,IS,j)                   ) * 3.0_RP &
                                + ( QTRC(KS,IS,j+1)+QTRC(KS,IS,j-1) ) * 2.0_RP &
                                + ( QTRC(KS,IS,j+2)+QTRC(KS,IS,j-2) ) &
                                + ( QTRC(KS+1,IS,j)                 ) * 2.0_RP &
                                ) / 11.0_RP
             qv_diff(KE,IS,j) = ( ( QTRC(KE,IS,j)                   ) * 3.0_RP &
                                + ( QTRC(KE,IS,j+1)+QTRC(KE,IS,j-1) ) * 2.0_RP &
                                + ( QTRC(KE,IS,j+2)+QTRC(KE,IS,j-2) ) &
                                + ( QTRC(KE-1,IS,j)                 ) * 2.0_RP &
                                ) / 11.0_RP
          end do

       else
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

       end if

       call PROF_rapend  ("NumFilter_Main", 3)

       call PROF_rapstart("NumFilter_Comm", 3)

       call COMM_vars8(qv_diff, 1)
       call COMM_wait (qv_diff, 1)

       call PROF_rapend  ("NumFilter_Comm", 3)

    end if

    if ( is_qv ) then

       if ( ND_USE_RS ) then

          call PROF_rapstart("NumFilter_Main", 3)

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JS-1, JE+2
          do i = max(IS-1,1), min(IE+2,IA)
          do k = KS, KE
             qv_diff(k,i,j) = QTRC(k,i,j) - REF_qv(k,i,j)
          enddo
          enddo
          enddo

          call PROF_rapend  ("NumFilter_Main", 3)

       endif

       call calc_numdiff( work, iwork,      & ! (out)
                          qv_diff,          & ! (in)
                          TwoD,             & ! (in)
                          ND_LAPLACIAN_NUM, & ! (in)
                          0, 0, 0, KE )

    else ! not qv

       call calc_numdiff( work, iwork,      & ! (out)
                          QTRC,             & ! (in)
                          TwoD,             & ! (in)
                          ND_LAPLACIAN_NUM, & ! (in)
                          0, 0, 0, KE )

    endif ! QV or not?


    call PROF_rapstart("NumFilter_Main", 3)


    !$omp parallel private(i,j,k)

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       num_diff_q(k,i,j,ZDIR) = work(k,i,j,ZDIR,iwork) * nd_coef_cdz(k) &
                              * 0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff_q(1:KS-2,i,j,ZDIR) = 0.0_RP
       num_diff_q(KS-1,i,j,ZDIR) = work(KS-1,i,j,ZDIR,iwork) * nd_coef_cdz(KS-1) &
                                 * DENS(KS,i,j)
       num_diff_q(KE  ,i,j,ZDIR) = work(KE  ,i,j,ZDIR,iwork) * nd_coef_cdz(KE  ) &
                                 * DENS(KE,i,j)
       num_diff_q(KE+1:KA,i,j,ZDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          num_diff_q(k,i,j,XDIR) = work(k,i,j,XDIR,iwork) * nd_coef_cdx(i) &
                                 * 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          num_diff_q(1:KS-1,i,j,XDIR) = 0.0_RP
          num_diff_q(KS  ,i,j,XDIR) = num_diff_q(KS  ,i,j,XDIR) * ND_SFC_FACT
          num_diff_q(KS+1,i,j,XDIR) = num_diff_q(KS+1,i,j,XDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
          num_diff_q(KE+1:KA,i,j,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait
    end if

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       num_diff_q(k,i,j,YDIR) = work(k,i,j,YDIR,iwork) * nd_coef_cdy(j) &
                              * 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff_q(1:KS-1,i,j,YDIR) = 0.0_RP
       num_diff_q(KS  ,i,j,YDIR) = num_diff_q(KS  ,i,j,YDIR) * ND_SFC_FACT
       num_diff_q(KS+1,i,j,YDIR) = num_diff_q(KS+1,i,j,YDIR) * (1.0_RP + ND_SFC_FACT) * 0.5_RP
       num_diff_q(KE+1:KA,i,j,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    !$omp end parallel
    call PROF_rapend  ("NumFilter_Main", 3)

    call PROF_rapstart("NumFilter_Comm", 3)

    call COMM_vars8( num_diff_q(:,:,:,ZDIR), I_COMM_QTRC_Z )
    if ( .not. TwoD ) &
    call COMM_vars8( num_diff_q(:,:,:,XDIR), I_COMM_QTRC_X )
    call COMM_vars8( num_diff_q(:,:,:,YDIR), I_COMM_QTRC_Y )

    call COMM_wait ( num_diff_q(:,:,:,ZDIR), I_COMM_QTRC_Z )
    if ( .not. TwoD ) &
    call COMM_wait ( num_diff_q(:,:,:,XDIR), I_COMM_QTRC_X )
    call COMM_wait ( num_diff_q(:,:,:,YDIR), I_COMM_QTRC_Y )

    call PROF_rapend  ("NumFilter_Comm", 3)

    return
  end subroutine ATMOS_DYN_fvm_numfilter_flux_q

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_fvm_numfilter_tend( &
       phi_t, &
       phi, &
       rdz, rdx, rdy, &
       TwoD, &
       KO, IO, JO )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    real(RP), intent(out) :: phi_t(KA,IA,JA)
    real(RP), intent(in ) :: phi  (KA,IA,JA)
    real(RP), intent(in ) :: rdz(:)
    real(RP), intent(in ) :: rdx(:)
    real(RP), intent(in ) :: rdy(:)
    logical,  intent(in ) :: TwoD
    integer , intent(in ) :: KO
    integer , intent(in ) :: IO
    integer , intent(in ) :: JO

    real(RP) :: flux(KA,IA,JA,3)

    integer :: k, i, j

    call calc_diff3( flux,      & ! (out)
                     phi,       & ! (in)
                     TwoD,      & ! (in)
                     KO, IO, JO )

    if ( .not. TwoD ) &
    call COMM_vars8( flux(:,:,:,XDIR), 1 )
    call COMM_vars8( flux(:,:,:,YDIR), 2 )
    if ( .not. TwoD ) &
    call COMM_wait ( flux(:,:,:,XDIR), 1 )
    call COMM_wait ( flux(:,:,:,YDIR), 2 )

    if ( TwoD ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          phi_t(k,IS,j) = ( flux(k+KO,IS,j,ZDIR) - flux(k-1+KO,IS,j,ZDIR) ) * RDZ(k) &
                        + ( flux(k,IS,j+JO,YDIR) - flux(k,IS,j-1+JO,YDIR) ) * RDY(j)
       end do
       end do
    else
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          phi_t(k,i,j) = ( flux(k+KO,i,j,ZDIR) - flux(k-1+KO,i,j,ZDIR) ) * RDZ(k) &
                       + ( flux(k,i+IO,j,XDIR) - flux(k,i-1+IO,j,XDIR) ) * RDX(i) &
                       + ( flux(k,i,j+JO,YDIR) - flux(k,i,j-1+JO,YDIR) ) * RDY(j)
       end do
       end do
       end do
    end if

    return
  end subroutine ATMOS_DYN_fvm_numfilter_tend

  subroutine ATMOS_DYN_fvm_apply_numfilter( &
      num_diff,                                          & ! (out)
      DENS, MOMZ, MOMX, MOMY, RHOT,                      & ! (inout)
      CDZ, CDX, CDY, FDZ, FDX, FDY,                      & ! (in)
      RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                & ! (in)
      TwoD, DT,                                          & ! (in)
      GSQRT, MAPF, REF_dens, REF_pott,                   & ! (in)
      ND_COEF, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS  ) ! (in)
   
      implicit none

      real(RP), intent(out) :: num_diff(KA,IA,JA,5,3)

      real(RP), intent(inout)  :: DENS(KA,IA,JA)
      real(RP), intent(inout)  :: MOMZ(KA,IA,JA)
      real(RP), intent(inout)  :: MOMX(KA,IA,JA)
      real(RP), intent(inout)  :: MOMY(KA,IA,JA)
      real(RP), intent(inout)  :: RHOT(KA,IA,JA)
      real(RP), intent(in)  :: CDZ(KA)
      real(RP), intent(in)  :: CDX(IA)
      real(RP), intent(in)  :: CDY(JA)
      real(RP), intent(in)  :: FDZ(KA-1)
      real(RP), intent(in)  :: FDX(IA-1)
      real(RP), intent(in)  :: FDY(JA-1)
      real(RP), intent(in)  :: RCDZ(KA)
      real(RP), intent(in)  :: RCDX(IA)
      real(RP), intent(in)  :: RCDY(JA)
      real(RP), intent(in)  :: RFDZ(KA-1)
      real(RP), intent(in)  :: RFDX(IA-1)
      real(RP), intent(in)  :: RFDY(JA-1)

      logical,  intent(in)  :: TwoD
      real(RP), intent(in)  :: dt

      real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
      real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor

      real(RP), intent(in)  :: REF_dens(KA,IA,JA)
      real(RP), intent(in)  :: REF_pott(KA,IA,JA)

      real(RP), intent(in)  :: ND_COEF
      integer,  intent(in)  :: ND_LAPLACIAN_NUM
      real(RP), intent(in)  :: ND_SFC_FACT
      logical,  intent(in)  :: ND_USE_RS
      !-----------------------------------------------------------

      call ATMOS_DYN_fvm_numfilter_flux( &
         num_diff,                                         &
         DENS, MOMZ, MOMX, MOMY, RHOT,                     &
         CDZ, CDX, CDY, FDZ, FDX, FDY,                     &
         TwoD, dt,                                         &
         REF_dens, REF_pott,                               &
         ND_COEF, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS )
      
      call fvm_add_FluxDiv_xyz( DENS, num_diff(:,:,:,I_DENS,:), RCDZ, RCDX, RCDY, GSQRT, MAPF, TwoD, dt )
      call fvm_add_FluxDiv_xyw( MOMZ, num_diff(:,:,:,I_MOMZ,:), RFDZ, RCDX, RCDY, GSQRT, MAPF, TwoD, dt )
      call fvm_add_FluxDiv_uyz( MOMX, num_diff(:,:,:,I_MOMX,:), RCDZ, RFDX, RCDY, GSQRT, MAPF, TwoD, dt )
      call fvm_add_FluxDiv_xvz( MOMY, num_diff(:,:,:,I_MOMY,:), RCDZ, RCDX, RFDY, GSQRT, MAPF, TwoD, dt )
      call fvm_add_FluxDiv_xyz( RHOT, num_diff(:,:,:,I_RHOT,:), RCDZ, RCDX, RCDY, GSQRT, MAPF, TwoD, dt )
      
      return
   end subroutine ATMOS_DYN_fvm_apply_numfilter

  !-- private subroutines ---------------------------------------------------------------------------

  subroutine calc_numdiff(&
       work, iwork,       &
       data,              &
       TwoD,              &
       nd_laplacian_num,  &
       KO, IO, JO, KEE    )
    use scale_prc, only: PRC_abort
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    
    real(RP), intent(out) :: work(KA,IA,JA,3,2)
    integer,  intent(out) :: iwork
    real(RP), intent(in)  :: data(KA,IA,JA)
    logical,  intent(in)  :: TwoD
    integer,  intent(in)  :: nd_laplacian_num
    integer,  intent(in)  :: KO
    integer,  intent(in)  :: IO
    integer,  intent(in)  :: JO
    integer,  intent(in)  :: KEE

    integer :: i_in, i_out, i_tmp
    integer :: no
    integer :: ho4_itr_num
    !-----------------------------------------------------------------------------

    call PROF_rapstart("NumFilter_Main", 3)

    if (mod(nd_laplacian_num,2)==0) then
      call calc_diff3( work(:,:,:,:,1), & ! (out)
                       data,           & ! (in)
                       TwoD,           & ! (in)
                       KO, IO, JO )      ! (in)
      ho4_itr_num = nd_laplacian_num / 2
    else if (mod(nd_laplacian_num,2)==1 .and. nd_laplacian_num <= 3) then
      call calc_diff1( work(:,:,:,:,1), & ! (out)
                       data,           & ! (in)
                       TwoD,           & ! (in)
                       KO, IO, JO )      ! (in)      
      ho4_itr_num = nd_laplacian_num / 2 + 1
    else
      LOG_ERROR('calc_numdiff',*) 'Numerical filter with nd_laplacian_num=', nd_laplacian_num
      LOG_ERROR_CONT(*) 'is not supported. Check!'
      call PRC_abort    
    end if

    call PROF_rapend  ("NumFilter_Main", 3)

    !###########################################################################
    ! High order coefficients
    !###########################################################################

    i_in  = 1
    i_out = 2

    do no = 2, ho4_itr_num

       call PROF_rapstart("NumFilter_Comm", 3)

       call COMM_vars8( work(:,:,:,ZDIR,i_in),  16 )
       if ( .not. TwoD ) &
       call COMM_vars8( work(:,:,:,XDIR,i_in),  17 )
       call COMM_vars8( work(:,:,:,YDIR,i_in),  18 )

       call COMM_wait ( work(:,:,:,ZDIR,i_in),  16 )
       if ( .not. TwoD ) &
       call COMM_wait ( work(:,:,:,XDIR,i_in),  17 )
       call COMM_wait ( work(:,:,:,YDIR,i_in),  18 )

       call PROF_rapend  ("NumFilter_Comm", 3)

       call PROF_rapstart("NumFilter_Main", 3)

       call calc_diff4( work(:,:,:,:,i_out), & ! (out)
                        work(:,:,:,:,i_in),  & ! (in)
                        CNZ4(:,:,1+KO),      & ! (in)
                        CNX4(:,:,1+IO),      & ! (in)
                        CNY4(:,:,1+JO),      & ! (in)
                        TwoD,                & ! (in)
                        KEE                  ) ! (in)

       call PROF_rapend  ("NumFilter_Main", 3)

       ! swap pointer target
       i_tmp = i_in
       i_in  = i_out
       i_out = i_tmp
    enddo

    iwork = i_in

    return
  end subroutine calc_numdiff

  !-----------------------------------------------------------------------------
  subroutine calc_diff1( &
    diff,                &
    phi,                 &
    TwoD,                &
    KO, IO, JO           )

    implicit none
   
    real(RP), intent(out) :: diff(KA,IA,JA,3)
    real(RP), intent(in ) :: phi(KA,IA,JA)
    logical,  intent(in ) :: TwoD
    integer , intent(in ) :: KO
    integer , intent(in ) :: IO
    integer , intent(in ) :: JO

    integer :: kee
    integer :: k, i, j
    !---------------------------------------------------------------------------

    KEE = KE-KO

    if ( KO == 0 ) then

      !$omp parallel default(none) private(i,j,k)  &
      !$omp shared(JS,JE,IS,IE,KS,KE,phi,diff,CNZ1)

      !$omp do OMP_SCHEDULE_ collapse(2)
      do j = JS, JE
      do i = IS, IE
      do k = KS, KE-1
#ifdef DEBUG
         call CHECK( __LINE__, phi(k+1,i,j) )
         call CHECK( __LINE__, phi(k  ,i,j) )
#endif
         diff(k,i,j,ZDIR) = CNZ1(k,1) * ( phi(k+1,i,j) - phi(k,i,j) )
      enddo
      enddo
      enddo
      !$omp end do nowait
      !$omp do OMP_SCHEDULE_ collapse(2)
      do j = JS, JE
      do i = IS, IE
         diff(KS-1,i,j,ZDIR) = - diff(KS  ,i,j,ZDIR)
         diff(KS-2,i,j,ZDIR) = - diff(KS+1,i,j,ZDIR)
         diff(KE  ,i,j,ZDIR) = - diff(KE-1,i,j,ZDIR)
         diff(KE+1,i,j,ZDIR) = - diff(KE-2,i,j,ZDIR)
         diff(KE+2,i,j,ZDIR) = 0.0_RP
      end do
      end do
      !$omp end do nowait
      !$omp end parallel 
    else ! K0=1

      !$omp parallel default(none) private(i,j,k) &
      !$omp shared(JS,JE,IS,IE,KS,KE,phi,diff,CNZ1)

      !$omp do OMP_SCHEDULE_ collapse(2) 
      do j = JS, JE
      do i = IS, IE
      do k = KS+1, KE-1
#ifdef DEBUG
         call CHECK( __LINE__, phi(k  ,i,j) )
         call CHECK( __LINE__, phi(k-1,i,j) )
#endif
         diff(k,i,j,ZDIR) = CNZ1(k,2) * ( phi(k,i,j) - phi(k-1,i,j) )
      enddo
      enddo
      enddo
      !$omp end do nowait
      !$omp do OMP_SCHEDULE_ collapse(2)
      do j = JS, JE
      do i = IS, IE
         diff(KS  ,i,j,ZDIR) = CNZ1(KS,2) * ( phi(KS,i,j) - phi(KS+1,i,j) )
         diff(KS-1,i,j,ZDIR) = - diff(KS+2,i,j,ZDIR)
         diff(KS-2,i,j,ZDIR) = - diff(KS+3,i,j,ZDIR)
         diff(KE  ,i,j,ZDIR) = CNZ1(KE,2) * (             - phi(KE-1,i,j) )
         diff(KE+1,i,j,ZDIR) = - diff(KE  ,i,j,ZDIR)
         diff(KE+2,i,j,ZDIR) = - diff(KE-1,i,j,ZDIR)
      end do
      end do
      !$omp end do nowait
      !$omp end parallel 
    end if


    !$omp parallel default(none) private(i,j,k)                   &
    !$omp shared(IO,IS,IE,JO,JS,JE,KS,KE,KA,KEE,phi,diff,CNX1,CNY1,TwoD)

    if ( .not. TwoD ) then

       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KEE
#ifdef DEBUG
          call CHECK( __LINE__, phi(k,i+1-IO,j) )
          call CHECK( __LINE__, phi(k,i  -IO,j) )
#endif
          diff(k,i,j,XDIR) = CNX1(i,1+IO) * ( phi(k,i+1-IO,j) - phi(k,i-IO,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          diff(   1:KS-1,i,j,XDIR) = 0.0_RP
          diff(KE+1:KA  ,i,j,XDIR) = 0.0_RP
       enddo
       enddo
       !$omp end do

    end if

    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KEE
#ifdef DEBUG
      call CHECK( __LINE__, phi(k,i,j+1-JO) )
      call CHECK( __LINE__, phi(k,i,j  -JO) )
#endif
      diff(k,i,j,YDIR) = CNY1(j,1+JO) * ( phi(k,i,j+1-JO) - phi(k,i,j-JO) )
    enddo
    enddo
    enddo
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
      diff(   1:KS-1,i,j,YDIR) = 0.0_RP
      diff(KE+1:KA  ,i,j,YDIR) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    !$omp end parallel

    return
  end subroutine calc_diff1

  !-----------------------------------------------------------------------------
  subroutine calc_diff3( &
       diff,             &
       phi,              &
       TwoD,             &
       KO, IO, JO        )
    implicit none
    
    real(RP), intent(out) :: diff(KA,IA,JA,3)
    real(RP), intent(in ) :: phi(KA,IA,JA)
    logical,  intent(in ) :: TwoD
    integer , intent(in ) :: KO
    integer , intent(in ) :: IO
    integer , intent(in ) :: JO

    integer :: kee
    integer :: k, i, j
    !---------------------------------------------------------------------------

    KEE = KE-KO

    if ( KO == 0 ) then

       !$omp parallel default(none) private(i,j,k)  &
       !$omp shared(JS,JE,IS,IE,KS,KE,phi,diff,CNZ3)

       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, phi(k+2,i,j) )
          call CHECK( __LINE__, phi(k+1,i,j) )
          call CHECK( __LINE__, phi(k  ,i,j) )
          call CHECK( __LINE__, phi(k-1,i,j) )
#endif
          diff(k,i,j,ZDIR) = ( + CNZ3(1,k+1,1) * phi(k+2,i,j) &
                               - CNZ3(2,k+1,1) * phi(k+1,i,j) &
                               + CNZ3(3,k+1,1) * phi(k  ,i,j) &
                               - CNZ3(1,k  ,1) * phi(k-1,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
#ifdef DEBUG
          call CHECK( __LINE__, phi(KS+2,i,j) )
          call CHECK( __LINE__, phi(KS+1,i,j) )
          call CHECK( __LINE__, phi(KS,i,j) )
          call CHECK( __LINE__, phi(KE,i,j) )
          call CHECK( __LINE__, phi(KE-1,i,j) )
          call CHECK( __LINE__, phi(KE-2,i,j) )
#endif
          diff(KS,i,j,ZDIR) = ( + CNZ3(1,KS+1,1) * phi(KS+2,i,j) &
                                - CNZ3(2,KS+1,1) * phi(KS+1,i,j) &
                                + CNZ3(3,KS+1,1) * phi(KS  ,i,j) &
                                - CNZ3(1,KS  ,1) * phi(KS+1,i,j) )
          diff(KS-1,i,j,ZDIR) = - diff(KS  ,i,j,ZDIR)
          diff(KS-2,i,j,ZDIR) = - diff(KS+1,i,j,ZDIR)
          diff(KE-1,i,j,ZDIR) = ( + CNZ3(1,KE  ,1) * phi(KE-1,i,j) &
                                  - CNZ3(2,KE  ,1) * phi(KE  ,i,j) &
                                  + CNZ3(3,KE  ,1) * phi(KE-1,i,j) &
                                  - CNZ3(1,KE-1,1) * phi(KE-2,i,j) )
          diff(KE  ,i,j,ZDIR) = - diff(KE-1,i,j,ZDIR)
          diff(KE+1,i,j,ZDIR) = - diff(KE-2,i,j,ZDIR)
          diff(KE+2,i,j,ZDIR) = 0.0_RP
       end do
       end do
       !$omp end do nowait
       !$omp end parallel 
    else ! K0=1

       !$omp parallel default(none) private(i,j,k) &
       !$omp shared(JS,JE,IS,IE,KS,KE,phi,diff,CNZ3)

       !$omp do OMP_SCHEDULE_ collapse(2) 
       do j = JS, JE
       do i = IS, IE
       do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, phi(k+1,i,j) )
          call CHECK( __LINE__, phi(k  ,i,j) )
          call CHECK( __LINE__, phi(k-1,i,j) )
          call CHECK( __LINE__, phi(k-2,i,j) )
#endif
          diff(k,i,j,ZDIR) = ( + CNZ3(1,k  ,2) * phi(k+1,i,j) &
                               - CNZ3(2,k  ,2) * phi(k  ,i,j) &
                               + CNZ3(3,k  ,2) * phi(k-1,i,j) &
                               - CNZ3(1,k-1,2) * phi(k-2,i,j) )
       enddo
       enddo
       enddo
       !$omp end do nowait
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
#ifdef DEBUG
          call CHECK( __LINE__, phi(KS+2,i,j) )
          call CHECK( __LINE__, phi(KS+1,i,j) )
          call CHECK( __LINE__, phi(KS,i,j) )
          call CHECK( __LINE__, phi(KS+1,i,j) )
          call CHECK( __LINE__, phi(KS  ,i,j) )
          call CHECK( __LINE__, phi(KE-1,i,j) )
          call CHECK( __LINE__, phi(KE-2,i,j) )
          call CHECK( __LINE__, phi(KE-3,i,j) )
#endif
          diff(KS+1,i,j,ZDIR) = ( + CNZ3(1,KS+1,2) * phi(KS+2,i,j) &
                                  - CNZ3(2,KS+1,2) * phi(KS+1,i,j) &
                                  + CNZ3(3,KS+1,2) * phi(KS  ,i,j) &
                                  - CNZ3(1,KS  ,2) * phi(KS+1,i,j) )
          diff(KS  ,i,j,ZDIR) = - diff(KS+1,i,j,ZDIR)
          diff(KS-1,i,j,ZDIR) = - diff(KS+2,i,j,ZDIR)
          diff(KS-2,i,j,ZDIR) = - diff(KS+3,i,j,ZDIR)
          diff(KE-1,i,j,ZDIR) = ( - CNZ3(2,KE-1,2) * phi(KE-1,i,j) &
                                  + CNZ3(3,KE-1,2) * phi(KE-2,i,j) &
                                  - CNZ3(1,KE-2,2) * phi(KE-3,i,j) )
          diff(KE  ,i,j,ZDIR) = ( + CNZ3(1,KE  ,2) * phi(KE-1,i,j) &
                                  + CNZ3(3,KE  ,2) * phi(KE-1,i,j) &
                                  - CNZ3(1,KE-1,2) * phi(KE-2,i,j) )
          diff(KE+1,i,j,ZDIR) = - diff(KE  ,i,j,ZDIR)
          diff(KE+2,i,j,ZDIR) = - diff(KE-1,i,j,ZDIR)
       end do
       end do
       !$omp end do nowait
       !$omp end parallel 
    end if

    if ( .not. TwoD ) then
       if ( IO == 0 ) then
          !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IE,KS,KEE,phi,diff,CNX3)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KEE
#ifdef DEBUG
             call CHECK( __LINE__, phi(k,i+2,j) )
             call CHECK( __LINE__, phi(k,i+1,j) )
             call CHECK( __LINE__, phi(k,i  ,j) )
             call CHECK( __LINE__, phi(k,i-1,j) )
#endif
             diff(k,i,j,XDIR) = ( + CNX3(1,i+1,1) * phi(k,i+2,j) &
                                  - CNX3(2,i+1,1) * phi(k,i+1,j) &
                                  + CNX3(3,i+1,1) * phi(k,i  ,j) &
                                  - CNX3(1,i  ,1) * phi(k,i-1,j) )
          enddo
          enddo
          enddo
       else
          !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IE,KS,KEE,phi,diff,CNX3)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KEE
#ifdef DEBUG
             call CHECK( __LINE__, phi(k,i+1,j) )
             call CHECK( __LINE__, phi(k,i  ,j) )
             call CHECK( __LINE__, phi(k,i-1,j) )
             call CHECK( __LINE__, phi(k,i-2,j) )
#endif
             diff(k,i,j,XDIR) = ( + CNX3(1,i  ,2) * phi(k,i+1,j) &
                                  - CNX3(2,i  ,2) * phi(k,i  ,j) &
                                  + CNX3(3,i  ,2) * phi(k,i-1,j) &
                                  - CNX3(1,i-1,2) * phi(k,i-2,j) )
          enddo
          enddo
          enddo
       end if

       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          diff(   1:KS-1,i,j,XDIR) = 0.0_RP
          diff(KE+1:KA  ,i,j,XDIR) = 0.0_RP
       enddo
       enddo
    end if

    if ( JO == 0 ) then
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ &
       !$omp shared(JS,JE,IS,IE,KS,KEE,phi,diff,CNY3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KEE
#ifdef DEBUG
          call CHECK( __LINE__, phi(k,i,j+2) )
          call CHECK( __LINE__, phi(k,i,j+1) )
          call CHECK( __LINE__, phi(k,i,j  ) )
          call CHECK( __LINE__, phi(k,i,j-1) )
#endif
          diff(k,i,j,YDIR) = ( + CNY3(1,j+1,1) * phi(k,i,j+2) &
                               - CNY3(2,j+1,1) * phi(k,i,j+1) &
                               + CNY3(3,j+1,1) * phi(k,i,j  ) &
                               - CNY3(1,j  ,1) * phi(k,i,j-1) )
       enddo
       enddo
       enddo
    else
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ &
       !$omp shared(JS,JE,IS,IE,KS,KEE,phi,diff,CNY3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KEE
#ifdef DEBUG
          call CHECK( __LINE__, phi(k,i,j+1) )
          call CHECK( __LINE__, phi(k,i,j  ) )
          call CHECK( __LINE__, phi(k,i,j-1) )
          call CHECK( __LINE__, phi(k,i,j-2) )
#endif
          diff(k,i,j,YDIR) = ( + CNY3(1,j  ,2) * phi(k,i,j+1) &
                               - CNY3(2,j  ,2) * phi(k,i,j  ) &
                               + CNY3(3,j  ,2) * phi(k,i,j-1) &
                               - CNY3(1,j-1,2) * phi(k,i,j-2) )
       enddo
       enddo
       enddo
    end if

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
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
       TwoD,         &
       k1            )
    implicit none

    real(RP), intent(out) :: num_diff_pt1(KA,IA,JA,3)
    real(RP), intent(in)  :: num_diff_pt0(KA,IA,JA,3)
    real(RP), intent(in)  :: CNZ4(5,KA)
    real(RP), intent(in)  :: CNX4(5,IA)
    real(RP), intent(in)  :: CNY4(5,JA)
    logical,  intent(in)  :: TwoD
    integer,  intent(in)  :: k1

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$omp parallel private(i,j,k) 

    !$omp do OMP_SCHEDULE_ collapse(2)
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
    !$omp end do nowait
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       num_diff_pt1(KS-1,i,j,ZDIR) = - num_diff_pt1(KS  ,i,j,ZDIR)
       num_diff_pt1(KS-2,i,j,ZDIR) = - num_diff_pt1(KS+1,i,j,ZDIR)
       num_diff_pt1(KE  ,i,j,ZDIR) = - num_diff_pt1(KE-1,i,j,ZDIR)
       num_diff_pt1(KE+1,i,j,ZDIR) = - num_diff_pt1(KE-2,i,j,ZDIR)
    enddo
    enddo
    !$omp end do nowait

    if ( .not. TwoD ) then
       !$omp do OMP_SCHEDULE_ collapse(2)
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
       !$omp end do nowait
    end if

    !$omp do OMP_SCHEDULE_ collapse(2)
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
    !$omp end do nowait

    !$omp end parallel

    return
  end subroutine calc_diff4

  !-----------------------------------------------------------------------------

  subroutine fvm_add_FluxDiv_xyz( &
    var, flux, RCDZ, RCDX, RCDY, GSQRT, MAPF, TwoD, dt )
   
    implicit none
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in) :: flux(KA,IA,JA,3)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)   !< map factor
    logical, intent(in) :: TwoD
    real(RP), intent(in) :: dt

    integer :: i, j, k
    real(RP) :: flux_tmp(KA,IA,JA,3)
    real(RP) :: tend_h, tend_v
    !---------------------------------------------

    !$omp parallel private(i,j,k,tend_h,tend_v)

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
      flux_tmp(k,i,j,ZDIR) = GSQRT(k,i,j,I_XYW) * flux(k,i,j,ZDIR) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) )
    enddo
    enddo
    enddo
    !$omp workshare
    flux_tmp(KS-1,:,:,ZDIR) = 0.0_RP
    flux_tmp(KE  ,:,:,ZDIR) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS-1, IE
    do k = KS, KE
      flux_tmp(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) * flux(k,i,j,XDIR)
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j = JS-1, JE
    do i = IS, IE
    do k = KS, KE
      flux_tmp(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) * flux(k,i,j,YDIR)
    enddo
    enddo
    enddo

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      tend_v =  - (flux_tmp(k,i,j,ZDIR) - flux_tmp(k-1,i,j,ZDIR)) * RCDZ(k)
      tend_h =  - (flux_tmp(k,i,j,YDIR) - flux_tmp(k,i,j-1,YDIR)) * RCDY(j)
      if (.not. TwoD) tend_h =  tend_h - (flux_tmp(k,i,j,XDIR) - flux_tmp(k,i-1,j,XDIR)) * RCDX(i)
      var(k,i,j) = var(k,i,j) + dt  * ( tend_v + tend_h ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
    enddo
    enddo
    enddo

    !$omp end parallel

    return
  end subroutine fvm_add_FluxDiv_xyz

  subroutine fvm_add_FluxDiv_xyw( var, flux, RFDZ, RCDX, RCDY, GSQRT, MAPF, TwoD, dt )
    implicit none
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in) :: flux(KA,IA,JA,3)
    real(RP), intent(in) :: RFDZ(KA-1)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)   !< map factor
    logical, intent(in) :: TwoD
    real(RP), intent(in) :: dt

    integer :: i, j, k
    real(RP) :: flux_tmp(KA,IA,JA,3)
    real(RP) :: tend_h, tend_v
    !---------------------------------------------

    !$omp parallel private(i,j,k,tend_h,tend_v)

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE-1
      flux_tmp(k-1,i,j,ZDIR) = GSQRT(k,i,j,I_XYZ) * flux(k,i,j,ZDIR) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) )
    enddo
    enddo
    enddo
    !$omp workshare
    flux_tmp(KS-1,:,:,ZDIR) = 0.0_RP
    flux_tmp(KE-1,:,:,ZDIR) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS-1, IE
    do k = KS, KE-1
      flux_tmp(k,i,j,XDIR) = GSQRT(k,i,j,I_UYW) / MAPF(i,j,2,I_UY) * flux(k,i,j,XDIR)
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j = JS-1, JE
    do i = IS, IE
    do k = KS, KE-1
      flux_tmp(k,i,j,YDIR) = GSQRT(k,i,j,I_XVW) / MAPF(i,j,1,I_XV) * flux(k,i,j,YDIR)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
      tend_v =  - (flux_tmp(k,i,j,ZDIR) - flux_tmp(k-1,i,j,ZDIR)) * RFDZ(k)
      tend_h =  - (flux_tmp(k,i,j,YDIR) - flux_tmp(k,i,j-1,YDIR)) * RCDY(j)
      if (.not. TwoD) tend_h =  tend_h - (flux_tmp(k,i,j,XDIR) - flux_tmp(k,i-1,j,XDIR)) * RCDX(i)
      var(k,i,j) = var(k,i,j) + dt  * ( tend_v + tend_h ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYW)
    enddo
    enddo
    enddo

    !$omp end parallel

    return
  end subroutine fvm_add_FluxDiv_xyw


  subroutine fvm_add_FluxDiv_uyz( var, flux, RCDZ, RFDX, RCDY, GSQRT, MAPF, TwoD, dt )
    implicit none
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in) :: flux(KA,IA,JA,3)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)   !< map factor
    logical, intent(in) :: TwoD
    real(RP), intent(in) :: dt

    integer :: i, j, k
    real(RP) :: flux_tmp(KA,IA,JA,3)
    real(RP) :: tend_h, tend_v
    !---------------------------------------------

    !$omp parallel private(i,j,k,tend_h,tend_v)

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
      flux_tmp(k,i,j,ZDIR) = GSQRT(k,i,j,I_UYW) * flux(k,i,j,ZDIR) / ( MAPF(i,j,1,I_UY)*MAPF(i,j,2,I_UY) )
    enddo
    enddo
    enddo
    !$omp workshare
    flux_tmp(KS-1,:,:,ZDIR) = 0.0_RP
    flux_tmp(KE  ,:,:,ZDIR) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE+1
    do k = KS, KE
      flux_tmp(k,i-1,j,XDIR) = GSQRT(k,i,j,I_XYZ) / MAPF(i,j,2,I_XY) * flux(k,i,j,XDIR)
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j = JS-1, JE
    do i = IS, IE
    do k = KS, KE
      flux_tmp(k,i,j,YDIR) = GSQRT(k,i,j,I_UVZ) / MAPF(i,j,1,I_UV) * flux(k,i,j,YDIR)
    enddo
    enddo
    enddo

    !$omp do collapse(2)
    do j= JS, JE
    do i= IS, min(IE,IEH)
    do k= KS, KE
      tend_v =  - (flux_tmp(k,i,j,ZDIR) - flux_tmp(k-1,i,j,ZDIR)) * RCDZ(k)
      tend_h =  - (flux_tmp(k,i,j,YDIR) - flux_tmp(k,i,j-1,YDIR)) * RCDY(j)
      if (.not. TwoD) tend_h =  tend_h - (flux_tmp(k,i,j,XDIR) - flux_tmp(k,i-1,j,XDIR)) * RFDX(i)
      var(k,i,j) = var(k,i,j) + dt  * ( tend_v + tend_h ) * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) / GSQRT(k,i,j,I_UYZ)
    enddo
    enddo
    enddo

    !$omp end parallel

    return
  end subroutine fvm_add_FluxDiv_uyz

  subroutine fvm_add_FluxDiv_xvz( var, flux, RCDZ, RCDX, RFDY, GSQRT, MAPF, TwoD, dt )
    implicit none
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in) :: flux(KA,IA,JA,3)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RFDY(JA-1)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)   !< map factor
    logical, intent(in) :: TwoD
    real(RP), intent(in) :: dt

    integer :: i, j, k
    real(RP) :: flux_tmp(KA,IA,JA,3)
    real(RP) :: tend_h, tend_v
    !---------------------------------------------
 
    !$omp parallel private(i,j,k,tend_h,tend_v)

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
      flux_tmp(k,i,j,ZDIR) = GSQRT(k,i,j,I_XVW) * flux(k,i,j,ZDIR) / ( MAPF(i,j,1,I_XV)*MAPF(i,j,2,I_XV) )
    enddo
    enddo
    enddo
    !$omp workshare
    flux_tmp(KS-1,:,:,ZDIR) = 0.0_RP
    flux_tmp(KE  ,:,:,ZDIR) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do j = JS, JE
    do i = IS-1, IE
    do k = KS, KE
      flux_tmp(k,i,j,XDIR) = GSQRT(k,i,j,I_UVZ) / MAPF(i,j,2,I_UV) * flux(k,i,j,XDIR)
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j = JS, JE+1
    do i = IS, IE
    do k = KS, KE
      flux_tmp(k,i,j-1,YDIR) = GSQRT(k,i,j,I_XYZ) / MAPF(i,j,1,I_XY) * flux(k,i,j,YDIR)
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j = JS, min(JE,JEH)
    do i = IS, IE
    do k = KS, KE
      tend_v =  - (flux_tmp(k,i,j,ZDIR) - flux_tmp(k-1,i,j,ZDIR)) * RCDZ(k)
      tend_h =  - (flux_tmp(k,i,j,YDIR) - flux_tmp(k,i,j-1,YDIR)) * RFDY(j)
      if (.not. TwoD) &
         tend_h =  tend_h - (flux_tmp(k,i,j,XDIR) - flux_tmp(k,i-1,j,XDIR)) * RCDX(i)
      var(k,i,j) = var(k,i,j) + dt  * ( tend_v + tend_h ) * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) / GSQRT(k,i,j,I_XVZ)
    enddo
    enddo
    enddo

    !$omp end parallel

    return
  end subroutine fvm_add_FluxDiv_xvz

end module scale_atmos_dyn_fvm_numfilter
