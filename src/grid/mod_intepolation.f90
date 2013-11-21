!-------------------------------------------------------------------------------
!> module INTERPOLATION
!!
!! @par Description
!!          spacial interpolation module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-08-17 (H.Yashiro)  [new]
!!
!<
module mod_interpolation
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: INTERP_setup
  public :: INTERP_vertical

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
  integer,  private, save :: INTERP_xi2z_idx (KA,IA,JA,2) !< index set   for vertical interpolation
  real(RP), private, save :: INTERP_xi2z_coef(KA,IA,JA,3) !< coefficient for vertical interpolation

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_setup
    use mod_grid, only: &
       GRID_CZ
    use mod_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    implicit none

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[INTERPOLATION]/Categ[COMMON]'

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( GRID_CZ(k) <= REAL_FZ(KS-1,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS   ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS   ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.D0
          INTERP_xi2z_coef(k,i,j,2) = 0.D0
          INTERP_xi2z_coef(k,i,j,3) = 1.D0 ! set UNDEF

       elseif( GRID_CZ(k) <= REAL_CZ(KS,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS   ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS
          INTERP_xi2z_coef(k,i,j,1) = 0.D0
          INTERP_xi2z_coef(k,i,j,2) = 1.D0
          INTERP_xi2z_coef(k,i,j,3) = 0.D0

       elseif( GRID_CZ(k) >= REAL_CZ(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE
          INTERP_xi2z_idx (k,i,j,2) = KE   ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 1.D0
          INTERP_xi2z_coef(k,i,j,2) = 0.D0
          INTERP_xi2z_coef(k,i,j,3) = 0.D0

       elseif( GRID_CZ(k) >= REAL_FZ(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE   ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KE   ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.D0
          INTERP_xi2z_coef(k,i,j,2) = 0.D0
          INTERP_xi2z_coef(k,i,j,3) = 1.D0 ! set UNDEF

       else

          do kk = KS+1, KE
             kp = kk
             if( GRID_CZ(k) < REAL_CZ(kk,i,j) ) exit
          enddo

          INTERP_xi2z_idx (k,i,j,1) = kp - 1
          INTERP_xi2z_idx (k,i,j,2) = kp
          INTERP_xi2z_coef(k,i,j,1) = ( REAL_CZ(kp,i,j) - GRID_CZ(k)        ) &
                                    / ( REAL_CZ(kp,i,j) - REAL_CZ(kp-1,i,j) )
          INTERP_xi2z_coef(k,i,j,2) = ( GRID_CZ(k)      - REAL_CZ(kp-1,i,j) ) &
                                    / ( REAL_CZ(kp,i,j) - REAL_CZ(kp-1,i,j) )
          INTERP_xi2z_coef(k,i,j,3) = 0.D0

       endif
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       INTERP_xi2z_idx ( 1:KS-1,i,j,1) = KS   ! dummmy
       INTERP_xi2z_idx ( 1:KS-1,i,j,2) = KS   ! dummmy
       INTERP_xi2z_coef( 1:KS-1,i,j,1) = 0.D0
       INTERP_xi2z_coef( 1:KS-1,i,j,2) = 0.D0
       INTERP_xi2z_coef( 1:KS-1,i,j,3) = 1.D0 ! set UNDEF

       INTERP_xi2z_idx (KE+1:KA,i,j,1) = KE   ! dummmy
       INTERP_xi2z_idx (KE+1:KA,i,j,2) = KE   ! dummmy
       INTERP_xi2z_coef(KE+1:KA,i,j,1) = 0.D0
       INTERP_xi2z_coef(KE+1:KA,i,j,2) = 0.D0
       INTERP_xi2z_coef(KE+1:KA,i,j,3) = 1.D0 ! set UNDEF
    enddo
    enddo

    return
  end subroutine INTERP_setup

  !-----------------------------------------------------------------------------
  !> Reset random seed
  subroutine INTERP_vertical( &
       var,  &
       var_Z )
    use mod_const, only: &
       CONST_UNDEF
    implicit none

    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       var_Z(k,i,j) = INTERP_xi2z_coef(k,i,j,1) * var(INTERP_xi2z_idx(k,i,j,1),i,j) &
                    + INTERP_xi2z_coef(k,i,j,2) * var(INTERP_xi2z_idx(k,i,j,2),i,j) &
                    + INTERP_xi2z_coef(k,i,j,3) * CONST_UNDEF

!       if ( i == IS .AND. j == JS ) then
!       if( IO_L ) write(IO_FID_LOG,*) k,i,j, &
!                                      var_Z(k,i,j), &
!                                      INTERP_xi2z_idx(k,i,j,1), &
!                                      INTERP_xi2z_idx(k,i,j,2), &
!                                      var(INTERP_xi2z_idx(k,i,j,1),i,j), &
!                                      var(INTERP_xi2z_idx(k,i,j,2),i,j)
!       endif
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_vertical

end module mod_interpolation
!-------------------------------------------------------------------------------
