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
module scale_interpolation
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: INTERP_setup
  public :: INTERP_vertical_xi2z
  public :: INTERP_vertical_z2xi

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: INTERP_available = .false. !< vertical interpolation is available? (have meaning?)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, allocatable :: INTERP_xi2z_idx (:,:,:,:) !< index set   for vertical interpolation (xi->z)
  real(RP), private, allocatable :: INTERP_xi2z_coef(:,:,:,:) !< coefficient for vertical interpolation (xi->z)
  integer,  private, allocatable :: INTERP_z2xi_idx (:,:,:,:) !< index set   for vertical interpolation (z->xi)
  real(RP), private, allocatable :: INTERP_z2xi_coef(:,:,:,:) !< coefficient for vertical interpolation (z->xi)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_setup
    use scale_grid, only: &
       GRID_CZ, &
       GRID_FZ
    use scale_topography, only: &
       TOPO_exist
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    implicit none

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[INTERPOLATION]/Categ[COMMON]'

    INTERP_available = TOPO_exist
    if( IO_L ) write(IO_FID_LOG,*) '*** Interpolation is available? : ', INTERP_available

    allocate( INTERP_xi2z_idx (KA,IA,JA,2) )
    allocate( INTERP_xi2z_coef(KA,IA,JA,3) )
    allocate( INTERP_z2xi_idx (KA,IA,JA,2) )
    allocate( INTERP_z2xi_coef(KA,IA,JA,3) )

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( GRID_CZ(k) <= REAL_FZ(KS-1,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS   ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS   ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( GRID_CZ(k) <= REAL_CZ(KS,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS   ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 1.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       elseif( GRID_CZ(k) > REAL_CZ(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE
          INTERP_xi2z_idx (k,i,j,2) = KE   ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 1.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       elseif( GRID_CZ(k) > REAL_FZ(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE   ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KE   ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS+1, KE
             kp = kk
             if( GRID_CZ(k) <= REAL_CZ(kk,i,j) ) exit
          enddo

          INTERP_xi2z_idx (k,i,j,1) = kp - 1
          INTERP_xi2z_idx (k,i,j,2) = kp
          INTERP_xi2z_coef(k,i,j,1) = ( REAL_CZ(kp,i,j) - GRID_CZ(k)        ) &
                                    / ( REAL_CZ(kp,i,j) - REAL_CZ(kp-1,i,j) )
          INTERP_xi2z_coef(k,i,j,2) = ( GRID_CZ(k)      - REAL_CZ(kp-1,i,j) ) &
                                    / ( REAL_CZ(kp,i,j) - REAL_CZ(kp-1,i,j) )
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( REAL_CZ(k,i,j) <= GRID_FZ(KS-1) ) then

          INTERP_z2xi_idx (k,i,j,1) = KS   ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KS   ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( REAL_CZ(k,i,j) <= GRID_CZ(KS) ) then

          INTERP_z2xi_idx (k,i,j,1) = KS   ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KS
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 1.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       elseif( REAL_CZ(k,i,j) > GRID_CZ(KE) ) then

          INTERP_z2xi_idx (k,i,j,1) = KE
          INTERP_z2xi_idx (k,i,j,2) = KE   ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 1.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       elseif( REAL_CZ(k,i,j) > GRID_FZ(KE) ) then

          INTERP_z2xi_idx (k,i,j,1) = KE   ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KE   ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS+1, KE
             kp = kk
             if( REAL_CZ(k,i,j) <= GRID_CZ(kk) ) exit
          enddo

          INTERP_z2xi_idx (k,i,j,1) = kp - 1
          INTERP_z2xi_idx (k,i,j,2) = kp
          INTERP_z2xi_coef(k,i,j,1) = ( GRID_CZ(kp)    - REAL_CZ(k,i,j) ) &
                                    / ( GRID_CZ(kp)    - GRID_CZ(kp-1)  )
          INTERP_z2xi_coef(k,i,j,2) = ( REAL_CZ(k,i,j) - GRID_CZ(kp-1)  ) &
                                    / ( GRID_CZ(kp)    - GRID_CZ(kp-1)  )
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA
       INTERP_xi2z_idx ( 1:KS-1,i,j,1) = KS   ! dummmy
       INTERP_xi2z_idx ( 1:KS-1,i,j,2) = KS   ! dummmy
       INTERP_xi2z_coef( 1:KS-1,i,j,1) = 0.0_RP
       INTERP_xi2z_coef( 1:KS-1,i,j,2) = 0.0_RP
       INTERP_xi2z_coef( 1:KS-1,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_xi2z_idx (KE+1:KA,i,j,1) = KE   ! dummmy
       INTERP_xi2z_idx (KE+1:KA,i,j,2) = KE   ! dummmy
       INTERP_xi2z_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_xi2z_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_xi2z_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_z2xi_idx ( 1:KS-1,i,j,1) = KS   ! dummmy
       INTERP_z2xi_idx ( 1:KS-1,i,j,2) = KS   ! dummmy
       INTERP_z2xi_coef( 1:KS-1,i,j,1) = 0.0_RP
       INTERP_z2xi_coef( 1:KS-1,i,j,2) = 0.0_RP
       INTERP_z2xi_coef( 1:KS-1,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_z2xi_idx (KE+1:KA,i,j,1) = KE   ! dummmy
       INTERP_z2xi_idx (KE+1:KA,i,j,2) = KE   ! dummmy
       INTERP_z2xi_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_z2xi_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_z2xi_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF
    enddo
    enddo

    return
  end subroutine INTERP_setup

  !-----------------------------------------------------------------------------
  !> Reset random seed
  subroutine INTERP_vertical_xi2z( &
       var,  &
       var_Z )
    use scale_const, only: &
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
  end subroutine INTERP_vertical_xi2z

  !-----------------------------------------------------------------------------
  !> Reset random seed
  subroutine INTERP_vertical_z2xi( &
       var,   &
       var_Xi )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       var_Xi(k,i,j) = INTERP_z2xi_coef(k,i,j,1) * var(INTERP_z2xi_idx(k,i,j,1),i,j) &
                     + INTERP_z2xi_coef(k,i,j,2) * var(INTERP_z2xi_idx(k,i,j,2),i,j) &
                     + INTERP_z2xi_coef(k,i,j,3) * CONST_UNDEF
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_vertical_z2xi

end module scale_interpolation
!-------------------------------------------------------------------------------
