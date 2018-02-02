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
#include "inc_openmp.h"
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
  public :: INTERP_vertical_xih2zh
  public :: INTERP_vertical_z2xi
  public :: INTERP_vertical_zh2xih
  public :: INTERP_setup_pres
  public :: INTERP_update_pres
  public :: INTERP_vertical_xi2p
  public :: INTERP_vertical_xih2p

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: INTERP_available = .false. !< topography exists & vertical interpolation has meaning?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! full level
  integer,  private, allocatable :: INTERP_xi2z_idx (:,:,:,:) !< index set   for vertical interpolation (xi->z)
  real(RP), private, allocatable :: INTERP_xi2z_coef(:,:,:,:) !< coefficient for vertical interpolation (xi->z)
  integer,  private, allocatable :: INTERP_z2xi_idx (:,:,:,:) !< index set   for vertical interpolation (z->xi)
  real(RP), private, allocatable :: INTERP_z2xi_coef(:,:,:,:) !< coefficient for vertical interpolation (z->xi)

  integer,  private, allocatable :: INTERP_xi2p_idx (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(RP), private, allocatable :: INTERP_xi2p_coef(:,:,:,:) !< coefficient for vertical interpolation (xi->p)

  ! half level
  integer,  private, allocatable :: INTERP_xih2zh_idx (:,:,:,:) !< index set   for vertical interpolation (xih->zh)
  real(RP), private, allocatable :: INTERP_xih2zh_coef(:,:,:,:) !< coefficient for vertical interpolation (xih->zh)
  integer,  private, allocatable :: INTERP_zh2xih_idx (:,:,:,:) !< index set   for vertical interpolation (zh->xih)
  real(RP), private, allocatable :: INTERP_zh2xih_coef(:,:,:,:) !< coefficient for vertical interpolation (zh->xih)

  integer,  private, allocatable :: INTERP_xih2p_idx (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(RP), private, allocatable :: INTERP_xih2p_coef(:,:,:,:) !< coefficient for vertical interpolation (xi->p)

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
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[INTERPOLATION] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** No namelists.'

    INTERP_available = TOPO_exist

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Topography exists & interpolation has meaning? : ', INTERP_available


    ! full level

    allocate( INTERP_xi2z_idx (KA,IA,JA,2) )
    allocate( INTERP_xi2z_coef(KA,IA,JA,3) )
    allocate( INTERP_z2xi_idx (KA,IA,JA,2) )
    allocate( INTERP_z2xi_coef(KA,IA,JA,3) )

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,GRID_CZ,REAL_FZ,REAL_CZ,INTERP_xi2z_idx,INTERP_xi2z_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( GRID_CZ(k) <= REAL_FZ(KS-1,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS     ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( GRID_CZ(k) <= REAL_CZ(KS,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 1.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       elseif( GRID_CZ(k) > REAL_FZ(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE     ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( GRID_CZ(k) > REAL_CZ(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE
          INTERP_xi2z_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 1.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

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

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,REAL_CZ,GRID_FZ,GRID_CZ,INTERP_z2xi_idx,INTERP_z2xi_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( REAL_CZ(k,i,j) <= GRID_FZ(KS-1) ) then

          INTERP_z2xi_idx (k,i,j,1) = KS     ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KS     ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( REAL_CZ(k,i,j) <= GRID_CZ(KS) ) then

          INTERP_z2xi_idx (k,i,j,1) = KS     ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KS
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 1.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       elseif( REAL_CZ(k,i,j) > GRID_FZ(KE) ) then

          INTERP_z2xi_idx (k,i,j,1) = KE     ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KE     ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( REAL_CZ(k,i,j) > GRID_CZ(KE) ) then

          INTERP_z2xi_idx (k,i,j,1) = KE
          INTERP_z2xi_idx (k,i,j,2) = KE     ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 1.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

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

    !$omp parallel do default(none) private(i,j) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,INTERP_xi2z_idx,INTERP_xi2z_coef,INTERP_z2xi_idx,INTERP_z2xi_coef,KS,KE,KA)
    do j = 1, JA
    do i = 1, IA
       INTERP_xi2z_idx ( 1:KS-1,i,j,1) = KS     ! dummmy
       INTERP_xi2z_idx ( 1:KS-1,i,j,2) = KS     ! dummmy
       INTERP_xi2z_coef( 1:KS-1,i,j,1) = 0.0_RP
       INTERP_xi2z_coef( 1:KS-1,i,j,2) = 0.0_RP
       INTERP_xi2z_coef( 1:KS-1,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_xi2z_idx (KE+1:KA,i,j,1) = KE     ! dummmy
       INTERP_xi2z_idx (KE+1:KA,i,j,2) = KE     ! dummmy
       INTERP_xi2z_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_xi2z_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_xi2z_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_z2xi_idx ( 1:KS-1,i,j,1) = KS     ! dummmy
       INTERP_z2xi_idx ( 1:KS-1,i,j,2) = KS     ! dummmy
       INTERP_z2xi_coef( 1:KS-1,i,j,1) = 0.0_RP
       INTERP_z2xi_coef( 1:KS-1,i,j,2) = 0.0_RP
       INTERP_z2xi_coef( 1:KS-1,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_z2xi_idx (KE+1:KA,i,j,1) = KE     ! dummmy
       INTERP_z2xi_idx (KE+1:KA,i,j,2) = KE     ! dummmy
       INTERP_z2xi_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_z2xi_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_z2xi_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j) &
    !$omp shared(JA,IA,INTERP_xi2z_idx,INTERP_xi2z_coef,INTERP_z2xi_idx,INTERP_z2xi_coef,KS,KE,KA)
    do j = 1, JA
    do i = 1, IA
       INTERP_xi2z_idx ( 1:KS-1,i,j,1) = KS     ! dummmy
       INTERP_xi2z_idx ( 1:KS-1,i,j,2) = KS     ! dummmy
       INTERP_xi2z_coef( 1:KS-1,i,j,1) = 0.0_RP
       INTERP_xi2z_coef( 1:KS-1,i,j,2) = 0.0_RP
       INTERP_xi2z_coef( 1:KS-1,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_xi2z_idx (KE+1:KA,i,j,1) = KE     ! dummmy
       INTERP_xi2z_idx (KE+1:KA,i,j,2) = KE     ! dummmy
       INTERP_xi2z_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_xi2z_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_xi2z_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_z2xi_idx ( 1:KS-1,i,j,1) = KS     ! dummmy
       INTERP_z2xi_idx ( 1:KS-1,i,j,2) = KS     ! dummmy
       INTERP_z2xi_coef( 1:KS-1,i,j,1) = 0.0_RP
       INTERP_z2xi_coef( 1:KS-1,i,j,2) = 0.0_RP
       INTERP_z2xi_coef( 1:KS-1,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_z2xi_idx (KE+1:KA,i,j,1) = KE     ! dummmy
       INTERP_z2xi_idx (KE+1:KA,i,j,2) = KE     ! dummmy
       INTERP_z2xi_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_z2xi_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_z2xi_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF
    enddo
    enddo


    ! half level

    allocate( INTERP_xih2zh_idx (0:KA,IA,JA,2) )
    allocate( INTERP_xih2zh_coef(0:KA,IA,JA,3) )
    allocate( INTERP_zh2xih_idx (0:KA,IA,JA,2) )
    allocate( INTERP_zh2xih_coef(0:KA,IA,JA,3) )

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,GRID_FZ,REAL_FZ,INTERP_xih2zh_idx,INTERP_xih2zh_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE
       if ( GRID_FZ(k) < REAL_FZ(KS-1,i,j) ) then

          INTERP_xih2zh_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xih2zh_idx (k,i,j,2) = KS     ! dummmy
          INTERP_xih2zh_coef(k,i,j,1) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,2) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( GRID_FZ(k) > REAL_FZ(KE,i,j) ) then

          INTERP_xih2zh_idx (k,i,j,1) = KE     ! dummmy
          INTERP_xih2zh_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xih2zh_coef(k,i,j,1) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,2) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS, KE
             kp = kk
             if( GRID_FZ(k) <= REAL_FZ(kk,i,j) ) exit
          enddo

          INTERP_xih2zh_idx (k,i,j,1) = kp - 1
          INTERP_xih2zh_idx (k,i,j,2) = kp
          INTERP_xih2zh_coef(k,i,j,1) = ( REAL_FZ(kp,i,j) - GRID_FZ(k)        ) &
                                      / ( REAL_FZ(kp,i,j) - REAL_FZ(kp-1,i,j) )
          INTERP_xih2zh_coef(k,i,j,2) = ( GRID_FZ(k)      - REAL_FZ(kp-1,i,j) ) &
                                      / ( REAL_FZ(kp,i,j) - REAL_FZ(kp-1,i,j) )
          INTERP_xih2zh_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,GRID_FZ,REAL_FZ,INTERP_zh2xih_idx,INTERP_zh2xih_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE
       if ( REAL_FZ(k,i,j) <= GRID_FZ(KS-1) ) then

          INTERP_zh2xih_idx (k,i,j,1) = KS     ! dummmy
          INTERP_zh2xih_idx (k,i,j,2) = KS     ! dummmy
          INTERP_zh2xih_coef(k,i,j,1) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,2) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( REAL_FZ(k,i,j) > GRID_FZ(KE) ) then

          INTERP_zh2xih_idx (k,i,j,1) = KE     ! dummmy
          INTERP_zh2xih_idx (k,i,j,2) = KE     ! dummmy
          INTERP_zh2xih_coef(k,i,j,1) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,2) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS, KE
             kp = kk
             if( REAL_FZ(k,i,j) <= GRID_FZ(kk) ) exit
          enddo

          INTERP_zh2xih_idx (k,i,j,1) = kp - 1
          INTERP_zh2xih_idx (k,i,j,2) = kp
          INTERP_zh2xih_coef(k,i,j,1) = ( GRID_FZ(kp)    - REAL_FZ(k,i,j) ) &
                                      / ( GRID_FZ(kp)    - GRID_FZ(kp-1)  )
          INTERP_zh2xih_coef(k,i,j,2) = ( REAL_FZ(k,i,j) - GRID_FZ(kp-1)  ) &
                                      / ( GRID_FZ(kp)    - GRID_FZ(kp-1)  )
          INTERP_zh2xih_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j) &
    !$omp shared(JA,IA,INTERP_xih2zh_idx,INTERP_xih2zh_coef,INTERP_zh2xih_idx,INTERP_zh2xih_coef,KS,KE,KA)
    do j = 1, JA
    do i = 1, IA
       INTERP_xih2zh_idx ( 1:KS-2,i,j,1) = KS     ! dummmy
       INTERP_xih2zh_idx ( 1:KS-2,i,j,2) = KS     ! dummmy
       INTERP_xih2zh_coef( 1:KS-2,i,j,1) = 0.0_RP
       INTERP_xih2zh_coef( 1:KS-2,i,j,2) = 0.0_RP
       INTERP_xih2zh_coef( 1:KS-2,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_xih2zh_idx (KE+1:KA,i,j,1) = KE     ! dummmy
       INTERP_xih2zh_idx (KE+1:KA,i,j,2) = KE     ! dummmy
       INTERP_xih2zh_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_xih2zh_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_xih2zh_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_zh2xih_idx ( 1:KS-2,i,j,1) = KS     ! dummmy
       INTERP_zh2xih_idx ( 1:KS-2,i,j,2) = KS     ! dummmy
       INTERP_zh2xih_coef( 1:KS-2,i,j,1) = 0.0_RP
       INTERP_zh2xih_coef( 1:KS-2,i,j,2) = 0.0_RP
       INTERP_zh2xih_coef( 1:KS-2,i,j,3) = 1.0_RP ! set UNDEF

       INTERP_zh2xih_idx (KE+1:KA,i,j,1) = KE     ! dummmy
       INTERP_zh2xih_idx (KE+1:KA,i,j,2) = KE     ! dummmy
       INTERP_zh2xih_coef(KE+1:KA,i,j,1) = 0.0_RP
       INTERP_zh2xih_coef(KE+1:KA,i,j,2) = 0.0_RP
       INTERP_zh2xih_coef(KE+1:KA,i,j,3) = 1.0_RP ! set UNDEF
    enddo
    enddo

    return
  end subroutine INTERP_setup

  !-----------------------------------------------------------------------------
  subroutine INTERP_vertical_xi2z( &
       var,  &
       var_Z )
    use scale_const, only: &
       CONST_UNDEF
    use scale_grid, only: &
       GRID_CZ
    use scale_grid_real, only: &
       REAL_CZ
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

    real(RP) :: FDZ(KA)
    real(RP) :: MD(KA)
    real(RP) :: U(KA)
    real(RP) :: V(KA)
    real(RP) :: c1, c2, c3
    real(RP) :: d
    integer  :: kk

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( KMAX == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xi2z_idx,INTERP_xi2z_coef,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = 1, KA
          var_Z(k,i,j) = INTERP_xi2z_coef(k,i,j,1) * var(INTERP_xi2z_idx(k,i,j,1),i,j) &
                       + INTERP_xi2z_coef(k,i,j,2) * var(INTERP_xi2z_idx(k,i,j,2),i,j) &
                       + INTERP_xi2z_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else

       ! cubic spline
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(FDZ,MD,U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,KMAX,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xi2z_idx,INTERP_xi2z_coef,GRID_CZ,REAL_CZ,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB

          do k = KS, KE-1
             FDZ(k) = REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j)
          end do

          MD(KS+1) = 2.0 * ( FDZ(KS) + FDZ(KS+1) ) + FDZ(KS)
          do k = KS+2, KE-2
             MD(k) = 2.0 * ( FDZ(k-1) + FDZ(k) )
          end do
          MD(KE-1) = 2.0 * ( FDZ(KE-2) + FDZ(KE-1) ) + FDZ(KE-1)

          do k = KS+1, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / FDZ(k) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / FDZ(k-1)
          end do

          call MATRIX_SOLVER_tridiagonal( &
               KMAX-2, &
               FDZ(KS+1:KE-1), MD(KS+1:KE-1), FDZ(KS+1:KE-1), V(KS+1:KE-1), & ! (in)
               U(KS+1:KE-1) ) ! (out)
          U(KS) = U(KS+1)
          U(KE) = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_xi2z_idx(k,i,j,1), KE-1)
             c3 = ( U(kk+1) - U(kk) ) / FDZ(kk)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / FDZ(kk) - ( U(kk) * 2.0_RP + U(kk+1) ) * FDZ(kk)
             d = GRID_CZ(k) - REAL_CZ(kk,i,j)
             var_Z(k,i,j) = INTERP_xi2z_coef(k,i,j,3) * CONST_UNDEF &
                          + ( 1.0_RP - INTERP_xi2z_coef(k,i,j,3) ) &
                          * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_vertical_xi2z

  !-----------------------------------------------------------------------------
  subroutine INTERP_vertical_z2xi( &
       var,   &
       var_Xi )
    use scale_const, only: &
       CONST_UNDEF
    use scale_grid, only: &
       GRID_CZ, &
       FDZ => GRID_FDZ
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_grid_real, only: &
       REAL_CZ
    implicit none

    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

    real(RP) :: MD(KA)
    real(RP) :: U(KA)
    real(RP) :: V(KA)
    real(RP) :: c1, c2, c3
    real(RP) :: d
    integer  :: kk

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( KMAX == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Xi) &
       !$omp shared(INTERP_z2xi_idx,INTERP_z2xi_coef,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = 1, KA
          var_Xi(k,i,j) = INTERP_z2xi_coef(k,i,j,1) * var(INTERP_z2xi_idx(k,i,j,1),i,j) &
                        + INTERP_z2xi_coef(k,i,j,2) * var(INTERP_z2xi_idx(k,i,j,2),i,j) &
                        + INTERP_z2xi_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else

       ! cubic spline

       MD(KS+1) = 2.0 * ( FDZ(KS) + FDZ(KS+1) ) + FDZ(KS)
       do k = KS+2, KE-2
          MD(k) = 2.0 * ( FDZ(k-1) + FDZ(k) )
       end do
       MD(KE-1) = 2.0 * ( FDZ(KE-2) + FDZ(KE-1) ) + FDZ(KE-1)

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,KMAX,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Xi,MD,FDZ) &
       !$omp shared(INTERP_z2xi_idx,INTERP_z2xi_coef,GRID_CZ,REAL_CZ,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB

          do k = KS+1, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / FDZ(k) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / FDZ(k-1)
          end do

          call MATRIX_SOLVER_tridiagonal( &
               KMAX-2, &
               FDZ(KS+1:KE-1), MD(KS+1:KE-1), FDZ(KS+1:KE-1), V(KS+1:KE-1), & ! (in)
               U(KS+1:KE-1) ) ! (out)
          U(KS) = U(KS+1)
          U(KE) = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_z2xi_idx(k,i,j,1), KE-1)
             c3 = ( U(kk+1) - U(kk) ) / FDZ(kk)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / FDZ(kk) - ( U(kk) * 2.0_RP + U(kk+1) ) * FDZ(kk)
             d = REAL_CZ(k,i,j) - GRID_CZ(kk)
             var_Xi(k,i,j) = INTERP_z2xi_coef(k,i,j,3) * CONST_UNDEF &
                           + ( 1.0_RP - INTERP_z2xi_coef(k,i,j,3) ) &
                           * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_vertical_z2xi

  !-----------------------------------------------------------------------------
  subroutine INTERP_vertical_xih2zh( &
       var,  &
       var_Z )
    use scale_const, only: &
       CONST_UNDEF
    use scale_grid, only: &
       GRID_FZ
    use scale_grid_real, only: &
       REAL_FZ
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

    real(RP) :: CDZ(KA)
    real(RP) :: MD(KA)
    real(RP) :: U(KA)
    real(RP) :: V(KA)
    real(RP) :: c1, c2, c3
    real(RP) :: d
    integer  :: kk

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( KMAX == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xih2zh_idx,INTERP_xih2zh_coef,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = 1, KA
          var_Z(k,i,j) = INTERP_xih2zh_coef(k,i,j,1) * var(INTERP_xih2zh_idx(k,i,j,1),i,j) &
                       + INTERP_xih2zh_coef(k,i,j,2) * var(INTERP_xih2zh_idx(k,i,j,2),i,j) &
                       + INTERP_xih2zh_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else

       ! cubic spline
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(CDZ,MD,U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,KMAX,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xih2zh_idx,INTERP_xih2zh_coef,GRID_FZ,REAL_FZ,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB

          do k = KS, KE
             CDZ(k) = REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j)
          end do

          MD(KS) = 2.0 * ( CDZ(KS) + CDZ(KS+1) ) + CDZ(KS)
          do k = KS+1, KE-2
             MD(k) = 2.0 * ( CDZ(k) + CDZ(k+1) )
          end do
          MD(KE-1) = 2.0 * ( CDZ(KE-1) + CDZ(KE) ) + CDZ(KE)

          do k = KS, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / CDZ(k+1) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / CDZ(k)
          end do

          call MATRIX_SOLVER_tridiagonal( &
               KMAX-1, &
               CDZ(KS+1:KE), MD(KS:KE-1), CDZ(KS+1:KE), V(KS:KE-1), & ! (in)
               U(KS:KE-1) ) ! (out)
          U(KS-1) = U(KS)
          U(KE)   = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_xih2zh_idx(k,i,j,1), KE-1)
             c3 = ( U(kk+1) - U(kk) ) / CDZ(kk+1)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / CDZ(kk+1) - ( U(kk) * 2.0_RP + U(kk+1) ) * CDZ(kk+1)
             d = GRID_FZ(k) - REAL_FZ(kk,i,j)
             var_Z(k,i,j) = INTERP_xih2zh_coef(k,i,j,3) * CONST_UNDEF &
                          + ( 1.0_RP - INTERP_xih2zh_coef(k,i,j,3) ) &
                          * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_vertical_xih2zh

  !-----------------------------------------------------------------------------
  subroutine INTERP_vertical_zh2xih( &
       var,   &
       var_Xi )
    use scale_const, only: &
       CONST_UNDEF
    use scale_grid, only: &
       GRID_FZ, &
       CDZ => GRID_CDZ
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_grid_real, only: &
       REAL_FZ
    implicit none

    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

    real(RP) :: MD(KA)
    real(RP) :: U(KA)
    real(RP) :: V(KA)
    real(RP) :: c1, c2, c3
    real(RP) :: d
    integer  :: kk

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( KMAX == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Xi) &
       !$omp shared(INTERP_zh2xih_idx,INTERP_zh2xih_coef,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB
       do k = 1, KA
          var_Xi(k,i,j) = INTERP_zh2xih_coef(k,i,j,1) * var(INTERP_zh2xih_idx(k,i,j,1),i,j) &
                        + INTERP_zh2xih_coef(k,i,j,2) * var(INTERP_zh2xih_idx(k,i,j,2),i,j) &
                        + INTERP_zh2xih_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else

       ! cubic spline

       MD(KS) = 2.0 * ( CDZ(KS) + CDZ(KS+1) ) + CDZ(KS)
       do k = KS+1, KE-2
          MD(k) = 2.0 * ( CDZ(k) + CDZ(k+1) )
       end do
       MD(KE-1) = 2.0 * ( CDZ(KE-1) + CDZ(KE) ) + CDZ(KE)

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,KMAX,ISB,IEB,JSB,JEB) &
       !$omp shared(var,var_Xi,MD,CDZ) &
       !$omp shared(INTERP_zh2xih_idx,INTERP_zh2xih_coef,GRID_FZ,REAL_FZ,CONST_UNDEF)
       do j = JSB, JEB
       do i = ISB, IEB

          do k = KS, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / CDZ(k+1) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / CDZ(k)
          end do

          call MATRIX_SOLVER_tridiagonal( &
               KMAX-1, &
               CDZ(KS+1:KE), MD(KS:KE-1), CDZ(KS+1:KE), V(KS:KE-1), & ! (in)
               U(KS:KE-1) ) ! (out)
          U(KS-1) = U(KS)
          U(KE)   = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_zh2xih_idx(k,i,j,1), KE-1)
             c3 = ( U(kk+1) - U(kk) ) / CDZ(kk+1)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / CDZ(kk+1) - ( U(kk) * 2.0_RP + U(kk+1) ) * CDZ(kk+1)
             d = REAL_FZ(k,i,j) - GRID_FZ(kk)
             var_Xi(k,i,j) = INTERP_zh2xih_coef(k,i,j,3) * CONST_UNDEF &
                           + ( 1.0_RP - INTERP_zh2xih_coef(k,i,j,3) ) &
                           * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_vertical_zh2xih

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_setup_pres( &
       Kpres )
    implicit none

    integer,  intent(in) :: Kpres
    !---------------------------------------------------------------------------

    allocate( INTERP_xi2p_idx (Kpres,IA,JA,2) )
    allocate( INTERP_xi2p_coef(Kpres,IA,JA,3) )

    return
  end subroutine INTERP_setup_pres

  !-----------------------------------------------------------------------------
  subroutine INTERP_update_pres( &
       Kpres,    &
       PRES,     &
       SFC_PRES, &
       Paxis     )
    implicit none

    integer,  intent(in) :: Kpres
    real(RP), intent(in) :: PRES    (KA,IA,JA) ! pressure in Xi coordinate [Pa]
    real(RP), intent(in) :: SFC_PRES(   IA,JA) ! surface pressure          [Pa]
    real(RP), intent(in) :: Paxis   (Kpres)    ! pressure level to output  [Pa]

    real(RP) :: LnPRES    (KA,IA,JA) ! (log) pressure in Xi coordinate [Pa]
    real(RP) :: LnSFC_PRES(   IA,JA) ! (log) surface pressure          [Pa]
    real(RP) :: LnPaxis   (Kpres)    ! (log) pressure level to output  [Pa]

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       LnPRES(k,i,j) = log( PRES(k,i,j) )
    enddo
    enddo
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
       LnSFC_PRES(i,j) = log( SFC_PRES(i,j) )
    enddo
    enddo

    LnPaxis(:) = log( Paxis(:) )

    do j = JSB, JEB
    do i = ISB, IEB
    do k = 1, Kpres
       if ( LnPaxis(k) >= LnSFC_PRES(i,j) ) then

          INTERP_xi2p_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xi2p_idx (k,i,j,2) = KS     ! dummmy
          INTERP_xi2p_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2p_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2p_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( LnPaxis(k) >= LnPRES(KS,i,j) ) then

          INTERP_xi2p_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xi2p_idx (k,i,j,2) = KS
          INTERP_xi2p_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2p_coef(k,i,j,2) = 1.0_RP
          INTERP_xi2p_coef(k,i,j,3) = 0.0_RP

       elseif( LnPaxis(k) < LnPRES(KE,i,j) ) then

          INTERP_xi2p_idx (k,i,j,1) = KE     ! dummmy
          INTERP_xi2p_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xi2p_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2p_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2p_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS+1, KE
             kp = kk
             if( LnPaxis(k) >= LnPRES(kk,i,j) ) exit
          enddo

          INTERP_xi2p_idx (k,i,j,1) = kp - 1
          INTERP_xi2p_idx (k,i,j,2) = kp
          INTERP_xi2p_coef(k,i,j,1) = ( LnPaxis(k)        - LnPRES (kp,i,j) ) &
                                    / ( LnPRES (kp-1,i,j) - LnPRES (kp,i,j) )
          INTERP_xi2p_coef(k,i,j,2) = ( LnPRES (kp-1,i,j) - LnPaxis(k)      ) &
                                    / ( LnPRES (kp-1,i,j) - LnPRES (kp,i,j) )
          INTERP_xi2p_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_update_pres

  !-----------------------------------------------------------------------------
  subroutine INTERP_vertical_xi2p( &
       Kpres, &
       var,   &
       var_P  )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in)  :: Kpres
    real(RP), intent(in)  :: var  (KA   ,IA,JA)
    real(RP), intent(out) :: var_P(Kpres,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = 1, Kpres
       var_P(k,i,j) = INTERP_xi2p_coef(k,i,j,1) * var(INTERP_xi2p_idx(k,i,j,1),i,j) &
                    + INTERP_xi2p_coef(k,i,j,2) * var(INTERP_xi2p_idx(k,i,j,2),i,j) &
                    + INTERP_xi2p_coef(k,i,j,3) * CONST_UNDEF
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_vertical_xi2p

  !-----------------------------------------------------------------------------
  subroutine INTERP_vertical_xih2p( &
       Kpres, &
       var,   &
       var_P  )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in)  :: Kpres
    real(RP), intent(in)  :: var  (KA   ,IA,JA)
    real(RP), intent(out) :: var_P(Kpres,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JSB, JEB
    do i = ISB, IEB
    do k = 1, Kpres
       var_P(k,i,j) = INTERP_xih2p_coef(k,i,j,1) * var(INTERP_xih2p_idx(k,i,j,1),i,j) &
                    + INTERP_xih2p_coef(k,i,j,2) * var(INTERP_xih2p_idx(k,i,j,2),i,j) &
                    + INTERP_xih2p_coef(k,i,j,3) * CONST_UNDEF
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_vertical_xih2p

end module scale_interpolation
