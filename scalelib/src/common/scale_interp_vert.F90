!-------------------------------------------------------------------------------
!> module INTERPOLATION vertical
!!
!! @par Description
!!          spacial interpolation module, vertical only (xi<->z,xi->pres)
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_interp_vert
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: INTERP_VERT_setcoef
  public :: INTERP_VERT_xi2z
  public :: INTERP_VERT_xih2zh
  public :: INTERP_VERT_z2xi
  public :: INTERP_VERT_zh2xih

  public :: INTERP_VERT_alloc_pres
  public :: INTERP_VERT_setcoef_pres
  public :: INTERP_VERT_xi2p
  public :: INTERP_VERT_xih2p

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
  integer,  private, allocatable :: INTERP_xi2z_idx   (:,:,:,:) !< index set   for vertical interpolation (xi->z)
  real(RP), private, allocatable :: INTERP_xi2z_coef  (:,:,:,:) !< coefficient for vertical interpolation (xi->z)
  integer,  private, allocatable :: INTERP_z2xi_idx   (:,:,:,:) !< index set   for vertical interpolation (z->xi)
  real(RP), private, allocatable :: INTERP_z2xi_coef  (:,:,:,:) !< coefficient for vertical interpolation (z->xi)

  integer,  private, allocatable :: INTERP_xi2p_idx   (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(RP), private, allocatable :: INTERP_xi2p_coef  (:,:,:,:) !< coefficient for vertical interpolation (xi->p)

  ! half level
  integer,  private, allocatable :: INTERP_xih2zh_idx (:,:,:,:) !< index set   for vertical interpolation (xih->zh)
  real(RP), private, allocatable :: INTERP_xih2zh_coef(:,:,:,:) !< coefficient for vertical interpolation (xih->zh)
  integer,  private, allocatable :: INTERP_zh2xih_idx (:,:,:,:) !< index set   for vertical interpolation (zh->xih)
  real(RP), private, allocatable :: INTERP_zh2xih_coef(:,:,:,:) !< coefficient for vertical interpolation (zh->xih)

  integer,  private, allocatable :: INTERP_xih2p_idx  (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(RP), private, allocatable :: INTERP_xih2p_coef (:,:,:,:) !< coefficient for vertical interpolation (xi->p)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_VERT_setcoef( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       TOPO_exist, &
       Xi, Xih,    &
       Z,  Zh      )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    logical,  intent(in)  :: TOPO_exist
    real(RP), intent(in)  :: Xi (  KA)
    real(RP), intent(in)  :: Xih(0:KA)
    real(RP), intent(in)  :: Z  (  KA,IA,JA)
    real(RP), intent(in)  :: Zh (0:KA,IA,JA)

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[INTERPOLATION] / Categ[ATMOS-RM GRID] / Origin[SCALElib]'
    LOG_INFO("INTERP_VERT_setcoef",*) 'No namelists.'

    INTERP_available = TOPO_exist

    LOG_NEWLINE
    LOG_INFO("INTERP_VERT_setcoef",*) 'Topography exists & interpolation has meaning? : ', INTERP_available


    ! full level

    allocate( INTERP_xi2z_idx (KA,IA,JA,2) )
    allocate( INTERP_xi2z_coef(KA,IA,JA,3) )
    allocate( INTERP_z2xi_idx (KA,IA,JA,2) )
    allocate( INTERP_z2xi_coef(KA,IA,JA,3) )

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,Xi,Zh,Z,INTERP_xi2z_idx,INTERP_xi2z_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( Xi(k) <= Zh(KS-1,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS     ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( Xi(k) <= Z(KS,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KS
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 1.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       elseif( Xi(k) > Zh(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE     ! dummmy
          INTERP_xi2z_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( Xi(k) > Z(KE,i,j) ) then

          INTERP_xi2z_idx (k,i,j,1) = KE
          INTERP_xi2z_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xi2z_coef(k,i,j,1) = 1.0_RP
          INTERP_xi2z_coef(k,i,j,2) = 0.0_RP
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       else

          do kk = KS+1, KE
             kp = kk
             if( Xi(k) <= Z(kk,i,j) ) exit
          enddo

          INTERP_xi2z_idx (k,i,j,1) = kp - 1
          INTERP_xi2z_idx (k,i,j,2) = kp
          INTERP_xi2z_coef(k,i,j,1) = ( Z (kp,i,j) - Xi(k)        ) &
                                    / ( Z (kp,i,j) - Z (kp-1,i,j) )
          INTERP_xi2z_coef(k,i,j,2) = ( Xi(k)      - Z (kp-1,i,j) ) &
                                    / ( Z (kp,i,j) - Z (kp-1,i,j) )
          INTERP_xi2z_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,Z,Xih,Xi,INTERP_z2xi_idx,INTERP_z2xi_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       if ( Z(k,i,j) <= Xih(KS-1) ) then

          INTERP_z2xi_idx (k,i,j,1) = KS     ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KS     ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( Z(k,i,j) <= Xi(KS) ) then

          INTERP_z2xi_idx (k,i,j,1) = KS     ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KS
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 1.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       elseif( Z(k,i,j) > Xih(KE) ) then

          INTERP_z2xi_idx (k,i,j,1) = KE     ! dummmy
          INTERP_z2xi_idx (k,i,j,2) = KE     ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( Z(k,i,j) > Xi(KE) ) then

          INTERP_z2xi_idx (k,i,j,1) = KE
          INTERP_z2xi_idx (k,i,j,2) = KE     ! dummmy
          INTERP_z2xi_coef(k,i,j,1) = 1.0_RP
          INTERP_z2xi_coef(k,i,j,2) = 0.0_RP
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       else

          do kk = KS+1, KE
             kp = kk
             if( Z(k,i,j) <= Xi(kk) ) exit
          enddo

          INTERP_z2xi_idx (k,i,j,1) = kp - 1
          INTERP_z2xi_idx (k,i,j,2) = kp
          INTERP_z2xi_coef(k,i,j,1) = ( Xi(kp)    - Z (k,i,j) ) &
                                    / ( Xi(kp)    - Xi(kp-1)  )
          INTERP_z2xi_coef(k,i,j,2) = ( Z (k,i,j) - Xi(kp-1)  ) &
                                    / ( Xi(kp)    - Xi(kp-1)  )
          INTERP_z2xi_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j) &
    !$omp shared(JA,IA,KS,KE,KA,INTERP_xi2z_idx,INTERP_xi2z_coef,INTERP_z2xi_idx,INTERP_z2xi_coef)
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
    !$omp shared(JA,IA,KS,KE,Xih,Zh,INTERP_xih2zh_idx,INTERP_xih2zh_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE
       if ( Xih(k) < Zh(KS-1,i,j) ) then

          INTERP_xih2zh_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xih2zh_idx (k,i,j,2) = KS     ! dummmy
          INTERP_xih2zh_coef(k,i,j,1) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,2) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( Xih(k) > Zh(KE,i,j) ) then

          INTERP_xih2zh_idx (k,i,j,1) = KE     ! dummmy
          INTERP_xih2zh_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xih2zh_coef(k,i,j,1) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,2) = 0.0_RP
          INTERP_xih2zh_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS, KE
             kp = kk
             if( Xih(k) <= Zh(kk,i,j) ) exit
          enddo

          INTERP_xih2zh_idx (k,i,j,1) = kp - 1
          INTERP_xih2zh_idx (k,i,j,2) = kp
          INTERP_xih2zh_coef(k,i,j,1) = ( Zh (kp,i,j) - Xih(k)        ) &
                                      / ( Zh (kp,i,j) - Zh (kp-1,i,j) )
          INTERP_xih2zh_coef(k,i,j,2) = ( Xih(k)      - Zh (kp-1,i,j) ) &
                                      / ( Zh (kp,i,j) - Zh (kp-1,i,j) )
          INTERP_xih2zh_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j,k,kk,kp) &
    !$omp shared(JA,IA,KS,KE,Xih,Zh,INTERP_zh2xih_idx,INTERP_zh2xih_coef)
    do j = 1, JA
    do i = 1, IA
    do k = KS-1, KE
       if ( Zh(k,i,j) <= Xih(KS-1) ) then

          INTERP_zh2xih_idx (k,i,j,1) = KS     ! dummmy
          INTERP_zh2xih_idx (k,i,j,2) = KS     ! dummmy
          INTERP_zh2xih_coef(k,i,j,1) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,2) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( Zh(k,i,j) > Xih(KE) ) then

          INTERP_zh2xih_idx (k,i,j,1) = KE     ! dummmy
          INTERP_zh2xih_idx (k,i,j,2) = KE     ! dummmy
          INTERP_zh2xih_coef(k,i,j,1) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,2) = 0.0_RP
          INTERP_zh2xih_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS, KE
             kp = kk
             if( Zh(k,i,j) <= Xih(kk) ) exit
          enddo

          INTERP_zh2xih_idx (k,i,j,1) = kp - 1
          INTERP_zh2xih_idx (k,i,j,2) = kp
          INTERP_zh2xih_coef(k,i,j,1) = ( Xih(kp)    - Zh (k,i,j) ) &
                                      / ( Xih(kp)    - Xih(kp-1)  )
          INTERP_zh2xih_coef(k,i,j,2) = ( Zh (k,i,j) - Xih(kp-1)  ) &
                                      / ( Xih(kp)    - Xih(kp-1)  )
          INTERP_zh2xih_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(i,j) &
    !$omp shared(JA,IA,KS,KE,KA,INTERP_xih2zh_idx,INTERP_xih2zh_coef,INTERP_zh2xih_idx,INTERP_zh2xih_coef)
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
  end subroutine INTERP_VERT_setcoef

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_xi2z( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       Xi,         &
       Z,          &
       var,        &
       var_Z       )
    use scale_const, only: &
       CONST_UNDEF
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Xi   (KA)
    real(RP), intent(in)  :: Z    (KA,IA,JA)
    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

    real(RP) :: FDZ(KA)
    real(RP) :: MD (KA)
    real(RP) :: U  (KA)
    real(RP) :: V  (KA)
    real(RP) :: c1, c2, c3, d
    integer  :: kmax

    integer  :: k, i, j, kk
    !---------------------------------------------------------------------------

    kmax = KE - KS + 1

    if ( kmax == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,IS,IE,JS,JE) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xi2z_idx,INTERP_xi2z_coef,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE
       do k = 1, KA
          var_Z(k,i,j) = INTERP_xi2z_coef(k,i,j,1) * var(INTERP_xi2z_idx(k,i,j,1),i,j) &
                       + INTERP_xi2z_coef(k,i,j,2) * var(INTERP_xi2z_idx(k,i,j,2),i,j) &
                       + INTERP_xi2z_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else ! cubic spline

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(FDZ,MD,U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,kmax,IS,IE,JS,JE) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xi2z_idx,INTERP_xi2z_coef,Xi,Z,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE

          do k = KS, KE-1
             FDZ(k) = Z(k+1,i,j) - Z(k,i,j)
          end do

          MD(KS+1) = 2.0 * ( FDZ(KS) + FDZ(KS+1) ) + FDZ(KS)
          do k = KS+2, KE-2
             MD(k) = 2.0 * ( FDZ(k-1) + FDZ(k) )
          end do
          MD(KE-1) = 2.0 * ( FDZ(KE-2) + FDZ(KE-1) ) + FDZ(KE-1)

          do k = KS+1, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / FDZ(k  ) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / FDZ(k-1)
          end do

          call MATRIX_SOLVER_tridiagonal( kmax-2,         & ! [IN]
                                          FDZ(KS+1:KE-1), & ! [IN]
                                          MD (KS+1:KE-1), & ! [IN]
                                          FDZ(KS+1:KE-1), & ! [IN]
                                          V  (KS+1:KE-1), & ! [IN]
                                          U  (KS+1:KE-1)  ) ! [OUT]

          U(KS) = U(KS+1)
          U(KE) = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_xi2z_idx(k,i,j,1),KE-1)

             c3 = ( U(kk+1) - U(kk) ) / FDZ(kk)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / FDZ(kk) - ( U(kk) * 2.0_RP + U(kk+1) ) * FDZ(kk)
             d  = Xi(k) - Z(kk,i,j)

             var_Z(k,i,j) = (        INTERP_xi2z_coef(k,i,j,3) ) * CONST_UNDEF &
                          + ( 1.0_RP-INTERP_xi2z_coef(k,i,j,3) ) * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_VERT_xi2z

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_z2xi( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       Z,          &
       Xi,         &
       var,        &
       var_Xi      )
    use scale_const, only: &
       CONST_UNDEF
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Z     (KA,IA,JA)
    real(RP), intent(in)  :: Xi    (KA)
    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

    real(RP) :: FDZ(KA)
    real(RP) :: MD (KA)
    real(RP) :: U  (KA)
    real(RP) :: V  (KA)
    real(RP) :: c1, c2, c3, d
    integer  :: kmax

    integer :: k, i, j, kk
    !---------------------------------------------------------------------------

    kmax = KE - KS + 1

    if ( kmax <= 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,IS,IE,JS,JE) &
       !$omp shared(var,var_Xi) &
       !$omp shared(INTERP_z2xi_idx,INTERP_z2xi_coef,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE
       do k = 1, KA
          var_Xi(k,i,j) = INTERP_z2xi_coef(k,i,j,1) * var(INTERP_z2xi_idx(k,i,j,1),i,j) &
                        + INTERP_z2xi_coef(k,i,j,2) * var(INTERP_z2xi_idx(k,i,j,2),i,j) &
                        + INTERP_z2xi_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else ! cubic spline

       do k = KS, KE-1
          FDZ(k) = Xi(k+1) - Xi(k)
       end do

       MD(KS+1) = 2.0 * ( FDZ(KS) + FDZ(KS+1) ) + FDZ(KS)
       do k = KS+2, KE-2
          MD(k) = 2.0 * ( FDZ(k-1) + FDZ(k) )
       end do
       MD(KE-1) = 2.0 * ( FDZ(KE-2) + FDZ(KE-1) ) + FDZ(KE-1)

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,kmax,IS,IE,JS,JE) &
       !$omp shared(var,var_Xi,MD,FDZ) &
       !$omp shared(INTERP_z2xi_idx,INTERP_z2xi_coef,Xi,Z,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE

          do k = KS+1, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / FDZ(k  ) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / FDZ(k-1)
          end do

          call MATRIX_SOLVER_tridiagonal( kmax-2,         & ! [IN]
                                          FDZ(KS+1:KE-1), & ! [IN]
                                          MD (KS+1:KE-1), & ! [IN]
                                          FDZ(KS+1:KE-1), & ! [IN]
                                          V  (KS+1:KE-1), & ! [IN]
                                          U  (KS+1:KE-1)  ) ! [OUT]

          U(KS) = U(KS+1)
          U(KE) = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_z2xi_idx(k,i,j,1),KE-1)

             c3 = ( U(kk+1) - U(kk) ) / FDZ(kk)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / FDZ(kk) - ( U(kk) * 2.0_RP + U(kk+1) ) * FDZ(kk)
             d  = Z(k,i,j) - Xi(kk)

             var_Xi(k,i,j) = (        INTERP_z2xi_coef(k,i,j,3) ) * CONST_UNDEF &
                           + ( 1.0_RP-INTERP_z2xi_coef(k,i,j,3) ) * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_VERT_z2xi

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_xih2zh( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       Xih,        &
       Zh,         &
       var,        &
       var_Z       )
    use scale_const, only: &
       CONST_UNDEF
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Xih  (0:KA)
    real(RP), intent(in)  :: Zh   (0:KA,IA,JA)
    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

    real(RP) :: CDZ(KA)
    real(RP) :: MD (KA)
    real(RP) :: U  (KA)
    real(RP) :: V  (KA)
    real(RP) :: c1, c2, c3, d
    integer  :: kmax

    integer  :: k, i, j, kk
    !---------------------------------------------------------------------------

    kmax = KE - KS + 1

    if ( kmax == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,IS,IE,JS,JE) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xih2zh_idx,INTERP_xih2zh_coef,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE
       do k = 1, KA
          var_Z(k,i,j) = INTERP_xih2zh_coef(k,i,j,1) * var(INTERP_xih2zh_idx(k,i,j,1),i,j) &
                       + INTERP_xih2zh_coef(k,i,j,2) * var(INTERP_xih2zh_idx(k,i,j,2),i,j) &
                       + INTERP_xih2zh_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else ! cubic spline

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(CDZ,MD,U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,kmax,IS,IE,JS,JE) &
       !$omp shared(var,var_Z) &
       !$omp shared(INTERP_xih2zh_idx,INTERP_xih2zh_coef,Xih,Zh,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE

          do k = KS, KE
             CDZ(k) = Zh(k,i,j) - Zh(k-1,i,j)
          end do

          MD(KS) = 2.0 * ( CDZ(KS) + CDZ(KS+1) ) + CDZ(KS)
          do k = KS+1, KE-2
             MD(k) = 2.0 * ( CDZ(k) + CDZ(k+1) )
          end do
          MD(KE-1) = 2.0 * ( CDZ(KE-1) + CDZ(KE) ) + CDZ(KE)

          do k = KS, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / CDZ(k+1) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / CDZ(k  )
          end do

          call MATRIX_SOLVER_tridiagonal( kmax-1,         & ! [IN]
                                          CDZ(KS+1:KE  ), & ! [IN]
                                          MD (KS  :KE-1), & ! [IN]
                                          CDZ(KS+1:KE  ), & ! [IN]
                                          V  (KS  :KE-1), & ! [IN]
                                          U  (KS  :KE-1)  ) ! [OUT]

          U(KS-1) = U(KS)
          U(KE)   = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_xih2zh_idx(k,i,j,1),KE-1)

             c3 = ( U(kk+1) - U(kk) ) / CDZ(kk+1)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / CDZ(kk+1) - ( U(kk) * 2.0_RP + U(kk+1) ) * CDZ(kk+1)
             d = Xih(k) - Zh(kk,i,j)

             var_Z(k,i,j) = (        INTERP_xih2zh_coef(k,i,j,3) ) * CONST_UNDEF &
                          + ( 1.0_RP-INTERP_xih2zh_coef(k,i,j,3) ) * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_VERT_xih2zh

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_zh2xih( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       Zh,         &
       Xih,        &
       var,        &
       var_Xi      )
    use scale_const, only: &
       CONST_UNDEF
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Zh   (0:KA,IA,JA)
    real(RP), intent(in)  :: Xih  (0:KA)
    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

    real(RP) :: CDZ(KA)
    real(RP) :: MD (KA)
    real(RP) :: U  (KA)
    real(RP) :: V  (KA)
    real(RP) :: c1, c2, c3, d
    integer  :: kmax

    integer  :: k, i, j, kk
    !---------------------------------------------------------------------------

    kmax = KE - KS + 1

    if ( kmax == 2 ) then

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp shared(KA,IS,IE,JS,JE) &
       !$omp shared(var,var_Xi) &
       !$omp shared(INTERP_zh2xih_idx,INTERP_zh2xih_coef,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE
       do k = 1, KA
          var_Xi(k,i,j) = INTERP_zh2xih_coef(k,i,j,1) * var(INTERP_zh2xih_idx(k,i,j,1),i,j) &
                        + INTERP_zh2xih_coef(k,i,j,2) * var(INTERP_zh2xih_idx(k,i,j,2),i,j) &
                        + INTERP_zh2xih_coef(k,i,j,3) * CONST_UNDEF
       enddo
       enddo
       enddo

    else ! cubic spline

       do k = KS, KE
          CDZ(k) = Xih(k) - Xih(k-1)
       end do

       MD(KS) = 2.0 * ( CDZ(KS) + CDZ(KS+1) ) + CDZ(KS)
       do k = KS+1, KE-2
          MD(k) = 2.0 * ( CDZ(k) + CDZ(k+1) )
       end do
       MD(KE-1) = 2.0 * ( CDZ(KE-1) + CDZ(KE) ) + CDZ(KE)

       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j) &
       !$omp private(U,V,kk,c1,c2,c3,d) &
       !$omp shared(KA,KS,KE,kmax,IS,IE,JS,JE) &
       !$omp shared(var,var_Xi,MD,CDZ) &
       !$omp shared(INTERP_zh2xih_idx,INTERP_zh2xih_coef,Xih,Zh,CONST_UNDEF)
       do j = JS, JE
       do i = IS, IE

          do k = KS, KE-1
             V(k) = ( var(k+1,i,j) - var(k  ,i,j) ) / CDZ(k+1) &
                  - ( var(k  ,i,j) - var(k-1,i,j) ) / CDZ(k  )
          end do

          call MATRIX_SOLVER_tridiagonal( kmax-1,         & ! [IN]
                                          CDZ(KS+1:KE  ), & ! [IN]
                                          MD (KS  :KE-1), & ! [IN]
                                          CDZ(KS+1:KE  ), & ! [IN]
                                          V  (KS  :KE-1), & ! [IN]
                                          U  (KS  :KE-1)  ) ! [OUT]

          U(KS-1) = U(KS)
          U(KE)   = U(KE-1)

          do k = 1, KA
             kk = min(INTERP_zh2xih_idx(k,i,j,1),KE-1)

             c3 = ( U(kk+1) - U(kk) ) / CDZ(kk+1)
             c2 = 3.0_RP * U(kk)
             c1 = ( var(kk+1,i,j) - var(kk,i,j) ) / CDZ(kk+1) - ( U(kk) * 2.0_RP + U(kk+1) ) * CDZ(kk+1)
             d = Zh(k,i,j) - Xih(kk)
             var_Xi(k,i,j) = INTERP_zh2xih_coef(k,i,j,3) * CONST_UNDEF &
                           + ( 1.0_RP - INTERP_zh2xih_coef(k,i,j,3) ) &
                           * ( ( ( c3 * d + c2 ) * d + c1 ) * d + var(kk,i,j) )
          end do

       end do
       end do

    end if

    return
  end subroutine INTERP_VERT_zh2xih

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine INTERP_VERT_alloc_pres( &
       Kpres, IA, JA )
    implicit none

    integer,  intent(in) :: Kpres
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    !---------------------------------------------------------------------------

    allocate( INTERP_xi2p_idx  (Kpres,IA,JA,2) )
    allocate( INTERP_xi2p_coef (Kpres,IA,JA,3) )
    allocate( INTERP_xih2p_idx (Kpres,IA,JA,2) )
    allocate( INTERP_xih2p_coef(Kpres,IA,JA,3) )

    return
  end subroutine INTERP_VERT_alloc_pres

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_setcoef_pres( &
       Kpres,      &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       PRES,       &
       PRESh,      &
       SFC_PRES,   &
       Paxis       )
    implicit none

    integer,  intent(in)  :: Kpres
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: PRES    (  KA,IA,JA) ! pressure in Xi coordinate [Pa]
    real(RP), intent(in)  :: PRESh   (0:KA,IA,JA) ! pressure in Xi coordinate [Pa], layer interface
    real(RP), intent(in)  :: SFC_PRES(     IA,JA) ! surface pressure          [Pa]
    real(RP), intent(in)  :: Paxis   (Kpres)      ! pressure level to output  [Pa]

    real(RP) :: LnPRES    (  KA,IA,JA) ! (log) pressure in Xi coordinate [Pa]
    real(RP) :: LnPRESh   (0:KA,IA,JA) ! (log) pressure in Xi coordinate [Pa], layer interface
    real(RP) :: LnSFC_PRES(     IA,JA) ! (log) surface pressure          [Pa]
    real(RP) :: LnPaxis   (Kpres)      ! (log) pressure level to output  [Pa]

    integer :: k, i, j, kk, kp
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       LnPRES(k,i,j) = log( PRES(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       LnPRESh(k,i,j) = log( PRESh(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       LnSFC_PRES(i,j) = log( SFC_PRES(i,j) )
    enddo
    enddo

    LnPaxis(:) = log( Paxis(:) )

    ! full level

    do j = JS, JE
    do i = IS, IE
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

    ! half level

    do j = JS, JE
    do i = IS, IE
    do k = 1, Kpres
       if ( LnPaxis(k) > LnSFC_PRES(i,j) ) then

          INTERP_xih2p_idx (k,i,j,1) = KS     ! dummmy
          INTERP_xih2p_idx (k,i,j,2) = KS     ! dummmy
          INTERP_xih2p_coef(k,i,j,1) = 0.0_RP
          INTERP_xih2p_coef(k,i,j,2) = 0.0_RP
          INTERP_xih2p_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       elseif( LnPaxis(k) < LnPRESh(KE,i,j) ) then

          INTERP_xih2p_idx (k,i,j,1) = KE     ! dummmy
          INTERP_xih2p_idx (k,i,j,2) = KE     ! dummmy
          INTERP_xih2p_coef(k,i,j,1) = 0.0_RP
          INTERP_xih2p_coef(k,i,j,2) = 0.0_RP
          INTERP_xih2p_coef(k,i,j,3) = 1.0_RP ! set UNDEF

       else

          do kk = KS, KE
             kp = kk
             if( LnPaxis(k) >= LnPRESh(kk,i,j) ) exit
          enddo

          INTERP_xih2p_idx (k,i,j,1) = kp - 1
          INTERP_xih2p_idx (k,i,j,2) = kp
          INTERP_xih2p_coef(k,i,j,1) = ( LnPaxis(k)        - LnPRESh(kp,i,j) ) &
                                     / ( LnPRESh(kp-1,i,j) - LnPRESh(kp,i,j) )
          INTERP_xih2p_coef(k,i,j,2) = ( LnPRESh(kp-1,i,j) - LnPaxis(k)      ) &
                                     / ( LnPRESh(kp-1,i,j) - LnPRESh(kp,i,j) )
          INTERP_xih2p_coef(k,i,j,3) = 0.0_RP

       endif
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_VERT_setcoef_pres

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_xi2p( &
       Kpres,      &
       KA,         &
       IA, IS, IE, &
       JA, JS, JE, &
       var,        &
       var_P       )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in)  :: Kpres
    integer,  intent(in)  :: KA
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: var  (KA   ,IA,JA)
    real(RP), intent(out) :: var_P(Kpres,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(k,i,j) &
    !$omp shared(Kpres,IS,IE,JS,JE) &
    !$omp shared(var,var_P) &
    !$omp shared(INTERP_xi2p_idx,INTERP_xi2p_coef,CONST_UNDEF)
    do j = JS, JE
    do i = IS, IE
    do k = 1, Kpres
       var_P(k,i,j) = INTERP_xi2p_coef(k,i,j,1) * var(INTERP_xi2p_idx(k,i,j,1),i,j) &
                    + INTERP_xi2p_coef(k,i,j,2) * var(INTERP_xi2p_idx(k,i,j,2),i,j) &
                    + INTERP_xi2p_coef(k,i,j,3) * CONST_UNDEF
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_VERT_xi2p

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_xih2p( &
       Kpres,      &
       KA,         &
       IA, IS, IE, &
       JA, JS, JE, &
       var,        &
       var_P       )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in)  :: Kpres
    integer,  intent(in)  :: KA
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: var  (KA   ,IA,JA)
    real(RP), intent(out) :: var_P(Kpres,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(k,i,j) &
    !$omp shared(Kpres,IS,IE,JS,JE) &
    !$omp shared(var,var_P) &
    !$omp shared(INTERP_xih2p_idx,INTERP_xih2p_coef,CONST_UNDEF)
    do j = JS, JE
    do i = IS, IE
    do k = 1, Kpres
       var_P(k,i,j) = INTERP_xih2p_coef(k,i,j,1) * var(INTERP_xih2p_idx(k,i,j,1),i,j) &
                    + INTERP_xih2p_coef(k,i,j,2) * var(INTERP_xih2p_idx(k,i,j,2),i,j) &
                    + INTERP_xih2p_coef(k,i,j,3) * CONST_UNDEF
    enddo
    enddo
    enddo

    return
  end subroutine INTERP_VERT_xih2p

end module scale_interp_vert
