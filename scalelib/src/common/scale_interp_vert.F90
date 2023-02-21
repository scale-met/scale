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
  use scale_io
  use scale_prof
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
  public :: INTERP_VERT_dealloc_pres
  public :: INTERP_VERT_setcoef_pres
  public :: INTERP_VERT_xi2p
  public :: INTERP_VERT_xih2p
  public :: INTERP_VERT_finalize

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
  real(RP), private, allocatable :: INTERP_xi2z_coef  (:,  :,:) !< coefficient for vertical interpolation (xi->z)
  integer,  private, allocatable :: INTERP_z2xi_idx   (:,:,:,:) !< index set   for vertical interpolation (z->xi)
  real(RP), private, allocatable :: INTERP_z2xi_coef  (:,  :,:) !< coefficient for vertical interpolation (z->xi)

  integer,  private, allocatable :: INTERP_xi2p_idx   (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(RP), private, allocatable :: INTERP_xi2p_coef  (:,  :,:) !< coefficient for vertical interpolation (xi->p)

  ! half level
  integer,  private, allocatable :: INTERP_xih2zh_idx (:,:,:,:) !< index set   for vertical interpolation (xih->zh)
  real(RP), private, allocatable :: INTERP_xih2zh_coef(:,  :,:) !< coefficient for vertical interpolation (xih->zh)
  integer,  private, allocatable :: INTERP_zh2xih_idx (:,:,:,:) !< index set   for vertical interpolation (zh->xih)
  real(RP), private, allocatable :: INTERP_zh2xih_coef(:,  :,:) !< coefficient for vertical interpolation (zh->xih)

  integer,  private, allocatable :: INTERP_xih2p_idx  (:,:,:,:) !< index set   for vertical interpolation (xi->p)
  real(RP), private, allocatable :: INTERP_xih2p_coef (:,  :,:) !< coefficient for vertical interpolation (xi->p)

  ! log pressure
  real(RP), private, allocatable :: LnPRES (:,:,:)
  real(RP), private, allocatable :: LnPRESh(:,:,:)
  real(RP), private, allocatable :: LnPaxis(:)


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
    use scale_interp, only: &
       INTERP_factor1d
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    logical,  intent(in)  :: TOPO_exist
    real(RP), intent(in)  :: Xi (  KA)
    real(RP), intent(in)  :: Xih(0:KA)
    real(RP), intent(in)  :: Z  (  KA,IA,JA)
    real(RP), intent(in)  :: Zh (0:KA,IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("INTERP_VERT_setcoef",*) 'Setup'
    LOG_INFO("INTERP_VERT_setcoef",*) 'No namelists.'

    INTERP_available = TOPO_exist

    LOG_NEWLINE
    LOG_INFO("INTERP_VERT_setcoef",*) 'Topography exists & interpolation has meaning? : ', INTERP_available

    ! full level

    allocate( INTERP_xi2z_idx (KA,2,IA,JA) )
    allocate( INTERP_xi2z_coef(KA,  IA,JA) )
    allocate( INTERP_z2xi_idx (KA,2,IA,JA) )
    allocate( INTERP_z2xi_coef(KA,  IA,JA) )
    !$acc enter data create(INTERP_xi2z_idx, INTERP_xi2z_coef, INTERP_z2xi_idx, INTERP_z2xi_coef)

    !$acc data copyin(Z, Xi)

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KA,KS,KE,Xi,Z,INTERP_xi2z_idx,INTERP_xi2z_coef)
    !$acc kernels
    !$acc loop independent
    do j = 1, JA
    !$acc loop independent
    do i = 1, IA
       call INTERP_factor1d( KA, KS, KE, KA, KS, KE, &
                             Z(:,i,j), Xi(:),           & ! (in)
                             INTERP_xi2z_idx (:,:,i,j), & ! (out)
                             INTERP_xi2z_coef(:,  i,j), & ! (out)
                             flag_extrap = .true.       ) ! (in)
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KA,KS,KE,Z,Xi,INTERP_z2xi_idx,INTERP_z2xi_coef)
    !$acc kernels
    !$acc loop independent
    do j = 1, JA
    !$acc loop independent
    do i = 1, IA
       call INTERP_factor1d( KA, KS, KE, KA, KS, KE, &
                             Xi(:), Z(:,i,j),           & ! (in)
                             INTERP_z2xi_idx (:,:,i,j), & ! (out)
                             INTERP_z2xi_coef(:,  i,j), & ! (out)
                             flag_extrap = .true.       ) ! (in)
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    ! half level

    allocate( INTERP_xih2zh_idx (KA,2,IA,JA) )
    allocate( INTERP_xih2zh_coef(KA,  IA,JA) )
    allocate( INTERP_zh2xih_idx (KA,2,IA,JA) )
    allocate( INTERP_zh2xih_coef(KA,  IA,JA) )
    !$acc enter data create(INTERP_xih2zh_idx, INTERP_xih2zh_coef, INTERP_zh2xih_idx, INTERP_zh2xih_coef)

    !$acc data copyin(Zh, Xih)

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KA,KS,KE,Xih,Zh,INTERP_xih2zh_idx,INTERP_xih2zh_coef)
    !$acc kernels
    !$acc loop independent
    do j = 1, JA
    !$acc loop independent
    do i = 1, IA
       call INTERP_factor1d( KA, KS-1, KE, KA, KS-1, KE, &
                             Zh(1:,i,j), Xih(1:),         & ! (in)
                             INTERP_xih2zh_idx (:,:,i,j), & ! (out)
                             INTERP_xih2zh_coef(:,  i,j), & ! (out)
                             flag_extrap = .false.        ) ! (in)
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KA,KS,KE,Xih,Zh,INTERP_zh2xih_idx,INTERP_zh2xih_coef)
    !$acc kernels
    !$acc loop independent
    do j = 1, JA
    !$acc loop independent
    do i = 1, IA
       call INTERP_factor1d( KA, KS-1, KE, KA, KS-1, KE, &
                             Xih(1:), Zh(1:,i,j),         & ! (in)
                             INTERP_zh2xih_idx (:,:,i,j), & ! (out)
                             INTERP_zh2xih_coef(:,  i,j), & ! (out)
                             flag_extrap = .true.         ) ! (in)
    enddo
    enddo
    !$acc end kernels

    !$acc end data

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
    use scale_interp, only: &
       INTERP_interp1d
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Xi   (KA)
    real(RP), intent(in)  :: Z    (KA,IA,JA)
    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

#ifdef _OPENACC
    real(RP) :: workr(KA,7)
    integer  :: worki(KA,2)
#endif

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE) &
    !$omp shared(Xi,Z,var,var_Z) &
    !$omp shared(INTERP_xi2z_idx,INTERP_xi2z_coef)
    !$acc kernels copyin(Xi, Z, var) copyout(var_Z)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent private(workr,worki)
    do i = IS, IE
       call INTERP_interp1d( KA, KS, KE, KA, KS, KE, &
#ifdef _OPENACC
                             workr(:,:), worki(:,:), &
#endif
                             INTERP_xi2z_idx (:,:,i,j), & ! (in)
                             INTERP_xi2z_coef(:,  i,j), & ! (in)
                             Z(:,i,j), Xi(:),           & ! (in)
                             var(:,i,j),                & ! (in)
                             var_Z(:,i,j)               ) ! (out)
    end do
    end do
    !$acc end kernels

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
    use scale_interp, only: &
       INTERP_interp1d
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Z     (KA,IA,JA)
    real(RP), intent(in)  :: Xi    (KA)
    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

#ifdef _OPENACC
    real(RP) :: workr(KA,7)
    integer  :: worki(KA,2)
#endif

    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE) &
    !$omp shared(Z,Xi,var,var_Xi) &
    !$omp shared(INTERP_z2xi_idx,INTERP_z2xi_coef)
    !$acc kernels copyin(Z, Xi, var) copyout(var_Xi)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent private(workr,worki)
    do i = IS, IE
       call INTERP_interp1d( KA, KS, KE, KA, KS, KE, &
#ifdef _OPENACC
                             workr(:,:), worki(:,:), &
#endif
                             INTERP_z2xi_idx (:,:,i,j), & ! (in)
                             INTERP_z2xi_coef(:,  i,j), & ! (in)
                             Xi(:), Z(:,i,j),           & ! (in)
                             var(:,i,j),                & ! (in)
                             var_Xi(:,i,j)              ) ! (out)
    end do
    end do
    !$acc end kernels

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
    use scale_interp, only: &
       INTERP_interp1d
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Xih  (0:KA)
    real(RP), intent(in)  :: Zh   (0:KA,IA,JA)
    real(RP), intent(in)  :: var  (KA,IA,JA)
    real(RP), intent(out) :: var_Z(KA,IA,JA)

#ifdef _OPENACC
    real(RP) :: workr(KA,7)
    integer  :: worki(KA,2)
#endif

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE) &
    !$omp shared(Xih,Zh,var,var_Z) &
    !$omp shared(INTERP_xih2zh_idx,INTERP_xih2zh_coef)
    !$acc kernels copyin(Xih, Zh, var) copyout(var_Z)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent private(workr,worki)
    do i = IS, IE
       call INTERP_interp1d( KA, KS-1, KE, KA, KS-1, KE, &
#ifdef _OPENACC
                             workr(:,:), worki(:,:), &
#endif
                             INTERP_xih2zh_idx (:,:,i,j), & ! (in)
                             INTERP_xih2zh_coef(:,  i,j), & ! (in)
                             Zh(1:,i,j), Xih(1:),         & ! (in)
                             var(:,i,j),                  & ! (in)
                             var_Z(:,i,j)                 ) ! (out)
    enddo
    enddo
    !$acc end kernels

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
    use scale_interp, only: &
       INTERP_interp1d
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: Zh   (0:KA,IA,JA)
    real(RP), intent(in)  :: Xih  (0:KA)
    real(RP), intent(in)  :: var   (KA,IA,JA)
    real(RP), intent(out) :: var_Xi(KA,IA,JA)

#ifdef _OPENACC
    real(RP) :: workr(KA,7)
    integer  :: worki(KA,2)
#endif

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE) &
    !$omp shared(Zh,Xih,var,var_Xi) &
    !$omp shared(INTERP_zh2xih_idx,INTERP_zh2xih_coef)
    !$acc kernels copyin(Zh, Xih, var) copyout(var_Xi)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent private(workr,worki)
    do i = IS, IE
       call INTERP_interp1d( KA, KS-1, KE, KA, KS-1, KE, &
#ifdef _OPENACC
                             workr(:,:), worki(:,:), &
#endif
                             INTERP_zh2xih_idx (:,:,i,j), & ! (in)
                             INTERP_zh2xih_coef(:,  i,j), & ! (in)
                             Xih(1:), Zh(1:,i,j),         & ! (in)
                             var(:,i,j),                  & ! (in)
                             var_Xi(:,i,j)                ) ! (out)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine INTERP_VERT_zh2xih

  !-----------------------------------------------------------------------------
  !> Setup for pressure coordinate
  subroutine INTERP_VERT_alloc_pres( &
       Kpres, KA, IA, JA )
    implicit none
    integer,  intent(in) :: Kpres
    integer,  intent(in) :: KA
    integer,  intent(in) :: IA
    integer,  intent(in) :: JA
    !---------------------------------------------------------------------------

    allocate( INTERP_xi2p_idx  (Kpres,2,IA,JA) )
    allocate( INTERP_xi2p_coef (Kpres,  IA,JA) )
    allocate( INTERP_xih2p_idx (Kpres,2,IA,JA) )
    allocate( INTERP_xih2p_coef(Kpres,  IA,JA) )
    !$acc enter data create(INTERP_xi2p_idx, INTERP_xi2p_coef, INTERP_xih2p_idx, INTERP_xih2p_coef)

    allocate( LnPRES (KA,IA,JA) )
    allocate( LnPRESh(KA,IA,JA) )
    allocate( LnPaxis(Kpres)    )
    !$acc enter data create(LnPRES, LnPRESh, LnPaxis)

    return
  end subroutine INTERP_VERT_alloc_pres

  !-----------------------------------------------------------------------------
  !> Finalize for pressure coordinate
  subroutine INTERP_VERT_dealloc_pres
    implicit none
    !---------------------------------------------------------------------------

    !$acc exit data delete(INTERP_xi2p_idx, INTERP_xi2p_coef, INTERP_xih2p_idx, INTERP_xih2p_coef)
    deallocate( INTERP_xi2p_idx   )
    deallocate( INTERP_xi2p_coef  )
    deallocate( INTERP_xih2p_idx  )
    deallocate( INTERP_xih2p_coef )

    !$acc exit data delete(LnPRES, LnPRESh, LnPaxis)
    deallocate( LnPRES  )
    deallocate( LnPRESh )
    deallocate( LnPaxis )

    return
  end subroutine INTERP_VERT_dealloc_pres

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
    use scale_interp, only: &
       INTERP_factor1d
    implicit none

    integer,  intent(in)  :: Kpres
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: PRES    (  KA,IA,JA) ! pressure in Xi coordinate [Pa]
    real(RP), intent(in)  :: PRESh   (0:KA,IA,JA) ! pressure in Xi coordinate [Pa], layer interface
    real(RP), intent(in)  :: SFC_PRES(     IA,JA) ! surface pressure          [Pa]
    real(RP), intent(in)  :: Paxis   (Kpres)      ! pressure level to output  [Pa]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data copyin(PRES, PRESh, SFC_PRES, Paxis)

    ! full level

!OCL SERIAL
    !$acc kernels
    do k = 1, Kpres
       LnPaxis(k) = - log( Paxis(k) )
    end do
    !$acc end kernels

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,Kpres,IS,IE,JS,JE) &
    !$omp shared(PRES,LnPRES,LnPaxis) &
    !$omp shared(INTERP_xi2p_idx,INTERP_xi2p_coef)
    !$acc kernels
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent
    do i = IS, IE
       do k = KS, KE
          LnPRES(k,i,j) = - log( PRES(k,i,j) )
       end do
       call INTERP_factor1d( KA, KS, KE, Kpres, 1, Kpres, &
                             LnPRES(:,i,j), LnPaxis(:), & ! (in)
                             INTERP_xi2p_idx (:,:,i,j), & ! (out)
                             INTERP_xi2p_coef(:,  i,j), & ! (out)
                             flag_extrap = .false.      ) ! (in)
    enddo
    enddo
    !$acc end kernels

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,Kpres,IS,IE,JS,JE) &
    !$omp shared(PRESh,LnPRESh,LnPaxis) &
    !$omp shared(INTERP_xih2p_idx,INTERP_xih2p_coef)
    !$acc kernels
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent
    do i = IS, IE
       do k = KS-1, KE
          LnPRESh(k,i,j) = - log( PRESh(k,i,j) )
       end do
       call INTERP_factor1d( KA, KS-1, KE, Kpres, 1, Kpres, &
                             LnPRESh(:,i,j), LnPaxis(:), & ! (in)
                             INTERP_xih2p_idx (:,:,i,j), & ! (out)
                             INTERP_xih2p_coef(:,  i,j), & ! (out)
                             flag_extrap = .false.       ) ! (in)
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine INTERP_VERT_setcoef_pres

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_xi2p( &
       Kpres,      &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       var,        &
       var_P       )
    use scale_interp, only: &
       INTERP_interp1d
    implicit none
    integer,  intent(in)  :: Kpres
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: var  (KA   ,IA,JA)
    real(RP), intent(out) :: var_P(Kpres,IA,JA)

#ifdef _OPENACC
    real(RP) :: workr(KA,7)
    integer  :: worki(KA,2)
#endif

    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,Kpres,IS,IE,JS,JE) &
    !$omp shared(LnPRES,LnPaxis,var,var_P) &
    !$omp shared(INTERP_xi2p_idx,INTERP_xi2p_coef)
    !$acc kernels copyin(var) copyout(var_P)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent private(workr,worki)
    do i = IS, IE
       call INTERP_interp1d( KA, KS, KE, Kpres, 1, Kpres, &
#ifdef _OPENACC
                             workr(:,:), worki(:,:), &
#endif
                             INTERP_xi2p_idx (:,:,i,j), & ! (in)
                             INTERP_xi2p_coef(:,  i,j), & ! (in)
                             LnPRES(:,i,j), LnPaxis(:), & ! (in)
                             var(:,i,j),                & ! (in)
                             var_P(:,i,j)               ) ! (out)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine INTERP_VERT_xi2p

  !-----------------------------------------------------------------------------
  subroutine INTERP_VERT_xih2p( &
       Kpres,      &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       var,        &
       var_P       )
    use scale_interp, only: &
       INTERP_interp1d
    implicit none

    integer,  intent(in)  :: Kpres
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: var  (KA   ,IA,JA)
    real(RP), intent(out) :: var_P(Kpres,IA,JA)

#ifdef _OPENACC
    real(RP) :: workr(KA,7)
    integer  :: worki(KA,2)
#endif

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,Kpres,IS,IE,JS,JE) &
    !$omp shared(LnPRESh,LnPaxis,var,var_P) &
    !$omp shared(INTERP_xih2p_idx,INTERP_xih2p_coef)
    !$acc kernels copyin(var) copyout(var_P)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent private(workr,worki)
    do i = IS, IE
       call INTERP_interp1d( KA, KS-1, KE, Kpres, 1, Kpres, &
#ifdef _OPENACC
                             workr(:,:), worki(:,:), &
#endif
                             INTERP_xih2p_idx (:,:,i,j), & ! (in)
                             INTERP_xih2p_coef(:,  i,j), & ! (in)
                             LnPRESh(:,i,j), LnPaxis(:), & ! (in)
                             var(:,i,j),                 & ! (in)
                             var_P(:,i,j)                ) ! (out)
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine INTERP_VERT_xih2p

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine INTERP_VERT_finalize
    implicit none
    !---------------------------------------------------------------------------

    !$acc exit data delete(INTERP_xi2z_idx, INTERP_xi2z_coef, INTERP_z2xi_idx, INTERP_z2xi_coef)
    deallocate( INTERP_xi2z_idx  )
    deallocate( INTERP_xi2z_coef )
    deallocate( INTERP_z2xi_idx  )
    deallocate( INTERP_z2xi_coef )

    !$acc exit data delete(INTERP_xih2zh_idx, INTERP_xih2zh_coef, INTERP_zh2xih_idx, INTERP_zh2xih_coef)
    deallocate( INTERP_xih2zh_idx  )
    deallocate( INTERP_xih2zh_coef )
    deallocate( INTERP_zh2xih_idx  )
    deallocate( INTERP_zh2xih_coef )

    return
  end subroutine INTERP_VERT_finalize

end module scale_interp_vert
