!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Common routines for turbulence
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-12-13 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_tb_common
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
  public :: ATMOS_PHY_TB_calc_strain_tensor
  public :: ATMOS_PHY_TB_diffusion_solver
  public :: ATMOS_PHY_TB_calc_tend_MOMZ
  public :: ATMOS_PHY_TB_calc_tend_MOMX
  public :: ATMOS_PHY_TB_calc_tend_MOMY
  public :: ATMOS_PHY_TB_calc_tend_phi
  public :: ATMOS_PHY_TB_calc_flux_phi

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_calc_strain_tensor( &
       S33_C, S11_C, S22_C, &
       S31_C, S12_C, S23_C, &
       S12_Z, S23_X, S31_Y, &
       S2,                  &
       DENS, MOMZ, MOMX, MOMY, &
       GSQRT, J13G, J23G, J33G, MAPF )
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY,  &
       I_UY,  &
       I_XV,  &
       I_UV
    use scale_grid, only: &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    implicit none

    real(RP), intent(out) :: S33_C (KA,IA,JA) ! (cell center)
    real(RP), intent(out) :: S11_C (KA,IA,JA)
    real(RP), intent(out) :: S22_C (KA,IA,JA)
    real(RP), intent(out) :: S31_C (KA,IA,JA)
    real(RP), intent(out) :: S12_C (KA,IA,JA)
    real(RP), intent(out) :: S23_C (KA,IA,JA)
    real(RP), intent(out) :: S12_Z (KA,IA,JA) ! (z edge or x-y plane)
    real(RP), intent(out) :: S23_X (KA,IA,JA) ! (x edge or y-z plane)
    real(RP), intent(out) :: S31_Y (KA,IA,JA) ! (y edge or z-x plane)
    real(RP), intent(out) :: S2    (KA,IA,JA) ! |S|^2

    real(RP), intent(in)  :: DENS  (KA,IA,JA)
    real(RP), intent(in)  :: MOMZ  (KA,IA,JA)
    real(RP), intent(in)  :: MOMX  (KA,IA,JA)
    real(RP), intent(in)  :: MOMY  (KA,IA,JA)

    real(RP), intent(in)  :: GSQRT (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G  (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G  (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G               !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF  (IA,JA,2,4)  !< map factor

    ! velocity
    real(RP) :: VELZ_C (KA,IA,JA)
    real(RP) :: VELZ_XY(KA,IA,JA)
    real(RP) :: VELX_C (KA,IA,JA)
    real(RP) :: VELX_YZ(KA,IA,JA)
    real(RP) :: VELY_C (KA,IA,JA)
    real(RP) :: VELY_ZX(KA,IA,JA)

    ! work space
    real(RP) :: WORK_V(KA,IA,JA) ! work space (vertex)
    real(RP) :: WORK_Z(KA,IA,JA) !            (z edge or x-y plane)
    real(RP) :: WORK_X(KA,IA,JA) !            (x edge or y-z plane)
    real(RP) :: WORK_Y(KA,IA,JA) !            (y edge or z-x plane)

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j

#ifdef DEBUG
    S33_C (:,:,:) = UNDEF
    S11_C (:,:,:) = UNDEF
    S22_C (:,:,:) = UNDEF
    S31_C (:,:,:) = UNDEF
    S12_C (:,:,:) = UNDEF
    S23_C (:,:,:) = UNDEF
    S12_Z (:,:,:) = UNDEF
    S23_X (:,:,:) = UNDEF
    S31_Y (:,:,:) = UNDEF
    S2    (:,:,:) = UNDEF

    VELZ_C (:,:,:) = UNDEF
    VELZ_XY(:,:,:) = UNDEF
    VELX_C (:,:,:) = UNDEF
    VELX_YZ(:,:,:) = UNDEF
    VELY_C (:,:,:) = UNDEF
    VELY_ZX(:,:,:) = UNDEF

    WORK_V(:,:,:) = UNDEF
    WORK_Z(:,:,:) = UNDEF
    WORK_X(:,:,:) = UNDEF
    WORK_Y(:,:,:) = UNDEF
#endif

   ! momentum -> velocity
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, MOMZ(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELZ_XY(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       VELZ_XY(KE,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMZ(k,i,j) )
       call CHECK( __LINE__, MOMZ(k-1,i,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELZ_C(k,i,j) = 0.5_RP * ( MOMZ(k,i,j) + MOMZ(k-1,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-2, JE+2
    do i = IS-2, IE+2
#ifdef DEBUG
       call CHECK( __LINE__, MOMZ(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i,j) )
#endif
       VELZ_C(KS,i,j) = 0.5_RP * MOMZ(KS,i,j) / DENS(KS,i,j) ! MOMZ(KS-1,i,j) = 0
    enddo

    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    do j = JS-1, JE+1
    do i = IS-2, IE+1
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMX(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELX_YZ(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+1
    do i = IS-2, IE+1
       VELX_YZ(KE+1,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-2, JE+2
    do i = IS-1, IE+2
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMX(k,i,j) )
       call CHECK( __LINE__, MOMX(k,i-1,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELX_C(k,i,j) = 0.5_RP * ( MOMX(k,i,j) + MOMX(k,i-1,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    !$omp parallel do default(none)                   &
    !$omp shared(JS,JE,IS,IE,KS,KE,MOMY,DENS,VELY_ZX) &
    !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMY(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELY_ZX(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-2, JE+1
    do i = IS-1, IE+1
       VELY_ZX(KE+1,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+2
    do i = IS-2, IE+2
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMY(k,i,j) )
       call CHECK( __LINE__, MOMY(k,i,j-1) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELY_C(k,i,j) = 0.5_RP * ( MOMY(k,i,j) + MOMY(k,i,j-1) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! w
       ! (x-y plane; x,y,w)
       ! WORK_Z = VELZ_XY
       ! (y-z plane; u,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i,j) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELZ_C(k,i+1,j) + VELZ_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          WORK_X(KE+1,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane; x,v,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i,j+1) )
       call CHECK( __LINE__, VELZ_C(k,i,j) )
#endif
          WORK_Y(k,i,j) = 0.5_RP * ( VELZ_C(k,i,j+1) + VELZ_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          WORK_Y(KE+1,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dw/dz
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S33_C(k,i,j) = ( VELZ_XY(k,i,j) - VELZ_XY(k-1,i,j) ) * RCDZ(k) &
                       * J33G / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S33_C(KS,i,j) = VELZ_XY(KS,i,j) * RCDZ(KS) & ! VELZ_XY(KS-1,i,j) == 0
                        * J33G / GSQRT(KS,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dx
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i-1,j) )
       call CHECK( __LINE__, GSQRT(k,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, J13G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J13G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S31_C(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i+1,j,I_XYZ)*VELZ_C(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*VELZ_C(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
               + ( J13G(k,i,j,I_XYW)*VELZ_XY(k,i,j) - J13G(k-1,i,j,I_XYW)*VELZ_XY(k-1,i,j) ) * RCDZ(k) &
               ) * MAPF(i,j,1,I_XY)

       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(KS,i+1,j) )
       call CHECK( __LINE__, VELZ_C(KS,i-1,j) )
       call CHECK( __LINE__, GSQRT(KS,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, J13G(KS,i,j,I_XYW) )
       call CHECK( __LINE__, VELZ_C(KE,i+1,j) )
       call CHECK( __LINE__, VELZ_C(KE,i-1,j) )
       call CHECK( __LINE__, GSQRT(KE,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KE,i,j) )
       call CHECK( __LINE__, J13G(KE,i,j,I_XYW) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S31_C(KS,i,j) = 0.5_RP * ( &
                 ( GSQRT(KS,i+1,j,I_XYZ)*VELZ_C(KS,i+1,j) - GSQRT(KS,i-1,j,I_XYZ)*VELZ_C(KS,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
               + ( J13G(KS,i,j,I_XYW)*VELZ_XY(KS,i,j) ) * RCDZ(KS) &
               ) * MAPF(i,j,1,I_XY)
          S31_C(KE,i,j) = 0.5_RP * ( &
                 ( GSQRT(KE,i+1,j,I_XYZ)*VELZ_C(KE,i+1,j) - GSQRT(KE,i-1,j,I_XYZ)*VELZ_C(KE,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
               - ( J13G(KE-1,i,j,I_XYW)*VELZ_XY(KE-1,i,j) ) * RCDZ(KE) &
               ) * MAPF(i,j,1,I_XY)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (y edge, u,y,w)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i+1,j) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S31_Y(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i+1,j,I_XYW)*VELZ_XY(k,i+1,j) - GSQRT(k,i,j,I_XYW)*VELZ_XY(k,i,j) ) * RFDX(i) &
               + ( J13G(k+1,i,j,I_UYZ)*WORK_X(k+1,i,j) - J13G(k,i,j,I_UYZ)*WORK_X (k,i,j)) * RFDZ(k) &
               ) * MAPF(i,j,1,I_UY)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dy
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i,j+1) )
       call CHECK( __LINE__, VELZ_C(k,i,j-1) )
       call CHECK( __LINE__, GSQRT(k,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, J23G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J23G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S23_C(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_XYZ)*VELZ_C(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*VELZ_C(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(k,i,j,I_XYW)*VELZ_XY(k,i,j) - J23G(k-1,i,j,I_XYW)*VELZ_XY(k-1,i,j) ) * RCDZ(k) &
               ) * MAPF(i,j,2,I_XY)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(KS,i,j+1) )
       call CHECK( __LINE__, VELZ_C(KS,i,j-1) )
       call CHECK( __LINE__, GSQRT(KS,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, J23G(KS,i,j,I_XYW) )
       call CHECK( __LINE__, VELZ_C(KE,i,j+1) )
       call CHECK( __LINE__, VELZ_C(KE,i,j-1) )
       call CHECK( __LINE__, GSQRT(KE,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KE,i,j) )
       call CHECK( __LINE__, J23G(KE,i,j,I_XYW) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S23_C(KS,i,j) = 0.5_RP * ( &
                 ( GSQRT(KS,i,j+1,I_XYZ)*VELZ_C(KS,i,j+1) - GSQRT(KS,i,j-1,I_XYZ)*VELZ_C(KS,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(KS,i,j,I_XYW)*VELZ_XY(KS,i,j) ) * RCDZ(KS) &
               ) * MAPF(i,j,2,I_XY)
          S23_C(KE,i,j) = 0.5_RP * ( &
                 ( GSQRT(KE,i,j+1,I_XYZ)*VELZ_C(KE,i,j+1) - GSQRT(KE,i,j-1,I_XYZ)*VELZ_C(KE,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               - ( J23G(KE-1,i,j,I_XYW)*VELZ_XY(KE-1,i,j) ) * RCDZ(KE) &
               ) * MAPF(i,j,2,I_XY)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (x edge; x,v,w)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j+1) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S23_X(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_XYW)*VELZ_XY(k,i,j+1) - GSQRT(k,i,j,I_XYW)*VELZ_XY(k,i,j) ) * RFDY(j) &
               + ( J23G(k+1,i,j,I_XVZ)*WORK_Y(k+1,i,j) - J23G(k,i,j,I_XVZ)*WORK_Y (k,i,j) ) * RFDZ(k) &
               ) * MAPF(i,j,2,I_XV)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! u
       ! (x-y plane; x,y,w)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * ( VELX_C(k+1,i,j) + VELX_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane; u,y,z)
       ! WORK_X = VELX_YZ
       ! (z-x plane; x,v,z)
       do j = JJS-1, JJE
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
#endif
          WORK_Y(k,i,j) = 0.5_RP * ( VELX_C(k,i,j+1) + VELX_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (vertex; u,v,w)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j+1) )
       call CHECK( __LINE__, J23G(k  ,i,j  ,I_UVZ) )
       call CHECK( __LINE__, J23G(k+1,i,j  ,I_UVZ) )
       call CHECK( __LINE__, J23G(k  ,i,j+1,I_UVZ) )
       call CHECK( __LINE__, J23G(k+1,i,j+1,I_UVZ) )
#endif
          WORK_V(k,i,j) = 0.25_RP &
               * ( J23G(k  ,i,j  ,I_UYZ)*VELX_YZ(k  ,i,j  ) &
                 + J23G(k+1,i,j  ,I_UYZ)*VELX_YZ(k+1,i,j  ) &
                 + J23G(k  ,i,j+1,I_UYZ)*VELX_YZ(k  ,i,j+1) &
                 + J23G(k+1,i,j+1,I_UYZ)*VELX_YZ(k+1,i,j+1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! du/dx
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i-1,j) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
       call CHECK( __LINE__, GSQRT(k,i-1,j,I_UYZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J13G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J13G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_C(k,i,j) = ( &
                 ( GSQRT(k,i,j,I_UYZ)*VELX_YZ(k,i,j) - GSQRT(k,i-1,j,I_UYZ)*VELX_YZ(k,i-1,j) ) * RCDX(i) &
               + ( J13G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J13G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) &
               ) * MAPF(i,j,1,I_XY) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(KS,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS,i-1,j) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_UYZ) )
       call CHECK( __LINE__, GSQRT(KS,i-1,j,I_UYZ) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, J13G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE,i-1,j) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_UYZ) )
       call CHECK( __LINE__, GSQRT(KE,i-1,j,I_UYZ) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, J13G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KE-1,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_C(KS,i,j) = ( &
                 ( GSQRT(KS,i,j,I_UYZ)*VELX_YZ(KS,i,j) - GSQRT(KS,i-1,j,I_UYZ)*VELX_YZ(KS,i-1,j) ) * RCDX(i) &
               + ( J13G(KS+1,i,j,I_XYZ)*VELX_C(KS+1,i,j) - J13G(KS,i,j,I_XYZ)*VELX_C(KS,i,j) ) * RFDZ(KS) &
               ) * MAPF(i,j,1,I_XY) / GSQRT(KS,i,j,I_XYZ)
          S11_C(KE,i,j) = ( &
                 ( GSQRT(KE,i,j,I_UYZ)*VELX_YZ(KE,i,j) - GSQRT(KE,i-1,j,I_UYZ)*VELX_YZ(KE,i-1,j) ) * RCDX(i) &
               + ( J13G(KE,i,j,I_XYZ)*VELX_C(KE,i,j) - J13G(KE-1,i,j,I_XYZ)*VELX_C(KE-1,i,j) ) * RFDZ(KE-1) &
               ) * MAPF(i,j,1,I_XY) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dz
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(k,i,j) )
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S31_C(k,i,j) = ( S31_C(k,i,j) & ! dw/dx
               + 0.5_RP * ( VELX_C(k+1,i,j) - VELX_C(k-1,i,j) ) * J33G / ( FDZ(k) + FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(KS,i,j) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S31_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
#endif
          S31_C(KS,i,j) = ( S31_C(KS,i,j) &
               + 0.5_RP * ( VELX_C(KS+1,i,j) - VELX_C(KS,i,j) ) * J33G * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S31_C(KE,i,j) = ( S31_C(KE,i,j) &
               + 0.5_RP * ( VELX_C(KE,i,j) - VELX_C(KE-1,i,j) ) * J33G * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge; u,y,w)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_Y(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S31_Y(k,i,j) = ( S31_Y(k,i,j) & ! dw/dx
               + 0.5_RP * ( VELX_YZ(k+1,i,j) - VELX_YZ(k,i,j) ) * J33G * RFDZ(k) &
               ) / GSQRT(k,i,j,I_UYW)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dy
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j-1) )
       call CHECK( __LINE__, GSQRT(k,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i,j-1,I_XYZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J23G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J23G(k-1,i,j,I_XYW) )
#endif
          S12_C(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_XYZ)*VELX_C(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*VELX_C(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J23G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) &
               ) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(KS,i,j+1) )
       call CHECK( __LINE__, VELX_C(KS,i,j-1) )
       call CHECK( __LINE__, GSQRT(KS,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, J23G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, VELX_C(KE,i,j+1) )
       call CHECK( __LINE__, VELX_C(KE,i,j-1) )
       call CHECK( __LINE__, GSQRT(KE,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, J23G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KE-1,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S12_C(KS,i,j) = 0.5_RP * ( &
                 ( GSQRT(KS,i,j+1,I_XYZ)*VELX_C(KS,i,j+1) - GSQRT(KS,i,j-1,I_XYZ)*VELX_C(KS,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(KS+1,i,j,I_XYZ)*VELX_C(KS+1,i,j) - J23G(KS,i,j,I_XYZ)*VELX_C(KS,i,j) ) * RFDZ(KS) &
               ) * MAPF(i,j,2,I_XY) / GSQRT(KS,i,j,I_XYZ)
          S12_C(KE,i,j) = 0.5_RP * ( &
                 ( GSQRT(KE,i,j+1,I_XYZ)*VELX_C(KE,i,j+1) - GSQRT(KE,i,j-1,I_XYZ)*VELX_C(KE,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(KE,i,j,I_XYZ)*VELX_C(KE,i,j) - J23G(KE-1,i,j,I_XYZ)*VELX_C(KE-1,i,j) ) * RFDZ(KE-1) &
               ) * MAPF(i,j,2,I_XY) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z edge; u,v,z)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S12_Z(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_UYZ)*VELX_YZ(k,i,j+1) - GSQRT(k,i,j,I_UYZ)*VELX_YZ(k,i,j) ) * RFDY(j) &
               + ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k) &
               ) * MAPF(i,j,2,I_UV)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(KS,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(KS,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j+1) )
       call CHECK( __LINE__, J23G(KS+1,i,j,I_UVZ) )
       call CHECK( __LINE__, J23G(KS  ,i,j,I_UVZ) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE-1,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE-1,i,j+1) )
       call CHECK( __LINE__, J23G(KE  ,i,j,I_UVZ) )
       call CHECK( __LINE__, J23G(KE-1,i,j,I_UVZ) )
#endif
          S12_Z(KS,i,j) = 0.25_RP * ( &
                 ( GSQRT(KS,i,j+1,I_UYZ)*VELX_YZ(KS,i,j+1) - GSQRT(KS,i,j,I_UYZ)*VELX_YZ(KS,i,j) ) * RFDY(j) &
               + ( J23G(KS+1,i,j,I_UVZ) * ( VELX_YZ(KS+1,i,j) + VELX_YZ(KS+1,i,j+1) ) &
                 - J23G(KS  ,i,j,I_UVZ) * ( VELX_YZ(KS  ,i,j) + VELX_YZ(KS  ,i,j+1) ) ) * RFDZ(KS) &
               ) * MAPF(i,j,2,I_UV)
          S12_Z(KE,i,j) = 0.25_RP * ( &
                 ( GSQRT(KE,i,j+1,I_UYZ)*VELX_YZ(KE,i,j+1) - GSQRT(KE,i,j,I_UYZ)*VELX_YZ(KE,i,j) ) * RFDY(j) &
               + ( J23G(KE  ,i,j,I_UVZ) * ( VELX_YZ(KE  ,i,j) + VELX_YZ(KE  ,i,j+1) ) &
                 - J23G(KE-1,i,j,I_UVZ) * ( VELX_YZ(KE-1,i,j) + VELX_YZ(KE-1,i,j+1) ) ) * RFDZ(KE-1) &
               ) * MAPF(i,j,2,I_UV)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! v
       ! (x-y plane; x,y,w)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * ( VELY_C(k+1,i,j) + VELY_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane; u,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELZ_C(k,i+1,j) + VELZ_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane; x,v,z)
       ! WORK_Y = VELY_ZX
       ! (vertex; u,v,w)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i+1,j) )
#endif
          WORK_V(k,i,j) = 0.25_RP &
               * ( J13G(k  ,i  ,j,I_XVZ)*VELY_ZX(k  ,i  ,j) &
                 + J13G(k+1,i  ,j,I_XVZ)*VELY_ZX(k+1,i  ,j) &
                 + J13G(k  ,i+1,j,I_XVZ)*VELY_ZX(k  ,i+1,j) &
                 + J13G(k+1,i+1,j,I_XVZ)*VELY_ZX(k+1,i+1,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dv/dy
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j-1) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
       call CHECK( __LINE__, GSQRT(k,i,j-1,I_XVZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J23G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J23G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S22_C(k,i,j) = ( &
                 ( GSQRT(k,i,j,I_XVZ)*VELY_ZX(k,i,j) - GSQRT(k,i,j-1,I_XVZ)*VELY_ZX(k,i,j-1) ) * RCDY(j) &
               + ( J23G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J23G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) &
               ) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j-1) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XVZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j-1,I_XVZ) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, J23G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDY(j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j-1) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XVZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XVZ) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, J23G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KE-1,i,j,I_XYZ) )
#endif
          S22_C(KS,i,j) = ( &
                 ( GSQRT(KS,i,j,I_XVZ)*VELY_ZX(KS,i,j) - GSQRT(KS,i,j-1,I_XVZ)*VELY_ZX(KS,i,j-1) ) * RCDY(j) &
               + ( J23G(KS+1,i,j,I_XYZ)*VELY_C(KS+1,i,j) - J23G(KS,i,j,I_XYZ)*VELY_C(KS,i,j) ) * RFDZ(KS) &
               ) * MAPF(i,j,2,I_XY) / GSQRT(KS,i,j,I_XYZ)
          S22_C(KE,i,j) = ( &
                 ( GSQRT(KE,i,j,I_XVZ)*VELY_ZX(KE,i,j) - GSQRT(KE,i,j-1,I_XVZ)*VELY_ZX(KE,i,j-1) ) * RCDY(j) &
               + ( J23G(KE,i,j,I_XYZ)*VELY_C(KE,i,j) - J23G(KE-1,i,j,I_XYZ)*VELY_C(KE-1,i,j) ) * RFDZ(KE-1) &
               ) * MAPF(i,j,2,I_XY) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dx
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S12_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k,i-1,j) )
       call CHECK( __LINE__, GSQRT(k,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i-1,j,I_XYZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J13G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J13G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_XYZ) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S12_C(k,i,j) = ( S12_C(k,i,j) & ! du/dy
               + 0.5_RP * ( &
                     ( GSQRT(k,i+1,j,I_XYZ)*VELY_C(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*VELY_C(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                   + ( J13G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J13G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) ) * MAPF(i,j,1,I_XY) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S12_C(KS,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i+1,j) )
       call CHECK( __LINE__, VELY_C(KS,i-1,j) )
       call CHECK( __LINE__, GSQRT(KS,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, J13G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, S12_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE,i+1,j) )
       call CHECK( __LINE__, VELY_C(KE,i-1,j) )
       call CHECK( __LINE__, GSQRT(KE,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, J13G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KE-1,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S12_C(KS,i,j) = ( S12_C(KS,i,j) & ! du/dy
               + 0.5_RP * ( &
                     ( GSQRT(KS,i+1,j,I_XYZ)*VELY_C(KS,i+1,j) - GSQRT(KS,i-1,j,I_XYZ)*VELY_C(KS,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                   + ( J13G(KS+1,i,j,I_XYZ)*VELY_C(KS+1,i,j) - J13G(KS,i,j,I_XYZ)*VELY_C(KS,i,j) ) * RFDZ(KS) ) &
                 * MAPF(i,j,1,I_XY) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S12_C(KE,i,j) = ( S12_C(KE,i,j) & ! du/dy
               + 0.5_RP * ( &
                     ( GSQRT(KE,i+1,j,I_XYZ)*VELY_C(KE,i+1,j) - GSQRT(KE,i-1,j,I_XYZ)*VELY_C(KE,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                   + ( J13G(KE,i,j,I_XYZ)*VELY_C(KE,i,j) - J13G(KE-1,i,j,I_XYZ)*VELY_C(KE-1,i,j) ) * RFDZ(KE-1) ) &
                 * MAPF(i,j,1,I_XY) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge; u,v,z)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S12_Z(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_Z(k,i,j) = ( S12_Z(k,i,j) &
               + 0.5_RP * ( &
                     ( GSQRT(k,i+1,j,I_XVZ)*VELY_ZX(k,i+1,j) - GSQRT(k,i,j,I_XVZ)*VELY_ZX(k,i,j) ) * RFDX(i) &
                   + ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k) ) * MAPF(i,j,1,I_UV) &
               ) / GSQRT(k,i,j,I_UVZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S12_Z(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i+1,j) )
       call CHECK( __LINE__, S12_Z(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i+1,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_Z(KS,i,j) = ( S12_Z(KS,i,j) &
               + 0.5_RP * ( &
                     ( GSQRT(KS,i+1,j,I_XVZ)*VELY_ZX(KS,i+1,j) - GSQRT(KS,i,j,I_XVZ)*VELY_ZX(KS,i,j) ) * RFDX(i) &
                   + ( J13G(KS+1,i,j,I_UVZ) * ( VELY_ZX(KS+1,i,j) + VELY_ZX(KS+1,i+1,j) ) &
                     - J13G(KS  ,i,j,I_UVZ) * ( VELY_ZX(KS  ,i,j) + VELY_ZX(KS  ,i+1,j) ) ) * RFDZ(KS) ) * MAPF(i,j,1,I_UV) &
               ) / GSQRT(KS,i,j,I_UVZ)
          S12_Z(KE,i,j) = ( S12_Z(KE,i,j) &
               + 0.5_RP * ( &
                     ( GSQRT(KE,i+1,j,I_XVZ)*VELY_ZX(KE,i+1,j) - GSQRT(KE,i,j,I_XVZ)*VELY_ZX(KE,i,j) ) * RFDX(i) &
                   + ( J13G(KE  ,i,j,I_UVZ) * ( VELY_ZX(KE  ,i,j) + VELY_ZX(KE  ,i+1,j) ) &
                     - J13G(KE-1,i,j,I_UVZ) * ( VELY_ZX(KE-1,i,j) + VELY_ZX(KE-1,i+1,j) ) ) * RFDZ(KE-1) ) * MAPF(i,j,1,I_UV) &
               ) / GSQRT(KE,i,j,I_UVZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dz
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S23_C(k,i,j) = ( S23_C(k,i,j) & ! dw/dy
               + 0.5_RP * ( VELY_C(k+1,i,j) - VELY_C(k-1,i,j) ) * J33G / ( FDZ(k) + FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(KS,i,j) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S23_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
#endif
          S23_C(KS,i,j) = ( S23_C(KS,i,j) &
               + 0.5_RP * ( VELY_C(KS+1,i,j) - VELY_C(KS,i,j) ) * J33G * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S23_C(KE,i,j) = ( S23_C(KE,i,j) &
               + 0.5_RP * ( VELY_C(KE,i,j) - VELY_C(KE-1,i,j) ) * J33G * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (x edge; x,v,w)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_X(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S23_X(k,i,j) = ( S23_X(k,i,j) &
               + 0.5_RP * ( VELY_ZX(k+1,i,j) - VELY_ZX(k,i,j) ) * J33G * RFDZ(k) &
               ) / GSQRT(k,i,j,I_XVW)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


       ! |S|^2 = 2*Sij*Sij
#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (cell center)
       !$omp parallel do default(none)                                            &
       !$omp shared(JJS,JJE,IIS,IIE,KS,KE,S11_C,S22_C,S33_C,S31_C,S12_C,S23_C,S2) &
       !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, S31_C(k,i,j) )
       call CHECK( __LINE__, S12_C(k,i,j) )
       call CHECK( __LINE__, S23_C(k,i,j) )
#endif
          S2(k,i,j) = max( 1e-10_RP, &
                 2.0_RP * ( S11_C(k,i,j)**2 + S22_C(k,i,j)**2 + S33_C(k,i,j)**2 ) &
               + 4.0_RP * ( S31_C(k,i,j)**2 + S12_C(k,i,j)**2 + S23_C(k,i,j)**2 ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_calc_strain_tensor

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_calc_flux_phi( &
       qflx_phi, &
       DENS, PHI, Kh, FACT, &
       GSQRT, J13G, J23G, J33G, MAPF, &
       a, b, c, dt, &
       implicit, &
       IIS, IIE, JJS, JJE )
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY,  &
       I_UY,  &
       I_XV,  &
       I_UV
    use scale_grid, only: &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    implicit none

    real(RP), intent(inout) :: qflx_phi(KA,IA,JA,3)
    real(RP), intent(in)    :: DENS    (KA,IA,JA)
    real(RP), intent(in)    :: PHI     (KA,IA,JA)
    real(RP), intent(in)    :: Kh      (KA,IA,JA)
    real(RP), intent(in)    :: FACT
    real(RP), intent(in)    :: GSQRT   (KA,IA,JA,7)
    real(RP), intent(in)    :: J13G    (KA,IA,JA,7)
    real(RP), intent(in)    :: J23G    (KA,IA,JA,7)
    real(RP), intent(in)    :: J33G
    real(RP), intent(in)    :: MAPF    (IA,JA,2,4)
    real(RP), intent(in)    :: a       (KA,IA,JA)
    real(RP), intent(in)    :: b       (KA,IA,JA)
    real(RP), intent(in)    :: c       (KA,IA,JA)
    real(DP), intent(in)    :: dt
    logical,  intent(in)    :: implicit
    integer,  intent(in)    :: IIS
    integer,  intent(in)    :: IIE
    integer,  intent(in)    :: JJS
    integer,  intent(in)    :: JJE

    real(RP) :: TEND(KA,IA,JA)
    real(RP) :: d(KA)

    integer :: k, i, j

    ! (x-y plane; x,y,w)
    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,DENS,Kh,PHI,qflx_phi,GSQRT,I_XYW,RFDZ,J33G) &
    !$omp shared(FDZ)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, Kh(k,i,j) )
       call CHECK( __LINE__, Kh(k+1,i,j) )
       call CHECK( __LINE__, PHI(k+1,i,j) )
       call CHECK( __LINE__, PHI(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
       qflx_phi(k,i,j,ZDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
               * ( Kh(k,i,j) + Kh(k+1,i,j) ) &
               * ( PHI(k+1,i,j)-PHI(k,i,j) ) * RFDZ(k) * J33G &
               / GSQRT(k,i,j,I_XYW)
    enddo
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JJS, JJE
    do i = IIS, IIE
       qflx_phi(KS-1,i,j,ZDIR) = 0.0_RP
       qflx_phi(KE  ,i,j,ZDIR) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    ! (y-z plane; u,y,z)
    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,DENS,Kh,PHI,qflx_phi,GSQRT,I_XYZ,RFDX,J13G,I_UYZ) &
    !$omp shared(FDZ)
    do j = JJS,   JJE
    do i = IIS-1, IIE
    do k = KS+1,  KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, Kh(k,i,j) )
       call CHECK( __LINE__, Kh(k,i+1,j) )
       call CHECK( __LINE__, PHI(k,i+1,j) )
       call CHECK( __LINE__, PHI(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
       qflx_phi(k,i,j,XDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
               * ( Kh(k,i,j) + Kh(k,i+1,j) ) &
               * ( &
                   ( GSQRT(k,i+1,j,I_XYZ) * PHI(k,i+1,j) &
                   - GSQRT(k,i  ,j,I_XYZ) * PHI(k,i  ,j) ) * RFDX(i) &
                 + ( J13G(k+1,i,j,I_UYZ) * ( PHI(k+1,i+1,j)+PHI(k+1,i,j) ) &
                   - J13G(k-1,i,j,I_UYZ) * ( PHI(k-1,i+1,j)+PHI(k-1,i,j) ) &
                   ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_UYZ)
    enddo
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JJS,   JJE
    do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i+1,j) )
       call CHECK( __LINE__, Kh(KS,i,j) )
       call CHECK( __LINE__, Kh(KS,i+1,j) )
       call CHECK( __LINE__, PHI(KS,i+1,j) )
       call CHECK( __LINE__, PHI(KS,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
       qflx_phi(KS,i,j,XDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(KS,i,j)+DENS(KS,i+1,j) ) &
               * ( Kh(KS,i,j) + Kh(KS,i+1,j) ) &
               * ( &
                   ( GSQRT(KS,i+1,j,I_XYZ) * PHI(KS,i+1,j) &
                   - GSQRT(KS,i  ,j,I_XYZ) * PHI(KS,i  ,j) ) * RFDX(i) &
                 + ( J13G(KS+1,i,j,I_UYZ) * ( PHI(KS+1,i+1,j)+PHI(KS+1,i,j) ) &
                   - J13G(KS  ,i,j,I_UYZ) * ( PHI(KS  ,i+1,j)+PHI(KS  ,i,j) ) &
                   ) * 0.5_RP * RFDZ(KS) &
               ) * MAPF(i,j,1,I_UY) / GSQRT(KS,i,j,I_UYZ)
       qflx_phi(KE,i,j,XDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(KE,i,j)+DENS(KE,i+1,j) ) &
               * ( Kh(KE,i,j) + Kh(KE,i+1,j) ) &
               * ( &
                   ( GSQRT(KE,i+1,j,I_XYZ) * PHI(KE,i+1,j) &
                   - GSQRT(KE,i  ,j,I_XYZ) * PHI(KE,i  ,j) ) * RFDX(i) &
                 + ( J13G(KE  ,i,j,I_UYZ) * ( PHI(KE  ,i+1,j)+PHI(KE  ,i,j) ) &
                   - J13G(KE-1,i,j,I_UYZ) * ( PHI(KE-1,i+1,j)+PHI(KE-1,i,j) ) &
                   ) * 0.5_RP * RFDZ(KE-1) &
               ) * MAPF(i,j,1,I_UY) / GSQRT(KE,i,j,I_UYZ)
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    ! (z-x plane; x,v,z)
    !$omp parallel do default(none)                                                          &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,Kh,PHI,RFDY,DENS,qflx_phi,GSQRT,I_XYZ,J23G,I_XVZ,FDZ) &
    !$omp shared(MAPF,I_XV)                                                                  &
    !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS,   IIE
    do k = KS+1,  KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, Kh(k,i,j) )
       call CHECK( __LINE__, Kh(k,i,j+1) )
       call CHECK( __LINE__, PHI(k,i,j+1) )
       call CHECK( __LINE__, PHI(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
       qflx_phi(k,i,j,YDIR) = - 0.25_RP &
               * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
               * ( Kh(k,i,j) + Kh(k,i,j+1) ) &
               * ( &
                     ( GSQRT(k,i,j+1,I_XYZ) * PHI(k,i,j+1) &
                     - GSQRT(k,i,j  ,I_XYZ) * PHI(k,i,j  ) ) * RFDY(j) &
                   + ( J23G(k+1,i,j,I_XVZ) * ( PHI(k+1,i,j+1)+PHI(k+1,i,j) ) &
                     - J23G(k-1,i,j,I_XVZ) * ( PHI(k-1,i,j+1)+PHI(k-1,i,j) ) &
                     ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
               ) * MAPF(i,j,2,I_XV) / GSQRT(k,i,j,I_XVZ)
    enddo
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JJS-1, JJE
    do i = IIS,   IIE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i,j+1) )
       call CHECK( __LINE__, Kh(KS,i,j) )
       call CHECK( __LINE__, Kh(KS,i,j+1) )
       call CHECK( __LINE__, PHI(KS,i,j+1) )
       call CHECK( __LINE__, PHI(KS,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
       qflx_phi(KS,i,j,YDIR) = - 0.25_RP &
               * ( DENS(KS,i,j)+DENS(KS,i,j+1) ) &
               * ( Kh(KS,i,j) + Kh(KS,i,j+1) ) &
               * ( &
                     ( GSQRT(KS,i,j+1,I_XYZ) * PHI(KS,i,j+1) &
                     - GSQRT(KS,i,j  ,I_XYZ) * PHI(KS,i,j  ) ) * RFDY(j) &
                   + ( J23G(KS+1,i,j,I_XVZ) * ( PHI(KS+1,i,j+1)+PHI(KS+1,i,j) ) &
                     - J23G(KS  ,i,j,I_XVZ) * ( PHI(KS  ,i,j+1)+PHI(KS  ,i,j) ) &
                     ) * 0.5_RP * RFDZ(KS) &
               ) * MAPF(i,j,2,I_XV) / GSQRT(KS,i,j,I_XVZ)
       qflx_phi(KE,i,j,YDIR) = - 0.25_RP &
               * ( DENS(KE,i,j)+DENS(KE,i,j+1) ) &
               * ( Kh(KE,i,j) + Kh(KE,i,j+1) ) &
               * ( &
                     ( GSQRT(KE,i,j+1,I_XYZ) * PHI(KE,i,j+1) &
                     - GSQRT(KE,i,j  ,I_XYZ) * PHI(KE,i,j  ) ) * RFDY(j) &
                   + ( J23G(KE  ,i,j,I_XVZ) * ( PHI(KE  ,i,j+1)+PHI(KE  ,i,j) ) &
                     - J23G(KE-1,i,j,I_XVZ) * ( PHI(KE-1,i,j+1)+PHI(KE-1,i,j) ) &
                     ) * 0.5_RP * RFDZ(KE-1) &
               ) * MAPF(i,j,2,I_XV) / GSQRT(KE,i,j,I_XVZ)
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    if ( implicit ) then
       call ATMOS_PHY_TB_calc_tend_phi( TEND, & ! (out)
                           qflx_phi, & ! (in)
                           GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                           IIS, IIE, JJS, JJE ) ! (in)

       do j = JJS, JJE
       do i = IIS, IIE

          do k = KS, KE
             d(k) = TEND(k,i,j)
          end do

          call ATMOS_PHY_TB_diffusion_solver( &
                  TEND(:,i,j),                & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE                               ) ! (in)

          do k = KS, KE-1
             qflx_phi(k,i,j,ZDIR) = qflx_phi(k,i,j,ZDIR) &
                     - 0.25_RP & ! 1/2/2
                     * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                     * ( Kh(k,i,j) + Kh(k+1,i,j) ) &
                     * dt * ( TEND(k+1,i,j)-TEND(k,i,j) ) * RFDZ(k) * J33G &
                     / GSQRT(k,i,j,I_XYW)
          end do

       end do
       end do

    end if

    return
  end subroutine ATMOS_PHY_TB_calc_flux_phi

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine ATMOS_PHY_TB_diffusion_solver( &
       phi, &
       a, b, c, d, &
       KE_TB )
    implicit none
    real(RP), intent(out) :: phi(KA)
    real(RP), intent(in)  :: a(KA)
    real(RP), intent(in)  :: b(KA)
    real(RP), intent(in)  :: c(KA)
    real(RP), intent(in)  :: d(KA)
    integer,  intent(in)  :: KE_TB
    real(RP) :: e(KA)
    real(RP) :: f(KA)
    real(RP) :: denom
    integer :: k

    e(KS) = - a(KS) / b(KS)
    f(KS) =   d(KS) / b(KS)
    do k = KS+1, KE_TB-1
       denom = b(k) + c(k)*e(k-1)
       e(k) = - a(k) / denom
       f(k) = ( d(k) - c(k)*f(k-1) ) / denom
    end do

    ! flux at the top boundary is zero
    phi(KE_TB) = ( d(KE_TB) - c(KE_TB)*f(KE_TB-1) ) / ( b(KE_TB) + c(KE_TB)*e(KE_TB-1) ) ! = f(KE_PBL)

    do k = KE_TB-1, KS, -1
       phi(k) = e(k) * phi(k+1) + f(k)
    end do
    do k = 1, KS-1
       phi(k) = 0.0_RP
    end do
    do k = KE_TB+1, KA
       phi(k) = 0.0_RP
    end do

    return
  end subroutine ATMOS_PHY_TB_diffusion_solver

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_calc_tend_MOMZ( &
       MOMZ_t_TB, &
       QFLX_MOMZ, &
       GSQRT, J13G, J23G, J33G, MAPF, &
       IIS, IIE, JJS, JJE )
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       CDZ  => GRID_CDZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_XY
    implicit none

    real(RP), intent(out) :: MOMZ_t_TB(KA,IA,JA)
    real(RP), intent(in)  :: QFLX_MOMZ(KA,IA,JA,3)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in)  :: J13G(KA,IA,JA,7)
    real(RP), intent(in)  :: J23G(KA,IA,JA,7)
    real(RP), intent(in)  :: J33G
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)
    integer , intent(in)  :: IIS
    integer , intent(in)  :: IIE
    integer , intent(in)  :: JJS
    integer , intent(in)  :: JJE

    integer :: k, i, j

    !$omp parallel do default(none)                                                          &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,MOMZ_t_TB,GSQRT,I_UYW,I_XVW,QFLX_MOMZ,RCDX,MAPF,I_XY) &
    !$omp shared(RCDY,J13G,I_XYZ,J23G,CDZ,J33G,RFDZ,I_XYW)                                   &
    !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS+1, KE-2
       MOMZ_t_TB(k,i,j) = &
            - ( ( GSQRT(k,i  ,j,I_UYW) * QFLX_MOMZ(k,i  ,j,XDIR) &
                - GSQRT(k,i-1,j,I_UYW) * QFLX_MOMZ(k,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XY) &
              + ( GSQRT(k,i,j  ,I_XVW) * QFLX_MOMZ(k,i,j  ,YDIR) &
                - GSQRT(k,i,j-1,I_XVW) * QFLX_MOMZ(k,i,j-1,YDIR) ) * RCDY(j) &
              + ( ( J13G (k+1,i,j,I_XYZ) * ( QFLX_MOMZ(k+1,i,j,XDIR) + QFLX_MOMZ(k+1,i-1,j,XDIR) ) &
                  - J13G (k-1,i,j,I_XYZ) * ( QFLX_MOMZ(k-1,i,j,XDIR) + QFLX_MOMZ(k-1,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XY) &
                + ( J23G (k+1,i,j,I_XYZ) * ( QFLX_MOMZ(k+1,i,j,YDIR) + QFLX_MOMZ(k+1,i,j-1,YDIR) ) &
                  - J23G (k-1,i,j,I_XYZ) * ( QFLX_MOMZ(k-1,i,j,YDIR) + QFLX_MOMZ(k-1,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_XY) &
                ) * 0.5_RP / ( CDZ(k+1)+CDZ(k) ) &
              + J33G * ( QFLX_MOMZ(k+1,i,j,ZDIR) - QFLX_MOMZ(k,i,j,ZDIR) ) * RFDZ(k) ) &
              / GSQRT(k,i,j,I_XYW)
    enddo
    enddo
    enddo

    do j = JJS, JJE
    do i = IIS, IIE
       MOMZ_t_TB(KS,i,j) = &
            - ( ( GSQRT(KS,i  ,j,I_UYW) * QFLX_MOMZ(KS,i  ,j,XDIR) &
                - GSQRT(KS,i-1,j,I_UYW) * QFLX_MOMZ(KS,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XY) &
              + ( GSQRT(KS,i,j  ,I_XVW) * QFLX_MOMZ(KS,i,j  ,YDIR) &
                - GSQRT(KS,i,j-1,I_XVW) * QFLX_MOMZ(KS,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_XY) &
              + ( ( J13G (KS+1,i,j,I_XYZ) * ( QFLX_MOMZ(KS+1,i,j,XDIR) + QFLX_MOMZ(KS+1,i-1,j,XDIR) ) &
                  - J13G (KS  ,i,j,I_XYZ) * ( QFLX_MOMZ(KS  ,i,j,XDIR) + QFLX_MOMZ(KS  ,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XY) &
                + ( J23G (KS+1,i,j,I_XYZ) * ( QFLX_MOMZ(KS+1,i,j,YDIR) + QFLX_MOMZ(KS+1,i,j-1,YDIR) ) &
                  - J23G (KS  ,i,j,I_XYZ) * ( QFLX_MOMZ(KS  ,i,j,YDIR) + QFLX_MOMZ(KS  ,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_XY) &
                ) * 0.5_RP * RCDZ(KS+1) &
              + J33G * ( QFLX_MOMZ(KS+1,i,j,ZDIR) - QFLX_MOMZ(KS,i,j,ZDIR) ) * RFDZ(KS) ) &
            / GSQRT(KS,i,j,I_XYW)
       MOMZ_t_TB(KE-1,i,j) = &
            - ( ( GSQRT(KE-1,i  ,j,I_UYW) * QFLX_MOMZ(KE-1,i  ,j,XDIR) &
                - GSQRT(KE-1,i-1,j,I_UYW) * QFLX_MOMZ(KE-1,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XY) &
              + ( GSQRT(KE-1,i,j  ,I_XVW) * QFLX_MOMZ(KE-1,i,j  ,YDIR) &
                - GSQRT(KE-1,i,j-1,I_XVW) * QFLX_MOMZ(KE-1,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_XY) &
              + ( ( J13G (KE-1,i,j,I_XYZ) * ( QFLX_MOMZ(KE-1,i,j,XDIR) + QFLX_MOMZ(KE-1,i-1,j,XDIR) ) &
                  - J13G (KE-2,i,j,I_XYZ) * ( QFLX_MOMZ(KE-2,i,j,XDIR) + QFLX_MOMZ(KE-2,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XY) &
                + ( J23G (KE-1,i,j,I_XYZ) * ( QFLX_MOMZ(KE-1,i,j,YDIR) + QFLX_MOMZ(KE-1,i,j-1,YDIR) ) &
                - J23G (KE-2,i,j,I_XYZ) * ( QFLX_MOMZ(KE-2,i,j,YDIR) + QFLX_MOMZ(KE-2,i,j-1,YDIR) ) &
                ) * MAPF(i,j,2,I_XY) &
              ) * 0.5_RP * RCDZ(KE-1) &
              + J33G * ( QFLX_MOMZ(KE,i,j,ZDIR) - QFLX_MOMZ(KE-1,i,j,ZDIR) ) * RFDZ(KE-1) ) &
            / GSQRT(KE-1,i,j,I_XYW)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_calc_tend_MOMZ

  subroutine ATMOS_PHY_TB_calc_tend_MOMX( &
       MOMX_t_TB, &
       QFLX_MOMX, &
       GSQRT, J13G, J23G, J33G, MAPF, &
       IIS, IIE, JJS, JJE )
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_UYW, &
       I_UYZ, &
       I_UVZ, &
       I_UY
    implicit none
    real(RP), intent(out) :: MOMX_t_TB(KA,IA,JA)
    real(RP), intent(in)  :: QFLX_MOMX(KA,IA,JA,3)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)
    real(RP), intent(in)  :: J13G(KA,IA,JA,7)
    real(RP), intent(in)  :: J23G(KA,IA,JA,7)
    real(RP), intent(in)  :: J33G
    integer , intent(in)  :: IIS
    integer , intent(in)  :: IIE
    integer , intent(in)  :: JJS
    integer , intent(in)  :: JJE

    integer :: k, i, j

    !$omp parallel do default(none)                                                          &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,MOMX_t_TB,GSQRT,I_UVZ,I_XYZ,QFLX_MOMX,RFDX,MAPF,I_UY) &
    !$omp shared(RCDY,J13G,I_UYW,J23G,FDZ,J33G,RCDZ,I_UYZ)                                   &
    !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS+1, KE-1
       MOMX_t_TB(k,i,j) = &
            - ( ( GSQRT(k,i+1,j,I_XYZ) * QFLX_MOMX(k,i+1,j,XDIR) &
                - GSQRT(k,i  ,j,I_XYZ) * QFLX_MOMX(k,i  ,j,XDIR) ) * RFDX(i) * MAPF(i,j,1,I_UY) &
              + ( GSQRT(k,i,j  ,I_UVZ) * QFLX_MOMX(k,i,j  ,YDIR) &
                - GSQRT(k,i,j-1,I_UVZ) * QFLX_MOMX(k,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_UY) &
              + ( ( J13G (k+1,i,j,I_UYW) * ( QFLX_MOMX(k+1,i+1,j,XDIR) + QFLX_MOMX(k+1,i,j  ,XDIR) ) &
                  - J13G (k-1,i,j,I_UYW) * ( QFLX_MOMX(k-1,i+1,j,XDIR) + QFLX_MOMX(k-1,i,j  ,XDIR) ) &
                  ) * MAPF(i,j,1,I_UY) &
                + ( J23G (k+1,i,j,I_UYW) * ( QFLX_MOMX(k+1,i  ,j,YDIR) + QFLX_MOMX(k+1,i,j-1,YDIR) ) &
                  - J23G (k-1,i,j,I_UYW) * ( QFLX_MOMX(k-1,i  ,j,YDIR) + QFLX_MOMX(k-1,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_UY) &
                ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
              + J33G * ( QFLX_MOMX(k,i,j,ZDIR) - QFLX_MOMX(k-1,i,j,ZDIR) ) * RCDZ(k) ) &
            / GSQRT(k,i,j,I_UYZ)
    enddo
    enddo
    enddo
    do j = JJS, JJE
    do i = IIS, IIE
       MOMX_t_TB(KS,i,j) = &
            - ( ( GSQRT(KS,i+1,j,I_XYZ) * QFLX_MOMX(KS,i+1,j,XDIR) &
                - GSQRT(KS,i  ,j,I_XYZ) * QFLX_MOMX(KS,i  ,j,XDIR) ) * RFDX(i) * MAPF(i,j,1,I_UY) &
              + ( GSQRT(KS,i,j  ,I_UVZ) * QFLX_MOMX(KS,i,j  ,YDIR) &
                - GSQRT(KS,i,j-1,I_UVZ) * QFLX_MOMX(KS,i,j-1,YDIR) ) * RCDY(j) &
              + ( ( J13G (KS+1,i,j,I_UYW) * ( QFLX_MOMX(KS+1,i+1,j,XDIR) + QFLX_MOMX(KS+1,i,j  ,XDIR) ) &
                  - J13G (KS  ,i,j,I_UYW) * ( QFLX_MOMX(KS  ,i+1,j,XDIR) + QFLX_MOMX(KS  ,i,j  ,XDIR) ) &
                  ) * MAPF(i,j,1,I_UY) &
              + ( J23G (KS+1,i,j,I_UYW) * ( QFLX_MOMX(KS+1,i  ,j,YDIR) + QFLX_MOMX(KS+1,i,j-1,YDIR) ) &
                - J23G (KS  ,i,j,I_UYW) * ( QFLX_MOMX(KS  ,i  ,j,YDIR) + QFLX_MOMX(KS  ,i,j-1,YDIR) ) &
                ) * MAPF(i,j,2,I_UY) &
              ) * 0.5_RP * RCDZ(KS) &
            + J33G * ( QFLX_MOMX(KS,i,j,ZDIR) ) * RFDZ(KS) ) &
          / GSQRT(KS,i,j,I_UYZ)
       MOMX_t_TB(KE,i,j) = &
            - ( ( GSQRT(KE,i+1,j,I_XYZ) * QFLX_MOMX(KE,i+1,j,XDIR) &
                - GSQRT(KE,i  ,j,I_XYZ) * QFLX_MOMX(KE,i  ,j,XDIR) ) * RFDX(i) * MAPF(i,j,1,I_UY) &
              + ( GSQRT(KE,i,j  ,I_UVZ) * QFLX_MOMX(KE,i,j  ,YDIR) &
                - GSQRT(KE,i,j-1,I_UVZ) * QFLX_MOMX(KE,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_UY)&
              + ( ( J13G (KE  ,i,j,I_UYW) * ( QFLX_MOMX(KE  ,i+1,j,XDIR) + QFLX_MOMX(KE    ,i,j  ,XDIR) ) &
                  - J13G (KE-1,i,j,I_UYW) * ( QFLX_MOMX(KE-1,i+1,j,XDIR) + QFLX_MOMX(KE-1  ,i,j  ,XDIR) ) &
                  ) * MAPF(i,j,1,I_UY) &
                + ( J23G (KE  ,i,j,I_UYW) * ( QFLX_MOMX(KE  ,i  ,j,YDIR) + QFLX_MOMX(KE-1+1,i,j-1,YDIR) ) &
                  - J23G (KE-1,i,j,I_UYW) * ( QFLX_MOMX(KE-1,i  ,j,YDIR) + QFLX_MOMX(KE-1  ,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_UY) &
                ) * 0.5_RP * RFDZ(KE-1) &
              - J33G * ( QFLX_MOMX(KE-1,i,j,ZDIR) ) * RCDZ(KE) ) &
            / GSQRT(KE,i,j,I_UYZ)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_calc_tend_MOMX

  subroutine ATMOS_PHY_TB_calc_tend_MOMY( &
       MOMY_t_TB, &
       QFLX_MOMY, &
       GSQRT, J13G, J23G, J33G, MAPF, &
       IIS, IIE, JJS, JJE )
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RFDZ => GRID_RFDZ, &
       RFDY => GRID_RFDY, &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XVW, &
       I_UVZ, &
       I_XV
    implicit none

    real(RP), intent(out) :: MOMY_t_TB(KA,IA,JA)
    real(RP), intent(in)  :: QFLX_MOMY(KA,IA,JA,3)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)
    real(RP), intent(in)  :: J13G(KA,IA,JA,7)
    real(RP), intent(in)  :: J23G(KA,IA,JA,7)
    real(RP), intent(in)  :: J33G
    integer , intent(in)  :: IIS
    integer , intent(in)  :: IIE
    integer , intent(in)  :: JJS
    integer , intent(in)  :: JJE

    integer :: k, i, j

    !$omp parallel do default(none)                                                          &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,MOMY_t_TB,GSQRT,I_UVZ,I_XYZ,QFLX_MOMY,RCDX,MAPF,I_XV) &
    !$omp shared(RFDY,J13G,I_XVW,J23G,FDZ,J33G,RCDZ)                                         &
    !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS+1, KE-1
       MOMY_t_TB(k,i,j) = &
            - ( ( GSQRT(k,i  ,j  ,I_UVZ) * QFLX_MOMY(k,i  ,j,XDIR) &
                - GSQRT(k,i-1,j  ,I_UVZ) * QFLX_MOMY(k,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XV) &
              + ( GSQRT(k,i  ,j+1,I_XYZ) * QFLX_MOMY(k,i,j+1,YDIR) &
                - GSQRT(k,i  ,j  ,I_XYZ) * QFLX_MOMY(k,i,j  ,YDIR) ) * RFDY(j) * MAPF(i,j,2,I_XV) &
              + ( ( J13G (k+1,i,j  ,I_XVW) * ( QFLX_MOMY(k+1,i,j  ,XDIR) + QFLX_MOMY(k+1,i-1,j,XDIR) ) &
                  - J13G (k-1,i,j  ,I_XVW) * ( QFLX_MOMY(k-1,i,j  ,XDIR) + QFLX_MOMY(k-1,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XV) &
                + ( J23G (k+1,i,j+1,I_XVW) * ( QFLX_MOMY(k+1,i,j+1,YDIR) + QFLX_MOMY(k+1,i  ,j,YDIR) ) &
                  - J23G (k-1,i,j+1,I_XVW) * ( QFLX_MOMY(k-1,i,j+1,YDIR) + QFLX_MOMY(k-1,i  ,j,YDIR) ) &
                  ) * MAPF(i,j,2,I_XV) &
                ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
              + J33G * ( QFLX_MOMY(k,i,j,ZDIR) - QFLX_MOMY(k-1,i,j,ZDIR) ) * RCDZ(k) ) &
            / GSQRT(k,i,j,I_XVW)
    enddo
    enddo
    enddo
    do j = JJS, JJE
    do i = IIS, IIE
       MOMY_t_TB(KS,i,j) = &
            - ( ( GSQRT(KS,i  ,j  ,I_UVZ) * QFLX_MOMY(KS,i  ,j,XDIR) &
                - GSQRT(KS,i-1,j  ,I_UVZ) * QFLX_MOMY(KS,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XV) &
              + ( GSQRT(KS,i  ,j+1,I_XYZ) * QFLX_MOMY(KS,i,j+1,YDIR) &
                - GSQRT(KS,i  ,j  ,I_XYZ) * QFLX_MOMY(KS,i,j  ,YDIR) ) * RFDY(j) * MAPF(i,j,2,I_XV) &
              + ( ( J13G (KS+1,i,j  ,I_XVW) * ( QFLX_MOMY(KS+1,i,j  ,XDIR) + QFLX_MOMY(KS+1,i-1,j,XDIR) ) &
                  - J13G (KS  ,i,j  ,I_XVW) * ( QFLX_MOMY(KS  ,i,j  ,XDIR) + QFLX_MOMY(KS  ,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XV) &
                + ( J23G (KS+1,i,j+1,I_XVW) * ( QFLX_MOMY(KS+1,i,j+1,YDIR) + QFLX_MOMY(KS+1,i  ,j,YDIR) ) &
                  - J23G (KS  ,i,j+1,I_XVW) * ( QFLX_MOMY(KS  ,i,j+1,YDIR) + QFLX_MOMY(KS  ,i  ,j,YDIR) ) &
                  ) * MAPF(i,j,2,I_XV) &
                ) * 0.5_RP * RFDZ(KS) &
              + J33G * ( QFLX_MOMY(KS,i,j,ZDIR) ) * RCDZ(KS) ) &
            / GSQRT(KS,i,j,I_XVW)
       MOMY_t_TB(KE,i,j) = &
            - ( ( GSQRT(KE,i  ,j  ,I_UVZ) * QFLX_MOMY(KE,i  ,j,XDIR) &
                - GSQRT(KE,i-1,j  ,I_UVZ) * QFLX_MOMY(KE,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XV) &
              + ( GSQRT(KE,i  ,j+1,I_XYZ) * QFLX_MOMY(KE,i,j+1,YDIR) &
                - GSQRT(KE,i  ,j  ,I_XYZ) * QFLX_MOMY(KE,i,j  ,YDIR) ) * RFDY(j) * MAPF(i,j,2,I_XV) &
              + ( ( J13G (KE  ,i,j  ,I_XVW) * ( QFLX_MOMY(KE  ,i,j  ,XDIR) + QFLX_MOMY(KE  ,i-1,j,XDIR) ) &
                  - J13G (KE-1,i,j  ,I_XVW) * ( QFLX_MOMY(KE-1,i,j  ,XDIR) + QFLX_MOMY(KE-1,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XV) &
                + ( J23G (KE  ,i,j+1,I_XVW) * ( QFLX_MOMY(KE  ,i,j+1,YDIR) + QFLX_MOMY(KE  ,i  ,j,YDIR) ) &
                  - J23G (KE-1,i,j+1,I_XVW) * ( QFLX_MOMY(KE-1,i,j+1,YDIR) + QFLX_MOMY(KE-1,i  ,j,YDIR) ) &
                  ) * MAPF(i,j,2,I_XV) &
                ) * 0.5_RP * RFDZ(KE-1) &
              - J33G * ( QFLX_MOMY(KE-1,i,j,ZDIR) ) * RCDZ(KE) ) &
            / GSQRT(KE,i,j,I_XVW)
    end do
    end do

    return
  end subroutine ATMOS_PHY_TB_calc_tend_MOMY

  subroutine ATMOS_PHY_TB_calc_tend_phi( &
       phi_t_TB, &
       QFLX_phi, &
       GSQRT, J13G, J23G, J33G, MAPF, &
       IIS, IIE, JJS, JJE )
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY
    implicit none

    real(RP), intent(out) :: phi_t_TB(KA,IA,JA)
    real(RP), intent(in)  :: QFLX_phi(KA,IA,JA,3)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)
    real(RP), intent(in)  :: J13G(KA,IA,JA,7)
    real(RP), intent(in)  :: J23G(KA,IA,JA,7)
    real(RP), intent(in)  :: J33G
    integer , intent(in)  :: IIS
    integer , intent(in)  :: IIE
    integer , intent(in)  :: JJS
    integer , intent(in)  :: JJE

    integer :: k, i, j

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JJS,JJE,IIS,IIE,KS,KE,phi_t_TB,GSQRT,I_UYZ,QFLX_phi,I_UVZ,RCDX,MAPF,I_XY) &
    !$omp shared(I_XVZ,J13G,I_XYW,J23G,FDZ,J33G,RCDZ,I_XYZ,RCDY)
    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS+1, KE-1
       phi_t_TB(k,i,j) = &
            - ( ( GSQRT(k,i  ,j,I_UYZ) * QFLX_phi(k,i  ,j,XDIR) &
                - GSQRT(k,i-1,j,I_UVZ) * QFLX_phi(k,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XY) &
              + ( GSQRT(k,i,j  ,I_XVZ) * QFLX_phi(k,i,j  ,YDIR) &
                - GSQRT(k,i,j-1,I_XVZ) * QFLX_phi(k,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_XY) &
              + ( ( J13G(k+1,i,j,I_XYW) * ( QFLX_phi(k+1,i,j,XDIR) + QFLX_phi(k+1,i-1,j,XDIR) ) &
                  - J13G(k-1,i,j,I_XYW) * ( QFLX_phi(k-1,i,j,XDIR) + QFLX_phi(k-1,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XY) &
                + ( J23G(k+1,i,j,I_XYW) * ( QFLX_phi(k+1,i,j,YDIR) + QFLX_phi(k+1,i,j-1,YDIR) ) &
                  - J23G(k-1,i,j,I_XYW) * ( QFLX_phi(k-1,i,j,YDIR) + QFLX_phi(k-1,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_XY) &
                ) * 0.5_RP / ( FDZ(k) + FDZ(k-1) ) &
              + J33G * ( QFLX_phi(k,i,j,ZDIR) - QFLX_phi(k-1,i,j,ZDIR) ) * RCDZ(k) ) &
            / GSQRT(k,i,j,I_XYZ)
    enddo
    enddo
    enddo
    do j = JJS, JJE
    do i = IIS, IIE
       phi_t_TB(KS,i,j) = &
            - ( ( GSQRT(KS,i  ,j,I_UYZ) * QFLX_phi(KS,i  ,j,XDIR) &
                - GSQRT(KS,i-1,j,I_UVZ) * QFLX_phi(KS,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XY) &
              + ( GSQRT(KS,i,j  ,I_XVZ) * QFLX_phi(KS,i,j  ,YDIR) &
                - GSQRT(KS,i,j-1,I_XVZ) * QFLX_phi(KS,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_XY) &
              + ( ( J13G(KS+1,i,j,I_XYW) * ( QFLX_phi(KS+1,i,j,XDIR) + QFLX_phi(KS+1,i-1,j,XDIR) ) &
                  - J13G(KS  ,i,j,I_XYW) * ( QFLX_phi(KS  ,i,j,XDIR) + QFLX_phi(KS  ,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XY) &
                + ( J23G(KS+1,i,j,I_XYW) * ( QFLX_phi(KS+1,i,j,YDIR) + QFLX_phi(KS+1,i,j-1,YDIR) ) &
                  - J23G(KS  ,i,j,I_XYW) * ( QFLX_phi(KS  ,i,j,YDIR) + QFLX_phi(KS  ,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_XY) &
                ) * 0.5_RP * RFDZ(KS) &
              + J33G * ( QFLX_phi(KS,i,j,ZDIR) ) * RCDZ(KS) ) &
            / GSQRT(KS,i,j,I_XYZ)
       phi_t_TB(KE,i,j) = &
            - ( ( GSQRT(KE,i  ,j,I_UYZ) * QFLX_phi(KE,i  ,j,XDIR) &
                - GSQRT(KE,i-1,j,I_UVZ) * QFLX_phi(KE,i-1,j,XDIR) ) * RCDX(i) * MAPF(i,j,1,I_XY) &
              + ( GSQRT(KE,i,j  ,I_XVZ) * QFLX_phi(KE,i,j  ,YDIR) &
                - GSQRT(KE,i,j-1,I_XVZ) * QFLX_phi(KE,i,j-1,YDIR) ) * RCDY(j) * MAPF(i,j,2,I_XY) &
              + ( ( J13G(KE  ,i,j,I_XYW) * ( QFLX_phi(KE  ,i,j,XDIR) + QFLX_phi(KE  ,i-1,j,XDIR) ) &
                  - J13G(KE-1,i,j,I_XYW) * ( QFLX_phi(KE-1,i,j,XDIR) + QFLX_phi(KE-1,i-1,j,XDIR) ) &
                  ) * MAPF(i,j,1,I_XY) &
                + ( J23G(KE  ,i,j,I_XYW) * ( QFLX_phi(KE  ,i,j,YDIR) + QFLX_phi(KE  ,i,j-1,YDIR) ) &
                  - J23G(KE-1,i,j,I_XYW) * ( QFLX_phi(KE-1,i,j,YDIR) + QFLX_phi(KE-1,i,j-1,YDIR) ) &
                  ) * MAPF(i,j,2,I_XY) &
                ) * 0.5_RP * RFDZ(KE-1) &
              - J33G * ( QFLX_phi(KE-1,i,j,ZDIR) ) * RCDZ(KE) ) &
            / GSQRT(KE,i,j,I_XYZ)
    end do
    end do

    return
  end subroutine ATMOS_PHY_TB_calc_tend_phi

end module scale_atmos_phy_tb_common

