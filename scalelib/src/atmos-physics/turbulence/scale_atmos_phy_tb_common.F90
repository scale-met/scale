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
  public :: ATMOS_PHY_TB_diffusion_solver
  public :: ATMOS_PHY_TB_calc_tend_MOMZ
  public :: ATMOS_PHY_TB_calc_tend_MOMX
  public :: ATMOS_PHY_TB_calc_tend_MOMY
  public :: ATMOS_PHY_TB_calc_tend_PHI

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
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY, &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ,                 &
       I_XYW,                 &
       I_UYW,                 &
       I_XVW,                 &
       I_UYZ,                 &
       I_XVZ,                 &
       I_UVZ,                 &
       I_XY,                  &
       I_UY,                  &
       I_XV,                  &
       I_UV
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
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY, &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ,                 &
       I_XYW,                 &
       I_UYW,                 &
       I_XVW,                 &
       I_UYZ,                 &
       I_XVZ,                 &
       I_UVZ,                 &
       I_XY,                  &
       I_UY,                  &
       I_XV,                  &
       I_UV
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
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY, &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ,                 &
       I_XYW,                 &
       I_UYW,                 &
       I_XVW,                 &
       I_UYZ,                 &
       I_XVZ,                 &
       I_UVZ,                 &
       I_XY,                  &
       I_UY,                  &
       I_XV,                  &
       I_UV
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
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY, &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ
    use scale_gridtrans, only: &
       I_XYZ,                 &
       I_XYW,                 &
       I_UYW,                 &
       I_XVW,                 &
       I_UYZ,                 &
       I_XVZ,                 &
       I_UVZ,                 &
       I_XY,                  &
       I_UY,                  &
       I_XV,                  &
       I_UV
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

    return
  end subroutine ATMOS_PHY_TB_diffusion_solver

end module scale_atmos_phy_tb_common

