!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)       [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_tb_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_TB_driver_setup
  public :: ATMOS_PHY_TB_driver

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
  !> Setup
  subroutine ATMOS_PHY_TB_driver_setup
    use scale_grid, only: &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY, &
       CZ  => GRID_CZ
    use mod_atmos_admin, only: &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_sw_phy_tb
    use scale_atmos_phy_tb, only: &
       ATMOS_PHY_TB_setup
    use mod_atmos_phy_tb_vars, only: &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_TB] / Origin[SCALE-LES]'

    if ( ATMOS_sw_phy_tb ) then

       ! setup library component
       call ATMOS_PHY_TB_setup( ATMOS_PHY_TB_TYPE, & ! [IN]
                                CDZ, CDX, CDY, CZ  ) ! [IN]

       ! run once (only for the diagnostic value)
       call ATMOS_PHY_TB_driver( .true., .false. )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    do j = JS, JE
    do i = IS, IE
       MOMZ_t_TB(KS-1,i,j) = 0.0_RP
       MOMZ_t_TB(KE  ,i,j) = 0.0_RP
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_driver_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_driver( update_flag, history_flag )
  !> Driver
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
       GSQRT => GTRANS_GSQRT, &
       J13G  => GTRANS_J13G,  &
       J23G  => GTRANS_J23G,  &
       J33G  => GTRANS_J33G,  &
       I_XYZ,                 &
       I_XYW,                 &
       I_UYW,                 &
       I_XVW,                 &
       I_UYZ,                 &
       I_XVZ,                 &
       I_UVZ
    use scale_time, only: &
       dt_TB => TIME_DTSEC_ATMOS_PHY_TB
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_tb, only: &
       ATMOS_PHY_TB
    use mod_atmos_vars, only: &
       DENS_av,           &
       MOMZ_av,           &
       MOMX_av,           &
       MOMY_av,           &
       RHOT_av,           &
       QTRC_av,           &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       QTRC_t => QTRC_tp
    use mod_atmos_phy_tb_vars, only: &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t, &
       MOMX_t_TB => ATMOS_PHY_TB_MOMX_t, &
       MOMY_t_TB => ATMOS_PHY_TB_MOMY_t, &
       RHOT_t_TB => ATMOS_PHY_TB_RHOT_t, &
       QTRC_t_TB => ATMOS_PHY_TB_QTRC_t, &
       TKE       => ATMOS_PHY_TB_TKE,    &
       NU        => ATMOS_PHY_TB_NU
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in) :: history_flag

    ! eddy viscosity/diffusion flux
    real(RP) :: QFLX_MOMZ(KA,IA,JA,3)
    real(RP) :: QFLX_MOMX(KA,IA,JA,3)
    real(RP) :: QFLX_MOMY(KA,IA,JA,3)
    real(RP) :: QFLX_RHOT(KA,IA,JA,3)
    real(RP) :: QFLX_QTRC(KA,IA,JA,QA,3)

    real(RP) :: Ri(KA,IA,JA)
    real(RP) :: Pr(KA,IA,JA)

    integer :: JJS, JJE
    integer :: IIS, IIE

    real(RP) :: RHOQ(KA,IA,JA)
    real(RP) :: total ! dummy

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call ATMOS_PHY_TB( QFLX_MOMZ, & ! [OUT]
                          QFLX_MOMX, & ! [OUT]
                          QFLX_MOMY, & ! [OUT]
                          QFLX_RHOT, & ! [OUT]
                          QFLX_QTRC, & ! [OUT]
                          TKE,       & ! [OUT]
                          NU,        & ! [OUT]
                          Ri,        & ! [OUT]
                          Pr,        & ! [OUT]
                          MOMZ_av,   & ! [IN]
                          MOMX_av,   & ! [IN]
                          MOMY_av,   & ! [IN]
                          RHOT_av,   & ! [IN]
                          DENS_av,   & ! [IN]
                          QTRC_av,   & ! [IN]
                          GSQRT,     & ! [IN]
                          J13G,      & ! [IN]
                          J23G,      & ! [IN]
                          J33G       ) ! [IN]

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-2
             MOMZ_t_TB(k,i,j) = - ( ( GSQRT(k,i  ,j,I_UYW) * QFLX_MOMZ(k,i  ,j,XDIR) &
                                    - GSQRT(k,i-1,j,I_UYW) * QFLX_MOMZ(k,i-1,j,XDIR) ) * RCDX(i) &
                                  + ( GSQRT(k,i,j  ,I_XVW) * QFLX_MOMZ(k,i,j  ,YDIR) &
                                    - GSQRT(k,i,j-1,I_XVW) * QFLX_MOMZ(k,i,j-1,YDIR) ) * RCDY(j) &
                                  + ( J13G (k+1,i,j,I_XYZ) * ( QFLX_MOMZ(k+1,i,j,XDIR) + QFLX_MOMZ(k+1,i-1,j,XDIR) ) &
                                    - J13G (k-1,i,j,I_XYZ) * ( QFLX_MOMZ(k-1,i,j,XDIR) + QFLX_MOMZ(k-1,i-1,j,XDIR) ) &
                                    + J23G (k+1,i,j,I_XYZ) * ( QFLX_MOMZ(k+1,i,j,YDIR) + QFLX_MOMZ(k+1,i,j-1,YDIR) ) &
                                    - J23G (k-1,i,j,I_XYZ) * ( QFLX_MOMZ(k-1,i,j,YDIR) + QFLX_MOMZ(k-1,i,j-1,YDIR) ) &
                                    ) * 0.5_RP / ( CDZ(k+1)+CDZ(k) ) &
                                  + J33G * ( QFLX_MOMZ(k+1,i,j,ZDIR) - QFLX_MOMZ(k,i,j,ZDIR) ) * RFDZ(k) &
                                  ) / GSQRT(k,i,j,I_XYW)
          enddo
          enddo
          enddo

          do j = JJS, JJE
          do i = IIS, IIE
             MOMZ_t_TB(KS,i,j) = - ( ( GSQRT(KS,i  ,j,I_UYW) * QFLX_MOMZ(KS,i  ,j,XDIR) &
                                     - GSQRT(KS,i-1,j,I_UYW) * QFLX_MOMZ(KS,i-1,j,XDIR) ) * RCDX(i) &
                                   + ( GSQRT(KS,i,j  ,I_XVW) * QFLX_MOMZ(KS,i,j  ,YDIR) &
                                     - GSQRT(KS,i,j-1,I_XVW) * QFLX_MOMZ(KS,i,j-1,YDIR) ) * RCDY(j) &
                                   + ( J13G (KS+1,i,j,I_XYZ) * ( QFLX_MOMZ(KS+1,i,j,XDIR) + QFLX_MOMZ(KS+1,i-1,j,XDIR) ) &
                                     - J13G (KS  ,i,j,I_XYZ) * ( QFLX_MOMZ(KS  ,i,j,XDIR) + QFLX_MOMZ(KS  ,i-1,j,XDIR) ) &
                                     + J23G (KS+1,i,j,I_XYZ) * ( QFLX_MOMZ(KS+1,i,j,YDIR) + QFLX_MOMZ(KS+1,i,j-1,YDIR) ) &
                                     - J23G (KS  ,i,j,I_XYZ) * ( QFLX_MOMZ(KS  ,i,j,YDIR) + QFLX_MOMZ(KS  ,i,j-1,YDIR) ) &
                                    ) * 0.5_RP * RCDZ(KS+1) &
                                  + J33G * ( QFLX_MOMZ(KS+1,i,j,ZDIR) - QFLX_MOMZ(KS,i,j,ZDIR) ) * RFDZ(KS) &
                                  ) / GSQRT(KS,i,j,I_XYW)
             MOMZ_t_TB(KE-1,i,j) = - ( ( GSQRT(KE-1,i  ,j,I_UYW) * QFLX_MOMZ(KE-1,i  ,j,XDIR) &
                                       - GSQRT(KE-1,i-1,j,I_UYW) * QFLX_MOMZ(KE-1,i-1,j,XDIR) ) * RCDX(i) &
                                     + ( GSQRT(KE-1,i,j  ,I_XVW) * QFLX_MOMZ(KE-1,i,j  ,YDIR) &
                                       - GSQRT(KE-1,i,j-1,I_XVW) * QFLX_MOMZ(KE-1,i,j-1,YDIR) ) * RCDY(j) &
                                     + ( J13G (KE-1,i,j,I_XYZ) * ( QFLX_MOMZ(KE-1,i,j,XDIR) + QFLX_MOMZ(KE-1,i-1,j,XDIR) ) &
                                       - J13G (KE-2,i,j,I_XYZ) * ( QFLX_MOMZ(KE-2,i,j,XDIR) + QFLX_MOMZ(KE-2,i-1,j,XDIR) ) &
                                       + J23G (KE-1,i,j,I_XYZ) * ( QFLX_MOMZ(KE-1,i,j,YDIR) + QFLX_MOMZ(KE-1,i,j-1,YDIR) ) &
                                       - J23G (KE-2,i,j,I_XYZ) * ( QFLX_MOMZ(KE-2,i,j,YDIR) + QFLX_MOMZ(KE-2,i,j-1,YDIR) ) &
                                      ) * 0.5_RP * RCDZ(KE-1) &
                                    + J33G * ( QFLX_MOMZ(KE,i,j,ZDIR) - QFLX_MOMZ(KE-1,i,j,ZDIR) ) * RFDZ(KE-1) &
                                    ) / GSQRT(KE-1,i,j,I_XYW)
          enddo
          enddo

          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-1
             MOMX_t_TB(k,i,j) = - ( ( GSQRT(k,i+1,j,I_XYZ) * QFLX_MOMX(k,i+1,j,XDIR) &
                                    - GSQRT(k,i  ,j,I_XYZ) * QFLX_MOMX(k,i  ,j,XDIR) ) * RFDX(i) &
                                  + ( GSQRT(k,i,j  ,I_UVZ) * QFLX_MOMX(k,i,j  ,YDIR) &
                                    - GSQRT(k,i,j-1,I_UVZ) * QFLX_MOMX(k,i,j-1,YDIR) ) * RCDY(j) &
                                  + ( J13G (k+1,i,j,I_UYW) * ( QFLX_MOMX(k+1,i+1,j,XDIR) + QFLX_MOMX(k+1,i,j  ,XDIR) ) &
                                    - J13G (k-1,i,j,I_UYW) * ( QFLX_MOMX(k-1,i+1,j,XDIR) + QFLX_MOMX(k-1,i,j  ,XDIR) ) &
                                    + J23G (k+1,i,j,I_UYW) * ( QFLX_MOMX(k+1,i  ,j,YDIR) + QFLX_MOMX(k+1,i,j-1,YDIR) ) &
                                    - J23G (k-1,i,j,I_UYW) * ( QFLX_MOMX(k-1,i  ,j,YDIR) + QFLX_MOMX(k-1,i,j-1,YDIR) ) &
                                    ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
                                  + J33G * ( QFLX_MOMX(k,i,j,ZDIR) - QFLX_MOMX(k-1,i,j,ZDIR) ) * RCDZ(k) &
                                  ) / GSQRT(k,i,j,I_UYZ)
          enddo
          enddo
          enddo
          do j = JJS, JJE
          do i = IIS, IIE
             MOMX_t_TB(KS,i,j) = - ( ( GSQRT(KS,i+1,j,I_XYZ) * QFLX_MOMX(KS,i+1,j,XDIR) &
                                     - GSQRT(KS,i  ,j,I_XYZ) * QFLX_MOMX(KS,i  ,j,XDIR) ) * RFDX(i) &
                                   + ( GSQRT(KS,i,j  ,I_UVZ) * QFLX_MOMX(KS,i,j  ,YDIR) &
                                     - GSQRT(KS,i,j-1,I_UVZ) * QFLX_MOMX(KS,i,j-1,YDIR) ) * RCDY(j) &
                                   + ( J13G (KS+1,i,j,I_UYW) * ( QFLX_MOMX(KS+1,i+1,j,XDIR) + QFLX_MOMX(KS+1,i,j  ,XDIR) ) &
                                     - J13G (KS  ,i,j,I_UYW) * ( QFLX_MOMX(KS  ,i+1,j,XDIR) + QFLX_MOMX(KS  ,i,j  ,XDIR) ) &
                                     + J23G (KS+1,i,j,I_UYW) * ( QFLX_MOMX(KS+1,i  ,j,YDIR) + QFLX_MOMX(KS+1,i,j-1,YDIR) ) &
                                     - J23G (KS  ,i,j,I_UYW) * ( QFLX_MOMX(KS  ,i  ,j,YDIR) + QFLX_MOMX(KS  ,i,j-1,YDIR) ) &
                                    ) * 0.5_RP * RCDZ(KS) &
                                   + J33G * ( QFLX_MOMX(KS,i,j,ZDIR) ) * RFDZ(KS) &
                                  ) / GSQRT(KS,i,j,I_UYZ)
             MOMX_t_TB(KE,i,j) = - ( ( GSQRT(KE,i+1,j,I_XYZ) * QFLX_MOMX(KE,i+1,j,XDIR) &
                                     - GSQRT(KE,i  ,j,I_XYZ) * QFLX_MOMX(KE,i  ,j,XDIR) ) * RFDX(i) &
                                   + ( GSQRT(KE,i,j  ,I_UVZ) * QFLX_MOMX(KE,i,j  ,YDIR) &
                                     - GSQRT(KE,i,j-1,I_UVZ) * QFLX_MOMX(KE,i,j-1,YDIR) ) * RCDY(j) &
                                   + ( J13G (KE  ,i,j,I_UYW) * ( QFLX_MOMX(KE  ,i+1,j,XDIR) + QFLX_MOMX(KE  ,i,j  ,XDIR) ) &
                                     - J13G (KE-1,i,j,I_UYW) * ( QFLX_MOMX(KE-1,i+1,j,XDIR) + QFLX_MOMX(KE-1  ,i,j  ,XDIR) ) &
                                     + J23G (KE  ,i,j,I_UYW) * ( QFLX_MOMX(KE  ,i  ,j,YDIR) + QFLX_MOMX(KE-1+1,i,j-1,YDIR) ) &
                                     - J23G (KE-1,i,j,I_UYW) * ( QFLX_MOMX(KE-1,i  ,j,YDIR) + QFLX_MOMX(KE-1  ,i,j-1,YDIR) ) &
                                    ) * 0.5_RP * RFDZ(KE-1) &
                                   - J33G * ( QFLX_MOMX(KE-1,i,j,ZDIR) ) * RCDZ(KE) &
                                  ) / GSQRT(KE,i,j,I_UYZ)
          enddo
          enddo

          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-1
             MOMY_t_TB(k,i,j) = - ( ( GSQRT(k,i  ,j  ,I_UVZ) * QFLX_MOMY(k,i  ,j,XDIR) &
                                    - GSQRT(k,i-1,j  ,I_UVZ) * QFLX_MOMY(k,i-1,j,XDIR) ) * RCDX(i) &
                                  + ( GSQRT(k,i  ,j+1,I_XYZ) * QFLX_MOMY(k,i,j+1,YDIR) &
                                    - GSQRT(k,i  ,j  ,I_XYZ) * QFLX_MOMY(k,i,j  ,YDIR) ) * RFDY(j) &
                                  + ( J13G (k+1,i,j  ,I_XVW) * ( QFLX_MOMY(k+1,i,j  ,XDIR) + QFLX_MOMY(k+1,i-1,j,XDIR) ) &
                                    - J13G (k-1,i,j  ,I_XVW) * ( QFLX_MOMY(k-1,i,j  ,XDIR) + QFLX_MOMY(k-1,i-1,j,XDIR) ) &
                                    + J23G (k+1,i,j+1,I_XVW) * ( QFLX_MOMY(k+1,i,j+1,YDIR) + QFLX_MOMY(k+1,i  ,j,YDIR) ) &
                                    - J23G (k-1,i,j+1,I_XVW) * ( QFLX_MOMY(k-1,i,j+1,YDIR) + QFLX_MOMY(k-1,i  ,j,YDIR) ) &
                                    ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
                                  + J33G * ( QFLX_MOMY(k,i,j,ZDIR) - QFLX_MOMY(k-1,i,j,ZDIR) ) * RCDZ(k) &
                                  ) / GSQRT(k,i,j,I_XVW)
          enddo
          enddo
          enddo
          do j = JJS, JJE
          do i = IIS, IIE
             MOMY_t_TB(KS,i,j) = - ( ( GSQRT(KS,i  ,j  ,I_UVZ) * QFLX_MOMY(KS,i  ,j,XDIR) &
                                     - GSQRT(KS,i-1,j  ,I_UVZ) * QFLX_MOMY(KS,i-1,j,XDIR) ) * RCDX(i) &
                                   + ( GSQRT(KS,i  ,j+1,I_XYZ) * QFLX_MOMY(KS,i,j+1,YDIR) &
                                     - GSQRT(KS,i  ,j  ,I_XYZ) * QFLX_MOMY(KS,i,j  ,YDIR) ) * RFDY(j) &
                                   + ( J13G (KS+1,i,j  ,I_XVW) * ( QFLX_MOMY(KS+1,i,j  ,XDIR) + QFLX_MOMY(KS+1,i-1,j,XDIR) ) &
                                     - J13G (KS  ,i,j  ,I_XVW) * ( QFLX_MOMY(KS  ,i,j  ,XDIR) + QFLX_MOMY(KS  ,i-1,j,XDIR) ) &
                                     + J23G (KS+1,i,j+1,I_XVW) * ( QFLX_MOMY(KS+1,i,j+1,YDIR) + QFLX_MOMY(KS+1,i  ,j,YDIR) ) &
                                     - J23G (KS  ,i,j+1,I_XVW) * ( QFLX_MOMY(KS  ,i,j+1,YDIR) + QFLX_MOMY(KS  ,i  ,j,YDIR) ) &
                                     ) * 0.5_RP * RFDZ(KS) &
                                   + J33G * ( QFLX_MOMY(KS,i,j,ZDIR) ) * RCDZ(KS) &
                                   ) / GSQRT(KS,i,j,I_XVW)
             MOMY_t_TB(KE,i,j) = - ( ( GSQRT(KE,i  ,j  ,I_UVZ) * QFLX_MOMY(KE,i  ,j,XDIR) &
                                     - GSQRT(KE,i-1,j  ,I_UVZ) * QFLX_MOMY(KE,i-1,j,XDIR) ) * RCDX(i) &
                                   + ( GSQRT(KE,i  ,j+1,I_XYZ) * QFLX_MOMY(KE,i,j+1,YDIR) &
                                     - GSQRT(KE,i  ,j  ,I_XYZ) * QFLX_MOMY(KE,i,j  ,YDIR) ) * RFDY(j) &
                                   + ( J13G (KE  ,i,j  ,I_XVW) * ( QFLX_MOMY(KE  ,i,j  ,XDIR) + QFLX_MOMY(KE  ,i-1,j,XDIR) ) &
                                     - J13G (KE-1,i,j  ,I_XVW) * ( QFLX_MOMY(KE-1,i,j  ,XDIR) + QFLX_MOMY(KE-1,i-1,j,XDIR) ) &
                                     + J23G (KE  ,i,j+1,I_XVW) * ( QFLX_MOMY(KE  ,i,j+1,YDIR) + QFLX_MOMY(KE  ,i  ,j,YDIR) ) &
                                     - J23G (KE-1,i,j+1,I_XVW) * ( QFLX_MOMY(KE-1,i,j+1,YDIR) + QFLX_MOMY(KE-1,i  ,j,YDIR) ) &
                                     ) * 0.5_RP * RFDZ(KE-1) &
                                   - J33G * ( QFLX_MOMY(KE-1,i,j,ZDIR) ) * RCDZ(KE) &
                                   ) / GSQRT(KE,i,j,I_XVW)
          enddo
          enddo

          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-1
             RHOT_t_TB(k,i,j) = - ( ( GSQRT(k,i  ,j,I_UYZ) * QFLX_RHOT(k,i  ,j,XDIR) &
                                    - GSQRT(k,i-1,j,I_UVZ) * QFLX_RHOT(k,i-1,j,XDIR) ) * RCDX(i) &
                                  + ( GSQRT(k,i,j  ,I_XVZ) * QFLX_RHOT(k,i,j  ,YDIR) &
                                    - GSQRT(k,i,j-1,I_XVZ) * QFLX_RHOT(k,i,j-1,YDIR) ) * RCDY(j) &
                                  + ( J13G(k+1,i,j,I_XYW) * ( QFLX_RHOT(k+1,i,j,XDIR) + QFLX_RHOT(k+1,i-1,j,XDIR) ) &
                                    - J13G(k-1,i,j,I_XYW) * ( QFLX_RHOT(k-1,i,j,XDIR) + QFLX_RHOT(k-1,i-1,j,XDIR) ) &
                                    + J23G(k+1,i,j,I_XYW) * ( QFLX_RHOT(k+1,i,j,YDIR) + QFLX_RHOT(k+1,i,j-1,YDIR) ) &
                                    - J23G(k-1,i,j,I_XYW) * ( QFLX_RHOT(k-1,i,j,YDIR) + QFLX_RHOT(k-1,i,j-1,YDIR) ) &
                                    ) * 0.5_RP / ( FDZ(k) + FDZ(k-1) ) &
                                  + J33G * ( QFLX_RHOT(k,i,j,ZDIR) - QFLX_RHOT(k-1,i,j,ZDIR) ) * RCDZ(k) &
                                  ) / GSQRT(k,i,j,I_XYZ)
          enddo
          enddo
          enddo
          do j = JJS, JJE
          do i = IIS, IIE
             RHOT_t_TB(KS,i,j) = - ( ( GSQRT(KS,i  ,j,I_UYZ) * QFLX_RHOT(KS,i  ,j,XDIR) &
                                     - GSQRT(KS,i-1,j,I_UVZ) * QFLX_RHOT(KS,i-1,j,XDIR) ) * RCDX(i) &
                                   + ( GSQRT(KS,i,j  ,I_XVZ) * QFLX_RHOT(KS,i,j  ,YDIR) &
                                     - GSQRT(KS,i,j-1,I_XVZ) * QFLX_RHOT(KS,i,j-1,YDIR) ) * RCDY(j) &
                                   + ( J13G(KS+1,i,j,I_XYW) * ( QFLX_RHOT(KS+1,i,j,XDIR) + QFLX_RHOT(KS+1,i-1,j,XDIR) ) &
                                     - J13G(KS  ,i,j,I_XYW) * ( QFLX_RHOT(KS  ,i,j,XDIR) + QFLX_RHOT(KS  ,i-1,j,XDIR) ) &
                                     + J23G(KS+1,i,j,I_XYW) * ( QFLX_RHOT(KS+1,i,j,YDIR) + QFLX_RHOT(KS+1,i,j-1,YDIR) ) &
                                     - J23G(KS  ,i,j,I_XYW) * ( QFLX_RHOT(KS  ,i,j,YDIR) + QFLX_RHOT(KS  ,i,j-1,YDIR) ) &
                                    ) * 0.5_RP * RFDZ(KS) &
                                   + J33G * ( QFLX_RHOT(KS,i,j,ZDIR) ) * RCDZ(KS) &
                                   ) / GSQRT(KS,i,j,I_XYZ)
             RHOT_t_TB(KE,i,j) = - ( ( GSQRT(KE,i  ,j,I_UYZ) * QFLX_RHOT(KE,i  ,j,XDIR) &
                                     - GSQRT(KE,i-1,j,I_UVZ) * QFLX_RHOT(KE,i-1,j,XDIR) ) * RCDX(i) &
                                   + ( GSQRT(KE,i,j  ,I_XVZ) * QFLX_RHOT(KE,i,j  ,YDIR) &
                                     - GSQRT(KE,i,j-1,I_XVZ) * QFLX_RHOT(KE,i,j-1,YDIR) ) * RCDY(j) &
                                   + ( J13G(KE  ,i,j,I_XYW) * ( QFLX_RHOT(KE  ,i,j,XDIR) + QFLX_RHOT(KE  ,i-1,j,XDIR) ) &
                                     - J13G(KE-1,i,j,I_XYW) * ( QFLX_RHOT(KE-1,i,j,XDIR) + QFLX_RHOT(KE-1,i-1,j,XDIR) ) &
                                     + J23G(KE  ,i,j,I_XYW) * ( QFLX_RHOT(KE  ,i,j,YDIR) + QFLX_RHOT(KE  ,i,j-1,YDIR) ) &
                                     - J23G(KE-1,i,j,I_XYW) * ( QFLX_RHOT(KE-1,i,j,YDIR) + QFLX_RHOT(KE-1,i,j-1,YDIR) ) &
                                    ) * 0.5_RP * RFDZ(KE-1) &
                                   - J33G * ( QFLX_RHOT(KE-1,i,j,ZDIR) ) * RCDZ(KE) &
                                   ) / GSQRT(KE,i,j,I_XYZ)
          enddo
          enddo

          do iq = 1, QA
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS+1, KE-1
               QTRC_t_TB(k,i,j,iq) = - ( ( GSQRT(k,i  ,j,I_UYZ) * QFLX_QTRC(k,i  ,j,iq,XDIR) &
                                         - GSQRT(k,i-1,j,I_UVZ) * QFLX_QTRC(k,i-1,j,iq,XDIR) ) * RCDX(i) &
                                       + ( GSQRT(k,i,j  ,I_XVZ) * QFLX_QTRC(k,i,j  ,iq,YDIR) &
                                         - GSQRT(k,i,j-1,I_XVZ) * QFLX_QTRC(k,i,j-1,iq,YDIR) ) * RCDY(j) &
                                       + ( J13G(k+1,i,j,I_XYW) * ( QFLX_QTRC(k+1,i,j,iq,XDIR) + QFLX_QTRC(k+1,i-1,j,iq,XDIR) ) &
                                         - J13G(k-1,i,j,I_XYW) * ( QFLX_QTRC(k-1,i,j,iq,XDIR) + QFLX_QTRC(k-1,i-1,j,iq,XDIR) ) &
                                         + J23G(k+1,i,j,I_XYW) * ( QFLX_QTRC(k+1,i,j,iq,YDIR) + QFLX_QTRC(k+1,i,j-1,iq,YDIR) ) &
                                         - J23G(k-1,i,j,I_XYW) * ( QFLX_QTRC(k-1,i,j,iq,YDIR) + QFLX_QTRC(k-1,i,j-1,iq,YDIR) ) &
                                         ) * 0.5_RP / ( FDZ(k) + FDZ(k-1) ) &
                                       + J33G * ( QFLX_QTRC(k,i,j,iq,ZDIR) - QFLX_QTRC(k-1,i,j,iq,ZDIR) ) * RCDZ(k) &
                                       ) / GSQRT(k,i,j,I_XYZ)
             enddo
             enddo
             enddo
             do j = JJS, JJE
             do i = IIS, IIE
               QTRC_t_TB(KS,i,j,iq) = - ( ( GSQRT(KS,i  ,j,I_UYZ) * QFLX_QTRC(KS,i  ,j,iq,XDIR) &
                                          - GSQRT(KS,i-1,j,I_UVZ) * QFLX_QTRC(KS,i-1,j,iq,XDIR) ) * RCDX(i) &
                                        + ( GSQRT(KS,i,j  ,I_XVZ) * QFLX_QTRC(KS,i,j  ,iq,YDIR) &
                                          - GSQRT(KS,i,j-1,I_XVZ) * QFLX_QTRC(KS,i,j-1,iq,YDIR) ) * RCDY(j) &
                                        + ( J13G(KS+1,i,j,I_XYW) * ( QFLX_QTRC(KS+1,i,j,iq,XDIR) + QFLX_QTRC(KS+1,i-1,j,iq,XDIR) ) &
                                          - J13G(KS  ,i,j,I_XYW) * ( QFLX_QTRC(KS  ,i,j,iq,XDIR) + QFLX_QTRC(KS  ,i-1,j,iq,XDIR) ) &
                                          + J23G(KS+1,i,j,I_XYW) * ( QFLX_QTRC(KS+1,i,j,iq,YDIR) + QFLX_QTRC(KS+1,i,j-1,iq,YDIR) ) &
                                          - J23G(KS  ,i,j,I_XYW) * ( QFLX_QTRC(KS,i,j,iq,YDIR) + QFLX_QTRC(KS,i,j-1,iq,YDIR) ) &
                                          ) * 0.5_RP * RFDZ(KS) &
                                        + J33G * ( QFLX_QTRC(KS,i,j,iq,ZDIR) ) * RCDZ(KS) &
                                         ) / GSQRT(KS,i,j,I_XYZ)
               QTRC_t_TB(KE,i,j,iq) = - ( ( GSQRT(KE,i  ,j,I_UYZ) * QFLX_QTRC(KE,i  ,j,iq,XDIR) &
                                          - GSQRT(KE,i-1,j,I_UVZ) * QFLX_QTRC(KE,i-1,j,iq,XDIR) ) * RCDX(i) &
                                        + ( GSQRT(KE,i,j  ,I_XVZ) * QFLX_QTRC(KE,i,j  ,iq,YDIR) &
                                          - GSQRT(KE,i,j-1,I_XVZ) * QFLX_QTRC(KE,i,j-1,iq,YDIR) ) * RCDY(j) &
                                        + ( J13G(KE  ,i,j,I_XYW) * ( QFLX_QTRC(KE  ,i,j,iq,XDIR) + QFLX_QTRC(KE  ,i-1,j,iq,XDIR) ) &
                                          - J13G(KE-1,i,j,I_XYW) * ( QFLX_QTRC(KE-1,i,j,iq,XDIR) + QFLX_QTRC(KE-1,i-1,j,iq,XDIR) ) &
                                          + J23G(KE  ,i,j,I_XYW) * ( QFLX_QTRC(KE  ,i,j,iq,YDIR) + QFLX_QTRC(KE  ,i,j-1,iq,YDIR) ) &
                                          - J23G(KE-1,i,j,I_XYW) * ( QFLX_QTRC(KE-1,i,j,iq,YDIR) + QFLX_QTRC(KE-1,i,j-1,iq,YDIR) ) &
                                          ) * 0.5_RP * RFDZ(KE-1) &
                                        - J33G * ( QFLX_QTRC(KE-1,i,j,iq,ZDIR) ) * RCDZ(KE) &
                                         ) / GSQRT(KE,i,j,I_XYZ)
             enddo
             enddo
          enddo

       enddo
       enddo

       if ( history_flag ) then

          call HIST_in( TKE(:,:,:), 'TKE', 'turburent kinetic energy', 'm2/s2', dt_TB )
          call HIST_in( NU (:,:,:), 'NU',  'eddy viscosity',           'm2/s',  dt_TB )
          call HIST_in( Ri (:,:,:), 'Ri',  'Richardson number',        'NIL',   dt_TB )
          call HIST_in( Pr (:,:,:), 'Pr',  'Prantle number',           'NIL',   dt_TB )

          call HIST_in( MOMZ_t_TB(:,:,:), 'MOMZ_t_TB', 'MOMZ tendency (TB)', 'kg/m2/s2', dt_TB )
          call HIST_in( MOMX_t_TB(:,:,:), 'MOMX_t_TB', 'MOMX tendency (TB)', 'kg/m2/s2', dt_TB )
          call HIST_in( MOMY_t_TB(:,:,:), 'MOMY_t_TB', 'MOMY tendency (TB)', 'kg/m2/s2', dt_TB )
          call HIST_in( RHOT_t_TB(:,:,:), 'RHOT_t_TB', 'RHOT tendency (TB)', 'kg/m2/s2', dt_TB )

          call HIST_in( QFLX_MOMZ(:,:,:,ZDIR), 'SGS_ZFLX_MOMZ', 'SGS Z FLUX of MOMZ', 'kg/m/s2', dt_TB)
          call HIST_in( QFLX_MOMZ(:,:,:,XDIR), 'SGS_XFLX_MOMZ', 'SGS X FLUX of MOMZ', 'kg/m/s2', dt_TB, xdim='half', zdim='half')
          call HIST_in( QFLX_MOMZ(:,:,:,YDIR), 'SGS_YFLX_MOMZ', 'SGS Y FLUX of MOMZ', 'kg/m/s2', dt_TB, ydim='half', zdim='half')

          call HIST_in( QFLX_MOMX(:,:,:,ZDIR), 'SGS_ZFLX_MOMX', 'SGS Z FLUX of MOMX', 'kg/m/s2', dt_TB, xdim='half', zdim='half')
          call HIST_in( QFLX_MOMX(:,:,:,XDIR), 'SGS_XFLX_MOMX', 'SGS X FLUX of MOMX', 'kg/m/s2', dt_TB)
          call HIST_in( QFLX_MOMX(:,:,:,YDIR), 'SGS_YFLX_MOMX', 'SGS Y FLUX of MOMX', 'kg/m/s2', dt_TB, xdim='half', ydim='half')

          call HIST_in( QFLX_MOMY(:,:,:,ZDIR), 'SGS_ZFLX_MOMY', 'SGS Z FLUX of MOMY', 'kg/m/s2', dt_TB, ydim='half', zdim='half')
          call HIST_in( QFLX_MOMY(:,:,:,XDIR), 'SGS_XFLX_MOMY', 'SGS X FLUX of MOMY', 'kg/m/s2', dt_TB, xdim='half', ydim='half')
          call HIST_in( QFLX_MOMY(:,:,:,YDIR), 'SGS_YFLX_MOMY', 'SGS Y FLUX of MOMY', 'kg/m/s2', dt_TB)

          call HIST_in( QFLX_RHOT(:,:,:,ZDIR), 'SGS_ZFLX_RHOT', 'SGS Z FLUX of RHOT', 'K*kg/m2/s', dt_TB, zdim='half')
          call HIST_in( QFLX_RHOT(:,:,:,XDIR), 'SGS_XFLX_RHOT', 'SGS X FLUX of RHOT', 'K*kg/m2/s', dt_TB, xdim='half')
          call HIST_in( QFLX_RHOT(:,:,:,YDIR), 'SGS_YFLX_RHOT', 'SGS Y FLUX of RHOT', 'K*kg/m2/s', dt_TB, ydim='half')

          if ( I_QV > 0 ) then
             call HIST_in( QFLX_QTRC(:,:,:,I_QV,ZDIR), 'SGS_ZFLX_QV', 'SGS Z FLUX of QV', 'kg/m2/s', dt_TB, zdim='half')
             call HIST_in( QFLX_QTRC(:,:,:,I_QV,XDIR), 'SGS_XFLX_QV', 'SGS X FLUX of QV', 'kg/m2/s', dt_TB, xdim='half')
             call HIST_in( QFLX_QTRC(:,:,:,I_QV,YDIR), 'SGS_YFLX_QV', 'SGS Y FLUX of QV', 'kg/m2/s', dt_TB, ydim='half')
          endif

       endif

    endif

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_TB(k,i,j)
       MOMX_t(k,i,j) = MOMX_t(k,i,j) + MOMX_t_TB(k,i,j)
       MOMY_t(k,i,j) = MOMY_t(k,i,j) + MOMY_t_TB(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_TB(k,i,j)
    enddo
    enddo
    enddo

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do iq = 1,  QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       QTRC_t(k,i,j,iq) = QTRC_t(k,i,j,iq) + QTRC_t_TB(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, MOMZ_t_TB(:,:,:), 'MOMZ_t_TB' )
       call STAT_total( total, MOMX_t_TB(:,:,:), 'MOMX_t_TB' )
       call STAT_total( total, MOMY_t_TB(:,:,:), 'MOMY_t_TB' )
       call STAT_total( total, RHOT_t_TB(:,:,:), 'RHOT_t_TB' )

       do iq = 1, QA
          RHOQ(:,:,:) = DENS_av(:,:,:) * QTRC_t_TB(:,:,:,iq)

          call STAT_total( total, RHOQ(:,:,:), trim(AQ_NAME(iq))//'_t_TB' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_TB_driver

end module mod_atmos_phy_tb_driver
