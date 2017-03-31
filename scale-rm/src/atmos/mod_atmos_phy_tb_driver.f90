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
  public :: ATMOS_PHY_TB_driver_config
  public :: ATMOS_PHY_TB_driver_setup
  public :: ATMOS_PHY_TB_driver_resume
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
  !> Config
  subroutine ATMOS_PHY_TB_driver_config
    use scale_atmos_phy_tb, only: &
       ATMOS_PHY_TB_config
    use mod_atmos_admin, only: &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_sw_phy_tb
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CONFIG] / Categ[ATMOS PHY_TB] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_tb ) then
       call ATMOS_PHY_TB_config( ATMOS_PHY_TB_TYPE )
    end if

    return
  end subroutine ATMOS_PHY_TB_driver_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_driver_setup
    use scale_grid, only: &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY
    use scale_grid_real, only: &
       CZ  => REAL_CZ
    use scale_atmos_phy_tb, only: &
       ATMOS_PHY_TB_setup
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_tb
    use mod_atmos_phy_tb_vars, only: &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_TB] / Origin[SCALE-RM]'

    ! initialize
    do j = JS, JE
    do i = IS, IE
       MOMZ_t_TB(KS-1,i,j) = 0.0_RP
       MOMZ_t_TB(KE  ,i,j) = 0.0_RP
    enddo
    enddo

    if ( ATMOS_sw_phy_tb ) then
       ! setup library component
       call ATMOS_PHY_TB_setup( CDZ, CDX, CDY, CZ ) ! [IN]
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_TB_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_PHY_TB_driver_resume
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_tb
    implicit none

    if ( ATMOS_sw_phy_tb ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Turbulence', 1)
       call ATMOS_PHY_TB_driver( update_flag = .true. )
       call PROF_rapend  ('ATM_Turbulence', 1)

    end if

    return
  end subroutine ATMOS_PHY_TB_driver_resume

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_TB_driver( update_flag )
    use scale_gridtrans, only: &
       GSQRT => GTRANS_GSQRT, &
       J13G  => GTRANS_J13G,  &
       J23G  => GTRANS_J23G,  &
       J33G  => GTRANS_J33G,  &
       MAPF  => GTRANS_MAPF
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_history, only: &
       HIST_in
    use scale_time, only: &
       dt_TB => TIME_DTSEC_ATMOS_PHY_TB
    use scale_atmos_phy_tb, only: &
       ATMOS_PHY_TB, &
       I_TKE
    use scale_atmos_phy_tb_common, only: &
       calc_tend_momz => ATMOS_PHY_TB_calc_tend_momz, &
       calc_tend_momx => ATMOS_PHY_TB_calc_tend_momx, &
       calc_tend_momy => ATMOS_PHY_TB_calc_tend_momy, &
       calc_tend_phi  => ATMOS_PHY_TB_calc_tend_phi
    use scale_atmos_hydrometeor, only: &
       I_QV
    use mod_atmos_vars, only: &
       DENS => DENS_av,   &
       MOMZ => MOMZ_av,   &
       MOMX => MOMX_av,   &
       MOMY => MOMY_av,   &
       RHOT => RHOT_av,   &
       QTRC => QTRC_av,   &
       N2   => N2,        &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_tb_vars, only: &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t, &
       MOMX_t_TB => ATMOS_PHY_TB_MOMX_t, &
       MOMY_t_TB => ATMOS_PHY_TB_MOMY_t, &
       RHOT_t_TB => ATMOS_PHY_TB_RHOT_t, &
       RHOQ_t_TB => ATMOS_PHY_TB_RHOQ_t
    use mod_atmos_phy_sf_vars, only: &
       SFLX_MW => ATMOS_PHY_SF_SFLX_MW, &
       SFLX_MU => ATMOS_PHY_SF_SFLX_MU, &
       SFLX_MV => ATMOS_PHY_SF_SFLX_MV, &
       SFLX_SH => ATMOS_PHY_SF_SFLX_SH, &
       SFLX_Q  => ATMOS_PHY_SF_SFLX_QTRC
    implicit none

    logical, intent(in) :: update_flag

    ! eddy viscosity/diffusion flux
    real(RP) :: QFLX_MOMZ(KA,IA,JA,3)
    real(RP) :: QFLX_MOMX(KA,IA,JA,3)
    real(RP) :: QFLX_MOMY(KA,IA,JA,3)
    real(RP) :: QFLX_RHOT(KA,IA,JA,3)
    real(RP) :: QFLX_RHOQ(KA,IA,JA,3,QA)

    real(RP) :: Nu(KA,IA,JA) ! eddy viscosity
    real(RP) :: Ri(KA,IA,JA) ! Richardson number
    real(RP) :: Pr(KA,IA,JA) ! Prandtl number

    real(RP) :: tend(KA,IA,JA)
    real(RP) :: total ! dummy

    integer  :: JJS, JJE
    integer  :: IIS, IIE

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       RHOQ_t_TB = 0.0_RP

       call ATMOS_PHY_TB( QFLX_MOMZ, QFLX_MOMX, QFLX_MOMY,        & ! [OUT]
                          QFLX_RHOT, QFLX_RHOQ,                   & ! [OUT]
                          RHOQ_t_TB,                              & ! [INOUT]
                          Nu, Ri, Pr,                             & ! [OUT]
                          MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2, & ! [IN]
                          SFLX_MW, SFLX_MU, SFLX_MV,              & ! [IN]
                          SFLX_SH, SFLX_Q,                        & ! [IN]
                          GSQRT, J13G, J23G, J33G, MAPF,          & ! [IN]
                          dt_TB                                   ) ! [IN]

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1
          call calc_tend_momz( MOMZ_t_TB,   & ! (out)
                               QFLX_MOMZ,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
       end do
       end do

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1
          call calc_tend_momx( MOMX_t_TB,   & ! (out)
                               QFLX_MOMX,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
       end do
       end do

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1
          call calc_tend_momy( MOMY_t_TB,   & ! (out)
                               QFLX_MOMY,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
       end do
       end do

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1
          call calc_tend_phi ( RHOT_t_TB,  & ! (out)
                               QFLX_RHOT,  & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
       end do
       end do

       do iq = 1, QA
          if ( iq == I_TKE .or. .not. TRACER_ADVC(iq) ) cycle

          do JJS = JS, JE, JBLOCK
          JJE = JJS+JBLOCK-1
          do IIS = IS, IE, IBLOCK
          IIE = IIS+IBLOCK-1
             call calc_tend_phi( tend(:,:,:), & ! (out)
                                 QFLX_RHOQ(:,:,:,:,iq),  & ! (in)
                                 GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                                 IIS, IIE, JJS, JJE ) ! (in)

             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                RHOQ_t_TB(k,i,j,iq) = RHOQ_t_TB(k,i,j,iq) + tend(k,i,j)
             end do
             end do
             end do

          end do
          end do
       end do

       call HIST_in( NU (:,:,:), 'NU',  'eddy viscosity',           'm2/s' , nohalo=.true. )
       call HIST_in( Ri (:,:,:), 'Ri',  'Richardson number',        'NIL'  , nohalo=.true. )
       call HIST_in( Pr (:,:,:), 'Pr',  'Prantle number',           'NIL'  , nohalo=.true. )

       call HIST_in( MOMZ_t_TB(:,:,:), 'MOMZ_t_TB', 'MOMZ tendency (TB)', 'kg/m2/s2',  nohalo=.true. )
       call HIST_in( MOMX_t_TB(:,:,:), 'MOMX_t_TB', 'MOMX tendency (TB)', 'kg/m2/s2',  nohalo=.true. )
       call HIST_in( MOMY_t_TB(:,:,:), 'MOMY_t_TB', 'MOMY tendency (TB)', 'kg/m2/s2',  nohalo=.true. )
       call HIST_in( RHOT_t_TB(:,:,:), 'RHOT_t_TB', 'RHOT tendency (TB)', 'K.kg/m3/s', nohalo=.true. )

       do iq = 1, QA
          call HIST_in( RHOQ_t_TB(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_TB',                      &
                        'RHO*'//trim(TRACER_NAME(iq))//' tendency (TB)', 'kg/m3/s', nohalo=.true. )
       enddo

       call HIST_in( QFLX_MOMZ(:,:,:,ZDIR), 'SGS_ZFLX_MOMZ', 'SGS Z FLUX of MOMZ', 'kg/m/s2', &
                     nohalo=.true.)
       call HIST_in( QFLX_MOMZ(:,:,:,XDIR), 'SGS_XFLX_MOMZ', 'SGS X FLUX of MOMZ', 'kg/m/s2', &
                     xdim='half', zdim='half', nohalo=.true.)
       call HIST_in( QFLX_MOMZ(:,:,:,YDIR), 'SGS_YFLX_MOMZ', 'SGS Y FLUX of MOMZ', 'kg/m/s2', &
                     ydim='half', zdim='half', nohalo=.true.)

       call HIST_in( QFLX_MOMX(:,:,:,ZDIR), 'SGS_ZFLX_MOMX', 'SGS Z FLUX of MOMX', 'kg/m/s2', &
                     xdim='half', zdim='half', nohalo=.true.)
       call HIST_in( QFLX_MOMX(:,:,:,XDIR), 'SGS_XFLX_MOMX', 'SGS X FLUX of MOMX', 'kg/m/s2', &
                     nohalo=.true.)
       call HIST_in( QFLX_MOMX(:,:,:,YDIR), 'SGS_YFLX_MOMX', 'SGS Y FLUX of MOMX', 'kg/m/s2', &
                     xdim='half', ydim='half', nohalo=.true.)

       call HIST_in( QFLX_MOMY(:,:,:,ZDIR), 'SGS_ZFLX_MOMY', 'SGS Z FLUX of MOMY', 'kg/m/s2', &
                     ydim='half', zdim='half', nohalo=.true.)
       call HIST_in( QFLX_MOMY(:,:,:,XDIR), 'SGS_XFLX_MOMY', 'SGS X FLUX of MOMY', 'kg/m/s2', &
                     xdim='half', ydim='half', nohalo=.true.)
       call HIST_in( QFLX_MOMY(:,:,:,YDIR), 'SGS_YFLX_MOMY', 'SGS Y FLUX of MOMY', 'kg/m/s2', &
                     nohalo=.true.)

       call HIST_in( QFLX_RHOT(:,:,:,ZDIR), 'SGS_ZFLX_RHOT', 'SGS Z FLUX of RHOT', 'K*kg/m2/s', &
                     zdim='half', nohalo=.true.)
       call HIST_in( QFLX_RHOT(:,:,:,XDIR), 'SGS_XFLX_RHOT', 'SGS X FLUX of RHOT', 'K*kg/m2/s', &
                     xdim='half', nohalo=.true.)
       call HIST_in( QFLX_RHOT(:,:,:,YDIR), 'SGS_YFLX_RHOT', 'SGS Y FLUX of RHOT', 'K*kg/m2/s', &
                     ydim='half', nohalo=.true.)


       do iq = 1, QA
          if ( iq == I_TKE .or. .not. TRACER_ADVC(iq) ) cycle

          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,iq), &
               'SGS_ZFLX_'//trim(TRACER_NAME(iq)), 'SGS Z FLUX of '//trim(TRACER_NAME(iq)), 'kg/m2/s', &
               zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,iq), &
               'SGS_XFLX_'//trim(TRACER_NAME(iq)), 'SGS X FLUX of '//trim(TRACER_NAME(iq)), 'kg/m2/s', &
               xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,iq), &
               'SGS_YFLX_'//trim(TRACER_NAME(iq)), 'SGS Y FLUX of '//trim(TRACER_NAME(ia)), 'kg/m2/s', &
               ydim='half', nohalo=.true.)
       end do

       if ( STATISTICS_checktotal ) then
          call STAT_total( total, MOMZ_t_TB(:,:,:), 'MOMZ_t_TB' )
          call STAT_total( total, MOMX_t_TB(:,:,:), 'MOMX_t_TB' )
          call STAT_total( total, MOMY_t_TB(:,:,:), 'MOMY_t_TB' )
          call STAT_total( total, RHOT_t_TB(:,:,:), 'RHOT_t_TB' )
          call STAT_total( total, Nu(:,:,:), 'Nu' )
          call STAT_total( total, Ri(:,:,:), 'Ri' )
          call STAT_total( total, Pr(:,:,:), 'Pr' )

          do iq = 1, QA
             call STAT_total( total, RHOQ_t_TB(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_TB' )
          enddo
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
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_TB(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_driver

end module mod_atmos_phy_tb_driver
