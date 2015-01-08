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
       CDY => GRID_CDY
    use scale_grid_real, only: &
       CZ  => REAL_CZ
    use scale_process, only: &
       PRC_MPIstop
    use mod_atmos_admin, only: &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_sw_phy_tb
    use scale_atmos_phy_tb, only: &
       ATMOS_PHY_TB_setup
    use mod_atmos_phy_tb_vars, only: &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t
    implicit none

!    NAMELIST / PARAM_ATMOS_PHY_TB / &

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_TB] / Origin[SCALE-LES]'

!    !--- read namelist
!    rewind(IO_FID_CONF)
!    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB,iostat=ierr)
!    if( ierr < 0 ) then !--- missing
!       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
!    elseif( ierr > 0 ) then !--- fatal error
!       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB. Check!'
!       call PRC_MPIstop
!    endif
!    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB)


    ! initialize
    do j = JS, JE
    do i = IS, IE
       MOMZ_t_TB(KS-1,i,j) = 0.0_RP
       MOMZ_t_TB(KE  ,i,j) = 0.0_RP
    enddo
    enddo

    if ( ATMOS_sw_phy_tb ) then

       ! setup library component
       call ATMOS_PHY_TB_setup( ATMOS_PHY_TB_TYPE, & ! [IN]
                                CDZ, CDX, CDY, CZ  ) ! [IN]

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM Turbulence', 1)
       call ATMOS_PHY_TB_driver( update_flag = .true. )
       call PROF_rapend  ('ATM Turbulence', 1)

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_TB_driver_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_driver( update_flag )
  !> Driver
    use scale_gridtrans, only: &
       GSQRT => GTRANS_GSQRT, &
       J13G  => GTRANS_J13G,  &
       J23G  => GTRANS_J23G,  &
       J33G  => GTRANS_J33G,  &
       MAPF  => GTRANS_MAPF
    use scale_statistics, only: &
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
       ATMOS_PHY_TB
    use scale_atmos_phy_tb_common, only: &
       calc_tend_momz => ATMOS_PHY_TB_calc_tend_momz, &
       calc_tend_momx => ATMOS_PHY_TB_calc_tend_momx, &
       calc_tend_momy => ATMOS_PHY_TB_calc_tend_momy, &
       calc_tend_phi  => ATMOS_PHY_TB_calc_tend_phi
    use mod_atmos_vars, only: &
       DENS => DENS_av,   &
       MOMZ => MOMZ_av,   &
       MOMX => MOMX_av,   &
       MOMY => MOMY_av,   &
       RHOT => RHOT_av,   &
       QTRC => QTRC_av,   &
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
       RHOQ_t_TB => ATMOS_PHY_TB_RHOQ_t, &
       TKE       => ATMOS_PHY_TB_TKE,    &
       NU        => ATMOS_PHY_TB_NU
    use mod_atmos_phy_sf_vars, only: &
       SFLX_MW => ATMOS_PHY_SF_SFLX_MW, &
       SFLX_MU => ATMOS_PHY_SF_SFLX_MU, &
       SFLX_MV => ATMOS_PHY_SF_SFLX_MV, &
       SFLX_SH => ATMOS_PHY_SF_SFLX_SH
    implicit none

    logical, intent(in) :: update_flag

    ! eddy viscosity/diffusion flux
    real(RP) :: QFLX_MOMZ(KA,IA,JA,3)
    real(RP) :: QFLX_MOMX(KA,IA,JA,3)
    real(RP) :: QFLX_MOMY(KA,IA,JA,3)
    real(RP) :: QFLX_RHOT(KA,IA,JA,3)
    real(RP) :: QFLX_RHOQ(KA,IA,JA,3,QA)

    real(RP) :: Ri(KA,IA,JA)
    real(RP) :: Pr(KA,IA,JA)

    integer :: JJS, JJE
    integer :: IIS, IIE

    real(RP) :: total ! dummy

    integer :: n, k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call ATMOS_PHY_TB( QFLX_MOMZ, & ! [OUT]
                          QFLX_MOMX, & ! [OUT]
                          QFLX_MOMY, & ! [OUT]
                          QFLX_RHOT, & ! [OUT]
                          QFLX_RHOQ, & ! [OUT]
                          TKE,       & ! [OUT]
                          NU,        & ! [OUT]
                          Ri,        & ! [OUT]
                          Pr,        & ! [OUT]
                          MOMZ,      & ! [IN]
                          MOMX,      & ! [IN]
                          MOMY,      & ! [IN]
                          RHOT,      & ! [IN]
                          DENS,      & ! [IN]
                          QTRC,      & ! [IN]
                          SFLX_MW,   & ! [IN]
                          SFLX_MU,   & ! [IN]
                          SFLX_MV,   & ! [IN]
                          SFLX_SH,   & ! [IN]
                          GSQRT,     & ! [IN]
                          J13G,      & ! [IN]
                          J23G,      & ! [IN]
                          J33G,      & ! [IN]
                          MAPF,      & ! [IN]
                          dt_TB      ) ! [IN]

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
          do JJS = JS, JE, JBLOCK
          JJE = JJS+JBLOCK-1
          do IIS = IS, IE, IBLOCK
          IIE = IIS+IBLOCK-1
             call calc_tend_phi ( RHOQ_t_TB(:,:,:,iq),    & ! (out)
                                  QFLX_RHOQ(:,:,:,:,iq),  & ! (in)
                                  GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                                  IIS, IIE, JJS, JJE ) ! (in)
          end do
          end do
       end do


       call COMM_vars8( TKE(:,:,:), 1 )
       call COMM_wait ( TKE(:,:,:), 1 )


       call HIST_in( TKE(:,:,:), 'TKE', 'turburent kinetic energy', 'm2/s2', nohalo=.true. )
       call HIST_in( NU (:,:,:), 'NU',  'eddy viscosity',           'm2/s' , nohalo=.true. )
       call HIST_in( Ri (:,:,:), 'Ri',  'Richardson number',        'NIL'  , nohalo=.true. )
       call HIST_in( Pr (:,:,:), 'Pr',  'Prantle number',           'NIL'  , nohalo=.true. )

       call HIST_in( MOMZ_t_TB(:,:,:), 'MOMZ_t_TB', 'MOMZ tendency (TB)', 'kg/m2/s2', nohalo=.true. )
       call HIST_in( MOMX_t_TB(:,:,:), 'MOMX_t_TB', 'MOMX tendency (TB)', 'kg/m2/s2', nohalo=.true. )
       call HIST_in( MOMY_t_TB(:,:,:), 'MOMY_t_TB', 'MOMY tendency (TB)', 'kg/m2/s2', nohalo=.true. )
       call HIST_in( RHOT_t_TB(:,:,:), 'RHOT_t_TB', 'RHOT tendency (TB)', 'K*kg/m3/s', nohalo=.true. )

       do iq = 1, QA
          call HIST_in( RHOQ_t_TB(:,:,:,iq), trim(AQ_NAME(iq))//'_t_TB', 'RHO*'//trim(AQ_NAME(iq))//' tendency (TB)', 'kg/m3/s', nohalo=.true. )
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

       if ( I_QV > 0 ) then
          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,I_QV), 'SGS_ZFLX_QV', 'SGS Z FLUX of QV', 'kg/m2/s', &
                        zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,I_QV), 'SGS_XFLX_QV', 'SGS X FLUX of QV', 'kg/m2/s', &
                        xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,I_QV), 'SGS_YFLX_QV', 'SGS Y FLUX of QV', 'kg/m2/s', &
                        ydim='half', nohalo=.true.)
       endif

#ifndef DRY
       if ( I_QC > 0 ) then
          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,I_QC), 'SGS_ZFLX_QC', 'SGS Z FLUX of QC', 'kg/m2/s', &
                        zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,I_QC), 'SGS_XFLX_QC', 'SGS X FLUX of QC', 'kg/m2/s', &
                        xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,I_QC), 'SGS_YFLX_QC', 'SGS Y FLUX of QC', 'kg/m2/s', &
                        ydim='half', nohalo=.true.)
       endif

       if ( I_QR > 0 ) then
          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,I_QR), 'SGS_ZFLX_QR', 'SGS Z FLUX of QR', 'kg/m2/s', &
                        zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,I_QR), 'SGS_XFLX_QR', 'SGS X FLUX of QR', 'kg/m2/s', &
                        xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,I_QR), 'SGS_YFLX_QR', 'SGS Y FLUX of QR', 'kg/m2/s', &
                        ydim='half', nohalo=.true.)
       endif

       if ( I_QI > 0 ) then
          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,I_QI), 'SGS_ZFLX_QI', 'SGS Z FLUX of QI', 'kg/m2/s', &
                        zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,I_QI), 'SGS_XFLX_QI', 'SGS X FLUX of QI', 'kg/m2/s', &
                        xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,I_QI), 'SGS_YFLX_QI', 'SGS Y FLUX of QI', 'kg/m2/s', &
                        ydim='half', nohalo=.true.)
       endif

       if ( I_QS > 0 ) then
          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,I_QS), 'SGS_ZFLX_QS', 'SGS Z FLUX of QS', 'kg/m2/s', &
                        zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,I_QS), 'SGS_XFLX_QS', 'SGS X FLUX of QS', 'kg/m2/s', &
                        xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,I_QS), 'SGS_YFLX_QS', 'SGS Y FLUX of QS', 'kg/m2/s', &
                        ydim='half', nohalo=.true.)
       endif

       if ( I_QG > 0 ) then
          call HIST_in( QFLX_RHOQ(:,:,:,ZDIR,I_QG), 'SGS_ZFLX_QG', 'SGS Z FLUX of QG', 'kg/m2/s', &
                        zdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,XDIR,I_QG), 'SGS_XFLX_QG', 'SGS X FLUX of QG', 'kg/m2/s', &
                        xdim='half', nohalo=.true.)
          call HIST_in( QFLX_RHOQ(:,:,:,YDIR,I_QG), 'SGS_YFLX_QG', 'SGS Y FLUX of QG', 'kg/m2/s', &
                        ydim='half', nohalo=.true.)
       endif
#endif

       if ( STATISTICS_checktotal ) then
          call STAT_total( total, MOMZ_t_TB(:,:,:), 'MOMZ_t_TB' )
          call STAT_total( total, MOMX_t_TB(:,:,:), 'MOMX_t_TB' )
          call STAT_total( total, MOMY_t_TB(:,:,:), 'MOMY_t_TB' )
          call STAT_total( total, RHOT_t_TB(:,:,:), 'RHOT_t_TB' )
          call STAT_total( total, TKE(:,:,:), 'TKE' )
          call STAT_total( total, Nu(:,:,:), 'Nu' )
          call STAT_total( total, Ri(:,:,:), 'Ri' )
          call STAT_total( total, Pr(:,:,:), 'Pr' )

          do iq = 1, QA
             call STAT_total( total, RHOQ_t_TB(:,:,:,iq), trim(AQ_NAME(iq))//'_t_TB' )
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
