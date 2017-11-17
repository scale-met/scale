!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_mp_driver
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
  public :: ATMOS_PHY_MP_driver_config
  public :: ATMOS_PHY_MP_driver_setup
  public :: ATMOS_PHY_MP_driver_resume
  public :: ATMOS_PHY_MP_driver

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
  subroutine ATMOS_PHY_MP_driver_config
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_config
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CONFIG] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_mp ) then
       call ATMOS_PHY_MP_config( ATMOS_PHY_MP_TYPE )
    end if

    return
  end subroutine ATMOS_PHY_MP_driver_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_driver_setup
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_sw_phy_mp
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_MP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_mp ) then

       ! setup library component
       call ATMOS_PHY_MP_setup

    else

       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
       if( IO_L ) write(IO_FID_LOG,*) '*** SFLX_rain and SFLX_snow is set to zero.'
       SFLX_rain(:,:) = 0.0_RP
       SFLX_snow(:,:) = 0.0_RP

    endif

    return
  end subroutine ATMOS_PHY_MP_driver_setup

  !-----------------------------------------------------------------------------
  !> resume
  subroutine ATMOS_PHY_MP_driver_resume
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_mp
    implicit none

    if ( ATMOS_sw_phy_mp ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Microphysics', 1)
       call ATMOS_PHY_MP_driver( update_flag = .true. )
       call PROF_rapend  ('ATM_Microphysics', 1)

    end if

    return
  end subroutine ATMOS_PHY_MP_driver_resume

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_MP_driver( update_flag )
    use scale_time, only: &
       dt_MP => TIME_DTSEC_ATMOS_PHY_MP
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP,                 &
       ATMOS_PHY_MP_EffectiveRadius, &
       QS_MP,                        &
       QE_MP
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       MOMZ   => MOMZ_av, &
       MOMX   => MOMX_av, &
       MOMY   => MOMY_av, &
       RHOT   => RHOT_av, &
       QTRC   => QTRC_av, &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp, &
       TEMP
    use mod_atmos_phy_mp_vars, only: &
       DENS_t_MP => ATMOS_PHY_MP_DENS_t,    &
       MOMZ_t_MP => ATMOS_PHY_MP_MOMZ_t,    &
       MOMX_t_MP => ATMOS_PHY_MP_MOMX_t,    &
       MOMY_t_MP => ATMOS_PHY_MP_MOMY_t,    &
       RHOT_t_MP => ATMOS_PHY_MP_RHOT_t,    &
       RHOQ_t_MP => ATMOS_PHY_MP_RHOQ_t,    &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE, &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_ae_vars, only: &
       CCN_t => ATMOS_PHY_AE_CCN_t
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: DENS0(KA,IA,JA)
    real(RP) :: MOMZ0(KA,IA,JA)
    real(RP) :: MOMX0(KA,IA,JA)
    real(RP) :: MOMY0(KA,IA,JA)
    real(RP) :: RHOT0(KA,IA,JA)
    real(RP) :: QTRC0(KA,IA,JA,QA)
    real(RP) :: CCN(KA,IA,JA)

    real(RP) :: Re    (KA,IA,JA,N_HYD) ! effective radius
    real(RP) :: precip(IA,JA)
    real(RP) :: total ! dummy

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

!OCL XFILL
       !$omp parallel do default(none) &
       !$omp shared(JA,IA,KA,DENS0,MOMZ0,MOMX0,MOMY0,RHOT0,DENS,MOMZ,MOMX,MOMY,RHOT) &
       !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = 1, JA
       do i  = 1, IA
       do k  = 1, KA
          DENS0(k,i,j) = DENS(k,i,j) ! save
          MOMZ0(k,i,j) = MOMZ(k,i,j) ! save
          MOMX0(k,i,j) = MOMX(k,i,j) ! save
          MOMY0(k,i,j) = MOMY(k,i,j) ! save
          RHOT0(k,i,j) = RHOT(k,i,j) ! save
       enddo
       enddo
       enddo

!OCL XFILL
       do iq = 1, QA
       do j  = 1, JA
       do i  = 1, IA
       do k  = 1, KA
          QTRC0(k,i,j,iq) = QTRC(k,i,j,iq) ! save
       enddo
       enddo
       enddo
       enddo

       CCN(:,:,:) = CCN_t(:,:,:) * dt_MP

       call ATMOS_PHY_MP( DENS0    (:,:,:),   & ! [INOUT]
                          MOMZ0    (:,:,:),   & ! [INOUT]
                          MOMX0    (:,:,:),   & ! [INOUT]
                          MOMY0    (:,:,:),   & ! [INOUT]
                          RHOT0    (:,:,:),   & ! [INOUT]
                          QTRC0    (:,:,:,:), & ! [INOUT]
                          CCN      (:,:,:),   & ! [IN]
                          EVAPORATE(:,:,:),   & ! [OUT]
                          SFLX_rain(:,:),     & ! [OUT]
                          SFLX_snow(:,:)      ) ! [OUT]

!OCL XFILL
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,KS,KE,DENS_t_MP,DENS0,DENS,MOMZ_t_MP,MOMZ0,MOMZ,MOMX_t_MP,MOMX0) &
       !$omp shared(MOMX,MOMY_t_MP,MOMY0,MOMY,RHOT_t_MP,RHOT0,RHOT,dt_MP)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          DENS_t_MP(k,i,j) = ( DENS0(k,i,j) - DENS(k,i,j) ) / dt_MP
          MOMZ_t_MP(k,i,j) = ( MOMZ0(k,i,j) - MOMZ(k,i,j) ) / dt_MP
          MOMX_t_MP(k,i,j) = ( MOMX0(k,i,j) - MOMX(k,i,j) ) / dt_MP
          MOMY_t_MP(k,i,j) = ( MOMY0(k,i,j) - MOMY(k,i,j) ) / dt_MP
          RHOT_t_MP(k,i,j) = ( RHOT0(k,i,j) - RHOT(k,i,j) ) / dt_MP
       enddo
       enddo
       enddo

!OCL XFILL
       do iq = QS_MP, QE_MP
       !$omp parallel do default(none) &
       !$omp shared(JS,JE,IS,IE,KS,KE,RHOQ_t_MP,iq,QTRC0,QTRC,DENS0,DENS,dt_MP) &
       !$omp private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t_MP(k,i,j,iq) = ( QTRC0(k,i,j,iq) * DENS0(k,i,j) &
                                - QTRC (k,i,j,iq) * DENS (k,i,j) ) / dt_MP
       enddo
       enddo
       enddo
       enddo

!OCL XFILL
       do j  = JS, JE
       do i  = IS, IE
          precip(i,j) = SFLX_rain(i,j) + SFLX_snow(i,j)
       end do
       end do

       call HIST_in( SFLX_rain(:,:),   'RAIN_MP',   'surface rain rate by MP',          'kg/m2/s',  nohalo=.true. )
       call HIST_in( SFLX_snow(:,:),   'SNOW_MP',   'surface snow rate by MP',          'kg/m2/s',  nohalo=.true. )
       call HIST_in( precip   (:,:),   'PREC_MP',   'surface precipitation rate by MP', 'kg/m2/s',  nohalo=.true. )
       call HIST_in( EVAPORATE(:,:,:), 'EVAPORATE', 'evaporated cloud number',          'num/m3/s', nohalo=.true. )

       call HIST_in( DENS_t_MP(:,:,:), 'DENS_t_MP', 'tendency DENS in MP', 'kg/m3/s'  , nohalo=.true. )
       call HIST_in( MOMZ_t_MP(:,:,:), 'MOMZ_t_MP', 'tendency MOMZ in MP', 'kg/m2/s2' , nohalo=.true. )
       call HIST_in( MOMX_t_MP(:,:,:), 'MOMX_t_MP', 'tendency MOMX in MP', 'kg/m2/s2' , nohalo=.true. )
       call HIST_in( MOMY_t_MP(:,:,:), 'MOMY_t_MP', 'tendency MOMY in MP', 'kg/m2/s2' , nohalo=.true. )
       call HIST_in( RHOT_t_MP(:,:,:), 'RHOT_t_MP', 'tendency RHOT in MP', 'K*kg/m3/s', nohalo=.true. )

       do iq = QS_MP, QE_MP
          call HIST_in( RHOQ_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP', &
                        'tendency rho*'//trim(TRACER_NAME(iq))//'in MP', 'kg/m3/s', nohalo=.true. )
       enddo

       call ATMOS_PHY_MP_EffectiveRadius( Re  (:,:,:,:), & ! [OUT]
                                          QTRC(:,:,:,:), & ! [IN]
                                          DENS(:,:,:),   & ! [IN]
                                          TEMP(:,:,:)    ) ! [IN]

       call HIST_in( Re(:,:,:,I_HC), 'Re_QC', 'effective radius cloud water', 'cm', nohalo=.true. )
       call HIST_in( Re(:,:,:,I_HR), 'Re_QR', 'effective radius rain',        'cm', nohalo=.true. )
       call HIST_in( Re(:,:,:,I_HI), 'Re_QI', 'effective radius ice water',   'cm', nohalo=.true. )
       call HIST_in( Re(:,:,:,I_HS), 'Re_QS', 'effective radius snow',        'cm', nohalo=.true. )
       call HIST_in( Re(:,:,:,I_HG), 'Re_QG', 'effective radius graupel',     'cm', nohalo=.true. )
       call HIST_in( Re(:,:,:,I_HH), 'Re_QH', 'effective radius hail',        'cm', nohalo=.true. )
    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,DENS_t,DENS_t_MP,MOMZ_t,MOMZ_t_MP,MOMX_t,MOMX_t_MP,MOMY_t) &
    !$omp shared(MOMY_t_MP,RHOT_t,RHOT_t_MP)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_MP(k,i,j)
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_MP(k,i,j)
       MOMX_t(k,i,j) = MOMX_t(k,i,j) + MOMX_t_MP(k,i,j)
       MOMY_t(k,i,j) = MOMY_t(k,i,j) + MOMY_t_MP(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_MP(k,i,j)
    enddo
    enddo
    enddo

    do iq = QS_MP, QE_MP
    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ &
    !$omp shared(JS,JE,IS,IE,KS,KE,RHOQ_t,iq,RHOQ_t_MP)
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_MP(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, DENS_t_MP(:,:,:), 'DENS_t_MP' )
       call STAT_total( total, MOMZ_t_MP(:,:,:), 'MOMZ_t_MP' )
       call STAT_total( total, MOMX_t_MP(:,:,:), 'MOMX_t_MP' )
       call STAT_total( total, MOMY_t_MP(:,:,:), 'MOMY_t_MP' )
       call STAT_total( total, RHOT_t_MP(:,:,:), 'RHOT_t_MP' )

       do iq = QS_MP, QE_MP
          call STAT_total( total, RHOQ_t_MP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_MP' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_MP_driver

end module mod_atmos_phy_mp_driver
