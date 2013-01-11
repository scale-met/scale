!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Kessler-type parametarization
!!          Reference: Kessler(1969)
!!                     Klemp and Wilhelmson(1978)
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-01-14 (Y.Miyamoto) [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2012-12-23 (H.Yashiro)  [mod] Reconstruction
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP
  public :: MP_kessler
  public :: MP_kessler_vterm

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public,  save :: vterm(KA,IA,JA,QA) ! terminal velocity of each tracer [m/s]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, save  :: MP_doreport_tendency = .false. ! report tendency of each process?
  logical, private, save  :: MP_donegative_fixer  = .true.  ! apply negative fixer?

  real(RP), private, save :: factor_vterm(KA) ! collection factor for terminal velocity of QR

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_comm, only: &
       COMM_horizontal_mean
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_MP, &
       DENS
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doreport_tendency, &
       MP_donegative_fixer

    real(RP) :: rho_prof(KA) ! averaged profile of rho

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** KESSLER-type parametarization'

    if ( trim(ATMOS_TYPE_PHY_MP) .ne. 'KESSLER' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_MP is not KESSLER. Check!'
       call PRC_MPIstop
    endif

    if (      I_QV <= 0 &
         .OR. I_QC <= 0 &
         .OR. I_QR <= 0 ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx KESSLER needs QV, QC, QR tracer. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    ! Calculate collection factor for terminal velocity of QR
    call COMM_horizontal_mean( rho_prof(:), DENS(:,:,:) )
    rho_prof(:) = rho_prof(:) * 1.E-3_RP ! [kg/m3]->[g/cc]

    do k = KS, KE
       factor_vterm(k) = sqrt( rho_prof(KS)/rho_prof(k) )
    enddo

    vterm(:,:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total,    &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    use mod_atmos_thermodyn, only: &
       THERMODYN_rhoe      => ATMOS_THERMODYN_rhoe,     &
       THERMODYN_rhot      => ATMOS_THERMODYN_rhot,     &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_history, only: &
       HIST_in
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    implicit none

    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)

    real(RP) :: temp (KA,IA,JA)
    real(RP) :: pres (KA,IA,JA)

    real(RP) :: flux_rain(KA,IA,JA)
    real(RP) :: flux_snow(KA,IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(kessler)'

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    call MP_kessler( QTRC_t(:,:,:,:), & ! [OUT]
                     RHOT  (:,:,:),   & ! [INOUT]
                     QTRC  (:,:,:,:), & ! [INOUT]
                     DENS  (:,:,:)    ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: after kessler main'
    call ATMOS_vars_total

    call MP_kessler_vterm( vterm(:,:,:,:), & ! [OUT]
                           DENS (:,:,:),   & ! [IN]
                           QTRC (:,:,:,:)  ) ! [IN]

    call THERMODYN_temp_pres( temp(:,:,:),  & ! [OUT]
                              pres(:,:,:),  & ! [OUT]
                              DENS(:,:,:),  & ! [IN]
                              RHOT(:,:,:),  & ! [IN]
                              QTRC(:,:,:,:) ) ! [IN]

    call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                         RHOT(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]
 
    call MP_precipitation( flux_rain(:,:,:), & ! [OUT]
                           flux_snow(:,:,:), & ! [OUT]
                           DENS (:,:,:),     & ! [INOUT]
                           MOMZ (:,:,:),     & ! [INOUT]
                           MOMX (:,:,:),     & ! [INOUT]
                           MOMY (:,:,:),     & ! [INOUT]
                           RHOE (:,:,:),     & ! [INOUT]
                           QTRC (:,:,:,:),   & ! [INOUT]
                           vterm(:,:,:,:),   & ! [IN]
                           temp (:,:,:)      ) ! [IN]

    call HIST_in( flux_rain(KS-1,:,:), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:), 'SNOW', 'surface snow rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:)+flux_snow(KS-1,:,:), 'PREC', 'surface precipitation rate', 'kg/m2/s', dt)


    call THERMODYN_rhot( RHOT(:,:,:),  & ! [OUT]
                         RHOE(:,:,:),  & ! [IN]
                         QTRC(:,:,:,:) ) ! [IN]

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: after precipitation'
    call ATMOS_vars_total

    call MP_saturation_adjustment( RHOT(:,:,:),   & ! [INOUT]
                                   QTRC(:,:,:,:), & ! [INOUT]
                                   DENS(:,:,:)    ) ! [IN]

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),  & ! [INOUT]
                               RHOT(:,:,:),  & ! [INOUT]
                               QTRC(:,:,:,:) ) ! [INOUT]
    endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: after mp'

    ! fill halo
    call ATMOS_vars_fillhalo

    ! log report total (optional)
    call ATMOS_vars_total

    if ( MP_doreport_tendency ) then
    endif

    return
  end subroutine ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics
  !-----------------------------------------------------------------------------
  subroutine MP_kessler( &
       QTRC_t, &
       RHOT0,  &
       QTRC0,  &
       DENS0   )
    use mod_const, only: &
       LHV00 => CONST_LH00
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use mod_atmos_thermodyn, only: &
       THERMODYN_rhoe      => ATMOS_THERMODYN_rhoe,     &
       THERMODYN_rhot      => ATMOS_THERMODYN_rhot,     &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    implicit none

    real(RP), intent(out)   :: QTRC_t(KA,IA,JA,QA)
    real(RP), intent(inout) :: RHOT0 (KA,IA,JA)
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)

    ! working
    real(RP) :: dens  ! density      [kg/m3]
    real(RP) :: rhot  ! density * PT [K*kg/m3]
    real(RP) :: q(QA) ! tracer Q [kg/kg]
    real(RP) :: temp  ! temperature [K]
    real(RP) :: pres  ! pressure [Pa]
    real(RP) :: rhoe  ! density * internal energy
    real(RP) :: qsat  ! saturated water vapor

    ! tendency
    real(RP) :: dq_evap ! tendency q (evaporation)
    real(RP) :: dq_auto ! tendency q (autoconversion)
    real(RP) :: dq_accr ! tendency q (accretion)
    real(RP) :: dqv, dqc, dqr, correct

    real(RP) :: vent_factor ! ventilation factor

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       dq_evap = 0.0_RP
       dq_auto = 0.0_RP
       dq_accr = 0.0_RP

       ! store to work
       dens = DENS0(k,i,j)
       rhot = RHOT0(k,i,j)
       do iq = I_QV, I_QR
          q(iq) = QTRC0(k,i,j,iq)
       enddo

       call THERMODYN_temp_pres( temp, & ! [OUT]
                                 pres, & ! [OUT]
                                 dens, & ! [IN]
                                 rhot, & ! [IN]
                                 q(:)  ) ! [IN]

       call THERMODYN_rhoe( rhoe, & ! [OUT]
                            rhot, & ! [IN]
                            q(:)  ) ! [IN]

       call SATURATION_dens2qsat_liq( qsat, & ! [OUT]
                                      temp, & ! [IN]
                                      dens  ) ! [IN]


       !##### Start Main #####

       ! Auto-conversion (QC->QR)
       dq_auto = 1.E-3_RP * max( q(I_QC)-1.E-3_RP, 0.0_RP )

       ! Accretion (QC->QR)
       dq_accr = 2.2_RP * q(I_QC) * q(I_QR)**0.875_RP

       ! Evaporation (QR->QV)
       vent_factor = 1.6_RP + 124.9_RP * ( dens * q(I_QR) )**0.2046_RP

       dq_evap = max( qsat-q(I_QV), 0.0_RP ) / qsat / dens  &
               * vent_factor * ( dens * q(I_QR) )**0.525_RP &
               / ( 5.4E5_RP + 2.55E8_RP / ( pres*qsat ) )

       ! tendency QV, QC, QR with limiter
       dqv = (                  dq_evap ) * dt
       dqc = ( -dq_auto-dq_accr         ) * dt
       dqr = (  dq_auto+dq_accr-dq_evap ) * dt
       correct = dqc + dqr

       dqc = max( -q(I_QC), dqc )
       dqr = max( -q(I_QR), dqr )
       dqv = dqv + correct - ( dqc + dqr )

       rhoe = rhoe - LHV00 * dqv

       !##### End Main #####


       ! update from 1D work

       call THERMODYN_rhot( rhot, & ! [OUT]
                            rhoe, & ! [IN]
                            q(:)  ) ! [IN]

       QTRC_t(k,i,j,I_QV) = dqv
       QTRC_t(k,i,j,I_QC) = dqc
       QTRC_t(k,i,j,I_QR) = dqr

       QTRC0(k,i,j,I_QV) = q(I_QV) + dqv
       QTRC0(k,i,j,I_QC) = q(I_QC) + dqc
       QTRC0(k,i,j,I_QR) = q(I_QR) + dqr

       RHOT0(k,i,j) = rhot

    enddo
    enddo
    enddo

    return
  end subroutine MP_kessler

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics (terminal velocity)
  !-----------------------------------------------------------------------------
  subroutine MP_kessler_vterm( &
       vterm, &
       DENS,  &
       QTRC   )
    implicit none

    real(RP), intent(inout) :: vterm(KA,IA,JA,QA)
    real(RP), intent(in)    :: DENS (KA,IA,JA)
    real(RP), intent(in)    :: QTRC (KA,IA,JA,QA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! only update QR
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       vterm(k,i,j,I_QR) = - 36.34_RP * ( DENS(k,i,j) * QTRC(k,i,j,I_QR) )**0.1364_RP * factor_vterm(k)
    enddo
    enddo
    enddo

    return
  end subroutine MP_kessler_vterm

end module mod_atmos_phy_mp 
