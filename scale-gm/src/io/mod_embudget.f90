!-------------------------------------------------------------------------------
!> Module budget monitoring
!!
!! @par Description
!!          This module is for monitoring the energy/mass budget
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_embudget
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: embudget_setup
  public :: embudget_monitor

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private :: MNT_ON   = .false.
  integer,  private :: MNT_INTV = 1
  integer,  private :: MNT_m_fid
  integer,  private :: MNT_e_fid

  real(RP), private :: rhoqd_sum_old   = 0.0_RP
  real(RP), private :: rhoqv_sum_old   = 0.0_RP
  real(RP), private :: rhoql_sum_old   = 0.0_RP
  real(RP), private :: rhoqi_sum_old   = 0.0_RP
  real(RP), private :: rhoqt_sum_old   = 0.0_RP
  real(RP), private :: rhophi_sum_old  = 0.0_RP
  real(RP), private :: rhoein_sum_old  = 0.0_RP
  real(RP), private :: rhokin_sum_old  = 0.0_RP
  real(RP), private :: rhoetot_sum_old = 0.0_RP

  real(RP), private :: Mass_budget_factor
  real(RP), private :: Energy_budget_factor

  logical,  private :: first = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine embudget_setup
    use scale_process, only: &
       PRC_IsMaster, &
       PRC_MPIstop
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI      => CONST_PI
    use mod_time, only: &
       TIME_DTL
    implicit none

    namelist / EMBUDGETPARAM / &
       MNT_INTV, &
       MNT_ON

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[embudget]/Category[nhm share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=EMBUDGETPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** EMBUDGETPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist EMBUDGETPARAM. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist EMBUDGETPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=EMBUDGETPARAM)

    if(.not.MNT_ON) return

    Mass_budget_factor   = 1.0_RP / ( TIME_DTL * real(MNT_INTV,kind=RP) * 4.0_RP * PI * RADIUS * RADIUS ) ! [kg/step] -> [kg/m2/s]
    Energy_budget_factor = 1.0_RP / ( TIME_DTL * real(MNT_INTV,kind=RP) * 4.0_RP * PI * RADIUS * RADIUS ) ! [J /step] -> [W/m2]
    if( IO_L ) write(IO_FID_LOG,*) "Mass_budget_factor   = ", Mass_budget_factor
    if( IO_L ) write(IO_FID_LOG,*) "Energy_budget_factor = ", Energy_budget_factor

    ! open budget.info file
    if ( PRC_IsMaster ) then
       MNT_m_fid  = IO_get_available_fid()
       open( unit   = MNT_m_fid,          &
             file   = 'MASS_BUDGET.info', &
             form   = 'formatted',        &
             status = 'unknown'           )

       MNT_e_fid = IO_get_available_fid()
       open( unit   = MNT_e_fid,            &
             file   = 'ENERGY_BUDGET.info', &
             form   = 'formatted',          &
             status = 'unknown'             )
    endif

    call diagnose_energy_mass

    return
  end subroutine embudget_setup

  !-----------------------------------------------------------------------------
  subroutine embudget_monitor
    use mod_time, only: &
       TIME_CSTEP
    implicit none
    !---------------------------------------------------------------------------

    if( .NOT. MNT_ON ) return

    if ( mod(TIME_CSTEP-1,MNT_INTV) == 0 ) then
       call diagnose_energy_mass
    endif

    return
  end subroutine embudget_monitor

  !-----------------------------------------------------------------------------
  subroutine diagnose_energy_mass
    use scale_process, only: &
       PRC_IsMaster
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI,     &
       CVdry  => CONST_CVdry
    use scale_atmos_hydrometeor, only: &
       LHV, &
       LHF
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_vmtr, only: &
       VMTR_getIJ_RGSGAM2, &
       VMTR_getIJ_PHI
    use mod_time, only: &
       TIME_CSTEP, &
       TIME_DTL
    use mod_gm_statistics, only: &
       GTL_global_sum, &
       GTL_global_sum_srf
    use mod_runconf, only: &
       TRC_vmax, &
       NQW_STR,  &
       NQW_END,  &
       I_QV,     &
       I_QC,     &
       I_QR,     &
       I_QI,     &
       I_QS,     &
       I_QG,     &
       CVW
    use mod_prgvar, only: &
       prgvar_get_withdiag
    use mod_cnvvar, only: &
       cnvvar_rhogkin
    use mod_thrmdyn, only: &
       THRMDYN_qd
    implicit none

    real(RP) :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq    (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: qd    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qd_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tmp   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhoq_sum    (TRC_vmax)
    real(RP) :: rhoein_q_sum(TRC_vmax)
    real(RP) :: rhoein_qd_sum

    real(RP) :: rhoqd_sum
    real(RP) :: rhoqv_sum
    real(RP) :: rhoql_sum
    real(RP) :: rhoqi_sum
    real(RP) :: rhoqt_sum
    real(RP) :: rhophi_sum
    real(RP) :: rhoein_sum
    real(RP) :: rhokin_sum
    real(RP) :: rhoetot_sum

    real(RP) :: rhoqd_sum_diff
    real(RP) :: rhoqv_sum_diff
    real(RP) :: rhoql_sum_diff
    real(RP) :: rhoqi_sum_diff
    real(RP) :: rhoqt_sum_diff
    real(RP) :: rhophi_sum_diff
    real(RP) :: rhoein_sum_diff
    real(RP) :: rhokin_sum_diff
    real(RP) :: rhoetot_sum_diff

    real(RP) :: VMTR_RGSGAM2   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_RGSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: VMTR_PHI       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_PHI_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: nq
    !---------------------------------------------------------------------------

    call VMTR_getIJ_RGSGAM2( VMTR_RGSGAM2, VMTR_RGSGAM2_pl )
    call VMTR_getIJ_PHI    ( VMTR_PHI,     VMTR_PHI_pl     )

    call prgvar_get_withdiag( rhog,   rhog_pl,   & ! [OUT]
                              rhogvx, rhogvx_pl, & ! [OUT]
                              rhogvy, rhogvy_pl, & ! [OUT]
                              rhogvz, rhogvz_pl, & ! [OUT]
                              rhogw,  rhogw_pl,  & ! [OUT]
                              rhoge,  rhoge_pl,  & ! [OUT]
                              rhogq,  rhogq_pl,  & ! [OUT]
                              rho,    rho_pl,    & ! [OUT]
                              pre,    pre_pl,    & ! [OUT]
                              tem,    tem_pl,    & ! [OUT]
                              vx,     vx_pl,     & ! [OUT]
                              vy,     vy_pl,     & ! [OUT]
                              vz,     vz_pl,     & ! [OUT]
                              w,      w_pl,      & ! [OUT]
                              q,      q_pl       ) ! [OUT]

    call THRMDYN_qd( ADM_gall,    & ! [IN]
                     ADM_kall,    & ! [IN]
                     ADM_lall,    & ! [IN]
                     q (:,:,:,:), & ! [IN]
                     qd(:,:,:)    ) ! [OUT]

    if ( ADM_have_pl ) then
       call THRMDYN_qd( ADM_gall_pl,    & ! [IN]
                        ADM_kall,       & ! [IN]
                        ADM_lall_pl,    & ! [IN]
                        q_pl (:,:,:,:), & ! [IN]
                        qd_pl(:,:,:)    ) ! [OUT]
    endif

    !----- Mass budget

    !--- total mass (dry air)
    tmp(:,:,:) = rho(:,:,:) * qd(:,:,:)
    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:) * qd_pl(:,:,:)
    endif
    rhoqd_sum = GTL_global_sum( tmp, tmp_pl )

    !--- total mass (each water category)
    do nq = NQW_STR, NQW_END
       tmp(:,:,:) = rho(:,:,:) * q(:,:,:,nq)
       if ( ADM_have_pl ) then
          tmp_pl(:,:,:) = rho_pl(:,:,:) * q_pl(:,:,:,nq)
       endif
       rhoq_sum(nq) = GTL_global_sum( tmp, tmp_pl )
    enddo

    !--- total mass (total/vapor/liquid/soild water)
    rhoqt_sum = 0.0_RP
    rhoqv_sum = 0.0_RP
    rhoql_sum = 0.0_RP
    rhoqi_sum = 0.0_RP
    do nq = NQW_STR, NQW_END
       rhoqt_sum = rhoqt_sum + rhoq_sum(nq)

       if    ( nq == I_QV ) then
          rhoqv_sum = rhoqv_sum + rhoq_sum(nq)
       elseif( nq == I_QC .OR. nq == I_QR ) then
          rhoql_sum = rhoql_sum + rhoq_sum(nq)
       elseif( nq == I_QI .OR. nq == I_QS  .OR. nq == I_QG ) then
          rhoqi_sum = rhoqi_sum + rhoq_sum(nq)
       endif
    enddo


    !----- Energy budget

    !--- potential energy
    tmp(:,:,:) = rho(:,:,:) * VMTR_PHI(:,:,:)
    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = rho_pl(:,:,:) * VMTR_PHI_pl(:,:,:)
    endif
    rhophi_sum = GTL_global_sum( tmp, tmp_pl )

    !--- internal energy (dry air)
    tmp = rho * qd * CVdry * tem
    if ( ADM_have_pl ) then
       tmp_pl = rho_pl * qd_pl * CVdry * tem_pl
    endif
    rhoein_qd_sum = GTL_global_sum( tmp, tmp_pl )

    !--- internal energy (each water category)
    do nq = NQW_STR,NQW_END

       tmp(:,:,:) = rho(:,:,:) * q(:,:,:,nq) * CVW(nq) * tem(:,:,:)
       if    ( nq == I_QV ) then
          tmp(:,:,:) = tmp(:,:,:) + rho(:,:,:) * q(:,:,:,nq) * LHV ! correct latent heat
       elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
          tmp(:,:,:) = tmp(:,:,:) - rho(:,:,:) * q(:,:,:,nq) * LHF ! correct latent heat
       endif

       if ( ADM_have_pl ) then
          tmp_pl(:,:,:) = rho_pl(:,:,:) * q_pl(:,:,:,nq) * CVW(nq) * tem_pl(:,:,:)
          if    ( nq == I_QV ) then
             tmp_pl(:,:,:) = tmp_pl(:,:,:) + rho_pl(:,:,:) * q_pl(:,:,:,nq) * LHV ! correct latent heat
          elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then
             tmp_pl(:,:,:) = tmp_pl(:,:,:) - rho_pl(:,:,:) * q_pl(:,:,:,nq) * LHF ! correct latent heat
          endif
       endif

       rhoein_q_sum(nq) = GTL_global_sum( tmp, tmp_pl )
    enddo

    !--- internal energy (total)
    rhoein_sum = rhoein_qd_sum
    do nq = NQW_STR,NQW_END
       rhoein_sum = rhoein_sum + rhoein_q_sum(nq)
    enddo

    !--- kinetic energy
    call cnvvar_rhogkin( rhog,   rhog_pl,   & ! [IN]
                         rhogvx, rhogvx_pl, & ! [IN]
                         rhogvy, rhogvy_pl, & ! [IN]
                         rhogvz, rhogvz_pl, & ! [IN]
                         rhogw,  rhogw_pl,  & ! [IN]
                         tmp,    tmp_pl     ) ! [OUT]

    tmp(:,:,:) = tmp(:,:,:) * VMTR_RGSGAM2(:,:,:)
    if ( ADM_have_pl ) then
       tmp_pl(:,:,:) = tmp_pl(:,:,:) * VMTR_RGSGAM2_pl(:,:,:)
    endif
    rhokin_sum = GTL_global_sum( tmp, tmp_pl )

    !--- total energy
    rhoetot_sum = rhoein_sum + rhophi_sum + rhokin_sum



    !##### File OUTPUT #####

    if ( first ) then
       ! [kg/m2], absolute value
       rhoqd_sum_diff   = rhoqd_sum   / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhoqv_sum_diff   = rhoqv_sum   / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhoql_sum_diff   = rhoql_sum   / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhoqi_sum_diff   = rhoqi_sum   / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhoqt_sum_diff   = rhoqt_sum   / ( 4.0_RP * PI * RADIUS * RADIUS )
       ! [J/m2], absolute value
       rhophi_sum_diff  = rhophi_sum  / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhoein_sum_diff  = rhoein_sum  / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhokin_sum_diff  = rhokin_sum  / ( 4.0_RP * PI * RADIUS * RADIUS )
       rhoetot_sum_diff = rhoetot_sum / ( 4.0_RP * PI * RADIUS * RADIUS )
    else
       ! [kg/m2/s], difference from previous step
       rhoqd_sum_diff   = ( rhoqd_sum   - rhoqd_sum_old   ) * Mass_budget_factor
       rhoqv_sum_diff   = ( rhoqv_sum   - rhoqv_sum_old   ) * Mass_budget_factor
       rhoql_sum_diff   = ( rhoql_sum   - rhoql_sum_old   ) * Mass_budget_factor
       rhoqi_sum_diff   = ( rhoqi_sum   - rhoqi_sum_old   ) * Mass_budget_factor
       rhoqt_sum_diff   = ( rhoqt_sum   - rhoqt_sum_old   ) * Mass_budget_factor
       ! [W/m2], difference from previous step
       rhophi_sum_diff  = ( rhophi_sum  - rhophi_sum_old  ) * Energy_budget_factor
       rhoein_sum_diff  = ( rhoein_sum  - rhoein_sum_old  ) * Energy_budget_factor
       rhokin_sum_diff  = ( rhokin_sum  - rhokin_sum_old  ) * Energy_budget_factor
       rhoetot_sum_diff = ( rhoetot_sum - rhoetot_sum_old ) * Energy_budget_factor
    endif

    if ( PRC_IsMaster ) then
       if ( first ) then
          write(MNT_m_fid,'(A6)' ,advance='no') '#STEP'
          write(MNT_m_fid,'(A16)',advance='no') 'Day'
          write(MNT_m_fid,'(A16)',advance='no') 'dry air mass '
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(g)'
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(l)'
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(s)'
          write(MNT_m_fid,'(A16)',advance='no') 'water mass(t)'
          write(MNT_m_fid,*)

          write(MNT_e_fid,'(A6)' ,advance='no') '#STEP'
          write(MNT_e_fid,'(A16)',advance='no') 'Day'
          write(MNT_e_fid,'(A16)',advance='no') 'energy(pot)'
          write(MNT_e_fid,'(A16)',advance='no') 'energy(int)'
          write(MNT_e_fid,'(A16)',advance='no') 'energy(kin)'
          write(MNT_e_fid,'(A16)',advance='no') 'energy(tot)'
          write(MNT_e_fid,*)
       endif

       ! mass budget
       write(MNT_m_fid,'(I6)'    ,advance='no') TIME_CSTEP
       write(MNT_m_fid,'(ES16.8)',advance='no') TIME_CSTEP * TIME_DTL / 86400.0_RP
       write(MNT_m_fid,'(ES16.8)',advance='no') rhoqd_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') rhoqv_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') rhoql_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') rhoqi_sum_diff
       write(MNT_m_fid,'(ES16.8)',advance='no') rhoqt_sum_diff
       write(MNT_m_fid,*)

       ! energy budget
       write(MNT_e_fid,'(I6)'    ,advance='no') TIME_CSTEP
       write(MNT_e_fid,'(ES16.8)',advance='no') TIME_CSTEP * TIME_DTL / 86400.0_RP
       write(MNT_e_fid,'(ES16.8)',advance='no') rhophi_sum_diff
       write(MNT_e_fid,'(ES16.8)',advance='no') rhoein_sum_diff
       write(MNT_e_fid,'(ES16.8)',advance='no') rhokin_sum_diff
       write(MNT_e_fid,'(ES16.8)',advance='no') rhoetot_sum_diff
       write(MNT_e_fid,*)
    endif

    rhoqd_sum_old   = rhoqd_sum
    rhoqv_sum_old   = rhoqv_sum
    rhoql_sum_old   = rhoql_sum
    rhoqi_sum_old   = rhoqi_sum
    rhoqt_sum_old   = rhoqt_sum
    rhophi_sum_old  = rhophi_sum
    rhoein_sum_old  = rhoein_sum
    rhokin_sum_old  = rhokin_sum
    rhoetot_sum_old = rhoetot_sum

    if ( first ) then
       first = .false.
    endif

    return
  end subroutine diagnose_energy_mass

end module mod_embudget
