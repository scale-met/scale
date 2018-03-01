!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Aerosol Microphysics driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_ae_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_driver_tracer_setup
  public :: ATMOS_PHY_AE_driver_setup
  public :: ATMOS_PHY_AE_driver_resume
  public :: ATMOS_PHY_AE_driver_tendency

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
  subroutine ATMOS_PHY_AE_driver_tracer_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE!, &
!       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
       IC_MIX, &
       IC_SEA, &
       IC_DUS, &
       NSIZ, &
       NKAP, &
       ATMOS_PHY_AE_kajino13_NAME, &
       ATMOS_PHY_AE_kajino13_DESC, &
       ATMOS_PHY_AE_kajino13_UNIT, &
       AE_CTG, &
       GAS_CTG, &
       IG_CGAS, &
       IG_H2SO4, &
       N_ATR, &
       QAEE, &
       QAES
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_phy_ae, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    implicit none

    integer, allocatable :: aero_idx(:,:,:,:)
    integer :: n_kap_max, n_siz_max, ncat_max
    integer :: NASIZ(3), NAKAP(3)
    character(len=H_SHORT) :: attribute, catego, aunit

    integer :: QS

    NAMELIST / PARAM_TRACER_KAJINO13 / &
       AE_CTG, &
       NASIZ,  &
       NAKAP

    integer :: m, ierr, ik, ic, ia0, is0

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_AE] / Origin[SCALE-RM]'

    ! note: tentatively, aerosol module should be called at all time. we need dummy subprogram.
    !if ( ATMOS_sw_phy_ae ) then

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          write(*,*) '### aerosol type(', ATMOS_PHY_AE_TYPE, '). is not recommended in current version!'

          ncat_max = max( IC_MIX, IC_SEA, IC_DUS )
   
          NASIZ(:) = 64
          NAKAP(:) = 1
   
          rewind(IO_FID_CONF)
          read(IO_FID_CONF,nml=PARAM_TRACER_KAJINO13,iostat=ierr)

          if( ierr < 0 ) then !--- missing
             if( IO_L ) write(IO_FID_LOG,*)  '*** Not found namelist. Default used.'
          elseif( ierr > 0 ) then !--- fatal error
             write(*,*) 'xxx Not appropriate names in namelist PARAM_TRACER_KAJINO13, Check!'
             call PRC_MPIstop
          end if
   
          if( IO_NML ) write(IO_FID_NML,nml=PARAM_TRACER_KAJINO13)
   
          if( AE_CTG > ncat_max ) then
             write(*,*) 'xxx AE_CTG should be smaller than', ncat_max+1, 'stop'
             call PRC_MPIstop
          endif
   
          allocate( NSIZ(AE_CTG) )
          allocate( NKAP(AE_CTG) )
   
          NKAP(1:AE_CTG) = NAKAP(1:AE_CTG)
          NSIZ(1:AE_CTG) = NASIZ(1:AE_CTG)
   
          if( maxval( NKAP ) /= 1 .OR. minval( NKAP ) /= 1 ) then
             write(*,*) 'xxx NKAP(:) /= 1 is not supported now, Stop!'
             call PRC_MPIstop
          end if
   
   !       do ia0 = 1, N_ATR
          do ic = 1, AE_CTG
          do ik = 1, NKAP(ic)
          do is0 = 1, NSIZ(ic)
             QA_AE   = QA_AE + N_ATR
          enddo
          enddo
          enddo
          QA_AE = QA_AE + GAS_CTG
   
          allocate( ATMOS_PHY_AE_kajino13_NAME(QA_AE) )
          allocate( ATMOS_PHY_AE_kajino13_DESC(QA_AE) )
          allocate( ATMOS_PHY_AE_kajino13_UNIT(QA_AE) )
   
          n_siz_max = 0
          n_kap_max = 0
          do ic = 1, AE_CTG
            n_siz_max = max(n_siz_max, NSIZ(ic))
            n_kap_max = max(n_kap_max, NKAP(ic))
          enddo
   
          allocate( aero_idx(N_ATR,AE_CTG,n_kap_max,n_siz_max) )
          m = 0
   !       do ia0 = 1, N_ATR
          do ic = 1, AE_CTG
          do ik = 1, NKAP(ic)
          do is0 = 1, NSIZ(ic)
          do ia0 = 1, N_ATR
            m = m+1
            aero_idx(ia0,ic,ik,is0) = m
          enddo
          enddo
          enddo
          enddo
   
          !-----------------------------------------------------------------------------
          !
          !++ calculate each category and aerosol
          !
          !-----------------------------------------------------------------------------
          ic = QA_AE-GAS_CTG+IG_H2SO4
          write(ATMOS_PHY_AE_kajino13_UNIT(ic),'(a)')  'kg/kg'
          ic = QA_AE-GAS_CTG+IG_CGAS
          write(ATMOS_PHY_AE_kajino13_UNIT(ic),'(a)')  'kg/kg'
   
   !       do ia0 = 1, N_ATR
   !         if( ia0 == 1 ) then
   !            write(attribute,'(a)') "Number"
   !            write(aunit,'(a)') "num/kg"
   !         elseif( ia0 == 2 ) then
   !            write(attribute,'(a)') "Section"
   !            write(aunit,'(a)') "m2/kg"
   !         elseif( ia0 == 3 ) then
   !            write(attribute,'(a)') "Volume"
   !            write(aunit,'(a)') "m3/kg"
   !         elseif( ia0 == 4 ) then
   !            write(attribute,'(a)') "Mass"
   !            write(aunit,'(a)') "kg/kg"
   !         elseif( ia0 == 5 ) then
   !            write(attribute,'(a)') "kpXmass"
   !            write(aunit,'(a)') "kg/kg"
   !         endif
          do ic = 1, AE_CTG       !aerosol category
          do ik = 1, NKAP(ic)   !kappa bin
          do is0 = 1, NSIZ(ic)
          do ia0 = 1, N_ATR
            if( ia0 == 1 ) then
               write(attribute,'(a)') "Number"
               write(aunit,'(a)') "num/kg"
            elseif( ia0 == 2 ) then
               write(attribute,'(a)') "Section"
               write(aunit,'(a)') "m2/kg"
            elseif( ia0 == 3 ) then
               write(attribute,'(a)') "Volume"
               write(aunit,'(a)') "m3/kg"
            elseif( ia0 == 4 ) then
               write(attribute,'(a)') "Mass"
               write(aunit,'(a)') "kg/kg"
            elseif( ia0 == 5 ) then
               write(attribute,'(a)') "kXm"
               write(aunit,'(a)') "kg/kg"
            endif
            if( ic == IC_MIX ) then
               write(catego,'(a)') "Sulf_"
            elseif( ic == IC_SEA ) then
               write(catego,'(a)') "Salt_"
            elseif( ic == IC_DUS ) then
               write(catego,'(a)') "Dust_"
            endif
            write(ATMOS_PHY_AE_kajino13_UNIT(aero_idx(ia0,ic,ik,is0)),'(a)')  trim(aunit)
            write(ATMOS_PHY_AE_kajino13_NAME(aero_idx(ia0,ic,ik,is0)),'(a,a,i0)') trim(catego), trim(attribute), is0
            write(ATMOS_PHY_AE_kajino13_DESC(aero_idx(ia0,ic,ik,is0)),'(a,a,a,i0)') trim(attribute), ' mixing radio of ', trim(catego), is0
          enddo
          enddo
          enddo
          enddo
          ic = QA_AE-GAS_CTG+IG_H2SO4
          write(ATMOS_PHY_AE_kajino13_NAME(ic),'(a)') 'H2SO4_Gas'
          ic = QA_AE-GAS_CTG+IG_CGAS
          write(ATMOS_PHY_AE_kajino13_NAME(ic),'(a)') 'Condensable_GAS'
   
          ic = QA_AE-GAS_CTG+IG_H2SO4
          write(ATMOS_PHY_AE_kajino13_DESC(ic),'(a)') 'Mixing ratio of H2SO4 Gas'
          ic = QA_AE-GAS_CTG+IG_CGAS
          write(ATMOS_PHY_AE_kajino13_DESC(ic),'(a)') 'Mixing ratio of Condensable GAS'
   
          deallocate(aero_idx)
   
          call TRACER_regist( QS,                         & ! [OUT]
                              QA_AE,                      & ! [IN]
                              ATMOS_PHY_AE_kajino13_NAME, & ! [IN]
                              ATMOS_PHY_AE_kajino13_DESC, & ! [IN]
                              ATMOS_PHY_AE_kajino13_UNIT  ) ! [IN]
   
          QA   = QA_AE
          QAES = QS
          QAEE = QS + QA_AE - 1

          QA_AE = QA_AE
          QS_AE = QAES
          QE_AE = QAEE
       case default
          write(*,*) 'xxx invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_MPIstop
       end select

    !else
       !if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    !endif

    return
  end subroutine ATMOS_PHY_AE_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_driver_setup
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE!, &
!       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
        ATMOS_PHY_AE_kajino13_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_AE] / Origin[SCALE-RM]'

    ! note: tentatively, aerosol module should be called at all time. we need dummy subprogram.
!    if ( ATMOS_sw_phy_ae ) then

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_setup
       case default
          write(*,*) 'xxx invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_MPIstop
       end select

    !else
       !if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    !endif

    return
  end subroutine ATMOS_PHY_AE_driver_setup


  !-----------------------------------------------------------------------------
  !> Resume
  subroutine ATMOS_PHY_AE_driver_resume
!     use mod_atmos_admin, only: &
!        ATMOS_sw_phy_ae
    implicit none

    ! note: tentatively, aerosol module should be called at all time. we need dummy subprogram.
!    if ( ATMOS_sw_phy_ae ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Aerosol', 1)
       call ATMOS_PHY_AE_driver_tendency( update_flag = .true. )
       call PROF_rapend  ('ATM_Aerosol', 1)

!    endif

    return
  end subroutine ATMOS_PHY_AE_driver_resume

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_AE_driver_tendency( update_flag )
    use scale_time, only: &
       dt_AE => TIME_DTSEC_ATMOS_PHY_AE
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_phy_ae, only: &
       ATMOS_PHY_AE, &
       QA_AE, &
       QS_AE, &
       QE_AE
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       MOMZ   => MOMZ_av, &
       MOMX   => MOMX_av, &
       MOMY   => MOMY_av, &
       RHOT   => RHOT_av, &
       QTRC   => QTRC_av, &
       RHOQ_t => RHOQ_tp, &
       QDRY, &
       PRES, &
       TEMP
    use mod_atmos_phy_ae_vars, only: &
       RHOQ_t_AE => ATMOS_PHY_AE_RHOQ_t, &
       CCN       => ATMOS_PHY_AE_CCN,   &
       CCN_t     => ATMOS_PHY_AE_CCN_t, &
       AE_EMIT   => ATMOS_PHY_AE_EMIT
    use mod_atmos_phy_mp_vars, only: &
       EVAPORATE => ATMOS_PHY_MP_EVAPORATE
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_tendency, &
       QAES, &
       QAEE
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: CN(KA,IA,JA)
    real(RP) :: NREG(KA,IA,JA)

    real(RP) :: total ! dummy

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

!OCL XFILL
       CCN      (:,:,:)   = 0.0_RP ! reset
!OCL XFILL
       RHOQ_t_AE(:,:,:,:) = 0.0_RP ! reset

       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          NREG(k,i,j) = EVAPORATE(k,i,j) * dt_AE
       enddo
       enddo
       enddo

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_KAJINO13_TENDENCY( KA, KS, KE,                 & ! [IN]
                                               IA, IS, IE,                 & ! [IN]
                                               JA, JS, JE,                 & ! [IN]
                                               QA_AE,                      & ! [IN]
                                               RHOT     (:,:,:),           & ! [IN]
                                               TEMP     (:,:,:),           & ! [IN]
                                               PRES     (:,:,:),           & ! [IN]
                                               QDRY     (:,:,:),           & ! [IN]
                                               NREG     (:,:,:),           & ! [IN]
                                               DENS     (:,:,:),           & ! [INOUT]
                                               MOMZ     (:,:,:),           & ! [INOUT]
                                               MOMX     (:,:,:),           & ! [INOUT]
                                               MOMY     (:,:,:),           & ! [INOUT]
                                               QTRC     (:,:,:,QAES:QAEE), & ! [INOUT]
                                               RHOQ_t_AE(:,:,:,QAES:QAEE), & ! [INOUT]
                                               AE_EMIT  (:,:,:,QAES:QAEE), & ! [INOUT]
                                               CN       (:,:,:),           & ! [OUT]
                                               CCN      (:,:,:) )            ! [OUT]
       case default
          write(*,*) 'xxx invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_MPIstop
       end select

       CCN_t(:,:,:) = CCN(:,:,:) / dt_AE

       call FILE_HISTORY_in( CN (:,:,:)*1.E-6_RP, 'CN',  'condensation nucrei',       'num/cc' )
       call FILE_HISTORY_in( CCN(:,:,:)*1.E-6_RP, 'CCN', 'cloud condensation nucrei', 'num/cc' )

    endif

    do iq = QS_AE, QE_AE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_AE(k,i,j,iq)
       enddo
       enddo
       enddo
    enddo

    if ( STATISTICS_checktotal ) then
       do iq = QS_AE, QE_AE
          call STAT_total( total, RHOQ_t_AE(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_AE' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_AE_driver_tendency

end module mod_atmos_phy_ae_driver
