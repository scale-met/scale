!-------------------------------------------------------------------------------
!> Module basic state
!!
!! @par Description
!!          This module is for the set of basic state for non-hydrostatic model
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_bsstate
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !
  !--- density
  real(RP),allocatable, public :: rho_bs(:,:,:)
  real(RP),allocatable, public :: rho_bs_pl(:,:,:)
  !
  !--- pressure
  real(RP),allocatable, public :: pre_bs(:,:,:)
  real(RP),allocatable, public :: pre_bs_pl(:,:,:)
  !
  !--- temperature
  real(RP),allocatable, public :: tem_bs(:,:,:)
  real(RP),allocatable, public :: tem_bs_pl(:,:,:)
  !
  !--- pot temperature
  real(RP),allocatable, public :: th_bs(:,:,:)
  real(RP),allocatable, public :: th_bs_pl(:,:,:)
  !
  !--- water vap.
  real(RP),allocatable, public :: qv_bs(:,:,:)
  real(RP),allocatable, public :: qv_bs_pl(:,:,:)
  !
  !--- geo-potential ( g X z )
  real(RP),allocatable, private :: phi(:,:,:)
  real(RP),allocatable, private :: phi_pl(:,:,:)
  !
  !--- Basic state type
  character(len=H_SHORT), public :: ref_type = 'NOBASE'
  !                                  ='TEM': temerature is given.
  !                                  ='TH' : potential temperature is given.
  !                                  ='NOBASE' : no basic state
  !                                  ='INPUT'  : input
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: bsstate_setup
  public :: bsstate_input_ref
  public :: bsstate_output_ref
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  !--- reference pressure at the ground
  real(RP), private :: pre_g = 101325.0_RP
  !
  !--- reference temperature at the ground
  real(RP), private :: tem_g = 300.0_RP
  !
  !--- reference pot. temperature at the ground
  real(RP), private :: th_g = 300.0_RP
  !
  !--- reference density at the ground ( calculated by using pre_g & tem_g )
  real(RP), private :: rho_g
  !
  !--- reference Brunt-Vaisala frequency ( used if ref_type = 'TH'. )
  real(RP), private :: BV_freq = 0.0_RP
  !
  !--- lapse rate ( used if ref_type = 'TEM'. )
  real(RP), private :: TGAMMA = 0.0_RP
  !
  !--- lower boundary of constant (potential) temperature
  real(RP), private :: ZT = 0.0_RP
  !
  !--- geopotential at the ground
  real(RP), parameter, public :: PHI0=0.0_RP
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !--- reference phi
  real(RP), allocatable, private :: phi_ref(:)
  !
  !--- reference density
  real(RP), allocatable, private :: rho_ref(:)
  !
  !--- reference pressure
  real(RP), allocatable, private :: pre_ref(:)
  !
  !--- reference temperature
  real(RP), allocatable, private :: tem_ref(:)
  !
  !--- water vapor
  real(RP), allocatable, private :: qv_ref(:)
  !
  !--- reference potential temperature
  real(RP), allocatable, private :: th_ref(:)
  !
  character(len=H_LONG), private :: ref_fname = 'ref.dat'
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: set_referencestate
  private :: set_basicstate
  private :: output_info
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine bsstate_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_kall,      &
       ADM_gall_pl,   &
       ADM_gall
    implicit none

    namelist / BSSTATEPARAM / &
         ref_type, & !--- type of basic state
         ZT,       & !--- if z>ZT, equi-temperature
         pre_g,    & !--- reference pressure
         tem_g,    & !--- reference temperature
         TGAMMA,   & !--- lapse rate when ref_type='TEM'.
         th_g,     & !--- reference potential temp. for when ref_type='TH'
         BV_freq,  & !--- Vaisala freq when ref_type='TH'
         ref_fname   !--- in/output base name if necessary.

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[basic state]/Category[nhm share]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=BSSTATEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** BSSTATEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist BSSTATEPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist BSSTATEPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=BSSTATEPARAM)

    !--- allocation of reference variables
    allocate(phi_ref(ADM_kall))
    allocate(rho_ref(ADM_kall))
    allocate(pre_ref(ADM_kall))
    allocate(tem_ref(ADM_kall))
    allocate(th_ref(ADM_kall))
    allocate(qv_ref(ADM_kall))
    !
    ! add by kgoto
    ! initialize
    phi_ref=0.0_RP
    rho_ref=0.0_RP
    pre_ref=0.0_RP
    tem_ref=0.0_RP
    th_ref=0.0_RP
    !--- allocation of the basic variables
    allocate(rho_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(rho_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl))
    allocate(pre_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(pre_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl))
    allocate(tem_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(tem_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl))
    allocate(th_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(th_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl))
    allocate(qv_bs(ADM_gall,ADM_kall,ADM_lall))
    allocate(qv_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl))
    allocate(phi(ADM_gall,ADM_kall,ADM_lall))
    allocate(phi_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl))
    !
    ! add by kgoto
    ! initialize
    rho_bs    = 0.0_RP
    rho_bs_pl = 0.0_RP
    pre_bs    = 0.0_RP
    pre_bs_pl = 0.0_RP
    tem_bs    = 0.0_RP
    tem_bs_pl = 0.0_RP
    th_bs    = 0.0_RP
    th_bs_pl = 0.0_RP
    qv_bs    = 0.0_RP
    qv_bs_pl = 0.0_RP
    !
    if ( ref_type=='INPUT') then
       !
       !---- input reference state
       call bsstate_input_ref(ref_fname)
       !
       !--- calculation of basic state
       call set_basicstate
    else !--- other type
       !
       !--- calculation of reference state
       call set_referencestate
       !
       !--- calculation of basic state
       call set_basicstate
    endif
    !
    !--- output the information
    call output_info

    return
  end subroutine bsstate_setup

  !-----------------------------------------------------------------------------
  subroutine output_info
    use mod_adm, only: &
       ADM_LOG_FID, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) 'Msg : Sub[nhm_bs_output_info]/Mod[basicstate]'
    write(ADM_LOG_FID,*) ' --- Basic state type         : ', trim(ref_type)
    write(ADM_LOG_FID,*) ' --- Reference pressure       : ', pre_g
    write(ADM_LOG_FID,*) ' --- Reference temperature    : ', tem_g
    write(ADM_LOG_FID,*) ' --- Reference density        : ', rho_g
    write(ADM_LOG_FID,*) ' --- Vaisala frequency        : ', BV_freq
    write(ADM_LOG_FID,*) ' --- Lapse rate of temperature: ', TGAMMA
    write(ADM_LOG_FID,*) ' --- Effective height         : ', ZT

    write(ADM_LOG_FID,*) '-------------------------------------------------------'
    write(ADM_LOG_FID,*) 'Level   Density  Pressure     Temp. Pot. Tem.        qv'
    do k=1,ADM_kall
       write(ADM_LOG_FID,'(I4,F12.4,F10.2,F10.2,F10.2,F10.7)') &
                         k, rho_ref(k), pre_ref(k), tem_ref(k),th_ref(k), qv_ref(k)
       if (  k == ADM_kmin-1 ) write(ADM_LOG_FID,*) '-------------------------------------------------------'
       if (  k == ADM_kmax   ) write(ADM_LOG_FID,*) '-------------------------------------------------------'
    enddo

    return
  end subroutine output_info

  !-----------------------------------------------------------------------------
  subroutine set_referencestate
    use scale_const, only: &
       CONST_GRAV,  &
       CONST_Rdry,  &
       CONST_Rvap,  &
       CONST_CPdry, &
       CONST_PRE00
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_kall,       &
       ADM_kmin,       &
       ADM_kmax
    use mod_grd, only :  &
         GRD_gz,         &
         GRD_dgz,        &
         GRD_gzh,        &
         GRD_afac,       &
         GRD_bfac
    implicit none

    real(RP) :: dpre_ref_k
    real(RP) :: pre_s, rho_s, total_mass0, total_mass, mass_diff_ratio

    real(RP) :: kappa

    integer :: k
    !---------------------------------------------------------------------------

    kappa = CONST_Rdry / CONST_CPdry

    !--- calculation of reference geopotential
    do k=1,ADM_kall
       phi_ref(k) = CONST_GRAV * GRD_gz(k)
    enddo
    !
    if ( ref_type=='NOBASE') then
       !
       phi_ref = 0.0_RP
       rho_ref = 0.0_RP
       pre_ref = 0.0_RP
       tem_ref = 0.0_RP
       th_ref  = 0.0_RP
       qv_ref  = 0.0_RP
       !
       pre_bs    = 0.0_RP
       pre_bs_pl = 0.0_RP
       tem_bs    = 0.0_RP
       tem_bs_pl = 0.0_RP
       th_bs    = 0.0_RP
       th_bs_pl = 0.0_RP
       qv_bs    = 0.0_RP
       qv_bs_pl = 0.0_RP
       rho_bs    = 0.0_RP
       rho_bs_pl = 0.0_RP
       !
       ! 04/12/25 M.Satoh add
       if ( tem_g /= 0.0_RP ) then
          rho_g = pre_g / CONST_Rdry / tem_g
       else
          rho_g = 0.0_RP
       endif
       !
    else if ( ref_type == 'TEM' ) then
       qv_ref = 0.0_RP
       !
       !---  calculation of reference temperature
       do k=1,ADM_kall
          if ( GRD_gz(k) <=ZT ) then
             tem_ref(k) = tem_g - TGAMMA * GRD_gz(k)
          else
             tem_ref(k) = tem_g - TGAMMA * ZT
          endif
       enddo
       !
       !--- calculation of density at the surface
       rho_g = pre_g / CONST_Rdry / tem_g
       !
       !--- calculation of reference pressure and density
       !--- just below the ground level
       pre_ref(ADM_kmin-1)                                &
            = pre_g                                        &
            + 0.5_RP                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_g
       rho_ref(ADM_kmin-1)        &
            = pre_ref(ADM_kmin-1) &
            / CONST_Rdry           &
            / tem_ref(ADM_kmin-1)
       !
       !--- Reference pressure and density at the first level
       pre_ref(ADM_kmin)                                        &
            = pre_ref(ADM_kmin-1)                               &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_g
       rho_ref(ADM_kmin)        &
            = pre_ref(ADM_kmin) &
            / CONST_Rdry         &
            / tem_ref(ADM_kmin)
       !
       !--- Reference pressure and density at k level
       !---    In this caluculation, the hydrostatic balance equation
       !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer
       !---    level ( k-1/2 ). RHO is obtained by extrapolation from
       !---    the values at lower level.
       !---
       !--- Consistent way(?) in the scheme.
       !---
       !--- pre_ref(k) - pre_ref(k-1)
       !--- = - 0.5_RP * ( GRD_afac(k) * rho_ref(k)
       !---              +GRD_bfac(k) * rho_ref(k-1)  )
       !---   * ( phi_ref(k) - phi_ref(k-1) )
       !---
       !--- rho_ref(k)*CONST_Rdry*tem_ref(k)  - pre_ref(k-1)
       !--- = - 0.5_RP * ( GRD_afac(k) * rho_ref(k)
       !---              +GRD_bfac(k) * rho_ref(k-1)  )
       !---   * ( phi_ref(k) - phi_ref(k-1) )
       !
       !--- rho_ref(k)*( CONST_Rdry*tem_ref(k)
       !---            + 0.5_RP * GRD_afac(k)
       !---            * ( phi_ref(k) - phi_ref(k-1) )
       !---            )
       !--- = pre_ref(k-1)
       !---            - 0.5_RP * GRD_bfac(k) * rho_ref(k-1)
       !---            * ( phi_ref(k) - phi_ref(k-1) )
       !---
       do k = ADM_kmin+1, ADM_kmax+1
          rho_ref(k) = &
               ( pre_ref(k-1) &
               - 0.5_RP * GRD_bfac(k) * rho_ref(k-1) &
               * ( phi_ref(k) - phi_ref(k-1) )&
               ) / &
               ( CONST_Rdry*tem_ref(k) &
               + 0.5_RP * GRD_afac(k)  &
               * ( phi_ref(k) - phi_ref(k-1) )  &
               )
          pre_ref(k) = rho_ref(k) * CONST_Rdry * tem_ref(k)
       enddo
       !
       !--- calculation of reference potential temperature
       do k = 1, ADM_kall
          th_ref(k)                                        &
               = tem_ref(k)                                &
               * ( CONST_PRE00 / pre_ref(k) )**kappa
       enddo
       !
    else if ( ref_type == 'RHO' ) then
       qv_ref = 0.0_RP
       !
       !---  calculation of reference density
       total_mass0 = CONST_PRE00/CONST_GRAV
       pre_s = CONST_PRE00
       do
          do k = ADM_kmin-1, ADM_kmax+1
             rho_ref(k) = pre_s/CONST_Rdry/tem_g / exp(CONST_GRAV*GRD_gz(k)/CONST_Rdry/tem_g)
          enddo
          total_mass = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             total_mass = total_mass + GRD_dgz(k)*rho_ref(k)
          enddo
          mass_diff_ratio = total_mass0/total_mass
          if ( abs(mass_diff_ratio-1.0_RP)<1.E-8_RP) then
             exit
          else
             pre_s = pre_s * mass_diff_ratio
          endif
       enddo
       !
       !--- calculation of density at the surface
       rho_s = pre_s /CONST_Rdry / tem_g
       !
       !--- calculation of reference pressure and density
       !--- just below the ground level
       pre_ref(ADM_kmin-1)                                &
            = pre_s                                       &
            + 0.5_RP                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_s
       !
       !--- Reference pressure and density at the first level
       pre_ref(ADM_kmin)                                        &
            = pre_ref(ADM_kmin-1)                               &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_s
       !
       !--- Reference pressure and density at k level
       do k = ADM_kmin+1, ADM_kmax+1
          pre_ref(k) = pre_ref(k-1) &
               - 0.5_RP * ( GRD_afac(k) * rho_ref(k)     &
                          +GRD_bfac(k) * rho_ref(k-1)  )&
               * ( phi_ref(k) - phi_ref(k-1) )
       enddo
       !
       !--- calculation of reference temperature & pot temperature
       do k = 1, ADM_kall
          tem_ref(k) = pre_ref(k) / rho_ref(k) / CONST_Rdry
          th_ref(k)                                        &
               = tem_ref(k)                                &
               * ( CONST_PRE00 / pre_ref(k) )**kappa
       enddo
       !
    else if ( ref_type == 'TH' ) then
       qv_ref = 0.0_RP
       !
       !--- calculation of reference pot. temp.
       !--- just below the surface level.
       th_ref(ADM_kmin-1)                                 &
            = th_g                                         &
            / exp ( BV_freq**2 / CONST_GRAV                     &
                   * ( GRD_gzh(ADM_kmin) - GRD_gz(ADM_kmin-1) ) )
       !
       !--- calculation of reference pot. temp. at the first level.
       th_ref(ADM_kmin)                                   &
            = th_g                                         &
            * exp ( BV_freq**2 / CONST_GRAV                     &
                   * ( GRD_gz(ADM_kmin) - GRD_gzh(ADM_kmin) ) )
       !
       !--- calculation of reference pot. temp. at k level
       do k = ADM_kmin+1, ADM_kmax+1
          if ( GRD_gz(k) <= ZT ) then
             th_ref(k) = th_ref(k-1)                     &
                  * exp ( BV_freq**2 / CONST_GRAV              &
                         * ( GRD_gz(k) - GRD_gz(k-1) ) )
          else
             th_ref(k) = th_ref(k-1)
          endif
       enddo
       !
       !--- calculation of density at the surface
       tem_g = th_g / ( CONST_PRE00 / pre_g )**kappa
       rho_g = pre_g / CONST_Rdry / tem_g
       !
       !--- calculation of reference pressure, temperature, density
       !--- just below the ground level.
       pre_ref(ADM_kmin-1)                                &
            = pre_g                                        &
            + 0.5_RP                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_g
       tem_ref(ADM_kmin-1)                                        &
            = th_ref(ADM_kmin-1)                                  &
            / ( CONST_PRE00 / pre_ref(ADM_kmin-1) )**kappa
       rho_ref(ADM_kmin-1)        &
            = pre_ref(ADM_kmin-1) &
            / CONST_Rdry           &
            / tem_ref(ADM_kmin-1)
       !
       !--- calculation of reference pressure and density
       !--- at the first level
       pre_ref(ADM_kmin)                                         &
            = pre_ref(ADM_kmin-1)                                &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_g
       tem_ref(ADM_kmin)                                       &
            = th_ref(ADM_kmin)                                 &
            / ( CONST_PRE00 / pre_ref(ADM_kmin) )**kappa
       rho_ref(ADM_kmin)                                       &
            = pre_ref(ADM_kmin) / CONST_Rdry / tem_ref(ADM_kmin)
       !
       !--- Reference pressure and density at k level
       !---    In this caluculation, the hydrostatic balance equation
       !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer
       !---    level ( k-1/2 ). RHO is obtained by extrapolation from
       !---    the values at lower level.
       !
       !--- fist guess
       do k = ADM_kmin+1, ADM_kmax+1
          pre_ref(k) = pre_ref(k-1)              &
               - ( phi_ref(k) - phi_ref(k-1) )   &
               * ( rho_ref(k-1)                  &
               + ( rho_ref(k-1) - rho_ref(k-2) ) &
               / ( GRD_gz(k-1)  - GRD_gz(k-2) )   &
               * ( GRD_gz(k)    - GRD_gz(k-1) ) * 0.5_RP )
          tem_ref(k)                                      &
               = th_ref(k)                                &
               / ( CONST_PRE00 / pre_ref(k) )**kappa
          rho_ref(k)        &
               = pre_ref(k) &
               / CONST_Rdry  &
               / tem_ref(k)
       enddo
       !--- hydro static balance adjustment
       do k = ADM_kmin+1, ADM_kmax+1
          do
             tem_ref(k) = th_ref(k)                                &
                  / ( CONST_PRE00 / pre_ref(k) )**kappa
             rho_ref(k) = &
                  ( pre_ref(k-1) &
                  - 0.5_RP * GRD_bfac(k) * rho_ref(k-1) &
                  * ( phi_ref(k) - phi_ref(k-1) )&
                  ) / &
                  ( CONST_Rdry*tem_ref(k) &
                  + 0.5_RP * GRD_afac(k)  &
                  * ( phi_ref(k) - phi_ref(k-1) )  &
                  )
             dpre_ref_k = rho_ref(k) * CONST_Rdry * tem_ref(k) - pre_ref(k)
             pre_ref(k) = pre_ref(k)+ dpre_ref_k
             if ( abs(dpre_ref_k) < 1.E-10_RP ) exit
             if ( ADM_have_pl ) then
                write(*,*) k,abs(dpre_ref_k)
             endif
          enddo
       enddo
    else if ( ref_type == 'TH-SP' ) then
       qv_ref = 0.0_RP
       !
       !--- calculation of reference pot. temp.
       !--- just below the surface level.
       do k=ADM_kmin-1,ADM_kmax+1
          if ( GRD_gz(k)<10000.0_RP) then
             th_ref(k) = th_g
          else
             th_ref(k) = th_g+10.0_RP/1000.0_RP*(GRD_gz(k)-10000.0_RP)
          endif
       enddo
       !
       !--- calculation of density at the surface
       tem_g = th_g / ( CONST_PRE00 / pre_g )**kappa
       rho_g = pre_g / CONST_Rdry / tem_g
       !
       !--- calculation of reference pressure, temperature, density
       !--- just below the ground level.
       pre_ref(ADM_kmin-1)                                &
            = pre_g                                        &
            + 0.5_RP                                       &
            * ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) &
            * rho_g
       tem_ref(ADM_kmin-1)                                        &
            = th_ref(ADM_kmin-1)                                  &
            / ( CONST_PRE00 / pre_ref(ADM_kmin-1) )**kappa
       rho_ref(ADM_kmin-1)        &
            = pre_ref(ADM_kmin-1) &
            / CONST_Rdry           &
            / tem_ref(ADM_kmin-1)
       !
       !--- calculation of reference pressure and density
       !--- at the first level
       pre_ref(ADM_kmin)                                         &
            = pre_ref(ADM_kmin-1)                                &
            - ( phi_ref(ADM_kmin) - phi_ref(ADM_kmin-1) ) * rho_g
       tem_ref(ADM_kmin)                                       &
            = th_ref(ADM_kmin)                                 &
            / ( CONST_PRE00 / pre_ref(ADM_kmin) )**kappa
       rho_ref(ADM_kmin)                                       &
            = pre_ref(ADM_kmin) / CONST_Rdry / tem_ref(ADM_kmin)
       !
       !--- Reference pressure and density at k level
       !---    In this caluculation, the hydrostatic balance equation
       !---    ( dP/dz=-RHO dPHI/dz ) is applied at the half integer
       !---    level ( k-1/2 ). RHO is obtained by extrapolation from
       !---    the values at lower level.
       !
       !--- fist guess
       do k = ADM_kmin+1, ADM_kmax+1
          pre_ref(k) = pre_ref(k-1)              &
               - ( phi_ref(k) - phi_ref(k-1) )   &
               * ( rho_ref(k-1)                  &
               + ( rho_ref(k-1) - rho_ref(k-2) ) &
               / ( GRD_gz(k-1)  - GRD_gz(k-2) )   &
               * ( GRD_gz(k)    - GRD_gz(k-1) ) * 0.5_RP )
          tem_ref(k)                                      &
               = th_ref(k)                                &
               / ( CONST_PRE00 / pre_ref(k) )**kappa
          rho_ref(k)        &
               = pre_ref(k) &
               / CONST_Rdry  &
               / tem_ref(k)
       enddo
       !--- hydro static balance adjustment
       do k = ADM_kmin+1, ADM_kmax+1
          do
             tem_ref(k) = th_ref(k)                                &
                  / ( CONST_PRE00 / pre_ref(k) )**kappa
             rho_ref(k) = &
                  ( pre_ref(k-1) &
                  - 0.5_RP * GRD_bfac(k) * rho_ref(k-1) &
                  * ( phi_ref(k) - phi_ref(k-1) )&
                  ) / &
                  ( CONST_Rdry*tem_ref(k) &
                  + 0.5_RP * GRD_afac(k)  &
                  * ( phi_ref(k) - phi_ref(k-1) )  &
                  )
             dpre_ref_k = rho_ref(k) * CONST_Rdry * tem_ref(k) - pre_ref(k)
             pre_ref(k) = pre_ref(k)+ dpre_ref_k
             if ( abs(dpre_ref_k) < 1.E-10_RP ) exit
             if ( ADM_have_pl ) then
                write(*,*) k,abs(dpre_ref_k)
             endif
          enddo
       enddo
    else if ( ref_type == 'OOYAMA' ) then
       call ooyama_reference
       do k = 1, ADM_kall
          pre_ref(k) = rho_ref(k) * tem_ref(k) &
               * ( (1.0_RP-qv_ref(k))*CONST_Rdry+qv_ref(k)*CONST_Rvap )
          th_ref(k)                                        &
               = tem_ref(k)                                &
               * ( CONST_PRE00 / pre_ref(k) )**kappa
       enddo
    endif

!   if ( ref_type == 'GCSS_CASE1' ) then
    if ( ref_type == 'GCSS_CASE1' .or. ref_type=='TWP-ICE') then ! 11/08/13 A.Noda [add]
       call gcss_reference
       do k = 1, ADM_kall
          tem_ref(k) = th_ref(k)                                &
               / ( CONST_PRE00 / pre_ref(k) )**kappa
          rho_ref(k) = pre_ref(k)/tem_ref(k) &
               / ( (1.0_RP-qv_ref(k))*CONST_Rdry+qv_ref(k)*CONST_Rvap )
       enddo
    endif


  end subroutine set_referencestate

  !-----------------------------------------------------------------------------
  subroutine gcss_reference
    implicit none

    integer :: fid

    fid = IO_get_available_fid()
    Open(fid,file=Trim(ref_fname),status='old',form='unformatted')
    read(fid) th_ref(:)
    read(fid) pre_ref(:)
    read(fid) qv_ref(:)
    close(fid)

  end subroutine gcss_reference

  !-----------------------------------------------------------------------------
  subroutine ooyama_reference
    use mod_adm, only: &
       ADM_kmin, &
       ADM_kmax
    use mod_grd, only: &
       GRD_gz
    implicit none

    character(len=H_LONG) :: fname = 'ooyama_profile.dat'

    real(RP), allocatable :: z_s  (:)
    real(RP), allocatable :: rho_s(:)
    real(RP), allocatable :: tem_s(:)
    real(RP), allocatable :: qv_s (:)

    integer :: fid
    integer :: kmax
    integer :: k, kk, kp

    real(RP) :: lag_intpl
    real(RP) :: z,z1,p1,z2,p2,z3,p3
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ((z-z2)*(z-z3))/((z1-z2)*(z1-z3))*p1 &
                                   + ((z-z1)*(z-z3))/((z2-z1)*(z2-z3))*p2 &
                                   + ((z-z1)*(z-z2))/((z3-z1)*(z3-z2))*p3
    !---------------------------------------------------------------------------

    !--- read sounding data ( ooyama(2001) )
    fid = IO_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          status = 'old',       &
          form   = 'unformatted')

       read(fid) kmax

       allocate( z_s  (kmax) )
       allocate( rho_s(kmax) )
       allocate( tem_s(kmax) )
       allocate( qv_s (kmax) )

       read(fid) z_s  (:)
       read(fid) rho_s(:)
       read(fid) tem_s(:)
       read(fid) qv_s (:)
    close(fid)

    !--- initialization
    do k = ADM_kmin, ADM_kmax
    do kk = 1, kmax-1
       if (       z_s(kk)   <= GRD_gz(k) &
            .AND. z_s(kk+1) >= GRD_gz(k) ) then
          if ( kk == 1 ) then
             kp = 2
          else
             kp = kk
          endif

          rho_ref(k) = lag_intpl(GRD_gz(k),z_s(kp+1),rho_s(kp+1),z_s(kp),rho_s(kp),z_s(kp-1),rho_s(kp-1))
          tem_ref(k) = lag_intpl(GRD_gz(k),z_s(kp+1),tem_s(kp+1),z_s(kp),tem_s(kp),z_s(kp-1),tem_s(kp-1))
          qv_ref (k) = lag_intpl(GRD_gz(k),z_s(kp+1),qv_s (kp+1),z_s(kp),qv_s (kp),z_s(kp-1),qv_s (kp-1))
          exit
       endif
    enddo
    enddo

    k  = ADM_kmin-1
    kp = 2

    rho_ref(k) = lag_intpl(GRD_gz(k),z_s(kp+1),rho_s(kp+1),z_s(kp),rho_s(kp),z_s(kp-1),rho_s(kp-1))
    tem_ref(k) = lag_intpl(GRD_gz(k),z_s(kp+1),tem_s(kp+1),z_s(kp),tem_s(kp),z_s(kp-1),tem_s(kp-1))
    qv_ref (k) = lag_intpl(GRD_gz(k),z_s(kp+1),qv_s (kp+1),z_s(kp),qv_s (kp),z_s(kp-1),qv_s (kp-1))

    k  = ADM_kmax+1
    kp = kmax-1

    rho_ref(k) = lag_intpl(GRD_gz(k),z_s(kp+1),rho_s(kp+1),z_s(kp),rho_s(kp),z_s(kp-1),rho_s(kp-1))
    tem_ref(k) = lag_intpl(GRD_gz(k),z_s(kp+1),tem_s(kp+1),z_s(kp),tem_s(kp),z_s(kp-1),tem_s(kp-1))
    qv_ref (k) = lag_intpl(GRD_gz(k),z_s(kp+1),qv_s (kp+1),z_s(kp),qv_s (kp),z_s(kp-1),qv_s (kp-1))

    deallocate( rho_s )
    deallocate( tem_s )
    deallocate( qv_s  )

    return
  end subroutine ooyama_reference

  !-----------------------------------------------------------------------------
  !> generation of basic state from reference state
  subroutine set_basicstate
    use scale_const, only: &
       CONST_GRAV
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_grd, only: &
       GRD_Z,     &
       GRD_vz,    &
       GRD_vz_pl
    use mod_vintrpl, only: &
       VINTRPL_zstar_level
    use mod_bndcnd, only: &
       bndcnd_thermo
    use mod_thrmdyn, only: &
       THRMDYN_qd,  &
       THRMDYN_rho, &
       THRMDYN_th
    use mod_runconf, only: &
       TRC_VMAX, &
       I_QV
    implicit none

    real(RP) :: q_bs    (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP) :: q_bs_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP) :: qd_bs   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qd_bs_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: k, l
    !---------------------------------------------------------------------------

    !--- calculation of geo-potential
    phi(:,:,:) = CONST_GRAV * GRD_vz(:,:,:,GRD_Z)
    if ( ADM_have_pl ) then
       phi_pl(:,:,:) = CONST_GRAV * GRD_vz_pl(:,:,:,GRD_Z)
    endif

    if( ref_type == 'NOBASE' ) return

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       pre_bs(:,k,l) = pre_ref(k)
       tem_bs(:,k,l) = tem_ref(k)
       qv_bs (:,k,l) = qv_ref (k)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          pre_bs_pl(:,k,l) = pre_ref(k)
          tem_bs_pl(:,k,l) = tem_ref(k)
          qv_bs_pl (:,k,l) = qv_ref (k)
       enddo
       enddo
    endif

    !-- from z-level to zstar-level
    call VINTRPL_zstar_level( pre_bs, pre_bs_pl, .false. )
    call VINTRPL_zstar_level( tem_bs, tem_bs_pl, .false. )
    call VINTRPL_zstar_level( qv_bs,  qv_bs_pl,  .false. )

    !--- Setting of mass concentration [Note] The basic state is "dry" and TKE=0
    q_bs(:,:,:,:)    = 0.0_RP
    q_bs(:,:,:,I_QV) = qv_bs(:,:,:)

    call THRMDYN_qd( ADM_gall,       & ! [IN]
                     ADM_kall,       & ! [IN]
                     ADM_lall,       & ! [IN]
                     q_bs (:,:,:,:), & ! [IN]
                     qd_bs(:,:,:)    ) ! [OUT]

    call THRMDYN_rho( ADM_gall,        & ! [IN]
                      ADM_kall,        & ! [IN]
                      ADM_lall,        & ! [IN]
                      tem_bs(:,:,:),   & ! [IN]
                      pre_bs(:,:,:),   & ! [IN]
                      qd_bs (:,:,:),   & ! [IN]
                      q_bs  (:,:,:,:), & ! [IN]
                      rho_bs(:,:,:)    ) ! [OUT]

    !--- set boundary conditions of basic state
    do l = 1, ADM_lall
       call bndcnd_thermo( ADM_gall,      & ! [IN]
                           tem_bs(:,:,l), & ! [INOUT]
                           rho_bs(:,:,l), & ! [INOUT]
                           pre_bs(:,:,l), & ! [INOUT]
                           phi   (:,:,l)  ) ! [IN]
    enddo

    call THRMDYN_th( ADM_gall,      & ! [IN]
                     ADM_kall,      & ! [IN]
                     ADM_lall,      & ! [IN]
                     tem_bs(:,:,:), & ! [IN]
                     pre_bs(:,:,:), & ! [IN]
                     th_bs (:,:,:)  ) ! [OUT]

    if ( ADM_have_pl ) then
       q_bs_pl(:,:,:,:)    = 0.0_RP
       q_bs_pl(:,:,:,I_QV) = qv_bs_pl(:,:,:)

       call THRMDYN_qd( ADM_gall_pl,       & ! [IN]
                        ADM_kall,          & ! [IN]
                        ADM_lall_pl,       & ! [IN]
                        q_bs_pl (:,:,:,:), & ! [IN]
                        qd_bs_pl(:,:,:)    ) ! [OUT]

       call THRMDYN_rho( ADM_gall_pl,        & ! [IN]
                         ADM_kall,           & ! [IN]
                         ADM_lall_pl,        & ! [IN]
                         tem_bs_pl(:,:,:),   & ! [IN]
                         pre_bs_pl(:,:,:),   & ! [IN]
                         qd_bs_pl (:,:,:),   & ! [IN]
                         q_bs_pl  (:,:,:,:), & ! [IN]
                         rho_bs_pl(:,:,:)    ) ! [OUT]


       do l = 1, ADM_lall_pl
          call bndcnd_thermo( ADM_gall_pl,      & ! [IN]
                              tem_bs_pl(:,:,l), & ! [INOUT]
                              rho_bs_pl(:,:,l), & ! [INOUT]
                              pre_bs_pl(:,:,l), & ! [INOUT]
                              phi_pl   (:,:,l)  ) ! [IN]
       enddo

       call THRMDYN_th( ADM_gall_pl,      & ! [IN]
                        ADM_kall,         & ! [IN]
                        ADM_lall_pl,      & ! [IN]
                        tem_bs_pl(:,:,:), & ! [IN]
                        pre_bs_pl(:,:,:), & ! [IN]
                        th_bs_pl (:,:,:)  ) ! [OUT]
    endif

    return
  end subroutine set_basicstate

  !-----------------------------------------------------------------------------
  subroutine bsstate_output_ref( basename )
    use mod_adm, only: &
       ADM_prc_me,     &
       ADM_prc_run_master
    implicit none

    Character(*), Intent(in) :: basename
    integer :: fid

    !--- output
    if ( ADM_prc_me==ADM_prc_run_master) then
       fid = IO_get_available_fid()
       Open(fid,file=Trim(basename),form='unformatted')
       Write(fid) pre_ref(:)
       Write(fid) tem_ref(:)
       Write(fid) qv_ref(:)
       Close(fid)
       !
    endif

    return
  end subroutine bsstate_output_ref

  !-----------------------------------------------------------------------------
  subroutine bsstate_input_ref( basename )
    use scale_const, only: &
         CONST_GRAV,  &
         CONST_Rdry,  &
         CONST_Rvap,  &
         CONST_CPdry, &
         CONST_PRE00
    use mod_adm, only :  &
         ADM_kall
    use mod_grd, only :  &
         GRD_gz
    implicit none

    Character(*), Intent(in) :: basename
    integer :: fid

    real(RP) :: kappa
    integer :: k

    kappa = CONST_Rdry / CONST_CPdry

    !--- input
    fid = IO_get_available_fid()
    Open(fid,file=Trim(basename),status='old',form='unformatted')
    read(fid) pre_ref(:)
    read(fid) tem_ref(:)
    read(fid) qv_ref(:)
    Close(fid)
    !
    !--- additional reference state.
    do k = 1, ADM_kall
       th_ref(k)                                        &
            = tem_ref(k)                                &
            * ( CONST_PRE00 / pre_ref(k) )**kappa
       phi_ref(k) = CONST_GRAV * GRD_gz(k)
       rho_ref(k) = pre_ref(k)/ ( (1.0_RP-qv_ref(k))*CONST_Rdry+qv_ref(k)*CONST_Rvap )/ tem_ref(k)
    enddo

  end subroutine bsstate_input_ref

end module mod_bsstate
