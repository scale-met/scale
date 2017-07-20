!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          kajino13 code
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-03-25 (Y.Sato)  [new]
!!
!<
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
module scale_atmos_phy_ae_kajino13
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_const, only: &
      mwair => CONST_Mdry, &              ! molecular weight for dry air
      mwwat => CONST_Mvap, &              ! mean molecular weight for water vapor [g/mol]
      dnwat => CONST_DWATR, &             ! water density   [kg/m3]
      rair  => CONST_Rdry, &              ! gas constant of dry air            [J/kg/K]
      rwat  => CONST_Rvap, &              ! gas constant of water              [J/kg/K]
      rgas  => CONST_R, &                 ! universal gas constant             [J/mol/K]
      stdatmpa =>  CONST_Pstd, &          ! standard pressure                   [Pa]
      stdtemp  =>  CONST_TEM00, &         ! standard temperature                [K]
      pi    => CONST_PI                   ! pi
  use scale_atmos_aerosol, only: &
      N_AE
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_kajino13_config
  public :: ATMOS_PHY_AE_kajino13_setup
  public :: ATMOS_PHY_AE_kajino13
  public :: ATMOS_PHY_AE_kajino13_EffectiveRadius
  public :: ATMOS_PHY_AE_kajino13_mkinit

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, private :: QA_AE  = 0

  character(len=H_SHORT), public, target, allocatable :: ATMOS_PHY_AE_kajino13_NAME(:)
  character(len=H_MID)  , public, target, allocatable :: ATMOS_PHY_AE_kajino13_DESC(:)
  character(len=H_SHORT), public, target, allocatable :: ATMOS_PHY_AE_kajino13_UNIT(:)

  real(RP), public, target :: ATMOS_PHY_AE_kajino13_DENS(N_AE) ! hydrometeor density [kg/m3]=[g/L]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: AE_CTG = 1
  integer, parameter :: GAS_CTG = 2
  integer, parameter :: N_ATR = 5
  integer, allocatable :: NKAP(:)
  integer, allocatable :: NSIZ(:)

  integer :: QAES
  integer :: QAEE

  integer, parameter :: IC_MIX   =  1
  integer, parameter :: IC_SEA   =  2
  integer, parameter :: IC_DUS   =  3

  integer, parameter :: IG_H2SO4 =  1
  integer, parameter :: IG_CGAS  =  2


  ! physical parameters
!  real(RP), parameter  :: mwair = 28.9628_RP              ! molecular weight for dry air
!  real(RP), parameter  :: mwwat = 18.0153_RP              ! mean molecular weight for water vapor [g/mol]
!  real(RP), parameter  :: dnwat = 1.E3_RP                 ! water density   [kg/m3]
!  real(RP), parameter  :: rair    = 287._RP               ! gas constant of dry air            [J/kg/K]
!  real(RP), parameter  :: rwat    = rair*mwair/mwwat      ! gas constant of water              [J/kg/K]
!  real(RP), parameter  :: rgas    = 8.31447_RP            ! universal gas constant             [J/mol/K]
!  real(RP), parameter  :: stdatmpa =  101325._RP          ! standard pressure                   [Pa]
!  real(RP), parameter  :: stdtemp  =  273.15_RP           ! standard temperature                [K]

  real(RP), parameter  :: avo     = 6.0221367e23_RP       ! avogadro number [#/mol]
  real(RP), parameter  :: boltz   = 1.38048e-23_RP        ! boltzmann constant [J K-1]
  real(RP), parameter  :: conv_n2m = 1.E6_RP/avo*1.E6_RP &!
                                   *98._RP                ! ug/m3 = conv_n2m * #/cm3
  real(RP), parameter  :: conv_m2n = 1._RP/conv_n2m       ! #/cm3 = conv_m2n * ug/m3
  real(RP), parameter  :: rhod_ae  = 1.83_RP              ! particle density [g/cm3] sulfate assumed
  real(RP), parameter  :: rho_kg   = rhod_ae*1.E3_RP      ! particle density [kg/m3] sulfate assumed
  real(RP), parameter  :: grav     = 9.80622_RP           ! mean gravitational acceleration [m/s2]
  real(RP), parameter  :: conv_ms_vl = 1.E-12_RP/rhod_ae  ! mass[ug/m3] to M3(volume)[m3/m3] rhod
  real(RP), parameter  :: conv_vl_ms = rhod_ae/1.E-12_RP  ! M3(volume)[m3/m3] to mass[ug/m3]
  real(RP), parameter  :: mwrat_s6   = 96._RP/98._RP      ! molecular weight ratio (part/gas) of sulfate
  real(RP), parameter  :: mwrat_s6_i = 1._RP/mwrat_s6     ! molecular weight ratio (gas/part) of sulfate
  real(RP), parameter  :: diffsulf =  9.36e-6_RP          ! std. molecular diffusivity of sulfuric acid [m2/s]
  real(RP), parameter  :: mwh2so4  =  98._RP              ! molecular weight for h2so4 gas      [g/mol]
  real(RP)             :: logk_aenucl = -12.4_RP          !constant coefficient for kinetic nucleation [-]
  real(RP)             :: pi6                             != pi / 6._RP              ! pi/6
  real(RP)             :: sixpi                           != 6._RP / pi              ! 6/pi
  real(RP)             :: forpi                           != 4._RP / pi              ! 4/pi
  ! parameters for condensation/coagulation
  real(RP), parameter  :: c1=1.425e-6_RP, c2=0.5039_RP, c3=108.3_RP   ! Perry's [Tableo 2-312]
  real(RP)             :: c4                                          ! Perry's [Tableo 2-312]
  integer              :: nbins_out = 1024

  !(attributes: zero chemistry version)
!  integer, parameter :: n_atr  = 5             !number of attributes
  integer, parameter :: ia_m0  = 1             !1. number conc        [#/m3]
  integer, parameter :: ia_m2  = 2             !2. 2nd mom conc       [m2/m3]
  integer, parameter :: ia_m3  = 3             !3. 3rd mom conc       [m3/m3]
  integer, parameter :: ia_ms  = 4             !4. mass conc          [ug/m3]
  integer, parameter :: ia_kp  = 5             !5. mean kappa * volume [-] ( ia_kp*ia_m3 nisuru)
  integer, parameter :: ik_out = 1
  !(set in subroutine aerosol_settings)
!  integer :: n_ctg   !number of category
  integer, allocatable :: n_siz(:)              !number of size bins           (n_ctg)
  real(RP),allocatable :: d_min(:)              !lower bound of 1st size bin   (n_ctg)
  real(RP),allocatable :: d_max(:)              !upper bound of last size bin  (n_ctg)
  integer, allocatable :: n_kap(:)              !number of kappa bins          (n_ctg)
  real(RP),allocatable :: k_min(:)              !lower bound of 1st kappa bin  (n_ctg)
  real(RP),allocatable :: k_max(:)              !upper bound of last kappa bin (n_ctg)
  !(diagnosed in subroutine aerosol_settings)

  !--- bin settings (lower bound, center, upper bound)
  real(RP),allocatable :: d_lw(:,:), d_ct(:,:), d_up(:,:)  !diameter [m]
  real(RP),allocatable :: k_lw(:,:), k_ct(:,:), k_up(:,:)  !kappa    [-]
  real(RP) :: dlogd, dk                                    !delta log(D), delta K
!  real(RP) :: deltt                                        ! dt

  !--- coagulation rule (i+j=k)
  integer, allocatable :: is_i(:), is_j(:), is_k(:)    !(mcomb)
  integer, allocatable :: ik_i(:), ik_j(:), ik_k(:)    !(mcomb)
  integer, allocatable :: ic_i(:), ic_j(:), ic_k(:)    !(mcomb)
  integer :: mcomb !combinations of sections
  integer :: is1, is2, mc, ik, is0, ic, ia0

  logical :: flag_npf = .false.
  logical :: flag_cond = .true.
  logical :: flag_coag = .true.
  logical :: flag_ccn_interactive = .true.
  logical :: flag_regeneration    = .true.
!  real(RP) :: m0_init = 0.E0      ! initial total num. conc. of modes (Atk,Acm,Cor) [#/m3]
!  real(RP) :: dg_init = 80.E-9    ! initial number equivalen diameters of modes     [m]
!  real(RP) :: sg_init = 1.6       ! initial standard deviation                      [-]
  real(RP) :: h2so4dt = 5.E-6_RP   ! h2so4 production rate (Temporal)                [ug/m3/s]
  real(RP) :: ocgasdt = 8.E-5_RP   ! other condensational bas production rate 16*h2so4dt (see Kajino et al. 2013)
!  real(RP) :: c_ratio = 16.0_RP    ! ratio of condensable mass to h2so4 (after NPF)  [-]
  real(RP) :: c_kappa = 0.3_RP     ! hygroscopicity of condensable mass              [-]
  real(RP) :: t_npf = 21600.0_RP   ! time until the gas to particle conversion occure[sec]
  integer :: n_ctg                            ! number of category
  integer :: n_trans                          ! number of total transport variables
  integer :: n_siz_max                        ! maximum number of size bins
  integer :: n_kap_max                        ! maximum number of kappa bins
  integer, allocatable :: it_procs2trans(:,:,:,:) !procs to trans conversion
  integer, allocatable :: ia_trans2procs(:) !trans to procs conversion
  integer, allocatable :: is_trans2procs(:) !trans to procs conversion
  integer, allocatable :: ik_trans2procs(:) !trans to procs conversion
  integer, allocatable :: ic_trans2procs(:) !trans to procs conversion
  real(RP), allocatable :: rnum_out(:)
  real(RP) :: dg_reg = 5.E-7_RP             ! dg of regenerated aerosol (droplet mode 500 nm) [m]
  real(RP) :: sg_reg = 1.6_RP               ! sg of regenerated aerosol [-]
  integer  :: is0_reg                       ! size bin of regenerated aerosol
  character(len=H_SHORT),allocatable :: ctg_name(:)
  real(RP), parameter :: cleannumber = 1.E-3 ! tiny number of aerosol per bin for neglectable mass,
                                             !  D =10 um resulted in 1.e-6 ug/m3
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_AE_kajino13_config( &
       AE_TYPE, &
       QA, QS )
    use scale_process, only: &
       PRC_MPIstop
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    character(len=*), intent(in)  :: AE_TYPE
    integer,          intent(out) :: QA
    integer,          intent(out) :: QS

    integer, allocatable :: aero_idx(:,:,:,:)
    integer :: n_kap_max, n_siz_max, ncat_max
    integer :: NASIZ(3), NAKAP(3)
    character(len=H_SHORT) :: attribute, catego, aunit

    NAMELIST / PARAM_TRACER_KAJINO13 / &
       AE_CTG, &
       NASIZ,  &
       NAKAP

    integer :: m, ierr, ik, ic, ia0, is0

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Aerosol Tracer] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Tracers for Kajino(2013) scheme'

    if ( AE_TYPE /= 'KAJINO13' .AND. AE_TYPE /= 'NONE' ) then
       write(*,*) 'xxx ATMOS_PHY_AE_TYPE is not KAJINO13. Check!'
       call PRC_MPIstop
    endif

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

!    do ia0 = 1, N_ATR
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
!    do ia0 = 1, N_ATR
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

!    do ia0 = 1, N_ATR
!      if( ia0 == 1 ) then
!         write(attribute,'(a)') "Number"
!         write(aunit,'(a)') "num/kg"
!      elseif( ia0 == 2 ) then
!         write(attribute,'(a)') "Section"
!         write(aunit,'(a)') "m2/kg"
!      elseif( ia0 == 3 ) then
!         write(attribute,'(a)') "Volume"
!         write(aunit,'(a)') "m3/kg"
!      elseif( ia0 == 4 ) then
!         write(attribute,'(a)') "Mass"
!         write(aunit,'(a)') "kg/kg"
!      elseif( ia0 == 5 ) then
!         write(attribute,'(a)') "kpXmass"
!         write(aunit,'(a)') "kg/kg"
!      endif
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

    return
  end subroutine ATMOS_PHY_AE_kajino13_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_kajino13_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_AE
    implicit none

    real(RP), allocatable :: d_min_inp(:)
    real(RP), allocatable :: d_max_inp(:)
    real(RP), allocatable :: k_min_inp(:)
    real(RP), allocatable :: k_max_inp(:)
    integer , allocatable :: n_kap_inp(:)

    real(RP), parameter :: d_min_def = 1.E-9_RP ! default lower bound of 1st size bin
    real(RP), parameter :: d_max_def = 1.E-5_RP ! upper bound of last size bin
    integer , parameter :: n_kap_def = 1        ! number of kappa bins
    real(RP), parameter :: k_min_def = 0.E0_RP  ! lower bound of 1st kappa bin
    real(RP), parameter :: k_max_def = 1.E0_RP  ! upper bound of last kappa bin

    NAMELIST / PARAM_ATMOS_PHY_AE_KAJINO13 / &
       h2so4dt, &
       ocgasdt, &
!       c_ratio, &
       c_kappa, &
       t_npf, &
!       d_min_inp, &
!       d_max_inp, &
!       k_min_inp, &
!       k_max_inp, &
!       n_kap_inp, &
       flag_npf, &
       flag_cond, &
       flag_coag, &
       flag_ccn_interactive, &
       flag_regeneration,    &
       dg_reg, &
       sg_reg, &
       logk_aenucl, &
       nbins_out

    integer :: it, ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Aerosol] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Kajino(2013) scheme'

    !--- setup parameter
    pi6   = pi / 6._RP              ! pi/6
    sixpi = 6._RP / pi              ! 6/pi
    forpi = 4._RP / pi              ! 4/pi


!    deltt = TIME_DTSEC_ATMOS_PHY_AE
    n_ctg = AE_CTG
    allocate( rnum_out(nbins_out) )
    allocate( n_siz(n_ctg) )
    allocate( d_min(n_ctg) )
    allocate( d_max(n_ctg) )
    allocate( n_kap(n_ctg) )
    allocate( k_min(n_ctg) )
    allocate( k_max(n_ctg) )
    allocate( d_min_inp(n_ctg) )
    allocate( d_max_inp(n_ctg) )
    allocate( n_kap_inp(n_ctg) )
    allocate( k_min_inp(n_ctg) )
    allocate( k_max_inp(n_ctg) )
    allocate( ctg_name(n_ctg) )

    n_siz(1:n_ctg) = NSIZ(1:n_ctg)       ! number of size bins
    d_min(1:n_ctg) = d_min_def ! lower bound of 1st size bin
    d_max(1:n_ctg) = d_max_def ! upper bound of last size bin
    n_kap(1:n_ctg) = n_kap_def ! number of kappa bins
    k_min(1:n_ctg) = k_min_def ! lower bound of 1st kappa bin
    k_max(1:n_ctg) = k_max_def ! upper bound of last kappa bin

    do it = 1, n_ctg
     if( n_ctg == 1 ) then
       write(ctg_name(it),'(a)') "Sulfate"
     elseif( n_ctg == 2 ) then
       write(ctg_name(it),'(a)') "Seasalt"
     elseif( n_ctg == 3 ) then
       write(ctg_name(it),'(a)') "Dust"
     endif
    enddo

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_AE_KAJINO13,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_AE_KAJINO13. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_AE_KAJINO13)

    !--- now only the default setting is supported
!    n_siz(1:n_ctg) = NSIZ(1:n_ctg)       ! number of size bins
!    d_min(1:n_ctg) = d_min_inp(1:n_ctg)  ! lower bound of 1st size bin
!    d_max(1:n_ctg) = d_max_inp(1:n_ctg)  ! upper bound of last size bin
!    n_kap(1:n_ctg) = n_kap_inp(1:n_ctg)  ! number of kappa bins
!    k_min(1:n_ctg) = k_min_inp(1:n_ctg)  ! lower bound of 1st kappa bin
!    k_max(1:n_ctg) = k_max_inp(1:n_ctg)  ! upper bound of last kappa bin

    !--- diagnose parameters (n_trans, n_siz_max, n_kap_max)
    n_trans   = 0
    n_siz_max = 0
    n_kap_max = 0
    do ic = 1, n_ctg
      n_trans   = n_trans + n_siz(ic) * n_kap(ic) * N_ATR
      n_siz_max = max(n_siz_max, n_siz(ic))
      n_kap_max = max(n_kap_max, n_kap(ic))
    enddo

    !--- bin settings
    allocate(d_lw(n_siz_max,n_ctg))
    allocate(d_ct(n_siz_max,n_ctg))
    allocate(d_up(n_siz_max,n_ctg))
    allocate(k_lw(n_kap_max,n_ctg))
    allocate(k_ct(n_kap_max,n_ctg))
    allocate(k_up(n_kap_max,n_ctg))
    d_lw(:,:) = 0.0_RP
    d_ct(:,:) = 0.0_RP
    d_up(:,:) = 0.0_RP
    k_lw(:,:) = 0.0_RP
    k_ct(:,:) = 0.0_RP
    k_up(:,:) = 0.0_RP

    do ic = 1, n_ctg
      dlogd = (log(d_max(ic)) - log(d_min(ic)))/float(n_siz(ic))
      do is0 = 1, n_siz(ic)  !size bin
        d_lw(is0,ic) = exp(log(d_min(ic))+dlogd* float(is0-1)      )
        d_ct(is0,ic) = exp(log(d_min(ic))+dlogd*(float(is0)-0.5_RP))
        d_up(is0,ic) = exp(log(d_min(ic))+dlogd* float(is0)        )
      enddo !is (1:n_siz(ic))

      dk  = (k_max(ic) - k_min(ic))/float(n_kap(ic))
      do ik = 1, n_kap(ic)  !size bin
        k_lw(ik,ic) = k_min(ic) + dk  * float(ik-1)
        k_ct(ik,ic) = k_min(ic) + dk  *(float(ik)-0.5_RP)
        k_up(ik,ic) = k_min(ic) + dk  * float(ik)
      enddo !ik (1:n_kap(ic))

    enddo !ic (1:n_ctg)

!find size bin of regenerated aerosols
    do is0 = 1, n_siz(ic_mix)
      if (dg_reg >= d_lw(is0,ic_mix) .AND. &
          dg_reg < d_up(is0,ic_mix) ) then
       is0_reg = is0
      endif !d_lw < dg_reg < d_up
    enddo

    !--- coagulation rule
    !   [ NOTE: current version has one category and single
    !     hygroscopicity bins and thus does not consider inter-category
    !     nor inter-hygroscopicity-section coagulation ]
    mcomb = 0
    do ic = 1, n_ctg     !=1
    do ik = 1, n_kap(ic) !=1
      do is2 = 1, n_siz(ic)
      do is1 = 1, n_siz(ic)
        if ( d_ct(is2,ic) >= d_ct(is1,ic) ) then
          mcomb = mcomb + 1
        endif !d_ct(is2) >= d_ct(is1)
      enddo !is1 (1:n_siz(ic)  )
      enddo !is2 (1:n_siz(ic)  )
    enddo !ik(1:n_kap(ic))
    enddo !ic(1:n_ctg)

    allocate(is_i(mcomb))
    allocate(is_j(mcomb))
    allocate(is_k(mcomb))
    allocate(ik_i(mcomb))
    allocate(ik_j(mcomb))
    allocate(ik_k(mcomb))
    allocate(ic_i(mcomb))
    allocate(ic_j(mcomb))
    allocate(ic_k(mcomb))

    mc = 0
    do ic = 1, n_ctg     !=1
    do ik = 1, n_kap(ic) !=1
      do is2 = 1, n_siz(ic)
      do is1 = 1, n_siz(ic)
        if ( d_ct(is2,ic) >= d_ct(is1,ic) ) then
          mc = mc + 1
          is_i(mc) = is1
          ik_i(mc) = ik
          ic_i(mc) = ic
          is_j(mc) = is2
          ik_j(mc) = ik
          ic_j(mc) = ic
          is_k(mc) = is2
          ik_k(mc) = ik
          ic_k(mc) = ic
        endif !d_ct(is2) >= d_ct(is1)
      enddo !is1 (1:n_siz(ic)  )
      enddo !is2 (1:n_siz(ic)  )
    enddo !ik(1:n_kap(ic))
    enddo !ic(1:n_ctg)

    !--- gas concentration
!    conc_h2so4 = 0.0_RP
!    conc_cgas = 0.0_RP

    allocate( it_procs2trans(N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( ia_trans2procs(n_trans) )
    allocate( is_trans2procs(n_trans) )
    allocate( ik_trans2procs(n_trans) )
    allocate( ic_trans2procs(n_trans) )

    it_procs2trans(:,:,:,:)= -999
    ia_trans2procs(:)      = 0
    is_trans2procs(:)      = 0
    ik_trans2procs(:)      = 0
    ic_trans2procs(:)      = 0

    !--- get pointer for trans2procs, procs2trans
    it = 0
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)  !size bin
    do ia0 = 1, N_ATR       !attributes
      it = it + 1
      it_procs2trans(ia0,is0,ik,ic)= it
      ia_trans2procs(it)         = ia0
      is_trans2procs(it)         = is0
      ik_trans2procs(it)         = ik
      ic_trans2procs(it)         = ic
    enddo !ia (1:n_atr_prog )
    enddo !is (1:n_siz(ic)  )
    enddo !ik (1:n_kap(ic)  )
    enddo !ic (1:n_ctg      )

    ATMOS_PHY_AE_kajino13_DENS(:) = rhod_ae

    return
  end subroutine ATMOS_PHY_AE_kajino13_setup

  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics
  subroutine ATMOS_PHY_AE_kajino13( &
       QQA,  &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       EMIT, &
       NREG, &
       QTRC, &
       CN,   &
       CCN,  &
       RHOQ_t_AE )
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       CONST_CPdry, &
       CONST_CVdry, &
       CONST_Rvap, &
       CONST_PRE00, &
       CONST_Rdry
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_atmos_saturation, only: &
       pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_AE
    use scale_history, only: &
       HIST_in
    implicit none
    integer,  intent(in)    :: QQA
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: EMIT(KA,IA,JA,QQA)
    real(RP), intent(in)    :: NREG(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(out)   :: CN(KA,IA,JA)
    real(RP), intent(out)   :: CCN(KA,IA,JA)
    real(RP), intent(inout) :: RHOQ_t_AE(KA,IA,JA,QA)

    real(RP) :: QTRC0(KA,IA,JA,QA)
    real(RP) :: QTRC1(KA,IA,JA,QA)

    !--- local
    real(RP) :: pres_ae(KA,IA,JA)
    real(RP) :: temp_ae(KA,IA,JA)
    real(RP) :: qv_ae(KA,IA,JA)
    real(RP) :: ssliq_ae(KA,IA,JA)
    real(RP) :: th(KA,IA,JA)
    real(RP) :: q(KA,IA,JA,QA)
    real(RP) :: qdry(KA,IA,JA)
    real(RP) :: rrhog(KA,IA,JA)
    real(RP) :: cva(KA,IA,JA)
    real(RP) :: cpa(KA,IA,JA)
    real(RP) :: t_ccn, t_cn
    real(RP) :: Rmoist
    real(RP) :: qsat_tmp
    real(RP) :: m0_reg, m2_reg, m3_reg      !regenerated aerosols [m^k/m3]
    real(RP) :: ms_reg                      !regenerated aerosol mass [ug/m3]
    real(RP) :: reg_factor_m2,reg_factor_m3 !to save cpu time for moment conversion
    !--- aerosol variables
    real(RP),allocatable :: aerosol_procs(:,:,:,:) !(n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP),allocatable :: aerosol_activ(:,:,:,:) !(n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP),allocatable :: emis_procs(:,:,:,:)    !(n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP),allocatable :: emis_gas(:)            !emission of gas
    real(RP) :: total_aerosol_mass(KA,IA,JA,n_ctg)
    real(RP) :: total_aerosol_number(KA,IA,JA,n_ctg)
    real(RP) :: total_emit_aerosol_mass(KA,IA,JA,n_ctg)
    real(RP) :: total_emit_aerosol_number(KA,IA,JA,n_ctg)
    !--- gas
!    real(RP),allocatable :: conc_h2so4(:,:,:,:)    !concentration [ug/m3]
!    real(RP) :: conc_gas(KA,IA,JA,GAS_CTG)    !concentration [ug/m3]
    real(RP) :: conc_gas(GAS_CTG)    !concentration [ug/m3]
    character(len=H_LONG) :: ofilename
    integer :: i, j, k, iq, it

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Aerosol(kajino13)'

    !--- Negative fixer
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
      do ic = 1, n_ctg       !aerosol category
      do ik = 1, n_kap(ic)   !kappa bin
      do is0 = 1, n_siz(ic)   !size bin
        if (QTRC(k,i,j,QAES-1+it_procs2trans(ia_m0,is0,ik,ic)) < 0.0_RP .or. &
            QTRC(k,i,j,QAES-1+it_procs2trans(ia_m2,is0,ik,ic)) < 0.0_RP .or. &
            QTRC(k,i,j,QAES-1+it_procs2trans(ia_m3,is0,ik,ic)) < 0.0_RP .or. &
            QTRC(k,i,j,QAES-1+it_procs2trans(ia_ms,is0,ik,ic)) < 0.0_RP .or. &
            QTRC(k,i,j,QAES-1+it_procs2trans(ia_kp,is0,ik,ic)) < 0.0_RP ) then
          QTRC(k,i,j,QAES-1+it_procs2trans(1:N_ATR,is0,ik,ic)) = 0.0_RP
        endif
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       QTRC0(k,i,j,iq) = QTRC(k,i,j,iq) ! save
    enddo
    enddo
    enddo
    enddo

    allocate( aerosol_procs (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( aerosol_activ (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( emis_procs    (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( emis_gas      (GAS_CTG)  )

    reg_factor_m2 = dg_reg**2._RP * exp( 2.0_RP *(log(sg_reg)**2._RP)) !m0_reg to m2_reg
    reg_factor_m3 = dg_reg**3._RP * exp( 4.5_RP *(log(sg_reg)**2._RP)) !m0_reg to m3_reg

    !--- convert SCALE variable to zerochem variable

    aerosol_procs(:,:,:,:) = 0.0_RP
    aerosol_activ(:,:,:,:) = 0.0_RP
    emis_procs(:,:,:,:) = 0.0_RP
    emis_gas(:) = 0.0_RP
    pres_ae(:,:,:) = 0.0_RP
    temp_ae(:,:,:) = 0.0_RP
    qv_ae(:,:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          rrhog(k,i,j) = 1.0_RP / DENS(k,i,j)
          th(k,i,j) = RHOT(k,i,j) * rrhog(k,i,j)
       enddo
       do k = KS, KE
          CALC_QDRY( qdry(k,i,j), QTRC, TRACER_MASS, k, i, j, iq )
       enddo
       do k = KS, KE
          CALC_CV( cva(k,i,j), qdry(k,i,j), QTRC, k, i, j, iq, CONST_CVdry, TRACER_CV )
       enddo
       do k = KS, KE
          CALC_R( Rmoist, qdry(k,i,j), QTRC, k, i, j, iq, CONST_Rdry, TRACER_R )
          cpa(k,i,j) = cva(k,i,j) + Rmoist
          CALC_PRE( pres_ae(k,i,j), DENS(k,i,j), th(k,i,j), Rmoist, cpa(k,i,j), CONST_PRE00 )
          temp_ae(k,i,j) = pres_ae(k,i,j) / ( DENS(k,i,j) * Rmoist )
          if ( I_QV > 0 ) then
             qv_ae(k,i,j) = QTRC(k,i,j,I_QV)
             !--- calculate super saturation of water
             call pres2qsat_liq( qsat_tmp,temp_ae(k,i,j),pres_ae(k,i,j) )
             ssliq_ae(k,i,j) = qv_ae(k,i,j)/qsat_tmp - 1.0_RP
          else
             ssliq_ae(k,i,j) = - 1.0_RP
          end if
       enddo

    enddo
    enddo

  ! tiny number, tiny mass
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
      do ic = 1, n_ctg       !aerosol category
      do ik = 1, n_kap(ic)   !kappa bin
      do is0 = 1, n_siz(ic)   !size bin
        if (QTRC0(k,i,j,QAES-1+it_procs2trans(ia_m0,is0,ik,ic))*DENS(k,i,j) < cleannumber) then
          do ia0 = 1, N_ATR
            QTRC0(k,i,j,QAES-1+it_procs2trans(ia0,is0,ik,ic)) = 0._RP !to save cpu time and avoid underflow
          enddo !ia0 (1:n_atr      )
        endif
      enddo !is (1:n_siz(ic)  )
      enddo !ik (1:n_kap(ic)  )
      enddo !ic (1:n_ctg      )
    enddo
    enddo
    enddo

    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
      QTRC1(k,i,j,iq) = QTRC0(k,i,j,iq) ! save
    enddo
    enddo
    enddo
    enddo

    !---- Calculate aerosol processs
    CN(:,:,:) = 0.0_RP
    CCN(:,:,:) = 0.0_RP
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE

       !--- aerosol_trans at initial time
       ! [xx/kg] -> [xx/m3]
       do it = 1, n_trans
            aerosol_procs(ia_trans2procs(it), &
                          is_trans2procs(it), &
                          ik_trans2procs(it), &
                          ic_trans2procs(it)) = QTRC1(k,i,j,QAES-1+it)*DENS(k,i,j)
            emis_procs(ia_trans2procs(it), &
                          is_trans2procs(it), &
                          ik_trans2procs(it), &
                          ic_trans2procs(it)) = EMIT(k,i,j,it)*DENS(k,i,j)
       enddo !it(1:n_trans)
       ! mixing ratio [kg/kg] -> concentration [ug/m3]
       conc_gas(1:GAS_CTG) &
               = QTRC1(k,i,j,QAEE-GAS_CTG+1:QAEE-GAS_CTG+GAS_CTG)*DENS(k,i,j)*1.E+9_RP

       emis_gas(1:GAS_CTG) = EMIT(k,i,j,QA_AE-GAS_CTG+IG_H2SO4:QA_AE-GAS_CTG+IG_CGAS)

       call aerosol_zerochem(           &
         TIME_DTSEC_ATMOS_PHY_AE,       & !--- in
         temp_ae(k,i,j),                & !--- in
         pres_ae(k,i,j),                & !--- in
         ssliq_ae(k,i,j),               & !--- in
         flag_npf, flag_cond, flag_coag,& !--- in
         aerosol_procs,                 & !--- inout
         conc_gas,                      & !--- inout
         emis_procs,                    & !--- out
         emis_gas,                      & !--- out
         aerosol_activ                  ) !--- out

! aerosol loss due to activation to cloud droplets
       if (flag_ccn_interactive) then
         do is0 = 1, n_siz(ic_mix)
         do ia0 = 1, N_ATR       !attributes
           aerosol_activ(ia0,is0,ik_out,ic_mix) = min(max(0._RP, aerosol_activ(ia0,is0,ik_out,ic_mix)), &
                                                                 aerosol_procs(ia0,is0,ik_out,ic_mix) )
           aerosol_procs(ia0,is0,ik_out,ic_mix) = &
           aerosol_procs(ia0,is0,ik_out,ic_mix) - aerosol_activ(ia0,is0,ik_out,ic_mix)
         enddo
         enddo
       endif !flag_ccn_interactive

! aerosol regeneration due to evaporation of cloud droplets
! using prescribed size parameters and to internal mixture category (ic_mix)
       if (flag_regeneration) then
         m0_reg = NREG(k,i,j)  !#/m3
!        m2_reg = m0_reg * dg_reg**2._RP * exp( 2.0_RP *(log(sg_reg)**2._RP)) !m2/m3
!        m3_reg = m0_reg * dg_reg**3._RP * exp( 4.5_RP *(log(sg_reg)**2._RP)) !m3/m3
         m2_reg = m0_reg * reg_factor_m2     !m2/m3
         m3_reg = m0_reg * reg_factor_m3     !m3/m3
         ms_reg = m3_reg * pi6 * conv_vl_ms  !ug/m3
         aerosol_procs(ia_m0,is0_reg,ik_out,ic_mix) = &
         aerosol_procs(ia_m0,is0_reg,ik_out,ic_mix) + m0_reg !#/m3
         aerosol_procs(ia_m2,is0_reg,ik_out,ic_mix) = &
         aerosol_procs(ia_m2,is0_reg,ik_out,ic_mix) + m2_reg !m2/m3
         aerosol_procs(ia_m3,is0_reg,ik_out,ic_mix) = &
         aerosol_procs(ia_m3,is0_reg,ik_out,ic_mix) + m3_reg !m3/m3
         aerosol_procs(ia_ms,is0_reg,ik_out,ic_mix) = &
         aerosol_procs(ia_ms,is0_reg,ik_out,ic_mix) + ms_reg !ug/m3
! additional attirbute to be added (ia_kp)
       endif !flag_regeneration

! diagnosed variables
       do is0 = 1, n_siz(ic_mix)
         CCN(k,i,j) = CCN(k,i,j) + aerosol_activ(ia_m0,is0,ik_out,ic_mix)
         CN(k,i,j)  = CN(k,i,j)  + aerosol_procs(ia_m0,is0,ik_out,ic_mix)
       enddo

!       call trans_ccn(aerosol_procs, aerosol_activ, t_ccn, t_cn,  &
!            n_ctg, n_kap_max, n_siz_max, N_ATR,         &
!            ic_mix, ia_m0, ia_m2, ia_m3, ik_out, n_siz, &
!            rnum_out, nbins_out)

!       CN(k,i,j) = t_cn
!       CCN(k,i,j) = t_ccn

       ! [xx/m3] -> [xx/kg]
       do ic = 1, n_ctg       !category
       do ik = 1, n_kap(ic)   !kappa bin
       do is0 = 1, n_siz(ic)   !size bin
       do ia0 = 1, N_ATR       !attributes
          QTRC1(k,i,j,QAES-1+it_procs2trans(ia0,is0,ik,ic)) = aerosol_procs(ia0,is0,ik,ic) / DENS(k,i,j)
       enddo !ia (1:N_ATR )
       enddo !is (1:n_siz(ic)  )
       enddo !ik (1:n_kap(ic)  )
       enddo !ic (1:n_ctg      )
       !  [ug/m3] -> mixing ratio [kg/kg]
       QTRC1(k,i,j,QAEE-GAS_CTG+1:QAEE-GAS_CTG+GAS_CTG) &
            = conc_gas(1:GAS_CTG) / DENS(k,i,j)*1.E-9_RP

       ! tiny number, tiny mass
       do ic = 1, n_ctg       !aerosol category
       do ik = 1, n_kap(ic)   !kappa bin
       do is0 = 1, n_siz(ic)   !size bin
         if (QTRC1(k,i,j,QAES-1+it_procs2trans(ia_m0,is0,ik,ic))*DENS(k,i,j) < cleannumber) then
           do ia0 = 1, N_ATR
             QTRC1(k,i,j,QAES-1+it_procs2trans(ia0,is0,ik,ic)) = 0._RP !to save cpu time and avoid underflow
           enddo !ia0 (1:n_atr      )
         endif
       enddo !is (1:n_siz(ic)  )
       enddo !ik (1:n_kap(ic)  )
       enddo !ic (1:n_ctg      )

       !--- Negative fixer
       do ic = 1, n_ctg       !aerosol category
       do ik = 1, n_kap(ic)   !kappa bin
       do is0 = 1, n_siz(ic)   !size bin
         if (QTRC(k,i,j,QAES-1+it_procs2trans(ia_m0,is0,ik,ic)) < 0.0_RP .or. &
             QTRC(k,i,j,QAES-1+it_procs2trans(ia_m2,is0,ik,ic)) < 0.0_RP .or. &
             QTRC(k,i,j,QAES-1+it_procs2trans(ia_m3,is0,ik,ic)) < 0.0_RP .or. &
             QTRC(k,i,j,QAES-1+it_procs2trans(ia_ms,is0,ik,ic)) < 0.0_RP .or. &
             QTRC(k,i,j,QAES-1+it_procs2trans(ia_kp,is0,ik,ic)) < 0.0_RP ) then
           QTRC(k,i,j,QAES-1+it_procs2trans(1:N_ATR,is0,ik,ic)) = 0.0_RP
         endif
       enddo
       enddo
       enddo

       ! for history
       total_aerosol_mass(k,i,j,:) = 0.0_RP
       total_aerosol_number(k,i,j,:) = 0.0_RP
       total_emit_aerosol_mass(k,i,j,:) = 0.0_RP
       total_emit_aerosol_number(k,i,j,:) = 0.0_RP
       do ic = 1, n_ctg
       do ik = 1, n_kap(ic)
       do is0 = 1, n_siz(ic)
           total_aerosol_mass       (k,i,j,ic) = total_aerosol_mass (k,i,j,ic) &
                                               + QTRC1(k,i,j,QAES-1+it_procs2trans(ia_ms,is0,ik,ic))
           total_aerosol_number     (k,i,j,ic) = total_aerosol_number     (k,i,j,ic) &
                                               + QTRC1(k,i,j,QAES-1+it_procs2trans(ia_m0,is0,ik,ic))
           total_emit_aerosol_mass  (k,i,j,ic) = total_emit_aerosol_mass  (k,i,j,ic) &
                                               + EMIT(k,i,j,it_procs2trans(ia_ms,is0,ik,ic))
           total_emit_aerosol_number(k,i,j,ic) = total_emit_aerosol_number(k,i,j,ic) &
                                               + EMIT(k,i,j,it_procs2trans(ia_m0,is0,ik,ic))
       enddo
       enddo
       enddo

    enddo
    enddo
    enddo

    do ic = 1, n_ctg
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'mass'
      call HIST_in( total_aerosol_mass  (:,:,:,ic), trim(ofilename), 'Total mass mixing ratio of aerosol', 'kg/kg' )
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'number'
      call HIST_in( total_aerosol_number(:,:,:,ic), trim(ofilename), 'Total number mixing ratio of aerosol', 'num/kg' )
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'mass_emit'
      call HIST_in( total_emit_aerosol_mass  (:,:,:,ic), trim(ofilename), 'Total mass mixing ratio of emitted aerosol', 'kg/kg' )
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'number_emit'
      call HIST_in( total_emit_aerosol_number(:,:,:,ic), trim(ofilename), 'Total number mixing ratio of emitted aerosol', 'num/kg' )
    enddo

    do ic = 1, n_ctg
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'mass'
      call HIST_in( total_aerosol_mass  (:,:,:,ic), trim(ofilename), 'Total mass mixing ratio of aerosol', 'kg/kg' )
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'number'
      call HIST_in( total_aerosol_number(:,:,:,ic), trim(ofilename), 'Total number mixing ratio of aerosol', 'num/kg' )
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'mass_emit'
      call HIST_in( total_emit_aerosol_mass  (:,:,:,ic), trim(ofilename), 'Total mass mixing ratio of emitted aerosol', 'kg/kg' )
      write(ofilename,'(a,a)') trim(ctg_name(ic)), 'number_emit'
      call HIST_in( total_emit_aerosol_number(:,:,:,ic), trim(ofilename), 'Total number mixing ratio of emitted aerosol', 'num/kg' )
    enddo

    call HIST_in( EMIT(:,:,:,QA_AE-GAS_CTG+IG_H2SO4), 'H2SO4_emit', 'Emission ratio of H2SO4 gas', 'ug/m3/s' )
    call HIST_in( EMIT(:,:,:,QA_AE-GAS_CTG+IG_CGAS),  'CGAS_emit',  'Emission ratio of Condensabule gas', 'ug/m3/s' )

    deallocate( aerosol_procs )
    deallocate( aerosol_activ )
    deallocate( emis_procs )
    deallocate( emis_gas )

    do iq = 1, QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
      RHOQ_t_AE(k,i,j,iq) = ( QTRC1(k,i,j,iq) - QTRC0(k,i,j,iq) ) * DENS(k,i,j) / TIME_DTSEC_ATMOS_PHY_AE
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_AE_kajino13
  !-----------------------------------------------------------------------------
  ! subroutine 4. aerosol_zerochem
  !-----------------------------------------------------------------------------
  subroutine aerosol_zerochem (        &
        deltt,                         & !--- in
        temp_k, pres_pa, super,        & !--- in
!        h2so4dt, c_ratio, c_kappa,     & !--- in
        flag_npf, flag_cond, flag_coag,& !--- in
        aerosol_procs,                 & !--- inout
        conc_gas,                      & !--- inout
        emis_procs,                    & !--- out
        emis_gas,                      & !--- out
        aerosol_activ                  )

    use scale_time, only: &
       TIME_NOWSEC, &
       TIME_STARTDAYSEC

    implicit none
    ! i/o variables
    real(DP),intent(in) :: deltt      ! delta t                                         [sec]
    real(RP),intent(in) :: temp_k     ! temperature                                     [K]
    real(RP),intent(in) :: pres_pa    ! pressure                                        [Pa]
    real(RP),intent(in) :: super      ! supersaturation                                 [-]
!    real(RP),intent(in) :: h2so4dt    ! h2so4 production rate                           [ug/m3/s]
!    real(RP),intent(in) :: c_ratio    ! ratio of condensable mass to h2so4 (after NPF)  [-]
!    real(RP),intent(in) :: c_kappa    ! hygroscopicity of condensable mass              [-]
    logical, intent(in) :: flag_npf   ! (on/off) new particle formation
    logical, intent(in) :: flag_cond  ! (on/off) condensation
    logical, intent(in) :: flag_coag  ! (on/off) coagulation
    real(RP),intent(inout) :: aerosol_procs(N_ATR,n_siz_max,n_kap_max,n_ctg)
    real(RP),intent(inout) :: conc_gas(GAS_CTG)
    real(RP),intent(out) :: aerosol_activ(N_ATR,n_siz_max,n_kap_max,n_ctg)
    real(RP),intent(in) :: emis_procs(N_ATR,n_siz_max,n_kap_max,n_ctg)
    real(RP),intent(in) :: emis_gas(GAS_CTG)
  ! local variables
    real(RP)            :: J_1nm      ! nucleation rate of 1nm particles [#/cm3/s]
    integer             :: ic_nuc     ! category  for 1nm new particles
    integer             :: ik_nuc     ! kappa bin for 1nm new particles
    integer             :: is_nuc     ! size bin  for 1nm new particles
    real(RP)            :: c_ratio, t_elaps
    real(RP)            :: chem_gas(GAS_CTG)

    chem_gas(IG_H2SO4) = h2so4dt
    chem_gas(IG_CGAS) = ocgasdt

    !--- convert unit of aerosol mass [ia=ia_ms]
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)  !size bin
      aerosol_procs(ia_ms,is0,ik,ic) = aerosol_procs(ia_ms,is0,ik,ic) * 1.E+9_RP ! [kg/m3] -> [ug/m3]
    enddo !is (1:n_siz(ic)  )
    enddo !ik (1:n_kap(ic)  )
    enddo !ic (1:n_ctg      )

  ! emission
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)   !size bin
    do ia0 = 1, N_ATR       !attributes (prognostic)
      aerosol_procs(ia0,is0,ik,ic) = aerosol_procs(ia0,is0,ik,ic) &
                                  + emis_procs(ia0,is0,ik,ic) * deltt
    enddo !ia0 (1:n_atr      )
    enddo !is (1:n_siz(ic)  )
    enddo !ik (1:n_kap(ic)  )
    enddo !ic (1:n_ctg      )

  ! tiny number, tiny mass
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)   !size bin
      if (aerosol_procs(ia_m0,is0,ik,ic) < cleannumber) then
        do ia0 = 1, N_ATR
          aerosol_procs(ia0,is0,ik,ic) = 0._RP !to save cpu time and avoid underflow
        enddo !ia0 (1:n_atr      )
      endif
    enddo !is (1:n_siz(ic)  )
    enddo !ik (1:n_kap(ic)  )
    enddo !ic (1:n_ctg      )

  ! update conc_h2so4
!    conc_h2so4 = conc_h2so4 + h2so4dt * deltt  ! [ug/m3]
    conc_gas(IG_H2SO4) = conc_gas(IG_H2SO4) + chem_gas(IG_H2SO4) * deltt  ! [ug/m3]
    conc_gas(IG_CGAS) = conc_gas(IG_CGAS) + chem_gas(IG_CGAS) * deltt  ! [ug/m3]

    conc_gas(IG_H2SO4) = conc_gas(IG_H2SO4) + emis_gas(IG_H2SO4) * deltt  ! [ug/m3]
    conc_gas(IG_CGAS) = conc_gas(IG_CGAS) + emis_gas(IG_CGAS) * deltt  ! [ug/m3]

    if( chem_gas(IG_H2SO4)+emis_gas(IG_H2SO4) /= 0.0_RP ) then
      c_ratio = ( chem_gas(IG_H2SO4)+emis_gas(IG_H2SO4)+chem_gas(IG_CGAS)+emis_gas(IG_CGAS) ) &
              / ( chem_gas(IG_H2SO4)+emis_gas(IG_H2SO4) )
    else
      c_ratio = 0.0_RP
    endif

  ! new particle formation
    t_elaps = TIME_NOWSEC - TIME_STARTDAYSEC
    J_1nm      = 0._RP
    ic_nuc     = ic_mix
    ik_nuc     = 1
    is_nuc     = 1
    if (flag_npf) then
      call PROF_rapstart('ATM_Aerosol_NPF',1)
      if(  t_elaps <= t_npf ) then !      call aerosol_nucleation(conc_h2so4, J_1nm)     !(i) conc_h2so4 (o) J_1nm
        call aerosol_nucleation(conc_gas(IG_H2SO4), J_1nm)     !(i) conc_h2so4 (o) J_1nm
      endif
      call PROF_rapend('ATM_Aerosol_NPF',1)
    else
      J_1nm      = 0._RP  !( new particle formation does not occur )
    endif !if flag_npf=.true.

  ! condensation
    if (flag_cond) then
      call PROF_rapstart('ATM_Aerosol_cond',1)
      call aerosol_condensation(J_1nm, temp_k, pres_pa, deltt,     &
             ic_nuc, ik_nuc, is_nuc, n_ctg, n_kap, n_siz, N_ATR,   &
             n_siz_max, n_kap_max, ia_m0, ia_m2, ia_m3, ia_ms,     &
  !          d_lw, d_ct, d_up, k_lw, k_ct, k_up,                   &
             d_lw, d_ct, d_up,                                     &
             conc_gas(IG_H2SO4), c_ratio, aerosol_procs            )
      call PROF_rapend('ATM_Aerosol_cond',1)
    endif !if flag_cond=.true.

  ! coagulation
    if (flag_coag) then
      call PROF_rapstart('ATM_Aerosol_coag',1)
      call aerosol_coagulation(deltt, temp_k, pres_pa,                &
             mcomb,is_i,is_j,is_k,ik_i,ik_j,ik_k,ic_i,ic_j,ic_k,      &
             N_ATR,n_siz_max,n_kap_max,n_ctg,ia_m0,ia_m2,ia_m3,ia_ms, &
             n_siz,n_kap,d_lw,d_ct,d_up,aerosol_procs)
      call PROF_rapend('ATM_Aerosol_coag',1)
    endif !if flag_coag=.true.

    call aerosol_activation(c_kappa, super, temp_k, ia_m0, ia_m2, ia_m3, &
                            N_ATR,n_siz_max,n_kap_max,n_ctg,n_siz,n_kap, &
                            d_ct,aerosol_procs, aerosol_activ)

!   conc_gas(IG_CGAS) = conc_gas(IG_H2SO4)*( c_ratio-1.0_RP )

    !--- convert unit of aerosol mass [ia=ia_ms]
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)   !size bin
      aerosol_procs(ia_ms,is0,ik,ic) = aerosol_procs(ia_ms,is0,ik,ic) * 1.E-9_RP ! [ug/m3] -> [kg/m3]
    enddo !is (1:n_siz(ic)  )
    enddo !ik (1:n_kap(ic)  )
    enddo !ic (1:n_ctg      )

    return
  end subroutine aerosol_zerochem
  !-----------------------------------------------------------------------------
  ! subroutine 1. aerosol_nucleation
  !----------------------------------------------------------------------------------------
  subroutine aerosol_nucleation(conc_h2so4, J_1nm)        ! (i) conc_h2so4 (i/o) J_1nm
    implicit none
    !i/o variables
    real(RP), intent(in)    :: conc_h2so4      !H2SO4 concentration   [ug/m3]
    real(RP), intent(inout) :: J_1nm           !nucleation rate of 1nm particles [#/cm3/s]
    !local variables
    real(RP)                :: conc_num_h2so4  !H2SO4 number concentration   [#/cm3]
  !  real(RP), parameter     :: logk = -12.4_RP !constant coefficient for kinetic nucleation [-]
  ! real(RP), parameter     :: logk = -11.4_RP !constant coefficient for kinetic nucleation [-]

    conc_num_h2so4 = conc_h2so4 * conv_m2n     ![ug/m3] -> [#/cm3]

    !emperical formula of new particle formation rate [Kuang et al., 2008]
    J_1nm  = (10._RP**(logk_aenucl)) * conc_num_h2so4 ** 2._RP

    return
  end subroutine aerosol_nucleation
  !----------------------------------------------------------------------------------------
  ! subroutine 2. aerosol_condensation
  !----------------------------------------------------------------------------------------
  subroutine aerosol_condensation(J_1nm, temp_k, pres_pa, deltt,     &
               ic_nuc, ik_nuc, is_nuc, n_ctg, n_kap, n_siz, n_atr,   &
               n_siz_max, n_kap_max, ia_m0, ia_m2, ia_m3, ia_ms,     &
  !              d_lw, d_ct, d_up, k_lw, k_ct, k_up,                   &
               d_lw, d_ct, d_up,                                     &
               conc_h2so4, c_ratio, aerosol_procs )

    implicit none
    !i/o variables
    real(RP),intent(in)    :: J_1nm      ! nucleation rate of 1nm particles [#/cm3/s]
    real(RP),intent(in)    :: temp_k     ! temperature                      [K]
    real(RP),intent(in)    :: pres_pa    ! pressure                         [Pa]
    real(DP),intent(in)    :: deltt      ! delta t                          [sec]
    integer, intent(in)    :: ic_nuc     ! category  for 1nm new particles
    integer, intent(in)    :: ik_nuc     ! kappa bin for 1nm new particles
    integer, intent(in)    :: is_nuc     ! size bin  for 1nm new particles
    integer, intent(in)    :: n_ctg      ! number of aerosol category
    integer, intent(in)    :: n_kap(n_ctg)      ! number of kappa bins
    integer, intent(in)    :: n_siz(n_ctg)      ! number of size bins
    integer, intent(in)    :: n_atr      ! number of aerosol category
    integer, intent(in)    :: n_siz_max  ! max of n_siz
    integer, intent(in)    :: n_kap_max  ! max of n_kap
    integer, intent(in)    :: ia_m0      !
    integer, intent(in)    :: ia_m2      !
    integer, intent(in)    :: ia_m3      !
    integer, intent(in)    :: ia_ms      !
    real(RP),intent(in)    :: c_ratio    ! ratio of condensable mass to h2so4 (after NPF)  [-]
    real(RP),intent(in)   :: d_lw(n_siz_max,n_ctg), d_ct(n_siz_max,n_ctg), d_up(n_siz_max,n_ctg)
  ! real(RP), dimension(n_kap_max,n_ctg), intent(in)    :: k_lw, k_ct, k_up
    real(RP),intent(inout) :: conc_h2so4 !concentration [ug/m3]
    real(RP),intent(inout) :: aerosol_procs(n_atr,n_siz_max,n_kap_max,n_ctg)
  !local variables
    real(RP), parameter  :: alpha = 0.1_RP       ! accomodation coefficient
    real(RP), parameter  :: cour  = 0.5_RP       ! courant number for condensation
    integer              :: isplt                ! number of time splitted
    real(RP)             :: dtsplt               ! splitted time step
    real(RP) :: sq_cbar_h2so4                   ! square of molecular speed of H2SO4 gas [m2/s2]
    real(RP) :: cbar_h2so4                      ! molecular speed of H2SO4 gas           [m/s]
    real(RP) :: dv_h2so4                        ! diffusivity of H2SO4 gas               [m2/s]
    real(RP) :: drive                           ! driving force of H2SO4 gas       [m3SO4/m3]
    real(RP) :: dm0dt_npf, dm2dt_npf, dm3dt_npf, dmsdt_npf     ! dMk/dt for new particle formation
    real(RP) :: dm0dt_cnd(n_siz_max,n_kap_max,n_ctg) ! dM0/dt for condensation/evaporation
    real(RP) :: dm2dt_cnd(n_siz_max,n_kap_max,n_ctg) ! dM2/dt for condensation/evaporation
    real(RP) :: dm3dt_cnd(n_siz_max,n_kap_max,n_ctg) ! dM3/dt for condensation/evaporation
    real(RP) :: dmsdt_cnd(n_siz_max,n_kap_max,n_ctg) ! dMs/dt for condensation/evaporation
    real(RP) :: rm0_pls(n_siz_max), rm0_mns(n_siz_max) ! used for moving center calc.
    real(RP) :: rm2_pls(n_siz_max), rm2_mns(n_siz_max) ! used for moving center calc.
    real(RP) :: rm3_pls(n_siz_max), rm3_mns(n_siz_max) ! used for moving center calc.
    real(RP) :: rms_pls(n_siz_max), rms_mns(n_siz_max) ! used for moving center calc.
    real(RP) :: m0t,m1t,m2t,m3t,dgt,sgt,dm2
    real(RP) :: gnc2,gnc3,gfm2,gfm3,harm2,harm3
    real(RP) :: lossrate,tmps6
    integer  :: ic, ik, is0, i, is1, is2

! nothing happens --> in future, c_ratio is removed and condensable gas should be used
    if (conc_h2so4 <= 0.0_RP) return
! nothing happens --> in future, c_ratio is removed and condensable gas should be used

    drive         = conc_h2so4 * mwrat_s6 * conv_ms_vl  ![ugH2SO4/m3]=>volume [m3SO4/m3]
    sq_cbar_h2so4 = 8.0_RP*rgas*temp_k/(pi*mwh2so4*1.E-3_RP)
    cbar_h2so4    = sqrt( sq_cbar_h2so4 )
    dv_h2so4      = diffsulf*(stdatmpa/pres_pa)*(temp_k/stdtemp)**1.75_RP

    dm0dt_npf = 0.0_RP
    dm2dt_npf = 0.0_RP
    dm3dt_npf = 0.0_RP
    dmsdt_npf = 0.0_RP
    dm0dt_cnd = 0.0_RP
    dm2dt_cnd = 0.0_RP
    dm3dt_cnd = 0.0_RP
    dmsdt_cnd = 0.0_RP

    !(npf rate)
    dm0dt_npf = J_1nm * 1.E6_RP                  ! [#/m3/s]
    dm2dt_npf = dm0dt_npf * 1.E-18_RP            ! [m2/m3/s]
    dm3dt_npf = dm0dt_npf * 1.E-27_RP            ! [m3/m3/s]
    dmsdt_npf = dm3dt_npf * pi6 * conv_vl_ms     ! [ug/m3/s]

   !(condensation rate)
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)   !size bin

      dm2dt_cnd(is0,ik,ic) = 0._RP
      dm3dt_cnd(is0,ik,ic) = 0._RP

      if (aerosol_procs(ia_m0,is0,ik,ic) &
         *aerosol_procs(ia_m2,is0,ik,ic) &
         *aerosol_procs(ia_m3,is0,ik,ic) > 0._RP) then
        m0t = aerosol_procs(ia_m0,is0,ik,ic)
        m2t = aerosol_procs(ia_m2,is0,ik,ic)
        m3t = aerosol_procs(ia_m3,is0,ik,ic)
        call diag_ds(m0t,m2t,m3t,dgt,sgt,dm2)
        if (dgt <= 0._RP) dgt = d_ct(is0,ic)
        if (sgt <= 0._RP) sgt = 1.3_RP
        m1t = m0t*dgt*exp(0.5_RP*(log(sgt)**2._RP)) ![m/m3]
        !=4kDvM_k-2
        gnc2 =  8._RP * dv_h2so4 * m0t           !m2/s*#/m3  = m2/m3/s
        gnc3 = 12._RP * dv_h2so4 * m1t           !m2/s*m/m3  = m3/m3/s
        != k/2acM_k-1
        gfm2 =          alpha * cbar_h2so4 * m1t !m/s *m/m3  = m2/m3/s
        gfm3 = 1.5_RP * alpha * cbar_h2so4 * m2t !m/s *m2/m3 = m3/m3/s
        != harmonic mean approach
        harm2 = gnc2 * gfm2 / ( gnc2 + gfm2 )    !m2/m3/s
        harm3 = gnc3 * gfm3 / ( gnc3 + gfm3 )    !m3/m3/s

        dm2dt_cnd(is0,ik,ic) = harm2*drive                         !m2/m3/s
        dm3dt_cnd(is0,ik,ic) = harm3*drive                         !m3/m3/s
        dmsdt_cnd(is0,ik,ic) = dm3dt_cnd(is0,ik,ic)*pi6*conv_vl_ms !ug/m3/s

      endif !aerosol number > 0

    enddo
    enddo
    enddo

  ! coagulation
    !=time split
    lossrate = dmsdt_npf

    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, n_siz(ic)   !size bin
      lossrate = lossrate + dmsdt_cnd(is0,ik,ic)
    enddo
    enddo
    enddo

    if (lossrate<=0._RP) return !nothing happens

    isplt  = int( (lossrate*deltt) / (cour*conc_h2so4) )+1
    dtsplt = deltt/float(isplt)

    !== NPF & condensation calculation
    loop_split: do i = 1, isplt
      tmps6      = conc_h2so4
      conc_h2so4 = conc_h2so4 - (lossrate * dtsplt) * mwrat_s6_i !ugH2SO4/m3

      if (conc_h2so4 < 0._RP) then
        dtsplt     = tmps6 * mwrat_s6 / lossrate
        conc_h2so4 = 0._RP
      endif !conc_h2so4 < 0

      !(npf)
      aerosol_procs(ia_m0,is_nuc,ik_nuc,ic_nuc) = & ! #/m3
          aerosol_procs(ia_m0,is_nuc,ik_nuc,ic_nuc) + dm0dt_npf * dtsplt
      aerosol_procs(ia_m2,is_nuc,ik_nuc,ic_nuc) = & ! m2/m3
          aerosol_procs(ia_m2,is_nuc,ik_nuc,ic_nuc) + dm2dt_npf * dtsplt
      aerosol_procs(ia_m3,is_nuc,ik_nuc,ic_nuc) = & ! m3/m3
          aerosol_procs(ia_m3,is_nuc,ik_nuc,ic_nuc) + dm3dt_npf * dtsplt
      aerosol_procs(ia_ms,is_nuc,ik_nuc,ic_nuc) = & ! ug/m3
          aerosol_procs(ia_ms,is_nuc,ik_nuc,ic_nuc) + dmsdt_npf * dtsplt

     !(condensation)
      do ic = 1, n_ctg       !aerosol category
      do ik = 1, n_kap(ic)   !kappa bin
      do is0 = 1, n_siz(ic)   !size bin

        aerosol_procs(ia_m2,is0,ik,ic) = & ! m2/m3
            aerosol_procs(ia_m2,is0,ik,ic) + dm2dt_cnd(is0,ik,ic) * dtsplt & !m2/m3
                                          * c_ratio                        !additional mass for condensation
        aerosol_procs(ia_m3,is0,ik,ic) = & ! m3/m3
            aerosol_procs(ia_m3,is0,ik,ic) + dm3dt_cnd(is0,ik,ic) * dtsplt & !m3/m3
                                          * c_ratio                        !additional mass for condensation
        aerosol_procs(ia_ms,is0,ik,ic) = & ! ug/m3
            aerosol_procs(ia_ms,is0,ik,ic) + dmsdt_cnd(is0,ik,ic) * dtsplt & !ug/m3
                                          * c_ratio                        !additional mass for condensation

      enddo
      enddo
      enddo

!      if (conc_h2so4 <= 0._RP) goto 777
      if (conc_h2so4 <= 0.0_RP) exit loop_split

    enddo loop_split

!    777 continue

    conc_h2so4 = max(conc_h2so4, 0._RP)

  !== moving center
    do ic = 1, n_ctg       !aerosol category
    do ik = 1, n_kap(ic)   !kappa bin

      rm0_pls = 0._RP
      rm2_pls = 0._RP
      rm3_pls = 0._RP
      rms_pls = 0._RP
      rm0_mns = 0._RP
      rm2_mns = 0._RP
      rm3_mns = 0._RP
      rms_mns = 0._RP

      do is1 = 1, n_siz(ic)   !size bin

        if (aerosol_procs(ia_m0,is1,ik,ic) &
           *aerosol_procs(ia_m2,is1,ik,ic) &
           *aerosol_procs(ia_m3,is1,ik,ic) > 0._RP) then
          m0t = aerosol_procs(ia_m0,is1,ik,ic)
          m2t = aerosol_procs(ia_m2,is1,ik,ic)
          m3t = aerosol_procs(ia_m3,is1,ik,ic)
          call diag_ds(m0t,m2t,m3t,dgt,sgt,dm2)
          if (dgt <= 0._RP) dgt = d_ct(is1,ic)
          if (sgt <= 0._RP) sgt = 1.3_RP

          do is2 = 1, n_siz(ic)
            if (dgt >= d_lw(is2,ic) .AND. dgt < d_up(is2,ic)) then !moving center
              rm0_pls(is2) = rm0_pls(is2) + aerosol_procs(ia_m0,is1,ik,ic) !is2 <=
              rm0_mns(is1) =                aerosol_procs(ia_m0,is1,ik,ic) !is1 =>
              rm2_pls(is2) = rm2_pls(is2) + aerosol_procs(ia_m2,is1,ik,ic) !is2 <=
              rm2_mns(is1) =                aerosol_procs(ia_m2,is1,ik,ic) !is1 =>
              rm3_pls(is2) = rm3_pls(is2) + aerosol_procs(ia_m3,is1,ik,ic) !is2 <=
              rm3_mns(is1) =                aerosol_procs(ia_m3,is1,ik,ic) !is1 =>
              rms_pls(is2) = rms_pls(is2) + aerosol_procs(ia_ms,is1,ik,ic) !is2 <=
              rms_mns(is1) =                aerosol_procs(ia_ms,is1,ik,ic) !is1 =>
              exit
            endif !d_lw(is2) < dgt < d_up(is2)
          enddo

        endif !aerosol number > 0

      enddo

      do is0 = 1, n_siz(ic)
        aerosol_procs(ia_m0,is0,ik,ic) = aerosol_procs(ia_m0,is0,ik,ic) &
                                      + rm0_pls(is0) - rm0_mns(is0)
        aerosol_procs(ia_m2,is0,ik,ic) = aerosol_procs(ia_m2,is0,ik,ic) &
                                      + rm2_pls(is0) - rm2_mns(is0)
        aerosol_procs(ia_m3,is0,ik,ic) = aerosol_procs(ia_m3,is0,ik,ic) &
                                      + rm3_pls(is0) - rm3_mns(is0)
        aerosol_procs(ia_ms,is0,ik,ic) = aerosol_procs(ia_ms,is0,ik,ic) &
                                      + rms_pls(is0) - rms_mns(is0)
      enddo !is(1:n_siz(ic))

    enddo
    enddo

    return
  end subroutine aerosol_condensation
  !----------------------------------------------------------------------------------------
  ! subroutine 3. aerosol_coagulation
  !----------------------------------------------------------------------------------------
  subroutine aerosol_coagulation(deltt, temp_k, pres_pa,                &
              mcomb,is_i,is_j,is_k,ik_i,ik_j,ik_k,ic_i,ic_j,ic_k,      &
              n_atr,n_siz_max,n_kap_max,n_ctg,ia_m0,ia_m2,ia_m3,ia_ms, &
              n_siz,n_kap,d_lw,d_ct,d_up,aerosol_procs)
      implicit none
      !i/o variables
      real(DP),intent(in) :: deltt      ! delta t        [sec]
      real(RP),intent(in) :: temp_k     ! temperature    [K]
      real(RP),intent(in) :: pres_pa    ! pressure       [Pa]
      integer, intent(in) :: mcomb !combinations of sections
      integer, intent(in) :: is_i(mcomb), is_j(mcomb), is_k(mcomb)
      integer, intent(in) :: ik_i(mcomb), ik_j(mcomb), ik_k(mcomb)
      integer, intent(in) :: ic_i(mcomb), ic_j(mcomb), ic_k(mcomb)
      integer, intent(in) :: n_ctg      ! number of aerosol category
      integer, intent(in) :: n_atr      ! number of aerosol category
      integer, intent(in) :: n_siz_max  ! max of n_siz
      integer, intent(in) :: n_kap_max  ! max of n_kap
      integer, intent(in) :: ia_m0      !
      integer, intent(in) :: ia_m2      !
      integer, intent(in) :: ia_m3      !
      integer, intent(in) :: ia_ms      !
      integer, intent(in) :: n_siz(n_ctg)      !
      integer, intent(in) :: n_kap(n_ctg)      !
      real(RP),intent(in):: d_lw(n_siz_max,n_ctg), d_ct(n_siz_max,n_ctg), d_up(n_siz_max,n_ctg)
      real(RP),intent(inout) :: aerosol_procs(n_atr,n_siz_max,n_kap_max,n_ctg)
      !local variables
      real(RP) :: dt_m0(n_siz_max,n_kap_max,n_ctg)
      real(RP) :: dt_m3(n_siz_max,n_kap_max,n_ctg)
      real(RP) :: dt_m6(n_siz_max,n_kap_max,n_ctg)
      real(RP) :: dt_ms(n_siz_max,n_kap_max,n_ctg)
      real(RP) :: sixth(n_siz_max,n_kap_max,n_ctg)
      real(RP) :: rm0_pls(n_siz_max), rm0_mns(n_siz_max) ! used for moving center calc.
      real(RP) :: rm2_pls(n_siz_max), rm2_mns(n_siz_max) ! used for moving center calc.
      real(RP) :: rm3_pls(n_siz_max), rm3_mns(n_siz_max) ! used for moving center calc.
      real(RP) :: rms_pls(n_siz_max), rms_mns(n_siz_max) ! used for moving center calc.
      real(RP) :: visair          ! viscosity of air [Pa s]=[kg/m/s]
      real(RP) :: lambda          ! mean free path of air molecules        [cm]
      real(RP) :: rkfm,rknc
      real(RP) :: m0i,m0j,m2i,m2j,m3i,m3j,msi,msj,dgi,dgj,sgi,sgj,dm2,m6i,m6j
      real(RP) :: dm0i,dm0j,dm0k,dm3i,dm3j,dm3k,dm6i,dm6j,dm6k,dmsi,dmsj,dmsk
      real(RP) :: m0t,m3t,m6t,dgt,sgt,m2t
      integer  :: mc, ic, ik, is0, is1, is2

      dt_m0(:,:,:) = 0._RP
      dt_m3(:,:,:) = 0._RP
      dt_m6(:,:,:) = 0._RP
      dt_ms(:,:,:) = 0._RP
      sixth(:,:,:) = 0._RP

      visair=c1*temp_k**c2/(1._RP+c3/temp_k)                ! viscosity of air [Pa s]=[kg/m/s]
      rkfm=(3._RP*boltz*temp_k /rho_kg)**0.5_RP             ![J/K*K*m3/kg]=[kg*m2/s2*m3/kg]=[m2.5/s]
      rknc= 2._RP*boltz*temp_k /(3._RP*visair)              ![J*m*s/kg]=[kg*m2/s2*m*s/kg]=[m3/s]
      c4 = pi*boltz*temp_k/(8._RP*mwair/avo*1.E-3_RP)       !
      lambda = visair/(0.499_RP*pres_pa)*c4**0.5_RP*1.E2_RP ! mean free path of air molecules   [cm]


    !--- 555 coagulation combination rule loop
        do mc = 1, mcomb
          m0i = aerosol_procs(ia_m0,is_i(mc),ik_i(mc),ic_i(mc))
          m0j = aerosol_procs(ia_m0,is_j(mc),ik_j(mc),ic_j(mc))
          m2i = aerosol_procs(ia_m2,is_i(mc),ik_i(mc),ic_i(mc))
          m2j = aerosol_procs(ia_m2,is_j(mc),ik_j(mc),ic_j(mc))
          m3i = aerosol_procs(ia_m3,is_i(mc),ik_i(mc),ic_i(mc))
          m3j = aerosol_procs(ia_m3,is_j(mc),ik_j(mc),ic_j(mc))
          msi = aerosol_procs(ia_ms,is_i(mc),ik_i(mc),ic_i(mc))
          msj = aerosol_procs(ia_ms,is_j(mc),ik_j(mc),ic_j(mc))

          if (m0i*m2i*m3i <= 0._RP .OR. m0j*m2j*m3j <= 0._RP) cycle

          call diag_ds(m0i,m2i,m3i,dgi,sgi,dm2)
          if (dgi <= 0._RP) dgi = d_ct(is_i(mc),ic_i(mc))
          if (sgi <= 0._RP) sgi = 1.3_RP
          call diag_ds(m0j,m2j,m3j,dgj,sgj,dm2)
          if (dgj <= 0._RP) dgj = d_ct(is_j(mc),ic_j(mc))
          if (sgj <= 0._RP) sgj = 1.3_RP

          m6i = m0i*dgi**6._RP*exp(18._RP*(log(sgi)**2._RP))
          m6j = m0j*dgj**6._RP*exp(18._RP*(log(sgj)**2._RP))
          sixth(is_i(mc),ik_i(mc),ic_i(mc)) = m6i
          sixth(is_j(mc),ik_j(mc),ic_j(mc)) = m6j

          !intra sectional coagulation
          if ( is_i(mc) == is_j(mc) .AND. &
               ik_i(mc) == ik_j(mc) .AND. &
               ic_i(mc) == ic_j(mc) ) then
            call aero_intra(m0i,  m3i,          & !i
                            dm0i, dm3i, dm6i,   & !o
                            dgi,  sgi,  lambda, & !i
                            rknc, rkfm, deltt   ) !i
            dm0j = dm0i
            dm3j = dm3i
            dm6j = dm6i
            dm0k = -dm0i
            dm3k = -dm3i
            dm6k = -dm6i
          !inter sectional coagulation
          else
            call aero_inter(m0i ,m3i ,m6i ,      & !input unchanged
                            m0j ,m3j ,m6j ,      & !input unchanged
                            dm0i,dm3i,dm6i,dmsi, & !output
                            dm0j,dm3j,dm6j,dmsj, & !output
                            dm0k,dm3k,dm6k,dmsk, & !output
                            dgi ,sgi, rhod_ae,   & !input unchanged
                            dgj ,sgj, rhod_ae,   & !input unchanged
                            rknc, rkfm, lambda,  & !input unchanged
                            deltt                ) !input unchanged
          endif

          dt_m0(is_i(mc),ik_i(mc),ic_i(mc)) = dt_m0(is_i(mc),ik_i(mc),ic_i(mc)) + dm0i
          dt_m0(is_j(mc),ik_j(mc),ic_j(mc)) = dt_m0(is_j(mc),ik_j(mc),ic_j(mc)) + dm0j
          dt_m0(is_k(mc),ik_k(mc),ic_k(mc)) = dt_m0(is_k(mc),ik_k(mc),ic_k(mc)) + dm0k
          dt_m3(is_i(mc),ik_i(mc),ic_i(mc)) = dt_m3(is_i(mc),ik_i(mc),ic_i(mc)) + dm3i
          dt_m3(is_j(mc),ik_j(mc),ic_j(mc)) = dt_m3(is_j(mc),ik_j(mc),ic_j(mc)) + dm3j
          dt_m3(is_k(mc),ik_k(mc),ic_k(mc)) = dt_m3(is_k(mc),ik_k(mc),ic_k(mc)) + dm3k
          dt_m6(is_i(mc),ik_i(mc),ic_i(mc)) = dt_m6(is_i(mc),ik_i(mc),ic_i(mc)) + dm6i
          dt_m6(is_j(mc),ik_j(mc),ic_j(mc)) = dt_m6(is_j(mc),ik_j(mc),ic_j(mc)) + dm6j
          dt_m6(is_k(mc),ik_k(mc),ic_k(mc)) = dt_m6(is_k(mc),ik_k(mc),ic_k(mc)) + dm6k
          dt_ms(is_i(mc),ik_i(mc),ic_i(mc)) = dt_ms(is_i(mc),ik_i(mc),ic_i(mc)) + dmsi
          dt_ms(is_j(mc),ik_j(mc),ic_j(mc)) = dt_ms(is_j(mc),ik_j(mc),ic_j(mc)) + dmsj
          dt_ms(is_k(mc),ik_k(mc),ic_k(mc)) = dt_ms(is_k(mc),ik_k(mc),ic_k(mc)) + dmsk

        enddo !555 continue !mc(1:mcomb)


        !--- redistribution
        do ic = 1, n_ctg       !aerosol category
        do ik = 1, n_kap(ic)   !kappa bin
        do is0 = 1, n_siz(ic)   !size bin
          if (aerosol_procs(ia_m0,is0,ik,ic) > 0._RP .AND. &                 !to avoid divided by zero
              dt_m0(is0,ik,ic)/aerosol_procs(ia_m0,is0,ik,ic) < 1._RP ) then !to avoid overflow in exp
            aerosol_procs(ia_m0,is0,ik,ic)=aerosol_procs(ia_m0,is0,ik,ic) &
                                          *exp( dt_m0(is0,ik,ic)/aerosol_procs(ia_m0,is0,ik,ic) )
          else !aerosol_procs = 0
            aerosol_procs(ia_m0,is0,ik,ic)=aerosol_procs(ia_m0,is0,ik,ic) + dt_m0(is0,ik,ic)
          endif
          if (aerosol_procs(ia_m3,is0,ik,ic) > 0._RP .AND. &                 !to avoid divided by zero
              dt_m3(is0,ik,ic)/aerosol_procs(ia_m3,is0,ik,ic) < 1._RP ) then !to avoid overflow in exp
            aerosol_procs(ia_m3,is0,ik,ic)=aerosol_procs(ia_m3,is0,ik,ic) &
                                          *exp( dt_m3(is0,ik,ic)/aerosol_procs(ia_m3,is0,ik,ic) )
          else !aerosol_procs = 0
            aerosol_procs(ia_m3,is0,ik,ic)=aerosol_procs(ia_m3,is0,ik,ic) + dt_m3(is0,ik,ic)
          endif
          if (sixth(is0,ik,ic)               > 0._RP .AND. &                 !to avoid divided by zero
              dt_m6(is0,ik,ic)/sixth(is0,ik,ic) < 1._RP               ) then !to avoid overflow in exp
            sixth(is0,ik,ic)              =sixth(is0,ik,ic) &
                                          *exp( dt_m6(is0,ik,ic)/sixth(is0,ik,ic)               )
          else !aerosol_procs = 0
            sixth(is0,ik,ic)              =sixth(is0,ik,ic)               + dt_m6(is0,ik,ic)
          endif
          if (aerosol_procs(ia_ms,is0,ik,ic) > 0._RP .AND. &                 !to avoid divided by zero
              dt_ms(is0,ik,ic)/aerosol_procs(ia_ms,is0,ik,ic) < 1._RP ) then !to avoid overflow in exp
            aerosol_procs(ia_ms,is0,ik,ic)=aerosol_procs(ia_ms,is0,ik,ic) &
                                          *exp( dt_ms(is0,ik,ic)/aerosol_procs(ia_ms,is0,ik,ic) )
          else !aerosol_procs = 0
            aerosol_procs(ia_ms,is0,ik,ic)=aerosol_procs(ia_ms,is0,ik,ic) + dt_ms(is0,ik,ic)
          endif
          m0t = aerosol_procs(ia_m0,is0,ik,ic)
          m3t = aerosol_procs(ia_m3,is0,ik,ic)
          m6t = sixth              (is0,ik,ic)
          call diag_ds6(m0t,m3t,m6t,   &
                        m2t,dgt,sgt,dm2)
          if (dgt <= 0._RP) dgt = d_ct(is0,ic)
          if (sgt <= 0._RP) sgt = 1.3_RP
          if (m3t == 0._RP) aerosol_procs(ia_ms,is0,ik,ic)=0._RP
          aerosol_procs(ia_m2,is0,ik,ic)=m2t
          ! Y.Sato added
          aerosol_procs(ia_m0,is0,ik,ic)=m0t
          aerosol_procs(ia_m3,is0,ik,ic)=m3t
        enddo
        enddo
        enddo

    !== moving center
        do ic = 1, n_ctg       !aerosol category
        do ik = 1, n_kap(ic)   !kappa bin

          rm0_pls = 0._RP
          rm2_pls = 0._RP
          rm3_pls = 0._RP
          rms_pls = 0._RP
          rm0_mns = 0._RP
          rm2_mns = 0._RP
          rm3_mns = 0._RP
          rms_mns = 0._RP

          do is1 = 1, n_siz(ic)   !size bin

            if (aerosol_procs(ia_m0,is1,ik,ic) &
               *aerosol_procs(ia_m2,is1,ik,ic) &
               *aerosol_procs(ia_m3,is1,ik,ic) > 0._RP) then
              m0t = aerosol_procs(ia_m0,is1,ik,ic)
              m2t = aerosol_procs(ia_m2,is1,ik,ic)
              m3t = aerosol_procs(ia_m3,is1,ik,ic)
              call diag_ds(m0t,m2t,m3t,dgt,sgt,dm2)
              if (dgt <= 0._RP) dgt = d_ct(is1,ic)
              if (sgt <= 0._RP) sgt = 1.3_RP

              do is2 = 1, n_siz(ic)
                if (dgt >= d_lw(is2,ic) .AND. dgt < d_up(is2,ic)) then !moving center
                  rm0_pls(is2) = rm0_pls(is2) + aerosol_procs(ia_m0,is1,ik,ic) !is2 <=
                  rm0_mns(is1) =                aerosol_procs(ia_m0,is1,ik,ic) !is1 =>
                  rm2_pls(is2) = rm2_pls(is2) + aerosol_procs(ia_m2,is1,ik,ic) !is2 <=
                  rm2_mns(is1) =                aerosol_procs(ia_m2,is1,ik,ic) !is1 =>
                  rm3_pls(is2) = rm3_pls(is2) + aerosol_procs(ia_m3,is1,ik,ic) !is2 <=
                  rm3_mns(is1) =                aerosol_procs(ia_m3,is1,ik,ic) !is1 =>
                  rms_pls(is2) = rms_pls(is2) + aerosol_procs(ia_ms,is1,ik,ic) !is2 <=
                  rms_mns(is1) =                aerosol_procs(ia_ms,is1,ik,ic) !is1 =>
                  exit
                endif !d_lw(is2) < dgt < d_up(is2)
              enddo

            endif !aerosol number > 0

          enddo

          do is0 = 1, n_siz(ic)
            aerosol_procs(ia_m0,is0,ik,ic) = aerosol_procs(ia_m0,is0,ik,ic) &
                                          + rm0_pls(is0) - rm0_mns(is0)
            aerosol_procs(ia_m2,is0,ik,ic) = aerosol_procs(ia_m2,is0,ik,ic) &
                                          + rm2_pls(is0) - rm2_mns(is0)
            aerosol_procs(ia_m3,is0,ik,ic) = aerosol_procs(ia_m3,is0,ik,ic) &
                                          + rm3_pls(is0) - rm3_mns(is0)
            aerosol_procs(ia_ms,is0,ik,ic) = aerosol_procs(ia_ms,is0,ik,ic) &
                                          + rms_pls(is0) - rms_mns(is0)
          enddo !is(1:n_siz(ic))

        enddo
        enddo


      return
  end subroutine aerosol_coagulation
  !----------------------------------------------------------------------------------------
  ! subroutine 4. aerosol_activation
  !   Abdul-Razzak et al.,   JGR, 1998 [AR98]
  !   Abdul-Razzak and Ghan, JGR, 2000 [AR00]
  !----------------------------------------------------------------------------------------
  subroutine aerosol_activation(c_kappa, super, temp_k, ia_m0, ia_m2, ia_m3, &
                                n_atr,n_siz_max,n_kap_max,n_ctg,n_siz,n_kap, &
                                d_ct, aerosol_procs, aerosol_activ)

      implicit none
    !i/o variables
      real(RP),intent(in) :: super      ! supersaturation                                 [-]
      real(RP),intent(in) :: c_kappa    ! hygroscopicity of condensable mass              [-]
      real(RP),intent(in) :: temp_k     ! temperature
      integer, intent(in) :: ia_m0, ia_m2, ia_m3
      integer, intent(in) :: n_atr
      integer, intent(in) :: n_siz_max
      integer, intent(in) :: n_kap_max
      integer, intent(in) :: n_ctg
      real(RP),intent(in) :: d_ct(n_siz_max,n_ctg)
      real(RP) :: aerosol_procs(n_atr,n_siz_max,n_kap_max,n_ctg)
      real(RP) :: aerosol_activ(n_atr,n_siz_max,n_kap_max,n_ctg)
      integer, intent(in) :: n_siz(n_ctg), n_kap(n_ctg)
    !local variables
      real(RP),parameter :: two3 = 2._RP/3._RP
      real(RP),parameter :: rt2  = sqrt(2._RP)
      real(RP),parameter :: twort2  = rt2       ! 2/sqrt(2) = sqrt(2)
      real(RP),parameter :: thrrt2  = 3._RP/rt2 ! 3/sqrt(2)
      real(RP) :: smax_inv                ! inverse
      real(RP) :: am,scrit_am,aa,tc,st,bb,ac
      real(RP) :: m0t,m2t,m3t,dgt,sgt,dm2
      real(RP) :: d_crit                  ! critical diameter
      real(RP) :: tmp1, tmp2, tmp3        !
      real(RP) :: ccn_frc,cca_frc,ccv_frc ! activated number,area,volume
      integer  :: is0, ik, ic

      aerosol_activ(:,:,:,:) = 0._RP

      if (super<=0._RP) return

      smax_inv = 1._RP / super

    !--- surface tension of water
      tc = temp_k - stdtemp
      if (tc >= 0._RP ) then
        st  = 75.94_RP-0.1365_RP*tc-0.3827e-3_RP*tc**2._RP !Gittens, JCIS, 69, Table 4.
      else !t[deg C]<0.
        st  = 75.93_RP                +0.115_RP*tc        &  !Eq.5-12, pp.130, PK97
            + 6.818e-2_RP*tc**2._RP+6.511e-3_RP*tc**3._RP &
            + 2.933e-4_RP*tc**4._RP+6.283e-6_RP*tc**5._RP &
            + 5.285e-8_RP*tc**6._RP
      endif
      st      = st * 1.E-3_RP                    ![J/m2]

    !-- Kelvin effect
    !          [J m-2]  [kg mol-1]       [m3 kg-1] [mol K J-1] [K-1]
      aa  = 2._RP * st * mwwat * 1.E-3_RP / (dnwat * rgas * temp_k ) ![m] Eq.5 in AR98

      do ic = 1, n_ctg
      do ik = 1, n_kap(ic)
      do is0 = 1, n_siz(ic)
        m0t = aerosol_procs(ia_m0,is0,ik,ic)
        m2t = aerosol_procs(ia_m2,is0,ik,ic)
        m3t = aerosol_procs(ia_m3,is0,ik,ic)
        call diag_ds(m0t,m2t,m3t,dgt,sgt,dm2)
        if (dgt <= 0._RP) dgt = d_ct(is0,ic)
        if (sgt <= 0._RP) sgt = 1.3_RP
        am  = dgt * 0.5_RP  !geometric dry mean radius [m]
        bb  = c_kappa
        if (bb > 0._RP .AND. am > 0._RP ) then
          scrit_am = 2._RP/sqrt(bb)*(aa/(3._RP*am))**1.5_RP !AR00 Eq.9
        else
          scrit_am = 0._RP
        endif
        ac     = am * (scrit_am * smax_inv)**two3      !AR00 Eq.12
        d_crit = ac * 2._RP
        tmp1   = log(d_crit) - log(dgt)
        tmp2   = 1._RP/(rt2*log(sgt))
        tmp3   = log(sgt)
        ccn_frc= 0.5_RP*(1._RP-erf(tmp1*tmp2))
        cca_frc= 0.5_RP*(1._RP-erf(tmp1*tmp2-twort2*tmp3))
        ccv_frc= 0.5_RP*(1._RP-erf(tmp1*tmp2-thrrt2*tmp3))
        ccn_frc=min(max(ccn_frc,0.0_RP),1.0_RP)
        cca_frc=min(max(cca_frc,0.0_RP),1.0_RP)
        ccv_frc=min(max(ccv_frc,0.0_RP),1.0_RP)
        aerosol_activ(ia_m0,is0,ik,ic) = ccn_frc * aerosol_procs(ia_m0,is0,ik,ic)
        aerosol_activ(ia_m2,is0,ik,ic) = cca_frc * aerosol_procs(ia_m2,is0,ik,ic)
        aerosol_activ(ia_m3,is0,ik,ic) = ccv_frc * aerosol_procs(ia_m3,is0,ik,ic)
        aerosol_activ(ia_ms,is0,ik,ic) = ccv_frc * aerosol_procs(ia_ms,is0,ik,ic)
      enddo !is(1:n_siz(ic))
      enddo !ik(1:n_kap(ic))
      enddo !ic(1:n_ctg)

      return
  end subroutine aerosol_activation
  !-----------------------------------------------------------------------------
  ! subroutine 5. diag_ds
  !----------------------------------------------------------------------------------------
  subroutine diag_ds(m0,m2,m3,  & !i
                     dg,sg,dm2)   !o
    implicit none
    real(RP)            :: m0,m2,m3,dg,sg,m3_bar,m2_bar
    real(RP)            :: m2_new,m2_old,dm2
    real(RP), parameter :: sgmax=2.5_RP
    real(RP), parameter :: rk1=2._RP
    real(RP), parameter :: rk2=3._RP
    real(RP), parameter :: ratio  =rk1/rk2
    real(RP), parameter :: rk1_hat=1._RP/(ratio*(rk2-rk1))
    real(RP), parameter :: rk2_hat=ratio/(rk1-rk2)
    real(DP), parameter :: tiny=1.E-50_DP

    dm2=0._RP

    if (m0 <= tiny .OR. m2 <= tiny .OR. m3 <= tiny) then
      m0=0._RP
      m2=0._RP
      m3=0._RP
      dg=-1._RP
      sg=-1._RP
      return
    endif

    m2_old = m2
    m3_bar = m3/m0
    m2_bar = m2/m0
    dg     = m2_bar**rk1_hat*m3_bar**rk2_hat

    if (m2_bar/m3_bar**ratio < 1._RP) then !stdev > 1.
      sg     = exp(sqrt(2._RP/(rk1*(rk1-rk2))  &
             * log(m2_bar/m3_bar**ratio) ))
    endif

    if (sg > sgmax) then
  !    print *,'sg=',sg
      sg = sgmax
      call diag_d2(m0,m3,sg, & !i
                   m2,dg     ) !o
  !   print *,'warning! modified as sg exceeded sgmax (diag_ds)'
    endif

    if (m2_bar/m3_bar**ratio >= 1._RP) then !stdev < 1.
      sg = 1._RP
      call diag_d2(m0,m3,sg, & !i
                   m2,dg     ) !o
  !   print *,'warning! modified as sg < 1. (diag_ds)'
    endif

    m2_new = m2
    dm2    = m2_old - m2_new !m2_pres - m2_diag

    return
  end subroutine diag_ds
  !----------------------------------------------------------------------------------------
  ! subroutine 6. diag_d2
  !----------------------------------------------------------------------------------------
  subroutine diag_d2(m0,m3,sg, & !i
                       m2,dg     ) !o
    implicit none
    real(RP)            :: dg,sg,m0,m2,m3,aaa
    real(RP), parameter :: one3=1._RP/3._RP
    aaa = m0               * exp( 4.5_RP * (log(sg)**2._RP) )
    dg  =(m3/aaa)**one3
    m2  = m0 * dg ** 2._RP * exp( 2.0_RP * (log(sg)**2._RP) )

    return
  end subroutine diag_d2
  !----------------------------------------------------------------------------------------
  ! subroutine 7. aero_intra            :: intra_sectional Brownian coagulation
  !----------------------------------------------------------------------------------------
  subroutine aero_intra(m0, m3,         & !input
                          dm0,dm3,dm6,    & !output
                          dg, sg, lambda, & !input
                          rknc,rkfm,dt)     !input

  !---intra-modal coagulation for spheres
    implicit none

  !=== Input
    real(RP), intent(in)   :: m0         ![/m3]
    real(RP), intent(in)   :: m3         ![m2/m3]
  ! real(RP), intent(in)   :: m6         ![m3/m3]
    real(RP), intent(out)  :: dm0        ![/m3]
    real(RP), intent(out)  :: dm3        ![m2/m3]
    real(RP), intent(out)  :: dm6        ![m3/m3]
    real(RP), intent(in)   :: dg         ![m]
    real(RP), intent(in)   :: sg         ![-]
    real(RP), intent(in)   :: lambda     ! mean free path of air molecules [cm]
    real(RP), intent(in)   :: rknc, rkfm ! constants for both regimes
    real(DP), intent(in)   :: dt         ! dt

  !=== Intermidiate quantities
    real(RP) :: mm2,mm1,m1,m4,m2,mm1p5,mm0p5,m0p5,m3p5,m1p5,m5,m2p5
    real(RP) :: dm0dt,dm6dt,dm0dt_nc,dm6dt_nc,dm0dt_fm,dm6dt_fm,bbb

  !=== Moment formulation for M0, M3 & M6 to extract Dg and Sg
    real(RP), parameter :: rk1=3._RP
    real(RP), parameter :: rk2=6._RP
    real(RP), parameter :: ratio=rk1/rk2
    real(RP), parameter :: rk1_hat=1._RP/(ratio*(rk2-rk1))
    real(RP), parameter :: rk2_hat=ratio/(rk1-rk2)

  !--- initialization
    dm0   = 0._RP
    dm3   = 0._RP
    dm6   = 0._RP
    dm0dt = 0._RP
    dm6dt = 0._RP

  !--- near continuum regime---------------------------------------------------------
  ! M(k)  = M(0) * dmean **   k  * exp( k**2.*0.5 *(dlog(stdev)**2.))
    mm2  =m0*dg**(-2._RP) *exp(2.0_RP *(log(sg)**2._RP)) !M-2 (dM0/dt)
    mm1  =m0*dg**(-1._RP) *exp(0.5_RP *(log(sg)**2._RP)) !M-1 (dM0/dt)
    m1   =m0*dg**  1._RP  *exp(0.5_RP *(log(sg)**2._RP)) !M1  (dM0/dt dM6/dt)
    m2   =m0*dg**  2._RP  *exp(2.0_RP *(log(sg)**2._RP)) !M2
    m4   =m0*dg**  4._RP  *exp(8.0_RP *(log(sg)**2._RP)) !M4  (dM6/dt)
  ! m6   =m0*dg**  6._RP  *exp(18._RP *(log(sg)**2._RP)) !M6

    dm0dt_nc=-1._RP*rknc*(m0*m0+m1*mm1+2.492e-2_RP*lambda*(m0*mm1+m1*mm2))
    dm6dt_nc= 2._RP*rknc*(m3*m3+m4*m2 +2.492e-2_RP*lambda*(m3*m2 +m4*m1 ))

  !--- free molecule regime----------------------------------------------------------
  ! M(k)  = M(0 * dmean **   k  * exp( k**2.*0.5 *(dlog(stdev)**2.) )
    mm1p5=m0*dg**(-1.5_RP)*exp(1.125_RP*(log(sg)**2._RP))!M-1.5(dM0/dt)
    mm0p5=m0*dg**(-0.5_RP)*exp(0.125_RP*(log(sg)**2._RP))!M-0.5(dM0/dt)
    m0p5 =m0*dg**  0.5_RP *exp(0.125_RP*(log(sg)**2._RP))!M0.5 (dM0/dt)
    m3p5 =m0*dg**( 3.5_RP)*exp(6.125_RP*(log(sg)**2._RP))!M3.5 (dM6/dt)
    m1p5 =m0*dg**( 1.5_RP)*exp(1.125_RP*(log(sg)**2._RP))!M1.5 (dM6/dt)
    m5   =m0*dg**  5._RP  *exp(12.50_RP*(log(sg)**2._RP))!M5   (dM6/dt)
    m2p5 =m0*dg**  2.5_RP *exp(3.125_RP*(log(sg)**2._RP))!M2.5 (dM6/dt)

    bbb     =1._RP+1.2_RP *exp(-2._RP*sg)-0.646_RP*exp(-0.35_RP*sg**2._RP) ! Eq.9 in Park et al. (1999) JAS

    dm0dt_fm= -1._RP*bbb*rkfm*(m0*m0p5+m2*mm1p5+2._RP*m1*mm0p5)
    dm6dt_fm=  2._RP*bbb*rkfm*(m3p5*m3+m1p5*m5 +2._RP*m2p5*m4 )

  !--- harmonic mean approach
    if (dm0dt_fm+dm0dt_nc /= 0._RP) dm0dt = dm0dt_fm*dm0dt_nc/(dm0dt_fm+dm0dt_nc)
    if (dm6dt_fm+dm6dt_nc /= 0._RP) dm6dt = dm6dt_fm*dm6dt_nc/(dm6dt_fm+dm6dt_nc)

    dm0    = dm0dt * dt
    dm3    = 0._RP
    dm6    = dm6dt * dt

    return
  end subroutine aero_intra
  !----------------------------------------------------------------------------------------
  ! subroutine 8. aero_inter            :: inter_sectional Brownian coagulation
  !----------------------------------------------------------------------------------------
  subroutine aero_inter(m0i ,m3i ,m6i ,      & !input unchanged
                        m0j ,m3j ,m6j ,      & !input unchanged
                        dm0i,dm3i,dm6i,dmsi, & !output
                        dm0j,dm3j,dm6j,dmsj, & !output
                        dm0k,dm3k,dm6k,dmsk, & !output
                        dgi ,sgi, rhoi,      & !input unchanged
                        dgj ,sgj, rhoj,      & !input unchanged
                        rknc,  rkfm, lambda, & !input unchanged
                        dtrest               ) !input unchanged

  !---inter-modal coagulation for spheres
    implicit none

  !=== input/output variables
    real(RP), intent(in)  :: m0i,m0j                 ! [/m3]
    real(RP), intent(in)  :: m3i,m3j                 ! [m3/m3]
    real(RP), intent(in)  :: m6i,m6j                 ! [m6/m3]
    real(RP), intent(out) :: dm0i,dm0j,dm0k          ! [/m3/s]
    real(RP), intent(out) :: dm3i,dm3j,dm3k          ! [m3/m3/s]
    real(RP), intent(out) :: dm6i,dm6j,dm6k          ! [m6/m3/s]
    real(RP), intent(out) :: dmsi,dmsj,dmsk          ! [ug/m3/s]
    real(RP), intent(in)  :: dgi,sgi,rhoi            ! [m][-][g/cm3]
    real(RP), intent(in)  :: dgj,sgj,rhoj            ! [m][-][g/cm3]
    real(RP), intent(in)  :: rknc, rkfm              ! constants for both regimes
    real(RP), intent(in)  :: lambda     ! mean free path of air molecules [cm]
    real(DP)              :: dtrest     ! time[s]

  !=== local variables
    real(RP) :: dm0dt_i,dm3dt_i,dm6dt_i
    real(RP) :: dm0dt_j,dm3dt_j,dm6dt_j
    real(RP) :: dm0dt_k,dm3dt_k,dm6dt_k
    real(RP) :: mm2i,mm1i,m1i,m2i,m4i,m5i,m7i,m8i
    real(RP) :: mm2j,mm1j,m1j,m2j,m4j,m5j,m7j,m8j
    real(RP) :: mm1p5i,mm0p5i,m0p5i,m1p5i,m2p5i,m3p5i,m4p5i,m5p5i,m6p5i
    real(RP) :: mm1p5j,mm0p5j,m0p5j,m1p5j,m2p5j,m3p5j,m4p5j,m5p5j,m6p5j
    real(RP) :: dm6dt_nc_i,dm3dt_nc_i,dm0dt_nc_i
    real(RP) :: dm6dt_nc_j,dm3dt_nc_j,dm0dt_nc_j
    real(RP) :: dm6dt_nc_k,dm3dt_nc_k,dm0dt_nc_k
    real(RP) :: dm6dt_fm_i,dm3dt_fm_i,dm0dt_fm_i
    real(RP) :: dm6dt_fm_j,dm3dt_fm_j,dm0dt_fm_j
    real(RP) :: dm6dt_fm_k,dm3dt_fm_k,dm0dt_fm_k
    real(RP) :: bbb,gamma1,gamma2,alpha,beta

  !=== initialization
    dm0dt_i = 0._RP
    dm3dt_i = 0._RP
    dm6dt_i = 0._RP
    dm0dt_j = 0._RP
    dm3dt_j = 0._RP
    dm6dt_j = 0._RP
    dm0dt_k = 0._RP
    dm3dt_k = 0._RP
    dm6dt_k = 0._RP
    dm0i    = 0._RP
    dm3i    = 0._RP
    dm6i    = 0._RP
    dmsi    = 0._RP
    dm0j    = 0._RP
    dm3j    = 0._RP
    dm6j    = 0._RP
    dmsj    = 0._RP
    dm0k    = 0._RP
    dm3k    = 0._RP
    dm6k    = 0._RP
    dmsk    = 0._RP

  !---near continuum regime----------------------------------------------------------------------
  ! (mode i)
    mm2i=m0i*dgi**(-2._RP)*exp( 2.0_RP*(log(sgi)**2._RP))!M-2
    mm1i=m0i*dgi**(-1._RP)*exp( 0.5_RP*(log(sgi)**2._RP))!M-1
    m1i =m0i*dgi**  1._RP *exp( 0.5_RP*(log(sgi)**2._RP))!M1
    m2i =m0i*dgi**  2._RP *exp( 2.0_RP*(log(sgi)**2._RP))!M2
  ! m3i =m3i
    m4i =m0i*dgi**  4._RP *exp( 8.0_RP*(log(sgi)**2._RP))!M4
    m5i =m0i*dgi**  5._RP *exp(12.5_RP*(log(sgi)**2._RP))!M5
  ! m6i =m6i
    m7i =m0i*dgi**  7._RP *exp(24.5_RP*(log(sgi)**2._RP))!M7
  ! (mode j)
    mm2j=m0j*dgj**(-2._RP)*exp( 2.0_RP*(log(sgj)**2._RP))!M-2
    mm1j=m0j*dgj**(-1._RP)*exp( 0.5_RP*(log(sgj)**2._RP))!M-1
    m1j =m0j*dgj**  1._RP *exp( 0.5_RP*(log(sgj)**2._RP))!M1
    m2j =m0j*dgj**  2._RP *exp( 2.0_RP*(log(sgj)**2._RP))!M2
  ! m3j =m3j
    m4j =m0j*dgj**  4._RP *exp( 8.0_RP*(log(sgj)**2._RP))!M4
    m5j =m0j*dgj**  5._RP *exp(12.5_RP*(log(sgj)**2._RP))!M5
  ! m6j =m6j
    m7j =m0j*dgj**  7._RP *exp(24.5_RP*(log(sgj)**2._RP))!M7

    dm0dt_nc_i = -rknc*(2._RP*m0i*m0j+ mm1i*m1j             &
                                      + mm1j*m1i            &
                   +2.492e-2_RP*lambda*(m0i *mm1j+m1i*mm2j  &
                                      + m0j *mm1i+m1j*mm2i) )
    dm3dt_nc_i = -rknc*(2._RP*m3i*m0j+ m2i *m1j             &
                                      + mm1j*m4i            &
                   +2.492e-2_RP*lambda*(m3i *mm1j+m4i*mm2j  &
                                      + m0j *m2i +m1j*m1i ) )
    dm6dt_nc_i = -rknc*(2._RP*m6i*m0j+ m5i *m1j             &
                                      + mm1j*m7i            &
                   +2.492e-2_RP*lambda*(m6i *mm1j+m7i*mm2j  &
                                      + m0j *m5i +m1j*m4i ) )
    dm0dt_nc_j = -rknc*(2._RP*m0i*m0j+ mm1i*m1j             & !=dm0dt_nc_i
                                      + mm1j*m1i            &
                   +2.492e-2_RP*lambda*(m0i *mm1j+m1i*mm2j  &
                                      + m0j *mm1i+m1j*mm2i) )
    dm3dt_nc_j = -rknc*(2._RP*m0i*m3j+ mm1i*m4j             &
                                      + m2j *m1i            &
                   +2.492e-2_RP*lambda*(m0i *m2j +m1i*m1j   &
                                      + m3j *mm1i+m4j*mm2i) )
    dm6dt_nc_j = -rknc*(2._RP*m0i*m6j+ mm1i*m7j             &
                                      + m5j *m1i            &
                   +2.492e-2_RP*lambda*(m0i *m5j +m1i*m4j   &
                                      + m6j *mm1i+m7j*mm2i) )
    dm0dt_nc_k = -dm0dt_nc_i
    dm3dt_nc_k = -dm3dt_nc_i-dm3dt_nc_j
    dm6dt_nc_k = -dm6dt_nc_i-dm6dt_nc_j                     &
            +2._RP*rknc*(2._RP*m3i*m3j+ m2i *m4j            &
                                      + m2j *m4i            &
                   +2.492e-2_RP*lambda*(m3i *m2j +m4i*m1j   &
                                      + m3j *m2i +m4j*m1i)  )

  !---free molecule regime-----------------------------------------------------------------------
  ! (mode i)
    mm1p5i=m0i*dgi**(-1.5_RP)*exp( 1.125_RP*(log(sgi)**2._RP) ) !M-1.5
    mm0p5i=m0i*dgi**(-0.5_RP)*exp( 0.125_RP*(log(sgi)**2._RP) ) !M-0.5
    m0p5i =m0i*dgi**  0.5_RP *exp( 0.125_RP*(log(sgi)**2._RP) ) !M0.5
    m1p5i =m0i*dgi**  1.5_RP *exp( 1.125_RP*(log(sgi)**2._RP) ) !M1.5
    m2p5i =m0i*dgi**  2.5_RP *exp( 3.125_RP*(log(sgi)**2._RP) ) !M2.5
    m3p5i =m0i*dgi**  3.5_RP *exp( 6.125_RP*(log(sgi)**2._RP) ) !M3.5
    m4p5i =m0i*dgi**  4.5_RP *exp(10.125_RP*(log(sgi)**2._RP) ) !M4.5
    m5p5i =m0i*dgi**  5.5_RP *exp(15.125_RP*(log(sgi)**2._RP) ) !M5.5
    m6p5i =m0i*dgi**  6.5_RP *exp(21.125_RP*(log(sgi)**2._RP) ) !M6.5
    m8i   =m0i*dgi**   8._RP *exp(32.000_RP*(log(sgi)**2._RP) ) !M8

  ! (mode j)
    mm1p5j=m0j*dgj**(-1.5_RP)*exp( 1.125_RP*(log(sgj)**2._RP) ) !M-1.5
    mm0p5j=m0j*dgj**(-0.5_RP)*exp( 0.125_RP*(log(sgj)**2._RP) ) !M-0.5
    m0p5j =m0j*dgj**  0.5_RP *exp( 0.125_RP*(log(sgj)**2._RP) ) !M0.5
    m1p5j =m0j*dgj**  1.5_RP *exp( 1.125_RP*(log(sgj)**2._RP) ) !M1.5
    m2p5j =m0j*dgj**  2.5_RP *exp( 3.125_RP*(log(sgj)**2._RP) ) !M2.5
    m3p5j =m0j*dgj**  3.5_RP *exp( 6.125_RP*(log(sgj)**2._RP) ) !M3.5
    m4p5j =m0j*dgj**  4.5_RP *exp(10.125_RP*(log(sgj)**2._RP) ) !M4.5
    m5p5j =m0j*dgj**  5.5_RP *exp(15.125_RP*(log(sgj)**2._RP) ) !M5.5
    m6p5j =m0j*dgj**  6.5_RP *exp(21.125_RP*(log(sgj)**2._RP) ) !M6.5
    m8j   =m0j*dgj**   8._RP *exp(32.000_RP*(log(sgj)**2._RP) ) !M8

  !---approximation function---------------------------------------------------------------------
    alpha   = dgj/dgi
    beta    =(1._RP-(sqrt(1._RP+alpha**3._RP)/(1._RP+sqrt(alpha**3._RP))))/(1._RP-1._RP/sqrt(2._RP))
    gamma1  =(sgi        +alpha*sgj       )/(1._RP+alpha)
    gamma2  =(sgi**2._RP +alpha*sgj**2._RP)/(1._RP+alpha)
    bbb     = 1._RP+1.2_RP*beta*dexp(real(-2._RP*gamma1,kind=DP))-0.646_RP*beta*dexp(real(-0.35_RP*gamma2,kind=DP))! Kajino (2011) JAS

    dm0dt_fm_i=-bbb*rkfm*(                                 &
                m0i *m0p5j +m2i *mm1p5j +2._RP*m1i *mm0p5j &
              + m0j *m0p5i +m2j *mm1p5i +2._RP*m1j *mm0p5i )
    dm3dt_fm_i=-bbb*rkfm*(                                 &
                m3i *m0p5j +m5i *mm1p5j +2._RP*m4i *mm0p5j &
              + m0j *m3p5i +m2j *m1p5i  +2._RP*m1j *m2p5i  )
    dm6dt_fm_i=-bbb*rkfm*(                                 &
                m6i *m0p5j +m8i *mm1p5j +2._RP*m7i *mm0p5j &
              + m0j *m6p5i +m2j *m4p5i  +2._RP*m1j *m5p5i  )
    dm0dt_fm_j= dm0dt_fm_i
    dm3dt_fm_j=-bbb*rkfm*(                                 &
                m0i *m3p5j +m2i *m1p5j  +2._RP*m1i *m2p5j  &
              + m3j *m0p5i +m5j *mm1p5i +2._RP*m4j *mm0p5i )
    dm6dt_fm_j=-bbb*rkfm*(                                 &
                m0i *m6p5j +m2i *m4p5j  +2._RP*m1i *m5p5j  &
              + m6j *m0p5i +m8j *mm1p5i +2._RP*m7j *mm0p5i )
    dm0dt_fm_k=-dm0dt_fm_i
    dm3dt_fm_k=-dm3dt_fm_i-dm3dt_fm_j
    dm6dt_fm_k=-dm6dt_fm_i-dm6dt_fm_j                      &
        +2._RP * bbb*rkfm*(                                &
                m3i *m3p5j +m5i *m1p5j  +2._RP*m4i *m2p5j  &
              + m3j *m3p5i +m5j *m1p5i  +2._RP*m4j *m2p5i  )

  !---harmonic mean approach---------------------------------------------------------------------
    if (dm0dt_fm_i+dm0dt_nc_i/=0._RP) &
    dm0dt_i  = dm0dt_fm_i*dm0dt_nc_i/(dm0dt_fm_i+dm0dt_nc_i)
    if (dm3dt_fm_i+dm3dt_nc_i/=0._RP) &
    dm3dt_i  = dm3dt_fm_i*dm3dt_nc_i/(dm3dt_fm_i+dm3dt_nc_i)
    if (dm6dt_fm_i+dm6dt_nc_i/=0._RP) &
    dm6dt_i  = dm6dt_fm_i*dm6dt_nc_i/(dm6dt_fm_i+dm6dt_nc_i)
    if (dm0dt_fm_j+dm0dt_nc_j/=0._RP) &
    dm0dt_j  = dm0dt_fm_j*dm0dt_nc_j/(dm0dt_fm_j+dm0dt_nc_j)
    if (dm3dt_fm_j+dm3dt_nc_j/=0._RP) &
    dm3dt_j  = dm3dt_fm_j*dm3dt_nc_j/(dm3dt_fm_j+dm3dt_nc_j)
    if (dm6dt_fm_j+dm6dt_nc_j/=0._RP) &
    dm6dt_j  = dm6dt_fm_j*dm6dt_nc_j/(dm6dt_fm_j+dm6dt_nc_j)
    if (dm0dt_fm_k+dm0dt_nc_k/=0._RP) &
    dm0dt_k  = dm0dt_fm_k*dm0dt_nc_k/(dm0dt_fm_k+dm0dt_nc_k)
  ! if (dm3dt_fm_k+dm3dt_nc_k/=0._RP) &
  ! dm3dt_k  = dm3dt_fm_k*dm3dt_nc_k/(dm3dt_fm_k+dm3dt_nc_k)
    if (dm6dt_fm_k+dm6dt_nc_k/=0._RP) &
    dm6dt_k  = dm6dt_fm_k*dm6dt_nc_k/(dm6dt_fm_k+dm6dt_nc_k)

    dm0i     = dm0dt_i * dtrest
    dm3i     = dm3dt_i * dtrest
    dm6i     = dm6dt_i * dtrest
    dmsi     = dm3i    * rhoi *pi6*1.E12_RP !3rd to mass !m3/m3=>ug/m3
    dm0j     = dm0dt_j * dtrest
    dm3j     = dm3dt_j * dtrest
    dm6j     = dm6dt_j * dtrest
    dmsj     = dm3j    * rhoj *pi6*1.E12_RP !3rd to mass !m3/m3=>ug/m3
    dm0k     = dm0dt_k * dtrest
  ! dm3k     = dm3dt_k * dtrest
    dm6k     = dm6dt_k * dtrest
    dm3k     = -dm3i   -dm3j
    dmsk     = -dmsi   -dmsj

    return
  end subroutine aero_inter
  !----------------------------------------------------------------------------------------
  ! subroutine 9. diag_ds6
  !----------------------------------------------------------------------------------------
  subroutine diag_ds6(m0,m3,m6,m2,dg,sg,dm2)
    implicit none
    real(RP)            :: m0,m2,m3,m6,dg,sg,m3_bar,m6_bar
    real(RP)            :: m2_new,m2_old,dm2
    real(RP), parameter :: sgmax=2.5_RP
    real(RP), parameter :: rk1=3._RP
    real(RP), parameter :: rk2=6._RP
    real(RP), parameter :: ratio  =rk1/rk2
    real(RP), parameter :: rk1_hat=1._RP/(ratio*(rk2-rk1))
    real(RP), parameter :: rk2_hat=ratio/(rk1-rk2)
    real(DP), parameter :: tiny=1.E-50_DP

    dm2=0._RP
    if (m0 <= tiny .OR. m3 <= tiny .OR. m6 <= tiny) then
      m0=0._RP
      m2=0._RP
      m3=0._RP
      dg=-1._RP
      sg=-1._RP
      return
    endif

    m2_old = m2
    m3_bar = m3/m0
    m6_bar = m6/m0
    dg     = m3_bar**rk1_hat*m6_bar**rk2_hat

    if (m3_bar/m6_bar**ratio < 1._RP) then !stdev > 1.
      sg     = exp(dsqrt(real(2._RP/(rk1*(rk1-rk2))  &
                  * log(m3_bar/m6_bar**ratio), kind=DP )))
      m2     = m0*dg**2._RP*exp(2._RP*(log(sg)**2._RP))
      m2_old = m2
    endif

    if (sg > sgmax) then
 !     print *,'sg=',sg
      sg = sgmax
      call diag_d2(m0,m3,sg, & !input
                   m2,dg     ) !output
  !   print *,'warning! modified as sg exceeded sgmax (diag_ds6)'
    endif

    if (m3_bar/m6_bar**ratio >= 1._RP) then !stdev < 1.
      sg = 1._RP
      call diag_d2(m0,m3,sg, & !input
                   m2,dg     ) !output
  !   print *,'warning! modified as sg < 1.(diag_ds6)'
    endif

    m2_new = m2
    dm2    = m2_old - m2_new !m2_pres - m2_diag


    return
  end subroutine diag_ds6
  !-----------------------------------------------------------------------------
  subroutine trans_ccn(aerosol_procs, aerosol_activ, t_ccn, t_cn,  &
                    n_ctg, n_kap_max, n_siz_max, n_atr,         &
                    ic_out, ia_m0, ia_m2, ia_m3, ik_out, n_siz, &
!                    rnum_out, nbins_out, dlog10d_out)
                    rnum_out, nbins_out)

    implicit none
    !i/o variables
    integer, intent(in) :: n_ctg, n_kap_max, n_siz_max, n_atr
    integer, intent(in) :: ic_out, ik_out
    integer, intent(in) :: ia_m0, ia_m2, ia_m3
    integer, intent(in) :: n_siz(n_ctg)
    real(RP),intent(in) :: aerosol_procs(n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP),intent(in) :: aerosol_activ(n_atr,n_siz_max,n_kap_max,n_ctg)
    integer, intent(in) :: nbins_out
    real(RP),intent(out):: rnum_out(nbins_out)
!    real(RP),intent(out):: dlog10d_out
    real(RP),intent(out):: t_ccn !total ccn number concentration [#/m3]
    real(RP),intent(out):: t_cn  !total cn number concentration  [#/m3]
    real(RP) :: dlog10d_out
    !local variables
    real(RP) :: m0t,m2t,m3t,dgt,sgt,dm2
    real(RP) :: d_lw2,d_up2,dg2,sg2,fnum0,fnum1
    real(RP) :: dlogd
    integer  :: is0, is_out
    real(RP), parameter :: d_max = 1.E-5_RP
    real(RP), parameter :: d_min = 1.E-9_RP
    real(RP) :: d_lw(nbins_out), d_up(nbins_out)

    rnum_out(:) = 0._RP
!
!    dlogd = (log(d_max) - log(d_min))/float(nbins_out)
!    dlog10d_out = (log10(d_max)-log10(d_min))/float(nbins_out)
!
!    do is_out = 1, nbins_out  !size bin
!      d_lw(is_out) = exp(log(d_min)+dlogd* float(is_out-1) )
!      d_up(is_out) = exp(log(d_min)+dlogd* float(is_out)   )
!    enddo !is_out
!
    t_ccn = 0._RP
    t_cn  = 0._RP
    do is0 = 1, n_siz(ic_out)
!      if (aerosol_procs(ia_m0,is0,ik_out,ic_out) > 0._RP) then
!        m0t = aerosol_procs(ia_m0,is0,ik_out,ic_out)
!        m2t = aerosol_procs(ia_m2,is0,ik_out,ic_out)
!        m3t = aerosol_procs(ia_m3,is0,ik_out,ic_out)
!        call diag_ds(m0t,m2t,m3t,dgt,sgt,dm2)
!        if (dgt <= 0._RP) dgt = d_ct(is0,ic_out)
!        if (sgt <= 0._RP) sgt = 1.3_RP
!
!        do is_out = 1, nbins_out
!          sgt   = max(sgt,1.0000001_RP) !to avoid floating divide by zero
!          d_lw2 = log(d_lw(is_out))
!          d_up2 = log(d_up(is_out))
!          dg2   = log(dgt)
!          sg2   = log(sgt)
!          fnum0 = m0t*0.5_RP*(1._RP+erf((d_lw2-dg2)/(sqrt(2.0_RP)*sg2)))
!          fnum1 = m0t*0.5_RP*(1._RP+erf((d_up2-dg2)/(sqrt(2.0_RP)*sg2)))
!          rnum_out(is_out) = rnum_out(is_out) + fnum1 - fnum0
!        enddo
!
!      endif !number>0

      t_ccn = t_ccn + aerosol_activ(ia_m0,is0,ik_out,ic_out)
      t_cn  = t_cn  + aerosol_procs(ia_m0,is0,ik_out,ic_out)

    enddo

  return
  end subroutine trans_ccn
  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AE_kajino13_EffectiveRadius( &
       Re,   &
       QTRC, &
       RH    )
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_aerosol, only: &
       N_AE
    implicit none
    real(RP), intent(out) :: Re  (KA,IA,JA,N_AE) ! effective radius
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)   ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: RH  (KA,IA,JA)      ! relative humidity         (0-1)
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = 0.0_RP

!    Re(:,:,:,I_ae_seasalt) = 2.E-4_RP
!    Re(:,:,:,I_ae_dust   ) = 4.E-6_RP
!    Re(:,:,:,I_ae_bc     ) = 4.E-8_RP
!    Re(:,:,:,I_ae_oc     ) = RH(:,:,:)
!    Re(:,:,:,I_ae_sulfate) = RH(:,:,:)

    return
  end subroutine ATMOS_PHY_AE_kajino13_EffectiveRadius

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_AE_kajino13_mkinit( &
       QTRC, CCN, &
       DENS, RHOT, &
       m0_init, dg_init, sg_init, &
       d_min_inp, d_max_inp, &
       k_min_inp, k_max_inp, &
       n_kap_inp )
    use scale_const, only: &
       PI => CONST_PI, &
       CONST_CVdry, &
       CONST_CPdry, &
       CONST_Rdry, &
       CONST_PRE00
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_atmos_saturation, only: &
       SATURATION_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq
    implicit none

    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(out)   :: CCN (KA,IA,JA)
    real(RP), intent(in)    :: DENS(KA,IA,JA)
    real(RP), intent(in)    :: RHOT(KA,IA,JA)
    real(RP), intent(in)    :: m0_init             ! initial total num. conc. of modes (Atk,Acm,Cor) [#/m3]
    real(RP), intent(in)    :: dg_init             ! initial number equivalen diameters of modes     [m]
    real(RP), intent(in)    :: sg_init             ! initial standard deviation                      [-]
    real(RP), intent(in)    :: d_min_inp(3)
    real(RP), intent(in)    :: d_max_inp(3)
    real(RP), intent(in)    :: k_min_inp(3)
    real(RP), intent(in)    :: k_max_inp(3)
    integer,  intent(in)    :: n_kap_inp(3)

    integer,  parameter :: ia_m0  = 1             ! 1. number conc        [#/m3]
    integer,  parameter :: ia_m2  = 2             ! 2. 2nd mom conc       [m2/m3]
    integer,  parameter :: ia_m3  = 3             ! 3. 3rd mom conc       [m3/m3]
    integer,  parameter :: ia_ms  = 4             ! 4. mass conc          [ug/m3]
    integer,  parameter :: ia_kp  = 5
    integer,  parameter :: ik_out = 1

    real(RP)            :: c_kappa     = 0.3_RP     ! hygroscopicity of condensable mass [-]
    real(RP), parameter :: cleannumber = 1.e-3_RP
    ! local variables
    real(RP), parameter :: rhod_ae    = 1.83_RP              ! particle density [g/cm3] sulfate assumed
    real(RP), parameter :: conv_vl_ms = rhod_ae/1.e-12_RP     ! M3(volume)[m3/m3] to mass[m3/m3]

    integer  :: n_trans
    real(RP) :: m0t, dgt, sgt, m2t, m3t, mst
    real(RP), allocatable :: aerosol_procs(:,:,:,:) ! (n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP), allocatable :: aerosol_activ(:,:,:,:) ! (n_atr,n_siz_max,n_kap_max,n_ctg)
    real(RP), allocatable :: emis_procs   (:,:,:,:) ! (n_atr,n_siz_max,n_kap_max,n_ctg)
    !--- gas
    real(RP) :: conc_gas(GAS_CTG)    !concentration [ug/m3]
    integer  :: n_siz_max, n_kap_max, n_ctg
!   integer,  allocatable :: it_procs2trans(:,:,:,:) !procs to trans conversion
!   integer,  allocatable :: ia_trans2procs(:) !trans to procs conversion
!   integer,  allocatable :: is_trans2procs(:) !trans to procs conversion
!   integer,  allocatable :: ik_trans2procs(:) !trans to procs conversion
!   integer,  allocatable :: ic_trans2procs(:)
    !--- bin  settings (lower bound, center, upper bound)
    real(RP), allocatable :: d_lw(:,:), d_ct(:,:), d_up(:,:)  !diameter [m]
    real(RP), allocatable :: k_lw(:,:), k_ct(:,:), k_up(:,:)  !kappa    [-]
    real(RP)  :: dlogd, dk                                    !delta log(D), delta K
    real(RP), allocatable :: d_min(:)              !lower bound of 1st size bin   (n_ctg)
    real(RP), allocatable :: d_max(:)              !upper bound of last size bin  (n_ctg)
    integer,  allocatable :: n_kap(:)              !number of kappa bins          (n_ctg)
    real(RP), allocatable :: k_min(:)              !lower bound of 1st kappa bin  (n_ctg)
    real(RP), allocatable :: k_max(:)

    real(RP) :: pott, qdry, pres
    real(RP) :: temp, cpa, cva, qsat_tmp, ssliq, Rmoist
    real(RP) :: pi6

    integer  :: ia0, ik, is0, ic, k, i, j, it, iq
    !---------------------------------------------------------------------------


    pi6   = pi / 6._RP
    n_ctg = AE_CTG

    allocate( d_min(n_ctg) )
    allocate( d_max(n_ctg) )
    allocate( n_kap(n_ctg) )
    allocate( k_min(n_ctg) )
    allocate( k_max(n_ctg) )


    d_min(1:n_ctg) = d_min_inp(1:n_ctg)  ! lower bound of 1st size bin
    d_max(1:n_ctg) = d_max_inp(1:n_ctg)  ! upper bound of last size bin
    n_kap(1:n_ctg) = n_kap_inp(1:n_ctg)  ! number of kappa bins
    k_min(1:n_ctg) = k_min_inp(1:n_ctg)  ! lower bound of 1st kappa bin
    k_max(1:n_ctg) = k_max_inp(1:n_ctg)  ! upper bound of last kappa bin

    !--- diagnose parameters (n_trans, n_siz_max, n_kap_max)
    n_trans   = 0
    n_siz_max = 0
    n_kap_max = 0
    do ic = 1, n_ctg
      n_trans   = n_trans + NSIZ(ic) * NKAP(ic) * N_ATR
      n_siz_max = max(n_siz_max, NSIZ(ic))
      n_kap_max = max(n_kap_max, NKAP(ic))
    enddo

    allocate( aerosol_procs (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( aerosol_activ (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    allocate( emis_procs    (N_ATR,n_siz_max,n_kap_max,n_ctg)  )
    aerosol_procs(:,:,:,:) = 0._RP
    aerosol_activ(:,:,:,:) = 0._RP
    emis_procs   (:,:,:,:) = 0._RP

!   allocate( it_procs2trans(N_ATR,n_siz_max,n_kap_max,n_ctg)  )
!   allocate( ia_trans2procs(n_trans) )
!   allocate( is_trans2procs(n_trans) )
!   allocate( ik_trans2procs(n_trans) )
!   allocate( ic_trans2procs(n_trans) )

    !bin setting
    allocate(d_lw(n_siz_max,n_ctg))
    allocate(d_ct(n_siz_max,n_ctg))
    allocate(d_up(n_siz_max,n_ctg))
    allocate(k_lw(n_kap_max,n_ctg))
    allocate(k_ct(n_kap_max,n_ctg))
    allocate(k_up(n_kap_max,n_ctg))
    d_lw(:,:) = 0._RP
    d_ct(:,:) = 0._RP
    d_up(:,:) = 0._RP
    k_lw(:,:) = 0._RP
    k_ct(:,:) = 0._RP
    k_up(:,:) = 0._RP

    do ic = 1, n_ctg

      dlogd = (log(d_max(ic)) - log(d_min(ic)))/real(NSIZ(ic),kind=RP)

      do is0 = 1, NSIZ(ic)  !size bin
        d_lw(is0,ic) = exp(log(d_min(ic))+dlogd* real(is0-1,kind=RP)        )
        d_ct(is0,ic) = exp(log(d_min(ic))+dlogd*(real(is0  ,kind=RP)-0.5_RP))
        d_up(is0,ic) = exp(log(d_min(ic))+dlogd* real(is0  ,kind=RP)        )
      enddo !is (1:n_siz(ic))
      dk    = (k_max(ic) - k_min(ic))/real(n_kap(ic),kind=RP)
      do ik = 1, n_kap(ic)  !size bin
        k_lw(ik,ic) = k_min(ic) + dk  * real(ik-1,kind=RP)
        k_ct(ik,ic) = k_min(ic) + dk  *(real(ik  ,kind=RP)-0.5_RP)
        k_up(ik,ic) = k_min(ic) + dk  * real(ik  ,kind=RP)
      enddo !ik (1:n_kap(ic))

    enddo !ic (1:n_ctg)
!    ik  = 1       !only one kappa bin

    m0t = m0_init !total M0 [#/m3]
    dgt = dg_init ![m]
    sgt = sg_init ![-]

    if ( m0t <= cleannumber ) then
       m0t = cleannumber
       dgt = 0.1E-6_RP
       sgt = 1.3_RP
       if( IO_L ) write(IO_FID_LOG,*) '*** WARNING! Initial aerosol number is set as ', cleannumber, '[#/m3]'
    endif

    m2t = m0t*dgt**(2.d0) *dexp(2.0d0 *(dlog(real(sgt,kind=DP))**2.d0)) !total M2 [m2/m3]
    m3t = m0t*dgt**(3.d0) *dexp(4.5d0 *(dlog(real(sgt,kind=DP))**2.d0)) !total M3 [m3/m3]
    mst = m3t*pi6*conv_vl_ms                              !total Ms [ug/m3]

    do ic = 1, n_ctg
    !aerosol_procs initial condition
    do ik = 1, n_kap(ic)   !kappa bin
    do is0 = 1, NSIZ(ic)
       if (dgt >= d_lw(is0,ic) .and. dgt < d_up(is0,ic)) then
         aerosol_procs(ia_m0,is0,ik,ic) = aerosol_procs(ia_m0,is0,ik,ic) + m0t          ![#/m3]
         aerosol_procs(ia_m2,is0,ik,ic) = aerosol_procs(ia_m2,is0,ik,ic) + m2t          ![m2/m3]
         aerosol_procs(ia_m3,is0,ik,ic) = aerosol_procs(ia_m3,is0,ik,ic) + m3t          ![m3/m3]
         aerosol_procs(ia_ms,is0,ik,ic) = aerosol_procs(ia_ms,is0,ik,ic) + mst*1.E-9_RP ![kg/m3]
       elseif (dgt < d_lw(1,ic)) then
         aerosol_procs(ia_m0,1 ,ik,ic) = aerosol_procs(ia_m0,1 ,ik,ic) + m0t          ![#/m3]
         aerosol_procs(ia_m2,1 ,ik,ic) = aerosol_procs(ia_m2,1 ,ik,ic) + m2t          ![m2/m3]
         aerosol_procs(ia_m3,1 ,ik,ic) = aerosol_procs(ia_m3,1 ,ik,ic) + m3t          ![m3/m3]
         aerosol_procs(ia_ms,1 ,ik,ic) = aerosol_procs(ia_ms,1 ,ik,ic) + mst*1.E-9_RP ![kg/m3]
       elseif (dgt >= d_up(NSIZ(ic),ic)) then
         aerosol_procs(ia_m0,NSIZ(ic),ik,ic) = aerosol_procs(ia_m0,NSIZ(ic),ik,ic) + m0t          ![#/m3]
         aerosol_procs(ia_m2,NSIZ(ic),ik,ic) = aerosol_procs(ia_m2,NSIZ(ic),ik,ic) + m2t          ![m2/m3]
         aerosol_procs(ia_m3,NSIZ(ic),ik,ic) = aerosol_procs(ia_m3,NSIZ(ic),ik,ic) + m3t          ![m3/m3]
         aerosol_procs(ia_ms,NSIZ(ic),ik,ic) = aerosol_procs(ia_ms,NSIZ(ic),ik,ic) + mst*1.E-9_RP ![kg/m3]
       endif
    enddo
    enddo
    enddo

    conc_gas(:) = 0.0_RP
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       it = 0
       do ic  = 1, n_ctg    ! category
       do ik  = 1, NKAP(ic) ! kappa bin
       do is0 = 1, NSIZ(ic) ! size bin
       do ia0 = 1, N_ATR    ! attributes
          QTRC(k,i,j,QAES+it) = aerosol_procs(ia0,is0,ik,ic) / DENS(k,i,j) !#,m2,m3,kg/m3 -> #,m2,m3,kg/kg
          it = it + 1
       enddo !ia0 (1:N_ATR )
       enddo !is (1:n_siz(ic)  )
       enddo !ik (1:n_kap(ic)  )
       enddo !ic (1:n_ctg      )

       do ic  = 1, GAS_CTG ! GAS category
          QTRC(k,i,j,QAES+it) = conc_gas(ic) / DENS(k,i,j) * 1.E-9_RP !mixing ratio [kg/kg]
          it = it + 1
       enddo
    enddo
    enddo
    enddo

    CCN(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       !--- calculate super saturation of water
       if ( I_QV > 0 ) then
          CALC_QDRY(qdry, QTRC, TRACER_MASS, k, i, j, iq)
          CALC_R   (Rmoist, qdry, QTRC, k, i, j, iq, CONST_Rdry,  TRACER_R )
          CALC_CV  (cva,    qdry, QTRC, k, i, j, iq, CONST_CVdry, TRACER_CV)
          CALC_CP  (cpa,    qdry, QTRC, k, i, j, iq, CONST_CPdry, TRACER_CP)

          pres = CONST_PRE00 * ( RHOT(k,i,j) * Rmoist / CONST_PRE00 )**(cpa/cva)
          temp = pres / ( DENS(k,i,j) * Rmoist )

          call SATURATION_pres2qsat_liq( qsat_tmp, temp, pres )

          ssliq = QTRC(k,i,j,I_QV) / qsat_tmp - 1.0_RP
       else
          ssliq = -1.0_RP
       endif

       call aerosol_activation( c_kappa, ssliq, temp, ia_m0, ia_m2, ia_m3,       &
                                N_ATR, n_siz_max, n_kap_max, n_ctg, NSIZ, n_kap, &
                                d_ct, aerosol_procs, aerosol_activ               )

       do is0 = 1, NSIZ(ic_mix)
          CCN(k,i,j) = CCN(k,i,j) + aerosol_activ(ia_m0,is0,ik_out,ic_mix)
       enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_AE_kajino13_mkinit

end module scale_atmos_phy_ae_kajino13
