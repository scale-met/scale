!-------------------------------------------------------------------------------
!> Module numerical filter
!!
!! @par Description
!!          This module contains subroutines for numerical smoothings or filters
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_numfilter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: numfilter_setup
  public :: numfilter_rayleigh_damping
  public :: numfilter_hdiffusion
  public :: numfilter_vdiffusion
  public :: numfilter_divdamp
  public :: numfilter_divdamp_2d

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: NUMFILTER_DOrayleigh            = .false. ! use rayleigh damping?
  logical, public :: NUMFILTER_DOhorizontaldiff      = .false. ! use horizontal diffusion?
  logical, public :: NUMFILTER_DOhorizontaldiff_lap1 = .false. ! use horizontal 1st-order damping? (for upper layer)
  logical, public :: NUMFILTER_DOverticaldiff        = .false. ! use vertical diffusion?
  logical, public :: NUMFILTER_DOdivdamp             = .false. ! use 3D divergence damping?
  logical, public :: NUMFILTER_DOdivdamp_v           = .false. ! use 3D divergence damping for vertical velocity?
  logical, public :: NUMFILTER_DOdivdamp_2d          = .false. ! use 2D divergence damping?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: numfilter_rayleigh_damping_setup
  private :: numfilter_hdiffusion_setup
  private :: numfilter_vdiffusion_setup
  private :: numfilter_divdamp_setup
  private :: numfilter_divdamp_2d_setup
  private :: numfilter_smooth_1var
  private :: height_factor

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_RHOG     = 1 ! Density x G^1/2 x gamma^2
  integer, private, parameter :: I_RHOGVX   = 2 ! Density x G^1/2 x gamma^2 x Horizontal velocity (X-direction)
  integer, private, parameter :: I_RHOGVY   = 3 ! Density x G^1/2 x gamma^2 x Horizontal velocity (Y-direction)
  integer, private, parameter :: I_RHOGVZ   = 4 ! Density x G^1/2 x gamma^2 x Horizontal velocity (Z-direction)
  integer, private, parameter :: I_RHOGW    = 5 ! Density x G^1/2 x gamma^2 x Vertical   velocity
  integer, private, parameter :: I_RHOGE    = 6 ! Density x G^1/2 x gamma^2 x Internal Energy
  integer, private, parameter :: I_RHOGETOT = 7 ! Density x G^1/2 x gamma^2 x Total Energy

  real(RP), public,  allocatable :: rayleigh_coef  (:)             ! Rayleigh damping coefficient at cell center
  real(RP), private, allocatable :: rayleigh_coef_h(:)             ! Rayleigh damping coefficient at cell wall
  logical, private              :: rayleigh_damp_only_w = .false. ! damp only w?

  real(RP), public,  allocatable :: Kh_coef   (:,:,:)              ! horizontal diffusion coefficient at cell center
  real(RP), private, allocatable :: Kh_coef_pl(:,:,:)
  integer, private              :: lap_order_hdiff = 2            ! laplacian order
  real(RP), private              :: hdiff_fact_rho  = 1.E-2_RP
  real(RP), private              :: hdiff_fact_q    = 0.0_RP
  real(RP), private              :: Kh_coef_minlim  = 0.E+00_RP
  real(RP), private              :: Kh_coef_maxlim  = 1.E+30_RP

  logical, private              :: hdiff_nonlinear = .false.
  real(RP), private              :: ZD_hdiff_nl     = 25000.0_RP     ! hight for decay of nonlinear diffusion

  real(RP), public,  allocatable :: Kh_coef_lap1   (:,:,:)         ! Kh_coef but 1st order laplacian
  real(RP), private, allocatable :: Kh_coef_lap1_pl(:,:,:)

  real(RP), public,  allocatable :: Kv_coef  (:)                   ! vertical diffusion coefficient at cell center
  real(RP), private, allocatable :: Kv_coef_h(:)                   ! vertical diffusion coefficient at cell wall

  real(RP), public,  allocatable :: divdamp_coef   (:,:,:)         ! divergence damping coefficient at cell center
  real(RP), private, allocatable :: divdamp_coef_pl(:,:,:)
  integer, private              :: lap_order_divdamp = 2          ! laplacian order
  real(RP), private              :: divdamp_coef_v    = 0.0_RP

  real(RP), public,  allocatable :: divdamp_2d_coef   (:,:,:)      ! divergence damping coefficient at cell center
  real(RP), private, allocatable :: divdamp_2d_coef_pl(:,:,:)
  integer, private              :: lap_order_divdamp_2d = 1       ! laplacian order

  logical, private              :: dep_hgrid = .false.            ! depend on the horizontal grid spacing?
  real(RP), private              :: AREA_ave                       ! averaged grid area

  logical, private              :: smooth_1var = .true.           ! should be false for stretched grid [add] S.Iga 20120721

  logical, private              :: deep_effect = .false.
  real(RP), private, allocatable :: Kh_deep_factor       (:)
  real(RP), private, allocatable :: Kh_deep_factor_h     (:)
  real(RP), private, allocatable :: Kh_lap1_deep_factor  (:)
  real(RP), private, allocatable :: Kh_lap1_deep_factor_h(:)
  real(RP), private, allocatable :: divdamp_deep_factor  (:)

  logical, private              :: debug = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine numfilter_setup
    use mod_adm, only: &
       ADM_CTL_FID,   &
       ADM_proc_stop, &
       ADM_GLEVEL,    &
       ADM_kall
    use mod_cnst, only: &
       PI     => CNST_PI, &
       RADIUS => CNST_ERADIUS
    use mod_grd, only: &
       GRD_gz,   &
       GRD_gzh
    implicit none

    ! rayleigh damping
    real(RP)                 :: alpha_r         = 0.0_RP                 ! coefficient for rayleigh damping
    real(RP)                 :: ZD              = 25000.0_RP             ! lower limit of rayleigh damping [m]
    ! horizontal diffusion
    character(len=ADM_NSYS) :: hdiff_type      = 'NONDIM_COEF'        ! diffusion type
    real(RP)                 :: gamma_h         = 1.0_RP / 16.0_RP / 10.0_RP ! coefficient    for horizontal diffusion
    real(RP)                 :: tau_h           = 160000.0_RP            ! e-folding time for horizontal diffusion [sec]
    ! horizontal diffusion (1st order laplacian)
    character(len=ADM_NSYS) :: hdiff_type_lap1 = 'DIRECT'             ! diffusion type
    real(RP)                 :: gamma_h_lap1    = 0.0_RP                 ! height-dependent gamma_h but 1st-order laplacian
    real(RP)                 :: tau_h_lap1      = 160000.0_RP            ! height-dependent tau_h   but 1st-order laplacian [sec]
    real(RP)                 :: ZD_hdiff_lap1   = 25000.0_RP             ! lower limit of horizontal diffusion [m]
    ! vertical diffusion
    real(RP)                 :: gamma_v         = 0.0_RP                 ! coefficient of vertical diffusion
    ! 3D divergence damping
    character(len=ADM_NSYS) :: divdamp_type    = 'NONDIM_COEF'        ! damping type
    real(RP)                 :: alpha_d         = 0.0_RP                 ! coefficient    for divergence damping
    real(RP)                 :: tau_d           = 132800.0_RP            ! e-folding time for divergence damping
    real(RP)                 :: alpha_dv        = 0.0_RP                 ! vertical coefficient
    ! 2D divergence damping
    character(len=ADM_NSYS) :: divdamp_2d_type = 'NONDIM_COEF'        ! damping type
    real(RP)                 :: alpha_d_2d      = 0.0_RP                 ! coefficient    for divergence damping
    real(RP)                 :: tau_d_2d        = 1328000.0_RP           ! e-folding time for divergence damping [sec]
    real(RP)                 :: ZD_d_2d         = 25000.0_RP             ! lower limit of divergence damping [m]

    namelist / NUMFILTERPARAM / &
       alpha_r,              &
       ZD,                   &
       rayleigh_damp_only_w, &
       hdiff_type,           &
       lap_order_hdiff,      &
       gamma_h,              &
       tau_h,                &
       ZD_hdiff_nl,          &
       hdiff_fact_rho,       &
       hdiff_fact_q,         &
       Kh_coef_minlim,       &
       Kh_coef_maxlim,       &
       hdiff_type_lap1,      &
       gamma_h_lap1,         &
       tau_h_lap1,           &
       ZD_hdiff_lap1,        &
       gamma_v,              &
       divdamp_type,         &
       lap_order_divdamp,    &
       alpha_d,              &
       tau_d,                &
       alpha_dv,             &
       divdamp_2d_type,      &
       lap_order_divdamp_2d, &
       alpha_d_2d,           &
       tau_d_2d,             &
       ZD_d_2d,              &
       dep_hgrid,            &
       smooth_1var,          &
       deep_effect,          &
       debug

    real(RP) :: global_area, global_grid

    integer :: k
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[numfilter]/Category[nhm dynamics]'
    rewind(ADM_CTL_FID)
    read(ADM_CTL_FID,nml=NUMFILTERPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(ADM_LOG_FID,*) '*** NUMFILTERPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,          *) 'xxx Not appropriate names in namelist NUMFILTERPARAM. STOP.'
       write(ADM_LOG_FID,*) 'xxx Not appropriate names in namelist NUMFILTERPARAM. STOP.'
       call ADM_proc_stop
    endif
    write(ADM_LOG_FID,nml=NUMFILTERPARAM)

    global_area = 4.0_RP * PI * RADIUS * RADIUS
    global_grid = 10.0_RP * 4.0_RP**ADM_GLEVEL
    AREA_ave = global_area / global_grid



    call numfilter_rayleigh_damping_setup( alpha_r, & ! [IN]
                                           ZD       ) ! [IN]

    call numfilter_hdiffusion_setup( hdiff_type,      & ! [IN]
                                     dep_hgrid,       & ! [IN]
                                     smooth_1var,     & ! [IN]
                                     lap_order_hdiff, & ! [IN]
                                     gamma_h,         & ! [IN]
                                     tau_h,           & ! [IN]
                                     hdiff_type_lap1, & ! [IN]
                                     gamma_h_lap1,    & ! [IN]
                                     tau_h_lap1,      & ! [IN]
                                     ZD_hdiff_lap1    ) ! [IN]

    call numfilter_vdiffusion_setup( gamma_v ) ! [IN]

    call numfilter_divdamp_setup( divdamp_type,      & ! [IN]
                                  dep_hgrid,         & ! [IN]
                                  smooth_1var,       & ! [IN]
                                  lap_order_divdamp, & ! [IN]
                                  alpha_d,           & ! [IN]
                                  tau_d,             & ! [IN]
                                  alpha_dv           ) ! [IN]

    call numfilter_divdamp_2d_setup( divdamp_2d_type,      & ! [IN]
                                     dep_hgrid,            & ! [IN]
                                     lap_order_divdamp_2d, & ! [IN]
                                     alpha_d_2d,           & ! [IN]
                                     tau_d_2d,             & ! [IN]
                                     ZD_d_2d               ) ! [IN]

    allocate( Kh_deep_factor       (ADM_kall) )
    allocate( Kh_deep_factor_h     (ADM_kall) )
    allocate( Kh_lap1_deep_factor  (ADM_kall) )
    allocate( Kh_lap1_deep_factor_h(ADM_kall) )
    allocate( divdamp_deep_factor  (ADM_kall) )
    Kh_deep_factor       (:) = 0.0_RP
    Kh_deep_factor_h     (:) = 0.0_RP
    Kh_lap1_deep_factor  (:) = 0.0_RP
    Kh_lap1_deep_factor_h(:) = 0.0_RP
    divdamp_deep_factor  (:) = 0.0_RP

    if ( deep_effect ) then
       write(ADM_LOG_FID,*) 'xxx this feature is tentatively suspended. stop.'
       call ADM_proc_stop
       do k = 1, ADM_kall
          Kh_deep_factor       (k) = ( (GRD_gz (k)+RADIUS) / RADIUS )**(2*lap_order_hdiff)
          Kh_deep_factor_h     (k) = ( (GRD_gzh(k)+RADIUS) / RADIUS )**(2*lap_order_hdiff)
          Kh_lap1_deep_factor  (k) = ( (GRD_gz (k)+RADIUS) / RADIUS )**2
          Kh_lap1_deep_factor_h(k) = ( (GRD_gzh(k)+RADIUS) / RADIUS )**2
          divdamp_deep_factor  (k) = ( (GRD_gz (k)+RADIUS) / RADIUS )**(2*lap_order_divdamp)
       enddo
    endif

    return
  end subroutine numfilter_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for rayleigh damping
  subroutine numfilter_rayleigh_damping_setup( &
       alpha, &
       zlimit )
    use mod_adm, only: &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       EPS => CNST_EPS_ZERO
    use mod_grd, only: &
       GRD_htop, &
       GRD_gz,   &
       GRD_gzh
    implicit none

    real(RP), intent(in) :: alpha  ! coefficient for rayleigh damping
    real(RP), intent(in) :: zlimit ! lower limit of rayleigh damping [m]

    real(RP) :: fact(ADM_kall)

    integer :: k
    !---------------------------------------------------------------------------

    if ( alpha > 0.0_RP ) NUMFILTER_DOrayleigh = .true.

    allocate( rayleigh_coef  (ADM_kall) )
    allocate( rayleigh_coef_h(ADM_kall) )

    call height_factor( ADM_kall, GRD_gz(:), GRD_htop, zlimit, fact(:) )

    rayleigh_coef(:) = alpha * fact(:)

    call height_factor( ADM_kall, GRD_gzh(:), GRD_htop, zlimit, fact(:) )

    rayleigh_coef_h(:) = alpha * fact(:)

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '-----   Rayleigh damping   -----'

    if ( NUMFILTER_DOrayleigh ) then
       if ( debug ) then
          write(ADM_LOG_FID,*) '    z[m]      ray.coef   e-time(2DX)'
          k = ADM_kmax + 1
          write(ADM_LOG_FID,'(1x,F8.2,3E14.6)') GRD_gzh(k), rayleigh_coef_h(k), 1.0_RP/( rayleigh_coef_h(k)+EPS )
          do k = ADM_kmax, ADM_kmin, -1
             write(ADM_LOG_FID,'(1x,F8.2,3E14.6)') GRD_gz (k), rayleigh_coef  (k), 1.0_RP/( rayleigh_coef  (k)+EPS )
             write(ADM_LOG_FID,'(1x,F8.2,3E14.6)') GRD_gzh(k), rayleigh_coef_h(k), 1.0_RP/( rayleigh_coef_h(k)+EPS )
          enddo
       else
          write(ADM_LOG_FID,*) '=> used.'
       endif
    else
       write(ADM_LOG_FID,*) '=> not used.'
    endif

    return
  end subroutine numfilter_rayleigh_damping_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for horizontal numerical diffusion
  subroutine numfilter_hdiffusion_setup( &
       hdiff_type,      &
       dep_hgrid,       &
       smooth_1var,     &
       lap_order,       &
       gamma,           &
       tau,             &
       hdiff_type_lap1, &
       gamma_lap1,      &
       tau_lap1,        &
       zlimit_lap1      )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       PI  => CNST_PI,       &
       EPS => CNST_EPS_ZERO
    use mod_grd, only: &
       GRD_htop, &
       GRD_gz
    use mod_gmtr, only: &
       GMTR_area,    &
       GMTR_area_pl
    use mod_time, only: &
       TIME_DTL
    use mod_gtl, only: &
       GTL_max_k, &
       GTL_min_k
    use mod_runconf, only: &
       DYN_DIV_NUM
    implicit none

    character(len=*), intent(in) :: hdiff_type      ! type of horizontal diffusion
    logical,          intent(in) :: dep_hgrid       ! depend on each horizontal grid?
    logical,          intent(in) :: smooth_1var     ! apply smoothing to coef?
    integer,          intent(in) :: lap_order       ! laplacian order
    real(RP),          intent(in) :: gamma           ! coefficient    for horizontal diffusion
    real(RP),          intent(in) :: tau             ! e-folding time for horizontal diffusion
    character(len=*), intent(in) :: hdiff_type_lap1 ! type of horizontal diffusion (lap1)
    real(RP),          intent(in) :: gamma_lap1      ! coefficient    for horizontal diffusion (lap1)
    real(RP),          intent(in) :: tau_lap1        ! e-folding time for horizontal diffusion (lap1)
    real(RP),          intent(in) :: zlimit_lap1     ! lower limit of horizontal diffusion (lap1) [m]

    real(RP) :: fact(ADM_kall)

    real(RP) :: e_fold_time   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: e_fold_time_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: coef_max, coef_min
    real(RP) :: eft_max,  eft_min

    real(RP) :: large_step_dt

    integer :: k, l
    !---------------------------------------------------------------------------

    allocate( Kh_coef   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    Kh_coef   (:,:,:) = 0.0_RP
    Kh_coef_pl(:,:,:) = 0.0_RP

    if ( hdiff_type == "DIRECT" ) then
       if( gamma > 0.0_RP ) NUMFILTER_DOhorizontaldiff = .true.

       ! gamma is an absolute value.
       Kh_coef   (:,:,:) = gamma
       Kh_coef_pl(:,:,:) = gamma

    elseif( hdiff_type == "NONDIM_COEF" ) then
       if( gamma > 0.0_RP ) NUMFILTER_DOhorizontaldiff = .true.

       large_step_dt = TIME_DTL / real(DYN_DIV_NUM,kind=RP)

       ! gamma is a non-dimensional number.
       if ( dep_hgrid ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
             Kh_coef(:,k,l) = gamma / large_step_dt * GMTR_area(:,l)**lap_order
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                Kh_coef_pl(:,k,l) = gamma / large_step_dt * GMTR_area_pl(:,l)**lap_order
             enddo
             enddo
          endif
       else
          Kh_coef   (:,:,:) = gamma / large_step_dt * AREA_ave**lap_order
          Kh_coef_pl(:,:,:) = gamma / large_step_dt * AREA_ave**lap_order
       endif

    elseif( hdiff_type == "E_FOLD_TIME" ) then
       if( tau > 0.0_RP ) NUMFILTER_DOhorizontaldiff = .true.

       ! tau is e-folding time for 2*dx waves.
       if ( dep_hgrid ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
             Kh_coef(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**(2*lap_order) / ( tau+EPS )
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                Kh_coef_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**(2*lap_order) / ( tau+EPS )
             enddo
             enddo
          endif
       else
          Kh_coef   (:,:,:) = ( sqrt(AREA_ave)/PI )**(2*lap_order) / ( tau+EPS )
          Kh_coef_pl(:,:,:) = ( sqrt(AREA_ave)/PI )**(2*lap_order) / ( tau+EPS )
       endif

    elseif( hdiff_type  == "NONLINEAR1" ) then
       NUMFILTER_DOhorizontaldiff = .true.
       hdiff_nonlinear            = .true.

       Kh_coef   (:,:,:) = -999.0_RP
       Kh_coef_pl(:,:,:) = -999.0_RP
    endif

    if (       hdiff_type /= "DIRECT"     &
         .AND. hdiff_type /= "NONLINEAR1" ) then

       if ( smooth_1var ) then ! iga 20120721 (add if)
          call numfilter_smooth_1var( Kh_coef(:,:,:), Kh_coef_pl(:,:,:) )
       endif
       Kh_coef(:,:,:) = max( Kh_coef(:,:,:), Kh_coef_minlim )

    endif

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '-----   Horizontal numerical diffusion   -----'
    if ( NUMFILTER_DOhorizontaldiff ) then
       if ( .NOT. hdiff_nonlinear ) then
          if ( debug ) then
             do l = 1, ADM_lall
             do k = 1, ADM_kall
                e_fold_time(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**(2*lap_order) &
                                   / ( Kh_coef(:,k,l)+EPS )
             enddo
             enddo

             if ( ADM_have_pl ) then
                do l = 1, ADM_lall_pl
                do k = 1, ADM_kall
                   e_fold_time_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**(2*lap_order) &
                                         / ( Kh_coef_pl(:,k,l)+EPS )
                enddo
                enddo
             endif

             write(ADM_LOG_FID,*) '    z[m]      max coef      min coef  max eft(2DX)  min eft(2DX)'
             do k = ADM_kmax, ADM_kmin, -1
                eft_max  = GTL_max_k( e_fold_time, e_fold_time_pl, k )
                eft_min  = GTL_min_k( e_fold_time, e_fold_time_pl, k )
                coef_max = GTL_max_k( Kh_coef, Kh_coef_pl, k )
                coef_min = GTL_min_k( Kh_coef, Kh_coef_pl, k )
                write(ADM_LOG_FID,'(1x,F8.2,4E14.6)') GRD_gz(k), coef_min, coef_max, eft_max, eft_min
             enddo
          else
             write(ADM_LOG_FID,*) '=> used.'
          endif
       else
          write(ADM_LOG_FID,*) '=> Nonlinear filter is used.'
       endif
    else
       write(ADM_LOG_FID,*) '=> not used.'
    endif



    allocate( Kh_coef_lap1   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( Kh_coef_lap1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    Kh_coef_lap1    = 0.0_RP
    Kh_coef_lap1_pl = 0.0_RP

    if ( hdiff_type_lap1 == "DIRECT" ) then
       if( gamma_lap1 > 0.0_RP ) NUMFILTER_DOhorizontaldiff_lap1 = .true.

       ! gamma is an absolute value.
       Kh_coef_lap1   (:,:,:) = gamma_lap1
       Kh_coef_lap1_pl(:,:,:) = gamma_lap1

    elseif( hdiff_type_lap1 == "NONDIM_COEF" ) then
       if( gamma_lap1 > 0.0_RP ) NUMFILTER_DOhorizontaldiff_lap1 = .true.

       large_step_dt = TIME_DTL / real(DYN_DIV_NUM,kind=RP)

       ! gamma is a non-dimensional number.
       if ( dep_hgrid ) then

          do l = 1, ADM_lall
          do k = 1, ADM_kall
             Kh_coef_lap1(:,k,l) = gamma_lap1 / large_step_dt * GMTR_area(:,l)
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                Kh_coef_lap1_pl(:,k,l) = gamma_lap1 / large_step_dt * GMTR_area_pl(:,l)
             enddo
             enddo
          endif

       else
          Kh_coef_lap1   (:,:,:) = gamma_lap1 / large_step_dt * AREA_ave
          Kh_coef_lap1_pl(:,:,:) = gamma_lap1 / large_step_dt * AREA_ave
       endif

    elseif( hdiff_type_lap1 == "E_FOLD_TIME" ) then
       if( tau_lap1 > 0.0_RP ) NUMFILTER_DOhorizontaldiff_lap1 = .true.

       ! tau is e-folding time for 2*dx waves.
       if ( dep_hgrid ) then

          do l = 1, ADM_lall
          do k = 1, ADM_kall
             Kh_coef_lap1(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**2 / ( tau_lap1+EPS )
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                Kh_coef_lap1_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**2 / ( tau_lap1+EPS )
             enddo
             enddo
          endif

       else
          Kh_coef_lap1   (:,:,:) = ( sqrt(AREA_ave)/PI )**2 / ( tau_lap1+EPS )
          Kh_coef_lap1_pl(:,:,:) = ( sqrt(AREA_ave)/PI )**2 / ( tau_lap1+EPS )
       endif

    endif

    call height_factor( ADM_kall, GRD_gz(:), GRD_htop, zlimit_lap1, fact(:) )

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       Kh_coef_lap1(:,k,l) = Kh_coef_lap1(:,k,l) * fact(k)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          Kh_coef_lap1_pl(:,k,l) = Kh_coef_lap1_pl(:,k,l) * fact(k)
       enddo
       enddo
    endif

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '-----   Horizontal numerical diffusion (1st order laplacian)   -----'
    if ( NUMFILTER_DOhorizontaldiff_lap1 ) then
       if ( debug ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
             e_fold_time(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**2 / ( Kh_coef_lap1(:,k,l)+EPS )
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                e_fold_time_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**2 / ( Kh_coef_lap1_pl(:,k,l)+EPS )
             enddo
             enddo
          endif

          write(ADM_LOG_FID,*) '    z[m]      max coef      min coef  max eft(2DX)  min eft(2DX)'
          do k = ADM_kmax, ADM_kmin, -1
             eft_max  = GTL_max_k( e_fold_time,  e_fold_time_pl,  k )
             eft_min  = GTL_min_k( e_fold_time,  e_fold_time_pl,  k )
             coef_max = GTL_max_k( Kh_coef_lap1, Kh_coef_lap1_pl, k )
             coef_min = GTL_min_k( Kh_coef_lap1, Kh_coef_lap1_pl, k )
             write(ADM_LOG_FID,'(1x,F8.2,4E14.6)') GRD_gz(k), coef_min, coef_max, eft_max, eft_min
          enddo
       else
          write(ADM_LOG_FID,*) '=> used.'
       endif
    else
       write(ADM_LOG_FID,*) '=> not used.'
    endif

    return
  end subroutine numfilter_hdiffusion_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for vertical numerical diffusion
  subroutine numfilter_vdiffusion_setup( &
       gamma )
    use mod_adm, only: &
       ADM_kall, &
       ADM_kmin, &
       ADM_kmax
    use mod_cnst, only: &
       PI  => CNST_PI,       &
       EPS => CNST_EPS_ZERO
    use mod_grd, only: &
       GRD_gz,   &
       GRD_gzh,  &
       GRD_dgz,  &
       GRD_dgzh
    use mod_time, only: &
       TIME_DTL
    use mod_runconf, only: &
       DYN_DIV_NUM
    implicit none

    real(RP), intent(in) :: gamma ! coefficient for vertical diffusion

    real(RP) :: large_step_dt

    integer :: k
    !---------------------------------------------------------------------------

    if ( gamma > 0.0_RP ) NUMFILTER_DOverticaldiff = .true.

    allocate( Kv_coef  (ADM_kall) )
    allocate( Kv_coef_h(ADM_kall) )

    large_step_dt = TIME_DTL / real(DYN_DIV_NUM,kind=RP)

    ! 6th order vertical numerical diffusion
    Kv_coef  (:) = gamma * GRD_dgz (:)**6 / large_step_dt
    Kv_coef_h(:) = gamma * GRD_dgzh(:)**6 / large_step_dt

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '-----   Vertical numerical diffusion   -----'
    if ( NUMFILTER_DOverticaldiff ) then
       if ( debug ) then
          write(ADM_LOG_FID,*) '    z[m]          coef   e-time(2DX)'
          k = ADM_kmax + 1
          write(ADM_LOG_FID,'(1x,F8.2,3E14.6)') GRD_gzh(k), Kv_coef_h(k), (GRD_dgzh(k)/PI)**6 / ( Kv_coef_h(k)+EPS )
          do k = ADM_kmax, ADM_kmin, -1
             write(ADM_LOG_FID,'(1x,F8.2,3E14.6)') GRD_gzh(k), Kv_coef_h(k), (GRD_dgzh(k)/PI)**6 / ( Kv_coef_h(k)+EPS )
             write(ADM_LOG_FID,'(1x,F8.2,3E14.6)') GRD_gz (k), Kv_coef  (k), (GRD_dgz (k)/PI)**6 / ( Kv_coef  (k)+EPS )
          enddo
       else
          write(ADM_LOG_FID,*) '=> used.'
       endif
    else
       write(ADM_LOG_FID,*) '=> not used.'
    endif

    return
  end subroutine numfilter_vdiffusion_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for 3D divergence damping
  subroutine numfilter_divdamp_setup( &
       divdamp_type, &
       dep_hgrid,    &
       smooth_1var,  &
       lap_order,    &
       alpha,        &
       tau,          &
       alpha_v       )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       PI    => CNST_PI,       &
       EPS   => CNST_EPS_ZERO, &
       SOUND => CNST_SOUND
    use mod_grd, only: &
       GRD_gz
    use mod_gmtr, only: &
       GMTR_area,    &
       GMTR_area_pl
    use mod_time, only: &
       TIME_DTS
    use mod_gtl, only: &
       GTL_max_k, &
       GTL_min_k
    use mod_runconf, only: &
       DYN_DIV_NUM
    implicit none

    character(len=*), intent(in) :: divdamp_type ! type of divergence damping
    logical,          intent(in) :: dep_hgrid    ! depend on each horizontal grid?
    logical,          intent(in) :: smooth_1var  ! apply smoothing to coef?
    integer,          intent(in) :: lap_order    ! laplacian order
    real(RP),          intent(in) :: alpha        ! coefficient    for divergence damping
    real(RP),          intent(in) :: tau          ! e-folding time for divergence damping
    real(RP),          intent(in) :: alpha_v      ! coefficient    for divergence damping

    real(RP) :: e_fold_time   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: e_fold_time_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: coef
    real(RP) :: coef_max, coef_min
    real(RP) :: eft_max,  eft_min

    real(RP) :: small_step_dt

    integer :: k, l
    !---------------------------------------------------------------------------

    allocate( divdamp_coef   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( divdamp_coef_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    divdamp_coef    = 0.0_RP
    divdamp_coef_pl = 0.0_RP

    if ( divdamp_type == "DIRECT") then
       if( alpha > 0.0_RP ) NUMFILTER_DOdivdamp = .true.

       ! alpha_d is an absolute value.
       coef = alpha

       divdamp_coef   (:,:,:) = coef
       divdamp_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "NONDIM_COEF" ) then
       if( alpha > 0.0_RP ) NUMFILTER_DOdivdamp = .true.

       small_step_dt = TIME_DTS / real(DYN_DIV_NUM,kind=RP)

       ! alpha_d is a non-dimensional number.
       ! alpha_d * (c_s)^p * dt^{2p-1}
       coef = alpha * ( SOUND * SOUND )**lap_order * small_step_dt**(2*lap_order-1)

       divdamp_coef   (:,:,:) = coef
       divdamp_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "E_FOLD_TIME" ) then
       if( tau > 0.0_RP ) NUMFILTER_DOdivdamp = .true.

       ! tau_d is e-folding time for 2*dx.
       if ( dep_hgrid ) then

          do l = 1, ADM_lall
          do k = 1, ADM_kall
             divdamp_coef(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**(2*lap_order) / ( tau+EPS )
          enddo
          enddo
          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                divdamp_coef_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**(2*lap_order) / ( tau+EPS )
             enddo
             enddo
          endif

       else

          coef = ( sqrt(AREA_ave)/PI )**(2*lap_order) / ( tau+EPS )

          divdamp_coef   (:,:,:) = coef
          divdamp_coef_pl(:,:,:) = coef

       endif
    endif

    if ( divdamp_type /= "DIRECT" ) then
       if ( smooth_1var ) then ! iga 20120721 (add if)
          call numfilter_smooth_1var( divdamp_coef(:,:,:), divdamp_coef_pl(:,:,:) )
       endif
       divdamp_coef(:,:,:) = max( divdamp_coef(:,:,:), Kh_coef_minlim )
    endif

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '-----   3D divergence damping   -----'
    if ( NUMFILTER_DOdivdamp ) then
       if ( debug ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
             e_fold_time(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**(2*lap_order) &
                                / ( divdamp_coef(:,k,l)+EPS )
          enddo
          enddo

          e_fold_time_pl(:,:,:) = 0.0_RP

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                e_fold_time_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**(2*lap_order) &
                                      / ( divdamp_coef_pl(:,k,l)+EPS )
             enddo
             enddo
          endif

          write(ADM_LOG_FID,*) '    z[m]      max coef      min coef  max eft(2DX)  min eft(2DX)'
          do k = ADM_kmax, ADM_kmin, -1
             eft_max  = GTL_max_k( e_fold_time,  e_fold_time_pl,  k )
             eft_min  = GTL_min_k( e_fold_time,  e_fold_time_pl,  k )
             coef_max = GTL_max_k( divdamp_coef, divdamp_coef_pl, k )
             coef_min = GTL_min_k( divdamp_coef, divdamp_coef_pl, k )
             write(ADM_LOG_FID,'(1x,F8.2,4E14.6)') GRD_gz(k), coef_min, coef_max, eft_max, eft_min
          enddo
       else
          write(ADM_LOG_FID,*) '=> used.'
       endif
    else
       write(ADM_LOG_FID,*) '=> not used.'
    endif

    if( alpha_v > 0.0_RP ) NUMFILTER_DOdivdamp_v = .true.

    small_step_dt = TIME_DTS / real(DYN_DIV_NUM,kind=RP)

    divdamp_coef_v = -alpha_v * SOUND * SOUND * small_step_dt

    return
  end subroutine numfilter_divdamp_setup

  !-----------------------------------------------------------------------------
  !> setup coefficient for vertical numerical diffusion
  subroutine numfilter_divdamp_2d_setup( &
       divdamp_type, &
       dep_hgrid,    &
       lap_order,    &
       alpha,        &
       tau,          &
       zlimit        )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       PI    => CNST_PI,       &
       EPS   => CNST_EPS_ZERO, &
       SOUND => CNST_SOUND
    use mod_grd, only: &
       GRD_htop, &
       GRD_gz
    use mod_gmtr, only: &
       GMTR_area,    &
       GMTR_area_pl
    use mod_time, only: &
       TIME_DTS
    use mod_gtl, only: &
       GTL_max_k, &
       GTL_min_k
    use mod_runconf, only: &
       DYN_DIV_NUM
    implicit none

    character(len=*), intent(in) :: divdamp_type ! type of divergence damping
    logical,          intent(in) :: dep_hgrid    ! depend on each horizontal grid?
    integer,          intent(in) :: lap_order    ! laplacian order
    real(RP),          intent(in) :: alpha        ! coefficient    for divergence damping
    real(RP),          intent(in) :: tau          ! e-folding time for divergence damping
    real(RP),          intent(in) :: zlimit       ! lower limit of divergence damping [m]

    real(RP) :: fact(ADM_kall)

    real(RP) :: e_fold_time   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: e_fold_time_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: coef
    real(RP) :: coef_max, coef_min
    real(RP) :: eft_max,  eft_min

    real(RP) :: small_step_dt

    integer :: k, l
    !---------------------------------------------------------------------------

    allocate( divdamp_2d_coef   (ADM_gall,   ADM_kall,ADM_lall   ) )
    allocate( divdamp_2d_coef_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    divdamp_2d_coef    = 0.0_RP
    divdamp_2d_coef_pl = 0.0_RP

    if ( divdamp_type == "DIRECT" ) then
       if( alpha > 0.0_RP ) NUMFILTER_DOdivdamp_2d = .true.

       ! alpha is the absolute value.
       coef = alpha

       divdamp_2d_coef   (:,:,:) = coef
       divdamp_2d_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "NONDIM_COEF" ) then
       if( alpha > 0.0_RP ) NUMFILTER_DOdivdamp_2d = .true.

       small_step_dt = TIME_DTS / real(DYN_DIV_NUM,kind=RP)

       ! alpha is the non-dimensional number.
       ! alpha * (c_s)^p * dt^{2p-1}
       coef = alpha * ( SOUND * SOUND )**lap_order * small_step_dt**(2*lap_order-1)

       divdamp_2d_coef   (:,:,:) = coef
       divdamp_2d_coef_pl(:,:,:) = coef

    elseif( divdamp_type == "E_FOLD_TIME" ) then
       if( tau > 0.0_RP ) NUMFILTER_DOdivdamp_2d = .true.

       ! tau is e-folding time for 2*dx.
       if ( dep_hgrid ) then

          do l = 1, ADM_lall
          do k = 1, ADM_kall
             divdamp_2d_coef(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**(2*lap_order) / ( tau+EPS )
          enddo
          enddo
          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                divdamp_2d_coef_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**(2*lap_order) / ( tau+EPS )
             enddo
             enddo
          endif

       else

          coef = ( sqrt(AREA_ave)/PI )**(2*lap_order) / ( tau+EPS )

          divdamp_2d_coef   (:,:,:) = coef
          divdamp_2d_coef_pl(:,:,:) = coef

       endif
    endif

    call height_factor( ADM_kall, GRD_gz(:), GRD_htop, zlimit, fact(:) )

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       divdamp_2d_coef(:,k,l) = divdamp_2d_coef(:,k,l) * fact(k)
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          divdamp_2d_coef_pl(:,k,l) = divdamp_2d_coef_pl(:,k,l) * fact(k)
       enddo
       enddo
    endif

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '-----   2D divergence damping   -----'
    if ( NUMFILTER_DOdivdamp_2d ) then
       if ( debug ) then
          do l = 1, ADM_lall
          do k = 1, ADM_kall
             e_fold_time(:,k,l) = ( sqrt(GMTR_area(:,l))/PI )**(2*lap_order_divdamp) &
                                / ( divdamp_2d_coef(:,k,l)+EPS )
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
                e_fold_time_pl(:,k,l) = ( sqrt(GMTR_area_pl(:,l))/PI )**(2*lap_order_divdamp) &
                                      / ( divdamp_2d_coef_pl(:,k,l)+EPS )
             enddo
             enddo
          else
             e_fold_time_pl(:,:,:) = 0.0_RP
          endif

          write(ADM_LOG_FID,*) '    z[m]      max coef      min coef  max eft(2DX)  min eft(2DX)'
          do k = ADM_kmax, ADM_kmin, -1
             eft_max  = GTL_max_k( e_fold_time,  e_fold_time_pl,  k )
             eft_min  = GTL_min_k( e_fold_time,  e_fold_time_pl,  k )
             coef_max = GTL_max_k( divdamp_coef, divdamp_coef_pl, k )
             coef_min = GTL_min_k( divdamp_coef, divdamp_coef_pl, k )
             write(ADM_LOG_FID,'(1x,F8.2,4E14.6)') GRD_gz(k), coef_min, coef_max, eft_max, eft_min
          enddo
       else
          write(ADM_LOG_FID,*) '=> used.'
       endif
    else
       write(ADM_LOG_FID,*) '=> not used.'
    endif

    return
  end subroutine numfilter_divdamp_2d_setup

  !-----------------------------------------------------------------------------
  !> Rayleigh damping
  subroutine numfilter_rayleigh_damping( &
       rhog,    rhog_pl,    &
       vx,      vx_pl,      &
       vy,      vy_pl,      &
       vz,      vz_pl,      &
       w,       w_pl,       &
       frhogvx, frhogvx_pl, &
       frhogvy, frhogvy_pl, &
       frhogvz, frhogvz_pl, &
       frhogw,  frhogw_pl   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    implicit none

    real(RP), intent(in)    :: rhog      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: rhog_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vx        (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vx_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vy        (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vy_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vz        (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vz_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: w         (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: w_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogw    (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: coef
    integer :: g, k, l
    !---------------------------------------------------------------------------

    if( .NOT. NUMFILTER_DOrayleigh ) return

    call DEBUG_rapstart('____numfilter_rayleigh_damping')

    if ( .NOT. rayleigh_damp_only_w ) then
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          coef = rayleigh_coef(k) * rhog(g,k,l)

          frhogvx(g,k,l) = frhogvx(g,k,l) - coef * vx(g,k,l)
          frhogvy(g,k,l) = frhogvy(g,k,l) - coef * vy(g,k,l)
          frhogvz(g,k,l) = frhogvz(g,k,l) - coef * vz(g,k,l)
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             coef = rayleigh_coef(k) * rhog_pl(g,k,l)

             frhogvx_pl(g,k,l) = frhogvx_pl(g,k,l) - coef * vx_pl(g,k,l)
             frhogvy_pl(g,k,l) = frhogvy_pl(g,k,l) - coef * vy_pl(g,k,l)
             frhogvz_pl(g,k,l) = frhogvz_pl(g,k,l) - coef * vz_pl(g,k,l)
          enddo
          enddo
          enddo
       endif
    endif

    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax+1
    do g = 1, ADM_gall
       frhogw(g,k,l) = frhogw(g,k,l) &
                     - rayleigh_coef_h(k) * w(g,k,l) * ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
                                                       + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall_pl
          frhogw_pl(g,k,l) = frhogw_pl(g,k,l) &
                           - rayleigh_coef_h(k) * w_pl(g,k,l) * ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
                                                                + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
       enddo
       enddo
       enddo
    endif

    call DEBUG_rapend('____numfilter_rayleigh_damping')

    return
  end subroutine numfilter_rayleigh_damping

  !-----------------------------------------------------------------------------
  !> horizontal numerical diffusion
  subroutine numfilter_hdiffusion( &
       rhog,       rhog_pl,      &
       rho,        rho_pl,       &
       vx,         vx_pl,        &
       vy,         vy_pl,        &
       vz,         vz_pl,        &
       w,          w_pl,         &
       tem,        tem_pl,       &
       q,          q_pl,         &
       tendency,   tendency_pl,  &
       tendency_q, tendency_q_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CVdry => CNST_CV
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_htop, &
       GRD_gz
    use mod_oprt, only: &
       OPRT_horizontalize_vec, &
       OPRT_laplacian,         &
       OPRT_diffusion
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_time, only: &
       TIME_DTL
    use mod_runconf, only: &
       TRC_VMAX,     &
       TRC_ADV_TYPE, &
       DYN_DIV_NUM
    use mod_bsstate, only: &
       rho_bs,    &
       rho_bs_pl, &
       tem_bs,    &
       tem_bs_pl
    implicit none

    real(RP), intent(in)  :: rhog         (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhog_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rho          (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rho_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx           (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy           (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz           (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: w            (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: w_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: tem          (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: tem_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: q            (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP), intent(in)  :: q_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP), intent(out) :: tendency     (ADM_gall,   ADM_kall,ADM_lall   ,7)
    real(RP), intent(out) :: tendency_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl,7)
    real(RP), intent(out) :: tendency_q   (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP), intent(out) :: tendency_q_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: KH_coef_h        (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: KH_coef_h_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: KH_coef_lap1_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: KH_coef_lap1_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: vtmp        (ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(RP) :: vtmp_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: vtmp2       (ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(RP) :: vtmp2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP) :: qtmp        (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP) :: qtmp_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP) :: qtmp2       (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP) :: qtmp2_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: vtmp_lap1   (ADM_gall,   ADM_kall,ADM_lall   ,6)
    real(RP) :: vtmp_lap1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: qtmp_lap1   (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP) :: qtmp_lap1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: wk       (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: wk_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhog_h   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhog_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), parameter :: cfact = 2.0_RP
    real(RP), parameter :: T0    = 300.0_RP

    real(RP) :: fact  (ADM_kall)
    real(RP) :: kh_max(ADM_kall)
    real(RP) :: d2T_dx2, coef

    real(RP) :: large_step_dt

    integer :: g, k, l, nq, p
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____numfilter_hdiffusion')


    if ( hdiff_nonlinear ) then
       call height_factor( ADM_kall, GRD_gz(:), GRD_htop, ZD_hdiff_nl, fact(:) )

       kh_max(:) = ( 1.0_RP - fact(:) ) * Kh_coef_maxlim &
                 + (        fact(:) ) * Kh_coef_minlim
    endif

    rhog_h(:,ADM_kmin-1,:) = 0.0_RP
    do l = 1, ADM_lall
    do k = ADM_kmin, ADM_kmax+1
    do g = 1, ADM_gall
       rhog_h(g,k,l) = ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
                       + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
    enddo
    enddo
    enddo

    rhog_h_pl(:,ADM_kmin-1,:) = 0.0_RP
    do l = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmax+1
    do g = 1, ADM_gall_pl
       rhog_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
                          + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
    enddo
    enddo
    enddo

    vtmp   (:,:,:,1) = vx    (:,:,:)
    vtmp   (:,:,:,2) = vy    (:,:,:)
    vtmp   (:,:,:,3) = vz    (:,:,:)
    vtmp   (:,:,:,4) = w     (:,:,:)
    vtmp   (:,:,:,5) = tem   (:,:,:) - tem_bs   (:,:,:)
    vtmp   (:,:,:,6) = rho   (:,:,:) - rho_bs   (:,:,:)

    vtmp_pl(:,:,:,1) = vx_pl (:,:,:)
    vtmp_pl(:,:,:,2) = vy_pl (:,:,:)
    vtmp_pl(:,:,:,3) = vz_pl (:,:,:)
    vtmp_pl(:,:,:,4) = w_pl  (:,:,:)
    vtmp_pl(:,:,:,5) = tem_pl(:,:,:) - tem_bs_pl(:,:,:)
    vtmp_pl(:,:,:,6) = rho_pl(:,:,:) - rho_bs_pl(:,:,:)

    ! copy beforehand
    if ( NUMFILTER_DOhorizontaldiff_lap1 ) then
       vtmp_lap1   (:,:,:,:) = vtmp   (:,:,:,:)
       vtmp_lap1_pl(:,:,:,:) = vtmp_pl(:,:,:,:)
    endif

    ! high order laplacian
    do p = 1, lap_order_hdiff
       ! for momentum
       call OPRT_laplacian( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & ! [OUT]
                            vtmp (:,:,:,1), vtmp_pl (:,:,:,1)  ) ! [IN]

       call OPRT_laplacian( vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & ! [OUT]
                            vtmp (:,:,:,2), vtmp_pl (:,:,:,2)  ) ! [IN]

       call OPRT_laplacian( vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & ! [OUT]
                            vtmp (:,:,:,3), vtmp_pl (:,:,:,3)  ) ! [IN]

       call OPRT_laplacian( vtmp2(:,:,:,4), vtmp2_pl(:,:,:,4), & ! [OUT]
                            vtmp (:,:,:,4), vtmp_pl (:,:,:,4)  ) ! [IN]

       ! for scalar
       if ( p == lap_order_hdiff ) then

          if ( hdiff_nonlinear ) then
             large_step_dt = TIME_DTL / real(DYN_DIV_NUM,kind=RP)

             do l = 1, ADM_lall
             do k = 1, ADM_kall
             do g = 1, ADM_gall
                d2T_dx2 = abs(vtmp(g,k,l,5)) / T0 * AREA_ave
                coef    = cfact * ( AREA_ave * AREA_ave ) / large_step_dt * d2T_dx2

                KH_coef(g,k,l) = max( min( coef, Kh_max(k) ), Kh_coef_minlim )
             enddo
             enddo
             enddo

             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                d2T_dx2 = abs(vtmp_pl(g,k,l,5)) / T0 * AREA_ave
                coef    = cfact * ( AREA_ave * AREA_ave ) / large_step_dt * d2T_dx2

                KH_coef_pl(g,k,l) = max( min( coef, Kh_max(k) ), Kh_coef_minlim )
             enddo
             enddo
             enddo

             do l = 1, ADM_lall
                do k = ADM_kmin+1, ADM_kmax
                   KH_coef_h(:,k,l) = 0.5_RP * ( KH_coef(:,k,l) + KH_coef(:,k-1,l) )
                enddo
                KH_coef_h(:,ADM_kmin,l) = 0.0_RP
             enddo

             do l = 1, ADM_lall_pl
                do k = ADM_kmin+1, ADM_kmax
                   KH_coef_h_pl(:,k,l) = 0.5_RP * ( KH_coef_pl(:,k,l) + KH_coef_pl(:,k-1,l) )
                enddo
                KH_coef_h_pl(:,ADM_kmin,l) = 0.0_RP
             enddo
          else
             KH_coef_h   (:,:,:) = KH_coef   (:,:,:)
             KH_coef_h_pl(:,:,:) = KH_coef_pl(:,:,:)
          endif ! nonlinear1

          wk   (:,:,:) = rhog   (:,:,:) * CVdry * KH_coef   (:,:,:)
          wk_pl(:,:,:) = rhog_pl(:,:,:) * CVdry * KH_coef_pl(:,:,:)

          call OPRT_diffusion( vtmp2(:,:,:,5), vtmp2_pl(:,:,:,5), & ! [OUT]
                               vtmp (:,:,:,5), vtmp_pl (:,:,:,5), & ! [IN]
                               wk   (:,:,:)  , wk_pl   (:,:,:)    ) ! [IN]

          wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_rho * KH_coef   (:,:,:)
          wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_rho * KH_coef_pl(:,:,:)

          call OPRT_diffusion( vtmp2(:,:,:,6), vtmp2_pl(:,:,:,6), & ! [OUT]
                               vtmp (:,:,:,6), vtmp_pl (:,:,:,6), & ! [IN]
                               wk   (:,:,:)  , wk_pl   (:,:,:)    ) ! [IN]
       else
          call OPRT_laplacian( vtmp2(:,:,:,5), vtmp2_pl(:,:,:,5), & ! [OUT]
                               vtmp (:,:,:,5), vtmp_pl (:,:,:,5)  ) ! [IN]

          call OPRT_laplacian( vtmp2(:,:,:,6), vtmp2_pl(:,:,:,6), & ! [OUT]
                               vtmp (:,:,:,6), vtmp_pl (:,:,:,6)  ) ! [IN]
       endif

       vtmp   (:,:,:,:) = -vtmp2   (:,:,:,:)
       vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

       call COMM_data_transfer( vtmp, vtmp_pl )

    enddo ! laplacian order loop

    !--- 1st order laplacian filter
    if ( NUMFILTER_DOhorizontaldiff_lap1 ) then

       KH_coef_lap1_h   (:,:,:) = KH_coef_lap1   (:,:,:)
       KH_coef_lap1_h_pl(:,:,:) = KH_coef_lap1_pl(:,:,:)

       call OPRT_laplacian( vtmp2    (:,:,:,1), vtmp2_pl    (:,:,:,1), & ! [OUT]
                            vtmp_lap1(:,:,:,1), vtmp_lap1_pl(:,:,:,1)  ) ! [IN]

       call OPRT_laplacian( vtmp2    (:,:,:,2), vtmp2_pl    (:,:,:,2), & ! [OUT]
                            vtmp_lap1(:,:,:,2), vtmp_lap1_pl(:,:,:,2)  ) ! [IN]

       call OPRT_laplacian( vtmp2    (:,:,:,3), vtmp2_pl    (:,:,:,3), & ! [OUT]
                            vtmp_lap1(:,:,:,3), vtmp_lap1_pl(:,:,:,3)  ) ! [IN]

       call OPRT_laplacian( vtmp2    (:,:,:,4), vtmp2_pl    (:,:,:,4), & ! [OUT]
                            vtmp_lap1(:,:,:,4), vtmp_lap1_pl(:,:,:,4)  ) ! [IN]

       wk   (:,:,:) = rhog   (:,:,:) * CVdry * KH_coef_lap1   (:,:,:)
       wk_pl(:,:,:) = rhog_pl(:,:,:) * CVdry * KH_coef_lap1_pl(:,:,:)

       call OPRT_diffusion( vtmp2    (:,:,:,5), vtmp2_pl    (:,:,:,5), & ! [OUT]
                            vtmp_lap1(:,:,:,5), vtmp_lap1_pl(:,:,:,5), & ! [IN]
                            wk       (:,:,:),   wk_pl       (:,:,:)    ) ! [IN]

       wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_rho * KH_coef_lap1   (:,:,:)
       wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_rho * KH_coef_lap1_pl(:,:,:)

       call OPRT_diffusion( vtmp2    (:,:,:,6), vtmp2_pl    (:,:,:,6), & ! [OUT]
                            vtmp_lap1(:,:,:,6), vtmp_lap1_pl(:,:,:,6), & ! [IN]
                            wk       (:,:,:),   wk_pl       (:,:,:)    ) ! [IN]

       vtmp_lap1   (:,:,:,:) = -vtmp2   (:,:,:,:)
       vtmp_lap1_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

       call COMM_data_transfer( vtmp_lap1, vtmp_lap1_pl )
    else
       KH_coef_lap1_h   (:,:,:) = 0.0_RP
       KH_coef_lap1_h_pl(:,:,:) = 0.0_RP

       vtmp_lap1   (:,:,:,:) = 0.0_RP
       vtmp_lap1_pl(:,:,:,:) = 0.0_RP
    endif

    !--- Update tendency
    do l = 1, ADM_lall
    do k = 1, ADM_kall
       do g = 1, ADM_gall
          tendency(g,k,l,I_RHOGVX) = - ( vtmp     (g,k,l,1) * KH_coef     (g,k,l) &
                                       + vtmp_lap1(g,k,l,1) * KH_coef_lap1(g,k,l) ) * rhog(g,k,l)
          tendency(g,k,l,I_RHOGVY) = - ( vtmp     (g,k,l,2) * KH_coef     (g,k,l) &
                                       + vtmp_lap1(g,k,l,2) * KH_coef_lap1(g,k,l) ) * rhog(g,k,l)
          tendency(g,k,l,I_RHOGVZ) = - ( vtmp     (g,k,l,3) * KH_coef     (g,k,l) &
                                       + vtmp_lap1(g,k,l,3) * KH_coef_lap1(g,k,l) ) * rhog(g,k,l)
          tendency(g,k,l,I_RHOGW ) = - ( vtmp     (g,k,l,4) * KH_coef_h     (g,k,l) &
                                       + vtmp_lap1(g,k,l,4) * KH_coef_lap1_h(g,k,l) ) * rhog_h(g,k,l)
       enddo

       do g = 1, ADM_gall
          tendency(g,k,l,I_RHOGE   ) = - ( vtmp(g,k,l,5) + vtmp_lap1(g,k,l,5) )
          tendency(g,k,l,I_RHOGETOT) = - ( vtmp(g,k,l,5) + vtmp_lap1(g,k,l,5) )
          tendency(g,k,l,I_RHOG    ) = - ( vtmp(g,k,l,6) + vtmp_lap1(g,k,l,6) )
       enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             tendency_pl(g,k,l,I_RHOGVX) = - ( vtmp_pl     (g,k,l,1) * KH_coef_pl     (g,k,l) &
                                             + vtmp_lap1_pl(g,k,l,1) * KH_coef_lap1_pl(g,k,l) ) * rhog_pl(g,k,l)
             tendency_pl(g,k,l,I_RHOGVY) = - ( vtmp_pl     (g,k,l,2) * KH_coef_pl     (g,k,l) &
                                             + vtmp_lap1_pl(g,k,l,2) * KH_coef_lap1_pl(g,k,l) ) * rhog_pl(g,k,l)
             tendency_pl(g,k,l,I_RHOGVZ) = - ( vtmp_pl     (g,k,l,3) * KH_coef_pl     (g,k,l) &
                                             + vtmp_lap1_pl(g,k,l,3) * KH_coef_lap1_pl(g,k,l) ) * rhog_pl(g,k,l)
             tendency_pl(g,k,l,I_RHOGW ) = - ( vtmp_pl     (g,k,l,4) * KH_coef_h_pl     (g,k,l) &
                                             + vtmp_lap1_pl(g,k,l,4) * KH_coef_lap1_h_pl(g,k,l) ) * rhog_h_pl(g,k,l)
          enddo

          do g = 1, ADM_gall_pl
             tendency_pl(g,k,l,I_RHOGE   ) = - ( vtmp_pl(g,k,l,5) + vtmp_lap1_pl(g,k,l,5) )
             tendency_pl(g,k,l,I_RHOGETOT) = - ( vtmp_pl(g,k,l,5) + vtmp_lap1_pl(g,k,l,5) )
             tendency_pl(g,k,l,I_RHOG    ) = - ( vtmp_pl(g,k,l,6) + vtmp_lap1_pl(g,k,l,6) )
          enddo
       enddo
       enddo
    else
       tendency_pl(:,:,:,:) = 0.0_RP
    endif

    call OPRT_horizontalize_vec( tendency(:,:,:,I_RHOGVX), tendency_pl(:,:,:,I_RHOGVX), & !--- [INOUT]
                                 tendency(:,:,:,I_RHOGVY), tendency_pl(:,:,:,I_RHOGVY), & !--- [INOUT]
                                 tendency(:,:,:,I_RHOGVZ), tendency_pl(:,:,:,I_RHOGVZ)  ) !--- [INOUT]

    !---------------------------------------------------------------------------
    ! For tracer
    !---------------------------------------------------------------------------
    ! 08/04/12 [Mod] T.Mitsui, hyper diffusion is needless for tracer if MIURA2004
    !                          because that is upwind-type advection(already diffusive)
    if ( TRC_ADV_TYPE /= 'MIURA2004' ) then

       qtmp   (:,:,:,:) = q   (:,:,:,:)
       qtmp_pl(:,:,:,:) = q_pl(:,:,:,:)

       ! copy beforehand
       if ( NUMFILTER_DOhorizontaldiff_lap1 ) then
          qtmp_lap1   (:,:,:,:) = qtmp   (:,:,:,:)
          qtmp_lap1_pl(:,:,:,:) = qtmp_pl(:,:,:,:)
       endif

       ! high order laplacian filter
       do p = 1, lap_order_hdiff

          if ( p == lap_order_hdiff ) then

             wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_q * KH_coef   (:,:,:)
             wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_q * KH_coef_pl(:,:,:)

             do nq = 1, TRC_VMAX
                call OPRT_diffusion( qtmp2(:,:,:,nq), qtmp2_pl(:,:,:,nq), & ! [OUT]
                                     qtmp (:,:,:,nq), qtmp_pl (:,:,:,nq), & ! [IN]
                                     wk   (:,:,:),    wk_pl   (:,:,:)     ) ! [IN]
             enddo
          else
             do nq = 1, TRC_VMAX
                call OPRT_laplacian( qtmp2(:,:,:,nq), qtmp2_pl(:,:,:,nq), & ! [OUT]
                                     qtmp (:,:,:,nq), qtmp_pl (:,:,:,nq)  ) ! [IN]
             enddo
          endif

          qtmp   (:,:,:,:) = -qtmp2   (:,:,:,:)
          qtmp_pl(:,:,:,:) = -qtmp2_pl(:,:,:,:)

          call COMM_data_transfer( qtmp, qtmp_pl )

       enddo ! laplacian order loop

       !--- 1st order laplacian filter
       if ( NUMFILTER_DOhorizontaldiff_lap1 ) then

          wk   (:,:,:) = rhog   (:,:,:) * hdiff_fact_q * KH_coef_lap1   (:,:,:)
          wk_pl(:,:,:) = rhog_pl(:,:,:) * hdiff_fact_q * KH_coef_lap1_pl(:,:,:)

          do nq = 1, TRC_VMAX
             call OPRT_diffusion( qtmp2    (:,:,:,nq), qtmp2_pl    (:,:,:,nq), & ! [OUT]
                                  qtmp_lap1(:,:,:,nq), qtmp_lap1_pl(:,:,:,nq), & ! [IN]
                                  wk       (:,:,:),    wk_pl       (:,:,:)     ) ! [IN]
          enddo

          qtmp_lap1   (:,:,:,:) = -qtmp2   (:,:,:,:)
          qtmp_lap1_pl(:,:,:,:) = -qtmp2_pl(:,:,:,:)

          call COMM_data_transfer( qtmp_lap1(:,:,:,:), qtmp_lap1_pl(:,:,:,:) )
       else
          qtmp_lap1   (:,:,:,:) = 0.0_RP
          qtmp_lap1_pl(:,:,:,:) = 0.0_RP
       endif

       do nq = 1, TRC_VMAX
       do l  = 1, ADM_lall
       do k  = ADM_kmin, ADM_kmax
          tendency_q(:,k,l,nq) = - ( qtmp(:,k,l,nq) + qtmp_lap1(:,k,l,nq) )
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do nq = 1, TRC_VMAX
          do l  = 1, ADM_lall_pl
          do k  = ADM_kmin, ADM_kmax
             tendency_q_pl(:,k,l,nq) = - ( qtmp_pl(:,k,l,nq) + qtmp_lap1_pl(:,k,l,nq) )
          enddo
          enddo
          enddo
       else
          tendency_q_pl(:,:,:,:) = 0.0_RP
       endif
    else
       tendency_q   (:,:,:,:) = 0.0_RP
       tendency_q_pl(:,:,:,:) = 0.0_RP
    endif ! apply filter to tracer?

    call DEBUG_rapend('____numfilter_hdiffusion')

    return
  end subroutine numfilter_hdiffusion

  !-----------------------------------------------------------------------------
  !> vertical numerical diffusion
  subroutine numfilter_vdiffusion( &
       rhog,       rhog_pl,      &
       rho,        rho_pl,       &
       vx,         vx_pl,        &
       vy,         vy_pl,        &
       vz,         vz_pl,        &
       w,          w_pl,         &
       tem,        tem_pl,       &
       q,          q_pl,         &
       tendency,   tendency_pl,  &
       tendency_q, tendency_q_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       CVdry => CNST_CV
    use mod_grd, only: &
       GRD_rdgz,  &
       GRD_rdgzh
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_GSGAM2H,    &
       VMTR_GSGAM2H_pl, &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    use mod_runconf, only: &
       TRC_VMAX
    use mod_bsstate, only: &
       rho_bs,    &
       rho_bs_pl, &
       tem_bs,    &
       tem_bs_pl
    implicit none

    real(RP), intent(in)    :: rhog         (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: rhog_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rho          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: rho_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vx           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vx_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vy           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vy_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vz           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vz_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: w            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: w_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: tem          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: tem_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: q            (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP), intent(in)    :: q_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)
    real(RP), intent(inout) :: tendency     (ADM_gall,   ADM_kall,ADM_lall   ,7)
    real(RP), intent(inout) :: tendency_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl,7)
    real(RP), intent(inout) :: tendency_q   (ADM_gall,   ADM_kall,ADM_lall   ,TRC_VMAX)
    real(RP), intent(inout) :: tendency_q_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    integer, parameter :: vmax  = 6
    integer, parameter :: I_VX  = 1
    integer, parameter :: I_VY  = 2
    integer, parameter :: I_VZ  = 3
    integer, parameter :: I_W   = 4
    integer, parameter :: I_TEM = 5
    integer, parameter :: I_RHO = 6

    real(RP) :: rhog_h   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_h_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: flux    (ADM_gall   ,ADM_kall,ADM_lall   ,vmax+TRC_VMAX)
    real(RP) :: flux_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax+TRC_VMAX)
    real(RP) :: vtmp0   (ADM_gall   ,ADM_kall,ADM_lall   ,vmax+TRC_VMAX)
    real(RP) :: vtmp0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax+TRC_VMAX)
    real(RP) :: vtmp1   (ADM_gall   ,ADM_kall,ADM_lall   ,vmax+TRC_VMAX)
    real(RP) :: vtmp1_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,vmax+TRC_VMAX)

    integer :: g, k, l, nq, p
    !---------------------------------------------------------------------------

    if( .NOT. NUMFILTER_DOverticaldiff ) return

    call DEBUG_rapstart('____numfilter_vdiffusion')

    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          rhog_h(g,k,l) = ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
                          + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
       enddo
       enddo
       rhog_h(:,ADM_kmin-1,l) = rhog_h(:,ADM_kmin,l)
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
                                + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
          enddo
          enddo
          rhog_h_pl(:,ADM_kmin-1,l) = rhog_h_pl(:,ADM_kmin,l)
       enddo
    endif

    vtmp0(:,:,:,I_VX ) = vx (:,:,:)
    vtmp0(:,:,:,I_VY ) = vy (:,:,:)
    vtmp0(:,:,:,I_VZ ) = vz (:,:,:)
    vtmp0(:,:,:,I_W  ) = w  (:,:,:)
    vtmp0(:,:,:,I_TEM) = tem(:,:,:) - tem_bs(:,:,:)
    vtmp0(:,:,:,I_RHO) = rho(:,:,:) - rho_bs(:,:,:)
    do nq = 1, TRC_VMAX
       vtmp0(:,:,:,vmax+nq) = rho(:,:,:) * q(:,:,:,nq)
    enddo

    !--- bottom boundary
    vtmp0(:,ADM_kmin-1,:,I_VX ) = vtmp0(:,ADM_kmin,:,I_VX)
    vtmp0(:,ADM_kmin-1,:,I_VY ) = vtmp0(:,ADM_kmin,:,I_VY)
    vtmp0(:,ADM_kmin-1,:,I_VZ ) = vtmp0(:,ADM_kmin,:,I_VZ)
    vtmp0(:,ADM_kmin  ,:,I_W  ) = 0.0_RP
    vtmp0(:,ADM_kmin-1,:,I_TEM) = 3.0_RP * vtmp0(:,ADM_kmin  ,:,I_TEM) &
                                - 3.0_RP * vtmp0(:,ADM_kmin+1,:,I_TEM) &
                                + 1.0_RP * vtmp0(:,ADM_kmin+2,:,I_TEM)
    vtmp0(:,ADM_kmin-1,:,I_RHO) = 3.0_RP * vtmp0(:,ADM_kmin  ,:,I_RHO) &
                                - 3.0_RP * vtmp0(:,ADM_kmin+1,:,I_RHO) &
                                + 1.0_RP * vtmp0(:,ADM_kmin+2,:,I_RHO)

    !--- top boundary
    vtmp0(:,ADM_kmax+1,:,I_VX ) = vtmp0(:,ADM_kmax,:,I_VX)
    vtmp0(:,ADM_kmax+1,:,I_VY ) = vtmp0(:,ADM_kmax,:,I_VY)
    vtmp0(:,ADM_kmax+1,:,I_VZ ) = vtmp0(:,ADM_kmax,:,I_VZ)
    vtmp0(:,ADM_kmax+1,:,I_W  ) = 0.0_RP
    vtmp0(:,ADM_kmax+1,:,I_TEM) = 3.0_RP * vtmp0(:,ADM_kmax  ,:,I_TEM) &
                                - 3.0_RP * vtmp0(:,ADM_kmax-1,:,I_TEM) &
                                + 1.0_RP * vtmp0(:,ADM_kmax-2,:,I_TEM)
    vtmp0(:,ADM_kmax+1,:,I_RHO) = 3.0_RP * vtmp0(:,ADM_kmax  ,:,I_RHO) &
                                - 3.0_RP * vtmp0(:,ADM_kmax-1,:,I_RHO) &
                                + 1.0_RP * vtmp0(:,ADM_kmax-2,:,I_RHO)

    do nq = 1, TRC_VMAX
       vtmp0(:,ADM_kmin-1,:,vmax+nq) = 3.0_RP * vtmp0(:,ADM_kmin  ,:,vmax+nq) &
                                     - 3.0_RP * vtmp0(:,ADM_kmin+1,:,vmax+nq) &
                                     + 1.0_RP * vtmp0(:,ADM_kmin+2,:,vmax+nq)
       vtmp0(:,ADM_kmax+1,:,vmax+nq) = 3.0_RP * vtmp0(:,ADM_kmax  ,:,vmax+nq) &
                                     - 3.0_RP * vtmp0(:,ADM_kmax-1,:,vmax+nq) &
                                     + 1.0_RP * vtmp0(:,ADM_kmax-2,:,vmax+nq)
    enddo

    do l = 1, ADM_lall
       do p = 1, 2
          do k = ADM_kmin, ADM_kmax
             vtmp1(:,k,l,I_VX ) = ( ( vtmp0(:,k+1,l,I_VX ) - vtmp0(:,k  ,l,I_VX ) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_VX ) - vtmp0(:,k-1,l,I_VX ) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp1(:,k,l,I_VY ) = ( ( vtmp0(:,k+1,l,I_VY ) - vtmp0(:,k  ,l,I_VY ) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_VY ) - vtmp0(:,k-1,l,I_VY ) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp1(:,k,l,I_VZ ) = ( ( vtmp0(:,k+1,l,I_VZ ) - vtmp0(:,k  ,l,I_VZ ) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_VZ ) - vtmp0(:,k-1,l,I_VZ ) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp1(:,k,l,I_TEM) = ( ( vtmp0(:,k+1,l,I_TEM) - vtmp0(:,k  ,l,I_TEM) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_TEM) - vtmp0(:,k-1,l,I_TEM) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             vtmp1(:,k,l,I_RHO) = ( ( vtmp0(:,k+1,l,I_RHO) - vtmp0(:,k  ,l,I_RHO) ) * GRD_rdgzh(k+1) &
                                  - ( vtmp0(:,k  ,l,I_RHO) - vtmp0(:,k-1,l,I_RHO) ) * GRD_rdgzh(k)   &
                                  ) * GRD_rdgz(k)
             do nq = 1, TRC_VMAX
                vtmp1(:,k,l,vmax+nq) = ( ( vtmp0(:,k+1,l,vmax+nq) - vtmp0(:,k  ,l,vmax+nq) ) * GRD_rdgzh(k+1) &
                                       - ( vtmp0(:,k  ,l,vmax+nq) - vtmp0(:,k-1,l,vmax+nq) ) * GRD_rdgzh(k)   &
                                       )  * GRD_rdgz(k)
             enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
             vtmp1(:,k,l,I_W) = ( ( vtmp0(:,k+1,l,I_W) - vtmp0(:,k  ,l,I_W) ) * GRD_rdgz(k)   &
                                - ( vtmp0(:,k  ,l,I_W) - vtmp0(:,k-1,l,I_W) ) * GRD_rdgz(k-1) &
                                ) * GRD_rdgzh(k)
          enddo

          if ( p == 1 ) then
             !--- bottom boundary
             vtmp1(:,ADM_kmin-1,l,I_VX ) = vtmp1(:,ADM_kmin  ,l,I_VX )
             vtmp1(:,ADM_kmin-1,l,I_VY ) = vtmp1(:,ADM_kmin  ,l,I_VY )
             vtmp1(:,ADM_kmin-1,l,I_VZ ) = vtmp1(:,ADM_kmin  ,l,I_VZ )
             vtmp1(:,ADM_kmin  ,l,I_W  ) = vtmp1(:,ADM_kmin+1,l,I_W  )
             vtmp1(:,ADM_kmin-1,l,I_TEM) = vtmp1(:,ADM_kmin  ,l,I_TEM) * 2.0_RP - vtmp1(:,ADM_kmin+1,l,I_TEM)
             vtmp1(:,ADM_kmin-1,l,I_RHO) = vtmp1(:,ADM_kmin  ,l,I_RHO) * 2.0_RP - vtmp1(:,ADM_kmin+1,l,I_RHO)

             !--- top boundary
             vtmp1(:,ADM_kmax+1,l,I_VX ) = vtmp1(:,ADM_kmax,l,I_VX )
             vtmp1(:,ADM_kmax+1,l,I_VY ) = vtmp1(:,ADM_kmax,l,I_VY )
             vtmp1(:,ADM_kmax+1,l,I_VZ ) = vtmp1(:,ADM_kmax,l,I_VZ )
             vtmp1(:,ADM_kmax+1,l,I_W  ) = vtmp1(:,ADM_kmax,l,I_W  )
             vtmp1(:,ADM_kmax+1,l,I_TEM) = vtmp1(:,ADM_kmax,l,I_TEM) * 2.0_RP - vtmp1(:,ADM_kmax-1,l,I_TEM)
             vtmp1(:,ADM_kmax+1,l,I_RHO) = vtmp1(:,ADM_kmax,l,I_RHO) * 2.0_RP - vtmp1(:,ADM_kmax-1,l,I_RHO)

             do nq = 1, TRC_VMAX
                vtmp1(:,ADM_kmin-1,l,vmax+nq) = 2.0_RP * vtmp1(:,ADM_kmin,l,vmax+nq) - vtmp1(:,ADM_kmin+1,l,vmax+nq)
                vtmp1(:,ADM_kmax+1,l,vmax+nq) = 2.0_RP * vtmp1(:,ADM_kmax,l,vmax+nq) - vtmp1(:,ADM_kmax-1,l,vmax+nq)
             enddo

             vtmp0(:,:,l,:) = vtmp1(:,:,l,:)
          elseif( p == 2 ) then
             !--- bottom boundary
             vtmp1(:,ADM_kmin-1,l,I_VX ) = vtmp1(:,ADM_kmin  ,l,I_VX )
             vtmp1(:,ADM_kmin-1,l,I_VY ) = vtmp1(:,ADM_kmin  ,l,I_VY )
             vtmp1(:,ADM_kmin-1,l,I_VZ ) = vtmp1(:,ADM_kmin  ,l,I_VZ )
             vtmp1(:,ADM_kmin  ,l,I_W  ) = vtmp1(:,ADM_kmin+1,l,I_W  )
             vtmp1(:,ADM_kmin-1,l,I_TEM) = vtmp1(:,ADM_kmin  ,l,I_TEM)
             vtmp1(:,ADM_kmin-1,l,I_RHO) = vtmp1(:,ADM_kmin  ,l,I_RHO)

             !--- top boundary
             vtmp1(:,ADM_kmax+1,l,I_VX ) = vtmp1(:,ADM_kmax,l,I_VX )
             vtmp1(:,ADM_kmax+1,l,I_VY ) = vtmp1(:,ADM_kmax,l,I_VY )
             vtmp1(:,ADM_kmax+1,l,I_VZ ) = vtmp1(:,ADM_kmax,l,I_VZ )
             vtmp1(:,ADM_kmax+1,l,I_W  ) = vtmp1(:,ADM_kmax,l,I_W  )
             vtmp1(:,ADM_kmax+1,l,I_TEM) = vtmp1(:,ADM_kmax,l,I_TEM)
             vtmp1(:,ADM_kmax+1,l,I_RHO) = vtmp1(:,ADM_kmax,l,I_RHO)

             do nq = 1, TRC_VMAX
                vtmp1(:,ADM_kmin-1,l,vmax+nq) = vtmp1(:,ADM_kmin,l,vmax+nq)
                vtmp1(:,ADM_kmax+1,l,vmax+nq) = vtmp1(:,ADM_kmax,l,vmax+nq)
             enddo

             vtmp0(:,:,l,:) = vtmp1(:,:,l,:)
          endif
       enddo

       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          flux(g,k,l,I_VX ) = Kv_coef_h(k) * ( vtmp0(g,k,l,I_VX )-vtmp0(g,k-1,l,I_VX ) ) * GRD_rdgzh(k) * rhog_h(g,k,l)
          flux(g,k,l,I_VY ) = Kv_coef_h(k) * ( vtmp0(g,k,l,I_VY )-vtmp0(g,k-1,l,I_VY ) ) * GRD_rdgzh(k) * rhog_h(g,k,l)
          flux(g,k,l,I_VZ ) = Kv_coef_h(k) * ( vtmp0(g,k,l,I_VZ )-vtmp0(g,k-1,l,I_VZ ) ) * GRD_rdgzh(k) * rhog_h(g,k,l)
          flux(g,k,l,I_TEM) = Kv_coef_h(k) * ( vtmp0(g,k,l,I_TEM)-vtmp0(g,k-1,l,I_TEM) ) * GRD_rdgzh(k) * rhog_h(g,k,l) * CVdry
          flux(g,k,l,I_RHO) = Kv_coef_h(k) * ( vtmp0(g,k,l,I_RHO)-vtmp0(g,k-1,l,I_RHO) ) * GRD_rdgzh(k) * VMTR_GSGAM2H(g,k,l)
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          flux(g,k,l,I_W) = Kv_coef(k) * ( vtmp0(g,k+1,l,I_W)-vtmp0(g,k,l,I_W) ) * GRD_rdgz(k) * rhog(g,k,l)
       enddo
       enddo

       do nq = 1, TRC_VMAX
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          flux(g,k,l,vmax+nq) = Kv_coef_h(k) * ( vtmp0(g,k,l,vmax+nq) - vtmp0(g,k-1,l,vmax+nq) ) * GRD_rdgzh(k)
       enddo
       enddo
       enddo

       !--- update tendency
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          tendency(g,k,l,I_RHOG    ) = tendency(g,k,l,I_RHOG    ) + ( flux(g,k+1,l,I_RHO) - flux(g,k,l,I_RHO) ) * GRD_rdgz(k)
          tendency(g,k,l,I_RHOGVX  ) = tendency(g,k,l,I_RHOGVX  ) + ( flux(g,k+1,l,I_VX ) - flux(g,k,l,I_VX ) ) * GRD_rdgz(k)
          tendency(g,k,l,I_RHOGVY  ) = tendency(g,k,l,I_RHOGVY  ) + ( flux(g,k+1,l,I_VY ) - flux(g,k,l,I_VY ) ) * GRD_rdgz(k)
          tendency(g,k,l,I_RHOGVZ  ) = tendency(g,k,l,I_RHOGVZ  ) + ( flux(g,k+1,l,I_VZ ) - flux(g,k,l,I_VZ ) ) * GRD_rdgz(k)
          tendency(g,k,l,I_RHOGE   ) = tendency(g,k,l,I_RHOGE   ) + ( flux(g,k+1,l,I_TEM) - flux(g,k,l,I_TEM) ) * GRD_rdgz(k)
          tendency(g,k,l,I_RHOGETOT) = tendency(g,k,l,I_RHOGETOT) + ( flux(g,k+1,l,I_TEM) - flux(g,k,l,I_TEM) ) * GRD_rdgz(k)
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          tendency(g,k,l,I_RHOGW) = tendency(g,k,l,I_RHOGW) + ( flux(g,k,l,I_W) - flux(g,k-1,l,I_W) ) * GRD_rdgzh(k)
       enddo
       enddo

       do nq = 1, TRC_VMAX
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          tendency_q(g,k,l,nq) = tendency_q(g,k,l,nq) + ( flux(g,k+1,l,vmax+nq) - flux(g,k,l,vmax+nq) ) * GRD_rdgz(k)
       enddo
       enddo
       enddo
    enddo

    if ( ADM_have_pl ) then

       vtmp0_pl(:,:,:,I_VX ) = vx_pl (:,:,:)
       vtmp0_pl(:,:,:,I_VY ) = vy_pl (:,:,:)
       vtmp0_pl(:,:,:,I_VZ ) = vz_pl (:,:,:)
       vtmp0_pl(:,:,:,I_W  ) = w_pl  (:,:,:)
       vtmp0_pl(:,:,:,I_TEM) = tem_pl(:,:,:) - tem_bs_pl(:,:,:)
       vtmp0_pl(:,:,:,I_RHO) = rho_pl(:,:,:) - rho_bs_pl(:,:,:)

       do nq = 1, TRC_VMAX
          vtmp0_pl(:,:,:,vmax+nq) = rho_pl(:,:,:) * q_pl(:,:,:,nq)
       enddo

       !--- bottom boundary
       vtmp0_pl(:,ADM_kmin-1,:,I_VX ) = vtmp0_pl(:,ADM_kmin,:,I_VX)
       vtmp0_pl(:,ADM_kmin-1,:,I_VY ) = vtmp0_pl(:,ADM_kmin,:,I_VY)
       vtmp0_pl(:,ADM_kmin-1,:,I_VZ ) = vtmp0_pl(:,ADM_kmin,:,I_VZ)
       vtmp0_pl(:,ADM_kmin  ,:,I_W  ) = 0.0_RP
       vtmp0_pl(:,ADM_kmin-1,:,I_TEM) = 3.0_RP * vtmp0_pl(:,ADM_kmin  ,:,I_TEM) &
                                      - 3.0_RP * vtmp0_pl(:,ADM_kmin+1,:,I_TEM) &
                                      + 1.0_RP * vtmp0_pl(:,ADM_kmin+2,:,I_TEM)
       vtmp0_pl(:,ADM_kmin-1,:,I_RHO) = 3.0_RP * vtmp0_pl(:,ADM_kmin  ,:,I_RHO) &
                                      - 3.0_RP * vtmp0_pl(:,ADM_kmin+1,:,I_RHO) &
                                      + 1.0_RP * vtmp0_pl(:,ADM_kmin+2,:,I_RHO)

       !--- top boundary
       vtmp0_pl(:,ADM_kmax+1,:,I_VX ) = vtmp0_pl(:,ADM_kmax,:,I_VX)
       vtmp0_pl(:,ADM_kmax+1,:,I_VY ) = vtmp0_pl(:,ADM_kmax,:,I_VY)
       vtmp0_pl(:,ADM_kmax+1,:,I_VZ ) = vtmp0_pl(:,ADM_kmax,:,I_VZ)
       vtmp0_pl(:,ADM_kmax+1,:,I_W  ) = 0.0_RP
       vtmp0_pl(:,ADM_kmax+1,:,I_TEM) = 3.0_RP * vtmp0_pl(:,ADM_kmax  ,:,I_TEM) &
                                      - 3.0_RP * vtmp0_pl(:,ADM_kmax-1,:,I_TEM) &
                                      + 1.0_RP * vtmp0_pl(:,ADM_kmax-2,:,I_TEM)
       vtmp0_pl(:,ADM_kmax+1,:,I_RHO) = 3.0_RP * vtmp0_pl(:,ADM_kmax  ,:,I_RHO) &
                                      - 3.0_RP * vtmp0_pl(:,ADM_kmax-1,:,I_RHO) &
                                      + 1.0_RP * vtmp0_pl(:,ADM_kmax-2,:,I_RHO)

       do nq = 1, TRC_VMAX
          vtmp0_pl(:,ADM_kmin-1,:,vmax+nq) = 3.0_RP * vtmp0_pl(:,ADM_kmin  ,:,vmax+nq) &
                                           - 3.0_RP * vtmp0_pl(:,ADM_kmin+1,:,vmax+nq) &
                                           + 1.0_RP * vtmp0_pl(:,ADM_kmin+2,:,vmax+nq)
          vtmp0_pl(:,ADM_kmax+1,:,vmax+nq) = 3.0_RP * vtmp0_pl(:,ADM_kmax  ,:,vmax+nq) &
                                           - 3.0_RP * vtmp0_pl(:,ADM_kmax-1,:,vmax+nq) &
                                           + 1.0_RP * vtmp0_pl(:,ADM_kmax-2,:,vmax+nq)
       enddo

       do l = 1, ADM_lall
          do p = 1, 2
             do k = ADM_kmin, ADM_kmax
                vtmp1_pl(:,k,l,I_VX ) = ( ( vtmp0_pl(:,k+1,l,I_VX )-vtmp0_pl(:,k  ,l,I_VX ) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k  ,l,I_VX )-vtmp0_pl(:,k-1,l,I_VX ) ) * GRD_rdgzh(k  ) &
                                        ) * GRD_rdgz(k)
                vtmp1_pl(:,k,l,I_VY ) = ( ( vtmp0_pl(:,k+1,l,I_VY )-vtmp0_pl(:,k  ,l,I_VY ) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k  ,l,I_VY )-vtmp0_pl(:,k-1,l,I_VY ) ) * GRD_rdgzh(k  ) &
                                        ) * GRD_rdgz(k)
                vtmp1_pl(:,k,l,I_VZ ) = ( ( vtmp0_pl(:,k+1,l,I_VZ )-vtmp0_pl(:,k  ,l,I_VZ ) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k  ,l,I_VZ )-vtmp0_pl(:,k-1,l,I_VZ ) ) * GRD_rdgzh(k  ) &
                                        ) * GRD_rdgz(k)
                vtmp1_pl(:,k,l,I_TEM) = ( ( vtmp0_pl(:,k+1,l,I_TEM)-vtmp0_pl(:,k  ,l,I_TEM) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k  ,l,I_TEM)-vtmp0_pl(:,k-1,l,I_TEM) ) * GRD_rdgzh(k  ) &
                                        ) * GRD_rdgz(k)
                vtmp1_pl(:,k,l,I_RHO) = ( ( vtmp0_pl(:,k+1,l,I_RHO)-vtmp0_pl(:,k  ,l,I_RHO) ) * GRD_rdgzh(k+1) &
                                        - ( vtmp0_pl(:,k  ,l,I_RHO)-vtmp0_pl(:,k-1,l,I_RHO) ) * GRD_rdgzh(k  ) &
                                        ) * GRD_rdgz(k)

                do nq = 1, TRC_VMAX
                   vtmp1_pl(:,k,l,vmax+nq) = ( ( vtmp0_pl(:,k+1,l,vmax+nq)-vtmp0_pl(:,k  ,l,vmax+nq) ) * GRD_rdgzh(k+1) &
                                             - ( vtmp0_pl(:,k  ,l,vmax+nq)-vtmp0_pl(:,k-1,l,vmax+nq) ) * GRD_rdgzh(k  ) &
                                             )  * GRD_rdgz(k)
                enddo
             enddo

             do k = ADM_kmin+1, ADM_kmax
                vtmp1_pl(:,k,l,I_W) = ( ( vtmp0_pl(:,k+1,l,I_W)-vtmp0_pl(:,k  ,l,I_W) ) * GRD_rdgz(k  ) &
                                      - ( vtmp0_pl(:,k  ,l,I_W)-vtmp0_pl(:,k-1,l,I_W) ) * GRD_rdgz(k-1) &
                                      ) * GRD_rdgzh(k)
             enddo

             if ( p == 1 ) then

                !--- bottom boundary
                vtmp1_pl(:,ADM_kmin-1,l,I_VX ) = vtmp1_pl(:,ADM_kmin  ,l,I_VX )
                vtmp1_pl(:,ADM_kmin-1,l,I_VY ) = vtmp1_pl(:,ADM_kmin  ,l,I_VY )
                vtmp1_pl(:,ADM_kmin-1,l,I_VZ ) = vtmp1_pl(:,ADM_kmin  ,l,I_VZ )
                vtmp1_pl(:,ADM_kmin  ,l,I_W  ) = vtmp1_pl(:,ADM_kmin+1,l,I_W  )
                vtmp1_pl(:,ADM_kmin-1,l,I_TEM) = vtmp1_pl(:,ADM_kmin  ,l,I_TEM) * 2.0_RP - vtmp1_pl(:,ADM_kmin+1,l,I_TEM)
                vtmp1_pl(:,ADM_kmin-1,l,I_RHO) = vtmp1_pl(:,ADM_kmin  ,l,I_RHO) * 2.0_RP - vtmp1_pl(:,ADM_kmin+1,l,I_RHO)

                !--- top boundary
                vtmp1_pl(:,ADM_kmax+1,l,I_VX ) = vtmp1_pl(:,ADM_kmax,l,I_VX )
                vtmp1_pl(:,ADM_kmax+1,l,I_VY ) = vtmp1_pl(:,ADM_kmax,l,I_VY )
                vtmp1_pl(:,ADM_kmax+1,l,I_VZ ) = vtmp1_pl(:,ADM_kmax,l,I_VZ )
                vtmp1_pl(:,ADM_kmax+1,l,I_W  ) = vtmp1_pl(:,ADM_kmax,l,I_W  )
                vtmp1_pl(:,ADM_kmax+1,l,I_TEM) = vtmp1_pl(:,ADM_kmax,l,I_TEM) * 2.0_RP - vtmp1_pl(:,ADM_kmax-1,l,I_TEM)
                vtmp1_pl(:,ADM_kmax+1,l,I_RHO) = vtmp1_pl(:,ADM_kmax,l,I_RHO) * 2.0_RP - vtmp1_pl(:,ADM_kmax-1,l,I_RHO)

                do nq = 1, TRC_VMAX
                   vtmp1_pl(:,ADM_kmin-1,l,vmax+nq) = 2.0_RP * vtmp1_pl(:,ADM_kmin  ,l,vmax+nq) &
                                                    - 1.0_RP * vtmp1_pl(:,ADM_kmin+1,l,vmax+nq)
                   vtmp1_pl(:,ADM_kmax+1,l,vmax+nq) = 2.0_RP * vtmp1_pl(:,ADM_kmax  ,l,vmax+nq) &
                                                    - 1.0_RP * vtmp1_pl(:,ADM_kmax-1,l,vmax+nq)
                enddo

             elseif( p == 2 ) then

                vtmp1_pl(:,ADM_kmin-1,l,I_VX ) = vtmp1_pl(:,ADM_kmin  ,l,I_VX )
                vtmp1_pl(:,ADM_kmin-1,l,I_VY ) = vtmp1_pl(:,ADM_kmin  ,l,I_VY )
                vtmp1_pl(:,ADM_kmin-1,l,I_VZ ) = vtmp1_pl(:,ADM_kmin  ,l,I_VZ )
                vtmp1_pl(:,ADM_kmin  ,l,I_W  ) = vtmp1_pl(:,ADM_kmin+1,l,I_W  )
                vtmp1_pl(:,ADM_kmin-1,l,I_TEM) = vtmp1_pl(:,ADM_kmin  ,l,I_TEM)
                vtmp1_pl(:,ADM_kmin-1,l,I_RHO) = vtmp1_pl(:,ADM_kmin  ,l,I_RHO)

                vtmp1_pl(:,ADM_kmax+1,l,I_VX ) = vtmp1_pl(:,ADM_kmax,l,I_VX )
                vtmp1_pl(:,ADM_kmax+1,l,I_VY ) = vtmp1_pl(:,ADM_kmax,l,I_VY )
                vtmp1_pl(:,ADM_kmax+1,l,I_VZ ) = vtmp1_pl(:,ADM_kmax,l,I_VZ )
                vtmp1_pl(:,ADM_kmax+1,l,I_W  ) = vtmp1_pl(:,ADM_kmax,l,I_W  )
                vtmp1_pl(:,ADM_kmax+1,l,I_TEM) = vtmp1_pl(:,ADM_kmax,l,I_TEM)
                vtmp1_pl(:,ADM_kmax+1,l,I_RHO) = vtmp1_pl(:,ADM_kmax,l,I_RHO)

                do nq = 1, TRC_VMAX
                   vtmp1_pl(:,ADM_kmin-1,l,vmax+nq) = vtmp1_pl(:,ADM_kmin,l,vmax+nq)
                   vtmp1_pl(:,ADM_kmax+1,l,vmax+nq) = vtmp1_pl(:,ADM_kmax,l,vmax+nq)
                enddo

             endif
          enddo
          vtmp0_pl(:,:,:,:) = vtmp1_pl(:,:,:,:)

          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall
             flux_pl(g,k,l,I_VX ) = Kv_coef_h(k) * ( vtmp0_pl(g,k,l,I_VX )-vtmp0_pl(g,k-1,l,I_VX ) ) &
                                  * GRD_rdgzh(k) * rhog_h_pl(g,k,l)
             flux_pl(g,k,l,I_VY ) = Kv_coef_h(k) * ( vtmp0_pl(g,k,l,I_VY )-vtmp0_pl(g,k-1,l,I_VY ) ) &
                                  * GRD_rdgzh(k) * rhog_h_pl(g,k,l)
             flux_pl(g,k,l,I_VZ ) = Kv_coef_h(k) * ( vtmp0_pl(g,k,l,I_VZ )-vtmp0_pl(g,k-1,l,I_VZ ) ) &
                                  * GRD_rdgzh(k) * rhog_h_pl(g,k,l)
             flux_pl(g,k,l,I_TEM) = Kv_coef_h(k) * ( vtmp0_pl(g,k,l,I_TEM)-vtmp0_pl(g,k-1,l,I_TEM) ) &
                                  * GRD_rdgzh(k) * rhog_h_pl(g,k,l) * CVdry
             flux_pl(g,k,l,I_RHO) = Kv_coef_h(k) * ( vtmp0_pl(g,k,l,I_RHO)-vtmp0_pl(g,k-1,l,I_RHO) ) &
                                  * GRD_rdgzh(k) * VMTR_GSGAM2H_pl(g,k,l)
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             flux_pl(g,k,l,I_W) = Kv_coef(k) * ( vtmp0_pl(g,k+1,l,I_W)-vtmp0_pl(g,k,l,I_W) ) &
                                * GRD_rdgz(k) * rhog_pl(g,k,l)
          enddo
          enddo

          do nq = 1, TRC_VMAX
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall
             flux_pl(g,k,l,vmax+nq) = Kv_coef_h(k) * ( vtmp0_pl(g,k,l,vmax+nq)-vtmp0_pl(g,k-1,l,vmax+nq) ) &
                                    * GRD_rdgzh(k) * rhog_h_pl(g,k,l)
          enddo
          enddo
          enddo

          !--- update tendency
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             tendency_pl(g,k,l,I_RHOG    ) = tendency_pl(g,k,l,I_RHOG    ) &
                                           + ( flux_pl(g,k+1,l,I_RHO) - flux_pl(g,k,l,I_RHO) ) * GRD_rdgz(k)
             tendency_pl(g,k,l,I_RHOGVX  ) = tendency_pl(g,k,l,I_RHOGVX  ) &
                                           + ( flux_pl(g,k+1,l,I_VX ) - flux_pl(g,k,l,I_VX ) ) * GRD_rdgz(k)
             tendency_pl(g,k,l,I_RHOGVY  ) = tendency_pl(g,k,l,I_RHOGVY  ) &
                                           + ( flux_pl(g,k+1,l,I_VY ) - flux_pl(g,k,l,I_VY ) ) * GRD_rdgz(k)
             tendency_pl(g,k,l,I_RHOGVZ  ) = tendency_pl(g,k,l,I_RHOGVZ  ) &
                                           + ( flux_pl(g,k+1,l,I_VZ ) - flux_pl(g,k,l,I_VZ ) ) * GRD_rdgz(k)
             tendency_pl(g,k,l,I_RHOGE   ) = tendency_pl(g,k,l,I_RHOGE   ) &
                                           + ( flux_pl(g,k+1,l,I_TEM) - flux_pl(g,k,l,I_TEM) ) * GRD_rdgz(k)
             tendency_pl(g,k,l,I_RHOGETOT) = tendency_pl(g,k,l,I_RHOGETOT) &
                                           + ( flux_pl(g,k+1,l,I_TEM) - flux_pl(g,k,l,I_TEM) ) * GRD_rdgz(k)
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall
             tendency_pl(g,k,l,I_RHOGW) = tendency_pl(g,k,l,I_RHOGW) &
                                        + ( flux_pl(g,k,l,I_W) - flux_pl(g,k-1,l,I_W) ) * GRD_rdgzh(k)
          enddo
          enddo

          do nq = 1, TRC_VMAX
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall
             tendency_q_pl(g,k,l,nq) = tendency_q_pl(g,k,l,nq) &
                                     + ( flux_pl(g,k+1,l,vmax+nq) - flux_pl(g,k,l,vmax+nq) ) * GRD_rdgz(k)
          enddo
          enddo
          enddo

       enddo
    endif

    call OPRT_horizontalize_vec( tendency(:,:,:,I_RHOGVX), tendency_pl(:,:,:,I_RHOGVX), & !--- [INOUT]
                                 tendency(:,:,:,I_RHOGVY), tendency_pl(:,:,:,I_RHOGVY), & !--- [INOUT]
                                 tendency(:,:,:,I_RHOGVZ), tendency_pl(:,:,:,I_RHOGVZ)  ) !--- [INOUT]

    call DEBUG_rapend('____numfilter_vdiffusion')

    return
  end subroutine numfilter_vdiffusion

  !-----------------------------------------------------------------------------
  !> 3D divergence damping
  subroutine numfilter_divdamp( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       gdx,    gdx_pl,    &
       gdy,    gdy_pl,    &
       gdz,    gdz_pl,    &
       gdvz,   gdvz_pl    )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only:  &
       GRD_rdgzh
    use mod_oprt, only: &
       OPRT_horizontalize_vec, &
       OPRT_divdamp
    use mod_oprt3d, only: &
       OPRT3D_divdamp
    use mod_src, only: &
       src_flux_convergence, &
       I_SRC_default
    implicit none

    real(RP), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vx ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*Vy ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw    (ADM_gall,   ADM_kall,ADM_lall   ) ! rho*w  ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: gdx      (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( G^1/2 x gam2 )
    real(RP), intent(out) :: gdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: gdy      (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( G^1/2 x gam2 )
    real(RP), intent(out) :: gdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: gdz      (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( G^1/2 x gam2 )
    real(RP), intent(out) :: gdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: gdvz     (ADM_gall,   ADM_kall,ADM_lall   ) ! (grad div)_x ( G^1/2 x gam2 )
    real(RP), intent(out) :: gdvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: vtmp    (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(RP) :: vtmp_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(RP) :: vtmp2   (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(RP) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    real(RP) :: cnv     (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: cnv_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: k, l, p
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____numfilter_divdamp')

    if ( .NOT. NUMFILTER_DOdivdamp ) then
       gdx    (:,:,:) = 0.0_RP
       gdx_pl (:,:,:) = 0.0_RP
       gdy    (:,:,:) = 0.0_RP
       gdy_pl (:,:,:) = 0.0_RP
       gdz    (:,:,:) = 0.0_RP
       gdz_pl (:,:,:) = 0.0_RP
       gdvz   (:,:,:) = 0.0_RP
       gdvz_pl(:,:,:) = 0.0_RP
       call DEBUG_rapend('____numfilter_divdamp')
       return
    endif

    !--- 3D divergence divdamp
    call OPRT3D_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & ! [OUT]
                         vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & ! [OUT]
                         vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & ! [OUT]
                         rhogvx(:,:,:),  rhogvx_pl(:,:,:),  & ! [IN]
                         rhogvy(:,:,:),  rhogvy_pl(:,:,:),  & ! [IN]
                         rhogvz(:,:,:),  rhogvz_pl(:,:,:),  & ! [IN]
                         rhogw (:,:,:),  rhogw_pl (:,:,:)   ) ! [IN]

    if ( lap_order_divdamp > 1 ) then
       do p = 1, lap_order_divdamp-1

          call COMM_data_transfer( vtmp2, vtmp2_pl )

          !--- note : sign changes
          vtmp   (:,:,:,:) = -vtmp2   (:,:,:,:)
          vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

          !--- 2D dinvergence divdamp
          call OPRT_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & ! [OUT]
                             vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & ! [OUT]
                             vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & ! [OUT]
                             vtmp (:,:,:,1), vtmp_pl (:,:,:,1), & ! [IN]
                             vtmp (:,:,:,2), vtmp_pl (:,:,:,2), & ! [IN]
                             vtmp (:,:,:,3), vtmp_pl (:,:,:,3)  ) ! [IN]
       enddo ! lap_order
    endif

    !--- X coeffcient
    gdx(:,:,:) = divdamp_coef(:,:,:) * vtmp2(:,:,:,1)
    gdy(:,:,:) = divdamp_coef(:,:,:) * vtmp2(:,:,:,2)
    gdz(:,:,:) = divdamp_coef(:,:,:) * vtmp2(:,:,:,3)

    if ( ADM_have_pl ) then
       gdx_pl(:,:,:) = divdamp_coef_pl(:,:,:) * vtmp2_pl(:,:,:,1)
       gdy_pl(:,:,:) = divdamp_coef_pl(:,:,:) * vtmp2_pl(:,:,:,2)
       gdz_pl(:,:,:) = divdamp_coef_pl(:,:,:) * vtmp2_pl(:,:,:,3)
    endif

    call OPRT_horizontalize_vec( gdx(:,:,:), gdx_pl(:,:,:), & !--- [INOUT]
                                 gdy(:,:,:), gdy_pl(:,:,:), & !--- [INOUT]
                                 gdz(:,:,:), gdz_pl(:,:,:)  ) !--- [INOUT]

    if ( NUMFILTER_DOdivdamp_v ) then

       call src_flux_convergence( rhogvx(:,:,:), rhogvx_pl(:,:,:), & ! [IN]
                                  rhogvy(:,:,:), rhogvy_pl(:,:,:), & ! [IN]
                                  rhogvz(:,:,:), rhogvz_pl(:,:,:), & ! [IN]
                                  rhogw (:,:,:), rhogw_pl (:,:,:), & ! [IN]
                                  cnv   (:,:,:), cnv_pl   (:,:,:), & ! [OUT]
                                  I_SRC_default                    ) ! [IN]

       do l = 1, ADM_lall
          do k = ADM_kmin+1, ADM_kmax
             gdvz(:,k,l) = divdamp_coef_v * ( cnv(:,k,l) - cnv(:,k-1,l) ) * GRD_rdgzh(k)
          enddo
          gdvz(:,ADM_kmin-1,l) = 0.0_RP
          gdvz(:,ADM_kmin  ,l) = 0.0_RP
          gdvz(:,ADM_kmax+1,l) = 0.0_RP
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin+1, ADM_kmax
                gdvz_pl(:,k,l) = divdamp_coef_v * ( cnv_pl(:,k,l) - cnv_pl(:,k-1,l) ) * GRD_rdgzh(k)
             enddo
             gdvz_pl(:,ADM_kmin-1,l) = 0.0_RP
             gdvz_pl(:,ADM_kmin  ,l) = 0.0_RP
             gdvz_pl(:,ADM_kmax+1,l) = 0.0_RP
          enddo
       endif

    else
       gdvz   (:,:,:) = 0.0_RP
       gdvz_pl(:,:,:) = 0.0_RP
    endif

    call DEBUG_rapend('____numfilter_divdamp')

    return
  end subroutine numfilter_divdamp

  !-----------------------------------------------------------------------------
  !> 2D dinvergence divdamp
  subroutine numfilter_divdamp_2d( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       gdx,    gdx_pl,    &
       gdy,    gdy_pl,    &
       gdz,    gdz_pl     )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_comm, only: &
       COMM_data_transfer
    use mod_oprt, only: &
       OPRT_horizontalize_vec, &
       OPRT_divdamp
    implicit none

    real(RP), intent(in)  :: rhogvx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(out) :: gdx      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: gdx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: gdy      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: gdy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: gdz      (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP), intent(out) :: gdz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: vtmp    (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(RP) :: vtmp_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(RP) :: vtmp2   (ADM_gall,   ADM_kall,ADM_lall   ,3)
    real(RP) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    integer :: p
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____numfilter_divdamp_2d')

    if ( .NOT. NUMFILTER_DOdivdamp_2d ) then
       gdx   (:,:,:) = 0.0_RP
       gdx_pl(:,:,:) = 0.0_RP
       gdy   (:,:,:) = 0.0_RP
       gdy_pl(:,:,:) = 0.0_RP
       gdz   (:,:,:) = 0.0_RP
       gdz_pl(:,:,:) = 0.0_RP
       call DEBUG_rapend('____numfilter_divdamp_2d')
       return
    endif

    !--- 2D dinvergence divdamp
    call OPRT_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & ! [OUT]
                       vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & ! [OUT]
                       vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & ! [OUT]
                       rhogvx(:,:,:),  rhogvx_pl(:,:,:),  & ! [IN]
                       rhogvy(:,:,:),  rhogvy_pl(:,:,:),  & ! [IN]
                       rhogvz(:,:,:),  rhogvz_pl(:,:,:)   ) ! [IN]

    if ( lap_order_divdamp_2d > 1 ) then
       do p = 1, lap_order_divdamp_2d-1

          call COMM_data_transfer(vtmp2,vtmp2_pl)

          !--- note : sign changes
          vtmp   (:,:,:,:) = -vtmp2   (:,:,:,:)
          vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

          !--- 2D dinvergence divdamp
          call OPRT_divdamp( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & ! [OUT]
                             vtmp2(:,:,:,2), vtmp2_pl(:,:,:,2), & ! [OUT]
                             vtmp2(:,:,:,3), vtmp2_pl(:,:,:,3), & ! [OUT]
                             vtmp (:,:,:,1), vtmp_pl (:,:,:,1), & ! [IN]
                             vtmp (:,:,:,2), vtmp_pl (:,:,:,2), & ! [IN]
                             vtmp (:,:,:,3), vtmp_pl (:,:,:,3)  ) ! [IN]

       enddo ! lap_order
    endif

    !--- X coeffcient
    gdx(:,:,:) = divdamp_2d_coef(:,:,:) * vtmp2(:,:,:,1)
    gdy(:,:,:) = divdamp_2d_coef(:,:,:) * vtmp2(:,:,:,2)
    gdz(:,:,:) = divdamp_2d_coef(:,:,:) * vtmp2(:,:,:,3)

    if ( ADM_have_pl ) then
       gdx_pl(:,:,:) = divdamp_2d_coef_pl(:,:,:) * vtmp2_pl(:,:,:,1)
       gdy_pl(:,:,:) = divdamp_2d_coef_pl(:,:,:) * vtmp2_pl(:,:,:,2)
       gdz_pl(:,:,:) = divdamp_2d_coef_pl(:,:,:) * vtmp2_pl(:,:,:,3)
    endif

    call OPRT_horizontalize_vec( gdx(:,:,:), gdx_pl(:,:,:), & !--- [INOUT]
                                 gdy(:,:,:), gdy_pl(:,:,:), & !--- [INOUT]
                                 gdz(:,:,:), gdz_pl(:,:,:)  ) !--- [INOUT]

    call DEBUG_rapend('____numfilter_divdamp_2d')

    return
  end subroutine numfilter_divdamp_2d

  !-----------------------------------------------------------------------------
  !> smoothing
  subroutine numfilter_smooth_1var( &
       s, s_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall
    use mod_comm, only: &
       COMM_data_transfer
    use mod_gmtr, only: &
       GMTR_area,    &
       GMTR_area_pl
    use mod_oprt, only: &
       OPRT_laplacian
    implicit none

    real(RP), intent(inout) :: s   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: s_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: vtmp    (ADM_gall   ,ADM_kall,ADM_lall   ,1)
    real(RP) :: vtmp_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,1)
    real(RP) :: vtmp2   (ADM_gall   ,ADM_kall,ADM_lall   ,1)
    real(RP) :: vtmp2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,1)

    real(RP), parameter :: ggamma_h = 1.0_RP / 16.0_RP / 10.0_RP
    integer, parameter :: itelim = 80

    integer :: p, ite
    integer :: k, l
    !---------------------------------------------------------------------------

    do ite = 1, itelim

       vtmp(:,:,:,1) = s(:,:,:)

       vtmp_pl(:,:,:,:) = 0.0_RP

       if ( ADM_have_pl ) then
          vtmp_pl(:,:,:,1) = s_pl(:,:,:)
       endif

       call COMM_data_transfer( vtmp, vtmp_pl )

       do p = 1, 2
          vtmp2   (:,:,:,:) = 0.0_RP
          vtmp2_pl(:,:,:,:) = 0.0_RP

          call OPRT_laplacian( vtmp2(:,:,:,1), vtmp2_pl(:,:,:,1), & ! [OUT]
                               vtmp (:,:,:,1), vtmp_pl (:,:,:,1)  ) ! [IN]

          vtmp   (:,:,:,:) = -vtmp2   (:,:,:,:)
          vtmp_pl(:,:,:,:) = -vtmp2_pl(:,:,:,:)

          call COMM_data_transfer( vtmp, vtmp_pl )
       enddo

       do l = 1, ADM_lall
       do k = 1, ADM_kall
          s(:,k,l) = s(:,k,l) - ggamma_h * GMTR_area(:,l)**2 * vtmp(:,k,l,1)
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
             s_pl(:,k,l) = s_pl(:,k,l) - ggamma_h * GMTR_area_pl(:,l)**2 * vtmp_pl(:,k,l,1)
          enddo
          enddo
       endif
    enddo

    vtmp   (:,:,:,1) = s   (:,:,:)
    vtmp_pl(:,:,:,1) = s_pl(:,:,:)

    call COMM_data_transfer( vtmp, vtmp_pl )

    s   (:,:,:) = vtmp   (:,:,:,1)
    s_pl(:,:,:) = vtmp_pl(:,:,:,1)

    return
  end subroutine numfilter_smooth_1var

  !-----------------------------------------------------------------------------
  !> calc height factor
  subroutine height_factor( &
       kdim,          &
       z,             &
       z_top,         &
       z_bottomlimit, &
       factor         )
    use mod_cnst, only: &
       PI => CNST_PI
    implicit none

    integer, intent(in)  :: kdim          ! number of vertical grid
    real(RP), intent(in)  :: z(kdim)       ! height [m]
    real(RP), intent(in)  :: z_top         ! height top [m]
    real(RP), intent(in)  :: z_bottomlimit ! bottom limit of the factor [m]
    real(RP), intent(out) :: factor(kdim)  ! height-dependent factor [0-1]

    real(RP) :: sw

    integer :: k
    !---------------------------------------------------------------------------

    do k = 1, kdim
       sw = 0.5_RP + sign( 0.5_RP, z(k)-z_bottomlimit )

       factor(k) = sw * 0.5_RP * ( 1.0_RP - cos( PI * (z(k)-z_bottomlimit) / (z_top-z_bottomlimit)) )
    enddo

    return
  end subroutine

end module mod_numfilter
