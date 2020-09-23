!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Boundary treatment
!!
!! @par Description
!!          Boundary treatment of model domain
!!          Additional forcing, Sponge layer, rayleigh dumping
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_bnd_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_BOUNDARY_driver_setup
  public :: ATMOS_BOUNDARY_driver_set
  public :: ATMOS_BOUNDARY_driver_finalize
  public :: ATMOS_BOUNDARY_driver_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: BND_QA    !> # of tracer at boundary
  integer,  public, allocatable :: BND_IQ(:) !> index of tracer

  real(RP), public, allocatable :: ATMOS_BOUNDARY_DENS(:,:,:)   !> reference DENS (with HALO)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_VELZ(:,:,:)   !> reference VELZ (with HALO)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_VELX(:,:,:)   !> reference VELX (with HALO)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_VELY(:,:,:)   !> reference VELY (with HALO)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_POTT(:,:,:)   !> reference POTT (with HALO)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_QTRC(:,:,:,:) !> reference QTRC (with HALO)

  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_DENS(:,:,:)   !> damping coefficient for DENS (0-1)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_VELZ(:,:,:)   !> damping coefficient for VELZ (0-1)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_VELX(:,:,:)   !> damping coefficient for VELX (0-1)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_VELY(:,:,:)   !> damping coefficient for VELY (0-1)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_POTT(:,:,:)   !> damping coefficient for POTT (0-1)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_QTRC(:,:,:,:) !> damping coefficient for QTRC (0-1)

  real(RP), public, allocatable :: ATMOS_BOUNDARY_MFLUX_OFFSET_X(:,:,:) !> mass flux offset (west, east)
  real(RP), public, allocatable :: ATMOS_BOUNDARY_MFLUX_OFFSET_Y(:,:,:) !> mass flux offset (south, north)

  real(RP), public              :: ATMOS_BOUNDARY_SMOOTHER_FACT  = 0.2_RP  !> fact for smoother to damping

  logical,  public              :: ATMOS_BOUNDARY_UPDATE_FLAG    = .false. !> switch for real case

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_BOUNDARY_var_fillhalo
  private :: ATMOS_BOUNDARY_alpha_fillhalo
  private :: ATMOS_BOUNDARY_ref_fillhalo
  private :: ATMOS_BOUNDARY_setalpha
  private :: ATMOS_BOUNDARY_setinitval
  private :: ATMOS_BOUNDARY_read
  private :: ATMOS_BOUNDARY_write
  private :: ATMOS_BOUNDARY_generate
  private :: ATMOS_BOUNDARY_initialize_file
  private :: ATMOS_BOUNDARY_initialize_online
  private :: ATMOS_BOUNDARY_update_file
  private :: ATMOS_BOUNDARY_update_online_parent
  private :: ATMOS_BOUNDARY_update_online_daughter
  private :: ATMOS_BOUNDARY_firstsend
  private :: ATMOS_BOUNDARY_send
  private :: ATMOS_BOUNDARY_recv

  abstract interface
     subroutine getbnd( &
          bnd_DENS, &
          bnd_VELZ, &
          bnd_VELX, &
          bnd_VELY, &
          bnd_POTT, &
          bnd_QTRC, &
          now_step, &
          update_step )
       use scale_precision
       implicit none

       real(RP), intent(out) :: bnd_DENS(:,:,:)
       real(RP), intent(out) :: bnd_VELZ(:,:,:)
       real(RP), intent(out) :: bnd_VELX(:,:,:)
       real(RP), intent(out) :: bnd_VELY(:,:,:)
       real(RP), intent(out) :: bnd_POTT(:,:,:)
       real(RP), intent(out) :: bnd_QTRC(:,:,:,:)
       integer,  intent(in)  :: now_step
       integer,  intent(in)  :: update_step
     end subroutine getbnd
  end interface

  procedure(getbnd), pointer :: get_boundary => NULL()
  private :: get_boundary
  private :: get_boundary_same_parent
  private :: get_boundary_nearest_neighbor
  private :: get_boundary_lerp_initpoint
  private :: get_boundary_lerp_midpoint

  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: ATMOS_BOUNDARY_TYPE         = 'NONE'
  character(len=H_LONG),  private :: ATMOS_BOUNDARY_IN_BASENAME  = ''
  logical,                private :: ATMOS_BOUNDARY_IN_CHECK_COORDINATES = .true.
  logical,                private :: ATMOS_BOUNDARY_IN_AGGREGATE
  character(len=H_LONG),  private :: ATMOS_BOUNDARY_OUT_BASENAME = ''
  character(len=H_MID),   private :: ATMOS_BOUNDARY_OUT_TITLE    = 'SCALE-RM BOUNDARY CONDITION'  !< title of the output file
  character(len=H_SHORT), private :: ATMOS_BOUNDARY_OUT_DTYPE    = 'DEFAULT'                      !< REAL4 or REAL8
  logical,                private :: ATMOS_BOUNDARY_OUT_AGGREGATE

  logical,               private :: ATMOS_BOUNDARY_USE_DENS     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELZ     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELX     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELY     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_PT       = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_QV       = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_QHYD     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_CHEM     = .false. ! read from file?

  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELZ   =   0.0_RP ! velocity w      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELX   =   0.0_RP ! velocity u      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELY   =   0.0_RP ! velocity v      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_PT     = 300.0_RP ! potential temp. at boundary, 300 [K]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_QTRC   =   0.0_RP ! tracer          at boundary, 0   [kg/kg]

  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_DENS = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_VELZ = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_VELX = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_VELY = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_PT   = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_QTRC = 1.0_RP ! alpha factor again default

  real(RP),              private :: ATMOS_BOUNDARY_FRACZ        =   1.0_RP ! fraction of boundary region for dumping (z) (0-1)
  real(RP),              private :: ATMOS_BOUNDARY_FRACX        =   1.0_RP ! fraction of boundary region for dumping (x) (0-1)
  real(RP),              private :: ATMOS_BOUNDARY_FRACY        =   1.0_RP ! fraction of boundary region for dumping (y) (0-1)
  real(RP),              private :: ATMOS_BOUNDARY_tauz                    ! maximum value for damping tau (z) [s]
  real(RP),              private :: ATMOS_BOUNDARY_taux                    ! maximum value for damping tau (x) [s]
  real(RP),              private :: ATMOS_BOUNDARY_tauy                    ! maximum value for damping tau (y) [s]

  real(DP), private              :: ATMOS_BOUNDARY_UPDATE_DT    =  0.0_DP ! inteval time of boudary data update [s]
  integer,  private              :: UPDATE_NSTEP

  logical,               private :: ATMOS_GRID_NUDGING_uv       = .false.  ! grid nudging
  logical,               private :: ATMOS_GRID_NUDGING_pt       = .false.  ! grid nudging
  logical,               private :: ATMOS_GRID_NUDGING_qv       = .false.  ! grid nudging
  real(RP),              private :: ATMOS_GRID_NUDGING_tau                 ! Damping tau for grid nudging [s]

  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_DENS(:,:,:,:)   ! reference DENS (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELZ(:,:,:,:)   ! reference VELZ (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELX(:,:,:,:)   ! reference VELX (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELY(:,:,:,:)   ! reference VELY (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_POTT(:,:,:,:)   ! reference POTT (with HALO)
  real(RP),              private, allocatable, target :: ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,:) ! reference QTRC (with HALO)

  ! work
  real(RP),              private, allocatable, target :: Q_WORK(:,:,:,:) ! QV + Qe


  character(len=H_SHORT), private :: ATMOS_BOUNDARY_interp_TYPE = 'lerp_initpoint' ! type of boundary interporation

  integer,               private :: ATMOS_BOUNDARY_START_DATE(6) = (/ -9999, 0, 0, 0, 0, 0 /) ! boundary initial date

  integer,               private :: ATMOS_BOUNDARY_fid = -1

  integer,               private :: now_step
  integer,               private :: boundary_timestep = 0
  logical,               private :: ATMOS_BOUNDARY_LINEAR_V = .false.  ! linear or non-linear profile of relax region
  logical,               private :: ATMOS_BOUNDARY_LINEAR_H = .true.   ! linear or non-linear profile of relax region
  real(RP),              private :: ATMOS_BOUNDARY_EXP_H    = 2.0_RP   ! factor of non-linear profile of relax region
  logical,               private :: ATMOS_BOUNDARY_ONLINE   = .false.  ! boundary online update by communicate inter-domain
  logical,               private :: ATMOS_BOUNDARY_ONLINE_MASTER = .false.  ! master domain in communicate inter-domain

  logical,               private :: ATMOS_BOUNDARY_DENS_ADJUST  = .false.
  real(RP),              private :: ATMOS_BOUNDARY_DENS_ADJUST_tau = -1.0_RP

  logical,               private :: do_parent_process       = .false.
  logical,               private :: do_daughter_process     = .false.
  logical,               private :: l_bnd = .false.

  real(DP),              private :: boundary_time_initdaysec

  integer,               private :: ref_size = 3
  integer,               private :: ref_old  = 1
  integer,               private :: ref_now  = 2
  integer,               private :: ref_new  = 3

  ! for mass flux offset
  real(DP),              private :: MASSTOT_now = 0.0_DP
  real(DP),              private :: MASSFLX_now = 0.0_DP
  real(RP), allocatable, private :: AREAZUY_W(:,:), AREAZUY_E(:,:)
  real(RP), allocatable, private :: OFFSET_TIME_FACT(:)
  real(RP), allocatable, private :: MFLUX_OFFSET_X(:,:,:,:)
  real(RP), allocatable, private :: MFLUX_OFFSET_Y(:,:,:,:)
  real(RP), allocatable, private, target :: zero_x(:,:), zero_y(:,:)


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_BOUNDARY_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_time, only: &
       DT => TIME_DTSEC
    use scale_file, only: &
       FILE_AGGREGATE
    use scale_comm_cartesC_nest, only: &
       ONLINE_BOUNDARY_DIAGQHYD, &
       ONLINE_BOUNDARY_USE_QHYD, &
       USE_NESTING,         &
       OFFLINE,             &
       ONLINE_IAM_PARENT,   &
       ONLINE_IAM_DAUGHTER, &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       I_QV
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_phy_ch_vars, only: &
       QS_CH, &
       QE_CH
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_AREAZUY_X
    implicit none

    namelist / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_TYPE,           &
       ATMOS_BOUNDARY_IN_BASENAME,    &
       ATMOS_BOUNDARY_IN_CHECK_COORDINATES, &
       ATMOS_BOUNDARY_IN_AGGREGATE,   &
       ATMOS_BOUNDARY_OUT_BASENAME,   &
       ATMOS_BOUNDARY_OUT_TITLE,      &
       ATMOS_BOUNDARY_OUT_DTYPE,      &
       ATMOS_BOUNDARY_OUT_AGGREGATE,  &
       ATMOS_BOUNDARY_USE_VELZ,       &
       ATMOS_BOUNDARY_USE_VELX,       &
       ATMOS_BOUNDARY_USE_VELY,       &
       ATMOS_BOUNDARY_USE_PT,         &
       ATMOS_BOUNDARY_USE_DENS,       &
       ATMOS_BOUNDARY_USE_QV,         &
       ATMOS_BOUNDARY_USE_QHYD,       &
       ATMOS_BOUNDARY_USE_CHEM,       &
       ATMOS_BOUNDARY_DENS_ADJUST,    &
       ATMOS_BOUNDARY_DENS_ADJUST_tau, &
       ATMOS_BOUNDARY_VALUE_VELZ,     &
       ATMOS_BOUNDARY_VALUE_VELX,     &
       ATMOS_BOUNDARY_VALUE_VELY,     &
       ATMOS_BOUNDARY_VALUE_PT,       &
       ATMOS_BOUNDARY_VALUE_QTRC,     &
       ATMOS_BOUNDARY_ALPHAFACT_DENS, &
       ATMOS_BOUNDARY_ALPHAFACT_VELZ, &
       ATMOS_BOUNDARY_ALPHAFACT_VELX, &
       ATMOS_BOUNDARY_ALPHAFACT_VELY, &
       ATMOS_BOUNDARY_ALPHAFACT_PT,   &
       ATMOS_BOUNDARY_ALPHAFACT_QTRC, &
       ATMOS_BOUNDARY_SMOOTHER_FACT,  &
       ATMOS_BOUNDARY_FRACZ,          &
       ATMOS_BOUNDARY_FRACX,          &
       ATMOS_BOUNDARY_FRACY,          &
       ATMOS_BOUNDARY_tauz,           &
       ATMOS_BOUNDARY_taux,           &
       ATMOS_BOUNDARY_tauy,           &
       ATMOS_BOUNDARY_UPDATE_DT,      &
       ATMOS_BOUNDARY_START_DATE,     &
       ATMOS_BOUNDARY_LINEAR_V,       &
       ATMOS_BOUNDARY_LINEAR_H,       &
       ATMOS_BOUNDARY_EXP_H,          &
       ATMOS_BOUNDARY_interp_TYPE,    &
       ATMOS_GRID_NUDGING_uv,         &
       ATMOS_GRID_NUDGING_pt,         &
       ATMOS_GRID_NUDGING_qv,         &
       ATMOS_GRID_NUDGING_tau

    integer :: k, i, j, iq
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_BOUNDARY_setup",*) 'Setup'


    ATMOS_BOUNDARY_tauz = DT * 10.0_RP
    ATMOS_BOUNDARY_taux = DT * 10.0_RP
    ATMOS_BOUNDARY_tauy = DT * 10.0_RP

    ATMOS_BOUNDARY_IN_AGGREGATE  = FILE_AGGREGATE
    ATMOS_BOUNDARY_OUT_AGGREGATE = FILE_AGGREGATE

    ATMOS_GRID_NUDGING_tau = 10.0_RP * 24.0_RP * 3600.0_RP   ! 10days [s]

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_BOUNDARY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_BOUNDARY_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_BOUNDARY_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_BOUNDARY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_BOUNDARY)

    call IO_filename_replace( ATMOS_BOUNDARY_IN_BASENAME, 'ATMOS_BOUNDARY_IN_BASENAME' )

    ! setting switches
    if( .NOT. USE_NESTING ) then
       ATMOS_BOUNDARY_ONLINE = .false.
    else
       if( OFFLINE ) then
          ATMOS_BOUNDARY_ONLINE = .false.
       else
          ATMOS_BOUNDARY_ONLINE = .true.
       endif
    endif
    do_parent_process   = .false.
    do_daughter_process = .false.
    ATMOS_BOUNDARY_ONLINE_MASTER = .false.
    if ( ATMOS_BOUNDARY_ONLINE ) then
       if ( ONLINE_IAM_PARENT ) then
          do_parent_process = .true.
          if ( .NOT. ONLINE_IAM_DAUGHTER ) then
             ATMOS_BOUNDARY_ONLINE_MASTER = .true.
          endif
       endif
       if ( ONLINE_IAM_DAUGHTER ) then
          do_daughter_process = .true.
          ATMOS_BOUNDARY_USE_QHYD = ONLINE_BOUNDARY_USE_QHYD
       endif
    endif

    allocate( BND_IQ(QA) )
    BND_IQ(:) = -1
    BND_QA = 0
    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       BND_QA = BND_QA + 1
       BND_IQ(I_QV) = BND_QA
       if( ATMOS_BOUNDARY_USE_QHYD ) then
          do iq = QS_MP+1, QE_MP
             BND_QA = BND_QA + 1
             BND_IQ(iq) = BND_QA
          end do
       end if
    end if
    if( ATMOS_BOUNDARY_USE_CHEM ) then
       do iq = QS_CH, QE_CH
          BND_QA = BND_QA + 1
          BND_IQ(iq) = BND_QA
       end do
    endif

    if ( ATMOS_BOUNDARY_DENS_ADJUST_tau <= 0.0_RP ) then
       ATMOS_BOUNDARY_DENS_ADJUST_tau = max( real(ATMOS_BOUNDARY_UPDATE_DT,kind=RP) / 6.0_RP, &
                                             ATMOS_BOUNDARY_taux, ATMOS_BOUNDARY_tauy )
    end if



    allocate( ATMOS_BOUNDARY_DENS(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_VELZ(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_VELX(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_VELY(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_POTT(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_QTRC(KA,IA,JA,BND_QA) )
    ATMOS_BOUNDARY_DENS(:,:,:)   = UNDEF
    ATMOS_BOUNDARY_VELZ(:,:,:)   = UNDEF
    ATMOS_BOUNDARY_VELX(:,:,:)   = UNDEF
    ATMOS_BOUNDARY_VELY(:,:,:)   = UNDEF
    ATMOS_BOUNDARY_POTT(:,:,:)   = UNDEF
    ATMOS_BOUNDARY_QTRC(:,:,:,:) = UNDEF

    allocate( ATMOS_BOUNDARY_alpha_DENS(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_alpha_VELZ(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_alpha_VELX(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_alpha_VELY(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_alpha_POTT(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_alpha_QTRC(KA,IA,JA,BND_QA) )
    ATMOS_BOUNDARY_alpha_DENS(:,:,:)   = 0.0_RP
    ATMOS_BOUNDARY_alpha_VELZ(:,:,:)   = 0.0_RP
    ATMOS_BOUNDARY_alpha_VELX(:,:,:)   = 0.0_RP
    ATMOS_BOUNDARY_alpha_VELY(:,:,:)   = 0.0_RP
    ATMOS_BOUNDARY_alpha_POTT(:,:,:)   = 0.0_RP
    ATMOS_BOUNDARY_alpha_QTRC(:,:,:,:) = 0.0_RP

    allocate( ATMOS_BOUNDARY_MFLUX_OFFSET_X(KA,JA,2) )
    allocate( ATMOS_BOUNDARY_MFLUX_OFFSET_Y(KA,IA,2) )
    ATMOS_BOUNDARY_MFLUX_OFFSET_X(:,:,:) = 0.0_RP
    ATMOS_BOUNDARY_MFLUX_OFFSET_Y(:,:,:) = 0.0_RP

    if ( ATMOS_BOUNDARY_TYPE == 'REAL' .OR. do_daughter_process ) then
       l_bnd = .true.
    else
       l_bnd = .false.
    end if

    if ( l_bnd ) then

       select case(ATMOS_BOUNDARY_interp_TYPE)
       case('same_parent')
          get_boundary => get_boundary_same_parent
       case('nearest_neighbor')
          get_boundary => get_boundary_nearest_neighbor
       case('lerp_initpoint')
          get_boundary => get_boundary_lerp_initpoint
       case('lerp_midpoint')
          get_boundary => get_boundary_lerp_midpoint
       case default
          LOG_ERROR("ATMOS_BOUNDARY_setup",*) 'Wrong parameter in ATMOS_BOUNDARY_interp_TYPE. Check!'
          call PRC_abort
       end select

       allocate( ATMOS_BOUNDARY_ref_DENS(KA,IA,JA,ref_size) )
       allocate( ATMOS_BOUNDARY_ref_VELZ(KA,IA,JA,ref_size) )
       allocate( ATMOS_BOUNDARY_ref_VELX(KA,IA,JA,ref_size) )
       allocate( ATMOS_BOUNDARY_ref_VELY(KA,IA,JA,ref_size) )
       allocate( ATMOS_BOUNDARY_ref_POTT(KA,IA,JA,ref_size) )
       allocate( ATMOS_BOUNDARY_ref_QTRC(KA,IA,JA,BND_QA,ref_size) )
       ATMOS_BOUNDARY_ref_DENS(:,:,:,:)   = UNDEF
       ATMOS_BOUNDARY_ref_VELZ(:,:,:,:)   = UNDEF
       ATMOS_BOUNDARY_ref_VELX(:,:,:,:)   = UNDEF
       ATMOS_BOUNDARY_ref_VELY(:,:,:,:)   = UNDEF
       ATMOS_BOUNDARY_ref_POTT(:,:,:,:)   = UNDEF
       ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,:) = UNDEF

       ! initialize boundary value (reading file or waiting parent domain)
       if ( do_daughter_process ) then
          call ATMOS_BOUNDARY_initialize_online
       else
          if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
             call ATMOS_BOUNDARY_initialize_file
          else
             LOG_ERROR("ATMOS_BOUNDARY_setup",*) 'You need specify ATMOS_BOUNDARY_IN_BASENAME'
             call PRC_abort
          endif
       endif

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .true.

       ! for mass flux offset
       allocate( AREAZUY_W(KA,JA), AREAZUY_E(KA,JA) )
       allocate( MFLUX_OFFSET_X(KA,JA,2,2), MFLUX_OFFSET_Y(KA,IA,2,2) )

       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          AREAZUY_W(k,j) = ATMOS_GRID_CARTESC_REAL_AREAZUY_X(k,IS-1,j)
          AREAZUY_E(k,j) = ATMOS_GRID_CARTESC_REAL_AREAZUY_X(k,IE  ,j)
       end do
       end do
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          MFLUX_OFFSET_X(k,j,:,:) = 0.0_RP
       end do
       end do
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          MFLUX_OFFSET_Y(k,i,:,:) = 0.0_RP
       end do
       end do

    elseif ( ATMOS_BOUNDARY_TYPE == 'NONE' ) then

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'CONST' ) then

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'INIT' ) then

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'OFFLINE' ) then

       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          call ATMOS_BOUNDARY_read
       else
          LOG_ERROR("ATMOS_BOUNDARY_setup",*) 'You need specify ATMOS_BOUNDARY_IN_BASENAME'
          call PRC_abort
       endif

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    else
       LOG_ERROR("ATMOS_BOUNDARY_setup",*) 'unsupported ATMOS_BOUNDARY_TYPE. Check!', trim(ATMOS_BOUNDARY_TYPE)
       call PRC_abort
    endif

    if ( USE_NESTING ) ATMOS_BOUNDARY_UPDATE_FLAG = .true.


    !----- report data -----
    LOG_NEWLINE
    LOG_INFO("ATMOS_BOUNDARY_setup",*) 'Atmospheric boundary parameters '
    LOG_INFO_CONT(*) 'Atmospheric boundary type                      : ', ATMOS_BOUNDARY_TYPE
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Is VELZ used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_VELZ
    LOG_INFO_CONT(*) 'Is VELX used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_VELX
    LOG_INFO_CONT(*) 'Is VELY used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_VELY
    LOG_INFO_CONT(*) 'Is PT used in atmospheric boundary?            : ', ATMOS_BOUNDARY_USE_PT
    LOG_INFO_CONT(*) 'Is DENS used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_DENS
    LOG_INFO_CONT(*) 'Is QV   used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_QV
    LOG_INFO_CONT(*) 'Is QHYD used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_QHYD
    LOG_INFO_CONT(*) 'Is CHEM used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_CHEM
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Atmospheric boundary VELZ values               : ', ATMOS_BOUNDARY_VALUE_VELZ
    LOG_INFO_CONT(*) 'Atmospheric boundary VELX values               : ', ATMOS_BOUNDARY_VALUE_VELX
    LOG_INFO_CONT(*) 'Atmospheric boundary VELY values               : ', ATMOS_BOUNDARY_VALUE_VELY
    LOG_INFO_CONT(*) 'Atmospheric boundary PT values                 : ', ATMOS_BOUNDARY_VALUE_PT
    LOG_INFO_CONT(*) 'Atmospheric boundary QTRC values               : ', ATMOS_BOUNDARY_VALUE_QTRC
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Atmospheric boundary smoother factor           : ', ATMOS_BOUNDARY_SMOOTHER_FACT
    LOG_INFO_CONT(*) 'Atmospheric boundary z-fraction                : ', ATMOS_BOUNDARY_FRACZ
    LOG_INFO_CONT(*) 'Atmospheric boundary x-fraction                : ', ATMOS_BOUNDARY_FRACX
    LOG_INFO_CONT(*) 'Atmospheric boundary y-fraction                : ', ATMOS_BOUNDARY_FRACY
    LOG_INFO_CONT(*) 'Atmospheric boundary z-relaxation time         : ', ATMOS_BOUNDARY_TAUZ
    LOG_INFO_CONT(*) 'Atmospheric boundary x-relaxation time         : ', ATMOS_BOUNDARY_TAUX
    LOG_INFO_CONT(*) 'Atmospheric boundary y-relaxation time         : ', ATMOS_BOUNDARY_TAUY
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Atmospheric boundary update dt                 : ', ATMOS_BOUNDARY_UPDATE_DT
    LOG_INFO_CONT(*) 'Atmospheric boundary start date                : ', ATMOS_BOUNDARY_START_DATE(:)
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Linear profile in vertically relax region      : ', ATMOS_BOUNDARY_LINEAR_V
    LOG_INFO_CONT(*) 'Linear profile in horizontally relax region    : ', ATMOS_BOUNDARY_LINEAR_H
    LOG_INFO_CONT(*) 'Non-linear factor in horizontally relax region : ', ATMOS_BOUNDARY_EXP_H
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Online nesting for lateral boundary            : ', ATMOS_BOUNDARY_ONLINE

    LOG_INFO_CONT(*) 'Does lateral boundary exist in this domain?    : ', l_bnd
    if ( l_bnd ) then
       LOG_INFO_CONT(*) 'Lateral boundary interporation type            : ', trim(ATMOS_BOUNDARY_interp_TYPE)
    endif
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Is grid nudging used for VELX & VELY?          : ', ATMOS_GRID_NUDGING_uv
    LOG_INFO_CONT(*) 'Is grid nudging used for POTT?                 : ', ATMOS_GRID_NUDGING_pt
    LOG_INFO_CONT(*) 'Is grid nudging used for QV?                   : ', ATMOS_GRID_NUDGING_qv
    LOG_INFO_CONT(*) 'Relaxation time for grid nudging               : ', ATMOS_GRID_NUDGING_tau

    LOG_INFO_CONT(*) 'Density adjustment                             : ', ATMOS_BOUNDARY_DENS_ADJUST
    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       LOG_INFO_CONT(*) 'Density relaxation time                        : ', ATMOS_BOUNDARY_DENS_ADJUST_tau
    end if

    if ( ONLINE_BOUNDARY_DIAGQHYD ) then
       allocate( Q_WORK(KA,IA,JA,NESTQA) )
    end if

    allocate( zero_x(KA,JA), zero_y(KA,IA) )
    !$omp parallel do
    do j = JS, JE
    do k = KS, KE
       zero_x(k,j) = 0.0_RP
    end do
    end do
    !$omp parallel do
    do i = IS, IE
    do k = KS, KE
       zero_y(k,i) = 0.0_RP
    end do
    end do



    return
  end subroutine ATMOS_BOUNDARY_driver_setup

  !-----------------------------------------------------------------------------
  !> set
  subroutine ATMOS_BOUNDARY_driver_set
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       QV,   &
       Qe
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none

    if ( do_parent_process ) then !online [parent]
       call ATMOS_BOUNDARY_firstsend( &
            DENS, MOMZ, MOMX, MOMY, RHOT, QTRC(:,:,:,QS_MP:QE_MP), QV, Qe )
    end if

    if ( l_bnd ) then

       ! initialize boundary value (reading file or waiting parent domain)
       if ( do_daughter_process ) then
          call ATMOS_BOUNDARY_set_online
       else
          if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
             call ATMOS_BOUNDARY_set_file
          endif
       endif

    elseif ( ATMOS_BOUNDARY_TYPE == 'CONST' ) then

       call ATMOS_BOUNDARY_generate

    elseif ( ATMOS_BOUNDARY_TYPE == 'INIT' ) then

       call ATMOS_BOUNDARY_setinitval( DENS, & ! [IN]
                                       MOMZ, & ! [IN]
                                       MOMX, & ! [IN]
                                       MOMY, & ! [IN]
                                       RHOT, & ! [IN]
                                       QTRC  ) ! [IN]
    endif

    if( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
       call ATMOS_BOUNDARY_write
    endif

    if ( ATMOS_BOUNDARY_UPDATE_FLAG ) then

       call history_bnd( &
            ATMOS_BOUNDARY_DENS, &
            ATMOS_BOUNDARY_VELZ, &
            ATMOS_BOUNDARY_VELX, &
            ATMOS_BOUNDARY_VELY, &
            ATMOS_BOUNDARY_POTT, &
            ATMOS_BOUNDARY_QTRC )
    end if

    return
  end subroutine ATMOS_BOUNDARY_driver_set

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_var_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       ATMOS_BOUNDARY_DENS(   1:KS-1,i,j) = ATMOS_BOUNDARY_DENS(KS,i,j)
       ATMOS_BOUNDARY_VELZ(   1:KS-1,i,j) = ATMOS_BOUNDARY_VELZ(KS,i,j)
       ATMOS_BOUNDARY_VELX(   1:KS-1,i,j) = ATMOS_BOUNDARY_VELX(KS,i,j)
       ATMOS_BOUNDARY_VELY(   1:KS-1,i,j) = ATMOS_BOUNDARY_VELY(KS,i,j)
       ATMOS_BOUNDARY_POTT(   1:KS-1,i,j) = ATMOS_BOUNDARY_POTT(KS,i,j)

       ATMOS_BOUNDARY_DENS(KE+1:KA,  i,j) = ATMOS_BOUNDARY_DENS(KE,i,j)
       ATMOS_BOUNDARY_VELZ(KE+1:KA,  i,j) = ATMOS_BOUNDARY_VELZ(KE,i,j)
       ATMOS_BOUNDARY_VELX(KE+1:KA,  i,j) = ATMOS_BOUNDARY_VELX(KE,i,j)
       ATMOS_BOUNDARY_VELY(KE+1:KA,  i,j) = ATMOS_BOUNDARY_VELY(KE,i,j)
       ATMOS_BOUNDARY_POTT(KE+1:KA,  i,j) = ATMOS_BOUNDARY_POTT(KE,i,j)

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(   1:KS-1,i,j,iq) = ATMOS_BOUNDARY_QTRC(KS,i,j,iq)
          ATMOS_BOUNDARY_QTRC(KE+1:KA,  i,j,iq) = ATMOS_BOUNDARY_QTRC(KE,i,j,iq)
       end do
    end do
    end do

    call COMM_vars8( ATMOS_BOUNDARY_DENS(:,:,:),   1 )
    call COMM_vars8( ATMOS_BOUNDARY_VELZ(:,:,:),   2 )
    call COMM_vars8( ATMOS_BOUNDARY_VELX(:,:,:),   3 )
    call COMM_vars8( ATMOS_BOUNDARY_VELY(:,:,:),   4 )
    call COMM_vars8( ATMOS_BOUNDARY_POTT(:,:,:),   5 )
    do iq = 1, BND_QA
       call COMM_vars8( ATMOS_BOUNDARY_QTRC(:,:,:,iq), 5+iq )
    end do

    call COMM_wait ( ATMOS_BOUNDARY_DENS(:,:,:),   1, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_VELZ(:,:,:),   2, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_VELX(:,:,:),   3, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_VELY(:,:,:),   4, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_POTT(:,:,:),   5, .false. )
    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_QTRC(:,:,:,iq), 5+iq, .false. )
    end do

    return
  end subroutine ATMOS_BOUNDARY_var_fillhalo

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_alpha_fillhalo
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       ATMOS_BOUNDARY_alpha_DENS(   1:KS-1,i,j) = ATMOS_BOUNDARY_alpha_DENS(KS,i,j)
       ATMOS_BOUNDARY_alpha_VELZ(   1:KS-1,i,j) = ATMOS_BOUNDARY_alpha_VELZ(KS,i,j)
       ATMOS_BOUNDARY_alpha_VELX(   1:KS-1,i,j) = ATMOS_BOUNDARY_alpha_VELX(KS,i,j)
       ATMOS_BOUNDARY_alpha_VELY(   1:KS-1,i,j) = ATMOS_BOUNDARY_alpha_VELY(KS,i,j)
       ATMOS_BOUNDARY_alpha_POTT(   1:KS-1,i,j) = ATMOS_BOUNDARY_alpha_POTT(KS,i,j)

       ATMOS_BOUNDARY_alpha_DENS(KE+1:KA,  i,j) = ATMOS_BOUNDARY_alpha_DENS(KE,i,j)
       ATMOS_BOUNDARY_alpha_VELZ(KE+1:KA,  i,j) = ATMOS_BOUNDARY_alpha_VELZ(KE,i,j)
       ATMOS_BOUNDARY_alpha_VELX(KE+1:KA,  i,j) = ATMOS_BOUNDARY_alpha_VELX(KE,i,j)
       ATMOS_BOUNDARY_alpha_VELY(KE+1:KA,  i,j) = ATMOS_BOUNDARY_alpha_VELY(KE,i,j)
       ATMOS_BOUNDARY_alpha_POTT(KE+1:KA,  i,j) = ATMOS_BOUNDARY_alpha_POTT(KE,i,j)

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_alpha_QTRC(   1:KS-1,i,j,iq) = ATMOS_BOUNDARY_alpha_QTRC(KS,i,j,iq)
          ATMOS_BOUNDARY_alpha_QTRC(KE+1:KA,  i,j,iq) = ATMOS_BOUNDARY_alpha_QTRC(KE,i,j,iq)
       end do
    enddo
    enddo

    call COMM_vars8( ATMOS_BOUNDARY_alpha_DENS(:,:,:),   1 )
    call COMM_vars8( ATMOS_BOUNDARY_alpha_VELZ(:,:,:),   2 )
    call COMM_vars8( ATMOS_BOUNDARY_alpha_VELX(:,:,:),   3 )
    call COMM_vars8( ATMOS_BOUNDARY_alpha_VELY(:,:,:),   4 )
    call COMM_vars8( ATMOS_BOUNDARY_alpha_POTT(:,:,:),   5 )
    do iq = 1, BND_QA
       call COMM_vars8( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq), 5+iq )
    end do

    call COMM_wait ( ATMOS_BOUNDARY_alpha_DENS(:,:,:),   1, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_VELZ(:,:,:),   2, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_VELX(:,:,:),   3, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_VELY(:,:,:),   4, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_POTT(:,:,:),   5, .false. )
    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq), 5+iq, .false. )
    end do

    return
  end subroutine ATMOS_BOUNDARY_alpha_fillhalo

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_ref_fillhalo( &
       ref_idx )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    ! arguments
    integer, intent(in) :: ref_idx

    ! works
    integer :: i, j, iq
    !---------------------------------------------------------------------------

    !$omp parallel do default(none) private(i,j,iq) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JSB,JEB,ISB,IEB,ATMOS_BOUNDARY_ref_DENS,ref_idx,ATMOS_BOUNDARY_ref_VELZ) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELX,ATMOS_BOUNDARY_ref_VELY,ATMOS_BOUNDARY_ref_POTT) &
    !$omp shared(KA,KS,BND_QA,ATMOS_BOUNDARY_ref_QTRC,KE)
    do j = JSB, JEB
    do i = ISB, IEB
       ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,ref_idx) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,ref_idx) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,ref_idx) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,ref_idx) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,ref_idx) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,ref_idx)

       ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,ref_idx) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,ref_idx) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,ref_idx) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,ref_idx) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,ref_idx)
       ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,ref_idx) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,ref_idx)

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,ref_idx) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,ref_idx)
          ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,ref_idx) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,ref_idx)
       end do
    end do
    end do

    call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,ref_idx), 1 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,ref_idx), 2 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,ref_idx), 3 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,ref_idx), 4 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,ref_idx), 5 )

    do iq = 1, BND_QA
       call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,ref_idx), 5+iq )
    end do

    call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,ref_idx), 1, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,ref_idx), 2, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,ref_idx), 3, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,ref_idx), 4, .false. )
    call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,ref_idx), 5, .false. )

    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,ref_idx), 5+iq, .false. )
    end do

    return
  end subroutine ATMOS_BOUNDARY_ref_fillhalo

  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine ATMOS_BOUNDARY_setalpha
    use scale_const, only: &
       EPS => CONST_EPS, &
       PI  => CONST_PI
    use scale_atmos_grid_cartesC, only: &
       CBFZ => ATMOS_GRID_CARTESC_CBFZ, &
       CBFX => ATMOS_GRID_CARTESC_CBFX, &
       CBFY => ATMOS_GRID_CARTESC_CBFY, &
       FBFZ => ATMOS_GRID_CARTESC_FBFZ, &
       FBFX => ATMOS_GRID_CARTESC_FBFX, &
       FBFY => ATMOS_GRID_CARTESC_FBFY
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ
    use scale_atmos_hydrometeor, only: &
       I_QV
    real(RP) :: coef_z, alpha_z1, alpha_z2
    real(RP) :: coef_x, alpha_x1, alpha_x2
    real(RP) :: coef_y, alpha_y1, alpha_y2
    real(RP) :: alpha_zm, alpha_xm, alpha_ym
    real(RP) :: ee1, ee2

    real(RP) :: alpha_nug ! grid nudging in inner domain

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    ! check invalid fraction
    ATMOS_BOUNDARY_FRACZ = max( min( ATMOS_BOUNDARY_FRACZ, 1.0_RP ), EPS )
    ATMOS_BOUNDARY_FRACX = max( min( ATMOS_BOUNDARY_FRACX, 1.0_RP ), EPS )
    ATMOS_BOUNDARY_FRACY = max( min( ATMOS_BOUNDARY_FRACY, 1.0_RP ), EPS )

    if ( ATMOS_BOUNDARY_tauz <= 0.0_RP ) then ! invalid tau
       coef_z = 0.0_RP
    else
       coef_z = 1.0_RP / ATMOS_BOUNDARY_tauz
    endif

    if ( ATMOS_BOUNDARY_taux <= 0.0_RP ) then ! invalid tau
       coef_x = 0.0_RP
    else
       coef_x = 1.0_RP / ATMOS_BOUNDARY_taux
    endif

    if ( ATMOS_BOUNDARY_tauy <= 0.0_RP ) then ! invalid tau
       coef_y = 0.0_RP
    else
       coef_y = 1.0_RP / ATMOS_BOUNDARY_tauy
    endif

    if ( ATMOS_GRID_NUDGING_tau <= 0.0_RP ) then ! invalid tau
       alpha_nug = 0.0_RP
    else
       alpha_nug = 1.0_RP / ATMOS_GRID_NUDGING_tau
    endif

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KS,KE,CBFZ,ATMOS_BOUNDARY_FRACZ,FBFZ,ATMOS_BOUNDARY_LINEAR_V,coef_z,CBFX)            &
    !$omp shared(ATMOS_BOUNDARY_FRACX,PI,FBFX,ATMOS_BOUNDARY_LINEAR_H,coef_x)     &
    !$omp shared(ATMOS_BOUNDARY_EXP_H,CBFY,ATMOS_BOUNDARY_FRACY,FBFY,coef_y,l_bnd)     &
    !$omp shared(do_daughter_process) &
    !$omp shared(ONLINE_USE_VELZ,ATMOS_BOUNDARY_USE_VELZ,ATMOS_BOUNDARY_alpha_VELZ,ATMOS_BOUNDARY_ALPHAFACT_VELZ) &
    !$omp shared(ATMOS_BOUNDARY_USE_DENS,ATMOS_BOUNDARY_alpha_DENS,ATMOS_BOUNDARY_ALPHAFACT_DENS)        &
    !$omp shared(ATMOS_BOUNDARY_USE_VELX,ATMOS_BOUNDARY_alpha_VELX,ATMOS_BOUNDARY_ALPHAFACT_VELX)        &
    !$omp shared(ATMOS_BOUNDARY_USE_VELY,ATMOS_BOUNDARY_alpha_VELY,ATMOS_BOUNDARY_ALPHAFACT_VELY)        &
    !$omp shared(ATMOS_BOUNDARY_USE_PT,ATMOS_BOUNDARY_alpha_POTT,ATMOS_BOUNDARY_ALPHAFACT_PT)        &
    !$omp shared(ATMOS_BOUNDARY_USE_QV,ATMOS_BOUNDARY_alpha_QTRC,ATMOS_BOUNDARY_ALPHAFACT_QTRC)          &
    !$omp shared(ATMOS_BOUNDARY_DENS_ADJUST,ATMOS_BOUNDARY_DENS_ADJUST_tau,ATMOS_BOUNDARY_UPDATE_DT) &
    !$omp shared(BND_QA,BND_IQ,I_QV) &
    !$omp shared(ATMOS_GRID_NUDGING_uv,ATMOS_GRID_NUDGING_pt,ATMOS_GRID_NUDGING_qv) &
    !$omp shared(alpha_nug) &
    !$omp private(i,j,k,iq) &
    !$omp private(ee1,ee2,alpha_z1,alpha_z2,alpha_x1,alpha_x2,alpha_y1,alpha_y2,alpha_zm,alpha_xm,alpha_ym)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ee1 = CBFZ(k)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACZ ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACZ ) / ATMOS_BOUNDARY_FRACZ
       endif

       ee2 = FBFZ(k)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACZ ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACZ ) / ATMOS_BOUNDARY_FRACZ
       endif

       if ( .not. ATMOS_BOUNDARY_LINEAR_V ) then
          if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
             ee1 = 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
             ee1 = 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          endif
          if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
             ee2 = 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
          elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
             ee2 = 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
          endif
       endif

       alpha_z1 = coef_z * ee1
       alpha_z2 = coef_z * ee2
       alpha_zm = ee1 / ATMOS_BOUNDARY_DENS_ADJUST_tau


       ee1 = CBFX(i)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACX ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACX ) / ATMOS_BOUNDARY_FRACX
       endif

       ee2 = FBFX(i)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACX ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACX ) / ATMOS_BOUNDARY_FRACX
       endif

       if ( .not. ATMOS_BOUNDARY_LINEAR_H ) then
          ee1 = ee1 * exp( -(1.0_RP-ee1) * ATMOS_BOUNDARY_EXP_H )
          ee2 = ee2 * exp( -(1.0_RP-ee2) * ATMOS_BOUNDARY_EXP_H )
       end if

       alpha_x1 = coef_x * ee1
       alpha_x2 = coef_x * ee2
       alpha_xm = ee1 / ATMOS_BOUNDARY_DENS_ADJUST_tau


       ee1 = CBFY(j)
       if ( ee1 <= 1.0_RP - ATMOS_BOUNDARY_FRACY ) then
          ee1 = 0.0_RP
       else
          ee1 = ( ee1 - 1.0_RP + ATMOS_BOUNDARY_FRACY ) / ATMOS_BOUNDARY_FRACY
       endif

       ee2 = FBFY(j)
       if ( ee2 <= 1.0_RP - ATMOS_BOUNDARY_FRACY ) then
          ee2 = 0.0_RP
       else
          ee2 = ( ee2 - 1.0_RP + ATMOS_BOUNDARY_FRACY ) / ATMOS_BOUNDARY_FRACY
       endif

       if ( .not. ATMOS_BOUNDARY_LINEAR_H ) then
          ee1 = ee1 * exp( -(1.0_RP-ee1) * ATMOS_BOUNDARY_EXP_H )
          ee2 = ee2 * exp( -(1.0_RP-ee2) * ATMOS_BOUNDARY_EXP_H )
       end if

       alpha_y1 = coef_y * ee1
       alpha_y2 = coef_y * ee2
       alpha_ym = ee1 / ATMOS_BOUNDARY_DENS_ADJUST_tau


       if ( l_bnd ) then
          if (      (       do_daughter_process .and. ONLINE_USE_VELZ ) &
               .or. ( .not. do_daughter_process .and. ATMOS_BOUNDARY_USE_VELZ ) ) then
             ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_z2, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELZ
          else
             ATMOS_BOUNDARY_alpha_VELZ(:,:,:) = 0.0_RP
          end if
          if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_zm, alpha_xm, alpha_ym )
          else
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_DENS
          end if
          ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_z1, alpha_x2, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELX
          ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_z1, alpha_x1, alpha_y2 ) * ATMOS_BOUNDARY_ALPHAFACT_VELY
          ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_PT
          do iq = 1, BND_QA
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
          end do
       else
          if ( ATMOS_BOUNDARY_USE_DENS ) then
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_DENS
          else
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = 0.0_RP
          end if
          if ( ATMOS_BOUNDARY_USE_VELZ ) then
             ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_z2, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELZ
          else
             ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = 0.0_RP
          end if
          if ( ATMOS_BOUNDARY_USE_VELX ) then
             ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_z1, alpha_x2, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELX
          else
             ATMOS_BOUNDARY_alpha_VELX(k,i,j) = 0.0_RP
          end if
          if ( ATMOS_BOUNDARY_USE_VELY ) then
             ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_z1, alpha_x1, alpha_y2 ) * ATMOS_BOUNDARY_ALPHAFACT_VELY
          else
             ATMOS_BOUNDARY_alpha_VELY(k,i,j) = 0.0_RP
          end if
          if ( ATMOS_BOUNDARY_USE_PT ) then
             ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_PT
          else
             ATMOS_BOUNDARY_alpha_POTT(k,i,j) = 0.0_RP
          end if
          do iq = 1, BND_QA
             if ( I_QV > 0 .and. iq == BND_IQ(I_QV) ) then
                if ( ATMOS_BOUNDARY_USE_QV ) then
                   ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
                else
                   ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = 0.0_RP
                endif
             else
                ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
             end if
          end do

          ! internal grid nudging
          if ( ATMOS_GRID_NUDGING_uv ) then
             ATMOS_BOUNDARY_alpha_VELX(k,i,j)   = max( ATMOS_BOUNDARY_alpha_VELX(k,i,j), alpha_nug )
             ATMOS_BOUNDARY_alpha_VELY(k,i,j)   = max( ATMOS_BOUNDARY_alpha_VELY(k,i,j), alpha_nug )
          endif
          if ( ATMOS_GRID_NUDGING_pt ) then
             ATMOS_BOUNDARY_alpha_POTT(k,i,j)   = max( ATMOS_BOUNDARY_alpha_POTT(k,i,j), alpha_nug )
          endif
          if ( ATMOS_GRID_NUDGING_qv ) then
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,1) = max( ATMOS_BOUNDARY_alpha_QTRC(k,i,j,1), alpha_nug )
          endif

       end if
    enddo
    enddo
    enddo


    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_setalpha

  !-----------------------------------------------------------------------------
  !> Read boundary data
  subroutine ATMOS_BOUNDARY_setinitval( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,BND_QA)

    integer :: i, j, k, iq, iqb
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_DENS(k,i,j) = DENS(k,i,j)
       ATMOS_BOUNDARY_VELZ(k,i,j) = MOMZ(k,i,j) / ( DENS(k,i,j)+DENS(k+1,i,  j  ) ) * 2.0_RP
       ATMOS_BOUNDARY_POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       do iq = 1, QA
          iqb = BND_IQ(iq)
          if ( iqb > 0 ) ATMOS_BOUNDARY_QTRC(k,i,j,iqb) = QTRC(k,i,j,iq)
       end do
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA-1
    do k = KS, KE
       ATMOS_BOUNDARY_VELX(k,i,j) = MOMX(k,i,j) / ( DENS(k,i,j)+DENS(k,  i+1,j  ) ) * 2.0_RP
    enddo
    enddo
    enddo
    do j = 1, JA
    do k = KS, KE
       ATMOS_BOUNDARY_VELX(k,IA,j) = MOMX(k,IA,j) / DENS(k,IA,j)
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_VELY(k,i,j) = MOMY(k,i,j) / ( DENS(k,i,j)+DENS(k,  i,  j+1) ) * 2.0_RP
    enddo
    enddo
    enddo
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_VELY(k,i,JA) = MOMY(k,i,JA) / DENS(k,i,JA)
    enddo
    enddo

    call ATMOS_BOUNDARY_var_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_setinitval

  !-----------------------------------------------------------------------------
  !> Read boundary data
  subroutine ATMOS_BOUNDARY_read
    use scale_prc, only: &
       PRC_abort
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_read, &
       FILE_CARTESC_close
    implicit none

    integer :: fid
    integer :: iq, iqb
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( ATMOS_BOUNDARY_IN_BASENAME, fid, aggregate=ATMOS_BOUNDARY_IN_AGGREGATE )

    if ( ATMOS_BOUNDARY_IN_CHECK_COORDINATES ) then
       call FILE_CARTESC_check_coordinates( fid, atmos=.true. )
    end if

    if (      ATMOS_BOUNDARY_USE_DENS &
         .OR. ATMOS_BOUNDARY_USE_VELZ &
         .OR. ATMOS_BOUNDARY_USE_VELX &
         .OR. ATMOS_BOUNDARY_USE_VELY &
         .OR. ATMOS_BOUNDARY_USE_PT &
         ) then
       call FILE_CARTESC_read( fid, 'DENS', 'ZXY', ATMOS_BOUNDARY_DENS(:,:,:) )
    end if
    if ( ATMOS_BOUNDARY_USE_DENS ) then
       call FILE_CARTESC_read( fid, 'ALPHA_DENS', 'ZXY', ATMOS_BOUNDARY_alpha_DENS(:,:,:) )
    endif

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FILE_CARTESC_read( fid, 'VELZ', 'ZHXY', ATMOS_BOUNDARY_VELZ(:,:,:) )
       call FILE_CARTESC_read( fid, 'ALPHA_VELZ', 'ZHXY', ATMOS_BOUNDARY_alpha_VELZ(:,:,:) )
    endif

    if ( ATMOS_BOUNDARY_USE_VELX ) then
       call FILE_CARTESC_read( fid, 'VELX', 'ZXHY', ATMOS_BOUNDARY_VELX(:,:,:) )
       call FILE_CARTESC_read( fid, 'ALPHA_VELX', 'ZXHY', ATMOS_BOUNDARY_alpha_VELX(:,:,:) )
    endif

    if ( ATMOS_BOUNDARY_USE_VELY ) then
       call FILE_CARTESC_read( fid, 'VELY', 'ZXYH', ATMOS_BOUNDARY_VELY(:,:,:) )
       call FILE_CARTESC_read( fid, 'ALPHA_VELY', 'ZXYH', ATMOS_BOUNDARY_alpha_VELY(:,:,:) )
    endif

    if ( ATMOS_BOUNDARY_USE_PT ) then
       call FILE_CARTESC_read( fid, 'PT', 'ZXY', ATMOS_BOUNDARY_POTT(:,:,:) )
       call FILE_CARTESC_read( fid, 'ALPHA_PT', 'ZXY', ATMOS_BOUNDARY_alpha_POTT(:,:,:) )
    endif

    do iq = 1, QA
       iqb = BND_IQ(iq)
       if ( iqb > 0 ) then
          call FILE_CARTESC_read( fid, TRACER_NAME(iq), 'ZXY', ATMOS_BOUNDARY_QTRC(:,:,:,iqb) )
          call FILE_CARTESC_read( fid, 'ALPHA_'//trim(TRACER_NAME(iq)), 'ZXY', ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iqb) )
       endif
    end do

    call FILE_CARTESC_close( fid )


    call ATMOS_BOUNDARY_var_fillhalo
    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_read

  !-----------------------------------------------------------------------------
  !> Write boundary data
  subroutine ATMOS_BOUNDARY_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_create,    &
       FILE_CARTESC_def_var,   &
       FILE_CARTESC_enddef,    &
       FILE_CARTESC_write_var, &
       FILE_CARTESC_close
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ
    implicit none

    integer :: fid
    integer :: vid_dens, vid_a_dens
    integer :: vid_velz, vid_a_velz
    integer :: vid_velx, vid_a_velx
    integer :: vid_vely, vid_a_vely
    integer :: vid_pott, vid_a_pott
    integer :: vid_qtrc(BND_QA), vid_a_qtrc(BND_QA)
    integer :: iq, iqb
    !---------------------------------------------------------------------------

    call FILE_CARTESC_create( ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, & ! [IN]
                              ATMOS_BOUNDARY_OUT_DTYPE,                              & ! [IN]
                              fid,                                                   & ! [OUT]
                              aggregate = ATMOS_BOUNDARY_OUT_AGGREGATE               ) ! [IN]

    if (      ATMOS_BOUNDARY_USE_DENS &
         .OR. ATMOS_BOUNDARY_USE_VELZ &
         .OR. ATMOS_BOUNDARY_USE_VELX &
         .OR. ATMOS_BOUNDARY_USE_VELY &
         .OR. ATMOS_BOUNDARY_USE_PT &
         .OR. l_bnd                   ) then
       call FILE_CARTESC_def_var( fid, 'DENS', 'Reference Density', 'kg/m3', 'ZXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_dens )
    else
       vid_dens = -1
    end if

    if ( ATMOS_BOUNDARY_USE_DENS .OR. l_bnd ) then
       call FILE_CARTESC_def_var( fid, 'ALPHA_DENS', 'Alpha for DENS', '1', 'ZXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_a_dens )
    else
       vid_a_dens = -1
    end if

    if ( ATMOS_BOUNDARY_USE_VELZ .OR. (l_bnd .AND. ONLINE_USE_VELZ) ) then
       call FILE_CARTESC_def_var( fid, 'VELZ', 'Reference Velocity w', 'm/s', 'ZHXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_velz )
       call FILE_CARTESC_def_var( fid, 'ALPHA_VELZ', 'Alpha for VELZ', '1', 'ZHXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_a_velz )
    else
       vid_velz = -1
    end if

    if ( ATMOS_BOUNDARY_USE_VELX .OR. l_bnd ) then
       call FILE_CARTESC_def_var( fid, 'VELX', 'Reference Velocity u', 'm/s', 'ZXHY', ATMOS_BOUNDARY_OUT_DTYPE, vid_velx )
       call FILE_CARTESC_def_var( fid, 'ALPHA_VELX', 'Alpha for VELX', '1', 'ZXHY', ATMOS_BOUNDARY_OUT_DTYPE, vid_a_velx )
    else
       vid_velx = -1
    end if

    if ( ATMOS_BOUNDARY_USE_VELY .OR. l_bnd ) then
       call FILE_CARTESC_def_var( fid, 'VELY', 'Reference Velocity y', 'm/s', 'ZXYH', ATMOS_BOUNDARY_OUT_DTYPE, vid_vely )
       call FILE_CARTESC_def_var( fid, 'ALPHA_VELY', 'Alpha for VELY', '1', 'ZXYH', ATMOS_BOUNDARY_OUT_DTYPE, vid_a_vely )
    else
       vid_vely = -1
    end if

    if ( ATMOS_BOUNDARY_USE_PT .OR. l_bnd ) then
       call FILE_CARTESC_def_var( fid, 'PT', 'Reference PT', 'K', 'ZXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_pott )
       call FILE_CARTESC_def_var( fid, 'ALPHA_PT', 'Alpha for PT', '1', 'ZXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_a_pott )
    else
       vid_pott = -1
    end if

    do iq = 1, QA
       iqb = BND_IQ(iq)
       if ( iqb > 0 ) then
          call FILE_CARTESC_def_var( fid, TRACER_NAME(iq), 'Reference '//trim(TRACER_NAME(iq)), TRACER_UNIT(iq), 'ZXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_qtrc(iqb) )
          call FILE_CARTESC_def_var( fid, 'ALPHA_'//trim(TRACER_NAME(iq)), 'Alpha for '//trim(TRACER_NAME(iq)), '1', 'ZXY', ATMOS_BOUNDARY_OUT_DTYPE, vid_a_qtrc(iqb) )
       end if
    end do


    call FILE_CARTESC_enddef( fid )


    if ( vid_dens > 0 ) then
       call FILE_CARTESC_write_var( fid, vid_dens, ATMOS_BOUNDARY_DENS(:,:,:), 'DENS', 'ZXY' )
    end if

    if ( vid_a_dens > 0 ) then
       call FILE_CARTESC_write_var( fid, vid_a_dens, ATMOS_BOUNDARY_alpha_DENS(:,:,:), 'ALPHA_DENS', 'ZXY' )
    end if

    if ( vid_velz > 0 ) then
       call FILE_CARTESC_write_var( fid, vid_velz, ATMOS_BOUNDARY_VELZ(:,:,:), 'VELZ', 'ZHXY' )
       call FILE_CARTESC_write_var( fid, vid_a_velz, ATMOS_BOUNDARY_alpha_VELZ(:,:,:), 'ALPHA_VELZ', 'ZHXY' )
    end if

    if ( vid_velx > 0 ) then
       call FILE_CARTESC_write_var( fid, vid_velx, ATMOS_BOUNDARY_VELX(:,:,:), 'VELX', 'ZXHY' )
       call FILE_CARTESC_write_var( fid, vid_a_velx, ATMOS_BOUNDARY_alpha_VELX(:,:,:), 'ALPHA_VELX', 'ZXHY' )
    endif

    if ( vid_vely > 0 ) then
       call FILE_CARTESC_write_var( fid, vid_vely, ATMOS_BOUNDARY_VELY(:,:,:), 'VELY', 'ZXYH' )
       call FILE_CARTESC_write_var( fid, vid_a_vely, ATMOS_BOUNDARY_alpha_VELY(:,:,:), 'ALPHA_VELY', 'ZXYH' )
    end if

    if ( vid_pott > 0 ) then
       call FILE_CARTESC_write_var( fid, vid_pott, ATMOS_BOUNDARY_POTT(:,:,:), 'PT', 'ZXY' )
       call FILE_CARTESC_write_var( fid, vid_a_pott, ATMOS_BOUNDARY_alpha_POTT(:,:,:), 'ALPHA_PT', 'ZXY' )
    end if

    do iqb = 1, BND_QA
       if ( vid_qtrc(iqb) > 0 ) then
          call FILE_CARTESC_write_var( fid, vid_qtrc(iqb), ATMOS_BOUNDARY_QTRC(:,:,:,iqb), TRACER_NAME(iqb), 'ZXY' )
          call FILE_CARTESC_write_var( fid, vid_a_qtrc(iqb), ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iqb), 'ALPHA_'//trim(TRACER_NAME(iqb)), 'ZXY' )
       end if
    end do


    call FILE_CARTESC_close( fid )


    return
  end subroutine ATMOS_BOUNDARY_write

  !-----------------------------------------------------------------------------
  !> generate boundary data
  subroutine ATMOS_BOUNDARY_generate
    use scale_atmos_refstate, only: &
         ATMOS_REFSTATE_dens
    implicit none

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_REFSTATE_DENS(k,i,j)
       ATMOS_BOUNDARY_VELZ(k,i,j) = ATMOS_BOUNDARY_VALUE_VELZ
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_VALUE_VELX
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_VALUE_VELY
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_VALUE_PT
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_VALUE_QTRC
       end do
    enddo
    enddo
    enddo

    call ATMOS_BOUNDARY_var_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_generate

  !-----------------------------------------------------------------------------
  !> Initialize boundary value for real case experiment
  subroutine ATMOS_BOUNDARY_initialize_file
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec, &
       CALENDAR_date2char
    use scale_time, only: &
       TIME_NOWDATE
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_check_coordinates
    implicit none

    integer           :: boundary_time_startday
    real(DP)          :: boundary_time_startsec
    real(DP)          :: boundary_time_startms
    integer           :: boundary_time_offset_year
    character(len=27) :: boundary_chardate

    if ( ATMOS_BOUNDARY_START_DATE(1) == -9999 ) then
       ATMOS_BOUNDARY_START_DATE = TIME_NOWDATE
    end if

    !--- calculate time of the initial step in boundary file [no offset]
    boundary_time_startms     = 0.0_DP
    boundary_time_offset_year = 0
    call CALENDAR_date2char( boundary_chardate,            & ! [OUT]
                             ATMOS_BOUNDARY_START_DATE(:), & ! [IN]
                             boundary_time_startms         ) ! [IN]

    call CALENDAR_date2daysec( boundary_time_startday,       & ! [OUT]
                               boundary_time_startsec,       & ! [OUT]
                               ATMOS_BOUNDARY_START_DATE(:), & ! [IN]
                               boundary_time_startms,        & ! [IN]
                               boundary_time_offset_year     ) ! [IN]

    boundary_time_initdaysec = CALENDAR_combine_daysec( boundary_time_startday, boundary_time_startsec )

    LOG_INFO("ATMOS_BOUNDARY_initialize_file",'(1x,A,A)') 'BOUNDARY START Date     : ', boundary_chardate

    call FILE_CARTESC_open( ATMOS_BOUNDARY_IN_BASENAME, ATMOS_BOUNDARY_fid )

    if ( ATMOS_BOUNDARY_IN_CHECK_COORDINATES ) then
       call FILE_CARTESC_check_coordinates( ATMOS_BOUNDARY_fid, atmos=.true. )
    end if

    return
  end subroutine ATMOS_BOUNDARY_initialize_file

  !-----------------------------------------------------------------------------
  !> Set boundary value for real case experiment
  subroutine ATMOS_BOUNDARY_set_file
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI => CONST_PI
    use scale_time, only: &
       TIME_NOWDATE,      &
       TIME_DTSEC
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec
    implicit none

    real(RP) :: bnd_DENS(KA,IA,JA)        ! damping coefficient for DENS (0-1)
    real(RP) :: bnd_VELZ(KA,IA,JA)        ! damping coefficient for VELZ (0-1)
    real(RP) :: bnd_VELX(KA,IA,JA)        ! damping coefficient for VELX (0-1)
    real(RP) :: bnd_VELY(KA,IA,JA)        ! damping coefficient for VELY (0-1)
    real(RP) :: bnd_POTT(KA,IA,JA)        ! damping coefficient for POTT (0-1)
    real(RP) :: bnd_QTRC(KA,IA,JA,BND_QA) ! damping coefficient for QTRC (0-1)

    integer  :: run_time_startdate(6)
    integer  :: run_time_startday
    real(DP) :: run_time_startsec
    real(DP) :: run_time_startms
    integer  :: run_time_offset_year
    real(DP) :: run_time_nowdaysec

    real(DP) :: boundary_diff_daysec
    real(RP) :: boundary_inc_offset
    integer  :: fillgaps_steps

    real(RP) :: total

    integer  :: i, j, k, n, iq
    !---------------------------------------------------------------------------

    if ( ATMOS_BOUNDARY_UPDATE_DT <= 0.0_DP ) then
       LOG_ERROR("ATMOS_BOUNDARY_set_file",*) 'You need specify ATMOS_BOUNDARY_UPDATE_DT as larger than 0.0'
       call PRC_abort
    endif
    UPDATE_NSTEP = nint( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
    if ( abs(UPDATE_NSTEP * TIME_DTSEC - ATMOS_BOUNDARY_UPDATE_DT) > 1E-10_DP ) then
       LOG_ERROR("ATMOS_BOUNDARY_set_file",*) 'ATMOS_BOUNDARY_UPDATE_DT is not multiple of DT'
       call PRC_abort
    end if

    !--- recalculate time of the run [no offset]
    run_time_startdate(:) = TIME_NOWDATE(:)
    run_time_startms      = 0.0_DP
    run_time_offset_year  = 0

    call CALENDAR_date2daysec( run_time_startday,     & ! [OUT]
                               run_time_startsec,     & ! [OUT]
                               run_time_startdate(:), & ! [IN]
                               run_time_startms,      & ! [IN]
                               run_time_offset_year   ) ! [IN]

    run_time_nowdaysec = CALENDAR_combine_daysec( run_time_startday, run_time_startsec )

    boundary_diff_daysec = run_time_nowdaysec - boundary_time_initdaysec
    boundary_timestep    = 1 + int( boundary_diff_daysec / ATMOS_BOUNDARY_UPDATE_DT )
    boundary_inc_offset  = mod( boundary_diff_daysec, ATMOS_BOUNDARY_UPDATE_DT )
    fillgaps_steps       = int( boundary_inc_offset / TIME_DTSEC )

    LOG_INFO("ATMOS_BOUNDARY_set_file",*) 'BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep
    LOG_INFO("ATMOS_BOUNDARY_set_file",*) 'BOUNDARY OFFSET:', boundary_inc_offset
    LOG_INFO("ATMOS_BOUNDARY_set_file",*) 'BOUNDARY FILLGAPS STEPS:', fillgaps_steps


    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       allocate( OFFSET_TIME_FACT(0:UPDATE_NSTEP) )
       total = 0.0_RP
       !$omp parallel do reduction(+:total)
       do n = 0, UPDATE_NSTEP
          OFFSET_TIME_FACT(n) = 1.0_RP - cos( 2.0_RP * PI * ( n - 1 ) / UPDATE_NSTEP )
          total = total + OFFSET_TIME_FACT(n)
       end do
       total = total / UPDATE_NSTEP
       !$omp parallel do
       do n = 0, UPDATE_NSTEP
          OFFSET_TIME_FACT(n) = OFFSET_TIME_FACT(n) / total
       end do
    end if

    ! read boundary data from input file
    call ATMOS_BOUNDARY_update_file( ref_now )

    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       call calc_mass( ref_now )
    end if

    boundary_timestep = boundary_timestep + 1
    call ATMOS_BOUNDARY_update_file( ref_new )

    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       call calc_mass( ref_new )
       call set_offset
    end if

    ! copy now to old
    !$omp parallel do default(none) private(i,j,k,iq) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KS,KE,ATMOS_BOUNDARY_ref_DENS,ref_old,ref_now,ATMOS_BOUNDARY_ref_VELX) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELY,ATMOS_BOUNDARY_ref_POTT,BND_QA,ATMOS_BOUNDARY_ref_QTRC)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now)
       ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now)
       ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now)
       ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_old) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now)
       end do
    end do
    end do
    end do

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ATMOS_BOUNDARY_ref_VELZ(k,i,j,:) = ATMOS_BOUNDARY_VALUE_VELZ
       end do
       end do
       end do
    end if

    now_step = fillgaps_steps

    ! set boundary data
    call set_boundary( ATMOS_BOUNDARY_USE_VELZ )

    return
  end subroutine ATMOS_BOUNDARY_set_file

  !-----------------------------------------------------------------------------
  !> Initialize boundary value for real case experiment [online daughter]
  subroutine ATMOS_BOUNDARY_initialize_online
    use scale_prc, only: &
       PRC_abort
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_recvwait_issue, &
       ONLINE_USE_VELZ,          &
       PARENT_DTSEC,             &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    implicit none

    ! parameters
    integer, parameter  :: handle = 2

    ATMOS_BOUNDARY_UPDATE_DT = PARENT_DTSEC(handle)

    if ( NESTQA > BND_QA ) then
       LOG_ERROR("ATMOS_BOUNDARY_initialize_online",*) 'NEST_BND_QA exceeds BND_QA'
       LOG_ERROR_CONT(*) 'This must not be occur.'
       LOG_ERROR_CONT(*) 'Please send your configuration file to SCALE develop member.'
       call PRC_abort
    end if

    call COMM_CARTESC_NEST_recvwait_issue( handle, NESTQA )

    return
  end subroutine ATMOS_BOUNDARY_initialize_online

  !-----------------------------------------------------------------------------
  !> Set boundary value for real case experiment [online daughter]
  subroutine ATMOS_BOUNDARY_set_online
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI => CONST_PI
    use scale_time, only: &
       TIME_DTSEC,        &
       TIME_NSTEP
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ,          &
       PARENT_NSTEP
    implicit none

    ! parameters
    integer, parameter  :: handle = 2

    real(RP) :: total

    ! works
    integer  :: i, j, k, n, iq
    !---------------------------------------------------------------------------

    ! import data from parent domain
    boundary_timestep = 1
    LOG_INFO("ATMOS_BOUNDARY_set_online",*) 'BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep

    call ATMOS_BOUNDARY_update_online_daughter( ref_now )

    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       call calc_mass( ref_now )
    end if

    boundary_timestep = boundary_timestep + 1
    LOG_INFO("ATMOS_BOUNDARY_set_online",*) 'BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep

    call ATMOS_BOUNDARY_update_online_daughter( ref_new )

    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       call calc_mass( ref_new )
       call set_offset
    end if

    ! copy now to old
    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now)
       ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now)
       ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now)
       ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_old) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_old) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now)
       end do
    end do
    end do
    end do

    UPDATE_NSTEP = nint( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
    if ( UPDATE_NSTEP * TIME_DTSEC /= ATMOS_BOUNDARY_UPDATE_DT ) then
       LOG_ERROR("ATMOS_BOUNDARY_set_online",*) 'DT of the parent is not multiple of the DT'
       call PRC_abort
    end if
    if ( UPDATE_NSTEP * PARENT_NSTEP(handle) /= TIME_NSTEP ) then
       LOG_ERROR("ATMOS_BOUNDARY_set_online",*) 'DURATION must be the same as that of the parent'
       call PRC_abort
    end if

    now_step = 0 ! should be set as zero in initialize process

    if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
       allocate( OFFSET_TIME_FACT(UPDATE_NSTEP) )
       total = 0.0_RP
       !$omp parallel do reduction(+:total)
       do n = 0, UPDATE_NSTEP
          OFFSET_TIME_FACT(n) = 1.0_RP - cos( 2.0_RP * PI * ( n - 1 ) / UPDATE_NSTEP )
          total = total + OFFSET_TIME_FACT(n)
       end do
       total = total / UPDATE_NSTEP
       !$omp parallel do
       do n = 0, UPDATE_NSTEP
          OFFSET_TIME_FACT(n) = OFFSET_TIME_FACT(n) / total
       end do
    end if

    ! set boundary data
    call set_boundary( ONLINE_USE_VELZ )

    return
  end subroutine ATMOS_BOUNDARY_set_online

  !-----------------------------------------------------------------------------
  !> First send boundary value
  subroutine ATMOS_BOUNDARY_firstsend( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, QV, Qe )
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    implicit none

    ! arguments
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(in) :: QV  (KA,IA,JA)
    real(RP), intent(in) :: Qe  (KA,IA,JA,N_HYD)
    !---------------------------------------------------------------------------

    ! send data at the first time
    if ( do_parent_process ) then !online [parent]
       ! issue send
       call ATMOS_BOUNDARY_send( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, QV, Qe )
    endif

    return
  end subroutine ATMOS_BOUNDARY_firstsend

  !-----------------------------------------------------------------------------
  !> Finalize boundary value
  subroutine ATMOS_BOUNDARY_driver_finalize
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_recvwait_issue, &
       COMM_CARTESC_NEST_recv_cancel, &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    use scale_file_cartesC, only: &
       FILE_CARTESC_close
    implicit none

    ! works
    integer :: handle
    !---------------------------------------------------------------------------

    if ( do_parent_process ) then !online [parent]
       handle = 1
       call COMM_CARTESC_NEST_recvwait_issue( handle, NESTQA )
    endif

    if ( do_daughter_process ) then !online [daughter]
       handle = 2
       call COMM_CARTESC_NEST_recv_cancel( handle )
    endif

    if ( ATMOS_BOUNDARY_fid > 0 ) then
       call FILE_CARTESC_close( ATMOS_BOUNDARY_fid )
       ATMOS_BOUNDARY_fid = -1
    end if

    return
  end subroutine ATMOS_BOUNDARY_driver_finalize

  !-----------------------------------------------------------------------------
  !> Update boundary value with a constant time boundary
  subroutine ATMOS_BOUNDARY_driver_update( &
       last_step )
    use scale_prc, only: &
       PRC_abort
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ,       &
       COMM_CARTESC_NEST_test
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC, &
       QV,   &
       Qe
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    implicit none

    logical, intent(in) :: last_step

    integer :: handle
    !---------------------------------------------------------------------------

    if ( do_parent_process ) then !online [parent]
       ! should be called every time step
       call ATMOS_BOUNDARY_update_online_parent( DENS,MOMZ,MOMX,MOMY,RHOT,QTRC(:,:,:,QS_MP:QE_MP), QV, Qe )
    endif

    if ( l_bnd ) then

       ! step boundary
       now_step = now_step + 1

       ! update referce vars
       if ( last_step ) then
          now_step = min( now_step, UPDATE_NSTEP-1 )
       else if ( now_step >= UPDATE_NSTEP ) then
          now_step          = 0
          boundary_timestep = boundary_timestep + 1

          call update_ref_index

          if ( do_daughter_process ) then !online [daughter]
             call ATMOS_BOUNDARY_update_online_daughter( ref_new )
          else
             call ATMOS_BOUNDARY_update_file( ref_new )
          end if

          if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
             call calc_mass( ref_new )
          end if

       end if

       call set_boundary( ONLINE_USE_VELZ )

       if ( ATMOS_BOUNDARY_DENS_ADJUST ) then
          call set_offset
       end if

    elseif ( do_parent_process ) then
       ! do nothing
    else
       LOG_ERROR("ATMOS_BOUNDARY_update",*) '[BUG] invalid path'
       call PRC_abort
    end if

    call history_bnd( ATMOS_BOUNDARY_DENS, &
                      ATMOS_BOUNDARY_VELZ, &
                      ATMOS_BOUNDARY_VELX, &
                      ATMOS_BOUNDARY_VELY, &
                      ATMOS_BOUNDARY_POTT, &
                      ATMOS_BOUNDARY_QTRC )

    ! To be enable to do asynchronous communicaton
    if ( do_parent_process ) then !online [parent]
       handle = 1
       call COMM_CARTESC_NEST_test( handle )
    endif
    if ( do_daughter_process ) then !online [daughter]
       handle = 2
       call COMM_CARTESC_NEST_test( handle )
    endif

    return
  end subroutine ATMOS_BOUNDARY_driver_update

  !-----------------------------------------------------------------------------
  !> Update reference boundary from file
  subroutine ATMOS_BOUNDARY_update_file( ref )
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    use mod_atmos_phy_mp_vars, only: &
       QS_MP
    implicit none

    integer, intent(in) :: ref

    integer :: fid, iq, iqb
    !---------------------------------------------------------------------------

    LOG_INFO("ATMOS_BOUNDARY_update_file",*) "Atmos Boundary: read from boundary file(timestep=", boundary_timestep, ")"

    fid = ATMOS_BOUNDARY_fid

    call FILE_CARTESC_read( fid, 'DENS', 'ZXY',  ATMOS_BOUNDARY_ref_DENS(:,:,:,ref), step=boundary_timestep )
    call FILE_CARTESC_read( fid, 'VELX', 'ZXHY', ATMOS_BOUNDARY_ref_VELX(:,:,:,ref), step=boundary_timestep )
    call FILE_CARTESC_read( fid, 'VELY', 'ZXYH', ATMOS_BOUNDARY_ref_VELY(:,:,:,ref), step=boundary_timestep )
    call FILE_CARTESC_read( fid, 'PT',   'ZXY',  ATMOS_BOUNDARY_ref_POTT(:,:,:,ref), step=boundary_timestep )
    do iq = 1, QA
       iqb = BND_IQ(iq)
       if ( iqb > 0 ) then
          call FILE_CARTESC_read( fid, TRACER_NAME(iq), 'ZXY', ATMOS_BOUNDARY_ref_QTRC(:,:,:,iqb,ref), step=boundary_timestep )
       end if
    end do

    call FILE_CARTESC_flush( fid )

    ! fill HALO in reference
    call ATMOS_BOUNDARY_ref_fillhalo( ref )

    return
  end subroutine ATMOS_BOUNDARY_update_file

  !-----------------------------------------------------------------------------
  !> Send reference boundary value to daughter domain by communicate
  subroutine ATMOS_BOUNDARY_update_online_parent( &
       DENS, & ! [in]
       MOMZ, & ! [in]
       MOMX, & ! [in]
       MOMY, & ! [in]
       RHOT, & ! [in]
       QTRC, & ! [in]
       QV,   & ! [in]
       Qe    ) ! [in]
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_recvwait_issue, &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    implicit none

    ! arguments
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(in) :: QV  (KA,IA,JA)
    real(RP), intent(in) :: Qe  (KA,IA,JA,N_HYD)

    integer, parameter :: handle = 1
    !---------------------------------------------------------------------------

    LOG_INFO("ATMOS_BOUNDARY_update_online_parent",*)"ATMOS BOUNDARY update online: PARENT"

    ! issue wait
    call COMM_CARTESC_NEST_recvwait_issue( handle, NESTQA )

    ! issue send
    call ATMOS_BOUNDARY_send( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, QV, Qe )

    return
  end subroutine ATMOS_BOUNDARY_update_online_parent

  !-----------------------------------------------------------------------------
  !> Update reference boundary by communicate with parent domain
  subroutine ATMOS_BOUNDARY_update_online_daughter( &
       ref ) ! [in]
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_recvwait_issue, &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    implicit none

    ! arguments
    integer,  intent(in) :: ref

    integer, parameter :: handle = 2
    !---------------------------------------------------------------------------

    LOG_INFO("ATMOS_BOUNDARY_update_online_daughter",'(1X,A,I5)') 'ATMOS BOUNDARY update online: DAUGHTER', boundary_timestep

    ! issue wait
    call ATMOS_BOUNDARY_recv( ref )

    ! fill HALO in reference
    call ATMOS_BOUNDARY_ref_fillhalo( ref )

    ! issue receive
    call COMM_CARTESC_NEST_recvwait_issue( handle, NESTQA )

    return
  end subroutine ATMOS_BOUNDARY_update_online_daughter

  !-----------------------------------------------------------------------------
  !> Send boundary value
  subroutine ATMOS_BOUNDARY_send( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, QV, Qe )
    use scale_comm_cartesC_nest, only: &
       ONLINE_BOUNDARY_DIAGQHYD,   &
       COMM_CARTESC_NEST_nestdown,    &
       DAUGHTER_KA,           &
       DAUGHTER_IA,           &
       DAUGHTER_JA,           &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    implicit none

    ! parameters
    integer, parameter  :: handle = 1

    ! arguments
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in), target :: QTRC(KA,IA,JA,QA_MP)
    real(RP), intent(in) :: QV  (KA,IA,JA)
    real(RP), intent(in) :: Qe  (KA,IA,JA,N_HYD)

    ! works
    real(RP), pointer :: Q(:,:,:,:)
    real(RP) :: dummy_d( DAUGHTER_KA(handle), DAUGHTER_IA(handle), DAUGHTER_JA(handle), NESTQA )

    integer :: iq
    !---------------------------------------------------------------------------

    if ( ONLINE_BOUNDARY_DIAGQHYD ) then
       Q => Q_WORK
       Q(:,:,:,1) = QV(:,:,:)
       do iq = 2, NESTQA
          Q(:,:,:,iq) = Qe(:,:,:,iq-1)
       end do
    else
       Q => QTRC
    end if

    call COMM_CARTESC_NEST_nestdown( handle,           &
                                     NESTQA,           &
                                     DENS(:,:,:),      &  !(KA,IA,JA)
                                     MOMZ(:,:,:),      &  !(KA,IA,JA)
                                     MOMX(:,:,:),      &  !(KA,IA,JA)
                                     MOMY(:,:,:),      &  !(KA,IA,JA)
                                     RHOT(:,:,:),      &  !(KA,IA,JA)
                                     Q   (:,:,:,:),    &  !(KA,IA,JA,NESTQA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,:)  )  !(KA,IA,JA,NESTQA)

    return
  end subroutine ATMOS_BOUNDARY_send

  !-----------------------------------------------------------------------------
  !> Recieve boundary value
  subroutine ATMOS_BOUNDARY_recv( &
       ref_idx )
    use scale_prc, only: &
       PRC_abort
    use scale_comm_cartesC_nest, only: &
       ONLINE_BOUNDARY_DIAGQHYD,   &
       COMM_CARTESC_NEST_nestdown, &
       PARENT_KA,                  &
       PARENT_IA,                  &
       PARENT_JA,                  &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
    implicit none

    ! parameters
    integer, parameter  :: handle = 2

    ! arguments
    integer, intent(in) :: ref_idx

    ! works
    real(RP), pointer :: Q(:,:,:,:)
    real(RP) :: dummy_p( PARENT_KA(handle), PARENT_IA(handle), PARENT_JA(handle), NESTQA )
    !---------------------------------------------------------------------------

!OCL XFILL
    dummy_p(:,:,:,:) = 0.0_RP

    if ( ONLINE_BOUNDARY_DIAGQHYD ) then
       Q => Q_WORK
    else
       Q => ATMOS_BOUNDARY_ref_QTRC(:,:,:,1:NESTQA,ref_idx)
    end if

    call COMM_CARTESC_NEST_nestdown( handle,                                 &
                                     NESTQA,                                 &
                                     dummy_p(:,:,:,1),                       & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                       & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                       & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                       & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                       & !(KA,IA,JA)
                                     dummy_p(:,:,:,1:NESTQA),                & !(KA,IA,JA,NESTQA)
                                     ATMOS_BOUNDARY_ref_DENS(:,:,:,ref_idx), & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_VELZ(:,:,:,ref_idx), & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_VELX(:,:,:,ref_idx), & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_VELY(:,:,:,ref_idx), & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_POTT(:,:,:,ref_idx), & !(KA,IA,JA)
                                     Q(:,:,:,:)                              ) !(KA,IA,JA,NESTQA)


    if ( ONLINE_BOUNDARY_DIAGQHYD ) then
       call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, 1, IA, JA, 1, JA, &
                                           Q(:,:,:,1), Q(:,:,:,2:),                 & ! (in)
                                           ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,ref_idx) ) ! (out)
    end if

    return
  end subroutine ATMOS_BOUNDARY_recv

  !-----------------------------------------------------------------------------
  subroutine set_boundary( use_velz )
    use scale_prc_cartesC, only: &
       PRC_HAS_W,   &
       PRC_HAS_E,   &
       PRC_HAS_S,   &
       PRC_HAS_N
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none
    logical, intent(in) :: use_velz

    real(RP) :: bnd_DENS(KA,IA,JA)        ! damping coefficient for DENS (0-1)
    real(RP) :: bnd_VELZ(KA,IA,JA)        ! damping coefficient for VELZ (0-1)
    real(RP) :: bnd_VELX(KA,IA,JA)        ! damping coefficient for VELX (0-1)
    real(RP) :: bnd_VELY(KA,IA,JA)        ! damping coefficient for VELY (0-1)
    real(RP) :: bnd_POTT(KA,IA,JA)        ! damping coefficient for POTT (0-1)
    real(RP) :: bnd_QTRC(KA,IA,JA,BND_QA) ! damping coefficient for QTRC (0-1)

    integer :: i, j, k, iq, iqb


    ! get boundaryal coefficients
    call get_boundary( bnd_DENS(:,:,:),   & ! [OUT]
                       bnd_VELZ(:,:,:),   & ! [OUT]
                       bnd_VELX(:,:,:),   & ! [OUT]
                       bnd_VELY(:,:,:),   & ! [OUT]
                       bnd_POTT(:,:,:),   & ! [OUT]
                       bnd_QTRC(:,:,:,:), & ! [OUT]
                       now_step,          & ! [IN]
                       UPDATE_NSTEP       ) ! [IN]

    ! update boundary vars
    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_DENS(k,i,j) = bnd_DENS(k,i,j)
       ATMOS_BOUNDARY_VELX(k,i,j) = bnd_VELX(k,i,j)
       ATMOS_BOUNDARY_VELY(k,i,j) = bnd_VELY(k,i,j)
       ATMOS_BOUNDARY_POTT(k,i,j) = bnd_POTT(k,i,j)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = bnd_QTRC(k,i,j,iq)
       end do
    end do
    end do
    end do
    if ( USE_VELZ ) then
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ATMOS_BOUNDARY_VELZ(k,i,j) = bnd_VELZ(k,i,j)
       end do
       end do
       end do
    end if

    ! fill HALO in western region
    if ( .NOT. PRC_HAS_W ) then
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JA,JS,IS,KS,KE,QA,DENS,MOMX,RHOT,QTRC) &
       !$omp shared(ATMOS_BOUNDARY_DENS,ATMOS_BOUNDARY_VELX) &
       !$omp shared(ATMOS_BOUNDARY_POTT,ATMOS_BOUNDARY_QTRC) &
       !$omp shared(BND_QA,BND_IQ) &
       !$omp private(i,j,k,iq,iqb)
       do j = 1, JA
       do i = 1, IS-1
       do k = KS, KE
          DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
          MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
          RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
          do iq = 1, QA
             iqb = BND_IQ(iq)
             if ( iqb > 0 ) then
                QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iqb)
             else
                QTRC(k,i,j,iq) = QTRC(k,IS,j,iq)
             end if
          end do
       end do
       end do
       end do
       !$omp parallel do
       do j = 1, JA-1
       do i = 1, IS-1
       do k = KS, KE
          MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IS-1
       do k = KS, KE
          MOMY(k,i,JA) = ATMOS_BOUNDARY_VELY(k,i,JA) &
                  * ATMOS_BOUNDARY_DENS(k,i,JA)
       end do
       end do
       if ( USE_VELZ ) then
          !$omp parallel do
          do j = 1, JA
          do i = 1, IS-1
          do k = KS, KE-1
             MOMZ(k,i,j) = ATMOS_BOUNDARY_VELZ(k,i,j) &
                     * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k+1,i,j) ) * 0.5_RP
          end do
          end do
          end do
       else
          !$omp parallel do
          do j = 1, JA
          do i = 1, IS-1
          do k = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,IS,j)
          end do
          end do
          end do
       end if
    end if

    ! fill HALO in eastern region
    if ( .NOT. PRC_HAS_E ) then
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JA,IE,IA,KS,KE,QA) &
       !$omp shared(DENS,RHOT,QTRC) &
       !$omp shared(ATMOS_BOUNDARY_DENS,ATMOS_BOUNDARY_VELX) &
       !$omp shared(ATMOS_BOUNDARY_POTT,ATMOS_BOUNDARY_QTRC) &
       !$omp shared(BND_QA,BND_IQ) &
       !$omp private(i,j,k,iq,iqb)
       do j = 1, JA
       do i = IE+1, IA
       do k = KS, KE
          DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
          RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
          do iq = 1, QA
             iqb = BND_IQ(iq)
             if ( iqb > 0 ) then
                QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iqb)
             else
                QTRC(k,i,j,iq) = QTRC(k,IE,j,iq)
             end if
          end do
       end do
       end do
       end do
       !$omp parallel do
       do j = 1, JA
       do i = IE, IA-1
       do k = KS, KE
          MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
       end do
       end do
       end do
       !$omp parallel do
       do j = 1, JA
       do k = KS, KE
          MOMX(k,IA,j) = ATMOS_BOUNDARY_VELX(k,IA,j) * ATMOS_BOUNDARY_DENS(k,IA,j)
       end do
       end do
       !$omp parallel do
       do j = 1, JA-1
       do i = IE+1, IA
       do k = KS, KE
          MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
       end do
       end do
       end do
       do i = IE+1, IA
       do k = KS, KE
          MOMY(k,i,JA) = ATMOS_BOUNDARY_VELY(k,i,JA) &
                  * ATMOS_BOUNDARY_DENS(k,i,JA)
       end do
       end do
       if ( USE_VELZ ) then
          !$omp parallel do
          do j = 1, JA
          do i = IE+1, IA
          do k = KS, KE-1
             MOMZ(k,i,j) = ATMOS_BOUNDARY_VELZ(k,i,j) &
                     * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k+1,i,j) ) * 0.5_RP
          end do
          end do
          end do
       else
          !$omp parallel do
          do j = 1, JA
          do i = IE+1, IA
          do k = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,IE,j)
          end do
          end do
          end do
       end if
    end if

    ! fill HALO in southern region
    if ( .NOT. PRC_HAS_S ) then
       do j = 1, JS-1
       do i = 1, IA
       do k = KS, KE
          DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
          MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
          RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
          do iq = 1, QA
             iqb = BND_IQ(iq)
             if ( iqb > 0 ) then
                QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iqb)
             else
                QTRC(k,i,j,iq) = QTRC(k,i,JS,iq)
             end if
          end do
       end do
       end do
       end do
       do j = 1, JS-1
       do i = 1, IA-1
       do k = KS, KE
          MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
       end do
       end do
       end do
       do j = 1, JS-1
       do k = KS, KE
          MOMX(k,IA,j) = ATMOS_BOUNDARY_VELX(k,IA,j) &
                  * ATMOS_BOUNDARY_DENS(k,IA,j)
       end do
       end do
       if ( USE_VELZ ) then
          do j = 1, JS-1
          do i = 1, IA
          do k = KS, KE-1
             MOMZ(k,i,j) = ATMOS_BOUNDARY_VELZ(k,i,j) &
                     * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k+1,i,j) ) * 0.5_RP
          end do
          end do
          end do
       else
          do j = 1, JS-1
          do i = 1, IA
          do k = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,i,JS)
          end do
          end do
          end do
       end if
    end if

    ! fill HALO in northern region
    if ( .NOT. PRC_HAS_N ) then
       do j = JE+1, JA
       do i = 1, IA
       do k = KS, KE
          DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
          RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
          do iq = 1, QA
             iqb = BND_IQ(iq)
             if ( iqb > 0 ) then
                QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iqb)
             else
                QTRC(k,i,j,iq) = QTRC(k,i,JE,iq)
             end if
          end do
       end do
       end do
       end do
       do j = JE, JA-1
       do i = 1, IA
       do k = KS, KE
          MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
       end do
       end do
       end do
       do i = 1, IA
       do k = KS, KE
          MOMY(k,i,JA) = ATMOS_BOUNDARY_VELY(k,i,JA) * ATMOS_BOUNDARY_DENS(k,i,JA)
       end do
       end do
       do j = JE+1, JA
       do i = 1, IA-1
       do k = KS, KE
          MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
       end do
       end do
       end do
       do j = JE+1, JA
       do k = KS, KE
          MOMX(k,IA,j) = ATMOS_BOUNDARY_VELX(k,IA,j) &
                  * ATMOS_BOUNDARY_DENS(k,IA,j)
       end do
       end do
       if ( USE_VELZ ) then
          do j = JE+1, JA
          do i = 1, IA
          do k = KS, KE-1
             MOMZ(k,i,j) = ATMOS_BOUNDARY_VELZ(k,i,j) &
                     * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k+1,i,j) ) * 0.5_RP
          end do
          end do
          end do
       else
          do j = JE+1, JA
          do i = 1, IA
          do k = KS, KE-1
             MOMZ(k,i,j) = MOMZ(k,i,JE)
          end do
          end do
          end do
       end if
    end if

    return
  end subroutine set_boundary
  !-----------------------------------------------------------------------------
  !> Get boundaryal coefficient with same parent data
  subroutine get_boundary_same_parent( &
       bnd_DENS, &
       bnd_VELZ, &
       bnd_VELX, &
       bnd_VELY, &
       bnd_POTT, &
       bnd_QTRC, &
       now_step, &
       update_step )
    implicit none

    ! arguments
    real(RP), intent(out) :: bnd_DENS(:,:,:)
    real(RP), intent(out) :: bnd_VELZ(:,:,:)
    real(RP), intent(out) :: bnd_VELX(:,:,:)
    real(RP), intent(out) :: bnd_VELY(:,:,:)
    real(RP), intent(out) :: bnd_POTT(:,:,:)
    real(RP), intent(out) :: bnd_QTRC(:,:,:,:)
    integer,  intent(in)  :: now_step
    integer,  intent(in)  :: update_step

    ! works
    integer :: i, j, k, iq
    integer :: ref
    !---------------------------------------------------------------------------

    if ( now_step == update_step ) then
       ref = ref_new
    else
       ref = ref_now
    end if

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       bnd_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref)
       bnd_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref)
       bnd_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref)
       bnd_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref)
       bnd_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref)
       do iq = 1, BND_QA
          bnd_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref)
       end do
    end do
    end do
    end do

    return
  end subroutine get_boundary_same_parent

  !-----------------------------------------------------------------------------
  !> Get boundaryal coefficient with nearest neighbor
  subroutine get_boundary_nearest_neighbor( &
       bnd_DENS, &
       bnd_VELZ, &
       bnd_VELX, &
       bnd_VELY, &
       bnd_POTT, &
       bnd_QTRC, &
       now_step, &
       update_step )
    implicit none

    ! parameters
    real(RP) :: EPS = 1.0E-4_RP

    ! arguments
    real(RP), intent(out) :: bnd_DENS(:,:,:)
    real(RP), intent(out) :: bnd_VELZ(:,:,:)
    real(RP), intent(out) :: bnd_VELX(:,:,:)
    real(RP), intent(out) :: bnd_VELY(:,:,:)
    real(RP), intent(out) :: bnd_POTT(:,:,:)
    real(RP), intent(out) :: bnd_QTRC(:,:,:,:)
    integer,  intent(in)  :: now_step
    integer,  intent(in)  :: update_step

    ! works
    integer :: i, j, k, iq
    integer :: ref_idx

    real(RP) :: real_nstep
    real(RP) :: half_nstep
    !---------------------------------------------------------------------------

    real_nstep = real( now_step, kind=RP )
    half_nstep = real( UPDATE_NSTEP, kind=RP ) * 0.5_RP

    ! this step before half of the parent step
    if( ( real_nstep - EPS ) < half_nstep ) then
      ref_idx = ref_now

    ! this step after half of the parent step
    else if( ( real_nstep - 1.0_RP + EPS ) > half_nstep ) then
      ref_idx = ref_new

    ! this step across half of the parent step
    else
      ref_idx = ref_now

    end if

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       bnd_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_idx)
       bnd_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_idx)
       bnd_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_idx)
       bnd_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_idx)
       bnd_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_idx)
       do iq = 1, BND_QA
          bnd_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_idx)
       end do
    end do
    end do
    end do

    return
  end subroutine get_boundary_nearest_neighbor

  !-----------------------------------------------------------------------------
  !> Get boundaryal coefficient with linear interpotation between initial points
  subroutine get_boundary_lerp_initpoint( &
       bnd_DENS, &
       bnd_VELZ, &
       bnd_VELX, &
       bnd_VELY, &
       bnd_POTT, &
       bnd_QTRC, &
       now_step, &
       update_step )
    implicit none

    ! arguments
    real(RP), intent(out) :: bnd_DENS(:,:,:)
    real(RP), intent(out) :: bnd_VELZ(:,:,:)
    real(RP), intent(out) :: bnd_VELX(:,:,:)
    real(RP), intent(out) :: bnd_VELY(:,:,:)
    real(RP), intent(out) :: bnd_POTT(:,:,:)
    real(RP), intent(out) :: bnd_QTRC(:,:,:,:)
    integer,  intent(in)  :: now_step
    integer,  intent(in)  :: update_step

    ! works
    integer :: i, j, k, iq

    real(RP) :: fact
    !---------------------------------------------------------------------------

    fact = ( now_step + 0.5_RP ) / update_step

    !$omp parallel do default(none) private(i,j,k,iq) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KS,KE,bnd_DENS,ATMOS_BOUNDARY_ref_DENS,ref_now,fact,ref_new,bnd_VELZ) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELZ,bnd_VELX,ATMOS_BOUNDARY_ref_VELX,bnd_VELY) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELY,bnd_POTT,ATMOS_BOUNDARY_ref_POTT,BND_QA,bnd_QTRC,ATMOS_BOUNDARY_ref_QTRC)
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       bnd_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now) * ( 1.0_RP-fact ) &
                       + ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_new) * fact
       bnd_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_now) * ( 1.0_RP-fact ) &
                       + ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_new) * fact
       bnd_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now) * ( 1.0_RP-fact) &
                       + ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_new) * fact
       bnd_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now) * ( 1.0_RP-fact) &
                       + ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_new) * fact
       bnd_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now) * ( 1.0_RP-fact ) &
                       + ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_new) * fact
       do iq = 1, BND_QA
          bnd_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now) * ( 1.0_RP-fact ) &
                             + ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_new) * fact
       end do
    end do
    end do
    end do

    return
  end subroutine get_boundary_lerp_initpoint

  !-----------------------------------------------------------------------------
  !> Get boundaryal coefficient with linear interpotation between mid-points
  subroutine get_boundary_lerp_midpoint( &
       bnd_DENS, &
       bnd_VELZ, &
       bnd_VELX, &
       bnd_VELY, &
       bnd_POTT, &
       bnd_QTRC, &
       now_step, &
       update_step )
    use scale_time, only: &
       TIME_DTSEC
    implicit none

    ! parameters
    real(RP) :: EPS = 1.0E-4_RP

    ! arguments
    real(RP), intent(out) :: bnd_DENS(:,:,:)
    real(RP), intent(out) :: bnd_VELZ(:,:,:)
    real(RP), intent(out) :: bnd_VELX(:,:,:)
    real(RP), intent(out) :: bnd_VELY(:,:,:)
    real(RP), intent(out) :: bnd_POTT(:,:,:)
    real(RP), intent(out) :: bnd_QTRC(:,:,:,:)
    integer,  intent(in)  :: now_step
    integer,  intent(in)  :: update_step

    ! works
    integer :: i, j, k, iq

    real(RP) :: real_nstep
    real(RP) :: half_nstep
    real(RP) :: t1, t2
    !---------------------------------------------------------------------------

    real_nstep = real( now_step, kind=RP )
    half_nstep = real( UPDATE_NSTEP, kind=RP ) * 0.5_RP

    ! this step before half of the parent step
    if( ( real_nstep - EPS ) < half_nstep ) then

       t1 = TIME_DTSEC / ATMOS_BOUNDARY_UPDATE_DT * ( real_nstep + half_nstep - 0.5_RP )

       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          bnd_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now) * t1              &
                          - ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_old) * ( t1 - 1.0_RP )
          bnd_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_now) * t1              &
                          - ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_old) * ( t1 - 1.0_RP )
          bnd_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now) * t1              &
                          - ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_old) * ( t1 - 1.0_RP )
          bnd_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now) * t1              &
                          - ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_old) * ( t1 - 1.0_RP )
          bnd_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now) * t1              &
                          - ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_old) * ( t1 - 1.0_RP )
          do iq = 1, BND_QA
             bnd_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now) * t1              &
                                - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_old) * ( t1 - 1.0_RP )
          end do
       end do
       end do
       end do

    ! this step after half of the parent step
    else if( ( real_nstep - 1.0_RP + EPS ) > half_nstep ) then

       t1 = TIME_DTSEC / ATMOS_BOUNDARY_UPDATE_DT * ( real_nstep - half_nstep - 0.5_RP )

       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          bnd_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_new) * t1              &
                          - ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now) * ( t1 - 1.0_RP )
          bnd_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_new) * t1              &
                          - ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_now) * ( t1 - 1.0_RP )
          bnd_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_new) * t1              &
                          - ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now) * ( t1 - 1.0_RP )
          bnd_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_new) * t1              &
                          - ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now) * ( t1 - 1.0_RP )
          bnd_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_new) * t1              &
                          - ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now) * ( t1 - 1.0_RP )
          do iq = 1, BND_QA
             bnd_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_new) * t1              &
                                - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now) * ( t1 - 1.0_RP )
          end do
       end do
       end do
       end do

    ! this step across half of the parent step
    else

       t1 = TIME_DTSEC / ATMOS_BOUNDARY_UPDATE_DT * ( real_nstep + half_nstep - 1.0_RP )
       t2 = TIME_DTSEC / ATMOS_BOUNDARY_UPDATE_DT * ( real_nstep - half_nstep )

       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          bnd_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_new) * t2 * 0.25_RP                   &
                          + ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now) * ( t1 - t2 + 3.0_RP ) * 0.25_RP &
                          - ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_old) * ( t1 - 1.0_RP ) * 0.25_RP
          bnd_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_new) * t2 * 0.25_RP                   &
                          + ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_now) * ( t1 - t2 + 3.0_RP ) * 0.25_RP &
                          - ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_old) * ( t1 - 1.0_RP ) * 0.25_RP
          bnd_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_new) * t2 * 0.25_RP                   &
                          + ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now) * ( t1 - t2 + 3.0_RP ) * 0.25_RP &
                          - ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_old) * ( t1 - 1.0_RP ) * 0.25_RP
          bnd_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_new) * t2 * 0.25_RP                   &
                          + ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now) * ( t1 - t2 + 3.0_RP ) * 0.25_RP &
                          - ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_old) * ( t1 - 1.0_RP ) * 0.25_RP
          bnd_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_new) * t2 * 0.25_RP                   &
                          + ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now) * ( t1 - t2 + 3.0_RP ) * 0.25_RP &
                          - ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_old) * ( t1 - 1.0_RP ) * 0.25_RP
          do iq = 1, BND_QA
             bnd_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_new) * t2 * 0.25_RP                   &
                                + ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now) * ( t1 - t2 + 3.0_RP ) * 0.25_RP &
                                - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_old) * ( t1 - 1.0_RP )
          end do
       end do
       end do
       end do

    end if

    return
  end subroutine get_boundary_lerp_midpoint

  !-----------------------------------------------------------------------------
  !> Update indices of array of boundary references
  subroutine update_ref_index
    implicit none

    ! works
    integer :: ref_tmp
    !---------------------------------------------------------------------------

    ref_tmp = ref_old
    ref_old = ref_now
    ref_now = ref_new
    ref_new = ref_tmp

    return
  end subroutine update_ref_index

  subroutine history_bnd( &
       ATMOS_BOUNDARY_DENS, &
       ATMOS_BOUNDARY_VELZ, &
       ATMOS_BOUNDARY_VELX, &
       ATMOS_BOUNDARY_VELY, &
       ATMOS_BOUNDARY_POTT, &
       ATMOS_BOUNDARY_QTRC )
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_atmos_phy_mp_vars, only: &
       QS_MP
    implicit none
    real(RP), intent(in) :: ATMOS_BOUNDARY_DENS(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_VELZ(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_VELX(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_VELY(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_POTT(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_QTRC(KA,IA,JA,BND_QA)

    integer :: iq, iqb

    call FILE_HISTORY_in( ATMOS_BOUNDARY_DENS(:,:,:), 'DENS_BND', 'Boundary Density',               'kg/m3'             )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELZ(:,:,:), 'VELZ_BND', 'Boundary velocity z-direction',  'm/s',  dim_type='ZHXY' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELX(:,:,:), 'VELX_BND', 'Boundary velocity x-direction',  'm/s',  dim_type='ZXHY' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELY(:,:,:), 'VELY_BND', 'Boundary velocity y-direction',  'm/s',  dim_type='ZXYH' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_POTT(:,:,:), 'PT_BND',   'Boundary potential temperature', 'K'                 )
    do iq = 1, QA
       iqb = BND_IQ(iq)
       if ( iqb > 0 ) then
          call FILE_HISTORY_in( ATMOS_BOUNDARY_QTRC(:,:,:,iqb), trim(TRACER_NAME(iq))//'_BND', &
                                trim(TRACER_NAME(iq))//' in boundary', 'kg/kg' )
       end if
    enddo

    return
  end subroutine history_bnd

  subroutine calc_mass( ref )
    use scale_prc_cartesC, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_S, &
       PRC_HAS_N
    use scale_statistics, only: &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       VOL          => ATMOS_GRID_CARTESC_REAL_VOL,          &
       TOTVOL       => ATMOS_GRID_CARTESC_REAL_TOTVOL,       &
       AREAZXV_Y    => ATMOS_GRID_CARTESC_REAL_AREAZXV_Y,    &
       TOTAREAZUY_X => ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X, &
       TOTAREAZXV_Y => ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y
    use scale_atmos_refstate, only: &
       DENS_ref => ATMOS_REFSTATE_dens
    use mod_atmos_vars, only: &
       DENS
    implicit none
    integer,  intent(in)  :: ref

    real(DP) :: masstot, masstot_current
    real(DP) :: massflx
    real(DP) :: offset_band, offset_bias
    real(DP) :: ref_tot
    real(DP) :: flx_w, flx_e, flx_s, flx_n
    real(DP) :: ref_w, ref_e, ref_s, ref_n

    real(RP), target :: work_x(KA,JA), work_y(KA,IA)
    real(RP), pointer :: ptr(:,:)

    integer :: k, i, j

    ! total mass
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                           ATMOS_BOUNDARY_ref_DENS(:,:,:,ref),      & ! (in)
                           "DENS_bnd",                              & ! (in)
                           VOL(:,:,:), TOTVOL,                      & ! (in)
                           log_suppress = .true., global = .true.,  & ! (in)
                           sum = masstot                            ) ! (out)

!!$    call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
!!$                           DENS(:,:,:), "DENS_bnd_update",         & ! (in)
!!$                           VOL(:,:,:), TOTVOL,                     & ! (in)
!!$                           log_suppress = .true., global = .true., & ! (in)
!!$                           sum = masstot_current                   ) ! (out)


    ! West
    if ( .NOT. PRC_HAS_W ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          work_x(k,j) = ATMOS_BOUNDARY_ref_VELX(k,IS-1,j,ref) &
                      * ( ATMOS_BOUNDARY_ref_DENS(k,IS-1,j,ref) + ATMOS_BOUNDARY_ref_DENS(k,IS,j,ref) ) * 0.5_RP
       end do
       end do
       ptr => work_x
    else
       ptr => zero_x
    end if
    call STATISTICS_total( KA, KS, KE, JA, JS, JE, &
                           ptr(:,:), "MFLUX_bnd_w",                & ! (in)
                           AREAZUY_W(:,:), TOTAREAZUY_X(IS-1),     & ! (in)
                           log_suppress = .true., global = .true., & ! (in)
                           sum = flx_w                             ) ! (out)
    if ( .NOT. PRC_HAS_W ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          work_x(k,j) = DENS_ref(k,IS,j)
       end do
       end do
    end if
    call STATISTICS_total( KA, KS, KE, JA, JS, JE, &
                           ptr(:,:), "DENS_ref_w",                 & ! (in)
                           AREAZUY_W(:,:), TOTAREAZUY_X(IS-1),     & ! (in)
                           log_suppress = .true., global = .true., & ! (in)
                           sum = ref_w                             ) ! (out)

    ! East
    if ( .NOT. PRC_HAS_E ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          work_x(k,j) = ATMOS_BOUNDARY_ref_VELX(k,IE,j,ref) &
                      * ( ATMOS_BOUNDARY_ref_DENS(k,IE,j,ref) + ATMOS_BOUNDARY_ref_DENS(k,IE+1,j,ref) ) * 0.5_RP
       end do
       end do
       ptr => work_x
    else
       ptr => zero_x
    end if
    call STATISTICS_total( KA, KS, KE, JA, JS, JE, &
                           ptr(:,:), "MFLUX_bnd_e",                & ! (in)
                           AREAZUY_E(:,:), TOTAREAZUY_X(IE),       & ! (in)
                           log_suppress = .true., global = .true., & ! (in)
                           sum = flx_e                             ) ! (out)
    if ( .NOT. PRC_HAS_E ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          work_x(k,j) = DENS_ref(k,IE,j)
       end do
       end do
    end if
    call STATISTICS_total( KA, KS, KE, JA, JS, JE, &
                           ptr(:,:), "DENS_ref_e",                 & ! (in)
                           AREAZUY_E(:,:), TOTAREAZUY_X(IE),       & ! (in)
                           log_suppress = .true., global = .true., & ! (in)
                           sum = ref_e                             ) ! (out)

    ! South
    if ( .NOT. PRC_HAS_S ) then
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          work_y(k,i) = ATMOS_BOUNDARY_ref_VELY(k,i,JS-1,ref) &
                      * ( ATMOS_BOUNDARY_ref_DENS(k,i,JS-1,ref) + ATMOS_BOUNDARY_ref_DENS(k,i,JS,ref) ) * 0.5_RP
       end do
       end do
       ptr => work_y
    else
       ptr => zero_y
    end if
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, &
                           ptr(:,:), "MFLUX_bnd_s",                 & ! (in)
                           AREAZXV_Y(:,:,JS-1), TOTAREAZXV_Y(JS-1), & ! (in)
                           log_suppress = .true., global = .true.,  & ! (in)
                           sum = flx_s                              ) ! (out)
    if ( .NOT. PRC_HAS_S ) then
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          work_y(k,i) = DENS_ref(k,i,JS)
       end do
       end do
    end if
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, &
                           ptr(:,:), "DENS_ref_s",                  & ! (in)
                           AREAZXV_Y(:,:,JS-1), TOTAREAZXV_Y(JS-1), & ! (in)
                           log_suppress = .true., global = .true.,  & ! (in)
                           sum = ref_s                              ) ! (out)

    ! North
    if ( .NOT. PRC_HAS_N ) then
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          work_y(k,i) = ATMOS_BOUNDARY_ref_VELY(k,i,JE,ref) &
                      * ( ATMOS_BOUNDARY_ref_DENS(k,i,JE,ref) + ATMOS_BOUNDARY_ref_DENS(k,i,JE+1,ref) ) * 0.5_RP
       end do
       end do
       ptr => work_y
    else
       ptr => zero_y
    end if
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, &
                           ptr(:,:), "MFLX_bnd_n",                 & ! (in)
                           AREAZXV_Y(:,:,JE), TOTAREAZXV_Y(JE),    & ! (in)
                           log_suppress = .true., global = .true., & ! (in)
                           sum = flx_n                             ) ! (out)
    if ( .NOT. PRC_HAS_N ) then
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          work_y(k,i) = DENS_ref(k,i,JE)
       end do
       end do
    end if
    call STATISTICS_total( KA, KS, KE, IA, IS, IE, &
                           ptr(:,:), "DENS_ref_n",                 & ! (in)
                           AREAZXV_Y(:,:,JE), TOTAREAZXV_Y(JE),    & ! (in)
                           log_suppress = .true., global = .true., & ! (in)
                           sum = ref_n                             ) ! (out)

    massflx = flx_w - flx_e + flx_s - flx_n

    offset_band = ( masstot - MASSTOT_now ) / ATMOS_BOUNDARY_UPDATE_DT &
                - ( massflx + MASSFLX_now ) * 0.5_DP
!    offset_bias = ( MASSTOT_now - masstot_current ) / ATMOS_BOUNDARY_UPDATE_DT

    LOG_INFO("ATMOS_BOUNDARY_calc_mass",*) "Offset_band is: ", offset_band, "(", masstot, masstot_now, massflx, massflx_now, ")"

    ref_tot = ref_w + ref_e + ref_s + ref_n
    offset_band = offset_band / ref_tot
!    offset_bias = offset_bias / ref_tot

    LOG_INFO_CONT(*) "          per dens  ", offset_band

    ! density of the reference state is used as weight
    if ( .not. PRC_HAS_W ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          MFLUX_OFFSET_X(k,j,1,1) = offset_band * DENS_ref(k,IS,j)
!          MFLUX_OFFSET_X(k,j,1,2) = offset_bias * DENS_ref(k,IS,j)
       end do
       end do
    end if
    if ( .not. PRC_HAS_E ) then
       !$omp parallel do
       do j = JS, JE
       do k = KS, KE
          MFLUX_OFFSET_X(k,j,2,1) = - offset_band * DENS_ref(k,IE,j)
!          MFLUX_OFFSET_X(k,j,2,2) = - offset_bias * DENS_ref(k,IE,j)
       end do
       end do
    end if
    if ( .not. PRC_HAS_S ) then
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          MFLUX_OFFSET_Y(k,i,1,1) = offset_band * DENS_ref(k,i,JS)
!          MFLUX_OFFSET_Y(k,i,1,2) = offset_bias * DENS_ref(k,i,JS)
       end do
       end do
    end if
    if ( .not. PRC_HAS_N ) then
       !$omp parallel do
       do i = IS, IE
       do k = KS, KE
          MFLUX_OFFSET_Y(k,i,2,1) = - offset_band * DENS_ref(k,i,JE)
!          MFLUX_OFFSET_Y(k,i,2,2) = - offset_bias * DENS_ref(k,i,JE)
       end do
       end do
    end if

    MASSTOT_now = masstot
    MASSFLX_now = massflx

    return
  end subroutine calc_mass


  subroutine set_offset
    integer :: k, i, j, n

    !$omp parallel do
    do j = JS, JE
    do k = KS, KE
    do n = 1, 2
       ATMOS_BOUNDARY_MFLUX_OFFSET_X(k,j,n) = MFLUX_OFFSET_X(k,j,n,1) * OFFSET_TIME_FACT(now_step) !&
!                                            + MFLUX_OFFSET_X(k,j,n,2)
    end do
    end do
    end do

    !$omp parallel do
    do i = IS, IE
    do k = KS, KE
    do n = 1, 2
       ATMOS_BOUNDARY_MFLUX_OFFSET_Y(k,i,n) = MFLUX_OFFSET_Y(k,i,n,1) * OFFSET_TIME_FACT(now_step) !&
!                                            + MFLUX_OFFSET_Y(k,i,n,2)
    end do
    end do
    end do

    return
  end subroutine set_offset

end module mod_atmos_bnd_driver
