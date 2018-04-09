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
module scale_atmos_boundary
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
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
  public :: ATMOS_BOUNDARY_setup
  public :: ATMOS_BOUNDARY_set
  public :: ATMOS_BOUNDARY_firstsend
  public :: ATMOS_BOUNDARY_finalize
  public :: ATMOS_BOUNDARY_update

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public              :: BND_QA !> # of tracer at boundary

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


  real(RP), public              :: ATMOS_BOUNDARY_SMOOTHER_FACT  =  0.2_RP ! fact for smoother to damping

  logical,  public              :: ATMOS_BOUNDARY_UPDATE_FLAG = .false. !> switch for real case

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
  character(len=H_LONG),  private :: ATMOS_BOUNDARY_OUT_BASENAME = ''
  character(len=H_MID),   private :: ATMOS_BOUNDARY_OUT_TITLE    = 'SCALE-RM BOUNDARY CONDITION'  !< title of the output file
  character(len=H_SHORT), private :: ATMOS_BOUNDARY_OUT_DTYPE    = 'DEFAULT'                      !< REAL4 or REAL8

  logical,               private :: ATMOS_BOUNDARY_USE_DENS     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELZ     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELX     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_VELY     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_POTT     = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_QV       = .false. ! read from file?
  logical,               private :: ATMOS_BOUNDARY_USE_QHYD     = .false. ! read from file?

  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELZ   =   0.0_RP ! velocity w      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELX   =   0.0_RP ! velocity u      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_VELY   =   0.0_RP ! velocity v      at boundary, 0   [m/s]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_POTT   = 300.0_RP ! potential temp. at boundary, 300 [K]
  real(RP),              private :: ATMOS_BOUNDARY_VALUE_QTRC   =   0.0_RP ! tracer          at boundary, 0   [kg/kg]

  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_DENS = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_VELZ = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_VELX = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_VELY = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_POTT = 1.0_RP ! alpha factor again default
  real(RP),              private :: ATMOS_BOUNDARY_ALPHAFACT_QTRC = 1.0_RP ! alpha factor again default

  real(RP),              private :: ATMOS_BOUNDARY_FRACZ        =   1.0_RP ! fraction of boundary region for dumping (z) (0-1)
  real(RP),              private :: ATMOS_BOUNDARY_FRACX        =   1.0_RP ! fraction of boundary region for dumping (x) (0-1)
  real(RP),              private :: ATMOS_BOUNDARY_FRACY        =   1.0_RP ! fraction of boundary region for dumping (y) (0-1)
  real(RP),              private :: ATMOS_BOUNDARY_tauz                    ! maximum value for damping tau (z) [s]
  real(RP),              private :: ATMOS_BOUNDARY_taux                    ! maximum value for damping tau (x) [s]
  real(RP),              private :: ATMOS_BOUNDARY_tauy                    ! maximum value for damping tau (y) [s]

  real(DP), private              :: ATMOS_BOUNDARY_UPDATE_DT    =  0.0_DP ! inteval time of boudary data update [s]
  integer,  private              :: UPDATE_NSTEP

  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_DENS(:,:,:,:)   ! reference DENS (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELZ(:,:,:,:)   ! reference VELZ (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELX(:,:,:,:)   ! reference VELX (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELY(:,:,:,:)   ! reference VELY (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_POTT(:,:,:,:)   ! reference POTT (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,:) ! reference QTRC (with HALO)

  character(len=H_LONG), private :: ATMOS_BOUNDARY_interp_TYPE = 'lerp_initpoint' ! type of boundary interporation

  integer,               private :: ATMOS_BOUNDARY_START_DATE(6) = (/ -9999, 0, 0, 0, 0, 0 /) ! boundary initial date

  integer,               private :: ATMOS_BOUNDARY_fid = -1

  integer,               private :: now_step
  integer,               private :: boundary_timestep = 0
  logical,               private :: ATMOS_BOUNDARY_LINEAR_V = .false.  ! linear or non-linear profile of relax region
  logical,               private :: ATMOS_BOUNDARY_LINEAR_H = .true.   ! linear or non-linear profile of relax region
  real(RP),              private :: ATMOS_BOUNDARY_EXP_H    = 2.0_RP   ! factor of non-linear profile of relax region
  logical,               private :: ATMOS_BOUNDARY_ONLINE   = .false.  ! boundary online update by communicate inter-domain
  logical,               private :: ATMOS_BOUNDARY_ONLINE_MASTER = .false.  ! master domain in communicate inter-domain

  integer,               private :: QA_MP
  integer,               private :: QS_MP
  integer,               private :: QE_MP

  logical,               private :: do_parent_process       = .false.
  logical,               private :: do_daughter_process     = .false.
  logical,               private :: l_bnd = .false.

  real(DP),              private :: boundary_time_initdaysec

  integer,               private :: ref_size = 3
  integer,               private :: ref_old  = 1
  integer,               private :: ref_now  = 2
  integer,               private :: ref_new  = 3

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_BOUNDARY_setup( QA_MP_in, QS_MP_in, QE_MP_in )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CONST_UNDEF
    use scale_time, only: &
       DT => TIME_DTSEC
    use scale_comm_cartesC_nest, only: &
       USE_NESTING,         &
       OFFLINE,             &
       ONLINE_IAM_PARENT,   &
       ONLINE_IAM_DAUGHTER
    implicit none
    integer, intent(in) :: QA_MP_in, QS_MP_in, QE_MP_in

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_TYPE,           &
       ATMOS_BOUNDARY_IN_BASENAME,    &
       ATMOS_BOUNDARY_IN_CHECK_COORDINATES, &
       ATMOS_BOUNDARY_OUT_BASENAME,   &
       ATMOS_BOUNDARY_OUT_TITLE,      &
       ATMOS_BOUNDARY_OUT_DTYPE,      &
       ATMOS_BOUNDARY_USE_VELZ,       &
       ATMOS_BOUNDARY_USE_VELX,       &
       ATMOS_BOUNDARY_USE_VELY,       &
       ATMOS_BOUNDARY_USE_POTT,       &
       ATMOS_BOUNDARY_USE_DENS,       &
       ATMOS_BOUNDARY_USE_QV,         &
       ATMOS_BOUNDARY_USE_QHYD,       &
       ATMOS_BOUNDARY_VALUE_VELZ,     &
       ATMOS_BOUNDARY_VALUE_VELX,     &
       ATMOS_BOUNDARY_VALUE_VELY,     &
       ATMOS_BOUNDARY_VALUE_POTT,     &
       ATMOS_BOUNDARY_VALUE_QTRC,     &
       ATMOS_BOUNDARY_ALPHAFACT_DENS, &
       ATMOS_BOUNDARY_ALPHAFACT_VELZ, &
       ATMOS_BOUNDARY_ALPHAFACT_VELX, &
       ATMOS_BOUNDARY_ALPHAFACT_VELY, &
       ATMOS_BOUNDARY_ALPHAFACT_POTT, &
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
       ATMOS_BOUNDARY_interp_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[BOUNDARY] / Categ[ATMOS-RM] / Origin[SCALElib]'


    QA_MP = QA_MP_in
    QS_MP = QS_MP_in
    QE_MP = QE_MP_in


    ATMOS_BOUNDARY_tauz = DT * 10.0_RP
    ATMOS_BOUNDARY_taux = DT * 10.0_RP
    ATMOS_BOUNDARY_tauy = DT * 10.0_RP

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
       endif
    endif

    if( ATMOS_BOUNDARY_USE_QHYD ) then
       BND_QA = QA_MP
    else if ( QA_MP > 0 ) then
       BND_QA = 1
    else
       BND_QA = 0
    end if

    allocate( ATMOS_BOUNDARY_DENS(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_VELZ(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_VELX(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_VELY(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_POTT(KA,IA,JA) )
    allocate( ATMOS_BOUNDARY_QTRC(KA,IA,JA,BND_QA) )
    ATMOS_BOUNDARY_DENS(:,:,:)   = CONST_UNDEF
    ATMOS_BOUNDARY_VELZ(:,:,:)   = CONST_UNDEF
    ATMOS_BOUNDARY_VELX(:,:,:)   = CONST_UNDEF
    ATMOS_BOUNDARY_VELY(:,:,:)   = CONST_UNDEF
    ATMOS_BOUNDARY_POTT(:,:,:)   = CONST_UNDEF
    ATMOS_BOUNDARY_QTRC(:,:,:,:) = CONST_UNDEF

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
       ATMOS_BOUNDARY_ref_DENS(:,:,:,:)   = CONST_UNDEF
       ATMOS_BOUNDARY_ref_VELZ(:,:,:,:)   = CONST_UNDEF
       ATMOS_BOUNDARY_ref_VELX(:,:,:,:)   = CONST_UNDEF
       ATMOS_BOUNDARY_ref_VELY(:,:,:,:)   = CONST_UNDEF
       ATMOS_BOUNDARY_ref_POTT(:,:,:,:)   = CONST_UNDEF
       ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,:) = CONST_UNDEF

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

    elseif ( ATMOS_BOUNDARY_TYPE == 'NONE' ) then

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'CONST' ) then

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'INIT' ) then

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'FILE' ) then

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
    LOG_INFO_CONT(*) 'Is POTT used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_POTT
    LOG_INFO_CONT(*) 'Is DENS used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_DENS
    LOG_INFO_CONT(*) 'Is QV   used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_QV
    LOG_INFO_CONT(*) 'Is QHYD used in atmospheric boundary?          : ', ATMOS_BOUNDARY_USE_QHYD
    LOG_NEWLINE
    LOG_INFO_CONT(*) 'Atmospheric boundary VELZ values               : ', ATMOS_BOUNDARY_VALUE_VELZ
    LOG_INFO_CONT(*) 'Atmospheric boundary VELX values               : ', ATMOS_BOUNDARY_VALUE_VELX
    LOG_INFO_CONT(*) 'Atmospheric boundary VELY values               : ', ATMOS_BOUNDARY_VALUE_VELY
    LOG_INFO_CONT(*) 'Atmospheric boundary POTT values               : ', ATMOS_BOUNDARY_VALUE_POTT
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
       LOG_INFO_CONT(*) 'Lateral boundary interporation type                : ', ATMOS_BOUNDARY_interp_TYPE
    endif

    return
  end subroutine ATMOS_BOUNDARY_setup

  !-----------------------------------------------------------------------------
  !> set
  subroutine ATMOS_BOUNDARY_set( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)


    call ATMOS_BOUNDARY_firstsend( &
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )

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
  end subroutine ATMOS_BOUNDARY_set

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_var_fillhalo
    use scale_comm, only: &
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
    use scale_comm, only: &
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
    use scale_comm, only: &
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

    real(RP) :: coef_z, alpha_z1, alpha_z2
    real(RP) :: coef_x, alpha_x1, alpha_x2
    real(RP) :: coef_y, alpha_y1, alpha_y2
    real(RP) :: ee1, ee2

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

    !$omp parallel do default(none) &
    !$omp shared(JA,IA,KA,CBFZ,ATMOS_BOUNDARY_FRACZ,FBFZ,ATMOS_BOUNDARY_LINEAR_V,coef_z,CBFX)            &
    !$omp shared(ATMOS_BOUNDARY_FRACX,PI,FBFX,ATMOS_BOUNDARY_LINEAR_H,coef_x)     &
    !$omp shared(ATMOS_BOUNDARY_EXP_H,CBFY,ATMOS_BOUNDARY_FRACY,FBFY,coef_y,l_bnd)     &
    !$omp shared(do_daughter_process,ONLINE_USE_VELZ,ATMOS_BOUNDARY_USE_VELZ,ATMOS_BOUNDARY_alpha_VELZ)  &
    !$omp shared(ATMOS_BOUNDARY_ALPHAFACT_VELZ,ATMOS_BOUNDARY_USE_DENS,ATMOS_BOUNDARY_alpha_DENS)        &
    !$omp shared(ATMOS_BOUNDARY_ALPHAFACT_DENS,ATMOS_BOUNDARY_USE_VELX,ATMOS_BOUNDARY_alpha_VELX)        &
    !$omp shared(ATMOS_BOUNDARY_USE_VELY,ATMOS_BOUNDARY_alpha_VELY,ATMOS_BOUNDARY_ALPHAFACT_VELY)        &
    !$omp shared(ATMOS_BOUNDARY_USE_POTT,ATMOS_BOUNDARY_alpha_POTT,ATMOS_BOUNDARY_ALPHAFACT_POTT)        &
    !$omp shared(ATMOS_BOUNDARY_USE_QV,ATMOS_BOUNDARY_alpha_QTRC,ATMOS_BOUNDARY_ALPHAFACT_QTRC)          &
    !$omp shared(ATMOS_BOUNDARY_USE_QHYD,BND_QA,ATMOS_BOUNDARY_ALPHAFACT_VELX) &
    !$omp private(i,j,k,ee1,ee2,alpha_z1,alpha_z2,alpha_x1,alpha_x2,alpha_y1,alpha_y2,iq) OMP_SCHEDULE_ collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

       alpha_z1 = 0.0_RP
       alpha_z2 = 0.0_RP
       if ( ATMOS_BOUNDARY_LINEAR_V ) then
          alpha_z1 = coef_z * ee1
          alpha_z2 = coef_z * ee2
       else
          if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
             alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
             alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          endif
          if    ( ee2 > 0.0_RP .AND. ee2 <= 0.5_RP ) then
             alpha_z2 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee2*PI ) )
          elseif( ee2 > 0.5_RP .AND. ee2 <= 1.0_RP ) then
             alpha_z2 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee2-0.5_RP)*PI ) )
          endif
       endif

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

       if ( ATMOS_BOUNDARY_LINEAR_H ) then
          alpha_x1 = coef_x * ee1
          alpha_x2 = coef_x * ee2
       else
          alpha_x1 = coef_x * ee1 * exp( -(1.0_RP-ee1) * ATMOS_BOUNDARY_EXP_H )
          alpha_x2 = coef_x * ee2 * exp( -(1.0_RP-ee2) * ATMOS_BOUNDARY_EXP_H )
       end if

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

       if ( ATMOS_BOUNDARY_LINEAR_H ) then
          alpha_y1 = coef_y * ee1
          alpha_y2 = coef_y * ee2
       else
          alpha_y1 = coef_y * ee1 * exp( -(1.0_RP-ee1) * ATMOS_BOUNDARY_EXP_H )
          alpha_y2 = coef_y * ee2 * exp( -(1.0_RP-ee2) * ATMOS_BOUNDARY_EXP_H )
       end if


       if ( l_bnd ) then
          if ( do_daughter_process ) then ! online
             if ( ONLINE_USE_VELZ ) then
                if( ATMOS_BOUNDARY_USE_VELZ ) then
                   ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_z2, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELZ
                else
                   ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELZ
                endif
             else
                if ( ATMOS_BOUNDARY_USE_VELZ ) then
                   ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = alpha_z2 * ATMOS_BOUNDARY_ALPHAFACT_VELZ
                end if
             end if
          else ! offline
             if ( ATMOS_BOUNDARY_USE_VELZ ) then
                ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_z2, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELZ
             endif
          end if
          if ( ATMOS_BOUNDARY_USE_DENS ) then
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_DENS
          else
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = 0.0_RP
!             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_DENS
          endif
          if ( ATMOS_BOUNDARY_USE_VELX ) then
             ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_z1, alpha_x2, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELX
          else
             ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_x2, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELX
          endif
          if ( ATMOS_BOUNDARY_USE_VELY ) then
             ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_z1, alpha_x1, alpha_y2 ) * ATMOS_BOUNDARY_ALPHAFACT_VELY
          else
             ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_x1, alpha_y2 ) * ATMOS_BOUNDARY_ALPHAFACT_VELY
          endif
          if ( ATMOS_BOUNDARY_USE_POTT ) then
             ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_POTT
          else
             ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_POTT
          endif
          if ( ATMOS_BOUNDARY_USE_QV   ) then
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,1) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
          else
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,1) = max( alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
          endif
          if ( ATMOS_BOUNDARY_USE_QHYD ) then
             do iq = 2, BND_QA
                ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
             end do
          else
             do iq = 2, BND_QA
                ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
             end do
          endif
       else
          ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_DENS
          ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_z2, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELZ
          ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_z1, alpha_x2, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_VELX
          ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_z1, alpha_x1, alpha_y2 ) * ATMOS_BOUNDARY_ALPHAFACT_VELY
          ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_POTT
          do iq = 1, BND_QA
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 ) * ATMOS_BOUNDARY_ALPHAFACT_QTRC
          end do
       end if
    enddo
    enddo
    enddo

    if ( l_bnd ) then
       if ( .NOT. ONLINE_USE_VELZ .AND. .NOT. ATMOS_BOUNDARY_USE_VELZ ) then
          ATMOS_BOUNDARY_alpha_VELZ(:,:,:) = 0.0_RP
       end if
    else
       if ( .NOT. ATMOS_BOUNDARY_USE_DENS ) then
          ATMOS_BOUNDARY_alpha_DENS(:,:,:) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_VELZ ) then
          ATMOS_BOUNDARY_alpha_VELZ(:,:,:) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_VELX ) then
          ATMOS_BOUNDARY_alpha_VELX(:,:,:) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_VELY ) then
          ATMOS_BOUNDARY_alpha_VELY(:,:,:) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_POTT ) then
          ATMOS_BOUNDARY_alpha_POTT(:,:,:) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_QV   ) then
          if ( BND_QA > 0 ) ATMOS_BOUNDARY_alpha_QTRC(:,:,:,1) = 0.0_RP
       end if
       if ( .NOT. ATMOS_BOUNDARY_USE_QHYD ) then
          do iq = 2, BND_QA
             ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq) = 0.0_RP
          end do
       end if
    end if


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

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       ATMOS_BOUNDARY_DENS(k,i,j) = DENS(k,i,j)
       ATMOS_BOUNDARY_VELZ(k,i,j) = MOMZ(k,i,j) / ( DENS(k,i,j)+DENS(k+1,i,  j  ) ) * 2.0_RP
       ATMOS_BOUNDARY_POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = QTRC(k,i,j,iq)
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
    integer :: iq
    !---------------------------------------------------------------------------

    call FILE_CARTESC_open( ATMOS_BOUNDARY_IN_BASENAME, fid )

    if ( ATMOS_BOUNDARY_IN_CHECK_COORDINATES ) then
       call FILE_CARTESC_check_coordinates( fid, atmos=.true. )
    end if

    if (      ATMOS_BOUNDARY_USE_DENS &
         .OR. ATMOS_BOUNDARY_USE_VELZ &
         .OR. ATMOS_BOUNDARY_USE_VELX &
         .OR. ATMOS_BOUNDARY_USE_VELY &
         .OR. ATMOS_BOUNDARY_USE_POTT &
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

    if ( ATMOS_BOUNDARY_USE_POTT ) then
       call FILE_CARTESC_read( fid, 'POTT', 'ZXY', ATMOS_BOUNDARY_POTT(:,:,:) )
       call FILE_CARTESC_read( fid, 'ALPHA_POTT', 'ZXY', ATMOS_BOUNDARY_alpha_POTT(:,:,:) )
    endif

    if ( ATMOS_BOUNDARY_USE_QV   ) then
       call FILE_CARTESC_read( fid, 'QV', 'ZXY', ATMOS_BOUNDARY_QTRC(:,:,:,1) )
       call FILE_CARTESC_read( fid, 'ALPHA_QV', 'ZXY', ATMOS_BOUNDARY_alpha_QTRC(:,:,:,1) )
    endif

    if ( ATMOS_BOUNDARY_USE_QHYD ) then
       do iq = 2, BND_QA
          call FILE_CARTESC_read( fid, TRACER_NAME(iq), 'ZXY', ATMOS_BOUNDARY_QTRC(:,:,:,iq) )
          call FILE_CARTESC_read( fid, 'ALPHA_'//trim(TRACER_NAME(iq)), 'ZXY', ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq) )
       end do
    endif

    call FILE_CARTESC_close( fid )


    call ATMOS_BOUNDARY_var_fillhalo
    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_read

  !-----------------------------------------------------------------------------
  !> Write boundary data
  subroutine ATMOS_BOUNDARY_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ
    implicit none

    integer :: iq
    !---------------------------------------------------------------------------

    if (      ATMOS_BOUNDARY_USE_DENS &
         .OR. ATMOS_BOUNDARY_USE_VELZ &
         .OR. ATMOS_BOUNDARY_USE_VELX &
         .OR. ATMOS_BOUNDARY_USE_VELY &
         .OR. ATMOS_BOUNDARY_USE_POTT &
         ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_DENS(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'DENS', 'Reference Density', 'kg/m3', 'ZXY',           &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    end if
    if ( ATMOS_BOUNDARY_USE_DENS .OR. l_bnd ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_DENS(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_DENS', 'Alpha for DENS', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELZ .OR. (l_bnd .AND. ONLINE_USE_VELZ) ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_VELZ(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELZ', 'Reference Velocity w', 'm/s', 'ZHXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_VELZ(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELZ', 'Alpha for VELZ', '1', 'ZHXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELX .OR. l_bnd ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_VELX(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELX', 'Reference Velocity u', 'm/s', 'ZXHY',         &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_VELX(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELX', 'Alpha for VELX', '1', 'ZXHY',           &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELY .OR. l_bnd ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_VELY(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELY', 'Reference Velocity y', 'm/s', 'ZXYH',         &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_VELY(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELY', 'Alpha for VELY', '1', 'ZXYH',           &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_POTT .OR. l_bnd ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_POTT(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'POTT', 'Reference POTT', 'K', 'ZXY',                  &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_POTT(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_POTT', 'Alpha for POTT', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_QV   .OR. l_bnd ) then
       call FILE_CARTESC_write( ATMOS_BOUNDARY_QTRC(:,:,:,1),                          &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'QV', 'Reference QV', 'kg/kg', 'ZXY',                  &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,1),                    &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_QV', 'Alpha for QV', '1', 'ZXY',                &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_QHYD ) then
       do iq = 2, BND_QA
          call FILE_CARTESC_write( ATMOS_BOUNDARY_QTRC(:,:,:,iq),                                    &
                             ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE,            &
                             TRACER_NAME(iq), 'Reference '//trim(TRACER_NAME(iq)), &
                             TRACER_UNIT(iq), 'ZXY', &
                             ATMOS_BOUNDARY_OUT_DTYPE                                          )
          call FILE_CARTESC_write( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq),                                      &
                             ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE,                    &
                             'ALPHA_'//trim(TRACER_NAME(iq)), 'Alpha for '//trim(TRACER_NAME(iq)), &
                             '1', 'ZXY', &
                             ATMOS_BOUNDARY_OUT_DTYPE                                                  )
       end do
    endif

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

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_REFSTATE_DENS(k,i,j)
       ATMOS_BOUNDARY_VELZ(k,i,j) = ATMOS_BOUNDARY_VALUE_VELZ
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_VALUE_VELX
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_VALUE_VELY
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_VALUE_POTT
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

    integer  :: i, j, k, iq
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

    ! read boundary data from input file
    call ATMOS_BOUNDARY_update_file( ref_now )

    boundary_timestep = boundary_timestep + 1
    call ATMOS_BOUNDARY_update_file( ref_new )

    ! copy now to old
    !$omp parallel do default(none) private(i,j,k,iq) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KA,ATMOS_BOUNDARY_ref_DENS,ref_old,ref_now,ATMOS_BOUNDARY_ref_VELX) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELY,ATMOS_BOUNDARY_ref_POTT,BND_QA,ATMOS_BOUNDARY_ref_QTRC)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

    ! set boundary data
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now)
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now)
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now)
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now)
       end do
    end do
    end do
    end do

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          ATMOS_BOUNDARY_VELZ(k,i,j) = ATMOS_BOUNDARY_VALUE_VELZ
       end do
       end do
       end do
    end if

    now_step = fillgaps_steps

    ! get time boundary
    call get_boundary( bnd_DENS(:,:,:),   & ! [OUT]
                       bnd_VELZ(:,:,:),   & ! [OUT]
                       bnd_VELX(:,:,:),   & ! [OUT]
                       bnd_VELY(:,:,:),   & ! [OUT]
                       bnd_POTT(:,:,:),   & ! [OUT]
                       bnd_QTRC(:,:,:,:), & ! [OUT]
                       now_step,          & ! [IN]
                       UPDATE_NSTEP       ) ! [IN]

    ! fill in gaps of the offset
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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
       LOG_ERROR_CONT(*) 'check consistency between'
       LOG_ERROR_CONT(*) '    ONLINE_BOUNDARY_USE_QHYD and ATMOS_BOUNDARY_USE_QHYD.'
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
    use scale_time, only: &
       TIME_DTSEC,        &
       TIME_NSTEP
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ,          &
       PARENT_NSTEP
    implicit none

    ! parameters
    integer, parameter  :: handle = 2

    ! works
    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    ! import data from parent domain
    boundary_timestep = 1
    LOG_INFO("ATMOS_BOUNDARY_set_online",*) 'BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep

    call ATMOS_BOUNDARY_update_online_daughter( ref_now )

    boundary_timestep = boundary_timestep + 1
    LOG_INFO("ATMOS_BOUNDARY_set_online",*) 'BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep

    call ATMOS_BOUNDARY_update_online_daughter( ref_new )

    ! copy now to old
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

    ! set boundary data
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,ref_now)
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,ref_now)
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,ref_now)
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,ref_now)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,ref_now)
       end do
    end do
    end do
    end do

    if ( ONLINE_USE_VELZ ) then
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          ATMOS_BOUNDARY_VELZ(k,i,j) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,ref_now)
       end do
       end do
       end do
    else if ( ATMOS_BOUNDARY_USE_VELZ ) then
       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
          ATMOS_BOUNDARY_VELZ(k,i,j) = ATMOS_BOUNDARY_VALUE_VELZ
       end do
       end do
       end do
    end if

    UPDATE_NSTEP = nint( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
    if ( UPDATE_NSTEP * PARENT_NSTEP(handle) /= TIME_NSTEP ) then
       LOG_ERROR("ATMOS_BOUNDARY_set_online",*) 'NSTEP is not multiple of PARENT_NSTEP'
       call PRC_abort
    end if

    now_step = 0 ! should be set as zero in initialize process

    return
  end subroutine ATMOS_BOUNDARY_set_online

  !-----------------------------------------------------------------------------
  !> First send boundary value
  subroutine ATMOS_BOUNDARY_firstsend( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    use scale_comm_cartesC_nest, only: &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    implicit none

    ! arguments
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA_MP)
    !---------------------------------------------------------------------------

    ! send data at the first time
    if ( do_parent_process ) then !online [parent]
       ! issue send
       call ATMOS_BOUNDARY_send( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    endif

    return
  end subroutine ATMOS_BOUNDARY_firstsend

  !-----------------------------------------------------------------------------
  !> Finalize boundary value
  subroutine ATMOS_BOUNDARY_finalize
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
  end subroutine ATMOS_BOUNDARY_finalize

  !-----------------------------------------------------------------------------
  !> Update boundary value with a constant time boundary
  subroutine ATMOS_BOUNDARY_update( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_HAS_W,   &
       PRC_HAS_E,   &
       PRC_HAS_S,   &
       PRC_HAS_N
    use scale_comm_cartesC_nest, only: &
       ONLINE_USE_VELZ,       &
       COMM_CARTESC_NEST_test
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP) :: bnd_DENS(KA,IA,JA)        ! damping coefficient for DENS (0-1)
    real(RP) :: bnd_VELZ(KA,IA,JA)        ! damping coefficient for VELZ (0-1)
    real(RP) :: bnd_VELX(KA,IA,JA)        ! damping coefficient for VELX (0-1)
    real(RP) :: bnd_VELY(KA,IA,JA)        ! damping coefficient for VELY (0-1)
    real(RP) :: bnd_POTT(KA,IA,JA)        ! damping coefficient for POTT (0-1)
    real(RP) :: bnd_QTRC(KA,IA,JA,BND_QA) ! damping coefficient for QTRC (0-1)

    integer :: handle
    integer :: i, j, k, iq, iqa
    !---------------------------------------------------------------------------

    if ( do_parent_process ) then !online [parent]
       ! should be called every time step
       call ATMOS_BOUNDARY_update_online_parent( DENS,MOMZ,MOMX,MOMY,RHOT,QTRC(:,:,:,QS_MP:QE_MP) )
    endif

    if ( l_bnd ) then
       ! update referce vars
       if ( now_step >= UPDATE_NSTEP ) then
          now_step          = 0
          boundary_timestep = boundary_timestep + 1

          call update_ref_index

          if ( do_daughter_process ) then !online [daughter]
             call ATMOS_BOUNDARY_update_online_daughter( ref_new )
          else
             call ATMOS_BOUNDARY_update_file( ref_new )
          end if
       end if

       ! step boundary
       now_step = now_step + 1

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
       do j  = 1, JA
       do i  = 1, IA
       do k  = 1, KA
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
       if ( ONLINE_USE_VELZ ) then
          do j  = 1, JA
          do i  = 1, IA
          do k  = 1, KA
             ATMOS_BOUNDARY_VELZ(k,i,j) = bnd_VELZ(k,i,j)
          end do
          end do
          end do
       end if

       ! fill HALO in western region
       if ( .NOT. PRC_HAS_W ) then
          !$omp parallel do default(none)                                               &
          !$omp shared(JS,IS,KA,DENS,ATMOS_BOUNDARY_DENS,MOMX,ATMOS_BOUNDARY_VELX,RHOT) &
          !$omp shared(ATMOS_BOUNDARY_POTT,ATMOS_BOUNDARY_QTRC,BND_QA,QS_MP,QTRC,QA,JA)    &
          !$omp private(i,j,k,iq,iqa) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = 1, IS-1
          do k = 1, KA
             DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
             MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
             RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
             do iq = 1, BND_QA
                iqa = iq + QS_MP - 1
                QTRC(k,i,j,iqa) = ATMOS_BOUNDARY_QTRC(k,i,j,iq)
             end do
             do iq = 1, QA
                if ( iq < QS_MP .or. iq >= BND_QA + QS_MP ) then
                   QTRC(k,i,j,iq) = QTRC(k,IS,j,iq) &
                               * ( 0.5_RP - sign(0.5_RP, ATMOS_BOUNDARY_VELX(k,IS-1,j)) )
                end if
             end do
          end do
          end do
          end do
          do j = 1, JA-1
          do i = 1, IS-1
          do k = 1, KA
             MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
          end do
          end do
          end do
          do i = 1, IS-1
          do k = 1, KA
             MOMY(k,i,JA) = ATMOS_BOUNDARY_VELY(k,i,JA) &
                  * ATMOS_BOUNDARY_DENS(k,i,JA)
          end do
          end do
          if ( ONLINE_USE_VELZ ) then
             do j = 1, JA
             do i = 1, IS-1
             do k = KS, KE-1
                MOMZ(k,i,j) = ATMOS_BOUNDARY_VELZ(k,i,j) &
                     * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k+1,i,j) ) * 0.5_RP
             end do
             end do
             end do
          else
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
          !$omp parallel do default(none)                                             &
          !$omp shared(JA,IE,IA,KA,DENS,ATMOS_BOUNDARY_DENS,ATMOS_BOUNDARY_VELX,RHOT) &
          !$omp shared(ATMOS_BOUNDARY_POTT,ATMOS_BOUNDARY_QTRC,BND_QA,QS_MP,QTRC,QA)  &
          !$omp private(i,j,k,iq,iqa) OMP_SCHEDULE_ collapse(2)
          do j = 1, JA
          do i = IE+1, IA
          do k = 1, KA
             DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
             RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
             do iq = 1, BND_QA
                iqa = iq + QS_MP - 1
                QTRC(k,i,j,iqa) = ATMOS_BOUNDARY_QTRC(k,i,j,iq)
             end do
             do iq = 1, QA
                if ( iq < QS_MP .or. iq >= BND_QA + QS_MP ) then
                   QTRC(k,i,j,iq) = QTRC(k,IE,j,iq) &
                               * ( 0.5_RP + sign(0.5_RP, ATMOS_BOUNDARY_VELX(k,IE,j)) )
                end if
             end do
          end do
          end do
          end do
          do j = 1, JA
          do i = IE, IA-1
          do k = 1, KA
             MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
          end do
          end do
          end do
          do j = 1, JA
          do k = 1, KA
             MOMX(k,IA,j) = ATMOS_BOUNDARY_VELX(k,IA,j) * ATMOS_BOUNDARY_DENS(k,IA,j)
          end do
          end do
          do j = 1, JA-1
          do i = IE+1, IA
          do k = 1, KA
             MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
          end do
          end do
          end do
          do i = IE+1, IA
          do k = 1, KA
             MOMY(k,i,JA) = ATMOS_BOUNDARY_VELY(k,i,JA) &
                  * ATMOS_BOUNDARY_DENS(k,i,JA)
          end do
          end do
          if ( ONLINE_USE_VELZ ) then
             do j = 1, JA
             do i = IE+1, IA
             do k = KS, KE-1
                MOMZ(k,i,j) = ATMOS_BOUNDARY_VELZ(k,i,j) &
                     * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k+1,i,j) ) * 0.5_RP
             end do
             end do
             end do
          else
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
          !$omp parallel do default(none)                                               &
          !$omp shared(JS,IA,KA,DENS,ATMOS_BOUNDARY_DENS,MOMY,ATMOS_BOUNDARY_VELY,RHOT) &
          !$omp shared(ATMOS_BOUNDARY_POTT,ATMOS_BOUNDARY_QTRC,BND_QA,QS_MP,QTRC,QA)    &
          !$omp private(i,j,k,iq,iqa) OMP_SCHEDULE_ collapse(2)
          do j = 1, JS-1
          do i = 1, IA
          do k = 1, KA
             DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
             MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
             RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
             do iq = 1, BND_QA
                iqa = iq + QS_MP - 1
                QTRC(k,i,j,iqa) = ATMOS_BOUNDARY_QTRC(k,i,j,iq)
             end do
             do iq = 1, QA
                if ( iq < QS_MP .or. iq >= BND_QA + QS_MP ) then
                   QTRC(k,i,j,iq) = QTRC(k,i,JS,iq) &
                               * ( 0.5_RP - sign(0.5_RP, ATMOS_BOUNDARY_VELY(k,i,JS-1)) )
                end if
             end do
          end do
          end do
          end do
          do j = 1, JS-1
          do i = 1, IA-1
          do k = 1, KA
             MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
          end do
          end do
          end do
          do j = 1, JS-1
          do k = 1, KA
             MOMX(k,IA,j) = ATMOS_BOUNDARY_VELX(k,IA,j) &
                  * ATMOS_BOUNDARY_DENS(k,IA,j)
          end do
          end do
          if ( ONLINE_USE_VELZ ) then
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
          !$omp parallel do default(none)                                             &
          !$omp shared(JE,JA,IA,KA,DENS,ATMOS_BOUNDARY_DENS,ATMOS_BOUNDARY_VELY,RHOT) &
          !$omp shared(ATMOS_BOUNDARY_POTT,ATMOS_BOUNDARY_QTRC,BND_QA,QS_MP,QTRC,QA)  &
          !$omp private(i,j,k,iq,iqa) OMP_SCHEDULE_ collapse(2)
          do j = JE+1, JA
          do i = 1, IA
          do k = 1, KA
             DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j)
             RHOT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) * ATMOS_BOUNDARY_DENS(k,i,j)
             do iq = 1, BND_QA
                iqa = iq + QS_MP - 1
                QTRC(k,i,j,iqa) = ATMOS_BOUNDARY_QTRC(k,i,j,iq)
             end do
             do iq = BND_QA+1, QA
                if ( iq < QS_MP .or. iq >= BND_QA + QS_MP ) then
                   QTRC(k,i,j,iq) = QTRC(k,i,JE,iq) &
                               * ( 0.5_RP + sign(0.5_RP, ATMOS_BOUNDARY_VELY(k,i,JE)) )
                end if
             end do
          end do
          end do
          end do
          do j = JE, JA-1
          do i = 1, IA
          do k = 1, KA
             MOMY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i,j+1) ) * 0.5_RP
          end do
          end do
          end do
          do i = 1, IA
          do k = 1, KA
             MOMY(k,i,JA) = ATMOS_BOUNDARY_VELY(k,i,JA) * ATMOS_BOUNDARY_DENS(k,i,JA)
          end do
          end do
          do j = JE+1, JA
          do i = 1, IA-1
          do k = 1, KA
             MOMX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) &
                  * ( ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_DENS(k,i+1,j) ) * 0.5_RP
          end do
          end do
          end do
          do j = JE+1, JA
          do k = 1, KA
             MOMX(k,IA,j) = ATMOS_BOUNDARY_VELX(k,IA,j) &
                  * ATMOS_BOUNDARY_DENS(k,IA,j)
          end do
          end do
          if ( ONLINE_USE_VELZ ) then
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
  end subroutine ATMOS_BOUNDARY_update

  !-----------------------------------------------------------------------------
  !> Update reference boundary from file
  subroutine ATMOS_BOUNDARY_update_file( ref )
    use scale_file_cartesC, only: &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush
    implicit none

    integer, intent(in) :: ref

    integer :: fid, iq
    !---------------------------------------------------------------------------

    LOG_INFO("ATMOS_BOUNDARY_update_file",*) "Atmos Boundary: read from boundary file(timestep=", boundary_timestep, ")"

    fid = ATMOS_BOUNDARY_fid

    call FILE_CARTESC_read( fid, 'DENS', 'ZXY',  ATMOS_BOUNDARY_ref_DENS(:,:,:,ref), step=boundary_timestep )
    call FILE_CARTESC_read( fid, 'VELX', 'ZXHY', ATMOS_BOUNDARY_ref_VELX(:,:,:,ref), step=boundary_timestep )
    call FILE_CARTESC_read( fid, 'VELY', 'ZXYH', ATMOS_BOUNDARY_ref_VELY(:,:,:,ref), step=boundary_timestep )
    call FILE_CARTESC_read( fid, 'POTT', 'ZXY',  ATMOS_BOUNDARY_ref_POTT(:,:,:,ref), step=boundary_timestep )
    do iq = 1, BND_QA
       call FILE_CARTESC_read( fid, TRACER_NAME(QS_MP+iq-1), 'ZXY', ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,ref), step=boundary_timestep )
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
       QTRC )  ! [in]
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_recvwait_issue, &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    implicit none

    ! arguments
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA_MP)

    integer, parameter :: handle = 1
    !---------------------------------------------------------------------------

    LOG_INFO("ATMOS_BOUNDARY_update_online_parent",*)"ATMOS BOUNDARY update online: PARENT"

    ! issue wait
    call COMM_CARTESC_NEST_recvwait_issue( handle, NESTQA )

    ! issue send
    call ATMOS_BOUNDARY_send( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )

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
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_nestdown,    &
       DAUGHTER_KA,           &
       DAUGHTER_IA,           &
       DAUGHTER_JA,           &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    implicit none

    ! parameters
    integer, parameter  :: handle = 1

    ! arguments
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,NESTQA)

    ! works
    real(RP) :: dummy_d( DAUGHTER_KA(handle), DAUGHTER_IA(handle), DAUGHTER_JA(handle), NESTQA )
    !---------------------------------------------------------------------------

!OCL XFILL
    dummy_d(:,:,:,:) = 0.0_RP

    call COMM_CARTESC_NEST_nestdown( handle,           &
                                     NESTQA,           &
                                     DENS(:,:,:),      &  !(KA,IA,JA)
                                     MOMZ(:,:,:),      &  !(KA,IA,JA)
                                     MOMX(:,:,:),      &  !(KA,IA,JA)
                                     MOMY(:,:,:),      &  !(KA,IA,JA)
                                     RHOT(:,:,:),      &  !(KA,IA,JA)
                                     QTRC(:,:,:,:),    &  !(KA,IA,JA,QA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,1), &  !(KA,IA,JA)
                                     dummy_d(:,:,:,:)  )  !(KA,IA,JA,QA)

    return
  end subroutine ATMOS_BOUNDARY_send

  !-----------------------------------------------------------------------------
  !> Recieve boundary value
  subroutine ATMOS_BOUNDARY_recv( &
       ref_idx )
    use scale_comm_cartesC_nest, only: &
       ONLINE_BOUNDARY_DIAGQNUM, &
       COMM_CARTESC_NEST_nestdown,       &
       PARENT_KA,                &
       PARENT_IA,                &
       PARENT_JA,                &
       NESTQA => COMM_CARTESC_NEST_BND_QA
    use scale_prc, only: &
       PRC_abort
    implicit none

    ! parameters
    integer, parameter  :: handle = 2

    ! arguments
    integer, intent(in) :: ref_idx

    ! works
    real(RP) :: dummy_p( PARENT_KA(handle), PARENT_IA(handle), PARENT_JA(handle), NESTQA )
    !---------------------------------------------------------------------------

!OCL XFILL
    dummy_p(:,:,:,:) = 0.0_RP

    call COMM_CARTESC_NEST_nestdown( handle,                                         &
                                     NESTQA,                                         &
                                     dummy_p(:,:,:,1),                               & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                               & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                               & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                               & !(KA,IA,JA)
                                     dummy_p(:,:,:,1),                               & !(KA,IA,JA)
                                     dummy_p(:,:,:,1:NESTQA),                        & !(KA,IA,JA,QA)
                                     ATMOS_BOUNDARY_ref_DENS(:,:,:,ref_idx),         & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_VELZ(:,:,:,ref_idx),         & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_VELX(:,:,:,ref_idx),         & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_VELY(:,:,:,ref_idx),         & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_POTT(:,:,:,ref_idx),         & !(KA,IA,JA)
                                     ATMOS_BOUNDARY_ref_QTRC(:,:,:,1:NESTQA,ref_idx) ) !(KA,IA,JA,QA)

    if ( ONLINE_BOUNDARY_DIAGQNUM ) then
!       call ATMOS_HYDROMETEOR_diagnose_number_concentration( ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,ref_idx) ) ! [INOUT]
       LOG_ERROR("ATMOS_BOUNDARY_recv",*) 'tentative disabled'
       call PRC_abort
    endif

    return
  end subroutine ATMOS_BOUNDARY_recv

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

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

    fact = REAL(now_step, kind=RP) / update_step

    !$omp parallel do default(none) private(i,j,k,iq) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JA,IA,KA,bnd_DENS,ATMOS_BOUNDARY_ref_DENS,ref_now,fact,ref_new,bnd_VELZ) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELZ,bnd_VELX,ATMOS_BOUNDARY_ref_VELX,bnd_VELY) &
    !$omp shared(ATMOS_BOUNDARY_ref_VELY,bnd_POTT,ATMOS_BOUNDARY_ref_POTT,BND_QA,bnd_QTRC,ATMOS_BOUNDARY_ref_QTRC)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
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

       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
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

       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
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

       do j = 1, JA
       do i = 1, IA
       do k = 1, KA
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
    implicit none
    real(RP), intent(in) :: ATMOS_BOUNDARY_DENS(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_VELZ(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_VELX(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_VELY(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_POTT(KA,IA,JA)
    real(RP), intent(in) :: ATMOS_BOUNDARY_QTRC(KA,IA,JA,BND_QA)

    integer :: iq

    call FILE_HISTORY_in( ATMOS_BOUNDARY_DENS(:,:,:), 'DENS_BND', 'Boundary Density',               'kg/m3'             )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELZ(:,:,:), 'VELZ_BND', 'Boundary velocity z-direction',  'm/s',  dim_type='ZHXY' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELX(:,:,:), 'VELX_BND', 'Boundary velocity x-direction',  'm/s',  dim_type='ZXHY' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_VELY(:,:,:), 'VELY_BND', 'Boundary velocity y-direction',  'm/s',  dim_type='ZXYH' )
    call FILE_HISTORY_in( ATMOS_BOUNDARY_POTT(:,:,:), 'POTT_BND', 'Boundary potential temperature', 'K'                 )
    do iq = 1, BND_QA
       call FILE_HISTORY_in( ATMOS_BOUNDARY_QTRC(:,:,:,iq), trim(TRACER_NAME(iq))//'_BND', &
                     'Boundary '//trim(TRACER_NAME(iq)), 'kg/kg' )
    enddo

    return
  end subroutine history_bnd

end module scale_atmos_boundary
