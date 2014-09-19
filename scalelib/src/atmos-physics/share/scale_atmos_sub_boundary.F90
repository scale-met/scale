!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Boundary treatment
!!
!! @par Description
!!          Boundary treatment of model domain
!!          Additional forcing, Sponge layer, rayleigh dumping
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-12-07 (Y.Miyamoto) [new]
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE-LES ver.3
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2014-05-18 (R.Yoshida)  [add] boudary read/update for real case
!! @li      2014-09-05 (R.Yoshida)  [add] boudary update by online communicate
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_boundary
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer

  use gtool_file_h, only: &
     File_REAL4,  &
     File_REAL8
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_BOUNDARY_setup
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

  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_DENS(:,:,:)   !> damping coefficient for DENS [0-1]
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_VELZ(:,:,:)   !> damping coefficient for VELZ [0-1]
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_VELX(:,:,:)   !> damping coefficient for VELX [0-1]
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_VELY(:,:,:)   !> damping coefficient for VELY [0-1]
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_POTT(:,:,:)   !> damping coefficient for POTT [0-1]
  real(RP), public, allocatable :: ATMOS_BOUNDARY_alpha_QTRC(:,:,:,:) !> damping coefficient for QTRC [0-1]


  real(DP), private             :: ATMOS_BOUNDARY_UPDATE_DT      =  0.0_DP ! inteval time of boudary data update [s]
  real(RP), public              :: ATMOS_BOUNDARY_SMOOTHER_FACT  =  0.2_RP ! fact for smoother to damping

  logical,  public              :: ATMOS_BOUNDARY_UPDATE_FLAG = .false. !> switch for real case

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_BOUNDARY_initialize
  private :: ATMOS_BOUNDARY_update_file
  private :: ATMOS_BOUNDARY_update_online
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: ATMOS_BOUNDARY_TYPE         = 'NONE'
  character(len=H_LONG), private :: ATMOS_BOUNDARY_IN_BASENAME  = ''
  character(len=H_LONG), private :: ATMOS_BOUNDARY_OUT_BASENAME = ''
  character(len=H_MID),  private :: ATMOS_BOUNDARY_OUT_TITLE    = 'SCALE-LES BOUNDARY CONDITION' !< title of the output file
  character(len=H_MID),  private :: ATMOS_BOUNDARY_OUT_DTYPE    = 'DEFAULT'                      !< REAL4 or REAL8

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

  real(RP),              private :: ATMOS_BOUNDARY_FRACZ        =   1.0_RP ! fraction of boundary region for dumping (z) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_FRACX        =   1.0_RP ! fraction of boundary region for dumping (x) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_FRACY        =   1.0_RP ! fraction of boundary region for dumping (y) [0-1]
  real(RP),              private :: ATMOS_BOUNDARY_tauz                    ! maximum value for damping tau (z) [s]
  real(RP),              private :: ATMOS_BOUNDARY_taux                    ! maximum value for damping tau (x) [s]
  real(RP),              private :: ATMOS_BOUNDARY_tauy                    ! maximum value for damping tau (y) [s]

  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_DENS(:,:,:,:)   ! reference DENS (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELZ(:,:,:,:)   ! reference VELZ (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELX(:,:,:,:)   ! reference VELX (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_VELY(:,:,:,:)   ! reference VELY (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_POTT(:,:,:,:)   ! reference POTT (with HALO)
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_ref_QTRC(:,:,:,:,:) ! reference QTRC (with HALO)

  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment_DENS(:,:,:)   ! damping coefficient for DENS [0-1]
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment_VELZ(:,:,:)   ! damping coefficient for VELZ [0-1]
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment_VELX(:,:,:)   ! damping coefficient for VELX [0-1]
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment_VELY(:,:,:)   ! damping coefficient for VELY [0-1]
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment_POTT(:,:,:)   ! damping coefficient for POTT [0-1]
  real(RP),              private, allocatable :: ATMOS_BOUNDARY_increment_QTRC(:,:,:,:) ! damping coefficient for QTRC [0-1]

  integer,               private :: ATMOS_BOUNDARY_START_DATE(6) = (/ -9999, 0, 0, 0, 0, 0 /) ! boundary initial date

  real(DP),              private :: last_updated      =  -999.0_DP
  real(DP),              private :: integrated_sec    = 0.0_DP
  integer,               private :: boundary_timestep = 0
  logical,               private :: ATMOS_BOUNDARY_LINEARZ  = .false.  ! linear or non-linear profile of relax region
  logical,               private :: ATMOS_BOUNDARY_ONLINE   = .false.  ! boundary online update by communicate inter-domain
  logical,               private :: ATMOS_BOUNDARY_ONLINE_MASTER = .false.  ! master domain in communicate inter-domain
  logical,               private :: do_parent_process       = .false.
  logical,               private :: do_daughter_process     = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_BOUNDARY_setup( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF
    use scale_time, only: &
       DT => TIME_DTSEC
    use scale_grid_nest, only: &
       USE_NESTING,     &
       OFFLINE,         &
       ONLINE_IAM_PARENT,   &
       ONLINE_IAM_DAUGHTER
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMZ(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    NAMELIST / PARAM_ATMOS_BOUNDARY / &
       ATMOS_BOUNDARY_TYPE,           &
       ATMOS_BOUNDARY_IN_BASENAME,    &
       ATMOS_BOUNDARY_OUT_BASENAME,   &
       ATMOS_BOUNDARY_OUT_TITLE,      &
       ATMOS_BOUNDARY_USE_VELZ,       &
       ATMOS_BOUNDARY_USE_VELX,       &
       ATMOS_BOUNDARY_USE_VELY,       &
       ATMOS_BOUNDARY_USE_POTT,       &
       ATMOS_BOUNDARY_USE_DENS,       &
       ATMOS_BOUNDARY_USE_QV,         &
       ATMOS_BOUNDARY_USE_QHYD,       &
       ATMOS_BOUNDARY_VALUE_VELZ,     &
       ATMOS_BOUNDARY_VALUE_VELY,     &
       ATMOS_BOUNDARY_VALUE_VELX,     &
       ATMOS_BOUNDARY_VALUE_POTT,     &
       ATMOS_BOUNDARY_VALUE_QTRC,     &
       ATMOS_BOUNDARY_SMOOTHER_FACT,  &
       ATMOS_BOUNDARY_FRACZ,          &
       ATMOS_BOUNDARY_FRACX,          &
       ATMOS_BOUNDARY_FRACY,          &
       ATMOS_BOUNDARY_tauz,           &
       ATMOS_BOUNDARY_taux,           &
       ATMOS_BOUNDARY_tauy,           &
       ATMOS_BOUNDARY_UPDATE_DT,      &
       ATMOS_BOUNDARY_START_DATE,     &
       ATMOS_BOUNDARY_LINEARZ

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Boundary]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_BOUNDARY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_BOUNDARY. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML) write(IO_FID_LOG,nml=PARAM_ATMOS_BOUNDARY)

    if( .NOT. USE_NESTING ) then
       ATMOS_BOUNDARY_ONLINE = .false.
    else
       if( OFFLINE ) then
          ATMOS_BOUNDARY_ONLINE = .false.
       else
          ATMOS_BOUNDARY_ONLINE = .true.
       endif
    endif
    if( IO_L ) write(IO_FID_LOG,*) '*** Online Nesting for Lateral Boundary:', ATMOS_BOUNDARY_ONLINE

    if( ATMOS_BOUNDARY_USE_QHYD ) then
       BND_QA = QA
    else
       BND_QA = I_QV
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

    ATMOS_BOUNDARY_tauz = DT * 10.0_RP
    ATMOS_BOUNDARY_taux = DT * 10.0_RP
    ATMOS_BOUNDARY_tauy = DT * 10.0_RP

    if ( ATMOS_BOUNDARY_TYPE == 'CONST' ) then

       call ATMOS_BOUNDARY_generate

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'INIT' ) then

       call ATMOS_BOUNDARY_setinitval( DENS, & ! [IN]
                                       MOMZ, & ! [IN]
                                       MOMX, & ! [IN]
                                       MOMY, & ! [IN]
                                       RHOT, & ! [IN]
                                       QTRC  ) ! [IN]

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'FILE' ) then

       if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
          call ATMOS_BOUNDARY_read
       else
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_IN_BASENAME'
          call PRC_MPIstop
       endif

       ATMOS_BOUNDARY_UPDATE_FLAG = .false.

    elseif ( ATMOS_BOUNDARY_TYPE == 'REAL' ) then

       allocate( ATMOS_BOUNDARY_ref_DENS(KA,IA,JA,2) )
       allocate( ATMOS_BOUNDARY_ref_VELZ(KA,IA,JA,2) )
       allocate( ATMOS_BOUNDARY_ref_VELX(KA,IA,JA,2) )
       allocate( ATMOS_BOUNDARY_ref_VELY(KA,IA,JA,2) )
       allocate( ATMOS_BOUNDARY_ref_POTT(KA,IA,JA,2) )
       allocate( ATMOS_BOUNDARY_ref_QTRC(KA,IA,JA,BND_QA,2) )

       allocate( ATMOS_BOUNDARY_increment_DENS(KA,IA,JA) )
       allocate( ATMOS_BOUNDARY_increment_VELZ(KA,IA,JA) )
       allocate( ATMOS_BOUNDARY_increment_VELX(KA,IA,JA) )
       allocate( ATMOS_BOUNDARY_increment_VELY(KA,IA,JA) )
       allocate( ATMOS_BOUNDARY_increment_POTT(KA,IA,JA) )
       allocate( ATMOS_BOUNDARY_increment_QTRC(KA,IA,JA,BND_QA) )

       ! setting switches
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

       ! initialize boundary value (reading file or waiting parent domain)
       if ( .NOT. USE_NESTING ) then ! without nesting
          if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
             call ATMOS_BOUNDARY_initialize
          else
             write(*,*) 'xxx You need specify ATMOS_BOUNDARY_IN_BASENAME'
             call PRC_MPIstop
          endif
       elseif ( ATMOS_BOUNDARY_ONLINE_MASTER ) then ! with nesting: master domain
          if ( ATMOS_BOUNDARY_IN_BASENAME /= '' ) then
             call ATMOS_BOUNDARY_initialize
          else
             write(*,*) 'xxx You need specify ATMOS_BOUNDARY_IN_BASENAME'
             call PRC_MPIstop
          endif
       else ! with nesting: not master domain
          call ATMOS_BOUNDARY_initialize_online
       endif

       call ATMOS_BOUNDARY_setalpha

       ATMOS_BOUNDARY_UPDATE_FLAG = .true.

       if ( .NOT. ONLINE_IAM_DAUGHTER .and. ATMOS_BOUNDARY_UPDATE_DT <= 0.0_DP ) then
          write(*,*) 'xxx You need specify ATMOS_BOUNDARY_UPDATE_DT as larger than 0.0'
          call PRC_MPIstop
       endif

    else
       write(*,*) 'xxx unsupported ATMOS_BOUNDARY_TYPE. Check!', trim(ATMOS_BOUNDARY_TYPE)
       call PRC_MPIstop
    endif

    if( ATMOS_BOUNDARY_OUT_BASENAME /= '' ) then
       call ATMOS_BOUNDARY_write
    endif

    return
  end subroutine ATMOS_BOUNDARY_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine ATMOS_BOUNDARY_var_fillhalo
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    integer :: i, j, iq
    !---------------------------------------------------------------------------

    do j  = 1, JA
    do i  = 1, IA
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

    call COMM_wait ( ATMOS_BOUNDARY_DENS(:,:,:),   1 )
    call COMM_wait ( ATMOS_BOUNDARY_VELZ(:,:,:),   2 )
    call COMM_wait ( ATMOS_BOUNDARY_VELX(:,:,:),   3 )
    call COMM_wait ( ATMOS_BOUNDARY_VELY(:,:,:),   4 )
    call COMM_wait ( ATMOS_BOUNDARY_POTT(:,:,:),   5 )
    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_QTRC(:,:,:,iq), 5+iq )
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

    do j  = 1, JA
    do i  = 1, IA
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

    call COMM_wait ( ATMOS_BOUNDARY_alpha_DENS(:,:,:),   1 )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_VELZ(:,:,:),   2 )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_VELX(:,:,:),   3 )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_VELY(:,:,:),   4 )
    call COMM_wait ( ATMOS_BOUNDARY_alpha_POTT(:,:,:),   5 )
    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq), 5+iq )
    end do

    return
  end subroutine ATMOS_BOUNDARY_alpha_fillhalo

  !-----------------------------------------------------------------------------
  !> Calc dumping coefficient alpha
  subroutine ATMOS_BOUNDARY_setalpha
    use scale_const, only: &
       EPS => CONST_EPS, &
       PI  => CONST_PI
    use scale_grid, only: &
       CBFZ => GRID_CBFZ, &
       CBFX => GRID_CBFX, &
       CBFY => GRID_CBFY, &
       FBFZ => GRID_FBFZ, &
       FBFX => GRID_FBFX, &
       FBFY => GRID_FBFY

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
       if ( ATMOS_BOUNDARY_LINEARZ ) then
          alpha_z1 = coef_z * ee1
       else
          if    ( ee1 > 0.0_RP .AND. ee1 <= 0.5_RP ) then
             alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP - cos( ee1*PI ) )
          elseif( ee1 > 0.5_RP .AND. ee1 <= 1.0_RP ) then
             alpha_z1 = coef_z * 0.5_RP * ( 1.0_RP + sin( (ee1-0.5_RP)*PI ) )
          endif
       endif

       alpha_z2 = 0.0_RP
       if ( ATMOS_BOUNDARY_LINEARZ ) then
          alpha_z2 = coef_z * ee2
       else
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

       alpha_x1 = coef_x * ee1
       alpha_x2 = coef_x * ee2

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

       alpha_y1 = coef_y * ee1
       alpha_y2 = coef_y * ee2


       if ( ATMOS_BOUNDARY_TYPE == 'REAL' ) then
          ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = alpha_z2
          if ( ATMOS_BOUNDARY_USE_DENS ) then
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 )
          else                        
             ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_x1, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_VELX ) then
             ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_z1, alpha_x2, alpha_y1 )
          else                        
             ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_x2, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_VELY ) then
             ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_z1, alpha_x1, alpha_y2 )
          else                        
             ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_x1, alpha_y2 )
          endif
          if ( ATMOS_BOUNDARY_USE_POTT ) then
             ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 )
          else                        
             ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_x1, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_QV   ) then
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,1) = max( alpha_z1, alpha_x1, alpha_y1 )
          else
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,1) = max( alpha_x1, alpha_y1 )
          endif
          if ( ATMOS_BOUNDARY_USE_QHYD ) then
             do iq = 2, BND_QA
                ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 )
             end do
          else
             do iq = 2, BND_QA
                ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_x1, alpha_y1 )
             end do
          endif
       else
          ATMOS_BOUNDARY_alpha_DENS(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 )
          ATMOS_BOUNDARY_alpha_VELZ(k,i,j) = max( alpha_z2, alpha_x1, alpha_y1 )
          ATMOS_BOUNDARY_alpha_VELX(k,i,j) = max( alpha_z1, alpha_x2, alpha_y1 )
          ATMOS_BOUNDARY_alpha_VELY(k,i,j) = max( alpha_z1, alpha_x1, alpha_y2 )
          ATMOS_BOUNDARY_alpha_POTT(k,i,j) = max( alpha_z1, alpha_x1, alpha_y1 )
          do iq = 1, BND_QA
             ATMOS_BOUNDARY_alpha_QTRC(k,i,j,iq) = max( alpha_z1, alpha_x1, alpha_y1 )
          end do
       end if
    enddo
    enddo
    enddo

    if ( .NOT. ATMOS_BOUNDARY_USE_VELZ ) then
       ATMOS_BOUNDARY_alpha_VELZ(:,:,:) = 0.0_RP
    end if
    if ( .NOT. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       if ( .NOT. ATMOS_BOUNDARY_USE_DENS ) then
          ATMOS_BOUNDARY_alpha_DENS(:,:,:) = 0.0_RP
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
          ATMOS_BOUNDARY_alpha_QTRC(:,:,:,1) = 0.0_RP
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
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ATMOS_BOUNDARY_DENS(k,i,j) = DENS(k,i,j)
       ATMOS_BOUNDARY_VELZ(k,i,j) = MOMZ(k,i,j) / ( DENS(k,i,j)+DENS(k+1,i,  j  ) ) * 2.0_RP
       ATMOS_BOUNDARY_VELX(k,i,j) = MOMX(k,i,j) / ( DENS(k,i,j)+DENS(k,  i+1,j  ) ) * 2.0_RP
       ATMOS_BOUNDARY_VELY(k,i,j) = MOMY(k,i,j) / ( DENS(k,i,j)+DENS(k,  i,  j+1) ) * 2.0_RP
       ATMOS_BOUNDARY_POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = QTRC(k,i,j,iq)
       end do
    enddo
    enddo
    enddo

    call ATMOS_BOUNDARY_var_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_setinitval

  !-----------------------------------------------------------------------------
  !> Read boundary data
  subroutine ATMOS_BOUNDARY_read
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: bname

    integer :: iq
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME

    if (      ATMOS_BOUNDARY_USE_DENS &
         .or. ATMOS_BOUNDARY_USE_VELZ &
         .or. ATMOS_BOUNDARY_USE_VELX &
         .or. ATMOS_BOUNDARY_USE_VELY &
         .or. ATMOS_BOUNDARY_USE_POTT &
         ) then
       call FileRead( reference_atmos(:,:,:), bname, 'DENS', 1, PRC_myrank )
       ATMOS_BOUNDARY_DENS(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    end if
    if ( ATMOS_BOUNDARY_USE_DENS ) then
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_DENS', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha_DENS(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELZ', 1, PRC_myrank )
       ATMOS_BOUNDARY_VELZ(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_VELZ', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha_VELZ(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELX ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELX', 1, PRC_myrank )
       ATMOS_BOUNDARY_VELX(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_VELX', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha_VELX(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_VELY ) then
       call FileRead( reference_atmos(:,:,:), bname, 'VELY', 1, PRC_myrank )
       ATMOS_BOUNDARY_VELY(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_VELY', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha_VELY(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_POTT ) then
       call FileRead( reference_atmos(:,:,:), bname, 'POTT', 1, PRC_myrank )
       ATMOS_BOUNDARY_POTT(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_POTT', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha_POTT(KS:KE,IS:IE,JS:JE) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_QV   ) then
       call FileRead( reference_atmos(:,:,:), bname, 'QV',   1, PRC_myrank )
       ATMOS_BOUNDARY_QTRC(KS:KE,IS:IE,JS:JE,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_QV', 1, PRC_myrank )
       ATMOS_BOUNDARY_alpha_QTRC(KS:KE,IS:IE,JS:JE,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    endif

    if ( ATMOS_BOUNDARY_USE_QHYD ) then
       do iq = 2, BND_QA
          call FileRead( reference_atmos(:,:,:), bname, AQ_NAME(iq), 1, PRC_myrank )
          ATMOS_BOUNDARY_QTRC(KS:KE,IS:IE,JS:JE,iq) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
          call FileRead( reference_atmos(:,:,:), bname, 'ALPHA_'//trim(AQ_NAME(iq)), 1, PRC_myrank )
          ATMOS_BOUNDARY_alpha_QTRC(KS:KE,IS:IE,JS:JE,iq) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
       end do
    endif

    call ATMOS_BOUNDARY_var_fillhalo
    call ATMOS_BOUNDARY_alpha_fillhalo

    return
  end subroutine ATMOS_BOUNDARY_read

  !-----------------------------------------------------------------------------
  !> Write boundary data
  subroutine ATMOS_BOUNDARY_write
    use scale_fileio, only: &
       FILEIO_write
    implicit none

    real(RP) :: buffer(KA,IA,JA)

    integer :: iq
    !---------------------------------------------------------------------------

    if (      ATMOS_BOUNDARY_USE_DENS &
         .or. ATMOS_BOUNDARY_USE_VELZ &
         .or. ATMOS_BOUNDARY_USE_VELX &
         .or. ATMOS_BOUNDARY_USE_VELY &
         .or. ATMOS_BOUNDARY_USE_POTT &
         ) then
       call FILEIO_write( ATMOS_BOUNDARY_DENS(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'DENS', 'Reference Density', 'kg/m3', 'ZXY',           &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    end if
    if ( ATMOS_BOUNDARY_USE_DENS .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_alpha_DENS(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_DENS', 'Alpha for DENS', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       call FILEIO_write( ATMOS_BOUNDARY_VELZ(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELZ', 'Reference Velocity w', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha_VELZ(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELZ', 'Alpha for VELZ', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELX .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_VELX(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELX', 'Reference Velocity u', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha_VELX(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELX', 'Alpha for VELX', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_VELY .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_VELY(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'VELY', 'Reference Velocity y', 'm/s', 'ZXY',          &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha_VELY(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_VELY', 'Alpha for VELY', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_POTT .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_POTT(:,:,:),                            &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'POTT', 'Reference POTT', 'K', 'ZXY',                  &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha_POTT(:,:,:),                      &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_POTT', 'Alpha for POTT', '1', 'ZXY',            &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_QV   .or. ATMOS_BOUNDARY_TYPE == 'REAL' ) then
       call FILEIO_write( ATMOS_BOUNDARY_QTRC(:,:,:,1),                          &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'QV', 'Reference QV', 'kg/kg', 'ZXY',                  &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
       call FILEIO_write( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,1),                    &
                          ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE, &
                          'ALPHA_QV', 'Alpha for QV', '1', 'ZXY',                &
                          ATMOS_BOUNDARY_OUT_DTYPE                               )
    endif

    if ( ATMOS_BOUNDARY_USE_QHYD ) then
       do iq = 2, BND_QA
          call FILEIO_write( ATMOS_BOUNDARY_QTRC(:,:,:,iq),                                    &
                             ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE,            &
                             AQ_NAME(iq), 'Reference '//trim(AQ_NAME(iq)), AQ_UNIT(iq), 'ZXY', &
                             ATMOS_BOUNDARY_OUT_DTYPE                                          )
          call FILEIO_write( ATMOS_BOUNDARY_alpha_QTRC(:,:,:,iq),                                      &
                             ATMOS_BOUNDARY_OUT_BASENAME, ATMOS_BOUNDARY_OUT_TITLE,                    &
                             'ALPHA_'//trim(AQ_NAME(iq)), 'Alpha for '//trim(AQ_NAME(iq)), '1', 'ZXY', &
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
  subroutine ATMOS_BOUNDARY_initialize()
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_time, only: &
       TIME_NOWDATE,      &
       TIME_OFFSET_YEAR,  &
       TIME_DTSEC,        &
       TIME_NOWDAYSEC
    use scale_calendar, only: &
       CALENDAR_date2daysec,    &
       CALENDAR_combine_daysec
    implicit none
    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    integer  :: run_time_startdate(6)
    integer  :: run_time_startday
    real(RP) :: run_time_startsec
    real(RP) :: run_time_startms
    real(RP) :: run_time_nowdaysec
    integer  :: boundary_time_startday
    real(RP) :: boundary_time_startsec
    real(RP) :: boundary_time_startms
    real(RP) :: boundary_time_initdaysec
    real(RP) :: boundary_diff_daysec
    real(RP) :: boundary_inc_offset
    integer  :: fillgaps_steps

    character(len=H_LONG) :: bname

    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    bname = ATMOS_BOUNDARY_IN_BASENAME
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I5.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') '*** BOUNDARY START Date     : ', &
         ATMOS_BOUNDARY_START_DATE(1),'/',ATMOS_BOUNDARY_START_DATE(2),'/',ATMOS_BOUNDARY_START_DATE(3),' ',  &
         ATMOS_BOUNDARY_START_DATE(4),':',ATMOS_BOUNDARY_START_DATE(5),':',ATMOS_BOUNDARY_START_DATE(6)

    if ( ATMOS_BOUNDARY_START_DATE(1) == -9999 ) then
       boundary_timestep = 1
       boundary_inc_offset = 0.D0
       fillgaps_steps = 0
    else
       !--- recalculate time of the run [no offset]
       run_time_startms = 0.0_DP
       run_time_startdate(:) = TIME_NOWDATE(:)
       run_time_startdate(1) = TIME_OFFSET_YEAR
       call CALENDAR_date2daysec( run_time_startday,     & ! [OUT]
                                  run_time_startsec,     & ! [OUT]
                                  run_time_startdate(:), & ! [IN]
                                  run_time_startms       ) ! [IN]
       run_time_nowdaysec  = CALENDAR_combine_daysec( run_time_startday, run_time_startsec )

       !--- calculate time of the initial step in boundary file [no offset]
       boundary_time_startms = 0.0_DP
       call CALENDAR_date2daysec( boundary_time_startday,       & ! [OUT]
                                  boundary_time_startsec,       & ! [OUT]
                                  ATMOS_BOUNDARY_START_DATE(:), & ! [IN]
                                  boundary_time_startms         ) ! [IN]
       boundary_time_initdaysec  = CALENDAR_combine_daysec( boundary_time_startday, boundary_time_startsec )

       boundary_diff_daysec = run_time_nowdaysec - boundary_time_initdaysec
       boundary_timestep = 1 + int( boundary_diff_daysec / ATMOS_BOUNDARY_UPDATE_DT )
       boundary_inc_offset = mod( boundary_diff_daysec, ATMOS_BOUNDARY_UPDATE_DT )
       fillgaps_steps = int( boundary_inc_offset / TIME_DTSEC )
    endif

    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep
    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY INCREMENT OFFSET:', boundary_inc_offset
    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY FILLGAPS STEPS:', fillgaps_steps

    ! read boundary data from input file
    call FileRead( reference_atmos(:,:,:), bname, 'DENS', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_DENS(KS:KE,IS:IE,JS:JE,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_VELX(KS:KE,IS:IE,JS:JE,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_VELY(KS:KE,IS:IE,JS:JE,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_POTT(KS:KE,IS:IE,JS:JE,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    do iq = 1, BND_QA
       call FileRead( reference_atmos(:,:,:), bname, AQ_NAME(iq), boundary_timestep, PRC_myrank )
       ATMOS_BOUNDARY_ref_QTRC(KS:KE,IS:IE,JS:JE,iq,1) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    end do

    boundary_timestep = boundary_timestep + 1
    call FileRead( reference_atmos(:,:,:), bname, 'DENS', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_DENS(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_VELX(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_VELY(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_POTT(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    do iq = 1, BND_QA
       call FileRead( reference_atmos(:,:,:), bname, AQ_NAME(iq), boundary_timestep, PRC_myrank )
       ATMOS_BOUNDARY_ref_QTRC(KS:KE,IS:IE,JS:JE,iq,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    end do

    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,1)
       ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,1)
       ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,1)
       ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,1)
       ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,1)

       ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,1)
       ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,1)
       ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,1)
       ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,1)
       ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,1)

       ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,2)

       ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,2)

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,1) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,1)
          ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,1) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,1)

          ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,2)
          ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,2)
       end do
    end do
    end do

    call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,1),  1 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,1),  2 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,1),  3 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,1),  4 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,1),  5 )

    call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,2),  6 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2),  7 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,2),  8 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,2),  9 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 10 )

    do iq = 1, BND_QA
       call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,1), 10+iq        )
       call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 10+iq+BND_QA )
    end do

    call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,1),  1 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,1),  2 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,1),  3 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,1),  4 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,1),  5 )

    call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,2),  6 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2),  7 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,2),  8 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,2),  9 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 10 )

    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,1), 10+iq        )
       call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 10+iq+BND_QA )
    end do

    ! set boundary data and time increment
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,1)
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,1)
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,1)
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,1)

       ATMOS_BOUNDARY_increment_DENS(k,i,j) = ( ATMOS_BOUNDARY_ref_DENS(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_DENS(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       ATMOS_BOUNDARY_increment_VELX(k,i,j) = ( ATMOS_BOUNDARY_ref_VELX(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_VELX(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       ATMOS_BOUNDARY_increment_VELY(k,i,j) = ( ATMOS_BOUNDARY_ref_VELY(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_VELY(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       ATMOS_BOUNDARY_increment_POTT(k,i,j) = ( ATMOS_BOUNDARY_ref_POTT(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_POTT(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1)

          ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq) = ( ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,2) &
                                                    - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1) ) &
                                                  / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       end do
    end do
    end do
    end do

    ! fill in gaps of the offset
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_increment_DENS(k,i,j) * dble(fillgaps_steps)
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) + ATMOS_BOUNDARY_increment_VELX(k,i,j) * dble(fillgaps_steps)
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) + ATMOS_BOUNDARY_increment_VELY(k,i,j) * dble(fillgaps_steps)
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) + ATMOS_BOUNDARY_increment_POTT(k,i,j) * dble(fillgaps_steps)
       do iq = 1, BND_QA
         ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iq) &
                                       + ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq) * dble(fillgaps_steps)
       end do
    end do
    end do
    end do

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       ATMOS_BOUNDARY_VELZ(:,:,:) = ATMOS_BOUNDARY_VALUE_VELZ
    end if

    last_updated = TIME_NOWDAYSEC - boundary_inc_offset

    return
  end subroutine ATMOS_BOUNDARY_initialize

  !-----------------------------------------------------------------------------
  !> Initialize boundary value for real case experiment [online daughter]
  subroutine ATMOS_BOUNDARY_initialize_online()
    use scale_process, only: &
       PRC_myrank
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_time, only: &
       TIME_DTSEC,        &
       TIME_NOWDAYSEC
    use scale_grid_nest, only: &
       NEST_COMM_nestdown, &
       PARENT_KA,          &
       PARENT_IA,          &
       PARENT_JA,          &
       DAUGHTER_KA,        &
       DAUGHTER_IA,        &
       DAUGHTER_JA,        &
       PRNT_KS,            &
       PRNT_KE,            &
       PRNT_IS,            &
       PRNT_IE,            &
       PRNT_JS,            &
       PRNT_JE,            &
       DATR_KS,            &
       DATR_KE,            &
       DATR_IS,            &
       DATR_IE,            &
       DATR_JS,            &
       DATR_JE,            &
       PARENT_DTSEC
    implicit none

    integer, parameter  :: handle = 2

    real(RP) :: dummy1_p(PARENT_KA(handle),  PARENT_IA(handle),  PARENT_JA(handle))
    real(RP) :: dummy2_p(PARENT_KA(handle),  PARENT_IA(handle),  PARENT_JA(handle),  BND_QA)

    integer  :: i, j, k, iq
    !---------------------------------------------------------------------------

    ATMOS_BOUNDARY_UPDATE_DT = PARENT_DTSEC(handle)

    dummy1_p(:,:,:)   = 0.0_RP
    dummy2_p(:,:,:,:) = 0.0_RP

    ! import data from parent domain
    boundary_timestep = 1
    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep
    if( IO_L ) write(IO_FID_LOG,*) '+++ waiting for communication accept from parent domain'

    call NEST_COMM_nestdown( handle,                                 &
                             BND_QA,                                 &
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy2_p(:,:,:,1:BND_QA),               &   !(KA,IA,JA,QA)
                             ATMOS_BOUNDARY_ref_DENS(:,:,:,1),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_VELX(:,:,:,1),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_VELY(:,:,:,1),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_POTT(:,:,:,1),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_QTRC(:,:,:,1:BND_QA,1) ) !(KA,IA,JA,QA)

    boundary_timestep = boundary_timestep + 1
    if( IO_L ) write(IO_FID_LOG,*) '+++ BOUNDARY TIMESTEP NUMBER FOR INIT:', boundary_timestep
    if( IO_L ) write(IO_FID_LOG,*) '+++ waiting for communication accept from parent domain'

    call NEST_COMM_nestdown( handle,                                 &
                             BND_QA,                                 &
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                             dummy2_p(:,:,:,1:BND_QA),               &   !(KA,IA,JA,QA)
                             ATMOS_BOUNDARY_ref_DENS(:,:,:,2),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_VELX(:,:,:,2),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_VELY(:,:,:,2),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_POTT(:,:,:,2),       &   !(KA,IA,JA)
                             ATMOS_BOUNDARY_ref_QTRC(:,:,:,1:BND_QA,2) ) !(KA,IA,JA,QA)

    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,1)
       ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,1)
       ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,1)
       ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,1)
       ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,1) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,1)

       ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,1)
       ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,1)
       ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,1)
       ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,1)
       ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,1) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,1)

       ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,2)

       ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,2)

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,1) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,1)
          ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,1) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,1)

          ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,2)
          ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,2)
       end do
    end do
    end do

    call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,1),  1 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,1),  2 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,1),  3 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,1),  4 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,1),  5 )

    call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,2),  6 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2),  7 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,2),  8 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,2),  9 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 10 )

    do iq = 1, BND_QA
       call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,1), 10+iq        )
       call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 10+iq+BND_QA )
    end do

    call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,1),  1 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,1),  2 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,1),  3 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,1),  4 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,1),  5 )

    call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,2),  6 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2),  7 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,2),  8 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,2),  9 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 10 )

    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,1), 10+iq        )
       call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 10+iq+BND_QA )
    end do

    ! set boundary data and time increment
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,1)
       ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,1)
       ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,1)
       ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,1)

       ATMOS_BOUNDARY_increment_DENS(k,i,j) = ( ATMOS_BOUNDARY_ref_DENS(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_DENS(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       ATMOS_BOUNDARY_increment_VELX(k,i,j) = ( ATMOS_BOUNDARY_ref_VELX(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_VELX(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       ATMOS_BOUNDARY_increment_VELY(k,i,j) = ( ATMOS_BOUNDARY_ref_VELY(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_VELY(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       ATMOS_BOUNDARY_increment_POTT(k,i,j) = ( ATMOS_BOUNDARY_ref_POTT(k,i,j,2) &
                                              - ATMOS_BOUNDARY_ref_POTT(k,i,j,1) ) &
                                            / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1)

          ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq) = ( ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,2) &
                                                    - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1) ) &
                                                  / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
       end do
    end do
    end do
    end do

    if ( ATMOS_BOUNDARY_USE_VELZ ) then
       ATMOS_BOUNDARY_VELZ(:,:,:) = ATMOS_BOUNDARY_VALUE_VELZ
    end if

    last_updated = TIME_NOWDAYSEC + ATMOS_BOUNDARY_UPDATE_DT

    return
  end subroutine ATMOS_BOUNDARY_initialize_online

  !-----------------------------------------------------------------------------
  !> Update boundary value with a constant time increment
  subroutine ATMOS_BOUNDARY_update ( &
       DENS, MOMX, MOMY, RHOT, QTRC  )
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_time, only: &
       TIME_DTSEC, &
       TIME_NOWDAYSEC
    use scale_grid_nest, only: &
       PARENT_DTSEC
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)

    integer :: handle
    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if ( ATMOS_BOUNDARY_TYPE == 'REAL' ) then

       if ( ATMOS_BOUNDARY_ONLINE ) then  !online
          if ( do_daughter_process ) then !online [daughter]

             !ATMOS_BOUNDARY_UPDATE_DT = PARENT_DTSEC
             handle = 2
             integrated_sec = TIME_NOWDAYSEC - last_updated
   
             if ( integrated_sec >= ATMOS_BOUNDARY_UPDATE_DT - EPS ) then
                boundary_timestep = boundary_timestep + 1

                call ATMOS_BOUNDARY_update_online( DENS,MOMX,MOMY,RHOT,QTRC,handle )

                do j  = 1, JA
                do i  = 1, IA
                do k  = 1, KA
                   ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,1)
                   ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,1)
                   ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,1)
                   ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,1)

                   ATMOS_BOUNDARY_increment_DENS(k,i,j) = ( ATMOS_BOUNDARY_ref_DENS(k,i,j,2) &
                                                          - ATMOS_BOUNDARY_ref_DENS(k,i,j,1) ) &
                                                        / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                   ATMOS_BOUNDARY_increment_VELX(k,i,j) = ( ATMOS_BOUNDARY_ref_VELX(k,i,j,2) &
                                                          - ATMOS_BOUNDARY_ref_VELX(k,i,j,1) ) &
                                                        / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                   ATMOS_BOUNDARY_increment_VELY(k,i,j) = ( ATMOS_BOUNDARY_ref_VELY(k,i,j,2) &
                                                          - ATMOS_BOUNDARY_ref_VELY(k,i,j,1) ) &
                                                        / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                   ATMOS_BOUNDARY_increment_POTT(k,i,j) = ( ATMOS_BOUNDARY_ref_POTT(k,i,j,2) &
                                                          - ATMOS_BOUNDARY_ref_POTT(k,i,j,1) ) &
                                                        / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                   do iq = 1, BND_QA
                      ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1)

                      ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq) = ( ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,2)   &
                                                                - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1) ) &
                                                              / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                   end do
                end do
                end do
                end do

                last_updated = last_updated + ATMOS_BOUNDARY_UPDATE_DT
             else
                do j  = 1, JA
                do i  = 1, IA
                do k  = 1, KA
                   ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_increment_DENS(k,i,j)
                   ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) + ATMOS_BOUNDARY_increment_VELX(k,i,j)
                   ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) + ATMOS_BOUNDARY_increment_VELY(k,i,j)
                   ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) + ATMOS_BOUNDARY_increment_POTT(k,i,j)
                   do iq = 1, BND_QA
                      ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iq) + ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq)
                   end do
                end do
                end do
                end do
             endif
          endif

          if ( do_parent_process ) then !online [parent]

             handle = 1
             call ATMOS_BOUNDARY_update_online( DENS,MOMX,MOMY,RHOT,QTRC,handle )

          endif

       endif !online

       if ( .not. do_daughter_process ) then !offline or online-parent

          integrated_sec = TIME_NOWDAYSEC - last_updated
   
          if ( integrated_sec >= ATMOS_BOUNDARY_UPDATE_DT - EPS ) then
             boundary_timestep = boundary_timestep + 1
   
             call ATMOS_BOUNDARY_update_file
   
             do j  = 1, JA
             do i  = 1, IA
             do k  = 1, KA
                ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_ref_DENS(k,i,j,1)
                ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_ref_VELX(k,i,j,1)
                ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_ref_VELY(k,i,j,1)
                ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_ref_POTT(k,i,j,1)
   
                ATMOS_BOUNDARY_increment_DENS(k,i,j) = ( ATMOS_BOUNDARY_ref_DENS(k,i,j,2) &
                                                       - ATMOS_BOUNDARY_ref_DENS(k,i,j,1) ) &
                                                     / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                ATMOS_BOUNDARY_increment_VELX(k,i,j) = ( ATMOS_BOUNDARY_ref_VELX(k,i,j,2) &
                                                       - ATMOS_BOUNDARY_ref_VELX(k,i,j,1) ) &
                                                     / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                ATMOS_BOUNDARY_increment_VELY(k,i,j) = ( ATMOS_BOUNDARY_ref_VELY(k,i,j,2) &
                                                       - ATMOS_BOUNDARY_ref_VELY(k,i,j,1) ) &
                                                     / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                ATMOS_BOUNDARY_increment_POTT(k,i,j) = ( ATMOS_BOUNDARY_ref_POTT(k,i,j,2) &
                                                       - ATMOS_BOUNDARY_ref_POTT(k,i,j,1) ) &
                                                     / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                do iq = 1, BND_QA
                   ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1)
   
                   ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq) = ( ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,2) &
                                                             - ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1) ) &
                                                           / ( ATMOS_BOUNDARY_UPDATE_DT / TIME_DTSEC )
                end do
             end do
             end do
             end do
   
             last_updated = last_updated + ATMOS_BOUNDARY_UPDATE_DT
          else
             do j  = 1, JA
             do i  = 1, IA
             do k  = 1, KA
                ATMOS_BOUNDARY_DENS(k,i,j) = ATMOS_BOUNDARY_DENS(k,i,j) + ATMOS_BOUNDARY_increment_DENS(k,i,j)
                ATMOS_BOUNDARY_VELX(k,i,j) = ATMOS_BOUNDARY_VELX(k,i,j) + ATMOS_BOUNDARY_increment_VELX(k,i,j)
                ATMOS_BOUNDARY_VELY(k,i,j) = ATMOS_BOUNDARY_VELY(k,i,j) + ATMOS_BOUNDARY_increment_VELY(k,i,j)
                ATMOS_BOUNDARY_POTT(k,i,j) = ATMOS_BOUNDARY_POTT(k,i,j) + ATMOS_BOUNDARY_increment_POTT(k,i,j)
                do iq = 1, BND_QA
                   ATMOS_BOUNDARY_QTRC(k,i,j,iq) = ATMOS_BOUNDARY_QTRC(k,i,j,iq) + ATMOS_BOUNDARY_increment_QTRC(k,i,j,iq)
                end do
             end do
             end do
             end do
          endif

       endif !offline or online-parent

    else
       write(*,*) 'xxx [BUG] unsupported type'
       call PRC_MPIstop
    end if

    return
  end subroutine ATMOS_BOUNDARY_update

  !-----------------------------------------------------------------------------
  !> Update reference boundary from file
  subroutine ATMOS_BOUNDARY_update_file()
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: reference_atmos(KMAX,IMAX,JMAX) !> restart file (no HALO)

    character(len=H_LONG) :: bname

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    if (IO_L) write(IO_FID_LOG,*)"*** Atmos Boundary: read from boundary file(timestep=", boundary_timestep, ")"

    bname = ATMOS_BOUNDARY_IN_BASENAME

    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       ATMOS_BOUNDARY_ref_DENS(k,i,j,1) = ATMOS_BOUNDARY_ref_DENS(k,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(k,i,j,1) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(k,i,j,1) = ATMOS_BOUNDARY_ref_VELX(k,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(k,i,j,1) = ATMOS_BOUNDARY_ref_VELY(k,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(k,i,j,1) = ATMOS_BOUNDARY_ref_POTT(k,i,j,2)
       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,2)
       end do
    end do
    end do
    end do

    call FileRead( reference_atmos(:,:,:), bname, 'DENS', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_DENS(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELX', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_VELX(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'VELY', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_VELY(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    call FileRead( reference_atmos(:,:,:), bname, 'POTT', boundary_timestep, PRC_myrank )
    ATMOS_BOUNDARY_ref_POTT(KS:KE,IS:IE,JS:JE,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    do iq = 1, BND_QA
       call FileRead( reference_atmos(:,:,:), bname, AQ_NAME(iq), boundary_timestep, PRC_myrank )
       ATMOS_BOUNDARY_ref_QTRC(KS:KE,IS:IE,JS:JE,iq,2) = reference_atmos(1:KMAX,1:IMAX,1:JMAX)
    end do

    do j  = 1, JA
    do i  = 1, IA
       ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,2)

       ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,2)
       ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,2)
       ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,2)

       do iq = 1, BND_QA
          ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,2)
          ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,2)
       end do
    end do
    end do

    call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,2), 1 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2), 2 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,2), 3 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,2), 4 )
    call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 5 )

    do iq = 1, BND_QA
       call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 5+iq )
    end do

    call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,2), 1 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2), 2 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,2), 3 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,2), 4 )
    call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 5 )

    do iq = 1, BND_QA
       call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 5+iq )
    end do

    return
  end subroutine ATMOS_BOUNDARY_update_file


  !-----------------------------------------------------------------------------
  !> Update reference boundary by communicate with parent domain
  subroutine ATMOS_BOUNDARY_update_online ( &
       DENS,   & ! [in]
       MOMX,   & ! [in]
       MOMY,   & ! [in]
       RHOT,   & ! [in]
       QTRC,   & ! [in]
       handle  ) ! [in]
    use scale_process, only: &
       PRC_myrank
    use scale_grid_nest, only: &
       NEST_COMM_nestdown, &
       PARENT_KA,          &
       PARENT_IA,          &
       PARENT_JA,          &
       DAUGHTER_KA,        &
       DAUGHTER_IA,        &
       DAUGHTER_JA,        &
       PRNT_KS,            &
       PRNT_KE,            &
       PRNT_IS,            &
       PRNT_IE,            &
       PRNT_JS,            &
       PRNT_JE,            &
       DATR_KS,            &
       DATR_KE,            &
       DATR_IS,            &
       DATR_IE,            &
       DATR_JS,            &
       DATR_JE
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: MOMX(KA,IA,JA)
    real(RP), intent(in) :: MOMY(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    integer,  intent(in) :: handle

    real(RP) :: dummy1_p(PARENT_KA(handle),  PARENT_IA(handle),  PARENT_JA(handle))
    real(RP) :: dummy1_d(DAUGHTER_KA(handle),DAUGHTER_IA(handle),DAUGHTER_JA(handle),2)
    real(RP) :: dummy2_p(PARENT_KA(handle),  PARENT_IA(handle),  PARENT_JA(handle),  BND_QA)
    real(RP) :: dummy2_d(DAUGHTER_KA(handle),DAUGHTER_IA(handle),DAUGHTER_JA(handle),BND_QA,2)

    integer :: i, j, k, iq
    !---------------------------------------------------------------------------

    dummy1_p(:,:,:)     = 0.0_RP
    dummy1_d(:,:,:,:)   = 0.0_RP
    dummy2_p(:,:,:,:)   = 0.0_RP
    dummy2_d(:,:,:,:,:) = 0.0_RP

    if ( handle == 1 .and. do_parent_process ) then ! [parent]
       if (IO_L) write(IO_FID_LOG,*)"*** NESTCOMM inter-domain as a PARENT"

       call NEST_COMM_nestdown( handle,                    &
                                BND_QA,                    &
                                DENS(:,:,:),               &  !(KA,IA,JA)
                                MOMX(:,:,:),               &  !(KA,IA,JA)
                                MOMY(:,:,:),               &  !(KA,IA,JA)
                                RHOT(:,:,:),               &  !(KA,IA,JA)
                                QTRC(:,:,:,1:BND_QA),      &  !(KA,IA,JA,QA)
                                dummy1_d(:,:,:,2),         &  !(KA,IA,JA)
                                dummy1_d(:,:,:,2),         &  !(KA,IA,JA)
                                dummy1_d(:,:,:,2),         &  !(KA,IA,JA)
                                dummy1_d(:,:,:,2),         &  !(KA,IA,JA)
                                dummy2_d(:,:,:,1:BND_QA,2) )  !(KA,IA,JA,QA)

    elseif ( handle == 2 .and. do_daughter_process ) then ! [daughter]
       if (IO_L) write(IO_FID_LOG,*)"*** NESTCOMM inter-domain as a DAUGHTER (step=", boundary_timestep, ")"

       do j  = 1, JA
       do i  = 1, IA
       do k  = 1, KA
          ATMOS_BOUNDARY_ref_DENS(k,i,j,1) = ATMOS_BOUNDARY_ref_DENS(k,i,j,2)
          ATMOS_BOUNDARY_ref_VELZ(k,i,j,1) = ATMOS_BOUNDARY_ref_VELZ(k,i,j,2)
          ATMOS_BOUNDARY_ref_VELX(k,i,j,1) = ATMOS_BOUNDARY_ref_VELX(k,i,j,2)
          ATMOS_BOUNDARY_ref_VELY(k,i,j,1) = ATMOS_BOUNDARY_ref_VELY(k,i,j,2)
          ATMOS_BOUNDARY_ref_POTT(k,i,j,1) = ATMOS_BOUNDARY_ref_POTT(k,i,j,2)
          do iq = 1, BND_QA
             ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,1) = ATMOS_BOUNDARY_ref_QTRC(k,i,j,iq,2)
          end do
       end do
       end do
       end do

       call NEST_COMM_nestdown( handle,                                 &
                                BND_QA,                                 &
                                dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                                dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                                dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                                dummy1_p(:,:,:),                        &   !(KA,IA,JA)
                                dummy2_p(:,:,:,1:BND_QA),               &   !(KA,IA,JA,QA)
                                ATMOS_BOUNDARY_ref_DENS(:,:,:,2),       &   !(KA,IA,JA)
                                ATMOS_BOUNDARY_ref_VELX(:,:,:,2),       &   !(KA,IA,JA)
                                ATMOS_BOUNDARY_ref_VELY(:,:,:,2),       &   !(KA,IA,JA)
                                ATMOS_BOUNDARY_ref_POTT(:,:,:,2),       &   !(KA,IA,JA)
                                ATMOS_BOUNDARY_ref_QTRC(:,:,:,1:BND_QA,2) ) !(KA,IA,JA,QA)

       do j  = 1, JA
       do i  = 1, IA
          ATMOS_BOUNDARY_ref_DENS(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_DENS(KS,i,j,2)
          ATMOS_BOUNDARY_ref_VELZ(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KS,i,j,2)
          ATMOS_BOUNDARY_ref_VELX(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELX(KS,i,j,2)
          ATMOS_BOUNDARY_ref_VELY(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_VELY(KS,i,j,2)
          ATMOS_BOUNDARY_ref_POTT(   1:KS-1,i,j,2) = ATMOS_BOUNDARY_ref_POTT(KS,i,j,2)

          ATMOS_BOUNDARY_ref_DENS(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_DENS(KE,i,j,2)
          ATMOS_BOUNDARY_ref_VELZ(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELZ(KE,i,j,2)
          ATMOS_BOUNDARY_ref_VELX(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELX(KE,i,j,2)
          ATMOS_BOUNDARY_ref_VELY(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_VELY(KE,i,j,2)
          ATMOS_BOUNDARY_ref_POTT(KE+1:KA,  i,j,2) = ATMOS_BOUNDARY_ref_POTT(KE,i,j,2)

          do iq = 1, BND_QA
             ATMOS_BOUNDARY_ref_QTRC(   1:KS-1,i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KS,i,j,iq,2)
             ATMOS_BOUNDARY_ref_QTRC(KE+1:KA,  i,j,iq,2) = ATMOS_BOUNDARY_ref_QTRC(KE,i,j,iq,2)
          end do
       end do
       end do

       call COMM_vars8( ATMOS_BOUNDARY_ref_DENS(:,:,:,2), 1 )
       call COMM_vars8( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2), 2 )
       call COMM_vars8( ATMOS_BOUNDARY_ref_VELX(:,:,:,2), 3 )
       call COMM_vars8( ATMOS_BOUNDARY_ref_VELY(:,:,:,2), 4 )
       call COMM_vars8( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 5 )

       do iq = 1, BND_QA
          call COMM_vars8( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 5+iq )
       end do

       call COMM_wait ( ATMOS_BOUNDARY_ref_DENS(:,:,:,2), 1 )
       call COMM_wait ( ATMOS_BOUNDARY_ref_VELZ(:,:,:,2), 2 )
       call COMM_wait ( ATMOS_BOUNDARY_ref_VELX(:,:,:,2), 3 )
       call COMM_wait ( ATMOS_BOUNDARY_ref_VELY(:,:,:,2), 4 )
       call COMM_wait ( ATMOS_BOUNDARY_ref_POTT(:,:,:,2), 5 )

       do iq = 1, BND_QA
          call COMM_wait ( ATMOS_BOUNDARY_ref_QTRC(:,:,:,iq,2), 5+iq )
       end do

    endif

    return
  end subroutine ATMOS_BOUNDARY_update_online

end module scale_atmos_boundary
