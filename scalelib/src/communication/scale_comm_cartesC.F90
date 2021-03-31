!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI Communication module for Cartesian C-grid
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_comm_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_tracer

  use scale_prc, only: &
     PRC_abort
  use scale_prc_cartesC, only: &
     PRC_next, &
     PRC_W,    &
     PRC_N,    &
     PRC_E,    &
     PRC_S,    &
     PRC_NW,   &
     PRC_NE,   &
     PRC_SW,   &
     PRC_SE,   &
     PRC_HAS_W, &
     PRC_HAS_N, &
     PRC_HAS_E, &
     PRC_HAS_S
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COMM_setup
  public :: COMM_regist
  public :: COMM_finalize
  public :: COMM_vars_init
  public :: COMM_vars8_init
  public :: COMM_vars
  public :: COMM_vars8
  public :: COMM_wait
  public :: COMM_gather
  public :: COMM_bcast

  interface COMM_vars
     module procedure COMM_vars_2D
     module procedure COMM_vars_3D
  end interface COMM_vars

  interface COMM_vars8
     module procedure COMM_vars8_2D
     module procedure COMM_vars8_3D
  end interface COMM_vars8

  interface COMM_wait
     module procedure COMM_wait_2D
     module procedure COMM_wait_3D
  end interface COMM_WAIT

  interface COMM_gather
     module procedure COMM_gather_2D
     module procedure COMM_gather_3D
  end interface COMM_gather

  interface COMM_bcast
     module procedure COMM_bcast_SCR
     module procedure COMM_bcast_1D
     module procedure COMM_bcast_2D
     module procedure COMM_bcast_3D
     module procedure COMM_bcast_4D
     module procedure COMM_bcast_INT_SCR
     module procedure COMM_bcast_INT_1D
     module procedure COMM_bcast_INT_2D
     module procedure COMM_bcast_LOGICAL_SCR
     module procedure COMM_bcast_CHARACTER
  end interface COMM_bcast

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: COMM_datatype   !< datatype of variable
  integer, public :: COMM_world      !< communication world ID

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private              :: COMM_vsize_max      !< # limit of communication variables at once
  integer,  private              :: COMM_vsize_max_pc      !< # limit of total communication variables for MPI PC

  logical,  private              :: COMM_IsAllPeriodic  !< periodic boundary condition?

  logical,  private              :: COMM_USE_MPI_PC = .true.

  type ginfo_t
     integer           :: KA
     integer           :: IA, IS, IE, IHALO
     integer           :: JA, JS, JE, JHALO
     integer           :: nreq_max          !< # limit of communication request at once
     integer           :: size2D_NS4        !< 2D data size (W/E    HALO, 4/8-direction comm.)
     integer           :: size2D_NS8        !< 2D data size (N/S    HALO,   4-direction comm.)
     integer           :: size2D_WE         !< 2D data size (N/S    HALO,   8-direction comm.)
     integer           :: size2D_4C         !< 2D data size (corner HALO,   8-direction comm.)
     integer           :: vars_id = 0       !< id of variables for persistent comm.
     real(RP), pointer :: recvpack_W2P(:,:) !< packing packet (receive, from W)
     real(RP), pointer :: recvpack_E2P(:,:) !< packing packet (receive, from E)
     real(RP), pointer :: sendpack_P2W(:,:) !< packing packet (send,    to W  )
     real(RP), pointer :: sendpack_P2E(:,:) !< packing packet (send,    to E  )
     integer,  pointer :: req_cnt (:)       !< request ID of each MPI send/recv
     integer,  pointer :: req_list(:,:)     !< request ID set of each variables
     integer,  pointer :: preq_cnt (:)      !< request ID of each MPI PC
     integer,  pointer :: preq_list(:,:)    !< request ID set of each variables for MPI PC
     integer,  pointer :: pseqid(:)         !< sequential ID of each variables for MPI PC
#ifdef DEBUG
     logical,  pointer :: use_packbuf(:)    !< using flag for packing buffer
#endif


  end type ginfo_t

  integer, private, parameter :: COMM_gid_max = 20
  integer, private            :: COMM_gid
  type(ginfo_t), private      :: ginfo(COMM_gid_max)

  logical, private :: initialized = .false.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine COMM_setup
    use scale_prc, only: &
       PRC_LOCAL_COMM_WORLD
    implicit none

    namelist / PARAM_COMM_CARTESC / &
       COMM_vsize_max, &
       COMM_vsize_max_pc, &
       COMM_USE_MPI_PC

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("COMM_setup",*) 'Setup'

    COMM_vsize_max = max( 10 + QA*2, 25 )
    COMM_vsize_max_pc = 50 + QA*2

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("COMM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("COMM_setup",*) 'Not appropriate names in namelist PARAM_COMM_CARTESC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_COMM_CARTESC)

    if ( PRC_HAS_N .AND. PRC_HAS_S .AND. PRC_HAS_W .AND. PRC_HAS_E ) then
       COMM_IsAllPeriodic = .true.
    else
       COMM_IsAllPeriodic = .false.
    endif

    if ( RP == kind(0.D0) ) then
       COMM_datatype = MPI_DOUBLE_PRECISION
    elseif( RP == kind(0.0) ) then
       COMM_datatype = MPI_REAL
    else
       LOG_ERROR("COMM_setup",*) 'precision is not supportd'
       call PRC_abort
    endif

    COMM_world = PRC_LOCAL_COMM_WORLD

    LOG_NEWLINE
    LOG_INFO("COMM_setup",*) 'Communication information'
    LOG_INFO_CONT(*)         'Maximum number of vars for one communication: ', COMM_vsize_max
    LOG_INFO_CONT(*)         'All side is periodic?                       : ', COMM_IsAllPeriodic

    initialized = .true.

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  !> Regist grid
  subroutine COMM_regist( &
       KA, IA, JA, IHALO, JHALO, &
       gid )
    implicit none

    integer, intent(in)  :: KA, IA, JA, IHALO, JHALO
    integer, intent(out) :: gid

    integer :: IMAX, JMAX
    integer :: nreq_NS, nreq_WE, nreq_4C

    if ( .not. initialized ) then
       LOG_ERROR("COMM_regist",*) 'COMM_setup must be called before calling COMM_regist'
       call PRC_abort
    end if

    COMM_gid = COMM_gid + 1
    if ( COMM_gid > COMM_gid_max ) then
       LOG_ERROR("COMM_regist",*) 'number of registed grid size exceeds the limit'
       call PRC_abort
    end if
    gid = COMM_gid

    if ( IA < IHALO * 3 ) then
       LOG_ERROR("COMM_regist",*) 'IA must be >= IHALO * 3'
       call PRC_abort
    end if
    if ( JA < JHALO * 3 ) then
       LOG_ERROR("COMM_regist",*) 'JA must be >= JHALO * 3'
       call PRC_abort
    end if

    IMAX = IA - IHALO * 2
    JMAX = JA - JHALO * 2

    ginfo(gid)%KA    = KA
    ginfo(gid)%IA    = IA
    ginfo(gid)%IS    = IHALO + 1
    ginfo(gid)%IE    = IA - IHALO
    ginfo(gid)%IHALO = IHALO
    ginfo(gid)%JA    = JA
    ginfo(gid)%JS    = JHALO + 1
    ginfo(gid)%JE    = JA - JHALO
    ginfo(gid)%JHALO = JHALO

    nreq_NS  = 2 * JHALO !--- send x JHALO, recv x JHALO
    nreq_WE  = 2         !--- send x 1    , recv x 1
    nreq_4C  = 2 * JHALO !--- send x JHALO, recv x JHALO

    if ( COMM_USE_MPI_PC ) then
       ginfo(gid)%nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C + 1
    else
       ginfo(gid)%nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C
    end if

    ginfo(gid)%size2D_NS4 = IA   * JHALO
    ginfo(gid)%size2D_NS8 = IMAX
    ginfo(gid)%size2D_WE  = JMAX * IHALO
    ginfo(gid)%size2D_4C  =        IHALO

    allocate( ginfo(gid)%recvpack_W2P(ginfo(gid)%size2D_WE * KA, COMM_vsize_max) )
    allocate( ginfo(gid)%recvpack_E2P(ginfo(gid)%size2D_WE * KA, COMM_vsize_max) )
    allocate( ginfo(gid)%sendpack_P2W(ginfo(gid)%size2D_WE * KA, COMM_vsize_max) )
    allocate( ginfo(gid)%sendpack_P2E(ginfo(gid)%size2D_WE * KA, COMM_vsize_max) )


#ifdef DEBUG
    allocate( ginfo(gid)%use_packbuf(COMM_vsize_max) )
    ginfo(gid)%use_packbuf(:) = .false.
#endif

    allocate( ginfo(gid)%req_cnt (                     COMM_vsize_max) )
    allocate( ginfo(gid)%req_list(ginfo(gid)%nreq_MAX, COMM_vsize_max) )
    ginfo(gid)%req_cnt (:)   = -1
    ginfo(gid)%req_list(:,:) = MPI_REQUEST_NULL

    if ( COMM_USE_MPI_PC ) then
       ginfo(gid)%vars_id = 0
       allocate( ginfo(gid)%preq_cnt (                      COMM_vsize_max_pc) )
       allocate( ginfo(gid)%preq_list(ginfo(gid)%nreq_MAX+1,COMM_vsize_max_pc) )
       ginfo(gid)%preq_cnt (:)   = -1
       ginfo(gid)%preq_list(:,:) = MPI_REQUEST_NULL

       allocate( ginfo(gid)%pseqid(COMM_vsize_max_pc) )
    end if

    LOG_NEWLINE
    LOG_INFO("COMM_regist",*) 'Register grid: id=', gid
    LOG_INFO_CONT(*)          'Data size of var (3D,including halo) [byte] : ', RP*KA*IA*JA
    LOG_INFO_CONT(*)          'Data size of halo                    [byte] : ', RP*KA*(2*IA*JHALO+2*JMAX*IHALO)
    LOG_INFO_CONT(*)          'Ratio of halo against the whole 3D grid     : ', real(2*IA*JHALO+2*JMAX*IHALO) / real(IA*JA)

    return
  end subroutine COMM_regist

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine COMM_finalize
    use mpi
    implicit none

    integer :: gid
    integer :: i, j, ierr
    !---------------------------------------------------------------------------

    do gid = 1, COMM_gid
       deallocate( ginfo(gid)%recvpack_W2P )
       deallocate( ginfo(gid)%recvpack_E2P )
       deallocate( ginfo(gid)%sendpack_P2W )
       deallocate( ginfo(gid)%sendpack_P2E )
#ifdef DEBUG
       deallocate( ginfo(gid)%use_packbuf )
#endif

       deallocate( ginfo(gid)%req_cnt )
       deallocate( ginfo(gid)%req_list )

       if ( COMM_USE_MPI_PC ) then
          do j = 1, COMM_vsize_max_pc
             do i = 1, ginfo(gid)%nreq_MAX+1
                if (ginfo(gid)%preq_list(i,j) .NE. MPI_REQUEST_NULL) &
                     call MPI_REQUEST_FREE(ginfo(gid)%preq_list(i,j), ierr)
             enddo
          enddo
          deallocate( ginfo(gid)%preq_cnt )
          deallocate( ginfo(gid)%preq_list )
          deallocate( ginfo(gid)%pseqid )
          ginfo(gid)%vars_id = 0
       end if

    end do

    COMM_gid = 0

    initialized = .false.

    return
  end subroutine COMM_finalize

  !-----------------------------------------------------------------------------
  !> Register variables
  subroutine COMM_vars_init( &
       varname, &
       var,     &
       vid,     &
       gid      )
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name
    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID

    integer,          intent(in), optional :: gid

    integer :: gid_
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_vars_init",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    if ( vid > COMM_vsize_max ) then
       LOG_ERROR("COMM_vars_init",*) 'vid exceeds max', vid, COMM_vsize_max, gid
       call PRC_abort
    end if

    if ( COMM_USE_MPI_PC ) then

       ginfo(gid_)%vars_id = ginfo(gid_)%vars_id + 1
       if ( ginfo(gid_)%vars_id > COMM_vsize_max_pc ) then
          LOG_ERROR("COMM_vars_init",*) 'number of variable for MPI PC exceeds max', ginfo(gid_)%vars_id, COMM_vsize_max_pc
          call PRC_abort
       end if

       call PROF_rapstart('COMM_init_pers', 2)
       call vars_init_mpi_pc(var, gid_, ginfo(gid_)%vars_id, vid)
       call PROF_rapend  ('COMM_init_pers', 2)

       vid = ginfo(gid_)%vars_id + COMM_vsize_max

       LOG_INFO("COMM_vars_init",'(1x,A,I3.3,A,I3.3,2A)') 'Initialize variable (grid ID = ', gid_, '): ID = ', vid, &
                                                                                       ', name = ', trim(varname)

    end if

    return
  end subroutine COMM_vars_init

  !-----------------------------------------------------------------------------
  !> Register variables
  subroutine COMM_vars8_init( &
       varname, &
       var,     &
       vid,     &
       gid      )
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name

    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID

    integer,          intent(in), optional :: gid

    integer :: gid_
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_vars8_init",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    if ( vid > COMM_vsize_max ) then
       LOG_ERROR("COMM_vars8_init",*) 'vid exceeds max', vid, COMM_vsize_max
       call PRC_abort
    end if

    if ( COMM_USE_MPI_PC ) then

       ginfo(gid_)%vars_id = ginfo(gid_)%vars_id + 1
       if ( ginfo(gid_)%vars_id > COMM_vsize_max_pc ) then
          LOG_ERROR("COMM_vars8_init",*) 'number of variable for MPI PC exceeds max', ginfo(gid_)%vars_id, COMM_vsize_max_pc
          call PRC_abort
       end if

       call PROF_rapstart('COMM_init_pers', 2)
       call vars8_init_mpi_pc(var, gid_, ginfo(gid_)%vars_id, vid)
       call PROF_rapend  ('COMM_init_pers', 2)

       vid = ginfo(gid_)%vars_id + COMM_vsize_max

       LOG_INFO("COMM_vars8_init",'(1x,A,I3.3,A,I3.3,2A)') 'Initialize variable (grid ID = ', gid_, '): ID = ', vid, &
                                                                                       ', name = ', trim(varname)

    end if

    return
  end subroutine COMM_vars8_init

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_3D(var, vid, gid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication

    integer,  intent(in)    :: vid        !< request ID

    integer,  intent(in), optional :: gid

    integer :: gid_
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_vars_3D",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    if ( vid > COMM_vsize_max ) then
       call PROF_rapstart('COMM_vars_pers', 2)
       call vars_3D_mpi_pc(var, gid_, vid-COMM_vsize_max)
       call PROF_rapend  ('COMM_vars_pers', 2)
    else
       call PROF_rapstart('COMM_vars', 2)
       call vars_3D_mpi(var, gid_, vid)
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_3D(var, vid, gid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)

    integer,  intent(in)    :: vid

    integer,  intent(in), optional :: gid

    integer :: gid_
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_vars8_3D",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    if ( vid > COMM_vsize_max ) then
       call PROF_rapstart('COMM_vars_pers', 2)
       call vars_3D_mpi_pc(var, gid_, vid-COMM_vsize_max)
       call PROF_rapend  ('COMM_vars_pers', 2)
    else
       call PROF_rapstart('COMM_vars', 2)
       call vars8_3D_mpi(var, gid_, vid)
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars8_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_3D(var, vid, FILL_BND, gid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)

    integer, intent(in)     :: vid

    logical, intent(in), optional :: FILL_BND
    integer, intent(in), optional :: gid

    logical :: FILL_BND_
    integer :: gid_
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_wait_3D",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    if ( vid > COMM_vsize_max ) then
       call PROF_rapstart('COMM_wait_pers', 2)
       call wait_3D_mpi_pc(var, gid_, vid-COMM_vsize_max)
       call PROF_rapend  ('COMM_wait_pers', 2)
    else
       call PROF_rapstart('COMM_wait', 2)
       call wait_3D_mpi(var, gid_, vid)
       call PROF_rapend  ('COMM_wait', 2)
    end if

    ! copy inner data to boundary
    if ( .NOT. COMM_IsAllPeriodic ) then
       if ( FILL_BND_ ) then
          call copy_boundary_3D(var, gid_)
       end if
    end if

    return
  end subroutine COMM_wait_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_2D(var, vid, gid)
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer,  intent(in)    :: vid

    integer,  intent(in), optional :: gid

    integer :: gid_
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_vars_2D",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    call PROF_rapstart('COMM_vars', 2)
    call vars_2D_mpi(var, gid_, vid)
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_2D(var, vid, gid)
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer,  intent(in)    :: vid

    integer,  intent(in), optional :: gid

    integer :: gid_
    !---------------------------------------------------------------------------

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_vars8_2D",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    call PROF_rapstart('COMM_vars', 2)
    call vars8_2D_mpi(var, gid_, vid)
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars8_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_2D(var, vid, FILL_BND, gid)
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer,  intent(in)    :: vid

    logical,  intent(in), optional :: FILL_BND
    integer,  intent(in), optional :: gid

    logical :: FILL_BND_
    integer :: gid_
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    gid_ = 1
    if ( present(gid) ) gid_ = gid
    if ( gid_ > COMM_gid_max ) then
       LOG_ERROR("COMM_wait_2D",*) 'gid is invalid', gid_, COMM_gid_max
       call PRC_abort
    end if

    call PROF_rapstart('COMM_wait', 2)
    call wait_2D_mpi(var, gid_, vid)
    call PROF_rapend  ('COMM_wait', 2)

    if( .NOT. COMM_IsAllPeriodic ) then
       if ( FILL_BND_ ) then
          call copy_boundary_2D(var, gid_)
       end if
    end if

    return
  end subroutine COMM_wait_2D

  !-----------------------------------------------------------------------------
  !> calculate horizontal mean (global total with communication) 2D
  subroutine COMM_horizontal_mean_2D( &
       IA, IS, IE, JA, JS, JE, &
       var,    &
       varmean )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: var(IA,JA) !< 2D value

    real(RP), intent(out) :: varmean  !< horizontal mean

    real(DP) :: stat(2)
    real(DP) :: allstat(2)
    real(DP) :: zerosw

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    stat(:) = 0.0_DP
    do j = JS, JE
    do i = IS, IE
       if ( abs(var(i,j)) < abs(CONST_UNDEF) ) then
          stat(1) = stat(1) + var(i,j)
          stat(2) = stat(2) + 1.0_DP
       endif
    enddo
    enddo

    ! All reduce
    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    call MPI_Allreduce( stat,                 &
                        allstat,              &
                        2,                    &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM,              &
                        COMM_world,           &
                        ierr                  )
    call PROF_rapend  ('COMM_Allreduce', 2)

    zerosw = 0.5_DP - sign(0.5_DP, allstat(1) - 1.E-12_DP )
    varmean = allstat(1) / ( allstat(2) + zerosw ) * ( 1.0_DP - zerosw )
    !LOG_INFO("COMM_horizontal_mean_2D",*) varmean, allstat(1), allstat(2)

    return
  end subroutine COMM_horizontal_mean_2D

  !-----------------------------------------------------------------------------
  !> calculate horizontal mean (global total with communication) 3D
  subroutine COMM_horizontal_mean_3D( &
       KA, IA, IS, IE, JA, JS, JE, &
       var,    &
       varmean )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer,  intent(in)  :: KA
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: var(KA,IA,JA) !< 3D value

    real(RP), intent(out) :: varmean(KA)   !< horizontal mean

    real(DP) :: stat   (KA,2)
    real(DP) :: allstat(KA,2)
    real(DP) :: zerosw

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    stat(:,:) = 0.0_DP
    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       if ( abs(var(k,i,j)) < abs(CONST_UNDEF) ) then
          stat(k,1) = stat(k,1) + var(k,i,j)
          stat(k,2) = stat(k,2) + 1.0_DP
       endif
    enddo
    enddo
    enddo

    ! All reduce
    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    call MPI_Allreduce( stat,                 &
                        allstat,              &
                        KA * 2,               &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM,              &
                        COMM_world,           &
                        ierr                  )
    call PROF_rapend  ('COMM_Allreduce', 2)

    do k = 1, KA
       zerosw = 0.5_DP - sign(0.5_DP, allstat(k,2) - 1.E-12_DP )
       varmean(k) = allstat(k,1) / ( allstat(k,2) + zerosw ) * ( 1.0_DP - zerosw )
       !LOG_INFO("COMM_horizontal_mean_3D",*) k, varmean(k), allstatval(k), allstatcnt(k)
    enddo

    return
  end subroutine COMM_horizontal_mean_3D

  !-----------------------------------------------------------------------------
  !> Get data from whole process value in 2D field
  subroutine COMM_gather_2D( &
       IA, JA, &
       send, &
       recv  )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)  :: IA, JA      !< dimension size
    real(RP), intent(in)  :: send(IA,JA) !< send buffer

    real(RP), intent(out) :: recv(:,:,:) !< receive buffer (IA,JA,nprcs)

    integer :: sendcounts, recvcounts
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = IA * JA
    recvcounts = IA * JA

    call MPI_GATHER( send(:,:),      &
                     sendcounts,     &
                     COMM_datatype,  &
                     recv(:,:,:),    &
                     recvcounts,     &
                     COMM_datatype,  &
                     PRC_masterrank, &
                     COMM_world,     &
                     ierr            )

    return
  end subroutine COMM_gather_2D

  !-----------------------------------------------------------------------------
  !> Get data from whole process value in 3D field
  subroutine COMM_gather_3D( &
       KA, IA, JA, &
       send, &
       recv  )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)  :: KA, IA, JA     !< dimension size
    real(RP), intent(in)  :: send(KA,IA,JA) !< send buffer

    real(RP), intent(out) :: recv(:,:,:,:) !< receive buffer(KA,IA,JA,nprcs)

    integer :: sendcounts, recvcounts
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = KA * IA * JA
    recvcounts = KA * IA * JA

    call MPI_GATHER( send(:,:,:),    &
                     sendcounts,     &
                     COMM_datatype,  &
                     recv(:,:,:,:),  &
                     recvcounts,     &
                     COMM_datatype,  &
                     PRC_masterrank, &
                     COMM_world,     &
                     ierr            )

    return
  end subroutine COMM_gather_3D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar field
  subroutine COMM_bcast_SCR( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(inout) :: var  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    COMM_datatype,  &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_SCR

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D field
  subroutine COMM_bcast_1D( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA       !< dimension size

    real(RP), intent(inout) :: var(IA)  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    call MPI_BCAST( var(:),         &
                    counts,         &
                    COMM_datatype,  &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_1D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 2D field
  subroutine COMM_bcast_2D( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA, JA     !< dimension size

    real(RP), intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    COMM_datatype,  &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_2D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 3D field
  subroutine COMM_bcast_3D( KA, IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA    !< dimension size

    real(RP), intent(inout) :: var(KA,IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA

    call MPI_BCAST( var(:,:,:),     &
                    counts,         &
                    COMM_datatype,  &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_3D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 4D field
  subroutine COMM_bcast_4D( KA, IA, JA, NT, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA, NT   !< dimension size

    real(RP), intent(inout) :: var(KA,IA,JA,NT) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA * NT
    if ( KA>0 .AND. IA>0 .AND. JA>0 .AND. NT>0 .AND. &
         counts < 0 ) then
       LOG_ERROR("COMM_bcast_4D",*) 'counts overflow'
       call PRC_abort
    end if

    call MPI_BCAST( var(:,:,:,:),   &
                    counts,         &
                    COMM_datatype,  &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_4D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar (integer)
  subroutine COMM_bcast_INT_SCR( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(inout) :: var !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_INT_SCR

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D field (integer)
  subroutine COMM_bcast_INT_1D( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(in)    :: IA      !< dimension size
    integer, intent(inout) :: var(IA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_INT_1D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 2D field (integer)
  subroutine COMM_bcast_INT_2D( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(in)    :: IA, JA     !< dimension size

    integer, intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_INT_2D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar (logical)
  subroutine COMM_bcast_LOGICAL_SCR( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    logical, intent(inout) :: var   !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_LOGICAL,    &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_LOGICAL_SCR

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in character
  subroutine COMM_bcast_CHARACTER( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    character(len=*), intent(inout) :: var   !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = len(var)

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_CHARACTER,  &
                    PRC_masterrank, &
                    COMM_world,     &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_CHARACTER

!-------------------------------------------------------------------------------
! private routines
!-------------------------------------------------------------------------------
  subroutine vars_init_mpi_pc(var, gid, vid, seqid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: gid
    integer,  intent(in) :: vid
    integer,  intent(in) :: seqid

    integer :: ireq, tag, ierr
    logical :: flag

    integer :: KA
    integer :: JA, JS, JE, JHALO

    integer :: nreq
    integer :: i

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

    KA    = ginfo(gid)%KA
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype,                 &
                        MPI_PROC_NULL, tag+ginfo(gid)%nreq_max+1, COMM_world, &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,vid), ierr )

    !--- From 4-Direction HALO communicate
    ! From S
    call MPI_RECV_INIT( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,                &
                        PRC_next(PRC_S), tag+1, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! From N
    call MPI_RECV_INIT( var(:,:,JE+1:JA), ginfo(gid)%size2D_NS4*KA, COMM_datatype,               &
                        PRC_next(PRC_N), tag+2, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    if ( .not. PRC_TwoD ) then
       ! From E
       call MPI_RECV_INIT( ginfo(gid)%recvpack_E2P(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,&
                           PRC_next(PRC_E), tag+3, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From W
       call MPI_RECV_INIT( ginfo(gid)%recvpack_W2P(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                           PRC_next(PRC_W), tag+4, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr  )
       ireq = ireq + 1
    end if

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
       ! To W HALO
       call MPI_SEND_INIT( ginfo(gid)%sendpack_P2W(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                           PRC_next(PRC_W), tag+3, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr  )
       ireq = ireq + 1
       ! To E HALO
       call MPI_SEND_INIT( ginfo(gid)%sendpack_P2E(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                           PRC_next(PRC_E), tag+4, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr  )
       ireq = ireq + 1
    end if
    ! To N HALO
    call MPI_SEND_INIT( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype,         &
                        PRC_next(PRC_N), tag+1, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! To S HALO
    call MPI_SEND_INIT( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,         &
                        PRC_next(PRC_S), tag+2, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
    ireq = ireq + 1

    ginfo(gid)%preq_cnt(vid) = ireq - 1
    ginfo(gid)%pseqid(vid) = seqid

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(vid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,vid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    return
  end subroutine vars_init_mpi_pc

  subroutine vars8_init_mpi_pc(var, gid, vid, seqid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: gid
    integer,  intent(in) :: vid
    integer,  intent(in) :: seqid

    integer :: ireq, tag, tagc
    integer :: ierr
    logical :: flag

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JA, JS, JE, JHALO

    integer :: nreq
    integer :: i, j

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype,                 &
                        MPI_PROC_NULL, tag+ginfo(gid)%nreq_max+1, COMM_world, &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,vid), ierr )


    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          tagc = 0
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                       &
                                 PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From SW
          tagc = 10
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                           &
                                  PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NE
          tagc = 20
          do j = JE+1, JA
             call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                       &
                                 PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NW
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                          &
                                 PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From E
          tagc = 60
          call MPI_RECV_INIT( ginfo(gid)%recvpack_E2P(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                              PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! From W
          tagc = 70
          call MPI_RECV_INIT( ginfo(gid)%recvpack_W2P(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                              PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
       end if
       ! From S
       tagc = 40
       do j = 1, JS-1
          call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                       &
                              PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JA
          call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                       &
                              PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !--- To 8-Direction HALO communicate
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                   &
                          PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                   &
                          PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       if ( .not. PRC_TwoD ) then
          ! To W HALO
          tagc = 60
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2W(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                              PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! To E HALO
          tagc = 70
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2E(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                              PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! To NW HALO
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                 PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To NE HALO
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                 &
                                 PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SW HALO
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                 PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SE HALO
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                 &
                                 PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       end if

    else ! non-periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                       &
                                    PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                    PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                    PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From SW
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                          &
                                    PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NE
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                       &
                                    PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                    PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                    PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NW
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                          &
                                    PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From E
          if ( PRC_HAS_E ) then
             tagc = 60
             call MPI_RECV_INIT( ginfo(gid)%recvpack_E2P(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                                 PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             tagc = 70
             call MPI_RECV_INIT( ginfo(gid)%recvpack_W2P(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                                 PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
       end if
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                   &
                             PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JA
             call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                   &
                             PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                   &
                             PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                   &
                             PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       if ( .not. PRC_TwoD ) then
          ! To W HALO
          if ( PRC_HAS_W ) then
             tagc = 60
             call MPI_SEND_INIT( ginfo(gid)%sendpack_P2W(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                                 PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             tagc = 70
             call MPI_SEND_INIT( ginfo(gid)%sendpack_P2E(:,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype,   &
                                 PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                    PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                 &
                                    PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                    PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                &
                                    PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                    PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                    PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                 &
                                    PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                    PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                &
                                    PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
       end if

    endif

    ginfo(gid)%preq_cnt(vid) = ireq - 1
    ginfo(gid)%pseqid(vid) = seqid

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(vid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,vid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    return
  end subroutine vars8_init_mpi_pc

  subroutine vars_3D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication
    integer,  intent(in)    :: gid        !< grid ID
    integer,  intent(in)    :: vid        !< request ID


    integer :: ireq, tag

    integer :: KA
    integer :: JA, JS, JE, JHALO

    integer :: ierr
    !---------------------------------------------------------------------------

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

    KA    = ginfo(gid)%KA
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(vid) ) then
       LOG_ERROR("vars_3D_mpi",*) 'packing buffer is already used', vid
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 4-Direction HALO communicate
       ! From S
       call MPI_IRECV( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,               &
                       PRC_next(PRC_S), tag+1, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From N
       call MPI_IRECV( var(:,:,JE+1:1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,               &
                       PRC_next(PRC_N), tag+2, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
       if ( .not. PRC_TwoD ) then
          ! From E
          call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                          PRC_next(PRC_E), tag+3, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! From W
          call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                          PRC_next(PRC_W), tag+4, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       end if

       !--- To 4-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          call pack_3D(var, gid, vid)

          ! To W HALO
          call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                          PRC_next(PRC_W), tag+3, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! To E HALO
          call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                          PRC_next(PRC_E), tag+4, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       end if
       ! To N HALO
       call MPI_ISEND( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype,        &
                       PRC_next(PRC_N), tag+1, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To S HALO
       call MPI_ISEND( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,        &
                       PRC_next(PRC_S), tag+2, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1

    else ! non-periodic condition

       !--- From 4-Direction HALO communicate
       ! From S
       if ( PRC_HAS_S ) then
          call MPI_IRECV( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,               &
                          PRC_next(PRC_S), tag+1, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From N
       if ( PRC_HAS_N ) then
          call MPI_IRECV( var(:,:,JE+1:JA), ginfo(gid)%size2D_NS4*KA, COMM_datatype,              &
                          PRC_next(PRC_N), tag+2, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       if ( .not. PRC_TwoD ) then
          ! From E
          if ( PRC_HAS_E ) then
             call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                             PRC_next(PRC_E), tag+3, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                             PRC_next(PRC_W), tag+4, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
       end if

       !--- To 4-Direction HALO communicate
       if ( .not. PRC_TwoD ) then

          call pack_3D(var, gid, vid)

          ! To W HALO
          if ( PRC_HAS_W ) then
             call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                             PRC_next(PRC_W), tag+3, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype, &
                             PRC_next(PRC_E), tag+4, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
       end if
       ! To N HALO
       if ( PRC_HAS_N ) then
          call MPI_ISEND( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype,        &
                          PRC_next(PRC_N), tag+1, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          call MPI_ISEND( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype,        &
                          PRC_next(PRC_S), tag+2, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

    endif

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars_3D_mpi

  subroutine vars8_3D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: ireq, tag, tagc

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JA, JS, JE, JHALO

    integer :: ierr
    integer :: j
    !---------------------------------------------------------------------------

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    tag  = vid * 100
    ireq = 1

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(vid) ) then
       LOG_ERROR("vars8_3D_mpi",*) 'packing buffer is already used', vid
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          tagc = 0
          do j = 1, JS-1
             call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From SW
          tagc = 10
          do j = 1, JS-1
             call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NE
          tagc = 20
          do j = JE+1, JA
             call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NW
          tagc = 30
          do j = JE+1, JA
             call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From E
          tagc = 60
          call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                          PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! From W
          tagc = 70
          call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                          PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       end if
       ! From S
       tagc = 40
       do j = 1, JS-1
          call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                          PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JA
          call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                          PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !--- To 8-Direction HALO communicate
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                          PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                          PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       if ( .not. PRC_TwoD ) then

          call pack_3D(var, gid, vid)

          ! To W HALO
          tagc = 60
          call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                          PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! To E HALO
          tagc = 70
          call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                          PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! To NW HALO
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To NE HALO
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SW HALO
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SE HALO
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       end if

    else ! non-periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                     &
                                PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                     &
                                PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From SW
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NE
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                      &
                                PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                     &
                                PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                     &
                                PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NW
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                         &
                                PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From E
          if ( PRC_HAS_E ) then
             tagc = 60
             call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                             PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             tagc = 70
             call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                             PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
       end if
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = 1, JS-1
             call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                             PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JA
             call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                             PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                             PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype,                      &
                             PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       if ( .not. PRC_TwoD ) then

          call pack_3D(var, gid, vid)

          ! To W HALO
          if ( PRC_HAS_W ) then
             tagc = 60
             call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                             PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             tagc = 70
             call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype,    &
                             PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif

          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_NW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                       &
                                PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                &
                                PRC_next(PRC_NE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                     &
                                PRC_next(PRC_N), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,               &
                                PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_SW), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                        &
                                PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                       &
                                PRC_next(PRC_W), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                &
                                PRC_next(PRC_SE), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,                     &
                                PRC_next(PRC_S), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype,               &
                                PRC_next(PRC_E), tag+tagc, COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
       end if

    endif

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_3D_mpi

  subroutine vars_2D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: JA, JS, JE, JHALO

    integer :: ireq, tag
    integer :: ierr
    !---------------------------------------------------------------------------

    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(vid) ) then
       LOG_ERROR("vars_2D_mpi",*) 'packing buffer is already used', vid
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 4-Direction HALO communicate
        ! From S
        call MPI_IRECV( var(:,1:JS-1), ginfo(gid)%size2D_NS4,           &
                        COMM_datatype, PRC_next(PRC_S), tag+1,          &
                        COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
        ireq = ireq + 1

        ! From N
        call MPI_IRECV( var(:,JE+1:JA), ginfo(gid)%size2D_NS4,          &
                        COMM_datatype, PRC_next(PRC_N), tag+2,          &
                        COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
        ireq = ireq + 1

        if ( .not. PRC_TwoD ) then
           ! From E
           call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_E), tag+3,                &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1

           ! From W
           call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_W), tag+4,                &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1
        end if


        !--- To 4-Direction HALO communicate

        if ( .not. PRC_TwoD ) then

           call pack_2D(var, gid, vid)

           ! To W HALO communicate
           call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_W), tag+3,                &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1

           ! To E HALO communicate
           call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_E), tag+4,                &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1
        end if

        ! To N HALO communicate
        call MPI_ISEND( var(:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4,    &
                        COMM_datatype, PRC_next(PRC_N), tag+1,          &
                        COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
        ireq = ireq + 1

        ! To S HALO communicate
        call MPI_ISEND( var(:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4,    &
                        COMM_datatype, PRC_next(PRC_S), tag+2,          &
                        COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
        ireq = ireq + 1

    else
    !--- non-periodic condition
        !--- From 4-Direction HALO communicate
        ! From S
        if ( PRC_HAS_S ) then
            call MPI_IRECV( var(:,1:JS-1), ginfo(gid)%size2D_NS4,           &
                            COMM_datatype, PRC_next(PRC_S), tag+1,          &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! From N
        if ( PRC_HAS_N ) then
            call MPI_IRECV( var(:,JE+1:JA), ginfo(gid)%size2D_NS4,          &
                            COMM_datatype, PRC_next(PRC_N), tag+2,          &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        if ( .not. PRC_TwoD ) then
           ! From E
           if ( PRC_HAS_E ) then
              call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_E), tag+3,                &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif

           ! From W
           if ( PRC_HAS_W ) then
              call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_W), tag+4,                &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif
        end if


        !--- To 4-Direction HALO communicate

        if ( .not. PRC_TwoD ) then

           call pack_2D(var, gid, vid)

           ! To W HALO communicate
           if ( PRC_HAS_W ) then
              call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_W), tag+3,                &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif

           ! To E HALO communicate
           if ( PRC_HAS_E ) then
              call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_E), tag+4,                &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif
        end if

        ! To N HALO communicate
        if ( PRC_HAS_N ) then
            call MPI_ISEND( var(:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4,    &
                            COMM_datatype, PRC_next(PRC_N), tag+1,          &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! To S HALO communicate
        if ( PRC_HAS_S ) then
            call MPI_ISEND( var(:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4,    &
                            COMM_datatype, PRC_next(PRC_S), tag+2,          &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

    endif

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars_2D_mpi

  subroutine vars8_2D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: IS, IE, IHALO
    integer :: JA, JS, JE, JHALO

    integer :: ireq, tag, tagc

    integer :: ierr
    integer :: j
    !---------------------------------------------------------------------------

    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(vid) ) then
       LOG_ERROR("vars8_2D_mpi",*) 'packing buffer is already used', vid
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 8-Direction HALO communicate
        if ( .not. PRC_TwoD ) then
           ! From SE
           tagc = 0
           do j = 1, JS-1
              call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                              COMM_datatype, PRC_next(PRC_SE), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From SW
           tagc = 10
           do j = 1, JS-1
              call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                              COMM_datatype, PRC_next(PRC_SW), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From NE
           tagc = 20
           do j = JE+1, JA
              call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                              COMM_datatype, PRC_next(PRC_NE), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From NW
           tagc = 30
           do j = JE+1, JA
              call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                              COMM_datatype, PRC_next(PRC_NW), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From E
           call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_E), tag+60,               &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1
           ! From W
           call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_W), tag+70,               &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1
        end if
        ! From S
        tagc = 40
        do j = 1, JS-1
            call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,               &
                            COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
        enddo
        ! From N
        tagc = 50
        do j = JE+1, JA
            call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,               &
                            COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo


        !--- To 8-Direction HALO communicate

        ! To N HALO communicate
        tagc = 40
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,               &
                            COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To S HALO communicate
        tagc = 50
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,               &
                            COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                            COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        if ( .not. PRC_TwoD ) then

           call pack_2D(var, gid, vid)

           ! To W HALO communicate
           call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_W), tag+60,               &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1

           ! To E HALO communicate
           call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE, &
                           COMM_datatype, PRC_next(PRC_E), tag+70,               &
                           COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
           ireq = ireq + 1

           ! To NW HALO communicate
           tagc = 0
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                              COMM_datatype, PRC_next(PRC_NW), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To NE HALO communicate
           tagc = 10
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                              COMM_datatype, PRC_next(PRC_NE), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To SW HALO communicate
           tagc = 20
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                              COMM_datatype, PRC_next(PRC_SW), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To SE HALO communicate
           tagc = 30
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                              COMM_datatype, PRC_next(PRC_SE), tag+tagc,      &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        end if
    else
    !--- non-periodic condition
        !--- From 8-Direction HALO communicate
        if ( .not. PRC_TwoD ) then
           ! From SE
           if ( PRC_HAS_S .AND. PRC_HAS_E ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                                 COMM_datatype, PRC_next(PRC_SE), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                                 COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                                 COMM_datatype, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From SW
           if ( PRC_HAS_S .AND. PRC_HAS_W ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                                 COMM_datatype, PRC_next(PRC_SW), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                                 COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                                 COMM_datatype, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From NE
           if ( PRC_HAS_N .AND. PRC_HAS_E ) then
              tagc = 20
              do j = JE+1, JE+JHALO
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                                 COMM_datatype, PRC_next(PRC_NE), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                                 COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,              &
                                 COMM_datatype, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From NW
           if ( PRC_HAS_N .AND. PRC_HAS_W ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                                 COMM_datatype, PRC_next(PRC_NW), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                                 COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                 &
                                 COMM_datatype, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From E
           if ( PRC_HAS_E ) then
              call MPI_IRECV( ginfo(gid)%recvpack_E2P(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_E), tag+60,               &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif

           ! From W
           if ( PRC_HAS_W ) then
              call MPI_IRECV( ginfo(gid)%recvpack_W2P(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_W), tag+70,               &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif

        end if

        ! From S
        if ( PRC_HAS_S ) then
           tagc = 40
           do j = 1, JS-1
              call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,               &
                              COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        ! From N
        if ( PRC_HAS_N ) then
           tagc = 50
           do j = JE+1, JA
              call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,               &
                              COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif


        !! RECEIVE

        ! To N HALO communicate
        if ( PRC_HAS_N ) then
           tagc = 40
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,               &
                              COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        ! To S HALO communicate
        if ( PRC_HAS_S ) then
           tagc = 50
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,               &
                              COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        if ( .not. PRC_TwoD ) then

           call pack_2D(var, gid, vid)

           ! To W HALO communicate
           if ( PRC_HAS_W ) then
              call MPI_ISEND( ginfo(gid)%sendpack_P2W(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_W), tag+60,               &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif

           ! To E HALO communicate
           if ( PRC_HAS_E ) then
              call MPI_ISEND( ginfo(gid)%sendpack_P2E(:,vid), ginfo(gid)%size2D_WE, &
                              COMM_datatype, PRC_next(PRC_E), tag+70,               &
                              COMM_world, ginfo(gid)%req_list(ireq,vid), ierr       )
              ireq = ireq + 1
           endif

           ! To NW HALO communicate
           if ( PRC_HAS_N .AND. PRC_HAS_W ) then
              tagc = 0
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype, PRC_next(PRC_NW), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 10
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 20
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To NE HALO communicate
           if ( PRC_HAS_N .AND. PRC_HAS_E ) then
              tagc = 10
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                                 COMM_datatype, PRC_next(PRC_NE), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 0
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                                 COMM_datatype, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 30
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                                 COMM_datatype, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To SW HALO communicate
           if ( PRC_HAS_S .AND. PRC_HAS_W ) then
              tagc = 20
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype, PRC_next(PRC_SW), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 30
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 0
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To SE HALO communicate
           if ( PRC_HAS_S .AND. PRC_HAS_E ) then
              tagc = 30
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                                 COMM_datatype, PRC_next(PRC_SE), tag+tagc,      &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 20
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                                 COMM_datatype, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 10
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,        &
                                 COMM_datatype, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif
        end if

    endif

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_2D_mpi

  subroutine vars_3D_mpi_pc(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid
    integer :: ierr
    !---------------------------------------------------------------------------

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(ginfo(gid)%pseqid(vid)) ) then
       LOG_ERROR("vars_3D_mpi_pc",*) 'packing buffer is already used', vid, ginfo(gid)%pseqid(vid)
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(ginfo(gid)%pseqid(vid)) = .true.
#endif

    if ( .not. PRC_TwoD ) call pack_3D(var, gid, ginfo(gid)%pseqid(vid))

    call MPI_STARTALL(ginfo(gid)%preq_cnt(vid), ginfo(gid)%preq_list(1:ginfo(gid)%preq_cnt(vid),vid), ierr)

    return
  end subroutine vars_3D_mpi_pc

  subroutine wait_3D_mpi(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call MPI_WAITALL( ginfo(gid)%req_cnt (vid),                           &
                      ginfo(gid)%req_list(1:ginfo(gid)%req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                &
                      ierr                                                )
    if ( .not. PRC_TwoD ) call unpack_3D(var, gid, vid)

#ifdef DEBUG
    ginfo(gid)%use_packbuf(vid) = .false.
#endif

    return
  end subroutine wait_3D_mpi

  subroutine wait_2D_mpi(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call MPI_WAITALL( ginfo(gid)%req_cnt(vid),                            &
                      ginfo(gid)%req_list(1:ginfo(gid)%req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                &
                      ierr                                                )
    if ( .not. PRC_TwoD ) call unpack_2D(var, gid, vid)

#ifdef DEBUG
    ginfo(gid)%use_packbuf(vid) = .false.
#endif

    return
  end subroutine wait_2D_mpi

  subroutine wait_3D_mpi_pc(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: ierr

    !--- wait packets
    call MPI_WAITALL( ginfo(gid)%preq_cnt (vid),                            &
                      ginfo(gid)%preq_list(1:ginfo(gid)%preq_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                  &
                      ierr                                                  )
    if ( .not. PRC_TwoD ) call unpack_3D(var, gid, ginfo(gid)%pseqid(vid))

#ifdef DEBUG
    ginfo(gid)%use_packbuf(ginfo(gid)%pseqid(vid)) = .false.
#endif

    return
  end subroutine wait_3D_mpi_pc

  subroutine pack_3D( var, gid, vid )
    implicit none

    real(RP), intent(in) :: var(:,:,:)
    integer,  intent(in) :: gid
    integer,  intent(in) :: vid

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JS, JE

    integer :: k, i, j, n

    call PROF_rapstart('COMM_pack', 3)

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- packing packets to West
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IS+IHALO-1
       do k = 1, KA
          n = (j-JS) * KA * IHALO &
            + (i-IS) * KA         &
            + k
          ginfo(gid)%sendpack_P2W(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo
       !--- packing packets to East
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE-IHALO+1, IE
       do k = 1, KA
          n = (j-JS)         * KA * IHALO &
            + (i-IE+IHALO-1) * KA         &
            + k
          ginfo(gid)%sendpack_P2E(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo

    else

       if ( PRC_HAS_W ) then
          !--- packing packets to West
          !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IHALO,KA,var,ginfo,gid,vid)
          do j = JS, JE
          do i = IS, IS+IHALO-1
          do k = 1, KA
             n = (j-JS) * KA * IHALO &
               + (i-IS) * KA         &
               + k
             ginfo(gid)%sendpack_P2W(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif
       if ( PRC_HAS_E ) then
          !--- packing packets to East
          !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IE,IHALO,KA,var,ginfo,gid,vid)
          do j = JS, JE
          do i = IE-IHALO+1, IE
          do k = 1, KA
             n = (j-JS)         * KA * IHALO &
               + (i-IE+IHALO-1) * KA         &
               + k
             ginfo(gid)%sendpack_P2E(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif

    end if

    call PROF_rapend('COMM_pack', 3)

    return
  end subroutine pack_3D

  subroutine pack_2D(var, gid, vid)
    implicit none

    real(RP), intent(in) :: var(:,:)
    integer,  intent(in) :: vid
    integer,  intent(in) :: gid

    integer :: IS, IE, IHALO
    integer :: JS, JE

    integer :: i, j, n

    call PROF_rapstart('COMM_pack', 3)

    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- To 4-Direction HALO communicate
       !--- packing packets to West
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IS+IHALO-1
          n = (j-JS) * IHALO &
            + (i-IS) + 1
          ginfo(gid)%sendpack_P2W(n,vid) = var(i,j)
       enddo
       enddo

       !--- packing packets to East
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE-IHALO+1, IE
          n = (j-JS)         * IHALO &
            + (i-IE+IHALO-1) + 1
          ginfo(gid)%sendpack_P2E(n,vid) = var(i,j)
       enddo
       enddo

    else

       !--- To 4-Direction HALO communicate
       !--- packing packets to West
       if ( PRC_HAS_W ) then
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IS+IHALO-1
             n = (j-JS) * IHALO &
               + (i-IS) + 1
             ginfo(gid)%sendpack_P2W(n,vid) = var(i,j)
          enddo
          enddo
       endif

       !--- packing packets to East
       if ( PRC_HAS_E ) then
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IE-IHALO+1, IE
             n = (j-JS)         * IHALO &
               + (i-IE+IHALO-1) + 1
             ginfo(gid)%sendpack_P2E(n,vid) = var(i,j)
          enddo
          enddo
       endif

    end if

    call PROF_rapend('COMM_pack', 3)

    return
  end subroutine pack_2D

  subroutine unpack_3D(var, gid, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JS, JE
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_unpack', 3)

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- unpacking packets from East
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE+1, IE+IHALO
       do k = 1, KA
          n = (j-JS)   * KA * IHALO &
            + (i-IE-1) * KA         &
            + k
          var(k,i,j) = ginfo(gid)%recvpack_E2P(n,vid)
       enddo
       enddo
       enddo

       !--- unpacking packets from West
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS-IHALO, IS-1
       do k = 1, KA
          n = (j-JS)       * KA * IHALO &
            + (i-IS+IHALO) * KA         &
            + k
          var(k,i,j) = ginfo(gid)%recvpack_W2P(n,vid)
       enddo
       enddo
       enddo

    else ! non-periodic condition

        if ( PRC_HAS_E ) then
           !--- unpacking packets from East
           !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
           !$omp shared(JS,JE,IE,IHALO,KA,var,ginfo,gid,vid)
           do j = JS, JE
           do i = IE+1, IE+IHALO
           do k = 1, KA
              n = (j-JS)   * KA * IHALO &
                + (i-IE-1) * KA         &
                + k
              var(k,i,j) = ginfo(gid)%recvpack_E2P(n,vid)
           enddo
           enddo
           enddo
        endif

        if ( PRC_HAS_W ) then
           !--- unpacking packets from West
           !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
           !$omp shared(JS,JE,IS,IHALO,KA,var,ginfo,gid,vid)
           do j = JS, JE
           do i = IS-IHALO, IS-1
           do k = 1, KA
              n = (j-JS)       * KA * IHALO &
                + (i-IS+IHALO) * KA         &
                + k
              var(k,i,j) = ginfo(gid)%recvpack_W2P(n,vid)
           enddo
           enddo
           enddo
        endif

    end if

    call PROF_rapend('COMM_unpack', 3)

    return
  end subroutine unpack_3D

  subroutine unpack_2D(var, gid, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: IS, IE, IHALO
    integer :: JS, JE

    integer :: i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_unpack', 3)

    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE

    if( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- unpacking packets from East
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IE+1, IE+IHALO
           n = (j-JS)   * IHALO &
             + (i-IE-1) + 1
           var(i,j) = ginfo(gid)%recvpack_E2P(n,vid)
        enddo
        enddo

        !--- unpacking packets from West
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IS-IHALO, IS-1
           n = (j-JS)       * IHALO &
             + (i-IS+IHALO) + 1
           var(i,j) = ginfo(gid)%recvpack_W2P(n,vid)
        enddo
        enddo

    else
    !--- non-periodic condition

       !--- unpacking packets from East / copy inner data to HALO(East)
       if( PRC_HAS_E ) then
          !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IE+1, IE+IHALO
             n = (j-JS)   * IHALO &
               + (i-IE-1) + 1
             var(i,j) = ginfo(gid)%recvpack_E2P(n,vid)
          enddo
          enddo
       end if


       !--- unpacking packets from West / copy inner data to HALO(West)
       if( PRC_HAS_W ) then
          !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS-IHALO, IS-1
             n = (j-JS)       * IHALO &
               + (i-IS+IHALO) + 1
             var(i,j) = ginfo(gid)%recvpack_W2P(n,vid)
          enddo
          enddo
       end if

    end if

    call PROF_rapend('COMM_unpack', 3)

    return
  end subroutine unpack_2D

  subroutine copy_boundary_3D(var, gid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid

    integer :: KA
    integer :: IS, IE, IHALO
    integer :: JS, JE, JHALO

    integer :: k, i, j
    !---------------------------------------------------------------------------

    KA    = ginfo(gid)%KA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$omp parallel

    !--- copy inner data to HALO(North)
    if ( .NOT. PRC_HAS_N ) then
       do j = JE+1, JE+JHALO
       !$omp do
       do i = IS, IE
       do k = 1, KA
          var(k,i,j) = var(k,i,JE)
       enddo
       enddo
       !$omp end do nowait
       enddo
    endif

    !--- copy inner data to HALO(South)
    if ( .NOT. PRC_HAS_S ) then
       do j = JS-JHALO, JS-1
       !$omp do
       do i = IS, IE
       do k = 1, KA
          var(k,i,j) = var(k,i,JS)
       enddo
       enddo
       !$omp end do nowait
       enddo
    endif

    if ( .not. PRC_TwoD ) then

       !--- copy inner data to HALO(East)
       if ( .NOT. PRC_HAS_E ) then
          !$omp do
          do j = JS, JE
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$omp end do nowait
       end if

       !--- copy inner data to HALO(West)
       if ( .NOT. PRC_HAS_W ) then
          !$omp do
          do j = JS, JE
          do i = IS-IHALO, IS-1
             var(:,i,j) = var(:,IS,j)
          enddo
          enddo
          !$omp end do nowait
       end if

       !--- copy inner data to HALO(NorthWest)
       if ( .NOT. PRC_HAS_N .AND. &
            .NOT. PRC_HAS_W ) then
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,JE)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_N ) then
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,i,JE)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_W ) then
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,j)
          enddo
          enddo
          enddo
       endif

       !--- copy inner data to HALO(SouthWest)
       if ( .NOT. PRC_HAS_S .AND. &
            .NOT. PRC_HAS_W ) then
          do j = JS-IHALO, JS-1
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,JS)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_S ) then
          do j = JS-IHALO, JS-1
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,i,JS)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_W ) then
          do j = JS-IHALO, JS-1
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,j)
          enddo
          enddo
          enddo
       endif

       !--- copy inner data to HALO(NorthEast)
       if ( .NOT. PRC_HAS_N .AND. &
            .NOT. PRC_HAS_E ) then
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,JE)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_N ) then
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,i,JE)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_E ) then
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
       endif

       !--- copy inner data to HALO(SouthEast)
       if ( .NOT. PRC_HAS_S .AND. &
            .NOT. PRC_HAS_E ) then
          do j = JS-IHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,JS)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_S ) then
          do j = JS-IHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,i,JS)
          enddo
          enddo
          enddo
       elseif( .NOT. PRC_HAS_E ) then
          do j = JS-IHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
       endif

    end if

    !$omp end parallel

    return
  end subroutine copy_boundary_3D

  subroutine copy_boundary_2D(var, gid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid

    integer :: IS, IE, IHALO
    integer :: JS, JE, JHALO

    integer :: i, j
    !---------------------------------------------------------------------------

    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$omp parallel

    !--- copy inner data to HALO(North)
    if( .NOT. PRC_HAS_N ) then
       do j = JE+1, JE+JHALO
       !$omp do
       do i = IS, IE
          var(i,j) = var(i,JE)
       enddo
       !$omp end do nowait
       enddo
    endif

    !--- copy inner data to HALO(South)
    if( .NOT. PRC_HAS_S ) then
       do j = JS-JHALO, JS-1
       !$omp do
       do i = IS, IE
          var(i,j) = var(i,JS)
       enddo
       !$omp end do nowait
       enddo
    endif

    if ( .not. PRC_TwoD ) then

       if( .NOT. PRC_HAS_E ) then
          !$omp do
          do j = JS, JE
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$omp end do nowait
       endif

       if( .NOT. PRC_HAS_W ) then
          !$omp do
          do j = JS, JE
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$omp end do nowait
       endif

       !--- copy inner data to HALO(NorthWest)
       if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_W ) then
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,JE)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_N ) then
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
             var(i,j) = var(i,JE)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_W ) then
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
       endif

       !--- copy inner data to HALO(SouthWest)
       if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_W ) then
          do j = JS-IHALO, JS-1
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,JS)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_S ) then
          do j = JS-IHALO, JS-1
          do i = IS-IHALO, IS-1
             var(i,j) = var(i,JS)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_W ) then
          do j = JS-IHALO, JS-1
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
       endif

       !--- copy inner data to HALO(NorthEast)
       if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_E ) then
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,JE)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_N ) then
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(i,JE)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_E ) then
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
       endif

       !--- copy inner data to HALO(SouthEast)
       if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_E ) then
          do j = JS-IHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,JS)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_S ) then
          do j = JS-IHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(i,JS)
          enddo
          enddo
       elseif( .NOT. PRC_HAS_E ) then
          do j = JS-IHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
       endif

    end if

    !$omp end parallel

    return
  end subroutine copy_boundary_2D

end module scale_comm_cartesC
