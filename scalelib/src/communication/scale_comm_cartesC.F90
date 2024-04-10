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
  use mpi_f08
  use iso_c_binding
  use scale_precision
  use scale_io
  use scale_prof
  use scale_tracer
#ifdef _OPENACC
  use openacc
#endif

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
     module procedure COMM_bcast_SCR_SP
     module procedure COMM_bcast_SCR_DP
     module procedure COMM_bcast_1D_SP
     module procedure COMM_bcast_1D_DP
     module procedure COMM_bcast_2D_SP
     module procedure COMM_bcast_2D_DP
     module procedure COMM_bcast_3D_SP
     module procedure COMM_bcast_3D_DP
     module procedure COMM_bcast_4D_SP
     module procedure COMM_bcast_4D_DP
     module procedure COMM_bcast_INT_SCR
     module procedure COMM_bcast_INT_1D
     module procedure COMM_bcast_INT_2D
     module procedure COMM_bcast_LOGICAL_SCR
     module procedure COMM_bcast_LOGICAL_1D
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
  integer,  private              :: COMM_vsize_max_pc   !< # limit of total communication variables for MPI PC

  logical,  private              :: COMM_IsAllPeriodic  !< periodic boundary condition?

  logical,  private              :: COMM_USE_MPI_PC       = .true.  !< MPI persistent communication
#ifdef __FUJITSU
  logical,  private              :: COMM_USE_MPI_ONESIDED = .true.  !< MPI one-sided communication
#else
  logical,  private              :: COMM_USE_MPI_ONESIDED = .false. !< MPI one-sided communication
#endif

  type(MPI_Datatype), public :: COMM_datatype_t
  type(MPI_Comm),     public :: COMM_world_t


#ifdef _OPENACC
  type ptr_t
     real(RP), pointer :: ptr(:,:,:)
  end type ptr_t
#endif
  type ginfo_t
     integer              :: KA
     integer              :: IA, IS, IE, IHALO
     integer              :: JA, JS, JE, JHALO
     integer              :: nreq_max             !< # limit of communication request at once
     integer              :: size2D_NS4           !< 2D data size (W/E    HALO, 4/8-direction comm.)
     integer              :: size2D_NS8           !< 2D data size (N/S    HALO,   4-direction comm.)
     integer              :: size2D_WE            !< 2D data size (N/S    HALO,   8-direction comm.)
     integer              :: size2D_4C            !< 2D data size (corner HALO,   8-direction comm.)
     integer              :: vars_num = 0         !< numbers of variables for persistent comm.
     real(RP),    pointer :: recvpack_WE2P(:,:,:) !< packing packet (receive, from W and E)
     real(RP),    pointer :: sendpack_P2WE(:,:,:) !< packing packet (send,    to   W and E)
     type(c_ptr), allocatable :: recvbuf_WE(:)    !< receive buffer for MPI_Put (from W and E)
     type(c_ptr), allocatable :: recvbuf_NS(:)    !< receive buffer for MPI_Put (from N and S)
     integer,           allocatable :: req_cnt (:)    !< request ID of each MPI send/recv
     type(MPI_Request), allocatable :: req_list(:,:)  !< request ID set of each variables
     integer,           allocatable :: preq_cnt (:)   !< request ID of each MPI PC
     type(MPI_Request), allocatable :: preq_list(:,:) !< request ID set of each variables for MPI PC
     integer,       allocatable :: packid(:)          !< ID of pack
     type(MPI_Win), allocatable :: win_packWE(:)      !< window ID for MPI onesided
     type(MPI_Win), allocatable :: win_packNS(:)      !< window ID for MPI onesided
#ifdef DEBUG
     logical,       allocatable :: use_packbuf(:)     !< using flag for packing buffer
#endif
#ifdef _OPENACC
     logical,     allocatable :: device_alloc(:)
     type(ptr_t), allocatable :: device_ptr(:)
#endif
  end type ginfo_t

  integer, private, parameter  :: COMM_gid_max = 20
  integer, private             :: COMM_gid
  type(ginfo_t), private       :: ginfo(COMM_gid_max)

  type(MPI_Group), private     :: group_packWE !< MPI_Group for pack
  type(MPI_Group), private     :: group_packNS !< MPI_Group for vars
  logical, private             :: group_packWE_created = .false.
  logical, private             :: group_packNS_created = .false.

  logical, private             :: initialized = .false.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine COMM_setup
    use scale_prc, only: &
       PRC_abort, &
       PRC_LOCAL_COMM_WORLD
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    namelist / PARAM_COMM_CARTESC / &
       COMM_vsize_max, &
       COMM_vsize_max_pc, &
       COMM_USE_MPI_PC, &
       COMM_USE_MPI_ONESIDED

    integer         :: ranks(8)
    type(MPI_Group) :: group

    integer :: n, m
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("COMM_setup",*) 'Setup'

    COMM_vsize_max = max( 10 + QA*2, 25 )
    COMM_vsize_max_pc = 50 + QA*2

#ifdef _OPENACC
    COMM_USE_MPI_ONESIDED = .false.
#endif

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
       COMM_datatype_t = MPI_DOUBLE_PRECISION
    elseif( RP == kind(0.0) ) then
       COMM_datatype_t = MPI_REAL
    else
       LOG_ERROR("COMM_setup",*) 'precision is not supportd'
       call PRC_abort
    endif
    COMM_datatype = COMM_datatype_t%MPI_VAL

    COMM_world = PRC_LOCAL_COMM_WORLD
    COMM_world_t%MPI_VAL = COMM_world

    COMM_gid = 0

#ifdef _OPENACC
    if ( COMM_USE_MPI_ONESIDED ) then
       LOG_WARN("COMM_setup",*) "Open MPI does not support one-sided APIs with CUDA-aware UCX"
    end if
#endif

    if ( COMM_USE_MPI_ONESIDED ) then

       COMM_USE_MPI_PC = .false.

       call MPI_Comm_group( COMM_world_t, group, ierr )

       n = 0
       if ( PRC_HAS_S ) then
          n = 1
          ranks(n) = PRC_next(PRC_S)
       end if
       if ( PRC_HAS_N ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_N) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_N)
          end if
       end if
       if ( PRC_HAS_N .and. PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_NW) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_NW)
          end if
       else if ( PRC_HAS_N ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_N) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_N)
          end if
       else if ( PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_W) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_W)
          end if
       end if
       if ( PRC_HAS_N .and. PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_NE) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_NE)
          end if
       else if ( PRC_HAS_N ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_N) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_N)
          end if
       else if ( PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_E) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_E)
          end if
       end if
       if ( PRC_HAS_S .and. PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_SW) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_SW)
          end if
       else if ( PRC_HAS_S ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_S) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_S)
          end if
       else if ( PRC_HAS_W ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_W) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_W)
          end if
       end if
       if ( PRC_HAS_S .and. PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_SE) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_SE)
          end if
       else if ( PRC_HAS_S ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_S) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_S)
          end if
       else if ( PRC_HAS_E ) then
          do m = 1, n
             if ( ranks(m) == PRC_next(PRC_E) ) exit
          end do
          if ( m == n + 1 ) then
             n = n + 1
             ranks(n) = PRC_next(PRC_E)
          end if
       end if
       if ( n > 0 ) then
          call MPI_Group_incl( group, n, ranks, group_packNS, ierr )
          group_packNS_created = .true.
       else
          group_packNS_created = .false.
       end if

       n = 0
       if ( .not. PRC_TwoD ) then
          if ( PRC_HAS_W ) then
             n = 1
             ranks(n) = PRC_next(PRC_W)
          end if
          if ( PRC_HAS_E ) then
             if ( n == 0 .or. ranks(1) .ne. PRC_next(PRC_E) ) then
                n = n + 1
                ranks(n) = PRC_next(PRC_E)
             end if
          end if
       end if
       if ( n > 0 ) then
          call MPI_Group_incl( group, n, ranks, group_packWE, ierr )
          group_packWE_created = .true.
       else
          group_packWE_created = .false.
       end if

       call MPI_Group_free( group, ierr )
    end if

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
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in)  :: KA, IA, JA, IHALO, JHALO
    integer, intent(out) :: gid

    integer :: IMAX, JMAX
    integer :: nreq_NS, nreq_WE, nreq_4C

    type(MPI_Info) :: win_info
    integer(kind=MPI_ADDRESS_KIND) :: size

    integer :: ierr
    integer :: n

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

    allocate( ginfo(gid)%sendpack_P2WE(ginfo(gid)%size2D_WE * KA, 2, COMM_vsize_max) )
    !$acc enter data create(ginfo(gid)%sendpack_P2WE)

#ifdef DEBUG
    allocate( ginfo(gid)%use_packbuf(COMM_vsize_max) )
    ginfo(gid)%use_packbuf(:) = .false.
#endif

#ifdef _OPENACC
    allocate( ginfo(gid)%device_alloc(COMM_vsize_max+COMM_vsize_max_pc) )
    allocate( ginfo(gid)%device_ptr(COMM_vsize_max+1:COMM_vsize_max_pc) )
    ginfo(gid)%device_alloc(:) = .false.
#endif

    if ( COMM_USE_MPI_ONESIDED ) then

       allocate( ginfo(gid)%recvbuf_WE(COMM_vsize_max) )
       allocate( ginfo(gid)%recvbuf_NS(COMM_vsize_max) )

       allocate( ginfo(gid)%win_packWE(COMM_vsize_max) )
       allocate( ginfo(gid)%win_packNS(COMM_vsize_max) )

       call MPI_Info_create(win_info, ierr)
       call MPI_Info_set(win_info, "no_locks", "true", ierr)
       call MPI_Info_set(win_info, "same_size", "true", ierr)
       call MPI_Info_set(win_info, "same_disp_unit", "true", ierr)

       do n = 1, COMM_vsize_max
          size = ginfo(gid)%size2D_WE * KA * 2 * RP
#ifdef _OPENACC
          block
            real(RP), pointer :: pack(:)
            call MPI_Alloc_mem(size, MPI_INFO_NULL, ginfo(gid)%recvbuf_WE(n), ierr)
            call C_F_pointer(ginfo(gid)%recvbuf_WE(n), pack, (/ size/RP /))
            !$acc enter data create(pack)
            !$acc host_data use_device(pack)
            call MPI_Win_create(pack, size, ginfo(gid)%size2D_WE*KA*RP, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%win_packWE(n), ierr)
            !$acc end host_data
          end block
#else
          call MPI_Win_allocate(size, ginfo(gid)%size2D_WE*KA*RP, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%recvbuf_WE(n), ginfo(gid)%win_packWE(n), ierr)
#endif
          size = ginfo(gid)%size2D_NS4 * KA * 2 * RP
#ifdef _OPENACC
          block
            real(RP), pointer :: pack(:)
            call MPI_Alloc_mem(size, MPI_INFO_NULL, ginfo(gid)%recvbuf_NS(n), ierr)
            call C_F_pointer(ginfo(gid)%recvbuf_NS(n), pack, (/ size/RP /))
            !$acc enter data create(pack)
            !$acc host_data use_device(pack)
            call MPI_Win_create(pack, size, RP, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%win_packNS(n), ierr)
            !$acc end host_data
          end block
#else
          call MPI_Win_allocate(size, RP, &
                                win_info, COMM_world_t, &
                                ginfo(gid)%recvbuf_NS(n), ginfo(gid)%win_packNS(n), ierr)
#endif
       end do

       call MPI_Info_free(win_info, ierr)

       do n = 1, COMM_vsize_max
          call MPI_Win_post( group_packWE, MPI_MODE_NOSTORE, ginfo(gid)%win_packWE(n), ierr )
          call MPI_Win_post( group_packNS, MPI_MODE_NOSTORE, ginfo(gid)%win_packNS(n), ierr )
       end do

       ginfo(gid)%vars_num = 0
       allocate( ginfo(gid)%packid(COMM_vsize_max_pc) )

    else

       allocate( ginfo(gid)%recvpack_WE2P(ginfo(gid)%size2D_WE * KA, 2, COMM_vsize_max) )
       !$acc enter data create(ginfo(gid)%recvpack_WE2P)

       allocate( ginfo(gid)%req_cnt (                     COMM_vsize_max) )
       allocate( ginfo(gid)%req_list(ginfo(gid)%nreq_MAX, COMM_vsize_max) )
       ginfo(gid)%req_cnt (:)   = -1
       ginfo(gid)%req_list(:,:) = MPI_REQUEST_NULL

       if ( COMM_USE_MPI_PC ) then
          ginfo(gid)%vars_num = 0
          allocate( ginfo(gid)%packid(COMM_vsize_max_pc) )
          allocate( ginfo(gid)%preq_cnt (                      COMM_vsize_max_pc) )
          allocate( ginfo(gid)%preq_list(ginfo(gid)%nreq_MAX+1,COMM_vsize_max_pc) )
          ginfo(gid)%preq_cnt (:)   = -1
          ginfo(gid)%preq_list(:,:) = MPI_REQUEST_NULL
       end if

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
    implicit none

    integer :: gid
    integer :: i, j, ierr
    !---------------------------------------------------------------------------

    do gid = 1, COMM_gid

       if ( COMM_USE_MPI_ONESIDED ) then

          do i = 1, COMM_vsize_max
             call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(i), ierr )
             call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(i), ierr )
         end do

          do i = 1, COMM_vsize_max
             call MPI_Win_complete( ginfo(gid)%win_packWE(i), ierr )
             call MPI_Win_complete( ginfo(gid)%win_packNS(i), ierr )
          end do

          do i = 1, COMM_vsize_max
             call MPI_Win_wait( ginfo(gid)%win_packWE(i), ierr )
             call MPI_Win_wait( ginfo(gid)%win_packNS(i), ierr )
          end do

          do i = 1, COMM_vsize_max
             call MPI_Win_free(ginfo(gid)%win_packWE(i), ierr)
             call MPI_Win_free(ginfo(gid)%win_packNS(i), ierr)
#ifdef _OPENACC
             block
               real(RP), pointer :: pack(:)
               integer :: KA
               KA = ginfo(gid)%KA
               call C_F_pointer( ginfo(gid)%recvbuf_WE(i), pack, (/ginfo(gid)%size2D_WE*KA*2/) )
               !$acc exit data delete(pack)
               call C_F_pointer( ginfo(gid)%recvbuf_NS(i), pack, (/ginfo(gid)%size2D_NS4*KA*2/) )
               !$acc exit data delete(pack)
             end block
             call MPI_Free_mem(ginfo(gid)%recvbuf_WE(i), ierr)
             call MPI_Free_mem(ginfo(gid)%recvbuf_NS(i), ierr)
#endif
          end do

          deallocate( ginfo(gid)%packid )
          ginfo(gid)%vars_num = 0

          deallocate( ginfo(gid)%win_packWE )
          deallocate( ginfo(gid)%win_packNS )

          deallocate( ginfo(gid)%recvbuf_WE )
          deallocate( ginfo(gid)%recvbuf_NS )

       else

          if ( COMM_USE_MPI_PC ) then

             do j = 1, COMM_vsize_max_pc
                do i = 1, ginfo(gid)%nreq_MAX+1
                   if (ginfo(gid)%preq_list(i,j) .NE. MPI_REQUEST_NULL) &
                        call MPI_REQUEST_FREE(ginfo(gid)%preq_list(i,j), ierr)
                enddo
#ifdef _OPENACC
                if ( ginfo(gid)%device_alloc(j+COMM_vsize_max) ) then
                   !$acc exit data delete(ginfo(gid)%device_ptr(j+COMM_vsize_max)%ptr)
                end if
#endif
             enddo
             deallocate( ginfo(gid)%preq_cnt )
             deallocate( ginfo(gid)%preq_list )
             deallocate( ginfo(gid)%packid )
             ginfo(gid)%vars_num = 0

          end if

          deallocate( ginfo(gid)%req_cnt )
          deallocate( ginfo(gid)%req_list )

          !$acc exit data delete(ginfo(gid)%recvpack_WE2P)
          deallocate( ginfo(gid)%recvpack_WE2P )

       end if

       !$acc exit data delete(ginfo(gid)%sendpack_P2WE)
       deallocate( ginfo(gid)%sendpack_P2WE )
#ifdef DEBUG
       deallocate( ginfo(gid)%use_packbuf )
#endif

    end do

    if ( COMM_USE_MPI_ONESIDED ) then
       if ( group_packWE_created ) then
          call MPI_Group_free(group_packWE, ierr)
          group_packWE_created = .false.
       end if
       if ( group_packNS_created ) then
          call MPI_Group_free(group_packNS, ierr)
          group_packNS_created = .false.
       end if
    end if


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
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name
    real(RP), target, intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID

    integer,          intent(in), optional :: gid

    integer :: gid_
    integer :: vars_id
    !---------------------------------------------------------------------------

    if ( .not. COMM_USE_MPI_PC ) return
#ifdef _OPENACC
    if ( .not. acc_is_present(var) ) return
#endif

    call PROF_rapstart('COMM_init_pers', 2)

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

    ginfo(gid_)%vars_num = ginfo(gid_)%vars_num + 1
    if ( ginfo(gid_)%vars_num > COMM_vsize_max_pc ) then
       LOG_ERROR("COMM_vars_init",*) 'number of variable for MPI PC exceeds max', ginfo(gid_)%vars_num, COMM_vsize_max_pc
       call PRC_abort
    end if

    vars_id = ginfo(gid_)%vars_num
    ginfo(gid_)%packid(vars_id) = vid

#ifdef _OPENACC
    if ( .not. acc_is_present(var) ) then
       ginfo(gid_)%device_alloc(vars_id+COMM_vsize_max) = .true.
       ginfo(gid_)%device_ptr(vars_id*COMM_vsize_max)%ptr => var
       !$acc enter data copyin(var)
    end if
#endif

    call vars_init_mpi_pc(var, gid_, vars_id, vid)

    vid = vars_id + COMM_vsize_max

    LOG_INFO("COMM_vars_init",'(1x,A,I3.3,A,I3.3,2A)') 'Initialize variable (grid ID = ', gid_, '): ID = ', vid, &
                                                                                       ', name = ', trim(varname)

    call PROF_rapend  ('COMM_init_pers', 2)

    return
  end subroutine COMM_vars_init

  !-----------------------------------------------------------------------------
  !> Register variables
  subroutine COMM_vars8_init( &
       varname, &
       var,     &
       vid,     &
       gid      )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name

    real(RP), target, intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID

    integer,          intent(in), optional :: gid

    integer :: gid_
    integer :: vars_id
    !---------------------------------------------------------------------------

    if ( .not. COMM_USE_MPI_PC ) return
#ifdef _OPENACC
    if ( .not. acc_is_present(var) ) return
#endif

    call PROF_rapstart('COMM_init_pers', 2)

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

    ginfo(gid_)%vars_num = ginfo(gid_)%vars_num + 1
    if ( ginfo(gid_)%vars_num > COMM_vsize_max_pc ) then
       LOG_ERROR("COMM_vars8_init",*) 'number of variable for MPI PC exceeds max', ginfo(gid_)%vars_num, COMM_vsize_max_pc
       call PRC_abort
    end if

    vars_id = ginfo(gid_)%vars_num
    ginfo(gid_)%packid(vars_id) = vid

#ifdef _OPENACC
    if ( .not. acc_is_present(var) ) then
       ginfo(gid_)%device_alloc(vars_id+COMM_vsize_max) = .true.
       ginfo(gid_)%device_ptr(vars_id+COMM_vsize_max)%ptr => var
       !$acc enter data copyin(var)
    end if
#endif

    call vars8_init_mpi_pc(var, gid_, vars_id, vid)

    vid = vars_id + COMM_vsize_max

    LOG_INFO("COMM_vars8_init",'(1x,A,I3.3,A,I3.3,2A)') 'Initialize variable (grid ID = ', gid_, '): ID = ', vid, &
                                                                                       ', name = ', trim(varname)

    call PROF_rapend  ('COMM_init_pers', 2)

    return
  end subroutine COMM_vars8_init

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_3D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
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
       if ( COMM_USE_MPI_ONESIDED ) then
          call vars_3D_mpi_onesided(var, gid_, vid)
       else
          call vars_3D_mpi(var, gid_, vid)
       end if
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_3D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
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
       if ( COMM_USE_MPI_ONESIDED ) then
          call vars8_3D_mpi_onesided(var, gid_, vid)
       else
          call vars8_3D_mpi(var, gid_, vid)
       end if
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars8_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_3D(var, vid, FILL_BND, gid)
    use scale_prc, only: &
       PRC_abort
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
       if ( COMM_USE_MPI_ONESIDED ) then
          call wait_3D_mpi_onesided(var, gid_, vid)
       else
          call wait_3D_mpi(var, gid_, vid)
       end if
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
    use scale_prc, only: &
       PRC_abort
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
    if ( COMM_USE_MPI_ONESIDED ) then
       call vars_2D_mpi_onesided(var, gid_, vid)
    else
       call vars_2D_mpi(var, gid_, vid)
    end if
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_2D(var, vid, gid)
    use scale_prc, only: &
       PRC_abort
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
    if ( COMM_USE_MPI_ONESIDED ) then
       call vars8_2D_mpi_onesided(var, gid_, vid)
    else
       call vars8_2D_mpi(var, gid_, vid)
    end if
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars8_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_2D(var, vid, FILL_BND, gid)
    use scale_prc, only: &
       PRC_abort
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
    if ( COMM_USE_MPI_ONESIDED ) then
       call wait_2D_mpi_onesided(var, gid_, vid)
    else
       call wait_2D_mpi(var, gid_, vid)
    end if
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
    real(DP) :: stat1, stat2
    real(DP) :: allstat(2)
    real(DP) :: zerosw

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    stat1 = 0.0_DP
    stat2 = 0.0_DP
    !$omp parallel do reduction(+:stat1,stat2)
    !$acc kernels if(acc_is_present(var))
    !$acc loop reduction(+:stat1,stat2)
    do j = JS, JE
    !$acc loop reduction(+:stat1,stat2)
    do i = IS, IE
       if ( abs(var(i,j)) < abs(CONST_UNDEF) ) then
          stat1 = stat1 + var(i,j)
          stat2 = stat2 + 1.0_DP
       endif
    enddo
    enddo
    !$acc end kernels

    stat(:) = (/stat1, stat2/)

    ! All reduce
    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    call MPI_Allreduce( stat,                 &
                        allstat,              &
                        2,                    &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM,              &
                        COMM_world_t,         &
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
#ifdef _OPENACC
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    !$acc data create(stat, allstat) if(flag_device)

    !$acc kernels if(flag_device)
    stat(:,:) = 0.0_DP
    !$acc end kernels
    !$acc kernels if(flag_device)
    !$acc loop independent
    do j = JS, JE
    !$acc loop independent
    do i = IS, IE
    do k = 1,  KA
       if ( abs(var(k,i,j)) < abs(CONST_UNDEF) ) then
          !$acc atomic update
          stat(k,1) = stat(k,1) + var(k,i,j)
          !$acc end atomic
          !$acc atomic update
          stat(k,2) = stat(k,2) + 1.0_DP
          !$acc end atomic
       endif
    enddo
    enddo
    enddo
    !$acc end kernels


    ! All reduce
    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    !$acc host_data use_device(stat, allstat) if(flag_device)
    call MPI_Allreduce( stat,                 &
                        allstat,              &
                        KA * 2,               &
                        MPI_DOUBLE_PRECISION, &
                        MPI_SUM,              &
                        COMM_world_t,         &
                        ierr                  )
    !$acc end host_data
    call PROF_rapend  ('COMM_Allreduce', 2)

    !$acc kernels if(flag_device)
    do k = 1, KA
       zerosw = 0.5_DP - sign(0.5_DP, allstat(k,2) - 1.E-12_DP )
       varmean(k) = allstat(k,1) / ( allstat(k,2) + zerosw ) * ( 1.0_DP - zerosw )
       !LOG_INFO("COMM_horizontal_mean_3D",*) k, varmean(k), allstatval(k), allstatcnt(k)
    enddo
    !$acc end kernels

    !$acc end data

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

    !$acc host_data use_device(send, recv) if(acc_is_present(send))
    call MPI_GATHER( send(:,:),       &
                     sendcounts,      &
                     COMM_datatype_t, &
                     recv(:,:,:),     &
                     recvcounts,      &
                     COMM_datatype_t, &
                     PRC_masterrank,  &
                     COMM_world_t,    &
                     ierr             )
    !$acc end host_data

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

    !$acc host_data use_device(send, recv) if(acc_is_present(send))
    call MPI_GATHER( send(:,:,:),     &
                     sendcounts,      &
                     COMM_datatype_t, &
                     recv(:,:,:,:),   &
                     recvcounts,      &
                     COMM_datatype_t, &
                     PRC_masterrank,  &
                     COMM_world_t,    &
                     ierr             )
    !$acc end host_data

    return
  end subroutine COMM_gather_3D

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in scalar field
  subroutine COMM_bcast_SCR_SP( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    real(SP), intent(inout) :: var  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_SCR_SP
  subroutine COMM_bcast_SCR_DP( var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    real(DP), intent(inout) :: var  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = 1

    call MPI_BCAST( var,            &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_SCR_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D field
  subroutine COMM_bcast_1D_SP( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA       !< dimension size

    real(SP), intent(inout) :: var(IA)  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_1D_SP
  subroutine COMM_bcast_1D_DP( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA       !< dimension size

    real(DP), intent(inout) :: var(IA)  !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_1D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 2D field
  subroutine COMM_bcast_2D_SP( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA, JA     !< dimension size

    real(SP), intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_2D_SP
  subroutine COMM_bcast_2D_DP( IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: IA, JA     !< dimension size

    real(DP), intent(inout) :: var(IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_2D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 3D field
  subroutine COMM_bcast_3D_SP( KA, IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA    !< dimension size

    real(SP), intent(inout) :: var(KA,IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:),     &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_3D_SP
  subroutine COMM_bcast_3D_DP( KA, IA, JA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA    !< dimension size

    real(DP), intent(inout) :: var(KA,IA,JA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = KA * IA * JA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:),     &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_3D_DP

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 4D field
  subroutine COMM_bcast_4D_SP( KA, IA, JA, NT, var )
    use scale_prc, only: &
       PRC_abort, &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA, NT   !< dimension size

    real(SP), intent(inout) :: var(KA,IA,JA,NT) !< broadcast buffer

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

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:,:),   &
                    counts,         &
                    MPI_REAL,       &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_4D_SP
  subroutine COMM_bcast_4D_DP( KA, IA, JA, NT, var )
    use scale_prc, only: &
       PRC_abort, &
       PRC_masterrank
    implicit none

    integer,  intent(in)    :: KA, IA, JA, NT   !< dimension size

    real(DP), intent(inout) :: var(KA,IA,JA,NT) !< broadcast buffer

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

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:,:,:),   &
                    counts,         &
                    MPI_DOUBLE_PRECISION, &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_4D_DP

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
                    COMM_world_t,   &
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
                    COMM_world_t,   &
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

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:,:),       &
                    counts,         &
                    MPI_INTEGER,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

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
                    COMM_world_t,   &
                    ierr            )

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_LOGICAL_SCR

  !-----------------------------------------------------------------------------
  !> Broadcast data for whole process value in 1D (logical)
  subroutine COMM_bcast_LOGICAL_1D( IA, var )
    use scale_prc, only: &
       PRC_masterrank
    implicit none

    integer, intent(in)    :: IA      !< dimension size
    logical, intent(inout) :: var(IA) !< broadcast buffer

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = IA

    !$acc host_data use_device(var) if(acc_is_present(var))
    call MPI_BCAST( var(:),         &
                    counts,         &
                    MPI_LOGICAL,    &
                    PRC_masterrank, &
                    COMM_world_t,   &
                    ierr            )
    !$acc end host_data

    call PROF_rapend('COMM_Bcast', 2)

    return
  end subroutine COMM_bcast_LOGICAL_1D

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
                    COMM_world_t,   &
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

#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
#endif

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

    KA    = ginfo(gid)%KA
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$acc host_data use_device(var)

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype_t,                 &
                        MPI_PROC_NULL, tag+ginfo(gid)%nreq_max+1, COMM_world_t, &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,vid), ierr   )

    !--- From 4-Direction HALO communicate
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_RECV_INIT( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,                &
                           PRC_next(PRC_S), tag+1, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
       ireq = ireq + 1
    end if
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_RECV_INIT( var(:,:,JE+1:JA), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,               &
                           PRC_next(PRC_N), tag+2, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
       ireq = ireq + 1
    end if
    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%recvpack_WE2P(:,:,seqid)
       !$acc host_data use_device(ptr)
#endif

       ! From E
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_RECV_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,2,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_E), tag+3, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
       end if
       ! From W
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_RECV_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,1,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_W), tag+4, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
          ireq = ireq + 1
       end if
       !$acc end host_data
    end if

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,seqid)
       !$acc host_data use_device(ptr)
#endif
       ! To W HALO
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_SEND_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,1,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_W), tag+3, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
          ireq = ireq + 1
       end if
       ! To E HALO
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_SEND_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,2,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_E), tag+4, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
          ireq = ireq + 1
       end if
       !$acc end host_data
    end if
    ! To N HALO
    if ( PRC_HAS_N ) then
       call MPI_SEND_INIT( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,         &
                           PRC_next(PRC_N), tag+1, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
       ireq = ireq + 1
    end if
    ! To S HALO
    if ( PRC_HAS_S ) then
       call MPI_SEND_INIT( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,         &
                           PRC_next(PRC_S), tag+2, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
       ireq = ireq + 1
    end if

    ginfo(gid)%preq_cnt(vid) = ireq - 1

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(vid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,vid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    !$acc end host_data

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

#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
#endif

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

    !$acc host_data use_device(var)

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype_t,                 &
                        MPI_PROC_NULL, tag+ginfo(gid)%nreq_max+1, COMM_world_t, &
                        ginfo(gid)%preq_list(ginfo(gid)%nreq_max+1,vid), ierr   )


    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          tagc = 0
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                       &
                                 PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From SW
          tagc = 10
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                           &
                                  PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NE
          tagc = 20
          do j = JE+1, JA
             call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                       &
                                 PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NW
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                          &
                                 PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From E
          tagc = 60
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,seqid)
          !$acc host_data use_device(ptr)
          call MPI_RECV_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,2,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#endif
                              PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! From W
          tagc = 70
#ifdef _OPENACC
          call MPI_RECV_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,1,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
          !$acc end host_data
          ireq = ireq + 1
       end if
       ! From S
       tagc = 40
       do j = 1, JS-1
          call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                       &
                              PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JA
          call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                       &
                              PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !--- To 8-Direction HALO communicate
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                   &
                          PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                   &
                          PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       if ( .not. PRC_TwoD ) then
          ! To W HALO
          tagc = 60
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,seqid)
          !$acc host_data use_device(ptr)
          call MPI_SEND_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,1,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
          ireq = ireq + 1
          ! To E HALO
          tagc = 70
#ifdef _OPENACC
          call MPI_SEND_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
          call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,2,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                              PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
          !$acc end host_data
          ireq = ireq + 1
          ! To NW HALO
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                 PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To NE HALO
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                 &
                                 PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SW HALO
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                 PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SE HALO
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                 &
                                 PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
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
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                       &
                                    PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                    PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                    PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From SW
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                          &
                                    PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NE
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                       &
                                    PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                    PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                    PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NW
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                          &
                                    PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_RECV_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,seqid)
          !$acc host_data use_device(ptr)
#endif
          ! From E
          if ( PRC_HAS_E ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_RECV_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,2,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                                 PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_RECV_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_RECV_INIT( ginfo(gid)%recvpack_WE2P(:,1,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                                 PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
             ireq = ireq + 1
          endif
          !$acc end host_data
       end if
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = 1, JS-1
             call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                   &
                             PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JA
             call MPI_RECV_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                   &
                             PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                   &
                             PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                   &
                             PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,seqid)
          !$acc host_data use_device(ptr)
#endif
          ! To W HALO
          if ( PRC_HAS_W ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_SEND_INIT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,1,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                                 PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_SEND_INIT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,   &
#else
             call MPI_SEND_INIT( ginfo(gid)%sendpack_P2WE(:,2,seqid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                                 PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr  )
             ireq = ireq + 1
          endif
          !$acc end host_data
          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                    PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                 &
                                    PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_SEND_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                    PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                &
                                    PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                    PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_SEND_INIT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                    PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                 &
                                    PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_SEND_INIT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                    PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_SEND_INIT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                &
                                    PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%preq_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
       end if

    endif

    ginfo(gid)%preq_cnt(vid) = ireq - 1

    ! to finish initial processes of MPI
    nreq = ginfo(gid)%preq_cnt(vid)
    do i = 1, 32
       call MPI_TESTALL( nreq, ginfo(gid)%preq_list(1:nreq,vid), &
                         flag, MPI_STATUSES_IGNORE, ierr         )
    enddo

    !$acc end host_data

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
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ierr
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    ireq = 1

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(vid) ) then
       LOG_ERROR("vars_3D_mpi",*) 'packing buffer is already used', vid
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(vid) = .true.
#endif

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    !$acc host_data use_device(var) if(flag_device)

    !--- From 4-Direction HALO communicate
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_IRECV( var(:,:,1:JS-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,               &
                       PRC_next(PRC_S), tag+1, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_IRECV( var(:,:,JE+1:JA), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,              &
                       PRC_next(PRC_N), tag+2, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif
       ! From E
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+3, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr    )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+4, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr    )
          ireq = ireq + 1
       endif
       !$acc end host_data
    end if

    !$acc end host_data

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
       call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)
    end if

    !$acc host_data use_device(var) if(flag_device)

    ! To N HALO
    if ( PRC_HAS_N ) then
       call MPI_ISEND( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,        &
                       PRC_next(PRC_N), tag+1, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    ! To S HALO
    if ( PRC_HAS_S ) then
       call MPI_ISEND( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t,        &
                       PRC_next(PRC_S), tag+2, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif

    !$acc end host_data

    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc wait
       !$acc host_data use_device(ptr) if(flag_device)
#endif
       ! To W HALO
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+3, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr    )
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+4, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr    )
          ireq = ireq + 1
       endif

       !$acc end host_data
    end if

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars_3D_mpi

  subroutine vars_3D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication
    integer,  intent(in)    :: gid        !< grid ID
    integer,  intent(in)    :: vid        !< request ID

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)

    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then
       call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)
    end if

    !$acc host_data use_device(var)

    ! To N HALO
    if ( PRC_HAS_N ) then
       disp = 0
       call MPI_PUT( var(:,:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
    endif
    ! To S HALO
    if ( PRC_HAS_S ) then
       disp = KA * IA * JHALO
       call MPI_PUT( var(:,:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS4*KA, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
    endif

    !$acc end host_data

    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc wait
       !$acc host_data use_device(ptr)
#endif

       ! To W HALO
       if ( PRC_HAS_W ) then
          disp = 1
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
       endif

       !$acc end host_data
    end if

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )

    !$acc end data

    return
  end subroutine vars_3D_mpi_onesided

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
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ierr
    integer :: j
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    tag  = ( (gid - 1) * COMM_vsize_max + vid ) * 100
    tag  = vid * 100
    ireq = 1

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(vid) ) then
       LOG_ERROR("vars8_3D_mpi",*) 'packing buffer is already used', vid
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(vid) = .true.
#endif

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !$acc host_data use_device(var) if(flag_device)

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          tagc = 0
          do j = 1, JS-1
             call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                             PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From SW
          tagc = 10
          do j = 1, JS-1
             call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NE
          tagc = 20
          do j = JE+1, JA
             call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                             PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! From NW
          tagc = 30
          do j = JE+1, JA
             call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                             PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
          !$acc host_data use_device(ptr) if(flag_device)
#endif
          ! From E
          tagc = 60
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! From W
          tagc = 70
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          !$acc end host_data
       end if
       ! From S
       tagc = 40
       do j = 1, JS-1
          call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                          PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JA
          call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                          PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !--- To 8-Direction HALO communicate
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                          PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                          PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)

          !$acc host_data use_device(var) if(flag_device)

          ! To NW HALO
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                             PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To NE HALO
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                &
                             PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SW HALO
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                             PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
          ! To SE HALO
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                &
                             PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc wait
          !$acc host_data use_device(ptr) if(flag_device)
#endif
          ! To W HALO
          tagc = 60
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
          ! To E HALO
          tagc = 70
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          !$acc end host_data
          ireq = ireq + 1

       end if

    else ! non-periodic condition

       !$acc host_data use_device(var) if(flag_device)

       !--- From 8-Direction HALO communicate
       if ( .not. PRC_TwoD ) then
          ! From SE
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                     &
                                PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                     &
                                PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From SW
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NE
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                      &
                                PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                     &
                                PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_IRECV( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                     &
                                PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! From NW
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                         &
                                PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_IRECV( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
#ifdef _OPENACC
          ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
          !$acc host_data use_device(ptr) if(flag_device)
#endif
          ! From E
          if ( PRC_HAS_E ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                             PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! From W
          if ( PRC_HAS_W ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                             PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          !$acc end host_data
       end if
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = 1, JS-1
             call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                             PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JA
             call MPI_IRECV( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                             PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                             PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t,                      &
                             PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)

          !$acc host_data use_device(var) if(flag_device)

          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_NW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 20
             do j = JE+1, JA
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                       &
                                PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             tagc = 10
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                &
                                PRC_next(PRC_NE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_N ) then
             tagc = 0
             do j = JE-JHALO+1, JE
                call MPI_ISEND( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                     &
                                PRC_next(PRC_N), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 30
             do j = JE+1, JA
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,               &
                                PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_SW), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                        &
                                PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_W ) then
             tagc = 0
             do j = 1, JS-1
                call MPI_ISEND( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                       &
                                PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             tagc = 30
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                &
                                PRC_next(PRC_SE), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_S ) then
             tagc = 20
             do j = JS, JS+JHALO-1
                call MPI_ISEND( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,                     &
                                PRC_next(PRC_S), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          else if ( PRC_HAS_E ) then
             tagc = 10
             do j = 1, JS-1
                call MPI_ISEND( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t,               &
                                PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                ireq = ireq + 1
                tagc = tagc + 1
             enddo
          endif

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc wait
          !$acc host_data use_device(ptr) if(flag_device)
#endif

          ! To W HALO
          if ( PRC_HAS_W ) then
             tagc = 60
#ifdef _OPENACC
             call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                             PRC_next(PRC_W), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             tagc = 70
#ifdef _OPENACC
             call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t,    &
#else
             call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                             PRC_next(PRC_E), tag+tagc, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
          endif
          !$acc end host_data

       end if

    endif

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_3D_mpi

  subroutine vars8_3D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
    integer :: j
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    KA    = ginfo(gid)%KA
    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)
    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !$acc host_data use_device(var)

       !--- To 8-Direction HALO communicate
       ! To N HALO
       do j = JE-JHALO+1, JE
          disp = KA * ( IHALO + IA * ( j - JE+JHALO-1 ) )
          call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packNS(vid), ierr )
       enddo
       ! To S HALO
       do j = JS, JS+JHALO-1
          disp = KA * ( IHALO + IA * ( j - JS + JHALO ) )
          call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packNS(vid), ierr )
       enddo

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)

          !$acc host_data use_device(var)
          ! To NW HALO
          do j = JE-JHALO+1, JE
             disp = KA * ( IE + IA * ( j - JE+JHALO-1 ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
          ! To NE HALO
          do j = JE-JHALO+1, JE
             disp = KA * ( IA * ( j - JE+JHALO-1 ) )
             call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
          ! To SW HALO
          do j = JS, JS+JHALO-1
             disp = KA * ( IE + IA * ( j - JS + JHALO ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
          ! To SE HALO
          do j = JS, JS+JHALO-1
             disp = KA * ( IA * ( j - JS + JHALO ) )
             call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc wait
          !$acc host_data use_device(ptr)
#endif
          ! To W HALO
          disp = 1
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          ! To E HALO
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                        PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          !$acc end host_data

       end if

    else ! non-periodic condition

       !$acc host_data use_device(var)

       !--- To 8-Direction HALO communicate
       ! To N HALO
       if ( PRC_HAS_N ) then
          do j = JE-JHALO+1, JE
             disp = KA * ( IHALO + IA * ( j - JE+JHALO-1 ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          do j = JS, JS+JHALO-1
             disp = KA * ( IHALO + IA * ( j - JS + JHALO ) )
             call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
       endif

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)

          !$acc host_data use_device(var)

          ! To NW HALO
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IE + IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_W ) then
             do j = JE+1, JA
                disp = KA * ( IE + IA * ( j - JE-1 + JHALO ) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif
          ! To NE HALO
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = KA * ( IE + IA * ( j - JE+JHALO-1 ) )
                call MPI_PUT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_E ) then
             do j = JE+1, JA
                disp = KA * IA * ( j - JE-1 + JHALO )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif
          ! To SW HALO
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IE + IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_W ) then
             do j = 1, JS-1
                disp = KA * ( IE + IA * (j-1) )
                call MPI_PUT( var(1,IS,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif
          ! To SE HALO
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = KA * ( IE + IA * ( j - JS + JHALO ) )
                call MPI_PUT( var(1,IE+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_E ) then
             do j = 1, JS-1
                disp = KA * ( IA * ( j - 1 ) )
                call MPI_PUT( var(1,IE-IHALO+1,j), ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C*KA, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif

          !$acc end host_data

#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc wait
          !$acc host_data use_device(ptr)
#endif

          ! To W HALO
          if ( PRC_HAS_W ) then
             disp = 1
#ifdef _OPENACC
             call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                           PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
          endif
          ! To E HALO
          if ( PRC_HAS_E ) then
             disp = 0
#ifdef _OPENACC
             call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
#endif
                           PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE*KA, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
          endif
          !$acc end host_data

       end if

    endif

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )

    !$acc end data

    return
  end subroutine vars8_3D_mpi_onesided

  subroutine vars_2D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ireq, tag
    integer :: ierr
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
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

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    !$acc host_data use_device(var) if(flag_device)

    !--- From 4-Direction HALO communicate
    ! From S
    if ( PRC_HAS_S ) then
       call MPI_IRECV( var(:,1:JS-1), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                       PRC_next(PRC_S), tag+1, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    ! From N
    if ( PRC_HAS_N ) then
       call MPI_IRECV( var(:,JE+1:JA), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                       PRC_next(PRC_N), tag+2, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif

    if ( .not. PRC_TwoD ) then
#ifdef _OPENACC
       ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif
       ! From E
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+3, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+4, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       !$acc end host_data
    end if

    !$acc end host_data

    !--- To 4-Direction HALO communicate
    if ( .not. PRC_TwoD ) then

       call packWE_2D( IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)

#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr) if(flag_device)
#endif

       ! To W HALO communicate
       if ( PRC_HAS_W ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                          PRC_next(PRC_W), tag+3, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To E HALO communicate
       if ( PRC_HAS_E ) then
#ifdef _OPENACC
          call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                          PRC_next(PRC_E), tag+4, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

       !$acc end host_data

    end if

    !$acc host_data use_device(var) if(flag_device)

    ! To N HALO communicate
    if ( PRC_HAS_N ) then
       call MPI_ISEND( var(:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
            PRC_next(PRC_N), tag+1, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif
    ! To S HALO communicate
    if ( PRC_HAS_S ) then
       call MPI_ISEND( var(:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
            PRC_next(PRC_S), tag+2, COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
       ireq = ireq + 1
    endif

    !$acc end host_data

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars_2D_mpi

  subroutine vars_2D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)

    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    !--- To 4-Direction HALO communicate

    if ( .not. PRC_TwoD ) then

       call packWE_2D( IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, vid)

#ifdef _OPENACC
       ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
       !$acc host_data use_device(ptr)
#endif

       ! To W HALO communicate
       if ( PRC_HAS_W ) then
          disp = 1
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
       endif
       ! To E HALO communicate
       if ( PRC_HAS_E ) then
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                        PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
       endif

       !$acc end host_data

    end if

    !$acc host_data use_device(var)

    ! To N HALO communicate
    if ( PRC_HAS_N ) then
       disp = 0
       call MPI_PUT( var(:,JE-JHALO+1:JE), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
    endif
    ! To S HALO communicate
    if ( PRC_HAS_S ) then
       disp = IA * JHALO
       call MPI_PUT( var(:,JS:JS+JHALO-1), ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS4, COMM_datatype_t, &
                     ginfo(gid)%win_packNS(vid), ierr )
    endif

    !$acc end host_data

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )

    !$acc end data

    return
  end subroutine vars_2D_mpi_onesided

  subroutine vars8_2D_mpi(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    integer :: ireq, tag, tagc

    integer :: ierr
    integer :: j
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
    logical :: flag_device
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
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

#ifdef _OPENACC
    flag_device = acc_is_present(var)
#endif

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 8-Direction HALO communicate
        !$acc host_data use_device(var) if(flag_device)

        if ( .not. PRC_TwoD ) then
           ! From SE
           tagc = 0
           do j = 1, JS-1
              call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                              COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From SW
           tagc = 10
           do j = 1, JS-1
              call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                              COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From NE
           tagc = 20
           do j = JE+1, JA
              call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                              COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
           ! From NW
           tagc = 30
           do j = JE+1, JA
              call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                              COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
#ifdef _OPENACC
           ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif
           ! From E
#ifdef _OPENACC
           call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
           call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_E), tag+60,                &
                           COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
           ireq = ireq + 1
           ! From W
#ifdef _OPENACC
           call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
           call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_W), tag+70,                &
                           COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
           ireq = ireq + 1
           !$acc end host_data
        end if
        ! From S
        tagc = 40
        do j = 1, JS-1
            call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                            COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                            COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
        enddo
        ! From N
        tagc = 50
        do j = JE+1, JA
            call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                            COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                            COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo


        !--- To 8-Direction HALO communicate

        ! To N HALO communicate
        tagc = 40
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                            COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                            COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To S HALO communicate
        tagc = 50
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                            COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                            COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        !$acc end host_data

        if ( .not. PRC_TwoD ) then

           call packWE_2D( IA, IS, IE, JA, JS, JE, &
                           IHALO, &
                           var, gid, vid)

           !$acc host_data use_device(var) if(flag_device)
#ifdef _OPENACC
           ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif

           ! To W HALO communicate
#ifdef _OPENACC
           call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
           call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_W), tag+60,                &
                           COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
           ireq = ireq + 1

           ! To E HALO communicate
#ifdef _OPENACC
           call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
           call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                           COMM_datatype_t, PRC_next(PRC_E), tag+70,                &
                           COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
           ireq = ireq + 1
           !$acc end host_data

           ! To NW HALO communicate
           tagc = 0
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                  &
                              COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To NE HALO communicate
           tagc = 10
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,          &
                              COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To SW HALO communicate
           tagc = 20
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                  &
                              COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           ! To SE HALO communicate
           tagc = 30
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,          &
                              COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,      &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo

           !$acc end host_data

        end if

    else
    !--- non-periodic condition
        !--- From 8-Direction HALO communicate

        !$acc host_data use_device(var) if(flag_device)

        if ( .not. PRC_TwoD ) then
           ! From SE
           if ( PRC_HAS_S .AND. PRC_HAS_E ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From SW
           if ( PRC_HAS_S .AND. PRC_HAS_W ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From NE
           if ( PRC_HAS_N .AND. PRC_HAS_E ) then
              tagc = 20
              do j = JE+1, JE+JHALO
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_IRECV( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! From NW
           if ( PRC_HAS_N .AND. PRC_HAS_W ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_IRECV( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

#ifdef _OPENACC
           ptr => ginfo(gid)%recvpack_WE2P(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif
           ! From E
           if ( PRC_HAS_E ) then
#ifdef _OPENACC
              call MPI_IRECV( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
              call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_E), tag+60,                &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
              ireq = ireq + 1
           endif

           ! From W
           if ( PRC_HAS_W ) then
#ifdef _OPENACC
              call MPI_IRECV( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
              call MPI_IRECV( ginfo(gid)%recvpack_WE2P(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_W), tag+70,                &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
              ireq = ireq + 1
           endif
           !$acc end host_data

        end if

        ! From S
        if ( PRC_HAS_S ) then
           tagc = 40
           do j = 1, JS-1
              call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                              COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        ! From N
        if ( PRC_HAS_N ) then
           tagc = 50
           do j = JE+1, JA
              call MPI_IRECV( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                              COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif


        !! RECEIVE

        ! To N HALO communicate
        if ( PRC_HAS_N ) then
           tagc = 40
           do j = JE-JHALO+1, JE
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                              COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        ! To S HALO communicate
        if ( PRC_HAS_S ) then
           tagc = 50
           do j = JS, JS+JHALO-1
              call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_NS8,                 &
                              COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
              ireq = ireq + 1
              tagc = tagc + 1
           enddo
        endif

        !$acc end host_data

        if ( .not. PRC_TwoD ) then

           call packWE_2D( IA, IS, IE, JA, JS, JE, &
                           IHALO, &
                           var, gid, vid)

           !$acc host_data use_device(var) if(flag_device)
#ifdef _OPENACC
           ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
           !$acc host_data use_device(ptr) if(flag_device)
#endif

           ! To W HALO communicate
           if ( PRC_HAS_W ) then
#ifdef _OPENACC
              call MPI_ISEND( ptr(:,1), ginfo(gid)%size2D_WE, &
#else
              call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_W), tag+60,                &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
              ireq = ireq + 1
           endif

           ! To E HALO communicate
           if ( PRC_HAS_E ) then
#ifdef _OPENACC
              call MPI_ISEND( ptr(:,2), ginfo(gid)%size2D_WE, &
#else
              call MPI_ISEND( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, &
#endif
                              COMM_datatype_t, PRC_next(PRC_E), tag+70,                &
                              COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr        )
              ireq = ireq + 1
           endif
           !$acc end host_data

           ! To NW HALO communicate
           if ( PRC_HAS_N .AND. PRC_HAS_W ) then
              tagc = 0
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                  &
                                 COMM_datatype_t, PRC_next(PRC_NW), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 10
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 20
              do j = JE+1, JA
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                  &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To NE HALO communicate
           if ( PRC_HAS_N .AND. PRC_HAS_E ) then
              tagc = 10
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,          &
                                 COMM_datatype_t, PRC_next(PRC_NE), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_N ) then
              tagc = 0
              do j = JE-JHALO+1, JE
                 call MPI_ISEND( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_N), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 30
              do j = JE+1, JA
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,          &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To SW HALO communicate
           if ( PRC_HAS_S .AND. PRC_HAS_W ) then
              tagc = 20
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                  &
                                 COMM_datatype_t, PRC_next(PRC_SW), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 30
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(1,j), ginfo(gid)%size2D_4C,                   &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_W ) then
              tagc = 0
              do j = 1, JS-1
                 call MPI_ISEND( var(IS,j), ginfo(gid)%size2D_4C,                  &
                                 COMM_datatype_t, PRC_next(PRC_W), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           ! To SE HALO communicate
           if ( PRC_HAS_S .AND. PRC_HAS_E ) then
              tagc = 30
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,          &
                                 COMM_datatype_t, PRC_next(PRC_SE), tag+tagc,      &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_S ) then
              tagc = 20
              do j = JS, JS+JHALO-1
                 call MPI_ISEND( var(IE+1,j), ginfo(gid)%size2D_4C,                &
                                 COMM_datatype_t, PRC_next(PRC_S), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           else if ( PRC_HAS_E ) then
              tagc = 10
              do j = 1, JS-1
                 call MPI_ISEND( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C,          &
                                 COMM_datatype_t, PRC_next(PRC_E), tag+tagc,       &
                                 COMM_world_t, ginfo(gid)%req_list(ireq,vid), ierr )
                 ireq = ireq + 1
                 tagc = tagc + 1
              enddo
           endif

           !$acc end host_data

        end if

    endif

    ginfo(gid)%req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_2D_mpi

  subroutine vars8_2D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: IA, IS, IE, IHALO
    integer :: JA, JS, JE, JHALO

    integer(kind=MPI_ADDRESS_KIND) :: disp

    integer :: ierr
    integer :: j
#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:)
#endif
    !---------------------------------------------------------------------------

    IA    = ginfo(gid)%IA
    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JA    = ginfo(gid)%JA
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$acc data copyin(var)

    call MPI_Win_start( group_packWE, 0, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_start( group_packNS, 0, ginfo(gid)%win_packNS(vid), ierr )

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition

       !--- To 8-Direction HALO communicate

       !$acc host_data use_device(var)

       ! To N HALO communicate
       do j = JE-JHALO+1, JE
          disp = IHALO + IA * ( j - JE+JHALO-1 )
          call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                        PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                        ginfo(gid)%win_packNS(vid), ierr )
       enddo
       ! To S HALO communicate
       do j = JS, JS+JHALO-1
          disp = IHALO + IA * ( j - JS + JHALO )
          call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
               PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
               ginfo(gid)%win_packNS(vid), ierr )
       enddo

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call packWE_2D( IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)

          !$acc host_data use_device(var)
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif

          ! To W HALO communicate
          disp = 1
#ifdef _OPENACC
          call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                        PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                        ginfo(gid)%win_packWE(vid), ierr )
          ! To E HALO communicate
          disp = 0
#ifdef _OPENACC
          call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
          call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
               PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
               ginfo(gid)%win_packWE(vid), ierr )
          !$acc end host_data
          ! To NW HALO communicate
          do j = JE-JHALO+1, JE
             disp = IE + IA * ( j - JE+JHALO-1 )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
          ! To NE HALO communicate
          do j = JE-JHALO+1, JE
             disp = IA * ( j - JE+JHALO-1 )
             call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
          ! To SW HALO communicate
          do j = JS, JS+JHALO-1
             disp = IE + IA * ( j - JS + JHALO )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
          ! To SE HALO communicate
          do j = JS, JS+JHALO-1
             disp = IA * ( j - JS + JHALO )
             call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo

          !$acc end host_data

       end if
    else
    !--- non-periodic condition

       !$acc host_data use_device(var)

       ! To N HALO communicate
       if ( PRC_HAS_N ) then
          do j = JE-JHALO+1, JE
             disp = IHALO + IA * ( j - JE+JHALO-1 )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           PRC_next(PRC_N), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
       endif
       ! To S HALO communicate
       if ( PRC_HAS_S ) then
          do j = JS, JS+JHALO-1
             disp = IHALO + IA * ( j - JS + JHALO )
             call MPI_PUT( var(IS,j), ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           PRC_next(PRC_S), disp, ginfo(gid)%size2D_NS8, COMM_datatype_t, &
                           ginfo(gid)%win_packNS(vid), ierr )
          enddo
       endif

       !$acc end host_data

       if ( .not. PRC_TwoD ) then

          call packWE_2D( IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, gid, vid)

          !$acc host_data use_device(var)
#ifdef _OPENACC
          ptr => ginfo(gid)%sendpack_P2WE(:,:,vid)
          !$acc host_data use_device(ptr)
#endif

          ! To W HALO communicate
          if ( PRC_HAS_W ) then
             disp = 1
#ifdef _OPENACC
             call MPI_PUT( ptr(:,1), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,1,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                           PRC_next(PRC_W), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
          endif
          ! To E HALO communicate
          if ( PRC_HAS_E ) then
             disp = 0
#ifdef _OPENACC
             call MPI_PUT( ptr(:,2), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#else
             call MPI_PUT( ginfo(gid)%sendpack_P2WE(:,2,vid), ginfo(gid)%size2D_WE, COMM_datatype_t, &
#endif
                           PRC_next(PRC_E), disp, ginfo(gid)%size2D_WE, COMM_datatype_t, &
                           ginfo(gid)%win_packWE(vid), ierr )
          endif
          !$acc end host_data
          ! To NW HALO communicate
          if ( PRC_HAS_N .AND. PRC_HAS_W ) then
             do j = JE-JHALO+1, JE
                disp = IE + IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_NW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_W ) then
             do j = JE+1, JA
                disp = IE + IA * ( j - JE-1 + JHALO )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif
          ! To NE HALO communicate
          if ( PRC_HAS_N .AND. PRC_HAS_E ) then
             do j = JE-JHALO+1, JE
                disp = IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_NE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_N ) then
             do j = JE-JHALO+1, JE
                disp = IE + IA * ( j - JE+JHALO-1 )
                call MPI_PUT( var(IE+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_N), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_E ) then
             do j = JE+1, JA
                disp = IA * ( j - JE-1 + JHALO )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif
          ! To SW HALO communicate
          if ( PRC_HAS_S .AND. PRC_HAS_W ) then
             do j = JS, JS+JHALO-1
                disp = IE + IA * ( j - JS + JHALO )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_SW), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = IA * ( j - JS + JHALO )
                call MPI_PUT( var(1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_W ) then
             do j = 1, JS-1
                disp = IE + IA * ( j - 1 )
                call MPI_PUT( var(IS,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_W), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif
          ! To SE HALO communicate
          if ( PRC_HAS_S .AND. PRC_HAS_E ) then
             do j = JS, JS+JHALO-1
                disp = IA * ( j - JS + JHALO )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_SE), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_S ) then
             do j = JS, JS+JHALO-1
                disp = IE + IA * ( j - JS + JHALO )
                call MPI_PUT( var(IE+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_S), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          else if ( PRC_HAS_E ) then
             do j = 1, JS-1
                disp = IA * ( j - 1 )
                call MPI_PUT( var(IE-IHALO+1,j), ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              PRC_next(PRC_E), disp, ginfo(gid)%size2D_4C, COMM_datatype_t, &
                              ginfo(gid)%win_packNS(vid), ierr )
             enddo
          endif

          !$acc end host_data
       end if

    endif

    call MPI_Win_complete( ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_complete( ginfo(gid)%win_packNS(vid), ierr )

    !$acc end data

    return
  end subroutine vars8_2D_mpi_onesided

  subroutine vars_3D_mpi_pc(var, gid, vid)
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO

    integer :: ierr
    !---------------------------------------------------------------------------

#ifdef DEBUG
    if ( ginfo(gid)%use_packbuf(ginfo(gid)%packid(vid)) ) then
       LOG_ERROR("vars_3D_mpi_pc",*) 'packing buffer is already used', vid, ginfo(gid)%packid(vid)
       call PRC_abort
    end if
    ginfo(gid)%use_packbuf(ginfo(gid)%packid(vid)) = .true.
#endif

#ifdef _OPENACC
    if ( ginfo(gid)%device_alloc(vid+COMM_vsize_max) ) then
       !$acc update device(var)
    end if
#endif

    if ( .not. PRC_TwoD ) then
       KA = ginfo(gid)%KA
       IA = ginfo(gid)%IA
       IS = ginfo(gid)%IS
       IE = ginfo(gid)%IE
       JA = ginfo(gid)%JA
       JS = ginfo(gid)%JS
       JE = ginfo(gid)%JE
       IHALO = ginfo(gid)%IHALO
       call packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                       IHALO, &
                       var, gid, ginfo(gid)%packid(vid))
       !$acc wait
    end if

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

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call MPI_WAITALL( ginfo(gid)%req_cnt (vid),                           &
                      ginfo(gid)%req_list(1:ginfo(gid)%req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                &
                      ierr                                                )
    if ( .not. PRC_TwoD ) then
       KA = ginfo(gid)%KA
       IA = ginfo(gid)%IA
       IS = ginfo(gid)%IS
       IE = ginfo(gid)%IE
       JA = ginfo(gid)%JA
       JS = ginfo(gid)%JS
       JE = ginfo(gid)%JE
       IHALO = ginfo(gid)%IHALO
       call unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, ginfo(gid)%recvpack_WE2P(:,:,vid) )
       !$acc wait
    end if

#ifdef DEBUG
    ginfo(gid)%use_packbuf(vid) = .false.
#endif

    return
  end subroutine wait_3D_mpi

  subroutine wait_3D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)     :: gid
    integer, intent(in)     :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    real(RP), pointer :: pack(:)

    integer :: ierr
    !---------------------------------------------------------------------------

    KA = ginfo(gid)%KA
    IA = ginfo(gid)%IA
    IS = ginfo(gid)%IS
    IE = ginfo(gid)%IE
    JA = ginfo(gid)%JA
    JS = ginfo(gid)%JS
    JE = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    call MPI_Win_wait( ginfo(gid)%win_packWE(vid), ierr )
    if ( .not. PRC_TwoD ) then
       call C_F_pointer( ginfo(gid)%recvbuf_WE(vid), pack, (/ginfo(gid)%size2D_WE*KA*2/) )
       call unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, pack )
    end if

    call MPI_Win_wait( ginfo(gid)%win_packNS(vid), ierr )
    call C_F_pointer( ginfo(gid)%recvbuf_NS(vid), pack, (/ginfo(gid)%size2D_NS4*KA*2/) )
    call unpackNS_3D( KA, IA, IS, IE, JA, JS, JE, &
                      JHALO, &
                      var, pack )

    !$acc wait

    call MPI_Win_post( group_packWE, MPI_MODE_NOSTORE, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_post( group_packNS, MPI_MODE_NOSTORE, ginfo(gid)%win_packNS(vid), ierr )

    return
  end subroutine wait_3D_mpi_onesided

  subroutine wait_2D_mpi(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call MPI_WAITALL( ginfo(gid)%req_cnt(vid),                            &
                      ginfo(gid)%req_list(1:ginfo(gid)%req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                &
                      ierr                                                )
    if ( .not. PRC_TwoD ) then
       KA = ginfo(gid)%KA
       IA = ginfo(gid)%IA
       IS = ginfo(gid)%IS
       IE = ginfo(gid)%IE
       JA = ginfo(gid)%JA
       JS = ginfo(gid)%JS
       JE = ginfo(gid)%JE
       IHALO = ginfo(gid)%IHALO
       call unpackWE_2D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, ginfo(gid)%recvpack_WE2P(:,:,vid) )
    end if

#ifdef DEBUG
    ginfo(gid)%use_packbuf(vid) = .false.
#endif

    return
  end subroutine wait_2D_mpi

  subroutine wait_2D_mpi_onesided(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO, JHALO

    real(RP), pointer :: pack(:)

    integer :: ierr
    !---------------------------------------------------------------------------

    KA = ginfo(gid)%KA
    IA = ginfo(gid)%IA
    IS = ginfo(gid)%IS
    IE = ginfo(gid)%IE
    JA = ginfo(gid)%JA
    JS = ginfo(gid)%JS
    JE = ginfo(gid)%JE
    IHALO = ginfo(gid)%IHALO
    JHALO = ginfo(gid)%JHALO

    call MPI_Win_wait( ginfo(gid)%win_packWE(vid), ierr )
    if ( .not. PRC_TwoD ) then
       call C_F_pointer( ginfo(gid)%recvbuf_WE(vid), pack, (/ginfo(gid)%size2D_WE*KA*2/)  )
       call unpackWE_2D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, pack )
    end if

    call MPI_Win_wait( ginfo(gid)%win_packNS(vid), ierr )
    call C_F_pointer( ginfo(gid)%recvbuf_NS(vid), pack, (/ginfo(gid)%size2D_NS4*KA*2/)  )
    call unpackNS_2D( IA, IS, IE, JA, JS, JE, &
                      JHALO, &
                      var, pack )

    call MPI_Win_post( group_packWE, MPI_MODE_NOSTORE, ginfo(gid)%win_packWE(vid), ierr )
    call MPI_Win_post( group_packNS, MPI_MODE_NOSTORE, ginfo(gid)%win_packNS(vid), ierr )

    return
  end subroutine wait_2D_mpi_onesided

  subroutine wait_3D_mpi_pc(var, gid, vid)
    use scale_prc_cartesC, only: &
       PRC_TwoD
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: gid
    integer,  intent(in)    :: vid

    integer :: KA
    integer :: IA, IS, IE
    integer :: JA, JS, JE
    integer :: IHALO

    integer :: pid
    integer :: ierr

    !--- wait packets
    call MPI_WAITALL( ginfo(gid)%preq_cnt (vid),                            &
                      ginfo(gid)%preq_list(1:ginfo(gid)%preq_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,                                  &
                      ierr                                                  )
    if ( .not. PRC_TwoD ) then
       KA = ginfo(gid)%KA
       IA = ginfo(gid)%IA
       IS = ginfo(gid)%IS
       IE = ginfo(gid)%IE
       JA = ginfo(gid)%JA
       JS = ginfo(gid)%JS
       JE = ginfo(gid)%JE
       IHALO = ginfo(gid)%IHALO
       pid = ginfo(gid)%packid(vid)
       call unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                         IHALO, &
                         var, ginfo(gid)%recvpack_WE2P(:,:,pid) )
       !$acc wait
    end if

#ifdef DEBUG
    ginfo(gid)%use_packbuf(ginfo(gid)%packid(vid)) = .false.
#endif

#ifdef _OPENACC
    if ( ginfo(gid)%device_alloc(vid+COMM_vsize_max) ) then
       !$acc update host(var)
    end if
#endif

    return
  end subroutine wait_3D_mpi_pc

  subroutine packWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                        IHALO, &
                        var, gid, vid)
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(in) :: var(KA,IA,JA)
    integer,  intent(in) :: gid
    integer,  intent(in) :: vid

    integer :: k, i, j, n

#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:,:)
    ptr => ginfo(gid)%sendpack_P2WE
#endif

    !$acc data copyin(var) if(acc_is_present(var))

    call PROF_rapstart('COMM_pack', 3)

    if ( PRC_HAS_W ) then
       !--- packing packets to West
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = IS, IS+IHALO-1
       !$acc loop independent vector
       do k = 1, KA
          n = (j-JS) * KA * IHALO &
            + (i-IS) * KA         &
            + k
#ifdef _OPENACC
          ptr(n,1,vid) = var(k,i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,1,vid) = var(k,i,j)
#endif
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    if ( PRC_HAS_E ) then
       !--- packing packets to East
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = IE-IHALO+1, IE
       !$acc loop independent vector
       do k = 1, KA
          n = (j-JS)         * KA * IHALO &
            + (i-IE+IHALO-1) * KA         &
            + k
#ifdef _OPENACC
          ptr(n,2,vid) = var(k,i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,2,vid) = var(k,i,j)
#endif
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    call PROF_rapend('COMM_pack', 3)

    !$acc end data

    return
  end subroutine packWE_3D

  subroutine packWE_2D( IA, IS, IE, JA, JS, JE, &
                        IHALO, &
                        var, gid, vid)
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(in) :: var(IA,JA)
    integer,  intent(in) :: vid
    integer,  intent(in) :: gid

    integer :: i, j, n

#ifdef _OPENACC
    real(RP), pointer :: ptr(:,:,:)
    ptr => ginfo(gid)%sendpack_P2WE
#endif
    !$acc data copyin(var) if(acc_is_present(var))

    call PROF_rapstart('COMM_pack', 3)

    if ( PRC_HAS_W ) then
       !--- To 4-Direction HALO communicate
       !--- packing packets to West
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
       !$acc loop independent
       do j = JS, JE
       !$acc loop independent
       do i = IS, IS+IHALO-1
          n = (j-JS) * IHALO &
            + (i-IS) + 1
#ifdef _OPENACC
          ptr(n,1,vid) = var(i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,1,vid) = var(i,j)
#endif
       enddo
       enddo
       !$acc end kernels
    end if

    if ( PRC_HAS_E ) then
       !--- packing packets to East
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
       !$acc loop independent
       do j = JS, JE
       !$acc loop independent
       do i = IE-IHALO+1, IE
          n = (j-JS)         * IHALO &
            + (i-IE+IHALO-1) + 1
#ifdef _OPENACC
          ptr(n,2,vid) = var(i,j)
#else
          ginfo(gid)%sendpack_P2WE(n,2,vid) = var(i,j)
#endif
       enddo
       enddo
       !$acc end kernels
    end if

    !$acc wait

    call PROF_rapend('COMM_pack', 3)

    !$acc end data

    return
  end subroutine packWE_2D

  subroutine unpackWE_3D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, buf )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in)    :: buf(KA,IHALO,JS:JE,2)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf) if(acc_is_present(var))

    call PROF_rapstart('COMM_unpack', 3)

    if ( PRC_HAS_E ) then
       !--- unpacking packets from East
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = IE+1, IA
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf(k,i-IE,j,2)
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    if ( PRC_HAS_W ) then
       !--- unpacking packets from West
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc parallel if(acc_is_present(var)) async
       !$acc loop collapse(2) gang
       do j = JS, JE
       do i = 1, IS-1
       !$acc loop vector
       do k = 1, KA
          var(k,i,j) = buf(k,i,j,1)
       enddo
       enddo
       enddo
       !$acc end parallel
    end if

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackWE_3D

  subroutine unpackWE_2D( KA, IA, IS, IE, JA, JS, JE, &
                          IHALO, &
                          var, buf )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: IHALO
    real(RP), intent(inout) :: var(IA,JA)
    real(RP), intent(in)    :: buf(IHALO,JS:JE,KA,2)

    integer :: i, j
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf) if(acc_is_present(var))

    call PROF_rapstart('COMM_unpack', 3)

    if( PRC_HAS_E ) then
        !--- unpacking packets from East
        !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
        do j = JS, JE
        do i = IE+1, IE+IHALO
           var(i,j) = buf(i-IE,j,1,2)
        enddo
        enddo
       !$acc end kernels
     end if

     if( PRC_HAS_W ) then
        !--- unpacking packets from West
        !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels if(acc_is_present(var)) async
        do j = JS, JE
        do i = IS-IHALO, IS-1
           var(i,j) = buf(i,j,1,1)
        enddo
        enddo
       !$acc end kernels
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackWE_2D

  subroutine unpackNS_3D( KA, IA, IS, IE, JA, JS, JE, &
                          JHALO, &
                          var, buf )
    implicit none
    integer, intent(in) :: KA
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: JHALO
    real(RP), intent(inout) :: var(KA,IA,JA)
    real(RP), intent(in)    :: buf(KA,IA,JHALO,2)

    integer :: i, j, k
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf)

    call PROF_rapstart('COMM_unpack', 3)

    if ( PRC_HAS_S ) then
       !--- unpacking packets from S, SW, and SE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels async
       do j = 1, JS-1
       do i = 1, IA
       do k = 1, KA
          var(k,i,j) = buf(k,i,j,1)
       enddo
       enddo
       enddo
       !$acc end kernels
    else
       if ( PRC_HAS_W ) then
          !--- unpacking packets from SW
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = 1, JS-1
          do i = 1, IS-1
          do k = 1, KA
             var(k,i,j) = buf(k,i,j,1)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from SE
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = 1, JS-1
          do i = IE+1, IA
          do k = 1, KA
             var(k,i,j) = buf(k,i,j,1)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    if ( PRC_HAS_N ) then
       !--- unpacking packets from N, NW, and NE
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels async
       do j = JE+1, JA
       do i = 1, IA
       do k = 1, KA
          var(k,i,j) = buf(k,i,j-JE,2)
       enddo
       enddo
       enddo
       !$acc end kernels
    else
       if ( PRC_HAS_W ) then
          !--- unpacking packets from NW
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = JE+1, JA
          do i = 1, IS-1
          do k = 1, KA
             var(k,i,j) = buf(k,i,j-JE,2)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from NE
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          !$acc kernels async
          do j = JE+1, JA
          do i = IE+1, IA
          do k = 1, KA
             var(k,i,j) = buf(k,i,j-JE,2)
          enddo
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackNS_3D

  subroutine unpackNS_2D( IA, IS, IE, JA, JS, JE, &
                          JHALO, &
                          var, buf )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer, intent(in) :: JHALO
    real(RP), intent(inout) :: var(IA,JA)
    real(RP), intent(in)    :: buf(IA,JHALO,2)

    integer :: i, j
    !---------------------------------------------------------------------------

    !$acc data copy(var) copyin(buf)

    call PROF_rapstart('COMM_unpack', 3)

    if ( PRC_HAS_S ) then
       !--- unpacking packets from S, SW, and SE
       !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels async
       do j = 1, JS-1
       do i = 1, IA
          var(i,j) = buf(i,j,1)
       enddo
       enddo
       !$acc end kernels
    else
       if ( PRC_HAS_W ) then
          !--- unpacking packets from SW
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = 1, JS-1
          do i = 1, IS-1
             var(i,j) = buf(i,j,1)
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from SE
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = 1, JS-1
          do i = IE+1, IA
             var(i,j) = buf(i,j,1)
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    if ( PRC_HAS_N ) then
       !--- unpacking packets from N, NW, and NE
       !$omp parallel do private(i,j) OMP_SCHEDULE_
       !$acc kernels async
       do j = JE+1, JA
       do i = 1, IA
          var(i,j) = buf(i,j-JE,2)
       enddo
       enddo
       !$acc end kernels
    else
       if ( PRC_HAS_W ) then
          !--- unpacking packets from NW
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = JE+1, JA
          do i = 1, IS-1
             var(i,j) = buf(i,j-JE,2)
          enddo
          enddo
          !$acc end kernels
       end if
       if ( PRC_HAS_E ) then
          !--- unpacking packets from NE
          !$omp parallel do private(i,j) OMP_SCHEDULE_
          !$acc kernels async
          do j = JE+1, JA
          do i = IE+1, IA
             var(i,j) = buf(i,j-JE,2)
          enddo
          enddo
          !$acc end kernels
       end if
    end if

    !$acc wait

    call PROF_rapend('COMM_unpack', 3)

    !$acc end data

    return
  end subroutine unpackNS_2D

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

    !$acc data copy(var)

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
       !$acc kernels async
       do j = JE+1, JE+JHALO
       !$omp do
       do i = IS, IE
       do k = 1, KA
          var(k,i,j) = var(k,i,JE)
       enddo
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    !--- copy inner data to HALO(South)
    if ( .NOT. PRC_HAS_S ) then
       !$acc kernels async
       !$acc loop independent
       do j = JS-JHALO, JS-1
       !$omp do
       do i = IS, IE
       do k = 1, KA
          var(k,i,j) = var(k,i,JS)
       enddo
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    if ( .not. PRC_TwoD ) then

       !--- copy inner data to HALO(East)
       if ( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          !$omp do
          do j = JS, JE
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$omp end do nowait
          !$acc end kernels
       end if

       !--- copy inner data to HALO(West)
       if ( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          !$omp do
          do j = JS, JE
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(:,i,j) = var(:,IS,j)
          enddo
          enddo
          !$omp end do nowait
          !$acc end kernels
       end if

       !--- copy inner data to HALO(NorthWest)
       if ( .NOT. PRC_HAS_N .AND. &
            .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,i,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthWest)
       if ( .NOT. PRC_HAS_S .AND. &
            .NOT. PRC_HAS_W ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,i,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
          do k = 1, KA
             var(k,i,j) = var(k,IS,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(NorthEast)
       if ( .NOT. PRC_HAS_N .AND. &
            .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,i,JE)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthEast)
       if ( .NOT. PRC_HAS_S .AND. &
            .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,i,JS)
          enddo
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
          do k = 1, KA
             var(k,i,j) = var(k,IE,j)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

    end if

    !$omp end parallel

    !$acc wait

    !$acc end data

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

    !$acc data copy(var)

    IS    = ginfo(gid)%IS
    IE    = ginfo(gid)%IE
    IHALO = ginfo(gid)%IHALO
    JS    = ginfo(gid)%JS
    JE    = ginfo(gid)%JE
    JHALO = ginfo(gid)%JHALO

    !$omp parallel

    !--- copy inner data to HALO(North)
    if( .NOT. PRC_HAS_N ) then
       !$acc kernels async
       do j = JE+1, JE+JHALO
       !$omp do
       do i = IS, IE
          var(i,j) = var(i,JE)
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    !--- copy inner data to HALO(South)
    if( .NOT. PRC_HAS_S ) then
       !$acc kernels async
       !$acc loop independent
       do j = JS-JHALO, JS-1
       !$omp do
       do i = IS, IE
          var(i,j) = var(i,JS)
       enddo
       !$omp end do nowait
       enddo
       !$acc end kernels
    endif

    if ( .not. PRC_TwoD ) then

       if( .NOT. PRC_HAS_E ) then
          !$omp do
          !$acc kernels async
          do j = JS, JE
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       endif

       if( .NOT. PRC_HAS_W ) then
          !$omp do
          !$acc kernels async
          do j = JS, JE
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       endif

       !--- copy inner data to HALO(NorthWest)
       if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IS-IHALO, IS-1
             var(i,j) = var(i,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthWest)
       if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_W ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IS-IHALO, IS-1
             var(i,j) = var(i,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_W ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          !$acc loop independent
          do i = IS-IHALO, IS-1
             var(i,j) = var(IS,j)
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(NorthEast)
       if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_N ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(i,JE)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JE+1, JE+JHALO
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$acc end kernels
       endif

       !--- copy inner data to HALO(SouthEast)
       if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_S ) then
          !$acc kernels async
          !$acc loop independent
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(i,JS)
          enddo
          enddo
          !$acc end kernels
       elseif( .NOT. PRC_HAS_E ) then
          !$acc kernels async
          do j = JS-JHALO, JS-1
          do i = IE+1, IE+IHALO
             var(i,j) = var(IE,j)
          enddo
          enddo
          !$acc end kernels
       endif

    end if

    !$omp end parallel

    !$acc wait

    !$acc end data

    return
  end subroutine copy_boundary_2D

end module scale_comm_cartesC
