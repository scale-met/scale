!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI Communication module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-11 (R.Yoshida)   [new]
!! @li      2011-11-11 (H.Yashiro)   [mod] Integrate to SCALE-LES ver.3
!! @li      2012-01-10 (Y.Ohno)      [mod] Nonblocking communication (MPI)
!! @li      2012-01-23 (Y.Ohno)      [mod] Self unpacking (MPI)
!! @li      2012-03-12 (H.Yashiro)   [mod] REAL4(MPI)
!! @li      2012-03-12 (Y.Ohno)      [mod] RDMA communication
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-03-27 (H.Yashiro)   [mod] Area/volume weighted total value report
!! @li      2014-06-13 (R.Yoshida)   [mod] gather data from whole processes
!! @li      2014-11-26 (S.Nishizawa) [mod] MPI persistent communication (MPI PC)
!!
!<
#include "inc_openmp.h"
module scale_comm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_MPIstop
  use scale_rm_process, only: &
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
  public :: COMM_vars_init
  public :: COMM_vars8_init
  public :: COMM_vars
  public :: COMM_vars8
  public :: COMM_wait
  public :: COMM_horizontal_mean
  public :: COMM_horizontal_max
  public :: COMM_horizontal_min
  public :: COMM_gather
  public :: COMM_bcast
  public :: COMM_cleanup

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

  interface COMM_horizontal_max
     module procedure COMM_horizontal_max_2D
     module procedure COMM_horizontal_max_3D
  end interface COMM_horizontal_max

  interface COMM_horizontal_min
     module procedure COMM_horizontal_min_2D
     module procedure COMM_horizontal_min_3D
  end interface COMM_horizontal_min

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
  end interface COMM_bcast

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: COMM_datatype   !< datatype of variable
  integer, public :: COMM_world      !< communication world ID
  logical, public :: COMM_FILL_BND = .true. !< switch whether fill boundary data

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private              :: COMM_nreq_max       !< # limit of communication request at once
  integer,  private              :: COMM_vsize_max      !< # limit of communication variables at once
  integer,  private              :: COMM_vsize_max_pc      !< # limit of total communication variables for MPI PC

  logical,  private              :: COMM_IsAllPeriodic  !< periodic boundary condition?

  integer,  private              :: COMM_size2D_NS4     !< 2D data size (W/E    HALO, 4/8-direction comm.)
  integer,  private              :: COMM_size2D_NS8     !< 2D data size (N/S    HALO,   4-direction comm.)
  integer,  private              :: COMM_size2D_WE      !< 2D data size (N/S    HALO,   8-direction comm.)
  integer,  private              :: COMM_size2D_4C      !< 2D data size (corner HALO,   8-direction comm.)

  integer,  private              :: COMM_vars_id = 0    !< id of variables

  logical,  private              :: COMM_USE_MPI_PC = .true.

  real(RP), private, allocatable :: recvpack_W2P(:,:)   !< packing packet (receive, from W)
  real(RP), private, allocatable :: recvpack_E2P(:,:)   !< packing packet (receive, from E)
  real(RP), private, allocatable :: sendpack_P2W(:,:)   !< packing packet (send,    to W  )
  real(RP), private, allocatable :: sendpack_P2E(:,:)   !< packing packet (send,    to E  )
#ifdef DEBUG
  logical,  private, allocatable :: use_packbuf(:)      !< using flag for packing buffer
#endif

  integer,  private, allocatable :: req_cnt (:)         !< request ID of each MPI send/recv
  integer,  private, allocatable :: req_list(:,:)       !< request ID set of each variables
  integer,  private, allocatable :: preq_cnt (:)        !< request ID of each MPI PC
  integer,  private, allocatable :: preq_list(:,:)      !< request ID set of each variables for MPI PC
  integer,  private, allocatable :: pseqid(:)           !< sequential ID of each variables for MPI PC

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine COMM_setup
    use scale_stdio, only: &
       IO_FID_CONF
    use scale_process, only: &
       PRC_LOCAL_COMM_WORLD
    implicit none

    NAMELIST / PARAM_COMM / &
       COMM_vsize_max, &
       COMM_vsize_max_pc, &
       COMM_USE_MPI_PC

    integer :: nreq_NS, nreq_WE, nreq_4C

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[COMM] / Categ[ATMOS-RM COMM] / Origin[SCALElib]'

    COMM_vsize_max = max( 10 + QA*2, 25 )
    COMM_vsize_max_pc = 50 + QA*2

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_COMM. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_COMM)

    nreq_NS  = 2 * JHALO !--- send x JHALO, recv x JHALO
    nreq_WE  = 2         !--- send x 1    , recv x 1
    nreq_4C  = 2 * JHALO !--- send x JHALO, recv x JHALO

    if ( COMM_USE_MPI_PC ) then
       COMM_nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C + 1
    else
       COMM_nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C
    end if

    COMM_size2D_NS4 = IA   * JHALO
    COMM_size2D_NS8 = IMAX
    COMM_size2D_WE  = JMAX * IHALO
    COMM_size2D_4C  =        IHALO

    allocate( recvpack_W2P(COMM_size2D_WE*KA,COMM_vsize_max) )
    allocate( recvpack_E2P(COMM_size2D_WE*KA,COMM_vsize_max) )
    allocate( sendpack_P2W(COMM_size2D_WE*KA,COMM_vsize_max) )
    allocate( sendpack_P2E(COMM_size2D_WE*KA,COMM_vsize_max) )
#ifdef DEBUG
    allocate( use_packbuf(COMM_vsize_max) )
    use_packbuf(:) = .false.
#endif

    allocate( req_cnt (              COMM_vsize_max) )
    allocate( req_list(COMM_nreq_MAX,COMM_vsize_max) )
    req_cnt (:)   = -1
    req_list(:,:) = MPI_REQUEST_NULL

    if ( COMM_USE_MPI_PC ) then
       allocate( preq_cnt (                COMM_vsize_max_pc) )
       allocate( preq_list(COMM_nreq_MAX+1,COMM_vsize_max_pc) )
       preq_cnt (:)   = -1
       preq_list(:,:) = MPI_REQUEST_NULL

       allocate( pseqid(COMM_vsize_max_pc) )
    end if

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
       write(*,*) 'xxx precision is not supportd'
       call PRC_MPIstop
    endif

    COMM_world = PRC_LOCAL_COMM_WORLD

#ifdef _USE_RDMA
    call rdma_setup( COMM_vsize_max_pc, &
                     IA,                &
                     JA,                &
                     KA,                &
                     IHALO,             &
                     JHALO,             &
                     IS,                &
                     IE,                &
                     JS,                &
                     JE,                &
                     PRC_next(PRC_W),   &
                     PRC_next(PRC_N),   &
                     PRC_next(PRC_E),   &
                     PRC_next(PRC_S)    )
#endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Maximum number of vars for one communication: ', &
                                   COMM_vsize_max
    if( IO_L ) write(IO_FID_LOG,*) '*** Data size of var (3D,including halo) [byte] : ', &
                                   RP*KA*IA*JA
    if( IO_L ) write(IO_FID_LOG,*) '*** Data size of halo                    [byte] : ', &
                                   RP*KA*(2*IA*JHALO+2*JMAX*IHALO)
    if( IO_L ) write(IO_FID_LOG,*) '*** Ratio of halo against the whole 3D grid     : ', &
                                   real(2*IA*JHALO+2*JMAX*IHALO) / real(IA*JA)
    if( IO_L ) write(IO_FID_LOG,*) '*** All side is periodic?                       : ', COMM_IsAllPeriodic

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  !> Register variables
  subroutine COMM_vars_init( &
       varname, &
       var,     &
       vid      )
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name
    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID
    !---------------------------------------------------------------------------

    if ( vid > COMM_vsize_max ) then
       write(*,*) 'xxx vid exceeds max', vid, COMM_vsize_max
       call PRC_MPIstop
    end if

    if ( COMM_USE_MPI_PC ) then

       COMM_vars_id = COMM_vars_id + 1
       if ( COMM_vars_id > COMM_vsize_max_pc ) then
          write(*,*) 'xxx number of variable for MPI PC exceeds max', COMM_vars_id, COMM_vsize_max_pc
          call PRC_MPIstop
       end if

#ifdef _USE_RDMA
       call PROF_rapstart('COMM_init_RDMA', 2)
       call set_rdma_variable(var, COMM_vars_id-1)
       call PROF_rapend  ('COMM_init_RDMA', 2)
#else
       call PROF_rapstart('COMM_init_pers', 2)
       call vars_init_mpi_pc(var, COMM_vars_id, vid)
       call PROF_rapend  ('COMM_init_pers', 2)
#endif

       vid = COMM_vars_id + COMM_vsize_max

       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,2A)') '*** [Pers.COMM] Initialize variable : ID = ', vid, &
                                                                                       ', name = ', trim(varname)

    end if

    return
  end subroutine COMM_vars_init

  !-----------------------------------------------------------------------------
  !> Register variables
  subroutine COMM_vars8_init( &
       varname, &
       var,     &
       vid      )
    implicit none

    character(len=*), intent(in)    :: varname    !< variable name
    real(RP),         intent(inout) :: var(:,:,:) !< variable array for register
    integer,          intent(inout) :: vid        !< variable ID
    !---------------------------------------------------------------------------

    if ( vid > COMM_vsize_max ) then
       write(*,*) 'xxx vid exceeds max', vid, COMM_vsize_max
       call PRC_MPIstop
    end if

    if ( COMM_USE_MPI_PC ) then

       COMM_vars_id = COMM_vars_id + 1
       if ( COMM_vars_id > COMM_vsize_max_pc ) then
          write(*,*) 'xxx number of variable for MPI PC exceeds max', COMM_vars_id, COMM_vsize_max_pc
          call PRC_MPIstop
       end if

#ifdef _USE_RDMA
       call PROF_rapstart('COMM_init_RDMA', 2)
       call set_rdma_variable(var, COMM_vars_id-1)
       call PROF_rapend  ('COMM_init_RDMA', 2)
#else
       call PROF_rapstart('COMM_init_pers', 2)
       call vars8_init_mpi_pc(var, COMM_vars_id, vid)
       call PROF_rapend  ('COMM_init_pers', 2)
#endif

       vid = COMM_vars_id + COMM_vsize_max

       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,2A)') '*** [Pers.COMM] Initialize variable : ID = ', vid, &
                                                                                       ', name = ', trim(varname)

    end if

    return
  end subroutine COMM_vars8_init

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_3D(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication
    integer,  intent(in)    :: vid        !< request ID
    !---------------------------------------------------------------------------

    if ( vid > COMM_vsize_max ) then
#ifdef _USE_RDMA
       call PROF_rapstart('COMM_vars_RDMA', 2)
       call rdma_put(vid-COMM_vsize_max-1, 1)
       call PROF_rapend  ('COMM_vars_RDMA', 2)
#else
       call PROF_rapstart('COMM_vars_pers', 2)
       call vars_3D_mpi_pc(var, vid-COMM_vsize_max)
       call PROF_rapend  ('COMM_vars_pers', 2)
#endif
    else
       call PROF_rapstart('COMM_vars', 2)
       call vars_3D_mpi(var, vid)
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_3D(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: vid
    !---------------------------------------------------------------------------

    if ( vid > COMM_vsize_max ) then
#ifdef _USE_RDMA
       call PROF_rapstart('COMM_vars_RDMA', 2)
       call rdma_put8(vid-COMM_vsize_max-1,1)
       call PROF_rapend  ('COMM_vars_RDMA', 2)
#else
       call PROF_rapstart('COMM_vars_pers', 2)
       call vars_3D_mpi_pc(var, vid-COMM_vsize_max)
       call PROF_rapend  ('COMM_vars_pers', 2)
#endif
    else
       call PROF_rapstart('COMM_vars', 2)
       call vars8_3D_mpi(var, vid)
       call PROF_rapend  ('COMM_vars', 2)
    end if

    return
  end subroutine COMM_vars8_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_3D(var, vid, FILL_BND)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid
    logical, intent(in), optional :: FILL_BND

    logical :: FILL_BND_
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    if ( vid > COMM_vsize_max ) then
#ifdef _USE_RDMA
       ! do nothing
#else
       call PROF_rapstart('COMM_wait_pers', 2)
       call wait_3D_mpi_pc(var, vid-COMM_vsize_max)
       call PROF_rapend  ('COMM_wait_pers', 2)
#endif
    else
       call PROF_rapstart('COMM_wait', 2)
       call wait_3D_mpi(var, vid)
       call PROF_rapend  ('COMM_wait', 2)
    end if

    ! copy inner data to boundary
    if ( .NOT. COMM_IsAllPeriodic ) then
       if ( FILL_BND_ ) then
          call copy_boundary_3D(var)
       end if
    end if

    return
  end subroutine COMM_wait_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_2D(var, vid)
    implicit none
    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: vid
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_vars', 2)
    call vars_2D_mpi(var, vid)
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_2D(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: vid
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_vars', 2)
    call vars8_2D_mpi(var, vid)
    call PROF_rapend  ('COMM_vars', 2)

    return
  end subroutine COMM_vars8_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_2D(var, vid, FILL_BND)
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: vid
    logical,  intent(in), optional :: FILL_BND

    logical :: FILL_BND_
    !---------------------------------------------------------------------------

    FILL_BND_ = .true.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    call PROF_rapstart('COMM_wait', 2)
    call wait_2D_mpi(var, vid)
    call PROF_rapend  ('COMM_wait', 2)

    if( .NOT. COMM_IsAllPeriodic ) then
       if ( FILL_BND_ ) then
          call copy_boundary_2D(var)
       end if
    end if

    return
  end subroutine COMM_wait_2D

  !-----------------------------------------------------------------------------
  !> calculate horizontal mean (global total with communication)
  subroutine COMM_horizontal_mean( varmean, var )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    real(RP), intent(out) :: varmean(KA)       !< horizontal mean
    real(RP), intent(in)  :: var    (KA,IA,JA) !< 3D value

    real(RP) :: statval   (KA)
    real(RP) :: statcnt   (KA)
    real(RP) :: allstatval(KA)
    real(RP) :: allstatcnt(KA)
    real(RP) :: zerosw

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    statval(:) = 0.0_RP
    statcnt(:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       if ( abs(var(k,i,j)) < abs(CONST_UNDEF) ) then
          statval(k) = statval(k) + var(k,i,j)
          statcnt(k) = statcnt(k) + 1.D0
       endif
    enddo
    enddo
    enddo

    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval(1),    &
                        allstatval(1), &
                        KA,            &
                        COMM_datatype, &
                        MPI_SUM,       &
                        COMM_world,    &
                        ierr           )
    ! All reduce
    call MPI_Allreduce( statcnt(1),    &
                        allstatcnt(1), &
                        KA,            &
                        COMM_datatype, &
                        MPI_SUM,       &
                        COMM_world,    &
                        ierr           )

    call PROF_rapend  ('COMM_Allreduce', 2)

    do k = 1, KA
       zerosw = 0.5_RP - sign(0.5_RP, allstatcnt(k) - 1.E-12_RP )
       varmean(k) = allstatval(k) / ( allstatcnt(k) + zerosw ) * ( 1.0_RP - zerosw )
       !if( IO_L ) write(IO_FID_LOG,*) k, varmean(k), allstatval(k), allstatcnt(k)
    enddo

    return
  end subroutine COMM_horizontal_mean

  !-----------------------------------------------------------------------------
  !> Get maximum value in horizontal area
  subroutine COMM_horizontal_max_2D( varmax, var )
    implicit none

    real(RP), intent(out) :: varmax     !< horizontal maximum
    real(RP), intent(in)  :: var(IA,JA) !< 2D value

    real(RP) :: statval
    real(RP) :: allstatval

    integer :: ierr
    !---------------------------------------------------------------------------

    statval = maxval(var(IS:IE,JS:JE))

    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval,       &
                        allstatval,    &
                        1,             &
                        COMM_datatype, &
                        MPI_MAX,       &
                        COMM_world,    &
                        ierr           )

    call PROF_rapend  ('COMM_Allreduce', 2)

    varmax = allstatval

    return
  end subroutine COMM_horizontal_max_2D

  !-----------------------------------------------------------------------------
  !> Get maximum value in 3D volume
  subroutine COMM_horizontal_max_3D( varmax, var )
    use scale_const, only: &
       CONST_HUGE
    implicit none

    real(RP), intent(out) :: varmax(KA)       !< horizontal maximum
    real(RP), intent(in)  :: var   (KA,IA,JA) !< 3D value

    real(RP) :: statval   (KA)
    real(RP) :: allstatval(KA)

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    statval(:) = -1.E19_RP
    do k = KS, KE
       statval(k) = maxval(var(k,IS:IE,JS:JE))
    enddo

    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval(1),    &
                        allstatval(1), &
                        KA,            &
                        COMM_datatype, &
                        MPI_MAX,       &
                        COMM_world,    &
                        ierr           )

    call PROF_rapend  ('COMM_Allreduce', 2)

    do k = KS, KE
       varmax(k) = allstatval(k)
    enddo
    varmax(   1:KS-1) = -CONST_HUGE
    varmax(KE+1:KA  ) = -CONST_HUGE

    return
  end subroutine COMM_horizontal_max_3D

  !-----------------------------------------------------------------------------
  !> Get minimum value in horizontal area
  subroutine COMM_horizontal_min_2D( varmin, var )
    implicit none

    real(RP), intent(out) :: varmin     !< horizontal minimum
    real(RP), intent(in)  :: var(IA,JA) !< 2D value

    real(RP) :: statval
    real(RP) :: allstatval

    integer :: ierr
    !---------------------------------------------------------------------------

    statval = minval(var(IS:IE,JS:JE))

    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval,       &
                        allstatval,    &
                        1,             &
                        COMM_datatype, &
                        MPI_MIN,       &
                        COMM_world,    &
                        ierr           )

    call PROF_rapend  ('COMM_Allreduce', 2)

    varmin = allstatval

    return
  end subroutine COMM_horizontal_min_2D

  !-----------------------------------------------------------------------------
  !> Get minimum value in 3D volume
  subroutine COMM_horizontal_min_3D( varmin, var )
    use scale_const, only: &
       CONST_HUGE
    implicit none

    real(RP), intent(out) :: varmin(KA)       !< horizontal minimum
    real(RP), intent(in)  :: var   (KA,IA,JA) !< 3D value

    real(RP) :: statval   (KA)
    real(RP) :: allstatval(KA)

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    statval(:) = -1.E19_RP
    do k = KS, KE
       statval(k) = minval(var(k,IS:IE,JS:JE))
    enddo

    ! [NOTE] always communicate globally
    call PROF_rapstart('COMM_Allreduce', 2)
    ! All reduce
    call MPI_Allreduce( statval(1),    &
                        allstatval(1), &
                        KA,            &
                        COMM_datatype, &
                        MPI_MIN,       &
                        COMM_world,    &
                        ierr           )

    call PROF_rapend  ('COMM_Allreduce', 2)

    do k = KS, KE
       varmin(k) = allstatval(k)
    enddo
    varmin(   1:KS-1) = CONST_HUGE
    varmin(KE+1:KA  ) = CONST_HUGE

    return
  end subroutine COMM_horizontal_min_3D

  !-----------------------------------------------------------------------------
  !> Get data from whole process value in 2D field
  subroutine COMM_gather_2D( recv, send, gIA, gJA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(out) :: recv(:,:) !< receive buffer (gIA,gJA)
    real(RP), intent(in)  :: send(:,:) !< send buffer (gIA,gJA)
    integer,  intent(in)  :: gIA       !< dimension size of x
    integer,  intent(in)  :: gJA       !< dimension size of y

    integer :: sendcounts, recvcounts
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = gIA * gJA
    recvcounts = gIA * gJA

    call MPI_GATHER( send(:,:),      &
                     sendcounts,     &
                     COMM_datatype,  &
                     recv(:,:),      &
                     recvcounts,     &
                     COMM_datatype,  &
                     PRC_masterrank, &
                     COMM_world,     &
                     ierr            )

    return
  end subroutine COMM_gather_2D

  !-----------------------------------------------------------------------------
  !> Get data from whole process value in 3D field
  subroutine COMM_gather_3D( recv, send, gIA, gJA, gKA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(out) :: recv(:,:,:) !< receive buffer(gIA,gJA,gKA)
    real(RP), intent(in)  :: send(:,:,:) !< send buffer   (gIA,gJA,gKA)
    integer,  intent(in)  :: gIA         !< dimension size of x
    integer,  intent(in)  :: gJA         !< dimension size of y
    integer,  intent(in)  :: gKA         !< dimension size of z

    integer :: sendcounts, recvcounts
    integer :: ierr
    !---------------------------------------------------------------------------

    sendcounts = gIA * gJA * gKA
    recvcounts = gIA * gJA * gKA

    call MPI_GATHER( send(:,:,:),    &
                     sendcounts,     &
                     COMM_datatype,  &
                     recv(:,:,:),    &
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
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(inout) :: var  !< broadcast buffer (gIA)

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
  subroutine COMM_bcast_1D( var, gIA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(inout) :: var(:)  !< broadcast buffer (gIA)
    integer,  intent(in)    :: gIA       !< dimension size of x

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = gIA

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
  subroutine COMM_bcast_2D( var, gIA, gJA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(inout) :: var(:,:)  !< broadcast buffer (gIA,gJA)
    integer,  intent(in)    :: gIA       !< dimension size of x
    integer,  intent(in)    :: gJA       !< dimension size of y

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = gIA * gJA

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
  subroutine COMM_bcast_3D( var, gIA, gJA, gKA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(inout) :: var(:,:,:)  !< broadcast buffer(gIA,gJA,gKA)
    integer,  intent(in)    :: gIA         !< dimension size of x
    integer,  intent(in)    :: gJA         !< dimension size of y
    integer,  intent(in)    :: gKA         !< dimension size of z

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = gIA * gJA * gKA

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
  subroutine COMM_bcast_4D( var, gIA, gJA, gKA, gTime )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    real(RP), intent(inout) :: var(:,:,:,:) !< broadcast buffer(gIA,gJA,gKA,gTime)
    integer,  intent(in)    :: gIA          !< dimension size of x
    integer,  intent(in)    :: gJA          !< dimension size of y
    integer,  intent(in)    :: gKA          !< dimension size of z
    integer,  intent(in)    :: gTime        !< dimension size of time

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = gIA * gJA * gKA * gTime
    if ( gIA>0 .AND. gJA>0 .AND. gKA>0 .AND. gTime>0 .AND. &
         counts < 0 ) then
       write(*,*) 'xxx counts overflow'
       call PRC_MPIstop
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
    use scale_process, only: &
       PRC_masterrank
    implicit none

    integer, intent(inout) :: var   !< broadcast buffer (gIA)

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
  subroutine COMM_bcast_INT_1D( var, gIA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    integer, intent(inout) :: var(:)   !< broadcast buffer (gIA)
    integer, intent(in)    :: gIA      !< dimension size of x

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = gIA

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
  subroutine COMM_bcast_INT_2D( var, gIA, gJA )
    use scale_process, only: &
       PRC_masterrank
    implicit none

    integer, intent(inout) :: var(:,:)  !< broadcast buffer (gIA,gJA)
    integer, intent(in)    :: gIA       !< dimension size of x
    integer, intent(in)    :: gJA       !< dimension size of y

    integer :: counts
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_Bcast', 2)

    counts = gIA * gJA

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
    use scale_process, only: &
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

!-------------------------------------------------------------------------------
! private routines
!-------------------------------------------------------------------------------
  subroutine vars_init_mpi_pc(var, vid, seqid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: vid
    integer,  intent(in) :: seqid

    integer :: ireq, tag, ierr
    logical :: flag

    integer :: kd
    integer :: i

    tag = vid * 100
    ireq = 1

    kd = size(var, 1)

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype,           &
                        MPI_PROC_NULL, tag+COMM_nreq_max+1, COMM_world, &
                        preq_list(COMM_nreq_max+1,vid), ierr )

    !--- From 4-Direction HALO communicate
    ! From S
    call MPI_RECV_INIT( var(:,:,JS-JHALO:JS-1), COMM_size2D_NS4*kd, COMM_datatype,      &
                        PRC_next(PRC_S), tag+1, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! From N
    call MPI_RECV_INIT( var(:,:,JE+1:JE+JHALO), COMM_size2D_NS4*kd, COMM_datatype,      &
                        PRC_next(PRC_N), tag+2, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! From E
    call MPI_RECV_INIT( recvpack_E2P(:,seqid), COMM_size2D_WE*kd,  COMM_datatype,      &
                        PRC_next(PRC_E), tag+3, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! From W
    call MPI_RECV_INIT( recvpack_W2P(:,seqid), COMM_size2D_WE*kd,  COMM_datatype,      &
                        PRC_next(PRC_W), tag+4, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1

    !--- To 4-Direction HALO communicate
    ! To W HALO
    call MPI_SEND_INIT( sendpack_P2W(:,seqid), COMM_size2D_WE*kd,  COMM_datatype,      &
                        PRC_next(PRC_W), tag+3, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! To E HALO
    call MPI_SEND_INIT( sendpack_P2E(:,seqid), COMM_size2D_WE*kd,  COMM_datatype,      &
                        PRC_next(PRC_E), tag+4, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! To N HALO
    call MPI_SEND_INIT( var(:,:,JE-JHALO+1:JE), COMM_size2D_NS4*kd, COMM_datatype,      &
                        PRC_next(PRC_N), tag+1, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1
    ! To S HALO
    call MPI_SEND_INIT( var(:,:,JS:JS+JHALO-1), COMM_size2D_NS4*kd, COMM_datatype,      &
                        PRC_next(PRC_S), tag+2, COMM_world, preq_list(ireq,vid), ierr )
    ireq = ireq + 1

    preq_cnt(vid) = ireq - 1
    pseqid(vid) = seqid

    ! to finish initial processes of MPIa
    do i = 1, 32
       call MPI_TESTALL( preq_cnt(vid), preq_list(1:preq_cnt(vid),vid), &
                         flag, MPI_STATUSES_IGNORE, ierr )
    enddo

    return
  end subroutine vars_init_mpi_pc

  subroutine vars8_init_mpi_pc(var, vid, seqid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in) :: vid
    integer,  intent(in) :: seqid

    integer :: ireq, tag, tagc
    integer :: ierr
    integer :: kd
    logical :: flag

    integer :: i, j

    kd = size(var, 1)

    tag  = vid * 100
    ireq = 1

    ! register whole array to inner table of MPI and/or lower library
    ! otherwise a lot of sub small segments would be registered
    call MPI_SEND_INIT( var(:,:,:), size(var), COMM_datatype,         &
                        MPI_PROC_NULL, tag+COMM_nreq_max+1, COMM_world, &
                        preq_list(COMM_nreq_max+1,vid), ierr )


    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       ! From SE
       tagc = 0
       do j = JS-JHALO, JS-1
          call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_SE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From SW
       tagc = 10
       do j = JS-JHALO, JS-1
          call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_SW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From NE
       tagc = 20
       do j = JE+1, JE+JHALO
          call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_NE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From NW
       tagc = 30
       do j = JE+1, JE+JHALO
          call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_NW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From S
       tagc = 40
       do j = JS-JHALO, JS-1
          call MPI_RECV_INIT( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                          PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JE+JHALO
          call MPI_RECV_INIT( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                          PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From E
       tagc = 60
       call MPI_RECV_INIT( recvpack_E2P(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From W
       tagc = 70
       call MPI_RECV_INIT( recvpack_W2P(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
       ireq = ireq + 1

       !--- To 8-Direction HALO communicate
       ! To W HALO
       tagc = 60
       call MPI_SEND_INIT( sendpack_P2W(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To E HALO
       tagc = 70
       call MPI_SEND_INIT( sendpack_P2E(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_SEND_INIT( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                          PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_SEND_INIT( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                          PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To NW HALO
       tagc = 0
       do j = JE-JHALO+1, JE
          call MPI_SEND_INIT( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_NW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To NE HALO
       tagc = 10
       do j = JE-JHALO+1, JE
          call MPI_SEND_INIT( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_NE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To SW HALO
       tagc = 20
       do j = JS, JS+JHALO-1
          call MPI_SEND_INIT( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_SW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To SE HALO
       tagc = 30
       do j = JS, JS+JHALO-1
          call MPI_SEND_INIT( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_SE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

    else ! non-periodic condition

       !--- From 8-Direction HALO communicate
       ! From SE
       if ( PRC_HAS_S .AND. PRC_HAS_E ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From SW
       if ( PRC_HAS_S .AND. PRC_HAS_W ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From NE
       if ( PRC_HAS_N .AND. PRC_HAS_E ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From NW
       if ( PRC_HAS_N .AND. PRC_HAS_W ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = JS-JHALO, JS-1
             call MPI_RECV_INIT( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JE+JHALO
             call MPI_RECV_INIT( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From E
       if ( PRC_HAS_E ) then
          tagc = 60
          call MPI_RECV_INIT( recvpack_E2P(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
          tagc = 70
          call MPI_RECV_INIT( recvpack_W2P(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

       !--- To 8-Direction HALO communicate
       ! To W HALO
       if ( PRC_HAS_W ) then
          tagc = 60
          call MPI_SEND_INIT( sendpack_P2W(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
          tagc = 70
          call MPI_SEND_INIT( sendpack_P2E(:,seqid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                             PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                             PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To NW HALO
       if ( PRC_HAS_N .AND. PRC_HAS_W ) then
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_SEND_INIT( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To NE HALO
       if ( PRC_HAS_N .AND. PRC_HAS_E ) then
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_SEND_INIT( var(1,IE+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_N), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To SW HALO
       if ( PRC_HAS_S .AND. PRC_HAS_W ) then
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IS-IHALO,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_SEND_INIT( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_W), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To SE HALO
       if ( PRC_HAS_S .AND. PRC_HAS_E ) then
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_SEND_INIT( var(1,IE+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_S), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_SEND_INIT( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_E), tag+tagc, COMM_world, preq_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

    endif

    preq_cnt(vid) = ireq - 1
    pseqid(vid) = seqid

    ! to finish initial processes of MPIa
    do i = 1, 32
       call MPI_TESTALL( preq_cnt(vid), preq_list(1:preq_cnt(vid),vid), &
                         flag, MPI_STATUSES_IGNORE, ierr )
    enddo

    return
  end subroutine vars8_init_mpi_pc

  subroutine vars_3D_mpi(var, vid)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< atmospheric 3D variable to communication
    integer,  intent(in)    :: vid        !< request ID


    integer :: ireq, tag

    integer :: kd
    integer :: ierr
    !---------------------------------------------------------------------------

    tag  = vid * 100
    ireq = 1

    kd = size(var, 1)

#ifdef DEBUG
    if ( use_packbuf(vid) ) then
       write(*,*) 'packing buffer is already used', vid
       call PRC_MPIstop
    end if
    use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 4-Direction HALO communicate
       ! From S
       call MPI_IRECV( var(:,:,JS-JHALO:JS-1), COMM_size2D_NS4*kd, COMM_datatype,      &
                       PRC_next(PRC_S), tag+1, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From N
       call MPI_IRECV( var(:,:,JE+1:JE+JHALO), COMM_size2D_NS4*kd, COMM_datatype,      &
                       PRC_next(PRC_N), tag+2, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From E
       call MPI_IRECV( recvpack_E2P(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                       PRC_next(PRC_E), tag+3, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From W
       call MPI_IRECV( recvpack_W2P(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                       PRC_next(PRC_W), tag+4, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1

       call pack_3D(var, vid)

       !--- To 4-Direction HALO communicate
       ! To W HALO
       call MPI_ISEND( sendpack_P2W(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                       PRC_next(PRC_W), tag+3, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To E HALO
       call MPI_ISEND( sendpack_P2E(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                       PRC_next(PRC_E), tag+4, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To N HALO
       call MPI_ISEND( var(:,:,JE-JHALO+1:JE), COMM_size2D_NS4*kd, COMM_datatype,      &
                       PRC_next(PRC_N), tag+1, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To S HALO
       call MPI_ISEND( var(:,:,JS:JS+JHALO-1), COMM_size2D_NS4*kd, COMM_datatype,      &
                       PRC_next(PRC_S), tag+2, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1

    else ! non-periodic condition

       !--- From 4-Direction HALO communicate
       ! From S
       if ( PRC_HAS_S ) then
          call MPI_IRECV( var(:,:,JS-JHALO:JS-1), COMM_size2D_NS4*kd, COMM_datatype,      &
                          PRC_next(PRC_S), tag+1, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From N
       if ( PRC_HAS_N ) then
          call MPI_IRECV( var(:,:,JE+1:JE+JHALO), COMM_size2D_NS4*kd, COMM_datatype,      &
                          PRC_next(PRC_N), tag+2, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From E
       if ( PRC_HAS_E ) then
          call MPI_IRECV( recvpack_E2P(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                          PRC_next(PRC_E), tag+3, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
          call MPI_IRECV( recvpack_W2P(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                          PRC_next(PRC_W), tag+4, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

       call pack_3D(var, vid)

       !--- To 4-Direction HALO communicate
       ! To W HALO
       if ( PRC_HAS_W ) then
          call MPI_ISEND( sendpack_P2W(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                          PRC_next(PRC_W), tag+3, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
          call MPI_ISEND( sendpack_P2E(:,vid),    COMM_size2D_WE*kd,  COMM_datatype,      &
                          PRC_next(PRC_E), tag+4, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To N HALO
       if ( PRC_HAS_N ) then
          call MPI_ISEND( var(:,:,JE-JHALO+1:JE), COMM_size2D_NS4*kd, COMM_datatype,      &
                          PRC_next(PRC_N), tag+1, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          call MPI_ISEND( var(:,:,JS:JS+JHALO-1), COMM_size2D_NS4*kd, COMM_datatype,      &
                          PRC_next(PRC_S), tag+2, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

    endif

    req_cnt(vid) = ireq - 1

    return
  end subroutine vars_3D_mpi

  subroutine vars8_3D_mpi(var, vid)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: vid

    integer :: ireq, tag, tagc

    integer :: kd

    integer :: ierr
    integer :: j
    !---------------------------------------------------------------------------

    tag  = vid * 100
    ireq = 1

    kd = size(var, 1)

#ifdef DEBUG
    if ( use_packbuf(vid) ) then
       write(*,*) 'packing buffer is already used', vid
       call PRC_MPIstop
    end if
    use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       ! From SE
       tagc = 0
       do j = JS-JHALO, JS-1
          call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From SW
       tagc = 10
       do j = JS-JHALO, JS-1
          call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From NE
       tagc = 20
       do j = JE+1, JE+JHALO
          call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From NW
       tagc = 30
       do j = JE+1, JE+JHALO
          call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                          PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From S
       tagc = 40
       do j = JS-JHALO, JS-1
          call MPI_IRECV( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                          PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JE+JHALO
          call MPI_IRECV( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                          PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From E
       tagc = 60
       call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From W
       tagc = 70
       call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1

       call pack_3D(var, vid)


       !--- To 8-Direction HALO communicate
       ! To W HALO
       tagc = 60
       call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To E HALO
       tagc = 70
       call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                       PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                          PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                          PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To NW HALO
       tagc = 0
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To NE HALO
       tagc = 10
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To SW HALO
       tagc = 20
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To SE HALO
       tagc = 30
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                          PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

    else ! non-periodic condition

       !--- From 8-Direction HALO communicate
       ! From SE
       if ( PRC_HAS_S .AND. PRC_HAS_E ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From SW
       if ( PRC_HAS_S .AND. PRC_HAS_W ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From NE
       if ( PRC_HAS_N .AND. PRC_HAS_E ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IE+1,j),     COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From NW
       if ( PRC_HAS_N .AND. PRC_HAS_W ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size2D_4C*kd, COMM_datatype,                &
                             PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From S
       if ( PRC_HAS_S ) then
          tagc = 40
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_HAS_N ) then
          tagc = 50
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IS,j),     COMM_size2D_NS8*kd, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From E
       if ( PRC_HAS_E ) then
          tagc = 60
          call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_HAS_W ) then
          tagc = 70
          call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

       call pack_3D(var, vid)

       !--- To 8-Direction HALO communicate
       ! To W HALO
       if ( PRC_HAS_W ) then
          tagc = 60
          call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_HAS_E ) then
          tagc = 70
          call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE*kd, COMM_datatype,             &
                          PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To N HALO
       if ( PRC_HAS_N ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_HAS_S ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), COMM_size2D_NS8*kd, COMM_datatype,                    &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To NW HALO
       if ( PRC_HAS_N .AND. PRC_HAS_W ) then
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,1,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_ISEND( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To NE HALO
       if ( PRC_HAS_N .AND. PRC_HAS_E ) then
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_N ) then
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IE+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To SW HALO
       if ( PRC_HAS_S .AND. PRC_HAS_W ) then
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,1,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_W ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_ISEND( var(1,IS,j),         COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To SE HALO
       if ( PRC_HAS_S .AND. PRC_HAS_E ) then
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_S ) then
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IE+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       else if ( PRC_HAS_E ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size2D_4C*kd, COMM_datatype,              &
                             PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

    endif

    req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_3D_mpi

  subroutine vars_2D_mpi(var, vid)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)    :: vid

    integer :: ireq, tag
    integer :: ierr
    !---------------------------------------------------------------------------

    tag = vid * 100
    ireq = 1

#ifdef DEBUG
    if ( use_packbuf(vid) ) then
       write(*,*) 'packing buffer is already used', vid
       call PRC_MPIstop
    end if
    use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 4-Direction HALO communicate
        ! From S
        call MPI_IRECV( var(:,JS-JHALO:JS-1), COMM_size2D_NS4,        &
                        COMM_datatype, PRC_next(PRC_S), tag+1, &
                        COMM_world, req_list(ireq,vid), ierr    )
        ireq = ireq + 1

        ! From N
        call MPI_IRECV( var(:,JE+1:JE+JHALO), COMM_size2D_NS4,        &
                        COMM_datatype, PRC_next(PRC_N), tag+2, &
                        COMM_world, req_list(ireq,vid), ierr    )
        ireq = ireq + 1

        ! From E
        call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,       &
                        COMM_datatype, PRC_next(PRC_E), tag+3, &
                        COMM_world, req_list(ireq,vid), ierr )
        ireq = ireq + 1

        ! From W
        call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,       &
                        COMM_datatype, PRC_next(PRC_W), tag+4, &
                        COMM_world, req_list(ireq,vid), ierr )
        ireq = ireq + 1

        call pack_2D(var, vid)

        ! To W HALO communicate
        call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,       &
                        COMM_datatype, PRC_next(PRC_W), tag+3, &
                        COMM_world, req_list(ireq,vid), ierr )
        ireq = ireq + 1

        ! To E HALO communicate
        call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,       &
                        COMM_datatype, PRC_next(PRC_E), tag+4, &
                        COMM_world, req_list(ireq,vid), ierr )
        ireq = ireq + 1

        ! To N HALO communicate
        call MPI_ISEND( var(:,JE-JHALO+1:JE), COMM_size2D_NS4,        &
                        COMM_datatype, PRC_next(PRC_N), tag+1, &
                        COMM_world, req_list(ireq,vid), ierr    )
        ireq = ireq + 1

        ! To S HALO communicate
        call MPI_ISEND( var(:,JS:JS+JHALO-1), COMM_size2D_NS4,        &
                        COMM_datatype, PRC_next(PRC_S), tag+2, &
                        COMM_world, req_list(ireq,vid), ierr    )
        ireq = ireq + 1

    else
    !--- non-periodic condition
        !--- From 4-Direction HALO communicate
        ! From S
        if ( PRC_HAS_S ) then
            call MPI_IRECV( var(:,JS-JHALO:JS-1), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_S), tag+1, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

        ! From N
        if ( PRC_HAS_N ) then
            call MPI_IRECV( var(:,JE+1:JE+JHALO), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_N), tag+2, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

        ! From E
        if ( PRC_HAS_E ) then
            call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_E), tag+3, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! From W
        if ( PRC_HAS_W ) then
            call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_W), tag+4, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        call pack_2D(var, vid)

        ! To W HALO communicate
        if ( PRC_HAS_W ) then
            call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_W), tag+3, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! To E HALO communicate
        if ( PRC_HAS_E ) then
            call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_E), tag+4, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! To N HALO communicate
        if ( PRC_HAS_N ) then
            call MPI_ISEND( var(:,JE-JHALO+1:JE), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_N), tag+1, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

        ! To S HALO communicate
        if ( PRC_HAS_S ) then
            call MPI_ISEND( var(:,JS:JS+JHALO-1), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_S), tag+2, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

    endif

    req_cnt(vid) = ireq - 1

    return
  end subroutine vars_2D_mpi

  subroutine vars8_2D_mpi(var, vid)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: vid

    integer :: ireq, tag, tagc

    integer :: ierr
    integer :: j
    !---------------------------------------------------------------------------

    tag   = vid * 100
    ireq = 1

#ifdef DEBUG
    if ( use_packbuf(vid) ) then
       write(*,*) 'packing buffer is already used', vid
       call PRC_MPIstop
    end if
    use_packbuf(vid) = .true.
#endif

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 8-Direction HALO communicate
        ! From SE
        tagc = 0
        do j = JS-JHALO, JS-1
            call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                            COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From SW
        tagc = 10
        do j = JS-JHALO, JS-1
            call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                            COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From NE
        tagc = 20
        do j = JE+1, JE+JHALO
            call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                            COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From NW
        tagc = 30
        do j = JE+1, JE+JHALO
            call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                            COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From S
        tagc = 40
        do j = JS-JHALO, JS-1
            call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr       )
             ireq = ireq + 1
             tagc = tagc + 1
        enddo
        ! From N
        tagc = 50
        do j = JE+1, JE+JHALO
            call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr       )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From E
        call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,           &
                        COMM_datatype, PRC_next(PRC_E), tag+60, &
                        COMM_world, req_list(ireq,vid), ierr     )
        ireq = ireq + 1
        ! From W
        call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,           &
                        COMM_datatype, PRC_next(PRC_W), tag+70, &
                        COMM_world, req_list(ireq,vid), ierr     )
        ireq = ireq + 1

        call pack_2D(var, vid)

        ! To W HALO communicate
        call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,            &
                        COMM_datatype, PRC_next(PRC_W), tag+60, &
                        COMM_world, req_list(ireq,vid), ierr     )
        ireq = ireq + 1

        ! To E HALO communicate
        call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,           &
                        COMM_datatype, PRC_next(PRC_E), tag+70, &
                        COMM_world, req_list(ireq,vid), ierr     )
        ireq = ireq + 1

        ! To N HALO communicate
        tagc = 40
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr       )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To S HALO communicate
        tagc = 50
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr       )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To NW HALO communicate
        tagc = 0
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                            COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To NE HALO communicate
        tagc = 10
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                            COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To SW HALO communicate
        tagc = 20
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                            COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To SE HALO communicate
        tagc = 30
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                            COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                            COMM_world, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
    else
    !--- non-periodic condition
        !--- From 8-Direction HALO communicate
        ! From SE
        if ( PRC_HAS_S .AND. PRC_HAS_E ) then
            tagc = 0
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_S ) then
            tagc = 0
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_E ) then
            tagc = 0
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_E), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From SW
        if ( PRC_HAS_S .AND. PRC_HAS_W ) then
            tagc = 10
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_S ) then
            tagc = 10
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_W ) then
            tagc = 10
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_W), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From NE
        if ( PRC_HAS_N .AND. PRC_HAS_E ) then
            tagc = 20
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_N ) then
            tagc = 20
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_E ) then
            tagc = 20
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_E), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From NW
        if ( PRC_HAS_N .AND. PRC_HAS_W ) then
            tagc = 30
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_N ) then
            tagc = 30
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_W ) then
            tagc = 30
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_W), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From S
        if ( PRC_HAS_S ) then
            tagc = 40
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr       )
                 ireq = ireq + 1
                 tagc = tagc + 1
            enddo
        endif

        ! From N
        if ( PRC_HAS_N ) then
            tagc = 50
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr       )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From E
        if ( PRC_HAS_E ) then
            call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,             &
                            COMM_datatype, PRC_next(PRC_E), tag+60, &
                            COMM_world, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        ! From W
        if ( PRC_HAS_W ) then
            call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,           &
                            COMM_datatype, PRC_next(PRC_W), tag+70, &
                            COMM_world, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        call pack_2D(var, vid)

        ! To W HALO communicate
        if ( PRC_HAS_W ) then
            call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,           &
                            COMM_datatype, PRC_next(PRC_W), tag+60, &
                            COMM_world, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        ! To E HALO communicate
        if ( PRC_HAS_E ) then
            call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,           &
                            COMM_datatype, PRC_next(PRC_E), tag+70, &
                            COMM_world, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        ! To N HALO communicate
        if ( PRC_HAS_N ) then
            tagc = 40
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr       )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To S HALO communicate
        if ( PRC_HAS_S ) then
            tagc = 50
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr       )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To NW HALO communicate
        if ( PRC_HAS_N .AND. PRC_HAS_W ) then
            tagc = 0
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_N ) then
            tagc = 10
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_W ) then
            tagc = 20
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_W), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To NE HALO communicate
        if ( PRC_HAS_N .AND. PRC_HAS_E ) then
            tagc = 10
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_N ) then
            tagc = 0
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_E ) then
            tagc = 30
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_E), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To SW HALO communicate
        if ( PRC_HAS_S .AND. PRC_HAS_W ) then
            tagc = 20
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_S ) then
            tagc = 30
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_W ) then
            tagc = 0
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_W), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To SE HALO communicate
        if ( PRC_HAS_S .AND. PRC_HAS_E ) then
            tagc = 30
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_S ) then
            tagc = 20
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
         else if ( PRC_HAS_E ) then
            tagc = 10
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_E), tag+tagc, &
                                COMM_world, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

    endif

    req_cnt(vid) = ireq - 1

    return
  end subroutine vars8_2D_mpi

  subroutine vars_3D_mpi_pc(var, vid)
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid
    integer :: ierr
    !---------------------------------------------------------------------------

#ifdef DEBUG
    if ( use_packbuf(pseqid(vid)) ) then
       write(*,*) 'packing buffer is already used', vid, pseqid(vid)
       call PRC_MPIstop
    end if
    use_packbuf(pseqid(vid)) = .true.
#endif

    call pack_3D(var, pseqid(vid))

    call MPI_STARTALL(preq_cnt(vid), preq_list(1:preq_cnt(vid),vid), ierr)

    return
  end subroutine vars_3D_mpi_pc

  subroutine wait_3D_mpi(var, vid)
    implicit none
    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call MPI_WAITALL( req_cnt (vid),                &
                      req_list(1:req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,          &
                      ierr                          )
    call unpack_3D(var, vid)

#ifdef DEBUG
    use_packbuf(vid) = .false.
#endif

    return
  end subroutine wait_3D_mpi

  subroutine wait_2D_mpi(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: vid

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- wait packets
    call MPI_WAITALL( req_cnt(vid), &
                      req_list(1:req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE, &
                      ierr )
    call unpack_2D(var, vid)

#ifdef DEBUG
    use_packbuf(vid) = .false.
#endif

    return
  end subroutine wait_2D_mpi

  subroutine wait_3D_mpi_pc(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: vid

    integer :: ierr

    !--- wait packets
    call MPI_WAITALL( preq_cnt (vid),                &
                      preq_list(1:preq_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,          &
                      ierr                          )
    call unpack_3D(var, pseqid(vid))

#ifdef DEBUG
    use_packbuf(pseqid(vid)) = .false.
#endif

    return
  end subroutine wait_3D_mpi_pc

  subroutine pack_3D(var, vid)
    implicit none

    real(RP), intent(in) :: var(:,:,:)
    integer,  intent(in) :: vid

    integer :: kd
    integer :: k, i, j, n

    kd = size(var, 1)

    call PROF_rapstart('COMM_pack', 3)

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- packing packets to West
!OCL NORECURRENCE(sendpack_P2W)
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IS+IHALO-1
       do k = 1, kd
          n = (j-JS) * kd * IHALO &
            + (i-IS) * kd         &
            + k
          sendpack_P2W(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo
       !--- packing packets to East
!OCL NORECURRENCE(sendpack_P2E)
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE-IHALO+1, IE
       do k = 1, kd
          n = (j-JS)         * kd * IHALO &
            + (i-IE+IHALO-1) * kd         &
            + k
          sendpack_P2E(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo

    else

       if ( PRC_HAS_W ) then
          !--- packing packets to West
!OCL NORECURRENCE(sendpack_P2W)
          !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IHALO,kd,var,sendpack_P2W,vid) 
          do j = JS, JE
          do i = IS, IS+IHALO-1
          do k = 1, kd
             n = (j-JS) * kd * IHALO &
               + (i-IS) * kd         &
               + k
             sendpack_P2W(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif
       if ( PRC_HAS_E ) then
          !--- packing packets to East
!OCL NORECURRENCE(sendpack_P2E)
          !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IE,IHALO,kd,var,sendpack_P2E,vid) 
          do j = JS, JE
          do i = IE-IHALO+1, IE
          do k = 1, kd
             n = (j-JS)         * kd * IHALO &
               + (i-IE+IHALO-1) * kd         &
               + k
             sendpack_P2E(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif

    end if

    call PROF_rapend('COMM_pack', 3)

    return
  end subroutine pack_3D

  subroutine pack_2D(var, vid)
    implicit none

    real(RP), intent(in) :: var(:,:)
    integer,  intent(in) :: vid

    integer :: i, j, n

    call PROF_rapstart('COMM_pack', 3)

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- To 4-Direction HALO communicate
       !--- packing packets to West
!OCL NORECURRENCE(sendpack_P2W)
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IS+IHALO-1
          n = (j-JS) * IHALO &
            + (i-IS) + 1
          sendpack_P2W(n,vid) = var(i,j)
       enddo
       enddo

       !--- packing packets to East
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
!OCL NORECURRENCE(sendpack_P2E)
       do j = JS, JE
       do i = IE-IHALO+1, IE
          n = (j-JS)         * IHALO &
            + (i-IE+IHALO-1) + 1
          sendpack_P2E(n,vid) = var(i,j)
       enddo
       enddo

    else

       !--- To 4-Direction HALO communicate
       !--- packing packets to West
       if ( PRC_HAS_W ) then
!OCL NORECURRENCE(sendpack_P2W)
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IS+IHALO-1
             n = (j-JS) * IHALO &
               + (i-IS) + 1
             sendpack_P2W(n,vid) = var(i,j)
          enddo
          enddo
       endif

       !--- packing packets to East
       if ( PRC_HAS_E ) then
!OCL NORECURRENCE(sendpack_P2E)
       !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IE-IHALO+1, IE
             n = (j-JS)         * IHALO &
               + (i-IE+IHALO-1) + 1
             sendpack_P2E(n,vid) = var(i,j)
          enddo
          enddo
       endif

    end if

    call PROF_rapend('COMM_pack', 3)

    return
  end subroutine pack_2D

  subroutine unpack_3D(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: vid

    integer :: kd
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    kd = size(var, 1)

    call PROF_rapstart('COMM_unpack', 3)

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- unpacking packets from East
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE+1, IE+IHALO
       do k = 1, kd
          n = (j-JS)   * kd * IHALO &
            + (i-IE-1) * kd         &
            + k
          var(k,i,j) = recvpack_E2P(n,vid)
       enddo
       enddo
       enddo

       !--- unpacking packets from West
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS-IHALO, IS-1
       do k = 1, kd
          n = (j-JS)       * kd * IHALO &
            + (i-IS+IHALO) * kd         &
            + k
          var(k,i,j) = recvpack_W2P(n,vid)
       enddo
       enddo
       enddo

    else ! non-periodic condition

        if ( PRC_HAS_E ) then
           !--- unpacking packets from East
           !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
           !$omp shared(JS,JE,IE,IHALO,kd,var,recvpack_E2P,vid) 
           do j = JS, JE
           do i = IE+1, IE+IHALO
           do k = 1, kd
              n = (j-JS)   * kd * IHALO &
                + (i-IE-1) * kd         &
                + k
              var(k,i,j) = recvpack_E2P(n,vid)
           enddo
           enddo
           enddo
        endif

        if ( PRC_HAS_W ) then
           !--- unpacking packets from West
           !$omp parallel do default(none) private(i,j,k,n) OMP_SCHEDULE_ collapse(2) &
           !$omp shared(JS,JE,IS,IHALO,kd,var,recvpack_W2P,vid) 
           do j = JS, JE
           do i = IS-IHALO, IS-1
           do k = 1, kd
              n = (j-JS)       * kd * IHALO &
                + (i-IS+IHALO) * kd         &
                + k
              var(k,i,j) = recvpack_W2P(n,vid)
           enddo
           enddo
           enddo
        endif

    end if

    call PROF_rapend('COMM_unpack', 3)

    return
  end subroutine unpack_3D

  subroutine unpack_2D(var, vid)
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer,  intent(in)    :: vid

    integer :: i, j, n
    !---------------------------------------------------------------------------

    call PROF_rapstart('COMM_unpack', 3)

    if( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- unpacking packets from East
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IE+1, IE+IHALO
           n = (j-JS)   * IHALO &
             + (i-IE-1) + 1
           var(i,j) = recvpack_E2P(n,vid)
        enddo
        enddo

        !--- unpacking packets from West
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IS-IHALO, IS-1
           n = (j-JS)       * IHALO &
             + (i-IS+IHALO) + 1
           var(i,j) = recvpack_W2P(n,vid)
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
             var(i,j) = recvpack_E2P(n,vid)
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
             var(i,j) = recvpack_W2P(n,vid)
          enddo
          enddo
       end if

    end if

    call PROF_rapend('COMM_unpack', 3)

    return
  end subroutine unpack_2D

  subroutine copy_boundary_3D(var)
    implicit none

    real(RP), intent(inout) :: var(:,:,:)

    integer :: i, j
    !---------------------------------------------------------------------------

    !--- copy inner data to HALO(North)
    if ( .NOT. PRC_HAS_N ) then
       do j = JE+1, JE+JHALO
       do i = IS, IE
          var(:,i,j) = var(:,i,JE)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(South)
    if ( .NOT. PRC_HAS_S ) then
       do j = JS-JHALO, JS-1
       do i = IS, IE
          var(:,i,j) = var(:,i,JS)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(East)
    if ( .NOT. PRC_HAS_E ) then
       do j = JS, JE
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,IE,j)
       enddo
       enddo
    end if

    !--- copy inner data to HALO(West)
    if ( .NOT. PRC_HAS_W ) then
       do j = JS, JE
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,IS,j)
       enddo
       enddo
    end if

    !--- copy inner data to HALO(NorthWest)
    if ( .NOT. PRC_HAS_N .AND. &
         .NOT. PRC_HAS_W ) then
       do j = JE+1, JE+JHALO
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,IS,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_N ) then
       do j = JE+1, JE+JHALO
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,i,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_W ) then
       do j = JE+1, JE+JHALO
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,IS,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(SouthWest)
    if ( .NOT. PRC_HAS_S .AND. &
         .NOT. PRC_HAS_W ) then
       do j = JS-IHALO, JS-1
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,IS,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_S ) then
       do j = JS-IHALO, JS-1
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,i,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_W ) then
       do j = JS-IHALO, JS-1
       do i = IS-IHALO, IS-1
          var(:,i,j) = var(:,IS,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(NorthEast)
    if ( .NOT. PRC_HAS_N .AND. &
         .NOT. PRC_HAS_E ) then
       do j = JE+1, JE+JHALO
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,IE,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_N ) then
       do j = JE+1, JE+JHALO
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,i,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_E ) then
       do j = JE+1, JE+JHALO
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,IE,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(SouthEast)
    if ( .NOT. PRC_HAS_S .AND. &
         .NOT. PRC_HAS_E ) then
       do j = JS-IHALO, JS-1
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,IE,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_S ) then
       do j = JS-IHALO, JS-1
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,i,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_E ) then
       do j = JS-IHALO, JS-1
       do i = IE+1, IE+IHALO
          var(:,i,j) = var(:,IE,j)
       enddo
       enddo
    endif

    return
  end subroutine copy_boundary_3D

  subroutine copy_boundary_2D(var)
    implicit none

    real(RP), intent(inout) :: var(:,:)

    integer :: i, j
    !---------------------------------------------------------------------------
    !--- copy inner data to HALO(North)
    if( .NOT. PRC_HAS_N ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IS, IE
          var(i,j) = var(i,JE)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(South)
    if( .NOT. PRC_HAS_S ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-JHALO, JS-1
       do i = IS, IE
          var(i,j) = var(i,JS)
       enddo
       enddo
    endif

    if( .NOT. PRC_HAS_E ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE+1, IE+IHALO
          var(i,j) = var(IE,j)
       enddo
       enddo
    endif

    if( .NOT. PRC_HAS_W ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS-IHALO, IS-1
          var(i,j) = var(IS,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(NorthWest)
    if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_W ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IS-IHALO, IS-1
          var(i,j) = var(IS,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_N ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IS-IHALO, IS-1
          var(i,j) = var(i,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_W ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IS-IHALO, IS-1
          var(i,j) = var(IS,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(SouthWest)
    if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_W ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-IHALO, JS-1
       do i = IS-IHALO, IS-1
          var(i,j) = var(IS,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_S ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-IHALO, JS-1
       do i = IS-IHALO, IS-1
          var(i,j) = var(i,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_W ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-IHALO, JS-1
       do i = IS-IHALO, IS-1
          var(i,j) = var(IS,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(NorthEast)
    if( .NOT. PRC_HAS_N .AND. .NOT. PRC_HAS_E ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IE+1, IE+IHALO
          var(i,j) = var(IE,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_N ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IE+1, IE+IHALO
          var(i,j) = var(i,JE)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_E ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JE+1, JE+JHALO
       do i = IE+1, IE+IHALO
          var(i,j) = var(IE,j)
       enddo
       enddo
    endif

    !--- copy inner data to HALO(SouthEast)
    if( .NOT. PRC_HAS_S .AND. .NOT. PRC_HAS_E ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-IHALO, JS-1
       do i = IE+1, IE+IHALO
          var(i,j) = var(IE,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_S ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-IHALO, JS-1
       do i = IE+1, IE+IHALO
          var(i,j) = var(i,JS)
       enddo
       enddo
    elseif( .NOT. PRC_HAS_E ) then
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS-IHALO, JS-1
       do i = IE+1, IE+IHALO
          var(i,j) = var(IE,j)
       enddo
       enddo
    endif

    return
  end subroutine copy_boundary_2D

  subroutine COMM_cleanup
    use mpi
    implicit none

    NAMELIST / PARAM_COMM / &
       COMM_vsize_max_pc, &
       COMM_USE_MPI_PC

    integer :: i, j, ierr
    !---------------------------------------------------------------------------

    deallocate( recvpack_W2P )
    deallocate( recvpack_E2P )
    deallocate( sendpack_P2W )
    deallocate( sendpack_P2E )
#ifdef DEBUG
    deallocate( use_packbuf )
#endif

    deallocate( req_cnt )
    deallocate( req_list )

    if ( COMM_USE_MPI_PC ) then
       do j=1, COMM_vsize_max_pc
          do i=1, COMM_nreq_MAX+1
             if (preq_list(i,j) .NE. MPI_REQUEST_NULL) &
                 call MPI_REQUEST_FREE(preq_list(i,j), ierr)
          enddo
       enddo
       deallocate( preq_cnt )
       deallocate( preq_list )
       deallocate( pseqid )
    end if
  end subroutine COMM_cleanup

end module scale_comm
