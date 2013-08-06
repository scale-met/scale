!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI Communication module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-11 (R.Yoshida)  [new]
!! @li      2011-11-11 (H.Yashiro)  [mod] Integrate to SCALE3
!! @li      2012-01-10 (Y.Ohno)     [mod] Nonblocking communication (MPI)
!! @li      2012-01-23 (Y.Ohno)     [mod] Self unpacking (MPI)
!! @li      2012-03-12 (H.Yashiro)  [mod] REAL4(MPI)
!! @li      2012-03-12 (Y.Ohno)     [mod] RDMA communication
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2012-03-27 (H.Yashiro)  [mod] Area/volume weighted total value report
!!
!<
#include "inc_openmp.h"
module mod_comm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use dc_types, only: &
     SP, &
     DP
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COMM_setup
  public :: COMM_vars
  public :: COMM_vars8
  public :: COMM_wait
#ifdef _USE_RDMA
  public :: COMM_set_rdma_variable
  public :: COMM_rdma_vars
  public :: COMM_rdma_vars8
#endif
  public :: COMM_stats
  public :: COMM_total
  public :: COMM_horizontal_mean

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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public, save :: COMM_total_doreport = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, save :: COMM_vsize_max       = 1024       !< # limit of communication variables at once
  logical,  private, save :: COMM_total_globalsum = .false.    !< calculate total with global communication?

  logical,  private, save :: COMM_IsAllPeriodic                !< periodic boundary condition?
  integer,  private, save :: COMM_datatype                     !< datatype of variable
  integer,  private, save :: COMM_world                        !< communication world ID

  integer,  private, save :: COMM_size3D_WE                    !< 3D data size (W/E    HALO, 4/8-direction comm.)
  integer,  private, save :: COMM_size3D_NS4                   !< 3D data size (N/S    HALO,   4-direction comm.)
  integer,  private, save :: COMM_size3D_NS8                   !< 3D data size (N/S    HALO,   8-direction comm.)
  integer,  private, save :: COMM_size3D_4C                    !< 3D data size (corner HALO,   8-direction comm.)

  integer,  private, save :: COMM_size2D_NS4                   !< 2D data size (W/E    HALO, 4/8-direction comm.)
  integer,  private, save :: COMM_size2D_NS8                   !< 2D data size (N/S    HALO,   4-direction comm.)
  integer,  private, save :: COMM_size2D_WE                    !< 2D data size (N/S    HALO,   8-direction comm.)
  integer,  private, save :: COMM_size2D_4C                    !< 2D data size (corner HALO,   8-direction comm.)

  real(RP), private, allocatable, save :: recvpack_W2P   (:,:) !< packing packet (receive, from W)
  real(RP), private, allocatable, save :: recvpack_E2P   (:,:) !< packing packet (receive, from E)
  real(RP), private, allocatable, save :: sendpack_P2W   (:,:) !< packing packet (send,    to W  )
  real(RP), private, allocatable, save :: sendpack_P2E   (:,:) !< packing packet (send,    to E  )

  integer,  private, allocatable, save :: req_cnt (:)          !< request ID of each MPI send/recv
  integer,  private, allocatable, save :: req_list(:,:)        !< request ID set of each variables

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine COMM_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop, &
       PRC_NEXT,    &
       PRC_W,       &
       PRC_N,       &
       PRC_E,       &
       PRC_S
    implicit none

    NAMELIST / PARAM_COMM / &
       COMM_vsize_max,      &
       COMM_total_doreport, &
       COMM_total_globalsum

    integer :: nreq_NS, nreq_WE, nreq_4C, nreq_MAX

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[COMM]/Categ[COMMON]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_COMM. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_COMM)

    ! only for register
    call TIME_rapstart('COMM vars MPI')
    call TIME_rapend  ('COMM vars MPI')
    call TIME_rapstart('COMM wait MPI')
    call TIME_rapend  ('COMM wait MPI')
    call TIME_rapstart('COMM Bcast MPI')
    call TIME_rapend  ('COMM Bcast MPI')
    call TIME_rapstart('COMM Allreduce MPI')
    call TIME_rapend  ('COMM Allreduce MPI')

    nreq_NS  = 2 * JHALO !--- send x JHALO, recv x JHALO
    nreq_WE  = 2         !--- send x 1    , recv x 1
    nreq_4C  = 2 * JHALO !--- send x JHALO, recv x JHALO
    nreq_MAX = 2 * nreq_NS + 2 * nreq_WE + 4 * nreq_4C

    COMM_size3D_NS4 = IA   * KA * JHALO
    COMM_size3D_NS8 = IMAX * KA
    COMM_size3D_WE  = JMAX * KA * IHALO
    COMM_size3D_4C  =        KA * IHALO

    COMM_size2D_NS4 = IA   * JHALO
    COMM_size2D_NS8 = IMAX
    COMM_size2D_WE  = JMAX * IHALO
    COMM_size2D_4C  =        IHALO

    allocate( recvpack_W2P   (COMM_size3D_WE,COMM_vsize_max) )
    allocate( recvpack_E2P   (COMM_size3D_WE,COMM_vsize_max) )
    allocate( sendpack_P2W   (COMM_size3D_WE,COMM_vsize_max) )
    allocate( sendpack_P2E   (COMM_size3D_WE,COMM_vsize_max) )

    allocate( req_cnt (COMM_vsize_max) )
    allocate( req_list(nreq_MAX,COMM_vsize_max) )
    req_cnt (:)   = 0
    req_list(:,:) = 0

    if (      PRC_NEXT(PRC_N) == MPI_PROC_NULL &
         .OR. PRC_NEXT(PRC_S) == MPI_PROC_NULL &
         .OR. PRC_NEXT(PRC_W) == MPI_PROC_NULL &
         .OR. PRC_NEXT(PRC_E) == MPI_PROC_NULL ) then
       COMM_IsAllPeriodic = .false.
    else
       COMM_IsAllPeriodic = .true.
    endif

    if ( RP == kind(0.D0) ) then
       COMM_datatype = MPI_DOUBLE_PRECISION
    elseif( RP == kind(0.0) ) then
       COMM_datatype = MPI_REAL
    else
       write(*,*) 'xxx precision is not supportd'
       call PRC_MPIstop
    endif

    COMM_world = MPI_COMM_WORLD

#ifdef _USE_RDMA
    call rdma_setup(COMM_vsize_max,  &
                    IA,              &
                    JA,              &
                    KA,              &
                    IHALO,           &
                    JHALO,           &
                    IS,              &
                    IE,              &
                    JS,              &
                    JE,              &
                    PRC_NEXT(PRC_W), &
                    PRC_NEXT(PRC_N), &
                    PRC_NEXT(PRC_E), &
                    PRC_NEXT(PRC_S)  )
#endif

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_3D(var, vid)
    use mod_process, only: &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    implicit none

    real(RP), intent(inout) :: var(:,:,:) !< 3D variable to communication
    integer,  intent(in)    :: vid        !< request ID

    integer :: ireq, tag

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    tag  = vid * 100
    ireq = 1

    call TIME_rapstart('COMM vars MPI')

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 4-Direction HALO communicate
       ! From S
       call MPI_IRECV( var(:,:,JS-JHALO:JS-1), COMM_size3D_NS4, COMM_datatype,      &
                       PRC_next(PRC_S), tag+1, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From N
       call MPI_IRECV( var(:,:,JE+1:JE+JHALO), COMM_size3D_NS4, COMM_datatype,      &
                       PRC_next(PRC_N), tag+2, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From E
       call MPI_IRECV( recvpack_E2P(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                       PRC_next(PRC_E), tag+3, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From W
       call MPI_IRECV( recvpack_W2P(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                       PRC_next(PRC_W), tag+4, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1

       ! packing packet to W HALO
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IS+IHALO-1
       do k = 1, KA
          n = (j-JS) * KA * IHALO &
            + (i-IS) * KA         &
            + k
          sendpack_P2W(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo
       ! packing packet to E HALO
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE-IHALO+1, IE
       do k = 1, KA
          n = (j-JS)         * KA * IHALO &
            + (i-IE+IHALO-1) * KA         &
            + k
          sendpack_P2E(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo

       !--- To 4-Direction HALO communicate
       ! To W HALO
       call MPI_ISEND( sendpack_P2W(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                       PRC_next(PRC_W), tag+3, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To E HALO
       call MPI_ISEND( sendpack_P2E(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                       PRC_next(PRC_E), tag+4, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To N HALO
       call MPI_ISEND( var(:,:,JE-JHALO+1:JE), COMM_size3D_NS4, COMM_datatype,      &
                       PRC_next(PRC_N), tag+1, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To S HALO
       call MPI_ISEND( var(:,:,JS:JS+JHALO-1), COMM_size3D_NS4, COMM_datatype,      &
                       PRC_next(PRC_S), tag+2, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1

    else ! non-periodic condition

       !--- From 4-Direction HALO communicate
       ! From S
       if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
          call MPI_IRECV( var(:,:,JS-JHALO:JS-1), COMM_size3D_NS4, COMM_datatype,      &
                          PRC_next(PRC_S), tag+1, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From N
       if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
          call MPI_IRECV( var(:,:,JE+1:JE+JHALO), COMM_size3D_NS4, COMM_datatype,      &
                          PRC_next(PRC_N), tag+2, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From E
       if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
          call MPI_IRECV( recvpack_E2P(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                          PRC_next(PRC_E), tag+3, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
          call MPI_IRECV( recvpack_W2P(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                          PRC_next(PRC_W), tag+4, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

       ! packing packet to W HALO
       if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
          !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IS+IHALO-1
          do k = 1, KA
             n = (j-JS) * KA * IHALO &
               + (i-IS) * KA         &
               + k
             sendpack_P2W(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif
       ! packing packet to E HALO
       if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
          !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IE-IHALO+1, IE
          do k = 1, KA
             n = (j-JS)         * KA * IHALO &
               + (i-IE+IHALO-1) * KA         &
               + k
             sendpack_P2E(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif

       !--- To 4-Direction HALO communicate
       ! To W HALO
       if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
          call MPI_ISEND( sendpack_P2W(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                          PRC_next(PRC_W), tag+3, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
          call MPI_ISEND( sendpack_P2E(:,vid),    COMM_size3D_WE,  COMM_datatype,      &
                          PRC_next(PRC_E), tag+4, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To N HALO
       if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
          call MPI_ISEND( var(:,:,JE-JHALO+1:JE), COMM_size3D_NS4, COMM_datatype,      &
                          PRC_next(PRC_N), tag+1, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To S HALO
       if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
          call MPI_ISEND( var(:,:,JS:JS+JHALO-1), COMM_size3D_NS4, COMM_datatype,      &
                          PRC_next(PRC_S), tag+2, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

    endif

    req_cnt(vid) = ireq - 1

    call TIME_rapend  ('COMM vars MPI')

    return
  end subroutine COMM_vars_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_3D(var, vid)
    use mod_process, only: &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S,    &
       PRC_NW,   &
       PRC_NE,   &
       PRC_SW,   &
       PRC_SE
    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer,  intent(in)    :: vid

    integer :: ireq, tag, tagc

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    tag  = vid * 100
    ireq = 1

    call TIME_rapstart('COMM vars MPI')

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- From 8-Direction HALO communicate
       ! From SE
       tagc = 0
       do j = JS-JHALO, JS-1
          call MPI_IRECV( var(1,IE+1,j),     COMM_size3D_4C, COMM_datatype,                &
                          PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From SW
       tagc = 10
       do j = JS-JHALO, JS-1
          call MPI_IRECV( var(1,IS-IHALO,j), COMM_size3D_4C, COMM_datatype,                &
                          PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From NE
       tagc = 20
       do j = JE+1, JE+JHALO
          call MPI_IRECV( var(1,IE+1,j),     COMM_size3D_4C, COMM_datatype,                &
                          PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From NW
       tagc = 30
       do j = JE+1, JE+JHALO
          call MPI_IRECV( var(1,IS-IHALO,j), COMM_size3D_4C, COMM_datatype,                &
                          PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From S
       tagc = 40
       do j = JS-JHALO, JS-1
          call MPI_IRECV( var(1,IS,j),     COMM_size3D_NS8, COMM_datatype,                &
                          PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From N
       tagc = 50
       do j = JE+1, JE+JHALO
          call MPI_IRECV( var(1,IS,j),     COMM_size3D_NS8, COMM_datatype,                &
                          PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! From E
       tagc = 60
       call MPI_IRECV( recvpack_E2P(:,vid), COMM_size3D_WE, COMM_datatype,             &
                       PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! From W
       tagc = 70
       call MPI_IRECV( recvpack_W2P(:,vid), COMM_size3D_WE, COMM_datatype,             &
                       PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1

       !--- packing packets to West
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IS+IHALO-1
       do k = 1, KA
          n = (j-JS) * KA * IHALO &
            + (i-IS) * KA         &
            + k
          sendpack_P2W(n,vid) = var(k,i,j)
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
          sendpack_P2E(n,vid) = var(k,i,j)
       enddo
       enddo
       enddo

       !--- To 8-Direction HALO communicate
       ! To W HALO
       tagc = 60
       call MPI_ISEND( sendpack_P2W(:,vid), COMM_size3D_WE, COMM_datatype,             &
                       PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To E HALO
       tagc = 70
       call MPI_ISEND( sendpack_P2E(:,vid), COMM_size3D_WE, COMM_datatype,             &
                       PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
       ireq = ireq + 1
       ! To N HALO
       tagc = 40
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j), COMM_size3D_NS8, COMM_datatype,                    &
                          PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To S HALO
       tagc = 50
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j), COMM_size3D_NS8, COMM_datatype,                    &
                          PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To NW HALO
       tagc = 0
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IS,j),         COMM_size3D_4C, COMM_datatype,              &
                          PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To NE HALO
       tagc = 10
       do j = JE-JHALO+1, JE
          call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size3D_4C, COMM_datatype,              &
                          PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To SW HALO
       tagc = 20
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IS,j),         COMM_size3D_4C, COMM_datatype,              &
                          PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo
       ! To SE HALO
       tagc = 30
       do j = JS, JS+JHALO-1
          call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size3D_4C, COMM_datatype,              &
                          PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
          tagc = tagc + 1
       enddo

    else ! non-periodic condition

       !--- From 8-Direction HALO communicate
       ! From SE
       if ( PRC_next(PRC_SE) /= MPI_PROC_NULL ) then
          tagc = 0
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IE+1,j),     COMM_size3D_4C, COMM_datatype,                &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From SW
       if ( PRC_next(PRC_SW) /= MPI_PROC_NULL ) then
          tagc = 10
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size3D_4C, COMM_datatype,                &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From NE
       if ( PRC_next(PRC_NE) /= MPI_PROC_NULL ) then
          tagc = 20
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IE+1,j),     COMM_size3D_4C, COMM_datatype,                &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From NW
       if ( PRC_next(PRC_NW) /= MPI_PROC_NULL ) then
          tagc = 30
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IS-IHALO,j), COMM_size3D_4C, COMM_datatype,                &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From S
       if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
          tagc = 40
          do j = JS-JHALO, JS-1
             call MPI_IRECV( var(1,IS,j),     COMM_size3D_NS8, COMM_datatype,                &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From N
       if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
          tagc = 50
          do j = JE+1, JE+JHALO
             call MPI_IRECV( var(1,IS,j),     COMM_size3D_NS8, COMM_datatype,                &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! From E
       if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
          tagc = 60
          call MPI_IRECV( recvpack_E2P(:,vid), COMM_size3D_WE, COMM_datatype,             &
                          PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! From W
       if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
          tagc = 70
          call MPI_IRECV( recvpack_W2P(:,vid), COMM_size3D_WE, COMM_datatype,             &
                          PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif

       if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
          !--- packing packets to West
          !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IS, IS+IHALO-1
          do k = 1, KA
             n = (j-JS) * KA * IHALO &
               + (i-IS) * KA         &
               + k
             sendpack_P2W(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif
       if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
          !--- packing packets to East
          !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
          do j = JS, JE
          do i = IE-IHALO+1, IE
          do k = 1, KA
             n = (j-JS)         * KA * IHALO &
               + (i-IE+IHALO-1) * KA         &
               + k
             sendpack_P2E(n,vid) = var(k,i,j)
          enddo
          enddo
          enddo
       endif

       !--- To 8-Direction HALO communicate
       ! To W HALO
       if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
          tagc = 60
          call MPI_ISEND( sendpack_P2W(:,vid), COMM_size3D_WE, COMM_datatype,             &
                          PRC_next(PRC_W), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To E HALO
       if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
          tagc = 70
          call MPI_ISEND( sendpack_P2E(:,vid), COMM_size3D_WE, COMM_datatype,             &
                          PRC_next(PRC_E), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
          ireq = ireq + 1
       endif
       ! To N HALO
       if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
          tagc = 40
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j), COMM_size3D_NS8, COMM_datatype,                    &
                             PRC_next(PRC_N), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To S HALO
       if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
          tagc = 50
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j), COMM_size3D_NS8, COMM_datatype,                    &
                             PRC_next(PRC_S), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To NW HALO
       if ( PRC_next(PRC_NW) /= MPI_PROC_NULL ) then
          tagc = 0
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IS,j),         COMM_size3D_4C, COMM_datatype,              &
                             PRC_next(PRC_NW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To NE HALO
       if ( PRC_next(PRC_NE) /= MPI_PROC_NULL ) then
          tagc = 10
          do j = JE-JHALO+1, JE
             call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size3D_4C, COMM_datatype,              &
                             PRC_next(PRC_NE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To SW HALO
       if ( PRC_next(PRC_SW) /= MPI_PROC_NULL ) then
          tagc = 20
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IS,j),         COMM_size3D_4C, COMM_datatype,              &
                             PRC_next(PRC_SW), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif
       ! To SE HALO
       if ( PRC_next(PRC_SE) /= MPI_PROC_NULL ) then
          tagc = 30
          do j = JS, JS+JHALO-1
             call MPI_ISEND( var(1,IE-IHALO+1,j), COMM_size3D_4C, COMM_datatype,              &
                             PRC_next(PRC_SE), tag+tagc, COMM_world, req_list(ireq,vid), ierr )
             ireq = ireq + 1
             tagc = tagc + 1
          enddo
       endif

    endif

    req_cnt(vid) = ireq - 1

    call TIME_rapend  ('COMM vars MPI')

    return
  end subroutine COMM_vars8_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_3D(var, vid)
    use mod_process, only: &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S

    implicit none

    real(RP), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM wait MPI')

    !--- wait packets
    call MPI_WAITALL( req_cnt (vid),                &
                      req_list(1:req_cnt(vid),vid), &
                      MPI_STATUSES_IGNORE,          &
                      ierr                          )

    if ( COMM_IsAllPeriodic ) then ! periodic condition

       !--- unpacking packets from East
       !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IE+1, IE+IHALO
       do k = 1, KA
          n = (j-JS)   * KA * IHALO &
            + (i-IE-1) * KA         &
            + k
          var(k,i,j) = recvpack_E2P(n,vid)
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
          var(k,i,j) = recvpack_W2P(n,vid)
       enddo
       enddo
       enddo

    else ! non-periodic condition

       if ( PRC_next(PRC_N) == MPI_PROC_NULL ) then
          !--- copy inner data to HALO(North)
          do j = JE+1, JE+JHALO
          do i = IS, IE
             var(:,i,j) = var(:,i,JE)
          enddo
          enddo
       endif

       if ( PRC_next(PRC_S) == MPI_PROC_NULL ) then
          !--- copy inner data to HALO(South)
          do j = JS-JHALO, JS-1
          do i = IS, IE
             var(:,i,j) = var(:,i,JS)
          enddo
          enddo
       endif

        if ( PRC_next(PRC_E) == MPI_PROC_NULL ) then
           !--- copy inner data to HALO(East)
           do j = JS, JE
           do i = IE+1, IE+IHALO
              var(:,i,j) = var(:,IE,j)
           enddo
           enddo
        else
           !--- unpacking packets from East
           !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
           do j = JS, JE
           do i = IE+1, IE+IHALO
           do k = 1, KA
              n = (j-JS)   * KA * IHALO &
                + (i-IE-1) * KA         &
                + k
              var(k,i,j) = recvpack_E2P(n,vid)
           enddo
           enddo
           enddo
        endif

        if ( PRC_next(PRC_W) == MPI_PROC_NULL ) then
           !--- copy inner data to HALO(West)
           do j = JS, JE
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,IS,j)
           enddo
           enddo
        else
           !--- unpacking packets from West
           !$omp parallel do private(i,j,k,n) OMP_SCHEDULE_ collapse(2)
           do j = JS, JE
           do i = IS-IHALO, IS-1
           do k = 1, KA
              n = (j-JS)       * KA * IHALO &
                + (i-IS+IHALO) * KA         &
                + k
              var(k,i,j) = recvpack_W2P(n,vid)
           enddo
           enddo
           enddo
        endif

        !--- copy inner data to HALO(NorthWest)
        if (       PRC_next(PRC_N) == MPI_PROC_NULL &
             .AND. PRC_next(PRC_W) == MPI_PROC_NULL ) then

           do j = JE+1, JE+JHALO
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,IS,JE)
           enddo
           enddo

        elseif( PRC_next(PRC_N) == MPI_PROC_NULL ) then

           do j = JE+1, JE+JHALO
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,i,JE)
           enddo
           enddo

        elseif( PRC_next(PRC_W) == MPI_PROC_NULL ) then

           do j = JE+1, JE+JHALO
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,IS,j)
           enddo
           enddo

        endif

        !--- copy inner data to HALO(SouthWest)
        if (       PRC_next(PRC_S) == MPI_PROC_NULL &
             .AND. PRC_next(PRC_W) == MPI_PROC_NULL ) then

           do j = JS-IHALO, JS-1
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,IS,JS)
           enddo
           enddo

        elseif( PRC_next(PRC_S) == MPI_PROC_NULL ) then

           do j = JS-IHALO, JS-1
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,i,JS)
           enddo
           enddo

        elseif( PRC_next(PRC_W) == MPI_PROC_NULL ) then

           do j = JS-IHALO, JS-1
           do i = IS-IHALO, IS-1
              var(:,i,j) = var(:,IS,j)
           enddo
           enddo

        endif

        !--- copy inner data to HALO(NorthEast)
        if (       PRC_next(PRC_N) == MPI_PROC_NULL &
             .AND. PRC_next(PRC_E) == MPI_PROC_NULL ) then

            do j = JE+1, JE+JHALO
            do i = IE+1, IE+IHALO
               var(:,i,j) = var(:,IE,JE)
            enddo
            enddo

        elseif( PRC_next(PRC_N) == MPI_PROC_NULL ) then

            do j = JE+1, JE+JHALO
            do i = IE+1, IE+IHALO
               var(:,i,j) = var(:,i,JE)
            enddo
            enddo

        elseif( PRC_next(PRC_E) == MPI_PROC_NULL ) then

            do j = JE+1, JE+JHALO
            do i = IE+1, IE+IHALO
               var(:,i,j) = var(:,IE,j)
            enddo
            enddo

        endif

        !--- copy inner data to HALO(SouthEast)
        if (       PRC_next(PRC_S) == MPI_PROC_NULL &
             .AND. PRC_next(PRC_E) == MPI_PROC_NULL ) then

            do j = JS-IHALO, JS-1
            do i = IE+1, IE+IHALO
               var(:,i,j) = var(:,IE,JS)
            enddo
            enddo

        elseif( PRC_next(PRC_S) == MPI_PROC_NULL ) then

            do j = JS-IHALO, JS-1
            do i = IE+1, IE+IHALO
               var(:,i,j) = var(:,i,JS)
            enddo
            enddo

        elseif( PRC_next(PRC_E) == MPI_PROC_NULL ) then

            do j = JS-IHALO, JS-1
            do i = IE+1, IE+IHALO
               var(:,i,j) = var(:,IE,j)
            enddo
            enddo

        endif

    endif

    call TIME_rapend  ('COMM wait MPI')

    return
  end subroutine COMM_wait_3D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_2D(var, vid)
    use mod_process, only: &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)    :: vid

    integer :: ireq, tag
    integer :: ierr
    integer :: i, j, n
    !---------------------------------------------------------------------------

    tag = vid * 100
    ireq = 1

    call TIME_rapstart('COMM vars MPI')

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

        !--- To 4-Direction HALO communicate
        !--- packing packets to West
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IS, IS+IHALO-1
            n =  (j-JS) * IHALO &
               + (i-IS) + 1
            sendpack_P2W(n,vid) = var(i,j)
        enddo
        enddo

        !--- packing packets to East
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IE-IHALO+1, IE
            n =  (j-JS)         * IHALO &
               + (i-IE+IHALO-1) + 1
            sendpack_P2E(n,vid) = var(i,j)
        enddo
        enddo

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
        if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
            call MPI_IRECV( var(:,JS-JHALO:JS-1), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_S), tag+1, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

        ! From N
        if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
            call MPI_IRECV( var(:,JE+1:JE+JHALO), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_N), tag+2, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

        ! From E
        if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_E), tag+3, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! From W
        if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_W), tag+4, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        !--- To 4-Direction HALO communicate
        !--- packing packets to West
        if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            do j = JS, JE
            do i = IS, IS+IHALO-1
                n =  (j-JS) * IHALO &
                   + (i-IS) + 1
                sendpack_P2W(n,vid) = var(i,j)
            enddo
            enddo
        endif

        !--- packing packets to East
        if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            do j = JS, JE
            do i = IE-IHALO+1, IE
                n =  (j-JS)         * IHALO &
                   + (i-IE+IHALO-1) + 1
                sendpack_P2E(n,vid) = var(i,j)
            enddo
            enddo
         endif

        ! To W HALO communicate
        if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_W), tag+3, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! To E HALO communicate
        if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,       &
                            COMM_datatype, PRC_next(PRC_E), tag+4, &
                            COMM_world, req_list(ireq,vid), ierr )
            ireq = ireq + 1
        endif

        ! To N HALO communicate
        if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
            call MPI_ISEND( var(:,JE-JHALO+1:JE), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_N), tag+1, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

        ! To S HALO communicate
        if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
            call MPI_ISEND( var(:,JS:JS+JHALO-1), COMM_size2D_NS4,        &
                            COMM_datatype, PRC_next(PRC_S), tag+2, &
                            COMM_world, req_list(ireq,vid), ierr    )
            ireq = ireq + 1
        endif

    endif

    req_cnt(vid) = ireq - 1

    call TIME_rapend  ('COMM vars MPI')

    return
  end subroutine COMM_vars_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_2D(var, vid)
    use mod_process, only: &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S,    &
       PRC_NW,   &
       PRC_NE,   &
       PRC_SW,   &
       PRC_SE
    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)    :: vid

    integer :: ireq, tag, tagc

    integer :: ierr
    integer :: i, j, n
    !---------------------------------------------------------------------------

    tag   = vid * 100
    ireq = 1

    call TIME_rapstart('COMM vars MPI')

    if ( COMM_IsAllPeriodic ) then
    !--- periodic condition
        !--- From 8-Direction HALO communicate
        ! From SE
        tagc = 0
        do j = JS-JHALO, JS-1
            call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                            COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From SW
        tagc = 10
        do j = JS-JHALO, JS-1
            call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                            COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From NE
        tagc = 20
        do j = JE+1, JE+JHALO
            call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                            COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From NW
        tagc = 30
        do j = JE+1, JE+JHALO
            call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                            COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From S
        tagc = 40
        do j = JS-JHALO, JS-1
            call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
             ireq = ireq + 1
             tagc = tagc + 1
        enddo
        ! From N
        tagc = 50
        do j = JE+1, JE+JHALO
            call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
        ! From E
        call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,           &
                        COMM_datatype, PRC_next(PRC_E), tag+60, &
                        MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
        ireq = ireq + 1
        ! From W
        call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,           &
                        COMM_datatype, PRC_next(PRC_W), tag+70, &
                        MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
        ireq = ireq + 1

        !--- To 8-Direction HALO communicate
        !--- packing packets to West
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IS, IS+IHALO-1
            n =  (j-JS) * IHALO &
               + (i-IS) + 1
            sendpack_P2W(n,vid) = var(i,j)
        enddo
        enddo

        !--- packing packets to East
        !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
        do j = JS, JE
        do i = IE-IHALO+1, IE
            n =  (j-JS)         * IHALO &
               + (i-IE+IHALO-1) + 1
            sendpack_P2E(n,vid) = var(i,j)
        enddo
        enddo

        ! To W HALO communicate
        call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,            &
                        COMM_datatype, PRC_next(PRC_W), tag+60, &
                        MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
        ireq = ireq + 1

        ! To E HALO communicate
        call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,           &
                        COMM_datatype, PRC_next(PRC_E), tag+70, &
                        MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
        ireq = ireq + 1

        ! To N HALO communicate
        tagc = 40
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To S HALO communicate
        tagc = 50
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                            COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To NW HALO communicate
        tagc = 0
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                            COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To NE HALO communicate
        tagc = 10
        do j = JE-JHALO+1, JE
            call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                            COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To SW HALO communicate
        tagc = 20
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                            COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo

        ! To SE HALO communicate
        tagc = 30
        do j = JS, JS+JHALO-1
            call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                            COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
            ireq = ireq + 1
            tagc = tagc + 1
        enddo
    else
    !--- non-periodic condition
        !--- From 8-Direction HALO communicate
        ! From SE
        if ( PRC_next(PRC_SE) /= MPI_PROC_NULL ) then
            tagc = 0
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From SW
        if ( PRC_next(PRC_SW) /= MPI_PROC_NULL ) then
            tagc = 10
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From NE
        if ( PRC_next(PRC_NE) /= MPI_PROC_NULL ) then
            tagc = 20
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IE+1,j), COMM_size2D_4C,                      &
                                COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From NW
        if ( PRC_next(PRC_NW) /= MPI_PROC_NULL ) then
            tagc = 30
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IS-IHALO,j), COMM_size2D_4C,                  &
                                COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From S
        if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
            tagc = 40
            do j = JS-JHALO, JS-1
                call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
                 ireq = ireq + 1
                 tagc = tagc + 1
            enddo
        endif

        ! From N
        if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
            tagc = 50
            do j = JE+1, JE+JHALO
                call MPI_IRECV( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! From E
        if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            call MPI_IRECV( recvpack_E2P(:,vid), COMM_size2D_WE,             &
                            COMM_datatype, PRC_next(PRC_E), tag+60, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        ! From W
        if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            call MPI_IRECV( recvpack_W2P(:,vid), COMM_size2D_WE,           &
                            COMM_datatype, PRC_next(PRC_W), tag+70, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        !--- To 8-Direction HALO communicate
        !--- packing packets to West
        if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
            do j = JS, JE
            do i = IS, IS+IHALO-1
                n =  (j-JS) * IHALO &
                   + (i-IS) + 1
                sendpack_P2W(n,vid) = var(i,j)
            enddo
            enddo
        endif

        !--- packing packets to East
        if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
            do j = JS, JE
            do i = IE-IHALO+1, IE
                n =  (j-JS)         * IHALO &
                   + (i-IE+IHALO-1) + 1
                sendpack_P2E(n,vid) = var(i,j)
            enddo
            enddo
        endif

        ! To W HALO communicate
        if ( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            call MPI_ISEND( sendpack_P2W(:,vid), COMM_size2D_WE,           &
                            COMM_datatype, PRC_next(PRC_W), tag+60, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        ! To E HALO communicate
        if ( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            call MPI_ISEND( sendpack_P2E(:,vid), COMM_size2D_WE,           &
                            COMM_datatype, PRC_next(PRC_E), tag+70, &
                            MPI_COMM_WORLD, req_list(ireq,vid), ierr     )
            ireq = ireq + 1
        endif

        ! To N HALO communicate
        if ( PRC_next(PRC_N) /= MPI_PROC_NULL ) then
            tagc = 40
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_N), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To S HALO communicate
        if ( PRC_next(PRC_S) /= MPI_PROC_NULL ) then
            tagc = 50
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IS,j), COMM_size2D_NS8,                      &
                                COMM_datatype, PRC_next(PRC_S), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr       )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To NW HALO communicate
        if ( PRC_next(PRC_NW) /= MPI_PROC_NULL ) then
            tagc = 0
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_NW), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To NE HALO communicate
        if ( PRC_next(PRC_NE) /= MPI_PROC_NULL ) then
            tagc = 10
            do j = JE-JHALO+1, JE
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_NE), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To SW HALO communicate
        if ( PRC_next(PRC_SW) /= MPI_PROC_NULL ) then
            tagc = 20
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IS,j), COMM_size2D_4C,                        &
                                COMM_datatype, PRC_next(PRC_SW), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

        ! To SE HALO communicate
        if ( PRC_next(PRC_SE) /= MPI_PROC_NULL ) then
            tagc = 30
            do j = JS, JS+JHALO-1
                call MPI_ISEND( var(IE-IHALO+1,j), COMM_size2D_4C,                &
                                COMM_datatype, PRC_next(PRC_SE), tag+tagc, &
                                MPI_COMM_WORLD, req_list(ireq,vid), ierr        )
                ireq = ireq + 1
                tagc = tagc + 1
            enddo
        endif

    endif

    req_cnt(vid) = ireq - 1

    call TIME_rapend  ('COMM vars MPI')

    return
  end subroutine COMM_vars8_2D

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_2D(var, vid)
    use mod_process, only: &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S

    implicit none

    real(RP), intent(inout) :: var(:,:)
    integer, intent(in)    :: vid

    integer :: ierr
    integer :: i, j, n
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM wait MPI')

    !--- wait packets
    call MPI_WAITALL(req_cnt(vid), req_list(1:req_cnt(vid),vid), MPI_STATUSES_IGNORE, ierr)

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

        !--- copy inner data to HALO(North)
        if( PRC_next(PRC_N) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IS, IE
               var(i,j) = var(i,JE)
            enddo
            enddo
        endif

        !--- copy inner data to HALO(South)
        if( PRC_next(PRC_S) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-JHALO, JS-1
            do i = IS, IE
               var(i,j) = var(i,JS)
            enddo
            enddo
        endif

        !--- unpacking packets from East / copy inner data to HALO(East)
        if( PRC_next(PRC_E) /= MPI_PROC_NULL ) then
            !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
            do j = JS, JE
            do i = IE+1, IE+IHALO
               n = (j-JS)   * IHALO &
                 + (i-IE-1) + 1
               var(i,j) = recvpack_E2P(n,vid)
            enddo
            enddo
        else
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS, JE
            do i = IE+1, IE+IHALO
               var(i,j) = var(IE,j)
            enddo
            enddo
        endif

        !--- unpacking packets from West / copy inner data to HALO(West)
        if( PRC_next(PRC_W) /= MPI_PROC_NULL ) then
            !$omp parallel do private(i,j,n) OMP_SCHEDULE_ collapse(2)
            do j = JS, JE
            do i = IS-IHALO, IS-1
               n = (j-JS)       * IHALO &
                 + (i-IS+IHALO) + 1
               var(i,j) = recvpack_W2P(n,vid)
            enddo
            enddo
        else
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS, JE
            do i = IS-IHALO, IS-1
               var(i,j) = var(IS,j)
            enddo
            enddo
        endif

        !--- copy inner data to HALO(NorthWest)
        if( PRC_next(PRC_N) == MPI_PROC_NULL .AND. PRC_next(PRC_W) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IS-IHALO, IS-1
               var(i,j) = var(IS,JE)
            enddo
            enddo
        elseif( PRC_next(PRC_N) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IS-IHALO, IS-1
               var(i,j) = var(i,JE)
            enddo
            enddo
        elseif( PRC_next(PRC_W) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IS-IHALO, IS-1
               var(i,j) = var(IS,j)
            enddo
            enddo
        endif

        !--- copy inner data to HALO(SouthWest)
        if( PRC_next(PRC_S) == MPI_PROC_NULL .AND. PRC_next(PRC_W) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-IHALO, JS-1
            do i = IS-IHALO, IS-1
               var(i,j) = var(IS,JS)
            enddo
            enddo
        elseif( PRC_next(PRC_S) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-IHALO, JS-1
            do i = IS-IHALO, IS-1
               var(i,j) = var(i,JS)
            enddo
            enddo
        elseif( PRC_next(PRC_W) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-IHALO, JS-1
            do i = IS-IHALO, IS-1
               var(i,j) = var(IS,j)
            enddo
            enddo
        endif

        !--- copy inner data to HALO(NorthEast)
        if( PRC_next(PRC_N) == MPI_PROC_NULL .AND. PRC_next(PRC_E) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IE+1, IE+IHALO
               var(i,j) = var(IE,JE)
            enddo
            enddo
        elseif( PRC_next(PRC_N) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IE+1, IE+IHALO
               var(i,j) = var(i,JE)
            enddo
            enddo
        elseif( PRC_next(PRC_E) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JE+1, JE+JHALO
            do i = IE+1, IE+IHALO
               var(i,j) = var(IE,j)
            enddo
            enddo
        endif

        !--- copy inner data to HALO(SouthEast)
        if( PRC_next(PRC_S) == MPI_PROC_NULL .AND. PRC_next(PRC_E) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-IHALO, JS-1
            do i = IE+1, IE+IHALO
               var(i,j) = var(IE,JS)
            enddo
            enddo
        elseif( PRC_next(PRC_S) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-IHALO, JS-1
            do i = IE+1, IE+IHALO
               var(i,j) = var(i,JS)
            enddo
            enddo
        elseif( PRC_next(PRC_E) == MPI_PROC_NULL ) then
            !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
            do j = JS-IHALO, JS-1
            do i = IE+1, IE+IHALO
               var(i,j) = var(IE,j)
            enddo
            enddo
        endif

    endif


    call TIME_rapend  ('COMM wait MPI')

    return
  end subroutine COMM_wait_2D

#ifdef _USE_RDMA
  !-----------------------------------------------------------------------------
  !> Register variables on the memory for direct memory access (RDMA)
  subroutine COMM_set_rdma_variable(var, vid)
    implicit none

    real(RP), intent(in) :: var(:,:,:) !< variable for register
    integer,  intent(in) :: vid        !< variable ID
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** set RDMA ID:', vid-1

    call set_rdma_variable(var, vid-1);

    return
  end subroutine

  !-----------------------------------------------------------------------------
  !> Communicate halo value with 4 directions (RDMA)
  subroutine COMM_rdma_vars(vid, num)
    implicit none

    integer, intent(in) :: vid !< variable ID
    integer, intent(in) :: num !< number of variables from vid
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM RDMA')

    !--- put data
    call rdma_put(vid-1, num)

    call TIME_rapend  ('COMM RDMA')

    return
  end subroutine COMM_rdma_vars

  !-----------------------------------------------------------------------------
  !> Communicate halo value with 8 directions (RDMA)
  subroutine COMM_rdma_vars8(vid, num)
    implicit none

    integer, intent(in) :: vid !< variable ID
    integer, intent(in) :: num !< number of variables from vid
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM RDMA')

    !--- put data
    call rdma_put8(vid-1, num)

    call TIME_rapend  ('COMM RDMA')

    return
  end subroutine COMM_rdma_vars8
#endif

  !-----------------------------------------------------------------------------
  !> Search global maximum & minimum value
  subroutine COMM_stats(var, varname)
    use mod_process, only: &
       PRC_nmax,   &
       PRC_myrank
    use mod_const, only: &
       CONST_UNDEF8, &
       CONST_UNDEF2
    implicit none

    real(RP),         intent(inout) :: var(:,:,:,:) !< values
    character(len=*), intent(in)    :: varname(:)   !< name of item

    logical :: halomask(KA,IA,JA)

    real(RP), allocatable :: statval   (:,:,:)
    integer,  allocatable :: statidx   (:,:,:,:)
    real(RP), allocatable :: allstatval(:,:)
    integer,  allocatable :: allstatidx(:,:,:)
    integer               :: vsize

    integer :: ierr
    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:,:,:,:),4)

    halomask(:,:,:) = .false.
    halomask(KS:KE,IS:IE,JS:JE) = .true.

    allocate( statval(  vsize,2,0:PRC_nmax-1) ); statval(:,:,:)   = CONST_UNDEF8
    allocate( statidx(3,vsize,2,0:PRC_nmax-1) ); statidx(:,:,:,:) = CONST_UNDEF2

    allocate( allstatval(  vsize,2) ); allstatval(:,:)   = CONST_UNDEF8
    allocate( allstatidx(1,vsize,2) ); allstatidx(:,:,:) = CONST_UNDEF2

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Variable Statistics ***'
    do v = 1, vsize
       statval(  v,1,PRC_myrank) = maxval(var(:,:,:,v),mask=halomask)
       statval(  v,2,PRC_myrank) = minval(var(:,:,:,v),mask=halomask)
       statidx(:,v,1,PRC_myrank) = maxloc(var(:,:,:,v),mask=halomask)
       statidx(:,v,2,PRC_myrank) = minloc(var(:,:,:,v),mask=halomask)

! statistics on each node
!       if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,3(I5,A))') '*** MAX = ', &
!                                             statval(  v,1,PRC_myrank),'(', &
!                                             statidx(1,v,1,PRC_myrank),',', &
!                                             statidx(2,v,1,PRC_myrank),',', &
!                                             statidx(3,v,1,PRC_myrank),')'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,3(I5,A))') '*** MIN = ', &
!                                             statval(  v,2,PRC_myrank),'(', &
!                                             statidx(1,v,2,PRC_myrank),',', &
!                                             statidx(2,v,2,PRC_myrank),',', &
!                                             statidx(3,v,2,PRC_myrank),')'
    enddo

    ! MPI broadcast
    call TIME_rapstart('COMM Bcast MPI')
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,1,p),       &
                       vsize*2,              &
                       COMM_datatype,          &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
       call MPI_Bcast( statidx(1,1,1,p),     &
                       3*vsize*2,            &
                       MPI_INTEGER,          &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
    enddo
    call TIME_rapend  ('COMM Bcast MPI')

    do v = 1, vsize
       allstatval(v,1)   = maxval(statval(v,1,:))
       allstatval(v,2)   = minval(statval(v,2,:))
       allstatidx(:,v,1) = maxloc(statval(v,1,:))-1
       allstatidx(:,v,2) = minloc(statval(v,2,:))-1
       if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,4(I5,A))') '  MAX =', &
                                                    allstatval(  v,1), '(', &
                                                    allstatidx(1,v,1), ',', &
                                      statidx(1,v,1,allstatidx(1,v,1)),',', &
                                      statidx(2,v,1,allstatidx(1,v,1)),',', &
                                      statidx(3,v,1,allstatidx(1,v,1)),')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,4(I5,A))') '  MIN =', &
                                                    allstatval(  v,2), '(', &
                                                    allstatidx(1,v,2), ',', &
                                      statidx(1,v,2,allstatidx(1,v,2)),',', &
                                      statidx(2,v,2,allstatidx(1,v,2)),',', &
                                      statidx(3,v,2,allstatidx(1,v,2)),')'
    enddo

    deallocate( statval )
    deallocate( statidx )

    deallocate( allstatval )
    deallocate( allstatidx )

    return
  end subroutine COMM_stats

  !-----------------------------------------------------------------------------
  !> Calc volume/area-weighted global sum
  subroutine COMM_total( allstatval, var, varname )
    use mod_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use mod_geometrics, only: &
       area => GEOMETRICS_area, &
       vol  => GEOMETRICS_vol
    implicit none

    real(RP),         intent(out) :: allstatval !< volume/area-weighted total
    real(RP),         intent(in)  :: var(:,:,:) !< 3D value
    character(len=*), intent(in)  :: varname    !< name of item

    real(RP) :: statval
    integer  :: ksize

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    ksize = size(var(:,:,:),1)

    statval = 0.0_RP
    if ( ksize == KA ) then ! 3D
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          statval = statval + var(k,i,j) * vol(k,i,j)
       enddo
       enddo
       enddo
    elseif( ksize == 1 ) then ! 2D
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2) reduction(+:statval)
       do j = JS, JE
       do i = IS, IE
          statval = statval + var(1,i,j) * area(1,i,j)
       enddo
       enddo
    endif

    if ( .not. ( statval > -1.0_RP .or. statval < 1.0_RP ) ) then ! must be NaN
       write(*,*) 'xxx NaN is detected for ', trim(varname), ' in rank ', PRC_myrank
       call PRC_MPIstop
    endif

    if ( COMM_total_globalsum ) then
       call TIME_rapstart('COMM Allreduce MPI')
       ! All reduce
       call MPI_Allreduce( statval,              &
                           allstatval,           &
                           1,                    &
                           COMM_datatype,          &
                           MPI_SUM,              &
                           MPI_COMM_WORLD,       &
                           ierr                  )

       call TIME_rapend  ('COMM Allreduce MPI')

       ! statistics over the all node
       if ( varname /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,1PE24.17)') &
                     '[', varname, '] SUM(global) =', allstatval
       endif
    else
       allstatval = statval

       ! statistics on each node
       if ( varname /= "" ) then ! if varname is empty, suppress output
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,1PE24.17)') &
                     '[', varname, '] SUM(local)  =', statval
       endif
    endif

    return
  end subroutine COMM_total

  !-----------------------------------------------------------------------------
  !> calculate horizontal mean (global total with communication)
  subroutine COMM_horizontal_mean( varmean, var )
    use mod_process, only: &
       PRC_nmax
    implicit none

    real(RP), intent(out) :: varmean(KA)       !< horizontal mean
    real(RP), intent(in)  :: var    (KA,IA,JA) !< 3D value

    real(RP) :: statval   (KA)
    real(RP) :: allstatval(KA)

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    statval(:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       statval(k) = statval(k) + var(k,i,j)
    enddo
    enddo
    enddo

    do k = KS, KE
       statval(k) = statval(k) / real(IMAX*JMAX,kind=RP)
    enddo

    ! [NOTE] always communicate globally
    call TIME_rapstart('COMM Allreduce MPI')
    ! All reduce
    call MPI_Allreduce( statval(1),           &
                        allstatval(1),        &
                        KA,                   &
                        COMM_datatype,          &
                        MPI_SUM,              &
                        MPI_COMM_WORLD,       &
                        ierr                  )

    call TIME_rapend  ('COMM Allreduce MPI')

    do k = KS, KE
       varmean(k) = allstatval(k) / real(PRC_nmax,kind=RP) 
    enddo
    varmean(   1:KS-1) = 0.0_RP
    varmean(KE+1:KA  ) = 0.0_RP

    return
  end subroutine COMM_horizontal_mean

end module mod_comm
!-------------------------------------------------------------------------------
