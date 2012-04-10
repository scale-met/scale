!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI module for SCALE3 (Communication Core)
!!
!! @author H.Tomita and SCALE developpers
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
!-------------------------------------------------------------------------------
module mod_comm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
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
  !++ Public procedure
  !
  public :: COMM_setup
  public :: COMM_vars
  public :: COMM_vars8
  public :: COMM_wait
  public :: COMM_vars_r4
  public :: COMM_vars8_r4
  public :: COMM_wait_r4
#ifdef _USE_RDMA
  public :: COMM_set_rdma_variable
  public :: COMM_rdma_vars
  public :: COMM_rdma_vars8
#endif
  public :: COMM_stats
  public :: COMM_total

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public, save :: COMM_total_doreport  = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, save :: COMM_vsize_max  = 30
  logical, private, save :: COMM_total_globalsum = .false.

  integer, private, save :: datasize_NS4
  integer, private, save :: datasize_NS8
  integer, private, save :: datasize_WE
  integer, private, save :: datasize_4C

  integer, private, save :: IREQ_CNT_NS
  integer, private, save :: IREQ_CNT_WE
  integer, private, save :: IREQ_CNT_4C
  integer, private, save :: IREQ_CNT_MAX

  real(8), private, allocatable, save :: recvpack_W2P(:,:)
  real(8), private, allocatable, save :: recvpack_E2P(:,:)
  real(8), private, allocatable, save :: sendpack_P2W(:,:)
  real(8), private, allocatable, save :: sendpack_P2E(:,:)
  real(4), private, allocatable, save :: recvpack_W2P_r4(:,:)
  real(4), private, allocatable, save :: recvpack_E2P_r4(:,:)
  real(4), private, allocatable, save :: sendpack_P2W_r4(:,:)
  real(4), private, allocatable, save :: sendpack_P2E_r4(:,:)

  integer, private, allocatable, save :: ireq_cnt(:)
  integer, private, allocatable, save :: ireq_list(:,:)

  !-----------------------------------------------------------------------------
contains

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

    IREQ_CNT_NS  = 2 * JHALO !--- sendxJHALO recvxJHALO
    IREQ_CNT_WE  = 2         !--- sendx1 recvx1
    IREQ_CNT_4C  = 2 * JHALO !--- sendxJHALO recvxJHALO
    IREQ_CNT_MAX = 2 * IREQ_CNT_NS + 2 * IREQ_CNT_WE + 4 * IREQ_CNT_4C

    datasize_NS4 = IA   * KA * JHALO
    datasize_NS8 = IMAX * KA
    datasize_WE  = JMAX * KA * IHALO
    datasize_4C  =        KA * IHALO

    allocate( recvpack_W2P(datasize_WE,COMM_vsize_max) )
    allocate( recvpack_E2P(datasize_WE,COMM_vsize_max) )
    allocate( sendpack_P2W(datasize_WE,COMM_vsize_max) )
    allocate( sendpack_P2E(datasize_WE,COMM_vsize_max) )

    allocate( recvpack_W2P_r4(datasize_WE,COMM_vsize_max) )
    allocate( recvpack_E2P_r4(datasize_WE,COMM_vsize_max) )
    allocate( sendpack_P2W_r4(datasize_WE,COMM_vsize_max) )
    allocate( sendpack_P2E_r4(datasize_WE,COMM_vsize_max) )

    allocate( ireq_cnt(COMM_vsize_max) ) ;              ireq_cnt(:)   = 0
    allocate( ireq_list(IREQ_CNT_MAX,COMM_vsize_max) ); ireq_list(:,:) = 0

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
  subroutine COMM_vars(var, vid)
    use mod_process, only : &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    implicit none

    real(8), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ireqc, tag
    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    tag = vid * 100
    ireqc = 1

    call TIME_rapstart('COMM vars MPI')

    !-- From 4-Direction HALO communicate
    ! From S
    call MPI_IRECV( var(:,:,JS-JHALO:JS-1), datasize_NS4,         &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag+1, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From N
    call MPI_IRECV( var(:,:,JE+1:JE+JHALO), datasize_NS4,         &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag+2, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From E
    call MPI_IRECV( recvpack_E2P(:,vid), datasize_WE,             &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_E), tag+3, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! From W
    call MPI_IRECV( recvpack_W2P(:,vid), datasize_WE,             &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_W), tag+4, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr )
    ireqc = ireqc + 1

    !-- To 4-Direction HALO communicate
    !--- packing packets to West
    do j = JS, JE
    do i = IS, IS+IHALO-1
    do k = 1, KA
        n =  (j-JS) * KA * IHALO &
           + (i-IS) * KA         &
           + k
        sendpack_P2W(n,vid) = var(k,i,j)
    enddo
    enddo
    enddo

    !--- packing packets to East
    do j = JS, JE
    do i = IE-IHALO+1, IE
    do k = 1, KA
        n =  (j-JS)         * KA * IHALO &
           + (i-IE+IHALO-1) * KA         &
           + k
        sendpack_P2E(n,vid) = var(k,i,j)
        n = n + 1
    enddo
    enddo
    enddo

    ! To W HALO communicate
    call MPI_ISEND( sendpack_P2W(:,vid), datasize_WE,             &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_W), tag+3, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! To E HALO communicate
    call MPI_ISEND( sendpack_P2E(:,vid), datasize_WE,             &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_E), tag+4, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! To N HALO communicate
    call MPI_ISEND( var(:,:,JE-JHALO+1:JE), datasize_NS4,         &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag+1, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! To S HALO communicate
    call MPI_ISEND( var(:,:,JS:JS+JHALO-1), datasize_NS4,         &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag+2, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ireq_cnt(vid) = ireqc - 1

    call TIME_rapend  ('COMM vars MPI')

    return
  end subroutine COMM_vars

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8(var, vid)
    use mod_process, only : &
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

    real(8), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ireqc, tag, tagc

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    tag   = vid * 100
    ireqc = 1

    call TIME_rapstart('COMM vars MPI')

    !-- From 8-Direction HALO communicate
    ! From SE
    tagc = 0
    do j = JS-JHALO, JS-1
        call MPI_IRECV( var(1,IE+1,j), datasize_4C,                       &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_SE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From SW
    tagc = 10
    do j = JS-JHALO, JS-1
        call MPI_IRECV( var(1,IS-IHALO,j), datasize_4C,                   &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_SW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From NE
    tagc = 20
    do j = JE+1, JE+JHALO
        call MPI_IRECV( var(1,IE+1,j), datasize_4C,                       &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_NE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From NW
    tagc = 30
    do j = JE+1, JE+JHALO
        call MPI_IRECV( var(1,IS-IHALO,j), datasize_4C,                   &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_NW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From S
    tagc = 40
    do j = JS-JHALO, JS-1
        call MPI_IRECV( var(1,IS,j), datasize_NS8,                       &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
         ireqc = ireqc + 1
         tagc  = tagc  + 1
    enddo
    ! From N
    tagc = 50
    do j = JE+1, JE+JHALO
        call MPI_IRECV( var(1,IS,j), datasize_NS8,                       &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From E
    call MPI_IRECV( recvpack_E2P(:,vid), datasize_WE,              &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_E), tag+60, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1
    ! From W
    call MPI_IRECV( recvpack_W2P(:,vid), datasize_WE,              &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_W), tag+70, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1

    !-- To 8-Direction HALO communicate
    !--- packing packets to West
    do j = JS, JE
    do i = IS, IS+IHALO-1
    do k = 1, KA
        n =  (j-JS) * KA * IHALO &
           + (i-IS) * KA         &
           + k
        sendpack_P2W(n,vid) = var(k,i,j)
    enddo
    enddo
    enddo

    !--- packing packets to East
    do j = JS, JE
    do i = IE-IHALO+1, IE
    do k = 1, KA
        n =  (j-JS)         * KA * IHALO &
           + (i-IE+IHALO-1) * KA         &
           + k
        sendpack_P2E(n,vid) = var(k,i,j)
        n = n + 1
    enddo
    enddo
    enddo

    ! To W HALO communicate
    call MPI_ISEND( sendpack_P2W(:,vid), datasize_WE,              &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_W), tag+60, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1

    ! To E HALO communicate
    call MPI_ISEND( sendpack_P2E(:,vid), datasize_WE,              &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_E), tag+70, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1

    ! To N HALO communicate
    tagc = 40
    do j = JE-JHALO+1, JE
        call MPI_ISEND( var(1,IS,j), datasize_NS8,                       &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To S HALO communicate
    tagc = 50
    do j = JS, JS+JHALO-1
        call MPI_ISEND( var(1,IS,j), datasize_NS8,                       &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To NW HALO communicate
    tagc = 0
    do j = JE-JHALO+1, JE
        call MPI_ISEND( var(1,IS,j), datasize_4C,                         &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_NW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To NE HALO communicate
    tagc = 10
    do j = JE-JHALO+1, JE
        call MPI_ISEND( var(1,IE-IHALO+1,j), datasize_4C,                 &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_NE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To SW HALO communicate
    tagc = 20
    do j = JS, JS+JHALO-1
        call MPI_ISEND( var(1,IS,j), datasize_4C,                         &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_SW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To SE HALO communicate
    tagc = 30
    do j = JS, JS+JHALO-1
        call MPI_ISEND( var(1,IE-IHALO+1,j), datasize_4C,                 &
                        MPI_DOUBLE_PRECISION, PRC_next(PRC_SE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ireq_cnt(vid) = ireqc - 1

    call TIME_rapend  ('COMM vars MPI')

    return
  end subroutine COMM_vars8

  !-----------------------------------------------------------------------------
  subroutine COMM_wait(var, vid)
    implicit none

    real(8), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM wait MPI')

    !--- wait packets
    call MPI_WAITALL(ireq_cnt(vid), ireq_list(1:ireq_cnt(vid),vid), MPI_STATUSES_IGNORE, ierr)

    !--- unpacking packets from East
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


    call TIME_rapend  ('COMM wait MPI')

    return
  end subroutine COMM_wait

  !-----------------------------------------------------------------------------
  subroutine COMM_vars_r4(var, vid)
    use mod_process, only : &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    implicit none

    real(4), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ireqc, tag

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    tag = vid * 100
    ireqc = 1

    call TIME_rapstart('COMM vars(real4) MPI')

    !-- From 4-Direction HALO communicate
    ! From S
    call MPI_IRECV( var(:,:,JS-JHALO:JS-1), datasize_NS4,         &
                    MPI_REAL            , PRC_next(PRC_S), tag+1, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From N
    call MPI_IRECV( var(:,:,JE+1:JE+JHALO), datasize_NS4,         &
                    MPI_REAL            , PRC_next(PRC_N), tag+2, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From E
    call MPI_IRECV( recvpack_E2P_r4(:,vid), datasize_WE,          &
                    MPI_REAL            , PRC_next(PRC_E), tag+3, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From W
    call MPI_IRECV( recvpack_W2P_r4(:,vid), datasize_WE,          &
                    MPI_REAL            , PRC_next(PRC_W), tag+4, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    !-- To 4-Direction HALO communicate
    !--- packing packets to West
    do j = JS, JE
    do i = IS, IS+IHALO-1
    do k = 1, KA
        n =  (j-JS) * KA * IHALO &
           + (i-IS) * KA         &
           + k
        sendpack_P2W_r4(n,vid) = var(k,i,j)
    enddo
    enddo
    enddo

    !--- packing packets to East
    do j = JS, JE
    do i = IE-IHALO+1, IE
    do k = 1, KA
        n =  (j-JS)         * KA * IHALO &
           + (i-IE+IHALO-1) * KA         &
           + k
        sendpack_P2E_r4(n,vid) = var(k,i,j)
        n = n + 1
    enddo
    enddo
    enddo

    ! To W HALO communicate
    call MPI_ISEND( sendpack_P2W_r4(:,vid), datasize_WE,          &
                    MPI_REAL            , PRC_next(PRC_W), tag+3, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! To E HALO communicate
    call MPI_ISEND( sendpack_P2E_r4(:,vid), datasize_WE,          &
                    MPI_REAL            , PRC_next(PRC_E), tag+4, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! To N HALO communicate
    call MPI_ISEND( var(:,:,JE-JHALO+1:JE), datasize_NS4,         &
                    MPI_REAL            , PRC_next(PRC_N), tag+1, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! To S HALO communicate
    call MPI_ISEND( var(:,:,JS:JS+JHALO-1), datasize_NS4,         &
                    MPI_REAL            , PRC_next(PRC_S), tag+2, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ireq_cnt(vid) = ireqc - 1

    call TIME_rapend  ('COMM vars(real4) MPI')

    return
  end subroutine COMM_vars_r4

  !-----------------------------------------------------------------------------
  subroutine COMM_vars8_r4(var, vid)
    use mod_process, only : &
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

    real(4), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ireqc, tag, tagc

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    tag   = vid * 100
    ireqc = 1

    call TIME_rapstart('COMM vars(real4) MPI')

    !-- From 8-Direction HALO communicate
    ! From SE
    tagc = 0
    do j = JS-JHALO, JS-1
        call MPI_IRECV( var(1,IE+1,j), datasize_4C,                       &
                        MPI_REAL            , PRC_next(PRC_SE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From SW
    tagc = 10
    do j = JS-JHALO, JS-1
        call MPI_IRECV( var(1,IS-IHALO,j), datasize_4C,                   &
                        MPI_REAL            , PRC_next(PRC_SW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From NE
    tagc = 20
    do j = JE+1, JE+JHALO
        call MPI_IRECV( var(1,IE+1,j), datasize_4C,                       &
                        MPI_REAL            , PRC_next(PRC_NE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From NW
    tagc = 30
    do j = JE+1, JE+JHALO
        call MPI_IRECV( var(1,IS-IHALO,j), datasize_4C,                   &
                        MPI_REAL            , PRC_next(PRC_NW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From S
    tagc = 40
    do j = JS-JHALO, JS-1
        call MPI_IRECV( var(1,IS,j), datasize_NS8,                       &
                        MPI_REAL            , PRC_next(PRC_S), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
         ireqc = ireqc + 1
         tagc  = tagc  + 1
    enddo
    ! From N
    tagc = 50
    do j = JE+1, JE+JHALO
        call MPI_IRECV( var(1,IS,j), datasize_NS8,                       &
                        MPI_REAL            , PRC_next(PRC_N), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo
    ! From E
    call MPI_IRECV( recvpack_E2P_r4(:,vid), datasize_WE,           &
                    MPI_REAL            , PRC_next(PRC_E), tag+60, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1
    ! From W
    call MPI_IRECV( recvpack_W2P_r4(:,vid), datasize_WE,           &
                    MPI_REAL            , PRC_next(PRC_W), tag+70, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1

    !-- To 8-Direction HALO communicate
    !--- packing packets to West
    do j = JS, JE
    do i = IS, IS+IHALO-1
    do k = 1, KA
        n =  (j-JS) * KA * IHALO &
           + (i-IS) * KA         &
           + k
        sendpack_P2W_r4(n,vid) = var(k,i,j)
    enddo
    enddo
    enddo

    !--- packing packets to East
    do j = JS, JE
    do i = IE-IHALO+1, IE
    do k = 1, KA
        n =  (j-JS)         * KA * IHALO &
           + (i-IE+IHALO-1) * KA         &
           + k
        sendpack_P2E_r4(n,vid) = var(k,i,j)
        n = n + 1
    enddo
    enddo
    enddo

    ! To W HALO communicate
    call MPI_ISEND( sendpack_P2W_r4(:,vid), datasize_WE,           &
                    MPI_REAL            , PRC_next(PRC_W), tag+60, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1

    ! To E HALO communicate
    call MPI_ISEND( sendpack_P2E_r4(:,vid), datasize_WE,           &
                    MPI_REAL            , PRC_next(PRC_E), tag+70, &
                    MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr     )
    ireqc = ireqc + 1

    ! To N HALO communicate
    tagc = 40
    do j = JE-JHALO+1, JE
        call MPI_ISEND( var(1,IS,j), datasize_NS8,                       &
                        MPI_REAL            , PRC_next(PRC_N), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To S HALO communicate
    tagc = 50
    do j = JS, JS+JHALO-1
        call MPI_ISEND( var(1,IS,j), datasize_NS8,                       &
                        MPI_REAL            , PRC_next(PRC_S), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr       )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To NW HALO communicate
    tagc = 0
    do j = JE-JHALO+1, JE
        call MPI_ISEND( var(1,IS,j), datasize_4C,                         &
                        MPI_REAL            , PRC_next(PRC_NW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To NE HALO communicate
    tagc = 10
    do j = JE-JHALO+1, JE
        call MPI_ISEND( var(1,IE-IHALO+1,j), datasize_4C,                 &
                        MPI_REAL            , PRC_next(PRC_NE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To SW HALO communicate
    tagc = 20
    do j = JS, JS+JHALO-1
        call MPI_ISEND( var(1,IS,j), datasize_4C,                         &
                        MPI_REAL            , PRC_next(PRC_SW), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ! To SE HALO communicate
    tagc = 30
    do j = JS, JS+JHALO-1
        call MPI_ISEND( var(1,IE-IHALO+1,j), datasize_4C,                 &
                        MPI_REAL            , PRC_next(PRC_SE), tag+tagc, &
                        MPI_COMM_WORLD, ireq_list(ireqc,vid), ierr        )
        ireqc = ireqc + 1
        tagc  = tagc  + 1
    enddo

    ireq_cnt(vid) = ireqc - 1

    call TIME_rapend  ('COMM vars(real4) MPI')

    return
  end subroutine COMM_vars8_r4

  !-----------------------------------------------------------------------------
  subroutine COMM_wait_r4(var, vid)
    implicit none

    real(4), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ierr
    integer :: i, j, k, n
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM wait(real4) MPI')

    !--- wait packets
    call MPI_WAITALL(ireq_cnt(vid), ireq_list(1:ireq_cnt(vid),vid), MPI_STATUSES_IGNORE, ierr)

    !--- unpacking packets from East
    do j = JS, JE
    do i = IE+1, IE+IHALO
    do k = 1,  KA
       n = (j-JS)   * KA * IHALO &
         + (i-IE-1) * KA         &
         + k

       var(k,i,j) = recvpack_E2P_r4(n,vid)
    enddo
    enddo
    enddo

    !--- unpacking packets from West
    do j = JS, JE
    do i = IS-IHALO, IS-1
    do k = 1,  KA
       n = (j-JS)       * KA * IHALO &
         + (i-IS+IHALO) * KA         &
         + k

       var(k,i,j) = recvpack_W2P_r4(n,vid)
    enddo
    enddo
    enddo

    call TIME_rapend  ('COMM wait(real4) MPI')

    return
  end subroutine COMM_wait_r4

#ifdef _USE_RDMA
  !-----------------------------------------------------------------------------
  subroutine COMM_set_rdma_variable(var, vid)
    implicit none

    real(8), intent(in) :: var(:,:,:)
    integer, intent(in) :: vid
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** set RDMA ID:', vid-1

    call set_rdma_variable(var, vid-1);

    return
  end subroutine
#endif

#ifdef _USE_RDMA
  !-----------------------------------------------------------------------------
  subroutine COMM_rdma_vars(vid, num)
    implicit none

    integer, intent(in) :: vid
    integer, intent(in) :: num
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM RDMA')

    !--- put data
    call rdma_put(vid-1, num)

    call TIME_rapend  ('COMM RDMA')

    return
  end subroutine COMM_rdma_vars
#endif

#ifdef _USE_RDMA
  !-----------------------------------------------------------------------------
  subroutine COMM_rdma_vars8(vid, num)
    implicit none

    integer, intent(in) :: vid
    integer, intent(in) :: num
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM RDMA')

    !--- put data
    call rdma_put8(vid-1, num)

    call TIME_rapend  ('COMM RDMA')

    return
  end subroutine COMM_rdma_vars8
#endif

  !-----------------------------------------------------------------------------
  subroutine COMM_stats(var, varname)
    use mod_process, only : &
       PRC_nmax,   &
       PRC_myrank
    use mod_const, only : &
       CONST_UNDEF8, &
       CONST_UNDEF2
    implicit none

    real(8),          intent(inout) :: var(:,:,:,:)
    character(len=*), intent(in)    :: varname(:)

    logical :: halomask(KA,IA,JA)

    real(8), allocatable :: statval   (:,:,:)
    integer, allocatable :: statidx   (:,:,:,:)
    real(8), allocatable :: allstatval(:,:)
    integer, allocatable :: allstatidx(:,:,:)
    integer              :: vsize

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
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,1,p),       &
                       vsize*2,              &
                       MPI_DOUBLE_PRECISION, &
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

    return
  end subroutine COMM_stats

  !-----------------------------------------------------------------------------
  subroutine COMM_total( var, varname )
    use mod_process, only : &
       PRC_nmax,   &
       PRC_myrank
    use mod_const, only : &
       CONST_UNDEF8
    use mod_geometrics, only: &
       area    => GEOMETRICS_area,    &
       vol     => GEOMETRICS_vol,     &
       totarea => GEOMETRICS_totarea, &
       totvol  => GEOMETRICS_totvol
    implicit none

    real(8),           intent(inout) :: var(KA,IA,JA)
    character(len=*),  intent(in)    :: varname

    real(8) :: statval(0:PRC_nmax-1)
    real(8) :: allstatval
    real(8) :: ksize

    integer :: ierr
    integer :: k, i, j, p
    !---------------------------------------------------------------------------

    if ( COMM_total_doreport ) then

       statval(:) = CONST_UNDEF8
       ksize = size(var(:,:,:),1)

       statval(PRC_myrank) = 0.D0
       if ( ksize == KA ) then ! 3D
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             statval(PRC_myrank) = statval(PRC_myrank) + var(k,i,j) * vol(k,i,j)
          enddo
          enddo
          enddo
          statval(PRC_myrank) = statval(PRC_myrank) / totvol
       elseif( ksize == 1 ) then ! 2D
          do j = JS, JE
          do i = IS, IE
             statval(PRC_myrank) = statval(PRC_myrank) + var(1,i,j) * area(1,i,j)
          enddo
          enddo
          statval(PRC_myrank) = statval(PRC_myrank) / totarea
       endif

       if ( COMM_total_globalsum ) then
          call TIME_rapstart('COMM Bcast MPI')
          ! MPI broadcast
          do p = 0, PRC_nmax-1
             call MPI_Bcast( statval(p),           &
                             1,                    &
                             MPI_DOUBLE_PRECISION, &
                             p,                    &
                             MPI_COMM_WORLD,       &
                             ierr                  )
          enddo

          allstatval = 0.D0
          do p = 0, PRC_nmax-1
             allstatval = allstatval + statval(p)
          enddo
          allstatval = allstatval / PRC_nmax

          if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,1PE24.17)') &
                     '[', varname, '] SUM(global) =', allstatval

          call TIME_rapend  ('COMM Bcast MPI')
       else
          ! statistics on each node
          if( IO_L ) write(IO_FID_LOG,'(1x,A,A8,A,1PE24.17)') &
                     '[', varname, '] SUM(local)  =', statval(PRC_myrank)
       endif
    endif

    return
  end subroutine COMM_total

end module mod_comm
